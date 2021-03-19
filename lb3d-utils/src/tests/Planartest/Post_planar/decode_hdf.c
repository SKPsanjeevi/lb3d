/* Read HDF from the specified file, write ASCII to stdout.
*  C Implementation: decode_hdf
*
* Description: 
*
*
* Author: Frank Raischel, 210 (+49 711 685-3594) <raischel@ica1.uni-stuttgart.de>, (C) 2007
*
* Copyright: See COPYING file that comes with this distribution
*
*/



#include<stdio.h>
#include<rpc/rpc.h>
#include<sys/types.h>
#include<sys/stat.h>
#include<hdf5.h>

#include "main.h"

extern off_t   file_size(char *);
extern int     extrema(Stats *, float *);


int decode_hdf(char *prefix, Params *params, unsigned long timestep,
	   float ***data, Stats *stats)
{                      
int     i, j, k, nx,ny,nz;
float   a=0.0;	
float   val=0.0;
float ***t;
   char    *hdffile;
   FILE	   *infile;
   FILE *READPIPE;
   
       hid_t       file_id, dataset_id;  /* identifiers */
       hid_t       datatype,dataspace;
       hid_t       class2;
      hid_t class;
      H5T_order_t order;
      size_t      size;
      hsize_t     dims[3];
       herr_t      status;
			 int         status_n;  

   char	   *hdftestfile;
   char    *filetest;
   char	    *alt_filename;
   
   hdffile = (char *)malloc(PATHNAMELEN * sizeof(char));
   hdftestfile = (char *)malloc(PATHNAMELEN * sizeof(char));
   alt_filename = (char *)malloc(PATHNAMELEN * sizeof(char));
   filetest = (char *)malloc((PATHNAMELEN +20) * sizeof(char));

   sprintf(hdffile, "%s/%s_%s_t%06lu.h5",  params->datapath, prefix, params->gr_out_file, timestep);

   sprintf(hdftestfile, "%s/%s_%s_t%06lu-*.h5", params->datapath, prefix, params->gr_out_file, timestep);	
   sprintf(filetest, "ls %s", hdftestfile);	


   if(file_size(hdffile) == 0) {
	   if( NULL==(READPIPE=popen(filetest, "r"))){
		   fprintf(stderr, "Cannot open pipe with command ' %s '!\n", filetest);
		   return(-1);
	   }
	   else{
			fscanf(READPIPE, "%s",alt_filename );
			if(file_size(alt_filename) == 0) {
				fprintf(stderr, "\n\tError: File %s has zero length\n", alt_filename);
				fflush(stderr);
				free(hdffile);
				/*  return(1) is decode() OK, 
				*  return(0) is zero filelength, 
				*  return(-1) is error */
				return(0);
				}
			else{
				sprintf(hdffile, "%s", alt_filename);
			}	
	   }
   }

	 //FR:
	 // Not needed, use H5Fopen() instead:
/*    if (NULL == (infile = fopen(hdffile,"r"))) { */
/*       fprintf(stderr, "\n\tError: Could not open data file %s\n",   hdffile); */
/*       fprintf(stderr, "\a\nError in reading size of file %s\n" */
/* 		      "\tIn function tools.c:file_size() called by decode()\n", */
/* 		      hdffile); */
/*       fflush(stderr); */
/*       free(hdffile); */
/*       return(-1); */
/*    } */

   printf("\n\tReading file %s", hdffile);

	 /* Swtup hdf stream */

    /* Open an existing hdf-file. */
       file_id = H5Fopen(hdffile, H5F_ACC_RDONLY, H5P_DEFAULT);
       /* Open an existing dataset. */
        dataset_id = H5Dopen1(file_id,"OutArray");

        datatype=H5Dget_type(dataset_id);
        class=H5Tget_class(datatype);
         if(class==H5T_FLOAT)printf("Data set has float type\n");

         class=H5Tget_native_type(datatype,H5T_DIR_DESCEND);
         //if(class!=H5T_NATIVE_FLOAT)printf("Data set not of type float\n");
         //if(class!=H5T_NATIVE_DOUBLE)printf("Data set not of type double\n");
         //if(class!=H5T_FLOAT)printf("Data set not of any type\n");

        order=H5Tget_order(datatype);
        size=H5Tget_size(datatype);
        printf("Data size is %lu bytes \n",size);
        dataspace=H5Dget_space(dataset_id);
        status_n=H5Sget_simple_extent_dims(dataspace,dims,NULL);

        printf("Dims: %lu x %lu x %lu\n",(unsigned long)(dims[0]),(unsigned long)(dims[1]), (unsigned long)(dims[2]));


   //pxx = f3tensor(1, dims[0], 1,dims[1], 1, dims[2]);	
        nx=dims[0];
	ny=dims[1];
	nz=dims[2];

	/* Allocate tensor (i.e. 3d array) and set pointers */


          t=(float ***)malloc(nx*sizeof(float **));
          if(!t){fprintf(stderr,"Error: (1) can't allocate memory \n");
          exit(1);}
          t[0]=(float **)malloc(nx*ny*sizeof(float *));
          if(!t[0]){fprintf(stderr,"Error: (2) can't allocate memory \n"); exit(1);}
          t[0][0]=(float *)calloc(nx*ny*nz,sizeof(float));
          if(!t[0][0]){fprintf(stderr,"Error: (3) can't allocate memory \n");exit(1);}
          
          for(j=1;j<ny;j++){t[0][j]=t[0][j-1]+nz;}

            for(i=1;i<nx;i++){
      	    t[i]=t[i-1]+ny; 
            t[i][0]=t[i-1][0]+ny*nz;
            for(j=1;j<ny;j++){
	      t[i][j]=t[i][j-1]+nz;
	    }
	    } 


				/* Read hdf data */

				//status = H5Dread(dataset_id,class, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0][0][0]);
			status = H5Dread(dataset_id,class, H5S_ALL, H5S_ALL, H5P_DEFAULT, &t[0][0][0]);  
				

				//!
				/* Extrema missing */
				//!
       /* Close the dataset. */
          status = H5Dclose(dataset_id); 

       /* Close the hdf file. */
         status = H5Fclose(file_id);
				 
				 // already closed		
				 //			 fclose(hdffile);


//now, find the extrema, and add here code for changing x,y,z , if necessary, since:
// CAVE: hdf coordinate ordering differs between C and Fortran!!!!

//this is a very ugly hack, but it works.
if(nx==params->nx && ny==params->ny && nz==params->nz){
 for(i=0;i< params->nz;++i){
    for(j=0;j< params->ny;++j){
     for(k=0;k< params->nx;++k){
       if(i!=0||j!=0||k!=0){
	val=t[i][j][k];
	extrema(stats, &val);  
	  }
	data[i+1][j+1][k+1]=t[i][j][k];
     }
    }
  }
}
else if(nz==params->nx && ny==params->ny && nx==params->nz){
 for(i=0;i< params->nx;++i){
    for(j=0;j< params->ny;++j){
     for(k=0;k< params->nz;++k){
       if(i!=0||j!=0||k!=0){
	val=t[i][j][k];
	extrema(stats, &val);  
	  }
	data[i+1][j+1][k+1]=t[k][j][i];
     }
    }
  }
}
else{
fprintf(stderr, "Unknown ordering in file %s !\n", hdffile);
return(-1);
}


return(1);
}//end decode_hdf()


