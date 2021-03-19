/* 
 *  Reading the existing filename for scalars handed over by xdrstats.c.
 *  Only for floats, not for doubles.
 * Finds out about the dimensions by itself, given it is a 3d array.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <rpc/rpc.h>
#include <hdf5.h>

void hdfread(char *FILE) {

   hid_t       file_id, dataset_id,type_id;  /* identifiers */
   hid_t       dataspace;
   herr_t      status;
   hsize_t     dims[3];
   int         i, j, k,status_n,nx,ny,nz;
#ifdef SGL
   float       phimax,phimin,phisum,phisum2;
   float       sd;
   float       ***t;
#else
   double      phimax,phimin,phisum,phisum2;
   double      sd;
   double      ***t;
#endif

   /* Open an existing file. */
   /* When error in code sometimes all *.h5 files were deleted */
   /* A copy of colour_sc02_t000010.h5 is called c10 in the same directory */
    file_id = H5Fopen(FILE, H5F_ACC_RDONLY, H5P_DEFAULT); 

   /* Open an existing dataset. */
      dataset_id = H5Dopen1(file_id,"OutArray"); 

      /* Get dimensions */
      dataspace=H5Dget_space(dataset_id);
      status_n=H5Sget_simple_extent_dims(dataspace,dims,NULL);
 
     printf("Dims: %lu x %lu x %lu \n",(unsigned long)(dims[0]),(unsigned long)(dims[1]), (unsigned long)(dims[2]) );    

      nx=dims[0];
      ny=dims[1];
      nz=dims[2];
 
        /* Allocate tensor (i.e. 3d array) and set pointers */
#ifdef SGL
          t=(float ***)malloc(nx*sizeof(float **));
          if(!t){fprintf(stderr,"Error: (1) can't allocate memory \n");
          exit(1);}
          t[0]=(float **)malloc(nx*ny*sizeof(float *));
          if(!t[0]){fprintf(stderr,"Error: (2) can't allocate memory \n"); exit(1);}
          t[0][0]=(float *)malloc(nx*ny*nz*sizeof(float));
          if(!t[0][0]){fprintf(stderr,"Error: (3) can't allocate memory \n");exit(1);}
#else
          t=(double ***)malloc(nx*sizeof(double **));
          if(!t){fprintf(stderr,"Error: (1) can't allocate memory \n");
          exit(1);}
          t[0]=(double **)malloc(nx*ny*sizeof(double *));
          if(!t[0]){fprintf(stderr,"Error: (2) can't allocate memory \n"); exit(1);}
          t[0][0]=(double *)malloc(nx*ny*nz*sizeof(double));
          if(!t[0][0]){fprintf(stderr,"Error: (3) can't allocate memory \n");exit(1);}
#endif
          
          for(j=1;j<ny;j++){t[0][j]=t[0][j-1]+nz;}

            for(i=1;i<nx;i++){
            t[i]=t[i-1]+ny;
            t[i][0]=t[i-1][0]+ny*nz;
            for(j=1;j<ny;j++){
              t[i][j]=t[i][j-1]+nz;
            }
            } 

#ifdef SGL
	    type_id=H5T_NATIVE_FLOAT;
#else
	    type_id=H5T_NATIVE_DOUBLE;
#endif

   /* Read  the dataset. */
   
      status = H5Dread(dataset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
                   &t[0][0][0]);
   
   /* Close the dataset. */
      status = H5Dclose(dataset_id); 

   /* Close the file. */
      status = H5Fclose(file_id);
 

   /* Do original xdrstats stuff. Easier to happen here for h5-files. */
phisum=phimax=phimin=t[0][0][0];
phisum2=phisum*phisum;

   for(i=0;i<nx;++i){
    for(j=0;j<ny;++j){
     for(k=0;k<nz;++k){
       if(i!=0||j!=0||k!=0){		
      phisum+=t[i][j][k];			
      phisum2+=t[i][j][k]*t[i][j][k];			
      if (t[i][j][k]>phimax) { phimax=t[i][j][k];}
      if (t[i][j][k]<phimin) { phimin=t[i][j][k];}
  } }
  }
 }

#ifdef SGL
	phisum/=(float)(nx*ny*nz);
	phisum2/=(float)(nx*ny*nz);
#else
	phisum/=(double)(nx*ny*nz);
	phisum2/=(double)(nx*ny*nz);
#endif

	sd=sqrt(phisum2-phisum*phisum);
   
	printf("Max %f Min %f Mean %f SD %f\n",
		phimax,phimin,phisum,sd);

        free((char *)(t[0][0]));
        free((char *)(t[0]));
        free((char *)(t)); 

return;
}
