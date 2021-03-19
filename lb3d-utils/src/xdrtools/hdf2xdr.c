/* 
 *  PROGRAM BY E. BREITMOSER, EPCC
 *
 * Reads hdf5 file and converts it to xdr file of same name but xdr-extension
 *
 * ./a.out colour_pl03_t000800 ; write hdf-file without .h5 extension!!!
 */

#include <stdio.h>
#include <stdlib.h>
#include <rpc/rpc.h>
#include <math.h>
#include <hdf5.h>

#define MAXSTRING 160

int main (int argc, char *argv[])
{
        hid_t       file_id, dataset_id;  /* identifiers */
        hid_t       datatype,dataspace;
        hid_t       class2;
        H5T_class_t class;
        H5T_order_t order;
        size_t      size;
        hsize_t     dims[3];
        herr_t      status;
        int         status_n;
	FILE	    *f1,*f2;
	int         i,j,k;
	XDR	    xdrs;
	char        *fname;
        char        extxdr[] = ".xdr",exthdf[] = ".h5";
        char        xdrfile[MAXSTRING],hdffile[MAXSTRING];
        char        command[MAXSTRING];
	int         nx,ny,nz;
#ifdef SGL
        float       ***t;
#else 
        double       ***t;
#endif

	if (argc!=2) {
		fprintf(stderr,"Syntax: %s <filename> \n",argv[0]);
		return -1;
	}

	fname=argv[1];
	/*	 nx=atoi(argv[2]);
	 ny=atoi(argv[3]);
	 nz=atoi(argv[4]); */
   


	    /*        for(i=0;i<nx;i++){
	  for(j=0;j<ny;j++){
	    for(k=0;k<nz;k++){
              t[i][j][k]=23.678;
               printf("Hallo %d %d %d \n",i,j,k); 
		  }}} */

        sprintf(hdffile,"%s%s",fname,exthdf);
        printf("Reading %s \n",hdffile);
	/*  if(NULL==(f1=fopen(hdffile,"r"))){
	  perror("Unable to open input file");
	  return -1;
	  }*/

       /* Open an existing hdf-file. */
          file_id = H5Fopen(hdffile, H5F_ACC_RDONLY, H5P_DEFAULT); 
       /* Open an existing dataset. */
          dataset_id = H5Dopen1(file_id,"OutArray"); 

        datatype=H5Dget_type(dataset_id);
	class=H5Tget_class(datatype);
         if(class==H5T_FLOAT)printf("Data set has float type\n");

	 class=H5Tget_native_type(datatype,H5T_DIR_DESCEND);
         if(class!=H5T_NATIVE_FLOAT)printf("Data set not of type float\n");
         if(class!=H5T_NATIVE_DOUBLE)printf("Data set not of type double\n");
         if(class!=H5T_FLOAT)printf("Data set not of any type\n");

	order=H5Tget_order(datatype);
	size=H5Tget_size(datatype);
        printf("Data size is %d bytes \n",size);
	dataspace=H5Dget_space(dataset_id);
	status_n=H5Sget_simple_extent_dims(dataspace,dims,NULL);
       
        printf("Dims: %lu x %lu x %lu\n",(unsigned long)(dims[0]),(unsigned long)(dims[1]), (unsigned long)(dims[2]));


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
          t[0][0]=(float *)calloc(nx*ny*nz,sizeof(float));
          if(!t[0][0]){fprintf(stderr,"Error: (3) can't allocate memory \n");exit(1);}
#else 
          t=(double ***)malloc(nx*sizeof(double **));
          if(!t){fprintf(stderr,"Error: (1) can't allocate memory \n");
          exit(1);}
          t[0]=(double **)malloc(nx*ny*sizeof(double *));
          if(!t[0]){fprintf(stderr,"Error: (2) can't allocate memory \n"); exit(1);}
          t[0][0]=(double *)calloc(nx*ny*nz,sizeof(double));
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

       /* Read  the dataset. */
   
	   status = H5Dread(dataset_id,class, H5S_ALL, H5S_ALL, 
	     H5P_DEFAULT,&t[0][0][0]); 

       /* Close the dataset. */
          status = H5Dclose(dataset_id); 

       /* Close the hdf file. */
         status = H5Fclose(file_id);

        sprintf(xdrfile,"%s%s",fname,extxdr);

	/* Create xdrfile, stop if it already exists/overwrite? */
        sprintf(command,"touch %s ",xdrfile);
        system(command);
        if(NULL==(f2=fopen(xdrfile,"w"))){
	  perror("Unable to open output file");
	  return -1;
	}        
        printf("%s \n",xdrfile);

	  xdrstdio_create(&xdrs,f2,XDR_ENCODE); 
    
        for(i=0;i<nx;++i){
	  for(j=0;j<ny;++j){
	    for(k=0;k<nz;++k){
	      /* printf("The end: %f %d %d %d\n",t[i][j][k],i,j,k);  */
#ifdef SGL  
              xdr_float(&xdrs,&t[i][j][k]);
#else
	      xdr_double(&xdrs,&t[i][j][k]); 
#endif
	}
	}
	}
	     xdr_destroy(&xdrs); 
       
	       fclose(f2); 
	       /* free((char *)(t[0][0]));
        free((char *)(t[0]));
        free((char *)(t)); */
	return 0;
}

