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
#define INDX(nx,ny,nz,d,i,j,k,l) i*ny*nz*d+j*nz*d+k*d+l

int main (int argc, char *argv[])
{
        hid_t       file_id, dataset_id,type_id;  /* identifiers */
        hid_t       datatype,dataspace,file_dataspace,mem_dataspace;
        H5T_class_t class;
        H5T_order_t order;
        size_t      size;
        hsize_t     dims[4];
        herr_t      status;
        int         status_n;
	FILE	    *f1,*f2;
	int         i,j,k,l;
        XDR         xdrs;
	char        *fname;
        char        extxdr[] = ".xdr",exthdf[] = ".h5";
        char        xdrfile[MAXSTRING],hdffile[MAXSTRING];
        char        command[MAXSTRING];
	int         nx,ny,nz,d;
#ifdef SGL
        float       *t;
#else
        double      *t;
#endif
	if (argc!=2) {
		fprintf(stderr,"Syntax: %s <filename> \n",argv[0]);
		return -1;
	}

	fname=argv[1];

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

	order=H5Tget_order(datatype);
	size=H5Tget_size(datatype);
        printf("Data size is %d bytes \n",size);
	dataspace=H5Dget_space(dataset_id);
	file_dataspace=mem_dataspace=dataspace;
	status_n=H5Sget_simple_extent_dims(dataspace,dims,NULL);
       
        printf("Dims: %lu x %lu x %lu x %lu\n",(unsigned long)(dims[0]),(unsigned long)(dims[1]), (unsigned long)(dims[2]), (unsigned long)(dims[3]));


	nx=dims[0];
	ny=dims[1];
	nz=dims[2];
        d=dims[3];

        type_id=H5Tget_native_type(datatype,H5T_DIR_DESCEND);
#ifdef SGL
	t=malloc(sizeof(float)*dims[0]*dims[1]*dims[2]*dims[3]);
#else
	t=malloc(sizeof(double)*dims[0]*dims[1]*dims[2]*dims[3]);
#endif

        H5Dread(dataset_id,type_id,mem_dataspace,file_dataspace,
                H5P_DEFAULT,t);

       /* Close the dataset. */
          status = H5Dclose(dataset_id); 
	  H5Tclose(datatype);
	  H5Sclose(dataspace);
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
    
        for(i=0;i<(nx*ny*nz*d);i++){
	  /* printf("The end: %f %d\n",t[i],i); */
	
#ifdef SGL
     xdr_float(&xdrs,&t[i]); 
#else
     xdr_double(&xdrs,&t[i]);
#endif
	}
	     xdr_destroy(&xdrs); 
       
	       fclose(f2); 
	       /* free((char *)(t[0][0]));
        free((char *)(t[0]));
        free((char *)(t)); */
	return 0;
}

