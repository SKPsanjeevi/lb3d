#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h>
#include "hdfread.h"

int hdfread(char *FILE,
	    char *prefix,
	    FLOATNUM ****data) {
	
	//by A. Sarkar, ICP, 2007
	// from lbe/version5/utils/xdrtoold/hdfread.c
	// takes in filename of hdf file to be opened, returns pointer *t with values

	char	    *basename;
	hid_t       file_id, dataset_id,type_id;  /* identifiers */
	hid_t       dataspace;
	herr_t      status;
	hsize_t     dims[3];
	int         i, j, k,status_n; 
	int	    nx, ny, nz;

	FLOATNUM ***t;  

	basename = strchr(FILE, 95);
	sprintf(FILE,"%s%s", prefix, basename);
	//printf("file:%s\n", FILE);

	/* Open an existing file. */
	/* When error in code sometimes all *.h5 files were deleted */
	/* A copy of colour_sc02_t000010.h5 is called c10 in the same directory */
	file_id = H5Fopen(FILE, H5F_ACC_RDONLY, H5P_DEFAULT); 

	/* Open an existing dataset. */
	dataset_id = H5Dopen(file_id,"OutArray"); 

	/* Get dimensions */
	dataspace=H5Dget_space(dataset_id);
	status_n=H5Sget_simple_extent_dims(dataspace,dims,NULL);
 
	//printf("Dims: %lu x %lu x %lu \n",(unsigned long)(dims[0]),(unsigned long)(dims[1]), (unsigned long)(dims[2]) );    
	//    printf(".");

	nx=dims[0];
	ny=dims[1];
	nz=dims[2];
 
	t=(FLOATNUM ***)malloc(nx*sizeof(FLOATNUM **));
	if(!t){fprintf(stderr,"Error: (1) can't allocate memory \n");
	exit(1);}
	t[0]=(FLOATNUM **)malloc(nx*ny*sizeof(FLOATNUM *));
	if(!t[0]){fprintf(stderr,"Error: (2) can't allocate memory \n"); exit(1);}
	t[0][0]=(FLOATNUM *)malloc(nx*ny*nz*sizeof(FLOATNUM));
	if(!t[0][0]){fprintf(stderr,"Error: (3) can't allocate memory \n");exit(1);}
      
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

      
	*data=t;
 

	/* Do original xdrstats stuff. Easier to happen here for h5-files. */
	//    printf(".done\n");
    
	return(0);
}

int hdfread2(char *FILE,
	     FLOATNUM ****data) {
	     
	
	//by A. Sarkar, ICP, 2007
	// from lbe/version5/utils/xdrtools/hdfread.c
	// takes in filename of hdf file to be opened and pointer to array
	// returns 0 on success

	hid_t       file_id, dataset_id,type_id;  /* identifiers */
	hid_t       dataspace;
	herr_t      status;
	hsize_t     dims[3];
	int         i, j, k,status_n; 
	int	    nx, ny, nz;

	FLOATNUM ***t;  

	//basename = strchr(FILE, 95);
	//sprintf(FILE,"%s%s", prefix, basename);
	printf("file:%s\n", FILE);

	/* Open an existing file. */
	/* When error in code sometimes all *.h5 files were deleted */
	/* A copy of colour_sc02_t000010.h5 is called c10 in the same directory */
	file_id = H5Fopen(FILE, H5F_ACC_RDONLY, H5P_DEFAULT); 

	/* Open an existing dataset. */
	dataset_id = H5Dopen(file_id,"OutArray"); 

	/* Get dimensions */
	dataspace=H5Dget_space(dataset_id);
	status_n=H5Sget_simple_extent_dims(dataspace,dims,NULL);
 
	//printf("Dims: %lu x %lu x %lu \n",(unsigned long)(dims[0]),(unsigned long)(dims[1]), (unsigned long)(dims[2]) );    
	//    printf(".");

	nx=dims[0];
	ny=dims[1];
	nz=dims[2];
 
	t=(FLOATNUM ***)malloc(nx*sizeof(FLOATNUM **));
	if(!t){fprintf(stderr,"Error: (1) can't allocate memory \n");
	exit(1);}
	t[0]=(FLOATNUM **)malloc(nx*ny*sizeof(FLOATNUM *));
	if(!t[0]){fprintf(stderr,"Error: (2) can't allocate memory \n"); exit(1);}
	t[0][0]=(FLOATNUM *)malloc(nx*ny*nz*sizeof(FLOATNUM));
	if(!t[0][0]){fprintf(stderr,"Error: (3) can't allocate memory \n");exit(1);}
      
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

	
	*data=t;
 

	/* Do original xdrstats stuff. Easier to happen here for h5-files. */
	//    printf(".done\n");
    
	return(0);
}

int hdfwrite(char *FILE,
	     FLOATNUM ****data) {
	     
	
	//by A. Sarkar, ICP, 2007
	// from lbe/version5/utils/xdrtools/hdfread.c
	// takes in filename of hdf file to be opened and pointer to array
	// returns 0 on success

	hid_t       file_id, dataset_id,type_id;  /* identifiers */
	hid_t       dataspace;
	herr_t      status;
	hsize_t     dims[3];
	int         i, j, k,status_n; 
	int	    nx, ny, nz;

	FLOATNUM ***t;  

	//basename = strchr(FILE, 95);
	//sprintf(FILE,"%s%s", prefix, basename);
	printf("file:%s\n", FILE);

	/* Open an existing file. */
	/* When error in code sometimes all *.h5 files were deleted */
	/* A copy of colour_sc02_t000010.h5 is called c10 in the same directory */
	file_id = H5Fopen(FILE, H5F_ACC_RDONLY, H5P_DEFAULT); 

	/* Open an existing dataset. */
	dataset_id = H5Dopen(file_id,"OutArray"); 

	/* Get dimensions */
	dataspace=H5Dget_space(dataset_id);
	status_n=H5Sget_simple_extent_dims(dataspace,dims,NULL);
 
	//printf("Dims: %lu x %lu x %lu \n",(unsigned long)(dims[0]),(unsigned long)(dims[1]), (unsigned long)(dims[2]) );    
	//    printf(".");

	nx=dims[0];
	ny=dims[1];
	nz=dims[2];
 
	t=(FLOATNUM ***)malloc(nx*sizeof(FLOATNUM **));
	if(!t){fprintf(stderr,"Error: (1) can't allocate memory \n");
	exit(1);}
	t[0]=(FLOATNUM **)malloc(nx*ny*sizeof(FLOATNUM *));
	if(!t[0]){fprintf(stderr,"Error: (2) can't allocate memory \n"); exit(1);}
	t[0][0]=(FLOATNUM *)malloc(nx*ny*nz*sizeof(FLOATNUM));
	if(!t[0][0]){fprintf(stderr,"Error: (3) can't allocate memory \n");exit(1);}
      
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

	
	*data=t;
 

	/* Do original xdrstats stuff. Easier to happen here for h5-files. */
	//    printf(".done\n");
    
	return(0);
}

int getDims(char *FILE, int *nx, int *ny, int *nz)
{
	char	    *basename;
	hid_t       file_id, dataset_id,type_id;  /* identifiers */
	hid_t       dataspace;
	herr_t      status;
	hsize_t     dims[3];
	int         i, j, k,status_n; 

	file_id = H5Fopen(FILE, H5F_ACC_RDONLY, H5P_DEFAULT); 

	/* Open an existing dataset. */
	dataset_id = H5Dopen(file_id,"OutArray"); 

	/* Get dimensions */
	dataspace=H5Dget_space(dataset_id);
	status_n=H5Sget_simple_extent_dims(dataspace,dims,NULL);

	*nx=dims[0];
	*ny=dims[1];
	*nz=dims[2];

	return(0);
}
