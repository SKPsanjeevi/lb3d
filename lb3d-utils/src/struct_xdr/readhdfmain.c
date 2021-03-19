/* 
 *  Reading the existing dataset for scalars handed over by xdrstats.c.
 *  Only for floats, not for doubles.
 */

#include <stdio.h>
#include <float.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <rpc/rpc.h>
#include <hdf5.h>
#include "main.h"

int hdfread(Params *params,Dat *dat,Stats *stats,char *FILE) {

   hid_t       file_id, dataset_id;  /* identifiers */
   herr_t      status;
   int         i, j, k;
   /* Use malloc so that array goes to heap, stack too small for 
      large array sizes; will produce a segmentation fault fllowed by
      a core dump */
   float       (*dset_data)[params->nx][params->ny]=(float(*)[params->nx][params->ny]) malloc(sizeof(int [params->nx][params->ny])*params->nz);; 

   /*     fprintf(stderr,"nx: %d %s \n",params->nx,FILE); */
   /* Open an existing file. */
   /* When error in code sometimes all *.h5 files were deleted */
   /* A copy of colour_sc02_t000010.h5 is called c10 in the same directory */
    file_id = H5Fopen(FILE, H5F_ACC_RDWR, H5P_DEFAULT); 

   /* Open an existing dataset. */
      dataset_id = H5Dopen1(file_id,"OutArray"); 

   /* Read  the dataset. */
   
      status = H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
                    dset_data);
   
   /* Close the dataset. */
      status = H5Dclose(dataset_id); 

   /* Close the file. */
      status = H5Fclose(file_id);
 

   fprintf(stderr,"hdf file closed \n");
   /* Do original decode_FLOAT stuff. Easier to happen here for h5-files. */

      stats->min = stats->min2 =  DBL_MAX;
      stats->max = stats->max2 = -DBL_MAX;

   for(i=0;i<=params->nx-1;++i){
    for(j=0;j<=params->ny-1;++j){
     for(k=0;k<=params->nz-1;++k){
      dat->data[i+1][j+1][k+1] =(double)dset_data[i][j][k];	       	
      if(dset_data[i][j][k] > stats->max) stats->max = dset_data[i][j][k];
      if(dset_data[i][j][k] < stats->min) stats->min = dset_data[i][j][k];
   }
  }
 }
 /*   fprintf(stderr,"stats_max %f,stats_min %f \n",stats->max,stats->min);*/
   
	free(dset_data);
return(1);
}
