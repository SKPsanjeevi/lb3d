#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <string.h>

#include "hdf5_helper.h"
#include "lbe_version.h"

using namespace std;

int main (int argc, char *argv[]) {
  // Dimensions
  int nx, ny, nz, nv;
  int n = 0;

  int dlen;

  // File names
  char *fname;
  char  fname_out[256];

  hid_t file_out, dcpl, ext_space, dset_ext;
  herr_t status;

  hsize_t offset[3] = {0,0,0};
  hsize_t stride[3] = {1,1,1};
  hsize_t count[3] = {1,1,1};
  
  char fstr[32];

  char dset_name[] = "OutArray";

  // Array pointers
  float *v;

  hsize_t *dims;

  hsize_t out_dims[3];
  // Flags
  bool fname_set = false;

  for (int i = 1; i < argc ; i++) {

    // fprintf(stdout,"  argc = %d, argv[%d] = <%s>\n",argc,i,argv[i]);

    if ( strcmp(argv[i],"-f") == 0 ) {
      if ( i+1 < argc ) {
	fname = argv[++i];
	//fprintf(stdout,"  File <%s>\n",fname);
	fname_set = true;
      }
      else {
	fprintf(stderr,"  Missing argument to -f flag. \n");
	exit(1);
      }
    }
    else {
      //fprintf(stdout,"  No matches for <%s>.\n",argv[i]);
    }
  }

  if (! (fname_set ) ) {
    fprintf(stdout,"arr-to-velz (%s)\n\n", GIT_DESC);
    fprintf(stderr,"File -f is mandatory.\n\n");
    fprintf(stderr,"  -f <filename>\n");
    fprintf(stderr,"\n");
    exit(1);
  }


  if (hdfReadAllVecComp(fname,&v,&dims,2) != 0) {
    fprintf(stderr,"ERROR\n");
    exit(1);
  }
  hdfToIntDims(dims,&nx,&ny,&nz,&nv);
  hdfFreeDims(dims);
  fprintf(stdout,"Found %d x %d x %d (%d) dataset.\n", nx, ny, nz, nv);

  // for(int k=0;k<nz;k++) {
  //   for(int j=0;j<ny;j++) {
  //     for(int i=0;i<nx;i++) {
  // 	fprintf(stdout,"%d %d %d\n",i,j,k);
  // 	fprintf(stdout,"         = %f\n",v[k*ny*nx+j*nx+i]);
  //     }
  //   }
  // }

  out_dims[0] = nz;
  out_dims[1] = ny;
  out_dims[2] = nx;

  sprintf(fname_out,"velz_%s",fname+4);
  fprintf(stdout,"Writing HDF5 file %s...\n", fname_out);
  file_out = H5Fcreate(fname_out, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  dcpl = H5Pcreate(H5P_DATASET_CREATE);
    // fprintf(stdout,"Setting chunks to %d x %d x %d\n", (int) chunk[0], (int) chunk[1], (int) chunk[2]);
    // status = H5Pset_chunk(dcpl, 3, chunk);
  ext_space = H5Screate_simple(3,out_dims, NULL);
  dset_ext = H5Dcreate(file_out, dset_name, H5T_NATIVE_FLOAT, ext_space, H5P_DEFAULT, dcpl, H5P_DEFAULT);
  status = H5Sselect_hyperslab(ext_space, H5S_SELECT_SET, offset, stride, count, out_dims);
  status = H5Dwrite(dset_ext, H5T_NATIVE_FLOAT, H5S_ALL, ext_space, H5P_DEFAULT, v);
  fprintf(stdout,"Finished writing HDF5 file\n");

  free(v);



  //hdfFree4DArray(v);

  return 0;
}
