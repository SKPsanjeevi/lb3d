#include <stdio.h>
#include <stdlib.h>
#include <rpc/rpc.h>
#include <rpc/xdr.h>
#include <math.h>
#include <time.h>
#include <set>

#ifdef OLD_HEADER_FILENAME
#include <iostream.h>
#else
#include <iostream>
#endif
#include <string>

#ifndef H5_NO_NAMESPACE
#ifndef H5_NO_STD
  using std::cout;
  using std::endl;
#endif  // H5_NO_STD
#endif

#include "H5Cpp.h"

#define DATASET "OutArray"

#define MAXSTRING 160

#define XAXIS 1
#define YAXIS 2
#define ZAXIS 3

#define RANK 3

#ifndef H5_NO_NAMESPACE
  using namespace H5;
#endif

using namespace std;

int main (int argc, char *argv[]) {
  char fname[MAXSTRING], fname_out[MAXSTRING];
  char *axis_char;
  int wall, iolet, axis;
  int nx, ny, nz;
  float ***pos;
  float w;

  //HDF variable declarations
  hid_t           file, file_out, org_space, ext_space, dset_org, dset_ext, dcpl;    // Handles
  herr_t          status;
  H5D_layout_t    layout;
  hsize_t         dims_org[RANK],
                  dims_ext[RANK],
                  offset[RANK],
                  offset_out[RANK] = {0, 0, 0},
                  stride[RANK] = {1, 1, 1},
                  count[RANK] = {1, 1, 1},
                  chunk[RANK] = {1, 1, 1}; // This is a default value and might be overwritten later if the original dataset was chunked

  fprintf(stdout,"=== ROck Sample EMbedder  === \n");

  if (argc!=3) {
    fprintf(stderr,"Syntax: %s <source filename> <destination filename> \n", argv[0]);
    exit(1);
  }

  sprintf(fname,"%s",argv[1]);

  fprintf(stdout,"Attempting to read hdf file %s \n", fname);
  file = H5Fopen (fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file < 0) {
    fprintf(stderr,"Invalid file, aborting \n");
    exit(1);
  }

  dset_org = H5Dopen(file, DATASET, H5P_DEFAULT);
  org_space = H5Dget_space(dset_org);
  H5Sget_simple_extent_dims(org_space,dims_org,NULL);
  fprintf (stdout,"Dataspace dimensions are %d x %d x %d \n", (int) dims_org[2], (int) dims_org[1], (int) dims_org[0]);

  dcpl = H5Dget_create_plist(dset_org);
  layout = H5Pget_layout(dcpl);
  printf("Storage layout for %s is: ", DATASET);
  switch (layout) {
      case H5D_COMPACT:
        printf ("H5D_COMPACT\n");
        break;
      case H5D_CONTIGUOUS:
        printf ("H5D_CONTIGUOUS\n");
        break;
      case H5D_CHUNKED:
        printf ("H5D_CHUNKED\n");
        H5Pget_chunk(dcpl,3,chunk);
        printf ("Chunk dimensions are %d x %d x %d\n", (int) chunk[2], (int) chunk[1], (int) chunk[0]);
        break;
      case H5D_LAYOUT_ERROR:
        printf("H5D_LAYOUT_ERROR\n");
        break;
      case H5D_NLAYOUTS:
        printf("H5D_NLAYOUTS\n");
        break;
  }

  fprintf(stdout,"Create tube in which direction [x|y|z]?\n");
  scanf("%s",axis_char);
  if (strcmp(axis_char,"x") == 0) {
    axis = XAXIS;
  }
  else if (strcmp(axis_char,"y") == 0) {
    axis = YAXIS;
  }
  else if (strcmp(axis_char,"z") == 0) {
    axis = ZAXIS;
  }
  else {
    fprintf(stderr,"Unknown axis, aborting\n");
    exit(1);
  }
  fprintf(stdout,"Enter wall width:\n");
  scanf("%d",&wall);
  fprintf(stdout,"Enter inlet/outlet size:\n");
  scanf("%d",&iolet);
  fprintf(stdout,"Enter wettability of the rock tube:\n");
  scanf("%f",&w);

  /* Depending on the chosen orientation of the rock tube, calculate offsets and dimensions */
  switch(axis) {
    case XAXIS:
      dims_ext[0] = dims_org[0] + 2*wall;
      dims_ext[1] = dims_org[1] + 2*wall;
      dims_ext[2] = dims_org[2] + 2*iolet;
      offset[0] = wall;
      offset[1] = wall;
      offset[2] = iolet;
      break;
    case YAXIS:
      dims_ext[0] = dims_org[0] + 2*wall;
      dims_ext[1] = dims_org[1] + 2*iolet;
      dims_ext[2] = dims_org[2] + 2*wall;
      offset[0] = wall;
      offset[1] = iolet;
      offset[2] = wall;
      break;
    case ZAXIS:
      dims_ext[0] = dims_org[0] + 2*iolet;
      dims_ext[1] = dims_org[1] + 2*wall;
      dims_ext[2] = dims_org[2] + 2*wall;
      offset[0] = iolet;
      offset[1] = wall;
      offset[2] = wall;
      break;
    default:
      fprintf(stderr,"Axis not well-defined. Aborting...\n");
      exit(1);
  }

  nz = dims_ext[0];
  ny = dims_ext[1];
  nx = dims_ext[2];

  fprintf(stdout,"New system size will be %d x %d x %d\n", nx, ny, nz);

  /* Allocate tensor (i.e. 3d array) and set pointers */
  fprintf(stdout,"Setting up 3d array... \n");
  pos = (float ***) malloc(nz*sizeof(float **));
  if (!pos) {
    fprintf(stderr,"Error: (1) can't allocate memory, aborting... \n");
    exit(1);
  }
  pos[0] = (float **) malloc(nz*ny*sizeof(float *));
  if (!pos[0]) { 
    fprintf(stderr,"Error: (2) can't allocate memory, aborting... \n");
    exit(1);
  }
  pos[0][0] = (float *) calloc(nz*ny*nx,sizeof(float));
  if (!pos[0][0]) {
    fprintf(stderr,"Error: (3) can't allocate memory, aborting... \n");
    exit(1);
  }

  /* Calculate vector pointers */
  for(int j=1;j<ny;j++) {
    pos[0][j] = pos[0][j-1]+nx;
  }
  for(int k=1;k<nz;k++) {
    pos[k]=pos[k-1]+ny;
    pos[k][0]=pos[k-1][0]+ny*nx;
    for(int j=1;j<ny;j++) {
      pos[k][j]=pos[k][j-1]+nx;
    }
  }

  fprintf(stdout,"Finished setting up 3d array \n");

  fprintf(stdout,"Creating tube... \n");

  /*
  We set up the array 'pos' in one go by setting zeroes or the wettability for all sites,
  depending on the requested orientation of the tube.
  */

  for(int k=0;k<nz;k++) {
    for(int j=0;j<ny;j++) {
      for(int i=0;i<nx;i++) {
        switch(axis) {
          case XAXIS:
            if ( (j < wall || j >= ny - wall) ||
                 (k < wall || k >= nz - wall) ) {
               pos[k][j][i] = w;
            }
            else {
              pos[k][j][i] = 0.0;
            }
            break;
          case YAXIS:
            if ( (i < wall || i >= nx - wall) ||
                 (k < wall || k >= nz - wall) ) {
               pos[k][j][i] = w;
            }
            else {
              pos[k][j][i] = 0.0;
            }
            break;
          case ZAXIS:
            if ( (i < wall || i >= nx - wall) ||
                 (j < wall || j >= ny - wall) ) {
               pos[k][j][i] = w;
            }
            else {
              pos[k][j][i] = 0.0;
            }
            break; 
          default:
            fprintf(stderr,"Axis not well-defined. Aborting...\n");
            exit(1);
        }
      }
    }
  }
  fprintf(stdout,"Finished creating tube\n");

  /*
  Array 'pos' has already been set up with zeroes and the rock tube;
  now we are going to place the original data into the same array with the correct offset.
  This will save us from holding a second (large) array in memory.
  */

  ext_space = H5Screate_simple(RANK,dims_ext, NULL);

  fprintf(stdout,"Reading original file into tube\n");
  status = H5Sselect_hyperslab(ext_space, H5S_SELECT_SET, offset, NULL, count, dims_org);
  status = H5Dread(dset_org, H5T_NATIVE_FLOAT, ext_space, org_space, H5P_DEFAULT, pos[0][0]);

  //sprintf(fname_out,"%s-emb%s-%d-%d.h5", argv[1], axis_char, wall, iolet);
  sprintf(fname_out,"%s",argv[2]);

  /*
  Array 'pos' now holds all the correct data, so we just need to set up chunking
  (copied from the chunks in the in-file), and select the entire dataset.
  Finally, the data is written to file.
  */

  fprintf(stdout,"Writing file %s...\n", fname_out);
  file_out = H5Fcreate(fname_out, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  dcpl = H5Pcreate(H5P_DATASET_CREATE);
  //status = H5Pset_chunk(dcpl, 3, chunk);
  dset_ext = H5Dcreate(file_out, DATASET, H5T_NATIVE_FLOAT, ext_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Sselect_hyperslab(ext_space, H5S_SELECT_SET, offset_out, stride, count, dims_ext);
  status = H5Dwrite(dset_ext, H5T_NATIVE_FLOAT, H5S_ALL, ext_space, H5P_DEFAULT, pos[0][0]);

  /* Close all the HDF5 handles */

  status = H5Pclose (dcpl);
  status = H5Dclose (dset_ext);
  status = H5Dclose (dset_org);
  status = H5Sclose (org_space);
  status = H5Sclose (ext_space);
  status = H5Fclose (file);
  status = H5Fclose (file_out);

  fprintf(stdout,"Done! \n");

  return 0;
}