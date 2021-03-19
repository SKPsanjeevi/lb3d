#include "h5converter.h"

#ifndef H5_NO_NAMESPACE
  using namespace H5;
#endif

using namespace std;

const int rank = 3;

void usage() {
  fprintf(stderr,"Flags -i and -o are required. \n");
  fprintf(stderr,"Either -s or -x, -y, and -z are required for XDR or RAW input. \n\n");
  fprintf(stderr,"  -d         <dataset name> (for HDF5 only, default: 'OutArray') \n");
  fprintf(stderr,"  -i         <inputfile>  <input file type> \n");
  fprintf(stderr,"  -o         <outputfile> <output file type> \n");
  fprintf(stderr,"  -s         <nx> <ny> <nz> \n");
  fprintf(stderr,"  -x -y -z   <n[x|y|z]> \n");
  fprintf(stderr,"\n");
  fprintf(stderr,"Valid file type specifications: 'hdf', 'xdr', 'raw'. \n");
  fprintf(stderr,"\n");
  exit(1);
}

int main (int argc, char *argv[]) {
  //  FILE  *f;
  int nx = 0;
  int ny = 0;
  int nz = 0;

  int cluster_id = 0;
  int reverse = 0;

  // File names
  char *fname_in, *fname_out;
  char *dset_name;

  bool fname_in_set = false;
  bool fname_out_set = false;

  int ***f_array;
  float ***o_array;

  //HDF variable declarations
  hid_t           file, file_out, org_space, ext_space, dset_org, dset_ext, dcpl;    // Handles
  herr_t          status;
  H5D_layout_t    layout;
  hsize_t         dims[rank],
                  offset[rank] = {0, 0, 0},
                  stride[rank] = {1, 1, 1},
                  count[rank] = {1, 1, 1},
                  chunk[rank] = {1, 1, 1}; // This is a default value and might be overwritten later

  fprintf(stdout,"\nH5 cluster splitter -- converts from HDF5 to HDF5. \n\n");

  dset_name = "OutArray";

  for (int i = 1; i < argc ; i++) {

    // fprintf(stdout,"  argc = %d, argv[%d] = '%s'. \n", argc, i, argv[i]);
    if ( strcmp(argv[i],"-?") == 0 ) {
      usage();
    }
    else if ( strcmp(argv[i],"-c") == 0 ) {
      if ( i+1 < argc ) {
	cluster_id = atoi(argv[++i]);
	//fprintf(stdout,"  Cluster id '%d'. \n", cluster_id);
      }
      else {
	fprintf(stderr,"  Missing arguments to -c <cluster id>. \n");
	exit(1);
      }
    }
    else if ( strcmp(argv[i],"-d") == 0 ) {
      if ( i+1 < argc ) {
	dset_name = argv[++i];
	//fprintf(stdout,"  Dataset '%s'. \n", dset_name);
      }
      else {
	fprintf(stderr,"  Missing arguments to -d <dsetname>. \n");
	exit(1);
      }
    }
    else if ( strcmp(argv[i],"-i") == 0 ) {
      if ( i+1 < argc ) {
	fname_in = argv[++i];
	//fprintf(stdout,"  Input file name '%s', type '%s'. \n", fname_in, ftype_in);
	fname_in_set = true;
      }
      else {
	fprintf(stderr,"  Missing arguments to -i <filename> <filetype>. \n");
	exit(1);
      }
    }
    else if ( strcmp(argv[i],"-o") == 0 ) {
      if ( i+1 < argc ) {
	fname_out = argv[++i];
	//fprintf(stdout,"  Output file name '%s', type '%s'. \n", fname_out, ftype_out);
	fname_out_set = true;
      }
      else {
	fprintf(stderr,"  Missing argument to -o <filename> <filetype>. \n");
	exit(1);
      }
    }
    else if (strcmp(argv[i],"-r") == 0) {
      reverse = 1;
    }
    else {
      //fprintf(stdout,"  No matches found for flag '%s'. \n", argv[i]);
    }
  }

  if ( !fname_in_set || !fname_out_set) {
    //fprintf(stderr,"No input or output file specified. Aborting.\n");
    usage();
  }

  // Setting name of dataset
  //fprintf(stdout,"Set dset_name to '%s'.\n",dset_name);
  
  // HDF prep
  file = H5Fopen (fname_in, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file < 0) {
    fprintf(stderr,"[HDF5] Invalid file '%s'. Aborting. \n\n", fname_in);
    exit(1);
  }
    
  dset_org = H5Dopen(file, dset_name, H5P_DEFAULT);
  org_space = H5Dget_space(dset_org);
  H5Sget_simple_extent_dims(org_space,dims,NULL);
  fprintf (stdout,"[HDF5] Dataspace dimensions are %d x %d x %d.\n", (int) dims[2], (int) dims[1], (int) dims[0]);
    
  nx = (int) dims[2];
  ny = (int) dims[1];
  nz = (int) dims[0];

  dcpl = H5Dget_create_plist(dset_org);
  layout = H5Pget_layout(dcpl);
  fprintf(stdout, "[HDF5] Storage layout for '%s' is: ", dset_name);
  switch (layout) {
  case H5D_COMPACT:
    fprintf(stdout,"H5D_COMPACT\n");
    break;
  case H5D_CONTIGUOUS:
    fprintf(stdout,"H5D_CONTIGUOUS\n");
    break;
  case H5D_CHUNKED:
    fprintf(stdout,"H5D_CHUNKED");
    H5Pget_chunk(dcpl,3,chunk);
    fprintf(stdout," - chunk dimensions are %d x %d x %d\n", (int) chunk[2], (int) chunk[1], (int) chunk[0]);
    break;
  case H5D_LAYOUT_ERROR:
    fprintf(stdout,"H5D_LAYOUT_ERROR\n");
    break;
  case H5D_NLAYOUTS:
    fprintf(stdout,"H5D_NLAYOUTS\n");
    break;
  }

  /* Allocate tensor (i.e. 3d array) and set pointers */
  fprintf(stdout,"Setting up 3d array... \n");
  f_array = (int ***) malloc(nz*sizeof(int **));
  if (!f_array) {
    fprintf(stderr,"Error: (1) cannot allocate memory. Aborting. \n\n");
    exit(1);
  }
  f_array[0] = (int **) malloc(nz*ny*sizeof(int *));
  if (!f_array[0]) { 
    fprintf(stderr,"Error: (2) cannot allocate memory. Aborting. \n\n");
    exit(1);
  }
  f_array[0][0] = (int *) calloc(nz*ny*nx,sizeof(int));
  if (!f_array[0][0]) {
    fprintf(stderr,"Error: (3) cannot allocate memory. Aborting. \n\n");
    exit(1);
  }

  /* Calculate vector pointers */
  for(int j=1;j<ny;j++) {
    f_array[0][j] = f_array[0][j-1]+nx;
  }
  for(int k=1;k<nz;k++) {
    f_array[k]=f_array[k-1]+ny;
    f_array[k][0]=f_array[k-1][0]+ny*nx;
    for(int j=1;j<ny;j++) {
      f_array[k][j]=f_array[k][j-1]+nx;
    }
  }

  fprintf(stdout,"Finished setting up input 3d array.\n");

  /* Allocate tensor (i.e. 3d array) and set pointers */
  fprintf(stdout,"Setting up 3d array... \n");
  o_array = (float ***) malloc(nz*sizeof(float **));
  if (!o_array) {
    fprintf(stderr,"Error: (1) cannot allocate memory. Aborting. \n\n");
    exit(1);
  }
  o_array[0] = (float **) malloc(nz*ny*sizeof(float *));
  if (!o_array[0]) { 
    fprintf(stderr,"Error: (2) cannot allocate memory. Aborting. \n\n");
    exit(1);
  }
  o_array[0][0] = (float *) calloc(nz*ny*nx,sizeof(float));
  if (!o_array[0][0]) {
    fprintf(stderr,"Error: (3) cannot allocate memory. Aborting. \n\n");
    exit(1);
  }

  /* Calculate vector pointers */
  for(int j=1;j<ny;j++) {
    o_array[0][j] = o_array[0][j-1]+nx;
  }
  for(int k=1;k<nz;k++) {
    o_array[k]=o_array[k-1]+ny;
    o_array[k][0]=o_array[k-1][0]+ny*nx;
    for(int j=1;j<ny;j++) {
      o_array[k][j]=o_array[k][j-1]+nx;
    }
  }

  fprintf(stdout,"Finished setting up output 3d array.\n");

  dims[0] = nz;
  dims[1] = ny;
  dims[2] = nx;
  ext_space = H5Screate_simple(rank,dims, NULL);

  status = H5Sselect_hyperslab(ext_space, H5S_SELECT_SET, offset, NULL, count, dims);
  status = H5Dread(dset_org, H5T_NATIVE_INT, ext_space, org_space, H5P_DEFAULT, f_array[0][0]);

  fprintf(stdout,"Reading complete -- proceeding to writing stage...\n");

  for (int i=0;i<nx;i++) {
    for (int j=0;j<ny;j++) {
      for (int k=0;k<nz;k++) {
	if (reverse == 1) {
	  if ( f_array[i][j][k] == 0 ) { 
	    o_array[i][j][k] = 0;
	  }
	  else {
	    if ( f_array[i][j][k] == cluster_id ) {
	      o_array[i][j][k] = 0;
	    }
	    else {
	      o_array[i][j][k] = 1;
	    }
	  }	      
	}
	else {
	  if ( f_array[i][j][k] == cluster_id ) {
	    o_array[i][j][k] = 1;
	  }
	  else {
	    o_array[i][j][k] = 0;
	  }
	}
      }
    }
  }

  fprintf(stdout,"[HDF5] Writing HDF5 file '%s'...\n", fname_out);
  file_out = H5Fcreate(fname_out, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  dcpl = H5Pcreate(H5P_DATASET_CREATE);
  // fprintf(stdout,"Setting chunks to %d x %d x %d\n", (int) chunk[0], (int) chunk[1], (int) chunk[2]);
  // status = H5Pset_chunk(dcpl, 3, chunk);
  dset_ext = H5Dcreate(file_out, dset_name, H5T_NATIVE_FLOAT, ext_space, H5P_DEFAULT, dcpl, H5P_DEFAULT);
  status = H5Sselect_hyperslab(ext_space, H5S_SELECT_SET, offset, stride, count, dims);
  status = H5Dwrite(dset_ext, H5T_NATIVE_FLOAT, H5S_ALL, ext_space, H5P_DEFAULT, o_array[0][0]);
  fprintf(stdout,"[HDF5] Finished writing HDF5 file (status = %d).\n", status);
  fprintf(stdout,"Running h5repack -f SHUF -f GZIP=9 %s <repacked file> is recommended.\n", fname_out);

  return 0;
}
