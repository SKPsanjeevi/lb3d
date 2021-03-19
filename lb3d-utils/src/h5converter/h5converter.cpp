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
  FILE  *f;
  XDR   xdrs_in, xdrs_out;
  DataTypes filetype_in = UNDEFINED;
  DataTypes filetype_out = UNDEFINED;
  int nx = 0;
  int ny = 0;
  int nz = 0;
  float xdr_f;

  // File names
  char *fname_in, *fname_out;
  char *ftype_in, *ftype_out;
  char *dset_name;

  bool fname_in_set = false;
  bool fname_out_set = false;

  float ***f_array;

  //HDF variable declarations
  hid_t           file, file_out, org_space, ext_space, dset_org, dset_ext, dcpl;    // Handles
  herr_t          status;
  H5D_layout_t    layout;
  hsize_t         dims[rank],
                  offset[rank] = {0, 0, 0},
                  stride[rank] = {1, 1, 1},
                  count[rank] = {1, 1, 1},
                  chunk[rank] = {1, 1, 1}; // This is a default value and might be overwritten later

  fprintf(stdout,"\nH5 converter -- converts from HDF5, XDRF, RAW to HDF5, XDRF. \n\n");

  dset_name = "OutArray";

  for (int i = 1; i < argc ; i++) {

    // fprintf(stdout,"  argc = %d, argv[%d] = '%s'. \n", argc, i, argv[i]);
    if ( strcmp(argv[i],"-?") == 0 ) {
      usage();
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
      if ( i+2 < argc ) {
	fname_in = argv[++i];
	ftype_in = argv[++i];
	//fprintf(stdout,"  Input file name '%s', type '%s'. \n", fname_in, ftype_in);
	fname_in_set = true;
      }
      else {
	fprintf(stderr,"  Missing arguments to -i <filename> <filetype>. \n");
	exit(1);
      }
    }
    else if ( strcmp(argv[i],"-o") == 0 ) {
      if ( i+2 < argc ) {
	fname_out = argv[++i];
	ftype_out = argv[++i];
	//fprintf(stdout,"  Output file name '%s', type '%s'. \n", fname_out, ftype_out);
	fname_out_set = true;
      }
      else {
	fprintf(stderr,"  Missing argument to -o <filename> <filetype>. \n");
	exit(1);
      }
    }
    else if ( strcmp(argv[i],"-s") == 0 ) {
      if ( i+3 < argc ) {
	nx = atoi(argv[++i]);
	ny = atoi(argv[++i]);
	nz = atoi(argv[++i]);
	//fprintf(stdout,"  NX '%d'. \n", nx);
	//fprintf(stdout,"  NY '%d'. \n", ny);
	//fprintf(stdout,"  NZ '%d'. \n", nz);
      }
      else {
	fprintf(stderr,"  Missing arguments to -s <nx> <ny> <nz> flag. \n");
	exit(1);
      }
    }
    else if ( strcmp(argv[i],"-x") == 0 ) {
      if ( i+1 < argc ) {
	nx = atoi(argv[++i]);
	//fprintf(stdout,"  NX '%d'. \n", nx);
      }
      else {
	fprintf(stderr,"  Missing argument to -x <nx> flag. \n");
	exit(1);
      }
    }
    else if ( strcmp(argv[i],"-y") == 0 ) {
      if ( i+1 < argc ) {
	ny = atoi(argv[++i]);
	//fprintf(stdout,"  NY '%d'. \n", ny);
      }
      else {
	fprintf(stderr,"  Missing argument to -y <ny> flag. \n");
	exit(1);
      }
    }
    else if ( strcmp(argv[i],"-z") == 0 ) {
      if ( i+1 < argc ) {
	nz = atoi(argv[++i]);
	//fprintf(stdout,"  NZ '%d'. \n", nz);
      }
      else {
	fprintf(stderr,"  Missing argument to -z <nz> flag. \n");
	exit(1);
      }
    }

    else {
      //fprintf(stdout,"  No matches found for flag '%s'. \n", argv[i]);
    }
  }

  if ( !fname_out_set || !fname_in_set ) {
    //fprintf(stderr,"No input or output file specified. Aborting.\n");
    usage();
  }

  // Setting input file type
  //fprintf(stdout,"Setting ftype_in ...\n");
  if ( strcmp(ftype_in, "xdr" ) == 0 ) {
    filetype_in = XDRF;
  }
  else if ( strcmp(ftype_in, "hdf" ) == 0 ) {
    filetype_in = HDF5;
  }
  else if ( strcmp(ftype_in, "raw" ) == 0 ) {
    filetype_in = UINT8;
  }
  else {
    printf("Don't know how to classify inputfile based on '%s'. Aborting. \n\n", ftype_in);
    exit(1);
  }

  // Setting output file type
  //fprintf(stdout,"Setting ftype_out ...\n");
  if ( strcmp(ftype_out, "xdr" ) == 0 ) {
    filetype_out = XDRF;
  }
  else if ( strcmp(ftype_out, "hdf" ) == 0 ) {
    filetype_out = HDF5;
  }
  else {
    printf("Don't know how to classify outputfile based on '%s'. Aborting. \n\n", ftype_out);
    exit(1);
  }

  // Setting name of dataset
  //fprintf(stdout,"Set dset_name to '%s'.\n",dset_name);

  // Verifying system size
  if ( filetype_in == XDRF || filetype_in == UINT8 ) {
    if (nx <= 0 || ny <= 0 || nz <= 0) {
      fprintf(stderr,"[XDRF] nx, ny, nz need to be > 0. Aborting.\n\n");
      exit(1);
    }
  }
  
  // HDF prep
  if ( filetype_in == HDF5 ) {
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
  }

  /* Allocate tensor (i.e. 3d array) and set pointers */
  fprintf(stdout,"Setting up 3d array... \n");
  f_array = (float ***) malloc(nz*sizeof(float **));
  if (!f_array) {
    fprintf(stderr,"Error: (1) cannot allocate memory. Aborting. \n\n");
    exit(1);
  }
  f_array[0] = (float **) malloc(nz*ny*sizeof(float *));
  if (!f_array[0]) { 
    fprintf(stderr,"Error: (2) cannot allocate memory. Aborting. \n\n");
    exit(1);
  }
  f_array[0][0] = (float *) calloc(nz*ny*nx,sizeof(float));
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

  fprintf(stdout,"Finished setting up 3d array.\n");

  // Start reading
  if (filetype_in == UINT8) {
    uint8_t *ibuf;
    if (NULL==(f=fopen(fname_in,"rb"))) {
      fprintf(stderr,"[RAW] Error opening file '%s'. Aborting.\n\n",fname_in);
      exit(1);
    }
    fprintf(stdout,"[RAW] Reading binary integer file '%s'.\n", fname_in);
    if (NULL==(ibuf = (uint8_t*) malloc(sizeof(uint8_t)*nx*ny*nz))) {
      fprintf(stderr,"[RAW] Error allocating integer array. Aborting. \n\n");
      exit(1);
    }
    fread(ibuf,1,nx*ny*nz,f);
    fprintf(stdout,"[RAW] Finished reading binary integer file.\n");

    fprintf(stdout,"[RAW] Copying data from 1D to 3D array...\n");
    for (int k = 1; k <= nz; k++) {
      for (int j = 1; j <= ny; j++) {
        for (int i = 1; i <= nx; i++) {
	  f_array[k-1][j-1][i-1] = ibuf[(k-1)*nx*ny + (j-1)*nx + (i-1)];
	}
      }
    }
    fprintf(stdout,"[RAW] Finished copying data to 3D array.\n");
    free(ibuf);
    fclose(f);
  }

  if (filetype_in == XDRF ) {
    if (NULL==(f=fopen(fname_in,"r"))) {
      fprintf(stderr,"[XDRF] Error opening file '%s'. Aborting. \n\n",fname_in);
      exit(1);
    }

    fprintf(stdout,"[XDRF]  Reading XDRF file '%s'.\n", fname_out);

    xdrstdio_create(&xdrs_in,f,XDR_DECODE);

    for (int k = 1; k <= nz; k++) {
      for (int j = 1; j <= ny; j++) {
        for (int i = 1; i <= nx; i++) {
          if (1!=xdr_float(&xdrs_in,&xdr_f)) {
            printf("%f\n",xdr_f);
            perror("[XDRF] xdr_float()");
            exit(1);
          }
	  f_array[k-1][j-1][i-1] = xdr_f;
        }
      }
    }
    xdr_destroy(&xdrs_in);
    fprintf(stdout,"[XDRF] Finished reading XDRF file.\n");
    fclose(f);
  }

  if (filetype_in == HDF5 || filetype_out == HDF5) {
    dims[0] = nz;
    dims[1] = ny;
    dims[2] = nx;
    ext_space = H5Screate_simple(rank,dims, NULL);
  }

  if (filetype_in == HDF5) {
    status = H5Sselect_hyperslab(ext_space, H5S_SELECT_SET, offset, NULL, count, dims);
    status = H5Dread(dset_org, H5T_NATIVE_FLOAT, ext_space, org_space, H5P_DEFAULT, f_array[0][0]);
  }

  fprintf(stdout,"Reading complete -- proceeding to writing stage...\n");

  
  if (filetype_out == HDF5) {
    /*
    fprintf(stdout,"For HDF5 input file: 0 chunk size means keep original chunk size\n");
    fprintf(stdout,"Enter chunk_x: \n");
    scanf("%d",&chunk_x);
    if (chunk_x != 0)
      chunk[0] = chunk_x;
    fprintf(stdout,"Enter chunk_y: \n");
    scanf("%d",&chunk_y);
    if (chunk_y != 0)
      chunk[1] = chunk_y;
    fprintf(stdout,"Enter chunk_z: \n");
    scanf("%d",&chunk_z);
    if (chunk_z != 0)
      chunk[2] = chunk_z;
    */

    fprintf(stdout,"[HDF5] Writing HDF5 file '%s'...\n", fname_out);
    file_out = H5Fcreate(fname_out, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    dcpl = H5Pcreate(H5P_DATASET_CREATE);
    // fprintf(stdout,"Setting chunks to %d x %d x %d\n", (int) chunk[0], (int) chunk[1], (int) chunk[2]);
    // status = H5Pset_chunk(dcpl, 3, chunk);
    dset_ext = H5Dcreate(file_out, dset_name, H5T_NATIVE_FLOAT, ext_space, H5P_DEFAULT, dcpl, H5P_DEFAULT);
    status = H5Sselect_hyperslab(ext_space, H5S_SELECT_SET, offset, stride, count, dims);
    status = H5Dwrite(dset_ext, H5T_NATIVE_FLOAT, H5S_ALL, ext_space, H5P_DEFAULT, f_array[0][0]);
    fprintf(stdout,"[HDF5] Finished writing HDF5 file (status = %d).\n", status);
    fprintf(stdout,"Running h5repack -f SHUF -f GZIP=9 %s <repacked file> is recommended.\n", fname_out);
  }

  if (filetype_out == XDRF) {
    if (NULL==(f=fopen(fname_out,"w"))) {
      fprintf(stderr,"[XDRF] Error opening file '%s'. Aborting. \n\n",fname_out);
      exit(1);
    }

    fprintf(stdout,"[XDRF] Writing XDRF file '%s'.\n",fname_out);

    xdrstdio_create(&xdrs_out,f,XDR_ENCODE);

    for (int k = 1; k <= nz; k++) {
      for (int j = 1; j <= ny; j++) {
        for (int i = 1; i <= nx; i++) {
          xdr_f = f_array[k-1][j-1][i-1];
          if (1!=xdr_float(&xdrs_out,&xdr_f)) {
            printf("%f\n",xdr_f);
            perror("xdr_float()");
            exit(1);
          }
        }
      }
    }
    xdr_destroy(&xdrs_out);
    fprintf(stdout,"[XDRF] Finished writing XDRF file.\n");
    fclose(f);
  }

  return 0;
}
