#include "hdf5_helper.h"

const int HDF_RANK = 3;
const char *HDF_DEFAULT_DATASET = "OutArray";

VerbosityLevel HDF_VERBOSITY = NONE;

int LB3dAttributes::getString(char* key, string* value) {
  map<string,string>::iterator it;

  it = values.find(key);
  if ( it == values.end() ) {
    return error;
  }
  // Strip first and last characters from the string
  // Since the whitespace has been stripped, they should be the quotes
  *value = values[key].substr( 1, values[key].length()-2 );
  return ok;
}

int LB3dAttributes::getBool(char* key, bool* value) {
  map<string,string>::iterator it;

  it = values.find(key);
  if ( it == values.end() ) {
    return error;
  }

  if (values[key]==".false.") {
    *value = false;
  }
  else if (values[key]==".true.") {
    *value = true;
  }
  else {
    return error;
  }
  return ok;
}

int LB3dAttributes::getDouble(char* key, double* value) {
  map<string,string>::iterator it;

  it = values.find(key);
  if ( it == values.end() ) {
    return error;
  }
  *value = atof(values[key].c_str());
  return ok;
}

int LB3dAttributes::getInteger(char* key, int* value) {
  map<string,string>::iterator it;

  it = values.find(key);
  if ( it == values.end() ) {
    return error;
  }
  *value = atoi(values[key].c_str());
  return ok;
}

void LB3dAttributes::trimSpaces(string& str) {
    // Trim Both leading and trailing spaces
    size_t startpos = str.find_first_not_of(" \t"); // Find the first character position after excluding leading blank spaces
    size_t endpos = str.find_last_not_of(" \t"); // Find the first character position from reverse af
 
    // if all spaces or empty return an empty string
    if(( string::npos == startpos ) || ( string::npos == endpos))
    {
        str = "";
    }
    else
        str = str.substr( startpos, endpos-startpos+1 );
}

void LB3dAttributes::parseLines() {
  size_t eqpos;
  string key, value, line;

  values.clear();
  // Start from the bottom so one hits differential input files first
  for (int i = lines.size() - 1; i >= 0; i--) {
    line = lines.at(i);
    eqpos = line.find("=");
    //fprintf(stdout,"%d ",i);
    if (eqpos!=string::npos) {
      key = line.substr(0,eqpos);
      LB3dAttributes::trimSpaces(key);
      value = line.substr(eqpos+1,line.length());
      LB3dAttributes::trimSpaces(value);
      //fprintf(stdout,"%d [%s] [%s]\n",i,key.c_str(),value.c_str());
      values.insert(pair<string,string>(key,value));
    }
  }

//     // showing contents:
//   cout << "values contains:\n";
//   for ( values_it=values.begin() ; values_it != values.end(); values_it++ )
//     cout << (*values_it).first << " => " << (*values_it).second << endl;

}

LB3dAttributes::LB3dAttributes(const char* filename) {
  //const int attrname_len = 256;
  //char      attrname[attrname_len];
  hid_t     file_id, dataset_id, attr_id, dspace_id;
  hid_t     atype;
  hsize_t   dims[1],sdim;
  char**    rdata;
  //string*   str;

  file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  dataset_id = H5Dopen(file_id,HDF_DEFAULT_DATASET, H5P_DEFAULT);
  //fprintf(stdout,"Dataset has %d attributes\n",H5Aget_num_attrs(dataset_id));
  attr_id = H5Aopen_by_idx(dataset_id, ".", H5_INDEX_NAME ,H5_ITER_NATIVE, 0, H5P_DEFAULT, H5P_DEFAULT);
  atype = H5Aget_type (attr_id);
  //H5Aget_name(attr_id,attrname_len,attrname);
  //fprintf(stdout,"Attribute name is '%s'\n",attrname);
  //attr_id = H5Aopen(dataset_id,attrname,H5P_DEFAULT);

  sdim = H5Tget_size (atype);
  sdim++;
  dspace_id = H5Aget_space (attr_id);
  H5Sget_simple_extent_dims (dspace_id, dims, NULL);
    /*
   * Allocate array of pointers to rows.
   */
  rdata = (char **) malloc (dims[0] * sizeof (char *));

  /*
   * Allocate space for integer data.
   */
  rdata[0] = (char *) malloc (dims[0] * sdim * sizeof (char));
   
  /*
   * Set the rest of the pointers to rows to the correct addresses.
   */
  for (int i=1; i<dims[0]; i++)
    rdata[i] = rdata[0] + i * sdim;

  /*
   * Create the memory datatype.
   */
  atype = H5Tcopy (H5T_C_S1);
  H5Tset_size (atype, sdim);

  /*
   * Read the data.
   */
  H5Aread (attr_id, atype, rdata[0]);
    /*
   * Close and release resources.
   */

  for (int i=0; i<dims[0]; i++) {
    //fprintf(stdout,"%s\n",rdata[i]);
    string str(rdata[i]);
    lines.push_back(str);
  }

  free (rdata[0]);
  free (rdata);

  H5Sclose(dspace_id);
  H5Aclose(attr_id);
  H5Dclose(dataset_id);
  H5Fclose(file_id);
  
  LB3dAttributes::parseLines();
}

void hdfToIntDims(hsize_t dims[], int *nx, int *ny, int *nz) {
  if (HDF_VERBOSITY > 0) fprintf (stdout,"Dataspace dimensions are %d x %d x %d \n", (int) dims[2], (int) dims[1], (int) dims[0]);
  *nx = (int) dims[2];
  *ny = (int) dims[1];
  *nz = (int) dims[0];
}

void hdfToIntDims(hsize_t dims[], int *nx, int *ny, int *nz, int *nv) {
  if (HDF_VERBOSITY > 0) fprintf (stdout,"Dataspace dimensions are %d x %d x %d x %d\n", (int) dims[3], (int) dims[2], (int) dims[1], (int) dims[0]);
  *nv = (int) dims[3];
  *nx = (int) dims[2];
  *ny = (int) dims[1];
  *nz = (int) dims[0];
}


int hdfGetDims(char* filename, hsize_t **pdims) {
  return hdfGetDims(filename, pdims, HDF_DEFAULT_DATASET);
}

int hdfGetDims(char* filename, hsize_t **dims, const char *dsetname) {
  hid_t     file_src, dspace_src, dset_src;

  *dims = (hsize_t*) malloc(HDF_RANK*sizeof(hsize_t));
  if (!dims) {
    fprintf(stderr,"ERROR: Can't allocate memory for dims, aborting... \n");
    return 1;
  }

  file_src = H5Fopen (filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_src < 0) {
    fprintf(stderr,"ERROR: %d - Invalid file %s, aborting \n",file_src,filename);
    return 1;
  }
  dset_src = H5Dopen(file_src, dsetname, H5P_DEFAULT);
  dspace_src = H5Dget_space(dset_src);
  H5Sget_simple_extent_dims(dspace_src, *dims, NULL);

  H5Dclose (dset_src);
  H5Sclose (dspace_src);
  H5Fclose (file_src);
  return 0;
}

int hdfReadAll(char* filename, float ****pdata, hsize_t **pdims) {
  return hdfReadAll(filename, pdata, pdims, HDF_DEFAULT_DATASET);
}

int hdfReadAll(char* filename, float ****pdata, hsize_t **pdims, const char *dsetname) {
  //HDF variable declarations
  hid_t           file_src, dspace_src, dspace_tar, dset_src;
  const hsize_t   offset[HDF_RANK] = {0, 0, 0},
                  count[HDF_RANK]  = {1, 1, 1};
  int             nx, ny, nz;
  float           ***data;
  hsize_t         *dims;

  dims = (hsize_t*) malloc(HDF_RANK*sizeof(hsize_t));
  if (!dims) {
    fprintf(stderr,"ERROR: Can't allocate memory for dims, aborting... \n");
    return 1;
  }

  file_src = H5Fopen (filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_src < 0) {
    fprintf(stderr,"ERROR: %d - Invalid file %s, aborting \n",file_src,filename);
    return 1;
  }

  dset_src = H5Dopen(file_src, dsetname, H5P_DEFAULT);
  dspace_src = H5Dget_space(dset_src);
  H5Sget_simple_extent_dims(dspace_src,dims,NULL);
  nx = (int) dims[2];
  ny = (int) dims[1];
  nz = (int) dims[0];

  //fprintf (stdout,"Dataspace dimensions are %d x %d x %d \n", nx, ny, nz);

  /* Allocate tensor (i.e. 3d array) and set pointers */
  //fprintf(stdout,"Setting up 3d array for file %s... \n",filename);
  data = (float ***) malloc(nz*sizeof(float **));
  if (!data) {
    fprintf(stderr,"ERROR: (1) can't allocate memory, aborting... \n");
    return 1;
  }
  data[0] = (float **) malloc(nz*ny*sizeof(float *));
  if (!data[0]) { 
    fprintf(stderr,"ERROR: (2) can't allocate memory, aborting... \n");
    return 1;
  }
  data[0][0] = (float *) calloc(nz*ny*nx,sizeof(float));
  if (!data[0][0]) {
    fprintf(stderr,"ERROR: (3) can't allocate memory, aborting... \n");
    return 1;
  }

  /* Calculate vector pointers */
  for(int j=1;j<ny;j++) {
    data[0][j] = data[0][j-1]+nx;
  }
  for(int k=1;k<nz;k++) {
    data[k]=data[k-1]+ny;
    data[k][0]=data[k-1][0]+ny*nx;
    for(int j=1;j<ny;j++) {
      data[k][j]=data[k][j-1]+nx;
    }
  }

  dspace_tar = H5Screate_simple(HDF_RANK, dims, NULL);
  H5Sselect_hyperslab(dspace_tar, H5S_SELECT_SET, offset, NULL, count, dims);
  H5Dread(dset_src, H5T_NATIVE_FLOAT, dspace_tar, dspace_src, H5P_DEFAULT, &data[0][0][0]);

  *pdata = data;
  *pdims = dims;

  H5Dclose (dset_src);
  H5Sclose (dspace_src);
  H5Sclose (dspace_tar);
  H5Fclose (file_src);

  return 0;
}

int hdfReadSlab(char* filename, float ****pdata, const int ox, const int oy, const int oz, const int dx, const int dy, const int dz) {
  return hdfReadSlab(filename, pdata, ox, oy, oz, dx, dy, dz, HDF_DEFAULT_DATASET );
}

int hdfReadSlab(char* filename, float ****pdata, const int ox, const int oy, const int oz, const int dx, const int dy, const int dz, const char *dsetname) {
  //HDF variable declarations
  hid_t           file_src, dspace_src, dspace_tar, dset_src;    // Handles
  const hsize_t   offset[HDF_RANK] = {oz, oy, ox},
                  count[HDF_RANK]  = {1, 1, 1};
  float           ***data;
  hsize_t         *dims;

  /* Allocate tensor (i.e. 3d array) and set pointers */
  //fprintf(stdout,"Setting up 3d array for file %s... \n",filename);
  data = (float ***) malloc(dz*sizeof(float **));
  if (!data) {
    fprintf(stderr,"ERROR: (1) can't allocate memory, aborting... \n");
    return 1;
  }
  data[0] = (float **) malloc(dz*dy*sizeof(float *));
  if (!data[0]) { 
    fprintf(stderr,"ERROR: (2) can't allocate memory, aborting... \n");
    return 1;
  }
  data[0][0] = (float *) calloc(dz*dy*dx,sizeof(float));
  if (!data[0][0]) {
    fprintf(stderr,"ERROR: (3) can't allocate memory, aborting... \n");
    return 1;
  }

  /* Calculate vector pointers */
  for(int j=1;j<dy;j++) {
    data[0][j] = data[0][j-1]+dx;
  }
  for(int k=1;k<dz;k++) {
    data[k]=data[k-1]+dy;
    data[k][0]=data[k-1][0]+dy*dx;
    for(int j=1;j<dy;j++) {
      data[k][j]=data[k][j-1]+dx;
    }
  }

  // Attempt to allocate memory for dataslab dimensions
  dims = (hsize_t*) malloc(HDF_RANK*sizeof(hsize_t));
  if (!dims) {
    fprintf(stderr,"ERROR: Can't allocate memory for dims, aborting... \n");
    return 1;
  }
  dims[0] = (hsize_t) dz;
  dims[1] = (hsize_t) dy;
  dims[2] = (hsize_t) dx;

  // Attempt to open file
  file_src = H5Fopen (filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_src < 0) {
    fprintf(stderr,"ERROR: %d  - Invalid file %s, aborting \n",file_src,filename);
    return 1;
  }
  dset_src = H5Dopen(file_src, dsetname, H5P_DEFAULT);
  dspace_src = H5Dget_space(dset_src);

  H5Sselect_hyperslab(dspace_src, H5S_SELECT_SET, offset, NULL, count, dims);
  dspace_tar = H5Screate_simple(HDF_RANK, dims, NULL);
  H5Dread(dset_src, H5T_NATIVE_FLOAT, dspace_tar, dspace_src, H5P_DEFAULT, &data[0][0][0]);

  *pdata = data;

  H5Dclose (dset_src);
  H5Sclose (dspace_src);
  H5Sclose (dspace_tar);
  H5Fclose (file_src);

  free(dims);

  return 0;
}

int hdfReadAllVecComp(char* filename, float **pdata, hsize_t **pdims, const int voffset) {
  return hdfReadAllVecComp(filename, pdata, pdims, voffset, HDF_DEFAULT_DATASET);
}

int hdfReadAllVecComp(char* filename, float **pdata, hsize_t **pdims, const int voffset, const char *dsetname) {
  //HDF variable declarations
  hid_t           file_src, dspace_src, dspace_tar, dset_src;
  const hsize_t   offset[4] = {0, 0, 0, voffset};
  const hsize_t   count[4]  = {1, 1, 1, 1};
  int             nx, ny, nz;
  float           *data;
  hsize_t         *dims;

  dims = (hsize_t*) malloc(4*sizeof(hsize_t));
  if (!dims) {
    fprintf(stderr,"ERROR: Can't allocate memory for dims, aborting... \n");
    return 1;
  }

  file_src = H5Fopen (filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_src < 0) {
    fprintf(stderr,"ERROR: %d - Invalid file %s, aborting \n",file_src,filename);
    return 1;
  }

  dset_src = H5Dopen(file_src, dsetname, H5P_DEFAULT);
  dspace_src = H5Dget_space(dset_src);
  H5Sget_simple_extent_dims(dspace_src,dims,NULL);
  nx = (int) dims[2];
  ny = (int) dims[1];
  nz = (int) dims[0];

  dims[3] = (hsize_t) 1;

  // fprintf (stdout,"Dataspace dimensions are %d x %d x %d \n", nx, ny, nz);

  // Allocate tensor (i.e. 3d array) and set pointers
  //fprintf(stdout,"Setting up 3d array for file %s... \n",filename);
  data = (float *) malloc(nz*ny*nx*sizeof(float ));
  if (!data) {
    fprintf(stderr,"ERROR: (1) can't allocate memory, aborting... \n");
    return 1;
  }

  H5Sselect_hyperslab(dspace_src, H5S_SELECT_SET, offset, NULL, count, dims);
  dspace_tar = H5Screate_simple(HDF_RANK, dims, NULL);
  H5Dread(dset_src, H5T_NATIVE_FLOAT, dspace_tar, dspace_src, H5P_DEFAULT, &data[0]);

  *pdata = data;
  *pdims = dims;

  H5Dclose (dset_src);
  H5Sclose (dspace_src);
  H5Sclose (dspace_tar);
  H5Fclose (file_src);

  return 0;
}

int hdfFree3DArray(float ***data) {
  free(data[0][0]);
  free(data[0]);
  free(data);
  return 0;
}

int hdfFreeDims(hsize_t *dims) {
  free(dims);
  return 0;
}

int hdfFree4DArray(float ****data) {
  free(data[0][0][0]);
  free(data[0][0]);
  free(data[0]);
  free(data);
  return 0;
}
