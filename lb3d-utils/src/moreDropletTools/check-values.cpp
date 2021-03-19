#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <string.h>

#include "hdf5_helper.h"
#include "lbe_version.h"

using namespace std;

int main (int argc, char *argv[]) {
  // Dimensions
  int nx, ny, nz;
  int n = 0;

  int dlen;

  // File names
  char *fname;

  char fstr[32];

  // Array pointers
  float ***f;

  hsize_t *dims;
  double cutoff = 0.0;

  // Flags
  bool cutoff_set = false;
  bool fname_set = false;

  for (int i = 1; i < argc ; i++) {

    // fprintf(stdout,"  argc = %d, argv[%d] = <%s>\n",argc,i,argv[i]);

    if ( strcmp(argv[i],"-c") == 0 || strcmp(argv[i],"-cutoff") == 0 ) {
      if ( i+1 < argc ) {
	cutoff = atof(argv[++i]);
	//fprintf(stdout,"  Cutoff <%f>\n",cutoff);
	cutoff_set = true;
      }
      else {
	fprintf(stderr,"  Missing argument to -c flag. \n");
	exit(1);
      }
       }
    else if ( strcmp(argv[i],"-f") == 0 ) {
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

  if (! (fname_set && cutoff_set) ) {
    fprintf(stdout,"check-values (%s)\n\n", GIT_DESC);
    fprintf(stderr,"File -f and cutoff -c are mandatory.\n\n");
    fprintf(stderr,"  -c <cutoff>\n");
    fprintf(stderr,"  -f <filename>\n");
    fprintf(stderr,"\n");
    exit(1);
  }


  if (hdfReadAll(fname,&f,&dims) != 0) {
    exit(1);
  }
  hdfToIntDims(dims,&nx,&ny,&nz);
  hdfFreeDims(dims);
  fprintf(stdout,"Found %d x %d x %d dataset. Cutoff = %f .\n", nx, ny, nz, cutoff);

  dlen = ceil(log10(max(nx,max(ny,nz))));

  sprintf(fstr,"%%%d.%dd %%%d.%dd %%%d.%dd: %%+6.6e | ", dlen, dlen, dlen, dlen, dlen, dlen);

  //fprintf(stdout,"%s",fstr);

  for(int i=0;i<nx;i++) {
    for(int j=0;j<ny;j++) {
      for(int k=0;k<nz;k++) {
	if (fabs(f[k][j][i]) > cutoff) {
	  n++;
	  fprintf(stdout,fstr,i,j,k,f[k][j][i]);
	  for (int a=0;a<floor(logl(fabs(f[k][j][i])/cutoff)) && a < 30;a++) {
	    fprintf(stdout,"#");
	  }
	  fprintf(stdout,"\n");
	}
      }
    }
  }

  fprintf(stdout,"Found %d values > cutoff = %f .\n", n, cutoff);


  hdfFree3DArray(f);

  return 0;
}
