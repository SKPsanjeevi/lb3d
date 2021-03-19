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

  int time = 0;

  // File names
  char *odfname;
  char *wdfname;

  char fstr[32];

  // Array pointers
  float ***od;
  float ***wd;

  double cutoff = 0.0;

  double oind, ooutd;
  double wind, woutd;

  hsize_t *dims;

  // Flags
  bool odfname_set = false;
  bool wdfname_set = false;

  for (int i = 1; i < argc ; i++) {

    // fprintf(stdout,"  argc = %d, argv[%d] = <%s>\n",argc,i,argv[i]);

    if ( strcmp(argv[i],"-c") == 0 ) {
      if ( i+1 < argc ) {
	cutoff = atof(argv[++i]);
	//fprintf(stdout,"  File <%s>\n",fname);
      }
      else {
	fprintf(stderr,"  Missing argument to -c flag. \n");
	exit(1);
      }
    }
    else if ( strcmp(argv[i],"-o") == 0 ) {
      if ( i+1 < argc ) {
	odfname = argv[++i];
	//fprintf(stdout,"  File <%s>\n",fname);
	odfname_set = true;
      }
      else {
	fprintf(stderr,"  Missing argument to -o flag. \n");
	exit(1);
      }
    }
    else if ( strcmp(argv[i],"-t") == 0 ) {
      if ( i+1 < argc ) {
	time = atoi(argv[++i]);
	//fprintf(stdout,"  File <%s>\n",fname);
      }
      else {
	fprintf(stderr,"  Missing argument to -t flag. \n");
	exit(1);
      }
    }
    else if ( strcmp(argv[i],"-w") == 0 ) {
      if ( i+1 < argc ) {
	wdfname = argv[++i];
	//fprintf(stdout,"  File <%s>\n",fname);
	wdfname_set = true;
      }
      else {
	fprintf(stderr,"  Missing argument to -w flag. \n");
	exit(1);
      }
    }
    else {
      //fprintf(stdout,"  No matches for <%s>.\n",argv[i]);
    }
  }

  if (! (odfname_set && wdfname_set) ) {
    fprintf(stdout,"check-minority (%s)\n\n", GIT_DESC);
    fprintf(stderr,"OD file -o and WD file -w are mandatory.\n\n");
    fprintf(stderr,"  -o <filename>\n");
    fprintf(stderr,"  -t <time>\n");
    fprintf(stderr,"  -w <filename>\n");
    fprintf(stderr,"\n");
    exit(1);
  }

  if (hdfReadAll(odfname,&od,&dims) != 0) {
    exit(1);
  }
  if (hdfReadAll(wdfname,&wd,&dims) != 0) {
    exit(1);
  }
  hdfToIntDims(dims,&nx,&ny,&nz);
  hdfFreeDims(dims);

  for(int i=0;i<nx;i++) {
    for(int j=0;j<ny;j++) {
      for(int k=0;k<nz;k++) {
	if (fabs(od[k][j][i]) > cutoff) {
	  n++;
	  fprintf(stdout,fstr,i,j,k,od[k][j][i]);
	  for (int a=0;a<floor(logl(fabs(od[k][j][i])/cutoff)) && a < 30;a++) {
	    fprintf(stdout,"#");
	  }
	  fprintf(stdout,"\n");
	}
      }
    }
  }

  fprintf(stdout,"%d, %E, %E, %E, %E\n",time,od[0][0][0],od[nz/2][ny/2][nx/2],wd[0][0][0],wd[nz/2][ny/2][nx/2]);

  hdfFree3DArray(od);
  hdfFree3DArray(wd);

  return 0;
}
