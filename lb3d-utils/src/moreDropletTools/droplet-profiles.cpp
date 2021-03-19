//#include <stdlib.h>
//#include <iostream>
#include <string.h>
#include <limits>
#include "hdf5_helper.h"
#include "lbe_version.h"
#include "droplet-profiles.h"

using namespace std;

int find_extentminz(float*** N,const int dx, const int dy, const int dz, const double cutoff) {
  for(int k=0;k<dz;k++) {
    for(int j=0;j<dy;j++) {
      for(int i=0;i<dx;i++) {
	if (N[k][j][i] > cutoff) {
	  return k;
	}
      }
    }
  }
  return -1;
}

int find_extentmaxz(float*** N,const int dx, const int dy, const int dz, const double cutoff) {
  for(int k=dz-1;k>=0;k--) {
    for(int j=0;j<dy;j++) {
      for(int i=0;i<dx;i++) {
	if (N[k][j][i] > cutoff) {
	  return k;
	}
      }
    }
  }
  return -1;
}


int find_topz(float*** N,const int dx, const int dy, const int dz, const double cutoff) {
  for(int i=dx-1;i>=0;i--) {
    for(int k=dz-1;k>=0;k--) {
      for(int j=0;j<dy;j++) {
	if (N[k][j][i] > cutoff) {
	  fprintf(stdout,"Bottom z: %d %d %d\n",i,j,k);
	  return k;
	}
      }
    }
  }
  return -1;
}

int find_botz(float*** N,const int dx, const int dy, const int dz, const double cutoff) {
  for(int i=0;i<dx;i++) {
    for(int k=0;k<dz;k++) {
      for(int j=0;j<dy;j++) {
	if (N[k][j][i] > cutoff) {
	  fprintf(stdout,"Bottom z: %d %d %d\n",i,j,k);
	  return k;
	}
      }
    }
  }
  return -1;
}

double gammadot(float***N,const int dx, const int j, const int dz, const int ox, const int oy, const int oz) {
  double gd;
  for(int k=0;k<dz;k++) {
    for(int i=1;i<dx;i++) {
      fprintf(stdout,"# %d %d %d %lf\n",i+ox,j+oy,k+oz,N[k][j][i]);
      if (N[k][j][i] < N[k][j][i-1]) {
	gd = (N[k][j][i] - N[k][j][0]) / (double) i;
	// fprintf(stdout,"Returning %lf\n",gd);
	return gd;
      }
    }
  }
  gd = (N[0][j][dx-1]-N[0][j][0]) / (double) (dx -1);
  // fprintf(stdout,"Returning %lf %lf %d %d %d %lf\n",N[0][j][dx-1],N[0][j][0],dx,j,dz,gd);
  return gd;
  //return numeric_limits<double>::quiet_NaN();;
}

int main (int argc, char *argv[]) {
  // Lattice dimensions
  int nx, ny, nz;
  // Offsets
  int ox, oy, oz;
  // Slice dimensions
  int dx, dy, dz;

  // Offsets
  int offset;

  // File names
  char *fname_od,*fname_vel;

  // Array pointers
  float ***N;

  double cutoff = 0.0;

  hsize_t *dims;

  CardinalDirection dir;

  // Flags
  bool fname_od_set = false;
  bool fname_vel_set = false;
  bool dir_set = false;
  // bool offset_set = false;

  bool header = false;
  // bool verbose = false;

  for (int i = 1; i < argc ; i++) {

    // fprintf(stdout,"  argc = %d, argv[%d] = <%s>\n",argc,i,argv[i]);

    if      ( strcmp(argv[i],"-c") == 0 ) {
      if ( i+1 < argc ) {
	cutoff = atof(argv[++i]);
	//fprintf(stdout,"  Cutoff <%f>\n",cutoff);
      }
      else {
	fprintf(stderr,"  Missing argument to -c flag. \n");
	exit(1);
      }
    }
    else if ( strcmp(argv[i],"-d") == 0 || strcmp(argv[i],"--direction") == 0   ) {
      if ( i+1 < argc ) {
	++i;
     	if      ( strcmp(argv[i],"x") == 0 || strcmp(argv[i],"X") == 0 ) {
	  dir_set = true;
	  dir = X;
	}
	else if ( strcmp(argv[i],"y") == 0 || strcmp(argv[i],"Y") == 0 ) {
	  dir_set = true;
	  dir = Y;
	}
	else if ( strcmp(argv[i],"z") == 0 || strcmp(argv[i],"Z") == 0 ) {
	  dir_set = true;
	  dir = Z;
	}
	else {
	  fprintf(stderr,"  Unknown argument to -d flag. \n");
	  exit(1);
	}
      }
      else {
	fprintf(stderr,"  Missing argument to -d flag. \n");
     	exit(1);
       }
    }
    else if ( strcmp(argv[i],"-f") == 0 ) {
      if ( i+1 < argc ) {
	fname_od = argv[++i];
	//fprintf(stdout,"  File <%s>\n",fname_od);
	fname_od_set = true;
      }
      else {
	fprintf(stderr,"  Missing argument to -f flag. \n");
	exit(1);
      }
    }
    else if ( strcmp(argv[i],"-h") == 0 || strcmp(argv[i],"--header") == 0    ) {
      header = true;
      //fprintf(stdout,"  Header\n");
    }
    else if ( strcmp(argv[i],"-o") == 0 || strcmp(argv[i],"--offset") == 0 ) {
      if ( i+1 < argc ) {
	// offset_set = true;
	offset = atoi(argv[++i]);
      }
      else {
	fprintf(stderr,"  Missing argument to -o flag. \n");
     	exit(1);
      }
    }
    else if ( strcmp(argv[i],"-v") == 0 ) {
      if ( i+1 < argc ) {
	fname_vel = argv[++i];
	//fprintf(stdout,"  File <%s>\n",fname_vel);
	fname_vel_set = true;
      }
      else {
	fprintf(stderr,"  Missing argument to -v flag. \n");
	exit(1);
      }
    }  }

  if (! (fname_od_set && fname_vel_set && dir_set )  ) {
    fprintf(stdout,"droplet-profiles (%s)\n\n", GIT_DESC);
    fprintf(stderr,"Flags -d -f -v are mandatory. \n\n");
    fprintf(stderr,"  -c <cutoff>\n");
    fprintf(stderr,"  -d <direction>\n");
    fprintf(stderr,"  -f <od file>\n");
    fprintf(stderr,"  -h \n");
    fprintf(stderr,"  -o <offset>\n");
    fprintf(stderr,"  -v <vel file>\n");
    fprintf(stderr,"\n");
    exit(1);
  }

  // Get lattice size
  if (hdfGetDims(fname_vel,&dims) != 0) {
    exit(1);
  }
  hdfToIntDims(dims,&nx,&ny,&nz);
  hdfFreeDims(dims);


  // Find droplet extent
  if (hdfReadAll(fname_od,&N,&dims) != 0) {
    exit(1);
  }
  hdfToIntDims(dims,&dx,&dy,&dz);
  hdfFreeDims(dims);

  int extentmin = find_extentminz(N,dx,dy,dz,cutoff);
  int extentmax = find_extentmaxz(N,dx,dy,dz,cutoff);

  fprintf(stdout,"# Drop extent: z-min = %d , z-max = %d\n",extentmin,extentmax);

  hdfFree3DArray(N);

  // Found extent

  if (header) {
    fprintf(stdout,"# z    calc     imp\n");
  }

  double gd_max = 0.0;
  double gd_avg_drop = 0.0;
  double gd_avg_all = 0.0;

  for (int o=0;o<nz;o++) {
 
    offset = o;

    switch(dir) {
    case X:
      ox = offset;
      oy = 0;
      oz = 0;
      dx = 1;
      dy = ny;
      dz = nz;
      break;
    case Y:
      ox = 0;
      oy = offset;
      oz = 0;
      dx = nx;
      dy = 1;
      dz = nz;
      break;
    case Z:
      ox = 0;
      oy = 0;
      oz = offset;
      dx = nx;
      dy = ny;
      dz = 1;
      break;
    }

    if (hdfReadSlab(fname_vel,&N,ox,oy,oz,dx,dy,dz) != 0) {
      exit(1);
    }

    int j = (dy/2)-1;

    double shear_u;

    // Get attributes
    LB3dAttributes attr = LB3dAttributes(fname_od);
    attr.getDouble("shear_u",&shear_u);
    //fprintf(stdout,"shear_u = %f\n",shear_u);

    double gd = gammadot(N,dx,j,dz,ox,oy,oz);
  
    if (gd > gd_max) gd_max = gd;
    if (o <= extentmax && o >= extentmin) gd_avg_drop += gd;
    gd_avg_all += gd;

    fprintf(stdout,"%d %lf %lf\n",o,gd,2.0*shear_u/(dx-1.0));

    hdfFree3DArray(N);

  }

  gd_avg_drop /= (double) (extentmax-extentmin+1);
  gd_avg_all /= nz;

  fprintf(stdout,"# gd_max = %lf, gd_avg_drop = %lf, gd_avg_all = %lf\n",gd_max,gd_avg_drop,gd_avg_all);

  return 0;
}


