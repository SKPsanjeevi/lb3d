#include "droplet_fit.h"
#include <vector>
#include <math.h>

#ifdef WRITEBMP
#include "write_bmp.h"
#endif

using namespace std;

EllipseFit getEllipseFitFromSlice(char* fname, const CardinalDirection dir, const int offset) {
  float ***c;
  vector< vector<double> > points;
  vector<double> pos;

  int nx, ny, nz;
  int dx, dy, dz;
  int ox, oy, oz;

#ifdef WRITEBMP
  int bmpx, bmpy;
  char* bmpname;
  double scale = 8.0;
#endif

  double x,y;

  hsize_t *dims;

  if (hdfGetDims(fname,&dims) != 0) {
    fprintf(stderr,"hdfGetDims returned error...\n");
    exit(1);
  }
  hdfToIntDims(dims,&nx,&ny,&nz);
  hdfFreeDims(dims);

  switch (dir) {

    case X:
      dx = 1;
      dy = ny;
      dz = nz;
#ifdef WRITEBMP
      bmpx = int(scale*ny);
      bmpy = int(scale*nz);
      bmpname = "testX.bmp";
#endif

      if ( offset < 0 ) {
        ox = nx/2;
      }
      else {
        ox = offset;
      }
      oy = 0;
      oz = 0;
      break;

    case Y:
      dx = nx;
      dy = 1;
      dz = nz;
#ifdef WRITEBMP
      bmpx = int(scale*nz);
      bmpy = int(scale*nx);
      bmpname = "testY.bmp";
#endif

      ox = 0;
      if ( offset < 0 ) {
        oy = ny/2;
      }
      else {
        oy = offset;
      }
      oz = 0;
      break;

    case Z:
      dx = nx;
      dy = ny;
      dz = 1;
#ifdef WRITEBMP
      bmpx = int(scale*nx);
      bmpy = int(scale*ny);
      bmpname = "testZ.bmp";
#endif

      ox = 0;
      oy = 0;
      if ( offset < 0 ) {
        oz = nz/2;
      }
      else {
        oz = offset;
      }

      break;
  }

  //fprintf(stdout,"Dataset size: %d x %d x %d.\n",nx,ny,nz);
  //fprintf(stdout,"hdfReadSlab: %d x %d x %d at (%d, %d, %d).\n",dx,dy,dz,ox,oy,oz);

  if (hdfReadSlab(fname,&c,ox,oy,oz,dx,dy,dz) != 0) {
    fprintf(stderr,"hdfReadSlab returned error...\n");
    exit(1);
  }

  switch (dir) {

    case X:
      for (int k = 0; k < dz; k++)
        for (int j = 1; j < dy; j++)
          for (int i = 0; i < dx; i++)
            if ( c[k][j][i] * c[k][j-1][i] < 0 ) {
              x = (j-1) + ( c[k][j-1][i] / ( c[k][j-1][i] - c[k][j][i] ) );
              y = k;
              pos.push_back(x);
              pos.push_back(y);
              points.push_back(pos);
              pos.clear();
            }
      for (int k = 1; k < dz; k++)
        for (int j = 0; j < dy; j++)
          for (int i = 0; i < dx; i++)
            if ( c[k][j][i] * c[k-1][j][i] < 0 ) {
              x = j;
              y = (k-1) + ( c[k-1][j][i] / ( c[k-1][j][i] - c[k][j][i] ) );
              pos.push_back(x);
              pos.push_back(y);
              points.push_back(pos);
              pos.clear();
            }
      break;

    case Y:
      for (int k = 1; k < dz; k++)
        for (int j = 0; j < dy; j++)
          for (int i = 0; i < dx; i++)
            if ( c[k][j][i] * c[k-1][j][i] < 0 ) {
              x = (k-1) + ( c[k-1][j][i] / ( c[k-1][j][i] - c[k][j][i] ) );
              y = i;
              pos.push_back(x);
              pos.push_back(y);
              points.push_back(pos);
              pos.clear();
            }
      for (int k = 0; k < dz; k++)
        for (int j = 0; j < dy; j++)
          for (int i = 1; i < dx; i++)
            if ( c[k][j][i] * c[k][j][i-1] < 0 ) {
              x = k;
              y = (i-1) + ( c[k][j][i-1] / ( c[k][j][i-1] - c[k][j][i] ) );
              pos.push_back(x);
              pos.push_back(y);
              points.push_back(pos);
              pos.clear();
            }
      break;

    case Z:
      for (int k = 0; k < dz; k++)
        for (int j = 0; j < dy; j++)
          for (int i = 1; i < dx; i++)
            if ( c[k][j][i] * c[k][j][i-1] < 0 ) {
              x = (i-1) + ( c[k][j][i-1] / ( c[k][j][i-1] - c[k][j][i] ) );
              y = j;
              pos.push_back(x);
              pos.push_back(y);
              points.push_back(pos);
              pos.clear();
            }
      for (int k = 0; k < dz; k++)
        for (int j = 1; j < dy; j++)
          for (int i = 0; i < dx; i++)
            if ( c[k][j][i] * c[k][j-1][i] < 0 ) {
              x = i;
              y = (j-1) + ( c[k][j-1][i] / ( c[k][j-1][i] - c[k][j][i] ) );
              pos.push_back(x);
              pos.push_back(y);
              points.push_back(pos);
              pos.clear();
            }
      break;
  }

  hdfFree3DArray(c);

  EllipseFit fit(points, NONE);
  fit.doFit();

#ifdef WRITEBMP

  vector<double> centre = fit.getCentre();
  double rotation = fit.getRotation();
  
  double **red = Make2DDoubleArray(bmpx,bmpy);
  double **green = Make2DDoubleArray(bmpx,bmpy);
  double **blue = Make2DDoubleArray(bmpx,bmpy);

  for (int i=0; i<bmpx; i++) {
    for (int j=0; j<bmpy; j++) {
      red[i][j] = 1.0;
      green[i][j] = 1.0;
      blue[i][j] = 1.0;
    }
  }

  for(int i = 0; i < points.size(); i++) {
    red  [ int(scale*points[i][0]) ][ int(scale*points[i][1]) ] = 1.0;
    green[ int(scale*points[i][0]) ][ int(scale*points[i][1]) ] = 0.0;
    blue [ int(scale*points[i][0]) ][ int(scale*points[i][1]) ] = 0.0;
  }

  for( int i =0; i < bmpx ; i++) {
    double ypos = (i - scale*centre[0])*tan(rotation) + scale*centre[1];
    int nypos = int (max(min(ypos,1.0*bmpy),0.0));
    red  [i][ nypos ] = 0.0;
    green[i][ nypos ] = 0.0;
    blue [i][ nypos ] = 1.0;
  }

  red  [ int(scale*centre[0]) ][ int(scale*centre[1]) ] = 0.0;
  green[ int(scale*centre[0]) ][ int(scale*centre[1]) ] = 1.0;
  blue [ int(scale*centre[0]) ][ int(scale*centre[1]) ] = 0.0;

  write_bmp(bmpname,bmpx,bmpy,red,green,blue);

  Free2DDoubleArray(red,nx);
  Free2DDoubleArray(green,nx);
  Free2DDoubleArray(blue,nx);

#endif

  return fit;
}

EllipseFit getEllipseFitFromSlice(char* fname, const CardinalDirection dir) {
  return getEllipseFitFromSlice(fname, dir, -1);
}
