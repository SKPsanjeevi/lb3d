#include "write_bmp.h"

int write_bmp(char fname[], int nx, int ny, double **red, double **green, double **blue) {

  FILE *f;
  unsigned char *img = NULL;

  int x, y, r, g, b;

  int yres = ny;

  int filesize = 54 + 3*nx*ny;  //w is your image width, h is image height, both int

  img = (unsigned char *)malloc(3*nx*ny);
  memset(img,0,sizeof(img));

  for(int i=0; i<nx; i++) {
    for(int j=0; j<ny; j++) {
      x=i; y=(yres-1)-j;
      r = int(red[i][j]*255);
      g = int(green[i][j]*255);
      b = int(blue[i][j]*255);
      if (r > 255) r=255;
      if (g > 255) g=255;
      if (b > 255) b=255;
      img[(x+y*nx)*3+2] = (unsigned char)(r);
      img[(x+y*nx)*3+1] = (unsigned char)(g);
      img[(x+y*nx)*3+0] = (unsigned char)(b);
    }
  }

  unsigned char bmpfileheader[14] = {'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0};
  unsigned char bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0};
  unsigned char bmppad[3] = {0,0,0};

  bmpfileheader[ 2] = (unsigned char)(filesize    );
  bmpfileheader[ 3] = (unsigned char)(filesize>> 8);
  bmpfileheader[ 4] = (unsigned char)(filesize>>16);
  bmpfileheader[ 5] = (unsigned char)(filesize>>24);

  bmpinfoheader[ 4] = (unsigned char)(      nx    );
  bmpinfoheader[ 5] = (unsigned char)(      nx>> 8);
  bmpinfoheader[ 6] = (unsigned char)(      nx>>16);
  bmpinfoheader[ 7] = (unsigned char)(      nx>>24);
  bmpinfoheader[ 8] = (unsigned char)(      ny    );
  bmpinfoheader[ 9] = (unsigned char)(      ny>> 8);
  bmpinfoheader[10] = (unsigned char)(      ny>>16);
  bmpinfoheader[11] = (unsigned char)(      ny>>24);

  f = fopen(fname,"wb");
  fwrite(bmpfileheader,1,14,f);
  fwrite(bmpinfoheader,1,40,f);
  for(int j=0; j<ny; j++) {
    fwrite(img+(nx*(ny-j-1)*3),3,nx,f);
    fwrite(bmppad,1,(4-(nx*3)%4)%4,f);
  }
  fclose(f);

  free(img);
  return 0;
}

double** Make2DDoubleArray(int arraySizeX, int arraySizeY) {
  double** theArray;
  theArray = (double**) malloc(arraySizeX*sizeof(double*));
  for (int i = 0; i < arraySizeX; i++)
     theArray[i] = (double*) malloc(arraySizeY*sizeof(double));
   return theArray;
} 


void Free2DDoubleArray(double **a, int nx) {
  for (int i = 0; i < nx; i++){
   free(a[i]);
  }
  free(a);
}
