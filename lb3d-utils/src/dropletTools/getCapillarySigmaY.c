/*includes*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "hdfread.h"

/* declarations */
int getSigma(FLOATNUM ***pxx,
		   FLOATNUM ***pzz,
		   int xpos,
		   int ypos,
		   FLOATNUM *sigma,
		   int xdelta,
		   int ydelta);

/* globals */
int		nx, ny, nz;
FLOATNUM	pi = 3.1415926535897932384626433832795028841971693993751058209749445923;
/* definitions */

int main(int argc,
	 char *argv[])
{
	int		retVal = 1; // 0 on success

	int		xPos,yPos;
	int		xDelta, yDelta;

	FLOATNUM	zPos1r, zPos2r; // right
	FLOATNUM	zPos1c, zPos2c; // center 
	FLOATNUM	zPos1l, zPos2l; // left
	FLOATNUM	***pzz;
	FLOATNUM	***pyy;
	char		*filename;
	char		*outfilename;

	int		direction;

	FLOATNUM	height1, height2;
	FLOATNUM	sigma;

	FLOATNUM	radius;

	filename = argv[1];
	outfilename = argv[2];

	if (filename) {
		getDims(filename, &nx, &ny, &nz);
		hdfread(filename, "pzz", &pzz);
		hdfread(filename, "pyy", &pyy);
		
		//	hdfwrite(outfilename, nx, ny, nz, colour);

		retVal = 0;
	}

	getSigma(pzz,pyy,nx/2,nz/2,&sigma,0,0);

	printf("sigma: %f\n",sigma);

}


int getSigma(FLOATNUM ***pzz,
	      FLOATNUM ***pyy,
	      int xpos,
	      int ypos,
	      FLOATNUM *sigma, 
	      int xdelta,
	      int ydelta)

{

	int		i,j,k,l;
	FLOATNUM	deltaPlocal[ny];
	FLOATNUM	sigma1 = 0.;

	printf("nx: %d, ny: %d, nz: %d\n",nx,ny,nz);
	l=0;

	for (j=ypos-ydelta;j<=ypos+ydelta;j++) {
		for (k=xpos-xdelta;k<=xpos+xdelta;k++) {
			deltaPlocal[0] = 0.;
			deltaPlocal[0] = ( pyy[j][0][k] ) - ( pzz[j][0][k] ); 
		}
	}

	for (i=0;i<ny-1;i++) {
		printf("%d\n",i); 
		for (j=ypos-ydelta;j<=ypos+ydelta;j++) {
			for (k=xpos-xdelta;k<=xpos+xdelta;k++) {
				deltaPlocal[i+1] = 0.;
				deltaPlocal[i+1] = ( pyy[j][i+1][k] ) - ( pzz[j][i+1][k] );
				l++;
				printf("%d %d %e\n",i,l, deltaPlocal[i]); 
				sigma1 = sigma1 + deltaPlocal[i];
			}
		}	
	}

	*sigma = sigma1 / 2.0;
	
	
	return(0);
}

