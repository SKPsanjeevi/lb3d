/*includes*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "hdfread.h"

/* declarations */
int getBase(FLOATNUM ***colour,
	    int offset,
	    FLOATNUM *radius);

/* globals */
int		nx, ny, nz;
FLOATNUM	pi = 3.1415926535897932384626433832795028841971693993751058209749445923;
/* definitions */

int main(int argc,
	 char *argv[])
{
	int		retVal = 1; // 0 on success

	FLOATNUM	***colour;
	char		*filename;

	int		offset;

	FLOATNUM	radius;

	filename = argv[1];

	if (filename) {
		getDims(filename, &nx, &ny, &nz);
		hdfread(filename, "colour", &colour);
		retVal = 0;
	}

	getBase(colour,3,&radius);

	printf("%f\n",radius);

}

int getBase(FLOATNUM ***colour,
	    int offset,
	    FLOATNUM *radius)
{

	int		i,j,k;
	FLOATNUM	deltaColour;
	FLOATNUM	virtPos;
	FLOATNUM	lowY, highY;


		
	for (j = 1; j < ny-1; j++) {

		deltaColour = colour[offset][j][k] - colour[offset][j-1][k];

		if ( ((colour[offset][j][k] <= 0.) && (colour[offset][j+1][k] > 0.) )) {

			virtPos = j - (colour[offset][j][k] ) / ( colour[offset][j+1][k] - colour[offset][j][k] );
			lowY = virtPos;
		}
		     if ( ((colour[offset][j][k] >= 0.) && (colour[offset][j+1][k] < 0.) )) {
			virtPos = (j + 1) + ( colour[offset][j+1][k] ) / ( colour[offset][j][k] - colour[offset][j+1][k]);
			highY = virtPos;
		}
	}
	*radius = highY - lowY;
	
	return(0);
}
