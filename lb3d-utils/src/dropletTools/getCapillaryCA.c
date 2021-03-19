/// getCapillaryCA.c 
/// Author: Sebastian Schmieschek
///
/// Geometrically determine the contact angle of a surface  
/// of two fluid components in a pseudo 2D channel
///
/// Expects a 3D hdf colour output file as argument
///
/// ATTENTION: hdfread returns array as array[z][y][x]
///
/// Finds height via sign change of colour field at x=2, x=nx/2 and x=nx-2
/// Expects base to be nx - 1
/// Calculates geometrical contact angle as 
/// theta = atan((base / 2), (radius - height))
/// Prints result to stdout (in radiants)

/*includes*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "hdfread.h"

/* declarations */
int getZeroColourZ(FLOATNUM ***colour,
		   int xpos,
		   int ypos,
		   FLOATNUM *z1,
		   FLOATNUM *z2,
		   int xdelta,
		   int ydelta);

int calcGeometricContactAngle(FLOATNUM base,
			      FLOATNUM height,
			      FLOATNUM *theta);
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
	FLOATNUM	***colour;
	char		*filename;
	char		*outfilename;

	int		direction;

	FLOATNUM	height1, height2;
	FLOATNUM	theta1,theta2;

	FLOATNUM	radius;


	// Read filename to array colour
	filename = argv[1];

	if (filename) {
		getDims(filename,&nx,&ny,&nz);
		hdfread2(filename,&colour);

		retVal = 0;
	}
	
	// Determine sign change of colour field
	getZeroColourZ(colour,2,ny/2,&zPos1r,&zPos2r,0,0);
	getZeroColourZ(colour,nz/2,ny/2,&zPos1c,&zPos2c,0,0);
	getZeroColourZ(colour,nz-2,ny/2,&zPos1l,&zPos2l,0,0);

	
	// Determine height as difference of left hand side to center and rhs to center
	height1 = zPos2l - zPos2c;
	height2 = zPos2r - zPos2c;

	//	printf("height1: %f\theight2: %f\n",height1, height2);

	// Calculate contact angle
	calcGeometricContactAngle(nz-1,height1,&theta1);
	calcGeometricContactAngle(nz-1,height2,&theta2);
	
	//	printf("theta1: %f\ttheta2: %f\n",theta1,theta2);
	//      printf("%f %f\n",theta1, theta1*(180/pi));
	printf("%f\n",(pi - theta2));
	free(colour);
	
	return(retVal);
}



/// Determine z position of sign change in colour field for given x and y position
/// z1 gives value for negative -> positive, z2 for positive -> negative transition
/// if xdelta and/or ydelta are not 0 mean values for an area +\- delta will be returned

int getZeroColourZ(FLOATNUM ***colour,
		   int xpos,
		   int ypos,
		   FLOATNUM *z1, 
		   FLOATNUM *z2,
		   int xdelta,
		   int ydelta)

{

	int		i,j,k;
	int		valCount1, valCount2;
	FLOATNUM	zPos1, zPos2;

	valCount1 = 0;
	valCount2 = 0;
	zPos1 = 0;
	zPos2 = 0;


	for (i=10;i<=nx-12;i++) {
		for (j=ypos-ydelta;j<=ypos+ydelta;j++) {
			for (k=xpos-xdelta;k<=xpos+xdelta;k++) {
				//printf("%f\t%d\n",colour[i][j][k],i); 
				if (colour[i][j][k] <= 0. && colour[i+1][j][k] > 0.) {
					zPos1 += i + (colour[i][j][k] ) / (colour[i+1][j][k] - colour[i][j][k]);
					valCount1++;
					//printf("bong\t%f\t%d\n",zPos1,valCount1);
				}
				if (colour[i][j][k] >= 0. && colour[i+1][j][k] < 0.) {
					zPos2 += (i + 1) + (colour[i+1][j][k] ) / ( colour[i][j][k] - colour[i+1][j][k]);
					valCount2++;
					//printf("bing\t%f\t%d\n",zPos2,valCount2);
				}
			}
		}	
	}
     
	*z1=zPos1/valCount1;
	*z2=zPos2/valCount2;

	return(0);
}


int calcGeometricContactAngle(FLOATNUM base,
			      FLOATNUM height,
			      FLOATNUM *theta)
{
	FLOATNUM	radius = 0.;

	radius = -((4*(pow(height,2)))+(pow(base,2)))/(8 * height);
	//		printf("base height radius: %f\n",radius);

	*theta = atan2((base / 2), (radius - height));// * (180 / pi);

	return(0);
}
