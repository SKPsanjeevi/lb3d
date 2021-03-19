/// getCapillaryCAnew.c 
/// Author: Qingguang Xie, Sebastian Schmieschek
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
   int      xstart, xend;
	FLOATNUM	zPos1r, zPos2r; // right
	FLOATNUM	zPos1c, zPos2c; // center 
	FLOATNUM	zPos1l, zPos2l; // left
	FLOATNUM	***colour;
	char		*filename;
	char		*outfilename;

	int		direction;

	FLOATNUM	base;
	FLOATNUM	theta;
   FLOATNUM height;
	FLOATNUM	radius;


	// Read filename to array colour
	filename = argv[1];

	if (filename) {
		getDims(filename,&nx,&ny,&nz);
		hdfread2(filename,&colour);

		retVal = 0;
	}
    printf("first point ?\n");
    scanf("%d",&xstart);	
    printf("second point ?\n");
    scanf("%d",&xend);
	// Determine sign change of colour field
	getZeroColourXstart(colour,xstart,ny/2,&zPos1l,&zPos2l,10,0);
   getZeroColourXend(colour,xend,ny/2,&zPos1r,&zPos2r,10,0);
	getZeroColourZ(colour,nz/2,ny/2,&zPos1c,&zPos2c,0,0);
    	
//getZeroColourZ(colour,nz-2,ny/2,&zPos1l,&zPos2l,0,0);

	
	// Determine height as difference of left hand side to center and rhs to center
	//height1 = zPos2l - zPos2c;
	//height2 = zPos2r - zPos2c;
   base  = zPos2r - zPos2l;
   //base2  = 
   height = zPos2c - 1;
   printf("height: %f\n",height);
   printf("zPos1l: %f\n",zPos1l);
   printf("zPos2l: %f\n",zPos2l);
   printf("zPos1r: %f\n",zPos1r);
   printf("zPos2r: %f\n",zPos2r);
   printf("zPos1c: %f\n",zPos1c);
   printf("zPos2c: %f\n",zPos2c);
   
   printf("base: %f\n", base);
	// Calculate contact angle
	calcGeometricContactAngle(base,height,&theta);
	//calcGeometricContactAngle(nz-1,height2,&theta2);
	
   printf("theta: %f\n",theta);
	printf("theta: %f\n",theta*(180/pi));
	//printf("%f\n",(pi - theta));
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


	for (i=5;i<=nx-5;i++) {
		for (j=ypos-ydelta;j<=ypos+ydelta;j++) {
			for (k=xpos-xdelta;k<=xpos+xdelta;k++) {
				//printf("%f\t%d\n",colour[i][j][k],i); 
				if (colour[i][j][k] <= 0. && colour[i+1][j][k] > 0.) {
					zPos1 += i + (colour[i][j][k] ) / (colour[i+1][j][k] - colour[i][j][k]);
					valCount1++;
				//	printf("bong\t%f\t%d\n",zPos1,valCount1);
				}
				if (colour[i][j][k] >= 0. && colour[i+1][j][k] < 0.) {
					zPos2 += (i + 1) + (colour[i+1][j][k] ) / ( colour[i][j][k] - colour[i+1][j][k]);
					valCount2++;
				//	printf("bing\t%f\t%d\n",zPos2,valCount2);
				}
			}
		}	
	}
   if (valCount1 != 0)  *z1=zPos1/valCount1;
	if (valCount2 != 0)	*z2=zPos2/valCount2;

	return(0);
}

/// Determine x position of sign change in colour field for given x and y position
/// z1 gives value for negative -> positive, z2 for positive -> negative transition
/// if xdelta and/or ydelta are not 0 mean values for an area +\- delta will be returned

int getZeroColourXstart(FLOATNUM ***colour,
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


	for (i=1;i<=1;i++) {
		for (j=ypos-ydelta;j<=ypos+ydelta;j++) {
			for (k=xpos+xdelta;k>=xpos-xdelta;k--) {
				//printf("%f\t%d\n",colour[i][j][k],i); 
				if (colour[i][j][k] <= 0. && colour[i][j][k-1] > 0.) {
					zPos1 += k - (colour[i][j][k] ) / (colour[i][j][k-1] - colour[i][j][k]);
					valCount1++;
					//printf("bong\t%f\t%d\n",zPos1,valCount1);
				}
				if (colour[i][j][k] >= 0. && colour[i][j][k-1] < 0.) {
					zPos2 += (k-1) - (colour[i][j][k-1] ) / ( colour[i][j][k] - colour[i][j][k-1]);
					valCount2++;
					//printf("bing\t%f\t%d\n",zPos2,valCount2);
				}
			}
		}	
	}
      
        if (valCount1 != 0)
	     {*z1=zPos1/valCount1;}
        if (valCount2 != 0)
	     {*z2=zPos2/valCount2;}

	return(0);
}
/// Determine x position of sign change in colour field for given x and y position
/// z1 gives value for negative -> positive, z2 for positive -> negative transition
/// if xdelta and/or ydelta are not 0 mean values for an area +\- delta will be returned

int getZeroColourXend(FLOATNUM ***colour,
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


	for (i=1;i<=1;i++) {
		for (j=ypos-ydelta;j<=ypos+ydelta;j++) {
			for (k=xpos-xdelta;k<=xpos+xdelta;k++) {
				//printf("%f\t%d\n",colour[i][j][k],i); 
				if (colour[i][j][k] <= 0. && colour[i][j][k+1] > 0.) {
					zPos1 += k + (colour[i][j][k] ) / (colour[i][j][k+1] - colour[i][j][k]);
					valCount1++;
					//printf("bong\t%f\t%d\n",zPos1,valCount1);
				}
				if (colour[i][j][k] >= 0. && colour[i][j][k+1] < 0.) {
					zPos2 += (k + 1) + (colour[i][j][k+1] ) / ( colour[i][j][k] - colour[i][j][k+1]);
					valCount2++;
					//printf("bing\t%f\t%d\n",zPos2,valCount2);
				}
			}
		}	
	}
      
        if (valCount1 != 0)
	     {*z1=zPos1/valCount1;}
        if (valCount2 != 0)
	     {*z2=zPos2/valCount2;}

	return(0);
}
int calcGeometricContactAngle(FLOATNUM base,
			      FLOATNUM height,
			      FLOATNUM *theta)
{
	FLOATNUM	radius = 0.;

	radius = ((4*(pow(height,2)))+(pow(base,2)))/(8 * height);
			printf(" radius: %f\n",radius);

	*theta = atan2((base / 2), (radius - height));// * (180 / pi);

	return(0);
}
