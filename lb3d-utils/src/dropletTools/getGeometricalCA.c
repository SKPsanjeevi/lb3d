// getContactAngle.c
//
// Filename:	getContactAngle.c
// Authors:	Sebastian Schmieschek, A. Sarkar
// Version:	0.2
// Last change:	20.03.2008
// 
// Purpose:	Geometrically determine the contact angle of a droplet 
//		which is above a surface at z = 2 by evaluating the colour field.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "hdfread.h"
#include "getContactAngle.h"


int main(int argc,
	 char *argv[])
{
	int	       retVal = 1; // 0 on success

	char	*filename;

	struct geometry_t CAresult;
	FLOATNUM	***colour;
	int		measureHeight;
	FLOATNUM	threshold;

	
	filename =		argv[1];
	measureHeight =		argv[2] ? atoi(argv[2]) : 5;
	threshold =		argv[3] ? atof(argv[3]) : 0.2;

	if (filename) {
		getDims(filename,&nx,&ny,&nz);
		hdfread(filename, "colour", &colour);
		retVal = 0;
	}

	calcGeometricalContactAngle(colour, measureHeight, threshold, &CAresult);

	if ( DEBUG == 1 ) {
		printf("AboveGround:\t%f\nheight:\t\t%f\t%f\t+%f\nbase:\t\t%f\t%f\t+%f\nradius:\t\t%f\t%f\t+%f\ntheta:\t\t%f\t%f\t+%f\n",
		       CAresult.aboveGround, CAresult.height, CAresult.minHeight, CAresult.maxHeight, CAresult.base, CAresult.maxBase, CAresult.minBase, CAresult.radius, CAresult.minRadius, CAresult.maxRadius, CAresult.theta, CAresult.minTheta, CAresult.maxTheta);
	}
	else {
		printf("%f\t%f\t%f\n",CAresult.theta, CAresult.minTheta, CAresult.maxTheta);
	}
}

int calcGeometricalContactAngle(FLOATNUM ***colour,
				int tmpHeightOffset,
				FLOATNUM colourThreshold,
				struct geometry_t *result )
{

	int		i,j,k,l;

	FLOATNUM	virtPos = 0.;
	FLOATNUM	lowZ = 0.;
	FLOATNUM	highZ = 0.;
	FLOATNUM	lowY = 0.;
	FLOATNUM	highY = 0.;
	FLOATNUM	aboveGround = 0.;

	FLOATNUM	tmpHeight = 0.;
	FLOATNUM	tmpBase = 0.;
	FLOATNUM	height = 0.;
	FLOATNUM	base = 0.;
	FLOATNUM	diam = 0.;
	FLOATNUM	radius = 0.;
	FLOATNUM	theta = 0.;
	FLOATNUM	deltaColour = 0.;
    

	int		minLowZ = 0, minLowY = 0;
	int		maxLowZ = 0, maxLowY = 0;
	int		minHighZ = 0, minHighY = 0;
	int		maxHighZ = 0, maxHighY = 0;

	FLOATNUM	minBase = 0.0, maxBase = 0.0;
	int		minHeight = 0, maxHeight = 0;
	int		tmpMinHeight = 0, tmpMaxHeight = 0;
	int		tmpMinBase = 0.0, tmpMaxBase = 0.0;
	FLOATNUM	minRadius = 0.0, maxRadius = 0.0;
	FLOATNUM	minTheta = 0.0, maxTheta = 0.0;

	j = ny / 2;
	k = nx / 2;
	for (i = 0; i < nz-3; i++) {
		//if (i > 1) {
		deltaColour = colour[i+1][j][k] - colour[i][j][k];
		//}
		//printf("%f\t%f\n",colour[i][j][k],deltaColour);
		if ( ((colour[i][j][k] <= 0.) && (colour[i+1][j][k] > 0.))) { 
			l = 1;
			while (fabs(colour[i-l][j][k]) < colourThreshold && (i-l) != 0) {
				l++;
			}
			
			minLowZ = i-l;
			l = 1;
			while (fabs(colour[i+l][j][k]) < colourThreshold) {
				l++;
			}
			maxLowZ = i+l;
			virtPos = i - (colour[i][j][k] ) / (colour[i+1][j][k] - colour[i][j][k]);
			lowZ = virtPos;
			
		}
		if ( ((colour[i][j][k] >= 0.) && (colour[i+1][j][k] < 0.))) { // && (colourThreshold == 0.))
			//   || ( deltaColour < -colourThreshold ) ) {
			l = 1;
			while (fabs(colour[i-l][j][k]) < colourThreshold && (i-l) != 0) {
				l++;
			}
			
			minHighZ = i-l;
			l = 1;
			while (fabs(colour[i+l][j][k]) < colourThreshold) {
				l++;
			}
			maxHighZ = i+l;
			virtPos = (i + 1) + (colour[i+1][j][k] ) / ( colour[i][j][k] - colour[i+1][j][k]);
			highZ = virtPos;
			//highZ = i;
		}
		//printf("lowZhighZ:%f\t%f\n",lowZ,highZ);
	}
	

	if (lowZ != 1.) {
		aboveGround = lowZ;
	}
	else {
		// Colour change from zero is found in rock-area, assuming the drop has contact
		// and set its base to the first non-rock z value.
		lowZ = 2.;
	}
	tmpHeight = highZ - lowZ - tmpHeightOffset; 
	
	i = 2 + tmpHeightOffset;
	for (j = 1; j < ny-1; j++) {
		//if (j > 1) {
		deltaColour = colour[i][j][k] - colour[i][j-1][k];
		//}
		//printf("%f\t%f\n",colour[i][j][k],deltaColour);
		if ( ((colour[i][j][k] <= 0.) && (colour[i][j+1][k] > 0.) )) {//&& (colourThreshold == 0.)) 
			//     || ( deltaColour > colourThreshold ) ) {
			l = 1;
			while (fabs(colour[i][j-l][k]) < colourThreshold && (j-l) != 0) {
				l++;
			}
			
			minLowY = j-l;
			l = 1;
			while (fabs(colour[i][j+l][k]) < colourThreshold) {
				l++;
			}
			maxLowY = j+l;
			virtPos = j - (colour[i][j][k] ) / ( colour[i][j+1][k] - colour[i][j][k] );
			lowY = virtPos;
		}
		if ( ((colour[i][j][k] >= 0.) && (colour[i][j+1][k] < 0.) )) { // && (colourThreshold == 0.)) 
			//    || ( deltaColour < -colourThreshold ) ) {
			l = 1;
			while (fabs(colour[i][j-l][k]) < colourThreshold && (j-l) != 0) {
				l++;
			}
			
			minHighY = j-l;
			l = 1;
			while (fabs(colour[i][j+l][k]) < colourThreshold) {
				l++;
			}
			maxHighY = j+l;
			virtPos = (j + 1) + ( colour[i][j+1][k] ) / ( colour[i][j][k] - colour[i][j+1][k]);
			highY = virtPos;
		}
		//printf("%f\t%f\n",lowY,highY);
	}

	tmpBase = highY - lowY;
	
	tmpMinHeight = minHighZ - maxLowZ - tmpHeightOffset;
	tmpMaxHeight = maxHighZ - minLowZ - tmpHeightOffset;
	tmpMinBase = minHighY - maxLowY;
	tmpMaxBase = maxHighY - minLowY;


	minHeight = tmpMinHeight + tmpHeightOffset;
	height = tmpHeight + tmpHeightOffset;
	maxHeight = tmpMaxHeight + tmpHeightOffset;

	if (aboveGround != 0.0 ) {
		minRadius = minHeight;
		radius = height;
		maxRadius = maxHeight;
		
		minBase = 0.0;
		base = 0.0;
		maxBase = 0.0;
	}
	else {
		minRadius = ((4*(pow(tmpMinHeight,2)))+(pow(tmpMaxBase,2)))/(8 * tmpMinHeight);
		radius = ((4*(pow(tmpHeight,2)))+(pow(tmpBase,2)))/(8 * tmpHeight);
		maxRadius = ((4*(pow(tmpMaxHeight,2)))+(pow(tmpMinBase,2)))/(8 * tmpMaxHeight);

		minBase = 2* sqrt(pow(minRadius,2)-pow((minHeight-minRadius),2));
		base = 2* sqrt(pow(radius,2)-pow((height-radius),2));
		maxBase = 2* sqrt(pow(maxRadius,2)-pow((maxHeight-maxRadius),2));
	}

	minTheta = atan2((minBase / 2), (minRadius - minHeight)) * (180 / pi);
	theta = atan2((base / 2), (radius - height)) * (180 / pi);
	maxTheta = atan2((maxBase / 2), (maxRadius - maxHeight)) * (180 / pi);

	if (DEBUG == 1) {
    
	printf("%d %f %d\n",minLowZ,lowZ,maxLowZ);
	printf("%d %f %d\n",minHighZ,highZ,maxHighZ);

	printf("%d %f %d\n",minLowY,lowY,maxLowY);
	printf("%d %f %d\n",minHighY,highY,maxHighY);

	printf("%d %d\n", minHeight, maxHeight);
	printf("%f %f\n", minBase, maxBase);
	printf("%f %f\n",minRadius,maxRadius);

	printf("%f %f\n",minTheta,maxTheta);

	}

	result->aboveGround = aboveGround;
	result->minHeight = minHeight - height;
	result->height = height;
	result->maxHeight = maxHeight - height;
	result->minBase = minBase - base;
	result->base = base;
	result->maxBase = maxBase - base;
	result->minRadius = minRadius - radius;
	result->radius = radius;
	result->maxRadius = maxRadius - radius;
	result->minTheta = minTheta - theta;
	result->theta = theta;
	result->maxTheta = maxTheta - theta;

	return(0);
}
