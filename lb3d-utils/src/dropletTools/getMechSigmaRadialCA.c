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

	struct sigma_t mechRadialSigmaResult;

	FLOATNUM	offset;
	FLOATNUM	radius;
	FLOATNUM	height;

	FLOATNUM	***redDens;
	FLOATNUM	***blueDens;

	FLOATNUM	***colour;

	FLOATNUM	***pxx;
	FLOATNUM	***pyy;
	FLOATNUM	***pzz;

	FLOATNUM	***pxy;
	FLOATNUM	***pxz;
	FLOATNUM	***pyz;

	filename = argv[1];


	if (filename) {
		getDims(filename,&nx,&ny,&nz);
		hdfread(filename, "colour", &colour);
		hdfread(filename, "pxx", &pxx);
		hdfread(filename, "pxy", &pxy);
		hdfread(filename, "pyy", &pyy);
		hdfread(filename, "pyz", &pyz);
		hdfread(filename, "pzz", &pzz);
		hdfread(filename, "pxz", &pxz);
		hdfread(filename, "wd", &blueDens);
		hdfread(filename, "od", &redDens);
		retVal = 0;
	}

	getGeometry(colour, 5, &offset, &radius, &height);

	calcMechSigmaRadial(pxx, pxy, pyz, pxz, pzz, offset, radius, height, &mechRadialSigmaResult);
    
	printf("%f\n",mechRadialSigmaResult.theta);
}




int calcMechSigmaRadial(FLOATNUM ***pxx,
			FLOATNUM ***pxy,
			FLOATNUM ***pyz,
			FLOATNUM ***pxz,
			FLOATNUM ***pzz,
			FLOATNUM offset,
			FLOATNUM radius,
			FLOATNUM height,
			struct sigma_t *result)
{
    
	int         i, j1, k1, j2, k2, j1center, k1center, j2center, k2center, j1delta, k1delta, j2delta, k2delta, dropCenterZ;
	int		SGvalCount=0;
	int		LGvalCount=0;
	int		SLvalCount=0;
    
	FLOATNUM	deltaPlocal1[nz];
	FLOATNUM	deltaPlocal2[nz];
    
	FLOATNUM	sigmaSL = 0.;
	FLOATNUM	sigmaSG = 0.;
	FLOATNUM	sigmaLG = 0.;
	FLOATNUM	sigmaLG2 = 0.;

	FLOATNUM	theta;

	dropCenterZ = (int) (offset + height - radius);
	//printf("centerZ:\t%u\n",dropCenterZ);
	j1center = ny / 2;
	k1center = nx / 2;
	j2center = 5;
	k2center = 5;

	j1delta = 0;
	k1delta = 0;
	j2delta = 0;
	k2delta = 0;

	for (j1 = j1center - j1delta; j1 <= j1center + j1delta; j1++) {
	for (k1 = k1center - k1delta; k1 <= k1center + k1delta; k1++) {
	for (j2 = j2center - j2delta; j2 <= j2center + j2delta; j2++) {
	for (k2 = k2center - k2delta; k2 <= k2center + k2delta; k2++) {

		deltaPlocal1[0] = 0.;
		deltaPlocal1[0] = ( pzz[0][j1][k1]  -  pxx[0][j1][k1] );
		deltaPlocal2[0] = 0.;
		deltaPlocal2[0] = ( pzz[0][j2][k2]  -  pxx[0][j2][k2] );
		for (i = 0; i < nz - 1; i++) {
			deltaPlocal1[i+1] = 0.;
			deltaPlocal1[i+1] = ( pzz[i+1][j1][k1]  - pxx[i+1][j1][k1] );
			deltaPlocal2[i+1] = 0.;
			deltaPlocal2[i+1] = ( pzz[i+1][j2][k2]  -  pxx[i+1][j2][k2] );

			if (i < 15) {

 				sigmaSL = sigmaSL + ( deltaPlocal1[i] ); 
 				sigmaSG = sigmaSG + ( deltaPlocal2[i] ); 

			}
			else if (i > (int) (dropCenterZ && i < (int) nz - 15)) {

				sigmaLG = sigmaLG + ( deltaPlocal1[i]  * (i/(int) (offset+height)) );

			}
			
		}
	}
	}
	}
	}       

	theta = acos( (sigmaSG -sigmaSL ) / sigmaLG ) * ( 180 / pi);
	result->sigmaSL = sigmaSL;
	result->sigmaSG = sigmaSG;
	result->sigmaLG = sigmaLG;
	result->theta = theta;
	if ( DEBUG ==1 ) {
		printf("sigmaSL:\t%f\nsigmaSG:\t%f\nsigmaLG:\t%f\ntheta:\t\t\t\t%f\n",sigmaSL,sigmaSG,sigmaLG, theta);
	}
	
		
		
	return(0);

}

int getGeometry(FLOATNUM ***colour,
		int tmpHeightOffset,
		FLOATNUM *offset,
		FLOATNUM *radius,
		FLOATNUM *height)
{

	int		i,j,k;
	
	FLOATNUM	virtPos = 0.;
	FLOATNUM	lowZ = 0.;
	FLOATNUM	highZ = 0.;
	FLOATNUM	lowY = 0.;
	FLOATNUM	highY = 0.;
	FLOATNUM	aboveGround = 0.;
	
	FLOATNUM	tmpHeight = 0.;
	FLOATNUM	tmpBase = 0.;
	FLOATNUM	base = 0.;
	FLOATNUM	deltaColour = 0.;

	j = ny / 2;
	k = nx / 2;
	for (i = 0; i < nz-3; i++) {
		deltaColour = colour[i+1][j][k] - colour[i][j][k];
		if ( ((colour[i][j][k] <= 0.) && (colour[i+1][j][k] > 0.))) { 	
			virtPos = i - (colour[i][j][k] ) / (colour[i+1][j][k] - colour[i][j][k]);
			lowZ = virtPos;
		}
		if ( ((colour[i][j][k] >= 0.) && (colour[i+1][j][k] < 0.))) { // && (colourThreshold == 0.))
			virtPos = (i + 1) + (colour[i+1][j][k] ) / ( colour[i][j][k] - colour[i+1][j][k]);
			highZ = virtPos;
		}
	}
	

	if (lowZ != 1.) {
		aboveGround = lowZ;
	}
	else {
		lowZ = 2.;
	}
	tmpHeight = highZ - lowZ - tmpHeightOffset; 
	
	i = 2 + tmpHeightOffset;
	for (j = 1; j < ny-1; j++) {
		deltaColour = colour[i][j][k] - colour[i][j-1][k];
		if ( ((colour[i][j][k] <= 0.) && (colour[i][j+1][k] > 0.) )) {
			virtPos = j - (colour[i][j][k] ) / ( colour[i][j+1][k] - colour[i][j][k] );
			lowY = virtPos;
		}
		if ( ((colour[i][j][k] >= 0.) && (colour[i][j+1][k] < 0.) )) {
			virtPos = (j + 1) + ( colour[i][j+1][k] ) / ( colour[i][j][k] - colour[i][j+1][k]);
			highY = virtPos;
		}
	}

	tmpBase = highY - lowY;
	
	*height = tmpHeight + tmpHeightOffset;

	if (aboveGround != 0.0 ) {
		*radius = *height;
		base = 0.0;
	}
	else {
		*radius = ((4*(pow(tmpHeight,2)))+(pow(tmpBase,2)))/(8 * tmpHeight);
		base = 2* sqrt(pow(*radius,2)-pow((height-radius),2));
	}

	*offset = aboveGround;

	return(0);
}
