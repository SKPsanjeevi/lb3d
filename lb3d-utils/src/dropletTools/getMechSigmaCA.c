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

	char		*filename;

	struct sigma_t	mechSigmaResult;

	FLOATNUM	***colour;
	FLOATNUM	***redDens;

	FLOATNUM	***pxx;
	FLOATNUM	***pyy;
	FLOATNUM	***pzz;

	FLOATNUM	offset = 0.0;
	FLOATNUM	height = 0.0;

	filename = argv[1];

	if (filename) {
		getDims(filename,&nx,&ny,&nz);
		hdfread(filename, "colour", &colour);
		hdfread(filename, "pxx", &pxx);
		hdfread(filename, "pyy", &pyy);
		hdfread(filename, "pzz", &pzz);
		hdfread(filename, "od", &redDens);
		retVal = 0;
	}
	
	getHeight(&offset, &height, colour);

	calcSigmaMech(pxx, pzz, 0, height, 0, &mechSigmaResult);

	printf("%f\n",mechSigmaResult.theta);

}
int calcSigmaMech(FLOATNUM ***pxx,
		  FLOATNUM ***pzz,
		  FLOATNUM offset,
		  FLOATNUM height,
		  FLOATNUM radius,
		  struct sigma_t *result)
{
  
	int         i, j1, k1, j2, k2, dropCenterZ;

	FLOATNUM	deltaPlocal1[nz];
	FLOATNUM	deltaPlocal2[nz];
    
	FLOATNUM	sigmaSL = 0.;
	FLOATNUM	sigmaSG = 0.;
	FLOATNUM	sigmaLG = 0.;

	FLOATNUM	theta;

	FLOATNUM	factor = 1.;
    

	j1 = ny / 2;
	k1 = nx / 2;
	j2 = 5;
	k2 = 5;
	deltaPlocal1[0] = 0.;
	deltaPlocal1[0] = ( pzz[0][j1][k1] ) - ( pxx[0][j1][k1] ); 
	deltaPlocal2[0] = 0.;
	deltaPlocal2[0] = ( pzz[0][j2][k2] ) - ( pxx[0][j2][k2] );    
	for (i = 0; i < nz - 1; i++) {
		deltaPlocal1[i+1] =0.;
		deltaPlocal1[i+1] = ( pzz[i+1][j1][k1] ) - ( pxx[i+1][j1][k1] );
		deltaPlocal2[i+1] =0.;
		deltaPlocal2[i+1] = ( pzz[i+1][j2][k2] ) - ( pxx[i+1][j2][k2] );

		if (i < 15 && i > 1) {
			
			sigmaSL = sigmaSL + (factor * deltaPlocal1[i]);
			sigmaSG = sigmaSG + (factor * deltaPlocal2[i]);
			
		}
		else if (i > (int) (height - 15) && i < (int) (offset + height + 15)) {
	    
			sigmaLG = sigmaLG + (factor * deltaPlocal1[i]); 

		}	
	}
	theta = acos( ( sigmaSG - sigmaSL ) / sigmaLG ) * ( 180 / pi);
	result->sigmaSL = sigmaSL;
	result->sigmaSL = sigmaSG;
	result->sigmaSL = sigmaLG;
	result->theta = theta;

	if ( DEBUG == 1 ) {
		printf("sigmaSL:\t%f\nsigmaSG:\t%f\nsigmaLG:\t%f\ntheta:\t\t\t\t%f\n",sigmaSL,sigmaSG,sigmaLG, theta);
	}
    
    

	return(0);

}

int getHeight(FLOATNUM *offset,
	      FLOATNUM *height, 
	      FLOATNUM ***colour)
{

	int		i,j,k;
	FLOATNUM	deltaColour;
	int		minLowZ = 0, lowZ = 0, maxLowZ = 0;
	int		minHighZ = 0, highZ = 0, maxHighZ = 0;

	FLOATNUM	virtPos;

	j = ny / 2;
	k = nx / 2;
	for (i = 0; i < nz-3; i++) {
		deltaColour = colour[i+1][j][k] - colour[i][j][k];
		if ( ((colour[i][j][k] <= 0.) && (colour[i+1][j][k] > 0.))) { 
			virtPos = i - (colour[i][j][k] ) / (colour[i+1][j][k] - colour[i][j][k]);
			lowZ = virtPos;
		}
		if ( ((colour[i][j][k] >= 0.) && (colour[i+1][j][k] < 0.))) { 
			virtPos = (i + 1) + (colour[i+1][j][k] ) / ( colour[i][j][k] - colour[i+1][j][k]);
			highZ = virtPos;
		}
	}
      	if (lowZ != 1.) {
		*offset = lowZ;
	}
	else {
		lowZ = 2.;
	}

	*height = highZ;

	return 0;
}
