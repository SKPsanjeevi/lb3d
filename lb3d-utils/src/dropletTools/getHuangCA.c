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

	struct sigma_t huangResult;

	FLOATNUM	gBR;
	FLOATNUM	rockC;

	FLOATNUM	***redDens;
	FLOATNUM	***blueDens;

	FLOATNUM	***colour;

	FLOATNUM	***pxx;
	FLOATNUM	***pyy;
	FLOATNUM	***pzz;

	FLOATNUM	***pxy;
	FLOATNUM	***pxz;
	FLOATNUM	***pyz;

	FLOATNUM	invRadius;

	filename = argv[1];
	gBR = argv[2] ? atof(argv[2]) : usage();
	rockC = argv[3] ? atof(argv[3]) : usage();

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

	calcSigmaHuang2(redDens, blueDens, gBR, rockC, &huangResult );
	printf("%f\n",huangResult.theta);
}

usage()
{
	printf("so nu ma nech mien jong.\n");
	exit(0);
}

calcSigmaHuang2(FLOATNUM ***redDens,
	       FLOATNUM ***blueDens,
	       FLOATNUM gBR,
	       FLOATNUM rockC,
	       struct sigma_t *result )
{
	int         i, j, k;

	FLOATNUM	sigmaSL = 0.;
	FLOATNUM	sigmaSG = 0.;
	FLOATNUM	sigmaLG = 0.;

	FLOATNUM	maxRedDens = 0.;
	FLOATNUM	minRedDens = 1.;

	FLOATNUM	maxBlueDens = 0.;
	FLOATNUM	minBlueDens = 1.;

	FLOATNUM	gROCK = gBR * (6./19. * (1-exp(-rockC)));

	FLOATNUM	theta;

	j = ny / 2;
	k = nx / 2;

	for (i = 3; i < nz; i++) {
		if (redDens[i][j][k] > maxRedDens) {
			maxRedDens = redDens[i][j][k];
		}
		if (redDens[i][j][k] < minRedDens) {
			minRedDens = redDens[i][j][k];
		}

		if (blueDens[i][j][k] > maxBlueDens) {
			maxBlueDens = blueDens[i][j][k];
		}
		if (blueDens[i][j][k] < minBlueDens) {
			minBlueDens = blueDens[i][j][k];
		}
	}
 
	sigmaSL = gROCK;

	sigmaSG = -gROCK;

 	sigmaLG = gBR * ( ( maxRedDens - minRedDens ) / 2 ); 

	theta = acos( ( sigmaSL - sigmaSG ) / sigmaLG ) * ( 180 / pi);
	result->sigmaSL = sigmaSL;
	result->sigmaSG = sigmaSG;
	result->sigmaLG = sigmaLG;
	result->theta = theta;
	
	if (DEBUG == 1) {
		printf("maxRedDens:\t%f\nminRedDens:\t%f\n",maxRedDens,minRedDens);
		printf("sigmaSL:\t%f\nsigmaSG:\t%f\nsigmaLG:\t%f\ntheta:\t\t\t\t%f\n", sigmaSL, sigmaSG, sigmaLG,theta);
	}
}
