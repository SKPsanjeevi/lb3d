// getMinMax.c
//
// Filename:	getMinMax.c
// Authors:	Sebastian Schmieschek
// Version:	0.1
// Last change:	21.08.2008
// 
// Purpose:	Read an LB3D hdf output file and determine minimum and maximum value

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "hdfread.h"
#include "getContactAngle.h"


int main(int argc,
	 char *argv[])
{
	int		retVal = 1; // 0 on success

	char		*filename;
	FLOATNUM	***data;

	FLOATNUM	threshold;
	FLOATNUM	sumDens;
	FLOATNUM	avgDens;

	int		latticeSites = 0.0;
	int		i=0,j=0,k=0;
	int		nx=0,ny=0,nz=0;

	filename = argv[1];
        threshold = atof(argv[2]);

	
	if (filename) {
		getDims(filename,&nx,&ny,&nz);
		hdfread2(filename, &data);
		retVal = 0;
	}

	for (i=0;i<nx;i++) {
		for (j=0;j<ny;j++) {
			for (k=0;k<nz;k++) {
			  if (data[i][j][k] > threshold) {
			  sumDens += data[i][j][k];
			  latticeSites++;
			  }
			}
		}
	}
	avgDens = sumDens / latticeSites;

	printf("dens: %f on %d sites (%f average)\n",sumDens,latticeSites,avgDens);

	return retVal;
   
}
