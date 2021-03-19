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

	FLOATNUM	maxVal=-5000.0;
	FLOATNUM	minVal=5000.0;

	int		i=0,j=0,k=0;
	int		nx=0,ny=0,nz=0;
	int		outtype;

	filename = argv[1];
	outtype = argv[2] ? atoi(argv[2]) : 2;
	
	if (filename) {
		getDims(filename,&nx,&ny,&nz);
		hdfread2(filename, &data);
		retVal = 0;
	}

	for (i=0;i<nx;i++) {
		for (j=0;j<ny;j++) {
			for (k=0;k<nz;k++) {
				if (abs(data[i][j][k]) > 5000.0) {
					printf("values out of range");
					return(-1);
				}
				else {
					if (data[i][j][k] > maxVal) {
						maxVal = data[i][j][k];
					}
					if (data[i][j][k] < minVal) {
						minVal = data[i][j][k];
					}
					if (data[i][j][k] < 0.0) {
					  //printf("%d,%d,%d: %f\n",i,j,k,data[i][j][k]);
					} 
				}

			}
		}
	}

	switch (outtype) {
	case 0:
		printf("%f\n",minVal);
		break;
	case 1:
		printf("%f\n",maxVal);
		break;
	default:
		printf("%f\t%f\n",minVal,maxVal);
	}
	return retVal;
   
}
