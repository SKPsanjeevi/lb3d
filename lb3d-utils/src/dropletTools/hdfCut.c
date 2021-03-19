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
	char		outfilename[70];

	FILE		*oFile;

	FLOATNUM	***data;

	FLOATNUM	maxVal=-5000.0;
	FLOATNUM	minVal=5000.0;

	int		i=0,j=0,k=0;
	int		nx=0,ny=0,nz=0;
	int		lx=-1,ly=-1,lz=-1;
	int		ux=-1,uy=-1,uz=-1;
	char		*x,*y,*z;

	if (argc != 5) usage();


	filename	= argv[1];// ? argv[1] : usage();
	x		= argv[2];// ? argv[2] : usage();
	y		= argv[3];// ? argv[3] : usage();
	z		= argv[4];// ? argv[4] : usage();
	

	//printf("%s %s %s %s\n", filename,  x, y, z);
	if (filename) {
		getDims(filename,&nz,&ny,&nx);
		hdfread2(filename, &data);
		retVal = 0;
		sprintf(outfilename,"%s",filename);
	}

	if (*x != 'x'&& isdigit(*x)) {
		lx = atoi(x);
		ux = atoi(x);
		sprintf(outfilename,"%s.X%d",outfilename, lx);
	}
	else {
		lx = 0;
		ux = nx-1;
		//lx = 1;
		//ux = nx/2;//nx-2;
	}
	if (*y != 'y' && isdigit(*y)) {
		ly = atoi(y);
		uy = atoi(y);
		sprintf(outfilename,"%s.Y%d",outfilename, ly);
	}
	else {
		ly = 0;
		uy = ny-1;
		// surface luwa256
		//ly = 1;
		//uy = ny-2;
	}
	if (*z != 'z' && isdigit(*z)) {
		lz = atoi(z);
		uz = atoi(z);
		sprintf(outfilename,"%s.Z%d",outfilename, lz);
	}
	else {
		lz = 0;
		uz = nz-1;
		// surface luwa256
		//lz = 20;
		//uz = 140; //nz-2;
	}

	sprintf(outfilename,"%s.dat",outfilename);

	printf("Writing to: %s\n",outfilename);

	oFile = fopen(outfilename, "w");


	for (i=lx;i<=ux;i++) {
		for (j=ly;j<=uy;j++) {
			for (k=lz;k<=uz;k++) {
				if (lx != ux) {fprintf(oFile,"%d ", i);}
				if (ly != uy) {fprintf(oFile,"%d ", j);}
				if (lz != uz) {fprintf(oFile,"%d ", k);}
				fprintf(oFile,"%f\n",data[k][j][i]);
			}
		}
	}

	fclose(oFile);

	//printf("%d %d %d %d %d %d\n",lx,ly,lz,ux,uy,uz);

	return retVal;
   
}

int usage()
{

	printf("usage: <file> <x> <y> <z>\n\n{x,y,z}:\tx,y,z if variable,\n\t\tinteger if fixed\n");
	exit(0);
}
