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
	}
	if (*y != 'y' && isdigit(*y)) {
		ly = atoi(y);
		uy = atoi(y);
		sprintf(outfilename,"%s.Y%d",outfilename, ly);
	}
	else {
		ly = 0;
		uy = ny-1;
	}
	if (*z != 'z' && isdigit(*z)) {
		lz = atoi(z);
		uz = atoi(z);
		sprintf(outfilename,"%s.Z%d",outfilename, lz);
	}
	else {
		lz = 0;
		uz = nz-1;
	}
 		
	sprintf(outfilename,"%s.gle.dat",outfilename);

	printf("Writing to: %s\n",outfilename);

	oFile = fopen(outfilename, "w");

	if (lx == ux && ly == uy && lz == uz) {
		// point
		fprintf(oFile,"%f\n", data[lz][ly][lx]);
	}

	if (lx != ux && ly == uy && lz == uz) {
		// x-line
		for (i=lx;i<=ux;i++) {
			fprintf(oFile,"%d %f\n", i, data[lz][ly][i]);
		}
	}

	if (lx == ux && ly != uy && lz == uz) {
		// y-line
		for (j=ly;j<=uy;j++) {
			fprintf(oFile,"%d %f\n", j, data[lz][j][lx]);
		}
	}

	if (lx == ux && ly == uy && lz != uz) {
		// z-line
		for (k=lz;k<=uz;k++) {
			fprintf(oFile,"%d %f\n", k,data[k][ly][lx]);
		}
	}

	if (lx != ux && ly != uy && lz == uz) {
		// xy-plane
		fprintf(oFile,"! NX %d NY %d\n", ux, uy);
		for (i=lx;i<=ux;i++) {
			for (j=ly;j<=uy;j++) {
				fprintf(oFile,"%f ", data[lz][j][i]);
			}
			fprintf(oFile,"\n");
		}
	}

	if (lx != ux && ly == uy && lz != uz) {
		// xz-plane
		//fprintf(oFile,"! NX %d NY %d\n", ux/2-1, 120);
		for (k=241;k<=261;k++) {	
			for (i=1;i<=ux/2;i++) {
				
				fprintf(oFile,"%f ", -data[k][ly][i]);
			}
			fprintf(oFile,"\n%d ",i);
		}
	}

	if (lx == ux && ly != uy && lz != uz) {
		// yz-plane
		fprintf(oFile,"! NX %d NY %d\n", uy, uz);
		for (j=ly;j<=uy;j++) {
			for (k=lz;k<=uz;k++) {
				fprintf(oFile,"%f ", data[k][j][lx]);
			}
			fprintf(oFile,"\n");
		}
	}

	if (lx != ux && ly != uy && lz != uz) {
		// xyz-space
		for (i=lx;i<=ux;i++) {
			for (j=ly;j<=uy;j++) {
				for (k=lz;k<=uz;k++) {
					fprintf(oFile,"%d %d %d %f\n", k, j, i, data[k][j][i]);
				}
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
