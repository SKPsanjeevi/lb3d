/*includes*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "hdfread.h"

/* declarations */
int getIsoColourZ(FLOATNUM ***colour,
		  int xpos,
		  int ypos,
		  FLOATNUM *z1,
		  FLOATNUM *z2,
		  int xdelta,
		  int ydelta,
		  FLOATNUM threshold);
/* globals */
int		nx, ny, nz;

/* definitions */

int main(int argc,
	 char *argv[])
{
	int		retVal = 1; // 0 on success

	int		xPos,yPos;
	int		xDelta, yDelta;

	FLOATNUM	zPos1, zPos2;
	FLOATNUM	***colour;
	char		*filename;
	int		direction;
	FLOATNUM	threshold;

	filename = argv[1];
	xPos = atoi(argv[2]);
	yPos = atoi(argv[3]);

	xDelta = atoi(argv[4]);
	yDelta = atoi(argv[5]);

	direction = atoi(argv[6]);

	threshold = atof(argv[7]);

	if (filename) {
		getDims(filename,&nx,&ny,&nz);
		//hdfread(filename, "colour",&colour);
		hdfread2(filename, &colour);
		retVal = 0;
	}
	
	getIsoColourZ(colour,xPos,yPos,&zPos1,&zPos2,xDelta,yDelta,threshold);
	zPos2 = zPos2; 
	if (direction == 0) {
		printf("%f\n",zPos2);
	}
	else if (direction == 1) {
		printf("%f %f\n",threshold,zPos1);
	}
	
	return(retVal);
}




int getIsoColourZ(FLOATNUM ***colour,
		   int xpos,
		   int ypos,
		   FLOATNUM *z1, 
		   FLOATNUM *z2,
		   int xdelta,
		  int ydelta,
		  FLOATNUM threshold)

{

	int		i,j,k;
	int		valCount1, valCount2;
	FLOATNUM	zPos1, zPos2;

	valCount1 = 0;
	valCount2 = 0;


	for (i=10;i<=nx-12;i++) {
		for (j=ypos-ydelta;j<=ypos+ydelta;j++) {
			for (k=xpos-xdelta;k<=xpos+xdelta;k++) {
				//printf("%f\t%d\n",colour[i][j][k],i); 
				if (colour[i][j][k] <= threshold && colour[i+1][j][k] > threshold) {
					zPos1 += (i+1) + (colour[i][j][k] ) / (colour[i+1][j][k] - colour[i][j][k]);
					valCount1++;
					//printf("bong\t%f\t%d\n",zPos1,valCount1);
				}
				if (colour[i][j][k] >= threshold && colour[i+1][j][k] < threshold) {
					zPos2 += (i) - (colour[i+1][j][k] ) / ( colour[i][j][k] - colour[i+1][j][k]);
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
