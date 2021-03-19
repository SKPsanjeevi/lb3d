/*includes*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "hdfread.h"

//#define DEBUG 0;

/* declarations */
int 
getCapillaryNumber(FLOATNUM mu,
		   FLOATNUM v,
		   FLOATNUM sigma,
		   FLOATNUM *Ca);

int 
getZeroColourZ(FLOATNUM ***colour,
		   int xpos,
		   int ypos,
		   FLOATNUM *z1,
		   FLOATNUM *z2,
		   int xdelta,
		   int ydelta);

int
calcGeometricalRadius(FLOATNUM base,
		      FLOATNUM height,
		      FLOATNUM *radius);


int 
calcGeometricContactAngle(FLOATNUM base,
			  FLOATNUM height,
			  FLOATNUM radius,
			  FLOATNUM *theta);

int 
getMechDefSigma(FLOATNUM ***pxx,
		FLOATNUM ***pzz,
		int xpos,
		int ypos,
		FLOATNUM *sigma,
		int xdelta,
		int ydelta);

void
getCapillaryLaplaceSigma(FLOATNUM ***pxx,
		FLOATNUM ***pyy,
		FLOATNUM ***pzz,
		FLOATNUM *sigma,
		int lowerZ,
		int higherZ,
		int offset,
		FLOATNUM height,
		FLOATNUM theta
		);

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

	FLOATNUM	zPos1r, zPos2r; // right
	FLOATNUM	zPos1c, zPos2c; // center 
	FLOATNUM	zPos1l, zPos2l; // left
	FLOATNUM	***colour;
	char		*filename;
	char		*outfilename;

	int		direction;

	FLOATNUM	height1, height2;
	FLOATNUM	theta1,theta2;

	FLOATNUM	radius, radius1, radius2;	

	FLOATNUM	***pxx, ***pyy, ***pzz;

	FLOATNUM	sigma, sigma2;

	/* Lucas Washburn variables */
	
	FLOATNUM	v_cap = 0.0;
	FLOATNUM	t_d = 0.0;
	FLOATNUM	z1 = 0.0;
	FLOATNUM	z2 = 0.0;
	FLOATNUM	rho = 0.0;
	FLOATNUM	theta = 0.0;
	FLOATNUM	mu = 0.0;
	int		timestep;
	int		maxTime=0;
	int		height=127;
	int		length=448;
	int		z_0=32;
	int DEBUG = 0;


	filename = argv[1];
	
	if (filename) {
		getDims(filename,&nx,&ny,&nz);
		//	hdfread(filename, "colour",&colour);
		hdfread2(filename, &colour);
		hdfread(filename, "pxx", &pxx);
		hdfread(filename, "pyy", &pyy);
		hdfread(filename, "pzz", &pzz);

		
		//	hdfwrite(outfilename, nx, ny, nz, colour);

		retVal = 0;
	}
	


	getZeroColourZ(colour,1,4,&zPos1r,&zPos2r,0,0);
	getZeroColourZ(colour,(height/2),4,&zPos1c,&zPos2c,0,0);
	getZeroColourZ(colour,(height-1),4,&zPos1l,&zPos2l,0,0);

	height1 = zPos1r - zPos1c;
	height2 = zPos2r - zPos2c;
 

 	getCapillaryLaplaceSigma(pxx,pyy,pzz,&sigma2,zPos1c,zPos2c,0,height,2.2);

		printf("sigma: %f\n",sigma2);
		//}
	
	//sigma = 0.07;

	return(retVal);
}

int getCapillaryNumber(FLOATNUM mu,
		       FLOATNUM v,
		       FLOATNUM sigma,
		       FLOATNUM *Ca)
{

	*Ca = ( mu * v ) / sigma;
	return (0);

}


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

	//	printf("blabla\n");

	for (i=10;i<=nx-12;i++) {
		for (j=ypos-ydelta;j<=ypos+ydelta;j++) {
			for (k=xpos-xdelta;k<=xpos+xdelta;k++) {
				//printf("%f\t%d\n",colour[i][j][k],i); 
				if (colour[i][j][k] <= 0. && colour[i+1][j][k] > 0.) {
					zPos1 += i + (colour[i][j][k] ) / (colour[i+1][j][k] - colour[i][j][k]);
					valCount1++;
					//printf("bong\t%f\t%d\n",zPos1,valCount1);
				}
				if (colour[i][j][k] >= 0. && colour[i+1][j][k] < 0.) {
					zPos2 += (i + 1) + (colour[i+1][j][k] ) / ( colour[i][j][k] - colour[i+1][j][k]);
					valCount2++;
					//printf("bing\t%f\t%d\n",zPos2,valCount2);
				}
			}
		}	
	}
     
	printf("%f %f\n",zPos1/valCount1,zPos2/valCount2);

	*z1=zPos1/valCount1;
	*z2=zPos2/valCount2;

	return(0);
}

int calcGeometricalRadius(FLOATNUM base,
			  FLOATNUM height,
			  FLOATNUM *radius)
{
	*radius = ((4*(pow(height,2)))+(pow(base,2)))/(8 * height);
	//	printf("base height radius: %f\n",radius);
	return(0);
}

int calcGeometricContactAngle(FLOATNUM base,
			      FLOATNUM height,
			      FLOATNUM radius,
			      FLOATNUM *theta)
{
	FLOATNUM	bla;
	*theta = atan2((base / 2), (radius - height)) ;//* (180 / pi);
	bla = atan2((base / 2), (radius - height)) * (180 / pi);
	printf("%f\n",bla);
	return(0);
}
int getMechDefSigma(FLOATNUM ***pxx,
	      FLOATNUM ***pzz,
	      int xpos,
	      int ypos,
	      FLOATNUM *sigma, 
	      int xdelta,
	      int ydelta)

{

	int		i,j,k,l;
	FLOATNUM	deltaPlocal[nx];
	FLOATNUM	sigma1 = 0.;

	//	printf("nx: %d, ny: %d, nz: %d\n",nx,ny,nz);
	l=0;

	for (j=ypos-ydelta;j<=ypos+ydelta;j++) {
		for (k=xpos-xdelta;k<=xpos+xdelta;k++) {
			deltaPlocal[0] = 0.;
			deltaPlocal[0] = ( pzz[0][j][k] ) - ( pxx[0][j][k] ); 
		}
	}

	for (i=0;i<nx-1;i++) {
		//		printf("%d\n",i); 
		for (j=ypos-ydelta;j<=ypos+ydelta;j++) {
			for (k=xpos-xdelta;k<=xpos+xdelta;k++) {
				deltaPlocal[i+1] = 0.;
				deltaPlocal[i+1] = ( pzz[i+1][j][k] ) - ( pxx[i+1][j][k] );
				l++;
				//printf("%d %d %f\n",i,l, deltaPlocal[i]); 
				sigma1 = sigma1 + deltaPlocal[i];
			}
		}	
	}

	*sigma = sigma1 / 2.0;
	
	//printf("%f\n",*sigma);

	return(0);
}




void
getCapillaryLaplaceSigma(FLOATNUM ***pxx,
			 FLOATNUM ***pyy,
			 FLOATNUM ***pzz,
			 FLOATNUM *sigma,
			 int lowerZ,
			 int higherZ,
			 int offset,
			 FLOATNUM height,
			 FLOATNUM theta)
{

	FLOATNUM	pDrop, pEnv;
	int		dropPointCount = 0;
	int		envPointCount = 0;
	FLOATNUM	deltaPdropEnv;

	int		i, j, k;
	int		dj, dk;

	j=ny/2;
	k=nz/2;
	for (i = 0; i < nx - 1; i++) {
		//		if ( ( i > lowerZ + offset ) && ( i < higherZ - offset )) {
		if (( i < higherZ - offset ) && (i > offset) ) {
			//			pDrop = pDrop + ( sqrt( pow(pxx[i][j][k],2) + pow(pyy[i][j][k],2) + pow(pzz[i][j][k],2) ));
			pDrop = pDrop + fabs((pxx[i][j][k] + pyy[i][j][k] + pzz[i][j][k])/3);
			//			for (dj=-5;dj<5;dj++) {
			//				for (dk=-5;dk<5;dk++) {
			//				pDrop = pDrop + fabs((pxx[i][j+dj][k+dk] + pyy[i][j+dj][k+dk] + pzz[i][j+dj][k+dk])/3);
							dropPointCount = dropPointCount + 1;
				//				}
				//			}
			printf("Dropi:\t%u\t%f\t%u\n", i, pDrop, dropPointCount);
		}
		//else if ( ( i < lowerZ - offset ) || ( i > higherZ + offset )) {
		else if (( i > higherZ + offset )  && (i < nx - offset) ) {
			//			pEnv = pEnv + ( sqrt( pow(pxx[i][j][k],2) + pow(pyy[i][j][k],2) + pow(pzz[i][j][k],2) ));
			pEnv = pEnv + fabs((pxx[i][j][k] + pyy[i][j][k] + pzz[i][j][k])/3);
			//			for (dj=-5;dj<5;dj++) {
			//				for (dk=-5;dk<5;dk++) {
			//				pEnv = pEnv + fabs((pxx[i][j+dj][k+dk] + pyy[i][j+dj][k+dk] + pzz[i][j+dj][k+dk])/3);
							envPointCount = envPointCount + 1;
			//				}
			//			}
			printf("Envi:\t%u\t%f\t%u\n", i, pEnv, envPointCount);
		}
	}
	
	printf("%f %d %f\n", pDrop, dropPointCount, ( pDrop / dropPointCount ));
	printf("%f %d %f\n", pEnv, envPointCount, ( pEnv / envPointCount ));

	deltaPdropEnv = ( pDrop / dropPointCount )  - ( pEnv / envPointCount ) ;
	if (deltaPdropEnv < 0) {deltaPdropEnv = -1 * deltaPdropEnv;}
	
		printf("%f\n",deltaPdropEnv);

	*sigma = ( deltaPdropEnv * height) / ( 2 * 0.244);//cos(theta) );
	return;
}
