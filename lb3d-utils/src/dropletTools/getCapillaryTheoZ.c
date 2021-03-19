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
	int		timeoff=0;
	int DEBUG = 0;

	FLOATNUM	v = 0.003;
	FLOATNUM	Ca;

	int		offset;

	filename = argv[1];
	z_0	 = atoi(argv[2]);
	maxTime  = atoi(argv[3]);
	rho	 = atof(argv[4]);
	height	 = atoi(argv[5]);
	length   = atoi(argv[6]);
	sigma = atof(argv[7]);
	theta2 = atof(argv[8]);// / (180 / pi);
	timeoff = argv[9] ? atoi(argv[9]) : 0;
	
	if (z_0 != 32) {
		offset = 32 - z_0;
		length +=offset;
	}

	if (filename) {
		getDims(filename,&nx,&ny,&nz);
		hdfread(filename, "colour",&colour);
		//hdfread2(filename, &colour);
		hdfread(filename, "pxx", &pxx);
		hdfread(filename, "pyy", &pyy);
		hdfread(filename, "pzz", &pzz);

		
		//	hdfwrite(outfilename, nx, ny, nz, colour);

		retVal = 0;
	}
	
	//	if (argc > 5) {

	//	}
			//	else {
			//		getMechDefSigma(pxx,pzz,64,4,&sigma,0,0);
			//	}


		

	

	getZeroColourZ(colour,1,4,&zPos1r,&zPos2r,0,0);
	getZeroColourZ(colour,(height/2),4,&zPos1c,&zPos2c,0,0);
	getZeroColourZ(colour,(height-1),4,&zPos1l,&zPos2l,0,0);

	

	//	printf("%f\t%f\t%f\t%f\t%f\t%f\n",zPos1r,zPos2r,zPos1c,zPos2c,zPos1l,zPos2l);


	height1 = zPos1r - zPos1c;
	height2 = zPos2r - zPos2c;
 
	//	printf("height1: %f\theight2: %f\n",height1, height2);

	//calcGeometricalRadius(127,height1,&radius1);
	//calcGeometricalRadius(127,height2,&radius2);


	//calcGeometricContactAngle(127,height1,radius1,&theta1);
	//calcGeometricContactAngle(127,height2,radius2,&theta2);
	
	//theta1=2.82743339;
	//theta2=0.314159265;
		


 	getCapillaryLaplaceSigma(pxx,pyy,pzz,&sigma2,zPos2c,zPos1c,30,height,theta2);
	//printf("theta1: %f\ttheta2: %f\n",theta1,theta2); 
	if (DEBUG == 1) {
		printf("sigma: %f\n",sigma);
	}
	
	//sigma = 0.07;
	//sigma = 0.053;
	
	mu = rho * (1./6.);
	
	getCapillaryNumber(mu,v,sigma,&Ca);

	//	Ca = 0.03

	//	theta1 = acos(cos(theta1) - (18*pow(Ca,1.2)));
	//theta2 = acos(cos(theta2) - (18*pow(Ca,1.2)));

	v_cap = sigma / mu;
	t_d = (pow(height,2) * rho) / (12 * mu);

	//	free(colour);
	//	free(pxx);
	//	free(pzz);

	
	//	theta1 = theta1;

	


	for (timestep = 0; timestep <= maxTime; timestep++) {
		z1 = ((v_cap * height * cos(theta1))/(6*length))*t_d*((exp(-timestep/t_d)+(timestep/t_d)-1)) + z_0;
		z2 = ((v_cap * height * cos(theta2))/(6*length))*t_d*((exp(-timestep/t_d)+(timestep/t_d)-1)) + z_0;
		if (timestep%100 == 0) {
			if (DEBUG ==1) {
				printf("%d\t%f\t%f\n",timestep,z1,z2);
			} else {
				printf("%d\t%f\n",timestep-timeoff,z2);
			}
		}
	}
	if (DEBUG == 1) {
		printf("theta2: %f\tcos(theta2): %f\n",theta2*(180/pi),cos(theta2));
	}
	
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

	j=ny/2;
	k=nz/2;
	for (i = 0; i < nx - 1; i++) {
		if ( ( i > lowerZ + offset ) && ( i < higherZ - offset )) {
			pDrop = pDrop + ( sqrt( pow(pxx[i][j][k],2) + pow(pyy[i][j][k],2) + pow(pzz[i][j][k],2) ));
			dropPointCount = dropPointCount + 1;
			//			printf("Dropi:\t%u\t%f\t%u\n", i, pDrop, dropPointCount);
		}
		else if ( ( i < lowerZ - offset ) || ( i > higherZ + offset )) {
			pEnv = pEnv + ( sqrt( pow(pxx[i][j][k],2) + pow(pyy[i][j][k],2) + pow(pzz[i][j][k],2) ));
			envPointCount = envPointCount + 1;
			//			printf("Envi:\t%u\t%f\t%u\n", i, pEnv, envPointCount);
		}
	}
	
	deltaPdropEnv = ( pDrop / dropPointCount )  - ( pEnv / envPointCount ) ;
	if (deltaPdropEnv < 0) {deltaPdropEnv = -1 * deltaPdropEnv;}
	
	//	printf("%f\n",deltaPdropEnv);

	*sigma = ( deltaPdropEnv * height) / ( 2 * cos(theta) );
	return;
}
