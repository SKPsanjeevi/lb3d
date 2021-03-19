// getContactAngle.c
//
// Filename:	getContactAngle.c
// Authors:	Sebastian Schmieschek, A. Sarkar
// Version:	0.2
// Last change:	20.03.2008
// 
// Purpose:	Determine the contact angle of a droplet using variuos methods,
//		This is mostly obsolete, since the methods were put into separate programs


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

	struct geometry_t CAresult;
	struct sigma_t mechSigmaResult;
	struct sigma_t mechRadialSigmaResult;
	struct laplace_t laplaceResult;
	struct sigma_t benziResult;
	struct sigma_t huangResult;
	//struct geometry_t huangResult;

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
	gBR = atof(argv[2]);
	rockC = atof(argv[3]);

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

 

	// Calculate geometrical values
	//printf("\n*** Geometry ***\n");
	calcGeometricalContactAngle(colour, 5, 0.01, &CAresult);

	invRadius = 1/CAresult.radius;

	printf("AboveGround:\t\t%fheigth:\t\t%f\nbase:\t\t%f\nradius:\t\t%f\ntheta:\t\t\t\t%f\n",
	       CAresult.aboveGround, CAresult.height, CAresult.base, CAresult.radius, CAresult.theta);
    
	//    getColourShift(colour, CAresult.aboveGround, CAresult.radius, CAresult.heigth);
	//printf("\n*** Mechanical Definition of surface Tension ***\n");

	calcSigmaMech(pxx, pzz, CAresult.aboveGround, CAresult.radius, CAresult.height, &mechSigmaResult);

	calcSigmaMechRadial(pxx, pxy, pyz, pxz, pzz, CAresult.aboveGround, CAresult.radius, CAresult.height, &mechRadialSigmaResult);
    
	//printf("\n*** Laplace Law ***\n");
	calcSigmaLaplace(pxx, pyy, pzz, CAresult.aboveGround, CAresult.radius, CAresult.height, &laplaceResult);

	calcSigmaLaplace2(colour, pxx, pyy, pzz, CAresult.aboveGround, CAresult.radius, CAresult.height, &laplaceResult);
    
	/*     printf("\n*** Benzi et al ***\n"); */
	/*     if (calcSigmaBenzi(redDens, blueDens, CAresult.aboveGround, CAresult.radius, CAresult.heigth, &benziResult ) == 0) { */

	/*     } */

	/*     printf("\n*** Benzi et al using Psi = 1 - exp(-n) ***\n"); */
	if (calcSigmaBenzi(redDens, blueDens, CAresult.aboveGround, CAresult.radius, CAresult.height, &benziResult ) == 0) {
	  printf("BenziOK\n");
	} 
	//printf("\n*** Huang et al ***\n");
	calcSigmaHuang(redDens, blueDens, gBR, rockC, CAresult.aboveGround, CAresult.radius, CAresult.height, &huangResult );

	//printf("geo\t\tmech\t\tmechR\t\thuang\n");
	printf("%f\t%f\t%f\t%f\t%f\n", CAresult.theta, mechSigmaResult.theta, mechRadialSigmaResult.theta, huangResult.theta, benziResult.sigmaSL);
	//printf("geo\t\tmechR\t\thuang\n");
	//printf("%f\n%f\t%f\n%f\t%f\n", CAresult.theta, mechRadialSigmaResult.theta, mechRadialSigmaResult.sigmaLG, huangResult.theta, huangResult.sigmaLG);
	//printf("%f\t%f\t%f\t%f\n", CAresult.radius, laplaceResult.sigmaLG, invRadius, laplaceResult.deltaP);

}

getColourShift(FLOATNUM ***colour,
	       FLOATNUM offset,
	       FLOATNUM radius,
	       FLOATNUM heigth)
{
	int i,j,k,l, dropCenterZ;
	FLOATNUM virtPos;
	FLOATNUM lowY, highY, lowZ, highZ;
    
	dropCenterZ = (int) (offset + heigth - radius);
	k = nx / 2;
 
	for (l = 0; l < 4; l++) {
		for (i = 0; i < nz-1; i++) {
			for (j = 0; j < ny-1; j++) {
				if (i < dropCenterZ && j >= 64){// || i < dropCenterZ && j < 64) {
					if ( (colour[i][j][k] <= 0.) && (colour[i+1][j][k] > 0.) ) {
						virtPos = i - (colour[i][j][k] ) / (colour[i+1][j][k] - colour[i][j][k]);
						lowZ = virtPos;
						if ( l == 3 && i < dropCenterZ && j >= 64 && j < 78) printf("%u %f\n",j,lowZ);
						if ( l == 0 && i < dropCenterZ && j < 64 && j>= 50) printf("%u %f\n",j,lowZ);
					}
				}
				if (i >= dropCenterZ && j >= 64){// || i >= dropCenterZ && j < 64) {
					if ( (colour[i][j][k] >= 0.) && (colour[i+1][j][k] < 0.) ) {
						virtPos = (i + 1) + (colour[i+1][j][k] ) / ( colour[i][j][k] - colour[i+1][j][k]);
						highZ = virtPos;
						if ( l == 2 && i >= dropCenterZ && j >= 64) printf("%u %f\n",j,highZ);
						if ( l == 1 && i >= dropCenterZ && j < 64) printf("%u %f\n",j,highZ);
	
					}
				}
				if (j < 50) {
					if ( (colour[i][j][k] <= 0.) && (colour[i][j+1][k] > 0.) ) {		
						virtPos = j - (colour[i][j][k] ) / ( colour[i][j+1][k] - colour[i][j][k] );
						lowY = virtPos;
						if (l == 0 ) printf("%f %u\n",lowY,i);
					}
				}
				if (j >= 78) {
					if ( (colour[i][j][k] >= 0.) && (colour[i][j+1][k] < 0.) ) {
						virtPos = (j + 1) + ( colour[i][j+1][k] ) / ( colour[i][j][k] - colour[i][j+1][k]);
						highY = virtPos;
						if (l== 3) printf("%f %u\n",highY,i);
					}
				}
			}
		}
	}
	return(0);
}

/* int getZeroColourZ(FLOATNUM ***colour, */
/* 		 int xpos, */
/* 		 int ypos, */
/* 		 FLOATNUM *result) */
/* { */


/* 	int		i,j1,k1,j2,k2; */
/* 	FLOATNUM	z1; */
/* 	FLOATNUM	z2; */


/* 	for (i=0;i<=nz-1;i++) { */
		
/* 		if (colour[i][ypos][xpos] <= 0. && colour[i+1][ypos][xpos] > 0.) { */
/* 			z1 = i - (colour[i][ypos][xpos] ) / (colour[i+1][ypos][xpos] - colour[i][ypos][xpos]); */
/* 		} */
/* 		else if (colour[i][ypos][xpos] >= 0. && colour[i+1][ypos][xpos] < 0.) { */
/* 			z1 = (i + 1) + (colour[i+1][ypos][xpos] ) / ( colour[i][ypos][xpos] - colour[i+1][ypos][xpos]); */
/* 		} */

/* } */



int calcGeometricalContactAngle(FLOATNUM ***colour,
				int tmpHeigthOffset,
				FLOATNUM colourThreshold,
				struct geometry_t *result )
{

	int		i,j,k;

	FLOATNUM	virtPos = 0.;
	FLOATNUM	lowZ = 0.;
	FLOATNUM	highZ = 0.;
	FLOATNUM	lowY = 0.;
	FLOATNUM	highY = 0.;
	FLOATNUM	aboveGround = 0.;

	FLOATNUM	tmpHeigth = 0.;
	FLOATNUM	tmpBase = 0.;
	FLOATNUM	heigth = 0.;
	FLOATNUM	base = 0.;
	FLOATNUM	diam = 0.;
	FLOATNUM	radius = 0.;
	FLOATNUM	theta = 0.;
	FLOATNUM	deltaColour = 0.;
    


	j = ny / 2;
	k = nx / 2;
	for (i = 0; i < nz-3; i++) {
		//if (i > 1) {
		deltaColour = colour[i+1][j][k] - colour[i][j][k];
			//}
		//printf("%f\t%f\n",colour[i][j][k],deltaColour);
		if ( ((colour[i][j][k] <= 0.) && (colour[i+1][j][k] > 0.))) { //)&& (colourThreshold == 0.)) 
		      //|| ( deltaColour > colourThreshold ) ) {
						virtPos = i - (colour[i][j][k] ) / (colour[i+1][j][k] - colour[i][j][k]);
						lowZ = virtPos;
			//lowZ = i;
		}
		if ( ((colour[i][j][k] >= 0.) && (colour[i+1][j][k] < 0.))) { // && (colourThreshold == 0.))
			//   || ( deltaColour < -colourThreshold ) ) {
						virtPos = (i + 1) + (colour[i+1][j][k] ) / ( colour[i][j][k] - colour[i+1][j][k]);
						highZ = virtPos;
			//highZ = i;
		}
		//printf("lowZhighZ:%f\t%f\n",lowZ,highZ);
	}
    
 
	if (lowZ != 1.) {
		aboveGround = lowZ;
	}
	else {
		// Colour change from zero is found in rock-area, assuming the drop has contact
		// and set its base to the first non-rock z value.
		lowZ = 2.;
	}
	tmpHeigth = highZ - lowZ - tmpHeigthOffset; 
   
	i = 7;
	for (j = 1; j < ny-1; j++) {
		//if (j > 1) {
		deltaColour = colour[i][j][k] - colour[i][j-1][k];
			//}
		//printf("%f\t%f\n",colour[i][j][k],deltaColour);
		if ( ((colour[i][j][k] <= 0.) && (colour[i][j+1][k] > 0.) )) {//&& (colourThreshold == 0.)) 
		     //     || ( deltaColour > colourThreshold ) ) {
			virtPos = j - (colour[i][j][k] ) / ( colour[i][j+1][k] - colour[i][j][k] );
			lowY = virtPos;
		}
		     if ( ((colour[i][j][k] >= 0.) && (colour[i][j+1][k] < 0.) )) { // && (colourThreshold == 0.)) 
		     //    || ( deltaColour < -colourThreshold ) ) {
			virtPos = (j + 1) + ( colour[i][j+1][k] ) / ( colour[i][j][k] - colour[i][j+1][k]);
			highY = virtPos;
		}
		     //printf("%f\t%f\n",lowY,highY);
	}
	tmpBase = highY - lowY;

	radius = ((4*(pow(tmpHeigth,2)))+(pow(tmpBase,2)))/(8 * tmpHeigth);

	heigth = tmpHeigth + tmpHeigthOffset;
	base = 2* sqrt(pow(radius,2)-pow((heigth-radius),2));
	theta = atan2((base / 2), (radius - heigth)) * (180 / pi);
    

	result->aboveGround = aboveGround;
	result->height = heigth;
	result->base = base;
	result->radius = radius;
	result->theta = theta;

	return(0);
}

int calcGeometricContactAngle(FLOATNUM base,
			      FLOATNUM height,
			      FLOATNUM *theta)
{
	FLOATNUM	radius = 0.;
	//	FLOATNUM	theta = 0.;

	radius = ((4*(pow(height,2)))+(pow(base,2)))/(8 * height);
	*theta = atan2((base / 2), (radius - height)) * (180 / pi);

	return(0);
}


int calcSigmaMech(FLOATNUM ***pxx,
		  FLOATNUM ***pzz,
		  FLOATNUM offset,
		  FLOATNUM radius,
		  FLOATNUM heigth,
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
    

	dropCenterZ = (int) (offset + heigth - radius);
	//printf("centerZ:\t%u\n",dropCenterZ);
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
			//if (offset == 0) {

				/* 		if((i == 0) || (i == (9))) factor = 3./8.; */
				/* 		else if((i == 1) || (i == (8))) factor = 7./6.; */
				/* 		else if((i == 2) || (i == (7))) factor = 23./24.; */
				/* 		else factor = 1.; */
		
				sigmaSL = sigmaSL + (factor * deltaPlocal1[i]);
				sigmaSG = sigmaSG + (factor * deltaPlocal2[i]);
		
				//				/* 		if (deltaPlocal1[i] < deltaPlocal1[i+1]) { */
				/* 		    sigmaSL = sigmaSL + ( ( ( deltaPlocal1[i+1] - deltaPlocal1[i] ) / 2 ) + deltaPlocal1[i+1] ); */
				/* 		} else if (deltaPlocal1[i] > deltaPlocal1[i+1]) { */
				/* 		    sigmaSL = sigmaSL + ( ( ( deltaPlocal1[i] - deltaPlocal1[i+1] ) / 2 ) + deltaPlocal1[i] ); */
				/* 		} */
				/* 		if (deltaPlocal2[i] < deltaPlocal2[i+1]) { */
				/* 		    sigmaSG = sigmaSG + ( ( ( deltaPlocal2[i+1] - deltaPlocal2[i] ) / 2 ) + deltaPlocal2[i+1] ); */
				/* 		} else if (deltaPlocal2[i] > deltaPlocal2[i+1]) { */
				/* 		    sigmaSG = sigmaSG + ( ( ( deltaPlocal2[i] - deltaPlocal2[i+1] ) / 2 ) + deltaPlocal2[i] ); */
				/* 		} */
				/* 	    } /\* else { *\/ */
				/* 		//sigmaSG = sigmaSG + deltaPlocal[i]; */
				/* 		if (deltaPlocal[i] < deltaPlocal[i+1]) { */
				/* 		    sigmaSG = sigmaSG + ( ( ( deltaPlocal1[i+1] - deltaPlocal1[i] ) / 2 ) + deltaPlocal1[i+1] ); */
				/* 		} else if (deltaPlocal[i] > deltaPlocal[i+1]) { */
				/* 		    sigmaSG = sigmaSG + ( ( ( deltaPlocal1[i] - deltaPlocal1[i+1] ) / 2 ) + deltaPlocal1[i] ); */
				/* 		} */
				/* 	    } */
				//}
	    
		}
		else if (i > (int) (offset + heigth - 15) && i < (int) (offset + heigth + 15)) {
	    
			/* 	    if((i == (int) (offset + heigth - 9)) || (i == (int) (offset + heigth + 9))) factor = 3./8.; */
			/* 	    else if((i == (int) (offset + heigth - 8)) || (i == (int) (offset + heigth + 8))) factor = 7./6.; */
			/* 	    else if((i == (int) (offset + heigth - 7)) || (i == (int) (offset + heigth + 7))) factor = 23./24.; */
			/* 	    else factor = 1.; */
			sigmaLG = sigmaLG + (factor * deltaPlocal1[i]); 
			/* 	    if (deltaPlocal1[i] < deltaPlocal1[i+1]) { */
			/* 		sigmaLG = sigmaLG + ( ( ( deltaPlocal1[i+1] - deltaPlocal1[i] ) / 2 ) + deltaPlocal1[i+1] ); */
			/* 	    } else if (deltaPlocal1[i] > deltaPlocal1[i+1]) { */
			/* 		sigmaLG = sigmaLG + ( ( ( deltaPlocal1[i] - deltaPlocal1[i+1] ) / 2 ) + deltaPlocal1[i] ); */
			/* 	    } */
		}
	
	}
	theta = acos( ( sigmaSG - sigmaSL ) / sigmaLG ) * ( 180 / pi);
	result->sigmaSL = sigmaSL;
	result->sigmaSL = sigmaSG;
	result->sigmaSL = sigmaLG;
	result->theta = theta;

	//	printf("sigmaSL:\t%f\nsigmaSG:\t%f\nsigmaLG:\t%f\ntheta:\t\t\t\t%f\n",sigmaSL,sigmaSG,sigmaLG, theta);
    
    

	return(0);

}


int calcSigmaMechRadial(FLOATNUM ***pxx,
			FLOATNUM ***pxy,
			FLOATNUM ***pyz,
			FLOATNUM ***pxz,
			FLOATNUM ***pzz,
			FLOATNUM offset,
			FLOATNUM radius,
			FLOATNUM heigth,
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

	dropCenterZ = (int) (offset + heigth - radius);
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
			//printf("%d\t%f\t%f\n",i,deltaPlocal1[i],deltaPlocal2[i]);	

			if (i < 8) {

/* 				sigmaSL = sigmaSL + ( deltaPlocal1[i] ); */
/* 				SLvalCount++; */
/* 				sigmaSG = sigmaSG + ( deltaPlocal2[i] ); */
/* 				SGvalCount++; */
				//if (offset == 0) {
				if (deltaPlocal1[i] < deltaPlocal1[i+1]) {
					sigmaSL = sigmaSL + ( ( ( deltaPlocal1[i+1] - deltaPlocal1[i] ) / 2 ) + deltaPlocal1[i+1] );
					SLvalCount++;
				} else if (deltaPlocal1[i] > deltaPlocal1[i+1]) {
					sigmaSL = sigmaSL + ( ( ( deltaPlocal1[i] - deltaPlocal1[i+1] ) / 2 ) + deltaPlocal1[i] );
					SLvalCount++;
				}
				if (deltaPlocal2[i] < deltaPlocal2[i+1]) {
					sigmaSG = sigmaSG + ( ( ( deltaPlocal2[i+1] - deltaPlocal2[i] ) / 2 ) + deltaPlocal2[i+1] );
					SGvalCount++;
				} else if (deltaPlocal2[i] > deltaPlocal2[i+1]) {
					sigmaSG = sigmaSG + ( ( ( deltaPlocal2[i] - deltaPlocal2[i+1] ) / 2 ) + deltaPlocal2[i] );
					SGvalCount++;
				}
				//} /* else { */
				/* 		//sigmaSG = sigmaSG + deltaPlocal[i]; */
			/* 		if (deltaPlocal[i] < deltaPlocal[i+1]) { */
			/* 		    sigmaSG = sigmaSG + ( ( ( deltaPlocal1[i+1] - deltaPlocal1[i] ) / 2 ) + deltaPlocal1[i+1] ); */
			/* 		} else if (deltaPlocal[i] > deltaPlocal[i+1]) { */
			/* 		    sigmaSG = sigmaSG + ( ( ( deltaPlocal1[i] - deltaPlocal1[i+1] ) / 2 ) + deltaPlocal1[i] ); */
			/* 		} */
			/* 	    } */
				
			}
			else if (i > (int) (offset + heigth - 10) && i < (int) (offset + heigth + 10)) {
				/* sigmaLG = sigmaLG + ( deltaPlocal1[i]  * (i/(int) (offset+heigth)) ); */
/* 				sigmaLG2 = sigmaLG2 + ( deltaPlocal1[i]  * pow((i/(int) (offset+heigth)),2) ); */
/* 				LGvalCount++; */
				if (deltaPlocal1[i] < deltaPlocal1[i+1]) {
					sigmaLG = sigmaLG + ( ( ( ( deltaPlocal1[i+1] - deltaPlocal1[i] ) / 2 ) + deltaPlocal1[i+1] ) * (i/(int) (offset+heigth)) );
					LGvalCount++;
				} else if (deltaPlocal1[i] > deltaPlocal1[i+1]) {
					sigmaLG = sigmaLG + ( ( ( ( deltaPlocal1[i] - deltaPlocal1[i+1] ) / 2 ) + deltaPlocal1[i] ) * (i/((int) (offset+heigth))) );
					LGvalCount++;
				}
			}
			
		}
	}
	}
	}
	}
	sigmaSL = sigmaSL ;// 25600.;// SLvalCount;
       	sigmaSG = sigmaSG ;// 25600.;// SGvalCount;
	sigmaLG = sigmaLG ;// 25600.;// LGvalCount;
       

	theta = acos( (sigmaSG -sigmaSL ) / sigmaLG ) * ( 180 / pi);
	result->sigmaSL = sigmaSL;
	result->sigmaSG = sigmaSG;
	result->sigmaLG = sigmaLG;
	result->theta = theta;
	printf("sigmaSL:\t%f\nsigmaSG:\t%f\nsigmaLG:\t%f\ntheta:\t\t\t\t%f\n",sigmaSL,sigmaSG,sigmaLG, theta);
	
		
		
	return(0);

}

int calcSigmaLaplace(FLOATNUM ***pxx,
		     FLOATNUM ***pyy,
		     FLOATNUM ***pzz,
		     FLOATNUM offset,
		     FLOATNUM radius,
		     FLOATNUM heigth,
		     struct laplace_t *result )
{

	int		i, j, k;
	int		dropPointCount = 0;
	int		envPointCount = 0;

	FLOATNUM	deltaPdropEnv = 0.;
	FLOATNUM	pDrop = 0.;
	FLOATNUM	pEnv = 0.;
	FLOATNUM	sigma = 0.;


	j = ny / 2;
	k = nx / 2;
    
	for (i = (int) ( offset + 10 ); i < ( nz - 10 ); i++) {//for (i=1;i<64;i++) {
		if ( i < (int) (offset + heigth - 10) ) {//	if (i == nz/2) {

			pDrop = pDrop + ( ( pxx[i][j][k] + pyy[i][j][k] + pzz[i][j][k] ) / 3 );
			dropPointCount = dropPointCount + 1;
			//printf("Dropi:\t%u\t%f\t%u\n", i, pDrop, dropPointCount);
		}
		else if ( i > (int) ( offset + heigth + 10 ) ) {//else if (i == 63) {

			pEnv = pEnv + ( ( pxx[i][j][k] + pyy[i][j][k] + pzz[i][j][k] ) / 3 );
			envPointCount = envPointCount + 1;
			//printf("Envi:\t%u\t%f\t%u\n", i, pEnv, envPointCount);
		}
	}

	deltaPdropEnv = ( pDrop / dropPointCount )  - ( pEnv / envPointCount ) ;
	if (deltaPdropEnv < 0) {deltaPdropEnv = -1 * deltaPdropEnv;}

	sigma = deltaPdropEnv / ( 2 / radius );
	result->sigmaLG = sigma;
	result->deltaP = deltaPdropEnv;
	//printf("deltaPdropEnv:\t%f\nsigmaLG:\t%f\n", deltaPdropEnv, sigma);
	return(0);

}

int calcSigmaLaplace2(FLOATNUM ***colour,
		      FLOATNUM ***pxx,
		      FLOATNUM ***pyy,
		      FLOATNUM ***pzz,
		      FLOATNUM offset,
		      FLOATNUM radius,
		      FLOATNUM heigth,
		      struct laplace_t *result )
{

	int		i, j, k;
	int		dropPointCount = 0;
	int		envPointCount = 0;

	int		myOffset = 5;
	int		lowerZ, higherZ;

	FLOATNUM	deltaPdropEnv = 0.;
	FLOATNUM	pDrop = 0.;
	FLOATNUM	pEnv = 0.;
	FLOATNUM	sigma = 0.;

	FLOATNUM	invRadius = 0.;

	j = ny / 2;
	k = nx / 2;
    
	for (i = 0; i < nz - 2; i++) {
		if ( ((colour[i][j][k] <= 0.) && (colour[i+1][j][k] > 0.))) {
			lowerZ = i;
		}
		if ( ((colour[i][j][k] > 0.) && (colour[i+1][j][k] <= 0.))) {
			higherZ = i;
		}
	}
	printf("%d\t%d\n",lowerZ,higherZ);
	for (i = 0; i < nz - 1; i++) {
		if ( ( i > lowerZ + myOffset ) && ( i < higherZ - myOffset )) {
			pDrop = pDrop + ( sqrt( pow(pxx[i][j][k],2) + pow(pyy[i][j][k],2) + pow(pzz[i][j][k],2) ));
			dropPointCount = dropPointCount + 1;
			//printf("Dropi:\t%u\t%f\t%u\n", i, pDrop, dropPointCount);
		}
		else if ( ( i < lowerZ - myOffset ) || ( i > higherZ + myOffset )) {
			pEnv = pEnv + ( sqrt( pow(pxx[i][j][k],2) + pow(pyy[i][j][k],2) + pow(pzz[i][j][k],2) ));
			envPointCount = envPointCount + 1;
			//printf("Envi:\t%u\t%f\t%u\n", i, pEnv, envPointCount);
		}
	}

	deltaPdropEnv = ( pDrop / dropPointCount )  - ( pEnv / envPointCount ) ;
	if (deltaPdropEnv < 0) {deltaPdropEnv = -1 * deltaPdropEnv;}

	sigma = deltaPdropEnv / ( 2 / radius );
	result->sigmaLG = sigma;
	result->deltaP = deltaPdropEnv;

	invRadius = 1/radius;

	printf("deltaPdropEnv:\t%f\nsigmaLG:\t%f\t%f\n", deltaPdropEnv, sigma, invRadius);
	return(0);

}

calcSigmaBenziPsi(FLOATNUM ***redDens,
		  FLOATNUM ***blueDens,
		  FLOATNUM offset,
		  FLOATNUM radius,
		  FLOATNUM heigth,
		  struct sigma_t *result )
{
	int         i, j1, k1, j2, k2;

	FLOATNUM	sigmaSLred = 0.;
	FLOATNUM	sigmaSGred = 0.;
	FLOATNUM	sigmaLGred = 0.;
	FLOATNUM	sigmaSLblue = 0.;
	FLOATNUM	sigmaSGblue = 0.;
	FLOATNUM	sigmaLGblue = 0.;
	FLOATNUM	sigmaSL = 0.;
	FLOATNUM	sigmaSG = 0.;
	FLOATNUM	sigmaLG = 0.;

	FLOATNUM	psiRed1[nz];
	FLOATNUM	psiBlue1[nz];
	FLOATNUM	psiRed2[nz];
	FLOATNUM	psiBlue2[nz];

	FLOATNUM	deltaRedDensLocal1[nz];
	FLOATNUM	deltaBlueDensLocal1[nz];
	FLOATNUM	deltaRedDensLocal2[nz];
	FLOATNUM	deltaBlueDensLocal2[nz];

	FLOATNUM	thetaRed;
	FLOATNUM	thetaBlue;
	FLOATNUM	theta;
    
	j1 = ny / 2;
	k1 = nx / 2;
	j2 = 5;
	k2 = 5;
	psiRed1[0] = 1 - exp( -redDens[0][j1][k1] );
	psiRed2[0] = 1 - exp( -redDens[0][j2][k2] );
	psiBlue1[0] = 1 - exp( -blueDens[0][j1][k1] );
	psiBlue2[0] = 1 - exp( -blueDens[0][j2][k2] );
	for (i = 0; i < nz - 1; i++) {

		deltaRedDensLocal1[i] = 0.;
		deltaBlueDensLocal1[i] = 0.;

		psiRed1[i+1] = 1 - exp( -redDens[i+1][j1][k1] );
		psiRed2[i+1] = 1 - exp( -redDens[i+1][j2][k2] );
		psiBlue1[i+1] = 1 - exp( -blueDens[i+1][j1][k1] );
		psiBlue2[i+1] = 1 - exp( -blueDens[i+1][j2][k2] );

		deltaRedDensLocal1[i] = psiRed1[i] - psiRed1[i+1];
		deltaRedDensLocal2[i] = psiRed2[i] - psiRed2[i+1];
		deltaBlueDensLocal1[i] = psiBlue1[i] - psiBlue1[i+1];
		deltaBlueDensLocal2[i] = psiBlue2[i] - psiBlue2[i+1];

		if ( i < 10  && i > 1) {
			if (offset == 0) {
				sigmaSLred = sigmaSLred + (deltaRedDensLocal1[i] * deltaRedDensLocal1[i]);
				sigmaSGblue = sigmaSGblue + (deltaBlueDensLocal1[i] * deltaBlueDensLocal1[i]);
				sigmaSGred = sigmaSGred + (deltaRedDensLocal2[i] * deltaRedDensLocal2[i]);
				sigmaSLblue = sigmaSLblue + (deltaBlueDensLocal2[i] * deltaBlueDensLocal2[i]);
				sigmaSL = sigmaSL + (sqrt(pow(deltaRedDensLocal1[i],2)) * sqrt(pow(deltaBlueDensLocal1[i],2)));
				sigmaSG = sigmaSG + (sqrt(pow(deltaRedDensLocal2[i],2)) * sqrt(pow(deltaBlueDensLocal2[i],2)));
			} else {
				sigmaSGred = sigmaSGred + (deltaRedDensLocal1[i] * deltaRedDensLocal1[i]);
				sigmaSLblue = sigmaSLblue + (deltaBlueDensLocal1[i] * deltaBlueDensLocal1[i]);
			}
		}
		else if (i > (int) (offset + heigth - 10) && i < (int) (offset + heigth + 10)) {
			sigmaLGred = sigmaLGred + (deltaRedDensLocal1[i] * deltaRedDensLocal1[i]);
			sigmaLGblue = sigmaLGblue + (deltaBlueDensLocal1[i] * deltaBlueDensLocal1[i]);
			sigmaLG = sigmaLG + (sqrt(pow(deltaRedDensLocal1[i],2)) * sqrt(pow(deltaBlueDensLocal1[i],2)));
		}
	}

	thetaRed = acos( ( -sigmaSLred - sigmaSGred ) / sigmaLGred ) * ( 180 / pi);
	thetaBlue = acos( ( -sigmaSLblue - sigmaSGblue ) / sigmaLGblue ) * ( 180 / pi);
	theta = acos( ( -sigmaSL - sigmaSG ) / sigmaLG ) * ( 180 / pi);

	//printf("sigmaSLR:\t%f\nsigmaSGR:\t%f\nsigmaLGR:\t%f\nthetaR:\t%f\nsigmaSLB:\t%f\nsigmaSGB:\t%f\nsigmaLGB:\t%f\nthetaB:\t%f\nsigmaSL:\t%f\nsigmaSG:\t%f\nsigmaLG:\t%f\ntheta:\t\t\t\t%f\n", sigmaSLred, sigmaSGred, sigmaLGred, thetaRed, sigmaSLblue, sigmaSGblue, sigmaLGblue, thetaBlue, sigmaSL, sigmaSG, sigmaLG, theta);
	result->sigmaSL = sigmaSL;
	result->sigmaSL = sigmaSG;
	result->sigmaSL = sigmaLG;
	result->theta = theta;
	return(0);
}

calcSigmaBenzi(FLOATNUM ***redDens,
	       FLOATNUM ***blueDens,
	       FLOATNUM offset,
	       FLOATNUM radius,
	       FLOATNUM heigth,
	       struct sigma_t *result )
{
	int         i, j1, k1, j2, k2;

	FLOATNUM	sigmaSLred = 0.;
	FLOATNUM	sigmaSGred = 0.;
	FLOATNUM	sigmaLGred = 0.;
	FLOATNUM	sigmaSLblue = 0.;
	FLOATNUM	sigmaSGblue = 0.;
	FLOATNUM	sigmaLGblue = 0.;
	FLOATNUM	sigmaSL = 0.;
	FLOATNUM	sigmaSG = 0.;
	FLOATNUM	sigmaLG = 0.;

	FLOATNUM	deltaRedDensLocal1[nz];
	FLOATNUM	deltaBlueDensLocal1[nz];
	FLOATNUM	deltaRedDensLocal2[nz];
	FLOATNUM	deltaBlueDensLocal2[nz];

	FLOATNUM	thetaRed;
	FLOATNUM	thetaBlue;
	FLOATNUM	theta;
    
	j1 = ny / 2;
	k1 = nx / 2;
	j2 = 5;
	k2 = 5;

	for (i = 0; i < nz - 1; i++) {

		deltaRedDensLocal1[i] = 0.;
		deltaBlueDensLocal1[i] = 0.;

		deltaRedDensLocal1[i] = redDens[i][j1][k1] - redDens[i+1][j1][k1];
		deltaBlueDensLocal1[i] = blueDens[i][j1][k1] - blueDens[i+1][j1][k1];
		deltaRedDensLocal2[i] = redDens[i][j2][k2] - redDens[i+1][j2][k2];
		deltaBlueDensLocal2[i] = blueDens[i][j2][k2] - blueDens[i+1][j2][k2];

		if ( i < 10 && i > 1) {
			if (offset == 0) {
				sigmaSLred = sigmaSLred + (deltaRedDensLocal1[i] * deltaRedDensLocal1[i]);
				sigmaSGblue = sigmaSGblue + (deltaBlueDensLocal1[i] * deltaBlueDensLocal1[i]);
				sigmaSGred = sigmaSGred + (deltaRedDensLocal2[i] * deltaRedDensLocal2[i]);
				sigmaSLblue = sigmaSLblue + (deltaBlueDensLocal2[i] * deltaBlueDensLocal2[i]);
				sigmaSL = sigmaSL + (sqrt(pow(deltaRedDensLocal1[i],2)) * sqrt(pow(deltaBlueDensLocal1[i],2)));
				sigmaSG = sigmaSG + (sqrt(pow(deltaRedDensLocal2[i],2)) * sqrt(pow(deltaBlueDensLocal2[i],2)));
			} else {
				sigmaSGred = sigmaSGred + (deltaRedDensLocal1[i] * deltaRedDensLocal1[i]);
				sigmaSLblue = sigmaSLblue + (deltaBlueDensLocal1[i] * deltaBlueDensLocal1[i]);
			}
		}
		else if (i > (int) (offset + heigth - 10) && i < (int) (offset + heigth + 10)) {
			sigmaLGred = sigmaLGred + (deltaRedDensLocal1[i] * deltaRedDensLocal1[i]);
			sigmaLGblue = sigmaLGblue + (deltaBlueDensLocal1[i] * deltaBlueDensLocal1[i]);
			sigmaLG = sigmaLG + (sqrt(pow(deltaRedDensLocal1[i],2)) * sqrt(pow(deltaBlueDensLocal1[i],2)));
		}
	}

	thetaRed = acos( ( -sigmaSLred - sigmaSGred ) / sigmaLGred ) * ( 180 / pi);
	thetaBlue = acos( ( -sigmaSLblue - sigmaSGblue ) / sigmaLGblue ) * ( 180 / pi);
	theta = acos( ( -sigmaSL - sigmaSG ) / sigmaLG ) * ( 180 / pi);
	result->sigmaSL = sigmaSL;
	result->sigmaSL = sigmaSG;
	result->sigmaSL = sigmaLG;
	result->theta = theta;
	//printf("sigmaSLR:\t%f\nsigmaSGR:\t%f\nsigmaLGR:\t%f\nthetaR:\t%f\nsigmaSLB:\t%f\nsigmaSGB:\t%f\nsigmaLGB:\t%f\nthetaB:\t%f\nsigmaSL:\t%f\nsigmaSG:\t%f\nsigmaLG:\t%f\ntheta:\t\t\t\t%f\n", sigmaSLred, sigmaSGred, sigmaLGred, thetaRed, sigmaSLblue, sigmaSGblue, sigmaLGblue, thetaBlue, sigmaSL, sigmaSG, sigmaLG, theta);

	return(0);
}

calcSigmaHuangPsi(FLOATNUM ***redDens,
		  FLOATNUM ***blueDens,
		  FLOATNUM offset,
		  FLOATNUM radius,
		  FLOATNUM heigth,
		  struct sigma_t *result )
{
	int         i, j, k;
    
	FLOATNUM	sigmaSLp = 0.;
	FLOATNUM	sigmaSGp = 0.;
	FLOATNUM	sigmaLGp = 0.;
 
	FLOATNUM	maxRedDens = 0.;
	FLOATNUM	minRedDens = 1.;

	FLOATNUM	maxBlueDens = 0.;
	FLOATNUM	minBlueDens = 1.;

	FLOATNUM	Gbr = 0.14;
	FLOATNUM	Grock = -0.3 * Gbr;
	//    printf("MEIN GBR IST %f\nMAIN ROCK ist %f\n",Gbr,Grock);

	FLOATNUM	thetap;

	FLOATNUM	psiRedMin;
	FLOATNUM	psiRedMax;

	j = ny / 2;
	k = nx / 2;

	for (i = 2; i < nz; i++) {
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
	psiRedMin = 1 - exp( -minRedDens );
	psiRedMax = 1 - exp( -maxRedDens );

	sigmaSLp = Grock * psiRedMax;
	sigmaSGp = Grock * psiRedMin;
	sigmaLGp = Gbr * ( ( psiRedMax - psiRedMin ) / 2 );

	thetap= acos( ( sigmaSLp - sigmaSGp ) / sigmaLGp ) * ( 180 / pi);
	result->sigmaSL = sigmaSLp;
	result->sigmaSL = sigmaSGp;
	result->sigmaSL = sigmaLGp;
	result->theta = thetap;
	// printf("maxRedDens:\t%f\nminRedDens:\t%f\n",maxRedDens,minRedDens);
	// printf("sigmaSL:\t%f\nsigmaSG:\t%f\nsigmaLG:\t%f\ntheta:\t\t\t\t%f\n", sigmaSLp, sigmaSGp, sigmaLGp,thetap);
}

calcSigmaHuang(FLOATNUM ***redDens,
	       FLOATNUM ***blueDens,
	       FLOATNUM gBR,
	       FLOATNUM rockC,
	       FLOATNUM offset,
	       FLOATNUM radius,
	       FLOATNUM heigth,
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

	//    FLOATNUM	Gbr = 0.18;
	FLOATNUM	gROCK = (((5. * rockC) + 14.) * gBR)/19.;
	//printf("MEIN GBR IST %f\nMAIN ROCK ist %f\n",gBR,gROCK);

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
 
 	sigmaSL = (((5. * rockC) + (14. * minRedDens) * gBR) / 19.) * maxRedDens;//gROCK * maxRedDens;//redDens[2][j][k]; *
 	sigmaSG = (((5. * rockC) + (14. * maxRedDens) * gBR) / 19.) * minRedDens;//gROCK * minRedDens;
 	sigmaLG = gBR * ( ( maxRedDens - minRedDens ) / 2 );

/* 	sigmaSL = -gROCK*(5./14.);//gROCK * maxRedDens;//redDens[2][j][k]; */
/* 	sigmaSG = 0.;// gROCK;//(((5. * rockC) + (14. * maxRedDens) * gBR) / 19.) * minRedDens;//gROCK * minRedDens; */
/* 	sigmaLG = gBR * ( ( maxRedDens - minRedDens ) / 2 ); */

	theta = acos( ( sigmaSL - sigmaSG ) / sigmaLG ) * ( 180 / pi);
	result->sigmaSL = sigmaSL;
	result->sigmaSG = sigmaSG;
	result->sigmaLG = sigmaLG;
	result->theta = theta;
	//printf("maxRedDens:\t%f\nminRedDens:\t%f\n",maxRedDens,minRedDens);
	//printf("sigmaSL:\t%f\nsigmaSG:\t%f\nsigmaLG:\t%f\ntheta:\t\t\t\t%f\n", sigmaSL, sigmaSG, sigmaLG,theta);
}


