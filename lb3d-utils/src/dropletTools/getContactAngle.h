#ifndef GETCONTACTANGLE_H
#define GETCONTACTANGLE_H

struct geometry_t{
	FLOATNUM	aboveGround;
	FLOATNUM	minHeight;
	FLOATNUM	height;
	FLOATNUM	maxHeight;
	FLOATNUM	minBase;
	FLOATNUM	base;
	FLOATNUM	maxBase;
	FLOATNUM	minRadius;
	FLOATNUM	radius;
	FLOATNUM	maxRadius;
	FLOATNUM	minTheta;
	FLOATNUM	theta;
	FLOATNUM	maxTheta;
};

struct sigma_t {
	FLOATNUM	sigmaSL;
	FLOATNUM	sigmaSG;
	FLOATNUM	sigmaLG;
	FLOATNUM	theta;
};

struct laplace_t {
	FLOATNUM	sigmaLG;
	FLOATNUM	deltaP;
};

// Function declarations

int getColourShift	(FLOATNUM ***colour,
			 FLOATNUM offset,
			 FLOATNUM radius,
			 FLOATNUM heigth);

int calcGeometricalContactAngle(FLOATNUM ***colour,
				int tmpHeigthOffset,
				FLOATNUM colourThreshold,
				struct geometry_t *result );

int calcSigmaMech	(FLOATNUM ***pxx,
			 FLOATNUM ***pzz,
			 FLOATNUM offset,
			 FLOATNUM radius,
			 FLOATNUM heigth,
			 struct sigma_t *result);

int calcSigmaMechRadial	(FLOATNUM ***pxx,
			 FLOATNUM ***pxy,
			 FLOATNUM ***pyz,
			 FLOATNUM ***pxz,
			 FLOATNUM ***pzz,
			 FLOATNUM offset,
			 FLOATNUM radius,
			 FLOATNUM heigth,
			 struct sigma_t *result);

int calcSigmaLaplace	(FLOATNUM ***pxx,
			 FLOATNUM ***pyy,
			 FLOATNUM ***pzz,
			 FLOATNUM offset,
			 FLOATNUM radius,
			 FLOATNUM heigth,
			 struct laplace_t *result );

int calcSigmaLaplace2	(FLOATNUM ***colour,
			 FLOATNUM ***pxx,
			 FLOATNUM ***pyy,
			 FLOATNUM ***pzz,
			 FLOATNUM offset,
			 FLOATNUM radius,
			 FLOATNUM heigth,
			 struct laplace_t *result );

int calcSigmaBenziPsi	(FLOATNUM ***redDens,
			 FLOATNUM ***blueDens,
			 FLOATNUM offset,
			 FLOATNUM radius,
			 FLOATNUM heigth,
			 struct sigma_t *result );

int calcSigmaBenzi	(FLOATNUM ***redDens,
			 FLOATNUM ***blueDens,
			 FLOATNUM offset,
			 FLOATNUM radius,
			 FLOATNUM heigth,
			 struct sigma_t *result );

int calcSigmaHuang	(FLOATNUM ***redDens,
			 FLOATNUM ***blueDens,
			 FLOATNUM gBR,
			 FLOATNUM rockC,
			 FLOATNUM offset,
			 FLOATNUM radius,
			 FLOATNUM heigth,
			 struct sigma_t *result );


FLOATNUM	pi = 3.1415926535897932384626433832795028841971693993751058209749445923;

int		nx, ny, nz;

int		DEBUG = 0;

#endif //GETCONTACTANGLE_H
