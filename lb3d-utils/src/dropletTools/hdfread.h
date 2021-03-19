#ifndef myFileFuncs_h
#define myFileFuncs_h
#include <hdf5.h>

#ifdef SGL
typedef float FLOATNUM;
#else
typedef double FLOATNUM;
#endif

// Function declarations
int hdfread	(char *,
		 char *,
		 FLOATNUM ****);

int getDims	(char *FILE,
		 int *nx,
		 int *ny,
		 int *nz);
#endif //myFileFuncs_h
