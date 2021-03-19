#include "readbin.h"

enum datatype {
	TYPE_UNDEFINED	= -1,
	TYPE_SCALAR	= 1,
	TYPE_2SCALAR	= 2,
	TYPE_3SCALAR	= 3,
	TYPE_VECTOR	= 4
};

typedef enum datatype datatype;
datatype parse_type(char *s) ;
int read_xdr(	datatype type,	double *,double *, double *,
			unsigned int , unsigned int , unsigned int ,
			char *filename);
int write_xdr(	datatype type,	double *,double *, double *,
			unsigned int , unsigned int , unsigned int ,
			char *filename);
int read_xdr_float(	datatype type,	float *,float *, float *,
			unsigned int , unsigned int , unsigned int ,
			char *filename);
int write_xdr_float(	datatype type,	float *,float *, float *,
			unsigned int , unsigned int , unsigned int ,
			char *filename);
int read_rawbin(	datatype type,	F_REAL *,F_REAL *, F_REAL *,
			unsigned int , unsigned int , unsigned int ,
			char *filename);
int write_rawbin(	datatype type,	F_REAL *,F_REAL *, F_REAL *,
			unsigned int , unsigned int , unsigned int ,
			char *filename);
int read_bin(		datatype type,	double *,double *, double *,
			unsigned int , unsigned int , unsigned int ,
			char *filename);
int write_bin(		datatype type,	double *,double *, double *,
			unsigned int , unsigned int , unsigned int ,
			char *filename);
