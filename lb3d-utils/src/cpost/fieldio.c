
/* Set of routines to read and write
 * XDR and fortran unformatted data.
 */


#include <stdio.h>
#include <stdlib.h>
#include <rpc/rpc.h>
#include <sys/stat.h>

#include "readbin.h"
#include "fieldio.h"

datatype parse_type(char *s) {
	if (0==strcmp("scalar",s)) {
		return TYPE_SCALAR;
	}

	if (0==strcmp("2scalar",s)) {
		return TYPE_2SCALAR;
	}

	if (0==strcmp("3scalar",s)) {
		return TYPE_3SCALAR;
	}

	if (0==strcmp("vector",s)) {
		return TYPE_VECTOR;
	}

	return TYPE_UNDEFINED;

}

/*
 * Read XDR data of the given type into the supplied arrays.
 * Scalar data only requires arr1, 2scalar arr1 and arr2, 3scalar all three.
 * Vector data is read into arr1.
 * Returns 0 if all went well.
 */

int read_xdr(	datatype type,	double *arr1,double *arr2, double *arr3,
			unsigned int nx, unsigned int ny, unsigned int nz,
			char *filename)
{
	FILE *infile;
	XDR	xdrs;
	double *p1,*p2,*p3;
	int i;

	if (NULL==(infile=fopen(filename,"r"))) {
		perror("fopen()");
		return -1;
	}

	xdrstdio_create(&xdrs,infile,XDR_DECODE);

	switch(type) {
		case TYPE_SCALAR:
			p1=arr1;
			for (i=0;i<(nx*ny*nz);i++) {
				xdr_double(&xdrs,p1++);
			}
			break;
	
		case TYPE_2SCALAR:
			p1=arr1;
			p2=arr2;
		
			for (i=0;i<(nx*ny*nz);i++) {
				xdr_double(&xdrs,p1++);
				xdr_double(&xdrs,p2++);
			}
			break;
	
		case TYPE_3SCALAR:
			p1=arr1;
			p2=arr2;
			p3=arr3;
		
			for (i=0;i<(nx*ny*nz);i++) {
				if (1!=xdr_double(&xdrs,p1++)) 
					perror("xdr_double(p1)");
				if (1!=xdr_double(&xdrs,p2++))
					perror("xdr_double(p2)");
				if (1!=xdr_double(&xdrs,p3++))
					perror("xdr_double(p3)");
			}
			break;
	
		case TYPE_VECTOR:
			p1=arr1;
		
			for (i=0;i<(nx*ny*nz);i++) {
				xdr_double(&xdrs,p1++);
				xdr_double(&xdrs,p1++);
				xdr_double(&xdrs,p1++);
			}
			break;
	
		case TYPE_UNDEFINED:
		default:
			return -1;
	}

	fclose(infile);
	return 0;

}

/*
 * Write XDR data of the given type from the supplied arrays.
 * Scalar data only requires arr1, 2scalar arr1 and arr2, 3scalar all three.
 * Vector data is read from arr1.
 * Returns 0 if all went well.
 */

int write_xdr(	datatype type,	double *arr1,double *arr2, double *arr3,
			unsigned int nx, unsigned int ny, unsigned int nz,
			char *filename)
{
	FILE *infile;
	XDR	xdrs;
	double *p1,*p2,*p3;
	int i;

	if (NULL==(infile=fopen(filename,"w"))) {
		perror("fopen()");
		return -1;
	}

	xdrstdio_create(&xdrs,infile,XDR_ENCODE);

	switch(type) {
		case TYPE_SCALAR:
			p1=arr1;
			for (i=0;i<(nx*ny*nz);i++) {
				xdr_double(&xdrs,p1++);
			}
			break;
	
		case TYPE_2SCALAR:
			p1=arr1;
			p2=arr2;
		
			for (i=0;i<(nx*ny*nz);i++) {
				xdr_double(&xdrs,p1++);
				xdr_double(&xdrs,p2++);
			}
			break;
	
		case TYPE_3SCALAR:
			p1=arr1;
			p2=arr2;
			p3=arr3;
		
			for (i=0;i<(nx*ny*nz);i++) {
				if (1!=xdr_double(&xdrs,p1++)) 
					perror("xdr_double(p1)");
				if (1!=xdr_double(&xdrs,p2++))
					perror("xdr_double(p2)");
				if (1!=xdr_double(&xdrs,p3++))
					perror("xdr_double(p3)");
			}
			break;
	
		case TYPE_VECTOR:
			p1=arr1;
		
			for (i=0;i<(nx*ny*nz);i++) {
				xdr_double(&xdrs,p1++);
				xdr_double(&xdrs,p1++);
				xdr_double(&xdrs,p1++);
			}
			break;
	
		case TYPE_UNDEFINED:
		default:
			return -1;
	}

	fclose(infile);
	return 0;

}



/*
 * Read XDR data of the FLOAT type into the supplied arrays.
 * Scalar data only requires arr1, 2scalar arr1 and arr2, 3scalar all three.
 * Vector data is read into arr1.
 * Returns 0 if all went well.
 */

int read_xdr_float(	datatype type,	float *arr1,float *arr2, float *arr3,
			unsigned int nx, unsigned int ny, unsigned int nz,
			char *filename)
{
	FILE *infile;
	XDR	xdrs;
	float *p1,*p2,*p3;
	int i;

	if (NULL==(infile=fopen(filename,"r"))) {
		perror("fopen()");
		return -1;
	}

	xdrstdio_create(&xdrs,infile,XDR_DECODE);

	switch(type) {
		case TYPE_SCALAR:
			p1=arr1;
			for (i=0;i<(nx*ny*nz);i++) {
				xdr_float(&xdrs,p1++);
			}
			break;
	
		case TYPE_2SCALAR:
			p1=arr1;
			p2=arr2;
		
			for (i=0;i<(nx*ny*nz);i++) {
				xdr_float(&xdrs,p1++);
				xdr_float(&xdrs,p2++);
			}
			break;
	
		case TYPE_3SCALAR:
			p1=arr1;
			p2=arr2;
			p3=arr3;
		
			for (i=0;i<(nx*ny*nz);i++) {
				if (1!=xdr_float(&xdrs,p1++)) 
					perror("xdr_float(p1)");
				if (1!=xdr_float(&xdrs,p2++))
					perror("xdr_float(p2)");
				if (1!=xdr_float(&xdrs,p3++))
					perror("xdr_float(p3)");
			}
			break;
	
		case TYPE_VECTOR:
			p1=arr1;
		
			for (i=0;i<(nx*ny*nz);i++) {
				xdr_float(&xdrs,p1++);
				xdr_float(&xdrs,p1++);
				xdr_float(&xdrs,p1++);
			}
			break;
	
		case TYPE_UNDEFINED:
		default:
			return -1;
	}

	fclose(infile);
	return 0;

}

/*
 * Write XDR data of the FLOAT type from the supplied arrays.
 * Scalar data only requires arr1, 2scalar arr1 and arr2, 3scalar all three.
 * Vector data is read from arr1.
 * Returns 0 if all went well.
 */

int write_xdr_float(	datatype type,	float *arr1,float *arr2, float *arr3,
			unsigned int nx, unsigned int ny, unsigned int nz,
			char *filename)
{
	FILE *infile;
	XDR	xdrs;
	float *p1,*p2,*p3;
	int i;

	if (NULL==(infile=fopen(filename,"w"))) {
		perror("fopen()");
		return -1;
	}

	xdrstdio_create(&xdrs,infile,XDR_ENCODE);

	switch(type) {
		case TYPE_SCALAR:
			p1=arr1;
			for (i=0;i<(nx*ny*nz);i++) {
				xdr_float(&xdrs,p1++);
			}
			break;
	
		case TYPE_2SCALAR:
			p1=arr1;
			p2=arr2;
		
			for (i=0;i<(nx*ny*nz);i++) {
				xdr_float(&xdrs,p1++);
				xdr_float(&xdrs,p2++);
			}
			break;
	
		case TYPE_3SCALAR:
			p1=arr1;
			p2=arr2;
			p3=arr3;
		
			for (i=0;i<(nx*ny*nz);i++) {
				if (1!=xdr_float(&xdrs,p1++)) 
					perror("xdr_float(p1)");
				if (1!=xdr_float(&xdrs,p2++))
					perror("xdr_float(p2)");
				if (1!=xdr_float(&xdrs,p3++))
					perror("xdr_float(p3)");
			}
			break;
	
		case TYPE_VECTOR:
			p1=arr1;
		
			for (i=0;i<(nx*ny*nz);i++) {
				xdr_float(&xdrs,p1++);
				xdr_float(&xdrs,p1++);
				xdr_float(&xdrs,p1++);
			}
			break;
	
		case TYPE_UNDEFINED:
		default:
			return -1;
	}

	fclose(infile);
	return 0;

}



/*
 * Read SGI-style Fortran unformatted data of the given type into the supplied arrays.
 * Scalar data only requires arr1, 2scalar arr1 and arr2, 3scalar all three.
 * Vector data is read from arr1.
 * Returns 0 if all went well.
 *
 * NB: this reads into arrays of F_REALs, which may differ from the C double type.
 * Use read_bin() if unsure.
 */

int read_rawbin(	datatype type,	F_REAL *arr1,F_REAL *arr2, F_REAL *arr3,
			unsigned int nx, unsigned int ny, unsigned int nz,
			char *filename)
{
	struct stat statbuf;
	FILE *infile;
	off_t infilesize;
#ifdef	F_USE_HEADER
	F_HEADER headword;
#endif /* F_USE_HEADER */
	char *errmsg;
	int  errmsglen;

	errmsglen = 26 + 1 + strlen(filename);
	errmsg = malloc(sizeof(char)*errmsglen);


	/* Check the length of the input file. */

	if (0!=stat(filename,&statbuf)) {
	   sprintf(errmsg,"Couldn't stat() field file %s", filename);
	   perror(errmsg);
	   return -1;
	}

	infilesize = sizeof(F_REAL) * nx * ny * nz;

	switch(type) {
		case TYPE_SCALAR:
			break;
		case TYPE_2SCALAR:
			infilesize *= 2;
			break;
		case TYPE_3SCALAR:
		case TYPE_VECTOR:
			infilesize *= 3;
			break;
		case TYPE_UNDEFINED:
		default:
			return -1;
	}

#ifdef	F_USE_HEADER
	infilesize += sizeof(F_HEADER)*2;
#endif /* F_USE_HEADER */

	/* If size is wrong, warn but continue. */

	if (infilesize != statbuf.st_size) {
		fprintf(stderr,"WARNING! %s has incorrect file size.\n",filename);
	}

	/* Open the input file */

	if (NULL==(infile=fopen(filename,"r"))) {
		perror("Could not open input file");
		return -1;
	}

#ifdef	F_USE_HEADER
	if (-1==(fread(&headword,sizeof(headword),1,infile))) {
		perror("fread()");
		return -1;
	}
#endif /* F_USE_HEADER */
	switch(type) {
		case TYPE_SCALAR:

			if (-1==fread(arr1,sizeof(F_REAL),nx*ny*nz,infile)) {
				perror("fread()");
				return -1;
			}
			break;

		case TYPE_2SCALAR:
			if (-1==fread(arr1,sizeof(F_REAL),nx*ny*nz,infile)) {
				perror("fread()");
				return -1;
			}
			if (-1==fread(arr2,sizeof(F_REAL),nx*ny*nz,infile)) {
				perror("fread()");
				return -1;
			}
			break;

		case TYPE_3SCALAR:

			if (-1==fread(arr1,sizeof(F_REAL),nx*ny*nz,infile)) {
				perror("fread()");
				return -1;
			}
			if (-1==fread(arr2,sizeof(F_REAL),nx*ny*nz,infile)) {
				perror("fread()");
				return -1;
			}
			if (-1==fread(arr3,sizeof(F_REAL),nx*ny*nz,infile)) {
				perror("fread()");
				return -1;
			}
			break;

		case TYPE_VECTOR:
			if (-1==fread(arr1,sizeof(F_REAL),nx*ny*nz*3,infile)) {
				perror("fread()");
				return -1;
			}
			break;

		default:
			fprintf(stderr,"Warning: bad type passed to read_bin()\n");
			return -1;

		/* Unreachable */
	} /* switch (type) */
	fclose(infile);

	return 0;
}

/*
 * Write SGI-style Fortran unformatted data of the given type from the supplied arrays.
 * Scalar data only requires arr1, 2scalar arr1 and arr2, 3scalar all three.
 * Vector data is read from arr1.
 * Returns 0 if all went well.
 * NB: this writes from arrays of F_REALs, which may differ from the C double type.
 * Use write_bin() if unsure.
 */

int write_rawbin(	datatype type,	F_REAL *arr1,F_REAL *arr2, F_REAL *arr3,
			unsigned int nx, unsigned int ny, unsigned int nz,
			char *filename)
{
	FILE *outfile;
#ifdef	F_USE_HEADER
	F_HEADER headword;
#endif /* F_USE_HEADER */
	off_t	outfilesize;


	/* Calculate the size of the data to be written. */
	/* This goes into the header. */

	outfilesize = sizeof(F_REAL) * nx * ny * nz;

	switch(type) {
		case TYPE_SCALAR:
			break;
		case TYPE_2SCALAR:
			outfilesize *= 2;
			break;
		case TYPE_3SCALAR:
		case TYPE_VECTOR:
			outfilesize *= 3;
			break;
		case TYPE_UNDEFINED:
		default:
			return -1;
	}


	/* Open the output file */

	if (NULL==(outfile=fopen(filename,"w"))) {
		perror("Could not open output file");
		return -1;
	}

#ifdef	F_USE_HEADER
	/* Write the header word */

	headword = outfilesize; /* FIXME - possible overflow in cast? */

	if (-1==(fwrite(&headword,sizeof(headword),1,outfile))) {
		perror("fwrite()");
		return -1;
	}
#endif /* F_USE_HEADER */

	switch(type) {
		case TYPE_SCALAR:

			if (-1==fwrite(arr1,sizeof(F_REAL),nx*ny*nz,outfile)) {
				perror("fwrite()");
				return -1;
			}
			break;

		case TYPE_2SCALAR:
			if (-1==fwrite(arr1,sizeof(F_REAL),nx*ny*nz,outfile)) {
				perror("fwrite()");
				return -1;
			}
			if (-1==fwrite(arr2,sizeof(F_REAL),nx*ny*nz,outfile)) {
				perror("fwrite()");
				return -1;
			}
			break;

		case TYPE_3SCALAR:

			if (-1==fwrite(arr1,sizeof(F_REAL),nx*ny*nz,outfile)) {
				perror("fwrite()");
				return -1;
			}
			if (-1==fwrite(arr2,sizeof(F_REAL),nx*ny*nz,outfile)) {
				perror("fwrite()");
				return -1;
			}
			if (-1==fwrite(arr3,sizeof(F_REAL),nx*ny*nz,outfile)) {
				perror("fwrite()");
				return -1;
			}
			break;

		case TYPE_VECTOR:
			if (-1==fwrite(arr1,sizeof(F_REAL),nx*ny*nz*3,outfile)) {
				perror("fwrite()");
				return -1;
			}
			break;

		default:
			fprintf(stderr,"Warning: bad type passed to read_bin()\n");
			return -1;

		/* Unreachable */
	} /* switch (type) */

#ifdef	F_USE_HEADER
	/* Write the footer word */

	if (-1==(fwrite(&headword,sizeof(headword),1,outfile))) {
		perror("fwrite()");
		return -1;
	}
#endif /* F_USE_HEADER */

	fclose(outfile);

	return 0;
}

/*
 * Read SGI-style Fortran unformatted data of the given type into the supplied arrays.
 * Scalar data only requires arr1, 2scalar arr1 and arr2, 3scalar all three.
 * Vector data is read from arr1.
 * Returns 0 if all went well.
 *
 * NB: all data is cast to the C double() type. This routine calls read_rawbin()
 * directly if sizeof(F_REAL) == sizeof(double).
 */

#ifdef F_REAL_IS_C_DOUBLE
int read_bin(		datatype type,	double *arr1,double *arr2, double *arr3,
			unsigned int nx, unsigned int ny, unsigned int nz,
			char *filename) {
	return read_rawbin(type,arr1,arr2,arr3,nx,ny,nz,filename);
}
#else

int read_bin(		datatype type,	double *arr1,double *arr2, double *arr3,
			unsigned int nx, unsigned int ny, unsigned int nz,
			char *filename)
{
	struct stat statbuf;
	F_REAL *tmparr1=NULL,*tmparr2=NULL,*tmparr3=NULL;
	FILE *infile;
	off_t infilesize;
#ifdef	F_USE_HEADER
	F_HEADER headword;
#endif /* F_USE_HEADER */
	double *cp;	/* Fortran-data pointers */
	F_REAL *fp;	/* C-data pointers */
	int i;


	/* Check the length of the input file. */

	if (0!=stat(filename,&statbuf)) {
		perror("Couldn't stat() input file");
		return -1;
	}

	infilesize = sizeof(F_REAL) * nx * ny * nz;

	switch(type) {
		case TYPE_SCALAR:
			break;
		case TYPE_2SCALAR:
			infilesize *= 2;
			break;
		case TYPE_3SCALAR:
		case TYPE_VECTOR:
			infilesize *= 3;
			break;
		case TYPE_UNDEFINED:
		default:
			return -1;
	}

	infilesize += sizeof(F_HEADER)*2;

	/* If size is wrong, warn but continue. */

	if (infilesize != statbuf.st_size) {
		fprintf(stderr,"WARNING! %s has incorrect file size.\n",filename);
	}

	/* Allocate temporary storage.
	 * (This is the "memory is cheap" method. An alternative
	 * is to fread() one datum at a time, cast it to double, and write
	 * it into the array; this is probably slower I/O wise, but more
	 * memory-efficient.)
	 */

	switch(type) {
		case TYPE_VECTOR:
			if (NULL==(tmparr1=malloc(sizeof(double)*nx*ny*nz*3))) {
				perror("malloc()");
				return -1;
			}
			break;

		case TYPE_3SCALAR:
			if (NULL==(tmparr3=malloc(sizeof(double)*nx*ny*nz))) {
				perror("malloc()");
				return -1;
			}

			/* Fall through */

		case TYPE_2SCALAR:
			if (NULL==(tmparr2=malloc(sizeof(double)*nx*ny*nz))) {
				perror("malloc()");
				return -1;
			}

			/* Fall through */

		case TYPE_SCALAR:
			if (NULL==(tmparr1=malloc(sizeof(double)*nx*ny*nz))) {
				perror("malloc()");
				return -1;
			}
			break;
		default:
			fprintf(stderr,"Can't happen!\n");
			return -1;
	}

	/* Open the input file */

	if (NULL==(infile=fopen(filename,"r"))) {
		perror("Could not open input file");
		return -1;
	}

#ifdef	F_USE_HEADER
	if (-1==(fread(&headword,sizeof(headword),1,infile))) {
		perror("fread()");
		return -1;
	}
#endif /* F_USE_HEADER */

	/* Read the data in */

	switch(type) {
		case TYPE_SCALAR:

			if (-1==fread(tmparr1,sizeof(F_REAL),nx*ny*nz,infile)) {
				perror("fread()");
				return -1;
			}
			break;

		case TYPE_2SCALAR:
			if (-1==fread(tmparr1,sizeof(F_REAL),nx*ny*nz,infile)) {
				perror("fread()");
				return -1;
			}
			if (-1==fread(tmparr2,sizeof(F_REAL),nx*ny*nz,infile)) {
				perror("fread()");
				return -1;
			}
			break;

		case TYPE_3SCALAR:

			if (-1==fread(tmparr1,sizeof(F_REAL),nx*ny*nz,infile)) {
				perror("fread()");
				return -1;
			}
			if (-1==fread(tmparr2,sizeof(F_REAL),nx*ny*nz,infile)) {
				perror("fread()");
				return -1;
			}
			if (-1==fread(tmparr3,sizeof(F_REAL),nx*ny*nz,infile)) {
				perror("fread()");
				return -1;
			}
			break;

		case TYPE_VECTOR:
			if (-1==fread(tmparr1,sizeof(F_REAL),nx*ny*nz*3,infile)) {
				perror("fread()");
				return -1;
			}
			break;

		default:
			fprintf(stderr,"Warning: bad type passed to read_bin()\n");
			return -1;

		/* Unreachable */
	} /* switch (type) */

	/* Convert to C double type, then free the temp storage */

	switch(type) {
		case TYPE_VECTOR:
			fp=tmparr1;
			cp=arr1;
			for (i=0;i<(3*nx*ny*nz);i++) {
				*cp++ = (double) (*fp++);
			}
			free(tmparr1);
			break;
		case TYPE_3SCALAR:
			fp=tmparr3;
			cp=arr3;
			for (i=0;i<(nx*ny*nz);i++) {
				*cp++ = (double) (*fp++);
			}
			free(tmparr3);
			/* Fall through */
		case TYPE_2SCALAR:
			fp=tmparr2;
			cp=arr2;
			for (i=0;i<(nx*ny*nz);i++) {
				*cp++ = (double) (*fp++);
			}
			free(tmparr2);
			/* Fall through */
		case TYPE_SCALAR:
			fp=tmparr1;
			cp=arr1;
			for (i=0;i<(nx*ny*nz);i++) {
				*cp++ = (double) (*fp++);
			}
			free(tmparr1);
			break;
		default:
			fprintf(stderr,"Can't happen.\n");
			return -1;
			
	}
	
	fclose(infile);

	return 0;
}

#endif

/*
 * Write SGI-style Fortran unformatted data of the given type from the supplied arrays.
 * Scalar data only requires arr1, 2scalar arr1 and arr2, 3scalar all three.
 * Vector data is read from arr1.
 * Returns 0 if all went well.
 * NB: this writes from arrays of F_REALs, which may differ from the C double type.
 * Use read_bin() if unsure.
 */

#ifdef F_REAL_IS_C_DOUBLE
int write_bin(		datatype type,	double *arr1,double *arr2, double *arr3,
			unsigned int nx, unsigned int ny, unsigned int nz,
			char *filename) {
	return write_rawbin(type,arr1,arr2,arr3,nx,ny,nz,filename);
}
#else 


int write_bin(		datatype type,	double *arr1,double *arr2, double *arr3,
			unsigned int nx, unsigned int ny, unsigned int nz,
			char *filename)
{
	F_REAL *tmparr1=NULL,*tmparr2=NULL,*tmparr3=NULL;
	FILE *outfile;
	off_t outfilesize;
#ifdef	F_USE_HEADER
	F_HEADER headword;
#endif /* F_USE_HEADER */
	double *cp;	/* Fortran-data pointers */
	F_REAL *fp;	/* C-data pointers */
	int i;


	outfilesize = sizeof(F_REAL) * nx * ny * nz;

	switch(type) {
		case TYPE_SCALAR:
			break;
		case TYPE_2SCALAR:
			outfilesize *= 2;
			break;
		case TYPE_3SCALAR:
		case TYPE_VECTOR:
			outfilesize *= 3;
			break;
		case TYPE_UNDEFINED:
		default:
			return -1;
	}

	/* Allocate temporary storage.
	 */

	switch(type) {
		case TYPE_VECTOR:
			if (NULL==(tmparr1=malloc(sizeof(double)*nx*ny*nz*3))) {
				perror("malloc()");
				return -1;
			}
			break;

		case TYPE_3SCALAR:
			if (NULL==(tmparr3=malloc(sizeof(double)*nx*ny*nz))) {
				perror("malloc()");
				return -1;
			}

			/* Fall through */

		case TYPE_2SCALAR:
			if (NULL==(tmparr2=malloc(sizeof(double)*nx*ny*nz))) {
				perror("malloc()");
				return -1;
			}

			/* Fall through */

		case TYPE_SCALAR:
			if (NULL==(tmparr1=malloc(sizeof(double)*nx*ny*nz))) {
				perror("malloc()");
				return -1;
			}
			break;
		default:
			fprintf(stderr,"Can't happen!\n");
			return -1;
	}

	/* Convert to C double type */

	switch(type) {
		case TYPE_VECTOR:
			fp=tmparr1;
			cp=arr1;
			for (i=0;i<(3*nx*ny*nz);i++) {
				*fp++ = (double) (*cp++);
			}
			break;
		case TYPE_3SCALAR:
			fp=tmparr3;
			cp=arr3;
			for (i=0;i<(nx*ny*nz);i++) {
				*fp++ = (double) (*cp++);
			}
			/* Fall through */
		case TYPE_2SCALAR:
			fp=tmparr2;
			cp=arr2;
			for (i=0;i<(nx*ny*nz);i++) {
				*fp++ = (double) (*cp++);
			}
			/* Fall through */
		case TYPE_SCALAR:
			fp=tmparr1;
			cp=arr1;
			for (i=0;i<(nx*ny*nz);i++) {
				*fp++ = (double) (*cp++);
			}
			break;
		default:
			fprintf(stderr,"Can't happen.\n");
			return -1;
			
	}

	/* Open the output file */

	if (NULL==(outfile=fopen(filename,"w"))) {
		perror("Could not open output file");
		return -1;
	}

#ifdef	F_USE_HEADER
	/* Write the header word */

	headword = outfilesize; /* FIXME - possible overflow in cast? */

	if (-1==(fwrite(&headword,sizeof(headword),1,outfile))) {
		perror("fwrite()");
		return -1;
	}
#endif /* F_USE_HEADER */


	/* Write the data out */

	switch(type) {
		case TYPE_SCALAR:

			if (-1==fwrite(tmparr1,sizeof(F_REAL),nx*ny*nz,outfile)) {
				perror("fwrite()");
				return -1;
			}
			free(tmparr1);
			break;

		case TYPE_2SCALAR:
			if (-1==fwrite(tmparr1,sizeof(F_REAL),nx*ny*nz,outfile)) {
				perror("fwrite()");
				return -1;
			}
			if (-1==fwrite(tmparr2,sizeof(F_REAL),nx*ny*nz,outfile)) {
				perror("fwrite()");
				return -1;
			}
			free(tmparr1);
			free(tmparr2);
			break;

		case TYPE_3SCALAR:

			if (-1==fwrite(tmparr1,sizeof(F_REAL),nx*ny*nz,outfile)) {
				perror("fwrite()");
				return -1;
			}
			if (-1==fwrite(tmparr2,sizeof(F_REAL),nx*ny*nz,outfile)) {
				perror("fwrite()");
				return -1;
			}
			if (-1==fwrite(tmparr3,sizeof(F_REAL),nx*ny*nz,outfile)) {
				perror("fwrite()");
				return -1;
			}
			free(tmparr1);
			free(tmparr2);
			free(tmparr3);
			break;

		case TYPE_VECTOR:
			if (-1==fwrite(tmparr1,sizeof(F_REAL),nx*ny*nz*3,outfile)) {
				perror("fwrite()");
				return -1;
			}
			free(tmparr1);
			break;

		default:
			fprintf(stderr,"Can't happen!\n");
			return -1;

		/* Unreachable */
	} /* switch (type) */

#ifdef	F_USE_HEADER
	/* Write the header word */

	if (-1==(fwrite(&headword,sizeof(headword),1,outfile))) {
		perror("fwrite()");
		return -1;
	}
#endif /* F_USE_HEADER */
	
	fclose(outfile);

	return 0;
}
#endif

