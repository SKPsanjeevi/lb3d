/* Read XDR format data from the specified file, write ASCII to stdout.
 *
 * Syntax: xdr2ascii nx ny nz filename type
 *
 * where type is one of scalar, 2scalar, 3scalar, vector
 *
 */

#include <stdio.h>
#include <rpc/rpc.h>

enum datatype {
	TYPE_UNDEFINED	= -1,
	TYPE_SCALAR	= 1,
	TYPE_2SCALAR	= 2,
	TYPE_3SCALAR	= 3,
	TYPE_VECTOR	= 4
};

typedef enum datatype datatype;

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

int main (int argc, char *argv[])
{
	int nx,ny,nz;
	int i;
	/* NB: type "double" should be at least 64 bits wide.
	 * if not, it's trivial to change this code to use long double
	 * for the buffers.
	 */

	double	a,b,c;
	FILE	*infile;
	XDR	xdrs;

	datatype type = TYPE_UNDEFINED;


	/* Ensure correct no of arguments */

	if (6!=argc) {
		fprintf(stderr,
			"Syntax: %s nx ny nz filename type\n"
			"nx, nx, nz = dimensions of data\n"
			"filename = name of input file\n"
			"type = scalar, 2scalar, 3scalar, vector\n",
			argv[0]);
		return -1;
	}

	/* Parse the dimension arguments */

	nx = atoi(argv[1]);
	ny = atoi(argv[2]);
	nz = atoi(argv[3]);

	/* Sanity check them */

	if ( (nx<=0) || (ny<=0) || (nz<=0) ) {
		fprintf(stderr,"Invalid size %d %d %d\n",nx,ny,nz);
		return -1;
	}

	/* Parse the type argument */

	if (TYPE_UNDEFINED == (type = parse_type(argv[5]))) {
		fprintf(stderr,"Bad type %s\n",argv[5]);
		return -1;
	}

	/* Open the output file */

	if (NULL==(infile=fopen(argv[4],"r"))) {
		perror("Could not open input file");
		return -1;
	}

	/* Set up an XDR stream */

	xdrstdio_create(&xdrs,infile,XDR_DECODE);

	/* Now read and parse stdin. */

	switch (type) {
		case TYPE_SCALAR:

			for (i=0;i<nx*ny*nz;i++) {

				/* An XDR double is 64 bits wide. */

				xdr_double(&xdrs,&a);
				printf("%24.20lf\n",a);
			}

			break;

		case TYPE_2SCALAR:
			for (i=0;i<nx*ny*nz;i++) {
				xdr_double(&xdrs,&a);
				xdr_double(&xdrs,&b);
				printf("%24.20lf %24.20lf\n",a,b);
			}
			break;

		case TYPE_3SCALAR:
			for (i=0;i<nx*ny*nz;i++) {
				xdr_double(&xdrs,&a);
				xdr_double(&xdrs,&b);
				xdr_double(&xdrs,&c);
				printf("%lf %lf %lf\n",a,b,c);
			}
			break;

		case TYPE_VECTOR:
			for (i=0;i<nx*ny*nz;i++) {
				xdr_double(&xdrs,&a);
				xdr_double(&xdrs,&b);
				xdr_double(&xdrs,&c);
				printf("%lf %lf %lf\n",a,b,c);
			}
			break;


	}

	/* Close the XDR stream */

	xdr_destroy(&xdrs);

	return 0; /* And they all lived happily ever after. */

}
