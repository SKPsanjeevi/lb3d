/* Read XDR from the specified file, write ASCII to stdout.
 *
 * Syntax: xdr2all nx ny nz ox oy oz type filename 
 *
 * where type is one of scalar, 2scalar, 3scalar, vector
 *
 */

#include <stdio.h>
#include <rpc/rpc.h>

#define ARGV_NX	1
#define ARGV_NY	2
#define ARGV_NZ	3
#define ARGV_OX 4
#define ARGV_OY 5
#define ARGV_OZ 6
#define ARGV_TYPE 7
#define ARGV_FILE 8

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
	int ox,oy,oz;
	int x,y,z;

	/* int i; */
	/* NB: type "double" should be at least 64 bits wide.
	 * if not, it's trivial to change this code to use long double
	 * for the buffers.
	 */

#ifdef SGL
	float	a,b,c;
#else
	double	a,b,c;
#endif

	FILE	*infile;
	XDR	xdrs;

	datatype type = TYPE_UNDEFINED;


	/* Ensure correct no of arguments */

	if (9!=argc) {
		fprintf(stderr,
			"Syntax: %s nx ny nz ox oy oz type filename\n"
			"nx, ny, nz = dimensions of ASCII data\n"
			"ox, oy, oz = offset\n"
			"filename = name of output file\n"
			"type = scalar, 2scalar, 3scalar, vector\n",
			argv[0]);
		return -1;
	}

	/* Parse the dimension arguments */

	nx = atoi(argv[ARGV_NX]);
	ny = atoi(argv[ARGV_NY]);
	nz = atoi(argv[ARGV_NZ]);

	/* Sanity check them */

	if ( (nx<=0) || (ny<=0) || (nz<=0) ) {
		fprintf(stderr,"Invalid size %d %d %d\n",nx,ny,nz);
		return -1;
	}

	/* Parse the offset arguments */

	ox = atoi(argv[ARGV_OX]);
	oy = atoi(argv[ARGV_OY]);
	oz = atoi(argv[ARGV_OZ]);

	/* Parse the type argument */

	if (TYPE_UNDEFINED == (type = parse_type(argv[ARGV_TYPE]))) {
		fprintf(stderr,"Bad type %s\n",argv[ARGV_TYPE]);
		return -1;
	}

	/* Open the output file */

	if (NULL==(infile=fopen(argv[ARGV_FILE],"r"))) {
		perror("Could not open output file");
		return -1;
	}

	/* Set up an XDR stream */

	xdrstdio_create(&xdrs,infile,XDR_DECODE);

	/* Now read and parse stdin. */

	switch (type) {
		case TYPE_SCALAR:

			for (z=1;z<=nz;z++) {
			 for (y=1;y<=ny;y++) {
			  for (x=1;x<=nx;x++) {
#ifdef SGL
				xdr_float(&xdrs,&a);
				printf("%d %d %d %e\n",
					x+ox,y+oy,z+oz,a);
#else
				xdr_double(&xdrs,&a);
				printf("%d %d %d %.20lf\n",
					x+ox,y+oy,z+oz,a);
#endif
			  }
			 }
			}

			break;

		case TYPE_2SCALAR:


			for (z=1;z<=nz;z++) {
			 for (y=1;y<=ny;y++) {
			  for (x=1;x<=nx;x++) {
#ifdef SGL
				xdr_float(&xdrs,&a);
				xdr_float(&xdrs,&b);
				printf("%d %d %d %e %e\n",
					x+ox,y+oy,z+oz,a,b);
#else
				xdr_double(&xdrs,&a);
				xdr_double(&xdrs,&b);
				printf("%d %d %d %.20lf %.20lf\n",
					x+ox,y+oy,z+oz,a,b);
#endif
			  }
			 }
			}


			break;

		case TYPE_3SCALAR:
		case TYPE_VECTOR:

			for (z=1;z<=nz;z++) {
			 for (y=1;y<=ny;y++) {
			  for (x=1;x<=nx;x++) {
#ifdef SGL
				xdr_float(&xdrs,&a);
				xdr_float(&xdrs,&b);
				xdr_float(&xdrs,&c);
				printf("%d %d %d %e %e %e\n",
					x+ox,y+oy,z+oz,a,b,c);
#else
				xdr_double(&xdrs,&a);
				xdr_double(&xdrs,&b);
				xdr_double(&xdrs,&c);
				printf("%d %d %d %.20lf %.20lf %.20lf\n",
					x+ox,y+oy,z+oz,a,b,c);
#endif
			  }
			 }
			}

			break;

	}

	/* Close the XDR stream */

	xdr_destroy(&xdrs);

	return 0; /* And they all lived happily ever after. */

}
