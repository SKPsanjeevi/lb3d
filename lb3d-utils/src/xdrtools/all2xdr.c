/* Read ASCII from stdin, write XDR format data do the specified file.
 *
 * Syntax: ascii2xdr nx ny nz ox oy oz type filename 
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
	int i;
	/* NB: type "double" should be at least 64 bits wide.
	 * if not, it's trivial to change this code to use long double
	 * for the buffers.
	 */

	double	*a,*b,*c;
	double *ap,*bp,*cp;
	double	ba,bb,bc;
	FILE	*outfile;
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

	if (NULL==(outfile=fopen(argv[ARGV_FILE],"w"))) {
		perror("Could not open output file");
		return -1;
	}

	/* Set up an XDR stream */

	xdrstdio_create(&xdrs,outfile,XDR_ENCODE);

	/* Allocate the first array (which will always be used) */

	if (NULL==(a=(double *)calloc(nx*ny*nz,sizeof(double)))) {
		perror("Unable to allocate array");
		return -1;
	}

	/* Now read and parse stdin. */

	switch (type) {
		case TYPE_SCALAR:

			for (i=0;i<nx*ny*nz;i++) {
				fscanf(stdin,"%d %d %d %lf \n",
					&x,&y,&z,&ba);

				x=x-ox-1;
				y=y-oy-1;
				z=z-oz-1;

				/* Sanity check */

				if (	   (x<0) || (x>=nx)
					|| (y<0) || (y>=ny)
					|| (z<0) || (z>=nz)
				) {
					fprintf(stderr,
						"Bad coords %d %d %d\n",
						x+ox+1,y+oy+1,z+oz+1
						);
					return -1;
				}

				/* Store the array in Fortran order */

				*(a + z*ny*nx + y*nx + x)=ba;
			}

			/* Note that this does not check for
			 * multiply-defined cells, and will write
			 * rubbish for those cells and undefined ones
			 * if the input contains them.
			 */

			ap = a;
			for (z=0;z<nz;z++) {
			 for (y=0;y<ny;y++) {
			  for (x=0;x<nx;x++) {
				xdr_double(&xdrs,ap++);
			  }
			 }
			}

			free(a);
			break;

		case TYPE_2SCALAR:

			if (NULL==(b=(double *)calloc(nx*ny*nz,sizeof(double)))) {
				perror("Unable to allocate array");
				return -1;
			}

			for (i=0;i<nx*ny*nz;i++) {
				fscanf(stdin,"%d %d %d %lf %lf \n",
					&x,&y,&z,&ba,&bb);

				x=x-ox-1;
				y=y-oy-1;
				z=z-oz-1;

				/* Sanity check */

				if (	   (x<0) || (x>=nx)
					|| (y<0) || (y>=ny)
					|| (z<0) || (z>=nz)
				) {
					fprintf(stderr,
						"Bad coords %d %d %d\n",
						x+ox+1,y+oy+1,z+oz+1);
					return -1;
				}

				/* Store the array in Fortran order */

				/* This is a nonoptimal memory access
				 * pattern. It could be improved. */

				*(a + z*ny*nx + y*nx + x)=ba;
				*(b + z*ny*nx + y*nx + x)=bb;
			}

			ap = a;
			bp = b;
			for (z=0;z<nz;z++) {
			 for (y=0;y<ny;y++) {
			  for (x=0;x<nx;x++) {
				xdr_double(&xdrs,ap++);
				xdr_double(&xdrs,bp++);
			  }
			 }
			}

			free(a);
			free(b);

			break;

		case TYPE_3SCALAR:
		case TYPE_VECTOR:

			if (NULL==(b=(double *)calloc(nx*ny*nz,sizeof(double)))) {
				perror("Unable to allocate array");
				return -1;
			}
			if (NULL==(c=(double *)calloc(nx*ny*nz,sizeof(double)))) {
				perror("Unable to allocate array");
				return -1;
			}

			for (i=0;i<nx*ny*nz;i++) {
				fscanf(stdin,"%d %d %d %lf %lf %lf \n",
					&x,&y,&z,&ba,&bb,&bc);
				fflush(stdout);

				x=x-ox-1;
				y=y-oy-1;
				z=z-oz-1;

				/* Sanity check */

				if (	   (x<0) || (x>=nx)
					|| (y<0) || (y>=ny)
					|| (z<0) || (z>=nz)
				) {
					fprintf(stderr,
						"Bad coords %d %d %d\n",
						x+ox+1,y+oy+1,z+oz+1);
					return -1;
				}

				/* Store the array in Fortran order */

				/* This is a nonoptimal memory access
				 * pattern. It could be improved. */

				*(a + z*ny*nx + y*nx + x)=ba;
				*(b + z*ny*nx + y*nx + x)=bb;
				*(c + z*ny*nx + y*nx + x)=bc;
			}

			ap = a;
			bp = b;
			cp = c;
			for (z=0;z<nz;z++) {
			 for (y=0;y<ny;y++) {
			  for (x=0;x<nx;x++) {
				xdr_double(&xdrs,ap++);
				xdr_double(&xdrs,bp++);
				xdr_double(&xdrs,cp++);
			  }
			 }
			}

			free(a);
			free(b);
			free(c);

			break;

	}

	/* Close the XDR stream */

	xdr_destroy(&xdrs);

	return 0; /* And they all lived happily ever after. */

}
