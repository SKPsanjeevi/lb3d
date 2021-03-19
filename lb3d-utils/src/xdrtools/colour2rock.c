/* Read a colour file in XDR format.
 * For each site, write a corresponding obstacle value in ALL format
 * to stdout.
 * If the site has -ve colour, write 0 (no obstacle), else 1.
 *
 * Syntax: colour2rock nx ny nz ox oy oz filename 
 *
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
#define ARGV_FILE 7


int main (int argc, char *argv[])
{
	int nx,ny,nz;
	int ox,oy,oz;
	int x,y,z;
	/* NB: type "double" should be at least 64 bits wide.
	 * if not, it's trivial to change this code to use long double
	 * for the buffers.
	 */

	double	a;
	FILE	*infile;
	XDR	xdrs;



	/* Ensure correct no of arguments */

	if (8!=argc) {
		fprintf(stderr,
			"Syntax: %s nx ny nz ox oy oz filename\n"
			"nx, ny, nz = dimensions of ASCII data\n"
			"ox, oy, oz = offset\n"
			"filename = name of output file\n",
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

	/* Open the output file */

	if (NULL==(infile=fopen(argv[ARGV_FILE],"r"))) {
		perror("Could not open output file");
		return -1;
	}

	/* Set up an XDR stream */

	xdrstdio_create(&xdrs,infile,XDR_DECODE);

	/* Now read and parse stdin. */

	for (z=1;z<=nz;z++) {
	 for (y=1;y<=ny;y++) {
	  for (x=1;x<=nx;x++) {
		xdr_double(&xdrs,&a);
		if (a<=0)
			printf("%d %d %d %d\n",
				x+ox,y+oy,z+oz,1);
		else
			printf("%d %d %d %d\n",
				x+ox,y+oy,z+oz,0);
	  }
	 }
	}

	/* Close the XDR stream */

	xdr_destroy(&xdrs);

	return 0; /* And they all lived happily ever after. */

}
