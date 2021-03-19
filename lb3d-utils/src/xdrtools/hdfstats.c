/*  3-D VERSION OF THE LIKE
 *  PROGRAM BY J. CHIN
 *
 *  version to read scalar hdf5 output files; E. Breitmoser, Oct. 2003
 */

extern void hdfread(char *);

#include <stdio.h>
#include <stdlib.h>
#include <rpc/rpc.h>
#include <math.h>

#define BUFSIZE 1024 /* So shoot me. */


int main (int argc, char *argv[])
{
	FILE	*f;
	int	i;
	XDR	xdrs;
	char *fname;
#ifdef SGL
	float	buf;
	float sd;
#else
	double	buf;
	double sd;
#endif


	if (argc!=2) {
		fprintf(stderr,"Syntax: %s <filename>\n",argv[0]);
		return -1;
	}

	fname=argv[1];

        if(NULL==(f=fopen(fname,"r"))){
	  perror("Unable to open input file");
	  return -1;
	}

	/* Read the hdf5 file and calculate max, min, mean, sd */
	     hdfread(fname); 

	return 0;
}
