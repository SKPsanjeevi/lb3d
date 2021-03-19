/*  3-D VERSION OF THE LIKE
 *  PROGRAM BY J. CHIN
 */


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
	int nx,ny,nz;
#ifdef SGL
	float	buf;
	float phimax,phimin,phisum,phisum2;
	float sd;
#else
	double	buf;
	double phimax,phimin,phisum,phisum2;
	double sd;
#endif


	if (argc!=5) {
		fprintf(stderr,"Syntax: %s <filename> <nx> <ny> <nz>\n",argv[0]);
		return -1;
	}

	fname=argv[1];
	nx=atoi(argv[2]);
	ny=atoi(argv[3]);
	nz=atoi(argv[4]);

	if (NULL==(f=fopen(fname,"r"))) {
		perror("Unable to open input file");
		return -1;
	}

	xdrstdio_create(&xdrs,f,XDR_DECODE);
	
#ifdef SGL
	xdr_float(&xdrs,&buf);
#else
	xdr_double(&xdrs,&buf);
#endif
	phisum=phimax=phimin=buf;
	phisum2=buf*buf;
	for (i=0;i<(nx*ny*nz-1);i++) {
#ifdef SGL
			xdr_float(&xdrs,&buf);
#else
			xdr_double(&xdrs,&buf);
#endif
			phisum+=buf;
			phisum2+=buf*buf;
			if (buf>phimax) { phimax=buf;}
			if (buf<phimin) { phimin=buf;}
	}
	xdr_destroy(&xdrs);

#ifdef SGL
	phisum/=(float)(nx*ny*nz);
	phisum2/=(float)(nx*ny*nz);
#else
	phisum/=(double)(nx*ny*nz);
	phisum2/=(double)(nx*ny*nz);
#endif
	sd=sqrt(phisum2-phisum*phisum);

	printf("Max %f Min %f Mean %f SD %f\n",
		phimax,phimin,phisum,sd);
	fclose(f);
	return 0;
}
