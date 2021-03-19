/* Downsize XDR data */
#include <stdio.h>
#include <stdlib.h>
#include <rpc/rpc.h>


	
int main (int argc, char *argv[])
{

	FILE *infile,*outfile;
	char *infname,*outfname;
	XDR xdrs;
	int nx,ny,nz,ox,oy,oz;
	int newx,newy,newz;
	int i,x,y,z;

#ifdef SGL	
	float *phi,*p;
	float sum;
#else
	double *phi,*p;
	double sum;
#endif	

	if (argc!=12) {
		fprintf(stderr,
	"Usage: %s <nx> <ny> <nz> <newx> <newy> <newz> <ox> <oy> <oz> <infile> <outfile>\n",
			argv[0]);
		return -1;
	}

	nx=atoi(argv[1]);
	ny=atoi(argv[2]);
	nz=atoi(argv[3]);

	newx=atoi(argv[4]);
	newy=atoi(argv[5]);
	newz=atoi(argv[6]);

	ox=atoi(argv[7]);
	oy=atoi(argv[8]);
	oz=atoi(argv[9]);

	infname=argv[10];
	outfname=argv[11];
	if ( (nx<=0) || (ny<=0) ) {
		fprintf(stderr,"Bad dimensions: %d x %d\n",nx,ny);
		return -1;
	}

#ifdef SGL	
	if (NULL==(phi=malloc(sizeof(float)*nx*ny*nz))) {
		perror("malloc()");
		return -1;
	}
#else
	if (NULL==(phi=malloc(sizeof(double)*nx*ny*nz))) {
		perror("malloc()");
		return -1;
	}
#endif

	if (NULL==(infile=fopen(infname,"r"))) {
		perror("fopen()");
		return -1;
	}
	if (NULL==(outfile=fopen(outfname,"w"))) {
		perror("fopen()");
		return -1;
	}

	/* Snarf the data into one large buffer, taking a sum as we go along.
	 */

	xdrstdio_create(&xdrs,infile,XDR_DECODE);

	p=phi;
	sum=0.0;
	for (i=0;i<nx*ny*nz;i++) {
#ifdef SGL	
		if (1!=xdr_float(&xdrs,p++)) {
			fprintf(stderr,"xdr_float() failed on datum %d.\n",i);
			return -1;
		}
#else
		if (1!=xdr_double(&xdrs,p++)) {
			fprintf(stderr,"xdr_double() failed on datum %d.\n",i);
			return -1;
		}
#endif
		sum+= *(p-1);
	}
	xdr_destroy(&xdrs);
	fclose(infile);

	xdrstdio_create(&xdrs,outfile,XDR_ENCODE);

	p=phi;
	for (z=0;z<newz;z++) {
	for (y=0;y<newy;y++) {
	for (x=0;x<newx;x++) {
#ifdef SGL
		float *val;
#else
		double *val;
#endif
		val=(phi + nx*ny*(z+oz) + nx*(y+oy) + (x+ox));
#ifdef SGL	
		xdr_float(&xdrs,val);
#else
		xdr_double(&xdrs,val);
#endif
	}
	}
	}

	fclose(outfile);
	free(phi);
	return 0;
}
