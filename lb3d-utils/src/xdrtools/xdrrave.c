/* Read in an XDR output file, and radially average the values */

#include <stdio.h>
#include <rpc/rpc.h>
#include <math.h>

int main (int argc, char *argv[])
{
	int	nx,ny,nz;
	int	*counter;
	int	nbins;
	FILE	*f;
	int	i,j,k;
#ifdef SGL
	float	cx,cy,cz;
	float	*bin;
	float	r;
	float	x,y,buf;
#else
	double	cx,cy,cz;
	double	*bin;
	double	r;
	double	x,y,buf;
#endif
	XDR	xdrs;

	if (argc!=5) {
		fprintf(stderr,"Syntax: %s nx ny nz filename\n",
			argv[0]);
		return -1;
	}

	nx=atoi(argv[1]);
	ny=atoi(argv[2]);
	nz=atoi(argv[3]);

	cx=nx/2.0;
	cy=ny/2.0;
	cz=nz/2.0;


	nbins=(int)ceil(sqrt(nx*nx+ny*ny));

#ifdef SGL
	if (NULL==(bin=(float *)calloc(nbins,sizeof(float)))) {
		perror("Unable to allocate binnery");
		return -1;
	}
#else
	if (NULL==(bin=(double *)calloc(nbins,sizeof(double)))) {
		perror("Unable to allocate binnery");
		return -1;
	}
#endif
	if (NULL==(counter=(int *)calloc(nbins,sizeof(int)))) {
		perror("Unable to allocate counters");
		return -1;
	}

	if (NULL==(f=fopen(argv[4],"r"))) {
		perror("Unable to open input file");
		return -1;
	}

	xdrstdio_create(&xdrs,f,XDR_DECODE);

	/* Initialise binnery and counters */

	for (i=0;i<nbins;i++) {
		bin[i]=0;
		counter[i]=0;
	}

	for (k=1;k<=nz;k++) {
		for (j=1;j<=ny;j++) {
			y=j-0.5;
			for (i=1;i<=nx;i++) {
				x=i-0.5;
				r=sqrt(	(x-cx)*(x-cx)+
					(y-cy)*(y-cy)
					);
				counter[(int)floor(r)]++;
#ifdef SGL
				xdr_float(&xdrs,&buf);
#else
				xdr_double(&xdrs,&buf);
#endif
				bin[(int)floor(r)]+=buf;
			}
		}
	}

	xdr_destroy(&xdrs);
	fclose(f);

	/* Now calculate and output the radial averages */

	for (i=0;i<nbins;i++) {
		if (counter[i]!=0)
#ifdef SGL
		printf("%f %f\n",(float)i+0.5,bin[i]/(float)counter[i]);
#else
		printf("%f %f\n",(double)i+0.5,bin[i]/(double)counter[i]);
#endif
	}

	return 0;
}
