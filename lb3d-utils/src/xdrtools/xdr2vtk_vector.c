#include <stdio.h>
#include <stdlib.h>
#include <rpc/rpc.h>


	
int main (int argc, char *argv[])
{

	FILE *infile,*outfile;
	char *infname,*outfname;
	XDR xdrs;
        float p2;
	int nx,ny,nz;
	int i,x,y,z,t;
#ifdef SGL
	float *phi,*p;
	float sum;
#else
	double *phi,*p;
	double sum;
#endif

	if (argc!=6) {
		fprintf(stderr,"Usage: %s <nx> <ny> <nz> <infile> <outfile>\n",
			argv[0]);
		return -1;
	}

	nx=atoi(argv[1]);
	ny=atoi(argv[2]);
	nz=atoi(argv[3]);
	infname=argv[4];
	outfname=argv[5];
	if ( (nx<=0) || (ny<=0) ) {
		fprintf(stderr,"Bad dimensions: %d x %d\n",nx,ny);
		return -1;
	}

#ifdef SGL
	if (NULL==(phi=malloc(sizeof(float)*nx*ny*nz*3))) {
		perror("malloc()");
		return -1;
	}
#else
	if (NULL==(phi=malloc(sizeof(double)*nx*ny*nz*3))) {
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
	for (i=0;i<nx*ny*nz*3;i++) {
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

	/* Set "sum" to contain the mean */

#ifdef SGL
	sum/=(float)(nx*ny);
#else
	sum/=(double)(nx*ny);
#endif

	/* Start writing vtk file */

	fprintf(outfile,"# vtk DataFile Version 2.0\nGenerated by xdr2vtk\nBINARY\n");
	fprintf(outfile,"DATASET STRUCTURED_POINTS\n");
	fprintf(outfile,"DIMENSIONS %d %d %d\n",nx,ny,nz);
	fprintf(outfile,"SPACING 1 1 1\n");
	fprintf(outfile,"ORIGIN 0 0 0\n");
	fprintf(outfile,"POINT_DATA %d\n",nx*ny*nz);
	fprintf(outfile,"VECTORS vectors float\n");
/*	fprintf(outfile,"LOOKUP_TABLE default\n");*/

	p=phi;
	for (z=0;z<nz;z++) {
	for (y=0;y<ny;y++) {
	for (x=0;x<nx;x++) {
	for (t=0;t<3;t++) {
                p2 = *p++;
		fwrite(&p2,sizeof(float),1,outfile);
	}
	}
	}
	}

	fclose(outfile);
	free(phi);
	return 0;
}
