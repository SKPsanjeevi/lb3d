#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>

/* Downsize XMT raw data.
 * This data are integers between 0 and 255.
 * A threshold is set to discriminate between pore sites (<=threshold) and rock sites (>threshold),
 * and only the rock sites and relative coordinates are written on file */


int main(int argc, char *argv[])
{

	int nx,ny,nz;
	int xmin,ymin,zmin,xmax,ymax,zmax;
	int x,y,z;
	int threshold=115;
	unsigned int ic;
	int fc;
	char *infilename,*outfilename;
	FILE *infile,*outfile;
	unsigned char *buffer=NULL;
	unsigned char *p=NULL;

	if (12!=argc) {
		fprintf(stderr,
	"Syntax: %s <nx> <ny> <nz> <xmin> <ymin> <zmin> <xmax> <ymax> <zmax> <file.raw> <outfile>\n",

		argv[0]);
		return -1;
	}

	nx=atoi(argv[1]);
	ny=atoi(argv[2]);
	nz=atoi(argv[3]);
	xmin=atoi(argv[4]);
	ymin=atoi(argv[5]);
	zmin=atoi(argv[6]);
	xmax=atoi(argv[7]);
	ymax=atoi(argv[8]);
	zmax=atoi(argv[9]);
	infilename=argv[10];
	outfilename=argv[11];


	if (NULL==(buffer=malloc(nx*ny*nz*sizeof(char)))) {
		perror("malloc()");
		return -1;
	}
	if (NULL==(infile=fopen(infilename,"r"))) {
		perror("fopen()");
		return -1;
	}

	if (0==fread(buffer,sizeof(char),nx*ny*nz,infile)) {
		perror("fread()");
		return -1;
	}
	fclose(infile);
		
	if (NULL==(outfile=fopen(outfilename,"w"))) {
		perror("fopen()");
		return -1;
	}

	p=buffer;
	for (z=1;z<=nz;z++) {
	for (y=1;y<=ny;y++) {
	for (x=1;x<=nx;x++) {
		ic = *p++;
		if (
			(x<=xmax) &&
			(y<=ymax) &&
			(z<=zmax) &&
			(x>=xmin) &&
			(y>=ymin) &&
			(z>=zmin) ) { 
		            if(ic>threshold)
		                fc=1;
		             else
		                fc=0;
                      /*  fprintf(outfile,"%d %d %d %d\n",x-xmin+1,y-ymin+1,z-zmin+1,fc); */
		     if(fc==1) fprintf(outfile,"%d %d %d %d\n",x-xmin+1,y-ymin+1,z-zmin+1,fc);
		}
			
	}}}
	fclose(outfile);
	
	nx = xmax-xmin+1;
	ny = ymax-ymin+1;
	nz = zmax-zmin+1;
	printf("nx= %d ny= %d nz= %d \n",nx,ny,nz);


	return 0;
}
