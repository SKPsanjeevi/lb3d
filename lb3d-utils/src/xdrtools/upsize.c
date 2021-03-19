#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* VTK upsizer
 *
 * Syntax: ./upsize <vtkfile> <n>
 *
 * Reads in a VTK-format file, and very crudely parses its header.
 * It then scales the dataset up by a factor of n, to generate
 * a larger dataset, which is written to the file "out.vtk".
 * The dataset is also coerced into a 1-byte-per-voxel format, where
 * 0 = minimum value of dataset, 255 = max value, and everything else
 * is linearly interpolated between max and min.
 *
 * XXX WARNING XXX : it assumes input files are 4-byte floating point
 * arrays.  If you write a VTK file from an SGI and try to read it in a
 * linux machine, this will fail due to different floating-point
 * formats. The code to upscale the dataset without changing to
 * byte-per-voxel still works, since it just shuffles data around
 * without trying to parse it.
 *
 * eg, if foo.xdr is 128^3, then `upsize foo.xdr 3` will write a 384^3
 * file to out.xdr.
 * 
 * Initial version October 2003 Jonathan Chin
 */

/* TODO:
 *
 * Allow byte-ification to be turned on or off from command line
 * Proper header parsing.
 */

/* The following symbols are defined to allow large (>2GB) file
 * support under linux.
 */
#ifndef _FILE_OFFSET_BITS
#define _FILE_OFFSET_BITS 64
#endif
#ifndef _LARGEFILE_SOURCE
#define _LARGEFILE_SOURCE
#endif

const char *cvsid="$Id$";

struct dataset {
	unsigned int nx,ny,nz;
	void *data;
};

/* Read a line from the given file descriptor.
 * Return 1 if the line matches the given string, else
 * print an error message and return 0.
 */
int matchline(FILE *f, char *str)
{
	char linebuf[1024]; /* Hatchet job */

	if (1!=fscanf(f,"%s",linebuf)) { 
		fprintf(stderr,"Premature end of file!\n");
		return 0;
	}

	if (0!=strcmp(str,linebuf)) { 
		fprintf(stderr,"Expected \"%s\", got \"%s\"\n",str,linebuf);
		return 0;
	}
	return 1;
}

/* Take the path of a VTK file containing a floating-point
 * structured points dataset. Read it in, and return
 * the corresponding dataset structure, or NULL on error.
 */

struct dataset *dataset_fromvtkfloats(char *infilename) 
{
	struct dataset *dset=NULL;
	int nx=128,ny=128,nz=128;
	float *data=NULL;
	FILE *f=NULL;
	int junk;

	if (NULL==(f=fopen(infilename,"r"))) {
		perror("fopen()"); goto error;
	}

	while (0x0a != fgetc(f)); /* Skip first line. */
	while (0x0a != fgetc(f)); /* Skip second line. */
	if (!matchline(f,"BINARY")) { goto misparse; }
	if (!matchline(f,"DATASET")) { goto misparse; }
	if (!matchline(f,"STRUCTURED_POINTS")) { goto misparse; }
	if (!matchline(f,"DIMENSIONS")) { goto misparse; }
	if (3!=fscanf(f,"%d %d %d",&nx,&ny,&nz)) { goto misparse; }
	if (!matchline(f,"SPACING")) { goto misparse; }
	if (3!=fscanf(f,"%d %d %d",&junk,&junk,&junk)) { goto misparse; }
	if (!matchline(f,"ORIGIN")) { goto misparse; }
	if (3!=fscanf(f,"%d %d %d",&junk,&junk,&junk)) { goto misparse; }
	if (!matchline(f,"POINT_DATA")) { goto misparse; }
	if (1!=fscanf(f,"%d",&junk)) { goto misparse; }
	if (junk != nx*ny*nz) {
		fprintf(stderr,"Expected %d, got %d\n",nx*ny*nz,junk);
		goto misparse;
	}
	if (!matchline(f,"SCALARS")) { goto misparse; }
	if (!matchline(f,"scalars")) { goto misparse; }
	if (!matchline(f,"float")) { goto misparse; }
	if (!matchline(f,"LOOKUP_TABLE")) { goto misparse; }
	if (!matchline(f,"default")) { goto misparse; }
	while (0x0a != fgetc(f)); /* Skip past end of last line. */

	/* input filehandle is at beginning of binary data */

	if (4!=sizeof(float)) {
		fprintf(stderr,"Your floats are the wrong size!\n"); goto error;
	}

	if (NULL==(data = malloc(4*nx*ny*nz))) {
		perror("malloc()"); goto error;
	}
	
	if (nx*ny*nz != fread(data,4,nx*ny*nz,f)) {
		perror("fread()"); goto error;
	}

	if (NULL==(dset=malloc(sizeof(struct dataset)))) {
		perror("malloc()"); goto error;
	}
	dset->nx=nx; dset->ny=ny; dset->nz=nz; dset->data = data;
	fclose(f);
	return dset;

misparse:
	fprintf(stderr,"Unable to parse input file headers.\n");
error:
	if (data) { free(data); }
	if (dset) { free(dset); }
	return NULL;
}

/* Take a dataset, assumed to be in float format.
 * Upscale and write to the named file.
 * Return 1 on success, else zero.
 */

int upscale_to_vtkfloat(struct dataset *dset, char *outfilename, int nreps)
{
	FILE *outfile=NULL;
	int x,y,xrep,yrep,zrep;
	int nx,ny,nz;
	float *data;

	nx = dset->nx;
	ny = dset->ny;
	nz = dset->nz;
	data = dset->data;

	if (NULL==(outfile=fopen(outfilename,"w"))) {
		perror("fopen()");
		return 0;
	}

	fprintf(outfile,"# vtk DataFile Version 2.0\n");
	fprintf(outfile,"Generated by grind.c\n");
	fprintf(outfile,"BINARY\n");
	fprintf(outfile,"DATASET STRUCTURED_POINTS\n");
	fprintf(outfile,"DIMENSIONS %d %d %d\n",nx*nreps,ny*nreps,nz*nreps);
	fprintf(outfile,"SPACING 1 1 1\n");
	fprintf(outfile,"ORIGIN 0 0 0\n");
	fprintf(outfile,"POINT_DATA %d\n",nx*ny*nz*nreps*nreps*nreps);
	fprintf(outfile,"SCALARS scalars float\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");

	for (xrep=0;xrep<nreps;xrep++) {
	for (x=0;x<nx;x++) {
		for (yrep=0;yrep<nreps;yrep++) {
		for (y=0;y<ny;y++) {
			for (zrep=0;zrep<nreps;zrep++) {
				fwrite(
					data + y*nz + x * ny * nz,
					sizeof(float),
					nz,
					outfile
				);
			} /* Z */
		}} /* Y */
	}} /* X */

	fclose(outfile);
	return 1;
}

/* Take a float dataset, return a new byte dataset, or NULL on error */
struct dataset *dataset_to_bytes(struct dataset *dset)
{
	unsigned int nx,ny,nz;
	unsigned char *outdata=NULL;
	struct dataset *outdset=NULL;
	float phi,phimax,phimin;
	float *p;
	unsigned char *q;
	int i;

	nx=dset->nx; ny=dset->ny; nz=dset->nz;

	if (NULL==(outdata=malloc(nx*ny*nz))) {
		perror("dataset_to_bytes: malloc()"); goto error;
	}

	if (NULL==(outdset=malloc(sizeof(struct dataset)))) {
		perror("dataset_to_bytes: malloc()"); goto error;
	}

	/* Get max and min values */

	p=dset->data;
	phimax = phimin = *p;
	for (i=0;i<nx*ny*nz;i++) {
		phi = *p++;
		phimax = (phi>phimax) ? phi : phimax;
		phimin = (phi<phimin) ? phi : phimin;
	}

	q=outdata; p=dset->data;

	for (i=0;i<nx*ny*nz;i++) {
		phi = *p++;
		phi = (phi-phimin)/(phimax-phimin)*255;
		if (phi>255) { phi=255; } if (phi<0) { phi=0; }
		*q++ = (unsigned char)phi;
	}
	outdset->nx = nx;
	outdset->ny = ny;
	outdset->nz = nz;
	outdset->data = outdata;
	return outdset;

error:
	if (outdata) { free(outdata); }
	if (outdset) { free(outdset); }
	return NULL;
}

/* Write a byte dataset to the named file. 1 on success, 0 on error. */

int upscale_to_vtkbyte(struct dataset *dset, char *outfilename, int nreps)
{
	FILE *outfile=NULL;
	int x,y,xrep,yrep,zrep;
	int nx,ny,nz;
	unsigned char *data;

	nx = dset->nx;
	ny = dset->ny;
	nz = dset->nz;
	data = dset->data;

	if (NULL==(outfile=fopen(outfilename,"w"))) {
		perror("fopen()");
		return 0;
	}

	fprintf(outfile,"# vtk DataFile Version 2.0\n");
	fprintf(outfile,"Generated by grind.c\n");
	fprintf(outfile,"BINARY\n");
	fprintf(outfile,"DATASET STRUCTURED_POINTS\n");
	fprintf(outfile,"DIMENSIONS %d %d %d\n",nx*nreps,ny*nreps,nz*nreps);
	fprintf(outfile,"SPACING 1 1 1\n");
	fprintf(outfile,"ORIGIN 0 0 0\n");
	fprintf(outfile,"POINT_DATA %d\n",nx*ny*nz*nreps*nreps*nreps);
	fprintf(outfile,"SCALARS scalars unsigned_char\n");
	fprintf(outfile,"LOOKUP_TABLE default\n");

	for (xrep=0;xrep<nreps;xrep++) {
	for (x=0;x<nx;x++) {
		for (yrep=0;yrep<nreps;yrep++) {
		for (y=0;y<ny;y++) {
			for (zrep=0;zrep<nreps;zrep++) {
				fwrite(
					data + y*nz + x * ny * nz,
					sizeof(unsigned char),
					nz,
					outfile
				);
			} /* Z */
		}} /* Y */
	}} /* X */

	fclose(outfile);
	return 1;
}

int main(int argc, char *argv[])
{

	char *infilename = "foo.vtk";
	char *outfilename="out.vtk";

	int nreps=2; /* x6 = 768 */
	struct dataset *dset=NULL;
	struct dataset *byteset=NULL;

	if (1!=argc) { infilename = argv[1];  }
	if (2<=argc) { nreps = atoi(argv[2]); }

	if (NULL==(dset=dataset_fromvtkfloats(infilename))) {
		fprintf(stderr,"Unable to read input data\n");
		return -1;
	}

	fprintf(stderr,"Read %dx%dx%d voxels.\n", dset->nx,dset->ny,dset->nz);
	fflush(stderr);
	fprintf(stderr,"Upscaling by factor of %d to %dx%dx%d...\n",
			nreps,dset->nx*nreps,dset->ny*nreps,dset->nz*nreps);

	if (NULL==(byteset=dataset_to_bytes(dset))) {
		fprintf(stderr,"Failed byte conversion\n");
		return -1;
	}
	fprintf(stderr,"ok\n");

	/* now fake up our own VTK file.. */
	if (!upscale_to_vtkbyte(byteset,outfilename,nreps)) {
		fprintf(stderr,"Failed to write dataset\n");
		return -1;
	}

	fprintf(stderr,"All done!\n");
	
	return 0;
}