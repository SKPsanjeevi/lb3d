
/* A C-based postprocessor for LBE. 		*/
/* Requires the C<-> Fortran I/O library.	*/
/*						*/

/*/////////////////////////////////////////// */
/*///		NOTES			   // */
/*/// In order to postprocess a job that   // */
/*/// been restarted from a past configur- // */
/*/// ation, just mofify ../code/input-file / */
/*/// such that 'sci_start' equals the     // */
/*/// restarted time step.		   // */
/*/////////////////////////////////////////// */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "namelist.h"
#include "fieldio.h"

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#define TIMECHARS 8        /* CHANGED TO 6 */
#define TIMECHARSTR "8"    /* CHANGED TO 6 */
#define CPUCHARS 6
#define CPUCHARSTR "6"
#define SQUARE(A)   ((A)*(A))
#define PI          (3.14159265358979323846264338327950288)
#define WIN_IN      0.1
#define WIN_OUT     2.0   /*  Condition: 
			   *  (WIN_OUT * if2) <= L
			   *  where L is a half-diagonal length:
			   *  L = nx * sqrt(3)/2
			   *  for a cubic lattice.
			   *
			   */

/* Information about the CPU topology */

int ncpus=0;
struct cpupos {
	unsigned int x;
	unsigned int y;
	unsigned int z;
};
struct cpupos *topology;
unsigned int sx,sy,sz; /* Size of subdomain */

/* Global variables related to the simulation */

char   *pathstub="./";
int    nx,ny,nz;
int    n_iteration,n_sci,sci_start;
int    sci_int,sci_od,sci_wd,sci_vel,sci_flo,sci_arrows,sci_dir,sci_sur;
int    sci_spress = 1, init_cond;
double fr1,fr2;   /* Sizes for vesicle */

/* Allocate a buffer, read a file into it, null-terminate the buffer,
 * return either the buffer or NULL on error.
 */

char *readfile(char *filename)
{
	struct stat sbuf;
	off_t filesize;
	char *buf;
	FILE *infile;

	if (0!=stat(filename,&sbuf)) {
		perror("stat()");
		return NULL;
	}

	filesize = sbuf.st_size ;

	if (NULL==(buf=malloc(filesize+1))) {
		perror("malloc()");
		return NULL;
	}

	if (NULL==(infile=fopen(filename,"r"))) {
		perror("fopen()");
		return NULL;
	}

	if (-1==fread(buf,1,filesize,infile)) {
		perror("fread()");
		return NULL;
	}
	fclose(infile);
	buf[filesize]=0;
	return buf;
}

int postprocess_taskfarm() {

	fprintf(stderr,"Taskfarming unimplemented\n");
	return -1;
}

/* Determine CPU topology from given filename.
 * If this is not possible, return zero, otherwise
 * return the number of CPUs.
 */


int find_topology(char *filename)
{
	char *buffer,*p;
	int i;
	int j,x,y,z;
	unsigned int maxx=0,maxy=0,maxz=0;

	ncpus=0;
	if (NULL==(buffer=readfile(filename))) {
		return -1;
	}
	/* Number of CRs in file == number of CPUs. */

	for (p=buffer;0!=*p;p++) {
		if ('\n'==*p) {
			ncpus++;
		}
	}

	if (NULL==(topology=malloc(sizeof(struct cpupos)*ncpus))) {
		perror("malloc()");
		return -1;
	}

	/* Parse in the CPU positions */

	for (i=0;i<ncpus;i++) {
		topology[i].x=-1;
	}

	p=buffer;
	for (i=0;i<ncpus;i++) {
		sscanf(p,"%d %d %d %d\n",&j,&x,&y,&z);

		if ( (x<0)||(y<0)||(z<0) ) {
			printf("Bad CPU coords %d %d %d\n",x,y,z);
			return -1;
		}

		if ( (j>=0) && (j<ncpus) ) {
			topology[j].x=x;
			topology[j].y=y;
			topology[j].z=z;
			while('\n'!=*p) p++;
			while('\n'==*p) p++;
			printf("CPU %d at (%d, %d, %d)\n",j,
				topology[j].x,
				topology[j].y,
				topology[j].z);
		} else {
			fprintf(stderr,"Bad CPU number %u\n",j);
			return -1;
		}
	}

	for (i=0;i<ncpus;i++) {
		if (topology[i].x==-1) {
			fprintf(stderr,"Incomplete CPU topology\n");
			return -1;
		}
	}

	/* Now find the size of a subdomain. */

	for (i=0;i<ncpus;i++) {
		if (topology[i].x > maxx) maxx=topology[i].x;
		if (topology[i].y > maxy) maxy=topology[i].y;
		if (topology[i].z > maxz) maxz=topology[i].z;
	}

	sx = nx/(1+maxx);
	sy = ny/(1+maxy);
	sz = nz/(1+maxz);

	if (nx!=(sx*(1+maxx))) {
		fprintf(stderr,"Strange number of CPUs in X direction.\n");
		return -1;
	}
	if (ny!=(sy*(1+maxy))) {
		fprintf(stderr,"Strange number of CPUs in Y direction.\n");
		return -1;
	}
	if (nz!=(sz*(1+maxz))) {
		fprintf(stderr,"Strange number of CPUs in Z direction.\n");
		return -1;
	}

	return 0;

}


int reassemble_field(datatype type, char *path, char *prefix, int t, char *suffix)
{
	char *infilename=NULL,*outfilename=NULL;
	unsigned int inflen,outflen;
	int rank,i;

	double *bigarr1=NULL,*bigarr2=NULL,*bigarr3=NULL;
	double *arr1=NULL,*arr2=NULL,*arr3=NULL;
	double *p1,*p2,*p3;

	unsigned int x,y,z,bx,by,bz;
	FILE *fldfile;
	unsigned int veclen=0;


	switch(type) {
		case TYPE_VECTOR:
			if (0==veclen) veclen=3;
			if (NULL==(bigarr1=malloc(sizeof(double)*3*nx*ny*nz))) {
				perror("malloc()");
				return -1;
			}
			if (NULL==(arr1=malloc(sizeof(double)*3*sx*sy*sz))) {
				perror("malloc()");
				free(bigarr1);
				return -1;
			}
			break;
		case TYPE_3SCALAR:
			if (0==veclen) veclen=3;
			if (NULL==(bigarr3=malloc(sizeof(double)*nx*ny*nz))) {
				perror("malloc()");
				return -1;
			}
			if (NULL==(arr3=malloc(sizeof(double)*sx*sy*sz))) {
				perror("malloc()");
				free(bigarr3);
				return -1;
			}
			/* Fall through */
		case TYPE_2SCALAR:
			if (0==veclen) veclen=2;
			if (NULL==(bigarr2=malloc(sizeof(double)*nx*ny*nz))) {
				perror("malloc()");
				if (NULL!=bigarr3) { free(bigarr3);free(arr3);}
				return -1;
			}
			if (NULL==(arr2=malloc(sizeof(double)*sx*sy*sz))) {
				perror("malloc()");
				if (NULL!=bigarr3) { free(bigarr3);free(arr3);}
				free(bigarr2);
				return -1;
			}
			/* Fall through */
		case TYPE_SCALAR:
			if (0==veclen) veclen=1;
			if (NULL==(bigarr1=malloc(sizeof(double)*nx*ny*nz))) {
				perror("malloc()");
				if (NULL!=bigarr3) { free(bigarr3);free(arr3);}
				if (NULL!=bigarr2) { free(bigarr2);free(arr2);}
				return -1;
			}
			if (NULL==(arr1=malloc(sizeof(double)*sx*sy*sz))) {
				perror("malloc()");
				if (NULL!=bigarr3) { free(bigarr3);free(arr3);}
				if (NULL!=bigarr2) { free(bigarr2);free(arr2);}
				free(bigarr1);
				return -1;
			}
			break;
		default:
			fprintf(stderr,"Bad type passed to reassemble_field()\n");
			return -1;
	}

	/* Set up the filename buffers */

	/*
	 * colour_ch_t002500_p0001.bin
	 * `---------'`----'`'`--'`--'
	 *       |     |     |  |    \- suffix
	 *       |     |     |   \----- CPUCHARS   
	 *       |     |      \-------- 2 extra chars      
	 *       |      \------------ TIMECHARS       
	 *        \------------------ prefix               
	 *                            (plus NUL at end)
	 *
	 * Transforms to:
	 *
	 * colour_ch_t002500.xdr
	 *
	 */

	inflen=strlen(path)+strlen(prefix)+TIMECHARS+2+CPUCHARS+strlen(suffix)+1;
	outflen=strlen(path)+strlen(prefix)+TIMECHARS+4+1;

	if (NULL==(infilename=malloc(inflen))) {
		perror("malloc()");
		if (NULL!=bigarr3) { free(bigarr3);free(arr3);}
		if (NULL!=bigarr2) { free(bigarr2);free(arr2);}
		free(bigarr1); free(arr1);
		return -1;
	}
	if (NULL==(outfilename=malloc(outflen))) {
		perror("malloc()");
		if (NULL!=bigarr3) { free(bigarr3);free(arr3);}
		if (NULL!=bigarr2) { free(bigarr2);free(arr2);}
		free(bigarr1); free(arr1);
		free(infilename);
		return -1;
	}

	/* Read in the single-CPU BIN files and reassemble */

	for (rank=0;rank<ncpus;rank++) {

		snprintf(infilename,inflen,
			"%s%s%0"TIMECHARSTR"d_p%0"CPUCHARSTR"d%s",
			path,prefix,t,rank,suffix);

		printf("Reading \"%s\"...",infilename);
		fflush(stdout); /* Otherwise console can get out-of-order */

		if (0!=read_bin(type,arr1,arr2,arr3,sx,sy,sz,infilename)) {
			goto bailout; /* Dijkstra probably (etc) */
		}

		printf(" OK\n");

		p1=arr1;
		p2=arr2;
		p3=arr3;
		for (z=0;z<sz;z++) {
		 for (y=0;y<sy;y++) {
		  for (x=0;x<sx;x++) {
			bx=x+(sx*topology[rank].x);
			by=y+(sy*topology[rank].y);
			bz=z+(sz*topology[rank].z);
			switch(type) {
				case (TYPE_VECTOR):
				 *(bigarr1+(bx*3)+(nx*3)*by+(nx*3)*ny*bz)=*p1++;
				 *(bigarr1+(bx*3)+(nx*3)*by+(nx*3)*ny*bz+1)=*p1++;
				 *(bigarr1+(bx*3)+(nx*3)*by+(nx*3)*ny*bz+2)=*p1++;
				break;
				case TYPE_3SCALAR:
				*(bigarr3 + bx + nx*by + nx*ny*bz)=*p3++;
				case TYPE_2SCALAR:
				*(bigarr2 + bx + nx*by + nx*ny*bz)=*p2++;
				case TYPE_SCALAR:
				*(bigarr1 + bx + nx*by + nx*ny*bz)=*p1++;
				break;
				default:
					fprintf(stderr,"1==2 in reassemble!\n");
					goto bailout;
			}
		  }
		 }
		}
	}


	snprintf(outfilename,outflen,"%s%s%0"TIMECHARSTR"d.xdr",path,prefix,t);
	printf("Writing \"%s\"...",outfilename);
	fflush(stdout);

	if (0!=write_xdr(type,bigarr1,bigarr2,bigarr3,nx,ny,nz,outfilename)) {
		goto bailout;
	}
	putchar('.');
	snprintf(outfilename,outflen,"%s%s%0"TIMECHARSTR"d.fld",path,prefix,t);

	if (NULL==(fldfile=fopen(outfilename,"w"))) {
		perror("fopen()");
		goto bailout;
	}

	fprintf(fldfile,"# AVS\nndim = 3\n");
	fprintf(fldfile,"dim1 = %d\ndim2 = %d\ndim3 = %d\n",nx,ny,nz);
	fprintf(fldfile,
		"nspace = 3\nveclen = %d\ndata = xdr_double\nfield = uniform\n",
		veclen);
	for (i=1;i<=veclen;i++) {
		fprintf(fldfile,
		 "variable %d file=%s"
		 "%0"TIMECHARSTR"d.xdr filetype=binary skip=%d stride=%u\n",
		 i,prefix,t,8*(i-1),veclen);
	}
	fclose(fldfile);

	printf(" OK\n");
	if (NULL!=bigarr3) { free(bigarr3);free(arr3);}
	if (NULL!=bigarr2) { free(bigarr2);free(arr2);}
	free(bigarr1); free(arr1);
	free(infilename);
	free(outfilename);
	return 0;


bailout:
	if (NULL!=bigarr3) { free(bigarr3);free(arr3);}
	if (NULL!=bigarr2) { free(bigarr2);free(arr2);}
	free(bigarr1); free(arr1);
	free(infilename);
	free(outfilename);

	return -1;

	
}


int reassemble_field_floatxdr(
   datatype type, char *path, char *prefix, int t, char *suffix
   )
{
	char *infilename=NULL,*outfilename=NULL;
	unsigned int inflen,outflen;
	int rank,i;

	float *bigarr1=NULL,*bigarr2=NULL,*bigarr3=NULL;
	float *arr1=NULL,*arr2=NULL,*arr3=NULL;
	float *p1,*p2,*p3;

	unsigned int x,y,z,bx,by,bz;
	FILE *fldfile;
	unsigned int veclen=0;


	switch(type) {
		case TYPE_VECTOR:
			if (0==veclen) veclen=3;
			if (NULL==(bigarr1=malloc(sizeof(float)*3*nx*ny*nz))) {
				perror("malloc()");
				return -1;
			}
			if (NULL==(arr1=malloc(sizeof(float)*3*sx*sy*sz))) {
				perror("malloc()");
				free(bigarr1);
				return -1;
			}
			break;
		case TYPE_3SCALAR:
			if (0==veclen) veclen=3;
			if (NULL==(bigarr3=malloc(sizeof(float)*nx*ny*nz))) {
				perror("malloc()");
				return -1;
			}
			if (NULL==(arr3=malloc(sizeof(float)*sx*sy*sz))) {
				perror("malloc()");
				free(bigarr3);
				return -1;
			}
			/* Fall through */
		case TYPE_2SCALAR:
			if (0==veclen) veclen=2;
			if (NULL==(bigarr2=malloc(sizeof(float)*nx*ny*nz))) {
				perror("malloc()");
				if (NULL!=bigarr3) { free(bigarr3);free(arr3);}
				return -1;
			}
			if (NULL==(arr2=malloc(sizeof(float)*sx*sy*sz))) {
				perror("malloc()");
				if (NULL!=bigarr3) { free(bigarr3);free(arr3);}
				free(bigarr2);
				return -1;
			}
			/* Fall through */
		case TYPE_SCALAR:
			if (0==veclen) veclen=1;
			if (NULL==(bigarr1=malloc(sizeof(float)*nx*ny*nz))) {
				perror("malloc()");
				if (NULL!=bigarr3) { free(bigarr3);free(arr3);}
				if (NULL!=bigarr2) { free(bigarr2);free(arr2);}
				return -1;
			}
			if (NULL==(arr1=malloc(sizeof(float)*sx*sy*sz))) {
				perror("malloc()");
				if (NULL!=bigarr3) { free(bigarr3);free(arr3);}
				if (NULL!=bigarr2) { free(bigarr2);free(arr2);}
				free(bigarr1);
				return -1;
			}
			break;
		default:
			fprintf(stderr,"Bad type passed to reassemble_field()\n");
			return -1;
	}

	/* Set up the filename buffers */

	/*
	 * colour_ch_t002500_p0001.bin
	 * `---------'`----'`'`--'`--'
	 *       |     |     |  |    \- suffix
	 *       |     |     |   \----- CPUCHARS   
	 *       |     |      \-------- 2 extra chars      
	 *       |      \------------ TIMECHARS       
	 *        \------------------ prefix               
	 *                            (plus NUL at end)
	 *
	 * Transforms to:
	 *
	 * colour_ch_t002500.xdr
	 *
	 */

	inflen=strlen(path)+strlen(prefix)+TIMECHARS+2+CPUCHARS+strlen(suffix)+1;
	outflen=strlen(path)+strlen(prefix)+TIMECHARS+4+1;

	if (NULL==(infilename=malloc(inflen))) {
		perror("malloc()");
		if (NULL!=bigarr3) { free(bigarr3);free(arr3);}
		if (NULL!=bigarr2) { free(bigarr2);free(arr2);}
		free(bigarr1); free(arr1);
		return -1;
	}
	if (NULL==(outfilename=malloc(outflen))) {
		perror("malloc()");
		if (NULL!=bigarr3) { free(bigarr3);free(arr3);}
		if (NULL!=bigarr2) { free(bigarr2);free(arr2);}
		free(bigarr1); free(arr1);
		free(infilename);
		return -1;
	}

	/* Read in the single-CPU XDR files and reassemble */

	for (rank=0;rank<ncpus;rank++) {

		snprintf(infilename,inflen,
			"%s%s%0"TIMECHARSTR"d_p%0"CPUCHARSTR"d%s",
			path,prefix,t,rank,suffix);

		printf("Reading \"%s\"...",infilename);
		fflush(stdout); /* Otherwise console can get out-of-order */

		if (0!=read_xdr_float(type,arr1,arr2,arr3,sx,sy,sz,infilename)) {
			goto bailout; /* Dijkstra probably (etc) */
		}

		printf(" OK\n");

		p1=arr1;
		p2=arr2;
		p3=arr3;
		for (z=0;z<sz;z++) {
		 for (y=0;y<sy;y++) {
		  for (x=0;x<sx;x++) {
			bx=x+(sx*topology[rank].x);
			by=y+(sy*topology[rank].y);
			bz=z+(sz*topology[rank].z);
			switch(type) {
				case (TYPE_VECTOR):
				 *(bigarr1+(bx*3)+(nx*3)*by+(nx*3)*ny*bz)=*p1++;
				 *(bigarr1+(bx*3)+(nx*3)*by+(nx*3)*ny*bz+1)=*p1++;
				 *(bigarr1+(bx*3)+(nx*3)*by+(nx*3)*ny*bz+2)=*p1++;
				break;
				case TYPE_3SCALAR:
				*(bigarr3 + bx + nx*by + nx*ny*bz)=*p3++;
				case TYPE_2SCALAR:
				*(bigarr2 + bx + nx*by + nx*ny*bz)=*p2++;
				case TYPE_SCALAR:
				*(bigarr1 + bx + nx*by + nx*ny*bz)=*p1++;
				break;
				default:
					fprintf(stderr,"1==2 in reassemble!\n");
					goto bailout;
			}
		  }
		 }
		}
	}


	snprintf(outfilename,outflen,"%s%s%0"TIMECHARSTR"d.xdr",path,prefix,t);
	printf("Writing \"%s\"...",outfilename);
	fflush(stdout);

	if (0!=write_xdr_float(type,bigarr1,bigarr2,bigarr3,nx,ny,nz,outfilename)) {
		goto bailout;
	}
	putchar('.');
	snprintf(outfilename,outflen,"%s%s%0"TIMECHARSTR"d.fld",path,prefix,t);

	if (NULL==(fldfile=fopen(outfilename,"w"))) {
		perror("fopen()");
		goto bailout;
	}

	fprintf(fldfile,"# AVS\nndim = 3\n");
	fprintf(fldfile,"dim1 = %d\ndim2 = %d\ndim3 = %d\n",nx,ny,nz);
	fprintf(fldfile,
		"nspace = 3\nveclen = %d\ndata = xdr_float\nfield = uniform\n",
		veclen);
	for (i=1;i<=veclen;i++) {
		fprintf(fldfile,
		 "variable %d file=%s"
		 "%0"TIMECHARSTR"d.xdr filetype=binary skip=%d stride=%u\n",
		 i,prefix,t,8*(i-1),veclen);
	}
	fclose(fldfile);

	printf(" OK\n");
	if (NULL!=bigarr3) { free(bigarr3);free(arr3);}
	if (NULL!=bigarr2) { free(bigarr2);free(arr2);}
	free(bigarr1); free(arr1);
	free(infilename);
	free(outfilename);
	return 0;


bailout:
	if (NULL!=bigarr3) { free(bigarr3);free(arr3);}
	if (NULL!=bigarr2) { free(bigarr2);free(arr2);}
	free(bigarr1); free(arr1);
	free(infilename);
	free(outfilename);

	return -1;

	
}


int reassemble_field_xdr(
   datatype type, char *path, char *prefix, int t, char *suffix
   )
{
	char *infilename=NULL,*outfilename=NULL;
	unsigned int inflen,outflen;
	int rank,i;

	double *bigarr1=NULL,*bigarr2=NULL,*bigarr3=NULL;
	double *arr1=NULL,*arr2=NULL,*arr3=NULL;
	double *p1,*p2,*p3;

	unsigned int x,y,z,bx,by,bz;
	FILE *fldfile;
	unsigned int veclen=0;


	switch(type) {
		case TYPE_VECTOR:
			if (0==veclen) veclen=3;
			if (NULL==(bigarr1=malloc(sizeof(double)*3*nx*ny*nz))) {
				perror("malloc()");
				return -1;
			}
			if (NULL==(arr1=malloc(sizeof(double)*3*sx*sy*sz))) {
				perror("malloc()");
				free(bigarr1);
				return -1;
			}
			break;
		case TYPE_3SCALAR:
			if (0==veclen) veclen=3;
			if (NULL==(bigarr3=malloc(sizeof(double)*nx*ny*nz))) {
				perror("malloc()");
				return -1;
			}
			if (NULL==(arr3=malloc(sizeof(double)*sx*sy*sz))) {
				perror("malloc()");
				free(bigarr3);
				return -1;
			}
			/* Fall through */
		case TYPE_2SCALAR:
			if (0==veclen) veclen=2;
			if (NULL==(bigarr2=malloc(sizeof(double)*nx*ny*nz))) {
				perror("malloc()");
				if (NULL!=bigarr3) { free(bigarr3);free(arr3);}
				return -1;
			}
			if (NULL==(arr2=malloc(sizeof(double)*sx*sy*sz))) {
				perror("malloc()");
				if (NULL!=bigarr3) { free(bigarr3);free(arr3);}
				free(bigarr2);
				return -1;
			}
			/* Fall through */
		case TYPE_SCALAR:
			if (0==veclen) veclen=1;
			if (NULL==(bigarr1=malloc(sizeof(double)*nx*ny*nz))) {
				perror("malloc()");
				if (NULL!=bigarr3) { free(bigarr3);free(arr3);}
				if (NULL!=bigarr2) { free(bigarr2);free(arr2);}
				return -1;
			}
			if (NULL==(arr1=malloc(sizeof(double)*sx*sy*sz))) {
				perror("malloc()");
				if (NULL!=bigarr3) { free(bigarr3);free(arr3);}
				if (NULL!=bigarr2) { free(bigarr2);free(arr2);}
				free(bigarr1);
				return -1;
			}
			break;
		default:
			fprintf(stderr,"Bad type passed to reassemble_field()\n");
			return -1;
	}

	/* Set up the filename buffers */

	/*
	 * colour_ch_t002500_p0001.bin
	 * `---------'`----'`'`--'`--'
	 *       |     |     |  |    \- suffix
	 *       |     |     |   \----- CPUCHARS   
	 *       |     |      \-------- 2 extra chars      
	 *       |      \------------ TIMECHARS       
	 *        \------------------ prefix               
	 *                            (plus NUL at end)
	 *
	 * Transforms to:
	 *
	 * colour_ch_t002500.xdr
	 *
	 */

	inflen=strlen(path)+strlen(prefix)+TIMECHARS+2+CPUCHARS+strlen(suffix)+1;
	outflen=strlen(path)+strlen(prefix)+TIMECHARS+4+1;

	if (NULL==(infilename=malloc(inflen))) {
		perror("malloc()");
		if (NULL!=bigarr3) { free(bigarr3);free(arr3);}
		if (NULL!=bigarr2) { free(bigarr2);free(arr2);}
		free(bigarr1); free(arr1);
		return -1;
	}
	if (NULL==(outfilename=malloc(outflen))) {
		perror("malloc()");
		if (NULL!=bigarr3) { free(bigarr3);free(arr3);}
		if (NULL!=bigarr2) { free(bigarr2);free(arr2);}
		free(bigarr1); free(arr1);
		free(infilename);
		return -1;
	}

	/* Read in the single-CPU XDR files and reassemble */

	for (rank=0;rank<ncpus;rank++) {

		snprintf(infilename,inflen,
			"%s%s%0"TIMECHARSTR"d_p%0"CPUCHARSTR"d%s",
			path,prefix,t,rank,suffix);

		printf("Reading \"%s\"...",infilename);
		fflush(stdout); /* Otherwise console can get out-of-order */

		if (0!=read_xdr(type,arr1,arr2,arr3,sx,sy,sz,infilename)) {
			goto bailout; /* Dijkstra probably (etc) */
		}

		printf(" OK\n");

		p1=arr1;
		p2=arr2;
		p3=arr3;
		for (z=0;z<sz;z++) {
		 for (y=0;y<sy;y++) {
		  for (x=0;x<sx;x++) {
			bx=x+(sx*topology[rank].x);
			by=y+(sy*topology[rank].y);
			bz=z+(sz*topology[rank].z);
			switch(type) {
				case (TYPE_VECTOR):
				 *(bigarr1+(bx*3)+(nx*3)*by+(nx*3)*ny*bz)=*p1++;
				 *(bigarr1+(bx*3)+(nx*3)*by+(nx*3)*ny*bz+1)=*p1++;
				 *(bigarr1+(bx*3)+(nx*3)*by+(nx*3)*ny*bz+2)=*p1++;
				break;
				case TYPE_3SCALAR:
				*(bigarr3 + bx + nx*by + nx*ny*bz)=*p3++;
				case TYPE_2SCALAR:
				*(bigarr2 + bx + nx*by + nx*ny*bz)=*p2++;
				case TYPE_SCALAR:
				*(bigarr1 + bx + nx*by + nx*ny*bz)=*p1++;
				break;
				default:
					fprintf(stderr,"1==2 in reassemble!\n");
					goto bailout;
			}
		  }
		 }
		}
	}


	snprintf(outfilename,outflen,"%s%s%0"TIMECHARSTR"d.xdr",path,prefix,t);
	printf("Writing \"%s\"...",outfilename);
	fflush(stdout);

	if (0!=write_xdr(type,bigarr1,bigarr2,bigarr3,nx,ny,nz,outfilename)) {
		goto bailout;
	}
	putchar('.');
	snprintf(outfilename,outflen,"%s%s%0"TIMECHARSTR"d.fld",path,prefix,t);

	if (NULL==(fldfile=fopen(outfilename,"w"))) {
		perror("fopen()");
		goto bailout;
	}

	fprintf(fldfile,"# AVS\nndim = 3\n");
	fprintf(fldfile,"dim1 = %d\ndim2 = %d\ndim3 = %d\n",nx,ny,nz);
	fprintf(fldfile,
		"nspace = 3\nveclen = %d\ndata = xdr_double\nfield = uniform\n",
		veclen);
	for (i=1;i<=veclen;i++) {
		fprintf(fldfile,
		 "variable %d file=%s"
		 "%0"TIMECHARSTR"d.xdr filetype=binary skip=%d stride=%u\n",
		 i,prefix,t,8*(i-1),veclen);
	}
	fclose(fldfile);

	printf(" OK\n");
	if (NULL!=bigarr3) { free(bigarr3);free(arr3);}
	if (NULL!=bigarr2) { free(bigarr2);free(arr2);}
	free(bigarr1); free(arr1);
	free(infilename);
	free(outfilename);
	return 0;


bailout:
	if (NULL!=bigarr3) { free(bigarr3);free(arr3);}
	if (NULL!=bigarr2) { free(bigarr2);free(arr2);}
	free(bigarr1); free(arr1);
	free(infilename);
	free(outfilename);

	return -1;

	
}




/*
 *
 *           path         prefix
 *     /~~~~~~~~~~~~~~~\/~~~~~~~~~\
 *     ../output/folder/colour_ch_t0123_p0456.bin
 *     `-------' `----'
 *      pathstub  folder
 *
 */



int postprocess_parallel(char *cmdfilename, int sglxdr, int dblxdr) {

	char *dotinputfile,*ifilename,*p;
	nmlist	nl;
	char *folder,*gr_out_file;
	char *path; unsigned int pathlen;
	char *topfname; unsigned int topflen;
	char *prefix; unsigned int preflen;

	int t;
	double radius=0.0;
	double com[3]={0,0,0};


	if (0!=strcmp(cmdfilename,".input-file")){ 
         ifilename = skip_whitespace(cmdfilename);
        } else { 
	 if (NULL==(dotinputfile=readfile(".input-file"))) {
		fprintf(stderr,"Could not read \".input-file\"\n");
		return -1;
	 }
	 ifilename=skip_whitespace(dotinputfile);
        } 

	p=skip_to_whitespace(ifilename);
	if (0!=*p) {
		*p=0;
	}

	/* ifilename now points to the null-term name of the input file.
	 * Try to parse it.
	 *
	 */
	
	/*if (0!=nmlist_parse_file(&nl,dotinputfile)) {     Jens 28/05/02*/
	if (0!=nmlist_parse_file(&nl,ifilename)) {
		fprintf(stderr,"Unable to parse input file \"%s\"\n",
			ifilename);
		free(dotinputfile);
		return -1;
	}

	/* The input file parsed OK. */

	free(dotinputfile);

	/* Now fetch values of interest. */

	/* Integers */

	if (0!=nmlist_lookup_int(&nl,"fixed_input","nx",&nx)) {
		fprintf(stderr,"Error: nx undefined in input file.\n");
		return -1;
	}
	if (0!=nmlist_lookup_int(&nl,"fixed_input","ny",&ny)) {
		fprintf(stderr,"Error: ny undefined in input file.\n");
		return -1;
	}
	if (0!=nmlist_lookup_int(&nl,"fixed_input","nz",&nz)) {
		fprintf(stderr,"Error: nz undefined in input file.\n");
		return -1;
	}
	if (0!=nmlist_lookup_int(&nl,"variable_input","n_iteration",
		&n_iteration)) {
		fprintf(stderr,"Error: n_iteration undefined in input file.\n");
		return -1;
	}
	if (0!=nmlist_lookup_int(&nl,"variable_input","n_sci",&n_sci)) {
		fprintf(stderr,"Error: n_sci undefined in input file.\n");
		return -1;
	}
	if (0!=nmlist_lookup_int(&nl,"variable_input","sci_start",&sci_start)) {
		sci_start=0;
	}
	if (0!=nmlist_lookup_int(&nl,"variable_input","init_cond",&init_cond)) {
	        fprintf(stderr,
			"Error: init_cond undefined in input file.\n");
		return -1;
	}

	/* Floats */
	if (0!=nmlist_lookup_double(&nl,"variable_input","fr1",&fr1)) {
	   fprintf(stderr,
		   "Error: fr1 undefined in input file.\n");
	   return -1;
	}
	if (0!=nmlist_lookup_double(&nl,"variable_input","fr2",&fr2)) {
	   fprintf(stderr,
		   "Error: fr2 undefined in input file.\n");
	   return -1;
	}

	/* Booleans */

	if (-1==(sci_int=nmlist_lookup_bool(&nl,"variable_input","sci_int"))) {
		sci_int=FALSE;
	}
	if (-1==(sci_sur=nmlist_lookup_bool(&nl,"variable_input","sci_sur"))) {
		sci_sur=FALSE;
	}
	if (-1==(sci_vel=nmlist_lookup_bool(&nl,"variable_input","sci_vel"))) {
		sci_vel=FALSE;
	}
	if (-1==(sci_flo=nmlist_lookup_bool(&nl,"variable_input","sci_flo"))) {
		sci_flo=FALSE;
	}
	if (-1==(sci_od=nmlist_lookup_bool(&nl,"variable_input","sci_od"))) {
		sci_od=FALSE;
	}
	if (-1==(sci_wd=nmlist_lookup_bool(&nl,"variable_input","sci_wd"))) {
		sci_wd=FALSE;
	}
	if (-1==(sci_arrows=nmlist_lookup_bool(&nl,"variable_input","sci_arrows"))){
		sci_arrows=FALSE;
	}
	if (-1==(sci_dir=nmlist_lookup_bool(&nl,"variable_input","sci_dir"))){
		printf("sci_dir not found\n");
		sci_dir=FALSE;
	}
	printf("sci_dir=%d\n",sci_dir);

	/* Strings */

	if (NULL==(folder=nmlist_lookup(&nl,"variable_input",
			"folder"))) {
		fprintf(stderr,"Error: folder undefined in input file\n");
		return -1;
	}
	if (NULL==(gr_out_file=nmlist_lookup(&nl,"variable_input",
			"gr_out_file"))) {
		fprintf(stderr,"Error: gr_out_file undefined in input file\n");
		return -1;
	}

	pathlen = strlen(pathstub)+1+strlen(folder)+1+1;
	if (NULL==(path=malloc(pathlen))) { 
		perror("malloc()");
		return -1;
	}
	snprintf(path,pathlen,"%s/%s/",pathstub,folder);

	topflen=strlen(pathstub)+1+strlen(folder)+1+strlen("coords_")
		+strlen(gr_out_file)+1;
	if (NULL==(topfname=malloc(topflen))) {
		perror("malloc()");
		return -1;
	}
	snprintf(topfname,topflen,"%s/%s/coords_%s",pathstub,folder,gr_out_file);

	preflen=strlen("heregoesthelongeststring_")+strlen(gr_out_file)+2+1;
	if (NULL==(prefix=malloc(preflen))) {
		perror("malloc()");
		return -1;
	}

	if (0!=find_topology(topfname)) {
		fprintf(stderr,"Can't find topology file \"%s\""
				"- unable to postprocess\n",topfname);
		free(path);
		free(prefix);
		free(topfname);
		return -1;
	}
 
 	/* Print out parameters:	*/
	printf("\nsci_int=%d",sci_int);
	printf("\nsci_sur=%d",sci_sur);
	printf("\nsci_vel=%d",sci_vel);
	printf("\nsci_flo=%d",sci_flo);
	printf("\nsci_arrows=%d",sci_arrows);
	printf("\nsci_dir=%d",sci_dir);
	printf("\nsci_od=%d",sci_od);
	printf("\nsci_wd=%d",sci_wd);
	printf("\nfolder=%s",folder);

	for (t=0;t<=n_iteration;t++) {
	   if ( (t>=sci_start) && (0==(t%n_sci)) ) {

	      
	      if (sci_int) {
		 snprintf(prefix,preflen,"colour_%s_t",gr_out_file);
		 if(sglxdr)
		    reassemble_field_floatxdr(TYPE_SCALAR,path,prefix,t,".xdr");
		 else if(dblxdr)
		    reassemble_field_xdr(TYPE_SCALAR,path,prefix,t,".xdr");
		 else
		    reassemble_field(TYPE_SCALAR,path,prefix,t,".bin");
	      }
	      
	      if (sci_sur) {
		 snprintf(prefix,preflen,"sur_%s_t",gr_out_file);
		 if(sglxdr)
		    reassemble_field_floatxdr(TYPE_SCALAR,path,prefix,t,".xdr");
		 else if(dblxdr)
		    reassemble_field_xdr(TYPE_SCALAR,path,prefix,t,".xdr");
		 else
		    reassemble_field(TYPE_SCALAR,path,prefix,t,".bin");
	      }
	      
	      if (sci_vel) {
		 snprintf(prefix,preflen,"vel_%s_t",gr_out_file);
		 if(sglxdr)
		    reassemble_field_floatxdr(TYPE_SCALAR,path,prefix,t,".xdr");
		 else if(dblxdr)
		    reassemble_field_xdr(TYPE_SCALAR,path,prefix,t,".xdr");
		 else
		    reassemble_field(TYPE_SCALAR,path,prefix,t,".bin");
	      }
	      
/* Old version, used when only one flo file (3 scalars) was produced */
/* In the new version two (if compiled with -DNOSURFACTANT) or three flo files (1 scalar) are produced, one for each phase) */ 
/*	      if (sci_flo) {
		 snprintf(prefix,preflen,"flo_%s_t",gr_out_file);
		 if(sglxdr)
		    reassemble_field_floatxdr(TYPE_3SCALAR,path,prefix,t,".xdr");
		 else if(dblxdr)
		    reassemble_field_xdr(TYPE_3SCALAR,path,prefix,t,".xdr");
		 else
		    reassemble_field(TYPE_3SCALAR,path,prefix,t,".bin");
	      }
*/

/* New version */
	      if (sci_flo) {
		 snprintf(prefix,preflen,"flooil_%s_t",gr_out_file);
		 if(sglxdr)
		    reassemble_field_floatxdr(TYPE_SCALAR,path,prefix,t,".xdr");
		 else if(dblxdr)
		    reassemble_field_xdr(TYPE_SCALAR,path,prefix,t,".xdr");
		 else
		    reassemble_field(TYPE_SCALAR,path,prefix,t,".bin");

		 snprintf(prefix,preflen,"flowater_%s_t",gr_out_file);
		 if(sglxdr)
		    reassemble_field_floatxdr(TYPE_SCALAR,path,prefix,t,".xdr");
		 else if(dblxdr)
		    reassemble_field_xdr(TYPE_SCALAR,path,prefix,t,".xdr");
		 else
		    reassemble_field(TYPE_SCALAR,path,prefix,t,".bin");

		 snprintf(prefix,preflen,"flosurf_%s_t",gr_out_file);
		 if(sglxdr)
		    reassemble_field_floatxdr(TYPE_SCALAR,path,prefix,t,".xdr");
		 else if(dblxdr)
		    reassemble_field_xdr(TYPE_SCALAR,path,prefix,t,".xdr");
		 else
		    reassemble_field(TYPE_SCALAR,path,prefix,t,".bin");
	      }


	      
	      if (sci_arrows) {
		 snprintf(prefix,preflen,"arr_%s_t",gr_out_file);
		 if(sglxdr)
		    reassemble_field_floatxdr(TYPE_VECTOR,path,prefix,t,".xdr");
		 else if(dblxdr)
		    reassemble_field_xdr(TYPE_VECTOR,path,prefix,t,".xdr");
		 else
		    reassemble_field(TYPE_VECTOR,path,prefix,t,".bin");
	      }
	      
	      if (sci_dir) {
		 snprintf(prefix,preflen,"dir_%s_t",gr_out_file);
		 if(sglxdr)
		    reassemble_field_floatxdr(TYPE_VECTOR,path,prefix,t,".xdr");
		 else if(dblxdr)
		    reassemble_field_xdr(TYPE_VECTOR,path,prefix,t,".xdr");
		 else
		    reassemble_field(TYPE_VECTOR,path,prefix,t,".bin");
	      }
	      
	      
	      if (sci_od) {
		 snprintf(prefix,preflen,"od_%s_t",gr_out_file);
		 if(sglxdr)
		    reassemble_field_floatxdr(TYPE_SCALAR,path,prefix,t,".xdr");
		 else if(dblxdr)
		    reassemble_field_xdr(TYPE_SCALAR,path,prefix,t,".xdr");
		 else
		    reassemble_field(TYPE_SCALAR,path,prefix,t,".bin");
	      }


	      if (sci_wd) {
		 snprintf(prefix,preflen,"wd_%s_t",gr_out_file);
		 if(sglxdr)
		    reassemble_field_floatxdr(TYPE_SCALAR,path,prefix,t,".xdr");
		 else if(dblxdr)
		    reassemble_field_xdr(TYPE_SCALAR,path,prefix,t,".xdr");
		 else
		    reassemble_field(TYPE_SCALAR,path,prefix,t,".bin");
	      }


	      /*
		if (sci_spress) {
		snprintf(prefix,preflen,"spressure_%s_t",gr_out_file);
		reassemble_field(TYPE_SCALAR,path,prefix,t,".bin");
		}
	      */
	      if((init_cond == 8)||
		 (init_cond == 3)||
		 (init_cond == -2)||
		 (init_cond == 4)) {
		 snprintf(prefix,preflen,"pxx_%s_t",gr_out_file);
		 if(sglxdr)
		    reassemble_field_floatxdr(TYPE_SCALAR,path,prefix,t,".xdr");
		 else if(dblxdr)
		    reassemble_field_xdr(TYPE_SCALAR,path,prefix,t,".xdr");
		 else
		    reassemble_field(TYPE_SCALAR,path,prefix,t,".bin");
		 snprintf(prefix,preflen,"pyy_%s_t",gr_out_file);
		 if(sglxdr)
		    reassemble_field_floatxdr(TYPE_SCALAR,path,prefix,t,".xdr");
		 else if(dblxdr)
		    reassemble_field_xdr(TYPE_SCALAR,path,prefix,t,".xdr");
		 else
		    reassemble_field(TYPE_SCALAR,path,prefix,t,".bin");
		 snprintf(prefix,preflen,"pzz_%s_t",gr_out_file);
		 if(sglxdr)
		    reassemble_field_floatxdr(TYPE_SCALAR,path,prefix,t,".xdr");
		 else if(dblxdr)
		    reassemble_field_xdr(TYPE_SCALAR,path,prefix,t,".xdr");
		 else
		    reassemble_field(TYPE_SCALAR,path,prefix,t,".bin");
		 snprintf(prefix,preflen,"pxy_%s_t",gr_out_file);
		 if(sglxdr)
		    reassemble_field_floatxdr(TYPE_SCALAR,path,prefix,t,".xdr");
		 else if(dblxdr)
		    reassemble_field_xdr(TYPE_SCALAR,path,prefix,t,".xdr");
		 else
		    reassemble_field(TYPE_SCALAR,path,prefix,t,".bin");
		 snprintf(prefix,preflen,"pyz_%s_t",gr_out_file);
		 if(sglxdr)
		    reassemble_field_floatxdr(TYPE_SCALAR,path,prefix,t,".xdr");
		 else if(dblxdr)
		    reassemble_field_xdr(TYPE_SCALAR,path,prefix,t,".xdr");
		 else
		    reassemble_field(TYPE_SCALAR,path,prefix,t,".bin");
		 snprintf(prefix,preflen,"pxz_%s_t",gr_out_file);
		 if(sglxdr)
		    reassemble_field_floatxdr(TYPE_SCALAR,path,prefix,t,".xdr");
		 else if(dblxdr)
		    reassemble_field_xdr(TYPE_SCALAR,path,prefix,t,".xdr");
		 else
		    reassemble_field(TYPE_SCALAR,path,prefix,t,".bin");
		 /*Next commented out lines were there  */
		 /*for debugging purposes.		*/
		 
		 /* snprintf(prefix,preflen,"popz_%s_t",gr_out_file);	*/
		 /* reassemble_field(TYPE_SCALAR,path,prefix,t,".bin");	*/
		 /* snprintf(prefix,preflen,"popx_%s_t",gr_out_file);	*/
		 /* reassemble_field(TYPE_SCALAR,path,prefix,t,".bin");	*/
		 /* snprintf(prefix,preflen,"p1xx_%s_t",gr_out_file);	*/
		 /* reassemble_field(TYPE_SCALAR,path,prefix,t,".bin");	*/
		 /* snprintf(prefix,preflen,"p1yy_%s_t",gr_out_file);	*/
		 /* reassemble_field(TYPE_SCALAR,path,prefix,t,".bin");	*/
		 /* snprintf(prefix,preflen,"p1zz_%s_t",gr_out_file);	*/
		 /* reassemble_field(TYPE_SCALAR,path,prefix,t,".bin");	*/
		 /* snprintf(prefix,preflen,"p2xx_%s_t",gr_out_file);	*/
		 /* reassemble_field(TYPE_SCALAR,path,prefix,t,".bin");	*/
		 /* snprintf(prefix,preflen,"p2yy_%s_t",gr_out_file);	*/	
		 /* reassemble_field(TYPE_SCALAR,path,prefix,t,".bin");	*/
		 /* snprintf(prefix,preflen,"p2zz_%s_t",gr_out_file);	*/
		 /* reassemble_field(TYPE_SCALAR,path,prefix,t,".bin"); */
	      }
	      if(init_cond == 1) {
		 /* This case is not finished yet,			*/
		 /* since osmotic_pressure() hasn't been		*/
		 /* adapted to the new SGLXDR case. --NGS, Apr2003.  	*/

		 snprintf(prefix,preflen,"od_%s_t",gr_out_file);
		 snprintf(prefix,preflen,"wd_%s_t",gr_out_file);
		 actual_radius(path,prefix,t,&radius,com);
		 snprintf(prefix,preflen,"scp_%s_t",gr_out_file);
		 if(sglxdr)
		    reassemble_field_floatxdr(TYPE_SCALAR,path,prefix,t,".xdr");
		 else if(dblxdr)
		    reassemble_field_xdr(TYPE_SCALAR,path,prefix,t,".xdr");
		 else
		    reassemble_field(TYPE_SCALAR,path,prefix,t,".bin");
		 osmotic_pressure(radius,com,path,prefix,t);
	      }


	      /********************************/

	      if(t==0) {
		 snprintf(prefix,preflen,"checkpoint_%s_t",gr_out_file);
		 reassemble_field(TYPE_SCALAR,path,prefix,t,".bin");
	      }		 

	      /********************************/

	   }
	}

	printf("\n\bDone!\n");
	return 0;
}




int
osmotic_pressure(double radius, double com[], char *path, char *prefix, int t)
{
   double         *spress=NULL, *a=NULL;
   double         osmotic_pressure=0.0;
   char           suffix[]="osmotic";


   if(NULL==(spress=malloc(sizeof(double)*nx*ny*nz))) {
      perror("malloc() in post.c:osmotic_pressure()");
      exit(-1);
   }

   decode_xdr(spress,a,TYPE_SCALAR,path,prefix,t);

   average_pressure(spress,radius,com,&osmotic_pressure,t);

   dump_scalar(osmotic_pressure,path,prefix,suffix,t);

   return 0;
}


double
min(double a, double b, double c)
{
   double min;

   min = a;
   if(b < min) min = b;
   if(c < min) min = c;

   return min;
}


int
decode_xdr(double *a1, double *a2, datatype type, 
	   char *path, char *prefix, int t)
{
   double         *a=NULL;
   unsigned int   inflen;
   char           *infilename=NULL;

   inflen=strlen(path)+strlen(prefix)+TIMECHARS+4+1;

   if (NULL==(infilename=malloc(inflen))) {
      perror("malloc()");
      return -1;
   }

   snprintf(infilename,inflen,"%s%s%0"TIMECHARSTR"d.xdr",path,prefix,t);
   printf("Reading \"%s\"...",infilename);
   fflush(stdout);

   switch(type) {
      case TYPE_SCALAR:

	 if(0!=read_xdr(TYPE_SCALAR,a1,a2,a,nx,ny,nz,infilename)) {
	    perror("post.c:decode_xdr():read_xdr()");
	    return -1;
	 }
	 printf(" OK\n");
	 break;

      case TYPE_2SCALAR:

	 if(0!=read_xdr(TYPE_2SCALAR,a1,a2,a,nx,ny,nz,infilename)) {
	    perror("post.c:decode_xdr():read_xdr()");
	    return -1;
	 }
	 printf(" OK\n");
	 break;
   }

   /*  Array spress is stored as */	
   /*  .xdr in Fortran order	 */
   /*  i.e., spress[x][y][z] 	 */
   /*  with x the fastest.       */

   return 0;
}


int
average_pressure(double *spress, double radius, double com[],
		 double *osmotic_pressure, int t)
{
   int       offset[3];
   int       x,y,z;
   double    spress_in=0.0, spress_out=0.0;
   double    if1, if2, rad, maxrad=0.0;
   double    in, out;
   double    field;
   FILE      *aux;

   /* spress_in  = -DBL_MAX; */
   /* spress_out = DBL_MAX; */
   in  = 0.0;
   out = 0.0;

   /* This first part as taken from
    * [lbe-lite code] lbe_init.F90:lbe_init_radial()
    *
    * fr1, fr2 and nx,ny,nz are global vars.
    *
    *
   if1 = fr1 * min((double)nx,(double)ny,(double)nz)/2.;
   if2 = fr2 * min((double)nx,(double)ny,(double)nz)/2.;
   offset[0] = -nx/2;
   offset[1] = -ny/2;
   offset[2] = -nz/2;

   */   
   

   /*  This new part takes the actual radius
    *  and centre of mass, rather than assuming 
    *  a droplet centred at the lattice's centre
    *  that keeps its initial radius.
    *
    */
   if1 = if2 = radius;
   offset[0] = -com[0];
   offset[1] = -com[1];
   offset[2] = -com[2];


   /*
   aux = fopen("../output/droplet1/radii.dat", "w");
   fprintf(aux, "%d\n", t);
   */

   for(z = 0; z < nz; z++) {
      for(y = 0; y < ny; y++) {
	 for(x = 0; x < nx; x++) {
	    
	    rad = sqrt(SQUARE(x+offset[0])+
	               SQUARE(y+offset[1])+
	               SQUARE(z+offset[2]));
	    field = *(spress + x + nx*y + nx*ny*z);

	    if(rad > maxrad) 
	       maxrad = rad;

	    /*  BEST METHOD FOR MEASURING PRESSURE:
	     *  Maximum pressure inside and
	     *  average pressure outside.
	     *  Best because it's the easiest way to
	     *  define the pressure inside and outside
	     *  considering pressure doesn't keep
	     *  constant in space.
	     *
	     */

	    if(rad <= WIN_IN * if1) { 
/*	       spress_in += field;
	       in++;
*/    
	       /* Determine the maximum scalar
		* pressure inside the
		* droplet. 
		* spress_in >= 0.0 always.
		*/

	       if(field > spress_in) 
		  spress_in = field;
  
	    } else if(rad >= WIN_OUT * if2) {
	       spress_out += field;
	       out++;
	       
               /* Determine the maximum scalar
		* pressure outside the
		* droplet. 
		* spress_out >= 0.0 always.
		*/
/*
               if(field < spress_out) 
		  spress_out = field;
*/
	       
	    }
            /* FOR DEBUGGING ONLY
	     *
	     */	    
/*	    if(t == 800) {
	       fprintf(aux,"%d\t%d\t%d\t%16.8f\t%16.8f\t%16.8f\t%16.8f\n",
		       x,y,z,rad,field,spress_in,spress_out);
	       fflush(aux);
	    }
*/

	 }
      }
   }

/*
   if(in!=0.0)
      spress_in /= in;
*/
   if(out!=0.0)
      spress_out /= out;


   *osmotic_pressure = spress_in - spress_out; 

   /* FOR DEBUGGING ONLY
    *
    *
   printf("\nt                = %d"
          "\nif1              = %16.8f"   
          "\nif2              = %16.8f"
	  "\nCentre of mass   = (%d,%d,%d)"
          "\nMax radius       = %16.8f"
          "\nin               = %16.8f"
	  "\nout              = %16.8f"
	  "\nspress_in        = %16.8f"
	  "\nspress_out       = %16.8f"
	  "\nosmotic_pressure = %16.8f", 
	  t,if1,if2,(int)com[0],(int)com[1],(int)com[2],
	  maxrad,in,out,spress_in,
	  spress_out,*osmotic_pressure);


   fclose(aux);
*/


   return 0;
}


int
dump_scalar(double scalar, char *path, char *prefix, char *suffix, int t)
{
   int        i,outflen,newstrlen;
   char       *outfilename, *newprefix;
   FILE       *fp;

   newstrlen = strlen(prefix)-2;
   if (NULL==(newprefix=malloc(newstrlen))) {
      perror("malloc()");
      return -1;
   }

   /* Remove trailing "-t" from prefix */
   for(i=0; i < newstrlen; i++)
      newprefix[i] = prefix[i];
   /*
   fprintf(stdout,"\npost.c:dump_osmotic(): newprefix=%s",newprefix);
   */
   fflush(stdout);

   outflen=strlen(path)+newstrlen+1+strlen(suffix)+4+1;
   /*
    * Format of osmotic pressure file is:
    *
    *           path       newprefix  suffix
    *     /~~~~~~~~~~~~~~~\/~~~~~~~\ /~~~~~\
    *     ../output/folder/colour_ch_osmotic.dat
    *     `-------' `----'
    *      pathstub  folder
    * 
    * +1 comes from connector _ between newprefix and suffix
    * +4 comes from the extension .dat
    * +1 comes from the \0 character
    */

   if (NULL==(outfilename=malloc(outflen))) {
      perror("malloc()");
      return -1;
   }

   snprintf(outfilename,outflen,"%s%s_%s.dat",path,newprefix,suffix);
   printf("\nWriting \"%s\"...",outfilename);
   fflush(stdout);

   if(t == 0) {
      if(remove(outfilename) == 0) {
         printf("\nAlready existing %s removed...\n",
                outfilename);
      }
   }
   if((fp = fopen(outfilename, "a"))==NULL) {
      printf("\nFile %s can't be opened\n", outfilename);
      exit(0);
   }
   fprintf(fp, "%d\t%16.8f\n", t, scalar);
   fclose(fp);

   printf(" OK\n");
   
   return 0;
}


double
min_field(double *field)
{
   int    x,y,z;
   double value;
   double min_field=DBL_MAX;

   for(z=0; z<nz; z++) {
      for(y=0; y<ny; y++) {
	 for(x=0; x<nx; x++) {
	    
	    value = *(field + x + y*nx + z*nx*ny);
	    if(value < min_field)
	       min_field = value;
	    
	 }
      }
   }
   return min_field;
}


double
max_field(double *field)
{
   int    x,y,z;
   double value;
   double max_field=-DBL_MAX;

   for(z=0; z<nz; z++) {
      for(y=0; y<ny; y++) {
	 for(x=0; x<nx; x++) {
	    
	    value = *(field + x + y*nx + z*nx*ny);
	    if(value > max_field)
	       max_field = value;
		       
	 }
      }
   }
   return max_field;
}


double
avg_field(double *field)
{
   int    x,y,z;
   double avg_field=0.0;

   for(z=0; z<nz; z++) {
      for(y=0; y<ny; y++) {
	 for(x=0; x<nx; x++) {
	    
	    avg_field += *(field + x + y*nx + z*nx*ny);
	
	 }
      }
   }
   return (avg_field / ((double)nx*ny*nz));
}



int
compute_radius(double *density, double *radius)
{
   double  max_dens, min_dens, avg_dens;

   max_dens = max_field(density);
   min_dens = min_field(density);
   avg_dens = avg_field(density);

   if(max_dens - min_dens == 0.0) {
      perror("Divide-by-zero in post.c:compute_radius()");
      exit(-1);
   }
   *radius = (avg_dens - min_dens)/(max_dens - min_dens);
   *radius *= 3 / (4 * PI);
   *radius = pow(*radius,1./3.) * nx;

   printf("\nmax_dens = %16.8f"
	  "\nmin_dens = %16.8f"
	  "\navg_dens = %16.8f"
          "\nradius   = %16.8f", 
	  max_dens,min_dens,avg_dens,*radius);

   return 0;
}


int
centre_of_mass(double *density, double com[])
{
   int    x,y,z;
   double dens, tot_dens;

   com[0] = 0.0;
   com[1] = 0.0;
   com[2] = 0.0;
   tot_dens = 0.0;

   for(z=0; z<nz; z++) {
      for(y=0; y<ny; y++) {
	 for(x=0; x<nx; x++) {
	  
	    dens      = *(density + x + y*nx + z*nx*ny);
	    com[0]   += dens * x;
	    com[1]   += dens * y;
	    com[2]   += dens * z;
	    tot_dens += dens;
	
	 }
      }
   }

   com[0] = floor(com[0]/tot_dens);
   com[1] = floor(com[1]/tot_dens);
   com[2] = floor(com[2]/tot_dens);

   return 0;
}


int
actual_radius(char *path, char *prefix, int t, double *radius, double com[])
{
   double   *red_dens=NULL, *blue_dens=NULL;
   char     suffix[]="radius";

   if(NULL==(red_dens=malloc(sizeof(double)*nx*ny*nz))) {
      perror("malloc() in post.c:actual_radius()");
      exit(-1);
   }
   if(NULL==(blue_dens=malloc(sizeof(double)*nx*ny*nz))) {
      perror("malloc() in post.c:actual_radius()");
      exit(-1);
   }

   decode_xdr(red_dens,blue_dens,TYPE_2SCALAR,path,prefix,t);

   compute_radius(red_dens,radius);

   centre_of_mass(red_dens,com);

   dump_scalar(*radius,path,prefix,suffix,t);


   return 0;
}





int 
main (int argc, char *argv[])
{
   int sglxdr;
   int dblxdr;
  
   if(2==argc) {
      /* Assume argument is <DBLBIN | SGLXDR | DBLXDR> */
      if(0==strcmp(argv[1],"DBLBIN")){
	 sglxdr=0;
	 dblxdr=0;}
      else if(0==strcmp(argv[1],"SGLXDR")){
	 sglxdr=1;
	 dblxdr=0;}
      else if(0==strcmp(argv[1],"DBLXDR")){
	 sglxdr=0;
	 dblxdr=1;}
      else 
	 goto syntax;
      return postprocess_parallel(".input-file",sglxdr,dblxdr);
   } 

   else if((4==argc) && (0==strcmp(argv[2],"-f"))) {
      if(0==strcmp(argv[1],"DBLBIN")){
         sglxdr=0;
	 dblxdr=0;}
      else if(0==strcmp(argv[1],"SGLXDR")){
         sglxdr=1;
	 dblxdr=0;}
      else if(0==strcmp(argv[1],"DBLXDR")){
	 sglxdr=0;
	 dblxdr=1;}
      return postprocess_parallel(argv[3],sglxdr,dblxdr);
   } 

   else {
syntax:
      fprintf(stderr,"Syntax: post <DBLBIN | SGLXDR | DBLXDR> "
	     	 "[-f input-file]\n"
		 "where you have to specify whether the data is\n"
		 "double precision Fortran binary or single precision XDR or double precision XDR\n");
      return -1;
   }
}



