#include <stdio.h>
#include <stdlib.h>

#include <inttypes.h>
#include "vol3d.h"

#define NREPS 16

#undef VOL3DTEST_MESSY

#define TESTDATA1 \
	{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23}
#define TEST1X 2
#define TEST1Y 3
#define TEST1Z 4
#define TESTDATA1SIZE (TEST1X*TEST1Y*TEST1Z)

#define declare_testarray1(dtype) \
dtype dtype ## _data1[TESTDATA1SIZE] = TESTDATA1;

declare_testarray1(uint8_t)
declare_testarray1(int8_t)
declare_testarray1(uint16_t)
declare_testarray1(int16_t)
declare_testarray1(uint32_t)
declare_testarray1(int32_t)
declare_testarray1(uint64_t)
declare_testarray1(int64_t)
declare_testarray1(float)
declare_testarray1(double)

struct vol3d *testvol1[VOL3D_NTYPES];


#define writebinfile(dtype) \
	if (NULL==(f=fopen( #dtype ".bin","w"))) { \
		perror("fopen()"); return -1; \
	} \
	fwrite( dtype ## _data1,sizeof(dtype),TESTDATA1SIZE,f); \
	fclose(f);

static int vol3dtest_setup(void)
{
	FILE *f=NULL;

	writebinfile(uint8_t)
	writebinfile(int8_t)
	writebinfile(uint16_t)
	writebinfile(int16_t)
	writebinfile(uint32_t)
	writebinfile(int32_t)
	writebinfile(uint64_t)
	writebinfile(int64_t)
	writebinfile(float)
	writebinfile(double)

	return 0;
}

#define test_read1(dtype) \
	if (0!=vol3d_read_raw(v,#dtype ".bin")) { \
		fprintf(stderr,"FAILED\n");return -1; \
	} \
	if (0!=memcmp(v->data, dtype ## _data1,sizeof( dtype ## _data1))) { \
		fprintf(stderr,"FAILED on memcmp\n"); return -1; \
	} \
	vol3d_destroy(v);

static int vol3dtest_loading(void)
{
	struct vol3d *v=NULL;

	fprintf(stderr,"Testing vol3d_read_raw...");

	v=vol3d_new(TEST1X,TEST1Y,TEST1Z,VOL3D_UINT8);
	test_read1(uint8_t)
	v=vol3d_new(TEST1X,TEST1Y,TEST1Z,VOL3D_INT8);
	test_read1(int8_t)
	v=vol3d_new(TEST1X,TEST1Y,TEST1Z,VOL3D_UINT16);
	test_read1(uint16_t)
	v=vol3d_new(TEST1X,TEST1Y,TEST1Z,VOL3D_INT16);
	test_read1(int16_t)
	v=vol3d_new(TEST1X,TEST1Y,TEST1Z,VOL3D_UINT32);
	test_read1(uint32_t)
	v=vol3d_new(TEST1X,TEST1Y,TEST1Z,VOL3D_INT32);
	test_read1(int32_t)
	v=vol3d_new(TEST1X,TEST1Y,TEST1Z,VOL3D_UINT64);
	test_read1(uint64_t)
	v=vol3d_new(TEST1X,TEST1Y,TEST1Z,VOL3D_INT64);
	test_read1(int64_t)
	
	fprintf(stderr,"OK\n");
	return 0;
}

struct streamtest_data {
	int i;
};
/* Generate a stream of int16 elements in the sequence 0, 1, 2, 3, ... */

static void newstreamtest(void *dst, void *data)
{
	struct streamtest_data *std;

	std = (struct streamtest_data *)data;
	*(int16_t *)dst = (int16_t)std->i;
	std->i++;
}

static int vol3dtest_new_stream(void)
{
	struct vol3d *v=NULL;

	struct streamtest_data std;

	std.i=0;

	fprintf(stderr,"Testing vol3d_new_stream...");

	if (NULL==(v=vol3d_new_stream(TEST1X,TEST1Y,TEST1Z,
					VOL3D_INT16,
					newstreamtest,
					&std))) {
		fprintf(stderr,"FAILED\n"); return -1;
	}

	if (0!=memcmp(v->data,int16_t_data1,sizeof(int16_t_data1))) {
		fprintf(stderr,"FAILED on memcmp\n");
		return -1;
	}
	vol3d_destroy(v);
	fprintf(stderr,"OK\n");
	return 0;

}

struct streamreader {
	int expected;
	int ok;
};

/* Function which will read in data from a uint16_t element stream,
 * and check that each element has the value given in 
 * ((streamreader *)data)->expected .
 */
static void streamreader(void *el, void *data)
{
	struct streamreader *sr=NULL;
	uint16_t elval;

	sr = (struct streamreader *)data;
	elval = *(uint16_t *)el;

	if (elval != sr->expected) {
		fprintf(stderr,"FAILED: expected %d, got %d\n",
				sr->expected,elval);
		sr->ok=0;
		return;
	}
	sr->expected++;
}

static int vol3dtest_stream(void)
{
	struct vol3d *v=NULL;
	struct streamreader sr;

	sr.expected=0;
	sr.ok=1;

	fprintf(stderr,"Testing vol3d_stream...");
	if (NULL==(v=vol3d_new(TEST1X,TEST1Y,TEST1Z,VOL3D_INT16))) {
		fprintf(stderr,"FAILED on vol3d_new\n"); return -1;
	}

	memcpy(v->data,int16_t_data1,sizeof(int16_t_data1)) ;

	vol3d_stream(v,streamreader,(void *)&sr);

	if (sr.ok) { 
		fprintf(stderr,"OK\n");
		return 0;
	} else {
		return -1;
	}

	
}

struct addfuncdata {
	struct vol3d *v1,*v2;
	int i;
};

/* Set el to be the sum of corresponding voxels in
 * v1 and v2 as referenced in addfuncdata structure
 */
static void addingstreamfunc(void *el, void *data)
{
	struct addfuncdata *mydata=NULL;
	int16_t x,y;

	mydata = (struct addfuncdata *)data;

	x = *((int16_t *) mydata->v1->data + mydata->i);
	y = *((int16_t *) mydata->v2->data + mydata->i);
	mydata->i++;
	*(int16_t *)el = x+y;

}

static int vol3dtest_streaminout(void)
{
	struct vol3d *v16a=NULL,*v16b=NULL,*v16c=NULL;
	int16_t *p;
	struct streamreader sr;
	int i;
	struct addfuncdata afd;

	sr.expected=0;
	sr.ok=1;


	fprintf(stderr,"Testing vol3d_stream plumbing...");
	if (NULL==(v16a=vol3d_new(TEST1X,TEST1Y,TEST1Z,VOL3D_INT16))) {
		fprintf(stderr,"FAILED on vol3d_new\n"); return -1;
	}

	memcpy(v16a->data,int16_t_data1,sizeof(int16_t_data1)) ;

	if (NULL==(v16b=vol3d_new_copy(v16a))) {
		fprintf(stderr,"FAILED on vol3d_new_copy\n");
		return -1;
	}

	p=v16b->data;
	for (i=0;i<TESTDATA1SIZE;i++) {
		int16_t phi;
		phi = *p;
		*p++ = 1-phi;
	}

	afd.i=0;
	afd.v1 = v16a;
	afd.v2 = v16b;

	if (NULL==(v16c=vol3d_new_stream(TEST1X,TEST1Y,TEST1Z,
					VOL3D_INT16,
					addingstreamfunc,
					&afd))) {
		fprintf(stderr,"FAILED on vol3d_new_stream\n");
		return -1;
	}

	p=v16c->data;
	for (i=0;i<TESTDATA1SIZE;i++) {
		if (1!=*p++) {
			fprintf(stderr,"FAILED at element %d:"
				"got %d, expected 1\n",
				i,*(p-1));
			return -1;
		}
	}

	vol3d_destroy(v16a);
	vol3d_destroy(v16b);
	vol3d_destroy(v16c);

	fprintf(stderr,"OK\n");
	return 0;
}

static int vol3dtest_new_copy(void)
{
	struct vol3d *v=NULL,*vc=NULL;

	fprintf(stderr,"Testing vol3d_new_copy...");
	if (NULL==(v=vol3d_new(TEST1X,TEST1Y,TEST1Z,VOL3D_INT16))) {
		fprintf(stderr,"FAILED on vol3d_new\n"); return -1;
	}

	memcpy(v->data,int16_t_data1,sizeof(int16_t_data1)) ;

	if (NULL==(vc=vol3d_new_copy(v))) {
		fprintf(stderr,"FAILED\n"); return -1;
	}

	if (
		(v->nx!=vc->nx) ||
		(v->ny!=vc->ny) ||
		(v->nz!=vc->nz) ) {
		fprintf(stderr,"FAILED on dimensions check\n");
		return -1;
	}

	if (v->datatype != vc->datatype) {
		fprintf(stderr,"FAILED on datatype check\n");
		return -1;
	}
	if (v->elsize != vc->elsize) {
		fprintf(stderr,"FAILED on elsize check\n");
		return -1;
	}

	if (0!=memcmp(v->data,vc->data,sizeof(int16_t_data1))) {
		fprintf(stderr,"FAILED data verification\n");
		return -1;
	}

	if ((0!=vol3d_destroy(vc))||(0!=vol3d_destroy(v))) {
		fprintf(stderr,"FAILED on destruction\n"); return -1;
	}
	fprintf(stderr,"OK\n");
	return 0;

}

static int vol3dtest_creation(void)
{
	struct vol3d *v=NULL;
	int i;

	fprintf(stderr,"Testing creation and destruction...");
	fflush(stderr);
	for (i=0;i<NREPS;i++) {
		if (NULL==(v=vol3d_new(12,34,56,VOL3D_UINT8))) {
			fprintf(stderr,"FAILED on vol3d_new");
			return -1;
		}
		if (0!=vol3d_destroy(v)) {
			fprintf(stderr,"FAILED on destroy");
			return -1;
		}
	}
	fprintf(stderr,"OK\n");
	return 0;
}

static int vol3dtest_creationtypes(void)
{
	struct vol3d *v=NULL;
	int i;
	int nx=3,ny=4,nz=5;

	fprintf(stderr,"Testing creation types...");
	fflush(stderr);

#ifdef VOL3DTEST_MESSY
	if (NULL!=vol3d_new(nx,ny,nz,0)) {
		fprintf(stderr,"FAILED: vol3d_new works on datatype=0\n");
		return -1;
	}
#endif /* VOL3DTEST_MESSY */

	for (i=1;i<VOL3D_NTYPES;i++) {
		if (NULL==(v=vol3d_new(nx,ny,nz,i))) {
			fprintf(stderr,"FAILED: vol3d_new(type=%d)\n",i);
			return -1;
		}
		if (i!=v->datatype) {
			fprintf(stderr,"FAILED: vol3d_new: bad datatype\n");
			return -1;
		}
		if (vol3d_get_datatype_size(i) != v->elsize) {
			fprintf(stderr,"FAILED: vol3d_new: bad elsize\n");
			return -1;
		}
		if ((nx!=v->nx)||(ny!=v->ny)||(nz!=v->nz)) {
			fprintf(stderr,"FAILED: vol3d_new: bad dimensions\n");
			return -1;
		}
		if (0!=vol3d_destroy(v)) {
			fprintf(stderr,"FAILED on vol3d_destroy()\n");
			return -1;
		}
	}
	fprintf(stderr,"OK\n");
	return 0;
}


static int (*testfunc[])() = { 
	vol3dtest_setup,
	vol3dtest_creation,
	vol3dtest_creationtypes,
	vol3dtest_loading,
	vol3dtest_new_stream,
	vol3dtest_new_copy,
	vol3dtest_stream,
	vol3dtest_streaminout,
NULL } ;

static int vol3d_test(void)
{
	int (*myfunc)();
	int i;
	int nfailed=0;

	fprintf(stderr,"\tTesting vol3d methods...\n");
	myfunc=testfunc[0];
	for (i=0;NULL!=myfunc;i++) {
		if (0!=myfunc()) {
			nfailed++;
		}
		myfunc=testfunc[i+1];
	}

	/* temp visual tests */
	
	return nfailed;
}


int main(int argc, char *argv[])
{
	fflush(stderr);
	return vol3d_test();
}
