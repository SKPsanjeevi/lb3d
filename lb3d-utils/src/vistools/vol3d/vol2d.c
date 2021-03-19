#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <netinet/in.h>

#include "config.h"
#include "vol2d.h"

/*

=head1 NAME

vol2d.c - 2-dimensional volume library

=head1 SYNOPSIS

 #include "vol2d.h"

 struct vol2d *v;
 int x,y;
 int nx=16,ny=32;
 
 if (NULL==(v=vol2d_new(16,32,VOL_TYPE_UINT8))) {
 	fprintf(stderr,"vol2d_new() failed\n");
 }
 
 for (y=0;y<ny;y++) {
 for (x=0;x<nx;x++) {
 	*vol2d_element(v,x,y) = x+y;
 }}
 
 vol2d_write_raw(v,"mydata.bin");
 
 vol2d_destroy(v);

=head1 DESCRIPTION

  This library deals with homogeneous two-dimensional arrays of
  elements of arbitrary size. It allows these arrays to be read
  from and written to files, and to be manipulated through
  streaming transfer functions.

=cut

=head1 ROUTINES

=over 1

=cut

*/



/*

=item B<vol2d_get_datatype_size>

  size_t vol2d_get_datatype_size(vol_datatype datatype)
  
  Returns the size, in bytes, of a single element of the given
  vol2d datatype. Returns zero if VOL_TYPE_UNKNOWN or an invalid
  type are given.

=cut

*/

size_t vol2d_get_datatype_size(vol_datatype datatype)
{
	switch(datatype)
	{
		case VOL_TYPE_UNKNOWN:
			fprintf(stderr,"WARNING: vol2d_get_datatype_size"
					"(VOL_TYPE_UNKNOWN)\n");
			return 0;
		case VOL_TYPE_UINT8:
			return 1;
		case VOL_TYPE_INT8:
			return 1;
		case VOL_TYPE_UINT16:
			return 2;
		case VOL_TYPE_INT16:
			return 2;
		case VOL_TYPE_UINT32:
			return 4;
		case VOL_TYPE_INT32:
			return 4;
		case VOL_TYPE_UINT64:
			return 8;
		case VOL_TYPE_INT64:
			return 8;
		case VOL_TYPE_FLOAT:
			return VOL_TYPE_FLOATSIZE;
		case VOL_TYPE_FLOAT2:
			return 2*VOL_TYPE_FLOATSIZE;
		case VOL_TYPE_FLOAT3:
			return 3*VOL_TYPE_FLOATSIZE;
		case VOL_TYPE_DOUBLE:
			return VOL_TYPE_DOUBLESIZE;
		case VOL_TYPE_DOUBLE2:
			return 2*VOL_TYPE_DOUBLESIZE;
		case VOL_TYPE_DOUBLE3:
			return 3*VOL_TYPE_DOUBLESIZE;
		default:
			fprintf(stderr,"WARNING: vol2d_get_datatype_size"
					"(unknown type %d)\n",datatype);
			return 0;
	}
}

/*

=item B<vol2d_new>

  struct vol2d *vol2d_new(
		unsigned int nx, unsigned int ny, 
		vol_datatype datatype)

  v = vol2d_new(12,34,VOL_TYPE_UINT8);

This routine allocates a new vol2d structure of specified size,
to contain elements of given type. The element data is uninitialized.

Returns C<NULL> on failure.

=cut

*/

struct vol2d *vol2d_new(
		unsigned int nx, unsigned int ny, 
		vol_datatype datatype)
{
	struct vol2d *v=NULL;
	size_t elsize=0;

	if (VOL_TYPE_UNKNOWN==datatype) {
		fprintf(stderr,"WARNING: vol2d_new(VOL_TYPE_UNKNOWN)\n");
		return NULL;
	} else {
		if (0==(elsize=vol2d_get_datatype_size(datatype))) {
			fprintf(stderr,"vol2d_new: vol2d_get_datatype_size=0\n"
			       );
			return NULL;
		}
	}

	if (NULL==(v=vol2d_new_elsize(nx,ny,elsize))) {
		fprintf(stderr,"vol2d_new: vol2d_new_elsize failed\n");
		return NULL;
	}
	v->datatype = datatype;
	return v;
	
}

/*

=item B<vol2d_new_elsize>

  struct vol2d *vol2d_new_elsize(
		unsigned int nx, unsigned int ny,
		size_t elsize)

  v = vol2d_new(12,34,1);

This routine allocates a new vol2d structure of specified size,
to contain elements of size C<elsize>. The element data is
uninitialized, and the volume data type is set to
C<VOL_TYPE_UNKNOWN>.

Returns C<NULL> on failure.

=cut

*/

struct vol2d *vol2d_new_elsize(
		unsigned int nx, unsigned int ny,
		size_t elsize)
{
	struct vol2d *v=NULL;

	if (NULL==(v=malloc(sizeof(struct vol2d)))) {
		perror("malloc()");
		return NULL;
	}

	if (NULL==(v->data=malloc(elsize*nx*ny))) {
		perror("malloc()");
		free(v);
		return NULL;
	}
	v->nx=nx;
	v->ny=ny;
	v->elsize=elsize;
	v->datatype = VOL_TYPE_UNKNOWN;
	return v;
}
/*

=item B<vol2d_new_copy>

  struct vol2d *vol2d_new_copy(struct vol2d *vold)

  vnew = vol2d_new_copy(vold);

This routine creates a new vol2d object, and sets its size, datatype,
and elements to be identical to the original one, vold.

Returns C<NULL> on failure.

=cut

*/

struct vol2d *vol2d_new_copy(struct vol2d *vold)
{
	struct vol2d *v=NULL;

	if (NULL==(v=vol2d_new_copytype(vold))) {
		fprintf(stderr,"vol2d_new_copy: vol2d_new_copytype failed\n");
		return NULL;
	}

	memcpy(v->data,vold->data,v->nx*v->ny*v->elsize);
	return v;
}


/*

=item B<vol2d_new_copytype>

  struct vol2d *vol2d_new_copytype(struct vol2d *vold)

  vnew = vol2d_new_copytype(vold);

This routine creates a new vol2d object, and sets its size and datatype
to be identical to the original one, but does not initialize any
elements.

Returns C<NULL> on failure.

=cut

*/
struct vol2d *vol2d_new_copytype(struct vol2d *vold)
{
	struct vol2d *v=NULL;

	if (NULL==(v=vol2d_new_elsize(vold->nx,vold->ny,
					vold->elsize))) {
		fprintf(stderr,"vol2d_new_copy: vol2d_new_elsize failed\n");
		return NULL;
	}
	v->datatype = vold->datatype;
	return v;
}
/*

=item B<vol2d_new_transfer_elsize>

  struct vol2d *vol2d_new_transfer(
	struct vol2d *vold,
	size_t elsize,
	void (*transferfunc)(void *, void *, void *),
	void *data
	);

  struct vol2d *v1,*v2;
  int nx=16,ny=16;


  void myfunc(void *src, void *dest, void *data) {
  	double outval;
	unsigned char inval;

	inval = *((unsigned char*)src);
	outval = (double) inval;
	*( (double *)dest) = outval;
  }

  ...
  
  v1 = vol2d_new(nx,ny,VOL_TYPE_UINT8);

  vol2d_read_raw(v1,"uchar.dat");
  v2 = vol2d_new_transfer_elsize(v1,sizeof(double),myfunc);


This routine creates a new vol2d object, of the same
dimensions as the given object. Each element in the new object
is created by passing the corresponding element of the old
object through a transfer function. The C<data> argument is passed
straight through as the third argument to the transfer function.

Returns C<NULL> on failure.

=cut

*/

struct vol2d *vol2d_new_transfer_elsize(
	struct vol2d *vold,
	size_t elsize,
	void (*transferfunc)(void *, void *, void *),
	void *data
	) {

	struct vol2d *v=NULL;
	int i;
	void *src,*dst;

	if (NULL==(v=vol2d_new_elsize(vold->nx,vold->ny,elsize))) {
		fprintf(stderr,"vol2d_new_transfer: vol2d_new_elsize failed\n");
		return NULL;
	}

	src = vold->data; dst=v->data;

	for (i=0;i<v->nx*v->ny;i++) {
		transferfunc(src,dst,data);
		src =(void *)
			((char *)src + vold->elsize);
		dst =(void *)
			((char *)dst + elsize);
	}
	return v;
}

/*

=item B<vol2d_new_transfer>

  struct vol2d *vol2d_new_transfer(
  	struct vol2d *vold,
  	vol_datatype datatype,
  	void (*transferfunc)(void *, void *, void *),
  	void *data


This routine creates a new vol2d object, of the same dimensions as the
given object, with elements of type C<datatype>. Each element in the new
object is created by passing the corresponding element of the old object
through a transfer function. The C<data> argument is passed straight
through as the third argument to the transfer function.

Returns C<NULL> on failure.

=cut

*/

struct vol2d *vol2d_new_transfer(
	struct vol2d *vold,
	vol_datatype datatype,
	void (*transferfunc)(void *, void *, void *),
	void *data
	) {

	struct vol2d *v=NULL;
	size_t elsize=0;

	if (VOL_TYPE_UNKNOWN==datatype) {
		fprintf(stderr,"vol2d_new_transfer(VOL_TYPE_UNKNOWN)\n");
		return NULL;
	} else {
		if (0==(elsize=vol2d_get_datatype_size(datatype))) {
			fprintf(stderr,"vol2d_new_transfer: "
					"vol2d_get_datatype_size=0\n"
			       );
			return NULL;
		}
	}

	if (NULL==(v=vol2d_new_transfer_elsize(vold,
					elsize,transferfunc,data))){
		fprintf(stderr,"vol2d_new_transfer: "
				"vol2d_new_transfer_elsize failed\n");
		return NULL;
	}

	v->datatype = datatype;
	return v;
}

/*

=item B<vol2d_destroy>

  int vol2d_destroy(struct vol2d *v);

  vol2d_destroy(v);

This routine deallocated the given vol2d object.
Returns zero on success, nonzero on failure.

=cut

*/

int vol2d_destroy(struct vol2d *v)
{
	if (NULL==v) {
		fprintf(stderr,"WARNING: vol2d_destroy(NULL)\n");
		return -1;
	}
	if (NULL==v->data) {
		fprintf(stderr,"WARNING: vol2d_destroy(data=NULL)\n");
		free(v);
		return -2;
	} else {
		free(v->data);
	}
	free(v);
	return 0;
}

/*

=item B<vol2d_write_raw_fhandle>

  int vol2d_write_raw_fhandle(struct vol2d *v, FILE *f)

This routine writes the raw data in the given vol2d object to
the given open file descriptor. Returns zero on success,
nonzero on failure.

=cut

*/

int vol2d_write_raw_fhandle(struct vol2d *v, FILE *f)
{
	size_t length;

	length = v->nx*v->ny ;

	if (length!=fwrite(v->data,v->elsize,length,f)) {
		perror("fwrite()");
		return -1;
	} else {
		return 0;
	}

}
/*

=item B<vol2d_write_raw>

  int vol2d_write_raw(struct vol2d *v, char *fname)

This routine writes the raw data in the given vol2d object to
the given file. Returns zero on success, nonzero on failure.

=cut

*/

int vol2d_write_raw(struct vol2d *v, char *fname)
{
	FILE *f;
	size_t length;

	if (NULL==(f=fopen(fname,"wb"))) {
		perror("fopen()");
		return -1;
	}

	length = v->nx*v->ny ;
	if (length!=fwrite(v->data,v->elsize,length,f)) {
		perror("fwrite()");
		fclose(f);
		return -1;
	} else {
		fclose(f);
		return 0;
	}
}

/*

=item B<vol2d_read_raw_fhandle>

  int vol2d_read_raw_fhandle(struct vol2d *v, FILE *f)

This routine reads the appropriate amount of raw data
from the given open filehandle, into the given vol2d
object. Returns zero on success, nonzero on failure.

=cut

*/

int vol2d_read_raw_fhandle(struct vol2d *v, FILE *f)
{
	size_t length;


	length = v->nx*v->ny ;

	if (length!=fread(v->data,v->elsize,length,f)) {
		perror("fread()");
		return -1;
	} else {
		return 0;
	}
}

/*

=item B<vol2d_read_raw>

  int vol2d_read_raw(struct vol2d *v, char *fname)

This routine reads the appropriate amount of raw data from the
given file, into the given vol2d object. Returns zero on
success, nonzero on failure.

=cut

*/

int vol2d_read_raw(struct vol2d *v, char *fname)
{
	FILE *f;
	size_t length;

	if (NULL==(f=fopen(fname,"rb"))) {
		perror("fopen()");
		return -1;
	}

	length = v->nx*v->ny ;
	if (length!=fread(v->data,v->elsize,length,f)) {
		perror("fread()");
		fclose(f);
		return -1;
	} else {
		fclose(f);
		return 0;
	}
}

/*

=item B<vol2d_element>

  void *vol2d_element(struct vol2d *v,
        unsigned int x, unsigned int y)

Returns a pointer to the element at (x,y) in the given volume.
Coordinates are zero-based, so the first element is at (0,0,0).
Performs no bounds-checking!

=cut

*/

void *vol2d_element(struct vol2d *v,
		unsigned int x, unsigned int y)
{
	return (void *) (
		((char *)v->data)
			+ (v->elsize)*(x + y*v->nx ));
}

/*

=item B<vol2d_element_bc>

  void *vol2d_element_bc(struct vol2d *v,
        unsigned int x, unsigned int y)

Returns a pointer to the element at (x,y) in the given volume,
or NULL if x or y are out of bounds.

=cut

*/

void *vol2d_element_bc(struct vol2d *v,
		unsigned int x, unsigned int y)
{
        if ( (x>=v->nx)||(y>=v->ny)) {
                fprintf(stderr,"vol2d_element_bc:(%d,%d) in (%d,%d)\n",
                                x,y,v->nx,v->ny);
                return NULL;
        }
	return (void *) (
		((char *)v->data)
			+ (v->elsize)*(x + y*v->nx ));
}

/*

=item B<vol2d_new_subvol>

  struct vol2d *vol2d_new_subvol(struct vol2d *vold, 
                  unsigned int x0, unsigned int y0,
                  unsigned int dx, unsigned int dy)

Creates a new vol2d object of size (dx,dy), created from the 
subvolume of that size whose lower-coordinate corner is at
(x0,y0). Returns a pointer to a new vol2d object on success, 
or C<NULL> on failure.

=cut

*/

struct vol2d *vol2d_new_subvol(struct vol2d *vold, 
	unsigned int x0, unsigned int y0,
	unsigned int dx, unsigned int dy)
{
	size_t spansize; /* Number of bytes in a span of dx */

	struct vol2d *v=NULL;
	void *p=NULL;
	int y;

	/* Sanity-check arguments */

	if (
		  ((x0+dx)>vold->nx)
		||((y0+dy)>vold->ny)
	   ) {
		fprintf(stderr,"vol2d_new_subvol called with silly arguments:"
			"\nr=(%u,%u) dr=(%u,%u)\n"
			"vol has dimensions (%u,%u)\n",
			x0,y0,dx,dy,vold->nx,vold->ny);
		return NULL;
	}

	if (NULL==(v=vol2d_new_elsize(dx,dy,vold->elsize))) {
		fprintf(stderr,"vol2d_new_subvol: vol2d_new_elsize failed\n");
		return NULL;
	}
	v->datatype = vold->datatype;

	p = v->data;

	spansize = dx * vold->elsize;

	for (y=y0;y<(y0+dy);y++) {
		memcpy(p,vol2d_element(vold,x0,y),spansize);
		p = (void *)(
			((char *)p)
			+ spansize
		);
	}

	return v;


}

/*
=item B<vol2d_set>

  void vol2d_set(struct vol2d *v, int c)

This routine sets all data in C<v> to the byte value C<c>. it is useful
for zeroing a freshly-allocated volume.

=cut
*/

void vol2d_set(struct vol2d *v, int c)
{
	memset(v->data,c,v->nx*v->ny*v->elsize);
}

/*

=item B<vol2d_subvol_readstream>


  int vol2d_subvol_readstream(struct vol2d *v, 
  	unsigned int x0, unsigned int y0,
  	unsigned int dx, unsigned int dy,
	void (*streamfunc)(void *, void *),
	void *data
	)

This routine reads data from the given stream function into the
subvolume of the C<vol2d> object C<v> of size (dx,dy) which starts at
(x0,y0). The subvolume must fit entirely into C<v>.

=cut
   
*/

int vol2d_subvol_readstream(struct vol2d *v, 
	unsigned int x0, unsigned int y0,
	unsigned int dx, unsigned int dy,
	void (*streamfunc)(void *, void *),
	void *data)
{

	void *p=NULL;
	int i,y;

	/* Sanity-check arguments */

	if (
		  ((x0+dx)>v->nx)
		||((y0+dy)>v->ny)
	   ) {
		fprintf(stderr,"vol2d_subvol_readstream "
			"called with silly arguments:"
			"\nr=(%u,%u) dr=(%u,%u)\n"
			"vol has dimensions (%u,%u)\n",
			x0,y0,dx,dy,v->nx,v->ny);
		return -1;
	}

	p = v->data;

	for (y=y0;y<(y0+dy);y++) {
		p = vol2d_element(v,x0,y);
		for (i=0;i<dx;i++) {
			streamfunc(p,data);
			p = (char *)p + v->elsize;
		}
	}

	return 0;


}
/*

=item B<vol2d_set_subvol>

  struct vol2d *vol2d_set_subvol(struct vol2d *v, 
                  unsigned int x0, unsigned int y0,
                  unsigned int dx, unsigned int dy,
		  void *val)

Takes a vol2d object, and sets all elements in the specified
subvolume to copies of the single supplied element C<val>.
Returns 0 on success, nonzero on failure.

=cut

*/
int vol2d_set_subvol(struct vol2d *v, 
                unsigned int x0, unsigned int y0,
                unsigned int dx, unsigned int dy,
		void *val)

{
	size_t spansize; /* Number of bytes in a span of dx */
	void *span=NULL;

	void *p=NULL;
	int i;
	int y;

	/* Sanity-check arguments */

	if (
		  ((x0+dx)>v->nx) /* Volume too large */
		||((y0+dy)>v->ny)
		||(0==dx*dy)   /* Volume of zero size */
	   ) {
		fprintf(stderr,"vol2d_set_subvol called with silly arguments:"
			"\nr=(%u,%u) dr=(%u,%u)\n"
			"vol has dimensions (%u,%u)\n",
			x0,y0,dx,dy,v->nx,v->ny);
		return -1;
	}

	spansize = dx * v->elsize;

	if (NULL==(span=malloc(spansize))) {
		perror("vol2d_set_subvol: malloc()");
		return -1;
	}

	/* Place dx element copies into span. */
	for (i=0;i<dx;i++) {
		memcpy((char *)span+(i*v->elsize),val,v->elsize);
	}

	p = v->data;


	for (y=y0;y<(y0+dy);y++) {
		memcpy(p,span,spansize);
		p = (char *)p +  spansize;
	}

	free(span);

	return 0;

}


/*

=item B<vol2d_new_stream_elsize>

  struct vol2d *vol2d_new_stream_elsize(
  	unsigned int nx, unsigned int ny,
  	size_t elsize,
  	void (*streamfunc)(void *, void *),
	void *data)

This function creates a new C<vol2d> object of given size, containing
elements of size C<elsize>. The function C<streamfunc> is then called
once per element, in Fortran order. The first argument to C<streamfunc>
is a pointer to the element to be written; the second argument is taken
from the C<data> argument to C<vol2d_new_stream_elsize>.  On success,
returns a pointer to a new C<vol2d> object; on failure, returns C<NULL>.

  void *mystream(void *dst, void *data)
  {
  	*(uint8_t *)dst = 0xff;
  }
 
  ...
 
  struct vol2d *v;
 
  v = vol2d_new_stream_elsize(8,8,sizeof(uint8_t),mystream);

=cut

*/

struct vol2d *vol2d_new_stream_elsize(
		unsigned int nx, unsigned int ny,
		size_t elsize,
		void (*streamfunc)(void *, void *),
		void *data)
{
	struct vol2d *v=NULL;
	void *p=NULL;
	int i;

	if (NULL==(v=vol2d_new_elsize(nx,ny,elsize))) {
		fprintf(stderr,"vol2d_new_stream_elsize: "
				"vol2d_new_elsize failed\n");
		return NULL;
	}

	p = v->data;
	for (i=0;i<nx*ny;i++) {
		streamfunc(p,data);
		p = (void *) ((char *)p + elsize);
	}
	return v;
}

/*

=item B<vol2d_new_stream>

  struct vol2d *vol2d_new_stream(
  	unsigned int nx, unsigned int ny,
  	vol_datatype datatype,
  	void (*streamfunc)(void *, void *),
  	void *data)

This function creates a new C<vol2d> object of given size, containing
elements of type C<datatype>. The function C<streamfunc> is then called
once per element, in Fortran order. The first argument to C<streamfunc>
is a pointer to the element to be written; the second argument is taken
from the C<data> argument to C<vol2d_new_stream_elsize>.  On success,
returns a pointer to a new C<vol2d> object; on failure, returns C<NULL>.

  void *mystream(void *dst, void *data)
  {
  	*(uint8_t *)dst = 0xff;
  }
 
  ...
 
  struct vol2d *v;
 
  v = vol2d_new_stream_elsize(8,8,VOL_TYPE_UINT8,mystream);

=cut

*/


struct vol2d *vol2d_new_stream(
		unsigned int nx, unsigned int ny,
		vol_datatype datatype,
		void (*streamfunc)(void *, void *),
		void *data)
{
	struct vol2d *v=NULL;
	size_t elsize;

	if (VOL_TYPE_UNKNOWN==datatype) {
		fprintf(stderr,"vol2d_new_stream(VOL_TYPE_UNKNOWN)\n");
		return NULL;
	} else {
		if (0==(elsize=vol2d_get_datatype_size(datatype))) {
			fprintf(stderr,"vol2d_new_stream: "
					"vol2d_get_datatype_size=0\n"
			       );
			return NULL;
		}
	}

	if (NULL==(v=vol2d_new_stream_elsize(nx,ny,elsize,streamfunc,data))){
		fprintf(stderr,"vol2d_new_stream: "
				"vol2d_new_stream_elsize failed\n");
		return NULL;
	}
	v->datatype = datatype;

	return v;
}

/*

=item B<vol2d_stream>

  void vol2d_stream(struct vol2d *v,
  	void (*streamfunc)(void *, void *),
	void *data);

  void mystream(void *el, void *data)
  {
  	printf("%u\n",(unsigned int) *( (uint8_t *)el);
  }

  ...

  struct vol2d *v;
  v = vol2d_new(3,3,3,VOL_TYPE_UINT8);

  ...

  vol2d_stream(v,streamfunc;


This function iterates through all elements in the given C<vol2d>
object, in Fortran order, and calls C<streamfunc> with two arguments: a
pointer to the element, and the C<data> argument passed to
C<vol2d_stream>.

=cut

*/

void vol2d_stream(struct vol2d *v,
	void (*streamfunc)(void *, void *),
	void *data)

{
	unsigned int i,n;
	void *p;
	size_t elsize;

	n=v->nx*v->ny;
	elsize = v->elsize;
	p = v->data;

	for (i=0;i<n;i++) {
		streamfunc(p,data);
		p = (void *) ((char *)p + elsize);
	}
	return;
}
