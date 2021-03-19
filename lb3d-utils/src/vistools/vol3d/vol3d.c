#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <netinet/in.h>

#include "vol3d.h"

/*

=head1 NAME

vol3d.c - 3-dimensional volume library

=head1 SYNOPSIS

 #include "vol3d.h"

 struct vol3d *v;
 int x,y,z;
 int nx=16,ny=32,nz=64;
 
 if (NULL==(v=vol3d_new(16,32,24,VOL_TYPE_UINT8))) {
 	fprintf(stderr,"vol3d_new() failed\n");
 }
 
 for (z=0;z<nz;z++) {
 for (y=0;y<ny;y++) {
 for (x=0;x<nx;x++) {
 	*vol3d_element(v,x,y,z) = x+y+z;
 }}}
 
 vol3d_write_raw(v,"mydata.bin");
 
 vol3d_destroy(v);

=head1 DESCRIPTION

  This library deals with homogeneous three-dimensional arrays of
  elements of arbitrary size. It allows these arrays to be read
  from and written to files, and to be manipulated through
  streaming transfer functions.

=cut

=head1 ROUTINES

=over 1

=cut

*/



/*

=item B<vol3d_get_datatype_size>

  size_t vol3d_get_datatype_size(vol_datatype datatype)
  
  Returns the size, in bytes, of a single element of the given
  vol3d datatype. Returns zero if VOL_TYPE_UNKNOWN or an invalid
  type are given.

=cut

*/

size_t vol3d_get_datatype_size(vol_datatype datatype)
{
	switch(datatype)
	{
		case VOL_TYPE_UNKNOWN:
			fprintf(stderr,"WARNING: vol3d_get_datatype_size"
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
			fprintf(stderr,"WARNING: vol3d_get_datatype_size"
					"(unknown type %d)\n",datatype);
			return 0;
	}
}

/*

=item B<vol3d_new>

  struct vol3d *vol3d_new(
		unsigned int nx, unsigned int ny, unsigned int nz,
		vol_datatype datatype)

  v = vol3d_new(12,34,56,VOL_TYPE_UINT8);

This routine allocates a new vol3d structure of specified size,
to contain elements of given type. The element data is uninitialized.

Returns C<NULL> on failure.

=cut

*/

struct vol3d *vol3d_new(
		unsigned int nx, unsigned int ny, unsigned int nz,
		vol_datatype datatype)
{
	struct vol3d *v=NULL;
	size_t elsize=0;

	if (VOL_TYPE_UNKNOWN==datatype) {
		fprintf(stderr,"WARNING: vol3d_new(VOL_TYPE_UNKNOWN)\n");
		return NULL;
	} else {
		if (0==(elsize=vol3d_get_datatype_size(datatype))) {
			fprintf(stderr,"vol3d_new: vol3d_get_datatype_size=0\n"
			       );
			return NULL;
		}
	}

	if (NULL==(v=vol3d_new_elsize(nx,ny,nz,elsize))) {
		fprintf(stderr,"vol3d_new: vol3d_new_elsize failed\n");
		return NULL;
	}
	v->datatype = datatype;
	return v;
	
}

/*

=item B<vol3d_new_elsize>

  struct vol3d *vol3d_new_elsize(
		unsigned int nx, unsigned int ny, unsigned int nz,
		size_t elsize)

  v = vol3d_new(12,34,56,1);

This routine allocates a new vol3d structure of specified size,
to contain elements of size C<elsize>. The element data is
uninitialized, and the volume data type is set to
C<VOL_TYPE_UNKNOWN>.

Returns C<NULL> on failure.

=cut

*/

struct vol3d *vol3d_new_elsize(
		unsigned int nx, unsigned int ny, unsigned int nz,
		size_t elsize)
{
	struct vol3d *v=NULL;

	if (NULL==(v=malloc(sizeof(struct vol3d)))) {
		perror("malloc()");
		return NULL;
	}

	if (NULL==(v->data=malloc(elsize*nx*ny*nz))) {
		perror("malloc()");
		free(v);
		return NULL;
	}
	v->nx=nx;
	v->ny=ny;
	v->nz=nz;
	v->elsize=elsize;
	v->datatype = VOL_TYPE_UNKNOWN;
	return v;
}
/*

=item B<vol3d_new_copy>

  struct vol3d *vol3d_new_copy(struct vol3d *vold)

  vnew = vol3d_new_copy(vold);

This routine creates a new vol3d object, and sets its size, datatype,
and elements to be identical to the original one, vold.

Returns C<NULL> on failure.

=cut

*/

struct vol3d *vol3d_new_copy(struct vol3d *vold)
{
	struct vol3d *v=NULL;

	if (NULL==(v=vol3d_new_copytype(vold))) {
		fprintf(stderr,"vol3d_new_copy: vol3d_new_copytype failed\n");
		return NULL;
	}

	memcpy(v->data,vold->data,v->nx*v->ny*v->nz*v->elsize);
	return v;
}


/*

=item B<vol3d_new_copytype>

  struct vol3d *vol3d_new_copytype(struct vol3d *vold)

  vnew = vol3d_new_copytype(vold);

This routine creates a new vol3d object, and sets its size and datatype
to be identical to the original one, but does not initialize any
elements.

Returns C<NULL> on failure.

=cut

*/
struct vol3d *vol3d_new_copytype(struct vol3d *vold)
{
	struct vol3d *v=NULL;

	if (NULL==(v=vol3d_new_elsize(vold->nx,vold->ny,vold->nz,
					vold->elsize))) {
		fprintf(stderr,"vol3d_new_copy: vol3d_new_elsize failed\n");
		return NULL;
	}
	v->datatype = vold->datatype;
	return v;
}
/*

=item B<vol3d_new_transfer_elsize>

  struct vol3d *vol3d_new_transfer(
	struct vol3d *vold,
	size_t elsize,
	void (*transferfunc)(void *, void *, void *),
	void *data
	);

  struct vol3d *v1,*v2;
  int nx=16,ny=16,nz=32;


  void myfunc(void *src, void *dest, void *data) {
  	double outval;
	unsigned char inval;

	inval = *((unsigned char*)src);
	outval = (double) inval;
	*( (double *)dest) = outval;
  }

  ...
  
  v1 = vol3d_new(nx,ny,nz,VOL_TYPE_UINT8);

  vol3d_read_raw(v1,"uchar.dat");
  v2 = vol3d_new_transfer_elsize(v1,sizeof(double),myfunc);


This routine creates a new vol3d object, of the same
dimensions as the given object. Each element in the new object
is created by passing the corresponding element of the old
object through a transfer function. The C<data> argument is passed
straight through as the third argument to the transfer function.

Returns C<NULL> on failure.

=cut

*/

struct vol3d *vol3d_new_transfer_elsize(
	struct vol3d *vold,
	size_t elsize,
	void (*transferfunc)(void *, void *, void *),
	void *data
	) {

	struct vol3d *v=NULL;
	int i;
	void *src,*dst;

	if (NULL==(v=vol3d_new_elsize(vold->nx,vold->ny,vold->nz,elsize))) {
		fprintf(stderr,"vol3d_new_transfer: vol3d_new_elsize failed\n");
		return NULL;
	}

	src = vold->data; dst=v->data;

	for (i=0;i<v->nx*v->ny*v->nz;i++) {
		transferfunc(src,dst,data);
		src =(void *)
			((char *)src + vold->elsize);
		dst =(void *)
			((char *)dst + elsize);
	}
	return v;
}

/*

=item B<vol3d_new_transfer>

  struct vol3d *vol3d_new_transfer(
  	struct vol3d *vold,
  	vol_datatype datatype,
  	void (*transferfunc)(void *, void *, void *),
  	void *data


This routine creates a new vol3d object, of the same dimensions as the
given object, with elements of type C<datatype>. Each element in the new
object is created by passing the corresponding element of the old object
through a transfer function. The C<data> argument is passed straight
through as the third argument to the transfer function.

Returns C<NULL> on failure.

=cut

*/

struct vol3d *vol3d_new_transfer(
	struct vol3d *vold,
	vol_datatype datatype,
	void (*transferfunc)(void *, void *, void *),
	void *data
	) {

	struct vol3d *v=NULL;
	size_t elsize=0;

	if (VOL_TYPE_UNKNOWN==datatype) {
		fprintf(stderr,"vol3d_new_transfer(VOL_TYPE_UNKNOWN)\n");
		return NULL;
	} else {
		if (0==(elsize=vol3d_get_datatype_size(datatype))) {
			fprintf(stderr,"vol3d_new_transfer: "
					"vol3d_get_datatype_size=0\n"
			       );
			return NULL;
		}
	}

	if (NULL==(v=vol3d_new_transfer_elsize(vold,
					elsize,transferfunc,data))){
		fprintf(stderr,"vol3d_new_transfer: "
				"vol3d_new_transfer_elsize failed\n");
		return NULL;
	}

	v->datatype = datatype;
	return v;
}

/*

=item B<vol3d_destroy>

  int vol3d_destroy(struct vol3d *v);

  vol3d_destroy(v);

This routine deallocates the given vol3d object.
Returns zero on success, nonzero on failure.

=cut

*/

int vol3d_destroy(struct vol3d *v)
{
	if (NULL==v) {
		fprintf(stderr,"WARNING: vol3d_destroy(NULL)\n");
		return -1;
	}
	if (NULL==v->data) {
		fprintf(stderr,"WARNING: vol3d_destroy(data=NULL)\n");
		free(v);
		return -2;
	} else {
		free(v->data);
	}
	free(v);
	return 0;
}

/*

=item B<vol3d_write_raw_fhandle>

  int vol3d_write_raw_fhandle(struct vol3d *v, FILE *f)

This routine writes the raw data in the given vol3d object to
the given open file descriptor. Returns zero on success,
nonzero on failure.

=cut

*/

int vol3d_write_raw_fhandle(struct vol3d *v, FILE *f)
{
	size_t length;

	length = v->nx*v->ny*v->nz ;

	if (length!=fwrite(v->data,v->elsize,length,f)) {
		perror("fwrite()");
		return -1;
	} else {
		return 0;
	}

}
/*

=item B<vol3d_write_raw>

  int vol3d_write_raw(struct vol3d *v, char *fname)

This routine writes the raw data in the given vol3d object to
the given file. Returns zero on success, nonzero on failure.

=cut

*/

int vol3d_write_raw(struct vol3d *v, char *fname)
{
	FILE *f;
	size_t length;

	if (NULL==(f=fopen(fname,"wb"))) {
		perror("fopen()");
		return -1;
	}

	length = v->nx*v->ny*v->nz ;
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

=item B<vol3d_read_raw_fhandle>

  int vol3d_read_raw_fhandle(struct vol3d *v, FILE *f)

This routine reads the appropriate amount of raw data
from the given open filehandle, into the given vol3d
object. Returns zero on success, nonzero on failure.

=cut

*/

int vol3d_read_raw_fhandle(struct vol3d *v, FILE *f)
{
	size_t length;


	length = v->nx*v->ny*v->nz ;

	if (length!=fread(v->data,v->elsize,length,f)) {
		perror("fread()");
		return -1;
	} else {
		return 0;
	}
}

/*

=item B<vol3d_read_raw>

  int vol3d_read_raw(struct vol3d *v, char *fname)

This routine reads the appropriate amount of raw data from the
given file, into the given vol3d object. Returns zero on
success, nonzero on failure.

=cut

*/

int vol3d_read_raw(struct vol3d *v, char *fname)
{
	FILE *f;
	size_t length;

	if (NULL==(f=fopen(fname,"rb"))) {
		perror("fopen()");
		return -1;
	}

	length = v->nx*v->ny*v->nz ;
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

=item B<vol3d_element>

  void *vol3d_element(struct vol3d *v,
        unsigned int x, unsigned int y, unsigned int z)

Returns a pointer to the element at (x,y,z) in the given volume.
Coordinates are zero-based, so the first element is at (0,0,0).
Performs no bounds-checking!

=cut

*/

void *vol3d_element(struct vol3d *v,
		unsigned int x, unsigned int y, unsigned int z)
{
	return (void *) (
		((char *)v->data)
			+ (v->elsize)*(x + y*v->nx + z*v->nx*v->ny));
}

/*

=item B<vol3d_element_wrap>

  void *vol3d_element_wrap(struct vol3d *v,
        unsigned int x, unsigned int y, unsigned int z)

Returns a pointer to the element at (x,y,z) in the given volume.
Coordinates are zero-based, so the first element is at (0,0,0).
Coordinates are wrapped-round to imnpose periodic boundary conditions,
so (-1,-1,-1) is the same as (nx-1,ny-1,nz-1).

=cut

*/

void  *vol3d_element_wrap(struct vol3d *v, int x, int y, int z)
{
        if (x<0) { x += v->nx; }
        if (y<0) { y += v->ny; }
        if (z<0) { z += v->nz; }

        if (x >= v->nx) { x -= v->nx; }
        if (y >= v->ny) { y -= v->ny; }
        if (z >= v->nz) { z -= v->nz; }

        return vol3d_element(v,x,y,z);

}

/*

=item B<vol3d_new_subvol>

  struct vol3d *vol3d_new_subvol(struct vol3d *vold, 
                  unsigned int x0, unsigned int y0, unsigned int x0,
                  unsigned int dx, unsigned int dy, unsigned int dz)

Creates a new vol3d object of size (dx,dy,dz), created from the 
subvolume of that size whose lower-coordinate corner is at
(x0,y0,z0). Returns a pointer to a new vol3d object on success, 
or C<NULL> on failure.

=cut

*/

struct vol3d *vol3d_new_subvol(struct vol3d *vold, 
	unsigned int x0, unsigned int y0, unsigned int z0,
	unsigned int dx, unsigned int dy, unsigned int dz)
{
	size_t spansize; /* Number of bytes in a span of dx */

	struct vol3d *v=NULL;
	void *p=NULL;
	int z,y;

	/* Sanity-check arguments */

	if (
		  ((x0+dx)>vold->nx)
		||((y0+dy)>vold->ny)
		||((z0+dz)>vold->nz)
	   ) {
		fprintf(stderr,"vol3d_new_subvol called with silly arguments:"
			"\nr=(%u,%u,%u) dr=(%u,%u,%u)\n"
			"vol has dimensions (%u,%u,%u)\n",
			x0,y0,z0,dx,dy,dz,vold->nx,vold->ny,vold->nz);
		return NULL;
	}

	if (NULL==(v=vol3d_new_elsize(dx,dy,dz,vold->elsize))) {
		fprintf(stderr,"vol3d_new_subvol: vol3d_new_elsize failed\n");
		return NULL;
	}
	v->datatype = vold->datatype;

	p = v->data;

	spansize = dx * vold->elsize;

	for (z=z0;z<(z0+dz);z++) {
	for (y=y0;y<(y0+dy);y++) {
		memcpy(p,vol3d_element(vold,x0,y,z),spansize);
		p = (void *)(
			((char *)p)
			+ spansize
		);
	}}

	return v;


}

/*

=item B<vol3d_new_swapzx>

  struct vol3d *vol3d_new_swapzx(struct vol3d *vold)

Copies the given volume into a new one, with the X and Z
coordinates reversed. This can be used to convert between
"Fortran order" (X varies fastest) and "C order" (Z fastest).
Returns a pointer to a new vol3d object on success, 
or C<NULL> on failure.

=cut

*/

struct vol3d *vol3d_new_swapzx(struct vol3d *vold)
{

	struct vol3d *vnew=NULL;
	void *p=NULL;
	int z,y,x;
        int newx,newy,newz;
        vol_datatype type;
        size_t size;

        newx = vold->nz;
        newy = vold->ny;
        newz = vold->nx;
        type = vold->datatype;
        size = vol3d_get_datatype_size(type);

	if (NULL==(vnew=vol3d_new_elsize(newx,newy,newz,size))) {
		fprintf(stderr,"vol3d_new: vol3d_new_elsize failed\n");
		return NULL;
	}
	vnew->datatype = type;

	p = vnew->data;

        for (z=0;z<vold->nz;z++) {
        for (y=0;y<vold->ny;y++) {
        for (x=0;x<vold->nx;x++) {
		memcpy(
                        vol3d_element(vnew,z,y,x),
                        vol3d_element(vold,x,y,z),
                        size
                      );
        }}}

	return vnew;


}

/*

=item B<vol3d_new_rotxyz:>

  struct vol3d *vol3d_new_rotxyz(struct vol3d *vold)

Copies the given volume into a new one, with X, Y, and Z
coordinates cyclically permuted to rotate about (1,1,1).
Returns a pointer to a new vol3d object on success, 
or C<NULL> on failure.

=cut

*/

struct vol3d *vol3d_new_rotxyz(struct vol3d *vold)
{

	struct vol3d *vnew=NULL;
	void *p=NULL;
	int z,y,x;
        int newx,newy,newz;
        vol_datatype type;
        size_t size;

        newx = vold->nz;
        newy = vold->ny;
        newz = vold->nx;
        type = vold->datatype;
        size = vol3d_get_datatype_size(type);

	if (NULL==(vnew=vol3d_new_elsize(newx,newy,newz,size))) {
		fprintf(stderr,"vol3d_new_rotxyz: vol3d_new_elsize failed\n");
		return NULL;
	}
	vnew->datatype = type;

	p = vnew->data;

        for (z=0;z<vold->nz;z++) {
        for (y=0;y<vold->ny;y++) {
        for (x=0;x<vold->nx;x++) {
		memcpy(
                        vol3d_element(vnew,z,x,y),
                        vol3d_element(vold,x,y,z),
                        size
                      );
        }}}

	return vnew;


}



/* This is a routine to test an algorithm out.. */
struct vol3d *vol3d_new_subvol2(struct vol3d *vold, 
	unsigned int x0, unsigned int y0, unsigned int z0,
	unsigned int dx, unsigned int dy, unsigned int dz)
{
	size_t spansize; /* Number of bytes in a span of dx */

	struct vol3d *v=NULL;
	void *p=NULL;
	int z,y;

	/* Sanity-check arguments */

	if (
		  ((x0+dx)>vold->nx)
		||((y0+dy)>vold->ny)
		||((z0+dz)>vold->nz)
	   ) {
		fprintf(stderr,"vol3d_new_subvol called with silly arguments:"
			"\nr=(%u,%u,%u) dr=(%u,%u,%u)\n"
			"vol has dimensions (%u,%u,%u)\n",
			x0,y0,z0,dx,dy,dz,vold->nx,vold->ny,vold->nz);
		return NULL;
	}

	if (NULL==(v=vol3d_new_elsize(dx,dy,dz,vold->elsize))) {
		fprintf(stderr,"vol3d_new_subvol: vol3d_new_elsize failed\n");
		return NULL;
	}
	v->datatype = vold->datatype;

	p = v->data;

	spansize = dx * vold->elsize;

	for (z=z0;z<(z0+dz);z++) {
	for (y=y0;y<(y0+dy);y++) {
		memcpy(p,vol3d_element(vold,x0,y,z),spansize);
		p = (void *)(
			((char *)p)
			+ spansize
		);
	}}

	return v;


}

/*
=item B<vol3d_set>

  void vol3d_set(struct vol3d *v, int c)

This routine sets all data in C<v> to the byte value C<c>. it is useful
for zeroing a freshly-allocated volume.

=cut
*/

void vol3d_set(struct vol3d *v, int c)
{
	memset(v->data,c,v->nx*v->ny*v->nz*v->elsize);
}

/*

=item B<vol3d_subvol_readstream>


  int vol3d_subvol_readstream(struct vol3d *v, 
  	unsigned int x0, unsigned int y0, unsigned int z0,
  	unsigned int dx, unsigned int dy, unsigned int dz,
	void (*streamfunc)(void *, void *),
	void *data
	)

This routine reads data from the given stream function into the
subvolume of the C<vol3d> object C<v> of size (dx,dy,dz) which starts at
(x0,y0,z0). The subvolume must fit entirely into C<v>.

=cut
   
*/

int vol3d_subvol_readstream(struct vol3d *v, 
	unsigned int x0, unsigned int y0, unsigned int z0,
	unsigned int dx, unsigned int dy, unsigned int dz,
	void (*streamfunc)(void *, void *),
	void *data)
{

	void *p=NULL;
	int i,z,y;

	/* Sanity-check arguments */

	if (
		  ((x0+dx)>v->nx)
		||((y0+dy)>v->ny)
		||((z0+dz)>v->nz)
	   ) {
		fprintf(stderr,"vol3d_subvol_readstream "
			"called with silly arguments:"
			"\nr=(%u,%u,%u) dr=(%u,%u,%u)\n"
			"vol has dimensions (%u,%u,%u)\n",
			x0,y0,z0,dx,dy,dz,v->nx,v->ny,v->nz);
		return -1;
	}

	p = v->data;

	for (z=z0;z<(z0+dz);z++) {
	for (y=y0;y<(y0+dy);y++) {
		p = vol3d_element(v,x0,y,z);
		for (i=0;i<dx;i++) {
			streamfunc(p,data);
			p = (char *)p + v->elsize;
		}
	}}

	return 0;


}
/*

=item B<vol3d_set_subvol>

  struct vol3d *vol3d_set_subvol(struct vol3d *v, 
                  unsigned int x0, unsigned int y0, unsigned int x0,
                  unsigned int dx, unsigned int dy, unsigned int dz,
		  void *val)

Takes a vol3d object, and sets all elements in the specified
subvolume to copies of the single supplied element C<val>.
Returns 0 on success, nonzero on failure.

=cut

*/
int vol3d_set_subvol(struct vol3d *v, 
                unsigned int x0, unsigned int y0, unsigned int z0,
                unsigned int dx, unsigned int dy, unsigned int dz,
		void *val)

{
	size_t spansize; /* Number of bytes in a span of dx */
	void *span=NULL;

	void *p=NULL;
	int i;
	int z,y;

	/* Sanity-check arguments */

	if (
		  ((x0+dx)>v->nx) /* Volume too large */
		||((y0+dy)>v->ny)
		||((z0+dz)>v->nz)
		||(0==dx*dy*dz)   /* Volume of zero size */
	   ) {
		fprintf(stderr,"vol3d_set_subvol called with silly arguments:"
			"\nr=(%u,%u,%u) dr=(%u,%u,%u)\n"
			"vol has dimensions (%u,%u,%u)\n",
			x0,y0,z0,dx,dy,dz,v->nx,v->ny,v->nz);
		return -1;
	}

	spansize = dx * v->elsize;

	if (NULL==(span=malloc(spansize))) {
		perror("vol3d_set_subvol: malloc()");
		return -1;
	}

	/* Place dx element copies into span. */
	for (i=0;i<dx;i++) {
		memcpy((char *)span+(i*v->elsize),val,v->elsize);
	}

	p = v->data;


	for (z=z0;z<(z0+dz);z++) {
	for (y=y0;y<(y0+dy);y++) {
		memcpy(p,span,spansize);
		p = (char *)p +  spansize;
	}}

	free(span);

	return 0;

}

/*

=item B<vol3d_new_upscale>

  struct vol3d *vol3d_new_upscale(struct vol3d *v, 
                  unsigned int nx, unsigned int ny, unsigned int nz)

Takes a vol3d object, of size (x,y,z), and creates a new one of size
(x*nx,y*ny,z*nz), composed of (nx,ny,nz) copies of the original volume.
nx,ny, and nz must all be greater than zero. Returns NULL on error.

=cut

*/

struct vol3d *vol3d_new_upscale(struct vol3d *v, 
                  unsigned int nx, unsigned int ny, unsigned int nz)
{

        struct vol3d *vnew=NULL;
        int x,y,z;
        int vx,vy,vz;
        size_t elsize;
        void *src,*dst;

        if ( (0==nx)||(0==ny)||(0==nx) ) {
                fprintf(stderr,"vol3d_new_upscale: bad args %u,%u,%u\n",
                                nx,ny,nz);
                return NULL;
        }

        vx=v->nx; vy=v->ny; vz=v->nz;
        elsize = vol3d_get_datatype_size(v->datatype);

        if (NULL==(vnew=vol3d_new(vx*nx,vy*ny,vz*nz,v->datatype))) {
                fprintf(stderr,"vol3d_new_upscale: vol3d_new failed.\n");
                return NULL;
        }

        /* New volume is allocated; now we just fill it in. */
        
        for (z=0;z<vz*nz;z++) {
        for (y=0;y<vy*ny;y++) {
        for (x=0;x<vx*nx;x++) {

                src = vol3d_element(v,x%vx,y%vy,z%vz);
                dst = vol3d_element(vnew,x,y,z);

                memcpy(dst,src,elsize);

        }}}

        return vnew;

}


/*

=item B<vol3d_readzslice_raw_fhandle>

  int vol3d_readzslice_raw_fhandle(struct vol3d *v, unsigned int z, FILE *f);

This routine takes an existing vol3d object, and a value of
the Z coordinate. Raw data is read from the given filehandle,
and written to the corresponding slice of the volume.

Returns 0 on success, nonzero on error.

=cut
   
*/

int vol3d_readzslice_raw_fhandle(struct vol3d *v, unsigned int z, FILE *f)
{
	struct vol3d vslice;
	if (z>=v->nz) {
		fprintf(stderr,"vol3d_readzslice_raw_fhandle: z=%u"
				"for system with nz=%u\n",z,v->nz);
		return -1;
	}


	/* Set up a dummy vol3d object, referring to the Z
	 * subslice of v
	 */

	vslice.nx = v->nx; vslice.ny = v->ny; vslice.nz = 1;
	vslice.elsize = v->elsize;
	vslice.data = vol3d_element(v,0,0,z);

	if (0!=vol3d_read_raw_fhandle(&vslice,f)) {
		fprintf(stderr,"vol3d_readzslice_raw_fhandle :"
				"vol3d_read_raw_fhandle failed.\n");
		return -1;
	}

	return 0;
}

/*

=item B<vol3d_readzslice_raw>

  int vol3d_readzslice_raw(struct vol3d *v, unsigned int z, char *fname);

This routine takes an existing vol3d object, and a value of
the Z coordinate. Raw data is read from the named file,
and written to the corresponding slice of the volume.

Returns 0 on success, nonzero on error.

=cut
   
*/

int vol3d_readzslice_raw(struct vol3d *v, unsigned int z,char *fname )
{
	FILE *f=NULL;

	if (NULL==(f=fopen(fname,"rb"))) {
		perror("vol3d_readzslice_raw: fopen");
		return -1;
	}

	if (0!=vol3d_readzslice_raw_fhandle(v,z,f)) {
		fprintf(stderr,"vol3d_readzslice_raw: "
			"vol3d_readzslice_raw_fhandle failed\n");
		fclose(f);
		return -1;
	}

	fclose(f);
	return 0;
}

/*

=item B<vol3d_new_stream_elsize>

  struct vol3d *vol3d_new_stream_elsize(
  	unsigned int nx, unsigned int ny, unsigned int nz,
  	size_t elsize,
  	void (*streamfunc)(void *, void *),
	void *data)

This function creates a new C<vol3d> object of given size, containing
elements of size C<elsize>. The function C<streamfunc> is then called
once per element, in Fortran order. The first argument to C<streamfunc>
is a pointer to the element to be written; the second argument is taken
from the C<data> argument to C<vol3d_new_stream_elsize>.  On success,
returns a pointer to a new C<vol3d> object; on failure, returns C<NULL>.

  void *mystream(void *dst, void *data)
  {
  	*(uint8_t *)dst = 0xff;
  }
 
  ...
 
  struct vol3d *v;
 
  v = vol3d_new_stream_elsize(8,8,8,sizeof(uint8_t),mystream);

=cut

*/

struct vol3d *vol3d_new_stream_elsize(
		unsigned int nx, unsigned int ny, unsigned int nz,
		size_t elsize,
		void (*streamfunc)(void *, void *),
		void *data)
{
	struct vol3d *v=NULL;
	void *p=NULL;
	int i;

	if (NULL==(v=vol3d_new_elsize(nx,ny,nz,elsize))) {
		fprintf(stderr,"vol3d_new_stream_elsize: "
				"vol3d_new_elsize failed\n");
		return NULL;
	}

	p = v->data;
	for (i=0;i<nx*ny*nz;i++) {
		streamfunc(p,data);
		p = (void *) ((char *)p + elsize);
	}
	return v;
}

/*

=item B<vol3d_new_stream>

  struct vol3d *vol3d_new_stream(
  	unsigned int nx, unsigned int ny, unsigned int nz,
  	vol_datatype datatype,
  	void (*streamfunc)(void *, void *),
  	void *data)

This function creates a new C<vol3d> object of given size, containing
elements of type C<datatype>. The function C<streamfunc> is then called
once per element, in Fortran order. The first argument to C<streamfunc>
is a pointer to the element to be written; the second argument is taken
from the C<data> argument to C<vol3d_new_stream_elsize>.  On success,
returns a pointer to a new C<vol3d> object; on failure, returns C<NULL>.

  void *mystream(void *dst, void *data)
  {
  	*(uint8_t *)dst = 0xff;
  }
 
  ...
 
  struct vol3d *v;
 
  v = vol3d_new_stream_elsize(8,8,8,VOL_TYPE_UINT8,mystream);

=cut

*/


struct vol3d *vol3d_new_stream(
		unsigned int nx, unsigned int ny, unsigned int nz,
		vol_datatype datatype,
		void (*streamfunc)(void *, void *),
		void *data)
{
	struct vol3d *v=NULL;
	size_t elsize;

	if (VOL_TYPE_UNKNOWN==datatype) {
		fprintf(stderr,"vol3d_new_stream(VOL_TYPE_UNKNOWN)\n");
		return NULL;
	} else {
		if (0==(elsize=vol3d_get_datatype_size(datatype))) {
			fprintf(stderr,"vol3d_new_stream: "
					"vol3d_get_datatype_size=0\n"
			       );
			return NULL;
		}
	}

	if (NULL==(v=vol3d_new_stream_elsize(nx,ny,nz,elsize,streamfunc,data))){
		fprintf(stderr,"vol3d_new_stream: "
				"vol3d_new_stream_elsize failed\n");
		return NULL;
	}
	v->datatype = datatype;

	return v;
}

/*

=item B<vol3d_stream>

  void vol3d_stream(struct vol3d *v,
  	void (*streamfunc)(void *, void *),
	void *data);

  void mystream(void *el, void *data)
  {
  	printf("%u\n",(unsigned int) *( (uint8_t *)el);
  }

  ...

  struct vol3d *v;
  v = vol3d_new(3,3,3,VOL_TYPE_UINT8);

  ...

  vol3d_stream(v,streamfunc);


This function iterates through all elements in the given C<vol3d>
object, in Fortran order, and calls C<streamfunc> with two arguments: a
pointer to the element, and the C<data> argument passed to
C<vol3d_stream>.

=cut

*/

void vol3d_stream(struct vol3d *v,
	void (*streamfunc)(void *, void *),
	void *data)

{
	unsigned int i,n;
	void *p;
	size_t elsize;

	n=v->nx*v->ny*v->nz;
	elsize = v->elsize;
	p = v->data;

	for (i=0;i<n;i++) {
		streamfunc(p,data);
		p = (void *) ((char *)p + elsize);
	}
	return;
}

/*
=item B<vol3d_slice_stream>

  int vol3d_slice_stream(struct vol3d *v,
        enum vol3d_direction dir,
        unsigned int index,
  	void (*streamfunc)(void *, void *),
	void *data);

  void mystream(void *el, void *data)
  {
  	printf("%u\n",(unsigned int) *( (uint8_t *)el);
  }

  ...

  struct vol3d *v;
  v = vol3d_new(3,3,3,VOL_TYPE_UINT8);

  ...

  vol3d_slice_stream(v,VOL3D_Z,2,streamfunc);


This function works like C<vol3d_stream()>, except that it only streams
through elements in a plane perpendicular to the stated direction.

Returns zero on success; nonzero indicates an error condition, in which case
there are no guarantees on the state of the volume.

=cut

*/

int vol3d_slice_stream(struct vol3d *v,
        enum vol3d_direction dir,
        unsigned int index,
	void (*streamfunc)(void *, void *),
	void *data)

{
        unsigned int imax=0;
        unsigned int x=0,y=0,z=0;

        switch(dir) {
                case VOL3D_X: imax=v->nx; break;
                case VOL3D_Y: imax=v->ny; break;
                case VOL3D_Z: imax=v->nz; break;
        }
        if (index>=imax) {
                fprintf(stderr,"vol3d_slice_stream: index out of range.\n");
                return -1;
        }

        switch(dir) {
                case VOL3D_X:
                        x=index;
                        for (z=0;z<v->nz;z++) {
                        for (y=0;y<v->ny;y++) {
                                streamfunc(vol3d_element(v,x,y,z), data);
                        }}
                        break;
                case VOL3D_Y:
                        y=index;
                        for (x=0;x<v->nx;x++) {
                        for (z=0;z<v->nz;z++) {
                                streamfunc(vol3d_element(v,x,y,z), data);
                        }}
                        break;
                case VOL3D_Z:
                        z=index;
                        for (y=0;y<v->ny;y++) {
                        for (x=0;x<v->nx;x++) {
                                streamfunc(vol3d_element(v,x,y,z), data);
                        }}
                        break;
        }

	return 0;
}

/*
=item B<vol2d_slice_2d>

  struct vol2d *vol3d_slice_2d(struct vol3d *v,
        enum vol3d_direction dir,
        unsigned int index);

  struct vol2d *slice=NULL;

  struct vol3d *v = ... ;

  slice = vol3d_slice_2d(v,VOL3D_X,3);

This function creates a C<vol2d> object from the corresponding slice
of the C<vol3d> volume. Returns a pointer to a C<struct vol3d>, or
C<NULL> on error.

=cut

*/

/* Streamfunction: first argument is a pointer to the element in
 * the 3d volume; the second argument is a pointer to a pointer to
 * the corresponding element of the 2d slice.
 */

static void vol3d_slice_2d_streamfunc(void *element, void *data)
{
        float **p=data;

        *((*p)++) = *(float *)element;

}

struct vol2d *vol3d_slice_2d(struct vol3d *v,
        enum vol3d_direction dir,
        unsigned int index)

{
        unsigned int slicex=0,slicey=0;
        unsigned int imax=0;
        struct vol2d *slice=NULL;
        float *p=NULL;

        switch(dir) {
                case VOL3D_Z: imax=v->nz; slicey=v->ny; slicex=v->nx; break;
                case VOL3D_Y: imax=v->ny; slicey=v->nx; slicex=v->nz; break;
                case VOL3D_X: imax=v->nx; slicey=v->nz; slicex=v->ny; break;
                default:
                        fprintf(stderr,"vol3d_slice_2d: Can't happen!\n");
                        return NULL;
        }
        if (index>=imax) {
                fprintf(stderr,"vol3d_slice_stream: index out of range.\n");
                return NULL;
        }

        if (NULL==(slice=vol2d_new(slicex,slicey,v->datatype))) {
                fprintf(stderr,"vol3d_slice_2d: vol2d_new failed\n");
                return NULL;
        }

        p = slice->data;

        if (0!=vol3d_slice_stream(v,dir,index,
                vol3d_slice_2d_streamfunc,&p)) {
                fprintf(stderr,"vol3d_slice_2d: vol3d_slice_stream failed.\n");
                vol2d_destroy(slice);
                return NULL;
        }

        return slice;

}

/*

=item B<vol3d_new_uc2rgba_lookup>

  struct vol3d *vol3d_new_uc2rgba_lookup(
         struct vol3d *v,
	 uint8 *r, uint8 *g, uint8 *b, uint8 *a);


  uint8 linear[256];
  int i;
  struct vol3d *vfrom,*vto;

  vfrom = vol3d_new(16,32,64,VOL_TYPE_UINT8);

  ...

  for (i=0;i<255;i++) {
  	linear[i] = i;
  }

  vto = vol3d_new_uc2rgba_lookup(vfrom,
  	linear,linear,linear,linear);


This routine will B<only> work on a vol3d object containing unsigned
8-bit char (ie byte) elements, and produces a vol3d object containing
unsigned 32-bit elements containing network-order RGBA tuples found
by looking the source byte value up in the supplied tables.

=cut


*/

struct vol3d *vol3d_new_uc2rgba_lookup(
	struct vol3d *v, uint8_t *rt, uint8_t *gt, uint8_t *bt, uint8_t *at)
{
	uint8_t *from;
	uint32_t *to;
	uint8_t data_in;
	unsigned int i;

	struct vol3d *vnew=NULL;

	if (1!=v->elsize) {
		fprintf(stderr,"vol3d_new_uc2rgba_lookup called with "
				"elsize=%u instead of 1\n",(unsigned)v->elsize);
		return NULL;
	}

	if ((NULL==rt)||(NULL==gt)||(NULL==bt)||(NULL==at)) {
		fprintf(stderr,"vol3d_new_uc2rgba_lookup called with "
				"NULL lookup table\n");
		return NULL;
	}

	if (NULL==(vnew=vol3d_new(v->nx,v->ny,v->nz,VOL_TYPE_UINT32))) {
		fprintf(stderr,"vol3d_new_uc2rgba_lookup: "
				"vol3d_new failed\n");
		return NULL;
	}

	from = (uint8_t *)v->data;
	to = (uint32_t *)vnew->data;
	for (i=0;i<v->nx*v->ny*v->nz;i++) {
		uint8_t r,g,b,a;
		unsigned long hostlong;

		data_in = *from++;

		r = rt[data_in];
		g = gt[data_in];
		b = bt[data_in];
		a = at[data_in];

		hostlong = r & 0xff;
		hostlong <<= 8;
		hostlong |= (g & 0xff);
		hostlong <<= 8;
		hostlong |= (b & 0xff);
		hostlong <<= 8;
		hostlong |= (a & 0xff);

		*to++  = (uint32_t) htonl(hostlong);
	}

	return vnew;


}

/*

=item B<vol3d_parse_dims_string>

  void vol3d_parse_dims_string(int *nx,int *ny,int *nz,char *str);

  char *string = "32,32,128";
  int nx,ny,nz;
  vol3d_parse_dims_string(&nx,&ny,&nz,string)

This routine decodes a string of the form <nx>,<ny>,<nz> into
three integers. Useful for decoding command-line arguments.

=cut

*/

void vol3d_parse_dims_string(unsigned int *nx,unsigned int *ny,unsigned int *nz,char *str)
{
        char *s=NULL,*p;
        size_t len;
        int i=0;

        len=strlen(str);
        if (NULL==(s=malloc(len+1))) {
                perror("vol3d_parse_dims_string: malloc");
                exit(-1);
        }

        strncpy(s,str,len+1);

        /* Read first digit */
        while (isdigit(s[i])) { i++; }
        if (',' != s[i]) {
                fprintf(stderr,"vol3d_parse_dims_string: Bad dimensions\n");
                exit(-1);
        }
        s[i++]=0;
        if (i>=len) {
                fprintf(stderr,"vol3d_parse_dims_string: Bad dimensions\n");
                exit(-1);
        }
        *nx = (unsigned int) atoi(s); p=s+i;

        /* Read second digit */
        while (isdigit(s[i])) { i++; }
        if (',' != s[i]) {
                fprintf(stderr,"vol3d_parse_dims_string: Bad dimensions\n");
                exit(-1);
        }
        s[i++]=0;
        if (i>=len) {
                fprintf(stderr,"vol3d_parse_dims_string: Bad dimensions\n");
                exit(-1);
        }
        *ny = (unsigned int) atoi(p); p=s+i;

        /* Read third digit */
        while (isdigit(s[i])) { i++; }
        if (0 != s[i]) {
                fprintf(stderr,"vol3d_parse_dims_string: Bad dimensions\n");
                exit(-1);
        }
        if (i!=len) {
                fprintf(stderr,"vol3d_parse_dims_string: Bad dimensions\n");
                exit(-1);
        }
        *nz = (unsigned int) atoi(p);

        free(s);
}

void vol3d_parse_floats_string(float *nx,float *ny,float *nz,char *str)
{
        char *s=NULL,*p;
        size_t len;
        int i=0;

        len=strlen(str);
        if (NULL==(s=malloc(len+1))) {
                perror("vol3d_parse_floats_string: malloc");
                exit(-1);
        }

        strncpy(s,str,len+1);

        /* Read first digit */
        while ((s[i]=='.') || isdigit(s[i])) { i++; }
        if (',' != s[i]) {
                fprintf(stderr,"vol3d_parse_floats_string: Bad input\n");
                exit(-1);
        }
        s[i++]=0;
        if (i>=len) {
                fprintf(stderr,"vol3d_parse_floats_string: Bad input\n");
                exit(-1);
        }
        *nx = (float) atof(s); p=s+i;

        /* Read second digit */
        while ((s[i]=='.') || isdigit(s[i])) { i++; }
        if (',' != s[i]) {
                fprintf(stderr,"vol3d_parse_floats_string: Bad input\n");
                exit(-1);
        }
        s[i++]=0;
        if (i>=len) {
                fprintf(stderr,"vol3d_parse_floats_string: Bad input\n");
                exit(-1);
        }
        *ny = (float) atof(p); p=s+i;

        /* Read third digit */
        while ((s[i]=='.') || isdigit(s[i])) { i++; }
        if (0 != s[i]) {
                fprintf(stderr,"vol3d_parse_floats_string: Bad input\n");
                exit(-1);
        }
        if (i!=len) {
                fprintf(stderr,"vol3d_parse_floats_string: Bad input\n");
                exit(-1);
        }
        *nz = (float) atof(p);

        free(s);
}

/*

=back

=cut

*/
