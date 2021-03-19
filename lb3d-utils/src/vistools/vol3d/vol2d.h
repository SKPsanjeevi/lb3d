#ifndef INCLUDED_VOL2D_H
#define INCLUDED_VOL2D_H

/* vol2d.c */

#define VOL_TYPE_FLOATSIZE	sizeof(float)
#define VOL_TYPE_DOUBLESIZE sizeof(double)
#define VOL2D_FLOATSIZE	VOL_TYPE_FLOATSIZE
#define VOL2D_DOUBLESIZE VOL_TYPE_DOUBLESIZE

enum vol_datatype {
	VOL_TYPE_UNKNOWN	= 0,
	VOL_TYPE_UINT8		= 1,
	VOL_TYPE_INT8		= 2,
	VOL_TYPE_UINT16		= 3,
	VOL_TYPE_INT16		= 4,
	VOL_TYPE_UINT32		= 5,
	VOL_TYPE_INT32		= 6,
	VOL_TYPE_UINT64		= 7,
	VOL_TYPE_INT64		= 8,
	VOL_TYPE_FLOAT		= 9,
	VOL_TYPE_FLOAT2		= 10,
	VOL_TYPE_FLOAT3		= 11,
	VOL_TYPE_DOUBLE		= 12,
	VOL_TYPE_DOUBLE2	= 13,
	VOL_TYPE_DOUBLE3	= 14
};

#define VOL_NTYPES	15
#define VOL2D_NTYPES   VOL_NTYPES

typedef enum vol_datatype vol_datatype;
typedef enum vol_datatype vol2d_datatype;




struct vol2d {
	unsigned int nx,ny;
	size_t elsize;
	vol_datatype datatype;
	void *data;
};


size_t vol2d_get_datatype_size(vol_datatype datatype);
struct vol2d *vol2d_new(unsigned int nx, unsigned int ny, vol_datatype datatype);
struct vol2d *vol2d_new_elsize(unsigned int nx, unsigned int ny, size_t elsize);
struct vol2d *vol2d_new_copy(struct vol2d *vold);
struct vol2d *vol2d_new_copytype(struct vol2d *vold);
struct vol2d *vol2d_new_transfer_elsize(struct vol2d *vold, size_t elsize, void (*transferfunc)(void *, void *, void *), void *data);
struct vol2d *vol2d_new_transfer(struct vol2d *vold, vol_datatype datatype, void (*transferfunc)(void *, void *, void *), void *data);
int vol2d_destroy(struct vol2d *v);
int vol2d_write_raw_fhandle(struct vol2d *v, FILE *f);
int vol2d_write_raw(struct vol2d *v, char *fname);
int vol2d_read_raw_fhandle(struct vol2d *v, FILE *f);
int vol2d_read_raw(struct vol2d *v, char *fname);
void *vol2d_element(struct vol2d *v, unsigned int x, unsigned int y);
void *vol2d_element_bc(struct vol2d *v, unsigned int x, unsigned int y);
struct vol2d *vol2d_new_subvol(struct vol2d *vold, unsigned int x0, unsigned int y0, unsigned int dx, unsigned int dy);
int vol2d_set_subvol(struct vol2d *v, unsigned int x0, unsigned int y0, unsigned int dx, unsigned int dy, void *val);
struct vol2d *vol2d_new_stream_elsize(unsigned int nx, unsigned int ny, size_t elsize, void (*streamfunc)(void *, void *), void *data);
struct vol2d *vol2d_new_stream(unsigned int nx, unsigned int ny, vol_datatype datatype, void (*streamfunc)(void *, void *), void *data);
void vol2d_stream(struct vol2d *v, void (*streamfunc)(void *, void *), void *data);
int vol2d_subvol_readstream(struct vol2d *v, 
	unsigned int x0, unsigned int y0,
	unsigned int dx, unsigned int dy,
	void (*streamfunc)(void *, void *),
	void *data);
void vol2d_set(struct vol2d *v, int c);

#endif /* INCLUDED_VOL2D_H */
