#ifndef INCLUDED_VOL3D_H
#define INCLUDED_VOL3D_H

#define VOL3D_FLOATSIZE	sizeof(float)
#define VOL3D_DOUBLESIZE sizeof(double)

#include <stdio.h>
#include <ctype.h>
#include <inttypes.h>

#include "vol2d.h"


typedef enum vol_datatype vol3d_datatype;

#define VOL3D_NTYPES	VOL_NTYPES

enum vol3d_direction {
        VOL3D_X = 0,
        VOL3D_Y = 1,
        VOL3D_Z = 2
};

struct vol3d {
	unsigned int nx,ny,nz;
	size_t elsize;
	vol_datatype datatype;
	void *data;
};


size_t vol3d_get_datatype_size(vol_datatype datatype);
struct vol3d *vol3d_new(unsigned int nx, unsigned int ny, unsigned int nz, vol_datatype datatype);
struct vol3d *vol3d_new_elsize(unsigned int nx, unsigned int ny, unsigned int nz, size_t elsize);
struct vol3d *vol3d_new_copy(struct vol3d *vold);
struct vol3d *vol3d_new_copytype(struct vol3d *vold);
struct vol3d *vol3d_new_transfer_elsize(struct vol3d *vold, size_t elsize, void (*transferfunc)(void *, void *, void *), void *data);
struct vol3d *vol3d_new_transfer(struct vol3d *vold, vol_datatype datatype, void (*transferfunc)(void *, void *, void *), void *data);
int vol3d_destroy(struct vol3d *v);
int vol3d_write_raw_fhandle(struct vol3d *v, FILE *f);
int vol3d_write_raw(struct vol3d *v, char *fname);
int vol3d_read_raw_fhandle(struct vol3d *v, FILE *f);
int vol3d_read_raw(struct vol3d *v, char *fname);
void *vol3d_element(struct vol3d *v, unsigned int x, unsigned int y, unsigned int z);
void  *vol3d_element_wrap(struct vol3d *v, int x, int y, int z);
struct vol3d *vol3d_new_subvol(struct vol3d *vold, unsigned int x0, unsigned int y0, unsigned int z0, unsigned int dx, unsigned int dy, unsigned int dz);
struct vol3d *vol3d_new_swapzx(struct vol3d *vold);
struct vol3d *vol3d_new_rotxyz(struct vol3d *vold);
int vol3d_set_subvol(struct vol3d *v, unsigned int x0, unsigned int y0, unsigned int z0, unsigned int dx, unsigned int dy, unsigned int dz, void *val);
int vol3d_readzslice_raw_fhandle(struct vol3d *v, unsigned int z, FILE *f);
int vol3d_readzslice_raw(struct vol3d *v, unsigned int z, char *fname);
struct vol3d *vol3d_new_stream_elsize(unsigned int nx, unsigned int ny, unsigned int nz, size_t elsize, void (*streamfunc)(void *, void *), void *data);
struct vol3d *vol3d_new_stream(unsigned int nx, unsigned int ny, unsigned int nz, vol_datatype datatype, void (*streamfunc)(void *, void *), void *data);
void vol3d_stream(struct vol3d *v, void (*streamfunc)(void *, void *), void *data);
struct vol3d *vol3d_new_uc2rgba_lookup(struct vol3d *v, uint8_t *rt, uint8_t *gt, uint8_t *bt, uint8_t *at);
int vol3d_subvol_readstream(struct vol3d *v, 
	unsigned int x0, unsigned int y0, unsigned int z0,
	unsigned int dx, unsigned int dy, unsigned int dz,
	void (*streamfunc)(void *, void *),
	void *data);
void vol3d_set(struct vol3d *v, int c);
void vol3d_parse_dims_string(unsigned int *nx,unsigned int *ny,unsigned int *nz,char *str);
void vol3d_parse_floats_string(float *nx,float *ny,float *nz,char *str);
int vol3d_slice_stream(struct vol3d *v,
        enum vol3d_direction dir,
        unsigned int index,
	void (*streamfunc)(void *, void *),
	void *data);
struct vol2d *vol3d_slice_2d(struct vol3d *v,
        enum vol3d_direction dir,
        unsigned int index);
struct vol3d *vol3d_new_upscale(struct vol3d *v, 
                  unsigned int nx, unsigned int ny, unsigned int nz);

#endif /* INCLUDED_VOL3D_H */
