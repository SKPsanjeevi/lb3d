#ifndef INCLUDED_VOL3DF_H
#define INCLUDED_VOL3DF_H

#include "vol3d.h"
#include "vector.h"
#include "config.h"

struct vol3df_stats {
        float max;
        float min;
        float mean;
        float variance;
};

struct vol3df_dist {
        unsigned int nbins;
        unsigned int *count;
        float dphi;
        struct vol3df_stats *stats;
};

enum vol3df_filetype {
        VOL3DF_FTYPE_UNKNOWN = -1,
        VOL3DF_FTYPE_RAW = 0,
        VOL3DF_FTYPE_XDR = 1,
        VOL3DF_FTYPE_HDF5= 2,
        VOL3DF_FTYPE_VTK = 3,
        VOL3DF_FTYPE_UCVTK = 4,
        VOL3DF_FTYPE_UC = 5,
        VOL3DF_FTYPE_ISO = 6,
        VOL3DF_FTYPE_UCNRRD = 7,
        VOL3DF_FTYPE_FNRRD = 8
};


struct vol3df_file_hints {
        enum vol3df_filetype type;
        char *filename;
        char *h5path;
        unsigned int nx,ny,nz;
};

extern const struct vol3df_file_hints vol3df_default_file_hints;

void vol3df_shape_operator(struct vol3d *v, float *S,
        float rx, float ry, float rz);
float vol3df_mean_curvature(struct vol3d *v, float rx, float ry, float rz);
float vol3df_gaussian_curvature(struct vol3d *v,float rx, float ry, float rz);
int vol3df_curvatures(struct vol3d *v,float rx, float ry, float rz,
        float *H,float *K) ;
vector vol3df_gradientvec(struct vol3d *v,int x, int y, int z);
vector vol3df_trilinear_gradient_vec(struct vol3d *v, float rx, float ry, float rz);
float vol3df_trilinear_point(struct vol3d *v, float x, float y, float z);
void vol3df_trilinear_gradient(struct vol3d *v,float rx, float ry, float rz,
                float *gxp, float *gyp, float *gzp);
void vol3df_gradient(struct vol3d *v,int x, int y, int z,
                float *dxp, float *dyp, float *dzp);
struct vol3d *vol3df_new(unsigned int nx, unsigned int ny, unsigned int nz);
struct vol3d *vol3df_new_zero(unsigned int nx,unsigned int ny, unsigned int nz);
void vol3df_xdr_streamfunc(void *el, void *data);
struct vol3d *vol3df_new_from_xdrfh(unsigned int nx, unsigned int ny, unsigned int nz, FILE *f);
int vol3df_write_xdrfh(struct vol3d *v, FILE *f);
struct vol3d *vol3df_new_from_rawfh(unsigned int nx, unsigned int ny, unsigned int nz, FILE *f);
struct vol3d *vol3df_new_from_isofh (FILE *f) ;
int vol3df_write_fh(struct vol3d *v, FILE *f);
int vol3df_write_isofh(struct vol3d *v, FILE *f);
struct vol3d *vol3df_new_fromhdf5(char *fname, char *path);
struct vol3d *vol3df_new_from_file_heuristic(
                char *filename,
                struct vol3df_file_hints *hints
                );
struct vol3d *vol3df_new_from_fh_heuristic(
                FILE *f,
                struct vol3df_file_hints *hints
                );

int vol3df_write_hdf5(struct vol3d *v, char *fname, char *path);
int vol3df_write_vtkfh(struct vol3d *v, FILE *f);
struct vol3df_stats *vol3df_getstats(struct vol3d *v);
struct vol3d *vol3df_new_erosion(struct vol3d *v);
struct vol3d *vol3df_new_dilation(struct vol3d *v);
float vol3df_total_vol_euler(struct vol3d *v,
                int (*testfunc)(float, void * ),
                void *data);
float vol3df_total_surf_euler(struct vol3d *v,
                int (*testfunc)(float, void * ),
                void *data);
int vol3df_gt(struct vol3d *v, float thresh);
int vol3df_lt(struct vol3d *v, float thresh);
int vol3df_range(struct vol3d *v, float pmax, float pmin);
struct vol3d *vol3df_new_nncount(struct vol3d *v,
                int (*testfunc)(float, void * ),
                void *data);
int vol3df_write_ucvtkfh(struct vol3d *v, FILE *f);
int vol3df_write_ucnrrdfh(struct vol3d *v, FILE *f);
int vol3df_write_fnrrdfh(struct vol3d *v, FILE *f);
int vol3df_write_ucfh (struct vol3d *v,FILE *f);
struct vol3df_dist *vol3df_getdist(struct vol3d *v, int nbins);
int vol3df_abs(struct vol3d *v);
struct vol3d *vol3df_new_smear(struct vol3d *v, int delta);
float vol3df_wrap(struct vol3d *v, int x, int y, int z);
float *vol3df_wrap_ptr(struct vol3d *v, int x, int y, int z);
void vol3df_fill_func(struct vol3d *v,
                float (*func)(int , int , int , void *), void *data);
struct vol3d *vol3df_new_pointwise_binaryop( struct vol3d *a, struct vol3d *b,
                float (*func)(float, float));
struct vol3d *vol3df_new_vol_euler(struct vol3d *v,
                int (*testfunc)(float, void * ),
                void *data);
struct vol3d *vol3df_new_surf_euler(struct vol3d *v,
                int (*testfunc)(float, void * ),
                void *data);
struct vol3d *vol3df_new_roll(struct vol3d *v, int dx, int dy, int dz);
struct vol3d *vol3df_new_vol_euler8(struct vol3d *v, float threshold);
struct vol3d *vol3df_new_surf_euler8(struct vol3d *v, float threshold);
enum vol3df_filetype vol3df_parse_ftype(char *s);

int vol3df_whitenoise(struct vol3d *v,float magnitude,float zeropoint);

#ifdef HAVE_FFTW
struct vol3d *vol3df_new_rfft(struct vol3d *v);
struct vol3d *vol3df_new_spectrum(struct vol3d *v);
struct binnery *vol3df_sfactor(struct vol3d *v);
struct vol3d *vol3df_autocorrelation(struct vol3d *v);
int vol3df_inverse_f_noise(struct vol3d *v);
#endif /* HAVE_FFTW */

int vol3df_logabs(struct vol3d *v);
int vol3df_zerobounds(struct vol3d *v);
int vol3df_rangefilter(struct vol3d *v, float pmax, float pmin);
enum vol3df_filetype vol3df_guess_fname_ftype(char *s);
#endif /* INCLUDED_VOL3DF_H */
struct vol3d *vol3df_new_downsample(struct vol3d *v, int subx, int suby,int subz);
#ifdef HAVE_PNG
int vol3df_write_pngs(struct vol3d *v, enum vol3d_direction dir,char **filenames );
int vol3df_write_png(struct vol3d *v, enum vol3d_direction dir,
                unsigned int index, char *filename, float phimax, float phimin);
int vol3df_write_png_normalize(struct vol3d *v, enum vol3d_direction dir,
                unsigned int index, char *fname);
#endif /* HAVE_PNG */
struct vol3d *vol3df_new_scale_interpolated(struct vol3d *v, int newx, int newy, int newz);
struct vol3d *vol3df_new_transformed(struct vol3d *v,float *mat);
struct vol3d *vol3df_new_H(struct vol3d *v);
int vol3df_normalise(struct vol3d *v,float amplitude);
matrix vol3df_trilinear_hessian_matrix(struct vol3d *v,
        float x, float y, float z);
void vol3df_hessian(struct vol3d *v, int x, int y, int z,
        float *hxx,float *hxy, float *hxz, float *hyy, float *hyz, float *hzz);
void vol3df_trilinear_hessian(struct vol3d *v, float x, float y, float z,
 float *hxxp,float *hxyp, float *hxzp, float *hyyp, float *hyzp, float *hzzp);
matrix vol3df_geometry_tensor(struct vol3d *v, float x, float y, float z);

float *vol3df_ortholine(struct vol3d *v, enum vol3d_direction dir, unsigned int coord0, unsigned int coord1);

