#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>

#include <unistd.h>

#ifdef HAVE_HDF5
#include <hdf5.h>
#endif

#include "vol3df.h"
#include "binnery.h"

#include "isosurface.h"






int ltfilter(float x,void *d) { 
        if (x<*(float *)d) { return 1.0; } else { return 0.0; }

}

int gtfilter(float x,void *d) { 
        if (x>*(float *)d) { return 1.0; } else { return 0.0; }
}


/* write the output. Special-case HDF5. */

int write_output(struct vol3d *v,char *outfilename,
                enum vol3df_filetype otype, char *h5path)
{

        FILE *outfile=NULL;

        if (VOL3DF_FTYPE_UNKNOWN == otype) { 
                /* Try to guess from the filename extension */
                if (VOL3DF_FTYPE_UNKNOWN==(otype=vol3df_guess_fname_ftype(outfilename))) {
                        fprintf(stderr,
                               "Please specify output file type with -g\n");
                        return -1;
                }
        }

        if (VOL3DF_FTYPE_HDF5 == otype) {
                if (0==strcmp("-",outfilename)) {
                       fprintf(stderr,"Filter operation disallowed for HDF5\n");
                       return -1;
                }
                if (NULL==h5path) {
                        /* Use default */
                        h5path = "/OutArray";
                }

                if (0!=vol3df_write_hdf5(v,outfilename,h5path)) {
                        fprintf(stderr,"vol3df_write_hdf5 failed\n");
                        return -1;
                }

                return 0;
        }

        /* Output filename of "-" -> stdout */

        if (0==strcmp("-",outfilename)) {
                outfile=stdout;
        } else {
                if (NULL==(outfile=fopen(outfilename,"wb"))) {
                        perror("fopen()");
                        return -1;
                }
        }

        /* outfile is now a valid file descriptor */



        switch(otype) {
                case VOL3DF_FTYPE_RAW:
                        if (0!=vol3df_write_fh(v,outfile)) {
                                fprintf(stderr,"vol3df_write_fh failed\n");
                        }
                        break;
                case VOL3DF_FTYPE_XDR:
                        if (0!=vol3df_write_xdrfh(v,outfile)) {
                                fprintf(stderr,"vol3df_write_xdrfh failed\n");
                        }
                        break;
                case VOL3DF_FTYPE_VTK:
                        if (0!=vol3df_write_vtkfh(v,outfile)) {
                                fprintf(stderr,"vol3df_write_vtkfh failed\n");
                        }
                        break;
                case VOL3DF_FTYPE_UCVTK:
                        if (0!=vol3df_write_ucvtkfh(v,outfile)) {
                                fprintf(stderr,"vol3df_write_ucvtkfh failed\n");
                        }
                        break;
                case VOL3DF_FTYPE_UC:
                        if (0!=vol3df_write_ucfh(v,outfile)) {
                                fprintf(stderr,"vol3df_write_ucfh failed\n");
                        }
                        break;
                case VOL3DF_FTYPE_ISO:
                        if (0!=vol3df_write_isofh(v,outfile)) {
                                fprintf(stderr,"vol3df_write_isofh failed\n");
                        }
                        break;
                case VOL3DF_FTYPE_UCNRRD:
                        if (0!=vol3df_write_ucnrrdfh(v,outfile)) {
                               fprintf(stderr,"vol3df_write_ucnrrdfh failed\n");
                        }
                        break;
                case VOL3DF_FTYPE_FNRRD:
                        if (0!=vol3df_write_fnrrdfh(v,outfile)) {
                               fprintf(stderr,"vol3df_write_fnrrdfh failed\n");
                        }
                        break;
                default:
                        fprintf(stderr,"Internal error: unknown output type\n");
                        return -1;
        }
        return 0;
}


int null_wrapper(int argc, char *argv[], int optind, struct vol3d **vp)
{
        return 0;
}

/* Today we are going to abuse the preprocessor in the name of laziness. */

#define wrap_simplefilter(name,func) \
int name ## _wrapper(int argc, char *argv[], int optind, struct vol3d **vp) \
{ \
        struct vol3d *vnew=NULL; \
        if (NULL==(vnew= func (*vp))) { return -1; } \
        vol3d_destroy(*vp); *vp = vnew; return 0; \
} \

wrap_simplefilter(erode,vol3df_new_erosion)
wrap_simplefilter(dilate,vol3df_new_dilation)
wrap_simplefilter(swapzx,vol3d_new_swapzx)
wrap_simplefilter(rotxyz,vol3d_new_rotxyz)
wrap_simplefilter(H,vol3df_new_H)

int areadensity_wrapper(int argc, char *argv[], int optind, struct vol3d **vp)
{
        float isoval=0.0;
        struct vol3d *vnew=NULL;

        if (optind!=argc) { isoval=(float)atof(argv[optind++]);}

        if (NULL==(vnew=vol3df_new_areadensity(*vp,isoval))){return -1; }
        vol3d_destroy(*vp); *vp=vnew; return 0;
}

int h2density_wrapper(int argc, char *argv[], int optind, struct vol3d **vp)
{
        float isoval=0.0;
        struct vol3d *vnew=NULL;

        if (optind!=argc) { isoval=(float)atof(argv[optind++]);}

        if (NULL==(vnew=vol3df_new_h2density(*vp,isoval))){return -1; }
        vol3d_destroy(*vp); *vp=vnew; return 0;
}

int nngt_wrapper(int argc, char *argv[], int optind, struct vol3d **vp)
{
        float thresh=0.0;
        struct vol3d *vnew=NULL;

        if (optind!=argc) { thresh=(float)atof(argv[optind++]);}

        if (NULL==(vnew=vol3df_new_nncount(*vp,gtfilter,&thresh))){return -1; }
        vol3d_destroy(*vp); *vp=vnew; return 0;
}

int nnlt_wrapper(int argc, char *argv[], int optind, struct vol3d **vp)
{
        float thresh=0.0;
        struct vol3d *vnew=NULL;

        if (optind!=argc) { thresh=(float)atof(argv[optind++]);}

        if (NULL==(vnew=vol3df_new_nncount(*vp,ltfilter,&thresh))){return -1; }
        vol3d_destroy(*vp); *vp=vnew; return 0;
}

int open_wrapper(int argc, char *argv[], int optind, struct vol3d **vp) {
        struct vol3d *vnew=NULL;

        if (NULL==(vnew=vol3df_new_erosion(*vp))) { return -1; }
        vol3d_destroy(*vp); *vp = vnew ; vnew = NULL;

        if (NULL==(vnew=vol3df_new_dilation(*vp))) { return -1; }
        vol3d_destroy(*vp); *vp = vnew;
        return 0;
}

int close_wrapper(int argc, char *argv[], int optind, struct vol3d **vp) {
        struct vol3d *vnew=NULL;

        if (NULL==(vnew=vol3df_new_dilation(*vp))) { return -1; }
        vol3d_destroy(*vp); *vp = vnew;

        if (NULL==(vnew=vol3df_new_erosion(*vp))) { return -1; }
        vol3d_destroy(*vp); *vp = vnew ; vnew = NULL;

        return 0;
}

int subvol_wrapper(int argc, char *argv[], int optind, struct vol3d **vp) {
        struct vol3d *vnew=NULL;
        struct vol3d *v=NULL;
        int x0=0,y0=0,z0=0;
        int dx=0,dy=0,dz=0;
        int nargs;

        nargs = argc-optind;

        v=*vp;


        switch(nargs) {
                case 6:
                        x0=atoi(argv[optind+3]);
                        y0=atoi(argv[optind+4]);
                        z0=atoi(argv[optind+5]);
                        /* Fall through */
                case 3:
                        dx=atoi(argv[optind]);
                        dy=atoi(argv[optind+1]);
                        dz=atoi(argv[optind+2]);
                        break;
                default:
                        fprintf(stderr,
                                "Bad number of arguments\n");
                        return -1;
        }

        if (NULL==(vnew=
                    vol3d_new_subvol(v,
                            x0,y0,z0,
                            dx,dy,dz))) {
                fprintf(stderr,"vol3d_new_subvol failed\n");
                return -1;
        }
        vol3d_destroy(v);
        *vp = vnew;
        return 0;
}

int scale_wrapper(int argc, char *argv[], int optind, struct vol3d **vp) {
        struct vol3d *vnew=NULL;
        struct vol3d *v=NULL;
        unsigned int nx,ny,nz;
        int nargs;

        nargs = argc-optind;

        v=*vp;


        if (nargs!=3) {
                        fprintf(stderr,
                                "Bad number of arguments\n");
                        return -1;
        }

        nx=atoi(argv[optind]);
        ny=atoi(argv[optind+1]);
        nz=atoi(argv[optind+2]);

        if (NULL==(vnew=vol3df_new_scale_interpolated(v,
                            nx,ny,nz))) {
                fprintf(stderr,"vol3d_new_scale_interpolated failed\n");
                return -1;
        }
        vol3d_destroy(v);
        *vp = vnew;
        return 0;
}

int whitenoise_wrapper(int argc, char *argv[], int optind, struct vol3d **vp) {
    int nargs;
    float magnitude=1.0;
    float zeropoint = 0.0;

    nargs = argc-optind;

    if (optind!=argc) { magnitude = atof(argv[optind++]); };
    if (optind!=argc) { zeropoint = atof(argv[optind++]); };
    if (optind!=argc) { fprintf(stderr,"Too many arguments for whitenoise\n");};

    vol3df_whitenoise(*vp,magnitude,zeropoint);
    return 0;
}

#ifdef HAVE_FFTW
int ifnoise_wrapper(int argc, char *argv[], int optind, struct vol3d **vp) {

    vol3df_inverse_f_noise(*vp);
    return 0;
}
#endif /* HAVE_FFTW */

int rot111_wrapper(int argc, char *argv[], int optind, struct vol3d **vp) {
        struct vol3d *vnew=NULL;
        struct vol3d *v=NULL;
        float mat[9];

        mat[0] = sqrt(2.0/3.0);
        mat[1] = 1.0/sqrt(3.0);
        mat[2] = 0;
        mat[3] = -1.0/sqrt(6.0);
        mat[4] = 1/sqrt(3.0);
        mat[5] = 1.0/sqrt(2.0);
        mat[6] = 1.0/sqrt(6.0);
        mat[7] = -1.0/sqrt(3.0);
        mat[8] = 1.0/sqrt(2.0);

        v=*vp;

        if (NULL==(vnew=vol3df_new_transformed(*vp,mat))) {
                fprintf(stderr,"vol3d_new_transformed failed\n");
                return -1;
        }
        vol3d_destroy(v);
        *vp = vnew;
        return 0;
}

int shear_wrapper(int argc, char *argv[], int optind, struct vol3d **vp) {
        struct vol3d *vnew=NULL;
        struct vol3d *v=NULL;
        float mat[9] = { 
                        1,   0, 0,
                        0,   1, 0,
                        0,   0, 1 };
        int nargs;

        nargs = argc-optind;

        if (nargs!=1) {
                        fprintf(stderr,
                                "shear takes one argument!\n");
                        return -1;
        }

        mat[1] = atof(argv[optind]);

        v=*vp;

        if (NULL==(vnew=vol3df_new_transformed(*vp,mat))) {
                fprintf(stderr,"vol3d_new_transformed failed\n");
                return -1;
        }
        vol3d_destroy(v);
        *vp = vnew;
        return 0;
}

int upscale_wrapper(int argc, char *argv[], int optind, struct vol3d **vp) {
        struct vol3d *vnew=NULL;
        struct vol3d *v=NULL;
        unsigned int nx,ny,nz;
        int nargs;

        nargs = argc-optind;

        v=*vp;


        if (nargs!=3) {
                        fprintf(stderr,
                                "Bad number of arguments\n");
                        return -1;
        }

        nx=atoi(argv[optind]);
        ny=atoi(argv[optind+1]);
        nz=atoi(argv[optind+2]);

        if (NULL==(vnew=vol3d_new_upscale(v,
                            nx,ny,nz))) {
                fprintf(stderr,"vol3d_new_upscale failed\n");
                return -1;
        }
        vol3d_destroy(v);
        *vp = vnew;
        return 0;
}

int gt_wrapper(int argc, char *argv[], int optind, struct vol3d **vp) {
        float thresh=0.0;
        if (optind!=argc) { thresh=(float)atof(argv[optind++]);}
        vol3df_gt(*vp,thresh); return 0;
}

int lt_wrapper(int argc, char *argv[], int optind, struct vol3d **vp) {
        float thresh=0.0;
        if (optind!=argc) { thresh=(float)atof(argv[optind++]);}
        vol3df_lt(*vp,thresh); return 0;
}

int range_wrapper(int argc, char *argv[], int optind, struct vol3d **vp) {
        float pmax,pmin;
        if (argc != optind + 2) { 
                fprintf(stderr,"Bad number of arguments\n");
                return 0;
        }

        pmax = (float)atof(argv[optind++]);
        pmin = (float)atof(argv[optind++]);

        vol3df_range(*vp,pmax,pmin); return 0;
}

int rangefilter_wrapper(int argc, char *argv[], int optind, struct vol3d **vp) {
        float pmax,pmin;
        if (argc != optind + 2) { 
                fprintf(stderr,"Bad number of arguments\n");
                return 0;
        }

        pmax = (float)atof(argv[optind++]);
        pmin = (float)atof(argv[optind++]);

        vol3df_rangefilter(*vp,pmax,pmin); return 0;
}
        
int abs_wrapper(int argc, char *argv[], int optind, struct vol3d **vp) {
        vol3df_abs(*vp); return 0;
}

int stats_wrapper(int argc, char *argv[], int optind, struct vol3d **vp) {
        struct vol3df_stats *stats=NULL;
        struct vol3d *v=NULL;

        stats=vol3df_getstats(v=*vp);
        printf("Total %f max %f min %f mean %f sd %f\n",
                        stats->mean*v->nx*v->ny*v->nz,
                        stats->max,
                        stats->min,
                        stats->mean,
                        sqrt(stats->variance));
        return 0;
}

int dist_wrapper(int argc, char *argv[], int optind, struct vol3d **vp) {

        struct vol3df_dist *dist=NULL;
        int i;
        int nbins=32;

        if (optind!=argc) { nbins=atoi(argv[optind++]);}

        if (NULL==(dist=vol3df_getdist(*vp,nbins))) { return -1; }

        for (i=0;i<dist->nbins;i++) {
              printf("%f %d\n",
                 dist->stats->min+dist->dphi*(float)i,
                 dist->count[i]
              );
        }
        return 0;
}

int smear_wrapper(int argc, char *argv[], int optind, struct vol3d **vp) {
        int delta=2; /* Maximum displacement */
        struct vol3d *vnew=NULL;

        if (optind!=argc) { delta=atoi(argv[optind++]);}
        if (NULL==(vnew=vol3df_new_smear(*vp,delta))){return -1; }
        vol3d_destroy(*vp); *vp=vnew; return 0;
}


int surf_euler_wrapper(int argc, char *argv[], int optind, struct vol3d **vp) {
        float threshold=0.0;
        struct vol3d *vnew=NULL;

        if (optind!=argc) { threshold=atof(argv[optind++]);}
        if (NULL==(vnew=vol3df_new_surf_euler(*vp,gtfilter,&threshold)))
                {return -1; }
        vol3d_destroy(*vp); *vp=vnew; return 0;
}

int vol_euler_wrapper(int argc, char *argv[], int optind, struct vol3d **vp) {
        float threshold=0.0;
        struct vol3d *vnew=NULL;

        if (optind!=argc) { threshold=atof(argv[optind++]);}
        if (NULL==(vnew=vol3df_new_vol_euler(*vp,gtfilter,&threshold)))
                {return -1; }
        vol3d_destroy(*vp); *vp=vnew; return 0;
}

int total_vol_euler_wrapper(int argc, char *argv[], int optind, struct vol3d **vp) {

        float threshold=0.0;
        float chi=0.0;

        if (optind!=argc) { threshold=atof(argv[optind++]);}

        chi = vol3df_total_vol_euler(*vp,gtfilter,&threshold);
        printf("%.10f\n",chi);

        return 0;
}

int total_surf_euler_wrapper(int argc, char *argv[], int optind, struct vol3d **vp) {

        float threshold=0.0;
        float chi=0.0;

        if (optind!=argc) { threshold=atof(argv[optind++]);}

        chi = vol3df_total_surf_euler(*vp,gtfilter,&threshold);
        printf("%.10f\n",chi);

        return 0;
}

int surf_euler8_wrapper(int argc, char *argv[], int optind, struct vol3d **vp)
{
        float thresh=0.0;
        struct vol3d *vnew=NULL;

        if (optind!=argc) { thresh=(float)atof(argv[optind++]);}

        if (NULL==(vnew=vol3df_new_surf_euler8(*vp,thresh))){return -1; }
        vol3d_destroy(*vp); *vp=vnew; return 0;
}

int vol_euler8_wrapper(int argc, char *argv[], int optind, struct vol3d **vp)
{
        float thresh=0.0;
        struct vol3d *vnew=NULL;

        if (optind!=argc) { thresh=(float)atof(argv[optind++]);}

        if (NULL==(vnew=vol3df_new_vol_euler8(*vp,thresh))){return -1; }
        vol3d_destroy(*vp); *vp=vnew; return 0;
}

/* Fill given volume with constant value (zero by default) */

int fill_wrapper(int argc, char *argv[], int optind, struct vol3d **vp) {
        int i,n;
        struct vol3d *v=NULL;
        float fillval=0.0;

        v=*vp;
        n=v->nx*v->ny*v->nz;

        if (optind!=argc) { fillval=atof(argv[optind++]);}
        for (i=0;i<n;i++) {
                *((float *)v->data + i) = fillval;
        }
        return 0;
}

int setpoint_wrapper(int argc, char *argv[], int optind, struct vol3d **vp) {
        struct vol3d *v=NULL;
        int nargs;
        int x,y,z;
        float phi;

        v=*vp;
        nargs = argc-optind;

        if (4!=nargs)
            { fprintf(stderr,"Bad number of arguments\n"); return -1; }

        x=atoi(argv[optind++]);
        y=atoi(argv[optind++]);
        z=atoi(argv[optind++]);

        if ((x>=v->nx)||(y>v->ny)||(z>v->nz)||(x<0)||(y<0)||(z<0)) {
            fprintf(stderr,"Point %d %d %d out of range.\n",x,y,z);
            return -1;
        }

        phi = atof(argv[optind++]);

        *(float *)vol3d_element(v,x,y,z) = phi;

        return 0;
}

#ifdef HAVE_FFTW

int rfft_wrapper(int argc, char *argv[], int optind, struct vol3d **vp) {
        struct vol3d *vnew=NULL;

        if (NULL==(vnew=vol3df_new_rfft(*vp))){return -1; }
        vol3d_destroy(*vp); *vp=vnew; return 0;
}

int spectrum_wrapper(int argc, char *argv[], int optind, struct vol3d **vp) {
        struct vol3d *vnew=NULL;

        if (NULL==(vnew=vol3df_new_spectrum(*vp))){return -1; }
        vol3d_destroy(*vp); *vp=vnew; return 0;
}

int autocorrelation_wrapper(int argc, char *argv[], int optind, struct vol3d **vp) {

        if (NULL==(vol3df_autocorrelation(*vp))) { return -1; }
        fprintf(stderr,"done wrapper\n");fflush(stderr);
        return 0;
}
#endif /* HAVE_FFTW */

int logabs_wrapper(int argc, char *argv[], int optind, struct vol3d **vp) {
        return vol3df_logabs(*vp);
}

int zerobounds_wrapper(int argc, char *argv[], int optind, struct vol3d **vp) {
        return vol3df_zerobounds(*vp);
}


struct wavevector { float kx,ky,kz; };

static float gyroid_func(int x, int y, int z, void *data)
{
        struct wavevector *gi;

        gi=data;
        return cos(x*gi->kx)*sin(y*gi->ky)
             + cos(y*gi->ky)*sin(z*gi->kz)
             + cos(z*gi->kz)*sin(x*gi->kx);
}

static float sine_func(int x, int y, int z, void *data)
{
        struct wavevector *k;

        k=data;
        return sin(x*k->kx)*sin(y*k->ky)*sin(z*k->kz);
}

static float psurface_func(int x, int y, int z, void *data)
{
        struct wavevector *gi;
        float xx,yy,zz;

        gi=data;
        xx = x*gi->kx; yy = y*gi->ky; zz = z*gi->kz; 

        return cos(xx) + cos(yy) + cos(zz);
}

static float dsurface_func(int x, int y, int z, void *data)
{
        struct wavevector *gi;
        float xx,yy,zz;

        gi=data;
        xx = x*gi->kx; yy = y*gi->ky; zz = z*gi->kz; 

        return    sin(xx)*sin(yy)*sin(zz)
                + sin(xx)*cos(yy)*cos(zz)
                + cos(xx)*sin(yy)*cos(zz)
                + cos(xx)*cos(yy)*sin(zz);
}

static float gaussian_func(int x, int y, int z, void *data)
{
        struct wavevector *gi;
        float xx,yy,zz;

        gi=data;
        xx = x*gi->kx; yy = y*gi->ky; zz = z*gi->kz; 

        return exp(-(xx*xx+yy*yy+zz*zz));

}


int wavefunc_wrapper(int argc, char *argv[], int optind, struct vol3d **vp,
                float (*func)(int , int , int , void *)) {
        struct vol3d *v=NULL;
        float gnx=1,gny=1,gnz=1;
        struct wavevector k;

        v=*vp;

        if (optind!=argc) { gnx=gny=gnz=atof(argv[optind++]); }
        if (optind!=argc) { gny=gnz=atof(argv[optind++]); }
        if (optind!=argc) { gnz=atof(argv[optind++]); }

        k.kx = 2.0*M_PI/((float)v->nx)*gnx;
        k.ky = 2.0*M_PI/((float)v->ny)*gny;
        k.kz = 2.0*M_PI/((float)v->nz)*gnz;

        vol3df_fill_func(v,func,&k);

        return 0;
}

int roll_wrapper(int argc, char *argv[], int optind, struct vol3d **vp) {

        struct vol3d *vnew=NULL;
        int dx=1,dy=1,dz=1;

        if (optind!=argc) { dx=atoi(argv[optind++]); }
        if (optind!=argc) { dy=atoi(argv[optind++]); }
        if (optind!=argc) { dz=atoi(argv[optind++]); }

        if (NULL==(vnew=vol3df_new_roll(*vp,dx,dy,dz))){return -1; }
        vol3d_destroy(*vp); *vp=vnew;

        return 0;
}


int gyroid_wrapper(int argc, char *argv[], int optind, struct vol3d **vp) {

        return wavefunc_wrapper(argc,argv,optind,vp,gyroid_func);
}

int sine_wrapper(int argc, char *argv[], int optind, struct vol3d **vp) {

        return wavefunc_wrapper(argc,argv,optind,vp,sine_func);
}

int psurface_wrapper(int argc, char *argv[], int optind, struct vol3d **vp) {

        return wavefunc_wrapper(argc,argv,optind,vp,psurface_func);
}

int dsurface_wrapper(int argc, char *argv[], int optind, struct vol3d **vp) {

        return wavefunc_wrapper(argc,argv,optind,vp,dsurface_func);
}

int gaussian_wrapper(int argc, char *argv[], int optind, struct vol3d **vp) {

        return wavefunc_wrapper(argc,argv,optind,vp,gaussian_func);
}

int func_wrapper(int argc, char *argv[], int optind, struct vol3d **vp,
                float (*func)(int , int , int , void *)) {
        struct vol3d *v=NULL;

        v=*vp;

        vol3df_fill_func(v,func,NULL);

        return 0;
}

static float odist_func(int x, int y, int z, void *data)
{
        float dx,dy,dz,*o;

        o = (float *)data;
        dx = x-o[0]; dy = y-o[1]; dz = z-o[2];

        return sqrt(dx*dx+dy*dy+dz*dz);
}

static float ogauss_func(int x, int y, int z, void *data)
{
        float dx,dy,dz,*o;

        o = (float *)data;
        dx = x-o[0]; dy = y-o[1]; dz = z-o[2];

        return exp(-(dx*dx+dy*dy+dz*dz));
}

int odist_wrapper(int argc, char *argv[], int optind, struct vol3d **vp) {
        struct vol3d *v=NULL;
        float o[3];

        v=*vp;

        o[0]=0.5*(float)v->nx; o[1]=0.5*(float)v->ny; o[2]=0.5*(float)v->nz;

        vol3df_fill_func(v,odist_func,o);

        return 0;
}

int ogauss_wrapper(int argc, char *argv[], int optind, struct vol3d **vp) {
        struct vol3d *v=NULL;
        float o[3];

        v=*vp;

        o[0]=0.5*(float)v->nx; o[1]=0.5*(float)v->ny; o[2]=0.5*(float)v->nz;

        vol3df_fill_func(v,ogauss_func,o);

        return 0;
}

static float cylrad_func(int x, int y, int z, void *data)
{
        float dx,dy;
        struct vol3d *v;

        v=(struct vol3d *)data;

        dx = y-0.5*v->ny;
        dy = z-0.5*v->nz;

        return sqrt(dx*dx+dy*dy);
}

static float xjump_func(int x, int y, int z, void *data)
{
        float threshold;

        threshold = *(float *)data;

        if (x<=threshold) 
        { return 1; } else { return 0; }

}

int xjump_wrapper(int argc, char *argv[], int optind, struct vol3d **vp) {
        struct vol3d *v=NULL;
        float threshold;

        v=*vp;
        threshold=0.5*v->nx;

        vol3df_fill_func(v,xjump_func,&threshold);

        return 0;
}

static float lamella_func(int x, int y, int z, void *data)
{
        int spacing,width;
        
        spacing=*(int *)data;
        width = *((int *)data + 1);

        x = x%spacing;
        if (x<width) { return 1.0; } else { return 0.0; }
}

int lamella_wrapper(int argc, char *argv[], int optind, struct vol3d **vp) {
        struct vol3d *v=NULL;
        int data[2];

        v=*vp;

        if ( (argc-optind) != 2 ) {
                fprintf(stderr,"lamella needs two arguments: spacing, width\n");
                return -1;
        }
        data[0]=atoi(argv[optind++]); /* spacing */
        data[1]=atoi(argv[optind++]); /* width */

        vol3df_fill_func(v,lamella_func,data);

        return 0;
}


int cylrad_wrapper(int argc, char *argv[], int optind, struct vol3d **vp) {


        vol3df_fill_func(*vp,cylrad_func,*vp);

        return 0;
}

int binaryop_wrapper(int argc, char *argv[], int optind, struct vol3d **vp,
                float (*func)(float,float)
                ) {

        char *bfilename=NULL;
        struct vol3d *a=NULL,*b=NULL,*c=NULL;
        struct vol3df_file_hints hints = vol3df_default_file_hints;

        a = *vp;

        if (optind!=argc) { bfilename=argv[optind++]; }
        else { fprintf(stderr,"Need name of file to subtract\n"); return -1;}

        hints.nx = a->nx;
        hints.ny = a->ny;
        hints.nz = a->nz;

        if (NULL==(b=vol3df_new_from_file_heuristic(bfilename,&hints))) {
                fprintf(stderr,"Failed to read %s\n",bfilename);
                return -1;
        }

        if (NULL==(c=vol3df_new_pointwise_binaryop(a,b,func))) {
                fprintf(stderr,"vol3df_new_pointwise_binaryop failed\n");
                return -1;
        }
        *vp=c;
        return 0;
}

float subtract_floats(float a, float b) { return a-b; }

int subtract_wrapper(int argc, char *argv[], int optind, struct vol3d **vp)
{
        return binaryop_wrapper(argc,argv,optind,vp,subtract_floats);
}

float add_floats(float a, float b) { return a+b; }

int add_wrapper(int argc, char *argv[], int optind, struct vol3d **vp)
{
        return binaryop_wrapper(argc,argv,optind,vp,add_floats);
}


#ifdef HAVE_FFTW
int sfactor_wrapper(int argc, char *argv[], int optind, struct vol3d **vp) {

        struct binnery *b=NULL;
        float k1,k1_err,l0;

        if (NULL==(b=vol3df_sfactor(*vp))) { return -1; }

        binnery_scaled_dump(b,stdout,2.0*M_PI/(*vp)->nx);
        binnery_scaled_1stmoment(b,2.0*M_PI/(*vp)->nx,&k1,&k1_err);
        l0=2.0*M_PI/k1;
        printf("# k1= %.10f err %.10f L= %.10f \n",k1,k1_err,l0);
        binnery_destroy(b);

        return 0;
}
#endif /* HAVE_FFTW */

/* Split domain across cdx x cdy x cdz CPUs. Set the
 * region corresponding to the CPU at (ccx,ccy,ccz) to one.
 */

int volutil_set_subdomain(struct vol3d *v,
                int cdx, int cdy, int cdz, int ccx, int ccy, int ccz)
{

        int x,y,z;
        int gx,gy,gz; /* global coords of low-xyz corner */
        int tnx,tny,tnz; /* global dimensions */
        int nx,ny,nz;    /* Size of subdomain */
        float *p=NULL;

        tnx=v->nx; tny=v->ny; tnz=v->nz;
        nx = tnx/cdx; ny = tny/cdy; nz = tnz/cdz;
        gx = ccx*nx; gy = ccy*ny; gz = ccz*nz;

        printf("set_subdomain: CPU array %dx%dx%d, CPU at (%d,%d,%d)\n",
                        cdx,cdy,cdz,
                        ccx,ccy,ccz);
        printf("tn=(%d,%d,%d) n=(%d,%d,%d), g=(%d,%d,%d)\n",
                        tnx,tny,tnz,
                        nx,ny,nz,
                        gx,gy,gz);

        if ( (0!=(tnx % cdx))
           ||(0!=(tny % cdy))
           ||(0!=(tnz % cdz))
           ) {
                fprintf(stderr,"volume not spanned.\n");
                return -1;
        }

        if (VOL_TYPE_FLOAT != v->datatype) {
                fprintf(stderr,"volume not float.\n");
                return -1;
        }

        /* Seek to low-(xyz) corner */
        p = (float *)v->data + gx + tny*gy + tnx*tny*gz;

        for (z=0;z<nz;z++) {
                for (y=0;y<ny;y++) {
                        /* Process an X-span */
                        for (x=0;x<nx;x++) { *p++ = 1.0; }

                        /* Skip tnx-nx sites to next span in XY slab */

                        p += (tnx-nx);
                }

                /* Skip (tnx*tny)-(nx*ny) sites to next XY slab in volume */

                p += ((tnx*tny)-(tnx*ny));
        }

        return 0;
        
}

int subdomain_wrapper(int argc, char *argv[], int optind, struct vol3d **vp) {

        int nargs;
        int cdx, cdy, cdz, ccx, ccy, ccz;

        nargs = argc-optind;

        if (nargs!=6) {
                fprintf(stderr,
                        "subdomain needs 6 args: cdx cdy cdz ccx ccy ccz\n");
                return -1;
        }

        cdx=atoi(argv[optind++]);
        cdy=atoi(argv[optind++]);
        cdz=atoi(argv[optind++]);
        ccx=atoi(argv[optind++]);
        ccy=atoi(argv[optind++]);
        ccz=atoi(argv[optind++]);

        return volutil_set_subdomain(*vp,cdx,cdy,cdz,ccx,ccy,ccz);

}

int downsample_wrapper(int argc, char *argv[], int optind, struct vol3d **vp) {

        int nargs;
        int subx,suby,subz;
        struct vol3d *vnew=NULL;

        nargs = argc-optind;

        if (nargs!=3) {
                fprintf(stderr,
                        "downsample needs 3 args: subx suby subz\n");
                return -1;
        }

        subx=atoi(argv[optind++]);
        suby=atoi(argv[optind++]);
        subz=atoi(argv[optind++]);

        if (NULL==(vnew=vol3df_new_downsample(*vp,subx,suby,subz))) { return -1; }
        vol3d_destroy(*vp); *vp = vnew ; vnew = NULL;

        return 0;

}


#ifdef HAVE_PNG


int pngslice_wrapper(int argc, char *argv[], int optind, struct vol3d **vp) {

        char dirchar;

        char *prefix=NULL;
        char *filename=NULL;

        unsigned int nslices=0,fnlen=0;
        enum vol3d_direction dir = VOL3D_X;
        int retval=-1;
        int coord=0;


        if (3!=(argc-optind)) {
                fprintf(stderr,
                        "pngslices needs 3 args: direction coord  prefix\n");
                return -1;
        }

        dirchar = tolower(*(argv[optind++]));

        switch(dirchar) {
                case 'x': dir = VOL3D_X; nslices=(*vp)->nx; break;
                case 'y': dir = VOL3D_Y; nslices=(*vp)->ny; break;
                case 'z': dir = VOL3D_Z; nslices=(*vp)->nz; break;
                default:
                        fprintf(stderr,"Unknown direction \"%c\".",dirchar);
                        return -1;
        }

        coord = atoi(argv[optind++]);

        prefix = argv[optind++];

        /* Filename = 
         *                prefix        strlen(prefix)
         *              + '.png'        4
         *              + NUL           1
         *
         */


        fnlen = strlen(prefix) + 5 ;
        filename = malloc(fnlen);
        snprintf(filename,fnlen,"%s.png",prefix);

        if (0!=vol3df_write_png_normalize(*vp,dir,coord,filename)) {
                fprintf(stderr,"vol3df_write_png_normalize failed.\n");
                retval=-1; goto finish;
        }

        retval=0; /* Successful completion */

finish:
        if (NULL!=filename) { free(filename); }
        return retval;


}


int pngslices_wrapper(int argc, char *argv[], int optind, struct vol3d **vp) {

        char dirchar;

        char *prefix=NULL;
        char *formatstring=NULL;
        char *filenames=NULL;
        char **fnamelist=NULL;

        unsigned int nslices=0,ndigits=0,fnlen=0;
        enum vol3d_direction dir = VOL3D_X;
        int retval=-1;
        int i;


        if (2!=(argc-optind)) {
                fprintf(stderr,
                        "pngslices needs 2 args: direction prefix\n");
                return -1;
        }

        dirchar = tolower(*(argv[optind++]));

        switch(dirchar) {
                case 'x': dir = VOL3D_X; nslices=(*vp)->nx; break;
                case 'y': dir = VOL3D_Y; nslices=(*vp)->ny; break;
                case 'z': dir = VOL3D_Z; nslices=(*vp)->nz; break;
                default:
                        fprintf(stderr,"Unknown direction \"%c\".",dirchar);
                        return -1;
        }

        prefix = argv[optind++];

        /* Filename = 
         *                prefix        strlen(prefix)
         *              + '_'           1
         *              + (digits)      ndigits
         *              + '.png'        4
         *              + NUL           1
         *
         * Format string = '%s_%0'      5
         *              + (digits)      log10(ndigits)  
         *              + 'd.png'       5
         *              + NUL           1
         */

        ndigits = (unsigned int) ceil(log10(nslices));
        ndigits = (ndigits==0) ? 1 : ndigits;

        fnlen = strlen(prefix) + 6 + ndigits;
        fprintf(stderr,"nslices=%d ndigits=%u fnlen=%d\n",nslices,ndigits,fnlen);

        /* Allocate a format string */

        if (NULL==(formatstring=malloc(11+ndigits)))
                {perror("malloc");return -1; }

        /* Build the format string */

        snprintf(formatstring,11+ndigits,"%%s_%%0%dd.png",ndigits);
        fprintf(stderr,"formatstring=<%s>\n",formatstring);

        /* Allocate space for an array of pointers to filenames */
        if (NULL==(fnamelist=malloc(sizeof(char *)*nslices))) {
                perror("malloc"); retval=-1; goto finish;
        }
        /* Allocate space for the filenames themselves */
        if (NULL==(filenames=malloc(fnlen*nslices))) {
                perror("malloc"); retval=-1; goto finish;
        }

        /* For each slice, create the filename for that slice, and
         * set up a pointer to the filename.
         */

        for (i=0;i<nslices;i++) {
                fnamelist[i]=filenames + i*fnlen; 
                snprintf(fnamelist[i],fnlen,formatstring, prefix, i);
                *(fnamelist[i]+fnlen-1)=0;
        }

        if (0!=vol3df_write_pngs(*vp,dir,fnamelist)) {
                fprintf(stderr,"vol3df_write_pngs failed.\n");
                retval=-1; goto finish;
        }

        retval=0; /* Successful completion */

finish: /* Dijkstra probably hates me, etc */
        if (NULL!=formatstring) { free(formatstring); }
        if (NULL!=fnamelist) { free(fnamelist); }
        if (NULL!=filenames) { free(filenames); }
        return retval;


        return 0;

}

#endif /* HAVE_PNG */

int ortholine_wrapper(int argc, char *argv[], int optind, struct vol3d **vp) {

        char dirchar;
        float *line=NULL;

        unsigned int nslices=0;
        unsigned int coord0; unsigned int coord1;
        enum vol3d_direction dir = VOL3D_X;
        int i;

        if (3!=(argc-optind)) {
                fprintf(stderr,
                        "orthoslice needs 3 args: direction coord0 coord1\n");
                return -1;
        }

        dirchar = tolower(*(argv[optind++]));

        switch(dirchar) {
                case 'x': dir = VOL3D_X; nslices=(*vp)->nx; break;
                case 'y': dir = VOL3D_Y; nslices=(*vp)->ny; break;
                case 'z': dir = VOL3D_Z; nslices=(*vp)->nz; break;
                default:
                        fprintf(stderr,"Unknown direction \"%c\".",dirchar);
                        return -1;
        }

        coord0 = (unsigned)atoi(argv[optind++]);
        coord1 = (unsigned)atoi(argv[optind++]);



        if (NULL==(line=vol3df_ortholine(*vp,dir,coord0,coord1))) {
            fprintf(stderr,"ortholine_wrapper() failed.\n");
            return -1;
        }

        for (i=0;i<nslices;i++)
            { fprintf(stdout,"%f\n",line[i]); }

        free(line);

        return 0;

}

struct operation {
        char *name;
        char *summary;
        int (*op_fn)(int , char *[], int, struct vol3d **);
        int need_input;
        int need_output;
};

struct operation oplist[] = {
        { "null",NULL,
                null_wrapper, 1, 1 },
        { "swapzx",NULL,
                swapzx_wrapper,  1, 1 },
        { "rotxyz",NULL,
                rotxyz_wrapper,  1, 1 },
        { "erode",NULL,
                erode_wrapper,  1, 1 },
        { "dilate",NULL,
                dilate_wrapper,  1, 1 },
        { "open",NULL,
                open_wrapper, 1, 1 },
        { "close",NULL,
                close_wrapper, 1, 1 },
        { "nngt","<threshold>",
                nngt_wrapper, 1, 1 },
        { "nnlt", "<threshold>",
                nnlt_wrapper, 1, 1 },
        { "subvol","<dx> <dy> <dz> [<x> <y> <z>]",
                subvol_wrapper, 1, 1 },
        { "upscale","<nx> <ny> <nz>",
                upscale_wrapper, 1, 1 },
        { "scale","<nx> <ny> <nz>",
                scale_wrapper, 1, 1 },
        { "gt","<threshold>",
                gt_wrapper, 1, 1 },
        { "lt","<threshold>",
                lt_wrapper, 1, 1 },
        { "range","<max> <min>",
                range_wrapper, 1, 1 },
        { "rangefilter","<max> <min>",
                rangefilter_wrapper, 1, 1 },
        { "abs",NULL,
                abs_wrapper, 1, 1 },
        { "stats",NULL,
                stats_wrapper, 1, 0 },
        { "ortholine",NULL,
                ortholine_wrapper, 1, 0 },
        { "dist","[<nbins>]",
                dist_wrapper, 1, 0 },
        { "smear","[<delta>]",
                smear_wrapper, 1, 1},
        { "iso","[<isoval>]",
                NULL, 1, 0},
        { "isovtk","[<isoval>]",
                NULL, 1, 0},
        { "fill","[<fillval>]",
                fill_wrapper, 0, 1},
        { "setpoint","<x> <y> <z> <val>",
                setpoint_wrapper, 0, 1},
        { "gyroid","[<nx>[ <ny>[ <nx>]]]",
                gyroid_wrapper, 0, 1},
        { "psurface","[<nx>[ <ny>[ <nx>]]]",
                psurface_wrapper, 0, 1},
        { "dsurface","[<nx>[ <ny>[ <nx>]]]",
                dsurface_wrapper, 0, 1},
        { "gaussian","[<nx>[ <ny>[ <nx>]]]",
                gaussian_wrapper, 0, 1},
        { "sine","[<nx>[ <ny>[ <nx>]]]",
                sine_wrapper, 0, 1},
        { "odist","[<nx>[ <ny>[ <nx>]]]",
                odist_wrapper, 0, 1},
        { "ogauss","[<nx>[ <ny>[ <nx>]]]",
                ogauss_wrapper, 0, 1},
        { "cylrad",NULL,
                cylrad_wrapper, 0, 1},
        { "subtract","<filename> [<filetype>]",
                subtract_wrapper,1,1},
        { "add","<filename> [<filetype>]",
                add_wrapper,1,1},
        { "total_surfeuler", "[<threshold>]",
                total_surf_euler_wrapper, 1, 0},
        { "total_voleuler", "[<threshold>]",
                total_vol_euler_wrapper, 1, 0},
        { "surfeuler", "[<threshold>]",
                surf_euler_wrapper, 1, 1},
        { "voleuler", "[<threshold>]",
                vol_euler_wrapper, 1, 1},
        { "roll", "[<dx> [<dy> [<dz>]]]",
                roll_wrapper, 1, 1},
        { "surf_euler8","[<threshold>]",
                surf_euler8_wrapper, 1, 1 },
        { "vol_euler8","[<threshold>]",
                vol_euler8_wrapper, 1, 1 },
        { "lamella","<spacing> <width>",
                lamella_wrapper, 0, 1},
        { "downsample","<subx> <suby> <subz>",
                downsample_wrapper, 1, 1},
        { "shear","",
                shear_wrapper, 1, 1 },
        { "rot111","",
                rot111_wrapper, 1, 1 },
        { "whitenoise","[<magnitude> [<zero>]]",
                whitenoise_wrapper, 0, 1},
#ifdef HAVE_PNG
        { "pngslices","<direction [xyz]> <prefix>",
                pngslices_wrapper, 1, 0 },
        { "pngslice","<direction [xyz]> <index> <prefix>",
                pngslice_wrapper, 1, 0 },
#endif /* HAVE_PNG */
#ifdef HAVE_FFTW
        { "rfft","",
                rfft_wrapper, 1, 1 },

        { "spectrum","",
                spectrum_wrapper, 1, 1 },
        { "sfactor","",
                sfactor_wrapper, 1, 0 },
        { "autocorrelation","",
                autocorrelation_wrapper,1,1},
        { "ifnoise","",
                ifnoise_wrapper, 0, 1},
#endif /* HAVE_FFTW */

        { "logabs","",
                logabs_wrapper, 1, 1 },
        { "zerobounds","",
                zerobounds_wrapper, 1, 1 },
        { "subdomain","<cdx> <cdy> <cdz> <ccx> <ccy> <ccz>",
                subdomain_wrapper, 1, 1 },
        { "H",NULL,
                H_wrapper,  1, 1 },
        { "areadensity","[<isoval>]",
                areadensity_wrapper, 1, 1 },
        { "h2density","[<isoval>]",
                h2density_wrapper, 1, 1 },
        { NULL, NULL, 0, 0 }
};

void usage(char *name)
{
        struct operation *op=NULL;

        fprintf(stderr,"Usage: %s <options> [<operation>]\n"
                "\t-o <output filename> (required)\n"
                "\t-i <input filename>  (required)\n"
                "\t-d <nx>,<ny>,<nz> -- specify dimensions\n"
                "\t-p <HDF5 path> -- specify path inside HDF5 file\n"
                "\t-f <intype> -- specify input file format (default xdr)\n"
                "\t-g <outtype> -- specify output file format (default xdr)\n"
                "\t-h,-? -- print this message\n"
                "\nPossible input types are "
#ifdef HAVE_HDF5
                "hdf5,"
#endif /* HAVE_HDF5 */
                "raw,iso,xdr\n"
                "Possible output types are raw,xdr,iso,"
#ifdef HAVE_HDF5
                "hdf5,"
#endif /* HAVE_HDF5 */
                "vtk,ucvtk,ucnrrd,fnrrd,uc\n"
                ,name
                );
        fprintf(stderr, "Possible operations are:\n");

        for (op=oplist; NULL!=op->name ; op++) {
                fprintf(stderr,"\t%s\t%s\n", op->name,
                                (NULL==op->summary) ? "" : op->summary);
        }

        return;
}

int main(int argc, char *argv[])
{

        int c;
        const char *optstring="ni:o:d:f:g:h?p:";
        char *outfilename=NULL; 
        char *infilename=NULL; 
        enum vol3df_filetype itype = VOL3DF_FTYPE_UNKNOWN;
        enum vol3df_filetype otype = VOL3DF_FTYPE_UNKNOWN;
        struct vol3d *v=NULL;
        unsigned int nx=0,ny=0,nz=0;
        char *h5path=NULL;
        char *opstr="null"; /* String describing operation to be performed */
        int need_output=0, need_input = 1 ;
        int (*op_fn)(int , char *[], int, struct vol3d **)=NULL;
        int calc_vnormals=0;

        struct vol3df_file_hints hints = vol3df_default_file_hints;



        if (1==argc) { usage(argv[0]); return -1; }

        while (-1!=(c=getopt(argc,argv,optstring))) {
                switch(c) {
                        case 'n':
                                calc_vnormals = 1;
                                break;
                        case 'i':
                                hints.filename = infilename=optarg;
                                break;
                        case 'o':
                                outfilename=optarg;
                                break;
                        case 'd':
                                vol3d_parse_dims_string(&nx,&ny,&nz,optarg);
                                hints.nx = nx;
                                hints.ny = ny;
                                hints.nz = nz;
                                break;
                        case 'f':
                                itype = vol3df_parse_ftype(optarg);
                                if (
                                        (VOL3DF_FTYPE_RAW != itype)
                                        && (VOL3DF_FTYPE_XDR != itype)
                                        && (VOL3DF_FTYPE_ISO != itype)
                                        && (VOL3DF_FTYPE_HDF5 != itype)
                                ) {
                                        fprintf(stderr,
                                              "Unsupported input file type\n");
                                        return -1;
                                }
                                hints.type = itype;
                                break;
                        case 'g':
                                otype = vol3df_parse_ftype(optarg);
                                if (VOL3DF_FTYPE_UNKNOWN == otype) {
                                        fprintf(stderr,
                                              "Unknown output file type\n");
                                        return -1;
                                }
                                break;
                        case 'p':
                                h5path = hints.h5path = optarg;
                                break;
                        case 'h':
                        case '?':
                                usage(argv[0]);
                                return -1;
                        default:
                                fprintf(stderr,"Unknown option: %c\n",optopt);
                                return -1;
                }
        }


        /* Now, look at the rest of the arguments. */

        if (optind == argc) { 
                /* No operation; do nothing. */
                op_fn = null_wrapper; need_input = need_output = 1;
        } else {
                /* Determine which operation to perform, whether
                 * it needs an input file to operate on, and whether
                 * it produces an output file.
                 */
                struct operation *op;

                opstr = argv[optind++];

                for (op=oplist; NULL!=op->name ; op++) {
                        if (0==strcmp(opstr,op->name)) {
                                op_fn = op->op_fn;
                                need_input = op->need_input;
                                need_output = op->need_output;
                                break;
                        }
                }
                if (NULL==op->name) {
                        fprintf(stderr,"Unknown operation %s\n",opstr);
                        return -1;
                }
        }

        /* op_fn, need_(in|out)put defined */

        /* Load the file even if need_input is not set, since the
         * dimensions for the output file can then be guessed from
         * the input file.
         */

        if ((need_input)||(NULL!=infilename)) {

                if (0==strcmp("-",infilename)) {
                        /* Read from stdin */
                        if (NULL==(v=vol3df_new_from_fh_heuristic(stdin,&hints)))
                        {
                                fprintf(stderr,"Failed to read input stream\n");
                                return -1;
                        }
                } else {

                if (NULL==(v=vol3df_new_from_file_heuristic(infilename,&hints)))
                {
                        fprintf(stderr,"Failed to read input file\n");
                        return -1;
                }
                }

                nx=v->nx; ny=v->ny; nz=v->nz;

                fprintf(stderr,"Loaded %dx%dx%d\n",nx,ny,nz);
        } else {
                if (0!=nx) {
                        /* We've been given dimensions; create a new
                         * volume of that size.
                         */
                        if (NULL==(v=vol3df_new(nx,ny,nz))) {
                                fprintf(stderr,"vol3df_new failed\n");
                                return -1;
                        }
                        nx=v->nx; ny=v->ny; nz=v->nz;

                } else {
                        fprintf(stderr,"Please specify dimensions with -d\n");
                        return -1;
                }


        }

        /* Perform operation */

        /* Isosurfacing is a special case. */

        if ((0==strcmp("isovtk",argv[optind-1])) 
         ||(0==strcmp("iso",argv[optind-1]))) {

                struct mcubes_state *mc=NULL;
                float isoval = 0.0;
                FILE *outfile=NULL;

                if (argc != optind) { isoval = atof(argv[optind++]); }
                if (NULL==(mc=mcubes_new(v,isoval,NULL,1))) {
                        fprintf(stderr,"mcubes_new failed.\n");
                        return -1;
                }
                if (0!=mcubes_calc_vertices(mc)) {
                        fprintf(stderr,"mcubes_calc_vertices failed.\n");
                        return -1;
                }


                mcubes_lewiner(mc);

                if(NULL==(outfile=fopen(outfilename,"wb")))
                        { perror("fopen"); return -1; }

                if (0==strcmp("isovtk",argv[optind-1])) {
                    if (0!=mcubes_write_vtkfh(mc,outfile,1)) {
                            fprintf(stderr,"mesh_write_vtkfh failed\n");
                            return -1;
                    }
                } else {
                    if (0!=mcubes_write_surf_fh(mc,outfile)) {
                        fprintf(stderr,"mcubes_write_surf failed\n");
                        return -1;
                    }
                }

                mcubes_destroy(mc);

        } else if (0!= op_fn(argc,argv,optind,&v)) {
                fprintf(stderr,"Operation \"%s\" failed.\n",opstr);
                return -1;
        }

        if (NULL==v) {
                fprintf(stderr,"Internal error: NULL v after operation\n");
                return -1;
        }

        /* At this point, v is defined and contains processed data.
         * Write this to a file if necessary.
         */

        if (need_output) {
                int retval;
                if (NULL==outfilename) {
                        fprintf(stderr,"Please specify output filename\n");
                        return -1;
                }
                retval = write_output(v,outfilename,otype,h5path);
                vol3d_destroy(v);
                return retval;
        } else {

                vol3d_destroy(v);
                exit(0);
                return 0;
        }

}

