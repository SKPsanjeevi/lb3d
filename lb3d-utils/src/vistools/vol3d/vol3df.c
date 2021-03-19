/* Subclass of vol3d for floating-point data */

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>


#include <rpc/rpc.h>

#ifdef HAVE_HDF5
#include <hdf5.h>
#endif

#ifdef HAVE_FFTW
#include <fftw3.h>
#endif

#include "vol3d.h"
#include "vol3df.h"

#include "vol2df.h"

#include "binnery.h"

/* Default set of file hints, representing the state of knowing
 * nothing about a file.
 */

const struct vol3df_file_hints vol3df_default_file_hints = {
        VOL3DF_FTYPE_UNKNOWN,   /* type */
        NULL,                   /* filename */
        NULL,                   /* h5path */
        0,0,0                   /* dimensions */
};


/* Lookup table for 8x the Euler number for the surface intersecting the
 * 1^3 region surrounding the central vertex of a 2x2x2 voxel cluster;
 * factor of eight kept for compatibility with volume Euler number.
 *
 * So, this table is the one you use if you want to calculate the Euler number
 * of a *surface* , the sum of V-E+F.
 */


int surf_euler_lookup[256] = {
0,2,2,0,2,0,0,-2,2,0,0,-2,0,-2,-2,0,2,0,0,-2,0,-2,2,-4,-4,-2,
-2,-4,-2,-4,0,-2,2,0,0,-2,-4,-2,-2,-4,0,2,-2,-4,-2,0,-4,-2,0,
-2,-2,0,-2,-4,0,-2,-2,0,-4,-2,0,-2,-2,0,2,0,-4,-2,0,-2,-2,-4,
0,2,-2,0,-2,-4,-4,-2,0,-2,-2,-4,-2,0,0,-2,-2,0,0,-2,-4,-2,-2,
0,0,2,-2,0,-2,0,0,-2,2,8,0,2,0,2,-2,0,-2,-4,-4,-2,-4,-2,-2,0,
0,2,-2,0,-2,0,-4,2,2,-4,0,-2,0,-2,2,0,0,-2,-2,-4,-2,-4,-4,-2,
0,-2,2,0,2,0,8,2,-2,0,0,-2,0,-2,2,0,0,-2,-2,-4,-2,0,0,-2,-2,
0,0,-2,-4,-2,-2,0,-2,-4,-4,-2,0,-2,2,0,-4,-2,-2,0,-2,-4,0,2,
0,-2,-2,0,-2,-4,0,-2,-2,0,-4,-2,0,-2,-2,0,-2,-4,0,-2,-4,-2,2,
0,-4,-2,-2,-4,-2,0,0,2,-2,0,-4,-2,-4,-2,-2,-4,-4,2,-2,0,-2,0,
0,2,0,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,0
};

/* Lookup table for 8x the Euler number for the volume intersecting the
 * 1^3 region surrounding the central vertex of a 2x2x2 voxel cluster;
 * factor of eight ensures that this is always integral.
 *
 * This table is the one you use if you want to calculate the Euler number of
 * a *volume*, the sum of V-E+F-C.
 */

const int vol_euler_lookup[256] = {

0,1,1,0,1,0,-2,-1,1,-2,0,-1,0,-1,-1,0,1,0,-2,-1,-2,-1,-1,-2,
-6,-3,-3,-2,-3,-2,0,-1,1,-2,0,-1,-6,-3,-3,-2,-2,-1,-1,-2,-3,
0,-2,-1,0,-1,-1,0,-3,-2,0,-1,-3,0,-2,-1,0,1,1,0,1,-2,-6,-3,
0,-1,-3,-2,-2,-1,-3,0,-1,-2,-2,-1,0,-1,-3,-2,-1,0,0,-1,-3,
0,0,1,-2,-1,1,0,-2,-1,-3,0,-3,0,0,1,-1,4,0,3,0,3,1,2,-1,-2,
-2,-1,-2,-1,1,0,0,3,1,2,1,2,2,1,1,-6,-2,-3,-2,-3,-1,0,0,-3,
-1,-2,-1,-2,-2,-1,-2,-3,-1,0,-1,0,4,3,-3,0,0,1,0,1,3,2,0,-3,
-1,-2,-3,0,0,1,-1,0,0,-1,-2,1,-1,0,-1,-2,-2,-1,0,1,3,2,-2,1,
-1,0,1,2,2,1,0,-3,-3,0,-1,-2,0,1,-1,0,-2,1,0,-1,-1,0,-1,-2,
0,1,-2,-1,3,2,-2,1,1,2,-1,0,2,1,-1,0,-2,1,-2,1,1,2,-2,3,-1,2,
-1,2,0,1,0,-1,-1,0,-1,0,2,1,-1,2,0,1,0,1,1,0

};


/* euler_delta[i] contains the offset vector for the ith voxel in
 * the 2x2x2 Euler-lookup cluster.
 */
const int euler_delta[8][3] = {
        { 0, 0, 0},
        { 1, 0, 0},
        { 0, 1, 0},
        { 1, 1, 0},
        { 0, 0, 1},
        { 1, 0, 1},
        { 0, 1, 1},
        { 1, 1, 1}
};


void vol3df_unsupported_mess(char *featurename)
{
        fprintf(stderr,"%s support was not compiled in.\n",featurename);
}

struct vol3d *vol3df_new(unsigned int nx, unsigned int ny, unsigned int nz)
{
        return vol3d_new(nx,ny,nz,VOL_TYPE_FLOAT);
}

struct vol3d *vol3df_new_zero(unsigned int nx, unsigned int ny, unsigned int nz)
{
        struct vol3d *v=NULL;
        float *p=NULL;
        int i=0;

        if (NULL==(v=vol3d_new(nx,ny,nz,VOL_TYPE_FLOAT))) {
                fprintf(stderr,"vol3df_new_zero: vol3df_new failed\n");
                return NULL;
        }
        p=(float *)v->data;

        for (i=0;i<nx*ny*nz;i++) { *p++=0.0; }

        return v;

}

void vol3df_xdr_streamfunc(void *el,void *data)
{

        xdr_float((XDR *)data,(float *)el);

}


/* Look at a string (eg from a file extension) and try to guess a filetype */

enum vol3df_filetype vol3df_parse_ftype(char *s)
{
        if (0==strcmp("raw",s)) {
                return VOL3DF_FTYPE_RAW;
        } else if (0==strcmp("xdr",s)) {
                return VOL3DF_FTYPE_XDR;
        } else if (0==strcmp("hdf5",s)) {
                return VOL3DF_FTYPE_HDF5;
        } else if (0==strcmp("h5",s)) {
                return VOL3DF_FTYPE_HDF5;
        } else if (0==strcmp("vtk",s)) {
                return VOL3DF_FTYPE_VTK;
        } else if (0==strcmp("ucvtk",s)) {
                return VOL3DF_FTYPE_UCVTK;
        } else if (0==strcmp("uc",s)) {
                return VOL3DF_FTYPE_UC;
        } else if (0==strcmp("iso",s)) {
                return VOL3DF_FTYPE_ISO;
        } else if (0==strcmp("ucnrrd",s)) {
                return VOL3DF_FTYPE_UCNRRD;
        } else if (0==strcmp("fnrrd",s)) {
                return VOL3DF_FTYPE_FNRRD;
        } else if (0==strcmp("nrrd",s)) {
                return VOL3DF_FTYPE_FNRRD; /* Assume floating-point nrrd */
        } else {
                fprintf(stderr,"vol3df_parse_ftype: don't recognize \"%s\".\n",
                               s);
                return VOL3DF_FTYPE_UNKNOWN;
        }

}

/* Look at a filename and guess its type based on
 * extension. Return VOL3DF_FTYPE_UNKNOWN if unable to guess.
 */

enum vol3df_filetype vol3df_guess_fname_ftype(char *s)
{
        size_t len;
        char *ext=NULL;
        int i;

        len=strlen(s);

        /* Find last period */

        for (i=len-1;i>=0;i--) {
                if ('.'==s[i]) { break; }
        }

        if (i<0) { return VOL3DF_FTYPE_UNKNOWN; } /* No extension */

        /* We have an extension */
        ext=s+i+1;

        return vol3df_parse_ftype(ext);

}

/* Returns gradient of field */
void vol3df_gradient(struct vol3d *v,int x, int y, int z,
                float *dxp, float *dyp, float *dzp)
{

        *dxp = 0.5*( vol3df_wrap(v,x+1,y  ,z  ) - vol3df_wrap(v,x-1,y  ,z  ));
        *dyp = 0.5*( vol3df_wrap(v,x  ,y+1,z  ) - vol3df_wrap(v,x  ,y-1,z  ));
        *dzp = 0.5*( vol3df_wrap(v,x  ,y  ,z+1) - vol3df_wrap(v,x  ,y  ,z-1));

}

void vol3df_hessian(struct vol3d *v, int x, int y, int z,
        float *hxx,float *hxy, float *hxz, float *hyy, float *hyz, float *hzz)
{
        float cent = -2.0*vol3df_wrap(v,x,y,z); /* -2*(centre value) */
        *hxx = cent + vol3df_wrap(v,x+1,y  ,z  )+vol3df_wrap(v,x-1,y  ,z  );
        *hyy = cent + vol3df_wrap(v,x  ,y+1,z  )+vol3df_wrap(v,x  ,y-1,z  );
        *hzz = cent + vol3df_wrap(v,x  ,y  ,z+1)+vol3df_wrap(v,x  ,y  ,z-1);
        *hxy = 0.25*(
                 vol3df_wrap(v,x+1,y+1,z  )
                -vol3df_wrap(v,x+1,y-1,z  )
                -vol3df_wrap(v,x-1,y+1,z  )
                +vol3df_wrap(v,x-1,y-1,z  )
                );
        *hxz = 0.25*(
                 vol3df_wrap(v,x+1,y  ,z+1)
                -vol3df_wrap(v,x+1,y  ,z-1)
                -vol3df_wrap(v,x-1,y  ,z+1)
                +vol3df_wrap(v,x-1,y  ,z-1)
                );
        *hyz = 0.25*(
                 vol3df_wrap(v,x  ,y+1,z+1)
                -vol3df_wrap(v,x  ,y+1,z-1)
                -vol3df_wrap(v,x  ,y-1,z+1)
                +vol3df_wrap(v,x  ,y-1,z-1)
                );
}

void vol3df_trilinear_hessian(struct vol3d *v, float x, float y, float z,
 float *hxxp,float *hxyp, float *hxzp, float *hyyp, float *hyzp, float *hzzp)
{
    int lx,ly,lz;
    int dxi,dyi,dzi;
    float phixc[2],phiyc[2],phizc[2];
    float dx,dy,dz;

    lx = (int) floor(x); ly = (int) floor(y); lz = (int) floor(z);
    dx = x-lx; dy=y-ly; dz=z-lz;

    phixc[0] = 1.0 - dx; phixc[1] = dx; 
    phiyc[0] = 1.0 - dy; phiyc[1] = dy; 
    phizc[0] = 1.0 - dz; phizc[1] = dz; 

    *hxxp=*hyyp=*hzzp=*hxyp=*hxzp=*hyzp=0.0;
    

    for (dzi=0;dzi<=1;dzi++) {
    for (dyi=0;dyi<=1;dyi++) {
    for (dxi=0;dxi<=1;dxi++) {
        float hxx,hyy,hzz,hxy,hxz,hyz;
        vol3df_hessian(v,lx+dxi,ly+dyi,lz+dzi,&hxx,&hxy,&hxz,&hyy,&hyz,&hzz);

        *hxxp += phixc[dxi]*phiyc[dyi]*phizc[dzi] * hxx;
        *hyyp += phixc[dxi]*phiyc[dyi]*phizc[dzi] * hyy;
        *hzzp += phixc[dxi]*phiyc[dyi]*phizc[dzi] * hzz;
        *hxyp += phixc[dxi]*phiyc[dyi]*phizc[dzi] * hxy;
        *hxzp += phixc[dxi]*phiyc[dyi]*phizc[dzi] * hxz;
        *hyzp += phixc[dxi]*phiyc[dyi]*phizc[dzi] * hyz;
    }}}


    return;
}

matrix vol3df_trilinear_hessian_matrix(struct vol3d *v,
        float x, float y, float z)
{
    matrix m;

    vol3df_trilinear_hessian(v,x,y,z,
            &m.el[0],&m.el[1],&m.el[2],&m.el[4],&m.el[5],&m.el[8]);
    m.el[3] = m.el[1];
    m.el[6] = m.el[2];
    m.el[7] = m.el[5];
    return m;

}

vector vol3df_gradientvec(struct vol3d *v,int x, int y, int z)
{
    vector vec;

    vec.x = 0.5*( vol3df_wrap(v,x+1,y  ,z  ) - vol3df_wrap(v,x-1,y  ,z  ));
    vec.y = 0.5*( vol3df_wrap(v,x  ,y+1,z  ) - vol3df_wrap(v,x  ,y-1,z  ));
    vec.z = 0.5*( vol3df_wrap(v,x  ,y  ,z+1) - vol3df_wrap(v,x  ,y  ,z-1));
    return vec;

}

vector vol3df_trilinear_gradient_vec(struct vol3d *v, float rx, float ry, float rz)
{
    vector vec;
    vol3df_trilinear_gradient(v,rx,ry,rz,&vec.x,&vec.y,&vec.z);
    return vec;
}

/* Assume that (rx,ry,rz) is a point on some surface of constant phi.
 * Return -0.5*del . ( del phi / ( | del phi | ) ) , equal
 * to the mean curvature of the surface at that point.
 * See, eg, Gozdz & Holyst, PRL 76, 15 p2726--2729 (1996).
 */

float vol3df_mean_curvature(struct vol3d *v, float rx, float ry, float rz)
{
    return 0.5*matrix_trace(vol3df_geometry_tensor(v,rx,ry,rz));
    /*
    float S[9];
    vol3df_shape_operator(v,S,rx,ry,rz);
    return -0.5*(S[0]+S[4]+S[8]);
    */
    /*
    float delta=0.01;
    float sum=0.0;
    float dx,dy,dz;

    vector C,C_x, C_y, C_z;

    / C(r) = del phi / |del phi| /

    C  =vector_normalized(vol3df_trilinear_gradient_vec(v,rx,ry,rz));
    C_x=vector_normalized(vol3df_trilinear_gradient_vec(v,rx+delta,ry,rz));
    C_y=vector_normalized(vol3df_trilinear_gradient_vec(v,rx,ry+delta,rz));
    C_z=vector_normalized(vol3df_trilinear_gradient_vec(v,rx,ry,rz+delta));

    dx = (C_x.x - C.x)/delta;
    dy = (C_y.y - C.y)/delta;
    dz = (C_z.z - C.z)/delta;

    sum = -0.5 * ( dx + dy + dz);
    return sum;
    */

}

/* Find the (approximate) matrix of negative derivatives of the normal field,
 * at a point.  The Gaussian curvature is the determinant of this matrix; the
 * mean curvature is half of the trace.  See
 * http://mathworld.wolfram.com/ShapeOperator.html
 */

void vol3df_shape_operator(struct vol3d *v, float *S,
        float rx, float ry, float rz)
{

    float H[3][3]; /* Hessian */
    float u[3],Hu[3]; /* gradient, Hessian dot gradient */
    float r; /* Magnitude of gradient */
    float hxx,hxy,hxz,hyy,hyz,hzz;
    float ir3,ir;

    int i,j;

    vol3df_trilinear_gradient(v,rx,ry,rz,u,u+1,u+2);

    r = sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
    ir = 1.0/r;
    ir3 = 1.0/(r*r*r);

    vol3df_trilinear_hessian(v,rx,ry,rz, &hxx,&hxy,&hxz,&hyy,&hyz,&hzz);
    H[0][0] = hxx;
    H[0][1] = H[1][0] = hxy;
    H[0][2] = H[2][0] = hxz;
    H[1][1] = hyy;
    H[1][2] = H[2][1] = hyz;
    H[2][2] = hzz;

    for (i=0;i<3;i++) {
        Hu[i] = H[i][0]*u[0] + H[i][1]*u[1] + H[i][2]*u[2];
    }

    for (i=0;i<3;i++) {
    for (j=0;j<3;j++) {
        *(S+3*i+j) = ir3 * u[i]*Hu[j] - ir * H[i][j];
    }}


}

/*
 * See http://www.cs.utah.edu/~gk/papers/vis03/talk/slide016.html
 *
 * Define a surface normal n.
 * Define the projection matrix P = I - nn
 * Let the Hessian be H.
 * Let the gradient vector be g
 * The geometry tensor is G = -PHP/|g|.
 * Tr(G) = c1+c2 , twice the mean curvature.
 * frob(G) = sqrt(k1*k1 + k2*k2)
 */

matrix vol3df_geometry_tensor(struct vol3d *v, float x, float y, float z)
{
    vector n,g;
    matrix P,H,G;
    float scalefactor;

    g = vol3df_trilinear_gradient_vec(v,x,y,z);
    scalefactor = 1.0/vector_modulus(g);
    n = vector_normalized(g);

    P = matrix_sub(matrix_identity(), matrix_outer_vectors(n,n));

    H = vol3df_trilinear_hessian_matrix(v,x,y,z);

    G = matrix_scale( matrix_product(P,matrix_product(H,P)),scalefactor);

    return G;
}

int vol3df_curvatures(struct vol3d *v,float rx, float ry, float rz,
        float *H,float *K) {
    matrix G;
    float Tr,Fr;

    G = vol3df_geometry_tensor(v,rx,ry,rz);
    Tr = matrix_trace(G);
    Fr = matrix_frob(G);

    *K = 0.5*( Tr*Tr - Fr*Fr);
    *H = 0.5*Tr;
    return 0;
}


float vol3df_gaussian_curvature(struct vol3d *v,float rx, float ry, float rz) {
    matrix G;
    float Tr,Fr;

    G = vol3df_geometry_tensor(v,rx,ry,rz);
    Tr = matrix_trace(G);
    Fr = matrix_frob(G);

    return 0.5*( Tr*Tr - Fr*Fr);
}


void vol3df_trilinear_gradient(struct vol3d *v,float x, float y, float z,
                float *gxp, float *gyp, float *gzp)
{

    int lx,ly,lz;
    int dxi,dyi,dzi;
    float phixc[2],phiyc[2],phizc[2];
    float dpx=0.0,dpy=0.0,dpz=0.0;
    float dx,dy,dz;

    lx = (int) floor(x); ly = (int) floor(y); lz = (int) floor(z);
    dx = x-lx; dy=y-ly; dz=z-lz;

    phixc[0] = 1.0 - dx; phixc[1] = dx; 
    phiyc[0] = 1.0 - dy; phiyc[1] = dy; 
    phizc[0] = 1.0 - dz; phizc[1] = dz; 

    for (dzi=0;dzi<=1;dzi++) {
    for (dyi=0;dyi<=1;dyi++) {
    for (dxi=0;dxi<=1;dxi++) {
        float phix,phiy,phiz;
        vol3df_gradient(v,x+dxi,y+dyi,z+dzi,&phix,&phiy,&phiz);
        dpx += phixc[dxi]*phiyc[dyi]*phizc[dzi] * phix;
        dpy += phixc[dxi]*phiyc[dyi]*phizc[dzi] * phiy;
        dpz += phixc[dxi]*phiyc[dyi]*phizc[dzi] * phiz;
    }}}

    *gxp=dpx;
    *gyp=dpy;
    *gzp=dpz;

    return;
}


float vol3df_trilinear_point(struct vol3d *v, float x, float y, float z)
{

    int lx,ly,lz;
    int dxi,dyi,dzi;
    float phi=0.0;
    float phixc[2],phiyc[2],phizc[2];

    lx = (int) floor(x); ly = (int) floor(y); lz = (int) floor(z);

    phixc[0] = 1.0 - (x-lx); phixc[1] = x-lx; 
    phiyc[0] = 1.0 - (y-ly); phiyc[1] = y-ly; 
    phizc[0] = 1.0 - (z-lz); phizc[1] = z-lz; 

    for (dzi=0;dzi<=1;dzi++) {
    for (dyi=0;dyi<=1;dyi++) {
    for (dxi=0;dxi<=1;dxi++) {
        phi += phixc[dxi]*phiyc[dyi]*phizc[dzi]
               *vol3df_wrap(v,lx+dxi,ly+dyi,lz+dzi);
    }}}

    return phi;
}

/* Take a filehandle and some hints about what it might contain.
 * Attempt to construct a volume from this information.
 */

struct vol3d *vol3df_new_from_fh_heuristic(
                FILE *f,
                struct vol3df_file_hints *hints
                )
{
        struct vol3df_file_hints myhints = vol3df_default_file_hints;

        if (NULL==hints) { hints = &myhints; }

        /* Is type known? */

        if (VOL3DF_FTYPE_UNKNOWN == hints->type) {
                /* No ; try to guess. Is filename known? */
                if (NULL==hints->filename) {
                        /* No. Return an error.
                         * FIXME: one could try to guess from, eg,
                         * the first 128 bytes, but this might require
                         * rewinding of filehandles etc.
                         */
                        fprintf(stderr,
                                "vol3df_new_from_fh_heuristic: "
                                "Unable to guess file type; no filename.\n");
                        return NULL;
                } else {
                        hints->type =vol3df_guess_fname_ftype(hints->filename);
                        if (VOL3DF_FTYPE_UNKNOWN==hints->type) {
                                fprintf(stderr,
                                        "vol3df_new_from_fh_heuristic: "
                                        "Unable to guess file type\n");
                                return NULL;
                        }
                }
        }

        /* Type is now known, or at least guessed. */

        if (VOL3DF_FTYPE_HDF5==hints->type) {
                if (NULL==(hints->filename)) {
                        fprintf(stderr,"vol3df_new_from_fh_heuristic: "
                                "need to know HDF5 filename\n");
                        return NULL;
                }


                /* HDF5, for all its sins, contains all the metadata
                 * you could ever want, so just go ahead and load it.
                 */

                return vol3df_new_fromhdf5(hints->filename,hints->h5path);
        } else if (VOL3DF_FTYPE_ISO==hints->type) {
                if (NULL==(hints->filename)) {
                        fprintf(stderr,"vol3df_new_from_fh_heuristic: "
                                "need to know ISO filename\n");
                        return NULL;
                }

                return vol3df_new_from_isofh(f);
        }

        /* Otherwise, we need to know the dimensions. */

        if (0==hints->nx*hints->ny*hints->nz) {
                /* Dimensions not specified.
                 * If the file is XDR or raw, and we know the
                 * filename, then we can try to guess from size.
                 */

                if ( (NULL!=hints->filename) && (
                                (VOL3DF_FTYPE_RAW==hints->type)
                             || (VOL3DF_FTYPE_XDR==hints->type)
                     )) {
                        /* We can try to guess from size. */

                        unsigned long nvoxels;
                        unsigned long csize;
                        unsigned long voxelsize;
                        struct stat statbuf;

                        if (VOL3DF_FTYPE_XDR==hints->type) {
                                voxelsize = 4; /* size of XDR float */
                        } else {
                                voxelsize = sizeof(float);
                        }

                        if (0!=stat(hints->filename,&statbuf)) {
                                perror("vol3df_new_from_fh_heuristic: stat");
                                return NULL;
                        }
                        /* We have filesize. */

                        nvoxels = (unsigned long) statbuf.st_size/voxelsize;

                        csize=(unsigned long)ceil(pow((double)nvoxels,1.0/3.0));

                        /* If the number of voxels is exactly a cube,
                         * then guess cubic dimensions. Else give up.
                         */

                        if (csize*csize*csize==nvoxels) {
                                fprintf(stderr,"vol3df_new_from_fh_heuristic:"
                                               " guessing %lux%lux%lu\n",
                                               csize,csize,csize);
                                hints->nx = hints->ny = hints->nz = csize;
                        } else {
                                fprintf(stderr,"vol3df_new_from_fh_heuristic:"
                                               " unable to guess file size;"
                                               " nvoxels=%lu.\n",nvoxels);
                                return NULL;
                        }

                } else {

                        fprintf(stderr,"vol3df_new_from_fh_heuristic: "
                                       "need to know dataset dimensions.\n");
                        return NULL;
                }
        }

        /* We know filetype and dimensions. Load it. */

        switch(hints->type) {
                case VOL3DF_FTYPE_RAW:
                        return vol3df_new_from_rawfh(
                                        hints->nx,hints->ny,hints->nz,f);
                case VOL3DF_FTYPE_XDR:
                        return vol3df_new_from_xdrfh(
                                        hints->nx,hints->ny,hints->nz,f);
                default:
                        fprintf(stderr,"vol3df_new_from_fh_heuristic: "
                                       "unable to load this filetype\n");
                        return NULL;
        }

        /* Unreachable */

}

/* Attempt to open given filename; accept hints about what it contains. */

struct vol3d *vol3df_new_from_file_heuristic(
                char *filename,
                struct vol3df_file_hints *hints
                )
{
        struct vol3df_file_hints myhints = vol3df_default_file_hints;
        FILE *f=NULL;
        struct vol3d *v=NULL;

        if (NULL==hints) { hints = &myhints; } /* Use default hints */

        hints->filename = filename;
        if (NULL==(f=fopen(filename,"rb"))) {
                perror("vol3df_new_from_file_heuristic: fopen");
                return NULL;
        }

        if (NULL==(v=vol3df_new_from_fh_heuristic(f,hints))) {
                fprintf(stderr,"vol3df_new_from_file_heuristic: "
                               "vol3df_new_from_fh_heuristic failed\n");
        }
        fclose(f);
        return v;
}

struct vol3d *vol3df_new_from_xdrfh (
                unsigned int nx, unsigned int ny, unsigned int nz,FILE *f)
{
        struct vol3d *v=NULL;
        XDR xdrs;

        if (NULL==(v=vol3df_new(nx,ny,nz))) {
                fprintf(stderr,"vol3df_new_from_xdrfh: vol3df_new failed\n");
                return NULL;
        }

        xdrstdio_create(&xdrs,f,XDR_DECODE);

        vol3d_stream(v,vol3df_xdr_streamfunc,&xdrs);

        xdr_destroy(&xdrs);

        return v;

}

int vol3df_write_xdrfh (struct vol3d *v,FILE *f)
{
        XDR xdrs;

        xdrstdio_create(&xdrs,f,XDR_ENCODE);

        vol3d_stream(v,vol3df_xdr_streamfunc,&xdrs);

        xdr_destroy(&xdrs);

        return 0;

}

struct ucdata {
        float max;
        float min;
        FILE *f;
};

void vol3df_uc_streamfunc(void *el,void *data)
{
        unsigned char u;
        float phi;
        struct ucdata *ucdata;

        phi = *(float *)el;
        ucdata = (struct ucdata *)data;

        phi = 255.0 * ( phi - ucdata->min )/(ucdata->max - ucdata->min);
        u = (unsigned char)floor(phi);
        fwrite(&u,1,1,ucdata->f);

}

int vol3df_write_ucfh (struct vol3d *v,FILE *f)
{
        struct vol3df_stats *stats=NULL;
        struct ucdata u;

        if (NULL==(stats=vol3df_getstats(v))) {
                fprintf(stderr,"vol3df_write_ucfh: vol3df_stats() failed\n");
                return -1;
        }
        u.max = stats->max; u.min = stats->min; u.f = f;

        vol3d_stream(v,vol3df_uc_streamfunc,&u);

        free(stats);

        return 0;

}

struct vol3d *vol3df_new_from_isofh (FILE *f) {
        struct vol3d *v=NULL;
        int nx,ny,nz;

        if (1!=fread(&nx,sizeof(int),1,f)) { perror("fread"); return NULL; }
        if (1!=fread(&ny,sizeof(int),1,f)) { perror("fread"); return NULL; }
        if (1!=fread(&nz,sizeof(int),1,f)) { perror("fread"); return NULL; }

        if (NULL==(v=vol3df_new(nx,ny,nz))) {
                fprintf(stderr,"vol3df_new_fromxdrstream: vol3df_new failed\n");
                return NULL;
        }

        if (0!=vol3d_read_raw_fhandle(v,f)) {
                fprintf(stderr,
                           "vol3df_new_from_isofh: vol3d_read_raw_fh failed\n");
                vol3d_destroy(v);
                return NULL;
        }

        return v;
}

struct vol3d *vol3df_new_from_rawfh (
                unsigned int nx, unsigned int ny, unsigned int nz,FILE *f)
{
        struct vol3d *v=NULL;

        if (NULL==(v=vol3df_new(nx,ny,nz))) {
                fprintf(stderr,"vol3df_new_fromxdrstream: vol3df_new failed\n");
                return NULL;
        }

        if (0!=vol3d_read_raw_fhandle(v,f)) {
                fprintf(stderr,
                           "vol3df_new_from_rawfh: vol3d_read_raw_fh failed\n");
                vol3d_destroy(v);
                return NULL;
        }

        return v;
}

int vol3df_write_isofh(struct vol3d *v, FILE *f)
{
        int nx,ny,nz;
        nx=v->nx;
        ny=v->ny;
        nz=v->nz;

        fwrite(&nx,sizeof(int),1,f);
        fwrite(&ny,sizeof(int),1,f);
        fwrite(&nz,sizeof(int),1,f);
        return vol3d_write_raw_fhandle(v,f);
}

int vol3df_write_fh(struct vol3d *v, FILE *f)
{
        return vol3d_write_raw_fhandle(v,f);
}


#ifdef HAVE_HDF5
/* Take a location and pathname for a dataset in an HDF5 file;
 * return 1 if the corresponding object is a 3d dataset.
 */

static int is_3d_dataset(hid_t loc_id,char *name)
{
        hid_t dset_id,dspace_id;
        int ndims;

        dset_id = H5Dopen(loc_id,name);
        dspace_id = H5Dget_space(dset_id);

        ndims = H5Sget_simple_extent_ndims(dspace_id);

        H5Sclose(dspace_id);
        H5Dclose(dset_id);

        if (3==ndims) { return 1; } else { return 0; }

}

/* This function is passed to H5Giterate() as a callback. H5Giterate
 * calls it once for every object in the "/" tree. If it returns zero,
 * H5Giterate will continute and call it for the next object; if it
 * returns positive or negative, then H5Giterate will return.
 *
 * Returns positive if a 3d voxel dataset is found, negative on error.
 */

static herr_t find_3d_callback (hid_t loc_id, const char *name, void *data)
{

        H5G_stat_t statbuf;
        char *path = NULL; /* Pointer to my copy of path string */

        /* Find out the object type */

        H5Gget_objinfo(loc_id,
                        name,
                        1,   /* follow links */
                        &statbuf);

        if (H5G_DATASET == statbuf.type) {
                /* If it's a 3d dataset, make a copy of the pathname
                 * and pass it back to the caller. */

                if (is_3d_dataset(loc_id,(char *)name)) {
                        size_t pathlen=0;
                        pathlen = strlen(name);
                        if (NULL==(path=malloc(1 + pathlen)))
                                { perror("malloc"); return -1; }
                        strncpy(path,name,pathlen+1);
                        *( (char **) data) = path;
                        return 1;
                }
        }

        return 0; /* Continue */

}

/* Take the ID of an HDF5 file; try to find a 3d dataset in the root
 * group. If found, return a (newly malloced) pointer to the pathname;
 * return NULL otherwise.
 */
static char *find_3d_dataset(hid_t file_id)
{
        char *pathname=NULL;
        if (1==H5Giterate(file_id,
                        "/",
                        NULL, /* Pointer to start integer */
                        find_3d_callback, /* Callback function */
                        &pathname /* Callback data */
                  )) {
                return pathname;
        } else { return NULL; }
}
#endif /* HAVE_HDF5 */

struct vol3d *vol3df_new_fromhdf5(char *fname,char *path)
#ifdef HAVE_HDF5
{
        hid_t file_id,dset_id,dspace_id,type_id;
        int ndims;
        hsize_t dims[3];
        struct vol3d *v=NULL;
        herr_t retval;

        if (NULL==fname) {
                fprintf(stderr,"vol3df_new_fromhdf5: NULL filename\n");
                return NULL;
        }

        if (0>(file_id=H5Fopen(fname,H5F_ACC_RDONLY,H5P_DEFAULT))) {
                fprintf(stderr,"vol3df_new_fromhdf5: H5Fopen failed\n");
                return NULL;
        }

        if (NULL==path) {
                /* No path given; see if we can find a dataset */
                if (NULL==(path=find_3d_dataset(file_id))) {
                        fprintf(stderr,
                                "vol3df_new_fromhdf5: no 3d dataset found\n");
                        return NULL;
                }
                /* Maybe it should be possible to suppress this message
                 * with a verbosity flag..
                 */

                fprintf(stderr,"Found dataset \"%s\"\n",path);
        }

        if (0>(dset_id=H5Dopen(file_id,path))) {
                fprintf(stderr,"vol3df_new_fromhdf5: H5Dopen failed\n");
                return NULL;
        }

        if (0>(dspace_id=H5Dget_space(dset_id))) {
                fprintf(stderr,"vol3df_new_fromhdf5: H5Dget_space failed\n");
                return NULL;
        }

        if (0>(type_id = H5Dget_type(dset_id))) {
                fprintf(stderr,"vol3df_new_fromhdf5: H5Dget_type failed\n");
                return NULL;
        }


        if (3!=(ndims=H5Sget_simple_extent_ndims(dspace_id))) {
                fprintf(stderr,"vol3df_new_fromhdf5: ndims=%d, !=3\n",
                                ndims);
                return NULL;
        }

        H5Sget_simple_extent_dims(dspace_id,dims,NULL);

        if (NULL==(v=vol3df_new(
                                (unsigned int)dims[0],
                                (unsigned int)dims[1],
                                (unsigned int)dims[2]
                             ))) {
                fprintf(stderr,"vol3df_new_fromhdf5: vol3df_new failed\n");
                return NULL;
        }

        if (0>(retval = H5Dread(dset_id,
                        H5T_NATIVE_FLOAT,      /* mem_type_id */
                        H5S_ALL,
                        H5S_ALL,
                        H5P_DEFAULT,            /* xfer_plist_id */
                        v->data
                        ))) {
                fprintf(stderr,"H5Dread() failed\n");
                vol3d_destroy(v);
                return NULL;
        }

        H5Sclose(dspace_id);
        H5Dclose(dset_id);
        H5Fclose(file_id);

        return v;

}
#else
{
        vol3df_unsupported_mess("hdf5");
        return NULL;
}
#endif /* HAVE_HDF5 */

int vol3df_write_hdf5(struct vol3d *v, char *fname, char *path)
#ifdef HAVE_HDF5
{
        hid_t file_id,dset_id,dspace_id;
        hsize_t dims[3];

        if (0>(file_id=H5Fcreate(fname,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT))){
                fprintf(stderr,"vol3df_write_hdf5: H5Fcreate failed\n");
                return -1;
        }

        dims[0]=v->nx;
        dims[1]=v->ny;
        dims[2]=v->nz;


        /* FIXME: this will probably fuck up if the path is nested;
         * it should create groups as required.
         */

        if (0>(dspace_id=H5Screate_simple(3,dims,NULL))) {
                fprintf(stderr,"vol3df_write_hdf5: H5Screate_simple failed\n");
                return -1;
        }

        if (0>(dset_id=H5Dcreate(file_id,path,H5T_IEEE_F32BE,
                                        dspace_id,H5P_DEFAULT))) {
                fprintf(stderr,"vol3df_write_hdf5: H5Dcreate failed\n");
                return -1;
        }

        if (0>H5Dwrite(dset_id,
                        H5T_NATIVE_FLOAT,      /* mem_type_id */
                        H5S_ALL,                /* mem_space_id */
                        dspace_id,              /* file_space_id */
                        H5P_DEFAULT,            /* xfer_plist_id */
                        v->data
                      )) {
                fprintf(stderr,"vol3df_write_hdf5: H5Dwrite failed\n");
                return -1;
        }

        H5Dclose(dset_id);
        H5Sclose(dspace_id);
        H5Fclose(file_id);

        return 0;

}
#else
{
        vol3df_unsupported_mess("hdf5");
        return -1;
}
#endif /* HAVE_HDF5 */

int vol3df_write_vtkfh(struct vol3d *v, FILE *f)
{
        fprintf(f,
          "# vtk DataFile Version 2.0\nGenerated by vol3df library\nBINARY\n");
        fprintf(f,"DATASET STRUCTURED_POINTS\n");
        fprintf(f,"DIMENSIONS %d %d %d\n",v->nx,v->ny,v->nz);
        fprintf(f,"ORIGIN 0 0 0\n");
        fprintf(f,"SPACING 1 1 1\n");
        fprintf(f,"POINT_DATA %d\n",v->nx*v->ny*v->nz);
        fprintf(f,"SCALARS scalars float 1\n");
        fprintf(f,"LOOKUP_TABLE default\n");

        return vol3df_write_xdrfh (v,f);
}

int vol3df_write_ucvtkfh(struct vol3d *v, FILE *f)
{
        fprintf(f,
          "# vtk DataFile Version 2.0\nGenerated by vol3df library\nBINARY\n");
        fprintf(f,"DATASET STRUCTURED_POINTS\n");
        fprintf(f,"DIMENSIONS %d %d %d\n",v->nx,v->ny,v->nz);
        fprintf(f,"ORIGIN 0 0 0\n");
        fprintf(f,"SPACING 1 1 1\n");
        fprintf(f,"POINT_DATA %d\n",v->nx*v->ny*v->nz);
        fprintf(f,"SCALARS scalars unsigned_char 1\n");
        fprintf(f,"LOOKUP_TABLE default\n");

        return vol3df_write_ucfh (v,f);
}

int vol3df_write_fnrrdfh(struct vol3d *v, FILE *f)
{
        fprintf(f,"NRRD0001\n");
        fprintf(f,"# Generated by vol3d\n");
        fprintf(f,"type: float\ndimension: 3\n");
        fprintf(f,"sizes: %d %d %d\n",v->nx,v->ny,v->nz);
        fprintf(f,"endian: big\nencoding: raw\n\n");

        return vol3df_write_xdrfh (v,f);
}

int vol3df_write_ucnrrdfh(struct vol3d *v, FILE *f)
{
        fprintf(f,"NRRD0001\n");
        fprintf(f,"# Generated by vol3d\n");
        fprintf(f,"type: unsigned char\ndimension: 3\n");
        fprintf(f,"sizes: %d %d %d\n",v->nx,v->ny,v->nz);
        fprintf(f,"encoding: raw\n\n");

        return vol3df_write_ucfh (v,f);
}


static void accum_streamfunc(void *el,void *data)
{

        struct vol3df_dist *dist;
        int i;
        float phi;

        phi = *(float *)el;
        dist = (struct vol3df_dist *)data;
        i=(int)floor( ( phi - dist->stats->min)/dist->dphi );

        if (i<0) {
                fprintf(stderr,"Warning: phi=%f < phimin=%f\n",
                        phi,dist->stats->min);
                i=0;
        }
        if (i>=dist->nbins) {
                fprintf(stderr,"Warning: phi=%f > phimax=%f\n",
                        phi,dist->stats->max);
                i=dist->nbins-1;
        }

        dist->count[i]++;
}

struct vol3df_dist *vol3df_getdist(struct vol3d *v, int nbins)
{

        struct vol3df_dist *dist=NULL;
        int i;


        if (NULL==(dist=malloc(sizeof(struct vol3df_dist)))) {
                perror("vol3df_getdist: malloc");
                return NULL;
        }

        dist->nbins = nbins;

        if (NULL==(dist->stats=vol3df_getstats(v))) {
                fprintf(stderr,"vol3df_dist: vol3df_getstats failed\n");
                return NULL;
        }

        /* Allocate bins */

        if (NULL==(dist->count=malloc(sizeof(unsigned int)*nbins))) {
                perror("malloc()");
                return NULL;
        }

        /* Zero all counts */

        for (i=0;i<dist->nbins;i++) {
                dist->count[i]=0;
        }

        dist->dphi = (dist->stats->max - dist->stats->min)/(float)(nbins-1);

        vol3d_stream(v,accum_streamfunc,dist);

        return dist;

}

struct vol3df_stats *vol3df_getstats(struct vol3d *v)
{
        struct vol3df_stats *stat=NULL;
        unsigned long i;
        float phi;

        if (VOL_TYPE_FLOAT!=v->datatype) {
                fprintf(stderr,"vol3df_stats: datatype is not float\n");
                return NULL;
        }

        if (NULL==(stat=malloc(sizeof(struct vol3df_stats)))) {
                perror("vol3df_stats: malloc");
                return NULL;
        }

        stat->max = stat->min = stat->mean = phi = *(float *)v->data;
        stat->variance = phi*phi;

        for (i=1;i<v->nx*v->ny*v->nz;i++) {
                phi = *( (float *)v->data + i);
                stat->mean += phi;
                stat->variance   += phi*phi;
                if (stat->min>phi) { stat->min=phi; }
                if (stat->max<phi) { stat->max=phi; }
        }

        stat->mean /= (float) (v->nx*v->ny*v->nz);
        stat->variance /= (float) (v->nx*v->ny*v->nz);
        stat->variance -= stat->mean*stat->mean;


        return stat;
}

int vol3df_normalise(struct vol3d *v,float amplitude)
{
    struct vol3df_stats *stats=NULL;
    unsigned long n,i;
    float *p=NULL;
    if (NULL==(stats=vol3df_getstats(v))) {
        fprintf(stderr,"vol3df_normalise: vol3df_getstats failed\n");
        return -1;
    }

    n=v->nx*v->ny*v->nz;

    p=(float *)v->data;

    for (i=0;i<n;i++) {
        p[i] = amplitude*(2.0*(p[i]-stats->min)/(stats->max-stats->min)-1.0);
    }

    return 0;

}

void abs_streamfunc(void *el, void *data)
{
        *(float *)el=fabs(*(float *)el);
        
}

/* This is where I start wishing that C had closures */

void gt_streamfunc(void *el, void *data)
{
        float phi,thresh;
        phi = *(float *)el; thresh = *(float *) data;
        if (phi > thresh) { *(float *)el = 1.0; } else { *(float *)el = 0.0; }
}

void lt_streamfunc(void *el, void *data)
{
        float phi,thresh;
        phi = *(float *)el; thresh = *(float *) data;
        if (phi < thresh) { *(float *)el = 1.0; } else { *(float *)el = 0.0; }
}

int vol3df_abs(struct vol3d *v)
{
        if (VOL_TYPE_FLOAT!=v->datatype) {
                fprintf(stderr,"vol3df_abs: datatype is not float\n");
                return -1;
        }

        vol3d_stream(v,abs_streamfunc,NULL);
        return 0;
}

void range_streamfunc(void *el, void *data)
{
        float phi,pmax,pmin;

        phi = *(float *)el;
        pmax = *(float *)data;
        pmin = *( (float *)data + 1);
        
        if ( (phi > pmin) && (phi < pmax) ) {
                *(float *)el = 1.0;
        } else {
                *(float *)el = 0.0;
        }
}

void rangefilter_streamfunc(void *el, void *data)
{
        float phi,pmax,pmin;

        phi = *(float *)el;
        pmax = *(float *)data;
        pmin = *( (float *)data + 1);
        
        if ( (phi > pmin) && (phi < pmax) ) {
                /* leave untouched */
        } else {
                *(float *)el = 0.0;
        }
}

int vol3df_range(struct vol3d *v, float pmax, float pmin)
{
        float p[2];

        if (VOL_TYPE_FLOAT!=v->datatype) {
                fprintf(stderr,"vol3df_range: datatype is not float\n");
                return -1;
        }

        p[0]=pmax; p[1]=pmin;

        vol3d_stream(v,range_streamfunc,p);
        return 0;
}

int vol3df_rangefilter(struct vol3d *v, float pmax, float pmin)
{
        float p[2];

        if (VOL_TYPE_FLOAT!=v->datatype) {
                fprintf(stderr,"vol3df_rangefilter: datatype is not float\n");
                return -1;
        }

        p[0]=pmax; p[1]=pmin;

        vol3d_stream(v,rangefilter_streamfunc,p);
        return 0;
}

int vol3df_gt(struct vol3d *v, float thresh)
{
        if (VOL_TYPE_FLOAT!=v->datatype) {
                fprintf(stderr,"vol3df_gt: datatype is not float\n");
                return -1;
        }

        vol3d_stream(v,gt_streamfunc,&thresh);
        return 0;
}

int vol3df_lt(struct vol3d *v, float thresh)
{
        if (VOL_TYPE_FLOAT!=v->datatype) {
                fprintf(stderr,"vol3df_lt: datatype is not float\n");
                return -1;
        }

        vol3d_stream(v,lt_streamfunc,&thresh);
        return 0;
}

int wrap(int x, int nx)
{
        if (x<0) { return x+nx; }
        if (x>=nx) { return x-nx; }
        return x;
}

float vol3df_wrap(struct vol3d *v, int x, int y, int z)
{ return *vol3df_wrap_ptr(v,x,y,z); }

float *vol3df_wrap_ptr(struct vol3d *v, int x, int y, int z)
{
        while (x<0) { x += v->nx; }
        while (y<0) { y += v->ny; }
        while (z<0) { z += v->nz; }

        while (x >= v->nx) { x -= v->nx; }
        while (y >= v->ny) { y -= v->ny; }
        while (z >= v->nz) { z -= v->nz; }

        return (
                        (float *)v->data
                        + z*v->nx*v->ny
                        + y*v->nx
                        + x
                );

}

static float vmin(float x, float y)
{
        return (x>y) ? y : x;
}

static float vmax(float x, float y)
{
        return (x>y) ? x : y;
}

/*
 * Hey now
 * Hey now now
 * Sing this erosion to me..
 */

struct vol3d *vol3df_new_erosion(struct vol3d *v)
{
        int x,y,z;
        int dx,dy,dz;
        float phi;
        struct vol3d *vnew=NULL;

        if (VOL_TYPE_FLOAT!=v->datatype) {
                fprintf(stderr,"vol3df_new_erosion: datatype is not float\n");
                return NULL;
        }

        if (NULL==(vnew=vol3d_new_copytype(v))) {
                fprintf(stderr,"vol3df_new_erosion: vol3d_copytype failed\n");
                return NULL;
        }

        for (z=0;z<v->nz;z++) {
        for (y=0;y<v->ny;y++) {
        for (x=0;x<v->nx;x++) {
                phi = *(float *)vol3d_element(v,x,y,z);

                for (dz=-1;dz<=1;dz++) {
                for (dy=-1;dy<=1;dy++) {
                for (dx=-1;dx<=1;dx++) {
                        phi = vmin(phi,vol3df_wrap(v,x+dx,y+dy,z+dz)); 
                }}}

                *(float *)vol3d_element(vnew,x,y,z) = phi;
        }}}

        return vnew;
}

struct vol3d *vol3df_new_dilation(struct vol3d *v)
{
        int x,y,z;
        int dx,dy,dz;
        float phi;
        struct vol3d *vnew=NULL;

        if (VOL_TYPE_FLOAT!=v->datatype) {
                fprintf(stderr,"vol3df_new_dilation: datatype is not float\n");
                return NULL;
        }

        if (NULL==(vnew=vol3d_new_copytype(v))) {
                fprintf(stderr,"vol3df_new_dilation: vol3d_copytype failed\n");
                return NULL;
        }

        for (z=0;z<v->nz;z++) {
        for (y=0;y<v->ny;y++) {
        for (x=0;x<v->nx;x++) {
                phi = *(float *)vol3d_element(v,x,y,z);

                for (dz=-1;dz<=1;dz++) {
                for (dy=-1;dy<=1;dy++) {
                for (dx=-1;dx<=1;dx++) {
                        phi = vmax(phi,vol3df_wrap(v,x+dx,y+dy,z+dz));
                }}}

                *(float *)vol3d_element(vnew,x,y,z) = phi;
        }}}

        return vnew;
}


/* Return the Euler number for the given volume */

#define VOL3D_EULER_SURFACE 0
#define VOL3D_EULER_VOLUME  1


static float vol3df_totaleuler(struct vol3d *v,
                int (*testfunc)(float, void * ),
                void *data, int eulertype)
{
        int x,y,z;
        long int sum=0; /* 8x Euler number is summed here */
        if (
                (eulertype!=VOL3D_EULER_SURFACE) && 
                (eulertype!=VOL3D_EULER_VOLUME) )  {
                fprintf(stderr,"vol3df_totaleuler: bad eulertype %d\n",
                                eulertype);
                return -1;
        }

        if (VOL_TYPE_FLOAT!=v->datatype) {
                fprintf(stderr,"vol3df_euler: datatype is not float\n");
                return -1;
        }

        for (z=0;z<v->nz;z++) {
        for (y=0;y<v->ny;y++) {
        for (x=0;x<v->nx;x++) {
                int mask=0; /* 1 bit set for each selected voxel */
                int bit=1;
                int i;

                for (i=0;i<8;i++) {
                        int xx,yy,zz;
                        xx = x + euler_delta[i][0];
                        yy = y + euler_delta[i][1];
                        zz = z + euler_delta[i][2];

                        if (testfunc(vol3df_wrap(v,xx,yy,zz),data)) {
                                mask |= bit;
                        }

                        bit <<=1;
                }

                /* "mask" now has 1 bit set for each corresponding
                 * selected voxel. Look up and sum 8x the Euler number.
                 */

                if (VOL3D_EULER_SURFACE==eulertype) {
                        sum += surf_euler_lookup[mask];
                } else {
                        sum += vol_euler_lookup[mask];
                }

        }}}

        return ((float)sum)/8.0;
}

float vol3df_total_surf_euler(struct vol3d *v,
                int (*testfunc)(float, void * ),
                void *data)
{
        return vol3df_totaleuler(v,testfunc,data,VOL3D_EULER_SURFACE);
}

float vol3df_total_vol_euler(struct vol3d *v,
                int (*testfunc)(float, void * ),
                void *data)
{
        return vol3df_totaleuler(v,testfunc,data,VOL3D_EULER_VOLUME);
}

/* At each voxel, calculate the Euler number of the 1x1x1
 * cube centred at the (1,1,1) vertex of the voxel, calculated from the
 * local 2x2x2 cluster. Store the Euler number in the new volume, which
 * is returned. Return NULL on error.
 */


static struct vol3d *vol3df_new_euler(struct vol3d *v,
                int (*testfunc)(float, void * ),
                void *data, int eulertype)
{
        int x,y,z;
        struct vol3d *vnew=NULL;

        if (
                (eulertype!=VOL3D_EULER_SURFACE) && 
                (eulertype!=VOL3D_EULER_VOLUME) )  {
                fprintf(stderr,"vol3df_new_euler: bad eulertype %d\n",
                                eulertype);
                return NULL;
        }

        if (VOL_TYPE_FLOAT!=v->datatype) {
                fprintf(stderr,"vol3df_euler: datatype is not float\n");
                return NULL;
        }

        if (NULL==(vnew=vol3d_new_copytype(v))) {
                fprintf(stderr,"vol3df_new_euler: vol3d_copytype failed\n");
                return NULL;
        }

        for (z=0;z<v->nz;z++) {
        for (y=0;y<v->ny;y++) {
        for (x=0;x<v->nx;x++) {
                int mask=0; /* 1 bit set for each selected voxel */
                int bit=1;
                int i;

                for (i=0;i<8;i++) {
                        int xx,yy,zz;
                        xx = x + euler_delta[i][0];
                        yy = y + euler_delta[i][1];
                        zz = z + euler_delta[i][2];

                        if (testfunc(vol3df_wrap(v,xx,yy,zz),data)) {
                                mask |= bit;
                        }

                        bit <<=1;
                }

                /* "mask" now has 1 bit set for each corresponding
                 * selected voxel. Look up and store the Euler number.
                 */


                if (VOL3D_EULER_SURFACE==eulertype) {
                        *(float *)vol3d_element(vnew,x,y,z)
                                = surf_euler_lookup[mask]/8.0;
                } else {
                        *(float *)vol3d_element(vnew,x,y,z)
                                = vol_euler_lookup[mask]/8.0;
                }

        }}}

        return vnew;
}

struct vol3d *vol3df_new_surf_euler(struct vol3d *v,
                int (*testfunc)(float, void * ),
                void *data) {
        return vol3df_new_euler(v,testfunc,data,VOL3D_EULER_SURFACE);
}
struct vol3d *vol3df_new_vol_euler(struct vol3d *v,
                int (*testfunc)(float, void * ),
                void *data) {
        return vol3df_new_euler(v,testfunc,data,VOL3D_EULER_VOLUME);
}


static struct vol3d *vol3df_new_euler8(struct vol3d *v, float threshold,int eulertype)
{
        struct vol3d *vnew=NULL;
        int vnx,vny,vnz;
        int bx,by,bz; /* block x,y,z: refers to element of vnew */
        const int bsize=8; /* Block size */

        if (
                (eulertype!=VOL3D_EULER_SURFACE) && 
                (eulertype!=VOL3D_EULER_VOLUME) )  {
                fprintf(stderr,"vol3df_new_euler8: bad eulertype %d\n",
                                eulertype);
                return NULL;
        }

        if ((0!=(v->nx % bsize))||(0!=(v->ny % bsize))||(0!=(v->nz % bsize)) ) {
                fprintf(stderr,
                    "vol3df_new_euler8: dimensions must be multiples of %d\n",
                    bsize);
                return NULL;
        }

        if (VOL_TYPE_FLOAT!=v->datatype) {
                fprintf(stderr,"vol3df_new_euler8: datatype is not float\n");
                return NULL;
        }

        vnx=v->nx/bsize; vny=v->ny/bsize; vnz=v->nz/bsize;

        if (NULL==(vnew=vol3df_new(vnx,vny,vnz))) {
                fprintf(stderr,"vol3df_new_euler8: vol3df_new() failed\n");
                return NULL;
        }

        for (bz=0;bz<vnz;bz++) {
        for (by=0;by<vny;by++) {
        for (bx=0;bx<vnx;bx++) {
                long int sum=0; /* sum 8*chi */
                int dx,dy,dz;

                for (dz=0;dz<bsize;dz++) {
                for (dy=0;dy<bsize;dy++) {
                for (dx=0;dx<bsize;dx++) {
                        int mask=0; /* 1 bit set for each selected voxel */
                        int bit=1;
                        int i;

                        int x = bx*bsize + dx;
                        int y = by*bsize + dy;
                        int z = bz*bsize + dz;

                        for (i=0;i<8;i++) {
                                int xx,yy,zz;
                                xx = x + euler_delta[i][0];
                                yy = y + euler_delta[i][1];
                                zz = z + euler_delta[i][2];

                                if (vol3df_wrap(v,xx,yy,zz)>threshold) {
                                        mask |= bit;
                                }

                                bit <<=1;
                        }

                        /* "mask" now has 1 bit set for each corresponding
                         * selected voxel. Look up and sum 8x the Euler number.
                         */

                        if (VOL3D_EULER_SURFACE==eulertype) {
                                sum += surf_euler_lookup[mask];
                        } else {
                                sum += vol_euler_lookup[mask];
                        }

                }}} /* dx,dy,dz */

                * (float *)vol3d_element(vnew,bx,by,bz) = ((float)sum)/8.0;

        }}} /* bx,by,bz */

        return vnew;


}

struct vol3d *vol3df_new_surf_euler8(struct vol3d *v, float threshold)
{ return vol3df_new_euler8(v,threshold,VOL3D_EULER_SURFACE); }

struct vol3d *vol3df_new_vol_euler8(struct vol3d *v, float threshold)
{ return vol3df_new_euler8(v,threshold,VOL3D_EULER_VOLUME); }

/* Go through each site, and apply the given boolean-valued function
 * to it. If true, then replace it by the number of nearest neighbours.
 * If false, set to zero.
 */

struct vol3d *vol3df_new_nncount(struct vol3d *v,
                int (*testfunc)(float, void * ),
                void *data)
{
        int x,y,z;
        float phi;
        struct vol3d *vnew=NULL;

        if (VOL_TYPE_FLOAT!=v->datatype) {
                fprintf(stderr,"vol3df_new_dilation: datatype is not float\n");
                return NULL;
        }

        if (NULL==(vnew=vol3d_new_copytype(v))) {
                fprintf(stderr,"vol3df_new_nncount: vol3d_copytype failed\n");
                return NULL;
        }

        for (z=0;z<v->nz;z++) {
        for (y=0;y<v->ny;y++) {
        for (x=0;x<v->nx;x++) {
                float phinew=0.0;
                int dx,dy,dz;

                phi = *(float *)vol3d_element(v,x,y,z);

                if (testfunc(phi,data)) {

                        for (dz=-1;dz<=1;dz++) {
                        for (dy=-1;dy<=1;dy++) {
                        for (dx=-1;dx<=1;dx++) {
                                phinew += testfunc(
                                        vol3df_wrap(v,x+dx,y+dy,z+dz),data
                                               );
                        }}}
                        phinew--; /* Don't count yourself */

                }

                *(float *)vol3d_element(vnew,x,y,z) = phinew;
        }}}

        return vnew;
}


/* Run through the lattice. At each voxel (x,y,z), sum up the neighbours
 * in the range (x,y,z) to (x+d-1,y+d-1,z+d-1). Replace the voxel with the
 * mean of its neighbourhood.
 */

struct vol3d *vol3df_new_smear(struct vol3d *v, int delta)
{
        int x,y,z;
        int nx,ny,nz;
        struct vol3d *vnew=NULL;
        float stencilsize; /* No of voxels we average over */

        nx=v->nx; ny=v->ny; nz=v->nz; 

        if ( (delta>=(nx/2)) || (delta>=(ny/2)) || (delta>=(nz/2)) ) {
                fprintf(stderr,
                "vol3df_new_smear: too large smearsize %d for %dx%dx%d\n",
                delta,nx,ny,nz); return NULL;
        }

        stencilsize = delta*delta*delta;

        if (NULL==(vnew=vol3d_new_copytype(v))) {
                fprintf(stderr,"vol3df_new_smear: vol3d_copytype failed\n");
                return NULL;
        }

        for (z=0;z<v->nz;z++) {
        for (y=0;y<v->ny;y++) {
        for (x=0;x<v->nx;x++) {
                int dx,dy,dz;
                float phi=0.0;

                        for (dz=0;dz<delta;dz++) {
                        for (dy=0;dy<delta;dy++) {
                        for (dx=0;dx<delta;dx++) {
                                phi += vol3df_wrap(v,x+dx,y+dy,z+dz);
                        }}}

                *(float *)vol3d_element(vnew,x,y,z) = phi/stencilsize;
        }}}

        return vnew;

}

/* Take a vol3d structure. Return a new one, in which each voxel
 * contains the value of the voxel of the previous volume, at the same
 * position plus a given vector. Wrap periodically.
 */
struct vol3d *vol3df_new_roll(struct vol3d *v, int dx, int dy, int dz)
{
        int x,y,z;
        struct vol3d *vnew=NULL;

        if (VOL_TYPE_FLOAT!=v->datatype) {
                fprintf(stderr,"vol3df_euler: datatype is not float\n");
                return NULL;
        }

        if (NULL==(vnew=vol3d_new_copytype(v))) {
                fprintf(stderr,"vol3df_new_smear: vol3d_copytype failed\n");
                return NULL;
        }

        for (z=0;z<v->nz;z++) {
        for (y=0;y<v->ny;y++) {
        for (x=0;x<v->nx;x++) {
                *(float *)vol3d_element(vnew,x,y,z) =
                        vol3df_wrap(v,x-dx,y-dy,z-dz);
        }}}

        return vnew;

}

/* Take a vol3df, and a pointer to a function. The function takes three
 * coordinate arguments, and a void pointer. The function is invoked
 * once for each voxel; the resulting float is stored at that voxel.
 * The "data" argument is passed straight through to the function.
 */

void vol3df_fill_func(struct vol3d *v,
                float (*func)(int , int , int , void *), void *data)
{
        int x,y,z;
        for (z=0;z<v->nz;z++) {
        for (y=0;y<v->ny;y++) {
        for (x=0;x<v->nx;x++) {
                *(float *)vol3d_element(v,x,y,z) = func(x,y,z,data);
        }}}

}

/* Take two vol3d arguments, a and b. a and b must have the same
 * dimensions. Create a new vol3d object, created by applying the
 * supplied function to elements of a and b.
 *
 * FIXME: should be generalized and moved into vol3d.c
 */

struct vol3d *vol3df_new_pointwise_binaryop( struct vol3d *a, struct vol3d *b,
                float (*func)(float, float))
{
        struct vol3d *c=NULL;
        int i,nels;
        float *pa,*pb,*pc;

        if ((VOL_TYPE_FLOAT!=a->datatype) || (VOL_TYPE_FLOAT!=b->datatype)) {
                fprintf(stderr,
                    "vol3df_new_pointwise_binaryop: datatype is not float\n");
                return NULL;
        }
        if ( (a->nx != b->nx) || (a->ny != b->ny) || (a->nz != b->nz)) {
                fprintf(stderr,
                        "vol3df_new_pointwise_binaryop: "
                        "input volumes must have same dimensions\n");
                return NULL;
        }

        if (NULL==(c=vol3d_new_copytype(a))) {
                fprintf(stderr,
                        "vol3df_new_pointwise_binaryop: "
                        "vol3d_new_copytype failed\n");
                return NULL;
        }

        pa=a->data; pb=b->data; pc=c->data;
        nels = a->nx * a->ny * a->nz;

        for (i=0;i<nels;i++) {
                *pc++ = func(*pa++,*pb++);
        }

        return c;

}

struct vol3d *vol3df_new_rfft(struct vol3d *v)
#ifdef HAVE_FFTW
{
        fftwf_plan plan;
        struct vol3d *vnew=NULL;
        int newx,newy,newz;

        /* Dimensions of new array are slightly larger */
        newx=v->nx+2; newy=v->ny; newz = v->nz;

        if (NULL==(vnew=vol3d_new(newx,newy,newz,VOL_TYPE_FLOAT)))
        { fprintf(stderr,"vol3df_new_fft: vol3d_new failed\n"); return NULL;}

        /* FFTW expects a C-order array whereas we use Fortran order,
         * so reverse X, y, and Z.
         */

        plan = fftwf_plan_dft_r2c_3d (v->nz,v->ny,v->nx,
                                (float *)v->data, 
                                (fftwf_complex *) vnew->data,
                                FFTW_ESTIMATE);

        fftwf_execute(plan);
        fftwf_destroy_plan(plan);

        return vnew;
        
}
#else
{
        vol3df_unsupported_mess("fftw");
        return NULL;
}
#endif /* HAVE_FFTW */

struct vol3d *vol3df_new_spectrum(struct vol3d *v)
#ifdef HAVE_FFTW
{
        fftwf_plan plan;
        struct vol3d *v2=NULL,*vnew=NULL;
        int x,y,z,newx,newy,newz,px,py,pz;

        /* Dimensions of new array are slightly larger */
        newx=2*(1+(v->nx/2)); newy=v->ny; newz = v->nz;
        px=newx/2; py=newy; pz=newz;

        if (NULL==(v2=vol3d_new(newx,newy,newz,VOL_TYPE_FLOAT)))
        { fprintf(stderr,"vol3df_new_spectrum: vol3d_new failed\n"); return NULL;}
        if (NULL==(vnew=vol3d_new(px,py,pz,VOL_TYPE_FLOAT)))
        { fprintf(stderr,"vol3df_new_spectrum: vol3d_new failed\n"); return NULL;}

        plan = fftwf_plan_dft_r2c_3d(v->nz,v->ny,v->nx,
                                (float *)v->data, 
                                (fftwf_complex *) v2->data,
                                FFTW_ESTIMATE);

        fftwf_execute(plan);
        fftwf_destroy_plan(plan);

        for (z=0;z<pz;z++) {
        for (y=0;y<py;y++) {
        for (x=0;x<px;x++) {
                float *p;
                float re,im;

                p = (float *)v2->data
                        + 2*( x + px*y + px*py*z );
                re = *p; im = *(p+1);
                *(float *)vol3d_element(vnew,x,y,z) = re*re+im*im;
                
        }}}

        vol3d_destroy(v2);


        return vnew;
        
}
#else
{
        vol3df_unsupported_mess("fftw");
        return NULL;
}
#endif /* HAVE_FFTW */

int is_mirrored_3d(struct vol3d *v, int x)
{
        if (0==(v->nx%2)) {
                if ( (0==x)||(x==(v->nx/2))) { return 0; }
                else { return 1; }
        } else {
                if (0==x) { return 0; } else { return 1; }
        }
}

struct binnery *vol3df_sfactor(struct vol3d *v)
#ifdef HAVE_FFTW
{
        struct vol3d *transform=NULL;
        int n,nbins;
        struct binnery *binnery=NULL;
        int i,j,k,r;
        unsigned int tnx,tny,tnz;

        if (VOL_TYPE_FLOAT != v->datatype)
        { fprintf(stderr,"vol3df_sfactor: not a float dataset\n");return NULL;}


        if (NULL==(transform=vol3df_new_rfft(v)))
        { fprintf(stderr,"vol3df_sfactor: new_rfft failed\n"); return NULL;}

        n= v->nx * v->ny * v->nz;
        tnx=transform->nx/2; tny=transform->ny; tnz=transform->nz;

        /* Set nbins equal to half the minimum of (nx,ny,nz) */

        if (v->nx<v->ny) { nbins=v->nx; } else { nbins=v->ny; }
        if (nbins>v->nz) { nbins=v->nz; }
        nbins /= 2;

        if (NULL==(binnery=binnery_new(nbins)))
        { fprintf(stderr,"vol2df_sfactor: binnery_new failed\n");return NULL;}

        for (k=0;k<tnz;k++) {
        for (j=0;j<tny;j++) {
        for (i=0;i<tnx;i++) {
                int ii,jj,kk;
                float re,im,s;

                re = *( (float *)transform->data + 2*(i + j*tnx + k*tnx*tny)  );
                im = *( (float *)transform->data + 2*(i + j*tnx + k*tnx*tny)+1);
                s = (re*re+im*im);

                ii=i;
                if (j>(v->ny/2)) { jj= j - v->ny; } else { jj=j; }
                if (k>(v->nz/2)) { kk= k - v->nz; } else { kk=k; }
                /* ii>0; jj,kk may be pos or neg; (ii,jj,kk) corresponds to the
                 * unwrapped frequency vector. See notes of 19/5/2004 in
                 * lab book.
                 */

                r = (int)floor(sqrt(ii*ii+jj*jj+kk*kk));
                binnery_add_point(binnery,s,r);

                /* If this point is the complex conjugate of another,
                 * redundant, point, then sum the conjugate point
                 * as well.
                 */
                if (is_mirrored_3d(v,i)) {
                        ii = v->nx - i;
                        jj = v->ny - j;
                        if (jj > (v->ny/2)) { jj -= v->ny; }
                        if (kk > (v->nz/2)) { kk -= v->nz; }
                        r = (int)floor(sqrt(ii*ii+jj*jj+kk*kk));
                        binnery_add_point(binnery,s,r);
                }
        }}}


        vol3d_destroy(transform);


        return binnery;
        
}
#else
{
        vol3df_unsupported_mess("fftw");
        return NULL;
}
#endif /* HAVE_FFTW */

#ifdef HAVE_FFTW
int vol3df_inverse_f_noise(struct vol3d *v)
{
        struct vol3d *transform=NULL;
        int i,j,k,n;
        unsigned int tnx,tny,tnz;
        fftwf_plan plan;

        if (VOL_TYPE_FLOAT != v->datatype)
        { fprintf(stderr,"vol3df_sfactor: not a float dataset\n");return -1;}

        /* Fill with white noise in (-1,1) */

        vol3df_whitenoise(v,1.0,0.0);


        if (NULL==(transform=vol3df_new_rfft(v)))
        { fprintf(stderr,"vol3df_sfactor: new_rfft failed\n"); return -1;}

        n= v->nx * v->ny * v->nz;
        tnx=transform->nx/2; tny=transform->ny; tnz=transform->nz;

        for (k=0;k<tnz;k++) {
        for (j=0;j<tny;j++) {
        for (i=0;i<tnx;i++) {
                int ii,jj,kk;
                int mag2;
                float scalefactor;

                ii=i;
                if (j>(v->ny/2)) { jj= j - v->ny; } else { jj=j; }
                if (k>(v->nz/2)) { kk= k - v->nz; } else { kk=k; }
                /* ii>0; jj,kk may be pos or neg; (ii,jj,kk) corresponds to the
                 * unwrapped frequency vector. See notes of 19/5/2004 in
                 * lab book.
                 */

                mag2 = ii*ii+jj*jj+kk*kk;

                if (0!=mag2) {
                    scalefactor = 1.0/sqrt((float)mag2);

                    *( (float *)transform->data + 2*(i + j*tnx + k*tnx*tny)  )
                        *= scalefactor;
                    *( (float *)transform->data + 2*(i + j*tnx + k*tnx*tny)+1)
                        *= scalefactor;
                }

        }}}


        /* Now inverse transform to get inverse-f noise */

        plan = fftwf_plan_dft_c2r_3d(v->nz,v->ny, v->nx,
                                (fftwf_complex *) transform->data,
                                (float *)v->data, 
                                FFTW_ESTIMATE);

        fftwf_execute(plan);

        fftwf_destroy_plan(plan);

        vol3d_destroy(transform);

        /* Normalise it */
        vol3df_normalise(v,1.0);

        return 0;
        
}
#endif /* HAVE_FFTW */


static float rand11();
static float rand11()
{
	return (((float)rand()/RAND_MAX)-0.5)*2.0;
}

/* Return a random float between zeropoint-magnitude and zeropoint+magnitude */
static float whitenoise(float magnitude,float zeropoint) {
    return zeropoint + magnitude*rand11();

}

struct floatpair  { float x; float y; };

static void vol3df_whitenoise_streamfunc(void *el, void *data)
{
    struct floatpair *fp;
    fp = (struct floatpair *) data;

    *(float *)el = whitenoise(fp->x,fp->y);
}

/* Fill the volume with noise of a flat distribution in the range (0,1).
 */

int vol3df_whitenoise(struct vol3d *v,float magnitude,float zeropoint)
{
    struct floatpair fp;
    if (VOL_TYPE_FLOAT!=v->datatype) {
            fprintf(stderr,"vol3df_whitenoise: datatype is not float\n");
            return -1;
    }
    fp.x = magnitude;
    fp.y = zeropoint;

    vol3d_stream(v,vol3df_whitenoise_streamfunc,&fp);
    return 0;
}


struct vol3d *vol3df_autocorrelation(struct vol3d *v)
#ifdef HAVE_FFTW
{
        struct vol3d *transform=NULL;
        fftwf_plan plan;
        unsigned int i,n=0;
        float *p=NULL;
        unsigned int nrealpoints=0;


        nrealpoints = v->nx*v->ny*v->nz;
        if (NULL==(transform=vol3df_new(2*(1+(v->nx/2)),v->ny,v->nz))) {
                fprintf(stderr,
                          "vol3df_new_autocorrelation: vol3df_new failed.\n");
                return NULL;
        }

        plan = fftwf_plan_dft_r2c_3d(v->nz,v->ny,v->nx,
                                (float *)v->data, 
                                (fftwf_complex *) transform->data,
                                FFTW_ESTIMATE);

        fftwf_execute(plan);
        fftwf_destroy_plan(plan);

        /* Run through the transform array; replace each complex
         * element with its real absolute square.
         * Also, divide by the number of points, to preserve normalization
         * after the inverse transform.
         */

        p=(float *)transform->data;
        n=(transform->nx/2)*transform->ny*transform->nz;

        for (i=0;i<n;i++) {
                float re,im;

                re=*p;
                im=*(p+1);

                *p = (re*re+im*im)/nrealpoints;
                *(p+1)=0.0;

                p+=2;
        }

        plan = fftwf_plan_dft_c2r_3d(v->nz,v->ny, v->nx,
                                (fftwf_complex *) transform->data,
                                (float *)v->data, 
                                FFTW_ESTIMATE);

        fftwf_execute(plan);

        fftwf_destroy_plan(plan);

        return v;


}
#else
{
        vol3df_unsupported_mess("fftw");
        return NULL;
}
#endif


int vol3df_logabs(struct vol3d *v)
{
        int i,n;
        float *p=NULL;

        if (VOL_TYPE_FLOAT!=v->datatype) {
                fprintf(stderr,"vol3df_logabs: datatype is not float\n");
                return -1;
        }
        n=v->nx*v->ny*v->nz;
        p=(float *)v->data;

        for (i=0;i<n;i++) {
                float phi;
                phi = *p;
                if (0==phi) { *p++=0; } else { *p++ = log(fabs(phi)); }
        }

        return 0;

}

int vol3df_zerobounds(struct vol3d *v)
{
        int x,y,z;

        if (VOL_TYPE_FLOAT!=v->datatype) {
                fprintf(stderr,"vol3df_logabs: datatype is not float\n");
                return -1;
        }

        for (z=0;z<v->nz;z++) {
        for (y=0;y<v->ny;y++) {
                *(float *)vol3d_element(v,0,y,z)=0.0;
                *(float *)vol3d_element(v,v->nx-1,y,z)=0.0;
        }}

        for (z=0;z<v->nz;z++) {
        for (x=0;x<v->nx;x++) {
                *(float *)vol3d_element(v,x,0,z)=0.0;
                *(float *)vol3d_element(v,x,v->ny-1,z)=0.0;
        }}

        for (y=0;y<v->ny;y++) {
        for (x=0;x<v->nx;x++) {
                *(float *)vol3d_element(v,x,y,0)=0.0;
                *(float *)vol3d_element(v,x,y,v->nz-1)=0.0;
        }}

        return 0;

}

/* Create a new vol3df object of given dimensions. Each point
 * in the new volume contains the value of the mean curvature at
 * the corresponding point in the original volume.
 */

struct vol3d *vol3df_new_H(struct vol3d *v)
{
        struct vol3d *vnew=NULL;
        int x,y,z;

        if (VOL_TYPE_FLOAT!=v->datatype) {
                fprintf(stderr,"vol3df_new_H: datatype is not float.\n");
                return NULL;
        }


        if (NULL==(vnew=vol3d_new_copytype(v))) {
                fprintf(stderr,
                      "vol3df_new_H: vol3d_new_copytype failed.\n");
                return NULL;
        }

        for (z=0;z<v->nz;z++) {
        for (y=0;y<v->ny;y++) {
        for (x=0;x<v->nx;x++) {
                *(float *)vol3d_element(vnew,x,y,z)
                    = vol3df_mean_curvature(v,(float)x,(float)y,(float)z);
        }}}

        return vnew;


}

/* Create a new vol3df object of given dimensions. The value at each point
 * in the new volume is found by trilinearly interpolating the data
 * in the original volume.
 */

struct vol3d *vol3df_new_scale_interpolated(struct vol3d *v, int newx, int newy, int newz)
{
        struct vol3d *vnew=NULL;
        int x,y,z;
        float rx,ry,rz;


        if (VOL_TYPE_FLOAT!=v->datatype) {
                fprintf(stderr,"vol3df_downsample: datatype is not float.\n");
                return NULL;
        }


        if (NULL==(vnew=vol3df_new(newx,newy,newz))) {
                fprintf(stderr,
                      "vol3df_new_scale_interpolated: vol3df_new failed.\n");
                return NULL;
        }


        for (z=0;z<newz;z++) {
        for (y=0;y<newy;y++) {
        for (x=0;x<newx;x++) {
                rx = ((float)x/(float)newx)*(float)(v->nx);
                ry = ((float)y/(float)newy)*(float)(v->ny);
                rz = ((float)z/(float)newz)*(float)(v->nz);
                *(float *)vol3d_element(vnew,x,y,z)
                        = vol3df_trilinear_point(v,rx,ry,rz);
        }}}

        return vnew;


}

/* Create a new vol3df object of the same dimensions as the original one.
 * Coordinates are mapped between the two objects by the given matrix.
 * Scalar values are trilinearly interpolated.
 */

struct vol3d *vol3df_new_transformed(struct vol3d *v,float *mat)
{
    struct vol3d *vnew=NULL;
    int x,y,z;
    int nx,ny,nz;
    float rx,ry,rz;


    if (VOL_TYPE_FLOAT!=v->datatype) {
            fprintf(stderr,"vol3df_downsample: datatype is not float.\n");
            return NULL;
    }


    if (NULL==(vnew=vol3d_new_copytype(v))) {
            fprintf(stderr,"vol3df_new_transformed: vol3d_copytype failed\n");
            return NULL;
    }
    nx=vnew->nx;
    ny=vnew->ny;
    nz=vnew->nz;

    for (z=0;z<nz;z++) {
    for (y=0;y<ny;y++) {
    for (x=0;x<nx;x++) {
        rx= mat[0]*x + mat[1]*y + mat[2]*z;
        ry= mat[3]*x + mat[4]*y + mat[5]*z;
        rz= mat[6]*x + mat[7]*y + mat[8]*z;
        *(float *)vol3d_element(vnew,x,y,z)
                    = vol3df_trilinear_point(v,rx,ry,rz);
    }}}

    return vnew;

}

struct vol3d *vol3df_new_downsample(struct vol3d *v, int subx, int suby, int subz)
{
        struct vol3d *vnew=NULL;
        unsigned int nx,ny,nz;
        unsigned int xs,ys,zs; /* Coords in subvolume */
        unsigned int x,y,z; /* Coords in original volume */

        if (VOL_TYPE_FLOAT!=v->datatype) {
                fprintf(stderr,"vol3df_downsample: datatype is not float.\n");
                return NULL;
        }

        nx=v->nx; ny=v->ny; nz=v->nz; 

        if ( (0!=(nx%subx)) || (0!=(ny%suby)) || (0!=(nz%subz)) ) {
                fprintf(stderr,
                "vol3df_downsample: volume must be multiple of subvolume.\n");
                return NULL;
        }

        if (NULL==(vnew=vol3df_new_zero(nx/subx,ny/suby,nz/subz))) {
                fprintf(stderr,"vol3df_downsample: vol3df_new_zero failed.\n");
                return NULL;
        }


        for (z=0;z<nz;z++) {
        for (y=0;y<ny;y++) {
        for (x=0;x<nx;x++) {
                xs = x/subx; ys = y/suby; zs = z/subz;
                *(float *)vol3d_element(vnew,xs,ys,zs)
                        +=
                *(float *)vol3d_element(v,x,y,z) ;
        }}}

        return vnew;


}

#ifdef HAVE_PNG
int vol3df_write_pngs(struct vol3d *v, enum vol3d_direction dir, char **fnames )
{
        struct vol3df_stats *stats=NULL;
        unsigned int i,imax=0;

        if (NULL==(stats=vol3df_getstats(v))) {
                fprintf(stderr,
                     "vol3df_write_pngs: vol3df_stats() failed\n");
                return -1;
        }

        switch(dir) {
                case VOL3D_Z: imax=v->nz; break;
                case VOL3D_Y: imax=v->ny; break;
                case VOL3D_X: imax=v->nx; break;
        }

        for (i=0;i<imax;i++) {
               if(0==vol3df_write_png(v,dir,i,fnames[i],stats->max,stats->min)){
                       fprintf(stderr,"Wrote \"%s\".\n",fnames[i]);
               } else {
                       fprintf(stderr,"Unable to build \"%s\".",fnames[i]);
               }
        }

        free(stats);
        return 0;

}

int vol3df_write_png(struct vol3d *v, enum vol3d_direction dir,
                unsigned int index, char *filename, float phimax, float phimin)
{
        struct vol2d *slice=NULL;

        if (NULL==(slice=vol3d_slice_2d(v,dir,index))) {
                fprintf(stderr,"vol3df_write_png: vol3d_slice_2d failed.\n");
                return -1;
        }

        if (0!=vol2df_write_png(slice,filename,phimax,phimin)) {
                fprintf(stderr,
                        "vol3df_write_png: vol2df_write_png failed.\n");
                vol2d_destroy(slice);
                return -1;
        }
        vol2d_destroy(slice);
        return 0;

}

int vol3df_write_png_normalize(struct vol3d *v, enum vol3d_direction dir,
                unsigned int index, char *fname)
{
        struct vol3df_stats *stats=NULL;

        if (NULL==(stats=vol3df_getstats(v))) {
                fprintf(stderr,
                     "vol3df_write_png_normalize: vol3df_stats() failed\n");
                return -1;
        }

        if (0!=vol3df_write_png(v,dir,index,fname,stats->max,stats->min)) {
                fprintf(stderr,
                     "vol3df_write_png_normalize: vol3df_write_png failed.\n");
                free(stats);
                return -1;
        }

        free(stats);
        return 0;

}
#endif /* HAVE_PNG */

/* Take a vol3d volume and a direction: X, Y, or Z. If the direction is X, 
 * (coord0,coord1)=(y,z) ; Y->(z,x) ; Z->(x,y).
 * Allocate an array of floats and fill it with the values along the
 * corresponding line; return this.
 * Return NULL on error.
 */

float *vol3df_ortholine(struct vol3d *v, enum vol3d_direction dir, unsigned int coord0, unsigned int coord1)
{
    unsigned int asize=0;
    float *arr=NULL;

    switch(dir) {
        case VOL3D_X: asize = v->nx; break;
        case VOL3D_Y: asize = v->ny; break;
        case VOL3D_Z: asize = v->nz; break;
    }

    if (NULL==(arr=malloc(sizeof(float)*asize)))
        { perror("vol3df_ortholine: malloc"); return NULL; }

    switch(dir) {
        int i;
        case VOL3D_X:
             for (i=0;i<asize;i++) { arr[i]=vol3df_wrap(v,i,coord0,coord1);}
             break;
        case VOL3D_Y:
             for (i=0;i<asize;i++) { arr[i]=vol3df_wrap(v,coord1,i,coord0);}
             break;
        case VOL3D_Z:
             for (i=0;i<asize;i++) { arr[i]=vol3df_wrap(v,coord0,coord1,i);}
             break;
    }
    return arr;
}

