#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include <rpc/rpc.h>

#include "vol3df.h"
#include "mcubes_tables.h"

#include "growarr.h"

#include "LookUpTable.h"
#include "vector.h"
#include "isosurface.h"


static int mcubes_xdr_uint32(XDR *xdrs,uint32_t u)
{
    unsigned uu;
    uu=u;
    return xdr_u_int(xdrs,&uu);
}



/* marching cubes */

void mcubes_geometry_dump(struct mcubes_state *mc, FILE *fh)
{
        int i;
        int nvertices,ntris;
        fprintf(fh,"Dumping geometry data (open wide..)\n");

        nvertices = growarr_nitems(mc->vertex);
        ntris = growarr_nitems(mc->triangle);

        fprintf(fh,"Vertices\n========\n");
        for (i=0;i<nvertices;i++) {
                struct vector v;
                v = *(struct vector *)growarr_get_nth(mc->vertex,i);
                fprintf(fh,"%d %.12f %.12f %.12f\n",i,v.x,v.y,v.z);
        }
        fprintf(fh,"Triangles\n========\n");
        for (i=0;i<ntris;i++) {
                struct tri t;
                t = *(struct tri *)growarr_get_nth(mc->triangle,i);
                fprintf(fh,"%d %u %u %u\n",i,t.a,t.b,t.c);
        }
        fprintf(fh,"Geometry dump done.\n\n");


}

/* Calculate Euler characteristic chi = V - E + F, and write
 * into "chi" pointer. Return 0 on success, nonzero on failure.
 *
 * For a manifold surface composed of triangles, each edge will correspond
 * to exactly two triangles; each triangle therefore contributes 3/2 edges.
 * Therefore E=(2/3)*F, so chi = V-F/2.
 */
int mcubes_chi(struct mcubes_state *mc,float *chi)
{
        int ntris=-1;

        if (NULL==(mc->triangle)) {
                fprintf(stderr,"mcubes_chi: No triangle data\n");
                return -1;
        } else {
                ntris=growarr_nitems(mc->triangle);
        }
        if (NULL==(mc->vertex)) {
                fprintf(stderr,"mcubes_chi: No vertex data\n");
                return -1;
        }
        *chi = (float)mc->nvertices-0.5*(float)ntris;
        return 0;

}

int mcubes_genus(struct mcubes_state *mc, float *g)
{
    float chi;
    if (0!=mcubes_chi(mc,&chi)) {
        fprintf(stderr,"mcubes_genus: mcubes_chi failed.\n"); return -1;
    }
    *g = 0.5*(2.0-chi);
    if (fabs( (*g)-floor(*g)) > FLT_EPSILON ) {
        fprintf(stderr,"*** WARNING genus=%f\n",*g);
    }
    return 0;
}

void mcubes_info(struct mcubes_state *mc,FILE *fh)
{
        int nverts=-1,ntris=-1;

        fprintf(fh,"isoval=%.12f\n",mc->isoval);
        fprintf(fh,"dims %dx%dx%d\n",mc->nx,mc->ny,mc->nz);
        if (NULL==mc->voxel) {
                fprintf(fh,"No voxel data\n");
        } else {
                fprintf(fh,"voxel data %dx%dx%d\n",
                                mc->voxel->nx,
                                mc->voxel->ny,
                                mc->voxel->nz);
        }

        if (NULL==mc->colvol) {
                fprintf(fh,"No colour volume\n");
        } else {
                fprintf(fh,"colour data %dx%dx%d\n",
                                mc->voxel->nx,
                                mc->voxel->ny,
                                mc->voxel->nz);
        }

        if (NULL==(mc->triangle)) {
                fprintf(fh,"No triangle data\n");
        } else {
                fprintf(fh,"%u triangles\n",ntris=growarr_nitems(mc->triangle));
        }
        if (NULL==(mc->vertex)) {
                fprintf(fh,"No vertex data\n");
        } else {
                fprintf(fh,"%u vertices\n",nverts=growarr_nitems(mc->vertex));
        }
#ifdef DEBUG
        if (NULL==(mc->alias_map)) {
                fprintf(fh,"No alias map data\n");
        } else {
                fprintf(fh,"%u alias maps\n",nverts=growarr_nitems(mc->alias_map));
        }
#endif /* DEBUG */
        if (NULL==(mc->normal)) {
                fprintf(fh,"No normal data\n");
        } else {
                fprintf(fh,"%u normals\n",growarr_nitems(mc->normal));
        }
        if (NULL==(mc->scalar)) {
                fprintf(fh,"No scalar data\n");
        } else {
                fprintf(fh,"%u scalars\n",growarr_nitems(mc->scalar));
        }
        if (NULL==(mc->scalarstats)) {
            fprintf(fh,"No scalar statistics\n");
        } else {
            fprintf(stderr,"Scalar max %f min %f mean %f var %f\n",
                    mc->scalarstats->max,
                    mc->scalarstats->min,
                    mc->scalarstats->mean,
                    mc->scalarstats->variance);
        }
        if (NULL==(mc->xvert)) {
                fprintf(fh,"No xvertex data\n");
        } else {
                fprintf(fh,"xvertex data %dx%dx%d\n",
                                mc->xvert->nx,
                                mc->xvert->ny,
                                mc->xvert->nz);
        }
        if (NULL==(mc->yvert)) {
                fprintf(fh,"No yvertex data\n");
        } else {
                fprintf(fh,"yvertex data %dx%dx%d\n",
                                mc->yvert->nx,
                                mc->yvert->ny,
                                mc->yvert->nz);
        }
        if (NULL==(mc->zvert)) {
                fprintf(fh,"No zvertex data\n");
        } else {
                fprintf(fh,"zvertex data %dx%dx%d\n",
                                mc->zvert->nx,
                                mc->zvert->ny,
                                mc->zvert->nz);
        }
        if ((nverts!=-1)&&(ntris!=-1)) {
                float chi;
                fprintf(stderr,"%u unique vertices\n",mc->nvertices);
                chi = (float)mc->nvertices-0.5*(float)ntris;
                fprintf(fh,"Euler number is %.12f, assuming manifold surface.\n",
                                chi);
                fprintf(fh,"Corresponding genus is %.12f\n",0.5*(2-chi));
        }
        if (mc->calc_area) {
                fprintf(stderr,"Accumulated area is %.12f\n",mc->total_area);
        }
        if (mc->calc_mean_curvature) {
            float meanH,meanH2;
            meanH = mc->sum_H_area / mc->total_area;
            meanH2= mc->sum_H2_area / mc->total_area;
                fprintf(stderr,"Averaged mean curvature is %.12f\n",meanH);
                fprintf(stderr,"Averaged squared mean curvature is %.12f\n",
                        meanH2);
        }
        if (mc->calc_gaussian_curvature) {
            float meanK,meanK2;
            meanK = mc->sum_K_area / mc->total_area;
            meanK2= mc->sum_K2_area / mc->total_area;
                fprintf(stderr,"Integrated Gaussian curvature is %.12f\n",
                        mc->sum_K_area);
                fprintf(stderr,"Averaged Gaussian curvature is %.12f\n",meanK);
                fprintf(stderr,"Averaged squared Gaussian curvature is %.12f\n",
                        meanK2);
        }

}

/* Call before triangulation to enable counting up the area */
void mcubes_count_area(struct mcubes_state *mc)
{
    mc->calc_area=1;
}

/* Enable scalar summing */
void mcubes_count_scalar(struct mcubes_state *mc)
{
    mcubes_count_area(mc);
    mc->calc_sum_scalar=1;
}

/* Enable curvature counting */

void mcubes_count_curvature(struct mcubes_state *mc)
{
    mcubes_count_area(mc);
    mc->calc_mean_curvature=1;
    mc->calc_gaussian_curvature=1;
}

int mcubes_scalar_curvature(struct mcubes_state *mc)
{
    fprintf(stderr,"Setting SC calculation ON\n");fflush(stderr);
    if (NULL==(mc->scalar=growarr_new(sizeof(float),1))) {
        fprintf(stderr,"mcubes_scalar_curvature: growarr_new failed\n");
        return -1;
    }

    mc->scalar_curvature = 1;
    return 0;
}

int mcubes_gaussian_curvature(struct mcubes_state *mc)
{
    mc->scalar_gaussian = 1;
    return 0;
}

void mcubes_destroy(struct mcubes_state *mc)
{
        if (NULL==mc) { return; }

        mcubes_release_input_data(mc);
        mcubes_purge_cache(mc);
        mcubes_purge_geometry(mc);
        free(mc);
        return;
}

struct mcubes_state *mcubes_new_fromsurf(char *filename)
{
    struct mcubes_state *mc=NULL;
    FILE *fh=NULL;

    if (NULL==(fh=fopen(filename,"rb"))) {
        perror("mcubes_new_fromsurf: fopen");
        return NULL;
    }

    if (NULL==(mc=mcubes_new_fromsurf_fh(fh))) {
        fprintf(stderr,"mcubes_new_fromsurf: mcubes_new_fromsurf_fh failed.\n");
        fclose(fh);
        return NULL;
    }
    fclose(fh);
    return mc;

}

struct mcubes_state *mcubes_new_fromsurf_fh(FILE *fh)
{
    struct mcubes_state *mc=NULL;
    uint32_t i;
    uint32_t nvertices,ntriangles,nnormals;
    unsigned u;
    XDR xdrs;

    if (NULL==(mc=mcubes_new_minimal())) {
        fprintf(stderr,"mcubes_new_fromsurf_fh: mcubes_new_minimal failed\n");
        return NULL;
    }

    xdrstdio_create(&xdrs,fh,XDR_DECODE);


    xdr_u_int(&xdrs,&u);
    if (SURF_BLOCKTYPE_FILE != u) {
        fprintf(stderr,"mcubes_new_fromsurf_fh: bad file blocktype %x\n",
                (unsigned) u);
        goto error;
    }

    /* Read block sizes */

    xdr_u_int(&xdrs,&u); nvertices=u;
    xdr_u_int(&xdrs,&u); ntriangles=u;
    xdr_u_int(&xdrs,&u); nnormals=u;


    if (0!=nnormals) {
        if (nvertices != nnormals) {
            fprintf(stderr,"mcubes_new_fromsurf_fh: "
                    "(nvertices=%u) != (nnormals=%u)\n",
                    nvertices,nnormals);
            goto error;
        }
    }

    /* Allocate storage */

    if (NULL==(mc->triangle=growarr_new(sizeof(struct tri),ntriangles)))
            { fprintf(stderr,"growarr_new failed\n"); goto error; }
    if (NULL==(mc->vertex=growarr_new(sizeof(vector),ntriangles)))
            { fprintf(stderr,"growarr_new failed\n"); goto error; }
    if (0!=nnormals) {
        if (NULL==(mc->normal=growarr_new(sizeof(vector),nnormals)))
                { fprintf(stderr,"growarr_new failed\n"); goto error; }
    }

    fprintf(stderr,"%u triangles %u vertices %u normals\n",
            ntriangles,nvertices,nnormals);fflush(stderr);

    /* Read number of unique vertices */
    xdr_u_int(&xdrs,&u); mc->nvertices = u;

    /* Read vertex block */

        xdr_u_int(&xdrs,&u);
        if (SURF_BLOCKTYPE_VERTICES != u) {
            fprintf(stderr,"mcubes_new_fromsurf_fh: bad vertex blocktype %x\n",
                    u);
            goto error;
        }

        /* Read vertices */

        for (i=0;i<nvertices;i++) {
            struct vector v;
            xdr_float(&xdrs,&(v.x));
            xdr_float(&xdrs,&(v.y));
            xdr_float(&xdrs,&(v.z));
            growarr_add(mc->vertex,&v);
        }

    /* Read triangle block */

        xdr_u_int(&xdrs,&u);
        if (SURF_BLOCKTYPE_TRIANGLES != u) {
            fprintf(stderr,"mcubes_new_fromsurf_fh: "
                    "bad triangle blocktype %x\n", u);
            goto error;
        }

        /* Read triangles */

        for (i=0;i<ntriangles;i++) {
                struct tri t;
                xdr_u_int(&xdrs,&u); t.a = (uint32_t) u;
                xdr_u_int(&xdrs,&u); t.b = (uint32_t) u;
                xdr_u_int(&xdrs,&u); t.c = (uint32_t) u;

                growarr_add(mc->triangle,&t);
        }

    if (NULL!=mc->normal) {
        xdr_u_int(&xdrs,&u);
        if (SURF_BLOCKTYPE_NORMALS != u) {
            fprintf(stderr,"mcubes_new_fromsurf_fh: bad vertex blocktype %x\n",
                    u);
            goto error;
        }

        /* Read vertices */

        for (i=0;i<nvertices;i++) {
            struct vector v;
            xdr_float(&xdrs,&(v.x));
            xdr_float(&xdrs,&(v.y));
            xdr_float(&xdrs,&(v.z));
            growarr_add(mc->normal,&v);
        }
    }

    xdr_destroy(&xdrs);

    return mc;


error:
    mcubes_destroy(mc);
    return NULL;
}

struct mcubes_state *mcubes_new_minimal(void)
{
    struct mcubes_state *mc=NULL;
    if (NULL==(mc=malloc(sizeof(struct mcubes_state))))
            { perror("mcubes_new_minimal: malloc"); return NULL; }

    /* All pointers should be NULL by default */
    mc->triangle = mc->vertex = mc->normal = mc->scalar = NULL;
    mc->xvert = mc->yvert = mc->zvert = NULL;

    mc->nx = 0.0; mc->ny =0.0; mc->nz = 0.0;

    mc->voxel = NULL;
    mc->isoval = 0.0;
    mc->nvertices = 0;

    mc->calc_area=0;
    mc->calc_mean_curvature=0;
    mc->calc_gaussian_curvature=0;
    mc->calc_sum_scalar=0;
    mc->total_area=0;
    mc->fatal_error=0;
    mc->sum_H_area=0;
    mc->sum_H2_area=0;
    mc->sum_H4_area=0;
    mc->sum_scalar_area=0;
    mc->sum_scalar2_area=0;
    mc->scalar_curvature=0;
    mc->scalar_gaussian=0;

    return mc;
}

struct mcubes_state *mcubes_new(struct vol3d *voxel, float isoval,
                struct vol3d *colvol, int calc_normals)
{
        struct mcubes_state *mc=NULL;

        if (NULL==(mc=mcubes_new_minimal())) {
            fprintf(stderr,"mcubes_new: mcubes_new_minimal failed.\n");
            return NULL; 
        }

        /* If we are going to calculate normals, then allocate space for them */
        if (calc_normals) {
                if (NULL==(mc->normal=growarr_new(sizeof(vector),1)))
                        {fprintf(stderr,"growarr_new failed\n"); goto bailout;}
        }

        mc->nx = voxel->nx; mc->ny =voxel->ny; mc->nz = voxel->nz;

        mc->voxel = voxel;
        mc->isoval = isoval;

        /* If a "colvol" argument is passed, then assume that it is a vol3df
         * object, whose values determine the colour at each vertex.
         * This means that we need to store an extra scalar, the interpolated
         * value of the "colvol" field, at each vertex.
         * This can later be turned into an RGBA tuple for colouring the
         * geometry.
         */

        mc->colvol = colvol;
        mc->scalarstats = NULL;
        if (NULL != colvol) {
                if ( 
                          (mc->nx != colvol->nx)
                        ||(mc->ny != colvol->ny)
                        ||(mc->nz != colvol->nz)) {
                        fprintf(stderr, "mcubes_new: "
                            "colour and voxels fields are different sizes!\n");
                        goto bailout;
                }

                if (NULL==(mc->scalar=growarr_new(sizeof(float),1)))
                    {fprintf(stderr,"growarr_new failed\n"); goto bailout;}
                if (NULL==(mc->scalarstats=vol3df_getstats(colvol)))
                    {fprintf(stderr,"vol3df_getstats failed\n");goto bailout;}

                    
        }

        /* Allocate storage for triangles and vertices */

        if (NULL==(mc->triangle=growarr_new(sizeof(struct tri),1)))
                { fprintf(stderr,"growarr_new failed\n"); goto bailout; }
        if (NULL==(mc->vertex=growarr_new(sizeof(vector),1)))
                { fprintf(stderr,"growarr_new failed\n"); goto bailout; }

#ifdef DEBUG
        if (NULL==(mc->alias_map=growarr_new(sizeof(vertindex),1)))
                { fprintf(stderr,"growarr_new failed\n"); goto bailout; }
#endif /* DEBUG */

        /* Allocate temporary storage for vertex indices. */

        /* NB vertex indices are currently unsigned 32-bit integers.
         * That's enough to handle a 1024^3 volume with three vertices
         * per voxel, but larger volumes would mean we need to move
         * to 64 bits, and almost certainly a less wasteful storage method
         */

        if (NULL==(mc->xvert=vol3d_new(1+mc->nx,1+mc->ny,1+mc->nz,
                                        VOL_TYPE_UINT32)))
         { fprintf(stderr,"mcubes_new: vol3d_new failed\n");goto bailout;}
        if (NULL==(mc->yvert=vol3d_new(1+mc->nx,1+mc->ny,1+mc->nz,
                                        VOL_TYPE_UINT32)))
         { fprintf(stderr,"mcubes_new: vol3d_new failed\n");goto bailout;}
        if (NULL==(mc->zvert=vol3d_new(1+mc->nx,1+mc->ny,1+mc->nz,
                                        VOL_TYPE_UINT32)))
         { fprintf(stderr,"mcubes_new: vol3d_new failed\n");goto bailout;}

        vol3d_set(mc->xvert,0); vol3d_set(mc->yvert,0); vol3d_set(mc->zvert,0);


#ifdef DEBUG
        /* Fill vertex volumes with an "invalid" value. This assumes that
         * we'll have no more than V_BADVAL vertices, but makes debugging
         * the cases where nonexistent vertices are referenced a lot easier.
         */
    
        {
             int x,y,z;
             for (z=0;z<=mc->nz;z++) {
             for (y=0;y<=mc->ny;y++) {
             for (x=0;x<=mc->nx;x++) {
                     *(vertindex *) vol3d_element(mc->xvert,x,y,z) = V_BADVAL;
                     *(vertindex *) vol3d_element(mc->yvert,x,y,z) = V_BADVAL;
                     *(vertindex *) vol3d_element(mc->zvert,x,y,z) = V_BADVAL;
             }}}
        }
#endif /* DEBUG */

        return mc;

bailout:
        mcubes_destroy(mc);
        return NULL;

}

void mcubes_cache_vertex(struct mcubes_state *mc,struct vol3d *vertvol,
                int ix, int iy, int iz,
                vertindex v)
{
        *(vertindex *)vol3d_element(vertvol,ix,iy,iz) = v;

#ifdef DEBUG
        /* Maintain alias map */
        if ((ix==mc->nx)||(iy==mc->ny)||(iz==mc->nz)) {
                /* in halo region */
            vertindex vreal;
                vreal = *(vertindex *)vol3d_element(vertvol,
                                ix%mc->nx,
                                iy%mc->ny,
                                iz%mc->nz);
                growarr_add(mc->alias_map,&vreal);
#ifdef SICKENINGLY_VERBOSE
                fprintf(stderr,"Aliasing %u (%d %d %d) to %u (%d %d %d)\n",
                                v,ix,iy,iz,
                                vreal,
                                ix%mc->nx,
                                iy%mc->ny,
                                iz%mc->nz);
#endif /* SICKENINGLY_VERBOSE */
        } else {
                growarr_add(mc->alias_map,&v);
        }
#endif /* DEBUG */

#ifdef SICKENINGLY_VERBOSE
        {
        char vi;
                if (vertvol == mc->xvert) {
                        vi='x';
                } else if (vertvol == mc->yvert) {
                        vi='y';
                } else if (vertvol == mc->zvert) {
                        vi='z';
                } else {
                        vi='?';
                }

        fprintf(stderr,"%d %d %d = %c vert %u\n", ix,iy,iz,vi,v);
        }
#endif /* SICKENINGLY_VERBOSE */
}


vertindex mcubes_add_vertex(struct mcubes_state *mc, 
                float x, float y, float z)
{
        vector v;
        vertindex index;
        int i;

        v.x=x; v.y=y;v.z=z;

        if (0>(i=growarr_add(mc->vertex,&v))) {
                fprintf(stderr,"mcubes_add_vertex: growarr_add failed.\n");
                return -1;
        }
        index=i;


        if (NULL!=mc->normal) {
                vector norm;
                vol3df_trilinear_gradient(mc->voxel,x,y,z,
                                &norm.x,&norm.y,&norm.z);
                norm = vector_normalized(norm);
                if (index!=growarr_add(mc->normal,&norm)) {
                        fprintf(stderr,
                            "mcubes_add_vertex: growarr_add normal failed.\n");

                        return -1;
                }
        } /* if (normals) */

        if (NULL!=mc->colvol) {
                float scalar;
                scalar = vol3df_trilinear_point(mc->colvol,x,y,z);
                if (index != growarr_add(mc->scalar,&scalar)) {
                        fprintf(stderr,
                             "mcubes_add_vertex: growarr_add scalar failed.\n");
                        return -1;
                }
        } /* if (scalar colouring) */

        if (mc->scalar_curvature) {
            /* Calculate and store curvature */

            float H;
            if (NULL==mc->scalar) {
                fprintf(stderr,"scalar_curvature but no mc->scalar!\n");
                fprintf(stderr,"Shouldn't happen; bailing out.\n");
                exit(-1);
            }

            if (mc->scalar_gaussian) {
                H = vol3df_gaussian_curvature(mc->voxel,v.x,v.y,v.z);
            } else {
                H = vol3df_mean_curvature(mc->voxel,v.x,v.y,v.z); 
            }

            if (0 > growarr_add(mc->scalar,&H)) {
                    fprintf(stderr,
                     "mcubes_add_vertex: growarr_add curvature failed.\n");
                    return -1;
            }

        }

#ifdef SICKENINGLY_VERBOSE
        {
            vertindex vi;
            vi = growarr_nitems(mc->vertex)-1;
        fprintf(stderr,"Vertex %d : %f %f %f\n", vi,x,y,z);
        fflush(stderr);

                }
#endif /* SICKENINGLY_VERBOSE */

        return index;
}

/* Calculate and store all the points where the isosurface crosses the
 * cuberille.
 *
 * Return 0 on success. If nonzero is returned, then it should be assumed
 * that the calculation has failed (eg due to lack of memory), and
 * the marching cubes calculation should be abandoned.
 */
int mcubes_calc_vertices(struct mcubes_state *mc)
{
        int x,y,z;
        float isoval;
        vertindex vert;   /* Game, cat? */

        isoval = mc->isoval;

        /* NB: each ordinate ranges from 0 to N, *not* 0 to N-1.
         * The vertex arrays are of size (nx+1)*(ny+1)*(nz+1) to account
         * for periodic boundaries.
         */
        for (z=0;z<1+mc->nz;z++) {
        for (y=0;y<1+mc->ny;y++) {
        for (x=0;x<1+mc->nx;x++) {
               float delta;
               float phi,phix,phiy,phiz;
               int p,px,py,pz;
               int xhalo=0,yhalo=0,zhalo=0;
               int is_in_halo=0; /* Set to 1 if (x,y,z) in halo region */

               if ((x==mc->nx)||(y==mc->ny)||(z==mc->nz)) {
                       is_in_halo=1;
                       if (x==mc->nx) { xhalo=1; }
                       if (y==mc->ny) { yhalo=1; }
                       if (z==mc->nz) { zhalo=1; }
               }

               phi  = vol3df_wrap(mc->voxel,x,y,z)   - isoval; 
               phix = vol3df_wrap(mc->voxel,x+1,y,z) - isoval; 
               phiy = vol3df_wrap(mc->voxel,x,y+1,z) - isoval; 
               phiz = vol3df_wrap(mc->voxel,x,y,z+1) - isoval; 

               /* If anything is hovering around the isovalue, coerce
                * it to a point above the isosurface.
                */

               if (fabs(phi ) < FLT_EPSILON) { phi  = FLT_EPSILON; }
               if (fabs(phix) < FLT_EPSILON) { phix = FLT_EPSILON; }
               if (fabs(phiy) < FLT_EPSILON) { phiy = FLT_EPSILON; }
               if (fabs(phiz) < FLT_EPSILON) { phiz = FLT_EPSILON; }

               p=(phi<0); px=(phix<0); py=(phiy<0); pz=(phiz<0);

               if ((!xhalo) && ( p != px)) { 
                       /* There's a crossing in the X direction */
                       delta = phi/(phi-phix);
                       vert=mcubes_add_vertex(mc,
                                        (float)x+delta,y,z);
                       mcubes_cache_vertex(mc,mc->xvert,x,y,z,vert);
                       if (!is_in_halo) { mc->nvertices++; }
               }
               if ((!yhalo)&&(p != py)) { 
                       /* There's a crossing in the Y direction */
                       delta = phi/(phi-phiy);
                       vert=mcubes_add_vertex(mc,
                                        x,(float)y+delta,z);
                       mcubes_cache_vertex(mc,mc->yvert,x,y,z,vert);
                       if (!is_in_halo) { mc->nvertices++; }
               }
               if ((!zhalo)&&(p != pz)) { 
                       /* There's a crossing in the Z direction */
                       delta = phi/(phi-phiz);
                       vert=mcubes_add_vertex(mc,
                                        x,y,(float)z+delta);
                       mcubes_cache_vertex(mc,mc->zvert,x,y,z,vert);
                       if (!is_in_halo) { mc->nvertices++; }
               }

        }}} /* xyz */

        if (mc->fatal_error) {
                fprintf(stderr,"mcubes_calc_vertices: fatal error detected.\n");
                return -1;
        }

        return 0;

}

void dump_verts(struct vol3d *v)
{
        int x,y,z;
        if (VOL_TYPE_UINT32 != v->datatype) {
                fprintf(stderr,"dump_verts: wrong datatype\n");
                return;
        }
        for (z=0;z<v->nz;z++) {
        for (y=0;y<v->ny;y++) {
        for (x=0;x<v->nx;x++) {
                uint32_t phi;
                phi = *(uint32_t *) vol3d_element(v,x,y,z);
                if (phi!=V_BADVAL) {
                        fprintf(stderr,"(%d,%d,%d) : %u\n",x,y,z,phi);
                }
        }}}
}

float mcubes_area(struct mcubes_state *mc)
{
    if (!mc->calc_area) { return -1; }
    return mc->total_area;
}

float mcubes_average_curvature(struct mcubes_state *mc)
{
    if (!mc->calc_mean_curvature) { return -1; }

    return mc->sum_H_area / mc->total_area;
}

float mcubes_average_gaussian_curvature(struct mcubes_state *mc)
{
    if (!mc->calc_gaussian_curvature) { return -1; }

    return mc->sum_K_area / mc->total_area;
}

/* Return the sum over all triangles of
 * (average value of scalar over triangle)*(area of triangle)
 */
float mcubes_sum_scalar_area(struct mcubes_state *mc)
{
    if (!mc->calc_sum_scalar) { return -1; }

    return mc->sum_scalar_area;
}

/* Return the sum over all triangles of
 * (average value of scalar^2 over triangle)*(area of triangle)
 */
float mcubes_sum_scalar2_area(struct mcubes_state *mc)
{
    if (!mc->calc_sum_scalar) { return -1; }

    return mc->sum_scalar2_area;
}


float mcubes_meansquared_curvature(struct mcubes_state *mc)
{
    if (!mc->calc_mean_curvature) { return -1; }

    return mc->sum_H2_area / mc->total_area;
}

float mcubes_meansquared_gaussian_curvature(struct mcubes_state *mc)
{
    if (!mc->calc_gaussian_curvature) { return -1; }

    return mc->sum_K2_area / mc->total_area;
}

float mcubes_meanfourth_curvature(struct mcubes_state *mc)
{
    if (!mc->calc_mean_curvature) { return -1; }

    return mc->sum_H4_area / mc->total_area;
}

int mcubes_add_triangle(struct mcubes_state *mc, struct tri *t)
{
        int trino;
        float area;

        if (0>(trino=growarr_add(mc->triangle,t))) {
                fprintf(stderr,"mcubes_add_triangle: "
                            "growarr_add triangle failed.\n");
                return -1;
        }

#ifdef SICKENINGLY_VERBOSE
        fprintf(stderr,"Added triangle %d: %u %u %u\n",
                        trino,t->a,t->b,t->c);
#endif /* SICKENINGLY_VERBOSE */
        if (mc->calc_area) {

            area = triangle_area(
                    *(struct vector *)growarr_get_nth(mc->vertex,t->a),
                    *(struct vector *)growarr_get_nth(mc->vertex,t->b),
                    *(struct vector *)growarr_get_nth(mc->vertex,t->c)
            );
            mc->total_area += area;

            if (mc->calc_mean_curvature) {
                float H=0.0,H2=0.0,H4=0.0;
                float Ha,Hb,Hc,Ka,Kb,Kc;
                vector va,vb,vc;
                va= *(struct vector *)growarr_get_nth(mc->vertex,t->a);
                vb= *(struct vector *)growarr_get_nth(mc->vertex,t->b);
                vc= *(struct vector *)growarr_get_nth(mc->vertex,t->c);

                Ha = vol3df_mean_curvature(mc->voxel,va.x,va.y,va.z);
                Hb = vol3df_mean_curvature(mc->voxel,vb.x,vb.y,vb.z);
                Hc = vol3df_mean_curvature(mc->voxel,vc.x,vc.y,vc.z);

                H=(Ha+Hb+Hc)/3.0;
                H2=(Ha*Ha + Hb*Hb + Hc*Hc)/3.0;
                H4=( Ha*Ha*Ha*Ha + Hb*Hb*Hb*Hb + Hc*Hc*Hc*Hc )/3.0;

                /* It's possible to just add up H at each vertex, and then
                 * divide by the number of vertices to obtain a measure of
                 * the average mean curvature. However, there will be more
                 * triangles (and therefore more vertices) in regions of high
                 * curvature, which will skew the statistics. Instead, weight
                 * each curvature measurement by the area of the triangle
                 * over which it was averaged.
                 */

                mc->sum_H_area  += H*area;
                mc->sum_H2_area += H2*area;
                mc->sum_H4_area += H4*area;

                if (mc->calc_gaussian_curvature) {
                    float K=0.0,K2=0.0;
                    Ka = vol3df_gaussian_curvature(mc->voxel,va.x,va.y,va.z);
                    Kb = vol3df_gaussian_curvature(mc->voxel,vb.x,vb.y,vb.z);
                    Kc = vol3df_gaussian_curvature(mc->voxel,vc.x,vc.y,vc.z);

                    K=(Ka+Kb+Kc)/3.0;
                    K2=(Ka*Ka + Kb*Kb + Kc*Kc)/3.0;

                    mc->sum_K_area  += K*area;
                    mc->sum_K2_area += K2*area;
                }

            }

            if (mc->calc_sum_scalar) {
                float sum_scalar=0.0;
                float sum_scalar2=0.0;
                float phi;

                phi = *(float *)growarr_get_nth(mc->scalar,t->a);
                sum_scalar += phi ; sum_scalar2 += phi*phi;
                phi = *(float *)growarr_get_nth(mc->scalar,t->b);
                sum_scalar += phi ; sum_scalar2 += phi*phi;
                phi = *(float *)growarr_get_nth(mc->scalar,t->c);
                sum_scalar += phi ; sum_scalar2 += phi*phi;

                sum_scalar /= 3.0;
                sum_scalar2 /= 3.0;

                mc->sum_scalar_area  += sum_scalar*area;
                mc->sum_scalar2_area += sum_scalar2*area;
            }
        }
        return 0;
}

vertindex mcubes_get_edge_vertex(struct mcubes_state *mc,
        int x, int y, int z, int edgeno)
{
        const int *edge_delta=NULL;
        struct vol3d *vertvol=NULL;

        edge_delta = cube_edge_delta[edgeno];

        /* Determine whether the vertex lies on a link running in the
         * X, Y, or Z direction.
         */
        switch(edge_delta[3]) {
                case 0: vertvol = mc->xvert; break;
                case 1: vertvol = mc->yvert; break;
                case 2: vertvol = mc->zvert; break;
                default: fprintf(stderr,
                               "mcubes_get_edge_vertex: bad edge table lookup "
                               "(%d,%d,%d) edge %d !\n",x,y,z,edgeno);
                         return 0;
        }
        return *(vertindex *)vol3d_element_wrap(vertvol,
                x + edge_delta[0],
                y + edge_delta[1],
                z + edge_delta[2]);

}


int mcubes_write_vtkfh(struct mcubes_state *mc, FILE *f, int binary)
{
        int i,nverts,ntris;
        XDR xdrs;

        nverts = growarr_nitems(mc->vertex);
        ntris = growarr_nitems(mc->triangle);

        if (binary) { xdrstdio_create(&xdrs,f,XDR_ENCODE); }

        /* Write header */

        fprintf(f,"# vtk DataFile Version 3.0\n");
        fprintf(f,"Generated from vol3d isosurface\n");
        fprintf(f, binary ? "BINARY\n" : "ASCII\n");
        fprintf(f,"DATASET POLYDATA\n");
        fprintf(f,"POINTS %d float\n",nverts);

        /* Write vertices */

        for (i=0;i<nverts;i++) {
                struct vector v;

                v = *(struct vector *)growarr_get_nth(mc->vertex,i);

                if (binary) {
                        xdr_float(&xdrs,&v.x);
                        xdr_float(&xdrs,&v.y);
                        xdr_float(&xdrs,&v.z);
                } else {
                        fprintf(f,"%.12f %.12f %.12f\n",v.x,v.y,v.z);
                }
        }

        /* Write triangles */
        fprintf(f,"\nPOLYGONS %d %d\n",ntris,ntris*4);
        for (i=0;i<ntris;i++) {
                struct tri t;
                t = *(struct tri *)growarr_get_nth(mc->triangle,i);

                if (binary) {
                        int n=3;
                        int ta,tb,tc;
                        ta=t.a;tb=t.b;tc=t.c;
                        xdr_int(&xdrs,&n);
                        xdr_int(&xdrs,&ta);
                        xdr_int(&xdrs,&tb);
                        xdr_int(&xdrs,&tc);
                } else {
                        fprintf(f,"3 %u %u %u\n", t.a,t.b,t.c);
                }
        }
        if (NULL!=mc->normal) {
                /* Write normals */
                int nnorms = growarr_nitems(mc->normal);
                fprintf(f,"\nPOINT_DATA %d\n",nnorms);
                /* NB assume normals in 1:1 corresp with vertices */
                fprintf(f,"NORMALS cell_normals float\n");
                for (i=0;i<nnorms;i++) {
                        struct vector v;

                        v = *(struct vector *)growarr_get_nth(mc->normal,i);

                        if (binary) {
                                xdr_float(&xdrs,&v.x);
                                xdr_float(&xdrs,&v.y);
                                xdr_float(&xdrs,&v.z);
                        } else {
                                fprintf(f,"%.12f %.12f %.12f\n",v.x,v.y,v.z);
                        }
                }
        }
        
        return 0;
}

/* ------------  Modified from Thomas Lewiner's C++ code --------------------*/

static void print_cube(float *v) 
{
        fprintf(stderr,"\t%.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f\n",
        v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7]) ;
}

/*
 * Test a face
 * if face>0 return true if the face contains a part of the surface
 */
static int test_face( int face, float *v )
{
  float A,B,C,D ;

  switch( face )
  {
  case -1 : case 1 :  A = v[0] ;  B = v[4] ;  C = v[5] ;  D = v[1] ;  break ;
  case -2 : case 2 :  A = v[1] ;  B = v[5] ;  C = v[6] ;  D = v[2] ;  break ;
  case -3 : case 3 :  A = v[2] ;  B = v[6] ;  C = v[7] ;  D = v[3] ;  break ;
  case -4 : case 4 :  A = v[3] ;  B = v[7] ;  C = v[4] ;  D = v[0] ;  break ;
  case -5 : case 5 :  A = v[0] ;  B = v[3] ;  C = v[2] ;  D = v[1] ;  break ;
  case -6 : case 6 :  A = v[4] ;  B = v[7] ;  C = v[6] ;  D = v[5] ;  break ;
  default : fprintf(stderr,"test_face: Invalid face code %d\n", face ) ;
            print_cube(v) ;
            A = B = C = D = 0 ;
  };

  return face * A * ( A*C - B*D ) >= 0  ;  /* face and A invert signs */
}

vertindex mcubes_get_x_vert(struct mcubes_state *mc, int x, int y, int z)
{ return *(vertindex *)vol3d_element(mc->xvert,x,y,z); }
vertindex mcubes_get_y_vert(struct mcubes_state *mc, int x, int y, int z)
{ return *(vertindex *)vol3d_element(mc->yvert,x,y,z); }
vertindex mcubes_get_z_vert(struct mcubes_state *mc, int x, int y, int z)
{ return *(vertindex *)vol3d_element(mc->zvert,x,y,z); }

/*
 * Add the triangles for a voxel.
 * The optional v12 argument is the index of the central vertex, if
 * required.
 *
 * Returns 0 on success. Nonzero indicates that an error occurred, and that
 * the triangulation should be abandoned.
 */
int mcubes_triangulate_voxel(struct mcubes_state *mc, const signed char* trig,
                char n, vertindex v12, int i, int j, int k )
{
  vertindex    tv[3] ;
  int t;

  for( t = 0 ; t < 3*n ; t++ )
  {

    switch( trig[t] )
    {
    case  0 : tv[ t % 3 ] = mcubes_get_x_vert(mc, i , j , k ) ; break ;
    case  1 : tv[ t % 3 ] = mcubes_get_y_vert(mc,i+1, j , k ) ; break ;
    case  2 : tv[ t % 3 ] = mcubes_get_x_vert(mc, i ,j+1, k ) ; break ;
    case  3 : tv[ t % 3 ] = mcubes_get_y_vert(mc, i , j , k ) ; break ;
    case  4 : tv[ t % 3 ] = mcubes_get_x_vert(mc, i , j ,k+1) ; break ;
    case  5 : tv[ t % 3 ] = mcubes_get_y_vert(mc,i+1, j ,k+1) ; break ;
    case  6 : tv[ t % 3 ] = mcubes_get_x_vert(mc, i ,j+1,k+1) ; break ;
    case  7 : tv[ t % 3 ] = mcubes_get_y_vert(mc, i , j ,k+1) ; break ;
    case  8 : tv[ t % 3 ] = mcubes_get_z_vert(mc, i , j , k ) ; break ;
    case  9 : tv[ t % 3 ] = mcubes_get_z_vert(mc,i+1, j , k ) ; break ;
    case 10 : tv[ t % 3 ] = mcubes_get_z_vert(mc,i+1,j+1, k ) ; break ;
    case 11 : tv[ t % 3 ] = mcubes_get_z_vert(mc, i ,j+1, k ) ; break ;
    case 12 : tv[ t % 3 ] = v12 ; break ;
    default : break ;
    }

    if( -1 == tv[t%3] )
    {
      fprintf(stderr,"mcubes_triangulate_voxel: invalid triangle\n") ;
      return -1;
    }
    
    if( 2== t%3 )
    {
        struct tri mytri;
        mytri.a = tv[0];
        mytri.b = tv[1];
        mytri.c = tv[2];
        if (0!=mcubes_add_triangle(mc,&mytri)) {
                fprintf(stderr,"mcubes_triangulate_voxel: "
                                "mcubes_add_triangle() failed.\n");
                return -1;
        }
    }
  } /* Loop over triangles for this case */
  return 0;
}

/*
 * In certain cases, it turns out to be best to add an extra vertex, inside
 * the voxel; see Chernyaev's paper.
 * This vertex is defined to be the average of the other vertices on the voxel.
 */

vertindex mcubes_add_c_vertex(struct mcubes_state *mc,
                int x, int y, int z, int cubestate )
{
        struct vector newvert;
        int i=0;
        const signed char *crossings=NULL;

        newvert.x = newvert.y = newvert.z = 0.0;

        crossings = casesClassic[cubestate];

        while (-1 != crossings[i] ) {
                struct vector crossvertex;
                vertindex vertex_index;

                vertex_index = mcubes_get_edge_vertex(mc,x,y,z,crossings[i]);
                crossvertex = *(struct vector *)
                        growarr_get_nth(mc->vertex,vertex_index);
                newvert = vector_add(newvert,crossvertex);
                i++;
                
        }
        /* i = number of crossings */
        newvert = vector_scale(newvert,1.0/(float)i);
        mc->nvertices++;
#ifdef DEBUG
        {
                vertindex foo;
                foo=mcubes_add_vertex(mc,newvert.x,newvert.y,newvert.z);
                growarr_add(mc->alias_map,&foo);
#ifdef SICKENINGLY_VERBOSE
                fprintf(stderr,"Added central vertex %u at %d %d %d\n",
                        foo,x,y,z);
#endif /* SICKENINGLY_VERBOSE */
                return foo;
        }
#else
        return mcubes_add_vertex(mc,newvert.x,newvert.y,newvert.z);
#endif /* DEBUG */

        

}

/*
 * Test the interior of a cube
 * if s == 7, return true  if the interior is empty
 * if s ==-7, return false if the interior is empty
 */
int test_interior( int s, float *_cube, int _case, int _config, int _subconfig)
{
    float t, At=0, Bt=0, Ct=0, Dt=0, a, b ;
    char  test =  0 ;
    char  edge = -1 ; /* reference edge of the triangulation */

    switch( _case ) {
        case  4 :
        case 10 :
            a = ( _cube[4] - _cube[0] ) 
              * ( _cube[6] - _cube[2] ) 
              - ( _cube[7] - _cube[3] ) 
              * ( _cube[5] - _cube[1] ) ;
            b =  _cube[2] * ( _cube[4] - _cube[0] ) 
               + _cube[0] * ( _cube[6] - _cube[2] )
               - _cube[1] * ( _cube[7] - _cube[3] ) 
               - _cube[3] * ( _cube[5] - _cube[1] ) ;
            t = - b / (2*a) ;

            if( t<0 || t>1 ) return s>0 ;

            At = _cube[0] + ( _cube[4] - _cube[0] ) * t ;
            Bt = _cube[3] + ( _cube[7] - _cube[3] ) * t ;
            Ct = _cube[2] + ( _cube[6] - _cube[2] ) * t ;
            Dt = _cube[1] + ( _cube[5] - _cube[1] ) * t ;
        break ;

        case  6 :
        case  7 :
        case 12 :
        case 13 :
            switch( _case )
            {
            case  6 : edge = test6 [_config][2] ; break ;
            case  7 : edge = test7 [_config][4] ; break ;
            case 12 : edge = test12[_config][3] ; break ;
            case 13 : edge = tiling13_5_1[_config][_subconfig][0] ; break ;
            }
        switch( edge ) {
            case  0 :
                t  = _cube[0] / ( _cube[0] - _cube[1] ) ;
                At = 0 ;
                Bt = _cube[3] + ( _cube[2] - _cube[3] ) * t ;
                Ct = _cube[7] + ( _cube[6] - _cube[7] ) * t ;
                Dt = _cube[4] + ( _cube[5] - _cube[4] ) * t ;
                break ;
            case  1 :
                t  = _cube[1] / ( _cube[1] - _cube[2] ) ;
                At = 0 ;
                Bt = _cube[0] + ( _cube[3] - _cube[0] ) * t ;
                Ct = _cube[4] + ( _cube[7] - _cube[4] ) * t ;
                Dt = _cube[5] + ( _cube[6] - _cube[5] ) * t ;
                break ;
            case  2 :
                t  = _cube[2] / ( _cube[2] - _cube[3] ) ;
                At = 0 ;
                Bt = _cube[1] + ( _cube[0] - _cube[1] ) * t ;
                Ct = _cube[5] + ( _cube[4] - _cube[5] ) * t ;
                Dt = _cube[6] + ( _cube[7] - _cube[6] ) * t ;
                break ;
            case  3 :
                t  = _cube[3] / ( _cube[3] - _cube[0] ) ;
                At = 0 ;
                Bt = _cube[2] + ( _cube[1] - _cube[2] ) * t ;
                Ct = _cube[6] + ( _cube[5] - _cube[6] ) * t ;
                Dt = _cube[7] + ( _cube[4] - _cube[7] ) * t ;
                break ;
            case  4 :
                t  = _cube[4] / ( _cube[4] - _cube[5] ) ;
                At = 0 ;
                Bt = _cube[7] + ( _cube[6] - _cube[7] ) * t ;
                Ct = _cube[3] + ( _cube[2] - _cube[3] ) * t ;
                Dt = _cube[0] + ( _cube[1] - _cube[0] ) * t ;
                break ;
            case  5 :
                t  = _cube[5] / ( _cube[5] - _cube[6] ) ;
                At = 0 ;
                Bt = _cube[4] + ( _cube[7] - _cube[4] ) * t ;
                Ct = _cube[0] + ( _cube[3] - _cube[0] ) * t ;
                Dt = _cube[1] + ( _cube[2] - _cube[1] ) * t ;
                break ;
            case  6 :
                t  = _cube[6] / ( _cube[6] - _cube[7] ) ;
                At = 0 ;
                Bt = _cube[5] + ( _cube[4] - _cube[5] ) * t ;
                Ct = _cube[1] + ( _cube[0] - _cube[1] ) * t ;
                Dt = _cube[2] + ( _cube[3] - _cube[2] ) * t ;
                break ;
            case  7 :
                t  = _cube[7] / ( _cube[7] - _cube[4] ) ;
                At = 0 ;
                Bt = _cube[6] + ( _cube[5] - _cube[6] ) * t ;
                Ct = _cube[2] + ( _cube[1] - _cube[2] ) * t ;
                Dt = _cube[3] + ( _cube[0] - _cube[3] ) * t ;
                break ;
            case  8 :
                t  = _cube[0] / ( _cube[0] - _cube[4] ) ;
                At = 0 ;
                Bt = _cube[3] + ( _cube[7] - _cube[3] ) * t ;
                Ct = _cube[2] + ( _cube[6] - _cube[2] ) * t ;
                Dt = _cube[1] + ( _cube[5] - _cube[1] ) * t ;
                break ;
            case  9 :
                t  = _cube[1] / ( _cube[1] - _cube[5] ) ;
                At = 0 ;
                Bt = _cube[0] + ( _cube[4] - _cube[0] ) * t ;
                Ct = _cube[3] + ( _cube[7] - _cube[3] ) * t ;
                Dt = _cube[2] + ( _cube[6] - _cube[2] ) * t ;
                break ;
            case 10 :
                t  = _cube[2] / ( _cube[2] - _cube[6] ) ;
                At = 0 ;
                Bt = _cube[1] + ( _cube[5] - _cube[1] ) * t ;
                Ct = _cube[0] + ( _cube[4] - _cube[0] ) * t ;
                Dt = _cube[3] + ( _cube[7] - _cube[3] ) * t ;
                break ;
            case 11 :
                t  = _cube[3] / ( _cube[3] - _cube[7] ) ;
                At = 0 ;
                Bt = _cube[2] + ( _cube[6] - _cube[2] ) * t ;
                Ct = _cube[1] + ( _cube[5] - _cube[1] ) * t ;
                Dt = _cube[0] + ( _cube[4] - _cube[0] ) * t ;
                break ;
            default :
                fprintf(stderr, "Invalid edge %d\n", edge ) ;
                print_cube(_cube) ;
                break ;
        }
        break ;

        default :
                fprintf(stderr,"Invalid ambiguous case %d\n", _case );
                print_cube(_cube) ;
                break ;

    } /* switch (_case) */

    if( At >= 0 ) test ++ ;
    if( Bt >= 0 ) test += 2 ;
    if( Ct >= 0 ) test += 4 ;
    if( Dt >= 0 ) test += 8 ;
    switch( test ) {
        case  0 : return s>0 ;
        case  1 : return s>0 ;
        case  2 : return s>0 ;
        case  3 : return s>0 ;
        case  4 : return s>0 ;
        case  5 : if( At * Ct <  Bt * Dt ) return s>0 ; break ;
        case  6 : return s>0 ;
        case  7 : return s<0 ;
        case  8 : return s>0 ;
        case  9 : return s>0 ;
        case 10 : if( At * Ct >= Bt * Dt ) return s>0 ; break ;
        case 11 : return s<0 ;
        case 12 : return s>0 ;
        case 13 : return s<0 ;
        case 14 : return s<0 ;
        case 15 : return s<0 ;
    }

    return s<0 ;
}


/* This routine corresponds roughly with (and is mostly taken from)
 * process_cube() in Lewiner's implementation.
 */

int mcubes_lewiner_cube(struct mcubes_state *mc, int x, int y, int z)
{
    float v[8];
    int cube=0,two_to_i=0,i=0;

    int _case,_config,_subconfig;
    float isoval;
    vertindex v12=0;

    isoval = mc->isoval;

    cube=0;

    v[0] = vol3df_wrap(mc->voxel,x  ,y  ,z  ) - isoval;
    v[1] = vol3df_wrap(mc->voxel,x+1,y  ,z  ) - isoval;
    v[2] = vol3df_wrap(mc->voxel,x+1,y+1,z  ) - isoval;
    v[3] = vol3df_wrap(mc->voxel,x  ,y+1,z  ) - isoval;
    v[4] = vol3df_wrap(mc->voxel,x  ,y  ,z+1) - isoval;
    v[5] = vol3df_wrap(mc->voxel,x+1,y  ,z+1) - isoval;
    v[6] = vol3df_wrap(mc->voxel,x+1,y+1,z+1) - isoval;
    v[7] = vol3df_wrap(mc->voxel,x  ,y+1,z+1) - isoval;

    /* Resolve any ambiguity about points very close to the surface. */
    for (i=0;i<8;i++) {
            if (fabs(v[i]) < FLT_EPSILON) { v[i] = FLT_EPSILON; }
    }

    /* For each point on the nearest-neighbour grid:
     * if the point is below the isosurface value,
     * set the corresponding bit in the cube.
     */

    two_to_i = 1;
    for (i=0;i<8;i++) {
            if (v[i] > 0) {
                    cube |= two_to_i;
            }
            two_to_i <<= 1;
    }

    /* cube is now a number between 0 and 255 inclusive,
     * representing the voxel state.
     */

    _case = cases[cube][0];
    _config = cases[cube][1];
    _subconfig = 0 ;

#ifdef SICKENINGLY_VERBOSE
        fprintf(stderr,"(%d %d %d) state %d case %d config %d\n",
                x,y,z,cube,_case,_config);
#endif /* SICKENINGLY_VERBOSE */

    switch( _case ) {

        case  0 :
            break ;

        case  1 :
            mcubes_triangulate_voxel(mc, 
                tiling1[_config], 1 ,v12,x,y,z) ;
            break ;

        case  2 :
            mcubes_triangulate_voxel(mc, 
                tiling2[_config], 2 ,v12,x,y,z) ;
            break ;

        case  3 :
            if( test_face( test3[_config],v) )
            mcubes_triangulate_voxel(mc, 
                tiling3_2[_config], 4 ,v12,x,y,z) ; /* 3.2 */
            else
            mcubes_triangulate_voxel(mc, 
                tiling3_1[_config], 2 ,v12,x,y,z) ; /* 3.1 */
            break ;

        case  4 :
            if( test_interior( test4[_config],v,_case,_config,_subconfig) )
            mcubes_triangulate_voxel(mc, 
                tiling4_1[_config], 2 ,v12,x,y,z) ; /* 4.1.1 */
            else
            mcubes_triangulate_voxel(mc, 
                tiling4_2[_config], 6 ,v12,x,y,z) ; /* 4.1.2 */
            break ;

        case  5 :
            mcubes_triangulate_voxel(mc, 
                tiling5[_config], 3 ,v12,x,y,z) ;
            break ;

        case  6 :
            if( test_face( test6[_config][0],v) )
            mcubes_triangulate_voxel(mc, 
                tiling6_2[_config], 5 ,v12,x,y,z) ; /* 6.2 */
            else
            {
            if( test_interior( test6[_config][1],v,_case,_config,_subconfig) )
            mcubes_triangulate_voxel(mc, 
                tiling6_1_1[_config], 3 ,v12,x,y,z) ; /* 6.1.1 */
            else
            mcubes_triangulate_voxel(mc, 
                tiling6_1_2[_config], 7 ,v12,x,y,z) ; /* 6.1.2 */
            }
            break ;

        case  7 :
            if( test_face( test7[_config][0] ,v) ) _subconfig +=  1 ;
            if( test_face( test7[_config][1] ,v) ) _subconfig +=  2 ;
            if( test_face( test7[_config][2] ,v) ) _subconfig +=  4 ;
            switch( _subconfig )
            {
        case 0 :
            mcubes_triangulate_voxel(mc, 
                tiling7_1[_config], 3 ,v12,x,y,z) ; break ;
        case 1 :
            mcubes_triangulate_voxel(mc, 
                tiling7_2[_config][0], 5 ,v12,x,y,z) ; break ;
        case 2 :
            mcubes_triangulate_voxel(mc, 
                tiling7_2[_config][1], 5 ,v12,x,y,z) ; break ;
        case 3 :
            v12 = mcubes_add_c_vertex(mc,x,y,z,cube) ;
            mcubes_triangulate_voxel(mc, 
                tiling7_3[_config][0], 9, v12,x,y,z) ; break ;
        case 4 :
            mcubes_triangulate_voxel(mc, 
                tiling7_2[_config][2], 5 ,v12,x,y,z) ; break ;
        case 5 :
            v12 = mcubes_add_c_vertex(mc,x,y,z,cube) ;
            mcubes_triangulate_voxel(mc, 
                tiling7_3[_config][1], 9, v12,x,y,z) ; break ;
        case 6 :
            v12 = mcubes_add_c_vertex(mc,x,y,z,cube) ;
            mcubes_triangulate_voxel(mc, 
                tiling7_3[_config][2], 9, v12,x,y,z) ; break ;
        case 7 :
            if( test_interior( test7[_config][3],v,_case,_config,_subconfig) )
            mcubes_triangulate_voxel(mc, 
                tiling7_4_2[_config], 9 ,v12,x,y,z) ;
            else
            mcubes_triangulate_voxel(mc, 
                tiling7_4_1[_config], 5 ,v12,x,y,z) ;
            break ;
            };
            break ;

        case  8 :
            mcubes_triangulate_voxel(mc, 
                tiling8[_config], 2 ,v12,x,y,z) ;
            break ;

        case  9 :
            mcubes_triangulate_voxel(mc, 
                tiling9[_config], 4 ,v12,x,y,z) ;
            break ;

        case 10 :
            if( test_face( test10[_config][0],v) )
            {
            if( test_face( test10[_config][1],v) )
            mcubes_triangulate_voxel(mc, 
                tiling10_1_1_[_config], 4 ,v12,x,y,z) ; /* 10.1.1 */
            else
            {
            v12 = mcubes_add_c_vertex(mc,x,y,z,cube) ;
            mcubes_triangulate_voxel(mc, 
                tiling10_2[_config], 8, v12,x,y,z) ; /* 10.2 */
            }
            }
            else
            {
            if( test_face( test10[_config][1],v) )
            {
            v12 = mcubes_add_c_vertex(mc,x,y,z,cube) ;
            mcubes_triangulate_voxel(mc, 
                tiling10_2_[_config], 8, v12,x,y,z) ; /* 10.2 */
            }
            else
            {
            if( test_interior( test10[_config][2],v,_case,_config,_subconfig) )
            mcubes_triangulate_voxel(mc, 
                tiling10_1_1[_config], 4 ,v12,x,y,z) ; /* 10.1.1 */
            else
            mcubes_triangulate_voxel(mc, 
                tiling10_1_2[_config], 8 ,v12,x,y,z) ; /* 10.1.2 */
            }
            }
            break ;

        case 11 :
            mcubes_triangulate_voxel(mc, 
                tiling11[_config], 4 ,v12,x,y,z) ;
            break ;

        case 12 :
            if( test_face( test12[_config][0],v) )
            {
            if( test_face( test12[_config][1],v) )
            mcubes_triangulate_voxel(mc, 
                tiling12_1_1_[_config], 4 ,v12,x,y,z) ; /* 12.1.1 */
            else
            {
            v12 = mcubes_add_c_vertex(mc,x,y,z,cube) ;
            mcubes_triangulate_voxel(mc, 
                tiling12_2[_config], 8, v12,x,y,z) ; /* 12.2 */
            }
            }
            else
            {
            if( test_face( test12[_config][1],v) )
            {
            v12 = mcubes_add_c_vertex(mc,x,y,z,cube) ;
            mcubes_triangulate_voxel(mc, 
                tiling12_2_[_config], 8, v12,x,y,z) ; /* 12.2 */
            }
            else
            {
            if( test_interior( test12[_config][2],v,_case,_config,_subconfig) )
            mcubes_triangulate_voxel(mc, 
                tiling12_1_1[_config], 4 ,v12,x,y,z) ; /* 12.1.1 */
            else
            mcubes_triangulate_voxel(mc, 
                tiling12_1_2[_config], 8 ,v12,x,y,z) ; /* 12.1.2 */
            }
            }
            break ;

        case 13 :
            if( test_face( test13[_config][0] ,v) ) _subconfig +=  1 ;
            if( test_face( test13[_config][1] ,v) ) _subconfig +=  2 ;
            if( test_face( test13[_config][2] ,v) ) _subconfig +=  4 ;
            if( test_face( test13[_config][3] ,v) ) _subconfig +=  8 ;
            if( test_face( test13[_config][4] ,v) ) _subconfig += 16 ;
            if( test_face( test13[_config][5] ,v) ) _subconfig += 32 ;
            switch( subconfig13[_subconfig] )
            {
        case 0 :/* 13.1 */
            mcubes_triangulate_voxel(mc, 
                tiling13_1[_config], 4 ,v12,x,y,z) ; break ;

        case 1 :/* 13.2 */
            mcubes_triangulate_voxel(mc, 
                tiling13_2[_config][0], 6 ,v12,x,y,z) ; break ;
        case 2 :/* 13.2 */
            mcubes_triangulate_voxel(mc, 
                tiling13_2[_config][1], 6 ,v12,x,y,z) ; break ;
        case 3 :/* 13.2 */
            mcubes_triangulate_voxel(mc, 
                tiling13_2[_config][2], 6 ,v12,x,y,z) ; break ;
        case 4 :/* 13.2 */
            mcubes_triangulate_voxel(mc, 
                tiling13_2[_config][3], 6 ,v12,x,y,z) ; break ;
        case 5 :/* 13.2 */
            mcubes_triangulate_voxel(mc, 
                tiling13_2[_config][4], 6 ,v12,x,y,z) ; break ;
        case 6 :/* 13.2 */
            mcubes_triangulate_voxel(mc, 
                tiling13_2[_config][5], 6 ,v12,x,y,z) ; break ;

        case 7 :/* 13.3 */
            v12 = mcubes_add_c_vertex(mc,x,y,z,cube) ;
            mcubes_triangulate_voxel(mc, 
                tiling13_3[_config][0], 10, v12,x,y,z) ; break ;
        case 8 :/* 13.3 */
            v12 = mcubes_add_c_vertex(mc,x,y,z,cube) ;
            mcubes_triangulate_voxel(mc, 
                tiling13_3[_config][1], 10, v12,x,y,z) ; break ;
        case 9 :/* 13.3 */
            v12 = mcubes_add_c_vertex(mc,x,y,z,cube) ;
            mcubes_triangulate_voxel(mc, 
                tiling13_3[_config][2], 10, v12,x,y,z) ; break ;
        case 10 :/* 13.3 */
            v12 = mcubes_add_c_vertex(mc,x,y,z,cube) ;
            mcubes_triangulate_voxel(mc, 
                tiling13_3[_config][3], 10, v12,x,y,z) ; break ;
        case 11 :/* 13.3 */
            v12 = mcubes_add_c_vertex(mc,x,y,z,cube) ;
            mcubes_triangulate_voxel(mc, 
                tiling13_3[_config][4], 10, v12,x,y,z) ; break ;
        case 12 :/* 13.3 */
            v12 = mcubes_add_c_vertex(mc,x,y,z,cube) ;
            mcubes_triangulate_voxel(mc, 
                tiling13_3[_config][5], 10, v12,x,y,z) ; break ;
        case 13 :/* 13.3 */
            v12 = mcubes_add_c_vertex(mc,x,y,z,cube) ;
            mcubes_triangulate_voxel(mc, 
                tiling13_3[_config][6], 10, v12,x,y,z) ; break ;
        case 14 :/* 13.3 */
            v12 = mcubes_add_c_vertex(mc,x,y,z,cube) ;
            mcubes_triangulate_voxel(mc, 
                tiling13_3[_config][7], 10, v12,x,y,z) ; break ;
        case 15 :/* 13.3 */
            v12 = mcubes_add_c_vertex(mc,x,y,z,cube) ;
            mcubes_triangulate_voxel(mc, 
                tiling13_3[_config][8], 10, v12,x,y,z) ; break ;
        case 16 :/* 13.3 */
            v12 = mcubes_add_c_vertex(mc,x,y,z,cube) ;
            mcubes_triangulate_voxel(mc, 
                tiling13_3[_config][9], 10, v12,x,y,z) ; break ;
        case 17 :/* 13.3 */
            v12 = mcubes_add_c_vertex(mc,x,y,z,cube) ;
            mcubes_triangulate_voxel(mc, 
                tiling13_3[_config][10], 10, v12,x,y,z) ; break ;
        case 18 :/* 13.3 */
            v12 = mcubes_add_c_vertex(mc,x,y,z,cube) ;
            mcubes_triangulate_voxel(mc, 
                tiling13_3[_config][11], 10, v12,x,y,z) ; break ;

        case 19 :/* 13.4 */
            v12 = mcubes_add_c_vertex(mc,x,y,z,cube) ;
            mcubes_triangulate_voxel(mc, 
                tiling13_4[_config][0], 12, v12,x,y,z) ; break ;
        case 20 :/* 13.4 */
            v12 = mcubes_add_c_vertex(mc,x,y,z,cube) ;
            mcubes_triangulate_voxel(mc, 
                tiling13_4[_config][1], 12, v12,x,y,z) ; break ;
        case 21 :/* 13.4 */
            v12 = mcubes_add_c_vertex(mc,x,y,z,cube) ;
            mcubes_triangulate_voxel(mc, 
                tiling13_4[_config][2], 12, v12,x,y,z) ; break ;
        case 22 :/* 13.4 */
            v12 = mcubes_add_c_vertex(mc,x,y,z,cube) ;
            mcubes_triangulate_voxel(mc, 
                tiling13_4[_config][3], 12, v12,x,y,z) ; break ;

        case 23 :/* 13.5 */
            _subconfig = 0 ;
            if( test_interior( test13[_config][6] ,v,_case,_config,_subconfig) )
            mcubes_triangulate_voxel(mc, 
                tiling13_5_1[_config][0], 6 ,v12,x,y,z) ;
            else
            mcubes_triangulate_voxel(mc, 
                tiling13_5_2[_config][0], 10 ,v12,x,y,z) ;
            break ;
        case 24 :/* 13.5 */
            _subconfig = 1 ;
            if( test_interior( test13[_config][6] ,v,_case,_config,_subconfig) )
            mcubes_triangulate_voxel(mc, 
                tiling13_5_1[_config][1], 6 ,v12,x,y,z) ;
            else
            mcubes_triangulate_voxel(mc, 
                tiling13_5_2[_config][1], 10 ,v12,x,y,z) ;
            break ;
        case 25 :/* 13.5 */
            _subconfig = 2 ;
            if( test_interior( test13[_config][6] ,v,_case,_config,_subconfig) )
            mcubes_triangulate_voxel(mc, 
                tiling13_5_1[_config][2], 6 ,v12,x,y,z) ;
            else
            mcubes_triangulate_voxel(mc, 
                tiling13_5_2[_config][2], 10 ,v12,x,y,z) ;
            break ;
        case 26 :/* 13.5 */
            _subconfig = 3 ;
            if( test_interior( test13[_config][6] ,v,_case,_config,_subconfig) )
            mcubes_triangulate_voxel(mc, 
                tiling13_5_1[_config][3], 6 ,v12,x,y,z) ;
            else
            mcubes_triangulate_voxel(mc, 
                tiling13_5_2[_config][3], 10 ,v12,x,y,z) ;
            break ;

        case 27 :/* 13.3 */
            v12 = mcubes_add_c_vertex(mc,x,y,z,cube) ;
            mcubes_triangulate_voxel(mc, 
                tiling13_3_[_config][0], 10, v12,x,y,z) ; break ;
        case 28 :/* 13.3 */
            v12 = mcubes_add_c_vertex(mc,x,y,z,cube) ;
            mcubes_triangulate_voxel(mc, 
                tiling13_3_[_config][1], 10, v12,x,y,z) ; break ;
        case 29 :/* 13.3 */
            v12 = mcubes_add_c_vertex(mc,x,y,z,cube) ;
            mcubes_triangulate_voxel(mc, 
                tiling13_3_[_config][2], 10, v12,x,y,z) ; break ;
        case 30 :/* 13.3 */
            v12 = mcubes_add_c_vertex(mc,x,y,z,cube) ;
            mcubes_triangulate_voxel(mc, 
                tiling13_3_[_config][3], 10, v12,x,y,z) ; break ;
        case 31 :/* 13.3 */
            v12 = mcubes_add_c_vertex(mc,x,y,z,cube) ;
            mcubes_triangulate_voxel(mc, 
                tiling13_3_[_config][4], 10, v12,x,y,z) ; break ;
        case 32 :/* 13.3 */
            v12 = mcubes_add_c_vertex(mc,x,y,z,cube) ;
            mcubes_triangulate_voxel(mc, 
                tiling13_3_[_config][5], 10, v12,x,y,z) ; break ;
        case 33 :/* 13.3 */
            v12 = mcubes_add_c_vertex(mc,x,y,z,cube) ;
            mcubes_triangulate_voxel(mc, 
                tiling13_3_[_config][6], 10, v12,x,y,z) ; break ;
        case 34 :/* 13.3 */
            v12 = mcubes_add_c_vertex(mc,x,y,z,cube) ;
            mcubes_triangulate_voxel(mc, 
                tiling13_3_[_config][7], 10, v12,x,y,z) ; break ;
        case 35 :/* 13.3 */
            v12 = mcubes_add_c_vertex(mc,x,y,z,cube) ;
            mcubes_triangulate_voxel(mc, 
                tiling13_3_[_config][8], 10, v12,x,y,z) ; break ;
        case 36 :/* 13.3 */
            v12 = mcubes_add_c_vertex(mc,x,y,z,cube) ;
            mcubes_triangulate_voxel(mc, 
                tiling13_3_[_config][9], 10, v12,x,y,z) ; break ;
        case 37 :/* 13.3 */
            v12 = mcubes_add_c_vertex(mc,x,y,z,cube) ;
            mcubes_triangulate_voxel(mc, 
                tiling13_3_[_config][10], 10, v12,x,y,z) ; break ;
        case 38 :/* 13.3 */
            v12 = mcubes_add_c_vertex(mc,x,y,z,cube) ;
            mcubes_triangulate_voxel(mc, 
                tiling13_3_[_config][11], 10, v12,x,y,z) ; break ;

        case 39 :/* 13.2 */
            mcubes_triangulate_voxel(mc, 
                tiling13_2_[_config][0], 6 ,v12,x,y,z) ; break ;
        case 40 :/* 13.2 */
            mcubes_triangulate_voxel(mc, 
                tiling13_2_[_config][1], 6 ,v12,x,y,z) ; break ;
        case 41 :/* 13.2 */
            mcubes_triangulate_voxel(mc, 
                tiling13_2_[_config][2], 6 ,v12,x,y,z) ; break ;
        case 42 :/* 13.2 */
            mcubes_triangulate_voxel(mc, 
                tiling13_2_[_config][3], 6 ,v12,x,y,z) ; break ;
        case 43 :/* 13.2 */
            mcubes_triangulate_voxel(mc, 
                tiling13_2_[_config][4], 6 ,v12,x,y,z) ; break ;
        case 44 :/* 13.2 */
            mcubes_triangulate_voxel(mc, 
                tiling13_2_[_config][5], 6 ,v12,x,y,z) ; break ;

        case 45 :/* 13.1 */
            mcubes_triangulate_voxel(mc, 
                tiling13_1_[_config], 4 ,v12,x,y,z) ; break ;

            default :
            printf("Marching Cubes: Impossible case 13?\n" ) ;  print_cube(v) ;
            }
            break ;

        case 14 :
            mcubes_triangulate_voxel(mc, 
                tiling14[_config], 4 ,v12,x,y,z) ;
            break ;

    };
    return 0;
}


/* Triangulate isosurface using topologically-correct approach of
 * Lewiner et al, in turn based on the analysis of Chernyaev.
 */

int mcubes_lewiner(struct mcubes_state *mc)
{
        int x,y,z;
        float genus;


#ifdef SICKENINGLY_VERBOSE
        fprintf(stderr,"Triangulating...\n");
#endif

        if (mc->calc_area) {
                /* Zero the area counter */
                mc->total_area=0.0;
            if (mc->calc_mean_curvature) {
                mc->sum_H_area = mc->sum_H2_area=mc->sum_H4_area=0;
            }
            if (mc->calc_gaussian_curvature) {
                mc->sum_K_area = mc->sum_K2_area=0;
            }
        }
        if (mc->calc_sum_scalar) {
            mc->sum_scalar_area = 0.0;
            mc->sum_scalar2_area = 0.0;
        }


        for (z=0;z<mc->nz;z++) {
        for (y=0;y<mc->ny;y++) {
        for (x=0;x<mc->nx;x++) {
                mcubes_lewiner_cube(mc,x,y,z);
        }}} /* xyz */
        if (mc->fatal_error) {
                fprintf(stderr,"mcubes_lewiner: fatal error detected!\n");
                return -1;
        }
#ifdef SICKENINGLY_VERBOSE
        fprintf(stderr,"Triangulated.\n");
#endif

        /* Check that the genus is integral. If this is not the case,
         * then the mesh is not manifold.
         */

        if (0!=mcubes_genus(mc,&genus)) {
            fprintf(stderr,"mcubes_lewiner: mcubes_genus failed.\n");
            return -1;
        }

        if (fabs(genus-floor(genus)) > FLT_EPSILON) {
            fprintf(stderr,"mcubes_lewiner: NON-MANIFOLD MESH GENERATED!\n");
            fprintf(stderr,"Bailing out.\n");
            exit(-1); /* A bit drastic, but I'd rather know about it. */
        }

        if (mc->scalar_curvature) {
            /* Calculate verticial curvature statistics */

            int i,imax;
            float sum=0.0,sum2=0.0;
            float phimax,phimin;
            struct vol3df_stats *stats=NULL;
            fprintf(stderr,"doing curvatures..\n");fflush(stderr);

            if (NULL==(stats=malloc(sizeof(struct vol3df_stats)))) {
                perror("mcubes_lewiner: malloc(vol3df_stats)");
                return -1;
            }
            if (NULL==mc->scalar) {
                fprintf(stderr,"mcubes_lewiner: scalar_curvature but no scalar!\n");
                return -1;
            }

            phimax=phimin=*(float *)growarr_get_nth(mc->scalar,0);
            imax = growarr_nitems(mc->scalar);
            for (i=0;i<imax;i++) {
                float phi;
                phi = *(float *)growarr_get_nth(mc->scalar,i);
                sum += phi ; sum2 += phi*phi;
                if (phi>phimax) { phimax=phi; }
                if (phi<phimin) { phimin=phi; }
            }

            stats->max = phimax;
            stats->min = phimin;

            sum/=imax; sum2/=imax;
            stats->mean = sum;
            stats->variance = sum2 - sum*sum;
            mc->scalarstats = stats;

            fprintf(stderr,"Done curvatures.\n");fflush(stderr);
        }


        return 0;
}

/* Triangulate isosurface using classic Lorensen-Cline Marching Cubes */

int mcubes_classic(struct mcubes_state *mc)
{
        int x,y,z;
        float grid[8];
        int cube=0,two_to_i=0,i=0;
        int nt;

        float isoval;


#ifdef SICKENINGLY_VERBOSE
        fprintf(stderr,"Triangulating...\n");
#endif

        isoval = mc->isoval;

        if (mc->calc_area) {
                /* Zero the area counter */
                mc->total_area=0.0;
        }
        if (mc->calc_mean_curvature) {
                /* Zero the area counter */
            mc->sum_H_area=mc->sum_H2_area=mc->sum_H4_area=0.0;
        }
        if (mc->calc_gaussian_curvature) {
            mc->sum_K_area=mc->sum_K2_area=0.0;
        }
        if (mc->calc_sum_scalar) {
            mc->sum_scalar_area = 0.0;
            mc->sum_scalar2_area = 0.0;
        }

        for (z=0;z<mc->nz;z++) {
        for (y=0;y<mc->ny;y++) {
        for (x=0;x<mc->nx;x++) {
                cube=0; nt=0;


                grid[0] = vol3df_wrap(mc->voxel,x  ,y  ,z  ) - isoval;
                grid[1] = vol3df_wrap(mc->voxel,x+1,y  ,z  ) - isoval;
                grid[2] = vol3df_wrap(mc->voxel,x+1,y+1,z  ) - isoval;
                grid[3] = vol3df_wrap(mc->voxel,x  ,y+1,z  ) - isoval;
                grid[4] = vol3df_wrap(mc->voxel,x  ,y  ,z+1) - isoval;
                grid[5] = vol3df_wrap(mc->voxel,x+1,y  ,z+1) - isoval;
                grid[6] = vol3df_wrap(mc->voxel,x+1,y+1,z+1) - isoval;
                grid[7] = vol3df_wrap(mc->voxel,x  ,y+1,z+1) - isoval;

                /* Resolve any ambiguity about points which are
                 * very close to the surface.
                 */
                for (i=0;i<8;i++) {
                        if (fabs(grid[i]) < FLT_EPSILON)
                                { grid[i] = FLT_EPSILON; }
                }

                /* For each point on the nearest-neighbour grid:
                 * if the point is below the isosurface value,
                 * set the corresponding bit in the cube.
                 */

                two_to_i = 1;
                for (i=0;i<8;i++) {
                        if (grid[i] > 0.0) {
                                cube |= two_to_i;
                        }
                        two_to_i <<= 1;
                }

                /* Let's cop out and do classic MC. */

                while (casesClassic[cube][3*nt] != -1) { nt++; }
                mcubes_triangulate_voxel(mc,
                        casesClassic[cube],nt,0,x,y,z);


        }}} /* xyz */
        fprintf(stderr,"Classic MC done.\n");
        return 0;
}

/* 
 * Calculate a list of edges, and determine the consistency of a triangulation.
 */

struct edge { int count; vertindex hivert; };
#define MAXEDGES 20 /* max no of edges for a given vertex */

/* The edge table holds one entry for each vertex. Each entry contains
 * a list of the edges owned by that vertex.
 *
 * An edge belongs to its lowest-numbered vertex.
 *
 */

static void register_edge(struct edge *edgetable, vertindex hi, vertindex lo)
{
        int i=0;
        struct edge *entry=NULL;

        entry = edgetable+MAXEDGES*lo;

        for (i=0;i<MAXEDGES;i++) {
                if (hi == entry[i].hivert) {
                        /* This edge has already been created in the table.
                         * Increment its reference count.
                         */
                        entry[i].count++;
                        return;
                }
                if (-1==entry[i].count) {
                        /* This edge has not been created. Create it and
                         * set refcount=1.
                         */
                        entry[i].hivert = hi;
                        entry[i].count=1;
                        return;
                }
        }
        fprintf(stderr,"register_edge failed! hi=%u lo=%u\n",hi,lo);
        exit(-1);

}

static vertindex remap_vertindex(struct mcubes_state *mc, vertindex vindex)
{


#ifdef DEBUG
        vertindex vnew;
        vnew = *(vertindex *)growarr_get_nth(mc->alias_map,vindex);

        return vnew;
#else
        return vindex;
#endif

}

static struct tri remap_triangle(struct mcubes_state *mc, struct tri tri)
{
        struct tri newtri;

        newtri.a = remap_vertindex(mc,tri.a);
        newtri.b = remap_vertindex(mc,tri.b);
        newtri.c = remap_vertindex(mc,tri.c);
        if ( (newtri.a!=tri.a) ||(newtri.a!=tri.a) ||(newtri.a!=tri.a) ) {
            /*
            fprintf(stderr,"mapping:\n%d -> %d\n%d -> %d\n%d -> %d\n",
                    tri.a,newtri.a,
                    tri.b,newtri.b,
                    tri.c,newtri.c);
            fflush(stderr);
            */
        }

        return newtri;
}

static void register_triangle(struct mcubes_state *mc,
                struct edge *edgetable, struct tri tri)
{
        vertindex hi,lo;

        /* Determine if any of the triangle's vertices are in the halo
         * region. If so, change them to the corresponding non-halo
         * vertices.
         */

        tri = remap_triangle(mc,tri);


        /* A -> B */
        if (tri.a > tri.b)
                { hi = tri.a ; lo = tri.b; } else { hi = tri.b ; lo = tri.a; }
        register_edge(edgetable,hi,lo);

        /* C -> B */
        if (tri.c > tri.b)
                { hi = tri.c ; lo = tri.b; } else { hi = tri.b ; lo = tri.c; }
        register_edge(edgetable,hi,lo);

        /* C -> A */
        if (tri.c > tri.a)
                { hi = tri.c ; lo = tri.a; } else { hi = tri.a ; lo = tri.c; }
        register_edge(edgetable,hi,lo);

}

/* Check that each vertex is referenced by at least one triangle. */

void mcubes_check_vertex_usage(struct mcubes_state *mc)
{ 
        int nvertices,i;
        int *ulist=NULL;
        int ntris;

        fprintf(stderr,"Checking vertex usage...\n");

        nvertices = growarr_nitems(mc->vertex);
        ntris = growarr_nitems(mc->triangle);

        if (NULL==(ulist=malloc(sizeof(int)*nvertices))) 
                { perror("mcubes_check_vertex_usage: malloc"); return; }

        for (i=0;i<nvertices;i++) { ulist[i]=0; }

        for (i=0;i<ntris;i++) {
                struct tri t;

                t=*(struct tri *)growarr_get_nth(mc->triangle,i);
                if (t.a>nvertices) {
                        fprintf(stderr,"Triangle %d has bad vertex ref %u!\n",
                                        i,t.a); exit(-1); }
                ulist[t.a]=1;
                if (t.b>nvertices) {
                        fprintf(stderr,"Triangle %d has bad vertex ref %u!\n",
                                        i,t.b); exit(-1); }
                ulist[t.b]=1;
                if (t.c>nvertices) {
                        fprintf(stderr,"Triangle %d has bad vertex ref %u!\n",
                                        i,t.c); exit(-1); }
                ulist[t.c]=1;
        }

        for (i=0;i<nvertices;i++) {
                if (0==ulist[i]) {
                        fprintf(stderr,"*** Unused vertex %d\n",i);
                }
        }
        fprintf(stderr,"Done vertex usage check\n");

        free(ulist);
        return;

}

/* Delete any triangles whose centroids lie outside the inclusive
 * range (x:x+dx,y:y+dy,z:z+dz).
 */

int mcubes_triangles_subvol(struct mcubes_state *mc,
        int x0, int y0, int z0,
        int dx, int dy, int dz)
{
    int ntris;
    struct growarr *tri_new=NULL;
    struct tri t;
    int i;

    if (NULL==(tri_new=growarr_new(sizeof(struct tri),dz*dy*dz))) {
        fprintf(stderr,"mcubes_triangles_subvol: growarr_new failed.\n");
        return -1;
    }

    ntris = growarr_nitems(mc->triangle);

    for (i=0;i<ntris;i++) {
        struct vector va,vb,vc,centroid;

        t  = *(struct tri *)growarr_get_nth(mc->triangle,i);
        va = *(struct vector *)growarr_get_nth(mc->vertex,t.a);
        vb = *(struct vector *)growarr_get_nth(mc->vertex,t.b);
        vc = *(struct vector *)growarr_get_nth(mc->vertex,t.c);
        centroid = vector_add(va,vb);
        centroid = vector_add(centroid,vc);
        centroid = vector_scale(centroid,1.0/3.0);

        if (
                   (centroid.x >= x0) && (centroid.x <= (x0+dx) )
                && (centroid.y >= y0) && (centroid.y <= (y0+dy) )
                && (centroid.z >= z0) && (centroid.z <= (z0+dz) )
           ) {
            /* Inside clipping region */
            growarr_add(tri_new,&t);
            fprintf(stderr,"Clip added  %d %d %d\n",t.a,t.b,t.c);
            fflush(stderr);
        }
    }

    /* tri_new contains new triangle array; insert into mcubes state
     * and delete the old one.
     */

    growarr_destroy(mc->triangle);
    mc->triangle = tri_new;
    return 0;

}

void mcubes_audit(struct mcubes_state *mc)
{
        struct edge *edgetable=NULL;
        struct edge *entry=NULL;
        int nvertices,ntris;
        int i;
        int nedges=0;
        int nbad=0;

        mcubes_check_vertex_usage(mc);

        nvertices = growarr_nitems(mc->vertex);
        ntris = growarr_nitems(mc->triangle);

        if (NULL==(edgetable=malloc(sizeof(struct edge)*nvertices*MAXEDGES)))
                { perror("malloc"); exit(-1); }

        fprintf(stderr,"Initialising edge table\n");
        for (i=0;i<nvertices*MAXEDGES;i++) {
                edgetable[i].count = -1;
                edgetable[i].hivert = 0;
        }

        /* Run through all triangles.
         * Each triangle has three vertices: A, B, and C.
         * There is an edge between A and B, B and C, and A and C.
         * Register each of these edges.
         */
        fprintf(stderr,"Traversing triangle list of %d\n",ntris);

        for (i=0;i<ntris;i++) {
                struct tri *trip=NULL;
                struct tri mytri;
                trip = (struct tri *)growarr_get_nth(mc->triangle,i);
                if (NULL==trip) {
                        fprintf(stderr,"NULL trip %d\n",i); exit(-1);
                }
                mytri = *trip;
                /*
                    fprintf(stderr,"debug tri %d : %d %d %d\n",
                            i,mytri.a,mytri.b,mytri.c);
                    fflush(stderr);
                    */
                register_triangle(mc,edgetable,mytri);
        }

        /*
         * Now run through the edge table.  Count the number of unique
         * edges.  Also count (and list) all edges which have a
         * reference count which is not equal to two.  An edge should
         * join exactly two triangles, and therefore have a reference
         * count of two. A refcount of one would imply that an edge is
         * hanging off in space somewhere, and that there is a hole in
         * the triangulation, which is therefore not a manifold surface.
         */
        fprintf(stderr,"Traversing edge table\n");

        for (i=0;i<nvertices;i++) {
                int j;
                entry = edgetable + MAXEDGES*i;
                j=0;
                while (entry[j].count !=-1) {
                        nedges++;
                        if ((entry[j].count!=2)&&(entry[j].count!=4)) {
#ifdef DEBUG
                                struct vector v;
                                vertindex alias;
                                alias = *(vertindex *)growarr_get_nth(mc->alias_map,i);
                                v = *(struct vector *)growarr_get_nth(mc->vertex,i);
                              fprintf(stderr,"Vertex %d (%u) edge %d refcount %d {%02.08f,%02.08f,%02.08f}\n",
                                i,alias,j,entry[j].count,
                                v.x,v.y,v.z);
                              {
                                  int k;
                                  for (k=0;k<MAXEDGES;k++) {
                                      fprintf(stderr,
                                              "edge %d: hv %d count %d\n",
                                              k,
                                              entry[k].hivert,
                                              entry[k].count);
                                  }

                              }
#endif
                              nbad++;
                        }
                        j++;
                }
        }
        fprintf(stderr,"%d edges.\n",nedges);
        fprintf(stderr,"%d bad edges.\n",nbad);
        fprintf(stderr,"V=%d E=%d F=%d chi=%d\n",
                        nvertices,nedges,ntris,nvertices-nedges+ntris);

        /*
        mcubes_geometry_dump(mc,stderr);
        */


}

/* Unreference the volume and colour volume data, if they exist. */
/* This DOES NOT deallocate the volumes! */
int mcubes_release_input_data(struct mcubes_state *mc)
{
    mc->voxel = mc->colvol = NULL;
    mc->scalarstats = NULL;
    return 0;
}

/* Unreference and deallocate the volume and colour volume data. */
int mcubes_purge_input_data(struct mcubes_state *mc)
{
    if (NULL!=mc->voxel) { vol3d_destroy(mc->voxel); }
    if (NULL!=mc->colvol) { vol3d_destroy(mc->colvol); }
    if (NULL!=mc->scalarstats) { free(mc->scalarstats); }
    return mcubes_release_input_data(mc);
}

int mcubes_purge_cache(struct mcubes_state *mc)
{
    if (NULL!=mc->xvert) { vol3d_destroy(mc->xvert); mc->xvert=NULL; }
    if (NULL!=mc->yvert) { vol3d_destroy(mc->yvert); mc->yvert=NULL; }
    if (NULL!=mc->zvert) { vol3d_destroy(mc->zvert); mc->zvert=NULL; }
    return 0;
}

int mcubes_purge_geometry(struct mcubes_state *mc)
{
#ifdef DEBUG
    if (NULL!=mc->alias_map)
        { growarr_destroy(mc->alias_map); mc->alias_map=NULL; }
    if (NULL!=mc->triangle)
        { growarr_destroy(mc->triangle); mc->triangle=NULL; }
    if (NULL!=mc->vertex)
        { growarr_destroy(mc->vertex); mc->vertex=NULL; }
    if (NULL!=mc->normal)
        { growarr_destroy(mc->normal); mc->normal=NULL; }
    if (NULL!=mc->scalar)
        { growarr_destroy(mc->scalar); mc->scalar=NULL; }
#endif /* DEBUG */
    return 0;
}



int mcubes_write_surf(struct mcubes_state *mc, char *fname)
{
    FILE *fh=NULL;

    if (NULL==(fh=fopen(fname,"wb")))
        { perror("mcubes_write_surf: fopen"); return -1; }

    if (0!=mcubes_write_surf_fh(mc,fh)) {
        fprintf(stderr,"mcubes_write_surf: mcubes_write_surf_fh failed\n");
        fclose(fh);
        return -1;
    }

    fclose(fh);

    return 0;
}


/* Write an isosurface in "surf" format to the given filehandle. */
int mcubes_write_surf_fh(struct mcubes_state *mc, FILE *fh)
{

    XDR xdrs;
    uint32_t i;
    uint32_t nvertices,ntriangles,nnormals=0;

    xdrstdio_create(&xdrs,fh,XDR_ENCODE);

    /* Write magic file header */

    mcubes_xdr_uint32(&xdrs,SURF_BLOCKTYPE_FILE);

    /* Write no of vertices, triangles, normals */

    mcubes_xdr_uint32(&xdrs,nvertices=(uint32_t) growarr_nitems(mc->vertex));
    mcubes_xdr_uint32(&xdrs,ntriangles=(uint32_t) growarr_nitems(mc->triangle));
    if (NULL!=mc->normal) {
            mcubes_xdr_uint32(&xdrs,nnormals=(uint32_t) growarr_nitems(mc->normal));
    } else {
            mcubes_xdr_uint32(&xdrs,0);
    }

    /* Write number of unique vertices */
    mcubes_xdr_uint32(&xdrs,(uint32_t) mc->nvertices);

    /* Write vertex block */

            mcubes_xdr_uint32(&xdrs,SURF_BLOCKTYPE_VERTICES);

            /* Write vertices */

            for (i=0;i<nvertices;i++) {
                struct vector v;
                v = *(struct vector *)growarr_get_nth(mc->vertex,i);
                xdr_float(&xdrs,&(v.x));
                xdr_float(&xdrs,&(v.y));
                xdr_float(&xdrs,&(v.z));
            }

    /* Write triangle block */

            mcubes_xdr_uint32(&xdrs,SURF_BLOCKTYPE_TRIANGLES);

            /* Write triangles */

            for (i=0;i<ntriangles;i++) {
                    struct tri t;
                    t = *(struct tri *)growarr_get_nth(mc->triangle,i);
                    mcubes_xdr_uint32(&xdrs,t.a);
                    mcubes_xdr_uint32(&xdrs,t.b);
                    mcubes_xdr_uint32(&xdrs,t.c);
            }

    if (NULL!=mc->normal) {
    /* Write normal block */

            mcubes_xdr_uint32(&xdrs,SURF_BLOCKTYPE_NORMALS);

            /* Write number of normals */

            mcubes_xdr_uint32(&xdrs,(uint32_t) growarr_nitems(mc->normal));

            /* Write vertices */

            for (i=0;i<nnormals;i++) {
                struct vector v;
                v = *(struct vector *)growarr_get_nth(mc->normal,i);
                xdr_float(&xdrs,&(v.x));
                xdr_float(&xdrs,&(v.y));
                xdr_float(&xdrs,&(v.z));
            }
    }


    xdr_destroy(&xdrs);

    return 0;

}


static void triangle_area_centroid(
        struct mcubes_state *mc, struct vol3d *v,
        uint32_t trino,float *area,vector *centroid)
{
        struct tri t;
        vector vertex[3];

        t = *(struct tri *)growarr_get_nth(mc->triangle,trino);

        vertex[0] = *(vector *)growarr_get_nth(mc->vertex,t.a);
        vertex[1] = *(vector *)growarr_get_nth(mc->vertex,t.b);
        vertex[2] = *(vector *)growarr_get_nth(mc->vertex,t.c);

        *centroid = vector_scale(
                vector_add(vertex[0],vector_add(vertex[1],vertex[2])),
                1.0/3.0);

        *area = triangle_area(vertex[0],vertex[1],vertex[2]);

        return;

}

void triangle_propertyarea_centroid(
        struct mcubes_state *mc, struct vol3d *v,
        uint32_t trino,float *proparea,vector *centroid, int property)
{
    struct tri t;
    vector vertex[3];
    float area;
    float Pa,Pb,Pc,P;

    t = *(struct tri *)growarr_get_nth(mc->triangle,trino);

    vertex[0] = *(vector *)growarr_get_nth(mc->vertex,t.a);
    vertex[1] = *(vector *)growarr_get_nth(mc->vertex,t.b);
    vertex[2] = *(vector *)growarr_get_nth(mc->vertex,t.c);

    *centroid = vector_scale(
            vector_add(vertex[0],vector_add(vertex[1],vertex[2])),
            1.0/3.0);

    area = triangle_area(vertex[0],vertex[1],vertex[2]);

    switch(property) {
        case SURF_PROPERTY_SCALAR:
           Pa=vol3df_trilinear_point(v,vertex[0].x,vertex[0].y,vertex[0].z);
           Pb=vol3df_trilinear_point(v,vertex[1].x,vertex[1].y,vertex[1].z);
           Pc=vol3df_trilinear_point(v,vertex[2].x,vertex[2].y,vertex[2].z);
           break;
        case SURF_PROPERTY_H:
          Pa = vol3df_mean_curvature(v,vertex[0].x,vertex[0].y,vertex[0].z);
          Pb = vol3df_mean_curvature(v,vertex[1].x,vertex[1].y,vertex[1].z);
          Pc = vol3df_mean_curvature(v,vertex[2].x,vertex[2].y,vertex[2].z);
          break;
        case SURF_PROPERTY_H2:
          Pa = vol3df_mean_curvature(v,vertex[0].x,vertex[0].y,vertex[0].z);
          Pb = vol3df_mean_curvature(v,vertex[1].x,vertex[1].y,vertex[1].z);
          Pc = vol3df_mean_curvature(v,vertex[2].x,vertex[2].y,vertex[2].z);
          Pa=Pa*Pa; Pb = Pb*Pb; Pc = Pc*Pc;
          break;
        case SURF_PROPERTY_K:
          Pa = vol3df_gaussian_curvature(v,vertex[0].x,vertex[0].y,vertex[0].z);
          Pb = vol3df_gaussian_curvature(v,vertex[1].x,vertex[1].y,vertex[1].z);
          Pc = vol3df_gaussian_curvature(v,vertex[2].x,vertex[2].y,vertex[2].z);
          break;
        case SURF_PROPERTY_K2:
          Pa = vol3df_gaussian_curvature(v,vertex[0].x,vertex[0].y,vertex[0].z);
          Pb = vol3df_gaussian_curvature(v,vertex[1].x,vertex[1].y,vertex[1].z);
          Pc = vol3df_gaussian_curvature(v,vertex[2].x,vertex[2].y,vertex[2].z);
          Pa=Pa*Pa; Pb = Pb*Pb; Pc = Pc*Pc;
          break;
        default:
          fprintf(stderr,"triangle_propertyarea_centroid: bad property %d\n",
                  property);
          return;
    }

    P = (Pa+Pb+Pc)/3.0;

    *proparea = P*area;

    return;

}

static void triangle_H2area_centroid(
        struct mcubes_state *mc, struct vol3d *v,
        uint32_t trino,float *h2area,vector *centroid)
{
        struct tri t;
        vector vertex[3];
        float area;
        float Ha,Hb,Hc,H;

        t = *(struct tri *)growarr_get_nth(mc->triangle,trino);

        vertex[0] = *(vector *)growarr_get_nth(mc->vertex,t.a);
        vertex[1] = *(vector *)growarr_get_nth(mc->vertex,t.b);
        vertex[2] = *(vector *)growarr_get_nth(mc->vertex,t.c);

        *centroid = vector_scale(
                vector_add(vertex[0],vector_add(vertex[1],vertex[2])),
                1.0/3.0);

        area = triangle_area(vertex[0],vertex[1],vertex[2]);

        Ha = vol3df_mean_curvature(v,vertex[0].x,vertex[0].y,vertex[0].z);
        Hb = vol3df_mean_curvature(v,vertex[1].x,vertex[1].y,vertex[1].z);
        Hc = vol3df_mean_curvature(v,vertex[2].x,vertex[2].y,vertex[2].z);

        H = (Ha+Hb+Hc)/3.0;

        *h2area = H*H*area;

        return;

}


struct vol3d *vol3df_new_areadensity(struct vol3d *v,float isoval)
{
    struct vol3d *areavol=NULL;
    struct mcubes_state *mc=NULL;
    int i;
    uint32_t ntriangles;

    if (NULL==(mc=mcubes_new(v,isoval,NULL,0))) {
            fprintf(stderr,"mcubes_new failed.\n");
            return NULL;
    }

    if (0!=mcubes_calc_vertices(mc)) {
            fprintf(stderr,"mcubes_calc_vertices failed.\n");
            return NULL;
    }

    mcubes_lewiner(mc);

    /* Free vertex cache */

    mcubes_purge_cache(mc);

    if (NULL==(areavol=vol3df_new_zero(v->nx,v->ny,v->nz))) {
        fprintf(stderr,"vol3df_new_areadensity: vol3df_new_zero failed\n");
        mcubes_destroy(mc);
        return NULL;
    }

    /* For each triangle on the surface, calculate its area. Find its centroid,
     * and add the area to the value in the voxel containing the centroid.
     */

    ntriangles = growarr_nitems(mc->triangle);

    for (i=0;i<ntriangles;i++) {
        float area;
        int x,y,z;
        vector centroid;
        float *p;

        triangle_area_centroid(mc,v,i,&area,&centroid);
        x=(int)floor(centroid.x);
        y=(int)floor(centroid.y);
        z=(int)floor(centroid.z);

        p = vol3df_wrap_ptr(areavol,x,y,z);
        *p += fabs(area);

    }

    mcubes_destroy(mc);
    return areavol;

}


struct vol3d *vol3df_new_h2density(struct vol3d *v,float isoval)
{
    struct vol3d *areavol=NULL,*h2vol=NULL;
    struct mcubes_state *mc=NULL;
    int i;
    uint32_t ntriangles;
    int x,y,z;

    if (NULL==(mc=mcubes_new(v,isoval,NULL,0))) {
            fprintf(stderr,"mcubes_new failed.\n");
            return NULL;
    }

    if (0!=mcubes_calc_vertices(mc)) {
            fprintf(stderr,"mcubes_calc_vertices failed.\n");
            return NULL;
    }

    mcubes_lewiner(mc);

    /* Free vertex cache */

    mcubes_purge_cache(mc);

    /* Find the area density. */

    if (NULL==(areavol=vol3df_new_zero(v->nx,v->ny,v->nz))) {
        fprintf(stderr,"vol3df_new_areadensity: vol3df_new_zero failed\n");
        mcubes_destroy(mc);
        return NULL;
    }


    ntriangles = growarr_nitems(mc->triangle);

    for (i=0;i<ntriangles;i++) {
        float area;
        vector centroid;
        float *p;

        triangle_area_centroid(mc,v,i,&area,&centroid);
        x=(int)floor(centroid.x);
        y=(int)floor(centroid.y);
        z=(int)floor(centroid.z);

        p = vol3df_wrap_ptr(areavol,x,y,z);
        *p += fabs(area);

    }

    /* Find H^2 density */

    if (NULL==(h2vol=vol3df_new_zero(v->nx,v->ny,v->nz))) {
        fprintf(stderr,"vol3df_new_h2density: vol3df_new_zero failed\n");
        mcubes_destroy(mc);
        return NULL;
    }


    ntriangles = growarr_nitems(mc->triangle);

    for (i=0;i<ntriangles;i++) {
        float h2area;
        vector centroid;
        float *p;

        triangle_H2area_centroid(mc,v,i,&h2area,&centroid);
        x=(int)floor(centroid.x);
        y=(int)floor(centroid.y);
        z=(int)floor(centroid.z);

        p = vol3df_wrap_ptr(h2vol,x,y,z);
        *p += h2area;

    }

    mcubes_destroy(mc);

    /* Now normalize */

    for (z=0;z<v->nz;z++) {
    for (y=0;y<v->ny;y++) {
    for (x=0;x<v->nx;x++) {
        float area,h2area;
        area = vol3df_wrap(areavol,x,y,z);
        h2area = vol3df_wrap(h2vol,x,y,z);
        if (0.0<area) {
                *(vol3df_wrap_ptr(h2vol,x,y,z)) = h2area/area;
        } else {
                *(vol3df_wrap_ptr(h2vol,x,y,z)) = 0.0;
        }
    }}}

    vol3d_destroy(areavol);

    return h2vol;

}

