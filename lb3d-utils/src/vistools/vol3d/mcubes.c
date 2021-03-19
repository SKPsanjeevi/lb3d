#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <rpc/rpc.h>

#include "mcubes_tables.h"
#include "vol3df.h"


struct vertex vertex_add (struct vertex a, struct vertex b)
{
        struct vertex c;
        c.x = a.x + b.x; c.y = a.y + b.y; c.z = a.z + b.z;
        return c;
}

struct vertex vertex_sub (struct vertex a, struct vertex b)
{
        struct vertex c;
        c.x = a.x - b.x; c.y = a.y - b.y; c.z = a.z - b.z;
        return c;
}

struct vertex vertex_cross(struct vertex a, struct vertex b)
{
        struct vertex c;
        c.x = a.y * b.z - a.z * b.y;
        c.y = a.z * b.x - a.x * b.z;
        c.z = a.x * b.y - a.y * b.x;
        return c;
}

float vertex_modulus(struct vertex v)
{
        return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

/* Return the area of the triangle defined by three vertices */

float tri_area(struct vertex a, struct vertex b, struct vertex c)
{
        struct vertex v1,v2;

        v1 = vertex_sub(c,a);
        v2 = vertex_sub(b,a);

        return 0.5*vertex_modulus(vertex_cross(v1,v2));
}

struct vertex vertex_scale(struct vertex a, float s)
{
        struct vertex b;
        b.x = s * a.x;
        b.y = s * a.y;
        b.z = s * a.z;
        return b;
}

/* Take two vertices, the associated scalars, and an isovalue.
 * Calculate, by linear interpolation, the point on the line
 * at which the scalar is equal to the isovalue. Return this point.
 */
struct vertex vertex_interpolate(struct vertex ra, float phia, 
                struct vertex rb, float phib,
                float isoval)
{
        float lambda;
        struct vertex delta;

        lambda = (isoval - phia)/(phib - phia);

        delta = vertex_sub(rb,ra);
        delta = vertex_scale(delta,lambda);
        return vertex_add(ra,delta);


}


int mesh_destroy(struct mesh *m)
{
        if (NULL!=m) {
                if (NULL==m->vert) { free(m->vert); }
                if (NULL==m->tri)  { free(m->tri) ; }
                free(m);
        }
        return 0;
}

/* Expand the max number of vertices to nnew, or
 * use heuristics if nnew=0
 */
int mesh_grow_vert(struct mesh *m, int nnew)
{
        struct vertex *vnew=NULL;
        struct vertex *vcnew=NULL;
        struct vertex *vnnew=NULL;

        if (0==nnew) { nnew = m->maxverts * 2; }

        if (NULL==(vnew=realloc(m->vert,sizeof(struct vertex)*nnew)))
                { perror("malloc"); return -1;}

        if (NULL!=m->vertcol) {
               if (NULL==(vcnew=realloc(m->vertcol,sizeof(struct vertex)*nnew)))
                        { perror("malloc"); return -1;}
                m->vertcol = vcnew;
        }
        if (NULL!=m->vertnorm) {
               if (NULL==(vnnew=realloc(m->vertnorm,
                                               sizeof(struct vertex)*nnew)))
                        { perror("malloc"); return -1;}
                m->vertnorm = vnnew;
        }

        m->maxverts = nnew;
        m->vert = vnew; 
        return 0;
}



int mesh_grow_tri(struct mesh *m, int nnew)
{
        struct triangle *tnew=NULL;

        if (0==nnew) { nnew = m->maxtris * 2; }

        if (NULL==(tnew=realloc(m->tri,sizeof(struct triangle)*nnew)))
                { perror("malloc"); return -1;}

        m->maxtris = nnew;
        m->tri = tnew; return 0;
}

int mesh_add_coloured_vert(struct mesh *m, struct vertex v, struct vol3d *vol,struct vol3df_stats *mystats, struct vertex *vnormal)
{
        float r0,g0,b0,r1,g1,b1;
        float phi,r,g,b;
        float phimax = 1.0;
        float phimin = 0.0;

        phimax = mystats->max;
        phimin = mystats->min;

        if (m->maxverts == m->nextvert) {
                if (0!=mesh_grow_vert(m,0)) {
                        fprintf(stderr,"mesh_add_cvert: grow_vert failed\n");
                        return -1;
                }
        }

        m->vert[m->nextvert].x = v.x;
        m->vert[m->nextvert].y = v.y;
        m->vert[m->nextvert].z = v.z;

        if ((NULL!=m->vertnorm) && (NULL!=vnormal)) {
                m->vertnorm[m->nextvert].x = vnormal->x;
                m->vertnorm[m->nextvert].y = vnormal->y;
                m->vertnorm[m->nextvert].z = vnormal->z;
        }

        phi = vol3df_trilinear_point(vol,v.x,v.y,v.z);

        phi = (phi-phimin)/(phimax-phimin);

        phi = phi*phi*phi;

        r0=1.0; g0=0.0; b0=0.0;

        r1=0.0; g1=1.0; b1=1.0;

        r = r0 + phi*(r1-r0);
        g = g0 + phi*(g1-g0);
        b = b0 + phi*(b1-b0);


        m->vertcol[m->nextvert].x = r;
        m->vertcol[m->nextvert].y = g;
        m->vertcol[m->nextvert].z = b;


        return m->nextvert++;
}

int mesh_add_vert(struct mesh *m, struct vertex v, struct vertex *vnormal)
{
        if (m->maxverts == m->nextvert) {
                if (0!=mesh_grow_vert(m,0)) {
                        fprintf(stderr,"mesh_add_vert: grow_vert failed\n");
                        return -1;
                }
        }

        m->vert[m->nextvert].x = v.x;
        m->vert[m->nextvert].y = v.y;
        m->vert[m->nextvert].z = v.z;

        if ((NULL!=m->vertnorm) && (NULL!=vnormal)) {
                m->vertnorm[m->nextvert].x = vnormal->x;
                m->vertnorm[m->nextvert].y = vnormal->y;
                m->vertnorm[m->nextvert].z = vnormal->z;
        }

        return m->nextvert++;
}


int mesh_add_tri(struct mesh *m, struct triangle t)
{
        if (m->maxtris == m->nexttri) {
                if (0!=mesh_grow_tri(m,0)) {
                        fprintf(stderr,"mesh_add_tri: grow_tri failed\n");
                        return -1;
                }
        }

        m->tri[m->nexttri++] = t;

        return 0;
}

struct mesh *mesh_new()
{
        struct mesh *m=NULL;
        const int nstart=1; /* Initial no of triangles */

        if (NULL==(m=malloc(sizeof(struct mesh))))
                { perror("malloc"); return NULL; }

        m->vert=NULL; m->tri=NULL; m->vertcol=NULL;
        m->nextvert = m->nexttri = 0;
        m->maxverts=3*nstart; m->maxtris=nstart;

        if (NULL==(m->vert=malloc(m->maxverts*sizeof(struct vertex))))
                { perror("malloc"); mesh_destroy(m); return NULL; }
        if (NULL==(m->tri=malloc(m->maxtris*sizeof(struct triangle))))
                { perror("malloc"); mesh_destroy(m); return NULL; }

        return m;

}

int mesh_alloc_vcol(struct mesh *m)
{
        if (NULL==(m->vertcol=malloc(m->maxverts*sizeof(struct vertex))))
                { perror("malloc"); return -1; }
        return 0;
}

int mesh_alloc_vnorm(struct mesh *m)
{
        if (NULL==(m->vertnorm=malloc(m->maxverts*sizeof(struct vertex))))
                { perror("malloc"); return -1; }
        return 0;
}

struct mesh *vol3df_marching_cubes(struct vol3d *v, float isoval, int calc_vnormals)
{
        return vol3df_marching_cubes_full(v,isoval,NULL,calc_vnormals);
}

struct mesh *vol3df_marching_cubes_full(struct vol3d *v, float isoval,
                struct vol3d *vcol, int calc_vnormals)
{
        struct mesh *m=NULL;
        int x,y,z;
        float grid[8]; /* neighbourhood of point */
        int i,two_to_i;
        struct vol3df_stats *stats=NULL;
        struct vertex vnormal;
        struct vertex *vnp=NULL;


        if (NULL==(m=mesh_new()))
                { fprintf(stderr,"mesh_new() failed\n"); return NULL; }

        if (NULL!=vcol) {
                if (0!=mesh_alloc_vcol(m)) {
                        fprintf(stderr,"unable to allocate vertex colours");
                        mesh_destroy(m);
                        return NULL;
                }
                stats = vol3df_getstats(vcol);
        }

        if (calc_vnormals) {
                if (0!=mesh_alloc_vnorm(m)) {
                        fprintf(stderr,"unable to allocate vertex normals");
                        mesh_destroy(m);
                        return NULL;
                }
                vnp=&vnormal;
        }

        for (z=0;z<v->nz;z++) {
        for (y=0;y<v->ny;y++) {
        for (x=0;x<v->nx;x++) {
                int index=0;
                struct vertex r;

                int tri_vert[12];

                r.x = x; r.y=y; r.z=z;

                grid[0] = vol3df_wrap(v,x  ,y  ,z  );
                grid[1] = vol3df_wrap(v,x+1,y  ,z  );
                grid[2] = vol3df_wrap(v,x+1,y+1,z  );
                grid[3] = vol3df_wrap(v,x  ,y+1,z  );
                grid[4] = vol3df_wrap(v,x  ,y  ,z+1);
                grid[5] = vol3df_wrap(v,x+1,y  ,z+1);
                grid[6] = vol3df_wrap(v,x+1,y+1,z+1);
                grid[7] = vol3df_wrap(v,x  ,y+1,z+1);


                /* For each point on the nearest-neighbour grid:
                 * if the point is below the isosurface value,
                 * set the corresponding bit in the index.
                 */
                two_to_i = 1;
                for (i=0;i<8;i++) {
                        if (grid[i] < isoval) {
                                index |= two_to_i;
                        }
                        two_to_i <<= 1;
                }

                /* For each edge on the nearest-neighbour grid:
                 * look in the edge table to see if there is a
                 * crossing on this edge. If so, calculate its
                 * coordinates by interpolation, and store them
                 * in the corresponding entry in the tri_vert table
                 */
                two_to_i = 1;

                for (i=0;i<12;i++) {
                        if (edgeTable[index] & two_to_i) {
                                /* There's a crossing on ith edge */
                                struct vertex ra,rb,rc;
                                int a,b; float phia, phib;

                                /* a and b are vertices joined by ith edge */
                                a = cube_edge[i][0];
                                b = cube_edge[i][1];

                                /* ra, rb are coords of a, b, etc */
                                ra = vertex_add(r,cube_vert[a]); phia = grid[a];
                                rb = vertex_add(r,cube_vert[b]); phib = grid[b];


                                rc = vertex_interpolate(ra,phia,rb,phib,isoval);

                                if (calc_vnormals) {
                                        vol3df_trilinear_gradient(v,
                                                        rc.x,rc.y,rc.z,
                                                        &vnormal.x,
                                                        &vnormal.y,
                                                        &vnormal.z);
                                        vnormal.x = - vnormal.x;
                                        vnormal.y = - vnormal.y;
                                        vnormal.z = - vnormal.z;
                                }

                                /* Add vertex to vertex list; get its index */
                                if (NULL==vcol) {
                                      if (0>(tri_vert[i]=mesh_add_vert(m,rc,vnp))) {
                                               fprintf(stderr,
                                                      "mesh_add_vert failed\n");
                                                mesh_destroy(m);
                                                return NULL;
                                      }
                                } else {
                                   if (0>(tri_vert[i]=mesh_add_coloured_vert(
                                                                   m,rc,vcol,stats,vnp)
                                      )) {
                                            fprintf(stderr,
                                                   "mesh_add_cvert failed\n");
                                             mesh_destroy(m);
                                             return NULL;
                                   }
                                }
                        }
                        two_to_i <<= 1;
                }

                /* For each edge i which contains a crossing, 
                 * tri_vert[i] now contains the index of the crossing vertex.
                 * Now create the triangles.
                 */

                for (i=0;triTable[index][i] != -1; i += 3) {
                        struct triangle t;

                        t.a = tri_vert[triTable[index][i]];
                        t.b = tri_vert[triTable[index][i+1]];
                        t.c = tri_vert[triTable[index][i+2]];
                        mesh_add_tri(m,t);
                }


        }}}

        if (NULL!=stats) { free(stats); }

        return m;
}



int mesh_write_vtkfh(struct mesh *m, FILE *f, int binary)
{
        int i,nverts,ntris;
        XDR xdrs;

        nverts = m->nextvert;
        ntris = m->nexttri;

        if (binary) { xdrstdio_create(&xdrs,f,XDR_ENCODE); }

        /* Write header */

        fprintf(f,"# vtk DataFile Version 3.0\n");
        fprintf(f,"Generated from marching cubes\n");
        fprintf(f, binary ? "BINARY\n" : "ASCII\n");
        fprintf(f,"DATASET POLYDATA\n");
        fprintf(f,"POINTS %d float\n",nverts);

        /* Write vertices */

        for (i=0;i<nverts;i++) {
                struct vertex v;

                v=m->vert[i];

                if (binary) {
                        xdr_float(&xdrs,&v.x);
                        xdr_float(&xdrs,&v.y);
                        xdr_float(&xdrs,&v.z);
                } else {
                        fprintf(f,"%f %f %f\n",v.x,v.y,v.z);
                }
        }

        /* Write triangles */
        fprintf(f,"\nPOLYGONS %d %d\n",ntris,ntris*4);
        for (i=0;i<ntris;i++) {
                struct triangle t;
                t = m->tri[i];
                if (binary) {
                        int n=3;
                        xdr_int(&xdrs,&n);
                        xdr_int(&xdrs,&t.a);
                        xdr_int(&xdrs,&t.b);
                        xdr_int(&xdrs,&t.c);
                } else {
                        fprintf(f,"3 %d %d %d\n", t.a,t.b,t.c);
                }
        }
        
        return 0;
}
