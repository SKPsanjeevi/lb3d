#ifndef INCLUDED_MCUBES_H
#define INCLUDED_MCUBES_H

#include "vol3d.h"
#include "vol3df.h"

struct vertex { float x,y,z; };
struct triangle { int a,b,c; }; /* Indices into vertex list */
struct mesh {
        struct vertex *vert;
        struct vertex *vertcol;
        struct vertex *vertnorm;
        struct triangle *tri;
        int nextvert; /* Index of next vertex to be added */
        int maxverts; /* Max number of vertices */
        int nexttri,maxtris; /* Similarly for triangles */
};

struct vertex vertex_cross(struct vertex a, struct vertex b);
struct vertex vertex_add(struct vertex a, struct vertex b);
struct vertex vertex_sub(struct vertex a, struct vertex b);
struct vertex vertex_scale(struct vertex a, float s);
struct vertex vertex_interpolate(struct vertex ra, float phia, struct vertex rb, float phib, float isoval);
int mesh_destroy(struct mesh *m);
int mesh_grow_vert(struct mesh *m, int nnew);
int mesh_grow_tri(struct mesh *m, int nnew);
int mesh_add_vert(struct mesh *m, struct vertex v,struct vertex *vnorm);
int mesh_add_tri(struct mesh *m, struct triangle t);
struct mesh *mesh_new(void);
struct mesh *vol3df_marching_cubes(struct vol3d *v, float isoval,int calc_vnormals);
struct mesh *vol3df_marching_cubes_full(struct vol3d *v, float isoval,
                struct vol3d *vcol, int calc_vnormals);
int mesh_write_vtkfh(struct mesh *m, FILE *f, int binary);
int mesh_add_coloured_vert(struct mesh *m, struct vertex v, struct vol3d *vol,struct vol3df_stats *mystats,struct vertex *vnp);

#endif /* INCLUDED_MCUBES_H */
