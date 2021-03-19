#ifndef  INCLUDED_ISOSURFACE_H 
#define  INCLUDED_ISOSURFACE_H 

#include "vector.h"

#undef DEBUG
#undef SICKENINGLY_VERBOSE

#define SURF_BLOCKTYPE_TRIANGLES  1
#define SURF_BLOCKTYPE_NORMALS    2
#define SURF_BLOCKTYPE_VERTICES   3
#define SURF_BLOCKTYPE_FILE 0x1234d00f

#define SURF_PROPERTY_SCALAR 1
#define SURF_PROPERTY_H 2
#define SURF_PROPERTY_H2 3
#define SURF_PROPERTY_K 4
#define SURF_PROPERTY_K2 5
#define SURF_PROPERTY_AREA 6





typedef uint32_t vertindex;
#define V_BADVAL 99999 /* bad vertex index (for debugging) */

struct tri {
        vertindex a,b,c;
};



struct mcubes_state {
        /* Universal data */
        float isoval;
        int nx,ny,nz;
        /* If the fatal_error flag is set, this indicates that an error has
         * occurred, the state of the mcubes object is probably invalid,
         * and any ongoing operation should be abandoned.
         */
        int fatal_error;

        /* Input data */
        struct vol3d *voxel; /* Voxel data to triangulate */
        struct vol3d *colvol; /* Voxel data for colours */
        struct vol3df_stats *scalarstats; /* Statistics for colour data */

        /* Geometry data */
#ifdef DEBUG
        struct growarr *alias_map;
#endif /* DEBUG */

        struct growarr *triangle,*vertex,*normal,*scalar;
        unsigned int nvertices; /* Number of unique (non-border) vertices */
        int scalar_curvature; /* If set, define scalar to be curvature */
        int scalar_gaussian; /* scalar curvature is Gaussian, not mean */

        /* Cache data */

        struct vol3d *xvert,*yvert,*zvert; /* Hold vertex information */

        /* Statistical data */

        int calc_area; /* Flag for area calculation */
        double total_area;
        int calc_mean_curvature; /* Flag for mean curvature calculation */
        int calc_gaussian_curvature; /* Flag for gaussian curv calculation */
        int calc_sum_scalar; /* Flag for total scalar*area calculation */
        double sum_H_area; /* Sum of (curvature * area) */
        double sum_K_area; /* Sum of (K * area) */
        double sum_K2_area; /* Sum of (K^2 * area) */
        double sum_H2_area; /* Sum of ( (curvature^2) * area ) */
        double sum_H4_area; /* Sum of ( (curvature^4) * area ) */
        double sum_scalar_area; /* Sum of (scalar * area) */
        double sum_scalar2_area; /* Sum of (scalar^2 * area) */

};


void mcubes_geometry_dump(struct mcubes_state *mc, FILE *fh);
void mcubes_info(struct mcubes_state *mc, FILE *fh);
int mcubes_gaussian_curvature(struct mcubes_state *mc);
void mcubes_destroy(struct mcubes_state *mc);
struct mcubes_state *mcubes_new(struct vol3d *voxel, float isoval, struct vol3d *colvol, int calc_normals);
struct mcubes_state *mcubes_new_minimal(void);
void mcubes_cache_vertex(struct mcubes_state *mc, struct vol3d *vertvol, int ix, int iy, int iz, vertindex v);
vertindex mcubes_add_vertex(struct mcubes_state *mc, float x, float y, float z);
int mcubes_calc_vertices(struct mcubes_state *mc);
void dump_verts(struct vol3d *v);
int mcubes_add_triangle(struct mcubes_state *mc, struct tri *t);
vertindex mcubes_get_edge_vertex(struct mcubes_state *mc, int x, int y, int z, int edgeno);
int mcubes_write_vtkfh(struct mcubes_state *mc, FILE *f, int binary);
vertindex mcubes_get_x_vert(struct mcubes_state *mc, int x, int y, int z);
vertindex mcubes_get_y_vert(struct mcubes_state *mc, int x, int y, int z);
vertindex mcubes_get_z_vert(struct mcubes_state *mc, int x, int y, int z);
int mcubes_triangulate_voxel(struct mcubes_state *mc, const signed char *trig, char n, vertindex v12, int i, int j, int k);
vertindex mcubes_add_c_vertex(struct mcubes_state *mc, int x, int y, int z, int cubestate);
int test_interior(int s, float *v, int m_case, int m_config, int m_subconfig);
int mcubes_lewiner_cube(struct mcubes_state *mc, int x, int y, int z);
int mcubes_lewiner(struct mcubes_state *mc);
int mcubes_classic(struct mcubes_state *mc);
void mcubes_audit(struct mcubes_state *mc);

int mcubes_release_input_data(struct mcubes_state *mc);
int mcubes_purge_input_data(struct mcubes_state *mc);
int mcubes_purge_cache(struct mcubes_state *mc);
int mcubes_purge_geometry(struct mcubes_state *mc);
int mcubes_triangles_subvol(struct mcubes_state *mc,
        int x0, int y0, int z0,
        int dx, int dy, int dz);
int mcubes_chi(struct mcubes_state *mc,float *chi);
int mcubes_genus(struct mcubes_state *mc, float *g);
float mcubes_average_curvature(struct mcubes_state *mc);
float mcubes_average_gaussian_curvature(struct mcubes_state *mc);
float mcubes_area(struct mcubes_state *mc);
float mcubes_meansquared_curvature(struct mcubes_state *mc);
float mcubes_meanfourth_curvature(struct mcubes_state *mc);
float mcubes_meansquared_gaussian_curvature(struct mcubes_state *mc);
void mcubes_count_curvature(struct mcubes_state *mc);
void mcubes_count_area(struct mcubes_state *mc);
int mcubes_scalar_curvature(struct mcubes_state *mc);
void mcubes_count_scalar(struct mcubes_state *mc);
float mcubes_sum_scalar_area(struct mcubes_state *mc);
float mcubes_sum_scalar2_area(struct mcubes_state *mc);
int mcubes_write_surf(struct mcubes_state *mc, char *fname);
int mcubes_write_surf_fh(struct mcubes_state *mc, FILE *fh);
struct mcubes_state *mcubes_new_fromsurf(char *filename);
struct mcubes_state *mcubes_new_fromsurf_fh(FILE *fh);
struct vol3d *vol3df_new_areadensity(struct vol3d *v,float isoval);
struct vol3d *vol3df_new_h2density(struct vol3d *v,float isoval);

#endif /* INCLUDED_ISOSURFACE_H */
