#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "vol3df.h"
#include "isosurface.h"
#include "growarr.h"

/* Structure containing the straightforward-to-measure values on a surface:
 * area, H integrated over area, H^2 integrated, ... .
 * Average of H can be found from areaH/area, etc.
 */

struct surface_stats {
    float area;
    float phi;
    float H;
    float H2;
    float H4;
    float K;
    float K2;
};

/* Return a surface_stats structure giving the area and average H, H2, etc
 * for a single triangle on the surface.
 */
struct surface_stats triangle_stats(struct mcubes_state *mc, struct vol3d *v,
        uint32_t trino)
{
        struct tri t;
        int j;
        vector vertex[3];
        float H=0.0,K=0.0,phi=0.0;
        struct surface_stats stats;

        t = *(struct tri *)growarr_get_nth(mc->triangle,trino);

        vertex[0] = *(vector *)growarr_get_nth(mc->vertex,t.a);
        vertex[1] = *(vector *)growarr_get_nth(mc->vertex,t.b);
        vertex[2] = *(vector *)growarr_get_nth(mc->vertex,t.c);

        stats.area = triangle_area(vertex[0],vertex[1],vertex[2]);

        for (j=0;j<3;j++) {
            float myH,myK;

            phi+=vol3df_trilinear_point(v,vertex[j].x,vertex[j].y,vertex[j].z);        

            vol3df_curvatures(v,vertex[j].x,vertex[j].y,vertex[j].z,&myH,&myK);
            H += myH; K += myK;
        }

        phi/=3.0; H/=3.0; K/=3.0;

        stats.phi = phi;
        stats.H = H; stats.H2 = H*H; stats.H4 = H*H*H*H;
        stats.K = K; stats.K2 = K*K;

        return stats;

}

/* Return a surface-stats structure containing the area, and the
 * values of H, H^2 etc INTEGRATED over the triangle.
 * ie area*H, area*H^2, etc.
 */

struct surface_stats surfavg_stats(struct mcubes_state *mc, struct vol3d *v)
{

    struct surface_stats sumstats;
    uint32_t i,ntriangles;

    sumstats.area = sumstats.H = sumstats.H2 = sumstats.H4 = 0.0;
    sumstats.K = sumstats.K2 = sumstats.phi = 0.0;


    ntriangles = growarr_nitems(mc->triangle);

    for (i=0;i<ntriangles;i++) {
        struct surface_stats tristats;
        float area;

        tristats = triangle_stats(mc,v,i);
        area = tristats.area;
        sumstats.area += area;
        sumstats.H += area*tristats.H;
        sumstats.H2 += area*tristats.H2;
        sumstats.H4 += area*tristats.H4;
        sumstats.K += area*tristats.K;
        sumstats.K2 += area*tristats.K2;
        sumstats.phi += area*tristats.phi;
    }

    return sumstats;

}


float surfavg_area(struct mcubes_state *mc, struct vol3d *v)
{
    uint32_t i,ntriangles;
    float areasum = 0.0;

    ntriangles = growarr_nitems(mc->triangle);

    for (i=0;i<ntriangles;i++) {
        struct tri t;

        t = *(struct tri *)growarr_get_nth(mc->triangle,i);
        areasum += triangle_area(
                *(struct vector *)growarr_get_nth(mc->vertex,t.a),
                *(struct vector *)growarr_get_nth(mc->vertex,t.b),
                *(struct vector *)growarr_get_nth(mc->vertex,t.c)
        );

    }

    return (float)areasum;

}

int main(int argc, char *argv[])
{
    
    char *surffilename=NULL,*volfilename=NULL,*op=NULL;
    struct mcubes_state *mc=NULL;
    struct vol3d *v=NULL;
    struct vol3df_file_hints hints = vol3df_default_file_hints;
    float chi,g;
    int retval=-1;

    if ((4!=argc)) {
            fprintf(stderr,"%s surffile volfile op\n",argv[0]);
            return -1;
    }
    surffilename=argv[1];
    hints.filename = volfilename=argv[2];
    op=argv[3];


    if (NULL==(mc=mcubes_new_fromsurf(surffilename))) {
        fprintf(stderr,"mcubes_new_fromsurf failed\n");
        return -1;
    }

    if (0!=mcubes_chi(mc,&chi))
        { fprintf(stderr,"mcubes_chi() failed.\n"); exit(-1); }
    if (0!=mcubes_genus(mc,&g))
        { fprintf(stderr,"mcubes_genus() failed.\n"); exit(-1); }

    if (NULL==(v=vol3df_new_from_file_heuristic(volfilename,&hints)))
    {
            fprintf(stderr,"Failed to read volume file %s\n",
                    volfilename);
            return -1;
    }

    /* Both volume and surface loaded; now decide what to do with them.. */

    if (0==strcmp("area",op)) {
        fprintf(stderr,"area= %.12f\n",surfavg_area(mc,v));
        retval=0;
    } else if (0==strcmp("stats",op)) {
        struct surface_stats stats;

        stats = surfavg_stats(mc,v);
        fprintf(stdout,"(area,H,H2,H4,K,K2,phi)= "
                "%.12f %.12f %.12f %.12f %.12f %.12f %.12f\n",
                stats.area,
                stats.H/stats.area,
                stats.H2/stats.area,
                stats.H4/stats.area,
                stats.K/stats.area,
                stats.K2/stats.area,
                stats.phi/stats.area);
    } else if (0==strcmp("genus",op)) {
        struct surface_stats stats;
        float chi,g;
        if (0!=mcubes_chi(mc,&chi))
            { fprintf(stderr,"mcubes_chi() failed.\n"); exit(-1); }
        if (0!=mcubes_genus(mc,&g))
            { fprintf(stderr,"mcubes_genus() failed.\n"); exit(-1); }

        stats = surfavg_stats(mc,v);
        printf("(chi,g,H,H2,H4,area,N,G,G2)= "
                "%.12f %.12f %.12f %.12f %.12f %.12f %u %.12f %.12f\n",
                chi,g,
                stats.H/stats.area,
                stats.H2/stats.area,
                stats.H4/stats.area,
                stats.area,
                growarr_nitems(mc->triangle),
                stats.K/stats.area,
                stats.K2/stats.area);
    } else {
        fprintf(stderr,"Unknown operation %s\n",op);
        return -1;
    }

    return retval;





    mcubes_info(mc,stdout);

    mcubes_destroy(mc);

    return 0;
}
