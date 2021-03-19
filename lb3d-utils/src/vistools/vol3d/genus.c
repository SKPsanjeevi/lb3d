#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "vol3df.h"
#include "isosurface.h"
#include "growarr.h"

int main(int argc, char *argv[])
{
        char *infilename=NULL;
        struct vol3df_file_hints hints = vol3df_default_file_hints;
        struct vol3d *v=NULL;
        struct mcubes_state *mc=NULL;
        float H,H2,H4,area;
        float K,K2;
        float chi,g;
        float isoval=0.0;

        if ((2!=argc)&&(3!=argc)) {
                fprintf(stderr,"%s infile [isoval]\n",argv[0]);
                return -1;
        }
        hints.filename = infilename=argv[1]; 

        if (3==argc) { isoval = atof(argv[2]); }

        if (NULL==(v=vol3df_new_from_file_heuristic(infilename,&hints)))
        {
                fprintf(stderr,"Failed to read input file\n");
                return -1;
        }

        if (NULL==(mc=mcubes_new(v,isoval,NULL,0))) {
                fprintf(stderr,"mcubes_new failed.\n");
                return -1;
        }
        mcubes_count_curvature(mc); /* Also counts area */

        if (0!=mcubes_calc_vertices(mc)) {
                fprintf(stderr,"mcubes_calc_vertices failed.\n");
                return -1;
        }


        mcubes_lewiner(mc);
        if (0!=mcubes_chi(mc,&chi))
            { fprintf(stderr,"mcubes_chi() failed.\n"); exit(-1); }
        if (0!=mcubes_genus(mc,&g))
            { fprintf(stderr,"mcubes_genus() failed.\n"); exit(-1); }

        H = mcubes_average_curvature(mc);
        H2 = mcubes_meansquared_curvature(mc);
        H4 = mcubes_meanfourth_curvature(mc);
        K = mcubes_average_gaussian_curvature(mc);
        K2 = mcubes_meansquared_gaussian_curvature(mc);
        area = mcubes_area(mc);
        

        printf("(chi,g,H,H2,H4,area,N,K,K2)= %.12f %.12f %.12f %.12f %.12f %.12f %u %.12f %.12f\n",chi,g,H,H2,H4,area,growarr_nitems(mc->triangle),K,K2);
        mcubes_destroy(mc);
        return 0;
}
