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
        float G,G2;
        float chi,g;
        float isoval=0.0;

        if ((3!=argc)&&(4!=argc)) {
                fprintf(stderr,"%s infile outfile [isoval]\n",argv[0]);
                return -1;
        }
        hints.filename = infilename=argv[1]; 

        if (4==argc) { isoval = atof(argv[3]); }

        if (NULL==(v=vol3df_new_from_file_heuristic(infilename,&hints)))
        {
                fprintf(stderr,"Failed to read input file\n");
                return -1;
        }

        if (NULL==(mc=mcubes_new(v,isoval,NULL,1))) {
                fprintf(stderr,"mcubes_new failed.\n");
                return -1;
        }

        if (0!=mcubes_calc_vertices(mc)) {
                fprintf(stderr,"mcubes_calc_vertices failed.\n");
                return -1;
        }


        mcubes_lewiner(mc);
        if (0!=mcubes_chi(mc,&chi))
            { fprintf(stderr,"mcubes_chi() failed.\n"); exit(-1); }
        if (0!=mcubes_genus(mc,&g))
            { fprintf(stderr,"mcubes_genus() failed.\n"); exit(-1); }

        mcubes_write_surf

        mcubes_destroy(mc);
        return 0;
}
