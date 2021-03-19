#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "vol3df.h"
#include "isosurface.h"
#include "growarr.h"

int main(int argc, char *argv[])
{
        char *infilename=NULL;
        struct mcubes_state *mc=NULL;
        float chi,g;

        if ((2!=argc)) {
                fprintf(stderr,"%s infile\n",argv[0]);
                return -1;
        }

        if (NULL==(mc=mcubes_new_fromsurf(argv[1]))) {
            fprintf(stderr,"mcubes_new_fromsurf failed\n");
            return -1;
        }

        if (0!=mcubes_chi(mc,&chi))
            { fprintf(stderr,"mcubes_chi() failed.\n"); exit(-1); }
        if (0!=mcubes_genus(mc,&g))
            { fprintf(stderr,"mcubes_genus() failed.\n"); exit(-1); }

        mcubes_info(mc,stdout);

        mcubes_destroy(mc);

        return 0;
}
