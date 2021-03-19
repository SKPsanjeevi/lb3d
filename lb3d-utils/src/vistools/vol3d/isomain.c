#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "vol3df.h"
#include "isosurface.h"

int main(int argc, char *argv[])
{
        char *infilename=NULL,*outfilename=NULL;
        struct vol3df_file_hints hints = vol3df_default_file_hints;
        struct vol3d *v=NULL;
        struct mcubes_state *mc=NULL;

        struct vol3d *vv=NULL;
        int x,y,z;

        if (2!=argc) {
                fprintf(stderr,"%s infile\n",argv[0]);
                return -1;
        }
        hints.filename = infilename=argv[1]; outfilename=argv[2];
        if (NULL==(v=vol3df_new_from_file_heuristic(infilename,&hints)))
        {
                fprintf(stderr,"Failed to read input file\n");
                return -1;
        }

        if (NULL==(mc=mcubes_new(v,0.0,NULL,1))) {
                fprintf(stderr,"mcubes_new failed.\n");
                return -1;
        }
        mcubes_count_curvature(mc); /* Also counts area */

        if (0!=mcubes_calc_vertices(mc)) {
                fprintf(stderr,"mcubes_calc_vertices failed.\n");
                return -1;
        }

        /* dodgy bits */

        vv=vol3df_new(mc->nx,mc->ny,mc->nz);

        for (z=0;z<mc->nz;z++) {
        for (y=0;y<mc->ny;y++) {
        for (x=0;x<mc->nx;x++) {
                vertindex vi;
                vi = *(vertindex *)vol3d_element(mc->xvert,x,y,z);
                *(float *)vol3d_element(vv,x,y,z)=(float)vi;
        }}}

        /*
        of = fopen("xvert.xdr","wb");
        vol3df_write_xdrfh(vv,of);
        fclose(of);

        for (z=0;z<mc->nz;z++) {
        for (y=0;y<mc->ny;y++) {
        for (x=0;x<mc->nx;x++) {
                vertindex vi;
                vi = *(vertindex *)vol3d_element(mc->yvert,x,y,z);
                *(float *)vol3d_element(vv,x,y,z)=(float)vi;
        }}}

        of = fopen("yvert.xdr","wb");
        vol3df_write_xdrfh(vv,of);
        fclose(of);
        for (z=0;z<mc->nz;z++) {
        for (y=0;y<mc->ny;y++) {
        for (x=0;x<mc->nx;x++) {
                vertindex vi;
                vi = *(vertindex *)vol3d_element(mc->zvert,x,y,z);
                *(float *)vol3d_element(vv,x,y,z)=(float)vi;
        }}}

        of = fopen("zvert.xdr","wb");
        vol3df_write_xdrfh(vv,of);
        fclose(of);
        */

        fflush(stderr);
        mcubes_lewiner(mc);
        mcubes_purge_cache(mc);
        mcubes_purge_input_data(mc); /* destroys v as well */

        fflush(stderr);
        mcubes_info(mc,stderr);
        fflush(stderr);
        mcubes_audit(mc); 
        fflush(stderr);

        /*
        of=fopen("geom.vtk","wb");
        mcubes_write_vtkfh(mc, of, 1);
        fclose(of);
        */

        return 0;

}
