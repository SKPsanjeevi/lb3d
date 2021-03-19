#include <stdio.h>
#include <stdlib.h>

#include "config.h"

#include "vol2du8.h"
#include "vol2df.h"
#include "vol2d.h"

static void truncfunc(void *el, void *data)
{
        float phi;

        phi = fabs(*(float *)el);
        if (phi > 10) { phi=1.0; }
        if (phi<=0) { phi=0.0; }
        *(float *)el = (phi);

}

int main (int argc, char *argv[])
{
        char *infilename=NULL,*outfilename=NULL;

        struct vol2d *in=NULL,*out=NULL,*pngdata=NULL,*outpng=NULL;

        if (3!=argc) {
                fprintf(stderr,"%s <in.png> <out.png>\n",argv[0]);
                return -1;
        }

        infilename=argv[1]; outfilename=argv[2];

        if (NULL==(pngdata=vol2du8_new_from_gspng(infilename))) {
                fprintf(stderr,"vol2du8_new_from_gspng failed\n");
                return -1;
        }

        if (NULL==(in=vol2df_new_from_2du8(pngdata))) {
                fprintf(stderr,"conversion from u8 failed\n");
                return -1;
        }

        if (NULL==(out=vol2df_new_fullspectrum(in))) {
                fprintf(stderr,"FFT failed\n");
                return -1;
        }
        printf("Input: ");
        vol2df_dumpstats(in,stdout);
        vol2d_destroy(in);

        vol2d_stream(out,truncfunc,NULL);

        if (NULL==(outpng=vol2du8_new_from_2df_normalize(out))) {
                fprintf(stderr,"Conversion to u8 failed\n");
                return -1;
        }

        vol2d_stream(out,truncfunc,NULL);
        printf("Output: ");
        vol2df_dumpstats(out,stdout);

        vol2df_dumpdist(out,32,stdout);
        vol2d_destroy(out);

        if (0!=vol2du8_write_png(outpng,outfilename)) {
                fprintf(stderr,"Output write failed\n");
                return -1;
        }

        vol2d_destroy(outpng);


        return 0;
}
