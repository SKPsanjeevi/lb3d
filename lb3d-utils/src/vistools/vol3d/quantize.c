#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <png.h>

#include "vol3d.h"
#include "vol3df.h"
#include "lut.h"
#include "lutWidget.h"

int vol3d_f_to_u8_ftable(struct vol3d *v,
        struct vol3d **quantp,struct u8_freqtable **tablep, float phimax,float phimin)
{
    struct vol3d *quantized=NULL;
    struct u8_freqtable *table=NULL;

    /* Pointers to float value and corresponding quantized byte */
    uint8_t *quantel=NULL;
    float *phip=NULL;
    float scale=0.0;

    unsigned long n,nmax;

    if (NULL==(quantized=vol3d_new(v->nx,v->ny,v->nz,VOL_TYPE_UINT8))) {
        fprintf(stderr,"vol3d_f_to_u8_ftable: vol3d_new failed\n");
        goto bailout;
    }

    if (NULL==(table=malloc(sizeof(struct u8_freqtable)))) {
        perror("vol3d_f_to_u8_ftable: malloc");
        goto bailout;
    }


    nmax = v->nx*v->ny*v->nz;
    phip = (float *)v->data;
    quantel = (uint8_t *)quantized->data;

    scale = 255.0/(phimax-phimin);
    for (n=0;n<nmax;n++) {
        float phi;
        uint8_t quant;

        phi = *phip++;
        phi = ( phi-phimin)*scale;
        if (phi<0.0) { phi=0.0; }
        if (phi>255.0) { phi=255.0; }
        quant = (uint8_t) phi;

        *quantel++ = quant;
        table->freq[quant]++;
    }

    /* All done! */

    *quantp = quantized;
    *tablep = table;
    return 0;

bailout:
    *quantp=NULL;
    *tablep=NULL;
    if (NULL!=quantized) { vol3d_destroy(quantized); }
    if (NULL!=table) { free(table); }
    return -1;
}

int rgba_write_png_slice(struct vol3d *v, char *filename, unsigned int z) 
{
        FILE *outfile=NULL;
        png_structp png_ptr = NULL;
        png_infop info_ptr = NULL;
        png_byte **row_pointers=NULL;
        const int interlace_flag=1; /* Set for adam7 interlacing */
        int nx,ny;
        int i;
        unsigned char *pixdata;

        if (VOL_TYPE_UINT32!=v->datatype) {
                fprintf(stderr,"rgba_write_png_slice: datatype is not u32\n");
                return -1;
        }


        nx=v->nx; ny=v->ny;

        pixdata = (unsigned char *)v->data + z*nx*ny*4;

        if (NULL==(outfile=fopen(filename,"wb")))
                { perror("fopen()");  return -1; }

        /* Set up the row pointers, because libpng likes to be fed
         * one scanline at a time.
         */

        if (NULL==(row_pointers=malloc(sizeof(png_byte *)*ny)))
                { perror("malloc()");  return -1; }
        for (i=0;i<ny;i++) {
                row_pointers[i] = pixdata + nx*i*4;
        }

        /* Let's do the libpng dance */

        if (NULL==(png_ptr = png_create_write_struct(
                PNG_LIBPNG_VER_STRING,NULL,NULL,NULL))) {
                fprintf(stderr,"png_create_write_struct failed\n");
                free(row_pointers);
                return -1;
        }

        if (NULL==(info_ptr=png_create_info_struct(png_ptr))) {
                png_destroy_write_struct(&png_ptr,NULL);
                free(row_pointers);
                return -1;
        }

        if (setjmp(png_ptr->jmpbuf)) {
                png_destroy_write_struct(&png_ptr,&info_ptr);
                fprintf(stderr,"PNG write error\n");
                free(row_pointers);
                return -1;
        }

        png_init_io(png_ptr,outfile);

        png_set_IHDR(png_ptr,info_ptr,
                        nx,ny,
                        8,      /* bits per channel */
                        PNG_COLOR_TYPE_RGB_ALPHA,
                        interlace_flag?PNG_INTERLACE_ADAM7 : PNG_INTERLACE_NONE,
                        PNG_COMPRESSION_TYPE_DEFAULT,
                        PNG_FILTER_TYPE_DEFAULT
                    );

        /* Write the header */

        png_write_info(png_ptr,info_ptr);
        if (interlace_flag) {
                int niters=png_set_interlace_handling(png_ptr);
                for (i=0;i<niters;i++) {
                        png_write_rows(png_ptr,row_pointers,ny);
                }
        } else {
                png_write_rows(png_ptr,row_pointers,ny);
        }

        png_write_end(png_ptr,info_ptr);
        png_destroy_write_struct(&png_ptr,&info_ptr);

        free(row_pointers);
        fclose(outfile);
        
        return 0;
}

void myCallback(struct lutWidget *self,unsigned char c, void *ptr)
{
    if (' '==c) {
        struct vol3d *quantized=NULL,*vrgba=NULL;
        int i;

        quantized = (struct vol3d *)ptr;

        if (NULL==(vrgba=vol3d_u8_lut_to_rgba(quantized,self->lut))) {
            fprintf(stderr,"vol3d_u8_lut_to_rgba failed.\n");
            return;
        }

        for (i=0;i<vrgba->nz;i++) {
            char buf[16];
            snprintf(buf,16,"slice-%04d.png",i);
            rgba_write_png_slice(vrgba,buf,i);
            printf("Wrote \"%s\"\n",buf);fflush(stdout);
        }

        vol3d_destroy(vrgba);
    }
}


int main(int argc, char *argv[])
{
    char *infilename=NULL;
    struct vol3df_file_hints hints = vol3df_default_file_hints;
    struct vol3d *v=NULL;
    struct u8_freqtable *ft=NULL;
    struct vol3d *quantized = NULL,*vrgba=NULL;
    struct vol3df_stats *stats=NULL;
    struct u8_rgba_lut *lut=NULL;
    struct lutWidget *lw=NULL;

    if ((2!=argc)) {
            fprintf(stderr,"%s infile\n",argv[0]);
            return -1;
    }
    hints.filename = infilename=argv[1]; 

    if (NULL==(v=vol3df_new_from_file_heuristic(infilename,&hints))) {
            fprintf(stderr,"Failed to read input file\n");
            return -1;
    }

    if (NULL==(stats=vol3df_getstats(v))) {
        fprintf(stderr,"vol3df_stats failed\n");
        return -1;
    }

    if (0!=vol3d_f_to_u8_ftable(v,&quantized,&ft,stats->max,stats->min)) {
        fprintf(stderr,"vol3d_f_to_u8_freqtable failed.\n");
        return -1;
    }

    if (NULL==(lut=u8_rgba_lut_default())) {
        fprintf(stderr,"u8_rgba_lut_default failed.\n");
        return -1;
    }

    if (NULL==(vrgba=vol3d_u8_lut_to_rgba(quantized,lut))) {
        fprintf(stderr,"vol3d_u8_lut_to_rgba failed.\n");
        return -1;
    }


    glutInit(&argc,argv);

    if (NULL==(lw=lutWidget_new(ft,lut))) {
        fprintf(stderr,"lutWidget_new failed.\n");
        return -1;
    }

    lutWidget_setKeyCallback(lw,myCallback,quantized);

    glutMainLoop();
    printf("# Done!\n");

    return 0;

}
