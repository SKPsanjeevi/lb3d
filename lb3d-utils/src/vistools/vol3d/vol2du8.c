#include <stdio.h>
#include <stdlib.h>
#include "config.h"

#ifdef HAVE_PNG
#include <png.h>
#endif /* HAVE_PNG */

#include "vol2d.h"
#include "vol2du8.h"
#include "vol2df.h"

struct vol2d *vol2du8_new(int nx, int ny) {
        struct vol2d *v=NULL;

        if (NULL==(v=vol2d_new(nx,ny,VOL_TYPE_UINT8))) {
                fprintf(stderr,"vol2du8_new: vol2d_new failed\n");
                return NULL;
        }
        return v;
}

#ifdef HAVE_PNG

#define PNGMAGICSIZE 4 /* Size of PNG magic number */





/* Create a new 8-bit unsigned integer vol2d object from
 * a greyscale PNG file
 */
struct vol2d *vol2du8_new_from_gspng(char *fname) 
{
	FILE *f;
	int i,j;
        struct vol2d *v=NULL;

	char pngheader[PNGMAGICSIZE];

	/* PNG metadata */
	png_structp png_ptr;
	png_infop infoptr;
	png_uint_32 pngwidth,pngheight;
	int bpp,coltype,intertype,comptype,filtype;
	int png_transforms = PNG_TRANSFORM_IDENTITY;
	png_bytep *rowptr;
	png_bytep q;
	unsigned char *field,*p;

	/* Open the png file */
	if (NULL==(f=fopen(fname,"rb"))) {
		perror("vol2du8_new_from_gspng: fopen()");
		return NULL;
	}

	/* Check the magic number to make sure it's a PNG */

	if (PNGMAGICSIZE != fread(pngheader, 1, PNGMAGICSIZE,f)) {
		perror("vol2du8_new_from_gspng: fread()");
		return NULL;
	}
	if (0!=png_sig_cmp((png_bytep)pngheader,(png_size_t)0,PNGMAGICSIZE)) {
		fprintf(stderr,
                          "vol2du8_new_from_gspng: Bad PNG magic number.\n");
		return NULL;
	}

	/* Create the PNG structure to hold state for libpng */

	if (NULL==(png_ptr=png_create_read_struct(
		PNG_LIBPNG_VER_STRING,
		NULL,NULL,NULL))) {
			fclose(f);
			fprintf(stderr,"png_create_read_struct() failed\n");
			return NULL;
	}

	/* Allocate structure for image information */

	if (NULL==(infoptr = png_create_info_struct(png_ptr))) {
		fclose(f);
		png_destroy_read_struct(&png_ptr,(png_infopp)NULL,
			(png_infopp)NULL);
		fprintf(stderr,"png_create_info_struct() failed.\n");
		return NULL;
	}

	/* Set up PNG error handling. */

	if (setjmp(png_jmpbuf(png_ptr))) {
		png_destroy_read_struct(&png_ptr,
			(png_infopp)NULL,(png_infopp)NULL);
		fclose(f);
		fprintf(stderr,"Unable to install PNG error handler\n");
		return NULL;
	}

	/* Set up PNG I/O handling */

	png_init_io(png_ptr,f);

	/* Tell it we've already checked the signature */
	png_set_sig_bytes(png_ptr,PNGMAGICSIZE);

	/* Read the goddamned image already */

	png_read_png(png_ptr,infoptr,png_transforms,NULL);

	/* Fetch image header data */

	png_get_IHDR(png_ptr,infoptr, &pngwidth,&pngheight,
		&bpp,&coltype,&intertype,&comptype,&filtype);

        if (NULL==(v=vol2du8_new(pngwidth,pngheight))) {
                fprintf(stderr,"vol2du8_new_form_gspng: vol2du8_new failed\n");
                return NULL;
        }
        field = (unsigned char *)v->data;

	if (8!=bpp) {
		fprintf(stderr,
			"vol2du8_new_from_gspng: PNG file should be 8bpp\n");
		png_destroy_read_struct(&png_ptr,&infoptr,(png_infopp)NULL);
		fclose(f);
                vol2d_destroy(v);
		return NULL;
	}

	/* Get image data */

	rowptr = png_get_rows(png_ptr,infoptr);

	if (PNG_COLOR_TYPE_GRAY != coltype) {
		fprintf(stderr,
		"vol2du8_new_from_gspng: image is not a grayscale PNG\n");
		png_destroy_read_struct(&png_ptr,&infoptr,(png_infopp)NULL);
		fclose(f);
                vol2d_destroy(v);
		return NULL;
	}

	

	/* Transfer image from png to my buffer */

	p=field;
	for (j=0;j<pngheight;j++) {
		q=rowptr[j];
		for (i=0;i<pngwidth;i++) {
			*p++=(unsigned char)(*q++);
		}
	}

	/* Clean up */

	png_destroy_read_struct(&png_ptr,&infoptr,(png_infopp)NULL);
	fclose(f);
	return v;

}

int vol2du8_write_png(struct vol2d *v, char *filename) 
{
        FILE *outfile=NULL;
        png_structp png_ptr = NULL;
        png_infop info_ptr = NULL;
        png_byte **row_pointers=NULL;
        const int interlace_flag=1; /* Set for adam7 interlacing */
        int nx,ny;
        int i;
        unsigned char *pixdata;

        if (VOL_TYPE_UINT8!=v->datatype) {
                fprintf(stderr,"vol2du8_save_png: datatype is not float\n");
                return -1;
        }

        pixdata = (unsigned char *)v->data;

        nx=v->nx; ny=v->ny;

        if (NULL==(outfile=fopen(filename,"wb")))
                { perror("fopen()");  return -1; }

        /* Set up the row pointers, because libpng likes to be fed
         * one scanline at a time.
         */

        if (NULL==(row_pointers=malloc(sizeof(png_byte *)*ny)))
                { perror("malloc()");  return -1; }
        for (i=0;i<ny;i++) {
                row_pointers[i] = pixdata + nx*i;
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
                        PNG_COLOR_TYPE_GRAY,
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

#endif /* HAVE_PNG */

struct vol2d *vol2du8_new_from_2df(struct vol2d *v,float pmax,float pmin)
{
        unsigned char *ucdata=NULL;
        float phi;
        int i;
        int nx,ny;
        struct vol2d *vnew=NULL;
        float *fdata=NULL;

        if (VOL_TYPE_FLOAT != v->datatype) {
                fprintf(stderr,"vol2du8_from_2df: not a float dataset\n");
                return NULL;
        }
        nx=v->nx; ny=v->ny;

        if (NULL==(vnew=vol2du8_new(nx,ny))) {
                fprintf(stderr,"vol2du8_from_2df: vol2du8_new failed\n");
                return NULL;
        }

        ucdata = (unsigned char *)vnew->data;
        fdata = (float *)v->data;

	for (i=0;i<nx*ny;i++) {
		phi = floor(255.0*(fdata[i]-pmin)/(pmax-pmin));
		ucdata[i] = (unsigned char) phi;

	}

        return vnew;

}


struct vol2d *vol2du8_new_from_2df_normalize(struct vol2d *v)
{
        float pmax,pmin,phi;
        int i;
        float *fdata=NULL;

        if (VOL_TYPE_FLOAT != v->datatype) {
                fprintf(stderr,"vol2du8_from_2df: not a float dataset\n");
                return NULL;
        }

        fdata = (float *)v->data;

        pmax=pmin=fdata[0];

        for (i=1;i<v->nx*v->ny;i++) {
                phi = fdata[i];
                if (phi>pmax) { pmax=phi; }
                if (phi<pmin) { pmin=phi; }
        }

        return vol2du8_new_from_2df(v,pmax,pmin);

}
