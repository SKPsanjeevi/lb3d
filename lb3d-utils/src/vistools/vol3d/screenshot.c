#include "config.h"
#include "screenshot.h"

int cGLerror(char *mess)
{
    GLenum err;
    char *str="(Unknown error)";

    err=glGetError();

    if (GL_NO_ERROR==err) { return 0; }

    switch(err) {
        case GL_INVALID_ENUM: str="GL_INVALID_ENUM"; break;
        case GL_INVALID_VALUE: str="GL_INVALID_VALUE"; break;
        case GL_INVALID_OPERATION: str="GL_INVALID_OPERATION"; break;
        case GL_STACK_OVERFLOW: str="GL_STACK_OVERFLOW"; break;
        case GL_STACK_UNDERFLOW: str="GL_STACK_UNDERFLOW"; break;
        case GL_OUT_OF_MEMORY: str="GL_OUT_OF_MEMORY"; break;
        case GL_TABLE_TOO_LARGE: str="GL_TABLE_TOO_LARGE"; break;
        default: break;
    }

    fprintf(stderr,"%s: GL error: %s\n",mess,str);

    return -1;

}

int gl_screenshot_write_png(unsigned char *pixdata, int nx, int ny, char *filename, int alphaFlag) 
#ifdef HAVE_PNG
{
        FILE *outfile=NULL;
        png_structp png_ptr = NULL;
        png_infop info_ptr = NULL;
        png_byte **row_pointers=NULL;
        const int interlace_flag=1; /* Set for adam7 interlacing */
        int i;
        int channels=3;

        if (alphaFlag) { channels=4; }

        if (NULL==(outfile=fopen(filename,"wb")))
                { perror("fopen()");  return -1; }

        fprintf(stderr,"opened ss file\n");fflush(stderr);

        /* Set up the row pointers, because libpng likes to be fed
         * one scanline at a time.
         */

        if (NULL==(row_pointers=malloc(sizeof(png_byte *)*ny)))
                { perror("malloc()");  return -1; }
        for (i=0;i<ny;i++) {
                row_pointers[i] = pixdata + nx*channels*(ny-i-1);
        }
        fprintf(stderr,"alloced pointers\n");fflush(stderr);

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
                        alphaFlag ? PNG_COLOR_TYPE_RGBA : PNG_COLOR_TYPE_RGB,
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

        fprintf(stderr,"wrote header\n");fflush(stderr);
        png_write_end(png_ptr,info_ptr);
        png_destroy_write_struct(&png_ptr,&info_ptr);

        fprintf(stderr,"wrote png\n");fflush(stderr);

        free(row_pointers);
        fclose(outfile);
        fprintf(stderr,"done png\n");fflush(stderr);
        
        return 0;
}
#else
{
        fprintf(stderr,
           "gl_screenshot_write_png failed: compiled without PNG support.\n");
        return -1;
}
#endif /* HAVE_PNG */



int gl_screenshot_png(int width, int height, char *filename, int alphaFlag)
#ifdef HAVE_PNG
{
    unsigned char *pixbuf=NULL;
    int channels=3;

    if (alphaFlag) { channels=4; }
        fprintf(stderr,"glsp(%dx%d , %s)\n",width,height,filename);fflush(stderr);

    if (NULL==(pixbuf=malloc(channels*width*height)))
    { perror("gl_screenshot_png: malloc"); return -1; }

        fprintf(stderr,"allocated buffer\n");fflush(stderr);

    cGLerror("pre-screenshot");
    fprintf(stderr,"glPixelStorei(GL_PACK_ROW_LENGTH,%d)..",width);
    fflush(stderr);

    glPixelStorei(GL_PACK_ROW_LENGTH,width);
    cGLerror("post");

    fprintf(stderr,"done.\n");
        fflush(stderr);
        fprintf(stderr,"glReadPixels(0,0,%d,%d,%d,GL_UNSIGNED_BYTE,%p)..",
                width,height,alphaFlag ? GL_RGBA : GL_RGB , pixbuf);
        fflush(stderr);
    cGLerror("pre");

    glReadPixels(0,0,width,height,
            alphaFlag ? GL_RGBA : GL_RGB,
            GL_UNSIGNED_BYTE,(GLvoid *)pixbuf);

    fprintf(stderr,"done.\n");
    cGLerror("post");
        fflush(stderr);
        fprintf(stderr,"done glReadPixels\n");fflush(stderr);

    if (0!=gl_screenshot_write_png(pixbuf,width,height,filename,alphaFlag)) {
        fprintf(stderr,"gl_screenshot_png: gl_screenshot_write_png failed.\n");
        free(pixbuf);
        return -1;
    }
        fprintf(stderr,"done gl_screenshot_png\n");fflush(stderr);

    free(pixbuf);

    return 0;
}
#else
{
        fprintf(stderr,
           "gl_screenshot_png failed: compiled without PNG support.\n");
        return -1;
}
#endif /* HAVE_PNG */

