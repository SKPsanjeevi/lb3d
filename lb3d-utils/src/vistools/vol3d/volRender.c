#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include <sys/types.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <unistd.h>

#include <rpc/rpc.h>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#ifdef HAVE_PNG
#include <png.h>
#endif

#include "vol3d.h"
#include "vol3df.h"
#include "lut.h"
#include "quat.h"
#include "lutWidget.h"
#include "glutViewer.h"
#include "adjFloatList.h"

#include "dirscan.h"

#include "screenshot.h"

#include "volRender.h"

static const int serializeVersion=0;
const char *defaultStateFilename="save.vrState";


#define NCUBEVECTORS 26
vector cubeVector[NCUBEVECTORS];

matrix matrix_fromGLarray(GLfloat *m)
{
    matrix mat;
    int i;
    for (i=0;i<9;i++) {
        mat.el[i] = m[i];
    }

    return mat;
}

void getNearestCubeVecs(vector *forwardp, vector *upp, GLfloat *m)
{

    matrix mat;
    vector upVec, forwardVec,cubeForward,cubeUp;
    GLfloat dp,max;
    int i;

    mat = matrix_fromGLarray(m);
    upVec = vector_normalized(vector_matrix_dot(mat,vector_new(0,1,0)));
    forwardVec = vector_normalized(vector_matrix_dot(mat,vector_new(0,0,-1)));

    /* Find the cubeVector with the highest dot product with the forwardVec. */
    /* This is the new forward vector, cubeForward. */

    max=0.0;

    for (i=0;i<NCUBEVECTORS;i++) {
        dp = vector_dot(cubeVector[i],forwardVec);
        if (dp>max) {
            max=dp;
            cubeForward = cubeVector[i];
        }
    }


    /* Now find the cubeVector with the highest DP with upVec, but which
     * is also perpendicular to the cubeForward vector. This is the new upVec.
     */

    max=0.0;

    for (i=0;i<NCUBEVECTORS;i++) {
        /* Skip vectors which are not normal to cubeForward. */
        if (fabs(vector_dot(cubeVector[i],cubeForward))>1e-6) { continue; }

        dp = vector_dot(cubeVector[i],upVec);
        if (dp>max) {
            max=dp;
            cubeUp = cubeVector[i];
        }
    }


    *forwardp = cubeForward;
    *upp = cubeUp;

    return;
}

void glVertex_vec(vector v)
{
    glVertex3f(v.x,v.y,v.z);
}

void snapToCube(struct glutViewer *gv)
{
    GLfloat modelview[16];
    GLfloat m[9];
    GLfloat mNew[9];
    int i,j;
    vector cubeForward,cubeUp,cubeRight;
    quat_t newQuat;

    glGetFloatv(GL_MODELVIEW_MATRIX,modelview);

    /* Extract rotation component */

    for (i=0;i<3;i++) {
    for (j=0;j<3;j++) {
        m[3*i+j] = modelview[4*i+j];
    }}

    getNearestCubeVecs(&cubeForward,&cubeUp,m);

    cubeRight = vector_cross(cubeForward,cubeUp);

    cubeForward = vector_normalized(cubeForward);
    cubeUp = vector_normalized(cubeUp);
    cubeRight = vector_normalized(cubeRight);

    mNew[0] = cubeRight.x;
    mNew[3] = cubeRight.y;
    mNew[6] = cubeRight.z;
    mNew[1] = cubeUp.x;
    mNew[4] = cubeUp.y;
    mNew[7] = cubeUp.z;
    mNew[2] =-cubeForward.x;
    mNew[5] =-cubeForward.y;
    mNew[8] =-cubeForward.z;

    newQuat = quatfromrotmatrix(mNew);

    gv->orquat = newQuat;

    glutViewer_postRedisplay(gv);

}


void initCubeVectors(void)
{
        cubeVector[ 0] = vector_new(1,0,0);
        cubeVector[ 1] = vector_new(0,1,0);
        cubeVector[ 2] = vector_new(0,0,1);
        cubeVector[ 3] = vector_new(-1,0,0);
        cubeVector[ 4] = vector_new(0,-1,0);
        cubeVector[ 5] = vector_new(0,0,-1);

        cubeVector[ 6] = vector_normalized(vector_new(-1,-1,-1));
        cubeVector[ 7] = vector_normalized(vector_new(-1,-1, 1));
        cubeVector[ 8] = vector_normalized(vector_new(-1, 1,-1));
        cubeVector[ 9] = vector_normalized(vector_new(-1, 1, 1));
        cubeVector[10] = vector_normalized(vector_new( 1,-1,-1));
        cubeVector[11] = vector_normalized(vector_new( 1,-1, 1));
        cubeVector[12] = vector_normalized(vector_new( 1, 1,-1));
        cubeVector[13] = vector_normalized(vector_new( 1, 1, 1));

        cubeVector[14] = vector_normalized(vector_new(0,-1,-1));
        cubeVector[15] = vector_normalized(vector_new(0,-1, 1));
        cubeVector[16] = vector_normalized(vector_new(0, 1,-1));
        cubeVector[17] = vector_normalized(vector_new(0, 1, 1));
        cubeVector[18] = vector_normalized(vector_new(-1,0,-1));
        cubeVector[19] = vector_normalized(vector_new(-1,0, 1));
        cubeVector[20] = vector_normalized(vector_new( 1,0,-1));
        cubeVector[21] = vector_normalized(vector_new( 1,0, 1));
        cubeVector[22] = vector_normalized(vector_new(-1,-1,0));
        cubeVector[23] = vector_normalized(vector_new(-1, 1,0));
        cubeVector[24] = vector_normalized(vector_new( 1,-1,0));
        cubeVector[25] = vector_normalized(vector_new( 1, 1,0));

}

int volRender_xdr(struct volRender *self, XDR *xdrp)
{
    int i,j;
    struct glutViewer *gv=NULL;

    if (NULL==(gv=self->gv)) {
        fprintf(stderr,"volRender_xdr: NULL gv!\n");
        fflush(stderr);
        return -1;
    }

    /* Serialize selected glutViewer state */

    xdr_int(xdrp,&gv->orthographic);
    xdr_float(xdrp,&gv->orthoSize);
    xdr_float(xdrp,&gv->fovy);
    xdr_float(xdrp,&gv->nearClip);
    xdr_float(xdrp,&gv->geomx);
    xdr_float(xdrp,&gv->geomy);
    xdr_float(xdrp,&gv->geomz);

    xdr_float(xdrp,&gv->orquat.q0);
    xdr_float(xdrp,&gv->orquat.q1);
    xdr_float(xdrp,&gv->orquat.q2);
    xdr_float(xdrp,&gv->orquat.q3);

    xdr_float(xdrp,&gv->sphererad);

    /* Serialize selected volRender state */

    /* Serialize RGBA LUT */

    for (i=0;i<256;i++) {
        for (j=0;j<4;j++) {
             xdr_u_char(xdrp,&(self->lut->ent[i].c[j]));
        }
    }
    xdr_int(xdrp,&self->slicesLow);
    xdr_int(xdrp,&self->slicesHigh);
    xdr_float(xdrp,&self->fogStart);
    xdr_float(xdrp,&self->fogEnd);
    xdr_float(xdrp,&self->fogDensity);
    xdr_float(xdrp,&self->globalAlpha);
    xdr_int(xdrp,&self->fogOn);
    xdr_int(xdrp,&self->wrapFrames);

    return 0;
}


int volRender_serializeFH(struct volRender *self,FILE *f,
        enum volRenderFileMode mode)
{
    XDR xdrs;

    /* Create the appropriate XDR structure, and write or
     * verify the state file version number.
     */

    if (VOLRENDER_READ==mode) {
       int version;

       xdrstdio_create(&xdrs,f,XDR_DECODE);
       xdr_int(&xdrs,&version);

       if (version!=serializeVersion) {
           fprintf(stderr,"volRender_serializeFH: version mismatch:\n"
                   "File version %d, program version %d\n",
                   version,serializeVersion);
           xdr_destroy(&xdrs);
           return -1;
       }

    } else if (VOLRENDER_WRITE==mode) {

        xdrstdio_create(&xdrs,f,XDR_ENCODE);
        xdr_int(&xdrs,(int *)&serializeVersion);


    } else {
        fprintf(stderr,"volRender_serializeFH: bad mode %d\n",mode);
        xdr_destroy(&xdrs);
        return -1;
    }

    /* Call XDR (de)serialization routine to read/write actual state */

    if (0!=volRender_xdr(self,&xdrs)) {
        fprintf(stderr,
                "volRender_serializeToFH: volRender_xdr failed.\n");
        xdr_destroy(&xdrs);
        return -1;
    }

    xdr_destroy(&xdrs);

    return 0;
}

int volRender_serializeFile(struct volRender *self, char *fname,
        enum volRenderFileMode mode)
{
    FILE *fh=NULL;
    char *fileMode;

    if (VOLRENDER_READ==mode) {
        fileMode = "rb";
    } else if (VOLRENDER_WRITE==mode) {
        fileMode = "wb";
    } else {
        fprintf(stderr,"volRender_serializeFile: bad mode \"%d\"\n",mode);
        return -1;
    }

    if (NULL==(fh=fopen(fname,fileMode))) {
        perror("volRender_serializeFile: fopen");
        return -1;
    }

    if (0!=volRender_serializeFH(self,fh,mode)) {
        fprintf(stderr,
                "volRender_serializeFile: volRender_serializeFH failed.\n");
        fclose(fh);
        return -1;
    }

    fclose(fh);
    return 0;
}

long gettime_ms(void)
{
    struct timeval tv;

    gettimeofday(&tv,NULL);

    return 1000*tv.tv_sec + (tv.tv_usec/1000);

}


int checkGLerror(char *mess)
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

    for (n=0;n<256;n++) { table->freq[n]=0; }

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
#ifdef DEBUG
    {
        int i;
        fprintf(stderr,"Rebuilt frequency table\n\n");fflush(stderr);
        for (i=0;i<256;i++) {
            fprintf(stderr,"%d %d\n",i,table->freq[i]);
        }
        fprintf(stderr,"\n");fflush(stderr);
    }
#endif /* DEBUG */
    return 0;

bailout:
    *quantp=NULL;
    *tablep=NULL;
    if (NULL!=quantized) { vol3d_destroy(quantized); }
    if (NULL!=table) { free(table); }
    return -1;
}

int rgba_write_png_slice(struct vol3d *v, char *filename, unsigned int z) 
#ifdef HAVE_PNG
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
#else 
{
    fprintf(stderr,
            "rgba_write_png_slice failed: compiled without PNG support\n");
    return -1;
}
#endif /* HAVE_PNG */

void vertex(GLfloat tx, GLfloat ty, GLfloat tz,
        GLfloat x, GLfloat y, GLfloat z)
{
    glTexCoord3f(tx,ty,tz);
    glVertex3f(x,y,z);

    /*
    fprintf(stderr,"tex (%f,%f,%f) v (%f,%f,%f)\n",
            tx,ty,tz,
            x,y,z);
    fflush(stderr);
    */
}

void myDisplay(struct glutViewer *gv, void *data)
{
    struct volRender *vr=NULL;
    GLfloat cubeSize=1.0; /* Size of volume cube */
    GLfloat spacing; /* Separation between polygons */
    GLfloat l;
    int i;
    GLfloat m[16],mi[16];
    quat_t orientation;
    GLfloat root3,iroot3;

    long hiRenderingTime=0; /* Time to render at high res */

    /*
    GLfloat ambient[]={1.0,1.0,1.0,1.0};
    GLfloat shininess[]={50};
    GLfloat position[]={1,1,1,1};
    GLfloat specular[]={1.0,1.0,1.0,1.0};
    */


    vr = (struct volRender *)data;

    hiRenderingTime = gettime_ms();

    root3=sqrt(3.0);
    iroot3=1.0/root3;

    cubeSize=1.0;
    l=0.5;
    spacing = cubeSize/(GLfloat)vr->slices;

    orientation = glutViewer_getQuat(gv);
    quatmatrixandinv(gv->orquat,m,mi);


    glDisable(GL_FOG);
    glDisable(GL_TEXTURE_3D);
    glEnable(GL_COLOR_MATERIAL);
    glDisable(GL_BLEND);
    glColor4f(1,1,1,1);
    if (vr->drawBox) {
        glutWireCube(1.0/sqrt(3.0));
    }

    if (0) { /* Render glyphs */


        glPushMatrix();

            glScalef(iroot3,iroot3,iroot3);
            glTranslatef(-0.5,-0.5,-0.5);

            glPointSize(10.0);
            glColor4f(1,0,1,0.5);

            glEnable(GL_POINT_SMOOTH);
            glDisable(GL_BLEND);

            glBegin(GL_POINTS);

                glColor4f(1,1,1,1);
                glVertex3f(0,0,0);

                glColor4f(1,0,1,1);
                glVertex3f(0.2,0.2,0.8);
                glVertex3f(0.2,0.8,0.2);
                glVertex3f(0.2,0.8,0.8);
                glVertex3f(0.8,0.2,0.2);
                glVertex3f(0.8,0.2,0.8);
                glVertex3f(0.8,0.8,0.2);
                glVertex3f(0.8,0.8,0.8);

            glEnd();

        glPopMatrix();

    }

    if (vr->fogOn) {
        glFogi(GL_FOG_MODE,GL_EXP2);
        glFogfv(GL_FOG_COLOR,gv->bgColor);
        glFogf(GL_FOG_DENSITY,vr->fogDensity);
        glFogf(GL_FOG_START,vr->fogStart);
        glFogf(GL_FOG_END,vr->fogEnd);
        glHint(GL_FOG,GL_NICEST);

        glEnable(GL_FOG);

    } else { glDisable(GL_FOG); }

    /* Move into billboard mode */
    glDisable(GL_COLOR_MATERIAL);
    glEnable(GL_TEXTURE_3D);

    glMatrixMode(GL_MODELVIEW); glPushMatrix();

    glGetFloatv(GL_MODELVIEW_MATRIX,m);

    m[0]=m[5]=m[10]=1;
    m[1]=m[2]=
    m[4]=m[6]=
    m[8]=m[9]=0;

    glLoadMatrixf(m);



    /* Now set the texture matrix to orient the volume */

    glMatrixMode(GL_TEXTURE); glPushMatrix();

    glLoadIdentity();
    glTranslatef(0.5,0.5,0.5);
    glScalef(root3,root3,root3);
    glMultMatrixf(mi);

    glEnable(GL_TEXTURE_3D);
    glDisable(GL_CULL_FACE);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    /*glBlendFunc(GL_CONSTANT_ALPHA,GL_ONE); */
    glColor4f(1,1,1,vr->globalAlpha);

    glBegin(GL_QUADS);
    for (i=0;i<vr->slices;i++) {
        GLfloat z,tz;

        z=i*spacing-l;
        tz=i*spacing-0.5;

        vertex( -0.5,  0.5, tz,-l, l, z);
        vertex(  0.5,  0.5, tz, l, l, z);
        vertex(  0.5, -0.5, tz, l,-l, z);
        vertex( -0.5, -0.5, tz,-l,-l, z);

    }
    glEnd();

    glPopMatrix(); /* GL_TEXTURE */

    glMatrixMode(GL_MODELVIEW); glPopMatrix();

}

/* Let's not bother with viewer-aligned slices.
 * Just set the border regions to alpha=0.
 * XXX This is a truly gruesome hack.
 */
void vrgba_cheat(struct vol3d *quantized, struct vol3d *vrgba)
{
    unsigned long x,y,z;

    /* +/- X planes */

    for (z=0;z<vrgba->nz;z++) {
    for (y=0;y<vrgba->ny;y++) {
        *((uint32_t *)vol3d_element(vrgba,0,y,z))=0;
        *((uint32_t *)vol3d_element(vrgba,vrgba->nx-1,y,z))=0;
    }}

    /* +/- Y planes */

    for (z=0;z<vrgba->nz;z++) {
    for (x=0;x<vrgba->nx;x++) {
        *((uint32_t *)vol3d_element(vrgba,x,0,z))=0;
        *((uint32_t *)vol3d_element(vrgba,x,vrgba->ny-1,z))=0;
    }}

    /* +/- Z planes */

    for (x=0;x<vrgba->nx;x++) {
    for (y=0;y<vrgba->ny;y++) {
        *((uint32_t *)vol3d_element(vrgba,x,y,0))=0;
        *((uint32_t *)vol3d_element(vrgba,x,y,vrgba->nz-1))=0;
    }}
}

int volRender_setTexture(struct volRender *self)
{
    int oldWindow;

    oldWindow = glutGetWindow();

    glutSetWindow(self->gv->win);

    vrgba_cheat(self->quantized[self->currVol],self->vrgba);

    checkGLerror("volRender_setTexture: pre-existing error");

    glTexImage3D(
        GL_TEXTURE_3D, /* target */
        0, /* level */
        GL_RGBA, /* internalFormat (no of channels) */
        self->quantized[self->currVol]->nx, /* width/height/depth */
        self->quantized[self->currVol]->ny,
        self->quantized[self->currVol]->nz,
        0, /* border */
        GL_RGBA, /* format */
        GL_UNSIGNED_BYTE, /* type */
        self->vrgba->data);
    checkGLerror("glTexImage3D");

    glutSetWindow(oldWindow);
    return 0;
}

/* Map the current volume through the LUT and redisplay */

int volRender_reload(struct volRender *self)
{
    struct vol3d *vrgba=NULL;
    if (NULL==(vrgba=vol3d_u8_lut_to_rgba(
                    self->quantized[self->currVol],
                    self->lw->lut))) {

        fprintf(stderr,"volRender_reload: vol3d_u8_lut_to_rgba failed.\n");
        return -1;
    }

    if (NULL!=self->vrgba) { vol3d_destroy(self->vrgba); }

    self->vrgba = vrgba;
    volRender_setTexture(self);
    glutViewer_postRedisplay(self->gv);

    return 0;
}

void myKeyCallback(struct lutWidget *self,unsigned char c, void *ptr)
{
    struct volRender *vr=NULL;
    if (' '==c) {

        fprintf(stderr,"Recalculating RGBA volume...");fflush(stderr);

        if (NULL==(vr=(struct volRender *)ptr)) {
            fprintf(stderr,"myKeyCallback: NULL data!\n");
            fflush(stderr);
            return;
        }
        volRender_reload(vr);
        fprintf(stderr,"done.\n");fflush(stderr);
    }
}

int volRender_destroy(struct volRender *self)
{
    if (NULL==self) { return -1; }
    if (NULL!=self->quantized) { free(self->quantized); }
    if (NULL!=self->ft) { free(self->ft); }
    if (NULL!=self->afl) { adjFloatList_destroy(self->afl); }
    /* Don't free the names array; let the user deal with it. */
    free(self);
    return 0;
}

struct volRender *volRender_new(int nvols)
{
    struct volRender *self=NULL;
    int i;

    if (NULL==(self=malloc(sizeof(struct volRender))))
        { perror("volRender_new: malloc"); return NULL; }

    if (NULL==(self->quantized=malloc(sizeof(struct vol3d *)*nvols))) {
        perror("volRender_new: malloc");
        volRender_destroy(self);
        return NULL;
    }

    if (NULL==(self->ft=malloc(sizeof(struct u8_freqtable *)*nvols))) {
        perror("volRender_new: malloc");
        volRender_destroy(self);
        return NULL;
    }

    if (NULL==(self->afl=adjFloatList_new())) {
        fprintf(stderr,"volRender_new: adjFloatList_new failed.\n");
        volRender_destroy(self);
        return NULL;
    }

    adjFloatList_wrapOn(self->afl);


    self->name=NULL;
    self->nVolumes = nvols;
    for (i=0;i<nvols;i++) { self->quantized[i]=NULL; self->ft[i]=NULL; }


    self->gv=NULL;
    self->vrgba = NULL;
    self->lut = NULL;
    self->lw = NULL;
    self->currVol=0;
    self->drawBox=1;

    self->fogStart = 1.0;
    self->fogEnd = 64.0;
    self->fogDensity = 0.0465;

    self->slicesLow = 256;
    self->slicesHigh = 512;
    self->slices = 512;

    self->globalAlpha = 0.1;

    self->lastClickTime=0;
    self->lodTime = 500;

    adjFloatList_addFloat(self->afl,
            &self->globalAlpha, 0.0, 1.0, 0.005, "globalAlpha");
    adjFloatList_addFloat(self->afl,
            &self->fogStart, 0.0, 100.0, 1.0, "fogStart");
    adjFloatList_addFloat(self->afl,
            &self->fogEnd, 0.0, 100.0, 1.0, "fogEnd");
    adjFloatList_addFloat(self->afl,
            &self->fogDensity, 0.0, 10.0, 0.001, "fogDensity");

    self->newestTime=0;
    self->minAge = 3; /* Default number of seconds a file has to sit still */
    self->prefix=NULL;
    self->rangeMax=1.0;
    self->rangeMin=0.0;
    self->useRange=0;

    return self;
}

void volRender_keyboardFunc(unsigned char key, int x, int y)
{
    struct volRender *self=NULL;
    struct glutViewer *gv=NULL;

    if (NULL==(gv=winReg_getPtr(0))) {
        fprintf(stderr,"volRender_keyboardFunc: winReg_getPtr failed\n");
        return;
    }

    if (NULL==(self = gv->data)) {
        fprintf(stderr,"volRender_keyboardFunc: NULL gv->data!\n");
        return;
    }

    if ('Q'==key) {
        fprintf(stderr,"Closing.\n");fflush(stderr);
        exit(0);
    }

    switch(key) {
        case 's': /* Screenshot */
        fprintf(stderr,"caught screenshot\n");fflush(stderr);
            {
                char *pngname=NULL;
                unsigned int pngnamelen=0;

                if (NULL!=self->name) {
                    /* pngname is the filename plus ".png" (plus nul) */
                    pngnamelen = 5+strlen(self->name[self->currVol]) ;
                    if (NULL!=(pngname=malloc(pngnamelen))) {
                        strncpy(pngname,self->name[self->currVol],pngnamelen);
                        strncat(pngname,".png",5);
                        pngname[pngnamelen-1]=0; /* Always null-terminate */
                    }
                }

                if (NULL==pngname) {
                    gl_screenshot_png(gv->width,gv->height,"screenshot.png",
                            gv->alphaBufferFlag);
                    fprintf(stderr,"Saved screenshot.png\n"); fflush(stderr);
                } else {
                    gl_screenshot_png(gv->width,gv->height,pngname,
                            gv->alphaBufferFlag);
                    fprintf(stderr,"Saved %s\n",pngname); fflush(stderr);
                    free(pngname);
                    pngname=NULL; pngnamelen=0;
                }
            }
            break;

        case 'S': /* Save state to disk */
        fprintf(stderr,"Saving viewer state...\n");fflush(stderr);
            if(0==volRender_serializeFile(self,
                        (char *)defaultStateFilename,VOLRENDER_WRITE)){
                fprintf(stderr,"Saved to %s\n",defaultStateFilename);
                fflush(stderr);
            } else {
                fprintf(stderr,"Failed to save state!\n"); fflush(stderr);
            }
            break;

        case 'p': /* Toggle orthographic projection */
            glutViewer_toggleOrtho(gv);
            fprintf(stderr,"Orthographic projection %s.\n",
                    gv->orthographic ? "on" : "off");
            fflush(stderr);
            break;

        case 'o': /* Toggle outline box */
            self->drawBox = !self->drawBox;
            fprintf(stderr,"Outline box %s.\n",
                    self->drawBox ? "on" : "off");
            fflush(stderr);
            break;
            
        case 'f': /* Toggle fogging */ 
            if (self->fogOn) { self->fogOn=0; } else { self->fogOn=1; }
            fprintf(stderr,"Fogging turned %s.\n",
                    self->fogOn ? "on" : "off");
            fflush(stderr);
            break;
        case '>': /* Next float value */
            adjFloatList_next(self->afl);
            fprintf(stderr,"%s = %f\n",
                    adjFloatList_getCurrentName(self->afl),
                    adjFloatList_getCurrentValue(self->afl)
                   );
            fflush(stderr);
            break;
        case '<': /* Prev float value */
            adjFloatList_prev(self->afl);
            fprintf(stderr,"%s = %f\n",
                    adjFloatList_getCurrentName(self->afl),
                    adjFloatList_getCurrentValue(self->afl)
                   );
            fflush(stderr);
            break;
        case '.': /* Increment float */
            adjFloatList_incFloat(self->afl);
            fprintf(stderr,"%s = %f\n",
                    adjFloatList_getCurrentName(self->afl),
                    adjFloatList_getCurrentValue(self->afl)
                   );
            break;
        case ',': /* Decrement float */
            adjFloatList_decFloat(self->afl);
            fprintf(stderr,"%s = %f\n",
                    adjFloatList_getCurrentName(self->afl),
                    adjFloatList_getCurrentValue(self->afl)
                   );
            break;
        case '`':
            snapToCube(gv);
            break;
        case '+':
            wheelUpFunc(gv,NULL);
            break;
        case '-':
            wheelDownFunc(gv,NULL);
            break;

        default:
            break;

    }
    glutViewer_postRedisplay(gv);
}

void volRender_advanceFrame(struct volRender *self, int delta)
{
    int newFrame;
    int oldWin;

    newFrame = self->currVol + delta;
    if (newFrame>=self->nVolumes) {
        if (self->wrapFrames) {
            newFrame = 0;
        } else {
            newFrame = self->currVol;
        }
    }

    if (newFrame<0) {
        if (self->wrapFrames) {
            newFrame = self->nVolumes-1;
        } else {
            newFrame = 0;
        }
    }
    self->currVol=newFrame;




    volRender_reload(self);
    lutWidget_setFT(self->lw,self->ft[self->currVol]);

    if (NULL!=self->name) {
        fprintf(stderr,"Displaying \"%s\"\n",self->name[self->currVol]);
        oldWin=glutGetWindow();
            glutSetWindowTitle(self->name[self->currVol]);
        glutSetWindow(oldWin);
    } else {
        fprintf(stderr,"Displaying volume %d\n",self->currVol);
    }

    fflush(stderr);

}


void wheelUpFunc(struct glutViewer *gv, void *data)
{
    struct volRender *self=NULL;

    if (NULL==(self=gv->data)) {
        fprintf(stderr,"wheelUpFunc: NULL gv->data !\n");
        fflush(stderr);
    }
    volRender_advanceFrame(self,+1);
}

void wheelDownFunc(struct glutViewer *gv, void *data)
{
    struct volRender *self=NULL;

    if (NULL==(self=gv->data)) {
        fprintf(stderr,"wheelUpFunc: NULL gv->data !\n");
        fflush(stderr);
    }
    volRender_advanceFrame(self,-1);
}

    

/* Quantize the given volume into the ith volume
 * in the volume renderer.
 */

int volRender_quantizeVolume(struct volRender *self,struct vol3d *vol, int i)
{
    struct vol3df_stats *stats=NULL;
    if (self->useRange) {
        if (0!=vol3d_f_to_u8_ftable(vol,
                    &(self->quantized[i]),
                    &(self->ft[i]),
                    self->rangeMax,
                    self->rangeMin)) {
            fprintf(stderr,"vol3d_f_to_u8_freqtable failed.\n");
            return -1;
        }
    } else {
        if (NULL==(stats=vol3df_getstats(vol))) {
            fprintf(stderr,"vol3df_stats failed\n");
            return -1;
        }

        if (0!=vol3d_f_to_u8_ftable(vol,
                    &(self->quantized[i]),
                    &(self->ft[i]),
                    stats->max,
                    stats->min)) {
            fprintf(stderr,"vol3d_f_to_u8_freqtable failed.\n");
            return -1;
        }
    }

    return 0;
}

/* Load and quantize the named file into the ith volume
 * in the volume renderer.
 */

int volRender_loadVolume(struct volRender *self, char *filename, int i)
{
    struct vol3df_file_hints hints = vol3df_default_file_hints;
    struct vol3d *vol=NULL;

    hints.filename = filename; 

    if (NULL==(vol=vol3df_new_from_file_heuristic(filename,&hints))) {
            fprintf(stderr,"Failed to read input file\n");
            return -1;
    }

    if (0!=volRender_quantizeVolume(self,vol,i)) {
        fprintf(stderr,
              "volRender_loadVolume: volRender_quantizeVolume failed.\n");
        vol3d_destroy(vol);
        return -1;
    }

    vol3d_destroy(vol);

    fprintf(stderr,"Loaded and quantized \"%s\"\n",filename);
    fflush(stderr);

    return 0;
}


/* Scanner timer callback.
 * See if there are any files which: (a) start with self->prefix,
 * (b) have been modified since self->newestTime, and (c) are older
 * than self->minAge. If so, load the newest one.
 */
void volRender_scannerTimerFunc(int win)
{
    struct glutViewer *gv=NULL;
    struct volRender *self=NULL;
    int winOld;

    char *newestFilename=NULL;

    struct vol3df_file_hints fileHints = vol3df_default_file_hints;
    struct vol3d *v=NULL;
    struct stat statbuf;

    /* Find out volrender/window context */

    winOld = glutGetWindow();
    glutSetWindow(win);

    if (NULL==(gv=winReg_getPtr(0))) {
        fprintf(stderr,"volRender_scannerTimerFunc: winReg_getPtr failed\n");
        glutSetWindow(winOld);
        return;
    }
    self=(struct volRender *)gv->data;

    fprintf(stderr,"Scanning..\n");fflush(stderr);
    glutTimerFunc(1000*self->minAge,volRender_scannerTimerFunc,win);

    if (NULL==(
     newestFilename=findRipeFile(self->prefix,self->newestTime,self->minAge)))
        { return; } /* Nothing found. */ 

    /* Update mtime */

    if (0!=stat(newestFilename,&statbuf)) { perror("stat"); return ; }
    self->newestTime = statbuf.st_mtime;
    fprintf(stderr,"Newest file is %s\n",newestFilename);

    /* Attempt to load it */       

    fileHints.filename = newestFilename;
    if (NULL==(v=vol3df_new_from_file_heuristic(
        fileHints.filename,&fileHints))) {
        fprintf(stderr,"Failed to load \"%s\".\n",newestFilename);
        return ; /* oh well, doesn't matter. */
    }
    fprintf(stderr,"Loaded %dx%dx%d voxels from %s\n",
            v->nx,v->ny,v->nz,fileHints.filename);
    fflush(stderr);

    if (0!=volRender_quantizeVolume(self,v,0)) {
        fprintf(stderr,
              "volRender_loadVolume: volRender_quantizeVolume failed.\n");
        vol3d_destroy(v);
        return ;
    }

    vol3d_destroy(v);

    fprintf(stderr,"Loaded and quantized \"%s\"\n",newestFilename);
    fflush(stderr);

    /* Free old "current" filename ; install new one. */
    if (NULL!=self->name) {
        free(self->name[0]);
        self->name[0] = newestFilename;
    }

    volRender_reload(self);
    lutWidget_setFT(self->lw,self->ft[self->currVol]);
}

/* LOD timer callback.
 * Check if more than lodTime milliseconds have elapsed since
 * lastClickTime. If so, set the LOD to high, and post a redisplay.
 */

void volRender_lodTimerFunc(int win)
{
    long now;
    struct glutViewer *gv=NULL;
    struct volRender *self=NULL;
    long delta=5;
    long elapsed_ms;
    int winOld;

    winOld = glutGetWindow();
    glutSetWindow(win);

    if (NULL==(gv=winReg_getPtr(0))) {
        fprintf(stderr,"volRender_lodTimerFunc: winReg_getPtr failed\n");
        glutSetWindow(winOld);
        return;
    }
    self=(struct volRender *)gv->data;

    now=gettime_ms();

    elapsed_ms = now - self->lastClickTime;

    if ( labs( elapsed_ms - self->lodTime) <= delta ) {
        self->slices = self->slicesHigh;
        glutViewer_postRedisplay(gv);
    } else {
       glutTimerFunc(labs(elapsed_ms-self->lodTime),volRender_lodTimerFunc,win);
    }

    glutSetWindow(winOld);
}

/* Called whenever the mouse is dragged.
 * Set the lastClickTime to the current time, set the LOD to low,
 * and set a timer to be called in lodTime milliseconds.
 * Then drop through to the glutViewer callback.
 */
void volRender_motionFunc(int x, int y)
{
    struct glutViewer *gv=NULL;
    struct volRender *self=NULL;

    if (NULL==(gv=winReg_getPtr(0))) {
        fprintf(stderr,"volRender_motionFunc: winReg_getPtr failed\n");
        return;
    }

    self=(struct volRender *)gv->data;
    self->lastClickTime = gettime_ms();

    if (self->slicesLow == self->slices) {
        glutViewer_motionFunc(x,y);
        return;
    }

    self->slices = self->slicesLow;


    glutTimerFunc(self->lodTime,volRender_lodTimerFunc,gv->win);

    glutViewer_motionFunc(x,y);
}


