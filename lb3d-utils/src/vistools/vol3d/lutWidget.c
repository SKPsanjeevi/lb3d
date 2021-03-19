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

#include "lutWidget.h"





/* Convert HSV colour specification to RGB.
 * Adapted from Foley et al, "Introduction to Computer Graphics",
 * 1990 Addison-Wesley, p415-416
 */
static void hsv2rgb(GLfloat *r,GLfloat *g,GLfloat *b,
	GLfloat h,GLfloat s,GLfloat v)
{
	GLfloat f,p,q,t;
	int i;

	if (0.0==s) {	/* Colour in black-and-white centre line */
		*r=*g=*b=v;
		return;
	}
	while (h>360.0) {
		h-=360.0;
	}
	/* h is now below 360 */
	h/=60.0;

	i=floor(h);	/* h now in [0,6] */
	f=h-i; /* f is fractional part of h */
	p=v*(1-s);
	q=v*(1-s*f);
	t=v*(1-s*(1-f));

	switch(i) {
		case 0:
			*r=v; *g=t; *b=p; break;
		case 1:
			*r=q; *g=v; *b=p; break;
		case 2:
			*r=p; *g=v; *b=t; break;
		case 3:
			*r=p; *g=q; *b=v; break;
		case 4:
			*r=t; *g=p; *b=v; break;
		case 5:
			*r=v; *g=p; *b=q; break;
	}
	
}

/* Return a new LUT with default values */

struct u8_rgba_lut *u8_rgba_lut_default(void)
{
        struct u8_rgba_lut *lut=NULL;
        uint8_t i; 

        if (NULL==(lut=malloc(sizeof(struct u8_rgba_lut)))) {
            perror("u8_rgba_lut_default: malloc");
            return NULL;
        }

        for (i=0;i<255;i++) {
            GLfloat h,s,v;
            GLfloat r,g,b;
            h = ((float)i)/255.0*360.0;
            s = 1.0;
            v = 1.0;
            hsv2rgb(&r,&g,&b,h,s,v);

            lut->ent[i].c[0] = 255*r;
            lut->ent[i].c[1] = 255*g;
            lut->ent[i].c[2] = 255*b;
            lut->ent[i].c[3] = i;
        }

        return lut;
}

/* Return true if x is a positive integer power of 2; false otherwise. */

int isPowerOfTwo(int x) {
        if (x<=0) { return 0; }
        while (0==(x & 1)) { x >>= 1; }
        return (1==x);
}

/* Take a positive integer; if it's a power of two, return it; else 
 * return the next highest power of two. */

static int nextPowerOfTwo(int x) {
    int xOriginal = x;
    int y=2;

    if (x<=0) { return 0; } /* Not a positive integer! */

    while (1!=x) { x>>=1; y<<=1; }

    if (xOriginal == (y>>1)) { return xOriginal; }

    return y;
}


struct vol3d *vol3d_u8_lut_to_rgba(struct vol3d *u8, 
        struct u8_rgba_lut *lut)
{
    struct vol3d *vrgba=NULL;
    uint8_t *rgba=NULL;
    uint8_t *p=NULL;
    int rgbx,rgby,rgbz;
    int ux,uy,uz;
    int i;

    /* Make sure the RGBA volume dimensions are all integer powers of 2. */
    rgbx = nextPowerOfTwo(u8->nx);
    rgby = nextPowerOfTwo(u8->ny);
    rgbz = nextPowerOfTwo(u8->nz);

    if (VOL_TYPE_UINT8 != u8->datatype) {
        fprintf(stderr,"vol3d_u8_lut_to_rgba: need a UINT8 volume\n");
        return NULL;
    }

    if (NULL==(vrgba=vol3d_new(rgbx,rgby,rgbz,VOL_TYPE_UINT32))) {
        fprintf(stderr,"vol3d_u8_lut_to_rgba: vol3d_new failed.\n");
        return NULL;
    }


    p=(uint8_t *)u8->data;

    for (uz=0;uz<u8->nz;uz++) {
        for (uy=0;uy<u8->ny;uy++) {
           for (ux=0;ux<u8->nx;ux++) {

               /* Read a byte from the volume */
               uint8_t pixval=*p++;

               /* Find address of corresponding texel */

               rgba = (uint8_t *)vrgba->data 
                    + 4*( ux + rgbx*uy + rgbx*rgby*uz);

               /* Look up the byte , write an RGBA tuple to the texture. */
               for (i=0;i<4;i++) { *rgba++ = lut->ent[pixval].c[i]; }
           }
        }
    }
    fprintf(stderr,"rgb=(%d,%d,%d) b=(%d,%d,%d)\n",
            rgbx,rgby,rgbz,u8->nx,u8->ny,u8->nz);fflush(stderr);
    
    return vrgba;

}

void lutWidget_postRedisplay(struct lutWidget *self)
{
    int oldWin;
    oldWin = glutGetWindow();
    glutSetWindow(self->win);
    glutPostRedisplay();
    glutSetWindow(oldWin);
}

/* Update the frequency table graph */
int lutWidget_setFT(struct lutWidget *self, struct u8_freqtable *freqtable)
{
    int oldWin;
    oldWin=glutGetWindow();
    glutSetWindow(self->win);

    fprintf(stderr,"Resetting frequency table\n");fflush(stderr);
    self->freqtable = freqtable;
    lutWidget_makeFreqGraph(self);
    lutWidget_postRedisplay(self);

    glutSetWindow(oldWin);
    return 0;
}




/*
 * Generate the frequency graph display list,
 * from the freqtable.
 */

int lutWidget_makeFreqGraph(struct lutWidget *self)
{
    GLuint dlist;
    uint32_t freqMax;
    int i;
    unsigned int voxelval;

    /* Create a display list if one has not already been set */
    if (0==(dlist=self->freqGraphList)) {
        if (0==(dlist=self->freqGraphList=glGenLists(1))) {
            fprintf(stderr,"lutWidget_makeFreqGraph: glGenLists");
            return -1;
        }
    }
   /*
    * Find the largest population from the frequency table, *apart* from zero,
    * which we assume has an enormous population corresponding to empty space.
    */

    freqMax = self->freqtable->freq[1];
    for (i=10;i<256;i++) {
        uint32_t f;
        f=self->freqtable->freq[i];
        if (f>freqMax) { freqMax = f; }
    }

    /* Build a line graph of population. */

    glNewList(dlist,GL_COMPILE);
        glBegin(GL_LINE_STRIP);
        for (voxelval=0;voxelval<256;voxelval++) {
            GLfloat x,y;
            /*
            glColor4ub(
                    self->lut->ent[voxelval].c[0],
                    self->lut->ent[voxelval].c[1],
                    self->lut->ent[voxelval].c[2],
                    self->lut->ent[voxelval].c[3]
                    );
                    */
            x = voxelval/255.0;
            y = self->freqtable->freq[voxelval]/(float)freqMax;
            glVertex3f(x,y,0);
        }
        glEnd();
    glEndList();

    return 0;
}


int lutWidget_makeTFGraph(struct lutWidget *self)
{
    int i;

    /* Define display lists if required */
    if (0==self->tfList[0]) {
        GLuint dl;
        if (0==(dl=glGenLists(4))) {
            fprintf(stderr,"lutWidget_makeTFGraph: glGenLists");
            return -1;
        }

        for (i=0;i<4;i++) {
            self->tfList[i]=dl+i;
        }
    }

    for (i=0;i<4;i++) {
        unsigned int voxelval;

        glNewList(self->tfList[i],GL_COMPILE);
            glBegin(GL_LINE_STRIP);
                for (voxelval=0;voxelval<256;voxelval++) {
                    GLfloat x,y;
                    GLuint component;

                    x = voxelval/255.0;
                    component = self->lut->ent[voxelval].c[i];
                    y = ((float)component)/255.0;

                    if (3==i) {
                        glColor4ub(
                                self->lut->ent[voxelval].c[0],
                                self->lut->ent[voxelval].c[1],
                                self->lut->ent[voxelval].c[2],
                                255 
                               );
                    }

                    glVertex3f(x,y,0);
                }
            glEnd();
        glEndList();
    }

    return 0;
}


int lutWidget_drawTFGraph(struct lutWidget *self) {

    GLuint dl;

    if (0==(dl=self->tfList[0])) {
        fprintf(stderr,"lutWidget_drawTFGraph: no display list\n");
        return -1;
    }


    if (self->tfFlags[0]) { glColor4f(1,0,0,1); glCallList(self->tfList[0]); }
    if (self->tfFlags[1]) { glColor4f(0,1,0,1); glCallList(self->tfList[1]); }
    if (self->tfFlags[2]) { glColor4f(0,0,1,1); glCallList(self->tfList[2]); }
    if (self->tfFlags[3]) { 
        glColor4f(0.5,0.5,0.5,1); glCallList(self->tfList[3]);
    }

    return 0;

}


int lutWidget_drawFreqGraph(struct lutWidget *self)
{
    GLuint dl;

    if (0==(dl=self->freqGraphList)) {
        fprintf(stderr,"lutWidget_drawFreqGraph: no display list defined\n");
        return -1;
    }

    glColor4f(1,1,1,1);
    glCallList(dl);
    return 0;

}

/* Convert window coordinates to floats between 0 and 1 */
void lutWidget_win2ortho(struct lutWidget *self,
        int xw, int yw, GLfloat *xgp, GLfloat *ygp)
{
    *xgp = xw/(float)self->width;
    *ygp = 1.0-(yw/(float)self->height);
}

/* Convert GL ortho coordinates for a window in [0,1] to 
 * integer pixel coordinates */

void lutWidget_ortho2win(struct lutWidget *self,
        GLfloat xg, GLfloat yg, int *xwp, int *ywp)
{
    *xwp = (int)(xg*(float)self->width);
    *ywp = (int)(((float)self->height)*(1.0-yg));
}

/* Ensure that pointed-to value is clipped to [0,255] */
static void clipByte(unsigned int *x)
{ if (*x>255) { *x=255; } }

void lutWidget_setTFPoint(struct lutWidget *self,unsigned int voxelval,unsigned int tfval)
{
    int i;

    clipByte(&voxelval);
    clipByte(&tfval);

    /* both voxelval and tfval between 0 and 255 inclusive */

    for (i=0;i<4;i++) {
        if (self->tfFlags[i]) {
            self->lut->ent[voxelval].c[i] = tfval;
        }
    }
}

void lutWidget_leftClick(struct lutWidget *self,GLfloat tx,GLfloat ty)
{
    self->clickx = tx;
    self->clicky = ty;

    if (tx<0) { tx=0; } if (ty<0) { ty=0; }
    lutWidget_setTFPoint(self,tx*255.0,ty*255.0);
    lutWidget_makeTFGraph(self);

    glutPostRedisplay();

}


void lutWidget_leftDrag(struct lutWidget *self, GLfloat tx, GLfloat ty)
{

    GLfloat prevx,prevy;
    unsigned int x1,x2,y1,y2;
    unsigned int dx,voxelval;

    if (ty<0) { ty=0; }

    prevx = self->clickx;
    prevy = self->clicky;
    self->clickx = self->cursx = tx;
    self->clicky = self->cursy = ty;



    /* Ensure x1<x2 */

    if (prevx<tx) {
        x1 = (unsigned int)(255.0*prevx);
        y1 = (unsigned int)(255.0*prevy);
        x2 = (unsigned int)(255.0*tx);
        y2 = (unsigned int)(255.0*ty);
    } else {
        x2 = (unsigned int)(255.0*prevx);
        y2 = (unsigned int)(255.0*prevy);
        x1 = (unsigned int)(255.0*tx);
        y1 = (unsigned int)(255.0*ty);
    }

    clipByte(&x1);
    clipByte(&x2);
    clipByte(&y1);
    clipByte(&y2);

    dx = x2-x1;


    if (0==dx) { lutWidget_leftClick(self,tx,ty); return; }

    for (voxelval=x1;voxelval<=x2;voxelval++) {
        GLfloat interp;
        unsigned tfval;

        interp = (voxelval-x1)/(float)dx;
        tfval=(unsigned int) ( (1.0-interp)*(float)y1 + interp*(float)y2);
        clipByte(&tfval);

        lutWidget_setTFPoint(self,voxelval,tfval);
    }

    lutWidget_makeTFGraph(self);

    glutPostRedisplay();

}

void lutWidget_motionFunc(int x, int y)
{
    struct lutWidget *self=NULL;
    GLfloat tx,ty;

    if (NULL==(self=winReg_getPtr(0))) {
        fprintf(stderr,"lutWidget_motionFunc: winReg_getPtr failed\n");
        return;
    }

    /* deal with floating-point coordinates wherever possible */
    lutWidget_win2ortho(self,x,y,&tx,&ty);


    if (self->buttonState[LUTWIDGET_LEFT_BUTTON]) {
        lutWidget_leftDrag(self,tx,ty);
    }

    glutPostRedisplay();

}



void lutWidget_passiveMotionFunc(int x, int y)
{
    GLfloat cx,cy;
    struct lutWidget *self=NULL;
    if (NULL==(self=winReg_getPtr(0))) {
        fprintf(stderr,"lutWidget_passiveMotionFunc: winReg_getPtr failed\n");
        return;
    }

    lutWidget_win2ortho(self,x,y,&cx,&cy);

    self->cursx = cx;
    self->cursy = cy;

    glutPostRedisplay(); 
}

void lutWidget_reshapeFunc(int w, int h)
{
    struct lutWidget *self=NULL;
    if (NULL==(self=winReg_getPtr(0))) {
        fprintf(stderr,"lutWidget_reshapeFunc: winReg_getPtr failed\n");
        return;
    }
    glViewport(0,0,w,h);
    self->width = w;
    self->height = h;
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0,1,0,1,0,1);
    glMatrixMode(GL_MODELVIEW);
}

void lutWidget_wheelUp(struct lutWidget *self)
{ return; }
void lutWidget_wheelDown(struct lutWidget *self)
{ return; }

void lutWidget_mouseFunc(int button, int state, int x, int y)
{
    int butno=0;
    struct lutWidget *self=NULL;

    if (NULL==(self=winReg_getPtr(0))) {
        fprintf(stderr,"lutWidget_mouseFunc: winReg_getPtr failed\n");
        return;
    }

    switch (button) {
        case GLUT_LEFT_BUTTON:
            butno=LUTWIDGET_LEFT_BUTTON; break;
        case GLUT_MIDDLE_BUTTON:
            butno=LUTWIDGET_MIDDLE_BUTTON; break;
        case GLUT_RIGHT_BUTTON:
            butno=LUTWIDGET_RIGHT_BUTTON; break;
        case WHEEL_UP:
            if (GLUT_DOWN==state) { lutWidget_wheelUp(self); }
            break;
        case WHEEL_DOWN:
            if (GLUT_DOWN==state) { lutWidget_wheelDown(self); }
            break;
        default:
            /* unknown button; do nothing but return silently. */
            return;
    }

    self->buttonState[butno] = (GLUT_DOWN==state) ? 1 : 0;

    if (self->buttonState[LUTWIDGET_LEFT_BUTTON]) {
        GLfloat tx,ty;

        lutWidget_win2ortho(self,x,y,&tx,&ty);
        lutWidget_leftClick(self,tx,ty);
        glutPostRedisplay();
    }
}

void lutWidget_drawCursor(struct lutWidget *self)
{

    glPushMatrix();
        glLoadIdentity();

        glLineWidth(1.0);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
        glTranslatef(self->cursx,self->cursy,0);

        if (self->cursorList) {
            glCallList(self->cursorList);
        }
    glPopMatrix();
}


void lutWidget_newCursor(struct lutWidget *self)
{
#define LUTWIDGET_CURSOR_ALPHA 1.0
    const GLfloat crossSize=0.02;
    const GLfloat rgbaWidth=0.01;
    int i,j;

    GLuint dlist;
    vector crossCoords[4];
    vector squareCoords[5];
    vector offset;
    GLfloat color[4][4] = {
            {1,0,0,LUTWIDGET_CURSOR_ALPHA},
            {0,1,0,LUTWIDGET_CURSOR_ALPHA},
            {0,0,1,LUTWIDGET_CURSOR_ALPHA},
            {1,1,1,LUTWIDGET_CURSOR_ALPHA}
    };

    crossCoords[0] = vector_new(-crossSize,0,0);
    crossCoords[1] = vector_new( crossSize,0,0);
    crossCoords[2] = vector_new(0,-crossSize,0);
    crossCoords[3] = vector_new(0, crossSize,0);

    squareCoords[0] = vector_new(0,0,0);
    squareCoords[1] = vector_new(rgbaWidth,0,0);
    squareCoords[2] = vector_new(rgbaWidth,-rgbaWidth,0);
    squareCoords[3] = vector_new(0,-rgbaWidth,0);
    squareCoords[4] = vector_new(0,0,0);

    offset = vector_new(crossSize,crossSize,0);

    if (0==(dlist=self->cursorList)) {
        if (0==(dlist=self->cursorList=glGenLists(1))) {
            fprintf(stderr,"lutWidget_newCursor: glGenLists()=0\n");
            return;
        }
    }

    glNewList(dlist,GL_COMPILE);

        glColor4f(1,1,1,LUTWIDGET_CURSOR_ALPHA);
        glBegin(GL_LINES);
            for (i=0;i<4;i++) {
                vector v;
                v=crossCoords[i];
                glVertex3f(v.x,v.y,v.z);
            }
        glEnd();

        for (i=0;i<4;i++) {
            glColor4fv(&(color[i][0]));

            if (self->tfFlags[i]) {
                glBegin(GL_POLYGON);
            } else {
                glBegin(GL_LINE_STRIP);
            }

            for (j=0;j<5;j++) {
                vector v;
                v = vector_add(squareCoords[j],offset);
                glVertex3f(v.x,v.y,v.z);
            }
            glEnd();

            offset.y -= rgbaWidth;
        }

    glEndList();
}

void lutWidget_changeCursor(struct lutWidget *self)
{
    lutWidget_newCursor(self);
    glutPostRedisplay();
}

void lutWidget_toggleTFComponent(struct lutWidget *self, int component)
{
    int prev;

    prev = self->tfFlags[component];
    if (prev) { self->tfFlags[component]=0; } 
    else { self->tfFlags[component]=1; }
    lutWidget_changeCursor(self);
}

void lutWidget_keyboardFunc(unsigned char key, int x, int y)
{
    struct lutWidget *self=NULL;

    if (NULL==(self=winReg_getPtr(0))) {
        fprintf(stderr,"lutWidget_keyboardFunc: winReg_getPtr failed\n");
        return;
    }


    switch(tolower(key)) {
        case 'r':
            lutWidget_toggleTFComponent(self,0);
            break;
        case 'g':
            lutWidget_toggleTFComponent(self,1);
            break;
        case 'b':
            lutWidget_toggleTFComponent(self,2);
            break;
        case 'a':
            lutWidget_toggleTFComponent(self,3);
            break;

    }

    if (NULL!=self->keyCallback) {
        self->keyCallback(self,key,self->keyCallbackData);
    }
}


void lutWidget_displayFunc(void)
{
    struct lutWidget *self=NULL;

    if (NULL==(self=winReg_getPtr(0))) {
        fprintf(stderr,"lutWidget_displayFunc: winReg_getPtr failed\n");
        return;
    }

    glClearColor(0,0,0,1);
    glClear(GL_COLOR_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0,1,0,1,0,1);
    glMatrixMode(GL_MODELVIEW);

    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_BLEND);
    glEnable(GL_LINE_SMOOTH);
    glLineWidth(2.0);

    glLoadIdentity();
    glPushMatrix();
        lutWidget_drawFreqGraph(self);
        lutWidget_drawTFGraph(self);
    glPopMatrix();


    glLineWidth(1.0); glDisable(GL_BLEND);
    if (self->cursorOn) {
        lutWidget_drawCursor(self);
    }
    
    glutSwapBuffers();
}

void lutWidget_entryFunc(int state)
{
    struct lutWidget *self=NULL;

    if (NULL==(self=winReg_getPtr(0))) {
        fprintf(stderr,"lutWidget_entryFunc: winReg_getPtr failed\n");
        return;
    }

    if (GLUT_ENTERED==state) {
        self->cursorOn = 1;
    } else {
        self->cursorOn = 0;
    }
    glutPostRedisplay();



}


int lutWidget_initWindow(struct lutWidget *self)
{

    glutInitDisplayMode(GLUT_RGBA|GLUT_DOUBLE);
    glutInitWindowSize(self->width,self->height);
    self->win = glutCreateWindow("LUT widget");
    glutSetCursor(GLUT_CURSOR_NONE);

    winReg_setPtr(0,self);

    glutMotionFunc(lutWidget_motionFunc);
    glutPassiveMotionFunc(lutWidget_passiveMotionFunc);
    glutReshapeFunc(lutWidget_reshapeFunc);
    glutMouseFunc(lutWidget_mouseFunc);
    glutKeyboardFunc(lutWidget_keyboardFunc);
    glutEntryFunc(lutWidget_entryFunc);
    glutDisplayFunc(lutWidget_displayFunc);

    return 0;
}

void lutWidget_init(struct lutWidget *widget)
{
    lutWidget_initWindow(widget);
    lutWidget_makeFreqGraph(widget);
    lutWidget_makeTFGraph(widget);
    lutWidget_newCursor(widget);
}


struct lutWidget *lutWidget_new(struct u8_freqtable *freqtable,struct u8_rgba_lut *lut)
{
    struct lutWidget *widget=NULL;
    int i;

    if (NULL==(widget=malloc(sizeof(struct lutWidget))))
        { perror("lutWidget_new: malloc"); return NULL; }

    widget->freqtable = freqtable;
    widget->lut = lut;
    widget->width = 768;
    widget->height= 384;
    widget->cursorOn = 1;
    widget->tfList[0]=0;
    widget->keyCallback = NULL;
    widget->keyCallbackData = NULL;

    for (i=0;i<2;i++) { widget->tfFlags[i]=0; }
    widget->tfFlags[3]=1; /* start off with only alpha showing */

    widget->clickx = widget->clicky = 0;
    widget->cursx = widget->cursy = 0;

    for (i=0;i<3;i++) { widget->buttonState[i]=0; }

    lutWidget_init(widget);

    return widget;
}

void lutWidget_setKeyCallback(struct lutWidget *self,
        void (*callback)(struct lutWidget *,unsigned char,void *),void *data)
{
    self->keyCallback = callback;
    self->keyCallbackData = data;
}

