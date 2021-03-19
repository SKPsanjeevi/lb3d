/* Compile me with:
 * cc -c -o glutViewer.o glutViewer.c  
 *
 * Wants to link with quat.c , winReg.c , and libglut.
 *
 * Boilerplate setup code to open a window using GLUT,
 * and render an object in it.
 *
 * 2002,2003,2005 Jonathan Chin <jon-src@earth.li>
 */


/*
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <unistd.h>
#include <math.h>
#include "quat.h"
#include "winReg.h"
#include "glutViewer.h"



/* Process mouse motion from x0,y0 to x1,y1 ; update orientation quaternion. */

void glutViewer_mouseRotate(struct glutViewer *self,
        int x0,int y0, int x1,int y1) {

        GLfloat mx,my,m,s;
        GLfloat theta;
        quat_t rotquat;

        s = self->sphererad;
        my = (GLfloat)(x1-x0);
        mx = (GLfloat)(y1-y0);
        m=sqrt(mx*mx+my*my);

        if ((m>0) && (m<s)) {
                theta = m/s;

                mx /= m;
                my /= m;

                rotquat = quatrotation(theta,mx,my,0.0);
                self->orquat = quatmultiply(rotquat,self->orquat);
        }

}

void glutViewer_mouseZoom(struct glutViewer *self,int dz)
{
        const GLfloat zoomscale=0.03;
        const GLfloat orthoScale = 0.005;
        if (self->orthographic) {
            self->orthoSize += orthoScale * dz;
        } else {
            self->geomz-=zoomscale*dz;
        }
        return;
}

void glutViewer_mouseTranslate(struct glutViewer *self,int dx,int dy)
{
        const GLfloat transscale=0.0025;
        self->geomx-=transscale*dx;
        self->geomy+=transscale*dy; /* GL uses left-handed coordinates */
        return;
}


void glutViewer_reshapeFunc(int width, int height)
{
    struct glutViewer *self=NULL;

    if (NULL==(self=winReg_getPtr(0))) {
        fprintf(stderr,"glutViewer_reshapeFunc: winReg_getPtr failed\n");
        return;
    }

    self->width = width;
    self->height = height;
    glViewport(0,0,self->width,self->height);
    if (self->width<self->height) {
            self->sphererad=0.5*self->height;
    } else {
            self->sphererad=0.5*self->width;
    }
}

void glutViewer_displayFunc(void)
{
    struct glutViewer *self=NULL;
    GLfloat m[16];

    if (NULL==(self=winReg_getPtr(0))) {
        fprintf(stderr,"glutViewer_displayFunc: winReg_getPtr failed\n");
        return;
    }
    
    if (self->stereoFlag) { glDrawBuffer(GL_BACK); }

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    if (self->orthographic) {
        GLfloat size=self->orthoSize;
      glOrtho(-size,size,-size,size,self->nearClip,64.0);
    } else {
        /* Set for perspective projection */

        gluPerspective(self->fovy,(float)self->width/(float)self->height,
                self->nearClip,1024.0 );
    }
    glMatrixMode(GL_MODELVIEW);

    glDepthMask(GL_TRUE);
    glClearColor(
            self->bgColor[0],
            self->bgColor[1],
            self->bgColor[2],
            self->bgColor[3]);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

    glLoadIdentity();

    if (self->stereoFlag) {
        quatmatrix(self->orquat,m);
        glPushMatrix();
            glDrawBuffer(GL_BACK_LEFT);
            glTranslatef(self->geomx,self->geomy,self->geomz);
            glTranslatef(self->eyesize,0,0);

            glMultMatrixf(m);

            if (self->displayCallback) {
                self->displayCallback(self,self->data);
            }
        glPopMatrix();
        glPushMatrix();
            glDrawBuffer(GL_BACK_RIGHT);
            glClear(GL_DEPTH_BUFFER_BIT);
            glTranslatef(self->geomx,self->geomy,self->geomz);
            glTranslatef(-self->eyesize,0,0);

            glMultMatrixf(m);

            if (self->displayCallback) {
                self->displayCallback(self,self->data);
            }
        glPopMatrix();

    } else {

        /* Position and orient geometry */
        glTranslatef(self->geomx,self->geomy,self->geomz);
        quatmatrix(self->orquat,m);
        glMultMatrixf(m);

        if (self->displayCallback) {
            self->displayCallback(self,self->data);
        }
    }
    glFinish();
    glutSwapBuffers();
}

void glutViewer_mouseFunc(int button, int state, int x, int y)
{
    struct glutViewer *self=NULL;

    if (NULL==(self=winReg_getPtr(0))) {
        fprintf(stderr,"glutViewer_mouseFunc: winReg_getPtr failed\n");
        return;
    }
    self->clickx = x;
    self->clicky = y;

#define WHEEL_UP 3      /* Experimental measurement.. */
#define WHEEL_DOWN 4

    switch (button) {
        case GLUT_LEFT_BUTTON:
            self->buttonState[0] = (GLUT_DOWN == state) ? 1 : 0;
            break;
        case GLUT_MIDDLE_BUTTON:
            self->buttonState[1] = (GLUT_DOWN == state) ? 1 : 0;
            break;
        case GLUT_RIGHT_BUTTON:
            self->buttonState[2] = (GLUT_DOWN == state) ? 1 : 0;
            break;
        case WHEEL_UP:
            if (GLUT_DOWN==state) {
                if (NULL!=self->wheelUpFunc) {
                    self->wheelUpFunc(self,self->wheelUpData);
                }
            }
            break;
        case WHEEL_DOWN:
            if (GLUT_DOWN==state) {
                if (NULL!=self->wheelDownFunc) {
                    self->wheelDownFunc(self,self->wheelDownData);
                }
            }
            break;
    }
}

void glutViewer_motionFunc(int x, int y)
{
    struct glutViewer *self=NULL;

    if (NULL==(self=winReg_getPtr(0))) {
        fprintf(stderr,"glutViewer_mouseFunc: winReg_getPtr failed\n");
        return;
    }

    if (self->buttonState[0]) {
        /* Left button */
            glutViewer_mouseRotate(self,self->clickx,self->clicky,x,y);
    } else if (self->buttonState[1]) {
        /* Middle button */
            glutViewer_mouseTranslate(self,self->clickx-x,self->clicky-y);
    } else {
        /* Assume right button */
            glutViewer_mouseZoom(self,y-self->clicky);
    }
    self->clickx=x;self->clicky=y;
    glutPostRedisplay();
}

int glutViewer_init(struct glutViewer *self)
{
    unsigned int flags = GLUT_RGBA|GLUT_DOUBLE|GLUT_DEPTH;
    char *title="GLUT viewer";

    if (self->alphaBufferFlag) { flags |= GLUT_ALPHA; }
    if (self->stereoFlag) { flags |= GLUT_STEREO; title = "GLUT stereo viewer";}
    glutInitDisplayMode(flags);
    glutInitWindowSize(self->width,self->height);
    glutCreateWindow(title);

    winReg_setPtr(0,self);
    self->win = glutGetWindow();

    glutMotionFunc(glutViewer_motionFunc);
    glutReshapeFunc(glutViewer_reshapeFunc);
    glutMouseFunc(glutViewer_mouseFunc);
    glutDisplayFunc(glutViewer_displayFunc);

    /* Set up OpenGL state */

    glViewport(0,0,self->width,self->height);

    glShadeModel(GL_SMOOTH);


    glClearDepth(1.0);

    glDepthFunc(GL_LESS);
    glDepthMask(GL_TRUE);   /* Enable depth buffer writes */
    glEnable(GL_DEPTH_TEST);

    /*
    glCullFace(GL_BACK);
    glEnable(GL_CULL_FACE);

    glDisable(GL_TEXTURE_2D);
    */

    return 0;
}


struct glutViewer *glutViewer_new_minimal(void)
{
    struct glutViewer *gv=NULL;
    int i;

    if (NULL==(gv=malloc(sizeof(struct glutViewer))))
        { perror("glutViewer_new: malloc"); return NULL; }

    gv->width = 512;
    gv->height = 512;
    gv->cursx = gv->cursy = gv->cursx = gv->cursy = 0;
    gv->orthographic=0;
    gv->orthoSize = 0.3;

    gv->fovy = 20.0;

    for (i=0;i<3;i++) { gv->buttonState[i]=0; }
    gv->stereoFlag= 0;
    gv->eyesize=0.1;

    gv->bgColor[0] = 0.3;
    gv->bgColor[1] = 0.3;
    gv->bgColor[2] = 0.5;
    gv->bgColor[4] = 0.0;

    gv->nearClip = 0.01;
    gv->geomx = gv->geomy = 0;
    gv->geomz = -5;

    gv->orquat = quatrotation(0,0,0,1);

    gv->sphererad = 256;

    gv->data = NULL;
    gv->displayCallback = NULL;
    gv->displayCallbackData = NULL;

    gv->alphaBufferFlag = 0;

    return gv;

}

struct glutViewer *glutViewer_new(void)
{
    struct glutViewer *gv=NULL;

    if (NULL==(gv=glutViewer_new_minimal())) {
        fprintf(stderr,"glutViewer_new: glutViewer_new_minimal failed.\n");
        return NULL;
    }
    glutViewer_init(gv);
    return gv;
}

struct glutViewer *glutViewer_new_stereo(void)
{
    struct glutViewer *gv=NULL;

    if (NULL==(gv=glutViewer_new_minimal())) {
        fprintf(stderr,"glutViewer_new_stereo: glutViewer_new_minimal failed.\n");
        return NULL;
    }
    gv->stereoFlag = 1;
    glutViewer_init(gv);
    return gv;
}

int glutViewer_setData(struct glutViewer *self, void *data)
{ self->data=data; return 0; }

void *glutViewer_getData(struct glutViewer *self)
{ return self->data; }

quat_t glutViewer_getQuat(struct glutViewer *self)
{ return self->orquat; }

void glutViewer_wheelUpFunc(struct glutViewer *self, void (*func)(struct glutViewer *, void *))
{
    self->wheelUpFunc = func;
}

void glutViewer_wheelUpData(struct glutViewer *self, void *data)
{
    self->wheelUpData = data;
}

void glutViewer_wheelDownFunc(struct glutViewer *self, void (*func)(struct glutViewer *, void *))
{
    self->wheelDownFunc = func;
}

void glutViewer_wheelDownData(struct glutViewer *self, void *data)
{
    self->wheelDownData = data;
}

void (*glutViewer_setDisplayCallback(struct glutViewer *self, void (*newCallback)(struct glutViewer *, void *)))(struct glutViewer *, void *)
{
    void (*old)(struct glutViewer *, void *);

    old = self->displayCallback;
    self->displayCallback = newCallback;
    return old;
}

void (*glutViewer_getDisplayCallback(
            struct glutViewer *self))(struct glutViewer *, void *)
{ return self->displayCallback; }

void *glutViewer_setDisplayCallbackData(struct glutViewer *self, void *data)
{
    void *old;

    old = self->data;
    self->data = data;
    return old;
}

void *glutViewer_getDisplayCallbackData(struct glutViewer *self, void *data)
{ return self->displayCallbackData; }

void glutViewer_postRedisplay(struct glutViewer *self)
{
    int oldWin;
    oldWin = glutGetWindow();
    glutSetWindow(self->win);
    glutPostRedisplay();
    glutSetWindow(oldWin);
}

void glutViewer_toggleOrtho(struct glutViewer *self)
{
    if (self->orthographic) {
        self->orthographic = 0;
    } else {
        self->orthographic = 1;
    }

    glutViewer_postRedisplay(self);
}
