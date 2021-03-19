#include <stdlib.h>
#include <stdio.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "quat.h"
#include "winReg.h"
#include "glutViewer.h"

void myDisplay(struct glutViewer *gv, void *data)
{
    if (NULL!=data) {
        printf("<%s>\n",(char *)data);fflush(stdout);
        glutViewer_setDisplayCallbackData(gv,NULL);
    }
    glClearColor(0,0,0,0);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    glDepthMask(GL_FALSE);
    glDisable(GL_DEPTH_TEST);
    glColor4f(0,1,0,0.9);
    glutWireCube(1.0);
    glColor4f(1,1,1,1);
    glPointSize(10.0);
    glBegin(GL_POINTS);
    glVertex3f(0,0,0);
    glEnd();
    glutSwapBuffers();
}

int main(int argc, char *argv[])
{

    struct glutViewer *gv=NULL;

    glutInit(&argc,argv);

    if (NULL==(gv=glutViewer_new()))
    { fprintf(stderr,"glutViewer_new() failed.\n"); return -1; }

    glutViewer_setDisplayCallback(gv,myDisplay);
    glutViewer_setDisplayCallbackData(gv,
            "This is the display callback callback.\n");

    glutMainLoop();

    return 0;
}       
