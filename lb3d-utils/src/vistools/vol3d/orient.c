#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include "quat.h"
#include "vector.h"
#include "winReg.h"
#include "glutViewer.h"

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

void myDisplay(struct glutViewer *gv, void *data)
{
    GLfloat modelview[16];
    GLfloat m[9];
    int i,j;
    vector cubeForward,cubeUp,cubeRight;

    glGetFloatv(GL_MODELVIEW_MATRIX,modelview);

    /* Extract rotation component */

    for (i=0;i<3;i++) {
    for (j=0;j<3;j++) {
        m[3*i+j] = modelview[4*i+j];
    }}

    getNearestCubeVecs(&cubeForward,&cubeUp,m);

    cubeRight = vector_cross(cubeForward,cubeUp);

    glClearColor(0,0,0,0);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

    glEnable(GL_BLEND);
    glDisable(GL_LIGHTING);
    glLineWidth(3.0);

    glColor4f(0,1,0,0.5);
    glutWireCube(1.0);

    glLineWidth(5.0);
    glBegin(GL_LINES);

        glColor4f(1,0,1,0.5);
        glVertex3f(0,0,0);
        glVertex_vec(cubeForward);

        glColor4f(0,0,1,0.5);
        glVertex3f(0,0,0);
        glVertex_vec(cubeUp);

        glColor4f(1,0,0,0.5);
        glVertex3f(0,0,0);
        glVertex_vec(cubeRight);

    glEnd();


    glutSwapBuffers();
}

void keyFunc(unsigned char key, int x, int y)
{
    struct glutViewer *gv=NULL;

    if (NULL==(gv=winReg_getPtr(0))) {
        fprintf(stderr,"volRender_keyboardFunc: winReg_getPtr failed\n");
        return;
    }

    if ('S'==key) {
        snapToCube(gv);
    }
}

int main(int argc, char *argv[])
{

    struct glutViewer *gv=NULL;

    initCubeVectors();

    glutInit(&argc,argv);

    if (NULL==(gv=glutViewer_new()))
    { fprintf(stderr,"glutViewer_new() failed.\n"); return -1; }

    gv->fovy=20.0;

    glutViewer_setDisplayCallback(gv,myDisplay);
    glutViewer_setDisplayCallbackData(gv,
            "This is the display callback callback.\n");

    glutKeyboardFunc(keyFunc);
    glutMainLoop();

    return 0;
}       
