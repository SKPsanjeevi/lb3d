#include <math.h>
#include <stdio.h>
#include <GL/glu.h>

#include "D_CGL.h"

#define DISPLAYMODE_CASE_FUNCTION void D_CGL::updateListCube
#define DISPLAYMODE_CASE(a,b) case a : updateListCube ## b(pHeader,pArray); break;

// 2.*360/(2*m_pi)
#define ACosToThetaGrad 114.59156

// 
// Anfang jedes Funktionsaufrufes
//
#define DISPLAYMODE_FUNCTION_BEGIN(colormode) \
void D_CGL::updateListCube ## colormode(D_Header_t *pHeader,D_Data_t *pArray)\
{\
\
 GLfloat OffsetX,OffsetY,OffsetZ;\
\
 OffsetX=0.5*(pHeader->xMax+pHeader->xMin);\
 OffsetY=0.5*(pHeader->yMax+pHeader->yMin);\
 OffsetZ=0.5*(pHeader->zMax+pHeader->zMin);\
    
//
// Hier wird Sourcecode fuer die Farbbestimmung eingefuegt
//

//
// Anfang der Arbeitschleife
//

#define DISPLAYMODE_FUNCTION_MIDDLE \
 glNewList(m_glList, GL_COMPILE );\
  glEnable(GL_COLOR_MATERIAL);\
  glColorMaterial(GL_FRONT,GL_AMBIENT_AND_DIFFUSE);\
  for(long int i=0;i<pHeader->Number;++i){\
   glPushMatrix();\

//
// Hier wird Sourcecode fuer die Farbbestimmung eingefuegt
//


//
// Ende der Arbeitschleife und Abschluss der Funktion
//
// s=0.80599598*... => W"urfel hat das selbe Volumen der Kugel
#define DISPLAYMODE_FUNCTION_END \
    float m[16];\
    float xx=pArray[i].qx*pArray[i].qx;\
    float yy=pArray[i].qy*pArray[i].qy;\
    float zz=pArray[i].qz*pArray[i].qz;\
    float wx=pArray[i].qw*pArray[i].qx;\
    float wy=pArray[i].qw*pArray[i].qy;\
    float wz=pArray[i].qw*pArray[i].qz;\
    float xy=pArray[i].qx*pArray[i].qy;\
    float xz=pArray[i].qx*pArray[i].qz;\
    float yz=pArray[i].qy*pArray[i].qz;\
    m[ 0]=1.-2.*(yy+zz);m[ 4]=2.*(xy-wz);m[ 8]=2.*(xz+wy);m[12]=pArray[i].x-OffsetX;\
    m[ 1]=2.*(xy+wz);m[ 5]=1.-2.*(xx+zz);m[ 9]=2.*(yz-wx);m[13]=pArray[i].y-OffsetY;\
    m[ 2]=2.*(xz-wy);m[ 6]=2.*(yz+wx);m[10]=1.-2.*(xx+yy);m[14]=pArray[i].z-OffsetZ;\
    m[ 3]=0.00;m[ 7]=0.00;m[11]=0.00;m[15]=1.00;\
    glMultMatrixf(m);\
    GLfloat s;\
    s=0.80599598*pArray[i].r;\
    glBegin(GL_QUADS);\
     glNormal3f(0.,0.,1.);\
     glVertex3f(-s,s,s);\
     glVertex3f(s,s,s);\
     glVertex3f(s,-s,s);\
     glVertex3f(-s,-s,s);\
     glNormal3f(0.,1.,0.);\
     glVertex3f(s,s,s);\
     glVertex3f(s,s,-s);\
     glVertex3f(-s,s,-s);\
     glVertex3f(-s,s,s);\
     glNormal3f(0.,0.,-1.);\
     glVertex3f(-s,s,-s);\
     glVertex3f(s,s,-s);\
     glVertex3f(s,-s,-s);\
     glVertex3f(-s,-s,-s);\
     glNormal3f(0.,-1.,0.);\
     glVertex3f(s,-s,s);\
     glVertex3f(s,-s,-s);\
     glVertex3f(-s,-s,-s);\
     glVertex3f(-s,-s,s);\
     glNormal3f(1.,0.,0.);\
     glVertex3f(s,-s,s);\
     glVertex3f(s,-s,-s);\
     glVertex3f(s,s,-s);\
     glVertex3f(s,s,s);\
     glNormal3f(-1.,0.,0.);\
     glVertex3f(-s,-s,s);\
     glVertex3f(-s,-s,-s);\
     glVertex3f(-s,s,-s);\
     glVertex3f(-s,s,s);\
    glEnd(); \
   glPopMatrix();\
  }\
  glDisable(GL_COLOR_MATERIAL);\
 glEndList();\
}\

#include "D_CGL_ColorModes_include.h"

