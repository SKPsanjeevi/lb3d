#include <math.h>
#include <stdio.h>
#include <GL/glu.h>

#include "D_CGL.h"

#define DISPLAYMODE_CASE_FUNCTION void D_CGL::updateListSphere
#define DISPLAYMODE_CASE(a,b) case a : updateListSphere ## b(pHeader,pArray); break;

//
// Anfang jedes Funktionsaufrufes 
//
#define DISPLAYMODE_FUNCTION_BEGIN(colormode) \
void D_CGL::updateListSphere ## colormode(D_Header_t *pHeader,D_Data_t *pArray)\
{\
\
 GLUquadricObj *pObj;\
 GLfloat OffsetX,OffsetY,OffsetZ;\
\
 OffsetX=0.5*(pHeader->xMax+pHeader->xMin);\
 OffsetY=0.5*(pHeader->yMax+pHeader->yMin);\
 OffsetZ=0.5*(pHeader->zMax+pHeader->zMin);\
 m_HeadOffsetX=OffsetX;\
 m_HeadOffsetY=OffsetY;\
 m_HeadOffsetZ=OffsetZ;\


//
// Hier wird Sourcecode fuer die Farbbestimmung eingefuegt
//
//glDepthMask(GL_FALSE);
//glDepthMask(GL_TRUE);

//
// Anfang der Arbeitschleife
//
//  glColorMaterial(GL_FRONT,GL_AMBIENT);

#define DISPLAYMODE_FUNCTION_MIDDLE \
 pObj=gluNewQuadric();\
 gluQuadricDrawStyle(pObj,(GLenum)GLU_FILL);\
 gluQuadricNormals(pObj,(GLenum)GLU_SMOOTH);\
 glNewList(m_glListParticle, GL_COMPILE );\
 gluSphere(pObj,1.,15,10);\
 glEndList();\
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
#define DISPLAYMODE_FUNCTION_END \
    glTranslatef(pArray[i].x-OffsetX,pArray[i].y-OffsetY,pArray[i].z-OffsetZ);\
    if(pArray[i].r<m_radius)gluSphere(pObj,pArray[i].r,8,4);\
    else gluSphere(pObj,pArray[i].r,15,10);\
   glPopMatrix();\
  }\
  glDisable(GL_COLOR_MATERIAL);\
 glEndList();\
 gluDeleteQuadric(pObj);\
}\


#include "D_CGL_ColorModes_include.h"

