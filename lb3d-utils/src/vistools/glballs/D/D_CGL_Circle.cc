#include <math.h>
#include <GL/glu.h>

#include "D_CGL.h"


#define D_GLCIRCLE(STEP,ALPHA) {\
 glBegin(GL_POLYGON);\
  for(int v=0;v<(STEP);v++)glVertex3f(pArray[i].r*std::sin(v*(ALPHA)),pArray[i].r*std::cos(v*(ALPHA)),0.);\
 glEnd();\
 glColor3f(0.,0.,0.);\
 glBegin(GL_LINE_LOOP);\
  for(int v=0;v<(STEP);v++)glVertex3f(pArray[i].r*std::sin(v*(ALPHA)),pArray[i].r*std::cos(v*(ALPHA)),0.0001);\
 glEnd();\
}

#define DISPLAYMODE_CASE_FUNCTION void D_CGL::updateListCircle
#define DISPLAYMODE_CASE(a,b) case a : updateListCircle ## b(pHeader,pArray); break;

//
// Anfang jedes Funktionsaufrufes
//
#define DISPLAYMODE_FUNCTION_BEGIN(colormode) \
void D_CGL::updateListCircle ## colormode(D_Header_t *pHeader,D_Data_t *pArray)\
{\
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
  glColorMaterial(GL_FRONT,GL_AMBIENT);\
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
    glRotatef( -m_zRot, 0.0, 0.0, 1.0 );\
    glRotatef( -m_yRot, 0.0, 1.0, 0.0 ); \
    glRotatef( -m_xRot, 1.0, 0.0, 0.0 ); \
    if(pArray[i].r<m_radius)D_GLCIRCLE(8,0.78539816)\
    else D_GLCIRCLE(16,0.39269908)\
   glPopMatrix();\
  }\
  glDisable(GL_COLOR_MATERIAL);\
 glEndList();\
}\


#include "D_CGL_ColorModes_include.h"

