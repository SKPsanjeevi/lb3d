#include "D_CGL.h"

/*
*/
void D_CGL::initBorders(int i)
{
 GL_CQuadIO quadio(m_glListFixedFilename[i]);
 int numberofquads=quadio.readnumberofquads();
 GL_CQuad backquad;
  
  glNewList(m_glListFixed[i], GL_COMPILE );
  
  glEnable(GL_COLOR_MATERIAL);
  glColorMaterial(GL_FRONT,GL_AMBIENT);

  glBegin(GL_QUADS);
  
  for (int i=0;i<numberofquads;++i){

    backquad=quadio.readquad();


    glColor3f(backquad.rgb[0],backquad.rgb[1],backquad.rgb[2]);

    glNormal3fv(backquad.norm);
    glVertex3f(backquad.vertex[0],backquad.vertex[1],backquad.vertex[2]);
    glVertex3f(backquad.vertex[3],backquad.vertex[4],backquad.vertex[5]);
    glVertex3f(backquad.vertex[6],backquad.vertex[7],backquad.vertex[8]);
    glVertex3f(backquad.vertex[9],backquad.vertex[10],backquad.vertex[11]);
  }
  glDisable(GL_COLOR_MATERIAL);
  
  glEnd();


  glEndList();
}
