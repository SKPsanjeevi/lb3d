#ifndef GL_CQUAD
#define GL_CQUAD

#include <cmath>
#include "T_CVector.h"

struct GL_CQuad{
  inline void set1(float x1,float y1, float z1) {vertex[0]=x1;vertex[1]=y1;vertex[2]=z1;}
  inline void set2(float x2,float y2, float z2) {vertex[3]=x2;vertex[4]=y2;vertex[5]=z2;}
  inline void set3(float x3,float y3, float z3) {vertex[6]=x3;vertex[7]=y3;vertex[8]=z3;}
  inline void set4(float x4,float y4, float z4) {vertex[9]=x4;vertex[10]=y4;vertex[11]=z4;}
  inline void setnorm(float x, float y, float z){norm[0]=x;norm[1]=y;norm[2]=z;}
  inline void set1(const T_CVector &CV) {vertex[0]=CV.getX();vertex[1]=CV.getY();vertex[2]=CV.getZ();}
  inline void set2(const T_CVector &CV) {vertex[0+3]=CV.getX();vertex[1+3]=CV.getY();vertex[2+3]=CV.getZ();}
  inline void set3(const T_CVector &CV) {vertex[0+6]=CV.getX();vertex[1+6]=CV.getY();vertex[2+6]=CV.getZ();}
  inline void set4(const T_CVector &CV) {vertex[0+9]=CV.getX();vertex[1+9]=CV.getY();vertex[2+9]=CV.getZ();}
  inline void setnorm(const T_CVector &CV){norm[0]=CV.getX();norm[1]=CV.getY();norm[2]=CV.getZ();}
  inline void setrgb(char r, char g, char b) {rgb[0]=r;rgb[1]=g;rgb[2]=b;}
  float vertex[12];
  float norm[3];
  char rgb[3];
  //    double r,g,b;
  
};

#endif
