#ifndef D_CGL_H
#define D_CGL_H

#include <math.h>
#include <GL/glu.h>

#include "T_CArg.h"
#include "D_CData.h"
#include "GL_CQuadIO.h"

#define DISPLAYMODE_SPHERE 0
#define DISPLAYMODE_TEXTURE 1
#define DISPLAYMODE_CUBE 2
#define DISPLAYMODE_CIRCLE 3
#define DISPLAYMODE_BLIP 4

#define COLORMODE_VELOCITY 0
#define COLORMODE_VELOCITYFIXED 1
#define COLORMODE_ANGULARVELOCITY 2
#define COLORMODE_VELOCITYNORM 3
#define COLORMODE_VELOCITYNORMFIXED 4
#define COLORMODE_TYPE 5
#define COLORMODE_ID 6
#define COLORMODE_COLOR 7
#define COLORMODE_FOREGROUND 8

#define NUMBER_OF_BORDERS 4

#define checkImageHeight 32
#define checkImageWidth 32
class D_CGL
{
 public:
  D_CGL(void);
  virtual ~D_CGL(void);
  void init(T_CArg *pCArg);
  void finish(T_CArg *pCArg);

  virtual void updateList(void);

  bool stop(void);
 private:
  void updateListLine(D_HeaderLine_t *pHeaderLine,D_DataLine_t *pArrayLine);
 private:
  void updateListSphere(D_Header_t *pHeader,D_Data_t *pArray,int m_ColorMode);
  void updateListTexture(D_Header_t *pHeader,D_Data_t *pArray,int m_ColorMode);
  void updateListCube(D_Header_t *pHeader,D_Data_t *pArray,int m_ColorMode);
  void updateListCircle(D_Header_t *pHeader,D_Data_t *pArray,int m_ColorMode);
  void updateListBlip(D_Header_t *pHeader,D_Data_t *pArray,int m_ColorMode);
  
  void updateListSphereVelocity(D_Header_t *pHeader,D_Data_t *pArray);
  void updateListSphereVelocityFixed(D_Header_t *pHeader,D_Data_t *pArray);
  void updateListSphereAngularVelocity(D_Header_t *pHeader,D_Data_t *pArray);
  void updateListSphereVelocityNorm(D_Header_t *pHeader,D_Data_t *pArray);
  void updateListSphereVelocityNormFixed(D_Header_t *pHeader,D_Data_t *pArray);
  void updateListSphereType(D_Header_t *pHeader,D_Data_t *pArray);
  void updateListSphereId(D_Header_t *pHeader,D_Data_t *pArray);
  void updateListSphereColor(D_Header_t *pHeader,D_Data_t *pArray);
  void updateListSphereForeground(D_Header_t *pHeader,D_Data_t *pArray);

  void updateListTextureVelocity(D_Header_t *pHeader,D_Data_t *pArray);
  void updateListTextureVelocityFixed(D_Header_t *pHeader,D_Data_t *pArray);
  void updateListTextureAngularVelocity(D_Header_t *pHeader,D_Data_t *pArray);
  void updateListTextureVelocityNorm(D_Header_t *pHeader,D_Data_t *pArray);
  void updateListTextureVelocityNormFixed(D_Header_t *pHeader,D_Data_t *pArray);
  void updateListTextureType(D_Header_t *pHeader,D_Data_t *pArray);
  void updateListTextureId(D_Header_t *pHeader,D_Data_t *pArray);
  void updateListTextureColor(D_Header_t *pHeader,D_Data_t *pArray);
  void updateListTextureForeground(D_Header_t *pHeader,D_Data_t *pArray);

  void updateListCubeVelocity(D_Header_t *pHeader,D_Data_t *pArray);
  void updateListCubeVelocityFixed(D_Header_t *pHeader,D_Data_t *pArray);
  void updateListCubeAngularVelocity(D_Header_t *pHeader,D_Data_t *pArray);
  void updateListCubeVelocityNorm(D_Header_t *pHeader,D_Data_t *pArray);
  void updateListCubeVelocityNormFixed(D_Header_t *pHeader,D_Data_t *pArray);
  void updateListCubeType(D_Header_t *pHeader,D_Data_t *pArray);
  void updateListCubeId(D_Header_t *pHeader,D_Data_t *pArray);
  void updateListCubeColor(D_Header_t *pHeader,D_Data_t *pArray);
  void updateListCubeForeground(D_Header_t *pHeader,D_Data_t *pArray);

  void updateListCircleVelocity(D_Header_t *pHeader,D_Data_t *pArray);
  void updateListCircleVelocityFixed(D_Header_t *pHeader,D_Data_t *pArray);
  void updateListCircleAngularVelocity(D_Header_t *pHeader,D_Data_t *pArray);
  void updateListCircleVelocityNorm(D_Header_t *pHeader,D_Data_t *pArray);
  void updateListCircleVelocityNormFixed(D_Header_t *pHeader,D_Data_t *pArray);
  void updateListCircleType(D_Header_t *pHeader,D_Data_t *pArray);
  void updateListCircleId(D_Header_t *pHeader,D_Data_t *pArray);
  void updateListCircleColor(D_Header_t *pHeader,D_Data_t *pArray);
  void updateListCircleForeground(D_Header_t *pHeader,D_Data_t *pArray);

  void updateListBlipVelocity(D_Header_t *pHeader,D_Data_t *pArray);
  void updateListBlipVelocityFixed(D_Header_t *pHeader,D_Data_t *pArray);
  void updateListBlipAngularVelocity(D_Header_t *pHeader,D_Data_t *pArray);
  void updateListBlipVelocityNorm(D_Header_t *pHeader,D_Data_t *pArray);
  void updateListBlipVelocityNormFixed(D_Header_t *pHeader,D_Data_t *pArray);
  void updateListBlipType(D_Header_t *pHeader,D_Data_t *pArray);
  void updateListBlipId(D_Header_t *pHeader,D_Data_t *pArray);
  void updateListBlipColor(D_Header_t *pHeader,D_Data_t *pArray);
  void updateListBlipForeground(D_Header_t *pHeader,D_Data_t *pArray);

 private:
  void initBorders(int Number);
 public:

  void initializeGL();
  void paintGL();
  void resizeGL(int w,int h);

  void reportId(char c,T_Int32_t id);

 public:
  // bestimmt den Modus in welcher Form die Partikel dargestellt werden sollen
  int m_DisplayMode;
  // bestimmt den Modus welche Bedeutung die Farben auf den Partikeln haben
  int m_ColorMode;

 public:
  // OpenGL Liste aller Partikel
  GLuint m_glList;
  // OpenGL Liste aller Linien
  GLuint m_glListLine;
  // OpenGL Liste fuer feststehende Waende etc
  GLuint m_glListFixed[NUMBER_OF_BORDERS];
  //
  char *m_glListFixedFilename[NUMBER_OF_BORDERS];

  GLuint m_glListParticle;

 public:
  GLfloat m_HeadOffsetX;
  GLfloat m_HeadOffsetY;
  GLfloat m_HeadOffsetZ;

 public:
  GLfloat m_ratio;
  GLfloat m_xRot;
  GLfloat m_yRot;
  GLfloat m_zRot;
  GLfloat m_scale;
  GLfloat m_viewscale;
  GLfloat m_offset;
  GLfloat m_offsetx;
  GLfloat m_offsety;
  GLfloat m_radius;
  bool m_bPerspective;
  GLfloat m_VelocityFixed;
  GLfloat m_NormX,m_NormY,m_NormZ;
  GLfloat m_NormFixedScale;
  bool m_bLight;
  GLfloat m_lightdistance;
  bool m_bTypeRed;
  bool m_bTypeGreen;
  bool m_bTypeBlue;
  T_Int32_t m_TypeRed;
  T_Int32_t m_TypeGreen;
  T_Int32_t m_TypeBlue;
  bool m_bIdRed;
  bool m_bIdGreen;
  bool m_bIdBlue;
  T_Int32_t m_IdRed;
  T_Int32_t m_IdGreen;
  T_Int32_t m_IdBlue;
  bool m_bColor;
  bool m_bglListLine;
  bool m_bglListFixed[NUMBER_OF_BORDERS];
  int m_pictureheight;
  D_Data_t *m_pArray;
  D_Header_t *m_pHeader;
  D_DataLine_t *m_pArrayLine;
  D_HeaderLine_t *m_pHeaderLine;
  bool m_bStop;
  char *m_pForeground;
  GLfloat m_paForeground[4];
  char *m_pBackground;

 private:
  struct {
   GLubyte checkImage[checkImageHeight][checkImageWidth][4];
   GLubyte otherImage[checkImageHeight][checkImageWidth][4];
  } m_texture;
};

#endif // D_CGL_H
