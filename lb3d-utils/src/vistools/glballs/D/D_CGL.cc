#include <math.h>
#include <stdio.h>
#include <GL/glu.h>

#include "D_CGL.h"

/*
*/
D_CGL::D_CGL(void)
{
 m_ratio=1.;
 m_pHeader=NULL;
 m_pArray=NULL;

 m_pHeaderLine=NULL;
 m_pArrayLine=NULL;
 m_glListLine=0;
 m_bglListLine=FALSE;

 m_glList=0;
 for(int i=0;i<NUMBER_OF_BORDERS;i++){
  m_glListFixed[i]=0;
  m_glListFixedFilename[i]=NULL;
  m_bglListFixed[i]=FALSE;
 }

 m_bStop=FALSE;
 m_NormFixedScale=1.;
 m_DisplayMode=DISPLAYMODE_SPHERE;
 m_ColorMode=COLORMODE_VELOCITY;

 m_HeadOffsetX=0.;
 m_HeadOffsetY=0.;
 m_HeadOffsetZ=0.;

 m_pForeground=NULL;
 m_pBackground=NULL;
}

/*
*/
D_CGL::~D_CGL(void)
{
}

/*
*/
void D_CGL::init(T_CArg *pCArg)
{
 // default object rotation
 pCArg->getParameter("display","ratio",&m_ratio);
 pCArg->getParameter("display","xRot",&m_xRot);
 pCArg->getParameter("display","yRot",&m_yRot);
 pCArg->getParameter("display","zRot",&m_zRot);
 while(m_xRot<0.)m_xRot+=360;
 while(m_yRot<0.)m_yRot+=360;
 while(m_zRot<0.)m_zRot+=360;
 while(m_xRot>360.)m_xRot-=360;
 while(m_yRot>360.)m_yRot-=360;
 while(m_zRot>360.)m_zRot-=360;
 // default object scale
 pCArg->getParameter("display","scale",&m_scale);			
 pCArg->getParameter("display","viewscale",&m_viewscale);			
 pCArg->getParameter("display","offset",&m_offset);
 pCArg->getParameter("display","offsetx",&m_offsetx);
 pCArg->getParameter("display","offsety",&m_offsety);
 pCArg->getParameter("display","radius",&m_radius);
 pCArg->getParameter("display","perspective",&m_bPerspective);

 char *pforeground;
 pCArg->getParameter("display","foreground",&pforeground);
 m_pForeground=strdup(pforeground);

 // Vordergrund sei Weiss
 {
  int r=0,g=0,b=0,a=1; 
  sscanf(m_pForeground,"x%2x%2x%2x%2x",&r,&g,&b,&a);
  m_paForeground[0]=(r%256)/256.; 
  m_paForeground[1]=(g%256)/256.; 
  m_paForeground[2]=(b%256)/256.; 
  m_paForeground[3]=(a%256)/256.; 
 }

 char *pbackground;
 pCArg->getParameter("display","background",&pbackground);
 m_pBackground=strdup(pbackground);

 // Waehlt den Darstellungsmode
 char *pMode;
 pCArg->getParameter("display","mode",&pMode);
 m_DisplayMode=DISPLAYMODE_SPHERE;
 if(!strcmp("sphere",pMode))m_DisplayMode=DISPLAYMODE_SPHERE;
 if(!strcmp("texture",pMode))m_DisplayMode=DISPLAYMODE_TEXTURE;
 if(!strcmp("cube",pMode))m_DisplayMode=DISPLAYMODE_CUBE;
 if(!strcmp("circle",pMode))m_DisplayMode=DISPLAYMODE_CIRCLE;
 if(!strcmp("blip",pMode))m_DisplayMode=DISPLAYMODE_BLIP;

 // Waehlt den Farbmode
 pCArg->getParameter("display","color",&pMode);
 m_ColorMode=COLORMODE_VELOCITY;
 if(!strcmp("velocity",pMode))m_ColorMode=COLORMODE_VELOCITY;
 if(!strcmp("velocityfixed",pMode))m_ColorMode=COLORMODE_VELOCITYFIXED;
 if(!strcmp("angularvelocity",pMode))m_ColorMode=COLORMODE_ANGULARVELOCITY;
 if(!strcmp("velocitynorm",pMode))m_ColorMode=COLORMODE_VELOCITYNORM;
 if(!strcmp("velocitynormfixed",pMode))m_ColorMode=COLORMODE_VELOCITYNORMFIXED;
 if(!strcmp("type",pMode))m_ColorMode=COLORMODE_TYPE;
 if(!strcmp("id",pMode))m_ColorMode=COLORMODE_ID;
 if(!strcmp("color",pMode))m_ColorMode=COLORMODE_COLOR;
 if(!strcmp("fg",pMode))m_ColorMode=COLORMODE_FOREGROUND;
 if(!strcmp("foreground",pMode))m_ColorMode=COLORMODE_FOREGROUND;

 // Parameter fuer die verschiedenen Farbmodi
 pCArg->getParameter("display","velocityfixed",&m_VelocityFixed);
 pCArg->getParameter("display","normx",&m_NormX);
 pCArg->getParameter("display","normy",&m_NormY);
 pCArg->getParameter("display","normz",&m_NormZ);
 pCArg->getParameter("display","normfixedscale",&m_NormFixedScale);
 pCArg->getParameter("display","typeRed",&m_bTypeRed); 
 pCArg->getParameter("display","typeGreen",&m_bTypeGreen); 
 pCArg->getParameter("display","typeBlue",&m_bTypeBlue); 
 pCArg->getParameter("display","typeRed",&m_TypeRed); 
 pCArg->getParameter("display","typeGreen",&m_TypeGreen); 
 pCArg->getParameter("display","typeBlue",&m_TypeBlue); 
 pCArg->getParameter("display","idRed",&m_bIdRed); 
 pCArg->getParameter("display","idGreen",&m_bIdGreen); 
 pCArg->getParameter("display","idBlue",&m_bIdBlue); 
 pCArg->getParameter("display","idRed",&m_IdRed); 
 pCArg->getParameter("display","idGreen",&m_IdGreen); 
 pCArg->getParameter("display","idBlue",&m_IdBlue); 
 if(!m_bTypeRed)m_TypeRed=-1;
 if(!m_bTypeGreen)m_TypeGreen=-1;
 if(!m_bTypeBlue)m_TypeBlue=-1;
 if(!m_bIdRed)m_IdRed=-1;
 if(!m_bIdGreen)m_IdGreen=-1;
 if(!m_bIdBlue)m_IdBlue=-1;

 pCArg->getParameter("display","light",&m_bLight);
 m_lightdistance=1.; // pCArg->getDoubleParameter("display","lightdistance");

 {
  char s[10];
  for(int i=0;i<NUMBER_OF_BORDERS;i++){
   sprintf(s,"border%i",i+1);
   if(pCArg->isValidParameter("display",s)){
    pCArg->getParameter("display",s,&m_glListFixedFilename[i]);
    m_bglListFixed[i]=TRUE;
   }
  }
 }
 pCArg->getParameter("display","lines",&m_bglListLine);

}

void D_CGL::finish(T_CArg *pCArg)
{
 char s[100];
 sprintf(s,"%f",m_ratio);
 pCArg->setParameter("display","ratio",s);
 sprintf(s,"%f",m_xRot);
 pCArg->setParameter("display","xRot",s);
 sprintf(s,"%f",m_yRot);
 m_yRot=pCArg->setParameter("display","yRot",s);
 sprintf(s,"%f",m_zRot);
 m_zRot=pCArg->setParameter("display","zRot",s);
 sprintf(s,"%f",m_scale);
 m_scale=pCArg->setParameter("display","scale",s);			
 sprintf(s,"%f",m_viewscale);
 m_viewscale=pCArg->setParameter("display","viewscale",s);
 sprintf(s,"%f",m_offset);
 m_offset=pCArg->setParameter("display","offset",s);
 sprintf(s,"%f",m_offsetx);
 m_offsetx=pCArg->setParameter("display","offsetx",s);
 sprintf(s,"%f",m_offsety);
 m_offsety=pCArg->setParameter("display","offsety",s);
 sprintf(s,"%f",m_radius);
 m_radius=pCArg->setParameter("display","radius",s);
 if(m_bPerspective)pCArg->setParameter("display","perspective","true");
 else pCArg->setParameter("display","perspective","false");

 switch(m_DisplayMode){
  case DISPLAYMODE_SPHERE: pCArg->setParameter("display","mode","sphere"); break;
  case DISPLAYMODE_TEXTURE: pCArg->setParameter("display","mode","texture"); break;
  case DISPLAYMODE_CUBE: pCArg->setParameter("display","mode","cube"); break;
  case DISPLAYMODE_CIRCLE: pCArg->setParameter("display","mode","circle"); break;
  case DISPLAYMODE_BLIP: pCArg->setParameter("display","mode","blip"); break;
 }

 switch(m_ColorMode){
  case COLORMODE_VELOCITY: pCArg->setParameter("display","color","velocity"); break;
  case COLORMODE_VELOCITYFIXED: pCArg->setParameter("display","color","velocityfixed"); break;
  case COLORMODE_ANGULARVELOCITY: pCArg->setParameter("display","color","angularvelocity"); break;
  case COLORMODE_VELOCITYNORM: pCArg->setParameter("display","color","velocitynorm"); break;
  case COLORMODE_VELOCITYNORMFIXED: pCArg->setParameter("display","color","velocitynormfixed"); break;
  case COLORMODE_TYPE: pCArg->setParameter("display","color","type"); break;
  case COLORMODE_ID: pCArg->setParameter("display","color","id"); break;
  case COLORMODE_COLOR: pCArg->setParameter("display","color","color"); break;
  case COLORMODE_FOREGROUND: pCArg->setParameter("display","color","foreground"); break;
 }

 sprintf(s,"%f",m_NormX);
 m_NormX=pCArg->setParameter("display","normx",s);
 sprintf(s,"%f",m_NormY);
 m_NormY=pCArg->setParameter("display","normy",s);
 sprintf(s,"%f",m_NormZ);
 m_NormZ=pCArg->setParameter("display","normz",s);
 sprintf(s,"%f",m_NormFixedScale);
 pCArg->setParameter("display","normfixedscale",s);
 if(m_bTypeRed){
  sprintf(s,"%i",m_TypeRed);
  pCArg->setParameter("display","typeRed",s);
 } else pCArg->setParameter("display","typeRed","false");
 if(m_bTypeGreen){
  sprintf(s,"%i",m_TypeGreen);
  pCArg->setParameter("display","typeGreen",s);
 } else pCArg->setParameter("display","typeGreen","false");
 if(m_bTypeBlue){
  sprintf(s,"%i",m_TypeBlue);
  pCArg->setParameter("display","typeBlue",s);
 } else pCArg->setParameter("display","typeBlue","false");
 if(m_bIdRed){
  sprintf(s,"%i",m_IdRed);
  pCArg->setParameter("display","idRed",s);
 } else pCArg->setParameter("display","idRed","false");
 if(m_bIdGreen){
  sprintf(s,"%i",m_IdGreen);
  pCArg->setParameter("display","idGreen",s);
 } else pCArg->setParameter("display","idGreen","false");
 if(m_bIdBlue){
  sprintf(s,"%i",m_IdBlue);
  pCArg->setParameter("display","idBlue",s);
 } else pCArg->setParameter("display","idBlue","false");

 if(m_bLight)pCArg->setParameter("display","light","true");
 else pCArg->setParameter("display","light","false");
// sprintf(s,"%f",m_lightdistance);
// m_lightdistance=pCArg->setParameter("display","lightdistance",s);

 glDeleteLists(m_glList,1);

 for(int i=0;i<NUMBER_OF_BORDERS;i++){
  glDeleteLists(m_glListFixed[i],1);
 }
 if(m_bglListLine)pCArg->setParameter("display","lines","true");
 else pCArg->setParameter("display","line","false");
}


/*!
  Set up the OpenGL rendering state, and define display list
*/

void D_CGL::initializeGL()
{
 GLfloat light_position[4];
 GLfloat model_ambientlight[] = { 0.8, 0.8, 0.8, 1.0 };
 GLfloat model_ambient[] = { 1.0, 1.0, 1.0, 1.0 };
 GLfloat mat_zero[] = { 0., 0., 0., 0. };
 GLfloat mat_ambient[] = { 0.2, 0.2, 0.2, 0.8 };
 GLfloat mat_diffuse[] = { 0.8, 0.8, 0.8, 0.8 };

 GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };

 GLfloat mat_shininess[] = { 100. };
 // setzt Position der ersten Lichtquelle
 light_position[0] = m_lightdistance*m_scale; 
 light_position[1] = m_lightdistance*m_scale; 
 light_position[2] = m_lightdistance*m_scale; 
 light_position[3] = 0.;
 glLightfv(GL_LIGHT0, GL_POSITION, light_position);
 // setzt ungerichte Hintergundbeleuchtung
 if(m_bLight) glLightModelfv(GL_LIGHT_MODEL_AMBIENT, model_ambientlight);
 else  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, model_ambient);
 // Reflexion auf Oberfl"achen

  glMaterialfv (GL_FRONT, GL_AMBIENT, mat_ambient);
  glMaterialfv (GL_FRONT, GL_DIFFUSE, mat_diffuse);
  glMaterialfv (GL_FRONT, GL_SPECULAR, mat_specular);
  glMaterialfv (GL_FRONT, GL_SHININESS, mat_shininess);
  glMaterialfv (GL_BACK, GL_AMBIENT, mat_zero);
  glMaterialfv (GL_BACK, GL_DIFFUSE, mat_zero);
  glMaterialfv (GL_BACK, GL_SPECULAR, mat_zero);
  glMaterialfv (GL_BACK, GL_SHININESS, mat_shininess);


 // schaltet Beleuchtung generell an
 glEnable(GL_LIGHTING);
 // schaltet erste Lichtquelle an
 if(m_bLight) glEnable(GL_LIGHT0);
 // verdeckte Pixel werden ausgeblendet
 glEnable(GL_DEPTH_TEST);


//DGTMP
//glEnable (GL_BLEND);
//glBlendFunc (GL_SRC_ALPHA, GL_ONE);
//glDepthMask(GL_FALSE);

 // With smooth shading, the color at each vertex is treated individually. 
 glShadeModel (GL_SMOOTH);

 // Hintergrund sei Schwarz
 {
  int r=0,g=0,b=0;
  sscanf(m_pBackground,"x%2x%2x%2x",&r,&g,&b);
  r=r%256; g=g%256; b=b%256;
  glClearColor((double)r/256.,(double)g/256.,(double)b/256.,0.);
 }

 int i, j, c;
 for(i=0;i<checkImageHeight;i++){
  for(j=0;j<checkImageWidth;j++){
   c=(((i&0x8)==0)^((j&0x8)==0))*150+105;
   m_texture.checkImage[i][j][0] = (GLubyte) c;
   m_texture.checkImage[i][j][1] = (GLubyte) c;
   m_texture.checkImage[i][j][2] = (GLubyte) c;
   c = (((i&0x10)==0)^((j&0x10)==0))*255;
   m_texture.otherImage[i][j][0] = (GLubyte) c;
   m_texture.otherImage[i][j][1] = (GLubyte) 0;
   m_texture.otherImage[i][j][2] = (GLubyte) 0;
   m_texture.otherImage[i][j][3] = (GLubyte) 255;
  }
 }
 glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
 glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
 glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
 glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  
 glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, checkImageWidth, checkImageHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, m_texture.checkImage);
 glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_COLOR, GL_DECAL);


 // alloziert eine Verwaltungliste f"ur gl Operationen
 if(m_glList==0) m_glList=glGenLists(1);
 if(m_bglListLine && m_glListLine==0) m_glListLine=glGenLists(1);
 for(int i=0;i<NUMBER_OF_BORDERS;i++){
  if(m_glListFixed[i]==0 && m_glListFixedFilename[i]){
   m_glListFixed[i]=glGenLists(1);
   initBorders(i); 
  }
 }
 // m_glListParticle=glGenLists(1); //DGTMP
 updateList();
}


/*!
  Set up the OpenGL view port, matrix mode, etc.
*/

void D_CGL::resizeGL( int w, int h )
{
 glViewport( 0, 0, (GLint)w, (GLint)h );
 glMatrixMode( GL_PROJECTION );
 glLoadIdentity();
 if(m_bPerspective){
  gluPerspective(40.*m_viewscale,m_ratio,1.,600);
 } else {
  glOrtho(-m_viewscale,+m_viewscale,-m_ratio*m_viewscale,+m_ratio*m_viewscale,1.,1000.);
 }
}


/*!
  Paint the . The actual openGL commands for drawing the box are
  performed here.
*/

void D_CGL::paintGL()
{
/*
 GLfloat light_position[4];
 // setzt Position der ersten Lichtquelle
 light_position[0] = m_lightdistance*m_scale; 
 light_position[1] = m_lightdistance*m_scale; 
 light_position[2] = m_lightdistance*m_scale; 
 light_position[3] = 0.;
 glLightfv(GL_LIGHT0, GL_POSITION, light_position);
*/
 glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

 glMatrixMode( GL_MODELVIEW );
 glLoadIdentity();
 glTranslatef( m_offsetx, m_offsety, m_offset );
 glScalef( m_scale, m_scale, m_scale );

 glRotatef( m_xRot, 1.0, 0.0, 0.0 ); 
 glRotatef( m_yRot, 0.0, 1.0, 0.0 ); 
 glRotatef( m_zRot, 0.0, 0.0, 1.0 );

 if(m_bglListLine){
  glDisable(GL_LIGHTING);
   glCallList(m_glListLine);
  glEnable(GL_LIGHTING);
 }

 glCallList(m_glList);

 for(int i=0;i<NUMBER_OF_BORDERS;i++){
  if(m_bglListFixed[i]){
   glTranslatef(-m_HeadOffsetX,-m_HeadOffsetY,-m_HeadOffsetZ);
   glCallList(m_glListFixed[i]);
  }
 }

 glFlush();
}

/*!
  Generate an OpenGL display list for the object to be shown, i.e. the 
*/
void D_CGL::updateList()
{	
 if((!m_pHeader) || (!m_pArray))return;

 if(m_DisplayMode==DISPLAYMODE_TEXTURE){
  glEnable(GL_TEXTURE_2D);
 } else {
  glDisable(GL_TEXTURE_2D);
 }

 if(m_bglListLine){
  updateListLine(m_pHeaderLine,m_pArrayLine);
 }

 switch(m_DisplayMode){
  case DISPLAYMODE_SPHERE:
   updateListSphere(m_pHeader,m_pArray,m_ColorMode);
   break;
  case DISPLAYMODE_TEXTURE:
   updateListTexture(m_pHeader,m_pArray,m_ColorMode);
   break;
  case DISPLAYMODE_CUBE:
   updateListCube(m_pHeader,m_pArray,m_ColorMode);
   break;
  case DISPLAYMODE_CIRCLE:
   updateListCircle(m_pHeader,m_pArray,m_ColorMode);
   break;
  case DISPLAYMODE_BLIP:
   updateListBlip(m_pHeader,m_pArray,m_ColorMode);
   break;
 }
}

/*!
  search for particles with this Id and print their parameters to stdout
*/
void D_CGL::reportId(char c,T_Int32_t id)
{	
 if((!m_pHeader) || (!m_pArray))return;
 
 for(long int i=0;i<m_pHeader->Number;++i){
  if(m_pArray[i].id==id){
   fprintf(stdout,"%c %"T_R3_PRINTF" %"T_R3_PRINTF" %"T_R3_PRINTF" %"T_R3_PRINTF" %"T_R3_PRINTF" %"T_R3_PRINTF" %"T_R3_PRINTF" %"T_R3_PRINTF" %"T_R3_PRINTF" %"T_R3_PRINTF" %"T_R3_PRINTF" %"T_R3_PRINTF" %"T_R3_PRINTF" %"T_R3_PRINTF" %"T_ID_PRINTF" %"T_PARTICLECLASSEXTERN_PRINTF" %"T_Int32_PRINTF"\n",
  c,
  m_pArray[i].x,m_pArray[i].y,m_pArray[i].z,  
  m_pArray[i].vx,m_pArray[i].vy,m_pArray[i].vz,  
  m_pArray[i].r,
  m_pArray[i].qw,m_pArray[i].qx,m_pArray[i].qy,m_pArray[i].qz,  
  m_pArray[i].wx,m_pArray[i].wy,m_pArray[i].wz,  
  m_pArray[i].id, m_pArray[i].type, m_pArray[i].color
  );
  }
 }
} 

void D_CGL::updateListLine(D_HeaderLine_t *pHeaderLine,D_DataLine_t *pArrayLine)
{
 GLfloat OffsetX,OffsetY,OffsetZ;
 OffsetX=0.5*(m_pHeader->xMax+m_pHeader->xMin);
 OffsetY=0.5*(m_pHeader->yMax+m_pHeader->yMin);
 OffsetZ=0.5*(m_pHeader->zMax+m_pHeader->zMin);

 if(pHeaderLine->Number){
  glNewList(m_glListLine, GL_COMPILE );
  glBegin(GL_LINES);
  for(long int i=0;i<pHeaderLine->Number;++i){
   glColor3f(
    0.00390625*(pArrayLine[i].color & 255),
    0.00390625*((pArrayLine[i].color >> 8)& 255),
    0.00390625*((pArrayLine[i].color >> 16)& 255)
   );
   glVertex3f(pArrayLine[i].x1-OffsetX,pArrayLine[i].y1-OffsetY,pArrayLine[i].z1-OffsetZ);
   glVertex3f(pArrayLine[i].x2-OffsetX,pArrayLine[i].y2-OffsetY,pArrayLine[i].z2-OffsetZ);
  }
  glEnd();
  glEndList();
 } else {
  glNewList(m_glListLine, GL_COMPILE );
  glEndList();
 }
}
