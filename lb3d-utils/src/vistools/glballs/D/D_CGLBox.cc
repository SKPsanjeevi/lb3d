#include <math.h>
#include <stdio.h>
#include <GL/glu.h>

#include "D_CGLBox.h"

/*!
  Create a D_CGLBox widget
*/
D_CGLBox::D_CGLBox( QWidget* parent, const char* name, const QGLWidget* shareWidget )
    : QGLWidget( parent, name, shareWidget )
{
 m_filenameCounter=0l;
 m_bStop=FALSE;
 m_keyPressEventColor=0;
 setFocusPolicy(QWidget::StrongFocus);
}


/*!
  Create a D_CGLBox widget
*/
D_CGLBox::D_CGLBox( const QGLFormat& format, QWidget* parent, const char* name, 
	      const QGLWidget* shareWidget )
    : QGLWidget( format, parent, name, shareWidget )
{
 m_filenameCounter=0l;
 m_bStop=FALSE;
}

void D_CGLBox::init(T_CArg *pCArg)
{
 m_CGL.init(pCArg);

 char *filename;
 pCArg->getParameter("output","file",&filename);
 m_filename=filename;
 char *value;
 pCArg->getParameter("output","format",&value);
 m_format=QString(value).upper();
 pCArg->getParameter("output","height",&m_pictureheight);

// setSizePolicy(QSizePolicy(QSizePolicy::Preferred,QSizePolicy::Preferred,TRUE));
 setSizePolicy(QSizePolicy(QSizePolicy::Preferred,QSizePolicy::Minimum,TRUE));
 updateGeometry();
}

void D_CGLBox::finish(T_CArg *pCArg)
{
 m_CGL.finish(pCArg);
}


/*!
  Release allocated resources
*/
D_CGLBox::~D_CGLBox()
{
}

void D_CGLBox::setInterface(D_Header_t *pHeader,D_Data_t *pArray,D_HeaderLine_t *pHeaderLine,D_DataLine_t *pArrayLine)
{
 m_CGL.m_pHeader=pHeader;
 m_CGL.m_pArray=pArray;
 m_CGL.m_pHeaderLine=pHeaderLine;
 m_CGL.m_pArrayLine=pArrayLine;
}

void D_CGLBox::updateList(void)
{
 m_CGL.updateList();
}


/*!
  Set up the OpenGL rendering state, and define display list
*/

void D_CGLBox::initializeGL()
{
 m_CGL.initializeGL();
}


/*!
  Set up the OpenGL view port, matrix mode, etc.
*/

void D_CGLBox::resizeGL( int w, int h )
{
 m_CGL.resizeGL(w,h);
}


/*!
  Paint the box. The actual openGL commands for drawing the box are
  performed here.
*/

void D_CGLBox::paintGL()
{
 m_CGL.paintGL();
}

/*!
  Set the rotation angle of the object to \e degrees around the X axis.
*/

void D_CGLBox::setXRotation( int degrees )
{
    m_CGL.m_xRot = (GLfloat)(degrees % 360);
    updateGL();
}


/*!
  Set the rotation angle of the object to \e degrees around the Y axis.
*/

void D_CGLBox::setYRotation( int degrees )
{
    m_CGL.m_yRot = (GLfloat)(degrees % 360);
    updateGL();
}


/*!
  Set the rotation angle of the object to \e degrees around the Z axis.
*/
void D_CGLBox::setZRotation( int degrees )
{
    m_CGL.m_zRot = (GLfloat)(degrees % 360);
    updateGL();
}

/*!
  Set the scale of the object 
*/
void D_CGLBox::setScale(int scale)
{
    m_CGL.m_scale = (GLfloat)0.2*std::exp(0.03*(float)(scale));
    updateGL();
}

/*!
  Set the viewscale of the object 
*/
void D_CGLBox::setViewScale(int scale)
{
    m_CGL.m_viewscale = (GLfloat)0.2*std::exp(0.03*(float)(scale));
    resizeGL(width(),height());
    updateGL();
}


/*!
  Set the offset of the object 
*/
void D_CGLBox::setOffset(int offset)
{
    m_CGL.m_offset = 0.1*(GLfloat)(offset);
    updateGL();
}

/*!
  Set the offset of the object 
*/
void D_CGLBox::setOffsetX(int x)
{
    m_CGL.m_offsetx = 0.1*(GLfloat)(x);
    updateGL();
}

/*!
  Set the offset of the object 
*/
void D_CGLBox::setOffsetY(int y)
{
    m_CGL.m_offsety = 0.1*(GLfloat)(y);
    updateGL();
}

/*
void D_CGLBox::mousePressEvent(QMouseEvent *e)
{
 m_dragPoint=e->pos();
}
void D_CGLBox::mouseReleaseEvent(QMouseEvent *e)
{
 if(e->button() & RightButton){
  m_CGL.m_offsetx=0.;
  m_CGL.m_offsety=0.;
 } else {
  m_CGL.m_offsetx+=10*(e->pos().x()-m_dragPoint.x())/(float)width();
  m_CGL.m_offsety-=10*(e->pos().y()-m_dragPoint.y())/(float)height();
 }
 updateGL();
}
*/

void D_CGLBox::slotStop()
{
 if(m_bStop)m_bStop=FALSE;
 else m_bStop=TRUE;
}

#include "qimage.h"
#include "qfiledialog.h"

void D_CGLBox::slotSave()
{
 QString filters;
 QString filename;

 m_bStop=TRUE;

 for(unsigned int i=0;i<QImageIO::outputFormats().count();i++){
  QString str=QString(QImageIO::outputFormats().at(i)).lower();
  filters+=QString("*.%1 ").arg(str);
 }
 filename=QFileDialog::getSaveFileName(QString::null,filters,this);
 if(filename.isEmpty()){
  m_bStop=FALSE;
  return;
 }
 for(unsigned int i=0;i<QImageIO::outputFormats().count();i++){
  QRegExp reg=QRegExp(QString(QImageIO::outputFormats().at(i)),FALSE,FALSE);
  if(reg.match(filename)!=-1){
   QPixmap Pixmap;
   Pixmap=renderPixmap((int)(m_CGL.m_ratio*m_pictureheight),m_pictureheight,TRUE);
   if(!Pixmap.isNull()){
    Pixmap.save(filename,QString(QImageIO::outputFormats().at(i)).latin1());
   }
   m_bStop=FALSE;
   return;
  }
 }
 m_bStop=FALSE;
}

void D_CGLBox::savePicture(void)
{
 //double v_width = width();
 //double v_height = height();

 // Make GL Context current
 makeCurrent();

 // Copy from OpenGL
 QImage *tempImage = new QImage( grabFrameBuffer() );

 if(!tempImage)return;

 QString filename=m_filename+QString().sprintf(".%.6i.",m_filenameCounter++)+m_format.lower();
 if(!tempImage->save(filename,m_format.latin1()))return;

 // Cleaning memory
 delete tempImage;

/*
 QPixmap Pixmap;
 makeCurrent();
 Pixmap=renderPixmap((int)(m_CGL.m_ratio*m_pictureheight),m_pictureheight,FALSE);
 if(!Pixmap.isNull()){
  QString filename=m_filename+QString().sprintf(".%.6i.",m_filenameCounter++)+m_format.lower();
  Pixmap.save(filename,m_format.latin1());
 }
*/
}

int D_CGLBox::heightForWidth(int w) const
{
 return (int)(w/m_CGL.m_ratio);
}

bool D_CGLBox::stop(void)
{
 return m_bStop;
}

void D_CGLBox::setDisplayMode(int mode){
 m_CGL.m_DisplayMode=mode;
 updateList();
 updateGL();
}

void D_CGLBox::setColorMode(int mode){
 m_CGL.m_ColorMode=mode;
 updateList();
 updateGL();
}

void D_CGLBox::keyPressEvent(QKeyEvent *e){
 if(!e->state()==Qt::NoButton)return;

 if(e->key()==Qt::Key_0 || e->key()==Qt::Key_1 || e->key()==Qt::Key_2 || e->key()==Qt::Key_3 || e->key()==Qt::Key_4 || e->key()==Qt::Key_5 || e->key()==Qt::Key_6 || e->key()==Qt::Key_7 || e->key()==Qt::Key_8 || e->key()==Qt::Key_9){
  m_keyPressEventNumber=m_keyPressEventNumber+e->text();
  return;
 }
 if(e->key()==Qt::Key_Right){
  m_keyPressEventColor=(m_keyPressEventColor+1)%3;
 }
 if(e->key()==Qt::Key_Left){
  m_keyPressEventColor=(m_keyPressEventColor+2)%3;
 }
 if(e->key()==Qt::Key_Up){
  if(m_CGL.m_ColorMode==COLORMODE_TYPE){
   switch(m_keyPressEventColor){
    case 0: m_CGL.m_bTypeRed=TRUE; m_CGL.m_TypeRed++; break;
    case 1: m_CGL.m_bTypeGreen=TRUE; m_CGL.m_TypeGreen++; break;
    case 2: m_CGL.m_bTypeBlue=TRUE; m_CGL.m_TypeBlue++; break;
   }
   updateList();
   updateGL();
  }
  if(m_CGL.m_ColorMode==COLORMODE_ID){
   switch(m_keyPressEventColor){
    case 0: m_CGL.m_bIdRed=TRUE; m_CGL.m_IdRed++; break;
    case 1: m_CGL.m_bIdGreen=TRUE; m_CGL.m_IdGreen++; break;
    case 2: m_CGL.m_bIdBlue=TRUE; m_CGL.m_IdBlue++; break;
   }
   updateList();
   updateGL();
  }
 }
 if(e->key()==Qt::Key_Down){
  if(m_CGL.m_ColorMode==COLORMODE_TYPE){
   switch(m_keyPressEventColor){
    case 0: m_CGL.m_TypeRed--; if(m_CGL.m_TypeRed<0)m_CGL.m_bTypeRed=FALSE; break;
    case 1: m_CGL.m_TypeGreen--; if(m_CGL.m_TypeGreen<0)m_CGL.m_bTypeGreen=FALSE; break;
    case 2: m_CGL.m_TypeBlue--; if(m_CGL.m_TypeBlue<0)m_CGL.m_bTypeBlue=FALSE; break;
   }
   updateList();
   updateGL();
  }
  if(m_CGL.m_ColorMode==COLORMODE_ID){
   switch(m_keyPressEventColor){
    case 0: m_CGL.m_IdRed--; if(m_CGL.m_IdRed<0)m_CGL.m_bIdRed=FALSE; break;
    case 1: m_CGL.m_IdGreen--; if(m_CGL.m_IdGreen<0)m_CGL.m_bIdGreen=FALSE; break;
    case 2: m_CGL.m_IdBlue--; if(m_CGL.m_IdBlue<0)m_CGL.m_bIdBlue=FALSE; break;
   }
   updateList();
   updateGL();
  }
 }
 if(e->key()==Qt::Key_Space){
  if(m_CGL.m_ColorMode==COLORMODE_ID){
   switch(m_keyPressEventColor){
    case 0: if(m_CGL.m_bIdRed)m_CGL.reportId('R',m_CGL.m_IdRed); break;
    case 1: if(m_CGL.m_bIdGreen)m_CGL.reportId('G',m_CGL.m_IdGreen); break;
    case 2: if(m_CGL.m_bIdBlue)m_CGL.reportId('B',m_CGL.m_IdBlue); break;
   }
  }
 }

 if(m_keyPressEventNumber.isEmpty())return;
 if(e->key()==Qt::Key_R){
  if(m_CGL.m_ColorMode==COLORMODE_TYPE){
   m_keyPressEventColor=0;
   m_CGL.m_bTypeRed=TRUE;
   m_CGL.m_TypeRed=m_keyPressEventNumber.toInt();  
   updateList();
   updateGL();
  }
  if(m_CGL.m_ColorMode==COLORMODE_ID){
   m_keyPressEventColor=0;
   m_CGL.m_bIdRed=TRUE;
   m_CGL.m_IdRed=m_keyPressEventNumber.toInt();  
   updateList();
   updateGL();
  }
 }
 if(e->key()==Qt::Key_G){
  if(m_CGL.m_ColorMode==COLORMODE_TYPE){
   m_keyPressEventColor=1;
   m_CGL.m_bTypeGreen=TRUE;
   m_CGL.m_TypeGreen=m_keyPressEventNumber.toInt();  
   updateList();
   updateGL();
  }
  if(m_CGL.m_ColorMode==COLORMODE_ID){
   m_keyPressEventColor=1;
   m_CGL.m_bIdGreen=TRUE;
   m_CGL.m_IdGreen=m_keyPressEventNumber.toInt();  
   updateList();
   updateGL();
  }
 }
 if(e->key()==Qt::Key_B){
  if(m_CGL.m_ColorMode==COLORMODE_TYPE){
   m_keyPressEventColor=2;
   m_CGL.m_bTypeBlue=TRUE;
   m_CGL.m_TypeBlue=m_keyPressEventNumber.toInt();  
   updateList();
   updateGL();
  }
  if(m_CGL.m_ColorMode==COLORMODE_ID){
   m_keyPressEventColor=2;
   m_CGL.m_bIdBlue=TRUE;
   m_CGL.m_IdBlue=m_keyPressEventNumber.toInt();  
   updateList();
   updateGL();
  }
 }
 m_keyPressEventNumber=QString();
}
