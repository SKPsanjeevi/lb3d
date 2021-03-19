#ifndef D_CGLBox_H
#define D_CGLBox_H

#include <math.h>
#include <qgl.h>
#include <qpoint.h>
#include <qpixmap.h>
#include <qregexp.h>

#include "T_CArg.h"
#include "D_CData.h"
#include "D_CLoader.h"
#include "D_CGL.h"

class D_CLoader;
class D_CGLBox : public QGLWidget
{
 Q_OBJECT
 public:
  D_CGLBox(QWidget *parent,const char *name,const QGLWidget *shareWidget=0);
  D_CGLBox(const QGLFormat &format,QWidget *parent,const char *name, 
	   const QGLWidget *shareWidget=0);
  virtual ~D_CGLBox();
  void init(T_CArg *pCArg);
  void finish(T_CArg *pCArg);

  void setInterface(D_Header_t *pHeader,D_Data_t *pArray,D_HeaderLine_t *pHeaderLine,D_DataLine_t *pArrayLine);
  void updateList(void);
  void savePicture(void);

  bool stop(void);
 private:
  D_CGL m_CGL;

 public:
  inline int getDisplayMode(void) const;
  inline int getColorMode(void) const;
  inline int getXRotation(void) const;
  inline int getYRotation(void) const;
  inline int getZRotation(void) const;
  inline int getScale(void) const;
  inline int getViewScale(void) const;
  inline int getOffset(void) const;
//  inline int getOffsetX(void) const;
//  inline int getOffsetY(void) const;

 public:
  void setDisplayMode(int mode);
  void setColorMode(int mode);
 public slots:
  void setXRotation(int degrees);
  void setYRotation(int degrees);
  void setZRotation(int degrees);
  void setScale(int scale);
  void setViewScale(int scale);
  void setOffset(int offset);
  void setOffsetX(int offset);
  void setOffsetY(int offset);
  void slotSave(void);
  void slotStop(void);


 protected:
  void keyPressEvent(QKeyEvent *e);

  void initializeGL();
  void paintGL();
  void resizeGL(int w,int h);
//  void mousePressEvent(QMouseEvent *e);
//  void mouseReleaseEvent(QMouseEvent *e);
  int heightForWidth(int w) const;
//  QSize sizeHint () const;

 private:
  QPoint m_dragPoint;
  QString m_filename;
  QString m_format;
  T_Int32_t m_filenameCounter;
  int m_pictureheight;
  bool m_bStop;
  QString m_keyPressEventNumber;
  int m_keyPressEventColor;
};

#include "D_CGLBox_include.h"

#endif // D_CGLBox_H
