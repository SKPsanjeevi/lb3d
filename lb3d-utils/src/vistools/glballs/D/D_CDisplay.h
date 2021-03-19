#if !defined(D_CDISPLAY_H)
#define D_CDISPLAY_H

#include <qwidget.h>
#include <qmainwindow.h>
#include "D_CLoader.h"
#include <math.h>

#include <qpushbutton.h>
#include <qslider.h>
#include <qlayout.h>
#include <qframe.h>
#include <qlabel.h>
#include <qmenubar.h>
#include <qpopupmenu.h>
#include <qapplication.h>
#include <qkeycode.h>
#include <qpixmap.h>
#include <qpainter.h>
#include "D_CGLBox.h"


class D_CGLBox;

class D_CDisplay : public QWidget
{
  Q_OBJECT
 public:
  D_CDisplay(QMainWindow *pMain,T_CArg *pCArg,D_CLoader *pCLoader,QWidget* parent = 0, const char* name = 0);
  void finish(T_CArg *pCArg);
  D_CGLBox *getGLBox(void){return m_pCGLBox;}
 private:
  D_CGLBox *m_pCGLBox;
  QMainWindow *m_pMain;

 public slots:
  void slotDisplayModeSphere(void);
  void slotDisplayModeTexture(void);
  void slotDisplayModeCube(void);
  void slotDisplayModeCircle(void);
  void slotDisplayModeBlip(void);
  void slotColorModeVelocity(void);
  void slotColorModeVelocityFixed(void);
  void slotColorModeAngularVelocity(void);
  void slotColorModeVelocityNorm(void);
  void slotColorModeVelocityNormFixed(void);
  void slotColorModeType(void);
  void slotColorModeId(void);
  void slotColorModeColor(void);
  void slotColorModeForeground(void);
 private :
  int OptionIDDisplayModeSphere;
  int OptionIDDisplayModeTexture;
  int OptionIDDisplayModeCube;
  int OptionIDDisplayModeCircle;
  int OptionIDDisplayModeBlip;
  int OptionIDColorModeVelocity;
  int OptionIDColorModeVelocityFixed;
  int OptionIDColorModeAngularVelocity;
  int OptionIDColorModeVelocityNorm;
  int OptionIDColorModeVelocityNormFixed;
  int OptionIDColorModeType;
  int OptionIDColorModeId;
  int OptionIDColorModeColor;
  int OptionIDColorModeForeground;
 bool m_bPictures;
};

#endif 
