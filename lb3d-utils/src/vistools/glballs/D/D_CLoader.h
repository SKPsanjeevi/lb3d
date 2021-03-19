
#if !defined(D_CLOADER_H)
#define D_CLOADER_H

#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdio.h>

#include "T_CArg.h"
#include "D_CData.h"

#if defined(WITH_QT)
#include <qstring.h>
#include <qtimer.h>
#include <qstatusbar.h>

#include "D_CGLBox.h"
class D_CGLBox;
#else
#include "D_CGLX.h"
class D_CGLX;
#endif

#if !defined(WITH_QT)
#include "fakeqt.h"
#endif

class D_CLoader : public QObject
{
 Q_OBJECT
 public:
  D_CLoader(void);
  ~D_CLoader(void);

  void initDataFormat(char *pFormat);
  void initHeaderFormat(char *pFormat);

 public slots:
  void load(void);
  void rewind(void);

#if defined(WITH_QT)
 public:
  bool init(QStatusBar *pStatusBar,D_CGLBox *pBox,T_CArg *pCArg);
 private:
  QTimer *m_pTimer;
  int m_timestep;
  D_CGLBox *m_pBox;
  QStatusBar *m_pStatusBar;
  bool m_bPictures;
#else 
 public:
  bool init(D_CGLX *pBox,T_CArg *pCArg);
 private:
  D_CGLX *m_pGLX;
#endif

 private:
  FILE *m_file;
  T_N_t m_DataFormatPosMax;
  T_N_t *m_pDataFormat;
  T_N_t m_HeaderFormatPosMax;
  T_N_t *m_pHeaderFormat;
  D_Data_t *m_pArray;
  T_Int32_t m_ArrayLength;
  D_Header_t *m_pHeader;
  D_DataLine_t *m_pArrayLine;
  T_Int32_t m_ArrayLineLength;
  D_HeaderLine_t *m_pHeaderLine;
  int m_format;
  int m_headerformat;
  bool m_btMin,m_btMax;
  double m_tMin,m_tMax;
  T_Int32_t m_Step;
  T_Int32_t m_it;
};

#endif
