#include "D_CDisplay.h"
#include <qapplication.h>
#include <qmainwindow.h>
#include <qgl.h>

/*
  The main program is here. 
*/

int main( int argc, char **argv )
{
 bool bNoError=TRUE;

 QApplication *pA;
 D_CDisplay *pDisplay;
 T_CArg *pCArg=NULL;
 D_CLoader *pCLoader=NULL;
 QMainWindow *pMainWindow=NULL;

 pCArg=new T_CArg("main");
 pCLoader=new D_CLoader();
 if(!pCArg || !pCLoader)bNoError=FALSE;
 
 if(bNoError){
  pCArg->registerParameter("main","file","f","-","filename");
  pCArg->registerParameter("main","format",NULL,"100","data format (7|8|100)");
  pCArg->registerParameter("main","header",NULL,"7","header format (7|8)");
  pCArg->registerParameter("main","timer","t","10","time waited between reading from file"); 
  pCArg->registerParameter("main","tmin",NULL,"false","plot for time >="); 
  pCArg->registerParameter("main","tmax",NULL,"false","plot for time <="); 
  pCArg->registerParameter("main","step",NULL,"1","plot every x steps"); 
  pCArg->registerParameter("display","width",NULL,"500","width of the main window");
  pCArg->registerParameter("display","height",NULL,"500","height of the main window");
  pCArg->registerParameter("display","ratio",NULL,"1.","ratio width : height");
  pCArg->registerParameter("display","xRot","x","0.","rotation arround x-axis");
  pCArg->registerParameter("display","yRot","y","0.","rotation arround y-axis");
  pCArg->registerParameter("display","zRot","z","0.","rotation arround z-axis");
  pCArg->registerParameter("display","scale","s","1.","scaling on the objects");
  pCArg->registerParameter("display","viewscale","vs","1.","scaling on the size of the viewing window");
  pCArg->registerParameter("display","offset","o","-3.","offset from the screen");
  pCArg->registerParameter("display","offsetx","ox","0.","offset from the screen");
  pCArg->registerParameter("display","offsety","oy","0.","offset from the screen");
  pCArg->registerParameter("display","radius","r","0.1","particles with smaller radiuses are painted as oktaeders");
  pCArg->registerParameter("display","perspective",NULL,"true","perspective viewing");
  pCArg->registerParameter("display","foreground","fg","xFFFFFF","foreground color");
  pCArg->registerParameter("display","background","bg","x000000","background color");
  pCArg->registerParameter("display","mode","m","sphere","controls the shape of particles");
  pCArg->registerParameter("display","color","c","velocity","use velocity as color");
  pCArg->registerParameter("display","velocityfixed","vyfd","1","max velocity for fixed displaymode");
  pCArg->registerParameter("display","normx","nmx","1","direction");
  pCArg->registerParameter("display","normy","nmy","0","direction");
  pCArg->registerParameter("display","normz","nmz","0","direction");
  pCArg->registerParameter("display","normfixedscale","nmfs","1","scale for direction");
  pCArg->registerParameter("display","typeRed","tr","false","draw paricles of this type in red");
  pCArg->registerParameter("display","typeGreen","tg","false","draw paricles of this type in green");
  pCArg->registerParameter("display","typeBlue","tb","false","draw paricles of this type in blue");
  pCArg->registerParameter("display","idRed","ir","false","draw paricles of this id in red");
  pCArg->registerParameter("display","idGreen","ig","false","draw paricles of this id in green");
  pCArg->registerParameter("display","idBlue","ib","false","draw paricles of this id in blue");
  pCArg->registerParameter("display","light",NULL,"true","with spot light");
//  pCArg->registerParameter("display","lightdistance",NULL,"1.","distance of the spot light");
  pCArg->registerParameter("display","lines",NULL,"false","draw lines");

  for(int i=0;i<NUMBER_OF_BORDERS;i++){
   char s[10];
   sprintf(s,"border%i",i+1);
   pCArg->registerParameter("display",s,NULL,"false","with borders from this file");
  }
  pCArg->registerParameter("output","flag",NULL,"false","output display as picture");
  pCArg->registerParameter("output","file",NULL,"output","picture basename");
  pCArg->registerParameter("output","format",NULL,"ppm","picture format");
  pCArg->registerParameter("output","height",NULL,"100","picture height");
 }
 if(bNoError){
  bNoError&=pCArg->parseArgs(argc,argv);
 }

 if(bNoError){
  QApplication::setColorSpec( QApplication::CustomColor );
  pA=new QApplication(argc,argv);			

  if ( !QGLFormat::hasOpenGL() ) {
   qWarning( "This system has no OpenGL support. Exiting." );
   return -1;
  }

  pMainWindow=new QMainWindow();
  pDisplay=new D_CDisplay(pMainWindow,pCArg,pCLoader,pMainWindow);
  if(!pDisplay||(!pMainWindow))bNoError=FALSE;
 }
 if(bNoError){
  pMainWindow->setCentralWidget(pDisplay);
  bNoError&=pCLoader->init(pMainWindow->statusBar(),pDisplay->getGLBox(),pCArg);
 }
 if(bNoError){
  int width,height;
  pCArg->getParameter("display","width",&width);
  pCArg->getParameter("display","height",&height);
  pMainWindow->resize(width,height);
  pA->setMainWidget(pMainWindow);
  pMainWindow->show();
  bNoError=pA->exec();
  (pDisplay->getGLBox())->finish(pCArg);
 }

 T_DELETE(pCArg);
 return bNoError;
}
