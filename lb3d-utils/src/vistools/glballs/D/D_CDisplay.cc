#include "D_CDisplay.h"

D_CDisplay::D_CDisplay(QMainWindow *pMain,T_CArg *pCArg,D_CLoader *pCLoader,QWidget *parent, const char *name )
    : QWidget( parent, name )
{
 pCArg->getParameter("output","flag",&m_bPictures);
 m_bPictures=FALSE;
    
 // Create an OpenGL widget
 m_pCGLBox=new D_CGLBox( this, "D_CGLBox1");
 if(m_bPictures)m_pCGLBox->hide();
 m_pCGLBox->init(pCArg);

 m_pMain=pMain;

 // Create a menu
 QPopupMenu *file=new QPopupMenu();
 if(!m_bPictures){
  file->insertItem("&Stop",m_pCGLBox, SLOT(slotStop()),CTRL+Key_S );
  file->insertItem("Save Picture",m_pCGLBox, SLOT(slotSave()),CTRL+Key_P );
  file->insertItem("&Rewind",pCLoader, SLOT(rewind()),CTRL+Key_R );
 }
 file->insertItem("Exit",qApp,SLOT(quit()),CTRL+Key_Q);

 QPopupMenu *options;
// if(!m_bPictures){
  // Create a menu
  options=new QPopupMenu();
  options->setCaption("DisplayMode");
  OptionIDDisplayModeSphere=options->insertItem("&Sphere",this, SLOT(slotDisplayModeSphere()),ALT+Key_S );
  OptionIDDisplayModeTexture=options->insertItem("Te&xture",this, SLOT(slotDisplayModeTexture()),ALT+Key_X );
  OptionIDDisplayModeCube=options->insertItem("C&ube",this, SLOT(slotDisplayModeCube()),ALT+Key_U );
  OptionIDDisplayModeCircle=options->insertItem("Ci&rcle",this, SLOT(slotDisplayModeCircle()),ALT+Key_R );
  OptionIDDisplayModeBlip=options->insertItem("&Blip",this, SLOT(slotDisplayModeBlip()),ALT+Key_B );
  options->insertSeparator();
  options->setCaption("ColorMode");
  OptionIDColorModeVelocity=options->insertItem("&Velocity",this, SLOT(slotColorModeVelocity()),ALT+Key_V );
  OptionIDColorModeVelocityFixed=options->insertItem("Velocity&Fixed",this, SLOT(slotColorModeVelocityFixed()),ALT+Key_F );
  OptionIDColorModeAngularVelocity=options->insertItem("&AngularVelocity",this, SLOT(slotColorModeAngularVelocity()),ALT+Key_A );
  OptionIDColorModeVelocityNorm=options->insertItem("Velocity&Norm",this, SLOT(slotColorModeVelocityNorm()),ALT+Key_N );
  OptionIDColorModeVelocityNormFixed=options->insertItem("VelocityN&ormFixed",this, SLOT(slotColorModeVelocityNormFixed()),ALT+Key_O );
  OptionIDColorModeType=options->insertItem("&Type",this, SLOT(slotColorModeType()),ALT+Key_T );
  OptionIDColorModeId=options->insertItem("&Id",this, SLOT(slotColorModeId()),ALT+Key_I );
  OptionIDColorModeColor=options->insertItem("&Color",this, SLOT(slotColorModeColor()),ALT+Key_C );
  OptionIDColorModeForeground=options->insertItem("Fore&ground",this, SLOT(slotColorModeForeground()),ALT+Key_G );
  options->setCheckable(TRUE); 

  switch(m_pCGLBox->getDisplayMode()){
   case DISPLAYMODE_SPHERE:
    options->setItemChecked(OptionIDDisplayModeSphere,TRUE);
    break;
   case DISPLAYMODE_TEXTURE:
    options->setItemChecked(OptionIDDisplayModeTexture,TRUE);
    break;
   case DISPLAYMODE_CUBE:
    options->setItemChecked(OptionIDDisplayModeCube,TRUE);
    break;
   case DISPLAYMODE_CIRCLE:
    options->setItemChecked(OptionIDDisplayModeCircle,TRUE);
    break;
  }
// }


 pMain->statusBar();
 // Create a menu bar
 QMenuBar *m=pMain->menuBar();
 //m->setSeparator( QMenuBar::InWindowsStyle );
 m->insertItem("&File",file);
 m->insertItem("&Options",options);

if(!m_bPictures){
    int value;
    int min,max;

    // Create the three sliders; one for each rotation axis
    value=m_pCGLBox->getXRotation();
    QSlider *x=new QSlider(0,360,60,value,QSlider::Vertical,this,"xsl");
    x->setTickmarks(QSlider::Left);
    connect(x,SIGNAL(valueChanged(int)),m_pCGLBox,SLOT(setXRotation(int)));

    value=m_pCGLBox->getYRotation();
    QSlider *y=new QSlider(0,360,60,value,QSlider::Vertical,this,"ysl");
    y->setTickmarks(QSlider::Left);
    connect(y,SIGNAL(valueChanged(int)),m_pCGLBox,SLOT(setYRotation(int)));

    value=m_pCGLBox->getZRotation();
    QSlider *z=new QSlider (0,360,60,value,QSlider::Vertical,this,"zsl");
    z->setTickmarks(QSlider::Left);
    connect(z,SIGNAL(valueChanged(int)),m_pCGLBox,SLOT(setZRotation(int)));

    value=m_pCGLBox->getScale();
    min=0;
    max=100;
    if(value<min)min=value; 
    if(value>max)max=value; 
    QSlider *SScale=new QSlider(min,max,10,value,QSlider::Vertical,this,"scalesl");
    SScale->setTickmarks(QSlider::Left);
    connect(SScale,SIGNAL(valueChanged(int)),m_pCGLBox,SLOT(setScale(int)));

    value=m_pCGLBox->getViewScale();
    min=0;
    max=100;
    if(value<min)min=value; 
    if(value>max)max=value; 
    QSlider *SViewScale=new QSlider(min,max,10,value,QSlider::Vertical,this,"scalesl");
    SViewScale->setTickmarks(QSlider::Left);
    connect(SViewScale,SIGNAL(valueChanged(int)),m_pCGLBox,SLOT(setViewScale(int)));

    value=m_pCGLBox->getOffset();
    min=-100;
    max=0;
    if(value<min)min=value; 
    if(value>max)max=value; 
    QSlider *SOffset=new QSlider(min,max,10,value,QSlider::Vertical,this,"scalesl");
    SOffset->setTickmarks(QSlider::Left);
    connect(SOffset,SIGNAL(valueChanged(int)),m_pCGLBox,SLOT(setOffset(int)));


    // Now that we have all the widgets, put them into a nice layout

    // Put the sliders on top of each other
    QVBoxLayout *vlayout1=new QVBoxLayout(20,"vlayout1");
    vlayout1->addWidget(x);
    vlayout1->addWidget(y);
    vlayout1->addWidget(z);

    // Put the GL widget inside the frame
    QVBoxLayout *flayout1=new QVBoxLayout(20,"flayout1");
    flayout1->addStretch();
    flayout1->addWidget(m_pCGLBox);
    flayout1->addStretch();

    QVBoxLayout *vlayout2=new QVBoxLayout(20,"vlayout2");
    vlayout2->addWidget(SScale,1);
    vlayout2->addWidget(SViewScale,1);
    vlayout2->addWidget(SOffset,1);

    // Top level layout, puts the sliders to the left of the frame/GL widget
    QHBoxLayout *hlayout=new QHBoxLayout(this,20,20,"hlayout");
    hlayout->addLayout(vlayout1);
    hlayout->addLayout(flayout1);
    hlayout->addLayout(vlayout2);
 } 
}

void D_CDisplay::finish(T_CArg *pCArg)
{
 m_pCGLBox->finish(pCArg); 
}

void D_CDisplay::slotDisplayModeSphere(void)
{
 if(!m_bPictures){
 QMenuBar *m=m_pMain->menuBar();
 m->setItemChecked(OptionIDDisplayModeSphere,TRUE);
 m->setItemChecked(OptionIDDisplayModeTexture,FALSE);
 m->setItemChecked(OptionIDDisplayModeCube,FALSE);
 m->setItemChecked(OptionIDDisplayModeCircle,FALSE);
 m->setItemChecked(OptionIDDisplayModeBlip,FALSE);
 m_pCGLBox->setDisplayMode(DISPLAYMODE_SPHERE);
}
}

void D_CDisplay::slotDisplayModeTexture(void)
{
 if(!m_bPictures){
 QMenuBar *m=m_pMain->menuBar();
 m->setItemChecked(OptionIDDisplayModeSphere,FALSE);
 m->setItemChecked(OptionIDDisplayModeTexture,TRUE);
 m->setItemChecked(OptionIDDisplayModeCube,FALSE);
 m->setItemChecked(OptionIDDisplayModeCircle,FALSE);
 m->setItemChecked(OptionIDDisplayModeBlip,FALSE);
 m_pCGLBox->setDisplayMode(DISPLAYMODE_TEXTURE);
}
}

void D_CDisplay::slotDisplayModeCube(void)
{
 if(!m_bPictures){
 QMenuBar *m=m_pMain->menuBar();
 m->setItemChecked(OptionIDDisplayModeSphere,FALSE);
 m->setItemChecked(OptionIDDisplayModeTexture,FALSE);
 m->setItemChecked(OptionIDDisplayModeCube,TRUE);
 m->setItemChecked(OptionIDDisplayModeCircle,FALSE);
 m->setItemChecked(OptionIDDisplayModeBlip,FALSE);
 m_pCGLBox->setDisplayMode(DISPLAYMODE_CUBE);
}
}

void D_CDisplay::slotDisplayModeCircle(void)
{
 if(!m_bPictures){
 QMenuBar *m=m_pMain->menuBar();
 m->setItemChecked(OptionIDDisplayModeSphere,FALSE);
 m->setItemChecked(OptionIDDisplayModeTexture,FALSE);
 m->setItemChecked(OptionIDDisplayModeCube,FALSE);
 m->setItemChecked(OptionIDDisplayModeCircle,TRUE);
 m->setItemChecked(OptionIDDisplayModeBlip,FALSE);
 m_pCGLBox->setDisplayMode(DISPLAYMODE_CIRCLE);
}
}

void D_CDisplay::slotDisplayModeBlip(void)
{
 if(!m_bPictures){
 QMenuBar *m=m_pMain->menuBar();
 m->setItemChecked(OptionIDDisplayModeSphere,FALSE);
 m->setItemChecked(OptionIDDisplayModeTexture,FALSE);
 m->setItemChecked(OptionIDDisplayModeCube,FALSE);
 m->setItemChecked(OptionIDDisplayModeCircle,FALSE);
 m->setItemChecked(OptionIDDisplayModeBlip,TRUE);
 m_pCGLBox->setDisplayMode(DISPLAYMODE_BLIP);
}
}


void D_CDisplay::slotColorModeVelocity(void)
{
 if(!m_bPictures){
 QMenuBar *m=m_pMain->menuBar();
 m->setItemChecked(OptionIDColorModeVelocity,TRUE);
 m->setItemChecked(OptionIDColorModeVelocityFixed,FALSE);
 m->setItemChecked(OptionIDColorModeAngularVelocity,FALSE);
 m->setItemChecked(OptionIDColorModeVelocityNorm,FALSE);
 m->setItemChecked(OptionIDColorModeVelocityNormFixed,FALSE);
 m->setItemChecked(OptionIDColorModeType,FALSE);
 m->setItemChecked(OptionIDColorModeId,FALSE);
 m->setItemChecked(OptionIDColorModeColor,FALSE);
 m->setItemChecked(OptionIDColorModeForeground,FALSE);
 m_pCGLBox->setColorMode(COLORMODE_VELOCITY);
}
}

void D_CDisplay::slotColorModeVelocityFixed(void)
{
 if(!m_bPictures){
 QMenuBar *m=m_pMain->menuBar();
 m->setItemChecked(OptionIDColorModeVelocity,FALSE);
 m->setItemChecked(OptionIDColorModeVelocityFixed,TRUE);
 m->setItemChecked(OptionIDColorModeAngularVelocity,FALSE);
 m->setItemChecked(OptionIDColorModeVelocityNorm,FALSE);
 m->setItemChecked(OptionIDColorModeVelocityNormFixed,FALSE);
 m->setItemChecked(OptionIDColorModeType,FALSE);
 m->setItemChecked(OptionIDColorModeId,FALSE);
 m->setItemChecked(OptionIDColorModeColor,FALSE);
 m->setItemChecked(OptionIDColorModeForeground,FALSE);
 m_pCGLBox->setColorMode(COLORMODE_VELOCITYFIXED);
}
}

void D_CDisplay::slotColorModeAngularVelocity(void)
{
 if(!m_bPictures){
 QMenuBar *m=m_pMain->menuBar();
 m->setItemChecked(OptionIDColorModeVelocity,FALSE);
 m->setItemChecked(OptionIDColorModeVelocityFixed,FALSE);
 m->setItemChecked(OptionIDColorModeAngularVelocity,TRUE);
 m->setItemChecked(OptionIDColorModeVelocityNorm,FALSE);
 m->setItemChecked(OptionIDColorModeVelocityNormFixed,FALSE);
 m->setItemChecked(OptionIDColorModeType,FALSE);
 m->setItemChecked(OptionIDColorModeId,FALSE);
 m->setItemChecked(OptionIDColorModeColor,FALSE);
 m->setItemChecked(OptionIDColorModeForeground,FALSE);
 m_pCGLBox->setColorMode(COLORMODE_ANGULARVELOCITY);
}
}

void D_CDisplay::slotColorModeVelocityNorm(void)
{
 if(!m_bPictures){
 QMenuBar *m=m_pMain->menuBar();
 m->setItemChecked(OptionIDColorModeVelocity,FALSE);
 m->setItemChecked(OptionIDColorModeVelocityFixed,FALSE);
 m->setItemChecked(OptionIDColorModeAngularVelocity,FALSE);
 m->setItemChecked(OptionIDColorModeVelocityNorm,TRUE);
 m->setItemChecked(OptionIDColorModeVelocityNormFixed,FALSE);
 m->setItemChecked(OptionIDColorModeType,FALSE);
 m->setItemChecked(OptionIDColorModeId,FALSE);
 m->setItemChecked(OptionIDColorModeColor,FALSE);
 m->setItemChecked(OptionIDColorModeForeground,FALSE);
 m_pCGLBox->setColorMode(COLORMODE_VELOCITYNORM);
}
}

void D_CDisplay::slotColorModeVelocityNormFixed(void)
{
 if(!m_bPictures){
 QMenuBar *m=m_pMain->menuBar();
 m->setItemChecked(OptionIDColorModeVelocity,FALSE);
 m->setItemChecked(OptionIDColorModeVelocityFixed,FALSE);
 m->setItemChecked(OptionIDColorModeAngularVelocity,FALSE);
 m->setItemChecked(OptionIDColorModeVelocityNorm,FALSE);
 m->setItemChecked(OptionIDColorModeVelocityNormFixed,TRUE);
 m->setItemChecked(OptionIDColorModeType,FALSE);
 m->setItemChecked(OptionIDColorModeId,FALSE);
 m->setItemChecked(OptionIDColorModeColor,FALSE);
 m->setItemChecked(OptionIDColorModeForeground,FALSE);
 m_pCGLBox->setColorMode(COLORMODE_VELOCITYNORMFIXED);
}
}

void D_CDisplay::slotColorModeType(void)
{
 if(!m_bPictures){
 QMenuBar *m=m_pMain->menuBar();
 m->setItemChecked(OptionIDColorModeVelocity,FALSE);
 m->setItemChecked(OptionIDColorModeVelocityFixed,FALSE);
 m->setItemChecked(OptionIDColorModeAngularVelocity,FALSE);
 m->setItemChecked(OptionIDColorModeVelocityNorm,FALSE);
 m->setItemChecked(OptionIDColorModeVelocityNormFixed,FALSE);
 m->setItemChecked(OptionIDColorModeType,TRUE);
 m->setItemChecked(OptionIDColorModeId,FALSE);
 m->setItemChecked(OptionIDColorModeColor,FALSE);
 m->setItemChecked(OptionIDColorModeForeground,FALSE);
 m_pCGLBox->setColorMode(COLORMODE_TYPE);
}
}

void D_CDisplay::slotColorModeId(void)
{
 if(!m_bPictures){
 QMenuBar *m=m_pMain->menuBar();
 m->setItemChecked(OptionIDColorModeVelocity,FALSE);
 m->setItemChecked(OptionIDColorModeVelocityFixed,FALSE);
 m->setItemChecked(OptionIDColorModeAngularVelocity,FALSE);
 m->setItemChecked(OptionIDColorModeVelocityNorm,FALSE);
 m->setItemChecked(OptionIDColorModeVelocityNormFixed,FALSE);
 m->setItemChecked(OptionIDColorModeType,FALSE);
 m->setItemChecked(OptionIDColorModeId,TRUE);
 m->setItemChecked(OptionIDColorModeColor,FALSE);
 m->setItemChecked(OptionIDColorModeForeground,FALSE);
 m_pCGLBox->setColorMode(COLORMODE_ID);
}
}

void D_CDisplay::slotColorModeColor(void)
{
 if(!m_bPictures){
 QMenuBar *m=m_pMain->menuBar();
 m->setItemChecked(OptionIDColorModeVelocity,FALSE);
 m->setItemChecked(OptionIDColorModeVelocityFixed,FALSE);
 m->setItemChecked(OptionIDColorModeAngularVelocity,FALSE);
 m->setItemChecked(OptionIDColorModeVelocityNorm,FALSE);
 m->setItemChecked(OptionIDColorModeVelocityNormFixed,FALSE);
 m->setItemChecked(OptionIDColorModeType,FALSE);
 m->setItemChecked(OptionIDColorModeId,FALSE);
 m->setItemChecked(OptionIDColorModeColor,TRUE);
 m->setItemChecked(OptionIDColorModeForeground,FALSE);
 m_pCGLBox->setColorMode(COLORMODE_COLOR);
}
}

void D_CDisplay::slotColorModeForeground(void)
{
 if(!m_bPictures){
 QMenuBar *m=m_pMain->menuBar();
 m->setItemChecked(OptionIDColorModeVelocity,FALSE);
 m->setItemChecked(OptionIDColorModeVelocityFixed,FALSE);
 m->setItemChecked(OptionIDColorModeAngularVelocity,FALSE);
 m->setItemChecked(OptionIDColorModeVelocityNorm,FALSE);
 m->setItemChecked(OptionIDColorModeVelocityNormFixed,FALSE);
 m->setItemChecked(OptionIDColorModeType,FALSE);
 m->setItemChecked(OptionIDColorModeId,FALSE);
 m->setItemChecked(OptionIDColorModeColor,FALSE);
 m->setItemChecked(OptionIDColorModeForeground,TRUE);
 m_pCGLBox->setColorMode(COLORMODE_FOREGROUND);
}
}

