#include "D_CLoader.h"

D_CLoader::D_CLoader(void)
{
 m_file=NULL;
#if defined(WITH_QT)
 m_pTimer=NULL;
#endif
 m_pArray=NULL;
 m_ArrayLength=0;
 m_pHeader=new D_CHeader();
 m_pHeader->Number=0;
 m_pHeader->Time=0.;
 m_pHeader->xMin=0.;
 m_pHeader->yMin=0.;
 m_pHeader->zMin=0.;
 m_pHeader->xMax=0.;
 m_pHeader->yMax=0.;
 m_pHeader->zMax=0.;
 m_pDataFormat=NULL;
 m_pHeaderFormat=NULL;
 m_pArrayLine=NULL;
 m_ArrayLineLength=0;
 m_pHeaderLine=new D_CHeaderLine();
 m_pHeaderLine->Number=0;
}

D_CLoader::~D_CLoader(void)
{
 if(m_file!=stdin)fclose(m_file);
#if defined(WITH_QT)
 if(m_pTimer){
  m_pTimer->stop();
  T_DELETE(m_pTimer)
 }
#endif
 T_DELETE(m_pArray)
 T_DELETE(m_pDataFormat)
 T_DELETE(m_pHeaderFormat)
}

#if defined(WITH_QT)
bool D_CLoader::init(QStatusBar *pStatusBar,D_CGLBox *pBox,T_CArg *pCArg)
#else
bool D_CLoader::init(D_CGLX *pGLX,T_CArg *pCArg)
#endif
{
 char *name;
 pCArg->getParameter("main","file",&name);
 if(!strcmp(name,"-")){
  m_file=stdin;
 } else {
  m_file=fopen(name,"r");
 }

 if(m_file){
#if defined(WITH_QT)
  m_pStatusBar=pStatusBar;
  m_pBox=pBox;
#else
  m_pGLX=pGLX;
#endif
  pCArg->getParameter("main","tmin",&m_btMin);
  pCArg->getParameter("main","tmax",&m_btMax);
  pCArg->getParameter("main","tmin",&m_tMin);
  pCArg->getParameter("main","tmax",&m_tMax);
  pCArg->getParameter("main","step",&m_Step);
  m_it=0l;
#if defined(WITH_QT)
  pCArg->getParameter("output","flag",&m_bPictures);

  m_pTimer=new QTimer();
  connect(m_pTimer,SIGNAL(timeout()),SLOT(load()));
  pCArg->getParameter("main","timer",&m_timestep);
  if(m_btMin && m_tMin>m_pHeader->Time){
   m_pTimer->start(1,FALSE);
  } else {
   m_pTimer->start(m_timestep,FALSE);
  }
#endif
  char *string;
  pCArg->getParameter("main","header",&string);
  initHeaderFormat(string);
  pCArg->getParameter("main","format",&string);
  initDataFormat(string);
  return TRUE;
 } else return FALSE;
}

void D_CLoader::initDataFormat(char *pFormat)
{
 int imax;
 int x,v,r,q,w,id,t;
 bool bFlag=TRUE;
 T_N_t pos;

 if(!strcmp(pFormat,"7")){
  m_format=7;
  return;
 }
 if(!strcmp(pFormat,"8")){
  m_format=8;
  return;
 }
 if(!strcmp(pFormat,"100")){
  m_format=100;
  return;
 }
 if(!strcmp(pFormat,"101")){
  m_format=101;
  return;
 }
 if(!strcmp(pFormat,"200")){
  m_format=200;
  return;
 }
 m_format=0;

 x=v=r=q=w=id=t=0;
 imax=strlen(pFormat);

 pos=0;
 for(int i=0;i<imax;++i){
  switch(pFormat[i]){
   case '-' :
    pos++;
    break;
   case 'x' :
    if(x<3){
     x++;
     pos++;
    }
    break;
   case 'v' : 
    if(v<3){
     v++;
     pos++;
    }
    break;
   case 'r' : 
    if(r<1){
     r++;
     pos++;
    }
    break;
   case 'q' : 
    if(q<4){
     q++;
     pos++;
    }
    break;
   case 'w' : 
    if(w<3){
     w++;
     pos++;
    }
    break;
  case 'i' : 
    if(id<1){
     id++;
     pos++;
    }
    break;
   case 't' : 
    if(t<1){
     t++;
     pos++;
    }
    break;
  }
 }
 m_DataFormatPosMax=pos;

 T_MALLOC(m_pDataFormat,m_DataFormatPosMax,T_N_t,bFlag)
 if(!bFlag)return;

 x=v=r=q=w=id=t=0;
 pos=0;
 for(int i=0;i<imax;++i){
  switch(pFormat[i]){
   case '-' :
    m_pDataFormat[pos]=0;
    pos++;
    break;
   case 'x' :
    if(x<3){
     m_pDataFormat[pos]=1+x;
     x++;
     pos++;
    }
    break;
   case 'v' : 
    if(v<3){
     m_pDataFormat[pos]=4+v;
     v++;
     pos++;
    }
    break;
   case 'r' : 
    if(r<1){
     m_pDataFormat[pos]=7;
     r++;
     pos++;
    }
    break;
   case 'q' : 
    if(q<4){
     m_pDataFormat[pos]=8+q;
     q++;
     pos++;
    }
    break;
   case 'w' : 
    if(w<3){
     m_pDataFormat[pos]=12+w;
     w++;
     pos++;
    }
    break;
   case 'i' : 
    if(id<1){
     m_pDataFormat[pos]=15;
     id++;
     pos++;
    }
    break;
   case 't' :
    if(t<1){ 
     m_pDataFormat[pos]=16;
     t++;
     pos++;
    }
    break;
  }
 }
}

void D_CLoader::initHeaderFormat(char *pFormat)
{
 int imax;
 int x,y,z,X,Y,Z,n,t;
 bool bFlag=TRUE;
 T_N_t pos;

 if(!strcmp(pFormat,"7")){
  m_headerformat=7;
  return;
 }
 if(!strcmp(pFormat,"8")){
  m_headerformat=8;
  return;
 }
 if(!strcmp(pFormat,"200")){
  m_headerformat=200;
  return;
 }
 m_headerformat=0;

 x=y=z=X=Y=Z=n=t=0;
 imax=strlen(pFormat);

 pos=0;
 for(int it=0;it<imax;++it){
  switch(pFormat[it]){
   case '-' :
    pos++;
    break;
   case 'n' : 
    if(n<1){
     n++;
     pos++;
    }
    break;
   case 't' : 
    if(t<1){
     t++;
     pos++;
    }
    break;
   case 'x' :
    if(x<1){
     x++;
     pos++;
    }
    break;
   case 'y' :
    if(y<1){
     y++;
     pos++;
    }
    break;
   case 'z' :
    if(z<1){
     z++;
     pos++;
    }
    break;
   case 'X' :
    if(X<1){
     X++;
     pos++;
    }
    break;
   case 'Y' :
    if(Y<1){
     Y++;
     pos++;
    }
    break;
   case 'Z' :
    if(Z<1){
     Z++;
     pos++;
    }
    break;
  }
 }
 m_HeaderFormatPosMax=pos;

 T_MALLOC(m_pHeaderFormat,m_HeaderFormatPosMax,T_N_t,bFlag)
 if(!bFlag)return;

 x=y=z=X=Y=Z=n=t=0;
 pos=0;
 for(int it=0;it<imax;++it){
  switch(pFormat[it]){
   case '-' :
    m_pHeaderFormat[pos]=0;
    pos++;
    break;
   case 'n' : 
    if(n<1){
     m_pHeaderFormat[pos]=1;
     n++;
     pos++;
    }
    break;
   case 't' : 
    if(t<1){
     m_pHeaderFormat[pos]=2;
     t++;
     pos++;
    }
    break;
   case 'x' :
    if(x<1){
     m_pHeaderFormat[pos]=3;
     x++;
     pos++;
    }
    break;
   case 'y' :
    if(y<1){
     m_pHeaderFormat[pos]=4;
     y++;
     pos++;
    }
    break;
   case 'z' :
    if(z<1){
     m_pHeaderFormat[pos]=5;
     z++;
     pos++;
    }
    break;
   case 'X' :
    if(X<1){
     m_pHeaderFormat[pos]=6;
     X++;
     pos++;
    }
    break;
   case 'Y' :
    if(Y<1){
     m_pHeaderFormat[pos]=7;
     Y++;
     pos++;
    }
    break;
   case 'Z' :
    if(Z<1){
     m_pHeaderFormat[pos]=8;
     Z++;
     pos++;
    }
    break;
  }
 }
}


void D_CLoader::load(void)
{
 T_R3_t dummy;
  T_N_t pos;

 char c,cc;
 int fd;
 int Result_GETFL,Result;

#if defined(WITH_QT)
 if(m_pBox->stop())return;
#endif

 fd=fileno(m_file);
 Result_GETFL=fcntl(fd,F_GETFL,0) ;
 if(Result_GETFL<0){
  fprintf(stderr,"Error: fcntl(..,F_GETFL,..)\n");
  return;
 }
 Result=fcntl(fd,F_SETFL,Result_GETFL | O_NONBLOCK) ;
 if(Result<0){
  fcntl(fd,F_SETFL,Result_GETFL) ;
  fprintf(stderr,"Error: fcntl(..,F_SETFL,..)\n");
  return;
 }

 c=fgetc(m_file);

 if(m_headerformat<200){
  // eliminiert Zeilen mit Kommentaren (nur fuer ascii formate sinnvoll)
  while(c=='#'){
   do{
    cc=fgetc(m_file);
   }while(cc!='\n' && cc!=EOF);
   if(cc==EOF)c=EOF;
   else c=fgetc(m_file);
  }
 }

 /*
  stellt den Stream wieder zur"uck
 */
 Result=fcntl(fd,F_SETFL,Result_GETFL) ;
 if(Result<0){
  fprintf(stderr,"Error: fcntl(..,F_SETFL,..)\n");
  return;
 }

 if(c==EOF){
  return;
 }
 ungetc(c,m_file);
 
 switch(m_headerformat){
  default:

   memset(m_pHeader,0,sizeof(D_Header_t));
   for(pos=0;pos<m_HeaderFormatPosMax;++pos){
    switch(m_pHeaderFormat[pos]){
     default:
     case 0:
      fscanf(m_file,"%lf",&dummy);
      break;
     case 1 :
      fscanf(m_file,"%i",&(m_pHeader->Number));
      break;
     case 2 :
      fscanf(m_file,"%lf",&(m_pHeader->Time));
      break;
     case 3 :
      fscanf(m_file,"%lf",&(m_pHeader->xMin));
      break;
     case 4 :
      fscanf(m_file,"%lf",&(m_pHeader->yMin));
      break;
     case 5 :
      fscanf(m_file,"%lf",&(m_pHeader->zMin));
      break;
     case 6 :
      fscanf(m_file,"%lf",&(m_pHeader->xMax));
      break;
     case 7 :
      fscanf(m_file,"%lf",&(m_pHeader->yMax));
      break;
     case 8 :
      fscanf(m_file,"%lf",&(m_pHeader->zMax));
      break;
    }
   }
   while(fgetc(m_file)!='\n');
  break;
 case 7:
   memset(m_pHeader,0,sizeof(D_Header_t));
   fscanf(m_file,"%i %lf %lf %lf %lf %lf %lf %lf",
    &m_pHeader->Number, &m_pHeader->Time,
    &m_pHeader->xMin, &m_pHeader->yMin, &m_pHeader->zMin,
    &m_pHeader->xMax, &m_pHeader->yMax, &m_pHeader->zMax);
   break;
  case 8:
   memset(m_pHeader,0,sizeof(D_Header_t));
   fscanf(m_file,"%i %lf %lf %lf %lf %lf",
    &m_pHeader->Number, &m_pHeader->Time,
    &m_pHeader->xMin, &m_pHeader->yMin,
    &m_pHeader->xMax, &m_pHeader->yMax);
   break;
  case 200:
   memset(m_pHeader,0,sizeof(D_Header_t));
   fread(m_pHeader,sizeof(D_Header_t),1,m_file);
//   fprintf(stderr,"%i %e\n",m_pHeader->Number,m_pHeader->Time);
   break;
 }

 if(m_format<200){
  // eliminiert restliche white-space (nur fuer ascii formate)
  while(fgetc(m_file)!='\n');
 }

  // alloziert wenn n"otig mehr Speicher
  if(m_pHeader->Number>m_ArrayLength){
   D_Data_t *pArray;
   pArray=(D_Data_t *)realloc(m_pArray,m_pHeader->Number*sizeof(D_Data_t));
   if(pArray){
    m_pArray=pArray;
    m_ArrayLength=m_pHeader->Number;
   } else if(m_pArray){
    free(m_pArray);
    m_ArrayLength=0;
   }
  }
  if(!m_pArray)return;
  // lie"st Daten ein

  if(m_format<200){
   // eliminiert Zeilen mit Kommentaren (nur fuer ascii Formate sinnvoll)
   c=fgetc(m_file);
   while(c=='#'){
    while(fgetc(m_file)!='\n');
    c=fgetc(m_file);
   }
   ungetc(c,m_file);
  }

  memset(m_pArray,0,m_pHeader->Number*sizeof(D_Data_t));
  switch(m_format){
   default:
    for(int i=0;i<m_pHeader->Number;++i){
     memset(&(m_pArray[i]),0,sizeof(D_Data_t));
     for(pos=0;pos<m_DataFormatPosMax;++pos){
      switch(m_pDataFormat[pos]){
       default:
       case 0:
        fscanf(m_file,"%lf",&dummy);
        break;
       case 1 :
        fscanf(m_file,"%lf",&(m_pArray[i].x));
        break;
       case 2 :
        fscanf(m_file,"%lf",&(m_pArray[i].y));
        break;
       case 3 :
        fscanf(m_file,"%lf",&(m_pArray[i].z));
        break;
       case 4 :
        fscanf(m_file,"%lf",&(m_pArray[i].vx));
        break;
       case 5 :
        fscanf(m_file,"%lf",&(m_pArray[i].vy));
        break;
       case 6 :
        fscanf(m_file,"%lf",&(m_pArray[i].vz));
        break;
       case 7 :
        fscanf(m_file,"%lf",&(m_pArray[i].r));
        break;
       case 8 :
        fscanf(m_file,"%lf",&(m_pArray[i].qw));
        break;
       case 9 :
        fscanf(m_file,"%lf",&(m_pArray[i].qx));
        break;
       case 10 :
        fscanf(m_file,"%lf",&(m_pArray[i].qy));
        break;
       case 11 :
        fscanf(m_file,"%lf",&(m_pArray[i].qz));
        break;
       case 12 :
        fscanf(m_file,"%lf",&(m_pArray[i].wx));
        break;
       case 13 :
        fscanf(m_file,"%lf",&(m_pArray[i].wy));
        break;
       case 14 :
        fscanf(m_file,"%lf",&(m_pArray[i].wz));
        break;
       case 15 :
        fscanf(m_file,"%i",&(m_pArray[i].id));
        break;
       case 16 :
        fscanf(m_file,"%i",&(m_pArray[i].type));
        break;
      }
     }
     while(fgetc(m_file)!='\n');
    }
    break;
   case 7:
    for(int i=0;i<m_pHeader->Number;++i){
     fscanf(m_file,"%lf %lf %lf %lf %lf %lf %lf",
      &(m_pArray[i].x), &(m_pArray[i].y), &(m_pArray[i].z),
      &(m_pArray[i].vx), &(m_pArray[i].vy), &(m_pArray[i].vz),
      &(m_pArray[i].r));
     while(fgetc(m_file)!='\n');
    }
    break;
   case 8:
    for(int i=0;i<m_pHeader->Number;++i){
     fscanf(m_file,"%lf %lf %lf %lf %lf",
      &(m_pArray[i].x), &(m_pArray[i].y),
      &(m_pArray[i].vx), &(m_pArray[i].vy),
      &(m_pArray[i].r));
     while(fgetc(m_file)!='\n');
    }
    break;
   case 100:
    for(int i=0;i<m_pHeader->Number;++i){
     fscanf(m_file,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %i %i %i",
      &(m_pArray[i].x), &(m_pArray[i].y), &(m_pArray[i].z),
      &(m_pArray[i].vx), &(m_pArray[i].vy), &(m_pArray[i].vz),
      &(m_pArray[i].r),
      &(m_pArray[i].qw), &(m_pArray[i].qx), &(m_pArray[i].qy), &(m_pArray[i].qz),
      &(m_pArray[i].wx), &(m_pArray[i].wy), &(m_pArray[i].wz),
      &(m_pArray[i].id),&(m_pArray[i].type),&(m_pArray[i].color));
     while(fgetc(m_file)!='\n');

    }
    break;
   case 101:
    for(int i=0;i<m_pHeader->Number;++i){
     fscanf(m_file,"%lf %lf %lf %lf %lf %lf %lf %i %i",
      &(m_pArray[i].x), &(m_pArray[i].y), &(m_pArray[i].z),
      &(m_pArray[i].vx), &(m_pArray[i].vy), &(m_pArray[i].vz),
      &(m_pArray[i].r),
      &(m_pArray[i].id),&(m_pArray[i].type));
     while(fgetc(m_file)!='\n');
    }
    break;
   case 200:
    fread(m_pArray,sizeof(D_Data_t),m_pHeader->Number,m_file);
//fprintf(stderr," 0: %e %i\n",m_pArray[0].x,m_pArray[0].id);
//fprintf(stderr,"10: %e %i\n",m_pArray[10].x,m_pArray[10].id);
    break;
  }


///////////////////////////////////////////////////////////////
 if(m_headerformat<200){

 fd=fileno(m_file);
 Result_GETFL=fcntl(fd,F_GETFL,0) ;
 if(Result_GETFL<0){
  fprintf(stderr,"Error: fcntl(..,F_GETFL,..)\n");
  return;
 }
 Result=fcntl(fd,F_SETFL,Result_GETFL | O_NONBLOCK) ;
 if(Result<0){
  fcntl(fd,F_SETFL,Result_GETFL) ;
  fprintf(stderr,"Error: fcntl(..,F_SETFL,..)\n");
  return;
 }

 c=fgetc(m_file);

 Result=fcntl(fd,F_SETFL,Result_GETFL) ;
 if(Result<0){
  fprintf(stderr,"Error: fcntl(..,F_SETFL,..)\n");
  return;
 }

 if(c==EOF || c!='#'){
  if(c!=EOF){
   ungetc(c,m_file);
   m_pHeaderLine->Number=0;
  }
 } else {
  c=fgetc(m_file);
  if(c!='L'){
   while(fgetc(m_file)!='\n');
   m_pHeaderLine->Number=0;
  } else {
   /// Linien einlesen
   fscanf(m_file,"%i",&(m_pHeaderLine->Number));
   while(fgetc(m_file)!='#');
   // alloziert wenn n"otig mehr Speicher
   if(m_pHeaderLine->Number>m_ArrayLineLength){
    D_DataLine_t *pArray;
    pArray=(D_DataLine_t *)realloc(m_pArrayLine,m_pHeaderLine->Number*sizeof(D_DataLine_t));
    if(pArray){
     m_pArrayLine=pArray;
     m_ArrayLineLength=m_pHeaderLine->Number;
    } else if(m_pArrayLine){
     free(m_pArrayLine);
     m_ArrayLineLength=0;
    }
   }
   if(!m_pArrayLine)return;
 
   for(int i=0;i<m_pHeaderLine->Number;++i){
    fscanf(m_file,"%lf %lf %lf %lf %lf %lf %i",
    &(m_pArrayLine[i].x1), &(m_pArrayLine[i].y1), &(m_pArrayLine[i].z1),
    &(m_pArrayLine[i].x2), &(m_pArrayLine[i].y2), &(m_pArrayLine[i].z2),
    &(m_pArrayLine[i].color));
    while(fgetc(m_file)!='#');
   }
///
  }
 }
 }
 

///////////////////////////////////////////////////////////////

#if defined(WITH_QT)
 QString s;
 s.sprintf("%i Particles, Time: %f",m_pHeader->Number,m_pHeader->Time);
 m_pStatusBar->message(s);
#endif

 if(m_it==0){
  if(!((m_btMin && m_tMin>m_pHeader->Time) 
   || (m_btMax && m_tMax<m_pHeader->Time))){
#if defined(WITH_QT)
   m_pBox->setInterface(m_pHeader,m_pArray,m_pHeaderLine,m_pArrayLine);
   if(m_bPictures){
    m_pBox->updateList();
    m_pBox->updateGL();
    m_pBox->savePicture();
   } else {
    m_pBox->updateList();
    m_pBox->updateGL();
   }
#else
   m_pGLX->setInterface(m_pHeader,m_pArray);
#endif
  }
 } 
 m_it++;
 if(m_it>=m_Step)m_it=0;

#if defined(WITH_QT)
 if(m_btMin && m_tMin<=m_pHeader->Time){
  m_pTimer->changeInterval(m_timestep);
 }
 if((m_btMax && m_tMax<=m_pHeader->Time)){
  m_pTimer->stop();
 }
#endif
} 

void D_CLoader::rewind(void)
{
 if(m_file!=stdin)std::rewind(m_file);
}

