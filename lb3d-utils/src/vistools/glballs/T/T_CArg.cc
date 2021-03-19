#include "T_CArg.h"

// Defaultkonstruktor
/*
*/
T_CArgItem::T_CArgItem(const char *Class,const char *Name,const char *ShortName,const char *Value,const char *Description)
{
 bool bFlag=TRUE;
 T_STRDUP(m_Class,Class,bFlag)
 T_STRDUP(m_Name,Name,bFlag)
 T_STRDUP(m_ShortName,ShortName,bFlag)
 T_STRDUP(m_Value,Value,bFlag)
 T_STRDUP(m_Description,Description,bFlag)
 if(!bFlag)fprintf(stderr,"Error : constructor T_CArgItem\n");
}

// Defaultdestruktor
/*
*/
T_CArgItem::~T_CArgItem(void)
{
 T_FREE(m_Class);
 T_FREE(m_Name);
 T_FREE(m_ShortName);
 T_FREE(m_Value);
 T_FREE(m_Description);
}

// setzt den Wert dieses Parameters
/*
*/
bool T_CArgItem::set(const char *Value)
{
 bool bFlag=TRUE;
 if(m_Value)T_FREE(m_Value);
 T_STRDUP(m_Value,Value,bFlag);
 return bFlag;
}

// gibt den Wert dieses Parameters aus
/*
*/
char *T_CArgItem::get(void) const
{
 return m_Value;
}

// schreibt den Wert dieses Parameters in einen Stream
/*
*/
void T_CArgItem::describe(FILE *out)
{
 fprintf(out,"-%s::%s",m_Class,m_Name);
 if(m_Value)fprintf(out,"=%s",m_Value);
 if(m_Description)fprintf(out,"\t\t%s",m_Description);
 if(m_ShortName)fprintf(out," (-%s)",m_ShortName);
 fprintf(out,"\n");
}

// schreibt die Kurzbeschreibung des Parameters in einen Stream
/*
*/
void T_CArgItem::dump(FILE *out)
{
 if(m_Description)fprintf(out,"# %s\n",m_Description);
 if(m_Value)fprintf(out,"%s::%s=%s\n",m_Class,m_Name,m_Value);
}

// pr"uft ob der Parameter auf die angegebene Klasse und den angegebenen Namen pa"st
/*!
 \sa T_CArgItem::hasName()
*/
bool T_CArgItem::isOption(const char *Class,const char *Name) const
{
 if(isClass(Class)){
  if(hasName(Name)){
   return TRUE;
  }
 }
 return FALSE;
}

// pr"uft ob der Parameter auf die angegebene Klasse pa"st
/*
*/
bool T_CArgItem::isClass(const char *Class) const
{
 return !strcmp(Class,m_Class);
}

// pr"uft ob der Parameter auf den angegebenen Namen pa"st
/*!
 Dabei wird sowohl der normale Name als auch der Kurzname verglichen
*/
bool T_CArgItem::hasName(const char *Name) const
{
 if((!m_Name) || (!Name))return FALSE;
 if(!strcmp(Name,m_Name))return TRUE;
 if(!m_ShortName)return FALSE;
 if(!strcmp(Name,m_ShortName))return TRUE;
 return FALSE;
}

/*************************************************************************/

// Defaultkonstruktor
/*!
 Der "ubergebene String wird als Defaultklasse ben"utzt. Also jedes mal
 wenn beim registrieren ioder parsen eines Parameters keine Klasse angegeben 
 wurde.
*/
T_CArg::T_CArg(const char *DefaultClass)
{
 bool bFlag=TRUE;
 if(DefaultClass){
  T_STRDUP(m_DefaultClass,DefaultClass,bFlag)
 }
 registerParameter("arg","dump",NULL,"false","dump parameters into file");
 registerParameter("arg","file",NULL,"false","use file instead of commandline");
 m_CCAI.setMaster();
 if(!bFlag)fprintf(stderr,"Error : constructor T_CArg\n");
}

// Defaultdestruktur
/*!
 wurde der parameter -arg::dump=Dateiname gesetzt, werden die gespeicherten 
 Parameter in die entsprechende Datei geschrieben.
*/
T_CArg::~T_CArg(void)
{
 dumpArgs();
 T_FREE(m_DefaultClass);
}

// registriert einen Parameter
/*!
 F"ugt einen neuen Parameter in die Liste ein. 
 Wird die Klasse ausgelassen, wird statt dessen der Defaultwert ben"utzt.
 \warning Es wird nicht auf Doppelbenennungen "uberpr"uft !
*/
void T_CArg::registerParameter(const char *Class,const char *Name,const char *ShortName,const char *Value,const char *Description)
{
 T_CArgItem *pCAI;
 if(!Class)Class=m_DefaultClass;
 if(!Name)return;
 T_NEW(pCAI,T_CArgItem,(Class,Name,ShortName,Value,Description));
 if(pCAI)pCAI->linkWith(&m_CCAI);
}

// parst die von der main() Funktion stammenden Kommandozeilenparameter
/*!
 Dabei wird zun"achst kontrolliert, ob eines der Argumente -?, -h, --help
 "ubergeben wurde. Falls dem so ist, wird eine Beschreibung der vorhandenen
 Parameter ausgegeben und dem Programm eine Aufforderung zum Abbruch
 zur"uckgemeldet. 

 Andernfalls wird erst "uber die Kommandozeile ermittelt, ob eine
 Konfigurationsdatei existiert und diese gegebenenfalls geparst.
 Erst dann werden die Kommandozeilenparameter endg"ultig "ubernommen.
*/
bool T_CArg::parseArgs(int argc,char *argv[])
{
 int i;

 for(i=0;i<argc;++i){
  if((!strcmp(argv[i],"-?")) || (!strcmp(argv[i],"-h")) || (!strcmp(argv[i],"--help"))){
   describe();
   return FALSE;
  } 
 }
 parseArgsInternal(argc,argv,TRUE);
 parseFile(TRUE);
 parseArgsInternal(argc,argv,FALSE);
// dumpArgs();
 return TRUE;
}

// parst eine Datei
/*
 \warning eine Zeile darf nicht mehr als 511 Zeichen haben
*/
void T_CArg::parseFile(bool Flag)
{
 char *Class,*Name,*Value;
 char *posLeft,*posRight;
 char *buffer;
 FILE *in;
 char *filename;
 bool bFlag=TRUE;

 buffer=NULL;
 Class=NULL;
 Name=NULL;
 Value=NULL;

 T_MALLOC(buffer,512,char,bFlag)
 if(!bFlag)return;

 if(!isValidParameter("arg","file")){
  T_FREE(buffer);
  return;
 }
 getParameter("arg","file",&filename);
 in=fopen(filename,"r");
 if(!in){
  T_FREE(buffer);
  return;
 }

 while((!feof(in))  && bFlag){
  fgets(buffer,512-1,in);
  if(buffer[0]=='#')continue;
  if(buffer[strlen(buffer)-1]=='\n')buffer[strlen(buffer)-1]='\0';
  Class=NULL;
  Name=NULL;
  Value=NULL;

  posLeft=buffer;
  /* suche */
  posRight=strstr(posLeft,"::");
  if((!posRight) || (posLeft==posRight)){
   Class=NULL;
   if(posLeft==posRight)posLeft=posRight+2;
  } else {
   T_MALLOC(Class,posRight-posLeft+1,char,bFlag)
   if(!bFlag)break;
   strncpy(Class,posLeft,posRight-posLeft);
   Class[posRight-posLeft]=0;
   posLeft=posRight+2;
  }
  posRight=strchr(posLeft,'=');
  if(posLeft!=posRight){
   if(!posRight){
    T_STRDUP(Name,posLeft,bFlag)
    T_STRDUP_STATIC(Value,"true",bFlag)
    if(!bFlag)break;
   } else {
    T_MALLOC(Name,posRight-posLeft+1,char,bFlag)
    if(!bFlag)break;
    strncpy(Name,posLeft,posRight-posLeft);
    Name[posRight-posLeft]=0;
    posLeft=posRight+1;
    if(strlen(posLeft)>0){
     T_STRDUP(Value,posLeft,bFlag)
     if(!bFlag)break;
    }
   }
  }
  if(Name && Value){
   if(!setParameter(Class,Name,Value)){
    if(Flag)fprintf(stderr,"ignoring unknown option : %s\n",buffer);
   }
  } else {
   if(Flag)fprintf(stderr,"ignoring illegal option : %s\n",buffer);
  }
  T_FREE(Class);
  T_FREE(Name);
  T_FREE(Value);
 }
 T_FREE(Class);
 T_FREE(Name);
 T_FREE(Value);
 T_FREE(buffer);
 fclose(in);
}
 
// parst die von der main() Funktion stammenden Kommandozeilenparameter
/*!
 \internal
*/
void T_CArg::parseArgsInternal(int argc,char *argv[],bool Flag)
{
 char *Class,*Name,*Value;
 char *posLeft,*posRight;
 int i;
 bool bFlag=TRUE;

 for(i=0;i<argc;++i){
  if(argv[i][0]=='-' || argv[i][0]=='+'){
   if(argv[i][1]=='-' && argv[i][2]=='\0'){
    i=argc;
   } else {
   Class=NULL;
   Name=NULL;
   Value=NULL;

   /* abschneiden der '-' chars */
   posLeft=argv[i];
   while(posLeft[0]=='-' || posLeft[0]=='+')posLeft++;
   /* suche */
   posRight=strstr(posLeft,"::");
   if((!posRight) || (posLeft==posRight)){
    Class=NULL;
    if(posLeft==posRight)posLeft=posRight+2;
   } else {
    T_MALLOC(Class,posRight-posLeft+1,char,bFlag)
    if(!bFlag)break;
    strncpy(Class,posLeft,posRight-posLeft);
    Class[posRight-posLeft]=0;
    posLeft=posRight+2;
   }
   posRight=strchr(posLeft,'=');
   if(posLeft!=posRight){
    if(!posRight){
     T_STRDUP(Name,posLeft,bFlag)
     if(argv[i][0]=='-'){
      T_STRDUP_STATIC(Value,"false",bFlag)
     } else {
      T_STRDUP_STATIC(Value,"true",bFlag)
     }
     if(!bFlag)break;
    } else {
     T_MALLOC(Name,posRight-posLeft+1,char,bFlag)
     if(!bFlag)break;
     strncpy(Name,posLeft,posRight-posLeft);
     Name[posRight-posLeft]=0;
     posLeft=posRight+1;
     if(strlen(posLeft)>0){
      T_STRDUP(Value,posLeft,bFlag)
      if(!bFlag)break;
     }
    }
   }
   if(Name && Value){
    if(!setParameter(Class,Name,Value)){
     if(Flag)fprintf(stderr,"ignoring unknown option : %s\n",argv[i]);
    }
   } else {
    if(Flag)fprintf(stderr,"ignoring illegal option : %s\n",argv[i]);
   }
   T_FREE(Class);
   T_FREE(Name);
   T_FREE(Value);
  }}
 }
}
 
// gibt die Beschreibung aller Parameter auf die Standartausgabe
/*
*/
void T_CArg::describe(void)
{
 T_CContainerIterator I;

 for(I.begin(m_CCAI);I.isValid();++I){
  ((T_CArgItem *)I.get())->describe(stdout);
 }
}

// schreibt alle Parameter in die "uber Kommandozeile angegebene Datei
/*
*/
bool T_CArg::dumpArgs(void)
{
 char *filename;
 FILE *out;
 if(isValidParameter("arg","dump")){
  getParameter("arg","dump",&filename);
  if(!strcmp(filename,"-"))dumpArgs(stdout);
  else {
   out=fopen(filename,"w");
   dumpArgs(out);
   fclose(out);
  }
  return TRUE;
 }
 return FALSE;
}

// schreibt alle Parameter in einen Stream
/*
*/
void T_CArg::dumpArgs(FILE *out)
{
 T_CContainerIterator I;
 for(I.begin(m_CCAI);I.isValid();++I){
  ((T_CArgItem *)I.get())->dump(out);
 }
}

// liefert den Wert eines Parameters als bool
bool T_CArg::isValidParameter(const char *Class,const char *Name)
{
 bool buffer;
 getParameter(Class,Name,&buffer);
 return buffer;
}

// setzt einen Parameter auf einen Wert
/*!
 Wird keine Klasse angegeben, wird zun"achst versucht, ob ein Parameter
 der Defaultklasse einen passenden Namen hat. Ist dies nicht der Fall
 wird nochmals unter allen Klassen nach einem passenden Parameter gesucht.
 Tretten bei der direkten Benennung der Klasse, der Defaultklasse oder
 dem Versuch Parameter einer beliebigen Klasse zu verwenden mehrere passende
 Variablen auf, so wird die Routine ohne den Wert gesetzt zu haben
 mit einem FALSE beendet.
*/
bool T_CArg::setParameter(const char *Class,const char *Name,const char *Value)
{
 T_CContainerIterator I;
 T_CArgItem *pCAI,*pCAI_Ok;
 int Counter;

 pCAI_Ok=NULL;

 if(Class){
  Counter=0;
  for(I.begin(m_CCAI);I.isValid();++I){
   pCAI=(T_CArgItem *)I.get();
   if(pCAI->isOption(Class,Name)){
    pCAI_Ok=pCAI;
    ++Counter;
   }
  }
  if(Counter==1){
   return pCAI_Ok->set(Value);
  } else {
   return FALSE;
  }
 } else {
  Class=m_DefaultClass;
  Counter=0;
  for(I.begin(m_CCAI);I.isValid();++I){
   pCAI=(T_CArgItem *)I.get();
   if(pCAI->isOption(Class,Name)){
    pCAI_Ok=pCAI;
    ++Counter;
   }
  }
  if(Counter>1){
   return FALSE;
  } 
  if(Counter==1){
   return pCAI_Ok->set(Value);
  }
  if(Counter==0){
   Counter=0;
   for(I.begin(m_CCAI);I.isValid();++I){
    pCAI=(T_CArgItem *)I.get();
    if(pCAI->hasName(Name)){
     pCAI_Ok=pCAI;
     ++Counter;
    }
   }
   if(Counter==1){
    return pCAI_Ok->set(Value);
   } else {
    return FALSE;
   } 
  }
 }
 return FALSE;
}

// liefert den Wert eines Parameters als String
/*
*/
void T_CArg::getParameter(const char *Class,const char *Name,char **pValue)
{
 T_CContainerIterator I;
 T_CArgItem *pCAI;

 for(I.begin(m_CCAI);I.isValid();++I){
  pCAI=(T_CArgItem *)I.get();
  if(pCAI->isOption(Class,Name)){
   *pValue=pCAI->get();
   return;
  }
 }
 *pValue=NULL;
}

// liefert den Wert eines Parameters als T_Real32_t
/*
*/
void T_CArg::getParameter(const char *Class,const char *Name,T_Real32_t *pValue)
{
 T_CContainerIterator I;
 T_CArgItem *pCAI;

 for(I.begin(m_CCAI);I.isValid();++I){
  pCAI=(T_CArgItem *)I.get();
  if(pCAI->isOption(Class,Name)){
   sscanf(pCAI->get(),"%"T_Real32_SCANF,pValue);
   return;
  }
 }
 *pValue=T_Real32_ZERO;
}

// liefert den Wert eines Parameters als T_Real64_t
/*
*/
void T_CArg::getParameter(const char *Class,const char *Name,T_Real64_t *pValue)
{
 T_CContainerIterator I;
 T_CArgItem *pCAI;

 for(I.begin(m_CCAI);I.isValid();++I){
  pCAI=(T_CArgItem *)I.get();
  if(pCAI->isOption(Class,Name)){
   sscanf(pCAI->get(),"%"T_Real64_SCANF,pValue);
   return;
  }
 }
 *pValue=T_Real64_ZERO;
}

// liefert den Wert eines Parameters als T_Int32_t
/*
*/
void T_CArg::getParameter(const char *Class,const char *Name,T_Int32_t *pValue)
{
 T_CContainerIterator I;
 T_CArgItem *pCAI;

 for(I.begin(m_CCAI);I.isValid();++I){
  pCAI=(T_CArgItem *)I.get();
  if(pCAI->isOption(Class,Name)){
   sscanf(pCAI->get(),"%"T_Int32_SCANF,pValue);
   return;
  }
 }
 *pValue=T_Int32_ZERO;
}

// liefert den Wert eines Parameters als T_Int64_t
/*
*/
void T_CArg::getParameter(const char *Class,const char *Name,T_Int64_t *pValue)
{
 T_CContainerIterator I;
 T_CArgItem *pCAI;

 for(I.begin(m_CCAI);I.isValid();++I){
  pCAI=(T_CArgItem *)I.get();
  if(pCAI->isOption(Class,Name)){
   sscanf(pCAI->get(),"%"T_Int64_SCANF,pValue);
   return;
  }
 }
 *pValue=T_Int64_ZERO;
}

// liefert den Wert eines Parameters als bool
/*
*/
void T_CArg::getParameter(const char *Class,const char *Name,bool *pValue)
{
 T_CContainerIterator I;
 T_CArgItem *pCAI;
 char *Value;

 for(I.begin(m_CCAI);I.isValid();++I){
  pCAI=(T_CArgItem *)I.get();
  if(pCAI->isOption(Class,Name)){
   Value=pCAI->get();
   if(!Value){
    *pValue=FALSE;
    return;
   }
   if((!strcmp(Value,"false")) || (!strcmp(Value,"False")) || (!strcmp(Value,"FALSE"))){
    *pValue=FALSE;
    return;
   }
   *pValue=TRUE;
   return;
  }
 }
 *pValue=FALSE;
 return;
}

// liefert den Wert eines Parameters Minimum..Maximum als T_Real64_t
/*
*/
void T_CArg::getParameter(const char *Class,const char *Name,T_Real64_t *pMin,T_Real64_t *pMax)
{
 T_CContainerIterator I;
 T_CArgItem *pCAI;
 char *pRegion,*pParameter;
 int Skip;

 for(I.begin(m_CCAI);I.isValid();++I){
  pCAI=(T_CArgItem *)I.get();
  if(pCAI->isOption(Class,Name)){
   pRegion=pCAI->get();
   sscanf(pRegion,"%"T_Real64_SCANF,pMin);
   pParameter=strstr(pRegion,"..");
   Skip=2;
   if(!pParameter){
    pParameter=strstr(pRegion,":");
    Skip=1;
   }
   if(!pParameter){
    pParameter=strstr(pRegion,",");
    Skip=1;
   }
   if(pParameter){
    sscanf(pParameter+Skip,"%"T_Real64_SCANF,pMax);
   } else {
    sscanf(pRegion,"%"T_Real64_SCANF,pMax);
   }
   if(*pMin>*pMax){
    T_Real64_t tmp;
    tmp=*pMin;
    *pMin=*pMax;
    *pMax=tmp;
   }
   return;
  }
 }
 *pMin=T_Real64_ZERO;
 *pMax=T_Real64_ZERO;
 return;
}

// liefert den Wert eines Parameters Minimum..Maximum als T_Int32_t
/*
*/
void T_CArg::getParameter(const char *Class,const char *Name,T_Int32_t *pMin,T_Int32_t *pMax)
{
 T_CContainerIterator I;
 T_CArgItem *pCAI;
 char *pRegion,*pParameter;
 int Skip;

 for(I.begin(m_CCAI);I.isValid();++I){
  pCAI=(T_CArgItem *)I.get();
  if(pCAI->isOption(Class,Name)){
   pRegion=pCAI->get();
   sscanf(pRegion,"%"T_Int32_SCANF,pMin);
   pParameter=strstr(pRegion,"..");
   Skip=2;
   if(!pParameter){
    pParameter=strstr(pRegion,":");
    Skip=1;
   }
   if(!pParameter){
    pParameter=strstr(pRegion,",");
    Skip=1;
   }
   if(pParameter){
    sscanf(pParameter+Skip,"%"T_Int32_SCANF,pMax);
   } else {
    sscanf(pRegion,"%"T_Int32_SCANF,pMax);
   }
   if(*pMin>*pMax){
    T_Int32_t tmp;
    tmp=*pMin;
    *pMin=*pMax;
    *pMax=tmp;
   }
   return;
  }
 }
 *pMin=T_Int32_ZERO;
 *pMax=T_Int32_ZERO;
 return;
}

// liefert den ersten Wert eines Parameterbereiches als T_Int64_t
/*
*/
void T_CArg::getParameter(const char *Class,const char *Name,T_Int64_t *pMin,T_Int64_t *pMax)
{
 T_CContainerIterator I;
 T_CArgItem *pCAI;
 char *pRegion,*pParameter;
 int Skip;

 for(I.begin(m_CCAI);I.isValid();++I){
  pCAI=(T_CArgItem *)I.get();
  if(pCAI->isOption(Class,Name)){
   pRegion=pCAI->get();
   sscanf(pRegion,"%"T_Int64_SCANF,pMin);
   pParameter=strstr(pRegion,"..");
   Skip=2;
   if(!pParameter){
    pParameter=strstr(pRegion,":");
    Skip=1;
   }
   if(!pParameter){
    pParameter=strstr(pRegion,",");
    Skip=1;
   }
   if(pParameter){
    sscanf(pParameter+Skip,"%"T_Int64_SCANF,pMax);
   } else {
    sscanf(pRegion,"%"T_Int64_SCANF,pMax);
   }
   if(*pMin>*pMax){
    T_Int64_t tmp;
    tmp=*pMin;
    *pMin=*pMax;
    *pMax=tmp;
   }
   return;
  }
 }
 *pMin=T_Int64_ZERO;
 *pMax=T_Int64_ZERO;
 return;
}
