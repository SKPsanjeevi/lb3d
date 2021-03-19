
// M : malloc
// S : strdup
// F : free
// N : new
// D : delete
// E : found an error 

// muessen im Programm auf sinnvolle Werte gesetzt werden
extern FILE *g_stdmemleak;

#if defined(T_ENABLE_MACROSMEMLEAK_C)
#undef T_MALLOC
#undef T_FREE
#undef T_STRDUP
#undef T_STRDUP_STATIC
#undef T_REGISTER

/*!
 Ersatz f"ur die Funktion malloc
*/
#define T_MALLOC(POINTER,NUMBER,TYPE,BOOL) {\
 POINTER=(TYPE*)malloc((NUMBER)*sizeof(TYPE));\
 if(!POINTER){\
  BOOL=FALSE;\
  fprintf(g_stdmemleak,"E malloc returns NULL : "#TYPE"[%i] %s:%i\n",NUMBER,__FILE__,__LINE__);\
  fflush(g_stdmemleak);\
 } else {\
  fprintf(g_stdmemleak,"M %p %i "#TYPE"[%i] %s:%i\n",POINTER,(NUMBER)*sizeof(TYPE),NUMBER,__FILE__,__LINE__);\
  fflush(g_stdmemleak);\
 }\
}

/*!
 Ersatz f"ur die Funktion strdup
*/
#define T_STRDUP(STRING,CONSTSTRING,BOOL) {\
 if(CONSTSTRING){\
  STRING=strdup(CONSTSTRING);\
  if(!STRING){\
   BOOL=FALSE;\
   fprintf(g_stdmemleak,"E strdup returns NULL char[%i] %s:%i\n",strlen(CONSTSTRING),__FILE__,__LINE__);\
   fflush(g_stdmemleak);\
  } else {\
   int i=strlen(STRING);\
   fprintf(g_stdmemleak,"S %p %i char[%i] %s:%i\n",STRING,i,i,__FILE__,__LINE__);\
   fflush(g_stdmemleak);\
  }\
 } else STRING=NULL;\
}
 
/*!
 Ersatz f"ur die Funktion strdup
*/
#define T_STRDUP_STATIC(STRING,CONSTSTRING,BOOL) {\
 STRING=strdup(CONSTSTRING);\
 if(!STRING){\
  BOOL=FALSE;\
  int i=strlen(CONSTSTRING);\
  fprintf(g_stdmemleak,"E new returns NULL : char[%i] %s:%i\n",i,__FILE__,__LINE__);\
  fflush(g_stdmemleak);\
 } else {\
  int i=strlen(STRING);\
  fprintf(g_stdmemleak,"S %p %i char[%i] %s:%i\n",STRING,i,i,__FILE__,__LINE__);\
  fflush(g_stdmemleak);\
 }\
}

/*!
 generiert ein zweidimensionales Array (Pointer[row][column])
*/
#define T_FREE(POINTER) {\
 if(POINTER){\
  fprintf(g_stdmemleak,"F %p %s:%i\n",POINTER,__FILE__,__LINE__);\
  fflush(g_stdmemleak);\
  free(POINTER);\
  POINTER=NULL;\
 }\
}

/*!
 registriert neu allozierten Speicher, z.B. der durch asprintf oder
 Routinen der libxml zurueckgegeben werden.
*/
#define T_REGISTER(POINTER,NUMBER,TYPE) {\
 if(!POINTER){\
  fprintf(g_stdmemleak,"E function returns NULL : "#TYPE"[%i] %s:%i\n",NUMBER,__FILE__,__LINE__);\
  fflush(g_stdmemleak);\
 } else {\
  fprintf(g_stdmemleak,"R %p %i "#TYPE"[%i] %s:%i\n",POINTER,(NUMBER)*sizeof(TYPE),NUMBER,__FILE__,__LINE__);\
  fflush(g_stdmemleak);\
 }\
}

#endif

#if defined(T_ENABLE_MACROSMEMLEAK_CC)
#undef T_NEW
#undef T_DELETE

/*!
 Ersatz f"ur den Operator new
*/
#define T_NEW(POINTER,CLASS,FUNC) {\
 (POINTER)=new CLASS FUNC;\
 if(POINTER)fprintf(g_stdmemleak,"N %p %i "#CLASS" %s:%i\n",POINTER,sizeof(CLASS),__FILE__,__LINE__);\
 else fprintf(g_stdmemleak,"E new returns NULL : "#CLASS" %s:%i\n",__FILE__,__LINE__);\
 fflush(g_stdmemleak);\
}

/*!
 Ersatz f"ur den Operator delete,
 es wird sichergestellt, da"s nur gesetzte Pointer
 freigegeben werden, und da"s nach der Freigabe
 der Pointer auf NULL gesetzt wird.
*/
#define T_DELETE(POINTER) {\
 if(POINTER){\
  fprintf(g_stdmemleak,"D %p %s:%i\n",POINTER,__FILE__,__LINE__);\
  fflush(g_stdmemleak);\
  delete POINTER;\
  POINTER=NULL;\
 }\
}
#endif
