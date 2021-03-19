#include "config.h"

#if !defined(T_Macro_H)
#define T_Macro_H

//#define _ISOC99_SOURCE
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <string.h> // strdup

/*!
 Ersatz f"ur die Funktion malloc
*/
#define T_MALLOC(POINTER,NUMBER,TYPE,BOOL) {\
 POINTER=(TYPE*)malloc((NUMBER)*sizeof(TYPE));\
 if(!POINTER)BOOL=FALSE;\
}

/*!
 generiert ein zweidimensionales Array (Pointer[row][column])
*/
#define T_MALLOC2D(POINTER,ROWS,COLUMNS,TYPE,BOOL) {\
 T_MALLOC(POINTER,ROWS,TYPE*,BOOL);\
 if(BOOL)T_MALLOC(POINTER[0],(ROWS)*(COLUMNS),TYPE,BOOL);\
 if(BOOL){\
  for(int malloc2d_i=1;malloc2d_i<(ROWS);++malloc2d_i){\
   POINTER[malloc2d_i]=POINTER[malloc2d_i-1]+(COLUMNS);\
  }\
 }\
 if(!BOOL)T_FREE(POINTER);\
}

/*!
 generiert ein dreidimensionales Array (Pointer[depth][height][width])
 die Belegung ist per Definition (Pointer[z][y][x])
*/
#define T_MALLOC3D(POINTER,DEPTH,HEIGHT,WIDTH,TYPE,BOOL) {\
 TYPE *malloc3d_p=NULL;\
 T_MALLOC(POINTER,DEPTH,TYPE**,BOOL);\
 if(BOOL){\
  for(int malloc3d_i=0;malloc3d_i<(DEPTH);++malloc3d_i){\
   T_MALLOC(POINTER[malloc3d_i],(HEIGHT),TYPE*,BOOL);\
 }}\
 if(BOOL)T_MALLOC(malloc3d_p,(DEPTH)*(HEIGHT)*(WIDTH),TYPE,BOOL);\
 if(BOOL){\
  for(int malloc3d_j=0;malloc3d_j<(DEPTH);++malloc3d_j){\
   for(int malloc3d_i=0;malloc3d_i<(HEIGHT);++malloc3d_i){\
    POINTER[malloc3d_j][malloc3d_i]=malloc3d_p;\
    malloc3d_p+=(WIDTH);\
  }}\
 } else {\
  T_FREE(malloc3d_p)\
 }\
 if(!BOOL){\
  if(POINTER){\
   for(int malloc3d_i=0;malloc3d_i<(DEPTH);++malloc3d_i){\
    T_FREE(POINTER[malloc3d_i]);\
   }\
  } else {\
   T_FREE(POINTER)\
}}}

/*!
 Ersatz f"ur die Funktion strdup
*/
#define T_STRDUP(STRING,CONSTSTRING,BOOL) {\
 if(CONSTSTRING){\
  STRING=strdup(CONSTSTRING);\
  if(!STRING)BOOL=FALSE;\
 } else STRING=NULL;\
}

/*!
 Ersatz f"ur die Funktion strdup
*/
#define T_STRDUP_STATIC(STRING,CONSTSTRING,BOOL) {\
 STRING=strdup(CONSTSTRING);\
 if(!STRING)BOOL=FALSE;\
}

/*!
 Ersatz f"ur die Funktion free
*/
#define T_FREE(POINTER) {if(POINTER){free(POINTER);POINTER=NULL;}}

/*!
 gibt ein mit T_MALLOC2D erzeugtes zweidimensionales Array frei
*/
#define T_FREE2D(POINTER) {\
 if(POINTER){\
  T_FREE(POINTER[0])\
  T_FREE(POINTER)\
 }\
}

/*!
 gibt ein mit T_MALLOC3D erzeugtes dreidimensionales Array frei
*/
#define T_FREE3D(POINTER,DEPTH) {\
 if(POINTER){\
  T_FREE(POINTER[0][0])\
  for(int free3d_i=0;free3d_i<(DEPTH);++free3d_i){\
   T_FREE(POINTER[free3d_i])\
  }\
  T_FREE(POINTER)\
 }\
}

/*!
 Ersatz f"ur die Funktion memset
*/
#define T_ZERO(POINTER,NUMBER,TYPE) {\
 memset(POINTER,0,(NUMBER)*sizeof(TYPE));\
}

/*!
 Ersatz f"ur den Operator new
*/
#define T_NEW(POINTER,CLASS,FUNC) {(POINTER)=new CLASS FUNC;}

/*!
 Ersatz f"ur den Operator delete,
 es wird sichergestellt, da"s nur gesetzte Pointer
 freigegeben werden, und da"s nach der Freigabe
 der Pointer auf NULL gesetzt wird.
*/
#define T_DELETE(POINTER) {if(POINTER){delete POINTER;POINTER=NULL;}}

/*!
 registriert neu allozierten Speicher, z.B. der durch asprintf oder
 Routinen der libxml zurueckgegeben werden.
*/
#define T_REGISTER(POINTER,NUMBER,TYPE) {}

/*!
 dump eines eindimensionalen Arrays (Pointer[length])
*/
#define T_DUMP(PCD,POINTER,LENGTH,TYPE,BOOL) {\
 BOOL&=PCD->Memory((T_MemoryPointer_t *)&POINTER,LENGTH*sizeof(TYPE));\
}

/*!
 dump eines zweidimensionales Arrays (Pointer[row][column])
*/
#define T_DUMP2D(PCD,POINTER,ROWS,COLUMNS,TYPE,BOOL) {\
 if(BOOL){\
  BOOL&=PCD->Memory((T_MemoryPointer_t*)&POINTER,(ROWS)*sizeof(TYPE*));\
  BOOL&=PCD->PointerArray((T_MemoryPointer_t*)POINTER,(ROWS));\
 }\
 if(BOOL)BOOL&=PCD->memory((T_MemoryPointer_t)POINTER[0],(ROWS)*(COLUMNS)*sizeof(TYPE));\
}

/*!
 dump eines dreidimensionales Arrays (Pointer[depth][height][width])
*/
#define T_DUMP3D(PCD,POINTER,DEPTH,HEIGHT,WIDTH,TYPE,BOOL) {\
 if(BOOL){\
  BOOL&=PCD->Memory((T_MemoryPointer_t*)&POINTER,(DEPTH)*sizeof(TYPE**));\
  BOOL&=PCD->PointerArray((T_MemoryPointer_t *)POINTER,(DEPTH));\
 }\
 if(BOOL){\
  for(int free3d_i=0;free3d_i<(DEPTH);++free3d_i){\
   BOOL&=PCD->memory((T_MemoryPointer_t)POINTER[free3d_i],(HEIGHT)*sizeof(TYPE*));\
   BOOL&=PCD->PointerArray((T_MemoryPointer_t *)POINTER[free3d_i],(HEIGHT));\
 }}\
 if(BOOL)BOOL&=PCD->memory((T_MemoryPointer_t)POINTER[0][0],(DEPTH)*(HEIGHT)*(WIDTH)*sizeof(TYPE));\
}

//! gibt den groesseren Wert von A und B zurueck
#define T_MAXIMUM(A,B) ((A) > (B) ? (A) : (B))
//! gibt den keineren Wert von A und B zurueck
#define T_MINIMUM(A,B) ((A) < (B) ? (A) : (B))
//! beschraenkt den Wertebereich auf MIN bis MAX
#define T_ADJUST_TO_RANGE(VAL,MIN,MAX) {if(VAL<MIN){VAL=MIN;}else {if(VAL>MAX)VAL=MAX;}}

#if defined(T_ENABLE_MACROSMEMLEAK_C) || defined(T_ENABLE_MACROSMEMLEAK_CC)
#include "T_MacrosMemLeak.h"
#endif

#endif
