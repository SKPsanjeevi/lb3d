#ifndef T_DEBUG_H
#define T_DEBUG_H 

//DGKEEP #ifdef DEBUG
#include <stdio.h>
/*!
 Dient dem Debuggen von Programmen. 
 Man f"ugt die Buchstaben DG an Stellen des Sources ein, bei
 denen man sich nicht sicher ist ob das Programm diese erreicht.
 Bei Programmablauf wird dann chronologisch ausgegeben welche
 Stellen im Source wann durchlaufen werden.

*/
#ifdef __GNUC__

#define DGP(VAR) fprintf(stderr,VAR);
#define DG {fflush(stdout);fprintf(stderr,"%s +%i  %s\n",__FILE__,__LINE__,__PRETTY_FUNCTION__);fflush(stderr);}
#define DGASM asm("nop;nop");
//asm("bis $31,$31,$31;bis $31,$31,$31");
#else
#define DG
#define DGASM 
#endif
/*
*/
//extern void DGB(char *file,int line,char *function);
//define DGB DGB(__FILE__,__LINE__,__PRETTY_FUNCTION__);


//DGKEEP #else

/*!
 Dient dem Debuggen von Programmen. 
 Man f"ugt die Buchstaben DG an Stellen des Sources ein, bei
 denen man sich nicht sicher ist ob das Programm diese erreicht.
 Bei Programmablauf wird dann chronologisch ausgegeben welche
 Stellen im Source wann durchlaufen werden.
*/
//DGKEEP #define DG
//DGKEEP #define DGASM 
//DGKEEP //define DGB
//DGKEEP #endif

#endif

