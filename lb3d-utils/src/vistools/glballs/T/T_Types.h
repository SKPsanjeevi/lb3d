#include "config.h"

#if !defined(T_Types_H)
#define T_Types_H

#include <float.h> // DBL_MAX ..
#if !defined(__STDC_LIMIT_MACROS)
#define __STDC_LIMIT_MACROS
#endif
//#include <stdint.h> // INT32_MAX ..


/*!
 Aussage ist wahr
*/
#define TRUE 1
/*!
 Aussage ist falsch
*/
#define FALSE 0

#if !defined(NULL)
/*!
 Null-pointer
*/
#define NULL 0
#endif

/*
 Datentyp f"ur 16 Bit unsigned Integer
 Definiert den Wert Null
*/
#define T_UInt16_ZERO 0
#define T_UInt16_MIN 0
#define T_UInt16_MAX (65535)
#if SIZEOF_UNSIGNED_SHORT_INT == 2
 typedef unsigned short int T_UInt16_t;
#define T_UInt16_SCANF "hi"
#define T_UInt16_PRINTF "hi"
#define T_UInt16_MPI MPI_UNSIGNED_SHORT
#else
#if SIZEOF_UNSIGNED_INT == 2
 typedef unsigned int T_UInt16_t;
#define T_UInt16_SCANF "i"
#define T_UInt16_PRINTF "i"
#define T_UInt16_MPI MPI_UNSIGNED
#else
//#warning has no 16 bit unsigned integer type
 typedef unsigned short int T_UInt16_t;
#define T_UInt16_SCANF "hi"
#define T_UInt16_PRINTF "hi"
#define T_UInt16_MPI MPI_UNSIGNED_SHORT
#endif
#endif

/*
 Datentyp f"ur 32 Bit Integer
 Definiert den Wert Null
*/
#define T_Int32_MIN (-2147483647-1)
#define T_Int32_MAX (2147483647)
#if SIZEOF_INT == 4
 typedef int T_Int32_t;
#define T_Int32_ZERO 0
#define T_Int32_SCANF "i"
#define T_Int32_PRINTF "i"
#define T_Int32_MPI MPI_INT
#else
#if SIZEOF_SHORT_INT == 4
 typedef short int T_Int32_t;
#define T_Int32_ZERO 0
#define T_Int32_SCANF "hi"
#define T_Int32_PRINTF "hi"
#define T_Int32_MPI MPI_SHORT
#else
#if SIZEOF_LONG_INT == 4
 typedef long int T_Int32_t;
#define T_Int32_ZERO 0l
#define T_Int32_SCANF "li"
#define T_Int32_PRINTF "li"
#define T_Int32_MPI MPI_LONG
#else
#error has no 32 bit integer type
#endif
#endif
#endif

/*
 Datentyp f"ur 64 Bit Integer
 Definiert den Wert Null
*/
#if SIZEOF_LONG_INT == 8
#define T_Int64_MIN (-9223372036854775807L-1)
#define T_Int64_MAX (9223372036854775807L)
 typedef long int T_Int64_t;
#define T_Int64_ZERO 0l
#define T_Int64_SCANF "li"
#define T_Int64_PRINTF "li"
#define T_Int64_MPI MPI_LONG
#else
#if SIZEOF_SHORT_INT == 8
#define T_Int64_MIN (-9223372036854775807-1)
#define T_Int64_MAX (9223372036854775807)
 typedef short int T_Int64_t;
#define T_Int64_ZERO 0
#define T_Int64_SCANF "hi"
#define T_Int64_PRINTF "hi"
#define T_Int64_MPI MPI_SHORT
#else
#if SIZEOF_INT == 8
#define T_Int64_MIN (-9223372036854775807-1)
#define T_Int64_MAX (9223372036854775807)
 typedef int T_Int64_t;
#define T_Int64_ZERO 0
#define T_Int64_SCANF "i"
#define T_Int64_PRINTF "i"
#define T_Int64_MPI MPI_INT
#else
#if SIZEOF_LONG_LONG_INT == 8
#define T_Int64_MIN (-9223372036854775807LL-1)
#define T_Int64_MAX (9223372036854775807LL)
 typedef long long int T_Int64_t;
#define T_Int64_ZERO 0ll
#define T_Int64_SCANF "lli"
#define T_Int64_PRINTF "lli"
#define T_Int64_MPI MPI_LONG_LONG_INT
#else
#error has no 64 bit integer type
#endif
#endif
#endif
#endif

/*
 Datentyp f"ur 32 Bit Real
 Definiert den Wert Null
*/
#if SIZEOF_FLOAT == 4
 typedef float T_Real32_t;
#define T_Real32_ZERO 0.
#define T_Real32_MIN -FLT_MAX
#define T_Real32_MAX FLT_MAX
#define T_Real32_SCANF "e"
#define T_Real32_PRINTF "e"
#define T_Real32_MPI MPI_FLOAT
#else
#if SIZEOF_DOUBLE == 4
 typedef double T_Real32_t;
#define T_Real32_ZERO 0.
#define T_Real32_MIN -DBL_MAX
#define T_Real32_MAX DBL_MAX
#define T_Real32_SCANF "le"
#define T_Real32_PRINTF "le"
#define T_Real32_MPI MPI_DOUBLE
#else
#error have no 32 bit real type
#endif
#endif

/*
 Datentyp f"ur 64 Bit Real
 Definiert den Wert Null
*/
#if SIZEOF_FLOAT == 8
 typedef float T_Real64_t;
#define T_Real64_ZERO 0.
#define T_Real64_MIN -FLT_MAX
#define T_Real64_MAX FLT_MAX
#define T_Real64_SCANF "e"
#define T_Real64_PRINTF "e"
#define T_Real64_MPI MPI_FLOAT
#else
#if SIZEOF_DOUBLE == 8
 typedef double T_Real64_t;
#define T_Real64_ZERO 0.
#define T_Real64_MIN -DBL_MAX
#define T_Real64_MAX DBL_MAX
#define T_Real64_SCANF "le"
#define T_Real64_PRINTF "e"
#define T_Real64_MPI MPI_DOUBLE
#else
#error have no 64 bit real type
#endif
#endif

#endif
