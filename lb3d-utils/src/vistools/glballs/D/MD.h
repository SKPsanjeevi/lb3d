#include "config.h"

#if !defined(MD_H)
#define MD_H

//#define MD_ENABLE_DISPLACEMENT
//#define MD_ENABLE_PARTICLECLUSTER

#include "T.h"

#if defined(MD_ENABLE_MPI)
#include "mpi.h"
extern MPI_Comm MPIComm;
#if defined(T_ENABLE_INTERCOMM)
extern MPI_Comm MPIInterComm;
#endif
extern char MPISimulationType;
#define MD_MPI_LEADER 0
#endif

typedef T_Int32_t MD_Particle_Counter_t;
#define MD_PARTICLE_COUNTER_ZERO T_Int32_ZERO
#define MD_PARTICLE_COUNTER_SCANF T_Int32_SCANF
#define MD_PARTICLE_COUNTER_PRINTF T_Int32_PRINTF
#define MD_PARTICLE_COUNTER_MPI T_Int32_MPI

typedef T_Int32_t MD_Cell_Counter_t;
#define MD_CELL_COUNTER_ZERO T_Int32_ZERO
#define MD_CELL_COUNTER_SCANF T_Int32_SCANF
#define MD_CELL_COUNTER_PRINTF T_Int32_PRINTF
#define MD_CELL_COUNTER_MPI T_Int32_MPI

#if defined(MD_ENABLE_CONTACT)
typedef T_Int32_t MD_Contact_Counter_t;
#define MD_CONTACT_COUNTER_ZERO T_Int32_ZERO
#define MD_CONTACT_COUNTER_SCANF T_Int32_SCANF
#define MD_CONTACT_COUNTER_PRINTF T_Int32_PRINTF
#define MD_CONTACT_COUNTER_MPI T_Int32_MPI
#endif

typedef T_Int32_t MD_UnitId_t;
#define MD_UNITID_ZERO T_Int32_ZERO
#define MD_UNITID_DEFAULT MD_UNITID_ZERO
#define MD_UNITID_SCANF T_Int32_SCANF
#define MD_UNITID_PRINTF T_Int32_PRINTF
#define MD_UNITID_MPI T_Int32_MPI

typedef T_Int32_t MD_Color_t;
#define MD_COLOR_ZERO T_Int32_ZERO
#define MD_COLOR_DEFAULT MD_COLOR_ZERO
#define MD_COLOR_SCANF T_Int32_SCANF
#define MD_COLOR_PRINTF T_Int32_PRINTF
#define MD_COLOR_MPI T_Int32_MPI

#endif
