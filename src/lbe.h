! Variable type
#define LBE_REAL MPI_REAL8

! IBM_PART implies XDRF and USEHDF
#ifdef IBM_PART
#ifndef USEXDRF
#define USEXDRF
#endif
#ifndef USEHDF
#define USEHDF
#endif
#endif

! IBM_PART always implies SINGLEFLUID, unless IBM_BINARYIBM is defined.
#ifdef IBM_PART
#ifndef IBM_BINARYIBM
#ifndef SINGLEFLUID
#define SINGLEFLUID
#endif
#else
#ifndef NOSURFACTANT
#define NOSURFACTANT
#endif
#endif
#endif

! IBM_PART and VARTAU imply IBM_INDEXFIELD
#if defined IBM_PART && defined VARTAU
#ifndef IBM_INDEXFIELD
#define IBM_INDEXFIELD
#endif
#endif

! IBM_PART and IBM_BINARYIBM imply IBM_INDEXFIELD
#if defined IBM_PART && defined IBM_BINARYIBM
#ifndef IBM_INDEXFIELD
#define IBM_INDEXFIELD
#endif
#endif

! Single fluid implies no surfactant
#ifdef SINGLEFLUID
#ifndef NOSURFACTANT
#define NOSURFACTANT
#endif
#endif

! Capture forbidden interpolated bounceback BC combinations
#if !defined(SINGLEFLUID) && defined(INTERPOLATEDBB)
#error "Interpolated bounceback (Bouzidi) BC is available only for single fluid."
#endif


! Capture forbidden IBM combinations
#if defined IBM_BINARYIBM && defined SINGLEFLUID
#error "IBM_BINARYIBM and SINGLEFLUID are not compatible."
#endif
#if defined IBM_BINARYIBM && defined NOFLUIDDYNAMICS
#error "IBM_BINARYIBM and NOFLUIDDYNAMICS are not compatible."
#endif
#if defined IBM_PART && defined MD
#error "IBM_PART and MD are not compatible."
#endif

! With the DEBUG_MPI flag set, show the debug output.
! Otherwise, replace the macros with nothing.
#ifdef DEBUG_MPI
#define DEBUG_MPI_MSG(msg) call log_mpi(msg)
#define DEBUG_CHECKMPI(ierror,msg) call checkmpi(ierror,msg)
#else
#define DEBUG_MPI_MSG(msg)
#define DEBUG_CHECKMPI(ierror,msg)
#endif

! Debug macro
#ifdef DEBUG
#define DEBUG_MSG(msg) call log_msg(msg)
#define DEBUG_MSG_ALL(msg) call log_msg(msg,.true.)
#define DEBUG_BARRIER(msg) call log_msg(msg,.true.) ; call MPI_Barrier(comm_cart) ; call log_msg("POST BARRIER", .true.)
#else
#define DEBUG_MSG(msg)
#define DEBUG_MSG_ALL(msg)
#define DEBUG_BARRIER(msg)
#endif

!#ifdef OLD_VEL
!call log_msg("Warning: You are using old version for velocity output!!! Shan-Chen force is not corrected added into velocity calculations.")
!#endif


