#include "lbe.h"

!> Defines the datatypes \c lbe_site, and \c radial_model.
!>
!> \warning If you change the lbe_site type, then you
!> MUST ALWAYS change the \c MPI_Type definitions in
!> \p lbe_parallel.f90, subroutine \c lbe_parallel_init().
module lbe_types_module

  use lbe_globals_module, only: rk

  implicit none
  public

  private rk

  type lbe_site
     sequence	! Contiguous in memory, please..
     real*8, dimension(19) :: n_r  ! No of red particles
        
#ifndef SINGLEFLUID
     real*8, dimension(19) :: n_b  ! No of blue particles
#endif
#ifndef NOSURFACTANT
     real*8, dimension(19) :: n_s  ! no of surf particles
     real*8, dimension(3)  :: d    ! Dipole moment vector
#endif

#ifdef INTERPOLATEDBB
     real*8 :: rho_0               ! Initial density
     !distance ratio between fluid node and the wall along the velocity directions
     real*8, dimension(19) :: delta_q  
#endif

#ifdef ELEC
     real(kind=rk) :: rho_p !< Scalar field for positive charges (plus)
     real(kind=rk) :: rho_m !< Scalar field for negative charges (minus)
     real(kind=rk) :: phi   !< Scalar field for electric potential
     real(kind=rk) :: eps   !< Local dielectric
     real(kind=rk), dimension(3) :: E
     real(kind=rk), dimension(3) :: elec_force !< electrokinetic force
#endif

     real*8	:: rock_state  ! (0,1) for rockp.
#ifdef DIST
     real*8	:: abst
#endif
#ifdef LOCALBC
     real*8	:: local_acccoef
#endif
#ifdef VARTAU
     real*8	:: taupos_r
#ifndef SINGLEFLUID
     real*8	:: taupos_b
#endif
#ifndef NOSURFACTANT
     real*8	:: taupos_s
#endif

#endif
     real*8	:: rock_colour
     real*8 :: rock_colour_r
     real*8 :: rock_colour_b
!#ifdef COMMON_VEL_FIX
#ifndef OLD_VEL
 real*8, dimension(19) :: n_r_pre  ! No of red particles pre collision
#ifndef SINGLEFLUID 
 real*8, dimension(19) :: n_b_pre  ! No of blue particles pre collision
#endif

#ifndef NOSURFACTANT
 real*8, dimension(19) :: n_s_pre  ! No of surfactant particles per collision
#endif

#endif


end type lbe_site

type radial_model
   real*8 :: n_r, n_b, n_s
   integer :: dip
end type radial_model

#ifdef LOCALBC
type bc_check
   sequence	! Contiguous in memory, please..
   integer, dimension(19) :: n
end type bc_check
#endif

end module lbe_types_module
