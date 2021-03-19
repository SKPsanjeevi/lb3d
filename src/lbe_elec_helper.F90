#include "lbe.h"

!> Helper routines for ELEC.
module lbe_elec_helper_module
#ifdef ELEC

  use lbe_elec_globals_module
  use lbe_helper_module, only: is_colloid, is_fluid, is_wall, is_rock, is_surface, every_n_time_steps, unixtime
  use lbe_log_module
  use lbe_parallel_module, only: check_allocate
  use mpi

  implicit none

  private

  public error_elec, log_msg_elec_ws, log_msg_elec, log_msg_elec_hdr ! Logging
  public get_local_eps, get_op, get_nb_phi  ! Wrappers
  public msgstr, stdin, stdout, stderr ! From log_module; expose these variables to ELEC as well.
  public is_colloid, is_fluid, is_wall, is_rock, is_surface, every_n_time_steps, unixtime ! From helper_module; expose these variables to ELEC as well.
  public check_allocate ! From lbe_parallel.

  character(len=4), parameter :: elec_log_prefix = "[EL]"

contains

! ----------------------------------------------------------------------------
!
!                                 LOGGING 
!
! ----------------------------------------------------------------------------

!> Log message with ELEC prefix and surrounding whitespace.
subroutine log_msg_elec_ws(msg, forall)
  implicit none

  character(len=*), intent(in)  :: msg    !< Text to be logged.
  logical, optional, intent(in) :: forall !< If supplied and set to true, all ranks will log.

  if ( present(forall) ) then
    call log_msg("", forall)
    call log_msg(elec_log_prefix // " " // trim(msg), forall)
    call log_msg("", forall)
  else
    call log_msg("")
    call log_msg(elec_log_prefix // " " // trim(msg))
    call log_msg("")
  end if

end subroutine log_msg_elec_ws

!> Log message with ELEC prefix.
subroutine log_msg_elec(msg, forall, error)
  implicit none

  character(len=*), intent(in)  :: msg    !< Text to be logged.
  logical, optional, intent(in) :: forall !< If supplied and set to true, all ranks will log.
  logical, optional, intent(in) :: error  !< If supplied and set to true, log an error.

  if ( present(forall) ) then
    if ( present(error) ) then
      call log_msg(elec_log_prefix // " " // trim(msg), forall = forall, error = error)
    else
      call log_msg(elec_log_prefix // " " // trim(msg), forall = forall)
    end if
  else
    if ( present(error) ) then
      call log_msg(elec_log_prefix // " " // trim(msg), error = error)
    else
      call log_msg(elec_log_prefix // " " // trim(msg))
    end if
  end if

end subroutine log_msg_elec

!> Log header message with ELEC prefix.
subroutine log_msg_elec_hdr(msg, forall)
  implicit none

  character(len=*), intent(in)  :: msg    !< Text to be logged.
  logical, optional, intent(in) :: forall !< If supplied and set to true, all ranks will log.

  if ( present(forall) ) then
    call log_msg_hdr(elec_log_prefix // " " // trim(msg), forall)
  else
    call log_msg_hdr(elec_log_prefix // " " // trim(msg))
  end if

end subroutine log_msg_elec_hdr

!> Log error message with ELEC prefix and exit.
subroutine error_elec(msg)
  implicit none

  character(len=*), intent(in)  :: msg    !< Text to be logged.

  integer :: mpierror

  call log_msg_elec("ERROR: "//trim(msg), forall = .true., error = .false.)
  call log_msg_elec("ERROR: "//trim(msg), forall = .true., error = .true. )
  call MPI_Abort(MPI_COMM_WORLD, -1, mpierror)
  stop

end subroutine error_elec

! ----------------------------------------------------------------------------
!
!                                WRAPPERS 
!
! ----------------------------------------------------------------------------

!> When using local_eps, returns local epsilon at a lattice site, depending on the use of fluid_on_elec and single- / multicomponent fluids.
!> This does not allow for the case of local_eps = .false., for performance reasons (another if-statement).
function get_local_eps(N, i, j, k)
  implicit none

  real(kind=rk) :: get_local_eps
  type(lbe_site), intent(in) :: N(1-halo_extent:, 1-halo_extent:, 1-halo_extent:)
  integer, intent(in) :: i, j, k

  if ( fluid_on_elec ) then
    if ( is_rock( N(i,j,k)%rock_state ) ) then
      get_local_eps = eps_wall
    else
#ifdef SINGLEFLUID
      get_local_eps = eps_r
#else
      get_local_eps = eps_avg * ( 1.0_rk - eps_gamma * get_op(N, i, j, k) )
#endif
    end if
  else
    get_local_eps = N(i,j,k)%eps
  end if
end function get_local_eps

!> Get local order parameter.
pure function get_op(N, i, j, k)
  implicit none

  real(kind=rk) :: get_op
  type(lbe_site), intent(in) :: N(1-halo_extent:, 1-halo_extent:, 1-halo_extent:)
  integer, intent(in) :: i, j, k

#ifndef SINGLEFLUID
  get_op = ( sum(N(i,j,k)%n_r(:)*g) - sum(N(i,j,k)%n_b(:)*g) ) / ( sum(N(i,j,k)%n_r(:)*g) + sum(N(i,j,k)%n_b(:)*g) )
#else
  get_op = 1.0_rk
#endif

end function get_op

!> Calculate neighbour position, taking into account boundary conditions for phi.
!> Also returns the potential on the neighbour site.
function get_nb_phi(N, veldir, i, j, k, ip, jp, kp)
  implicit none

  real(kind=rk) :: get_nb_phi
  type(lbe_site), intent(inout) :: N(1-halo_extent:, 1-halo_extent:, 1-halo_extent:)
  integer, intent(in)  :: veldir, i,  j,  k
  integer, intent(out) :: ip, jp, kp

  ip = i + cx(veldir)
  jp = j + cy(veldir)
  kp = k + cz(veldir)

  if ( boundary_phi_id == boundary_phi_neumannx ) then
    ! For Neumann boundary, enforce grad = 0, so take the neighbour to be the site itself.
    if ( ip + ccoords(1)*nx .lt. 1 ) then
      ip = 1
    end if
    if ( ip + ccoords(1)*nx .gt. tnx ) then
      ip = nx
    end if
    get_nb_phi = N(ip,jp,kp)%phi
  else if ( boundary_phi_id == boundary_phi_neumannz ) then
    ! For Neumann boundary, enforce grad = 0, so take the neighbour to be the site itself.
    if ( kp + ccoords(3)*nz .lt. 1 ) then
      kp = 1
    end if
    if ( kp + ccoords(3)*nz .gt. tnz ) then
      kp = nz
    end if
    get_nb_phi = N(ip,jp,kp)%phi
  else if ( boundary_phi_id == boundary_phi_dropz ) then
    get_nb_phi = N(ip,jp,kp)%phi
    ! For potential drop, shift the potential by phi_dropz.
    if ( kp + ccoords(3)*nz .lt. 1 ) then
      get_nb_phi = N(ip,jp,kp)%phi + phi_dropz
    end if
    if ( kp + ccoords(3)*nz .gt. tnz ) then
      get_nb_phi = N(ip,jp,kp)%phi - phi_dropz
    end if
  else
    get_nb_phi = N(ip,jp,kp)%phi
  end if

end function get_nb_phi

#endif
end module lbe_elec_helper_module

