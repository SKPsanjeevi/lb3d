#include "lbe.h"

!> Output for ELEC part of the code
module lbe_elec_output_module
#ifdef ELEC
  use lbe_elec_globals_module
  use lbe_elec_helper_module
  use lbe_elec_timer_module
  use lbe_parms_module, only: dump_format, lbeversion, lbeflags, lbeplatform
  use lbe_io_helper_module, only: lbe_make_filename_rank
#ifdef USEHDF
  use lbe_io_module, only: dump_scalar
#endif

  implicit none
  private

  public elec_postprocess, elec_dump_all_asc

contains

!> Check when to dump ELEC data, and then make the appropriate calls.
subroutine elec_postprocess(N)
  implicit none

  type(lbe_site), intent(inout) :: N(1-halo_extent:, 1-halo_extent:, 1-halo_extent:)

  if (     ( every_n_time_steps( n_sci_rho_p ) ) &
       .or.( every_n_time_steps( n_sci_rho_m ) ) &
       .or.( every_n_time_steps( n_sci_phi   ) ) &
       .or.( every_n_time_steps( n_sci_eps   ) ) &
       .or.( every_n_time_steps( n_sci_E   ) ) ) then

    call log_msg_elec("Started postprocessing:")

    if ( index(dump_format, 'hdf') .gt. 0 ) then
      call log_msg_elec("Dumping HDF5 data:")
#ifdef USEHDF
      call elec_dump_data(N)
#else
      call error_elec("FATAL ERROR: HDF5 support is switched off but HDF5 output is requested. Aborting...")
#endif
      call log_msg_elec("Finished dumping HDF5 data.")
    else ! not phdf5
      call log_msg_elec("Postprocessing, no HDF5. Not implemented.")
    end if
  end if

  if ( every_n_time_steps( n_sci_elec ) ) then
    call elec_dump_all_asc(N)
  end if

end subroutine elec_postprocess

!> This function will dump all elec-specific fields to disc in ascii format.
!> Will write one file per rank.
subroutine elec_dump_all_asc(N, prefix)
  implicit none

  type(lbe_site), intent(inout) :: N(1-halo_extent:, 1-halo_extent:, 1-halo_extent:)
  character(len=*), intent(in), optional :: prefix

  integer, parameter :: funit = 10

  character(len=256) :: filename
  integer :: nxi, nyi, nzi
  integer :: ti, tj, tk
  integer :: i, j, k
  integer :: imin, imax, jmin, jmax, kmin, kmax

  real(kind=rk) :: eps

  call start_timer(ti_elec_io)

  nxi = size(N,1) - 2*halo_extent
  nyi = size(N,2) - 2*halo_extent
  nzi = size(N,3) - 2*halo_extent

  ! This hack is needed to get output within the equilibration loops.
  ! Prefix will then contain, for example, substep and number of iterations.
  if ( present(prefix) ) then
    call lbe_make_filename_rank(filename, 'elec-'//prefix, '.asc', nt, myrankc)
  else
    call lbe_make_filename_rank(filename, 'elec', '.asc', nt, myrankc)
  end if

  ! Extend sizes when wanting to dump the halo (usually for debugging only).
  if ( dump_elec_halo ) then
    imin = 1 - halo_extent
    imax = nxi + halo_extent
    jmin = 1 - halo_extent
    jmax = nyi + halo_extent
    kmin = 1 - halo_extent
    kmax = nzi + halo_extent
  else
    imin = 1
    imax = nxi
    jmin = 1
    jmax = nyi
    kmin = 1
    kmax = nzi
  end if

  write(msgstr,"('Writing ELEC data to <',A,'>')") trim(filename) ; call log_msg_elec(msgstr)

  open(unit = funit, file = filename)

  ! Write header information
  write(funit, "('# ',A)") trim(lbeversion)
  write(funit, "('# ',A)") trim(lbeflags)
  write(funit, "('# ',A)") trim(lbeplatform)
  write(funit, "('#')")
  write(funit, "('# Raw ELEC data ; t = ',I8.8,X,A)") nt, prefix
  write(funit, "('#')")
  write(funit, "('# myrankc = ',I0)") myrankc
  write(funit, "('# tnx     = ',I0)") tnx
  write(funit, "('# tny     = ',I0)") tny
  write(funit, "('# tnz     = ',I0)") tnz
  write(funit, "('#')")
  write(funit, "(3(A7,X),8(A15,X),A5)") "#?    x", "y", "z", "rho_p", "rho_m", "phi", "eps", "Ex", "Ey", "Ez", "op", "rock"

  ! Write fields
  do k = kmin, kmax
    tk = k + ccoords(3)*nz
    do j = jmin, jmax
      tj = j + ccoords(2)*ny
      do i = imin, imax
        ti = i + ccoords(1)*nx

        if ( local_eps ) then
          eps = get_local_eps(N,i,j,k)
        else
          eps = eps_uniform
        end if

        write(funit,"(3(SP,I7.6,X),S,8(ES15.8,X),F5.2)") &
             ti, tj, tk, N(i,j,k)%rho_p, N(i,j,k)%rho_m, N(i,j,k)%phi, eps, N(i,j,k)%E(:), get_op(N,i,j,k), N(i,j,k)%rock_state
      end do
    end do
  end do

  close(unit = funit)

  call stop_timer(ti_elec_io)

end subroutine elec_dump_all_asc

#ifdef USEHDF
!> This subroutine just calls the specific routines to dump the data to disk
subroutine elec_dump_data(N)
  implicit none

  type(lbe_site), intent(inout) :: N(1-halo_extent:, 1-halo_extent:, 1-halo_extent:)

  call start_timer(ti_elec_io)

  if ( every_n_time_steps( n_sci_rho_p ) ) then
    call log_msg_elec("  Dumping RHO_P ...")
    call dump_rho_p(N)
  end if
  if ( every_n_time_steps( n_sci_rho_m ) ) then
    call log_msg_elec("  Dumping RHO_M ...")
    call dump_rho_m(N)
  end if
  if ( every_n_time_steps( n_sci_phi   ) ) then
    call log_msg_elec("  Dumping PHI ...")
    call dump_phi(N)
  end if
  if ( every_n_time_steps( n_sci_eps   ) ) then
    call log_msg_elec("  Dumping EPS ...")
    call dump_eps(N)
  end if
  if ( every_n_time_steps( n_sci_E     ) ) then
    call log_msg_elec("  Dumping E ...")
    call dump_E(N)
  end if

  call stop_timer(ti_elec_io)

end subroutine elec_dump_data

!> Copy-pasted function to dump rho_p into a file.
subroutine dump_rho_p(N)
  implicit none

  type(lbe_site), intent(inout) :: N(1-halo_extent:, 1-halo_extent:, 1-halo_extent:)

  real(kind=rk), dimension(:, :, :), allocatable :: rho_p
  integer :: nxi, nyi, nzi
  integer :: x, y, z, ierror

  nxi = size(N,1) - 2*halo_extent
  nyi = size(N,2) - 2*halo_extent
  nzi = size(N,3) - 2*halo_extent

  allocate( rho_p(1:nxi, 1:nyi, 1:nzi ), stat = ierror)
  call check_allocate(ierror, "Unable to allocate scalar buffer for rho_p output.")

  do z = 1, nzi
    do y = 1, nyi
      do x = 1, nxi
        rho_p(x,y,z) = N(x,y,z)%rho_p
      end do
    end do
  end do
  call dump_scalar(rho_p, 'rho_p')

  deallocate(rho_p, stat = ierror)
  call check_allocate(ierror, "Unable to deallocate scalar buffer for rho_p output.")

end subroutine dump_rho_p

!> Copy-pasted function to dump rho_m into a file.
subroutine dump_rho_m(N)
  implicit none

  type(lbe_site), intent(inout) :: N(1-halo_extent:, 1-halo_extent:, 1-halo_extent:)

  real(kind=rk), dimension(:, :, :), allocatable :: rho_m
  integer :: nxi, nyi, nzi
  integer :: x, y, z, ierror

  nxi = size(N,1) - 2*halo_extent
  nyi = size(N,2) - 2*halo_extent
  nzi = size(N,3) - 2*halo_extent

  allocate( rho_m(1:nxi, 1:nyi, 1:nzi ), stat = ierror)
  call check_allocate(ierror, "Unable to allocate scalar buffer for rho_m output.")

  do z = 1, nzi
    do y = 1, nyi
      do x = 1, nxi
        rho_m(x,y,z) = N(x,y,z)%rho_m
      end do
    end do
  end do
  call dump_scalar(rho_m, 'rho_m')

  deallocate(rho_m, stat = ierror)
  call check_allocate(ierror, "Unable to deallocate scalar buffer for rho_m output.")

end subroutine dump_rho_m

!> Copy-pasted function to dump phi into a file.
subroutine dump_phi(N)
  implicit none

  type(lbe_site), intent(inout) :: N(1-halo_extent:, 1-halo_extent:, 1-halo_extent:)

  real(kind=rk), dimension(:, :, :), allocatable :: phi
  integer :: nxi, nyi, nzi
  integer :: x, y, z, ierror

  nxi = size(N,1) - 2*halo_extent
  nyi = size(N,2) - 2*halo_extent
  nzi = size(N,3) - 2*halo_extent

  allocate( phi(1:nxi, 1:nyi, 1:nzi ), stat = ierror)
  call check_allocate(ierror, "Unable to allocate scalar buffer for phi output.")

  do z = 1, nzi
    do y = 1, nyi
      do x = 1, nxi
        phi(x,y,z) = N(x,y,z)%phi
      end do
    end do
  end do
  call dump_scalar(phi, 'phi')

  deallocate(phi, stat = ierror)
  call check_allocate(ierror, "Unable to deallocate scalar buffer for phi output.")

end subroutine dump_phi

!> Copy-pasted function to dump eps into a file.
subroutine dump_eps(N)
  implicit none

  type(lbe_site), intent(inout) :: N(1-halo_extent:, 1-halo_extent:, 1-halo_extent:)

  real(kind=rk), dimension(:, :, :), allocatable :: eps
  integer :: nxi, nyi, nzi
  integer :: x, y, z, ierror

  nxi = size(N,1) - 2*halo_extent
  nyi = size(N,2) - 2*halo_extent
  nzi = size(N,3) - 2*halo_extent

  allocate( eps(1:nxi, 1:nyi, 1:nzi ), stat = ierror)
  call check_allocate(ierror, "Unable to allocate scalar buffer for eps output.")

  do z = 1, nzi
    do y = 1, nyi
      do x = 1, nxi

        if ( local_eps ) then
          eps(x,y,z) = get_local_eps(N,x,y,z)
        else
          eps(x,y,x) = eps_uniform
        end if

      end do
    end do
  end do
  call dump_scalar(eps, 'eps')

  deallocate(eps, stat = ierror)
  call check_allocate(ierror, "Unable to deallocate scalar buffer for eps output.")

end subroutine dump_eps

!> Copy-pasted function to dump E into a triplet of files.
subroutine dump_E(N)
  implicit none

  type(lbe_site), intent(inout) :: N(1-halo_extent:, 1-halo_extent:, 1-halo_extent:)

  real(kind=rk), dimension(:, :, :), allocatable :: E
  integer :: nxi, nyi, nzi
  integer :: x, y, z, ierror

  nxi = size(N,1) - 2*halo_extent
  nyi = size(N,2) - 2*halo_extent
  nzi = size(N,3) - 2*halo_extent

  allocate( E(1:nxi, 1:nyi, 1:nzi ), stat = ierror)
  call check_allocate(ierror, "Unable to allocate scalar buffer for E output.")

  do z = 1, nzi
    do y = 1, nyi
      do x = 1, nxi
        E(x,y,z) = N(x,y,z)%E(1)
      end do
    end do
  end do
  call dump_scalar(E, 'Ex')

  do z = 1, nzi
    do y = 1, nyi
      do x = 1, nxi
        E(x,y,z) = N(x,y,z)%E(2)
      end do
    end do
  end do
  call dump_scalar(E, 'Ey')


  do z = 1, nzi
    do y = 1, nyi
      do x = 1, nxi
        E(x,y,z) = N(x,y,z)%E(3)
      end do
    end do
  end do
  call dump_scalar(E, 'Ez')

  deallocate(E, stat = ierror)
  call check_allocate(ierror, "Unable to deallocate scalar buffer for E output.")

end subroutine dump_E

! endif USEHDF
#endif

#endif
end module lbe_elec_output_module

