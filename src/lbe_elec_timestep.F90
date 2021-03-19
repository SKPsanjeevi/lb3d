#include "lbe.h"

!> Contains the main ELEC timestep.
module lbe_elec_timestep_module
#ifdef ELEC

  use lbe_elec_fluxes_module
  use lbe_elec_globals_module
  use lbe_elec_poisson_solver_module
  use lbe_elec_parallel_module

  implicit none

  private

  public elec_timestep

contains

!> One call calculates fluxes, moves charges accordingly and then calculates the new potential.
subroutine elec_timestep(N, ch_eq, ch_loops)
  implicit none

  type(lbe_site), intent(inout)  :: N(1-halo_extent:, 1-halo_extent:, 1-halo_extent:)
  logical, intent(out), optional :: ch_eq
  integer, intent(in),  optional :: ch_loops

  ! With advection we need to allocate a halo of one to avoid complicated code.
  real(kind=rk), dimension(1-flux_halo_extent:nx+flux_halo_extent, &
                         & 1-flux_halo_extent:ny+flux_halo_extent, &
                         & 1-flux_halo_extent:nz+flux_halo_extent) :: flux_site_plus
  real(kind=rk), dimension(1-flux_halo_extent:nx+flux_halo_extent, &
                         & 1-flux_halo_extent:ny+flux_halo_extent, &
                         & 1-flux_halo_extent:nz+flux_halo_extent) :: flux_site_minus

  ! Initialize flux arrays for summing
  flux_site_plus( 1-flux_halo_extent:, 1-flux_halo_extent:, 1-flux_halo_extent:) = 0.0_rk
  flux_site_minus(1-flux_halo_extent:, 1-flux_halo_extent:, 1-flux_halo_extent:) = 0.0_rk

  call halo_exchange(N, elec_halo)

  ! Calculate the charge fluxes over the lattice links due to diffusion.
  ! call log_msg_elec("Calling calc_fluxes_diffusion(...)")
  call calc_fluxes_diffusion(N, flux_site_plus, flux_site_minus)

  if ( fluid_on_elec ) then
    ! Do not advect during initial equilibration.
    if ( nt > 0 ) then    
      ! Calculate the charge fluxes over the lattice links due to advection.
      ! call log_msg_elec("Calling calc_fluxes_advection(...)")
      call calc_fluxes_advection(N, flux_site_plus, flux_site_minus)
    end if
  end if

  ! Apply the calculated fluxes by moving the charges.
  ! call log_msg_elec("Calling apply_fluxes(...)")
  if ( present( ch_eq ) .and. present(ch_loops) ) then
    ! Call this during initial equilibration; it keeps iterating until ch_eq = .true.
    call apply_fluxes(N, flux_site_plus, flux_site_minus, ch_eq, ch_loops)
  else
    ! Call this at later stages.
    call apply_fluxes(N, flux_site_plus, flux_site_minus)
  end if

  call halo_exchange(N, elec_halo)

  ! Solve the Poisson equation
  ! call log_msg_elec("Calling solve_poisson(...)")
  call solve_poisson(N)

  if ( E_solver_id == E_solver_fd ) then
    ! Calculate the electric field from the current potential.
    ! call log_msg_elec("Calling calc_E_fd(...)")
    call calc_E_fd(N)
  end if

end subroutine elec_timestep

#endif
end module lbe_elec_timestep_module

