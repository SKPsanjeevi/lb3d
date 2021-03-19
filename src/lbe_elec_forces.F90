#include "lbe.h"

!> Force couplings for ion species -> fluids.
module lbe_elec_forces_module
#ifdef ELEC

  use lbe_elec_helper_module
  use lbe_elec_globals_module
  use lbe_force_interface_module, only: add_force_to_all
  use lbe_globals_module, only: lbe_force

  implicit none
  private

  public elec_add_forces, elec_force_reset
  public lbe_force

  contains

!> Calculate forces caused by the electrolytes and apply them to the fluids. Cf. manual for details.
subroutine elec_add_forces(N)
  implicit none

  type(lbe_site), intent(inout) :: N(1-halo_extent:, 1-halo_extent:, 1-halo_extent:) !> The lattice

  integer :: i, j, k, ip, jp, kp, veldir  !> Dummy loop variables
  real(kind=rk) :: Fx, Fy, Fz             !> Forces

  real(kind=rk) :: grad_rho(3)
  real(kind=rk) :: grad_E2(3)
  real(kind=rk) :: grad_solv(3)
  real(kind=rk) :: nb_E2, nb_solv, diff_rho

  real(kind=rk) :: F_static(3), F_diel(3), F_solv(3), F_ideal(3)

  if ( elec_on_fluid ) then

    ! Charge advection requires velocity calculation in the halo, which requires the forces to be set correctly in the halo,
    ! so we loop over the physical domain + 1 here.
    do i = 0, nx+1

      if ( inv_fluid .eq. 17 ) then
        ! On-site fixed velocity boundary condition - no forces at BCs, and upper part of the system.
        if ( ( i + ccoords(1)*nx .le. 1 ) .or. ( i + ccoords(1)*nx .ge. tnx - noforce_offset_ux ) ) then
          cycle
        end if
      end if
      
      do j = 0, ny+1
        do k = 0, nz+1

#ifndef ELEC_NEWFORCE
          ! Electrostatic term
          F_static(:) = ec * ( N(i,j,k)%rho_p - N(i,j,k)%rho_m ) * N(i,j,k)%E(:)

          ! Calculate contribution from local dieletric properties.
#ifdef ELEC_NNONREST
          call error_elec("ELEC_NNONREST not implemented for elec_add_forces")
#else
          grad_rho(:) = 0.0_rk
          grad_E2(:) = 0.0_rk
          grad_solv(:) = 0.0_rk

          do veldir = 1, nnn
            ip = i + cx(veldir)
            jp = j + cy(veldir)
            kp = k + cz(veldir)
            ! Ideal pressure term [US 24/02/2015]
            ! This term needs to be added since the charged species
            ! are not included in the LB algorithm which consequently
            ! does not capture their ideal pressure
            diff_rho = N(ip,jp,kp)%rho_p - N(i,j,k)%rho_p + N(ip,jp,kp)%rho_m - N(i,j,k)%rho_m
            grad_rho(1) = grad_rho(1) + cx(veldir) * diff_rho
            grad_rho(2) = grad_rho(2) + cy(veldir) * diff_rho
            grad_rho(3) = grad_rho(3) + cz(veldir) * diff_rho
#ifndef SINGLEFLUID
            ! Dielectrophoretic term
            nb_E2 = 0.5_rk * ( N(ip,jp,kp)%E(1)**2 + N(ip,jp,kp)%E(2)**2 + N(ip,jp,kp)%E(3)**2 )
            grad_E2(1) = grad_E2(1) + ( cx(veldir) * nb_E2 )
            grad_E2(2) = grad_E2(2) + ( cy(veldir) * nb_E2 )
            grad_E2(3) = grad_E2(3) + ( cz(veldir) * nb_E2 )
            ! Solvation term
            nb_solv = 0.5_rk * get_op(N,ip,jp,kp) * ( N(ip,jp,kp)%rho_p * delta_mu_plus + N(ip,jp,kp)%rho_m * delta_mu_minus )
            grad_solv(1) = grad_solv(1) + ( cx(veldir) * nb_solv )
            grad_solv(2) = grad_solv(2) + ( cy(veldir) * nb_solv )
            grad_solv(3) = grad_solv(3) + ( cz(veldir) * nb_solv )
#endif
          end do
#endif
          ! Ideal pressure term
          F_ideal = 0.5_rk * grad_rho

#ifndef SINGLEFLUID
          ! Dielectrophoretic term
          F_diel = 0.5_rk * ( get_local_eps(N,i,j,k) - eps_avg ) * grad_E2(:)

          ! Solvation term
          F_solv = 0.5_rk * grad_solv(:)
#endif

          Fx = F_ideal(1) + F_static(1) + F_diel(1) + F_solv(1)
          Fy = F_ideal(2) + F_static(2) + F_diel(2) + F_solv(2)
          Fz = F_ideal(3) + F_static(3) + F_diel(3) + F_solv(3)
#else
          Fx = N(i,j,k)%elec_force(1)
          Fy = N(i,j,k)%elec_force(2)
          Fz = N(i,j,k)%elec_force(3)
#endif
          call add_force_to_all(Fx, Fy, Fz, i, j, k)

        end do
      end do
    end do

  end if

end subroutine elec_add_forces

!> Like lbe_force_reset in lbe_force.F90, except also the halo is blanked out.
subroutine elec_force_reset
  implicit none

  lbe_force(:, 0:n_spec, 1-force_halo_extent:nx+force_halo_extent, &
                       & 1-force_halo_extent:ny+force_halo_extent, &
                       & 1-force_halo_extent:nz+force_halo_extent) = 0.0_rk

end subroutine elec_force_reset

#endif
end module lbe_elec_forces_module
