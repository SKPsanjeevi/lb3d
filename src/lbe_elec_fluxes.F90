#include "lbe.h"

!> Main functions pertaining to evolution of charge densities on the lattice.
module lbe_elec_fluxes_module
#ifdef ELEC

  use lbe_elec_forces_module, only: lbe_force
  use lbe_elec_globals_module
  use lbe_elec_helper_module
  use lbe_elec_timer_module
  use lbe_io_module, only: corrected_velocity
  use mpi

  implicit none

  private

  public calc_fluxes_diffusion, calc_fluxes_advection, apply_fluxes, flux_halo_extent

  ! The halo of the temporary flux_site_ has to be set. Charge advection requires us to have correct
  ! fluxes in the halo, so for that we need flux_halo_extent = 1.
  integer, parameter :: flux_halo_extent = 1

contains

!> Calculate the charge fluxes caused by diffusion and store them in the flux_site_ arrays.
subroutine calc_fluxes_diffusion(N, flux_site_plus, flux_site_minus)
  implicit none

  type(lbe_site), intent(inout) :: N(1-halo_extent:, 1-halo_extent:, 1-halo_extent:)
  real(kind=rk), dimension(1-flux_halo_extent:, 1-flux_halo_extent:, 1-flux_halo_extent:), intent(inout) :: flux_site_plus
  real(kind=rk), dimension(1-flux_halo_extent:, 1-flux_halo_extent:, 1-flux_halo_extent:), intent(inout) :: flux_site_minus

  real(kind=rk) :: phi_tot_cur, phi_tot_nn
  real(kind=rk) :: exp_dphi, inv_exp_dphi
  real(kind=rk) :: flux_link_plus, flux_link_minus

  integer :: i, j, k, ip, jp, kp, tip, tjp, tkp, veldir

#ifndef SINGLEFLUID
  ! For fluid coupling (1-way)
  real(kind=rk) :: exp_delta_plus, inv_exp_delta_plus
  real(kind=rk) :: exp_delta_minus, inv_exp_delta_minus
  real(kind=rk) :: op_cur
#endif

  call start_timer(ti_elec_calc_fluxes_diffusion)

#ifdef ELEC_NEWFORCE
  ! elec forces are also needed in the halo?
  do i = 1, nx
    do j = 1, ny
      do k = 1, nz
#else
  ! Diffusive fluxes are required and thus calculated in the physical region only.
  do i = 1, nx
    do j = 1, ny
      do k = 1, nz
#endif
        ! reset elec_force on this node
        N(i,j,k)%elec_force = 0.0_rk
        if ( is_mobile_charge(N, i, j, k) ) then
          phi_tot_cur = N(i,j,k)%phi + ( Ex*real(i, kind=rk) + Ey*real(j, kind=rk) + Ez*real(k, kind=rk)  )
#ifndef SINGLEFLUID
          op_cur = get_op(N,i,j,k)
#endif
#ifdef ELEC_NNONREST
          do veldir = 1, nnonrest
            ip = i + cx(veldir)
            jp = j + cy(veldir)
            kp = k + cz(veldir)
            if ( is_mobile_charge(N, ip, jp, kp) ) then
              ! We probably do not want to overwrite the neighbour location, so use temps tip, tjp, tkp to store.
              phi_tot_nn = get_nb_phi(N, veldir, i, j,k, tip,tjp,tkp) + ( Ex*real(ip, kind=rk) + Ey*real(jp, kind=rk) + Ez*real(kp, kind=rk) )
#ifndef SINGLEFLUID
              if ( fluid_on_elec ) then
                call error_elec("Combination -DELEC_NNONREST + fluid_on_elec not implemented.")
              else
#endif
                exp_dphi = exp( beta * ec * ( phi_tot_nn - phi_tot_cur ) )
                exp_min_dphi = 1.0_rk / exp_dphi

                flux_link_plus  = -0.5_rk * D_plus * (1.0_rk + exp_min_dphi) * ( N(ip,jp,kp)%rho_p * exp_dphi - N(i,j,k)%rho_p )
                flux_link_minus = -0.5_rk * D_minus * (1.0_rk + exp_dphi) * ( N(ip,jp,kp)%rho_m * exp_min_dphi - N(i,j,k)%rho_m )

                ! WARNING: The 1 + 2 sqrt(2) prefactor is not sure to be correct, although it reproduces the correct results in the 1D liquid junction.
                flux_site_plus(i,j,k)  = flux_site_plus(i,j,k)  + flux_link_plus  / ( (1.0_rk + 2.0_rk * sqrt(2.0_rk) ) * flux_link_weight(veldir) )
                flux_site_minus(i,j,k) = flux_site_minus(i,j,k) + flux_link_minus / ( (1.0_rk + 2.0_rk * sqrt(2.0_rk) ) * flux_link_weight(veldir) )
#ifndef SINGLEFLUID
              end if ! fluid coupling
#endif
#ifdef ELEC_NEWFORCE
              ! site force to be applied in LB update
              N(i,j,k)%elec_force = N(i,j,k)%elec_force - 3.0_rk * w(veldir) * c(veldir,:) * (flux_link_plus/D_plus + flux_link_minus/D_minus)
#endif
#else
          do veldir = 1, nnn
            ip = i + cx(veldir)
            jp = j + cy(veldir)
            kp = k + cz(veldir)
            if ( is_mobile_charge(N, ip, jp, kp) ) then
              ! We probably do not want to overwrite the neighbour location, so use temps tip, tjp, tkp to store.
              phi_tot_nn = get_nb_phi(N, veldir, i, j, k,tip,tjp,tkp) + ( Ex*real(ip, kind=rk) + Ey*real(jp, kind=rk) + Ez*real(kp, kind=rk) )
#ifndef SINGLEFLUID
              if ( fluid_on_elec ) then
                exp_delta_plus  = exp( beta * ( ec * ( phi_tot_nn - phi_tot_cur ) + 0.5_rk * delta_mu_plus  * ( get_op(N,ip,jp,kp) - op_cur ) ) )
                inv_exp_delta_plus  = 1.0_rk / exp_delta_plus

                exp_delta_minus = exp( beta * ( ec * ( phi_tot_nn - phi_tot_cur ) + 0.5_rk * delta_mu_minus * ( get_op(N,ip,jp,kp) - op_cur ) ) )
                inv_exp_delta_minus = 1.0_rk / exp_delta_minus

                flux_link_plus  = -0.5_rk * D_plus  * (1.0_rk + inv_exp_delta_plus) * ( N(ip,jp,kp)%rho_p * exp_delta_plus - N(i,j,k)%rho_p )
                flux_link_minus = -0.5_rk * D_minus * (1.0_rk + exp_delta_minus) * ( N(ip,jp,kp)%rho_m * inv_exp_delta_minus - N(i,j,k)%rho_m )

                flux_site_plus(i,j,k)  = flux_site_plus(i,j,k)  + flux_link_plus
                flux_site_minus(i,j,k) = flux_site_minus(i,j,k) + flux_link_minus
              else
#endif
                exp_dphi = exp( beta * ec * ( phi_tot_nn - phi_tot_cur ) )
                inv_exp_dphi = 1.0_rk / exp_dphi

                flux_link_plus  = -0.5_rk * D_plus  * (1.0_rk + inv_exp_dphi) * ( N(ip,jp,kp)%rho_p * exp_dphi - N(i,j,k)%rho_p )
                flux_link_minus = -0.5_rk * D_minus * (1.0_rk + exp_dphi) * ( N(ip,jp,kp)%rho_m * inv_exp_dphi - N(i,j,k)%rho_m )

                flux_site_plus(i,j,k)  = flux_site_plus(i,j,k)  + flux_link_plus
                flux_site_minus(i,j,k) = flux_site_minus(i,j,k) + flux_link_minus
#ifndef SINGLEFLUID
              end if ! fluid coupling
#endif
#ifdef ELEC_NEWFORCE
              ! site force to be applied in LB update
              ! The factor 1/2 stems from w_i/c_s^2 for D3Q6 and D3Q7
              N(i,j,k)%elec_force = N(i,j,k)%elec_force - 0.5_rk * c(veldir,:) * (flux_link_plus/D_plus + flux_link_minus/D_minus)
#endif
#endif
            end if ! end norock at nn
          end do ! end nnonrest loop
        end if ! end norock at current
      end do ! end k loop
    end do ! ind j loop
  end do ! end i loop

  call stop_timer(ti_elec_calc_fluxes_diffusion)

end subroutine calc_fluxes_diffusion

!> Calculate the charge fluxes caused by advection and store them in the flux_site_ arrays.
subroutine calc_fluxes_advection(N, flux_site_plus, flux_site_minus)
  implicit none

  type(lbe_site), intent(inout) :: N(1-halo_extent:, 1-halo_extent:, 1-halo_extent:)
  real(kind=rk), dimension(1-flux_halo_extent:, 1-flux_halo_extent:, 1-flux_halo_extent:), intent(inout) :: flux_site_plus
  real(kind=rk), dimension(1-flux_halo_extent:, 1-flux_halo_extent:, 1-flux_halo_extent:), intent(inout) :: flux_site_minus

  integer :: i, j, k, ip, jp, kp, veldir
  real(kind=rk) :: vel(3), corr_vel(3,0:n_spec), delta_rho_plus, delta_rho_minus

  call start_timer(ti_elec_calc_fluxes_advection)

  ! To ensure charge conservation we might have to add fluxes moving from the halo into the physical domain, 
  ! so we loop over the physical domain + 1.
  do i = 0, nx+1
    do j = 0, ny+1
      do k = 0, nz+1
        if ( is_mobile_charge(N, i, j, k) ) then
          call corrected_velocity(N(i,j,k), lbe_force(:,0:n_spec,i,j,k), corr_vel(:,0:n_spec))
          vel = corr_vel(:,0)

#ifdef ELEC_NNONREST
          call error_elec("ELEC_NNONREST not implemented for calc_fluxes_advection.")
#else
          do veldir = 1, nnn
            ip = i + cx(veldir)
            jp = j + cy(veldir)
            kp = k + cz(veldir)
            if ( is_mobile_charge(N, ip, jp, kp) ) then
              delta_rho_plus  = 0.0_rk
              delta_rho_minus = 0.0_rk

              ! For one particular neighbour, we only have one of these if statements being true,
              ! so we can use else if blocks.
              ! The vel * c[xyz] product in the if-statements is taken to determine if the velocity
              ! and the current flux link are pointing in the same direction
              if ( vel(1) * cx(veldir) .gt. 0.0_rk ) then
                delta_rho_plus  = vel(1) * cx(veldir) * N(i,j,k)%rho_p
                delta_rho_minus = vel(1) * cx(veldir) * N(i,j,k)%rho_m
              else if ( vel(2) * cy(veldir) .gt. 0.0_rk ) then
                delta_rho_plus  = vel(2) * cy(veldir) * N(i,j,k)%rho_p
                delta_rho_minus = vel(2) * cy(veldir) * N(i,j,k)%rho_m
              else if ( vel(3) * cz(veldir) .gt. 0.0_rk ) then
                delta_rho_plus  = vel(3) * cz(veldir) * N(i,j,k)%rho_p
                delta_rho_minus = vel(3) * cz(veldir) * N(i,j,k)%rho_m
              end if

              flux_site_plus( i, j, k ) = flux_site_plus( i ,j ,k ) + delta_rho_plus
              flux_site_minus(i, j, k ) = flux_site_minus(i ,j ,k ) + delta_rho_minus
              ! We don't care what would happen in a second halo layer, so just don't update those sites.
              ! This avoids having to allocate an even larger halo just to solve an addressing problem.
              ! When summing over physical sites, the charge is still exactly conserved.
              if ( ip .ge. 1-flux_halo_extent .and. &
                 & jp .ge. 1-flux_halo_extent .and. &
                 & kp .ge. 1-flux_halo_extent .and. &
                 & ip .le. nx+flux_halo_extent .and. &
                 & jp .le. ny+flux_halo_extent .and. &
                 & kp .le. nz+flux_halo_extent ) then
                flux_site_plus( ip,jp,kp) = flux_site_plus( ip,jp,kp) - delta_rho_plus
                flux_site_minus(ip,jp,kp) = flux_site_minus(ip,jp,kp) - delta_rho_minus
              end if

            end if ! end norock at nn
          end do ! end nnonrest loop
#endif
        end if ! end norock at current
      end do ! end k loop
    end do ! ind j loop
  end do ! end i loop

  call stop_timer(ti_elec_calc_fluxes_advection)

end subroutine calc_fluxes_advection

!> Move charges on the lattice according to the flux_site_ arrays.
subroutine apply_fluxes(N, flux_site_plus, flux_site_minus, ch_eq_out, ch_loops)
  implicit none

  type(lbe_site), intent(inout) :: N(1-halo_extent:, 1-halo_extent:, 1-halo_extent:)
  real(kind=rk), dimension(1-flux_halo_extent:,1-flux_halo_extent:,1-flux_halo_extent:), intent(in) :: flux_site_plus
  real(kind=rk), dimension(1-flux_halo_extent:,1-flux_halo_extent:,1-flux_halo_extent:), intent(in) :: flux_site_minus
  logical, intent(out), optional :: ch_eq_out
  integer, intent(in), optional :: ch_loops

  integer :: mpierror

  real(kind=rk) :: tot_plus, tot_minus
  real(kind=rk) :: maxfluxl, maxflux_pl, maxflux_ml
  real(kind=rk) :: max_accuracy_cpl, max_accuracy_cml, max_accuracy_cp, max_accuracy_cm
  real(kind=rk) :: mf_try, accuracy_cpl, accuracy_cml, tot_plusl, tot_minusl
  integer :: i, j, k, ip, jp, kp, im, jm, km
  logical :: ch_eq

  call start_timer(ti_elec_apply_fluxes)

  ch_eq = .true.
  tot_plusl  = 0.0_rk
  tot_minusl = 0.0_rk
  maxfluxl   = 0.0_rk
  maxflux_pl = 0.0_rk
  maxflux_ml = 0.0_rk
  max_accuracy_cpl = 0.0_rk
  max_accuracy_cml = 0.0_rk

  ip = 0
  jp = 0
  kp = 0
  im = 0
  jm = 0
  km = 0

  ! Fluxes need only be applied in the physical region, so we loop over the physical region only.
  do i = 1, nx
    do j = 1, ny
      do k = 1, nz
        if ( is_mobile_charge(N, i, j, k) ) then

          ! Move charges
          N(i,j,k)%rho_p = N(i,j,k)%rho_p - flux_site_plus(i,j,k)
          N(i,j,k)%rho_m = N(i,j,k)%rho_m - flux_site_minus(i,j,k)

          ! Calculate accuracies
          if (abs(flux_site_plus(i,j,k)) > abs(flux_site_minus(i,j,k) ) ) then
            mf_try = abs(flux_site_plus(i,j,k))
          else
            mf_try = abs(flux_site_minus(i,j,k))
          end if

          if ( mf_try > maxfluxl) then
            maxfluxl = mf_try
            maxflux_pl = abs(flux_site_plus(i,j,k))
            maxflux_ml = abs(flux_site_minus(i,j,k))
          end if

          accuracy_cpl = abs(flux_site_plus(i,j,k) / N(i,j,k)%rho_p )
          if (max_accuracy_cpl < accuracy_cpl) then
            max_accuracy_cpl = accuracy_cpl
            ip = i
            jp = j
            kp = k
          end if

          accuracy_cml = abs(flux_site_minus(i,j,k) / N(i,j,k)%rho_m )
          if (max_accuracy_cml < accuracy_cml) then
            max_accuracy_cml = accuracy_cml
            im = i
            jm = j
            km = k
          end if

          tot_plusl  = tot_plusl  + flux_site_plus(i,j,k)
          tot_minusl = tot_minusl + flux_site_minus(i,j,k)
        end if ! fluid
      end do ! k 
    end do ! j
  end do ! i

  ! Check accuracies.
  call MPI_Allreduce(max_accuracy_cpl, max_accuracy_cp, 1, LBE_REAL, MPI_MAX, Comm_cart, mpierror)
  call MPI_Allreduce(max_accuracy_cml, max_accuracy_cm, 1, LBE_REAL, MPI_MAX, Comm_cart, mpierror)

  if ( max_accuracy_cp > acc_fluxes .or. max_accuracy_cm > acc_fluxes ) then
    ch_eq = .false.
  end if

  if ( present(ch_eq_out) .and. present(ch_loops) ) then
    ch_eq_out = ch_eq

    if ( n_show_eq > 0 ) then
      if( ( mod(ch_loops, n_show_eq) .eq. 0 ) ) then
        if ( .not. ch_eq ) then
          write(msgstr,"('  Fluxes not yet small enough: cp = ',ES15.8,' , cm = ',ES15.8, ' acc_fluxes = ', ES15.8)") max_accuracy_cp, max_accuracy_cm, acc_fluxes
          call log_msg_elec(msgstr)
        end if
      end if
    end if
  end if

  if ( dbg_check_zero_flux ) then
    ! Check for zero total flux.
    call MPI_Allreduce(tot_plusl , tot_plus , 1, LBE_REAL, MPI_SUM, Comm_cart, mpierror)
    call MPI_Allreduce(tot_minusl, tot_minus, 1, LBE_REAL, MPI_SUM, Comm_cart, mpierror)

    if ( abs(tot_plus) > meps .or. abs(tot_minus) > meps ) then
      write(msgstr,"('Nonzero total flux: ',2(ES15.8,X), ' Local fluxes: ',2(ES15.8,X))") tot_plus, tot_minus, tot_plusl, tot_minusl
      call error_elec(msgstr)
    else
      ! write(msgstr,"('  Total fluxes: ',2(ES15.8,X))") tot_plus, tot_minus
      ! call log_msg_elec(msgstr)
      ! write(msgstr,"('  Total flux: ',2(ES15.8,X), ' Local fluxes: ',2(ES15.8,X))") tot_plus, tot_minus, tot_plusl, tot_minusl
      ! call log_msg_elec(msgstr,.true.)
    end if
  end if

  call stop_timer(ti_elec_apply_fluxes)

end subroutine apply_fluxes

logical function is_mobile_charge(N, i, j, k)
  implicit none

  type(lbe_site), intent(in) :: N(1-halo_extent:, 1-halo_extent:, 1-halo_extent:)
  integer, intent(in) :: i, j, k

  if ( inv_fluid .eq. 17 ) then
    ! In combination with on-site fixed velocity boundary condition - no moving charges at top and bottom BCs.
    if ( ( i + ccoords(1)*nx .eq. 1 ) .or. ( i + ccoords(1)*nx .eq. tnx ) ) then
      is_mobile_charge = .false.
    else
      is_mobile_charge = is_fluid( N(i,j,k)%rock_state )
    end if
  else
    is_mobile_charge = is_fluid( N(i,j,k)%rock_state )
  end if
end function is_mobile_charge

#endif
end module lbe_elec_fluxes_module

