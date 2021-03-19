#include "lbe.h"

!> Contains functions used for initialising the ELEC part of the system, and destruction and cleanup afterwards.
module lbe_elec_init_module
#ifdef ELEC

  use lbe_elec_fluxes_module
  use lbe_elec_helper_module
  use lbe_elec_globals_module
  use lbe_elec_output_module, only: elec_dump_all_asc
  use lbe_elec_parallel_module
  use lbe_elec_poisson_solver_module
  use lbe_elec_timer_module
  use lbe_elec_timestep_module
  use lbe_globals_module, only: rock_value
  use lbe_helper_module, only: request_halo_extent, n_sites, n_sites_fluid, n_sites_rock, n_sites_surface, n_sites_particle
#ifdef MD
  use lbe_md_helper_module, only: count_particles_all
#endif
  use mpi

  implicit none

  private

  public elec_init_system, elec_restore_init_system, elec_initial_equilibration, elec_shutdown

contains

!> This is the general initialization routine. All ELEC-specific fields should be set up here.
subroutine elec_init_system(N)
  implicit none

  type(lbe_site), intent(inout) :: N(1-halo_extent:, 1-halo_extent:, 1-halo_extent:)

  call start_timer(ti_elec_init)

  call log_msg_elec_hdr("Initializing ELEC fields")

#ifdef P3M
  if ( poisson_solver_id == poisson_solver_p3m .and. local_eps ) then
    call log_msg_elec("Allocating and initializing E_prev array.")
    call lbe_elec_allocate_E_prev(N)
    call log_msg_elec("Successfully allocated and initialized E_prev.")
  end if
#endif

  ! Initialize zero potential field.
  call log_msg_elec("Initializing phi and E ...")
  call lbe_elec_init_phi_zero(N)
  call log_msg_elec("Successfully initialized phi and E!")

  ! Initialize charge fields.
  call log_msg_elec("Initializing rho ...")

  select case( trim(rho_init) )
    case( "capacitor" )
      call lbe_elec_init_rho_capacitor(N)
    case( "liquidjunction" )
      call lbe_elec_init_rho_liquid_junction(N)
    case( "slipbc" )
      call lbe_elec_init_rho_slipbc(N)
    case( "uniform" )
      call lbe_elec_init_rho_uniform(N)
    case( "surface" )
       call lbe_elec_init_rho_surface(N)
    case default
      write(msgstr,"('  Invalid rho_init value <',A,'> .')") trim(rho_init)
      call error_elec(msgstr)
  end select

  call log_msg_elec("Successfully initialized rho!")

  ! If we are using local permittivity, initialize the eps field.
  if ( local_eps .and. ( .not. fluid_on_elec ) ) then

    call log_msg_elec("Initializing local epsilon ...")

    select case( trim(eps_init) )
      case( "uniform" )
        call lbe_elec_init_eps_uniform(N)
      case( "capacitor" )
        call lbe_elec_init_eps_capacitor(N)
      case( "colloid" )
        call lbe_elec_init_eps_colloid(N)
      case default
        write(msgstr,"('  Invalid eps_init value <',A,'> .')") trim(eps_init)
        call error_elec(msgstr)
    end select

    call log_msg_elec("Successfully initialized local epsilon!")

  end if

  call log_msg_elec("Initializing rock state ...")

  ! If we want to modify existing rock state (as set in the normal LB3D rock
  ! initialization), overwrite the rock field.
  select case( trim(rock_init) )
    case ( "none" )
      call log_msg_elec("  Initializing rock with option <none>.")
    case( "capacitor" )
      call lbe_elec_init_rock_capacitor(N)
    case default
      write(msgstr,"('  Invalid rock_init value <',A,'> .')") trim(rock_init)
      call error_elec(msgstr)
  end select

  call log_msg_elec("Successfully initialized ELEC rock modifications!")

  ! Dump to see if all went well.
  call elec_dump_all_asc(N, "post-init")

  call stop_timer(ti_elec_init)

  call log_msg_elec_ws("Finished initializing ELEC fields.")

end subroutine elec_init_system

subroutine elec_restore_init_system(N)
  implicit none

  type(lbe_site), intent(inout) :: N(1-halo_extent:, 1-halo_extent:, 1-halo_extent:)

  call halo_exchange(N, elec_halo)
  call solve_poisson(N)
  call calc_E_fd(N)

end subroutine elec_restore_init_system

#ifdef P3M
! ----------------------------------------------------------------------------
!
!                              INIT E_PREV
!
! ----------------------------------------------------------------------------

subroutine lbe_elec_allocate_E_prev(N)
  implicit none

  type(lbe_site), intent(in) :: N(1-halo_extent:, 1-halo_extent:, 1-halo_extent:)

  integer :: ierror

  allocate( E_prev(-1:size(N,1)+2, -1:size(N,2)+2, -1:size(N,3)+2, 3 ), stat = ierror)
  call check_allocate(ierror, "Unable to allocate E_prev.")

  E_prev(-1:,-1:,-1:,:) = 0.0_rk

end subroutine lbe_elec_allocate_E_prev

#endif

! ----------------------------------------------------------------------------
!
!                                INIT PHI
!
! ----------------------------------------------------------------------------

!> This just zeroes out the potential field (will be calculated from the charges
!> on the first equilibration step).
subroutine lbe_elec_init_phi_zero(N)
  implicit none

  type(lbe_site), intent(inout) :: N(1-halo_extent:, 1-halo_extent:, 1-halo_extent:)

  integer :: i, j, k

  call log_msg_elec("  Initializing phi and E to zero.")

  do i = 1, nx
    do j = 1, ny
      do k = 1, nz
        N(i,j,k)%phi = 0.0_rk
        N(i,j,k)%E(:) = 0.0_rk
      end do
    end do
  end do
end subroutine lbe_elec_init_phi_zero

! ----------------------------------------------------------------------------
!
!                                INIT RHO
!
! ----------------------------------------------------------------------------

!> Initialize the charge field in such a way that charge is distributed evenly over
!> walls and colloids (according to Q_wall and Q_colloid), with charge to compensate
!> distributed evenly over fluid sites. This also takes into account salt concentration.
subroutine lbe_elec_init_rho_uniform(N)
  implicit none

  type(lbe_site), intent(inout) :: N(1-halo_extent:, 1-halo_extent:, 1-halo_extent:)

  integer :: x, y, z
  integer :: colloidsites,  wallsites,  fluidsites,  totalsites  !< global count of site types
  integer :: colloidsitesl, wallsitesl, fluidsitesl, totalsitesl !< local count of site types
  integer :: mpierror

  real(kind=rk) :: colloidcharge, wallcharge, totalcharge, verifycharge, verifychargel
  real(kind=rk) :: con_salt

  real(kind=rk) :: rho_colloid, rho_p_colloid, rho_m_colloid
  real(kind=rk) :: rho_wall,    rho_p_wall,    rho_m_wall
  real(kind=rk) :: rho_fluid,   rho_p_fluid,   rho_m_fluid

  call log_msg_elec("  Initializing rho with option <uniform>.")

#ifdef MD
  call count_particles_all(n_colloids)
  write(msgstr,"('  Using ',I0,' colloid(s) from MD.')") n_colloids
  call log_msg_elec(msgstr)
#endif

  ! Count site types
  colloidsitesl = 0
  wallsitesl    = 0
  fluidsitesl   = 0
  totalsitesl   = 0
  do x = 1, nx
    do y = 1, ny
      do z = 1, nz
        if ( is_colloid( N(x,y,z)%rock_state ) ) then
          colloidsitesl = colloidsitesl + 1
        else if ( is_wall( N(x,y,z)%rock_state ) ) then
          wallsitesl = wallsitesl + 1
        else
          fluidsitesl = fluidsitesl + 1
        end if
        totalsitesl = totalsitesl + 1
      end do
    end do
  end do

  call MPI_Allreduce(colloidsitesl, colloidsites, 1, MPI_INTEGER, MPI_SUM, Comm_cart, mpierror)
  call MPI_Allreduce(wallsitesl, wallsites, 1, MPI_INTEGER, MPI_SUM, Comm_cart, mpierror)
  call MPI_Allreduce(fluidsitesl, fluidsites, 1, MPI_INTEGER, MPI_SUM, Comm_cart, mpierror)
  call MPI_Allreduce(totalsitesl, totalsites, 1, MPI_INTEGER, MPI_SUM, Comm_cart, mpierror)

  write(msgstr,"('  Fluid sites: ',I0,', colloid sites: ',I0,', wall sites: ',I0,', total: ',I0)") fluidsites, colloidsites, wallsites, totalsites
  call log_msg_elec(msgstr)

  ! Calculate total charges
  colloidcharge = real(n_colloids, kind=rk)*Q_colloid

  if ( wallsites > 0 ) then
    wallcharge = Q_wall
  else
    wallcharge = 0.0_rk
  end if

  totalcharge = colloidcharge + wallcharge

  write(msgstr,"('  Number of colloids: ',I0,', charge per colloid: ',ES15.8,', total charge on colloids: ',ES15.8)") n_colloids, Q_colloid, colloidcharge
  call log_msg_elec(msgstr)
  write(msgstr,"('  Wall charge: ',ES15.8,', total charge: ',ES15.8)") wallcharge, totalcharge
  call log_msg_elec(msgstr)

  ! Calculate charge densities to be placed on fluid sites.
  rho_fluid = -totalcharge / real(fluidsites, kind=rk)

  write(msgstr,"('  Counterion concentration: ',ES15.8)") rho_fluid
  call log_msg_elec(msgstr)

  ! Add salt, if debye_length is not zero.
  if ( debye_length .ne. 0.0_rk ) then
     con_salt = 1.0_rk / ( 4.0_rk * pi * bjerrum_length * debye_length * debye_length )
  else
     con_salt = 0.0_rk
  end if

  write(msgstr,"('  Salt concentration: ',ES15.8)") con_salt
  call log_msg_elec(msgstr)

  ! Calculate charge densities to be placed on colloid sites.
  if ( colloidsites > 0 ) then
    rho_colloid = colloidcharge / real(colloidsites, kind=rk)
    rho_p_colloid =  0.5_rk * rho_colloid
    rho_m_colloid = -0.5_rk * rho_colloid
  else
    rho_p_colloid = 0.0_rk
    rho_m_colloid = 0.0_rk
  end if

  ! Calculate charge densities to be placed on wall sites.
  if ( wallsites > 0 ) then
    rho_wall = Q_wall / real(wallsites, kind=rk)
    rho_p_wall =  0.5_rk*rho_wall
    rho_m_wall = -0.5_rk*rho_wall
  else
    rho_p_wall = 0.0_rk
    rho_m_wall = 0.0_rk
  end if

  if (rho_fluid >= 0.0_rk) then
     rho_p_fluid = 0.5_rk * con_salt + rho_fluid
     rho_m_fluid = 0.5_rk * con_salt
  else
     rho_p_fluid = 0.5_rk * con_salt
     rho_m_fluid = 0.5_rk * con_salt - rho_fluid
  end if

  write(msgstr,"('  rho_p_fluid   = ', ES15.8,' , rho_m_fluid   = ', ES15.8)") rho_p_fluid, rho_m_fluid
  call log_msg_elec(msgstr)
  write(msgstr,"('  rho_p_colloid = ', ES15.8,' , rho_m_colloid = ', ES15.8)") rho_p_colloid, rho_m_colloid
  call log_msg_elec(msgstr)
  write(msgstr,"('  rho_p_wall    = ', ES15.8,' , rho_m_wall    = ', ES15.8)") rho_p_wall, rho_m_wall
  call log_msg_elec(msgstr)

  ! Now loop over the lattice and assign the charge densities.
  call log_msg_elec("  Initializing system ...")
  do x = 1, nx
    do y = 1, ny
      do z = 1, nz
        if ( is_colloid( N(x,y,z)%rock_state ) ) then
          N(x,y,z)%rho_p = rho_p_colloid
          N(x,y,z)%rho_m = rho_m_colloid
        else if ( is_wall( N(x,y,z)%rock_state ) ) then
          N(x,y,z)%rho_p = rho_p_wall
          N(x,y,z)%rho_m = rho_m_wall
        else
          N(x,y,z)%rho_p = rho_p_fluid
          N(x,y,z)%rho_m = rho_m_fluid
        end if
      end do
    end do
  end do

  ! Check for overall neutrality of the system.
  if ( acc_neutrality == 0.0_rk ) then
    call log_msg_elec("  WARNING: acc_neutrality identically zero. Please choose a value larger than zero to perform a neutrality check.")
  else
    call log_msg_elec("  Checking for neutrality ...")

    verifychargel = 0.0_rk
    do x = 1, nx
      do y = 1, ny
        do z = 1, nz
          verifychargel = verifychargel + N(x,y,z)%rho_p - N(x,y,z)%rho_m
        end do
      end do
    end do

    call MPI_Allreduce(verifychargel, verifycharge, 1, LBE_REAL, MPI_SUM, Comm_cart, mpierror)
    ! write(msgstr,"('  verifychargel = ', ES15.8,' , verifycharge = ',ES15.8)") verifychargel, verifycharge
    ! call log_msg_elec(msgstr)
    ! write(msgstr,"('  verifycharge = ',ES15.8, ' , acc_neutrality = ',ES15.8 )") verifycharge, acc_neutrality
    ! call log_msg_elec(msgstr)

    if ( myrankc == 0 ) then
      if ( abs(verifycharge) > acc_neutrality ) then
        call error_elec("Neutrality violated, aborting...")
      else
        call log_msg_elec("    System is neutral up to acc_neutrality.")
      end if
    end if
  end if

end subroutine lbe_elec_init_rho_uniform

!> Initialize the charges such that Q_wall is distributed over surface nodes
subroutine lbe_elec_init_rho_surface(N)
  implicit none

  type(lbe_site), intent(inout) :: N(1-halo_extent:, 1-halo_extent:, 1-halo_extent:)

  integer :: x, y, z
  integer :: colloidsites, wallsites, fluidsites, surfacesites, totalsites
  integer :: mpierror

  real(kind=rk) :: colloidcharge, wallcharge, totalcharge, verifycharge, verifychargel
  real(kind=rk) :: con_salt

  real(kind=rk) :: rho_colloid, rho_p_colloid, rho_m_colloid
  real(kind=rk) :: rho_wall,    rho_p_wall,    rho_m_wall
  real(kind=rk) :: rho_fluid,   rho_p_fluid,   rho_m_fluid

  fluidsites   = n_sites_fluid
  colloidsites = n_sites_particle
  wallsites    = n_sites_rock
  surfacesites = n_sites_surface
  totalsites   = n_sites

  write(msgstr,"('  Fluid sites: ',I0,', colloid sites: ',I0,', wall sites: ',I0,', surface sites: ',I0,', total: ',I0)") fluidsites, colloidsites, wallsites, surfacesites, totalsites
  call log_msg_elec(msgstr)

  ! Calculate total charges
  colloidcharge = real(n_colloids, kind=rk)*Q_colloid

  if ( wallsites > 0 ) then
    wallcharge = Q_wall
  else
    wallcharge = 0.0_rk
  end if

  totalcharge = colloidcharge + wallcharge

  write(msgstr,"('  Number of colloids: ',I0,', charge per colloid: ',ES15.8,', total charge on colloids: ',ES15.8)") n_colloids, Q_colloid, colloidcharge
  call log_msg_elec(msgstr)
  write(msgstr,"('  Wall charge: ',ES15.8,', total charge: ',ES15.8)") wallcharge, totalcharge
  call log_msg_elec(msgstr)

  ! Calculate charge densities to be placed on fluid sites.
  rho_fluid = -totalcharge / real(fluidsites, kind=rk)

  write(msgstr,"('  Counterion concentration: ',ES15.8)") rho_fluid
  call log_msg_elec(msgstr)

  ! Add salt, if debye_length is not zero.
  if ( debye_length .ne. 0.0_rk ) then
     con_salt = 1.0_rk / ( 4.0_rk * pi * bjerrum_length * debye_length * debye_length )
  else
     con_salt = 0.0_rk
  end if

  write(msgstr,"('  Salt concentration: ',ES15.8)") con_salt
  call log_msg_elec(msgstr)

  ! Calculate charge densities to be placed on colloid sites.
  if ( colloidsites > 0 ) then
    rho_colloid = colloidcharge / real(colloidsites, kind=rk)
    rho_p_colloid =  0.5_rk * rho_colloid
    rho_m_colloid = -0.5_rk * rho_colloid
  else
    rho_p_colloid = 0.0_rk
    rho_m_colloid = 0.0_rk
  end if

  ! Calculate charge densities to be placed on wall sites.
  if ( wallsites > 0 ) then
    rho_wall = Q_wall / real(surfacesites, kind=rk)
    rho_p_wall =  0.5_rk*rho_wall
    rho_m_wall = -0.5_rk*rho_wall
  else
    rho_p_wall = 0.0_rk
    rho_m_wall = 0.0_rk
  end if

  if (rho_fluid >= 0.0_rk) then
     rho_p_fluid = 0.5_rk * con_salt + rho_fluid
     rho_m_fluid = 0.5_rk * con_salt
  else
     rho_p_fluid = 0.5_rk * con_salt
     rho_m_fluid = 0.5_rk * con_salt - rho_fluid
  end if

  write(msgstr,"('  rho_p_fluid   = ', ES15.8,' , rho_m_fluid   = ', ES15.8)") rho_p_fluid, rho_m_fluid
  call log_msg_elec(msgstr)
  write(msgstr,"('  rho_p_colloid = ', ES15.8,' , rho_m_colloid = ', ES15.8)") rho_p_colloid, rho_m_colloid
  call log_msg_elec(msgstr)
  write(msgstr,"('  rho_p_wall    = ', ES15.8,' , rho_m_wall    = ', ES15.8)") rho_p_wall, rho_m_wall
  call log_msg_elec(msgstr)

  ! Now loop over the lattice and assign the charge densities.
  call log_msg_elec("  Initializing system ...")
  do x = 1, nx
    do y = 1, ny
      do z = 1, nz
        if ( is_colloid( N(x,y,z)%rock_state ) ) then
          N(x,y,z)%rho_p = rho_p_colloid
          N(x,y,z)%rho_m = rho_m_colloid
        else if ( is_rock( N(x,y,z)%rock_state ) ) then
           if ( is_surface( N(x,y,z)%rock_state ) ) then
              N(x,y,z)%rho_p = rho_p_wall
              N(x,y,z)%rho_m = rho_m_wall
           else
              N(x,y,z)%rho_p = 0.0_rk
              N(x,y,z)%rho_m = 0.0_rk
           end if
        else
          N(x,y,z)%rho_p = rho_p_fluid
          N(x,y,z)%rho_m = rho_m_fluid
        end if
      end do
    end do
  end do

  ! Check for overall neutrality of the system.
  if ( acc_neutrality == 0.0_rk ) then
    call log_msg_elec("  WARNING: acc_neutrality identically zero. Please choose a value larger than zero to perform a neutrality check.")
  else
    call log_msg_elec("  Checking for neutrality ...")

    verifychargel = 0.0_rk
    do x = 1, nx
      do y = 1, ny
        do z = 1, nz
          verifychargel = verifychargel + N(x,y,z)%rho_p - N(x,y,z)%rho_m
        end do
      end do
    end do

    call MPI_Allreduce(verifychargel, verifycharge, 1, LBE_REAL, MPI_SUM, Comm_cart, mpierror)
    ! write(msgstr,"('  verifychargel = ', ES15.8,' , verifycharge = ',ES15.8)") verifychargel, verifycharge
    ! call log_msg_elec(msgstr)
    ! write(msgstr,"('  verifycharge = ',ES15.8, ' , acc_neutrality = ',ES15.8 )") verifycharge, acc_neutrality
    ! call log_msg_elec(msgstr)

    if ( myrankc == 0 ) then
      if ( abs(verifycharge) > acc_neutrality ) then
        call error_elec("Neutrality violated, aborting...")
      else
        call log_msg_elec("    System is neutral up to acc_neutrality.")
      end if
    end if
  end if

end subroutine lbe_elec_init_rho_surface

!> Initialize the charge fields to simulate a liquid junction. All sites will be neutral, but
!> the densities of both ion species will be offset by delta_rho_lj on the top and bottom z-half of the system.
subroutine lbe_elec_init_rho_liquid_junction(N)
  implicit none

  type(lbe_site), intent(inout) :: N(1-halo_extent:, 1-halo_extent:, 1-halo_extent:)

  integer :: x, y, z

  call log_msg_elec("  Initializing rho with option <liquidjunction>.")

  ! Loop over the lattice and initialize the charge densities.
  do x = 1, nx
    do y = 1, ny
      do z = 1, nz
        if ( z + ccoords(3)*nz .le. tnz/2 ) then
          N(x,y,z)%rho_p = rho_lj + delta_rho_lj
          N(x,y,z)%rho_m = rho_lj + delta_rho_lj
        else
          N(x,y,z)%rho_p = rho_lj - delta_rho_lj
          N(x,y,z)%rho_m = rho_lj - delta_rho_lj
        end if
      end do
    end do
  end do

end subroutine lbe_elec_init_rho_liquid_junction

!> Initialize the charge fields to simulate a capacitor: initialize two charged x-planes
!> with charge +1 and -1, respectively.
!> WARNING: should probably be rewritten to match the normal LB3D rock init boundary_cond = 6
!> for placing charges.
subroutine lbe_elec_init_rho_capacitor(N)
  implicit none

  type(lbe_site), intent(inout) :: N(1-halo_extent:, 1-halo_extent:, 1-halo_extent:)

  integer :: x, y, z
  real(kind=rk) :: sigma

  sigma = Q_wall / ( real(tnx,kind=rk) * real(tny,kind=rk) )

  call log_msg_elec("  Initializing rho with option <capacitor>.")
  write(msgstr,"('   Q_wall = ',E16.8, ', sigma = ',E16.8)") Q_wall, sigma
  call log_msg_elec(msgstr)

  do x = 1, nx
    do y = 1, ny
      do z = 1, nz
        N(x,y,z)%rho_p = 0.0_rk
        N(x,y,z)%rho_m = 0.0_rk

        if ( z + ccoords(3)*nz .eq. 1*tnz/4 ) then
          N(x,y,z)%rho_p = -0.5_rk * sigma
          N(x,y,z)%rho_m = +0.5_rk * sigma
        end if

        if ( z + ccoords(3)*nz .eq. 3*tnz/4 + 1 ) then
          N(x,y,z)%rho_p = +0.5_rk * sigma
          N(x,y,z)%rho_m = -0.5_rk * sigma
        end if

      end do
    end do
  end do

end subroutine lbe_elec_init_rho_capacitor

!> Initialize the charge field in such a way that charge density Q_slip is set on slip areas and
!> Q_noslip is set on noslip areas (matching inv_fluid == 17), for x = 1. For x = tnx charge is set to zero.
!> Compentsating charges to attain neutrality are distributed evenly over fluid sites. 
!> This initialization also takes into account salt concentration.
!> No other obstacles are currently supported.
subroutine lbe_elec_init_rho_slipbc(N)
  implicit none

  type(lbe_site), intent(inout) :: N(1-halo_extent:, 1-halo_extent:, 1-halo_extent:)

  integer :: x, y, z, ty
  integer :: boundarysites, fluidsites, totalsites  !< global count of site types
  integer :: mpierror

  real(kind=rk) :: boundarycharge, boundarychargel, verifycharge, verifychargel
  real(kind=rk) :: con_salt

  real(kind=rk) :: rho_boundary,  rho_p_boundary, rho_m_boundary
  real(kind=rk) :: rho_fluid, rho_p_fluid, rho_m_fluid
  real(kind=rk) :: Q_stripe

  integer :: ipb, ipg

  call log_msg_elec("  Initializing rho with option <slipbc>.")

  ! Convert to integer
  ipb = pb
  ipg = pg

  ! Count site types - this assumes no other obstacles are present.
  boundarysites = 2*tny*tnz
  fluidsites = (tnx-2)*tny*tnz
  totalsites = boundarysites + fluidsites

  write(msgstr,"('  Fluid sites: ',I0,', boundary sites: ',I0,', total: ',I0)") fluidsites, boundarysites, totalsites
  call log_msg_elec(msgstr)

  ! Set charges on boundary layers x = 1 and x = tnx, according to stripe width parameters pb and pg,
  ! and charge densities Q_slip and Q_noslip.
  boundarychargel = 0.0_rk
  do x = 1, nx
    if ( x + ccoords(1)*nx .eq. 1 ) then
      Q_stripe = Q_slip
      do y = 1, ny
        ty = ccoords(2) * ny + y
        if ( mod(ty,ipb) == 0) then
          Q_stripe = Q_slip
        end if
        if ( mod(ty,ipb) >= ipg) then
          Q_stripe = Q_noslip
        end if

        do z = 1, nz
          boundarychargel = boundarychargel + Q_stripe
          N(x,y,z)%rho_p =  0.5_rk*Q_stripe
          N(x,y,z)%rho_m = -0.5_rk*Q_stripe
        end do
      end do
    else if ( x + ccoords(1)*nx .eq. tnx ) then
      do y = 1, ny
        do z = 1, nz
          N(x,y,z)%rho_p = 0.0_rk
          N(x,y,z)%rho_m = 0.0_rk
        end do
      end do
    end if
  end do

  ! Make all ranks aware of total charge placed on the boundary.
  call MPI_Allreduce(boundarychargel, boundarycharge, 1, LBE_REAL, MPI_SUM, Comm_cart, mpierror)

  write(msgstr,"('  Boundary charge: ',ES15.8)") boundarycharge
  call log_msg_elec(msgstr)

  ! Add salt, if debye_length is not zero.
  if ( debye_length .ne. 0.0_rk ) then
     con_salt = 1.0_rk / ( 4.0_rk * pi * bjerrum_length * debye_length * debye_length )
  else
     con_salt = 0.0_rk
  end if

  ! Calculate charge densities to be placed on fluid sites.
  rho_fluid = -boundarycharge / real(fluidsites, kind=rk)

  write(msgstr,"('  Counterion concentration: ',ES15.8)") rho_fluid
  call log_msg_elec(msgstr)

  write(msgstr,"('  Salt concentration: ',ES15.8)") con_salt
  call log_msg_elec(msgstr)

  if (rho_fluid >= 0.0_rk) then
     rho_p_fluid = 0.5_rk * con_salt + rho_fluid
     rho_m_fluid = 0.5_rk * con_salt
  else
     rho_p_fluid = 0.5_rk * con_salt
     rho_m_fluid = 0.5_rk * con_salt - rho_fluid
  end if

  write(msgstr,"('  rho_p_fluid   = ', ES15.8,' , rho_m_fluid   = ', ES15.8)") rho_p_fluid, rho_m_fluid
  call log_msg_elec(msgstr)

  ! Now loop over the lattice and assign the charge densities. Do not touch x = 1 and x = tnx again.
  call log_msg_elec("  Initializing system ...")
  do x = 1, nx
    if ( ( x + ccoords(1)*nx .gt. 1 ) .and. ( x + ccoords(1)*nx .lt. tnx ) ) then
      do y = 1, ny
        do z = 1, nz
          N(x,y,z)%rho_p = rho_p_fluid
          N(x,y,z)%rho_m = rho_m_fluid
        end do
      end do
    end if
  end do

  ! Check for overall neutrality of the system.
  if ( acc_neutrality == 0.0_rk ) then
    call log_msg_elec("  WARNING: acc_neutrality identically zero. Please choose a value larger than zero to perform a neutrality check.")
  else
    call log_msg_elec("  Checking for neutrality ...")

    verifychargel = 0.0_rk
    do x = 1, nx
      do y = 1, ny
        do z = 1, nz
          verifychargel = verifychargel + N(x,y,z)%rho_p - N(x,y,z)%rho_m
        end do
      end do
    end do

    call MPI_Allreduce(verifychargel, verifycharge, 1, LBE_REAL, MPI_SUM, Comm_cart, mpierror)
    ! write(msgstr,"('  verifychargel = ', ES15.8,' , verifycharge = ',ES15.8)") verifychargel, verifycharge
    ! call log_msg_elec(msgstr)
    ! write(msgstr,"('  verifycharge = ',ES15.8, ' , acc_neutrality = ',ES15.8 )") verifycharge, acc_neutrality
    ! call log_msg_elec(msgstr)

    if ( myrankc == 0 ) then
      if ( abs(verifycharge) > acc_neutrality ) then
        call error_elec("Neutrality violated, aborting...")
      else
        call log_msg_elec("    System is neutral up to acc_neutrality.")
      end if
    end if
  end if

end subroutine lbe_elec_init_rho_slipbc

! ----------------------------------------------------------------------------
!
!                                INIT EPS
!
! ----------------------------------------------------------------------------

!> Uniformly initialize the local eps.
!> Setting negative value for eps_uniform should be equivalent to using local_eps = .false. (legacy behaviour).
subroutine lbe_elec_init_eps_uniform(N)
  implicit none

  type(lbe_site), intent(inout) :: N(1-halo_extent:, 1-halo_extent:, 1-halo_extent:)

  integer :: i, j, k

  call log_msg_elec("  Initializing eps with option <uniform>, setting eps = eps_uniform everywhere.")

  ! Set eps on the lattice
  do i = 1, nx
    do j = 1, ny
      do k = 1, nz
        N(i,j,k)%eps = eps_uniform
      end do
    end do
  end do

end subroutine lbe_elec_init_eps_uniform

!> Uniformly initialize the local eps, but with different values for fluid and colloid sites.
!> WARNING: hardcoded variables exist here, this should probably be rewritten to be more flexible.
subroutine lbe_elec_init_eps_colloid(N)
  implicit none

  type(lbe_site), intent(inout) :: N(1-halo_extent:, 1-halo_extent:, 1-halo_extent:)

  integer :: i, j, k
  real(kind=rk) :: eps_colloid, eps_fluid

  eps_colloid = 1.0_rk / ( pi * bjerrum_length)
  eps_fluid   = 0.8_rk / ( pi * bjerrum_length)

  call log_msg_elec("  Initializing eps with option <colloid>.")
  write(msgstr,"('  eps in colloid = ', ES15.8,' , eps in fluid = ', ES15.8)") eps_colloid, eps_fluid
  call log_msg_elec(msgstr)

  do i = 1, nx
    do j = 1, ny
      do k = 1, nz
        if ( is_colloid(N(i,j,k)%rock_state) ) then
          N(i,j,k)%eps = eps_colloid
        else
          N(i,j,k)%eps = eps_fluid
        end if
      end do
    end do
  end do

end subroutine lbe_elec_init_eps_colloid

!> Initialize the local eps in lamellae to create a capacitor test system.
subroutine lbe_elec_init_eps_capacitor(N)
  implicit none

  type(lbe_site), intent(inout) :: N(1-halo_extent:, 1-halo_extent:, 1-halo_extent:)

  integer :: i, j, k
  real(kind=rk) :: eps_top, eps_bottom

  eps_top    = 2.0_rk * eps_uniform
  eps_bottom = 1.0_rk * eps_uniform

  call log_msg_elec("  Initializing eps with option <capacitor>.")
  write(msgstr,"('  eps on top = ', ES15.8,' , eps at bottom = ',ES15.8)") eps_top, eps_bottom
  call log_msg_elec(msgstr)

  do i = 1, nx
    do j = 1, ny
      do k = 1, nz
        if ( k + ccoords(3)*nz .gt. tnz / 2 ) then
          N(i,j,k)%eps = eps_top
        else
          N(i,j,k)%eps = eps_bottom
        end if
      end do
    end do
  end do

end subroutine lbe_elec_init_eps_capacitor

! ----------------------------------------------------------------------------
!
!                                INIT ROCK
!
! ----------------------------------------------------------------------------

!> Initialize the rock state to match the charge initialization for a capacitor.
subroutine lbe_elec_init_rock_capacitor(N)
  implicit none

  type(lbe_site), intent(inout) :: N(1-halo_extent:, 1-halo_extent:, 1-halo_extent:)

  integer :: i, j, k

  call log_msg_elec("  Initializing rock with option <capacitor>.")

  do i = 1, nx
    do j = 1, ny
      do k = 1, nz
        if ( ( k + ccoords(3)*nz .eq. 1*tnz/4 ) .or. ( k + ccoords(3)*nz .eq. 3*tnz/4 + 1 ) ) then
          N(i,j,k)%rock_state = rock_value
        end if
      end do
    end do
  end do

end subroutine lbe_elec_init_rock_capacitor

! ----------------------------------------------------------------------------
!
!                              EQUILIBRATION
!
! ----------------------------------------------------------------------------

!> At the start of the simulation the electrolytes have to be brought into balance.
subroutine elec_initial_equilibration(N)
  implicit none

  type(lbe_site), intent(inout) :: N(1-halo_extent:, 1-halo_extent:, 1-halo_extent:)

  real(kind=rk) :: Ex_tmp, Ey_tmp, Ez_tmp
  logical :: ch_equilibrated
  integer :: ch_loops

  character(len=12) :: dump_prefix

  ! Store the electric field values so we can disable it for now.
  Ex_tmp = Ex
  Ey_tmp = Ey
  Ez_tmp = Ez

  call start_timer(ti_elec_equilibration)

  call log_msg_elec_hdr("Initial equilibration")

  call log_msg_elec("Suppressing electric field for now.")
  ! Temporarily disable external electric field.
  Ex = 0.0_rk
  Ey = 0.0_rk
  Ez = 0.0_rk

  call log_msg_elec("Solving initial Poisson equation ...")
  ! Solve Poisson equation for current charge distribution.
  call solve_poisson(N)

  if ( E_solver_id == E_solver_fd ) then
    ! Calculate the electric field from the current potential.
    ! call log_msg_elec("Calling calc_E_fd(...)")
    call calc_E_fd(N)
  end if

  call elec_dump_all_asc(N, "post-sor")

  call log_msg_elec_ws("Equilibrating charges without field ...")

  ! Now allow the charges to move and equilibrate more.
  ch_equilibrated = .false.
  ch_loops = 0
  do while ( ( .not. ch_equilibrated .and. ( n_eq_noE < 0 ) ) .or. ( ch_loops < n_eq_noE ) )
    ch_loops = ch_loops + 1
    if ( n_show_eq > 0 ) then
      if ( ( mod(ch_loops, n_show_eq) .eq. 0 ) ) then
        write(msgstr,"('Performing nofield equilibration iteration ',I0,' .')") ch_loops
        call log_msg_elec(msgstr)
      end if
    end if
    call elec_timestep(N, ch_equilibrated, ch_loops)
    if ( n_dump_eq > 0 ) then
      if ( mod(ch_loops, n_dump_eq) .eq. 0 ) then
        write(dump_prefix,"('NOE',I8.8)") ch_loops
        call elec_dump_all_asc(N, trim(dump_prefix))
      end if
    end if
  end do
  if ( n_eq_noE < 0 ) then
    write(msgstr,"('  Charge diffusion without field converged with requested accuracy ',ES15.8,' after ',I0,' iterations.')") acc_fluxes, ch_loops
    call log_msg_elec(msgstr)
  else
    write(msgstr,"('  Charge diffusion without field stopped after ',I0,' requested iterations.')") ch_loops
    call log_msg_elec(msgstr)
  end if

  call elec_dump_all_asc(N, "post-noE")

  call log_msg_elec_ws("Equilibrating charges with external field enabled ... ")

  ! Restore the external electric field.
  Ex = Ex_tmp
  Ey = Ey_tmp
  Ez = Ez_tmp

  write(msgstr,"('Restored electric field to ',ES15.8,' ',ES15.8,' ',ES15.8,' .')") Ex, Ey, Ez
  call log_msg_elec(msgstr)

  ! Now allow the charges to move and equilibrate with the external electric field enabled.
  ch_equilibrated = .false.
  ch_loops = 0
  do while ( ( .not. ch_equilibrated .and. ( n_eq_E < 0 ) ) .or. ( ch_loops < n_eq_E ) )
    ch_loops = ch_loops + 1
    if ( n_show_eq > 0 ) then
      if( ( mod(ch_loops, n_show_eq) .eq. 0 ) ) then
        write(msgstr,"('Performing equilibration iteration ',I0,' .')") ch_loops
        call log_msg_elec(msgstr)
      end if
    end if
    call elec_timestep(N, ch_equilibrated, ch_loops)
    if ( n_dump_eq > 0 ) then
      if ( mod(ch_loops, n_dump_eq) .eq. 0 ) then
        write(dump_prefix,"('E',I8.8)") ch_loops
        call elec_dump_all_asc(N, trim(dump_prefix))
      end if
    end if
  end do
  if ( n_eq_E < 0 ) then
    write(msgstr,"('  Charge diffusion with field converged with requested accuracy ',ES15.8,' after ',I0,' iterations.')") acc_fluxes, ch_loops
    call log_msg_elec(msgstr)
  else
    write(msgstr,"('  Charge diffusion with field stopped after ',I0,' requested iterations.')") ch_loops
    call log_msg_elec(msgstr)
  end if

  call elec_dump_all_asc(N, "post-E")

  call log_msg_elec_ws("Finished initial equilibration.")

  call stop_timer(ti_elec_equilibration)

end subroutine elec_initial_equilibration

! ----------------------------------------------------------------------------
!
!                                 CLEANUP   
!
! ----------------------------------------------------------------------------

subroutine elec_shutdown()
  implicit none

  integer :: ierror

#ifdef P3M
  if ( poisson_solver_id == poisson_solver_p3m .and. local_eps ) then
    deallocate(E_prev, stat = ierror)
    call check_allocate(ierror, "Unable to deallocate E_prev.")
  end if
#endif

end subroutine elec_shutdown

#endif
end module lbe_elec_init_module

