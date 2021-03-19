#include "lbe.h"

!> code related to the advection and collision stages of each
!> timestep.
module lbe_collision_module
  use lbe_bdist_module
  use lbe_globals_module
  use lbe_helper_module, only: is_wall, is_fluid, is_colloid, CS_EOS
  use lbe_log_module
  use lbe_leesedwards_module
  use lbe_MRT
  use lbe_parms_module, only: BGK_id,collisiontype_id&
       &,MRT_id,SCMP,ZFOSdiag,acccoef,amass_b,amass_r,amass_s&
       &,bcsel,beta,boundary_cond,d_0,g_accn_fluid_r,g_accn_fluid_b,g_accn_fluid_s&
       &,g_accn,g_accn_x,g_accn_y,g_accn_max&
       &,g_accn_max_x,g_accn_max_y,g_accn_min,g_accn_min_x,g_accn_min_y,g_br&
       &,g_bs,g_rr,g_ss,g_wb,g_wr,kbT,Tcs,omega_d&
       &,omegabulk_b,omegabulk_r,omegabulk_s,psifunc,s02_r,s03_b,s03_r,s03_s,s05_b&
       &,s05_r,s05_s,s11_b,s11_r,s11_s,s17_b,s17_r,s17_s&
       &,tau_wb,tau_wr,inv_fluid,nx,ny,nz, &
       & omega_r, omega_b, omega_s, &
       & get_tau_r, get_omega_r, &
       & interpolation_scheme, interpolation_order, &
       & rock_colour_double, g_value_double, gw_double_wet, gw_wet_global
#ifndef SINGLEFLUID
  use lbe_parms_module, only: get_tau_b, get_omega_b, tau_b
#endif
#ifndef NOSURFACTANT
  use lbe_parms_module, only: get_tau_s, get_omega_s, tau_s
#endif
  use lbe_parallel_module, only: start,ccoords,tnz, Abend
#ifndef NOSURFACTANT
  use lbe_parallel_module, only: lbe_all_dipole_exchange
#endif
  use lbe_timer_module
  use lbe_types_module, only: lbe_site
#ifdef LOCALBC
  use lbe_types_module, only: bc_check
#endif 
#ifdef IBM_BINARYIBM
  use lsuperstruct_data_module, only: g_ri, g_bi, ibm_colour, noslipcoeff
  use lsuperstruct_interface_module
  use lbe_force_interface_module, only : add_force_to_n_r, add_force_to_n_b, add_force_to_all
#endif
#if AXISYM
  use lbe_axisym_module, only: lbe_axisym_correction_sc_forces, lbe_get_psi, lbe_gradient1, lbe_gradient2
#endif
  implicit none
  private

#ifndef SINGLEFLUID
  public lbe_advection,lbe_interactions,lbe_boundary_condition
#else
  public lbe_interactions,lbe_boundary_condition
#endif
  public lbe_calculate_sc_forces, md_calculate_sc_forces, zero_force_offset

  type(lbe_site), dimension(:,:,:), allocatable :: Nbuf

contains


  ! #ifndef SINGLEFLUID
  !> Performs advection on the particles in array \c N, writing the advected
  !> dipole moments for each lattice vector into the array \c d_adv.
  !>
  !> Also implements the bounce-back boundary conditions for rock sites.
  !>
  !> \todo Optimize the code (see FIXME below)!
  subroutine lbe_advection(N,d_adv)
    implicit none
    integer :: x, y, z
    integer :: xa, ya, za, xb, yb, zb
    integer :: bcz
    integer :: s,i
    type(lbe_site), dimension(0:,0:,0:),intent(inout) :: N
#ifndef NOSURFACTANT
    real*8, dimension(1:,1:,0:,0:,0:),intent(inout) :: d_adv
#else
    integer,intent(inout) :: d_adv
#endif
    real*8 :: init, finl, diff
    integer :: xx, yy, zz

    real*8 :: pred, pblue, psurf, qred, qblue, qsurf
    integer :: xp, yp, zp, xq, yq, zq, sp, sq
    integer :: fooflag=0

#ifdef LOCALBC
    type(bc_check),dimension(:,:,:),allocatable :: useddirections
#endif

#ifdef NOEDGESTEP
    integer, dimension(19,2) :: projection
    integer xaa, yaa, zaa
    integer xab, yab, zab
    integer xba, yba, zba
    integer xbb, ybb, zbb

    !projections of diagonal velocities:
    projection(1,:)  = (/ 19,19 /)
    projection(2,:)  = (/ 19,19 /)
    projection(3,:)  = (/ 19,19 /)
    projection(4,:)  = (/ 19,19 /)
    projection(5,:)  = (/ 19,19 /)
    projection(6,:)  = (/ 19,19 /)
    projection(7,:)  = (/ 1, 3 /)
    projection(8,:)  = (/ 1, 4 /)
    projection(9,:)  = (/ 1, 5 /)
    projection(10,:) = (/ 1, 6 /)
    projection(11,:) = (/ 2, 3 /)
    projection(12,:) = (/ 2, 4 /)
    projection(13,:) = (/ 2, 5 /)
    projection(14,:) = (/ 2, 6 /)
    projection(15,:) = (/ 3, 5 /)
    projection(16,:) = (/ 3, 6 /)
    projection(17,:) = (/ 4, 5 /)
    projection(18,:) = (/ 4, 6 /)
    projection(19,:) = (/ 19,19 /)
#endif

    ! FIXME there must be some way of improving this:
    ! preferably using one of the more terse
    ! array-transfer expressions provided by F90.

    if (.not. allocated(Nbuf)) allocate(Nbuf(0:nx+1,0:ny+1,-1:1))
#ifdef LOCALBC
    if (.not.allocated(useddirections)) allocate(useddirections(1:nx,1:ny,1:nz))
    do i=1,19
       useddirections(:,:,:)%n(i)=0
    end do
#endif
    Nbuf(0:nx+1,0:ny+1,-1) = N(0:nx+1,0:ny+1,0)
    Nbuf(0:nx+1,0:ny+1,0) = N(0:nx+1,0:ny+1,1)

    do z=1,nz
       Nbuf(0:nx+1,0:ny+1,1) = N(0:nx+1,0:ny+1,z+1)
       do y=1,ny
          do x=1,nx
             if( bcsel .ge. 0 )then
#ifdef SINGLEFLUID
#ifndef LOCALBC
                call lbe_boundary_condition(N,d_adv,x,y,z,Nbuf)
#endif
#ifdef LOCALBC
                call lbe_boundary_condition(N,d_adv,x,y,z,Nbuf,useddirections)
#endif
#else
#ifndef LOCALBC
                call lbe_boundary_condition(N,d_adv,x,y,z)
#endif
#ifdef LOCALBC
                call lbe_boundary_condition(N,d_adv,x,y,z,useddirections)
#endif
#endif
             end if
             ! Loop over each site in N.
             ! For each site, loop over each non-rest vector, and
             ! move the particles which wouuld travel *into* that site,
             ! via that vector.

             do s=1,nnonrest
                xb = x - cx(s)
                yb = y - cy(s)
                zb = z - cz(s)
                bcz = - cz(s)
#ifdef NOEDGESTEP
                xba = x - cx(projection(s,1))
                yba = y - cy(projection(s,1))
                zba = z - cz(projection(s,1))
                xbb = x - cx(projection(s,2))
                ybb = y - cy(projection(s,2))
                zbb = z - cz(projection(s,2))
#endif

! Remember not to advect on top of the
                ! bounced vectors..
                if ( bcsel .ne. 2 ) then
                   if ( N(x,y,z)%rock_state == 0. ) then
                      if ( N(xb,yb,zb)%rock_state == 0.) then
#ifdef NOEDGESTEP
                         if( .not. ( s .gt. 6 .and. &
                              N(xba,yba,zba)%rock_state .ne. 0. .and. &
                              N(xbb,ybb,zbb)%rock_state .ne. 0. ) ) then
#endif
                            N(x,y,z)%n_r(s) = Nbuf(xb,yb,bcz)%n_r(s)

#ifndef SINGLEFLUID
                            N(x,y,z)%n_b(s) = Nbuf(xb,yb,bcz)%n_b(s)

#endif
#ifndef NOSURFACTANT
                            N(x,y,z)%n_s(s) = Nbuf(xb,yb,bcz)%n_s(s)
                            d_adv(:,s,x,y,z) = Nbuf(xb,yb,bcz)%d
#endif
#ifdef NOEDGESTEP
                         endif
#endif
                      endif
                   endif
                else
                   N(x,y,z)%n_r(s) = Nbuf(xb,yb,bcz)%n_r(s)
#ifndef SINGLEFLUID
                   N(x,y,z)%n_b(s) = Nbuf(xb,yb,bcz)%n_b(s)
#endif
#ifndef NOSURFACTANT
                   N(x,y,z)%n_s(s) = Nbuf(xb,yb,bcz)%n_s(s)
                   d_adv(:,s,x,y,z) = Nbuf(xb,yb,bcz)%d
#endif
                endif
             end do  !s
#ifndef NOSURFACTANT
             d_adv(:,restvec,x,y,z) = Nbuf(x,y,0)%d
             ! Rest particles stay the same.
#endif
          end do     !x
       end do        !y
       Nbuf(0:nx+1,0:ny+1,-1)=Nbuf(0:nx+1,0:ny+1,0)
       Nbuf(0:nx+1,0:ny+1,0)=Nbuf(0:nx+1,0:ny+1,1)
    end do           !z

#ifndef NOSURFACTANT
    !Ensure d_adv halo is defined, or there can be random problems
    d_adv(:,:,0,:,:) = 0.d0
    d_adv(:,:,nx+1,:,:) = 0.d0
    d_adv(:,:,:,0,:) = 0.d0
    d_adv(:,:,:,ny+1,:) = 0.d0
    d_adv(:,:,:,:,0) = 0.d0
    d_adv(:,:,:,:,nz+1) = 0.d0
#endif

#ifdef LOCALBC
    deallocate(useddirections)
#endif
  end subroutine lbe_advection

  !> Collision step.
  !>
  !> Takes a subdomain array \c N / \c whole_N and an array of
  !> advected dipole moments \c d_adv, and collides and redistributes
  !> the particles at each site.
  !>
  !> If this is ever run in a situation where a Fortran 95 compiler
  !> is available, it could very possibly benefit from being rephrased
  !> in terms of the \c forall statement in F95.
  !>
  !> \param[in] whole_N local lattice chunk with halo of extent \c
  !> halo_extent
  subroutine lbe_interactions(N,whole_N,d_adv)
    type(lbe_site), dimension(0:,0:,0:) :: N
    type(lbe_site),intent(in)&
         & :: whole_N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
#ifndef NOSURFACTANT
    real*8, dimension(1:,1:,0:,0:,0:) :: d_adv
#else
    integer :: d_adv
#endif

    real(kind=rk) :: g_accn_local, g_accn_local_x, g_accn_local_y
    integer :: x, y, z, s
    integer :: tz, tx, ty
    integer :: xa, ya, za, ss
    logical :: minz, maxz

    real*8, dimension(3) :: h_c
    real*8 :: hc_m ! Magnitude of h_c
    ! h_s is a vector involving the surfactant part of the Hamiltonian.
    real*8, dimension(3) :: h_s

    real*8 :: hs_dot

    ! h is the sum of h_s and h_c
    real*8, dimension(3) :: h
    real*8 :: h_m  ! Magnitude of h, times beta

    ! sd_eq is the equilibrium surfactant density.
    real*8, dimension(3) :: sd_eq

    ! These three are used in some more of Hudong's code
    ! which I don't quite follow.
    real*8, dimension(3) :: h_bar
    real*8 :: sd_mag
    real*8 :: exh

    ! sN is the surfactant density at a site.
    real*8 :: sN,sN_tmp
    real*8, dimension(3) :: sd_tmp

    real*8 :: psi_b, psi_r

    real*8, dimension(3) :: pas = (/ 0.d0, 0.d0, 0.d0 /)
    real*8, dimension(19) :: tmp_Nn_r,tmp_Nn_b,tmp_Nn_s

    ! Nett forces
    real*8, dimension(nx,ny,nz,3) :: f_b
    real*8, dimension(nx,ny,nz,3) :: f_r
    real*8, dimension(nx,ny,nz,3) :: f_s

    ! Pointwise forces (used only for IBM plug-in)
#ifdef IBM_BINARYIBM
    real(kind=rk), dimension(3) :: f_loc_r
    real(kind=rk), dimension(3) :: f_loc_b
    real(kind=rk), dimension(3) :: f_loc_i
#endif

    ! Velocities
    real*8, dimension(3) :: p_r, p_b, p_s, p_tilde, u_tilde

    ! Number densities
    real*8 :: bN_s, rN_s, sN_s

    ! Weighted masses
    real*8 :: rho_b, rho_r, rho_s, rho_tilde

    ! Total fluid density
    real*8 :: rho_tot

    ! Acceleration due to force acting on total fluid
    real*8, dimension(3) :: acc_tot

    ! Boltzmann distributions
    real*8, dimension(nvecs) :: boltz_b, boltz_r, boltz_s

    ! Equilibrium distributions at a site
    real*8, dimension(3) :: bN_eq(nvecs), rN_eq(nvecs), sN_eq(nvecs)
    real*8 :: pon_r, pon_s, pon_b
    real*8, dimension(3) :: du_r, du_s, du_b
    real*8, dimension(3) :: dd_r, dd_s, dd_b

!!!!!!!!!!!!!!!!MRT VARIABLES START
    real*8, dimension(19) :: S_hat_r
    real*8, dimension(19) :: S_hat_b
    real*8, dimension(19) :: S_hat_s


    ! temporary variable for force calculation
    real*8 :: forceterm

#ifdef AXISYM
    real*8,dimension(:,:,:,:),allocatable :: psi_mat

    real*8,dimension(:,:,:),allocatable :: dy_psi_r,dyz_psi_r,dyy_psi_r
#ifndef SINGLEFLUID
    real*8,dimension(:,:,:),allocatable :: dy_psi_b,dyz_psi_b,dyy_psi_b
#endif
#ifndef NOSURFACTANT
    real*8,dimension(:,:,:),allocatable :: dy_psi_s,dyz_psi_s,dyy_psi_s
#endif
#endif



#ifdef AXISYM
    if (.not. allocated(psi_mat)) allocate(psi_mat(1:nx,1:ny,1:nz,1:3))

    if (.not. allocated(dy_psi_r)) allocate(dy_psi_r(1:nx,1:ny,1:nz))
    if (.not. allocated(dyz_psi_r)) allocate(dyz_psi_r(1:nx,1:ny,1:nz))
    if (.not. allocated(dyy_psi_r)) allocate(dyy_psi_r(1:nx,1:ny,1:nz))
#ifndef SINGLEFLUID
    if (.not. allocated(dy_psi_b)) allocate(dy_psi_b(1:nx,1:ny,1:nz))
    if (.not. allocated(dyz_psi_b)) allocate(dyz_psi_b(1:nx,1:ny,1:nz))
    if (.not. allocated(dyy_psi_b)) allocate(dyy_psi_b(1:nx,1:ny,1:nz))
#endif
#ifndef NOSURFACTANT
    if (.not. allocated(dy_psi_s)) allocate(dy_psi_s(1:nx,1:ny,1:nz))
    if (.not. allocated(dyz_psi_s)) allocate(dyz_psi_s(1:nx,1:ny,1:nz))
    if (.not. allocated(dyy_psi_s)) allocate(dyy_psi_s(1:nx,1:ny,1:nz))
#endif
#endif

    call start_timer(ti_intf)

    ! Diagonal collision matrix S_hat, given in Phil Trans R Soc Lond A 360, pp 437-451,
    ! transformed to match LB3D lattice definition,
    ! further simplifiable by choosing s10 = s14.
    ! The kinematic viscosity then is given by ( 1 / 3 ) * ( ( 1 / s10 ) - ( 1 / 2 ) ).
    ! The bulk viscosity is given by ( ( 5 - 9 * c_s^2 ) / 9 ) * ( ( 1 / s2 ) - ( 1 / 2 ) ) ; c_s^2 = 1 / 3.
    ! TODO:
    ! - What the hell is this definition doing here? If it is constant, it should be put somewhere else.
    !   If it is non-constant, it must be computed after the viscosity is set.
    !   In any case, it should not be here. (Added by Timm (2012-02-24))

    if ( collisiontype_id == MRT_id) then

!      write(*,"('Inside lbe_interactions')")
!      S_hat_r(1) = 0.d0  ! rho           0.00d0
!      S_hat_r(2) = s02_r !omegabulk_r    ! e       1.19d0
!      S_hat_r(3) = s03_r    ! epsilon       1.40d0
!      S_hat_r(4) = 0.d0  ! j_x           0.00d0
!      S_hat_r(5) = s05_r    ! q_x           1.20d0
!      S_hat_r(6) = 0.d0  ! j_y           0.00d0
!      S_hat_r(7) = s05_r    ! q_y           1.20d0
!      S_hat_r(8) = 0.d0  ! j_z           0.00d0
!      S_hat_r(9) = s05_r    ! q_z           1.20d0
!      S_hat_r(10)= omega_r    ! 3p_xx         1.00d0
!      S_hat_r(11)= s11_r    ! 3pi_xx        1.40d0
!      S_hat_r(12)= omega_r    ! p_ww          1.00d0
!      S_hat_r(13)= s11_r    ! pi_ww         1.40d0
!      S_hat_r(14)= omega_r    ! p_xy          1.00d0
!      S_hat_r(15)= omega_r    ! p_yz          1.00d0
!      S_hat_r(16)= omega_r    ! p_xz          1.00d0
!      S_hat_r(17)= s17_r    ! m_y           1.98d0
!      S_hat_r(18)= s17_r    ! m_x           1.98d0
!      S_hat_r(19)= s17_r    ! m_z           1.98d0


      S_hat_r(1) = 1.d0  ! rho           0.00d0
      S_hat_r(2) = omegabulk_r    ! e       1.19d0
      S_hat_r(3) = s03_r    ! epsilon       1.40d0
      S_hat_r(4) = 1.d0  ! j_x           0.00d0
      S_hat_r(5) = s05_r    ! q_x           1.20d0
      S_hat_r(6) = 1.d0  ! j_y           0.00d0
      S_hat_r(7) = s05_r    ! q_y           1.20d0
      S_hat_r(8) = 1.d0  ! j_z           0.00d0
      S_hat_r(9) = s05_r    ! q_z           1.20d0
      S_hat_r(10)= omega_r    ! 3p_xx         1.00d0
      S_hat_r(11)= s11_r    ! 3pi_xx        1.40d0
      S_hat_r(12)= omega_r    ! p_ww          1.00d0
      S_hat_r(13)= s11_r    ! pi_ww         1.40d0
      S_hat_r(14)= omega_r    ! p_xy          1.00d0
      S_hat_r(15)= omega_r    ! p_yz          1.00d0
      S_hat_r(16)= omega_r    ! p_xz          1.00d0
      S_hat_r(17)= s17_r    ! m_y           1.98d0
      S_hat_r(18)= s17_r    ! m_x           1.98d0
      S_hat_r(19)= s17_r    ! m_z           1.98d0

      S_hat_b(1) = 1.d0  ! rho           0.00d0
      S_hat_b(2) = omegabulk_b    ! e             1.19d0
      S_hat_b(3) = s03_b    ! epsilon       1.40d0
      S_hat_b(4) = 1.d0  ! j_x           0.00d0
      S_hat_b(5) = s05_b    ! q_x           1.20d0
      S_hat_b(6) = 1.d0  ! j_y           0.00d0
      S_hat_b(7) = s05_b    ! q_y           1.20d0
      S_hat_b(8) = 1.d0  ! j_z           0.00d0
      S_hat_b(9) = s05_b    ! q_z           1.20d0
      S_hat_b(10)= omega_b    ! 3p_xx         1.00d0
      S_hat_b(11)= s11_b    ! 3pi_xx        1.40d0
      S_hat_b(12)= omega_b    ! p_ww          1.00d0
      S_hat_b(13)= s11_b    ! pi_ww         1.40d0
      S_hat_b(14)= omega_b    ! p_xy          1.00d0
      S_hat_b(15)= omega_b    ! p_yz          1.00d0
      S_hat_b(16)= omega_b    ! p_xz          1.00d0
      S_hat_b(17)= s17_b    ! m_y           1.98d0
      S_hat_b(18)= s17_b    ! m_x           1.98d0
      S_hat_b(19)= s17_b    ! m_z           1.98d0

      S_hat_s(1) = 1.d0  ! rho           0.00d0
      S_hat_s(2) = omegabulk_s    ! e             1.19d0
      S_hat_s(3) = s03_s    ! epsilon       1.40d0
      S_hat_s(4) = 1.d0  ! j_x           0.00d0
      S_hat_s(5) = s05_s    ! q_x           1.20d0
      S_hat_s(6) = 1.d0  ! j_y           0.00d0
      S_hat_s(7) = s05_s    ! q_y           1.20d0
      S_hat_s(8) = 1.d0  ! j_z           0.00d0
      S_hat_s(9) = s05_s    ! q_z           1.20d0
      S_hat_s(10)= omega_s    ! 3p_xx         1.00d0
      S_hat_s(11)= s11_s    ! 3pi_xx        1.40d0
      S_hat_s(12)= omega_s    ! p_ww          1.00d0
      S_hat_s(13)= s11_s    ! pi_ww         1.40d0
      S_hat_s(14)= omega_s    ! p_xy          1.00d0
      S_hat_s(15)= omega_s    ! p_yz          1.00d0
      S_hat_s(16)= omega_s    ! p_xz          1.00d0
      S_hat_s(17)= s17_s    ! m_y           1.98d0
      S_hat_s(18)= s17_s    ! m_x           1.98d0
      S_hat_s(19)= s17_s    ! m_z           1.98d0
    end if

!!!!!!!!!!!!!!!!MRT VARIABLES END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Surfactant computations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifndef NOSURFACTANT
    ! Loop over all lattice sites (icky)

    do z=1,nz
       do y=1,ny
          do x=1,nx
             h_c = 0.d0
             h_s = 0.d0
             sN = 0.d0

             ! Loop over all vectors, and calculate h_s:
             do s=1,nvecs
                sN = sN + N(x,y,z)%n_s(s)*g(s)
                h_s = h_s + g(s)*N(x,y,z)%n_s(s)*d_adv(:,s,x,y,z)
             end do ! s
             ! Use a normalizing factor, and ensure that it is nonzero.
             ! FNORD: I'm not sure what the point of this is, exactly.

             sN_tmp = max(sN,dble(10.e-5))
             sd_tmp = h_s/sN_tmp

             ! Now loop over only the nonrest vectors.

             do s=1,nnonrest
                psi_b = 0.d0
                psi_r = 0.d0
                pas = 0.d0

                xa = x + cx(s)
                ya = y + cy(s)
                za = z + cz(s)

                do ss=1,nvecs
                   psi_b = psi_b + N(xa,ya,za)%n_b(ss)*g(ss)
                   psi_r = psi_r + N(xa,ya,za)%n_r(ss)*g(ss)
                   pas = pas + N(xa,ya,za)%n_s(ss)*d_adv(:,ss,xa,ya,za)*g(ss)
                end do
                ! psi_b, psi_r give the blue, red densities at (xa,ya,za).
                ! Calculate magnitude of h_c
                hc_m = - sign(dble(1.0),g_bs)*(psi_b - psi_r)*g(s)

                ! Make it into a vector.
                h_c = h_c + hc_m * c(s,:)

                !hs_dot = D/2.*( cx(s)*pas(1) + cy(s)*pas(2) + cz(s)*pas(3) )
                hs_dot = D/2.d0*sum(c(s,:)*pas)

                ! Make h_s into a vector
                h_s = h_s + (pas - c(s,:)*hs_dot) * g(s)
             end do ! s over non-rest vectors

             h = h_c + h_s

             ! h_m = beta * magnitude of h
             h_m = beta*sqrt(sum(h*h))

             sd_eq = 1.0_8/3.0_8*beta*h

             ! What's going on here? Looks like a clipping routine
             ! to avoid divide-by-zeroes, but I'm not sure.
             if (h_m > 10.d-3) then
                exh = 1.d0/tanh(h_m)
                h_bar = beta*h/h_m
                sd_mag = exh - 1.d0/h_m
                sd_eq = h_bar*sd_mag
             endif

             ! This line was added in a hurry to allow d_0 to be varied
             sd_eq = sd_eq*d_0

             N(x,y,z)%d = omega_d*sd_eq + ( 1.d0 - omega_d )*sd_tmp
          end do ! x
       end do ! y
    end do ! z

    if ( ( boundary_cond .ne. 0 ) .and. &
         ( inv_fluid .eq. 5 .or. inv_fluid .eq. 6 ) ) then
       call le_all_dipole_exchange(d_adv,N,whole_N)
    else
       call lbe_all_dipole_exchange(d_adv,N)
    endif
#endif

    ! The equilibrium values and Hamiltonians have now been calculated.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Compute Shan-Chen forces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Revised by Timm (2012-02-24).
    ! The forces are first stored in temporary arrays f_r, f_b, and f_s (depending on the number of species).
    ! Afterwards, the forces are added to lbe_force where they will survive the entire time step.
    ! This way, the forces are available for velocity correction without recomputation.
    do z = 1, nz
      do y = 1, ny
        do x = 1, nx
          ! Since the forces for collision are computed here, only non-rock sites are considered.
          ! If Shan-Chen forces are required on rock sites, e.g., for a force or torque balance,
          ! the corresponding module has to call lbe_calculate_sc_forces on its own.

          if(N(x, y, z)%rock_state == 0.d0) then
            ! The following code is executed only in a very special case:
            ! Using IBM particles with IBM index field and exactly two fluid components.
            ! This should be very uninteresting for most users.
#ifdef IBM_BINARYIBM
            call ibm_calculate_sc_forces(N, x, y, z, f_loc_r, f_loc_b, f_loc_i)
            call add_force_to_n_r(f_loc_r(1), f_loc_r(2), f_loc_r(3), x, y, z)
            call add_force_to_n_b(f_loc_b(1), f_loc_b(2), f_loc_b(3), x, y, z)
            call add_force_to_all(f_loc_i(1), f_loc_i(2), f_loc_i(3), x, y, z)
#else
	    ! This is the standard code everybody is used to:
#ifdef SINGLEFLUID
            call lbe_calculate_sc_forces(N, x, y, z, f_r)
            lbe_force(1:3, 1, x, y, z) = lbe_force(1:3, 1, x, y, z) + f_r(x, y, z, 1:3)
#else
#ifdef NOSURFACTANT
            call lbe_calculate_sc_forces(N, x, y, z, f_b, f_r)
            lbe_force(1:3, 1, x, y, z) = lbe_force(1:3, 1, x, y, z) + f_r(x, y, z, 1:3)
            lbe_force(1:3, 2, x, y, z) = lbe_force(1:3, 2, x, y, z) + f_b(x, y, z, 1:3)
#else
            call lbe_calculate_sc_forces(N, x, y, z, f_b, f_r, f_s)
            lbe_force(1:3, 1, x, y, z) = lbe_force(1:3, 1, x, y, z) + f_r(x, y, z, 1:3)
            lbe_force(1:3, 2, x, y, z) = lbe_force(1:3, 2, x, y, z) + f_b(x, y, z, 1:3)
            lbe_force(1:3, 3, x, y, z) = lbe_force(1:3, 3, x, y, z) + f_s(x, y, z, 1:3)
#endif 
#endif
#endif 

         endif
        end do
      end do
    end do


! sudhir
#ifdef AXISYM

    call lbe_get_psi(N,psi_mat)
!    psi_mat(1,:,:,:) = psi_mat(3,:,:,:)
!    psi_mat(4,:,:,:) = psi_mat(2,:,:,:)

#ifdef SINGLEFLUID
    call lbe_gradient1(psi_mat,dy_psi_r)
    call lbe_gradient2(dy_psi_r,dyz_psi_r,dyy_psi_r)
#else
#ifdef NOSURFACTANT
    call lbe_gradient1(psi_mat,dy_psi_r,dy_psi_b)
    call lbe_gradient2(dy_psi_r,dyz_psi_r,dyy_psi_r,&
         & dy_psi_b,dyz_psi_b,dyy_psi_b)
#else
    call lbe_gradient1(psi_mat,dy_psi_r,dy_psi_b,dy_psi_s)
    call lbe_gradient2(dy_psi_r,dyz_psi_r,dyy_psi_r,&
         & dy_psi_b,dyz_psi_b,dyy_psi_b,&
         & dy_psi_s,dyz_psi_s,dyy_psi_s)
#endif
#endif
!!$

    do z = 1,nz
      do y = 1,ny
        do x = 2,nx-1
           if(N(x, y, z)%rock_state == 0.d0) then

#ifdef SINGLEFLUID
              call lbe_axisym_correction_sc_forces(psi_mat,&
                   & dy_psi_r, dyz_psi_r, dyy_psi_r,&
                   & x, y, z, f_r)
              f_r(x,y,z,1) = 0.0_rk
              lbe_force(1:3, 1, x, y, z) = lbe_force(1:3, 1, x, y, z) + f_r(x, y, z, 1:3)
#else
#ifdef NOSURFACTANT
              call lbe_axisym_correction_sc_forces(psi_mat,&
                   & dy_psi_r, dyz_psi_r, dyy_psi_r,&
                   & dy_psi_b, dyz_psi_b, dyy_psi_b,&
                   & x, y, z, f_b, f_r)
              f_r(x,y,z,1) = 0.0_rk
              f_b(x,y,z,1) = 0.0_rk
              lbe_force(2:3, 1, x, y, z) = lbe_force(2:3, 1, x, y, z) + f_r(x, y, z, 2:3)
              lbe_force(2:3, 2, x, y, z) = lbe_force(2:3, 2, x, y, z) + f_b(x, y, z, 2:3)
#else
              call lbe_axisym_correction_sc_forces(psi_mat,&
                   & dy_psi_r, dyz_psi_r, dyy_psi_r,&
                   & dy_psi_b, dyz_psi_b, dyy_psi_b,&
                   & dy_psi_s, dyz_psi_s, dyy_psi_s,&
                   & x, y, z, f_b, f_r, f_s)
              lbe_force(1:3, 1, x, y, z) = lbe_force(1:3, 1, x, y, z) + f_r(x, y, z, 1:3)
              lbe_force(1:3, 2, x, y, z) = lbe_force(1:3, 2, x, y, z) + f_b(x, y, z, 1:3)
              lbe_force(1:3, 3, x, y, z) = lbe_force(1:3, 3, x, y, z) + f_s(x, y, z, 1:3)
#endif
#endif

              if((abs(f_r(x,y,z,1))) .ge. 1e-10) then
                 print *,x,y,z,dyz_psi_r(x,y,z),dyy_psi_r(x,y,z)
                 print *,x,y,z,dy_psi_r(x,y,z)
                 print *,'-------------------------------------'
                 stop
              end if
              
           end if

        end do
      end do
    end do
#endif ! AXISYM

#ifdef IBM_BINARYIBM
!     call compute_multicomponent_forces(N)
#endif

    call stop_timer(ti_intf)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Collision stage
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call start_timer(ti_intc)

    ! Loop over all lattice sites.
    do z = 1,nz
      do y = 1,ny
        do x = 1,nx
          ! Collision is only performed on fluid sites.
          if(N(x, y, z)%rock_state == 0.d0) then

            ! ================================
            ! COMPUTE DENSITIES AND VELOCITIES
            ! ================================

            ! Compute species densities and velocities:
            ! - rN_s = particle density of red species etc.
            ! - p_r = velocity of red species etc.
            ! TODO
            ! - Wouldn't it be simpler to use Fortran's sum function here?
            p_b = 0.d0
            p_r = 0.d0
            p_s = 0.d0
            bN_s = 0.d0
            rN_s = 0.d0
            sN_s = 0.d0

            do s = 1, nvecs
              p_r = p_r + N(x, y, z)%n_r(s) * g(s) * c(s, :)
              rN_s = rN_s + N(x, y, z)%n_r(s) * g(s)
#ifndef SINGLEFLUID
              p_b = p_b + N(x, y, z)%n_b(s) * g(s) * c(s, :)
              bN_s = bN_s + N(x, y, z)%n_b(s) * g(s)
#endif
#ifndef NOSURFACTANT
              p_s = p_s + N(x, y, z)%n_s(s) * g(s) * c(s, :)
              sN_s = sN_s + N(x, y, z)%n_s(s) * g(s)
#endif
            end do

            ! Calculate weighted total momentum.
            p_tilde = amass_r * get_omega_r(N, x, y, z) * p_r
#ifndef SINGLEFLUID
            p_tilde = p_tilde + amass_b * get_omega_b(N, x, y, z) * p_b
#endif
#ifndef NOSURFACTANT
            p_tilde = p_tilde + amass_s * get_omega_s(N, x, y, z) * p_s
#endif

            ! Calculate mass densities of species.
            rho_r = amass_r * rN_s
            rho_b = amass_b * bN_s
            rho_s = amass_s * sN_s

            ! Calculate relaxation time averaged mass density.
            rho_tilde = rho_r * get_omega_r(N, x, y, z)
#ifndef SINGLEFLUID
            rho_tilde = rho_tilde + rho_b * get_omega_b(N, x, y, z)
#endif
#ifndef NOSURFACTANT
            rho_tilde = rho_tilde + rho_s * get_omega_s(N, x, y, z)
#endif

            ! Calculate total fluid density.
            rho_tot = rho_r + rho_b + rho_s

            ! Calculate relaxation time averaged velocity.
            ! Use a cutoff to avoid division by very small numbers.
            u_tilde = p_tilde / max(rho_tilde, 10.d-9)

            ! Calculate inverse particle densities.
            if(psifunc .eq. 0) then
              ! This expression is only for the original funny clipped version of the code.
              ! Evaluates to 1 unless densities go over 1.
              pon_r = min(1.d0 / max(rN_s, 10.d-9), 1.d0)
              pon_b = min(1.d0 / max(bN_s, 10.d-9), 1.d0)
              pon_s = min(1.d0 / max(sN_s, 10.d-9), 1.d0)
            else
              ! Sane code that does what the literature says.
              pon_r = 1.d0 / max(rN_s, 10.d-9)
              pon_b = 1.d0 / max(bN_s, 10.d-9)
              pon_s = 1.d0 / max(sN_s, 10.d-9)
            end if

            ! ==========================
            ! GRAVITATIONAL ACCELERATION
            ! ==========================

            ! Add the force due to g_accn to the 0-component of lbe_force
            ! since it acts on the total fluid rather than on the species.
            ! The force equals the acceleration times the total fluid density.
            ! NOTES:
            ! - As in the previous revision of the code, g_accn etc. act only on fluid sites, not on rock sites.
            ! TODO:
            ! - This functionality should be moved to some forcing module.

            ! Compute global position of the local site.
            tx = x + ccoords(1) * nx
            ty = y + ccoords(2) * ny
            tz = z + ccoords(3) * nz

            ! Check whether this site is within the externally driven chunk.
            if((tz .ge. g_accn_min) .and. (tz .le. g_accn_max) .and. &
                & (tx .ge. g_accn_min_x) .and. (tx .le. g_accn_max_x) .and. &
                & (ty .ge. g_accn_min_y) .and. (ty .le. g_accn_max_y)) then
              g_accn_local_x = g_accn_x
              g_accn_local_y = g_accn_y
              g_accn_local = g_accn
            else
              g_accn_local = 0.d0
              g_accn_local_y = 0.d0
              g_accn_local_x = 0.d0
            endif

            ! Consider inv_fluid settings and add the acceleration as a force to lbe_force.
            if((inv_fluid == 9) .or. (inv_fluid == 20) .or. (inv_fluid == 21)) then
              minz = (start(3) == 1)
              maxz = (start(3) >= (tnz - nz))

              if((minz) .and. (z .eq. 1)) then
                lbe_force(1, 0, x, y, z) = lbe_force(1, 0, x, y, z) + g_accn_local_x * rho_tot
                lbe_force(2, 0, x, y, z) = lbe_force(2, 0, x, y, z) + g_accn_local_y * rho_tot
                lbe_force(3, 0, x, y, z) = lbe_force(3, 0, x, y, z) + g_accn_local * rho_tot
              endif

              ! Only add force at z = tnz if inv_fluid = 9 (I know it is ugly)
              if((inv_fluid == 9) .and. (maxz) .and. (z .eq. tnz)) then
                lbe_force(1, 0, x, y, z) = lbe_force(1, 0, x, y, z) + g_accn_local_x * rho_tot
                lbe_force(2, 0, x, y, z) = lbe_force(2, 0, x, y, z) + g_accn_local_y * rho_tot
                lbe_force(3, 0, x, y, z) = lbe_force(3, 0, x, y, z) + g_accn_local * rho_tot
              endif
           else
               ! In all other cases, the gravitational acceleration is added to fluid sites without additional check.
#ifdef SINGLEFLUID
              lbe_force(1, 1, x, y, z) = lbe_force(1, 1, x, y, z) + g_accn_local_x * rho_r
              lbe_force(2, 1, x, y, z) = lbe_force(2, 1, x, y, z) + g_accn_local_y * rho_r
              lbe_force(3, 1, x, y, z) = lbe_force(3, 1, x, y, z) + g_accn_local * rho_r
#else
              if (g_accn_fluid_r) then
                 lbe_force(1, 1, x, y, z) = lbe_force(1, 1, x, y, z) + g_accn_local_x * rho_r
                 lbe_force(2, 1, x, y, z) = lbe_force(2, 1, x, y, z) + g_accn_local_y * rho_r
                 lbe_force(3, 1, x, y, z) = lbe_force(3, 1, x, y, z) + g_accn_local * rho_r
              endif
              if (g_accn_fluid_b) then
                 lbe_force(1, 2, x, y, z) = lbe_force(1, 2, x, y, z) + g_accn_local_x * rho_b
                 lbe_force(2, 2, x, y, z) = lbe_force(2, 2, x, y, z) + g_accn_local_y * rho_b
                 lbe_force(3, 2, x, y, z) = lbe_force(3, 2, x, y, z) + g_accn_local * rho_b
              endif
#endif
#ifndef NOSURFACTANT        
              if (g_accn_fluid_s) then
                 lbe_force(1, 3, x, y, z) = lbe_force(1, 2, x, y, z) + g_accn_local_x * rho_s
                 lbe_force(2, 3, x, y, z) = lbe_force(2, 2, x, y, z) + g_accn_local_y * rho_s
                 lbe_force(3, 3, x, y, z) = lbe_force(3, 2, x, y, z) + g_accn_local * rho_s
              endif
#endif
            endif

!!!Test:Store distribution function for Velocity Output before Collision
!#ifdef COMMON_VEL_FIX
#ifndef OLD_VEL
    N(x, y, z)%n_r_pre = N(x, y, z)%n_r
#ifndef SINGLEFLUID
	 N(x, y, z)%n_b_pre = N(x, y, z)%n_b
#ifndef NOSURFACTANT
	N(x, y, z)%n_s_pre = N(x, y, z)%n_s
#endif
#endif
#endif
            ! =================================
            ! PREPARE FLUID FORCE FOR COLLISION
            ! =================================

            ! All force contributions have been added to lbe_force now.
            ! The corresponding fluid acceleration can be computed.
            acc_tot(:) = lbe_force(:, 0, x, y, z) / rho_tot

            ! ===========================
            ! COLLISION FOR RED COMPONENT
            ! ===========================

            ! Compute velocity shift due to forces for the red component.
            du_r(:) = get_tau_r(N, x, y, z) * (pon_r * lbe_force(:, 1, x, y, z) / amass_r +  acc_tot)

            ! The dd_r array (Phys. Rev. E 65,(2002): Guo) used in the previous revision will be deprecated soon
            ! since all forces will be coupled in the same way in the future (user's choice: either Guo or Shan-Chen).
            ! For a temporary measure, dd_r is set to zero.
            dd_r(:) = 0.d0

            ! Perform collision.
            ! TODO
            ! - redefine boltz_dist by taking out dd_r later.
            ! - check: Is rN_eq and rN_s right input for equilibrium distribution and density
            if ( collisiontype_id == MRT_id ) then
              call mrt_dist(u_tilde + du_r - dd_r, dd_r / (get_tau_r(N, x, y, z) - 0.5d0), N(x, y, z)%n_r, S_hat_r)
            else if ( collisiontype_id == BGK_id) then
              call boltz_dist(u_tilde(1), u_tilde(2), u_tilde(3), du_r(1), du_r(2), du_r(3), dd_r(1), dd_r(2), dd_r(3), boltz_r)
              rN_eq = rN_s * boltz_r
              N(x, y, z)%n_r = N(x, y, z)%n_r - get_omega_r(N, x, y, z) * (N(x, y, z)%n_r - rN_eq)
            end if

            ! ============================
            ! COLLISION FOR BLUE COMPONENT
            ! ============================

#ifndef SINGLEFLUID
            ! Compute velocity shift due to forces for the blue component.
            du_b(:) = get_tau_b(N, x, y, z) * (pon_b * lbe_force(:, 2, x, y, z) / amass_b +  acc_tot)

            ! The dd_b array (Phys. Rev. E 65,(2002): Guo) used in the previous revision will be deprecated soon
            ! since all forces will be coupled in the same way in the future (user's choice: either Guo or Shan-Chen).
            ! For a temporary measure, dd_b is set to zero.
            dd_b(:) = 0.d0

            ! Perform collision.
            ! TODO
            ! - redefine boltz_dist by taking out dd_b later.
            ! - check: Is bN_eq and bN_s right input for equilibrium distribution and density
            if ( collisiontype_id == MRT_id) then
              call mrt_dist(u_tilde + du_b - dd_b, dd_b / (tau_b - 0.5d0), N(x, y, z)%n_b, S_hat_b)
            else if( collisiontype_id == BGK_id) then
              call boltz_dist(u_tilde(1), u_tilde(2), u_tilde(3), du_b(1), du_b(2), du_b(3), dd_b(1), dd_b(2), dd_b(3), boltz_b)
              bN_eq = bN_s * boltz_b
              N(x, y, z)%n_b = N(x, y, z)%n_b - get_omega_b(N, x, y, z) * (N(x, y, z)%n_b - bN_eq)
            end if
#endif

            ! ========================
            ! COLLISION FOR SURFACTANT
            ! ========================

#ifndef NOSURFACTANT
            ! Compute velocity shift due to forces for the surfactant.
            du_s(:) = get_tau_s(N, x, y, z) * (pon_s * lbe_force(:, 3, x, y, z) / amass_s +  acc_tot)

            ! The dd_s array (Phys. Rev. E 65,(2002): Guo) used in the previous revision will be deprecated soon
            ! since all forces will be coupled in the same way in the future (user's choice: either Guo or Shan-Chen).
            ! For a temporary measure, dd_s is set to zero.
            dd_s(:) = 0.d0

            ! Perform collision.
            ! TODO
            ! - redefine boltz_dist by taking out dd_s later.
            ! - check: Is sN_eq and sN_s right input for equilibrium distribution and density
            if( collisiontype_id == MRT_id ) then
              call mrt_dist(u_tilde + du_s - dd_s, dd_s / (get_tau_s(N, x, y, z) - 0.5d0), N(x, y, z)%n_s, S_hat_s)
            else if( collisiontype_id == BGK_id) then
              call boltz_dist(u_tilde(1),u_tilde(2),u_tilde(3), du_s(1), du_s(2), du_s(3), dd_s(1), dd_s(2), dd_s(3), boltz_s)
              sN_eq = sN_s * boltz_s
              N(x, y, z)%n_s = N(x, y, z)%n_s - get_omega_s(N, x, y, z) * (N(x, y, z)%n_s - sN_eq)
            end if
#endif

          ! End of code which is executed on fluid sites.
          else

            ! ==========================
            ! BOUNCE-BACK FOR ROCK SITES
            ! ==========================

            ! Mid-grid bounce-back no-slip boundary condition: if in rock, turn around.
            if(bcsel == 2) then
              ! tmp_Nn_r = 0.d0
              ! tmp_Nn_b = 0.d0
              ! tmp_Nn_s = 0.d0
              tmp_Nn_r = N(x, y, z)%n_r
#ifndef SINGLEFLUID
              tmp_Nn_b = N(x, y, z)%n_b
#endif
#ifndef NOSURFACTANT
              tmp_Nn_s = N(x, y, z)%n_s
#endif
              bb_bc_nonrest: do s = 1, nnonrest
                N(x, y, z)%n_r(bounce(s)) = tmp_Nn_r(s)
#ifndef SINGLEFLUID
                N(x, y, z)%n_b(bounce(s)) = tmp_Nn_b(s)
#endif
#ifndef NOSURFACTANT
                N(x, y, z)%n_s(bounce(s)) = tmp_Nn_s(s)

                ! Handle dipoles. Copied this from the on-grid bounce-back no slip (bcsel==0).
                ! Since this is a shift from pre-advection to post-collision, it should be equivalent.
                ! I (schmie) am not 100 percent sure, though.
                xa = x + cx(s)
                ya = y + cy(s)
                za = z + cz(s)
                d_adv(:, bounce(s), x, y, z) = N(xa, ya, za)%d
#endif
              end do bb_bc_nonrest
            endif
          endif
        end do
      end do
    end do
    ! Finished running over the lattice.

    call stop_timer(ti_intc)

#ifdef AXISYM
    deallocate(psi_mat)

    deallocate(dy_psi_r)
    deallocate(dyz_psi_r)
    deallocate(dyy_psi_r)
#ifndef SINGLEFLUID
    deallocate(dy_psi_b)
    deallocate(dyz_psi_b)
    deallocate(dyy_psi_b)
#endif
#ifndef NOSURFACTANT
    deallocate(dy_psi_s)
    deallocate(dyz_psi_s)
    deallocate(dyy_psi_s)
#endif
#endif
end subroutine lbe_interactions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> subroutine in which the Boundary condition is implemented
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! if started from lbe_collission_simple.F90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef SINGLEFLUID
#ifndef LOCALBC
subroutine lbe_boundary_condition(N,d_adv,x,y,z,Nbuffer)
#endif
#ifdef LOCALBC
subroutine lbe_boundary_condition(N,d_adv,x,y,z,Nbuffer,useddirections)
#endif
    implicit none
!#ifndef INTERPOLATEDBB    
!type(lbe_site), dimension(0:,0:,-1:),intent(inout) :: Nbuffer
!#else
type(lbe_site), dimension(1-halo_extent:,1-halo_extent:,-halo_extent:),intent(inout) :: Nbuffer
!#endif
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! if started from lbe_collission.F90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifndef SINGLEFLUID
#ifndef LOCALBC
  subroutine lbe_boundary_condition(N,d_adv,x,y,z)
#endif
#ifdef LOCALBC
subroutine lbe_boundary_condition(N,d_adv,x,y,z,useddirections)
#endif
implicit none
#endif

    integer, intent(in) :: x,y,z
    integer :: xa,ya,za,xb,yb,zb,bcz,xa1,xa2,ya1,ya2,za1,za2
    integer :: s, iter, g, Istat, iterate, counter,e,i,index1,index2
!#ifndef INTPEROLATEDBB    
!    type(lbe_site), dimension(0:,0:,0:), intent(inout) :: N
!#else
    type(lbe_site), dimension(1-halo_extent:,1-halo_extent:,1-halo_extent:), intent(inout) :: N
!#endif

#ifdef SINGLEFLUID
    integer :: d_adv
#ifdef INTERPOLATEDBB
    real(kind=rk) :: dist,q,q1,q2
    integer :: iPop,iPop2,xb1,yb1,zb1
#endif
#else
#ifndef NOSURFACTANT
    real*8, dimension(1:,1:,0:,0:,0:), intent(inout) :: d_adv
#else
    integer, intent(inout) :: d_adv
#endif
#endif
#ifdef LOCALBC
        type(bc_check),dimension(1:,1:,1:), intent(inout) :: useddirections
        integer :: slip_dir
#endif
    real*8 :: init,finl,diff,acccoef2
    integer :: xx,yy,zz

    integer :: xp,yp,zp,xq,yq,zq,sp,sq
    integer :: fooflag=0
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer, dimension(6,4) :: neig        !diagonal neighbour cells
    integer, dimension(6,4,2) :: refl
    integer, dimension(18) :: oppos
    integer, dimension(12,4) :: edge
    integer, dimension(1:3) :: dir1, dir2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#ifdef NOEDGESTEP
    integer, dimension(19,2) :: projection
    integer xaa,yaa,zaa,xab,yab,zab,xba,yba,zba,xbb,ybb,zbb

    !projections of diagonal velocities:
    projection(1,:)  = (/ 19,19 /)
    projection(2,:)  = (/ 19,19 /)
    projection(3,:)  = (/ 19,19 /)
    projection(4,:)  = (/ 19,19 /)
    projection(5,:)  = (/ 19,19 /)
    projection(6,:)  = (/ 19,19 /)
    projection(7,:)  = (/ 1, 3 /)
    projection(8,:)  = (/ 1, 4 /)
    projection(9,:)  = (/ 1, 5 /)
    projection(10,:) = (/ 1, 6 /)
    projection(11,:) = (/ 2, 3 /)
    projection(12,:) = (/ 2, 4 /)
    projection(13,:) = (/ 2, 5 /)
    projection(14,:) = (/ 2, 6 /)
    projection(15,:) = (/ 3, 5 /)
    projection(16,:) = (/ 3, 6 /)
    projection(17,:) = (/ 4, 5 /)
    projection(18,:) = (/ 4, 6 /)
    projection(19,:) = (/ 19,19 /)

#endif


        oppos(1) = 2   ! index of the vector in the opposite direction
        oppos(2) = 1
        oppos(3) = 4
        oppos(4) = 3
        oppos(5) = 6
        oppos(6) = 5
        oppos(7) = 12
        oppos(8) = 11
        oppos(9) = 14
        oppos(10) = 13
        oppos(11) = 8
        oppos(12) = 7
        oppos(13) = 10
        oppos(14) = 9
        oppos(15) = 18
        oppos(16) = 17
        oppos(17) = 16
        oppos(18) = 15

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        neig(1,:) = (/ 7,8,9,10 /) !for the first 6 velocities (along axis)
        neig(2,:) = (/ 11,12,13,14 /) !the 4 "diagonal" vectors that point in
        neig(3,:) = (/ 7,11,15,16 /)  !the same direction
        neig(4,:) = (/ 8,12,17,18 /)
        neig(5,:) = (/ 15,17,13,9 /)
        neig(6,:) = (/ 14,10,16,18 /)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    refl(1,1,1) = 3                ! first index the 6 "axis-velocities"
    refl(1,1,2) = 8   ! 2. index the 4 diagonal velocities from neig()
    refl(1,2,1) = 4 !1 -->  index of axis-vector associated with diagonal
    refl(1,2,2) = 7 !velocities 2 --> diagonal velocities
    refl(1,3,1) = 5
    refl(1,3,2) = 10
    refl(1,4,1) = 6
    refl(1,4,2) = 9

    refl(2,1,1) = 3
    refl(2,1,2) = 12
    refl(2,2,1) = 4
    refl(2,2,2) = 11
    refl(2,3,1) = 5
    refl(2,3,2) = 14
    refl(2,4,1) = 6
    refl(2,4,2) = 13

    refl(3,1,1) = 1
    refl(3,1,2) = 11
    refl(3,2,1) = 2
    refl(3,2,2) = 7
    refl(3,3,1) = 5
    refl(3,3,2) = 16
    refl(3,4,1) = 6
    refl(3,4,2) = 15

    refl(4,1,1) = 1
    refl(4,1,2) = 12
    refl(4,2,1) = 2
    refl(4,2,2) = 8
    refl(4,3,1) = 5
    refl(4,3,2) = 18
    refl(4,4,1) = 6
    refl(4,4,2) = 17

    refl(5,1,1) = 3
    refl(5,1,2) = 17
    refl(5,2,1) = 4
    refl(5,2,2) = 15
    refl(5,3,1) = 2
    refl(5,3,2) = 9
    refl(5,4,1) = 1
    refl(5,4,2) = 13

    refl(6,1,1) = 2
    refl(6,1,2) = 10
    refl(6,2,1) = 1
    refl(6,2,2) = 14
    refl(6,3,1) = 3
    refl(6,3,2) = 18
    refl(6,4,1) = 4
    refl(6,4,2) = 16
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        edge(1,:) = (/ 1,3,7,12 /)
        edge(2,:) = (/ 1,4,8,11 /)
        edge(3,:) = (/ 4,2,12,7 /)
        edge(4,:) = (/ 2,3,11,8 /)
        edge(5,:) = (/ 3,5,15,18 /)
        edge(6,:) = (/ 5,4,17,16 /)
        edge(7,:) = (/ 4,6,18,15 /)
        edge(8,:) = (/ 6,3,16,17 /)
        edge(9,:) = (/ 5,1,9,14 /)
        edge(10,:) = (/ 1,6,10,13 /)
        edge(11,:) = (/ 6,2,14,9 /)
        edge(12,:) = (/ 2,5,13,10 /)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     if (bcsel == 2) then       ! Choose BC ; 0 --> BB
! Mid-grid bounce-back is handled in lbe_interactions
! Doing nothing here.
#ifdef INTERPOLATEDBB
!        write(msgstr&
!             &,"('Call before actual interpolation..')")
!        write(msgstr&
!             &,"('(',I0,',',I0,',',I0,')')")&
!             & x,y,z
!        call log_msg(msgstr)
        do s=1,nnonrest
           xa = x + cx(s)
           ya = y + cy(s)
           za = z + cz(s)

           if (N(xa,ya,za)%rock_state .ne. 0) then
              if (Nbuffer(x,y,0)%rock_state ==0) then
                 dist = N(x,y,z)%delta_q(s)
                 if(dist < 0 .or. dist > 1) then
                   write(msgstr&
                        &,"('Invalid dist at (',I0,',',I0,',',I0,'), s = ',I0, F16.10)")&
                        & x,y,z,s,dist
                   call log_msg(msgstr)
                   write(msgstr&
                        &,"('Rockstate at (',I0,',',I0,',',I0,') = ' F16.10)")&
                        & xa,ya,za,N(xa,ya,za)%rock_state
                   call log_msg(msgstr)
                   write(msgstr&
                        &,"('Rockstate at (',I0,',',I0,',',I0,') = ' F16.10)")&
                        & x,y,z,Nbuffer(x,y,0)%rock_state
                   call log_msg(msgstr)
                   call Abend
                 end if

                 ! bouzidi interpolation scheme
                 ! Original paper bouzidi et al. (2001). For implementation
                 ! refer to Lallemand and Luo (2003)
!                 if(interpolation_scheme.eq.1) then
                   if(interpolation_order.eq.1) then ! linear
                     if(dist<0.5) then
                       xb = x
                       yb = y
                       zb = 0
                       q = 2.0_rk*dist
                       iPop2 = s
                     else
                       xb = x-cx(s)
                       yb = y-cy(s)
                       zb = -cz(s)
                       q = 1.0_rk/(2.0_rk*dist)
                       iPop2 = bounce(s)
                     endif
                     ! if the node in wall's opposite direction of fluid is also solid,
                     ! the code would crash.
                     ! a hack to recover midgrid bounceback for such case
                     if(N(x-cx(s),y-cy(s),z-cz(s))%rock_state .ne. 0) then
                       xb = x       
                       yb = y
                       zb = 0
                       q = 1._rk
                       iPop2 = s
                     endif
                     N(x,y,z)%n_r(bounce(s)) = q*Nbuffer(xa,ya,cz(s))%n_r(s)&
                                         &+(1.0_rk-q)*Nbuffer(xb,yb,zb)%n_r(iPop2)
                  ! write(msgstr&
                  !      &,"('After interpolation')")
                  ! call log_msg(msgstr)
!                   ! Mixed quadratic and linear interpolation (only one layer halo cells required)                                         
!                     if(dist<0.5) then
!                       N(x,y,z)%n_r(bounce(s)) = &
!                             &dist*(1.0_rk+2.0_rk*dist)*Nbuffer(xa,ya,cz(s))%n_r(s)&
!                             &+ (1.0_rk-4.0_rk*dist**2) * Nbuffer(x,y,0)%n_r(s)&
!                             &- dist*(1.0_rk-2.0_rk*dist)*Nbuffer(x-cx(s),y-cy(s),-cz(s))%n_r(s)
!                     else
!                       xb = x-cx(s)
!                       yb = y-cy(s)
!                       zb = -cz(s)
!                       q = 1.0_rk/(2.0_rk*dist)
!                       iPop2 = bounce(s)
!                       N(x,y,z)%n_r(bounce(s)) = q*Nbuffer(xa,ya,cz(s))%n_r(s)&
!                                              &+(1.0_rk-q)*Nbuffer(xb,yb,zb)%n_r(iPop2)
!                     endif
                   else if(interpolation_order.eq.2) then ! quadratic
                     if(dist<0.5) then
                       xb = x
                       yb = y
                       zb = 0
                       xb1 = x-cx(s)
                       yb1 = y-cy(s)
                       zb1 = -cz(s)
                       q  = dist
                       q1 = 1.0_rk+2.0_rk*dist
                       q2 = 1.0_rk-2.0_rk*dist
                       iPop2 = s
                     else
                       xb = x-cx(s)
                       yb = y-cy(s)
                       zb = -cz(s)
                       xb1 = x-2*cx(s)
                       yb1 = y-2*cy(s)
                       zb1 = -2*cz(s)
                       q  = 1.0_rk/(1.0_rk+2.0_rk*dist)
                       q1 = 1.0_rk/dist
                       q2 = 2.0_rk*dist-1.0_rk
                       iPop2 = bounce(s)
                     endif
                     N(x,y,z)%n_r(bounce(s)) = q*q1 *Nbuffer(xa,ya,cz(s))%n_r(s) + & 
                                             &q1*q2*Nbuffer(xb,yb,zb)%n_r(iPop2) - &
                                             &q2*q *Nbuffer(xb1,yb1,zb1)%n_r(iPop2)
                   endif ! linear or quadratic interpolation
!
!                 ! Yu interpolation scheme. Refer to Yu et al. (2003) titled
!                 ! "A Unified Boundary Treatment in Lattice Boltzmann Method" or
!                 ! "viscous flow computations with the method of lattice
!                 ! boltzmann equation"
!                 else if(interpolation_scheme.eq.2) then
!                   if(interpolation_order.eq.1) then ! linear
!                     xb = x-cx(s)
!                     yb = y-cy(s)
!                     zb = -cz(s)
!                     N(x,y,z)%n_r(bounce(s)) = 1/(1+dist)*((1-dist)*Nbuffer(x,y,0)%n_r(s)&
!                        &+dist*Nbuffer(xa,ya,cz(s))%n_r(s)+dist*Nbuffer(xb,yb,zb)%n_r(bounce(s)))
!                   else if(interpolation_order.eq.2) then ! quadratic
!                     ! To be implemented
!                   end if
!                 end if ! Linear or quadratic interpolation

 !                  if(dist<0.5) then
!!                   write(msgstr&
!!                        &,"('dist less than 0.5 (',I0,',',I0,',',I0,'), s = ',I0, F16.10)")&
!!                        & x,y,z,s,dist
!!                   call log_msg(msgstr)
 !                    xb = x
 !                    yb = y
 !                    zb = 0
 !                    q = 2*dist
 !                    iPop2 = s
 !                  else
!!                   write(msgstr&
!!                        &,"('dist more than 0.5 (',I0,',',I0,',',I0,'), s = ',I0, F16.10)")&
!!                        & x,y,z,s,dist
!!                   call log_msg(msgstr)
 !                    xb = x-cx(s)
 !                    yb = y-cy(s)
 !                    zb = -cz(s)
 !                    q = 1/(2*dist)
 !                    iPop2 = bounce(s)
 !                  endif

!!           write(msgstr&
!!                &,"('Call before actual interpolation..')")
!!           call log_msg(msgstr)
!!            if(xa.eq.6.and.ya.eq.11.and.za.eq.6) then
!!                   write(msgstr&
!!                        &,"('Nbuffer xb yb zb at (',I0,',',I0,',',I0,'), s = ',I0, F16.10)")&
!!                        & xb,yb,zb,s,Nbuffer(xb,yb,zb)%n_r(iPop2)
!!                   call log_msg(msgstr)
!!                   write(msgstr&
!!                        &,"('Nbuffer xa ya za at (',I0,',',I0,',',I0,'), s = ',I0, F16.10)")&
!!                        & xa,ya,za,s,Nbuffer(xa,ya,za)%n_r(s)
!!                   call log_msg(msgstr)
 !                  N(x,y,z)%n_r(bounce(s))=q*Nbuffer(xa,ya,cz(s))%n_r(s)+(1-q)*Nbuffer(xb,yb,zb)%n_r(iPop2)
 !                  !N(x,y,z)%n_r(bounce(s))=Nbuffer(xa,ya,cz(s))%n_r(s)
!!             endif
!!!                   N(x,y,z)%n_r(bounce(s))=Nbuffer(x,y,0)%n_r(s)
!!
!!           write(msgstr&
!!                &,"('Call after actual interpolation..')")
!!           call log_msg(msgstr)

            endif
          endif
        end do
#endif
!#ifdef INTERPOLATEDBB
!!        write(msgstr&
!!             &,"('Call before actual interpolation..')")
!!        call log_msg(msgstr)
!        do s=1,nnonrest
!           xa = x + cx(s)
!           ya = y + cy(s)
!           za = z + cz(s)
!
!           if (N(xa,ya,za)%rock_state .ne. 0) then
!
!              if (Nbuffer(x,y,0)%rock_state ==0) then
!                 dist = N(x,y,z)%delta_q(s)
!                 if(dist < 0 .or. dist > 1) then
!                   write(msgstr&
!                        &,"('Invalid dist at (',I0,',',I0,',',I0,'), s = ',I0, F16.10)")&
!                        & x,y,z,s,dist
!                   call log_msg(msgstr)
!                   write(msgstr&
!                        &,"('Rockstate at (',I0,',',I0,',',I0,') = ' F16.10)")&
!                        & xa,ya,za,N(xa,ya,za)%rock_state
!                   call log_msg(msgstr)
!                   write(msgstr&
!                        &,"('Rockstate at (',I0,',',I0,',',I0,') = ' F16.10)")&
!                        & x,y,z,Nbuffer(x,y,0)%rock_state
!                   call log_msg(msgstr)
!                   call Abend
!                 end if
!                 
!                 if(dist<0.5) then
!                   xb = x-cx(s)
!                   yb = y-cy(s)
!                   zb = -cz(s)
!                   q = 2*dist
!                   iPop2 = s
!                 else
!                   xb = x
!                   yb = y
!                   zb = 0
!                   q = 1/(2*dist)
!                   iPop2 = bounce(s)
!                 endif
!
!                 N(x,y,z)%n_r(bounce(s))=q*Nbuffer(x,y,0)%n_r(s)+(1-q)*Nbuffer(xb,yb,zb)%n_r(iPop2)
!              endif
!          endif
!        end do
!#endif
     elseif (bcsel == 0) then       ! Choose BC ; 0 --> BB
        do s=1,nnonrest
           xa = x + cx(s)
           ya = y + cy(s)
           za = z + cz(s)

#ifdef NOEDGESTEP
           xaa = x + cx(projection(s,1))
           yaa = y + cy(projection(s,1))
           zaa = z + cz(projection(s,1))
           xab = x + cx(projection(s,2))
           yab = y + cy(projection(s,2))
           zab = z + cz(projection(s,2))
#endif

           if (N(xa,ya,za)%rock_state .ne. 0) then

#ifdef SINGLEFLUID
              if (Nbuffer(x,y,0)%rock_state ==0) then
                 N(x,y,z)%n_r(bounce(s))=Nbuffer(x,y,0)%n_r(s)
#else
              if (Nbuf(x,y,0)%rock_state ==0) then
                 N(x,y,z)%n_r(bounce(s))=Nbuf(x,y,0)%n_r(s)
                 N(x,y,z)%n_b(bounce(s))=Nbuf(x,y,0)%n_b(s)
#endif
#ifndef NOSURFACTANT
                 N(x,y,z)%n_s(bounce(s))=Nbuf(x,y,0)%n_s(s)
                 d_adv(:,bounce(s),x,y,z)=Nbuf(xa,ya,cz(s))%d
#endif

#ifdef NOEDGESTEP
              else if (N(xaa,yaa,zaa)%rock_state.ne.0&
                   &.and.N(xab,yab,zab)%rock_state.ne.0) then
#ifdef SINGLEFLUID
                 if(Nbuffer(x,y,0)%rock_state ==0 .and. s .gt. 6 ) then
                    N(x,y,z)%n_r(bounce(s))=Nbuffer(x,y,0)%n_r(s)
#else
                 if(Nbuf(x,y,0)%rock_state ==0 .and. s .gt. 6 ) then
                    N(x,y,z)%n_r(bounce(s))=Nbuf(x,y,0)%n_r(s)
                    N(x,y,z)%n_b(bounce(s))=Nbuf(x,y,0)%n_b(s)
#endif
#ifndef NOSURFACTANT
                    N(x,y,z)%n_s(bounce(s))=Nbuf(x,y,0)%n_s(s)
                    d_adv(:,bounce(s),x,y,z)=Nbuf(xa,ya,cz(s))%d
#endif
                 endif
#endif
              endif
          endif
        end do
     else if (bcsel == 1) then ! BC selection
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!bcsel=1 --> slip bc (might be wrong for the dipole moment of the surfactant)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifndef LOCALBC

        if(N(x,y,z)%rock_state==0) then        !start check for rock in xyz

        !!!!!!!!!!!!!!!
        !start of g loop
        !!!!!!!!!!!!!!!
        do g=1,6
                                xa = x + cx(g)
                                ya = y + cy(g)
                                za = z + cz(g)


        !!!!!!!!!!!!!!!!!
        !check for orthogonal rocks and reflect n_r(g)=1-6
        !!!!!!!!!!!!!!!!

        if (N(xa,ya,za)%rock_state .ne. 0) then


#ifdef SINGLEFLUID
          N(x,y,z)%n_r(bounce(g))=Nbuffer(x,y,0)%n_r(g)
#endif
#ifndef SINGLEFLUID
          N(x,y,z)%n_r(bounce(g))=Nbuf(x,y,0)%n_r(g)
            N(x,y,z)%n_b(bounce(g))=Nbuf(x,y,0)%n_b(g)
#endif
#ifndef NOSURFACTANT
          N(x,y,z)%n_s(bounce(g))=Nbuf(x,y,0)%n_s(g)
                          d_adv(:,bounce(g),x,y,z)=Nbuf(xa,ya,cz(g))%d
#endif

        !!!!!!!!!!!!!!!!
        !start iter loop
        !!!!!!!!!!!!!!!!
        do iter=1,4
        !!!!!!!!!!!!!!!!
        !check for diagonally posed rock fields
        !!!!!!!!!!!!!!!!

        if(N(x+cx(neig(g,iter)),y+cy(neig(g,iter)),z+cz(neig(g,iter)))%rock_state .ne. 0) then


        !!!!!!!!!!!!!!!!
        !check for orthogonal rock fields
        !!!!!!!!!!!!!!!!
        if( N(x+cx(refl(g,iter,1)),y+cy(refl(g,iter,1)),z+cz(refl(g,iter,1)))%rock_state==0 ) then



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if((acccoef .ge. 0 ) .and. (acccoef .le. 1))then


        !!!!!!!!!!!!!!!!
        !Add diffusive part
        !!!!!!!!!!!!!!!!

#ifdef SINGLEFLUID
    N(x,y,z)%n_r(bounce(g)) = N(x,y,z)%n_r(bounce(g)) + 0.5*acccoef&
         & * Nbuffer(x+cx(refl(g,iter,1)),y+cy(refl(g,iter,1)),cz(refl(g,iter,1)))%n_r(refl(g,iter,2))

    N(x,y,z)%n_r(oppos(neig(g,iter))) = (1-acccoef)&
         & * Nbuffer(x+cx(refl(g,iter,1)),y+cy(refl(g,iter,1)),cz(refl(g,iter,1)))%n_r(refl(g,iter,2))
#endif

#ifndef SINGLEFLUID

    N(x,y,z)%n_r(bounce(g)) = N(x,y,z)%n_r(bounce(g)) + 0.5*acccoef&
         & * Nbuf(x+cx(refl(g,iter,1)),y+cy(refl(g,iter,1)),cz(refl(g,iter,1)))%n_r(refl(g,iter,2))

    N(x,y,z)%n_r(oppos(neig(g,iter))) = (1-acccoef)&
         & * Nbuf(x+cx(refl(g,iter,1)),y+cy(refl(g,iter,1)),cz(refl(g,iter,1)))%n_r(refl(g,iter,2))

    N(x,y,z)%n_b(bounce(g)) = N(x,y,z)%n_b(bounce(g)) + 0.5*acccoef&
         & * Nbuf(x+cx(refl(g,iter,1)),y+cy(refl(g,iter,1)),cz(refl(g,iter,1)))%n_b(refl(g,iter,2))
    N(x,y,z)%n_b(oppos(neig(g,iter))) = (1-acccoef)&
         & * Nbuf(x+cx(refl(g,iter,1)),y+cy(refl(g,iter,1)),cz(refl(g,iter,1)))%n_b(refl(g,iter,2))
#endif
#ifndef NOSURFACTANT
    N(x,y,z)%n_s(bounce(g)) = N(x,y,z)%n_s(bounce(g)) + 0.5*acccoef&
         & * Nbuf(x+cx(refl(g,iter,1)),y+cy(refl(g,iter,1)),cz(refl(g,iter,1)))%n_s(refl(g,iter,2))

    d_adv(:,bounce(g),x,y,z)=d_adv(:,bounce(g),x,y,z)+0.5*acccoef&
         & *Nbuf(x+cx(refl(g,iter,1)),y+cy(refl(g,iter,1)),cz(refl(g,iter,1)))%d

    N(x,y,z)%n_s(oppos(neig(g,iter))) = (1-acccoef)&
         & * Nbuf(x+cx(refl(g,iter,1)),y+cy(refl(g,iter,1)),cz(refl(g,iter,1)))%n_s(refl(g,iter,2))

    d_adv(:,bounce(g),x,y,z)=d_adv(:,bounce(g),x,y,z)+(1-acccoef)&
         & * Nbuf(x+cx(refl(g,iter,1)),y+cy(refl(g,iter,1)),cz(refl(g,iter,1)))%d

#endif

        elseif((acccoef > 1 ).and.(acccoef .le. 2))then

        acccoef2 = 2 - acccoef

#ifdef SINGLEFLUID
 N(x,y,z)%n_r(bounce(neig(g,iter))) = N(x,y,z)%n_r(bounce(neig(g,iter))) + 0.5*acccoef2 * Nbuffer(x,y,0)%n_r(neig(g,iter))


 N(x,y,z)%n_r(bounce(neig(g,iter))) &
      &= (1-acccoef2) * Nbuffer(x,y,0)%n_r(neig(g,iter))
#endif

#ifndef SINGLEFLUID

  N(x,y,z)%n_r(bounce(neig(g,iter))) = N(x,y,z)%n_r(bounce(neig(g,iter))) + 0.5*acccoef2 * Nbuf(x,y,0)%n_r(neig(g,iter))


  N(x,y,z)%n_r(bounce(neig(g,iter))) &
       &= (1-acccoef2) * Nbuf(x,y,0)%n_r(neig(g,iter))

  N(x,y,z)%n_b(bounce(neig(g,iter))) = N(x,y,z)%n_b(bounce(neig(g,iter))) + 0.5*acccoef2 * Nbuf(x,y,0)%n_b(neig(g,iter))

  N(x,y,z)%n_b(bounce(neig(g,iter))) &
       &= (1-acccoef2) * Nbuf(x,y,0)%n_b(neig(g,iter))
#endif
#ifndef NOSURFACTANT
  N(x,y,z)%n_s(bounce(neig(g,iter))) = N(x,y,z)%n_s(bounce(neig(g,iter))) + 0.5*acccoef2 * Nbuf(x,y,0)%n_s(neig(g,iter))

  N(x,y,z)%n_s(bounce(neig(g,iter))) &
       &= (1-acccoef2) * Nbuf(x,y,0)%n_s(neig(g,iter))

  d_adv(:,bounce(s),x,y,z)=Nbuf(xa,ya,cz(s))%d

  d_adv(:,bounce(neig(g,iter)),x,y,z)&
       &=d_adv(:,bounce(neig(g,iter)),x,y,z)+0.5*acccoef2&
       & *Nbuf(x+cx(neig(g,iter)),y+cy(neig(g,iter)),cz(neig(g,iter)))%d

  d_adv(:,bounce(neig(g,iter)),x,y,z)&
       &=d_adv(:,bounce(neig(g,iter)),x,y,z)+(1-acccoef2)&
       & * Nbuf(x+cx(neig(g,iter)),y+cy(neig(g,iter)),cz(neig(g,iter)))%d
#endif

else
   print*,"acccoef has an invalid value"
end if
endif !end check for orthogonal rock fields
endif !check diag r fields
        !!!!!!!!!!!!!!!!
        !end iter loop
        !!!!!!!!!!!!!!!!
end do
endif        !end of check fields in s=1-6 directions

!!!!!!!!!!!!!!!!
!end of g loop
!!!!!!!!!!!!!!!!
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do e=1,12
        xa = x + cx(edge(e,1))
        ya = y + cy(edge(e,1))
        za = z + cz(edge(e,1))
        xa1 = x + cx(edge(e,2))
        ya1 = y + cy(edge(e,2))
        za1 = z + cz(edge(e,2))
        xa2 = x + cx(edge(e,3))
        ya2 = y + cy(edge(e,3))
        za2 = z + cz(edge(e,3))
        if((N(xa,ya,za)%rock_state == 0 .and. N(xa1,ya1,za1)%rock_state == 0&
      & .and. N(xa2,ya2,za2)%rock_state .ne. 0)&
      &.or.((N(xa,ya,za)%rock_state .ne. 0 .and. N(xa1,ya1,za1)%rock_state .ne. 0 .and. N(xa2,ya2,za2)%rock_state .ne. 0)) )&
      & then !if there is just an edge or a corner, reflect diagonal

#ifdef SINGLEFLUID
        N(x,y,z)%n_r(bounce(edge(e,3)))=Nbuffer(x,y,0)%n_r(edge(e,3))
#endif
#ifndef SINGLEFLUID
        N(x,y,z)%n_r(bounce(edge(e,3)))=Nbuf(x,y,0)%n_r(edge(e,3))
        N(x,y,z)%n_b(bounce(edge(e,3)))=Nbuf(x,y,0)%n_b(edge(e,3))
#ifndef NOSURFACTANT
        N(x,y,z)%n_s(bounce(edge(e,3)))=Nbuf(x,y,0)%n_s(edge(e,3))
        d_adv(:,bounce(edge(e,3)),x,y,z)= Nbuf(x,y,0)%d
#endif
#endif
        endif

        end do !end of e loop

        !end of check for rock in xyz
        endif


#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! In the following part the local BC is implemented
! The general approch is:
!reflect the main Axis (1<=g<7)
!for the diagonal directions there are 4 cases a corner,
! 2 types of straight walls and an edge
!the g determines the direction in which it is checked for rocks
!the incoming directions at the site x,y,z are dependent on this g
!and on the local acccoef which ranges from 0 to 1 for slip to diffuse reflection
!and from 1 to 2 for diffuse reflection to BB-BC
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#ifdef LOCALBC


        !start check for rock in xyz
        rockxyzcheck: if(N(x,y,z)%rock_state==0) then
           directionloop: do g=1,nnonrest
              xa = x + cx(g)
              ya = y + cy(g)
              za = z + cz(g)
              rockdirectioncheck: if (N(xa,ya,za)%rock_state .ne. 0) then
                 reflectaxes: if(g.le.6)then
                    usedcheck1: if(useddirections(x,y,z)%n(bounce(g))==0)then
#ifdef SINGLEFLUID
                       N(x,y,z)%n_r(bounce(g))=Nbuffer(x,y,0)%n_r(g)
#endif
#ifndef SINGLEFLUID
                       N(x,y,z)%n_r(bounce(g))=Nbuf(x,y,0)%n_r(g)
        N(x,y,z)%n_b(bounce(g))=Nbuf(x,y,0)%n_b(g)
#endif
#ifndef NOSURFACTANT
        N(x,y,z)%n_s(bounce(g))=Nbuf(x,y,0)%n_s(g)
                          d_adv(:,bounce(g),x,y,z)=Nbuf(xa,ya,cz(g))%d
#endif

useddirections(x,y,z)%n(bounce(g)) = useddirections(x,y,z)%n(bounce(g))+1
else usedcheck1

#ifdef SINGLEFLUID
        N(x,y,z)%n_r(bounce(g))=Nbuffer(x,y,0)%n_r(g)
#endif
#ifndef SINGLEFLUID
        N(x,y,z)%n_r(bounce(g))=Nbuf(x,y,0)%n_r(g)
        N(x,y,z)%n_b(bounce(g))=Nbuf(x,y,0)%n_b(g)
#endif
#ifndef NOSURFACTANT
        N(x,y,z)%n_s(bounce(g))=Nbuf(x,y,0)%n_s(g)
                          d_adv(:,bounce(g),x,y,z)=Nbuf(xa,ya,cz(g))%d
#endif

useddirections(x,y,z)%n(bounce(g)) = useddirections(x,y,z)%n(bounce(g))+1
end if usedcheck1


else reflectaxes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Directions not matching the main Axis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        call lbe_get_direction(g,dir1,dir2,index1,index2)

getgeometry:if(N(x+dir1(1),y+dir1(2),z+dir1(3))%rock_state.ne.0 .and. N(x+dir2(1),y+dir2(2),z+dir2(3))%rock_state.ne.0)then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Corner
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

usedcheck2: if(useddirections(x,y,z)%n(bounce(g))==0)then
#ifdef SINGLEFLUID
        N(x,y,z)%n_r(bounce(g))=Nbuffer(x,y,0)%n_r(g)
#endif
#ifndef SINGLEFLUID
        N(x,y,z)%n_r(bounce(g))=Nbuf(x,y,0)%n_r(g)
        N(x,y,z)%n_b(bounce(g))=Nbuf(x,y,0)%n_b(g)
#endif
#ifndef NOSURFACTANT
        N(x,y,z)%n_s(bounce(g))=Nbuf(x,y,0)%n_s(g)
                          d_adv(:,bounce(g),x,y,z)=Nbuf(xa,ya,cz(g))%d
#endif

useddirections(x,y,z)%n(bounce(g)) = useddirections(x,y,z)%n(bounce(g))+1
else usedcheck2

#ifdef SINGLEFLUID
        N(x,y,z)%n_r(bounce(g))=Nbuffer(x,y,0)%n_r(g)
#endif
#ifndef SINGLEFLUID
        N(x,y,z)%n_r(bounce(g))=Nbuf(x,y,0)%n_r(g)
        N(x,y,z)%n_b(bounce(g))=Nbuf(x,y,0)%n_b(g)
#endif
#ifndef NOSURFACTANT
        N(x,y,z)%n_s(bounce(g))=Nbuf(x,y,0)%n_s(g)
                          d_adv(:,bounce(g),x,y,z)=Nbuf(xa,ya,cz(g))%d
#endif
useddirections(x,y,z)%n(bounce(g)) = useddirections(x,y,z)%n(bounce(g))+1
end if usedcheck2



else if(N(x+dir1(1),y+dir1(2),z+dir1(3))%rock_state==0 .and. N(x+dir2(1),y+dir2(2),z+dir2(3))%rock_state.ne.0)then getgeometry
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!straight wall 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        acccoef=N(xa,ya,za)%local_acccoef
acccoefcheck1:        if((acccoef .ge. 0 ) .and. (acccoef .le. 1))then


usedcheck3: if(useddirections(x,y,z)%n(bounce(g))==0)then
#ifdef SINGLEFLUID
        N(x,y,z)%n_r(bounce(g))=0
#endif
#ifndef SINGLEFLUID
        N(x,y,z)%n_r(bounce(g))=0
        N(x,y,z)%n_b(bounce(g))=0
#endif
#ifndef NOSURFACTANT
        N(x,y,z)%n_s(bounce(g))=0
        d_adv(:,bounce(g),x,y,z)=(/0,0,0/)
#endif
useddirections(x,y,z)%n(bounce(g)) = useddirections(x,y,z)%n(bounce(g))+1
else usedcheck3
useddirections(x,y,z)%n(bounce(g)) = useddirections(x,y,z)%n(bounce(g))+1
end if usedcheck3

usedcheck4: if(useddirections(x,y,z)%n(bounce(index2))==0)then
#ifdef SINGLEFLUID
        N(x,y,z)%n_r(bounce(index2))=0
#endif
#ifndef SINGLEFLUID
        N(x,y,z)%n_r(bounce(index2))=0
        N(x,y,z)%n_b(bounce(index2))=0
#endif
#ifndef NOSURFACTANT
        N(x,y,z)%n_s(bounce(index2))=0
        d_adv(:,bounce(index2),x,y,z)=(/0,0,0/)
#endif
useddirections(x,y,z)%n(bounce(index2)) = useddirections(x,y,z)%n(bounce(index2))+1
else usedcheck4
useddirections(x,y,z)%n(bounce(index2)) = useddirections(x,y,z)%n(bounce(index2))+1
end if usedcheck4


elseif((acccoef .gt. 1 ).and.(acccoef .le. 2))then        acccoefcheck1
        acccoef2 = 2 - acccoef



usedcheck5: if(useddirections(x,y,z)%n(bounce(g))==0)then
#ifdef SINGLEFLUID
        N(x,y,z)%n_r(bounce(g))=(1-acccoef2)*Nbuffer(x,y,0)%n_r(g)
#endif
#ifndef SINGLEFLUID
        N(x,y,z)%n_r(bounce(g))=(1-acccoef2)*Nbuf(x,y,0)%n_r(g)
        N(x,y,z)%n_b(bounce(g))=(1-acccoef2)*Nbuf(x,y,0)%n_b(g)
#endif
#ifndef NOSURFACTANT
        N(x,y,z)%n_s(bounce(g))=(1-acccoef2)*Nbuf(x,y,0)%n_s(g)
        d_adv(:,bounce(g),x,y,z)=(1-acccoef2)*Nbuf(xa,ya,cz(g))%d
#endif
useddirections(x,y,z)%n(bounce(g)) = useddirections(x,y,z)%n(bounce(g))+1
else usedcheck5
#ifdef SINGLEFLUID
        N(x,y,z)%n_r(bounce(g))=N(x,y,z)%n_r(bounce(g))+(1-acccoef2)*Nbuffer(x,y,0)%n_r(g)
#endif
#ifndef SINGLEFLUID
        N(x,y,z)%n_r(bounce(g))=N(x,y,z)%n_r(bounce(g))+(1-acccoef2)*Nbuf(x,y,0)%n_r(g)
        N(x,y,z)%n_b(bounce(g))=N(x,y,z)%n_b(bounce(g))+(1-acccoef2)*Nbuf(x,y,0)%n_b(g)
#endif
#ifndef NOSURFACTANT
        N(x,y,z)%n_s(bounce(g))=N(x,y,z)%n_s(bounce(g))+(1-acccoef2)*Nbuf(x,y,0)%n_s(g)
        d_adv(:,bounce(g),x,y,z)=d_adv(:,bounce(g),x,y,z)+(1-acccoef2)*Nbuf(xa,ya,cz(g))%d
#endif

useddirections(x,y,z)%n(bounce(g)) = useddirections(x,y,z)%n(bounce(g))+1
end if usedcheck5



usedcheck6: if(useddirections(x,y,z)%n(bounce(index2))==0)then
#ifdef SINGLEFLUID
        N(x,y,z)%n_r(bounce(index2))=0.5*acccoef2*Nbuffer(x,y,0)%n_r(g)
#endif
#ifndef SINGLEFLUID
        N(x,y,z)%n_r(bounce(index2))=0.5*acccoef2*Nbuf(x,y,0)%n_r(g)
        N(x,y,z)%n_b(bounce(index2))=0.5*acccoef2*Nbuf(x,y,0)%n_b(g)
#endif
#ifndef NOSURFACTANT
        N(x,y,z)%n_s(bounce(index2))=0.5*acccoef2*Nbuf(x,y,0)%n_s(g)
        d_adv(:,bounce(index2),x,y,z)=0.5*acccoef2*Nbuf(xa,ya,cz(g))%d
#endif
useddirections(x,y,z)%n(bounce(index2)) = useddirections(x,y,z)%n(bounce(index2))+1
else usedcheck6

#ifdef SINGLEFLUID
        N(x,y,z)%n_r(bounce(index2))=N(x,y,z)%n_r(bounce(index2))+0.5*acccoef2*Nbuffer(x,y,0)%n_r(g)
#endif
#ifndef SINGLEFLUID
        N(x,y,z)%n_r(bounce(index2))=N(x,y,z)%n_r(bounce(index2))+0.5*acccoef2*Nbuf(x,y,0)%n_r(g)
        N(x,y,z)%n_b(bounce(index2))=N(x,y,z)%n_b(bounce(index2))+0.5*acccoef2*Nbuf(x,y,0)%n_b(g)
#endif
#ifndef NOSURFACTANT
        N(x,y,z)%n_s(bounce(index2))=N(x,y,z)%n_s(bounce(index2))+0.5*acccoef2*Nbuf(x,y,0)%n_s(g)
        d_adv(:,bounce(index2),x,y,z)=d_adv(:,bounce(index2),x,y,z)+0.5*acccoef2*Nbuf(xa,ya,cz(g))%d
#endif

useddirections(x,y,z)%n(bounce(index2)) = useddirections(x,y,z)%n(bounce(index2))+1
end if usedcheck6


end if acccoefcheck1




        acccoef=N(x+dir2(1),y+dir2(2),z+dir2(3))%local_acccoef
call lbe_get_direction2(dir2,dir1,slip_dir)

acccoefcheck2:        if((acccoef .ge. 0 ) .and. (acccoef .le. 1))then


usedcheck7: if(useddirections(x,y,z)%n(bounce(g))==0)then
#ifdef SINGLEFLUID
        N(x,y,z)%n_r(bounce(g))=(1-acccoef)*Nbuffer(x+dir1(1),y+dir1(2),dir1(3))%n_r(slip_dir)
#endif
#ifndef SINGLEFLUID
        N(x,y,z)%n_r(bounce(g))=(1-acccoef)*Nbuf(x+dir1(1),y+dir1(2),dir1(3))%n_r(slip_dir)
        N(x,y,z)%n_b(bounce(g))=(1-acccoef)*Nbuf(x+dir1(1),y+dir1(2),dir1(3))%n_b(slip_dir)
#endif
#ifndef NOSURFACTANT
        N(x,y,z)%n_s(bounce(g))=(1-acccoef)*Nbuf(x+dir1(1),y+dir1(2),dir1(3))%n_s(slip_dir)
        d_adv(:,bounce(g),x,y,z)=(1-acccoef)*Nbuf(xa,ya,cz(g))%d
#endif
useddirections(x,y,z)%n(bounce(g)) = useddirections(x,y,z)%n(bounce(g))+1
else usedcheck7
#ifdef SINGLEFLUID
        N(x,y,z)%n_r(bounce(g))=N(x,y,z)%n_r(bounce(g))+(1-acccoef)*Nbuffer(x+dir1(1),y+dir1(2),dir1(3))%n_r(slip_dir)
#endif
#ifndef SINGLEFLUID
        N(x,y,z)%n_r(bounce(g))=N(x,y,z)%n_r(bounce(g))+(1-acccoef)*Nbuf(x+dir1(1),y+dir1(2),dir1(3))%n_r(slip_dir)
        N(x,y,z)%n_b(bounce(g))=N(x,y,z)%n_b(bounce(g))+(1-acccoef)*Nbuf(x+dir1(1),y+dir1(2),dir1(3))%n_b(slip_dir)
#endif
#ifndef NOSURFACTANT
        N(x,y,z)%n_s(bounce(g))=N(x,y,z)%n_s(bounce(g))+(1-acccoef)*Nbuf(x+dir1(1),y+dir1(2),dir1(3))%n_s(slip_dir)
        d_adv(:,bounce(g),x,y,z)=d_adv(:,bounce(g),x,y,z)+(1-acccoef)*Nbuf(xa,ya,cz(g))%d
#endif
useddirections(x,y,z)%n(bounce(g)) = useddirections(x,y,z)%n(bounce(g))+1
end if usedcheck7



usedcheck8: if(useddirections(x,y,z)%n(bounce(index2))==0)then
#ifdef SINGLEFLUID
        N(x,y,z)%n_r(bounce(index2))=0.5*acccoef*Nbuffer(x+dir1(1),y+dir1(2),dir1(3))%n_r(slip_dir)
#endif
#ifndef SINGLEFLUID
        N(x,y,z)%n_r(bounce(index2))=0.5*acccoef*Nbuf(x+dir1(1),y+dir1(2),dir1(3))%n_r(slip_dir)
        N(x,y,z)%n_b(bounce(index2))=0.5*acccoef*Nbuf(x+dir1(1),y+dir1(2),dir1(3))%n_b(slip_dir)
#endif
#ifndef NOSURFACTANT
        N(x,y,z)%n_s(bounce(index2))=0.5*acccoef*Nbuf(x+dir1(1),y+dir1(2),dir1(3))%n_s(slip_dir)
        d_adv(:,bounce(index2),x,y,z)=0.5*acccoef*Nbuf(xa,ya,cz(g))%d
#endif
useddirections(x,y,z)%n(bounce(index2)) = useddirections(x,y,z)%n(bounce(index2))+1
else usedcheck8
#ifdef SINGLEFLUID
        N(x,y,z)%n_r(bounce(index2))=N(x,y,z)%n_r(bounce(index2))+0.5*acccoef*Nbuffer(x+dir1(1),y+dir1(2),dir1(3))%n_r(slip_dir)
#endif
#ifndef SINGLEFLUID
        N(x,y,z)%n_r(bounce(index2))=N(x,y,z)%n_r(bounce(index2))+0.5*acccoef*Nbuf(x+dir1(1),y+dir1(2),dir1(3))%n_r(slip_dir)
        N(x,y,z)%n_b(bounce(index2))=N(x,y,z)%n_b(bounce(index2))+0.5*acccoef*Nbuf(x+dir1(1),y+dir1(2),dir1(3))%n_b(slip_dir)
#endif
#ifndef NOSURFACTANT
        N(x,y,z)%n_s(bounce(index2))=N(x,y,z)%n_s(bounce(index2))+0.5*acccoef*Nbuf(x+dir1(1),y+dir1(2),dir1(3))%n_s(slip_dir)
        d_adv(:,bounce(index2),x,y,z)=d_adv(:,bounce(index2),x,y,z)+0.5*acccoef*Nbuf(xa,ya,cz(g))%d
#endif
useddirections(x,y,z)%n(bounce(index2)) = useddirections(x,y,z)%n(bounce(index2))+1
end if usedcheck8



elseif((acccoef .gt. 1 ).and.(acccoef .le. 2))then        acccoefcheck2
        acccoef2 = 2 - acccoef

usedcheck9: if(useddirections(x,y,z)%n(bounce(g))==0)then
#ifdef SINGLEFLUID
        N(x,y,z)%n_r(bounce(g))=0
#endif
#ifndef SINGLEFLUID
        N(x,y,z)%n_r(bounce(g))=0
        N(x,y,z)%n_b(bounce(g))=0
#endif
#ifndef NOSURFACTANT
        N(x,y,z)%n_s(bounce(g))=0
        d_adv(:,bounce(g),x,y,z)=(/0,0,0/)
#endif
useddirections(x,y,z)%n(bounce(g)) = useddirections(x,y,z)%n(bounce(g))+1
else usedcheck9
useddirections(x,y,z)%n(bounce(g)) = useddirections(x,y,z)%n(bounce(g))+1
end if usedcheck9


usedcheck10: if(useddirections(x,y,z)%n(bounce(index2))==0)then
#ifdef SINGLEFLUID
        N(x,y,z)%n_r(bounce(index2))=0
#endif
#ifndef SINGLEFLUID
        N(x,y,z)%n_r(bounce(index2))=0
        N(x,y,z)%n_b(bounce(index2))=0
#endif
#ifndef NOSURFACTANT
        N(x,y,z)%n_s(bounce(index2))=0
        d_adv(:,bounce(index2),x,y,z)=0
#endif
useddirections(x,y,z)%n(bounce(index2)) = useddirections(x,y,z)%n(bounce(index2))+1
else usedcheck10
useddirections(x,y,z)%n(bounce(index2)) = useddirections(x,y,z)%n(bounce(index2))+1
end if usedcheck10

end if acccoefcheck2



else if(N(x+dir1(1),y+dir1(2),z+dir1(3))%rock_state.ne.0 .and. N(x+dir2(1),y+dir2(2),z+dir2(3))%rock_state==0)then getgeometry
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!straight wall 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        acccoef=N(xa,ya,za)%local_acccoef
acccoefcheck3:        if((acccoef .ge. 0 ) .and. (acccoef .le. 1))then


usedcheck11: if(useddirections(x,y,z)%n(bounce(g))==0)then
#ifdef SINGLEFLUID
        N(x,y,z)%n_r(bounce(g))=0
#endif
#ifndef SINGLEFLUID
        N(x,y,z)%n_r(bounce(g))=0
        N(x,y,z)%n_b(bounce(g))=0
#endif
#ifndef NOSURFACTANT
        N(x,y,z)%n_s(bounce(g))=0
        d_adv(:,bounce(g),x,y,z)=(/0,0,0/)
#endif
useddirections(x,y,z)%n(bounce(g)) = useddirections(x,y,z)%n(bounce(g))+1
else usedcheck11
useddirections(x,y,z)%n(bounce(g)) = useddirections(x,y,z)%n(bounce(g))+1
end if usedcheck11

usedcheck12: if(useddirections(x,y,z)%n(bounce(index1))==0)then
#ifdef SINGLEFLUID
        N(x,y,z)%n_r(bounce(index1))=0
#endif
#ifndef SINGLEFLUID
        N(x,y,z)%n_r(bounce(index1))=0
        N(x,y,z)%n_b(bounce(index1))=0
#endif
#ifndef NOSURFACTANT
        N(x,y,z)%n_s(bounce(index1))=0
        d_adv(:,bounce(index1),x,y,z)=(/0,0,0/)
#endif
useddirections(x,y,z)%n(bounce(index1)) = useddirections(x,y,z)%n(bounce(index1))+1
else usedcheck12
useddirections(x,y,z)%n(bounce(index1)) = useddirections(x,y,z)%n(bounce(index1))+1
end if usedcheck12


elseif((acccoef .gt. 1 ).and.(acccoef .le. 2))then        acccoefcheck3
        acccoef2 = 2 - acccoef



usedcheck13: if(useddirections(x,y,z)%n(bounce(g))==0)then
#ifdef SINGLEFLUID
        N(x,y,z)%n_r(bounce(g))=(1-acccoef2)*Nbuffer(x,y,0)%n_r(g)
#endif
#ifndef SINGLEFLUID
        N(x,y,z)%n_r(bounce(g))=(1-acccoef2)*Nbuf(x,y,0)%n_r(g)
        N(x,y,z)%n_b(bounce(g))=(1-acccoef2)*Nbuf(x,y,0)%n_b(g)
#endif
#ifndef NOSURFACTANT
        N(x,y,z)%n_s(bounce(g))=(1-acccoef2)*Nbuf(x,y,0)%n_s(g)
        d_adv(:,bounce(g),x,y,z)=(1-acccoef2)*Nbuf(xa,ya,cz(g))%d
#endif
useddirections(x,y,z)%n(bounce(g)) = useddirections(x,y,z)%n(bounce(g))+1
else usedcheck13
#ifdef SINGLEFLUID
        N(x,y,z)%n_r(bounce(g))=N(x,y,z)%n_r(bounce(g))+(1-acccoef2)*Nbuffer(x,y,0)%n_r(g)
#endif
#ifndef SINGLEFLUID
        N(x,y,z)%n_r(bounce(g))=N(x,y,z)%n_r(bounce(g))+(1-acccoef2)*Nbuf(x,y,0)%n_r(g)
        N(x,y,z)%n_b(bounce(g))=N(x,y,z)%n_b(bounce(g))+(1-acccoef2)*Nbuf(x,y,0)%n_b(g)
#endif
#ifndef NOSURFACTANT
        N(x,y,z)%n_s(bounce(g))=N(x,y,z)%n_s(bounce(g))+(1-acccoef2)*Nbuf(x,y,0)%n_s(g)
        d_adv(:,bounce(g),x,y,z)=d_adv(:,bounce(g),x,y,z)+(1-acccoef2)*Nbuf(xa,ya,cz(g))%d
#endif

useddirections(x,y,z)%n(bounce(g)) = useddirections(x,y,z)%n(bounce(g))+1
end if usedcheck13



usedcheck14: if(useddirections(x,y,z)%n(bounce(index1))==0)then
#ifdef SINGLEFLUID
        N(x,y,z)%n_r(bounce(index1))=0.5*acccoef2*Nbuffer(x,y,0)%n_r(g)
#endif
#ifndef SINGLEFLUID
        N(x,y,z)%n_r(bounce(index1))=0.5*acccoef2*Nbuf(x,y,0)%n_r(g)
        N(x,y,z)%n_b(bounce(index1))=0.5*acccoef2*Nbuf(x,y,0)%n_b(g)
#endif
#ifndef NOSURFACTANT
        N(x,y,z)%n_s(bounce(index1))=0.5*acccoef2*Nbuf(x,y,0)%n_s(g)
        d_adv(:,bounce(index1),x,y,z)=0.5*acccoef2*Nbuf(xa,ya,cz(g))%d
#endif
useddirections(x,y,z)%n(bounce(index1)) = useddirections(x,y,z)%n(bounce(index1))+1
else usedcheck14

#ifdef SINGLEFLUID
        N(x,y,z)%n_r(bounce(index1))=N(x,y,z)%n_r(bounce(index1))+0.5*acccoef2*Nbuffer(x,y,0)%n_r(g)
#endif
#ifndef SINGLEFLUID
        N(x,y,z)%n_r(bounce(index1))=N(x,y,z)%n_r(bounce(index1))+0.5*acccoef2*Nbuf(x,y,0)%n_r(g)
        N(x,y,z)%n_b(bounce(index1))=N(x,y,z)%n_b(bounce(index1))+0.5*acccoef2*Nbuf(x,y,0)%n_b(g)
#endif
#ifndef NOSURFACTANT
        N(x,y,z)%n_s(bounce(index1))=N(x,y,z)%n_s(bounce(index1))+0.5*acccoef2*Nbuf(x,y,0)%n_s(g)
        d_adv(:,bounce(index1),x,y,z)=d_adv(:,bounce(index1),x,y,z)+0.5*acccoef2*Nbuf(xa,ya,cz(g))%d
#endif

useddirections(x,y,z)%n(bounce(index1)) = useddirections(x,y,z)%n(bounce(index1))+1
end if usedcheck14


end if acccoefcheck3




        acccoef=N(x+dir1(1),y+dir1(2),z+dir1(3))%local_acccoef
call lbe_get_direction2(dir1,dir2,slip_dir)

acccoefcheck4:        if((acccoef .ge. 0 ) .and. (acccoef .le. 1))then


usedcheck15: if(useddirections(x,y,z)%n(bounce(g))==0)then
#ifdef SINGLEFLUID
        N(x,y,z)%n_r(bounce(g))=(1-acccoef)*Nbuffer(x+dir2(1),y+dir2(2),dir2(3))%n_r(slip_dir)
#endif
#ifndef SINGLEFLUID
        N(x,y,z)%n_r(bounce(g))=(1-acccoef)*Nbuf(x+dir2(1),y+dir2(2),dir2(3))%n_r(slip_dir)
        N(x,y,z)%n_b(bounce(g))=(1-acccoef)*Nbuf(x+dir2(1),y+dir2(2),dir2(3))%n_b(slip_dir)
#endif
#ifndef NOSURFACTANT
        N(x,y,z)%n_s(bounce(g))=(1-acccoef)*Nbuf(x+dir2(1),y+dir2(2),dir2(3))%n_s(slip_dir)
        d_adv(:,bounce(g),x,y,z)=(1-acccoef)*Nbuf(xa,ya,cz(g))%d
#endif
useddirections(x,y,z)%n(bounce(g)) = useddirections(x,y,z)%n(bounce(g))+1
else usedcheck15
#ifdef SINGLEFLUID
        N(x,y,z)%n_r(bounce(g))=N(x,y,z)%n_r(bounce(g))+(1-acccoef)*Nbuffer(x+dir2(1),y+dir2(2),dir2(3))%n_r(slip_dir)
#endif
#ifndef SINGLEFLUID
        N(x,y,z)%n_r(bounce(g))=N(x,y,z)%n_r(bounce(g))+(1-acccoef)*Nbuf(x+dir2(1),y+dir2(2),dir2(3))%n_r(slip_dir)
        N(x,y,z)%n_b(bounce(g))=N(x,y,z)%n_b(bounce(g))+(1-acccoef)*Nbuf(x+dir2(1),y+dir2(2),dir2(3))%n_b(slip_dir)
#endif
#ifndef NOSURFACTANT
        N(x,y,z)%n_s(bounce(g))=N(x,y,z)%n_s(bounce(g))+(1-acccoef)*Nbuf(x+dir2(1),y+dir2(2),dir2(3))%n_s(slip_dir)
        d_adv(:,bounce(g),x,y,z)=d_adv(:,bounce(g),x,y,z)+(1-acccoef)*Nbuf(xa,ya,cz(g))%d
#endif
useddirections(x,y,z)%n(bounce(g)) = useddirections(x,y,z)%n(bounce(g))+1
end if usedcheck15


usedcheck16: if(useddirections(x,y,z)%n(bounce(index1))==0)then
#ifdef SINGLEFLUID
        N(x,y,z)%n_r(bounce(index1))=0.5*acccoef*Nbuffer(x+dir2(1),y+dir2(2),dir2(3))%n_r(slip_dir)
#endif
#ifndef SINGLEFLUID
        N(x,y,z)%n_r(bounce(index1))=0.5*acccoef*Nbuf(x+dir2(1),y+dir2(2),dir2(3))%n_r(slip_dir)
        N(x,y,z)%n_b(bounce(index1))=0.5*acccoef*Nbuf(x+dir2(1),y+dir2(2),dir2(3))%n_b(slip_dir)
#endif
#ifndef NOSURFACTANT
        N(x,y,z)%n_s(bounce(index1))=0.5*acccoef*Nbuf(x+dir2(1),y+dir2(2),dir2(3))%n_s(slip_dir)
        d_adv(:,bounce(index1),x,y,z)=0.5*acccoef*Nbuf(xa,ya,cz(g))%d
#endif
useddirections(x,y,z)%n(bounce(index1)) = useddirections(x,y,z)%n(bounce(index1))+1
else usedcheck16
#ifdef SINGLEFLUID
        N(x,y,z)%n_r(bounce(index1))=N(x,y,z)%n_r(bounce(index1))+0.5*acccoef*Nbuffer(x+dir2(1),y+dir2(2),dir2(3))%n_r(slip_dir)
#endif
#ifndef SINGLEFLUID
        N(x,y,z)%n_r(bounce(index1))=N(x,y,z)%n_r(bounce(index1))+0.5*acccoef*Nbuf(x+dir2(1),y+dir2(2),dir2(3))%n_r(slip_dir)
        N(x,y,z)%n_b(bounce(index1))=N(x,y,z)%n_b(bounce(index1))+0.5*acccoef*Nbuf(x+dir2(1),y+dir2(2),dir2(3))%n_b(slip_dir)
#endif
#ifndef NOSURFACTANT
        N(x,y,z)%n_s(bounce(index1))=N(x,y,z)%n_s(bounce(index1))+0.5*acccoef*Nbuf(x+dir2(1),y+dir2(2),dir2(3))%n_s(slip_dir)
        d_adv(:,bounce(index1),x,y,z)=d_adv(:,bounce(index1),x,y,z)+0.5*acccoef*Nbuf(xa,ya,cz(g))%d
#endif
useddirections(x,y,z)%n(bounce(index1)) = useddirections(x,y,z)%n(bounce(index1))+1
end if usedcheck16



elseif((acccoef .gt. 1 ).and.(acccoef .le. 2))then        acccoefcheck4
        acccoef2 = 2 - acccoef

usedcheck17: if(useddirections(x,y,z)%n(bounce(g))==0)then
#ifdef SINGLEFLUID
        N(x,y,z)%n_r(bounce(g))=0
#endif
#ifndef SINGLEFLUID
        N(x,y,z)%n_r(bounce(g))=0
        N(x,y,z)%n_b(bounce(g))=0
#endif
#ifndef NOSURFACTANT
        N(x,y,z)%n_s(bounce(g))=0
        d_adv(:,bounce(g),x,y,z)=(/0,0,0/)
#endif
useddirections(x,y,z)%n(bounce(g)) = useddirections(x,y,z)%n(bounce(g))+1
else usedcheck17
useddirections(x,y,z)%n(bounce(g)) = useddirections(x,y,z)%n(bounce(g))+1
end if usedcheck17


usedcheck18: if(useddirections(x,y,z)%n(bounce(index1))==0)then
#ifdef SINGLEFLUID
        N(x,y,z)%n_r(bounce(index1))=0
#endif
#ifndef SINGLEFLUID
        N(x,y,z)%n_r(bounce(index1))=0
        N(x,y,z)%n_b(bounce(index1))=0
#endif
#ifndef NOSURFACTANT
        N(x,y,z)%n_s(bounce(index1))=0
        d_adv(:,bounce(index1),x,y,z)=0
#endif
useddirections(x,y,z)%n(bounce(index1)) = useddirections(x,y,z)%n(bounce(index1))+1
else usedcheck18
useddirections(x,y,z)%n(bounce(index1)) = useddirections(x,y,z)%n(bounce(index1))+1
end if usedcheck18


end if acccoefcheck4

else if(N(x+dir1(1),y+dir1(2),z+dir1(3))%rock_state==0 .and. N(x+dir2(1),y+dir2(2),z+dir2(3))%rock_state==0)then getgeometry
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Edge
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

usedchecklast: if(useddirections(x,y,z)%n(bounce(g))==0)then
#ifdef SINGLEFLUID
        N(x,y,z)%n_r(bounce(g))=Nbuffer(x,y,0)%n_r(g)
#endif
#ifndef SINGLEFLUID
        N(x,y,z)%n_r(bounce(g))=Nbuf(x,y,0)%n_r(g)
        N(x,y,z)%n_b(bounce(g))=Nbuf(x,y,0)%n_b(g)
#endif
#ifndef NOSURFACTANT
        N(x,y,z)%n_s(bounce(g))=Nbuf(x,y,0)%n_s(g)
                          d_adv(:,bounce(g),x,y,z)=Nbuf(xa,ya,cz(g))%d
#endif

useddirections(x,y,z)%n(bounce(g)) = useddirections(x,y,z)%n(bounce(g))+1
else usedchecklast

#ifdef SINGLEFLUID
        N(x,y,z)%n_r(bounce(g))=N(x,y,z)%n_r(bounce(g))+Nbuffer(x,y,0)%n_r(g)
#endif
#ifndef SINGLEFLUID
        N(x,y,z)%n_r(bounce(g))=N(x,y,z)%n_r(bounce(g))+Nbuf(x,y,0)%n_r(g)
        N(x,y,z)%n_b(bounce(g))=N(x,y,z)%n_b(bounce(g))+Nbuf(x,y,0)%n_b(g)
#endif
#ifndef NOSURFACTANT
        N(x,y,z)%n_s(bounce(g))=N(x,y,z)%n_s(bounce(g))+Nbuf(x,y,0)%n_s(g)
                          d_adv(:,bounce(g),x,y,z)=d_adv(:,bounce(g),x,y,z)+Nbuf(xa,ya,cz(g))%d
#endif
useddirections(x,y,z)%n(bounce(g)) = useddirections(x,y,z)%n(bounce(g))+1
end if usedchecklast


end if getgeometry

end if reflectaxes

end if rockdirectioncheck

end do directionloop

end if rockxyzcheck

#endif
!endif of the ifdef LOCALBC

        !end of BC selection
        endif


end subroutine lbe_boundary_condition


#ifdef LOCALBC
subroutine lbe_get_direction(g,dir1,dir2,index1,index2)
integer :: i,g,index1,index2
integer, dimension(1:3) :: dir1, dir2
if(cx(g).ne.0 .and. cy(g).ne.0)then
dir1(:)= (/ cx(g),0,0 /)
dir2(:)= (/ 0,cy(g),0 /)

do i=1,6
if(cx(i)==cx(g))index1=i
if(cy(i)==cy(g))index2=i
end do
end if

if(cx(g).ne.0 .and. cz(g).ne.0)then
dir1(:)= (/ cx(g),0,0 /)
dir2(:)= (/ 0,0,cz(g) /)
do i=1,6
if(cx(i)==cx(g))index1=i
if(cz(i)==cz(g))index2=i
end do
end if

if(cy(g).ne.0 .and. cz(g).ne.0)then
dir1(:)= (/ 0,cy(g),0 /)
dir2(:)= (/ 0,0,cz(g) /)
do i=1,6
if(cy(i)==cy(g))index1=i
if(cz(i)==cz(g))index2=i
end do
end if

end subroutine lbe_get_direction


subroutine lbe_get_direction2(dir_rock,dir_free,direction)
integer :: i,direction
integer, dimension(1:3) :: dir_rock,dir_free,dir_sum
dir_sum = dir_rock-dir_free
do i=1,nnonrest
if(cx(i)==dir_sum(1).and.cy(i)==dir_sum(2).and.cz(i)==dir_sum(3))direction=i
end do

end subroutine lbe_get_direction2
#endif

!> Calculates the forces at a site (x,y,z) and writes them to f_b,f_r[,f_s]
#ifdef SINGLEFLUID
subroutine lbe_calculate_sc_forces(N,x,y,z,f_r)
#else
#ifdef NOSURFACTANT
subroutine lbe_calculate_sc_forces(N,x,y,z,f_b,f_r)
#endif
#ifndef NOSURFACTANT
subroutine lbe_calculate_sc_forces(N,x,y,z,f_b,f_r,f_s)
#endif
#endif
  type(lbe_site),dimension(0:,0:,0:),intent(in) :: N
  integer,intent(in) :: x, y, z
  integer :: xa, ya, za, s, ss
  ! Nett forces
#ifndef SINGLEFLUID
  real*8, dimension(1:,1:,1:,1:),intent(out) :: f_b
#endif
  real*8, dimension(1:,1:,1:,1:),intent(out) :: f_r
#ifndef NOSURFACTANT
  real*8, dimension(1:,1:,1:,1:),intent(out) :: f_s

  real*8, dimension(3) :: da
  real*8 :: d_dotc, di_dot_d, di_dotc
#endif

  real*8,  dimension(3) :: f_rr, f_bb
  real*8,  dimension(3) :: f_br, f_rb
  real*8,  dimension(3) :: f_cs, f_sc

  !Independent wall interactions
  real*8,  dimension(3) :: f_wr, f_wb

#ifndef NOSURFACTANT
  real*8,  dimension(3) :: f_ss
  real*8 :: psi_s
#endif
  real*8 :: psi_b, psi_r

  !Independent wall interactions
  real*8 :: psi_wr, psi_wb
  real*8 :: rhowall

  ! Green's Functions
  real*8 green_b, green_r
  real*8 greenc_s, green_ss

  ! Independent wall interactions
  real*8 :: green_wr, green_wb

#ifdef BUGGYIFORT11
  character(len=256) :: bugbuf
#endif

#ifndef NOSURFACTANT
  f_ss = 0.d0
  psi_s = 0.d0
#endif
  psi_b = 0.d0
  psi_r = 0.d0
  f_br = 0.d0
  f_rb = 0.d0
  f_cs = 0.d0
  f_sc = 0.d0
  f_rr = 0.d0
  f_bb = 0.d0

  ! Independent wall interactions
  f_wr = 0.d0
  f_wb = 0.d0

  f_r(x,y,z,:) = 0.0_rk
#ifndef SINGLEFLUID
  f_b(x,y,z,:) = 0.0_rk
#ifndef NOSURFACTANT
  f_s(x,y,z,:) = 0.0_rk
#endif
#endif

  do s=1,nnonrest
     xa = x + cx(s)
     ya = y + cy(s)
     za = z + cz(s)

     if ( is_colloid( N(xa,ya,za)%rock_state ) ) then ! colloid node
        psi_r = N(xa,ya,za)%n_r(restvec)
        psi_wr = 0.0d0
#ifndef SINGLEFLUID
        psi_b = N(xa,ya,za)%n_b(restvec)
        psi_wb = 0.0d0
#ifndef NOSURFACTANT
        psi_s = N(xa,ya,za)%n_s(restvec)
#endif
#endif
     else if ( is_wall( N(xa,ya,za)%rock_state ) ) then ! wall node
        psi_r = N(xa,ya,za)%n_r(restvec)
        if (rock_colour_double) then
           psi_wr = N(xa,ya,za)%rock_colour_r
        else if (g_value_double) then 
           psi_wr = N(xa,ya,za)%rock_colour
        else if (gw_double_wet) then 
            psi_wr = psi_r
        else 
            psi_wr = N(xa,ya,za)%rock_colour
        end if
#ifndef SINGLEFLUID 
        psi_b = N(xa,ya,za)%n_b(restvec)
        if (rock_colour_double) then 
           psi_wb = N(xa,ya,za)%rock_colour_b
        else if (g_value_double) then
           psi_wb = N(xa,ya,za)%rock_colour
       else if (gw_double_wet) then
           psi_wb = psi_b
       else 
           psi_wb = N(xa,ya,za)%rock_colour
        end if
#ifndef NOSURFACTANT
        psi_s = N(xa,ya,za)%n_s(restvec)
#endif
#endif
!! sudhir: somthing wrong, caution ya is local to mpi node
#ifdef AXISYM ! neutral wetting
     else if ( ya .eq. 1) then ! density fluid  node
        psi_r = sum( N(xa,ya+1,za)%n_r(:)*g(:) )
        psi_wr = 0.d0
#ifndef SINGLEFLUID
        psi_b = sum( N(xa,ya+1,za)%n_b(:)*g(:) )
        psi_wb = 0.d0
#ifndef NOSURFACTANT
        psi_s = sum( N(xa,ya+1,za)%n_s(:)*g(:) )
#endif
#endif
     else if (ya .eq. ny) then ! density fluid  node
        psi_r = sum( N(xa,ya-1,za)%n_r(:)*g(:) )
        psi_wr = 0.d0
#ifndef SINGLEFLUID
        psi_b = sum( N(xa,ya-1,za)%n_b(:)*g(:) )
        psi_wb = 0.d0
#ifndef NOSURFACTANT
        psi_s = sum( N(xa,ya-1,za)%n_s(:)*g(:) )
#endif
#endif
#endif ! AXISYM

     else  ! density fluid  node
        psi_r = sum( N(xa,ya,za)%n_r(:)*g(:) )
        psi_wr = 0.d0
#ifndef SINGLEFLUID
        psi_b = sum( N(xa,ya,za)%n_b(:)*g(:) )
        psi_wb = 0.d0
#ifndef NOSURFACTANT
        psi_s = sum( N(xa,ya,za)%n_s(:)*g(:) )
#endif
#endif
     end if

     select case(psifunc)
     case (0)
        ! Original code
        ! Another funny clipping routine:-
        ! extracted from eff_mass() in Hudong's code.
        psi_r = min(1.d0, dble(psi_r))
#ifndef SINGLEFLUID
        psi_b = min(1.d0, dble(psi_b))
#ifndef NOSURFACTANT
        psi_s = min(1.d0, dble(psi_s))
#endif
#endif

        ! Independent wall interactions
        ! For now both fluid components interact with 'red wall'
        psi_wr = min(1.d0, dble(psi_wr))
#ifndef SINGLEFLUID
        psi_wb = min(1.d0, dble(psi_wb))
#endif
     case (1)
        ! psi = n, with no clipping
        ! nothing further to do.
     case (2)
        ! psi = 1 - exp(-n)
        ! no clipping
        psi_r = 1.d0 - exp(-psi_r)
#ifndef SINGLEFLUID
        psi_b = 1.d0 - exp(-psi_b)
#ifndef NOSURFACTANT
        psi_s = 1.d0 - exp(-psi_s)
#endif
#endif
        !Independent wall interactions
        !For now both fluid components interact with 'red wall'
        if (g_value_double) then
           ! do nothing, psi_wr = psi_wr
        else
           psi_wr = 1.d0 - exp(-psi_wr)
#ifndef SINGLEFLUID
           psi_wb = 1.d0 - exp(-psi_wb)
#endif
        end if
     case(3)
        ! Carnahan-Starling EOS
        ! psi = sqrt(((2*psi*Tcs*(1+psi+psi**2-psi**3)/(1-psi)**3)-4/3*psi**2)/(6*g_br))
        psi_r = CS_EOS(psi_r, g_br)
#ifndef SINGLEFLUID
        psi_b = CS_EOS(psi_b, g_br)
#ifndef NOSURFACTANT
        psi_s = CS_EOS(psi_s, g_br)
#endif
#endif
        !Independent wall interactions
        !For now both fluid components interact with 'red wall'
        psi_wr = CS_EOS(psi_wr, g_wr)
#ifndef SINGLEFLUID
        psi_wb = CS_EOS(psi_wb, g_wb)
#endif

     case default
        call log_msg("ERROR: Unknown psi functional, aborting ...")
        call Abend
     end select

     !Now, psi_s = psi^s(x+c_i). Calculate the Green's Functions.
     ! This is where referring to the surfactants as "s" rather
     ! than "g" becomes an advantage.
#ifndef SINGLEFLUID
     green_b = psi_r*g(s)
     green_r = psi_b*g(s)
     ! f += green_b * c_i
     ! ^              ^

     f_rb = f_rb + green_r*c(s,:)
     f_br = f_br + green_b*c(s,:)

     !Independent wall interactions
     green_wr = psi_wr*g(s)
     f_wr = f_wr + green_wr*c(s,:)
     green_wb = psi_wb*g(s)
     f_wb = f_wb + green_wb*c(s,:)

#ifdef MCMP
     f_rr = f_rr + green_r*c(s,:)
     f_bb = f_bb + green_b*c(s,:)
#endif
#else
     if (SCMP) then
        f_rr = f_rr + psi_r*g(s)*c(s,:)
        green_wr = psi_wr*g(s)
        f_wr = f_wr + green_wr*c(s,:)
     end if
#endif
     ! That's the colour-colour forces done for now.

#ifndef NOSURFACTANT
     ! Now think about the force on colours due to surfactants.
     ! da = surf direction vector = d^eq(x + c_i)*psi_s(x + c_i)

     da = N(xa,ya,za)%d*psi_s

     ! Dot with the vector c_i
     ! Remember that d (c c) == (d.c) c
     !               ^  ^ ^      ^ ^  ^

     d_dotc = sum(da*c(s,:))
     f_cs = f_cs + (da - D/2.d0 *d_dotc*c(s,:) )*g(s)

     ! That's the force on colours due to surfs done.
     ! Now think about the force on surfs due to colours.
     ! Dot c_i with the direction of the surfactants.

     di_dotc = sum( N(x,y,z)%d*c(s,:) )
     greenc_s = (psi_b - psi_r )*g(s)
     f_sc = f_sc + greenc_s*(N(x,y,z)%d-D/2.*di_dotc*c(s,:))

     ! Now think about the surf-surf forces.

     di_dot_d = sum(da*N(x,y,z)%d)
     green_ss = 2.d0*D*g(s)
     f_ss = f_ss + green_ss*( &
          + (di_dot_d - D/2.d0*di_dotc*d_dotc)*c(s,:)&
          + da*di_dotc + d_dotc*N(x,y,z)%d )
#endif
  end do ! s over nonrest

#ifndef SINGLEFLUID
  f_br = -g_br*f_br
  f_rb = -g_br*f_rb
#ifdef MCMP
  f_rr = -g_rr*f_rr
  f_bb = -g_bb*f_bb
#endif
  
  ! Independent wall interactions
  if (gw_wet_global) then 
  ! do nothing
  else if (g_value_double .or. gw_double_wet) then
     g_wr = N(xa,ya,za)%rock_colour_r ! rock_colour_r is equal to g_wr.
     g_wb = N(xa,ya,za)%rock_colour_b ! rock_colour_b is equal to g_wb.
 ! Test for Dennis .h5 rock files
 !if (xa==64.AND.ya ==64 .AND. za < 3) then
 !  print*, "g_wr,g_wb \n", za, g_wr, g_wb
 ! READ*
 !end if 
  end if
  f_wr = -g_wr*f_wr
  f_wb = -g_wb*f_wb
  
#ifndef NOSURFACTANT
  f_sc = 2.d0*g_bs*f_sc
  f_cs = -2.d0*g_bs*f_cs
  f_ss = -g_ss*f_ss
  
#endif
#else
  if (SCMP) then
     f_rr = -g_rr*f_rr
     f_wr = -g_wr*f_wr
  end if
#endif
  
  ! Find the nett force on each particle species at the site (x,y,z)
  select case(psifunc)
  case(0)
#ifndef SINGLEFLUID
     ! Original code
     f_b(x,y,z,:) = f_br + f_cs
     f_r(x,y,z,:) = f_rb - f_cs
#else
     if (SCMP) then
        f_r(x,y,z,:) = f_r(x,y,z,:) + f_rr
     end if
#endif
     
     f_r(x,y,z,:) = f_r(x,y,z,:) + f_wr/get_tau_r(N, x, y, z)*tau_wr
#ifndef SINGLEFLUID
     f_b(x,y,z,:) = f_b(x,y,z,:) + f_wb/get_tau_b(N, x, y, z)*tau_wb
#endif
     
#ifdef MCMP
     f_b(x,y,z,:) = f_b(x,y,z,:) + f_bb
     f_r(x,y,z,:) = f_r(x,y,z,:) + f_rr
#endif
     
#ifndef NOSURFACTANT
     f_s(x,y,z,:) = f_sc + f_ss
#endif
     
  case(1)
#ifndef SINGLEFLUID
     ! Original psi=n, but no clipping
     f_b(x,y,z,:) = (f_br + f_cs) * sum(N(x,y,z)%n_b(:)*g(:))
     f_r(x,y,z,:) = (f_rb - f_cs) * sum(N(x,y,z)%n_r(:)*g(:))
#else
     if (SCMP) then
        f_r(x,y,z,:) = f_rr*sum(N(x,y,z)%n_r(:)*g(:))
     end if
#endif

     !Independent wall interactions
       f_r(x,y,z,:) = f_r(x,y,z,:) + f_wr/get_tau_r(N, x, y, z)*tau_wr*sum(N(x,y,z)%n_r(:)*g(:))
#ifndef SINGLEFLUID
       f_b(x,y,z,:) = f_b(x,y,z,:) + f_wb/get_tau_b(N, x, y, z)*tau_wb*sum(N(x,y,z)%n_b(:)*g(:))
#endif

#ifdef MCMP
     f_b(x,y,z,:) = f_b(x,y,z,:) + f_bb*sum(N(x,y,z)%n_b(:)*g(:))
     f_r(x,y,z,:) = f_r(x,y,z,:) + f_rr*sum(N(x,y,z)%n_r(:)*g(:))
#endif

#ifndef NOSURFACTANT
     f_s(x,y,z,:) = (f_sc + f_ss) * sum(N(x,y,z)%n_s(:)*g(:))
#endif

  case(2)
#ifndef SINGLEFLUID
     ! psi = 1-exp(-n), no clipping
     f_b(x,y,z,:) = (f_br + f_cs) * (1.d0 - exp(sum(-N(x,y,z)%n_b(:)*g(:))))
     f_r(x,y,z,:) = (f_rb - f_cs) * (1.d0 - exp(sum(-N(x,y,z)%n_r(:)*g(:))))
#else
     if (SCMP) then
        f_r(x,y,z,1) = f_rr(1)*(1.d0 - exp(sum(-N(x,y,z)%n_r(:)*g(:))))
        f_r(x,y,z,2) = f_rr(2)*(1.d0 - exp(sum(-N(x,y,z)%n_r(:)*g(:))))
        f_r(x,y,z,3) = f_rr(3)*(1.d0 - exp(sum(-N(x,y,z)%n_r(:)*g(:))))
     end if
#endif

#ifdef BUGGYIFORT11
     write(bugbuf,"(3(E10.2,X))") f_r(x,y,z,:)
#endif

       f_r(x,y,z,:) = f_r(x,y,z,:) &
          + f_wr/get_tau_r(N, x, y, z)*tau_wr*(1.d0 - exp(sum(-N(x,y,z)%n_r(:)*g(:))))
#ifndef SINGLEFLUID
       f_b(x,y,z,:) = f_b(x,y,z,:) &
          + f_wb/get_tau_b(N, x, y, z)*tau_wb*(1.d0 - exp(sum(-N(x,y,z)%n_b(:)*g(:))))
#endif

#ifdef MCMP
     f_b(x,y,z,:) = f_b(x,y,z,:) + f_bb*(1.d0 - exp(sum(-N(x,y,z)%n_b(:)*g(:))))
     f_r(x,y,z,:) = f_r(x,y,z,:) + f_rr*(1.d0 - exp(sum(-N(x,y,z)%n_r(:)*g(:))))
#endif

#ifndef NOSURFACTANT
     f_s(x,y,z,:) = (f_sc + f_ss) * (1.d0 - exp(sum(-N(x,y,z)%n_s(:)*g(:))))
#endif

case(3)
#ifndef SINGLEFLUID
     ! psi = 1-exp(-n), no clipping
     f_b(x,y,z,:) = (f_br + f_cs) * CS_EOS(sum(-N(x,y,z)%n_b(:)*g(:)), g_br)
     f_r(x,y,z,:) = (f_rb - f_cs) * CS_EOS(sum(-N(x,y,z)%n_r(:)*g(:)), g_br)
#else
     if (SCMP) then
        f_r(x,y,z,1) = f_rr(1)*CS_EOS(sum(-N(x,y,z)%n_r(:)*g(:)),g_rr)
        f_r(x,y,z,2) = f_rr(2)*CS_EOS(sum(-N(x,y,z)%n_r(:)*g(:)),g_rr)
        f_r(x,y,z,3) = f_rr(3)*CS_EOS(sum(-N(x,y,z)%n_r(:)*g(:)),g_rr)
     end if
#endif

#ifdef BUGGYIFORT11
     write(bugbuf,"(3(E10.2,X))") f_r(x,y,z,:)
#endif

       f_r(x,y,z,:) = f_r(x,y,z,:) &
          + f_wr/get_tau_r(N, x, y, z)*tau_wr*CS_EOS(sum(-N(x,y,z)%n_r(:)*g(:)),g_wr)
#ifndef SINGLEFLUID
       f_b(x,y,z,:) = f_b(x,y,z,:) &
          + f_wb/get_tau_b(N, x, y, z)*tau_wb*CS_EOS(sum(-N(x,y,z)%n_b(:)*g(:)),g_wb)
#endif

#ifdef MCMP
       f_b(x,y,z,:) = f_b(x,y,z,:) + f_bb*CS_EOS(sum(-N(x,y,z)%n_b(:)*g(:)),g_br)
       f_r(x,y,z,:) = f_r(x,y,z,:) + f_rr*CS_EOS(sum(-N(x,y,z)%n_r(:)*g(:)),g_br)
#endif


     ! Not yet work for surfactant
     !#ifndef NOSURFACTANT
     !     f_s(x,y,z,:) = (f_sc + f_ss)*CS_EOS(sum(-N(x,y,z)%n_s(:)*g(:)),g_sc)
     !#endif
     
  case default
     call log_msg("ERROR: Unknown psi functional, aborting ...")
     call Abend
  end select

end subroutine lbe_calculate_sc_forces

! Calculates the forces at a site (x,y,z) and writes them to
! f_b,f_r[,f_s] .  Does the same as lbe_calculate_sc_forces() but
! operates on the whole lattice (whole_N instead of N in
! lbe.F90). Thus, f_r, f_b, and f_s need to contain halo nodes as well
! and halo_extent needs to be 2 or more. In consequence it is possible
! with this routine to calculate SC forces for lattice sites within a
! halo of extent 1.
!
#ifdef SINGLEFLUID
subroutine md_calculate_sc_forces(N,x,y,z,f_r)
#else
#ifdef NOSURFACTANT
subroutine md_calculate_sc_forces(N,x,y,z,f_b,f_r)
#endif
#ifndef NOSURFACTANT
subroutine md_calculate_sc_forces(N,x,y,z,f_b,f_r,f_s)
#endif
#endif
  type(lbe_site), dimension(1-halo_extent:,1-halo_extent:,1-halo_extent:) :: N
  integer,intent(in) :: x, y, z
  integer :: xa, ya, za, s, ss
  ! Nett forces
#ifndef SINGLEFLUID
  real*8, dimension(0:nx+1,0:ny+1,0:nz+1,3) :: f_b
#endif
  real*8, dimension(0:nx+1,0:ny+1,0:nz+1,3) :: f_r
#ifndef NOSURFACTANT
  real*8, dimension(0:nx+1,0:ny+1,0:nz+1,3) :: f_s

  real*8,dimension(3) :: da
  real*8 :: d_dotc, di_dot_d, di_dotc
#endif

  real*8,  dimension(3) :: f_rr, f_bb
  real*8,  dimension(3) :: f_br, f_rb
  real*8,  dimension(3) :: f_cs, f_sc

  !Independent wall interactions
  real*8,  dimension(3) :: f_wr, f_wb

#ifndef NOSURFACTANT
  real*8,  dimension(3) :: f_ss
  real*8 :: psi_s
#endif
 real*8 :: psi_b, psi_r

 !Independent wall interactions
 real*8 :: psi_wr, psi_wb
 real*8 :: rhowall

 ! Green's Functions
 real*8 green_b, green_r
 real*8 greenc_s, green_ss

!Independent wall interactions
 real*8 :: green_wr, green_wb

#ifndef NOSURFACTANT
  f_ss = 0.d0
  psi_s = 0.d0
#endif
  psi_b = 0.d0
  psi_r = 0.d0
  f_br = 0.d0
  f_rb = 0.d0
  f_cs = 0.d0
  f_sc = 0.d0
  f_rr = 0.d0
  f_bb = 0.d0

  !Independent wall interactions
  f_wr = 0.d0
  f_wb = 0.d0

  do s=1,nnonrest
     xa = x + cx(s)
     ya = y + cy(s)
     za = z + cz(s)

     if ( is_colloid( N(xa,ya,za)%rock_state ) ) then ! colloid node
        psi_r = N(xa,ya,za)%n_r(restvec)
        psi_wr = 0.0d0
#ifndef SINGLEFLUID
        psi_b = N(xa,ya,za)%n_b(restvec)
        psi_wb = 0.0d0
#ifndef NOSURFACTANT
        psi_s = N(xa,ya,za)%n_s(restvec)
#endif
#endif
     else if ( is_wall( N(xa,ya,za)%rock_state ) ) then ! wall node
        psi_r = N(xa,ya,za)%n_r(restvec)
       if (rock_colour_double) then
           psi_wr = N(xa,ya,za)%rock_colour_r
        else if (g_value_double) then 
           psi_wr = N(xa,ya,za)%rock_colour
        else if (gw_double_wet) then 
            psi_wr = psi_r
        else 
            psi_wr = N(xa,ya,za)%rock_colour
        end if
     
#ifndef SINGLEFLUID
        psi_b = N(xa,ya,za)%n_b(restvec)
      if (rock_colour_double) then 
           psi_wb = N(xa,ya,za)%rock_colour_b
        else if (g_value_double) then
           psi_wb = N(xa,ya,za)%rock_colour
       else if (gw_double_wet) then
           psi_wb = psi_b
       else 
           psi_wb = N(xa,ya,za)%rock_colour
        end if
#ifndef NOSURFACTANT
        psi_s = N(xa,ya,za)%n_s(restvec)
#endif
#endif
     else  ! fluid  node
        psi_r = sum( N(xa,ya,za)%n_r(:)*g(:) )
        psi_wr = 0.d0
#ifndef SINGLEFLUID
        psi_b = sum( N(xa,ya,za)%n_b(:)*g(:) )
        psi_wb = 0.d0
#ifndef NOSURFACTANT
        psi_s = sum( N(xa,ya,za)%n_s(:)*g(:) )
#endif
#endif
     end if

     select case(psifunc)
     case (0)
        ! Original code

        ! Another funny clipping routine:-
        ! extracted from eff_mass() in Hudong's code.
        psi_r = min(1.d0, dble(psi_r))
#ifndef SINGLEFLUID
        psi_b = min(1.d0, dble(psi_b))
#ifndef NOSURFACTANT
        psi_s = min(1.d0, dble(psi_s))
#endif
#endif

        !Independent wall interactions
        !For now both fluid components interact with 'red wall'
        psi_wr = min(1.d0, dble(psi_wr))
#ifndef SINGLEFLUID
        psi_wb = min(1.d0, dble(psi_wb))
#endif
     case (1)
        ! psi = n, with no clipping
        ! nothing further to do.
     case (2)
        ! psi = 1 - exp(-n)
        ! no clipping
        psi_r = 1.d0 - exp(-psi_r)
#ifndef SINGLEFLUID
        psi_b = 1.d0 - exp(-psi_b)
#ifndef NOSURFACTANT
        psi_s = 1.d0 - exp(-psi_s)
#endif
#endif
        !Independent wall interactions
        !For now both fluid components interact with 'red wall'
if (g_value_double) then
           ! do nothing, psi_wr = psi_wr
        else
        psi_wr = 1.d0 - exp(-psi_wr)
#ifndef SINGLEFLUID
        psi_wb = 1.d0 - exp(-psi_wb)
#endif
end if
     case (3)
        ! Carnahan-Starling EOS
        ! psi = sqrt(((2*psi*Tcs*(1+psi+psi**2-psi**3)/(1-psi)**3)-4/3*psi**2)/(6*g_br))
        psi_r = CS_EOS(psi_r, g_br)
#ifndef SINGLEFLUID
        psi_b = CS_EOS(psi_b, g_br)
#ifndef NOSURFACTANT
        psi_s = CS_EOS(psi_s, g_bs) 
#endif
#endif
        !Independent wall interactions
        !For now both fluid components interact with 'red wall'
        psi_wr = CS_EOS(psi_wr, g_wr)
#ifndef SINGLEFLUID
        psi_wb = CS_EOS(psi_wb, g_wb) 
#endif

     case default
        call log_msg("ERROR: Unknown psi functional, aborting ...")
        call Abend
     end select

      !Now, psi_s = psi^s(x+c_i). Calculate the Green's Functions.
      ! This is where referring to the surfactants as "s" rather
      ! than "g" becomes an advantage.
#ifndef SINGLEFLUID
      green_b = psi_r*g(s)
      green_r = psi_b*g(s)
      ! f += green_b * c_i
      ! ^              ^

      f_br = f_br + green_b*c(s,:)
      f_rb = f_rb + green_r*c(s,:)

     !Independent wall interactions
     green_wr = psi_wr*g(s)
     f_wr = f_wr + green_wr*c(s,:)
     green_wb = psi_wb*g(s)
     f_wb = f_wb + green_wb*c(s,:)

#ifdef MCMP
     f_rr = f_rr + green_r*c(s,:)
     f_bb = f_bb + green_b*c(s,:)
#endif
#else
     if (SCMP) then
        f_rr = f_rr + psi_r*g(s)*c(s,:)
        green_wr = psi_wr*g(s)
        f_wr = f_wr + green_wr*c(s,:)
     end if
#endif
     ! That's the colour-colour forces done for now.

#ifndef NOSURFACTANT
     ! Now think about the force on colours due to surfactants.
     ! da = surf direction vector = d^eq(x + c_i)*psi_s(x + c_i)

     da = N(xa,ya,za)%d*psi_s

     ! Dot with the vector c_i
     ! Remember that d (c c) == (d.c) c
     !               ^  ^ ^      ^ ^  ^

     d_dotc = sum(da*c(s,:))
     f_cs = f_cs + (da - D/2.d0 *d_dotc*c(s,:) )*g(s)

     ! That's the force on colours due to surfs done.
     ! Now think about the force on surfs due to colours.
     ! Dot c_i with the direction of the surfactants.

     di_dotc = sum( N(x,y,z)%d*c(s,:) )
     greenc_s = (psi_b - psi_r )*g(s)
     f_sc = f_sc + greenc_s*(N(x,y,z)%d-D/2.*di_dotc*c(s,:))

     ! Now think about the surf-surf forces.

     di_dot_d = sum(da*N(x,y,z)%d)
     green_ss = 2.d0*D*g(s)
     f_ss = f_ss + green_ss*( &
          + (di_dot_d - D/2.d0*di_dotc*d_dotc)*c(s,:) &
          + da*di_dotc + d_dotc*N(x,y,z)%d )
#endif
  end do ! s over nonrest

#ifndef SINGLEFLUID
  f_br = -g_br*f_br
  f_rb = -g_br*f_rb
#ifdef MCMP
  f_rr = -g_rr*f_rr
  f_bb = -g_bb*f_bb
#endif

 ! Independent wall interactions
if (gw_wet_global) then
    ! do nothing
    else if (g_value_double .or. gw_double_wet) then
       g_wr = N(xa,ya,za)%rock_colour_r ! rock_colour_r is equal to g_wr.
       g_wb = N(xa,ya,za)%rock_colour_b ! rock_colour_b is equal to g_wb.  
end if

  f_wr = -g_wr*f_wr
  f_wb = -g_wb*f_wb

#ifndef NOSURFACTANT
  f_sc =  2.d0*g_bs*f_sc
  f_cs = -2.d0*g_bs*f_cs
  f_ss = -g_ss*f_ss
#endif
#else
  if (SCMP) then
     f_rr = -g_rr*f_rr
     f_wr = -g_wr*f_wr
  end if
#endif

  ! Find the nett force on each particle species at the site (x,y,z)
  select case(psifunc)
  case(0)
#ifndef SINGLEFLUID
     ! Original code
     f_b(x,y,z,:) = f_br + f_cs
     f_r(x,y,z,:) = f_rb - f_cs
#else
     if (SCMP) then
        f_r(x,y,z,:) = f_r(x,y,z,:) + f_rr
     end if
#endif

        f_r(x,y,z,:) = f_r(x,y,z,:) + f_wr/get_tau_r(N, x, y, z)*tau_wr
#ifndef SINGLEFLUID
        f_b(x,y,z,:) = f_b(x,y,z,:) + f_wb/get_tau_b(N, x, y, z)*tau_wb
#endif

#ifdef MCMP
     f_b(x,y,z,:) = f_b(x,y,z,:) + f_bb
     f_r(x,y,z,:) = f_r(x,y,z,:) + f_rr
#endif

#ifndef NOSURFACTANT
     f_s(x,y,z,:) = f_sc + f_ss
#endif

  case(1)
#ifndef SINGLEFLUID
     ! Original psi=n, but no clipping
     f_b(x,y,z,:) = (f_br + f_cs)*sum(N(x,y,z)%n_b(:)*g(:))
     f_r(x,y,z,:) = (f_rb - f_cs)*sum(N(x,y,z)%n_r(:)*g(:))
#else
     if (SCMP) then
        f_r(x,y,z,:) = f_rr*sum(N(x,y,z)%n_r(:)*g(:))
     end if
#endif

        f_r(x,y,z,:) = f_r(x,y,z,:) + f_wr/get_tau_r(N, x, y, z)*tau_wr*sum(N(x,y,z)%n_r(:)*g(:))
#ifndef SINGLEFLUID
        f_b(x,y,z,:) = f_b(x,y,z,:) + f_wb/get_tau_b(N, x, y, z)*tau_wb* sum(N(x,y,z)%n_b(:)*g(:))
#endif

#ifdef MCMP
     f_b(x,y,z,:) = f_b(x,y,z,:) + f_bb*sum(N(x,y,z)%n_b(:)*g(:))
     f_r(x,y,z,:) = f_r(x,y,z,:) + f_rr*sum(N(x,y,z)%n_r(:)*g(:))
#endif

#ifndef NOSURFACTANT
     f_s(x,y,z,:) = (f_sc + f_ss)*sum(N(x,y,z)%n_s(:)*g(:))
#endif

  case(2)
#ifndef SINGLEFLUID
     ! psi = 1-exp(-n), no clipping
     f_b(x,y,z,:) = (f_br + f_cs)*(1.d0 - exp(sum(-N(x,y,z)%n_b(:)*g(:))))
     f_r(x,y,z,:) = (f_rb - f_cs)*(1.d0 - exp(sum(-N(x,y,z)%n_r(:)*g(:))))
#else
     if (SCMP) then
        f_r(x,y,z,1) = f_rr(1)*(1.d0 - exp(sum(-N(x,y,z)%n_r(:)*g(:))))
        f_r(x,y,z,2) = f_rr(2)*(1.d0 - exp(sum(-N(x,y,z)%n_r(:)*g(:))))
        f_r(x,y,z,3) = f_rr(3)*(1.d0 - exp(sum(-N(x,y,z)%n_r(:)*g(:))))
     end if
#endif

        f_r(x,y,z,:) = f_r(x,y,z,:) &
             + f_wr/get_tau_r(N, x, y, z)*tau_wr*(1.d0 - exp(sum(-N(x,y,z)%n_r(:)*g(:))))
#ifndef SINGLEFLUID
        f_b(x,y,z,:) = f_b(x,y,z,:) &
             + f_wb/get_tau_b(N, x, y, z)*tau_wb*(1.d0 - exp(sum(-N(x,y,z)%n_b(:)*g(:))))
#endif

#ifdef MCMP
     f_b(x,y,z,:) = f_b(x,y,z,:) + f_bb*(1.d0 - exp(sum(-N(x,y,z)%n_b(:)*g(:))))
     f_r(x,y,z,:) = f_r(x,y,z,:) + f_rr*(1.d0 - exp(sum(-N(x,y,z)%n_r(:)*g(:))))
#endif

#ifndef NOSURFACTANT
     f_s(x,y,z,:) = (f_sc + f_ss)*(1.d0 - exp(sum(-N(x,y,z)%n_s(:)*g(:))))
#endif


case(3)
#ifndef SINGLEFLUID
  
     f_b(x,y,z,:) = (f_br + f_cs) * CS_EOS(sum(-N(x,y,z)%n_b(:)*g(:)), g_br)
     f_r(x,y,z,:) = (f_rb - f_cs) * CS_EOS(sum(-N(x,y,z)%n_r(:)*g(:)), g_br)
#else
     if (SCMP) then
        f_r(x,y,z,1) = f_rr(1)*CS_EOS(sum(-N(x,y,z)%n_r(:)*g(:)),g_rr)
        f_r(x,y,z,2) = f_rr(2)*CS_EOS(sum(-N(x,y,z)%n_r(:)*g(:)),g_rr)
        f_r(x,y,z,3) = f_rr(3)*CS_EOS(sum(-N(x,y,z)%n_r(:)*g(:)),g_rr)
     end if
#endif


       f_r(x,y,z,:) = f_r(x,y,z,:) &
          + f_wr/get_tau_r(N, x, y, z)*tau_wr*CS_EOS(sum(-N(x,y,z)%n_r(:)*g(:)),g_wr)
#ifndef SINGLEFLUID
       f_b(x,y,z,:) = f_b(x,y,z,:) &
          + f_wb/get_tau_b(N, x, y, z)*tau_wb*CS_EOS(sum(-N(x,y,z)%n_b(:)*g(:)),g_wb)
#endif

#ifdef MCMP
     f_b(x,y,z,:) = f_b(x,y,z,:) + f_bb*CS_EOS(sum(-N(x,y,z)%n_b(:)*g(:)),g_br)
     f_r(x,y,z,:) = f_r(x,y,z,:) + f_rr*CS_EOS(sum(-N(x,y,z)%n_r(:)*g(:)),g_br)
#endif

!Not yet work for surfactant 
!#ifndef NOSURFACTANT
!     f_s(x,y,z,:) = (f_sc + f_ss) * CS_EOS(sum(-N(x,y,z)%n_s(:)*g(:)),g_sc)
!#endif
  case default
    call log_msg("ERROR: Unknown psi functional, aborting ...")
    call Abend
  end select

end subroutine md_calculate_sc_forces

!> Shan-Chen forces combined with immersed boundaries
!>
!> Computes the Shan-Chen forces f_loc_* at the specified site (x, y, z).
!> This subroutine will only be called if IBM_BINARYIBM is active.

#ifdef IBM_BINARYIBM
subroutine ibm_calculate_sc_forces(N, x, y, z, f_loc_r, f_loc_b, f_loc_i)
!   type(lbe_site), dimension(0:,0:,0:), intent(in) :: N !< lattice
!   integer, intent(in) :: x, y, z !< coordinates of current site
!   real(kind=rk), dimension(3), intent(out) :: f_loc_r, f_loc_b, f_loc_i !< local force components

  !!!!!!!!!!!!!!!!!!!!!!
  !!! MODEL NUMBER 1 !!!
  !!!!!!!!!!!!!!!!!!!!!!

!   ! Declare variables.
!   integer :: i ! velocity index
!   integer :: xn, yn, zn ! coordinates of neighboring site
!   real(kind=rk) :: den_loc_r, den_loc_b, den_loc_tot ! local mass densities
!   real(kind=rk) :: den_eff_loc_r, den_eff_loc_b, den_eff_loc_i ! local effective mass densities
!   real(kind=rk) :: psi_loc_r, psi_loc_b, psi_loc_i ! local potentials
!   real(kind=rk) :: den_nb_r, den_nb_b, den_nb_tot ! neighbor mass densities
!   real(kind=rk) :: den_eff_nb_r, den_eff_nb_b, den_eff_nb_i ! neighbor effective mass densities
!   real(kind=rk) :: psi_nb_r, psi_nb_b, psi_nb_i ! neighbor potentials
!   real(kind=rk) :: index_loc, index_nb ! local and neighbor IBM index field
! 
!   ! Compute densities on current site.
!   den_loc_r = sum(N(x, y, z)%n_r(:) * g(:))
!   den_loc_b = sum(N(x, y, z)%n_b(:) * g(:))
!   den_loc_tot = den_loc_r + den_loc_b
! 
!   ! Compute effective densities on current site based on the value of the index field.
!   index_loc = get_ibm_index(x, y, z)
!   den_eff_loc_r = den_loc_r * (1.0d0 - index_loc)
!   den_eff_loc_b = den_loc_b * (1.0d0 - index_loc)
!   den_eff_loc_i = den_loc_tot * index_loc
! 
!   ! Compute potentials on current site.
!   psi_loc_r = 1.d0 - exp(-den_eff_loc_r)
!   psi_loc_b = 1.d0 - exp(-den_eff_loc_b)
!   psi_loc_i = 1.d0 - exp(-den_eff_loc_i)
! 
!   ! Reset forces.
!   f_loc_r(:) = 0.0d0
!   f_loc_b(:) = 0.0d0
!   f_loc_i(:) = 0.0d0
! 
!   ! Loop over all neighbor sites
!   do i = 1, nnonrest
!     ! Identify neighbor coordinates.
!     xn = x + cx(i)
!     yn = y + cy(i)
!     zn = z + cz(i)
! 
!     ! Compute densities on neighbor site.
!     den_nb_r = sum(N(xn, yn, zn)%n_r(:) * g(:))
!     den_nb_b = sum(N(xn, yn, zn)%n_b(:) * g(:))
!     den_nb_tot = den_nb_r + den_nb_b
! 
!     ! Compute effective densities on neighbor site based on the value of the index field.
!     if(xn >= 1 .and. xn <= nx .and. yn >= 1 .and. yn <= ny .and. zn >= 1 .and. zn <= nz) then
!       index_nb = get_ibm_index(xn, yn, zn)
!     else
!       index_nb = 0.0d0
!     end if
! 
!     den_eff_nb_r = den_nb_r * (1.0d0 - index_nb)
!     den_eff_nb_b = den_nb_b * (1.0d0 - index_nb)
!     den_eff_nb_i = den_nb_tot * index_nb
! 
!     ! Compute potentials on neighbor site.
!     psi_nb_r = 1.d0 - exp(-den_eff_nb_r)
!     psi_nb_b = 1.d0 - exp(-den_eff_nb_b)
!     psi_nb_i = 1.d0 - exp(-den_eff_nb_i)
! 
!     ! Add force contributions.
!     f_loc_r(:) = f_loc_r(:) - (g_br * psi_nb_b + g_ri * psi_nb_i) * c(i, :) * g(i)
!     f_loc_b(:) = f_loc_b(:) - (g_br * psi_nb_r + g_bi * psi_nb_i) * c(i, :) * g(i)
!     f_loc_i(:) = f_loc_i(:) - (g_ri * psi_nb_r + g_bi * psi_nb_b) * c(i, :) * g(i)
!   end do
! 
!   ! Take into account the current potentials.
!   ! The results are the total SC forces acting on the red, blue, and index components of the current site.
!   f_loc_r(:) = f_loc_r(:) * psi_loc_r
!   f_loc_b(:) = f_loc_b(:) * psi_loc_b
!   f_loc_i(:) = f_loc_i(:) * psi_loc_i

  !!!!!!!!!!!!!!!!!!!!!!
  !!! MODEL NUMBER 2 !!!
  !!!!!!!!!!!!!!!!!!!!!!

  type(lbe_site), dimension(0:,0:,0:), intent(in) :: N !< lattice
  integer, intent(in) :: x, y, z !< coordinates of current site
  real(kind=rk), dimension(3), intent(out) :: f_loc_r, f_loc_b, f_loc_i !< local force components

  ! Declare variables.
  integer :: i ! velocity index
  integer :: xn, yn, zn ! coordinates of neighboring site
  real(kind=rk) :: den_loc_r, den_loc_b, den_loc_tot ! local mass densities
  real(kind=rk) :: den_eff_loc_r, den_eff_loc_b ! local effective mass densities
  real(kind=rk) :: psi_loc_r, psi_loc_b ! local potentials
  real(kind=rk) :: psi_eff_loc_r, psi_eff_loc_b ! local effective potentials
  real(kind=rk) :: den_nb_r, den_nb_b, den_nb_tot ! neighbor mass densities
  real(kind=rk) :: den_eff_nb_r, den_eff_nb_b ! neighbor effective mass densities
  real(kind=rk) :: psi_nb_r, psi_nb_b ! neighbor potentials
  real(kind=rk) :: psi_eff_nb_r, psi_eff_nb_b ! neighbor effective potentials
  real(kind=rk) :: index_loc, index_nb ! local and neighbor IBM index field
  real(kind=rk) :: colour_loc
  real(kind=rk), dimension(3) :: force_temp

  ! Compute densities on current site.
  den_loc_r = sum(N(x, y, z)%n_r(:) * g(:))
  den_loc_b = sum(N(x, y, z)%n_b(:) * g(:))
  den_loc_tot = den_loc_r + den_loc_b
  colour_loc = den_loc_r - den_loc_b

  ! Compute local true and effective densities and potentials.
  index_loc = get_ibm_index(x, y, z)
!   index_loc = floor(get_ibm_index(x, y, z) + 0.5d0)
  den_eff_loc_r = (den_loc_r + ibm_colour * index_loc) * (1.0d0 - index_loc)
  den_eff_loc_b = (den_loc_b - ibm_colour * index_loc) * (1.0d0 - index_loc)
  psi_loc_r = 1.d0 - exp(-den_loc_r)
  psi_loc_b = 1.d0 - exp(-den_loc_b)
  psi_eff_loc_r = 1.d0 - exp(-den_eff_loc_r)
  psi_eff_loc_b = 1.d0 - exp(-den_eff_loc_b)

  ! Reset forces.
  f_loc_r(:) = 0.0d0
  f_loc_b(:) = 0.0d0
  f_loc_i(:) = 0.0d0

  ! Loop over all neighbor sites
  do i = 1, nnonrest
    ! Identify neighbor coordinates.
    xn = x + cx(i)
    yn = y + cy(i)
    zn = z + cz(i)

    ! Compute effective densities on neighbor site based on the value of the index field.
    index_nb = get_ibm_index(xn, yn, zn)
!     index_nb = floor(get_ibm_index(x, y, z) + 0.5d0)

    ! Compute neighbour true and effective densities and potentials.
    den_nb_r = sum(N(xn, yn, zn)%n_r(:) * g(:))
    den_nb_b = sum(N(xn, yn, zn)%n_b(:) * g(:))
    den_nb_tot = den_nb_r + den_nb_b
    den_eff_nb_r = (den_nb_r + ibm_colour * index_nb) * (1.0d0 - index_nb)
    den_eff_nb_b = (den_nb_b - ibm_colour * index_nb) * (1.0d0 - index_nb)

    psi_nb_r = 1.d0 - exp(-den_nb_r)
    psi_nb_b = 1.d0 - exp(-den_nb_b)
    psi_eff_nb_r = 1.d0 - exp(-den_eff_nb_r)
    psi_eff_nb_b = 1.d0 - exp(-den_eff_nb_b)

    ! Compute force contributions.
    f_loc_r(:) = f_loc_r(:) - g_br * psi_eff_loc_r * psi_eff_nb_b * c(i, :) * g(i)
    f_loc_b(:) = f_loc_b(:) - g_br * psi_eff_loc_b * psi_eff_nb_r * c(i, :) * g(i)
!     force_temp(:) = 0.d0
    force_temp(:) = - noslipcoeff * (psi_loc_b * psi_nb_b - psi_loc_r * psi_nb_r) &
                    & * c(i, :) * g(i) * index_loc * colour_loc
    f_loc_r(:) = f_loc_r(:) + force_temp(:)
    f_loc_b(:) = f_loc_b(:) - force_temp(:)
  end do

!   !!!!!!!!!!!!!!!!!!!!!!
!   !!! MODEL NUMBER 3 !!!
!   !!!!!!!!!!!!!!!!!!!!!!

!   type(lbe_site), dimension(0:,0:,0:), intent(in) :: N !< lattice
!   integer, intent(in) :: x, y, z !< coordinates of current site
!   real(kind=rk), dimension(3), intent(out) :: f_loc_r, f_loc_b, f_loc_i !< local force components
!     
!   ! Declare variables.
!   integer :: i ! velocity index
!   integer :: xn, yn, zn ! coordinates of neighboring site
!   real(kind=rk) :: den_loc_r, den_loc_b ! local mass densities
!   real(kind=rk) :: psi_loc_r, psi_loc_b ! local potentials
!   real(kind=rk) :: den_nb_r, den_nb_b ! neighbor mass densities
!   real(kind=rk) :: psi_nb_r, psi_nb_b ! neighbor potentials
!   real(kind=rk) :: index_loc, index_nb ! local and neighbor IBM index field
! 
!   ! Reset forces.
!   f_loc_r(:) = 0.0d0
!   f_loc_b(:) = 0.0d0
!   f_loc_i(:) = 0.0d0
! 
!   ! Compute densities and potentials on current site if outside.
!   index_loc = get_ibm_index(x, y, z)
!   
!   if(index_loc < 0.5d0) then
!     den_loc_r = sum(N(x, y, z)%n_r(:) * g(:))
!     den_loc_b = sum(N(x, y, z)%n_b(:) * g(:))
!     psi_loc_r = 1.d0 - exp(-den_loc_r)
!     psi_loc_b = 1.d0 - exp(-den_loc_b)
!   
!     ! Loop over all neighbor sites
!     do i = 1, nnonrest
!       ! Identify neighbor coordinates.
!       xn = x + cx(i)
!       yn = y + cy(i)
!       zn = z + cz(i)
!     
!       index_nb = get_ibm_index(xn, yn, zn)
!       
!       if(index_nb < 0.5d0) then
!         ! Compute neighbour true and effective densities and potentials.
!         den_nb_r = sum(N(xn, yn, zn)%n_r(:) * g(:))
!         den_nb_b = sum(N(xn, yn, zn)%n_b(:) * g(:))
!         psi_nb_r = 1.d0 - exp(-den_nb_r)
!         psi_nb_b = 1.d0 - exp(-den_nb_b)
!       else
!         den_nb_r = den_loc_r
!         den_nb_b = den_loc_b
!         psi_nb_r = psi_loc_r
!         psi_nb_b = psi_loc_b
!       end if
!       
!       ! Compute force contributions.
!       f_loc_r(:) = f_loc_r(:) - g_br * psi_loc_r * psi_nb_b * c(i, :) * g(i)
!       f_loc_b(:) = f_loc_b(:) - g_br * psi_loc_b * psi_nb_r * c(i, :) * g(i)
!     end do
!   else
!     ! Loop over all neighbor sites
!     do i = 1, nnonrest
!       ! Identify neighbor coordinates.
!       xn = x + cx(i)
!       yn = y + cy(i)
!       zn = z + cz(i)
!     
!       index_nb = get_ibm_index(xn, yn, zn)
!       
!       if(index_nb < 0.5d0) then
!         ! Compute neighbour true and effective densities and potentials.
!         den_nb_r = sum(N(xn, yn, zn)%n_r(:) * g(:))
!         den_nb_b = sum(N(xn, yn, zn)%n_b(:) * g(:))
!         psi_nb_r = 1.d0 - exp(-den_nb_r)
!         psi_nb_b = 1.d0 - exp(-den_nb_b)
!         den_loc_r = den_nb_r
!         den_loc_b = den_nb_b
!         psi_loc_r = psi_nb_r
!         psi_loc_b = psi_nb_b
!       
!         ! Compute force contributions.
!         f_loc_r(:) = f_loc_r(:) - g_br * psi_loc_r * psi_nb_b * c(i, :) * g(i)
!         f_loc_b(:) = f_loc_b(:) - g_br * psi_loc_b * psi_nb_r * c(i, :) * g(i)
!       end if
!     end do
!   end if

  !!!!!!!!!!!!!!!!!!!!!!
  !!! MODEL NUMBER 4 !!!
  !!!!!!!!!!!!!!!!!!!!!!

!   type(lbe_site), dimension(0:,0:,0:), intent(in) :: N !< lattice
!   integer, intent(in) :: x, y, z !< coordinates of current site
!   real(kind=rk), dimension(3), intent(out) :: f_loc_r, f_loc_b, f_loc_i !< local force components
!     
!   ! Declare variables.
!   integer :: i ! velocity index
!   integer :: xn, yn, zn ! coordinates of neighboring site
!   real(kind=rk) :: den_loc_r, den_loc_b ! local mass densities
!   real(kind=rk) :: psi_loc_r, psi_loc_b ! local potentials
!   real(kind=rk) :: den_nb_r, den_nb_b ! neighbor mass densities
!   real(kind=rk) :: psi_nb_r, psi_nb_b ! neighbor potentials
!   real(kind=rk) :: index_loc, index_nb ! local and neighbor IBM index field
! 
!   ! Reset forces.
!   f_loc_r(:) = 0.0d0
!   f_loc_b(:) = 0.0d0
!   f_loc_i(:) = 0.0d0
! 
!   ! Compute densities and potentials on current site if outside.
!   index_loc = get_ibm_index(x, y, z)
!   den_loc_r = sum(N(x, y, z)%n_r(:) * g(:))
!   den_loc_b = sum(N(x, y, z)%n_b(:) * g(:))
!   psi_loc_r = 1.d0 - exp(-den_loc_r)
!   psi_loc_b = 1.d0 - exp(-den_loc_b)
!   
!   ! Loop over all neighbor sites
!   do i = 1, nnonrest
!     ! Identify neighbor coordinates.
!     xn = x + cx(i)
!     yn = y + cy(i)
!     zn = z + cz(i)
!     
!     ! Compute neighbour true and effective densities and potentials.
!     index_nb = get_ibm_index(xn, yn, zn)
!     den_nb_r = sum(N(xn, yn, zn)%n_r(:) * g(:))
!     den_nb_b = sum(N(xn, yn, zn)%n_b(:) * g(:))
!     psi_nb_r = 1.d0 - exp(-den_nb_r)
!     psi_nb_b = 1.d0 - exp(-den_nb_b)
!       
!     if(index_loc .lt. 0.5d0 .or. index_nb .lt. 0.5d0) then
!       ! Compute force contributions.
!       f_loc_r(:) = f_loc_r(:) - g_br * psi_loc_r * psi_nb_b * c(i, :) * g(i)
!       f_loc_b(:) = f_loc_b(:) - g_br * psi_loc_b * psi_nb_r * c(i, :) * g(i)
!     end if
!   end do

  ! Take into account the current potentials.
  ! The results are the total SC forces acting on the red, blue, and index components of the current site.
!   f_loc_r(:) = f_loc_r(:) * psi_eff_loc_r
!   f_loc_b(:) = f_loc_b(:) * psi_eff_loc_b
end subroutine ibm_calculate_sc_forces
#endif

subroutine zero_force_offset(N)
  integer :: x, y, z, xa, ya, za, s
  real(kind=rk) :: f_r_zero_offset, f_b_zero_offset, f_s_zero_offset, nfc
  type(lbe_site), dimension(1-halo_extent:,1-halo_extent:,1-halo_extent:) :: N

  do y = 0,(ny+1)
     do x = 0,(nx+1)
        do z = 0,(nz+1)
           if ( is_wall(N(x,y,z)%rock_state) ) then
              f_r_zero_offset = 0.d0
#ifndef SINGLEFLUID
              f_b_zero_offset = 0.d0
#ifndef NOSURFACTANT
              f_s_zero_offset = 0.d0
#endif
#endif
              nfc = 0.0d0
              if (ZFOSdiag) then
                 do s=1,nnonrest
                    xa = x + cx(s)
                    ya = y + cy(s)
                    za = z + cz(s)
                    if ( is_fluid(N(xa,ya,za)%rock_state) ) then
                       f_r_zero_offset = f_r_zero_offset + sum(N(xa,ya,za)%n_r(:)*g(:))
#ifndef SINGLEFLUID
                       f_b_zero_offset = f_b_zero_offset + sum(N(xa,ya,za)%n_b(:)*g(:))
#ifndef NOSURFACTANT
                       f_s_zero_offset = f_s_zero_offset + sum(N(xa,ya,za)%n_s(:)*g(:))
#endif
#endif
                       nfc = nfc + 1.0d0
                    end if
                 end do
              else
                 do s=1,6
                    xa = x + cx(s)
                    ya = y + cy(s)
                    za = z + cz(s)
                    if ( is_fluid(N(xa,ya,za)%rock_state) ) then
                       f_r_zero_offset = f_r_zero_offset + sum(N(xa,ya,za)%n_r(:)*g(:))
#ifndef SINGLEFLUID
                       f_b_zero_offset = f_b_zero_offset + sum(N(xa,ya,za)%n_b(:)*g(:))
#ifndef NOSURFACTANT
                       f_s_zero_offset = f_s_zero_offset + sum(N(xa,ya,za)%n_s(:)*g(:))
#endif
#endif
                       nfc = nfc + 1.0d0
                    end if
                 end do
              end if

              if ( nfc .ne. 0 ) then
                 N(x,y,z)%n_r(restvec) = f_r_zero_offset/nfc
#ifndef SINGLEFLUID
                 N(x,y,z)%n_b(restvec) = f_b_zero_offset/nfc
#ifndef NOSURFACTANT
                 N(x,y,z)%n_s(restvec) = f_s_zero_offset/nfc
#endif
#endif
              end if
           end if
        end do
     end do
  end do
end subroutine zero_force_offset

end module lbe_collision_module
