#include "lbe.h"

!> Single-fluid advection and collision without Shan-Chen interaction
!>
!> Contains all code related to the advection and collision stages of each timestep.
!> This version is for singlefluid only and does not contain any S-C forces.

module lbe_collision_simple_module

  use lbe_bdist_module
  use lbe_collision_module
  use lbe_globals_module
  use lbe_parallel_module, only: ccoords, start, tnz
  use lbe_parms_module
  use lbe_timer_module
  use lbe_types_module, only: lbe_site
  use lbe_helper_module, only: n_sites_fluid, lbe_count_sites  
#ifdef VELCONTROL
  use lbe_helper_module, only : n_sites
#endif
#ifdef MD
  use lbe_md_module, only: ramping, fix_v, forcebalance
#endif
#ifdef LOCALBC
  use lbe_types_module, only: bc_check,lbe_site
#endif
  use lbe_log_module
  implicit none
  private
#ifdef VELCONTROL  
  include 'mpif.h'
#endif  

  public lbe_advection_simple&
#ifdef INTERPOLATEDBB
         &,lbe_interpolated_dist&
#endif
         &, lbe_interactions_simple

  type(lbe_site), dimension(:,:,:), allocatable :: Nbuf
#ifdef INTERPOLATEDBB
  type(lbe_site), dimension(:,:,:), allocatable :: lbe_Nbuf
#endif  

contains

  !> Performs advection on the particles in array \c N
  !> Also implements the bounce-back boundary conditions for rock sites.
  !>
  !> \todo Optimize the code (see FIXME below)!
  subroutine lbe_advection_simple(N,d_adv)
    implicit none

    integer                             :: x,y,z,xa,ya,za,xb,yb,zb,bcz
    integer                             :: s,i
    type(lbe_site), dimension(0:,0:,0:) :: N
    integer                             :: d_adv
    real*8                              :: init,finl,diff
    integer                             :: xx,yy,zz

    integer                             :: xp,yp,zp,xq,yq,zq,sp,sq
    integer                             :: fooflag=0
#ifdef LOCALBC
	type(bc_check),dimension(:,:,:),allocatable :: useddirections
#endif

 ! FIXME there must be some way of improving this: preferably using
 ! one of the more terse array-transfer expressions provided by F90..

    if (.not. allocated(Nbuf)) allocate(Nbuf(0:nx+1,0:ny+1,-1:1))
#ifdef LOCALBC
	if (.not. allocated(useddirections)) allocate(useddirections(1:nx,1:ny,1:nz))
#endif
#ifdef LOCALBC
	do i=1,19
	useddirections(:,:,:)%n(i)=0
	end do
#endif

!#ifdef INTERPOLATEDBB
!  do s=1,19
!     write(*, "('s = ',I3,' n_r = ',F16.10)") s, N(6,10,6)%n_r(s)
!  enddo
!#endif

! Initialize buffer, including halo
    Nbuf(0:nx+1,0:ny+1,-1)=N(0:nx+1,0:ny+1,0)
    Nbuf(0:nx+1,0:ny+1,0)=N(0:nx+1,0:ny+1,1)

    do z=1,nz
       Nbuf(0:nx+1,0:ny+1,1)=N(0:nx+1,0:ny+1,z+1)
       do y=1,ny
          do x=1,nx
#ifndef INTERPOLATEDBB
if (bcsel.ge.0.d0) then
#ifdef SINGLEFLUID
#ifndef LOCALBC
             call lbe_boundary_condition(N,d_adv,x,y,z,Nbuf)
#endif
#ifdef LOCALBC
	     call lbe_boundary_condition(N,d_adv,x,y,z,Nbuf,useddirections)
#endif
#endif
end if
#endif !interpolatedbb
          ! Loop over each site in N.
          ! For each site, loop over each non-rest vector, and
          ! move the particles which would travel *into* that site,
          ! via that vector.


             do s=1,nnonrest
                xb = x - cx(s)
                yb = y - cy(s)
                zb = z - cz(s)
                bcz = - cz(s)

                        if (bcsel.ne.2)then
                                if (N(x,y,z)%rock_state == 0) then
                                        if(N(xb,yb,zb)%rock_state ==0) then
                                                N(x,y,z)%n_r(s)=Nbuf(xb,yb,bcz)%n_r(s)
                                        endif
                                endif
                        else
                                N(x,y,z)%n_r(s)=Nbuf(xb,yb,bcz)%n_r(s)
                        endif
             end do  !s
       end do     !x
    end do        !y
    Nbuf(0:nx+1,0:ny+1,-1)=Nbuf(0:nx+1,0:ny+1,0)
    Nbuf(0:nx+1,0:ny+1,0)=Nbuf(0:nx+1,0:ny+1,1)
 end do           !z

#ifdef LOCALBC
	deallocate(useddirections)
#endif

!#ifdef INTERPOLATEDBB
! if (.not. allocated(Nbuf)) allocate(Nbuf(0:nx+1,0:ny+1,-1:1))
!! Initialize buffer, including halo
!    Nbuf(0:nx+1,0:ny+1,-1)=N(0:nx+1,0:ny+1,0)
!    Nbuf(0:nx+1,0:ny+1,0)=N(0:nx+1,0:ny+1,1)
! do z=1,nz
!    Nbuf(0:nx+1,0:ny+1,1)=N(0:nx+1,0:ny+1,z+1)
!    do y=1,ny
!       do x=1,nx
!!         if (bcsel.ge.0.d0) then
!#ifdef SINGLEFLUID
!#ifndef LOCALBC
!             call lbe_boundary_condition(N,d_adv,x,y,z,Nbuf)
!#endif
!#endif
!!         end if
!       end do     !x
!    end do        !y
!    Nbuf(0:nx+1,0:ny+1,-1)=Nbuf(0:nx+1,0:ny+1,0)
!    Nbuf(0:nx+1,0:ny+1,0)=Nbuf(0:nx+1,0:ny+1,1)
! end do           !z
!#endif !interpolatedbb

end subroutine lbe_advection_simple



#ifdef INTERPOLATEDBB
! compute missing reflected distributions by interpolation
! on fluid node
subroutine lbe_interpolated_dist(N)
!  type(lbe_site),dimension(0:,0:,0:) :: N
  type(lbe_site),intent(inout) :: &
       &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
  integer          :: d_adv = 0
  integer          :: x,y,z,s,tz,tx,ty,he,i,endindex
  real(kind=rk)    :: fb(3)
  integer          :: xb,yb,zb
  logical          :: fexist
  character(len=1024) cfg_file_name
  integer      :: cfg_file_unit=12
  fb = 0.0_rk
  cfg_file_name = 'force.asc'
  he = halo_extent
  if (.not. allocated(lbe_Nbuf)) allocate(lbe_Nbuf(1-he:nx+he,1-he:ny+he,-he:he))
!  if (.not. allocated(lbe_Nbuf)) allocate(lbe_Nbuf(0:nx+1,0:ny+1,-1:1))
!  if (.not. allocated(Nbuf)) allocate(Nbuf(0:nx+1,0:ny+1,-1:1))
 endindex = 2*he-1 
        DEBUG_MSG("initialize halo")
! Initialize buffer, including halo
 do i=0,endindex
    lbe_Nbuf(1-he:nx+he,1-he:ny+he,-he+i)=N(1-he:nx+he,1-he:ny+he,-he+i+1)
 end do
!    lbe_Nbuf(0:nx+1,0:ny+1,-1)=N(0:nx+1,0:ny+1,0)
!    lbe_Nbuf(0:nx+1,0:ny+1,0)=N(0:nx+1,0:ny+1,1)
!    Nbuf(0:nx+1,0:ny+1,-1)=N(0:nx+1,0:ny+1,0)
!    Nbuf(0:nx+1,0:ny+1,0)=N(0:nx+1,0:ny+1,1)
        DEBUG_MSG("after initialize halo")
 do z=1,nz
! do z=1,nz
!    Nbuf(0:nx+1,0:ny+1,1)=N(0:nx+1,0:ny+1,z+1)
!        DEBUG_MSG("before passing")
!    lbe_Nbuf(0:nx+1,0:ny+1,1)=N(0:nx+1,0:ny+1,z+1)
    lbe_Nbuf(1-he:nx+he,1-he:ny+he,he)=N(1-he:nx+he,1-he:ny+he,z+he)
    do y=1,ny
       do x=1,nx
!         if (bcsel.ge.0.d0) then
#ifdef SINGLEFLUID
#ifndef LOCALBC
             call lbe_boundary_condition(N,d_adv,x,y,z,lbe_Nbuf)
!             call lbe_boundary_condition(N,d_adv,x,y,z,Nbuf)
#endif
#endif
!         end if
       end do     !x
    end do        !y
    do i=0,endindex
       lbe_Nbuf(1-he:nx+he,1-he:ny+he,-he+i)=lbe_Nbuf(1-he:nx+he,1-he:ny+he,-he+i+1)
    end do
!    lbe_Nbuf(0:nx+1,0:ny+1,-1)=lbe_Nbuf(0:nx+1,0:ny+1,0)
!    lbe_Nbuf(0:nx+1,0:ny+1,0)=lbe_Nbuf(0:nx+1,0:ny+1,1)
!    Nbuf(0:nx+1,0:ny+1,-1)=Nbuf(0:nx+1,0:ny+1,0)
!    Nbuf(0:nx+1,0:ny+1,0)=Nbuf(0:nx+1,0:ny+1,1)
 end do           !z

!    ! Total force test on single particle (all solid-fluid links basically):
!    ! works only for serial version
!    if(mod(nt,10).eq.0) then
!    do z=1,nz
!      do y=1,ny
!        do x=1,nx
!          if(N(x,y,z)%rock_state.eq.0) then
!            do s=1,nnonrest
!              xb = x + cx(s)
!              yb = y + cy(s)
!              zb = z + cz(s)
!              if(N(xb,yb,zb)%rock_state.ne.0) then
!                fb(:) = fb(:)+c(s,:)*(N(x,y,z)%n_r(bounce(s))+N(xb,yb,zb)%n_r(s))*g(s)
!!                counter=counter+1
!              endif
!            end do
!          end if
!        end do
!      end do
!    end do
!                 inquire(file=cfg_file_name,exist=fexist)
!                 if (fexist) then
!                    open(unit=cfg_file_unit,file=cfg_file_name,status='OLD'&
!                         &,position='APPEND',recl=650)
!                 else
!                    open(unit=cfg_file_unit,file=cfg_file_name,status='NEW'&
!                         &,position='APPEND',recl=650)
!                 endif
!!    write(msgstr&
!!    write(*&
!         !&,"('Force on sphere before collision : ' 3F16.10,', counter ='I0)") fb, counter
!    write(unit=cfg_file_unit&
!         &,fmt='(I10, 3F16.10)') nt,fb
!    close(cfg_file_unit)
!    end if
end subroutine lbe_interpolated_dist
#endif !interpolatedbb




!> Takes a subdomain array \c N and collides and redistributes the particles at each site.
!>
!> If this is ever run in a situation where a Fortran 95 compiler
!> is available, it could very possibly benefit from being rephrased
!> in terms of the \c forall statement in F95.

subroutine lbe_interactions_simple(N)
  type(lbe_site),dimension(0:,0:,0:) :: N
  integer                            :: x,y,z,s,tz,tx,ty
  logical                            :: minz,maxz
  real*8, dimension(19)              :: Ntemp

  ! Velocities
  real*8, dimension(3)               :: p_r, u_tilde

  ! Number density
  real*8                             :: rN_s

  ! Total fluid density
  real*8 :: rho_tot

  ! Acceleration due to force acting on total fluid
  real*8, dimension(3) :: acc_tot

  ! Boltzmann distributions
  real*8, dimension(nvecs)           :: boltz_r

  ! Equilibrium distributions at a site
  real*8, dimension(3)               :: rN_eq(nvecs)
  real*8                             :: pon_r
  real*8, dimension(3)               :: du_r
  real*8, dimension(3)               :: dd_r

!!!!!!!!!!!!!!!!MRT VARIABLES START
    real*8, dimension(19) :: S_hat_r

  ! temporary variable for force calculation
  real(kind=rk) :: g_accn_local,g_accn_local_x,g_accn_local_y
  real(kind=rk) :: forceterm
  real(kind=rk), dimension(3)        :: force_per_cell

  integer i,ii

#ifdef VELCONTROL
  real(kind=rk), dimension(3)        :: u_mean  = 0._rk
  real(kind=rk), dimension(3)        :: sum_vel = 0._rk
  real(kind=rk), dimension(3)        :: u_ref
  integer                            :: ierror, counter = 0
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Collision stage
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call start_timer(ti_intc)

  call lbe_count_sites(N)

#ifdef MD
  force_per_cell = totalforce/n_sites_fluid
#endif

#ifdef VELCONTROL
  ! Loop over all lattice sites.
  do z = 1, nz
    do y = 1, ny
      do x = 1, nx
        if (N(x,y,z)%rock_state == 0.d0) then
          rN_s = 0.d0
          p_r = 0.d0

          do s = 1, nvecs
            rN_s = rN_s + N(x, y, z)%n_r(s) * g(s)
            p_r = p_r + N(x, y, z)%n_r(s) * g(s) * c(s, :)
          end do

          ! Calculate weighted velocity.
          ! Use a cutoff to avoid division by very small numbers.
          u_tilde = p_r / max(rN_s, 10.d-9)
          sum_vel = sum_vel + u_tilde
        end if
      end do
    end do
  end do

  !write(*, "('sum_vel = ',3F16.10)") sum_vel
  call mpi_allreduce(sum_vel,u_mean,3,LBE_REAL,MPI_SUM,MPI_COMM_WORLD,ierror)  
  u_mean = u_mean/real(n_sites, kind=rk) !average for whole domain
  !u_mean = sum_vel/real(nx*ny*nz, kind=rk)
  u_ref = (/u_refx, u_refy, u_refz/)
  err_0 = u_ref - u_mean
  write(msgstr, "('Error_0 = ',3F16.10)") err_0/u_ref*100._rk
  call log_msg(msgstr)
  !write(*, "('u_mean = ',3F16.10)") u_mean
  !write(*, "('u_ref = ',3F16.10)") u_ref 
  !write(*, "('K_p, K_i, K_d = ',3F16.10)") k_p, k_i, k_d
  !write(*, "('u_mean = ',3F16.10)") u_mean
  !write(*, "('n_sites = ',I10)") n_sites
  !write(*, "('Error = ',3F16.10)") err_0
  sum_vel = 0.0_rk
  !gravity_new = gravity_old + (k_p+k_i+k_d)*err_0 - (k_p+2._rk*k_d)*err_1+&
  !              &k_d*err_2
  gravity_new = gravity_old + err_0/k_p**2.0_rk
  gravity_old = gravity_new
  err_2 = err_1
  err_1 = err_0
  write(msgstr, "('Gravity = ',3F16.10)") gravity_new
  call log_msg(msgstr)

#endif  !VELCONTROL


!write(*,"('Inside interactions_simple')")
!write(*,"('After interpolated BC')")


  ! Loop over all lattice sites.
  do z = 1, nz
    do y = 1, ny
      do x = 1, nx
        ! Collision is only performed on fluid sites.
        if (N(x,y,z)%rock_state == 0.d0) then
          ! ================================
          ! COMPUTE DENSITIES AND VELOCITIES
          ! ================================

          ! Compute species densities and velocities:
          ! - rN_s = particle density of red species
          ! - p_r = velocity of red species
          ! TODO
          ! - Wouldn't it be simpler to use Fortran's sum function here?
          rN_s = 0.d0
          p_r = 0.d0

          do s = 1, nvecs
            rN_s = rN_s + N(x, y, z)%n_r(s) * g(s)
            p_r = p_r + N(x, y, z)%n_r(s) * g(s) * c(s, :)
          end do

          ! Calculate weighted velocity.
          ! Use a cutoff to avoid division by very small numbers.
          u_tilde = p_r / max(rN_s, 10.d-9)

          ! Calculate total fluid density.
          rho_tot = amass_r * rN_s

          ! Calculate inverse particle density.
          if(psifunc .eq. 0) then
            ! This expression is only for the original funny clipped version of the code.
            ! Evaluates to 1 unless densities go over 1.
            pon_r = min(1.d0 / max(rN_s, 10.d-9), 1.d0)
          else
            ! Sane code that does what the literature says.
            pon_r = 1.d0 / max(rN_s, 10.d-9)
          end if

          ! ==========================
          ! GRAVITATIONAL ACCELERATION
          ! ==========================

          ! Add the force due to g_accn to the 0-component of lbe_force.
          ! The force equals the acceleration times the total fluid density.
          ! NOTES:
          ! - As in the previous revision of the code, g_accn etc. act only on fluid sites, not on rock sites.
          ! TODO:
          ! - This functionality should be moved to some forcing module.

          ! Compute global position of the local site.
          tx = x + ccoords(1) * nx
          ty = y + ccoords(2) * ny
          tz = z + ccoords(3) * nz

#ifndef MRT_SINGLEPARTICLE
#ifdef MD
          if(forcebalance) then
            lbe_force(1, 0, x, y, z) = lbe_force(1, 0, x, y, z) + force_per_cell(1) 
            lbe_force(2, 0, x, y, z) = lbe_force(2, 0, x, y, z) + force_per_cell(2)
            lbe_force(3, 0, x, y, z) = lbe_force(3, 0, x, y, z) + force_per_cell(3)
          end if
#endif
          
#ifdef VELCONTROL 
          g_accn_x = gravity_new(1)
          g_accn_y = gravity_new(2)
          g_accn   = gravity_new(3)
#endif          
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
            minz = (start(3)==1)
            maxz = (start(3)>=(tnz-nz))

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
            lbe_force(1, 0, x, y, z) = lbe_force(1, 0, x, y, z) + g_accn_local_x * rho_tot
            lbe_force(2, 0, x, y, z) = lbe_force(2, 0, x, y, z) + g_accn_local_y * rho_tot
            lbe_force(3, 0, x, y, z) = lbe_force(3, 0, x, y, z) + g_accn_local * rho_tot
          endif
#endif          
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

#ifndef MRT_SINGLEPARTICLE
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
          du_r(:) = get_tau_r(N, x, y, z) * (pon_r * lbe_force(:, 1, x, y, z) / amass_r + acc_tot)
#endif          

          !du_r(:) = 0.0_rk
!          du_r(:) = (/0.0_rk, 0.0_rk, 0.0_rk/)

!          do ii=1,3
!!            if(du_r(ii).ne.0.0_rk) then
!              write(*, "('pos = ',3I3,' du_r = ',3F16.10)") x,y,z,du_r
!!            end if
!          end do

          ! The dd_r array (Phys. Rev. E 65,(2002): Guo) used in the previous revision will be deprecated soon
          ! since all forces will be coupled in the same way in the future (user's choice: either Guo or Shan-Chen).
          ! For a temporary measure, dd_r is set to zero.
          dd_r(:) = 0.0_rk !0.d0

!#ifndef INTERPOLATEDBB
!          ! Perform collision.
!          ! TODO
!          ! - redefine boltz_dist by taking out dd_r later.
!          call boltz_dist(u_tilde(1), u_tilde(2), u_tilde(3), du_r(1), du_r(2), du_r(3), dd_r(1), dd_r(2), dd_r(3), boltz_r)
!          rN_eq = rN_s * boltz_r
!#else
          ! Required for interpolated bounce back scheme to prevent mass leakage
          ! in the system. Refer to He and Luo (1997) for this scheme
          ! Note: densities are already multiplied inside the function
!          call boltz_dist_const_density(u_tilde+du_r, rN_s, N(x,y,z)%rho_0, rN_eq)
!          call boltz_dist_const_density(u_tilde+du_r, rN_s, rN_s, rN_eq)
!#endif
!          N(x, y, z)%n_r = N(x, y, z)%n_r - get_omega_r(N, x, y, z) * (N(x, y, z)%n_r - rN_eq)
    if ( collisiontype_id == MRT_id) then
      S_hat_r(1) = 0.d0  ! rho           0.00d0
      S_hat_r(2) = s02_r !omegabulk_r    ! e       1.19d0
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
!    endif
!
!
!!   write(*,"(' pos =',3I3,' Rho = ',F16.10)") x,y,z,rN_s
!
!!         if(x.eq.6.and.y.eq.6.and.z.eq.6) then
!!           write(*,"('pos = ',3I3,' dist  = ',19F16.10)") x,y,z,N(x,y,z)%n_r(:)
!!         end if
!
!
!    if ( collisiontype_id == MRT_id ) then
      call mrt_dist(u_tilde + du_r - dd_r, dd_r / (get_tau_r(N, x, y, z) - 0.5d0), N(x, y, z)%n_r, S_hat_r)
    else if ( collisiontype_id == BGK_id) then
      !call boltz_dist_constrho(u_tilde(1), u_tilde(2), u_tilde(3), du_r(1), du_r(2), du_r(3), dd_r(1), dd_r(2), dd_r(3), rN_s, rN_eq)
      !call boltz_dist(u_tilde(1), u_tilde(2), u_tilde(3), 0.0_rk,0.0_rk,0.0_rk, 0.0_rk,0.0_rk,0.0_rk, boltz_r)
      call boltz_dist(u_tilde(1), u_tilde(2), u_tilde(3), du_r(1), du_r(2), du_r(3), dd_r(1), dd_r(2), dd_r(3), boltz_r)
      rN_eq = rN_s * boltz_r
      N(x, y, z)%n_r = N(x, y, z)%n_r - get_omega_r(N, x, y, z) * (N(x, y, z)%n_r - rN_eq)
    end if

#ifndef INTERPOLATEDBB
        ! End of code which is executed on fluid sites.
        else

          ! ==========================
          ! BOUNCE-BACK FOR ROCK SITES
          ! ==========================

          ! Rocks stand still.
          if(bcsel == 2) then
            Ntemp = N(x, y, z)%n_r

            do s = 1, 19
              N(x, y, z)%n_r(bounce(s)) = Ntemp(s)
            enddo
          endif
#endif !interpolatedbb
        endif
      end do
    end do
  end do
  ! Finished running over the lattice.
  call stop_timer(ti_intc)
end subroutine lbe_interactions_simple

end module lbe_collision_simple_module
