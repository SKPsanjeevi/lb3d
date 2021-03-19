#ifdef AXISYM 
#include "lbe.h"

!> code related to the axiysmmetric + multicomponent part.
!> For details of the model see http://link.aps.org/doi/10.1103/PhysRevE.88.013309
!> y is used for radial coordinares, and z is used for the axial coordinates

module lbe_axisym_module
  use lbe_globals_module
  use lbe_helper_module, only: is_wall, is_fluid, is_colloid, CS_EOS
  use lbe_log_module
  use lbe_parms_module, only: BGK_id,collisiontype_id,fluctuating_MRT_id&
       &,MRT_id,SCMP,ZFOSdiag,acccoef,amass_b,amass_r,amass_s&
       &,bcsel,beta,boundary_cond,d_0,g_accn_fluid_r,g_accn_fluid_b,g_accn_fluid_s&
       &,g_accn,g_accn_x,g_accn_y,g_accn_max&
       &,g_accn_max_x,g_accn_max_y,g_accn_min,g_accn_min_x,g_accn_min_y,g_br&
       &,g_bs,g_rr,g_ss,g_wb,g_wr,kbT,Tcs,omega_d&
       &,omegabulk_b,omegabulk_r,omegabulk_s,psifunc,s03_b,s03_r,s03_s,s05_b&
       &,s05_r,s05_s,s11_b,s11_r,s11_s,s17_b,s17_r,s17_s&
       &,tau_wb,tau_wr,inv_fluid,nx,ny,nz, &
       & omega_r, omega_b, omega_s, &
       & get_tau_r, get_omega_r, rock_colour_double, g_value_double
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
#if !defined(SINGLEFLUID) && defined(NOSURFACTANT) && defined(IBM_PART) && defined(IBM_INDEXFIELD)
  use lsuperstruct_data_module, only: g_ri, g_bi, ibm_colour
  use lsuperstruct_interface_module
  use lbe_force_interface_module, only : add_force_to_n_r, add_force_to_n_b, add_force_to_all
#endif
  implicit none
  private

  real(kind = rk), dimension(9) :: cz2d = (/ 1, 0, -1, 0, 1, -1, -1,  1, 0 /)
  real(kind = rk), dimension(9) :: cy2d = (/ 0, 1, 0, -1, 1, 1, -1,  -1, 0 /)
  real(kind = rk), parameter :: weight(nvecs) = &
    (/ 1.0d0/18.0d0,1.0d0/18.0d0, 1.0d0/18.0d0, 1.0d0/18.0d0, 1.0d0/18.0d0, 1.0d0/18.0d0, &
       1.0d0/36.0d0,1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, &
       1.0d0/36.0d0,1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, &
       1.0d0/3.0d0 &
    /)
  real(kind = rk) :: cs2 = 1.0d0/3.0d0
  real(kind = rk) :: cs4 = 1.0d0/9.0d0


  public lbe_axisym_axis_bc,lbe_axisym_correction_cty,lbe_axisym_correction_mom, &
       & lbe_axisym_correction_sc_forces, lbe_gradient1, lbe_gradient2, lbe_get_psi, &
       & lbe_get_hydrovar
  
  type(lbe_site), dimension(:,:,:), allocatable :: Nbuf
  
contains

  subroutine lbe_get_hydrovar(N,u,rho)
    implicit none
    integer :: x, y, z
    type(lbe_site), dimension(0:,0:,0:),intent(in) :: N

    real*8, dimension(0:,0:,0:,0:)   ,intent(out) :: rho
    real*8, dimension(0:,0:,0:,0:,0:),intent(out) :: u
                     !(x/y/z, x,y,z, r/b/s)!  
      do z=1,nz
         do y=1,ny
            do x=1,nx

               rho(x,y,z,1) = sum(N(x,y,z)%n_r(:)*g(:))
               u(1,x,y,z,1) = sum(N(x,y,z)%n_r(:)*g(:)*cx(:))
               u(2,x,y,z,1) = sum(N(x,y,z)%n_r(:)*g(:)*cy(:))
               u(3,x,y,z,1) = sum(N(x,y,z)%n_r(:)*g(:)*cz(:))
#ifndef SINGLEFLUID
               rho(x,y,z,2) = sum(N(x,y,z)%n_b(:)*g(:))
               u(1,x,y,z,2) = sum(N(x,y,z)%n_b(:)*g(:)*cx(:))
               u(2,x,y,z,2) = sum(N(x,y,z)%n_b(:)*g(:)*cy(:))
               u(3,x,y,z,2) = sum(N(x,y,z)%n_b(:)*g(:)*cz(:))
#ifndef NOSURFACTANT
               rho(x,y,z,3) = sum(N(x,y,z)%n_s(:)*g(:))
               u(1,x,y,z,3) = sum(N(x,y,z)%n_s(:)*g(:)*cx(:))
               u(2,x,y,z,3) = sum(N(x,y,z)%n_s(:)*g(:)*cy(:))
               u(3,x,y,z,3) = sum(N(x,y,z)%n_s(:)*g(:)*cz(:))
#endif
#endif
               end do
            end do
         end do
    
  end subroutine lbe_get_hydrovar


  subroutine lbe_get_psi(N,psi_mat)

    implicit none
    integer :: x, y, z
    type(lbe_site), dimension(0:,0:,0:),intent(inout) :: N

    real*8 :: rho_r
    real*8, dimension(1:,1:,1:,1:),intent(out) :: psi_mat

#ifndef SINGLEFLUID
    real*8 :: rho_b
#endif
#ifndef NOSURFACTANT
    real*8 :: rho_s
#endif


    do z=1,nz
       do y=1,ny
          do x=2,nx-1

             rho_r = sum(N(x,y,z)%n_r(:)*g(:))
#ifndef SINGLEFLUID
             rho_b = sum(N(x,y,z)%n_b(:)*g(:))
#ifndef NOSURFACTANT
             rho_s = sum(N(x,y,z)%n_s(:)*g(:))
#endif
#endif
             select case(psifunc)
             case (0)
                psi_mat(x,y,z,1) = min(1.d0, dble(rho_r))
#ifndef SINGLEFLUID
                psi_mat(x,y,z,2) = min(1.d0, dble(rho_b))
#ifndef NOSURFACTANT
                psi_mat(x,y,z,3) = min(1.d0, dble(rho_s))
#endif
#endif
             case (1)
                psi_mat(x,y,z,1) = rho_r
#ifndef SINGLEFLUID
                psi_mat(x,y,z,2) = rho_b
#ifndef NOSURFACTANT
                psi_mat(x,y,z,3) = rho_s
#endif
#endif
             case (2)
                ! psi = 1 - exp(-n)
                psi_mat(x,y,z,1) = 1.d0 - exp(-rho_r)
#ifndef SINGLEFLUID
                psi_mat(x,y,z,2) = 1.d0 - exp(-rho_b)
#ifndef NOSURFACTANT
                psi_mat(x,y,z,3) = 1.d0 - exp(-rho_s)
#endif
#endif
             case default
                call log_msg("ERROR: Unknown psi functional, aborting ...")
                call Abend
             end select

!             psi_mat(x,y,z,1) = (y-1.5)*(y-1.5)*(y-1.5)*(y-1.5)
          end do
       end do
    end do


    psi_mat(1,:,:,:) = psi_mat(nx-1,:,:,:)
    psi_mat(nx,:,:,:) = psi_mat(2,:,:,:)

  end subroutine lbe_get_psi


#ifdef SINGLEFLUID
    subroutine lbe_gradient1(psi_mat, dy_psi_r)
#else
#ifdef NOSURFACTANT
    subroutine lbe_gradient1(psi_mat, dy_psi_r, dy_psi_b)
#else
    subroutine lbe_gradient1(psi_mat, dy_psi_r, dy_psi_b, dy_psi_s)
#endif
#endif
      implicit none
      integer :: x,y,z
      integer :: xa,ya,za,x2a,y2a,z2a
      integer :: s

     
      real*8, dimension(1:,1:,1:,1:), intent(in) :: psi_mat

#ifdef SINGLEFLUID      
      real*8, dimension(1:,1:,1:),intent(out) :: dy_psi_r
#else
#ifdef NOSURFACTANT
      real*8, dimension(1:,1:,1:),intent(out) :: dy_psi_r
      real*8, dimension(1:,1:,1:),intent(out) :: dy_psi_b
#else
      real*8, dimension(1:,1:,1:),intent(out) :: dy_psi_r
      real*8, dimension(1:,1:,1:),intent(out) :: dy_psi_b
      real*8, dimension(1:,1:,1:),intent(out) :: dy_psi_s
#endif
#endif

#ifdef SINGLEFLUID      
      real*8, dimension(:,:,:), allocatable  :: psi_r
#else
#ifdef NOSURFACTANT
      real*8, dimension(:,:,:), allocatable  :: psi_r
      real*8, dimension(:,:,:), allocatable  :: psi_b
#else
      real*8, dimension(:,:,:), allocatable  :: psi_r
      real*8, dimension(:,:,:), allocatable  :: psi_b
      real*8, dimension(:,:,:), allocatable  :: psi_s
#endif
#endif



#ifdef SINGLEFLUID      
      if (.not. allocated(psi_r)) allocate(psi_r(1:nx,-1:ny+2,-1:nz+2))
      psi_r(1:nx,1:ny,1:nz) = psi_mat(:,:,:,1)
#else
#ifdef NOSURFACTANT
      if (.not. allocated(psi_r)) allocate(psi_r(1:nx,-1:ny+2,-1:nz+2))
      if (.not. allocated(psi_b)) allocate(psi_b(1:nx,-1:ny+2,-1:nz+2))
      psi_r(1:nx,1:ny,1:nz) = psi_mat(:,:,:,1)
      psi_b(1:nx,1:ny,1:nz) = psi_mat(:,:,:,2)
#else
      if (.not. allocated(psi_r)) allocate(psi_r(1:nx,-1:ny+2,-1:nz+2))
      if (.not. allocated(psi_b)) allocate(psi_b(1:nx,-1:ny+2,-1:nz+2))
      if (.not. allocated(psi_s)) allocate(psi_s(1:nx,-1:ny+2,-1:nz+2))
      psi_r(1:nx,1:ny,1:nz) = psi_mat(:,:,:,1)
      psi_b(1:nx,1:ny,1:nz) = psi_mat(:,:,:,2)
      psi_s(1:nx,1:ny,1:nz) = psi_mat(:,:,:,3)
#endif
#endif

! axis is at y = 1.5 [1-2]
#ifdef SINGLEFLUID      
      psi_r(:,1,:) = psi_r(:,2,:)
      psi_r(:,0,:) = psi_r(:,3,:)
      psi_r(:,-1,:) = psi_r(:,4,:)

      psi_r(:,ny,:) = psi_r(:,ny-1,:)
      psi_r(:,ny+1,:) = psi_r(:,ny-1,:)
      psi_r(:,ny+2,:) = psi_r(:,ny-1,:)

      !z-axis pbc
      psi_r(:,:,-1) = psi_r(:,:,nz-3)
      psi_r(:,:,0) = psi_r(:,:,nz-2)
      psi_r(:,:,1) = psi_r(:,:,nz-1)
      psi_r(:,:,nz) = psi_r(:,:,2)
      psi_r(:,:,nz+1) = psi_r(:,:,3)
      psi_r(:,:,nz+2) = psi_r(:,:,4)
     
!!$      call xcomm(psi_r)
!!$      call zcomm(psi_r)
#else
#ifdef NOSURFACTANT
      psi_r(:,1,:) = psi_r(:,2,:)
      psi_r(:,0,:) = psi_r(:,3,:)
      psi_r(:,-1,:) = psi_r(:,4,:)

      psi_r(:,ny+1,:) = psi_r(:,ny,:)
      psi_r(:,ny+2,:) = psi_r(:,ny,:)

      psi_b(:,1,:) = psi_b(:,2,:)
      psi_b(:,0,:) = psi_b(:,3,:)
      psi_b(:,-1,:) = psi_b(:,4,:)

      psi_b(:,ny+1,:) = psi_b(:,ny,:)
      psi_b(:,ny+2,:) = psi_b(:,ny,:)


      !z-axis pbc

      psi_r(:,:,-1) = psi_r(:,:,nz-3)
      psi_r(:,:,0) = psi_r(:,:,nz-2)
      psi_r(:,:,1) = psi_r(:,:,nz-1)
      psi_r(:,:,nz) = psi_r(:,:,2)
      psi_r(:,:,nz+1) = psi_r(:,:,3)
      psi_r(:,:,nz+2) = psi_r(:,:,4)

      psi_b(:,:,-1) = psi_b(:,:,nz-3)
      psi_b(:,:,0) = psi_b(:,:,nz-2)
      psi_b(:,:,1) = psi_b(:,:,nz-1)
      psi_b(:,:,nz) = psi_b(:,:,2)
      psi_b(:,:,nz+1) = psi_b(:,:,3)
      psi_b(:,:,nz+2) = psi_b(:,:,4)

#else
      ! not implemented

#endif
#endif

      do x = 2, nx-1
         do y = 1, ny
            do z = 1, nz
               
#ifdef SINGLEFLUID      
               dy_psi_r(x, y, z) = 0.0_rk
#else
#ifdef NOSURFACTANT
               dy_psi_r(x, y, z) = 0.0_rk
               dy_psi_b(x, y, z) = 0.0_rk
#else
               dy_psi_r(x, y, z) = 0.0_rk
               dy_psi_b(x, y, z) = 0.0_rk
               dy_psi_s(x, y, z) = 0.0_rk
#endif
#endif
               
               do s = 1, 8               
                  xa = x
                  ya = y + cy2d(s)
                  za = z + cz2d(s)
                  
                  x2a = x
                  y2a = y + 2.0*cy2d(s)
                  z2a = z + 2.0*cz2d(s)
#ifdef SINGLEFLUID      
                  dy_psi_r(x, y, z) = dy_psi_r(x, y, z) + (1.0/36.0)*(8.0*psi_r(xa,ya,za) &
                       &- psi_r(x2a,y2a,z2a))*cy2d(s)
#else
#ifdef NOSURFACTANT
                  dy_psi_r(x, y, z) = dy_psi_r(x, y, z) + (1.0/36.0)*(8.0*psi_r(xa,ya,za) &
                       &- psi_r(x2a,y2a,z2a))*cy2d(s)
                  dy_psi_b(x, y, z) = dy_psi_b(x, y, z) + (1.0/36.0)*(8.0*psi_b(xa,ya,za) &
                       &- psi_b(x2a,y2a,z2a))*cy2d(s)

#else
                  dy_psi_r(x, y, z) = dy_psi_r(x, y, z) + (1.0/36.0)*(8.0*psi_r(xa,ya,za) &
                       &- psi_r(x2a,y2a,z2a))*cy2d(s)
                  dy_psi_b(x, y, z) = dy_psi_b(x, y, z) + (1.0/36.0)*(8.0*psi_b(xa,ya,za) &
                       &- psi_b(x2a,y2a,z2a))*cy2d(s)
                  dy_psi_s(x, y, z) = dy_psi_s(x, y, z) + (1.0/36.0)*(8.0*psi_s(xa,ya,za) &
                       &- psi_s(x2a,y2a,z2a))*cy2d(s)
#endif
#endif

               end do
!               print *,x,y,z,4*(y-1.5)*(y-1.5)*(y-1.5),dy_psi_r(x,y,z)
            end do ! z
!!$            if (y .gt. 4) then
!!$               stop
!!$            end if
         end do ! y
      end do ! x

#ifdef SINGLEFLUID      
      dy_psi_r(1,:,:) = dy_psi_r(3,:,:)
      dy_psi_r(4,:,:) = dy_psi_r(2,:,:)
#else
#ifdef NOSURFACTANT
      dy_psi_r(1,:,:) = dy_psi_r(3,:,:)
      dy_psi_r(4,:,:) = dy_psi_r(2,:,:)
      dy_psi_b(1,:,:) = dy_psi_b(3,:,:)
      dy_psi_b(4,:,:) = dy_psi_b(2,:,:)
#else
      dy_psi_r(1,:,:) = dy_psi_r(3,:,:)
      dy_psi_r(4,:,:) = dy_psi_r(2,:,:)
      dy_psi_b(1,:,:) = dy_psi_b(3,:,:)
      dy_psi_b(4,:,:) = dy_psi_b(2,:,:)
      dy_psi_s(1,:,:) = dy_psi_s(3,:,:)
      dy_psi_s(4,:,:) = dy_psi_s(2,:,:)
#endif
#endif

#ifdef SINGLEFLUID      
      deallocate(psi_r)
#else
#ifdef NOSURFACTANT
      deallocate(psi_r)
      deallocate(psi_b)
#else
      deallocate(psi_r)
      dealloca te(psi_b)
      deallocate(psi_s)
#endif
#endif
 
    end subroutine lbe_gradient1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef SINGLEFLUID
    subroutine lbe_gradient2(dy_psi_r,dyz_psi_r,dyy_psi_r)
#else
#ifdef NOSURFACTANT
    subroutine lbe_gradient2(dy_psi_r,dyz_psi_r,dyy_psi_r,&
         & dy_psi_b,dyz_psi_b,dyy_psi_b)
#else
    subroutine lbe_gradient2(dy_psi_r,dyz_psi_r,dyy_psi_r,&
         & dy_psi_b,dyz_psi_b,dyy_psi_b,&
         & dy_psi_s,dyz_psi_s,dyy_psi_s)
#endif
#endif
      implicit none
      integer :: x,y,z
      integer :: xa,ya,za,x2a,y2a,z2a
      integer :: s


#ifdef SINGLEFLUID      
      real*8, dimension(1:,1:,1:),intent(in)  :: dy_psi_r
      real*8, dimension(1:,1:,1:),intent(out) :: dyz_psi_r,dyy_psi_r
#else
#ifdef NOSURFACTANT
      real*8, dimension(1:,1:,1:),intent(in)  :: dy_psi_r
      real*8, dimension(1:,1:,1:),intent(in)  :: dy_psi_b
      real*8, dimension(1:,1:,1:),intent(out) :: dyz_psi_r,dyy_psi_r
      real*8, dimension(1:,1:,1:),intent(out) :: dyz_psi_b,dyy_psi_b
#else
      real*8, dimension(1:,1:,1:),intent(in)  :: dy_psi_r
      real*8, dimension(1:,1:,1:),intent(in)  :: dy_psi_b
      real*8, dimension(1:,1:,1:),intent(in)  :: dy_psi_s
      real*8, dimension(1:,1:,1:),intent(out) :: dyz_psi_r,dyy_psi_r
      real*8, dimension(1:,1:,1:),intent(out) :: dyz_psi_b,dyy_psi_b
      real*8, dimension(1:,1:,1:),intent(out) :: dyz_psi_s,dyy_psi_s
#endif
#endif

#ifdef SINGLEFLUID      
      real*8, dimension(:,:,:), allocatable  :: dy_psi_r_buff
#else
#ifdef NOSURFACTANT
      real*8, dimension(:,:,:), allocatable  :: dy_psi_r_buff
      real*8, dimension(:,:,:), allocatable  :: dy_psi_b_buff
#else
      real*8, dimension(:,:,:), allocatable  :: dy_psi_r_buff
      real*8, dimension(:,:,:), allocatable  :: dy_psi_b_buff
      real*8, dimension(:,:,:), allocatable  :: dy_psi_s_buff
#endif
#endif


#ifdef SINGLEFLUID      
      if (.not. allocated(dy_psi_r_buff)) allocate(dy_psi_r_buff(1:nx,-1:ny+2,-1:nz+2))
      dy_psi_r_buff(1:nx,1:ny,1:nz) = dy_psi_r(:,:,:)
#else
#ifdef NOSURFACTANT
      if (.not. allocated(dy_psi_r_buff)) allocate(dy_psi_r_buff(1:nx,-1:ny+2,-1:nz+2))
      if (.not. allocated(dy_psi_b_buff)) allocate(dy_psi_b_buff(1:nx,-1:ny+2,-1:nz+2))
      dy_psi_r_buff(1:nx,1:ny,1:nz) = dy_psi_r(:,:,:)
      dy_psi_b_buff(1:nx,1:ny,1:nz) = dy_psi_b(:,:,:)
#else
      if (.not. allocated(dy_psi_r_buff)) allocate(dy_psi_r_buff(1:nx,-1:ny+2,-1:nz+2))
      if (.not. allocated(dy_psi_b_buff)) allocate(dy_psi_b_buff(1:nx,-1:ny+2,-1:nz+2))
      if (.not. allocated(dy_psi_s_buff)) allocate(dy_psi_s_buff(1:nx,-1:ny+2,-1:nz+2))
      dy_psi_r_buff(1:nx,1:ny,1:nz) = dy_psi_r(:,:,:)
      dy_psi_b_buff(1:nx,1:ny,1:nz) = dy_psi_b(:,:,:)
      dy_psi_s_buff(1:nx,1:ny,1:nz) = dy_psi_s(:,:,:)
#endif
#endif

#ifdef SINGLEFLUID      
      dy_psi_r_buff(:,-1,:) = -dy_psi_r_buff(:,4,:)
      dy_psi_r_buff(:,0,:) = -dy_psi_r_buff(:,3,:)
      dy_psi_r_buff(:,1,:) = -dy_psi_r_buff(:,2,:)

      dy_psi_r_buff(:,ny,:)   = dy_psi_r_buff(:,ny-1,:)
      dy_psi_r_buff(:,ny+1,:) = dy_psi_r_buff(:,ny-1,:)
      dy_psi_r_buff(:,ny+2,:) = dy_psi_r_buff(:,ny-1,:)

      dy_psi_r_buff(:,:,nz)   = dy_psi_r_buff(:,:,2)
      dy_psi_r_buff(:,:,nz+1) = dy_psi_r_buff(:,:,3)
      dy_psi_r_buff(:,:,nz+2) = dy_psi_r_buff(:,:,4)
      dy_psi_r_buff(:,:,-1)   = dy_psi_r_buff(:,:,nz-3)
      dy_psi_r_buff(:,:,0) = dy_psi_r_buff(:,:,nz-2)
      dy_psi_r_buff(:,:,1) = dy_psi_r_buff(:,:,nz-1)


#else
#ifdef NOSURFACTANT
      dy_psi_r_buff(:,-1,:) = -dy_psi_r_buff(:,4,:)
      dy_psi_r_buff(:,0,:) = -dy_psi_r_buff(:,3,:)
      dy_psi_r_buff(:,1,:) = -dy_psi_r_buff(:,2,:)

      dy_psi_r_buff(:,ny,:)   = dy_psi_r_buff(:,ny-1,:)
      dy_psi_r_buff(:,ny+1,:) = dy_psi_r_buff(:,ny,:)
      dy_psi_r_buff(:,ny+2,:) = dy_psi_r_buff(:,ny,:)

      dy_psi_r_buff(:,:,nz)   = dy_psi_r_buff(:,:,2)
      dy_psi_r_buff(:,:,nz+1) = dy_psi_r_buff(:,:,3)
      dy_psi_r_buff(:,:,nz+2) = dy_psi_r_buff(:,:,4)
      dy_psi_r_buff(:,:,-1)   = dy_psi_r_buff(:,:,nz-3)
      dy_psi_r_buff(:,:,0) = dy_psi_r_buff(:,:,nz-2)
      dy_psi_r_buff(:,:,1) = dy_psi_r_buff(:,:,nz-1)


      dy_psi_b_buff(:,-1,:) = -dy_psi_b_buff(:,4,:)
      dy_psi_b_buff(:,0,:) = -dy_psi_b_buff(:,3,:)
      dy_psi_b_buff(:,1,:) = -dy_psi_b_buff(:,2,:)

      dy_psi_b_buff(:,ny,:)   = dy_psi_b_buff(:,ny-1,:)
      dy_psi_b_buff(:,ny+1,:) = dy_psi_b_buff(:,ny,:)
      dy_psi_b_buff(:,ny+2,:) = dy_psi_b_buff(:,ny,:)

      dy_psi_b_buff(:,:,nz)   = dy_psi_b_buff(:,:,2)
      dy_psi_b_buff(:,:,nz+1) = dy_psi_b_buff(:,:,3)
      dy_psi_b_buff(:,:,nz+2) = dy_psi_b_buff(:,:,4)
      dy_psi_b_buff(:,:,-1)   = dy_psi_b_buff(:,:,nz-3)
      dy_psi_b_buff(:,:,0) = dy_psi_b_buff(:,:,nz-2)
      dy_psi_b_buff(:,:,1) = dy_psi_b_buff(:,:,nz-1)


#else
      dy_psi_r_buff(:,-1,:) = -dy_psi_r_buff(:,4,:)
      dy_psi_r_buff(:,1,:) = -dy_psi_r_buff(:,2,:)
      dy_psi_r_buff(:,0,:) = -dy_psi_r_buff(:,3,:)

      dy_psi_r_buff(:,ny,:)   = dy_psi_r_buff(:,ny-1,:)
      dy_psi_r_buff(:,ny+1,:) = dy_psi_r_buff(:,ny,:)
      dy_psi_r_buff(:,ny+2,:) = dy_psi_r_buff(:,ny,:)

      dy_psi_r_buff(:,:,nz)   = dy_psi_r_buff(:,:,2)
      dy_psi_r_buff(:,:,nz+1) = dy_psi_r_buff(:,:,3)
      dy_psi_r_buff(:,:,nz+2) = dy_psi_r_buff(:,:,4)
      dy_psi_r_buff(:,:,-1)   = dy_psi_r_buff(:,:,nz-3)
      dy_psi_r_buff(:,:,0) = dy_psi_r_buff(:,:,nz-2)
      dy_psi_r_buff(:,:,1) = dy_psi_r_buff(:,:,nz-1)


      dy_psi_b_buff(:,-1,:) = -dy_psi_b_buff(:,4,:)
      dy_psi_b_buff(:,1,:) = -dy_psi_b_buff(:,2,:)
      dy_psi_b_buff(:,0,:) = -dy_psi_b_buff(:,3,:)

      dy_psi_b_buff(:,ny,:)   = dy_psi_b_buff(:,ny-1,:)
      dy_psi_b_buff(:,ny+1,:) = dy_psi_b_buff(:,ny,:)
      dy_psi_b_buff(:,ny+2,:) = dy_psi_b_buff(:,ny,:)

      dy_psi_b_buff(:,:,nz)   = dy_psi_b_buff(:,:,2)
      dy_psi_b_buff(:,:,nz+1) = dy_psi_b_buff(:,:,3)
      dy_psi_b_buff(:,:,nz+2) = dy_psi_b_buff(:,:,4)
      dy_psi_b_buff(:,:,-1)   = dy_psi_b_buff(:,:,nz-3)
      dy_psi_b_buff(:,:,0) = dy_psi_b_buff(:,:,nz-2)
      dy_psi_b_buff(:,:,1) = dy_psi_b_buff(:,:,nz-1)


      dy_psi_s_buff(:,-1,:) = -dy_psi_s_buff(:,4,:)
      dy_psi_s_buff(:,1,:) = -dy_psi_s_buff(:,2,:)
      dy_psi_s_buff(:,0,:) = -dy_psi_s_buff(:,3,:)

      dy_psi_s_buff(:,ny,:)   = dy_psi_s_buff(:,ny-1,:)
      dy_psi_s_buff(:,ny+1,:) = dy_psi_s_buff(:,ny,:)
      dy_psi_s_buff(:,ny+2,:) = dy_psi_s_buff(:,ny,:)

      dy_psi_s_buff(:,:,nz)   = dy_psi_s_buff(:,:,2)
      dy_psi_s_buff(:,:,nz+1) = dy_psi_s_buff(:,:,3)
      dy_psi_s_buff(:,:,nz+2) = dy_psi_s_buff(:,:,4)
      dy_psi_s_buff(:,:,-1)   = dy_psi_s_buff(:,:,nz-3)
      dy_psi_s_buff(:,:,0) = dy_psi_s_buff(:,:,nz-2)
      dy_psi_s_buff(:,:,1) = dy_psi_s_buff(:,:,nz-1)
#endif
#endif


      do x = 2, nx-1
         do y = 1, ny
            do z = 1, nz

#ifdef SINGLEFLUID      
               dyz_psi_r(x, y, z) = 0.0
               dyy_psi_r(x, y, z) = 0.0
#else
#ifdef NOSURFACTANT
               dyz_psi_r(x, y, z) = 0.0
               dyy_psi_r(x, y, z) = 0.0

               dyz_psi_b(x, y, z) = 0.0
               dyy_psi_b(x, y, z) = 0.0
#else
               dyz_psi_r(x, y, z) = 0.0
               dyy_psi_r(x, y, z) = 0.0

               dyz_psi_b(x, y, z) = 0.0
               dyy_psi_b(x, y, z) = 0.0

               dyz_psi_s(x, y, z) = 0.0
               dyy_psi_s(x, y, z) = 0.0

#endif
#endif         
               do s = 1, 8               
                  xa = x
                  ya = y + cy2d(s)
                  za = z + cz2d(s)

                  x2a = x
                  y2a = y + 2.0*cy2d(s)
                  z2a = z + 2.0*cz2d(s)
#ifdef SINGLEFLUID      
                  dyz_psi_r(x, y, z) = dyz_psi_r(x, y, z) + (1.0/36.0)*(8.0*dy_psi_r_buff(xa,ya,za) &
                       &- dy_psi_r_buff(x2a,y2a,z2a))*cz2d(s)
                  dyy_psi_r(x, y, z) = dyy_psi_r(x, y, z) + (1.0/36.0)*(8.0*dy_psi_r_buff(xa,ya,za) &
                       &- dy_psi_r_buff(x2a,y2a,z2a))*cy2d(s)
#else
#ifdef NOSURFACTANT
                  dyz_psi_r(x, y, z) = dyz_psi_r(x, y, z) + (1.0/36.0)*(8.0*dy_psi_r_buff(xa,ya,za) &
                       &- dy_psi_r_buff(x2a,y2a,z2a))*cz2d(s)
                  dyz_psi_b(x, y, z) = dyz_psi_b(x, y, z) + (1.0/36.0)*(8.0*dy_psi_b_buff(xa,ya,za) &
                       &- dy_psi_b_buff(x2a,y2a,z2a))*cz2d(s)

                  dyy_psi_r(x, y, z) = dyy_psi_r(x, y, z) + (1.0/36.0)*(8.0*dy_psi_r_buff(xa,ya,za) &
                       &- dy_psi_r_buff(x2a,y2a,z2a))*cy2d(s)
                  dyy_psi_b(x, y, z) = dyy_psi_b(x, y, z) + (1.0/36.0)*(8.0*dy_psi_b_buff(xa,ya,za) &
                       &- dy_psi_b_buff(x2a,y2a,z2a))*cy2d(s)
#else
                  dyz_psi_r(x, y, z) = dyz_psi_r(x, y, z) + (1.0/36.0)*(8.0*dy_psi_r_buff(xa,ya,za) &
                       &- dy_psi_r_buff(x2a,y2a,z2a))*cz2d(s)
                  dyz_psi_b(x, y, z) = dyz_psi_b(x, y, z) + (1.0/36.0)*(8.0*dy_psi_b_buff(xa,ya,za) &
                       &- dy_psi_b_buff(x2a,y2a,z2a))*cz2d(s)
                  dyz_psi_s(x, y, z) = dyz_psi_s(x, y, z) + (1.0/36.0)*(8.0*dy_psi_s_buff(xa,ya,za) &
                       &- dy_psi_s_buff(x2a,y2a,z2a))*cz2d(s)


                  dyy_psi_r(x, y, z) = dyy_psi_r(x, y, z) + (1.0/36.0)*(8.0*dy_psi_r_buff(xa,ya,za) &
                       &- dy_psi_r_buff(x2a,y2a,z2a))*cy2d(s)
                  dyy_psi_b(x, y, z) = dyy_psi_b(x, y, z) + (1.0/36.0)*(8.0*dy_psi_b_buff(xa,ya,za) &
                       &- dy_psi_b_buff(x2a,y2a,z2a))*cy2d(s)
                  dyy_psi_s(x, y, z) = dyy_psi_s(x, y, z) + (1.0/36.0)*(8.0*dy_psi_s_buff(xa,ya,za) &
                       &- dy_psi_s_buff(x2a,y2a,z2a))*cy2d(s)

#endif
#endif
               end do

!               print *,y,4*(y-1.5)*(y-1.5)*(y-1.5),dy_psi_r(x,y,z),12*(y-1.5)*(y-1.5),dyy_psi_r(x, y, z)
!               print *,y,dyz_psi_r(x, y, z)

            end do ! z
!!$         if (y .gt. 4) then
!!$            stop
!!$         end if
         end do ! y
      end do !x

#ifdef SINGLEFLUID      
      dyy_psi_r(1,:,:) = dyy_psi_r(nx-1,:,:)
      dyy_psi_r(nx,:,:) = dyy_psi_r(2,:,:)

      dyz_psi_r(1,:,:) = dyz_psi_r(nx-1,:,:)
      dyz_psi_r(nx,:,:) = dyz_psi_r(2,:,:)
#else
#ifdef NOSURFACTANT
      dyy_psi_r(1,:,:) = dyy_psi_r(nx-1,:,:)
      dyy_psi_r(nx,:,:) = dyy_psi_r(2,:,:)
      dyy_psi_b(1,:,:) = dyy_psi_b(nx-1,:,:)
      dyy_psi_b(nx,:,:) = dyy_psi_b(2,:,:)

      dyz_psi_r(1,:,:) = dyz_psi_r(nx-1,:,:)
      dyz_psi_r(nx,:,:) = dyz_psi_r(2,:,:)
      dyz_psi_b(1,:,:) = dyz_psi_b(nx-1,:,:)
      dyz_psi_b(nx,:,:) = dyz_psi_b(2,:,:)
#else
      dyy_psi_r(1,:,:) = dyy_psi_r(nx-1,:,:)
      dyy_psi_r(nx,:,:) = dyy_psi_r(2,:,:)
      dyy_psi_b(1,:,:) = dyy_psi_b(nx-1,:,:)
      dyy_psi_b(nx,:,:) = dyy_psi_b(2,:,:)
      dyy_psi_s(1,:,:) = dyy_psi_s(nx-1,:,:)
      dyy_psi_s(nx,:,:) = dyy_psi_s(2,:,:)

      dyz_psi_r(1,:,:) = dyz_psi_r(nx-1,:,:)
      dyz_psi_r(nx,:,:) = dyz_psi_r(2,:,:)
      dyz_psi_b(1,:,:) = dyz_psi_b(nx-1,:,:)
      dyz_psi_b(nx,:,:) = dyz_psi_b(2,:,:)
      dyz_psi_s(1,:,:) = dyz_psi_s(nx-1,:,:)
      dyz_psi_s(nx,:,:) = dyz_psi_s(2,:,:)
#endif
#endif


#ifdef SINGLEFLUID      
      deallocate(dy_psi_r_buff)
#else
#ifdef NOSURFACTANT
      deallocate(dy_psi_r_buff)
      deallocate(dy_psi_b_buff)
#else
      deallocate(dy_psi_r_buff)
      deallocate(dy_psi_b_buff)
      deallocate(dy_psi_s_buff)
#endif
#endif
 
    end subroutine lbe_gradient2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



#ifdef SINGLEFLUID
  subroutine lbe_axisym_correction_sc_forces(psi_mat,&
       & dy_psi_r,dyz_psi_r,dyy_psi_r, &
       & x, y, z, f_r)
#else
#ifdef NOSURFACTANT
  subroutine lbe_axisym_correction_sc_forces(psi_mat,&
       & dy_psi_r,dyz_psi_r,dyy_psi_r, &
       & dy_psi_b,dyz_psi_b,dyy_psi_b, &
       & x, y, z, f_b, f_r)
#else
  subroutine lbe_axisym_correction_sc_forces(psi_mat,&
       & dy_psi_r,dyz_psi_r,dyy_psi_r, &
       & dy_psi_b,dyz_psi_b,dyy_psi_b, &
       & dy_psi_s,dyz_psi_s,dyy_psi_s, &
       & x, y, z, f_b, f_r, f_s)
#endif
#endif

    implicit none

    integer, intent(in) :: x, y, z
    integer :: xa, ya, za, s, ss

    real*8 :: psi_r, psi_b, psi_s
    real*8, dimension(1:,1:,1:,1:), intent(in) :: psi_mat

#ifdef SINGLEFLUID
    real*8, dimension(1:,1:,1:,1:),intent(inout) :: f_r
    real*8, dimension(1:,1:,1:), intent(in) :: dy_psi_r, dyz_psi_r, dyy_psi_r
#else
#ifdef NOSURFACTANT
    real*8, dimension(1:,1:,1:,1:),intent(inout) :: f_r
    real*8, dimension(1:,1:,1:,1:),intent(inout) :: f_b
    real*8, dimension(1:,1:,1:), intent(in) :: dy_psi_r, dyz_psi_r, dyy_psi_r
    real*8, dimension(1:,1:,1:), intent(in) :: dy_psi_b, dyz_psi_b, dyy_psi_b
#else
    real*8, dimension(1:,1:,1:,1:),intent(inout) :: f_r
    real*8, dimension(1:,1:,1:,1:),intent(inout) :: f_b
    real*8, dimension(1:,1:,1:,1:),intent(inout) :: f_s
    real*8, dimension(1:,1:,1:), intent(in) :: dy_psi_r, dyz_psi_r, dyy_psi_r
    real*8, dimension(1:,1:,1:), intent(in) :: dy_psi_b, dyz_psi_b, dyy_psi_b
    real*8, dimension(1:,1:,1:), intent(in) :: dy_psi_s, dyz_psi_s, dyy_psi_s
#endif
#endif
    
    real*8,  dimension(3) :: f_rr, f_bb
    real*8,  dimension(3) :: f_br, f_rb
    real*8,  dimension(3) :: f_cs, f_sc

#ifndef NOSURFACTANT
    real*8,  dimension(3) :: f_ss
#endif

    real*8 :: invy 
    real*8 :: g_rb, g_bb
    
    g_rb = g_br
    g_bb = g_rr
    invy = 1.0/(dfloat(y)-1.5)


#ifndef NOSURFACTANT
    f_ss = 0.d0
    psi_s = 0.d0
#endif

    f_rr = 0.d0
#ifndef SINGLEFLUID
    f_bb = 0.d0
    f_br = 0.d0
    f_rb = 0.d0
#endif

    psi_r = psi_mat(x,y,z,1)


    f_rr(1) = 0.0_rk
    f_rr(2) = (2.0)*-g_rr*psi_r*(invy*dyy_psi_r(x,y,z) - invy*invy*dy_psi_r(x,y,z))
    f_rr(3) = (2.0)*-g_rr*psi_r*invy*dyz_psi_r(x,y,z)
    
    f_r(x,y,z,:) =  f_r(x,y,z,:) + f_rr

#ifndef SINGLEFLUID
    psi_b = psi_mat(x,y,z,2)

    f_rb(1) = 0.0_rk
    f_rb(2) = (2.0)*-g_rb*psi_r*(invy*dyy_psi_b(x,y,z) - invy*invy*dy_psi_b(x,y,z))
    f_rb(3) = (2.0)*-g_rb*psi_r*invy*dyz_psi_b(x,y,z)

    f_bb(1) = 0.0_rk
    f_bb(2) = (2.0)*-g_bb*psi_b*(invy*dyy_psi_b(x,y,z) - invy*invy*dy_psi_b(x,y,z))
    f_bb(3) = (2.0)*-g_bb*psi_b*invy*dyz_psi_b(x,y,z)

    f_br(1) = 0.0_rk
    f_br(2) = (2.0)*-g_br*psi_b*(invy*dyy_psi_r(x,y,z) - invy*invy*dy_psi_r(x,y,z))
    f_br(3) = (2.0)*-g_br*psi_b*invy*dyz_psi_r(x,y,z)

    f_r(x,y,z,:) =  f_r(x,y,z,:) + f_rb
    f_b(x,y,z,:) =  f_b(x,y,z,:) + f_bb + f_br
#ifndef NOSURFACTANT
    psi_s = psi_mat(x,y,z,3)
    !! not implemented
    f_s(x,y,z,:) =  f_s(x,y,z,:)
#endif
#endif
    
  end subroutine lbe_axisym_correction_sc_forces


!!$!!! Boundary conditions tailored for axisymmetric module
  subroutine lbe_axisym_axis_bc(N)
    implicit none
    type(lbe_site), dimension(0:,0:,0:),intent(inout) :: N
    integer :: x, y, z
    integer :: xa, ya, za, xb, yb, zb

!!$    integer :: bcz
    integer :: s,i

!!$    real*8 :: init, finl, diff
    integer :: xx, yy, zz
    real*8, dimension(19) :: tmp_Nn_r,tmp_Nn_b,tmp_Nn_s

    !! there is a possibility of change on N(x,y=1,ny,z)
    !! if bcsel is not 0,1,2
    if (bcsel .eq. 3) then

    y = 1
    yy =2 
    do z=1,nz
       do x=1,nx
          bb_bc_nonrest_1: do s = 1, nnonrest
             N(x, y, z)%n_r(s) = N(x,yy,z)%n_r(reflect_y(s)) 
#ifndef SINGLEFLUID
             N(x, y, z)%n_b(s) = N(x,yy,z)%n_b(reflect_y(s)) 
#endif
#ifndef NOSURFACTANT
             N(x, y, z)%n_s(s) = N(x,yy,z)%n_s(reflect_y(s)) 
#endif
          end do bb_bc_nonrest_1
       end do
    end do

    y = ny
    yy =ny-1
    do z=1,nz
       do x=1,nx
          bb_bc_nonrest_ny: do s = 1, nnonrest
             xx = x + cx(s) 
             zz = z + cz(s)
             N(x, y, z)%n_r(s) = N(xx,yy,zz)%n_r(bounce(s)) 
#ifndef SINGLEFLUID
             N(x, y, z)%n_b(s) = N(xx,yy,zz)%n_b(bounce(s)) 
#endif
#ifndef NOSURFACTANT
             N(x, y, z)%n_s(s) = N(xx,yy,zz)%n_s(bounce(s)) 
#endif
          end do bb_bc_nonrest_ny
       end do
    end do

  end if ! bcsel .eq. 3

  end subroutine lbe_axisym_axis_bc

  subroutine lbe_calculate_axisym_sc_forces(N)
    implicit none
    integer :: x, y, z
!!$    integer :: xa, ya, za, xb, yb, zb
!!$    integer :: bcz
!!$    integer :: s,i
    type(lbe_site), dimension(0:,0:,0:),intent(inout) :: N
!!$    real*8 :: init, finl, diff
!!$    integer :: xx, yy, zz
!!$  
  end subroutine lbe_calculate_axisym_sc_forces


  subroutine lbe_axisym_correction_mom(N,u,rho)
    implicit none
    integer :: x, y, z
    integer :: s
    type(lbe_site), dimension(0:,0:,0:),intent(inout) :: N
    real*8, dimension(0:,0:,0:,0:,0:),intent(in) :: u
    real*8, dimension(0:,0:,0:,0:),intent(in) :: rho

    ! Number densities
    real*8 :: rN_r, rN_b, rN_s
    real*8 :: invy, visco_r, visco_b, visco_s 
    real*8 :: t1, t2, t3, t4
    real*8 :: dz_w, dy_w, dz_v, dy_v
    real*8, dimension(3) :: p_r, p_b, p_s
    real*8, dimension(3) :: mom_corr
  
   
    ! Velocities
    real*8, allocatable, dimension(:,:,:,:) :: u_r, u_b, u_s


    ! Weighted masses
    real*8, allocatable, dimension(:,:,:) :: rho_r, rho_b, rho_s
    real*8, allocatable, dimension(:,:,:,:) :: F_r, F_b, F_s


    if (.not. allocated(u_r)) allocate(u_r(1:3,1:nx,-1:ny+2,-1:nz+2))
    if (.not. allocated(rho_r)) allocate(rho_r(1:nx,1:ny,1:nz))
    if (.not. allocated(F_r)) allocate(F_r(1:19,1:nx,1:ny,1:nz))
#ifndef SINGLEFLUID
    if (.not. allocated(u_b)) allocate(u_b(1:3,1:nx,-1:ny+2,-1:nz+2))
    if (.not. allocated(rho_b)) allocate(rho_b(1:nx,1:ny,1:nz))
    if (.not. allocated(F_b)) allocate(F_b(1:19,1:nx,1:ny,1:nz))
#endif
#ifndef NOSURFACTANT
    if (.not. allocated(u_s)) allocate(u_s(1:3,1:nx,-1:ny+2,-1:nz+2))
    if (.not. allocated(rho_s)) allocate(rho_s(1:nx,1:ny,1:nz))
    if (.not. allocated(F_s)) allocate(F_s(1:19,1:nx,1:ny,1:nz))
#endif

    !pbc x
    N(1,:,:) = N(3,:,:)
    N(4,:,:) = N(2,:,:)


    do z=1,nz
       do y=1,ny
          do x=1,nx

             p_r = 0.d0
             p_b = 0.d0
             p_s = 0.d0
             rN_r = 0.d0
             rN_b = 0.d0
             rN_s = 0.d0
             
             F_r(:,x,y,z) = N(x,y,z)%n_r(:)*g(:)
#ifndef SINGLEFLUID
             F_b(:,x,y,z) = N(x,y,z)%n_b(:)*g(:)
#endif
#ifndef NOSURFACTANT
             F_s(:,x,y,z) = N(x,y,z)%n_s(:)*g(:)
#endif

             do s = 1, nvecs
                p_r = p_r + F_r(s,x, y, z)*c(s, :)
                rN_r = rN_r + F_r(s,x, y, z)
#ifndef SINGLEFLUID
                p_b = p_b +   F_b(s,x, y, z) * c(s, :)
                rN_b = rN_b +  F_b(s,x, y, z)
#endif
#ifndef NOSURFACTANT
                p_s = p_s +   F_s(s,x, y, z) * c(s, :)
                rN_s = rN_s +  F_(s,x, y, z)
#endif
             end do

             !             rho_r(x,y,z) = amass_r * rN_r
             !             u_r(:,x,y,z) = p_r/rho_r(x,y,z)

             rho_r(x,y,z) = 0.5*(amass_r * rN_r + rho(x,y,z,2))
             u_r(:,x,y,z) = 0.5*(p_r/(amass_r * rN_r) + u(:,x,y,z,1))
#ifndef SINGLEFLUID
             !             rho_b(x,y,z) = amass_b * rN_b
             !             u_b(:,x,y,z) = p_b/rho_b(x,y,z)
             
             rho_b(x,y,z) = 0.5*(amass_b * rN_b + rho(x,y,z,2))
             u_b(:,x,y,z) = 0.5*(p_b/(amass_b * rN_b) + u(:,x,y,z,2))

#endif
#ifndef NOSURFACTANT
             !             rho_s(x,y,z) = amass_s * rN_s
             !             u_s(:,x,y,z) = p_s/rho_s(x,y,z)

             rho_s(x,y,z) = 0.5*(amass_s * rN_s + rho(x,y,z,3))
             u_s(:,x,y,z) = 0.5*(p_s/(amass_s * rN_s) + u(:,x,y,z,3))

#endif

          end do
       end do
    end do

    ! boundary condition based velocities at ghost nodes
    ! bc - along y-axis
    u_r(1,:,-1,:) =  u_r(1,:,4,:)
    u_r(2,:,-1,:) = -u_r(2,:,4,:)
    u_r(3,:,-1,:) =  u_r(3,:,4,:)

    u_r(1,:,0,:) =  u_r(1,:,3,:)
    u_r(2,:,0,:) = -u_r(2,:,3,:)
    u_r(3,:,0,:) =  u_r(3,:,3,:)

    u_r(1,:,1,:) =  u_r(1,:,2,:)
    u_r(2,:,1,:) = -u_r(2,:,2,:)
    u_r(3,:,1,:) =  u_r(3,:,2,:)

    u_r(:,:,ny,:)   = u_r(:,:,ny-1,:)
    u_r(:,:,ny+1,:) = u_r(:,:,ny-1,:)
    u_r(:,:,ny+2,:) = u_r(:,:,ny-1,:)

    rho_r(:,1,:) = rho_r(:,2,:)
    rho_r(:,ny,:) = rho_r(:,ny-1,:)


    !z-axis pbc
    u_r(:,:,:,-1) = u_r(:,:,:,nz-3)
    u_r(:,:,:,0) = u_r(:,:,:,nz-2)
    u_r(:,:,:,1) = u_r(:,:,:,nz-1)
    u_r(:,:,:,nz) = u_r(:,:,:,2)
    u_r(:,:,:,nz+1) = u_r(:,:,:,3)
    u_r(:,:,:,nz+2) = u_r(:,:,:,4)

    rho_r(:,:,1) = rho_r(:,:,nz-1)
    rho_r(:,:,nz) = rho_r(:,:,2)

#ifndef SINGLEFLUID
    !y-axis
    u_b(1,:,-1,:) =  u_b(1,:,4,:)
    u_b(2,:,-1,:) = -u_b(2,:,4,:)
    u_b(3,:,-1,:) =  u_b(3,:,4,:)

    u_b(1,:,0,:) =  u_b(1,:,3,:)
    u_b(2,:,0,:) = -u_b(2,:,3,:)
    u_b(3,:,0,:) =  u_b(3,:,3,:)

    u_b(1,:,1,:) =  u_b(1,:,2,:)
    u_b(2,:,1,:) = -u_b(2,:,2,:)
    u_b(3,:,1,:) =  u_b(3,:,2,:)

    u_b(:,:,ny+1,:) = u_b(:,:,ny,:)
    u_b(:,:,ny+2,:) = u_b(:,:,ny,:)

    rho_b(:,1,:) = rho_b(:,2,:)
    rho_b(:,ny,:) = rho_b(:,ny-1,:)

    !z-axis
    u_b(:,:,:,-1) = u_b(:,:,:,nz-3)
    u_b(:,:,:,0) = u_b(:,:,:,nz-2)
    u_b(:,:,:,1) = u_b(:,:,:,nz-1)
    u_b(:,:,:,nz) = u_b(:,:,:,2)
    u_b(:,:,:,nz+1) = u_b(:,:,:,3)
    u_b(:,:,:,nz+2) = u_b(:,:,:,4)

    rho_b(:,:,1) = rho_b(:,:,nz-1)
    rho_b(:,:,nz) = rho_b(:,:,2)

#endif

    
    do z=1,nz
       do y=1,ny
          do x=1,nx
             invy = 1.0/(dfloat(y)-1.5)

             t1=(8.0*u_r(3,x,y,z+1) - u_r(3,x,y,z+2));
             t2=(8.0*u_r(3,x,y+1,z) - u_r(3,x,y+2,z));
             t3=(8.0*u_r(3,x,y,z-1) - u_r(3,x,y,z-2));
             t4=(8.0*u_r(3,x,y-1,z) - u_r(3,x,y-2,z));

             !u_r(3:)  = w
             dz_w = (1.0/12.0)*(t1-t3)
             dy_w = (1.0/12.0)*(t2-t4)

             t1=(8.0*u_r(2,x,y,z+1) - u_r(2,x,y,z+2));
             t2=(8.0*u_r(2,x,y+1,z) - u_r(2,x,y+2,z));
             t3=(8.0*u_r(2,x,y,z-1) - u_r(2,x,y,z-2));
             t4=(8.0*u_r(2,x,y-1,z) - u_r(2,x,y-2,z));

             dz_v = (1.0/12.0)*(t1-t3)
             dy_v = (1.0/12.0)*(t2-t4)


             !! do the same for blue fluid

             visco_r = cs2*(get_tau_r(N,x,y,z)-0.5)
             mom_corr(1) = 0.0
             mom_corr(2) = invy*rho_r(x,y,z)*(2.0*visco_r*(dy_v - invy*u_r(2,x,y,z)) -u_r(2,x,y,z)*u_r(2,x,y,z))
             mom_corr(3) = invy*rho_r(x,y,z)*(visco_r*(dy_w + dz_v) - u_r(3,x,y,z)*u_r(2,x,y,z))
#ifndef SINGLEFLUID
             t1=(8.0*u_b(3,x,y,z+1) - u_b(3,x,y,z+2));
             t2=(8.0*u_b(3,x,y+1,z) - u_b(3,x,y+2,z));
             t3=(8.0*u_b(3,x,y,z-1) - u_b(3,x,y,z-2));
             t4=(8.0*u_b(3,x,y-1,z) - u_b(3,x,y-2,z));


             dz_w = (1.0/12.0)*(t1-t3)
             dy_w = (1.0/12.0)*(t2-t4)

             t1=(8.0*u_b(2,x,y,z+1) - u_b(2,x,y,z+2));
             t2=(8.0*u_b(2,x,y+1,z) - u_b(2,x,y+2,z));
             t3=(8.0*u_b(2,x,y,z-1) - u_b(2,x,y,z-2));
             t4=(8.0*u_b(2,x,y-1,z) - u_b(2,x,y-2,z));

             dz_v = (1.0/12.0)*(t1-t3)
             dy_v = (1.0/12.0)*(t2-t4)

             visco_b = (get_tau_b(N,x,y,z)-0.5)/3.0
             mom_corr(1) = 0.0
             mom_corr(2) = invy*rho_b(x,y,z)*(2.0*visco_b*(dy_v - invy*u_b(2,x,y,z)) -u_b(2,x,y,z)*u_b(2,x,y,z))
             mom_corr(3) = invy*rho_b(x,y,z)*(visco_b*(dy_w + dz_v) - u_b(3,x,y,z)*u_b(2,x,y,z))
#endif
#ifndef NOSURFACTANT
             visco_s = (get_tau_s(N,x,y,z)-0.5)/3.0
             mom_corr(1) = 0.0
             mom_corr(2) = invy*rho_s(x,y,z)*(2.0*visco_s*(dy_v - invy*u_s(2,x,y,z)) -u_s(2,x,y,z)*u_s(2,x,y,z))
             mom_corr(3) = invy*rho_s(x,y,z)*(visco_s*(dy_w + dz_v) - u_s(3,x,y,z)*u_s(2,x,y,z))
#endif

             if ( N(x,y,z)%rock_state == 0. ) then


                do s = 1, nvecs 
!!!                   N(x,y,z)%n_r(s) = N(x,y,z)%n_r(s) + cs2*sum(c(s,:)*mom_corr)
                   F_r(s,x,y,z) = F_r(s,x,y,z) + weight(s)*4.0*sum(c(s,:)*mom_corr)/cs2
                   ! caution : unknown prefactor 4 ?
#ifndef SINGLEFLUID
!!!                   N(x,y,z)%n_b(s) = N(x,y,z)%n_b(s) + cs2*sum(c(s,:)*mom_corr)
                   F_b(s,x,y,z) = F_b(s,x,y,z) + weight(s)*4.0*sum(c(s,:)*mom_corr)/cs2
#endif
#ifndef NOSURFACTANT
!!!                   N(x,y,z)%n_s(s) = N(x,y,z)%n_s(s) + cs2*sum(c(s,:)*mom_corr)
                   F_s(s,x,y,z) = F_s(s,x,y,z) + weight(s)*4.0*sum(c(s,:)*mom_corr)/cs2
#endif
                end do

                N(x,y,z)%n_r(:) = F_r(:,x,y,z)/g(:)
#ifndef SINGLEFLUID
                N(x,y,z)%n_b(:) = F_b(:,x,y,z)/g(:)
#endif
#ifndef NOSURFACTANT
                N(x,y,z)%n_s(:) = F_s(:,x,y,z)/g(:)
#endif

             endif
          end do
       end do
    end do


    deallocate(u_r)
    deallocate(rho_r)
    deallocate(F_r)
#ifndef SINGLEFLUID
    deallocate(u_b)
    deallocate(rho_b)
    deallocate(F_b)
#endif
#ifndef NOSURFACTANT
    deallocate(u_s)
    deallocate(rho_s)
    deallocate(F_s)
#endif


  end subroutine lbe_axisym_correction_mom


  subroutine lbe_axisym_correction_cty(N,u,rho)
    implicit none
    integer :: x, y, z
    integer :: xa, ya, za, xb, yb, zb
    integer :: bcz
    integer :: s,i
    type(lbe_site), dimension(0:,0:,0:),intent(inout) :: N

    real*8, dimension(0:,0:,0:,0:)   ,intent(in) :: rho
    real*8, dimension(0:,0:,0:,0:,0:),intent(in) :: u
                     !(x/y/z, x,y,z, r/b/s)!  


    real*8, allocatable, dimension(:,:,:,:) :: F_r, F_b, F_s

    ! Velocities
    real*8, dimension(3) :: p_r, p_b, p_s
    real*8, dimension(3) :: u_r, u_b, u_s

    ! Number densities
    real*8 :: bN_s, rN_s, sN_s
 
    ! Weighted masses
    real*8 :: rho_b, rho_r, rho_s
    
    ! inverse of distance from axis (r == 0)
    real*8 :: invy

    !pbc x ! better to do in lbe.F90 
    N(1,:,:) = N(3,:,:)
    N(4,:,:) = N(2,:,:)


    do z=1,nz
       do y=1,ny
          do x=1,nx
             
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
             

             !             rho_r = amass_r * rN_r
             !             u_r = p_r/rho_r(x,y,z)
             
             rho_r = 0.5*(amass_r * rN_s + rho(x,y,z,1))
             u_r = 0.5*(p_r/(amass_r * rN_s) + u(:,x,y,z,1))
#ifndef SINGLEFLUID
             !             rho_b = amass_b * bN_s
             !             u_b = p_b/rho_b
             
             rho_b = 0.5*(amass_b * bN_s + rho(x,y,z,2))
             u_b = 0.5*(p_b/(amass_b * bN_s) + u(:,x,y,z,2))
             
#endif
#ifndef NOSURFACTANT
             !             rho_s(x,y,z) = amass_s * sN_s
             !             u_s(:,x,y,z) = p_s/rho_s(x,y,z)

             rho_s = 0.5*(amass_s * sN_s + rho(x,y,z,3))
             u_s = 0.5*(p_s/(amass_s * sN_s) + u(:,x,y,z,3))

#endif


             invy = 1.0/(dfloat(y)-1.5)

             if ( N(x,y,z)%rock_state == 0. ) then
                do s = 1, nvecs 

                   N(x,y,z)%n_r(s) = N(x,y,z)%n_r(s) - invy*rho_r*u_r(2)/sum(g)
                   
#ifndef SINGLEFLUID
                   N(x,y,z)%n_b(s) = N(x,y,z)%n_b(s) - invy*rho_b*u_b(2)/sum(g)
#endif
#ifndef NOSURFACTANT
                   N(x,y,z)%n_s(s) = N(x,y,z)%n_s(s) - invy*rho_s*u_s(2)/sum(g)
#endif
                end do

             endif

          end do     !x
       end do        !y
    end do           !z

  end subroutine lbe_axisym_correction_cty

end module lbe_axisym_module
#endif ! AXISYM
