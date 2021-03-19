#include "lbe.h"

!> Code responsible for initialising the system.
!>
!> Which option is taken depends on the \c init_cond line in the input file.
!> Valid values are:
!>
!> \li -5 = \c VEL Z	- gives each site the occupation specified
!> \li -4 = \c SIMPLE	- gives each site the occupation specified
!> \li -3 = \c INIT_RATIO	- gives each site the occupation specified
!>                                in fr, fg, and fb.
!> \li -1 = \c INIT_RAND
!> \li  1 = \c INIT_LAMX	- Oil-surf monolayer-water-surf monolayer in X dir
!> \li  2 = \c INIT_LAMY	- Ditto, but with lamellae perpendicular to Y dir
!> \li  3 = \c INIT_LAMZ	- Ditto, but with lamellae perpendicular to Z dir
!> \li  7 = \c INIT_RESTORE    - Restores the entire lattice state from one
!>                              previously saved.
!> \li  9 = \c INIT_UPSCALE    - Upscale from smaller system
!> \li 10 = \c INIT_CUTOUT    - put previously saved state at any point in the
!>                              system
!> \li 11 = \c INIT_BI_LAM_X
!> \li 12 = \c INIT_BI_LAM_Y
!> \li 13 = \c INIT_BI_LAM_Z
!> \li 14 = \c Binary droplet
!> \li 16 = \c INIT_BI_SIN_LAM_X

module lbe_init_module
  use lbe_globals_module
  use lbe_helper_module, only: is_restoring
  use lbe_io_checkpoint_module, only: restore_checkpoint
  use lbe_io_module, only: lbe_get_lbe_input
#ifdef USEXDRF
  use lbe_io_xdrf_module, only: restore_cutout,restore_upscale
#endif
  use lbe_log_module
#ifdef MD
  use lbe_md_globals_module, only: provide_uid2i
#endif
  use lbe_parallel_module, only: abend,ccoords,check_allocate,start,tnx,tny,tnz&
       &,comm_cart
  use lbe_parms_module
  use lbe_timer_module, only: register_timer
  ! use lbe_collision_module ! Only needed for INIT_HUDONG
  ! use lbe_collision_simple_module
  use lbe_bdist_module, only: boltz_dist, mrt_init_dist
  use lbe_types_module
  use lbe_init_rock_module, only: lbe_init_rock
  use luxury_rng_module, only: ranlux,rluxgo

! Intermediate Knudsen initialisation requires single halo exchange
#ifdef LOCALBC
  use lbe_parallel_module, only: halo, halo_exchange, acccoef_halo
#endif
#ifdef VARTAU 
  use lbe_parallel_module, only: halo, halo_exchange, vartau_halo
#endif

  implicit none
  include 'mpif.h'

  private

  public lbe_init_system, lbe_init_timers, lbe_unique_seeds, lbe_setup_fluxz
  public lbe_unique_id

contains

!> Globally set the unique ID for the simulation.
subroutine lbe_unique_id()
  implicit none
  integer :: ierror

  if ( myrankc == 0 ) then
    call system_clock(chk_uid)
    write(msgstr,"('Generated UID <',I10.10,'>.')") chk_uid
    call log_msg(msgstr)
    
    open (unit=2, file="chk_uid-info", status="unknown",position="append",action="write")
    write(2,'(I10.10)') chk_uid
    close(2)
  end if

  call MPI_Bcast(chk_uid,1,MPI_INTEGER,0,comm_cart,ierror)

end subroutine lbe_unique_id

!> This is the general initialization routine.
!>
!> It calls one of the other init routines, depending on the value of
!> the variable \c init_cond, specified in the input-file. The state
!> array \c lbe_N / \c whole_N , on return is filled with the initial
!> values.
!>
!> \note The halo region is zero'd out but still \b not defined after
!> this function returns - it must be filled in by a halo exchange
!> later.
!>
!> \param[in,out] lbe_N local chunk of the lattice with halo extent 1
!> (old LB3D style, should be removed sooner or later)
!>
!> \param[in,out] whole_N local chunk of the lattice with full halo of
!> depth \c halo_extent
subroutine lbe_init_system(lbe_N,whole_N)
  type(lbe_site),intent(inout) :: lbe_N(0:,0:,0:)
  type(lbe_site),intent(inout) :: &
       &whole_N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
  integer :: x,y,z
  type(radial_model), dimension(3) :: model
#ifdef INTERPOLATEDBBbla
  real(kind=rk) :: px(3)
  real(kind=rk) :: ppx(3,14),nxhf,nyhf,nzhf,nxreal,nyreal,nzreal
  real(kind=rk) :: ltemp(3), ltempsq
  integer :: xa,ya,za
  real(kind=rk) :: l(3),d(3),rop(2)
  real(kind=rk) :: lsq,msq,dist,dnorm,sub
  integer :: s,np
#endif
  call log_msg_hdr("Initializing system")

  ! Zero out the whole lattice, including the full halo
  do x=1-halo_extent,nx+halo_extent
     do y=1-halo_extent,ny+halo_extent
        do z=1-halo_extent,nz+halo_extent
           whole_N(x,y,z)%n_r = 0.0_rk
#ifndef SINGLEFLUID
           whole_N(x,y,z)%n_b = 0.0_rk
#endif
           whole_N(x,y,z)%rock_state = 0.0_rk
#ifndef NOSURFACTANT
           whole_N(x,y,z)%n_s = 0.0_rk
           whole_N(x,y,z)%d = 0.0_rk
#endif

#ifdef INTERPOLATEDBB
           whole_N(x,y,z)%delta_q = 0.5_rk
#endif

#ifdef VARTAU
           whole_N(x,y,z)%taupos_r = tau_r
#ifndef SINGLEFLUID
           whole_N(x,y,z)%taupos_b = tau_b
#endif
#ifndef NOSURFACTANT
           whole_N(x,y,z)%taupos_s = tau_s
#endif
#endif

!#ifdef COMMON_VEL_FIX
#ifndef OLD_VEL
				 whole_N(x,y,z)%n_r_pre = 0.0_rk
#ifndef SINGLEFLUID
             whole_N(x,y,z)%n_b_pre = 0.0_rk
#endif
#ifndef NOSURFACTANT
             whole_N(x,y,z)%n_s_pre = 0.0_rk
#endif
#endif

        end do
     end do
  end do

!  px(1)=nx/2.0
!  px(2)=ny/2.0
!  px(3)=nz/2.0
!  rop = (/4.5, 4.5/)
!   do x=1,nx
!     do y=1,ny
!        do z=1,nz
!          !if((x-px(1))**2+(z-px(3))**2>rop(1)**2) then
!          !if((x-px(1))**2+(y-px(2))**2>rop(1)**2) then
!          if((y-px(2))**2>rop(1)**2) then
!          !if((x-px(1))**2>rop(1)**2.or.(y-px(2))**2>rop(1)**2) then
!          !if((x-px(1))**2<rop(1)**2.and.(y-px(2))**2<rop(1)**2.and.(z-px(3))**2<rop(1)**2) then
!          !if(((x-px(1))**2+(y-px(2))**2+(z-px(3))**2).le.rop(1)**2) then
!            whole_N(x,y,z)%rock_state = 5.0
!          end if
!        end do
!     end do
!  end do
#ifdef INTERPOLATEDBBbla


! Circular channel flow
  px(1)=nx/2.0
  px(2)=ny/2.0
  px(3)=nz/2.0
  rop = (/R_test, R_test/)  ! (/R_orth,R_para/)
!  rop = (/9.5, 9.5/)  ! (/R_orth,R_para/)
!  rop = (/9.00, 9.000/)  ! (/R_orth,R_para/)
!  rop = (/18.4282362335, 18.4282362335/)  ! (/R_orth,R_para/)
!  rop = (/18.1402950424, 18.1402950424/)  ! (/R_orth,R_para/)
!  rop = (/45,45/)  ! (/R_orth,R_para/)
!  rop = (/45,45/)  ! (/R_orth,R_para/)
!  rop = (/15.0005,15.0005/)  ! (/R_orth,R_para/)
!  nxhf = (real(nx,kind=rk)+1.0)/2.0
!  nyhf = (real(ny,kind=rk)+1.0)/2.0
!  nzhf = (real(nz,kind=rk)+1.0)/2.0
!  nxreal = real(nx,kind=rk)
!  nyreal = real(ny,kind=rk)
!  nzreal = real(nz,kind=rk)
!
!  ppx(1,:) = (/nxhf,nxhf,nxhf,nxhf,1.0_rk,nxreal&
!  &,1.0_rk,1.0_rk,nxreal,nxreal,1.0_rk,1.0_rk,nxreal,nxreal/)
!  ppx(2,:) = (/nyhf,nyhf,1.0_rk,nyreal,nyhf,nyhf&
!  &,1.0_rk,nyreal,1.0_rk,nyreal,1.0_rk,nyreal,1.0_rk,nyreal/)
!  ppx(3,:) = (/1.0_rk,nzreal,nzhf,nzhf,nzhf,nzhf&
!  &,1.0_rk,nzreal,nzreal,1.0_rk,nzreal,1.0_rk,1.0_rk,nzreal/)

!  do np=1,14
!   px = ppx(:,np)
!   do x=1,nx
!     do y=1,ny
!        do z=1,nz
!          !if((x-px(1))**2+(z-px(3))**2>rop(1)**2) then
!          !if((x-px(1))**2+(y-px(2))**2>rop(1)**2) then
!          !if((y-px(2))**2>rop(1)**2) then
!          !if((x-px(1))**2>rop(1)**2.or.(y-px(2))**2>rop(1)**2.or.(z-px(3))**2>rop(1)**2) then
!          !if((x-px(1))**2<rop(1)**2.and.(y-px(2))**2<rop(1)**2.and.(z-px(3))**2<rop(1)**2) then
!          if(((x-px(1))**2+(y-px(2))**2+(z-px(3))**2).lt.rop(1)**2) then
!            whole_N(x,y,z)%rock_state = 5.0
!          end if
!        end do
!     end do
!  end do
!  end do


!  do np=1,14
!  px = ppx(:,np)
!   write(*,"(2F16.10)") rop
  do z=1,nz
     do y=1,ny
        do x=1,nx
            l = px-real((/x,y,z/),kind=rk)
            lsq = dot_product(l,l)
          do s=1,nnonrest
            xa = x + cx(s)
            ya = y + cy(s)
            za = z + cz(s)
          !  ! For sphere comment from here
          !  l = px-real((/xa,ya,za/),kind=rk)
          !  !l(2) = 0.0_rk
          !  l(3) = 0.0_rk
          !  lsq = dot_product(l,l)
          !  ! till here
            
            if (whole_N(xa,ya,za)%rock_state .ne. 0) then
              if (whole_N(x,y,z)%rock_state ==0) then
          !    ltemp = px-real((/xa,ya,za/),kind=rk)
          !    ltempsq = dot_product(ltemp,ltemp)
          !    if(ltempsq<rop(1)**2) then
             !!  d = real((/cx(s),cy(s),cz(s)/),kind=rk)
             !!  dnorm  = sqrt(dot_product(d,d))
             !!  d = d/dnorm
             !!  dist = dot_product(l,d)
             !!  msq = lsq-dist**2
             !!  sub = sqrt(rop(1)**2-msq) ! rop(1) implies sphere without orientation
               whole_N(x,y,z)%delta_q(s) = 0.9_rk !(dist-sub)/dnorm
      
             !   !d = real((/-cx(s),0,-cz(s)/),kind=rk)
             !   d = real((/-cx(s),-cy(s),0/),kind=rk)
             !   dnorm  = sqrt(dot_product(d,d))
             !   d = d/dnorm
             !   dist = dot_product(l,d)
             !   msq = lsq-dist**2
             !   sub = sqrt(rop(1)**2-msq) ! rop(1) implies sphere without orientation
             !   whole_N(x,y,z)%delta_q(s) = 1-((dist-sub)/dnorm)
             !   !whole_N(x,y,z)%delta_q(s) = 0.95_rk !0.00005_rk ! 1-((dist-sub)/dnorm)
      
                
      !          N(x,y,z)%delta_q(s) = 0.5_rk
      !          write(msgstr&
      !          write(*&
      !               &,"('x = ',I0,', y = ',I0,', z = ',I0,', s=',I0,' dist = ',F16.10)")&
      !               & x,y,z,s,whole_N(x,y,z)%delta_q(s)
      !          call log_msg(msgstr)
                if(whole_N(x,y,z)%delta_q(s)<0.or.whole_N(x,y,z)%delta_q(s)>1.or.isnan(whole_N(x,y,z)%delta_q(s))) then
                   write(*&
                       &,"('pos = ',3F16.10,', lddist = ',F16.10,', msq = ',F16.10,', sub=',F16.10,' dist = 'F16.10)")&
                       & l,dist,msq,sub,whole_N(x,y,z)%delta_q(s)
                  whole_N(x,y,z)%delta_q(s) = 0 
                end if
             ! end if
              end if
            end if
          end do
        end do
     end do
  end do
!  end do

!  write(*, "('Distances')")
!  do z=1,nz
!     do y=1,ny
!        do x=1,nx
!          do s=1,nnonrest
!            xa = x + cx(s)
!            ya = y + cy(s)
!            za = z + cz(s)
!            ! For sphere comment from here
!      !      l = px-real((/xa,ya,za/),kind=rk)
!      !      l(3) = 0.0_rk
!      !      lsq = dot_product(l,l)
!            ! till here
!            if (whole_N(xa,ya,za)%rock_state .ne. 0) then
!              if (whole_N(x,y,z)%rock_state ==0) then
!                l=real((/x,y,z/),kind=rk)+whole_N(x,y,z)%delta_q(s)*real((/cx(s),cy(s),cz(s)/),kind=rk)
!!                write(*,"(3F16.10)") l
!                l=l-px
!                write(*,"(F16.10)") dot_product(l,l)-rop(1)**2
!              end if
!            end if
!          end do 
!        end do
!     end do
!  end do

!  rop = (/4.5,4.5/)  ! (/R_orth,R_para/)
!
!!  do np=1,14
!!   px = ppx(:,np)
!   do x=1,nx
!     do y=1,ny
!        do z=1,nz
!          !if((x-px(1))**2+(y-px(2))**2>rop(1)**2) then
!          if(((x-px(1))**2+(y-px(2))**2+(z-px(3))**2).le.rop(1)**2) then
!            whole_N(x,y,z)%rock_state = 5.0
!          end if
!        end do
!     end do
!  end do
!!  end do
!   write(*,"(2F16.10)") rop
!  do z=1,nz
!     do y=1,ny
!        do x=1,nx
! ! do z=nz/2-6,nz/2+6
! !    do y=ny/2-6,ny/2+6
! !       do x=nx/2-6,nx/2+6
!            l = px-real((/x,y,z/),kind=rk)
!            lsq = dot_product(l,l)
!          do s=1,nnonrest
!            xa = x + cx(s)
!            ya = y + cy(s)
!            za = z + cz(s)
!            ! For sphere comment from here
!          !  l = px-real((/xa,ya,za/),kind=rk)
!          !  l(3) = 0.0_rk
!          !  lsq = dot_product(l,l)
!            ! till here
!            
!            if (whole_N(xa,ya,za)%rock_state .ne. 0) then
!              if (whole_N(x,y,z)%rock_state ==0) then
!              ltemp = px-real((/xa,ya,za/),kind=rk)
!              ltempsq = dot_product(ltemp,ltemp)
!              if(ltempsq<rop(1)**2) then
!                d = real((/cx(s),cy(s),cz(s)/),kind=rk)
!                dnorm  = sqrt(dot_product(d,d))
!                d = d/dnorm
!                dist = dot_product(l,d)
!                msq = lsq-dist**2
!                sub = sqrt(rop(1)**2-msq) ! rop(1) implies sphere without orientation
!                whole_N(x,y,z)%delta_q(s) = (dist-sub)/dnorm
!      
!          !      d = real((/-cx(s),-cy(s),0/),kind=rk)
!          !      dnorm  = sqrt(dot_product(d,d))
!          !      d = d/dnorm
!          !      dist = dot_product(l,d)
!          !      msq = lsq-dist**2
!          !      sub = sqrt(rop(1)**2-msq) ! rop(1) implies sphere without orientation
!          !      whole_N(x,y,z)%delta_q(s) = 1-((dist-sub)/dnorm)
!               ! whole_N(x,y,z)%delta_q(s) = 0.35_rk !0.00005_rk ! 1-((dist-sub)/dnorm)
!      
!                
!      !          N(x,y,z)%delta_q(s) = 0.5_rk
!      !          write(msgstr&
!      !          write(*&
!      !               &,"('x = ',I0,', y = ',I0,', z = ',I0,', s=',I0,' dist = ',F16.10)")&
!      !               & x,y,z,s,whole_N(x,y,z)%delta_q(s)
!      !          call log_msg(msgstr)
!                if(whole_N(x,y,z)%delta_q(s)<0.or.whole_N(x,y,z)%delta_q(s)>1.or.isnan(whole_N(x,y,z)%delta_q(s))) then
!      !            write(msgstr&
!      !             l = real((/x,y,z/),kind=rk)
!                   write(*&
!                       &,"('pos = ',3F16.10,', lddist = ',F16.10,', msq = ',F16.10,', sub=',F16.10,' dist = 'F16.10)")&
!                       & l,dist,msq,sub,whole_N(x,y,z)%delta_q(s)
!                  whole_N(x,y,z)%delta_q(s) = 0 
!                end if
!              end if
!              end if
!            end if
!          end do
!        end do
!     end do
!  end do
#endif

  ! And now fill it with something interesting
  if ( is_restoring() ) then
     call log_msg("Initializing system using INIT_RESTORE.")
     call lbe_init_restore(lbe_N)
  else
     select case(init_cond)
        
     case (-6)
        call log_msg("Initializing system using B_DIST with Z velocity (Poiseuille).")
        call lbe_init_velz(lbe_N)

     case (-5)
        call log_msg("Initializing system using B_DIST with Z velocity.")
        call lbe_init_velz(lbe_N)
        
     case (-4)
        call log_msg("Initializing system using B_DIST.")
        call lbe_init_simple(lbe_N)

     case (-3)
        call log_msg("Initializing system using INIT_RATIO.")
        call lbe_init_ratio(lbe_N)
        
     case (-1)
        call log_msg("Initializing system using INIT_RAND.")
        call lbe_init_rand(lbe_N)
        
     case (0)
        call log_msg("Initializing system using INIT_FRAC.")
        call lbe_init_frac(lbe_N)
        
     case (1)
        call log_msg("Initializing system using INIT_LAMX.")
        call lbe_init_lamx(lbe_N)
        
     case (2)
        call log_msg("Initializing system using INIT_LAMY.")
        call lbe_init_lamy(lbe_N)
        
     case (3)
        call log_msg("Initializing system using INIT_LAMZ.")
        call lbe_init_lamz(lbe_N)
        
    case (6)
       call log_msg("Initializing system using 'old' INIT_BI_LAM_Y from version 5.")
       call lbe_init_bi_lam_y_fromv5(lbe_N)
       
    case (7)
      ! This statement should not be reached, as is_restoring will have taken care of it.
       call error("Error in is_restoring()")
       
    case (9)
       call log_msg("Initializing system using INIT_UPSCALE.")
       call lbe_init_upscale(lbe_N)
       
    case (10)
       call log_msg("Initializing system using INIT_RATIO and INIT_CUTOUT")
       call lbe_init_ratio(lbe_N)
#ifdef USEXDRF
       call restore_cutout(lbe_N)
#else
       call log_msg("XDRF is disabled, can't perform INIT_CUTOUT")
       call Abend
#endif
       
    case (11)
       call log_msg("Initializing system using INIT_BI_LAM_X.")
       call lbe_init_bi_lam_x(lbe_N)
       
    case (12)
       call log_msg("Initializing system using INIT_BI_LAM_Y.")
       call lbe_init_bi_lam_y(lbe_N)
       
    case (13)
       call log_msg("Initializing system using INIT_BI_LAM_Z.")
       call lbe_init_bi_lam_z(lbe_N)
       
    case (14)
       call log_msg("Initializing system using Binary Droplet.")
       model(rad_inner)%n_r = fr
       model(rad_inner)%n_b = fb
       model(rad_inner)%n_s = fg
       model(rad_inner)%dip = fd
       
       model(rad_middle)%n_r = qr
       model(rad_middle)%n_b = qb
       model(rad_middle)%n_s = qg
       model(rad_middle)%dip = qd
       
       model(rad_outer)%n_r = pr
       model(rad_outer)%n_b = pb
       model(rad_outer)%n_s = pg
       model(rad_outer)%dip = pd
       call lbe_init_radial(lbe_N, model)
       
    case (16)
       call log_msg("Initializing system using INIT_BI_SIN_LAM_X.")
       call lbe_init_bi_sin_lam_x(lbe_N)
       
    case (19)
       call log_msg("Initializing system using INIT_CYLINDER_X.")
       call lbe_init_cylinder_x(lbe_N)
       
    case (21)
       call log_msg("Initializing system using INIT_CYLINDER_Z.")
       call lbe_init_cylinder_z(lbe_N)

    case (24)
       !call log_msg("Initializing system using INIT_SINUS.")
       !call lbe_init_sinus(lbe_N)
       !call log_msg("Initializing system using INIT_ratio_plusdisturb.")
       !call lbe_init_ratio_plusdisturb(lbe_N)
       call log_msg("Initializing system using '	old' INIT_FRAC with random number generator from version 5.")
       call lbe_init_frac_fromv5(lbe_N)
       
    case (25)
       call log_msg("Initializing system using INIT_HK_TEST.")
       call lbe_init_hk_test(lbe_N)
       
    case default
       write(msgstr,"('FATAL ERROR: Unknown init_cond <',I0,'>. Aborting...')") init_cond
       call log_msg(msgstr)
       call Abend
    end select
    
    ! Perhaps do some additional things, which are not to be done for restores.
    
    ! Add a perturbation if required.
    if (perturbation .ne. 0.0) then
       write(msgstr,"('Perturbing by fraction ',F16.10)") perturbation
       call log_msg(msgstr)
       call lbe_init_perturb(lbe_N)
    end if
    
    ! If shearing then remove any net x and z momentum
    if ( (boundary_cond .ne. 0) .and. &
         (inv_fluid .eq. 5 .or. inv_fluid .eq. 6) ) then
       call lbe_init_remove_momentum(lbe_N)
    end if

    ! Now initialize the rock if necessary.
    call lbe_init_rock(whole_N)
#ifdef VARTAU
    call halo_exchange(whole_N,vartau_halo)
#endif
#ifdef LOCALBC
    call halo_exchange(whole_N,acccoef_halo)
#endif  
 end if
 
 call log_msg("Successfully initialized system!")
end subroutine lbe_init_system

!>Restores the entire lattice state from one previously saved.
subroutine lbe_init_restore(N)
  implicit none
  type(lbe_site), dimension(0:,0:,0:), intent(inout) :: N

  call restore_checkpoint(N)
end subroutine lbe_init_restore

!>Builds a new lattice state by cloning copies of a smaller
!>saved system. 
subroutine lbe_init_upscale(N)
  implicit none
  type(lbe_site),dimension(0:,0:,0:) :: N

#ifdef USEXDRF
  call restore_upscale(n_restore, chk_uid, N)
#else
  call log_msg("XDRF is disabled, can't perform INIT_UPSCALE")
#endif
end subroutine lbe_init_upscale

!> Initializes fluid velocity in z-direction.
!> Supported initial conditions:
!> -4: constant velocitz
!> -5: Poiseuille profile (calculated from simulation parameters)

subroutine lbe_init_velz(N)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N
  real(kind=rk), dimension(3) :: vel_r, vel_b, vel_g ! initial velocities
  real(kind=rk) :: x_gl ! global x-position
  integer :: channel_diameter ! diameter of the channel (without rock nodes)
  real*8 :: theta, phi, r
  integer :: x, y, z

  do x = 1,nx
    do y = 1,ny
      do z = 1,nz
        ! Set velocity.
        if (init_cond == -6) then ! Poiseuille profile
          channel_diameter = tnx - 2 * boundary_width ! subtract rocks at both sides
          x_gl = x + ccoords(1) * nx - boundary_width - 0.5d0 - 0.5d0 * channel_diameter ! x_gl shall be zero at the channel center
          vel_r(1) = 0.0d0
          vel_r(2) = 0.0d0
          vel_r(3) = vel_poiseuille_max * (1.0d0 - (2.d0 * x_gl / channel_diameter)**2)
          vel_b(1) = 0.0d0
          vel_b(2) = 0.0d0
          vel_b(3) = vel_r(3)
          vel_g(1) = 0.0d0
          vel_g(2) = 0.0d0
          vel_g(3) = vel_r(3)
        else ! init_cond = -5 => constant velocity
          vel_r(1) = 0.0d0
          vel_r(2) = 0.0d0
          vel_r(3) = pr
          vel_b(1) = 0.0d0
          vel_b(2) = 0.0d0
          vel_b(3) = pb
          vel_g(1) = 0.0d0
          vel_g(2) = 0.0d0
          vel_g(3) = pg
        end if

        if (collisiontype_id .eq. MRT_id) then
          CALL mrt_init_dist((/vel_r(1), vel_r(2), vel_r(3)/), fr, N(x,y,z)%n_r)
        else
          CALL boltz_dist(vel_r(1), vel_r(2), vel_r(3), 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, N(x,y,z)%n_r(:))
          N(x,y,z)%n_r(:) = N(x,y,z)%n_r(:) * fr
        end if
#ifdef INTERPOLATEDBB
         ! calculating equillibrium distriubtion functions for interpolated
         ! bounce back scheme requires initial density information. The initial
         ! density of the system is stored to rho_0
         N(x,y,z)%rho_0 = sum(N(x,y,z)%n_r(:)*g(:))
#endif

#ifndef SINGLEFLUID
        if (collisiontype_id .eq. MRT_id) then
          CALL mrt_init_dist((/vel_b(1), vel_b(2), vel_b(3)/), fb, N(x,y,z)%n_b)
        else
          CALL boltz_dist(vel_b(1), vel_b(2), vel_b(3), 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, N(x,y,z)%n_b(:))
          N(x,y,z)%n_b(:) = N(x,y,z)%n_b(:) * fb
        end if
#endif
#ifndef NOSURFACTANT
        if (collisiontype_id .eq. MRT_id) then
          CALL mrt_init_dist((/vel_g(1), vel_g(2), vel_g(3)/), fg, N(x,y,z)%n_s)
        else
          CALL boltz_dist(vel_g(1), vel_g(2), vel_g(3), 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, N(x,y,z)%n_s(:))
          N(x,y,z)%n_s(:) = N(x,y,z)%n_s(:) * fg
        end if

        if (fd .eq. 0) then
          CALL random_number(theta)
          CALL random_number(phi)
          CALL random_number(r)
          theta = theta * pi ! 0 <= theta <= pi
          phi = phi * 2.0d0 * pi ! 0 <=  phi  <= 2pi
          N(x,y,z)%d(1) = r * sin(theta) * cos(phi)
          N(x,y,z)%d(2) = r * sin(theta) * sin(phi)
          N(x,y,z)%d(3) = r * cos(theta)
        else
          N(x,y,z)%d(:) = (/0.0d0, 1.0d0, 0.0d0/)
        end if
#endif
      end do
    end do
  end do
end subroutine lbe_init_velz

subroutine lbe_init_simple(N)
  type(lbe_site), dimension(0:,0:,0:) :: N
  real*8 :: theta, phi, r, F(19)
  integer :: x, y, z

  ! Now initialize my subdomain.

  if (collisiontype_id .eq. MRT_id) then
    CALL mrt_init_dist((/0.0_8,0._8,0.0_8/),1.0d0,F(:))
  else
    call boltz_dist(0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,F)
  end if


  do x = 1, nx
    do y = 1, ny
      do z = 1, nz
        if ( inv_fluid .eq. 11 ) then
          ! calculating a lineal interpolation of fr and pr for
          ! the initial condition for the density
          N(x,y,z)%n_r(:) = F(:)*( (pr - fr)*( start(3) + z - 2 )/(tnz - 1) + fr )
        else
          N(x,y,z)%n_r(:) = F(:)*fr
        end if
#ifdef INTERPOLATEDBB
         ! calculating equillibrium distriubtion functions for interpolated
         ! bounce back scheme requires initial density information. The initial
         ! density of the system is stored to rho_0
         N(x,y,z)%rho_0 = sum(N(x,y,z)%n_r(:)*g(:))
    !     if(x.eq.6.and.y.eq.6.and.z.eq.6) then
    !       write(*,"('pos = ',3I3,' dist  = ',19F16.10)") x,y,z,N(x,y,z)%n_r(:)
    !     end if


#endif
#ifndef SINGLEFLUID
        N(x,y,z)%n_b(:) = F(:)*fb
#endif
#ifndef NOSURFACTANT
        N(x,y,z)%n_s(:) = F(:)*fg

        if ( fd .eq. 0 ) then
          CALL random_number(theta)
          CALL random_number(phi)
          CALL random_number(r)
          theta = theta*pi ! 0 <= theta <= pi
          phi = phi*2.d0*pi ! 0 <=  phi  <= 2pi
          N(x,y,z)%d(1) = r*sin(theta)*cos(phi)
          N(x,y,z)%d(2) = r*sin(theta)*sin(phi)
          N(x,y,z)%d(3) = r*cos(theta)
        else
          N(x,y,z)%d(:) = (/ 0., 1., 0. /)
        end if
#endif
      end do
    end do
  end do
end subroutine lbe_init_simple

subroutine lbe_init_ratio(N)
  implicit none
  type(lbe_site),dimension(0:,0:,0:) :: N
  real*8 :: theta, phi, r, F(19)
  integer :: x, y, z
  
  if (collisiontype_id .eq. MRT_id) then
     CALL mrt_init_dist((/0.0_8,0._8,0.0_8/),1.0d0,F(:))
  else
     CALL boltz_dist(0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,F(:))
  end if

  do x = 1, nx
     do y = 1, ny
        do z = 1, nz
           N(x,y,z)%n_r(:) = F(:)*fr
#ifndef SINGLEFLUID
           N(x,y,z)%n_b(:) = F(:)*fb
#endif
#ifndef NOSURFACTANT
           N(x,y,z)%n_s(:) = F(:)*fg

           CALL random_number(theta)
           CALL random_number(phi)
           CALL random_number(r)
           theta = theta*pi ! 0 <= theta <= pi
           phi = phi * 2.d0 * pi ! 0 <=  phi  <= 2pi
           N(x,y,z)%d(1) = r*sin(theta)*cos(phi)
           N(x,y,z)%d(2) = r*sin(theta)*sin(phi)
           N(x,y,z)%d(3) = r*cos(theta)
#endif
        end do
     end do
  end do
end subroutine lbe_init_ratio

subroutine lbe_init_ratio_plusdisturb(N)
  implicit none
  type(lbe_site),dimension(0:,0:,0:) :: N
  real*8 :: theta, phi, r, F(19)
  integer :: x, y, z
  
  if (collisiontype_id .eq. MRT_id) then
     CALL mrt_init_dist((/0.0_8,0._8,0.0_8/),1.0d0,F(:))
  else
     CALL boltz_dist(0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,F(:))
  end if

  do x = 1, nx
     do y = 1, ny
        do z = 1, nz
           N(x,y,z)%n_r(:) = F(:)*fr
           if(x .eq. 20) then
               if(y .eq. 20) then
                   if(z .eq. 20) then
                       N(x,y,z)%n_r(:) = N(x,y,z)%n_r(:)*1.2
                   end if
               end if
           end if
           if(x .eq. 20) then
               if(y .eq. 16) then
                   if(z .eq. 8) then
                       N(x,y,z)%n_r(:) = N(x,y,z)%n_r(:)*1.2
                   end if
               end if
           end if
           N(x,y,z)%n_r(:) = F(:)*fr
           if(x .eq. 17) then
               if(y .eq. 17) then
                   if(z .eq. 17) then
                       N(x,y,z)%n_r(:) = N(x,y,z)%n_r(:)*0.8
                   end if
               end if
           end if
           if(x .eq. 17) then
               if(y .eq. 13) then
                   if(z .eq. 9) then
                       N(x,y,z)%n_r(:) = N(x,y,z)%n_r(:)*0.8
                   end if
               end if
           end if
#ifndef SINGLEFLUID
           N(x,y,z)%n_b(:) = F(:)*fb
           if(x .eq. 20) then
               if(y .eq. 20) then
                   if(z .eq. 20) then
                       N(x,y,z)%n_b(:) = N(x,y,z)%n_b(:)*1.2
                   end if
               end if
           end if
           if(x .eq. 20) then
               if(y .eq. 16) then
                   if(z .eq. 8) then
                       N(x,y,z)%n_b(:) = N(x,y,z)%n_b(:)*1.2
                   end if
               end if
           end if
           N(x,y,z)%n_r(:) = F(:)*fr
           if(x .eq. 17) then
               if(y .eq. 17) then
                   if(z .eq. 17) then
                       N(x,y,z)%n_b(:) = N(x,y,z)%n_b(:)*0.8
                   end if
               end if
           end if
           if(x .eq. 17) then
               if(y .eq. 13) then
                   if(z .eq. 9) then
                       N(x,y,z)%n_b(:) = N(x,y,z)%n_b(:)*0.8
                   end if
               end if
           end if
#endif
        end do
     end do
  end do
end subroutine lbe_init_ratio_plusdisturb

subroutine lbe_init_sinus(N)
  implicit none
  type(lbe_site),dimension(0:,0:,0:) :: N
  real*8 :: theta, phi, r, F(19)
  integer :: x, y, z
  
  if (collisiontype_id .eq. MRT_id) then
     CALL mrt_init_dist((/0.0_8,0._8,0.0_8/),1.0d0,F(:))
  else
     CALL boltz_dist(0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,F(:))
  end if

  do x = 1, nx
     do y = 1, ny
        do z = 1, nz
           N(x,y,z)%n_r(:) = F(:)*fr*(1+fr2*sin(2.0_rk*pi*x/fr1)*sin(2.0_rk*pi*y/fr1)*sin(2.0_rk*pi*z/fr1))
#ifndef SINGLEFLUID
           N(x,y,z)%n_b(:) = F(:)*fb*(1-fr2*sin(2.0_rk*pi*x/fr1)*sin(2.0_rk*pi*y/fr1)*sin(2.0_rk*pi*z/fr1))
#endif
        end do
     end do
  end do
end subroutine lbe_init_sinus


subroutine lbe_init_rand(N)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N
  real*8 :: theta, phi, r
  integer :: x, y, z
  real*8, dimension(nvecs) :: g_inv
  g_inv = 1./g

  ! First, read the LBE parameters from input-file.

  call lbe_get_lbe_input()
  ! Now initialize my subdomain.

  call log_msg("Initializing subdomain.")
  do x = 1,nx
    do y = 1,ny
      do z = 1,nz
        CALL random_number(N(x,y,z)%n_r(:))
        N(x,y,z)%n_r(:) = fr*N(x,y,z)%n_r(:)*g_inv/nvecs
#ifndef SINGLEFLUID
        CALL random_number(N(x,y,z)%n_b(:))
        N(x,y,z)%n_b(:) = fb*N(x,y,z)%n_b(:)*g_inv/nvecs
#endif

#ifndef NOSURFACTANT
        CALL random_number(N(x,y,z)%n_s(:))
        N(x,y,z)%n_s(:) = fg*N(x,y,z)%n_s(:)*g_inv/nvecs

        CALL random_number(theta)
        CALL random_number(phi)
        CALL random_number(r)	! 0 <=   r   <= 1
        theta = theta*pi	! 0 <= theta <= pi
        phi = phi*2.d0*pi	! 0 <=  phi  <= 2pi
        N(x,y,z)%d(1) = r*sin(theta)*cos(phi)
        N(x,y,z)%d(2) = r*sin(theta)*sin(phi)
        N(x,y,z)%d(3) = r*cos(theta)
#endif
      end do
    end do
  end do
end subroutine lbe_init_rand

subroutine lbe_init_frac(N)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N
  real*8 :: theta, phi, r, s, F(19)
  integer :: x, y, z
  real*8, dimension(nvecs) :: g_inv

  g_inv = 1.d0/g

  ! First, read the LBE parameters from input-file.
  CALL lbe_get_lbe_input()

  if (collisiontype_id .eq. MRT_id) then
    CALL mrt_init_dist((/0.0_8,0._8,0.0_8/),1.0d0,F(:))
  else
    CALL boltz_dist(0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,F(:))
  end if

  do x = 1,nx
    do y = 1,ny
      do z = 1,nz
        CALL random_number(N(x,y,z)%n_r(:))
        N(x,y,z)%n_r(:) = N(x,y,z)%n_r(:)/sum( N(x,y,z)%n_r(:)*g(:) )
        CALL random_number(s)
        N(x,y,z)%n_r(:) = fr*s*N(x,y,z)%n_r(:)
#ifndef SINGLEFLUID
        CALL random_number(N(x,y,z)%n_b(:))
        N(x,y,z)%n_b(:) = N(x,y,z)%n_b(:)/sum( N(x,y,z)%n_b(:)*g(:) )
        CALL random_number(s)
        N(x,y,z)%n_b(:) = fb*s*N(x,y,z)%n_b(:)
#ifndef NOSURFACTANT
        CALL random_number(N(x,y,z)%n_s(:))
        N(x,y,z)%n_s(:) = N(x,y,z)%n_s(:)/sum( N(x,y,z)%n_s(:)*g(:) )
        CALL random_number(s)
        N(x,y,z)%n_s(:) = fg*s*N(x,y,z)%n_s(:)

        CALL random_number(theta)
        CALL random_number(phi)
        CALL random_number(r)
        theta = theta*pi ! 0 <= theta <= pi
        phi = phi*2.d0*pi ! 0 <=  phi  <= 2pi
        N(x,y,z)%d(1) = r*sin(theta)*cos(phi)
        N(x,y,z)%d(2) = r*sin(theta)*sin(phi)
        N(x,y,z)%d(3) = r*cos(theta)
#endif
#endif
      end do
    end do
  end do
end subroutine lbe_init_frac

subroutine lbe_init_frac_fromv5(N)
  implicit none
  type(lbe_site), intent(inout), dimension(0:,0:,0:) :: N
  real(kind=rk) :: theta, phi, r
  integer :: x, y, z, s
  real(kind=rk), dimension(nvecs) :: g_inv
  integer :: irndcount
  real(kind=rk), dimension(:), allocatable :: rndtmp

  g_inv = 1.0_rk/g

  ! Generate all necessary random numbers
  allocate(rndtmp(nx*ny*nz*6*nvecs))
  call ranlux(rndtmp,nx*ny*nz*6*nvecs)
  irndcount = 1

  do x=1,nx
    do y=1,ny
      do z=1,nz
        do s=1,nvecs
          N(x,y,z)%n_r(s) = rndtmp(irndcount)
#ifndef SINGLEFLUID
          N(x,y,z)%n_b(s) = rndtmp(irndcount+1)
#endif
#ifndef NOSURFACTANT
          N(x,y,z)%n_s(s) = rndtmp(irndcount+2)
#endif
          irndcount = irndcount + 3
        end do
#ifndef NOSURFACTANT
        theta = rndtmp(irndcount)*pi      ! 0<=theta<=pi
        phi   = rndtmp(irndcount+1)*2.0_rk*pi ! 0<=phi<= 2pi
        r     = rndtmp(irndcount+2)       ! 0<= r <= 1
        irndcount = irndcount + 3

        N(x,y,z)%d(1) = r*sin(theta)*cos(phi)
        N(x,y,z)%d(2) = r*sin(theta)*sin(phi)
        N(x,y,z)%d(3) = r*cos(theta)

        ! Now distribute things.
        N(x,y,z)%n_s(:) = N(x,y,z)%n_s(:)*g_inv*fg/nvecs
#endif
        N(x,y,z)%n_r(:) = N(x,y,z)%n_r(:)*g_inv*fr/nvecs
#ifndef SINGLEFLUID
        N(x,y,z)%n_b(:) = N(x,y,z)%n_b(:)*g_inv*fb/nvecs
#endif
      end do
    end do
  end do
  deallocate(rndtmp)
end subroutine lbe_init_frac_fromv5

subroutine lbe_init_frac_determ(N)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N
  real*8 :: theta, phi, r, s, F(19)
  integer :: x, y, z
  real*8, dimension(nvecs) :: g_inv

  g_inv = 1.d0/g

  ! First, read the LBE parameters from input-file.
  CALL lbe_get_lbe_input()

  if (collisiontype_id .eq. MRT_id) then
    CALL mrt_init_dist((/0.0_8,0._8,0.0_8/),1.0d0,F(:))
  else
    CALL boltz_dist(0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,F(:))
  end if

  do x = 1,nx
    do y = 1,ny
      do z = 1,nz
        CALL random_number(N(x,y,z)%n_r(:))
        N(x,y,z)%n_r(:) = N(x,y,z)%n_r(:)/sum( N(x,y,z)%n_r(:)*g(:) )
        CALL random_number(s)
        N(x,y,z)%n_r(:) = fr*s*N(x,y,z)%n_r(:)
#ifndef SINGLEFLUID
        CALL random_number(N(x,y,z)%n_b(:))
        N(x,y,z)%n_b(:) = N(x,y,z)%n_b(:)/sum( N(x,y,z)%n_b(:)*g(:) )
        CALL random_number(s)
        N(x,y,z)%n_b(:) = fb*s*N(x,y,z)%n_b(:)
#endif
      end do
    end do
  end do
  do x = 1, nx
     do y = 1, ny
        do z = 1, nz
           N(x,y,z)%n_r(:) = F(:)*fr
#ifndef SINGLEFLUID
           N(x,y,z)%n_b(:) = F(:)*fb
#endif
        end do
     end do
  end do

end subroutine lbe_init_frac_determ

!> The mean x and z momentum are forced to be 0 - for Lees Edwards code to work.
subroutine lbe_init_remove_momentum(N)
  type(lbe_site), dimension(0:,0:,0:) :: N
  integer :: s
  real*8 :: sumposx,sumnegx,sumposz,sumnegz

  sumposx = 0.0_rk
  sumnegx = 0.0_rk
  sumposz = 0.0_rk
  sumnegz = 0.0_rk

  ! For Lees-Edwards x-momentum must sum to zero
  ! This is a mess - but it only runs once & I have other problems

  do s = 1,nnp
    sumnegx = sumnegx + sum(N%n_r(negx(s))*g(negx(s)))
    sumposx = sumposx + sum(N%n_r(posx(s))*g(posx(s)))
  enddo
  ! Clip to avoid division by zero.
  sumposx = max(dble(1e-9),sumposx)
  do s = 1,nnp
    N%n_r(posx(s)) = N%n_r(posx(s))*sumnegx/sumposx
  enddo

  sumnegx = 0.d0
  sumposx = 0.d0

#ifndef NOSURFACTANT
  do s = 1,nnp
    sumnegx = sumnegx + sum(N%n_s(negx(s))*g(negx(s)))
    sumposx = sumposx + sum(N%n_s(posx(s))*g(posx(s)))
  enddo
  sumposx = max(dble(1e-9),sumposx)
  do s = 1,nnp
    N%n_s(posx(s)) = N%n_s(posx(s))*sumnegx/sumposx
  enddo
#endif

#ifndef SINGLEFLUID
  sumnegx = 0.d0
  sumposx = 0.d0

  do s=1,nnp
    sumnegx = sumnegx + sum(N%n_b(negx(s))*g(negx(s)))
    sumposx = sumposx + sum(N%n_b(posx(s))*g(posx(s)))
  enddo
  sumposx = max(dble(1e-9),sumposx)
  do s = 1,nnp
    N%n_b(posx(s)) = N%n_b(posx(s))*sumnegx/sumposx
  enddo
#endif
  ! It is nice if mean z-momentum is zero as then graph
  ! intercepts at +/- shear_u

  sumnegz = 0.d0
  sumposz = 0.d0

  do s = 1,nnp
    sumnegz = sumnegz + sum(N%n_r(negz(s))*g(negz(s)))
    sumposz = sumposz + sum(N%n_r(posz(s))*g(posz(s)))
  enddo
  sumposz = max(dble(1e-9),sumposz)
  do s = 1,nnp
    N%n_r(posz(s)) = N%n_r(posz(s))*sumnegz/sumposz
  enddo

#ifndef NOSURFACTANT
  sumnegz = 0.d0
  sumposz = 0.d0

  do s = 1,nnp
    sumnegz = sumnegz + sum(N%n_s(negz(s))*g(negz(s)))
    sumposz = sumposz + sum(N%n_s(posz(s))*g(posz(s)))
  enddo
  sumposz = max(dble(1e-9),sumposz)
  do s = 1,nnp
    N%n_s(posz(s)) = N%n_s(posz(s))*sumnegz/sumposz
  enddo
#endif

#ifndef SINGLEFLUID
  sumnegz = 0.d0
  sumposz = 0.d0

  do s=1,nnp
    sumnegz = sumnegz + sum(N%n_b(negz(s))*g(negz(s)))
    sumposz = sumposz + sum(N%n_b(posz(s))*g(posz(s)))
  enddo
  sumposz = max(dble(1e-9),sumposz)
  do s = 1,nnp
    N%n_b(posz(s)) = N%n_b(posz(s))*sumnegz/sumposz
  enddo
#endif
end subroutine lbe_init_remove_momentum

!> \{
!> \name Init subroutines for ternary systems
!> lbe_init_lamx
!> lbe_init_lamy
!> lbe_init_lamz
!> lbe_init_radial
!>
!> N.B. Be careful with the precompiler ifdefs if adding/removing
!> subroutines from here.

!> Sets up the system to have lamellae perpendicular to the X-axis.
!> Format is: \c fr1 site widths of oil, then 1 surfactant,
!> then \c fr2 sites of water, then 1 surfactant, repeated.
!> Note that (\c fr1+\c fr2+2) should divide the lattice size evenly.
subroutine lbe_init_lamx(N)
  implicit none
  type(lbe_site),dimension(0:,0:,0:) :: N
  integer :: x, y, z, i
  real*8 :: theta, phi, r, F(19)
  integer :: lamwidth	! Width of one layer.
  integer :: tx	! x-coordinate in global system.

  if (collisiontype_id .eq. MRT_id) then
    CALL mrt_init_dist((/0.0_8,0._8,0.0_8/),1.0d0,F(:))
  else
    CALL boltz_dist(0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,F(:))
  end if

  ! Dipoles point towards water.
  ! FIXME Check that's the right way round.

  lamwidth = fr1 + fr2 + 2

  if (myrankc == 0) then
    if (mod(tnx,lamwidth) .ne. 0) then
      print*,'******* WARNING ******'
      print*,'* Truncated lamellae *'
      print*,'**********************'
    endif
  endif

  ! Note that tnx is zero-based.
  do x=1,nx
    tx = x-1 + ccoords(1)*nx
    i=mod(tx,lamwidth)

    if (i<fr1) then
      ! Make oil happen.
      do z = 1,nz
        do y = 1,ny
          N(x,y,z)%n_r(:) = F(:)*fr
#ifndef SINGLEFLUID
          N(x,y,z)%n_b(:) = F(:)*fb
#endif
#ifndef NOSURFACTANT
          N(x,y,z)%n_s(:) = F(:)*fg

          if ( fd .eq. 0 ) then
            CALL random_number(theta)
            CALL random_number(phi)
            CALL random_number(r)
            theta = theta*pi ! 0 <= theta <= pi
            phi = phi*2.d0*pi ! 0 <=  phi  <= 2pi
            N(x,y,z)%d(1) = r*sin(theta)*cos(phi)
            N(x,y,z)%d(2) = r*sin(theta)*sin(phi)
            N(x,y,z)%d(3) = r*cos(theta)
          else
            N(x,y,z)%d(:) = (/ 0., 1., 0. /)
          end if
#endif
        end do
      end do
    else if ( (i == fr1) .or. (i == (fr1+fr2+1)) ) then
      ! Make surf happen.
      do z = 1,nz
        do y = 1,ny
          N(x,y,z)%n_r(:) = F(:)*qr
#ifndef SINGLEFLUID
          N(x,y,z)%n_b(:) = F(:)*qb
#endif
#ifndef NOSURFACTANT            
          N(x,y,z)%n_s(:) = F(:)*qg

          if ( qd .eq. 0 ) then
            CALL random_number(theta)
            CALL random_number(phi)
            CALL random_number(r)
            theta = theta*pi ! 0 <= theta <= pi
            phi = phi*2.d0*pi ! 0 <=  phi  <= 2pi
            N(x,y,z)%d(1) = r*sin(theta)*cos(phi)
            N(x,y,z)%d(2) = r*sin(theta)*sin(phi)
            N(x,y,z)%d(3) = r*cos(theta)
          else
            N(x,y,z)%d(:) = (/ 0., 1., 0. /)
          end if
#endif
        end do
      end do
    elseif (i < (fr1 + 1 + fr2) ) then
      ! Make water happen.
      do z = 1,nz
        do y = 1,ny
          N(x,y,z)%n_r(:) = F(:)*pr
#ifndef SINGLEFLUID
          N(x,y,z)%n_b(:) = F(:)*pb
#endif
#ifndef NOSURFACTANT
          N(x,y,z)%n_s(:) = F(:)*pg

          if ( pd .eq. 0 ) then
            CALL random_number(theta)
            CALL random_number(phi)
            CALL random_number(r)
            theta = theta*pi ! 0 <= theta <= pi
            phi = phi*2.d0*pi ! 0 <=  phi  <= 2pi
            N(x,y,z)%d(1) = r*sin(theta)*cos(phi)
            N(x,y,z)%d(2) = r*sin(theta)*sin(phi)
            N(x,y,z)%d(3) = r*cos(theta)
          else
            N(x,y,z)%d(:) = (/ 0., 1., 0. /)
          end if
#endif
        end do
      end do
    else
      print*,"Can't happen in lam-x!"
      call Abend
    end if
  end do
end subroutine lbe_init_lamx

!> Sets up the system to have lamellae perpendicular to the Y-axis.
!> Format is: \c fr1 site widths of oil, then 1 surfactant,
!> then \c fr2 sites of water, then 1 surfactant, repeated.
!> Note that (\c fr1+\c fr2+2) should divide the lattice size evenly.
subroutine lbe_init_lamy(N)
  implicit none
  type(lbe_site),dimension(0:,0:,0:) :: N
  integer :: x, y ,z, i
  real*8 :: theta, phi, r, F(19)
  integer :: lamwidth	! Width of one layer.
  integer :: ty	! y-coordinate in global system.

  ! Dipoles point towards water.
  ! FIXME Check that's the right way round.

  lamwidth = fr1 + fr2 + 2

  if (myrankc == 0) then
    if (mod(tny,lamwidth) .ne. 0) then
      print*,'******* WARNING ******'
      print*,'* Truncated lamellae *'
      print*,'**********************'
    endif
  endif

  if (collisiontype_id .eq. MRT_id) then
    CALL mrt_init_dist((/0.0_8,0._8,0.0_8/),1.0d0,F(:))
  else
    CALL boltz_dist(0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,F(:))
  end if

  ! Note that tny is zero-based.
  do y = 1,ny
    ty = y - 1 + ccoords(2)*ny
    i = mod(ty,lamwidth)

    if ( i<fr1 ) then
      ! Make oil happen.
      do x = 1,nx
        do z = 1,nz
          N(x,y,z)%n_r(:) = F(:)*fr
#ifndef SINGLEFLUID
          N(x,y,z)%n_b(:) = F(:)*fb
#endif
#ifndef NOSURFACTANT
          N(x,y,z)%n_s(:) = F(:)*fg

          if ( fd .eq. 0 ) then
            CALL random_number(theta)
            CALL random_number(phi)
            CALL random_number(r)
            theta = theta*pi ! 0 <= theta <= pi
            phi = phi*2.d0*pi ! 0 <=  phi  <= 2pi
            N(x,y,z)%d(1) = r*sin(theta)*cos(phi)
            N(x,y,z)%d(2) = r*sin(theta)*sin(phi)
            N(x,y,z)%d(3) = r*cos(theta)
          else
            N(x,y,z)%d(:) = (/ 0., 1., 0. /)
          end if
#endif
        end do
      end do
    elseif ( (i == fr1) .or. (i==(fr1+fr2+1)) ) then
      ! Make surf happen.
      do x = 1,nx
        do z = 1,nz
          N(x,y,z)%n_r(:) = F(:)*qr
#ifndef SINGLEFLUID
          N(x,y,z)%n_b(:) = F(:)*qb
#endif
#ifndef NOSURFACTANT
          N(x,y,z)%n_s(:) = F(:)*qg

          if ( qd .eq. 0 ) then
            CALL random_number(theta)
            CALL random_number(phi)
            CALL random_number(r)
            theta = theta*pi ! 0 <= theta <= pi
            phi = phi*2.d0*pi ! 0 <=  phi  <= 2pi
            N(x,y,z)%d(1) = r*sin(theta)*cos(phi)
            N(x,y,z)%d(2) = r*sin(theta)*sin(phi)
            N(x,y,z)%d(3) = r*cos(theta)
          else
            N(x,y,z)%d(:) = (/ 0., 1., 0. /)
          end if
#endif
        end do
      end do
    elseif ( i < (fr1 + 1 + fr2) ) then
      ! Make water happen.
      do x = 1,nx
        do z = 1,nz
          N(x,y,z)%n_r(:) = F(:)*pr
#ifndef SINGLEFLUID
          N(x,y,z)%n_b(:) = F(:)*pb
#endif
#ifndef NOSURFACTANT
          N(x,y,z)%n_s(:) = F(:)*pg

          if ( pd .eq. 0 ) then
            CALL random_number(theta)
            CALL random_number(phi)
            CALL random_number(r)
            theta = theta*pi ! 0 <= theta <= pi
            phi = phi*2.d0*pi ! 0 <=  phi  <= 2pi
            N(x,y,z)%d(1) = r*sin(theta)*cos(phi)
            N(x,y,z)%d(2) = r*sin(theta)*sin(phi)
            N(x,y,z)%d(3) = r*cos(theta)
          else
            N(x,y,z)%d(:) = (/ 0., 1., 0. /)
          end if
#endif
        end do
      end do
    else
      print*,"Can't happen in lam-y!"
      call Abend
    end if
  end do
end subroutine lbe_init_lamy

!> 
!> Sets up the system to have lamellae perpendicular to the Y-axis.
!> Format is: \c fr1 site widths of oil, then 1 surfactant,
!> then \c fr2 sites of water, then 1 surfactant, repeated.
!> Note that (\c fr1+\c fr2+2) should divide the lattice size evenly.
subroutine lbe_init_lamz(N)
  implicit none
  type(lbe_site),dimension(0:,0:,0:) :: N
  integer :: x, y, z, i
  real*8 :: theta, phi, r, F(19)
  integer :: lamwidth	! Width of one layer.
  integer :: tz	! z-coordinate in global system.

  ! Dipoles point towards water.
  ! FIXME Check that's the right way round.

  lamwidth = fr1 + fr2 + 2

  if (myrankc == 0) then
    if (mod(tnz,lamwidth) .ne. 0) then
      print*,'******* WARNING ******'
      print*,'* Truncated lamellae *'
      print*,'**********************'
    endif
  endif

  if (collisiontype_id .eq. MRT_id) then
    CALL mrt_init_dist((/0.0_8,0._8,0.0_8/),1.0d0,F(:))
  else
    CALL boltz_dist(0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,F(:))
  end if

  ! Note that tnz is zero-based.
  do z=1,nz
    tz = z-1 + ccoords(3)*nz
    i=mod(tz,lamwidth)
    print*,'z=',z,'tz=',tz,'i=',i

    if ( i<fr1 ) then
      ! Make oil happen.
      print*,' oil'
      do x = 1,nx
        do y = 1,ny
          N(x,y,z)%n_r(:) = F(:)*fr
#ifndef SINGLEFLUID
          N(x,y,z)%n_b(:) = F(:)*fb
#endif
#ifndef NOSURFACTANT
          N(x,y,z)%n_s(:) = F(:)*fg

          if ( fd .eq. 0 ) then
            CALL random_number(theta)
            CALL random_number(phi)
            CALL random_number(r)
            theta = theta*pi ! 0 <= theta <= pi
            phi = phi*2.d0*pi ! 0 <=  phi  <= 2pi
            N(x,y,z)%d(1) = r*sin(theta)*cos(phi)
            N(x,y,z)%d(2) = r*sin(theta)*sin(phi)
            N(x,y,z)%d(3) = r*cos(theta)
          else
            N(x,y,z)%d(:) = (/ 0., 1., 0. /)
          end if
#endif
        end do
      end do
    elseif ((i == fr1) .or. (i==(fr1+fr2+1))) then
      ! Make surf happen.
      print*,' surf'
      do x = 1,nx
        do y = 1,ny
          N(x,y,z)%n_r(:) = F(:)*qr
#ifndef SINGLEFLUID
          N(x,y,z)%n_b(:) = F(:)*qb
#endif
#ifndef NOSURFACTANT
          N(x,y,z)%n_s(:) = F(:)*qg

          if ( qd .eq. 0 ) then
            CALL random_number(theta)
            CALL random_number(phi)
            CALL random_number(r)
            theta = theta*pi ! 0 <= theta <= pi
            phi = phi*2.d0*pi ! 0 <=  phi  <= 2pi
            N(x,y,z)%d(1) = r*sin(theta)*cos(phi)
            N(x,y,z)%d(2) = r*sin(theta)*sin(phi)
            N(x,y,z)%d(3) = r*cos(theta)
          else
            N(x,y,z)%d(:) = (/ 0., 1., 0. /)
          end if
#endif
        end do
      end do
    elseif (i < (fr1 + 1 + fr2) ) then
      ! Make water happen.
      print*,' water'
      do x=1,nx
        do y=1,ny
          N(x,y,z)%n_r(:) = F(:)*pr
#ifndef SINGLEFLUID
          N(x,y,z)%n_b(:) = F(:)*pb
#endif
#ifndef NOSURFACTANT
          N(x,y,z)%n_s(:) = F(:)*pg

          if ( pd .eq. 0 ) then
            CALL random_number(theta)
            CALL random_number(phi)
            CALL random_number(r)
            theta = theta*pi ! 0 <= theta <= pi
            phi = phi*2.d0*pi ! 0 <=  phi  <= 2pi
            N(x,y,z)%d(1) = r*sin(theta)*cos(phi)
            N(x,y,z)%d(2) = r*sin(theta)*sin(phi)
            N(x,y,z)%d(3) = r*cos(theta)
          else
            N(x,y,z)%d(:) = (/ 0., 1., 0. /)
          end if
#endif
        end do
      end do
    else
      print*,"Can't happen in lam-z!"
      call Abend
    end if
  end do
end subroutine lbe_init_lamz

!> This is a general routine for initialising droplike or vesicular
!> initial systems.  \c model contains an \c lbe_site type for each of
!> the three regions xE<lt>fr1, fr1E<lt>xE<gt>fr2, and xE<gt>fr2.
!>
!> Each site is initialised to the appropriate value.
!>
!>However, the dipole moment is ignored.
!>
!>The dipole moment of each site is set to have unit value, and
!>its orientation depends on the rock_state value for the region.
!>If \c rock_state is negative, the dipole points towards the centre of
!>the system (ie -r); if positive, it points away from the centre.
!>If zero, its orientation is chosen randomly.
!>
!>The ME3D initial conditions are all special cases of this model.
!>
!>Currently, the values are averaged over the eight points of the unit
!>cube. This sucks.
subroutine lbe_init_radial(N,model)
  type(lbe_site), dimension(0:,0:,0:) :: N
  integer :: x, y, z, s
  integer, dimension(3) :: centre     ! Global coords of centre of system.
  integer, dimension(3) :: base	      ! Global coords of corner of subdomain
  integer, dimension(3) :: offset
  integer, dimension(8,3) :: vertices ! Vertices of a unit cube
  real*8 :: if1, if2	              ! Real radii

  real*8, dimension(3) :: r	! Global coords of a given point
  real*8 :: rad			! Radius of point
  integer :: i		! Loop variable over vertices

  real*8 :: n_r
#ifndef SINGLEFLUID
  real*8 :: n_b
#endif
#ifndef NOSURFACTANT
  real*8 :: n_s
  real*8, dimension(3) :: d, dd
#endif
  real*8, parameter :: A = 1./8.	
  real*8 :: F(19)
  type(radial_model) :: model(3)  
  integer :: rad_index
#ifndef NOSURFACTANT
  real*8 :: theta, phi
#endif

  if (collisiontype_id .eq. MRT_id) then
    CALL mrt_init_dist((/0.0_8,0._8,0.0_8/),1.0d0,F(:))
  else
    CALL boltz_dist(0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,F(:))
  end if

  ! All coordinates in this routine are in the global
  ! coordinate system.

  ! Define the actual radii of the spheres.
  if1 = fr1*min(tnx, tny, tnz)/2.d0
  if2 = fr2*min(tnx, tny, tnz)/2.d0

  !Define the centre of the system.
  centre = (/ tnx, tny, tnz /)/2 + (/ drop_xshift, drop_yshift, drop_zshift /)

  ! Define the minimum xyz corner of my subdomain.
  base = ccoords*(/ nx, ny, nz /)
  ! Subtract the two.
  offset = base - centre

  ! Define the vectors of a unit cube.
  vertices(1,:) = (/ 0, 0, 0 /)
  vertices(2,:) = (/ 0, 0, 1 /)
  vertices(3,:) = (/ 0, 1, 0 /)
  vertices(4,:) = (/ 0, 1, 1 /)
  vertices(5,:) = (/ 1, 0, 0 /)
  vertices(6,:) = (/ 1, 0, 1 /)
  vertices(7,:) = (/ 1, 1, 0 /)
  vertices(8,:) = (/ 1, 1, 1 /)

  write(msgstr,"('if1 = ',F16.10,' if2 = ',F16.10)") if1, if2
  call log_msg(msgstr)

  ! Note that r is zero-based.

  do x = 1,nx
    do y = 1,ny
      do z = 1,nz
        n_r = 0.d0
#ifndef SINGLEFLUID
        n_b = 0.d0
#endif
#ifndef NOSURFACTANT
        n_s = 0.d0
        d = 0.d0
#endif
        do i = 1,8 ! Average over each corner of the site cube.
          r = (/ x - 1, y - 1, z - 1 /) + offset + vertices(i,:) 
          rad = sqrt(sum(r*r))

          ! Find the appropriate index.
          if ( rad < if1 ) then
            rad_index = rad_inner
          elseif ( (rad >= if1) .and. (rad <= if2) ) then
            rad_index = rad_middle
          else
            rad_index = rad_outer
          endif

          ! Set appropriate values.
          n_r = n_r + model(rad_index)%n_r
#ifndef SINGLEFLUID
          n_b = n_b + model(rad_index)%n_b
#endif
#ifndef NOSURFACTANT
          n_s = n_s + model(rad_index)%n_s

          ! Handle the dipole moment.
          ! (but only if there are actually any dipoles..)              
          if ( model(rad_index)%n_s > 0.d0 ) then
            if ( model(rad_index)%dip < 0.d0 ) then                  
              ! Dipole points towards centre.                
              d = d - r/rad                    
            elseif ( model(rad_index)%dip == 0) then
              ! Make a random dipole of unit size.
              call random_number(theta)
              call random_number(phi)
              theta = theta*pi	! 0 <= theta <= pi
              phi = phi*2.0d0*pi	! 0 <= phi <= 2pi
              dd(1) = sin(theta)*cos(phi)
              dd(2) = sin(theta)*sin(phi)
              dd(3) = cos(theta)

              ! Add it in to the averaging process.
              d = d + dd
            else
              ! Dipole points away from centre
              d = d + r/rad
            endif
          endif
#endif
        end do

        ! n_r, n_s, n_b are the required average particle densities for one site.

        N(x,y,z)%n_r = A*n_r*F(:)
#ifndef SINGLEFLUID
        N(x,y,z)%n_b = A*n_b*F(:)
#endif
#ifndef NOSURFACTANT
        N(x,y,z)%n_s = A*n_s*F(:)
        N(x,y,z)%d = A*d
#endif
        ! N(x,y,z)%rock_state=rad_index ! Useful for debugging.
      end do
    end do
  end do
end subroutine lbe_init_radial

subroutine lbe_init_cylinder_x(N)
  implicit none

  type(lbe_site),dimension(0:,0:,0:),intent(inout) :: N

  integer, dimension(8,3) :: vertices ! Vertices of a unit cube
  real(kind=rk) :: F(19)
  real(kind=rk) :: r_c2, r2, tr, tb

  integer, dimension(2) :: centre , base
  integer :: x,y,z,ty,tz,i

#ifndef NOSURFACTANT
  call log_msg("WARNING: cylinder initialization does not support surfactant.")
#endif

  if (collisiontype_id .eq. MRT_id) then
    call mrt_init_dist((/0.0_rk,0._rk,0.0_rk/),1.0_rk,F(:))
  else
    call boltz_dist(0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,F(:))
  end if

  r_c2 = ( fr1*min(tny, tnz)/2.0_rk )**2

  write(msgstr,"('r_c2 = ',F16.10)") r_c2
  call log_msg(msgstr)
  
  ! Define the vectors of a unit cube.
  vertices(1,:) = (/ 0, 0, 0 /)
  vertices(2,:) = (/ 0, 0, 1 /)
  vertices(3,:) = (/ 0, 1, 0 /)
  vertices(4,:) = (/ 0, 1, 1 /)
  vertices(5,:) = (/ 1, 0, 0 /)
  vertices(6,:) = (/ 1, 0, 1 /)
  vertices(7,:) = (/ 1, 1, 0 /)
  vertices(8,:) = (/ 1, 1, 1 /)

  !Define the centre of the system.
  ! centre = (/ tny, tnz /)/2 
  centre = (/ tny, tnz /)/2 + (/ drop_yshift, drop_zshift /)

  ! Define the minimum xyz corner of my subdomain.
  base = ccoords(2:3) * (/ ny, nz /)

  do x = 1,nx
    do y = 1,ny
      do z = 1,nz
        tr = 0.0_rk
        tb = 0.0_rk
        do i = 1,8
          ty = y + base(1) + vertices(i,2)
          tz = z + base(2) + vertices(i,3)
          r2 = ( ty - centre(1) )**2 + ( tz - centre(2) )**2
          if ( r2 < r_c2 ) then
            tr = tr + fr
#ifndef SINGLEFLUID
            tb = tb + fb
#endif
          else
            tr = tr + pr
#ifndef SINGLEFLUID
            tb = tb + pb
#endif
          end if
        end do
        N(x,y,z)%n_r = tr*F(:)/8.0_rk
#ifndef SINGLEFLUID
        N(x,y,z)%n_b = tb*F(:)/8.0_rk
#endif
      end do
    end do
  end do
end subroutine lbe_init_cylinder_x

subroutine lbe_init_cylinder_z(N)
  implicit none

  type(lbe_site),dimension(0:,0:,0:),intent(inout) :: N

  integer, dimension(8,3) :: vertices ! Vertices of a unit cube
  real(kind=rk) :: F(19)
  real(kind=rk) :: r_c2, r2, tr, tb

  integer, dimension(2) :: centre , base
  integer :: x,y,z,tx,ty,i

#ifndef NOSURFACTANT
  call log_msg("WARNING: cylinder initialization does not support surfactant.")
#endif

  if (collisiontype_id .eq. MRT_id) then
    call mrt_init_dist((/0.0_rk,0._rk,0.0_rk/),1.0_rk,F(:))
  else
    call boltz_dist(0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,F(:))
  end if

  r_c2 = ( fr1*min(tnx, tny)/2.0_rk )**2

  write(msgstr,"('r_c2 = ',F16.10)") r_c2
  call log_msg(msgstr)
  
  ! Define the vectors of a unit cube.
  vertices(1,:) = (/ 0, 0, 0 /)
  vertices(2,:) = (/ 0, 0, 1 /)
  vertices(3,:) = (/ 0, 1, 0 /)
  vertices(4,:) = (/ 0, 1, 1 /)
  vertices(5,:) = (/ 1, 0, 0 /)
  vertices(6,:) = (/ 1, 0, 1 /)
  vertices(7,:) = (/ 1, 1, 0 /)
  vertices(8,:) = (/ 1, 1, 1 /)

  !Define the centre of the system.
  centre = (/ tnx, tny /)/2

  ! Define the minimum xyz corner of my subdomain.
  base = ccoords(1:2) * (/ nx, ny /)

  do x = 1,nx
    do y = 1,ny
      do z = 1,nz
        tr = 0.0_rk
        tb = 0.0_rk
        do i = 1,8
          tx = x + base(1) + vertices(i,1)
          ty = y + base(2) + vertices(i,2)
          r2 = ( tx - centre(1) )**2 + ( ty - centre(2) )**2
          if ( r2 < r_c2 ) then
            tr = tr + fr
#ifndef SINGLEFLUID
            tb = tb + fb
#endif
          else
            tr = tr + pr
#ifndef SINGLEFLUID
            tb = tb + pb
#endif
          end if
        end do
        N(x,y,z)%n_r = tr*F(:)/8.0_rk
#ifndef SINGLEFLUID
        N(x,y,z)%n_b = tb*F(:)/8.0_rk
#endif
      end do
    end do
  end do
end subroutine lbe_init_cylinder_z

!> Sets up layers of \c fr1 sites of oil, \c fr2 sites of water, 
!> with the lamellae perpendicular to the y direction.
subroutine lbe_init_bi_lam_x(N)
  implicit none
  type(lbe_site),dimension(0:,0:,0:) :: N
  integer :: x, y, z, i
  real*8 :: theta, phi, r, F(19)
  integer :: lamwidth	! Width of one layer.
  integer :: tx	! y-coordinate in global system.


  lamwidth = fr1 + fr2

  if (myrankc == 0) then
    if (mod(tny,lamwidth) .ne. 0) then
      print*,'******* WARNING ******'
      print*,'* Truncated lamellae *'
      print*,'**********************'
    endif
  endif

  if (collisiontype_id .eq. MRT_id) then
    CALL mrt_init_dist((/0.0_8,0._8,0.0_8/),1.0d0,F(:))
  else
    CALL boltz_dist(0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,F(:))
  end if

  ! Note that ty is zero-based.
  do x=1,nx
    tx = x - 1 + ccoords(1)*nx
    i = mod(tx,lamwidth)

    if ( i<fr1 ) then
      ! Make oil happen.
      do y = 1,ny
        do z = 1,nz
          N(x,y,z)%n_r(:) = F(:)*fr
#ifndef SINGLEFLUID
          N(x,y,z)%n_b(:) = F(:)*fb
#endif
#ifndef NOSURFACTANT
          N(x,y,z)%n_s(:) = F(:)*fg

          if ( fd .eq. 0 ) then
            CALL random_number(theta)
            CALL random_number(phi)
            CALL random_number(r)
            theta = theta*pi ! 0 <= theta <= pi
            phi = phi*2.d0*pi ! 0 <=  phi  <= 2pi
            N(x,y,z)%d(1) = r*sin(theta)*cos(phi)
            N(x,y,z)%d(2) = r*sin(theta)*sin(phi)
            N(x,y,z)%d(3) = r*cos(theta)
          else
            N(x,y,z)%d(:) = (/ 0., 1., 0. /)
          end if
#endif
        end do
      end do
    else 
      ! Make water happen.
      do y = 1,ny
        do z = 1,nz
          N(x,y,z)%n_r(:) = F(:)*pr
#ifndef SINGLEFLUID
          N(x,y,z)%n_b(:) = F(:)*pb
#endif
#ifndef NOSURFACTANT
          N(x,y,z)%n_s(:) = F(:)*pg

          if ( qd .eq. 0 ) then
            CALL random_number(theta)
            CALL random_number(phi)
            CALL random_number(r)
            theta = theta*pi ! 0 <= theta <= pi
            phi = phi*2.d0*pi ! 0 <=  phi  <= 2pi
            N(x,y,z)%d(1) = r*sin(theta)*cos(phi)
            N(x,y,z)%d(2) = r*sin(theta)*sin(phi)
            N(x,y,z)%d(3) = r*cos(theta)
          else
            N(x,y,z)%d(:) = (/ 0., 1., 0. /)
          end if
#endif
        end do
      end do
    end if
  end do
end subroutine lbe_init_bi_lam_x

!> Sets up layers of \c fr1 sites of oil, \c fr2 sites of water, 
!> with the lamellae perpendicular to the y direction.
subroutine lbe_init_bi_lam_y(N)
  implicit none
  type(lbe_site),dimension(0:,0:,0:) :: N
  integer :: x, y, z, i
  real*8 :: theta, phi, r, F(19)
  integer :: lamwidth	! Width of one layer.
  integer :: ty	! y-coordinate in global system.

  lamwidth = fr1 + fr2

  if ( myrankc == 0 ) then
    if (mod(tny,lamwidth) .ne. 0) then
      print*,'******* WARNING ******'
      print*,'* Truncated lamellae *'
      print*,'**********************'
    endif
  endif

  if (collisiontype_id .eq. MRT_id) then
    CALL mrt_init_dist((/0.0_8,0._8,0.0_8/),1.0d0,F(:))
  else
    CALL boltz_dist(0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,F(:))
  end if

  ! Note that ty is zero-based.
  do y = 1,ny
    ty = y-1 + ccoords(2)*ny
    i = mod(ty,lamwidth)

    if (i<fr1) then
      ! Make oil happen.
      do x = 1,nx
        do z = 1,nz
          N(x,y,z)%n_r(:) = F(:)*fr
#ifndef SINGLEFLUID
          N(x,y,z)%n_b(:) = F(:)*fb
#endif
#ifndef NOSURFACTANT
          N(x,y,z)%n_s(:) = F(:)*fg

          if ( fd .eq. 0 ) then
            CALL random_number(theta)
            CALL random_number(phi)
            CALL random_number(r)
            theta = theta*pi ! 0 <= theta <= pi
            phi = phi*2.d0*pi ! 0 <=  phi  <= 2pi
            N(x,y,z)%d(1) = r*sin(theta)*cos(phi)
            N(x,y,z)%d(2) = r*sin(theta)*sin(phi)
            N(x,y,z)%d(3) = r*cos(theta)
          else
            N(x,y,z)%d(:) = (/ 0., 1., 0. /)
          end if
#endif
        end do
      end do
    else 
      ! Make water happen.
      do x = 1,nx
        do z = 1,nz
          N(x,y,z)%n_r(:) = F(:)*pr
#ifndef SINGLEFLUID
          N(x,y,z)%n_b(:) = F(:)*pb
#endif
#ifndef NOSURFACTANT
          N(x,y,z)%n_s(:) = F(:)*pg

          if ( pd .eq. 0 ) then
            CALL random_number(theta)
            CALL random_number(phi)
            CALL random_number(r)
            theta = theta*pi ! 0 <= theta <= pi
            phi = phi*2.d0*pi ! 0 <=  phi  <= 2pi
            N(x,y,z)%d(1) = r*sin(theta)*cos(phi)
            N(x,y,z)%d(2) = r*sin(theta)*sin(phi)
            N(x,y,z)%d(3) = r*cos(theta)
          else
            N(x,y,z)%d(:) = (/ 0., 1., 0. /)
          end if
#endif
        end do
      end do
    end if
  end do
end subroutine lbe_init_bi_lam_y

!> Sets up layers of \c fr1 sites of oil, \c fr2 sites of water,
!> with the lamellae perpendicular to the y direction.
subroutine lbe_init_bi_lam_y_fromv5(N)
  implicit none
  type(lbe_site), intent(inout), dimension(0:,0:,0:) :: N
  integer :: x,y,z,i
  real(kind=rk), dimension(nvecs) :: g_inv
  integer :: lamwidth	! Width of one layer.
  integer :: ty	! y-coordinate in global system.

  g_inv = 1.0_rk/g

  lamwidth = fr1 + fr2
  write(msgstr,"('lamwidth = ',I0)") lamwidth
  call log_msg(msgstr)

  if (myrankc == 0) then
    if (mod(tny,lamwidth) .ne. 0) then
      call log_msg("WARNING: truncated lamellae.")
    end if
  end if

  ! Note that ty is zero-based.
  do y=1,ny
    ty = y-1 + ccoords(2)*ny
    i = mod(ty,lamwidth)

    if (i<fr1) then
      ! Make oil happen.
      do x=1,nx
        do z=1,nz
          N(x,y,z)%n_r(:) = fr*g_inv/nvecs
#ifndef SINGLEFLUID
          N(x,y,z)%n_b(:) = 0.0_rk
#endif
#ifndef NOSURFACTANT
          N(x,y,z)%n_s(:) = 0.0_rk
          N(x,y,z)%d(:) = 0.0_rk
#endif
        end do
      end do
    else 
      ! Make water happen.
      do x=1,nx
        do z=1,nz
          N(x,y,z)%n_r(:) = 0.0_rk
#ifndef SINGLEFLUID
          N(x,y,z)%n_b(:) = fb*g_inv/nvecs
#endif

#ifndef NOSURFACTANT
          N(x,y,z)%n_s(:) = 0.0_rk
          N(x,y,z)%d(:) = 0.0_rk
#endif
        end do
      end do
    end if
  end do
end subroutine lbe_init_bi_lam_y_fromv5

!> Sets up a series of non-touching domains for HK benchmarking.
subroutine lbe_init_hk_test(N)
  implicit none
  type(lbe_site),dimension(0:,0:,0:) :: N
  integer :: x, y, z
  real(kind=rk) :: theta, phi, r, b, g, F(19)
  integer :: tz, ty, tx, mz, my, mx
  
  integer :: sieve_width, sieve_mod
 
  sieve_width = int(fr1)
  sieve_mod = sieve_width / 3
  
  if (collisiontype_id .eq. MRT_id) then
    call mrt_init_dist((/0.0_8,0._8,0.0_8/),1.0d0,F(:))
  else
    call boltz_dist(0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,F(:))
  end if

  ! Note that ty is zero-based.
  do z = 1,nz
    do y = 1,ny
      do x = 1,nx
        tz = z - 1 + ccoords(3)*nz
        ty = y - 1 + ccoords(2)*ny
        tx = x - 1 + ccoords(1)*nx

        mz = mod(tz, sieve_width)
        my = mod(ty, sieve_width)
        mx = mod(tx, sieve_width)
        if ( mz .lt. sieve_mod ) then
          if ( my .lt. sieve_mod ) then
            if ( mx .lt. sieve_mod ) then
              r = fr
              b = fb
              g = fg
            else
              r = pr
              b = pb
              g = pg
            end if
          else
            if ( mx .lt. sieve_mod ) then
              r = pr
              b = pb
              g = pg
            else
              r = fr
              b = fb
              g = fg
            end if
          end if
        else
          if ( my .lt. sieve_mod ) then
            if ( mx .lt. sieve_mod ) then
              r = pr
              b = pb
              g = pg
            else
              r = fr
              b = fb
              g = fg
            end if
          else
            if ( mx .lt. sieve_mod ) then
              r = fr
              b = fb
              g = fg
            else
              r = pr
              b = pb
              g = pg
            end if
          end if
        end if

        N(x,y,z)%n_r(:) = F(:)*r
#ifndef SINGLEFLUID
        N(x,y,z)%n_b(:) = F(:)*b
#endif
#ifndef NOSURFACTANT
        N(x,y,z)%n_s(:) = F(:)*g

        if ( fd .eq. 0 ) then
           call random_number(theta)
           call random_number(phi)
           call random_number(r)
           theta = theta*pi ! 0 <= theta <= pi
           phi = phi*2.d0*pi ! 0 <=  phi  <= 2pi
           N(x,y,z)%d(1) = r*sin(theta)*cos(phi)
           N(x,y,z)%d(2) = r*sin(theta)*sin(phi)
           N(x,y,z)%d(3) = r*cos(theta)
        else
           N(x,y,z)%d(:) = (/ 0., 1., 0. /)
        end if
#endif
        end do
      end do
   end do
end subroutine lbe_init_hk_test

!> Sets up layers of \c fr1 sites of oil, \c fr2 sites of water,
!> with the lamellae perpendicular to the z direction.
subroutine lbe_init_bi_lam_z(N)
  implicit none
  type(lbe_site),dimension(0:,0:,0:) :: N
  integer :: x, y, z, i
  real*8 :: theta, phi, r, F(19)
  integer :: lamwidth	! Width of one layer.
  integer :: tz	! z-coordinate in global system.

  lamwidth = fr1 + fr2

  if (myrankc == 0) then
    if (mod(tnz,lamwidth) .ne. 0) then
      print*,'******* WARNING ******'
      print*,'* Truncated lamellae *'
      print*,'**********************'
    endif
  endif

  if (collisiontype_id .eq. MRT_id) then
    CALL mrt_init_dist((/0.0_8,0._8,0.0_8/),1.0d0,F(:))
  else
    CALL boltz_dist(0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,F(:))
  end if

  ! Note that ty is zero-based.
  do z = 1,nz
    tz = z - 1 + ccoords(3)*nz
    i = mod(tz, lamwidth)

    if ( i < fr1 ) then
      ! Make oil happen.
      do x = 1,nx
        do y = 1,ny
          N(x,y,z)%n_r(:) = F(:)*fr
#ifndef SINGLEFLUID
          N(x,y,z)%n_b(:) = F(:)*fb
#endif
#ifndef NOSURFACTANT
          N(x,y,z)%n_s(:) = F(:)*fg

          if ( fd .eq. 0 ) then
            CALL random_number(theta)
            CALL random_number(phi)
            CALL random_number(r)
            theta = theta*pi ! 0 <= theta <= pi
            phi = phi*2.d0*pi ! 0 <=  phi  <= 2pi
            N(x,y,z)%d(1) = r*sin(theta)*cos(phi)
            N(x,y,z)%d(2) = r*sin(theta)*sin(phi)
            N(x,y,z)%d(3) = r*cos(theta)
          else
            N(x,y,z)%d(:) = (/ 0., 1., 0. /)
          end if
#endif
        end do
      end do
    else 
      ! Make water happen.
      do x = 1,nx
        do y = 1,ny
          N(x,y,z)%n_r(:) = F(:)*pr
#ifndef SINGLEFLUID
          N(x,y,z)%n_b(:) = F(:)*pb
#endif
#ifndef NOSURFACTANT
          N(x,y,z)%n_s(:) = F(:)*pg

          if ( pd .eq. 0 ) then
            CALL random_number(theta)
            CALL random_number(phi)
            CALL random_number(r)
            theta = theta*pi ! 0 <= theta <= pi
            phi = phi*2.d0*pi ! 0 <=  phi  <= 2pi
            N(x,y,z)%d(1) = r*sin(theta)*cos(phi)
            N(x,y,z)%d(2) = r*sin(theta)*sin(phi)
            N(x,y,z)%d(3) = r*cos(theta)
          else
            N(x,y,z)%d(:) = (/ 0., 1., 0. /)
          end if
#endif
        end do
      end do
    end if
  end do
end subroutine lbe_init_bi_lam_z

!> Sets up layers of \c fr1 sites of oil, \c fr2 sites of water, 
!> with the lamellae perpendicular to the z direction.
subroutine lbe_init_bi_sin_lam_x(N)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N
  integer :: x, y, z, i
  real*8  :: theta, phi, r, F(19) 
  integer :: lamwidth	! Width of one layer.
  integer :: tx, tz	! coordinates in global system.
  real*8  :: sinAmplitude, sinPeriods ! parameters of sinusoidal disturbance

  sinAmplitude = 5
  sinPeriods = 1

  lamwidth = fr1 + fr2

  if (myrankc == 0) then
    if (mod(tnz,lamwidth) .ne. 0) then
      print*,'******* WARNING ******'
      print*,'* Truncated lamellae *'
      print*,'**********************'
    endif
  endif

  if (collisiontype_id .eq. MRT_id) then
    CALL mrt_init_dist((/0.0_8,0._8,0.0_8/),1.0d0,F(:))
  else
    CALL boltz_dist(0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,F(:))
  end if

  ! Note that ty is zero-based.
  do z = 1,nz 
    tz = z - 1 + ccoords(3)*nz
    i = mod(tz, lamwidth)
    do x = 1,nx
      tx = x - 1 + ccoords(1)*nx
      do y = 1,ny
        if ( i > ( fr1 + sinAmplitude*sin( 1.0d0/sinPeriods*(2.d0*pi/nx)*x ) ) ) then
          ! Make oil happen.
          N(x,y,z)%n_r(:) = F(:)*fr
#ifndef SINGLEFLUID
          N(x,y,z)%n_b(:) = F(:)*fb
#endif
#ifndef NOSURFACTANT
          N(x,y,z)%n_s(:) = F(:)*fg

          if ( fd .eq. 0 ) then
            CALL random_number(theta)
            CALL random_number(phi)
            CALL random_number(r)
            theta = theta*pi ! 0 <= theta <= pi
            phi = phi*2.d0*pi ! 0 <=  phi  <= 2pi
            N(x,y,z)%d(1) = r*sin(theta)*cos(phi)
            N(x,y,z)%d(2) = r*sin(theta)*sin(phi)
            N(x,y,z)%d(3) = r*cos(theta)
          else
            N(x,y,z)%d(:) = (/ 0., 1., 0. /)
          end if
#endif
        else 
          N(x,y,z)%n_r(:) = F(:)*pr
#ifndef SINGLEFLUID
          N(x,y,z)%n_b(:) = F(:)*pb
#endif

#ifndef NOSURFACTANT
          N(x,y,z)%n_s(:) = F(:)*pg

          if ( pd .eq. 0 ) then
            CALL random_number(theta)
            CALL random_number(phi)
            CALL random_number(r)
            theta = theta*pi ! 0 <= theta <= pi
            phi = phi*2.d0*pi ! 0 <=  phi  <= 2pi
            N(x,y,z)%d(1) = r*sin(theta)*cos(phi)
            N(x,y,z)%d(2) = r*sin(theta)*sin(phi)
            N(x,y,z)%d(3) = r*cos(theta)
          else
            N(x,y,z)%d(:) = (/ 0., 1., 0. /)
          end if
#endif
        end if
      end do
    end do
  end do
end subroutine lbe_init_bi_sin_lam_x

!> This routine perturbs the system by up to a fraction \c perturbation.
!> It may be used to avoid metastable states.
!>
!> The occupation number on each site is multiplied by a random number
!> from a flat distrubution in the range
!> (1-perturbation)..(1+perturbation).
subroutine lbe_init_perturb(N)
  real*8, dimension(nvecs) :: wiggle
  integer :: x,y,z
  type(lbe_site),dimension(0:,0:,0:) :: N

  do z=1,nz
    do y=1,ny
      do x=1,nx
        call random_number(wiggle)
        wiggle = 1.d0 + perturbation*( 2.d0*wiggle - 1.d0 )
        N(x,y,z)%n_r(:) = N(x,y,z)%n_r(:)*wiggle
#ifndef SINGLEFLUID
        call random_number(wiggle)
        wiggle = 1.d0 + perturbation*( 2.d0*wiggle - 1.d0 )
        N(x,y,z)%n_b(:) = N(x,y,z)%n_b(:)*wiggle
#endif
        call random_number(wiggle)
#ifndef NOSURFACTANT
        wiggle = 1.d0 + perturbation*( 2.d0*wiggle - 1.d0 )
        N(x,y,z)%n_s(:) = N(x,y,z)%n_s(:)*wiggle
#endif
      end do
    end do
  end do

end subroutine lbe_init_perturb

!> This is a general routine for adding drops or vesicles to the system.
!> \c model contains an \c lbe_site type for each of the
!> three regions xE<lt>fr1, fr1E<lt>xE<gt>fr2, and xE<gt>fr2.
!> \c di contain the coordinates of the centre of teh drop/vesicle.
!> 
!> Each site is initialised to the appropriate value.
!> 
!>However, the dipole moment is ignored.
!>
!>The dipole moment of each site is set to have unit value, and
!>its orientation depends on the rock_state value for the region.
!>If \c rock_state is negative, the dipole points towards the centre of
!>the system (ie -r); if positive, it points away from the centre.
!>If zero, its orientation is chosen randomly.
!>
!>Currently, the values are averaged over the eight points of the unit
!>cube. This sucks.
subroutine lbe_add_radial(N,model,dx,dy,dz)
  implicit none
  type(lbe_site),dimension(0:,0:,0:) :: N
  integer :: x,y,z,s,dx,dy,dz
  integer, dimension(3) :: centre ! Global coords of centre of system.
  integer, dimension(3) :: base	! Global coords of corner of subdomain
  integer, dimension(3) :: offset
  integer, dimension(8,3) :: vertices ! Vertices of a unit cube
  real*8 :: if1, if2	! Real radii

  real*8, dimension(3) :: r		! Global coords of a given point
  real*8 :: rad			! Radius of point
  integer	:: i		! Loop variable over vertices

  real*8 :: n_r,n_b,n_s
  real*8, dimension(3) :: d,dd
  real*8, parameter :: A = 1./8.	
  real*8, dimension(nvecs) :: g_inv
  type(radial_model) :: model(2)  
  integer :: rad_index
  real*8 :: theta,phi

  g_inv = 1.d0/g

  ! All coordinates in this routine are in the global
  ! coordinate system.

  ! Define the actual radii of the spheres.
  if1 = fr1*min(tnx,tny,tnz)/2.d0
  if2 = fr2*min(tnx,tny,tnz)/2.d0
  if1 = fr1
  if2 = fr2

  !Define the centre of the system.
  centre = (/ dx, dy, dz /)
  ! Define the minimum xyz corner of my subdomain.
  base = ccoords* (/ nx, ny, nz /)
  ! Subtract the two.
  offset = base - centre

  ! Define the vectors of a unit cube.

  vertices(1,:) = (/ 0, 0, 0 /)
  vertices(2,:) = (/ 0, 0, 1 /)
  vertices(3,:) = (/ 0, 1, 0 /)
  vertices(4,:) = (/ 0, 1, 1 /)
  vertices(5,:) = (/ 1, 0, 0 /)
  vertices(6,:) = (/ 1, 0, 1 /)
  vertices(7,:) = (/ 1, 1, 0 /)
  vertices(8,:) = (/ 1, 1, 1 /)

  write(msgstr,"('if1 = ',F16.10,' if2 = ',F16.10)") if1, if2
  call log_msg(msgstr)

  ! 	if (myrankc == 0) then
  ! 		print*,'if12=',if1,if2
  ! 	endif

  ! Note that r is zero-based.

  do x=1,nx
    do y=1,ny
      do z=1,nz
        n_r=0.d0
        n_b=0.d0
        n_s=0.d0
        d=0.d0
        do i=1,8 ! Average over each corner of the site cube.
          r = (/ x-1,y-1,z-1 /)+offset + vertices(i,:) 
          rad = sqrt(sum(r*r))

          ! Find the appropriate index.

          if (rad < if1 ) then
            rad_index = rad_inner
          elseif ( (rad >= if1) .and. (rad <= if2) ) then
            rad_index = rad_middle
          else
            rad_index = rad_outer
          endif

          if (rad_index.ne.rad_outer) then

            ! Set appropriate values.

            n_r = n_r + model(rad_index)%n_r
            n_s = n_s + model(rad_index)%n_s
            n_b = n_b + model(rad_index)%n_b

            ! Handle the dipole moment.
            ! (but only if there are actually any dipoles..)

            if ( model(rad_index)%n_s > 0. ) then
              if ( model(rad_index)%dip < 0 ) then

                ! Dipole points towards centre.

                d = d - r/rad

              elseif ( model(rad_index)%dip == 0) then
                ! Make a random dipole of unit size.

                call random_number(theta)
                call random_number(phi)
                theta = theta*pi	! 0 <= theta <= pi
                phi = phi * 2.d0 * pi	! 0 <=  phi  <= 2pi
                dd(1) = sin(theta)*cos(phi)
                dd(2) = sin(theta)*sin(phi)
                dd(3) = cos(theta)

                ! Add it in to the averaging process.

                d = d + dd
              else
                ! Dipole points away from centre

                d = d + r/rad
              endif
            endif
          endif
        end do

        ! n_r, n_s, n_b are the required average particle densities
        ! for one site.

        if (rad_index.ne.rad_outer) then
          N(x,y,z)%n_r = A*n_r*g_inv/nvecs
          !N(x,y,z)%n_b = A*n_b*g_inv/nvecs
          !N(x,y,z)%n_s = A*n_s*g_inv/nvecs
          !N(x,y,z)%d = A*d
          !N(x,y,z)%rock_state=rad_index ! Useful for debugging.
        endif
      end do
    end do
  end do
end subroutine lbe_add_radial

!> register timers used in main LB code without MD
subroutine lbe_init_timers()

  call log_msg_ws("Registering timers ...")

  call register_timer('Total',ti_total)
  call register_timer('Advection',ti_adv)
  call register_timer('Int:Force',ti_intf)
  call register_timer('Int:Coll',ti_intc)
  call register_timer('HaloExchange',ti_halo)
  call register_timer('Invasion',ti_inv)
  call register_timer('Dump',ti_dump)

#ifdef IBM_PART
  call log_msg_ws("Registering IBM timers ...")
  call register_timer('IBM_init', ti_IBM_init)
  call register_timer('IBM_intspread', ti_IBM_intspread)
  call register_timer('IBM_forces', ti_IBM_forces)
  call register_timer('IBM_MPI', ti_IBM_MPI)
  call register_timer('IBM_update', ti_IBM_update)
  call register_timer('IBM_dump', ti_IBM_dump)
#endif

end subroutine lbe_init_timers

!> initialize default names for fluxz regions
subroutine lbe_setup_fluxz()
  integer :: i

#ifdef MD
  provide_uid2i = .true.
#endif

  do i=1,fluxz_regions
    write (fluxz_name(i),fmt='(I0)') i
  end do
end subroutine lbe_setup_fluxz

!> Initialize unique seeds
subroutine lbe_unique_seeds()
  integer :: isize,stat
  integer,allocatable :: seedarray(:)

  ! Now make sure that the initial random seed given to each processor
  ! is unique. BEWARE this is only a temporary measure - I know not
  ! if this will produce any hidden correlations between the random
  ! numbers. I leave this to the experts - it should not be difficult
  ! to modify the routine below.

  if ( seed .eq. 0 ) then
    call log_msg("WARNING: seed is set to zero - this sets equal seeds (0) on all ranks and might cause undefined behaviour.")
  end if

  ! ensure the random number seed remains in range
  seed = mod(seed * (myrankc + 1) * 10101, 300000)

  ! seed the builtin RNG
  call random_seed(size=isize)
  allocate (seedarray(isize), stat=stat)
  call check_allocate(stat, 'lbe_unique_seeds(): seedarray')
  seedarray = seed
  call random_seed(put=seedarray(1:isize))
  deallocate (seedarray)

  ! also seed the luxury_rng
  call rluxgo(4,seed,0,0)
end subroutine lbe_unique_seeds

end module lbe_init_module
