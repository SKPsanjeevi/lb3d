#include "lbe.h"

!> Invasive flow
!>
!>This file contains code which deals with invasive flow: routines to set
!>up appropriate buffering regions on startup, and routines to recolour
!>particles after advection during each timestep.
module lbe_invasion_module
  use lbe_globals_module
  use lbe_helper_module
  use lbe_log_module
  use lbe_parallel_module
  use lbe_parms_module, only: amass_b,amass_r,amass_s,fb,fg,fr,fr1,fr2,in_evp&
       &,inv_type,m_evp,m_evp_gr,m_evp_gb,m_evp_freq_f,m_evp_freq_a,nt,out_evp,m_evp_set_density,pb,pg,pr,pd,shear_u,SCMP,tau_b,tau_r, inv_fluid,nx,ny,nz,nt
  use lbe_bdist_module, only: boltz_dist
  use lbe_collision_module, only: lbe_calculate_sc_forces,zero_force_offset
  use lbe_types_module, only: lbe_site

  implicit none
  private

  public lbe_invade,lbe_invasion_adjust_halo,lbe_invade_after_evaporation 

contains

!> This is a wrapper routine which calls the invasion routine appropriate
!> to the value of \c inv_fluid.
subroutine lbe_invade(N)
  implicit none
  type(lbe_site), dimension(0:,0:,0:), intent(inout) :: N
  select case (inv_fluid)
  case (0)
    call lbe_invade_sampled(N)
  case (1)
    call lbe_invade_mixture(N,1.0_rk,0.0_rk,0.0_rk)
  case (2)
#ifndef SINGLEFLUID
    call lbe_invade_mixture(N,0.0_rk,0.0_rk,1.0_rk)
#else
    call error("Water invasion not supported with -DSINGLEFLUID.")
#endif
  case (3)
#ifndef NOSURFACTANT
    call lbe_invade_mixture(N,0.0_rk,1.0_rk,0.0_rk)
#else
    call error("Surfactant invasion not supported without surf.")
#endif
  case (4)
#ifndef NOSURFACTANT
    call lbe_invade_mixture(N,pr,pg,pb)
#else
    call lbe_invade_mixture(N,pr,0.0_rk,pb)
#endif
  case (5)
    !This is shear - do nothing here.
  case (6)
    !This is rotated bc's - do nothing here.
#ifdef NOSURFACTANT
  case (7)
    call lbe_invade_binary(N,1.0_rk,0.0_rk)
  case (8)
    call lbe_invade_binary(N,0.0_rk,1.0_rk)
#endif
  case (9)
    call lbe_invade_constmixture(N,fr,fg,fb,pr,pg,pb)
  case (10)
    call lbe_invade_constmixooutfr(N,pr,pg,pb,fr1,0.0_rk,fr2)
  case (11) 
    call lbe_invade_zouandhe(N,fr,fb,pr,pb,inv_type)
  case (12)
    call lbe_invade_zouandheflux(N,pr,pg,pb,inv_type) 
  case (13)
    call lbe_invade_zouandhe_angle(N,pr,pg,pb,fr1,fr2,inv_type)
! Martin: Hybrid-boundary condition: currently in a seperate branch version5-hybrid
! case (14)
!   call lbe_invade_hybrid(N)
! Martin: experimental boundaries contain some UNPHYSICAL cases and are not yet
! ready for use in production runs
  case (15)
    call lbe_invade_equilibrium(N,fr,fg,fb,pr,pg,pb,fr1,fr2,inv_type)
  case (16)
    call lbe_invade_poiseuille(N,inv_type)
  case (17)
    call lbe_invade_with_slip(N,pr,pb,pg,shear_u,inv_type)
  case (18)
    call lbe_invade_channel(N,fr,pr,pb,pg,inv_type)
  case (20)
    call lbe_invade_constmixoout(N,pr,pg,pb)
  case (21)
    call lbe_invade_constmixooutfr(N,pr,pg,pb,fr1,0.0_rk,fr2)
  case (22)
    call lbe_invade_binary_lamellae(N,1.0_rk,1.0_rk,fr1,fr2)
  case (23)
    call lbe_invade_zouandhemix(N,fr,fg,fb,pr,pg,pb)
  case (26)
    call lbe_invade_evaporation(N,m_evp,m_evp_gr,m_evp_gb,m_evp_freq_f,m_evp_freq_a,in_evp,out_evp,m_evp_set_density)
  case default
    write(msgstr,"('Illegal value: inv_fluid = ',I0)") inv_fluid
    call error(msgstr)
  end select

end subroutine lbe_invade

!> modifies a newly exchanged halo if required by the current setting
!> of \c inv_fluid
!>
!> \param[in,out] whole_N local chunk of the lattice with full
!> halo of depth \c halo_extent (as introduced in version5-md)
subroutine lbe_invasion_adjust_halo(whole_N)
  implicit none
  type(lbe_site),intent(inout) :: &
       &whole_N(1-halo_extent:,1-halo_extent:,1-halo_extent:)

  select case (inv_fluid)
  case (10)
    call extrapolate_maxz(whole_N)
  case (11)
    call mirror_rock_minz(whole_N)
    call mirror_rock_maxz(whole_N)
  case (12)
    call mirror_rock_minz(whole_N)
    call mirror_rock_maxz(whole_N)
  case (13)
    call mirror_rock_minz(whole_N)
    call mirror_rock_maxz(whole_N)
  case (20)
    call extrapolate_maxz(whole_N)
  case (21)
    call extrapolate_maxz(whole_N)
  case (23)
    call mirror_rock_minz(whole_N)
    call mirror_rock_maxz(whole_N)
  end select
end subroutine lbe_invasion_adjust_halo

!> At boundary regions, recolour all particles entering the boundary region
!> to be in the ratio \c mr red, \c ms surfactant, \c mb blue.
!> 
!> This currently only handles the Z boundaries. Adding the X and Y
!> boundaries, if desired, is a matter of cut, paste, and a few substitutions.
subroutine lbe_invade_mixture(N,mr,ms,mb)
	implicit none
	type(lbe_site),dimension(0:,0:,0:) :: N
	real*8 :: mr,ms,mb

	integer :: x,y,z
	logical minx,maxx,miny,maxy,minz,maxz
	real*8	:: pop
	integer :: i,j

	! Determine whether I sit on any boundaries.

	minx=(start(1)==1)
	miny=(start(2)==1)
	minz=(start(3)==1)
	maxx=(start(1)>=(tnx-nx))
	maxy=(start(2)>=(tny-ny))
	maxz=(start(3)>=(tnz-nz))

	if (minz) then
		do y=1,ny
		 do x=1,nx
			if (N(x,y,1)%rock_state == 0) then
				do i=1,nnp
					j=posz(i)

					pop=	N(x,y,1)%n_r(j) 
#ifndef SINGLEFLUID
					pop = pop + N(x,y,1)%n_b(j)
#endif
#ifndef NOSURFACTANT
					pop=	pop +                   &
						N(x,y,1)%n_s(j) 
					N(x,y,1)%n_s(j)=ms*pop
#endif
					N(x,y,1)%n_r(j)=mr*pop
#ifndef SINGLEFLUID
					N(x,y,1)%n_b(j)=mb*pop
#endif
				end do
			endif
			if (N(x,y,2)%rock_state == 0) then
				do i=1,nnp
		 		j=negz(i)
					pop=	N(x,y,2)%n_r(j)
#ifndef SINGLEFLUID
					pop = pop + N(x,y,2)%n_b(j)
#endif
#ifndef NOSURFACTANT
					pop=	pop             +       &
						N(x,y,2)%n_s(j)
					N(x,y,2)%n_s(j)=ms*pop
#endif
					N(x,y,2)%n_r(j)=mr*pop
#ifndef SINGLEFLUID
					N(x,y,2)%n_b(j)=mb*pop
#endif
				end do
			endif
		 end do
		end do
	endif ! minz
!	if (maxz) then
!		do y=1,ny
!		 do x=1,nx
!			if (N(x,y,tnz-1)%rock_state == 0) then
!				do i=1,nnp
!					j=posz(i)
!					pop=	N(x,y,tnz-1)%n_r(j) +	&
!						N(x,y,tnz-1)%n_s(j) +	&
!						N(x,y,tnz-1)%n_b(j)
!					N(x,y,tnz-1)%n_r(j)=mr*pop
!					N(x,y,tnz-1)%n_b(j)=mb*pop
!					N(x,y,tnz-1)%n_s(j)=ms*pop
!				end do
!			endif
!			if (N(x,y,tnz)%rock_state == 0) then
!				do i=1,nnp
!					j=negz(i)
!					pop=	N(x,y,tnz)%n_r(j) +	&
!						N(x,y,tnz)%n_s(j) +	&
!						N(x,y,tnz)%n_b(j)
!					N(x,y,tnz)%n_r(j)=mr*pop
!					N(x,y,tnz)%n_b(j)=mb*pop
!					N(x,y,tnz)%n_s(j)=ms*pop
!				end do
!			endif
!		 end do
!		end do
!	endif ! maxz
!	if (miny) then
!		do z=1,nz
!		 do x=1,nx
!			if (N(x,1,z)%rock_state == 0) then
!				do i=1,nnp
!					j=posy(i)
!					pop=	N(x,1,z)%n_r(j) +	&
!						N(x,1,z)%n_s(j) +	&
!						N(x,1,z)%n_b(j)
!					N(x,1,z)%n_r(j)=mr*pop
!					N(x,1,z)%n_b(j)=mb*pop
!					N(x,1,z)%n_s(j)=ms*pop
!				end do
!			endif
!			if (N(x,2,z)%rock_state == 0) then
!				do i=1,nnp
!					j=negy(i)
!					pop=	N(x,2,z)%n_r(j) +	&
!						N(x,2,z)%n_s(j) +	&
!						N(x,2,z)%n_b(j)
!					N(x,2,z)%n_r(j)=mr*pop
!					N(x,2,z)%n_b(j)=mb*pop
!					N(x,2,z)%n_s(j)=ms*pop
!				end do
!			endif
!		 end do
!		end do
!	endif ! miny
!	if (minx) then
!		do y=1,ny
!		 do z=1,nz
!			if (N(1,y,z)%rock_state == 0) then
!				do i=1,nnp
!					j=posx(i)
!					pop=	N(1,y,z)%n_r(j) +	&
!						N(1,y,z)%n_s(j) +	&
!						N(1,y,z)%n_b(j)
!					N(1,y,z)%n_r(j)=mr*pop
!					N(1,y,z)%n_b(j)=mb*pop
!					N(1,y,z)%n_s(j)=ms*pop
!				end do
!			endif
!			if (N(2,y,z)%rock_state == 0) then
!				do i=1,nnp
!					j=negx(i)
!					pop=	N(2,y,z)%n_r(j) +	&
!						N(2,y,z)%n_s(j) +	&
!						N(2,y,z)%n_b(j)
!					N(2,y,z)%n_r(j)=mr*pop
!					N(2,y,z)%n_b(j)=mb*pop
!					N(2,y,z)%n_s(j)=ms*pop
!				end do
!			endif
!		 end do
!		end do
!	endif ! minx

end subroutine lbe_invade_mixture

!> Constantly set densities at z=0,1 to to be \c inr red, \c ing surfactant,
!> \c inb blue and at z=tnz,tnz+1 to be \c outr red, \c outg surfactant, 
!> \c outb blue. This is useful for pressure driven flow.
subroutine lbe_invade_constmixture(N,inr,ing,inb,outr,outg,outb)
	implicit none
	type(lbe_site),dimension(0:,0:,0:) :: N
	real*8 :: inr,ing,inb,outr,outg,outb

	integer :: x,y,z
	logical minx,maxx,miny,maxy,minz,maxz
	real*8	:: pop
	integer :: i,j,k

	! Determine whether I sit on any boundaries.

	minx=(start(1)==1)
	miny=(start(2)==1)
	minz=(start(3)==1)
	maxx=(start(1)>=(tnx-nx))
	maxy=(start(2)>=(tny-ny))
	maxz=(start(3)>=(tnz-nz))

	if (minz) then
		do y=1,ny
		 do x=1,nx
			if (N(x,y,0)%rock_state == 0) then
!				do i=1,nnp
!					j=posz(i)

#ifndef NOSURFACTANT
					N(x,y,0)%n_s(:)=ing/nvecs/g(:)
#endif
					N(x,y,0)%n_r(:)=inr/nvecs/g(:)
#ifndef SINGLEFLUID
					N(x,y,0)%n_b(:)=inb/nvecs/g(:)
#endif
!				end do
			endif
			if (N(x,y,1)%rock_state == 0) then
!				do i=1,nnp
!		 		j=negz(i)
#ifndef NOSURFACTANT
					N(x,y,1)%n_s(:)=ing/nvecs/g(:)
#endif
					N(x,y,1)%n_r(:)=inr/nvecs/g(:)
#ifndef SINGLEFLUID
					N(x,y,1)%n_b(:)=inb/nvecs/g(:)
#endif
!				end do
			endif
		 end do
		end do
	endif ! minz

	if (maxz) then
		do y=1,ny
		 do x=1,nx
                        if (N(x,y,nz+1)%rock_state == 0) then
!                               do i=1,nnp
!                                       j=negz(i)
#ifndef NOSURFACTANT
                                        N(x,y,nz+1)%n_s(:)=outg/nvecs/g(:)
#endif
                                        N(x,y,nz+1)%n_r(:)=outr/nvecs/g(:)
#ifndef SINGLEFLUID
                                        N(x,y,nz+1)%n_b(:)=outb/nvecs/g(:)
#endif
!                               end do
                        endif
                        if (N(x,y,nz)%rock_state == 0) then
!                               do i=1,nnp
!                               j=posz(i)
#ifndef NOSURFACTANT
                                        N(x,y,nz)%n_s(:)=outg/nvecs/g(:)
#endif
                                        N(x,y,nz)%n_r(:)=outr/nvecs/g(:)
#ifndef SINGLEFLUID
                                        N(x,y,nz)%n_b(:)=outb/nvecs/g(:)
#endif
!                               end do
                        endif
                 end do
                end do
        endif ! maxz


end subroutine lbe_invade_constmixture

!> Constantly set densities at z=0,1 to to be \c inr red, \c ing surfactant,
!> \c inb blue. At z=tnz,tnz+1 an
!> infinitely long system is simulated by interpolating values at z=tnz+1
!> from values at z=tnz,tnz-1. 
!> This is useful for pressure driven flow.
subroutine lbe_invade_constmixoout(N,inr,ing,inb)
	implicit none
	type(lbe_site),dimension(0:,0:,0:) :: N
	real*8 :: inr,ing,inb
	real*8 :: theta,phi,r

	integer :: x,y,z
	logical minx,maxx,miny,maxy,minz,maxz
	real*8	:: pop
	integer :: i,j,k

	! Determine whether I sit on any boundaries.

	minx=(start(1)==1)
	miny=(start(2)==1)
	minz=(start(3)==1)
	maxx=(start(1)>=(tnx-nx))
	maxy=(start(2)>=(tny-ny))
	maxz=(start(3)>=(tnz-nz))

	if (minz) then
		do y=1,ny
		 do x=1,nx
			if (N(x,y,0)%rock_state == 0) then
				do i=1,nnp
					j=posz(i)

#ifndef NOSURFACTANT
                                call random_number(theta)
                                call random_number(phi)
                                call random_number(r)   ! 0 <=   r   <= 1
                                theta = theta*pi        ! 0 <= theta <= pi
                                phi = phi * 2. * pi     ! 0 <=  phi  <= 2pi
                                N(x,y,0)%d(1) = r*sin(theta)*cos(phi)
                                N(x,y,0)%d(2) = r*sin(theta)*sin(phi)
                                N(x,y,0)%d(3) = r*cos(theta)

                                        N(x,y,0)%n_s(j)=ing/nvecs/g(j)
#endif
                                        N(x,y,0)%n_r(j)=inr/nvecs/g(j)
#ifndef SINGLEFLUID
                                        N(x,y,0)%n_b(j)=inb/nvecs/g(j)
#endif
                                        k=negz(i)

#ifndef NOSURFACTANT
                                        N(x,y,0)%n_s(k)=0.d0
#endif
                                        N(x,y,0)%n_r(k)=0.d0
#ifndef SINGLEFLUID
                                        N(x,y,0)%n_b(k)=0.d0
#endif
                                end do
                        endif
                        if (N(x,y,1)%rock_state == 0) then
                                do i=1,nnp
                                j=negz(i)
#ifndef NOSURFACTANT
                                call random_number(theta)
                                call random_number(phi)
                                call random_number(r)   ! 0 <=   r   <= 1
                                theta = theta*pi        ! 0 <= theta <= pi
                                phi = phi * 2.d0 * pi     ! 0 <=  phi  <= 2pi
                                N(x,y,1)%d(1) = r*sin(theta)*cos(phi)
                                N(x,y,1)%d(2) = r*sin(theta)*sin(phi)
                                N(x,y,1)%d(3) = r*cos(theta)

                                N(x,y,1)%n_s(j)=ing/nvecs/g(j)
#endif
                                        N(x,y,1)%n_r(j)=inr/nvecs/g(j)
#ifndef SINGLEFLUID
                                        N(x,y,1)%n_b(j)=inb/nvecs/g(j)
#endif
                                end do
                        endif
                 end do
                end do
        endif ! minz

end subroutine lbe_invade_constmixoout


!> Constantly set densities at z=0,1 to to be \c inr red, \c ing surfactant,
!> \c inb blue. frr,frg,frg give x pos for different fluids. At z=tnz,tnz+1 an
!> infinitely long system is simulated by interpolating values at z=tnz+1
!> from values at z=tnz,tnz-1. 
!> This is useful for pressure driven flow.
subroutine lbe_invade_constmixooutfr(N,inr,ing,inb,frr,frg,frb)
	implicit none
	type(lbe_site),dimension(0:,0:,0:) :: N
	real*8 :: inr,ing,inb,frr,frg,frb

	integer :: x,y,z
	logical minx,maxx,miny,maxy,minz,maxz
	real*8	:: pop
	integer :: i,j,k,l
	integer :: lamwidth     ! Width of one layer.
	integer :: ty   ! y-coordinate in global system.

	! Determine whether I sit on any boundaries.

	minx=(start(1)==1)
	miny=(start(2)==1)
	minz=(start(3)==1)
	maxx=(start(1)>=(tnx-nx))
	maxy=(start(2)>=(tny-ny))
	maxz=(start(3)>=(tnz-nz))

	if (minz) then

	  lamwidth = frr + frg + frb
	  if ((myrankc == 0).AND.(nt==0)) then
                if (mod(tny,lamwidth) .ne. 0) then
                        print*,'******* WARNING ******'
                        print*,'* Truncated lamellae *'
                        print*,'**********************'
                endif
          endif



		do y=1,ny

                   ty = y-1 + ccoords(2)*ny
                   l=mod(ty,lamwidth)

		   if (l.lt.frr) then
		    ! Make oil happen.
		    do x=1,nx
			if (N(x,y,0)%rock_state == 0) then
				do i=1,nnp
					j=posz(i)

#ifndef NOSURFACTANT
					N(x,y,0)%n_s(j)=0.d0
					N(x,y,0)%d(:)=0.d0
#endif
					N(x,y,0)%n_r(j)=inr/nvecs/g(j)
#ifndef SINGLEFLUID
					N(x,y,0)%n_b(j)=0.d0
#endif

					k=negz(i)
#ifndef NOSURFACTANT
					N(x,y,0)%n_s(k)=0.d0
#endif
					N(x,y,0)%n_r(k)=0.d0
#ifndef SINGLEFLUID
					N(x,y,0)%n_b(k)=0.d0
#endif
				end do
			endif
			if (N(x,y,1)%rock_state == 0) then
				do i=1,nnp
		 		j=negz(i)
#ifndef NOSURFACTANT
					N(x,y,1)%n_s(j)=0.d0
					N(x,y,1)%d(:)=0.d0
#endif
					N(x,y,1)%n_r(j)=inr/nvecs/g(j)
#ifndef SINGLEFLUID
					N(x,y,1)%n_b(j)=0.d0
#endif
				end do
			endif
		 end do
		elseif ((l == frr) .or. (l.lt.(frr+frg))) then
		   ! Make surf happen.
		    do x=1,nx
			if (N(x,y,0)%rock_state == 0) then
				do i=1,nnp
					j=posz(i)

#ifndef NOSURFACTANT
					N(x,y,0)%n_s(j)=ing/nvecs/g(j)
					N(x,y,0)%d(:)= (/ 0.d0, 1.d0, 0.d0 /)
#endif
					N(x,y,0)%n_r(j)=0.d0
#ifndef SINGLEFLUID
					N(x,y,0)%n_b(j)=0.d0
#endif
					k=negz(i)

#ifndef NOSURFACTANT
					N(x,y,0)%n_s(k)=0.d0
#endif
					N(x,y,0)%n_r(k)=0.d0
#ifndef SINGLEFLUID
					N(x,y,0)%n_b(k)=0.d0
#endif
				end do
			endif
			if (N(x,y,1)%rock_state == 0) then
				do i=1,nnp
		 		j=negz(i)
#ifndef NOSURFACTANT
					N(x,y,1)%n_s(j)=ing/nvecs/g(j)
					N(x,y,1)%d(:)= (/ 0.d0, 1.d0, 0.d0  /)
#endif
					N(x,y,1)%n_r(j)=0.d0
#ifndef SINGLEFLUID
					N(x,y,1)%n_b(j)=0.d0
#endif
				end do
			endif
		 end do
		elseif ((l==frr+frg) .or. (l.lt.(frr+frg+frb))) then
		   ! Make water happen.
		    do x=1,nx
			if (N(x,y,0)%rock_state == 0) then
				do i=1,nnp
					j=posz(i)

#ifndef NOSURFACTANT
					N(x,y,0)%n_s(j)=0.d0
					N(x,y,0)%d(:)= 0.d0
#endif
					N(x,y,0)%n_r(j)=0.d0
#ifndef SINGLEFLUID
					N(x,y,0)%n_b(j)=inb/nvecs/g(j)
#endif
					k=negz(i)

#ifndef NOSURFACTANT
					N(x,y,0)%n_s(k)=0.d0
#endif
					N(x,y,0)%n_r(k)=0.d0
#ifndef SINGLEFLUID
					N(x,y,0)%n_b(k)=0.d0
#endif
				end do
			endif
			if (N(x,y,1)%rock_state == 0) then
				do i=1,nnp
		 		j=negz(i)
#ifndef NOSURFACTANT
					N(x,y,1)%n_s(j)=0.d0
					N(x,y,1)%d(:)= 0.d0
#endif
					N(x,y,1)%n_r(j)=0.d0
#ifndef SINGLEFLUID
					N(x,y,1)%n_b(j)=inb/nvecs/g(j)
#endif
				end do
			endif
		 end do
		endif

		end do
	endif ! minz

! I had to put the maxz stuff in lbe_parallel directly and know that that is
! bad... Jens, 14.03.05


end subroutine lbe_invade_constmixooutfr

!> If the CPU calling this routine is sitting at an edge of the global
!> system, it examines the buffer region within two sites of the edge.
!> All particles which have just entered this region are recoloured to match
!> the ratio of the adjacent inward lattice site.
!> 
!> Warning: this code has not been extensively tested. It does not conserve
!> momentum if non-unit masses are used.
subroutine lbe_invade_sampled(N)
  implicit none
  type(lbe_site),dimension(0:,0:,0:), intent(inout) :: N
  integer :: x, y, z, i, j
  logical minx, maxx, miny, maxy, minz, maxz
  real(kind=rk)	:: nr, pop
#ifndef SINGLEFLUID
  real(kind=rk) :: nb
#ifndef NOSURFACTANT
  real(kind=rk) :: ns
#endif
#endif

  ! Determine whether I sit on any boundaries.
  minx=(start(1)==1)
  miny=(start(2)==1)
  minz=(start(3)==1)
  maxx=(start(1)>=(tnx-nx))
  maxy=(start(2)>=(tny-ny))
  maxz=(start(3)>=(tnz-nz))

  if (minz) then
    do y=1,ny
      do x=1,nx
        if ( is_fluid(N(x,y,3)%rock_state) ) then
          nr = sum(N(x,y,3)%n_r(:)*g)
          pop = nr
#ifndef SINGLEFLUID
          nb = sum(N(x,y,3)%n_b(:)*g)
          pop = pop + nb
#ifndef NOSURFACTANT
          ns = sum(N(x,y,3)%n_s(:)*g)
          pop = pop + ns
#endif
#endif
        else
          pop = 0.0_rk
        end if

        ! nr,nb,ns end up containing the ratio of
        ! red to green to surf in the inward adjacent site.
        if (pop .eq. 0.0_rk) then
          ! Assign 1:1:1 ratio if adjacent site is empty.
          ! (or if it's a rock site).
#ifndef NOSURFACTANT
          nr = 1.0_rk/3.0_rk
          ns = 1.0_rk/3.0_rk
          nb = 1.0_rk/3.0_rk
#else
          nr = 1.0_rk/2.0_rk
#ifndef SINGLEFLUID
          nb = 1.0_rk/2.0_rk
#endif
#endif
        else
          nr = nr / pop
#ifndef SINGLEFLUID
          nb = nb / pop
#ifndef NOSURFACTANT
          ns = ns / pop
#endif
#endif
        end if
        do i=1, nnp
          ! Recolour the particles which are headed in the -Z direction.
          if ( is_fluid(N(x,y,2)%rock_state) ) then
            j = negz(i)
            pop = N(x,y,2)%n_r(j)
#ifndef SINGLEFLUID
            pop = pop + N(x,y,2)%n_b(j)
#ifndef NOSURFACTANT
            pop = pop + N(x,y,2)%n_s(j)
#endif
#endif
            ! Note that particle number should be conserved
            ! give or take a bit of roundoff error.
            ! (but we're using 64-bit reals, so that doesn't matter.
            ! Right? Right? FIXME!)
            N(x,y,2)%n_r(j) = pop * nr
#ifndef SINGLEFLUID
            N(x,y,2)%n_b(j) = pop * nb
#ifndef NOSURFACTANT
            N(x,y,2)%n_s(j) = pop * ns
#endif
#endif
          end if
          ! Now recolour the ones headed in the +Z direction.
          if ( is_fluid(N(x,y,1)%rock_state) ) then 
            j = posz(i)
            pop = N(x,y,1)%n_r(j)
#ifndef SINGLEFLUID
            pop = pop + N(x,y,1)%n_b(j)
#ifndef NOSURFACTANT
            pop = pop + N(x,y,1)%n_s(j)
#endif
#endif
            N(x,y,1)%n_r(j) = pop * nr
#ifndef SINGLEFLUID
            N(x,y,1)%n_b(j) = pop * nb
#ifndef NOSURFACTANT
            N(x,y,1)%n_s(j) = pop * ns
#endif
#endif
          end if
        end do
      end do
    end do
  end if ! minz

  if (maxz) then
    do y=1,ny
      do x=1,nx
        if ( is_fluid(N(x,y,nz-2)%rock_state) ) then
          nr = sum(N(x,y,nz-2)%n_r(:)*g)
          pop = nr
#ifndef SINGLEFLUID
          nb = sum(N(x,y,nz-2)%n_b(:)*g)
          pop = nr + nb
#ifndef NOSURFACTANT
          ns = sum(N(x,y,nz-2)%n_s(:)*g)
          pop = pop + ns
#endif
#endif
        else
          pop = 0.0_rk
        end if

        ! nr,nb,ns end up containing the ratio of
        ! red to green to surf in the inward adjacent site.
        if (pop .eq. 0.0_rk) then
          ! Assign 1:1:1 ratio if adjacent site is empty.
          ! (or if it's a rock site).
#ifndef NOSURFACTANT
          nr = 1.0_rk/3.0_rk
          ns = 1.0_rk/3.0_rk
          nb = 1.0_rk/3.0_rk
#else
          nr = 1.0_rk/2.0_rk
#ifndef SINGLEFLUID
          nb = 1.0_rk/2.0_rk
#endif
#endif
        else
          nr = nr / pop
#ifndef SINGLEFLUID
          nb = nb / pop
#ifndef NOSURFACTANT
          ns = ns / pop
#endif
#endif
        end if

        do i=1,nnp
          ! Recolour the particles which are headed in the +Z direction.
          if ( is_fluid(N(x,y,nz-1)%rock_state) ) then
            j = posz(i)
            pop = N(x,y,nz-1)%n_r(j)
#ifndef SINGLEFLUID
            pop = pop + N(x,y,nz-1)%n_b(j)
#ifndef NOSURFACTANT
            pop = pop + N(x,y,nz-1)%n_s(j)
#endif
#endif

            N(x,y,nz-1)%n_r(j) = pop * nr
#ifndef SINGLEFLUID
            N(x,y,nz-1)%n_b(j) = pop * nb
#ifndef NOSURFACTANT
            N(x,y,nz-1)%n_s(j) = pop * ns
#endif
#endif
          end if
          ! Now recolour the ones headed in the -Z direction.
          if ( is_fluid(N(x,y,nz)%rock_state) ) then 
            j = negz(i)
            pop = N(x,y,nz)%n_r(j)
#ifndef SINGLEFLUID
            pop = pop + N(x,y,nz)%n_b(j)
#ifndef NOSURFACTANT
            pop = pop + N(x,y,nz)%n_s(j)
#endif
#endif
            N(x,y,nz)%n_r(j)=pop*nr
#ifndef SINGLEFLUID
            N(x,y,nz)%n_b(j)=pop*nb
#ifndef NOSURFACTANT
            N(x,y,nz)%n_s(j)=pop*ns
#endif
#endif
          end if
        end do
      end do
    end do
  end if ! maxz

  if (maxy) then
    do z=1,nz
      do x=1,nx
        if ( is_fluid(N(x,ny-2,z)%rock_state) ) then
          nr = sum(N(x,ny-2,z)%n_r(:)*g)
          pop = nr
#ifndef SINGLEFLUID
          nb = sum(N(x,ny-2,z)%n_b(:)*g)
          pop = pop + nb
#ifndef NOSURFACTANT
          ns = sum(N(x,ny-2,z)%n_s(:)*g)
          pop = pop + ns
#endif
#endif
        else
          pop = 0.0_rk
        end if
        ! nr,nb,ns end up containing the ratio of
        ! red to green to surf in the inward adjacent site.

        if (pop .eq. 0.0_rk) then
          ! Assign 1:1:1 ratio if adjacent site is empty.
          ! (or if it's a rock site).
#ifndef NOSURFACTANT
          nr = 1.0_rk/3.0_rk
          ns = 1.0_rk/3.0_rk
          nb = 1.0_rk/3.0_rk
#else
          nr = 1.0_rk/2.0_rk
#ifndef SINGLEFLUID
          nb = 1.0_rk/2.0_rk
#endif
#endif
        else
          nr = nr / pop
#ifndef SINGLEFLUID
          nb = nb / pop
#ifndef NOSURFACTANT
          ns = ns / pop
#endif
#endif
        end if
        do i=1,nnp
          ! Recolour the particles which are headed in the +Y
          ! direction.
          if ( is_fluid(N(x,ny-1,z)%rock_state) ) then
            j = posy(i)
            pop = N(x,ny-1,z)%n_r(j)
#ifndef SINGLEFLUID
            pop = pop + N(x,ny-1,z)%n_b(j)
#ifndef NOSURFACTANT
            pop = pop + N(x,ny-1,z)%n_s(j)

#endif
#endif
            N(x,ny-1,z)%n_r(j) = pop * nr
#ifndef SINGLEFLUID
            N(x,ny-1,z)%n_b(j) = pop * nb
#ifndef NOSURFACTANT
            N(x,ny-1,z)%n_s(j) = pop * ns
#endif
#endif
          end if
          ! Now recolour the ones headed in the -Y direction.
   
          if ( is_fluid(N(x,ny,z)%rock_state) ) then 
            j = negy(i)
            pop = N(x,ny,z)%n_r(j)
#ifndef SINGLEFLUID
            pop = pop + N(x,ny,z)%n_b(j)
#ifndef NOSURFACTANT
            pop = pop + N(x,ny,z)%n_s(j)
#endif
#endif
            N(x,ny,z)%n_r(j) = pop * nr
#ifndef SINGLEFLUID
            N(x,ny,z)%n_b(j) = pop * nb
#ifndef NOSURFACTANT
            N(x,ny,z)%n_s(j) = pop * ns
#endif
#endif
          endif
          
        end do
      end do
    end do
  end if !maxy

  if (miny) then
    do z=1,nz
      do x=1,nx
        if ( is_fluid(N(x,3,z)%rock_state) ) then
          nr = sum(N(x,3,z)%n_r(:)*g)
          pop = nr
#ifndef SINGLEFLUID
          nb = sum(N(x,3,z)%n_b(:)*g)
          pop = pop + nb
#ifndef NOSURFACTANT
          ns = sum(N(x,3,z)%n_s(:)*g)
          pop = pop + ns   
#endif
#endif
        else
          pop = 0.0_rk
        end if
        ! nr,nb,ns end up containing the ratio of
        ! red to green to surf in the inward adjacent site.
        if (pop .eq. 0.0_rk) then
          ! Assign 1:1:1 ratio if adjacent site is empty.
          ! (or if it's a rock site).
#ifndef NOSURFACTANT
          nr = 1.0_rk/3.0_rk
          ns = 1.0_rk/3.0_rk
          nb = 1.0_rk/3.0_rk
#else
          nr = 1.0_rk/2.0_rk
#ifndef SINGLEFLUID
          nb = 1.0_rk/2.0_rk
#endif
#endif
        else
          nr = nr / pop
#ifndef SINGLEFLUID
          nb = nb / pop
#ifndef NOSURFACTANT
          ns = ns / pop
#endif
#endif
        end if
        do i=1,nnp
          ! Recolour the particles which are headed in the -Y direction.
          if ( is_fluid(N(x,2,z)%rock_state) ) then
            j = negy(i)
            pop = N(x,2,z)%n_r(j)
#ifndef SINGLEFLUID
            pop = pop + N(x,2,z)%n_b(j)
#ifndef NOSURFACTANT
            pop = pop +	N(x,2,z)%n_s(j)
#endif
#endif
            ! Note that particle number should be conserved
            ! give or take a bit of roundoff error.
            ! (but we're using 64-bit reals, so that doesn't matter.
            ! Right? Right? FIXME!)
            N(x,2,z)%n_r(j) = pop * nr
#ifndef SINGLEFLUID
            N(x,2,z)%n_b(j) = pop * nb
#ifndef NOSURFACTANT
            N(x,2,z)%n_s(j) = pop * ns
#endif
#endif
          end if
          ! Now recolour the ones headed in the +Y direction.
          if ( is_fluid(N(x,1,z)%rock_state) ) then 
            j = posy(i)
            pop = N(x,1,z)%n_r(j)
#ifndef SINGLEFLUID
            pop = pop + N(x,1,z)%n_b(j)
#ifndef NOSURFACTANT
            pop = pop + N(x,1,z)%n_s(j)
#endif
#endif

            N(x,1,z)%n_r(j) = pop * nr
#ifndef SINGLEFLUID
            N(x,1,z)%n_b(j) = pop * nb
#ifndef NOSURFACTANT
            N(x,1,z)%n_s(j) = pop * ns
#endif
#endif
          end if
        end do
      end do
    end do
  end if ! miny

  if (maxx) then
    do z=1,nz
      do y=1,ny
        if ( is_fluid(N(nx-2,y,z)%rock_state) ) then
          nr = sum(N(nx-2,y,z)%n_r(:)*g)
          pop = nr
#ifndef SINGLEFLUID
          nb = sum(N(nx-2,y,z)%n_b(:)*g)
          pop = pop + nb
#ifndef NOSURFACTANT
          ns = sum(N(nx-2,y,z)%n_s(:)*g)
          pop = pop + ns   
#endif
#endif
        else
          pop = 0.0_rk
        end if
        ! nr,nb,ns end up containing the ratio of
        ! red to green to surf in the inward adjacent site.
        if (pop .eq. 0.0_rk) then
          ! Assign 1:1:1 ratio if adjacent site is empty.
          ! (or if it's a rock site).
#ifndef NOSURFACTANT
          nr = 1.0_rk/3.0_rk
          ns = 1.0_rk/3.0_rk
          nb = 1.0_rk/3.0_rk
#else
          nr = 1.0_rk/2.0_rk
#ifndef SINGLEFLUID
          nb = 1.0_rk/2.0_rk
#endif
#endif
        else
          nr = nr / pop
#ifndef SINGLEFLUID
          nb = nb / pop
#ifndef NOSURFACTANT
          ns = ns / pop
#endif
#endif
        end if

        do i=1,nnp
          ! Recolour the particles which are headed in the +X
          ! direction.
          if ( is_fluid(N(nx-1,y,z)%rock_state) ) then
            j = posx(i)
            pop = N(nx-1,y,z)%n_r(j)
#ifndef SINGLEFLUID
            pop = pop + N(nx-1,y,z)%n_b(j)
#ifndef NOSURFACTANT
            pop = pop + N(nx-1,y,z)%n_s(j)
#endif
#endif
            N(nx-1,y,z)%n_r(j) = pop * nr
#ifndef SINGLEFLUID
            N(nx-1,y,z)%n_b(j) = pop * nb
#ifndef NOSURFACTANT
            N(nx-1,y,z)%n_s(j) = pop * ns
#endif
#endif
          end if
          ! Now recolour the ones headed in the -X direction.
          if ( is_fluid(N(nx,y,z)%rock_state) ) then 
            j = negx(i)
            pop = N(nx,y,z)%n_r(j)
#ifndef SINGLEFLUID
            pop = pop + N(nx,y,z)%n_b(j)
#endif
#ifndef NOSURFACTANT
            pop = pop + N(nx,y,z)%n_s(j)

#endif
            N(nx,y,z)%n_r(j) = pop * nr
#ifndef SINGLEFLUID
            N(nx,y,z)%n_b(j) = pop * nb
#ifndef NOSURFACTANT
            N(nx,y,z)%n_s(j) = pop * ns
#endif
#endif
          end if
        end do
      end do
    end do
  end if !maxy

  if (minx) then
    do z=1,nz
      do y=1,ny
        if ( is_fluid(N(3,y,z)%rock_state) ) then
          nr = sum(N(3,y,z)%n_r(:)*g)
          pop = nr
#ifndef SINGLEFLUID
          nb = sum(N(3,y,z)%n_b(:)*g)
          pop = pop + nb
#ifndef NOSURFACTANT
          ns = sum(N(3,y,z)%n_s(:)*g)
          pop = pop + ns  
#endif
#endif
        else
          pop = 0.0_rk
        end if
        ! nr,nb,ns end up containing the ratio of
        ! red to green to surf in the inward adjacent site.
        if (pop .eq. 0.0_rk) then
          ! Assign 1:1:1 ratio if adjacent site is empty.
          ! (or if it's a rock site).
#ifndef NOSURFACTANT
          nr = 1.0_rk/3.0_rk
          ns = 1.0_rk/3.0_rk
          nb = 1.0_rk/3.0_rk
#else
          nr = 1.0_rk/2.0_rk
#ifndef SINGLEFLUID
          nb = 1.0_rk/2.0_rk
#endif
#endif
        else
          nr = nr / pop
#ifndef SINGLEFLUID
          nb = nb / pop
#ifndef NOSURFACTANT
          ns = ns / pop
#endif
#endif
        end if
        do i=1,nnp
          ! Recolour the particles which are headed in the -X direction.
          if ( is_fluid(N(2,y,z)%rock_state) ) then
            j = negx(i)
            pop = N(2,y,z)%n_r(j)
#ifndef SINGLEFLUID
            pop = pop + N(2,y,z)%n_b(j)
#ifndef NOSURFACTANT
            pop = pop +	N(2,y,z)%n_s(j)
#endif
#endif
            ! Note that particle number should be conserved
            ! give or take a bit of roundoff error.
            ! (but we're using 64-bit reals, so that doesn't matter.
            ! Right? Right? FIXME!)
            N(2,y,z)%n_r(j) = pop * nr
#ifndef SINGLEFLUID
            N(2,y,z)%n_b(j) = pop * nb
#ifndef NOSURFACTANT
            N(2,y,z)%n_s(j) = pop * ns
#endif
#endif
          end if

          ! Now recolour the ones headed in the +X direction.
          if ( is_fluid(N(1,y,z)%rock_state) ) then 
            j = posx(i)
            pop = N(1,y,z)%n_r(j)
#ifndef SINGLEFLUID
            pop = pop + N(1,y,z)%n_b(j)
#ifndef NOSURFACTANT
            pop = pop + N(1,y,z)%n_s(j)
#endif
#endif
            N(1,y,z)%n_r(j) = pop * nr
#ifndef SINGLEFLUID
            N(1,y,z)%n_b(j) = pop * nb
#ifndef NOSURFACTANT
            N(1,y,z)%n_s(j) = pop * ns
#endif
#endif
          end if
        end do
      end do
    end do
  end if ! minx

end subroutine lbe_invade_sampled

!>At boundary regions, recolour all particles in the bottom boundary
!>region to be red=\c mr and blue=\c mb, and all particles in the top
!>boundary region to be red=\c mb and blue=\c mr.
!>
!> Added by Maddalena
!>
!>Only used for binary red/blue mixtures
!>
!>This only handles 2 layers of lattice sites in the Z boundaries.
subroutine lbe_invade_binary(N,mr,mb)
	implicit none
	type(lbe_site),dimension(0:,0:,0:) :: N
	real*8 :: mr,mb
        real*8, dimension(nvecs) :: g_inv

	integer :: x,y,z
	logical minz,maxz
	real*8	:: pop
	integer :: i

         g_inv = 1.d0/g


	! Determine whether I sit on any boundaries.

	minz=(start(3)==1)
	maxz=(start(3)>=(tnz-nz))

	if (minz) then
! Recolour bottom
		do y=1,ny
		 do x=1,nx
			if (N(x,y,1)%rock_state == 0) then
				do i=1,nvecs

					pop=	N(x,y,1)%n_r(i)
#ifndef SINGLEFLUID
					pop = pop + N(x,y,1)%n_b(i)
#endif

! Recolour with present densities
 				        N(x,y,1)%n_r(i)=mr*pop
#ifndef SINGLEFLUID
 				        N(x,y,1)%n_b(i)=mb*pop
#endif

				end do
			endif
			if (N(x,y,2)%rock_state == 0) then
				do i=1,nvecs
					pop=	N(x,y,2)%n_r(i)
#ifndef SINGLEFLUID
					pop = pop + N(x,y,2)%n_b(i)
#endif

! Recolour with present densities
 				         N(x,y,2)%n_r(i)=mr*pop
#ifndef SINGLEFLUID
 				         N(x,y,2)%n_b(i)=mb*pop
#endif

				end do
			endif
		 end do
		end do
	endif ! minz
 	if (maxz) then
! Recolour top
 		do y=1,ny
 		 do x=1,nx
 			if (N(x,y,nz-1)%rock_state == 0) then
 				do i=1,nvecs
 					pop=	N(x,y,nz-1)%n_r(i)
#ifndef SINGLEFLUID
 					pop = pop + N(x,y,nz-1)%n_b(i)
#endif
! Recolour with present densities
 				        N(x,y,nz-1)%n_r(i)=mb*pop
#ifndef SINGLEFLUID
 				        N(x,y,nz-1)%n_b(i)=mr*pop
#endif

 				end do
 			endif
 			if (N(x,y,nz)%rock_state == 0) then
 				do i=1,nvecs
 					pop=	N(x,y,nz)%n_r(i)
#ifndef SINGLEFLUID
 					pop = pop + N(x,y,nz)%n_b(i)
#endif

! Recolour with present densities
#ifndef SINGLEFLUID
 				        N(x,y,nz)%n_b(i)=mr*pop
#endif
 				        N(x,y,nz)%n_r(i)=mb*pop

 				end do
 			endif
 		 end do
 		end do
 	endif ! maxz
end subroutine lbe_invade_binary

!>At the bottom, recolour all particles the bottom boundary region to
!>be red=\c mr and blue=\c mb, in two lamellae of width frr,fbb
!>
!> Added by Jens, 12.08.09
!>
!>Only used for binary red/blue mixtures
!>
!>This only handles 2 layers of lattice sites in the Z boundaries.
subroutine lbe_invade_binary_lamellae(N,mr,mb,frr,frb)
        implicit none
        type(lbe_site),dimension(0:,0:,0:) :: N
	real*8 :: mr,mb,frr,frb
        real*8, dimension(nvecs) :: g_inv
        integer :: x,y,z
        logical minx,maxx,miny,maxy,minz,maxz
	real*8	:: pop
        integer :: i,j,k,l
        integer :: lamwidth     ! Width of one layer.
        integer :: ty   ! y-coordinate in global system.

        g_inv = 1.d0/g

        ! Determine whether I sit on any boundaries.

        minx=(start(1)==1)
        miny=(start(2)==1)
        minz=(start(3)==1)
        maxx=(start(1)>=(tnx-nx))
        maxy=(start(2)>=(tny-ny))
        maxz=(start(3)>=(tnz-nz))

        if (minz) then

          lamwidth = frr + frb
          if ((myrankc == 0).AND.(nt==0)) then
                if (mod(tny,lamwidth) .ne. 0) then
                        print*,'******* WARNING ******'
                        print*,'* Truncated lamellae *'
                        print*,'**********************'
                endif
          endif



! Recolour bottom
                 do y=1,ny
                   ty = y-1 + ccoords(2)*ny
                   l=mod(ty,lamwidth)
                  do x=1,nx
                   if (N(x,y,1)%rock_state == 0) then
                     do i=1,nvecs

                      pop = N(x,y,1)%n_r(i)
#ifndef SINGLEFLUID
                      pop = pop + N(x,y,1)%n_b(i)
#endif

! Recolour with present densities
                      if (l.lt.frr) then
                        N(x,y,1)%n_r(i)=mr*pop
#ifndef SINGLEFLUID
                        N(x,y,1)%n_b(i)=0.d0
                      else
                        N(x,y,1)%n_r(i)=0.d0
                        N(x,y,1)%n_b(i)=mb*pop
#endif
                      endif

                     end do
                   endif
                   if (N(x,y,2)%rock_state == 0) then
                    do i=1,nvecs
                      pop = N(x,y,2)%n_r(i)
#ifndef SINGLEFLUID
                      pop = pop + N(x,y,2)%n_b(i)
#endif

! Recolour with present densities
                      if (l.lt.frr) then
                        N(x,y,2)%n_r(i)=mr*pop
#ifndef SINGLEFLUID
                        N(x,y,2)%n_b(i)=0.d0
                      else
                        N(x,y,2)%n_r(i)=0.d0
                        N(x,y,2)%n_b(i)=mb*pop
#endif
                      endif
                    end do
                   endif
                 end do
                end do
           endif ! minz

end subroutine lbe_invade_binary_lamellae

!>  in and outflux with boundary conditions according to Zou and He, 
!>  Phys.Fluids 9, 1591, (1997), or, applied and formulated for D3Q19
!>  in Kutary et al. Computers and Geotechnics, 33, 381 (2006)
!>  Attention: They define the lattice vectors different from here.
!> 
!>  their notation  1 2 3 4 5 6 7 8  9  10 11 12 13 14 15 16 17 18 19
!>  our notation    1 2 3 4 5 6 7 11 12  8  9 10 14 13 15 16 18 17 19
!>
!>
!> There are two different variants:
!> itype = 0 is the original Zou and He boundary, and itype = 1 
!> includes the correct diagonal contributions
subroutine lbe_invade_zouandhe(N,inr,inb,outr,outb,itype)
  implicit none
  type(lbe_site), dimension(0:,0:,0:), intent(inout) :: N
  real*8, intent(in) :: inr, inb, outr, outb
  integer, intent(in) :: itype
  real*8 :: uz_r, uz_b, rho_r, rho_b, ux, uy
  real*8 :: F_r(19), F_b(19), fSCxr, fSCyr, fSCxb, fSCyb
  integer :: x, y, z
  logical minx, maxx, miny, maxy, minz, maxz
  integer :: i, j, k, tx, ty, tz
  real*8, dimension(nx,ny,nz,3) :: fSC_r
#ifndef SINGLEFLUID
  real*8, dimension(nx,ny,nz,3) :: fSC_b
#endif
  
  ! in the original paper flux was assumed to be perpendicular to the boundary
  ux = 0.
  uy = 0.
  
  ! Determine whether I sit on any boundaries.
  
  minx=(start(1)==1)
  miny=(start(2)==1)
  minz=(start(3)==1)
  maxx=(start(1)>=(tnx-nx))
  maxy=(start(2)>=(tny-ny))
  maxz=(start(3)>=(tnz-nz))

  if (SCMP) then
     if (start(3)==1) then
        do y = 0,(ny+1)
            do x = 0,(nx+1)
                N(x,y,0)%n_r(:) = N(x,y,2)%n_r(:)
            end do
        end do
     end if
     if (start(3)>=(tnz-nz)) then
        do y = 0,(ny+1)
           do x = 0,(nx+1)
              N(x,y,nz+1)%n_r(:) = N(x,y,nz-1)%n_r(:)
           end do
        end do
     end if
  end if
   

  if (minz) then
     do y=0,(ny+1)
        do x=0,(nx+1)
           if (N(x,y,1)%rock_state == 0) then             
#ifndef NOSURFACTANT
              ! not working
#endif
              F_r = N(x,y,1)%n_r(:)*g(:) 
              
              rho_r = inr
              
              if ( rho_r > 0.0 ) then
                 uz_r = 1.0d0 &
                      - ( F_r(1) + F_r(2) + F_r(3) + F_r(4) + F_r(7) + F_r(8) + F_r(11) + F_r(12) + F_r(19) &
                      + 2.0d0*( F_r(6) + F_r(10) + F_r(14) + F_r(16) + F_r(18) ) ) * amass_r/rho_r
              else 
                 uz_r = 0.0d0
              endif
              
              F_r(5)  = F_r(6) + rho_r*uz_r/amass_r/3.0d0
              if ( itype==0 ) then
                 F_r(15) = F_r(18) - ( F_r(3) - F_r(4) )/4.0d0 + rho_r*uz_r/amass_r/6.0d0
                 F_r(17) = F_r(16) + ( F_r(3) - F_r(4) )/4.0d0 + rho_r*uz_r/amass_r/6.0d0
                 F_r(9)  = F_r(14) - ( F_r(1) - F_r(2) )/4.0d0 + rho_r*uz_r/amass_r/6.0d0
                 F_r(13) = F_r(10) + ( F_r(1) - F_r(2) )/4.0d0 + rho_r*uz_r/amass_r/6.0d0
              else if (itype==1) then
                 F_r(9)  = uz_r*rho_r/amass_r/6.0d0 &
                      + F_r(14) + ( - F_r(1) + F_r(2) - F_r(7) - F_r(8) + F_r(11) + F_r(12) )/2.0d0
                 F_r(13) = uz_r*rho_r/amass_r/6.0d0 &
                      + F_r(10) + ( F_r(1) - F_r(2) + F_r(7) + F_r(8) - F_r(11) - F_r(12) )/2.0d0
                 F_r(15) = uz_r*rho_r/amass_r/6.0d0 &
                      + F_r(18) + ( - F_r(3) + F_r(4) - F_r(7) + F_r(8) - F_r(11) + F_r(12) )/2.0d0
                 F_r(17) = uz_r*rho_r/amass_r/6.0d0 &
                      + F_r(16) + ( F_r(3) - F_r(4) + F_r(7) - F_r(8) + F_r(11) - F_r(12) )/2.0d0 
              end if

              N(x,y,1)%n_r(:) = F_r(:)/g(:) 
#ifndef SINGLEFLUID
              F_r = N(x,y,1)%n_r(:)*g(:) 
              F_b = N(x,y,1)%n_b(:)*g(:) 

              rho_r = inr
              rho_b = inb

              fSCxr = 0.0d0
              fSCyr = 0.0d0
              fSCxb = 0.0d0
              fSCyb = 0.0d0

              F_r(5) = ( -2.0d0*F_r(10) - F_r(11) - F_r(12) - 2.0d0*F_r(14) - 2.0d0*F_r(16) &
                   - 2.0d0*F_r(18) - F_r(19) - F_r(1) - F_r(2) - F_r(3) - F_r(4) + F_r(6) &
                   - F_r(7) - F_r(8) )/3.0d0 + rho_r/3.0d0/amass_r

              F_r(9)  = ( -F_r(10) + F_r(11) + F_r(12) + 2.0d0*F_r(14) - F_r(16) - F_r(18) &
                   - F_r(19)/2.0d0 - 2.0d0*F_r(1) + F_r(2) - F_r(3)/2.0d0 - F_r(4)/2.0d0 &
                   - F_r(6) - 2.0d0*F_r(7) - 2.0d0*F_r(8) )/3.0d0 + rho_r/6.0d0/amass_r &
                   - fSCxr*tau_r/4.0d0/amass_r

              F_r(13) = ( 2.0d0*F_r(10) - 2.0d0*F_r(11) - 2.0d0*F_r(12) - F_r(14)  - F_r(16) &
                   - F_r(18) - F_r(19)/2.0d0 + F_r(1) - 2.0d0*F_r(2) - F_r(3)/2.0d0 &
                   - F_r(4)/2.0d0 - F_r(6) + F_r(7) + F_r(8) )/3.0d0 + rho_r/6.0d0/amass_r &
                   + fSCxr*tau_r/4.0d0/amass_r

              F_r(15) = ( -F_r(10) - 2.0d0*F_r(11) + F_r(12) - F_r(14) - F_r(16) + 2.0d0*F_r(18) &
                   - F_r(19)/2.0d0 - F_r(1)/2.0d0 - F_r(2)/2.0d0 - 2.0d0*F_r(3) + F_r(4) - F_r(6) &
                   - 2.0d0*F_r(7) + F_r(8) )/3.0d0 + rho_r/6.0d0/amass_r - fSCyr*tau_r/4.0d0/amass_r

              F_r(17) = ( -F_r(10) + F_r(11) - 2.0d0*F_r(12) - F_r(14) + 2.0d0*F_r(16) - F_r(18) &
                   - F_r(19)/2.0d0 - F_r(1)/2.0d0 - F_r(2)/2.0d0 + F_r(3) - 2.0d0*F_r(4) - F_r(6) &
                   + F_r(7) - 2.0d0*F_r(8) )/3.0d0 + rho_r/6.0d0/amass_r + fSCyr*tau_r/4.0d0/amass_r

  
              F_b(5) = ( -2.0d0*F_b(10) - F_b(11) - F_b(12) - 2.0d0*F_b(14) - 2.0d0*F_b(16) &
                   - 2.0d0*F_b(18) - F_b(19) - F_b(1) - F_b(2) - F_b(3) - F_b(4) + F_b(6) &
                   - F_b(7) - F_b(8) )/3.0d0 + rho_b/3.0d0/amass_b

              F_b(9)  = ( -F_b(10) + F_b(11) + F_b(12) + 2.0d0*F_b(14) - F_b(16) - F_b(18) &
                   - F_b(19)/2.0d0 - 2.0d0*F_b(1) + F_b(2) - F_b(3)/2.0d0 - F_b(4)/2.0d0 &
                   - F_b(6) - 2.0d0*F_b(7) - 2.0d0*F_b(8) )/3.0d0 + rho_b/6.0d0/amass_b &
                   - fSCxb*tau_b/4.0d0/amass_b

              F_b(13) = ( 2.0d0*F_b(10) - 2.0d0*F_b(11) - 2.0d0*F_b(12) - F_b(14) - F_b(16) &
                   - F_b(18) - F_b(19)/2.0d0 + F_b(1) - 2.0d0*F_b(2) - F_b(3)/2.0d0 &
                   - F_b(4)/2.0d0 - F_b(6) + F_b(7) + F_b(8) )/3.0d0 + rho_b/6.0d0/amass_r &
                   + fSCxb*tau_b/4.0d0/amass_b

              F_b(15) = ( -F_b(10) - 2.0d0*F_b(11) + F_b(12) - F_b(14) - F_b(16) + 2.0d0*F_b(18) &
                   - F_b(19)/2.0d0 - F_b(1)/2.0d0 - F_b(2)/2.0d0 - 2.0d0*F_b(3) + F_b(4) - F_b(6) &
                   - 2.0d0*F_b(7) + F_b(8) )/3.0d0 + rho_b/6.0d0/amass_b - fSCyb*tau_b/4.0d0/amass_b

              F_b(17) = ( -F_b(10) + F_b(11) - 2.0d0*F_b(12) - F_b(14) + 2.0d0*F_b(16) - F_b(18) &
                   - F_b(19)/2.0d0 - F_b(1)/2.0d0 - F_b(2)/2.0d0 + F_b(3) - 2.0d0*F_b(4) - F_b(6) &
                   + F_b(7) - 2.0d0*F_b(8) )/3.0d0 + rho_b/6.0d0/amass_b + fSCyb*tau_b/4.0d0/amass_b

              N(x,y,1)%n_r(:) = F_r(:)/g(:)
              N(x,y,1)%n_b(:) = F_b(:)/g(:) 
#endif
           endif
        end do
     end do
               

#ifdef SINGLEFLUID
     if (SCMP) then
        do y=0,(ny+1)
           do x=0,(nx+1)
              N(x,y,0)%n_r(:) = N(x,y,2)%n_r(:) 
              N(x,y,0)%rock_colour = N(x,y,2)%rock_colour
              N(x,y,0)%rock_state = N(x,y,2)%rock_state
           end do
        end do
        
        rho_r = inr
        
        do y=1,ny
           do x=1,nx
              if ( N(x,y,1)%rock_state == 0 ) then
                 CALL lbe_calculate_sc_forces(N,x,y,1,fSC_r)
                 
                 F_r = N(x,y,1)%n_r(:)*g(:)
                 
                 if ( rho_r > 0.0 ) then
                    uz_r = 1.0d0 &
                         - ( F_r(1) + F_r(2) + F_r(3) + F_r(4) + F_r(7) + F_r(8) + F_r(11) + F_r(12) + F_r(19) &
                         + 2.0d0*( F_r(6) + F_r(10) + F_r(14) + F_r(16) + F_r(18) ) )*amass_r/rho_r
                 else 
                    uz_r = 0.0d0
                 endif
                 
                 F_r(5)  = F_r(6) + rho_r*uz_r/amass_r/3.0d0
                 F_r(9)  = uz_r*rho_r/amass_r/6.0d0 - tau_r*fSC_r(x,y,1,1)/amass_r/4.0d0 &
                      + F_r(14) + ( - F_r(1) + F_r(2) - F_r(7) - F_r(8) + F_r(11) + F_r(12) )/2.0d0
                 F_r(13) = uz_r*rho_r/amass_r/6.0d0 + tau_r*fSC_r(x,y,1,1)/amass_r/4.0d0 &
                      + F_r(10) + ( F_r(1) - F_r(2) + F_r(7) + F_r(8) - F_r(11) - F_r(12) )/2.0d0
                 F_r(15) = uz_r*rho_r/amass_r/6.0d0 - tau_r*fSC_r(x,y,1,2)*amass_r/4.0d0 &
                      + F_r(18) + ( - F_r(3) + F_r(4) - F_r(7) + F_r(8) - F_r(11) + F_r(12) )/2.0d0
                 F_r(17) = uz_r*rho_r/amass_r/6.0d0 + tau_r*fSC_r(x,y,1,2)*amass_r/4.0d0 &
                      + F_r(16) + ( F_r(3) - F_r(4) + F_r(7) - F_r(8) + F_r(11) - F_r(12) )/2.0d0 
                 
                 N(x,y,1)%n_r(:) = F_r(:)/g(:) 
              endif
           end do
        end do
     endif
#endif
  endif ! minz
  
  if (maxz) then
     do y=0,(ny+1)
        do x=0,(nx+1)
           if (N(x,y,nz)%rock_state == 0) then
#ifndef NOSURFACTANT
              ! not working
#endif                
              F_r = N(x,y,nz)%n_r(:)*g(:) 
              
              rho_r = outr
              
              if ( rho_r > 0.0 ) then          
                 uz_r = -1.0d0 &
                      + ( F_r(1) + F_r(2) + F_r(3) + F_r(4) + F_r(7)  + F_r(8) + F_r(11) + F_r(12) + F_r(19) &
                      + 2.0d0*( F_r(5) + F_r(9) + F_r(13) + F_r(15) + F_r(17) ) )*amass_r/rho_r
              else 
                 uz_r = 0.0d0
              endif
              
              F_r(6) = F_r(5) - uz_r*rho_r/amass_r/3.0d0
              if ( itype==0 ) then
                 F_r(18) = F_r(15) + ( F_r(3) - F_r(4) )/4.0d0 - rho_r*uz_r/amass_r/6.0d0
                 F_r(16) = F_r(17) - ( F_r(3) - F_r(4) )/4.0d0 - rho_r*uz_r/amass_r/6.0d0
                 F_r(14) = F_r(9)  + ( F_r(1) - F_r(2) )/4.0d0 - rho_r*uz_r/amass_r/6.0d0
                 F_r(10) = F_r(13) - ( F_r(1) - F_r(2) )/4.0d0 - rho_r*uz_r/amass_r/6.0d0
              else if (itype==1) then
                 F_r(10) = - uz_r*rho_r/amass_r/6.0d0 &
                      + F_r(13) + ( - F_r(1) + F_r(2) - F_r(7) - F_r(8) + F_r(11) + F_r(12) )/2.0d0 
                 F_r(14) = - uz_r*rho_r/amass_r/6.0d0 &
                      + F_r(9) + ( F_r(1) - F_r(2) + F_r(7) + F_r(8) - F_r(11) - F_r(12) )/2.0d0
                 F_r(16) = - uz_r*rho_r/amass_r/6.0d0 &
                      + F_r(17) + ( - F_r(3) + F_r(4) - F_r(7) + F_r(8) - F_r(11) + F_r(12) )/2.0d0 
                 F_r(18) = - uz_r*rho_r/amass_r/6.0d0 &
                      + F_r(15) + ( F_r(3) - F_r(4) + F_r(7) - F_r(8) + F_r(11) - F_r(12) )/2.0d0 
              end if
              
              N(x,y,nz)%n_r(:) = F_r(:)/g(:) 
#ifndef SINGLEFLUID
              F_r = N(x,y,nz)%n_r(:)*g(:) 
              F_b = N(x,y,nz)%n_b(:)*g(:) 
              
              rho_r = outr
              rho_b = outb
              
              fSCxr = 0.0d0
              fSCyr = 0.0d0
              fSCxb = 0.0d0
              fSCyb = 0.0d0


              F_r(6) = ( -F_r(11) - F_r(12) - 2.0d0*F_r(13) - 2.0d0*F_r(15) - 2.0d0*F_r(17) &
                   - F_r(19) -  F_r(1) - F_r(2) - F_r(3) - F_r(4) + F_r(5) - F_r(7) - F_r(8) &
                   - 2.0d0*F_r(9) )/3.0d0 + rho_r/3.0d0/amass_r

              F_r(10)  = ( F_r(11) + F_r(12) + 2.0d0*F_r(13) - F_r(15) - F_r(17) - F_r(19)/2.0d0 &
                   - 2.0d0*F_r(1) +  F_r(2) - F_r(3)/2.0d0 - F_r(4)/2.0d0 - F_r(5) - 2.0d0*F_r(7) &
                   - 2.0d0*F_r(8) - F_r(9) )/3.0d0 + rho_r/6.0d0/amass_r - fSCxr*tau_r/4.0d0/amass_r

              F_r(14) = ( -2.0d0*F_r(11) - 2.0d0*F_r(12) - F_r(13) - F_r(15) - F_r(17) - F_r(19)/2.0d0 &
                   + F_r(1) - 2.0d0*F_r(2) - F_r(3)/2.0d0 - F_r(4)/2.0d0 - F_r(5) + F_r(7) + F_r(8) &
                   + 2.0d0*F_r(9) )/3.0d0 + rho_r/6.0d0/amass_r + fSCxr*tau_r/4.0d0/amass_r

              F_r(16) = ( -2.0d0*F_r(11) + F_r(12) - F_r(13) - F_r(15) + 2.0d0*F_r(17) - F_r(19)/2.0d0 &
                   - F_r(1)/2.0d0 - F_r(2)/2.0d0 - 2.0d0*F_r(3) + F_r(4) - F_r(5) - 2.0d0*F_r(7) + F_r(8) &
                   - F_r(9) )/3.0d0 + rho_r/6.0d0/amass_r - fSCyr*tau_r/4.0d0/amass_r

              F_r(18) = ( F_r(11) - 2.0d0*F_r(12) - F_r(13) + 2.0d0*F_r(15) - F_r(17) - F_r(19)/2.0d0 &
                   - F_r(1)/2.0d0 - F_r(2)/2.0d0 + F_r(3) - 2.0d0*F_r(4) - F_r(5) + F_r(7) - 2.0d0*F_r(8) &
                   - F_r(9) )/3.0d0 + rho_r/6.0d0/amass_r + fSCyr*tau_r/4.0d0/amass_r

  
              F_b(6) = ( -F_b(11) - F_b(12) - 2.0d0*F_b(13) - 2.0d0*F_b(15) - 2.0d0*F_b(17) &
                   - F_b(19) - F_b(1) - F_b(2) - F_b(3) - F_b(4) + F_b(5) - F_b(7) - F_b(8) &
                   - 2.0d0*F_b(9) )/3.0d0 + rho_b/3.0d0/amass_b

              F_b(10)  = ( F_b(11) + F_b(12) + 2.0d0*F_b(13) - F_b(15) - F_b(17) - F_b(19)/2.0d0 &
                   - 2.0d0*F_b(1) + F_b(2) - F_b(3)/2.0d0 - F_b(4)/2.0d0 - F_b(5) - 2.0d0*F_b(7) &
                   - 2.0d0*F_b(8) - F_b(9) )/3.0d0 + rho_b/6.0d0/amass_b - fSCxb*tau_b/4.0d0/amass_b

              F_b(14) = ( -2.0d0*F_b(11) - 2.0d0*F_b(12) - F_b(13) - F_b(15) - F_b(17) - F_b(19)/2.0d0 &
                   + F_b(1) - 2.0d0*F_b(2) - F_b(3)/2.0d0 - F_b(4)/2.0d0 - F_b(5) + F_b(7) + F_b(8) &
                   + 2.0d0*F_b(9) )/3.0d0 + rho_b/6.0d0/amass_b + fSCxb*tau_b/4.0d0/amass_b

              F_b(16) = ( -2.0d0*F_b(11) + F_b(12) - F_b(13) - F_b(15) + 2.0d0*F_b(17) - F_b(19)/2.0d0 &
                   - F_b(1)/2.0d0 - F_b(2)/2.0d0 - 2.0d0*F_b(3) + F_b(4) - F_b(5) - 2.0d0*F_b(7) + F_b(8) &
                   - F_b(9) )/3.0d0 + rho_b/6.0d0/amass_b - fSCyb*tau_b/4.0d0/amass_b

              F_b(18) = ( F_b(11) - 2.0d0*F_b(12) - F_b(13) + 2.0d0*F_b(15) - F_b(17) - F_b(19)/2.0d0 &
                   - F_b(1)/2.0d0 - F_b(2)/2.0d0 + F_b(3) - 2.0d0*F_b(4) - F_b(5) + F_b(7) - 2.0d0*F_b(8) &
                   - F_b(9) )/3.0d0 + rho_b/6.0d0/amass_b + fSCyb*tau_b/4.0d0/amass_b

              
              N(x,y,nz)%n_r(:) = F_r(:)/g(:) 
              N(x,y,nz)%n_b(:) = F_b(:)/g(:) 
#endif
           endif
        end do
     end do

#ifdef SINGLEFLUID          
     if (SCMP) then
        do y=0,(ny+1)
           do x=0,(nx+1)
              N(x,y,nz+1)%n_r(:) = N(x,y,nz-1)%n_r(:)
              N(x,y,nz+1)%rock_colour = N(x,y,nz-1)%rock_colour
              N(x,y,nz+1)%rock_state = N(x,y,nz-1)%rock_state
           end do
        end do
        
        rho_r = outr
                   
        do y=1,ny
           do x=1,nx
              if (N(x,y,nz)%rock_state == 0) then
                 CALL lbe_calculate_sc_forces(N,x,y,nz,fSC_r)
                 
                 F_r = N(x,y,nz)%n_r(:)*g(:)

                 if ( rho_r > 0.0 ) then          
                    uz_r = -1.0d0 &
                         + ( F_r(1) + F_r(2) + F_r(3) + F_r(4) + F_r(7)  + F_r(8) + F_r(11) + F_r(12) + F_r(19) &
                         + 2.0d0*( F_r(5) + F_r(9) + F_r(13) + F_r(15) + F_r(17) ) )*amass_r/rho_r
                 else 
                    uz_r = 0.0d0
                 endif

                 F_r(6) = F_r(5) - uz_r*rho_r/amass_r/3.0d0
                 F_r(10) = - uz_r*rho_r/amass_r/6.0d0 - tau_r*fSC_r(x,y,nz,1)/amass_r/4.0d0 &
                      + F_r(13) + ( - F_r(1) + F_r(2) - F_r(7) - F_r(8) + F_r(11) + F_r(12) )/2.0d0 
                 F_r(14) = - uz_r*rho_r/amass_r/6.0d0 + tau_r*fSC_r(x,y,nz,1)/amass_r/4.0d0 &
                      + F_r(9) + ( F_r(1) - F_r(2) + F_r(7) + F_r(8) - F_r(11) - F_r(12) )/2.0d0
                 F_r(16) = - uz_r*rho_r/amass_r/6.0d0 - tau_r*fSC_r(x,y,nz,2)/amass_r/4.0d0 &
                      + F_r(17) + ( - F_r(3) + F_r(4) - F_r(7) + F_r(8) - F_r(11) + F_r(12) )/2.0d0 
                 F_r(18) = - uz_r*rho_r/amass_r/6.0d0 + tau_r*fSC_r(x,y,nz,2)/amass_r/4.0d0 &
                      + F_r(15) + ( F_r(3) - F_r(4) + F_r(7) - F_r(8) + F_r(11) - F_r(12) )/2.0d0 
                 
                 N(x,y,nz)%n_r(:) = F_r(:)/g(:) 
              endif
           end do
        end do
     end if
#endif
  endif ! maxz

end subroutine lbe_invade_zouandhe


subroutine lbe_invade_zouandhemix(N,inr,ing,inb,outr,outg,outb)
	implicit none
	type(lbe_site),dimension(0:,0:,0:), intent(inout) :: N
        real*8, intent(in) :: inr, ing, inb, outr, outg, outb
	real*8 :: uz, rho, ux, uy
        real*8 :: F(19)
	integer :: x, y, z
	logical minx, maxx, miny, maxy, minz,  maxz
	real*8 :: pop
	integer :: i, j, k, tx, ty, tz
        real*8, dimension(nx,ny,nz,3) :: f_r
        
        ! in the original paper flux was assumed to be perpendicular to the boundary
        ux = 0.
        uy = 0.

	! Determine whether I sit on any boundaries.

	minx=(start(1)==1)
	miny=(start(2)==1)
	minz=(start(3)==1)
	maxx=(start(1)>=(tnx-nx))
	maxy=(start(2)>=(tny-ny))
	maxz=(start(3)>=(tnz-nz))

  if (SCMP) then
     if (start(3)==1) then
        do y = 0,(ny+1)
            do x = 0,(nx+1)
                N(x,y,0)%n_r(:) = N(x,y,2)%n_r(:)
            end do
        end do
     end if
     if (start(3)>=(tnz-nz)) then
        do y = 0,(ny+1)
           do x = 0,(nx+1)
              N(x,y,nz+1)%n_r(:) = N(x,y,nz-1)%n_r(:)
           end do
        end do
     end if
  end if

	if (minz) then
              do y=0,(ny+1)
		 do x=0,(nx+1)
			if (N(x,y,1)%rock_state == 0) then

! influx, original index system:
! uz = 1 - (F19 + F1+F2+F3+F4 +  F7+F8+F9+F10 +  2*(F6+F12+F13+F16+F17))/rho
!                           F5= F6+rho*uz/3.
!                           F15 = F17 - (F3-F4)/4. + rho*uz/6.
!                           F18 = F16 + (F3-F4)/4. + rho*uz/6.
!                           F11 = F13 - (F1-F2)/4. + rho*uz/6.
!                           F14 = F12 + (F1-F2)/4. + rho*uz/6.

#ifndef NOSURFACTANT
                        F = N(x,y,1)%n_s(:)*g(:) 
                      
                        rho = ing
                        if(rho>0.) then
                        uz = 1. - (F(19) + &
                                       F(1) + F(2) + F(3) + F(4) +  &
                                       F(7) + F(11) + F(12) + F(8) +  &
                                       2.*( F(6) + F(10) + F(14) + F(16) + F(18)))/rho
                        else 
                          uz = 0.
                        endif 
                                   
                        F(5)  = F(6) + rho*uz/3.                    
                        F(15) = F(18) - (F(7)+F(11)+F(3)-F(8)-F(12)-F(4))/2. &
                             + uz*rho/6. + uy*rho/2.
                        F(17) = F(16) - (F(8)+F(12)+F(4)-F(7)-F(11)-F(3))/2. &
                             + uz*rho/6. - uy*rho/2.
                        F(9)  = F(14) - (F(7)+F(8)+F(1)-F(11)-F(12)-F(2))/2. &
                             + uz*rho/6. + ux*rho/2.
                        F(13) = F(10) - (F(11)+F(12)+F(2)-F(7)-F(8)-F(1))/2. &
                             + uz*rho/6. - ux*rho/2.

                        N(x,y,1)%n_s(:) = F(:)/g(:) 
#endif

                        F = N(x,y,1)%n_r(:)*g(:) 
                       
                        rho = inr
                        ux = 0.0d0
                        uy = 0.0d0
                        uz = ing

                        if ( ( start(1) + x - 1.0d0 - tnx/2.0d0)**2 &
                             + ( start(2) + y - 1.0d0 - tny/2.0d0)**2 .le. inb**2 ) then 
!!!   density BC
                            if(rho>0.) then
                               uz = 1. - (F(19) + &
                                    F(1) + F(2) + F(3) + F(4) +  &
                                    F(7) + F(11) + F(12) + F(8) +  &
                                    2.*( F(6) + F(10) + F(14) + F(16) + F(18)))/rho
                            else 
                               uz = 0.
                            endif
                           
                            F(5)  = F(6) + rho*uz/3.
                           
                               F(15) = F(18) - (F(7)+F(11)+F(3)-F(8)-F(12)-F(4))/2. &
                                    + uz*rho/6. + uy*rho/2.
                               F(17) = F(16) - (F(8)+F(12)+F(4)-F(7)-F(11)-F(3))/2. &
                                    + uz*rho/6. - uy*rho/2.
                               F(9)  = F(14) - (F(7)+F(8)+F(1)-F(11)-F(12)-F(2))/2. &
                                    + uz*rho/6. + ux*rho/2.
                               F(13) = F(10) - (F(11)+F(12)+F(2)-F(7)-F(8)-F(1))/2. &
                                    + uz*rho/6. - ux*rho/2.
                         else
!!!  velocity  BC uz, ux=0, uy=0                                                    
                           rho = - ( F(1) + F(2) + F(3) + F(4) + F(7) + F(8) + F(11) + F(12) + F(19) &
                                + 2.0d0*( F(6) + F(10) + F(14)  + F(16) + F(18) ) )/(uz - 1.0d0)

                              F(5) = uz*rho/3.0d0 + F(6)
                              F(9) = (  3.0d0*ux + uz)*rho/6.0d0 & 
                                   + F(14) + ( - F(1) + F(2) - F(7) - F(8) + F(11) + F(12) )/2.0d0
                              F(13) = (- 3.0d0*ux + uz)*rho/6.0d0 &
                                   + F(10) + ( F(1) - F(2) + F(7) + F(8) - F(11) - F(12) )/2.0d0
                              F(15) = (  3.0d0*uy + uz)*rho/6.0d0 &
                                   + F(18) + ( - F(3) + F(4) - F(7) + F(8) - F(11) + F(12) )/2.0d0
                              F(17) = (- 3.0d0*uy + uz)*rho/0.6D1 &
                                   + F(16) + ( F(3) - F(4) + F(7) - F(8) + F(11) - F(12) )/2.0d0
                        end if
                        N(x,y,1)%n_r(:) = F(:)/g(:)
#ifndef SINGLEFLUID

                        F = N(x,y,1)%n_b(:)*g(:) 
                        
                        rho = inb
                        if(rho>0.) then
                        uz = 1. - (F(19) + &
                                       F(1) + F(2) + F(3) + F(4) +  &
                                       F(7) + F(11) + F(12) + F(8) +  &
                                       2.*( F(6) + F(10) + F(14) + F(16) + F(18)))/rho
                        else 
                          uz = 0.
                        endif                                       
                        F(5)  = F(6) + rho*uz/3.
                        
                          F(15) = F(18) - (F(7)+F(11)+F(3)-F(8)-F(12)-F(4))/2. &
                                  + uz*rho/6. + uy*rho/2.
                          F(17) = F(16) - (F(8)+F(12)+F(4)-F(7)-F(11)-F(3))/2. &
                                  + uz*rho/6. - uy*rho/2.
                          F(9)  = F(14) - (F(7)+F(8)+F(1)-F(11)-F(12)-F(2))/2. &
                                  + uz*rho/6. + ux*rho/2.
                          F(13) = F(10) - (F(11)+F(12)+F(2)-F(7)-F(8)-F(1))/2. &
                                  + uz*rho/6. - ux*rho/2.

                        N(x,y,1)%n_b(:) = F(:)/g(:) 
#endif
                     endif
                  end do
               end do
        
#ifdef SINGLEFLUID
               if (SCMP) then
                  do y=0,(ny+1)
                     do x=0,(nx+1)
                        N(x,y,0)%n_r(:) = N(x,y,2)%n_r(:) 
                     end do
                  end do

                  do y=1,ny
                     do x=1,nx

                        rho = inr
                        ux = 0.0d0
                        uy = 0.0d0
                        uz = ing

			if (N(x,y,1)%rock_state == 0) then
                            CALL lbe_calculate_sc_forces(N,x,y,1,f_r)

                            F = N(x,y,1)%n_r(:)*g(:)

                            if ( ( start(1) + x - 1.0d0 - tnx/2.0d0)**2 &                              
                                 + ( start(2) + y - 1.0d0 - tny/2.0d0)**2 .le. inb**2 ) then 
!!   density BC                                                           
                                 uz = 1.0d0 - ( F(1) + F(2) + F(3) + F(4) + F(7) + F(8) + F(11) + F(12) + F(19) &
                                      + 2.0d0*( F(6) + F(10) + F(14) + F(16) + F(18) ) )/rho
                               
                                 F(5) = uz*rho/3.0d0 + F(6)
                                 F(9) = uz*rho/6.0d0 + f_r(x,y,1,1)*rho/4.0d0 &
                                      + F(14) + ( - F(1) + F(2) - F(7) - F(8) + F(11) + F(12) )/2.0d0
                                 F(13) = uz*rho/6.0d0 - f_r(x,y,1,1)*rho/4.0d0 &
                                      + F(10) + ( F(1) - F(2) + F(7) + F(8) - F(11) - F(12) )/2.0d0
                                 F(15) = uz*rho/6.0d0 + f_r(x,y,1,2)*rho/4.0d0 &
                                      + F(18) + ( - F(3) + F(4) - F(7) + F(8) - F(11) + F(12) )/2.0d0
                                 F(17) = uz*rho/6.0d0 - f_r(x,y,1,2)*rho/4.0d0 &
                                      + F(16) + ( F(3) - F(4) + F(7) - F(8) + F(11) - F(12) )/2.0d0 
                             else
!!!  velocity  BC uz, ux=0, uy=0                                                        S                           
                              rho = - ( F(1) + F(2) + F(3) + F(4) + F(7) + F(8) + F(11) + F(12) + F(19) &
                                   + 2.0d0*( F(6) + F(10) + F(14)  + F(16) + F(18) ) )/(uz - 1.0d0)
                         
                              F(5) = uz*rho/3.0d0 + F(6)
                              F(9) = (   6.0d0*ux + 2.0d0*uz + 3.0d0*f_r(x,y,1,1) )*rho/12.0d0 & 
                                   + F(14) + ( - F(1) + F(2) - F(7) - F(8) + F(11) + F(12) )/2.0d0
                              F(13) = ( - 6.0d0*ux + 2.0d0*uz - 3.0d0*f_r(x,y,1,1) )*rho/12.0d0 &
                                   + F(10) + ( F(1) - F(2) + F(7) + F(8) - F(11) - F(12) )/2.0d0
                              F(15) = (  6.0d0*uy + 2.0d0*uz + 3.0d0*f_r(x,y,1,2) )*rho/12.0d0 &
                                   + F(18) + ( - F(3) + F(4) - F(7) + F(8) - F(11) + F(12) )/2.0d0
                              F(17) = ( - 6.0d0*uy + 2.0d0*uz - 3.0d0*f_r(x,y,1,2) )*rho/12.0d0 &
                                   + F(16) + ( F(3) - F(4) + F(7) - F(8) + F(11) - F(12) )/2.0d0
                            end if
                            N(x,y,1)%n_r(:) = F(:)/g(:)
                        endif
                      end do
                   end do
                endif
#endif
             endif ! minz
             
            if (maxz) then
            do y=0,(ny+1)
               do x=0,(nx+1)
                  if (N(x,y,nz)%rock_state == 0) then

#ifndef NOSURFACTANT
                     F = N(x,y,nz)%n_s(:)*g(:) 
                     
                     ux = 0.0d0
                     uy = 0.0d0
                     uz = outg
                            
                     rho = ( F(1) + F(2) + F(3) + F(4) + F(7) + F(8) + F(11) + F(12) + F(19) + &
                          2.0d0*( F(5) + F(9) + F(13) + F(15) + F(17) ) )/(uz + 1.0d0)

                        F(6) = - uz*rho/3.0d0 + F(5)
                        F(10) = (  3.0d0*ux - uz)*rho/6.0d0 & 
                             + F(13) + ( - F(1) + F(2) - F(7) - F(8) + F(11) + F(12) )/2.0d0
                        F(14) = (- 3.0d0*ux - uz)*rho/6.0d0 &
                             + F(9) + (  F(1) - F(2) + F(7) + F(8) - F(11) - F(12) )/2.0d0                
                        F(16) = (  3.0d0*uy - uz)*rho/6.0d0 &
                             + F(17) + ( - F(3) + F(4) - F(7) + F(8) - F(11) + F(12) )/2.0d0
                        F(18) = (- 3.0d0*uy - uz)*rho/0.6D1 &
                             + F(15) + ( F(3) - F(4) + F(7) - F(8) + F(11) - F(12) )/2.0d0

                     N(x,y,nz)%n_s(:) = F(:)/g(:) 
#endif
                    
                     F = N(x,y,nz)%n_r(:)*g(:) 
                     
                     ux = 0.0d0
                     uy = 0.0d0
                     uz = outr
                     
                     rho = ( F(1) + F(2) + F(3) + F(4) + F(7) + F(8) + F(11) + F(12) + F(19) &
                          + 2.0d0*( F(5) + F(9) + F(13) + F(15) + F(17) ) )/(uz + 1.0d0)

                        F(6) = - uz*rho/3.0d0 + F(5)
                        F(10) = (  3.0d0*ux - uz)*rho/6.0d0 & 
                            + F(13) + ( - F(1) + F(2) - F(7) - F(8) + F(11) + F(12) )/2.0d0
                        F(14) = (- 3.0d0*ux - uz)*rho/6.0d0 &
                             + F(9) + ( F(1) - F(2) + F(7) + F(8) - F(11) - F(12) )/2.0d0                
                        F(16) = (  3.0d0*uy - uz)*rho/6.0d0 &
                             + F(17) + ( - F(3) + F(4) - F(7) + F(8) - F(11) + F(12) )/2.0d0
                        F(18) = (- 3.0d0*uy - uz)*rho/0.6D1 &
                             + F(15) + ( F(3) - F(4) + F(7) - F(8) + F(11) - F(12) )/2.0d0
                     
                     N(x,y,nz)%n_r(:) = F(:)/g(:) 
#ifndef SINGLEFLUID
                        
                     F = N(x,y,nz)%n_b(:)*g(:) 
                     
                     ux = 0.0d0
                     uy = 0.0d0
                     uz = outb

                     rho = ( F(1) + F(2) + F(3) + F(4) + F(7) + F(8) + F(11) + F(12) + F(19) &
                          + 2.0d0*( F(5) + F(9) + F(13) + F(15) + F(17) ) )/(uz + 1.0d0)

                     F(6) = F(5) - uz*rho/3.0d0  
                     F(10) = ( 3.0d0*ux- uz)*rho/6.0d0 & 
                          + F(13) + ( - F(1) + F(2) - F(7) - F(8) + F(11) + F(12) )/2.0d0
                     F(14) = (- 3.0d0*ux- uz)*rho/6.0d0 &
                          + F(9) + ( F(1) - F(2) + F(7) + F(8) - F(11) - F(12) )/2.0d0                
                     F(16) = (3.0d0*uy- uz)*rho/6.0d0 &
                          + F(17) + ( - F(3) + F(4) - F(7) + F(8) - F(11) + F(12) )/2.0d0
                     F(18) = (- 3.0d0*uy- uz)*rho/0.6D1 &
                          + F(15) + ( F(3) - F(4) + F(7) - F(8) + F(11) - F(12) )/2.0d0
                       
                     N(x,y,nz)%n_b(:) = F(:)/g(:) 
#endif
                    endif
                 end do
                end do

#ifdef SINGLEFLUID            
                if (SCMP) then
                   
                   do y=0,(ny+1)
                      do x=0,(nx+1)
                         N(x,y,nz+1)%n_r(:) = N(x,y,nz-1)%n_r(:) 
                      end do
                   end do
                                    
                   ux = 0.0d0
                   uy = 0.0d0
                   uz = outr

                   do y=1,ny
                      do x=1,nx
                         if (N(x,y,nz)%rock_state == 0) then
                            
                            CALL lbe_calculate_sc_forces(N,x,y,nz,f_r)
                   
!                            ux = sum(N(x,y,nz-1)%n_r(:)*cx(:)*g(:))/sum(N(x,y,nz-1)%n_r(:)*g(:))
!                            uy = sum(N(x,y,nz-1)%n_r(:)*cy(:)*g(:))/sum(N(x,y,nz-1)%n_r(:)*g(:))
!                            uz = sum(N(x,y,nz-1)%n_r(:)*cz(:)*g(:))/sum(N(x,y,nz-1)%n_r(:)*g(:))
                         
                            F = N(x,y,nz)%n_r(:)*g(:)
                            
                            rho = ( F(1) + F(2) + F(3) + F(4) + F(7) + F(8) + F(11) + F(12) + F(19) &
                                 + 2.0d0*( F(5) + F(9) + F(13) + F(15) + F(17) ) )/(uz + 1.0d0)
                            
                            F(6) = - uz*rho/3.0d0 + F(5)
                            F(10) = rho*(  6.0d0*ux - 2.0d0*uz + 3.0d0*f_r(x,y,nz,1) )/12.0d0 & 
                                 + F(13) + ( - F(1) + F(2) - F(7) - F(8) + F(11) + F(12) )/2.0d0
                            F(14) = rho*(- 6.0d0*ux - 2.0d0*uz - 3.0d0*f_r(x,y,nz,1) )/12.0d0 &
                                 + F(9) + (  F(1) - F(2) + F(7) + F(8) - F(11) - F(12) )/2.0d0
                            F(16) = rho*(  6.0d0*uy - 2.0d0*uz + 3.0d0*f_r(x,y,nz,2) )/12.0d0 &
                                 + F(17) + ( - F(3) + F(4) - F(7) + F(8) - F(11) + F(12) )/2.0d0
                            F(18) = rho*(- 6.0d0*uy - 2.0d0*uz - 3.0d0*f_r(x,y,nz,2) )/12.0d0 &
                                 + F(15) + ( F(3) - F(4) + F(7) - F(8) + F(11) - F(12) )/2.0d0

                            N(x,y,nz)%n_r(:) = F(:)/g(:)
                         endif
                      end do
                   end do
                end if
#endif
          endif ! maxz

end subroutine lbe_invade_zouandhemix


!> the same as lbe_invade_zouandhe, except that we set u_z instead of
!> the density to a given value
!> itype = 0 is the original Zou and He boundary, and itype = 1
!> includes the correct diagonal contributions
!> negative itype values denote that u is the mass flow and not the velocity
subroutine lbe_invade_zouandheflux(N,uz_r,uz_g,uz_b,itype)
	implicit none
	type(lbe_site),dimension(0:,0:,0:), intent(inout)  :: N
        real*8, intent(in) :: uz_r,uz_g,uz_b
        integer, intent(in) :: itype
	real*8 :: rho, uz, ux, uy
        real*8 :: F(19)
	integer :: x,y,z
	logical minx,maxx,miny,maxy,minz,maxz
	real*8	:: pop
	integer :: i,j,k
        
        ! in the original paper flux was assumed to be perpendicular to the boundary
        ux = 0.
        uy = 0.

	! Determine whether I sit on any boundaries.

	minx=(start(1)==1)
	miny=(start(2)==1)
	minz=(start(3)==1)
	maxx=(start(1)>=(tnx-nx))
	maxy=(start(2)>=(tny-ny))
	maxz=(start(3)>=(tnz-nz))

	if (minz) then
		do y=1,ny
		 do x=1,nx
			if (N(x,y,1)%rock_state == 0) then


#ifndef NOSURFACTANT
                        F = N(x,y,1)%n_s(:)*g(:) 
                        
                        uz = uz_g
                        
                        if(itype<0) then
                          rho =  F(19) + F(1) + F(2) + F(3) + F(4) +  &
                                  F(7) + F(11) + F(12) + F(8) +  &
                                  2.*( F(6) + F(10) + F(14) + F(16) + F(18)) + uz
                          if (rho>0.) then
                          uz = uz / rho
                          else
                            uz = 0.
                          endif
                        else
                          rho = (F(19) + F(1) + F(2) + F(3) + F(4) +  &
                                  F(7) + F(11) + F(12) + F(8) +  &
                                  2.*( F(6) + F(10) + F(14) + F(16) + F(18)))/(1.-uz)
                        end if

                        F(5)  = F(6) + uz*rho/3.
                        
                        if (itype==0) then
                           F(15) = F(18) - (F(3)-F(4))/4. + rho*uz/6.
                           F(17) = F(16) + (F(3)-F(4))/4. + rho*uz/6.
                           F(9)  = F(14) - (F(1)-F(2))/4. + rho*uz/6.
                           F(13) = F(10) + (F(1)-F(2))/4. + rho*uz/6.
                        else if ((itype==1).or.(itype==-1)) then
                          F(15) = F(18) - (F(7)+F(11)+F(3)-F(8)-F(12)-F(4))/2. &
                                  + uz*rho/6. + uy*rho/2.
                          F(17) = F(16) - (F(8)+F(12)+F(4)-F(7)-F(11)-F(3))/2. &
                                  + uz*rho/6. - uy*rho/2.
                          F(9)  = F(14) - (F(7)+F(8)+F(1)-F(11)-F(12)-F(2))/2. &
                                  + uz*rho/6. + ux*rho/2.
                          F(13) = F(10) - (F(11)+F(12)+F(2)-F(7)-F(8)-F(1))/2. &
                                  + uz*rho/6. - ux*rho/2.
                        end if

                        N(x,y,1)%n_s(:) = F(:)/g(:) 
#endif
                        F = N(x,y,1)%n_r(:)*g(:) 
                        
                        uz = uz_r

                        if(itype<0) then
                          rho =  F(19) + F(1) + F(2) + F(3) + F(4) +  &
                                  F(7) + F(11) + F(12) + F(8) +  &
                                  2.*( F(6) + F(10) + F(14) + F(16) + F(18)) + uz
                          if (rho>0.) then
                          uz = uz / rho
                          else
                            uz = 0.
                          endif
                        else
                          rho = (F(19) + F(1) + F(2) + F(3) + F(4) +  &
                                  F(7) + F(11) + F(12) + F(8) +  &
                                  2.*( F(6) + F(10) + F(14) + F(16) + F(18)))/(1.-uz)
                        end if

                        F(5)  = F(6) + uz*rho/3.
                        
                        if (itype==0) then
                           F(15) = F(18) - (F(3)-F(4))/4. + rho*uz/6.
                           F(17) = F(16) + (F(3)-F(4))/4. + rho*uz/6.
                           F(9)  = F(14) - (F(1)-F(2))/4. + rho*uz/6.
                           F(13) = F(10) + (F(1)-F(2))/4. + rho*uz/6.
                        else if ((itype==1).or.(itype==-1)) then
                          F(15) = F(18) - (F(7)+F(11)+F(3)-F(8)-F(12)-F(4))/2. &
                                  + uz*rho/6. + uy*rho/2.
                          F(17) = F(16) - (F(8)+F(12)+F(4)-F(7)-F(11)-F(3))/2. &
                                  + uz*rho/6. - uy*rho/2.
                          F(9)  = F(14) - (F(7)+F(8)+F(1)-F(11)-F(12)-F(2))/2. &
                                  + uz*rho/6. + ux*rho/2.
                          F(13) = F(10) - (F(11)+F(12)+F(2)-F(7)-F(8)-F(1))/2. &
                                  + uz*rho/6. - ux*rho/2.
                        end if

                        N(x,y,1)%n_r(:) = F(:)/g(:) 
#ifndef SINGLEFLUID
                        F = N(x,y,1)%n_b(:)*g(:) 
                        
                        uz = uz_b

                        if(itype<0) then
                          rho =  F(19) + F(1) + F(2) + F(3) + F(4) +  &
                                  F(7) + F(11) + F(12) + F(8) +  &
                                  2.*( F(6) + F(10) + F(14) + F(16) + F(18)) + uz
                          if (rho>0.) then
                          uz = uz / rho
                          else
                            uz = 0.
                          endif
                        else
                          rho = (F(19) + F(1) + F(2) + F(3) + F(4) +  &
                                  F(7) + F(11) + F(12) + F(8) +  &
                                  2.*( F(6) + F(10) + F(14) + F(16) + F(18)))/(1.-uz)
                        end if

                        F(5)  = F(6) + uz*rho/3.
                        
                        if (itype==0) then
                           F(15) = F(18) - (F(3)-F(4))/4. + rho*uz/6.
                           F(17) = F(16) + (F(3)-F(4))/4. + rho*uz/6.
                           F(9)  = F(14) - (F(1)-F(2))/4. + rho*uz/6.
                           F(13) = F(10) + (F(1)-F(2))/4. + rho*uz/6.
                        else if ((itype==1).or.(itype==-1)) then
                          F(15) = F(18) - (F(7)+F(11)+F(3)-F(8)-F(12)-F(4))/2. &
                                  + uz*rho/6. + uy*rho/2.
                          F(17) = F(16) - (F(8)+F(12)+F(4)-F(7)-F(11)-F(3))/2. &
                                  + uz*rho/6. - uy*rho/2.
                          F(9)  = F(14) - (F(7)+F(8)+F(1)-F(11)-F(12)-F(2))/2. &
                                  + uz*rho/6. + ux*rho/2.
                          F(13) = F(10) - (F(11)+F(12)+F(2)-F(7)-F(8)-F(1))/2. &
                                  + uz*rho/6. - ux*rho/2.
                        end if

                        N(x,y,1)%n_b(:) = F(:)/g(:) 
#endif
			endif
		 end do
		end do
	endif ! minz

	if (maxz) then
		do y=1,ny
		 do x=1,nx
                        if (N(x,y,nz)%rock_state == 0) then
#ifndef NOSURFACTANT
                        F = N(x,y,nz)%n_s(:)*g(:) 
                        
                        uz = uz_g

                        if(itype<0) then
                          rho =  F(19) + F(1) + F(2) + F(3) + F(4) +  &
                                 F(7) + F(11) + F(12) + F(8) +  &
                                 2.*( F(5) + F(9) + F(13) + F(15) + F(17)) - uz
                          if (rho>0.) then
                          uz = uz / rho
                          else
                            uz = 0.
                          endif
                        else
                          rho = (F(19) + F(1) + F(2) + F(3) + F(4) +  &
                                 F(7) + F(11) + F(12) + F(8) +  &
                                 2.*( F(5) + F(9) + F(13) + F(15) + F(17)))/(uz+1.)
                        end if

                        F(6) = F(5)-uz*rho/3.
                        
                        if (itype==0) then
                          F(18) = F(15) + (F(3)-F(4))/4. - rho*uz/6.
                          F(16) = F(17) - (F(3)-F(4))/4. - rho*uz/6.
                          F(14) = F(9)  + (F(1)-F(2))/4. - rho*uz/6.
                          F(10) = F(13) - (F(1)-F(2))/4. - rho*uz/6.
                        else if ((itype==1).or.(itype==-1)) then
                          F(18) = F(15) + (F(7)+F(11)+F(3)-F(8)-F(12)-F(4))/2. &
                                  - uz*rho/6. - uy*rho/2.
                          F(16) = F(17) + (F(8)+F(12)+F(4)-F(7)-F(11)-F(3))/2. &
                                  - uz*rho/6. + uy*rho/2.
                          F(14)  = F(9) + (F(7)+F(8)+F(1)-F(11)-F(12)-F(2))/2. &
                                  - uz*rho/6. - ux*rho/2.
                          F(10) = F(13) + (F(11)+F(12)+F(2)-F(7)-F(8)-F(1))/2. &
                                  - uz*rho/6. + ux*rho/2.
                        end if

                        N(x,y,nz)%n_s(:) = F(:)/g(:) 
#endif
                        F = N(x,y,nz)%n_r(:)*g(:) 
                        
                        uz = uz_r

                        if(itype<0) then
                          rho =  F(19) + F(1) + F(2) + F(3) + F(4) +  &
                                 F(7) + F(11) + F(12) + F(8) +  &
                                 2.*( F(5) + F(9) + F(13) + F(15) + F(17)) - uz
                          if (rho>0.) then
                          uz = uz / rho
                          else
                            uz = 0.
                          endif
                        else
                          rho = (F(19) + F(1) + F(2) + F(3) + F(4) +  &
                                 F(7) + F(11) + F(12) + F(8) +  &
                                 2.*( F(5) + F(9) + F(13) + F(15) + F(17)))/(uz+1.)
                        end if

                        F(6) = F(5)-uz*rho/3.
                        
                        if (itype==0) then
                          F(18) = F(15) + (F(3)-F(4))/4. - rho*uz/6.
                          F(16) = F(17) - (F(3)-F(4))/4. - rho*uz/6.
                          F(14) = F(9)  + (F(1)-F(2))/4. - rho*uz/6.
                          F(10) = F(13) - (F(1)-F(2))/4. - rho*uz/6.
                        else if ((itype==1).or.(itype==-1)) then
                          F(18) = F(15) + (F(7)+F(11)+F(3)-F(8)-F(12)-F(4))/2. &
                                  - uz*rho/6. - uy*rho/2.
                          F(16) = F(17) + (F(8)+F(12)+F(4)-F(7)-F(11)-F(3))/2. &
                                  - uz*rho/6. + uy*rho/2.
                          F(14)  = F(9) + (F(7)+F(8)+F(1)-F(11)-F(12)-F(2))/2. &
                                  - uz*rho/6. - ux*rho/2.
                          F(10) = F(13) + (F(11)+F(12)+F(2)-F(7)-F(8)-F(1))/2. &
                                  - uz*rho/6. + ux*rho/2.
                        end if

                        N(x,y,nz)%n_r(:) = F(:)/g(:) 
#ifndef SINGLEFLUID
                        F = N(x,y,nz)%n_b(:)*g(:) 
                        
                        uz = uz_b

                        if(itype<0) then
                          rho =  F(19) + F(1) + F(2) + F(3) + F(4) +  &
                                 F(7) + F(11) + F(12) + F(8) +  &
                                 2.*( F(5) + F(9) + F(13) + F(15) + F(17)) - uz
                          if (rho>0.) then
                          uz = uz / rho
                          else
                            uz = 0.
                          endif
                        else
                          rho = (F(19) + F(1) + F(2) + F(3) + F(4) +  &
                                 F(7) + F(11) + F(12) + F(8) +  &
                                 2.*( F(5) + F(9) + F(13) + F(15) + F(17)))/(uz+1.)
                        end if

                        F(6) = F(5)-uz*rho/3.
                        
                        if (itype==0) then
                          F(18) = F(15) + (F(3)-F(4))/4. - rho*uz/6.
                          F(16) = F(17) - (F(3)-F(4))/4. - rho*uz/6.
                          F(14) = F(9)  + (F(1)-F(2))/4. - rho*uz/6.
                          F(10) = F(13) - (F(1)-F(2))/4. - rho*uz/6.
                        else if ((itype==1).or.(itype==-1)) then
                          F(18) = F(15) + (F(7)+F(11)+F(3)-F(8)-F(12)-F(4))/2. &
                                  - uz*rho/6. - uy*rho/2.
                          F(16) = F(17) + (F(8)+F(12)+F(4)-F(7)-F(11)-F(3))/2. &
                                  - uz*rho/6. + uy*rho/2.
                          F(14)  = F(9) + (F(7)+F(8)+F(1)-F(11)-F(12)-F(2))/2. &
                                  - uz*rho/6. - ux*rho/2.
                          F(10) = F(13) + (F(11)+F(12)+F(2)-F(7)-F(8)-F(1))/2. &
                                  - uz*rho/6. + ux*rho/2.
                        end if

                        N(x,y,nz)%n_b(:) = F(:)/g(:) 
#endif
                        endif
                 end do
                end do
        endif ! maxz


end subroutine lbe_invade_zouandheflux

subroutine lbe_invade_after_evaporation(N,m_evp,m_evp_gr,m_evp_gb,m_evp_freq_f,m_evp_freq_a,in,out,m_evp_set_density)
  implicit none
  type(lbe_site),dimension(0:,0:,0:), intent(inout) :: N
  real*8, intent(in) :: m_evp,m_evp_gr,m_evp_gb,m_evp_freq_f,m_evp_freq_a
  real*8 m_evp_use
  logical, intent(in) :: in(3), out(3)
  real*8 :: tmp_r, tmp_b, tmp_s
  real*8 :: tmp_r1, tmp_b1, tmp_bz
  integer :: x, y, z
  real*8 :: F(19)
  !real*8,dimension(:,:,:,:),allocatable :: F_r, F_b
  real*8, dimension(nx,ny,nz,3) :: F_r, F_b,F_s,F_r2, F_b2,F_s2
  logical :: ismin(3), ismax(3)
  integer :: iin(3), iout(3)
  logical, intent(in) :: m_evp_set_density
  real(kind=rk),dimension(1:3) :: momentum,vel

  ismin(1)=(start(1)==1)
  ismin(2)=(start(2)==1)
  ismin(3)=(start(3)==1)
  ismax(1)=(start(1)>=(tnx-nx))
  ismax(2)=(start(2)>=(tny-ny))
  ismax(3)=(start(3)>=(tnz-nz))
  iin(1) = 0
  iin(2) = 0
  iin(3) = 0
  iout(1) = 0
  iout(2) = 0
  iout(3) = 0
  if (.NOT.m_evp_set_density) then
  if (in(1)) iin(1)=1
  if (in(2)) iin(2)=1
  if (in(3)) iin(3)=1
  if (out(1)) iout(1)=1
  if (out(2)) iout(2)=1
  if (out(3)) iout(3)=1
  endif
  call boltz_dist(0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, F(:)) 

  m_evp_use=m_evp+m_evp_freq_a*DSIN(m_evp_freq_f*nt)
#ifndef SINGLEFLUID

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Set boundary density evaporation method
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (m_evp_set_density) then
!x
  if ( ismin(1) .and. in(1) ) then ! min x
     do y=1,ny
        do z=1,nz
           if (N(1,y,z)%rock_state == 0) then
              tmp_r = sum( N(1,y,z)%n_r(:)*g(:) )
              tmp_b = sum( N(1,y,z)%n_b(:)*g(:) )
              N(1,y,z)%n_r(:) = m_evp_use/amass_r * F(:)
              N(1,y,z)%n_b(:) = (tmp_b+tmp_r-m_evp_use)/amass_b * F(:)
           endif
        end do
     end do
  endif 

  if ( ismax(1) .and. out(1) ) then ! max x
     do y=1,ny
        do z=1,nz
           if (N(nx,y,z)%rock_state == 0) then
              tmp_r = sum( N(nx,y,z)%n_r(:)*g(:) )
              tmp_b = sum( N(nx,y,z)%n_b(:)*g(:) )
              N(nx,y,z)%n_r(:) = m_evp_use/amass_r * F(:)
              N(nx,y,z)%n_b(:) = (tmp_b+tmp_r-m_evp_use)/amass_b * F(:)
           endif
        end do
     end do
  endif 

!y
  if ( ismin(2) .and. in(2) ) then ! min y
     do x=(1+iin(1)),(nx-iout(1))
        do z=1,nz
           if (N(x,1,z)%rock_state == 0) then
              tmp_r = sum( N(x,1,z)%n_r(:)*g(:) )
              tmp_b = sum( N(x,1,z)%n_b(:)*g(:) )
              N(x,1,z)%n_r(:) = m_evp_use/amass_r * F(:)
              N(x,1,z)%n_b(:) = (tmp_b+tmp_r-m_evp_use)/amass_b * F(:)
           endif
        end do
     end do
  endif 

  if ( ismax(2) .and. out(2) ) then ! max y
     do x=(1+iin(1)),(nx-iout(1))
        do z=1,nz
           if (N(x,ny,z)%rock_state == 0) then
              tmp_r = sum( N(x,ny,z)%n_r(:)*g(:) )
              tmp_b = sum( N(x,ny,z)%n_b(:)*g(:) )
              N(x,ny,z)%n_r(:) = m_evp_use/amass_r * F(:)
              N(x,ny,z)%n_b(:) = (tmp_b+tmp_r-m_evp_use)/amass_b * F(:)
           endif
        end do
     end do
  endif 

!z
  if ( ismin(3) .and. in(3) ) then ! min z
     do x=(1+iin(1)),(nx-iout(1))
        do y=(1+iin(2)),(ny-iout(2))
           if (N(x,y,1)%rock_state == 0) then
              tmp_r = sum( N(x,y,1)%n_r(:)*g(:) )
              tmp_b = sum( N(x,y,1)%n_b(:)*g(:) )
              N(x,y,1)%n_r(:) = m_evp_use/amass_r * F(:)
              N(x,y,1)%n_b(:) = (tmp_b+tmp_r-m_evp_use)/amass_b * F(:)
           endif
        end do
     end do
  endif 


 if ( ismax(3) .and. out(3) ) then ! max z
     do x=(1+iin(1)),(nx-iout(1))
        do y=(1+iin(2)),(ny-iout(2))
           if (N(x,y,nz)%rock_state == 0) then
              tmp_r = sum( N(x,y,nz)%n_r(:)*g(:) )
              tmp_b = sum( N(x,y,nz)%n_b(:)*g(:) )
              !N(x,y,nz)%n_r(:) = m_evp_use/amass_r * F(:)
              !N(x,y,nz)%n_b(:) = (tmp_b+tmp_r-m_evp_use)/amass_b * F(:)

!! Test new boundary model with suitable velocity pu=-1/2F
! calculate shan-chen force
#ifdef SINGLEFLUID
              call lbe_calculate_sc_forces(N, x, y, nz, F_r)
              call lbe_calculate_sc_forces(N, x, y, nz-1, F_r2)
#else
#ifdef NOSURFACTANT
              call lbe_calculate_sc_forces(N, x, y, nz, F_b, F_r)
              call lbe_calculate_sc_forces(N, x, y, nz-1, F_b2, F_r2)
#else
              call lbe_calculate_sc_forces(N, x, y, nz, F_b, F_r, F_s)
              call lbe_calculate_sc_forces(N, x, y, nz-1, F_b2, F_r2, F_s2)
#endif 
#endif
 ! if (x==2 .AND. y==2 ) then
 !print*, "nz \n", nz
 !print*, "mom \n",   velocity(1,0), velocity(2,0),velocity(3,0)
 !print*, "force_r: \n",  F_r(x,y,nz,1), F_r(x,y,nz,2),F_r(x,y,nz,3)
 !print*, "force_b: \n",  F_b(x,y,nz,1), F_b(x,y,nz,2),F_b(x,y,nz,3)
 !print*, "momentum 1 2 x , y, z: \n",rho,  momentum(1),
 !momentum(2),momentum(3), p_r(1), p_r(2),p_r(3)
!do s = 1, nvecs 
!print*, s, local_distribution%n_r(s)*g(s), local_distribution%n_r_pre(s)*g(s),
!cx(s), cy(s), cz(s)
!end do
! READ*
!endif
          momentum(1)=amass_r * sum(N(x,y,nz)%n_r(:) * g(:) * c(:,1))
          momentum(2)=amass_r * sum(N(x,y,nz)%n_r(:) * g(:) * c(:,2))
          momentum(3)=amass_r * sum(N(x,y,nz)%n_r(:) * g(:) * c(:,3))

vel(:) = (momentum(:)-0.5_rk*F_r(x,y,nz,:)-0.5_rk*F_r2(x,y,nz-1,:)) / tmp_r
call boltz_dist(vel(1), vel(2), vel(3), 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, F(:))
N(x,y,nz)%n_r(:) = tmp_r/amass_r * F(:)
            momentum(1)=amass_b * sum(N(x,y,nz)%n_b(:) * g(:) * c(:,1))
            momentum(2)=amass_b * sum(N(x,y,nz)%n_b(:) * g(:) * c(:,2))
            momentum(3)=amass_b * sum(N(x,y,nz)%n_b(:) * g(:) * c(:,3))

vel(:) = (momentum(:)-0.5_rk*F_b(x,y,nz,:)-0.5_rk*F_b2(x,y,nz-1,:)) / tmp_b
  call boltz_dist(vel(1), vel(2), vel(3), 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,0.0d0, F(:))
N(x,y,nz)%n_b(:) = tmp_b/amass_b * F(:)

 
           endif
        end do
     end do
  endif 

endif
#endif
!^endif of #ifndef SINGLEFLUID

end subroutine lbe_invade_after_evaporation



subroutine lbe_invade_evaporation(N,m_evp,m_evp_gr,m_evp_gb,m_evp_freq_f,m_evp_freq_a,in,out,m_evp_set_density)
  implicit none
  type(lbe_site),dimension(0:,0:,0:), intent(inout) :: N
  real*8, intent(in) :: m_evp,m_evp_gr,m_evp_gb,m_evp_freq_f,m_evp_freq_a
  real*8 m_evp_use, m_evp_sin
  logical, intent(in) :: in(3), out(3)
  real*8 :: tmp_r, tmp_b, tmp_s
  real*8 :: tmp_r1, tmp_b1, tmp_bz
  integer :: x, y, z
  real*8 :: F(19)
  !real*8,dimension(:,:,:,:),allocatable :: F_r, F_b
  real*8, dimension(nx,ny,nz,3) :: F_r, F_b,F_s,F_r2, F_b2,F_s2
  logical :: ismin(3), ismax(3)
  integer :: iin(3), iout(3)
  logical, intent(in) :: m_evp_set_density
  real(kind=rk),dimension(1:3) :: momentum, momentum_r, momentum_b,velr,velb

  ismin(1)=(start(1)==1)
  ismin(2)=(start(2)==1)
  ismin(3)=(start(3)==1)
  ismax(1)=(start(1)>=(tnx-nx))
  ismax(2)=(start(2)>=(tny-ny))
  ismax(3)=(start(3)>=(tnz-nz))
  iin(1) = 0
  iin(2) = 0
  iin(3) = 0
  iout(1) = 0
  iout(2) = 0
  iout(3) = 0
  !if (in(1)) iin(1)=1
  !if (in(2)) iin(2)=1
  !if (in(3)) iin(3)=1
  !if (out(1)) iout(1)=1
  !if (out(2)) iout(2)=1
  !if (out(3)) iout(3)=1

  call boltz_dist(0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, F(:)) 

  !m_evp_use=m_evp+m_evp_freq_a*SIND(m_evp_freq_f*nt) !SIND is actually surprisingly incompatible, hence the manual version
  m_evp_use=m_evp+m_evp_freq_a*SIN(m_evp_freq_f*nt*0.0174532925199432957692369076848861271344287188854172)

#ifndef SINGLEFLUID

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Set boundary density evaporation method
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (m_evp_set_density) then
!x
  if ( ismin(1) .and. in(1) ) then ! min x
     do y=1,ny
        do z=1,nz
           if (N(1,y,z)%rock_state == 0) then
              tmp_r = sum( N(1,y,z)%n_r(:)*g(:) )
              tmp_b = sum( N(1,y,z)%n_b(:)*g(:) )
              N(1,y,z)%n_r(:) = m_evp_use/amass_r * F(:)
              N(1,y,z)%n_b(:) = (tmp_b+tmp_r-m_evp_use)/amass_b * F(:)
           endif
        end do
     end do
  endif 

  if ( ismax(1) .and. out(1) ) then ! max x
     do y=1,ny
        do z=1,nz
           if (N(nx,y,z)%rock_state == 0) then
              tmp_r = sum( N(nx,y,z)%n_r(:)*g(:) )
              tmp_b = sum( N(nx,y,z)%n_b(:)*g(:) )
              N(nx,y,z)%n_r(:) = m_evp_use/amass_r * F(:)
              N(nx,y,z)%n_b(:) = (tmp_b+tmp_r-m_evp_use)/amass_b * F(:)
           endif
        end do
     end do
  endif 

!y
  if ( ismin(2) .and. in(2) ) then ! min y
     do x=(1+iin(1)),(nx-iout(1))
        do z=1,nz
           if (N(x,1,z)%rock_state == 0) then
              tmp_r = sum( N(x,1,z)%n_r(:)*g(:) )
              tmp_b = sum( N(x,1,z)%n_b(:)*g(:) )
              N(x,1,z)%n_r(:) = m_evp_use/amass_r * F(:)
              N(x,1,z)%n_b(:) = (tmp_b+tmp_r-m_evp_use)/amass_b * F(:)
           endif
        end do
     end do
  endif 

  if ( ismax(2) .and. out(2) ) then ! max y
     do x=(1+iin(1)),(nx-iout(1))
        do z=1,nz
           if (N(x,ny,z)%rock_state == 0) then
              tmp_r = sum( N(x,ny,z)%n_r(:)*g(:) )
              tmp_b = sum( N(x,ny,z)%n_b(:)*g(:) )
              N(x,ny,z)%n_r(:) = m_evp_use/amass_r * F(:)
              N(x,ny,z)%n_b(:) = (tmp_b+tmp_r-m_evp_use)/amass_b * F(:)
           endif
        end do
     end do
  endif 

!z
  if ( ismin(3) .and. in(3) ) then ! min z
     do x=(1+iin(1)),(nx-iout(1))
        do y=(1+iin(2)),(ny-iout(2))
           if (N(x,y,1)%rock_state == 0) then
              tmp_r = sum( N(x,y,1)%n_r(:)*g(:) )
              tmp_b = sum( N(x,y,1)%n_b(:)*g(:) )
              N(x,y,1)%n_r(:) = m_evp_use/amass_r * F(:)
              N(x,y,1)%n_b(:) = (tmp_b+tmp_r-m_evp_use)/amass_b * F(:)
           endif
        end do
     end do
  endif 



  if ( ismax(3) .and. out(3) ) then ! max z
     do x=(1+iin(1)),(nx-iout(1))
        do y=(1+iin(2)),(ny-iout(2))
           if (N(x,y,nz)%rock_state == 0) then
              tmp_r = sum( N(x,y,nz)%n_r(:)*g(:) )
              tmp_b = sum( N(x,y,nz)%n_b(:)*g(:) )
              N(x,y,nz)%n_r(:) = m_evp_use/amass_r * F(:)
              N(x,y,nz)%n_b(:) = (tmp_b+tmp_r-m_evp_use)/amass_b * F(:)
!		
! Test slow start up
!       if (nt < 1000000) then
!				m_evp_sin = 0.0129609228-(0.0129609228-m_evp_use)*(nt-500000)/500000 
!		 else
!				m_evp_sin = m_evp_use
!		 endif 
!              N(x,y,nz)%n_r(:) = m_evp_sin/amass_r * F(:)
!              N(x,y,nz)%n_b(:) = (tmp_b+tmp_r-m_evp_sin)/amass_b * F(:)


!! Test new boundary model with suitable velocity pu=-1/2F
! calculate shan-chen force
!        momentum(1)=amass_r * sum(N(x,y,nz)%n_r(:) * g(:) * c(:,1))
!        momentum(2)=amass_r * sum(N(x,y,nz)%n_r(:) * g(:) * c(:,2))
!        momentum(3)=amass_r * sum(N(x,y,nz)%n_r(:) * g(:) * c(:,3))

!#ifdef SINGLEFLUID
!              call lbe_calculate_sc_forces(N, x, y, nz, F_r)
!              call lbe_calculate_sc_forces(N, x, y, nz-1, F_r2)
!#else
!#ifdef NOSURFACTANT
!              call lbe_calculate_sc_forces(N, x, y, nz, F_b, F_r)
!              call lbe_calculate_sc_forces(N, x, y, nz-1, F_b2, F_r2)
!#else
!              call lbe_calculate_sc_forces(N, x, y, nz, F_b, F_r, F_s)
!              call lbe_calculate_sc_forces(N, x, y, nz-1, F_b2, F_r2, F_s2)
!#endif 
!#endif
 ! if (x==2 .AND. y==2 ) then
 !print*, "nz \n", nz
 !print*, "mom \n",   velocity(1,0), velocity(2,0),velocity(3,0)
 !print*, "force_r: \n",  F_r(x,y,nz,1), F_r(x,y,nz,2),F_r(x,y,nz,3)
 !print*, "force_b: \n",  F_b(x,y,nz,1), F_b(x,y,nz,2),F_b(x,y,nz,3)
 !print*, "momentum 1 2 x , y, z: \n",rho,  momentum(1),
 !momentum(2),momentum(3), p_r(1), p_r(2),p_r(3)
!do s = 1, nvecs 
!print*, s, local_distribution%n_r(s)*g(s), local_distribution%n_r_pre(s)*g(s),
!cx(s), cy(s), cz(s)
!end do
! READ*
!endif

!vel(:) = (momentum(:)+0.5_rk*F_r(x,y,nz,:)+0.5_rk*F_r2(x,y,nz-1,:)) / m_evp_use
!velr(:) = momentum(:)/tmp_r
!momentum_b(:) = momentum(:) * (tmp_r - m_evp_use)/tmp_r 
!call boltz_dist(velr(1), velr(2), velr(3), 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, F(:))
!N(x,y,nz)%n_r(:) = m_evp_use/amass_r * F(:)
!          momentum(1)=amass_b * sum(N(x,y,nz)%n_b(:) * g(:) * c(:,1))
!          momentum(2)=amass_b * sum(N(x,y,nz)%n_b(:) * g(:) * c(:,2))
!          momentum(3)=amass_b * sum(N(x,y,nz)%n_b(:) * g(:) * c(:,3))
!vel(:) = (momentum(:) +0.5_rk*F_b(x,y,nz,:)+0.5_rk*F_b2(x,y,nz-1,:)) / (tmp_b+tmp_r-m_evp_use)
!velb(:) = (momentum(:) + momentum_b(:))/(tmp_b+tmp_r-m_evp_use)
!  call boltz_dist(velb(1), velb(2), velb(3), 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,0.0d0, F(:))
!N(x,y,nz)%n_b(:) = (tmp_b+tmp_r-m_evp_use)/amass_b * F(:)

           endif
        end do
     end do
  endif 

!And now as we are done with the red blue conversion, the same again with green to red and then green to blue
!m_evp_gr:
!x
#ifndef NOSURFACTANT
  if ( ismin(1) .and. in(1) ) then ! min x
     do y=1,ny
        do z=1,nz
           if (N(1,y,z)%rock_state == 0) then
              tmp_s = sum( N(1,y,z)%n_s(:)*g(:) )
              tmp_r = sum( N(1,y,z)%n_r(:)*g(:) )
              N(1,y,z)%n_s(:) = m_evp_gr/amass_s * F(:)
              N(1,y,z)%n_r(:) = (tmp_r+tmp_s-m_evp_gr)/amass_r * F(:)
           endif
        end do
     end do
  endif 

  if ( ismax(1) .and. out(1) ) then ! max x
     do y=1,ny
        do z=1,nz
           if (N(nx,y,z)%rock_state == 0) then
              tmp_s = sum( N(nx,y,z)%n_s(:)*g(:) )
              tmp_r = sum( N(nx,y,z)%n_r(:)*g(:) )
              N(nx,y,z)%n_s(:) = m_evp_gr/amass_s * F(:)
              N(nx,y,z)%n_r(:) = (tmp_r+tmp_s-m_evp_gr)/amass_r * F(:)
           endif
        end do
     end do
  endif 

!y
  if ( ismin(2) .and. in(2) ) then ! min y
     do x=(1+iin(1)),(nx-iout(1))
        do z=1,nz
           if (N(x,1,z)%rock_state == 0) then
              tmp_s = sum( N(x,1,z)%n_s(:)*g(:) )
              tmp_r = sum( N(x,1,z)%n_r(:)*g(:) )
              N(x,1,z)%n_s(:) = m_evp_gr/amass_s * F(:)
              N(x,1,z)%n_r(:) = (tmp_r+tmp_s-m_evp_gr)/amass_r * F(:)
           endif
        end do
     end do
  endif 

  if ( ismax(2) .and. out(2) ) then ! max y
     do x=(1+iin(1)),(nx-iout(1))
        do z=1,nz
           if (N(x,ny,z)%rock_state == 0) then
              tmp_s = sum( N(x,ny,z)%n_s(:)*g(:) )
              tmp_r = sum( N(x,ny,z)%n_r(:)*g(:) )
              N(x,ny,z)%n_s(:) = m_evp_gr/amass_s * F(:)
              N(x,ny,z)%n_r(:) = (tmp_r+tmp_s-m_evp_gr)/amass_r * F(:)
           endif
        end do
     end do
  endif 

!z
  if ( ismin(3) .and. in(3) ) then ! min z
     do x=(1+iin(1)),(nx-iout(1))
        do y=(1+iin(2)),(ny-iout(2))
           if (N(x,y,1)%rock_state == 0) then
              tmp_s = sum( N(x,y,1)%n_s(:)*g(:) )
              tmp_r = sum( N(x,y,1)%n_r(:)*g(:) )
              N(x,y,1)%n_s(:) = m_evp_gr/amass_s * F(:)
              N(x,y,1)%n_r(:) = (tmp_r+tmp_s-m_evp_gr)/amass_r * F(:)
           endif
        end do
     end do
  endif 

  if ( ismax(3) .and. out(3) ) then ! max z
     do x=(1+iin(1)),(nx-iout(1))
        do y=(1+iin(2)),(ny-iout(2))
           if (N(x,y,nz)%rock_state == 0) then
              tmp_s = sum( N(x,y,nz)%n_s(:)*g(:) )
              tmp_r = sum( N(x,y,nz)%n_r(:)*g(:) )
              N(x,y,nz)%n_s(:) = m_evp_gr/amass_s * F(:)
              N(x,y,nz)%n_r(:) = (tmp_r+tmp_s-m_evp_gr)/amass_r * F(:)
           endif
        end do
     end do
  endif 

!And now the evaporation green to blue
!m_evp_gb:
!x
  if ( ismin(1) .and. in(1) ) then ! min x
     do y=1,ny
        do z=1,nz
           if (N(1,y,z)%rock_state == 0) then
              tmp_s = sum( N(1,y,z)%n_s(:)*g(:) )
              tmp_b = sum( N(1,y,z)%n_b(:)*g(:) )
              N(1,y,z)%n_s(:) = m_evp_gb/amass_s * F(:)
              N(1,y,z)%n_b(:) = (tmp_b+tmp_s-m_evp_gb)/amass_b * F(:)
           endif
        end do
     end do
  endif 

  if ( ismax(1) .and. out(1) ) then ! max x
     do y=1,ny
        do z=1,nz
           if (N(nx,y,z)%rock_state == 0) then
              tmp_s = sum( N(nx,y,z)%n_s(:)*g(:) )
              tmp_b = sum( N(nx,y,z)%n_b(:)*g(:) )
              N(nx,y,z)%n_s(:) = m_evp_gb/amass_s * F(:)
              N(nx,y,z)%n_b(:) = (tmp_b+tmp_s-m_evp_gb)/amass_b * F(:)
           endif
        end do
     end do
  endif 

!y
  if ( ismin(2) .and. in(2) ) then ! min y
     do x=(1+iin(1)),(nx-iout(1))
        do z=1,nz
           if (N(x,1,z)%rock_state == 0) then
              tmp_s = sum( N(x,1,z)%n_s(:)*g(:) )
              tmp_b = sum( N(x,1,z)%n_b(:)*g(:) )
              N(x,1,z)%n_s(:) = m_evp_gb/amass_s * F(:)
              N(x,1,z)%n_b(:) = (tmp_b+tmp_s-m_evp_gb)/amass_b * F(:)
           endif
        end do
     end do
  endif 

  if ( ismax(2) .and. out(2) ) then ! max y
     do x=(1+iin(1)),(nx-iout(1))
        do z=1,nz
           if (N(x,ny,z)%rock_state == 0) then
              tmp_s = sum( N(x,ny,z)%n_s(:)*g(:) )
              tmp_b = sum( N(x,ny,z)%n_b(:)*g(:) )
              N(x,ny,z)%n_s(:) = m_evp_gb/amass_s * F(:)
              N(x,ny,z)%n_b(:) = (tmp_b+tmp_s-m_evp_gb)/amass_b * F(:)
           endif
        end do
     end do
  endif 

!z
  if ( ismin(3) .and. in(3) ) then ! min z
     do x=(1+iin(1)),(nx-iout(1))
        do y=(1+iin(2)),(ny-iout(2))
           if (N(x,y,1)%rock_state == 0) then
              tmp_s = sum( N(x,y,1)%n_s(:)*g(:) )
              tmp_b = sum( N(x,y,1)%n_b(:)*g(:) )
              N(x,y,1)%n_s(:) = m_evp_gb/amass_s * F(:)
              N(x,y,1)%n_b(:) = (tmp_b+tmp_s-m_evp_gb)/amass_b * F(:)
           endif
        end do
     end do
  endif 

  if ( ismax(3) .and. out(3) ) then ! max z
     do x=(1+iin(1)),(nx-iout(1))
        do y=(1+iin(2)),(ny-iout(2))
           if (N(x,y,nz)%rock_state == 0) then
              tmp_s = sum( N(x,y,nz)%n_s(:)*g(:) )
              tmp_b = sum( N(x,y,nz)%n_b(:)*g(:) )
              N(x,y,nz)%n_s(:) = m_evp_gb/amass_s * F(:)
              N(x,y,nz)%n_b(:) = (tmp_b+tmp_s-m_evp_gb)/amass_b * F(:)
           endif
        end do
     end do
  endif 
#endif
!^endif of #ifndef NOSURFACTANT




else
!above is the else of if (m_evp_set_density) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Relative evaporation method
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!x
  if ( ismin(1) .and. in(1) ) then ! min x
     do y=1,ny
        do z=1,nz
           if (N(1,y,z)%rock_state == 0) then
              tmp_r = sum( N(1,y,z)%n_r(:)*g(:) ) - m_evp_use/amass_r
              tmp_b = sum( N(1,y,z)%n_b(:)*g(:) ) + m_evp_use/amass_b              
              if ( tmp_r .gt. 0.0d0 .and. tmp_b .gt. 0.0d0 ) then
                 N(1,y,z)%n_r(:) = tmp_r * F(:)
                 N(1,y,z)%n_b(:) = tmp_b * F(:)
              else
                 N(1,y,z)%n_r(:) = sum( N(1,y,z)%n_r(:)*g(:) ) * F(:)
                 N(1,y,z)%n_b(:) = sum( N(1,y,z)%n_b(:)*g(:) ) * F(:)
              end if
           endif
        end do
     end do
  endif 

  if ( ismax(1) .and. out(1) ) then ! max x
     do y=1,ny
        do z=1,nz
           if (N(nx,y,z)%rock_state == 0) then
              tmp_r = sum( N(nx,y,z)%n_r(:)*g(:) ) - m_evp_use/amass_r
              tmp_b = sum( N(nx,y,z)%n_b(:)*g(:) ) + m_evp_use/amass_b                 
              if ( tmp_r .gt. 0.0d0 .and. tmp_b .gt. 0.0d0 ) then
                 N(nx,y,z)%n_r(:) = tmp_r * F(:)
                 N(nx,y,z)%n_b(:) = tmp_b * F(:)
              else
                 N(nx,y,z)%n_r(:) = sum( N(nx,y,z)%n_r(:)*g(:) ) * F(:)
                 N(nx,y,z)%n_b(:) = sum( N(nx,y,z)%n_b(:)*g(:) ) * F(:)
              end if
           endif
        end do
     end do
  endif 

!y
  if ( ismin(2) .and. in(2) ) then ! min y
     do x=(1+iin(1)),(nx-iout(1))
        do z=1,nz
           if (N(x,1,z)%rock_state == 0) then
              tmp_r = sum( N(x,1,z)%n_r(:)*g(:) ) - m_evp_use/amass_r
              tmp_b = sum( N(x,1,z)%n_b(:)*g(:) ) + m_evp_use/amass_b              
              if ( tmp_r .gt. 0.0d0 .and. tmp_b .gt. 0.0d0 ) then
                 N(x,1,z)%n_r(:) = tmp_r * F(:)
                 N(x,1,z)%n_b(:) = tmp_b * F(:)
              else
                 N(x,1,z)%n_r(:) = sum( N(x,1,z)%n_r(:)*g(:) ) * F(:)
                 N(x,1,z)%n_b(:) = sum( N(x,1,z)%n_b(:)*g(:) ) * F(:)
              end if
           endif
        end do
     end do
  endif 

  if ( ismax(2) .and. out(2) ) then ! max y
     do x=(1+iin(1)),(nx-iout(1))
        do z=1,nz
           if (N(x,ny,z)%rock_state == 0) then
              tmp_r = sum( N(x,ny,z)%n_r(:)*g(:) ) - m_evp_use/amass_r
              tmp_b = sum( N(x,ny,z)%n_b(:)*g(:) ) + m_evp_use/amass_b                 
              if ( tmp_r .gt. 0.0d0 .and. tmp_b .gt. 0.0d0 ) then
                 N(x,ny,z)%n_r(:) = tmp_r * F(:)
                 N(x,ny,z)%n_b(:) = tmp_b * F(:)
              else
                 N(x,ny,z)%n_r(:) = sum( N(x,ny,z)%n_r(:)*g(:) ) * F(:)
                 N(x,ny,z)%n_b(:) = sum( N(x,ny,z)%n_b(:)*g(:) ) * F(:)
              end if
           endif
        end do
     end do
  endif 

!z
  if ( ismin(3) .and. in(3) ) then ! min z
     do x=(1+iin(1)),(nx-iout(1))
        do y=(1+iin(2)),(ny-iout(2))
           if (N(x,y,1)%rock_state == 0) then
              tmp_r = sum( N(x,y,1)%n_r(:)*g(:) ) - m_evp_use/amass_r
              tmp_b = sum( N(x,y,1)%n_b(:)*g(:) ) + m_evp_use/amass_b              
              if ( tmp_r .gt. 0.0d0 .and. tmp_b .gt. 0.0d0 ) then
                 N(x,y,1)%n_r(:) = tmp_r * F(:)
                 N(x,y,1)%n_b(:) = tmp_b * F(:)
              else
                 N(x,y,1)%n_r(:) = sum( N(x,y,1)%n_r(:)*g(:) ) * F(:)
                 N(x,y,1)%n_b(:) = sum( N(x,y,1)%n_b(:)*g(:) ) * F(:)
              end if
           endif
        end do
     end do
  endif 

  if ( ismax(3) .and. out(3) ) then ! max z
     do x=(1+iin(1)),(nx-iout(1))
        do y=(1+iin(2)),(ny-iout(2))
           if (N(x,y,nz)%rock_state == 0) then
!              tmp_r = sum( N(x,y,nz)%n_r(:)*g(:) ) - m_evp_use/amass_r
!              tmp_b = sum( N(x,y,nz)%n_b(:)*g(:) ) + m_evp_use/amass_b     
              tmp_r = sum( N(x,y,nz)%n_r(:)*g(:) ) - m_evp_use/amass_r
              tmp_b = sum( N(x,y,nz)%n_b(:)*g(:) ) + tmp_r             
              if ( tmp_r .gt. 0.0d0 .and. tmp_b .gt. 0.0d0 ) then
                 N(x,y,nz)%n_r(:) = tmp_r * F(:)
                 N(x,y,nz)%n_b(:) = tmp_b * F(:)
              else
                 N(x,y,nz)%n_r(:) = sum( N(x,y,nz)%n_r(:)*g(:) ) * F(:)
                 N(x,y,nz)%n_b(:) = sum( N(x,y,nz)%n_b(:)*g(:) ) * F(:)
              end if
           endif
        end do
     end do
  endif 




!And now as we are done with the red blue conversion, the same again with green to red and then green to blue
!m_evp_gr:
!x
#ifndef NOSURFACTANT
  if ( ismin(1) .and. in(1) ) then ! min x
     do y=1,ny
        do z=1,nz
           if (N(1,y,z)%rock_state == 0) then
              tmp_s = sum( N(1,y,z)%n_s(:)*g(:) ) - m_evp_gr/amass_s
              tmp_r = sum( N(1,y,z)%n_r(:)*g(:) ) + m_evp_gr/amass_r              
              if ( tmp_r .gt. 0.0d0 .and. tmp_b .gt. 0.0d0 ) then
                 N(1,y,z)%n_s(:) = tmp_s * F(:)
                 N(1,y,z)%n_r(:) = tmp_r * F(:)
              else
                 N(1,y,z)%n_s(:) = sum( N(1,y,z)%n_s(:)*g(:) ) * F(:)
                 N(1,y,z)%n_r(:) = sum( N(1,y,z)%n_r(:)*g(:) ) * F(:)
              end if
           endif
        end do
     end do
  endif 

  if ( ismax(1) .and. out(1) ) then ! max x
     do y=1,ny
        do z=1,nz
           if (N(nx,y,z)%rock_state == 0) then
              tmp_s = sum( N(nx,y,z)%n_s(:)*g(:) ) - m_evp_gr/amass_s
              tmp_r = sum( N(nx,y,z)%n_r(:)*g(:) ) + m_evp_gr/amass_r                 
              if ( tmp_r .gt. 0.0d0 .and. tmp_b .gt. 0.0d0 ) then
                 N(nx,y,z)%n_s(:) = tmp_s * F(:)
                 N(nx,y,z)%n_r(:) = tmp_r * F(:)
              else
                 N(nx,y,z)%n_s(:) = sum( N(nx,y,z)%n_s(:)*g(:) ) * F(:)
                 N(nx,y,z)%n_r(:) = sum( N(nx,y,z)%n_r(:)*g(:) ) * F(:)
              end if
           endif
        end do
     end do
  endif 

!y
  if ( ismin(2) .and. in(2) ) then ! min y
     do x=(1+iin(1)),(nx-iout(1))
        do z=1,nz
           if (N(x,1,z)%rock_state == 0) then
              tmp_s = sum( N(x,1,z)%n_s(:)*g(:) ) - m_evp_gr/amass_s
              tmp_r = sum( N(x,1,z)%n_r(:)*g(:) ) + m_evp_gr/amass_r              
              if ( tmp_r .gt. 0.0d0 .and. tmp_b .gt. 0.0d0 ) then
                 N(x,1,z)%n_s(:) = tmp_s * F(:)
                 N(x,1,z)%n_r(:) = tmp_r * F(:)
              else
                 N(x,1,z)%n_s(:) = sum( N(x,1,z)%n_s(:)*g(:) ) * F(:)
                 N(x,1,z)%n_r(:) = sum( N(x,1,z)%n_r(:)*g(:) ) * F(:)
              end if
           endif
        end do
     end do
  endif 

  if ( ismax(2) .and. out(2) ) then ! max y
     do x=(1+iin(1)),(nx-iout(1))
        do z=1,nz
           if (N(x,ny,z)%rock_state == 0) then
              tmp_s = sum( N(x,ny,z)%n_s(:)*g(:) ) - m_evp_gr/amass_s
              tmp_r = sum( N(x,ny,z)%n_r(:)*g(:) ) + m_evp_gr/amass_r                 
              if ( tmp_r .gt. 0.0d0 .and. tmp_b .gt. 0.0d0 ) then
                 N(x,ny,z)%n_s(:) = tmp_s * F(:)
                 N(x,ny,z)%n_r(:) = tmp_r * F(:)
              else
                 N(x,ny,z)%n_s(:) = sum( N(x,ny,z)%n_s(:)*g(:) ) * F(:)
                 N(x,ny,z)%n_r(:) = sum( N(x,ny,z)%n_r(:)*g(:) ) * F(:)
              end if
           endif
        end do
     end do
  endif 

!z
  if ( ismin(3) .and. in(3) ) then ! min z
     do x=(1+iin(1)),(nx-iout(1))
        do y=(1+iin(2)),(ny-iout(2))
           if (N(x,y,1)%rock_state == 0) then
              tmp_s = sum( N(x,y,1)%n_s(:)*g(:) ) - m_evp_gr/amass_s
              tmp_r = sum( N(x,y,1)%n_r(:)*g(:) ) + m_evp_gr/amass_r              
              if ( tmp_r .gt. 0.0d0 .and. tmp_b .gt. 0.0d0 ) then
                 N(x,y,1)%n_s(:) = tmp_s * F(:)
                 N(x,y,1)%n_r(:) = tmp_r * F(:)
              else
                 N(x,y,1)%n_s(:) = sum( N(x,y,1)%n_s(:)*g(:) ) * F(:)
                 N(x,y,1)%n_r(:) = sum( N(x,y,1)%n_r(:)*g(:) ) * F(:)
              end if
           endif
        end do
     end do
  endif 

  if ( ismax(3) .and. out(3) ) then ! max z
     do x=(1+iin(1)),(nx-iout(1))
        do y=(1+iin(2)),(ny-iout(2))
           if (N(x,y,nz)%rock_state == 0) then
              tmp_s = sum( N(x,y,nz)%n_s(:)*g(:) ) - m_evp_gr/amass_s
              tmp_r = sum( N(x,y,nz)%n_r(:)*g(:) ) + m_evp_gr/amass_r                 
              if ( tmp_r .gt. 0.0d0 .and. tmp_b .gt. 0.0d0 ) then
                 N(x,y,nz)%n_s(:) = tmp_s * F(:)
                 N(x,y,nz)%n_r(:) = tmp_r * F(:)
              else
                 N(x,y,nz)%n_s(:) = sum( N(x,y,nz)%n_s(:)*g(:) ) * F(:)
                 N(x,y,nz)%n_r(:) = sum( N(x,y,nz)%n_r(:)*g(:) ) * F(:)
              end if
           endif
        end do
     end do
  endif 

!And now the evaporation green to blue
!m_evp_gb:
!x
  if ( ismin(1) .and. in(1) ) then ! min x
     do y=1,ny
        do z=1,nz
           if (N(1,y,z)%rock_state == 0) then
              tmp_s = sum( N(1,y,z)%n_s(:)*g(:) ) - m_evp_gb/amass_s
              tmp_b = sum( N(1,y,z)%n_b(:)*g(:) ) + m_evp_gb/amass_b              
              if ( tmp_s .gt. 0.0d0 .and. tmp_b .gt. 0.0d0 ) then
                 N(1,y,z)%n_s(:) = tmp_s * F(:)
                 N(1,y,z)%n_b(:) = tmp_b * F(:)
              else
                 N(1,y,z)%n_s(:) = sum( N(1,y,z)%n_s(:)*g(:) ) * F(:)
                 N(1,y,z)%n_b(:) = sum( N(1,y,z)%n_b(:)*g(:) ) * F(:)
              end if
           endif
        end do
     end do
  endif 

  if ( ismax(1) .and. out(1) ) then ! max x
     do y=1,ny
        do z=1,nz
           if (N(nx,y,z)%rock_state == 0) then
              tmp_s = sum( N(nx,y,z)%n_s(:)*g(:) ) - m_evp_gb/amass_s
              tmp_b = sum( N(nx,y,z)%n_b(:)*g(:) ) + m_evp_gb/amass_b                 
              if ( tmp_r .gt. 0.0d0 .and. tmp_b .gt. 0.0d0 ) then
                 N(nx,y,z)%n_s(:) = tmp_s * F(:)
                 N(nx,y,z)%n_b(:) = tmp_b * F(:)
              else
                 N(nx,y,z)%n_s(:) = sum( N(nx,y,z)%n_s(:)*g(:) ) * F(:)
                 N(nx,y,z)%n_b(:) = sum( N(nx,y,z)%n_b(:)*g(:) ) * F(:)
              end if
           endif
        end do
     end do
  endif 

!y
  if ( ismin(2) .and. in(2) ) then ! min y
     do x=(1+iin(1)),(nx-iout(1))
        do z=1,nz
           if (N(x,1,z)%rock_state == 0) then
              tmp_s = sum( N(x,1,z)%n_s(:)*g(:) ) - m_evp_gb/amass_s
              tmp_b = sum( N(x,1,z)%n_b(:)*g(:) ) + m_evp_gb/amass_b              
              if ( tmp_r .gt. 0.0d0 .and. tmp_b .gt. 0.0d0 ) then
                 N(x,1,z)%n_s(:) = tmp_s * F(:)
                 N(x,1,z)%n_b(:) = tmp_b * F(:)
              else
                 N(x,1,z)%n_s(:) = sum( N(x,1,z)%n_s(:)*g(:) ) * F(:)
                 N(x,1,z)%n_b(:) = sum( N(x,1,z)%n_b(:)*g(:) ) * F(:)
              end if
           endif
        end do
     end do
  endif 

  if ( ismax(2) .and. out(2) ) then ! max y
     do x=(1+iin(1)),(nx-iout(1))
        do z=1,nz
           if (N(x,ny,z)%rock_state == 0) then
              tmp_s = sum( N(x,ny,z)%n_s(:)*g(:) ) - m_evp_gb/amass_s
              tmp_b = sum( N(x,ny,z)%n_b(:)*g(:) ) + m_evp_gb/amass_b                 
              if ( tmp_r .gt. 0.0d0 .and. tmp_b .gt. 0.0d0 ) then
                 N(x,ny,z)%n_s(:) = tmp_s * F(:)
                 N(x,ny,z)%n_b(:) = tmp_b * F(:)
              else
                 N(x,ny,z)%n_s(:) = sum( N(x,ny,z)%n_s(:)*g(:) ) * F(:)
                 N(x,ny,z)%n_b(:) = sum( N(x,ny,z)%n_b(:)*g(:) ) * F(:)
              end if
           endif
        end do
     end do
  endif 

!z
  if ( ismin(3) .and. in(3) ) then ! min z
     do x=(1+iin(1)),(nx-iout(1))
        do y=(1+iin(2)),(ny-iout(2))
           if (N(x,y,1)%rock_state == 0) then
              tmp_s = sum( N(x,y,1)%n_s(:)*g(:) ) - m_evp_gb/amass_s
              tmp_b = sum( N(x,y,1)%n_b(:)*g(:) ) + m_evp_gb/amass_b              
              if ( tmp_r .gt. 0.0d0 .and. tmp_b .gt. 0.0d0 ) then
                 N(x,y,1)%n_s(:) = tmp_s * F(:)
                 N(x,y,1)%n_b(:) = tmp_b * F(:)
              else
                 N(x,y,1)%n_s(:) = sum( N(x,y,1)%n_s(:)*g(:) ) * F(:)
                 N(x,y,1)%n_b(:) = sum( N(x,y,1)%n_b(:)*g(:) ) * F(:)
              end if
           endif
        end do
     end do
  endif 

  if ( ismax(3) .and. out(3) ) then ! max z
     do x=(1+iin(1)),(nx-iout(1))
        do y=(1+iin(2)),(ny-iout(2))
           if (N(x,y,nz)%rock_state == 0) then
              tmp_s = sum( N(x,y,nz)%n_s(:)*g(:) ) - m_evp_gb/amass_s
              tmp_b = sum( N(x,y,nz)%n_b(:)*g(:) ) + m_evp_gb/amass_b                 
              if ( tmp_r .gt. 0.0d0 .and. tmp_b .gt. 0.0d0 ) then
                 N(x,y,nz)%n_s(:) = tmp_s * F(:)
                 N(x,y,nz)%n_b(:) = tmp_b * F(:)
              else
                 N(x,y,nz)%n_s(:) = sum( N(x,y,nz)%n_s(:)*g(:) ) * F(:)
                 N(x,y,nz)%n_b(:) = sum( N(x,y,nz)%n_b(:)*g(:) ) * F(:)
              end if
           endif
        end do
     end do
  endif 
#endif
!^endif of #ifndef NOSURFACTANT

endif
!above is the endif of if (m_evp_set_density) then

#endif
!^endif of #ifndef SINGLEFLUID

end subroutine lbe_invade_evaporation







!> here we have a general way to set the velocity on the z-boundaries:
!> no restriction to the direction of the flow. Additionally, an oscillating
!> infolw direction on a restricted area in the center of the bottom-xy-plane
!> can be specified. This should be a demonstrtion case to generate vortices.
!> work is ongoing on this boundary condition. Contact Martin for the latest
!> information
!> currently, there are the following variantsS:
!> itype = 1 : the outflux boundary is exactly the same as the influx, including
!>             width and frequency of the jet, which makes only sense
!>             if omega=0 or width = 0
!> itype = 2 : the outflux is parallel to the z-direction and averaged over the 
!>             influx area
!> itype = 3 : the outflux is parallel to z and restricted to the central area, 
!>             like on the influx plande
!> negative itype values denote that u is the mass flow and not the velocity
subroutine lbe_invade_zouandhe_angle(N,ux_i,uy_i,uz_i,width,omega,itype)
	implicit none
	type(lbe_site),dimension(0:,0:,0:), intent(inout)  :: N   ! lattice
        integer, intent(in) :: itype
        real*8, intent(in) :: ux_i,uy_i,uz_i,width,omega
	real*8 :: rho,ux,uy,uz,uintot
        real*8 :: F(19)
	integer :: x,y,z
	integer :: i,j,k
        
        ! itype = 1 : outflux = influx
        ! itype = 2 : outflux || z
        ! itype = 3 : outflux || z und schmal
        ! itype = 4 : parabolisch
        ! itype = 5 : parabolisch, aber Nx = Ny = 0
        ! itype = 6 : Scherfluss
        ! itype = 7 : Scherfluss new, width and omega are ignored
        ! itype < 0 : u sets momentum flux instead of velocity


        
        ! calculate total influx
        uintot = 0.
        z=1
        
        if ((width > 0).and.(itype.ne.7)) then
            do y=1,tny
                do x=1,tnx
                   uz=uz_i
                   if (width>0.) then
                        uz=uz*(((0.5-0.5*cos(y*2.*pi/tny)) &
                              *(0.5-0.5*cos(x*2.*pi/tnx)))**width)
                   endif
                   uintot = uintot + uz
                enddo
            enddo
        else
          uintot  = uz_i*(tnx*tny)
        endif          

        
        if(start(3)==1) then
           do y=1,ny
             do x=1,nx
               if (N(x,y,z)%rock_state == 0) then
                  ux=ux_i
                  uy=uy_i
                  uz=uz_i
                  if (width>0.) then

                    uz=uz*(((0.5-0.5*cos((y+start(2)-1)*2.*pi/tny)) &
                          *(0.5-0.5*cos((x+start(1)-1)*2.*pi/tnx)))**width)
                    ux=ux*(((0.5-0.5*cos((y+start(2)-1)*2.*pi/tny)) &
                          *(0.5-0.5*cos((x+start(1)-1)*2.*pi/tnx)))**width)
                    uy=uy*(((0.5-0.5*cos((y+start(2)-1)*2.*pi/tny)) &
                          *(0.5-0.5*cos((x+start(1)-1)*2.*pi/tnx)))**width)
                  
                    if ((itype == 4).or.(itype == 5)) then
                       uz=uz_i*&
                         (1.0-(((x+start(1)-1)-(tnx/2.0+0.5-omega))/width)*&
                              (((x+start(1)-1)-(tnx/2.0+0.5-omega))/width))
                       if (uz.lt.0.0) uz=0.
                       ux=ux_i*&
                         (1.0-(((x+start(1)-1)-(tnx/2.0+0.5-omega))/width)*&
                              (((x+start(1)-1)-(tnx/2.0+0.5-omega))/width))
                       if (uz.le.0.0) ux=0.
                       uy=uy_i*&
                         (1.0-(((x+start(1)-1)-(tnx/2.0+0.5-omega))/width)*&
                              (((x+start(1)-1)-(tnx/2.0+0.5-omega))/width))
                       if (uz.le.0.0) uy=0.
                    endif
                  endif

                  if ((omega>0.).and.(itype.ne.7)) then
                    if ((itype .ne. 4).and.(itype .ne. 5)) then
                       ux=ux*cos(timesteps_count*omega)
                       uy=uy*sin(timesteps_count*omega)
                    endif
                  endif


#ifndef NOSURFACTANT
                        F = N(x,y,1)%n_s(:)*g(:) 

                        if(itype<0) then
                          rho =  F(19) + F(1) + F(2) + F(3) + F(4) +  &
                                  F(7) + F(11) + F(12) + F(8) +  &
                                  2.*( F(6) + F(10) + F(14) + F(16) + F(18)) + uz
                          if (rho>0.) then
                            ux = ux / rho
                            uy = uy / rho
                            uz = uz / rho
                          else
                            uz = 0.
                            ux = 0.
                            uy = 0.
                          endif
                        else
                          rho = (F(19) + F(1) + F(2) + F(3) + F(4) +  &
                                  F(7) + F(11) + F(12) + F(8) +  &
                                  2.*( F(6) + F(10) + F(14) + F(16) + F(18)))/(1.-uz)
                        end if

                        F(5)  = F(6) + uz*rho/3.
                        F(15) = F(18) - (F(7)+F(11)+F(3)-F(8)-F(12)-F(4))/2. &
                                + uz*rho/6. + uy*rho/2.
                        F(17) = F(16) - (F(8)+F(12)+F(4)-F(7)-F(11)-F(3))/2. &
                                + uz*rho/6. - uy*rho/2.
                        F(9)  = F(14) - (F(7)+F(8)+F(1)-F(11)-F(12)-F(2))/2. &
                                + uz*rho/6. + ux*rho/2.
                        F(13) = F(10) - (F(11)+F(12)+F(2)-F(7)-F(8)-F(1))/2. &
                                + uz*rho/6. - ux*rho/2.
                        
                       if (itype==5) then
                           F(15) = F(18) - (F(3)-F(4))/4. + rho*uz/6.+ uy*rho/6.
                           F(17) = F(16) + (F(3)-F(4))/4. + rho*uz/6.- uy*rho/6.
                           F(9)  = F(14) - (F(1)-F(2))/4. + rho*uz/6.+ ux*rho/6.
                           F(13) = F(10) + (F(1)-F(2))/4. + rho*uz/6.- ux*rho/6.
                       endif

                        N(x,y,1)%n_s(:) = F(:)/g(:) 
#endif
                        F = N(x,y,1)%n_r(:)*g(:) 

                        if(itype<0) then
                          rho =  F(19) + F(1) + F(2) + F(3) + F(4) +  &
                                  F(7) + F(11) + F(12) + F(8) +  &
                                  2.*( F(6) + F(10) + F(14) + F(16) + F(18)) + uz
                          if (rho>0.) then
                            ux = ux / rho
                            uy = uy / rho
                            uz = uz / rho
                          else
                            uz = 0.
                            ux = 0.
                            uy = 0.
                          endif
                        else
                          rho = (F(19) + F(1) + F(2) + F(3) + F(4) +  &
                                  F(7) + F(11) + F(12) + F(8) +  &
                                  2.*( F(6) + F(10) + F(14) + F(16) + F(18)))/(1.-uz)
                        end if

                        F(5)  = F(6) + uz*rho/3.
                        F(15) = F(18) - (F(7)+F(11)+F(3)-F(8)-F(12)-F(4))/2. &
                                + uz*rho/6. + uy*rho/2.
                        F(17) = F(16) - (F(8)+F(12)+F(4)-F(7)-F(11)-F(3))/2. &
                                + uz*rho/6. - uy*rho/2.
                        F(9)  = F(14) - (F(7)+F(8)+F(1)-F(11)-F(12)-F(2))/2. &
                                + uz*rho/6. + ux*rho/2.
                        F(13) = F(10) - (F(11)+F(12)+F(2)-F(7)-F(8)-F(1))/2. &
                                + uz*rho/6. - ux*rho/2.

                       if (itype==5) then
                           F(15) = F(18) - (F(3)-F(4))/4. + rho*uz/6.+ uy*rho/6.
                           F(17) = F(16) + (F(3)-F(4))/4. + rho*uz/6.- uy*rho/6.
                           F(9)  = F(14) - (F(1)-F(2))/4. + rho*uz/6.+ ux*rho/6.
                           F(13) = F(10) + (F(1)-F(2))/4. + rho*uz/6.- ux*rho/6.
                       endif

                        N(x,y,1)%n_r(:) = F(:)/g(:) 
#ifndef SINGLEFLUID
                        F = N(x,y,1)%n_b(:)*g(:) 

                        if(itype<0) then
                          rho =  F(19) + F(1) + F(2) + F(3) + F(4) +  &
                                  F(7) + F(11) + F(12) + F(8) +  &
                                  2.*( F(6) + F(10) + F(14) + F(16) + F(18)) + uz
                          if (rho>0.) then
                            ux = ux / rho
                            uy = uy / rho
                            uz = uz / rho
                          else
                            uz = 0.
                            ux = 0.
                            uy = 0.
                          endif
                        else
                          rho = (F(19) + F(1) + F(2) + F(3) + F(4) +  &
                                  F(7) + F(11) + F(12) + F(8) +  &
                                  2.*( F(6) + F(10) + F(14) + F(16) + F(18)))/(1.-uz)
                        end if

                        F(5)  = F(6) + uz*rho/3.
                        F(15) = F(18) - (F(7)+F(11)+F(3)-F(8)-F(12)-F(4))/2. &
                                + uz*rho/6. + uy*rho/2.
                        F(17) = F(16) - (F(8)+F(12)+F(4)-F(7)-F(11)-F(3))/2. &
                                + uz*rho/6. - uy*rho/2.
                        F(9)  = F(14) - (F(7)+F(8)+F(1)-F(11)-F(12)-F(2))/2. &
                                + uz*rho/6. + ux*rho/2.
                        F(13) = F(10) - (F(11)+F(12)+F(2)-F(7)-F(8)-F(1))/2. &
                                + uz*rho/6. - ux*rho/2.

                       if (itype==5) then
                           F(15) = F(18) - (F(3)-F(4))/4. + rho*uz/6.+ uy*rho/6.
                           F(17) = F(16) + (F(3)-F(4))/4. + rho*uz/6.- uy*rho/6.
                           F(9)  = F(14) - (F(1)-F(2))/4. + rho*uz/6.+ ux*rho/6.
                           F(13) = F(10) + (F(1)-F(2))/4. + rho*uz/6.- ux*rho/6.
                       endif

                        N(x,y,1)%n_b(:) = F(:)/g(:) 
#endif
			endif
	            end do
        	 end do
      	endif 

	if (start(3)>=(tnz-nz)) then
            z=nz
            do y=1,ny
               do x=1,nx
               
                  ux=ux_i
                  uy=uy_i
                  uz=uz_i

                 if ((width>0.).and.((itype==3).or.(itype==-3).or.(itype==1).or.(itype==-1))) then
                    uz=uz*(((0.5-0.5*cos((y+start(2)-1)*2.*pi/tny)) &
                          *(0.5-0.5*cos((x+start(1)-1)*2.*pi/tnx)))**width)
                    ux=ux*(((0.5-0.5*cos((y+start(2)-1)*2.*pi/tny)) &
                          *(0.5-0.5*cos((x+start(1)-1)*2.*pi/tnx)))**width)
                    uy=uy*(((0.5-0.5*cos((y+start(2)-1)*2.*pi/tny)) &
                          *(0.5-0.5*cos((x+start(1)-1)*2.*pi/tnx)))**width)
                  else
                     uz = uintot/(tnx*tny)
                  endif
                  
                  if ((itype>1).or.(itype<-1)) then
                     if((itype.eq.6).or.(itype.eq.-6).or.(itype.eq.7).or.(itype.eq.-7)) then
                       ux = -1.0 * ux;
                       uy = -1.0 * uy;
                     else
                       ux = 0.
                       uy = 0.
                     endif
                  endif
                  if ((itype == 4).or.(itype == 5)) then
                     uz=uz_i*&
                         (1.0-(((x+start(1)-1)-(tnx/2.0+0.5+omega))/width)*&
                              (((x+start(1)-1)-(tnx/2.0+0.5+omega))/width))
                     if (uz.lt.0.0) uz=0.
                     ux=ux_i*&
                         (1.0-(((x+start(1)-1)-(tnx/2.0+0.5+omega))/width)*&
                              (((x+start(1)-1)-(tnx/2.0+0.5+omega))/width))
                     if (uz.le.0.0) ux=0.
                     uy=uy_i*&
                         (1.0-(((x+start(1)-1)-(tnx/2.0+0.5+omega))/width)*&
                              (((x+start(1)-1)-(tnx/2.0+0.5+omega))/width))
                     if (uz.le.0.0) uy=0.
                  endif

                  !print *,"itype = ",itype

                  if ((omega>0.).and.(itype.ne.7)) then
                    if ((itype .ne. 4).and.(itype .ne. 5)) then
                       ux=ux*cos(timesteps_count*omega)
                       uy=uy*sin(timesteps_count*omega)
                    endif
                  endif
                  if (N(x,y,z)%rock_state == 0) then
#ifndef NOSURFACTANT
                        F = N(x,y,z)%n_s(:)*g(:) 

                        if(itype<0) then
                          rho =  F(19) + F(1) + F(2) + F(3) + F(4) +  &
                                 F(7) + F(11) + F(12) + F(8) +  &
                                 2.*( F(5) + F(9) + F(13) + F(15) + F(17)) - uz
                          if (rho>0.) then
                            ux = ux / rho
                            uy = uy / rho
                            uz = uz / rho
                          else
                            uz = 0.
                            ux = 0.
                            uy = 0.
                          endif
                        else
                          rho = (F(19) + F(1) + F(2) + F(3) + F(4) +  &
                                 F(7) + F(11) + F(12) + F(8) +  &
                                 2.*( F(5) + F(9) + F(13) + F(15) + F(17)))/(uz+1.)
                        end if

                        F(6) = F(5)-uz*rho/3.
                        F(18) = F(15) + (F(7)+F(11)+F(3)-F(8)-F(12)-F(4))/2. &
                                - uz*rho/6. - uy*rho/2.
                        F(16) = F(17) + (F(8)+F(12)+F(4)-F(7)-F(11)-F(3))/2. &
                                - uz*rho/6. + uy*rho/2.
                        F(14)  = F(9) + (F(7)+F(8)+F(1)-F(11)-F(12)-F(2))/2. &
                                - uz*rho/6. - ux*rho/2.
                        F(10) = F(13) + (F(11)+F(12)+F(2)-F(7)-F(8)-F(1))/2. &
                                - uz*rho/6. + ux*rho/2.

                        
                        if (itype==5) then
                          F(18) = F(15) + (F(3)-F(4))/4. - rho*uz/6. - uy*rho/6.
                          F(16) = F(17) - (F(3)-F(4))/4. - rho*uz/6. + uy*rho/6.
                          F(14) = F(9)  + (F(1)-F(2))/4. - rho*uz/6. - ux*rho/6.
                          F(10) = F(13) - (F(1)-F(2))/4. - rho*uz/6. + ux*rho/6.
                        endif

                        N(x,y,nz)%n_s(:) = F(:)/g(:) 
#endif
                        F = N(x,y,nz)%n_r(:)*g(:) 

                        if(itype<0) then
                          rho =  F(19) + F(1) + F(2) + F(3) + F(4) +  &
                                 F(7) + F(11) + F(12) + F(8) +  &
                                 2.*( F(5) + F(9) + F(13) + F(15) + F(17)) - uz
                          if (rho>0.) then
                            ux = ux / rho
                            uy = uy / rho
                            uz = uz / rho
                          else
                            uz = 0.
                            ux = 0.
                            uy = 0.
                          endif
                        else
                          rho = (F(19) + F(1) + F(2) + F(3) + F(4) +  &
                                 F(7) + F(11) + F(12) + F(8) +  &
                                 2.*( F(5) + F(9) + F(13) + F(15) + F(17)))/(uz+1.)
                        end if

                        F(6) = F(5)-uz*rho/3.
                        F(18) = F(15) + (F(7)+F(11)+F(3)-F(8)-F(12)-F(4))/2. &
                                - uz*rho/6. - uy*rho/2.
                        F(16) = F(17) + (F(8)+F(12)+F(4)-F(7)-F(11)-F(3))/2. &
                                - uz*rho/6. + uy*rho/2.
                        F(14)  = F(9) + (F(7)+F(8)+F(1)-F(11)-F(12)-F(2))/2. &
                                - uz*rho/6. - ux*rho/2.
                        F(10) = F(13) + (F(11)+F(12)+F(2)-F(7)-F(8)-F(1))/2. &
                                - uz*rho/6. + ux*rho/2.
                        
                        if (itype==5) then
                          F(18) = F(15) + (F(3)-F(4))/4. - rho*uz/6. - uy*rho/6.
                          F(16) = F(17) - (F(3)-F(4))/4. - rho*uz/6. + uy*rho/6.
                          F(14) = F(9)  + (F(1)-F(2))/4. - rho*uz/6. - ux*rho/6.
                          F(10) = F(13) - (F(1)-F(2))/4. - rho*uz/6. + ux*rho/6.
                        endif

                        N(x,y,nz)%n_r(:) = F(:)/g(:) 
#ifndef SINGLEFLUID
                        F = N(x,y,nz)%n_b(:)*g(:) 
                        if(itype<0) then
                          rho =  F(19) + F(1) + F(2) + F(3) + F(4) +  &
                                 F(7) + F(11) + F(12) + F(8) +  &
                                 2.*( F(5) + F(9) + F(13) + F(15) + F(17)) - uz
                          if (rho>0.) then
                            ux = ux / rho
                            uy = uy / rho
                            uz = uz / rho
                          else
                            uz = 0.
                            ux = 0.
                            uy = 0.
                          endif
                        else
                          rho = (F(19) + F(1) + F(2) + F(3) + F(4) +  &
                                 F(7) + F(11) + F(12) + F(8) +  &
                                 2.*( F(5) + F(9) + F(13) + F(15) + F(17)))/(uz+1.)
                        end if

                        F(6) = F(5)-uz*rho/3.
                        F(18) = F(15) + (F(7)+F(11)+F(3)-F(8)-F(12)-F(4))/2. &
                                - uz*rho/6. - uy*rho/2.
                        F(16) = F(17) + (F(8)+F(12)+F(4)-F(7)-F(11)-F(3))/2. &
                                - uz*rho/6. + uy*rho/2.
                        F(14)  = F(9) + (F(7)+F(8)+F(1)-F(11)-F(12)-F(2))/2. &
                                - uz*rho/6. - ux*rho/2.
                        F(10) = F(13) + (F(11)+F(12)+F(2)-F(7)-F(8)-F(1))/2. &
                                - uz*rho/6. + ux*rho/2.
                        
                        if (itype==5) then
                          F(18) = F(15) + (F(3)-F(4))/4. - rho*uz/6. - uy*rho/6.
                          F(16) = F(17) - (F(3)-F(4))/4. - rho*uz/6. + uy*rho/6.
                          F(14) = F(9)  + (F(1)-F(2))/4. - rho*uz/6. - ux*rho/6.
                          F(10) = F(13) - (F(1)-F(2))/4. - rho*uz/6. + ux*rho/6.
                        endif

                        N(x,y,nz)%n_b(:) = F(:)/g(:) 
#endif
                        endif
                 end do
                end do
        endif 

end subroutine lbe_invade_zouandhe_angle

!> Primarily, this boundary condition is thought to be a playground to
!> investigate different boundary conditions.  Attention: Some cases
!> lead to UNPHYSICAL BEHAVIOR consult Latt et al., Phys. Rev. E, 77,
!> 056703 (2008) - they give some severe arguments against the way to
!> set boundary conditions just to the equilibrium values Work is
!> ongoing on this type of boundary conditions. When it can be used
!> for production runs, this state will be announced here!  Ask Martin
!> for the latest information about this boundary condition.
!>
!> inv_type = 0 : not yet implemented: BC3 in Latt et al.
!> inv_type = 1 : Bounce back of non-equilibrium-part
!> inv_type = 2 : Set only the unknown distribution to the equilibrium value
!> inv_type = 3 : Set all of them to the equilibrium value (definitely unphysical!)
!> inv_type = 4 : like 3 but with predefined density (even more unphysical)
subroutine lbe_invade_equilibrium(N,fr,fg,fb,ux_i,uy_i,uz_i,ux_o,uy_o,itype)
	implicit none
	type(lbe_site),dimension(0:,0:,0:), intent(inout)  :: N   ! lattice
        real*8, intent(in) :: ux_i,uy_i,uz_i,ux_o,uy_o,fr,fg,fb
        integer, intent(in) :: itype
	real*8 :: rho,ux,uy,uz,uintot
        real*8 :: F(19), F_eq(19),F_new(19)
	integer :: x,y,z
	integer :: i,j,k

        if(start(3)==1) then
          ! Cray compiler complains about undefined variables if z isn't set, so set it.
          ! A value of one seems to make sense because we are looking at the boundary in
          ! the z-direction (cf. similar construct after this if-block, where we look at
          ! the top of z). Currently this function is still set to abort anyway
          ! (cf. 45def50), so when this is used again, please check it before 
          ! the error line below is removed.
          z = 1
           do y=1,ny
             do x=1,nx
               call error(&
                    &'z is not set in lbe_invade_equilibrium()---fix the code!')
               if (N(x,y,z)%rock_state == 0) then
                  ux=ux_i
                  uy=uy_i
                  uz=uz_i

#ifndef NOSURFACTANT
                        F = N(x,y,1)%n_s(:)*g(:) 
                        rho = fg
                        
                        if(itype/=4) then
                        rho = (F(19) + F(1) + F(2) + F(3) + F(4) +  &
                                F(7) + F(11) + F(12) + F(8) +  &
                                2.*( F(6) + F(10) + F(14) + F(16) + F(18)))/(1.-uz)
                        endif

                
                        ! get equilibrium distribution:
                        call boltz_dist(ux,uy,uz,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,F_eq)
                        F_eq = F_eq(:)*g(:) * rho
             
                        !bounce the non-eqillibrium - part
                        do j = 1,19
                           F_new(j)=F(j)
                           if((itype==3).or.(itype==4)) then
                             F_new(j) = F_eq(j)
                           else
                             if (cz(j)==1) then
                                i = bounce(j)
                                if (itype==2) then
                                  F_new(j) = F_eq(j)
                                else 
                                  F_new(j) = F_eq(j) + F(i) - F_eq(i)
                                endif 
                             endif
                          endif
                        enddo 
                          
                        N(x,y,1)%n_s(:) = F_new(:)/g(:) 
#endif
                        F = N(x,y,1)%n_r(:)*g(:) 
                        
                        rho = fr

                        if(itype/=4) then
                        rho = (F(19) + F(1) + F(2) + F(3) + F(4) +  &
                                F(7) + F(11) + F(12) + F(8) +  &
                                2.*( F(6) + F(10) + F(14) + F(16) + F(18)))/(1.-uz)
                        endif

                
                        ! get equilibrium distribution:
                        call boltz_dist(ux,uy,uz,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,F_eq)
                        F_eq = F_eq(:)*g(:) * rho
             
                        !bounce the non-eqillibrium - part
                        do j = 1,19
                           F_new(j)=F(j)
                           if((itype==3).or.(itype==4)) then
                             F_new(j) = F_eq(j)
                           else
                             if (cz(j)==1) then
                                i = bounce(j)
                                if (itype==2) then
                                  F_new(j) = F_eq(j)
                                else 
                                  F_new(j) = F_eq(j) + F(i) - F_eq(i)
                                endif 
                             endif
                          endif
                        enddo 

                        N(x,y,1)%n_r(:) = F_new(:)/g(:) 
#ifndef SINGLEFLUID
                        F = N(x,y,1)%n_b(:)*g(:) 

                        rho = fb

                        if(itype/=4) then
                        rho = (F(19) + F(1) + F(2) + F(3) + F(4) +  &
                                F(7) + F(11) + F(12) + F(8) +  &
                                2.*( F(6) + F(10) + F(14) + F(16) + F(18)))/(1.-uz)
                        endif

                
                        ! get equilibrium distribution:
                        call boltz_dist(ux,uy,uz,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,F_eq)
                        F_eq = F_eq(:)*g(:) * rho
             
                        !bounce the non-eqillibrium - part
                        do j = 1,19
                           F_new(j)=F(j)
                           if((itype==3).or.(itype==4)) then
                             F_new(j) = F_eq(j)
                           else
                             if (cz(j)==1) then
                                i = bounce(j)
                                if (itype==2) then
                                  F_new(j) = F_eq(j)
                                else 
                                  F_new(j) = F_eq(j) + F(i) - F_eq(i)
                                endif 
                             endif
                          endif
                        enddo 

                        N(x,y,1)%n_b(:) = F_new(:)/g(:) 
#endif
			endif
	            end do
        	 end do
      	endif 
	if (start(3)>=(tnz-nz)) then
              z=nz
		do y=1,ny
		 do x=1,nx
                   uz=uz_i
                   ux=ux_o
                   uy=uy_o

                        if (N(x,y,z)%rock_state == 0) then
#ifndef NOSURFACTANT
                        F = N(x,y,z)%n_s(:)*g(:) 
                        
                        rho = fg

                        if(itype/=4) then
                        rho = (F(19) + F(1) + F(2) + F(3) + F(4) +  &
                                F(7) + F(11) + F(12) + F(8) +  &
                                2.*( F(6) + F(10) + F(14) + F(16) + F(18)))/(1.-uz)
                        endif

                        ! get equilibrium distribution:
                        call boltz_dist(ux,uy,uz,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,F_eq)
                        F_eq = F_eq(:)*g(:) * rho
             
                        !bounce the non-eqillibrium - part
                        do j = 1,19
                           F_new(j)=F(j)
                           if((itype==3).or.(itype==4)) then
                             F_new(j) = F_eq(j)
                           else
                             if (cz(j)==-1) then
                                i = bounce(j)
                                if (itype==2) then
                                  F_new(j) = F_eq(j)
                                else 
                                  F_new(j) = F_eq(j) + F(i) - F_eq(i)
                                endif 
                             endif
                          endif
                        enddo 

                        N(x,y,nz)%n_s(:) = F_new(:)/g(:) 
#endif
                        F = N(x,y,nz)%n_r(:)*g(:) 
                        rho = fr

                        if(itype/=4) then
                        rho = (F(19) + F(1) + F(2) + F(3) + F(4) +  &
                                F(7) + F(11) + F(12) + F(8) +  &
                                2.*( F(6) + F(10) + F(14) + F(16) + F(18)))/(1.-uz)
                        endif

                        ! get equilibrium distribution:
                        call boltz_dist(ux,uy,uz,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,F_eq)
                        F_eq = F_eq(:)*g(:) * rho
             
                        !bounce the non-eqillibrium - part
                        do j = 1,19
                           F_new(j)=F(j)
                           if((itype==3).or.(itype==4)) then
                             F_new(j) = F_eq(j)
                           else
                             if (cz(j)==-1) then
                                i = bounce(j)
                                if (itype==2) then
                                  F_new(j) = F_eq(j)
                                else 
                                  F_new(j) = F_eq(j) + F(i) - F_eq(i)
                                endif 
                             endif
                          endif
                        enddo 


                        N(x,y,nz)%n_r(:) = F_new(:)/g(:) 
#ifndef SINGLEFLUID
                        F = N(x,y,nz)%n_b(:)*g(:) 
                        
                        rho = fb

                        if(itype/=4) then
                        rho = (F(19) + F(1) + F(2) + F(3) + F(4) +  &
                                F(7) + F(11) + F(12) + F(8) +  &
                                2.*( F(6) + F(10) + F(14) + F(16) + F(18)))/(1.-uz)
                        endif

                        ! get equilibrium distribution:
                        call boltz_dist(ux,uy,uz,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,0.0_8,F_eq)
                        F_eq = F_eq(:)*g(:) * rho
             
                        !bounce the non-eqillibrium - part
                        do j = 1,19
                           F_new(j)=F(j)
                           if((itype==3).or.(itype==4)) then
                             F_new(j) = F_eq(j)
                           else
                             if (cz(j)==-1) then
                                i = bounce(j)
                                if (itype==2) then
                                  F_new(j) = F_eq(j)
                                else 
                                  F_new(j) = F_eq(j) + F(i) - F_eq(i)
                                endif 
                             endif
                          endif
                        enddo 


                        N(x,y,nz)%n_b(:) = F_new(:)/g(:) 
#endif
                        endif
                 end do
                end do
        endif 

end subroutine lbe_invade_equilibrium

!> currently, at x=1 and nx the velocity is set to zero. Nothing more
subroutine lbe_invade_poiseuille(N,itype)
	implicit none
	type(lbe_site),dimension(0:,0:,0:), intent(inout)  :: N   ! lattice
        integer, intent(in) :: itype
	real*8 :: rho,ux,uy,uz,uintot
        real*8 :: F(19)
	integer :: x,y,z
	integer :: i,j,k
        
        x=1
        if(start(1)==1) then
           do y=1,ny
             do z=1,nz
               if (N(x,y,z)%rock_state == 0) then
                  ux=0.
                  uy=0.
                  uz=0.


#ifndef NOSURFACTANT
                        F = N(x,y,z)%n_s(:)*g(:) 

                        rho =  F(19) + F(5) + F(6) + F(3) + F(4) +  &
                                  F(15) + F(16) + F(17) + F(18) +  &
                                  2.*( F(2) + F(13) + F(14) + F(11) + F(12)) + uz

                        F(1)  = F(2) + ux*rho/3.
                        F(9) = F(14) - (F(17)+F(15)+F(5)-F(18)-F(16)-F(6))/2. &
                                + ux*rho/6. + uz*rho/2.
                        F(10) = F(13) + (F(17)+F(15)+F(5)-F(18)-F(16)-F(6))/2. &
                                + ux*rho/6. - uz*rho/2.
                        F(7)  = F(12) - (F(15)+F(16)+F(3)-F(17)-F(18)-F(4))/2. &
                                + ux*rho/6. + uy*rho/2.
                        F(8) =  F(11) + (F(15)+F(16)+F(3)-F(17)-F(18)-F(4))/2. &
                                + ux*rho/6. - uy*rho/2.
                        if (itype==1) then
                          F(1)  = F(2) + ux*rho/3.
                          F(9) = F(14) + ux*rho/6. + uz*rho/6.
                          F(10) = F(13) + ux*rho/6. - uz*rho/6.
                          F(7)  = F(12) + ux*rho/6. + uy*rho/6.
                          F(8) =  F(11) + ux*rho/6. - uy*rho/6.
                        endif 
                          

                        N(x,y,z)%n_s(:) = F(:)/g(:) 
#endif
                        F = N(x,y,z)%n_r(:)*g(:) 

                        rho =  F(19) + F(5) + F(6) + F(3) + F(4) +  &
                                  F(15) + F(16) + F(17) + F(18) +  &
                                  2.*( F(2) + F(13) + F(14) + F(11) + F(12)) + uz

                        F(1)  = F(2) + ux*rho/3.
                        F(9) = F(14) - (F(17)+F(15)+F(5)-F(18)-F(16)-F(6))/2. &
                                + ux*rho/6. + uz*rho/2.
                        F(10) = F(13) + (F(17)+F(15)+F(5)-F(18)-F(16)-F(6))/2. &
                                + ux*rho/6. - uz*rho/2.
                        F(7)  = F(12) - (F(15)+F(16)+F(3)-F(17)-F(18)-F(4))/2. &
                                + ux*rho/6. + uy*rho/2.
                        F(8) =  F(11) + (F(15)+F(16)+F(3)-F(17)-F(18)-F(4))/2. &
                                + ux*rho/6. - uy*rho/2.

                        if (itype==1) then
                          F(1)  = F(2) + ux*rho/3.
                          F(9) = F(14) + ux*rho/6. + uz*rho/6.
                          F(10) = F(13) + ux*rho/6. - uz*rho/6.
                          F(7)  = F(12) + ux*rho/6. + uy*rho/6.
                          F(8) =  F(11) + ux*rho/6. - uy*rho/6.
                        endif 

                        N(x,y,z)%n_r(:) = F(:)/g(:) 
#ifndef SINGLEFLUID
                        F = N(x,y,z)%n_b(:)*g(:) 

                        rho =  F(19) + F(5) + F(6) + F(3) + F(4) +  &
                                  F(15) + F(16) + F(17) + F(18) +  &
                                  2.*( F(2) + F(13) + F(14) + F(11) + F(12)) + uz

                        F(1)  = F(2) + ux*rho/3.
                        F(9) = F(14) - (F(17)+F(15)+F(5)-F(18)-F(16)-F(6))/2. &
                                + ux*rho/6. + uz*rho/2.
                        F(10) = F(13) + (F(17)+F(15)+F(5)-F(18)-F(16)-F(6))/2. &
                                + ux*rho/6. - uz*rho/2.
                        F(7)  = F(12) - (F(15)+F(16)+F(3)-F(17)-F(18)-F(4))/2. &
                                + ux*rho/6. + uy*rho/2.
                        F(8) =  F(11) + (F(15)+F(16)+F(3)-F(17)-F(18)-F(4))/2. &
                                + ux*rho/6. - uy*rho/2.

                        if (itype==1) then
                          F(1)  = F(2) + ux*rho/3.
                          F(9) = F(14) + ux*rho/6. + uz*rho/6.
                          F(10) = F(13) + ux*rho/6. - uz*rho/6.
                          F(7)  = F(12) + ux*rho/6. + uy*rho/6.
                          F(8) =  F(11) + ux*rho/6. - uy*rho/6.
                        endif 

                        N(x,y,z)%n_b(:) = F(:)/g(:) 
#endif
			endif
	            end do
        	 end do
      	endif 

	if (start(1)>=(tnx-nx)) then
            x=nx
            do y=1,ny
               do z=1,nz
               
                  ux=0.
                  uy=0.
                  uz=0.

                  if (N(x,y,z)%rock_state == 0) then
#ifndef NOSURFACTANT
                        F = N(x,y,z)%n_s(:)*g(:) 

                        rho =  F(19) + F(5) + F(6) + F(3) + F(4) +  &
                                  F(15) + F(16) + F(17) + F(18) +  &
                                  2.*( F(1) + F(7) + F(8) + F(9) + F(10)) + uz

                        F(2)  = F(1) - ux*rho/3.
                        F(14) = F(9) + (F(17)+F(15)+F(5)-F(18)-F(16)-F(6))/2. &
                                - ux*rho/6. - uz*rho/2.
                        F(13) = F(10) - (F(17)+F(15)+F(5)-F(18)-F(16)-F(6))/2. &
                                - ux*rho/6. + uz*rho/2.
                        F(12)  = F(7) + (F(15)+F(16)+F(3)-F(17)-F(18)-F(4))/2. &
                                - ux*rho/6. - uy*rho/2.
                        F(11) = F(8) - (F(15)+F(16)+F(3)-F(17)-F(18)-F(4))/2. &
                                - ux*rho/6. + uy*rho/2.

                        if (itype==1) then
                        F(2)  = F(1) - ux*rho/3.
                        F(14) = F(9) - ux*rho/6. - uz*rho/6.
                        F(13) = F(10) - ux*rho/6. + uz*rho/6.
                        F(12)  = F(7) - ux*rho/6. - uy*rho/6.
                        F(11) = F(8) - ux*rho/6. + uy*rho/6.
                        endif 

                        N(x,y,z)%n_s(:) = F(:)/g(:) 
#endif
                        F = N(x,y,z)%n_r(:)*g(:) 

                        rho =  F(19) + F(5) + F(6) + F(3) + F(4) +  &
                                  F(15) + F(16) + F(17) + F(18) +  &
                                  2.*( F(1) + F(7) + F(8) + F(9) + F(10)) + uz

                        F(2)  = F(1) - ux*rho/3.
                        F(14) = F(9) + (F(17)+F(15)+F(5)-F(18)-F(16)-F(6))/2. &
                                - ux*rho/6. - uz*rho/2.
                        F(13) = F(10) - (F(17)+F(15)+F(5)-F(18)-F(16)-F(6))/2. &
                                - ux*rho/6. + uz*rho/2.
                        F(12)  = F(7) + (F(15)+F(16)+F(3)-F(17)-F(18)-F(4))/2. &
                                - ux*rho/6. - uy*rho/2.
                        F(11) = F(8) - (F(15)+F(16)+F(3)-F(17)-F(18)-F(4))/2. &
                                - ux*rho/6. + uy*rho/2.

                        if (itype==1) then
                        F(2)  = F(1) - ux*rho/3.
                        F(14) = F(9) - ux*rho/6. - uz*rho/6.
                        F(13) = F(10) - ux*rho/6. + uz*rho/6.
                        F(12)  = F(7) - ux*rho/6. - uy*rho/6.
                        F(11) = F(8) - ux*rho/6. + uy*rho/6.
                        endif 

                        N(x,y,z)%n_r(:) = F(:)/g(:) 
#ifndef SINGLEFLUID
                        F = N(x,y,z)%n_b(:)*g(:) 
                        rho =  F(19) + F(5) + F(6) + F(3) + F(4) +  &
                                  F(15) + F(16) + F(17) + F(18) +  &
                                  2.*( F(1) + F(7) + F(8) + F(9) + F(10)) + uz

                        F(2)  = F(1) - ux*rho/3.
                        F(14) = F(9) + (F(17)+F(15)+F(5)-F(18)-F(16)-F(6))/2. &
                                - ux*rho/6. - uz*rho/2.
                        F(13) = F(10) - (F(17)+F(15)+F(5)-F(18)-F(16)-F(6))/2. &
                                - ux*rho/6. + uz*rho/2.
                        F(12)  = F(7) + (F(15)+F(16)+F(3)-F(17)-F(18)-F(4))/2. &
                                - ux*rho/6. - uy*rho/2.
                        F(11) = F(8) - (F(15)+F(16)+F(3)-F(17)-F(18)-F(4))/2. &
                                - ux*rho/6. + uy*rho/2.

                        if (itype==1) then
                        F(2)  = F(1) - ux*rho/3.
                        F(14) = F(9) - ux*rho/6. - uz*rho/6.
                        F(13) = F(10) - ux*rho/6. + uz*rho/6.
                        F(12)  = F(7) - ux*rho/6. - uy*rho/6.
                        F(11) = F(8) - ux*rho/6. + uy*rho/6.
                        endif 

                        N(x,y,z)%n_b(:) = F(:)/g(:) 
#endif
                        endif
                 end do
                end do
        endif 
end subroutine lbe_invade_poiseuille


!> partial slip boundary Khalid Ahmed and Martin Hecht
!>
!> pr gives the slip parameter between 0 and 1
!> pb is the mean width of stripes (if any)
!> if pg is non-zero, pb+pg is the periodicity of the pattern and pg the width of one 
!> stripe, pb the one of the other one
!>    i.e. pr=1 pb=6 and pb=2 creates stripes of the witdh 2 and 6 with full- and noslip
!> shear_u the shear velocity
!> inv_type=0: homogeneous walls
!> inv_type=1: stripes in y-direction
!> inv_type=2: stripes in z-direction
!> inv_type=3: diagonal stripes 
!> inv_type=4: horizontal stripes with asymmetric slip
!> inv_type=5: vertical stripes with asymmetric slip
!> inv_type=6: continuously variing horizontal stripes 
!> inv_type=7: horizontal stripes alternating between 0 and pr
!> inv_type=8: vertical stripes alternating between 0 and pr

!> inv_type=10: horizontal (z) stripes alternating between 0 and pr, 
!> with a unit cell width pb and slip stripe width pg

!> inv_type=11: vertical (y) stripes alternating between 0 and pr, 
!> with a unit cell width pb and slip stripe width pg

!> inv_type=12: horizontal (z) stripes with a cosine slip modulation
!> according to b = pb + 2 * pb cos (2*\pi*z)

!> inv_type=13: vertical (y) stripes with a cosine slip modulation
!> according to b = pb + 2 * pb cos (2*\pi*y)

!> inv_type=13: vertical (y) stripes as in inv_type=11 with a 
!> full-slip boundary at x=tnx
 


subroutine lbe_invade_with_slip(N,pr_in,pb_in,pg_in,shear_u,w_type)
	implicit none
	type(lbe_site),dimension(0:,0:,0:), intent(inout)  :: N   ! lattice
        real*8, intent(in) :: pr_in
        real*8, intent(in) :: pb_in
        real*8, intent(in) :: pg_in
	real*8, intent(in) :: shear_u 
        integer :: pb
        integer :: pg
        real*8 :: pr

	real*8 :: shear_u_low_x = 0.d0

	integer, intent(in) :: w_type
	real*8 :: rho,ux,uy,uz,uintot
        real*8 :: F(19)
	integer :: x,y,z
	integer :: i,j,k
        integer :: tx,ty,tz
        real*8  :: b0, b1
        real*8  :: period,bLoc,bParm

	pr=pr_in
	pb=pb_in
	pg=pg_in
        x=1

        if(start(1)==1) then
                
                do y=1,ny
                        ty = ccoords(2) * ny + y
                        select case (w_type)
                        case (2)
                                if ( mod(ty,pb+pg) == 0) pr = 1-pr
                                if (pg.ne.0) then
                                        if ( mod(ty-pg,pb+pg) == 0) pr=1-pr
                                endif
                                
                        case (3)
                                if ( mod(ty,pb+pg) == 0)   pr = 1-pr
                                if (pg.ne.0) then
                                        if ( mod(ty-pg,pb+pg) == 0) pr=1-pr
                                endif
                                
                        case (5)
                                if ( mod(ty,(pb*2)) == 0) then
                                        pr=pg_in
                                elseif ( mod(ty,(pb*2)) == pb) then
                                        pr=pr_in
                                endif
                                
                        case (8)
                                if ( mod(ty,pb+pg) == 0) pr = 0
                                if (pg.ne.0) then
                                        if ( mod(ty-pg,pb+pg) == 0) pr=pr_in
                                endif
                                
                        case (11)
                                if ( mod(ty,pb) == 0) then
                                        pr=pr_in
                                endif
                                if ( mod(ty,pb) >= pg) then
                                        pr=0.
                                endif
                        case (13)
                                b0 = pr_in
                                b1 = pb_in
                                period = real(pg,kind=rk)/real(tny,kind=rk)

                                bLoc = real(b0,kind=rk) + 2.0_rk*real(b1,kind=rk) * cos(2.0_rk*real(pi,kind=rk)*real(period,kind=rk)*real(ty,kind=rk))
                                bParm = (3*bLoc) / ( (3*bLoc) + 1 )
                                pr = bParm

                        case (14)
                                if ( mod(ty,pb) == 0) then
                                        pr=pr_in
                                endif
                                if ( mod(ty,pb) >= pg) then
                                        pr=0.
                                endif

                        end select
                        do z=1,nz
                                tz = ccoords(3) * nz + z
                                select case (w_type)
                                        
                                case (1)
                                        if ( mod(tz,pb+pg) == 0) pr=1-pr
                                        if (pg.ne.0) then
                                                if ( mod(tz-pg,pb+pg) == 0) pr=1-pr
                                        endif
                                        
                                case (3)
                                        if ( mod(ty+tz,pb+pg) == 0)   pr = 1-pr
                                        if (pg.ne.0) then
                                                if ( mod(ty+tz-pg,pb+pg) == 0) pr=1-pr
                                        endif
                                        
                                case (4)
                                        if ( mod(tz,(pb*2)) == 0) then
                                                pr=pg_in
                                        elseif ( mod(tz,(pb*2)) == pb) then
                                                pr=pr_in
                                        endif

                                case (6)
                                        pr = 0.5*(pr_in-pg_in) * cos(tz/pb_in*2.*pi) + 0.5*(pg_in+pr_in)
                                        
                                case (7)
                                        if ( mod(tz,pb+pg) == 0) pr=0
                                        if (pg.ne.0) then
                                                if ( mod(tz-pg,pb+pg) == 0) pr=pr_in
                                        endif
                                        
                                case (10)
                                        if ( mod(tz,pb) == 0) then
                                                pr=pr_in
                                        endif
                                        if ( mod(tz,pb) >= pg) then
                                                pr=0.
                                        endif
                                        
                                case (12)
                                        b0 = pr_in
                                        b1 = pb_in
                                        period = real(pg,kind=rk)/real(tnz,kind=rk)
                                        bLoc = real(b0,kind=rk) + 2.0_rk*real(b1,kind=rk) * cos(2.0_rk*real(pi,kind=rk)*real(period,kind=rk)*real(tz,kind=rk))                 
                                        bParm = (3*bLoc) / ( (3*bLoc) + 1 )
                                        pr = bParm
                                end select
                                
                                
                                if (N(x,y,z)%rock_state == 0) then
                                        ux=0.
                                        uy=0.
                                        uz=0.
                                        
                                        
                                        F = N(x,y,z)%n_r(:)*g(:) 
                                        
                                        rho =  F(19) + F(5) + F(6) + F(3) + F(4) +  &
                                                F(15) + F(16) + F(17) + F(18) +  &
                                                2.*( F(2) + F(13) + F(14) + F(11) + F(12))
                                        
                                        F(1)  = F(2) + ux*rho/3.
                                        F(9) = (1-pr)*(F(14) - (F(17)+F(15)+F(5)-F(18)-F(16)-F(6))/2.) &
                                                + pr*(F(13)) &
                                                + (1-pr)*(ux*rho/6. + shear_u_low_x*rho/2.)
                                        F(10) = (1-pr)*(F(13) + (F(17)+F(15)+F(5)-F(18)-F(16)-F(6))/2.) &
                                                + pr*(F(14)) &
                                                + (1-pr)*(ux*rho/6. - shear_u_low_x*rho/2.)
                                        F(7)  = (1-pr)*(F(12) - (F(15)+F(16)+F(3)-F(17)-F(18)-F(4))/2.) &
                                                + pr*(F(11)) &
                                                + (1-pr)*(ux*rho/6. + uy*rho/2.)
                                        F(8) =  (1-pr)*(F(11) + (F(15)+F(16)+F(3)-F(17)-F(18)-F(4))/2.) &
                                                + pr*(F(12)) &
                                                + (1-pr)*(ux*rho/6. - uy*rho/2.)
                                        
                                        N(x,y,z)%n_r(:) = F(:)/g(:) 


#ifndef SINGLEFLUID
                                        F = N(x,y,z)%n_b(:)*g(:) 
                                        
                                        rho =  F(19) + F(5) + F(6) + F(3) + F(4) +  &
                                                F(15) + F(16) + F(17) + F(18) +  &
                                                2.*( F(2) + F(13) + F(14) + F(11) + F(12))
                                        
                                        F(1)  = F(2) + ux*rho/3.
                                        F(9) = (1-pr)*(F(14) - (F(17)+F(15)+F(5)-F(18)-F(16)-F(6))/2.) &
                                                + pr*(F(13)) &
                                                + (1-pr)*(ux*rho/6. + shear_u_low_x*rho/2.)
                                        F(10) = (1-pr)*(F(13) + (F(17)+F(15)+F(5)-F(18)-F(16)-F(6))/2.) &
                                                + pr*(F(14)) &
                                                + (1-pr)*(ux*rho/6. - shear_u_low_x*rho/2.)
                                        F(7)  = (1-pr)*(F(12) - (F(15)+F(16)+F(3)-F(17)-F(18)-F(4))/2.) &
                                                + pr*(F(11)) &
                                                + (1-pr)*(ux*rho/6. + uy*rho/2.)
                                        F(8) =  (1-pr)*(F(11) + (F(15)+F(16)+F(3)-F(17)-F(18)-F(4))/2.) &
                                                + pr*(F(12)) &
                                                + (1-pr)*(ux*rho/6. - uy*rho/2.)
                                        
                                        N(x,y,z)%n_b(:) = F(:)/g(:) 

#endif
                                endif
			end do
                end do
	endif


        if (start(1)>=(tnx-nx)) then
                x=nx

                do y=1,ny
                        ty = ccoords(2) * ny + y
                        select case (w_type)
                        case (0)

                        case (2)
                                if ( mod(ty,pb+pg) == 0) pr = 1-pr
                                if (pg.ne.0) then
                                        if ( mod(ty-pg,pb+pg) == 0) pr=1-pr
                                endif
                        case (3)
                                if ( mod(ty,pb+pg) == 0 )    pr = 1-pr
                                if (pg.ne.0) then
                                        if ( mod(ty-pg,pb+pg) == 0) pr=1-pr
                                endif
                        case (5)
                                if ( mod(ty,(pb*2)) == 0) then
                                        pr=pg_in
                                elseif ( mod(ty,(pb*2)) == pb) then
                                        pr=pr_in
                                endif
                        case (8)
                                if ( mod(ty,pb+pg) == 0) pr = 0
                                if (pg.ne.0) then
                                        if ( mod(ty-pg,pb+pg) == 0) pr=pr_in
                                endif
                        case (11)
                                pr = 0.
                        case (13)
                                pr = 0.
                        case (14)
                                pr = 1.
                        end select
                        do z=1,nz
                                tz = ccoords(3) * nz + z
                                select case (w_type)
                                case (0)

                                case (1)
                                        if ( mod(tz,pb+pg) == 0) pr=1-pr
                                        if (pg.ne.0) then
                                                if ( mod(tz-pg,pb+pg) == 0) pr=1-pr
                                        endif
                                case (3)
                                        if ( mod(ty+tz,pb+pg)==0 )   pr = 1-pr
                                        if (pg.ne.0) then
                                                if ( mod(ty+tz-pg,pb+pg) == 0) pr=1-pr
                                        endif
                                case (4)
                                        if ( mod(tz,(pb*2)) == 0) then
                                                pr=pg_in
                                        elseif ( mod(tz,(pb*2)) == pb) then
                                                pr=pr_in
                                        endif
                                case (6)
                                        pr = 0.5*(pr_in-pg_in) * cos(tz/pb_in*2.*pi) + 0.5*(pg_in+pr_in)

                                case (7)
                                        if ( mod(tz,pb+pg) == 0) pr=0
                                        if (pg.ne.0) then
                                                if ( mod(tz-pg,pb+pg) == 0) pr=pr_in
                                        endif

                                case (10)
                                        pr = 0.
                                case (12)
                                        pr = 0.

                                end select

                                ux=0.
                                uy=0.
                                uz=0.

                                if (N(x,y,z)%rock_state == 0) then
                                        F = N(x,y,z)%n_r(:)*g(:) 
                                        
                                        rho =  F(19) + F(5) + F(6) + F(3) + F(4) +  &
                                                F(15) + F(16) + F(17) + F(18) +  &
                                                2.*( F(1) + F(7) + F(8) + F(9) + F(10))
                                        
                                        F(2)  = F(1) - ux*rho/3.
                                        F(14) = (1-pr)*(F(9) + (F(17)+F(15)+F(5)-F(18)-F(16)-F(6))/2.) &
                                                + pr*(F(10)) &
                                                + (1-pr)*(- ux*rho/6. + shear_u*rho/2.)
                                        F(13) = (1-pr)*(F(10) - (F(17)+F(15)+F(5)-F(18)-F(16)-F(6))/2.) &
                                                + pr*(F(9)) &
                                                + (1-pr)*(- ux*rho/6. - shear_u*rho/2.)
                                        F(12)  = (1-pr)*(F(7) + (F(15)+F(16)+F(3)-F(17)-F(18)-F(4))/2.) &
                                                + pr*(F(8)) &
                                                + (1-pr)*(- ux*rho/6. - uy*rho/2.)
                                        F(11) = (1-pr)*(F(8) - (F(15)+F(16)+F(3)-F(17)-F(18)-F(4))/2.) &
                                                + pr*(F(7)) &
                                                + (1-pr)*(- ux*rho/6. + uy*rho/2.)
                                        
                                        N(x,y,z)%n_r(:) = F(:)/g(:) 
#ifndef SINGLEFLUID
                                        F = N(x,y,z)%n_b(:)*g(:) 
                                        
                                        rho =  F(19) + F(5) + F(6) + F(3) + F(4) +  &
                                                F(15) + F(16) + F(17) + F(18) +  &
                                                2.*( F(1) + F(7) + F(8) + F(9) + F(10))
                                        
                                        F(2)  = F(1) - ux*rho/3.
                                        F(14) = (1-pr)*(F(9) + (F(17)+F(15)+F(5)-F(18)-F(16)-F(6))/2.) &
                                                + pr*(F(10)) &
                                                + (1-pr)*(- ux*rho/6. + shear_u*rho/2.)
                                        F(13) = (1-pr)*(F(10) - (F(17)+F(15)+F(5)-F(18)-F(16)-F(6))/2.) &
                                                + pr*(F(9)) &
                                                + (1-pr)*(- ux*rho/6. - shear_u*rho/2.)
                                        F(12)  = (1-pr)*(F(7) + (F(15)+F(16)+F(3)-F(17)-F(18)-F(4))/2.) &
                                                + pr*(F(8)) &
                                                + (1-pr)*(- ux*rho/6. - uy*rho/2.)
                                        F(11) = (1-pr)*(F(8) - (F(15)+F(16)+F(3)-F(17)-F(18)-F(4))/2.) &
                                                + pr*(F(7)) &
                                                + (1-pr)*(- ux*rho/6. + uy*rho/2.)
                                        
                                        N(x,y,z)%n_b(:) = F(:)/g(:) 
#endif
                                endif
                        end do
                end do
        end if

end subroutine lbe_invade_with_slip



subroutine lbe_invade_channel(N,fr,pr,pb,pg,itype)
	implicit none
	type(lbe_site),dimension(0:,0:,0:), intent(inout)  :: N   ! lattice
        real*8, intent(in) :: fr,pr,pb,pg
        integer, intent(in) :: itype
	real*8 :: rho,ux,uy,uz,uintot
        real*8 :: F(19),F1(19),F2(19)
	integer :: x,y,z
	integer :: i,j,k
        integer :: startx,starty,startz,endx,endy,endz
        logical :: minx,miny,minz,maxx,maxy,maxz

	minx=(start(1)==1)
	miny=(start(2)==1)
	minz=(start(3)==1)
	maxx=(start(1)>=(tnx-nx))
	maxy=(start(2)>=(tny-ny))
	maxz=(start(3)>=(tnz-nz))


        startx = 1
        starty = 1
        startz = 1
        endx = nx
        endy = ny
        endz = nz

        if(minx) startx = 2
        if(miny) starty = 2
        if(minz) startz = 2
        if(maxx) endx = nx-1
        if(maxy) endy = ny-1
        if(maxz) endz = nz-1

        
        if(itype.eq.3) then
           call lbe_invade_slip_channel(N,pr,pb,pg)
           return
        endif
        
        ! no slip boundaries

        if(minx) then
        x=1
!           do y=starty,endy
!             do z=startz,endz
           do y=1,ny
             do z=1,nz
                  F = N(x,y,z)%n_r(:)*g(:) 


                  F(1)  = F(2) 
                  F(9)  = F(13)
                  F(10) = F(14)
                  F(7)  = F(11)
                  F(8)  = F(12)
                  
                  ! F=F*(1.0/SUM(F))
                  ! F(19) = F(19)-(SUM(F)-1.)

                 N(x,y,z)%n_r(:) = F(:)/g(:) 
	    end do
         end do
         endif

         if(maxx) then
         x=nx
!         do y=starty,endy
!            do z=startz,endz
           do y=1,ny
             do z=1,nz

                F = N(x,y,z)%n_r(:)*g(:) 

               F(2)  = F(1) 
               F(14) = F(10)
               F(13) = F(9)
               F(12) = F(8)
               F(11) = F(7)
                  
               ! F=F*(1.0/SUM(F))
               ! F(19) = F(19)-(SUM(F)-1.)

               N(x,y,z)%n_r(:) = F(:)/g(:) 
            end do
         end do
         endif


         if(miny) then
         y=1
!           do x=startx,endx
!             do z=startz,endz
           do x=1,nx
            do z=1,nz
                  F = N(x,y,z)%n_r(:)*g(:) 

                  F(3)  = F(4) 
                  F(7)  = F(8)
                  F(11) = F(12)
                  F(15) = F(17)
                  F(16) = F(18)
                  
                  ! F=F*(1.0/SUM(F))
                  ! F(19) = F(19)-(SUM(F)-1.)


                 N(x,y,z)%n_r(:) = F(:)/g(:) 
	    end do
         end do
         endif


         if(maxy) then
         y=ny
!         do x=startx,endx
!            do z=startz,endz
         do x=1,nx
            do z=1,nz
               F = N(x,y,z)%n_r(:)*g(:) 

                  F(4)  = F(3) 
                  F(12) = F(11)
                  F(8)  = F(7)
                  F(18) = F(16)
                  F(17) = F(15)
                  
                  ! F=F*(1.0/SUM(F))
                  ! F(19) = F(19)-(SUM(F)-1.)

               N(x,y,z)%n_r(:) = F(:)/g(:) 
            end do
         end do
         end if



        ! influx

        ux = 0.
        uy = 0.
        uz = pr
        
        if(minz) then
        z=1

	do y=starty,endy
	 do x=startx,endx
!	do y=1,ny
!	 do x=1,nx
                F = N(x,y,z)%n_r(:)*g(:) 
                
                        
!              if (itype==0) then

               uz=pr


               if (uz.lt.0.0) uz=0.

               rho = (F(19) + F(1) + F(2) + F(3) + F(4) +  &
                          F(7) + F(11) + F(12) + F(8) +  &
                          2.*( F(6) + F(10) + F(14) + F(16) + F(18)))/(1.-uz)

                F(5)  = F(6) + uz*rho/3.
                F(15) = F(18) - (F(7)+F(11)+F(3)-F(8)-F(12)-F(4))/2. &
                        + uz*rho/6. + uy*rho/2.
                F(17) = F(16) - (F(8)+F(12)+F(4)-F(7)-F(11)-F(3))/2. &
                        + uz*rho/6. - uy*rho/2.
                F(9)  = F(14) - (F(7)+F(8)+F(1)-F(11)-F(12)-F(2))/2. &
                        + uz*rho/6. + ux*rho/2.
                F(13) = F(10) - (F(11)+F(12)+F(2)-F(7)-F(8)-F(1))/2. &
                        + uz*rho/6. - ux*rho/2.

!              else if (itype==1) then
!
!                rho = fr
!
!                if(rho>0.) then
!                uz = 1. - (F(19) + &
!                               F(1) + F(2) + F(3) + F(4) +  &
!                               F(7) + F(11) + F(12) + F(8) +  &
!                               2.*( F(6) + F(10) + F(14) + F(16) + F(18)))/rho
!                else 
!                  uz = 0.
!                endif                                       
!                F(5)  = F(6) + rho*uz/3.
!
!                  F(15) = F(18) - (F(7)+F(11)+F(3)-F(8)-F(12)-F(4))/2. &
!                          + uz*rho/6. + uy*rho/2.
!                  F(17) = F(16) - (F(8)+F(12)+F(4)-F(7)-F(11)-F(3))/2. &
!                          + uz*rho/6. - uy*rho/2.
!                  F(9)  = F(14) - (F(7)+F(8)+F(1)-F(11)-F(12)-F(2))/2. &
!                          + uz*rho/6. + ux*rho/2.
!                  F(13) = F(10) - (F(11)+F(12)+F(2)-F(7)-F(8)-F(1))/2. &
!                          + uz*rho/6. - ux*rho/2.
!                end if


                N(x,y,z)%n_r(:) = F(:)/g(:) 
	 end do
	end do

!	do y=starty,endy
!	 do x=startx,endx
	do y=1,ny
	 do x=1,nx
          N(x,y,z+1)%n_r(:) = N(x,y,z)%n_r(:)
       !    N(x,y,z+1)%n_r(5) = N(x,y,z)%n_r(5)
       !    N(x,y,z+1)%n_r(15) = N(x,y,z)%n_r(15)
       !    N(x,y,z+1)%n_r(17) = N(x,y,z)%n_r(17)
       !    N(x,y,z+1)%n_r(9) = N(x,y,z)%n_r(9)
       !    N(x,y,z+1)%n_r(13) = N(x,y,z)%n_r(13)
         end do
        end do

        endif

        ! outflux
        if(maxz) then
        z=nz
        
	do y=starty,endy
	 do x=startx,endx
!	do y=1,ny
!	 do x=1,nx
!         N(x,y,z-1)%n_r(:) =2.*N(x,y,z-2)%n_r(:)-N(x,y,z-3)%n_r(:)
!         N(x,y,z)%n_r(:) = N(x,y,z-1)%n_r(:)
                F = N(x,y,z)%n_r(:)*g(:) 
                F1 = N(x,y,z-1)%n_r(:)*g(:) 
                F2 = N(x,y,z-2)%n_r(:)*g(:) 
!                N(x,y,z)%n_r(5) = 2.0_rk*N(x,y,z-1)%n_r(5)  - N(x,y,z-2)%n_r(5)
!                N(x,y,z)%n_r(9) = 2.0_rk*N(x,y,z-1)%n_r(9)  - N(x,y,z-2)%n_r(9)
!                N(x,y,z)%n_r(11) = 2.0_rk*N(x,y,z-1)%n_r(11)  - N(x,y,z-2)%n_r(11)
!                N(x,y,z)%n_r(13) = 2.0_rk*N(x,y,z-1)%n_r(13)  - N(x,y,z-2)%n_r(13)
!                N(x,y,z)%n_r(15) = 2.0_rk*N(x,y,z-1)%n_r(15)  - N(x,y,z-2)%n_r(15)
!                N(x,y,z)%n_r(:) = 2.0_rk*N(x,y,z-1)%n_r(:)  - N(x,y,z-2)%n_r(:)
!                N(x,y,z)%n_r(:) = N(x,y,z-1)%n_r(:)
!                F = 2.0_rk*F1-F2
                F = F1
!                F(6) = 2.0_rk*F1(6)-F2(6)
!                F(18) = 2.0_rk*F1(18)-F2(18)
!                F(16) = 2.0_rk*F1(16)-F2(16)
!                F(14) = 2.0_rk*F1(14)-F2(14)
!                F(10) = 2.0_rk*F1(10)-F2(10)
!                F(6)  = F1(6) 
!                F(18) = F1(18)
!                F(16) = F1(16)
!                F(14) = F1(14)
!                F(10) = F1(10)


!
!
!             if (itype==0) then
!                uz=pr
!
!               if (uz.lt.0.0) uz=0.
!
!
!
!                rho = (F(19) + F(1) + F(2) + F(3) + F(4) +  &
!                         F(7) + F(11) + F(12) + F(8) +  &
!                         2.*( F(5) + F(9) + F(13) + F(15) + F(17)))/(uz+1.)
!
!                F(6) = F(5)-uz*rho/3.
!                F(18) = F(15) + (F(7)+F(11)+F(3)-F(8)-F(12)-F(4))/2. &
!                        - uz*rho/6. - uy*rho/2.
!                F(16) = F(17) + (F(8)+F(12)+F(4)-F(7)-F(11)-F(3))/2. &
!                        - uz*rho/6. + uy*rho/2.
!                F(14)  = F(9) + (F(7)+F(8)+F(1)-F(11)-F(12)-F(2))/2. &
!                        - uz*rho/6. - ux*rho/2.
!                F(10) = F(13) + (F(11)+F(12)+F(2)-F(7)-F(8)-F(1))/2. &
!                        - uz*rho/6. + ux*rho/2.
!            else if (itype==1) then
!               F = N(x,y,z)%n_r(:)*g(:) 
!               rho = 1.0_rk !pr
!
!               if(rho>0.) then
!               uz = -1. + (F(19) + &
!                                F(1) + F(2) + F(3) + F(4) +  &
!                                F(7) + F(11) + F(12) + F(8) +  &
!                                2.*( F(5) + F(9) + F(13) + F(15) + F(17)))/rho
!               else 
!                 uz = 0.
!               endif                                       
!
!               F(6) = F(5)-rho*uz/3.
!                 F(18) = F(15) + (F(7)+F(11)+F(3)-F(8)-F(12)-F(4))/2. &
!                         - uz*rho/6. - uy*rho/2.
!                 F(16) = F(17) + (F(8)+F(12)+F(4)-F(7)-F(11)-F(3))/2. &
!                         - uz*rho/6. + uy*rho/2.
!                 F(14)  = F(9) + (F(7)+F(8)+F(1)-F(11)-F(12)-F(2))/2. &
!                         - uz*rho/6. - ux*rho/2.
!                 F(10) = F(13) + (F(11)+F(12)+F(2)-F(7)-F(8)-F(1))/2. &
!                         - uz*rho/6. + ux*rho/2.
!               end if
               N(x,y,z)%n_r(:) = F(:)/g(:) 
         end do
        end do
        endif
        
        ! edges...
!       if(minx.and.miny) then
!       x=1
!       y=1
!
!!       do z=startz,endz
!       do z=1,nz
!
!          F = N(x,y,z)%n_r(:)*g(:) 
!
!          F(1)  = F(2) 
!          F(3)  = F(4) 
!          F(7)  = F(12)
!!          F(9)  = F(14) - (F(5)-F(6))/4. 
!!          F(15) = F(18) - (F(5)-F(6))/4. 
!!          F(10) = F(13) + (F(5)-F(6))/4. 
!!          F(16) = F(17) + (F(5)-F(6))/4.
!          F(9)  = F(13)
!          F(15) = F(17)
!          F(10) = F(14)
!          F(16) = F(18)
!          F(8)  = 0.
!          F(11) = 0.
!          F(19) = 0.
!          F(8) = SUM(F)/22.
!          F(11) = F(8)
!          F(19) = F(8)*12.
!          
!          N(x,y,z)%n_r(:) = F(:)/g(:) 
!       end do
!       endif
!
!
!       if(minx.and.maxy) then
!       x=1
!       y=ny
!
!!       do z=startz,endz
!       do z=1,nz
!          F = N(x,y,z)%n_r(:)*g(:) 
!          F(1)  = F(2) 
!          F(4)  = F(3) 
!          F(8)  = F(11)
!!          F(9)  = F(14) - (F(5)-F(6))/4. 
!!          F(17) = F(16) - (F(5)-F(6))/4. 
!!          F(10) = F(13) + (F(5)-F(6))/4. 
!!          F(18) = F(15) + (F(5)-F(6))/4. 
!          F(9)  = F(13)
!          F(17) = F(15)
!          F(10) = F(14)
!          F(18) = F(16)
!          F(7)  = 0.
!          F(12) = 0.
!          F(19) = 0.
!          F(7) = SUM(F)/22.
!          F(12) = F(7)
!          F(19) = F(7)*12.
!
!          N(x,y,z)%n_r(:) = F(:)/g(:) 
!       end do
!       endif
!
!
!       if(maxx.and.miny) then
!       x=nx
!       y=1
!
!!       do z=startz,endz
!       do z=1,nz
!          F = N(x,y,z)%n_r(:)*g(:) 
!          F(2)  = F(1) 
!          F(3)  = F(4) 
!          F(11) = F(8)
!!          F(13) = F(10) - (F(5)-F(6))/4. 
!!          F(15) = F(18) - (F(5)-F(6))/4. 
!!          F(16) = F(17) + (F(5)-F(6))/4. 
!!          F(14) = F(9)  + (F(5)-F(6))/4. 
!          F(13) = F(9)
!          F(15) = F(17)
!          F(16) = F(18)
!          F(14) = F(10)
!          F(7)  = 0.
!          F(12) = 0.
!          F(19) = 0.
!          F(7) = SUM(F)/22.
!          F(12) = F(7)
!          F(19) = F(7)*12.
!
!          N(x,y,z)%n_r(:) = F(:)/g(:) 
!       end do
!       endif
!
!
!       if(maxx.and.maxy) then
!       x=nx
!       y=ny
!
!!       do z=startz,endz
!       do z=1,nz
!          F = N(x,y,z)%n_r(:)*g(:) 
!          F(2)  = F(1) 
!          F(4)  = F(3) 
!          F(12) = F(7)
!!          F(13) = F(10) - (F(5)-F(6))/4. 
!!          F(17) = F(16) - (F(5)-F(6))/4. 
!!          F(18) = F(15) + (F(5)-F(6))/4. 
!!          F(14) = F(9)  + (F(5)-F(6))/4. 
!          F(13) = F(9)
!          F(17) = F(15)
!          F(18) = F(16)
!          F(14) = F(10)
!          F(8)  = 0.
!          F(11) = 0.
!          F(19) = 0.
!          F(8) = SUM(F)/22.
!          F(11) = F(8)
!          F(19) = F(8)*12.
!
!          N(x,y,z)%n_r(:) = F(:)/g(:) 
!       end do
!       endif


!       if(miny.and.minz) then
!       y=1
!       z=1
!
!       do x=startx,endx
!          F = N(x,y,z)%n_r(:)*g(:) 
!          F(3)  = F(4) 
!          F(5)  = F(6) 
!          F(15) = F(18)
!!          F(7)  = F(12) - (F(1)-F(2))/4. 
!!          F(9)  = F(14) - (F(1)-F(2))/4. 
!!          F(11) = F(8)  + (F(1)-F(2))/4. 
!!          F(13) = F(10) + (F(1)-F(2))/4. 
!          F(7)  = F(8) 
!          F(9)  = F(10)
!          F(11) = F(12)
!          F(13) = F(14)
!          F(16)  = 0.
!          F(17) = 0.
!          F(19) = 0.
!          if (itype==1) then
!            F(16) = (fr-SUM(F))/14.
!          else
!            F(16) = SUM(F)/22.
!          endif
!          F(17) = F(16)
!          F(19) = F(16)*12.
!
!          N(x,y,z)%n_r(:) = F(:)/g(:) 
!       end do
!       endif
!
!
!       if(maxy.and.minz) then
!       y=ny
!       z=1
!
!       do x=startx,endx
!          F = N(x,y,z)%n_r(:)*g(:) 
!          F(4)  = F(3) 
!          F(5)  = F(6) 
!          F(17) = F(16)
!!          F(8)  = F(11) - (F(1)-F(2))/4. 
!!          F(9)  = F(14) - (F(1)-F(2))/4. 
!!          F(12) = F(7)  + (F(1)-F(2))/4. 
!!          F(13) = F(10) + (F(1)-F(2))/4. 
!          F(8)  = F(7) 
!          F(9)  = F(10)
!          F(12) = F(11)
!          F(13) = F(14)
!          F(18) = 0.
!          F(15) = 0.
!          F(19) = 0.
!          if (itype==1) then
!            F(18) = (fr-SUM(F))/14.
!          else
!            F(18) = SUM(F)/22.
!          endif
!          F(15) = F(18)
!          F(19) = F(18)*12.
!
!          N(x,y,z)%n_r(:) = F(:)/g(:) 
!       end do
!       endif
!
!
!       if(miny.and.maxz) then
!       y=1
!       z=nz
!
!       do x=startx,endx
!          F = N(x,y,z)%n_r(:)*g(:) 
!          F(3)  = F(4) 
!          F(6)  = F(5) 
!          F(16) = F(17)
!!          F(7)  = F(12) - (F(1)-F(2))/4. 
!!          F(10) = F(13) - (F(1)-F(2))/4.  
!!          F(11) = F(8)  + (F(1)-F(2))/4. 
!!          F(14) = F(9)  + (F(1)-F(2))/4. 
!          F(7)  = F(8) 
!          F(10) = F(9)
!          F(11) = F(12)
!          F(14) = F(13)
!          F(18) = 0.
!          F(15) = 0.
!          F(19) = 0.
!          if (itype==1) then
!            F(18) = (pr-SUM(F))/14.
!          else
!            F(18) = SUM(F)/22.
!          endif
!          F(15) = F(18)
!          F(19) = F(18)*12.
!
!          N(x,y,z)%n_r(:) = F(:)/g(:) 
!       end do
!       endif
!
!
!       if(maxy.and.maxz) then
!       y=ny
!       z=nz
!
!       do x=startx,endx
!          F = N(x,y,z)%n_r(:)*g(:) 
!          F(4)  = F(3) 
!          F(6)  = F(5) 
!          F(18) = F(15)
!!          F(8)  = F(11) - (F(1)-F(2))/4. 
!!          F(10) = F(13) - (F(1)-F(2))/4.  
!!          F(12) = F(7)  + (F(1)-F(2))/4.  
!!          F(14) = F(9)  + (F(1)-F(2))/4.  
!          F(8)  = F(7) 
!          F(10) = F(9)
!          F(12) = F(11)
!          F(14) = F(13)
!          F(16)  = 0.
!          F(17) = 0.
!          F(19) = 0.
!          if (itype==1) then
!            F(16) = (pr-SUM(F))/14.
!          else
!            F(16) = SUM(F)/22.
!          endif
!          F(17) = F(16)
!          F(19) = F(16)*12.
!
!          N(x,y,z)%n_r(:) = F(:)/g(:) 
!       end do
!       endif
!
!       
!       if(minx.and.minz) then
!       x=1
!       z=1
!
!       do y=starty,endy
!          F = N(x,y,z)%n_r(:)*g(:) 
!          F(1)  = F(2) 
!          F(5)  = F(6) 
!          F(9)  = F(14)
!!          F(7)  = F(12) - (F(3)-F(4))/4. 
!!          F(15) = F(18) - (F(3)-F(4))/4. 
!!          F(8)  = F(11) + (F(3)-F(4))/4. 
!!          F(17) = F(16) + (F(3)-F(4))/4. 
!          F(7)  = F(11) 
!          F(15) = F(16)
!          F(8)  = F(12)
!          F(17) = F(18)
!          F(10)  = 0.
!          F(13) = 0.
!          F(19) = 0.
!          if (itype==1) then
!            F(10) = (fr-SUM(F))/14.
!          else
!            F(10) = SUM(F)/22.
!          endif
!          F(13) = F(10)
!          F(19) = F(10)*12.
!
!          N(x,y,z)%n_r(:) = F(:)/g(:) 
!       end do
!       endif
!
!
!       if(maxx.and.minz) then
!       x=nx
!       z=1
!
!       do y=starty,endy
!          F = N(x,y,z)%n_r(:)*g(:) 
!          F(2)  = F(1) 
!          F(5)  = F(6) 
!          F(13) = F(10)
!!          F(11) = F(8)  - (F(3)-F(4))/4. 
!!          F(15) = F(18) - (F(3)-F(4))/4.  
!!          F(12) = F(7)  + (F(3)-F(4))/4. 
!!          F(17) = F(16) + (F(3)-F(4))/4. 
!          F(11) = F(7) 
!          F(15) = F(16)
!          F(12) = F(8)
!          F(17) = F(18)
!          F(9) = 0.
!          F(14) = 0.
!          F(19) = 0.
!          if (itype==1) then
!            F(9) = (fr-SUM(F))/14.
!          else
!            F(9) = SUM(F)/22.
!          endif
!          F(14) = F(9)
!          F(19) = F(9)*12.
!
!          N(x,y,z)%n_r(:) = F(:)/g(:) 
!       end do
!       endif
!
!
!       if(minx.and.maxz) then
!       x=1
!       z=nz
!
!       do y=starty,endy
!          F = N(x,y,z)%n_r(:)*g(:) 
!          F(1)  = F(2) 
!          F(6)  = F(5) 
!          F(10) = F(13)
!!          F(7)  = F(12) - (F(3)-F(4))/4. 
!!          F(16) = F(17) - (F(3)-F(4))/4. 
!!          F(8)  = F(11) + (F(3)-F(4))/4. 
!!          F(18) = F(15) + (F(3)-F(4))/4. 
!          F(7)  = F(11) 
!          F(16) = F(15)
!          F(8)  = F(12)
!          F(18) = F(17)
!          F(9)  = 0.
!          F(14) = 0.
!          F(19) = 0.
!          if (itype==1) then
!            F(9) = (pr-SUM(F))/14.
!          else
!            F(9) = SUM(F)/22.
!          endif
!          F(14) = F(9)
!          F(19) = F(9)*12.
!
!          N(x,y,z)%n_r(:) = F(:)/g(:) 
!       end do
!       endif
!
!       
!       if(maxx.and.maxz) then
!       x=nx
!       z=nz
!
!       do y=starty,endy
!          F = N(x,y,z)%n_r(:)*g(:) 
!          F(2)  = F(1) 
!          F(6)  = F(5) 
!          F(14) = F(9)
!!          F(11) = F(8)  - (F(3)-F(4))/4. 
!!          F(16) = F(17) - (F(3)-F(4))/4.  
!!          F(12) = F(7)  + (F(3)-F(4))/4.  
!!          F(18) = F(15) + (F(3)-F(4))/4. 
!          F(11) = F(7) 
!          F(16) = F(15)
!          F(12) = F(8)
!          F(18) = F(17)
!          F(10)  = 0.
!          F(13) = 0.
!          F(19) = 0.
!          if (itype==1) then
!            F(10) = (pr-SUM(F))/14.
!          else
!            F(10) = SUM(F)/22.
!          endif
!          F(13) = F(10)
!          F(19) = F(10)*12.
!
!          N(x,y,z)%n_r(:) = F(:)/g(:) 
!       end do
!       endif
!       
!       ! and now corners
!       if(minx.and.miny.and.minz) then
!       x=1
!       y=1
!       z=1
!
!       F = N(x,y,z)%n_r(:)*g(:) 
!       F(1)  = F(2) 
!       F(3)  = F(4) 
!       F(5)  = F(6)
!       F(7)  = F(12)
!       F(9)  = F(14)
!       F(15) = F(18)
!
!       F(8)  = 0.
!       F(11) = 0.
!       F(10) = 0.
!       F(13) = 0.
!       F(16) = 0.
!       F(17) = 0.
!       F(19) = 0.
!
!       if (itype==1) then
!         F(8) = (fr-SUM(F))/18.
!       else
!         F(8) = SUM(F)/18.
!       endif
!       F(11) = F(8)
!       F(10) = F(8)
!       F(13) = F(8)
!       F(16) = F(8)
!       F(17) = F(8)
!       F(19) = F(8)*12.
!       N(x,y,z)%n_r(:) = F(:)/g(:) 
!       endif
!       
!
!       if(maxx.and.miny.and.minz) then
!       x=nx
!       y=1
!       z=1
!
!       F = N(x,y,z)%n_r(:)*g(:) 
!       F(2)  = F(1) 
!       F(3)  = F(4) 
!       F(5)  = F(6)
!       F(11) = F(8)
!       F(13) = F(10)
!       F(15) = F(18)
!
!       F(7)  = 0.
!       F(9)  = 0.
!       F(12) = 0.
!       F(14) = 0.
!       F(16) = 0.
!       F(17) = 0.
!       F(19) = 0.
!
!       if (itype==1) then
!         F(7) = (fr-SUM(F))/18.
!       else
!         F(7) = SUM(F)/18.
!       endif
!       F(9)  = F(7)
!       F(12) = F(7)
!       F(14) = F(7)
!       F(16) = F(7)
!       F(17) = F(7)
!       F(19) = F(7)*12.
!       N(x,y,z)%n_r(:) = F(:)/g(:) 
!       endif
!
!
!       if(minx.and.maxy.and.minz) then
!       x=1
!       y=ny
!       z=1
!
!       F = N(x,y,z)%n_r(:)*g(:) 
!       F(1)  = F(2) 
!       F(4)  = F(3) 
!       F(5)  = F(6)
!       F(8)  = F(11)
!       F(9)  = F(14)
!       F(17) = F(16)
!
!       F(7)  = 0.
!       F(10) = 0.
!       F(12) = 0.
!       F(13) = 0.
!       F(15) = 0.
!       F(18) = 0.
!       F(19) = 0.
!
!       if (itype==1) then
!         F(7) = (fr-SUM(F))/18.
!       else
!         F(7) = SUM(F)/18.
!       endif
!       F(10) = F(7)
!       F(12) = F(7)
!       F(13) = F(7)
!       F(15) = F(7)
!       F(18) = F(7)
!       F(19) = F(7)*12.
!       N(x,y,z)%n_r(:) = F(:)/g(:) 
!       endif
!       
!
!       if(maxx.and.maxy.and.minz) then
!       x=nx
!       y=ny
!       z=1
!
!       F = N(x,y,z)%n_r(:)*g(:) 
!       F(2)  = F(1) 
!       F(4)  = F(3) 
!       F(5)  = F(6)
!       F(12) = F(7)
!       F(13) = F(10)
!       F(17) = F(16)
!
!       F(8)  = 0.
!       F(9)  = 0.
!       F(11) = 0.
!       F(14) = 0.
!       F(15) = 0.
!       F(18) = 0.
!       F(19) = 0.
!
!       if (itype==1) then
!         F(8) = (fr-SUM(F))/18.
!       else
!         F(8) = SUM(F)/18.
!       endif
!       F(9)  = F(8)
!       F(11) = F(8)
!       F(14) = F(8)
!       F(15) = F(8)
!       F(18) = F(8)
!       F(19) = F(8)*12.
!       N(x,y,z)%n_r(:) = F(:)/g(:) 
!       endif
!
!
!       if(minx.and.miny.and.maxz) then
!       x=1
!       y=1
!       z=nz
!
!       F = N(x,y,z)%n_r(:)*g(:) 
!       F(1)  = F(2) 
!       F(3)  = F(4) 
!       F(6)  = F(5)
!       F(10) = F(13)
!       F(7)  = F(12)
!       F(16) = F(17)
!
!       F(8)  = 0.
!       F(9)  = 0.
!       F(11) = 0.
!       F(15) = 0.
!       F(14) = 0.
!       F(18) = 0.
!       F(19) = 0.
!
!       if (itype==1) then
!         F(8) = (pr-SUM(F))/18.
!       else
!         F(8) = SUM(F)/18.
!       endif
!       F(9)  = F(8)
!       F(11) = F(8)
!       F(14) = F(8)
!       F(15) = F(8)
!       F(18) = F(8)
!       F(19) = F(8)*12.
!       N(x,y,z)%n_r(:) = F(:)/g(:) 
!       endif
!       
!
!       if(maxx.and.miny.and.maxz) then
!       x=nx
!       y=1
!       z=nz
!
!       F = N(x,y,z)%n_r(:)*g(:) 
!       F(2)  = F(1) 
!       F(3)  = F(4) 
!       F(6)  = F(5)
!       F(11) = F(8)
!       F(14) = F(9)
!       F(16) = F(17)
!
!       F(7)  = 0.
!       F(10) = 0.
!       F(12) = 0.
!       F(13) = 0.
!       F(15) = 0.
!       F(18) = 0.
!       F(19) = 0.
!
!       if (itype==1) then
!         F(7) = (pr-SUM(F))/18.
!       else
!         F(7) = SUM(F)/18.
!       endif
!       F(10)  = F(7)
!       F(12) = F(7)
!       F(13) = F(7)
!       F(15) = F(7)
!       F(18) = F(7)
!       F(19) = F(7)*12.
!       N(x,y,z)%n_r(:) = F(:)/g(:) 
!       endif
!
!
!       if(minx.and.maxy.and.maxz) then
!       x=1
!       y=ny
!       z=nz
!
!       F = N(x,y,z)%n_r(:)*g(:) 
!       F(1)  = F(2) 
!       F(4)  = F(3) 
!       F(6)  = F(5)
!       F(8)  = F(11)
!       F(10) = F(13)
!       F(18) = F(15)
!
!       F(7)  = 0.
!       F(9)  = 0.
!       F(12) = 0.
!       F(14) = 0.
!       F(16) = 0.
!       F(17) = 0.
!       F(19) = 0.
!
!       if (itype==1) then
!         F(7) = (pr-SUM(F))/18.
!       else
!         F(7) = SUM(F)/18.
!       endif
!       F(9)  = F(7)
!       F(12) = F(7)
!       F(14) = F(7)
!       F(16) = F(7)
!       F(17) = F(7)
!       F(19) = F(7)*12.
!       N(x,y,z)%n_r(:) = F(:)/g(:) 
!       endif
!
!
!            
!       if(maxx.and.maxy.and.maxz) then
!       x=nx
!       y=ny
!       z=nz
!
!       F = N(x,y,z)%n_r(:)*g(:) 
!       F(2)  = F(1) 
!       F(4)  = F(3) 
!       F(6)  = F(5)
!       F(12) = F(7)
!       F(14) = F(9)
!       F(18) = F(15)
!
!       F(8)  = 0.
!       F(10) = 0.
!       F(11) = 0.
!       F(13) = 0.
!       F(16) = 0.
!       F(17) = 0.
!       F(19) = 0.
!
!       if (itype==1) then
!         F(8) = (pr-SUM(F))/18.
!       else
!         F(8) = SUM(F)/18.
!       endif
!       F(10)  = F(8)
!       F(11) = F(8)
!       F(13) = F(8)
!       F(16) = F(8)
!       F(17) = F(8)
!       F(19) = F(8)*12.
!       N(x,y,z)%n_r(:) = F(:)/g(:)
!       endif
end subroutine lbe_invade_channel

!> a separate subroutine for inv_fluid=18 and inv_type=3
!>
!> the flow is driven by gravity, the channel walls are with
!> diagonal stripes of smoothly changing slip parameter.
!> \param[in,out] N local fluid populations
!> \param[in] pr is the magnitude of the k-vector for the pattern
!> \param[in] pb is the magnitude of the changes
!> \param[in] pg is the mean value of the slip parameter
subroutine lbe_invade_slip_channel(N,pr,pb,pg)
	implicit none
        real*8, intent(in) :: pr,pb,pg
	type(lbe_site),dimension(0:,0:,0:), intent(inout)  :: N   ! lattice
        real*8 :: F(19)
        real*8 :: s
	integer :: x,y,z

        ! partial slip boundaries

        x=1
           do y=2,ny-1
             do z=1,nz
                  F = N(x,y,z)%n_r(:)*g(:)

                  if(x>=y) then
                    s=cos(2.*pi*real(2-x-y+z)/pr)*pb*0.5+pg
                  else
                    s=cos(2.*pi*real(y+x-2+z)/pr)*pb*0.5+pg
                  endif

                  F(1)  = F(2) 
                  F(9) = (F(14) - (F(17)+F(15)+F(5)-F(18)-F(16)-F(6))/2.)*(1.-s)+s*F(13)
                  F(10) = (F(13) + (F(17)+F(15)+F(5)-F(18)-F(16)-F(6))/2.)*(1.-s)+s*F(14)
                  F(7)  = (F(12) - (F(15)+F(16)+F(3)-F(17)-F(18)-F(4))/2.)*(1.-s)+s*F(11)
                  F(8) =  (F(11) + (F(15)+F(16)+F(3)-F(17)-F(18)-F(4))/2.)*(1.-s)+s*F(12)
                  
                 N(x,y,z)%n_r(:) = F(:)/g(:) 
	    end do
         end do

         x=nx
         do y=2,ny-1
            do z=1,nz

                F = N(x,y,z)%n_r(:)*g(:) 

                  if(x>=y) then
                    s=cos(2.*pi*real(2-x-y+z)/pr)*pb*0.5+pg
                    else 
                    s=cos(2.*pi*real(y+x-2+z)/pr)*pb*0.5+pg
                  endif
               F(2)  = F(1) 
               F(14) = (F(9) + (F(17)+F(15)+F(5)-F(18)-F(16)-F(6))/2.)*(1.-s)+s*F(10)
               F(13) = (F(10) - (F(17)+F(15)+F(5)-F(18)-F(16)-F(6))/2.)*(1.-s)+s*F(9)
               F(12) = (F(7) + (F(15)+F(16)+F(3)-F(17)-F(18)-F(4))/2.)*(1.-s)+s*F(8)
               F(11) = (F(8) - (F(15)+F(16)+F(3)-F(17)-F(18)-F(4))/2.)*(1.-s)+s*F(7)
                  
               N(x,y,z)%n_r(:) = F(:)/g(:) 
            end do
         end do

        y=1
           do x=2,nx-1
             do z=1,nz
                  F = N(x,y,z)%n_r(:)*g(:) 

                  if(x>=y) then
                    s=cos(2.*pi*real(2-x-y+z)/pr)*pb*0.5+pg
                    else 
                    s=cos(2.*pi*real(y+x-2+z)/pr)*pb*0.5+pg
                  endif

                  F(3)  = F(4) 
                  F(7) = (F(12) - (F(1)+F(9)+F(10)-F(2)-F(13)-F(14))/2.)*(1.-s)+s*F(8)
                  F(11)  = (F(8) - (F(2)+F(13)+F(14)-F(1)-F(9)-F(10))/2.)*(1.-s)+s*F(12)
                  F(15) = (F(18) - (F(5)+F(9)+F(13)-F(6)-F(14)-F(10))/2.)*(1.-s)+s*F(17)
                  F(16) =  (F(17) - (F(6)+F(10)+F(14)-F(5)-F(9)-F(13))/2.)*(1.-s)+s*F(18)

                 N(x,y,z)%n_r(:) = F(:)/g(:) 
	    end do
         end do

         y=ny
         do x=2,nx-1
            do z=1,nz
               F = N(x,y,z)%n_r(:)*g(:) 

                  if(x>=y) then
                    s=cos(2.*pi*real(2-x-y+z)/pr)*pb*0.5+pg
                    else 
                    s=cos(2.*pi*real(y+x-2+z)/pr)*pb*0.5+pg
                  endif

                  F(4)  = F(3) 
                  F(12) = (F(7) + (F(1)+F(9)+F(10)-F(2)-F(13)-F(14))/2. )*(1.-s)+s*F(11)
                  F(8)  = (F(11) + (F(2)+F(13)+F(14)-F(1)-F(9)-F(10))/2.)*(1.-s)+s*F(7)
                  F(18) = (F(15) + (F(5)+F(9)+F(13)-F(6)-F(14)-F(10))/2.)*(1.-s)+s*F(16)
                  F(17) = (F(16) + (F(6)+F(10)+F(14)-F(5)-F(9)-F(13))/2.)*(1.-s)+s*F(15)
                  
               N(x,y,z)%n_r(:) = F(:)/g(:) 
            end do
         end do

        
        ! edges...

       x=1
       y=1

       do z=1,nz

          F = N(x,y,z)%n_r(:)*g(:) 

                  if(x>=y) then
                    s=cos(2.*pi*real(2-x-y+z)/pr)*pb*0.5+pg
                    else 
                    s=cos(2.*pi*real(y+x-2+z)/pr)*pb*0.5+pg
                  endif

          F(1)  = F(2) 
          F(3)  = F(4) 
          F(7)  = F(12)
          F(9)  = (F(14) - (F(5)-F(6))/4.)*(1.-s)+s*F(13)
          F(15) = (F(18) - (F(5)-F(6))/4.)*(1.-s)+s*F(17)
          F(10) = (F(13) + (F(5)-F(6))/4.)*(1.-s)+s*F(14)
          F(16) = (F(17) + (F(5)-F(6))/4.)*(1.-s)+s*F(18)
          F(8)  = 0.
          F(11) = 0.
          F(19) = 0.
          F(8) = SUM(F)/22.
          F(11) = F(8)
          F(19) = F(8)*12.
          
          N(x,y,z)%n_r(:) = F(:)/g(:) 
       end do


       x=1
       y=ny

       do z=1,nz
          F = N(x,y,z)%n_r(:)*g(:) 

                  if(x>=y) then
                    s=cos(2.*pi*real(2-x-y+z)/pr)*pb*0.5+pg
                    else 
                    s=cos(2.*pi*real(y+x-2+z)/pr)*pb*0.5+pg
                  endif

          F(1)  = F(2) 
          F(4)  = F(3) 
          F(8)  = F(11)
          F(9)  = (F(14) - (F(5)-F(6))/4. )*(1.-s)+s*F(13)
          F(17) = (F(16) - (F(5)-F(6))/4. )*(1.-s)+s*F(15)
          F(10) = (F(13) + (F(5)-F(6))/4. )*(1.-s)+s*F(14)
          F(18) = (F(15) + (F(5)-F(6))/4. )*(1.-s)+s*F(16)
          F(7)  = 0.
          F(12) = 0.
          F(19) = 0.
          F(7) = SUM(F)/22.
          F(12) = F(7)
          F(19) = F(7)*12.

          N(x,y,z)%n_r(:) = F(:)/g(:) 
       end do

       x=nx
       y=1

       do z=1,nz

                  if(x>=y) then
                    s=cos(2.*pi*real(2-x-y+z)/pr)*pb*0.5+pg
                    else 
                    s=cos(2.*pi*real(y+x-2+z)/pr)*pb*0.5+pg
                  endif

          F = N(x,y,z)%n_r(:)*g(:) 
          F(2)  = F(1) 
          F(3)  = F(4) 
          F(11) = F(8)
          F(13) = (F(10) - (F(5)-F(6))/4. )*(1.-s)+s*F(9)
          F(15) = (F(18) - (F(5)-F(6))/4. )*(1.-s)+s*F(17)
          F(16) = (F(17) + (F(5)-F(6))/4. )*(1.-s)+s*F(18)
          F(14) = (F(9)  + (F(5)-F(6))/4. )*(1.-s)+s*F(10)
          F(7)  = 0.
          F(12) = 0.
          F(19) = 0.
          F(7) = SUM(F)/22.
          F(12) = F(7)
          F(19) = F(7)*12.

          N(x,y,z)%n_r(:) = F(:)/g(:) 
       end do

       x=nx
       y=ny

       do z=1,nz
          F = N(x,y,z)%n_r(:)*g(:) 

                  if(x>=y) then
                    s=cos(2.*pi*real(2-x-y+z)/pr)*pb*0.5+pg
                    else 
                    s=cos(2.*pi*real(y+x-2+z)/pr)*pb*0.5+pg
                  endif

          F(2)  = F(1) 
          F(4)  = F(3) 
          F(12) = F(7)
          F(13) = (F(10) - (F(5)-F(6))/4. )*(1.-s)+s*F(9)
          F(17) = (F(16) - (F(5)-F(6))/4. )*(1.-s)+s*F(15)
          F(18) = (F(15) + (F(5)-F(6))/4. )*(1.-s)+s*F(16)
          F(14) = (F(9)  + (F(5)-F(6))/4. )*(1.-s)+s*F(10)
          F(8)  = 0.
          F(11) = 0.
          F(19) = 0.
          F(8) = SUM(F)/22.
          F(11) = F(8)
          F(19) = F(8)*12.

          N(x,y,z)%n_r(:) = F(:)/g(:)
       end do
end subroutine lbe_invade_slip_channel

!> fills the sites in the first halo layer at the global maximum
!> z-boundary with an extrapolation based on the two last real layers
!>
!> The rock state is copied from the last real layer. Also the first
!> layer of the x- and y-halo are extrapolated.
!>
!> \param[in,out] whole_N local chunk of the lattice with full
!> halo of depth \c halo_extent
subroutine extrapolate_maxz(whole_N)
    type(lbe_site),intent(inout) :: &
         &whole_N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
    integer :: x,y

    if (start(3)>=tnz-nz) then
       do x=0,nx+1
          do y=0,ny+1
             whole_N(x,y,nz+1)%n_r &
                  &= 2.0_rk*whole_N(x,y,nz)%n_r-whole_N(x,y,nz-1)%n_r
#ifndef SINGLEFLUID
             whole_N(x,y,nz+1)%n_b &
                  &= 2.0_rk*whole_N(x,y,nz)%n_b-whole_N(x,y,nz-1)%n_b
#endif
#ifndef NOSURFACTANT
             whole_N(x,y,nz+1)%n_s &
                  &= 2.0_rk*whole_N(x,y,nz)%n_s-whole_N(x,y,nz-1)%n_s
             whole_N(x,y,nz+1)%d &
                  &= 2.0_rk*whole_N(x,y,nz)%d-whole_N(x,y,nz-1)%d
#endif
             whole_N(x,y,nz+1)%rock_state = whole_N(x,y,nz)%rock_state
          end do
       end do
    end if
end subroutine extrapolate_maxz

!> fills the rock state in the first halo layer at the global maximum
!> z-boundary with the values from the second last real layer
!>
!> Also the first layer of the x- and y-halo is filled.
!>
!> \param[in,out] whole_N local chunk of the lattice with full
!> halo of depth \c halo_extent
subroutine mirror_rock_maxz(whole_N)
    type(lbe_site),intent(inout) :: &
         &whole_N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
    integer :: x,y

    if (start(3)>=tnz-nz) then
       do x=0,nx+1
          do y=0,ny+1
             whole_N(x,y,nz+1)%rock_state = whole_N(x,y,nz-1)%rock_state
          end do
       end do
    end if
end subroutine mirror_rock_maxz

!> fills the rock state in the first halo layer at the global minimum
!> z-boundary with the values from the second real layer
!>
!> Also the first layer of the x- and y-halo is filled.
!>
!> \param[in,out] whole_N local chunk of the lattice with full
!> halo of depth \c halo_extent
subroutine mirror_rock_minz(whole_N)
    type(lbe_site),intent(inout) :: &
         &whole_N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
    integer :: x,y

    if (start(3)==1) then
       do x=0,nx+1
          do y=0,ny+1
             whole_N(x,y,0)%rock_state = whole_N(x,y,2)%rock_state
          end do
       end do
    end if
end subroutine mirror_rock_minz

!> \name per-node on-site velocity and density boundary condition
!> \{
!> The implementation closely follows Hecht&Harting:
!> J. Stat. Mech. 2010, P01018 (2010) which is based on
!> Zou&He. Different from the earlier implementation the routines here
!> can be applied to single vectors of distribution components.
!>
!> The two-letter suffix of the names of the subroutines specifies the
!> axis normal to the associated boundary plane (x, y, z) and the
!> direction along this axis towards which the boundary is considered
!> open (p, n) for ("positive", "negative").

!> Prescribes the velocity for a single set of distributions. The
!> boundary is assumed to be open in negative z direction.
!>
!> \param[in,out] n local distribution vector, only unknown elements
!> are modified beyond numerical accuracy
!>
!> \param[in] u desired velocity vector
subroutine zou_he_hecht_velocity_dist_zn(n,u)
    real(kind=rk),intent(in) :: u(3)
    real(kind=rk),intent(inout) :: n(nvecs)
    real(kind=rk) :: rho

    n = n*g                     ! sanitize distribution space
    rho = (n(1)+n(2)+n(3)+n(4)+n(7)+n(8)+n(11)+n(12)+n(19)&
         &+2.0_rk*(n(6)+n(10)+n(14)+n(16)+n(18)))/(1.0_rk-u(3))

    call zou_he_hecht_dist_zn(n,rho*u)
    n = n/g                     ! de-sanitize again...
end subroutine zou_he_hecht_velocity_dist_zn

!> Prescribes the mass flow for a single set of distributions. The
!> boundary is assumed to be open in negative z direction.
!>
!> \param[in,out] n local distribution vector, only unknown elements
!> are modified beyond numerical accuracy
!>
!> \param[in] p desired mass flow vector
subroutine zou_he_hecht_massflow_dist_zn(n,p)
    real(kind=rk),intent(in) :: p(3)
    real(kind=rk),intent(inout) :: n(nvecs)

    n = n*g                     ! sanitize distribution space
    call zou_he_hecht_dist_zn(n,p)
    n = n/g                     ! de-sanitize again...
end subroutine zou_he_hecht_massflow_dist_zn

!> Prescribes the mass flow for a single set of distributions. The
!> boundary is assumed to be open in positive z direction.
!>
!> \param[in,out] n local distribution vector, only unknown elements
!> are modified beyond numerical accuracy
!>
!> \param[in] p desired mass flow vector
subroutine zou_he_hecht_massflow_dist_zp(n,p)
    real(kind=rk),intent(in) :: p(3)
    real(kind=rk),intent(inout) :: n(nvecs)

    n = n*g                     ! sanitize distribution space
    call zou_he_hecht_dist_zp(n,p)
    n = n/g                     ! de-sanitize again...
end subroutine zou_he_hecht_massflow_dist_zp

!> Prescribes the velocity for a single set of distributions. The
!> boundary is assumed to be open in positive z direction.
!>
!> \param[in,out] n local distribution vector, only unknown elements
!> are modified beyond numerical accuracy
!>
!> \param[in] u desired velocity vector
subroutine zou_he_hecht_velocity_dist_zp(n,u)
    real(kind=rk),intent(in) :: u(3)
    real(kind=rk),intent(inout) :: n(nvecs)
    real(kind=rk) :: rho

    n = n*g                     ! sanitize distribution space
    rho = (n(1)+n(2)+n(3)+n(4)+n(7)+n(8)+n(11)+n(12)+n(19)&
         &+2.0_rk*(n(5)+n(9)+n(13)+n(15)+n(17)))/(1.0_rk+u(3))

    call zou_he_hecht_dist_zp(n,rho*u)
    n = n/g                     ! de-sanitize again...
end subroutine zou_he_hecht_velocity_dist_zp

!> Prescribes the density for a single set of distributions. The
!> boundary is assumed to be open in negative z direction.
!>
!> The velocity in the xy plane is assumed to be zero.
!>
!> \param[in,out] n local distribution vector, only unknown elements
!> are modified beyond numerical accuracy
!>
!> \param[in] u desired velocity vector
subroutine zou_he_hecht_density_dist_zn(n,rho)
    real(kind=rk),intent(in) :: rho
    real(kind=rk),intent(inout) :: n(nvecs)
    real(kind=rk) :: rho_u(3)

    n = n*g                     ! sanitize distribution space
    rho_u(1:2) = 0.0_rk
    rho_u(3) = rho - (n(1)+n(2)+n(3)+n(4)+n(7)+n(8)+n(11)+n(12)+n(19)&
         &+2.0_rk*(n(6)+n(10)+n(14)+n(16)+n(18)))

    call zou_he_hecht_dist_zn(n,rho_u)
    n = n/g                     ! de-sanitize again...
end subroutine zou_he_hecht_density_dist_zn

!> Prescribes the density for a single set of distributions. The
!> boundary is assumed to be open in positive z direction.
!>
!> The velocity in the xy plane is assumed to be zero.
!>
!> \param[in,out] n local distribution vector, only unknown elements
!> are modified beyond numerical accuracy
!>
!> \param[in] u desired velocity vector
subroutine zou_he_hecht_density_dist_zp(n,rho)
    real(kind=rk),intent(in) :: rho
    real(kind=rk),intent(inout) :: n(nvecs)
    real(kind=rk) :: rho_u(3)

    n = n*g                     ! sanitize distribution space
    rho_u(1:2) = 0.0_rk
    rho_u(3) = n(1)+n(2)+n(3)+n(4)+n(7)+n(8)+n(11)+n(12)+n(19)&
         &+2.0_rk*(n(5)+n(9)+n(13)+n(15)+n(17)) - rho

    call zou_he_hecht_dist_zp(n,rho_u)
    n = n/g                     ! de-sanitize again...
end subroutine zou_he_hecht_density_dist_zp

!> Actually creates the unknown sanitized populations. The boundary is
!> assumed to be open in negative z direction.
!>
!> Number density \c rho and velocity \c u enter this subroutine as
!> product as they do in the expressions in the article of Hecht&Harting
!> (2010).
!>
!> \warning Both input and output distributions are considered to be
!> without the degeneracy vector \c g, such as can be obtained from
!> \c N(x,y,z)%n_r*g.
!>
!> \param[in,out] n local distribution vector, only unknown elements
!> are modified beyond numerical accuracy. These must be sane
!> distributions without \c g!
!>
!> \param[in] rho_u product of desired density and desired velocity vector
subroutine zou_he_hecht_dist_zn(n,rho_u)
    real(kind=rk),intent(in) :: rho_u(3)
    real(kind=rk),intent(inout) :: n(nvecs)
    real(kind=rk) :: nxzp,nyzp

    nxzp = 0.5_rk*(n(1)+n(7)+n(8)-n(2)-n(11)-n(12))
    nyzp = 0.5_rk*(n(3)+n(7)+n(11)-n(4)-n(8)-n(12))

    n(5) = n(6)+rho_u(3)/3.0_rk
    n(9) = n(14)+rho_u(3)/6.0_rk+0.5_rk*rho_u(1)-nxzp
    n(13) = n(10)+rho_u(3)/6.0_rk-0.5_rk*rho_u(1)+nxzp
    n(15) = n(18)+rho_u(3)/6.0_rk+0.5_rk*rho_u(2)-nyzp
    n(17) = n(16)+rho_u(3)/6.0_rk-0.5_rk*rho_u(2)+nyzp
end subroutine zou_he_hecht_dist_zn

!> Actually creates the unknown sanitized populations. The boundary is
!> assumed to be open in positive z direction.
!>
!> Number density \c rho and velocity \c u enter this subroutine as
!> product as they do in the expressions in the article of Hecht&Harting
!> (2010).
!>
!> \warning Both input and output distributions are considered to be
!> without the degeneracy vector \c g, such as can be obtained from
!> \c N(x,y,z)%n_r*g.
!>
!> \param[in,out] n local distribution vector, only unknown elements
!> are modified beyond numerical accuracy. These must be sane
!> distributions without \c g!
!>
!> \param[in] rho_u product of desired density and desired velocity vector
subroutine zou_he_hecht_dist_zp(n,rho_u)
    real(kind=rk),intent(in) :: rho_u(3)
    real(kind=rk),intent(inout) :: n(nvecs)
    real(kind=rk) :: nxzp,nyzp

    nxzp = 0.5_rk*(n(1)+n(7)+n(8)-n(2)-n(11)-n(12))
    nyzp = 0.5_rk*(n(3)+n(7)+n(11)-n(4)-n(8)-n(12))

    n(6) = n(5)-rho_u(3)/3.0_rk
    n(10) = n(13)-rho_u(3)/6.0_rk+0.5_rk*rho_u(1)-nxzp
    n(14) = n(9)-rho_u(3)/6.0_rk-0.5_rk*rho_u(1)+nxzp
    n(16) = n(17)-rho_u(3)/6.0_rk+0.5_rk*rho_u(2)-nyzp
    n(18) = n(15)-rho_u(3)/6.0_rk-0.5_rk*rho_u(2)+nyzp
end subroutine zou_he_hecht_dist_zp
!> \}

end module lbe_invasion_module
