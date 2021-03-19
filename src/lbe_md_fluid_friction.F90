#include "lbe.h"

!> coupling of particles and fluid by a friction force like that by
!> Ahlrichs/Dunweg in Int. J. Mod. Phys. C 9, 1429 (1998)
module lbe_md_fluid_friction_module
#ifdef MD
    use lbe_force_interface_module, only: add_force_to_all
    use lbe_globals_module, only: g,halo_extent,lp_sur,n_lp,opp_lp_sur,pi&
         &,input_dfile_unit,use_lbe_force,myrankc
    use lbe_helper_module, only: local_coordinates
    use lbe_log_module
    use lbe_md_globals_module
    use lbe_md_helper_module, only: fluid_velocity_and_viscosity&
         &,log_msg_md,error_md,log_msg_md_hdr
    use lbe_parallel_module, only: comm_cart
    use lbe_parms_module, only: amass_r,amass_b,amass_s,inp_file&
         &,arg_input_dfile_set,arg_input_dfile
    use lbe_types_module, only: lbe_site

    implicit none
    include 'mpif.h'
    private

    public fluid_f_interaction_friction,input_fluid_friction&
         &,setup_fluid_friction

    !> assumed Stokes radius for drag force calculation
    real(kind=rk),save :: R_stokes=0.05_rk
    !> save computation: \f$6\pi R_\mathrm{stokes}\f$
    real(kind=rk),save :: Rpi6

    namelist /md_fluid_friction/ R_stokes

contains

    !> read  \c /md_fluid_friction/
    subroutine input_fluid_friction
        integer ierror

        call log_msg_md_hdr("Reading MD F Friction input")

        if (myrankc.eq.0) then
           open (unit=md_input_file_unit,file=trim(inp_file)//'.md',err=100)
           read (unit=md_input_file_unit,nml=md_fluid_friction,err=100)
           close (unit=md_input_file_unit,err=100)
           !call log_msg_md('read /md_fluid_friction/ from file "'//trim(inp_file)&
           !     &//'.md"',.false.)
           !write (6,nml=md_fluid_friction)
        end if

        if ( arg_input_dfile_set ) then
          call log_msg_md("  Getting differential input...")
          open(UNIT = input_dfile_unit, FILE = arg_input_dfile, STATUS = 'UNKNOWN')
          read(UNIT = input_dfile_unit, NML = md_fluid_friction, IOSTAT = ierror)
          if (ierror .ne. 0) then
            call log_msg_md("    WARNING: Differential namelist not found or errors encountered.")
          endif
          close(UNIT = input_dfile_unit)
          call log_ws()
        end if

        write(msgstr,"('R_stokes           = ',F16.10)") R_stokes
        call log_msg(msgstr)
        call log_ws()

        call MPI_Bcast(R_stokes,1,MPI_REAL8,0,comm_cart,ierror)

        return
100     continue
        call error_md('Error reading md input file "'//trim(inp_file)//'.md"')
    end subroutine input_fluid_friction

    !> set \c use_lbe_force and some constants
    subroutine setup_fluid_friction
        use_lbe_force = .true.
        halo_exchange_before_md_run = .true.

        ! 6*pi*R_stokes
        Rpi6 = 6.0_rk*pi*R_stokes
    end subroutine setup_fluid_friction

    !> add dissipative forces in \c P(:)%f(:) and \c lbe_force(:,:,:,:,:)
    subroutine fluid_f_interaction_friction(N)
        type(lbe_site),intent(in) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        real(kind=rk) :: u(3) ! fluid velocity field at particle position
        real(kind=rk) :: mu ! fluid dynamic viscosity at particle position
        real(kind=rk) :: f(3)   ! friction force on md particle
        integer i,ii
        logical stuck

        i = atompnt
        particles: do ii = 1,nlocal
           ! get  u(:)  and  mu
           call fluid_velocity_and_viscosity&
                &(N,(/1.0_rk,1.0_rk,1.0_rk/),P(i)%x,u,mu,stuck)

           ! abort on the first stuck particle
           if (stuck) then
              print '(A,I10,A,3ES10.2,A,3ES10.2,A)'&
                   &,'WARNING: lost a particle in rock: uid='&
                   &,P(i)%uid,',x=(',P(i)%x(:),'),v=(',P(i)%v(:),')'
              call error_md('program exits...')
           end if

           ! Stokes friction: F=6 pi R mu deltav
           f(:) = Rpi6*mu*(u(:)-P(i)%v(:))

           ! force on particle
           P(i)%f(:) = P(i)%f(:) + f(:)

           ! force on fluid
           call apply_lbe_force(N,P(i)%x,-f)

           i = list(i)
        enddo particles
    end subroutine fluid_f_interaction_friction

    !> Apply the force \c f at (global) position \c x of \c
    !> lbe_force(:,:,:,:,:).
    !>
    !> \c x must be at some place that is
    !> covered by the local \c N(:,:,:) , however, this subroutine is
    !> smart enough to deal with pbc.
    subroutine apply_lbe_force(N,x,f)
        type(lbe_site),intent(in) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        real(kind=rk),intent(in) :: x(3),f(3)
        real(kind=rk) :: xx(3),pweight,ff(3)
        integer l,lx(3),lp(3),opp_lp(3)

        call local_coordinates(x(:),xx(:))
        lx(:) = floor(xx(:))

        ! distribute the force among surrounding lattice points according to
        ! trilinear interpolation
        lattice_points: do l=1,n_lp
           lp(:) = lx(:) + lp_sur(:,l)
           ! wall can be though of having infinite mass, so a force has no
           ! effect on it - this saves some operations
           lp_no_rock: if (N(lp(1),lp(2),lp(3))%rock_state==0.0_rk) then
              ! The weight for  lp(:)  is the volume of the hypercube spanned
              ! by  x(:)  and the opposite lattice point  opp_lp(:) . Because
              ! the total volume in between all surrounding lattice points is
              ! 1 no normalization is needed.
              opp_lp(:) = lx(:) + opp_lp_sur(:,l)
              pweight = abs(product(real(opp_lp(:),kind=rk)-xx(:)))

              ! Distribute the force among lbe particle species according to
              ! the total mass of all particles of each species at  lp(:) .
              ! This way, the acceleration for every species is the same.
              ff = f*pweight
              call add_force_to_all(ff(1),ff(2),ff(3),lp(1),lp(2),lp(3))
           end if lp_no_rock
        end do lattice_points
    end subroutine apply_lbe_force

#endif
end module lbe_md_fluid_friction_module
