#include "lbe.h"

!> Lennard-Jones potential between particles and rock sites at the
!> surface of walls
module lbe_md_rock_friction_module
#ifdef MD
    use lbe_bc_module, only: periodically_wrapped
    use lbe_globals_module, only: halo_extent,tsize, input_dfile_unit,myrankc,pi
    use lbe_helper_module, only: local_coordinates
    use lbe_log_module
    use lbe_md_boundary_condition_module, only: local_chunks,local_chunk_type&
         &,n_max_local_chunks
    use lbe_md_fluid_ladd_module, only: lubrication_force_rock&
         &,rc_lubrication_rock
    use lbe_md_fluid_ladd_parms_module, only: lubrication
    use lbe_md_globals_module
    use lbe_md_helper_module, only: log_msg_md,error_md,log_msg_md_hdr&
        &, fluid_velocity_and_viscosity
    use lbe_parallel_module, only: comm_cart
    use lbe_parms_module, only: inp_file,arg_input_dfile_set,arg_input_dfile
    use lbe_types_module, only: lbe_site

    implicit none
    include 'mpif.h'
    private

    public input_rock_friction&
         &,rock_f_interaction_friction,setup_rock_friction

    ! lennard jones parameters
    real(kind=rk),save :: sigma=2.0_rk

    real(kind=rk),save :: r_cut         ! lennard jones cutoff radius
    real(kind=rk),save :: r_cut_sq      ! r_cut**2
     !> assumed Stokes radius for drag force calculation
    real(kind=rk),save :: R_stokes=0.05_rk
    !> Threshold velocity when particle begin to move
    real(kind=rk),save :: threshold_v=0.0_rk
        !> Threshold velocity when particle begin to move
    real(kind=rk),save :: threshold_f=0.0_rk
    !> save computation: \f$6\pi R_\mathrm{stokes}\f$
    real(kind=rk),save :: Rpi6
    logical,save,public :: md_rock_fric_x=.false.
    logical,save,public :: md_rock_fric_y=.false.
    logical,save,public :: md_rock_fric_z=.false.
    logical,save,public :: freeze_all=.false.
 
    namelist /md_rock_friction/ R_stokes, sigma, md_rock_fric_x&
				&,md_rock_fric_y, md_rock_fric_z,threshold_v, threshold_f&
            &, freeze_all

contains

    !> read namelist  /md_rock_friction/  from input file
    subroutine input_rock_friction
        integer ierror

        call log_msg_md_hdr("Reading MD R friction input")

        ! These features are essential for this rock module. They are
        ! enabled here because in  setup_rock_friction()  it would be too late.
        collect_forces = .true.

        if (myrankc.eq.0) then
           open (unit=md_input_file_unit,file=trim(inp_file)//'.md',err=100)
           read (unit=md_input_file_unit,nml=md_rock_friction,err=100)
           close (unit=md_input_file_unit,err=100)
           !call log_msg_md('read /md_rock_friction/ from file "'//trim(inp_file)&
           !     &//'.md"',.false.)
           !write (6,nml=md_rock_friction)
        end if

        if ( arg_input_dfile_set ) then
          call log_msg_md("  Getting differential input...")
          open(UNIT = input_dfile_unit, FILE = arg_input_dfile, STATUS = 'UNKNOWN')
          read(UNIT = input_dfile_unit, NML = md_rock_friction, IOSTAT = ierror)
          if (ierror .ne. 0) then
              call log_msg_md("    WARNING: Differential namelist not found or errors encountered.")
          end if
          close(UNIT = input_dfile_unit)
          call log_ws()
        end if

		  write(msgstr,"('R_stokes            = ',F16.10)") R_stokes
          call log_msg(msgstr)
        write(msgstr,"('threshold_v         = ',F16.10)") threshold_v
            call log_msg(msgstr)
        write(msgstr,"('threshold_f         = ',F16.10)") threshold_f
            call log_msg(msgstr)
        write(msgstr,"('sigma              = ',F16.10)") sigma
        call log_msg(msgstr)
        !> Friction force direction
        write(msgstr,"('md_rock_fric_x      = ',L1)") md_rock_fric_x
            call log_msg(msgstr)
        write(msgstr,"('md_rock_fric_y      = ',L1)") md_rock_fric_y
            call log_msg(msgstr)
        write(msgstr,"('md_rock_fric_z      = ',L1)") md_rock_fric_z
              call log_msg(msgstr)
        write(msgstr,"('freeze_all      = ',L1)") freeze_all
            call log_msg(msgstr)
 
       call log_ws()

        call MPI_Bcast(R_stokes,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(threshold_v,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(threshold_f,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(sigma,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(md_rock_fric_x,1,MPI_LOGICAL,0,comm_cart,ierror)
        call MPI_Bcast(md_rock_fric_y,1,MPI_LOGICAL,0,comm_cart,ierror)
        call MPI_Bcast(md_rock_fric_z,1,MPI_LOGICAL,0,comm_cart,ierror) 
        call MPI_Bcast(freeze_all,1,MPI_LOGICAL,0,comm_cart,ierror) 
        r_cut = 1.0_rk*sigma

        return
100     continue
        call error_md('Error reading md input file "'//trim(inp_file)//'.md"')
    end subroutine input_rock_friction

    subroutine rock_f_interaction_friction(N)
        type(lbe_site),intent(in) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
	     real(kind=rk) :: u(3) ! fluid velocity field at particle position
        real(kind=rk) :: mu ! fluid dynamic viscosity at particle position
        real(kind=rk) :: f(3)   ! friction force on md particle
	logical stuck
        integer i,ii,j,x,y,z
        real(kind=rk) :: dist(3),dist_sq, v_sq,f_sq
        integer n_c
        type(local_chunk_type) :: c(n_max_local_chunks)

        i = atompnt
        particles: do ii = 1,nlocal+nother
           call local_chunks(P(i),r_cut,0,n_c,c)
           chunks: do j=1,n_c
              x_loop: do x=c(j)%minx(1),c(j)%maxx(1)
                 y_loop: do y=c(j)%minx(2),c(j)%maxx(2)
                    z_loop: do z=c(j)%minx(3),c(j)%maxx(3)
                       surface: if (N(x,y,z)%rock_state==-2.0_rk) then
                          dist = c(j)%xc-real((/x,y,z/),kind=rk)
                          dist_sq = dot_product(dist,dist)
                          cutoff: if (dist_sq<r_cut_sq) then
	!	 		 call fluid_velocity_and_viscosity&
         !       &(N,(/1.0_rk,1.0_rk,1.0_rk/),P(i)%x,u,mu,stuck)
           ! ! abort on the first stuck particle
           !if (stuck) then
           !   print '(A,I10,A,3ES10.2,A,3ES10.2,A)'&
           !        &,'WARNING: lost a particle in rock: uid='&
           !        &,P(i)%uid,',x=(',P(i)%x(:),'),v=(',P(i)%v(:),')'
           !   call error_md('program exits...')
           !end if
				if (freeze_all) then 

						P(i)%freeze(1) = 1.0_rk
                 P(i)%freeze(2) = 1.0_rk
                 P(i)%freeze(3) = 1.0_rk
				else

           ! Stokes friction: F=6 pi R mu deltav
           f(:) = R_stokes*P(i)%v(:)
          !  print '(A,3ES10.2,A,3ES10.2,A,3ES10.2,A)'& 
           !   &,'force: fx='&
           !   &,f(1),',fy=(',f(2),'),fz',f(3),')'
           ! force on particle
           if(md_rock_fric_x .AND. md_rock_fric_y) then 
             v_sq = sqrt(P(i)%v(1)*P(i)%v(1) + P(i)%v(2)*P(i)%v(2))
             f_sq = sqrt(P(i)%f(1)*P(i)%f(1) + P(i)%f(2)*P(i)%f(2))
             if (v_sq < threshold_v .or. f_sq < threshold_f) then         
               P(i)%freeze(1) = 1.0_rk
               P(i)%freeze(2) = 1.0_rk
               P(i)%freeze(3) = 1.0_rk
       		 else
               P(i)%freeze(1) = 0.0_rk
               P(i)%freeze(2) = 0.0_rk
	       P(i)%f(1) = P(i)%f(1) - threshold_f * cos(atan2(P(i)%v(2),P(i)%v(1)))
	       P(i)%f(2) = P(i)%f(2) - threshold_f * sin(atan2(P(i)%v(2),P(i)%v(1)))
			    end if
		     end if
           if(md_rock_fric_x .AND. md_rock_fric_z) then 
             v_sq = sqrt(P(i)%v(1)*P(i)%v(1) + P(i)%v(3)*P(i)%v(3))
             f_sq = sqrt(P(i)%f(1)*P(i)%f(1) + P(i)%f(3)*P(i)%f(3))
             if (v_sq < threshold_v .or. f_sq < threshold_f) then         
               P(i)%freeze(1) = 1.0_rk
               P(i)%freeze(3) = 1.0_rk
               P(i)%freeze(2) = 1.0_rk
       		 else
               P(i)%freeze(1) = 0.0_rk
               P(i)%freeze(3) = 0.0_rk
	       P(i)%f(1) = P(i)%f(1) - threshold_f * cos(atan2(P(i)%v(3),P(i)%v(1)))
	       P(i)%f(3) = P(i)%f(3) - threshold_f * sin(atan2(P(i)%v(3),P(i)%v(1)))
			    end if
		     end if
           if(md_rock_fric_z .AND. md_rock_fric_y) then 
             v_sq = sqrt(P(i)%v(3)*P(i)%v(3) + P(i)%v(2)*P(i)%v(2))
             f_sq = sqrt(P(i)%f(3)*P(i)%f(3) + P(i)%f(2)*P(i)%f(2))
             if (v_sq < threshold_v .or. f_sq < threshold_f) then         
               P(i)%freeze(3) = 1.0_rk
               P(i)%freeze(2) = 1.0_rk
               P(i)%freeze(1) = 1.0_rk
       		 else
               P(i)%freeze(3) = 0.0_rk
               P(i)%freeze(2) = 0.0_rk
	       P(i)%f(3) = P(i)%f(3) - threshold_f * sin(atan2(P(i)%v(3),P(i)%v(2)))
	       P(i)%f(2) = P(i)%f(2) - threshold_f * cos(atan2(P(i)%v(3),P(i)%v(2)))
			    end if
		     end if


end if ! freeze_all
           !P(i)%f(1) = P(i)%f(1) - f(1)
           !if(md_rock_fric_y) P(i)%f(2) = P(i)%f(2) - f(2)
           !if(md_rock_fric_z) P(i)%f(3) = P(i)%f(3) - f(3)

                             if (lubrication)&
                                  & call lubrication_force_rock(dist_sq,dist,i)
                          end if cutoff
                       end if surface
                    end do z_loop
                 end do y_loop
              end do x_loop
           end do chunks
           if (ii<=nlocal) then
              i = list(i)
           else
              i = i+1
           endif
        end do particles
    end subroutine rock_f_interaction_friction

    subroutine setup_rock_friction
        ! save expensive instructions in the force loop
        r_cut_sq = r_cut*r_cut
        
        !6*pi*R_stokes
        !Rpi6 = 6.0_rk*pi*R_stokes
        ! this check was added due to the dependence on
        ! rc_lubrication_rock below and because
        ! lubrication_force_rock() is not compatible with
        ! polydispersity yet
        if (polydispersity) call error_md("rock='friction' does not yet "&
             &//"support polydisperse particles---disable polydispersity "&
             &//"or select rock/='friction'!")

        !if (lubrication.and.r_cut<rc_lubrication_rock) call error_md(&
        !     &'r_cut too small! '&
        !     &//'(r_cut<rc_lubrication_rock but lubrication==.true.')
    end subroutine setup_rock_friction


#endif
end module lbe_md_rock_friction_module
