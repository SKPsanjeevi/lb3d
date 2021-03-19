#include "lbe.h"

!> Growing stage for stable and efficient initialization of stiff
!> particles at high volume concentrations
!>
!> After whatever particle initialization, particles are shrunk and
!> then grow back to their original size while feeling only potential
!> interactions before, after an equilibration phase, the actual
!> simulation starts. This is intended for efficient and stable
!> particle initialization, especially at high volume fractions of
!> stiff particles. Also growing beyond the original size is possible
!> to achieve additional particle separation.
module lbe_md_growing_stage_module
#ifdef MD
    use lbe_globals_module, only: input_dfile_unit,maxpos,minpos,tsize,myrankc
    use lbe_init_rock_module, only: rock_is_present
    use lbe_log_module
    use lbe_md_globals_module
    use lbe_md_helper_module, only: error_md,log_msg_md,log_msg_md_hdr
    use lbe_md_module, only: borders,boundary_width,calculate_all_orientations&
         &,communicate,exchange,fix_shear_boundary,integrate_q,integrate_v&
         &,integrate_w,integrate_x,list_still_valid,neighbor&
         &,reset_forces_and_torques
    use lbe_md_output_module, only: dump_configuration
    use lbe_md_potential_module, only: force_and_torque&
         &,set_potential_growth_factor
    use lbe_parallel_module, only: comm_cart
    use lbe_parms_module, only: arg_input_dfile,arg_input_dfile_set,inp_file

    implicit none
    private
    include 'mpif.h'

    public growing_stage_run,input_growing_stage,setup_growing_stage

    !> switches growing_stage on/off
    logical,save,public :: growing_stage=.false.

    !> start growth at this fraction of particle lengths
    real(kind=rk),save :: initial_length_factor=0.7_rk

    !> grow particles up to this multiple of particle lengths.
    !>
    !> Values >1 allow to start the simulation with a minimum gap
    !> between particles. After the growth stage, lengths are always
    !> reset to a factor of 1.
    real(kind=rk),save :: final_length_factor=1.02_rk

    !> damping for translational and rotational velocities during
    !> growth time steps
    real(kind=rk),save :: damping=0.0001_rk

    !> time step size during growth stage
    real(kind=rk),save :: dt=1.0_rk

    !> particle growth is stretch over that many time steps
    real(kind=rk),save :: growth_time=100000.0_rk

    !> additional time to let particles equilibrate after \c growth_time
    real(kind=rk),save :: equilibration_time=500000.0_rk

    !> dump every this many time steps during growth and equilibration
    !> (0: never)
    integer,save :: n_dump=0

    namelist /md_growing_stage/ damping,dt,equilibration_time&
         &,final_length_factor,initial_length_factor,growth_time,n_dump

contains

    !> calculate the scaling factor for a given point in time during
    !> growth stage
    !>
    !> \param[in] t current time
    !>
    !> The resulting scaling factors lead to a linear growth of volume
    !> over time (Thanks to TK for this hint!) rather than a linear
    !> growth of (linear) particle size.
    real(kind=rk) function current_scaling_factor(t)
        real(kind=rk),intent(in) :: t
        real(kind=rk) :: ff3,fi3

        ff3 = final_length_factor**3
        fi3 = initial_length_factor**3

        current_scaling_factor = (fi3+t*(ff3-fi3)/growth_time)**(1.0_rk/3.0_rk)
    end function current_scaling_factor

    !> perform whole particle growth stage including possible
    !> equilibration and reset particle lengths afterwards to the
    !> values read from the input file
    subroutine growing_stage_run
        real(kind=rk) :: f,t
        integer :: ts
        character(len=64) :: prefix

        call log_msg_md('Entering particle growing stage...')

        call set_potential_growth_factor(initial_length_factor)

        call exchange
        call borders
        call neighbor

        write (msgstr&
             &,"(' starting growth phase at length factor        ',F16.10)") &
             &initial_length_factor
        call log_msg_md(msgstr)

        ! actual growth phase
        t = 0.0_rk
        ts = 0
        do while (t<growth_time)
           ts = ts + 1
           if (n_dump/=0) then
              write (prefix,fmt='("growth-",F6.3)') current_scaling_factor(t)
              if (mod(ts,n_dump)==0) call dump_configuration(prefix=prefix,t=ts)
           end if

           call growth_stage_time_step

           t = t+dt
           f = current_scaling_factor(t)
           call set_potential_growth_factor(f)
        end do
        write (msgstr&
             &,"(' starting equilibration phase at length factor ',F16.10)") f
        call log_msg_md(msgstr)

        ! equilibration phase
        t = 0.0_rk
        ts = 0
        do while (t<equilibration_time)
           ts = ts + 1
           if (n_dump/=0) then
              write (prefix,fmt='("equil-",F6.3)') f
              if (mod(ts,n_dump)==0) call dump_configuration(prefix=prefix,t=ts)
           end if

           call growth_stage_time_step
           t = t+dt
        end do

        write (msgstr&
             &,"(' resetting potential length scaling factor to  ',F16.10)") 1.0
        call log_msg_md(msgstr)

        call set_potential_growth_factor(1.0_rk)
    end subroutine growing_stage_run

    !> minimal time step containing only integration, minimal
    !> communication, potential evaluation and call to \c
    !> fix_shear_boundary()
    subroutine growth_stage_time_step
        integer :: i,ii,k
        real(kind=rk),parameter :: small=1.0E-9_rk

        i = atompnt
        do ii = 1,nlocal
           if (use_rotation) call integrate_q(P(i),0.5_rk*dt)
#ifdef RWALK
           call error_md('Growing stage not implemented for RWALK---'&
                &//'set growing_stage=.false. or compile without RWALK!')
#else
           call integrate_x(P(i),dt)
#endif

           ! impose minimal (periodic) boundary conditions here
           do k=1,3
              if (P(i)%x(k)>=maxpos(k)) P(i)%x(k) = P(i)%x(k)-tsize(k)
              if (P(i)%x(k)<minpos(k)) P(i)%x(k) = P(i)%x(k)+tsize(k)
              if (P(i)%x(k)==maxpos(k)) P(i)%x(k) = P(i)%x(k)-small
           end do

           i = list(i)
        end do

        ! minimal communication, don't provide uid2i or rotation_s,
        ! they are not required by simple pair potentials
        if (list_still_valid()) then
           call communicate
        else
           call exchange
           call borders
           call neighbor
        end if

        if (calculate_orientations) call calculate_all_orientations

        ! the 1 pretends this would be the first (possibly) only
        ! substep of an LB step, it does not matter here
        call reset_forces_and_torques(1)

        call force_and_torque

        i = atompnt
        do ii = 1,nlocal
           call integrate_v(P(i),dt)
           if (use_rotation) call integrate_w(P(i),dt)

           P(i)%vnew = P(i)%vnew*(1.0_rk-damping)
           P(i)%wnew = P(i)%wnew*(1.0_rk-damping)

           i = list(i)
        enddo

        if (boundary_width/=0.0_rk) call fix_shear_boundary()
    end subroutine growth_stage_time_step

    !> read  \c /md_growing_stage/  from the md input file
    subroutine input_growing_stage
        integer ierror

        call log_msg_md_hdr("Reading MD growing stage input")

        if (myrankc==0) then
           open (unit=md_input_file_unit,file=trim(inp_file)//'.md',err=100)
           read (unit=md_input_file_unit,nml=md_growing_stage,err=100)
           close (unit=md_input_file_unit,err=100)
        end if

        if ( arg_input_dfile_set ) then
          call log_msg_md("  Getting differential input...")
          open (unit=input_dfile_unit,file=arg_input_dfile,status='UNKNOWN')
          read (unit=input_dfile_unit,nml=md_growing_stage,iostat=ierror)
          if (ierror/=0) call log_msg_md('    WARNING: Differential namelist '&
               &//'not found or errors encountered.')
          close (unit=input_dfile_unit)
          call log_ws()
        end if

        write (msgstr,"('initial_length_factor = ',F16.6)") &
             &initial_length_factor
        call log_msg(msgstr)
        write (msgstr,"('final_length_factor   = ',F16.6)") final_length_factor
        call log_msg(msgstr)
        write (msgstr,"('damping               = ',F16.6)") damping
        call log_msg(msgstr)
        write (msgstr,"('dt                    = ',F16.6)") dt
        call log_msg(msgstr)
        write (msgstr,"('growth_time           = ',F16.6)") growth_time
        call log_msg(msgstr)
        write (msgstr,"('equilibration_time    = ',F16.6)") equilibration_time
        call log_msg(msgstr)
        write (msgstr,"('n_dump                = ',I0)") n_dump
        call log_msg(msgstr)
        call log_ws()

        call MPI_Bcast(initial_length_factor,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(final_length_factor,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(damping,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(dt,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(growth_time,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(equilibration_time,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(n_dump,1,MPI_INTEGER,0,comm_cart,ierror)

        return
100     continue
        call error_md('Error reading md input file "'//trim(inp_file)//'.md"')
    end subroutine input_growing_stage

    !> setup stage for \c lbe_md_growing_stage_module
    subroutine setup_growing_stage
        if (rock_is_present()) call error_md('md_growing_stage does not '&
             &//'support rock geometries yet---set growing_stage=.false. or '&
             &//'disable the geometry (obs_file=''empty.dat'' and/or change '&
             &//'boundary_cond)')
    end subroutine setup_growing_stage
#endif
end module lbe_md_growing_stage_module
