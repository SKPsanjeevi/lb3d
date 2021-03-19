#include "lbe.h"

!> ought to provide the only \c TRACER subroutines that are called
!> directly from \c lbe.F90
module lbe_tracer_interface_module
#ifdef TRACER

    use lbe_globals_module, only: halo_extent
    use lbe_helper_module, only: check_dump_now,every_n_time_steps,is_restoring,request_halo_extent
    use lbe_log_module
    use lbe_parms_module, only: n_checkpoint,n_iteration,n_restore,nt
    use lbe_timer_module,only: start_timer,stop_timer
    use lbe_tracer_globals_module
    use lbe_tracer_helper_module, only: error_tracer,log_msg_tracer&
         &, log_msg_tracer_hdr
    use lbe_tracer_init_module, only: place_tracers,setup_general&
         &,setup_parallel,t_placement
    use lbe_tracer_input_module, only: tracer_restore_checkpoint
    use lbe_tracer_memory_module, only: setup_memory
    use lbe_tracer_module, only: tracer_time_step
    use lbe_tracer_output_module, only: dump_configuration,dump_tracer_profile&
         &,input_dump,n_dump,n_flow_rate&
	 &,n_msdx,n_msdx_sample,sample_msdx&
	 &,n_msdy,n_msdy_sample,sample_msdy&
	 &,n_msdz,n_msdz_sample,sample_msdz&
         &,n_profile,report_tracer_flow_rate&
         &,setup_dump,tracer_dump_checkpoint
    use lbe_tracer_input_module, only: read_tracer_input
    use lbe_types_module, only: lbe_site
#ifdef MD
    use lbe_md_globals_module, only: polydispersity
#endif

    implicit none
    private

    public tracer_halo_extent,tracer_init,tracer_input,tracer_run&
         &,tracer_shift_velocities,tracer_shutdown

contains

    !> requests  halo_extent  to match the requirements of \c TRACER
    subroutine tracer_halo_extent
      call request_halo_extent(1,"TRACER")
    end subroutine tracer_halo_extent

    !> initialization of \c TRACER code and the tracer particles
    !> themselves according to their initial condition
    !>
    !> This has to be run AFTER mpi cartesian topology is set up. \c
    !> MD \c interaction=='ladd' particles should be placed on the
    !> fluid as well, if present.
    !>
    !> \param[in] N local lattice chunk including full halo
    subroutine tracer_init(N)
        type(lbe_site),intent(in) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        integer i,ii

#ifndef SINGLEFLUID
        call error_tracer('TRACER not implemented yet for multi-component case'&
             &//'---compile with SINGLEFLUID!')
#endif
#ifdef MD
        if (polydispersity) call error_tracer('TRACER not implemented yet for '&
             &//'polydisperse MD particles---compile without MD or disable MD '&
             &//'polydispersity!')
#endif

        call log_msg_tracer_hdr('Initializing system')
        call setup_general
        call log_msg_tracer('Finished setup_general.')
        call setup_memory
        call log_msg_tracer('Finished setup_memory.')
        call setup_parallel
        call log_msg_tracer('Finished setup_parallel.')
        call setup_dump
        call log_msg_tracer('Finished setup_dump.')

        if ( is_restoring() ) then
           if (n_iteration<n_restore) call error_tracer(&
                &'n_iteration is smaller than the timestep to restore from')
           call tracer_restore_checkpoint()
        else
           if (t_placement==0) then
              call place_tracers(N)
              call log_msg_tracer('Created tracers.')
           else
              write (msgstr&
                   &,"('Creation of tracers postponed to time step ',I0,'.')") &
                   &t_placement
              call log_msg_tracer(msgstr)
           end if
        end if

        call log_msg_tracer('Initialization complete.')

        if (every_n_time_steps(n_dump)) call dump_configuration()
        if (every_n_time_steps(n_profile)) call dump_tracer_profile(N)

        ! don't sample when restoring, otherwise the checkpointed
        ! configuration would be sampled twice
        if (n_msdx/=0.and.every_n_time_steps(n_msdx_sample)&
             &.and..not.is_restoring()) call sample_msdx()
        if (n_msdy/=0.and.every_n_time_steps(n_msdy_sample)&
             &.and..not.is_restoring()) call sample_msdy()
	if (n_msdz/=0.and.every_n_time_steps(n_msdz_sample)&
	     &.and..not.is_restoring()) call sample_msdz()
    end subroutine tracer_init

    !> read input file, must be called before \c tracer_halo_extent() .
    subroutine tracer_input()
        call read_tracer_input()

        call input_dump()
    end subroutine tracer_input

    !> runs everything necessary for \c TRACER during one LB time step
    !>
    !> \param[in] N local lattice chunk including full halo
    !>
    !> This routine relies on an updated halo and according \c MD \c
    !> interaction=='ladd' particle postions (if present).
    subroutine tracer_run(N)
        type(lbe_site),intent(in) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)

        ! possible postponed placement of tracers
        if (nt==t_placement) call place_tracers(N)

        call tracer_time_step(N)

        if (n_msdx/=0.and.every_n_time_steps(n_msdx_sample)) call sample_msdx()
        if (n_msdy/=0.and.every_n_time_steps(n_msdy_sample)) call sample_msdy()
	if (n_msdz/=0.and.every_n_time_steps(n_msdz_sample)) call sample_msdz()

        ! when dumping now, the dumped tracer velocities are actually
        ! the velocities at the previous position. This errors seems
        ! acceptable, the alternative would be to split
        ! tracer_time_step() or extend tracer_particle_type with
        ! additional data elements
        if (every_n_time_steps(n_flow_rate)) then
           call start_timer(ti_tracer_dump)
           call report_tracer_flow_rate()
           call stop_timer(ti_tracer_dump)
        end if
        if (every_n_time_steps(n_dump)) then
           call start_timer(ti_tracer_dump)
           call dump_configuration()
           call stop_timer(ti_tracer_dump)
        end if
        if (every_n_time_steps(n_profile)) then
           call start_timer(ti_tracer_dump)
           call dump_tracer_profile(N)
           call stop_timer(ti_tracer_dump)
        end if
    end subroutine tracer_run

    !> shift all velocities by a constant offset
    !>
    !> \param[in] dv velocity offset
    subroutine tracer_shift_velocities(dv)
        real(kind=rk),intent(in) :: dv(3)
        integer :: i,ii

        i = atompnt
        do ii = 1,nlocal
           T(i)%v = T(i)%v + dv

           i = list(i)
        end do
    end subroutine tracer_shift_velocities

    !> clean up or final statistics for \c TRACER could be placed here
    subroutine tracer_shutdown()
    end subroutine tracer_shutdown

#endif
end module lbe_tracer_interface_module
