#include "lbe.h"

!> ought to provide the only md subroutines that are called directly
!> from the lb part of the code
module lbe_md_interface_module
#ifdef MD
    use lbe_bc_module, only: lbe_bc_complete_halo_exchange
    use lbe_globals_module, only: bounce,chunksize,cz,g,halo_extent,myrankc &
        &,forcesum,totalforce
#ifdef RWALK
    use lbe_helper_module, only: request_halo_extent
#endif
    use lbe_log_module
    use lbe_md_bc_leesedwards_module, only: input_md_leesedwards&
         &,setup_md_leesedwards
    use lbe_md_dynamic_module, only: set_dissolution_limits_z,setup_dynamic
    use lbe_md_fluid_module, only: fluid_force_halo_extent,fluid_halo_extent&
         &,fluid_post,fluid_middle,init_fluid,input_fluid,setup_fluid&
         &,shutdown_fluid
    use lbe_md_globals_module
    use lbe_md_growing_stage_module, only: input_growing_stage,growing_stage&
         &,growing_stage_run,setup_growing_stage
    use lbe_md_helper_module, only: error_md,log_msg_md,log_msg_md_hdr&
         &,orientation
    use lbe_md_init_module, only: setup_particles,setup_general,setup_parallel&
         &,setup_uid2i
    use lbe_md_input_module, only: input,md_restore_checkpoint
    use lbe_md_lyapunov_module, only: input_lyapunov,lyapunov,lyapunov_run&
         &,setup_lyapunov,shutdown_lyapunov
    use lbe_md_memory_module, only: setup_memory,shutdown_memory
    use lbe_md_module, only: borders,boundary_width,calculate_all_orientations&
         &,calculate_rotations_s,check_shear_boundary,create_uid2i,exchange&
         &,neighbor,integrate_vw,integrate_xq,statistics,setup_computations&
         &,status
    use lbe_md_output_module, only: dump_configuration,dump_rvz_profile&
         &,input_dump,n_dump,n_dump_substep,n_dump_substep_start,dump_summary&
	 &,n_msdx,n_msdx_sample,sample_msdx&
         &,n_msdy,n_msdy_sample,sample_msdy&
	 &,n_msdz,n_msdz_sample,sample_msdz&
	 &,n_rvz_profile_sample&
         &,setup_dump,summary
    use lbe_md_potential_module, only: input_potential,setup_potential
    use lbe_md_rock_module, only: input_rock,rock_halo_extent,setup_rock
    use lbe_timer_module,only: start_timer,stop_timer
    use lbe_parallel_module, only: comm_cart,start,tnx
    use lbe_parms_module, only: n_checkpoint,n_iteration,nt,n_restore,nx,ny,nz,inv_fluid
    use lbe_types_module, only: lbe_site
    use lbe_helper_module, only: check_dump_now,every_n_time_steps,is_restoring
    use lbe_md_magnetic_module, only: md_magnetic, input_magnetic
#ifdef DEBUG_MPI
    use lbe_parallel_module, only: log_mpi
#endif
    use lbe_md_fluid_ladd_mc_module, only: pfr, pfb, pfg, n_renew_rho, set_initial_average_density
    implicit none
    private

    public md_force_halo_extent,md_halo_extent,md_init,md_input&
         &,md_lbe_mass_scaler_apply,md_run_post,md_run_middle&
         &,md_run_pre_after_bc,md_run_pre_before_bc,md_shift_velocities&
         &,md_shutdown,md_calculate_dfz

    include 'mpif.h'

contains

    !> initializes  halo_extent  due to the chosen simulation parameters
    subroutine md_halo_extent
        call fluid_halo_extent
        call rock_halo_extent
#ifdef RWALK
        call request_halo_extent(4,"RWALK")
#endif
    end subroutine md_halo_extent

    !> initializes \c force_halo_extent due to the chosen simulation parameters
    subroutine md_force_halo_extent
        call fluid_force_halo_extent
    end subroutine md_force_halo_extent

    !> has to be run AFTER mpi cartesian topology is set up!
    !>
    !> \param[in,out] lbe_N local chunk of the lattice with halo
    !> extent 1 (old LB3D style)
    !>
    !> \param[in,out] N local lattice chunk including full halo
    !>
    !> \param[in,out] d_adv at the moment just passed to \c
    !> lbe_bc_complete_halo_exchange()
    subroutine md_init(lbe_N,N,d_adv)
        type(lbe_site),intent(inout) :: lbe_N(0:,0:,0:)
        type(lbe_site),intent(inout) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
#ifndef NOSURFACTANT
        real(kind=rk),dimension(1:,1:,0:,0:,0:),intent(inout) :: d_adv
#else
        integer,intent(inout) :: d_adv
#endif

        call log_msg_md_hdr('Initializing system')
        call setup_general
        call log_msg_md('Finished setup_general.')
        call setup_computations
        call log_msg_md('Finished setup_computations.')
        call setup_memory
        call log_msg_md('Finished setup_memory.')
        call setup_parallel
        call log_msg_md('Finished setup_parallel.')
        call setup_md_leesedwards(N)
        call log_msg_md('Finished setup_md_leesedwards.')
        call setup_potential
        call log_msg_md('Finished setup_potential.')
        call setup_fluid(N)
        call log_msg_md('Finished setup_fluid.')
        call setup_rock(N)
        call log_msg_md('Finished setup_rock.')
        call setup_dump
        call log_msg_md('Finished setup_dump.')
        call md_setup_dfz
        call log_msg_md('Finished md_setup_dfz.')
        if (growing_stage) then
           call setup_growing_stage
           call log_msg_md('Finished setup_growing_stage.')
        end if
        if (lyapunov) then
           call setup_lyapunov(N)
           call log_msg_md('Finished setup_lyapunov.')
        end if

        ! Initialize additional rock_state halo if requested by one of
        ! the fluid interaction modules. (As long as rock_state does
        ! not change during the simulation, no further communication
        ! of the md halo is necessary.)
        if (halo_extent>1) then
           ! Be aware that in case of Lees-Edwards boundary conditions
           ! only the 1st halo layer will be communicated but in this
           ! case there probably should be no walls anyway.
           call lbe_bc_complete_halo_exchange(lbe_N,N,d_adv)
           call log_msg_md('Initialized additional rock_state halo.')
        end if

        if ( is_restoring() ) then
           if (n_iteration<n_restore) call error_md(&
                &'n_iteration is smaller than the timestep to restore from')
           call md_restore_checkpoint()
        else
           call setup_particles(N)
           if (growing_stage) call growing_stage_run
        end if
        call log_msg_md('Created particles.')

        call setup_dynamic

        if (provide_uid2i) call setup_uid2i
        call exchange
        if (communicate_rotations_s) call calculate_rotations_s
        call borders
        if (provide_uid2i) call create_uid2i
        call neighbor
        call log_msg_md('Initialized neighbor lists.')

        ! calculate orientations for own and haloed particles - they
        ! might be needed for status calculations and fluid
        ! initialization
        if (calculate_orientations) call calculate_all_orientations

        call init_fluid(N)
        call log_msg_md('Initialized fluid.')

        call log_msg_md('Initialization complete.')

        if (n_dump/=0) call dump_configuration
        if (n_rvz_profile_sample/=0) call dump_rvz_profile(N)

        if (n_stat/=0) call status(N)

        ! don't sample when restoring, otherwise the checkpointed
        ! configuration would be sampled twice
        if (n_msdx/=0.and.every_n_time_steps(n_msdx_sample)&
             &.and..not.is_restoring()) call sample_msdx()
        if (n_msdy/=0.and.every_n_time_steps(n_msdy_sample)&
             &.and..not.is_restoring()) call sample_msdy()
	if (n_msdz/=0.and.every_n_time_steps(n_msdz_sample)&
	     &.and..not.is_restoring()) call sample_msdz()

        ! check only newly created configurations
        if (.not.is_restoring().and.boundary_width/=0.0_rk) &
             &call check_shear_boundary()

        md_initialized = .true.
    end subroutine md_init

    !> read input file, must be called before  md_halo_extent() .
    subroutine md_input
        call input
        call input_potential
        call input_fluid
        call input_rock
        call input_dump
        call input_md_leesedwards()
        if (md_magnetic) call input_magnetic
        if (growing_stage) call input_growing_stage
        if (lyapunov) call input_lyapunov
    end subroutine md_input

    !> rescales whatever MD quantities should be rescaled in sync with
    !> \c lbe_mass_scaler_module
    !>
    !> \param[in] factor mass rescaling factor as obtained from \c
    !> mass_scaler_apply()
    subroutine md_lbe_mass_scaler_apply(factor)
        real(kind=rk),intent(in) :: factor
        integer i,ii

#ifdef LADD_SURR_RHOF
        ! rescale the surrounding fluid densities. As they are in
        ! direct interplay with N%n_r, this seems much more important
        ! than worrying about the particle masses. To be on the safe
        ! side, rescale both owned and halo'ed particles and both the
        ! accumulation buffer and the rhof itself as it is not clear
        ! yet where exactly in the time step this subroutine call will
        ! end up when things are shuffled around...
        i = atompnt
        do ii=1,nlocal+nother
           P(i)%rhof = factor*P(i)%rhof
           P(i)%rhof_acc = factor*P(i)%rhof_acc

           if (ii<=nlocal) then
              i = list(i)
           else
              i = i+1
           end if
        end do
#endif
    end subroutine md_lbe_mass_scaler_apply

    !> called between lbe streaming and collision
    !>
    !> \param[in,out] lbe_N local chunk of the lattice with halo
    !> extent 1 (old LB3D style)
    !>
    !> \param[in,out] whole_N local chunk of the lattice with full
    !> halo of depth \c halo_extent
    !>
    !> \param[in,out] d_adv at the moment just passed to \c
    !> lbe_bc_complete_halo_exchange()
    subroutine md_run_middle(lbe_N,whole_N,d_adv)
        type(lbe_site),intent(inout) :: lbe_N(0:,0:,0:)
        type(lbe_site),intent(inout) :: &
             &whole_N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
#ifndef NOSURFACTANT
        real(kind=rk),dimension(1:,1:,0:,0:,0:),intent(inout) :: d_adv
#else
        integer,intent(inout) :: d_adv
#endif

        call start_timer(ti_md_fluid)
        DEBUG_MPI_MSG("Entering fluid_middle(...)")
        call fluid_middle(whole_N)
        DEBUG_MPI_MSG("Returned from fluid_middle(...)")
        call stop_timer(ti_md_fluid)

        ! Another halo exchange if requested by MD interaction
        if (halo_exchange_after_md_run_middle) &
             &call lbe_bc_complete_halo_exchange(lbe_N,whole_N,d_adv)
    end subroutine md_run_middle

    !> called after lbe streaming and collision
    subroutine md_run_post(N)
        type(lbe_site),intent(inout) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        
        integer i,ierror

        forcesum(:) = 0.0_rk
        totalforce(:) = 0.0_rk

        !for evaporation model, the inital density for mass correction is
        !calculated in every n_recalculate_rho timestep
          if (inv_fluid == 26) then 
             if (every_n_time_steps(n_renew_rho))   call set_initial_average_density(N)
          end if
        ! perform more than one integration step? The first
        !  integrate_xq()  is always called in  md_run_pre() .
        steps: do i=1,steps_per_lbe_step - 1
           call integrate_vw(N,i)
           if (n_dump_substep/=0.and.nt>=n_dump_substep_start) then
              if (mod(i,n_dump_substep)==0) call dump_configuration(i)
           end if

           call md_start_substep(i+1)
           call integrate_xq(N,i+1)
        end do steps

        DEBUG_MPI_MSG("Entering integrate_vw(...)")
        call integrate_vw(N,steps_per_lbe_step)
        DEBUG_MPI_MSG("Returned from integrate_vw(...)")

        if (lyapunov) call lyapunov_run(N)

        if (n_msdx/=0.and.every_n_time_steps(n_msdx_sample)) call sample_msdx()
        if (n_msdy/=0.and.every_n_time_steps(n_msdy_sample)) call sample_msdy()
	if (n_msdz/=0.and.every_n_time_steps(n_msdz_sample)) call sample_msdz()

        if (n_dump/=0) then
           if (mod(nt,n_dump)==0) then
              call start_timer(ti_md_dump)
              call dump_configuration
              call stop_timer(ti_md_dump)
           end if
        end if

        if (n_rvz_profile_sample/=0) then
           if (mod(nt,n_rvz_profile_sample)==0) then
              call start_timer(ti_md_dump)
              call dump_rvz_profile(N)
              call stop_timer(ti_md_dump)
           end if
        end if

        if (n_stat/=0) then
           if (mod(nt,n_stat)==0) call status(N)
        end if

        call start_timer(ti_md_fluid)
        DEBUG_MPI_MSG("Entering fluid_post(...)")
        call fluid_post(N)
        DEBUG_MPI_MSG("Returned from fluid_post(...)")
        call stop_timer(ti_md_fluid)
    
        call mpi_allreduce(forcesum,totalforce,3,MPI_REAL8,MPI_SUM,comm_cart,ierror)
                
    end subroutine md_run_post

    !> needs to be called before lbe streaming and collision
    !>
    !> \param[in] N local lattice chunk including full halo
    subroutine md_run_pre(N)
        type(lbe_site),intent(inout) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)

        ! announce start of first (and possible only) substep
        call md_start_substep(1)

        if (every_n_time_steps(n_dfz)) call md_reset_dfz

        DEBUG_MPI_MSG("Entering integrate_xq(...)")
        call integrate_xq(N,1)
        DEBUG_MPI_MSG("Returned from integrate_xq(...)")
    end subroutine md_run_pre

    !> is called between lbe streaming and the halo exchange/boundary
    !> condition update before
    !>
    !> \param[in] N local lattice chunk including full halo
    subroutine md_run_pre_after_bc(N)
        type(lbe_site),intent(inout) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)

        ! only if specially requested
        if (halo_exchange_before_md_run) call md_run_pre(N)

        if (every_n_time_steps(n_dfz)) call md_calculate_dfz(N)
    end subroutine md_run_pre_after_bc

    !> is called before the halo exchange/boundary condition update
    !> which comes before lbe streaming
    !>
    !> \param[in] N local lattice chunk including full halo
    subroutine md_run_pre_before_bc(N)
        type(lbe_site),intent(inout) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)

        ! this is the default
        if (.not.halo_exchange_before_md_run) call md_run_pre(N)
    end subroutine md_run_pre_before_bc

    !> shift all velocities by a constant offset
    !>
    !> The routine assumes to be called after the last integrate_vw()
    !> or before the first integrate_xq() of an LB step.
    !>
    !> \param[in] dv velocity offset
    subroutine md_shift_velocities(dv)
        real(kind=rk),intent(in) :: dv(3)
        integer :: i,ii

        i = atompnt
        do ii=1,nlocal+nother
           ! since this is not performance-critical shift all velocity
           ! elements of md_particle_type for both owned and halo'd
           ! particles to be on the safe side though some of them
           ! might be overwritten or not even used.
           P(i)%v = P(i)%v + dv
           P(i)%vnew = P(i)%vnew + dv
           P(i)%v_fluid = P(i)%v_fluid + dv
           P(i)%v_fluid_avg = P(i)%v_fluid_avg + dv

           ! This assumes that there were steps_per_lbe_step-1
           ! velocities accumulated in v_fluid_acc which is the case
           ! at the end of an LB step.
           P(i)%v_fluid_acc = P(i)%v_fluid_acc + &
                &real(steps_per_lbe_step-1,kind=rk)*dv

#ifdef RWALK
           P(i)%v_r = P(i)%v_r + dv
#endif

           if (ii<=nlocal) then
              i = list(i)
           else
              i = i+1
           end if
        end do
    end subroutine md_shift_velocities

    subroutine md_shutdown
        if (lyapunov) call shutdown_lyapunov
        call shutdown_fluid
        call statistics
        if ( dump_summary ) then
          call summary
        end if
        call shutdown_memory
    end subroutine md_shutdown

    !> announces the start of a new substep
    !>
    !> \param[in] s new substep index (\c 1..steps_per_lbe_step)
    subroutine md_start_substep(s)
        integer,intent(in) :: s

        nt_substep = s
    end subroutine md_start_substep

    subroutine md_calculate_dfz(N)
        type(lbe_site),intent(in) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        real(kind=rk) :: dfz(2),dfz_f_sum(2),dfz_p_sum(2),dfz_fp_sum(2)&
             &,dfz_pp_sum(2)
        integer ierror

        !  dfz_p ,  dfz_fp , and dfz_pp  were alreade summed up in force
        ! calculation, integration, and fluid coupling routines
        call md_local_fluid_momentum_transfer(N,dfz_f)

        ! The momentum transfer through each of both planes is caused
        ! by forces from both sides. So the values for both planes are
        ! multiplied by one half here in order to let their sum be the
        ! force from both sides times one timestep. The values in
        !  dfz_fp  and  dfz_pp  are already only the half of the
        ! momentum transferred across each plane, so they don't have
        ! to be halved here.
        dfz_f(:) = dfz_f(:)*0.5_rk
        dfz_p(:) = dfz_p(:)*0.5_rk

        call mpi_reduce(dfz_f,dfz_f_sum,2,MPI_REAL8,MPI_SUM,0,comm_cart,ierror)
        call mpi_reduce(dfz_p,dfz_p_sum,2,MPI_REAL8,MPI_SUM,0,comm_cart,ierror)
        call mpi_reduce(dfz_fp,dfz_fp_sum,2,MPI_REAL8,MPI_SUM,0,comm_cart&
             &,ierror)
        call mpi_reduce(dfz_pp,dfz_pp_sum,2,MPI_REAL8,MPI_SUM,0,comm_cart&
             &,ierror)

        if (myrankc==0) then
           dfz = dfz_f_sum+dfz_p_sum+dfz_fp_sum+dfz_pp_sum
           write (unit=6,fmt='(SP,2(A,ES15.8))') &
                &'dfz(1)   =',dfz(1),' dfz(2)   =',dfz(2)

           write (unit=6,fmt='(SP,2(A,ES15.8))') &
                &'dfz_f(1) =',dfz_f_sum(1),' dfz_f(2) =',dfz_f_sum(2)
           write (unit=6,fmt='(SP,2(A,ES15.8))') &
                &'dfz_p(1) =',dfz_p_sum(1),' dfz_p(2) =',dfz_p_sum(2)
           write (unit=6,fmt='(SP,2(A,ES15.8))') &
                &'dfz_fp(1)=',dfz_fp_sum(1),' dfz_fp(2)=',dfz_fp_sum(2)
           write (unit=6,fmt='(SP,2(A,ES15.8))') &
                &'dfz_pp(1)=',dfz_pp_sum(1),' dfz_pp(2)=',dfz_pp_sum(2)

           ! this is convenient for shear experiments
           write (unit=6,fmt='(A,SS,I9,A,SP,ES15.8)') &
                &'nt=',nt,' dfz_shear=',-dfz(1)+dfz(2)
        end if
    end subroutine md_calculate_dfz

    !> Adds to \c lfmt(1:2) the momentum transfer through the layers
    !> \c x==dfz_minx and \c x==dfz_maxx (as far as part of the local
    !> chunk) caused by the fluid itself.
    !>
    !> \todo This does neigher work for
    !> \c amass_(r|b|s)/=1.0 nor without \c SINGLEFLUID!
    subroutine md_local_fluid_momentum_transfer(N,lfmt)
        type(lbe_site),intent(in) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        real(kind=rk),intent(inout) :: lfmt(2)
        ! vectors that are advected in x direction AND carry momentum in z
        ! direction (p:positive x-direction, m:negative x-direction)
        integer,parameter,dimension(2) :: pvecs=(/9,10/),mvecs=(/13,14/)
        integer local_dfz_minx,local_dfz_maxx,r,s,t,y,z

        local_dfz_minx = dfz_minx+1-start(1)
        if (1<=local_dfz_minx.and.local_dfz_minx<=nx) then
           ! loop through lattice nodes on plane x=dfz_minx
           do y=1,ny
              do z=1,nz
                 do t=1,size(pvecs)
                    s = pvecs(t)
                    r = bounce(s)
                    ! both source and target node have to be fluid nodes
                    if (N(local_dfz_minx-1,y,z-cz(s))%rock_state==0.0_rk&
                         &.and.N(local_dfz_minx,y,z)%rock_state==0.0_rk) then
                       ! sum incoming and outgoing momentum
                       lfmt(1) = lfmt(1)&
                            &+N(local_dfz_minx-1,y,z-cz(s))%n_r(s)*cz(s)*g(s)&
                            &-N(local_dfz_minx,y,z)%n_r(r)*cz(r)*g(r)
                    end if
                 end do
              end do
           end do
        end if

        local_dfz_maxx = dfz_maxx+1-start(1)
        if (1<=local_dfz_maxx.and.local_dfz_maxx<=nx) then
           ! loop through lattice nodes on plane x=dfz_maxx
           do y=1,ny
              do z=1,nz
                 do t=1,size(pvecs)
                    s = mvecs(t)
                    r = bounce(s)
                    ! both source and target node have to be fluid nodes
                    if (N(local_dfz_maxx+1,y,z-cz(s))%rock_state==0.0_rk&
                         &.and.N(local_dfz_maxx,y,z)%rock_state==0.0_rk) then
                       ! sum incoming and outgoing momentum
                       lfmt(2) = lfmt(2)&
                            &+N(local_dfz_maxx+1,y,z-cz(s))%n_r(s)*cz(s)*g(s)&
                            &-N(local_dfz_maxx,y,z)%n_r(r)*cz(r)*g(r)
                    end if
                 end do
              end do
           end do
        end if
    end subroutine md_local_fluid_momentum_transfer

    !> reset dfz variables to zero
    subroutine md_reset_dfz
        dfz_f(:) = 0.0_rk
        dfz_p(:) = 0.0_rk
        dfz_fp(:) = 0.0_rk
        dfz_pp(:) = 0.0_rk
    end subroutine md_reset_dfz

    !> setup dfz boundary coordinates
    subroutine md_setup_dfz
        if (dfz_minx<0) dfz_minx = 1
        if (dfz_maxx<0) dfz_maxx = tnx
    end subroutine md_setup_dfz

#endif
end module lbe_md_interface_module
