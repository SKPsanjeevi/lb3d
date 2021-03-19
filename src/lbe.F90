
#ifndef LBE_RS

#include "lbe.h"

!> Reimplementation of LBE codes.
!> July 1999 Jonathan Chin
!>
!>This file contains the main program itself, but does little other than
!>call other routines at appropriate points. It serves to link the
!>whole program together.
program lbe
  use lbe_bc_module, only: lbe_bc_after_advection,lbe_bc_before_advection&
       &,lbe_bc_complete_halo_exchange,lbe_bc_halo_exchange,lbe_bc_init&
       &,lbe_bc_init_step,lbe_bc_input,lbe_bc_shutdown
  use lbe_collision_module
  use lbe_force_module, only: lbe_add_force_halo,lbe_force_apply,lbe_force_init&
       &,lbe_force_input,lbe_force_reset
  use lbe_galilean_stabilizer_module, only: lbe_galilean_stabilizer_init&
       &,lbe_galilean_stabilizer_input,lbe_galilean_stabilizer_run
  use lbe_globals_module
  use lbe_helper_module, only: report_check_NaN_function,unixtime,is_restoring, every_n_time_steps, set_halo_extent, request_halo_extent
  use lbe_init_module
  use lbe_io_interface_module, only: dump_parallel,lbe_io_halo_extent&
       &,lbe_io_init
  use lbe_io_module
  use lbe_io_checkpoint_module
  use lbe_io_stress_module, only: dump_stress
  use lbe_log_module
  use lbe_mass_scaler_module, only: lbe_mass_scaler_init,lbe_mass_scaler_input&
       &,lbe_mass_scaler_run
  use lbe_mean_square_displacement_mod, only: msd_init,msd_run
  use lbe_parms_module
  use lbe_parallel_module
  use lbe_timer_module, only: start_timer,stop_timer,sum_timer
  use lbe_types_module, only: lbe_site
#ifdef ELEC
  use lbe_elec_interface_module
#endif
#ifdef MD
  use lbe_md_interface_module, only: md_force_halo_extent, md_halo_extent, &
       &md_init, md_input, md_run_post, md_run_middle, md_run_pre_after_bc, &
       &md_run_pre_before_bc, md_shutdown
#endif
#ifdef SINGLEFLUID
  use lbe_collision_simple_module
#endif
#ifdef IBM_PART
  use lsuperstruct_module
#endif
#ifdef TRACER
  use lbe_tracer_interface_module, only: tracer_halo_extent,tracer_init&
       &,tracer_input,tracer_run,tracer_shutdown
#endif
#ifdef USEHDF
  use lbe_io_hdf5_module, only: lbe_io_init_hdf5, lbe_io_shutdown_hdf5,&
       lbe_add_inputfile_metadata_hdf5
#endif
#ifdef AXISYM
  use lbe_axisym_module
#endif
  use lbe_invasion_module, only: lbe_invade_after_evaporation
  implicit none

  integer :: ierror
  integer :: i, k
  integer :: tstart
  logical :: insanity

  ! This might be a queer hack but allows to widen the halo region just
  ! for md without need to touch anything outside the head of this file.

  ! Now  whole_N  is the whole lattice...
  type(lbe_site),dimension(:,:,:),allocatable,target :: whole_N
  ! ... N  is just a pointer to a part of whole_N .
  type(lbe_site),pointer :: N(:,:,:)


#ifndef NOSURFACTANT
  real(kind=rk), allocatable, dimension(:,:,:,:,:) :: d_adv
#else
  integer :: d_adv = 0
#endif

#ifdef AXISYM
  real*8, dimension(:,:,:,:)  ,allocatable :: rho
  real*8, dimension(:,:,:,:,:),allocatable :: u
  !(x/y/z, x,y,z, r/b/s)!  
#endif


  ! Initialise a few variables.
  insanity = .false.

  do i = 1, nvecs
    c(i,:) = (/ cx(i), cy(i), cz(i) /)
  end do

  ! Calculate relative coordinates for lattice points
  ! surrounding an arbitrary point
  combinations: do i = 1, n_lp
    dimensions: do k = 1, 3
      lp_sur(k,i) = mod((i-1)/2**(3-k),2)
    end do dimensions
  end do combinations
  opp_lp_sur(:,:) = 1 - lp_sur(:,:)

  ! Initialise MPI and set up the topology.
  call InitializeMPIcomms()

  write(msgstr,"('Starting ',A)") trim(lbeversion)
  call log_msg_hdr(msgstr)
  call log_msg(lbeplatform)
  call log_msg(lbeflags)

  ! Do not call this before MPI has been initialized
  call lbe_detect_flags()
  call lbe_report_flags()

  ! Report on which NaN check will be used, warn if none is available
  call report_check_NaN_function(msgstr)
  call log_msg_ws(msgstr)

#ifdef USEHDF
  ! Set up HDF5 if required.
  call lbe_io_init_hdf5()
#endif

  ! Parse command-line arguments
  call lbe_parse_arguments()

  ! Set up performance profiling timers for LB
  call lbe_init_timers()

! Read the fixed input.
call lbe_get_fixed_input()
! Generate unique seeds for every processor (this needs to be looked at)
call lbe_unique_seeds()
! Generate unique ID
call lbe_unique_id()

! Change from MPI_COMM_WORLD to a Cartesian grid
call ReorderMPIcomms()
! Divide the system into blocks
call lbe_divide_grid()
! initialize additional variables holding system size, etc.
call lbe_parallel_init()

! Read the lbe_input namelist.
call lbe_get_lbe_input()
! Read the boundary condition namelist.
call lbe_bc_input()

! Read the variable input.
call lbe_get_variable_input()
! Call this after lbe_get_lbe_input, lbe_get_variable_input and lbe_parse_arguments
call set_restore_parameters()

call lbe_io_halo_extent()
#ifdef MD
! read md input file - necessary to get appropriate  halo_extent and
! force_halo_extent afterwards.
call md_input()
call md_halo_extent()
call md_force_halo_extent()
#endif

if ( zeroforceoffset ) then
#ifdef SINGLEFLUID
if (SCMP) call request_halo_extent(2,"SINGLEFLUID with SCMP")
#else
call request_halo_extent(2,"NOSINGLEFLUID")
#ifdef NOSURFACTANT
!!$  if (SCMULTIRANGE) call request_halo_extent(2,"SCMULTIRANGE")
#endif
#endif
end if

#ifdef IBM_PART
call read_ibm_input()
call request_halo_extent(2, "IBM_PART") ! IBM requires a rock node halo of width >=2.
#endif

#ifdef ELEC
! Read the elec input file and process its parameters. Invalid combinations
! of parameters are also checked here.
call elec_read_input()
call request_halo_extent(3, "ELEC") ! ELEC requires a halo of width >=3.
use_lbe_force = .true.
#endif

#ifdef TRACER
call tracer_input()
call tracer_halo_extent()
#endif

  call lbe_force_input()
  call lbe_galilean_stabilizer_input()
  call lbe_mass_scaler_input()

  call set_halo_extent()

  ! Now that halo_extent has been set correctly, we can try to
  ! allocate the arena, but first make sure the halo is not too large
  ! to be filled by nearest-neighbor communication.
  call log_msg_hdr("Allocating lattice")

  if (any((/nx,ny,nz/)<halo_extent)) then
     call error('Requested halo_extent is larger than the local process domain'&
          &//'---try a more cubic decomposition or less CPUs!')
  end if
  allocate(whole_N(&
             &1-halo_extent:nx+halo_extent,&
             &1-halo_extent:ny+halo_extent,&
             &1-halo_extent:nz+halo_extent),stat=ierror)
  N => whole_N(0:nx+1,0:ny+1,0:nz+1)
  !> \warning I found no way to set \c lbound(N) to 0. Instead they
  !> are 1. However, this imposes no problem as long as nowhere in the
  !> main program there is a direct access to an element of \c
  !> N. Subroutines that receive \c N as an argument specify the
  !> starting index anyway.
  call check_allocate(ierror,'Unable to allocate arena. Aborting...')

#ifndef NOSURFACTANT
  allocate(d_adv(3,nvecs,0:nx+1,0:ny+1,0:nz+1), stat = ierror)
  call check_allocate(ierror,'Unable to allocate d_adv. Aborting...')
#endif

#ifdef AXISYM
  allocate(u(1:3,1:nx,1:ny,1:nz,1:3), stat = ierror)
  call check_allocate(ierror,'Unable to allocate u. Aborting...')
  allocate(rho(1:3,1:nx,1:ny,1:nz), stat = ierror)
  call check_allocate(ierror,'Unable to allocate rho. Aborting...')
#endif

  write(msgstr,"('Succesfully allocated lattice with halo_extent = ',I0,' .')") halo_extent
  call log_msg(msgstr)

  call msd_init()

  ! Initialise the collision code
  ! call lbe_collision_init()

  ! initialize fluxz region names before they are read from input file
  call lbe_setup_fluxz()

  call lbe_io_init()
#ifdef USEHDF
  ! If we use HDF5 we add the contents of additional input files to the metadata.
  if (myrankc .eq. 0) then
#ifdef MD
    call lbe_add_inputfile_metadata_hdf5(trim(inp_file)//".md","MD")
#endif
#ifdef ELEC
    call lbe_add_inputfile_metadata_hdf5(trim(inp_file)//".elec","ELEC")
#endif
  end if
#endif
#ifdef VARTAU
 call build_vartau_site_mpitype(vartau_mpitype)
 call build_all_chunk_mpitypes(N, vartau_halo, vartau_mpitype, halo_extent)
#endif
#ifdef LOCALBC
 call build_local_acccoef_site_mpitype(acccoef_mpitype)
 call build_all_chunk_mpitypes(N, acccoef_halo, acccoef_mpitype, halo_extent)
#endif

  ! Initialise the system.
  call lbe_init_system(N,whole_N)

  ! Write out the processor topology
  if ((myrankc == 0).and.(.not.(post))) then
    call lbe_write_topology()
    call log_msg("Wrote topology")
  end if

  ! If restoring, then don't start from timestep zero.
  if ( is_restoring() ) then
    tstart = n_restore
  else
    tstart = 0
  endif
  nt = tstart

  call lbe_bc_init(N,whole_N)
  call lbe_galilean_stabilizer_init()
  call lbe_mass_scaler_init()

#ifdef MD
  call md_init(N,whole_N,d_adv)

  ! since N might have been changed by md_init(), update the halo to
  ! be safe. This is probably needed only for MD Lees-Edwards where a
  ! proper halo exchange is not possible before MD was set up.
#ifndef NOFLUIDDYNAMICS
  call lbe_bc_complete_halo_exchange(N,whole_N,d_adv)
#endif

!!$!!$  ! This neutralizes lbe_halo_exchange_new() called by md_init() as far as
!!$!!$  ! it concerns the periodic z-boundaries of the periodic
!!$!!$  ! sub-volume. This is a hack that will become unnecessary as soon as
!!$!!$  ! the rock potentials do not require larger halos anymore. All other
!!$!!$  ! current features requiring larger halos are incompatible with
!!$!!$  ! periodic_inflow anyway.
!!$!!$  if (boundary=='periodic_inflow') &
!!$!!$       &call bc_before_advection_periodic_inflow(whole_N)
#endif

  call lbe_print_site_counts(N)

#ifdef TRACER
  call tracer_init(whole_N)
#endif

  ! The array lbe_force is allocated and initialized to 0 everywhere.
  ! After this call, lbe_force is ready to use.
  call lbe_force_init()

#ifdef ELEC
  ! Start up ELEC. The force arrays need to have been initialized,
  ! i.e. lbe_force_init() must have been called.
  ! Details of this sequence of calls can be found in the manual.
  call elec_init_timers()
  call elec_setup_parallel(whole_N)
  if ( .not. is_restoring() ) then
    call postprocess(N)
    call elec_init_system(whole_N)
    call halo_exchange(whole_N, elec_halo)
    call elec_initial_equilibration(whole_N)
    call elec_report_timers()
    call elec_reset_timers()
  else
    call elec_restore_init_system(whole_N)
  end if
#endif

  ! Initialization of the IBM plugin.
#ifdef IBM_PART
  call lagr_init(N,whole_N)
#ifdef IBM_BINARYIBM
  if(.not. restore) then
    call lagr_spread_forces()
    call lagr_recolour_interior(N)
    call lagr_dump_data(N) ! TODO: check this!
  end if
#endif


  ! The following routines have to be called if restarting from a checkpoint.
  if(restore) then
    call lbe_force_apply(N, whole_N)
    call lagr_compute_forces_int(nt)
    call lagr_distribute_nodes()
    call lagr_update_lut()
    call lagr_compute_forces_ext(N)
    call lagr_combine_forces()
    call lagr_spread_forces()
#ifdef IBM_INDEXFIELD
    call lagr_rebuild_index(N, .true.)
    call lagr_exchange_index_halo()
    call lagr_selective_wall_forces()
#endif
#ifdef IBM_BINARYIBM
    call lagr_multicomponent_forces(N)
    call lagr_recolour_interior(N)
#endif
  end if
#endif

  if ( nt >= n_sci_start ) then
     call start_timer(ti_dump)
     if (post) then
        call postprocess(N)
#ifdef ELEC
        call elec_postprocess(whole_N)
#endif
     else
        call dump_data(N)
     end if
     call dump_parallel(N,whole_N)
     call stop_timer(ti_dump)
  end if

  if (check_dump_now(sci_stress, n_sci_stress) ) then
    call log_msg("Dumping STRESSES...")
    call dump_stress(whole_N)
  end if

  if (.not.is_restoring()) call msd_run(whole_N)

  if ( n_sanity_check /= 0 ) call lbe_sanity_check(N, insanity)

  if ( insanity .or. every_n_time_steps(n_checkpoint) ) then
    call start_timer(ti_dump)
    call process_checkpoints(N, insanity)
    call stop_timer(ti_dump)
  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Let's do time loop, yeah..
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  nt = tstart + 1
  start_time_t = unixtime()
  call start_timer(ti_total)

  call log_msg_hdr("Starting time loop")

  do while (nt.le.n_iteration)
!    if(mod(nt,100).eq.0) then
      write(msgstr,"('Starting timestep ',i0)") nt
      call log_msg(msgstr)
!    end if

    ! DEBUG: Check I'm conserving everything I should do
    ! These functions are currently commented out in lbe_io and lbe_parallel
    ! call ins_total_species(N)
    ! call ins_total_momentum(N)
    ! call ins_rock(N)

#ifndef NOFLUIDDYNAMICS
    call lbe_bc_init_step()
#endif

#ifdef MD
    DEBUG_MSG("md_run_pre_before_bc")
    call md_run_pre_before_bc(whole_N)
#endif

#ifndef NOFLUIDDYNAMICS
    DEBUG_MSG("lbe_bc_before_advection")
    call lbe_bc_before_advection(N,whole_N)
#endif

#ifdef MD
    DEBUG_MSG("md_run_pre_after_bc")
    call md_run_pre_after_bc(whole_N)
#endif
    if (check_dump_now(sci_stress, n_sci_stress) ) then
      call log_msg("Dumping STRESSES...")
      call dump_stress(whole_N)
    end if
!! This is moved to the position after bc_after_advection
!    if (zeroforceoffset) then
!#ifdef SINGLEFLUID
!      if (SCMP) then
!        DEBUG_MSG("zero_force_offset")
!        call zero_force_offset(whole_N)
!      end if
!#else
!      DEBUG_MSG("zero_force_offset")
!      call zero_force_offset(whole_N)
!#endif
!    end if

    ! Propagation
    call start_timer(ti_adv)
#ifndef NOFLUIDDYNAMICS
#ifdef SINGLEFLUID
    DEBUG_MSG("lbe_advection_simple")
    call lbe_advection_simple(N, d_adv)
#else
    DEBUG_MSG("lbe_advection")
    call lbe_advection(N,d_adv)
#endif
#endif
    call stop_timer(ti_adv)

#ifndef NOFLUIDDYNAMICS
    DEBUG_MSG("lbe_bc_after_advection")
    call lbe_bc_after_advection(N,whole_N,d_adv)
#ifdef INTERPOLATEDBB
    ! compute missing reflected distributions by interpolation
    ! on fluid node
    DEBUG_MSG("lbe_interpolated_dist")
    call lbe_interpolated_dist(whole_N)
    !call lbe_interpolated_dist(lbe_N)
#endif        
#endif

!!! We should call zeroforceoffset here.
if (zeroforceoffset) then
#ifdef SINGLEFLUID
      if (SCMP) then
        DEBUG_MSG("zero_force_offset")
        call zero_force_offset(whole_N)
      end if
#else 
      DEBUG_MSG("zero_force_offset")
      call zero_force_offset(whole_N)
#endif
    end if

#ifdef TRACER
    DEBUG_MSG("tracer_run")
    call tracer_run(whole_N)
#endif

#ifdef MD
    ! give MD the opportunity to adjust to fluid configuration after streaming
    DEBUG_MSG("md_run_middle")
    call md_run_middle(N,whole_N,d_adv)
#endif

    ! Execute routines related to the IBM plug-in.
#ifdef IBM_PART
    DEBUG_MSG("IBM")
    call lagr_compute_vel_phys(N)
    call lagr_exchange_velocity_halo()
    call lagr_interpolate_velocities()
    call lagr_collect_nodes()
    call lagr_update_particles_pre()
    call lagr_exchange_particles()
    call lagr_update_particles_post()
    call lagr_dump_data(N)
#endif

    ! Reset lbe_forces.
    ! All forces added to lbe_force before this call will be overwritten.
#ifdef ELEC
    ! As ELEC adds force after lbe_add_force_halo, we need to still reset the halo as well.
    ! This is the only difference between elec_force_reset and the usual lbe_force_reset.
    DEBUG_MSG("elec_force_reset")
    call elec_force_reset()
#else
    DEBUG_MSG("lbe_force_reset")
    call lbe_force_reset()
#endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Add forces to lbe_force only after this comment. !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Add constant and/or Kolmogorov forces to lbe_force.
    DEBUG_MSG("lbe_force_apply")
    call lbe_force_apply(N, whole_N)

    ! Execute routines related to the IBM plug-in.
#ifdef IBM_PART
    DEBUG_MSG("IBM 2")
    call lagr_compute_forces_int(nt)
    call lagr_distribute_nodes()
    call lagr_update_lut()
    call lagr_compute_forces_ext(N)
    call lagr_combine_forces()
    call lagr_spread_forces()
#ifdef IBM_INDEXFIELD
    call lagr_rebuild_index(N, .false.)
    call lagr_exchange_index_halo()
    call lagr_selective_wall_forces()
#ifdef IBM_BINARYIBM
    !call lagr_recolour_interior(N)
    call lagr_multicomponent_forces(N)
#endif
#endif
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Add forces to lbe_force only before this comment. !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Exchange force halo.
    ! Add forces from the halo of neighboring processes to the physical region of lbe_force.
#ifndef NOFLUIDDYNAMICS
    DEBUG_MSG("lbe_add_force_halo")
    call lbe_add_force_halo()
#endif

#ifdef ELEC
#ifndef ELEC_NEWFORCE
    ! Contrary to the comment above, ELEC doesn't want its contribution to the 
    ! force halo exchanged and summed, so forces are added only here, *after* lbe_add_force_halo.
    DEBUG_MSG("elec_add_forces")
    call elec_add_forces(whole_N)

    ! Now we can start touching the charges.
    DEBUG_MSG("elec_timestep")
    call elec_timestep(whole_N)
#else
    ! ELEC forces do not need to be halo exchanged and summed, so they
    ! are added here *after* lbe_add_force_halo.
    ! The elec forces are calculated from the link fluxes during the
    ! elec update, so the elec timestep has to be performed first
    ! before adding the forces.
    DEBUG_MSG("elec_timestep")
    call elec_timestep(whole_N)
    DEBUG_MSG("elec_add_forces")
    call elec_add_forces(whole_N)
#endif
#endif

!!! Apply velocity to evaporation boundary before collision
!call lbe_invade_after_evaporation(N,m_evp,m_evp_gr,m_evp_gb,m_evp_freq_f,m_evp_freq_a,in_evp,out_evp,m_evp_set_density)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Shan-Chen forces and collision. !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef AXISYM
    DEBUG_MSG("lbe_get_hydrovar")
    call lbe_get_hydrovar(N,u,rho)
#endif

#ifndef NOFLUIDDYNAMICS
#ifdef SINGLEFLUID
    if ( ( collisiontype_id .ne. BGK_id ) .or. SCMP) then
      ! Ariels SINGLEFLUID MRT implemented in lbe_interactions_single.
      ! Now SINGLEFLUID MRT is reintegrated to lbe_interactions as special case of multicomponent MRT.
      DEBUG_MSG("lbe_interactions")
      call lbe_interactions_simple(N) !, whole_N, d_adv)
    else
      DEBUG_MSG("lbe_interactions_simple")
      call lbe_interactions_simple(N)
    end if
#else
    DEBUG_MSG("lbe_interactions")
    call lbe_interactions(N, whole_N, d_adv)
#endif
#endif

#ifdef AXISYM
!!! pre-post collision averaged hydrovar should be used 
    DEBUG_MSG("lbe_axisym_correction_cty")
    call lbe_axisym_correction_cty(N,u,rho)
    DEBUG_MSG("lbe_axisym_correction_mom")
    call lbe_axisym_correction_mom(N,u,rho)
!!!  boundary conditions implemented after collision
    DEBUG_MSG("lbe_axisym_axis_bc")
    call lbe_axisym_axis_bc(N)
#endif

#ifdef MD
    DEBUG_MSG("md_run_post")
    call md_run_post(whole_N)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Output of science data if required. !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    DEBUG_MSG("msd_run")
    call msd_run(whole_N)

    call start_timer(ti_dump)

    if (nt >= n_sci_start) then
      if ( post ) then
        DEBUG_MSG("postprocess")
        call postprocess(N)
#ifdef ELEC
        DEBUG_MSG("elec_postprocess")
        call elec_postprocess(whole_N)
#endif
      else
        DEBUG_MSG("dump_data")
        call dump_data(N)
      end if

      DEBUG_MSG("dump_parallel")
      call dump_parallel(N,whole_N)
    end if

    call stop_timer(ti_dump)

!!!!!!!!!!!!!!!!!!!!!
    ! Sanity check. !
!!!!!!!!!!!!!!!!!!!!!

    if ( every_n_time_steps(n_sanity_check) ) then
      DEBUG_MSG("lbe_sanity_check")
      call lbe_sanity_check(N, insanity)
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Potentially remove center of mass motion !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call lbe_galilean_stabilizer_run(whole_N)
    call lbe_mass_scaler_run(whole_N)

!!!!!!!!!!!!!!!!!!!!!!!!
    ! Checkpoint dump. !
!!!!!!!!!!!!!!!!!!!!!!!!

    ! Dump checkpoint.
    if ( insanity .or. every_n_time_steps(n_checkpoint) ) then
      call start_timer(ti_dump)
      DEBUG_MSG("process_checkpoints")
      call process_checkpoints(N, insanity)
      call stop_timer(ti_dump)
    end if

    ! Abort simulation if it is not sane.
    if (insanity) call Abend

    ! Increase the time step.
    nt = nt + 1

    ! timesteps_count always counts the no of timesteps performed
    ! in *this* run of the simulation, rather than measuring
    ! the time coordinate from simulation start.
    timesteps_count = timesteps_count + 1
  end do

  ! CALL set_endtime()
  end_time_t = unixtime()
  call stop_timer(ti_total)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Finished - now tidy up.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef ELEC
  call elec_shutdown()
#endif
#ifdef MD
  call md_shutdown()
#endif
#ifdef TRACER
  call tracer_shutdown()
#endif
  call sum_timer( (/ stdout /)) ! dump timer evaluation to standard output

#ifdef ELEC
  call elec_report_timers()
#endif

  call lbe_bc_shutdown()

#ifdef USEHDF
  call lbe_io_shutdown_hdf5()
#endif
  call FinalizeMPIcomms()

  call log_msg_ws("Exiting ...")

end program lbe
#endif
