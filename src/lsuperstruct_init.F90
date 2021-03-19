!> Lagrangian superstructure initialization module
!>
!> This module contains all subroutines responsible for initializing the Lagrangian superstructure.

#include "lbe.h"

module lsuperstruct_init_module

#ifdef IBM_PART

  ! Include external modules.
  use lbe_globals_module, only : halo_extent, rk, use_lbe_force, myrankc
  use lbe_helper_module, only : norm, n_sites_fluid, check_nan, is_restoring
  use lbe_init_rock_module, only: inside_rock_global, rock_is_present
  use lbe_io_helper_module, only : lbe_make_filename_restore_rank
  use lbe_log_module
  use lbe_parallel_module, only : Abend, check_allocate, comm_cart, gather_rock_state, tnx, tny, tnz
  use lbe_parms_module, only : nt, nx, ny, nz, inp_file, boundary_cond
  use lbe_types_module, only : lbe_site
  use lsuperstruct_data_module
  use lsuperstruct_dump_module, only : write_particles_vtk
  use lsuperstruct_IBM_module, only : IBM_spread_forces
  use lsuperstruct_interface_module, only : set_rock_grad
  use lsuperstruct_parallel_module
  use lsuperstruct_helper_module, only : distance_vector, add_particle
#ifdef IBM_INDEXFIELD
  use lsuperstruct_timeloop_module, only : rebuild_interior_index
  use lextobj_module, only : update_node_normals
#endif
  use lmesh_module, only : meshes, init_mesh_master, init_mesh_slave, allocate_mesh_memory
  use lextobj_module, only : particles, part_ind, init_particle, set_radius, set_center_position, rotate_particle, &
        & forces_compute, forces_friction, update_face_areas, update_node_areas
  use lsuperstruct_timeloop_module, only : compute_particle_forces_ext, combine_forces, &
        & update_particles_pre, update_particles_post, compute_particle_particle_forces, &
        & compute_wall_gradient_forces, update_node_interaction_LUT
  use lbe_io_xdrf_module

  implicit none
  include 'mpif.h'

  private
  public :: init_read_parameters, init_lsuperstruct, check_parameter_conflicts

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read parameters from input file
  !>
  !> All parameters which are not already contained in the "standard" LB3D parameter files are read here.
  !> Rank 0 reads the data and broadcasts the parameters.

  subroutine init_read_parameters()
    ! Declare variables.
    integer :: ierror ! error code

    ! Read parameters.
    if(myrankc == 0) then
      call log_msg("-----------( Reading IBM input )-----------", .false.)
      open(unit=42, file=inp_file, status='UNKNOWN')
      read(unit=42, nml=ibm_input)
#ifdef IBM_SWIMMER
      read(unit=42, nml=ibm_swimmers)
#endif 
      close(unit=42)

      call log_msg("", .false.)
      write(msgstr, "('IBM_RANGE = ', i0)") IBM_RANGE
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('mode_check_rock_pos = ', l1)") mode_check_rock_pos
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('mode_polydispersity = ', l1)") mode_polydispersity
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('mode_manual_pos = ', l1)") mode_manual_pos
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('mode_grow_particles = ', l1)") mode_grow_particles
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('rel_initial_radius = ', f0.3)") rel_initial_radius
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('num_steps_growth = ', i0)") num_steps_growth
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('growth_bending_modulus = ', f0.3)") growth_bending_modulus
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('dump_vtk = ', l1)") dump_vtk
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('time_step_dump_vtk = ', i0)") time_step_dump_vtk
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('dump_profiles = ', l1)") dump_profiles
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('time_step_dump_profiles = ', i0)") time_step_dump_profiles
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('dump_particles = ', l1)") dump_particles
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('time_step_dump_particles = ', i0)") time_step_dump_particles
      call log_msg(trim(msgstr), .false.)
#ifdef IBM_INDEXFIELD
      write(msgstr, "('dump_indexfield = ', l1)") dump_indexfield
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('time_step_dump_indexfield = ', i0)") time_step_dump_indexfield
      call log_msg(trim(msgstr), .false.)
#endif
      write(msgstr, "('num_meshes = ', i0)") num_meshes
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('int_strength_part_part = ', f0.3)") int_strength_part_part
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('int_strength_part_wall = ', f0.3)") int_strength_part_wall
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('int_strength_gradient = ', f0.3)") int_strength_gradient
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('friction = ', f0.3)") friction
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('IBM_massdensity = ', f0.3)") IBM_massdensity
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('surften = ', f0.3)") surften
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('contactangle = ', f0.3)") contactangle
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('width = ', f0.3)") width
      call log_msg(trim(msgstr), .false.)
#ifdef IBM_INDEXFIELD
      write(msgstr, "('viscosity_contrast = ', f0.3)") viscosity_contrast
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('time_step_rebuild_index = ', i0)") time_step_rebuild_index
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('selective_wall_force_magnitude = ', f0.3)") selective_wall_force_magnitude
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('selective_wall_force_range = ', f0.3)") selective_wall_force_range
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('g_ri = ', f0.3)") g_ri
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('g_bi = ', f0.3)") g_bi
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('noslipcoeff = ', f0.3)") noslipcoeff
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('ibm_colour = ', f0.3)") ibm_colour
      call log_msg(trim(msgstr), .false.)
#endif
#ifdef IBM_SWIMMER
      write(msgstr, "('swimmer_d12 = ', f0.3)") swimmer_d12
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('swimmer_d23 = ', f0.3)") swimmer_d23
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('swimmer_modulus12 = ', f0.3)") swimmer_modulus12
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('swimmer_modulus23 = ', f0.3)") swimmer_modulus23
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('swimmer_kb = ', f0.3)") swimmer_kb
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('swimmer_ka = ', f0.3)") swimmer_ka
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('swimmer_anchor = ', l1)") swimmer_anchor
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('swimmer_force_mag1 = ', f0.3)") swimmer_force_mag1
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('swimmer_force_mag2 = ', f0.3)") swimmer_force_mag2
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('swimmer_force_mag3 = ', f0.3)") swimmer_force_mag3
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('swimmer_period1 = ', f0.3)") swimmer_period1
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('swimmer_period2 = ', f0.3)") swimmer_period2
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('swimmer_period3 = ', f0.3)") swimmer_period3
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('swimmer_phase1 = ', f0.3)") swimmer_phase1
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('swimmer_phase2 = ', f0.3)") swimmer_phase2
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('swimmer_phase3 = ', f0.3)") swimmer_phase3
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('swimmer_manualpos = ', l1)") swimmer_manualpos
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('swimmer_pos = (',F16.10,',',F16.10,',',F16.10,')')") swimmer_pos(1), swimmer_pos(2), swimmer_pos(3)
      call log_msg(trim(msgstr), .false.)
      write(msgstr, "('swimmer_axis = (',F16.10,',',F16.10,',',F16.10,')')") swimmer_axis(1), swimmer_axis(2), swimmer_axis(3)
      call log_msg(trim(msgstr), .false.)
#endif
      call log_msg("", .false.)
    endif

    ! Broadcast parameters to remaining processes.
    call MPI_Bcast(IBM_RANGE, 1, MPI_INTEGER, 0, comm_cart, ierror)
    call MPI_Bcast(mode_check_rock_pos, 1, MPI_LOGICAL, 0, comm_cart, ierror)
    call MPI_Bcast(mode_polydispersity, 1, MPI_LOGICAL, 0, comm_cart, ierror)
    call MPI_Bcast(mode_manual_pos, 1, MPI_LOGICAL, 0, comm_cart, ierror)
    call MPI_Bcast(mode_grow_particles, 1, MPI_LOGICAL, 0, comm_cart, ierror)
    call MPI_Bcast(rel_initial_radius, 1, MPI_REAL8, 0, comm_cart, ierror)
    call MPI_Bcast(num_steps_growth, 1, MPI_INTEGER, 0, comm_cart, ierror)
    call MPI_Bcast(growth_bending_modulus, 1, MPI_REAL8, 0, comm_cart, ierror)
    call MPI_Bcast(dump_vtk, 1, MPI_LOGICAL, 0, comm_cart, ierror)
    call MPI_Bcast(time_step_dump_vtk, 1, MPI_INTEGER, 0, comm_cart, ierror)
    call MPI_Bcast(dump_profiles, 1, MPI_LOGICAL, 0, comm_cart, ierror)
    call MPI_Bcast(time_step_dump_profiles, 1, MPI_INTEGER, 0, comm_cart, ierror)
    call MPI_Bcast(dump_particles, 1, MPI_LOGICAL, 0, comm_cart, ierror)
    call MPI_Bcast(time_step_dump_particles, 1, MPI_INTEGER, 0, comm_cart, ierror)
    call MPI_Bcast(num_meshes, 1, MPI_INTEGER, 0, comm_cart, ierror)
    call MPI_Bcast(int_strength_part_part, 1, MPI_REAL8, 0, comm_cart, ierror)
    call MPI_Bcast(int_strength_part_wall, 1, MPI_REAL8, 0, comm_cart, ierror)
    call MPI_Bcast(int_strength_gradient, 1, MPI_REAL8, 0, comm_cart, ierror)
    call MPI_Bcast(friction, 1, MPI_REAL8, 0, comm_cart, ierror)
    call MPI_Bcast(IBM_massdensity, 1, MPI_REAL8, 0, comm_cart, ierror)
    call MPI_Bcast(surften, 1, MPI_REAL8, 0, comm_cart, ierror)
    call MPI_Bcast(contactangle, 1, MPI_REAL8, 0, comm_cart, ierror)
    call MPI_Bcast(width, 1, MPI_REAL8, 0, comm_cart, ierror)
#ifdef IBM_INDEXFIELD
    call MPI_Bcast(viscosity_contrast, 1, MPI_REAL8, 0, comm_cart, ierror)
    call MPI_Bcast(time_step_rebuild_index, 1, MPI_INTEGER, 0, comm_cart, ierror)
    call MPI_Bcast(selective_wall_force_magnitude, 1, MPI_REAL8, 0, comm_cart, ierror)
    call MPI_Bcast(selective_wall_force_range, 1, MPI_REAL8, 0, comm_cart, ierror)
    call MPI_Bcast(g_ri, 1, MPI_REAL8, 0, comm_cart, ierror)
    call MPI_Bcast(g_bi, 1, MPI_REAL8, 0, comm_cart, ierror)
    call MPI_Bcast(noslipcoeff, 1, MPI_REAL8, 0, comm_cart, ierror)
    call MPI_Bcast(ibm_colour, 1, MPI_REAL8, 0, comm_cart, ierror)
    call MPI_Bcast(dump_indexfield, 1, MPI_LOGICAL, 0, comm_cart, ierror)
    call MPI_Bcast(time_step_dump_indexfield, 1, MPI_INTEGER, 0, comm_cart, ierror)
#endif
#ifdef IBM_SWIMMER
    call MPI_Bcast(swimmer_d12, 1, MPI_REAL8, 0, comm_cart, ierror)
    call MPI_Bcast(swimmer_d23, 1, MPI_REAL8, 0, comm_cart, ierror)
    call MPI_Bcast(swimmer_modulus12, 1, MPI_REAL8, 0, comm_cart, ierror)
    call MPI_Bcast(swimmer_modulus23, 1, MPI_REAL8, 0, comm_cart, ierror)
    call MPI_Bcast(swimmer_kb, 1, MPI_REAL8, 0, comm_cart, ierror)
    call MPI_Bcast(swimmer_ka, 1, MPI_REAL8, 0, comm_cart, ierror)
    call MPI_Bcast(swimmer_anchor, 1, MPI_LOGICAL, 0, comm_cart, ierror)
    call MPI_Bcast(swimmer_force_mag1, 1, MPI_REAL8, 0, comm_cart, ierror)
    call MPI_Bcast(swimmer_force_mag2, 1, MPI_REAL8, 0, comm_cart, ierror)
    call MPI_Bcast(swimmer_force_mag3, 1, MPI_REAL8, 0, comm_cart, ierror)
    call MPI_Bcast(swimmer_period1, 1, MPI_REAL8, 0, comm_cart, ierror)
    call MPI_Bcast(swimmer_period2, 1, MPI_REAL8, 0, comm_cart, ierror)
    call MPI_Bcast(swimmer_period3, 1, MPI_REAL8, 0, comm_cart, ierror)
    call MPI_Bcast(swimmer_phase1, 1, MPI_REAL8, 0, comm_cart, ierror)
    call MPI_Bcast(swimmer_phase2, 1, MPI_REAL8, 0, comm_cart, ierror)
    call MPI_Bcast(swimmer_phase3, 1, MPI_REAL8, 0, comm_cart, ierror)
    call MPI_Bcast(swimmer_manualpos, 1, MPI_LOGICAL, 0, comm_cart, ierror)
    call MPI_Bcast(swimmer_pos, 3, MPI_REAL8, 0, comm_cart, ierror)
    call MPI_Bcast(swimmer_axis, 3, MPI_REAL8, 0, comm_cart, ierror)
#endif
  end subroutine init_read_parameters

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialization of the Lagrangian superstructure
  !>
  !> \param[in,out] lbe_N local chunk of the lattice with halo extent
  !> 1 (old LB3D style)
  !>
  !> \param[in] whole_N local chunk of the lattice with full halo of
  !> depth \c halo_extent
  !>
  !> All required initialization routines are called from here before the simulation loop starts.
  !>  1) The look-up table for node-node interactions and array for the lattice velocity are initialized.
  !>  2) Derived datatypes for communication are created (part 1)
  !>  3) Polydispersity is initialized (if required).
  !>  4) Meshes are initialized.
  !>  5) Derived datatypes for communication are created (part 2)
  !>  6) Particles are initialized (memory allocation)
  !>  7) Particles are positioned (if no restore) or restored (if restore)
  !>  8) The volume fraction is computed.
  !>  9) Particles are grown to their full size (if required).
  !> 10) Particles are manipulated if required for benchmarking or specific simulations.
  !> 11) If in IBM_INDEXFIELD mode, the index field has to be computed.
  !> This subroutine is the only one visible from the outside.
  !> It is called by \c lbe.F90 during simulation initialization.

  subroutine init_lsuperstruct(lbe_N, whole_N)
    type(lbe_site), dimension(0:, 0:, 0:), intent(inout) :: lbe_N !< lattice
    type(lbe_site),intent(in) :: whole_N(1-halo_extent:,1-halo_extent:,1-halo_extent:)

    ! Set the halo width required for IBM communication.
    ! Note that IBM_RANGE is an integer which may be odd. In this case, the integer division still is correct.
    IBM_HALO = IBM_RANGE / 2

    ! Activate lattice force inclusion in LB3D.
    ! This is of central importance since the IBM force coupling would not work otherwise.
    use_lbe_force = .true.

    ! Call initialization routines.
    call init_lattice_quantities() ! initialize lattice-based quantities (physical velocity, look-up table)
    call create_datatypes_pre() ! create derived MPI datatypes

    if( is_restoring() .eqv. .false.) then
      call set_rock_grad(lbe_N) ! find gradient of rock nodes
    end if

    call init_polydispersity() ! initialize polydispersity (if required)
    call init_meshes() ! initialize meshes
    call create_datatypes_post() ! create derived MPI datatypes
    call init_particles() ! allocate memory for particles

    if(is_restoring() .eqv. .false.) then
      call position_particles(whole_N) ! position particles in domain
    else
      call restore_particles() ! restore particles from a checkpoint
    end if

    call compute_volume_fraction() ! compute volume fraction of the particles
    call grow_particles(lbe_N) ! grow particles to their full size (if required)
    call manipulate_particles() ! manipulate particles

#ifdef IBM_INDEXFIELD
    call distribute_nodes()
    call IBM_spread_forces(.false.)
    call rebuild_interior_index(lbe_N, .true.)
#endif
  end subroutine init_lsuperstruct

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check conflicts of known parameter/flag combinations
  !>
  !> This subroutine is supposed to help the user to identify potential problems caused by
  !> improper combination of simulation parameters and flags.
  !> The developer is requested to keep this up to date.
  !> TODO: This list should be extended at some point.

  subroutine check_parameter_conflicts()
    ! Declare variables.
    character(len=200) :: message ! message string

  ! Selective wall force and IBM_INDEXFIELD
#ifdef IBM_INDEXFIELD
  if((selective_wall_force_magnitude .ne. 0.00) .and. (boundary_cond .ne. 6)) then
    write(message, "('Parameter conflict found:')")
    call log_msg(trim(message), .false.)
    write(message, "('  boundary_cond = ', i0, '; selective_wall_force_magnitude = ', F0.3)") boundary_cond, &
        & selective_wall_force_magnitude
    call log_msg(trim(message), .false.)
    write(message, "('  Selective wall force would have no effect. boundary_cond = 6 required.')")
    call log_msg(trim(message), .false.)
    call Abend()
  end if
#endif

  ! Viscosity contrast and VARTAU
#if defined IBM_INDEXFIELD && !defined VARTAU
  if(viscosity_contrast .ne. 1.0d0) then
    write(message, "('Parameter conflict found:')")
    call log_msg(trim(message), .false.)
    write(message, "('  viscosity_contrast = ', F0.3, ' found, but VARTAU not defined')") viscosity_contrast
    call log_msg(trim(message), .false.)
    write(message, "('  Viscosity contrast would have no effect. viscosity_contrast = 1 required.')")
    call log_msg(trim(message), .false.)
    call Abend()
  end if
#endif

  end subroutine check_parameter_conflicts

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialization of lattice-based quantities
  !>
  !> All lattice-based quantities are initialized in this subroutine.
  !> - physical lattice velocity
  !>   This velocity is the true velocity of the fluid and is required for the IBM coupling.
  !>   It is the first moment of the populations, corrected by half of the forces.
  !>   NOTES:
  !>   A halo is required to enable proper velocity treatment during IBM interpolations.
  !>   The halo size depends on the interaction of the IBM interpolations.
  !>   For the 2- and 3-point interpolations, the halo width is 1, for the 4-point interpolation, it is 2.
  !> - look-up table for node-node interactions
  !>   The look-up table consists of two parts:
  !>   1) The number of nodes in a subbox is stored in the 0-entry.
  !>   2) The global node indices in the subbox are stored in the subsequent entries.
  !>   NOTES:
  !>   A halo is required to enable proper force computation near process boundaries.
  !>   The halo size depends on the interaction force range.
  !>   For a range of 1 lattice unit, two halo nodes along each direction are sufficient.
  !>   MAX_NODES_SUBBOX is the maximum number of nodes which can be located in one subbox.
  !>   Lower halo nodes are stored in array indices <1.
  !> TODO: This subroutine should be revised after the remaining code has been cleaned.

  subroutine init_lattice_quantities()
    ! Declare variables.
    integer :: status ! status for checking allocation

    ! Allocate memory for the physical velocity.
    ! The first index is for the three components of the velocity.
    ! The remaining three indices are for the lattice dimensions.
    allocate(vel_phys(1:3, 1 - IBM_HALO:nx + IBM_HALO, 1 - IBM_HALO:ny + IBM_HALO, 1 - IBM_HALO:nz + IBM_HALO), stat=status)

    ! Check allocation status.
    if(status .ne. 0) then
      call error("memory for physical velocity could not be allocated")
    endif

    vel_phys(:, :, :, :) = 0.0d0

#ifdef IBM_INDEXFIELD
    ! Allocate memory for the interior/exterior index.
    ! The three indices are for the lattice dimensions.
    allocate(interior_index(1 - IBM_HALO:nx + IBM_HALO, 1 - IBM_HALO:ny + IBM_HALO, 1 - IBM_HALO:nz + IBM_HALO), stat=status)

    ! Check allocation status.
    if(status .ne. 0) then
      call error("memory for interior index could not be allocated")
    endif

    interior_index(:, :, :) = 0.0d0

    ! Allocate memory for the normal distance.
    ! The three indices are for the lattice dimensions.
    allocate(normal_distance(1:nx, 1:ny, 1:nz), stat=status)

    ! Check allocation status.
    if(status .ne. 0) then
      call error("memory for normal distance could not be allocated")
    endif

    normal_distance(:, :, :) = 1000.0d0

    ! Allocate memory for the check status array.
    ! The three indices are for the lattice dimensions.
    allocate(checked_index(1:nx, 1:ny, 1:nz), stat=status)

    ! Check allocation status.
    if(status .ne. 0) then
      call error("memory for status check array could not be allocated")
    endif

    checked_index(:, :, :) = 0.0d0
#endif

#ifdef IBM_BINARYIBM
    allocate(interior_index_old(1:nx, 1:ny, 1:nz), stat=status)

    ! Check allocation status.
    if(status .ne. 0) then
      call error("memory for old interior index could not be allocated")
    endif

    interior_index_old(:, :, :) = 0.0d0

    allocate(vel_red(1:3, 1 - IBM_HALO:nx + IBM_HALO, 1 - IBM_HALO:ny + IBM_HALO, 1 - IBM_HALO:nz + IBM_HALO), stat=status)

    ! Check allocation status.
    if(status .ne. 0) then
      call error("memory for red velocity could not be allocated")
    endif

    vel_red(:, :, :, :) = 0.0d0

    allocate(vel_blue(1:3, 1 - IBM_HALO:nx + IBM_HALO, 1 - IBM_HALO:ny + IBM_HALO, 1 - IBM_HALO:nz + IBM_HALO), stat=status)

    ! Check allocation status.
    if(status .ne. 0) then
      call error("memory for blue velocity could not be allocated")
    endif

    vel_blue(:, :, :, :) = 0.0d0

    allocate(vel_diff(1:3, 1 - IBM_HALO:nx + IBM_HALO, 1 - IBM_HALO:ny + IBM_HALO, 1 - IBM_HALO:nz + IBM_HALO), stat=status)

    ! Check allocation status.
    if(status .ne. 0) then
      call error("memory for velocity difference could not be allocated")
    endif

    vel_diff(:, :, :, :) = 0.0d0
#endif

    ! Allocate memory for the IBM force.
    ! The first index is for the three components of the velocity.
    ! The remaining three indices are for the lattice dimensions (no halo required).
    allocate(force_IBM(1:3, 1:nx, 1:ny, 1:nz), stat=status)

    ! Check allocation status.
    if(status .ne. 0) then
      call error("memory for IBM force could not be allocated")
    endif

    force_IBM(:, :, :, :) = 0.0d0

    ! Allocate memory for look-up table.
    allocate(LUT_nodes_in_range(1 - HALO_WIDTH:nx + HALO_WIDTH, 1 - HALO_WIDTH:ny + HALO_WIDTH, 1 - HALO_WIDTH:nz + HALO_WIDTH, &
      & 0:MAX_NODES_SUBBOX), stat=status)

    ! Check allocation status.
    if(status .ne. 0) then
      call error("memory for node look-up table could not be allocated")
    endif

    LUT_nodes_in_range(:, :, :, :) = 0

    ! Allocate memory for rock gradient.
    allocate(rock_grad(1 - IBM_HALO:nx + IBM_HALO, 1 - IBM_HALO:ny + IBM_HALO, 1 - IBM_HALO:nz + IBM_HALO), stat=status)

    ! Check allocation status.
    if(status .ne. 0) then
      call error("memory for rock gradient could not be allocated")
    endif

    rock_grad(:, :, :) = 0

    ! Allocate memory for rock gradient.
    allocate(gradient(3, 0:nx, 0:ny, 0:nz), stat=status)

    ! Check allocation status.
    if(status .ne. 0) then
      call error("memory for force gradient could not be allocated")
    endif

    gradient(:, :, :, :) = 0.0d0
  end subroutine init_lattice_quantities

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialization of polydispersity
  !>
  !> If in polydispersity mode, some data from the parameter files are overwritten:
  !> - The number of meshes is computed based on the desired average radius and polydispersity of the particles.
  !> - The number of particles is computed based on the desired volume fraction.
  !> In any case, the following connections have to be clarified:
  !> - Each mesh has to know the number of particles belonging to it.
  !> - Each mesh has to know its filename.
  !> - Each particle has to obtain its mesh index.
  !> - Each particle has to obtain its radius.
  !> TODO:
  !> - Add polydispersity functionality.

  subroutine init_polydispersity()
    ! Declare variables.
    integer :: m_i ! mesh index
    integer :: c_i ! particle index
    integer :: n ! counter
    integer :: ierror ! MPI error code
    integer :: error_stat ! switch for successful line reading
    integer, parameter :: init_file_unit = 13 ! input file unit
    logical :: file_exist ! switch for file existence
    character(len=1000) :: line_buffer ! buffer for reading lines from file

    ! Report.
    call log_msg("creating relations between particles and meshes...", .false.)

    ! Allocate arrays for particle-mesh membership.
    ! These arrays are only required for rank 0 at the beginning of the simulation.
    ! They should be deallocated at the end of the initialization routines.
    if(myrankc .eq. 0) then
      ! Allocate and set arrays depending on number of meshes.
      allocate(mesh_data(num_meshes))

      ! Switch off polydispersity mode.
      ! TODO: take out this temporary measure when polydispersity functionality is added
      if(mode_polydispersity .eqv. .true.) mode_polydispersity = .false.

      ! If not in polydispersity mode, the mesh data has to be read from a parameter file.
      ! TODO: add polydispersity functionality
      if(mode_polydispersity .eqv. .false.) then
        ! Report.
        call log_msg("  reading mesh data from parm_meshes.dat...", .false.)

	if(num_meshes .ne. 0) then
	  ! Check whether mesh file exists or not.
	  inquire(file="parm_meshes.dat", exist=file_exist)

	  ! If it does not exist, terminate the simulation.
	  if(file_exist .eqv. .false.) then
	    close(init_file_unit)
	    call error("file parm_meshes.dat not found, terminating simulation")
	  end if

	  ! Open the mesh parameter file.
	  open(init_file_unit, file="parm_meshes.dat", status='old', action='read')

	  ! Skip header line.
	  read(init_file_unit, '(a999)') line_buffer

	  ! Run over all meshes.
	  do m_i = 1, num_meshes
#ifndef IBM_FIXED
	    read(init_file_unit, fmt=*, iostat=error_stat) mesh_data(m_i)%filename, mesh_data(m_i)%num_particles, &
	      & mesh_data(m_i)%radius, mesh_data(m_i)%k_v, mesh_data(m_i)%k_at, mesh_data(m_i)%k_s, &
	      & mesh_data(m_i)%k_al, mesh_data(m_i)%k_b
#else
	    read(init_file_unit, fmt=*, iostat=error_stat) mesh_data(m_i)%filename, mesh_data(m_i)%num_particles, &
	      & mesh_data(m_i)%radius, mesh_data(m_i)%k_v, mesh_data(m_i)%k_at, mesh_data(m_i)%k_s, &
	      & mesh_data(m_i)%k_al, mesh_data(m_i)%k_b, mesh_data(m_i)%k_anchor
#endif

	    ! If the data for the mesh cannot be read from the file, terminate the simulation.
	    if(error_stat .ne. 0) then
	      close(init_file_unit)
	      call error("file parm_meshes.dat cannot be read or incorrectly set up, terminating simulation")
	    end if
	  end do

	  close(init_file_unit)
	end if
      end if

      ! Report.
      call log_msg("  counting particles...", .false.)

      ! Find number of particles.
      num_particles_gl = 0

      do m_i = 1, num_meshes
        num_particles_gl = num_particles_gl + mesh_data(m_i)%num_particles
      end do

      ! Allocate and set array depending on number of particles.
      allocate(particle_data(num_particles_gl))
    end if

    ! Broadcast global number of particles and meshes.
    call MPI_Bcast(num_particles_gl, 1, MPI_INTEGER, 0, comm_cart, ierror)

    ! Report.
    call log_msg("  setting up particle properties...", .false.)

    ! Assign meshes to particles.
    ! Rank 0 is responsible for assigning each particle with its properties (mesh, radius, elasticity).
    ! This step can be skipped when restoring from a checkpoint because the corresponding data is contained in the checkpoint.
    ! TODO: add polydispersity functionality
    if((myrankc .eq. 0) .and. (mode_polydispersity .eqv. .false.) .and. ( is_restoring() .eqv. .false.)) then
      ! Run over all particles and assign their properties.
      c_i = 1

      do m_i = 1, num_meshes
        do n = 1, mesh_data(m_i)%num_particles
          particle_data(c_i)%mesh_type = m_i
          particle_data(c_i)%radius = mesh_data(m_i)%radius
          particle_data(c_i)%k_v = mesh_data(m_i)%k_v
          particle_data(c_i)%k_at = mesh_data(m_i)%k_at
          particle_data(c_i)%k_s = mesh_data(m_i)%k_s
          particle_data(c_i)%k_al = mesh_data(m_i)%k_al
          particle_data(c_i)%k_b = mesh_data(m_i)%k_b
#ifdef IBM_FIXED
          particle_data(c_i)%k_anchor = mesh_data(m_i)%k_anchor
					if( mesh_data(m_i)%k_anchor .lt. 0.0d0 ) then
						call error("Negative spring constant read from parm_meshes.dat: ensure that all mesh-types have a non-negative spring constant defined")
					end if
#endif
          ! Check particle radius: If too large, abort the simulation.
          if(2.0d0 * particle_data(c_i)%radius > min(nx, ny, nz)) then
            call error("at least one particle radius is too large for subdomain size: recheck process distribution along the axes")
          end if

          c_i = c_i + 1
        end do
      end do
    end if
  end subroutine init_polydispersity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialization of the meshes
  !>
  !> The meshes are initialized.
  !> 1) The total number of particles is computed. For polydispersity mode, it is already known.
  !> 2) The mesh array is allocated, and the meshes are created.
  !> 3) The meshes are written to the disk for restart.

  subroutine init_meshes()
    ! Declare variables.
    integer :: c_i ! particle index
    integer :: m_i ! mesh index
    integer :: n_i ! node index
    integer :: f_i ! face index
    integer :: n ! counter
    integer :: ierror ! MPI error code
    integer :: num_nodes_max ! maximum number of nodex in any mesh
    integer :: num_faces_max ! maximum number of faces in any mesh
    integer, dimension(:), allocatable :: node_face_buffer ! buffer for sender face-node-LUT
    integer, dimension(:), allocatable :: num_nodes_in_mesh ! number of nodes in each mesh
    integer, dimension(:), allocatable :: num_faces_in_mesh ! number of faces in each mesh
    real(kind=rk), dimension(:), allocatable :: node_pos_buffer ! buffer for sending node positions
    character(len=200) :: message ! message string

    ! Report progress.
    write(message, "('initializing ', i0, ' Lagrangian mesh(es)...')") num_meshes
    call log_msg(trim(message), .false.)

    ! Allocate mesh array.
    ! The number of meshes is already known, either from the parameter file or from init_polydispersity().
    allocate(meshes(num_meshes))

    ! Report number of particles in each mesh and total number of particles.
    if(myrankc .eq. 0) then
      do m_i = 1, num_meshes
        write(message, "('  number of particles for mesh ', i0, ': ', i0)") m_i, mesh_data(m_i)%num_particles
        call log_msg(trim(message), .false.)
      end do

      write(message, "('  total number of particles: ', i0)") num_particles_gl
      call log_msg(trim(message), .false.)
    end if

#ifdef IBM_SWIMMER
    if(num_particles_gl .ne. 3) then
      write(message, "('Parameter conflict found:')")
      call log_msg(trim(message), .false.)
      write(message, "('  exactly 3 particles required if IBM_SWIMMER is used, but found ', i0)") num_particles_gl
      call log_msg(trim(message), .false.)
      call Abend()
    end if
#endif

    ! Only rank 0 reads the data from the mesh files.
    ! This way, excessive HDD use is avoided.
    ! Later on, rank 0 broadcasts the data to the remaining ranks.
    if(myrankc .eq. 0) then
      do m_i = 1, num_meshes
        ! Correct the mesh file name.
        ! All mesh files are expected to by located in the relative directory ../msh.
        mesh_data(m_i)%filename = "../msh/" // mesh_data(m_i)%filename

        ! Ask the mesh module to initialize the mesh.
        call init_mesh_master(mesh_data(m_i)%filename, meshes(m_i), m_i)
      end do
    end if

    ! At this point, rank 0 knows everything about all meshes.
    ! The next step is to tell the other processes about the sizes of the meshes (number of nodes and faces in each mesh).
    ! For this, two arrays are used which store the size of each mesh.
    ! Rank 0 fills these arrays with data and broadcasts them.
    ! The remaining processes can then allocate the memory within the mesh structure.
    call log_msg("broadcasting meshes...", .false.)
    allocate(num_nodes_in_mesh(num_meshes))
    allocate(num_faces_in_mesh(num_meshes))

    if(myrankc .eq. 0) then
      do m_i = 1, num_meshes
        num_nodes_in_mesh(m_i) = meshes(m_i)%num_nodes
        num_faces_in_mesh(m_i) = meshes(m_i)%num_faces
      end do
    end if

    call MPI_Bcast(num_nodes_in_mesh, num_meshes, MPI_INTEGER, 0, comm_cart, ierror)
    call MPI_Bcast(num_faces_in_mesh, num_meshes, MPI_INTEGER, 0, comm_cart, ierror)

    if(myrankc .ne. 0) then
      do m_i = 1, num_meshes
        call allocate_mesh_memory(meshes(m_i), num_nodes_in_mesh(m_i), num_faces_in_mesh(m_i))
      end do
    end if

    ! The next step is to broadcast the mesh data.
    ! The simplest way is to have one broadcast for each single mesh.
    ! The mesh buffer can be reused, but it must be large enough for the largest mesh.
    ! Rank 0 identifies the largest number of nodes and faces and broadcasts the results.
    ! These values will be used to allocate memory for the communication buffer.
    if(myrankc .eq. 0) then
      num_nodes_max = 0
      num_faces_max = 0

      do m_i = 1, num_meshes
        if(meshes(m_i)%num_nodes .gt. num_nodes_max) num_nodes_max = meshes(m_i)%num_nodes
        if(meshes(m_i)%num_faces .gt. num_faces_max) num_faces_max = meshes(m_i)%num_faces
      end do
    end if

    call MPI_Bcast(num_nodes_max, 1, MPI_INTEGER, 0, comm_cart, ierror)
    call MPI_Bcast(num_faces_max, 1, MPI_INTEGER, 0, comm_cart, ierror)

    allocate(node_pos_buffer(3 * num_nodes_max))
    allocate(node_face_buffer(3 * num_faces_max))

    ! The third step is to
    ! - fill the buffer with the relevant data (rank 0)
    ! - broadcast the buffers
    ! - copy the buffer to the meshes (rank not 0)
    ! - initialize the meshes based on the received data (rank not 0)
    do m_i = 1, num_meshes
      if(myrankc .eq. 0) then
        do n_i = 1, num_nodes_in_mesh(m_i)
          node_pos_buffer(3 * (n_i - 1) + 1) = meshes(m_i)%pos(1, n_i)
          node_pos_buffer(3 * (n_i - 1) + 2) = meshes(m_i)%pos(2, n_i)
          node_pos_buffer(3 * (n_i - 1) + 3) = meshes(m_i)%pos(3, n_i)
        end do

        do f_i = 1, num_faces_in_mesh(m_i)
          node_face_buffer(3 * (f_i - 1) + 1) = meshes(m_i)%neighbor_face_node(f_i, 1)
          node_face_buffer(3 * (f_i - 1) + 2) = meshes(m_i)%neighbor_face_node(f_i, 2)
          node_face_buffer(3 * (f_i - 1) + 3) = meshes(m_i)%neighbor_face_node(f_i, 3)
        end do
      end if

      call MPI_Bcast(node_pos_buffer, 3 * num_nodes_in_mesh(m_i), MPI_REAL8, 0, comm_cart, ierror)
      call MPI_Bcast(node_face_buffer, 3 * num_faces_in_mesh(m_i), MPI_INTEGER, 0, comm_cart, ierror)

      if(myrankc .ne. 0) then
        do n_i = 1, num_nodes_in_mesh(m_i)
          meshes(m_i)%pos(1, n_i) = node_pos_buffer(3 * (n_i - 1) + 1)
          meshes(m_i)%pos(2, n_i) = node_pos_buffer(3 * (n_i - 1) + 2)
          meshes(m_i)%pos(3, n_i) = node_pos_buffer(3 * (n_i - 1) + 3)
        end do

        do f_i = 1, num_faces_in_mesh(m_i)
          meshes(m_i)%neighbor_face_node(f_i, 1) = node_face_buffer(3 * (f_i - 1) + 1)
          meshes(m_i)%neighbor_face_node(f_i, 2) = node_face_buffer(3 * (f_i - 1) + 2)
          meshes(m_i)%neighbor_face_node(f_i, 3) = node_face_buffer(3 * (f_i - 1) + 3)
        end do

        call init_mesh_slave(meshes(m_i), m_i)
      end if
    end do

    ! Deallocate temporary meshes.
    deallocate(node_pos_buffer)
    deallocate(node_face_buffer)

    ! Report.
    call log_msg("  done", .false.)
  end subroutine init_meshes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialization of the particles
  !>
  !> The memory for the particles is allocated.
  !> The initial radius of the particles is set.

  subroutine init_particles()
    ! Report beginning of the initialization.
    call log_msg("allocating memory for particles...", .false.)

    ! Allocate arrays and initialize values.
    ! Allocate memory for particles in the subdomain.
    ! There are NUM_PART_LOC_MAX particles per subdomain allowed.
    allocate(particles(NUM_PART_LOC_MAX))
    allocate(part_ind(NUM_PART_LOC_MAX))
    particles(:)%currently_used = 0 ! set all particles to inactive
    part_ind(:) = -1 ! set all indices to invalid

    ! Set relative initial particle radius.
    ! If the particles shall be grown, their relative radius has already been taken from the parameter file,
    ! and nothing has to be done here.
    ! If the particles shall not be grown, their relative radius has to be set to unity.
    ! This is also the case if restoring from a checkpoint since checkpoints are always written after particle growth.
    if((mode_grow_particles .eqv. .false.) .or. ( is_restoring() .eqv. .true.)) then
      rel_initial_radius = 1
    end if

    call log_msg("  done", .false.)
  end subroutine init_particles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Positioning the particles
  !>
  !> \param[in] whole_N local chunk of the lattice with full halo of depth \c halo_extent
  !>
  !> The particles are positioned in the computational domain.
  !> This subroutine is only called when not restoring from a checkpoint because a checkpoint is only written
  !> after particles have been successfully positioned.
  !> Currently, there are the following options for particle positioning:
  !> 1) Position particles manually by using the initial positions from a file.
  !> 2) Position particless automatically and randomly with mutual overlap check (no check for overlap with walls!).
  !> The root process is responsible for initializing the particles and broadcasting relevant information.
  !> 1) Root checks if a new particle position overlaps with any of the already positioned particles.
  !>    NOTE: The overlap check is deactivated for manual particle positioning, so be careful then!
  !> 2) Once all positions are clear, the particle data is broadcast.
  !> 3) Each process creates the particles which are located in its physical domain.
  !> TODO: This subroutine should be revised.
  
  subroutine position_particles(whole_N)
    type(lbe_site),intent(in) :: whole_N(1-halo_extent:,1-halo_extent:,1-halo_extent:)

    ! Declare variables.
    type(init_particle_info), dimension(num_particles_gl) :: init_particle ! particle data to be sent by root to other processes
    type(init_particle_info) :: temp_particle ! temporary particle to store data
    integer, parameter :: init_file_unit = 13 ! input file unit
    integer :: error_stat ! switch for successful line reading
    integer :: counter_attempt ! counter for the number of attempts to position a particle without overlap
    integer :: c_i, c_j ! particle indices
    integer :: n_i ! node index
    integer :: ierror ! MPI error code
    integer :: rank_responsible ! rank responsible for a given particle
    integer :: num_particles_gl_check ! global number of particles after particle distribution
    integer :: stat
    real(kind=rk), dimension(3) :: dist_vec ! distance vector between particles
    real(kind=rk),allocatable,dimension(:,:,:) :: rock_state
    logical :: file_exist ! switch for file existence
    logical :: overlap ! switch for overlap detection
    logical :: do_check_rock_pos ! avoid random placement straight into rock?
    character(len=1000) :: line_buffer ! buffer for reading lines from file
    character(len=200) :: message ! message string

    ! Report beginning of the initialization.
    call log_msg("positioning particles...", .false.)

    ! only check for rock if desired by the user and if there is rock
    do_check_rock: if (rock_is_present().and.mode_check_rock_pos) then
       do_check_rock_pos = .true.
       if (myrankc==0) then
          allocate (rock_state&
               &(1-halo_extent:tnx+halo_extent&
               &,1-halo_extent:tny+halo_extent&
               &,1-halo_extent:tnz+halo_extent)&
               &,stat=stat)
          call check_allocate(stat,'position_particles(): rock_state')
       end if
       call gather_rock_state(whole_N,rock_state)
    else do_check_rock
       do_check_rock_pos = .false.
    end if do_check_rock

    ! Prepare particles.
    ! This is performed only by the root process.
    if(myrankc .eq. 0) then
      ! If manual particle positioning is specified by the user, the input file has to be checked for existence.
      if(mode_manual_pos .eqv. .true.) then
        ! Check whether position file exists or not.
        inquire(file="parm_positions.dat", exist=file_exist)

        ! If it does not exist, switch to automatic positioning and give a warning.
        if(file_exist .eqv. .false.) then
          mode_manual_pos = .false.
          call log_msg("  WARNING: input file parm_positions.dat not found, switching to automatic particle positioning", .false.)
        end if
      end if

      ! Open the particle position file if in manual particle positioning mode.
      if(mode_manual_pos .eqv. .true.) then
        ! Open file.
        open(init_file_unit, file="parm_positions.dat", status='old', action='read')

        ! Skip header line.
        read(init_file_unit, '(a999)') line_buffer
      end if

      ! Run over all particles.
      do c_i = 1, num_particles_gl
        ! Set the particle's global index, mesh_type, and radius
        temp_particle%particle_gl = c_i
        temp_particle%mesh_type = particle_data(c_i)%mesh_type
        temp_particle%radius = particle_data(c_i)%radius
        temp_particle%k_v = particle_data(c_i)%k_v
        temp_particle%k_at = particle_data(c_i)%k_at
        temp_particle%k_s = particle_data(c_i)%k_s
        temp_particle%k_al = particle_data(c_i)%k_al
        temp_particle%k_b = particle_data(c_i)%k_b
#ifdef IBM_DRAG
        temp_particle%force_const = particle_data(c_i)%force_const
#endif
#ifdef IBM_FIXED
        temp_particle%k_anchor = particle_data(c_i)%k_anchor
#endif

        ! If in manual mode, another line from the position file is read.
        if(mode_manual_pos .eqv. .true.) then
#ifdef IBM_DRAG
#ifdef IBM_FIXED
					call log_msg("ERROR: compiler flag IBM_DRAG and IBM_FIXED are not compatible!");
					call Abend
#endif 
					read(init_file_unit, fmt=*, iostat=error_stat) temp_particle%pos(:), temp_particle%angle, temp_particle%axis(:), temp_particle%force_const(:)
#else
					read(init_file_unit, fmt=*, iostat=error_stat) temp_particle%pos(:), temp_particle%angle, temp_particle%axis(:)
#endif
          ! If the data for the particle cannot be read from the file, give a warning and switch to automatic positioning.
          if(error_stat .ne. 0) then
            mode_manual_pos = .false.
            write(message, "('  WARNING: could not read line from file parm_positions.dat for particle ', i0 , ',' // &
            & ' switching to automatic particle positioning')") c_i
            call log_msg(trim(message), .false.)
            close(init_file_unit)
          end if
        end if

        ! Overwrite particle positions for swimmers.
#ifdef IBM_SWIMMER
        ! Perform some checks.
        if(swimmer_manualpos .eqv. .true.) then
          if(swimmer_anchor .eqv. .true.) then
            write(message, "('Swimmer anchor and manual positioning must not be switched on simultaneously.')")
            call error(trim(message))
          end if
          if(norm(swimmer_axis) .lt. 1.e-10) then
            write(message, "('Please use a non-zero swimmer alignment axis.')")
            call error(trim(message))
          end if
        end if
        
        ! Position beads.
        if(swimmer_manualpos .eqv. .false.) then
          if(c_i .eq. 1) then
            temp_particle%pos(1) = 0.5 * tnx - swimmer_d12
          end if

          if(c_i .eq. 2) then
            temp_particle%pos(1) = 0.5 * tnx
          end if
        
          if(c_i .eq. 3) then
            temp_particle%pos(1) = 0.5 * tnx + swimmer_d23
          end if
        
          temp_particle%pos(2) = 0.5 * tny
          temp_particle%pos(3) = 0.5 * tnz
        else
          if(c_i .eq. 1) then
            temp_particle%pos = swimmer_pos + swimmer_d12 * swimmer_axis / norm(swimmer_axis)
          end if

          if(c_i .eq. 2) then
            temp_particle%pos = swimmer_pos
          end if
          
          if(c_i .eq. 3) then
            temp_particle%pos = swimmer_pos - swimmer_d23 * swimmer_axis / norm(swimmer_axis)
          end if
        end if
#endif
        
        if(mode_manual_pos .eqv. .false.) then
          ! Obtain random numbers for particle orientation.
          call random_number(temp_particle%axis(:))
          temp_particle%axis(:) = (temp_particle%axis(:) - 0.5) * 2.0
          call random_number(temp_particle%angle)
          temp_particle%angle = temp_particle%angle * 360.0

          ! Position particle in such a way that no overlap occurs.
          ! Start with the assumption that an overlap exists.
          ! Count the number of attempts to position a particle.
          overlap = .true.
          counter_attempt = 0

          ! While the overlap condition is true, guess a new particle position and test again.
          try_positions: do while(overlap .eqv. .true.)
            ! Increase the counter and terminate the simulation if it becomes too large.
            counter_attempt = counter_attempt + 1

            if(counter_attempt .gt. 1000) then
              write(message, "('unsuccessful positioning of particle ', i0 , ' after 1000 attempts, volume fraction too large')") &
                & c_i
              call error(trim(message))
            end if

            ! Guess a new position.
            call random_number(temp_particle%pos(:))
            temp_particle%pos(1) = (temp_particle%pos(1) * (tnx - 2 * temp_particle%radius * rel_initial_radius)) &
                & + 0.5 + temp_particle%radius * rel_initial_radius
            temp_particle%pos(2) = (temp_particle%pos(2) * (tny - 2 * temp_particle%radius * rel_initial_radius)) &
                & + 0.5 + temp_particle%radius * rel_initial_radius
            temp_particle%pos(3) = (temp_particle%pos(3) * (tnz - 2 * temp_particle%radius * rel_initial_radius)) &
                & + 0.5 + temp_particle%radius * rel_initial_radius

            ! Assume that no overlap exists.
            overlap = .false.

            if (do_check_rock_pos) then
               ! reject positions that are completely surrounded by rock sites
               if (inside_rock_global(rock_state,temp_particle%pos)) then
                  overlap = .true.
                  cycle try_positions
               end if
            end if

            c_j = 0

            ! Check for overlap with all other particles already positioned.
            do while(c_j + 1 .lt. c_i)
              c_j = c_j + 1

              ! Compute distance vector between the particles.
              dist_vec = distance_vector(temp_particle%pos, init_particle(c_j)%pos)

              ! Check the distance and compare with the radii of the particles.
              ! Use a safety distance of 1 lattice node.
              if(norm(dist_vec) .lt. temp_particle%radius * rel_initial_radius &
                & + init_particle(c_j)%radius * rel_initial_radius + 1.0) then
                overlap = .true.
                exit
              end if
            end do
         end do try_positions
        end if

        ! Add the temporary particle to the particle initialization list.
        init_particle(c_i) = temp_particle
      end do

      ! Close the particle position file.
      if(mode_manual_pos .eqv. .true.) then
        close(init_file_unit)
      end if
    end if

    ! Broadcast particle list to all other processes.
    call MPI_Bcast(init_particle, num_particles_gl, MPI_PARTICLE_INIT, 0, comm_cart, ierror)

    ! Each rank runs over this list and takes only the particles it is responsible for.
    do c_i = 1, num_particles_gl
      ! If the rank is responsible, the particle is locally created and initialized.
      rank_responsible = rank_responsible_for_pos(init_particle(c_i)%pos)

      if(rank_responsible .eq. myrankc) then
        call add_particle(meshes(init_particle(c_i)%mesh_type))
        particles(part_ind(num_particles_loc))%particle_index_gl = c_i
        call set_radius(particles(part_ind(num_particles_loc)), init_particle(c_i)%radius, .true.)
        call set_center_position(particles(part_ind(num_particles_loc)), init_particle(c_i)%pos)
        call rotate_particle(particles(part_ind(num_particles_loc)), init_particle(c_i)%angle, init_particle(c_i)%axis)
        particles(part_ind(num_particles_loc))%k_v = init_particle(c_i)%k_v
        particles(part_ind(num_particles_loc))%k_at = init_particle(c_i)%k_at
        particles(part_ind(num_particles_loc))%k_s = init_particle(c_i)%k_s
        particles(part_ind(num_particles_loc))%k_al = init_particle(c_i)%k_al
        particles(part_ind(num_particles_loc))%k_b = init_particle(c_i)%k_b
#ifdef IBM_DRAG
        particles(part_ind(num_particles_loc))%force_const = init_particle(c_i)%force_const
#endif
#ifdef IBM_FIXED
				particles(part_ind(num_particles_loc))%k_anchor = init_particle(c_i)%k_anchor
#endif
      end if

      ! If no rank is responsible for a particle, the simulation has to be terminated.
      if(rank_responsible == -1) then
        write(message, "('rank ', i0, ' found particle ', i0, ' without responsible rank')") myrankc, c_i
        call error(trim(message))
      end if
    end do

    ! Check global number of particles.
    ! The sum of all particles in each process must be the same as the global number of particles before distribution.
    call MPI_Allreduce(num_particles_loc, num_particles_gl_check, 1, MPI_INTEGER, MPI_SUM, comm_cart, ierror)
    write(message, "('  consistency check: distributed ', i0, ' particles')") num_particles_gl_check
    call log_msg(trim(message), .false.)

    if(num_particles_gl_check .ne. num_particles_gl) then
      write(message, "('rank ', i0, ' reports ', i0, ' particles rather than ', i0)") myrankc, num_particles_gl_check, &
        & num_particles_gl
      call error(trim(message))
    end if
#ifdef IBM_FIXED
    do c_i = 1, num_particles_loc
      do n_i = 1, particles(part_ind(c_i))%num_nodes
				particles(part_ind(c_i))%node(n_i)%pos_anchor(:) = particles(part_ind(c_i))%node(n_i)%pos(:)
			end do
		end do
#endif
    ! Check sanity of particles.
    ! If any particle node is not defined, the simulation has to be terminated.
    do c_i = 1, num_particles_loc
      do n_i = 1, particles(part_ind(c_i))%num_nodes
        if(check_nan(particles(part_ind(c_i))%node(n_i)%pos(1)) .or. &
          & check_nan(particles(part_ind(c_i))%node(n_i)%pos(2)) .or. &
          & check_nan(particles(part_ind(c_i))%node(n_i)%pos(3))) then
          write(message, "('rank ', i0, ' detected NaN in position of particle ', i0)") myrankc, c_i
          call error(trim(message))
        end if
      end do
    end do

    ! Update areas.
    do c_i = 1, num_particles_loc
      call update_face_areas(particles(part_ind(c_i)))
      call update_node_areas(particles(part_ind(c_i)))
#ifdef IBM_INDEXFIELD
      call update_node_normals(particles(part_ind(c_i)))
#endif
    end do

    ! Deallocate arrays for mesh and particle initialization.
    ! They are not required anymore.
    if(myrankc .eq. 0) then
      deallocate(mesh_data)
      deallocate(particle_data)
      if (do_check_rock_pos) deallocate (rock_state)
    end if

    ! Report success of the initialization.
    call log_msg("  done", .false.)
  end subroutine position_particles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Restoring the particles from a checkpoint
  !>
  !> The particles are restored from a checkpoint which has been written in binary format (XDRF).
	!> TODO Add anchor spring constant and initial position to this

  subroutine restore_particles()
    ! Declare variables.
    integer :: file_unit ! input file unit
    integer :: ierror ! rror code
    integer :: c_i ! particle index
    integer :: n_i ! node index
    integer :: temp_num_particles_loc, temp_mesh_type ! temporary variables for reading checkpoint files
    character(len=100) :: filename ! output filename

    ! Set variables.
    file_unit = 42

    ! Report beginning of the initialization.
    call MPI_Barrier(comm_cart, ierror) ! just to make sure that all processes have reached this point.
    call log_msg("restoring particles...", .false.)

    ! Create and open checkpoint files.
    call lbe_make_filename_restore_rank(filename, 'checkpoint_IBM', '.xdr', myrankc)
    call xdrfopen(file_unit, filename, "r", ierror)
    call check_xdrfopen(ierror, filename)

    ! Read number of local particles.
    call xdrfint(file_unit, temp_num_particles_loc, ierror)

    ! Read all data for each local particle.
    do c_i = 1, temp_num_particles_loc
      ! Read particle mesh type and create new particle.
      call xdrfint(file_unit, temp_mesh_type, ierror)
      call add_particle(meshes(temp_mesh_type))

      ! Read particle data.
      call xdrfint(file_unit, particles(part_ind(c_i))%particle_index_gl, ierror)
      call xdrfint(file_unit, particles(part_ind(c_i))%num_jumps(1), ierror)
      call xdrfint(file_unit, particles(part_ind(c_i))%num_jumps(2), ierror)
      call xdrfint(file_unit, particles(part_ind(c_i))%num_jumps(3), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%radius_0, ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%radius, ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%volume_0, ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%volume, ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%surface_0, ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%surface, ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%k_s, ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%k_al, ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%k_b, ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%k_v, ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%k_at, ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%center(1), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%center(2), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%center(3), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%center_old(1), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%center_old(2), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%center_old(3), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%force_total(1), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%force_total(2), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%force_total(3), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%torque_total(1), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%torque_total(2), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%torque_total(3), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%linear_velocity(1), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%linear_velocity(2), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%linear_velocity(3), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%linear_momentum(1), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%linear_momentum(2), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%linear_momentum(3), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%angular_velocity(1), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%angular_velocity(2), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%angular_velocity(3), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%angular_momentum(1), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%angular_momentum(2), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%angular_momentum(3), ierror)
#ifdef IBM_FIXED
      call xdrfdouble(file_unit, particles(part_ind(c_i))%k_anchor, ierror)
#endif

      ! Read node data.
      do n_i = 1, particles(part_ind(c_i))%num_nodes
        call xdrfdouble(file_unit, particles(part_ind(c_i))%node(n_i)%pos(1), ierror)
        call xdrfdouble(file_unit, particles(part_ind(c_i))%node(n_i)%pos(2), ierror)
        call xdrfdouble(file_unit, particles(part_ind(c_i))%node(n_i)%pos(3), ierror)
        call xdrfdouble(file_unit, particles(part_ind(c_i))%node(n_i)%pos_old(1), ierror)
        call xdrfdouble(file_unit, particles(part_ind(c_i))%node(n_i)%pos_old(2), ierror)
        call xdrfdouble(file_unit, particles(part_ind(c_i))%node(n_i)%pos_old(3), ierror)
        call xdrfdouble(file_unit, particles(part_ind(c_i))%node(n_i)%vel(1), ierror)
        call xdrfdouble(file_unit, particles(part_ind(c_i))%node(n_i)%vel(2), ierror)
        call xdrfdouble(file_unit, particles(part_ind(c_i))%node(n_i)%vel(3), ierror)
#ifdef IBM_FIXED
        	call xdrfdouble(file_unit, particles(part_ind(c_i))%node(n_i)%pos_anchor(1), ierror)
        	call xdrfdouble(file_unit, particles(part_ind(c_i))%node(n_i)%pos_anchor(2), ierror)
        	call xdrfdouble(file_unit, particles(part_ind(c_i))%node(n_i)%pos_anchor(3), ierror)
#endif
      end do

      ! Compute dependent data.
      call update_face_areas(particles(part_ind(c_i)))
      call update_node_areas(particles(part_ind(c_i)))
#ifdef IBM_INDEXFIELD
      call update_node_normals(particles(part_ind(c_i)))
#endif
    end do

    ! Close the checkpoint file.
    call xdrfclose(file_unit, ierror)

    ! Report success of the initialization.
    call MPI_Barrier(comm_cart, ierror) ! just to make sure that all processes have reached this point.
    call log_msg("all particles restored", .false.)
  end subroutine restore_particles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computation of volume fraction of the particles
  !>
  !> The volume fraction in the simulation is computed and reported.
  !> It is defined as the total volume covered by particles divided by the total number of fluid lattice nodes.

  subroutine compute_volume_fraction()
    ! Declare variables
    integer :: num_fluid_nodes_gl ! global number of fluid nodes
    integer :: c_i ! particle index
    integer :: ierror ! MPI error code
    real(kind=rk) :: volume_loc ! volume filled by particles in subdomain
    real(kind=rk) :: volume_gl ! volume filled by particles in total domain
    real(kind=rk) :: volume_fraction ! volume fraction of particles
    character(len=200) :: message ! message

    ! Compute global number of fluid nodes.
    ! Here it is assumed that the number of fluid nodes (n_sites_fluid) has already been updated.
    num_fluid_nodes_gl = n_sites_fluid

    ! Compute local particle volume.
    volume_loc = 0

    do c_i = 1, num_particles_loc
      volume_loc = volume_loc + particles(part_ind(c_i))%volume_0
    end do

    ! Compute total particle volume.
    call MPI_Allreduce(volume_loc, volume_gl, 1, MPI_REAL8, MPI_SUM, comm_cart, ierror)

    ! Compute volume fraction.
    volume_fraction = volume_gl / num_fluid_nodes_gl

    ! Report result.
    write(message, "('particle volume fraction: ', F0.3, '%')") volume_fraction * 100d0
    call log_msg(trim(message), .false.)
  end subroutine compute_volume_fraction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Growth of particles to their full size
  !>
  !> If the particles have been positioned with a decreased size, they have to be grown to their full size again.
  !> The particles are successively inflated using simplified dynamics to reach the desired volume fraction.
  !> If particle growth is deactivated or no particles are used, this subroutine has no effect.
  !> If the simulation is restarted from a checkpoint, the particle growth is deactivated.
  !> The reason is that no checkpoints during particle growth are written.
  !> It is always expected that the system is within the time loop whenever a checkpoint is written.

  subroutine grow_particles(N)
    type(lbe_site), dimension(0:, 0:, 0:), intent(in) :: N !< lattice

    ! Declare variables.
    integer :: time ! time index for growing particles
    integer :: c_i ! particle index
    integer :: n_i ! node index
    real(kind=rk), parameter :: mass = 20.0 ! artificial mass for nodes and MD update.
    real(kind=rk) :: radius_new ! new particle radius
    real(kind=rk) :: vel_mag ! maximum node velocity
    character(len=200) :: message ! message string

    ! Check execution conditions.
    ! If growth is not required, skip this routine.
    if((mode_grow_particles .eqv. .false.) .or. (num_meshes .eq. 0) .or. ( is_restoring() .eqv. .true.)) then
      call log_msg("skipping particle growth (deactivated or restarting from checkpoint)", .false.)
      return
    end if

    ! Report start of growing.
    call log_msg("growing particles...", .false.)

    ! Execute growth loop.
    do time = 1, num_steps_growth
      ! Report progress.
      write(message, "('  ', F0.3, '% complete')") (time * 100.0d0) / num_steps_growth
      call log_msg(trim(message), .false.)

      ! Inflate particles and computen internal forces.
      do c_i = 1, num_particles_loc
        ! Inflate particles.
        radius_new = particles(part_ind(c_i))%radius_0 * (rel_initial_radius**3.0 + &
                      & ((1.0 - (rel_initial_radius)**3.0) * time) / num_steps_growth )**(1.0 / 3.0)
        call set_radius(particles(part_ind(c_i)), radius_new, .false.)

        ! Compute internal forces.
        call forces_compute(particles(part_ind(c_i)), surften, contactangle, width, friction, 0, k_b = growth_bending_modulus)
      end do

      ! Distribute nodes to neighboring processes and update node look-up table.
      call distribute_nodes()
      call update_node_interaction_LUT()

      ! Compute external forces and combine forces.
!      call compute_particle_forces_ext(N) ! currently not needed
      call compute_particle_particle_forces()
      call compute_wall_gradient_forces()
      call combine_forces()

      ! Collect nodes from neighboring processes.
      call collect_nodes()

      ! Compute friction forces.
      ! These forces are used to damp the kinetic motion of the particles which is constantly pumped into the system.
      ! Otherwise, particle deformation becomes irratic and the simulation unstable.
      do c_i = 1, num_particles_loc
        call forces_friction(particles(part_ind(c_i)), friction)
      end do

      ! Update node velocities and positions.
      ! A simplifies MD approach is used for time integration.
      do c_i = 1, num_particles_loc
        do n_i = 1, particles(part_ind(c_i))%num_nodes
          particles(part_ind(c_i))%node(n_i)%vel = particles(part_ind(c_i))%node(n_i)%vel &
            & + particles(part_ind(c_i))%node(n_i)%force_tot / mass

          ! Limit maximum node velocity to 0.1 in order to avoid strong deformations in a short time.
          vel_mag = norm(particles(part_ind(c_i))%node(n_i)%vel)
          if(vel_mag > 0.1) then
            particles(part_ind(c_i))%node(n_i)%vel = particles(part_ind(c_i))%node(n_i)%vel * 0.1 / vel_mag
          end if
        end do
      end do

#ifdef IBM_FIXED
			! Update anchor positions of nodes
			do c_i = 1, num_particles_loc
				do n_i = 1, particles(part_ind(c_i))%num_nodes
					particles(part_ind(c_i))%node(n_i)%pos_anchor = particles(part_ind(c_i))%node(n_i)%pos
				end do
			end do
#endif

      ! Update particles and exchange with neighboring processes.
      call update_particles_pre()
      call exchange_particles()
      call update_particles_post()
      call write_particles_vtk("vtk_lagrange_grow", time)
    end do

    ! Report end of growth process.
    call log_msg("  growth process complete", .false.)
  end subroutine grow_particles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Manipulate particles
  !>
  !> The developer may add some operations to manipulate the particles.
  !> This could be deformation, applying forces on single nodes etc.
  !> The idea is to make the initialization more flexible.
  !> This subroutine is currently only used for benchmarking the particle code.
  !> The subroutine is skipped if restarting from a checkpoint since the idea is
  !> to manipulate the particles only directly after initialization.

  subroutine manipulate_particles()
		! Declare variables.
 	
		! Check execution conditions.
    if(is_restoring() .eqv. .true.) then
      call log_msg("skipping particle manipulation (restarting from checkpoint)", .false.)
      return
    end if

    ! Nothing is happening right now.

  end subroutine manipulate_particles

! endif IBM_PART
#endif

end module lsuperstruct_init_module
