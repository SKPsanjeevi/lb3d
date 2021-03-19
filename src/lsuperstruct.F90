!> Lagrangian superstructure module
!>
!> This module contains all interface subroutines visible from \c lbe.F90
!> for the initialization, timeloop, and dumping related to the Lagrangian particles.
!> All subroutines are public, and this module should only be included in \c lbe.F90.
!> The task of this module is to relay the requests from \c lbe.F90 to the corresponding
!> modules for the Lagrangian particles and the IBM plug-in.
!> The runtime for each subroutine call is measured for time statistics.

#include "lbe.h"

module lsuperstruct_module

#ifdef IBM_PART

  ! Include external modules.
  use lbe_globals_module, only : halo_extent, ti_IBM_init, ti_IBM_forces, ti_IBM_intspread, ti_IBM_MPI, ti_IBM_update, ti_IBM_dump
  use lbe_parms_module, only : nt
  use lbe_types_module, only : lbe_site
  use lbe_timer_module
  use lsuperstruct_IBM_module
  use lsuperstruct_init_module
  use lsuperstruct_interface_module
  use lsuperstruct_timeloop_module
  use lsuperstruct_parallel_module
  use lsuperstruct_dump_module
  use lextobj_module

  implicit none

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read IBM input from parameter files

  subroutine read_ibm_input()
    call start_timer(ti_IBM_init)
    call init_read_parameters() ! => lsuperstruct_init_module
    call check_parameter_conflicts() ! => lsuperstruct_init_module
    call stop_timer(ti_IBM_init)
  end subroutine read_ibm_input

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize particle modules
  !>
  !> \param[in,out] lbe_N local chunk of the lattice with halo extent 1 (old LB3D style)
  !> \param[in] whole_N local chunk of the lattice with full halo of depth \c halo_extent

  subroutine lagr_init(lbe_N, whole_N)
    type(lbe_site), dimension(0:, 0:, 0:), intent(inout) :: lbe_N !< lattice
    type(lbe_site), intent(in) :: whole_N(1-halo_extent:,1-halo_extent:,1-halo_extent:)

    call start_timer(ti_IBM_init)
    call init_lsuperstruct(lbe_N, whole_N) ! => lsuperstruct_init_module
    call stop_timer(ti_IBM_init)
  end subroutine lagr_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute internal particle forces

  subroutine lagr_compute_forces_int(time)
    integer, intent(in) :: time !< simulation time
  
    call start_timer(ti_IBM_forces)
    call compute_particle_forces_int(time) ! => lsuperstruct_timeloop_module
    call stop_timer(ti_IBM_forces)
  end subroutine lagr_compute_forces_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Distribute nodes to neighboring processes (MPI)

  subroutine lagr_distribute_nodes()
    call start_timer(ti_IBM_MPI)
    call distribute_nodes() ! => lsuperstruct_parallel_module
    call stop_timer(ti_IBM_MPI)
  end subroutine lagr_distribute_nodes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update node look-up tables

  subroutine lagr_update_lut()
    call start_timer(ti_IBM_update)
#ifndef NOFLUIDDYNAMICS
    call update_node_interaction_LUT() ! => lsuperstruct_timeloop_module
#endif
    call stop_timer(ti_IBM_update)
  end subroutine lagr_update_lut

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute external particle forces

  subroutine lagr_compute_forces_ext(N)
    type(lbe_site), dimension(0:, 0:, 0:), intent(in) :: N !< lattice

    call start_timer(ti_IBM_forces)
    call compute_particle_forces_ext(N) ! => lsuperstruct_timeloop_module
    call stop_timer(ti_IBM_forces)
  end subroutine lagr_compute_forces_ext

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Add up internal and external force contributions

  subroutine lagr_combine_forces()
    call start_timer(ti_IBM_forces)
    call combine_forces() ! => lsuperstruct_timeloop_module
    call stop_timer(ti_IBM_forces)
  end subroutine lagr_combine_forces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Spread forces from the Lagrangian mesh to the Eulerian lattice

  subroutine lagr_spread_forces()
    call start_timer(ti_IBM_intspread)
#ifndef NOFLUIDDYNAMICS
    call IBM_spread_forces(.true.) ! => lsuperstruct_IBM_module
#endif
    call stop_timer(ti_IBM_intspread)
  end subroutine lagr_spread_forces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Rebuild interior/exterior index field
  !>
  !> This functionality is only enabled in IBM_INDEXFIELD mode.
  
#ifdef IBM_INDEXFIELD
  subroutine lagr_rebuild_index(N, force_rebuild)
    type(lbe_site), dimension(0:, 0:, 0:), intent(inout) :: N !< lattice
    logical, intent(in) :: force_rebuild !< whether to force a complete rebuild

    call start_timer(ti_IBM_intspread)
    call rebuild_interior_index(N, force_rebuild) ! => lsuperstruct_timeloop_module
    call stop_timer(ti_IBM_intspread)
  end subroutine lagr_rebuild_index
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Exchange index field halo

#ifdef IBM_INDEXFIELD
  subroutine lagr_exchange_index_halo()
    call start_timer(ti_IBM_MPI)
    call exchange_index_halo() ! => lsuperstruct_parallel_module
    call stop_timer(ti_IBM_MPI)
  end subroutine lagr_exchange_index_halo
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute selective wall forces (repulsion of external fluid by walls)
  !>
  !> This functionality is only enabled in IBM_INDEXFIELD mode.
  !> TODO: This is a highly experimental functionality which may not survive in the near future.

#ifdef IBM_INDEXFIELD
  subroutine lagr_selective_wall_forces()
    call start_timer(ti_IBM_intspread)
    call compute_selective_wall_forces() ! => lsuperstruct_interface_module
    call stop_timer(ti_IBM_intspread)
  end subroutine lagr_selective_wall_forces
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Recolour the interior region of the deformable particles.
  !>
  !> This functionality is only enabled in IBM_INDEXFIELD and NOSURFACTANT mode.
  !> TODO: This is a highly experimental functionality which may not survive in the near future.

#if IBM_BINARYIBM
   subroutine lagr_recolour_interior(N)
     type(lbe_site), dimension(0:, 0:, 0:), intent(inout) :: N !< lattice
 
     call start_timer(ti_IBM_intspread)
     call recolour_interior_sites(N) ! => lsuperstruct_interface_module
     call stop_timer(ti_IBM_intspread)
   end subroutine lagr_recolour_interior
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute multicomponent forces (extended no-slip force at particle surface)
  !>
  !> This functionality is only enabled in IBM_INDEXFIELD and NOSURFACTANT mode.
  !> TODO: This is a highly experimental functionality which may not survive in the near future.

#ifdef IBM_BINARYIBM
   subroutine lagr_multicomponent_forces(N)
     type(lbe_site), dimension(0:, 0:, 0:), intent(in) :: N !< lattice
 
     call start_timer(ti_IBM_intspread)
     call compute_multicomponent_forces(N) ! => lsuperstruct_interface_module
     call stop_timer(ti_IBM_intspread)
   end subroutine lagr_multicomponent_forces
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute physical velocity on lattice

  subroutine lagr_compute_vel_phys(N)
    type(lbe_site), dimension(0:, 0:, 0:), intent(in) :: N !< lattice

    call start_timer(ti_IBM_update)
#ifndef NOFLUIDDYNAMICS
    call compute_fluid_velocity(N) ! => lsuperstruct_interface_module
#endif
    call stop_timer(ti_IBM_update)
  end subroutine lagr_compute_vel_phys

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Exchange physical velocity halo

  subroutine lagr_exchange_velocity_halo()
    call start_timer(ti_IBM_MPI)
#ifndef NOFLUIDDYNAMICS
    call exchange_velocity_halo() ! => lsuperstruct_parallel_module
#endif
    call stop_timer(ti_IBM_MPI)
  end subroutine lagr_exchange_velocity_halo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Interpolate velocities from the Eulerian lattice to the Lagrangian mesh

  subroutine lagr_interpolate_velocities()
    call start_timer(ti_IBM_intspread)
#ifndef NOFLUIDDYNAMICS
    call IBM_interpolate_velocities() ! => lsuperstruct_IBM_module
#endif
    call stop_timer(ti_IBM_intspread)
  end subroutine lagr_interpolate_velocities

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Collect nodes from neighboring processes (MPI)

  subroutine lagr_collect_nodes()
    call start_timer(ti_IBM_MPI)
    call collect_nodes() ! => lsuperstruct_parallel_module
    call stop_timer(ti_IBM_MPI)
  end subroutine lagr_collect_nodes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update particles (part 1)

  subroutine lagr_update_particles_pre()
    call start_timer(ti_IBM_update)
    call update_particles_pre() ! => lsuperstruct_timeloop_module
    call stop_timer(ti_IBM_update)
  end subroutine lagr_update_particles_pre

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Exchange particles with neighboring processes (MPI)

  subroutine lagr_exchange_particles()
    call start_timer(ti_IBM_MPI)
    call exchange_particles() ! => lsuperstruct_parallel_module
    call stop_timer(ti_IBM_MPI)
  end subroutine lagr_exchange_particles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update particles (part 2)

  subroutine lagr_update_particles_post()
    call start_timer(ti_IBM_update)
    call update_particles_post() ! => lsuperstruct_timeloop_module
    call stop_timer(ti_IBM_update)
  end subroutine lagr_update_particles_post

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Dump all desired data

  subroutine lagr_dump_data(N)
    type(lbe_site), dimension(0:, 0:, 0:), intent(in) :: N !< lattice

    call start_timer(ti_IBM_dump)
    call write_particles_vtk("vtk_lagrange", nt) ! => lsuperstruct_dump_module
    call dump_lattice_profiles(N) ! => lsuperstruct_dump_module
    call dump_static_particle_data() ! => lsuperstruct_dump_module
    call write_particles_dat() ! => lsuperstruct_dump_module
#ifdef IBM_INDEXFIELD
    call dump_indexfield_lattice() ! => lsuperstruct_dump_module
#endif
    call stop_timer(ti_IBM_dump)
  end subroutine lagr_dump_data

! endif IBM_PART
#endif

end module lsuperstruct_module
