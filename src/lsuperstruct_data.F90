!> Lagrangian superstructure data module
!>
!> This module contains all data required for the Lagrangian superstructure.
!> It should be included in all submodules of \c lsuperstruct.F90 but not in \c lsuperstruct.F90 itself.

#include "lbe.h"

module lsuperstruct_data_module

#ifdef IBM_PART

  ! Include external modules.
  use lbe_globals_module, only : rk

  implicit none

  ! Hard-coded constants
  ! Be careful changing these values as this may deside about the fate of sanity.
  ! TODO:
  ! - Some of these should not be hard-coded in the future.
  integer, parameter :: NUM_PART_LOC_MAX = 300 !< maximum number of particles allowed in each subdomain
  integer, parameter :: MAX_NUM_NODES_SEND = 10000 !< maximum number of nodes which can be sent in one direction
  integer, parameter :: MAX_NUM_NODES = 300000 !< maximum number of nodes which can be located in one subdomain
  integer, parameter :: MAX_NODES_SUBBOX = 25 !< maximum number of nodes which can be located in a LUT subbox
  integer, parameter :: HALO_WIDTH = 2 !< halo width for sending particle nodes to neighboring processes
  real(kind=rk), parameter :: INT_RANGE_PART_PART = 1 !< range of particle node-particle node interaction force
  real(kind=rk), parameter :: INT_RANGE_WALL_PART = 1.5 !< range of particle node-wall node interaction force
  real(kind=rk), parameter :: DIST_CUTOFF = 0.5 !< cutoff distance for defining interior and exterior lattice nodes

  !> Definition of the superstructure node structure
  !>
  !> It contains all data for the nodes residing in the superstructure.

  type superstruct_node
    integer :: particle_gl !< global parent particle
    integer :: particle_loc !< global parent particle
    integer :: node !< node index in parent particle
    integer :: physical !< switch if node is in physical domain or not
    integer, dimension(3) :: send !< indicator for sending along the three axes
    integer, dimension(3) :: halo !< indicator for halo along the three axes
    real(kind = rk), dimension(3) :: pos !< physical node position
    real(kind = rk), dimension(3) :: vel !< physical node velocity
    real(kind = rk), dimension(3) :: force_tot !< total force on node (except for interaction forces)
    real(kind = rk), dimension(3) :: force_int !< interaction force on node
#ifdef IBM_FIXED
		real(kind = rk), dimension(3) :: pos_anchor !< initial physical node position 
#endif
    real(kind = rk), dimension(3) :: normal !< node normal vector
    real(kind = rk) :: curv_radius !< curvature radius
!#ifdef IBM_FIXED
		!real(kind = rk), dimension(3) :: force_anchor !< initial physical node position 
!#endif
  end type superstruct_node

  !> Definition of the particle information structure
  !>
  !> It contains all data for the particles to be set up during initialization.

  type init_particle_info
    integer :: particle_gl !< global particle index
    integer :: mesh_type !< mesh type
    real(kind=rk) :: radius !< particle radius
    real(kind=rk) :: angle !< rotation angle
    real(kind=rk) :: k_v !< rotation angle
    real(kind=rk) :: k_at !< rotation angle
    real(kind=rk) :: k_s !< rotation angle
    real(kind=rk) :: k_al !< rotation angle
    real(kind=rk) :: k_b !< rotation angle
#ifdef IBM_DRAG
    real(kind=rk), dimension(3) :: force_const !< constant force to be applied to particle
#endif
#ifdef IBM_FIXED
    real(kind=rk) :: k_anchor !< spring constant for anchor
#endif
    real(kind=rk), dimension(3) :: pos !< particle center position
    real(kind=rk), dimension(3) :: axis !< rotation axis
  end type init_particle_info

  !> Definition of the mesh information structure
  !>
  !> It contains all data for the meshes to be set up during initialization.

  type init_mesh_info
    integer :: num_particles
    real(kind=rk) :: radius
    real(kind=rk) :: k_v
    real(kind=rk) :: k_at
    real(kind=rk) :: k_s
    real(kind=rk) :: k_al
    real(kind=rk) :: k_b
#ifdef IBM_FIXED
		real(kind=rk) :: k_anchor
#endif
    character(len=50) :: filename
  end type init_mesh_info

  !> Definition of the particle information structure
  !>
  !> It contains all data for the particles to be positioned during initialization.

  type exchange_particle_info
    integer :: num_particles !< number of particles to send
    integer, dimension(NUM_PART_LOC_MAX) :: index_gl !< global particle index
    integer, dimension(NUM_PART_LOC_MAX) :: index_loc !< local particle index
    integer, dimension(NUM_PART_LOC_MAX) :: mesh_type !< mesh type
  end type exchange_particle_info

  !> Definition of the static particle data structure
  !>
  !> It contains all particle data required for dumping.

  type dump_static_particle_info
    integer :: particle_index_gl !< global particle index (unique for each particle, also in parallel simulation)
    real(kind=rk) :: volume !< current volume
    real(kind=rk) :: surface !< current surface
    real(kind=rk) :: erg_tot !< total energy
    real(kind=rk) :: erg_s !< strain energy
    real(kind=rk) :: erg_b !< bending energy
    real(kind=rk) :: erg_v !< volume energy
    real(kind=rk) :: erg_at !< surface energy
    real(kind=rk) :: erg_int !< interaction energy
		!TODO add anchor energy
    real(kind=rk) :: stresslet_xx !< particle stresslet (xx-component)
    real(kind=rk) :: stresslet_xy !< particle stresslet (xy-component)
    real(kind=rk) :: stresslet_xz !< particle stresslet (xz-component)
    real(kind=rk) :: stresslet_yy !< particle stresslet (yy-component)
    real(kind=rk) :: stresslet_yz !< particle stresslet (yz-component)
    real(kind=rk) :: stresslet_zz !< particle stresslet (zz-component)
    real(kind=rk), dimension(3) :: center !< current particle center position
    real(kind=rk), dimension(3) :: force_total !< total force acting on particle
    real(kind=rk), dimension(3) :: force_anchor !< total force acting on particle
    real(kind=rk), dimension(3) :: torque_total !< total torque acting on particle
    real(kind=rk), dimension(3) :: linear_velocity !< linear particle velocity
    real(kind=rk), dimension(3) :: linear_momentum !< linear particle momentum
    real(kind=rk), dimension(3) :: angular_velocity !< angular particle velocity
    real(kind=rk), dimension(3) :: angular_momentum !< angular particle momentum
    real(kind=rk), dimension(3) :: iner_ell_axes !< semiaxes of inertia ellipsoid
    real(kind=rk), dimension(3, 3) :: iner_ell_vecs !< orientations of inertia ellipsoid
  end type dump_static_particle_info

  ! Variables contained in the namelist and the parameter files
  ! These variables do not have to be saved in a checkpoint.
  ! All values given here are standard values which are replaced by the values in the parameter files.
  ! TODO: Check this list again after cleaning the code.
  integer :: IBM_RANGE = 2 !< range of the IBM interpolations
  logical :: mode_check_rock_pos = .false. !< avoid putting particles into rocks; costs a lot of memory on root => disabled on default
  logical :: mode_polydispersity = .false. !< switch for polydispersity mode
  logical :: mode_manual_pos = .false. !< switch for manual particle positioning
  logical :: mode_grow_particles = .true. !< switch for particle growth
  real(kind=rk) :: rel_initial_radius = 5.0d-1 !< initial relative radius of the particles
  integer :: num_steps_growth = 3000 !< number of time steps to grow particles
  real(kind=rk) :: growth_bending_modulus = -1.0_rk !< bending modulus for growing (negative value: use particle's own k_b)
  logical :: dump_vtk = .false. !< whether to dump particle VTK files
  integer :: time_step_dump_vtk = 100 !< number of time steps after which VTK files for the particles are dumped
  logical :: dump_profiles = .false. !< whether to dump lattice profiles
  integer :: time_step_dump_profiles = 100 !< number of time steps after which lattice profiles are dumped
  logical :: dump_particles = .false. !< whether to dump particle header data
  integer :: time_step_dump_particles = 100 !< number of time steps after which particle header data are dumped
  logical :: dump_indexfield = .false. !< whether to dump interior/exterior index field
  integer :: time_step_dump_indexfield = 100 !< number of time steps after which interior/exterior index field data are dumped
  integer :: num_meshes = 0 !< number of meshes
  real(kind=rk) :: int_strength_part_part = 1.0d-2 !< interaction strength (particle-particle)
  real(kind=rk) :: int_strength_part_wall = 1.0d-2 !< interaction strength (particle-wall)
  real(kind=rk) :: int_strength_gradient = 1.0d-2 !< interaction strength (rock gradient)
  real(kind=rk) :: friction = 0.75 !< friction coefficient for reducing maximum node velocity
  real(kind=rk) :: IBM_massdensity = 1.0 !< mass area density of nodes
  real(kind=rk) :: viscosity_contrast = 1.0d0 !< viscosity contrast between interior and exterior fluids
  real(kind=rk) :: selective_wall_force_magnitude = 0.0d0 !< force magnitude for selective particle-wall interaction
  real(kind=rk) :: selective_wall_force_range = 0.0d0 !< force magnitude for selective particle-wall interaction
  real(kind=rk) :: ibm_colour = 0.0d0 !< color of IBM particles
  real(kind=rk) :: g_ri = 0.0d0 !< red-index field interaction strength
  real(kind=rk) :: g_bi = 0.0d0 !< red-index field interaction strength
  real(kind=rk) :: noslipcoeff = 0.0d0 !< friction coefficient to enforce multicomponent no-slip
  real(kind=rk) :: surften = 0.0d0 !< surface tension of emulated interface
  real(kind=rk) :: contactangle = 90 !< contact angle (degrees)
  real(kind=rk) :: width = 0.5d0 !< width of emulated interface
  real(kind=rk) :: penalty_red_blue = 0.0d0 !< penalty modulus for elastic red-blue force
  integer :: time_step_rebuild_index = 10000 !< time after which the interior index is completely rebuilt
#ifdef IBM_SWIMMER
  real(kind=rk) :: swimmer_d12 = 20. !< distance between swimmer 1 and 2
  real(kind=rk) :: swimmer_d23 = 20. !< distance between swimmer 2 and 3
  real(kind=rk) :: swimmer_modulus12 = 0. !< distance between swimmer 1 and 2
  real(kind=rk) :: swimmer_modulus23 = 0. !< distance between swimmer 2 and 3
  real(kind=rk) :: swimmer_kb = 0. !< bending modulus between swimmer springs
  real(kind=rk) :: swimmer_ka = 0. !< anchor modulus
  logical :: swimmer_anchor = .false. !< use anchor force
  real(kind=rk) :: swimmer_force_mag1 = 0. !< force magnitude of swimmer 1
  real(kind=rk) :: swimmer_force_mag2 = 0. !< force magnitude of swimmer 2
  real(kind=rk) :: swimmer_force_mag3 = 0. !< force magnitude of swimmer 3
  real(kind=rk) :: swimmer_period1 = 0. !< oscillation period of swimmer 1
  real(kind=rk) :: swimmer_period2 = 0. !< oscillation period of swimmer 2
  real(kind=rk) :: swimmer_period3 = 0. !< oscillation period of swimmer 3
  real(kind=rk) :: swimmer_phase1 = 0. !< phase shift of swimmer 1
  real(kind=rk) :: swimmer_phase2 = 0. !< phase shift of swimmer 2
  real(kind=rk) :: swimmer_phase3 = 0. !< phase shift of swimmer 3
  logical :: swimmer_manualpos = .false. !< position swimmer manually
  real(kind=rk) :: swimmer_pos(3) = (/0., 0., 0./) !< initial swimmer position (central bead)
  real(kind=rk) :: swimmer_axis(3) = (/1., 0., 0./) !< initial swimmer axis
#endif

  ! Standard IBM namelist.
  namelist /IBM_INPUT/ IBM_RANGE, mode_check_rock_pos, mode_polydispersity, mode_manual_pos, mode_grow_particles, &
            & rel_initial_radius, num_steps_growth, growth_bending_modulus, &
            & dump_vtk, time_step_dump_vtk, dump_profiles, time_step_dump_profiles, dump_particles, time_step_dump_particles, &
            & dump_indexfield, time_step_dump_indexfield, &
            & num_meshes, int_strength_part_part, int_strength_part_wall, int_strength_gradient, friction, IBM_massdensity, &
            & viscosity_contrast, time_step_rebuild_index, selective_wall_force_magnitude, selective_wall_force_range, &
            & ibm_colour, g_ri, g_bi, noslipcoeff, surften, contactangle, width, penalty_red_blue

  ! Swimmer namelist.
#ifdef IBM_SWIMMER
  namelist /IBM_SWIMMERS/ swimmer_d12, swimmer_d23, swimmer_modulus12, swimmer_modulus23, &
            & swimmer_kb, swimmer_ka, swimmer_anchor, swimmer_force_mag1, swimmer_force_mag2, swimmer_force_mag3, &
            & swimmer_period1, swimmer_period2, swimmer_period3, swimmer_phase1, swimmer_phase2, swimmer_phase3, &
            & swimmer_manualpos, swimmer_pos, swimmer_axis
#endif
            
  ! Variables used during initialization and the time loop
  ! Some of these variables have to be saved in a checkpoint.
  integer :: IBM_HALO !< halo width for IBM communication
  integer :: num_particles_gl !< total number of particles in simulation
  integer :: num_particles_loc !< number of particles in subdomain (CHECKPOINT!)
  integer :: num_nodes_loc_vtk !< total number of nodes of all particles in subdomain (used for VTK output)
  integer :: num_faces_loc_vtk !< total number of faces of all particles in subdomain (used for VTK output)
  integer :: dist_min_gl !< minimum distance between two nodes in global domain
  integer :: dist_min_loc !< minimum distance between two nodes in local subdomain
  integer :: num_nodes_local !< number of nodes in local node list
  type(superstruct_node), dimension(MAX_NUM_NODES) :: nodes_local !< list of locally available nodes
  type(dump_static_particle_info), dimension(NUM_PART_LOC_MAX) :: static_particle_data_loc !< local static particle data to dump
  real(kind=rk), allocatable, dimension(:, :, :, :) :: vel_phys !< physical lattice velocity
  real(kind=rk), allocatable, dimension(:, :, :, :) :: force_IBM !< force due to IBM (used to compute the particle stress)
  integer, allocatable, dimension(:, :, :, :) :: LUT_nodes_in_range !< look-up table for nearby nodes
#ifdef IBM_INDEXFIELD
  real(kind=rk), allocatable, dimension(:, :, :) :: checked_index !< status array for index check
  real(kind=rk), allocatable, dimension(:, :, :) :: normal_distance !< normal distance to membrane
  real(kind=rk), allocatable, dimension(:, :, :) :: interior_index !< index indicating interior/exterior (current value)
  real(kind=rk), allocatable, dimension(:, :, :) :: interior_index_old !< index indicating interior/exterior (previous value)
#endif
#ifdef IBM_BINARYIBM
  real(kind=rk), allocatable, dimension(:, :, :, :) :: vel_red !< red fluid velocity
  real(kind=rk), allocatable, dimension(:, :, :, :) :: vel_blue !< blue fluid velocity
  real(kind=rk), allocatable, dimension(:, :, :, :) :: vel_diff !< difference of red and blue fluid velocities
#endif

  ! MPI derived datatypes
  ! These variables do not have to be saved in a checkpoint.
  integer :: MPI_NODE_DATA !< datatype for distributing nodes to neighboring processes
  integer :: MPI_LAGRANGE_NODE_DATA !< datatype for sending particle node data to a neighboring process
  integer :: MPI_LAGRANGE_FACE_DATA !< datatype for sending particle face data to a neighboring process
  integer :: MPI_LAGRANGE_PARTICLE_DATA !< datatype for sending particle header data to a neighboring process
  integer :: MPI_PARTICLE_INIT !< datatype for broadcasting particle positions during initialization
  integer :: MPI_STATIC_PARTICLE_DATA !< datatype for collecting static particle data
  integer :: datatype_send_vel_xpos, datatype_send_vel_xneg !< datatypes for physical velocity halo exchange
  integer :: datatype_send_vel_ypos, datatype_send_vel_yneg
  integer :: datatype_send_vel_zpos, datatype_send_vel_zneg
  integer :: datatype_recv_vel_xpos, datatype_recv_vel_xneg
  integer :: datatype_recv_vel_ypos, datatype_recv_vel_yneg
  integer :: datatype_recv_vel_zpos, datatype_recv_vel_zneg
  integer :: datatype_send_rock_grad_xpos, datatype_send_rock_grad_xneg !< datatypes for rock gradient halo exchange
  integer :: datatype_send_rock_grad_ypos, datatype_send_rock_grad_yneg
  integer :: datatype_send_rock_grad_zpos, datatype_send_rock_grad_zneg
  integer :: datatype_recv_rock_grad_xpos, datatype_recv_rock_grad_xneg
  integer :: datatype_recv_rock_grad_ypos, datatype_recv_rock_grad_yneg
  integer :: datatype_recv_rock_grad_zpos, datatype_recv_rock_grad_zneg
  integer :: datatype_send_index_xpos, datatype_send_index_xneg !< datatypes for IBM index halo exchange
  integer :: datatype_send_index_ypos, datatype_send_index_yneg
  integer :: datatype_send_index_zpos, datatype_send_index_zneg
  integer :: datatype_recv_index_xpos, datatype_recv_index_xneg
  integer :: datatype_recv_index_ypos, datatype_recv_index_yneg
  integer :: datatype_recv_index_zpos, datatype_recv_index_zneg

  ! Variables related to initialization or MPI communication
  ! These variables do not have to be saved in a checkpoint.
  type(superstruct_node), dimension(MAX_NUM_NODES_SEND) :: & !< buffers for sending nodes to neighbors
    & send_in_neg_x, send_in_pos_x, send_in_neg_y, send_in_pos_y, send_in_neg_z, send_in_pos_z
  type(superstruct_node), dimension(MAX_NUM_NODES_SEND) :: & !< buffers for receiving nodes from neighbors
    & recv_from_neg_x, recv_from_pos_x, recv_from_neg_y, recv_from_pos_y, recv_from_neg_z, recv_from_pos_z
  type(init_particle_info), allocatable, dimension(:) :: particle_data !< array for broadcasting particle data during initialization
  type(init_mesh_info), allocatable, dimension(:) :: mesh_data !< array for broadcasting mesh data during initialization
  integer, allocatable, dimension(:, :, :) :: rock_grad !< distance of rock node to fluid
  integer :: rock_grad_max !< maximum value of rock gradient
  real(kind=rk), allocatable, dimension(:, :, :, :) :: gradient !< gradient direction for force

! endif IBM_PART
#endif

end module lsuperstruct_data_module
