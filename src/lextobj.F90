!> Lagrangian extended object module
!>
!> This module contains all structures, variables, and subroutines responsible
!> for initializing and maintaining the Lagrangian extended objects (particles).

#include "lbe.h"

module lextobj_module

#ifdef IBM_PART

  ! Include external modules.
  use lbe_globals_module, only : pi, rk
  use lbe_parallel_module, only : tnx
  use lbe_helper_module, only : cross_product, norm, check_NaN, n_sites_fluid, diag_symm_matrix, present_and_non_negative
  use lbe_log_module
  use lmesh_module, only : meshes, MAX_NODE_NEIGHBORS, mesh_struct, node_diff, nodes_same

  implicit none
  private
  public :: extobj_struct, particles, part_ind, init_particle, deinit_particle, forces_compute, update_node_positions, &
            & update_face_areas, update_node_areas, update_particle_center, update_momenta, set_radius, &
            & set_center_position, rotate_particle, particle_node, particle_face, forces_friction, compute_angular_velocity, &
            & compute_particle_stresslet, compute_inertia_tensor, compute_inertia_ellipsoid, compute_volume_slice, &
            & compute_area_slice
#ifdef IBM_INDEXFIELD
  public :: update_node_normals
#endif

  !> Definition of the particle node structure
  !>
  !> It contains all properties which are unique for a particle node,
  !> e.g., its position, velocity, or force.

  type particle_node
    real(kind=rk) :: area !< node area (equal to Voronoi area)
    real(kind=rk), dimension(3) :: pos !< current node position
    real(kind=rk), dimension(3) :: pos_old !< former node position
    real(kind=rk), dimension(3) :: vel !< node velocity
    real(kind=rk), dimension(3) :: force_tot !< total force
    real(kind=rk), dimension(3) :: force_s !< strain force
    real(kind=rk), dimension(3) :: force_b !< bending force
    real(kind=rk), dimension(3) :: force_v !< volume force
    real(kind=rk), dimension(3) :: force_at !< total surface force
    real(kind=rk), dimension(3) :: force_int !< interaction force
    real(kind=rk), dimension(3) :: normal !< node normal vector
    real(kind=rk) :: curv_radius !< node normal vector
#ifdef IBM_DRAG
    real(kind=rk), dimension(3) :: force_c !< constant node forcing
#endif
#ifdef IBM_FIXED
    real(kind=rk), dimension(3) :: pos_anchor !< anchor node position
    real(kind=rk), dimension(3) :: force_anchor !< spring force anchoring particle
#endif
  end type particle_node

  !> Definition of the particle face structure
  !>
  !> It contains all properties which are unique for a particle face,
  !> e.g., its area and normal vector.

  type particle_face
    real(kind=rk) :: area !< face area
    real(kind=rk), dimension(3) :: normal !< face unit normal vector
  end type particle_face

  !> Definition of the extended object structure
  !>
  !> It contains all properties which are unique for a particle,
  !> e.g., its position or deformation state.

  type extobj_struct
    ! Variables with known size.
    integer :: num_nodes !< number of nodes
    integer :: num_faces !< number of faces
    integer :: mesh_type !< mesh belonging to particle
    integer :: particle_index_gl !< global particle index (unique for each particle, also in parallel simulation)
    integer :: currently_used = 0 !< switch indicating the particle usage
    integer, dimension(3) :: num_jumps !< number of jumps through periodic boundaries
    real(kind=rk) :: radius_0 !< reference radius
    real(kind=rk) :: radius !< current radius
    real(kind=rk) :: volume_0 !< reference volume
    real(kind=rk) :: volume !< current volume
    real(kind=rk) :: surface_0 !< reference surface
    real(kind=rk) :: surface !< current surface
    real(kind=rk) :: k_s !< strain modulus
    real(kind=rk) :: k_al !< local area modulus
    real(kind=rk) :: k_b !< bending modulus
    real(kind=rk) :: k_v !< volume modulus
    real(kind=rk) :: k_at !< surface modulus
    real(kind=rk) :: erg_tot !< total energy
    real(kind=rk) :: erg_s !< strain energy
    real(kind=rk) :: erg_b !< bending energy
    real(kind=rk) :: erg_v !< volume energy
    real(kind=rk) :: erg_at !< surface energy
    real(kind=rk) :: erg_int !< interaction energy
    real(kind=rk) :: stresslet_xx !< stresslet (xx-component)
    real(kind=rk) :: stresslet_xy !< stresslet (xy-component)
    real(kind=rk) :: stresslet_xz !< stresslet (xz-component)
    real(kind=rk) :: stresslet_yy !< stresslet (yy-component)
    real(kind=rk) :: stresslet_yz !< stresslet (yz-component)
    real(kind=rk) :: stresslet_zz !< stresslet (zz-component)
    real(kind=rk), dimension(3) :: center !< current particle center position
    real(kind=rk), dimension(3) :: center_old !< previous particle center position
    real(kind=rk), dimension(3) :: force_total !< total force acting on particle
    real(kind=rk), dimension(3) :: torque_total !< total torque acting on particle
    real(kind=rk), dimension(3) :: linear_velocity !< linear particle velocity
    real(kind=rk), dimension(3) :: linear_momentum !< linear particle momentum
    real(kind=rk), dimension(3) :: angular_velocity !< angular particle velocity
    real(kind=rk), dimension(3) :: angular_momentum !< angular particle momentum
    real(kind=rk), dimension(3, 3) :: inertia_tensor !< inertia tensor
    real(kind=rk), dimension(3) :: inertia_ell_axes !< semiaxes of inertia ellipsoid
    real(kind=rk), dimension(3, 3) :: inertia_ell_orient !< orientation vectors of inertia ellipsoid
#ifdef IBM_FIXED
    real(kind=rk) :: k_anchor !< anchor spring constant
		real(kind=rk), dimension(3) :: force_anchor !< total force acting on particle
#endif
#ifdef IBM_DRAG
    real(kind=rk), dimension(3) :: force_const !< constant force enforces onto mesh-nodes
#endif
    


    ! Variables with variable size.
    type(particle_node), allocatable, dimension(:) :: node !< nodes
    type(particle_face), allocatable, dimension(:) :: face !< faces
  end type extobj_struct

  ! Variable declarations
  ! These variables have to be saved in a checkpoint.
  type(extobj_struct), allocatable, dimension(:) :: particles !< array containing all particles required during simulation
  integer, allocatable, dimension(:) :: part_ind !< array for local particle indices (user interface)

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialization of a new particle
  !>
  !> The numbers of nodes and faces are set. The allocatable arrays are resized.
  !> Particle properties (elasticities) are set. The particle is resized and positioned temporarily.

  subroutine init_particle(particle, mesh)
    type(extobj_struct), intent(inout) :: particle !< element of the particle array
    type(mesh_struct), intent(in) :: mesh !< element of the mesh array

    ! Declare variables.
    integer :: n_i
    real(kind=rk), dimension(3) :: origin ! position of the origin

    ! Activate particle.
    particle%currently_used = 1

    ! Set mesh index and number of nodes/faces.
    particle%particle_index_gl = -1
    particle%mesh_type = mesh%mesh_type
    particle%num_nodes = mesh%num_nodes
    particle%num_faces = mesh%num_faces

    ! Resize node and face arrays.
    allocate(particle%node(particle%num_nodes))
    allocate(particle%face(particle%num_faces))

    ! Initialize other values.
    particle%num_jumps(:) = 0
    particle%k_s = 0.d0
    particle%k_al = 0.d0
    particle%k_b = 0.d0
    particle%k_v = 0.d0
    particle%k_at = 0.d0
    particle%erg_tot = 0.d0
    particle%erg_s = 0.d0
    particle%erg_b = 0.d0
    particle%erg_v = 0.d0
    particle%erg_at = 0.d0
    particle%erg_int = 0.d0
    particle%center(:) = 0.d0
    particle%center_old(:) = 0.d0
    particle%force_total(:) = 0.d0
    particle%torque_total(:) = 0.d0
    particle%linear_velocity(:) = 0.d0
    particle%linear_momentum(:) = 0.d0
    particle%angular_velocity(:) = 0.d0
    particle%angular_momentum(:) = 0.d0
    particle%inertia_tensor(:, :) = 0.d0

    ! Temporarily position particle nodes according to mesh nodes with zero velocity.
    do n_i = 1, particle%num_nodes
      particle%node(n_i)%pos(:) = mesh%pos(:, n_i)
      particle%node(n_i)%vel(:) = 0.0d0
    end do

    ! Set particle center to the origin.
    origin(:) = 0.d0
    call set_center_position(particle, origin)

    ! Temporarily set particle radius to unity.
    particle%radius = 1.d0
    particle%volume = 0.d0
    particle%surface = 0.d0
    call set_radius(particle, 1.d0, .true.)
  end subroutine init_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Deinitialization of a new particle
  !>
  !> The particle is marked inactive.
  !> Memory is deallocated.

  subroutine deinit_particle(particle)
    type(extobj_struct), intent(inout) :: particle !< element of the particle array

    ! Deactivate particle.
    particle%currently_used = 0
    particle%particle_index_gl = -1 ! use invalid particle index
    particle%mesh_type = -1 ! use invalid mesh index
    particle%num_nodes = 0
    particle%num_faces = 0

    ! Deallocate node and face arrays.
    deallocate(particle%node)
    deallocate(particle%face)
  end subroutine deinit_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update of particle node positions
  !>
  !> The node positions are updated according to Forward Euler.
  !> In principle, other time integration schemes may be used.
  !> The old positions are stored for data analysis (not required for the simulation itself).
  !> NOTE: The periodicity of the simulation is not taken into account here.
  !> It is checked whether the node velocity is reasonable.
  !> If the velocity is too large, first a warning, then an error is reported.

  subroutine update_node_positions(particle, IBM_massdensity, friction)
    type(extobj_struct), intent(inout) :: particle !< element of the particle array
    real(kind=rk), intent(inout) :: IBM_massdensity !< node mass density
    real(kind=rk), intent(inout) :: friction !< friction

    ! Declare variables.
    integer :: n_i ! node index
    real(kind=rk) :: mass ! mass of node
    real(kind=rk) :: C1, C2 ! integrators
    real(kind=rk), dimension(3) :: pos_temp ! half-time-step position
    character(len=200) :: message

    ! Update node positions and check velocities.
    do n_i = 1, particle%num_nodes
      ! If IBM is switched off, node velocities first have to be obtained from forces.
      ! The easiest way is forward Euler.
#ifdef NOFLUIDDYNAMICS
!       mass = IBM_massdensity
      mass = IBM_massdensity * particle%node(n_i)%area
!       C1 = exp(-friction / mass)
!       C2 = mass / friction * (1.0 - exp(-friction / mass))
!       particle%node(n_i)%pos_old = particle%node(n_i)%pos
!       pos_temp = particle%node(n_i)%pos_old + 0.5 * particle%node(n_i)%vel
!       particle%node(n_i)%vel = C1 * particle%node(n_i)%vel + C2 * particle%node(n_i)%force_tot / mass
!       particle%node(n_i)%pos = pos_temp + 0.5 * particle%node(n_i)%vel
      particle%node(n_i)%pos_old = particle%node(n_i)%pos
      pos_temp = particle%node(n_i)%pos_old + 0.5 * particle%node(n_i)%vel
      particle%node(n_i)%vel = (1.0 - friction) * particle%node(n_i)%vel + particle%node(n_i)%force_tot / mass
      particle%node(n_i)%pos = pos_temp + 0.5 * particle%node(n_i)%vel
#else
      ! Store old node position.
      particle%node(n_i)%pos_old = particle%node(n_i)%pos
      
      ! Update node position.
      particle%node(n_i)%pos = particle%node(n_i)%pos + particle%node(n_i)%vel
#endif

      ! Check node velocity.
      ! If the velocity is large (>0.5), a warning is printed along with the force information for that node.
      ! If the velocity is even larger (>1.0), the simulation is terminated due to instability problems.
      if(norm(particle%node(n_i)%vel) .gt. 0.5) then
        write(message, "('velocity of node ', i0, ' in particle ', i0, ' of mesh type ', i0, ' exceeds 0.5')") n_i, &
            & particle%particle_index_gl, particle%mesh_type
        call log_msg(trim(message), .true.)
        write(message, "('pos: (', f0.3, ', ', f0.3, ', ', f0.3, ')')") particle%node(n_i)%pos
        call log_msg(trim(message), .true.)
        write(message, "('vel: (', f0.3, ', ', f0.3, ', ', f0.3, ')')") particle%node(n_i)%vel
        call log_msg(trim(message), .true.)
        write(message, "('force_s: (', f0.3, ', ', f0.3, ', ', f0.3, ')')") particle%node(n_i)%force_s
        call log_msg(trim(message), .true.)
        write(message, "('force_b: (', f0.3, ', ', f0.3, ', ', f0.3, ')')") particle%node(n_i)%force_b
        call log_msg(trim(message), .true.)
        write(message, "('force_v: (', f0.3, ', ', f0.3, ', ', f0.3, ')')") particle%node(n_i)%force_v
        call log_msg(trim(message), .true.)
        write(message, "('force_at: (', f0.3, ', ', f0.3, ', ', f0.3, ')')") particle%node(n_i)%force_at
        call log_msg(trim(message), .true.)
        write(message, "('force_int: (', f0.3, ', ', f0.3, ', ', f0.3, ')')") particle%node(n_i)%force_int
        call log_msg(trim(message), .true.)
#ifdef IBM_DRAG
        write(message, "('force_c: (', f0.3, ', ', f0.3, ', ', f0.3, ')')") particle%node(n_i)%force_c
        call log_msg(trim(message), .true.)
#endif
#ifdef IBM_FIXED
        write(message, "('force_anchor: (', f0.3, ', ', f0.3, ', ', f0.3, ')')") particle%node(n_i)%force_anchor
        call log_msg(trim(message), .true.)
#endif
        if(norm(particle%node(n_i)%vel) .gt. 1.0) then
          call error("node velocity too large, simulation unstable")
        end if
      end if

      if(check_NaN(particle%node(n_i)%pos(1)) .or. check_NaN(particle%node(n_i)%pos(2)) &
        & .or. check_NaN(particle%node(n_i)%pos(3))) then
        write(message, "('Nan detected for node ', i0, ' in particle ', i0, ' of mesh type ', i0)") n_i, &
            & particle%particle_index_gl, particle%mesh_type
        call log_msg(trim(message), .true.)
        write(message, "('pos: (', f0.3, ', ', f0.3, ', ', f0.3, ')')") particle%node(n_i)%pos
        call log_msg(trim(message), .true.)
        write(message, "('vel: (', f0.3, ', ', f0.3, ', ', f0.3, ')')") particle%node(n_i)%vel
        call log_msg(trim(message), .true.)
        write(message, "('force_s: (', f0.3, ', ', f0.3, ', ', f0.3, ')')") particle%node(n_i)%force_s
        call log_msg(trim(message), .true.)
        write(message, "('force_b: (', f0.3, ', ', f0.3, ', ', f0.3, ')')") particle%node(n_i)%force_b
        call log_msg(trim(message), .true.)
        write(message, "('force_v: (', f0.3, ', ', f0.3, ', ', f0.3, ')')") particle%node(n_i)%force_v
        call log_msg(trim(message), .true.)
        write(message, "('force_at: (', f0.3, ', ', f0.3, ', ', f0.3, ')')") particle%node(n_i)%force_at
        call log_msg(trim(message), .true.)
        write(message, "('force_int: (', f0.3, ', ', f0.3, ', ', f0.3, ')')") particle%node(n_i)%force_int
        call log_msg(trim(message), .true.)
#ifdef IBM_DRAG
        write(message, "('force_c: (', f0.3, ', ', f0.3, ', ', f0.3, ')')") particle%node(n_i)%force_c
        call log_msg(trim(message), .true.)
#endif
#ifdef IBM_FIXED
        write(message, "('force_anchor: (', f0.3, ', ', f0.3, ', ', f0.3, ')')") particle%node(n_i)%force_anchor
        call log_msg(trim(message), .true.)
#endif
        call error("NaN detected, simulation unstable")
      end if
    end do
  end subroutine update_node_positions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update of particle face areas, normals, surface, and volume
  !>
  !> Based on the new node positions, the face normals and their areas are updated.
  !> The total particle surface and volume are updated as well.
  !> NOTE: The normal vectors point in outward direction. This is assured by the ordering of the nodes.

  subroutine update_face_areas(particle)
    type(extobj_struct), intent(inout) :: particle !< element of the particle array

    ! Declare variables.
    integer :: f_i ! face index
    integer :: n_1, n_2, n_3 ! node indices
    real(kind=rk) :: normal_length ! length of normal vector
    real(kind=rk) :: area ! face area
    real(kind=rk), dimension(3) :: normal ! face normal
    real(kind=rk), dimension(3) :: pos_1, pos_2, pos_3 ! node positions
    real(kind=rk), dimension(3) :: face_centroid ! centroid of face

    ! Reset particle surface and volume.
    particle%surface = 0.d0
    particle%volume = 0.d0

    ! Update particle face areas, normals, surface, and volume.
    do f_i = 1, particle%num_faces
      ! Get node indices and positions.
      n_1 = meshes(particle%mesh_type)%neighbor_face_node(f_i, 1)
      n_2 = meshes(particle%mesh_type)%neighbor_face_node(f_i, 2)
      n_3 = meshes(particle%mesh_type)%neighbor_face_node(f_i, 3)
      pos_1(:) = particle%node(n_1)%pos(:)
      pos_2(:) = particle%node(n_2)%pos(:)
      pos_3(:) = particle%node(n_3)%pos(:)

      ! Compute face normal and normalize to unit length.
      normal = cross_product((pos_1 - pos_2), (pos_3 - pos_2))
      normal_length = norm(normal)
      normal(:) = normal(:) / normal_length

      ! Compute face centroid and area.
      face_centroid = (pos_1 + pos_2 + pos_3) / 3.d0 - particle%center
      area = 0.5d0 * normal_length

      ! Update particle properties.
      particle%face(f_i)%normal(:) = normal(:)
      particle%face(f_i)%area = area
      particle%surface = particle%surface + area
      particle%volume = particle%volume + (area * dot_product(normal, face_centroid)) / 3.d0
    end do
  end subroutine update_face_areas

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update of particle node areas
  !>
  !> The node area is proportional to the Voronoi area of the node.
  !> In order to find the area of a node, all neighboring face areas have to be considered.

  subroutine update_node_areas(particle)
    type(extobj_struct), intent(inout) :: particle !< element of the particle array

    ! Declare variables.
    integer :: f_i ! face index
    integer :: n_i, n_1, n_2, n_3 ! node indices

    do n_i = 1, particle%num_nodes
      particle%node(n_i)%area = 0.0d0
    end do

    ! Run over faces and distribute area to nodes.
    do f_i = 1, particle%num_faces
      ! Identify nodes in face.
      n_1 = meshes(particle%mesh_type)%neighbor_face_node(f_i, 1)
      n_2 = meshes(particle%mesh_type)%neighbor_face_node(f_i, 2)
      n_3 = meshes(particle%mesh_type)%neighbor_face_node(f_i, 3)

      ! Update node areas.
      particle%node(n_1)%area = particle%node(n_1)%area + particle%face(f_i)%area / 3.d0
      particle%node(n_2)%area = particle%node(n_2)%area + particle%face(f_i)%area / 3.d0
      particle%node(n_3)%area = particle%node(n_3)%area + particle%face(f_i)%area / 3.d0
    end do
  end subroutine update_node_areas

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update of particle node normal vector and curvature radius
  !>
  !> The node normal and curvature radius are used for the identification of interior/exterior fluid nodes.
  !> This functionality is only required in the IBM_INDEXFIELD mode.
  !> The node normal is assumed to be the weighted average of the normal vector of all neighboring faces.
  !> The curvature is obtained by finding a spherical cap approximating the surface locally.

#ifdef IBM_INDEXFIELD
  subroutine update_node_normals(particle)
    type(extobj_struct), intent(inout) :: particle !< element of the particle array

    ! Declare variables.
    integer :: n_i, n_j ! node indices
    integer :: f_j ! face index
    integer :: i ! counter
    real(kind=rk), dimension(3) :: center ! center position of neighbors
    real(kind=rk) :: base_radius
    real(kind=rk), dimension(3) :: height

    do n_i = 1, particle%num_nodes
      ! Reset center position.
      center(:) = 0.0d0
      particle%node(n_i)%normal(:) = 0.0d0
      particle%node(n_i)%curv_radius = 0.0d0

      do i = 1, MAX_NODE_NEIGHBORS
        ! Identify node neighbor.
        n_j = meshes(particle%mesh_type)%neighbor_node_node(n_i, i)
        f_j = meshes(particle%mesh_type)%neighbor_node_face(n_i, i)

        ! Update center and normal.
        ! The normal is weighted by the area of the corresponding face.
        if(n_j > 0) then
          center(:) = center(:) + particle%node(n_j)%pos(:)
          particle%node(n_i)%normal(:) = particle%node(n_i)%normal(:) + particle%face(f_j)%normal(:) * particle%face(f_j)%area
        else
          exit
        end if
      end do

      ! Normalize center and normal.
      center(:) = center(:) / (i - 1)
      particle%node(n_i)%normal(:) = particle%node(n_i)%normal(:) / norm(particle%node(n_i)%normal(:))

      ! Compute base radius.
      base_radius = 0.0d0

      do i = 1, MAX_NODE_NEIGHBORS
        ! Identify node neighbor.
        n_j = meshes(particle%mesh_type)%neighbor_node_node(n_i, i)

        ! Update center and normal.
        if(n_j > 0) then
          base_radius = base_radius + norm(particle%node(n_j)%pos(:) - center(:))
        else
          exit
        end if
      end do

      base_radius = base_radius / (i - 1)

      ! Compute curvature radius.
      height(:) = particle%node(n_i)%pos(:) - center(:)
      particle%node(n_i)%curv_radius = (norm(height)**2 + base_radius**2) / (2.0d0 * norm(height))

      if(dot_product(height, particle%node(n_i)%normal) < 0.0d0) then
        particle%node(n_i)%curv_radius = -particle%node(n_i)%curv_radius
      end if
    end do
  end subroutine update_node_normals
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update of particle center
  !>
  !> The location of the particle center is updated.
  !> Periodicity is not taken into account here.
  !> The center location is computed as area-weighted average of all node locations.
  !> The old center position is stored for data analysis (not required for the simulation itself).

  subroutine update_particle_center(particle)
    type(extobj_struct), intent(inout) :: particle !< element of the particle array

    ! Declare variables.
    integer :: n_i ! node index

    ! Update particle center.
    particle%center_old(:) = particle%center(:)
    particle%center(:) = 0.d0

    do n_i = 1, particle%num_nodes
      particle%center = particle%center + (particle%node(n_i)%pos * particle%node(n_i)%area)
    end do

    particle%center = particle%center / particle%surface
  end subroutine update_particle_center

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set particle center to desired position
  !>
  !> The node positions are updated accordingly.
  !> Periodicity is not taken into account here: The user has to make sure that new center position is valid.

  subroutine set_center_position(particle, center_new)
    type(extobj_struct), intent(inout) :: particle !< element of the particle array
    real(kind=rk), dimension(3), intent(in) :: center_new

    ! Declare variables.
    integer :: n_i ! node index

    ! Shift each node according to the new center.
    do n_i = 1, particle%num_nodes
      particle%node(n_i)%pos = particle%node(n_i)%pos - particle%center + center_new
    end do

    ! Set new particle center.
    particle%center = center_new
  end subroutine set_center_position

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Resize particle shape
  !>
  !> The particle is resized radially by a factor (unity corresponds to no change).
  !> The center and the equilibrium properties of the membrane are not modified.

  subroutine resize_particle(particle, factor)
    type(extobj_struct), intent(inout) :: particle !< element of the particle array
    real(kind=rk), intent(in) :: factor !< scaling factor for the radius

    ! Declare variables.
    integer :: n_i ! node index

    ! Change node positions, but keep the center position invariant.
    do n_i = 1, particle%num_nodes
      particle%node(n_i)%pos = (particle%node(n_i)%pos - particle%center) * factor + particle%center
    end do
  end subroutine resize_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set particle radius
  !>
  !> The radius of the particle is set to a new value.
  !> Consequently, the equilibrium surface and volume of the membrane are updated.
  !> The user has to decide whether the reference radius 'radius_0' has to be changed as well or not.
  !> This is particularly important for the growth of the particles at the beginning of a simulation.
  !> TODO: I am not sure that the modifications below are correct. Has to be checked in a benchmark simulation with particle growth.

  subroutine set_radius(particle, radius_new, switch_radius_0)
    type(extobj_struct), intent(inout) :: particle !< element of the particle array
    real(kind=rk), intent(in) :: radius_new !< new particle radius
    logical, intent(in) :: switch_radius_0 !< whether to change reference radius or not

    ! Declare variables.
    real(kind=rk) :: radius_temp ! temporary radius

    ! Change radius
    radius_temp = particle%radius
    particle%radius = radius_new
    call resize_particle(particle, radius_new / radius_temp)

    if(switch_radius_0) then
      particle%radius_0 = radius_new
      particle%surface_0 = meshes(particle%mesh_type)%surface * (radius_new / radius_temp)**2
      particle%volume_0 = meshes(particle%mesh_type)%volume * (radius_new / radius_temp)**3
    end if
      
    call update_face_areas(particle)
    call update_node_areas(particle)
  end subroutine set_radius

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update particle momenta
  !>
  !> The momenta of the membrane are updated:
  !> - linear and angular momentum (computed from new positions and velocities)
  !> - linear and angular velocity (computed from new positions and velocities)
  !> - total membrane force (computed from old forces)
  !> - total membrane torque (computed from old positions and forces)
  !> These quantities are not required for the computations,
  !> but they are used for monitoring the simulation and data analysis.
  !> The surface takes the role of the mass.
  !> There is no functionality to compute the angular velocity of a deformable particle.
  !> This is rather a conceptual problem since a strict mathematical definition seems to be lacking.

  subroutine update_momenta(particle)
    type(extobj_struct), intent(inout) :: particle !< element of the particle array

    ! Declare variables.
    integer :: n_i ! node index

    ! Reset momenta.
    particle%linear_momentum(:) = 0.d0
    particle%angular_momentum(:) = 0.d0
    particle%force_total(:) = 0.d0
    particle%torque_total(:) = 0.d0

    ! Update momenta.
    do n_i = 1, particle%num_nodes
      particle%linear_momentum = particle%linear_momentum + (particle%node(n_i)%vel * particle%node(n_i)%area)
      particle%angular_momentum = particle%angular_momentum &
        + cross_product((particle%node(n_i)%pos - particle%center), particle%node(n_i)%vel * particle%node(n_i)%area)
      particle%force_total = particle%force_total + particle%node(n_i)%force_tot
      particle%torque_total = particle%torque_total &
        + cross_product((particle%node(n_i)%pos_old - particle%center_old), particle%node(n_i)%force_tot)
    end do

    particle%linear_velocity = particle%linear_momentum / particle%surface
  end subroutine update_momenta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Rotate particle
  !>
  !> The orientation of the particle is rotated about an axis by an angle.
  !> The angle is given in degrees and has to be converted to radiants.
  !> A rotation matrix is computed from the rotation axis and the angle.
  !> The particle center and other properties are not changed.

  subroutine rotate_particle(particle, angle, axis)
    type(extobj_struct), intent(inout) :: particle !< element of the particle array
    real(kind=rk), intent(in) :: angle !< rotation angle in degrees
    real(kind=rk), dimension(3), intent(in) :: axis !< rotation axis

    ! Declare variables.
    integer :: n_i ! node index
    real(kind=rk), dimension(3, 3) :: rot_matrix ! rotation matrix
    real(kind=rk), dimension(3) :: axis_norm ! normalized axis
    real(kind=rk) :: angle_rad ! angle in radiants

    ! Normalize rotation axis to unit length.
    axis_norm = axis / norm(axis)

    ! Convert angle to radiants
    angle_rad = angle * pi / 180.d0

    ! Create rotation matrix
    rot_matrix = reshape( &
      & (/cos(angle_rad) + axis_norm(1) * axis_norm(1) * (1.d0 - cos(angle_rad)), &
      & axis_norm(1) * axis_norm(2) * (1.d0 - cos(angle_rad)) - axis_norm(3) * sin(angle_rad), &
      & axis_norm(1) * axis_norm(3) * (1.d0 - cos(angle_rad)) + axis_norm(2) * sin(angle_rad), &
      & axis_norm(2) * axis_norm(1) * (1.d0 - cos(angle_rad)) + axis_norm(3) * sin(angle_rad), &
      & cos(angle_rad) + axis_norm(2) * axis_norm(2) * (1.d0 - cos(angle_rad)), &
      & axis_norm(2) * axis_norm(3) * (1.d0 - cos(angle_rad)) - axis_norm(1) * sin(angle_rad), &
      & axis_norm(3) * axis_norm(1) * (1.d0 - cos(angle_rad)) - axis_norm(2) * sin(angle_rad), &
      & axis_norm(3) * axis_norm(2) * (1.d0 - cos(angle_rad)) + axis_norm(1) * sin(angle_rad), &
      & cos(angle_rad) + axis_norm(3) * axis_norm(3) * (1.d0 - cos(angle_rad))/) &
      &, (/3,3/))

    ! Perform matrix multiplication.
    do n_i = 1, particle%num_nodes
      particle%node(n_i)%pos = matmul(rot_matrix, (particle%node(n_i)%pos - particle%center)) + particle%center
    end do
  end subroutine rotate_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute angle between two normal vectors
  !>
  !> The angle between two face unit normal vectors is computed.
  !> To increase numerical stability, extremely small angles are neglected:
  !> - If the cosine is near +1, the angle is set to 0.
  !> - If the cosine is near -1, the angle is set to pi.

  function angle_normals(particle, f_1, f_2)
    real(kind=rk) :: angle_normals !< return value
    type(extobj_struct), intent(in) :: particle !< element of the particle array
    integer :: f_1, f_2 !< face indices

    ! Declare variables.
    real(kind=rk) :: cosine ! cosine of angle

    ! Compute cosine of angle between face normals.
    cosine = dot_product(particle%face(f_1)%normal, particle%face(f_2)%normal)

    ! Compute angle from cosine.
    if(cosine .ge. 0.999999d0) then
      angle_normals = 0
    else if(cosine .le. -0.999999d0) then
      angle_normals = pi
    else
      angle_normals = acos(cosine)
    end if
  end function angle_normals

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute angle between two edges
  !>
  !> The angle between two edges is computed.
  !> The edges are specified by three nodes where the angle is taken at the second node.
  !> To increase numerical stability, extremely small angles are neglected.
  !> - If the cosine is near +1, the angle is set to 0.
  !> - If the cosine is near -1, the angle is set to pi.

  function angle_edges(particle, n_1, n_2, n_3)
    real(kind=rk) :: angle_edges !< return value
    type(extobj_struct), intent(in) :: particle !< element of the particle array
    integer :: n_1, n_2, n_3 !< node indices

    ! Declare variables.
    real(kind=rk), dimension(3) :: edge_1 ! first edge
    real(kind=rk), dimension(3) :: edge_2 ! second edge
    real(kind=rk) :: cosine ! cosine of angle between edges

    ! Compute edges and cosine of angle between them.
    edge_1 = particle%node(n_3)%pos - particle%node(n_2)%pos
    edge_2 = particle%node(n_1)%pos - particle%node(n_2)%pos
    cosine = dot_product(edge_1, edge_2) / (norm(edge_1) * norm(edge_2))

    ! Compute angle from cosine.
    if(cosine .ge. 0.999999d0) then
      angle_edges = 0
    else if(cosine .le. -0.999999d0) then
      angle_edges = pi
    else
      angle_edges = acos(cosine)
    end if
  end function angle_edges

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute internal particle forces
  !>
  !> The force computation subroutines are called.
  !> NOTE: \c forces_purge must always be called at the beginning in order to reset the membrane forces.

  subroutine forces_compute(particle, surften, contactangle, width, friction, time, k_b)
    type(extobj_struct), intent(inout) :: particle !< element of the particle array
    real(kind=rk), intent(in) :: surften !< emulated surface tension
    real(kind=rk), intent(in) :: contactangle !< contact angle
    real(kind=rk), intent(in) :: width !< width of emulated interface
    real(kind=rk), intent(in) :: friction !< friction coefficient
    integer, intent(in) :: time !< current time step
    real(kind=rk), intent(in), optional :: k_b !< overrides \c particle%k_b if specified and non-negative

    ! Call subroutines for force component computations.
    call forces_purge(particle) ! reset force arrays
    call forces_volume(particle) ! volume force
    call forces_surface(particle) ! surface forces
    call forces_strain(particle) ! strain forces
!     call forces_area(particle) ! area forces
    call forces_bending(particle, k_b) ! bending forces
#ifdef IBM_DRAG
		call forces_const_node(particle) ! constant node forcing
#endif
#ifdef IBM_FIXED	
		call forces_anchor(particle)
#endif
#ifdef IBM_EMULATETENSION
    call forces_interface_emulation(particle, surften, contactangle, width)
!     call forces_friction(particle, friction)
#endif

  end subroutine forces_compute

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Reset force arrays
  !>
  !> All forces and energies are set to zero.
  !> NOTE: This subroutine must be called at each time step before new forces and energies are computed.

  subroutine forces_purge(particle)
    type(extobj_struct), intent(inout) :: particle !< element of the particle array

    ! Declare variables.
    integer :: n_i ! node index

    ! Reset force arrays.
    do n_i = 1, particle%num_nodes
      particle%node(n_i)%force_tot(:) = 0.d0
      particle%node(n_i)%force_s(:) = 0.d0
      particle%node(n_i)%force_b(:) = 0.d0
      particle%node(n_i)%force_v(:) = 0.d0
      particle%node(n_i)%force_at(:) = 0.d0
      particle%node(n_i)%force_int(:) = 0.d0
#ifdef IBM_DRAG
			particle%node(n_i)%force_c(:) = 0.d0
#endif
#ifdef IBM_FIXED
      particle%node(n_i)%force_anchor(:) = 0.d0
#endif
    end do

    ! Reset energies.
    particle%erg_tot = 0.d0
    particle%erg_s = 0.d0
    particle%erg_b = 0.d0
    particle%erg_v = 0.d0
    particle%erg_at = 0.d0
    particle%erg_int = 0.d0
  end subroutine forces_purge

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute volume forces
  !>
  !> The volume forces are computed. They are proportional to the relative volume deviation.
  !> The force computation is based on the volume energy (principle of virtual work).
  !> If the volume modulus is quasi zero, the computation is skipped.
  !> The volume force is the sum of all face contributions (details in my dissertation, p. 48 and p. 151).

  subroutine forces_volume(particle)
    type(extobj_struct), intent(inout) :: particle !< element of the particle array

    ! Declare variables.
    real(kind=rk) :: strength ! force strength
    integer :: f_i ! face index
    integer :: n_1, n_2, n_3 ! node indices

    ! Skip computation if the volume modulus is small.
    if(particle%k_v .lt. 1.d-15) return

    ! Compute force strength.
    strength = particle%k_v * (particle%volume / particle%volume_0 - 1.d0) / 6.d0

    ! Compute volume forces.
    do f_i = 1, particle%num_faces
      n_1 = meshes(particle%mesh_type)%neighbor_face_node(f_i, 1)
      n_2 = meshes(particle%mesh_type)%neighbor_face_node(f_i, 2)
      n_3 = meshes(particle%mesh_type)%neighbor_face_node(f_i, 3)
      particle%node(n_1)%force_v = particle%node(n_1)%force_v + &
        & cross_product(particle%node(n_2)%pos - particle%center, particle%node(n_3)%pos - particle%center) * strength
      particle%node(n_2)%force_v = particle%node(n_2)%force_v + &
        & cross_product(particle%node(n_3)%pos - particle%center, particle%node(n_1)%pos - particle%center) * strength
      particle%node(n_3)%force_v = particle%node(n_3)%force_v + &
        & cross_product(particle%node(n_1)%pos - particle%center, particle%node(n_2)%pos - particle%center) * strength
    end do

    ! Compute volume energy and update total particle energy.
    particle%erg_v = 0.5_rk * particle%k_v * (particle%volume - particle%volume_0)**2 / particle%volume_0
    particle%erg_tot = particle%erg_tot + particle%erg_v
  end subroutine forces_volume

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute surface forces
  !>
  !> The surface forces are computed. They are proportional to the relative surface deviation.
  !> The force computation is based on the surface energy (principle of virtual work).
  !> If the surface modulus is quasi zero, the computation is skipped.
  !> The surface force is the sum of all face contributions (details in my dissertation, p. 47 and p. 150).

  subroutine forces_surface(particle)
    type(extobj_struct), intent(inout) :: particle !< element of the particle array

    ! Declare variables.
    real(kind=rk) :: strength ! force strength
    integer :: f_i ! face index
    integer :: n_1, n_2, n_3 ! node indices

    ! Skip computation if the surface modulus is small.
    if(particle%k_at .lt. 1d-15) return

    ! Compute force strength.
    strength = particle%k_at * (particle%surface / particle%surface_0 - 1.d0) / 2.d0

    ! Compute surface forces.
    do f_i = 1, particle%num_faces
      n_1 = meshes(particle%mesh_type)%neighbor_face_node(f_i, 1)
      n_2 = meshes(particle%mesh_type)%neighbor_face_node(f_i, 2)
      n_3 = meshes(particle%mesh_type)%neighbor_face_node(f_i, 3)
      particle%node(n_1)%force_at = particle%node(n_1)%force_at - &
        & cross_product((particle%node(n_3)%pos - particle%node(n_2)%pos), particle%face(f_i)%normal) * strength
      particle%node(n_2)%force_at = particle%node(n_2)%force_at - &
        & cross_product((particle%node(n_1)%pos - particle%node(n_3)%pos), particle%face(f_i)%normal) * strength
      particle%node(n_3)%force_at = particle%node(n_3)%force_at - &
        & cross_product((particle%node(n_2)%pos - particle%node(n_1)%pos), particle%face(f_i)%normal) * strength
    end do

    ! Compute surface energy and update total particle energy.
    particle%erg_at = 0.5d0 * particle%k_at * (particle%surface - particle%surface_0)**2 / particle%surface_0
    particle%erg_tot = particle%erg_tot + particle%erg_at
  end subroutine forces_surface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute area forces
  !>
  !> The area forces are computed.

  subroutine forces_area(particle)
    type(extobj_struct), intent(inout) :: particle !< element of the particle array

    ! Declare variables.
    real(kind=rk) :: strength ! force strength
    real(kind=rk) :: area_0 ! reference area
    integer :: f_i ! face index
    integer :: n_1, n_2, n_3 ! node indices

    ! Skip computation if the surface modulus is small.
    if(particle%k_al .lt. 1d-15) return

    ! Compute force strength.

    ! Compute surface forces.
    do f_i = 1, particle%num_faces
      area_0 = meshes(particle%mesh_type)%area(f_i) * (particle%radius)**2
      strength = particle%k_al * (particle%face(f_i)%area / area_0 * particle%surface_0 / particle%surface - 1.d0) / 2.d0
      n_1 = meshes(particle%mesh_type)%neighbor_face_node(f_i, 1)
      n_2 = meshes(particle%mesh_type)%neighbor_face_node(f_i, 2)
      n_3 = meshes(particle%mesh_type)%neighbor_face_node(f_i, 3)
      particle%node(n_1)%force_s = particle%node(n_1)%force_s - &
        & cross_product((particle%node(n_3)%pos - particle%node(n_2)%pos), particle%face(f_i)%normal) * strength
      particle%node(n_2)%force_s = particle%node(n_2)%force_s - &
        & cross_product((particle%node(n_1)%pos - particle%node(n_3)%pos), particle%face(f_i)%normal) * strength
      particle%node(n_3)%force_s = particle%node(n_3)%force_s - &
        & cross_product((particle%node(n_2)%pos - particle%node(n_1)%pos), particle%face(f_i)%normal) * strength
    end do

    ! Compute surface energy and update total particle energy.
!     particle%erg_at = 0.5d0 * particle%k_at * (particle%surface - particle%surface_0)**2 / particle%surface_0
!     particle%erg_tot = particle%erg_tot + particle%erg_at
  end subroutine forces_area

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute bending forces
  !>
  !> The bending forces are computed. The current angle between neighboring faces is compared to its equilibrium state.
  !> The restoring force is proportional to this deviation.
  !> If the bending modulus is quasi zero, the computation is skipped.
  !> The bending force is the sum of all face pair contributions (details in my dissertation, p. 46 and p. 149).
  !> TODO: It would be a good idea to add functionality to switch between bending models,
  !> e.g. zero or initial reference angles.

  subroutine forces_bending(particle, k_b)
    type(extobj_struct), intent(inout) :: particle !< element of the particle array
    real(kind=rk),intent(in),optional :: k_b !< overrides \c particle%k_b if specified and non-negative

    ! Declare variables.
    real(kind=rk) :: lambda ! effective bending rigidity
    integer :: f_i, f_j ! face indices
    integer :: n ! counter
    real(kind=rk) :: scalar_product ! scalar product of normal vectors
    real(kind=rk), dimension(3) :: normal_ij, normal_ji ! new normal vectors
    real(kind=rk) :: inv_area_i, inv_area_j ! inverse areas (for faster computations)
    integer :: n_1, n_2, n_3, n_4 ! node indices
    real(kind=rk) :: angle_current ! current angle between faces
    integer :: convex ! indicator for convex/concave surface patches
    real(kind=rk) :: strength ! force strength
    real(kind=rk) :: kk_b ! bending force actually used (either original or overwritten value)

    ! k_b overrides the particle's own k_b if it is specified and zero or larger
    if(present_and_non_negative(k_b)) then
       kk_b = k_b
    else
       kk_b = particle%k_b
    end if

    ! Skip computation if the bending modulus is small.
    if(kk_b .lt. 1.d-15) return

    ! Compute effective bending rigidity.
    lambda = kk_b * sqrt(3.d0)

    ! Compute bending forces
    do f_i = 1, particle%num_faces - 1
      do n = 1, 3
        f_j = meshes(particle%mesh_type)%neighbor_face_face(f_i, n)

        ! Each pair of faces is considered only once, hence the force is twice as large.
        ! This additional factor 2 is cancalled with the factor 2 in the denominator of the derivative.
        if(f_i .lt. f_j) then
          ! Compute scalar product of normals and new normal vectors.
          scalar_product = dot_product(particle%face(f_i)%normal, particle%face(f_j)%normal)
          normal_ij = particle%face(f_i)%normal - particle%face(f_j)%normal * scalar_product
          normal_ji = particle%face(f_j)%normal - particle%face(f_i)%normal * scalar_product
          normal_ij = normal_ij / norm(normal_ij)
          normal_ji = normal_ji / norm(normal_ji)

          ! Compute inverse face areas.
          inv_area_i = 1.d0 / particle%face(f_i)%area
          inv_area_j = 1.d0 / particle%face(f_j)%area

          ! Identify common and single nodes: Common nodes are located in both faces, single nodes in only one of both.
          ! n_2 is node in face f_i, but not in f_j.
          ! n_4 is node in face f_j, but not in f_i.
          ! n_1 and n_3 are the two common nodes in this specific order.
          n_2 = node_diff(meshes(particle%mesh_type), f_i, f_j)
          n_4 = node_diff(meshes(particle%mesh_type), f_j, f_i)
          call nodes_same(meshes(particle%mesh_type), f_i, f_j, n_1, n_3)

          ! Compute current angle between faces including the sign of the angle.
          ! convex is an integer for convex/concave angles, convex = 1 means that the pair of faces is convex.
          angle_current = angle_normals(particle, f_i, f_j)

          if(dot_product(cross_product(particle%face(f_i)%normal, particle%face(f_j)%normal), &
            & (particle%node(n_1)%pos - particle%node(n_3)%pos)) .ge. 0) then
            convex = 1
          else
            convex = -1
          end if

          ! Compute force strength.
!           strength = lambda * (angle_current - convex * meshes(particle%mesh_type)%angle(f_i, n)) ! with spontaneous curvature
          strength = lambda * (angle_current) ! without spontaneous curvature

          ! Update forces.
          particle%node(n_1)%force_b = particle%node(n_1)%force_b + &
            & ((cross_product(particle%node(n_2)%pos - particle%node(n_3)%pos, normal_ji) * inv_area_i + &
            & cross_product(particle%node(n_3)%pos - particle%node(n_4)%pos, normal_ij) * inv_area_j) * strength)
          particle%node(n_2)%force_b = particle%node(n_2)%force_b + &
            & (cross_product(particle%node(n_3)%pos - particle%node(n_1)%pos, normal_ji) * inv_area_i * strength)
          particle%node(n_3)%force_b = particle%node(n_3)%force_b + &
            & ((cross_product(particle%node(n_4)%pos - particle%node(n_1)%pos, normal_ij) * inv_area_j + &
            & cross_product(particle%node(n_1)%pos - particle%node(n_2)%pos, normal_ji) * inv_area_i) * strength)
          particle%node(n_4)%force_b = particle%node(n_4)%force_b + &
            & (cross_product(particle%node(n_1)%pos - particle%node(n_3)%pos, normal_ij) * inv_area_j * strength)

          ! Update bending energy.
!           particle%erg_b = particle%erg_b + (lambda * (convex * angle_current - meshes(particle%mesh_type)%angle(f_i, n))**2)
          particle%erg_b = particle%erg_b + lambda * (angle_current)**2
        end if
      end do
    end do

    ! Update total energy.
    particle%erg_tot = particle%erg_tot + particle%erg_b
  end subroutine forces_bending

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute strain forces
  !>
  !> The strain forces are computed. It is based on the constitutive model of the membrane.
  !> If the strain and area moduli are quasi zero, the computation is skipped.
  !> The strain force is the sum of all face contributions (details in my dissertation, p. 43 and p. 145).
  !> TODO: at some point it would be nice to add functionality to switch the constitutive model.

  subroutine forces_strain(particle)
    type(extobj_struct), intent(inout) :: particle !< element of the particle array

    ! Declare variables.
    integer :: f_i ! face index
    integer :: n_1, n_2, n_3 ! node indices
    real(kind=rk) :: l1_eq, l2_eq, cos_eq, inv_sin_eq ! reference geometry
    real(kind=rk) :: l1_def, l2_def, phi_def, cos_def, sin_def ! deformed geometry
    real(kind=rk) :: b0, a1, b1 ! shape functions
    real(kind=rk) :: Dxx, Dyy, Dxy, Gxx, Gyy, Gxy ! displacement gradient tensor
    real(kind=rk) :: I1, I2 ! principal strains
    real(kind=rk) :: w, dw_dI1, dw_dI2 ! energy and its derivatives
    real(kind=rk) :: dI1_dGxx, dI1_dGyy, dI2_dGxx, dI2_dGyy, dI2_dGxy ! derivatives of the strain invariants
    real(kind=rk) :: dGxx_du1x, dGxy_du0x, dGxy_du1x, dGxy_du1y ! derivatives of the squared deformation tensor #1
    real(kind=rk) :: dGyy_du0x, dGyy_du0y, dGyy_du1x, dGyy_du1y ! derivatives of the squared deformation tensor #2
    real(kind=rk) :: force1x, force1y, force2x, force2y ! force components in common plane
    real(kind=rk), dimension(3) :: ex, ey, ez ! coordinate system
    real(kind=rk), dimension(3) :: force_1, force_2, force_3 ! final forces

    ! Skip computation if the bending modulus is small.
    if((particle%k_s .lt. 1d-15) .and. (particle%k_al .lt. 1d-15)) return

    ! Compute strain forces.
    do f_i = 1, particle%num_faces
      n_1 = meshes(particle%mesh_type)%neighbor_face_node(f_i, 1)
      n_2 = meshes(particle%mesh_type)%neighbor_face_node(f_i, 2)
      n_3 = meshes(particle%mesh_type)%neighbor_face_node(f_i, 3)

      ! Compute reference geometry.
      ! The angles phi and phi_eq shall be located at node n_2.
      ! The aligned edges l1 and l1_eq shall be between nodes n_2 and n_3.
      ! The edges l2 and l2_eq shall be between nodes n_2 and n_1.
      l1_eq = meshes(particle%mesh_type)%l1(f_i) * particle%radius
      l2_eq = meshes(particle%mesh_type)%l2(f_i) * particle%radius
      cos_eq = meshes(particle%mesh_type)%cos(f_i)
      inv_sin_eq = meshes(particle%mesh_type)%inv_sin(f_i)

      ! Compute deformed geometry.
      l1_def = norm(particle%node(n_3)%pos - particle%node(n_2)%pos)
      l2_def = norm(particle%node(n_1)%pos - particle%node(n_2)%pos)
      phi_def = angle_edges(particle, n_1, n_2, n_3)
      cos_def = cos(phi_def)
      sin_def = sin(phi_def)

      ! Set shape functions.
      ! NOTE: In the current version of the code, the shape functions are stored in memory instead of recomputing them.
      b0 = meshes(particle%mesh_type)%b0(f_i) * particle%radius
      a1 = meshes(particle%mesh_type)%a1(f_i) * particle%radius
      b1 = meshes(particle%mesh_type)%b1(f_i) * particle%radius

      ! Compute deformation gradient tensor.
      Dxx = l1_def / l1_eq
      Dyy = (l2_def / l2_eq) * sin_def * inv_sin_eq
      Dxy = ((l2_def / l2_eq) * cos_def - (l1_def / l1_eq) * cos_eq) * inv_sin_eq
      Gxx = Dxx**2
      Gyy = Dxy**2 + Dyy**2
      Gxy = Dxx * Dxy

      ! Compute principal stresses.
      I1 = Gxx + Gyy - 2.0
      I2 = Gxx * Gyy - Gxy**2 - 1.d0

      ! Use Skalak's constitutive law.
      w = (particle%k_s / 12.d0) * (I1**2 + 2.d0 * I1 - 2.d0 * I2) + (particle%k_al / 12.d0) * I2**2
      dw_dI1 = (particle%k_s / 6.d0) * (I1 + 1.d0)
      dw_dI2 = -particle%k_s / 6.d0 + (particle%k_al / 6.d0) * I2

      ! Compute derivatives of the strain invariants.
      dI1_dGxx = 1.d0
      dI1_dGyy = 1.d0
      dI2_dGxx = Gyy
      dI2_dGyy = Gxx
      dI2_dGxy = -2.d0 * Gxy

      ! Compute derivatives of the squared deformation tensor.
      dGxx_du1x = 2.d0 * a1 * Dxx
      dGxy_du0x = b0 * Dxx
      dGxy_du1x = a1 * Dxy + b1 * Dxx
      dGxy_du1y = a1 * Dyy
      dGyy_du0x = 2.d0 * b0 * Dxy
      dGyy_du0y = 2.d0 * b0 * Dyy
      dGyy_du1x = 2.d0 * b1 * Dxy
      dGyy_du1y = 2.d0 * b1 * Dyy

      ! Compute force components in common plane.
      force1x = dw_dI1 * (dI1_dGyy * dGyy_du0x) + dw_dI2 * (dI2_dGyy * dGyy_du0x + dI2_dGxy * dGxy_du0x)
      force1y = dw_dI1 * (dI1_dGyy * dGyy_du0y) + dw_dI2 * (dI2_dGyy * dGyy_du0y)
      force2x = dw_dI1 * (dI1_dGxx * dGxx_du1x + dI1_dGyy * dGyy_du1x) + &
              & dw_dI2 * (dI2_dGxx * dGxx_du1x + dI2_dGyy * dGyy_du1x + dI2_dGxy * dGxy_du1x)
      force2y = dw_dI1 * (dI1_dGyy * dGyy_du1y) + dw_dI2 * (dI2_dGyy * dGyy_du1y + dI2_dGxy * dGxy_du1y)

      ! Set coordinate system.
      ex = particle%node(n_3)%pos - particle%node(n_2)%pos
      ex = ex / norm(ex)
      ez = cross_product((particle%node(n_3)%pos - particle%node(n_2)%pos), (particle%node(n_1)%pos - particle%node(n_2)%pos))
      ez = ez / norm(ez)
      ey = cross_product(ez, ex)

      ! Compute final forces.
      force_1 = ex * force1x + ey * force1y
      force_2 = ex * force2x + ey * force2y
      force_3 = -(force_1 + force_2)

      ! Update forces.
      ! The required face area for the forces is already contained in the shape functions.
      particle%node(n_1)%force_s = particle%node(n_1)%force_s - force_1
      particle%node(n_2)%force_s = particle%node(n_2)%force_s - force_2
      particle%node(n_3)%force_s = particle%node(n_3)%force_s - force_3

      ! Update strain energy.
      particle%erg_s = particle%erg_s + (w * meshes(particle%mesh_type)%area(f_i) * (particle%radius)**2)
    end do

    ! Update total energy.
    particle%erg_tot = particle%erg_tot + particle%erg_s
  end subroutine forces_strain

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute friction forces
  !>
  !> Friction forces are computed. They are proportional to the node velocity.
  !> The friction forces are only used for particle growth in order to reduce the kinetic energy.

  subroutine forces_friction(particle, friction)
    type(extobj_struct), intent(inout) :: particle !< element of the particle array
    real(kind=rk), intent(in) :: friction !< friction coefficient 

    ! Declare variables.
    integer :: n_i ! node index

    ! Compute friction forces.
    do n_i = 1, particle%num_nodes
      particle%node(n_i)%force_s = particle%node(n_i)%force_s - particle%node(n_i)%vel * friction
    end do
  end subroutine forces_friction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  
  !> Compute surface tension forces
  !>
  !> The fluid-fluid interface is emulated by pulling the surface within a virtual interface plane.
  !> This approach is only reasonable for a single-fluid simulation.
  
#ifdef IBM_EMULATETENSION
  subroutine forces_interface_emulation(particle, surften, contactangle, width)
    type(extobj_struct), intent(inout) :: particle !< element of the particle array
    real(kind=rk), intent(in) :: surften !< surface tension
    real(kind=rk), intent(in) :: contactangle !< contact angle
    real(kind=rk), intent(in) :: width !< interface thickness
  
    ! Declare variables.
    integer :: f_i ! face index
    integer :: n_i, n_1, n_2, n_3 ! node indices
    real(kind=rk) :: area ! exposed area of face
    real(kind=rk) :: x_int ! interface location along x-axis (interface is in y-z plane)
    real(kind=rk) :: force ! force on face
    real(kind=rk) :: area_tot ! total interface area
    real(kind=rk) :: cosangle ! cosine of current contact angle
    real(kind=rk) :: coscontact ! cosine of desired contact angle
    real(kind=rk) :: force_contact ! force due to contact angle
    real(kind=rk), dimension(3) :: normal_av ! average normal vector
    real(kind=rk), dimension(3) :: direction ! force direction vector
    real(kind=rk), dimension(3) :: origin1, origin2, normal ! intersecting plane properties
   
    ! Set plane intersecting properties.
    x_int = 0.5 * tnx
    origin1 = (/x_int - 0.5 * width, 0.0d0, 0.0d0/)
    origin2 = (/x_int + 0.5 * width, 0.0d0, 0.0d0/)
    normal = (/1.0d0, 0.0d0, 0.0d0/)
    
    ! Compute average normal vector and total interface area.
    normal_av = (/0.0, 0.0, 0.0/)
    area_tot = 0.0
    
    ! Loop over all face elements
    do f_i = 1, particle%num_faces
      ! Identify nodes in face and their positions.
      n_1 = meshes(particle%mesh_type)%neighbor_face_node(f_i, 1)
      n_2 = meshes(particle%mesh_type)%neighbor_face_node(f_i, 2)
      n_3 = meshes(particle%mesh_type)%neighbor_face_node(f_i, 3)

      ! Compute disposed area of each face.
      area = compute_area_slice(particle, f_i, origin1, normal) - compute_area_slice(particle, f_i, origin2, normal)
      
      ! Update average normal vector and total interface area.
      area_tot = area_tot + area
      normal_av = normal_av + particle%face(f_i)%normal * area

      ! Compute force on face.
      ! Area / width is effective length of interface on face.
      force = surften * area / width

      ! Compute direction vector.
      direction = (/0.0d0, particle%face(f_i)%normal(2), particle%face(f_i)%normal(3)/)

      ! Update forces.
      particle%node(n_1)%force_b = particle%node(n_1)%force_b + direction * force / 3.0d0
      particle%node(n_2)%force_b = particle%node(n_2)%force_b + direction * force / 3.0d0
      particle%node(n_3)%force_b = particle%node(n_3)%force_b + direction * force / 3.0d0
    end do
    
    ! Normalise average normal vector.
    normal_av = normal_av / area_tot
    
    ! Compute contact angle force.
    cosangle = normal_av(1)
    coscontact = cos(contactangle * 3.14159265 / 180.0)
    force_contact = (cosangle - coscontact)
    
    do n_i = 1, particle%num_nodes
      particle%node(n_i)%force_b(1) = particle%node(n_i)%force_b(1) + force_contact / particle%num_nodes
    end do
  end subroutine forces_interface_emulation
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute constant node force
  !>
	!> Distributes the total force on a particle homogeneously over the particle mesh nodes taking
	!> into account the total area corresponding to a node. 

#ifdef IBM_DRAG
  subroutine forces_const_node(particle)
		type(extobj_struct), intent(inout) :: particle !< element of the particle array

		! Declare variables.
		integer :: n_i ! node index
		
		! When no force added skip entire loop
		if (particle%force_const(1) .eq. 0.d0 .and. particle%force_const(2) .eq. 0.d0 .and. particle%force_const(3) .eq. 0.d0) return

		! Normalise force with node area.
		do n_i = 1, particle%num_nodes 
				particle%node(n_i)%force_c(:) = particle%force_const(:) * particle%node(n_i)%area / particle%surface
		end do

	end subroutine forces_const_node
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute anchor forces
  !>
  !> Computes force due to Hookean springs connecting each mesh-node to its anchor point.

#ifdef IBM_FIXED
  subroutine forces_anchor(particle)
    type(extobj_struct), intent(inout) :: particle !< element of the particle array

    ! Declare variables.
    integer :: n_i ! node index
		
		particle%force_anchor(:) = 0.0d0
    ! Compute friction forces.
    do n_i = 1, particle%num_nodes
      particle%node(n_i)%force_anchor(:) = particle%k_anchor * (particle%node(n_i)%pos_anchor(:) - particle%node(n_i)%pos(:))

			! Add total anchor force to particle structure for output
			particle%force_anchor(:) = particle%force_anchor(:) + particle%node(n_i)%force_anchor(:)
    end do
  end subroutine forces_anchor
#endif 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute angular velocity
  !>
  !> The surface-averaged angular velocity of the particle is computed.
  !> It is assumed that the rotation takes place in the xz-plane, thus only the y-component is computed.
  !> Nodes which are not too close to the rotation axis are considered and weighted by their area.
  !> Nodes too close to the axis may lead to strong noise since the result is divided by the square of this distance.

  subroutine compute_angular_velocity(particle)
    type(extobj_struct), intent(inout) :: particle !< element of the particle array

    ! Declare variables.
    integer :: n_i ! node index
    real(kind=rk) :: dist ! node distance to the rotation axis
    real(kind=rk) :: surface ! weighting surface
    real(kind=rk), dimension(3) :: vel ! node velocity relative to center
    real(kind=rk), dimension(3) :: pos ! node position relative to center
    real(kind=rk), dimension(3) :: omega ! node angular velocity

    ! Reset averaged angular velocity and weighting surface.
    particle%angular_velocity(:) = 0.0d0
    surface = 0.0d0

    ! Compute angular velocity.
    do n_i = 1, particle%num_nodes
      ! Use node velocity relative to particle center and project onto xy-plane.
      vel(:) = particle%node(n_i)%vel(:) - particle%linear_velocity(:)
      vel(2) = 0.0d0

      ! Use distance of node to particle center projected onto xy-plane.
      pos(:) = particle%node(n_i)%pos(:) - particle%center(:)
      pos(2) = 0.0d0
      dist = norm(pos)

      ! Check distance of node from rotation axis.
      ! If it is too small, neglect its contribution to avoid strong noise.
      if(dist < 1.0d0) then
        cycle
      end if

      ! Compute angular velocity of node.
      omega = cross_product(pos, vel) / dist**2

      ! Update averaged angular velocity and weighting surface.
      particle%angular_velocity(2) = particle%angular_velocity(2) + omega(2) * particle%node(n_i)%area
      surface = surface + particle%node(n_i)%area
    end do

    ! Normalize angular velocity by weighting surface.
    particle%angular_velocity(2) = particle%angular_velocity(2) / surface
  end subroutine compute_angular_velocity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute particle volume inertia tensor
  !>
  !> The inertia tensor of the particle is computed.
  !> This is the instantaneous inertia tensor based on the particle volume, not particle surface.
  !> This is the common approach to find the inertia ellipsoid of the particle.
  !> The particle is treated as a volume with constant density.
  !> The computation is based on the faces, not the nodes.

  subroutine compute_inertia_tensor(particle)
    type(extobj_struct), intent(inout) :: particle !< element of the particle array

    ! Declare variables.
    integer :: f_i ! face index
    integer :: n_1, n_2, n_3 ! node indices
    real(kind=rk), dimension(3) :: fc ! face centroid
    real(kind=rk) :: sp ! scalar product of pv and normal vector
    real(kind=rk) :: red_area ! reduced face area

    ! Reset inertia tensor.
    particle%inertia_tensor(:, :) = 0.0d0

    ! Compute independent components of the inertia tensor.
    do f_i = 1, particle%num_faces
      n_1 = meshes(particle%mesh_type)%neighbor_face_node(f_i, 1)
      n_2 = meshes(particle%mesh_type)%neighbor_face_node(f_i, 2)
      n_3 = meshes(particle%mesh_type)%neighbor_face_node(f_i, 3)
      fc(:) = (particle%node(n_1)%pos(:) + particle%node(n_2)%pos(:) + particle%node(n_3)%pos(:)) / 3.0d0 - particle%center(:)
      sp = dot_product(fc, particle%face(f_i)%normal)
      red_area = particle%face(f_i)%area / 5.0d0
      particle%inertia_tensor(1, 1) = particle%inertia_tensor(1, 1) + (red_area * ((dot_product(fc, fc) * sp) - fc(1)**2 * sp))
      particle%inertia_tensor(2, 2) = particle%inertia_tensor(2, 2) + (red_area * ((dot_product(fc, fc) * sp) - fc(2)**2 * sp))
      particle%inertia_tensor(3, 3) = particle%inertia_tensor(3, 3) + (red_area * ((dot_product(fc, fc) * sp) - fc(3)**2 * sp))
      particle%inertia_tensor(1, 2) = particle%inertia_tensor(1, 2) + (red_area * (-fc(1) * fc(2) * sp))
      particle%inertia_tensor(1, 3) = particle%inertia_tensor(1, 3) + (red_area * (-fc(1) * fc(3) * sp))
      particle%inertia_tensor(2, 3) = particle%inertia_tensor(2, 3) + (red_area * (-fc(2) * fc(3) * sp))
    end do

    ! Compute dependent components of the symmetric tensor.
    particle%inertia_tensor(2, 1) = particle%inertia_tensor(1, 2)
    particle%inertia_tensor(3, 1) = particle%inertia_tensor(1, 3)
    particle%inertia_tensor(3, 2) = particle%inertia_tensor(2, 3)
  end subroutine compute_inertia_tensor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute particle inertia ellipsoid
  !>
  !> The inertia ellpsoid of the particle is computed from its inertia tensor.
  !> The semiaxes and orientation vectors are sorted (largest axis first).

  subroutine compute_inertia_ellipsoid(particle)
    type(extobj_struct), intent(inout) :: particle !< element of the particle array

    ! Declare variables.
    real(kind=rk), dimension(3, 3) :: input_matrix ! matrix to diagonalize
    real(kind=rk), dimension(3, 3) :: output_vectors ! eigenvectors
    real(kind=rk), dimension(3) :: output_values ! eigenvalues
    real(kind=rk), dimension(3) :: temp_vector ! temporary vector for swapping
    real(kind=rk) :: temp_value ! temporary value for swapping

    ! Initialize matrix with inertia tensor.
    input_matrix(:, :) = particle%inertia_tensor(:, :)

    ! Call helper function to diagonalize the matrix.
    call diag_symm_matrix(input_matrix, output_vectors, output_values)

    ! Sort eigenvalues (smaller first) and eigenvectors.
    if(output_values(2) < output_values(1) .and. output_values(2) < output_values(3)) then
      temp_value = output_values(1)
      temp_vector(:) = output_vectors(1, :)
      output_values(1) = output_values(2)
      output_vectors(1, :) = output_vectors(2, :)
      output_values(2) = temp_value
      output_vectors(2, :) = temp_vector(:)
    else if(output_values(3) < output_values(1) .and. output_values(3) < output_values(2)) then
      temp_value = output_values(1)
      temp_vector(:) = output_vectors(1, :)
      output_values(1) = output_values(3)
      output_vectors(1, :) = output_vectors(3, :)
      output_values(3) = temp_value
      output_vectors(3, :) = temp_vector(:)
    end if

    if(output_values(3) < output_values(2)) then
      temp_value = output_values(2)
      temp_vector(:) = output_vectors(2, :)
      output_values(2) = output_values(3)
      output_vectors(2, :) = output_vectors(3, :)
      output_values(3) = temp_value
      output_vectors(3, :) = temp_vector(:)
    end if

    ! Copy eigenvectors to particle properties.
    particle%inertia_ell_orient(:, :) = output_vectors(:, :)

    ! Find semiaxes of inertia ellipsoid.
    particle%inertia_ell_axes(1) = sqrt(2.5d0 * (output_values(2) + output_values(3) - output_values(1)) / particle%volume)
    particle%inertia_ell_axes(2) = sqrt(2.5d0 * (output_values(3) + output_values(1) - output_values(2)) / particle%volume)
    particle%inertia_ell_axes(3) = sqrt(2.5d0 * (output_values(1) + output_values(2) - output_values(3)) / particle%volume)
  end subroutine compute_inertia_ellipsoid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute particle stresslet
  !>
  !> The stresslet of the particle is computed.
  !> For this, all internal particle forces are used. They obey total momentum and angular momentum conservation.
  !> External forces (particle-particle or particle-wall interactions) are not taken into account here.

  subroutine compute_particle_stresslet(particle)
    type(extobj_struct), intent(inout) :: particle !< element of the particle array

    ! Declare variables.
    integer :: n_i ! node index
    real(kind=rk), dimension(3) :: node_pos ! node position
    real(kind=rk), dimension(3) :: node_force ! node force

    ! Reset particle stresslet.
    particle%stresslet_xx = 0.0d0
    particle%stresslet_xy = 0.0d0
    particle%stresslet_xz = 0.0d0
    particle%stresslet_yy = 0.0d0
    particle%stresslet_yz = 0.0d0
    particle%stresslet_zz = 0.0d0

    ! Compute particle stresslet.
    ! The old node position (before particle update) has to be used for this.
    do n_i = 1, particle%num_nodes
      node_pos(:) = particle%node(n_i)%pos_old(:)
      node_force(:) = particle%node(n_i)%force_tot(:) - particle%node(n_i)%force_int(:)
      particle%stresslet_xx = particle%stresslet_xx - node_pos(1) * node_force(1)
      particle%stresslet_xy = particle%stresslet_xy - node_pos(1) * node_force(2)
      particle%stresslet_xz = particle%stresslet_xz - node_pos(1) * node_force(3)
      particle%stresslet_yy = particle%stresslet_yy - node_pos(2) * node_force(2)
      particle%stresslet_yz = particle%stresslet_yz - node_pos(2) * node_force(3)
      particle%stresslet_zz = particle%stresslet_zz - node_pos(3) * node_force(3)
    end do

    ! Normalize stresslet by volume of the fluid.
    particle%stresslet_xx = particle%stresslet_xx / n_sites_fluid
    particle%stresslet_xy = particle%stresslet_xy / n_sites_fluid
    particle%stresslet_xz = particle%stresslet_xz / n_sites_fluid
    particle%stresslet_yy = particle%stresslet_yy / n_sites_fluid
    particle%stresslet_yz = particle%stresslet_yz / n_sites_fluid
    particle%stresslet_zz = particle%stresslet_zz / n_sites_fluid
  end subroutine compute_particle_stresslet

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute volume of particle separated by a plane
  !>
  !> Any plane separates a particle into two parts (where one part can have zero size).
  !> The plane is defined by an origin and its normal vector.
  !> The volume of the part of the particle in normal direction of the plane is computed.
  !> This function is useful for finding the volume of a particle between two parallel planes.
  !> This allows to compute the local volume fraction in a planar geometry.
  !> Finding the correct volume is basically a tedious case-by-case analysis running over all particle faces.

  function compute_volume_slice(particle, slice_origin, slice_normal)
    real(kind=rk) :: compute_volume_slice !< return value
    type(extobj_struct), intent(in) :: particle !< element of the particle array
    real(kind=rk), dimension(3), intent(in) :: slice_origin !< origin of the slicing plane
    real(kind=rk), dimension(3), intent(in) :: slice_normal !< normal of the slicing plane

    ! Declare variables.
    integer :: f_i ! face index
    integer :: n_1, n_2, n_3 ! node indices
    real(kind=rk) :: dist_1, dist_2, dist_3 ! node distances from slicing plane
    real(kind=rk), dimension(3) :: slice_norm ! normalized slice normal vector
    real(kind=rk), dimension(3) :: pos_1, pos_2, pos_3 ! node positions
    real(kind=rk) :: face_area ! area of a face
    real(kind=rk) :: face_area_eff ! effective area of a face (if intersected)
    real(kind=rk), dimension(3) :: face_normal ! normal vector of a face
    real(kind=rk), dimension(3) :: face_center ! center position of a face
    real(kind=rk), dimension(3) :: intsec_12, intsec_13, intsec_21, intsec_23, intsec_31, intsec_32 ! plane intersections
    real(kind=rk) :: volume_slice ! volume in normal direction of the plane

    ! Normalize normal vector to unit length.
    slice_norm = slice_normal / norm(slice_normal)

    ! Start with zero volume.
    volume_slice = 0.d0
    compute_volume_slice = 0.d0

    ! Run over all faces and compute the volume contributions.
    do f_i = 1, particle%num_faces
      ! Find nodes in face.
      n_1 = meshes(particle%mesh_type)%neighbor_face_node(f_i, 1)
      n_2 = meshes(particle%mesh_type)%neighbor_face_node(f_i, 2)
      n_3 = meshes(particle%mesh_type)%neighbor_face_node(f_i, 3)

      ! Alias node/face data for more compact code.
      pos_1 = particle%node(n_1)%pos
      pos_2 = particle%node(n_2)%pos
      pos_3 = particle%node(n_3)%pos
      face_area = particle%face(f_i)%area
      face_normal = particle%face(f_i)%normal

      ! Compute signed distance of each node from the plane.
      dist_1 = dot_product(pos_1 - slice_origin, slice_norm)
      dist_2 = dot_product(pos_2 - slice_origin, slice_norm)
      dist_3 = dot_product(pos_3 - slice_origin, slice_norm)

      ! 1) If all three nodes are in normal direction of the plane, the face volume is completely taken into account.
      if(dist_1 > 0.d0 .and. dist_2 > 0.d0 .and. dist_3 > 0.d0) then
        ! Compute centroid of the face.
        face_center = (pos_1 + pos_2 + pos_3) / 3.d0

        ! Update sliced volume.
        volume_slice = volume_slice + dot_product(face_center - slice_origin, face_normal) * face_area / 3.d0

      ! 2) If only one node is in normal direction of the plane, the part of the face volume in normal direction is taken.
      ! 2a) Only node 1 is in normal direction.
      else if(dist_1 > 0.d0 .and. dist_2 <= 0.d0 .and. dist_3 <= 0.d0) then
        ! Compute intersection points.
        intsec_12 = pos_1 + (pos_2 - pos_1) * dist_1 / (dist_1 - dist_2)
        intsec_13 = pos_1 + (pos_3 - pos_1) * dist_1 / (dist_1 - dist_3)

        ! Compute effective face area.
        face_area_eff = norm(cross_product(intsec_12 - pos_1, intsec_13 - pos_1)) / 2.d0

        ! Update sliced volume.
        volume_slice = volume_slice + dot_product(pos_1 - slice_origin, face_normal) * face_area_eff / 3.d0

      ! 2b) Only node 2 is in normal direction.
      else if(dist_1 <= 0.d0 .and. dist_2 > 0.d0 .and. dist_3 <= 0.d0) then
        ! Compute intersection points.
        intsec_21 = pos_2 + (pos_1 - pos_2) * dist_2 / (dist_1 - dist_2)
        intsec_23 = pos_2 + (pos_3 - pos_2) * dist_2 / (dist_2 - dist_3)

        ! Compute effective face area.
        face_area_eff = norm(cross_product(intsec_21 - pos_2, intsec_23 - pos_2)) / 2.d0

        ! Update sliced volume.
        volume_slice = volume_slice + dot_product(pos_2 - slice_origin, face_normal) * face_area_eff / 3.d0

      ! 2c) Only node 3 is in normal direction.
      else if(dist_1 <= 0.d0 .and. dist_2 <= 0.d0 .and. dist_3 > 0.d0) then
        ! Compute intersection points.
        intsec_31 = pos_3 + (pos_1 - pos_3) * dist_3 / (dist_1 - dist_3)
        intsec_32 = pos_3 + (pos_2 - pos_3) * dist_3 / (dist_2 - dist_3)

        ! Compute effective face area.
        face_area_eff = norm(cross_product(intsec_31 - pos_3, intsec_32 - pos_3)) / 2.d0

        ! Update sliced volume.
        volume_slice = volume_slice + dot_product(pos_3 - slice_origin, face_normal) * face_area_eff / 3.d0

      ! 3) If two nodes are in normal direction of the plane, the part of the face volume in normal direction is taken
      !    as complement to the part in negative normal direction.
      ! 3a) Only node 3 is not in normal direction.
      else if(dist_1 > 0.d0 .and. dist_2 > 0.d0 .and. dist_3 <= 0.d0) then
        ! Compute intersection points.
        intsec_31 = pos_3 + (pos_1 - pos_3) * dist_3 / (dist_1 - dist_3)
        intsec_32 = pos_3 + (pos_2 - pos_3) * dist_3 / (dist_2 - dist_3)

        ! Compute effective face area.
        face_area_eff = face_area - norm(cross_product(intsec_31 - pos_3, intsec_32 - pos_3)) / 2.d0

        ! Update sliced volume.
        volume_slice = volume_slice + dot_product(pos_3 - slice_origin, face_normal) * face_area_eff / 3.d0

      ! 3b) Only node 2 is not in normal direction.
      else if(dist_1 > 0.d0 .and. dist_2 <= 0.d0 .and. dist_3 > 0.d0) then
        ! Compute intersection points.
        intsec_21 = pos_2 + (pos_1 - pos_2) * dist_2 / (dist_1 - dist_2)
        intsec_23 = pos_2 + (pos_3 - pos_2) * dist_2 / (dist_2 - dist_3)

        ! Compute effective face area.
        face_area_eff = face_area - norm(cross_product(intsec_21 - pos_2, intsec_23 - pos_2)) / 2.d0

        ! Update sliced volume.
        volume_slice = volume_slice + dot_product(pos_2 - slice_origin, face_normal) * face_area_eff / 3.d0

      ! 3c) Only node 1 is not in normal direction.
      else if(dist_1 <= 0.d0 .and. dist_2 > 0.d0 .and. dist_3 > 0.d0) then
        ! Compute intersection points.
        intsec_12 = pos_1 + (pos_2 - pos_1) * dist_1 / (dist_1 - dist_2)
        intsec_13 = pos_1 + (pos_3 - pos_1) * dist_1 / (dist_1 - dist_3)

        ! Compute effective face area.
        face_area_eff = face_area - norm(cross_product(intsec_12 - pos_1, intsec_13 - pos_1)) / 2.d0

        ! Update sliced volume.
        volume_slice = volume_slice + dot_product(pos_1 - slice_origin, face_normal) * face_area_eff / 3.d0

      ! Else, all three nodes are not in normal direction, and there is no volume contribution.
      end if
    end do

    compute_volume_slice = volume_slice
  end function compute_volume_slice

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute area of face separated by a plane
  !>
  !> Any plane separates a face into two parts (where one part can have zero size).
  !> The plane is defined by an origin and its normal vector.
  !> The area of the part of the face in normal direction of the plane is computed.
  !> This function is useful for finding the area of a face between two parallel planes.
  !> Finding the correct area is basically a tedious case-by-case analysis running over all particle faces.

  function compute_area_slice(particle, f_i, slice_origin, slice_normal)
    real(kind=rk) :: compute_area_slice !< return value
    type(extobj_struct), intent(in) :: particle !< element of the particle array
    integer, intent (in) :: f_i !< index of face to be evaluated
    real(kind=rk), dimension(3), intent(in) :: slice_origin !< origin of the slicing plane
    real(kind=rk), dimension(3), intent(in) :: slice_normal !< normal of the slicing plane

    ! Declare variables.
    integer :: n_1, n_2, n_3 ! node indices
    real(kind=rk) :: dist_1, dist_2, dist_3 ! node distances from slicing plane
    real(kind=rk), dimension(3) :: slice_norm ! normalized slice normal vector
    real(kind=rk), dimension(3) :: pos_1, pos_2, pos_3 ! node positions
    real(kind=rk) :: face_area ! area of a face
    real(kind=rk), dimension(3) :: intsec_12, intsec_13, intsec_21, intsec_23, intsec_31, intsec_32 ! plane intersections
    real(kind=rk) :: surface_slice ! volume in normal direction of the plane

    ! Normalize normal vector to unit length.
    slice_norm = slice_normal / norm(slice_normal)

    ! Find nodes in face.
    n_1 = meshes(particle%mesh_type)%neighbor_face_node(f_i, 1)
    n_2 = meshes(particle%mesh_type)%neighbor_face_node(f_i, 2)
    n_3 = meshes(particle%mesh_type)%neighbor_face_node(f_i, 3)

    ! Alias node/face data for more compact code.
    pos_1 = particle%node(n_1)%pos
    pos_2 = particle%node(n_2)%pos
    pos_3 = particle%node(n_3)%pos
    face_area = particle%face(f_i)%area

    ! Compute signed distance of each node from the plane.
    dist_1 = dot_product(pos_1 - slice_origin, slice_norm)
    dist_2 = dot_product(pos_2 - slice_origin, slice_norm)
    dist_3 = dot_product(pos_3 - slice_origin, slice_norm)

    ! 1) If all three nodes are in normal direction of the plane, the face surface is completely taken into account.
    if(dist_1 > 0.d0 .and. dist_2 > 0.d0 .and. dist_3 > 0.d0) then
      ! Update sliced volume.
      surface_slice = face_area

    ! 2) If only one node is in normal direction of the plane, the part of the face volume in normal direction is taken.
    ! 2a) Only node 1 is in normal direction.
    else if(dist_1 > 0.d0 .and. dist_2 <= 0.d0 .and. dist_3 <= 0.d0) then
      ! Compute intersection points.
      intsec_12 = pos_1 + (pos_2 - pos_1) * dist_1 / (dist_1 - dist_2)
      intsec_13 = pos_1 + (pos_3 - pos_1) * dist_1 / (dist_1 - dist_3)

      ! Compute effective face area.
      surface_slice = norm(cross_product(intsec_12 - pos_1, intsec_13 - pos_1)) / 2.d0

    ! 2b) Only node 2 is in normal direction.
    else if(dist_1 <= 0.d0 .and. dist_2 > 0.d0 .and. dist_3 <= 0.d0) then
      ! Compute intersection points.
      intsec_21 = pos_2 + (pos_1 - pos_2) * dist_2 / (dist_1 - dist_2)
      intsec_23 = pos_2 + (pos_3 - pos_2) * dist_2 / (dist_2 - dist_3)

      ! Compute effective face area.
      surface_slice = norm(cross_product(intsec_21 - pos_2, intsec_23 - pos_2)) / 2.d0

    ! 2c) Only node 3 is in normal direction.
    else if(dist_1 <= 0.d0 .and. dist_2 <= 0.d0 .and. dist_3 > 0.d0) then
      ! Compute intersection points.
      intsec_31 = pos_3 + (pos_1 - pos_3) * dist_3 / (dist_1 - dist_3)
      intsec_32 = pos_3 + (pos_2 - pos_3) * dist_3 / (dist_2 - dist_3)

      ! Compute effective face area.
      surface_slice = norm(cross_product(intsec_31 - pos_3, intsec_32 - pos_3)) / 2.d0

    ! 3) If two nodes are in normal direction of the plane, the part of the face volume in normal direction is taken
    !    as complement to the part in negative normal direction.
    ! 3a) Only node 3 is not in normal direction.
    else if(dist_1 > 0.d0 .and. dist_2 > 0.d0 .and. dist_3 <= 0.d0) then
      ! Compute intersection points.
      intsec_31 = pos_3 + (pos_1 - pos_3) * dist_3 / (dist_1 - dist_3)
      intsec_32 = pos_3 + (pos_2 - pos_3) * dist_3 / (dist_2 - dist_3)

      ! Compute effective face area.
      surface_slice = face_area - norm(cross_product(intsec_31 - pos_3, intsec_32 - pos_3)) / 2.d0

    ! 3b) Only node 2 is not in normal direction.
    else if(dist_1 > 0.d0 .and. dist_2 <= 0.d0 .and. dist_3 > 0.d0) then
      ! Compute intersection points.
      intsec_21 = pos_2 + (pos_1 - pos_2) * dist_2 / (dist_1 - dist_2)
      intsec_23 = pos_2 + (pos_3 - pos_2) * dist_2 / (dist_2 - dist_3)

      ! Compute effective face area.
      surface_slice = face_area - norm(cross_product(intsec_21 - pos_2, intsec_23 - pos_2)) / 2.d0

    ! 3c) Only node 1 is not in normal direction.
    else if(dist_1 <= 0.d0 .and. dist_2 > 0.d0 .and. dist_3 > 0.d0) then
      ! Compute intersection points.
      intsec_12 = pos_1 + (pos_2 - pos_1) * dist_1 / (dist_1 - dist_2)
      intsec_13 = pos_1 + (pos_3 - pos_1) * dist_1 / (dist_1 - dist_3)

      ! Compute effective face area.
      surface_slice = face_area - norm(cross_product(intsec_12 - pos_1, intsec_13 - pos_1)) / 2.d0

    ! Else, all three nodes are not in normal direction, and there is no volume contribution.
    else
      surface_slice = 0.0d0
    end if

    compute_area_slice = surface_slice
  end function compute_area_slice
  
! endif IBM_PART
#endif

end module lextobj_module
