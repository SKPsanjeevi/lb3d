!> Lagrangian mesh module
!>
!> This module contains all structures, variables, and subroutines responsible
!> for initializing and maintaining the Lagrangian meshes.

#include "lbe.h"

module lmesh_module

#ifdef IBM_PART

  ! Include external modules.
  use lbe_globals_module, only : pi, rk
  use lbe_helper_module, only : cross_product, norm
  use lbe_log_module

  implicit none
  private
  public :: mesh_struct, meshes, MAX_NODE_NEIGHBORS, init_mesh_master, init_mesh_slave, &
            & nodes_same, allocate_mesh_memory, node_diff

  !> Definition of the mesh structure
  !>
  !> It contains all properties which are common for all particles of the same mesh,
  !> e.g., the number of nodes or the neighbor lists (if the node connectivity is fixed).

  type mesh_struct
    ! scalars
    integer mesh_type !< mesh index
    integer num_nodes !< number of nodes in mesh
    integer num_faces !< number of faces in mesh
    real(kind=rk) surface !< mesh surface area
    real(kind=rk) volume !< mesh volume

    ! arrays
    real(kind=rk), allocatable, dimension(:,:) :: pos !< node positions
    real(kind=rk), allocatable, dimension(:) :: area !< face areas
    real(kind=rk), allocatable, dimension(:,:) :: angle !< angles between faces
    integer, allocatable, dimension(:,:) :: neighbor_face_face !< neighbor list, face to face
    integer, allocatable, dimension(:,:) :: neighbor_node_node !< neighbor list, node to node
    integer, allocatable, dimension(:,:) :: neighbor_node_face !< neighbor list, node to face
    integer, allocatable, dimension(:,:) :: neighbor_face_node !< neighbor list, face to node
    real(kind=rk), allocatable, dimension(:) :: l1, l2, phi, cos, sin, inv_sin !< face reference geometry
    real(kind=rk), allocatable, dimension(:) :: b0, a1, b1, a2, b2 !< face shape functions (for FEM)
  end type mesh_struct

  ! Variable declarations
  ! These variables do not have to be saved in a checkpoint.
  type(mesh_struct), allocatable, dimension(:) :: meshes !< array containing all meshes required during simulation
  integer, parameter :: MAX_NODE_NEIGHBORS = 7 !< maximum number of neighbors a node can have

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize mesh (master)
  !>
  !> This subroutine is called once for each mesh which has to be set up.
  !> It is only called by rank 0 which is responsible for reading the mesh files.

  subroutine init_mesh_master(filename, mesh, m_type)
    character(len=*), intent(in) :: filename !< name of the mesh file
    type(mesh_struct), intent(inout) :: mesh !< element of the mesh array
    integer, intent(in) :: m_type !< mesh index

    ! Declare variables.
    character(len=200) :: message ! message

    ! Report progress.
    write(message, "('initializing mesh from file ', a, '...')") trim(filename)
    call log_msg(trim(message), .false.)

    ! Set mesh index
    mesh%mesh_type = m_type

    ! Call subroutines for mesh initialization.
    call read_mesh_from_file(filename, mesh) ! read mesh from file and allocate memory
    call find_node_face_neighbors(mesh) ! set up look-up tables for node and face neighbors
    call set_face_normals(mesh) ! compute normal vectors and have them point in outward direction
    call set_reference_values(mesh) ! compute refence values of shape functions and geometry
    call report_mesh_statistics(mesh) ! report mesh properties to the user
  end subroutine init_mesh_master

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize mesh (slave)
  !>
  !> This subroutine is called once for each mesh which has to be set up.
  !> It is only called by ranks other than 0 which shall not read the mesh files.

  subroutine init_mesh_slave(mesh, m_type)
    type(mesh_struct), intent(inout) :: mesh !< element of the mesh array
    integer, intent(in) :: m_type !< mesh index

    ! Set mesh index
    mesh%mesh_type = m_type

    ! Call subroutines for mesh initialization.
    call find_node_face_neighbors(mesh) ! set up look-up tables for node and face neighbors
    call set_face_normals(mesh) ! compute normal vectors and have them point in outward direction
    call set_reference_values(mesh) ! compute refence values of shape functions and geometry
    call report_mesh_statistics(mesh) ! report mesh properties to the user
  end subroutine init_mesh_slave

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> A mesh file is read, and a corresponding mesh is created in memory.
  !>
  !> The initial node positions are set, and the face membership is defined.
  !> NOTE: Right now, only closed triangulated 2D meshes with fixed node connectivity are supported.
  !> The file format corresponds to that of Gmsh which allows inspection of mesh files without running this code.

  subroutine read_mesh_from_file(filename, mesh)
    character(len=*), intent(in) :: filename !< name of the mesh file
    type(mesh_struct), intent(inout) :: mesh !< element of the mesh array

    ! Declare variables.
    character(len=200) :: message ! message
    integer, parameter :: init_file_unit = 13 ! input file unit
    character(len=1000) :: line_buffer ! buffer for reading lines
    integer :: int_buf ! buffer for reading integers
    integer :: i ! counter
    integer :: num_nodes, num_faces ! number of nodes and faces
    logical :: file_exist ! switch for file existence

    ! Report progress.
    write(message, "('reading mesh from file ', a, '...')") trim(filename)
    call log_msg(trim(message), .false.)

    ! Check for existence of file.
    inquire(file=filename, exist=file_exist)

    ! If it does not exist, switch to automatic positioning and give a warning.
    if(file_exist .eqv. .false.) then
      call error("mesh file not found, expected in directory ../msh, please check directory and file name")
    end if

    ! Open file.
    open(init_file_unit, file=filename, status='old', action='read')

    ! Skip 4 header lines.
    do i = 1, 4
      read(init_file_unit, '(a999)') line_buffer
    end do

    ! Read number of nodes.
    read(init_file_unit, '(i10)') num_nodes

    ! Allocate arrays.
    mesh%num_nodes = num_nodes
    allocate(mesh%pos(3, num_nodes))
    allocate(mesh%neighbor_node_node(num_nodes, MAX_NODE_NEIGHBORS))
    allocate(mesh%neighbor_node_face(num_nodes, MAX_NODE_NEIGHBORS))

    ! Read node data
    do i = 1, num_nodes
      read(init_file_unit, fmt=*) int_buf, mesh%pos(:, i)
    end do

    ! Skip 2 header lines.
    do i = 1, 2
      read(init_file_unit, '(a)') line_buffer
    end do

    ! Read number of faces.
    read(init_file_unit, '(i10)') num_faces

    ! Allocate arrays.
    mesh%num_faces = num_faces
    allocate(mesh%area(num_faces))
    allocate(mesh%angle(num_faces, 3))
    allocate(mesh%neighbor_face_face(num_faces, 3))
    allocate(mesh%neighbor_face_node(num_faces, 3))
    allocate(mesh%l1(num_faces), mesh%l2(num_faces), mesh%phi(num_faces))
    allocate(mesh%cos(num_faces), mesh%sin(num_faces), mesh%inv_sin(num_faces))
    allocate(mesh%b0(num_faces), mesh%a1(num_faces), mesh%b1(num_faces))
    allocate(mesh%a2(num_faces), mesh%b2(num_faces))

    ! Read face data.
    do i = 1, num_faces
      read(init_file_unit, fmt=*) int_buf, int_buf, int_buf, int_buf, int_buf, int_buf, mesh%neighbor_face_node(i, :)
    end do

    ! Close file.
    close(init_file_unit)
  end subroutine read_mesh_from_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Memory allocation for mesh structures
  !>
  !> A slave process (rank other than 0) does not read the mesh files.
  !> Since the memory allocation operations are performed during reading the files,
  !> another version for memory allocation has to be provided without reading the files.

  subroutine allocate_mesh_memory(mesh, num_nodes, num_faces)
    type(mesh_struct), intent(inout) :: mesh !< element of the mesh array
    integer, intent(in) :: num_nodes !< number of nodes for the mesh
    integer, intent(in) :: num_faces !< number of faces for the mesh

    ! Set number of nodes and faces.
    mesh%num_nodes = num_nodes
    mesh%num_faces = num_faces

    ! Allocate memory for mesh arrays.
    allocate(mesh%pos(3, num_nodes))
    allocate(mesh%neighbor_node_node(num_nodes, MAX_NODE_NEIGHBORS))
    allocate(mesh%neighbor_node_face(num_nodes, MAX_NODE_NEIGHBORS))
    allocate(mesh%area(num_faces))
    allocate(mesh%angle(num_faces, 3))
    allocate(mesh%neighbor_face_face(num_faces, 3))
    allocate(mesh%neighbor_face_node(num_faces, 3))
    allocate(mesh%l1(num_faces), mesh%l2(num_faces), mesh%phi(num_faces))
    allocate(mesh%cos(num_faces), mesh%sin(num_faces), mesh%inv_sin(num_faces))
    allocate(mesh%b0(num_faces), mesh%a1(num_faces), mesh%b1(num_faces))
    allocate(mesh%a2(num_faces), mesh%b2(num_faces))
  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Node and face neighbors are identified.
  !>
  !> Three look-up tables for neighboring nodes and faces are created for an existing mesh:
  !> - face to neighboring faces
  !> - node to neighboring faces
  !> - node to neighboring nodes
  !> NOTE: The LUT \c neighbor_face_node has already been set in \c read_mesh_from_file.

  subroutine find_node_face_neighbors(mesh)
    type(mesh_struct), intent(inout) :: mesh !< element of the mesh array

    ! Declare variables.
    integer :: f_i, f_j ! face indices
    integer :: n_i, n_1, n_2, n_3 ! node indices
    integer :: m, n ! counters
    integer, dimension(3 * MAX_NODE_NEIGHBORS) :: temp_array ! auxilliary array for finding neighbors

    ! Report progress.
    call log_msg("setting up look-up tables for nodes and faces...", .false.)

    ! Valid LUT elements are non-negative integers.
    ! Use a negative value for the LUT elements at the beginning.
    ! This way, it is easy to detect invalid or unused entries later on.
    mesh%neighbor_face_face = -1
    mesh%neighbor_node_face = -1
    mesh%neighbor_node_node = -1

    ! Start with the LUT for face->face.
    ! Run over all faces and find other faces sharing exactly two nodes.
    do f_i = 1, mesh%num_faces
      n_1 = mesh%neighbor_face_node(f_i, 1)
      n_2 = mesh%neighbor_face_node(f_i, 2)
      n_3 = mesh%neighbor_face_node(f_i, 3)

      ! There are three possibilities for two common nodes.
      ! Each one uniquely corresponds to another neighboring face.
      do f_j = 1, mesh%num_faces
        if(f_j .ne. f_i) then
          if(is_node_in_face(mesh, n_1, f_j) .and. is_node_in_face(mesh, n_2, f_j)) then
            mesh%neighbor_face_face(f_i, 1) = f_j
          end if

          if(is_node_in_face(mesh, n_2, f_j) .and. is_node_in_face(mesh, n_3, f_j)) then
            mesh%neighbor_face_face(f_i, 2) = f_j
          end if

          if(is_node_in_face(mesh, n_3, f_j) .and. is_node_in_face(mesh, n_1, f_j)) then
            mesh%neighbor_face_face(f_i, 3) = f_j
          end if
        end if
      end do
    end do

    ! Continue with the LUT for node->face.
    ! Run over all nodes and find other nodes being member of the same face.
    do n_i = 1, mesh%num_nodes
      do f_i = 1, mesh%num_faces
        if(is_node_in_face(mesh, n_i, f_i)) then
          n = 1

          do while(mesh%neighbor_node_face(n_i, n) .ne. -1)
            n = n + 1
          end do

          mesh%neighbor_node_face(n_i, n) = f_i
        end if
      end do
    end do

    ! Conclude with the LUT for node->node.
    ! A temporary array is used to simplify the search for neighboring nodes.
    ! This array carries 'candidates' for a next neighbor which are checked in a second sweep.
    temp_array = -1

    do n_i = 1, mesh%num_nodes
      ! Set up the temporary array.
      do n = 1, MAX_NODE_NEIGHBORS
        f_i = mesh%neighbor_node_face(n_i, n)

        if(f_i .ne. -1) then
          temp_array(3 * (n - 1) + 1) = mesh%neighbor_face_node(f_i, 1)
          temp_array(3 * (n - 1) + 2) = mesh%neighbor_face_node(f_i, 2)
          temp_array(3 * (n - 1) + 3) = mesh%neighbor_face_node(f_i, 3)
        end if
      end do

      ! Filter candidates by avoiding double pairs.
      n = 1

      do m = 1, 3 * MAX_NODE_NEIGHBORS
        if((temp_array(m) .ne. -1) .and.&
          &(temp_array(m) .ne. n_i) .and.&
          &(.not.is_already_in_list(mesh, temp_array(m), n_i))) then
          mesh%neighbor_node_node(n_i, n) = temp_array(m)
          n = n + 1
        end if
      end do
    end do
  end subroutine find_node_face_neighbors

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Face normal vectors are oriented.
  !>
  !> By definition, normal vectors of closed membranes have to point in outward direction.
  !> This can be insured by storing the nodes for a face in a right-handed order.
  !> The order of the nodes in the look-up table \c neighbor_face_node is changed in such a way that normals point outwards.

  subroutine set_face_normals(mesh)
    type(mesh_struct), intent(inout) :: mesh !< element of the mesh array

    ! Declare variables.
    logical, dimension(mesh%num_faces) :: is_checked ! auxilliary array for storing check state of faces
    real(kind=rk), dimension(3) :: center ! temporary center position
    real(kind=rk), dimension(3) :: centroid ! center of mesh face
    real(kind=rk), dimension(3) :: normal, normal_i, normal_j ! normal vectors of faces
    real(kind=rk), dimension(3) :: pos_1, pos_2, pos_3 ! node positions
    integer :: n ! counter
    integer :: temp_node ! temporary node index
    integer :: n_i, n_1, n_2, n_3 ! node indices
    integer :: f_i, f_j, f_1, f_2, f_3 ! face indices

    ! Report progress.
    call log_msg("setting face normal vectors...", .false.)

    ! Start with assumption that no face has been checked.
    is_checked = .false.

    !!! Compute temporary center position of the mesh

    center = 0

    do n_i = 1, mesh%num_nodes
      center = center + mesh%pos(:,n_i)
    end do

    center = center / mesh%num_nodes

    !!! Start with first face (face index 1)
    ! It is assumed that the first face has a centroid pointing outwards.

    ! Identify nodes in first face
    n_1 = mesh%neighbor_face_node(1, 1)
    n_2 = mesh%neighbor_face_node(1, 2)
    n_3 = mesh%neighbor_face_node(1, 3)

    ! Compute face centroid and normal vector
    pos_1 = mesh%pos(:, n_1)
    pos_2 = mesh%pos(:, n_2)
    pos_3 = mesh%pos(:, n_3)
    centroid = (pos_1 + pos_2 + pos_3) / 3 - center
    normal = cross_product(pos_1 - pos_2, pos_3 - pos_2)

    ! Swap node indices if the normal vector points in different direction than the face center
    if(dot_product(centroid, normal) .lt. 0) then
      temp_node = mesh%neighbor_face_node(1, 1)
      mesh%neighbor_face_node(1, 1) = mesh%neighbor_face_node(1, 3)
      mesh%neighbor_face_node(1, 3) = temp_node
    end if

    is_checked(1) = .true.

    !!! Successively check all neighbors by running over all faces.
    ! For each face, the angle with respect to each of its three neighbors is computed.
    ! If the angle is too large, the order of nodes is changed so that the normal vector points outwards.
    ! After checking, the status of the face is changed to is_checked = true.

    do while(not_all_checked(mesh, is_checked))
      do f_i = 1, mesh%num_faces
        if(is_checked(f_i)) then
          ! Identify neighboring faces
          f_1 = mesh%neighbor_face_face(f_i, 1)
          f_2 = mesh%neighbor_face_face(f_i, 2)
          f_3 = mesh%neighbor_face_face(f_i, 3)

          ! Check normal vectors of pairs of faces
          if(.not.is_checked(f_1)) then
            call check_normals(mesh, f_1, f_i)
            is_checked(f_1) = .true.
          end if

          if(.not.is_checked(f_2)) then
            call check_normals(mesh, f_2, f_i)
            is_checked(f_2) = .true.
          end if

          if(.not.is_checked(f_3)) then
            call check_normals(mesh, f_3, f_i)
            is_checked(f_3) = .true.
          end if
        end if
      end do
    end do

    ! Check whether all angles are small.
    ! If at least one pair of neighboring faces is found with a large angle between their normals,
    ! something has gone wrong.

    do f_i = 1, mesh%num_faces
      do n = 1, 3
        ! Find neighboring face
        f_j = mesh%neighbor_face_face(f_i, n)

        ! Check angle between neigboring normals
        if(dot_product(normal_vector(mesh, f_i), normal_vector(mesh, f_j)) .lt. 0) then
          call error("Problem with mesh normal vector setup. Could not set all normal vectors. Check mesh file.")
        end if
      end do
    end do
  end subroutine set_face_normals

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Reference values for the undeformed mesh are computed.
  !>
  !> Some reference values (edge lengths, face angles) are required for force computation at each time step,
  !> but they do not change during the simulation.
  !> It is more efficient to keep these values in memory instead of recomputing them.
  !> NOTE: The efficiency advantage was obvious in a serial C++ code. It may be different here and could be tested again.

  subroutine set_reference_values(mesh)
    type(mesh_struct), intent(inout) :: mesh !< element of the mesh array

    ! Declare variables
    integer :: f_i, f_j ! face indices
    integer :: n_i, n_1, n_2, n_3 ! node indices
    integer :: i ! counter
    real(kind=rk), dimension(3) :: pos_1, pos_2, pos_3 ! node positions
    real(kind=rk), dimension(3) :: edge_23, edge_21 ! face edges
    real(kind=rk), dimension(3) :: center ! temporary center position
    real(kind=rk), dimension(3) :: normal ! face normal
    real(kind=rk), dimension(3) :: centroid ! face centroid
    real(kind=rk) :: angle_normals ! normal angle
    real(kind=rk) :: cosine ! cosine of normal angle
    real(kind=rk) :: angle_sign ! sign of angle

    ! Report progress.
    call log_msg("setting reference values...", .false.)

    !!! Compute temporary center position of the mesh

    center = 0

    do n_i = 1, mesh%num_nodes
      center = center + mesh%pos(:,n_i)
    end do

    center = center / mesh%num_nodes

    !!! Compute geometry and shape functions

    mesh%surface = 0
    mesh%volume = 0

    do f_i = 1, mesh%num_faces
      ! Identify nodes in the face
      n_1 = mesh%neighbor_face_node(f_i, 1)
      n_2 = mesh%neighbor_face_node(f_i, 2)
      n_3 = mesh%neighbor_face_node(f_i, 3)

      ! Compute edges
      pos_1 = mesh%pos(:,n_1)
      pos_2 = mesh%pos(:,n_2)
      pos_3 = mesh%pos(:,n_3)
      edge_23 = pos_3 - pos_2
      edge_21 = pos_1 - pos_2

      ! Compute normal vector and face centroid
      normal = cross_product(edge_21, edge_23)
      centroid = (pos_1 + pos_2 + pos_3) / 3 - center

      ! compute face area and update surface and volume
      mesh%area(f_i) = norm(normal) / 2
      mesh%surface = mesh%surface + mesh%area(f_i)
      mesh%volume = mesh%volume + (dot_product(normal, centroid) / 6)

      ! Compute edges
      pos_1 = mesh%pos(:,n_1)
      pos_2 = mesh%pos(:,n_2)
      pos_3 = mesh%pos(:,n_3)
      edge_23 = pos_3 - pos_2
      edge_21 = pos_1 - pos_2

      ! Set geometry
      mesh%l1(f_i) = norm(edge_23)
      mesh%l2(f_i) = norm(edge_21)
      mesh%phi(f_i) = acos(dot_product(edge_23, edge_21) / (norm(edge_23) * norm(edge_21)))
      mesh%cos(f_i) = cos(mesh%phi(f_i))
      mesh%sin(f_i) = sin(mesh%phi(f_i))
      mesh%inv_sin(f_i) = 1 / mesh%sin(f_i)

      ! Set shape functions
      mesh%b0(f_i) = mesh%l1(f_i) / 2
      mesh%a1(f_i) = -mesh%l2(f_i) * mesh%sin(f_i) / 2
      mesh%b1(f_i) = (mesh%l2(f_i) * mesh%cos(f_i) - mesh%l1(f_i)) / 2
      mesh%a2(f_i) = mesh%l2(f_i) * mesh%sin(f_i) / 2
      mesh%b2(f_i) = -mesh%l2(f_i) * mesh%cos(f_i) / 2

      ! Compute angles between faces
      do i = 1, 3
        ! Set neighboring face
        f_j = mesh%neighbor_face_face(f_i, i)

        ! Identify common nodes of both faces
        call nodes_same(mesh, f_i, f_j, n_1, n_2)
        pos_1 = mesh%pos(:,n_1)
        pos_2 = mesh%pos(:,n_2)

        ! Compute cosine of normal angle.
        cosine = dot_product(normal_vector(mesh, f_i), normal_vector(mesh, f_j))

        ! Compute angle from cosine.
        if(cosine .ge. 0.999999d0) then
          angle_normals = 0
        else if(cosine .le. -0.999999d0) then
          angle_normals = pi
        else
          angle_normals = acos(cosine)
        end if

        ! Compute sign of angle (can be concave or convex, has to be tested explicitly)
        angle_sign = dot_product(cross_product(normal_vector(mesh, f_i), normal_vector(mesh, f_j)), pos_1 - pos_2)

        if(angle_sign .gt. 0) then
          mesh%angle(f_i, i) = angle_normals
        else if(angle_sign .lt. 0) then
          mesh%angle(f_i, i) = -angle_normals
        else
          mesh%angle(f_i, i) = 0
        end if
      end do
    end do
  end subroutine set_reference_values

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Report mesh properties and statistics.
  !>
  !> After the mesh has been completely set up, its properties are reported to the user.
  !> The report contains the number of nodes and faces and the averages and variances of edgle lengths and face angles,
  !> as well as some other data which is not relevant for the further simulation.
  !> This subroutine is intented as a pure report tool which is useful for debugging.

  subroutine report_mesh_statistics(mesh)
    type(mesh_struct), intent(inout) :: mesh !< element of the mesh array

    ! Declare variables.
    integer :: f_i ! face index
    integer :: n_i, n_j ! node indices
    real(kind=rk) :: radius_eff ! effective radius
    real(kind=rk) :: volume_red ! reduced volume
    real(kind=rk) :: area_min, area_max, area_mean ! minimum, maximum, and mean areas
    real(kind=rk) :: area_var ! area variance
    real(kind=rk) :: angle_min, angle_max, angle_mean ! minimum, maximum, and mean angles
    real(kind=rk) :: angle_mag ! angle magnitude
    real(kind=rk) :: angle_var ! angle variance
    integer :: angle_num ! number of considered angles
    integer :: i ! counter
    integer :: num_min, num_max, num_curr ! minimum, maximum, and current number of neighboring nodes
    real(kind=rk) :: length_min, length_max, length_mean, length_curr ! minimum, maximum, mean, and current edge lengths
    real(kind=rk) :: length_var ! edge length variance
    integer :: length_num ! number of edges
    real(kind=rk) angle_edge_min, angle_edge_max, angle_edge_mean ! minimum, maximum, and mean edge angles
    real(kind=rk) angle_edge_mag ! edge angle magnitude
    real(kind=rk) angle_edge_var ! edge angle variance
    real(kind=rk), dimension(3) :: pos_1, pos_2, pos_3 ! node positions
    real(kind=rk), dimension(3) :: edge_1, edge_2 ! edges for edge angle computation
    character(len=200) :: message ! message string

    ! Compute effective radius and reduced volume.
    radius_eff = sqrt(mesh%surface / (4 * pi));
    volume_red = 3 * mesh%volume / (4 * pi * radius_eff**3);

    !!! Compute minimum, maximum, and mean areas.

    area_min = 10000 ! use some sufficiently large initial value
    area_max = 0
    area_mean = 0

    do f_i = 1, mesh%num_faces
      area_mean = area_mean + mesh%area(f_i)

      if(mesh%area(f_i) .lt. area_min) then
        area_min = mesh%area(f_i)
      end if

      if(mesh%area(f_i) .gt. area_max) then
        area_max = mesh%area(f_i)
      end if
    end do

    area_mean = area_mean / mesh%num_faces

    !!! Compute area variance.

    area_var = 0

    do f_i = 1, mesh%num_faces
      area_var = area_var + (area_mean - mesh%area(f_i))**2
    end do

    area_var = sqrt(area_var / mesh%num_faces) / area_mean;

    !!! Compute minimum, maximum, and mean angles between normals.

    angle_min = pi ! start with a sufficiently large value
    angle_max = 0
    angle_mean = 0
    angle_num = 0;

    do f_i = 1, mesh%num_faces
      do i = 1, 3
        angle_mag = abs(mesh%angle(f_i, i))
        angle_mean = angle_mean + angle_mag

        if(angle_mag .lt. angle_min) then
          angle_min = angle_mag
        end if

        if(angle_mag .gt. angle_max) then
          angle_max = angle_mag
        end if

        angle_num = angle_num + 1
      end do
    end do

    angle_mean = angle_mean / angle_num

    !!! Compute normal angle variance.

    angle_var = 0

    do f_i = 1, mesh%num_faces
      do i = 1, 3
        angle_var = angle_var + (angle_mean - mesh%angle(f_i, i))**2
      end do
    end do

    angle_var = sqrt(angle_var / angle_num) / angle_mean

    !!! Compute minimum and maximum edge angles and edge angle variance

    angle_edge_min = pi ! start with sufficiently large value
    angle_edge_mean = pi / 3 ! mean angle must be 60°
    angle_edge_max = 0
    angle_edge_var = 0

    do f_i = 1, mesh%num_faces
      do i = 1, 3
        pos_1 = mesh%pos(:, mesh%neighbor_face_node(f_i, mod(i    , 3) + 1))
        pos_2 = mesh%pos(:, mesh%neighbor_face_node(f_i, mod(i + 1, 3) + 1))
        pos_3 = mesh%pos(:, mesh%neighbor_face_node(f_i, mod(i + 2, 3) + 1))
        edge_1 = pos_1 - pos_2
        edge_2 = pos_3 - pos_2
        angle_edge_mag = abs(acos(dot_product(edge_1, edge_2) / (norm(edge_1) * norm(edge_1))))

        if(angle_edge_mag .lt. angle_edge_min) then
          angle_edge_min = angle_edge_mag
        end if

        if(angle_edge_mag .gt. angle_edge_max) then
          angle_edge_max = angle_edge_mag
        end if

        angle_edge_var = angle_edge_var + (angle_edge_mag - angle_edge_mean)**2
      end do
    end do

    angle_edge_var = sqrt(angle_edge_var / (3 * mesh%num_faces)) / angle_edge_mean

    !!! Compute minimum, maximum and mean edges

    length_min = 10000 ! start with sufficiently large value
    length_max = 0
    length_mean = 0
    length_num = 0

    do n_i = 1, mesh%num_nodes
      do n_j = n_i + 1, mesh%num_nodes - 1
        if(are_neighbor_nodes(mesh, n_i, n_j)) then
          length_curr = norm(mesh%pos(:, n_i) - mesh%pos(:, n_j))
          length_mean = length_mean + length_curr

          if(length_curr .le. length_min) then
            length_min = length_curr
          end if

          if(length_curr .ge. length_max) then
            length_max = length_curr
          end if

          length_num = length_num + 1
        end if
      end do
    end do

    length_mean = length_mean / length_num

    !!! Compute edge variance

    length_var = 0

    do n_i = 1, mesh%num_nodes
      do n_j = n_i + 1, mesh%num_nodes - 1
        if(are_neighbor_nodes(mesh, n_i, n_j)) then
          length_curr = norm(mesh%pos(:, n_i) - mesh%pos(:, n_j))
          length_var = length_var + (length_curr - length_mean)**2
        end if
      end do
    end do

    length_var = sqrt(length_var / length_num) / length_mean

    !!! Compute number of neighboring nodes

    num_min = MAX_NODE_NEIGHBORS ! start with sufficiently large number
    num_max = 0
    num_curr = 0

    do n_i = 1, mesh%num_nodes
      i = 1

      do while(mesh%neighbor_node_face(n_i, i) .ne. -1)
        num_curr = i
        i = i + 1
      end do

      if(num_curr .lt. num_min) then
        num_min = num_curr
      end if

      if(num_curr .gt. num_max) then
        num_max = num_curr
      end if
    end do

    !!! Report results.
    call log_msg("mesh properties:", .false.)
    write(message, "('  volume: ', f0.3)") mesh%volume
    call log_msg(trim(message), .false.)
    write(message, "('  surface: ', f0.3)") mesh%surface
    call log_msg(trim(message), .false.)
    write(message, "('  eff. radius: ', f0.3)") radius_eff
    call log_msg(trim(message), .false.)
    write(message, "('  red. volume: ', f0.3)") volume_red
    call log_msg(trim(message), .false.)
    write(message, "('  node neighbors: [', i0, ', ', i0, ']')") num_min, num_max
    call log_msg(trim(message), .false.)
    write(message, "('  edges [min, max]: [', f0.3, ', ', f0.3, ']')") length_min, length_max
    call log_msg(trim(message), .false.)
    write(message, "('  edges (mean, rel. var.): ', f0.3, ', ', f0.3)") length_mean, length_var
    call log_msg(trim(message), .false.)
    write(message, "('  areas [min, max]: [', f0.3, ', ', f0.3, ']')") area_min, area_max
    call log_msg(trim(message), .false.)
    write(message, "('  areas (mean, rel. var.): ', f0.3, ', ', f0.3)") area_mean, area_var
    call log_msg(trim(message), .false.)
    write(message, "('  normal angles [min, max]: [', f0.3, ', ', f0.3, ']')") angle_min, angle_max
    call log_msg(trim(message), .false.)
    write(message, "('  normal angles (mean, rel. var.): ', f0.3, ', ', f0.3)") angle_mean, angle_var
    call log_msg(trim(message), .false.)
    write(message, "('  edge angles [min, max]: [', f0.3, ', ', f0.3, ']')") angle_edge_min, angle_edge_max
    call log_msg(trim(message), .false.)
    write(message, "('  edge angles (mean, rel. var.): ', f0.3, ', ', f0.3)") angle_edge_mean, angle_edge_var
    call log_msg(trim(message), .false.)
  end subroutine report_mesh_statistics

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> The angle between two neighboring faces is checked.
  !>
  !> If the angle is larger than 90°, the orientation of the first face has to be changed.
  !> For doing so, the first and third elements in the look-up table neighbor_face_node are swapped,
  !> which changes the sign of the corresponding normal.

  subroutine check_normals(mesh, face_check, face_ref)
    type(mesh_struct), intent(inout) :: mesh !< element of the mesh array
    integer, intent(in) :: face_check !< face to be checked
    integer, intent(in) :: face_ref !< reference face

    ! Declare variables.
    integer :: temp_face ! temporary face index

    ! Swap first and third nodes if necessary.
    if(dot_product(normal_vector(mesh, face_check), normal_vector(mesh, face_ref)) .lt. 0) then
      temp_face = mesh%neighbor_face_node(face_check, 1)
      mesh%neighbor_face_node(face_check, 1) = mesh%neighbor_face_node(face_check, 3)
      mesh%neighbor_face_node(face_check, 3) = temp_face
    end if
  end subroutine check_normals

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute unit normal vector of a given face in a given mesh
  !>
  !> The orientation of the normal vector depends on the ordering of the nodes belonging to the face.
  !> The normal vector is normalized to unit length.

  function normal_vector(mesh, face)
    real(kind=rk), dimension(3) :: normal_vector !< return value
    type(mesh_struct), intent(in) :: mesh !< element of the mesh array
    integer, intent(in) :: face !< face index

    ! Declare variables.
    integer :: n_1, n_2, n_3 ! nodes in face
    real(kind=rk), dimension(3) :: pos_1, pos_2, pos_3 ! node positions

    ! Identify nodes in face.
    n_1 = mesh%neighbor_face_node(face, 1)
    n_2 = mesh%neighbor_face_node(face, 2)
    n_3 = mesh%neighbor_face_node(face, 3)

    ! Compute normal vector.
    pos_1 = mesh%pos(:, n_1)
    pos_2 = mesh%pos(:, n_2)
    pos_3 = mesh%pos(:, n_3)
    normal_vector = cross_product(pos_1 - pos_2, pos_3 - pos_2) ! note the order of the nodes
    normal_vector = normal_vector / norm(normal_vector) ! normalize to unit length
  end function normal_vector

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Find both common nodes of two neighboring faces
  !>
  !> Both returned nodes are member of both faces.
  !> The user has to make sure that both faces are neighbors.
  !> The order of the nodes carries information about the relative alignment of the faces.
  !> NOTE: The cycling order of nodes is \c n_1 -> \c n_2 -> \c n_3 -> \c n_1.

  subroutine nodes_same(mesh, face_i, face_j, node_1, node_2)
    type(mesh_struct), intent(in) :: mesh !< element of the mesh array
    integer, intent(in) :: face_i !< first input face index
    integer, intent(in) :: face_j !< second input face index
    integer, intent(out) :: node_1 !< first common node
    integer, intent(out) :: node_2 !< second common node

    ! Declare variables.
    integer :: n_1, n_2, n_3 ! nodes in first face

    ! Find nodes in first face.
    n_1 = mesh%neighbor_face_node(face_i, 1)
    n_2 = mesh%neighbor_face_node(face_i, 2)
    n_3 = mesh%neighbor_face_node(face_i, 3)

    ! Find common nodes
    if(n_1 .eq. node_diff(mesh, face_i, face_j)) then
      node_1 = n_2
      node_2 = n_3
    else if(n_2 .eq. node_diff(mesh, face_i, face_j)) then
      node_1 = n_3
      node_2 = n_1
    else
      node_1 = n_1
      node_2 = n_2
    end if
  end subroutine nodes_same

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> The node in one given face but not in the other given face is returned.
  !>
  !> NOTE: The user has to make sure that both faces are neighbors.

  function node_diff(mesh, face_1, face_2)
    integer :: node_diff !< return value
    type(mesh_struct), intent(in) :: mesh !< element of the mesh array
    integer, intent(in) :: face_1 !< first input face index
    integer, intent(in) :: face_2 !< second input face index

    ! Declare variables.
    integer :: n_1, n_2, n_3 ! nodes in first face

    ! Find nodes in first face.
    n_1 = mesh%neighbor_face_node(face_1, 1)
    n_2 = mesh%neighbor_face_node(face_1, 2)
    n_3 = mesh%neighbor_face_node(face_1, 3)

    ! Identify single node
    if(.not.is_node_in_face(mesh, n_1, face_2)) then
      node_diff = n_1
    else if(.not.is_node_in_face(mesh, n_2, face_2)) then
      node_diff = n_2
    else
      node_diff = n_3
    end if
  end function node_diff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Checks whether a given node belongs to a given face in a given mesh.
  !>
  !> The return value is \c true if the node belongs to the face and \c false otherwise.

  function is_node_in_face(mesh, node, face)
    logical :: is_node_in_face !< return value
    type(mesh_struct), intent(in) :: mesh !< element of the mesh array
    integer, intent(in) :: node !< node to check for membership
    integer, intent(in) :: face !< face to check for membership

    if((node .eq. mesh%neighbor_face_node(face, 1)) .or.&
      &(node .eq. mesh%neighbor_face_node(face, 2)) .or.&
      &(node .eq. mesh%neighbor_face_node(face, 3))) then
      is_node_in_face = .true.
    else
      is_node_in_face = .false.
    end if
  end function is_node_in_face

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check whether two given nodes are neighbors
  !>
  !> Two nodes are neighbors when the first is member of the look-up table \c neighbor_node_node of the second.

  function are_neighbor_nodes(mesh, node_1, node_2)
    logical :: are_neighbor_nodes !< return value
    type(mesh_struct), intent(in) :: mesh !< element of the mesh array
    integer, intent(in) :: node_1 !< first node index
    integer, intent(in) :: node_2 !< second node index

    ! Declare variables.
    integer :: i ! counter

    ! Start with false assumption
    are_neighbor_nodes = .false.

    ! Check for neighborhood
    do i = 1, MAX_NODE_NEIGHBORS
      if(node_1 .eq. mesh%neighbor_node_node(node_2, i)) then
        are_neighbor_nodes = .true.
      end if
    end do
  end function are_neighbor_nodes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Checks whether a given node has already been found as neighbor for another given node in a given mesh.
  !>
  !> The return value is \c true if the node to be checked has been found before and \c false otherwise.

  function is_already_in_list(mesh, node_check, node_ref)
    logical :: is_already_in_list !< return value
    type(mesh_struct), intent(in) :: mesh !< element of the mesh array
    integer, intent(in) :: node_check !< node to check for being in list
    integer, intent(in) :: node_ref !< reference node

    ! Declare variables.
    integer :: i ! counter

    is_already_in_list = .false. ! assume false at the beginning

    do i = 1, MAX_NODE_NEIGHBORS
      if(mesh%neighbor_node_node(node_ref, i) .eq. node_check) then
        is_already_in_list = .true.
      end if
    end do
  end function is_already_in_list

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Asks whether all faces have already been checked.
  !>
  !> It is asked whether all faces in the mesh have already been tested for their normal vector orientation.
  !> The function returns true if at least one face has not been checked yet and false otherwise.

  function not_all_checked(mesh, is_checked)
    logical :: not_all_checked !< return value
    type(mesh_struct), intent(in) :: mesh !< element of the mesh array
    logical, dimension(:), intent(in) :: is_checked !< auxilliary array storing the check status of the faces

    ! Declare variables.
    integer :: f_i ! face index

    ! Set standard return value.
    not_all_checked = .false.

    ! If at least one face has not been checked, return true.
    do f_i = 1, mesh%num_faces
      if(.not.is_checked(f_i)) then
        not_all_checked = .true.
        return
      end if
    end do
  end function not_all_checked

! endif IBM_PART
#endif

end module lmesh_module
