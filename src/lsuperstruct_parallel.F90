!> Lagrangian superstructure parallelization module
!>
!> This module contains all parallelization functionality required for the Lagrangian superstructure.

#include "lbe.h"

module lsuperstruct_parallel_module

#ifdef IBM_PART

  ! Include external modules.
  use lbe_globals_module, only: myrankc
  use lbe_log_module
  use lbe_parallel_module, only : tnx, tny, tnz, nnprocs, cdims, ccoords, comm_cart
  use lbe_parms_module, only : nx, ny, nz
  use lmesh_module, only : meshes
  use lsuperstruct_data_module
  use lsuperstruct_helper_module, only : add_particle, remove_particle
  use lextobj_module, only : particles, part_ind, particle_node, particle_face, extobj_struct, init_particle, &
                              & set_center_position, update_face_areas
  use lmesh_module, only : meshes

  implicit none
  include 'mpif.h'
  private
  public :: rank_responsible_for_pos, create_datatypes_pre, create_datatypes_post, distribute_nodes, collect_nodes, &
              & exchange_particles, exchange_velocity_halo, exchange_rock_grad_halo
#ifdef IBM_INDEXFIELD
  public :: exchange_index_halo
#endif

  ! variable declarations

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Rank responsible for position
  !>
  !> The process rank responsible for a position vector is returned.
  !> This process corresponds to the subdomain for which the position vector is physical.

  function rank_responsible_for_pos(pos)
    integer :: rank_responsible_for_pos !< return value
    real(kind=rk), dimension(3), intent(in) :: pos !< input position vector

    ! Declare variables.
    integer, dimension(3) :: coords_pos ! process position along the axes
    integer :: dummy_rank
    integer :: ierror ! MPI error code

    ! Identify process position along the axes.
    coords_pos(1) = floor((pos(1) - 0.5d0) / nx)
    coords_pos(2) = floor((pos(2) - 0.5d0) / ny)
    coords_pos(3) = floor((pos(3) - 0.5d0) / nz)

    ! Get rank of the process.
    rank_responsible_for_pos = -1 ! start with invalid value
    call MPI_Cart_rank(comm_cart, coords_pos, dummy_rank, ierror)
    rank_responsible_for_pos = dummy_rank
  end function rank_responsible_for_pos

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Create derived datatypes
  !>
  !> All subroutine calls for the creation of the 'pre' individual derived datatype are called.
  !> 'pre' means that these datatypes can be defined before the meshes are created.

  subroutine create_datatypes_pre()
    call create_datatype_superstruct_node()
    call create_datatype_lextobj_node()
    call create_datatype_lextobj_face()
    call create_datatype_particle_init()
    call create_datatype_static_particle()

    ! Create datatypes for sending velocities.
    call build_velocity_chunk_mpitype((/nx - IBM_HALO + 1, nx/),       (/1 - IBM_HALO, ny + IBM_HALO/), &
        & (/1 - IBM_HALO, nz + IBM_HALO/), datatype_send_vel_xpos)
    call build_velocity_chunk_mpitype((/1, 1 + IBM_HALO - 1/),         (/1 - IBM_HALO, ny + IBM_HALO/), &
        & (/1 - IBM_HALO, nz + IBM_HALO/), datatype_send_vel_xneg)
    call build_velocity_chunk_mpitype((/1 - IBM_HALO, nx + IBM_HALO/), (/ny - IBM_HALO + 1, ny/),       &
        & (/1 - IBM_HALO, nz + IBM_HALO/), datatype_send_vel_ypos)
    call build_velocity_chunk_mpitype((/1 - IBM_HALO, nx + IBM_HALO/), (/1, 1 + IBM_HALO - 1/),         &
        & (/1 - IBM_HALO, nz + IBM_HALO/), datatype_send_vel_yneg)
    call build_velocity_chunk_mpitype((/1 - IBM_HALO, nx + IBM_HALO/), (/1 - IBM_HALO, ny + IBM_HALO/), &
        & (/nz - IBM_HALO + 1, nz/),       datatype_send_vel_zpos)
    call build_velocity_chunk_mpitype((/1 - IBM_HALO, nx + IBM_HALO/), (/1 - IBM_HALO, ny + IBM_HALO/), &
        & (/1, 1 + IBM_HALO - 1/),         datatype_send_vel_zneg)
    call build_velocity_chunk_mpitype((/nx + 1, nx + IBM_HALO/),       (/1 - IBM_HALO, ny + IBM_HALO/), &
        & (/1 - IBM_HALO, nz + IBM_HALO/), datatype_recv_vel_xpos)
    call build_velocity_chunk_mpitype((/1 - IBM_HALO, 0/),             (/1 - IBM_HALO, ny + IBM_HALO/), &
        & (/1 - IBM_HALO, nz + IBM_HALO/), datatype_recv_vel_xneg)
    call build_velocity_chunk_mpitype((/1 - IBM_HALO, nx + IBM_HALO/), (/ny + 1, ny + IBM_HALO/),       &
        & (/1 - IBM_HALO, nz + IBM_HALO/), datatype_recv_vel_ypos)
    call build_velocity_chunk_mpitype((/1 - IBM_HALO, nx + IBM_HALO/), (/1 - IBM_HALO, 0/),             &
        & (/1 - IBM_HALO, nz + IBM_HALO/), datatype_recv_vel_yneg)
    call build_velocity_chunk_mpitype((/1 - IBM_HALO, nx + IBM_HALO/), (/1 - IBM_HALO, ny + IBM_HALO/), &
        & (/nz + 1, nz + IBM_HALO/),       datatype_recv_vel_zpos)
    call build_velocity_chunk_mpitype((/1 - IBM_HALO, nx + IBM_HALO/), (/1 - IBM_HALO, ny + IBM_HALO/), &
        & (/1 - IBM_HALO, 0/),             datatype_recv_vel_zneg)

    ! Create datatypes for sending rock gradient.
    call build_rock_grad_chunk_mpitype((/nx - IBM_HALO + 1, nx/),       (/1 - IBM_HALO, ny + IBM_HALO/), &
        &(/1 - IBM_HALO, nz + IBM_HALO/), datatype_send_rock_grad_xpos)
    call build_rock_grad_chunk_mpitype((/1, 1 + IBM_HALO - 1/),         (/1 - IBM_HALO, ny + IBM_HALO/), &
        & (/1 - IBM_HALO, nz + IBM_HALO/), datatype_send_rock_grad_xneg)
    call build_rock_grad_chunk_mpitype((/1 - IBM_HALO, nx + IBM_HALO/), (/ny - IBM_HALO + 1, ny/),       &
        & (/1 - IBM_HALO, nz + IBM_HALO/), datatype_send_rock_grad_ypos)
    call build_rock_grad_chunk_mpitype((/1 - IBM_HALO, nx + IBM_HALO/), (/1, 1 + IBM_HALO - 1/),         &
        & (/1 - IBM_HALO, nz + IBM_HALO/), datatype_send_rock_grad_yneg)
    call build_rock_grad_chunk_mpitype((/1 - IBM_HALO, nx + IBM_HALO/), (/1 - IBM_HALO, ny + IBM_HALO/), &
        & (/nz - IBM_HALO + 1, nz/),       datatype_send_rock_grad_zpos)
    call build_rock_grad_chunk_mpitype((/1 - IBM_HALO, nx + IBM_HALO/), (/1 - IBM_HALO, ny + IBM_HALO/), &
        & (/1, 1 + IBM_HALO - 1/),         datatype_send_rock_grad_zneg)
    call build_rock_grad_chunk_mpitype((/nx + 1, nx + IBM_HALO/),       (/1 - IBM_HALO, ny + IBM_HALO/), &
        & (/1 - IBM_HALO, nz + IBM_HALO/), datatype_recv_rock_grad_xpos)
    call build_rock_grad_chunk_mpitype((/1 - IBM_HALO, 0/),             (/1 - IBM_HALO, ny + IBM_HALO/), &
        & (/1 - IBM_HALO, nz + IBM_HALO/), datatype_recv_rock_grad_xneg)
    call build_rock_grad_chunk_mpitype((/1 - IBM_HALO, nx + IBM_HALO/), (/ny + 1, ny + IBM_HALO/),       &
        & (/1 - IBM_HALO, nz + IBM_HALO/), datatype_recv_rock_grad_ypos)
    call build_rock_grad_chunk_mpitype((/1 - IBM_HALO, nx + IBM_HALO/), (/1 - IBM_HALO, 0/),             &
        & (/1 - IBM_HALO, nz + IBM_HALO/), datatype_recv_rock_grad_yneg)
    call build_rock_grad_chunk_mpitype((/1 - IBM_HALO, nx + IBM_HALO/), (/1 - IBM_HALO, ny + IBM_HALO/), &
        & (/nz + 1, nz + IBM_HALO/),       datatype_recv_rock_grad_zpos)
    call build_rock_grad_chunk_mpitype((/1 - IBM_HALO, nx + IBM_HALO/), (/1 - IBM_HALO, ny + IBM_HALO/), &
        & (/1 - IBM_HALO, 0/),             datatype_recv_rock_grad_zneg)

    ! Create datatypes for sending index field.
#ifdef IBM_INDEXFIELD
    call build_index_chunk_mpitype((/nx - IBM_HALO + 1, nx/),       (/1 - IBM_HALO, ny + IBM_HALO/), &
        &(/1 - IBM_HALO, nz + IBM_HALO/), datatype_send_index_xpos)
    call build_index_chunk_mpitype((/1, 1 + IBM_HALO - 1/),         (/1 - IBM_HALO, ny + IBM_HALO/), &
        & (/1 - IBM_HALO, nz + IBM_HALO/), datatype_send_index_xneg)
    call build_index_chunk_mpitype((/1 - IBM_HALO, nx + IBM_HALO/), (/ny - IBM_HALO + 1, ny/),       &
        & (/1 - IBM_HALO, nz + IBM_HALO/), datatype_send_index_ypos)
    call build_index_chunk_mpitype((/1 - IBM_HALO, nx + IBM_HALO/), (/1, 1 + IBM_HALO - 1/),         &
        & (/1 - IBM_HALO, nz + IBM_HALO/), datatype_send_index_yneg)
    call build_index_chunk_mpitype((/1 - IBM_HALO, nx + IBM_HALO/), (/1 - IBM_HALO, ny + IBM_HALO/), &
        & (/nz - IBM_HALO + 1, nz/),       datatype_send_index_zpos)
    call build_index_chunk_mpitype((/1 - IBM_HALO, nx + IBM_HALO/), (/1 - IBM_HALO, ny + IBM_HALO/), &
        & (/1, 1 + IBM_HALO - 1/),         datatype_send_index_zneg)
    call build_index_chunk_mpitype((/nx + 1, nx + IBM_HALO/),       (/1 - IBM_HALO, ny + IBM_HALO/), &
        & (/1 - IBM_HALO, nz + IBM_HALO/), datatype_recv_index_xpos)
    call build_index_chunk_mpitype((/1 - IBM_HALO, 0/),             (/1 - IBM_HALO, ny + IBM_HALO/), &
        & (/1 - IBM_HALO, nz + IBM_HALO/), datatype_recv_index_xneg)
    call build_index_chunk_mpitype((/1 - IBM_HALO, nx + IBM_HALO/), (/ny + 1, ny + IBM_HALO/),       &
        & (/1 - IBM_HALO, nz + IBM_HALO/), datatype_recv_index_ypos)
    call build_index_chunk_mpitype((/1 - IBM_HALO, nx + IBM_HALO/), (/1 - IBM_HALO, 0/),             &
        & (/1 - IBM_HALO, nz + IBM_HALO/), datatype_recv_index_yneg)
    call build_index_chunk_mpitype((/1 - IBM_HALO, nx + IBM_HALO/), (/1 - IBM_HALO, ny + IBM_HALO/), &
        & (/nz + 1, nz + IBM_HALO/),       datatype_recv_index_zpos)
    call build_index_chunk_mpitype((/1 - IBM_HALO, nx + IBM_HALO/), (/1 - IBM_HALO, ny + IBM_HALO/), &
        & (/1 - IBM_HALO, 0/),             datatype_recv_index_zneg)
#endif
  end subroutine create_datatypes_pre

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Create derived datatypes
  !>
  !> All subroutine calls for the creation of the 'post' individual derived datatype are called.
  !> 'post' means that these datatypes can be defined only after the meshes are created.

  subroutine create_datatypes_post()
    call create_datatype_lextobj()
  end subroutine create_datatypes_post

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Create derived datatype for superstructure nodes
  !>
  !> The datatype for the superstructure node is created.
  !> It is used to distribute/collect nodes from/to particles to/from the superstructure.

  subroutine create_datatype_superstruct_node()
    ! Declare variables.
    integer :: ierror ! MPI error code
    type(superstruct_node) :: temp_node ! example superstruct_node so that MPI knows about its memory layout
    integer, dimension(2) :: old_type
    integer(kind=MPI_ADDRESS_KIND), dimension(2) :: offset
    integer(kind=MPI_ADDRESS_KIND) :: base
    integer(kind=MPI_ADDRESS_KIND), dimension(2) :: address
    integer, dimension(2) :: blocklen

    ! Create derived datatype for communication of particle nodes.
    ! This datatype carries along all information required for proper node handling:
    ! 1) 8 integers (particle index, node index, 3 send switches, 3 halo switches)
    ! 2) 12 reals (3 position, 3 velocity, 3 force, 3 interaction force)
    old_type = (/MPI_INTEGER, MPI_REAL8/)
    call MPI_Get_address(temp_node, base, ierror)
    call MPI_Get_address(temp_node%particle_gl, address(1), ierror)
    call MPI_Get_address(temp_node%pos(1), address(2), ierror)
#ifndef IBM_FIXED
    blocklen = (/10, 16/)
#else
    blocklen = (/10, 20/) ! 3 more REALS for position anchor
#endif
    offset = address(:) - base
    call MPI_Type_create_struct(2, blocklen, offset, old_type, MPI_NODE_DATA, ierror)
    call MPI_Type_commit(MPI_NODE_DATA, ierror)
  end subroutine create_datatype_superstruct_node

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Create derived datatype for Lagrangian extended object nodes
  !>
  !> The datatype for the Lagrangian extended object node is created.
  !> It is used to create the datatype for a total particle later on.

  subroutine create_datatype_lextobj_node()
    ! Declare variables.
    integer :: ierror ! MPI error code
    type(particle_node) :: temp_node ! example particle_node so that MPI knows about its memory layout
    integer, dimension(1) :: old_type
    integer(kind=MPI_ADDRESS_KIND), dimension(1) :: offset
    integer(kind=MPI_ADDRESS_KIND) :: base
    integer(kind=MPI_ADDRESS_KIND), dimension(1) :: address
    integer, dimension(1) :: blocklen

    ! Create derived datatype for communication of particle nodes.
    ! This datatype carries along all information required for proper node handling: 28 reals
    old_type = (/MPI_REAL8/)
    call MPI_Get_address(temp_node, base, ierror)
    call MPI_Get_address(temp_node%area, address(1), ierror)
#ifdef IBM_DRAG
		blocklen = (/35/)
#elif IBM_FIXED
		blocklen = (/36/)
#else
		blocklen = (/33/)
#endif
    offset = address(:) - base
    call MPI_Type_create_struct(1, blocklen, offset, old_type, MPI_LAGRANGE_NODE_DATA, ierror)
    call MPI_Type_commit(MPI_LAGRANGE_NODE_DATA, ierror)
  end subroutine create_datatype_lextobj_node

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Create derived datatype for Lagrangian extended object faces
  !>
  !> The datatype for the Lagrangian extended object face is created.
  !> It is used to create the datatype for a total particle later on.

  subroutine create_datatype_lextobj_face()
    ! Declare variables.
    integer :: ierror ! MPI error code
    type(particle_face) :: temp_face ! example particle_face so that MPI knows about its memory layout
    integer, dimension(2) :: old_type
    integer(kind=MPI_ADDRESS_KIND), dimension(2) :: offset
    integer(kind=MPI_ADDRESS_KIND) :: base
    integer(kind=MPI_ADDRESS_KIND), dimension(2) :: address
    integer, dimension(2) :: blocklen

    ! Create derived datatype for communication of particle faces.
    ! This datatype carries along all information required for proper face handling: 4 reals
    old_type = (/MPI_REAL8, MPI_REAL8/)
    call MPI_Get_address(temp_face, base, ierror)
    call MPI_Get_address(temp_face%area, address(1), ierror)
    call MPI_Get_address(temp_face%normal, address(2), ierror)
    blocklen = (/1, 3/)
    offset = address(:) - base
    call MPI_Type_create_struct(2, blocklen, offset, old_type, MPI_LAGRANGE_FACE_DATA, ierror)
    call MPI_Type_commit(MPI_LAGRANGE_FACE_DATA, ierror)
  end subroutine create_datatype_lextobj_face

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Create derived datatype for initialization of particle configuration
  !>
  !> The datatype for the initialization of the particle configuration is created.

  subroutine create_datatype_particle_init()
    ! Declare variables.
    integer :: ierror ! MPI error code
    type(init_particle_info) :: temp_particle ! example init_particle_info so that MPI knows about its memory layout
    integer, dimension(2) :: old_type
    integer(kind=MPI_ADDRESS_KIND), dimension(2) :: offset
    integer(kind=MPI_ADDRESS_KIND) :: base
    integer(kind=MPI_ADDRESS_KIND), dimension(2) :: address
    integer, dimension(2) :: blocklen

    ! Create derived datatype for communication of initial particle configuration.
    old_type = (/MPI_INTEGER, MPI_REAL8/)
    call MPI_Get_address(temp_particle, base, ierror)
    call MPI_Get_address(temp_particle%particle_gl, address(1), ierror)
    call MPI_Get_address(temp_particle%radius, address(2), ierror)
#ifdef IBM_DRAG
    blocklen = (/2, 16/)
#else
    blocklen = (/2, 13/)
#endif
    offset = address(:) - base
    call MPI_Type_create_struct(2, blocklen, offset, old_type, MPI_PARTICLE_INIT, ierror)
    call MPI_Type_commit(MPI_PARTICLE_INIT, ierror)
  end subroutine create_datatype_particle_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Create derived datatype for Lagrangian extended objects
  !>
  !> The datatype for the Lagrangian extended object is created.
  !> It is used to communicate particles between processes.

  subroutine create_datatype_lextobj()
    ! Declare variables.
    integer :: ierror ! MPI error code
    type(extobj_struct) :: temp_particle ! example extobj_struct so that MPI knows about its memory layout
    integer, dimension(2) :: old_type
    integer(kind=MPI_ADDRESS_KIND), dimension(2) :: offset
    integer(kind=MPI_ADDRESS_KIND) :: base
    integer(kind=MPI_ADDRESS_KIND), dimension(2) :: address
    integer, dimension(2) :: blocklen

    ! Create derived datatype for communication of static particle data.
    old_type = (/MPI_INTEGER, MPI_REAL8/)
    call MPI_Get_address(temp_particle, base, ierror)
    call MPI_Get_address(temp_particle%num_nodes, address(1), ierror)
    call MPI_Get_address(temp_particle%radius_0, address(2), ierror)
#ifdef IBM_DRAG
    blocklen = (/8, 71/) ! number of integers and reals in the static part of the particle
#else
    blocklen = (/8, 68/) ! number of integers and reals in the static part of the particle
#endif
    offset = address(:) - base
    call MPI_Type_create_struct(2, blocklen, offset, old_type, MPI_LAGRANGE_PARTICLE_DATA, ierror)
    call MPI_Type_commit(MPI_LAGRANGE_PARTICLE_DATA, ierror)
  end subroutine create_datatype_lextobj

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Create derived datatype for dumping static particle data
  !>
  !> The datatype for dumping the static particle data is created.
  !> It is used to gather particle data before it is written to the disk.

  subroutine create_datatype_static_particle()
    ! Declare variables.
    integer :: ierror ! MPI error code
    type(dump_static_particle_info) :: temp_particle ! example dump_static_particle_info so that MPI knows about its memory layout
    integer, dimension(2) :: old_type
    integer(kind=MPI_ADDRESS_KIND), dimension(2) :: offset
    integer(kind=MPI_ADDRESS_KIND) :: base
    integer(kind=MPI_ADDRESS_KIND), dimension(2) :: address
    integer, dimension(2) :: blocklen

    ! Create derived datatype for communication of static particle data.
    old_type = (/MPI_INTEGER, MPI_REAL8/)
    call MPI_Get_address(temp_particle, base, ierror)
    call MPI_Get_address(temp_particle%particle_index_gl, address(1), ierror)
    call MPI_Get_address(temp_particle%volume, address(2), ierror)
    blocklen = (/1, 47/) ! number of integers and reals
    offset = address(:) - base
    call MPI_Type_create_struct(2, blocklen, offset, old_type, MPI_STATIC_PARTICLE_DATA, ierror)
    call MPI_Type_commit(MPI_STATIC_PARTICLE_DATA, ierror)
  end subroutine create_datatype_static_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Create derived datatype for sending lattice velocity
  !>
  !> This subroutine builds a custom MPI datatype that represents the part of the passed array specified by
  !> the coordinate intervals ext_x(2), ext_y(2), ext_z(2) and stores it in datatype.

  subroutine build_velocity_chunk_mpitype(ext_x, ext_y, ext_z, datatype_out)
    integer,intent(in) :: ext_x(2), ext_y(2), ext_z(2) !< chunk ranges for each dimension
    integer,intent(out) :: datatype_out !< MPI derived datatype to be built

    ! Declare variables.
    integer :: xrow, xyplane, xyzchunk ! temporary MPI datatypes
    integer :: cnt, blocklength, ierror, lengths(1) ! MPI variables
    integer(kind=MPI_ADDRESS_KIND) addr1, addr2, stride, base, offset, displs(1) ! addresses

    ! MPI datatype for columns along the x-axis
    cnt = 1 + ext_x(2) - ext_x(1)
    blocklength = 3 ! three contiguous reals for the velocity components
    call MPI_Get_Address(vel_phys(1, 1, 1, 1), addr1, ierror)
    call MPI_Get_Address(vel_phys(1, 2, 1, 1), addr2, ierror)
    stride = addr2 - addr1
    call MPI_Type_create_hvector(cnt, blocklength, stride, MPI_REAL8, xrow, ierror)

    ! MPI datatype for slices in the xy-plane
    cnt = 1 + ext_y(2) - ext_y(1)
    blocklength = 1
    call MPI_Get_Address(vel_phys(1, 1, 1, 1), addr1, ierror)
    call MPI_Get_Address(vel_phys(1, 1, 2, 1), addr2, ierror)
    stride = addr2 - addr1
    call MPI_Type_create_hvector(cnt, blocklength, stride, xrow, xyplane, ierror)

    ! MPI datatype for whole chunk in xyz-volume
    cnt = 1 + ext_z(2) - ext_z(1)
    blocklength = 1
    call MPI_Get_Address(vel_phys(1, 1, 1, 1), addr1, ierror)
    call MPI_Get_Address(vel_phys(1, 1, 1, 2), addr2, ierror)
    stride = addr2 - addr1
    call MPI_Type_create_hvector(cnt, blocklength, stride, xyplane, xyzchunk, ierror)

    ! Position of the beginning of the chunk relative to the beginning of vel_phys
    call MPI_Get_Address(vel_phys, base, ierror)
    call MPI_Get_Address(vel_phys(1, ext_x(1), ext_y(1), ext_z(1)), offset, ierror)

    ! Shift the memory relative to the base of vel_array
    cnt = 1
    lengths = (/1/)
    displs = (/offset - base/)
    call MPI_Type_create_hindexed(cnt, lengths, displs, xyzchunk, datatype_out, ierror)
    call MPI_Type_commit(datatype_out, ierror)
  end subroutine build_velocity_chunk_mpitype

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Create derived datatype for sending lattice rock gradient
  !>
  !> This subroutine builds a custom MPI datatype that represents the part of the passed array specified by
  !> the coordinate intervals ext_x(2), ext_y(2), ext_z(2) and stores it in datatype.

  subroutine build_rock_grad_chunk_mpitype(ext_x, ext_y, ext_z, datatype_out)
    integer,intent(in) :: ext_x(2), ext_y(2), ext_z(2) !< chunk ranges for each dimension
    integer,intent(out) :: datatype_out !< MPI derived datatype to be built

    ! Declare variables.
    integer :: xrow, xyplane, xyzchunk ! temporary MPI datatypes
    integer :: cnt, blocklength, ierror, lengths(1) ! MPI variables
    integer(kind=MPI_ADDRESS_KIND) addr1, addr2, stride, base, offset, displs(1) ! addresses

    ! MPI datatype for columns along the x-axis
    cnt = 1 + ext_x(2) - ext_x(1)
    blocklength = 1 ! one integer for the rock gradient
    call MPI_Get_Address(rock_grad(1, 1, 1), addr1, ierror)
    call MPI_Get_Address(rock_grad(2, 1, 1), addr2, ierror)
    stride = addr2 - addr1
    call MPI_Type_create_hvector(cnt, blocklength, stride, MPI_INTEGER, xrow, ierror)

    ! MPI datatype for slices in the xy-plane
    cnt = 1 + ext_y(2) - ext_y(1)
    blocklength = 1
    call MPI_Get_Address(rock_grad(1, 1, 1), addr1, ierror)
    call MPI_Get_Address(rock_grad(1, 2, 1), addr2, ierror)
    stride = addr2 - addr1
    call MPI_Type_create_hvector(cnt, blocklength, stride, xrow, xyplane, ierror)

    ! MPI datatype for whole chunk in xyz-volume
    cnt = 1 + ext_z(2) - ext_z(1)
    blocklength = 1
    call MPI_Get_Address(rock_grad(1, 1, 1), addr1, ierror)
    call MPI_Get_Address(rock_grad(1, 1, 2), addr2, ierror)
    stride = addr2 - addr1
    call MPI_Type_create_hvector(cnt, blocklength, stride, xyplane, xyzchunk, ierror)

    ! Position of the beginning of the chunk relative to the beginning of rock_gradient
    call MPI_Get_Address(rock_grad, base, ierror)
    call MPI_Get_Address(rock_grad(ext_x(1), ext_y(1), ext_z(1)), offset, ierror)

    ! Shift the memory relative to the base of vel_array
    cnt = 1
    lengths = (/1/)
    displs = (/offset - base/)
    call MPI_Type_create_hindexed(cnt, lengths, displs, xyzchunk, datatype_out, ierror)
    call MPI_Type_commit(datatype_out, ierror)
  end subroutine build_rock_grad_chunk_mpitype

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Create derived datatype for sending IBM index field
  !>
  !> This subroutine builds a custom MPI datatype that represents the part of the passed array specified by
  !> the coordinate intervals ext_x(2), ext_y(2), ext_z(2) and stores it in datatype.

#ifdef IBM_INDEXFIELD
  subroutine build_index_chunk_mpitype(ext_x, ext_y, ext_z, datatype_out)
    integer,intent(in) :: ext_x(2), ext_y(2), ext_z(2) !< chunk ranges for each dimension
    integer,intent(out) :: datatype_out !< MPI derived datatype to be built

    ! Declare variables.
    integer :: xrow, xyplane, xyzchunk ! temporary MPI datatypes
    integer :: cnt, blocklength, ierror, lengths(1) ! MPI variables
    integer(kind=MPI_ADDRESS_KIND) addr1, addr2, stride, base, offset, displs(1) ! addresses

    ! MPI datatype for columns along the x-axis
    cnt = 1 + ext_x(2) - ext_x(1)
    blocklength = 1 ! one integer for the rock gradient
    call MPI_Get_Address(interior_index(1, 1, 1), addr1, ierror)
    call MPI_Get_Address(interior_index(2, 1, 1), addr2, ierror)
    stride = addr2 - addr1
    call MPI_Type_create_hvector(cnt, blocklength, stride, MPI_REAL8, xrow, ierror)

    ! MPI datatype for slices in the xy-plane
    cnt = 1 + ext_y(2) - ext_y(1)
    blocklength = 1
    call MPI_Get_Address(interior_index(1, 1, 1), addr1, ierror)
    call MPI_Get_Address(interior_index(1, 2, 1), addr2, ierror)
    stride = addr2 - addr1
    call MPI_Type_create_hvector(cnt, blocklength, stride, xrow, xyplane, ierror)

    ! MPI datatype for whole chunk in xyz-volume
    cnt = 1 + ext_z(2) - ext_z(1)
    blocklength = 1
    call MPI_Get_Address(interior_index(1, 1, 1), addr1, ierror)
    call MPI_Get_Address(interior_index(1, 1, 2), addr2, ierror)
    stride = addr2 - addr1
    call MPI_Type_create_hvector(cnt, blocklength, stride, xyplane, xyzchunk, ierror)

    ! Position of the beginning of the chunk relative to the beginning of rock_gradient
    call MPI_Get_Address(interior_index, base, ierror)
    call MPI_Get_Address(interior_index(ext_x(1), ext_y(1), ext_z(1)), offset, ierror)

    ! Shift the memory relative to the base of vel_array
    cnt = 1
    lengths = (/1/)
    displs = (/offset - base/)
    call MPI_Type_create_hindexed(cnt, lengths, displs, xyzchunk, datatype_out, ierror)
    call MPI_Type_commit(datatype_out, ierror)
  end subroutine build_index_chunk_mpitype
#endif  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Exchange halo of physical velocity
  !>
  !> After computing the physical velocity on the lattice, the data has to be exchanged
  !> and stored in the halo of the neighboring processes so that IBM can access it.

  subroutine exchange_velocity_halo()
    ! Declare variables
    integer :: ierror ! MPI error
    integer, parameter :: cnt = 1 ! message count

    ! Send velocities along x-axis.
    call MPI_Sendrecv(vel_phys, cnt, datatype_send_vel_xpos, nnprocs(1,2), 1, &
          & vel_phys, cnt, datatype_recv_vel_xneg, nnprocs(1,1), 1, comm_cart, MPI_STATUS_IGNORE, ierror)
    call MPI_Sendrecv(vel_phys, cnt, datatype_send_vel_xneg, nnprocs(1,1), 2, &
          & vel_phys, cnt, datatype_recv_vel_xpos, nnprocs(1,2), 2, comm_cart, MPI_STATUS_IGNORE, ierror)

    ! Send velocities along y-axis.
    call MPI_Sendrecv(vel_phys, cnt, datatype_send_vel_ypos, nnprocs(2,2), 3, &
          & vel_phys, cnt, datatype_recv_vel_yneg, nnprocs(2,1), 3, comm_cart, MPI_STATUS_IGNORE, ierror)
    call MPI_Sendrecv(vel_phys, cnt, datatype_send_vel_yneg, nnprocs(2,1), 4, &
          & vel_phys, cnt, datatype_recv_vel_ypos, nnprocs(2,2), 4, comm_cart, MPI_STATUS_IGNORE, ierror)

    ! Send velocities along z-axis.
    call MPI_Sendrecv(vel_phys, cnt, datatype_send_vel_zpos, nnprocs(3,2), 5, &
          & vel_phys, cnt, datatype_recv_vel_zneg, nnprocs(3,1), 5, comm_cart, MPI_STATUS_IGNORE, ierror)
    call MPI_Sendrecv(vel_phys, cnt, datatype_send_vel_zneg, nnprocs(3,1), 6, &
          & vel_phys, cnt, datatype_recv_vel_zpos, nnprocs(3,2), 6, comm_cart, MPI_STATUS_IGNORE, ierror)
  end subroutine exchange_velocity_halo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Exchange halo of rock gradient
  !>
  !> In order to find the correct rock gradient, the data has to be exchanged and stored in the halo
  !> of the neighboring processes.

  subroutine exchange_rock_grad_halo()
    ! Declare variables
    integer :: ierror ! MPI error
    integer, parameter :: cnt = 1 ! message count

    ! Send velocities along x-axis.
    call MPI_Sendrecv(rock_grad, cnt, datatype_send_rock_grad_xpos, nnprocs(1,2), 1, &
          & rock_grad, cnt, datatype_recv_rock_grad_xneg, nnprocs(1,1), 1, comm_cart, MPI_STATUS_IGNORE, ierror)
    call MPI_Sendrecv(rock_grad, cnt, datatype_send_rock_grad_xneg, nnprocs(1,1), 2, &
          & rock_grad, cnt, datatype_recv_rock_grad_xpos, nnprocs(1,2), 2, comm_cart, MPI_STATUS_IGNORE, ierror)

    ! Send velocities along y-axis.
    call MPI_Sendrecv(rock_grad, cnt, datatype_send_rock_grad_ypos, nnprocs(2,2), 3, &
          & rock_grad, cnt, datatype_recv_rock_grad_yneg, nnprocs(2,1), 3, comm_cart, MPI_STATUS_IGNORE, ierror)
    call MPI_Sendrecv(rock_grad, cnt, datatype_send_rock_grad_yneg, nnprocs(2,1), 4, &
          & rock_grad, cnt, datatype_recv_rock_grad_ypos, nnprocs(2,2), 4, comm_cart, MPI_STATUS_IGNORE, ierror)

    ! Send velocities along z-axis.
    call MPI_Sendrecv(rock_grad, cnt, datatype_send_rock_grad_zpos, nnprocs(3,2), 5, &
          & rock_grad, cnt, datatype_recv_rock_grad_zneg, nnprocs(3,1), 5, comm_cart, MPI_STATUS_IGNORE, ierror)
    call MPI_Sendrecv(rock_grad, cnt, datatype_send_rock_grad_zneg, nnprocs(3,1), 6, &
          & rock_grad, cnt, datatype_recv_rock_grad_zpos, nnprocs(3,2), 6, comm_cart, MPI_STATUS_IGNORE, ierror)
  end subroutine exchange_rock_grad_halo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Exchange index field
  !>
  !> The index field for interior/exterior IBM regions is exchanged and stored in the halo of the neighboring processes.

#ifdef IBM_INDEXFIELD
  subroutine exchange_index_halo()
    ! Declare variables
    integer :: ierror ! MPI error
    integer, parameter :: cnt = 1 ! message count

    ! Send velocities along x-axis.
    call MPI_Sendrecv(interior_index, cnt, datatype_send_index_xpos, nnprocs(1,2), 1, &
          & interior_index, cnt, datatype_recv_index_xneg, nnprocs(1,1), 1, comm_cart, MPI_STATUS_IGNORE, ierror)
    call MPI_Sendrecv(interior_index, cnt, datatype_send_index_xneg, nnprocs(1,1), 2, &
          & interior_index, cnt, datatype_recv_index_xpos, nnprocs(1,2), 2, comm_cart, MPI_STATUS_IGNORE, ierror)

    ! Send velocities along y-axis.
    call MPI_Sendrecv(interior_index, cnt, datatype_send_index_ypos, nnprocs(2,2), 3, &
          & interior_index, cnt, datatype_recv_index_yneg, nnprocs(2,1), 3, comm_cart, MPI_STATUS_IGNORE, ierror)
    call MPI_Sendrecv(interior_index, cnt, datatype_send_index_yneg, nnprocs(2,1), 4, &
          & interior_index, cnt, datatype_recv_index_ypos, nnprocs(2,2), 4, comm_cart, MPI_STATUS_IGNORE, ierror)

    ! Send velocities along z-axis.
    call MPI_Sendrecv(interior_index, cnt, datatype_send_index_zpos, nnprocs(3,2), 5, &
          & interior_index, cnt, datatype_recv_index_zneg, nnprocs(3,1), 5, comm_cart, MPI_STATUS_IGNORE, ierror)
    call MPI_Sendrecv(interior_index, cnt, datatype_send_index_zneg, nnprocs(3,1), 6, &
          & interior_index, cnt, datatype_recv_index_zpos, nnprocs(3,2), 6, comm_cart, MPI_STATUS_IGNORE, ierror)
  end subroutine exchange_index_halo
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Distribute nodes to neighboring processes
  !>
  !> The particle nodes are distributed to the neighboring processes.
  !> Communication is executed in three steps:
  !> 1) along the x-axis
  !> 2) along the y-axis
  !> 3) along the z-axis
  !> This way, all nodes reach their target, also in diagonal directions.
  !> A process sends all nodes for which at least one of the following conditions is true:
  !> 1) The node is located outside of the physical domain of the master process.
  !>    In this case, another process is responsible for the node and the node is a halo node for the original process.
  !> 2) The node is inside the master's physical domain but close to the process boundary (i.e., in the halo of the neighbor).
  !>    In this case, the master is still responsible (the node is not a halo node),
  !>    but the node info is required by the neighbors as well.
  !> NOTE: It is assumed that nodes only propagate to the next neighbor.

  subroutine distribute_nodes()
    ! Declare variables.
    integer :: num_nodes_send_in_neg_x, num_nodes_send_in_pos_x ! number of nodes to send along the x-axis
    integer :: num_nodes_send_in_neg_y, num_nodes_send_in_pos_y ! number of nodes to send along the y-axis
    integer :: num_nodes_send_in_neg_z, num_nodes_send_in_pos_z ! number of nodes to send along the z-axis
    integer :: c_i ! particle index
    integer :: n_i ! node index
    integer :: status_neg_x(mpi_status_size), status_pos_x(mpi_status_size) ! MPI status for sending nodes along the x-axis
    integer :: status_neg_y(mpi_status_size), status_pos_y(mpi_status_size) ! MPI status for sending nodes along the y-axis
    integer :: status_neg_z(mpi_status_size), status_pos_z(mpi_status_size) ! MPI status for sending nodes along the z-axis
    integer :: ierror ! MPI error code
    integer :: num_nodes_recv_from_pos_x, num_nodes_recv_from_neg_x ! number of nodes received along the x-axis
    integer :: num_nodes_recv_from_pos_y, num_nodes_recv_from_neg_y ! number of nodes received along the y-axis
    integer :: num_nodes_recv_from_pos_z, num_nodes_recv_from_neg_z ! number of nodes received along the z-axis
    type(superstruct_node) :: temp_node ! temporary superstructure node

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("starting subroutine distribute_nodes...", .false.)
#endif

    ! Reset all counters.
    num_nodes_local = 0
    num_nodes_send_in_neg_x = 0
    num_nodes_send_in_pos_x = 0
    num_nodes_send_in_neg_y = 0
    num_nodes_send_in_pos_y = 0
    num_nodes_send_in_neg_z = 0
    num_nodes_send_in_pos_z = 0

    ! Run over all particles in the subdomain and about all of their nodes.
    ! For each of these nodes, the x-position is tested.
    ! If the node is near the boundary or beyond, it is added to the corresponding send buffer.
    ! In any case, it remains in the local node array (even if it is far in the neighboring domain).
    do c_i = 1, num_particles_loc
      do n_i = 1, particles(part_ind(c_i))%num_nodes
        ! Copy relevant node data to temporary node.
        temp_node%particle_gl = particles(part_ind(c_i))%particle_index_gl
        temp_node%particle_loc = c_i
        temp_node%node = n_i
        temp_node%pos = particles(part_ind(c_i))%node(n_i)%pos
        temp_node%vel = particles(part_ind(c_i))%node(n_i)%vel
        temp_node%force_int = 0
#ifdef IBM_DRAG
        temp_node%force_tot = particles(part_ind(c_i))%node(n_i)%force_s + particles(part_ind(c_i))%node(n_i)%force_b + &
          & particles(part_ind(c_i))%node(n_i)%force_v + particles(part_ind(c_i))%node(n_i)%force_at + particles(part_ind(c_i))%node(n_i)%force_c 
#else
        temp_node%force_tot = particles(part_ind(c_i))%node(n_i)%force_s + particles(part_ind(c_i))%node(n_i)%force_b + &
          & particles(part_ind(c_i))%node(n_i)%force_v + particles(part_ind(c_i))%node(n_i)%force_at
#endif

#ifdef IBM_FIXED
        temp_node%force_tot = temp_node%force_tot + particles(part_ind(c_i))%node(n_i)%force_anchor
#endif
!#ifdef IBM_INDEXFIELD
        temp_node%normal = particles(part_ind(c_i))%node(n_i)%normal ! only required in IBM_INDEXFIELD mode
        temp_node%curv_radius = particles(part_ind(c_i))%node(n_i)%curv_radius ! only required in IBM_INDEXFIELD mode
!#endif

        temp_node%send(1) = floor((temp_node%pos(1) - 0.5d0 - ccoords(1) * nx - HALO_WIDTH) / (nx - 2.d0 * HALO_WIDTH))
        temp_node%halo(1) = floor((temp_node%pos(1) - 0.5d0 - ccoords(1) * nx) / nx)

        ! Check if node has to be put into one of the send buffers.
        if(temp_node%send(1) .eq. -1) then
          num_nodes_send_in_neg_x = num_nodes_send_in_neg_x + 1
          send_in_neg_x(num_nodes_send_in_neg_x) = temp_node
        else if(temp_node%send(1) .eq. 1) then
          num_nodes_send_in_pos_x = num_nodes_send_in_pos_x + 1
          send_in_pos_x(num_nodes_send_in_pos_x) = temp_node
        end if

        ! In any case, the node resides in the local domain.
        num_nodes_local = num_nodes_local + 1
        nodes_local(num_nodes_local) = temp_node
      end do
    end do

    ! Communicate nodes along the x-axis.
    ! Two MPI calls (send-receive) are required to send the nodes in positive and negative direction.
    call MPI_Sendrecv(send_in_neg_x(1), num_nodes_send_in_neg_x, MPI_NODE_DATA, nnprocs(1,1), 1, &
      & recv_from_pos_x(1), MAX_NUM_NODES_SEND, MPI_NODE_DATA, nnprocs(1,2), 1, comm_cart, status_neg_x, ierror)
    call MPI_Sendrecv(send_in_pos_x(1), num_nodes_send_in_pos_x, MPI_NODE_DATA, nnprocs(1,2), 2, &
      & recv_from_neg_x(1), MAX_NUM_NODES_SEND, MPI_NODE_DATA, nnprocs(1,1), 2, comm_cart, status_pos_x, ierror)

    ! Identify the number of received nodes.
    call MPI_Get_count(status_neg_x, MPI_NODE_DATA, num_nodes_recv_from_pos_x, ierror)
    call MPI_Get_count(status_pos_x, MPI_NODE_DATA, num_nodes_recv_from_neg_x, ierror)

    ! Copy received nodes to local buffer.
    ! Periodicity has to be considered here if nodes are sent through a periodic boundary.
    do n_i = 1, num_nodes_recv_from_pos_x
      ! Check for periodicity effect.
      if(ccoords(1) .eq. (cdims(1) - 1)) then
        recv_from_pos_x(n_i)%pos(1) = recv_from_pos_x(n_i)%pos(1) + tnx
      end if

      ! Add node to local list.
      num_nodes_local = num_nodes_local + 1
      nodes_local(num_nodes_local) = recv_from_pos_x(n_i)
    end do

    do n_i = 1, num_nodes_recv_from_neg_x
      ! Check for periodicity effect.
      if(ccoords(1) .eq. 0) then
        recv_from_neg_x(n_i)%pos(1) = recv_from_neg_x(n_i)%pos(1) - tnx
      end if

      ! Add node to local list.
      num_nodes_local = num_nodes_local + 1
      nodes_local(num_nodes_local) = recv_from_neg_x(n_i)
    end do

    ! Run over all nodes in the local buffer.
    ! For each of these nodes, the y-position is tested.
    ! If the node is near the boundary or beyond, it is added to the corresponding send buffer.
    ! In any case, it remains in the local node array (even if it is far in the neighboring domain).
    do n_i = 1, num_nodes_local
      nodes_local(n_i)%send(2) = floor((nodes_local(n_i)%pos(2) - 0.5d0 - ccoords(2) * ny - HALO_WIDTH) / (ny - 2.d0 * HALO_WIDTH))
      nodes_local(n_i)%halo(2) = floor((nodes_local(n_i)%pos(2) - 0.5d0 - ccoords(2) * ny) / ny)

      ! Check if node has to be put into one of the send buffers.
      if(nodes_local(n_i)%send(2) .eq. -1) then
        num_nodes_send_in_neg_y = num_nodes_send_in_neg_y + 1
        send_in_neg_y(num_nodes_send_in_neg_y) = nodes_local(n_i)
      else if(nodes_local(n_i)%send(2) .eq. 1) then
        num_nodes_send_in_pos_y = num_nodes_send_in_pos_y + 1
        send_in_pos_y(num_nodes_send_in_pos_y) = nodes_local(n_i)
      end if
    end do

    ! Communicate nodes along the y-axis.
    ! Two MPI calls (send-receive) are required to send the nodes in positive and negative direction.
    call MPI_Sendrecv(send_in_neg_y(1), num_nodes_send_in_neg_y, MPI_NODE_DATA, nnprocs(2,1), 3, &
      & recv_from_pos_y(1), MAX_NUM_NODES_SEND, MPI_NODE_DATA, nnprocs(2,2), 3, comm_cart, status_neg_y, ierror)
    call MPI_Sendrecv(send_in_pos_y(1), num_nodes_send_in_pos_y, MPI_NODE_DATA, nnprocs(2,2), 4, &
      & recv_from_neg_y(1), MAX_NUM_NODES_SEND, MPI_NODE_DATA, nnprocs(2,1), 4, comm_cart, status_pos_y, ierror)

    ! Identify the number of received nodes.
    call MPI_Get_count(status_neg_y, MPI_NODE_DATA, num_nodes_recv_from_pos_y, ierror)
    call MPI_Get_count(status_pos_y, MPI_NODE_DATA, num_nodes_recv_from_neg_y, ierror)

    ! Copy received nodes to local buffer.
    ! Periodicity has to be considered here if nodes are sent through a periodic boundary.
    do n_i = 1, num_nodes_recv_from_pos_y
      ! Check for periodicity effect.
      if(ccoords(2) .eq. (cdims(2) - 1)) then
        recv_from_pos_y(n_i)%pos(2) = recv_from_pos_y(n_i)%pos(2) + tny
      end if

      ! Add node to local list.
      num_nodes_local = num_nodes_local + 1
      nodes_local(num_nodes_local) = recv_from_pos_y(n_i)
    end do

    do n_i = 1, num_nodes_recv_from_neg_y
      ! Check for periodicity effect.
      if(ccoords(2) .eq. 0) then
        recv_from_neg_y(n_i)%pos(2) = recv_from_neg_y(n_i)%pos(2) - tny
      end if

      ! Add node to local list.
      num_nodes_local = num_nodes_local + 1
      nodes_local(num_nodes_local) = recv_from_neg_y(n_i)
    end do

    ! Run over all nodes in the local buffer.
    ! For each of these nodes, the z-position is tested.
    ! If the node is near the boundary or beyond, it is added to the corresponding send buffer.
    ! In any case, it remains in the local node array (even if it is far in the neighboring domain).
    do n_i = 1, num_nodes_local
      nodes_local(n_i)%send(3) = floor((nodes_local(n_i)%pos(3) - 0.5d0 - ccoords(3) * nz - HALO_WIDTH) / (nz - 2.d0 * HALO_WIDTH))
      nodes_local(n_i)%halo(3) = floor((nodes_local(n_i)%pos(3) - 0.5d0 - ccoords(3) * nz) / nz)

      ! Check if node has to be put into one of the send buffers.
      if(nodes_local(n_i)%send(3) .eq. -1) then
        num_nodes_send_in_neg_z = num_nodes_send_in_neg_z + 1
        send_in_neg_z(num_nodes_send_in_neg_z) = nodes_local(n_i)
      else if(nodes_local(n_i)%send(3) .eq. 1) then
        num_nodes_send_in_pos_z = num_nodes_send_in_pos_z + 1
        send_in_pos_z(num_nodes_send_in_pos_z) = nodes_local(n_i)
      end if
    end do

    ! Communicate nodes along the z-axis.
    ! Two MPI calls (send-receive) are required to send the nodes in positive and negative direction.
    call MPI_Sendrecv(send_in_neg_z(1), num_nodes_send_in_neg_z, MPI_NODE_DATA, nnprocs(3,1), 5, &
      & recv_from_pos_z(1), MAX_NUM_NODES_SEND, MPI_NODE_DATA, nnprocs(3,2), 5, comm_cart, status_neg_z, ierror)
    call MPI_Sendrecv(send_in_pos_z(1), num_nodes_send_in_pos_z, MPI_NODE_DATA, nnprocs(3,2), 6, &
      & recv_from_neg_z(1), MAX_NUM_NODES_SEND, MPI_NODE_DATA, nnprocs(3,1), 6, comm_cart, status_pos_z, ierror)

    ! Identify the number of received nodes.
    call MPI_Get_count(status_neg_z, MPI_NODE_DATA, num_nodes_recv_from_pos_z, ierror)
    call MPI_Get_count(status_pos_z, MPI_NODE_DATA, num_nodes_recv_from_neg_z, ierror)

    ! Copy received nodes to local buffer.
    ! Periodicity has to be considered here if nodes are sent through a periodic boundary.
    do n_i = 1, num_nodes_recv_from_pos_z
      ! Check for periodicity effect.
      if(ccoords(3) .eq. (cdims(3) - 1)) then
        recv_from_pos_z(n_i)%pos(3) = recv_from_pos_z(n_i)%pos(3) + tnz
      end if

      ! Add node to local list.
      num_nodes_local = num_nodes_local + 1
      nodes_local(num_nodes_local) = recv_from_pos_z(n_i)
    end do

    do n_i = 1, num_nodes_recv_from_neg_z
      ! Check for periodicity effect.
      if(ccoords(3) .eq. 0) then
        recv_from_neg_z(n_i)%pos(3) = recv_from_neg_z(n_i)%pos(3) - tnz
      end if

      ! Add node to local list.
      num_nodes_local = num_nodes_local + 1
      nodes_local(num_nodes_local) = recv_from_neg_z(n_i)
    end do

    ! Reset all local interaction forces and check whether node is in the physical domain or not.
    do n_i = 1, num_nodes_local
      ! Reset local interaction force.
      nodes_local(n_i)%force_int = 0

      ! Check if node is in physical domain.
      if(nodes_local(n_i)%pos(1) - 0.5d0 - ccoords(1) * nx .ge. 0 .and. &
        & nodes_local(n_i)%pos(1) - 0.5d0  - ccoords(1) * nx .lt. nx .and. &
        & nodes_local(n_i)%pos(2) - 0.5d0 - ccoords(2) * ny .ge. 0 .and. &
        & nodes_local(n_i)%pos(2) - 0.5d0  - ccoords(2) * ny .lt. ny .and. &
        & nodes_local(n_i)%pos(3) - 0.5d0 - ccoords(3) * nz .ge. 0 .and. &
        & nodes_local(n_i)%pos(3) - 0.5d0  - ccoords(3) * nz .lt. nz) then
        nodes_local(n_i)%physical = 1
      else
        nodes_local(n_i)%physical = 0
      end if
    end do

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("finished subroutine distribute_nodes", .false.)
#endif
  end subroutine distribute_nodes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Collect nodes from neighboring processes
  !>
  !> The particle nodes are collected from the neighboring processes.
  !> Communication is executed in three steps:
  !> 1) along the x-axis
  !> 2) along the y-axis
  !> 3) along the z-axis
  !> This way, all nodes are collected, also from diagonal directions.
  !> A process sends all nodes for which at least one of the following conditions is true:
  !> NOTE: It is assumed that nodes only propagate to the next neighbor.

  subroutine collect_nodes()
    ! Declare variables.
    integer :: num_nodes_send_in_neg_x, num_nodes_send_in_pos_x ! number of nodes to send along the x-axis
    integer :: num_nodes_send_in_neg_y, num_nodes_send_in_pos_y ! number of nodes to send along the y-axis
    integer :: num_nodes_send_in_neg_z, num_nodes_send_in_pos_z ! number of nodes to send along the z-axis
    integer :: n_i, n_j ! node indices
    integer :: c_j ! particle index
    integer :: status_neg_x(mpi_status_size), status_pos_x(mpi_status_size) ! MPI status for sending nodes along the x-axis
    integer :: status_neg_y(mpi_status_size), status_pos_y(mpi_status_size) ! MPI status for sending nodes along the y-axis
    integer :: status_neg_z(mpi_status_size), status_pos_z(mpi_status_size) ! MPI status for sending nodes along the z-axis
    integer :: ierror ! MPI error code
    integer :: num_nodes_recv_from_pos_x, num_nodes_recv_from_neg_x ! number of nodes received along the x-axis
    integer :: num_nodes_recv_from_pos_y, num_nodes_recv_from_neg_y ! number of nodes received along the y-axis
    integer :: num_nodes_recv_from_pos_z, num_nodes_recv_from_neg_z ! number of nodes received along the z-axis

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("starting subroutine collect_nodes...", .false.)
#endif

    ! Reset all counters.
    num_nodes_send_in_neg_x = 0
    num_nodes_send_in_pos_x = 0
    num_nodes_send_in_neg_y = 0
    num_nodes_send_in_pos_y = 0
    num_nodes_send_in_neg_z = 0
    num_nodes_send_in_pos_z = 0

    ! Check nodes for send request along z-axis and copy to send buffers.
    ! Only bulk nodes are treated in the following. Halo nodes are not required for velocity interpolation.
    ! There are three options:
    ! 1) If the node obeys halo(3) == -1, it has to be sent in positive z-direction where it originally came from.
    ! 2) If the node obeys halo(3) ==  1, it has to be sent in negative z-direction where it originally came from.
    ! 3) If the node obeys halo(3) ==  0, there are two additional options:
    !    a) The particle has to be sent in x- or y-direction in a later stage of the algorithm.
    !    b) The particle is already in its original domain. The corresponding particle is updated with the interpolated velocity.
    do n_i = 1, num_nodes_local
      ! Only consider nodes in the physical regime.
      if(nodes_local(n_i)%physical .eq. 1) then
        ! Check z-membership of node
        if(nodes_local(n_i)%halo(3) .eq. -1) then
          num_nodes_send_in_pos_z = num_nodes_send_in_pos_z + 1
          send_in_pos_z(num_nodes_send_in_pos_z) = nodes_local(n_i)
        else if(nodes_local(n_i)%halo(3) .eq. 1) then
          num_nodes_send_in_neg_z = num_nodes_send_in_neg_z + 1
          send_in_neg_z(num_nodes_send_in_neg_z) = nodes_local(n_i)
        else
          ! Check y-membership of node
          if(nodes_local(n_i)%halo(2) .eq. -1) then
            num_nodes_send_in_pos_y = num_nodes_send_in_pos_y + 1
            send_in_pos_y(num_nodes_send_in_pos_y) = nodes_local(n_i)
          else if(nodes_local(n_i)%halo(2) .eq. 1) then
            num_nodes_send_in_neg_y = num_nodes_send_in_neg_y + 1
            send_in_neg_y(num_nodes_send_in_neg_y) = nodes_local(n_i)
          else
            ! Check x-membership of node
            if(nodes_local(n_i)%halo(1) .eq. -1) then
              num_nodes_send_in_pos_x = num_nodes_send_in_pos_x + 1
              send_in_pos_x(num_nodes_send_in_pos_x) = nodes_local(n_i)
            else if(nodes_local(n_i)%halo(1) .eq. 1) then
              num_nodes_send_in_neg_x = num_nodes_send_in_neg_x + 1
              send_in_neg_x(num_nodes_send_in_neg_x) = nodes_local(n_i)
            else
              c_j = nodes_local(n_i)%particle_loc
              n_j = nodes_local(n_i)%node
              particles(part_ind(c_j))%node(n_j)%vel = nodes_local(n_i)%vel
              particles(part_ind(c_j))%node(n_j)%force_int = nodes_local(n_i)%force_int
              particles(part_ind(c_j))%node(n_j)%force_tot = nodes_local(n_i)%force_tot
            end if
          end if
        end if
      end if
    end do

    ! Communicate nodes along the z-axis.
    ! Two MPI calls (send-receive) are required to send the nodes in positive and negative direction.
    call MPI_Sendrecv(send_in_neg_z(1), num_nodes_send_in_neg_z, MPI_NODE_DATA, nnprocs(3,1), 5, &
      & recv_from_pos_z(1), MAX_NUM_NODES_SEND, MPI_NODE_DATA, nnprocs(3,2), 5, comm_cart, status_neg_z, ierror)
    call MPI_Sendrecv(send_in_pos_z(1), num_nodes_send_in_pos_z, MPI_NODE_DATA, nnprocs(3,2), 6, &
      & recv_from_neg_z(1), MAX_NUM_NODES_SEND, MPI_NODE_DATA, nnprocs(3,1), 6, comm_cart, status_pos_z, ierror)

    ! Identify the number of received nodes.
    call MPI_Get_count(status_neg_z, MPI_NODE_DATA, num_nodes_recv_from_pos_z, ierror)
    call MPI_Get_count(status_pos_z, MPI_NODE_DATA, num_nodes_recv_from_neg_z, ierror)

    ! Check nodes for send request along y-axis and copy to send buffers.
    do n_i = 1, num_nodes_recv_from_pos_z
      if(recv_from_pos_z(n_i)%halo(2) .eq. -1) then
        num_nodes_send_in_pos_y = num_nodes_send_in_pos_y + 1
        send_in_pos_y(num_nodes_send_in_pos_y) = recv_from_pos_z(n_i)
      else if(recv_from_pos_z(n_i)%halo(2) .eq. 1) then
        num_nodes_send_in_neg_y = num_nodes_send_in_neg_y + 1
        send_in_neg_y(num_nodes_send_in_neg_y) = recv_from_pos_z(n_i)
      else
        if(recv_from_pos_z(n_i)%halo(1) .eq. -1) then
          num_nodes_send_in_pos_x = num_nodes_send_in_pos_x + 1
          send_in_pos_x(num_nodes_send_in_pos_x) = recv_from_pos_z(n_i)
        else if(recv_from_pos_z(n_i)%halo(1) .eq. 1) then
          num_nodes_send_in_neg_x = num_nodes_send_in_neg_x + 1
          send_in_neg_x(num_nodes_send_in_neg_x) = recv_from_pos_z(n_i)
        else
          c_j = recv_from_pos_z(n_i)%particle_loc
          n_j = recv_from_pos_z(n_i)%node
          particles(part_ind(c_j))%node(n_j)%vel = recv_from_pos_z(n_i)%vel
          particles(part_ind(c_j))%node(n_j)%force_int = recv_from_pos_z(n_i)%force_int
          particles(part_ind(c_j))%node(n_j)%force_tot = recv_from_pos_z(n_i)%force_tot
        end if
      end if
    end do

    do n_i = 1, num_nodes_recv_from_neg_z
      if(recv_from_neg_z(n_i)%halo(2) .eq. -1) then
        num_nodes_send_in_pos_y = num_nodes_send_in_pos_y + 1
        send_in_pos_y(num_nodes_send_in_pos_y) = recv_from_neg_z(n_i)
      else if(recv_from_neg_z(n_i)%halo(2) .eq. 1) then
        num_nodes_send_in_neg_y = num_nodes_send_in_neg_y + 1
        send_in_neg_y(num_nodes_send_in_neg_y) = recv_from_neg_z(n_i)
      else
        if(recv_from_neg_z(n_i)%halo(1) .eq. -1) then
          num_nodes_send_in_pos_x = num_nodes_send_in_pos_x + 1
          send_in_pos_x(num_nodes_send_in_pos_x) = recv_from_neg_z(n_i)
        else if(recv_from_neg_z(n_i)%halo(1) .eq. 1) then
          num_nodes_send_in_neg_x = num_nodes_send_in_neg_x + 1
          send_in_neg_x(num_nodes_send_in_neg_x) = recv_from_neg_z(n_i)
        else
          c_j = recv_from_neg_z(n_i)%particle_loc
          n_j = recv_from_neg_z(n_i)%node
          particles(part_ind(c_j))%node(n_j)%vel = recv_from_neg_z(n_i)%vel
          particles(part_ind(c_j))%node(n_j)%force_int = recv_from_neg_z(n_i)%force_int
          particles(part_ind(c_j))%node(n_j)%force_tot = recv_from_neg_z(n_i)%force_tot
        end if
      end if
    end do

    ! Communicate nodes along the y-axis.
    ! Two MPI calls (send-receive) are required to send the nodes in positive and negative direction.
    call MPI_Sendrecv(send_in_neg_y(1), num_nodes_send_in_neg_y, MPI_NODE_DATA, nnprocs(2,1), 3, &
      & recv_from_pos_y(1), MAX_NUM_NODES_SEND, MPI_NODE_DATA, nnprocs(2,2), 3, comm_cart, status_neg_y, ierror)
    call MPI_Sendrecv(send_in_pos_y(1), num_nodes_send_in_pos_y, MPI_NODE_DATA, nnprocs(2,2), 4, &
      & recv_from_neg_y(1), MAX_NUM_NODES_SEND, MPI_NODE_DATA, nnprocs(2,1), 4, comm_cart, status_pos_y, ierror)

    ! Identify the number of received nodes.
    call MPI_Get_count(status_neg_y, MPI_NODE_DATA, num_nodes_recv_from_pos_y, ierror)
    call MPI_Get_count(status_pos_y, MPI_NODE_DATA, num_nodes_recv_from_neg_y, ierror)

    ! Check nodes for send request along x-axis and copy to send buffers.
    do n_i = 1, num_nodes_recv_from_pos_y
      if(recv_from_pos_y(n_i)%halo(1) .eq. -1) then
        num_nodes_send_in_pos_x = num_nodes_send_in_pos_x + 1
        send_in_pos_x(num_nodes_send_in_pos_x) = recv_from_pos_y(n_i)
      else if(recv_from_pos_y(n_i)%halo(1) .eq. 1) then
        num_nodes_send_in_neg_x = num_nodes_send_in_neg_x + 1
        send_in_neg_x(num_nodes_send_in_neg_x) = recv_from_pos_y(n_i)
      else
        c_j = recv_from_pos_y(n_i)%particle_loc
        n_j = recv_from_pos_y(n_i)%node
        particles(part_ind(c_j))%node(n_j)%vel = recv_from_pos_y(n_i)%vel
        particles(part_ind(c_j))%node(n_j)%force_int = recv_from_pos_y(n_i)%force_int
        particles(part_ind(c_j))%node(n_j)%force_tot = recv_from_pos_y(n_i)%force_tot
      end if
    end do

    do n_i = 1, num_nodes_recv_from_neg_y
      if(recv_from_neg_y(n_i)%halo(1) .eq. -1) then
        num_nodes_send_in_pos_x = num_nodes_send_in_pos_x + 1
        send_in_pos_x(num_nodes_send_in_pos_x) = recv_from_neg_y(n_i)
      else if(recv_from_neg_y(n_i)%halo(1) .eq. 1) then
        num_nodes_send_in_neg_x = num_nodes_send_in_neg_x + 1
        send_in_neg_x(num_nodes_send_in_neg_x) = recv_from_neg_y(n_i)
      else
        c_j = recv_from_neg_y(n_i)%particle_loc
        n_j = recv_from_neg_y(n_i)%node
        particles(part_ind(c_j))%node(n_j)%vel = recv_from_neg_y(n_i)%vel
        particles(part_ind(c_j))%node(n_j)%force_int = recv_from_neg_y(n_i)%force_int
        particles(part_ind(c_j))%node(n_j)%force_tot = recv_from_neg_y(n_i)%force_tot
      end if
    end do

    ! Communicate nodes along the x-axis.
    ! Two MPI calls (send-receive) are required to send the nodes in positive and negative direction.
    call MPI_Sendrecv(send_in_neg_x(1), num_nodes_send_in_neg_x, MPI_NODE_DATA, nnprocs(1,1), 1, &
      & recv_from_pos_x(1), MAX_NUM_NODES_SEND, MPI_NODE_DATA, nnprocs(1,2), 1, comm_cart, status_neg_x, ierror)
    call MPI_Sendrecv(send_in_pos_x(1), num_nodes_send_in_pos_x, MPI_NODE_DATA, nnprocs(1,2), 2, &
      & recv_from_neg_x(1), MAX_NUM_NODES_SEND, MPI_NODE_DATA, nnprocs(1,1), 2, comm_cart, status_pos_x, ierror)

    ! Identify the number of received nodes.
    call MPI_Get_count(status_neg_x, MPI_NODE_DATA, num_nodes_recv_from_pos_x, ierror)
    call MPI_Get_count(status_pos_x, MPI_NODE_DATA, num_nodes_recv_from_neg_x, ierror)

    ! Update velocities
    do n_i = 1, num_nodes_recv_from_pos_x
      c_j = recv_from_pos_x(n_i)%particle_loc
      n_j = recv_from_pos_x(n_i)%node
      particles(part_ind(c_j))%node(n_j)%vel = recv_from_pos_x(n_i)%vel
      particles(part_ind(c_j))%node(n_j)%force_int = recv_from_pos_x(n_i)%force_int
      particles(part_ind(c_j))%node(n_j)%force_tot = recv_from_pos_x(n_i)%force_tot
    end do

    do n_i = 1, num_nodes_recv_from_neg_x
      c_j = recv_from_neg_x(n_i)%particle_loc
      n_j = recv_from_neg_x(n_i)%node
      particles(part_ind(c_j))%node(n_j)%vel = recv_from_neg_x(n_i)%vel
      particles(part_ind(c_j))%node(n_j)%force_int = recv_from_neg_x(n_i)%force_int
      particles(part_ind(c_j))%node(n_j)%force_tot = recv_from_neg_x(n_i)%force_tot
    end do

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("finished subroutine collect_nodes", .false.)
#endif
  end subroutine collect_nodes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Exchange particles between processes
  !>
  !> At the end of each time step, particles have to be sent to neighboring processes if they have crossed a domain boundary.
  !> The following actions will be performed:
  !> 1) It is checked whether the particle is still in its master's x-range. If not, it will be marked to be sent in x-direction.
  !> 2) After sending particles in x-direction (first in negative, then in positive direction),
  !>    it is checked whether particles have to be sent in y-direction.
  !> 3) After sending particles in y-direction (first in negative, then in positive direction),
  !>    it is checked whether particles have to be sent in z-direction.
  !> 4) Finally, particles are sent in z-direction (first in negative, then in positive direction)

  subroutine exchange_particles()
    ! Declare variables
    integer :: c_i ! particle index
    real(kind=rk), dimension(3) :: new_center ! new particle center position
    type(exchange_particle_info) :: send_pos_x, send_neg_x ! particle info to send along x-axis
    type(exchange_particle_info) :: recv_pos_x, recv_neg_x ! particle info to recv along x-axis
    type(exchange_particle_info) :: send_pos_y, send_neg_y ! particle info to send along y-axis
    type(exchange_particle_info) :: recv_pos_y, recv_neg_y ! particle info to recv along y-axis
    type(exchange_particle_info) :: send_pos_z, send_neg_z ! particle info to send along z-axis
    type(exchange_particle_info) :: recv_pos_z, recv_neg_z ! particle info to recv along z-axis
    integer :: status ! MPI status
    integer, dimension(NUM_PART_LOC_MAX) :: request_send_1, request_send_2 ! MPI request
    integer, dimension(NUM_PART_LOC_MAX) :: request_recv_1, request_recv_2 ! MPI request
    integer :: ierror ! MPI error
    character(len=200) :: message ! message string
#ifdef IBM_DEBUG
    integer :: num_particles_total ! total number of particles
#endif

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("starting subroutine exchange_particles...", .false.)
#endif

    !!! Send particles in negative x-direction

    ! Reset counter for number of particles to be sent in negative x-direction.
    send_neg_x%num_particles = 0

    ! Identify particles to be sent in negative x-direction and update the corresponding buffer.
    ! NOTE: The physical domain starts at 0.5.
    do c_i = 1, num_particles_loc
      if(particles(part_ind(c_i))%center(1) .lt. (ccoords(1) * nx + 0.5)) then
        send_neg_x%num_particles = send_neg_x%num_particles + 1
        send_neg_x%index_gl(send_neg_x%num_particles) = particles(part_ind(c_i))%particle_index_gl
        send_neg_x%index_loc(send_neg_x%num_particles) = c_i
        send_neg_x%mesh_type(send_neg_x%num_particles) = particles(part_ind(c_i))%mesh_type

        ! Correct for periodicity.
        if(ccoords(1) .eq. 0) then
          new_center = particles(part_ind(c_i))%center
          new_center(1) = new_center(1) + tnx
          call set_center_position(particles(part_ind(c_i)), new_center)
          particles(part_ind(c_i))%num_jumps(1) = particles(part_ind(c_i))%num_jumps(1) - 1
        end if
      end if
    end do

    ! Exchange mesh-type buffers with neighboring processes.
    call MPI_Sendrecv(send_neg_x, 3 * NUM_PART_LOC_MAX + 1, MPI_INTEGER, nnprocs(1,1), 1, &
            & recv_pos_x, 3 * NUM_PART_LOC_MAX + 1, MPI_INTEGER, nnprocs(1,2), 1, comm_cart, status, ierror)

    ! Send particles in negative x-direction.
    do c_i = 1, send_neg_x%num_particles
      call MPI_Isend(particles(part_ind(send_neg_x%index_loc(c_i)))%num_nodes, 1, MPI_LAGRANGE_PARTICLE_DATA, &
            & nnprocs(1,1), c_i, comm_cart, request_send_1(c_i), ierror)
      call MPI_Isend(particles(part_ind(send_neg_x%index_loc(c_i)))%node, &
            & particles(part_ind(send_neg_x%index_loc(c_i)))%num_nodes, MPI_LAGRANGE_NODE_DATA, nnprocs(1,1), &
            & c_i + NUM_PART_LOC_MAX, comm_cart, request_send_2(c_i), ierror)
    end do

    ! Receive particles from positive x-direction.
    do c_i = 1, recv_pos_x%num_particles
      call add_particle(meshes(recv_pos_x%mesh_type(c_i)))
      call MPI_Irecv(particles(part_ind(num_particles_loc))%num_nodes, 1, MPI_LAGRANGE_PARTICLE_DATA, nnprocs(1,2), &
            & c_i, comm_cart, request_recv_1(c_i), ierror)
      call MPI_Irecv(particles(part_ind(num_particles_loc))%node, particles(part_ind(num_particles_loc))%num_nodes, &
            & MPI_LAGRANGE_NODE_DATA, nnprocs(1,2), c_i + NUM_PART_LOC_MAX, comm_cart, request_recv_2(c_i), ierror)
    end do

    ! Wait for all communications to be completed.
    call MPI_Waitall(send_neg_x%num_particles, request_send_1, MPI_STATUSES_IGNORE, ierror)
    call MPI_Waitall(send_neg_x%num_particles, request_send_2, MPI_STATUSES_IGNORE, ierror)
    call MPI_Waitall(recv_pos_x%num_particles, request_recv_1, MPI_STATUSES_IGNORE, ierror)
    call MPI_Waitall(recv_pos_x%num_particles, request_recv_2, MPI_STATUSES_IGNORE, ierror)

    ! The send buffers are usable again, remove the particles which have been sent.
    do c_i = send_neg_x%num_particles, 1, -1
      call remove_particle(send_neg_x%index_loc(c_i))
    end do

    !!! Send particles in positive x-direction

    ! Reset counter for number of particles to be sent in positive x-direction.
    send_pos_x%num_particles = 0

    ! Identify particles to be sent in positive x-direction and update the corresponding buffer.
    ! NOTE: The physical domain starts at 0.5.
    do c_i = 1, num_particles_loc
      if(particles(part_ind(c_i))%center(1) .ge. ((ccoords(1) + 1) * nx + 0.5)) then
        send_pos_x%num_particles = send_pos_x%num_particles + 1
        send_pos_x%index_gl(send_pos_x%num_particles) = particles(part_ind(c_i))%particle_index_gl
        send_pos_x%index_loc(send_pos_x%num_particles) = c_i
        send_pos_x%mesh_type(send_pos_x%num_particles) = particles(part_ind(c_i))%mesh_type

        ! Correct for periodicity.
        if(ccoords(1) .eq. cdims(1) - 1) then
          new_center = particles(part_ind(c_i))%center
          new_center(1) = new_center(1) - tnx
          call set_center_position(particles(part_ind(c_i)), new_center)
          particles(part_ind(c_i))%num_jumps(1) = particles(part_ind(c_i))%num_jumps(1) + 1
        end if
      end if
    end do

    ! Exchange mesh-type buffers with neighboring processes.
    call MPI_Sendrecv(send_pos_x, 3 * NUM_PART_LOC_MAX + 1, MPI_INTEGER, nnprocs(1,2), 2, &
            & recv_neg_x, 3 * NUM_PART_LOC_MAX + 1, MPI_INTEGER, nnprocs(1,1), 2, comm_cart, status, ierror)

    ! Send particles in positive x-direction.
    do c_i = 1, send_pos_x%num_particles
      call MPI_Isend(particles(part_ind(send_pos_x%index_loc(c_i)))%num_nodes, 1, MPI_LAGRANGE_PARTICLE_DATA, &
            & nnprocs(1,2), c_i, comm_cart, request_send_1(c_i), ierror)
      call MPI_Isend(particles(part_ind(send_pos_x%index_loc(c_i)))%node, &
            & particles(part_ind(send_pos_x%index_loc(c_i)))%num_nodes, MPI_LAGRANGE_NODE_DATA, nnprocs(1,2), &
            & c_i + NUM_PART_LOC_MAX, comm_cart, request_send_2(c_i), ierror)
    end do

    ! Receive particles from positive x-direction.
    do c_i = 1, recv_neg_x%num_particles
      call add_particle(meshes(recv_neg_x%mesh_type(c_i)))
      call MPI_Irecv(particles(part_ind(num_particles_loc))%num_nodes, 1, MPI_LAGRANGE_PARTICLE_DATA, nnprocs(1,1), &
            & c_i, comm_cart, request_recv_1(c_i), ierror)
      call MPI_Irecv(particles(part_ind(num_particles_loc))%node, particles(part_ind(num_particles_loc))%num_nodes, &
            & MPI_LAGRANGE_NODE_DATA, nnprocs(1,1), c_i + NUM_PART_LOC_MAX, comm_cart, request_recv_2(c_i), ierror)
    end do

    ! Wait for all communications to be completed.
    call MPI_Waitall(send_pos_x%num_particles, request_send_1, MPI_STATUSES_IGNORE, ierror)
    call MPI_Waitall(send_pos_x%num_particles, request_send_2, MPI_STATUSES_IGNORE, ierror)
    call MPI_Waitall(recv_neg_x%num_particles, request_recv_1, MPI_STATUSES_IGNORE, ierror)
    call MPI_Waitall(recv_neg_x%num_particles, request_recv_2, MPI_STATUSES_IGNORE, ierror)

    ! The send buffers are usable again, remove the particles which have been sent.
    do c_i = send_pos_x%num_particles, 1, -1
      call remove_particle(send_pos_x%index_loc(c_i))
    end do

    !!! Send particles in negative y-direction

    ! Reset counter for number of particles to be sent in negative y-direction.
    send_neg_y%num_particles = 0

    ! Identify particles to be sent in negative y-direction and update the corresponding buffer.
    ! NOTE: The physical domain starts at 0.5.
    do c_i = 1, num_particles_loc
      if(particles(part_ind(c_i))%center(2) .lt. (ccoords(2) * ny + 0.5)) then
        send_neg_y%num_particles = send_neg_y%num_particles + 1
        send_neg_y%index_gl(send_neg_y%num_particles) = particles(part_ind(c_i))%particle_index_gl
        send_neg_y%index_loc(send_neg_y%num_particles) = c_i
        send_neg_y%mesh_type(send_neg_y%num_particles) = particles(part_ind(c_i))%mesh_type

        ! Correct for periodicity.
        if(ccoords(2) .eq. 0) then
          new_center = particles(part_ind(c_i))%center
          new_center(2) = new_center(2) + tny
          call set_center_position(particles(part_ind(c_i)), new_center)
          particles(part_ind(c_i))%num_jumps(2) = particles(part_ind(c_i))%num_jumps(2) - 1
        end if
      end if
    end do

    ! Exchange mesh-type buffers with neighboring processes.
    call MPI_Sendrecv(send_neg_y, 3 * NUM_PART_LOC_MAX + 1, MPI_INTEGER, nnprocs(2,1), 3, &
            & recv_pos_y, 3 * NUM_PART_LOC_MAX + 1, MPI_INTEGER, nnprocs(2,2), 3, comm_cart, status, ierror)

    ! Send particles in negative y-direction.
    do c_i = 1, send_neg_y%num_particles
      call MPI_Isend(particles(part_ind(send_neg_y%index_loc(c_i)))%num_nodes, 1, MPI_LAGRANGE_PARTICLE_DATA, &
            & nnprocs(2,1), c_i, comm_cart, request_send_1(c_i), ierror)
      call MPI_Isend(particles(part_ind(send_neg_y%index_loc(c_i)))%node, &
            & particles(part_ind(send_neg_y%index_loc(c_i)))%num_nodes, MPI_LAGRANGE_NODE_DATA, nnprocs(2,1), &
            & c_i + NUM_PART_LOC_MAX, comm_cart, request_send_2(c_i), ierror)
    end do

    ! Receive particles from positive y-direction.
    do c_i = 1, recv_pos_y%num_particles
      call add_particle(meshes(recv_pos_y%mesh_type(c_i)))
      call MPI_Irecv(particles(part_ind(num_particles_loc))%num_nodes, 1, MPI_LAGRANGE_PARTICLE_DATA, nnprocs(2,2), &
            & c_i, comm_cart, request_recv_1(c_i), ierror)
      call MPI_Irecv(particles(part_ind(num_particles_loc))%node, particles(part_ind(num_particles_loc))%num_nodes, &
            & MPI_LAGRANGE_NODE_DATA, nnprocs(2,2), c_i + NUM_PART_LOC_MAX, comm_cart, request_recv_2(c_i), ierror)
    end do

    ! Wait for all communications to be completed.
    call MPI_Waitall(send_neg_y%num_particles, request_send_1, MPI_STATUSES_IGNORE, ierror)
    call MPI_Waitall(send_neg_y%num_particles, request_send_2, MPI_STATUSES_IGNORE, ierror)
    call MPI_Waitall(recv_pos_y%num_particles, request_recv_1, MPI_STATUSES_IGNORE, ierror)
    call MPI_Waitall(recv_pos_y%num_particles, request_recv_2, MPI_STATUSES_IGNORE, ierror)

    ! The send buffers are usable again, remove the particles which have been sent.
    do c_i = send_neg_y%num_particles, 1, -1
      call remove_particle(send_neg_y%index_loc(c_i))
    end do

    !!! Send particles in positive y-direction

    ! Reset counter for number of particles to be sent in positive y-direction.
    send_pos_y%num_particles = 0

    ! Identify particles to be sent in negative y-direction and update the corresponding buffer.
    ! NOTE: The physical domain starts at 0.5.
    do c_i = 1, num_particles_loc
      if(particles(part_ind(c_i))%center(2) .ge. ((ccoords(2) + 1) * ny + 0.5)) then
        send_pos_y%num_particles = send_pos_y%num_particles + 1
        send_pos_y%index_gl(send_pos_y%num_particles) = particles(part_ind(c_i))%particle_index_gl
        send_pos_y%index_loc(send_pos_y%num_particles) = c_i
        send_pos_y%mesh_type(send_pos_y%num_particles) = particles(part_ind(c_i))%mesh_type

        ! Correct for periodicity.
        if(ccoords(2) .eq. cdims(2) - 1) then
          new_center = particles(part_ind(c_i))%center
          new_center(2) = new_center(2) - tny
          call set_center_position(particles(part_ind(c_i)), new_center)
          particles(part_ind(c_i))%num_jumps(2) = particles(part_ind(c_i))%num_jumps(2) + 1
        end if
      end if
    end do

    ! Exchange mesh-type buffers with neighboring processes.
    call MPI_Sendrecv(send_pos_y, 3 * NUM_PART_LOC_MAX + 1, MPI_INTEGER, nnprocs(2,2), 4, &
            & recv_neg_y, 3 * NUM_PART_LOC_MAX + 1, MPI_INTEGER, nnprocs(2,1), 4, comm_cart, status, ierror)

    ! Send particles in positive y-direction.
    do c_i = 1, send_pos_y%num_particles
      call MPI_Isend(particles(part_ind(send_pos_y%index_loc(c_i)))%num_nodes, 1, MPI_LAGRANGE_PARTICLE_DATA, &
            & nnprocs(2,2), c_i, comm_cart, request_send_1(c_i), ierror)
      call MPI_Isend(particles(part_ind(send_pos_y%index_loc(c_i)))%node, &
            & particles(part_ind(send_pos_y%index_loc(c_i)))%num_nodes, MPI_LAGRANGE_NODE_DATA, nnprocs(2,2), &
            & c_i + NUM_PART_LOC_MAX, comm_cart, request_send_2(c_i), ierror)
    end do

    ! Receive particles from positive y-direction.
    do c_i = 1, recv_neg_y%num_particles
      call add_particle(meshes(recv_neg_y%mesh_type(c_i)))
      call MPI_Irecv(particles(part_ind(num_particles_loc))%num_nodes, 1, MPI_LAGRANGE_PARTICLE_DATA, nnprocs(2,1), &
            & c_i, comm_cart, request_recv_1(c_i), ierror)
      call MPI_Irecv(particles(part_ind(num_particles_loc))%node, particles(part_ind(num_particles_loc))%num_nodes, &
            & MPI_LAGRANGE_NODE_DATA, nnprocs(2,1), c_i + NUM_PART_LOC_MAX, comm_cart, request_recv_2(c_i), ierror)
    end do

    ! Wait for all communications to be completed.
    call MPI_Waitall(send_pos_y%num_particles, request_send_1, MPI_STATUSES_IGNORE, ierror)
    call MPI_Waitall(send_pos_y%num_particles, request_send_2, MPI_STATUSES_IGNORE, ierror)
    call MPI_Waitall(recv_neg_y%num_particles, request_recv_1, MPI_STATUSES_IGNORE, ierror)
    call MPI_Waitall(recv_neg_y%num_particles, request_recv_2, MPI_STATUSES_IGNORE, ierror)

    ! The send buffers are usable again, remove the particles which have been sent.
    do c_i = send_pos_y%num_particles, 1, -1
      call remove_particle(send_pos_y%index_loc(c_i))
    end do

    !!! Send particles in negative z-direction

    ! Reset counter for number of particles to be sent in negative z-direction.
    send_neg_z%num_particles = 0

    ! Identify particles to be sent in negative z-direction and update the corresponding buffer.
    ! NOTE: The physical domain starts at 0.5.
    do c_i = 1, num_particles_loc
      if(particles(part_ind(c_i))%center(3) .lt. (ccoords(3) * nz + 0.5)) then
        send_neg_z%num_particles = send_neg_z%num_particles + 1
        send_neg_z%index_gl(send_neg_z%num_particles) = particles(part_ind(c_i))%particle_index_gl
        send_neg_z%index_loc(send_neg_z%num_particles) = c_i
        send_neg_z%mesh_type(send_neg_z%num_particles) = particles(part_ind(c_i))%mesh_type

        ! Correct for periodicity.
        if(ccoords(3) .eq. 0) then
          new_center = particles(part_ind(c_i))%center
          new_center(3) = new_center(3) + tnz
          call set_center_position(particles(part_ind(c_i)), new_center)
          particles(part_ind(c_i))%num_jumps(3) = particles(part_ind(c_i))%num_jumps(3) - 1
        end if
      end if
    end do

    ! Exchange mesh-type buffers with neighboring processes.
    call MPI_Sendrecv(send_neg_z, 3 * NUM_PART_LOC_MAX + 1, MPI_INTEGER, nnprocs(3,1), 5, &
            & recv_pos_z, 3 * NUM_PART_LOC_MAX + 1, MPI_INTEGER, nnprocs(3,2), 5, comm_cart, status, ierror)

    ! Send particles in negative z-direction.
    do c_i = 1, send_neg_z%num_particles
      call MPI_Isend(particles(part_ind(send_neg_z%index_loc(c_i)))%num_nodes, 1, MPI_LAGRANGE_PARTICLE_DATA, &
            & nnprocs(3,1), c_i, comm_cart, request_send_1(c_i), ierror)
      call MPI_Isend(particles(part_ind(send_neg_z%index_loc(c_i)))%node, &
            & particles(part_ind(send_neg_z%index_loc(c_i)))%num_nodes, MPI_LAGRANGE_NODE_DATA, nnprocs(3,1), &
            & c_i + NUM_PART_LOC_MAX, comm_cart, request_send_2(c_i), ierror)
    end do

    ! Receive particles from positive z-direction.
    do c_i = 1, recv_pos_z%num_particles
      call add_particle(meshes(recv_pos_z%mesh_type(c_i)))
      call MPI_Irecv(particles(part_ind(num_particles_loc))%num_nodes, 1, MPI_LAGRANGE_PARTICLE_DATA, nnprocs(3,2), &
            & c_i, comm_cart, request_recv_1(c_i), ierror)
      call MPI_Irecv(particles(part_ind(num_particles_loc))%node, particles(part_ind(num_particles_loc))%num_nodes, &
            & MPI_LAGRANGE_NODE_DATA, nnprocs(3,2), c_i + NUM_PART_LOC_MAX, comm_cart, request_recv_2(c_i), ierror)
    end do

    ! Wait for all communications to be completed.
    call MPI_Waitall(send_neg_z%num_particles, request_send_1, MPI_STATUSES_IGNORE, ierror)
    call MPI_Waitall(send_neg_z%num_particles, request_send_2, MPI_STATUSES_IGNORE, ierror)
    call MPI_Waitall(recv_pos_z%num_particles, request_recv_1, MPI_STATUSES_IGNORE, ierror)
    call MPI_Waitall(recv_pos_z%num_particles, request_recv_2, MPI_STATUSES_IGNORE, ierror)

    ! The send buffers are usable again, remove the particles which have been sent.
    do c_i = send_neg_z%num_particles, 1, -1
      call remove_particle(send_neg_z%index_loc(c_i))
    end do

    !!! Send particles in positive z-direction

    ! Reset counter for number of particles to be sent in positive z-direction.
    send_pos_z%num_particles = 0

    ! Identify particles to be sent in positive z-direction and update the corresponding buffer.
    ! NOTE: The physical domain starts at 0.5.
    do c_i = 1, num_particles_loc
      if(particles(part_ind(c_i))%center(3) .ge. ((ccoords(3) + 1) * nz + 0.5)) then
        send_pos_z%num_particles = send_pos_z%num_particles + 1
        send_pos_z%index_gl(send_pos_z%num_particles) = particles(part_ind(c_i))%particle_index_gl
        send_pos_z%index_loc(send_pos_z%num_particles) = c_i
        send_pos_z%mesh_type(send_pos_z%num_particles) = particles(part_ind(c_i))%mesh_type

        ! Correct for periodicity.
        if(ccoords(3) .eq. cdims(3) - 1) then
          new_center = particles(part_ind(c_i))%center
          new_center(3) = new_center(3) - tnz
          call set_center_position(particles(part_ind(c_i)), new_center)
          particles(part_ind(c_i))%num_jumps(3) = particles(part_ind(c_i))%num_jumps(3) + 1
        end if
      end if
    end do

    ! Exchange mesh-type buffers with neighboring processes.
    call MPI_Sendrecv(send_pos_z, 3 * NUM_PART_LOC_MAX + 1, MPI_INTEGER, nnprocs(3,2), 6, &
            & recv_neg_z, 3 * NUM_PART_LOC_MAX + 1, MPI_INTEGER, nnprocs(3,1), 6, comm_cart, status, ierror)

    ! Send particles in positive z-direction.
    do c_i = 1, send_pos_z%num_particles
      call MPI_Isend(particles(part_ind(send_pos_z%index_loc(c_i)))%num_nodes, 1, MPI_LAGRANGE_PARTICLE_DATA, &
            & nnprocs(3,2), c_i, comm_cart, request_send_1(c_i), ierror)
      call MPI_Isend(particles(part_ind(send_pos_z%index_loc(c_i)))%node, &
            & particles(part_ind(send_pos_z%index_loc(c_i)))%num_nodes, MPI_LAGRANGE_NODE_DATA, nnprocs(3,2), &
            & c_i + NUM_PART_LOC_MAX, comm_cart, request_send_2(c_i), ierror)
    end do

    ! Receive particles from positive z-direction.
    do c_i = 1, recv_neg_z%num_particles
      call add_particle(meshes(recv_neg_z%mesh_type(c_i)))
      call MPI_Irecv(particles(part_ind(num_particles_loc))%num_nodes, 1, MPI_LAGRANGE_PARTICLE_DATA, nnprocs(3,1), &
            & c_i, comm_cart, request_recv_1(c_i), ierror)
      call MPI_Irecv(particles(part_ind(num_particles_loc))%node, particles(part_ind(num_particles_loc))%num_nodes, &
            & MPI_LAGRANGE_NODE_DATA, nnprocs(3,1), c_i + NUM_PART_LOC_MAX, comm_cart, request_recv_2(c_i), ierror)
    end do

    ! Wait for all communications to be completed.
    call MPI_Waitall(send_pos_z%num_particles, request_send_1, MPI_STATUSES_IGNORE, ierror)
    call MPI_Waitall(send_pos_z%num_particles, request_send_2, MPI_STATUSES_IGNORE, ierror)
    call MPI_Waitall(recv_neg_z%num_particles, request_recv_1, MPI_STATUSES_IGNORE, ierror)
    call MPI_Waitall(recv_neg_z%num_particles, request_recv_2, MPI_STATUSES_IGNORE, ierror)

    ! The send buffers are usable again, remove the particles which have been sent.
    do c_i = send_pos_z%num_particles, 1, -1
      call remove_particle(send_pos_z%index_loc(c_i))
    end do

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Reduce(num_particles_loc, num_particles_total, 1, MPI_INTEGER, MPI_SUM, 0, comm_cart, ierror) 
    write(message, "('counting ', i0, ' particles in total after particle exchange')") num_particles_total
    call log_msg(trim(message), .false.)
#endif

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("finished subroutine exchange_particles", .false.)
#endif
  end subroutine exchange_particles

! endif IBM_PART
#endif

end module lsuperstruct_parallel_module
