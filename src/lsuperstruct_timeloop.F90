!> Lagrangian superstructure timeloop module
!>
!> This module contains all subroutines responsible for maintaining the Lagrangian superstructure timeloop.

#include "lbe.h"

module lsuperstruct_timeloop_module

#ifdef IBM_PART

  ! Include external modules.
  use hoshen_kopelman_module
  use lbe_force_interface_module, only : add_force_to_all
  use lbe_helper_module, only : norm, cross_product
  use lbe_log_module
  use lbe_parallel_module, only : comm_cart, ccoords, tny, tnz
  use lbe_parms_module, only : nx, ny, nz, nt, tau_r
  use lbe_types_module, only : lbe_site
  use lsuperstruct_data_module
  use lsuperstruct_helper_module, only : distance_vector, set_interior_index
#ifdef IBM_INDEXFIELD
  use lextobj_module, only : update_node_normals
#endif
#if defined IBM_INDEXFIELD && defined VARTAU
  use lsuperstruct_interface_module, only: set_local_tau_r
#endif
#ifdef IBM_BINARYIBM
!   use lsuperstruct_interface_module, only: recolor_fluid, get_neighboring_color_ratio
#endif
  use lextobj_module, only : particles, part_ind, forces_compute, update_node_positions, update_face_areas, &
      & update_node_areas, update_particle_center, update_momenta
      

  implicit none
  include 'mpif.h'

  private
  public :: compute_particle_forces_int, compute_particle_forces_ext, combine_forces, update_particles_pre, update_particles_post, &
            & update_node_interaction_LUT, compute_wall_gradient_forces, compute_particle_particle_forces, compute_particle_wall_forces
#ifdef IBM_INDEXFIELD
  public :: rebuild_interior_index
#endif

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute internal particle forces
  !>
  !> The superstructure instructs the Lagrangian extended object module to compute the internal forces.
  !> These are forces which are defined for individual particles (e.g., strain, volume).
  !> External forces (e.g., particle interactions with other particles or walls) are not contained here,
  !> they are rather treated in \c compute_particle_forces_ext.
  !> TODO: Revise location of calculation of swimmer forces.

  subroutine compute_particle_forces_int(time)
    integer, intent(in) :: time !< simulation time
  
    ! Declare variables.
    integer :: c_i ! particle index
    integer :: ierror ! MPI errir code

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("starting subroutine compute_particle_forces_int...", .false.)
#endif

    ! Call internal force routine in particle module.
    do c_i = 1, num_particles_loc
      call forces_compute(particles(part_ind(c_i)), surften, contactangle, width, friction, time)
    end do

    ! Compute swimmer forces.
#ifdef IBM_SWIMMER
    call compute_swimmer_spring_forces(time) ! forces on swimmers due to spring interactions
#endif

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("finished subroutine compute_particle_forces_int", .false.)
#endif

  end subroutine compute_particle_forces_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute external particle forces
  !>
  !> The external particle forces (e.g., interactions) are computed.
  !> Internal forces (e.g., strain or bending) are not contained here,
  !> they are rather treated in \c compute_particle_forces_int.

  subroutine compute_particle_forces_ext(N)
    type(lbe_site), dimension(0:, 0:, 0:), intent(in) :: N !< lattice

    ! Call external forcing subroutines.
    call compute_particle_particle_forces() ! particle-particle interaction forces
    call compute_particle_wall_forces(N) ! wall interaction forces
  end subroutine compute_particle_forces_ext

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Combine internal and external forces
  !>
  !> The internal and external forces are combined in order to find the total force acting on each node.
  !> These forces are then required by the immersed boundary method.
  !> NOTE: Make sure that all force contributions are considered in the summation.

  subroutine combine_forces()
    ! Declare variables.
    integer :: n_i ! node index
    integer :: ierror ! MPI error

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("starting subroutine combine_forces...", .false.)
#endif

    ! Combine forces.
    do n_i = 1, num_nodes_local
      nodes_local(n_i)%force_tot = nodes_local(n_i)%force_tot + nodes_local(n_i)%force_int
    end do

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("finished subroutine combine_forces", .false.)
#endif

  end subroutine combine_forces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update particle positions (before particle exchange)
  !>
  !> The node positions are updated.

  subroutine update_particles_pre()
    ! Declare variables.
    integer :: c_i ! particle index

    ! Call update routines.
    do c_i = 1, num_particles_loc
      call update_node_positions(particles(part_ind(c_i)), IBM_massdensity, friction)
    end do
  end subroutine update_particles_pre

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update particle positions (after particle exchange)
  !>
  !> The face and node areas, particle center positions, and momenta are updated.

  subroutine update_particles_post()
    ! Declare variables.
    integer :: c_i ! particle index
    integer :: ierror ! MPI errir code

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("starting subroutine update_particles_post...", .false.)
#endif

    ! Reset number of nodes and faces in subdomain.
    num_nodes_loc_vtk = 0
    num_faces_loc_vtk = 0

    ! Call update routines.
    do c_i = 1, num_particles_loc
      num_nodes_loc_vtk = num_nodes_loc_vtk + particles(part_ind(c_i))%num_nodes
      num_faces_loc_vtk = num_faces_loc_vtk + particles(part_ind(c_i))%num_faces
      call update_face_areas(particles(part_ind(c_i)))
      call update_node_areas(particles(part_ind(c_i)))
#ifdef IBM_INDEXFIELD
      call update_node_normals(particles(part_ind(c_i)))
#endif
      call update_particle_center(particles(part_ind(c_i)))
      call update_momenta(particles(part_ind(c_i)))
    end do

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("finished subroutine update_particles_post", .false.)
#endif
  end subroutine update_particles_post

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Particle-particle interaction forces
  !>
  !> The node-node interaction forces are computed.
  !> For this, the look-up table for the node locations is required.
  !> The minimum distance between nodes is tracked in the process, it can be used for evaluation of the data.
  !> The algorithm consists of the following steps:
  !> 1) Run over all nodes.
  !> 2) For each node, scan its neighborhood using the look-up table.
  !> 3) Interaction forces will only be computed if the neighboring node
  !>    - is member of another particle than the original node,
  !>    - has not been considered before,
  !>    - is within interaction range.
  !> 4) The minimum distance between nodes is computed (start value is INTERACTION_RANGE).

  subroutine compute_particle_particle_forces()
    ! Declare variables.
    integer :: range_int ! numerical interaction range (number of required nodes)
    integer :: n_i, n_j ! node indices
    integer :: pos_x, pos_y, pos_z ! position of node on lattice
    integer :: c_i, c_j ! particle indices
    integer :: X, Y, Z ! positions of neighboring subboxes
    real(kind=rk), dimension(3) :: dist_vec ! distance vector between two node positions
    real(kind=rk) :: dist ! distance between two node positions
    real(kind=rk), dimension(3) :: force ! force between nodes
    integer :: ierror ! MPI error code
    character(len=200) :: message ! message string

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("starting subroutine compute_particle_particle_forces...", .false.)
#endif

    ! Check evaluation condition.
    ! If less than two particles are located in the flow, the interaction force is not computed.
    ! The minimum distance is set to zero.
    if((num_particles_gl .lt. 2) .or. (int_strength_part_part .eq. 0.0d0)) then
      ! Debugging routine
#ifdef IBM_DEBUG
      call MPI_Barrier(comm_cart, ierror)
      call log_msg("finished subroutine compute_particle_particle_forces", .false.)
#endif
      return
    end if

    ! Set interaction range.
    ! range_int is the required numerical interaction range.
    ! Initially set minimum distance to interaction range.
    range_int = ceiling(INT_RANGE_PART_PART)
    dist_min_loc = INT_RANGE_PART_PART

    ! Compute interaction forces.
    ! Run over all nodes in the subdomain.
    do n_i = 1, num_nodes_local
      ! Find position of node on lattice.
      ! For this, the global coordinate has to be converted to a local coordinate.
      ! NOTE: The same convention as in update_node_interaction_LUT() is used.
      pos_x = floor(nodes_local(n_i)%pos(1) + 0.5 - ccoords(1) * nx)
      pos_y = floor(nodes_local(n_i)%pos(2) + 0.5 - ccoords(2) * ny)
      pos_z = floor(nodes_local(n_i)%pos(3) + 0.5 - ccoords(3) * nz)

      ! Check node position.
      ! Only proceed if the node is located within the LUT region.
      if((pos_x .ge. 1 - HALO_WIDTH) .and. (pos_x .le. nx + HALO_WIDTH) .and. &
        & (pos_y .ge. 1 - HALO_WIDTH) .and. (pos_y .le. ny + HALO_WIDTH) .and. &
        & (pos_z .ge. 1 - HALO_WIDTH) .and. (pos_z .le. nz + HALO_WIDTH)) then
        ! Obtain the global particle index for this node.
        c_i = nodes_local(n_i)%particle_gl

        ! Loop over all x-positions.
        do X = pos_x - range_int, pos_x + range_int
          ! Check if neighboring box is still in the LUT region.
          if((X .lt. 1 - HALO_WIDTH) .or. (X .gt. nx + HALO_WIDTH)) cycle

          ! Loop over all y-positions.
          do Y = pos_y - range_int, pos_y + range_int
            ! Check if neighboring box is still in the LUT region.
            if((Y .lt. 1 - HALO_WIDTH) .or. (Y .gt. ny + HALO_WIDTH)) cycle

            ! Loop over all z-positions.
            do Z = pos_z - range_int, pos_z + range_int
              ! Check if neighboring box is still in the LUT region.
              if((Z .lt. 1 - HALO_WIDTH) .or. (Z .gt. nz + HALO_WIDTH)) cycle

              ! Loop over all nodes located in the subbox.
              do n_j = 1, LUT_nodes_in_range(X, Y, Z, 0)
                ! Check node ordering.
                ! An interaction force is only computed if the global particle index c_i is larger than c_j.
                ! This makes sure that interaction is considered only once per pair and not within particles.
                c_j = nodes_local(LUT_nodes_in_range(X, Y, Z, n_j))%particle_gl
                if(c_j .le. c_i) cycle

                ! Compute distance between the nodes.
                ! Neglect all node pairs with distances larger than the interaction range.
                ! Update the minimum distance between nodes.
                dist_vec = distance_vector(nodes_local(LUT_nodes_in_range(X, Y, Z, n_j))%pos, nodes_local(n_i)%pos)
                dist = norm(dist_vec)
                if(dist .ge. INT_RANGE_PART_PART) cycle
                if(dist .lt. dist_min_loc) dist_min_loc = dist

                ! Compute the interaction force and distribute it to the pair of nodes.
                force = (dist_vec / dist) * int_strength_part_part * (dist**(-2) - INT_RANGE_PART_PART**(-2))
                nodes_local(n_i)%force_int = nodes_local(n_i)%force_int + force
                nodes_local(LUT_nodes_in_range(X, Y, Z, n_j))%force_int = &
                  & nodes_local(LUT_nodes_in_range(X, Y, Z, n_j))%force_int - force
              end do
            end do
          end do
        end do
      end if
    end do

    ! Distribute the minimum distance to all processes
    call MPI_Allreduce(dist_min_loc, dist_min_gl, 1, LBE_REAL, MPI_MIN, comm_cart, ierror)

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("finished subroutine compute_particle_particle_forces", .false.)
#endif

  end subroutine compute_particle_particle_forces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Particle-wall interaction forces
  !>
  !> The node-rock node interaction forces are computed.
  !> The minimum distance between nodes and rock nodes is tracked in the process, it can be used for evaluation of the data.
  !> The algorithm consists of the following steps:
  !> 1) Run over all nodes.
  !> 2) For each node, scan its neighborhood for rock nodes.
  !> 3) Interaction forces will only be computed if the neighboring lattice node
  !>    - is a rock node (i.e., not a fluid node),
  !>    - is within interaction range.
  !> 4) The minimum distance between particle and rock nodes is computed.

  subroutine compute_particle_wall_forces(N)
    type(lbe_site), dimension(0:, 0:, 0:), intent(in) :: N !< lattice

    ! Declare variables.
    integer :: range_int ! numerical interaction range (number of required nodes)
    integer :: n_i ! node index
    integer :: pos_x, pos_y, pos_z ! position of node on lattice
    integer :: X, Y, Z ! positions of neighboring subboxes
    real(kind=rk), dimension(3) :: dist_vec ! distance vector between two node positions
    real(kind=rk) :: dist ! distance between two node positions
    real(kind=rk), dimension(3) :: latt_node_pos ! position of lattice node
    real(kind=rk), dimension(3) :: force ! force between nodes
    integer :: ierror ! MPI error code

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("starting subroutine compute_particle_wall_forces...", .false.)
#endif

    ! Set interaction range.
    ! range_int is the required numerical interaction range.
    ! Initially set minimum distance to interaction range.
    range_int = ceiling(INT_RANGE_WALL_PART - 0.5)
    dist_min_loc = INT_RANGE_WALL_PART

    ! Check execution condition.
    if(INT_RANGE_WALL_PART .eq. 0.0d0 .or. int_strength_part_wall .eq. 0.0d0) return

    ! Compute interaction forces.
    ! Run over all nodes in the subdomain.
    do n_i = 1, num_nodes_local
      ! Find position of node on lattice.
      ! For this, the global coordinate has to be converted to a local coordinate.
      ! NOTE: The same convention as in update_node_interaction_LUT() is used.
      pos_x = floor(nodes_local(n_i)%pos(1) + 0.5 - ccoords(1) * nx)
      pos_y = floor(nodes_local(n_i)%pos(2) + 0.5 - ccoords(2) * ny)
      pos_z = floor(nodes_local(n_i)%pos(3) + 0.5 - ccoords(3) * nz)

      ! Check node position.
      ! Only proceed if the node is located within the LUT region.
      if((pos_x .ge. 1 - HALO_WIDTH) .and. (pos_x .le. nx + HALO_WIDTH) .and. &
        & (pos_y .ge. 1 - HALO_WIDTH) .and. (pos_y .le. ny + HALO_WIDTH) .and. &
        & (pos_z .ge. 1 - HALO_WIDTH) .and. (pos_z .le. nz + HALO_WIDTH)) then

        if(N(pos_x, pos_y, pos_z)%rock_state .ne. 0.0d0) cycle

        ! Loop over all x-positions.
        do X = pos_x - range_int, pos_x + range_int
          ! Check if neighboring box is still in the LUT region.
          if((X .lt. 1 - HALO_WIDTH) .or. (X .gt. nx + HALO_WIDTH)) cycle

          ! Loop over all y-positions.
          do Y = pos_y - range_int, pos_y + range_int
            ! Check if neighboring box is still in the LUT region.
            if((Y .lt. 1 - HALO_WIDTH) .or. (Y .gt. ny + HALO_WIDTH)) cycle

            ! Loop over all z-positions.
            do Z = pos_z - range_int, pos_z + range_int
              ! Check if neighboring box is still in the LUT region.
              if((Z .lt. 1 - HALO_WIDTH) .or. (Z .gt. nz + HALO_WIDTH)) cycle

              ! Check whether target node is rock node.
              ! If it is a fluid node, cycle to next iteration.
              if(N(X, Y, Z)%rock_state .eq. 0.0d0) cycle

              ! Compute distance between node and lattice node.
              ! Neglect all node pairs with distances larger than the interaction range.
              latt_node_pos = (/X + ccoords(1) * nx, Y + ccoords(2) * ny, Z + ccoords(3) * nz/)
              dist_vec = distance_vector(latt_node_pos, nodes_local(n_i)%pos)
              dist = norm(dist_vec)
              if(dist .ge. INT_RANGE_WALL_PART) cycle

              ! Compute the interaction force.
              ! Use cutoff for the force so that it does not diverge for small distances.
              force = (dist_vec / dist) * int_strength_part_wall * (dist**(-2) - INT_RANGE_WALL_PART**(-2))
              if(norm(force) > 0.1) force = force * 0.1 / norm(force)
              nodes_local(n_i)%force_int = nodes_local(n_i)%force_int + force
              call add_force_to_all(-force(1), -force(2), -force(3), X, Y, Z)
            end do
          end do
        end do
      end if
    end do

    ! Distribute the minimum distance to all processes
!    call MPI_Allreduce(dist_min_loc, dist_min_gl, 1, LBE_REAL, MPI_MIN, comm_cart, ierror)

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("finished subroutine compute_particle_wall_forces", .false.)
#endif
  end subroutine compute_particle_wall_forces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Wall gradient forces
  !>
  !> During their growth, particles may still be located inside a wall.
  !> In order to repell them towards the fluid region, a wall gradient force is acting.
  !> This force points in the direction of the nearest fluid region.

  subroutine compute_wall_gradient_forces()
    ! Declare variables.
    integer :: n_i ! node index
    integer :: pos_x, pos_y, pos_z ! position of node on lattice
    integer :: X, Y, Z ! positions of neighboring subboxes
    integer :: ierror ! MPI error code
    real(kind=rk), dimension(3) :: force ! force between nodes

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("starting subroutine compute_wall_gradient_forces...", .false.)
#endif

    ! Compute gradient forces.
    ! Run over all nodes in the subdomain.
    do n_i = 1, num_nodes_local
      ! Find position of node on lattice.
      ! For this, the global coordinate has to be converted to a local coordinate.
      ! NOTE: The same convention as in update_node_interaction_LUT() is used.
      pos_x = floor(nodes_local(n_i)%pos(1) + 0.5 - ccoords(1) * nx)
      pos_y = floor(nodes_local(n_i)%pos(2) + 0.5 - ccoords(2) * ny)
      pos_z = floor(nodes_local(n_i)%pos(3) + 0.5 - ccoords(3) * nz)

      ! Check node position.
      ! Only proceed if the node is located within the LUT region.
      if((pos_x .ge. 1) .and. (pos_x .le. nx) .and. &
        & (pos_y .ge. 1) .and. (pos_y .le. ny) .and. &
        & (pos_z .ge. 1) .and. (pos_z .le. nz)) then
        ! Loop over all x-positions.
        force(:) = -gradient(:, pos_x, pos_y, pos_z) * int_strength_gradient
        nodes_local(n_i)%force_int = nodes_local(n_i)%force_int + force
      end if
    end do

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("finished subroutine compute_wall_gradient_forces", .false.)
#endif
  end subroutine compute_wall_gradient_forces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Swimmer spring forces
  !>
  !> The swimmer forces are computed.
  !> This is a very quick implementation which will not scale for a large number of swimmers:
  !> * Swimmer positions are collected by rank 0.
  !> * Rank 0 computes forces and broadcasts them.
  !> * Responsible ranks add forces to interaction forces.
  !> The following forces are considered:
  !> * External driving forces (specified by user)
  !> * Elastic spring forces according to specified spring constants
  !> * Bending rigidity to enforce linearity of springer
  !> * If desired, anchor force to restrict swimmer to original axis

#ifdef IBM_SWIMMER
  subroutine compute_swimmer_spring_forces(time)
    integer, intent(in) :: time !< simulation time
  
    ! Declare variables.
    integer :: n_i ! node index
    integer :: c_i ! cell index
    real(kind=rk), dimension(3) :: pos1, pos2, pos3 ! position of swimmer beads
    real(kind=rk), dimension(3) :: pos1_red, pos2_red, pos3_red ! position of swimmer beads (rank 0)
    real(kind=rk), dimension(3) :: pos_ref ! reference position for anchor force
    real(kind=rk), dimension(3) :: d_anchor1, d_anchor2, d_anchor3 ! deviations from reference position
    real(kind=rk), dimension(3) :: dvec12, dvec23, dvec13 ! distance vectors between bead pairs
    real(kind=rk) :: d12, d23, d13 ! distances between bead pairs
    real(kind=rk), dimension(3) :: cross ! cross product of distance vectors (indicating bending magnitude)
    real(kind=rk), dimension(3) :: force1_d, force2_d, force3_d ! driving forces
    real(kind=rk), dimension(3) :: force1_s, force2_s, force3_s ! spring forces
    real(kind=rk), dimension(3) :: force1_b, force2_b, force3_b ! bending forces
    real(kind=rk), dimension(3) :: force1_a, force2_a, force3_a ! anchor forces
    real(kind=rk), dimension(3) :: force1, force2, force3 ! total forces
    integer :: ierror ! MPI error code
    character(len=200) :: message ! message string

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("starting subroutine compute_swimmer_spring_forces...")
#endif

    ! Set particle positions.
    pos1 = 0.
    pos2 = 0.
    pos3 = 0.

    do c_i = 1, num_particles_loc
      if(particles(part_ind(c_i))%particle_index_gl .eq. 1) then
        pos1 = particles(part_ind(c_i))%center
      else if(particles(part_ind(c_i))%particle_index_gl .eq. 2) then
        pos2 = particles(part_ind(c_i))%center
      else if(particles(part_ind(c_i))%particle_index_gl .eq. 3) then
        pos3 = particles(part_ind(c_i))%center
      end if
    end do
    
    ! Reduce position data so that every process knows all swimmer locations.
    call MPI_Allreduce(pos1, pos1_red, 3, MPI_REAL8, MPI_SUM, comm_cart, ierror)
    call MPI_Allreduce(pos2, pos2_red, 3, MPI_REAL8, MPI_SUM, comm_cart, ierror)
    call MPI_Allreduce(pos3, pos3_red, 3, MPI_REAL8, MPI_SUM, comm_cart, ierror)

    ! Compute geometry.
    dvec12 = distance_vector(pos1_red, pos2_red)
    dvec23 = distance_vector(pos2_red, pos3_red)
    dvec13 = distance_vector(pos1_red, pos3_red)
    d12 = norm(dvec12)
    d23 = norm(dvec23)
    d13 = norm(dvec13)
    cross = cross_product(dvec12, dvec23) / (d12 * d23)
    
    ! Compute driving forces.
    ! It is assumed that all these forces act along one axis.
    ! This axis is defined by the positions of beads 1 and 3.
    ! Therefore, this model is only correct when the swimmer is roughly linear.
    do c_i = 1, num_particles_loc
      if(particles(part_ind(c_i))%particle_index_gl .eq. 1) then
        force1_d = -swimmer_force_mag1 * cos(6.28318530718 * time / swimmer_period1 + swimmer_phase1) * dvec13 / d13
      else if(particles(part_ind(c_i))%particle_index_gl .eq. 2) then
        force2_d = -swimmer_force_mag2 * cos(6.28318530718 * time / swimmer_period2 + swimmer_phase2) * dvec13 / d13
      else if(particles(part_ind(c_i))%particle_index_gl .eq. 3) then
        force3_d = -swimmer_force_mag3 * cos(6.28318530718 * time / swimmer_period3 + swimmer_phase3) * dvec13 / d13
      end if
    end do
    
    ! Compute spring forces.
    force1_s = swimmer_modulus12 * (d12 - swimmer_d12) * dvec12 / d12
    force3_s = -swimmer_modulus23 * (d23 - swimmer_d23) * dvec23 / d23
    force2_s = -force1_s - force3_s
    
    ! Compute bending forces.
    force1_b = swimmer_kb * cross_product(dvec12, cross) / d12
    force3_b = swimmer_kb * cross_product(dvec23, cross) / d23
    force2_b = -force1_b - force3_b

    ! Compute anchor forces.
    force1_a(:) = 0.
    force2_a(:) = 0.
    force3_a(:) = 0.
    
    if(swimmer_anchor .eq. .true.) then
      pos_ref(2:3) = (/ 0.5 * tny, 0.5 * tnz /)
      d_anchor1(2:3) = pos1(2:3) - pos_ref(2:3)
      d_anchor2(2:3) = pos2(2:3) - pos_ref(2:3)
      d_anchor3(2:3) = pos3(2:3) - pos_ref(2:3)
      force1_a(2:3) = -swimmer_ka * d_anchor1(2:3)
      force2_a(2:3) = -swimmer_ka * d_anchor2(2:3)
      force3_a(2:3) = -swimmer_ka * d_anchor3(2:3)
    end if
    
    ! Compute total forces.
    force1 = force1_s + force1_b + force1_d + force1_a
    force2 = force2_s + force2_b + force2_d + force2_a
    force3 = force3_s + force3_b + force3_d + force3_a
!     force1 = force1_s + force1_b + force1_a
!     force2 = force2_s + force2_b + force2_a
!     force3 = force3_s + force3_b + force3_a

    ! Add forces to beads.
    ! Forces are equally distributed to all nodes in a bead.
    ! For simplicity, all forces are added to the strain forces (force_s) of the nodes.
    do c_i = 1, num_particles_loc
      if(particles(part_ind(c_i))%particle_index_gl .eq. 1) then
        do n_i = 1, particles(part_ind(c_i))%num_nodes
          particles(part_ind(c_i))%node(n_i)%force_s = particles(part_ind(c_i))%node(n_i)%force_s + force1 / particles(part_ind(c_i))%num_nodes
        end do
      else if(particles(part_ind(c_i))%particle_index_gl .eq. 2) then
        do n_i = 1, particles(part_ind(c_i))%num_nodes
          particles(part_ind(c_i))%node(n_i)%force_s = particles(part_ind(c_i))%node(n_i)%force_s + force2 / particles(part_ind(c_i))%num_nodes
        end do
      else if(particles(part_ind(c_i))%particle_index_gl .eq. 3) then
        do n_i = 1, particles(part_ind(c_i))%num_nodes
          particles(part_ind(c_i))%node(n_i)%force_s = particles(part_ind(c_i))%node(n_i)%force_s + force3 / particles(part_ind(c_i))%num_nodes
        end do
      end if
    end do

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("finished subroutine compute_swimmer_spring_forces")
#endif
  end subroutine compute_swimmer_spring_forces
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update of look-up table for node-node interactions
  !>
  !> The look-up table for the node positions is updated.
  !> This table is used for node-node interaction forces.
  !> Each node is sorted into a subbox in which it is located.
  !> A halo region is used for efficient force computation.

  subroutine update_node_interaction_LUT()
    ! Declare variables.
    integer :: n_i ! node index
    integer :: pos_x, pos_y, pos_z ! position of node on lattice
    integer :: ierror ! MPI error

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("starting subroutine update_node_interaction_LUT...", .false.)
#endif

    ! Reset look-up table before updating it.
    LUT_nodes_in_range = 0

    ! Update table.
    ! The integer positions are defined in such a way that
    ! 1) The node is in the physical region if 1 .le. pos_x .le. nx
    ! 2) The node is in one of the halo regions if pos_x .lt. 1 or pos_x .gt. nx
    ! The same holds for the other directions.
    ! NOTE: The physical domain starts at 0.5.
    do n_i = 1, num_nodes_local
      pos_x = floor(nodes_local(n_i)%pos(1) + 0.5 - ccoords(1) * nx)
      pos_y = floor(nodes_local(n_i)%pos(2) + 0.5 - ccoords(2) * ny)
      pos_z = floor(nodes_local(n_i)%pos(3) + 0.5 - ccoords(3) * nz)

      ! If the node is located in the physical+halo region of the subdomain, the counter is increased and the node index is stored.
      if((pos_x .ge. 1 - HALO_WIDTH) .and. (pos_x .le. nx + HALO_WIDTH) .and. &
        & (pos_y .ge. 1 - HALO_WIDTH) .and. (pos_y .le. ny + HALO_WIDTH) .and. &
        & (pos_z .ge. 1 - HALO_WIDTH) .and. (pos_z .le. nz + HALO_WIDTH)) then
        LUT_nodes_in_range(pos_x, pos_y, pos_z, 0) = LUT_nodes_in_range(pos_x, pos_y, pos_z, 0) + 1

        ! Check if the LUT is sufficiently large to accept another node.
        ! If not, the simulation is terminated.
        if(LUT_nodes_in_range(pos_x, pos_y, pos_z, 0) .gt. MAX_NODES_SUBBOX) then
          call error("too many nodes in lattice subbox: particle radius or MAX_NODES_SUBBOX too small or particle overlap")
        end if

        LUT_nodes_in_range(pos_x, pos_y, pos_z, LUT_nodes_in_range(pos_x, pos_y, pos_z, 0)) = n_i
      end if
    end do

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("finished subroutine update_node_interaction_LUT", .false.)
#endif

  end subroutine update_node_interaction_LUT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Rebuild interior index
  !>
  !> If in IBM_INDEXFIELD mode, the index field has to be initialized
  !> once at the beginning or restart of the simulation.
  !> It may also be rebuilt from time to time in order to avoid numerical artifacts accumulating.
  !> The index field denotes fluid nodes as being inside or outside of a particle.
  !> It takes the value 0 if the site is completely outside and 1 if it is completely inside.
  !> Sites in the particle surface neighborhood have interpolated values between 0 and 1.

#ifdef IBM_INDEXFIELD
  subroutine rebuild_interior_index(N, force_rebuild)
    type(lbe_site), dimension(0:, 0:, 0:), intent(inout) :: N !< lattice
    logical, intent(in) :: force_rebuild !< enforces rebuild even if time step does not match

    ! Declare variables.
    integer :: x, y, z ! lattice coordinates
    integer :: i ! counter
    integer :: ierror ! MPI error
    logical :: recompute_field ! if index field has to be recomputed completely
    integer, allocatable, dimension(:) :: counter ! used to count the number of neighboring sites of a cluster
    integer, allocatable, dimension(:) :: cluster_indices ! list of cluster indices
    integer :: cluster_num ! number of clusters
    real(kind=rk), allocatable, dimension(:) :: surface_index ! used to identify inside/outside state of a cluster
    integer, allocatable, dimension(:, :, :) :: clusters ! cluster field
    integer :: num_clusters_max ! maximum number of clusters allowed
    real(kind=rk) :: visc_ext, visc_int ! external and internal viscosities
    real(kind=rk) :: tau_local ! local relaxation parameter

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("starting subroutine rebuild_interior_index...", .false.)
#endif

    ! Store old index field.
#ifdef NOSURFACTANT
  interior_index_old(:, :, :) = interior_index(:, :, :)
#endif

    ! Decide whether to update or to recompute.
    if((mod(nt, time_step_rebuild_index) == 0) .or. (force_rebuild .eqv. .true.)) then
      recompute_field = .true.
    else
      recompute_field = .false.
    end if

    ! It is assumed that the IBM spreading subroutine has been performed before.
    ! Note: This assumption is not checked here.
    ! During that routine, all distances of the lattice nodes close to the particle surfaces have been computed.
    ! From this distance field, the interior index field has to be computed first.
    ! This can only be done at lattice nodes which are close to particle surfaces.
    do x = 1, nx
      do y = 1, ny
        do z = 1, nz
          if(nint(checked_index(x, y, z)) == 1) then
            interior_index(x, y, z) = set_interior_index(normal_distance(x, y, z))
          end if
        end do
      end do
    end do

    ! Now, all lattice nodes which are not in direct vicinity of the particle surfaces have to be considered.
    ! There are two possibilities:
    ! 1) The interior index field is updated from a former valid field.
    ! 2) The interior index field is recomputed from scratch.

    ! Case 1: Update from a former field.
    ! Here, it is assumed that the particle surfaces move slowly (less than 0.5 lattice nodes per time step).
    ! Therefore, the new interior index takes the value
    ! 0 if the index of this lattice node was <  0.5 before,
    ! 1 if the index of this lattice node was >= 0.5 before.
    if(recompute_field .eqv. .false.) then
      do x = 1, nx
        do y = 1, ny
          do z = 1, nz
            if(nint(checked_index(x, y, z)) == 0) then
              if(interior_index(x, y, z) < 0.5d0) then
                interior_index(x, y, z) = 0.0d0
              else
                interior_index(x, y, z) = 1.0d0
              end if
            end if
          end do
        end do
      end do
    end if

    ! Case 2: Recompute from scratch.
    ! In this case, first all connected clusters of fluid sites are identified.
    ! This is done with an extended Hoshen-Kopelman algorithm.
    ! Once connected cluster have been identified, it has to be checked whether a given cluster is inside or outside.

    ! Perform Hoshen-Kopelman algorithm to identify clusters of connected lattice sites and count the number of clusters.
    if(recompute_field .eqv. .true.) then
      ! The maximum number of clusters is the number of lattice nodes.
      ! Note: It will never be as many clusters. But this way, not unexpected problems are caused.
      num_clusters_max = nx * ny * nz
      allocate(clusters(1:nx, 1:ny, 1:nz))
      clusters(:, :, :) = 0
      call hoshen_kopelman(checked_index, clusters, -0.5d0, 0.5d0)

      ! Count the number of clusters other than the original cluster.
      allocate(cluster_indices(num_clusters_max))
      cluster_indices(:) = 0
      cluster_num = 0

      do x = 1, nx
        do y = 1, ny
          do z = 1, nz
            ! Consider only clusters other than the original cluster
            if(clusters(x, y, z) > 0) then
              ! Check if the cluster index at the current site has already been detected.
              ! Otherwise, remember it and proceed with the next lattice site.
              ! Eventually, the number of clusters is identified.
              do i = 1, num_clusters_max
                if(cluster_indices(i) == clusters(x, y, z)) then
                  exit
                else if(cluster_indices(i) == 0) then
                  cluster_num = cluster_num + 1
                  cluster_indices(cluster_num) = clusters(x, y, z)
                  exit
                end if
              end do
            end if
          end do
        end do
      end do

      ! The cluster indices are not required any more, memory can be deallocated again.
      deallocate(cluster_indices)
    end if

    ! Each cluster has to be checked for its location (inside or outside).
    ! If a cluster is outside, all its lattice sites receive interior_index = 0,
    ! if a cluster is inside, all its lattice sites receive interior_index = 1.
    ! The following algorithm detects common boundaries of the clusters with the original cluster
    ! (i.e., the particle surface neighborhood).
    ! The average value of the index field along the cluster surface is computed.
    ! If the average is closer to 0, the cluster is outside.
    ! If it is closer to 1, the cluster is inside.
    if(recompute_field .eqv. .true.) then
      allocate(surface_index(num_clusters_max))
      allocate(counter(num_clusters_max))
      surface_index(:) = 0.0d0
      counter(:) = 0

      do x = 1, nx
        do y = 1, ny
          do z = 1, nz
            if((clusters(x, y, z) > 0)) then
              if((x - 1 >= 1) .and. (clusters(x - 1, y, z) == 0)) then
                counter(clusters(x, y, z)) = counter(clusters(x, y, z)) + 1
                surface_index(clusters(x, y, z)) = surface_index(clusters(x, y, z)) + interior_index(x - 1, y, z)
              end if
              if((x + 1 <= nx) .and. (clusters(x + 1, y, z) == 0)) then
                counter(clusters(x, y, z)) = counter(clusters(x, y, z)) + 1
                surface_index(clusters(x, y, z)) = surface_index(clusters(x, y, z)) + interior_index(x + 1, y, z)
              end if
              if((y - 1 >= 1) .and. (clusters(x, y - 1, z) == 0)) then
                counter(clusters(x, y, z)) = counter(clusters(x, y, z)) + 1
                surface_index(clusters(x, y, z)) = surface_index(clusters(x, y, z)) + interior_index(x, y - 1, z)
              end if
              if((y + 1 <= ny) .and. (clusters(x, y + 1, z) == 0)) then
                counter(clusters(x, y, z)) = counter(clusters(x, y, z)) + 1
                surface_index(clusters(x, y, z)) = surface_index(clusters(x, y, z)) + interior_index(x, y + 1, z)
              end if
              if((z - 1 >= 1) .and. (clusters(x, y, z - 1) == 0)) then
                counter(clusters(x, y, z)) = counter(clusters(x, y, z)) + 1
                surface_index(clusters(x, y, z)) = surface_index(clusters(x, y, z)) + interior_index(x, y, z - 1)
              end if
              if((z + 1 <= nz) .and. (clusters(x, y, z + 1) == 0)) then
                counter(clusters(x, y, z)) = counter(clusters(x, y, z)) + 1
                surface_index(clusters(x, y, z)) = surface_index(clusters(x, y, z)) + interior_index(x, y, z + 1)
              end if
            end if
          end do
        end do
      end do

      ! Normalize surface index.
      do i = 1, num_clusters_max
        if(counter(i) > 0) then
          surface_index(i) = surface_index(i) / counter(i)
        end if
      end do

      ! Set the clusters to inside or outside.
      ! If the counter is zero, there is no surface and the cluster has to be outside.
      ! Reason: Particles are not allowed to be larger than the subdomain size.
      ! Therefore, a subdomain without particle surface has to be outside of any particle.
      ! The same holds if the number of clusters is unity.
      do x = 1, nx
        do y = 1, ny
          do z = 1, nz
            if(clusters(x, y, z) > 0) then
              if(cluster_num == 1) then
                interior_index(x, y, z) = 0.0d0
              else
                if(surface_index(clusters(x, y, z)) > 0.5) then
                  interior_index(x, y, z) = 1.0d0
                else
                  interior_index(x, y, z) = 0.0d0
                end if
              end if
            end if
          end do
        end do
      end do

      ! Clean up memory.
      deallocate(surface_index)
      deallocate(counter)
      deallocate(clusters)
    end if

    ! Update relaxation parameter of the fluid.
#ifdef VARTAU
    visc_ext = tau_r - 0.5d0
    visc_int = viscosity_contrast * visc_ext

    do x = 1, nx
      do y = 1, ny
        do z = 1, nz
          if(interior_index(x, y, z) .eq. 0.0d0) then
            tau_local = tau_r
          else if (interior_index(x, y, z) .eq. 1.0d0) then
            tau_local = visc_int + 0.5d0
          else
            tau_local = visc_ext * (1.0d0 - interior_index(x, y, z)) + visc_int * interior_index(x, y, z) + 0.5d0
          end if

          call set_local_tau_r(N, x, y, z, tau_local)
        end do
      end do
    end do
#endif

!     ! Recolor fluid if index field changes.
! #ifdef NOSURFACTANT
!     do x = 1, nx
!       do y = 1, ny
!         do z = 1, nz
!           if((interior_index(x, y, z) .lt. 0.5d0) .and. (interior_index_old(x, y, z) .gt. 0.5d0)) then
!             call recolor_fluid(N, x, y, z, get_neighboring_color_ratio(N, x, y, z))
!           else if ((interior_index(x, y, z) .gt. 0.5d0) .and. (interior_index_old(x, y, z) .lt. 0.5d0)) then
! 	        call recolor_fluid(N, x, y, z, 0.5d0)
!           end if
!         end do
!       end do
!     end do
! #endif

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("finished subroutine rebuild_interior_index", .false.)
#endif
  end subroutine rebuild_interior_index
! endif IBM_INDEXFIELD
#endif

! endif IBM_PART
#endif

end module lsuperstruct_timeloop_module
