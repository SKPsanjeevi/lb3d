!> Lagrangian superstructure immmersed boundary method module
!>
!> This module contains all subroutines responsible for the immersed boundary method.
!> Additionally, the IBM spreading is used to set up the index field indicating interior/exterior regions.

#include "lbe.h"

module lsuperstruct_IBM_module

#ifdef IBM_PART

  ! Include external modules.
  use lbe_force_interface_module, only : add_force_to_all
  use lbe_helper_module, only : norm
  use lbe_log_module
  use lbe_parallel_module, only : comm_cart, ccoords
  use lbe_parms_module, only : nx, ny, nz
  use lbe_types_module, only : lbe_site
  use lsuperstruct_data_module
  use lextobj_module, only : particles, part_ind

  implicit none
  private
  public :: IBM_interpolate_velocities, IBM_spread_forces
  
  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Interpolate velocities
  !>
  !> The particle node velocities are interpolated from the fluid velocities.
  !> TODO: Tidy up after testing! This is a mess.

  subroutine IBM_interpolate_velocities()
    ! Declare variables.
    integer :: n_i ! node index
    integer :: X, Y, Z ! counters
    integer, dimension(3) :: fluid_pos ! fluid node position
    integer, dimension(3) :: fluid_low ! lowest fluid node position
    real(kind=rk), dimension(3) :: node_pos ! particle node position
    real(kind=rk), dimension(3) :: dist_low ! distance between lattice and particle node
    real(kind=rk), dimension(3) :: vel_intpol ! interpolated node velocity
    real(kind=rk), dimension(3, IBM_RANGE) :: weight ! interpolation weights
    integer :: ierror ! MPI error

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("starting subroutine IBM_interpolate_velocities...", .false.)
#endif

    ! Interpolate velocities.
    ! Run over all nodes which are located in the physical range of the process.
    do n_i = 1, num_nodes_local
      ! Only consider nodes in the physical regime.
      if(nodes_local(n_i)%physical .eq. 1) then
        ! Compute the node position in the local process.
        node_pos(1) = nodes_local(n_i)%pos(1) - ccoords(1) * nx
        node_pos(2) = nodes_local(n_i)%pos(2) - ccoords(2) * ny
        node_pos(3) = nodes_local(n_i)%pos(3) - ccoords(3) * nz

        ! Find lowest fluid node neighbor.
        ! This is the fluid node with the smallest coordinates still in interpolation range of the node.
        ! For fluid nodes in the halo, the values are 0 or negative.
        fluid_low(:) = floor(node_pos(:) - 0.5 * IBM_RANGE + 1)

        ! Compute distance of particle node to lowest neighbor.
        dist_low(:) = fluid_low(:) - node_pos(:)

        ! Compute interpolation weigths.
        ! The lowest distance is passed, and the results are stored in the weight array.
        call interpolation_weights(dist_low, weight)

        ! Reset interpolation velocity.
        vel_intpol(:) = 0

        ! Interpolate velocities from neighbors.
        do X = 1, IBM_RANGE
          do Y = 1, IBM_RANGE
            do Z = 1, IBM_RANGE
              ! Identify the position of the corresponding fluid node.
              fluid_pos(1) = fluid_low(1) + X - 1
              fluid_pos(2) = fluid_low(2) + Y - 1
              fluid_pos(3) = fluid_low(3) + Z - 1

              ! Get fluid velocity.
              ! The velocity contribution is added to the interpolated velocity via a wrapper subroutine.
              ! This subroutine also checks whether the target fluid node is fluid or rock.
              vel_intpol = vel_intpol + vel_phys(1:3, fluid_pos(1), fluid_pos(2), fluid_pos(3)) &
                * weight(1, X) * weight(2, Y) * weight(3, Z)
            end do
          end do
        end do

        ! Update interpolated velocity.
        nodes_local(n_i)%vel = vel_intpol
      end if
    end do

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("finished subroutine IBM_interpolate_velocities", .false.)
#endif

  end subroutine IBM_interpolate_velocities

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Spread forces and build interior index
  !>
  !> The particle node forces are spread to the fluid forces.
  !> At the same time, the index indicating interior/exterior regions of the particles is updated.
  !> Combining both approaches avoids additional lattice sweeps.
  !> TODO: Check relevance of force_IBM.

  subroutine IBM_spread_forces(spread_force)
    logical, intent(in) :: spread_force !< whether to actually spread the force

    ! Declare variables.
    integer :: n_i ! node index
    integer :: X, Y, Z ! counters
    integer, dimension(3) :: pos ! fluid node position
    integer, dimension(3) :: fluid_low ! lowest fluid node position
    real(kind=rk) :: dist ! distance between particle and lattice node
    real(kind=rk), dimension(3) :: node_pos ! particle node position
    real(kind=rk), dimension(3) :: dist_low ! distance between lowest fluid node and particle node
    real(kind=rk), dimension(3) :: force ! spread node force
    real(kind=rk), dimension(3, IBM_RANGE) :: weight ! interpolation weights
    real(kind=rk), dimension(3) :: center
    integer :: ierror ! MPI error

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("starting subroutine IBM_spread_forces...", .false.)
#endif

    ! Reset the IBM force array and the distance field to the membrane.
    ! The lattice force is initialized as zero force.
    ! The distance field is set to IBM_RANGE.
    force_IBM(:, :, :, :) = 0.0d0
#ifdef IBM_INDEXFIELD
    normal_distance(:, :, :) = 1000.0d0
    checked_index(:, :, :) = 0.0d0
#endif

    ! Spread forces.
    ! Run over all nodes which are located in the physical range of the process.
    do n_i = 1, num_nodes_local
      ! Compute the node position in the local process.
      node_pos(1) = nodes_local(n_i)%pos(1) - ccoords(1) * nx
      node_pos(2) = nodes_local(n_i)%pos(2) - ccoords(2) * ny
      node_pos(3) = nodes_local(n_i)%pos(3) - ccoords(3) * nz

      ! Only consider nodes in the physical+halo regime.
      if(node_pos(1) .ge. 0.5 - HALO_WIDTH .and. node_pos(1) .lt. nx + 0.5d0 + HALO_WIDTH .and. &
        & node_pos(2) .ge. 0.5 - HALO_WIDTH .and. node_pos(2) .lt. ny + 0.5d0 + HALO_WIDTH .and. &
        & node_pos(3) .ge. 0.5 - HALO_WIDTH .and. node_pos(3) .lt. nz + 0.5d0 + HALO_WIDTH) then
        ! Find lowest fluid node neighbor.
        ! This is the fluid node with the smallest coordinates still in interpolation range of the node.
        ! For fluid nodes in the halo, the values are 0 or negative.
        fluid_low(:) = floor(node_pos(:) - 0.5 * IBM_RANGE + 1)

        ! Compute distance of particle node to lowest neighbor.
        dist_low(:) = fluid_low(:) - node_pos(:)

        ! Compute interpolation weigths.
        ! The lowest distance is passed, and the results are stored in the weight array.
        call interpolation_weights(dist_low, weight)

        ! Spread force to neighbors.
        ! Make sure that only physical fluid nodes receive the force
        do X = 1, IBM_RANGE
          ! Compute x-component of target fluid node.
          pos(1) = fluid_low(1) + X - 1

          ! Continue only if target fluid node is physical.
          if((pos(1) .lt. 1) .or. (pos(1) .gt. nx)) cycle

          do Y = 1, IBM_RANGE
            ! Compute y-component of target fluid node.
            pos(2) = fluid_low(2) + Y - 1

            ! Continue only if target fluid node is physical.
            if((pos(2) .lt. 1) .or. (pos(2) .gt. ny)) cycle

            do Z = 1, IBM_RANGE
              ! Compute z-component of target fluid node.
              pos(3) = fluid_low(3) + Z - 1

              ! Continue only if target fluid node is physical.
              if((pos(3) .lt. 1) .or. (pos(3) .gt. nz)) cycle

              ! Compute the force contribution.
              ! It is the total particle node force weighted by the interpolation stencil for the corresponding fluid node.
              force = nodes_local(n_i)%force_tot * weight(1, X) * weight(2, Y) * weight(3, Z)

              ! Update the lattice force.
              ! The force contribution is added to the lattice force via a wrapper subroutine.
              if(spread_force .eqv. .true.) then
                call add_force_to_all(force(1), force(2), force(3), pos(1), pos(2), pos(3))
                force_IBM(1:3, pos(1), pos(2), pos(3)) = force_IBM(1:3, pos(1), pos(2), pos(3)) + force(1:3)
              end if

              ! Compute normal distance between particle and lattice node and update the distance field.
              ! Update the distance if its magnitude is smaller than the one already in memory.
              ! This way, the shortest distance is found eventually.
              ! \c dist is negative if the lattice node is inside the particle, positive otherwise.
#ifdef IBM_INDEXFIELD
              ! Use the local membrane curvature to find a better estimate for the distance.
              ! If the curvature radius is small, assume a locally spherical surface.
              if(abs(nodes_local(n_i)%curv_radius) < 10.0d0) then
                ! This is the center position of the sphere.
                center(:) = node_pos(:) - nodes_local(n_i)%curv_radius * nodes_local(n_i)%normal(:)

                ! Distinguish between a locally concave and convex surface.
                if(nodes_local(n_i)%curv_radius > 0.) then
                  dist = norm(pos - center) - nodes_local(n_i)%curv_radius ! convex
                else
                  dist = -norm(pos - center) - nodes_local(n_i)%curv_radius ! concave
                end if
              ! Otherwise, assume a flat interface.
              else
                dist = dot_product((pos - node_pos), nodes_local(n_i)%normal)
              end if

              ! Update the distance if its magnitude is smaller than before.
              if(abs(dist) < abs(normal_distance(pos(1), pos(2), pos(3)))) then
                normal_distance(pos(1), pos(2), pos(3)) = dist
                checked_index(pos(1), pos(2), pos(3)) = 1.0d0
              end if
#endif
            end do
          end do
        end do
      end if
    end do

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("finished subroutine IBM_spread_forces", .false.)
#endif

  end subroutine IBM_spread_forces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Interpolation weight
  !>
  !> The IBM interpolation weight is computed based on a given distance between particle and fluid node.
  !> The interpolation functions are symmetric: the magnitude of the distance is sufficient.
  !> The shape of the interpolation function depends on its range, there are three different algorithms implemented:
  !> 1) range 2
  !> 2) range 3
  !> 3) range 4
  !> The subroutine has been designed in such a way that it is numerically efficient (in terms of if-statements)
  !> and flexible (in terms of selecting the interpolation range as a parameter which is not hard-coded).

  subroutine interpolation_weights(dist_low, weight)
    real(kind=rk), dimension(3), intent(in) :: dist_low !< distance between particle node and lower fluid node
    real(kind=rk), dimension(3, IBM_RANGE), intent(out) :: weight !< interpolation weights

    ! Declare variables.
    integer :: coord ! coordinate index (for x, y, z)
    integer :: range ! range index (between 1 and IBM_RANGE)
    real(kind=rk) :: dist ! current distance along one axis

    ! Compute interpolation weights if the range is 2.
    if(IBM_RANGE .eq. 2) then
      ! Run over all spatial directions and all nodes within interpolation range
      do coord = 1, 3
        do range = 1, IBM_RANGE
          ! Compute the actual distance along one axis.
          dist = abs(dist_low(coord) + range - 1)

          ! Compute interpolation weights.
          if(dist .le. 1) then
            weight(coord, range) = 1 - dist
          else
            call log_msg("detected invalid distance in IBM interpolation", .true.)
            write(msgstr, "(f0.15)") dist
            call log_msg(trim(msgstr), .true.)
            weight(coord, range) = 0
          end if
        end do
      end do
    end if

    ! Compute interpolation weights if the range is 3.
    if(IBM_RANGE .eq. 3) then
      ! Run over all spatial directions and all nodes within interpolation range
      do coord = 1, 3
        do range = 1, IBM_RANGE
          ! Compute the actual distance along one axis.
          dist = abs(dist_low(coord) + range - 1)

          ! Compute interpolation weights.
          if(dist .le. 0.5) then
            weight(coord, range) = (1 + sqrt(1 - 3 * dist**2)) / 3.0
          else if(dist .le. 1.5) then
            weight(coord, range) = (5 - 3 * dist - sqrt(-2 + 6 * dist - 3 * dist**2)) / 6.0
          else
            call log_msg("detected invalid distance in IBM interpolation", .true.)
            write(msgstr, "(f0.15)") dist
            call log_msg(trim(msgstr), .true.)
            weight(coord, range) = 0
          end if
        end do
      end do
    end if

    ! Compute interpolation weights if the range is 4.
    if(IBM_RANGE .eq. 4) then
      ! Run over all spatial directions and all nodes within interpolation range
      do coord = 1, 3
        do range = 1, IBM_RANGE
          ! Compute the actual distance along one axis.
          dist = abs(dist_low(coord) + range - 1)

          ! Compute interpolation weights.
          if(dist .le. 1) then
            weight(coord, range) = (3 - 2 * dist + sqrt(1 + 4 * dist - 4 * dist**2)) / 8.0
          else if(dist .le. 2) then
            weight(coord, range) = (5 - 2 * dist - sqrt(-7 + 12 * dist - 4 * dist**2)) / 8.0
          else
            call log_msg("detected invalid distance in IBM interpolation", .true.)
            write(msgstr, "(f0.15)") dist
            call log_msg(trim(msgstr), .true.)
            weight(coord, range) = 0
          end if
        end do
      end do
    end if
  end subroutine interpolation_weights

! endif IBM_PART
#endif

end module lsuperstruct_IBM_module
