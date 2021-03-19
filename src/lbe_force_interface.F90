#include "lbe.h"

!> wrapper functions for applying and reading forces locally on the
!> system
module lbe_force_interface_module
    use lbe_globals_module

    implicit none
    private

    public add_force_to_all, &
          & add_force_to_n_r, &
#ifndef SINGLEFLUID
          & add_force_to_n_b, &
#endif
#ifndef NOSURFACTANT
          & add_force_to_n_s, &
#endif
          & get_total_force

contains

    !> Wrapper function for adding force to a lattice node
    !>
    !> Adds a force (fx, fy, fz) to \c lbe_force(1:3, 0, x, y, z) at position (x, y, z).
    !> This force affects everything located at the lattice node, irrespective of its rock state or color density.
    subroutine add_force_to_all(fx, fy, fz, x, y, z)
      real(kind=rk), intent(in) :: fx, fy, fz !< force to be added to lbe_force
      integer, intent(in) :: x, y, z !< position of lattice node where the force is added

      ! Add force to lbe_force.
      lbe_force(1, 0, x, y, z) = lbe_force(1, 0, x, y, z) + fx
      lbe_force(2, 0, x, y, z) = lbe_force(2, 0, x, y, z) + fy
      lbe_force(3, 0, x, y, z) = lbe_force(3, 0, x, y, z) + fz
    end subroutine add_force_to_all

    !> Wrapper function for adding force to the red component of a lattice node
    !>
    !> Adds a force (fx, fy, fz) to \c lbe_force(1:3, 1, x, y, z) at position (x, y, z).
    !> This force affects only the red component of the lattice node, irrespective of its rock state.
    subroutine add_force_to_n_r(fx, fy, fz, x, y, z)
      real(kind=rk), intent(in) :: fx, fy, fz !< force to be added to lbe_force
      integer, intent(in) :: x, y, z !< position of lattice node where the force is added

      ! Add force to lbe_force.
      lbe_force(1, 1, x, y, z) = lbe_force(1, 1, x, y, z) + fx
      lbe_force(2, 1, x, y, z) = lbe_force(2, 1, x, y, z) + fy
      lbe_force(3, 1, x, y, z) = lbe_force(3, 1, x, y, z) + fz
    end subroutine add_force_to_n_r

    !> Wrapper function for adding force to the blue component of a lattice node
    !>
    !> Adds a force (fx, fy, fz) to \c lbe_force(1:3, 2, x, y, z) at position (x, y, z).
    !> This force affects only the blue component of the lattice node, irrespective of its rock state.
    !> This function is only defined if SINGLEFLUID is not selected.
#ifndef SINGLEFLUID
    subroutine add_force_to_n_b(fx, fy, fz, x, y, z)
      real(kind=rk), intent(in) :: fx, fy, fz !< force to be added to lbe_force
      integer, intent(in) :: x, y, z !< position of lattice node where the force is added

      ! Add force to lbe_force.
      lbe_force(1, 2, x, y, z) = lbe_force(1, 2, x, y, z) + fx
      lbe_force(2, 2, x, y, z) = lbe_force(2, 2, x, y, z) + fy
      lbe_force(3, 2, x, y, z) = lbe_force(3, 2, x, y, z) + fz
    end subroutine add_force_to_n_b
#endif

    !> Wrapper function for adding force to the surfactant of a lattice node
    !>
    !> Adds a force (fx, fy, fz) to \c lbe_force(1:3, 3, x, y, z) at position (x, y, z).
    !> This force affects only the surfactant of the lattice node, irrespective of its rock state.
    !> This function is only defined if NOSURFACTANT is not selected.
#ifndef NOSURFACTANT
    subroutine add_force_to_n_s(fx, fy, fz, x, y, z)
      real(kind=rk), intent(in) :: fx, fy, fz !< force to be added to lbe_force
      integer, intent(in) :: x, y, z !< position of lattice node where the force is added

      ! Add force to lbe_force.
      lbe_force(1, 3, x, y, z) = lbe_force(1, 3, x, y, z) + fx
      lbe_force(2, 3, x, y, z) = lbe_force(2, 3, x, y, z) + fy
      lbe_force(3, 3, x, y, z) = lbe_force(3, 3, x, y, z) + fz
    end subroutine add_force_to_n_s
#endif

    !> Wrapper function for accessing the total force acting on a lattice node
    !>
    !> Writes the total force acting on a lattice node at position (x, y, z) into (fx, fy, fz).
    !> The total force is the sum of the forces stored in the (n_spec + 1) components of \c lbe_force(1:3, 0:n_spec, x, y, z).
    subroutine get_total_force(fx, fy, fz, x, y, z)
      real(kind=rk), intent(out) :: fx, fy, fz !< force acting on the lattice node
      integer, intent(in) :: x, y, z !< position of lattice node where the force is taken from

      ! Compute total force.
      fx = sum(lbe_force(1, 0:n_spec, x, y, z))
      fy = sum(lbe_force(2, 0:n_spec, x, y, z))
      fz = sum(lbe_force(3, 0:n_spec, x, y, z))
    end subroutine get_total_force

end module lbe_force_interface_module
