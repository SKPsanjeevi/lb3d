!> Lagrangian superstructure interface module
!>
!> This module contains all subroutines required for explicit interaction of the superstructure and LB3D.

#include "lbe.h"

module lsuperstruct_interface_module

#ifdef IBM_PART

  ! Include external modules.
#ifdef IBM_BINARYIBM
   use lbe_bdist_module, only: boltz_dist
   use lbe_force_interface_module, only: add_force_to_n_r, add_force_to_n_b
#endif
  use lbe_force_interface_module, only: add_force_to_all, get_total_force
  use lbe_log_module
  use lbe_parallel_module, only: ccoords, tnx, tny, tnz, comm_cart
  use lbe_parms_module, only: nx, ny, nz, boundary_width, boundary_cond
  use lbe_globals_module, only: nvecs, g, c, lbe_force
  use lbe_types_module, only: lbe_site
  use lsuperstruct_data_module
  use lsuperstruct_parallel_module, only: exchange_rock_grad_halo

  implicit none
  include 'mpif.h'

  private
  public :: compute_fluid_velocity, set_rock_grad
#ifdef IBM_INDEXFIELD
  public :: compute_selective_wall_forces, get_ibm_index
#endif
#if defined IBM_INDEXFIELD && defined VARTAU
  public :: set_local_tau_r
#endif
#ifdef IBM_BINARYIBM
   public :: recolour_interior_sites
   public :: compute_multicomponent_forces
   public :: recolor_fluid
   public :: get_neighboring_color_ratio
#endif

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get local red and blue densities
  !>
  !> The densities are computed from the populations.
  !> TODO: These functions should be completely redundant eventually (take precomputed values or LB3D functionality).

  function get_density_red(N, x, y, z)
    real(kind=rk) :: get_density_red !< return value
    type(lbe_site), dimension(0:, 0:, 0:), intent(in) :: N !< lattice
    integer :: x, y, z !< lattice coordinates

    ! Compute red density.
    get_density_red = sum(N(x, y, z)%n_r(:) * g(:))
  end function get_density_red

#ifndef SINGLEFLUID
  function get_density_blue(N, x, y, z)
    real(kind=rk) :: get_density_blue !< return value
    type(lbe_site), dimension(0:, 0:, 0:), intent(in) :: N !< lattice
    integer :: x, y, z !< lattice coordinates

    ! Compute blue density.
    get_density_blue = sum(N(x, y, z)%n_b(:) * g(:))
  end function get_density_blue
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute physical fluid velocity on the entire lattice
  !>
  !> The physical fluid velocity is computed on the entire lattice.
  !> The fluid velocity is the first moment of the populations, corrected by half of the force.
  !> It is checked at compiletime if SINGLEFLUID is defined:
  !> - if true: One component is used, the normal velocity computation is executed.
  !> - if false: Two components are used, the physical velocity has to be computed according to the Shan-Chen approach.
  !> TODO: This subroutine should be completely redundant eventually (take precomputed values or LB3D functionality).

  subroutine compute_fluid_velocity(N)
    type(lbe_site), dimension(0:, 0:, 0:), intent(in) :: N !< lattice

    ! Declare variables.
    integer :: x, y, z ! lattice coordinates
    real(kind=rk) :: den, den_red, den_blue ! densities
    real(kind=rk), dimension(3) :: mom, mom_bare_red, mom_bare_blue ! momenta
    real(kind=rk), dimension(3) :: force, force_red, force_blue ! force correction for velocity
		
    do z = 1, nz
      do y = 1, ny
        do x = 1, nx
          ! Consider fluid sites first.
          if(N(x, y, z)%rock_state .eq. 0.0d0) then
            ! Compute densities, momenta, and velocities.
#ifdef SINGLEFLUID
            call get_total_force(force(1), force(2), force(3), x, y, z)
            den = get_density_red(N, x, y, z)
            mom(1) = sum(N(x, y, z)%n_r(:) * g(:) * c(:, 1)) + 0.5d0 * force(1)
            mom(2) = sum(N(x, y, z)%n_r(:) * g(:) * c(:, 2)) + 0.5d0 * force(2)
            mom(3) = sum(N(x, y, z)%n_r(:) * g(:) * c(:, 3)) + 0.5d0 * force(3)
            vel_phys(:, x, y, z) = mom(:) / den
#else
#ifdef IBM_BINARYIBM
            den_red = get_density_red(N, x, y, z)
            den_blue = get_density_blue(N, x, y, z)
            den = den_red + den_blue
            force_red(:) = lbe_force(:, 1, x, y, z) + lbe_force(:, 0, x, y, z) * den_red / den
            force_blue(:) = lbe_force(:, 2, x, y, z) + lbe_force(:, 0, x, y, z) * den_blue / den
            mom_bare_red(1) = sum(N(x, y, z)%n_r(:) * g(:) * c(:, 1))
            mom_bare_red(2) = sum(N(x, y, z)%n_r(:) * g(:) * c(:, 2))
            mom_bare_red(3) = sum(N(x, y, z)%n_r(:) * g(:) * c(:, 3))
            mom_bare_blue(1) = sum(N(x, y, z)%n_b(:) * g(:) * c(:, 1))
            mom_bare_blue(2) = sum(N(x, y, z)%n_b(:) * g(:) * c(:, 2))
            mom_bare_blue(3) = sum(N(x, y, z)%n_b(:) * g(:) * c(:, 3))
            vel_phys(:, x, y, z) = (mom_bare_red(:) + mom_bare_blue(:) + 0.5d0 * (force_red(:) + force_blue(:))) / den
            vel_red(:, x, y, z) = mom_bare_red(:) / den_red
            vel_blue(:, x, y, z) = mom_bare_blue(:) / den_blue
#endif
#endif

            ! Compute the difference velocity of blue and red in the two-component case.
#ifdef IBM_BINARYIBM
            vel_diff(:, x, y, z) = vel_red(:, x, y, z) - vel_blue(:, x, y, z)
#endif
          ! Consider rock sites.
          else
            vel_phys(:, x, y, z) = 0.0d0
#ifdef IBM_BINARYIBM
            vel_red(:, x, y, z) = 0.0d0
            vel_blue(:, x, y, z) = 0.0d0
            vel_diff(:, x, y, z) = 0.0d0
#endif
          end if
        end do
      end do
    end do
  end subroutine compute_fluid_velocity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Recolour fluid sites inside IBM particles
  !>
  !> All fluid sites are recoloured according to the index function:
  !> - inside (index = 1): densities are redistributed to achieve zero colour
  !> - outside (index = 0): densities are not touched
  !> - boundary (0 < index < 1): only a fraction of the densities is recoloured
  !> TODO: This function is likely to get kicked out soon.

#ifdef IBM_BINARYIBM
   subroutine recolour_interior_sites(N)
     type(lbe_site), dimension(0:, 0:, 0:), intent(inout) :: N !< lattice
 
     ! Declare variables.
     integer :: x, y, z ! lattice coordinates
     real(kind=rk), dimension(nvecs) :: temp_r, temp_b ! buffer for populations
     integer :: ierror ! MPI error code
		 ! Declare variables for 2nd part
		 !integer :: x, y, z ! lattice coordinates
     real (kind=rk) :: local_index
     real (kind=rk) :: local_color

 
     ! Debugging routine
#ifdef IBM_DEBUG
     call MPI_Barrier(comm_cart, ierror)
     call log_msg("starting subroutine recolour_interior_sites...", .false.)
#endif
     call log_msg("starting subroutine recolour_interior_sites...", .false.)
 
     do z = 1, nz
       do y = 1, ny
         do x = 1, nx
           temp_r(:) = N(x, y, z)%n_r(:)
           temp_b(:) = N(x, y, z)%n_b(:)
           N(x, y, z)%n_r(:) = temp_r(:) * (1.0d0 - interior_index(x, y, z)) + 0.5d0 * (temp_r(:) + temp_b(:)) * interior_index(x, y, z)
           N(x, y, z)%n_b(:) = temp_b(:) * (1.0d0 - interior_index(x, y, z)) + 0.5d0 * (temp_r(:) + temp_b(:)) * interior_index(x, y, z)
         end do
       end do
     end do
 
     ! Debugging routine
#ifdef IBM_DEBUG
     call MPI_Barrier(comm_cart, ierror)
     call log_msg("finished subroutine recolour_interior_sites", .false.)
#endif
     call log_msg("finished subroutine recolour_interior_sites", .false.)
 
     !do z = 1, nz
     !  do y = 1, ny
     !    do x = 1, nx
     !      ! Get interior index.
     !      local_index = interior_index(x, y, z)
     !      
     !      if(local_index > 0.5d0) then
     !        local_color = get_neighboring_color_ratio(N, x, y, z)
     !      end if
     !      
     !      if(local_color .ge. 0.0d0) then
     !        call recolor_fluid(N, x, y, z, local_color)
     !      end if
     !    end do
     !  end do
     !end do
 
  end subroutine recolour_interior_sites
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set rock gradient
  !>
  !> For distributing the particles throughout the numerical domain, possible rock nodes pose a problem.
  !> Instead of avoiding rock nodes in the first place, one can distribute the particles assuming that there are no rocks.
  !> Afterwards, the particles are pushed out of the rock regions.
  !> For this, a rock node 'gradient' comes in handy. It defines the direction where the particles have to move
  !> in order to reach the fluid.
  !> A simple gradient can be constructed based on the onion principle:
  !> An index indicates how many rock layers are between the current node and the closest fluid region.
  !> This index can be found iteratively and in parallel.

  subroutine set_rock_grad(N)
    type(lbe_site), dimension(0:, 0:, 0:), intent(in) :: N !< lattice

    ! Declare variables.
    integer :: num_iter ! number of iteration steps
    integer :: max_iter_loc ! local maximum number of iterations required
    integer :: i ! counter
    integer :: x, y, z ! coordinates
    integer :: sx, sy, sz ! coordinate shifts for identifying neighbors
    integer :: ierror ! MPI error code
    logical :: is_in_bulk ! helper logical for identifying bulk state of the lattice node
    logical :: any_increase ! helper logical for identifying if a new layer has been identified
    character(len=200) :: message ! message string

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("starting subroutine set_rock_grad...", .false.)
#endif

    ! Report.
    call log_msg("finding rock gradient...", .false.)

    ! Set number of iterations.
    ! The maximum number of iterations is given by half of the largest system extension.
    num_iter = max(tnx, tny, tnz) / 2
    max_iter_loc = 0

    ! Set initial configuration: 0 for fluid, 2 for rock.
    ! The 2 is important since the value 1 will be used for the outermost fluid layer as additional buffer zone.
    do x = 1 - IBM_HALO, nx + IBM_HALO
      do y = 1 - IBM_HALO, ny + IBM_HALO
        do z = 1 - IBM_HALO, nz + IBM_HALO
          if(N(x, y, z)%rock_state .eq. 0.0d0) then
            rock_grad(x, y, z) = 0
          else
            rock_grad(x, y, z) = 2
            max_iter_loc = 2
          end if
        end do
      end do
    end do

    ! Use one additional fluid layer as buffer zone.
    ! This is to make sure that the particles still see a gradient when they already in the fluid.
    ! If at least one next neighbor has rock index 2, the node is set to 1.
    do x = 1, nx
      do y = 1, ny
        do z = 1, nz
          if(rock_grad(x, y, z) .eq. 0) then
            do sx = -1, 1
              do sy = -1, 1
                do sz = -1, 1
                  if(rock_grad(x + sx, y + sy, z + sz) .eq. 2) then
                    rock_grad(x, y, z) = 1
                  end if
                end do
              end do
            end do
          end if
        end do
      end do
    end do

    ! Exchange rock_grad halo with neighboring processes.
    call exchange_rock_grad_halo()

    ! Iterate and identify remaining layers.
    do i = 2, num_iter
      ! Assume first that no new layer has been found.
      any_increase = .false.

      ! Run over the lattice and identify new layers.
      do x = 1, nx
        do y = 1, ny
          do z = 1, nz
            ! Skip the current lattice node if its index is too low.
            if(rock_grad(x, y, z) .lt. i) cycle

            ! Assume that all neighboring nodes have an index which is not lower than the iteration index.
            is_in_bulk = .true.

            ! If at least one neighbor has an index which is too small, the node is not part of a new layer.
            ! In that case, the logical is set to false.
            do sx = -1, 1
              do sy = -1, 1
                do sz = -1, 1
                  if(rock_grad(x + sx, y + sy, z + sz) .lt. i) is_in_bulk = .false.
                end do
              end do
            end do

            ! If the node is in a new layer, increase its index.
            if(is_in_bulk) then
              rock_grad(x, y, z) = rock_grad(x, y, z) + 1
              any_increase = .true.
            end if
          end do
        end do
      end do

      ! Increase the counter of layers found so far.
      if(any_increase) max_iter_loc = i + 1

      ! Exchange the rock gradient.
      call exchange_rock_grad_halo()
    end do

    ! Find the globally maximum layer index of the rock region.
    call MPI_Allreduce(max_iter_loc, rock_grad_max, 1, MPI_INTEGER, MPI_MAX, comm_cart, ierror)

    ! Compute the rock gradient based on a simple FD scheme.
    do x = 1, nx
      do y = 1, ny
        do z = 1, nz
          gradient(1, x, y, z) = 0.5d0 * (rock_grad(x + 1, y, z) - rock_grad(x - 1, y, z))
          gradient(2, x, y, z) = 0.5d0 * (rock_grad(x, y + 1, z) - rock_grad(x, y - 1, z))
          gradient(3, x, y, z) = 0.5d0 * (rock_grad(x, y, z + 1) - rock_grad(x, y, z - 1))
        end do
      end do
    end do

    ! Report.
    write(message, "('  found maximum layer = ', i0)") rock_grad_max
    call log_msg(trim(message), .false.)

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("finished subroutine set_rock_grad", .false.)
#endif
  end subroutine set_rock_grad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute selective wall forces
  !>
  !> If in IBM_INDEXFIELD mode, there may be a selective wall force acting on the fluid.
  !> The force shall repell or attract only the fluid outside the particles.
  !> Therefore, the force is sensitive to the value of the interior/exterior index field.
  !> It is locally weighted by (1 - interior_index).
  !> TODO:
  !> - Currently, the force is only set up for a situation with two planar walls (boundary_cond = 6).
  !>   In all other cases, the subroutine will not have any effect.

#ifdef IBM_INDEXFIELD
  subroutine compute_selective_wall_forces()

    ! Declare variables.
    integer :: x, y, z ! local lattice coordinates
    integer :: x_gl ! global lattice coordinates
    integer :: ierror ! MPI error code
    real(kind=rk) :: x_dist_bot, x_dist_top ! x-distances between lattice site and bottom and top walls
    real(kind=rk) :: force_strength ! strength of the force
    real(kind=rk) :: force ! force

    ! Check execution condition.
    ! There will not be any effect if the geometry is not Poieuille-like (boundary_cond = 6).
    ! For efficiency, this subroutine is also skipped if either the force magnitude or the range is zero.
    if((boundary_cond .ne. 6) .or. (selective_wall_force_magnitude .eq. 0.0d0) .or. (selective_wall_force_range .eq. 0.0d0)) return

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("starting subroutine compute_selective_wall_forces...", .false.)
#endif

    ! Run over the lattice and find selective force.
    do x = 1, nx
      ! Compute global x-position and distances to bottom and top walls.
      x_gl = x + ccoords(1) * nx
      x_dist_bot = x_gl - boundary_width - 0.5d0
      x_dist_top = tnx - boundary_width - x_gl + 0.5d0

      ! Consider bottom wall region.
      ! Only lattice sites within the interaction range are taken into account.
      if((x_dist_bot >= 0.0d0) .and. (x_dist_bot <= selective_wall_force_range)) then
        force_strength = selective_wall_force_magnitude * 1.5819767d0 * (exp(-x_dist_bot / selective_wall_force_range) &
            & - 0.3678794d0)

        do y = 1, ny
          do z = 1, nz
            force = force_strength * (1.0d0 - interior_index(x, y, z))
            call add_force_to_all(force, 0.0d0, 0.0d0, x, y, z)
          end do
        end do
      end if

      ! Consider top wall region.
      ! Only lattice sites within the interaction range are taken into account.
      if((x_dist_top >= 0.0d0) .and. (x_dist_top <= selective_wall_force_range)) then
        force_strength = selective_wall_force_magnitude * 1.5819767d0 * (exp(-x_dist_top / selective_wall_force_range) &
            & - 0.3678794d0)

        do y = 1, ny
          do z = 1, nz
            force = force_strength * (1.0d0 - interior_index(x, y, z))
            call add_force_to_all(-force, 0.0d0, 0.0d0, x, y, z)
          end do
        end do
      end if
    end do

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("finished subroutine compute_selective_wall_forces", .false.)
#endif
  end subroutine compute_selective_wall_forces
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute no-slip forces for multicomponent fluid inside IBM particles
  !>
  !> A friction force between red and blue components is computed.
  !> This force is designed in such a way that it adds momentum to both components resulting in a relative no-slip of blue and red.
  !> The force is only active in regions inside the particle, i.e., it is weighted by the index field.
  !> TODO: This functionality is likely to be strongly modified.

#ifdef IBM_BINARYIBM 
	subroutine compute_multicomponent_forces(N)
     type(lbe_site), dimension(0:, 0:, 0:), intent(in) :: N !< lattice
 
     ! Declare variables.
     integer :: x, y, z ! lattice coordinates
     real(kind=rk), dimension(1:3) :: force ! no-slip force
     real(kind=rk) :: den_red, den_blue ! component densities
     real(kind=rk) :: colour ! fluid colour
     real(kind=rk), dimension(3) :: colour_grad ! fluid colour gradient
 
     do x = 1, nx
       do y = 1, ny
         do z = 1, nz
           ! Compute friction force
           den_red = get_density_red(N, x, y, z)
           den_blue = get_density_blue(N, x, y, z)
           force(:) = vel_diff(:, x, y, z) * interior_index(x, y, z) * den_red * den_blue / (den_red + den_blue)
           !call add_force_to_n_r(-force(1), -force(2), -force(3), x, y, z)
           !call add_force_to_n_b(+force(1), +force(2), +force(3), x, y, z)
           
           ! Compute colour force
           !if(interior_index(x, y, z) .eq. 1.d0) then
           !  call compute_colour_gradient(N, x, y, z, colour_grad)
           !  colour = den_red - den_blue
           !  force(:) = noslipcoeff * abs(colour) * interior_index(x, y, z) * colour_grad(:)
           !  call add_force_to_n_r(-force(1), -force(2), -force(3), x, y, z)
           !  call add_force_to_n_b(+force(1), +force(2), +force(3), x, y, z)
           !end if
         end do
       end do
     end do
 
!     ! Declare variables.
!     integer :: x, y, z ! lattice coordinates
!     integer :: i ! velocity index
!     real(kind=rk) :: den_r, den_b, den ! component densities and total density
!     real(kind=rk), dimension(1:3) :: force_r, force_b ! component forces
!     real(kind=rk), dimension(1:3) :: momentum_r, momentum_b ! red and blue momenta
!     real(kind=rk), dimension(1:3) :: correction_r, correction_b ! correction forces
!    
!     do x = 1, nx
!       do y = 1, ny
!         do z = 1, nz
!           ! Get densities.
!           den_r = get_density_red(N, x, y, z)
!           den_b = get_density_blue(N, x, y, z)
!           den = den_r + den_b
! 
!           ! Get corrent pre-collision momenta.
!           momentum_r(:) = 0.d0
!           momentum_b(:) = 0.d0
! 
!           do i = 1, nvecs
!             momentum_r(:) = momentum_r(:) + N(x, y, z)%n_r(i) * g(i) * c(i, :)
!             momentum_b(:) = momentum_b(:) + N(x, y, z)%n_b(i) * g(i) * c(i, :)
!           end do
! 
!           ! Obtain component forces.
!           force_r(:) = lbe_force(:, 1, x, y, z) + lbe_force(:, 0, x, y, z) * den_r / den
!           force_b(:) = lbe_force(:, 2, x, y, z) + lbe_force(:, 0, x, y, z) * den_b / den
! 
!           ! Compute post-collision momenta.
!           momentum_r(:) = momentum_r(:) + force_r(:)
!           momentum_b(:) = momentum_b(:) + force_b(:)
! 
!           ! Compute correction forces.
!           correction_r(:) = (den_r * momentum_b(:) - den_b * momentum_r(:)) / den * interior_index(x, y, z)
!           correction_b(:) = -correction_r(:)
!           ! Add correction forces.
!           call add_force_to_n_r(correction_r(1), correction_r(2), correction_r(3), x, y, z)
!           call add_force_to_n_b(correction_b(1), correction_b(2), correction_b(3), x, y, z)
!         end do
!       end do
!     end do
   end subroutine compute_multicomponent_forces
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the gradient of the colour field
  !>
  !> TODO: Will probably be removed soon.

! #if defined IBM_INDEXFIELD && !defined SINGLEFLUID && defined NOSURFACTANT
!   subroutine compute_colour_gradient(N, x, y, z, colour_grad)
!     type(lbe_site), dimension(0:, 0:, 0:), intent(in) :: N !< lattice
!     integer, intent(in) :: x, y, z !< lattice coordinates
!     real(kind=rk), dimension(3), intent(out) :: colour_grad !< colour gradient
!     
!     colour_grad(1) = 0.5d0 * (get_density_red(N, x + 1, y, z) - get_density_blue(N, x + 1, y, z) &
!                     & - (get_density_red(N, x - 1, y, z) - get_density_blue(N, x - 1, y, z)))
!     colour_grad(2) = 0.5d0 * (get_density_red(N, x, y + 1, z) - get_density_blue(N, x, y + 1, z) &
!                     & - (get_density_red(N, x, y - 1, z) - get_density_blue(N, x, y - 1, z)))
!     colour_grad(3) = 0.5d0 * (get_density_red(N, x, y, z + 1) - get_density_blue(N, x, y, z + 1) &
!                     & - (get_density_red(N, x, y, z - 1) - get_density_blue(N, x, y, z - 1)))
!   end subroutine compute_colour_gradient
! #endif
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update local relaxation time in LB3D
  !>
  !> This is a wrapper function to overwrite the local relaxation time of the red fluid.
  !> This functionality is required for spatially and temporally varying viscosities.

#if defined IBM_INDEXFIELD && defined VARTAU
  subroutine set_local_tau_r(N, x, y, z, tau_local)
    type(lbe_site), dimension(0:, 0:, 0:), intent(inout) :: N !< lattice
    integer, intent(in) :: x, y, z !< lattice site coordinates
    real(kind=rk), intent(in) :: tau_local !< input relaxation time

    ! Update relaxation time.
    N(x, y, z)%taupos_r = tau_local
  end subroutine set_local_tau_r
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Recolor fluid
  !>
  !> The red and blue components of a fluid site are recolored according to the desired ratio:
  !> 'ratio' is the fraction of the red component, values between 0 and 1 are allowed.
  !> The total density and momentum are conserved.
  !> TODO: This functionality is likely to be removed or strongly modified.

#ifdef IBM_BINARYIBM
   subroutine recolor_fluid(N, x, y, z, ratio)
     type(lbe_site), dimension(0:, 0:, 0:), intent(inout) :: N !< lattice
     integer, intent(in) :: x, y, z !< lattice site coordinates
     real(kind=rk), intent(in) :: ratio !< desired fluid ratio
 
     ! Declare variables.
     real(kind=rk) :: pop(nvecs) ! populations
 			
		 !type(lbe_site), dimension(0:, 0:, 0:), intent(inout) :: N !< lattice
     !integer, intent(in) :: x, y, z !< lattice site coordinates
     !real(kind=rk), intent(in) :: ratio !< desired fluid ratio
 
     ! Declare variables.
     real(kind=rk) :: den_r, den_b, den ! local densities
     real(kind=rk), dimension(1:3) :: vel ! local velocity
     !real(kind=rk) :: pop(nvecs) ! populations

     ! Compute total populations.
     pop(:) = N(x, y, z)%n_r(:) + N(x, y, z)%n_b(:)
 
     ! Overwrite red and blue populations.
     N(x, y, z)%n_r(:) = pop(:) * ratio
     N(x, y, z)%n_b(:) = pop(:) * (1.0d0 - ratio)
     
      
     ! Compute densities.
     den_r = sum(N(x, y, z)%n_r(:) * g(:))
     den_b = sum(N(x, y, z)%n_b(:) * g(:))
     den = den_r + den_b
 
     ! Compute bare velocity.
     vel(1) = sum((N(x, y, z)%n_r(:) + N(x, y, z)%n_b(:)) * g(:) * c(:, 1)) / den
     vel(2) = sum((N(x, y, z)%n_r(:) + N(x, y, z)%n_b(:)) * g(:) * c(:, 2)) / den
     vel(3) = sum((N(x, y, z)%n_r(:) + N(x, y, z)%n_b(:)) * g(:) * c(:, 3)) / den
 
     ! Compute equilibrium.
     call boltz_dist(vel(1), vel(2), vel(3), 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, pop)
 
     ! Overwrite red and blue populations.
     N(x, y, z)%n_r(:) = pop(:) * den * ratio
     N(x, y, z)%n_b(:) = pop(:) * den * (1.0d0 - ratio)
 
   end subroutine recolor_fluid
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get neighboring fluid color ratio
  !>
  !> TODO: This functionality is likely to be removed.

#ifdef IBM_BINARYIBM
  function get_neighboring_color_ratio(N, x, y, z)
    real(kind=rk) :: get_neighboring_color_ratio !< return value
    type(lbe_site), dimension(0:, 0:, 0:), intent(in) :: N !< lattice
    integer, intent(in) :: x, y, z !< lattice site coordinates

    ! Declare variables.
    integer :: xn, yn, zn ! neighboring lattice coordinates
    integer :: counter ! number of neighboring fluid sites taken for averaging
    real(kind=rk) :: den_r, den_b ! local red and blue density

    ! Compute average neighboring fluid ratio and count number of sites.
    get_neighboring_color_ratio = 0.0d0
    counter = 0

    do xn = x - 1, x + 1
      do yn = y - 1, y + 1
        do zn = z - 1, z + 1
          ! Valid site: interior index is more fluid than particle.
          if(interior_index(xn, yn, zn) < 0.5d0) then
            counter = counter + 1
            den_r = sum(N(xn, yn, zn)%n_r(:) * g(:))
            den_b = sum(N(xn, yn, zn)%n_b(:) * g(:))
            get_neighboring_color_ratio = get_neighboring_color_ratio + den_r / (den_r + den_b)
          end if
        end do
      end do
    end do

    ! Normalize by number of sites.
    if(counter .gt. 0) then
      get_neighboring_color_ratio = get_neighboring_color_ratio / counter
    else
      get_neighboring_color_ratio = -1
    end if
  end function get_neighboring_color_ratio
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> returns local IBM inside/outside index
  !>
  !> \param[in] x lattice x position (local coordinates)
  !> \param[in] y lattice y position (local coordinates)
  !> \param[in] z lattice z position (local coordinates)
  !>
  !> Returns value between \c 0.0 (outside of all IBM particles) and \c 1.0 (inside of an IBM particle).
  !> If compiled without \c IBM_INDEXFIELD, the return value is always \c 0.0 (outside).

  function get_ibm_index(x, y, z)
    real(kind=rk) :: get_ibm_index
    integer,intent(in) :: x, y, z

#ifdef IBM_INDEXFIELD
    get_ibm_index = interior_index(x, y, z)
#else
    get_ibm_index = 0.0_rk
#endif
  end function get_ibm_index

! endif IBM_PART
#endif

end module lsuperstruct_interface_module
