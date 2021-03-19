!> Lagrangian superstructure helper module
!>
!> This module contains all subroutines with useful tools for the Lagrangian superstructure.

#include "lbe.h"

module lsuperstruct_helper_module

#ifdef IBM_PART

  ! Include external modules.
  use lbe_bdist_module, only : equilibrium_distribution
  use lbe_globals_module, only : rk, myrankc, nvecs, g, c, lbe_force
  use lbe_log_module
  use lbe_parallel_module, only : tnx, tny, tnz, ccoords, comm_cart
  use lbe_parms_module, only : get_tau_r, nx, ny, nz
  use lbe_types_module, only : lbe_site
  use lsuperstruct_data_module
  use lextobj_module, only : particles, part_ind, init_particle, deinit_particle, compute_volume_slice
  use lmesh_module, only : mesh_struct, meshes

  implicit none

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute shortest distance between two points
  !>
  !> For two given position vectors, the shortest connection is computed.
  !> The connection is evaluated as a vector pointing from the first to the second point.
  !> Periodicity is taken into account.
  !> NOTE: The position vectors have to be global positions.
  !> TODO: Periodicty may have to be checked.

  function distance_vector(pos_1, pos_2)
    real(kind=rk), dimension(3) :: distance_vector !< return value
    real(kind=rk), dimension(3), intent(in) :: pos_1, pos_2 !< input position vectors

    ! Compute distance vector.
    distance_vector = pos_2 - pos_1

    ! Correct for periodicity in the simulation.
!    if(periodicity(1) .eqv. .true.) then
      if(distance_vector(1) .gt. 0.5 * tnx) then
        distance_vector(1) = distance_vector(1) - tnx
      else if(distance_vector(1) .lt. -0.5 * tnx) then
        distance_vector(1) = distance_vector(1) + tnx
      end if
!    end if

!    if(periodicity(2) .eqv. .true.) then
      if(distance_vector(2) .gt. 0.5 * tny) then
        distance_vector(2) = distance_vector(2) - tny
      else if(distance_vector(2) .lt. -0.5 * tny) then
        distance_vector(2) = distance_vector(2) + tny
      end if
!    end if

!    if(periodicity(3) .eqv. .true.) then
      if(distance_vector(3) .gt. 0.5 * tnz) then
        distance_vector(3) = distance_vector(3) - tnz
      else if(distance_vector(3) .lt. -0.5 * tnz) then
        distance_vector(3) = distance_vector(3) + tnz
      end if
!    end if
  end function distance_vector

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Add local particle
  !>
  !> A particle is added to the local particle array.
  !> 1) It is checked whether enough memory is available. Abort simulation if not.
  !> 2) The number of particles is increased.
  !> 3) The user interface is updated.
  !> 4) The particle is initialized.
  !> After returning, the particle is ready to use.

  subroutine add_particle(mesh)
    type(mesh_struct), intent(in) :: mesh !< element of the mesh array

    ! Declare variables.
    integer :: c_i ! particle index
    integer :: c_new ! particle index of new particle
    character(len=200) :: message ! message string

    ! Report creation of particle if in IBM debug mode.
#ifdef IBM_DEBUG
    write(message, "('adding particle of mesh ', i0)") mesh%mesh_type
    call log_msg(trim(message), .true.)
#endif

    ! Check memory.
    if(num_particles_loc .eq. NUM_PART_LOC_MAX) then
      call error("too many particles in at least one subdomain")
    end if

    ! Increase the number of local particles.
    num_particles_loc = num_particles_loc + 1

    ! Find the first empty slot in the particle array.
    ! The corresponding particle index is added to the user interface.
    do c_i = 1, NUM_PART_LOC_MAX
      if(particles(c_i)%currently_used .eq. 0) then
        c_new = c_i
        part_ind(num_particles_loc) = c_new
        exit
      end if
    end do

    ! Initialize the particle.
    call init_particle(particles(c_new), mesh)
  end subroutine add_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Remove local particle
  !>
  !> A particle is removed from the local particle array.
  !> 1) It is checked whether the particle is active. If not, nothing happens.
  !> 2) The particle is deinitialized.
  !> 3) The user interface is updated.
  !> 4) The number of particles is decreased.
  !> After returning, the particle is not available any more.

  subroutine remove_particle(c_rem)
    integer, intent(in) :: c_rem !< index of particle to be removed

    ! Declare variables.
    character(len=200) :: message ! message string

    ! Report destruction of particle if in IBM debug mode.
#ifdef IBM_DEBUG
    write(message, "('removing particle of mesh ', i0)") particles(part_ind(c_rem))%mesh_type
    call log_msg(trim(message), .true.)
#endif

    ! Return if there is no particle to remove.
    if(part_ind(c_rem) .eq. -1) then
      return
    end if

    ! Deinitialize the particle.
    call deinit_particle(particles(part_ind(c_rem)))

    ! Update the user interface.
    ! The index at the latest position in the interface is shifted to that of the removed particle.
    ! This way, there is no gap.
    part_ind(c_rem) = part_ind(num_particles_loc)
    part_ind(num_particles_loc) = -1

    ! Decrease the number of local particles.
    num_particles_loc = num_particles_loc - 1
  end subroutine remove_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compactify local particle list
  !>
  !> After removing local particles, there may be gaps in the local particle list.
  !> These gaps are removed by redistributing the particles in the local list.
  !> WARNING: This subroutine is experimental and currently unused.

  subroutine compactify_particle_list
    ! Declare variables.
    integer :: first_gap ! location of first gap found
    integer :: last_particle ! location of last particle found
    integer :: c_i ! counter
    logical :: swaps_left ! whether swaps are still required

    ! Start with initial values.
    first_gap = NUM_PART_LOC_MAX + 1
    last_particle = 0
    swaps_left = .true.

    do while(swaps_left .eqv. .true.)
      ! Find the first gap.
      do c_i = 1, NUM_PART_LOC_MAX
	if(part_ind(c_i) .lt. 0) then
	  first_gap = c_i
	  exit
	end if
      end do

      ! Find the last particle.
      do c_i = NUM_PART_LOC_MAX, 1, -1
	if(part_ind(c_i) .gt. 0) then
	  last_particle = c_i
	  exit
	end if
      end do

      if(last_particle .gt. first_gap) then
	! Swap the first gap and the last particle and mark the new gap.
	part_ind(first_gap) = part_ind(last_particle)
	part_ind(last_particle) = -1
      else
	swaps_left = .false.
      end if
    end do
  end subroutine compactify_particle_list

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute volume fraction in bins
  !>
  !> The entire geometry is subdivided into \tnx bins (indices in interval [1, tnx]) with unit thickness along the x-axis.
  !> For each of these bins, the particle volume fraction, averaged over the y- and z-axes is computed.
  !> Two neighboring bins (indices n, n+1) are separated by a slice (index n+1).
  !> Slice 1 is identical to the plane bounding the simulation volume at the lower end of the x-axis.

  subroutine calculate_particle_volume_fraction(vol_frac_loc)
    real(kind=rk), dimension(1:tnx), intent(inout) :: vol_frac_loc !< volume fraction profile

    ! Declare variables.
    integer :: c_i ! particle index
    integer :: n ! counter
    integer :: bin_start, bin_end ! first and last bin in which a particle is located.
    integer :: ierror ! MPI error code
    real(kind=rk), dimension(3) :: slice_normal = (/1.d0, 0.d0, 0.d0/) ! normal of the plane cutting the particle
    real(kind=rk), dimension(3) :: slice_origin ! origin of the plane cutting the particle
    real(kind=rk), dimension(tnx) :: slice_volume ! volume of particles in normal direction of the slice

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("starting subroutine calculate_particle_volume_fraction...", .false.)
#endif

    ! Reset volume fraction.
    vol_frac_loc(:) = 0.d0

    ! Loop over all particles in the subdomain.
    do c_i = 1, num_particles_loc
      ! Reset slice volume.
      slice_volume(:) = 0.d0

      ! Identify the lateral extension of the particle along the x-axis.
      ! \c bin_start and \c bin_end are the first and last bins, respectively, in which the particle is located.
      ! Both \c bin_start and \c bin_end are located within the interval [1, tnx].
      bin_start = 1 ! TODO: find value
      bin_end = tnx ! TODO: find value

      ! Compute the slice volumes of the particle for each slice cutting the particle.
      ! This gives the part of the volume of the particle which is located in direction \c slice_normal
      ! of the slice located at \c slice_origin.
      ! Slices are located between the lattice sites and, therefore, have half-integer positions.
      do n = bin_start, bin_end
        slice_origin = (/n - 0.5d0, 0.d0, 0.d0/)
        slice_volume(n) = compute_volume_slice(particles(part_ind(c_i)), slice_origin, slice_normal)
      end do

      ! Compute volume of the particle in every bin and update the bin volume fraction.
      ! The bin volume is identical to the volume of the particle between two adjacent slices which in turn
      ! is the difference of two successive slice volumes, except for the last bin:
      ! In that case, the bin volume equals the slice volume.
      do n = bin_start, bin_end - 1
        vol_frac_loc(n) = vol_frac_loc(n) + (slice_volume(n) - slice_volume(n+1))
      end do

      vol_frac_loc(bin_end) = vol_frac_loc(bin_end) + slice_volume(bin_end)
    end do

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("finished subroutine calculate_particle_volume_fraction", .false.)
#endif
  end subroutine calculate_particle_volume_fraction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute fluid stress profile
  !>
  !> This subroutine computes the local fluid deviatoric stress profile (along the x-axis).
  !> The algorithm bases on the nonequilibrium distributions.
  !> The resulting stress profile is written to \c str_fluid_loc.
  !> TODO:
  !> - Currently, only the stress for the red fluid is computed, based on the assumption
  !>   that the fluid is a single-component fluid.
  !>   It would be nice to generalize this to multi-component systems.
  !>   However, it is not clear at this point, how the stress tensor is defined in auch a situation.

  subroutine calculate_fluid_stresses(str_fluid_loc, N)
    real(kind=rk), dimension(1:9, 1:tnx), intent(inout) :: str_fluid_loc !< fluid stress profile
    type(lbe_site), dimension(0:, 0:, 0:), intent(in) :: N !< lattice

    ! Declare variables.
    integer :: i ! index
    integer :: x_loc, y_loc, z_loc ! local coordinates
    integer :: x_gl ! global coordinate
    integer :: ierror ! MPI error code
    real(kind=rk) :: den ! local fluid density (only red component!)
    real(kind=rk), dimension(3) :: vel ! local fluid velocity
    real(kind=rk), dimension(3) :: force ! local force density
    real(kind=rk), dimension(nvecs) :: f_eq, f_neq ! equilibrium and non-equilibrium populations (only red component!)
    real(kind=rk), dimension(3, 3) :: stress, force_correction ! temporary stress tensor and force correction
    character(len=200) :: message ! message string

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("starting subroutine calculate_fluid_stresses...", .false.)
#endif

    ! Compatibility warning
#ifndef SINGLEFLUID
    write(message, "('lsuperstruct_helper=>calculate_fluid_stresses: warning, currently only working properly for SINGLEFLUID!')")
    call log_msg(trim(message), .false.)
#endif

    ! Reset fluid stress array.
    str_fluid_loc(:, :) = 0.d0

    ! Loop over local lattice and compute averages.
    do x_loc = 1, nx
      ! Find global x-position.
      x_gl = x_loc + nx * ccoords(1)

      do y_loc = 1, ny
        do z_loc = 1, nz
          ! Reset temporary stress.
          stress(:, :) = 0.0d0

          ! Compute local fluid density.
          den = 0.0d0

          do i = 1, nvecs
            den = den + N(x_loc, y_loc, z_loc)%n_r(i) * g(i)
          end do

          ! Compute local equilibrium distribution.
          vel(:) = vel_phys(:, x_loc, y_loc, z_loc)
          force(:) = lbe_force(:, 0, x_loc, y_loc, z_loc)
          call equilibrium_distribution(den, vel_phys(:, x_loc, y_loc, z_loc), f_eq)
          f_neq(:) = N(x_loc, y_loc, z_loc)%n_r(:) * g(:) - f_eq(:)

          ! Compute local stress tensor.
          do i = 1, nvecs
            stress(1, 1) = stress(1, 1) + f_neq(i) * c(i, 1) * c(i, 1)
            stress(1, 2) = stress(1, 2) + f_neq(i) * c(i, 1) * c(i, 2)
            stress(1, 3) = stress(1, 3) + f_neq(i) * c(i, 1) * c(i, 3)
            stress(2, 1) = stress(2, 1) + f_neq(i) * c(i, 2) * c(i, 1)
            stress(2, 2) = stress(2, 2) + f_neq(i) * c(i, 2) * c(i, 2)
            stress(2, 3) = stress(2, 3) + f_neq(i) * c(i, 2) * c(i, 3)
            stress(3, 1) = stress(3, 1) + f_neq(i) * c(i, 3) * c(i, 1)
            stress(3, 2) = stress(3, 2) + f_neq(i) * c(i, 3) * c(i, 2)
            stress(3, 3) = stress(3, 3) + f_neq(i) * c(i, 3) * c(i, 3)
          end do

          ! Compute force correction.
          force_correction(1, 1) = 0.5d0 * (vel(1) * force(1) + force(1) * vel(1))
          force_correction(1, 2) = 0.5d0 * (vel(1) * force(2) + force(1) * vel(2))
          force_correction(1, 3) = 0.5d0 * (vel(1) * force(3) + force(1) * vel(3))
          force_correction(2, 1) = 0.5d0 * (vel(2) * force(1) + force(2) * vel(1))
          force_correction(2, 2) = 0.5d0 * (vel(2) * force(2) + force(2) * vel(2))
          force_correction(2, 3) = 0.5d0 * (vel(2) * force(3) + force(2) * vel(3))
          force_correction(3, 1) = 0.5d0 * (vel(3) * force(1) + force(3) * vel(1))
          force_correction(3, 2) = 0.5d0 * (vel(3) * force(2) + force(3) * vel(2))
          force_correction(3, 3) = 0.5d0 * (vel(3) * force(3) + force(3) * vel(3))

          ! Normalize stress.
          stress(:, :) = (stress(:, :) + force_correction(:, :)) * (0.5d0 / get_tau_r(N, x_loc, y_loc, z_loc) - 1.0d0) / (tny * tnz)

          ! Update stress profile.
          str_fluid_loc(1, x_gl) = str_fluid_loc(1, x_gl) + stress(1, 1)
          str_fluid_loc(2, x_gl) = str_fluid_loc(2, x_gl) + stress(1, 2)
          str_fluid_loc(3, x_gl) = str_fluid_loc(3, x_gl) + stress(1, 3)
          str_fluid_loc(4, x_gl) = str_fluid_loc(4, x_gl) + stress(2, 1)
          str_fluid_loc(5, x_gl) = str_fluid_loc(5, x_gl) + stress(2, 2)
          str_fluid_loc(6, x_gl) = str_fluid_loc(6, x_gl) + stress(2, 3)
          str_fluid_loc(7, x_gl) = str_fluid_loc(7, x_gl) + stress(3, 1)
          str_fluid_loc(8, x_gl) = str_fluid_loc(8, x_gl) + stress(3, 2)
          str_fluid_loc(9, x_gl) = str_fluid_loc(9, x_gl) + stress(3, 3)
        end do
      end do
    end do

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("finished subroutine calculate_fluid_stresses", .false.)
#endif
  end subroutine calculate_fluid_stresses

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Find interior index from normal distance
  !>
  !> The interior index depends on the distance.
  !> - For distances > +DIST_CUTOFF, the index is set to 0 (exterior).
  !> - For distances < -DIST_CUTOFF, the index is set to 1 (interior).
  !> - For all other distances, the index is linearly interpolated.

  function set_interior_index(distance)
    real(kind=rk) :: set_interior_index !< return value
    real(kind=rk), intent(in) :: distance !< input distance

    if(distance < -DIST_CUTOFF) then
      set_interior_index = 1.0d0
    else if(distance > DIST_CUTOFF) then
      set_interior_index = 0.0d0
    else
      set_interior_index = 0.5d0 * (DIST_CUTOFF - distance) / DIST_CUTOFF
    end if
  end function set_interior_index

! endif IBM_PART
#endif

end module lsuperstruct_helper_module
