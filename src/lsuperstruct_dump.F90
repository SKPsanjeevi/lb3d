!> Lagrangian superstructure dump module
!>
!> This module contains all structures, variables, and subroutines responsible
!> for dumping data related to the Lagrangian superstructure.

#include "lbe.h"

module lsuperstruct_dump_module

#ifdef IBM_PART

  ! Include external modules.
  use lbe_globals_module, only: myrankc
  use lbe_helper_module, only: density
  use lbe_io_helper_module, only: lbe_delete_file, lbe_make_filename_cp_rank, lbe_make_filename_output
  use lbe_io_module, only: dump_scalar
  use lbe_io_xdrf_module
  use lbe_log_module
  use lbe_parallel_module, only: ccoords, check_allocate, comm_cart, tnx, tny, tnz, nprocs
  use lbe_parms_module, only: nt, nx, ny, nz
  use lbe_types_module, only: lbe_site
  use lsuperstruct_data_module
  use lsuperstruct_helper_module, only: calculate_particle_volume_fraction, calculate_fluid_stresses
  use lmesh_module, only: meshes
  use lextobj_module, only: particles, part_ind, compute_angular_velocity, compute_particle_stresslet, &
      & compute_inertia_tensor, compute_inertia_ellipsoid

  implicit none
  include 'mpif.h'
  private
  public :: dump_IBM_checkpoint, IBM_delete_checkpoint, write_particles_vtk, dump_lattice_profiles, &
            & dump_static_particle_data, write_particles_dat
#ifdef IBM_INDEXFIELD
  public :: dump_indexfield_lattice
#endif


  ! variable declarations

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write particle checkpoint
  !>
  !> The particle checkpoint is written in binary format (XDRF)

  subroutine dump_IBM_checkpoint()
    ! Declare variables.
    integer :: file_unit ! output file unit
    integer :: c_i ! particle index
    integer :: n_i ! node index
    integer :: ierror ! node index
    character(len=100) :: filename ! output filename

    ! Set variables.
    file_unit = 42

    ! Report beginning of the initialization.
    call MPI_Barrier(comm_cart, ierror) ! just to make sure that all processes have reached this point.
    call log_msg("  Dumping IBM checkpoint ...")

    ! Create and open checkpoint files.
    call lbe_make_filename_cp_rank(filename, 'checkpoint_IBM', '.xdr', nt, myrankc)
    call xdrfopen(file_unit, filename, "w", ierror)
    call check_xdrfopen(ierror, filename)

    ! Write number of local particles.
    call xdrfint(file_unit, num_particles_loc, ierror)

    ! Write all data for each local particle.
    do c_i = 1, num_particles_loc
      call xdrfint(file_unit, particles(part_ind(c_i))%mesh_type, ierror)
      call xdrfint(file_unit, particles(part_ind(c_i))%particle_index_gl, ierror)
      call xdrfint(file_unit, particles(part_ind(c_i))%num_jumps(1), ierror)
      call xdrfint(file_unit, particles(part_ind(c_i))%num_jumps(2), ierror)
      call xdrfint(file_unit, particles(part_ind(c_i))%num_jumps(3), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%radius_0, ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%radius, ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%volume_0, ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%volume, ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%surface_0, ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%surface, ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%k_s, ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%k_al, ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%k_b, ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%k_v, ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%k_at, ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%center(1), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%center(2), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%center(3), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%center_old(1), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%center_old(2), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%center_old(3), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%force_total(1), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%force_total(2), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%force_total(3), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%torque_total(1), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%torque_total(2), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%torque_total(3), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%linear_velocity(1), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%linear_velocity(2), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%linear_velocity(3), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%linear_momentum(1), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%linear_momentum(2), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%linear_momentum(3), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%angular_velocity(1), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%angular_velocity(2), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%angular_velocity(3), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%angular_momentum(1), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%angular_momentum(2), ierror)
      call xdrfdouble(file_unit, particles(part_ind(c_i))%angular_momentum(3), ierror)
#ifdef IBM_FIXED
      call xdrfdouble(file_unit, particles(part_ind(c_i))%k_anchor, ierror)
#endif
      ! Write node data.
      do n_i = 1, particles(part_ind(c_i))%num_nodes
        call xdrfdouble(file_unit, particles(part_ind(c_i))%node(n_i)%pos(1), ierror)
        call xdrfdouble(file_unit, particles(part_ind(c_i))%node(n_i)%pos(2), ierror)
        call xdrfdouble(file_unit, particles(part_ind(c_i))%node(n_i)%pos(3), ierror)
        call xdrfdouble(file_unit, particles(part_ind(c_i))%node(n_i)%pos_old(1), ierror)
        call xdrfdouble(file_unit, particles(part_ind(c_i))%node(n_i)%pos_old(2), ierror)
        call xdrfdouble(file_unit, particles(part_ind(c_i))%node(n_i)%pos_old(3), ierror)
        call xdrfdouble(file_unit, particles(part_ind(c_i))%node(n_i)%vel(1), ierror)
        call xdrfdouble(file_unit, particles(part_ind(c_i))%node(n_i)%vel(2), ierror)
        call xdrfdouble(file_unit, particles(part_ind(c_i))%node(n_i)%vel(3), ierror)
#ifdef IBM_FIXED
        	call xdrfdouble(file_unit, particles(part_ind(c_i))%node(n_i)%pos_anchor(1), ierror)
        	call xdrfdouble(file_unit, particles(part_ind(c_i))%node(n_i)%pos_anchor(2), ierror)
        	call xdrfdouble(file_unit, particles(part_ind(c_i))%node(n_i)%pos_anchor(3), ierror)
#endif
      end do
    end do

    ! Close the checkpoint file.
    call xdrfclose(file_unit, ierror)

    ! Report success of the dump.
    call MPI_Barrier(comm_cart, ierror) ! just to make sure that all processes have reached this point.
    call log_msg("  Finished dumping IBM checkpoint.")
  end subroutine dump_IBM_checkpoint

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Delete old IBM checkpoint files
  !>
  !> Depending on the number of desired checkpoints, the older checkpoints are deleted.

subroutine IBM_delete_checkpoint(last)
  integer, intent(in) :: last !< index of last checkpoint

  ! Declare variables.
  character(len=256) :: filename

  ! Report.
  call log_msg("  Deleting IBM checkpoint files ...")

  ! Checkpoint formal is always 'xdr'.
  call lbe_make_filename_cp_rank(filename, 'checkpoint_IBM', '.xdr', last, myrankc)
  call lbe_delete_file(filename)  
end subroutine IBM_delete_checkpoint

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write particle data as VTK file
  !>
  !> The particle states are written as VTK files (ASCII).
  !> Each process writes one VTK file containing the particles in its local domain.
  !> TODO:
  !> - The current solution is a bit inefficient for two reasons:
  !>   1) The output is ASCII (large files).
  !>   2) Each process writes into its own file (inefficient for large process numbers).
  !>   The idea is to prepare a parallel non-VTK output. The resulting file can then be converted to VTK.

  subroutine write_particles_vtk(directory, time)
    character(len=*), intent(in) :: directory !< folder to write data into
    integer, intent(in) :: time !< time index for VTK output

    ! Declare variables.
    integer :: file_unit ! output file unit
    integer :: c_i ! particle index
    integer :: n_i ! node index
    integer :: f_i ! face index
    integer :: offset ! offset for getting correct node-face-relations
    integer :: ierror ! MPI error code
    character(len=100) :: filename ! output filename

    ! Check time step and execution condition.
    if((mod(time, time_step_dump_vtk) .ne. 0) .or. (dump_vtk .eqv. .false.)) return

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("starting subroutine write_particles_vtk...", .false.)
#endif

    ! Set variables.
    file_unit = 42

    ! Create filename and open file.
    write(filename, "(a, '/output_p', i0, '_t', i0, '.vtk')") directory, myrankc, time
    open(unit = file_unit, file = filename)

    ! Write header.
    write(file_unit, "('# vtk DataFile Version 2.0')")
    write(file_unit, "('Particles')")
    write(file_unit, "('ASCII')")
    write(file_unit, "('DATASET POLYDATA')")

    ! Write nodes.
    write(file_unit, "('POINTS ', i0, ' float')") num_nodes_loc_vtk

    do c_i = 1, num_particles_loc
      do n_i = 1, particles(part_ind(c_i))%num_nodes
        write(file_unit, "(F0.4, ' ', F0.4, ' ', F0.4)") particles(part_ind(c_i))%node(n_i)%pos(:)
      end do
    end do

    ! Write faces.
    write(file_unit, "('POLYGONS ', i0, ' ', i0)") num_faces_loc_vtk, 4 * num_faces_loc_vtk

    offset = 0

    do c_i = 1, num_particles_loc
      do f_i = 1, particles(part_ind(c_i))%num_faces
        write(file_unit, "('3 ', i0, ' ', i0, ' ', i0)") &
          & meshes(particles(part_ind(c_i))%mesh_type)%neighbor_face_node(f_i, 1) - 1 + offset, &
          & meshes(particles(part_ind(c_i))%mesh_type)%neighbor_face_node(f_i, 2) - 1 + offset, &
          & meshes(particles(part_ind(c_i))%mesh_type)%neighbor_face_node(f_i, 3) - 1 + offset
      end do

      offset = offset + particles(part_ind(c_i))%num_nodes
    end do

    ! Write cell (face) data.
    write(file_unit, "('CELL_DATA ', i0)") num_faces_loc_vtk

    ! Write area deviations.
    write(file_unit, "('SCALARS area_deviation float')")
    write(file_unit, "('LOOKUP_TABLE default')")

    do c_i = 1, num_particles_loc
      do f_i = 1, particles(part_ind(c_i))%num_faces
        write(file_unit, "(E10.3)") particles(part_ind(c_i))%face(f_i)%area / &
          & (meshes(particles(part_ind(c_i))%mesh_type)%area(f_i) * (particles(part_ind(c_i))%radius)**2) - 1.0
      end do
    end do

    ! Write particle index.
    write(file_unit, "('SCALARS cell_index integer')")
    write(file_unit, "('LOOKUP_TABLE default')")

    do c_i = 1, num_particles_loc
      do f_i = 1, particles(part_ind(c_i))%num_faces
        write(file_unit, "(i0)") particles(part_ind(c_i))%particle_index_gl
      end do
    end do

    ! Close file.
    close(file_unit)

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("finished subroutine write_particles_vtk", .false.)
#endif
  end subroutine write_particles_vtk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write particle contour
  !>
  !> The particle contours are written as VTK files (ASCII).
  !> Each process writes one DAT file containing the particle nodes in its local domain.

  subroutine write_particles_dat
    ! Declare variables.
    integer :: file_unit ! output file unit
    integer :: c_i ! particle index
    integer :: n_i ! node index
    integer :: ierror ! MPI error code
    character(len=100) :: filename ! output filename

    ! Check time step and execution condition.
    if((mod(nt, time_step_dump_particles) .ne. 0) .or. (dump_particles .eqv. .false.)) return

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("starting subroutine write_particles_dat...", .false.)
#endif

    ! Set variables.
    file_unit = 42

    ! Create filename and open file.
    write(filename, "('data/particle_nodes_p', i0, '_t', i0, '.vtk')") myrankc, nt 
    open(unit = file_unit, file = filename)

    ! Write header.
    write(file_unit, "('# x y z')")

    ! Write nodes.
    do c_i = 1, num_particles_loc
      do n_i = 1, particles(part_ind(c_i))%num_nodes
        write(file_unit, "(F0.6, ' ', F0.6, ' ', F0.6)") particles(part_ind(c_i))%node(n_i)%pos(:)
      end do
    end do

    ! Close file.
    close(file_unit)

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("finished subroutine write_particles_dat", .false.)
#endif
  end subroutine write_particles_dat
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Dump interior/exterior index field
  !>
  !> The interior/exterior index field is dumped if the time step matches.
  !> This subroutine is only defined if in IBM_INDEXFIELD mode.

#ifdef IBM_INDEXFIELD
  subroutine dump_indexfield_lattice()
    ! Declare variables.
    integer :: ierror ! MPI error code

    ! Check time step and execution condition.
    if((mod(nt, time_step_dump_indexfield) .ne. 0) .or. (dump_indexfield .eqv. .false.)) return

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("starting subroutine dump_indexfield_lattice...", .false.)
#endif

    call dump_scalar(interior_index(1:nx, 1:ny, 1:nz), 'indexfield')

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("finished subroutine dump_indexfield_lattice", .false.)
#endif
  end subroutine dump_indexfield_lattice
! endif IBM_INDEXFIELD
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Dump lattice based profiles
  !>
  !> All quantities defined on the lattice are dumped as function of x, averaged over the entire yz-plane.
  !> The following quantities are dumped:
  !> - fluid velocity
  !> - fluid stress
  !> - particle stress
  !> - particle volume fraction

  subroutine dump_lattice_profiles(N)
    type(lbe_site), dimension(0:, 0:, 0:), intent(in) :: N !< lattice

    ! Declare variables.
    integer :: x_loc, y_loc, z_loc ! local coordinates
    integer :: x_gl ! global coordinate
    integer :: ierror ! MPI error
    integer :: filename_unit ! file index
    real(kind=rk), dimension(1:tnx) :: den_fluid_loc, den_fluid_gl ! local and global fluid density profile
    real(kind=rk), dimension(3, 1:tnx) :: vel_fluid_loc, vel_fluid_gl ! local and global fluid velocity profile
    real(kind=rk), dimension(9, 1:tnx) :: str_fluid_loc, str_fluid_gl ! local and global fluid stress profile
    real(kind=rk), dimension(3, 1:tnx) :: force_part_loc, force_part_gl ! local and global particle force profile
    real(kind=rk), dimension(3, 1:tnx) :: str_part_gl ! global particle stress profile
    real(kind=rk), dimension(1:tnx) :: vol_frac_loc, vol_frac_gl ! local and global particle volume fraction
    character(len=256) filename ! dump file name
    character(len=1024) header ! header line

    ! Check time step and execution condition.
    if((mod(nt, time_step_dump_profiles) .ne. 0) .or. (dump_profiles .eqv. .false.)) return

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("starting subroutine dump_lattice_profiles...", .false.)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Prepare fluid density profile. !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! It is assumed that the velocity \c vel_phys is already up to date.

    ! Reset velocity array.
    den_fluid_loc(:) = 0.d0

    ! Loop over local lattice and compute averages.
    do x_loc = 1, nx
      ! Find global x-position
      x_gl = x_loc + nx * ccoords(1)

      do y_loc = 1, ny
        do z_loc = 1, nz
          den_fluid_loc(x_gl) = den_fluid_loc(x_gl) + density(N(x_loc, y_loc, z_loc))
        end do
      end do
    end do

    ! Reduce fluid velocity profile and store in memory of root.
    call MPI_Reduce(den_fluid_loc, den_fluid_gl, tnx, MPI_REAL8, MPI_SUM, 0, comm_cart, ierror)

    ! Normalize density velocity profile.
    if(myrankc == 0) then
      den_fluid_gl(:) = den_fluid_gl(:) / (tny * tnz)
    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Prepare fluid velocity profile. !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! It is assumed that the velocity \c vel_phys is already up to date.

    ! Reset velocity array.
    vel_fluid_loc(:, :) = 0.d0

    ! Loop over local lattice and compute averages.
    do x_loc = 1, nx
      ! Find global x-position
      x_gl = x_loc + nx * ccoords(1)

      do y_loc = 1, ny
        do z_loc = 1, nz
          vel_fluid_loc(1:3, x_gl) = vel_fluid_loc(1:3, x_gl) + vel_phys(1:3, x_loc, y_loc, z_loc)
        end do
      end do
    end do

    ! Reduce fluid velocity profile and store in memory of root.
    call MPI_Reduce(vel_fluid_loc, vel_fluid_gl, 3 * tnx, MPI_REAL8, MPI_SUM, 0, comm_cart, ierror)

    ! Normalize fluid velocity profile.
    if(myrankc == 0) then
      vel_fluid_gl(:, :) = vel_fluid_gl(:, :) / (tny * tnz)
    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Prepare fluid stress profile. !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! The fluid stresses are not computed on a regular basis.
    ! They have to be computed first.

    ! Compute fluid stresses.
    call calculate_fluid_stresses(str_fluid_loc, N)

    ! Reduce fluid stress profile and store in memory of root.
    call MPI_Reduce(str_fluid_loc, str_fluid_gl, 9 * tnx, MPI_REAL8, MPI_SUM, 0, comm_cart, ierror)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Prepare particle stress profile. !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! The particle stresses are not computed on a regular basis.
    ! They have to be computed first.

    ! Reset particle force array.
    force_part_loc(:, :) = 0.d0

    ! Loop over local lattice and compute averages.
    do x_loc = 1, nx
      ! Find global x-position
      x_gl = x_loc + nx * ccoords(1)

      do y_loc = 1, ny
        do z_loc = 1, nz
          force_part_loc(1:3, x_gl) = force_part_loc(1:3, x_gl) + force_IBM(1:3, x_loc, y_loc, z_loc)
        end do
      end do
    end do

    ! Reduce particle force profile and store in memory of root.
    call MPI_Reduce(force_part_loc, force_part_gl, 3 * tnx, MPI_REAL8, MPI_SUM, 0, comm_cart, ierror)

    ! Compute particle stress profile from particle force profile
    if(myrankc == 0) then
      ! Compute particle stress in the first bin, assuming that there is no force beyond the first bin.
      str_part_gl(1:3, 1) = force_part_gl(1:3, 1) / (2.d0 * tny * tnz)

      ! Compute particle stress in the remaining bins.
      do x_gl = 2, tnx
        str_part_gl(1:3, x_gl) = str_part_gl(1:3, x_gl - 1) &
            & + (force_part_gl(1:3, x_gl) + force_part_gl(1:3, x_gl - 1)) / (2.d0 * tny * tnz)
      end do
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Prepare particle volume fraction profile. !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! The particle volume fraction is not computed on a regular basis.
    ! It has to be computed first.

    ! Compute particle volume fraction.
    call calculate_particle_volume_fraction(vol_frac_loc)

    ! Reduce particle volume fraction profile and store in memory of root.
    call MPI_Reduce(vol_frac_loc, vol_frac_gl, tnx, MPI_REAL8, MPI_SUM, 0, comm_cart, ierror)

    ! Normalize particle volume fraction profile.
    if(myrankc == 0) then
      vol_frac_gl = vol_frac_gl / (tny * tnz)
    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Dump data to the disk. !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Root is responsible for dumping the data.
    if(myrankc == 0) then
      ! Prepare filename.
      filename_unit = 42
      call lbe_make_filename_output(filename, 'lattice_profile', '.asc', nt)

      ! Open file.
      open(unit=filename_unit, file=filename, status='REPLACE', action='WRITE')

      ! Write header lines.
      write(header, "('#lattice profile averaged over yz-planes at time ', i0)") nt
      write(unit=filename_unit, fmt='(A)') trim(header)
      write(header, "('bin    den_fl      vel_x        vel_y        vel_z        ' // &
      & 'str_f_xx     str_f_xy     str_f_xz     str_f_yx     str_f_yy     str_f_yz     ' // &
      & 'str_f_zx     str_f_zy     str_f_zz     str_p_xx     str_p_xy     str_p_xz     density')")
      write(unit=filename_unit, fmt='(A)') trim(header)

      ! Write data, one line for one bin.
      do x_gl = 1, tnx
        write(unit=filename_unit, fmt='(I5.5, X, 17(ES12.5, X))') &
          & x_gl, &
          & den_fluid_gl(x_gl), &
          & vel_fluid_gl(1:3, x_gl), &
          & str_fluid_gl(1:9, x_gl), &
          & str_part_gl(1:3, x_gl), &
          & vol_frac_gl(x_gl)
      end do

      ! Close file.
      close(filename_unit)
    endif

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("finished subroutine dump_lattice_profiles", .false.)
#endif
  end subroutine dump_lattice_profiles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Dump static particle data to the disk
  !>
  !> The static particle data is all data that does not concern the nodes and faces, e.g.,
  !> - center position and velocity
  !> - total force and torque acting on the particle
  !> - deformation energies
  !> These data are
  !> 1) computed for each individual particle by the hosting rank,
  !> 2) gathered via MPI communication, and
  !> 3) written to the disk by rank 0.

  subroutine dump_static_particle_data()
    ! Declare variables.
    integer :: c_i ! particle index
    integer :: i ! counter
    integer :: ierror ! MPI error code
    integer :: stat   ! allocate status
    integer :: filename_unit ! file index
    integer, dimension(1:nprocs) :: num_part_in_rank ! number of particles in each rank
    integer, dimension(1:nprocs) :: displacements ! displacements for memory access
    character(len=256) :: filename ! dump file name
    character(len=1024) :: header ! header line
    type(dump_static_particle_info), dimension(:), allocatable :: static_particle_data_gl ! global static particle data

    ! Check time step and execution condition.
    if((mod(nt, time_step_dump_particles) .ne. 0) .or. (dump_particles .eqv. .false.)) return

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("starting subroutine dump_static_particle_data...", .false.)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Prepare data for dumping. !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Compute all required observables which are not in memory yet.
    do c_i = 1, num_particles_loc
      call compute_particle_stresslet(particles(part_ind(c_i)))
      call compute_inertia_tensor(particles(part_ind(c_i)))
      call compute_inertia_ellipsoid(particles(part_ind(c_i)))
      call compute_angular_velocity(particles(part_ind(c_i)))
    end do

    ! Copy particle data to send buffer.
    do c_i = 1, num_particles_loc
      static_particle_data_loc(c_i)%particle_index_gl = particles(part_ind(c_i))%particle_index_gl
      static_particle_data_loc(c_i)%volume = particles(part_ind(c_i))%volume
      static_particle_data_loc(c_i)%surface = particles(part_ind(c_i))%surface
      static_particle_data_loc(c_i)%center(1) = particles(part_ind(c_i))%center(1) &
       & + real(particles(part_ind(c_i))%num_jumps(1) * tnx, kind=rk)
      static_particle_data_loc(c_i)%center(2) = particles(part_ind(c_i))%center(2) &
       & + real(particles(part_ind(c_i))%num_jumps(2) * tny, kind=rk)
      static_particle_data_loc(c_i)%center(3) = particles(part_ind(c_i))%center(3) &
       & + real(particles(part_ind(c_i))%num_jumps(3) * tnz, kind=rk)
      static_particle_data_loc(c_i)%linear_velocity(:) = particles(part_ind(c_i))%linear_velocity(:)
      static_particle_data_loc(c_i)%linear_momentum(:) = particles(part_ind(c_i))%linear_momentum(:)
      static_particle_data_loc(c_i)%angular_velocity(:) = particles(part_ind(c_i))%angular_velocity(:)
      static_particle_data_loc(c_i)%angular_momentum(:) = particles(part_ind(c_i))%angular_momentum(:)
      static_particle_data_loc(c_i)%force_total(:) = particles(part_ind(c_i))%force_total(:)
      static_particle_data_loc(c_i)%torque_total(:) = particles(part_ind(c_i))%torque_total(:)
      static_particle_data_loc(c_i)%erg_tot = particles(part_ind(c_i))%erg_tot
      static_particle_data_loc(c_i)%erg_s = particles(part_ind(c_i))%erg_s
      static_particle_data_loc(c_i)%erg_b = particles(part_ind(c_i))%erg_b
      static_particle_data_loc(c_i)%erg_v = particles(part_ind(c_i))%erg_v
      static_particle_data_loc(c_i)%erg_at = particles(part_ind(c_i))%erg_at
      static_particle_data_loc(c_i)%erg_int = particles(part_ind(c_i))%erg_int
      static_particle_data_loc(c_i)%stresslet_xx = particles(part_ind(c_i))%stresslet_xx
      static_particle_data_loc(c_i)%stresslet_xy = particles(part_ind(c_i))%stresslet_xy
      static_particle_data_loc(c_i)%stresslet_xz = particles(part_ind(c_i))%stresslet_xz
      static_particle_data_loc(c_i)%stresslet_yy = particles(part_ind(c_i))%stresslet_yy
      static_particle_data_loc(c_i)%stresslet_yz = particles(part_ind(c_i))%stresslet_yz
      static_particle_data_loc(c_i)%stresslet_zz = particles(part_ind(c_i))%stresslet_zz
      static_particle_data_loc(c_i)%iner_ell_axes(:) = particles(part_ind(c_i))%inertia_ell_axes(:)
      static_particle_data_loc(c_i)%iner_ell_vecs(:, :) = particles(part_ind(c_i))%inertia_ell_orient(:, :)
#ifdef IBM_FIXED
      static_particle_data_loc(c_i)%force_anchor(:) = particles(part_ind(c_i))%force_anchor(:)
#endif
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Collect data via MPI. !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Allocate memory for the global static particle data.
    ! If MPI_Allgather is used, the memory has to be allocated by each rank,
    ! otherwise it should only be allocated by root.
#ifdef MPI_ALLGV_FASTER_THAN_GV
    allocate(static_particle_data_gl(num_particles_gl),stat=stat)
    call check_allocate(stat&
         &,'dump_static_particle_data(): static_particle_data_gl')
#else
    if(myrankc == 0) then
      allocate(static_particle_data_gl(num_particles_gl),stat=stat)
      call check_allocate(stat&
           &,'dump_static_particle_data(): static_particle_data_gl')
    endif
#endif

    ! Count number of particles in each rank and tell the root/the other ranks about it.
#ifdef MPI_ALLGV_FASTER_THAN_GV
    call MPI_Allgather(num_particles_loc, 1, MPI_INTEGER, num_part_in_rank, 1, MPI_INTEGER, comm_cart, ierror)
#else
    call MPI_Gather(num_particles_loc, 1, MPI_INTEGER, num_part_in_rank, 1, MPI_INTEGER, 0, comm_cart, ierror)
#endif

    ! Compute displacements.
    displacements(1) = 0

    do i = 2, nprocs
      displacements(i) = displacements(i - 1) + num_part_in_rank(i - 1)
    end do

    ! Gather static particle data.
    ! Use MPI_Allgatherv if MPI_ALLGV_FASTER_THAN_GV is defined, MPI_Gatherv otherwise.
    ! The motivation is that the same approach for MPI_Allreduce and MPI_Reduce increases efficiency on Jugene.
    ! TODO: It has to be verified whether this also holds for MPI_Allgatherv and MPI_Gatherv.
#ifdef MPI_ALLGV_FASTER_THAN_GV
    call MPI_Allgatherv(static_particle_data_loc, num_particles_loc, MPI_STATIC_PARTICLE_DATA, &
            & static_particle_data_gl, num_part_in_rank, displacements, MPI_STATIC_PARTICLE_DATA, comm_cart, ierror)
#else
    call MPI_Gatherv(static_particle_data_loc, num_particles_loc, MPI_STATIC_PARTICLE_DATA, &
            & static_particle_data_gl, num_part_in_rank, displacements, MPI_STATIC_PARTICLE_DATA, 0, comm_cart, ierror)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Dump data to the disk. !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Root is responsible for dumping the data.
    if(myrankc == 0) then
      ! Prepare filename.
      filename_unit = 42
      call lbe_make_filename_output(filename, 'particle_data', '.asc', nt)

      ! Open file.
      open(unit=filename_unit, file=filename, status='REPLACE', action='WRITE')

      ! Write header lines.
      write(header, "('#particle data at time ', i0)") nt
      write(unit=filename_unit, fmt='(A)') trim(header)
#ifndef IBM_FIXED
      write(header, "('part   volume       surface      pos_x           pos_y           pos_z           &
      &vel_x        vel_y        vel_z        lmom_x       lmom_y       lmom_z       omega_x      omega_y      omega_z      &
      &amom_x       amom_y       amom_z       force_x      force_y      force_z      torque_x     torque_y     torque_z     &
      &erg_tot      erg_s        erg_b        erg_v        erg_at       erg_int      strlet_xx    strlet_xy    strlet_xz    &
      &strlet_yy    strlet_yz    strlet_zz    axis_a       axis_b       axis_c       evec_ax      evec_ay      evec_az      &
      &evec_bx      evec_by      evec_bz      evec_cx      evec_cy      evec_cz')")
#else
      write(header, "('part   volume       surface      pos_x           pos_y           pos_z           &
      &vel_x        vel_y        vel_z        lmom_x       lmom_y       lmom_z       omega_x      omega_y      omega_z      &
      &amom_x       amom_y       amom_z       force_x      force_y      force_z      torque_x     torque_y     torque_z     &
      &erg_tot      erg_s        erg_b        erg_v        erg_at       erg_int      strlet_xx    strlet_xy    strlet_xz    &
      &strlet_yy    strlet_yz    strlet_zz    axis_a       axis_b       axis_c       evec_ax      evec_ay      evec_az      &
      &evec_bx      evec_by      evec_bz      evec_cx      evec_cy      evec_cz      f_anch_x     f_anch_y     f_anch_z')")
#endif
      write(unit=filename_unit, fmt='(A)') trim(header)

      ! Write data, one line for one particle.
      do c_i = 1, num_particles_gl
        write(unit=filename_unit, &
#ifndef IBM_FIXED
        & fmt='((I5.5, X), 2(ES12.5, X), 3(ES15.8, X), 18(ES12.5, X), 6(ES12.5, X), 6(ES12.5, X), 12(ES12.5, X))') &
#else
				& fmt='((I5.5, X), 2(ES12.5, X), 3(ES15.8, X), 18(ES12.5, X), 6(ES12.5, X), 6(ES12.5, X), 15(ES12.5, X))') &
#endif
          & static_particle_data_gl(c_i)%particle_index_gl, &
          & static_particle_data_gl(c_i)%volume, &
          & static_particle_data_gl(c_i)%surface, &
          & static_particle_data_gl(c_i)%center(:), &
          & static_particle_data_gl(c_i)%linear_velocity(:), &
          & static_particle_data_gl(c_i)%linear_momentum(:), &
          & static_particle_data_gl(c_i)%angular_velocity(:), &
          & static_particle_data_gl(c_i)%angular_momentum(:), &
          & static_particle_data_gl(c_i)%force_total(:), &
          & static_particle_data_gl(c_i)%torque_total(:), &
          & static_particle_data_gl(c_i)%erg_tot, &
          & static_particle_data_gl(c_i)%erg_s, &
          & static_particle_data_gl(c_i)%erg_b, &
          & static_particle_data_gl(c_i)%erg_v, &
          & static_particle_data_gl(c_i)%erg_at, &
          & static_particle_data_gl(c_i)%erg_int, &
          & static_particle_data_gl(c_i)%stresslet_xx, &
          & static_particle_data_gl(c_i)%stresslet_xy, &
          & static_particle_data_gl(c_i)%stresslet_xz, &
          & static_particle_data_gl(c_i)%stresslet_yy, &
          & static_particle_data_gl(c_i)%stresslet_yz, &
          & static_particle_data_gl(c_i)%stresslet_zz, &
          & static_particle_data_gl(c_i)%iner_ell_axes(:), &
          & static_particle_data_gl(c_i)%iner_ell_vecs(1, :), &
          & static_particle_data_gl(c_i)%iner_ell_vecs(2, :), &
#ifndef IBM_FIXED
          & static_particle_data_gl(c_i)%iner_ell_vecs(3, :)
#else
          & static_particle_data_gl(c_i)%iner_ell_vecs(3, :), &
					& static_particle_data_gl(c_i)%force_anchor
#endif

      end do

      ! Close file.
      close(filename_unit)
    endif

#ifdef MPI_ALLGV_FASTER_THAN_GV
    deallocate (static_particle_data_gl)
#else
    if (myrankc==0) deallocate (static_particle_data_gl)
#endif

    ! Debugging routine
#ifdef IBM_DEBUG
    call MPI_Barrier(comm_cart, ierror)
    call log_msg("finished subroutine dump_static_particle_data", .false.)
#endif
  end subroutine dump_static_particle_data

! endif IBM_PART
#endif

end module lsuperstruct_dump_module
