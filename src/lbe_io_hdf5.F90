#include "lbe.h"

!> hdf5 IO
module lbe_io_hdf5_module
#ifdef USEHDF

  use HDF5
  use lbe_globals_module
  use lbe_helper_module, only: rotate_coordinates
  use lbe_log_module
  use lbe_parallel_module
  use lbe_parms_module, only: arg_input_dfile, arg_input_dfile_set, inp_file, &
       lbeversion, lbeplatform, flag_list, &
       dump_double, hdf_use_ibm_largeblock_io, hdf_use_independent_io,nx,ny,nz,nt, &
       dbg_report_hdf5, dbg_report_hdf5_timing
  use lbe_io_helper_module, only: lbe_make_filename_output

  implicit none
  include 'mpif.h'

  integer, parameter :: hdf5_metadata_len = 80
  character(len=hdf5_metadata_len), allocatable :: hdf5_metadata(:)
  character(len=32) :: hdfversion

  private

  public lbe_io_init_hdf5, lbe_init_metadata_hdf5, lbe_io_shutdown_hdf5
  public read_iscalar_phdf5, read_scalar_phdf5, read_vector_phdf5
  public dump_iscalar_phdf5, dump_scalar_phdf5, dump_vector_phdf5
  public lbe_read_cp_param_phdf5, lbe_write_cp_param_phdf5
  public lbe_add_inputfile_metadata_hdf5

contains

!>This routine calls \c h5open_f() to initialise HDF5. It should
!>only be called once during a given program run.
subroutine lbe_io_init_hdf5()
  implicit none
  integer :: err
  integer :: majnum, minnum, relnum

  call log_msg_hdr("Initializing HDF5")
  call h5open_f(err)
  call h5get_libversion_f(majnum, minnum, relnum, err)
  write(hdfversion,"(I0,'.',I0,'.',I0)") majnum, minnum, relnum
  write(msgstr,"('Using HDF5 version ',I0,'.',I0,'.',I0,' .')") majnum, minnum, relnum
  call log_msg(msgstr)

end subroutine lbe_io_init_hdf5

!>This routine calls \c h5close_f() to initialise HDF5. It should
!>only be called once during a given program run.
subroutine lbe_io_shutdown_hdf5()
  implicit none
  integer :: err

  call log_msg_ws("Shutting down HDF5...")
  call h5close_f(err)
end subroutine lbe_io_shutdown_hdf5

!> Reads an HDF5 file containig an integer scalar field in parallel.
!>
!> All CPUs open the file, but each CPU only loads its own subsection
!> of the file. The whole size of the file is supposed to be \c
!> tnx*tny*tnz data points.
!>
!> \param[in,out] iscalar The local chunk of the integer scalar field
!> is stored here, so the array ranges must be \c (1:nx,1:ny,1:nz) .
!>
!> \param[in] filename complete path and name of the HDF file to read from
subroutine read_iscalar_phdf5(iscalar,filename)
  implicit none

  integer, intent(inout) :: iscalar(1:,1:,1:)
  character(len=*), intent(in) :: filename
  character(len=8), parameter :: dsetname= 'OutArray' !Dataset name
  integer, parameter :: ndim = 3 ! Dataset rank

  integer :: info                ! MPI

  ! HDF variables
  integer(hid_t) :: file_id      ! File identifier
  integer(hid_t) :: plist_id     ! Property list identifier
  integer(hid_t) :: dset_id      ! Dataset identifier
  integer(hid_t) :: dsp_id       ! Dataspace identifier
  integer(hid_t) :: mdsp_id      ! Dataspace identifier
  integer        :: err          ! Capture HDF errors

  integer(hsize_t), dimension(3) :: stride = (/1,1,1/)    ! Stride
  integer(hsize_t), dimension(3) :: count = (/1,1,1/)     ! Number of blocks read from file
  integer(hsize_t), dimension(3) :: blocksize             ! Size of a block
  integer(hsize_t), dimension(3) :: offset                ! Offset of the block read from file
  integer(hsize_t), dimension(3) :: offset_m = (/0,0,0/)  ! Offset of the block in memory

  integer, dimension(nx,ny,nz) :: tmp_iscalar ! temporary storage
  integer(hsize_t), dimension(3) :: dims      ! dims of tmp_iscalar

  info = MPI_INFO_NULL ! MPI
  if (hdf_use_ibm_largeblock_io) then
    if (dbg_report_hdf5) call log_msg("HDF using IBM_largeblock_io")
    call MPI_Info_create(info, err)
    call MPI_Info_set(info, "IBM_largeblock_io", "true", err)
  end if

  ! Since we just want to take a block of data from the file and use
  ! it straightaway, these are all equal:
  dims = (/nx,ny,nz/)
  blocksize = (/nx,ny,nz/)

  ! Calculate which block to take
  offset = ccoords*dims

  if (dbg_report_hdf5) call log_msg_hdr("HDF scalar integer read debug started")

  write(msgstr,"('HDF attempting to read from file <',A,'>')") trim(filename)
  if (dbg_report_hdf5) call log_msg(msgstr)

  if (dbg_report_hdf5) call log_msg("HDF opening file")
  call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, err)
  call h5pset_fapl_mpio_f(plist_id, Comm_Cart, info, err)
  call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,err,access_prp = plist_id)
  call h5pclose_f(plist_id, err)

  if (dbg_report_hdf5) call log_msg("HDF accessing dataset")
  ! Access the dataset dsetname ('OutArray') from file and select its dataspace
  call h5dopen_f(file_id, dsetname, dset_id, err)
  call h5dget_space_f(dset_id, dsp_id, err)

  if (dbg_report_hdf5) call log_msg("HDF creating hyperslab 1")
  ! Select the data we want from the dataspace
  call h5sselect_hyperslab_f(dsp_id,H5S_SELECT_SET_F,offset,count,err,block=blocksize)

  ! Create a new dataspace to dump the data to - same shape as
  ! iscalar, hence the offset is zero
  if (dbg_report_hdf5) call log_msg("HDF creating hyperslab 2")
  call h5screate_simple_f(ndim, dims, mdsp_id, err)
  call h5sselect_hyperslab_f(mdsp_id,H5S_SELECT_SET_F,offset_m,count,err,block=blocksize)

  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, err)
  if (hdf_use_independent_io) then
    if (dbg_report_hdf5) call log_msg("HDF using H5FD_MPIO_INDEPENDENT_F")
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, err)
  else
    if (dbg_report_hdf5) call log_msg("HDF using H5FD_MPIO_COLLECTIVE_F")
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, err)
  end if

  ! read the data into tmp_iscalar
  if (dbg_report_hdf5) call log_msg("HDF reading dataset")
  call H5dread_f(dset_id,H5T_NATIVE_INTEGER,tmp_iscalar,dims,err&
       &,file_space_id=dsp_id,mem_space_id=mdsp_id,xfer_prp=plist_id)

  ! Close both dataspaces (s), the dataset (d) and finally the file (f)
  if (dbg_report_hdf5) call log_msg("HDF closing")
  call h5pclose_f(plist_id, err)
  call h5sclose_f(dsp_id, err)
  call h5sclose_f(mdsp_id, err)
  call h5dclose_f(dset_id, err)
  call h5fclose_f(file_id, err)
  if (dbg_report_hdf5) call log_msg("HDF finished closing handles")
  ! All done with HDF for now - all that remains is to move the data
  ! from the scalar array to the rock_state in a particular fashion.

  if (dbg_report_hdf5) call log_msg_hdr("HDF scalar integer read debug finished")

  iscalar(1:nx,1:ny,1:nz) = tmp_iscalar
end subroutine read_iscalar_phdf5

!> Reads an HDF5 file into the rock matrix.
!> All CPUs open the file, but each CPU only loads its own
!> subsection of the rock.
subroutine read_scalar_phdf5(scalar, filename, opt_rot)
  implicit none
  real(kind=8), intent(inout) :: scalar(1:,1:,1:) ! This is not using rk because of the HDF type.
  character(len=*), intent(in) :: filename
  character(len=*), intent(in), optional :: opt_rot

  character(len=8), parameter :: dsetname= 'OutArray' !Dataset name
  integer, parameter :: ndim = 3 ! Dataset rank

  integer :: info                ! MPI

  ! HDF variables
  integer(hid_t) :: file_id      ! File identifier
  integer(hid_t) :: plist_id     ! Property list identifier
  integer(hid_t) :: dset_id      ! Dataset identifier
  integer(hid_t) :: dsp_id       ! Dataspace identifier
  integer(hid_t) :: mdsp_id      ! Dataspace identifier
  integer        :: err          ! Capture HDF errors

  integer(hsize_t), dimension(3) :: stride = (/1,1,1/)    ! Stride
  integer(hsize_t), dimension(3) :: count = (/1,1,1/)     ! Number of blocks read from file
  integer(hsize_t), dimension(3) :: blocksize             ! Size of a block
  integer(hsize_t), dimension(3) :: offset                ! Offset of the block read from file
  integer(hsize_t), dimension(3) :: offset_m = (/0,0,0/)  ! Offset of the block in memory

  integer(hsize_t), dimension(3) :: dims                  ! Dims of scalar

  integer rnx, rny, rnz ! Rotated coordinates
  integer, dimension(3) :: rccoords
  character(len=3) :: rot

  ! Default: no rotation ('xyz' -> 'xyz')
  if (present(opt_rot)) then
    rot = opt_rot
  else
    rot = "xyz"
  end if

  call rotate_coordinates(rot, nx, ny, nz, rnx, rny, rnz)
  call rotate_coordinates(rot, ccoords(1), ccoords(2), ccoords(3), rccoords(1), rccoords(2), rccoords(3))

  info = MPI_INFO_NULL ! MPI
  if (hdf_use_ibm_largeblock_io) then
    if (dbg_report_hdf5) call log_msg('HDF using IBM_largeblock_io')
    call MPI_Info_create(info, err)
    call MPI_Info_set(info, "IBM_largeblock_io", "true", err)
  end if

  ! Since we just want to take a block of data from the file and use it straightaway, these are all equal:
  dims = (/rnx,rny,rnz/)
  blocksize = (/rnx,rny,rnz/)

  ! Calculate which block to take
  offset(1) = rccoords(1)*dims(1)
  offset(2) = rccoords(2)*dims(2)
  offset(3) = rccoords(3)*dims(3)

  ! Specify that we'll be using MPIO to open the datafile and open it

  if (dbg_report_hdf5) call log_msg_hdr("HDF scalar real read debug started")

  write(msgstr,"('HDF attempting to read from file <',A,'>')") trim(filename)
  if (dbg_report_hdf5) call log_msg(msgstr)

  if (dbg_report_hdf5) call log_msg("HDF opening file")
  call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, err)
  call h5pset_fapl_mpio_f(plist_id, Comm_Cart, info, err)
  call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, err, access_prp = plist_id)
  call h5pclose_f(plist_id, err)

  if (dbg_report_hdf5) call log_msg("HDF accessing dataset")
  ! Access the dataset dsetname ('OutArray') from file and select its dataspace
  call h5dopen_f(file_id, dsetname, dset_id, err)
  call h5dget_space_f(dset_id, dsp_id, err)

  if (dbg_report_hdf5) call log_msg("HDF creating hyperslab 1")
  ! Select the data we want from the dataspace
  call h5sselect_hyperslab_f(dsp_id, H5S_SELECT_SET_F, offset, count, err, block = blocksize)

  ! Create a new dataspace to dump the data to - same shape as the 'scalar' array, hence the offset is zero
  if (dbg_report_hdf5) call log_msg("HDF creating hyperslab 2")
  call h5screate_simple_f(ndim, dims, mdsp_id, err)
  call h5sselect_hyperslab_f(mdsp_id, H5S_SELECT_SET_F, offset_m, count, err, block = blocksize)

  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, err)
  if (hdf_use_independent_io) then
    if (dbg_report_hdf5) call log_msg("HDF using H5FD_MPIO_INDEPENDENT_F")
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, err)
  else
    if (dbg_report_hdf5) call log_msg("HDF using H5FD_MPIO_COLLECTIVE_F")
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, err)
  end if

  ! Finally read the data into 'scalar', yay!
  if (dbg_report_hdf5) call log_msg("HDF reading dataset")
  call H5dread_f(dset_id, H5T_NATIVE_DOUBLE, scalar, dims, err, file_space_id = dsp_id, mem_space_id = mdsp_id, xfer_prp = plist_id)

  ! Close both dataspaces (s), the dataset (d) and finally the file (f)
  if (dbg_report_hdf5) call log_msg("HDF closing")
  call h5pclose_f(plist_id, err)
  call h5sclose_f(dsp_id, err)
  call h5sclose_f(mdsp_id, err)
  call h5dclose_f(dset_id, err)
  call h5fclose_f(file_id, err)
  if (dbg_report_hdf5) call log_msg("HDF finished closing handles")
  ! All done with HDF for now - all that remains is to move the data from the scalar array to the rock_state in a particular fashion.

  if (dbg_report_hdf5) call log_msg_hdr("HDF scalar real read debug finished")
end subroutine read_scalar_phdf5

!> Reads an HDF5 grid vector file.
!> All CPUs open the files, but each CPU only loads its own
!> subsection of the grid.
subroutine read_vector_phdf5(vector,filename,vdim)
  implicit none
  real(kind=8), dimension(:,:,:,:), intent(inout) :: vector ! This is not using rk because of the HDF type.
  integer, intent(in) :: vdim
  character(len=*), intent(in) :: filename

  character(len=8), parameter :: dsetname= 'OutArray' !Dataset name
  integer, parameter :: ndim = 4 ! Dataset rank

  integer :: info                ! MPI

  ! HDF variables
  integer(hid_t) :: file_id      ! File identifier
  integer(hid_t) :: plist_id     ! Property list identifier
  integer(hid_t) :: dset_id      ! Dataset identifier
  integer(hid_t) :: dsp_id       ! Dataspace identifier in file
  integer(hid_t) :: mdsp_id      ! Dataspace identifier in memory
  integer        :: err          ! Capture HDF errors

  integer(hsize_t), dimension(4) :: stride = (/1,1,1,1/)    ! Stride
  integer(hsize_t), dimension(4) :: count = (/1,1,1,1/)     ! Number of blocks read from file
  integer(hsize_t), dimension(4) :: blocksize               ! Size of a block
  integer(hsize_t), dimension(4) :: offset                ! Offset of the block read from file
  integer(hsize_t), dimension(4) :: offset_m = (/0,0,0,0/)  ! Offset of the block in memory

  integer(hsize_t), dimension(4) :: dims                  ! Dims of scalar
  integer                        :: i,j,k                 ! Dummy loop variables

  info = MPI_INFO_NULL ! MPI
  if (hdf_use_ibm_largeblock_io) then
    if (dbg_report_hdf5) call log_msg("HDF using IBM_largeblock_io")
    call MPI_Info_create(info, err)
    call MPI_Info_set(info, "IBM_largeblock_io", "true", err)
  end if

#ifndef HDF5_FLIP
  ! Since we just want to take a block of data from the file and use it straightaway, these are all equal:
  dims = (/vdim,nx,ny,nz/)
  blocksize = (/vdim,nx,ny,nz/)

  ! Calculate which block to take
  offset(1) = 0
  offset(2) = ccoords(1)*dims(2)
  offset(3) = ccoords(2)*dims(3)
  offset(4) = ccoords(3)*dims(4)
#else
  ! Since we just want to take a block of data from the file and use it straightaway, these are all equal:
  dims = (/nx,ny,nz,vdim/)
  blocksize = (/nx,ny,nz,vdim/)

  ! Calculate which block to take
  offset(1) = ccoords(1)*dims(1)
  offset(2) = ccoords(2)*dims(2)
  offset(3) = ccoords(3)*dims(3)
  offset(4) = 0
#endif

  ! Specify that we'll be using MPIO to open the datafile and open it

  if (dbg_report_hdf5) call log_msg("HDF vector real read debug started")

  write(msgstr,"('HDF attempting to read from file <',A,'>')") trim(filename)
  if (dbg_report_hdf5) call log_msg(msgstr)

  if (dbg_report_hdf5) call log_msg("HDF opening file")
  call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, err)
  call h5pset_fapl_mpio_f(plist_id, Comm_Cart, info, err)
  call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, err, access_prp = plist_id)
  call h5pclose_f(plist_id, err)

  if (dbg_report_hdf5) call log_msg("HDF accessing dataset")
  ! Access the dataset dsetname ('OutArray') from file and select its dataspace
  call h5dopen_f(file_id, dsetname, dset_id, err)
  call h5dget_space_f(dset_id, dsp_id, err)

  if (dbg_report_hdf5) call log_msg("HDF creating hyperslab 1")
  ! Select the data we want from the dataspace
  call h5sselect_hyperslab_f(dsp_id, H5S_SELECT_SET_F, offset, count, err, block = blocksize)

  ! Create a new dataspace to dump the data to - same shape as the
  ! 'vector' array, hence the offset is zero
  if (dbg_report_hdf5) call log_msg("HDF creating hyperslab 2")
  call h5screate_simple_f(ndim, dims, mdsp_id, err)
  call h5sselect_hyperslab_f(mdsp_id, H5S_SELECT_SET_F, offset_m, count, err, block = blocksize)

  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, err)
  if (hdf_use_independent_io) then
    if (dbg_report_hdf5) call log_msg("HDF using H5FD_MPIO_INDEPENDENT_F")
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, err)
  else
    if (dbg_report_hdf5) call log_msg("HDF using H5FD_MPIO_COLLECTIVE_F")
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, err)
  end if

  ! Finally read the data into 'vector', yay!
  if (dbg_report_hdf5) call log_msg("HDF reading dataset")
  call H5dread_f(dset_id, H5T_NATIVE_DOUBLE, vector, dims, err, file_space_id = dsp_id, mem_space_id = mdsp_id, xfer_prp = plist_id)

  ! Close both dataspaces (s), the dataset (d) and finally the file (f)
  if (dbg_report_hdf5) call log_msg("HDF closing")
  call h5pclose_f(plist_id, err)
  call h5sclose_f(dsp_id, err)
  call h5sclose_f(mdsp_id, err)
  call h5dclose_f(dset_id, err)
  call h5fclose_f(file_id, err)
  if (dbg_report_hdf5) call log_msg("HDF finished closing handles")

  if (dbg_report_hdf5) call log_msg_hdr("HDF vector real read debug finished")
end subroutine read_vector_phdf5

!>  Called by subroutine dump_scalar.
!>  For HDF5 formatted IO. Only available for postprocessing (post=.true.)
subroutine dump_scalar_phdf5(scalar,name)
  implicit none
  real(kind=rk), dimension(1:, 1:, 1:  ),intent(inout) :: scalar

  character(len=*), intent(in) :: name
  character(len=1024) :: filename
  character(len=1024) :: attrname

  character(len=8), parameter :: dsetname = 'OutArray' !Dataset name

  integer(hid_t) :: file_id       ! File identifier
  integer(hid_t) :: dset_id       ! Dataset identifier
  integer(hid_t) :: filespace     ! Dataspace identifier in file
  integer(hid_t) :: memspace      ! Dataspace identifier in memory
  integer(hid_t) :: plist_id      ! Property list identifier
  integer(hid_t) :: type_id       ! Datatype id for array (real or double)

  integer(hsize_t),  dimension(3) :: dimsf
  integer(hsize_t),  dimension(3) :: dims

  integer(hsize_t),  dimension(3) :: count
  integer(hssize_t), dimension(3) :: offset
  ! integer(hsize_t),  dimension(3) :: chunk_dims

  integer, parameter :: ndim = 3 ! Dataset rank

  integer :: err ! Error flags
  integer :: info ! for MPI
  integer :: comm ! for MPI

  ! For debugging purposes, we have a lot of variables storing time, some MPI_Barrier calls and some MPI data gathering
  double precision :: t_s, t_f, t_init_s, t_init_f, t_write_s, t_write_f
  double precision :: t_close_s, t_close_f, t_closef_s, t_closef_f
  double precision :: t_wait, t_init, t_write, t_close, t_closef, t_total_nowait, t_total, p_wait
  integer :: tag = 1
  integer :: i
  integer status(MPI_STATUS_SIZE) 

  if (dbg_report_hdf5_timing) then
    call MPI_Barrier(Comm_cart,err)
    t_s = MPI_Wtime()
  end if

  call lbe_make_filename_output(filename, trim(name), '.h5', nt)
  attrname = filename ! attrname has to be unique, can't have it twice

  if (dbg_report_hdf5) call log_msg_hdr("HDF scalar real write debug started")

  write(msgstr,"('HDF attempting to write to file <',A,'>')") trim(filename)
  if (dbg_report_hdf5) call log_msg(msgstr)

  info = MPI_INFO_NULL

  if (hdf_use_ibm_largeblock_io) then
    if (dbg_report_hdf5) call log_msg("HDF using IBM_largeblock_io")
    call MPI_Info_create(info, err)
    call MPI_Info_set(info, "IBM_largeblock_io", "true", err)
  end if

  dimsf  = (/tnx, tny, tnz/)
  dims   = (/nx, ny, nz/)
  count  = (/nx, ny, nz/)
  ! chunk_dims = (/nx, ny, nz/)

  ! Note: the above is in fact equivalent to the count/stride/block formulation:
  ! count  = (/1, 1, 1/)
  ! stride = (/1, 1, 1/)
  ! block  = (/nx, ny, nz/)
  ! (with matching parameters in h5sselect_hyperslab_f

  offset(1) = ccoords(1)*nx
  offset(2) = ccoords(2)*ny
  offset(3) = ccoords(3)*nz

  if (dump_double) then
    if (dbg_report_hdf5) call log_msg("HDF datatype is H5T_NATIVE_DOUBLE")
    type_id = H5T_NATIVE_DOUBLE
  else
    if (dbg_report_hdf5) call log_msg("HDF datatype is H5T_NATIVE_REAL")
    type_id = H5T_NATIVE_REAL
  end if

  if (dbg_report_hdf5_timing) then
    call MPI_Barrier(Comm_cart,err)
    t_init_s = MPI_Wtime()
  end if

  ! Setup file access property list with parallel I/O access.
  if (dbg_report_hdf5) call log_msg("HDF creating file")
  call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, err)
  call h5pset_fapl_mpio_f(plist_id, Comm_cart, info, err)
  ! Create the file collectively.
  call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, err, access_prp = plist_id)
  if (dbg_report_hdf5) call log_msg("HDF closing property list handle")
  call h5pclose_f(plist_id, err)

  ! Create the data space for the dataset.
  if (dbg_report_hdf5) call log_msg("HDF creating filespace")
  call h5screate_simple_f(ndim, dimsf, filespace, err)

  ! Create chunked dataset.
  ! This should hopefully be needed nevermore
  ! if (dbg_report_hdf5) call log_msg("HDF creating chunked dataset")
  ! call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, err)
  ! call h5pset_chunk_f(plist_id, ndim, chunk_dims, err)
  ! call h5dcreate_f(file_id, dsetname, type_id, filespace, dset_id, err, plist_id)
  ! call h5pclose_f(plist_id, err)

  ! Create continuous dataset.
  if (dbg_report_hdf5) call log_msg("HDF creating continuous dataset")
  call h5dcreate_f(file_id, dsetname, type_id, filespace, dset_id, err)

  call h5sclose_f(filespace, err)

  ! Each process defines dataset in memory and writes it to the hyperslab in the file.
  if (dbg_report_hdf5) call log_msg("HDF creating memspace")
  call h5screate_simple_f(ndim, dims, memspace, err)

  ! Select hyperslab in the file.
  if (dbg_report_hdf5) call log_msg("HDF selecting hyperslab")
  call h5dget_space_f(dset_id, filespace, err)
  call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, err)

  ! Create property list for collective dataset write
  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, err)
  if (hdf_use_independent_io) then
    if (dbg_report_hdf5) call log_msg("HDF using H5FD_MPIO_INDEPENDENT_F")
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, err)
  else
    if (dbg_report_hdf5) call log_msg("HDF using H5FD_MPIO_COLLECTIVE_F")
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, err)
  end if

  if (dbg_report_hdf5_timing) then
    t_init_f = MPI_Wtime()
    call MPI_Barrier(Comm_cart,err)
    t_write_s = MPI_Wtime()
  end if

  ! Different write calls for double or single data (have to convert real*8 scalar to real*4)
  if (dump_double) then
    if (dbg_report_hdf5) call log_msg("HDF writing double data")
    call h5dwrite_f(dset_id, type_id, scalar, dims, err, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
  else
    if (dbg_report_hdf5) call log_msg("HDF writing single data")
    call h5dwrite_f(dset_id, type_id, real(scalar,4), dims, err, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
  end if

  if (dbg_report_hdf5_timing) then
    t_write_f = MPI_Wtime()
    call MPI_Barrier(Comm_cart,err)
    t_close_s = MPI_Wtime()
  end if

  ! Close dataspaces.
  if (dbg_report_hdf5) call log_msg("HDF closing filespace handle")
  call h5sclose_f(filespace, err)
  if (dbg_report_hdf5) call log_msg("HDF closing memspace handle")
  call h5sclose_f(memspace, err)
  ! Close the dataset.
  if (dbg_report_hdf5) call log_msg("HDF closing dataset handle")
  call h5dclose_f(dset_id, err)
  ! Close the property list.
  if (dbg_report_hdf5) call log_msg("HDF closing property list handle")
  call h5pclose_f(plist_id, err)

  if (dbg_report_hdf5_timing) then
    t_close_f = MPI_Wtime()
    call MPI_Barrier(Comm_cart,err)
    t_closef_s = MPI_Wtime()
  end if

  ! Close the file.
  if (dbg_report_hdf5) call log_msg("HDF closing file handle")
  call h5fclose_f(file_id, err)
  if (dbg_report_hdf5) call log_msg("HDF finished closing handles")

  if (dbg_report_hdf5_timing) then
    t_closef_f = MPI_Wtime()
  end if

  ! Call subroutine which adds the metadata, now the raw dataset exists
  ! Only possible from one processor
  if (myrankc .eq. 0 ) then
    if (dbg_report_hdf5) call log_msg("HDF writing metadata")
    call lbe_write_attr_phdf5(filename, dsetname, attrname)
    if (dbg_report_hdf5) call log_msg("HDF finished writing metadata")
  end if

  if (dbg_report_hdf5) call log_msg_hdr("HDF scalar real write debug finished")

  if (dbg_report_hdf5_timing) then
    ! This is a lot of debugging/timing stuff
    ! All processors create a string with their timings, then send it to rank 0.
    ! Rank 0 can then display them all in correct order
    call MPI_Barrier(Comm_cart,err)
    t_f = MPI_Wtime()

    call log_msg_hdr("HDF debug timer output started")
    if (dbg_report_hdf5) call log_msg("  RANK             INIT            WRITE            CLOSE       CLOSE FILE       WORK TOTAL             WAIT            TOTAL           %-WAIT")
    t_wait = t_write_s - t_init_f + t_close_s - t_write_f + t_closef_s - t_close_f + t_f - t_closef_f
    t_init = t_init_f - t_init_s
    t_write = t_write_f - t_write_s
    t_close = t_close_f - t_close_s
    t_closef = t_closef_f - t_closef_s

    t_total_nowait = t_init + t_write + t_close + t_closef
    t_total = t_f - t_s
    p_wait = t_wait / t_total

    ! Magic number '256' corresponds to the length of msgstr of course
    if ( myrankc .gt. 0 ) then
      write(msgstr,'(I6.6,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10)') myrankc, t_init, t_write, t_close, t_closef, t_total_nowait, t_wait, t_total, 100.0*p_wait
      call MPI_Send(msgstr, 256, MPI_CHARACTER, 0, tag, Comm_cart, err)
    else
      write(msgstr,'(I6.6,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10)') myrankc, t_init, t_write, t_close, t_closef, t_total_nowait, t_wait, t_total, 100.0*p_wait
      if (dbg_report_hdf5) call log_msg(msgstr)
      do i=1, nprocs-1
        call MPI_Recv(msgstr, 256, MPI_CHARACTER, i, tag, Comm_cart, status, err)
        if (dbg_report_hdf5) call log_msg(msgstr)
      end do
    end if

    call log_msg_hdr("HDF debug timer output finished")
  end if

end subroutine dump_scalar_phdf5

!> Called by subroutine dump_iscalar. For HDF5 formatted IO. Dumps
!> integer scalar fields.
subroutine dump_iscalar_phdf5(iscalar,name)
  implicit none
  integer,dimension(1:,1:,1:),intent(inout) :: iscalar
  character(len=*),intent(in) :: name

  character(len=1024) :: filename
  character(len=1024) :: attrname

  character(len=8), parameter :: dsetname = 'OutArray' !Dataset name

  integer(hid_t) :: file_id       ! File identifier
  integer(hid_t) :: dset_id       ! Dataset identifier
  integer(hid_t) :: filespace     ! Dataspace identifier in file
  integer(hid_t) :: memspace      ! Dataspace identifier in memory
  integer(hid_t) :: plist_id      ! Property list identifier
  integer(hid_t) :: type_id       ! Datatype id for array (real or double)

  integer(hsize_t),  dimension(3) :: dimsf
  integer(hsize_t),  dimension(3) :: dims

  integer(hsize_t),  dimension(3) :: count
  integer(hssize_t), dimension(3) :: offset
  ! integer(hsize_t),  dimension(3) :: chunk_dims

  integer, parameter :: ndim = 3 ! Dataset rank

  integer :: err ! Error flags
  integer :: info ! for MPI
  integer :: comm ! for MPI

  ! For debugging purposes, we have a lot of variables storing time,
  ! some MPI_Barrier calls and some MPI data gathering
  double precision :: t_s,t_f,t_init_s,t_init_f,t_write_s,t_write_f
  double precision :: t_close_s,t_close_f,t_closef_s,t_closef_f
  double precision :: t_wait,t_init,t_write,t_close,t_closef,t_total_nowait,t_total,p_wait
  integer, parameter :: tag=1
  integer :: i
  integer status(MPI_STATUS_SIZE)

  if ( dbg_report_hdf5_timing ) then
    call MPI_Barrier(Comm_cart,err)
    t_s = MPI_Wtime()
  end if

  call lbe_make_filename_output(filename, trim(name), '.h5', nt)
  attrname = filename ! attrname has to be unique, can't have it twice

  if (dbg_report_hdf5) call log_msg_hdr("HDF scalar integer write debug started")

  write(msgstr,"('HDF attempting to write to file <',A,'>')") trim(filename)
  if (dbg_report_hdf5) call log_msg(msgstr)

  info = MPI_INFO_NULL
  if (hdf_use_ibm_largeblock_io) then
    if (dbg_report_hdf5) call log_msg("HDF using IBM_largeblock_io")
    call MPI_Info_create(info, err)
    call MPI_Info_set(info, "IBM_largeblock_io", "true", err)
  end if

  dimsf  = (/tnx, tny, tnz/)
  dims   = (/nx, ny, nz/)
  count  = (/nx, ny, nz/)
  ! chunk_dims = (/nx, ny, nz/)

  ! Note: the above is in fact equivalent to the count/stride/block formulation:
  ! count  = (/1, 1, 1/)
  ! stride = (/1, 1, 1/)
  ! block  = (/nx, ny, nz/)
  ! (with matching parameters in h5sselect_hyperslab_f

  offset(1) = ccoords(1)*nx
  offset(2) = ccoords(2)*ny
  offset(3) = ccoords(3)*nz

  if (dbg_report_hdf5) call log_msg("HDF datatype is H5T_NATIVE_INTEGER")
  type_id = H5T_NATIVE_INTEGER

  if ( dbg_report_hdf5_timing ) then
    call MPI_Barrier(Comm_cart,err)
    t_init_s = MPI_Wtime()
  end if

  ! Setup file access property list with parallel I/O access.
  if (dbg_report_hdf5) call log_msg("HDF creating file")
  call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, err)
  call h5pset_fapl_mpio_f(plist_id, Comm_cart, info, err)
  ! Create the file collectively.
  call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,err,access_prp=plist_id)
  if (dbg_report_hdf5) call log_msg("HDF closing property list handle")
  call h5pclose_f(plist_id, err)

  ! Create the data space for the dataset.
  if (dbg_report_hdf5) call log_msg("HDF creating filespace")
  call h5screate_simple_f(ndim, dimsf, filespace, err)

  ! Create chunked dataset.
  ! This should hopefully be needed nevermore
  ! if (dbg_report_hdf5) call log_msg("HDF creating chunked dataset")
  ! call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, err)
  ! call h5pset_chunk_f(plist_id, ndim, chunk_dims, err)
  ! call h5dcreate_f(file_id, dsetname, type_id, filespace, dset_id, err, plist_id)
  ! call h5pclose_f(plist_id, err)

  ! Create continuous dataset.
  if (dbg_report_hdf5) call log_msg("HDF creating continuous dataset")
  call h5dcreate_f(file_id, dsetname, type_id, filespace, dset_id, err)

  call h5sclose_f(filespace, err)

  ! Each process defines dataset in memory and writes it to the hyperslab in the file.
  if (dbg_report_hdf5) call log_msg("HDF creating memspace")
  call h5screate_simple_f(ndim, dims, memspace, err)

  ! Select hyperslab in the file.
  if (dbg_report_hdf5) call log_msg("HDF selecting hyperslab")
  call h5dget_space_f(dset_id, filespace, err)
  call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, err)

  ! Create property list for collective dataset write
  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, err)
  if (hdf_use_independent_io) then
    if (dbg_report_hdf5) call log_msg("HDF using H5FD_MPIO_INDEPENDENT_F")
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, err)
  else
    if (dbg_report_hdf5) call log_msg("HDF using H5FD_MPIO_COLLECTIVE_F")
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, err)
  end if

  if ( dbg_report_hdf5_timing ) then
    t_init_f = MPI_Wtime()
    call MPI_Barrier(Comm_cart,err)
    t_write_s = MPI_Wtime()
  end if

  if (dbg_report_hdf5) call log_msg("HDF writing integer data")
  call h5dwrite_f(dset_id,type_id,iscalar,dims,err,file_space_id=filespace&
       &,mem_space_id=memspace,xfer_prp=plist_id)

  if ( dbg_report_hdf5_timing) then
    t_write_f = MPI_Wtime()
    call MPI_Barrier(Comm_cart,err)
    t_close_s = MPI_Wtime()
  end if

  ! Close dataspaces.
  if (dbg_report_hdf5) call log_msg("HDF closing filespace handle")
  call h5sclose_f(filespace, err)
  if (dbg_report_hdf5) call log_msg("HDF closing memspace handle")
  call h5sclose_f(memspace, err)
  ! Close the dataset.
  if (dbg_report_hdf5) call log_msg("HDF closing dataset handle")
  call h5dclose_f(dset_id, err)
  ! Close the property list.
  if (dbg_report_hdf5) call log_msg("HDF closing property list handle")
  call h5pclose_f(plist_id, err)

  if ( dbg_report_hdf5_timing ) then
    t_close_f = MPI_Wtime()
    call MPI_Barrier(Comm_cart,err)
    t_closef_s = MPI_Wtime()
  end if

  ! Close the file.
  if (dbg_report_hdf5) call log_msg("HDF closing file handle")
  call h5fclose_f(file_id, err)
  if (dbg_report_hdf5) call log_msg("HDF finished closing handles")

  if ( dbg_report_hdf5_timing ) then
    t_closef_f = MPI_Wtime()
  end if

  ! Call subroutine which adds the metadata, now the raw dataset exists
  ! Only possible from one processor
  if (myrankc .eq. 0 ) then
    if (dbg_report_hdf5) call log_msg("HDF writing metadata")
    call lbe_write_attr_phdf5(filename, dsetname, attrname)
    if (dbg_report_hdf5) call log_msg("HDF finished writing metadata")
  end if

  if (dbg_report_hdf5) call log_msg_hdr("HDF scalar integer write debug finished")

  if (dbg_report_hdf5_timing) then
    ! This is a lot of debugging/timing stuff
    ! All processors create a string with their timings, then send it to rank 0.
    ! Rank 0 can then display them all in correct order
    call MPI_Barrier(Comm_cart,err)
    t_f = MPI_Wtime()

    call log_msg_hdr("HDF debug timer output started")
    call log_msg("  RANK             INIT            WRITE            CLOSE       CLOSE FILE       WORK TOTAL             WAIT            TOTAL           %-WAIT")

    t_wait = t_write_s - t_init_f + t_close_s - t_write_f + t_closef_s - t_close_f + t_f - t_closef_f
    t_init = t_init_f - t_init_s
    t_write = t_write_f - t_write_s
    t_close = t_close_f - t_close_s
    t_closef = t_closef_f - t_closef_s

    t_total_nowait = t_init + t_write + t_close + t_closef
    t_total = t_f - t_s
    p_wait = t_wait / t_total

    ! Magic number '256' corresponds to the length of msgstr of course
    if ( myrankc .gt. 0 ) then
      write(msgstr,'(I6.6,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10)') myrankc, t_init, t_write, t_close, t_closef, t_total_nowait, t_wait, t_total, 100.0*p_wait
      call MPI_Send(msgstr, 256, MPI_CHARACTER, 0, tag, Comm_cart, err)
    else
      write(msgstr,'(I6.6,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10)') myrankc, t_init, t_write, t_close, t_closef, t_total_nowait, t_wait, t_total, 100.0*p_wait
      call log_msg(msgstr)
      do i=1, nprocs-1
        call MPI_Recv(msgstr, 256, MPI_CHARACTER, i, tag, Comm_cart, status, err)
        call log_msg(msgstr)
      end do
    end if

    call log_msg_hdr("HDF debug timer output finished")
  end if
end subroutine dump_iscalar_phdf5

!> Postprocessed IO done with PHDF5.
!> Dump an array of vectors in hdf5 format.
subroutine dump_vector_phdf5(vector,filename)
  implicit none
  real(kind=rk), dimension(1:,1:,1:,1:),intent(inout) :: vector
  character(len=*), intent(in) :: filename
  character(len=1024) :: attrname
  integer :: vecdim
  integer :: nxi,nyi,nzi
  !HDF5 stuff
  character(len=8), parameter :: dsetname='OutArray' !Dataset name

  integer(hid_t) :: file_id       ! File identifier
  integer(hid_t) :: dset_id       ! Dataset identifier
  integer(hid_t) :: filespace     ! Dataspace identifier in file
  integer(hid_t) :: memspace      ! Dataspace identifier in memory
  integer(hid_t) :: plist_id      ! Property list identifier
  integer(hid_t) :: type_id

  integer(hsize_t), dimension(4) :: dims
  integer(hsize_t), dimension(4) :: dimsf

  integer(hsize_t),  dimension(4) :: count
  integer(hssize_t), dimension(4) :: offset
  ! integer(hsize_t),  dimension(4) :: chunk_dims

  integer , parameter :: ndim = 4 ! Dataset rank

  integer :: err ! Error flags
  integer :: info ! for MPI
  integer :: comm ! for MPI

  ! For debugging purposes, we have a lot of variables storing time, some MPI_Barrier calls and some MPI data gathering
  double precision :: t_s, t_f, t_init_s, t_init_f, t_write_s, t_write_f
  double precision :: t_close_s, t_close_f, t_closef_s, t_closef_f
  double precision :: t_wait, t_init, t_write, t_close, t_closef, t_total_nowait, t_total, p_wait
  integer :: tag = 1
  integer :: i
  integer status(MPI_STATUS_SIZE) 

  if ( dbg_report_hdf5_timing ) then
    call MPI_Barrier(Comm_cart,err)
    t_s = MPI_Wtime()
  end if

#ifndef HDF5_FLIP
  vecdim = size(vector,1)
  nxi    = size(vector,2)
  nyi    = size(vector,3)
  nzi    = size(vector,4)
#else
  nxi    = size(vector,1)
  nyi    = size(vector,2)
  nzi    = size(vector,3)
  vecdim = size(vector,4)
#endif

  attrname = filename

  if (dbg_report_hdf5) call log_msg_hdr("HDF vector real write debug started")

  write(msgstr,"('HDF attempting to write to file <',A,'>')") trim(filename)
  if (dbg_report_hdf5) call log_msg(msgstr)

  write(msgstr,"('nxi = ',I0,', nyi = ',I0,', nzi = ',I0,', vecdim = ',I0)") nxi, nyi, nzi, vecdim
  if (dbg_report_hdf5) call log_msg(msgstr)

  info = MPI_INFO_NULL
  if (hdf_use_ibm_largeblock_io) then
    if (dbg_report_hdf5) call log_msg("HDF using IBM_largeblock_io")
    call MPI_Info_create(info, err)
    call MPI_Info_set(info, "IBM_largeblock_io", "true", err)
  end if

#ifndef HDF5_FLIP
  dimsf  = (/vecdim, tnx, tny, tnz/)
  dims   = (/vecdim, nx, ny, nz/)
  count  = (/vecdim, nx, ny, nz/)
  ! chunk_dims = (/vecdim, nx, ny, nz/)

  offset(1) = 0
  offset(2) = ccoords(1)*nx
  offset(3) = ccoords(2)*ny
  offset(4) = ccoords(3)*nz
#else
  dimsf  = (/tnx, tny, tnz, vecdim/)
  dims   = (/nx, ny, nz, vecdim/)
  count  = (/nx, ny, nz, vecdim/)
  ! chunk_dims = (/nx, ny, nz, vecdim/)

  offset(1) = ccoords(1)*nx
  offset(2) = ccoords(2)*ny
  offset(3) = ccoords(3)*nz
  offset(4) = 0

#endif

  if (dump_double) then
    if (dbg_report_hdf5) call log_msg('HDF datatype is H5T_NATIVE_DOUBLE')
    type_id = H5T_NATIVE_DOUBLE
  else
    if (dbg_report_hdf5) call log_msg('HDF datatype is H5T_NATIVE_REAL')
    type_id = H5T_NATIVE_REAL
  end if

  if ( dbg_report_hdf5_timing ) then
    call MPI_Barrier(Comm_cart,err)
    t_init_s = MPI_Wtime()
  end if

  ! Setup file access property list with parallel I/O access.
  if (dbg_report_hdf5) call log_msg("HDF creating file")
  call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, err)
  call h5pset_fapl_mpio_f(plist_id, Comm_Cart, info, err)
  ! Create the file collectively.
  call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, err, access_prp = plist_id)
  if (dbg_report_hdf5) call log_msg("HDF closing property list handle")
  call h5pclose_f(plist_id, err)

  ! Create the data space for the dataset.
  if (dbg_report_hdf5) call log_msg("HDF creating filespace")
  call h5screate_simple_f(ndim, dimsf, filespace, err)

  ! Create chunked dataset.
  ! This should hopefully be needed nevermore
  ! if (dbg_report_hdf5) call log_msg("HDF creating chunked dataset")
  ! call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, err)
  ! call h5pset_chunk_f(plist_id, ndim, chunk_dims, err)
  ! call h5dcreate_f(file_id, dsetname, type_id, filespace, dset_id, err, plist_id)
  ! call h5pclose_f(plist_id, err)

  ! Create continuous dataset.
  if (dbg_report_hdf5) call log_msg("HDF creating continuous dataset")
  call h5dcreate_f(file_id, dsetname, type_id, filespace, dset_id, err)

  call h5sclose_f(filespace, err)

  ! Each process defines dataset in memory and writes it to the hyperslab in the file.
  if (dbg_report_hdf5) call log_msg("HDF creating memspace")
  call h5screate_simple_f(ndim, dims, memspace, err)  

  ! Select hyperslab in the file.
  if (dbg_report_hdf5) call log_msg("HDF selecting hyperslab")
  call h5dget_space_f(dset_id, filespace, err)
  call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, err)

  ! Create property list for collective dataset write
  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, err)
  if (hdf_use_independent_io) then
    if (dbg_report_hdf5) call log_msg("HDF using H5FD_MPIO_INDEPENDENT_F")
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, err)
  else
    if (dbg_report_hdf5) call log_msg("HDF using H5FD_MPIO_COLLECTIVE_F")
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, err)
  end if

  if ( dbg_report_hdf5_timing ) then
    t_init_f = MPI_Wtime()
    call MPI_Barrier(Comm_cart,err)
    t_write_s = MPI_Wtime()
  end if

  ! Different write calls for double or single data (have to convert real*8 scalar to real*4
  if (dump_double) then
    if (dbg_report_hdf5) call log_msg("HDF writing double data")
    call h5dwrite_f(dset_id, type_id, vector, dims, err, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
  else
    if (dbg_report_hdf5) call log_msg("HDF writing single data")
    call h5dwrite_f(dset_id, type_id, real(vector,4), dims, err, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
  end if

  if ( dbg_report_hdf5_timing ) then
    t_write_f = MPI_Wtime()
    call MPI_Barrier(Comm_cart,err)
    t_close_s = MPI_Wtime()
  end if

  ! Close dataspaces.
  if (dbg_report_hdf5) call log_msg("HDF closing filespace handle")
  call h5sclose_f(filespace, err)
  if (dbg_report_hdf5) call log_msg("HDF closing memspace handle")
  call h5sclose_f(memspace, err)
  ! Close the dataset.
  if (dbg_report_hdf5) call log_msg("HDF closing dataset handle")
  call h5dclose_f(dset_id, err)
  ! Close the property list.
  if (dbg_report_hdf5) call log_msg("HDF closing property list handle")
  call h5pclose_f(plist_id, err)

  if ( dbg_report_hdf5_timing ) then
    t_close_f = MPI_Wtime()
    call MPI_Barrier(Comm_cart,err)
    t_closef_s = MPI_Wtime()
  end if

  ! Close the file.
  if (dbg_report_hdf5) call log_msg("HDF closing file handle")
  call h5fclose_f(file_id, err)
  if (dbg_report_hdf5) call log_msg("HDF finished closing handles")

  if ( dbg_report_hdf5_timing ) then
    t_closef_f = MPI_Wtime()
  end if

  ! Call subroutine which adds the metadata, now the raw dataset exists
  ! Only possible from one processor
  if (myrankc .eq. 0 ) then
    if (dbg_report_hdf5) call log_msg("HDF writing metadata")
    call lbe_write_attr_phdf5(filename, dsetname, attrname)
    if (dbg_report_hdf5) call log_msg("HDF finished writing metadata")
  end if

  if (dbg_report_hdf5) call log_msg_hdr("HDF vector real write debug finished")

  if ( dbg_report_hdf5_timing ) then
    ! This is a lot of debugging/timing stuff
    ! All processors create a string with their timings, then send it to rank 0.
    ! Rank 0 can then display them all in correct order
    call MPI_Barrier(Comm_cart,err)
    t_f = MPI_Wtime()

    call log_msg_hdr("HDF debug timer output started")
    call log_msg("  RANK             INIT            WRITE            CLOSE       CLOSE FILE       WORK TOTAL             WAIT            TOTAL           %-WAIT")

    t_wait = t_write_s - t_init_f + t_close_s - t_write_f + t_closef_s - t_close_f + t_f - t_closef_f
    t_init = t_init_f - t_init_s
    t_write = t_write_f - t_write_s
    t_close = t_close_f - t_close_s
    t_closef = t_closef_f - t_closef_s

    t_total_nowait = t_init + t_write + t_close + t_closef
    t_total = t_f - t_s
    p_wait = t_wait / t_total

    ! Magic number '256' corresponds to the length of msgstr of course
    if ( myrankc .gt. 0 ) then
      write(msgstr,'(I6.6,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10)') myrankc, t_init, t_write, t_close, t_closef, t_total_nowait, t_wait, t_total, 100.0*p_wait
      call MPI_Send(msgstr, 256, MPI_CHARACTER, 0, tag, Comm_cart, err)
    else
      write(msgstr,'(I6.6,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10,X,F16.10)') myrankc, t_init, t_write, t_close, t_closef, t_total_nowait, t_wait, t_total, 100.0*p_wait
      call log_msg(msgstr)
      do i=1, nprocs-1
        call MPI_Recv(msgstr, 256, MPI_CHARACTER, i, tag, Comm_cart, status, err)
        call log_msg(msgstr)
      end do
    end if

    call log_msg_hdr("HDF debug timer output finished")
  end if
end subroutine dump_vector_phdf5

!>This subroutine adds the checkpointing parameters to a file.
!>
!>It is called from root (rank 0) of the checkpointing subroutines.
subroutine lbe_write_cp_param_phdf5(data,filename,aname)
  implicit none

  character(len=*),intent(in) :: filename  ! File name
  character(len=*),intent(in) :: aname     ! Attribute / parameter name

  character(len=8), parameter :: dsetname='OutArray' ! Dataset name

  integer(hid_t)     :: file_id       ! File identifier 
  integer(hid_t)     :: dset_id       ! Dataset identifier 
  integer(hid_t)     :: attr_id       ! Attribute identifier 
  integer(hid_t)     :: aspace_id     ! Attribute Dataspace identifier 
  integer            :: arank = 1                        ! Attribute rank
  integer(hsize_t), dimension(1) :: adims           ! Attribute dimension
  integer            :: error ! Error flag

  ! HDF wants an array, so we make a 1d array of length 1
  real(kind=rk), intent(inout) :: data
  real(kind=8), dimension(1) :: data_array

  adims(1) = 1
  data_array(1) = data

  if (dbg_report_hdf5) call log_msg_hdr("HDF params write debug started")

  write(msgstr,"('HDF attempting to write to file <',A,'>')") trim(filename)
  if (dbg_report_hdf5) call log_msg(msgstr)

  if (dbg_report_hdf5) call log_msg("HDF opening file")
  call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,error)
  if (dbg_report_hdf5) call log_msg("HDF opening dataset")
  call h5dopen_f(file_id,dsetname,dset_id,error)
  ! Create scalar data space for the attribute.
  if (dbg_report_hdf5) call log_msg("HDF creating attribute dataspace")

  write(msgstr,"('  Param length: ',I0)") adims(1)
  if (dbg_report_hdf5) call log_msg(msgstr)

  call h5screate_simple_f(arank, adims, aspace_id, error)
  ! Create datatype for the attribute.
  ! call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
  ! call h5tset_size_f(atype_id, attrlen, error)
  ! Create dataset attribute.
  if (dbg_report_hdf5) call log_msg("HDF creating attribute")
  call h5acreate_f(dset_id, trim(aname), H5T_NATIVE_DOUBLE, aspace_id, attr_id, error)

  ! Write the attribute data.
  write(msgstr,"('HDF writing attribute <',A,'> = ',F16.10)") trim(aname), data_array(1)
  if (dbg_report_hdf5) call log_msg(msgstr)
  call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, data_array, adims, error)

  ! Close the attribute. 
  if (dbg_report_hdf5) call log_msg("HDF closing attribute handle")
  call h5aclose_f(attr_id, error)
  if (dbg_report_hdf5) call log_msg("HDF closing attribute space handle")
  call h5sclose_f(aspace_id, error)
  if (dbg_report_hdf5) call log_msg("HDF closing dataset handle")
  call h5dclose_f(dset_id,error)
  if (dbg_report_hdf5) call log_msg("HDF closing file handle")
  call h5fclose_f(file_id,error)

  if (dbg_report_hdf5) call log_msg_hdr("HDF params write debug finished")

end subroutine lbe_write_cp_param_phdf5

!>This subroutine reads the checkpointing parameters from a file.
!>
!>It is called from root (rank 0) of the checkpointing subroutines.
subroutine lbe_read_cp_param_phdf5(data,filename,aname)
  implicit none

  character(len=*), intent(in) :: filename  ! File name
  character(len=*), intent(in) :: aname  ! Attribute name

  character(len=8), parameter :: dsetname='OutArray' !Dataset name

  integer(hid_t)     :: file_id       ! File identifier 
  integer(hid_t)     :: dset_id       ! Dataset identifier 
  integer(hid_t)     :: attr_id       ! Attribute identifier 
  integer(hid_t)     :: aspace_id     ! Attribute Dataspace identifier 
  integer(hid_t)     :: type_id      ! Attribute Dataspace identifier 
  integer            :: arank = 1                        ! Attribute rank
  integer(hsize_t), dimension(1) :: adims           ! Attribute dimension
  integer            :: error ! Error flag

  real(kind=rk), intent(out) :: data
  real(kind=8), dimension(1) :: data_array

  ! Because HDF5 wants arrays
  adims(1) = 1

  if (dbg_report_hdf5) call log_msg_hdr("HDF params read debug started")

  write(msgstr,"('HDF attempting to read from file <',A,'>')") trim(filename)
  if (dbg_report_hdf5) call log_msg(msgstr)

  if (dbg_report_hdf5) call log_msg("HDF opening file")
  call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,error)
  if (dbg_report_hdf5) call log_msg("HDF opening dataset")
  call h5dopen_f(file_id,dsetname,dset_id,error)

  ! Create datatype for the attribute.
  ! call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
  ! call h5tset_size_f(atype_id, attrlen, error)
  ! Create dataset attribute.
  write(msgstr,"('HDF opening attribute with name <',A,'>')") trim(aname)
  if (dbg_report_hdf5) call log_msg(msgstr)
  call h5aopen_by_name_f(dset_id, ".", trim(aname), attr_id, error)
  if (dbg_report_hdf5) call log_msg("HDF opening attribute dataspace")
  call h5aget_space_f(attr_id, aspace_id, error)

  ! Read the attribute data.
  if (dbg_report_hdf5) call log_msg("HDF reading attribute")
  call h5aread_f(attr_id, H5T_NATIVE_DOUBLE, data_array, adims, error)
  data = data_array(1)

  write(msgstr,"('HDF read attribute <',A,'> = ',F16.10)") trim(aname), data
  if (dbg_report_hdf5) call log_msg(msgstr)

  ! Close the attribute. 
  if (dbg_report_hdf5) call log_msg("HDF closing attribute handle")
  call h5aclose_f(attr_id, error)
  if (dbg_report_hdf5) call log_msg("HDF closing attribute space handle")
  call h5sclose_f(aspace_id, error)
  if (dbg_report_hdf5) call log_msg("HDF closing dataset handle")
  call h5dclose_f(dset_id,error)
  if (dbg_report_hdf5) call log_msg("HDF closing file handle")
  call h5fclose_f(file_id,error)

  if (dbg_report_hdf5) call log_msg_hdr("HDF params read debug finished")

end subroutine lbe_read_cp_param_phdf5

!> This subroutine adds the HDF metadata to each postprocessed file.
!> It is called from root (rank 0) of the subroutines dump_..._phdf5.
subroutine lbe_write_attr_phdf5(filename,dsetname,aname)
  implicit none

  character(len=*), intent(in) :: filename  ! File name
  character(len=*), intent(in) :: dsetname  ! Dataset name
  character(len=*), intent(in) :: aname  ! Attribute name

  integer(hid_t)     :: file_id       ! File identifier 
  integer(hid_t)     :: dset_id       ! Dataset identifier 
  integer(hid_t)     :: attr_id       ! Attribute identifier 
  integer(hid_t)     :: aspace_id     ! Attribute Dataspace identifier 
  integer(hid_t)     :: atype_id      ! Attribute Dataspace identifier 
  integer            :: arank = 1                        ! Attribute rank
  integer(hsize_t), dimension(1) :: adims           ! Attribute dimension
  integer(size_t)    :: attrlen = 80   ! Length of the attribute string
  integer(hsize_t), dimension(7) :: data_dims
  integer            :: error ! Error flag

  adims = size(hdf5_metadata)
  data_dims(1) = size(hdf5_metadata)

  call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,error)
  call h5dopen_f(file_id,dsetname,dset_id,error)
  ! Create scalar data space for the attribute.
  call h5screate_simple_f(arank, adims, aspace_id, error)
  ! Create datatype for the attribute.
  call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
  call h5tset_size_f(atype_id, attrlen, error)
  ! Create dataset attribute.
  call h5acreate_f(dset_id, aname, atype_id, aspace_id, attr_id, error)
  ! Write the attribute data.
  call h5awrite_f(attr_id, atype_id, hdf5_metadata, data_dims, error)
  ! Close the attribute.
  call h5aclose_f(attr_id, error)
  ! Terminate access to the data space.
  call h5sclose_f(aspace_id, error)
  call h5dclose_f(dset_id,error)
  call h5fclose_f(file_id,error)

end subroutine lbe_write_attr_phdf5

!> This subroutine sets up the common metadata and basic input files.
subroutine lbe_init_metadata_hdf5()
  implicit none

  call log_msg_ws("Creating HDF5 metadata ...")
  call lbe_add_common_metadata_hdf5()

  call lbe_add_inputfile_metadata_hdf5(inp_file, "main")
  if ( arg_input_dfile_set ) then
    call lbe_add_inputfile_metadata_hdf5(arg_input_dfile, "diff")
  end if

end subroutine lbe_init_metadata_hdf5

!> This subroutine adds generic metadata, flag information and MPI information
!> to the HDF5 metadata.
subroutine lbe_add_common_metadata_hdf5()
  implicit none
  character(len=32)  :: username
  character(len=32) :: ccdims
  character(len=8)  :: cnprocs
  character(len=128) :: decomposition
  integer :: n_fixed, n_flags, n_mpi
  integer :: i, ierror, offset

  character(len=hdf5_metadata_len), allocatable :: metadata(:)
  character(len=hdf5_metadata_len) :: hdr

  ! Get the username
  call getenv('USER',username)

  n_fixed = 6
  if ( arg_input_dfile_set ) then
    n_fixed = n_fixed + 1
  end if
  n_flags = size(flag_list)
  n_mpi = 6

  allocate(metadata(n_fixed + n_mpi + n_flags + 1),stat=ierror)
  call check_allocate(ierror, "lbe_add_common_metadata_hdf5() : metadata")

  call make_hdf5_metadata_hdr("Metadata", hdr)
  metadata(1) = hdr
  metadata(2) = trim(lbeversion)
  metadata(3) = trim(lbeplatform)
  metadata(4) = '  HDF version: '//trim(hdfversion)
  metadata(5) = 'User: '//username
  metadata(6) = 'Input-file: '//trim(inp_file)
  if (arg_input_dfile_set ) then
    metadata(7) = 'Differential input-file: '//trim(arg_input_dfile)
  end if

  offset = n_fixed

  call make_hdf5_metadata_hdr("Compiler flags used", hdr)
  metadata(offset+1) = hdr
  do i=1,n_flags
    metadata(offset+1+i) = trim(flag_list(i))
  end do

  write(decomposition,"(I0,' x ',I0,' x ', I0)") nx, ny, nz
  write(ccdims,"(I0, ' x ',I0, ' x ', I0)") cdims(1), cdims(2), cdims(3)
  write(cnprocs,'(I0)') nprocs

  offset = offset + n_flags + 1

  call make_hdf5_metadata_hdr("MPI data", hdr)
  metadata(offset+1) = hdr
  metadata(offset+2) = 'Hostname: '//trim(hname)
  metadata(offset+3) = 'Simulation started: '//trim(startsimul)
  metadata(offset+4) = 'Total number of processors: '//trim(cnprocs)
  metadata(offset+5) = 'Processors using a '//trim(ccdims)//' grid.'
  metadata(offset+6) = 'Decomposition: '//trim(decomposition)

  call lbe_append_metadata_hdf5(metadata)

  deallocate(metadata)

end subroutine lbe_add_common_metadata_hdf5

!> Adds the contents of an input file to the HDF5 metadata and puts
!> 'tag' in the header and footer lines.
subroutine lbe_add_inputfile_metadata_hdf5(filename, tag)
  implicit none

  character(len=*), intent(in) :: filename, tag

  integer, parameter :: iunit = 10
  integer, parameter :: n_meta = 4

  character(len=hdf5_metadata_len), allocatable :: metadata(:)
  character(len=hdf5_metadata_len) :: counterchar, hdr
  integer :: n_lines
  integer :: i, eof, ierror

  eof = 0
  n_lines = 0
  open (iunit,file=filename,status='unknown',form='formatted')
  do while (eof .eq. 0)
    read (iunit,'(a80)',iostat = eof) counterchar
    if (eof .lt. 0) exit
    n_lines = n_lines + 1
  end do
  close(iunit)

  allocate(metadata(n_meta + n_lines),stat=ierror)
  call check_allocate(ierror, "lbe_add_inputfile_metadata_hdf5() : metadata")

  metadata(1) = ""
  call make_hdf5_metadata_hdr("Start input file ["//trim(tag)//"]", hdr)
  metadata(2) = hdr

  open(iunit,file=filename,status='unknown',form='formatted')
  do i = 3,n_meta+n_lines-2
    read (iunit,'(a80)') metadata(i)
  end do
  close(iunit)

  call make_hdf5_metadata_hdr("End input file ["//trim(tag)//"]", hdr)
  metadata(n_meta + n_lines - 1) = hdr
  metadata(n_meta + n_lines) = ""

  call lbe_append_metadata_hdf5(metadata)

  deallocate(metadata)

end subroutine lbe_add_inputfile_metadata_hdf5

!> Append an array of strings to the hdf5_metadata array.
subroutine lbe_append_metadata_hdf5(newdata)
  implicit none

  character(len=hdf5_metadata_len), intent(in) :: newdata(:)
  character(len=hdf5_metadata_len), allocatable :: tmp(:)

  integer :: ierror
  
  if ( .not. allocated(hdf5_metadata) ) then
    allocate (hdf5_metadata(0),stat=ierror)
    call check_allocate(ierror,'lbe_append_metadata_hdf5(): hdf5_metadata')
  end if

  allocate(tmp(size(hdf5_metadata)), stat=ierror)
  call check_allocate(ierror, "lbe_append_metadata_hdf5(): tmp")
  tmp(:) = hdf5_metadata(:)
  deallocate(hdf5_metadata)
  
  allocate (hdf5_metadata( size(tmp) + size(newdata) ), stat=ierror)
  call check_allocate(ierror,'lbe_append_metadata_hdf5(): hdf5_metadata')
  if ( size(tmp) .gt. 0) then
    hdf5_metadata(1:size(tmp)) = tmp(:)
  end if
  hdf5_metadata(size(tmp)+1:size(tmp)+size(newdata) ) = newdata(:)

  deallocate(tmp)

end subroutine lbe_append_metadata_hdf5

!> Make a header string of the correct length for HDF5 metadata.
subroutine make_hdf5_metadata_hdr(msg, hdr)
  implicit none

  character(len=*), intent(in)  :: msg
  character(len=hdf5_metadata_len), intent(out) :: hdr

  integer :: len, i

  len = len_trim(msg)
  len = hdf5_metadata_len - len

  hdr = "( " // trim(msg) // " )"

  do i = 1, len
    if ( mod(i,2) == 0 ) then
      hdr = "=" // trim(hdr)
    else
      hdr = trim(hdr) // "="
    end if
  end do

end subroutine make_hdf5_metadata_hdr

#endif
end module lbe_io_hdf5_module
