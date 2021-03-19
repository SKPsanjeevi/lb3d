!> Hoshen-Kopelman auxilliary module
!>
!> This module contains all subroutines required for cluster identification on a lattice.
!> Started by Timm Krüger, May 2012

#include "lbe.h"

module hoshen_kopelman_module

  use lbe_globals_module, only : myrankc, rk
  use lbe_log_module
  use lbe_parallel_module
  use lbe_parms_module, only: dbg_report_hk_timing

  implicit none
  include 'mpif.h'
  private

  public hoshen_kopelman

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Hoshen-Kopelman algorithm
  !>
  !> This is the main Hoshen-Kopelman (HK) algorithm.
  !> For general information, refer to the original algorithm:
  !> - J. Hoshen and R. Kopelman, "Percolation and cluster distribution.
  !>   I. Cluster multiple labeling technique and critical concentration algorithm",
  !>   Phys. Rev. B 14 (8), 1976
  !> - http://www.ocf.berkeley.edu/~fricke/projects/hoshenkopelman/hoshenkopelman.html
  !> The HK algorithm operates on a real scalar input array and returns the clusters in form of a real scalar output array.
  !> The user can specify a threshold region:
  !> - thr_top: scalar values below this threshold are considered in the cluster search.
  !> - thr_bot: scalar values above this threshold are considered in the cluster search.
  !> Example 1:
  !>   A region of scalar values between 0.3 and 0.4 shall be defined as the initial cluster,
  !>   i.e., all sites with scalar values <0.3 and >0.4 shall be sorted into clusters.
  !>   In this case, 'thr_top = 0.3' and 'thr_bot = 0.4'.
  !> Example 2:
  !>   A region of scalar values between 10 and 20 shall be sorted into clusters,
  !>   i.e., all sites with scalar values <10 and >20 shall be defined as the initial cluster.
  !>   In this case, 'thr_bot = 10' and 'thr_top = 20'.
  !> Example 3:
  !>   All scalar values above 0 shall be sorted into clusters.
  !>   In this case, 'thr_bot = 0'. 'thr_top' has to be set sufficiently large.
  !> Example 4:
  !>   All scalar values below 3.14159 shall be sorted into clusters.
  !>   In this case, 'thr_top = 3.14159'. 'thr_bot' has to be set sufficiently small.
  !> Convention for the output array:
  !> - All sites belonging to the initial cluster have value 0.
  !> - All sites subject to cluster identification have values 0 < n <= N
  !>   where n is the cluster index and N is the maximum number of clusters possible.
  !>   Note that N is not necessarily the number of clusters identified.
  !> - During cluster search, non-identified sites have value -1.
  !>   These values should have vanished after successfully completing the algorithm.

  subroutine hoshen_kopelman(input_array, output_array, thr_bot, thr_top, parallelize, n_clusters)
    implicit none

    real(kind=rk), dimension(1:,1:,1:), intent(in) :: input_array !< scalar input array (will not be overwritten)
    integer, dimension(1:,1:,1:), intent(out) :: output_array !< cluster output array (will be overwritten)
    real(kind=rk), intent(in) :: thr_bot, thr_top !< scalar thresholds for defining the initial cluster
    logical, intent(in), optional :: parallelize !< whether to sync the resulting cluster indices with neighboring processes or not
    integer, intent(out), optional :: n_clusters !< return the number of clusters found
    real(kind=rk) :: tstart !< timer

    ! Declare variables.
    integer, dimension(:), allocatable :: indices, indices_compact ! list of cluster indices
    integer :: x, y, z ! lattice coordinates
    integer :: dim_x, dim_y, dim_z ! array dimensions
    integer :: largest_index ! largest cluster label found so far
    integer :: n_x, n_y, n_z ! neighbor indices
    integer :: num_cluster_max ! maximum number of clusters allowed
    integer :: n ! counter
    integer :: ierror

    ! Start timing.
    tstart = MPI_Wtime()
    
    ! Find array dimensions.
    dim_x = size(input_array, 1)
    dim_y = size(input_array, 2)
    dim_z = size(input_array, 3)

    ! Check whether input and output arrays have the same size.
    ! Return if they do not match.
    if((size(output_array, 1) .ne. dim_x) .or. (size(output_array, 2) .ne. dim_y) .or. (size(output_array, 3) .ne. dim_z)) then
      call log_msg("WARNING (hoshen_kopelman): Input and output arrays have different size, aborting HK algorithm...")
      return
    end if

    ! Allocate memory for cluster indices.
    ! In the worst case, there are as many clusters as lattice sites to check.
    ! This will probably never happen, but it is better to be prepared...
    num_cluster_max = dim_x * dim_y * dim_z
    allocate(indices(1:num_cluster_max), stat=ierror)
    call check_allocate(ierror, "Unable to allocate indices buffer for HK.")

    ! Fill the index array with initial values.
    do n = 1, num_cluster_max
      indices(n) = n
    end do

    ! Set initial values of the output array.
    ! All sites belonging to the user-defined cluster are set to 0.
    ! All other sites have to be identified during the HK algorithm and are initially set to -1.
    do x = 1, dim_x
      do y = 1, dim_y
        do z = 1, dim_z
          if(thr_bot < thr_top) then
            if((input_array(x, y, z) < thr_bot) .or. (input_array(x, y, z) > thr_top)) then
              output_array(x, y, z) = 0
            else
              output_array(x, y, z) = -1
            end if
          else
            if((input_array(x, y, z) < thr_top) .or. (input_array(x, y, z) > thr_bot)) then
              output_array(x, y, z) = -1
            else
              output_array(x, y, z) = 0
            end if
          endif
        end do
      end do
    end do

    ! Perform the actual HK algorithm.
    ! The largest cluster index found so far is 0.
    largest_index = 0

    do x = 1, dim_x
      do y = 1, dim_y
        do z = 1, dim_z
          ! Only do anything if the site has to be processed.
          if(output_array(x, y, z) .eq. -1) then
            ! Get indices of neighbor sites which have already been processed.
            ! Take into account that the array is bounded.
            if(x > 1) then
              n_x = output_array(x - 1, y, z)
            else
              n_x = 0
            end if

            if(y > 1) then
              n_y = output_array(x, y - 1, z)
            else
              n_y = 0
            end if

            if(z > 1) then
              n_z = output_array(x, y, z - 1)
            else
              n_z = 0
            end if

            ! If all neighboring sites belong to the initial cluster,
            ! the current site receives a new cluster index.
            if((n_x .eq. 0) .and. (n_y .eq. 0) .and. (n_z .eq. 0)) then
              largest_index = largest_index + 1
              output_array(x, y, z) = largest_index
            ! If exactly one neighboring site belongs to a non-initial cluster,
            ! the current site receives the same index.
            else if((n_x .eq. 0) .and. (n_y .eq. 0) .and. (n_z .gt. 0)) then
              output_array(x, y, z) = find(n_z, indices)
            else if((n_x .eq. 0) .and. (n_y .gt. 0) .and. (n_z .eq. 0)) then
              output_array(x, y, z) = find(n_y, indices)
            else if((n_x .gt. 0) .and. (n_y .eq. 0) .and. (n_z .eq. 0)) then
              output_array(x, y, z) = find(n_x, indices)
            ! If exactly two neighboring sites belong to non-initial clusters,
            ! their indices have to be united and the current site receives their union index.
            else if((n_x .eq. 0) .and. (n_y .gt. 0) .and. (n_z .gt. 0)) then
              call union(n_y, n_z, indices)
              output_array(x, y, z) = find(n_y, indices)
            else if((n_x .gt. 0) .and. (n_y .eq. 0) .and. (n_z .gt. 0)) then
              call union(n_x, n_z, indices)
              output_array(x, y, z) = find(n_x, indices)
            else if((n_x .gt. 0) .and. (n_y .gt. 0) .and. (n_z .eq. 0)) then
              call union(n_x, n_y, indices)
              output_array(x, y, z) = find(n_x, indices)
            ! If all neighboring sites belong to non-initial clusters,
            ! their indices have to be unified and the current site receives their index.
            else
              call union(n_y, n_z, indices)
              call union(n_x, n_y, indices)
              output_array(x, y, z) = find(n_x, indices)
            end if
          end if
        end do
      end do
    end do

    allocate(indices_compact(1:num_cluster_max), stat=ierror)
    call check_allocate(ierror, "Unable to allocate compact indices buffer for HK.")
    indices_compact(:) = -1

    call compactify_mapping(indices, indices_compact, largest_index, largest_index)

    ! Perform a final sweep reducing all cluster indices to the merged and compactified index.
    do x = 1, dim_x
      do y = 1, dim_y
        do z = 1, dim_z
          if ( output_array(x, y, z) .gt. 0 ) then
            output_array(x, y, z) = indices_compact(output_array(x, y, z))
          end if
        end do
      end do
    end do

    if (dbg_report_hk_timing) then
       write(msgstr,"('HK: Serial HK took ',F16.10,' seconds.')") MPI_Wtime()-tstart
       call log_msg(msgstr)
    end if

    ! Sync the resulting cluster indices with neighboring processes.
    if (present(parallelize) ) then
      if (parallelize) then
        if (present(n_clusters) ) then
          call hoshen_kopelman_parallel(output_array, indices_compact, largest_index, n_clusters)
        else
          call hoshen_kopelman_parallel(output_array, indices_compact, largest_index)
        end if
      end if
    end if

    deallocate( indices_compact, stat=ierror )
    call check_allocate(ierror, "Failed deallocation of indices_compact")
    deallocate( indices, stat=ierror )
    call check_allocate(ierror, "Failed deallocation of indices")

    if (dbg_report_hk_timing) then
      write(msgstr,"('HK: HK took ',F16.10,' seconds and found ',I0,' clusters.')") MPI_Wtime()-tstart, n_clusters
      call log_msg(msgstr)
    end if

  end subroutine hoshen_kopelman

  subroutine hoshen_kopelman_parallel(index_array, indices, max_index_l, n_clusters)
    implicit none
    integer, dimension(1:,1:,1:), intent(inout) :: index_array !< cluster output array (will be overwritten)
    integer, dimension(:), intent(in) :: indices
    integer, intent(in) :: max_index_l
    integer, intent(out), optional :: n_clusters

    integer :: x, y, z ! lattice coordinates
    integer :: dim_x, dim_y, dim_z ! array dimensions
    integer, parameter :: tag = 1
    integer, dimension(:,:), allocatable :: sendbuf, recvbuf
    integer, dimension(:), allocatable :: local_id, remote_id, local_map, global_map, merge_map, global_map_compact, num_map_g, displs, sendcnts, gathered_local, gathered_remote
    integer :: i, dir
    integer :: mpierror, status, ierror
    integer :: local, remote, local_enc, remote_enc
    integer :: max_index_g
    integer :: num_map_l, cur_map, next_map, total_num_map, lr_id_len
    integer :: gl, gr
    integer :: mergel, merger, mergemin, merged, nextmap, curmap
    integer :: full_num_map

    logical :: insert_success

    real(kind=rk) :: tstart !< timer

    if (dbg_report_hk_timing) then
       ! Start timing.
       tstart = MPI_Wtime()
       write(msgstr,"('HK: Starting parallel HK: ',F16.10,' seconds.')") MPI_Wtime()-tstart
       call log_msg(msgstr)
     end if

    ! Find dimensions of the provided array index_array.
    dim_x = size(index_array, 1)
    dim_y = size(index_array, 2)
    dim_z = size(index_array, 3)

    ! Get largest global cluster ID.
    call MPI_Allreduce(max_index_l, max_index_g, 1, MPI_INTEGER, MPI_MAX, comm_cart, mpierror)
    ! write(msgstr,"('HK: Globally largest raw cluster index: ',I0)") max_index_g
    ! call log_msg(msgstr)

    ! Edge case: if there are *no* clusters, don't do any of this.
    if ( max_index_g .eq. 0 ) then
      if ( present(n_clusters) ) then
        n_clusters = 0
      end if
      return ! Nothing else to be done - exit subroutine.
    end if

    if (dbg_report_hk_timing) then
       write(msgstr,"('HK: Initializing lists: ',F16.10,' seconds.')") MPI_Wtime()-tstart
       call log_msg(msgstr)
    end if

    ! Allocate lists
    call init_id_lists(local_id, remote_id, dim_x, dim_y, dim_z, max_index_g, lr_id_len, ierror)
    call check_allocate(ierror, "Unable to allocate local_id or remote_id buffer for HK.")
    ! write(msgstr,"('Length of ID lists: ',I0)") lr_id_len
    ! call log_msg(msgstr)

    ! Start comparing the overlapping values, in all three directions, and both up and down.
    ! All values including the 'side-halos' are considered.
    ! "Positive" and "negative" directions are assigned by which way the send buffer will go.
    ! E.g. sendbuf(y,z) = index_array(1,y,z) is taking x = 1, which will be sent downward -> negative.
    ! Unique MPI tags are assigned by 2*dir for negative and 2*dir + 1 for positive directions.

    ! Compare in X-direction
    ! call MPI_Barrier(comm_cart,mpierror)
    ! call log_msg("HK: Starting communication: x-axis")
    allocate( sendbuf(1:dim_y,1:dim_z), stat=ierror )
    call check_allocate(ierror, "Unable to allocate sendbuf x buffer for HK.")
    allocate( recvbuf(1:dim_y,1:dim_z), stat=ierror )
    call check_allocate(ierror, "Unable to allocate recvbuf x buffer for HK.")
    dir = 1

    ! NEGATIVE X
    ! Send buffer will contain the haloed sites
    do y = 1, dim_y
      do z = 1, dim_z
        sendbuf(y,z) = index_array(1,y,z)
      end do
    end do

    call MPI_Sendrecv(sendbuf, dim_y*dim_z, MPI_INTEGER, nnprocs(dir,1), 2*dir, &
         recvbuf, dim_y*dim_z, MPI_INTEGER, nnprocs(dir,2), 2*dir, &
         comm_cart, MPI_STATUS_IGNORE, mpierror)

    ! Compare the recv buffer (which is the halo of the neighbour) to the physical region at max
    do y = 1, dim_y
      do z = 1, dim_z
        local = index_array(dim_x-1,y,z)
        remote = recvbuf(y,z)
        if ( ( local .eq. 0 .and. remote .ne. 0 ) .or. &
             ( local .ne. 0 .and. remote .eq. 0 ) ) then
          write(msgstr, "('Unmatched zero in HK x- : locally ',I0,', but received ', I0, ' from ',I0,' at y = ',I0,', z = ',I0 )") local, remote, nnprocs(dir,2), y, z
          call error(msgstr)
        end if
        if ( local .ne. 0 ) then
          local_enc = encode_cluster(local, myrankc, max_index_g )
          remote_enc = encode_cluster(remote, nnprocs(dir,2), max_index_g )
          call insert_pair_if_new(local_id, remote_id, local_enc, remote_enc, lr_id_len, insert_success )
          if ( .not. insert_success ) call error("Failed to insert mapping pair for HK x- ; array out of bounds.")
        end if
      end do
    end do

    ! POSITIVE X
    ! Send buffer will contain the haloed sites
    do y = 1, dim_y
      do z = 1, dim_z
        sendbuf(y,z) = index_array(dim_x,y,z)
      end do
    end do

    call MPI_Sendrecv(sendbuf, dim_y*dim_z, MPI_INTEGER, nnprocs(dir,2), 2*dir+1, &
         recvbuf, dim_y*dim_z, MPI_INTEGER, nnprocs(dir,1), 2*dir+1, &
         comm_cart, MPI_STATUS_IGNORE, mpierror)

    ! Compare the recv buffer (which is the halo of the neighbour) to the physical region at max
    do y = 1, dim_y
      do z = 1, dim_z
        local = index_array(2,y,z)
        remote = recvbuf(y,z)
        if ( ( local .eq. 0 .and. remote .ne. 0 ) .or. &
             ( local .ne. 0 .and. remote .eq. 0 ) ) then
          write(msgstr, "('Unmatched zero in HK x+ : locally ',I0,', but received ', I0, ' from ',I0,' at y = ',I0,', z = ',I0 )") local, remote, nnprocs(dir,1), y, z
          call error(msgstr)
        end if
        if ( local .ne. 0 ) then
          local_enc = encode_cluster(local, myrankc, max_index_g )
          remote_enc = encode_cluster(remote, nnprocs(dir,1), max_index_g )
          call insert_pair_if_new(local_id, remote_id, local_enc, remote_enc, lr_id_len, insert_success )
          if ( .not. insert_success ) call error("Failed to insert mapping pair for HK x+ ; array out of bounds.")
        end if
      end do
    end do

    deallocate(sendbuf, stat=ierror)
    call check_allocate(ierror, "Failed deallocation of sendbuf x")
    deallocate(recvbuf, stat=ierror)
    call check_allocate(ierror, "Failed deallocation of recvbuf x")

    ! Compare in Y-direction
    ! call MPI_Barrier(comm_cart,mpierror)
    ! call log_msg("HK: Starting communication: y-axis")
    allocate( sendbuf(1:dim_x,1:dim_z), stat=ierror )
    call check_allocate(ierror, "Unable to allocate sendbuf y buffer for HK.")
    allocate( recvbuf(1:dim_x,1:dim_z), stat=ierror )
    call check_allocate(ierror, "Unable to allocate recvbuf y buffer for HK.")
    dir = 2

    ! NEGATIVE Y
    ! Send buffer will contain the haloed sites
    do x = 1, dim_x
      do z = 1, dim_z
        sendbuf(x,z) = index_array(x,1,z)
      end do
    end do

    call MPI_Sendrecv(sendbuf, dim_x*dim_z, MPI_INTEGER, nnprocs(dir,1), 2*dir, &
         recvbuf, dim_x*dim_z, MPI_INTEGER, nnprocs(dir,2), 2*dir, &
         comm_cart, MPI_STATUS_IGNORE, mpierror)

    ! Compare the recv buffer (which is the halo of the neighbour) to the physical region at max
    do x = 1, dim_x
      do z = 1, dim_z
        local = index_array(x,dim_y-1,z)
        remote = recvbuf(x,z)
        if ( ( local .eq. 0 .and. remote .ne. 0 ) .or. &
             ( local .ne. 0 .and. remote .eq. 0 ) ) then
          write(msgstr, "('Unmatched zero in HK y- : locally ',I0,', but received ', I0, ' from ',I0,' at x = ',I0,', z = ',I0 )") local, remote, nnprocs(dir,2), x, z
          call error(msgstr)
        end if
        if ( local .ne. 0 ) then
          local_enc = encode_cluster(local, myrankc, max_index_g )
          remote_enc = encode_cluster(remote, nnprocs(dir,2), max_index_g )
          call insert_pair_if_new(local_id, remote_id, local_enc, remote_enc, lr_id_len, insert_success )
          if ( .not. insert_success ) call error("Failed to insert mapping pair for HK y- ; array out of bounds.")
        end if
      end do
    end do

    ! POSITIVE Y
    ! Send buffer will contain the haloed sites
    do x = 1, dim_x
      do z = 1, dim_z
        sendbuf(x,z) = index_array(x,dim_y,z)
      end do
    end do

    call MPI_Sendrecv(sendbuf, dim_x*dim_z, MPI_INTEGER, nnprocs(dir,2), 2*dir+1, &
         recvbuf, dim_x*dim_z, MPI_INTEGER, nnprocs(dir,1), 2*dir+1, &
         comm_cart, MPI_STATUS_IGNORE, mpierror)

    ! Compare the recv buffer (which is the halo of the neighbour) to the physical region at max
    do x = 1, dim_x
      do z = 1, dim_z
        local = index_array(x,2,z)
        remote = recvbuf(x,z)
        if ( ( local .eq. 0 .and. remote .ne. 0 ) .or. &
             ( local .ne. 0 .and. remote .eq. 0 ) ) then
          write(msgstr, "('Unmatched zero in HK y+ : locally ',I0,', but received ', I0, ' from ',I0,' at x = ',I0,', z = ',I0 )") local, remote, nnprocs(dir,1), x, z
          call error(msgstr)
        end if
        if ( local .ne. 0 ) then
          local_enc = encode_cluster(local, myrankc, max_index_g )
          remote_enc = encode_cluster(remote, nnprocs(dir,1), max_index_g )
          call insert_pair_if_new(local_id, remote_id, local_enc, remote_enc, lr_id_len, insert_success )
          if ( .not. insert_success ) call error("Failed to insert mapping pair for HK y+ ; array out of bounds.")
        end if
      end do
    end do

    deallocate(sendbuf, stat=ierror)
    call check_allocate(ierror, "Failed deallocation of sendbuf y")
    deallocate(recvbuf, stat=ierror)
    call check_allocate(ierror, "Failed deallocation of recvbuf y")

    ! Compare in Z-direction
    ! call MPI_Barrier(comm_cart,mpierror)
    ! call log_msg("HK: Starting communication: z-axis")
    allocate( sendbuf(1:dim_x,1:dim_y), stat=ierror )
    call check_allocate(ierror, "Unable to allocate sendbuf z buffer for HK.")
    allocate( recvbuf(1:dim_x,1:dim_y), stat=ierror )
    call check_allocate(ierror, "Unable to allocate recvbuf z buffer for HK.")
    dir = 3

    ! NEGATIVE Z
    ! Send buffer will contain the haloed sites
    do x = 1, dim_x
      do y = 1, dim_y
        sendbuf(x,y) = index_array(x,y,1)
      end do
    end do

    call MPI_Sendrecv(sendbuf, dim_x*dim_y, MPI_INTEGER, nnprocs(dir,1), 2*dir, &
         recvbuf, dim_x*dim_y, MPI_INTEGER, nnprocs(dir,2), 2*dir, &
         comm_cart, MPI_STATUS_IGNORE, mpierror)

    ! Compare the recv buffer (which is the halo of the neighbour) to the physical region at max
    do x = 1, dim_x
      do y = 1, dim_y
        local = index_array(x,y,dim_z-1)
        remote = recvbuf(x,y)
        if ( ( local .eq. 0 .and. remote .ne. 0 ) .or. &
             ( local .ne. 0 .and. remote .eq. 0 ) ) then
          write(msgstr, "('Unmatched zero in HK z- : locally ',I0,', but received ', I0, ' from ',I0,' at x = ',I0,', y = ',I0 )") local, remote, nnprocs(dir,2), x, y
          call error(msgstr)
        end if
        if ( local .ne. 0 ) then
          local_enc = encode_cluster(local, myrankc, max_index_g )
          remote_enc = encode_cluster(remote, nnprocs(dir,2), max_index_g )
          call insert_pair_if_new(local_id, remote_id, local_enc, remote_enc, lr_id_len, insert_success)
          if ( .not. insert_success ) call error("Failed to insert mapping pair for HK z- ; array out of bounds.")
        end if
      end do
    end do

    ! POSITIVE Z
    ! Send buffer will contain the haloed sites
    do x = 1, dim_x
      do y = 1, dim_y
        sendbuf(x,y) = index_array(x,y,dim_z)
      end do
    end do

    call MPI_Sendrecv(sendbuf, dim_x*dim_y, MPI_INTEGER, nnprocs(dir,2), 2*dir+1, &
         recvbuf, dim_x*dim_y, MPI_INTEGER, nnprocs(dir,1), 2*dir+1, &
         comm_cart, MPI_STATUS_IGNORE, mpierror)

    ! Compare the recv buffer (which is the halo of the neighbour) to the physical region at max
    do x = 1, dim_x
      do y = 1, dim_y
        local = index_array(x,y,2)
        remote = recvbuf(x,y)
        if ( ( local .eq. 0 .and. remote .ne. 0 ) .or. &
             ( local .ne. 0 .and. remote .eq. 0 ) ) then
          write(msgstr, "('Unmatched zero in HK z+ : locally ',I0,', but received ', I0, ' from ',I0,' at x = ',I0,', y = ',I0 )") local, remote, nnprocs(dir,1), x, y
          call error(msgstr)
        end if
        if ( local .ne. 0 ) then
          local_enc = encode_cluster(local, myrankc, max_index_g )
          remote_enc = encode_cluster(remote, nnprocs(dir,1), max_index_g )
          call insert_pair_if_new(local_id, remote_id, local_enc, remote_enc, lr_id_len, insert_success)
          if ( .not. insert_success ) call error("Failed to insert mapping pair for HK z+ ; array out of bounds.")
        end if
      end do
    end do

    deallocate(sendbuf, stat=ierror)
    call check_allocate(ierror, "Failed deallocation of sendbuf z")
    deallocate(recvbuf, stat=ierror)
    call check_allocate(ierror, "Failed deallocation of recvbuf z")

    if (dbg_report_hk_timing) then
       write(msgstr,"('HK: Finding leftover indices: ',F16.10,' seconds.')") MPI_Wtime()-tstart
       call log_msg(msgstr)
    end if

    ! Finally, loop over the lattice and map every extant cluster id to itself, if it was not present yet.
    do x = 1, dim_x
      do y = 1, dim_y
        do z = 1, dim_z
          if ( index_array(x, y, z) .gt. 0 ) then
            call insert_pair_if_local_new(local_id, remote_id, encode_cluster(index_array(x, y, z), myrankc, max_index_g), encode_cluster(index_array(x, y, z), myrankc, max_index_g), lr_id_len, insert_success )
            if ( .not. insert_success ) call error("Failed to insert mapping pair for whole lattice - array out of bounds.")
          end if
        end do
      end do
    end do

    ! call show_id_lists(local_id, remote_id, lr_id_len)
    ! call MPI_Barrier(comm_cart,mpierror)

    ! Get sizes of lists on the various processors.
    if (dbg_report_hk_timing) then
       write(msgstr,"('HK: Getting list sizes: ',F16.10,' seconds.')") MPI_Wtime()-tstart
       call log_msg(msgstr)
    end if

    num_map_l = get_list_size(local_id, lr_id_len)

    allocate( num_map_g(0:nprocs - 1), stat=ierror )
    call check_allocate(ierror, "Unable to allocate list_sizes buffer for HK.")

    call MPI_Allgather(num_map_l, 1, MPI_INTEGER, num_map_g, 1, MPI_INTEGER, comm_cart, mpierror)

    ! Now that we know the list sizes, gather the two id lists from all processors

    if (dbg_report_hk_timing) then
       write(msgstr,"('HK: Gather lists: ',F16.10,' seconds.')") MPI_Wtime()-tstart
       call log_msg(msgstr)
    end if

    total_num_map = sum(num_map_g)
    ! write(msgstr,"('Using total_num_map = sum(num_map_g) = ',I0,' for allocation of gathered_local and gathered_remote.')") total_num_map
    ! call log_msg(msgstr)
    allocate( gathered_local(1:total_num_map), stat=ierror )
    call check_allocate(ierror, "Unable to allocate gathered_local buffer for HK.")
    allocate( gathered_remote(1:total_num_map), stat=ierror )
    call check_allocate(ierror, "Unable to allocate gathered_remote buffer for HK.")
    allocate( displs(0:nprocs - 1), stat=ierror)
    call check_allocate(ierror, "Unable to allocate displs buffer for HK.")

    ! Calculate displacements to let MPI write the local arrays into the 'gathered' array on rank 0.
    displs(0) = 0
    do i=1, nprocs - 1
      displs(i) = displs(i-1) + num_map_g(i-1)
    end do

    ! Write all the local id mappings into the 'gathered' arrays.
    call MPI_Allgatherv(local_id, num_map_l, MPI_INTEGER, gathered_local, num_map_g, displs, MPI_INTEGER, comm_cart, mpierror)
    call MPI_Allgatherv(remote_id, num_map_l, MPI_INTEGER, gathered_remote, num_map_g, displs, MPI_INTEGER, comm_cart, mpierror)

    ! All the mapping data is now present on the root process
    full_num_map = max_index_g * nprocs

    ! Initialize the global map with negative numbers to see which ones are actually used after assignment
    if (dbg_report_hk_timing) then
       write(msgstr,"('HK: Initialize global map: ',F16.10,' seconds.')") MPI_Wtime()-tstart
       call log_msg(msgstr)
    end if
    ! write(msgstr,"('Using full_num_map = max_index_g * nprocs = ',I0,' * ',I0,' = ',I0,' for allocation of global_map and global_map_compact')") max_index_g, nprocs, full_num_map
    ! call log_msg(msgstr)
    allocate( global_map(1:full_num_map), stat=ierror )
    call check_allocate(ierror, "Unable to allocate global_map buffer for HK.")
    global_map(:) = -1

    allocate( merge_map(1:full_num_map), stat=ierror )
    call check_allocate(ierror, "Unable to allocate merge_map buffer for HK.")
    merge_map(:) = -1

    ! call log_msg("HK: Initialize compact global map")
    allocate( global_map_compact(1:full_num_map), stat=ierror )
    call check_allocate(ierror, "Unable to allocate global_map_compact buffer for HK.")
    global_map_compact(:) = -1

    if ( myrankc == 0 ) then
      ! Rank zero now gets to process all received mappings.
      ! We loop over the whole array of pairs.
      if (dbg_report_hk_timing) then
        write(msgstr,"('HK: Process lists: ',F16.10,' seconds.')") MPI_Wtime()-tstart
        call log_msg(msgstr)
      end if
      do i = 1, total_num_map
        if ( global_map(gathered_local(i) ) .lt. 0 ) then
          if ( global_map(gathered_remote(i) ) .lt. 0 ) then
            ! Neither id has been assigned so far:
            ! Set both to the lowest of the pair
            if ( gathered_local(i) .lt. gathered_remote(i) ) then
              global_map(gathered_local(i)) = gathered_local(i)
              global_map(gathered_remote(i)) = gathered_local(i)
            else
              global_map(gathered_local(i)) = gathered_remote(i)
              global_map(gathered_remote(i)) = gathered_remote(i)
            end if
          else
            ! Only the remote id has been assigned already:
            ! Set the local to the remote.
            global_map(gathered_local(i)) = global_map(gathered_remote(i))
          end if
        else
          if ( global_map(gathered_remote(i) ) .lt. 0 ) then
            ! Only the local id has been assigned already:
            ! Set the remote to the local.
            global_map(gathered_remote(i)) = global_map(gathered_local(i))
          else
            ! Now both values are already in the global map:
            ! These are both > 0 because of the above conditions.
            gl = global_map(gathered_local(i) )
            gr = global_map(gathered_remote(i) )

            if ( gl .ne. gr ) then
              ! If their mappings are already equal, we don't carry
              ! additional information and can skip the loop. Otherwise,
              ! find out what everything should merge into.

              ! Some checks are probably not required as they are guaranteed
              ! by the algorithm, but the gains are expected to be minor.
              ! Most importantly, this implementation gets rid of O(N^2)
              ! behaviour.

              mergel = gl
              do while ( ( merge_map(mergel) .gt. 0 ) )
                mergel = merge_map(mergel)
              end do

              merger = gr
              do while ( ( merge_map(merger) .gt. 0 ) )
                merger = merge_map(merger)
              end do

              if ( mergel .lt. merger ) then
                ! Update the merge map for the current local value;
                ! the rest of the chain already points to the right final value.
                if ( mergel .ne. gl ) then
                  merge_map(gl) = mergel
                end if
                ! Update the merge map for the current remote value;
                ! need to travel down the chain to the head of the list
                ! and as such we might as well update all references on the way...
                if ( mergel .ne. gr ) then
                  merge_map(gr) = mergel
                  curmap = gr
                  do while ( ( merge_map(curmap) .gt. 0 ) )
                    nextmap = merge_map(curmap)
                    merge_map(curmap) = mergel
                    curmap = nextmap
                  end do
                end if
                ! And we can set the values in question directly to the merge
                ! value as well.
                global_map(gathered_local(i)) = mergel
                global_map(gathered_remote(i)) = mergel
              else if ( merger .lt. mergel ) then
                ! Update the merge map for the current local value;
                ! the rest of the chain already points to the right final value.
                if ( merger .ne. gr ) then
                  merge_map(gr) = merger
                end if
                ! Update the merge map for the current remote value;
                ! need to travel down the chain to the head of the list
                ! and as such we might as well update all references on the way...
                if ( merger .ne. gl ) then
                  merge_map(gl) = merger
                  curmap = gl
                  do while ( ( merge_map(curmap) .gt. 0 ) )
                    nextmap = merge_map(curmap)
                    merge_map(curmap) = merger
                    curmap = nextmap
                  end do
                end if
                ! And we can set the values in question directly to the merge
                ! value as well.
                global_map(gathered_local(i)) = merger
                global_map(gathered_remote(i)) = merger
              else
                ! If both are equal, we don't need to update either chain, we
                ! just set both tails if necessary.
                if ( mergel .ne. gl ) then
                  merge_map(gl) = mergel
                end if
                if ( merger .ne. gr ) then
                  merge_map(gr) = merger
                end if
                global_map(gathered_local(i)) = mergel
                global_map(gathered_remote(i)) = merger
              end if
            end if
          end if
        end if
      end do

      !> Apply merge map
      if (dbg_report_hk_timing) then
        write(msgstr,"('HK: Apply merge map: ',F16.10,' seconds.')") MPI_Wtime()-tstart
        call log_msg(msgstr)
      end if

      do i = 1, full_num_map
        merged = global_map(i)
        if ( merged .gt. 0 ) then
          do while ( merge_map(merged) .gt. 0 )
            merged = merge_map(merged)
          end do
        end if
        if ( merged .gt. 0 ) then
          global_map(i) = merged
        end if
      end do

      !> Rescale indices
      if (dbg_report_hk_timing) then
        write(msgstr,"('HK: Rescale indices: ',F16.10,' seconds.')") MPI_Wtime()-tstart
        call log_msg(msgstr)
      end if
      call compactify_mapping( global_map, global_map_compact, full_num_map, next_map )
      ! call log_msg("HK: Report indices")
      ! call show_global_map( global_map, global_map_compact, full_num_map )
    end if

    call MPI_Barrier(comm_cart,mpierror)

    deallocate( gathered_local, stat=ierror )
    call check_allocate(ierror, "Failed deallocation of gathered_local")
    deallocate( gathered_remote, stat=ierror )
    call check_allocate(ierror, "Failed deallocation of gathered_remote")
    deallocate( num_map_g, stat=ierror )
    call check_allocate(ierror, "Failed deallocation of num_map_g")

    !> Prepare to scatter mappings
    if (dbg_report_hk_timing) then
      write(msgstr,"('HK: Scatter mappings: ',F16.10,' seconds.')") MPI_Wtime()-tstart
      call log_msg(msgstr)
    end if
    ! call log_msg("HK: Scatter mappings")
    displs(0) = 0
    do i=1, nprocs - 1
      displs(i) = displs(i-1) + max_index_g
    end do

    allocate( sendcnts(0:nprocs-1), stat=ierror )
    call check_allocate(ierror, "Unable to allocate sendcnts buffer for HK.")
    sendcnts(:) = max_index_g

    allocate( local_map(1:max_index_g), stat=ierror )
    call check_allocate(ierror, "Unable to allocate local_map buffer for HK.")
    call MPI_Scatterv( global_map_compact, sendcnts, displs, MPI_INTEGER, local_map, max_index_g, MPI_INTEGER, 0, comm_cart, mpierror)
    ! call MPI_Barrier(comm_cart,mpierror)

    deallocate( sendcnts, stat=ierror )
    call check_allocate(ierror, "Failed deallocation of sendcnts")
    deallocate( displs, stat=ierror )
    call check_allocate(ierror, "Failed deallocation of displs")

    !> Redo the indexing
    if (dbg_report_hk_timing) then
      write(msgstr,"('HK: Redo lattice sweep: ',F16.10,' seconds.')") MPI_Wtime()-tstart
      call log_msg(msgstr)
    end if
    ! call log_msg("HK: Redo lattice sweep")
    do x = 1, dim_x
      do y = 1, dim_y
        do z = 1, dim_z
          ! Ignore zero cluster, do not remap
          if ( index_array(x, y, z) .gt. 0) then
            ! If the target is negative, do not remap
            if ( local_map( index_array(x, y, z) ) .gt. 0 ) then
              index_array(x, y, z) = local_map( index_array(x, y, z) )
            end if
          end if
        end do
      end do
    end do

    ! Cleanup of parallel arrays
    deallocate( local_map, stat=ierror )
    call check_allocate(ierror, "Failed deallocation of local_map")
    deallocate( global_map_compact, stat=ierror )
    call check_allocate(ierror, "Failed deallocation of global_map_compact")
    deallocate( merge_map, stat=ierror )
    call check_allocate(ierror, "Failed deallocation of merge_map")
    deallocate( global_map, stat=ierror )
    call check_allocate(ierror, "Failed deallocation of global_map")
    deallocate( local_id, stat=ierror )
    call check_allocate(ierror, "Failed deallocation of local_id")
    deallocate( remote_id, stat=ierror )
    call check_allocate(ierror, "Failed deallocation of remote_id")

    !> Broadcast the number of clusters found
    if ( present(n_clusters) ) then
      if (dbg_report_hk_timing) then
        write(msgstr,"('HK: Broadcast clusters: ',F16.10,' seconds.')") MPI_Wtime()-tstart
        call log_msg(msgstr)
      end if
      n_clusters = next_map
      call MPI_Bcast( n_clusters, 1, MPI_INTEGER, 0, comm_cart, mpierror)
    end if

    ! call MPI_Barrier(comm_cart,mpierror)
    if (dbg_report_hk_timing) then
      write(msgstr,"('HK: Finished parallel HK: ',F16.10,' seconds.')") MPI_Wtime()-tstart
      call log_msg(msgstr)
    end if

  end subroutine hoshen_kopelman_parallel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Find cluster index
  !>
  !> This is a helper function for the Hoshen-Kopelman algorithm.
  !> It finds the lowest cluster index of united clusters.

  pure function find(n, indices)
    integer :: find !< return value
    integer, intent(in) :: n !< index to find
    integer, dimension(1:), intent(in) :: indices !< list of indices

    find = n

    do
      if(indices(find) .eq. find) exit
      find = indices(find)
    end do
  end function find

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unite cluster indices
  !>
  !> This is a helper function for the Hoshen-Kopelman algorithm.
  !> It links two cluster indices.

  pure subroutine union(index_1, index_2, indices)
    integer, intent(in) :: index_1, index_2 ! indices to unify
    integer, dimension(1:), intent(inout) :: indices ! cluster index list

    ! Declare variables.
    integer :: a, b

    a = find(index_1, indices)
    b = find(index_2, indices)

    if(a > b) then
      indices(a) = b
    else
      indices(b) = a
    end if
  end subroutine union

  !> Encode a cluster ID based on rank and largest index.
  pure integer function encode_cluster(id, rank, max_index)
    implicit none
    integer, intent(in) :: id, rank, max_index
    encode_cluster = id + rank * max_index
    ! write(msgstr,"(4(I0,x))") id, rank, max_index, encode_cluster
    ! call log_msg(msgstr)
  end function encode_cluster

  !> Decode a cluster ID based on rank and largest index.
  pure subroutine decode_cluster(id, max_index, id_local, rank)
    implicit none
    integer, intent(in) :: id, max_index
    integer, intent(out) :: id_local, rank

    rank = 0
    id_local = id
    do while ( id_local > max_index ) 
      id_local = id_local - max_index
      rank = rank + 1
    end do
  end subroutine decode_cluster

  !> Allocate and reset both ID lists to contain only -1.
  pure subroutine init_id_lists(local_list, remote_list, nx, ny, nz, max_index, arrlen, ierror)
    implicit none
    integer, dimension(:), intent(inout), allocatable :: local_list, remote_list
    integer, intent(in)  :: nx, ny, nz, max_index
    integer, intent(out) :: arrlen, ierror

    integer :: i

    ! Array size upper bounds:
    ! * If all sites on the domain edges contain different cluster ids (not actually
    !   possible because of the needed 'checkerboard' pattern) the maximum number of
    !   mappings is the number of lattice sites (first argument).
    ! * There can be at most mappings from max_index to max_index, but then also w.r.t.
    !   to all neighbouring domains. In 3D, this means a prefactor of 3^3 - 1 = 26.
    ! The space needed is then the smallest upper bound.
    arrlen = min( 2*( (nx*ny) + (ny*nz) + (nx*nz) ), 26 * max_index * max_index )

    allocate( local_list(1:arrlen), stat=ierror )
    if ( ierror .ne. 0 ) return
    allocate( remote_list(1:arrlen), stat=ierror )
    local_list(:) = -1
    remote_list(:) = -1
  end subroutine init_id_lists

  !> Insert a new mapping pair into the lists if the pair doesn't exist yet.
  !> Success is defined as either finding the pair exists and returning, or inserting the pair.
  !> Failure is then going out of bounds of the array.
  pure subroutine insert_pair_if_new(local_list, remote_list, local, remote, max_list_size, success)
    implicit none

    integer, dimension(:), intent(inout) :: local_list, remote_list
    integer, intent(in)  :: local, remote, max_list_size
    logical, intent(out) :: success

    integer :: i

    success = .true.
    i = 1
    do while ( local_list(i) .ge. 0 )
      if ( i .gt. max_list_size ) then
        success = .false.
        return
      end if
      if ( local_list(i) .eq. local .and. remote_list(i) .eq. remote ) return
      i = i +1
    end do
    ! write(msgstr,"('Inserting at ',I0,' : ',2(I0,X))") i, local, remote
    ! call log_msg(msgstr,.true.)

    local_list(i) = local
    remote_list(i) = remote
  end subroutine insert_pair_if_new

  !> Insert a new mapping pair into the lists if the local ID doesn't exist yet.
  !> Success is defined as either finding the local ID exists and returning, or inserting the pair.
  !> Failure is then going out of bounds of the array.
  pure subroutine insert_pair_if_local_new(local_list, remote_list, local, remote, max_list_size, success)
    implicit none

    integer, dimension(:), intent(inout) :: local_list, remote_list
    integer, intent(in)  :: local, remote, max_list_size
    logical, intent(out) :: success

    integer :: i

    success = .true.
    i = 1
    do while ( local_list(i) .ge. 0 )
      if ( i .gt. max_list_size ) then
        success = .false.
        return
      end if
      if ( local_list(i) .eq. local ) return
      i = i +1
    end do
    ! write(msgstr,"('Inserting local new at ',I0,' : ',2(I0,X))") i, local, remote
    ! call log_msg(msgstr,.true.)

    local_list(i) = local
    remote_list(i) = remote
  end subroutine insert_pair_if_local_new

  !> Find number of actual relevant elements in a cluster index map.
  pure integer function get_list_size(local_list, max_list_size)
    implicit none

    integer, dimension(:), intent(in) :: local_list
    integer, intent(in) :: max_list_size
    integer :: i

    i = 1
    do while ( ( local_list(i) .ge. 0 ) .and. ( i .le. max_list_size ) )
      i = i + 1
    end do
    get_list_size = i - 1
  end function get_list_size

  !> Will remove all ID gaps from an array to map cluster indices to eachother.
  subroutine compactify_mapping( map, map_compact, max_index_in, max_index_out )
    implicit none

    integer, dimension(:), intent(in) :: map            !< Unprocessed mapping data
    integer, dimension(:), intent(inout) :: map_compact !< Mapping data without any gaps
    integer, intent(in) :: max_index_in                 !< Number of elements to loop over
    integer, intent(out) :: max_index_out               !< Maximum index after removing gaps
    integer, dimension(:), allocatable  :: map_remap

    integer :: i, ierror
    integer :: next_map

    allocate( map_remap(1:max_index_in), stat=ierror )
    call check_allocate(ierror, "Unable to allocate map_remap buffer for HK.")
    map_remap(:) = -1

    ! It is assumed that previous steps of the algorithm supply an array where
    ! the first occurence of any id happends before the first occurence of any larger id
    ! (this is indeed the case!). We now fill an array mapping uncompactified ids
    ! to their compactified versions...
    next_map = 0
    do i = 1, max_index_in
      if ( map_remap(map(i)) .eq. -1 ) then
        next_map = next_map + 1
        map_remap(map(i)) = next_map
      end if
    end do
    ! ... and fill the compactified array with the mapped ids.
    do i = 1, max_index_in
      map_compact(i) = map_remap(map(i))
    end do

    max_index_out = next_map

    deallocate(map_remap, stat=ierror)
    call check_allocate(ierror, "Failed deallocation of map_remap")
  end subroutine compactify_mapping

  subroutine show_id_lists(local_list, remote_list, max_list_size)
    implicit none
    integer, dimension(:), intent(inout) :: local_list, remote_list
    integer, intent(in) :: max_list_size
    integer :: i

    i = 1
    do while ( ( local_list(i) .ge. 0 ) .and. ( i .le. max_list_size ) )
      write(msgstr,"(I4, ' local = ',I12,' , remote = ',I0)") i, local_list(i), remote_list(i)
      call log_msg(msgstr,.true.)
      i = i +1
    end do
  end subroutine show_id_lists

  subroutine show_global_map(global_map, global_map_compact, length)
    implicit none
    integer, dimension(:), intent(inout) :: global_map, global_map_compact
    integer, intent(in) :: length
    integer :: i

    do i = 1, length
      write(msgstr,"('global_map(',I12,') = ',I12, ' global_map_compact(',I12,') = ',I12)") i, global_map(i), i, global_map_compact(i)
      call log_msg(msgstr)
    end do
  end subroutine show_global_map

end module hoshen_kopelman_module

