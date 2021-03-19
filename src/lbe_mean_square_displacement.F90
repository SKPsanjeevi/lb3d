#include "lbe.h"

!> routines to calculate and dump the mean square displacement of
!> discrete particles
!>
!> The interface of this module is pretty concise. \c
!> msd_dump_checkpoint() , \c msd_delete_checkpoint() , \c msd_init()
!> , and \c msd_run() need to be run by the LB3D checkpointing
!> routines, the startup code, and during the main loop. File output
!> is done in \c msd_run() , checkpoint restoration in \c msd_init() .
!>
!> The only additional calls necessary to apply this module for the
!> analysis of a new type of mean square displacement calculations are
!> \c register_msd() to obtain a handle to a new \c
!> mean_square_displacement_type object, \c sample() to sample object
!> positions and \c close_trace() to mark an object as not being
!> supposed to be traced anymore.
!>
!> Reading huge checkpoint files in \c restore_single_checkpoint()
!> might fail due to insuffient memory. In that case it might help to
!> reduce the local buffer size parameter \c n_read_max .
module lbe_mean_square_displacement_mod
    use lbe_analysis_module, only: cached_avg_total_velocity
    use lbe_globals_module, only: halo_extent,myrankc,rk
    use lbe_helper_module, only: every_n_time_steps,is_restoring
    use lbe_io_helper_module, only: lbe_delete_file,lbe_make_filename_cp&
         &,lbe_make_filename_output,lbe_make_filename_restore
#ifdef USEXDRF
    use lbe_io_xdrf_module, only: check_xdrfopen
#endif
    use lbe_log_module
    use lbe_parallel_module, only: calculate_displacements,check_allocate&
         &,checkmpi,comm_cart,nprocs
    use lbe_parms_module, only: nt
    use lbe_types_module, only: lbe_site
    use map_module, only: Mii_clear,Mii_commit,Mii_exist,Mii_init,Mii_map&
         &,Mii_preinsert,Mii_provide_capacity,Mii_type

    implicit none
    private

    public close_trace,msd_dump_checkpoint,msd_delete_checkpoint,msd_init&
         &,msd_run,register_msd,sample

    include 'mpif.h'

    !> accumulators for one data point in a time difference--mean
    !> square displacement curve
    type point_type
       real(kind=rk) :: cum     !< sum of square displacements
       real(kind=rk) :: cum2    !< sum of squared square displacements
       integer :: n             !< number of samples
    end type point_type

    !> one sample of a 1D position of a Lagrangian object
    type sample_type
       integer :: uid           !< object id
       real(kind=rk) :: x       !< 1D object position
    end type sample_type

    !> a trace of the former positions of a Lagrangian object
    type trace_type
       integer :: uid           !< object id
       !> actual number of valid entries at the beginning of \c pos
       integer :: len
       !> history of former object positions, will be allocated with
       !> size \c mean_square_displacement_type%len but only first \c
       !> trace_type%len entries are valid.
       real(kind=rk),allocatable :: pos(:)
    end type trace_type

    !> contains everything needed to calculate the mean square
    !> displacement for one type of objects
    type mean_square_displacement_type
       character(len=32) :: name !< human-readable name for this object type
       integer :: n_dump         !< file output interval (LB steps, 0: never)
       integer :: n_sample       !< sampling interval (LB steps, 0: never)
       integer :: len            !< history length in units of \c n_sample

       !> \c .true. if displacement direction is periodic
       logical :: periodic

       !> periodic dimension of the displacement direction
       real(kind=rk) :: d_periodic

       !> only take into account complete traces of length \c len
       logical :: complete_only

       !> only take into account traces that start between \c
       !> min_startpos and \c max_startpos
       logical :: central_startpos_only

       !> minimum accepted start position for a trace with \c
       !> central_startpos_only
       real(kind=rk) :: min_startpos

       !> maximum accepted start position for a trace with \c
       !> central_startpos_only
       real(kind=rk) :: max_startpos

       !> if set to 1, 2, or 3, adjust sampled positions according to
       !> the x-, y-, or z-component of a possible drift of the whole
       !> system determined by the average system velocity
       integer :: subtract_drift

       !> accumulation buffer for time--mean square displacement
       !> curve. Allocated to size \c len.
       type(point_type),allocatable :: data(:)

       integer :: n_traces      !< current number of object traces
       type(trace_type),allocatable :: traces(:) !< buffer containing traces
       type(Mii_type) :: uid2trace !< maps object ids to index in \c traces

       integer :: n_samples     !< current number of samples
       type(sample_type),allocatable :: samples(:) !< buffer containing samples
    end type mean_square_displacement_type

    !> contains all objects of type \c mean_square_displacement_type
    type(mean_square_displacement_type),allocatable,save :: msds(:)

    !> used to prevent creation of new entries in \c msds after \c
    !> msd_init() was called
    logical,save :: initialized=.false.

    !> custom mpi data type representing \c point_type
    integer,save :: point_mpitype

    !> custom mpi data type representing \c sample_type
    integer,save :: sample_mpitype

    !> custom mpi data type representing \c trace_type except for the
    !> allocatable \c pos array.
    !>
    !> This assumes that the offsets to the elements in \c trace_type
    !> always are the same independently on how large \c pos is which
    !> is probably true if an allocatable element is a pointer,
    !> internally. In case this goes wrong when having \c
    !> mean_square_displacement_type objects with different \c len, \c
    !> mean_square_displacement_type would need its own \c
    !> trace_mpitype.
    integer,save :: trace_mpitype

    !> custom MPI reduction operation to sum up objects of type \c
    !> point_type
    integer,save :: sum_point_mpiop

contains

    !> allocates memory for a single msd object for later usage
    !>
    !> \param[in,out] m msd object
    subroutine allocate_msd(m)
        type(mean_square_displacement_type),intent(inout) :: m
        integer :: stat

        allocate (m%data(m%len),stat=stat)
        call check_allocate(stat,'allocate_msd(): m%data')
        call reset_data(m)

        m%n_samples = 0
        ! this is enlarged afterwards, if necessary
        allocate (m%samples(1),stat=stat)
        call check_allocate(stat,'allocate_msd(): m%samples')

        m%n_traces = 0
        ! this is enlarged afterwards, if necessary
        call allocate_traces(m%traces,1,m%len)
        call reset_traces(m)
        call Mii_init(m%uid2trace)
    end subroutine allocate_msd

    !> allocates an array of traces and the position arrays contained
    !>
    !> \param[in,out] t vector of traces to allocate
    !>
    !> \param[in] n number of traces to create
    !>
    !> \param[in] l length of position arrays
    subroutine allocate_traces(t,n,l)
        type(trace_type),allocatable,intent(inout) :: t(:)
        integer,intent(in) :: n,l
        integer :: i,stat

        allocate (t(n),stat=stat)
        call check_allocate(stat,'allocate_traces(): t')

        do i=1,n
           allocate (t(i)%pos(l),stat=stat)
           call check_allocate(stat,'allocate_traces(): t(i)%pos')
        end do
    end subroutine allocate_traces

    !> append content of recv buffer to the traces
    !>
    !> \param[in,out] m msd object
    !>
    !> \param[in] n number of samples to append
    !>
    !> \param[in,out] buf buffer containing the samples, will be
    !> freed by \c append_samples_to_traces at the end
    !>
    !> \param[in] whole_N local chunk of the lattice with full halo of
    !> depth \c halo_extent
    subroutine append_samples_to_traces(m,n,buf,whole_N)
        type(mean_square_displacement_type),intent(inout) :: m
        integer,intent(in) :: n
        type(sample_type),allocatable,intent(inout) :: buf(:)
        type(lbe_site),intent(in) :: &
             &whole_N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        integer :: i,j,n_new_traces,n_wrap,t
        real(kind=rk) :: drift,dx,new_x,v_avg(3)

        if (any(m%subtract_drift==(/1,2,3/))) then
           v_avg = cached_avg_total_velocity(whole_N)
           ! assume that the drift stayed constant since the last
           ! sample and that possible fluctuations in v_avg are
           ! smaller than a possible drift and smaller than the
           ! displacement recorded here. Anything more elaborate would
           ! require its own checkpointing.
           drift = v_avg(m%subtract_drift)*m%n_sample
        else
           drift = 0.0_rk
        end if

        n_new_traces = 0

        ! it is not possible that more than n traces will be added
        call check_allocate(Mii_provide_capacity(m%uid2trace,m%n_traces+n)&
             &,'append_samples_to_traces(): '&
             &//'m%uid2trace (Mii_provide_capacity())')

        append_buf: do i=1,n
           exist_in_traces: if (Mii_exist(m%uid2trace,buf(i)%uid)) then
              t = Mii_map(m%uid2trace,buf(i)%uid)
              invalid: if (buf(i)%x<0.0_rk) then
                 ! invalidate trace
                 m%traces(t)%len = 0
              else invalid
                 ! shift existing positions and prepend new position
                 ! as the new head of the trace
                 do j=min(m%traces(t)%len,m%len-1),1,-1
                    m%traces(t)%pos(j+1) = m%traces(t)%pos(j)+drift
                 end do

                 new_x = buf(i)%x
                 if (m%periodic.and.m%traces(t)%len>=1) then
                    ! wrap new position until its distance to previous
                    ! position (if existing) is minimal
                    dx = new_x-m%traces(t)%pos(2)
                    n_wrap = nint(dx/m%d_periodic)
                    new_x = new_x - real(n_wrap,kind=rk)*m%d_periodic
                 end if
                 m%traces(t)%pos(1) = new_x
                 if (m%traces(t)%len<m%len) m%traces(t)%len = m%traces(t)%len+1
              end if invalid
           else exist_in_traces
              ! if trace is or became invalid then don't do anything,
              ! otherwise create new trace
              valid: if (buf(i)%x>=0.0_rk) then
                 if (m%n_traces==size(m%traces)) call boost_traces(m)

                 m%n_traces = m%n_traces+1
                 n_new_traces = n_new_traces+1

                 m%traces(m%n_traces)%uid = buf(i)%uid
                 m%traces(m%n_traces)%pos(1) = buf(i)%x
                 m%traces(m%n_traces)%len = 1
              end if valid
           end if exist_in_traces
        end do append_buf

        ! update uid2trace
        do i=m%n_traces+1-n_new_traces,m%n_traces
           call Mii_preinsert(m%uid2trace,m%traces(i)%uid,i)
        end do
        call Mii_commit(m%uid2trace)

        deallocate (buf)
    end subroutine append_samples_to_traces

    !> increases the capacity of the sample buffer of an msd object
    !>
    !> \param[in,out] m msd object
    subroutine boost_samples(m)
        type(mean_square_displacement_type),intent(inout) :: m
        type(sample_type),allocatable :: tmp(:)
        integer :: stat

        allocate (tmp(size(m%samples)),stat=stat)
        call check_allocate(stat,'boost_samples(): tmp')

        tmp = m%samples
        deallocate (m%samples)

        allocate (m%samples(2*size(tmp)),stat=stat)
        call check_allocate(stat,'boost_samples(): m%samples')

        m%samples(1:size(tmp)) = tmp
        deallocate (tmp)
    end subroutine boost_samples

    !> increases the capacity of the trace buffer of an msd object
    !>
    !> \param[in,out] m msd object
    subroutine boost_traces(m)
        type(mean_square_displacement_type),intent(inout) :: m
        type(trace_type),allocatable :: tmp(:)
        integer :: stat

        call allocate_traces(tmp,size(m%traces),m%len)
        tmp = m%traces
        call deallocate_traces(m%traces)
        call allocate_traces(m%traces,2*size(tmp),m%len)
        m%traces(1:size(tmp)) = tmp
        call deallocate_traces(tmp)
    end subroutine boost_traces

    !> initialize custom mpi data type \c point_mpitype that
    !> represents \c point_type
    subroutine build_point_mpitype()
        integer,parameter :: n_blocks=5
        integer lengths(n_blocks),types(n_blocks)
        integer(kind=MPI_ADDRESS_KIND) :: base,addrs(n_blocks),displs(n_blocks)
        type(point_type) :: sample(2)
        integer :: ierror

        lengths(1) = 1          ! start of point in memory
        types(1) = MPI_LB
        call mpi_get_address(sample(1),addrs(1),ierror)

        lengths(2) = 1          ! cum
        types(2) = MPI_REAL8
        call mpi_get_address(sample(1)%cum,addrs(2),ierror)

        lengths(3) = 1          ! cum2
        types(3) = MPI_REAL8
        call mpi_get_address(sample(1)%cum2,addrs(3),ierror)

        lengths(4) = 1          ! n
        types(4) = MPI_INTEGER
        call mpi_get_address(sample(1)%n,addrs(4),ierror)

        lengths(5) = 1          ! next point
        types(5) = MPI_UB
        call mpi_get_address(sample(2),addrs(5),ierror)

        call mpi_get_address(sample(1),base,ierror) ! base address
        displs = addrs - base

        call mpi_type_create_struct(n_blocks,lengths,displs,types&
             &,point_mpitype,ierror)
        call checkmpi(ierror&
             &,'build_point_mpitype(): mpi_type_create_struct() failed')

        call mpi_type_commit(point_mpitype,ierror)
        call checkmpi(ierror&
             &,'build_point_mpitype(): mpi_type_commit() failed')
    end subroutine build_point_mpitype

    !> initialize custom mpi data type \c sample_mpitype that
    !> represents \c sample_type
    subroutine build_sample_mpitype()
        integer,parameter :: n_blocks=4
        integer lengths(n_blocks),types(n_blocks)
        integer(kind=MPI_ADDRESS_KIND) :: base,addrs(n_blocks),displs(n_blocks)
        type(sample_type) :: sample(2)
        integer :: ierror

        lengths(1) = 1          ! start of sample in memory
        types(1) = MPI_LB
        call mpi_get_address(sample(1),addrs(1),ierror)

        lengths(2) = 1          ! uid
        types(2) = MPI_INTEGER
        call mpi_get_address(sample(1)%uid,addrs(2),ierror)

        lengths(3) = 1          ! x
        types(3) = MPI_REAL8
        call mpi_get_address(sample(1)%x,addrs(3),ierror)

        lengths(4) = 1          ! next sample
        types(4) = MPI_UB
        call mpi_get_address(sample(2),addrs(4),ierror)

        call mpi_get_address(sample(1),base,ierror) ! base address
        displs = addrs - base

        call mpi_type_create_struct(n_blocks,lengths,displs,types&
             &,sample_mpitype,ierror)
        call checkmpi(ierror&
             &,'build_sample_mpitype(): mpi_type_create_struct() failed')

        call mpi_type_commit(sample_mpitype,ierror)
        call checkmpi(ierror&
             &,'build_sample_mpitype(): mpi_type_commit() failed')
    end subroutine build_sample_mpitype

    !> initialize custom mpi data type \c trace_mpitype that
    !> represents \c trace_type without the allocatable \c pos element
    subroutine build_trace_mpitype()
        integer,parameter :: n_blocks=4
        integer lengths(n_blocks),types(n_blocks)
        integer(kind=MPI_ADDRESS_KIND) :: base,addrs(n_blocks),displs(n_blocks)
        type(trace_type) :: sample(2)
        integer :: ierror

        lengths(1) = 1          ! start of trace in memory
        types(1) = MPI_LB
        call mpi_get_address(sample(1),addrs(1),ierror)

        lengths(2) = 1          ! uid
        types(2) = MPI_INTEGER
        call mpi_get_address(sample(1)%uid,addrs(2),ierror)

        lengths(3) = 1          ! len
        types(3) = MPI_INTEGER
        call mpi_get_address(sample(1)%len,addrs(3),ierror)

        lengths(4) = 1          ! next trace
        types(4) = MPI_UB
        call mpi_get_address(sample(2),addrs(4),ierror)

        ! the pos field is not included here since it's allocatable,
        ! so the offset to the other elements would not be constant in
        ! general

        call mpi_get_address(sample(1),base,ierror) ! base address
        displs = addrs - base

        call mpi_type_create_struct(n_blocks,lengths,displs,types&
             &,trace_mpitype,ierror)
        call checkmpi(ierror&
             &,'build_trace_mpitype(): mpi_type_create_struct() failed')

        call mpi_type_commit(trace_mpitype,ierror)
        call checkmpi(ierror&
             &,'build_trace_mpitype(): mpi_type_commit() failed')
    end subroutine build_trace_mpitype

    !> stop tracing a specific object in case it was sampled before
    !>
    !> \param[in] m handle to msd object
    !>
    !> \param[in] uid unique object id
    subroutine close_trace(m,uid)
        integer,intent(in) :: m,uid

        call sample(m,uid,-1.0_rk)
    end subroutine close_trace

    !> deallocates an array of traces and the position arrays contained
    !>
    !> \param[in,out] t vector of traces to deallocate
    subroutine deallocate_traces(t)
        type(trace_type),allocatable,intent(inout) :: t(:)
        integer :: i

        do i=1,size(t)
           deallocate (t(i)%pos)
        end do
        deallocate (t)
    end subroutine deallocate_traces

    !> delete a single msd checkpoint
    !>
    !> \param[in] m msd object
    !>
    !> \param[in] last timestep of the checkpoint to delete
    subroutine delete_single_checkpoint(m,last)
        type(mean_square_displacement_type),intent(in) :: m
        integer,intent(in) :: last
        character(len=1024) :: filename

        call lbe_make_filename_cp(filename,'msd-checkpoint-'//m%name,'.xdr'&
             &,last)
        call lbe_delete_file(filename)
    end subroutine delete_single_checkpoint

    !> send content of sample buffer to the appropriate processes
    !>
    !> \param[in,out] m msd object
    !>
    !> \param[out] n number of samples returned in \c rbuf
    !>
    !> \param[in,out] rbuf buffer to receive the samples on each
    !> process, allocated by \c distribute_samples but must be freed
    !> by the caller
    subroutine distribute_samples(m,n,rbuf)
        type(mean_square_displacement_type),intent(inout) :: m
        integer,intent(out) :: n
        type(sample_type),allocatable,intent(inout) :: rbuf(:)
        type(sample_type),allocatable :: sbuf(:)
        integer,allocatable :: rcounts(:),rdispls(:),scounts(:),sdispls(:)
        integer :: i,ierror,j,p,stat

        allocate (rcounts(0:nprocs-1),rdispls(0:nprocs-1),sbuf(m%n_samples)&
             &,scounts(0:nprocs-1),sdispls(0:nprocs-1),stat=stat)
        call check_allocate(stat,'distribute_samples(): '&
             &//'rcounts,rdispls,sbuf,scounts,sdispls')

        j = 0
        scounts = 0
        do p=0,nprocs-1
           do i=1,m%n_samples
              if (mod(m%samples(i)%uid,nprocs)==p) then
                 j = j+1
                 sbuf(j) = m%samples(i)
                 scounts(p) = scounts(p)+1
              end if
           end do
        end do

        if (j/=m%n_samples) call error('distribute_samples(): '&
             &//'failed to uniquely assign all samples to processes')

        call MPI_Alltoall(scounts,1,MPI_INTEGER&
             &,rcounts,1,MPI_INTEGER,MPI_COMM_WORLD,ierror)

        n = sum(rcounts)
        allocate (rbuf(n),stat=stat)
        call check_allocate(stat,'distribute_samples(): rbuf')

        call calculate_displacements(scounts,sdispls)
        call calculate_displacements(rcounts,rdispls)

        call MPI_Alltoallv(sbuf,scounts,sdispls,sample_mpitype&
             &,rbuf,rcounts,rdispls,sample_mpitype,MPI_COMM_WORLD,ierror)

        deallocate (rcounts,rdispls,sbuf,scounts,sdispls)
        m%n_samples = 0
    end subroutine distribute_samples

    !> reduce accumulated msd data from all nodes and dump it to file,
    !> reset data buffers
    !>
    !> \param[in] m handle to msd object
    subroutine dump_msd(m)
        type(mean_square_displacement_type),intent(inout) :: m
        integer,parameter :: mfu=22
        character(len=1024) :: mfn
        type(point_type),allocatable :: datasum(:)
        integer :: i,ierror,stat
        real(kind=rk) :: mean,n,stdev

        if (myrankc==0) then
           allocate (datasum(m%len),stat=stat)
           call check_allocate(stat,'dump_msd(): datasum')
        end if

        call MPI_Reduce(m%data,datasum,m%len,point_mpitype&
             &,sum_point_mpiop,0,comm_cart,ierror)

        if (myrankc==0) then
           call lbe_make_filename_output(mfn,'msd-'//m%name,'.asc',nt)
           open (unit=mfu,file=mfn,status='REPLACE',action='WRITE',recl=160)
           write (unit=mfu,fmt=&
                &'("# t [ts]             msd      stdev(msd)   #samples")')
           do i=1,m%len
              n = real(datasum(i)%n,kind=rk)
              mean = datasum(i)%cum/n
              stdev = sqrt((datasum(i)%cum2 - n*mean*mean) / (n - 1.0_rk))
!!! Replace the above with the following two lines and change
!!! update_msd_from_traces() accordingly (see comment there) and you
!!! have Welford's method for the standard deviation which numerically
!!! is more stable but not (yet) parallelized, so results will be
!!! wrong unless the code is run serially.
!!$              mean = datasum(i)%cum
!!$              stdev = sqrt(datasum(i)%cum2/(n-1.0_rk))
              write (unit=mfu,fmt='(I8,2(X,ES15.8),X,I10)') &
                   &(i-1)*m%n_sample,mean,stdev,datasum(i)%n
           end do
           close (mfu)
           deallocate (datasum)
        end if

        call reset_data(m)
        call reset_traces(m)
    end subroutine dump_msd

    !> writes the necessary information for one msd object into a
    !> checkpoint file
    !>
    !> \param[in] m msd object
    subroutine dump_single_checkpoint(m)
        type(mean_square_displacement_type),intent(in) :: m
        integer,parameter :: tag=1
        type(point_type),allocatable :: datasum(:)
        type(trace_type),allocatable :: rbuf(:),sbuf(:)
        real(kind=rk),allocatable :: prbuf(:),psbuf(:)
        character(len=1024) :: fn
        integer :: file_id,i,ierror,j,k,n_active,n_active_max,n_active_sum&
             &,n_pos,n_pos_max,n_recv,p,stat,status(MPI_STATUS_SIZE)

#ifndef USEXDRF
        call error('lbe_mean_square_displacement_mod: '&
             &//'dump_single_checkpoint() requires still XDRF')
#else
        if (myrankc==0) then
           allocate (datasum(m%len),stat=stat)
           call check_allocate(stat,'dump_single_checkpoint(): datasum')
        end if

        call MPI_Reduce(m%data,datasum,m%len,point_mpitype&
             &,sum_point_mpiop,0,comm_cart,ierror)

        if (myrankc==0) then
           call lbe_make_filename_cp(fn,'msd-checkpoint-'//m%name,'.xdr',nt)
           call xdrfopen(file_id,fn,'w',ierror)
           call check_xdrfopen(ierror,fn)

           ! dump data
           data: do i=1,m%len
              call lbe_xdrfpoint_type(file_id,datasum(i))
           end do data
           deallocate (datasum)
        end if

        ! now dump traces---first find global total number of active
        ! traces (determines file size) and maximum number of active
        ! traces per process (determines recv buffer size)
        n_active = count(m%traces(1:m%n_traces)%len>0)
        call MPI_Reduce(n_active,n_active_sum,1,MPI_INTEGER,MPI_SUM,0,comm_cart&
             &,ierror)
        call MPI_Reduce(n_active,n_active_max,1,MPI_INTEGER,MPI_MAX,0,comm_cart&
             &,ierror)

        allocate (sbuf(n_active),stat=stat)
        call check_allocate(stat,'dump_single_checkpoint(): sbuf')

        ! fill send buffer and count the local number of valid
        ! positions in all traces and the maximum of it among all
        ! tasks. This is needed for the buffers which are used to
        ! communicate the positions separately since they are not part
        ! of trace_mpitype.
        j = 0
        n_pos = 0
        do i=1,m%n_traces
           if (m%traces(i)%len>0) then
              j = j+1
              sbuf(j)%uid = m%traces(i)%uid
              sbuf(j)%len = m%traces(i)%len

              n_pos = n_pos+m%traces(i)%len
           end if
        end do

        allocate (psbuf(n_pos),stat=stat)
        call check_allocate(stat,'dump_single_checkpoint(): psbuf')

        ! fill position send buffer
        k = 0
        do i=1,m%n_traces
           do j=1,m%traces(i)%len
              k = k+1
              psbuf(k) = m%traces(i)%pos(j)
           end do
        end do

        call MPI_Reduce(n_pos,n_pos_max,1,MPI_INTEGER,MPI_MAX,0,comm_cart&
             &,ierror)

        if (myrankc==0) then
           ! total number of traces that will follow in the file
           call xdrfint(file_id,n_active_sum,ierror)

           allocate (rbuf(n_active_max),prbuf(n_pos_max),stat=stat)
           call check_allocate(stat,'dump_single_checkpoint(): rbuf,prbuf')
        end if

        ! send all traces to the root process, don't use collective
        ! communication in order to save memory.
        tasks: do p=0,nprocs-1
           ! prevent all non-root ranks from sending at the same time
           call MPI_Barrier(MPI_COMM_WORLD,ierror)

           rank0: if (myrankc/=0) then
              if (p==myrankc) then
                 call MPI_Send(sbuf,n_active,trace_mpitype,0,tag,comm_cart&
                      &,ierror)
                 call MPI_Send(psbuf,n_pos,LBE_REAL,0,tag,comm_cart,ierror)
              end if
           else rank0
              if (p==myrankc) then
                 rbuf(1:size(sbuf)) = sbuf
                 prbuf(1:size(psbuf)) = psbuf
                 n_recv = size(sbuf)
              else
                 call MPI_Recv(rbuf,n_active_max,trace_mpitype,p,tag,comm_cart&
                      &,status,ierror)
                 call MPI_Get_count(status,trace_mpitype,n_recv,ierror)
                 call MPI_Recv(prbuf,n_pos_max,LBE_REAL,p,tag,comm_cart&
                      &,MPI_STATUS_IGNORE,ierror)
              end if

              ! dump the received traces
              k = 0
              do i=1,n_recv
                 call xdrfint(file_id,rbuf(i)%uid,ierror)
                 call xdrfint(file_id,rbuf(i)%len,ierror)
                 do j=1,rbuf(i)%len
                    k = k+1
                    call xdrfdouble(file_id,prbuf(k),ierror)
                 end do
              end do
           end if rank0
        end do tasks

        if (myrankc==0) then
           call xdrfclose(file_id,ierror)
           deallocate (rbuf,prbuf)
        end if
        deallocate (sbuf,psbuf)
#endif
    end subroutine dump_single_checkpoint

    !> build and register custom mpi data types and custom reduction
    !> operations use by \c lbe_mean_square_displacement_mod
    subroutine init_msd_mpi()
        integer :: ierror

        call build_point_mpitype()
        call build_sample_mpitype()
        call build_trace_mpitype()

        call MPI_Op_create(sum_point,.true.,sum_point_mpiop,ierror)
    end subroutine init_msd_mpi

#ifdef USEXDRF
    !> reads or writes a msd data point from/to an XDR file
    !>
    !> \param[inout] handle XDR file identifier as obtained from \c xdrfopen()
    !>
    !> \param[inout] point msd \c point_type object
    !>
    !> Whether data is read or written is determined by how \c handle
    !> was opened. In this sense, \c lbe_xdrflbe_site() works the same
    !> way as the xdrf routines for builtin types (for instance \c
    !> xdrfdouble()).
    subroutine lbe_xdrfpoint_type(handle,point)
        integer,intent(inout) :: handle
        type(point_type),intent(inout) :: point
        integer :: ierror

        call xdrfdouble(handle,point%cum,ierror)
        call xdrfdouble(handle,point%cum2,ierror)
        call xdrfint(handle,point%n,ierror)
   end subroutine lbe_xdrfpoint_type
#endif

    !> delete all files of an msd checkpoint
    !>
    !> \param[in] last timestep of the checkpoint to delete
    subroutine msd_delete_checkpoint(last)
        integer,intent(in) :: last
        integer :: i

        if (size(msds)==0) return

        call log_msg("  Deleting msd checkpoint files ...")
        do i=1,size(msds)
           if (myrankc==0) call delete_single_checkpoint(msds(i),last)
        end do
    end subroutine msd_delete_checkpoint

    !> writes the necessary information of all msd objects into
    !> checkpoint files
    subroutine msd_dump_checkpoint
        integer :: i

        if (size(msds)==0) return

        call log_msg("  Writing msd checkpoint files ...")
        do i=1,size(msds)
           call log_msg('   <'//trim(msds(i)%name)//'>...')
           call dump_single_checkpoint(msds(i))
        end do
    end subroutine msd_dump_checkpoint

    !> initializes all msd objects and common global data
    !>
    !> This routine needs to be called after all msd objects were
    !> registered and before any of them is used.
    subroutine msd_init
        integer :: i,stat

        ! make sure the code further on works in a defined way also if
        ! msds is still unallocated
        if (.not.allocated(msds)) then
           allocate (msds(0),stat=stat)
           call check_allocate(stat,'msd_init(): msds')
        end if

        do i=1,size(msds)
           call allocate_msd(msds(i))
        end do

        call init_msd_mpi()

        initialized = .true.

        if (is_restoring()) call msd_restore_checkpoint()
    end subroutine msd_init

    !> restore all registered msd object from the respective
    !> checkpoint files
    subroutine msd_restore_checkpoint
        integer :: i

        if (size(msds)==0) return

        call log_msg("  Restoring msd checkpoint files ...")
        do i=1,size(msds)
           call log_msg('   -'//trim(msds(i)%name)//'...')
           call restore_single_checkpoint(msds(i))
        end do
    end subroutine msd_restore_checkpoint

    !> does everything necessary to keep msd calculation and output
    !> running during one time step
    !>
    !> \param[in] whole_N local chunk of the lattice with full halo of
    !> depth \c halo_extent
    subroutine msd_run(whole_N)
        type(lbe_site),intent(in) :: &
             &whole_N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        integer :: i

        if (size(msds)==0) return

        do i=1,size(msds)
           if (every_n_time_steps(msds(i)%n_sample)) then
              call log_msg('Mean square displacement <'//trim(msds(i)%name)&
                   &//'> : accumulating samples')
              call reduce_samples(msds(i),whole_N)
           end if
           if (every_n_time_steps(msds(i)%n_dump)) then
              call log_msg('Mean square displacement <'//trim(msds(i)%name)&
                   &//'> : dumping')
              call dump_msd(msds(i))
           end if
        end do
    end subroutine msd_run

    !> ensures that an msd object has at least a given capacity to
    !> store traces
    !>
    !> \param[in,out] m msd object
    !>
    !> \param[in] n requested minimum capacity
    subroutine provide_trace_capacity(m,n)
        type(mean_square_displacement_type),intent(inout) :: m
        integer,intent(in) :: n

        do while (size(m%traces)<n)
           call boost_traces(m)
        end do
    end subroutine provide_trace_capacity

    !> sort and append content of sample buffer into trace buffers of
    !> the appropriate processes and update the cumulative mean square
    !> displacement data with this data
    !>
    !> \param[in] m msd object
    !>
    !> \param[in] whole_N local chunk of the lattice with full halo of
    !> depth \c halo_extent
    subroutine reduce_samples(m,whole_N)
        type(mean_square_displacement_type),intent(inout) :: m
        type(lbe_site),intent(in) :: &
             &whole_N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        integer :: n_append
        type(sample_type),allocatable :: buf(:)

        call distribute_samples(m,n_append,buf)
        call append_samples_to_traces(m,n_append,buf,whole_N)
        call update_msd_from_traces(m)
    end subroutine reduce_samples

    !> obtain a handle to a newly created msd object
    !>
    !> \param[in] name human readable name, also used in file names
    !>
    !> \param[in] n_dump time interval for dumping
    !>
    !> \param[in] n_sample time interval for sampling positions
    !>
    !> \param[in] len time to record means square displacement in
    !> units of \c n_sample LB time steps
    !>
    !> \param[out] m returns handle to new msd object. A valid handle
    !> is a non-zero positive integer number.
    !>
    !> \param[in] periodic \c .true. if displacement direction is
    !> periodic [optional, defaults to \c .false. ]
    !>
    !> \param[in] d_periodic periodic dimension of displacement
    !> direction [optional, defaults to \c huge(d_periodic) ]
    !>
    !> \param[in] complete_only ignore non-complete traces [optional,
    !> defaults to \c .false. ]
    !>
    !> \param[in] central_startpos_only ignore traces that don not
    !> start within \c [min_startpos,max_startpos] [optional, defaults
    !> to \c .false. ]
    !>
    !> \param[in] min_startpos minimum accepted start position with \c
    !> central_startpos_only [optional, defaults to \c
    !> huge(min_startpos) ]
    !>
    !> \param[in] max_startpos maximum accepted start position with \c
    !> central_startpos_only [optional, defaults to \c
    !> huge(max_startpos) ]
    !>
    !> \param[in] subtract_drift if set to 1, 2, or 3, adjust sampled
    !> positions according to the x-, y-, or z-component of a possible
    !> drift of the whole system [optional, defaults to \c -1 ]
    subroutine register_msd(name,n_dump,n_sample,length,m,periodic,d_periodic&
         &,complete_only,central_startpos_only,min_startpos,max_startpos&
         &,subtract_drift)
        character(len=*),intent(in) :: name
        integer,intent(in) :: n_dump,n_sample,length
        integer,intent(out) :: m
        logical,intent(in),optional :: periodic
        real(kind=rk),intent(in),optional :: d_periodic
        logical,intent(in),optional :: complete_only
        logical,intent(in),optional :: central_startpos_only
        real(kind=rk),intent(in),optional :: min_startpos,max_startpos
        integer,intent(in),optional :: subtract_drift
        type(mean_square_displacement_type),allocatable :: tmp(:)
        integer :: stat

        if (initialized) call error('register_msd(): msd_init() was called '&
             &//'already, further msd objects must not be created anymore. ')

        ! make sure the code further on works in a defined way also if
        ! msds is still unallocated
        if (.not.allocated(msds)) then
           allocate (msds(0),stat=stat)
           call check_allocate(stat,'register_msd(): msds')
        end if

        ! append new msd object to msds
        allocate (tmp(size(msds)),stat=stat)
        call check_allocate(stat,'register_msd(): tmp')
        tmp = msds
        deallocate (msds)
        allocate (msds(size(tmp)+1),stat=stat)
        call check_allocate(stat,'register_msd(): msds')
        msds(1:size(tmp)) = tmp
        deallocate (tmp)

        ! initialize new msd object
        m = size(msds)
        msds(m)%name = name
        msds(m)%n_dump = n_dump
        msds(m)%n_sample = n_sample
        msds(m)%len = length

        ! initialize optional elements
        msds(m)%periodic = .false.
        if (present(periodic)) msds(m)%periodic = periodic

        msds(m)%d_periodic = huge(msds(m)%d_periodic)
        if (present(d_periodic)) msds(m)%d_periodic = d_periodic

        msds(m)%complete_only = .false.
        if (present(complete_only)) msds(m)%complete_only = complete_only

        msds(m)%central_startpos_only = .false.
        if (present(central_startpos_only)) &
             &msds(m)%central_startpos_only = central_startpos_only

        msds(m)%min_startpos = -huge(msds(m)%min_startpos)
        if (present(min_startpos)) msds(m)%min_startpos = min_startpos

        msds(m)%max_startpos = huge(msds(m)%max_startpos)
        if (present(max_startpos)) msds(m)%max_startpos = max_startpos

        msds(m)%subtract_drift = -1
        if (present(subtract_drift)) msds(m)%subtract_drift = subtract_drift
    end subroutine register_msd

    !> clear all accumulated msd data
    !>
    !> \param[in,out] m msd object
    subroutine reset_data(m)
        type(mean_square_displacement_type),intent(inout) :: m

        m%data(:)%cum = 0.0_rk
        m%data(:)%cum2 = 0.0_rk
        m%data(:)%n = 0
    end subroutine reset_data

    !> clear all accumulated trace data
    !>
    !> \param[in,out] m msd object
    !>
    !> The traces are not removes, just all existing position data is
    !> invalidated. If the set of uids that lead to valid samples
    !> changes a lot during a simulation, it is maybe more efficient
    !> to really remove the traces themselves.
    subroutine reset_traces(m)
        type(mean_square_displacement_type),intent(inout) :: m
        integer :: i

        do i=1,m%n_traces
           m%traces(i)%len = 0
        end do
    end subroutine reset_traces

    !> read the necessary information for one msd object from a
    !> checkpoint file
    !>
    !> \param[in] m msd object
    subroutine restore_single_checkpoint(m)
        type(mean_square_displacement_type),intent(inout) :: m
        integer,parameter :: n_read_max=10
        character(len=1024) :: fn
        type(trace_type),allocatable :: read_buf(:),rbuf(:),sbuf(:)
        real(kind=rk),allocatable :: prbuf(:),psbuf(:)
        integer :: c,file_id,i,ierror,j,n_chunks,n_recv,n_traces_total,p,pj&
             &,pn_recv,pn_send,stat,t,t_end
        integer :: counts(0:nprocs-1),pcounts(0:nprocs-1)
        integer :: displs(0:nprocs-1),pdispls(0:nprocs-1)

#ifndef USEXDRF
        call error('lbe_mean_square_displacement_mod: '&
             &//'restore_single_checkpoint() requires still XDRF')
#else
        allocate (prbuf(1),stat=stat)
        call check_allocate(stat,'restore_single_checkpoint(): prbuf')

        if (myrankc==0) then
           allocate (psbuf(1),stat=stat)
           call check_allocate(stat,'restore_single_checkpoint(): psbuf')

           call allocate_traces(read_buf,n_read_max,m%len)
           call allocate_traces(sbuf,n_read_max,0) ! no positions needed here

           call lbe_make_filename_restore(fn,'msd-checkpoint-'//m%name,'.xdr')
           call xdrfopen(file_id,fn,'r',ierror)
           call check_xdrfopen(ierror,fn)

           ! restore data
           data: do i=1,m%len
              call lbe_xdrfpoint_type(file_id,m%data(i))
           end do data

           ! total number of traces that will follow in the file
           call xdrfint(file_id,n_traces_total,ierror)
        else
           ! read msd data on root process only and reset it
           ! elsewhere, since later it is accumulated on root, anyway
           call reset_data(m)
        end if

        call MPI_Bcast(n_traces_total,1,MPI_INTEGER,0,comm_cart,ierror)
        n_chunks = (n_traces_total+n_read_max-1)/n_read_max

        ! now restore traces
        m%n_traces = 0
        chunks: do c=1,n_chunks
           rank0: if (myrankc==0) then
              t_end = min(n_read_max,n_traces_total-(c-1)*n_read_max)

              ! fill read buffer
              pn_send = 0       ! number of positions read and to send
              traces: do t=1,t_end
                 call xdrfint(file_id,read_buf(t)%uid,ierror)
                 call xdrfint(file_id,read_buf(t)%len,ierror)
                 positions: do i=1,read_buf(t)%len
                    call xdrfdouble(file_id,read_buf(t)%pos(i),ierror)
                 end do positions
                 pn_send = pn_send+read_buf(t)%len
              end do traces

              ! position send buffer
              if (size(psbuf)<pn_send) then
                 deallocate (psbuf)
                 allocate (psbuf(pn_send),stat=stat)
                 call check_allocate(stat,'restore_single_checkpoint(): psbuf')
              end if

              ! sort traces into owning processes
              j = 0                ! last trace pointer
              pj = 0               ! last position pointer
              counts = 0           ! number of traces per task
              pcounts = 0          ! number of positions per task
              processes: do p=0,nprocs-1
                 do t=1,t_end
                    if (mod(read_buf(t)%uid,nprocs)==p) then
                       j = j+1
                       sbuf(j)%uid = read_buf(t)%uid
                       sbuf(j)%len = read_buf(t)%len
                       do i=1,read_buf(t)%len
                          pj = pj+1
                          psbuf(pj) = read_buf(t)%pos(i)
                       end do
                       counts(p) = counts(p)+1
                       pcounts(p) = pcounts(p)+read_buf(t)%len
                    end if
                 end do
              end do processes
              call calculate_displacements(counts,displs)
              call calculate_displacements(pcounts,pdispls)
           end if rank0

           ! scatter traces
           call MPI_Scatter(counts,1,MPI_INTEGER,n_recv,1,MPI_INTEGER,0&
                &,comm_cart,ierror)
           call provide_trace_capacity(m,m%n_traces+n_recv)
           call MPI_Scatterv(sbuf,counts,displs,trace_mpitype&
                &,m%traces(m%n_traces+1:m%n_traces+n_recv),n_recv,trace_mpitype&
                &,0,comm_cart,ierror)

           ! scatter positions separately
           call MPI_Scatter(pcounts,1,MPI_INTEGER,pn_recv,1,MPI_INTEGER,0&
                &,comm_cart,ierror)
           if (size(prbuf)<pn_recv) then
              deallocate (prbuf)
              allocate (prbuf(pn_recv),stat=stat)
              call check_allocate(stat,'restore_single_checkpoint(): prbuf')
           end if
           call MPI_Scatterv(psbuf,pcounts,pdispls,LBE_REAL&
                &,prbuf,pn_recv,LBE_REAL,0,comm_cart,ierror)

           ! combine positions with traces
           pj = 0
           do t=m%n_traces+1,m%n_traces+n_recv
              do i=1,m%traces(t)%len
                 pj = pj+1
                 m%traces(t)%pos(i) = prbuf(pj)
              end do
           end do

           m%n_traces = m%n_traces+n_recv
        end do chunks

        if (myrankc==0) then
           call xdrfclose(file_id,ierror)
           call deallocate_traces(read_buf)
           call deallocate_traces(sbuf)
           deallocate (psbuf)
        end if
        deallocate (prbuf)

        ! build uid2trace
        call Mii_clear(m%uid2trace)
        call check_allocate(Mii_provide_capacity(m%uid2trace,m%n_traces)&
             &,'restore_single_checkpoint(): m%uid2trace'&
             &//' (Mii_provide_capacity())')
        do t=1,m%n_traces
           call Mii_preinsert(m%uid2trace,m%traces(t)%uid,t)
        end do
        call Mii_commit(m%uid2trace)
#endif
    end subroutine restore_single_checkpoint

    !> sample an 1D object position for later calculation of its mean
    !> square displacement
    !>
    !> \param[in] m handle to msd object
    !>
    !> \param[in] uid unique object id
    !>
    !> \param[in] x 1D object position
    !>
    !> Globally, all \c uid values passed in between two calls to \c
    !> reduce_samples() must be unique.
    subroutine sample(m,uid,x)
        integer,intent(in) :: m,uid
        real(kind=rk),intent(in) :: x

        if (msds(m)%n_samples==size(msds(m)%samples)) &
             &call boost_samples(msds(m))

        msds(m)%n_samples = msds(m)%n_samples + 1
        msds(m)%samples(msds(m)%n_samples)%uid = uid
        msds(m)%samples(msds(m)%n_samples)%x = x
    end subroutine sample

    !> custom mpi reduction operation to sum objects of type \c
    !> point_type
    !>
    !> \param[in] invec vector holding input data
    !>
    !> \param[in,out] invec vector holding other input data and
    !> expected to hold result on exit
    !>
    !> \param[in] len number of \c point_type elements to combine
    !>
    !> \param[in] type data type handle, passed by MPI but ignored
    !> here since the routine always works on \c point_type
    !>
    !> According to the requirements of MPI custom reuction
    !> operations, this routine combines the vectors \c invec and \c
    !> inoutvec element-wise and store the result in \c inoutvec.
    subroutine sum_point(invec,inoutvec,length,type)
        type(point_type),intent(in) :: invec(length)
        type(point_type),intent(inout) :: inoutvec(length)
        integer,intent(in) :: length,type
        integer i

        do i=1,length
           inoutvec(i)%cum = invec(i)%cum + inoutvec(i)%cum
           inoutvec(i)%cum2 = invec(i)%cum2 + inoutvec(i)%cum2
           inoutvec(i)%n = invec(i)%n + inoutvec(i)%n
        end do
    end subroutine sum_point

    !> accumulate data from trace list in \c cum, \c cum2, and \c n
    !>
    !> \param[in,out] m msd object
    subroutine update_msd_from_traces(m)
        type(mean_square_displacement_type),intent(inout) :: m
        real(kind=rk) :: sd
        integer :: i,t
!!! needed only for Welford's method, see comment in dump_msd()
!!$        real(kind=rk) :: delta

        do i=1,m%n_traces
           if (m%complete_only.and.m%traces(i)%len/=m%len) cycle
           if (m%central_startpos_only.and.&
                &(m%traces(i)%pos(m%traces(i)%len)<m%min_startpos&
                &.or.m%traces(i)%pos(m%traces(i)%len)>=m%max_startpos)) cycle

           do t=1,m%traces(i)%len
              sd = (m%traces(i)%pos(t)-m%traces(i)%pos(1))**2
              m%data(t)%n = m%data(t)%n + 1
              m%data(t)%cum = m%data(t)%cum + sd
              m%data(t)%cum2 = m%data(t)%cum2 + sd*sd
!!! For Welford's method replace the above two lines with the
!!! following three, read comment in dump_msd() first!
!!$              delta = sd - m%data(t)%cum
!!$              m%data(t)%cum = m%data(t)%cum + delta/real(m%data(t)%n,kind=rk)
!!$              m%data(t)%cum2 = m%data(t)%cum2 + delta*(sd-m%data(t)%cum)
           end do
        end do
    end subroutine update_msd_from_traces

end module lbe_mean_square_displacement_mod
