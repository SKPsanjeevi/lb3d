#include "lbe.h"

!>Code related to the parallelization of the calculation
module lbe_parallel_module
  use lbe_globals_module
  use lbe_helper_module, only: density,makedatestr,present_and_true
  use lbe_log_module
  use lbe_parms_module, only: boundary,dbg_n_start_mpi_debug&
       &,dbg_report_topology,inv_fluid,nt,nx,ny,nz
  use lbe_types_module

  implicit none
  include 'mpif.h'
  private

  public halo, nprocs, myrankw, cdims, nnprocs, comm_cart, ccoords, start, tnx, tny, tnz
  public tag_px, tag_py, tag_pz, tag_mx, tag_my, tag_mz
  public ls_mpitype, us_mpitype, lr_mpitype, ur_mpitype
  public start_time_t, end_time_t
#ifdef USEHDF
  public hname, startsimul
#endif

  public calculate_displacements, gather_rock_state, stats_rk
  public avg_fluid_density_facez, get_send_partner, owning_rankc
  public InitializeMPIcomms, ReorderMPIcomms, lbe_divide_grid, FinalizeMPIcomms, lbe_parallel_init ! For lbe.F90
  public lbe_setup_halo_exchange_new, lbe_halo_exchange_new
  public find_topology, recv_rock_par
  public checkmpi, check_allocate, Abend
  public send_final_lattice, recv_final_lattice, set_cpcoords, send_rock_bit, send_rock_release ! For lbe_io.F90
  public halo_exchange, build_all_chunk_mpitypes
  public rock_halo
#ifndef NOSURFACTANT
  public lbe_adv_dipole_exchange, lbe_all_dipole_exchange
#endif
#ifdef DEBUG_MPI
  public log_mpi
#endif
#ifdef VARTAU
  public vartau_halo,vartau_mpitype,build_vartau_site_mpitype
#endif
#ifdef LOCALBC
  public acccoef_halo,acccoef_mpitype,build_local_acccoef_site_mpitype
#endif     

  !> This type contains all data needed to do a call halo_exchange(N, h) type halo exchange.
  type halo
    sequence
    integer, dimension(3) :: ls      !< Lower send MPI type
    integer, dimension(3) :: us      !< Upper send MPI type
    integer, dimension(3) :: lr      !< Lower receive MPI type
    integer, dimension(3) :: ur      !< Upper receive MPI type
    integer               :: mpitype !< Base MPI type
    integer               :: extent  !< Halo extent
  end type halo

  type(halo), save :: rock_halo

  integer, save :: start_time_t, end_time_t

  !> \{
  !> \name Variables holding info about the processor topology
  !>
  !> Set to the number of CPUs taking part in the simulation.
  integer, save :: nprocs
  integer, save :: myrankw            !< rank of the CPU (MPI_COMM_WORLD)
  ! integer, save :: myrankc            !< rank of the CPU (Comm_Cart)
  integer, save :: cdims(nd) = 0      !< dimensions of the CPU topology.
  integer, save :: nnprocs(nd, 2)     !< ranks of the nearest-neighbour CPUs

  integer, save :: Comm_Cart      !< The Cartesian CPU topology communicator
  integer, save :: ccoords(nd) = 0    !< My coords in the Cartesian grid.
  integer, save :: start(nd) !< global coordinates of the start of the CPU's subdomain
  integer, save :: tnx,tny,tnz        !< Total size of lattice.
  !> \}

  !> Determines whether or not the Cartesian topology is periodic in
  !> each given dimension.
  !>
  !> \note This may change for different boundary conditions.
  logical, dimension(3), save :: periodic_p = .true.

  !> \{
  !> \name MPI Derived Datatypes
  !>
  !> These datatypes describe the surface of interface between
  !> two processors' domains.
  ! integer, save :: interf_px,interf_mx
  ! integer, save :: interf_py,interf_my
  ! integer, save :: interf_pz,interf_mz
  !> \}

  !> MPI tag for messages related to initialisation
  integer, parameter :: tag_init = 1
  !> \{
  !> \name MPI tags used in the halo exchange
  integer, parameter :: tag_px = 2,tag_mx = 3
  integer, parameter :: tag_py = 4,tag_my = 5
  integer, parameter :: tag_pz = 6,tag_mz = 7
  !> \}
  integer, parameter :: tag_rock = 8
  !> MPI tag for messages related to postprocessing
  integer, parameter :: tag_post = 15

  !> \{
  !> \name EYB HDF5 added
  character(len=MPI_MAX_PROCESSOR_NAME), save :: hname !< For hostname
  character(len=24), save :: startsimul                !< For starttime
  !> \}

  !> \{
  !> \name MPI data types---also read comment to \c lbe_halo_exchange_new()
  !> mpi datatypes representing the lower (l)/upper (u) part of the halo that
  !> is sent (s)/received (r) during the exchange for each direction (indices
  !> 1/2/3 represent x/y/z). Datatype definitions are relative to \c N as a
  !> whole.
  integer, save :: ls_mpitype(3),us_mpitype(3),lr_mpitype(3),ur_mpitype(3)
  !> \}

!!$  !> \{
!!$  !> \name for BC periodic_inflow (should be moved there as soon as possible)
!!$  integer,save :: last_periodic_z = 20 !< last z position of periodic part
!!$  real(kind=rk),save :: tsize_pi(3) !< dimensions of periodic sub-volume
!!$  real(kind=rk),save :: maxpos_pi(3) !< upper bounds of periodic part
!!$  real(kind=rk),save :: border_npi(2,3) !< dim's local non-periodic lo/hi bounds
!!$  real(kind=rk),save :: border_pi(2,3) !< dim's local periodic lo/hi bounds
!!$  !> \}

  ! Intermediate Knudsen adjustment additional datasets
#ifdef VARTAU
type(halo) :: vartau_halo
integer :: vartau_mpitype
#endif
#ifdef LOCALBC
type(halo) :: acccoef_halo
integer :: acccoef_mpitype
#endif     





contains

#ifdef DEBUG_MPI
subroutine log_mpi(msg)
  character(len=*), intent(in) :: msg
  character(len=256) :: msgbuf
  integer :: mpierr, mpistatus, i
  integer :: tag = 0

  if ( nt .ge. dbg_n_start_mpi_debug )  then
    write(msgbuf,'(A)') trim(msg)
    if ( myrankc .gt. 0 ) then
      ! Magic number '256' corresponds to the length of msgstr of course
      call MPI_Send(msgbuf, 256, MPI_CHARACTER, 0, tag, Comm_cart, mpierr)
    else
      write(msgstr,'("DEBUG_MPI @ t = ",I0," - Received from ",I6.6," : ",A)') nt, 0, trim(msgbuf)
      call log_msg(msgstr)
      do i=1, nprocs-1
        ! Magic number '256' corresponds to the length of msgstr of course
        call MPI_Recv(msgbuf, 256, MPI_CHARACTER, i, tag, Comm_cart, mpistatus, mpierr)
        write(msgstr,'("DEBUG_MPI @ t = ",I0," - Received from ",I6.6," : ",A)') nt, i, trim(msgbuf)
        call log_msg(msgstr)
      enddo
    endif
    call log_msg_ws("DEBUG_MPI BARRIER")
    call MPI_Barrier(Comm_cart,mpierr)
  endif
end subroutine log_mpi
#endif

    !> calculates the mass density averaged over all nodes of a global z-plane
    !>
    !> \param[in] N local lattice chunk ("lbe" N)
    !>
    !> \param[in] z global z coordinate of averaging plane
    !>
    !> \returns mass density
    function avg_fluid_density_facez(N,z)
        real(kind=rk) :: avg_fluid_density_facez
        type(lbe_site),intent(in) :: N(0:,0:,0:)
        integer,intent(in) :: z
        integer :: cnt,global_count,lx,ly,lz,ierror
        real(kind=rk) :: lsum,global_sum

        lz = z+1-start(3)
        cnt = 0
        lsum = 0.0_rk
        if (1<=lz.and.lz<=nz) then
           do ly=1,ny
              do lx=1,nx
                 if (N(lx,ly,lz)%rock_state==0.0_rk) then
                    lsum = lsum+density(N(lx,ly,lz))
                    cnt = cnt+1
                 end if
              end do
           end do
        end if

        call mpi_allreduce(lsum,global_sum,1,LBE_REAL,MPI_SUM,MPI_COMM_WORLD&
             &,ierror)
        call mpi_allreduce(cnt,global_count,1,MPI_INTEGER,MPI_SUM&
             &,MPI_COMM_WORLD,ierror)
        avg_fluid_density_facez = global_sum/real(global_count,kind=rk)
    end function avg_fluid_density_facez

    !> Builds a custom mpi datatype that represents the part of  N  specified by
    !> the coordinate intervals  xr(2),yr(2),zr(2)  and stores it in  lcmt .
    subroutine build_lattice_chunk_mpitype(N,xr,yr,zr,lsmt,lcmt)
        type(lbe_site),intent(in) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        integer,intent(in) :: xr(2),yr(2),zr(2) ! chunk ranges for each dim.
        integer,intent(in) :: lsmt  ! type used for each site
        integer,intent(out) :: lcmt ! type to be built
        integer :: xrow_mt,xyplane_mt,xyzchunk_mt ! temporary mpi data types
        integer :: blocks,lengths(1),ierror
        integer(kind=MPI_ADDRESS_KIND) :: adr1,adr2,base,offset,stride
        integer(kind=MPI_ADDRESS_KIND) :: displs(1)

        ! mpi datatype for slices of N like  N(xr(1):xr(2),y,z)
        blocks = 1 + xr(2) - xr(1)
        call mpi_get_address(N(1,1,1),adr1,ierror)
        call mpi_get_address(N(2,1,1),adr2,ierror)
        stride = adr2 - adr1
        call mpi_type_create_hvector(blocks,1,stride,lsmt,xrow_mt,ierror)

        ! mpi datatype for slices of N like  N(xr(1):xr(2),yr(1):yr(2),z)
        blocks = 1 + yr(2) - yr(1)
        call mpi_get_address(N(1,1,1),adr1,ierror)
        call mpi_get_address(N(1,2,1),adr2,ierror)
        stride = adr2 - adr1
        call mpi_type_create_hvector(blocks,1,stride,xrow_mt,xyplane_mt,ierror)

        ! mpi datatype for whole chunk  N(xr(1):xr(2),yr(1):yr(2),zr(1):zr(2))
        blocks = 1 + zr(2) - zr(1)
        call mpi_get_address(N(1,1,1),adr1,ierror)
        call mpi_get_address(N(1,1,2),adr2,ierror)
        stride = adr2 - adr1
        call mpi_type_create_hvector(blocks,1,stride,xyplane_mt,xyzchunk_mt,ierror)

        ! position of the beginning of the chunk relative to the beginning of N
        call mpi_get_address(N,base,ierror)
        call mpi_get_address(N(xr(1),yr(1),zr(1)),offset,ierror)
        offset = offset - base

        ! lcmt becomes a datatype like  xyzchunk_mt  but relative to  base
        lengths = (/1/)
        displs = (/ offset /)
        call mpi_type_create_hindexed(1,lengths,displs,xyzchunk_mt,lcmt,ierror)
        call mpi_type_commit(lcmt,ierror)
    end subroutine build_lattice_chunk_mpitype

    !> Build a custom mpi datatype that represents \c lbe_site (or parts of it)
    !>
    !> \param[in] idirs directions with value \c .true. of \c
    !> n_(r|b|s) are included
    !>
    !> \param[in] id \c d is included if \c .true. (only available
    !> without \c NOSURFACTANT)
    !>
    !> \param[in] irs \c rock_state is included if \c .true.
    subroutine build_lbe_site_mpitype(idirs&
#ifndef NOSURFACTANT
         &,id&
#endif
         &,irs,lsmt)
        logical,intent(in) :: idirs(19)&
#ifndef NOSURFACTANT
             &,id&
#endif
             &,irs
        integer,intent(out) :: lsmt ! type to be built
!#ifdef SINGLEFLUID
!  integer,parameter :: n_blocks_max = 7        
!#else
!#ifdef NOSURFACTANT
!integer,parameter :: n_blocks_max = 8
!#else
!integer,parameter :: n_blocks_max = 10
!#endif
!#endif

!#ifdef COMMON_VEL_FIX
#ifndef OLD_VEL
#ifdef SINGLEFLUID
    integer,parameter :: n_blocks_max = 8
#else
#ifdef NOSURFACTANT
     integer,parameter :: n_blocks_max = 9
#else
     integer,parameter :: n_blocks_max = 12
#endif
#endif
#else

#ifdef SINGLEFLUID
    integer,parameter :: n_blocks_max = 7
#else
#ifdef NOSURFACTANT
  integer,parameter :: n_blocks_max = 8
#else
  integer,parameter :: n_blocks_max = 10
#endif
#endif

#endif




        integer lengths(n_blocks_max),types(n_blocks_max)
        integer dcount,dlengths(19),ddispls(19),dirs_mt
        integer(kind=MPI_ADDRESS_KIND) :: addrs(n_blocks_max),displs(n_blocks_max)
        integer i,ierror,n
        type(lbe_site) :: sample(2)

        dlengths = 1
        dcount = 0
        do i=1,19
           if (idirs(i)) then
              dcount = dcount+1
              ddispls(dcount) = i-1
           end if
        end do
        if (dcount>0) call MPI_Type_indexed(&
             &dcount,dlengths,ddispls,MPI_REAL8,dirs_mt,ierror)

        n=0

        n = n+1
        lengths(n) = 1
        types(n) = MPI_LB
        call MPI_Get_address(sample(1),addrs(n),ierror)

        directions: if (dcount>0) then
           n = n+1
           lengths(n) = 1          ! n_r
           types(n) = dirs_mt
           call mpi_get_address(sample(1)%n_r(1),addrs(n),ierror)

#ifndef SINGLEFLUID
           n = n+1
           lengths(n) = 1          ! n_b
           types(n) = dirs_mt
           call mpi_get_address(sample(1)%n_b(1),addrs(n),ierror)

#ifndef NOSURFACTANT
           n = n+1
           lengths(n) = 1          ! n_s
           types(n) = dirs_mt
           call mpi_get_address(sample(1)%n_s(1),addrs(n),ierror)
#endif
#endif
        end if directions

#ifndef SINGLEFLUID
#ifndef NOSURFACTANT
        if (id) then
           n = n+1
           lengths(n) = 3          ! d
           types(n) = MPI_REAL8
           call mpi_get_address(sample(1)%d(1),addrs(n),ierror)
        end if
#endif
#endif

        if (irs) then
           n = n+1
           lengths(n) = 1 ! rock_state
           types(n) = MPI_REAL8
           call mpi_get_address(sample(1)%rock_state,addrs(n),ierror)

           n = n+1
           lengths(n) = 1 ! rock_colour
           types(n) = MPI_REAL8
           call mpi_get_address(sample(1)%rock_colour,addrs(n),ierror)
           
           n = n+1
             lengths(n) = 1 ! rock_colour_r
             types(n) = MPI_REAL8
             call mpi_get_address(sample(1)%rock_colour_r,addrs(n),ierror)
           n = n+1
             lengths(n) = 1 ! rock_colour_b
             types(n) = MPI_REAL8
             call mpi_get_address(sample(1)%rock_colour_b,addrs(n),ierror)
        end if
!#ifdef COMMON_VEL_FIX
#ifndef OLD_VEL
  direction: if (dcount>0) then
             n = n+1
             lengths(n) = 1          ! n_r_pre
             types(n) = dirs_mt
             call mpi_get_address(sample(1)%n_r_pre(1),addrs(n),ierror)
  
#ifndef SINGLEFLUID
             n = n+1
             lengths(n) = 1          ! n_b_pre
             types(n) = dirs_mt
             call mpi_get_address(sample(1)%n_b_pre(1),addrs(n),ierror)
  
#ifndef NOSURFACTANT
             n = n+1
             lengths(n) = 1          ! n_s_pre
             types(n) = dirs_mt
             call mpi_get_address(sample(1)%n_s_pre(1),addrs(n),ierror)
#endif
#endif
          end if direction
#endif

        n = n+1
        lengths(n) = 1   ! next lattice site
        types(n) = MPI_UB
        call mpi_get_address(sample(2),addrs(n),ierror)

        displs(2:n) = addrs(2:n) - addrs(1)
        displs(1) = 0
        call mpi_type_create_struct(n,lengths,displs,types,lsmt,ierror)
        call mpi_type_commit(lsmt,ierror)
    end subroutine build_lbe_site_mpitype

#ifdef VARTAU
subroutine build_vartau_site_mpitype(mpitype)
  implicit none

  integer, intent(out) :: mpitype ! type to be built

#ifdef SINGLEFLUID
  integer, parameter :: count = 1
#else
#ifdef NOSURFACTANT
  integer, parameter :: count = 2
#else
  integer, parameter :: count = 3
#endif
#endif

  type(lbe_site) :: sample
  integer(kind = MPI_ADDRESS_KIND) :: addrs(count), base, displs(count)
  integer :: blocklengths(count), types(count)
  integer :: mpierror

  ! Get base address of an lbe_site
  call MPI_Get_address(sample, base, mpierror)
  DEBUG_CHECKMPI(mpierror, "build_vartau_site_mpitype, MPI_Get_address : base")

  ! Component local relaxation time

  blocklengths(1) = 1
  types(1) = LBE_REAL
  call MPI_Get_address(sample%taupos_r,addrs(1), mpierror)
  DEBUG_CHECKMPI(mpierror, "build_vartau_site_mpitype, MPI_Get_address : 1")
#ifndef SINGLEFLUID
  blocklengths(2) = 1
  types(2) = LBE_REAL
  call MPI_Get_address(sample%taupos_b,addrs(2), mpierror)
  DEBUG_CHECKMPI(mpierror, "build_vartau_site_mpitype, MPI_Get_address : 2")
#ifndef NOSURFACTANT
  blocklengths(3) = 1
  types(3) = LBE_REAL
  call MPI_Get_address(sample%taupos_s,addrs(3), mpierror)
  DEBUG_CHECKMPI(mpierror, "build_vartau_site_mpitype, MPI_Get_address : 3")
#endif
#endif

  ! Calculate displacements
  displs(:) = addrs(:) - base

  call MPI_Type_create_struct(count, blocklengths, displs, types, mpitype, mpierror)
  DEBUG_CHECKMPI(mpierror, "build_vartau_site_mpitype, MPI_Type_create_struct")
  call MPI_Type_commit(mpitype, mpierror)
  DEBUG_CHECKMPI(mpierror, "build_vartau_site_mpitype, MPI_Type_commit")

end subroutine build_vartau_site_mpitype
#endif

#ifdef LOCALBC
subroutine build_local_acccoef_site_mpitype(mpitype)
  implicit none

  integer, intent(out) :: mpitype ! type to be built

  integer, parameter :: count = 1 ! only local_acccoef

  type(lbe_site) :: sample
  integer(kind = MPI_ADDRESS_KIND) :: addrs(count), base, displs(count)
  integer :: blocklengths(count), types(count)
  integer :: mpierror

  ! Get base address of an lbe_site
  call MPI_Get_address(sample, base, mpierror)
  DEBUG_CHECKMPI(mpierror, "build_local_acccoef_site_mpitype, MPI_Get_address : base")

  ! local_acccoef
  blocklengths(1) = 1
  types(1) = LBE_REAL
  call MPI_Get_address(sample%local_acccoef,addrs(1), mpierror)
  DEBUG_CHECKMPI(mpierror, "build_local_acccoef_site_mpitype, MPI_Get_address : 1")

  ! Calculate displacements
  displs(:) = addrs(:) - base

  call MPI_Type_create_struct(count, blocklengths, displs, types, mpitype, mpierror)
  DEBUG_CHECKMPI(mpierror, "build_local_acccoef_site_mpitype, MPI_Type_create_struct")
  call MPI_Type_commit(mpitype, mpierror)
  DEBUG_CHECKMPI(mpierror, "build_local_acccoef_site_mpitype, MPI_Type_commit")

end subroutine build_local_acccoef_site_mpitype
#endif


!> generate displacement vector as needed by \c mpi_gatherv() or
!> \c mpi_scatterv() within \c comm_cart
!>
!> \param[in] counts number of elements for each rank
!>
!> \param[out] displs resulting displacement vector
subroutine calculate_displacements(counts,displs)
    integer,intent(in) :: counts(0:nprocs-1)
    integer,intent(out) :: displs(0:nprocs-1)
    integer i

    displs(0) = 0
    do i=0,nprocs-2
       displs(i+1) = displs(i) + counts(i)
    end do
end subroutine calculate_displacements

    !> root gathers  N(:,:,:)%rock_state  of the local lattice
    !> chunks  N(:,:,:)  from all processes into one global array
    !>  rs(:,:,:) . A halo of  halo_extent  lu is added to  rs .
    !> Thus  rs  must have been allocated with dimension
    !>  1-halo_extent:(/tnx,tny,tnz/)+halo_extent .
    subroutine gather_rock_state(N,rs)
        type(lbe_site),intent(in) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        real(kind=rk),intent(inout) :: &
             &rs(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        real(kind=rk),allocatable,dimension(:) :: sbuf,rbuf
        integer pcoords(3,0:nprocs-1) ! cartesian topology coordinates
        integer os(3)           ! lattice offset for different processors
        integer i,j,k,p,r,ierror,stat,he,h1

        allocate (sbuf(nx*ny*nz),stat=stat)
        call check_allocate(stat,'gather_rock_state(): sbuf')
        if (myrankc==0) then
           allocate (rbuf(tnx*tny*tnz),stat=stat)
           call check_allocate(stat,'gather_rock_state(): rbuf')
        end if

        ! fill  sbuf
        p = 0
        do i=1,nx
           do j=1,ny
              do k=1,nz
                 p = p + 1
                 sbuf(p) = N(i,j,k)%rock_state
              end do
           end do
        end do

        call mpi_gather(sbuf,nx*ny*nz,MPI_REAL8&
             &,rbuf,nx*ny*nz,MPI_REAL8,0,comm_cart,ierror)
        deallocate (sbuf)

        rank0: if (myrankc==0) then
           ! copy data from  rbuf  to  rs
           call find_topology(pcoords)
           p = 0
           do r=0,nprocs-1
              os(:) = pcoords(:,r)*(/nx,ny,nz/)
              do i=os(1)+1,os(1)+nx
                 do j=os(2)+1,os(2)+ny
                    do k=os(3)+1,os(3)+nz
                       p = p + 1
                       rs(i,j,k) = rbuf(p)
                    end do
                 end do
              end do
           end do
           deallocate (rbuf)

           ! create halo
           he = halo_extent
           h1 = halo_extent - 1

           rs(  -h1:     0,1:tny,1:tnz) = rs(tnx-h1:tnx,1:tny,1:tnz)
           rs(tnx+1:tnx+he,1:tny,1:tnz) = rs(     1: he,1:tny,1:tnz)

           rs(-h1:tnx+he,  -h1:     0,1:tnz) = rs(-h1:tnx+he,tny-h1:tny,1:tnz)
           rs(-h1:tnx+he,tny+1:tny+he,1:tnz) = rs(-h1:tnx+he,     1: he,1:tnz)

           rs(-h1:tnx+he,-h1:tny+he,-h1:0)=rs(-h1:tnx+he,-h1:tny+he,tnz-h1:tnz)
           rs(-h1:tnx+he,-h1:tny+he,tnz+1:tnz+he)=rs(-h1:tnx+he,-h1:tny+he,1:he)
        end if rank0
    end subroutine gather_rock_state

    !> finds the rank of the communication partner for a given
    !> dimension and (pseudo-)direction
    !>
    !> \param[in] cdim Cartesian dimension (1-3 representing x,y,z)
    !>
    !> \param[in] dir direction: {1,3}: negative, {2,4}: positive; <3:
    !> periodic, >=3: non-periodic (see explanation below)
    !>
    !> \param[in] leesedwards consider global x boundaries to be
    !> Lees-Edwards boundaries (optional, defaults to \c .false.)
    !>
    !> \returns rank in \c comm_cart, can be \c MPI_PROC_NULL if no
    !> send is required in the given direction
    !>
    !> The pseudo-directions 3 and 4 are used in the case of BC
    !> periodic_inflow to designate the special case of z-swaps within
    !> the non-periodic part of the system. In this case, directions 1
    !> and 2 are still used for swapping within the periodic
    !> sub-volume.
    integer function get_send_partner(cdim,dir,leesedwards)
        integer,intent(in) :: cdim,dir
        logical,intent(in),optional :: leesedwards
        integer fnc,n
        real(kind=rk) rlpz

        if (boundary/='periodic_inflow'.or.cdim<3) then
           if (dir>2.or.cdims(cdim)==1) then
              ! Don't send with dir>2 in this case; also don't
              ! exchange if only 1 box in a dimension, Lees-Edwards
              ! directions are set up anew in each (sub-)step
              get_send_partner = MPI_PROC_NULL
           else if (present_and_true(leesedwards).and.cdim==1.and.&
                &((ccoords(1)==0.and.dir==1)&
                &.or.(ccoords(1)==cdims(1)-1.and.dir==2))) then
              ! no conventional downward/upward send across
              ! Lees-Edwards plane
              get_send_partner = MPI_PROC_NULL
           else
              get_send_partner = nnprocs(cdim,dir)
           end if

!!$        else if (dir<3) then
!!$           rlpz = real(last_periodic_z,kind=rk)
!!$           ! dir is one of {1,2} --- z-swap with BC periodic_inflow
!!$           ! within periodic sub-volume
!!$
!!$           if (start(3)>last_periodic_z) then
!!$              get_send_partner = MPI_PROC_NULL
!!$           else if (last_periodic_z<=nz) then
!!$              ! don't exchange if only 1 box in a dimension
!!$              get_send_partner = MPI_PROC_NULL
!!$           else
!!$              if (dir==1) then
!!$                 if (ccoords(3)==0) then
!!$                    get_send_partner = &
!!$                         &owning_rankc((/start(1),start(2),last_periodic_z/))
!!$                 else
!!$                    get_send_partner = nnprocs(3,1)
!!$                 end if
!!$              else              ! dir==2 (positive z)
!!$                 if (border(2,3)>rlpz) then
!!$                    get_send_partner = owning_rankc((/start(1),start(2),1/))
!!$                 else
!!$                    get_send_partner = nnprocs(3,2)
!!$                 end if
!!$              end if
!!$           end if
!!$        else
!!$           ! dir is one of {3,4} --- z-swap with BC periodic_inflow
!!$           ! within non-periodic sub-volume
!!$           if (border(2,3)<enslavement_threshold) then
!!$              get_send_partner = MPI_PROC_NULL
!!$           else
!!$              if (dir==3) then  ! negative z
!!$                 if (border(1,3)<enslavement_threshold) then
!!$                    get_send_partner = MPI_PROC_NULL ! first non-periodic cpu
!!$                 else
!!$                    get_send_partner = nnprocs(3,1)
!!$                 end if
!!$              else              ! positive z (dir==4)
!!$                 if (border(2,3)>real(tnz,kind=rk)) then
!!$                    get_send_partner = MPI_PROC_NULL
!!$                 else
!!$                    get_send_partner = nnprocs(3,2)
!!$                 end if
!!$              end if
!!$           end if
        end if
    end function get_send_partner

!> finds the process that owns a certain lattice position
!>
!> \param[in] x integer lattice position in global coordinates (halo nodes
!> are not allowed here)
!>
!> \returns owning rank in Cartesian communicator
integer function owning_rankc(x)
    integer,intent(in) :: x(3)
    integer ierror,process_coordinates(3),rank

    process_coordinates = (x-1)/(/nx,ny,nz/)
    call MPI_Cart_rank(comm_cart,process_coordinates,rank,ierror)
    owning_rankc = rank
end function owning_rankc

!> Find ave/min/max and histogram of one or more real(kind=rk) values
!> from all processors
subroutine stats_rk(data,ave,xmax,xmin,histo,histotmp,nhisto)
    real(kind=rk),intent(in) :: data
    real(kind=rk),intent(out) :: ave,xmax,xmin
    integer,intent(inout) :: histo(*),histotmp(*)
    integer,intent(in) :: nhisto
    integer j,ierror
    real(kind=rk) aave,del

    call mpi_allreduce(data,aave,1,MPI_REAL8,MPI_SUM,comm_cart,ierror)
    ave = aave/nprocs
    call mpi_allreduce(data,xmax,1,MPI_REAL8,MPI_MAX,comm_cart,ierror)
    call mpi_allreduce(data,xmin,1,MPI_REAL8,MPI_MIN,comm_cart,ierror)

    histotmp(1:nhisto) = 0

    del = xmax-xmin
    if (del.eq.0.0) then
       j = 1
    else
       j = (data-xmin)/del * nhisto + 1
       if (j.gt.nhisto) j = nhisto
    endif
    histotmp(j) = histotmp(j) + 1

    call mpi_allreduce(histotmp,histo,nhisto,MPI_INTEGER,MPI_SUM,&
         &comm_cart,ierror)
end subroutine stats_rk

!>This routine was taken from ME3D and modified.
!>It initialises the MPI libraries, determines the rank
!>of each CPU, and sets up the CPU topology.
subroutine InitializeMPIcomms
  implicit none
  integer :: ierror, reslen ! reslen is just a dummy variable to be supplied as an argument

  ! Initialize MPI
  call MPI_Init(ierror)
  if (ierror .ne. MPI_SUCCESS) then
    call log_msg('MPI_Init() failed. Aborting...')
    call Abend
  endif

  ! Find out my rank
  call MPI_Comm_rank(MPI_COMM_WORLD, myrankw, ierror)
  if (myrankw == 0) then
    call log_msg_ws("Initialized MPI_COMM_WORLD, starting output.")
  endif

  ! Find out the total number of processors.
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierror)
  if (myrankw == 0) then
    write(msgstr,"('Requested ',i0,' processors.')") nprocs
    call log_msg(msgstr)
  endif
  ! Added by Elena to be included in metadata
  call MPI_Get_processor_name(hname, reslen, ierror)

  ! Temporarily set Comm_Cart to world so we don't have to change the broadcast interface
  ! to exchange data before ReorderComms is called. Also myrankc = myrankw for now.
  Comm_Cart = MPI_COMM_WORLD
  myrankc = myrankw

end subroutine InitializeMPIcomms

subroutine ReorderMPIcomms
  integer :: ierror,i


  call log_msg_hdr("Creating MPI Cartesian grid")
  ! Create the virtual topology.
  ! This call returns the dimensions of an appropriate
  ! Cartesian grid into which the processors may be divided.
  ! Would be nice to get a 'good' grid automatically
  CALL MPI_Dims_create(nprocs, nd, cdims, ierror)
  if (ierror .ne. MPI_SUCCESS) then
    CALL log_msg("MPI_Dims_create() failed. Aborting...")
    call Abend
  endif

  if (myrankw == 0) then
    CALL makedatestr(startsimul) ! For HDF5 output
    write (msgstr,"('Processors using a ',i0,' x ',i0,' x ',i0,' grid.')") cdims
    CALL log_msg(msgstr)
  endif

  ! Now create a communicator to make it easy for each processor
  ! to talk only to its immediate neighbours.
  call MPI_Cart_create(MPI_COMM_WORLD& ! Make from all processors.
                     &,nd&             ! Dimension of lattice
                     &,cdims&          ! Dimensions of lattice
                     &,periodic_p&     ! Whether or not periodic
                     &,.true.&         ! Reorder ranks?
                     &,comm_cart&      ! New Cartesian comunicator.
                     &,ierror)

  if (ierror .ne. MPI_SUCCESS) then
    CALL log_msg("MPI_Cart_create() failed. Aborting...")
    call Abend
  endif

  ! Find out my new rank:- it may change if MPI_Cart_create()
  ! reorders the rankings so that the virtual topology corresponds
  ! more accurately to the physical topology.
  CALL MPI_Comm_rank(MPI_COMM_WORLD, myrankw, ierror)
  CALL MPI_Comm_rank(Comm_Cart, myrankc, ierror)

  ! Now fill in the ccoords array with my position in the Cartesian lattice.
  CALL MPI_Cart_get( Comm_Cart, nd, cdims, periodic_p, ccoords, ierror)

  ! Determine who my nearest neighbours are.
  do i = 1, nd
    CALL MPI_Cart_shift( Comm_Cart, i-1, 1, nnprocs(i,1), nnprocs(i,2), ierror )
  end do

  if ( dbg_report_topology ) then
    write(msgstr,"('DEBUG: comm_cart = ',I0)") comm_cart
    call log_msg(msgstr)
    call debug_report_ccoords
    call debug_report_hostnames
  endif
end subroutine ReorderMPIcomms

!> logs MPI cartesian coordinates and neighbors
!>
!> \warning This is likely to fail for extreme core counts.
subroutine debug_report_ccoords
  integer i,ierror,recv_ccoords(nd),recv_nnprocs(2*nd),status(MPI_STATUS_SIZE)
  integer :: tag = 1 ! MPI tag

  ! Report on the coordinates.
  if ( myrankc .gt. 0 ) then
    CALL MPI_Send(ccoords, nd , MPI_INTEGER, 0, tag, Comm_cart, ierror)
  else
    call log_msg("debug_report_ccoords(): reporting ccoords for each rank.")
    write (msgstr&
         &,"('  Rank ',I6.6,' has coordinates (',i0, ',', i0, ',', i0, ').')") &
         &myrankc, ccoords(1), ccoords(2), ccoords(3)
    call log_msg(msgstr)
    do i=1,nprocs-1
      call mpi_recv(recv_ccoords,nd,MPI_INTEGER,i,tag,Comm_cart,status,ierror)
      write (msgstr&
           &,"('  Rank ',I6.6,' has coordinates (',i0, ',', i0, ',', i0, ').')"&
           &) i, recv_ccoords(1), recv_ccoords(2), recv_ccoords(3)
      call log_msg(msgstr)
    enddo
  endif

  ! Report on nearest neighbours
  if ( myrankc .gt. 0 ) then
    CALL MPI_Send(nnprocs, 2*nd , MPI_INTEGER, 0, tag, Comm_cart, ierror)
  else
    call log_msg('debug_report_ccoords(): reporting nearest neighbours for '&
         &//'each rank (cx-1 cx+1 cy-1 cy+1 cz-1 cz+1).')
    write (msgstr,"('  Rank ',I6.6,' has neighbours ',6(I6.6,:,' '),'.')") &
         &myrankc,nnprocs(1,1:2),nnprocs(2,1:2),nnprocs(3,1:2)
    call log_msg(msgstr)
    do i=1, nprocs - 1
      call mpi_recv(recv_nnprocs,2*nd,MPI_INTEGER,i,tag,Comm_cart,status,ierror)
      write (msgstr,"('  Rank ',I6.6,' has neighbours ',6(I6.6,:,' '),'.')") &
           &i,recv_nnprocs(1),recv_nnprocs(4),recv_nnprocs(2),recv_nnprocs(5)&
           &,recv_nnprocs(3),recv_nnprocs(6)
      call log_msg(msgstr)
    enddo
  endif
end subroutine debug_report_ccoords

!> logs hostnames
!>
!> \warning This is likely to fail for extreme core counts.
subroutine debug_report_hostnames
  integer i,ierror,status(MPI_STATUS_SIZE)
  integer :: tag = 1 ! MPI tag
  character(len=MPI_MAX_PROCESSOR_NAME) :: recv_hname

  if ( myrankc .gt. 0 ) then
    call mpi_send(hname,MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER,0,tag,Comm_cart&
         &,ierror)
  else
    call log_msg("debug_report_hostnames(): reporting hostnames for each rank")
    write (msgstr,"('  Rank ',I6.6,' has been assigned to host <',A,'>.')") &
         &myrankc,trim(hname)
    call log_msg(msgstr)
    do i=1, nprocs-1
      call mpi_recv(recv_hname,MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER,i,tag&
           &,Comm_cart,status,ierror)
      write (msgstr,"('  Rank ',I6.6,' has been assigned to host <',A,'>.')") &
           &i,trim(recv_hname)
      call log_msg(msgstr)
    enddo
  endif
end subroutine debug_report_hostnames

!>  Divides the tnx x tny x tnz grid into nx x ny x nz blocks assigned
!>  to the processors
subroutine lbe_divide_grid()
  ! Must have the nx, ny, nz integer divisible by the
  ! number of processors assigned in that cartesian direction.

  if ( (MOD(nx,cdims(1))/= 0) .or. &
       (MOD(ny,cdims(2))/= 0) .or. &
       (MOD(nz,cdims(3))/= 0) ) then
    call log_msg("FATAL ERROR: Grid does not divide evenly by the number of processors. Aborting...")
    call Abend
  endif

  ! Remember the size of the original domain.
  tnx = nx; tny = ny; tnz = nz

  ! Now divide these digits so that a processor will only ever
  ! assign the amount of memory required for its own data.
  nx = nx/cdims(1)
  ny = ny/cdims(2)
  nz = nz/cdims(3)

  write(msgstr,"('Grid decomposition in blocks of size: ',i0,' x ',i0,' x ',i0,'.')") nx, ny, nz
  call log_msg(msgstr)

  ! Find out where the local data lies in relation to the global data
  ! set. The array `start' will identify the start of the local data
  ! relative to the global data set.
  start(1) = ccoords(1)*nx + 1
  start(2) = ccoords(2)*ny + 1
  start(3) = ccoords(3)*nz + 1
end subroutine lbe_divide_grid

!> Reports on calculation performance at the end of a simulation.
subroutine calculate_performance(delta_t)
  integer*8              :: nsites, nupdates
  integer                :: delta_t
  real*8                 :: updaterate, ratepercpu

  nsites = int(tnx,kind=8)*int(tny,kind=8)*int(tnz,kind=8)
  nupdates = timesteps_count * nsites
  updaterate = real(nupdates) / real(delta_t)
  ratepercpu = updaterate / real(nprocs)

  write(msgstr,"('Updated ',i0,' lattice sites for ',i0,' time steps -> ',i0,' site updates.')") &
    nsites, timesteps_count, nupdates
  call log_msg(msgstr)
  write(msgstr,"('This took ',i0,' seconds ->')") delta_t
  call log_msg(msgstr)
  write(msgstr,"('    ',F16.2,' updates per second')") updaterate
  call log_msg(msgstr)
  write(msgstr,"('    ',F16.2,' updates per second per CPU')") ratepercpu
  call log_msg(msgstr)

end subroutine calculate_performance

!>This routine is called just before program termination - it
!>tells the MPI libraries to clean up and finish.
subroutine FinalizeMPIcomms
  implicit none
  integer :: ierror

  call log_msg_hdr("Calling MPI_Finalize()")
  CALL MPI_Finalize(ierror)
  if (myrankc == 0) then
    CALL calculate_performance(end_time_t - start_time_t)
  endif
end subroutine FinalizeMPIcomms

!> initialize stuff related to communication
subroutine lbe_parallel_init()
  tsize = real((/ tnx,tny,tnz /),kind=rk)
  chunksize = real((/ nx,ny,nz /),kind=rk)

  maxpos = tsize + 0.5_rk
  minpos = 0.5_rk

  border(1,:) = start - 0.5_rk
  border(2,:) = start + chunksize - 0.5_rk
end subroutine lbe_parallel_init

!> sets up data types used by \c lbe_halo_exchange_new() ---also read
!> comment there
subroutine lbe_setup_halo_exchange_new(N)
    type(lbe_site),intent(inout) :: &
         &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
    integer :: he,h1,i,lbe_site_mpitype

    ! just abbreviations
    he = halo_extent
    h1 = halo_extent - 1

    ! setup custom mpi datatypes that simplify the halo exchange vastly

    call build_lbe_site_mpitype((/(.true.,i=1,19)/)&
#ifndef NOSURFACTANT
         &,.true.&
#endif
         &,.true.,lbe_site_mpitype)

    !! This type setup causes a full halo exchange, so lbe_halo_exchange()
    !! may be replaced by lbe_halo_exchange_new() .

    ! x swaps (still restricted y/z-range: data outside is out-dated)
    call build_lattice_chunk_mpitype&
         &(N,(/1,      he/),(/1,      ny/),(/1,      nz/)&
         &,lbe_site_mpitype,ls_mpitype(1))
    call build_lattice_chunk_mpitype&
         &(N,(/nx-h1,  nx/),(/1,      ny/),(/1,      nz/)&
         &,lbe_site_mpitype,us_mpitype(1))
    call build_lattice_chunk_mpitype&
         &(N,(/-h1,     0/),(/1,      ny/),(/1,      nz/)&
         &,lbe_site_mpitype,lr_mpitype(1))
    call build_lattice_chunk_mpitype&
         &(N,(/nx+1,nx+he/),(/1,      ny/),(/1,      nz/)&
         &,lbe_site_mpitype,ur_mpitype(1))

    ! y swaps (full x-range - up-to-date data was received in x-swaps)
    call build_lattice_chunk_mpitype&
         &(N,(/-h1, nx+he/),(/1,      he/),(/1,      nz/)&
         &,lbe_site_mpitype,ls_mpitype(2))
    call build_lattice_chunk_mpitype&
         &(N,(/-h1, nx+he/),(/ny-h1,  ny/),(/1,      nz/)&
         &,lbe_site_mpitype,us_mpitype(2))
    call build_lattice_chunk_mpitype&
         &(N,(/-h1, nx+he/),(/-h1,     0/),(/1,      nz/)&
         &,lbe_site_mpitype,lr_mpitype(2))
    call build_lattice_chunk_mpitype&
         &(N,(/-h1, nx+he/),(/ny+1,ny+he/),(/1,      nz/)&
         &,lbe_site_mpitype,ur_mpitype(2))

    ! z swaps (even full y-range - after z swap all data is complete)
    call build_lattice_chunk_mpitype&
         &(N,(/-h1, nx+he/),(/-h1, ny+he/),(/1,      he/)&
         &,lbe_site_mpitype,ls_mpitype(3))
    call build_lattice_chunk_mpitype&
         &(N,(/-h1, nx+he/),(/-h1, ny+he/),(/nz-h1,  nz/)&
         &,lbe_site_mpitype,us_mpitype(3))
    call build_lattice_chunk_mpitype&
         &(N,(/-h1, nx+he/),(/-h1, ny+he/),(/-h1,     0/)&
         &,lbe_site_mpitype,lr_mpitype(3))
    call build_lattice_chunk_mpitype&
         &(N,(/-h1, nx+he/),(/-h1, ny+he/),(/nz+1,nz+he/)&
         &,lbe_site_mpitype,ur_mpitype(3))


    !! This type setup causes lbe_halo_exchange_new() to exchange only
    !! that part of the lattice that is not exchanged by
    !! lbe_halo_exchange() already.

!    ! x swaps (still restricted y/z-range: data outside is out-dated)
!    call build_lattice_chunk_mpitype&
!         &(N,(/2,      he/),(/0,    ny+1/),(/0,    nz+1/)&
!         &,lbe_site_mpitype,ls_mpitype(1))
!    call build_lattice_chunk_mpitype&
!         &(N,(/nx-h1,nx-1/),(/0,    ny+1/),(/0,    nz+1/)&
!         &,lbe_site_mpitype,us_mpitype(1))
!    call build_lattice_chunk_mpitype&
!         &(N,(/-h1,    -1/),(/0,    ny+1/),(/0,    nz+1/)&
!         &,lbe_site_mpitype,lr_mpitype(1))
!    call build_lattice_chunk_mpitype&
!         &(N,(/nx+2,nx+he/),(/0,    ny+1/),(/0,    nz+1/)&
!         &,lbe_site_mpitype,ur_mpitype(1))
!
!    ! y swaps (full x-range - up-to-date data was received in x-swaps)
!    call build_lattice_chunk_mpitype&
!         &(N,(/-h1, nx+he/),(/2,      he/),(/0,    nz+1/)&
!         &,lbe_site_mpitype,ls_mpitype(2))
!    call build_lattice_chunk_mpitype&
!         &(N,(/-h1, nx+he/),(/ny-h1,ny-1/),(/0,    nz+1/)&
!         &,lbe_site_mpitype,us_mpitype(2))
!    call build_lattice_chunk_mpitype&
!         &(N,(/-h1, nx+he/),(/-h1,    -1/),(/0,    nz+1/)&
!         &,lbe_site_mpitype,lr_mpitype(2))
!    call build_lattice_chunk_mpitype&
!         &(N,(/-h1, nx+he/),(/ny+2,ny+he/),(/0,    nz+1/)&
!         &,lbe_site_mpitype,ur_mpitype(2))
!
!    ! z swaps (even full y-range - after z swap all data is complete)
!    call build_lattice_chunk_mpitype&
!         &(N,(/-h1, nx+he/),(/-h1, ny+he/),(/2,      he/)&
!         &,lbe_site_mpitype,ls_mpitype(3))
!    call build_lattice_chunk_mpitype&
!         &(N,(/-h1, nx+he/),(/-h1, ny+he/),(/nz-h1,nz-1/)&
!         &,lbe_site_mpitype,us_mpitype(3))
!    call build_lattice_chunk_mpitype&
!         &(N,(/-h1, nx+he/),(/-h1, ny+he/),(/-h1,    -1/)&
!         &,lbe_site_mpitype,lr_mpitype(3))
!    call build_lattice_chunk_mpitype&
!         &(N,(/-h1, nx+he/),(/-h1, ny+he/),(/nz+2,nz+he/)&
!         &,lbe_site_mpitype,ur_mpitype(3))
end subroutine lbe_setup_halo_exchange_new

!> yet another halo exchange routine
!>
!> This one was previously called \c md_halo_exchange() and it was
!> located in \c lbe_md_fluid_module. It was renamed and moved here
!> since it is required also by other modules outside \c lbe_md* . At
!> the moment, it is used exclusively to communicate data in the outer
!> halo layers (everything beyond the 1st layer that is already
!> exchanged by \c lbe_halo_exchange() ). Therefore the mpi types
!> hardcoded here are, at the moment, defined in a way that skips the
!> 1st halo layer. This could be easily changed. This could be
!> probably replaced by calling Stefan's more general \c
!> halo_exchange() routine with a newly defined halo type
!> instance. This should be done some time!
subroutine lbe_halo_exchange_new(N)
    type(lbe_site),intent(inout) :: &
         &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
    integer k,ierror
    integer status(MPI_STATUS_SIZE)

    do k = 1,3
       call mpi_sendrecv&   ! send "downward"
            &(N,1,ls_mpitype(k),nnprocs(k,1),0&
            &,N,1,ur_mpitype(k),nnprocs(k,2),0&
            &,comm_cart,status,ierror)
       call mpi_sendrecv&   ! send "upward"
            &(N,1,us_mpitype(k),nnprocs(k,2),0&
            &,N,1,lr_mpitype(k),nnprocs(k,1),0&
            &,comm_cart,status,ierror)
    end do
end subroutine lbe_halo_exchange_new

!!$!>Performs the halo exchange.
!!$!>
!!$!> WARNING: RELTIME is obsolete (changed to VARTAU); features have to be checked and maintained before use.
!!$!>Naming convention: there are two phases to the exchange in each dimension.
!!$!>When exchanging in the X direction, first data is propagated in the +x
!!$!>direction, then in the -x direction. "p" refers to *data* being propagated
!!$!>in the +ve direction, "m" refers to data propagating in the -ve direction;
!!$!>hence:
!!$!>pinbuf = buffer for data which will have travelled in the +ve direction
!!$!>poutbuf = buffer for data which will travel  in the +ve direction.
!!$!>s/p/m/ for -ve direction, etc.
!!$subroutine lbe_halo_exchange(N)
!!$	type(lbe_site),dimension(0:,0:,0:) :: N
!!$
!!$	real*8,dimension(:,:,:), allocatable :: pinbuf
!!$	real*8,dimension(:,:,:), allocatable :: minbuf
!!$	real*8,dimension(:,:,:), allocatable :: poutbuf
!!$	real*8,dimension(:,:,:), allocatable :: moutbuf
!!$
!!$	integer :: ibuf
!!$
!!$	integer, parameter :: pin=1,min=2,pout=3,mout=4
!!$
!!$	integer, dimension(4) :: requests
!!$	integer, dimension(MPI_STATUS_SIZE,4) :: statuses
!!$	integer :: ierror
!!$	integer :: i,j
!!$	integer :: added
!!$	integer :: x,y,z
!!$        character(LEN=80) :: bla
!!$
!!$        logical :: minz,maxz
!!$
!!$added =0
!!$#ifndef NOSURFACTANT
!!$        ibuf = 61
!!$#else
!!$        ibuf = 39
!!$#endif
!!$#ifdef SINGLEFLUID
!!$	ibuf = 20
!!$#endif
!!$
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$#ifndef RELTIME
!!$#ifdef DIST
!!$added = 1
!!$#ifndef NOSURFACTANT
!!$        ibuf = 62
!!$#else
!!$        ibuf = 40
!!$#endif
!!$#ifdef SINGLEFLUID
!!$	ibuf = 21
!!$#endif
!!$#endif
!!$#endif
!!$
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$#ifndef DIST
!!$#ifdef RELTIME
!!$#ifndef NOSURFACTANT
!!$        ibuf = 64
!!$	added = 3
!!$#else
!!$        ibuf = 41
!!$	added = 2
!!$#endif
!!$#ifdef SINGLEFLUID
!!$	ibuf = 21
!!$	added =1
!!$#endif
!!$#endif
!!$#endif
!!$
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$#ifdef DIST
!!$#ifdef RELTIME
!!$#ifndef NOSURFACTANT
!!$        ibuf = 65
!!$	added = 4
!!$#else
!!$        ibuf = 42
!!$	added = 3
!!$#endif
!!$#ifdef SINGLEFLUID
!!$	ibuf = 22
!!$	added = 2
!!$#endif
!!$#endif
!!$#endif
!!$
!!$#ifdef LOCALBC
!!$	ibuf = ibuf+1
!!$	added = added +1
!!$#endif
!!$
!!$	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$	! Swap in the X direction
!!$	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$! Elena: change array order of pinbuf, poutbuf, minbuf, moutbuf
!!$! Still, not all arrays will be contiguous in memory...
!!$
!!$        allocate(pinbuf(ibuf,0:ny+1,0:nz+1))
!!$        allocate(poutbuf(ibuf,0:ny+1,0:nz+1))
!!$        allocate(minbuf(ibuf,0:ny+1,0:nz+1))
!!$        allocate(moutbuf(ibuf,0:ny+1,0:nz+1))
!!$
!!$	!!!!!!
!!$	! Begin asynchronous receives
!!$	!!!!!!
!!$
!!$	call MPI_Irecv(		pinbuf,		& ! buf
!!$				(ny+2)*(nz+2)*ibuf,	& ! count
!!$				LBE_REAL,	& ! datatype
!!$				nnprocs(1,1),	& ! source
!!$				tag_px,		& ! tag
!!$				Comm_Cart,	& ! communicator
!!$				requests(pin),	& ! request
!!$				ierror)
!!$
!!$        bla = 'Bad +x async receive'
!!$	call checkmpi(ierror,bla)
!!$
!!$	call MPI_Irecv(		minbuf,		& ! buf
!!$				(ny+2)*(nz+2)*ibuf,	& ! count
!!$				LBE_REAL,	& ! datatype
!!$				nnprocs(1,2),	& ! source
!!$				tag_mx,		& ! tag
!!$				Comm_Cart,	& ! communicator
!!$				requests(min),	& ! request
!!$				ierror)
!!$
!!$        bla = 'Bad -x async receive'
!!$	call checkmpi(ierror,bla)
!!$
!!$	do j=0,nz+1
!!$	 do i=0,ny+1
!!$
!!$
!!$               poutbuf(1:19,i,j) = N(nx,i,j)%n_r(:)
!!$#ifndef SINGLEFLUID
!!$               poutbuf(20:38,i,j) = N(nx,i,j)%n_b(:)
!!$#endif
!!$#ifndef NOSURFACTANT
!!$               poutbuf(39:57,i,j) = N(nx,i,j)%n_s(:)
!!$               poutbuf(58:60,i,j) = N(nx,i,j)%d(:)
!!$#endif
!!$#ifdef LOCALBC
!!$		poutbuf(ibuf-added,i,j) = N(nx,i,j)%local_acccoef
!!$#endif
!!$
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$
!!$#ifndef DIST
!!$#ifndef RELTIME
!!$               poutbuf(ibuf,i,j) = N(nx,i,j)%rock_state
!!$#endif
!!$#endif
!!$
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$#ifndef RELTIME
!!$#ifdef DIST
!!$		poutbuf(ibuf-1,i,j) = N(nx,i,j)%rock_state
!!$		poutbuf(ibuf,i,j) = N(nx,i,j)%abst
!!$#endif
!!$#endif
!!$
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$
!!$#ifndef DIST
!!$#ifdef RELTIME
!!$		poutbuf(ibuf-1,i,j) = N(nx,i,j)%rock_state
!!$		poutbuf(ibuf,i,j) = N(nx,i,j)%taupos_r
!!$#ifndef SINGLEFLUID
!!$		poutbuf(ibuf-2,i,j) = N(nx,i,j)%taupos_b
!!$#endif
!!$#ifndef NOSURFACTANT
!!$		poutbuf(ibuf-2,i,j) = N(nx,i,j)%taupos_b
!!$		poutbuf(ibuf-3,i,j) = N(nx,i,j)%taupos_s
!!$#endif
!!$
!!$#endif
!!$#endif
!!$
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$
!!$#ifdef RELTIME
!!$#ifdef DIST
!!$		poutbuf(ibuf-2,i,j) = N(nx,i,j)%rock_state
!!$		poutbuf(ibuf-1,i,j) = N(nx,i,j)%abst
!!$		poutbuf(ibuf,i,j) = N(nx,i,j)%taupos_r
!!$
!!$#ifndef SINGLEFLUID
!!$		poutbuf(ibuf-3,i,j) = N(nx,i,j)%taupos_b
!!$#endif
!!$#ifndef NOSURFACTANT
!!$		poutbuf(ibuf-3,i,j) = N(nx,i,j)%taupos_b
!!$		poutbuf(ibuf-4,i,j) = N(nx,i,j)%taupos_s
!!$#endif
!!$#endif
!!$#endif
!!$
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!
!!$!
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$
!!$               moutbuf(1:19,i,j) = N(1,i,j)%n_r(:)
!!$#ifndef SINGLEFLUID
!!$               moutbuf(20:38,i,j) = N(1,i,j)%n_b(:)
!!$#endif
!!$#ifndef NOSURFACTANT
!!$               moutbuf(39:57,i,j) = N(1,i,j)%n_s(:)
!!$               moutbuf(58:60,i,j) = N(1,i,j)%d(:)
!!$#endif
!!$
!!$
!!$#ifdef LOCALBC
!!$		moutbuf(ibuf-added,i,j) = N(1,i,j)%local_acccoef
!!$#endif
!!$
!!$
!!$#ifndef DIST
!!$#ifndef RELTIME
!!$               moutbuf(ibuf,i,j) = N(1,i,j)%rock_state
!!$#endif
!!$#endif
!!$
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$#ifdef DIST
!!$#ifndef RELTIME
!!$		moutbuf(ibuf-1,i,j) = N(1,i,j)%rock_state
!!$		moutbuf(ibuf,i,j) = N(1,i,j)%abst
!!$#endif
!!$#endif
!!$
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$#ifndef DIST
!!$#ifdef RELTIME
!!$		moutbuf(ibuf-1,i,j) = N(1,i,j)%rock_state
!!$		moutbuf(ibuf,i,j) = N(1,i,j)%taupos_r
!!$
!!$#ifndef SINGLEFLUID
!!$		moutbuf(ibuf-2,i,j) = N(1,i,j)%taupos_b
!!$#endif
!!$#ifndef NOSURFACTANT
!!$		moutbuf(ibuf-2,i,j) = N(1,i,j)%taupos_b
!!$		moutbuf(ibuf-3,i,j) = N(1,i,j)%taupos_s
!!$#endif
!!$#endif
!!$#endif
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$#ifdef RELTIME
!!$#ifdef DIST
!!$		moutbuf(ibuf-2,i,j) = N(1,i,j)%rock_state
!!$		moutbuf(ibuf-1,i,j) = N(1,i,j)%abst
!!$		moutbuf(ibuf,i,j) = N(1,i,j)%taupos_r
!!$
!!$#ifndef SINGLEFLUID
!!$		moutbuf(ibuf-3,i,j) = N(1,i,j)%taupos_b
!!$#endif
!!$#ifndef NOSURFACTANT
!!$		moutbuf(ibuf-3,i,j) = N(1,i,j)%taupos_b
!!$		moutbuf(ibuf-4,i,j) = N(1,i,j)%taupos_s
!!$#endif
!!$#endif
!!$#endif
!!$
!!$	 end do
!!$	end do
!!$
!!$
!!$
!!$	!!!!!!
!!$	! Begin asynchronous sends
!!$	!!!!!!
!!$
!!$	call MPI_ISend(		poutbuf,	& ! buf
!!$				(ny+2)*(nz+2)*ibuf,	& ! count
!!$				LBE_REAL,	& ! datatype
!!$				nnprocs(1,2),	& ! dest
!!$				tag_px,		& ! tag
!!$				Comm_Cart,	& ! communicator
!!$				requests(pout),	& ! request
!!$				ierror)
!!$
!!$        bla = 'Bad +x async send'
!!$	call checkmpi(ierror,bla)
!!$
!!$	call MPI_ISend(		moutbuf,	& ! buf
!!$				(ny+2)*(nz+2)*ibuf,	& ! count
!!$				LBE_REAL,	& ! datatype
!!$				nnprocs(1,1),	& ! dest
!!$				tag_mx,		& ! tag
!!$				Comm_Cart,	& ! communicator
!!$				requests(mout),	& ! request
!!$				ierror)
!!$
!!$        bla = 'Bad -x async send'
!!$	call checkmpi(ierror,bla)
!!$
!!$	!!!!!!
!!$	! Now wait for all I/O to complete
!!$	!!!!!!
!!$
!!$	call MPI_Waitall(4,requests,statuses,ierror)
!!$
!!$        bla = 'MPI_Waitall() failed in X direction'
!!$	call checkmpi(ierror,bla)
!!$
!!$        do j=0,nz+1
!!$         do i=0,ny+1
!!$                 N(0,i,j)%n_r(:) = pinbuf(1:19,i,j)
!!$#ifndef SINGLEFLUID
!!$                 N(0,i,j)%n_b(:) = pinbuf(20:38,i,j)
!!$#endif
!!$#ifndef NOSURFACTANT
!!$                 N(0,i,j)%n_s(:) = pinbuf(39:57,i,j)
!!$                 N(0,i,j)%d(:)   = pinbuf(58:60,i,j)
!!$#endif
!!$
!!$
!!$#ifdef LOCALBC
!!$		N(0,i,j)%local_acccoef = pinbuf(ibuf-added,i,j)
!!$#endif
!!$
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$#ifndef DIST
!!$#ifndef RELTIME
!!$               N(0,i,j)%rock_state = pinbuf(ibuf,i,j)
!!$#endif
!!$#endif
!!$
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$#ifndef RELTIME
!!$#ifdef DIST
!!$		N(0,i,j)%rock_state = pinbuf(ibuf-1,i,j)
!!$		N(0,i,j)%abst = pinbuf(ibuf,i,j)
!!$#endif
!!$#endif
!!$
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$#ifndef DIST
!!$#ifdef RELTIME
!!$		N(0,i,j)%rock_state = pinbuf(ibuf-1,i,j)
!!$		N(0,i,j)%taupos_r = pinbuf(ibuf,i,j)
!!$
!!$#ifndef SINGLEFLUID
!!$		N(0,i,j)%taupos_b = pinbuf(ibuf-2,i,j)
!!$#endif
!!$#ifndef NOSURFACTANT
!!$		N(0,i,j)%taupos_b = pinbuf(ibuf-2,i,j)
!!$		N(0,i,j)%taupos_s = pinbuf(ibuf-3,i,j)
!!$#endif
!!$#endif
!!$#endif
!!$
!!$
!!$#ifdef RELTIME
!!$#ifdef DIST
!!$		N(0,i,j)%rock_state = pinbuf(ibuf-2,i,j)
!!$		N(0,i,j)%abst = pinbuf(ibuf-1,i,j)
!!$		N(0,i,j)%taupos_r = pinbuf(ibuf,i,j)
!!$
!!$#ifndef SINGLEFLUID
!!$		N(0,i,j)%taupos_b = pinbuf(ibuf-3,i,j)
!!$#endif
!!$#ifndef NOSURFACTANT
!!$		N(0,i,j)%taupos_b = pinbuf(ibuf-3,i,j)
!!$		N(0,i,j)%taupos_s = pinbuf(ibuf-4,i,j)
!!$#endif
!!$#endif
!!$#endif
!!$
!!$                 N(nx+1,i,j)%n_r(:) = minbuf(1:19,i,j)
!!$#ifndef SINGLEFLUID
!!$                 N(nx+1,i,j)%n_b(:) = minbuf(20:38,i,j)
!!$#endif
!!$#ifndef NOSURFACTANT
!!$                 N(nx+1,i,j)%n_s(:) = minbuf(39:57,i,j)
!!$                 N(nx+1,i,j)%d(:)   = minbuf(58:60,i,j)
!!$#endif
!!$
!!$#ifdef LOCALBC
!!$		N(nx+1,i,j)%local_acccoef = minbuf(ibuf-added,i,j)
!!$#endif
!!$
!!$#ifndef DIST
!!$#ifndef RELTIME
!!$               N(nx+1,i,j)%rock_state = minbuf(ibuf,i,j)
!!$#endif
!!$#endif
!!$
!!$#ifndef RELTIME
!!$#ifdef DIST
!!$		N(nx+1,i,j)%rock_state = minbuf(ibuf-1,i,j)
!!$		N(nx+1,i,j)%abst = minbuf(ibuf,i,j)
!!$#endif
!!$#endif
!!$
!!$#ifndef DIST
!!$#ifdef RELTIME
!!$		N(nx+1,i,j)%rock_state = minbuf(ibuf-1,i,j)
!!$		N(nx+1,i,j)%taupos_r = minbuf(ibuf,i,j)
!!$
!!$#ifndef SINGLEFLUID
!!$		N(nx+1,i,j)%taupos_b = minbuf(ibuf-2,i,j)
!!$#endif
!!$#ifndef NOSURFACTANT
!!$		N(nx+1,i,j)%taupos_b = minbuf(ibuf-2,i,j)
!!$		N(nx+1,i,j)%taupos_s = minbuf(ibuf-3,i,j)
!!$#endif
!!$#endif
!!$#endif
!!$
!!$#ifdef RELTIME
!!$#ifdef DIST
!!$		N(nx+1,i,j)%rock_state = minbuf(ibuf-2,i,j)
!!$		N(nx+1,i,j)%abst = minbuf(ibuf-1,i,j)
!!$		N(nx+1,i,j)%taupos_r = minbuf(ibuf,i,j)
!!$
!!$#ifndef SINGLEFLUID
!!$		N(nx+1,i,j)%taupos_b = minbuf(ibuf-3,i,j)
!!$#endif
!!$#ifndef NOSURFACTANT
!!$		N(nx+1,i,j)%taupos_b = minbuf(ibuf-3,i,j)
!!$		N(nx+1,i,j)%taupos_s = minbuf(ibuf-4,i,j)
!!$#endif
!!$#endif
!!$#endif
!!$         end do
!!$        end do
!!$
!!$	deallocate(pinbuf)
!!$	deallocate(poutbuf)
!!$	deallocate(minbuf)
!!$	deallocate(moutbuf)
!!$
!!$	!!!!
!!$	! X swap done.
!!$	!!!!
!!$
!!$	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$	! Swap in the Y direction
!!$	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$       allocate(pinbuf(ibuf,0:nx+1,0:nz+1))
!!$       allocate(poutbuf(ibuf,0:nx+1,0:nz+1))
!!$       allocate(moutbuf(ibuf,0:nx+1,0:nz+1))
!!$       allocate(minbuf(ibuf,0:nx+1,0:nz+1))
!!$
!!$	!!!!!!
!!$	! Begin asynchronous receives
!!$	!!!!!!
!!$
!!$	call MPI_Irecv(		pinbuf,		& ! buf
!!$				(nx+2)*(nz+2)*ibuf,	& ! count
!!$				LBE_REAL,	& ! datatype
!!$				nnprocs(2,1),	& ! source
!!$				tag_py,		& ! tag
!!$				Comm_Cart,	& ! communicator
!!$				requests(pin),	& ! request
!!$				ierror)
!!$        bla = 'Bad +y async receive'
!!$	call checkmpi(ierror,bla)
!!$
!!$	call MPI_Irecv(		minbuf,		& ! buf
!!$				(nx+2)*(nz+2)*ibuf,	& ! count
!!$				LBE_REAL,	& ! datatype
!!$				nnprocs(2,2),	& ! source
!!$				tag_my,		& ! tag
!!$				Comm_Cart,	& ! communicator
!!$				requests(min),	& ! request
!!$				ierror)
!!$        bla = 'Bad -y async receive'
!!$	call checkmpi(ierror,bla)
!!$
!!$        do j=0,nz+1
!!$         do i=0,nx+1
!!$                 poutbuf(1:19,i,j) = N(i,ny,j)%n_r(:)
!!$#ifndef SINGLEFLUID
!!$                 poutbuf(20:38,i,j) = N(i,ny,j)%n_b(:)
!!$#endif
!!$#ifndef NOSURFACTANT
!!$                 poutbuf(39:57,i,j) = N(i,ny,j)%n_s(:)
!!$                 poutbuf(58:60,i,j) = N(i,ny,j)%d(:)
!!$#endif
!!$
!!$#ifdef LOCALBC
!!$		 poutbuf(ibuf-added,i,j) = N(i,ny,j)%local_acccoef
!!$#endif
!!$
!!$
!!$#ifndef DIST
!!$#ifndef RELTIME
!!$               poutbuf(ibuf,i,j) = N(i,ny,j)%rock_state
!!$#endif
!!$#endif
!!$
!!$
!!$#ifndef RELTIME
!!$#ifdef DIST
!!$		poutbuf(ibuf-1,i,j) = N(i,ny,j)%rock_state
!!$		poutbuf(ibuf,i,j) = N(i,ny,j)%abst
!!$#endif
!!$#endif
!!$
!!$#ifndef DIST
!!$#ifdef RELTIME
!!$		poutbuf(ibuf-1,i,j) = N(i,ny,j)%rock_state
!!$		poutbuf(ibuf,i,j) = N(i,ny,j)%taupos_r
!!$
!!$#ifndef SINGLEFLUID
!!$		poutbuf(ibuf-2,i,j) = N(i,ny,j)%taupos_b
!!$#endif
!!$#ifndef NOSURFACTANT
!!$		poutbuf(ibuf-2,i,j) = N(i,ny,j)%taupos_b
!!$		poutbuf(ibuf-3,i,j) = N(i,ny,j)%taupos_s
!!$#endif
!!$#endif
!!$#endif
!!$
!!$#ifdef RELTIME
!!$#ifdef DIST
!!$		poutbuf(ibuf-2,i,j) = N(i,ny,j)%rock_state
!!$		poutbuf(ibuf-1,i,j) = N(i,ny,j)%abst
!!$		poutbuf(ibuf,i,j) = N(i,ny,j)%taupos_r
!!$
!!$#ifndef SINGLEFLUID
!!$		poutbuf(ibuf-3,i,j) = N(i,ny,j)%taupos_b
!!$#endif
!!$#ifndef NOSURFACTANT
!!$		poutbuf(ibuf-3,i,j) = N(i,ny,j)%taupos_b
!!$		poutbuf(ibuf-4,i,j) = N(i,ny,j)%taupos_s
!!$#endif
!!$#endif
!!$#endif
!!$
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$                 moutbuf(1:19,i,j) = N(i,1,j)%n_r(:)
!!$#ifndef SINGLEFLUID
!!$                 moutbuf(20:38,i,j) = N(i,1,j)%n_b(:)
!!$#endif
!!$#ifndef NOSURFACTANT
!!$                 moutbuf(39:57,i,j) = N(i,1,j)%n_s(:)
!!$                 moutbuf(58:60,i,j) = N(i,1,j)%d(:)
!!$#endif
!!$
!!$#ifdef LOCALBC
!!$		 moutbuf(ibuf-added,i,j) = N(i,1,j)%local_acccoef
!!$#endif
!!$
!!$#ifndef DIST
!!$#ifndef RELTIME
!!$               moutbuf(ibuf,i,j) = N(i,1,j)%rock_state
!!$#endif
!!$#endif
!!$
!!$#ifndef RELTIME
!!$#ifdef DIST
!!$		moutbuf(ibuf-1,i,j) = N(i,1,j)%rock_state
!!$		moutbuf(ibuf,i,j) = N(i,1,j)%abst
!!$#endif
!!$#endif
!!$
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$#ifndef DIST
!!$#ifdef RELTIME
!!$		moutbuf(ibuf-1,i,j) = N(i,1,j)%rock_state
!!$		moutbuf(ibuf,i,j) = N(i,1,j)%taupos_r
!!$#ifndef SINGLEFLUID
!!$		moutbuf(ibuf-2,i,j) = N(i,1,j)%taupos_b
!!$#endif
!!$#ifndef NOSURFACTANT
!!$		moutbuf(ibuf-2,i,j) = N(i,1,j)%taupos_b
!!$		moutbuf(ibuf-3,i,j) = N(i,1,j)%taupos_s
!!$#endif
!!$#endif
!!$#endif
!!$
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$#ifdef RELTIME
!!$#ifdef DIST
!!$		moutbuf(ibuf-2,i,j) = N(i,1,j)%rock_state
!!$		moutbuf(ibuf-1,i,j) = N(i,1,j)%abst
!!$		moutbuf(ibuf,i,j) = N(i,1,j)%taupos_r
!!$#ifndef SINGLEFLUID
!!$		moutbuf(ibuf-3,i,j) = N(i,1,j)%taupos_b
!!$#endif
!!$#ifndef NOSURFACTANT
!!$		moutbuf(ibuf-3,i,j) = N(i,1,j)%taupos_b
!!$		moutbuf(ibuf-4,i,j) = N(i,1,j)%taupos_s
!!$#endif
!!$#endif
!!$#endif
!!$         end do
!!$        end do
!!$
!!$	!!!!!!
!!$	! Begin asynchronous sends
!!$	!!!!!!
!!$
!!$	call MPI_ISend(		poutbuf,	& ! buf
!!$				(nx+2)*(nz+2)*ibuf,	& ! count
!!$				LBE_REAL,	& ! datatype
!!$				nnprocs(2,2),	& ! dest
!!$				tag_py,		& ! tag
!!$				Comm_Cart,	& ! communicator
!!$				requests(pout),	& ! request
!!$				ierror)
!!$        bla = 'Bad +y async send'
!!$	call checkmpi(ierror,bla)
!!$
!!$	call MPI_ISend(		moutbuf,	& ! buf
!!$				(nx+2)*(nz+2)*ibuf,	& ! count
!!$				LBE_REAL,	& ! datatype
!!$				nnprocs(2,1),	& ! dest
!!$				tag_my,		& ! tag
!!$				Comm_Cart,	& ! communicator
!!$				requests(mout),	& ! request
!!$				ierror)
!!$        bla = 'Bad -y async send'
!!$	call checkmpi(ierror,bla)
!!$
!!$	!!!!!!
!!$	! Now wait for all I/O to complete
!!$	!!!!!!
!!$
!!$	call MPI_Waitall(4,requests,statuses,ierror)
!!$        bla = 'MPI_Waitall() failed in Y direction'
!!$	call checkmpi(ierror,bla)
!!$
!!$        do j=0,nz+1
!!$         do i=0,nx+1
!!$                 N(i,0,j)%n_r(:) = pinbuf(1:19,i,j)
!!$#ifndef SINGLEFLUID
!!$                 N(i,0,j)%n_b(:) = pinbuf(20:38,i,j)
!!$#endif
!!$#ifndef NOSURFACTANT
!!$                 N(i,0,j)%n_s(:) = pinbuf(39:57,i,j)
!!$                 N(i,0,j)%d(:)   = pinbuf(58:60,i,j)
!!$#endif
!!$
!!$#ifdef LOCALBC
!!$		N(i,0,j)%local_acccoef   = pinbuf(ibuf-added,i,j)
!!$#endif
!!$
!!$#ifndef DIST
!!$#ifndef RELTIME
!!$               N(i,0,j)%rock_state = pinbuf(ibuf,i,j)
!!$#endif
!!$#endif
!!$
!!$#ifndef RELTIME
!!$#ifdef DIST
!!$		N(i,0,j)%rock_state = pinbuf(ibuf-1,i,j)
!!$		N(i,0,j)%abst = pinbuf(ibuf,i,j)
!!$#endif
!!$#endif
!!$
!!$#ifndef DIST
!!$#ifdef RELTIME
!!$		N(i,0,j)%rock_state = pinbuf(ibuf-1,i,j)
!!$		N(i,0,j)%taupos_r = pinbuf(ibuf,i,j)
!!$#ifndef SINGLEFLUID
!!$		N(i,0,j)%taupos_b = pinbuf(ibuf-2,i,j)
!!$#endif
!!$#ifndef NOSURFACTANT
!!$		N(i,0,j)%taupos_b = pinbuf(ibuf-2,i,j)
!!$		N(i,0,j)%taupos_s = pinbuf(ibuf-3,i,j)
!!$#endif
!!$#endif
!!$#endif
!!$
!!$#ifdef RELTIME
!!$#ifdef DIST
!!$		N(i,0,j)%rock_state = pinbuf(ibuf-2,i,j)
!!$		N(i,0,j)%abst = pinbuf(ibuf-1,i,j)
!!$		N(i,0,j)%taupos_r = pinbuf(ibuf,i,j)
!!$#ifndef SINGLEFLUID
!!$		N(i,0,j)%taupos_b = pinbuf(ibuf-3,i,j)
!!$#endif
!!$#ifndef NOSURFACTANT
!!$		N(i,0,j)%taupos_b = pinbuf(ibuf-3,i,j)
!!$		N(i,0,j)%taupos_s = pinbuf(ibuf-4,i,j)
!!$#endif
!!$#endif
!!$#endif
!!$
!!$!!!!!!!!!!!!!!!!!!!!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!$
!!$                 N(i,ny+1,j)%n_r(:) = minbuf(1:19,i,j)
!!$#ifndef SINGLEFLUID
!!$                 N(i,ny+1,j)%n_b(:) = minbuf(20:38,i,j)
!!$#endif
!!$#ifndef NOSURFACTANT
!!$                 N(i,ny+1,j)%n_s(:) = minbuf(39:57,i,j)
!!$                 N(i,ny+1,j)%d(:)   = minbuf(58:60,i,j)
!!$#endif
!!$
!!$#ifdef LOCALBC
!!$		N(i,ny+1,j)%local_acccoef   = minbuf(ibuf-added,i,j)
!!$#endif
!!$
!!$#ifndef DIST
!!$#ifndef RELTIME
!!$               N(i,ny+1,j)%rock_state = minbuf(ibuf,i,j)
!!$#endif
!!$#endif
!!$
!!$#ifndef RELTIME
!!$#ifdef DIST
!!$		N(i,ny+1,j)%rock_state = minbuf(ibuf-1,i,j)
!!$		N(i,ny+1,j)%abst = minbuf(ibuf,i,j)
!!$#endif
!!$#endif
!!$
!!$#ifndef DIST
!!$#ifdef RELTIME
!!$		N(i,ny+1,j)%rock_state = minbuf(ibuf-1,i,j)
!!$		N(i,ny+1,j)%taupos_r = minbuf(ibuf,i,j)
!!$#ifndef SINGLEFLUID
!!$		N(i,ny+1,j)%taupos_b = minbuf(ibuf-2,i,j)
!!$#endif
!!$#ifndef NOSURFACTANT
!!$		N(i,ny+1,j)%taupos_b = minbuf(ibuf-2,i,j)
!!$		N(i,ny+1,j)%taupos_s = minbuf(ibuf-3,i,j)
!!$#endif
!!$#endif
!!$#endif
!!$
!!$#ifdef RELTIME
!!$#ifdef DIST
!!$		N(i,ny+1,j)%rock_state = minbuf(ibuf-2,i,j)
!!$		N(i,ny+1,j)%abst = minbuf(ibuf-1,i,j)
!!$		N(i,ny+1,j)%taupos_r = minbuf(ibuf,i,j)
!!$#ifndef SINGLEFLUID
!!$		N(i,ny+1,j)%taupos_b = minbuf(ibuf-3,i,j)
!!$#endif
!!$#ifndef NOSURFACTANT
!!$		N(i,ny+1,j)%taupos_b = minbuf(ibuf-3,i,j)
!!$		N(i,ny+1,j)%taupos_s = minbuf(ibuf-4,i,j)
!!$#endif
!!$#endif
!!$#endif
!!$         end do
!!$        end do
!!$
!!$	deallocate(pinbuf)
!!$	deallocate(poutbuf)
!!$	deallocate(minbuf)
!!$	deallocate(moutbuf)
!!$
!!$	!!!!
!!$	! Y swap done.
!!$	!!!!
!!$
!!$	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$	! Swap in the Z direction
!!$	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$       allocate(pinbuf(ibuf,0:nx+1,0:ny+1))
!!$       allocate(poutbuf(ibuf,0:nx+1,0:ny+1))
!!$       allocate(moutbuf(ibuf,0:nx+1,0:ny+1))
!!$       allocate(minbuf(ibuf,0:nx+1,0:ny+1))
!!$
!!$       minz=(start(3)==1)
!!$       maxz=(start(3)>=(tnz-nz))
!!$
!!$	!!!!!!
!!$	! Begin asynchronous receives
!!$	!!!!!!
!!$
!!$	call MPI_Irecv(		pinbuf,		& ! buf
!!$				(nx+2)*(ny+2)*ibuf,	& ! count
!!$				LBE_REAL,	& ! datatype
!!$				nnprocs(3,1),	& ! source
!!$				tag_pz,		& ! tag
!!$				Comm_Cart,	& ! communicator
!!$				requests(pin),	& ! request
!!$				ierror)
!!$        bla = 'Bad +z async receive'
!!$	call checkmpi(ierror,bla)
!!$
!!$	call MPI_Irecv(		minbuf,		& ! buf
!!$				(nx+2)*(ny+2)*ibuf,	& ! count
!!$				LBE_REAL,	& ! datatype
!!$				nnprocs(3,2),	& ! source
!!$				tag_mz,		& ! tag
!!$				Comm_Cart,	& ! communicator
!!$				requests(min),	& ! request
!!$				ierror)
!!$        bla = 'Bad -z async receive'
!!$	call checkmpi(ierror,bla)
!!$
!!$        do j=0,ny+1
!!$         do i=0,nx+1
!!$                poutbuf(1:19,i,j) = N(i,j,nz)%n_r(:)
!!$#ifndef SINGLEFLUID
!!$                 poutbuf(20:38,i,j) = N(i,j,nz)%n_b(:)
!!$#endif
!!$#ifndef NOSURFACTANT
!!$                 poutbuf(39:57,i,j) = N(i,j,nz)%n_s(:)
!!$                 poutbuf(58:60,i,j) = N(i,j,nz)%d(:)
!!$#endif
!!$
!!$#ifdef LOCALBC
!!$		poutbuf(ibuf-added,i,j) = N(i,j,nz)%local_acccoef
!!$#endif
!!$
!!$#ifndef DIST
!!$#ifndef RELTIME
!!$               poutbuf(ibuf,i,j) = N(i,j,nz)%rock_state
!!$#endif
!!$#endif
!!$
!!$#ifndef RELTIME
!!$#ifdef DIST
!!$		poutbuf(ibuf-1,i,j) = N(i,j,nz)%rock_state
!!$		poutbuf(ibuf,i,j) = N(i,j,nz)%abst
!!$#endif
!!$#endif
!!$
!!$#ifndef DIST
!!$#ifdef RELTIME
!!$		poutbuf(ibuf-1,i,j) = N(i,j,nz)%rock_state
!!$		poutbuf(ibuf,i,j) = N(i,j,nz)%taupos_r
!!$#ifndef SINGLEFLUID
!!$		poutbuf(ibuf-2,i,j) = N(i,j,nz)%taupos_b
!!$#endif
!!$#ifndef NOSURFACTANT
!!$		poutbuf(ibuf-2,i,j) = N(i,j,nz)%taupos_b
!!$		poutbuf(ibuf-3,i,j) = N(i,j,nz)%taupos_s
!!$#endif
!!$#endif
!!$#endif
!!$
!!$#ifdef RELTIME
!!$#ifdef DIST
!!$		poutbuf(ibuf-2,i,j) = N(i,j,nz)%rock_state
!!$		poutbuf(ibuf-1,i,j) = N(i,j,nz)%abst
!!$		poutbuf(ibuf,i,j) = N(i,j,nz)%taupos_r
!!$#ifndef SINGLEFLUID
!!$		poutbuf(ibuf-3,i,j) = N(i,j,nz)%taupos_b
!!$#endif
!!$#ifndef NOSURFACTANT
!!$		poutbuf(ibuf-3,i,j) = N(i,j,nz)%taupos_b
!!$		poutbuf(ibuf-4,i,j) = N(i,j,nz)%taupos_s
!!$#endif
!!$#endif
!!$#endif
!!$
!!$                 moutbuf(1:19,i,j) = N(i,j,1)%n_r(:)
!!$#ifndef SINGLEFLUID
!!$                 moutbuf(20:38,i,j) = N(i,j,1)%n_b(:)
!!$#endif
!!$#ifndef NOSURFACTANT
!!$                 moutbuf(39:57,i,j) = N(i,j,1)%n_s(:)
!!$                 moutbuf(58:60,i,j) = N(i,j,1)%d(:)
!!$#endif
!!$
!!$#ifdef LOCALBC
!!$		moutbuf(ibuf-added,i,j) = N(i,j,1)%local_acccoef
!!$#endif
!!$
!!$#ifndef DIST
!!$#ifndef RELTIME
!!$               moutbuf(ibuf,i,j) = N(i,j,1)%rock_state
!!$#endif
!!$#endif
!!$
!!$#ifndef RELTIME
!!$#ifdef DIST
!!$		moutbuf(ibuf-1,i,j) = N(i,j,1)%rock_state
!!$		moutbuf(ibuf,i,j) = N(i,j,1)%abst
!!$#endif
!!$#endif
!!$
!!$#ifndef DIST
!!$#ifdef RELTIME
!!$		moutbuf(ibuf-1,i,j) = N(i,j,1)%rock_state
!!$		moutbuf(ibuf,i,j) = N(i,j,1)%taupos_r
!!$#ifndef SINGLEFLUID
!!$		moutbuf(ibuf-2,i,j) = N(i,j,1)%taupos_b
!!$#endif
!!$#ifndef NOSURFACTANT
!!$		moutbuf(ibuf-2,i,j) = N(i,j,1)%taupos_b
!!$		moutbuf(ibuf-3,i,j) = N(i,j,1)%taupos_s
!!$#endif
!!$#endif
!!$#endif
!!$
!!$#ifdef RELTIME
!!$#ifdef DIST
!!$		moutbuf(ibuf-2,i,j) = N(i,j,1)%rock_state
!!$		moutbuf(ibuf-1,i,j) = N(i,j,1)%abst
!!$		moutbuf(ibuf,i,j) = N(i,j,1)%taupos_r
!!$#ifndef SINGLEFLUID
!!$		moutbuf(ibuf-3,i,j) = N(i,j,1)%taupos_b
!!$#endif
!!$#ifndef NOSURFACTANT
!!$		moutbuf(ibuf-3,i,j) = N(i,j,1)%taupos_b
!!$		moutbuf(ibuf-4,i,j) = N(i,j,1)%taupos_s
!!$#endif
!!$#endif
!!$#endif
!!$         end do
!!$        end do
!!$
!!$	!!!!!!
!!$	! Begin asynchronous sends
!!$	!!!!!!
!!$
!!$	call MPI_ISend(		poutbuf,	& ! buf
!!$				(nx+2)*(ny+2)*ibuf,	& ! count
!!$				LBE_REAL,	& ! datatype
!!$				nnprocs(3,2),	& ! dest
!!$				tag_pz,		& ! tag
!!$				Comm_Cart,	& ! communicator
!!$				requests(pout),	& ! request
!!$				ierror)
!!$        bla = 'Bad +z async send'
!!$	call checkmpi(ierror,bla)
!!$
!!$	call MPI_ISend(		moutbuf,	& ! buf
!!$				(nx+2)*(ny+2)*ibuf,	& ! count
!!$				LBE_REAL,	& ! datatype
!!$				nnprocs(3,1),	& ! dest
!!$				tag_mz,		& ! tag
!!$				Comm_Cart,	& ! communicator
!!$				requests(mout),	& ! request
!!$				ierror)
!!$        bla = 'Bad -z async send'
!!$	call checkmpi(ierror,bla)
!!$
!!$	!!!!!!
!!$	! Now wait for all I/O to complete
!!$	!!!!!!
!!$
!!$	call MPI_Waitall(4,requests,statuses,ierror)
!!$        bla = 'MPI_Waitall() failed in Z direction'
!!$	call checkmpi(ierror,bla)
!!$
!!$! We do not want to exchange halos at z=0,tnz if we use pressure
!!$! driven flow.
!!$
!!$        do j=0,ny+1
!!$         do i=0,nx+1
!!$!!$#ifdef PERIODIC_INFLOW
!!$!!$            ! disable periodic boundaries in z direction and a
!!$!!$            ! possible halo exchange across the interface between
!!$!!$            ! periodic and non-periodic sub-volume (which is treated
!!$!!$            ! specially elsewhere)
!!$!!$            if (start(3)/=1+last_periodic_z.and..not.minz) then
!!$!!$#endif
!!$               N(i,j,0)%n_r(:) = pinbuf(1:19,i,j)
!!$#ifndef SINGLEFLUID
!!$               N(i,j,0)%n_b(:) = pinbuf(20:38,i,j)
!!$#endif
!!$#ifndef NOSURFACTANT
!!$               N(i,j,0)%n_s(:) = pinbuf(39:57,i,j)
!!$               N(i,j,0)%d(:)   = pinbuf(58:60,i,j)
!!$#endif
!!$
!!$#ifdef LOCALBC
!!$               N(i,j,0)%local_acccoef = pinbuf(ibuf-added,i,j)
!!$#endif
!!$
!!$#ifndef DIST
!!$#ifndef RELTIME
!!$
!!$  minz=(start(3)==1)
!!$  if (minz.and.(inv_fluid.eq.11.or.inv_fluid.eq.12.or.inv_fluid.eq.13.or.inv_fluid.eq.23)) then
!!$     N(i,j,0)%rock_state = N(i,j,2)%rock_state
!!$  else
!!$     N(i,j,0)%rock_state = pinbuf(ibuf,i,j)
!!$  end if
!!$
!!$#endif
!!$#endif
!!$
!!$#ifndef RELTIME
!!$#ifdef DIST
!!$               N(i,j,0)%rock_state = pinbuf(ibuf-1,i,j)
!!$               N(i,j,0)%abst = pinbuf(ibuf,i,j)
!!$#endif
!!$#endif
!!$
!!$#ifndef DIST
!!$#ifdef RELTIME
!!$               N(i,j,0)%rock_state = pinbuf(ibuf-1,i,j)
!!$               N(i,j,0)%taupos_r = pinbuf(ibuf,i,j)
!!$#ifndef SINGLEFLUID
!!$               N(i,j,0)%taupos_b = pinbuf(ibuf-2,i,j)
!!$#endif
!!$#ifndef NOSURFACTANT
!!$               N(i,j,0)%taupos_b = pinbuf(ibuf-2,i,j)
!!$               N(i,j,0)%taupos_s = pinbuf(ibuf-3,i,j)
!!$#endif
!!$#endif
!!$#endif
!!$
!!$#ifdef RELTIME
!!$#ifdef DIST
!!$               N(i,j,0)%rock_state = pinbuf(ibuf-2,i,j)
!!$               N(i,j,0)%abst = pinbuf(ibuf-1,i,j)
!!$               N(i,j,0)%taupos_r = pinbuf(ibuf,i,j)
!!$#ifndef SINGLEFLUID
!!$               N(i,j,0)%taupos_b = pinbuf(ibuf-3,i,j)
!!$#endif
!!$#ifndef NOSURFACTANT
!!$               N(i,j,0)%taupos_b = pinbuf(ibuf-3,i,j)
!!$               N(i,j,0)%taupos_s = pinbuf(ibuf-4,i,j)
!!$#endif
!!$#endif
!!$#endif
!!$!!$#ifdef PERIODIC_INFLOW
!!$!!$            end if
!!$!!$            ! disable periodic boundaries in z direction and a
!!$!!$            ! possible halo exchange across the interface between
!!$!!$            ! periodic and non-periodic sub-volume (which is treated
!!$!!$            ! specially elsewhere)
!!$!!$            if (nz+start(3)-1/=last_periodic_z.and..not.maxz) then
!!$!!$#endif
!!$               ! For pressure driven flow we interpolate densities at the end.
!!$               if (maxz.AND.&
!!$                    &((inv_fluid==20).OR.(inv_fluid==21).OR.(inv_fluid ==10))) &
!!$                    &then
!!$                  N(i,j,nz+1)%n_r(:) = 2*N(i,j,nz)%n_r(:)-N(i,j,nz-1)%n_r(:)
!!$#ifndef SINGLEFLUID
!!$                  N(i,j,nz+1)%n_b(:) = 2*N(i,j,nz)%n_b(:)-N(i,j,nz-1)%n_b(:)
!!$#endif
!!$#ifndef NOSURFACTANT
!!$                  N(i,j,nz+1)%n_s(:) = 2*N(i,j,nz)%n_s(:)-N(i,j,nz-1)%n_s(:)
!!$                  N(i,j,nz+1)%d(:) = 2*N(i,j,nz)%d(:)-N(i,j,nz-1)%d(:)
!!$#endif
!!$                  N(i,j,nz+1)%rock_state = N(i,j,nz)%rock_state
!!$               else
!!$                  N(i,j,nz+1)%n_r(:) = minbuf(1:19,i,j)
!!$#ifndef SINGLEFLUID
!!$                  N(i,j,nz+1)%n_b(:) = minbuf(20:38,i,j)
!!$#endif
!!$#ifndef NOSURFACTANT
!!$                  N(i,j,nz+1)%n_s(:) = minbuf(39:57,i,j)
!!$                  N(i,j,nz+1)%d(:)   = minbuf(58:60,i,j)
!!$#endif
!!$
!!$#ifdef LOCALBC
!!$                  N(i,j,nz+1)%local_acccoef = minbuf(ibuf-added,i,j)
!!$#endif
!!$
!!$#ifndef DIST
!!$#ifndef RELTIME
!!$                 maxz=(start(3)>=(tnz-nz))
!!$                 if (maxz.and.(inv_fluid.eq.11.or.inv_fluid.eq.12.or.inv_fluid.eq.13.or.inv_fluid.eq.23)) then
!!$                    N(i,j,nz+1)%rock_state = N(i,j,nz-1)%rock_state
!!$                 else
!!$                    N(i,j,nz+1)%rock_state = minbuf(ibuf,i,j)
!!$                 end if
!!$#endif
!!$#endif
!!$
!!$#ifndef RELTIME
!!$#ifdef DIST
!!$                  N(i,j,nz+1)%rock_state = minbuf(ibuf-1,i,j)
!!$                  N(i,j,nz+1)%abst = minbuf(ibuf,i,j)
!!$#endif
!!$#endif
!!$
!!$#ifndef DIST
!!$#ifdef RELTIME
!!$                  N(i,j,nz+1)%rock_state = minbuf(ibuf-1,i,j)
!!$                  N(i,j,nz+1)%taupos_r = minbuf(ibuf,i,j)
!!$#ifndef SINGLEFLUID
!!$                  N(i,j,nz+1)%taupos_b = minbuf(ibuf-2,i,j)
!!$#endif
!!$#ifndef NOSURFACTANT
!!$                  N(i,j,nz+1)%taupos_b = minbuf(ibuf-2,i,j)
!!$                  N(i,j,nz+1)%taupos_s = minbuf(ibuf-3,i,j)
!!$#endif
!!$#endif
!!$#endif
!!$
!!$#ifdef RELTIME
!!$#ifdef DIST
!!$                  N(i,j,nz+1)%rock_state = minbuf(ibuf-2,i,j)
!!$                  N(i,j,nz+1)%abst = minbuf(ibuf-1,i,j)
!!$                  N(i,j,nz+1)%taupos_r = minbuf(ibuf,i,j)
!!$#ifndef SINGLEFLUID
!!$                  N(i,j,nz+1)%taupos_b = minbuf(ibuf-3,i,j)
!!$#endif
!!$#ifndef NOSURFACTANT
!!$                  N(i,j,nz+1)%taupos_b = minbuf(ibuf-3,i,j)
!!$                  N(i,j,nz+1)%taupos_s = minbuf(ibuf-4,i,j)
!!$#endif
!!$#endif
!!$#endif
!!$               end if
!!$!!$#ifdef PERIODIC_INFLOW
!!$!!$            end if
!!$!!$#endif
!!$         end do
!!$        end do
!!$
!!$	deallocate(pinbuf)
!!$	deallocate(poutbuf)
!!$	deallocate(minbuf)
!!$	deallocate(moutbuf)
!!$
!!$	!!!!
!!$	! Z swap done.
!!$	!!!!
!!$
!!$	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$	! Halo exchange done.
!!$	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$end subroutine lbe_halo_exchange

#ifndef NOSURFACTANT

!> Like a halo exchange, but just swap the advected dipole moment vectors.
subroutine lbe_adv_dipole_exchange(d_adv)

	real*8,dimension(:,:,:,:),allocatable :: pinbuf
	real*8,dimension(:,:,:,:),allocatable :: poutbuf
	real*8,dimension(:,:,:,:),allocatable :: minbuf
	real*8,dimension(:,:,:,:),allocatable :: moutbuf
	real*8,dimension(1:,1:,0:,0:,0:) :: d_adv

	integer, parameter :: pin=1,mmin=2,pout=3,mout=4

	integer, dimension(4) :: requests
	integer, dimension(MPI_STATUS_SIZE,4) :: statuses
	integer :: ierror

	integer :: x,y,z,s,i
	real*8 :: t1,t2,t3

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Swap in the X direction
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	allocate(pinbuf(3,19,0:ny+1,0:nz+1))
	allocate(poutbuf(3,19,0:ny+1,0:nz+1))
	allocate(moutbuf(3,19,0:ny+1,0:nz+1))
	allocate(minbuf(3,19,0:ny+1,0:nz+1))

	!poutbuf = d_adv(:,:,nx,:,:) ! This feels like voodoo, but
	!moutbuf = d_adv(:,:,1,:,:)  ! it appears to work.

	!!!!!!
	! Begin asynchronous receives
	!!!!!!

	call MPI_Irecv(		pinbuf,		& ! buf
				(ny+2)*(nz+2)*57,&! count
				LBE_REAL,	& ! datatype
				nnprocs(1,1),	& ! source
				tag_px,		& ! tag
				Comm_Cart,	& ! communicator
				requests(pin),	& ! request
				ierror)
	call checkmpi(ierror,'(d)Bad +x async receive')

	call MPI_Irecv(		minbuf,		& ! buf
				(ny+2)*(nz+2)*57,	& ! count
				LBE_REAL,	& ! datatype
				nnprocs(1,2),	& ! source
				tag_mx,		& ! tag
				Comm_Cart,	& ! communicator
				requests(mmin),	& ! request
				ierror)
	call checkmpi(ierror,'(d)Bad -x async receive')

	do z=0,nz+1
	 do y=0,ny+1
		poutbuf(:,:,y,z)=d_adv(:,:,nx,y,z)
		moutbuf(:,:,y,z)=d_adv(:,:, 1,y,z)
	 end do
	end do

	!!!!!!
	! Begin asynchronous sends
	!!!!!!

	call MPI_ISend(		poutbuf,	& ! buf
				(ny+2)*(nz+2)*57,	& ! count
				LBE_REAL,	& ! datatype
				nnprocs(1,2),	& ! dest
				tag_px,		& ! tag
				Comm_Cart,	& ! communicator
				requests(pout),	& ! request
				ierror)
	call checkmpi(ierror,'(d)Bad +x async send')

	call MPI_ISend(		moutbuf,	& ! buf
				(ny+2)*(nz+2)*57,	& ! count
				LBE_REAL,	& ! datatype
				nnprocs(1,1),	& ! dest
				tag_mx,		& ! tag
				Comm_Cart,	& ! communicator
				requests(mout),	& ! request
				ierror)
	call checkmpi(ierror,'(d)Bad -x async send')

	!!!!!!
	! Now wait for all I/O to complete
	!!!!!!

	call MPI_Waitall(4,requests,statuses,ierror)
!	print*,'X statuses:',statuses
	call checkmpi(ierror,'(d)MPI_Waitall() failed in X direction')


	!d_adv(:,:,nx+1,:,:) = minbuf
	!d_adv(:,:,0,:,:) = pinbuf
	do z=0,nz+1
	 do y=0,ny+1
		d_adv(:,:,   0,y,z)=pinbuf(:,:,y,z)
		d_adv(:,:,nx+1,y,z)=minbuf(:,:,y,z)
	 end do
	end do

	!print*,poutbuf(:,:,:,:)-pinbuf(:,:,:,:)


	deallocate(pinbuf)
	deallocate(poutbuf)
	deallocate(minbuf)
	deallocate(moutbuf)

!	print*,'Rank ',myrankc,'completed X exchange.'

	!!!!
	! X swap done.
	!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Swap in the Y direction
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	allocate(pinbuf(3,19,0:nx+1,0:nz+1))
	allocate(poutbuf(3,19,0:nx+1,0:nz+1))
	allocate(moutbuf(3,19,0:nx+1,0:nz+1))
	allocate(minbuf(3,19,0:nx+1,0:nz+1))

	!poutbuf = d_adv(:,:,:,ny,:)
	!moutbuf = d_adv(:,:,:,1,:)



	!!!!!!
	! Begin asynchronous receives
	!!!!!!

	call MPI_Irecv(		pinbuf,		& ! buf
				(nx+2)*(nz+2)*57,	& ! count
				LBE_REAL,	& ! datatype
				nnprocs(2,1),	& ! source
				tag_py,		& ! tag
				Comm_Cart,	& ! communicator
				requests(pin),	& ! request
				ierror)
	call checkmpi(ierror,'(d)Bad +y async receive')

	call MPI_Irecv(		minbuf,		& ! buf
				(nx+2)*(nz+2)*57,	& ! count
				LBE_REAL,	& ! datatype
				nnprocs(2,2),	& ! source
				tag_my,		& ! tag
				Comm_Cart,	& ! communicator
				requests(mmin),	& ! request
				ierror)
	call checkmpi(ierror,'(d)Bad -y async receive')

	do z=0,nz+1
	 do x=0,nx+1
		poutbuf(:,:,x,z)=d_adv(:,:,x,ny,z)
		moutbuf(:,:,x,z)=d_adv(:,:,x, 1,z)
	 end do
	end do

	!!!!!!
	! Begin asynchronous sends
	!!!!!!

	call MPI_ISend(		poutbuf,	& ! buf
				(nx+2)*(nz+2)*57,	& ! count
				LBE_REAL,	& ! datatype
				nnprocs(2,2),	& ! dest
				tag_py,		& ! tag
				Comm_Cart,	& ! communicator
				requests(pout),	& ! request
				ierror)
	call checkmpi(ierror,'(d)Bad +y async send')

	call MPI_ISend(		moutbuf,	& ! buf
				(nx+2)*(nz+2)*57,	& ! count
				LBE_REAL,	& ! datatype
				nnprocs(2,1),	& ! dest
				tag_my,		& ! tag
				Comm_Cart,	& ! communicator
				requests(mout),	& ! request
				ierror)
	call checkmpi(ierror,'(d)Bad -y async send')

	!!!!!!
	! Now wait for all I/O to complete
	!!!!!!

	call MPI_Waitall(4,requests,statuses,ierror)
!	print*,'Y statuses:',statuses
	call checkmpi(ierror,'(d)MPI_Waitall() failed in Y direction')


	!d_adv(:,:,:,ny+1,:) = minbuf
	!d_adv(:,:,:,0,:) = pinbuf
	do z=0,nz+1
	 do x=0,nx+1
		d_adv(:,:,x,0,z)=pinbuf(:,:,x,z)
		d_adv(:,:,x,ny+1,z)=minbuf(:,:,x,z)
	 end do
	end do

	deallocate(pinbuf)
	deallocate(poutbuf)
	deallocate(minbuf)
	deallocate(moutbuf)

!	print*,'Rank ',myrankc,'completed Y exchange.'

	!!!!
	! Y swap done.
	!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Swap in the Z direction
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	allocate(pinbuf(3,19,0:nx+1,0:ny+1))
	allocate(poutbuf(3,19,0:nx+1,0:ny+1))
	allocate(moutbuf(3,19,0:nx+1,0:ny+1))
	allocate(minbuf(3,19,0:nx+1,0:ny+1))

	!poutbuf = d_adv(:,:,:,:,nz)
	!moutbuf = d_adv(:,:,:,:,1)


	!print*,'Rank ',myrankc,'Starting Z exchange.'

	!!!!!!
	! Begin asynchronous receives
	!!!!!!

	call MPI_Irecv(		pinbuf,		& ! buf
				(nx+2)*(ny+2)*57,	& ! count
				LBE_REAL,	& ! datatype
				nnprocs(3,1),	& ! source
				tag_pz,		& ! tag
				Comm_Cart,	& ! communicator
				requests(pin),	& ! request
				ierror)
	call checkmpi(ierror,'(d)Bad +z async receive')

	call MPI_Irecv(		minbuf,		& ! buf
				(nx+2)*(ny+2)*57,	& ! count
				LBE_REAL,	& ! datatype
				nnprocs(3,2),	& ! source
				tag_mz,		& ! tag
				Comm_Cart,	& ! communicator
				requests(mmin),	& ! request
				ierror)
	call checkmpi(ierror,'(d)Bad -z async receive')

	do y=0,ny+1
	 do x=0,nx+1
		poutbuf(:,:,x,y)=d_adv(:,:,x,y,nz)
		moutbuf(:,:,x,y)=d_adv(:,:,x,y, 1)
	 end do
	end do
	!print*,'send'

	!!!!!!
	! Begin asynchronous sends
	!!!!!!

	call MPI_ISend(		poutbuf,	& ! buf
				(nx+2)*(ny+2)*57,	& ! count
				LBE_REAL,	& ! datatype
				nnprocs(3,2),	& ! dest
				tag_pz,		& ! tag
				Comm_Cart,	& ! communicator
				requests(pout),	& ! request
				ierror)
	call checkmpi(ierror,'(d)Bad +z async send')

	call MPI_ISend(		moutbuf,	& ! buf
				(nx+2)*(ny+2)*57,	& ! count
				LBE_REAL,	& ! datatype
				nnprocs(3,1),	& ! dest
				tag_mz,		& ! tag
				Comm_Cart,	& ! communicator
				requests(mout),	& ! request
				ierror)
	call checkmpi(ierror,'(d)Bad -z async send')

	!!!!!!
	! Now wait for all I/O to complete
	!!!!!!
	!print*,'Waitall'

	call MPI_Waitall(4,requests,statuses,ierror)
!	print*,'Z statuses:',statuses
	call checkmpi(ierror,'(d)MPI_Waitall() failed in Z direction')


	!d_adv(:,:,:,:,nz+1) = minbuf
	!d_adv(:,:,:,:,0) = pinbuf

	!print*,'debuffer'

	do y=0,ny+1
	 do x=0,nx+1
		d_adv(:,:,x,y,   0)=pinbuf(:,:,x,y)
		d_adv(:,:,x,y,nz+1)=minbuf(:,:,x,y)
	 end do
	end do


	!print*,'deallocate'
	deallocate(pinbuf)
	deallocate(poutbuf)
	deallocate(minbuf)
	deallocate(moutbuf)

!	print*,'Rank ',myrankc,'completed Z exchange.'

	!!!!
	! Z swap done.
	!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Halo exchange done.
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end subroutine lbe_adv_dipole_exchange

!> Exchange both advected and real dipoles, the advected dipoles
!> being kept in the \c d_adv array.
subroutine lbe_all_dipole_exchange(d_adv,N)

	real*8,dimension(:,:,:,:),allocatable :: pinbuf
	real*8,dimension(:,:,:,:),allocatable :: poutbuf
	real*8,dimension(:,:,:,:),allocatable :: minbuf
	real*8,dimension(:,:,:,:),allocatable :: moutbuf
	real*8,dimension(1:,1:,0:,0:,0:) :: d_adv
	type(lbe_site),dimension(0:,0:,0:) :: N

	integer, parameter :: pin=1,mmin=2,pout=3,mout=4

	integer, dimension(4) :: requests
	integer, dimension(MPI_STATUS_SIZE,4) :: statuses
	integer :: ierror

	integer :: x,y,z,s,i
	real*8 :: t1,t2,t3

!	print*,'Processor ',myrankc,'beginning exchange'

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Swap in the X direction
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	allocate(pinbuf(3,nvecs+1,0:ny+1,0:nz+1))
	allocate(poutbuf(3,nvecs+1,0:ny+1,0:nz+1))
	allocate(moutbuf(3,nvecs+1,0:ny+1,0:nz+1))
	allocate(minbuf(3,nvecs+1,0:ny+1,0:nz+1))

	!poutbuf = d_adv(:,:,nx,:,:) ! This feels like voodoo, but
	!moutbuf = d_adv(:,:,1,:,:)  ! it appears to work.


	!!!!!!
	! Begin asynchronous receives
	!!!!!!

	call MPI_Irecv(		pinbuf,		& ! buf
				(ny+2)*(nz+2)*60,&! count
				LBE_REAL,	& ! datatype
				nnprocs(1,1),	& ! source
				tag_px,		& ! tag
				Comm_Cart,	& ! communicator
				requests(pin),	& ! request
				ierror)
	call checkmpi(ierror,'(d)Bad +x async receive')

	call MPI_Irecv(		minbuf,		& ! buf
				(ny+2)*(nz+2)*60,	& ! count
				LBE_REAL,	& ! datatype
				nnprocs(1,2),	& ! source
				tag_mx,		& ! tag
				Comm_Cart,	& ! communicator
				requests(mmin),	& ! request
				ierror)
	call checkmpi(ierror,'(d)Bad -x async receive')

	do z=0,nz+1
	 do y=0,ny+1
		poutbuf(:,1:nvecs,y,z)=d_adv(:,:,nx,y,z)
		moutbuf(:,1:nvecs,y,z)=d_adv(:,:, 1,y,z)
		poutbuf(:,nvecs+1,y,z)=N(nx,y,z)%d
		moutbuf(:,nvecs+1,y,z)=N(1,y,z)%d
	 end do
	end do

	!!!!!!
	! Begin asynchronous sends
	!!!!!!

	call MPI_ISend(		poutbuf,	& ! buf
				(ny+2)*(nz+2)*60,	& ! count
				LBE_REAL,	& ! datatype
				nnprocs(1,2),	& ! dest
				tag_px,		& ! tag
				Comm_Cart,	& ! communicator
				requests(pout),	& ! request
				ierror)
	call checkmpi(ierror,'(d)Bad +x async send')

	call MPI_ISend(		moutbuf,	& ! buf
				(ny+2)*(nz+2)*60,	& ! count
				LBE_REAL,	& ! datatype
				nnprocs(1,1),	& ! dest
				tag_mx,		& ! tag
				Comm_Cart,	& ! communicator
				requests(mout),	& ! request
				ierror)
	call checkmpi(ierror,'(d)Bad -x async send')

	!!!!!!
	! Now wait for all I/O to complete
	!!!!!!

	call MPI_Waitall(4,requests,statuses,ierror)
!	print*,'X statuses:',statuses
	call checkmpi(ierror,'(d)MPI_Waitall() failed in X direction')


	!d_adv(:,:,nx+1,:,:) = minbuf
	!d_adv(:,:,0,:,:) = pinbuf
	do z=0,nz+1
	 do y=0,ny+1
		d_adv(:,:,   0,y,z)=pinbuf(:,1:nvecs,y,z)
		d_adv(:,:,nx+1,y,z)=minbuf(:,1:nvecs,y,z)
		N(0,y,z)%d=pinbuf(:,nvecs+1,y,z)
		N(nx+1,y,z)%d=minbuf(:,nvecs+1,y,z)
	 end do
	end do

	!print*,poutbuf(:,:,:,:)-pinbuf(:,:,:,:)


	deallocate(pinbuf)
	deallocate(poutbuf)
	deallocate(minbuf)
	deallocate(moutbuf)

!	print*,'Rank ',myrankc,'completed X exchange.'

	!!!!
	! X swap done.
	!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Swap in the Y direction
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	allocate(pinbuf(3,nvecs+1,0:nx+1,0:nz+1))
	allocate(poutbuf(3,nvecs+1,0:nx+1,0:nz+1))
	allocate(moutbuf(3,nvecs+1,0:nx+1,0:nz+1))
	allocate(minbuf(3,nvecs+1,0:nx+1,0:nz+1))

	!poutbuf = d_adv(:,:,:,ny,:)
	!moutbuf = d_adv(:,:,:,1,:)



	!!!!!!
	! Begin asynchronous receives
	!!!!!!

	call MPI_Irecv(		pinbuf,		& ! buf
				(nx+2)*(nz+2)*60,	& ! count
				LBE_REAL,	& ! datatype
				nnprocs(2,1),	& ! source
				tag_py,		& ! tag
				Comm_Cart,	& ! communicator
				requests(pin),	& ! request
				ierror)
	call checkmpi(ierror,'(d)Bad +y async receive')

	call MPI_Irecv(		minbuf,		& ! buf
				(nx+2)*(nz+2)*60,	& ! count
				LBE_REAL,	& ! datatype
				nnprocs(2,2),	& ! source
				tag_my,		& ! tag
				Comm_Cart,	& ! communicator
				requests(mmin),	& ! request
				ierror)
	call checkmpi(ierror,'(d)Bad -y async receive')

	do z=0,nz+1
	 do x=0,nx+1
		poutbuf(:,1:nvecs,x,z)=d_adv(:,:,x,ny,z)
		moutbuf(:,1:nvecs,x,z)=d_adv(:,:,x, 1,z)
		poutbuf(:,nvecs+1,x,z)=N(x,ny,z)%d
		moutbuf(:,nvecs+1,x,z)=N(x, 1,z)%d
	 end do
	end do

	!!!!!!
	! Begin asynchronous sends
	!!!!!!

	call MPI_ISend(		poutbuf,	& ! buf
				(nx+2)*(nz+2)*60,	& ! count
				LBE_REAL,	& ! datatype
				nnprocs(2,2),	& ! dest
				tag_py,		& ! tag
				Comm_Cart,	& ! communicator
				requests(pout),	& ! request
				ierror)
	call checkmpi(ierror,'(d)Bad +y async send')

	call MPI_ISend(		moutbuf,	& ! buf
				(nx+2)*(nz+2)*60,	& ! count
				LBE_REAL,	& ! datatype
				nnprocs(2,1),	& ! dest
				tag_my,		& ! tag
				Comm_Cart,	& ! communicator
				requests(mout),	& ! request
				ierror)
	call checkmpi(ierror,'(d)Bad -y async send')

	!!!!!!
	! Now wait for all I/O to complete
	!!!!!!

	call MPI_Waitall(4,requests,statuses,ierror)
!	print*,'Y statuses:',statuses
	call checkmpi(ierror,'(d)MPI_Waitall() failed in Y direction')


	!d_adv(:,:,:,ny+1,:) = minbuf
	!d_adv(:,:,:,0,:) = pinbuf
	do z=0,nz+1
	 do x=0,nx+1
		d_adv(:,:,x,0,z)=pinbuf(:,1:nvecs,x,z)
		d_adv(:,:,x,ny+1,z)=minbuf(:,1:nvecs,x,z)
		N(x,0,z)%d=pinbuf(:,nvecs+1,x,z)
		N(x,ny+1,z)%d=minbuf(:,nvecs+1,x,z)
	 end do
	end do

	deallocate(pinbuf)
	deallocate(poutbuf)
	deallocate(minbuf)
	deallocate(moutbuf)

!	print*,'Rank ',myrankc,'completed Y exchange.'

	!!!!
	! Y swap done.
	!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Swap in the Z direction
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	allocate(pinbuf(3,nvecs+1,0:nx+1,0:ny+1))
	allocate(poutbuf(3,nvecs+1,0:nx+1,0:ny+1))
	allocate(moutbuf(3,nvecs+1,0:nx+1,0:ny+1))
	allocate(minbuf(3,nvecs+1,0:nx+1,0:ny+1))

	!poutbuf = d_adv(:,:,:,:,nz)
	!moutbuf = d_adv(:,:,:,:,1)


	!print*,'Rank ',myrankc,'Starting Z exchange.'

	!!!!!!
	! Begin asynchronous receives
	!!!!!!

	call MPI_Irecv(		pinbuf,		& ! buf
				(nx+2)*(ny+2)*60,	& ! count
				LBE_REAL,	& ! datatype
				nnprocs(3,1),	& ! source
				tag_pz,		& ! tag
				Comm_Cart,	& ! communicator
				requests(pin),	& ! request
				ierror)
	call checkmpi(ierror,'(d)Bad +z async receive')

	call MPI_Irecv(		minbuf,		& ! buf
				(nx+2)*(ny+2)*60,	& ! count
				LBE_REAL,	& ! datatype
				nnprocs(3,2),	& ! source
				tag_mz,		& ! tag
				Comm_Cart,	& ! communicator
				requests(mmin),	& ! request
				ierror)
	call checkmpi(ierror,'(d)Bad -z async receive')

	do y=0,ny+1
	 do x=0,nx+1
		poutbuf(:,1:nvecs,x,y)=d_adv(:,:,x,y,nz)
		moutbuf(:,1:nvecs,x,y)=d_adv(:,:,x,y, 1)
		poutbuf(:,nvecs+1,x,y)=N(x,y,nz)%d
		moutbuf(:,nvecs+1,x,y)=N(x,y, 1)%d
	 end do
	end do

	!!!!!!
	! Begin asynchronous sends
	!!!!!!

	call MPI_ISend(		poutbuf,	& ! buf
				(nx+2)*(ny+2)*60,	& ! count
				LBE_REAL,	& ! datatype
				nnprocs(3,2),	& ! dest
				tag_pz,		& ! tag
				Comm_Cart,	& ! communicator
				requests(pout),	& ! request
				ierror)
	call checkmpi(ierror,'(d)Bad +z async send')

	call MPI_ISend(		moutbuf,	& ! buf
				(nx+2)*(ny+2)*60,	& ! count
				LBE_REAL,	& ! datatype
				nnprocs(3,1),	& ! dest
				tag_mz,		& ! tag
				Comm_Cart,	& ! communicator
				requests(mout),	& ! request
				ierror)
	call checkmpi(ierror,'(d)Bad -z async send')

	!!!!!!
	! Now wait for all I/O to complete
	!!!!!!

	call MPI_Waitall(4,requests,statuses,ierror)
	call checkmpi(ierror,'(d)MPI_Waitall() failed in Z direction')


	!d_adv(:,:,:,:,nz+1) = minbuf
	!d_adv(:,:,:,:,0) = pinbuf

	do y=0,ny+1
	 do x=0,nx+1
		d_adv(:,:,x,y,   0)=pinbuf(:,1:nvecs,x,y)
		d_adv(:,:,x,y,nz+1)=minbuf(:,1:nvecs,x,y)
		N(x,y,   0)%d=pinbuf(:,nvecs+1,x,y)
		N(x,y,nz+1)%d=minbuf(:,nvecs+1,x,y)
	 end do
	end do

	deallocate(pinbuf)
	deallocate(poutbuf)
	deallocate(minbuf)
	deallocate(moutbuf)

!	print*,'Rank ',myrankc,'completed Z exchange.'

	!!!!
	! Z swap done.
	!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Halo exchange done.
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end subroutine lbe_all_dipole_exchange
#endif

!> This routine fills in an array, \c pcoords, with the coordinates
!> of each processor in the Cartesian topology.
!> Note that the array is zero-based, so that pcoords(:,i) holds
!> the coords of the ith processor.
subroutine find_topology(pcoords)

	implicit none
	integer :: np,pcoords(1:,0:) ! Processor rank, and coords in topology
	integer	:: sx,sy,sz	! Coords of start of a proc subdomain
	integer :: ierror

	do np=0,nprocs-1
		call MPI_cart_coords(	Comm_cart,		&
					np,			&
					nd,			&
					pcoords(:,np),		&
					ierror)
		call checkmpi(ierror,'failed to find processor coords')
	end do

end subroutine find_topology


!> This routine fills in an array, \c cpcoords, with the coordinates
!> of each domain starting point
subroutine set_cpcoords(cpcoords)
	integer, dimension(0:,:) :: cpcoords
	integer :: np,ierror,pcoords(nd)
        ! Create Array containing a cpu's SX SY SZ

	cpcoords = 0
        do np=nprocs-1,1,-1
         call MPI_cart_coords(   Comm_cart,              &
                                 np,                     &
                                 nd,                     &
                                 pcoords,                &
                                 ierror)
            call checkmpi(ierror,'failed to find processor coords')

                !Get sx,sy,sz
                cpcoords(np,1) = pcoords(1)*nx
                cpcoords(np,2) = pcoords(2)*ny
                cpcoords(np,3) = pcoords(3)*nz
	enddo

end subroutine set_cpcoords

! !> Take an array Nm, containing data for the entire domain.
! !> Cut it into subdomain-sized chunks, and feed these to the
! !> appropriate processors.
! !> Should be called only by rank 0.
! !> Returns with rank 0's subdomain completed as well.
! subroutine send_init_lattice(Nm,N)
! 	implicit none
! 	type(lbe_site),dimension(:,:,:) :: Nm
! 	type(lbe_site),dimension(0:,0:,0:) :: N
! 	integer :: np,pcoords(nd) ! Processor rank, and coords in topology
! 	integer	:: sx,sy,sz	! Coords of start of a proc subdomain
! 	integer :: ierror,i,j,k
! 	real*8,dimension(:,:,:,:),allocatable :: buffer
!         character(LEN=80) :: bla
! 	integer :: ibuf

! #ifndef NOSURFACTANT
!         ibuf = 61
! #else
!         ibuf = 39
! #endif
! #ifdef SINGLEFLUID
! 	ibuf = 20
! #endif

! 	allocate(buffer(nx,ny,nz,ibuf))

! 	do np=nprocs-1,1,-1

! 		call MPI_cart_coords(	Comm_cart,		&
! 					np,			&
! 					nd,			&
! 					pcoords,		&
! 					ierror)
!                 bla = 'failed to find processor coords'
! 		call checkmpi(ierror,bla)

! 		sx=pcoords(1)*nx
! 		sy=pcoords(2)*ny
! 		sz=pcoords(3)*nz

! 		! Fill local array with
! 		! appropriate bits.

! 		N(1:nx,1:ny,1:nz) =		&
! 		Nm(	sx+1:sx+nx,		&
! 			sy+1:sy+ny,		&
! 			sz+1:sz+nz)

! 	        do i=1,nx
!         	 do j=1,ny
! 	          do k=1,nz
!                    buffer(i,j,k, 1:19) = N(i,j,k)%n_r(:)
! #ifndef SINGLEFLUID
!                    buffer(i,j,k,20:38) = N(i,j,k)%n_b(:)
! #endif
! #ifndef NOSURFACTANT
!                    buffer(i,j,k,39:57) = N(i,j,k)%n_s(:)
!                    buffer(i,j,k,58:60) = N(i,j,k)%d(:)
! #endif
!                    buffer(i,j,k,ibuf)  = N(i,j,k)%rock_state
!          	  end do
! 	         end do
! 	        end do

! 		! Send to the appropriate
! 		! processor.

! 		print*,'Sending init data to processor ',np
! 		call MPI_send(	buffer,			& ! buf
! 				nx*ny*nz*ibuf,		& ! length
! 				LBE_REAL,		& ! datatype
! 				np,			& ! dest
! 				tag_init,		& ! tag
! 				Comm_Cart,		& ! communicator
! 				ierror)
!                 bla = 'failed to send init block'
! 		call checkmpi(ierror,bla)
! 	end do

!         deallocate(buffer)
! 	! Now fetch mine.

! 	sx=ccoords(1)*nx
! 	sy=ccoords(2)*ny
! 	sz=ccoords(3)*nz

! 	! Fill local array with
! 	! appropriate bits.

! 	N(1:nx,1:ny,1:nz) =		&
! 	Nm(	sx+1:sx+nx,		&
! 		sy+1:sy+ny,		&
! 		sz+1:sz+nz)

! end subroutine send_init_lattice

! !> Take the state array N, and wait for rank 0 to feed
! !> me the correct chunk, which is placed into N.
! subroutine recv_init_lattice(N)
! 	type(lbe_site),dimension(0:,0:,0:) :: N
! 	integer :: ierror,stat(MPI_STATUS_SIZE)
!         character(LEN=80) :: bla
! 	real*8,dimension(:,:,:,:),allocatable :: buffer
! 	integer i,j,k
! 	integer :: ibuf

! #ifndef NOSURFACTANT
!         ibuf = 61
! #else
!         ibuf = 39
! #endif
! #ifdef SINGLEFLUID
! 	ibuf = 20
! #endif

! 	allocate(buffer(nx,ny,nz,ibuf))
! 	call MPI_Recv(	buffer,				& ! buf
! 			nx*ny*nz*ibuf,			& ! length
! 			LBE_REAL,			& ! datatype
! 			0,				& ! source
! 			tag_init,			& ! tag
! 			Comm_Cart,			& ! communicator
! 			stat,				& ! status
! 			ierror)
!         bla = 'failed to receive init block'
! 	call checkmpi(ierror,bla)

! 	do i=1,nx
! 	 do j=1,ny
! 	  do k=1,nz
! 		N(i,j,k)%n_r(:) = buffer(i,j,k, 1:19)
! #ifndef SINGLEFLUID
! 		N(i,j,k)%n_b(:) = buffer(i,j,k,20:38)
! #endif
! #ifndef NOSURFACTANT
! 		N(i,j,k)%n_s(:) = buffer(i,j,k,39:57)
! 		N(i,j,k)%d(:)  = buffer(i,j,k,58:60)
! #endif
! 		N(i,j,k)%rock_state = buffer(i,j,k,ibuf)
! 	  end do
! 	 end do
! 	end do
! 	deallocate(buffer)

! end subroutine recv_init_lattice


!> Called by CPU of zero rank. Sends the
!> appropriate bit of rock data over to other CPUs
subroutine send_rock_bit(N,x,y,z,stone,cpcoords)
  type(lbe_site), dimension(0:,0:,0:) :: N
  integer, dimension(0:,:) :: cpcoords
  integer :: stat,ierror,buffer(4),x,y,z,np
  integer :: stone
  if ((x.le.(cpcoords(0,1)+nx)).and.(y.le.(cpcoords(0,2)+ny)) &
       .and.(z.le.(cpcoords(0,3)+nz))) then
    N(x-cpcoords(0,1),y-cpcoords(0,2),z-cpcoords(0,3))%rock_state=stone
  else ! Other CPUs:
    ! Find which CPU is the right one:
    do np=nprocs-1,1,-1

      if ((x.gt.cpcoords(np,1)).and.(x.le.(cpcoords(np,1)+nx)).and. &
           (y.gt.cpcoords(np,2)).and.(y.le.(cpcoords(np,2)+ny)).and. &
           (z.gt.cpcoords(np,3)).and.(z.le.(cpcoords(np,3)+nz))) then
        ! send it:
        buffer(1) = x-cpcoords(np,1)
        buffer(2) = y-cpcoords(np,2)
        buffer(3) = z-cpcoords(np,3)
        buffer(4) = int(stone)

        call MPI_send(  buffer,                 & ! buf
             4,                      & ! length
             MPI_INTEGER,            & ! datatype
             np,                     & ! dest
             tag_rock,               & ! tag
             Comm_cart,              & ! communicator
             ierror)

        call checkmpi(ierror,'Failed to send rock block')

        ! quit loop
        exit
      endif
    enddo
  endif

end subroutine send_rock_bit

!> Called by CPU of zero rank. Sends the 'release' command.
subroutine send_rock_release()
  integer :: buffer(4),ierror,np

  ! SEND 99999 TO RELEASE WAITING PROCESSORS
  buffer = 99999
  do np=nprocs-1,1,-1
    call MPI_send(  &
      buffer,       & ! buf
      4,            & ! length
      MPI_INTEGER,  & ! datatype
      np,           & ! dest
      tag_rock,     & ! tag
      Comm_cart,    & ! communicator
      ierror)
    call checkmpi(ierror,'Failed to send release processor command')
  enddo
end subroutine send_rock_release

!> Called by CPUs of nonzero rank. Waits for rank 0 to send the
!> appropriate bit of rock data over, then fills in \c N%rock_state.
subroutine recv_rock_par(N)
  type(lbe_site),dimension(0:,0:,0:) :: N
  integer :: stat(MPI_STATUS_SIZE),ierror,buffer(4)

  buffer = 0
  ! MAGIC NUMBER 99999, see send_rock_release()
  do while (buffer(4).ne.99999)
    call MPI_Recv(  &
      buffer,       & ! buf
      4,            & ! length
      MPI_INTEGER,  & ! datatype
      0,            & ! source
      tag_rock,     & ! tag
      Comm_Cart,    & ! communicator
      stat,         & ! status
      ierror)
    call checkmpi(ierror,'MPI_Recv() in recv_rock_par failed')
    if (buffer(4).ne.99999) then
      N(buffer(1),buffer(2),buffer(3))%rock_state=dble(buffer(4))
    endif
  enddo
end subroutine recv_rock_par

!> Call to check the MPI return variable \c ierror.
!> Should probably be inlined.
!> If ierror is not MPI_SUCCESS, prints \c string and abends.
subroutine checkmpi(ierror, msg)
  implicit none

  integer, intent(in) :: ierror
  character(len=*), intent(in) :: msg
  character(len=MPI_MAX_ERROR_STRING) :: error_string
  integer :: ierror_internal, error_string_len, error_class
  
  if (ierror .ne. MPI_SUCCESS) then
    call MPI_Error_class(ierror, error_class)
    call MPI_Error_string(ierror, error_string, error_string_len, ierror_internal)
    if (ierror_internal .ne. MPI_SUCCESS) then
      error_string = "unknown"
    end if
    write(msgstr,"('MPI returned error code ', I0, ' of class ', I0, '(',A,') : ', A)") ierror, error_class, trim(error_string), trim(msg)
    call error(msgstr)
  end if

end subroutine checkmpi

!> check if  stat , which should have been given to an  allocate -command
!> as the stat-argument, indicates an allocation failure; if yes, print
!>  msg  and stop the program.
subroutine check_allocate(stat,msg)
  integer,intent(in) :: stat
  character*(*),intent(in) :: msg

  if (stat/=0) then
     call log_msg('FATAL ERROR: allocate failed: '//msg,.true.)
     call Abend
  endif
end subroutine check_allocate

!> Abnormal ending
!>
!> This routine is called if something goes so disastrously wrong that
!> the simulation must be halted. At present, it's just a wrapper for
!> \c MPI_Abort(), but any required cleaning-up or checkpoint-on-error
!> code could be put here.
subroutine Abend
    integer ierror

    call log_msg("Abnormal ending...",.true.)
    call MPI_Abort(MPI_COMM_WORLD,-1,ierror)
    stop
end subroutine Abend

!> This is a debug routine, and probably best kept out of production code.
!> It calculates the total density of each species for the CPU's
!> subdomain, not including the halo regions, and then calls
!> \c MPI_Reduce() to find the global particle densities.
!>
!> Note that if this routine gets called, ALL CPUs must call it: if only
!> some do, then they will hang in the \c MPI-Reduce() call waiting for
!>others to call it.
!>
!>Also, \c MPI_Reduce() is extremely slow.
!>
!> \warning This is the only piece of code I have written which
!>forces variables to be \c real*8 - this seems to avoid some
!>rounding errors in the \c MPI_Reduce() call. Be aware that
!>further code added here should be careful not to mix
!>up \c real and \c real*8 variables in calls, which can lead
!>to bugs on non-total-64bit systems which are a complete
!>arse to track down.
!>
!> Jens 10.05.05: The real*8 causes segmentation faults on the NEC SX8.
!> Therefore, I have changed the according lines to use "real" and LB_REAL only.
!>
!> DEBUG function - this is never called at the moment
! subroutine ins_total_species(N)
!   type(lbe_site), dimension(0:,0:,0:) :: N
!   !	real*8 :: nred=0.,nblue=0.,nsurf=0., sumred=0.,sumblue=0.,sumsurf=0.
!   real*8 :: nred=0.,nblue=0.,nsurf=0., sumred=0.,sumblue=0.,sumsurf=0.
!   !	real*8 :: least = 0.
!   real*8 :: least = 0.
!   integer :: ierror
!
!   integer :: x,y,z,s
!
!   ! Add up all the densities
!   ! There's got to be a more elegant way of doing this
!   ! summation, but I can't be arsed thinking what.
!
!   !nred=0.
!   !nblue=0.
!   !nsurf=0.
!   do x=1,nx
!     do y=1,ny
!       do z=1,nz
!         if (N(x,y,z)%rock_state == 0.0) then
!           do s=1,nvecs
!             nred = nred + N(x,y,z)%n_r(s)*g(s)
! #ifndef SINGLEFLUID
!             nblue = nblue + N(x,y,z)%n_b(s)*g(s)
! #endif
! #ifndef NOSURFACTANT
!             nsurf = nsurf + N(x,y,z)%n_s(s)*g(s)
! #endif
!           end do
!         endif
!       end do
!     end do
!   end do
!
!   ! Now, add up all the subdomains' values of nred, and place
!   ! the sum in rank 0's bigsum.
!
!   CALL MPI_Reduce(nred, & ! sendbuf
!                   sumred, & ! recvbuf
!                   1, & ! count
!                   LBE_REAL, & ! datatype
!                   MPI_SUM, & ! operation
!                   0, & ! root
!                   Comm_Cart, & ! communicator
!                   ierror)
!   msgstr = 'MPI_Reduce() of red failed'
!   CALL checkmpi(ierror,msgstr)
!
!   ! Now, add up all the subdomains' values of nblue, and place
!   ! the sum in rank 0's bigsum.
!
! #ifndef SINGLEFLUID
!   CALL MPI_Reduce(nblue, & ! sendbuf
!                   sumblue, & ! recvbuf
!                   1, & ! count
!                   !MPI_REAL8, & ! datatype
!                   LBE_REAL, & ! datatype
!                   MPI_SUM, & ! operation
!                   0, & ! root
!                   Comm_Cart, & ! communicator
!                   ierror)
!   msgstr = 'MPI_Reduce() of blue failed'
!   CALL checkmpi(ierror,msgstr)
! #endif
!
! #ifndef NOSURFACTANT
!   ! Now, add up all the subdomains' values of nsurf, and place
!   ! the sum in rank 0's bigsum.
!
!   CALL MPI_Reduce(nsurf, & ! sendbuf
!                   sumsurf, & ! recvbuf
!                   1, & ! count
!                   !MPI_REAL8, & ! datatype
!                   LBE_REAL, & ! datatype
!                   MPI_SUM, & ! operation
!                   0, & ! root
!                   Comm_Cart, & ! communicator
!                   ierror)
!   msgstr = 'MPI_Reduce() of surf failed'
!   CALL checkmpi(ierror,msgstr)
! #endif
!
!
!   if (myrankc == 0) then
!     write(msgstr,"('Red, blue, surf densities: ',F16.10,' ',F16.10,' ',F16.10)") sumred, sumblue, sumsurf
!     call log_msg(msgstr)
!
!     least = min(sumred,sumblue,sumsurf)
!
!     write(msgstr,"('Red, blue, surf ratios: ',F16.10,' ',F16.10,' ',F16.10)") sumred/least, sumblue/least, sumsurf/least
!     call log_msg(msgstr)
!
!     sumred=sumred+sumblue+sumsurf
!
!     write(msgstr,"('Total particle density: ',F16.10)") sumred
!     call log_msg(msgstr)
!     ! write(msgstr,"('Epsilon: ',F16.10)") epsilon(sumred)
!     ! call log_msg(msgstr)
!   endif
!
! end subroutine ins_total_species

!> Sum total momentum and print.
!> More or less the same as \c ins_total_species - same caveats apply.
!>
!> DEBUG function - this is never called at the moment
! subroutine ins_total_momentum(N)
!   type(lbe_site), dimension(0:,0:,0:) :: N
!   integer :: ierror
!   integer :: x,y,z,s
!
!   !real*8	:: px,py,pz,sumx,sumy,sumz
!   real*8	:: px,py,pz,sumx,sumy,sumz
!
!   px=0.
!   py=0.
!   pz=0.
!   sumx=0.
!   sumy=0.
!   sumz=0.
!
!   do x=1,nx
!     do y=1,ny
!       do z=1,nz
!         do s=1,nvecs
!           px = px +       amass_r*N(x,y,z)%n_r(s)*cx(s)*g(s)
!           py = py +       amass_r*N(x,y,z)%n_r(s)*cy(s)*g(s)
!           pz = pz +       amass_r*N(x,y,z)%n_r(s)*cz(s)*g(s)
! #ifndef SINGLEFLUID
!           px = px +       amass_b*N(x,y,z)%n_b(s)*cx(s)*g(s)
!           py = py +       amass_b*N(x,y,z)%n_b(s)*cy(s)*g(s)
!           pz = pz +       amass_b*N(x,y,z)%n_b(s)*cz(s)*g(s)
! #endif
! #ifndef NOSURFACTANT
!           px = px +       amass_s*N(x,y,z)%n_s(s)*cx(s)*g(s)
!           py = py +       amass_s*N(x,y,z)%n_s(s)*cy(s)*g(s)
!           pz = pz +       amass_s*N(x,y,z)%n_s(s)*cz(s)*g(s)
! #endif
!         end do
!       end do
!     end do
!   end do
!
!   CALL MPI_Reduce(px, & ! sendbuf
!                   sumx, & ! recvbuf
!                   1, & ! count
!                   LBE_REAL, & ! datatype
!                   MPI_SUM, & ! operation
!                   0, & ! root
!                   Comm_Cart, & ! communicator
!                   ierror)
!   msgstr = "MPI_Reduce() of X momentum failed"
!   CALL checkmpi(ierror,msgstr)
!
!   CALL MPI_Reduce(py, & ! sendbuf
!                   sumy, & ! recvbuf
!                   1, & ! count
!                   !MPI_REAL8, & ! datatype
!                   LBE_REAL, & ! datatype
!                   MPI_SUM,& ! operation
!                   0, & ! root
!                   Comm_Cart, & ! communicator
!                   ierror)
!   msgstr = "MPI_Reduce() of Y momentum failed"
!   CALL checkmpi(ierror,msgstr)
!
!   CALL MPI_Reduce(pz, & ! sendbuf
!                   sumz, & ! recvbuf
!                   1, & ! count
!                   !MPI_REAL8, & ! datatype
!                   LBE_REAL, & ! datatype
!                   MPI_SUM, & ! operation
!                   0, & ! root
!                   Comm_Cart, & ! communicator
!                   ierror)
!   msgstr = "MPI_Reduce() of Z momentum failed"
!   CALL checkmpi(ierror,msgstr)
!
!   write(msgstr,"('Total momentum: ',F16.10,', ',F16.10,', ',F16.10)") px, py, pz
!   call log_msg(msgstr,.true.)
!   if (myrankc == 0) then
!     write(msgstr,"('Total momentum: ',F16.10,', ',F16.10,', ',F16.10)") sumx, sumy, sumz
!     call log_msg(msgstr)
!   endif
!
! end subroutine ins_total_momentum

!> Take the state array N, and send it to rank 0.
!> Added 19.06.02 by Jens
subroutine send_final_lattice(N)
        type(lbe_site),dimension(0:,0:,0:) :: N
        integer :: ierror,stat,i,j,k
	real*8,dimension(:,:,:,:),allocatable :: buffer
        integer :: ibuf

#ifndef NOSURFACTANT
        ibuf = 61
#else
        ibuf = 39
#endif
#ifdef SINGLEFLUID
	ibuf = 20
#endif

	allocate(buffer(ibuf,nx,ny,nz))

                do k=1,nz
                 do j=1,ny
                  do i=1,nx
                   buffer(1:19,i,j,k) = N(i,j,k)%n_r(:)
#ifndef SINGLEFLUID
                   buffer(20:38,i,j,k) = N(i,j,k)%n_b(:)
#endif
#ifndef NOSURFACTANT
                   buffer(39:57,i,j,k) = N(i,j,k)%n_s(:)
                   buffer(58:60,i,j,k) = N(i,j,k)%d(:)
#endif
                   buffer(ibuf,i,j,k)  = N(i,j,k)%rock_state
                  end do
                 end do
                end do

        call MPI_Send(  buffer,		      & ! buf
                        nx*ny*nz*ibuf,	      & ! length
                        LBE_REAL,             & ! datatype
                        0,                    & ! dest
                        tag_post,             & ! tag
                        Comm_Cart,       & ! communicator
                        ierror)
        call checkmpi(ierror,'failed to send final data block')
	deallocate(buffer)

end subroutine send_final_lattice

!> Collects the whole lattice on \c myrankc==0, only to be called by rank 0
!>
!> \param[in] N local lattice (of rank 0)
!>
!> \param[out] Nm global lattice
!>
!> This should be called only by rankw 0, all other ranks should call
!> \c send_final_lattice(). Also rank 0's subdomain \c N will be
!> copied into \c Nm. Also the halo (extent 1) of Nm will be filled
!> accordingly.
subroutine recv_final_lattice(N,Nm)
    type(lbe_site),dimension(0:,0:,0:),intent(in) :: N
    type(lbe_site),dimension(0:,0:,0:),intent(out) :: Nm
    real*8,dimension(:,:,:,:),allocatable :: buffer
    integer :: np,pcoords(nd) ! Processor rank, and coords in topology
    integer	:: sx,sy,sz	! Coords of start of a proc subdomain
    integer :: ierror,stat(MPI_STATUS_SIZE),i,j,k
    ! ibuf must have parameter attribute or its value changes after
    ! mpi_recv() . Obviously, ifort performs some optimizations and
    ! overlooks the side effects of mpi_recv() .
#ifdef SINGLEFLUID
    integer,parameter :: ibuf = 20
#else
#ifdef NOSURFACTANT
    integer,parameter :: ibuf = 39
#else
    integer,parameter :: ibuf = 61
#endif
#endif
    ! Same problem here: mpi_recv() / mpi_cart_coords() change their
    ! source/rank argument, so the loop counter must not be passed
    ! directly to these routines (Florian; ifort 9.1, mpich 1.0.4)
    integer tmpnp

    allocate(buffer(ibuf,nx,ny,nz))

    ! Fill big array with my own bits.
    sx=ccoords(1)*nx
    sy=ccoords(2)*ny
    sz=ccoords(3)*nz
    Nm(sx+1:sx+nx,sy+1:sy+ny,sz+1:sz+nz) = N(1:nx,1:ny,1:nz)

    ! collect other bits
    do np=nprocs-1,1,-1
       ! This works with ifort. Maybe outwitting other compilers needs
       ! something more tricky, like  tmpnp = nprocs-np .
       tmpnp = np
       call MPI_cart_coords(Comm_cart,tmpnp,nd,pcoords,ierror)
       call checkmpi(ierror,'failed to find processor coords')

       ! Get data from the appropriate processor.
       call MPI_Recv&
            &(buffer&                    ! buf
            &,nx*ny*nz*ibuf&             ! length
            &,LBE_REAL&                  ! datatype
            &,tmpnp&                     ! source
            &,tag_post&                  ! tag
            &,comm_cart&                 ! communicator
            &,stat&                      ! status
            &,ierror)
       call checkmpi(ierror,'failed to receive final data block')

       ! Fill big array with appropriate bits.
       sx=pcoords(1)*nx
       sy=pcoords(2)*ny
       sz=pcoords(3)*nz
       do k=1,nz
          do j=1,ny
             do i=1,nx
                Nm(sx+i,sy+j,sz+k)%n_r = buffer(1:19,i,j,k)
#ifndef SINGLEFLUID
                Nm(sx+i,sy+j,sz+k)%n_b = buffer(20:38,i,j,k)
#endif
#ifndef NOSURFACTANT
                Nm(sx+i,sy+j,sz+k)%n_s = buffer(39:57,i,j,k)
                Nm(sx+i,sy+j,sz+k)%d = buffer(58:60,i,j,k)
#endif
                Nm(sx+i,sy+j,sz+k)%rock_state = buffer(ibuf,i,j,k)
             end do
          end do
       end do
    end do
    deallocate(buffer)

    call lbe_halo_exchange_local(Nm)
end subroutine recv_final_lattice

!> fills the halo (extent 1) of a local chunk according to periodic
!> boundaries within itself
!>
!> \param[in,out] N local chunk of arbitrary size but expected to have
!> starting indices 0 and a halo of extent 1
subroutine lbe_halo_exchange_local(N)
    type(lbe_site),dimension(0:,0:,0:),intent(inout) :: N
    integer lx,ly,lz

    lx = size(N,1)-2
    ly = size(N,2)-2
    lz = size(N,3)-2

    N(     0,  1:ly,1:lz) = N(    lx,  1:ly,1:lz)
    N(  lx+1,  1:ly,1:lz) = N(     1,  1:ly,1:lz)
    N(0:lx+1,     0,1:lz) = N(0:lx+1,    ly,1:lz)
    N(0:lx+1,  ly+1,1:lz) = N(0:lx+1,     1,1:lz)
    N(0:lx+1,0:ly+1,   0) = N(0:lx+1,0:ly+1,  lz)
    N(0:lx+1,0:ly+1,lz+1) = N(0:lx+1,0:ly+1,   1)
end subroutine lbe_halo_exchange_local

subroutine build_all_chunk_mpitypes(N, h, mpitype, extent)
  implicit none

  integer,       intent(in)    :: extent, mpitype
  type(halo),    intent(inout) :: h
  type(lbe_site),intent(inout) :: N(1-halo_extent:,1-halo_extent:,1-halo_extent:)

  integer :: he, h1

  ! Set base MPI type of the halo
  h%mpitype = mpitype

  ! Set halo extent of the halo
  h%extent = extent

  ! just abbreviations
  he = h%extent
  h1 = h%extent - 1

  ! x swaps (still restricted y/z-range: data outside is out-dated)
  call build_chunk_mpitype(N,(/1,      he/),(/1,      ny/),(/1,      nz/),h,h%ls(1))
  call build_chunk_mpitype(N,(/nx-h1,  nx/),(/1,      ny/),(/1,      nz/),h,h%us(1))
  call build_chunk_mpitype(N,(/-h1,     0/),(/1,      ny/),(/1,      nz/),h,h%lr(1))
  call build_chunk_mpitype(N,(/nx+1,nx+he/),(/1,      ny/),(/1,      nz/),h,h%ur(1))

  ! y swaps (full x-range - up-to-date data was received in x-swaps)
  call build_chunk_mpitype(N,(/-h1, nx+he/),(/1,      he/),(/1,      nz/),h,h%ls(2))
  call build_chunk_mpitype(N,(/-h1, nx+he/),(/ny-h1,  ny/),(/1,      nz/),h,h%us(2))
  call build_chunk_mpitype(N,(/-h1, nx+he/),(/-h1,     0/),(/1,      nz/),h,h%lr(2))
  call build_chunk_mpitype(N,(/-h1, nx+he/),(/ny+1,ny+he/),(/1,      nz/),h,h%ur(2))

  ! z swaps (even full y-range - after z swap all data is complete)
  call build_chunk_mpitype(N,(/-h1, nx+he/),(/-h1, ny+he/),(/1,      he/),h,h%ls(3))
  call build_chunk_mpitype(N,(/-h1, nx+he/),(/-h1, ny+he/),(/nz-h1,  nz/),h,h%us(3))
  call build_chunk_mpitype(N,(/-h1, nx+he/),(/-h1, ny+he/),(/-h1,     0/),h,h%lr(3))
  call build_chunk_mpitype(N,(/-h1, nx+he/),(/-h1, ny+he/),(/nz+1,nz+he/),h,h%ur(3))

end subroutine build_all_chunk_mpitypes

subroutine build_chunk_mpitype(N, xr, yr, zr, h, cmt)
  implicit none

  type(halo),     intent(in)  :: h
  type(lbe_site), intent(in)  :: N(1-halo_extent:, 1-halo_extent:, 1-halo_extent:)
  integer,        intent(in)  :: xr(2), yr(2), zr(2) ! chunk ranges for each dim.
  integer,        intent(out) :: cmt ! type to be built

  integer(kind=MPI_ADDRESS_KIND) :: addr1, addr2, base, offset, stride
  integer(kind=MPI_ADDRESS_KIND) :: displs(1)
  integer :: bmt ! base MPI type
  integer :: xrow_mt, xyplane_mt, xyzchunk_mt ! temporary MPI data types
  integer :: cnt
  integer :: lengths(1) ! for  mpi_type_create_hindexed()
  integer :: mpierror

  integer, parameter :: blocklength = 1

  bmt = h%mpitype

  ! mpi datatype for slices of N like  N(xr(1):xr(2),y,z)
  cnt = 1 + xr(2) - xr(1)
  call MPI_Get_address(N(1,1,1), addr1, mpierror)
  DEBUG_CHECKMPI(mpierror, "build_chunk_mpitype, MPI_Get_address : x, 1, 1, 1")
  call MPI_Get_address(N(2,1,1), addr2, mpierror)
  DEBUG_CHECKMPI(mpierror, "build_chunk_mpitype, MPI_Get_address : x, 2, 1, 1")
  stride = addr2 - addr1
  call MPI_Type_create_hvector(cnt, blocklength, stride, bmt, xrow_mt, mpierror)
  DEBUG_CHECKMPI(mpierror, "build_chunk_mpitype: MPI_Type_create_hvector : x") 

  ! mpi datatype for slices of N like  N(xr(1):xr(2),yr(1):yr(2),z)
  cnt = 1 + yr(2) - yr(1)
  call MPI_Get_address(N(1,1,1), addr1, mpierror)
  DEBUG_CHECKMPI(mpierror, "build_chunk_mpitype, MPI_Get_address : y, 1, 1, 1")
  call MPI_Get_address(N(1,2,1), addr2, mpierror)
  DEBUG_CHECKMPI(mpierror, "build_chunk_mpitype, MPI_Get_address : y, 1, 2, 1")
  stride = addr2 - addr1
  call MPI_Type_create_hvector(cnt, blocklength, stride, xrow_mt, xyplane_mt, mpierror)
  DEBUG_CHECKMPI(mpierror, "build_chunk_mpitype: MPI_Type_create_hvector : y") 

  ! mpi datatype for whole chunk  N(xr(1):xr(2),yr(1):yr(2),zr(1):zr(2))
  cnt = 1 + zr(2) - zr(1)
  call MPI_Get_address(N(1,1,1), addr1, mpierror)
  DEBUG_CHECKMPI(mpierror, "build_chunk_mpitype, MPI_Get_address : z, 1, 1, 1")
  call MPI_Get_address(N(1,1,2), addr2, mpierror)
  DEBUG_CHECKMPI(mpierror, "build_chunk_mpitype, MPI_Get_address : z, 1, 1, 2")
  stride = addr2 - addr1
  call MPI_Type_create_hvector(cnt, blocklength, stride, xyplane_mt, xyzchunk_mt, mpierror)
  DEBUG_CHECKMPI(mpierror, "build_chunk_mpitype: MPI_Type_create_hvector : z") 

  ! position of the beginning of the chunk relative to the beginning of N
  call MPI_Get_address(N, base, mpierror)
  DEBUG_CHECKMPI(mpierror, "build_chunk_mpitype, MPI_Get_address : base")
  call MPI_Get_address(N(xr(1),yr(1),zr(1)), offset, mpierror)
  DEBUG_CHECKMPI(mpierror, "build_chunk_mpitype, MPI_Get_address : offset")
  offset = offset - base

  ! lcmt becomes a datatype like xyzchunk_mt but relative to base
  cnt   = 1
  lengths = (/ 1 /)
  displs  = (/ offset /)
  call MPI_Type_create_hindexed(cnt, lengths, displs, xyzchunk_mt, cmt, mpierror)
  DEBUG_CHECKMPI(mpierror, "build_chunk_mpitype, MPI_Type_create_hindexed")
  call MPI_Type_commit(cmt, mpierror)
  DEBUG_CHECKMPI(mpierror, "build_chunk_mpitype, MPI_Type_commit")

end subroutine build_chunk_mpitype

subroutine halo_exchange(N, h)
  implicit none

  type(lbe_site), intent(inout) :: N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
  type(halo), intent(inout)     :: h

  integer, parameter :: sendcount = 1
  integer, parameter :: sendtag   = 0
  integer, parameter :: recvcount = 1
  integer, parameter :: recvtag   = 0

  integer :: k, mpierror
  integer :: status(MPI_STATUS_SIZE)

  do k = 1,3
    call MPI_Sendrecv&   ! send "downward"
         &( N, sendcount, h%ls(k), nnprocs(k,1), sendtag&
         &, N, recvcount, h%ur(k), nnprocs(k,2), recvtag&
         &, comm_cart, status, mpierror)
    if(mpierror /= 0) then
       write(msgstr, "('Halo exchange MPI_Sendrecv for k = ',I0,' , downward.')") k
       DEBUG_CHECKMPI(mpierror, msgstr)
    end if
    call MPI_Sendrecv&   ! send "upward"
         &( N, sendcount, h%us(k), nnprocs(k,2), sendtag&
         &, N, recvcount, h%lr(k), nnprocs(k,1), recvtag&
         &, comm_cart, status, mpierror)
    if(mpierror /= 0) then
       write(msgstr, "('Halo exchange MPI_Sendrecv for k = ',I0,' , upward.')") k
       DEBUG_CHECKMPI(mpierror, msgstr)
    end if
  end do

end subroutine halo_exchange

!!$subroutine xcomm(phi)
!!$  implicit none
!!$  real*8, intent(inout) :: psi(-1:,-1:,-1:)
!!$  integer, parameter :: sendcount = 1
!!$  integer, parameter :: sendtag   = 0
!!$  integer, parameter :: recvcount = 1
!!$  integer, parameter :: recvtag   = 0
!!$
!!$  integer :: k, mpierror
!!$  integer :: status(MPI_STATUS_SIZE)
!!$
!!$  call MPI_Sendrecv&   ! send "downward"
!!$       &( psi, sendcount, h%ls(k), nnprocs(k,1), sendtag&
!!$       &, psi, recvcount, h%ur(k), nnprocs(k,2), recvtag&
!!$       &, comm_cart, status, mpierror)
!!$  if(mpierror /= 0) then
!!$     write(msgstr, "('Halo exchange MPI_Sendrecv for k = ',I0,' , downward.')") k
!!$     DEBUG_CHECKMPI(mpierror, msgstr)
!!$  end if
!!$  call MPI_Sendrecv&   ! send "upward"
!!$       &( psi, sendcount, h%us(k), nnprocs(k,2), sendtag&
!!$       &, psi, recvcount, h%lr(k), nnprocs(k,1), recvtag&
!!$       &, comm_cart, status, mpierror)
!!$  if(mpierror /= 0) then
!!$     write(msgstr, "('Halo exchange MPI_Sendrecv for k = ',I0,' , upward.')") k
!!$     DEBUG_CHECKMPI(mpierror, msgstr)
!!$  end if
!!$
!!$
!!$end subroutine xcomm
!!$
!!$subroutine ycomm(phi)
!!$  implicit none
!!$  real*8, intent(inout) :: psi(-1:,-1:,-1:)
!!$
!!$
!!$end subroutine ycomm
!!$
!!$subroutine zcomm(phi)
!!$  implicit none
!!$  real*8, intent(inout) :: psi(-1:,-1:,-1:)
!!$
!!$
!!$end subroutine zcomm

end module lbe_parallel_module
