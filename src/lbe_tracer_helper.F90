#include "lbe.h"

!> helper routines for \c TRACER
module lbe_tracer_helper_module
#ifdef TRACER

    use lbe_collision_module, only:  md_calculate_sc_forces
    use lbe_globals_module,only: c,cx,cy,cz,g,halo_extent,lp_sur,n_lp&
         &,nnonrest,opp_lp_sur,pi,tsize,myrankc
    use lbe_helper_module, only: is_fluid,local_coordinates,present_and_true&
         &,velocity
    use lbe_log_module
#ifdef MD
    use lbe_md_fluid_module, only: md_velocity
#endif
    use lbe_parallel_module, only: calculate_displacements,ccoords&
         &,check_allocate,checkmpi,comm_cart,nprocs&
         &,Abend
    use lbe_parms_module, only: amass_r,amass_b,amass_s,g_accn_min,g_accn_max&
         &,g_accn_min_x,g_accn_max_x,g_accn_min_y,g_accn_max_y,g_accn,g_accn_x&
         &,g_accn_y,nx,ny,nz
    use lbe_tracer_globals_module
    use lbe_types_module, only: lbe_site

    implicit none
    include 'mpif.h'
    private

    public build_tracer_mpitype,count_tracers,count_tracers_all,error_tracer&
         &,fluid_velocity,gather_tracers,gather_tracers_provide_rbuf&
         &,log_msg_tracer,log_msg_tracer_ws,log_msg_tracer_hdr,passed_plane

    character(len=8), parameter :: tracer_log_prefix = "[TRACER]"

contains

!> Log message with TRACER prefix and surrounding whitespace.
subroutine log_msg_tracer_ws(msg, forall)
  implicit none

  character(len=*), intent(in)  :: msg
  logical, optional, intent(in) :: forall

  if ( present(forall) ) then
    call log_msg("", forall)
    call log_msg(tracer_log_prefix // " " // trim(msg), forall)
    call log_msg("", forall)
  else
    call log_msg("")
    call log_msg(tracer_log_prefix // " " // trim(msg))
    call log_msg("")
  end if

end subroutine log_msg_tracer_ws

!> Log message with TRACER prefix.
subroutine log_msg_tracer(msg, forall, error)
  implicit none

  character(len=*), intent(in)  :: msg    !< Text to be logged.
  logical, optional, intent(in) :: forall !< If supplied and set to true, all ranks will log.
  logical, optional, intent(in) :: error  !< If supplied and set to true, log an error.

  if ( present(forall) ) then
    if ( present(error) ) then
      call log_msg(tracer_log_prefix // " " // trim(msg), forall = forall, error = error)
    else
      call log_msg(tracer_log_prefix // " " // trim(msg), forall = forall)
    end if
  else
    if ( present(error) ) then
      call log_msg(tracer_log_prefix // " " // trim(msg), error = error)
    else
      call log_msg(tracer_log_prefix // " " // trim(msg))
    end if
  end if

end subroutine log_msg_tracer

!> Log header message with TRACER prefix.
subroutine log_msg_tracer_hdr(msg, forall)
  implicit none

  character(len=*), intent(in)  :: msg    !< Text to be logged.
  logical, optional, intent(in) :: forall !< If supplied and set to true, all ranks will log.

  if ( present(forall) ) then
    call log_msg_hdr(tracer_log_prefix // " " // trim(msg), forall)
  else
    call log_msg_hdr(tracer_log_prefix // " " // trim(msg))
  end if

end subroutine log_msg_tracer_hdr

!> Log error message with TRACER prefix and exit.
subroutine error_tracer(msg)
  implicit none

  character(len=*), intent(in)  :: msg    !< Text to be logged.

  integer :: mpierror

  call log_msg_tracer("ERROR: "//trim(msg), forall = .true., error = .false.)
  call log_msg_tracer("ERROR: "//trim(msg), forall = .true., error = .true. )
  call MPI_Abort(MPI_COMM_WORLD, -1, mpierror)
  stop

end subroutine error_tracer

    !> build a custom mpi datatype representing a selection of the
    !> elements of \c tracer_particle_type controlled by the optional
    !> arguments
    !>
    !> \param[out] ptt custom mpi type
    !>
    !> \param[in] x include \c x ?
    !>
    !> \param[in] v include \c v ?
    !>
    !> \param[in] uid include \c uid ?
    !>
    !> \param[in] kind include \c kind ?
  subroutine build_tracer_mpitype(ptt,x,v,uid,akind)
      integer,intent(out) :: ptt
      logical,intent(in),optional :: x,v,uid,akind
      integer,parameter :: n_blocks_max=6
      integer n_blocks,ierror
      integer lengths(n_blocks_max),types(n_blocks_max)
      integer(kind=MPI_ADDRESS_KIND) :: base,addrs(n_blocks_max)&
           &,displs(n_blocks_max)
      type(tracer_particle_type) :: sample(2)

      n_blocks = 1
      lengths(n_blocks) = 1   ! start of particle in memory
      types(n_blocks) = MPI_LB
      call mpi_get_address(sample(1),addrs(n_blocks),ierror)

      if (present_and_true(x)) then
         n_blocks = n_blocks+1
         lengths(n_blocks) = 3          ! x
         types(n_blocks) = MPI_REAL8
         call mpi_get_address(sample(1)%x(1),addrs(n_blocks),ierror)
      end if

      if (present_and_true(v)) then
         n_blocks = n_blocks+1
         lengths(n_blocks) = 3          ! v
         types(n_blocks) = MPI_REAL8
         call mpi_get_address(sample(1)%v(1),addrs(n_blocks),ierror)
      end if

      if (present_and_true(uid)) then
         n_blocks = n_blocks+1
         lengths(n_blocks) = 1          ! uid
         types(n_blocks) = MPI_INTEGER
         call mpi_get_address(sample(1)%uid,addrs(n_blocks),ierror)
      end if

      if (present_and_true(akind)) then
         n_blocks = n_blocks+1
         lengths(n_blocks) = 1          ! kind
         types(n_blocks) = MPI_INTEGER
         call mpi_get_address(sample(1)%kind,addrs(n_blocks),ierror)
      end if

      n_blocks = n_blocks+1
      lengths(n_blocks) = 1   ! next particle
      types(n_blocks) = MPI_UB
      call mpi_get_address(sample(2),addrs(n_blocks),ierror)

      call mpi_get_address(sample(1),base,ierror) ! base address
      displs(1:n_blocks) = addrs(1:n_blocks) - base

      call mpi_type_create_struct(n_blocks,lengths,displs,types,ptt,ierror)
      call checkmpi(ierror&
           &,'build_tracer_mpitype(): mpi_type_create_struct() failed')

      call mpi_type_commit(ptt,ierror)
      call checkmpi(ierror&
           &,'build_tracer_mpitype(): mpi_type_commit() failed')
  end subroutine build_tracer_mpitype

  !> returns in \c n on process \c dst in \c comm_cart the total
  !> number of \c TRACER particles on all processes.
  subroutine count_tracers(n,dst,mod_uid)
      integer,intent(out) :: n
      integer,intent(in) :: dst
      integer,intent(in),optional :: mod_uid
      integer ierror,nl

      nl = nlocal_mod_uid(mod_uid)
      call mpi_reduce(nl,n,1,MPI_INTEGER,MPI_SUM,dst,comm_cart,ierror)
  end subroutine count_tracers

  !> returns in \c n the total number of \c TRACER particles on all
  !> processes and distributes it to all processes.
  subroutine count_tracers_all(n,mod_uid)
      integer,intent(out) :: n
      integer,intent(in),optional :: mod_uid
      integer ierror,nl

      nl = nlocal_mod_uid(mod_uid)
      call mpi_allreduce(nl,n,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierror)
  end subroutine count_tracers_all

  !> return the interpolated fluid velocity at global postion  x  in  v  and
  !>  .true.  in  rock  if  x  is enclosed by rock sites.  x  must be at
  !> some place that is covered by the local  N , however, this subroutine
  !> is smart enough to deal with pbc. Velocities of the lbe particle
  !> species r/b/s are weighted as  mask(1/2/3) .
  !> Set  halo_extent>=2  if you use this subroutine for arbitrary positions
  !> inside each local domain!
  subroutine fluid_velocity(N,x,v)
      type(lbe_site),intent(in) :: N(&
           &1-halo_extent:,1-halo_extent:,1-halo_extent:)
      real(kind=rk),intent(in) :: x(3)
      real(kind=rk),intent(out) :: v(3)
      integer l
      integer lx(3),lp(3),opp_lp(3) ! lattice point coordinates
      real(kind=rk) :: weight,xx(3)

      call local_coordinates(x,xx)
      lx = floor(xx)

      v = 0.0_rk
      lattice_points: do l=1,n_lp
         lp = lx + lp_sur(:,l)
         ! weight each lattice point with the volume of the
         ! hypercube spanned by the particle position and the
         ! opposite lattice point
         opp_lp = lx + opp_lp_sur(:,l)
         weight = abs(product(real(opp_lp,kind=rk)-xx))
#ifdef MD
         ! this assumes zero velocity at solid rock and rigid body
         ! motion at md particle sites
         v = v + weight*md_velocity(N,lp(1),lp(2),lp(3))
#else
         ! this assumes just zero velocity on rock sites
         if (is_fluid(N(lp(1),lp(2),lp(3))%rock_state)) &
              &v = v + weight*velocity(N(lp(1),lp(2),lp(3)))
#endif
      end do lattice_points
  end subroutine fluid_velocity

  !> Gathers data of all tracers from all processors on the root process.
  !>
  !> \param[out] tp array to hold tracer data, allocated by caller,
  !> only required on root process unless compiled with
  !> MPI_ALLGV_FASTER_THAN_GV
  !>
  !> \param[in] selection custom mpi type controlling which data out
  !> of \c tracer_particle_type is actually copied [optional,
  !> defaults to \c tracer_exch_mpitype ]
  !>
  !> \param[in] mod_uid gather only tracers for which \c
  !> mod(uid,mod_uid)==0 [optional, defaults to 1]
  subroutine gather_tracers(tp,selection,mod_uid)
      type(tracer_particle_type),intent(out) :: tp(1:)
      integer,intent(in),optional :: selection
      integer,intent(in),optional :: mod_uid
      integer selected_mpitype
      integer i,ii,iii,nl,scount,stype,ierror,stat
      ! for  mpi_type_indexed() :
      integer,allocatable,dimension(:),save :: ibuf,lengths
      ! only for root:
      integer counts(0:nprocs-1),displs(0:nprocs-1)

      if (present(selection)) then
         selected_mpitype = selection
      else
         selected_mpitype = tracer_exch_mpitype
      end if

      nl = nlocal_mod_uid(mod_uid)

      ! allocate  ibuf  and  lengths  on first invocation
      if (.not.allocated(ibuf)) then
         allocate (ibuf(nl),lengths(nl),stat=stat)
         call check_allocate(stat,'gather_tracers(): ibuf,lengths')
         lengths = 1        !  lengths(:)  remains  1  for all times
      end if

      ! nl can change during the simulation
      if (nl>size(ibuf)) then
         deallocate (ibuf,lengths)
         allocate (ibuf(2*nl),lengths(2*nl),stat=stat)
         call check_allocate(stat,'gather_tracers(): ibuf,lengths')
         lengths = 1        !  lengths(:)  remains  1  for all times
      end if

      if (nl>0) then
         ! create one type that refers to all local tracers
         i = atompnt
         if (present(mod_uid)) then
            iii = 0
            do ii=1,nlocal
               if (mod(T(i)%uid,mod_uid)==0) then
                  iii = iii+1
                  ibuf(iii) = i
               end if
               i = list(i)
            end do
         else
            do ii=1,nlocal
               ibuf(ii) = i
               i = list(i)
            end do
         end if

         call mpi_type_indexed(nl,lengths,ibuf,selected_mpitype,stype&
              &,ierror)
         call checkmpi(ierror,'gather_tracers(): mpi_type_indexed() failed')

         call mpi_type_commit(stype,ierror)
         call checkmpi(ierror,'gather_tracers(): mpi_type_commit() failed')
         scount = 1
      else
         ! sending one element of an empty type confuses
         ! mpi_gatherv() - send zero elements of a dummy type instead
         stype = MPI_INTEGER
         scount = 0
      end if

      ! first communicate tracer counts, then actual tracer data.
      ! Because -1 was not subtracted when filling ibuf, use T(0) as
      ! sendbuf instead of T(1).
#ifdef MPI_ALLGV_FASTER_THAN_GV
      ! MPI_allgatherv() produces unnecessary communication compared
      ! to MPI_gatherv(), but on jugene it is significantly faster.
      call mpi_allgather(nl,1,MPI_INTEGER&
           &,counts,1,MPI_INTEGER,MPI_COMM_WORLD,ierror)
      call checkmpi(ierror,'gather_tracers(): mpi_allgather() failed')
      call calculate_displacements(counts,displs)
      call mpi_allgatherv(T(0),scount,stype&
           &,tp,counts,displs,selected_mpitype,MPI_COMM_WORLD,ierror)
      call checkmpi(ierror,'gather_tracers(): mpi_allgatherv() failed')
#else
      call mpi_gather(nl,1,MPI_INTEGER,counts,1,MPI_INTEGER,0,comm_cart&
           &,ierror)
      call checkmpi(ierror,'gather_tracers(): mpi_gather() failed')
      if (myrankc==0) call calculate_displacements(counts,displs)
      call mpi_gatherv(T(0),scount,stype&
           &,tp,counts,displs,selected_mpitype,0,comm_cart,ierror)
      call checkmpi(ierror,'gather_tracers(): mpi_gatherv() failed')
#endif

      if (nl>0) then
         call mpi_type_free(stype,ierror)
         call checkmpi(ierror,'gather_tracers(): mpi_type_free() failed')
      end if
  end subroutine gather_tracers

  !> makes sure that a recv buffer for tracer data is allocated
  !> and sufficiently large to be used by \c gather_tracers()
  !>
  !> Depending on the actual implementation of \c
  !> gather_tracers(), the buffer is allocated on rank zero or on
  !> all ranks.
  !>
  !> \param[in,out] tp recv buffer
  !>
  !> \param[in] errstr string to be used in possible error messages
  !>
  !> \param[in] mod_uid take only into account tracers for which \c
  !> mod(uid,mod_uid)==0 [optional, defaults to 1]
  subroutine gather_tracers_provide_rbuf(tp,errstr,mod_uid)
      type(tracer_particle_type),allocatable,intent(inout) :: tp(:)
      character*(*),intent(in) :: errstr
      integer,intent(in),optional :: mod_uid
      integer :: capacity,n,stat

#ifdef MPI_ALLGV_FASTER_THAN_GV
      ! If gather_tracers() uses MPI_allgatherv() instead of
      ! MPI_gatherv() all tasks need to allocate a recv buffer
      call count_tracers_all(n,mod_uid)
#else
      call count_tracers(n,0,mod_uid)
      if (myrankc/=0) return
#endif
      capacity = n
      if (allocated(tp).and.n>size(tp)) then
         ! allow number of particles to change during a simulation
         deallocate (tp)
         capacity = capacity*2
      end if

      if (.not.allocated(tp)) then
         allocate (tp(capacity),stat=stat)
         call check_allocate(stat&
              &,'gather_tracers_provide_rbuf(): tp ('//errstr//')')
      end if
  end subroutine gather_tracers_provide_rbuf

    !> computes the number of local tracers for which \c
    !> mod(uid,mod_uid)==0
    !>
    !> \param[in] mod_uid defines the tracer uid filter, defaults to 1
    !> which means that all tracers are counted
    !>
    !> \returns local number of tracers with according uid
    integer function nlocal_mod_uid(mod_uid)
        integer,intent(in),optional :: mod_uid
        integer i,ii,n

        if (.not.present(mod_uid)) then
           nlocal_mod_uid = nlocal
        else if (mod_uid==1) then
           nlocal_mod_uid = nlocal
        else
           n = 0
           i = atompnt
           do ii=1,nlocal
              if (mod(T(i)%uid,mod_uid)==0) n = n+1

              i = list(i)
           end do
           nlocal_mod_uid = n
        end if
    end function nlocal_mod_uid

    !> checks whether two positions lie on opposite sides of a plane
    !> perpendicular to the x-direction
    !>
    !> \param[in] tx first position perpendicular to the plane
    !>
    !> \param[in] px second position perpendicular to the plane
    !>
    !> \param[in] plane_x position of the plane
    pure function passed_plane(tx,px,plane_x)
        logical :: passed_plane
        real(kind=rk),intent(in) :: tx,px,plane_x
        real(kind=rk) :: dtx,dpx

        ! relative position now and before
        dtx = tx-plane_x
        dpx = px-plane_x

        passed_plane = dtx*dpx<0.0_rk
    end function passed_plane

#endif
end module lbe_tracer_helper_module
