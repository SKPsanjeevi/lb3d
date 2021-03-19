#include "lbe.h"

!> commonly used parallelization routines for md particles
module lbe_md_parallel_module
#ifdef MD
    use lbe_globals_module
    use lbe_md_bc_leesedwards_module, only: md_leesedwards_adapt_after_gather
    use lbe_md_globals_module
    use lbe_md_helper_module, only: error_md
#ifdef MPI_ALLGV_FASTER_THAN_GV
    use lbe_md_helper_module, only: count_particles_all
#else
    use lbe_md_helper_module, only: count_particles
#endif
    use lbe_md_memory_module, only: boost_npmax
    use lbe_parallel_module, only: calculate_displacements,checkmpi&
         &,check_allocate,comm_cart,find_topology,nprocs

    implicit none
    private
    public gather_particles,gather_particles_provide_rbuf,scatter_particles

    include 'mpif.h'

contains

    !> Gathers data of all particles from all processors into the elements
    !> of  tp  on the root process. The optional custom mpi type  selection
    !> controls what data is copied, default is  particle_exch_mpitype .
    subroutine gather_particles(tp,selection)
        type(md_particle_type),intent(out) :: tp(1:)
        integer,intent(in),optional :: selection
        integer selected_mpitype
        integer i,ii,scount,stype,ierror,stat
        ! for  mpi_type_indexed() :
        integer,allocatable,dimension(:),save :: ibuf,lengths
        ! only for root:
        integer counts(0:nprocs-1),displs(0:nprocs-1)

        ! default selection of particle data to gather is  particle_exch_mpitype
        if (present(selection)) then
           selected_mpitype = selection
        else
           selected_mpitype = particle_exch_mpitype
        end if

        ! allocate  ibuf  and  lengths  on first invocation
        if (.not.allocated(ibuf)) then
           allocate (ibuf(nlocal),lengths(nlocal),stat=stat)
           call check_allocate(stat,'gather_particles(): ibuf,lengths')
           lengths = 1        !  lengths(:)  remains  1  for all times
        end if

        ! allow  nlocal  to alter during the simulation
        if (nlocal>size(ibuf)) then
           deallocate (ibuf,lengths)
           allocate (ibuf(2*nlocal),lengths(2*nlocal),stat=stat)
           call check_allocate(stat,'gather_particles(): ibuf,lengths')
           lengths = 1        !  lengths(:)  remains  1  for all times
        end if

        if (nlocal>0) then
           ! create one type that refers to all local particles
           i = atompnt
           do ii=1,nlocal
              ibuf(ii) = i
              i = list(i)
           end do

           call mpi_type_indexed(nlocal,lengths,ibuf,selected_mpitype,stype&
                &,ierror)
           call checkmpi(ierror,'gather_particles(): mpi_type_indexed() failed')

           call mpi_type_commit(stype,ierror)
           call checkmpi(ierror,'gather_particles(): mpi_type_commit() failed')
           scount = 1
        else
           ! sending one element of an empty type confuses
           ! mpi_gatherv() - send zero elements of a dummy type instead
           stype = MPI_INTEGER
           scount = 0
        end if

        ! first communicate particle counts, then actual particle data
        ! - because -1 was not subtracted when filling ibuf, use P(0)
        ! as sendbuf instead of P(1).
#ifdef MPI_ALLGV_FASTER_THAN_GV
        ! MPI_allgatherv() produces unnecessary communication compared
        ! to MPI_gatherv(), but on jugene it is significantly faster.
        call mpi_allgather(nlocal,1,MPI_INTEGER&
             &,counts,1,MPI_INTEGER,comm_cart,ierror)
        call checkmpi(ierror,'gather_particles(): mpi_allgather() failed')
        call calculate_displacements(counts,displs)
        call mpi_allgatherv(P(0),scount,stype&
             &,tp,counts,displs,selected_mpitype,comm_cart,ierror)
        call checkmpi(ierror,'gather_particles(): mpi_allgatherv() failed')
#else
        call mpi_gather(nlocal,1,MPI_INTEGER,counts,1,MPI_INTEGER,0,comm_cart&
             &,ierror)
        call checkmpi(ierror,'gather_particles(): mpi_gather() failed')
        if (myrankc==0) call calculate_displacements(counts,displs)
        call mpi_gatherv(P(0),scount,stype&
             &,tp,counts,displs,selected_mpitype,0,comm_cart,ierror)
        call checkmpi(ierror,'gather_particles(): mpi_gatherv() failed')
#endif

        if (nlocal>0) then
           call mpi_type_free(stype,ierror)
           call checkmpi(ierror,'gather_particles(): mpi_type_free() failed')
        end if

        if (md_leesedwards) call md_leesedwards_adapt_after_gather(tp,counts&
             &,displs)
    end subroutine gather_particles

    !> makes sure that a recv buffer for particle data is allocated
    !> and sufficiently large to be used by \c gather_particles()
    !>
    !> Depending on the actual implementation of \c
    !> gather_particles(), the buffer is allocated on rank zero or on
    !> all ranks.
    !>
    !> \param[inout] tp recv buffer
    !> \param[in] errstr string to be used in possible error messages
    subroutine gather_particles_provide_rbuf(tp,errstr)
        type(md_particle_type),allocatable,intent(inout) :: tp(:)
        character*(*),intent(in) :: errstr
        integer :: capacity,n,stat

#ifdef MPI_ALLGV_FASTER_THAN_GV
        ! If gather_particles() uses MPI_allgatherv() instead of
        ! MPI_gatherv() all tasks need to allocate a recv buffer
        call count_particles_all(n)
#else
        call count_particles(n,0)
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
                &,'gather_particles_provide_rbuf(): tp ('//errstr//')')
        end if
    end subroutine gather_particles_provide_rbuf

    !> distributes particles from the root process to the others
    !> depending on their positions
    !>
    !> \param[in] n number of particles to send (only read by root)
    !>
    !> \param[in] tp array of particles to send (only read by root)
    !>
    !> \param[in] selection [optional, defaults to \c
    !> particle_exch_mpitype] custom mpi type that determines which
    !> elements are send for each particle.
    subroutine scatter_particles(n,tp,selection)
        integer,intent(in) :: n
        type(md_particle_type),intent(in) :: tp(:)
        integer,intent(in),optional :: selection
        integer selected_mpitype
        ! for  mpi_type_indexed() :
        integer,allocatable,dimension(:),save :: ibuf,lengths
        integer i,j,r,ierror,stat,status(MPI_STATUS_SIZE)
        integer pcoords(3,0:nprocs-1) ! cartesian topology coordinates
        real(kind=rk) :: pmin(3),pmax(3)
        integer counts(0:nprocs-1) ! # of particles sent to each processor
        ! mpi indexed type sent to each processor except root (-> nproc-1  el.)
        integer stype(1:nprocs-1)
        ! index in  ibuf  where the chunk for each processor starts:
        integer pstart(0:nprocs-1)
        integer,parameter :: tag=0
        integer itmp,new_nlocal

        ! default selection of particle data to scatter is  particle_exch_mpitype
        if (present(selection)) then
           selected_mpitype = selection
        else
           selected_mpitype = particle_exch_mpitype
        end if

        rank_0a: if (myrankc==0) then
           ! allocate  ibuf  and  lengths  on first invocation
           if (.not.allocated(ibuf)) then
              allocate (ibuf(n),lengths(n),stat=stat)
              call check_allocate(stat,'scatter_particles(): ibuf,lengths')
              lengths(:) = 1  !  lengths(:)  remains  1  for all times
           end if

           ! allow  natoms  to alter during the simulation
           if (n>size(ibuf)) then
!!$              call log_msg_md&
!!$                   &('Index buffer is too small, allocating more memory...'&
!!$                   &,.true.)
              deallocate (ibuf,lengths)
              allocate (ibuf(2*n),lengths(2*n),stat=stat)
              call check_allocate(stat,'scatter_particles(): ibuf,lengths')
              lengths(:) = 1  !  lengths(:)  remains  1  for all times
           end if

           ! create a type for particles to be sent to rank  r  in  stype(r)
           call find_topology(pcoords)
           counts(:) = 0
           j = 0
           rank: do r=0,nprocs-1
              ! start index in  ibuf  for rank  r
              pstart(r) = j + 1
              ! cell boundaries for rank  r
              pmin(:) = pcoords(:,r)*chunksize(:) + 1.0_rk - 0.5_rk
              pmax(:) = pmin(:) + chunksize(:)
              ! iterate through all and index the right particles...
              particle: do i=1,n
                 if (all(tp(i)%x(:)>=pmin(:)).and.all(tp(i)%x(:)<pmax(:))) then
                    j = j + 1
                    ibuf(j) = i
                    counts(r) = counts(r) + 1
                 end if
              end do particle
              ! now really create the mpi type (except for root)
              if (r>0) then
                 call mpi_type_indexed(counts(r)&
                      &,lengths(pstart(r)),ibuf(pstart(r))&
                      &,selected_mpitype,stype(r),ierror)
                 call mpi_type_commit(stype(r),ierror)
              end if
           end do rank

           if (sum(counts)/=n) then
              write (6,*) 'global number of particles natoms=',n,', should be'&
                   &,sum(counts)
              call error_md('conservation of particle count violated in '&
                   &//'scatter_particles()')
           end if
        end if rank_0a

        call mpi_scatter(counts,1,MPI_INTEGER,new_nlocal,1,MPI_INTEGER,0&
             &,comm_cart,ierror)

        if (new_nlocal>npmax) call boost_npmax(new_nlocal)
        nlocal = new_nlocal

        rank_0b: if (myrankc==0) then
           ! these particles are already at the correct processor
           do j=1,counts(0)
              P(j) = tp(ibuf(j)) ! pstart(0)==1 -> no explicit pstart needed
           end do

           ! rest of the particles needs to be sent to other processors
           do r=1,nprocs-1
              call mpi_send(tp(0),1,stype(r),r,tag,comm_cart,ierror)
              call mpi_type_free(stype(r),ierror)
           end do
        else
           ! non-root processors receive particles from root
           call mpi_recv(P,nlocal,selected_mpitype,0,tag,comm_cart,status&
                &,ierror)
        end if rank_0b

        ! build linked list for the  nlocal  particles...
        do i=1,nlocal
           itmp = atompnt
           atompnt = freepnt
           freepnt = list(freepnt)
           list(atompnt) = itmp
        end do
    end subroutine scatter_particles

#endif
end module lbe_md_parallel_module
