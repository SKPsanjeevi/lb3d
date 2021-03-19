#include "lbe.h"

!> commonly used parallelization routines for \c TRACER particles
module lbe_tracer_parallel_module
#ifdef TRACER
    use lbe_globals_module, only: chunksize,myrankc
    use lbe_parallel_module, only: check_allocate,comm_cart,find_topology&
         &,nprocs
    use lbe_tracer_globals_module
    use lbe_tracer_helper_module, only: error_tracer
    use lbe_tracer_memory_module, only: boost_npmax

    implicit none
    private
    public scatter_tracers

    include 'mpif.h'

contains

    !> distributes tracers from the root process to the others
    !> depending on their positions
    !>
    !> \param[in] n number of tracers to send (only read by root)
    !>
    !> \param[in] tp array of tracers to send (only read by root)
    !>
    !> Only the elements of \c tracer_exch_mpitype are read/written.
    subroutine scatter_tracers(n,tp)
        integer,intent(in) :: n
        type(tracer_particle_type),intent(in) :: tp(:)
        ! for  mpi_type_indexed() :
        integer,allocatable,dimension(:),save :: ibuf,lengths
        integer i,j,r,ierror,stat,status(MPI_STATUS_SIZE)
        integer pcoords(3,0:nprocs-1) ! cartesian topology coordinates
        real(kind=rk) :: pmin(3),pmax(3)
        integer counts(0:nprocs-1) ! # of tracers sent to each processor
        ! mpi indexed type sent to each processor except root (-> nproc-1  el.)
        integer stype(1:nprocs-1)
        ! index in  ibuf  where the chunk for each processor starts:
        integer pstart(0:nprocs-1)
        integer,parameter :: tag=0
        integer itmp,new_nlocal

        rank_0a: if (myrankc==0) then
           ! allocate  ibuf  and  lengths  on first invocation
           if (.not.allocated(ibuf)) then
              allocate (ibuf(n),lengths(n),stat=stat)
              call check_allocate(stat,'scatter_tracers(): ibuf,lengths')
              lengths(:) = 1  !  lengths(:)  remains  1  for all times
           end if

           ! allow  natoms  to alter during the simulation
           if (n>size(ibuf)) then
              deallocate (ibuf,lengths)
              allocate (ibuf(2*n),lengths(2*n),stat=stat)
              call check_allocate(stat,'scatter_tracers(): ibuf,lengths')
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
                      &,tracer_exch_mpitype,stype(r),ierror)
                 call mpi_type_commit(stype(r),ierror)
              end if
           end do rank

           if (sum(counts)/=n) then
              write (6,*) 'global number of tracers natoms=',n,', should be'&
                   &,sum(counts)
              call error_tracer('conservation of tracer count violated in '&
                   &//'scatter_tracers()')
           end if
        end if rank_0a

        call mpi_scatter(counts,1,MPI_INTEGER,new_nlocal,1,MPI_INTEGER,0&
             &,comm_cart,ierror)

        if (new_nlocal>npmax) call boost_npmax(new_nlocal)
        nlocal = new_nlocal

        rank_0b: if (myrankc==0) then
           ! these tracersare already at the correct processor
           do j=1,counts(0)
              T(j) = tp(ibuf(j)) ! pstart(0)==1 -> no explicit pstart needed
           end do

           ! rest of the tracers needs to be sent to other processors
           do r=1,nprocs-1
              call mpi_send(tp(0),1,stype(r),r,tag,comm_cart,ierror)
              call mpi_type_free(stype(r),ierror)
           end do
        else
           ! non-root processors receive tracers from root
           call mpi_recv(T,nlocal,tracer_exch_mpitype,0,tag,comm_cart,status&
                &,ierror)
        end if rank_0b

        ! build linked list for the  nlocal  particles...
        do i=1,nlocal
           itmp = atompnt
           atompnt = freepnt
           freepnt = list(freepnt)
           list(atompnt) = itmp
        end do
    end subroutine scatter_tracers

#endif
end module lbe_tracer_parallel_module
