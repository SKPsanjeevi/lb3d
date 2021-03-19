#include "lbe.h"

!> memory management for \c TRACER part
module lbe_tracer_memory_module
#ifdef TRACER
    use lbe_globals_module, only: tsize
    use lbe_helper_module, only: is_restoring
    use lbe_log_module
    use lbe_parallel_module, only: check_allocate,nprocs,tnx,tny,tnz
    use lbe_parms_module, only: n_iteration,n_restore,nt,nx,ny,nz
    use lbe_tracer_globals_module
    use lbe_tracer_helper_module, only: log_msg_tracer

    implicit none
    private

    !> timestep during which npmax was changed the last time
    integer,save,public :: latest_boost_npmax

    public boost_npmax,setup_memory,shutdown_memory
contains

    !> enlarges the maximum number of owned tracers npmax (and \c T
    !> and \c list accordingly).
    !>
    !> Owned tracers in T are kept. Also the arrangement of tracers in
    !> T can change. The subroutine relies on a sane state of the
    !> global variables \c atompnt , \c list , \c nlocal , \c npmax ,
    !> and \c T .
    subroutine boost_npmax(need)
        integer,intent(in) :: need ! minimum necessary npmax
        type(tracer_particle_type),allocatable,dimension(:) :: tmp_T
        integer i,ii,itmp,stat

        ! Cannot use log_msg_tracer due to circular dependencies
        call log_msg_tracer('Boosting npmax...',.true.)

        do while (npmax<need)
           npmax = npmax*2
        end do

        allocate (tmp_T(nlocal),stat=stat)
        call check_allocate(stat,'boost_npmax(): tmp_T(nlocal)')

        i = atompnt
        do ii = 1,nlocal
           tmp_T(ii) = T(i)
           i = list(i)
        end do

        deallocate (T,list)
        allocate (T(npmax),list(npmax),stat=stat)
        call check_allocate(stat,'boost_npmax(): T(npmax),list(npmax)')

        T(1:nlocal) = tmp_T

        ! reset list of owned tracers according to status of 0 tracers...
        freepnt = 1
        do i = 1,npmax-1
           list(i) = i + 1
        enddo
        list(npmax) = 0
        atompnt = npmax + 1

        ! ...and then to  nlocal  tracers.
        do i=1,nlocal
           itmp = atompnt
           atompnt = freepnt
           freepnt = list(freepnt)
           list(atompnt) = itmp
        end do

        deallocate (tmp_T)

        latest_boost_npmax = nt
    end subroutine boost_npmax

    !> calculate array dimensions and allocate arrays
    !>
    !> \note These are just estimates.
    subroutine setup_memory
        integer stat

        npmax = 2*product((/nx,ny,nz/))
        nfmax = 2*2*(nx*ny+nx*nz+ny*nz)

        write(msgstr,"('npmax = ',I0)") npmax
        call log_msg_tracer(msgstr)
        write(msgstr,"('nfmax = ',I0)") nfmax
        call log_msg_tracer(msgstr)

        allocate (T(npmax),tbuf(nfmax),stat=stat)
        call check_allocate(stat,'T,tbuf')

        allocate (slens(nfmax),sidxs(nfmax),stat=stat)
        call check_allocate(stat,'slens,sidxs')
        slens = 1 ! remains always 1, needed for  mpi_type_indexed()

        allocate (list(npmax),stat=stat)
        call check_allocate(stat,'list')

        latest_boost_npmax = nt
    end subroutine setup_memory

    !> deallocate arrays
    subroutine shutdown_memory
        deallocate (T,tbuf)
        deallocate (slens,sidxs)
        deallocate (list)
    end subroutine shutdown_memory

#endif
end module lbe_tracer_memory_module
