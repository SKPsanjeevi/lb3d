#include "lbe.h"

!> memory management for md part
module lbe_md_memory_module
#ifdef MD
    use lbe_globals_module, only: tsize
    use lbe_helper_module, only: is_restoring
    use lbe_log_module
    use lbe_md_globals_module
    use lbe_md_helper_module, only: error_md,log_msg_md

    use lbe_parallel_module, only: check_allocate,nprocs,tnx,tny,tnz
    use lbe_parms_module, only: n_iteration,n_restore,nt, nx,ny,nz

    implicit none
    private

    !> timestep during which npmax was changed the last time
    integer,save,public :: latest_boost_npmax

    public boost_npmax,setup_memory,shutdown_memory
contains

    !> enlarges the maximum number of owned particles  npmax  (and  P ,  list ,
    !>  nnlist , and  bin  accordingly).
    !>
    !> Owned particles in P are kept but all information about
    !> neighbor lists, bins, and halo'ed particles is lost. Also the
    !> arrangement of particles in P can change. The subroutine relies
    !> on a sane state of the global variables atompnt , list , nlocal
    !> , nomax , npmax , and P .
    subroutine boost_npmax(need)
        integer,intent(in) :: need ! minimum necessary npmax
        type(md_particle_type),allocatable,dimension(:) :: tmp_P
        integer i,ii,itmp,stat

        ! Cannot use log_msg_md due to circular dependencies
        call log_msg_md('Boosting npmax...',.true.)

        do while (npmax<need)
           npmax = npmax*2
        end do

        allocate (tmp_P(nlocal),stat=stat)
        call check_allocate(stat,'boost_npmax(): tmp_P(nlocal)')

        i = atompnt
        do ii = 1,nlocal
           tmp_P(ii) = P(i)
           i = list(i)
        end do

        deallocate (P,list,nlist,nnlist,bin)
        allocate (P(npmax+nomax),list(npmax),nlist((npmax+1)*nnmax)&
             &,nnlist(npmax+1),bin(npmax+nomax),stat=stat)
        call check_allocate(stat,'boost_npmax(): '&
             &//'P(npmax+nomax),list(npmax),nlist((npmax+1)*nnmax)'&
             &//',nnlist(npmax+1),bin(npmax+nomax)')

        P(1:nlocal) = tmp_P

        ! reset list of owned particles according to status of 0 particles...
        freepnt = 1
        do i = 1,npmax-1
           list(i) = i + 1
        enddo
        list(npmax) = 0
        atompnt = npmax + 1

        ! ...and then to  nlocal  particles.
        do i=1,nlocal
           itmp = atompnt
           atompnt = freepnt
           freepnt = list(freepnt)
           list(atompnt) = itmp
        end do

        deallocate (tmp_P)

        latest_boost_npmax = nt
    end subroutine boost_npmax

    !> calculate array dimensions and allocate arrays
    !>
    !> \note These are just estimates.
    subroutine setup_memory
        integer namax           ! maximum number of owned and halo'ed particles
        integer stat
        select case (initial_placing)
        case ('fazex')
           npmax = max(10.0_rk,(1.2_rk*tny*tnz/alat**2)/nprocs)
        case ('fazey')
           npmax = max(10.0_rk,(1.2_rk*tnx*tnz/alat**2)/nprocs)
        case ('fazez')
           npmax = max(10.0_rk,(1.2_rk*tnx*tny/alat**2)/nprocs)
        case ('fcc')
           npmax = max(10.0_rk,(2.0_rk*4*product(tsize)/alat**3)/nprocs)
        case ('lyapunov')
           npmax = max(10.0_rk,(1.2_rk*5*tnx*tny/alat**2)/nprocs)
        case ('random')
           npmax = max(10.0_rk,(1.2_rk*5*product(tsize)/alat**3)/nprocs)
        case ('sc')
           npmax = max(10.0_rk,(1.5_rk*product(tsize)/alat**3)/nprocs)
        case ('sphere')
           npmax = max(10.0_rk,(1.0_rk*np_sphere))
        case default
           npmax = max(10.0_rk,(1.2_rk*product(tsize)/alat**3)/nprocs)
        end select

        nomax = 1.5_rk*npmax * ((1+2*rs/nx)*(1+2*rs/ny)*(1+2*rs/nz) - 1)
        namax = nomax + npmax
        nemax = 2*nomax
        nfmax = nomax

        write(msgstr,"('npmax = ',I0,' , nomax = ',I0,' , namax = ',I0)") &
             &npmax, nomax, namax
        call log_msg_md(msgstr)
        write(msgstr,"('nemax = ',I0,' , nfmax = ',I0)") nemax, nfmax
        call log_msg_md(msgstr)

        if (n_stat/=0) then
           if ( is_restoring() ) then
              ! the 2nd  + 1  is for the case when the quotient is zero.
              ntmax = 1 + (n_iteration-n_restore)/n_stat + 1
           else
              ntmax = 1 + n_iteration/n_stat
           end if
        end if

        allocate (P(namax),pbuf(nfmax),stat=stat)
        call check_allocate(stat,'P,pbuf')

        ! first real allocation will be done in  borders()  - if
        !  collect_forces  is set
        allocate (ftbuf(0),stat=stat)
        call check_allocate(stat,'ftbuf')

        allocate (slens(nemax),sidxs(nemax),stat=stat)
        call check_allocate(stat,'slens,sidxs')
        slens = 1 ! remains always 1, needed for  mpi_type_indexed()

        allocate (list(npmax),nnlist(npmax+1),stat=stat)
        call check_allocate(stat,'list,nnlist')

        allocate (nlist((npmax+1)*nnmax),stat=stat)
        call check_allocate(stat,'nlist')

        allocate (bin(namax),stat=stat)
        call check_allocate(stat,'bin')

        if (n_stat/=0) then
           allocate (tmparr(ntmax),engarr(ntmax),rpotarr(ntmax),prsarr(ntmax)&
                &,conarr(ntmax),momentumarr(3,ntmax),stat=stat)
           call check_allocate(stat&
                &,'tmparr,engarr,rpotarr,prsarr,conarr,momentumarr')

           allocate (e_trans(ntmax),e_rot(ntmax),stat=stat)
           call check_allocate(stat,'e_trans,e_rot')
        end if

        latest_boost_npmax = nt
    end subroutine setup_memory

    !> deallocate arrays
    subroutine shutdown_memory
        deallocate (P,pbuf)
        deallocate (slens,sidxs)
        deallocate (list,nlist,nnlist)
        deallocate (bin)

        if (n_stat/=0) then
           deallocate (tmparr,engarr,rpotarr,prsarr,conarr,momentumarr)
           deallocate (e_trans,e_rot)
        end if

        ! allocated in  setup_parallel()
        if (ineigh==1) deallocate (binpnt)
    end subroutine shutdown_memory

#endif
end module lbe_md_memory_module
