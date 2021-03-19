#include "lbe.h"

!> main molecular dynamics module
module lbe_tracer_module
#ifdef TRACER
    use lbe_globals_module, only: border,halo_extent,maxpos,minpos,tsize,myrankc
    use lbe_parallel_module, only: comm_cart
    use lbe_timer_module, only: start_timer,stop_timer
    use lbe_tracer_globals_module
    use lbe_tracer_helper_module, only: error_tracer,fluid_velocity
    use lbe_tracer_memory_module, only: boost_npmax
    use lbe_tracer_output_module, only: count_tracer_flow,n_flow_rate
    use lbe_tracer_recoloring_module, only: recolor_tracer_at_planes
    use lbe_types_module, only: lbe_site
#ifdef MD
    use lbe_md_globals_module, only: interaction
    use lbe_tracer_md_fluid_ladd_module, only: tracer_md_fluid_ladd_interaction
#endif

    implicit none
    private

    public exchange,setup_computations,setup_list,tracer_time_step

    include 'mpif.h'

contains

    !> one \c TRACER time step
    !>
    !> \param[in] N local lattice chunk including full halo
    !>
    !> interpolate velocities, move tracers accordingly, and exchange
    !> the ones leaving the local domain
    subroutine tracer_time_step(N)
        type(lbe_site),intent(in) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        real(kind=rk),parameter :: small = 1.0e-12_rk
        real(kind=rk) :: v_old(3),x_old(3)
        integer :: i,ii,k
#ifdef MD
        logical :: interaction_ladd

        interaction_ladd = interaction=='ladd'
#endif

        call start_timer(ti_tracer_move)
        i = atompnt
        move_tracers: do ii = 1,nlocal
           x_old = T(i)%x
           v_old = T(i)%v

           ! set  T(i)%v  just to the fluid velocity at point  T(i)%x
           call fluid_velocity(N,T(i)%x,T(i)%v)

           T(i)%x = T(i)%x + 1.5_rk*T(i)%v - 0.5_rk*v_old

#ifdef MD
           ! adjust tracer positions here, so the effect of particle
           ! interactions on tracer positions will directly enter the
           ! recoloring and counting routines
           if (interaction_ladd) call tracer_md_fluid_ladd_interaction(T(i))
#endif
           call recolor_tracer_at_planes(T(i),x_old)

           if (n_flow_rate/=0) call count_tracer_flow(T(i),x_old)

           ! periodic boundary conditions (having particles beyond
           ! maxpos would confuse exchange() )
           do k=1,3
              if (T(i)%x(k) >= maxpos(k)) T(i)%x(k) = T(i)%x(k) - tsize(k)
              if (T(i)%x(k) < minpos(k)) T(i)%x(k) = T(i)%x(k) + tsize(k)

              ! enforce that all coordinates are inside intervals with
              ! boundaries of type [minpos|maxpos[ . Tracers exactly
              ! at maxpos cause problems!
              if (T(i)%x(k)==maxpos(k)) T(i)%x(k) = T(i)%x(k) - small
           end do

           i = list(i)
        end do move_tracers
        call stop_timer(ti_tracer_move)

        call start_timer(ti_tracer_comm)
        call exchange
        call stop_timer(ti_tracer_comm)
    end subroutine tracer_time_step

    !> send out tracers that have left my box, receive ones entering
    !> my box for all 3 dimensions and all respective directions
    subroutine exchange
        integer k

        DEBUG_MPI_MSG("Entered exchange")

        k_loop: do k=1,3
!!$           if (boundary=='periodic_inflow'.and.k==3) then
!!$              call exchange_dir(k,sproc(k,1),sproc(k,2)&
!!$                   &,border_pi(1,k),border_pi(2,k),.true.)
!!$              call exchange_dir(k,sproc(k,3),sproc(k,4)&
!!$                   &,border_npi(1,k),border_npi(2,k),.false.)
!!$           else
              call exchange_dir(k,sproc(k,1),sproc(k,2),border(1,k),border(2,k))
!!$           end if
        end do k_loop

        DEBUG_MPI_MSG("Returning from exchange")
    end subroutine exchange

    !> actually send and recv tracers in a given dimension
    !>
    !> \param[in] dim dimension
    !>
    !> \param[in] parlo rankc of lower communication partner
    !>
    !> \param[in] parhi rankc of upper communication partner
    !>
    !> \param[in] bndlo lower chunk boundary
    !>
    !> \param[in] bndhi upper chunk boundary
    !>
    !> \param[in] periodic (optional) send only particles \c p with \c
    !> periodic==is_periodic(p), if not specified all particles are sent
    !>
    !> \note The check for \c freepnt!=0 should not be necessary...
    subroutine exchange_dir(curdim,parlo,parhi,bndlo,bndhi,periodic)
        integer,intent(in) :: curdim,parlo,parhi
        real(kind=rk),intent(in) :: bndlo,bndhi
        real(kind=rk) :: blo,bhi
        logical,intent(in),optional :: periodic
        integer :: actrcount,i,ierror,ii,iprev,itmp,j,maxrcount,ndelete&
             &,totrcount,status(MPI_STATUS_SIZE),stype
        integer,parameter :: tag=0
        logical ignore_per,per

        if (parlo==MPI_PROC_NULL.and.parhi==MPI_PROC_NULL) return

        ! don't sort out tracers that no cpu would take (Actually,
        ! one might think about deleting such particles... !!!)
        if (parlo==MPI_PROC_NULL) then
           blo = -huge(blo)
        else
           blo = bndlo
        end if
        if (parhi==MPI_PROC_NULL) then
           bhi = huge(bhi)
        else
           bhi = bndhi
        end if

        if (present(periodic)) then
           ignore_per = .false.
           per = periodic
        else
           ignore_per = .true.
           per = .true.         ! set this just to please debugging tools
        end if

        ndelete = 0
        iprev = 0
        j = 0
        i = atompnt

        ! index tracers leaving my box, update local list
        do ii = 1,nlocal
!!$           if ((ignore_per.or.(is_periodic(T(i)).eqv.per))&
           if ((ignore_per)&
                &.and.(T(i)%x(curdim)<blo.or.T(i)%x(curdim)>=bhi)) then
              j = j + 1
              if (j<=nfmax) sidxs(j) = i

              ndelete = ndelete + 1
              if (iprev.eq.0) then
                 atompnt = list(i)
              else
                 list(iprev) = list(i)
              endif
              itmp = list(i)
              list(i) = freepnt
              freepnt = i
              i = itmp
           else
              iprev = i
              i = list(i)
           endif
        enddo
        nlocal = nlocal-ndelete

        if (j>nfmax) then
           write (6,*) myrankc,curdim,blo,bhi,'j=',j,',nfmax=',nfmax
           call error_tracer('Sending too many exchange tracers. '&
                &//'(Boost nfmax to values>j !)')
        endif

        call mpi_type_indexed(j,slens,sidxs,tracer_exch_mpitype,stype,ierror)
        call mpi_type_commit(stype,ierror)

        ! send them out in both directions (if neighboring nodes are
        ! different)
        totrcount = 0
        maxrcount = nfmax
        call mpi_sendrecv(T(0),1,stype,parlo,tag&
             &,tbuf(1),maxrcount,tracer_exch_mpitype,parhi,tag&
             &,comm_cart,status,ierror)
        call mpi_get_count(status,tracer_exch_mpitype,actrcount,ierror)
        totrcount = totrcount + actrcount

        if (parlo/=parhi) then
           maxrcount = nfmax - totrcount
           call mpi_sendrecv(T(0),1,stype,parhi,tag&
                &,tbuf(totrcount+1),maxrcount,tracer_exch_mpitype,parlo,tag&
                &,comm_cart,status,ierror)
           call mpi_get_count(status,tracer_exch_mpitype,actrcount,ierror)
           totrcount = totrcount + actrcount
        end if

        call mpi_type_free(stype,ierror)

        ! check incoming tracers to see if they are in my box (could be
        ! in node's box on other side of the sender)
        sort_incoming: do j = 1,totrcount
           if (tbuf(j)%x(curdim)>=blo.and.tbuf(j)%x(curdim)<bhi) then
              if (nlocal+1>npmax) call boost_npmax(nlocal+1)
              nlocal = nlocal + 1
              ! The case freepnt==0 should never happen thanks
              ! to boost_npmax()
              if (freepnt.ne.0) then
                 itmp = atompnt
                 atompnt = freepnt
                 freepnt = list(freepnt)
                 list(atompnt) = itmp
                 ! assign only elements from  tracer_exch_mpitype ???
                 T(atompnt) = tbuf(j)
              endif
           endif
        enddo sort_incoming
    end subroutine exchange_dir

    !> initialize constants that simplify the computations
    subroutine setup_computations
    end subroutine setup_computations

    !> prepare list structure
    subroutine setup_list
        integer i

        freepnt = 1
        do i = 1,npmax
           list(i) = i + 1
        enddo
        list(npmax) = 0
        nlocal = 0
        atompnt = npmax + 1
    end subroutine setup_list

#endif
end module lbe_tracer_module
