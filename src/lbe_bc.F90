#include "lbe.h"

!> Branches off to different types of boundaries conditions (Zou/He,
!> "invasion", periodic, \c 'periodic_inflow' (in the future)).
module lbe_bc_module
!!$    use lbe_bc_periodic_inflow_module, only: bc_after_advection_periodic_inflow&
!!$         &,bc_before_advection_periodic_inflow,bc_init_periodic_inflow&
!!$         &,bc_input_periodic_inflow
    use lbe_helper_module, only: is_restoring
    use lbe_globals_module
    use lbe_log_module
    use lbe_invasion_module, only: lbe_invade,lbe_invasion_adjust_halo
    use lbe_leesedwards_module, only: le_adapt_halo_x,le_cleanup&
         &,le_halo_exchange_x,le_halo_exchange_yz,le_init,le_init_step
!#ifdef INTERPOLATEDBB 
!    use lbe_collision_simple_module, only: lbe_interpolated_dist
!#endif    
#ifndef NOSURFACTANT
    use lbe_leesedwards_module, only: le_adv_dipole_exchange
#endif
#ifdef MD
    use lbe_md_bc_leesedwards_module, only: md_before_le_adapt_halo_x
#endif
    use lbe_parallel_module, only: lbe_halo_exchange_new&
         &,lbe_setup_halo_exchange_new,tnx,tny,tnz !,tsize_pi,last_periodic_z
#ifndef NOSURFACTANT
    use lbe_parallel_module, only: lbe_adv_dipole_exchange
#endif
    use lbe_parms_module, only: boundary,boundary_cond,inv_fluid,nt,SCMP

    use lbe_timer_module, only: start_timer,stop_timer,sum_timer
    use lbe_types_module, only: lbe_site

    implicit none
    private

    public lbe_bc_after_advection,lbe_bc_before_advection&
         &,lbe_bc_complete_halo_exchange,lbe_bc_halo_exchange,lbe_bc_init&
         &,lbe_bc_init_step,lbe_bc_input,lbe_bc_le_halo_exchange&
         &,lbe_bc_shutdown,periodically_wrapped

contains

    !> supposed to contain all boundary-condition-related code to be
    !> called between LB advection and collision stage
    !>
    !> \param[in,out] lbe_N local chunk of the lattice with halo
    !> extent 1 (old LB3D style)
    !>
    !> \param[in,out] whole_N local chunk of the lattice with full
    !> halo of depth \c halo_extent
    !>
    !> \param[in,out] d_adv at the moment just passed to \c
    !> lbe_bc_complete_halo_exchange()
    subroutine lbe_bc_after_advection(lbe_N,whole_N,d_adv)
        type(lbe_site),intent(inout) :: lbe_N(0:,0:,0:)
        type(lbe_site),intent(inout) :: &
             &whole_N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
#ifndef NOSURFACTANT
        real(kind=rk),dimension(1:,1:,0:,0:,0:),intent(inout) :: d_adv
#else
        integer,intent(inout) :: d_adv
#endif

#ifdef MD
        ! [ don't know whether this is still necessary and/or correct,
        ! why is not also lbe_bc_halo_exchange() called here? ]
        if ( inv_fluid.eq.26 ) call lbe_halo_exchange_new(whole_N)
#endif

        ! [ don't know whether this is still necessary and/or correct ]
        ! important for x and y direction
        if ((inv_fluid.eq.11.or.inv_fluid.eq.23).and.SCMP) &
             &call lbe_bc_halo_exchange(whole_N)

        select case (boundary)
!!$        case ('periodic_inflow')
!!$           call bc_after_advection_periodic_inflow(lbe_N,whole_N)
        case ('periodic')
           ! nop
        case default
           call error('unknown boundary condition: boundary="'//boundary&
                &//'"')
        end select

        ! Modified by Maddalena, 19-Oct-2004
        ! do not invade if inv_fluid<0
        ! Useful to have boundary_cond.ne.0  and disabled invasion
        if (boundary_cond/=0.and.inv_fluid>=0) then
           call start_timer(ti_inv)
           DEBUG_MSG("lbe_invasion")
           call lbe_invade(lbe_N)
           call stop_timer(ti_inv)
        end if

        ! Another halo exchange
        DEBUG_MSG("lbe_bc_complete_halo_exchange")
        call lbe_bc_complete_halo_exchange(lbe_N,whole_N,d_adv)

!#ifdef INTERPOLATEDBB
!        ! compute missing reflected distributions by interpolation
!        ! on fluid node
!        DEBUG_MSG("before lbe_interpolated_dist")
!        call lbe_interpolated_dist(whole_N)
!        !call lbe_interpolated_dist(lbe_N)
!        DEBUG_MSG("after lbe_interpolated_dist")
!#endif        
    end subroutine lbe_bc_after_advection

    !> supposed to contain all boundary-condition-related code to be
    !> called between LB collision and advection stage
    !>
    !> \param[in,out] lbe_N local chunk of the lattice with halo
    !> extent 1 (old LB3D style)
    !>
    !> \param[in,out] whole_N local chunk of the lattice with full
    !> halo of depth \c halo_extent
    subroutine lbe_bc_before_advection(lbe_N,whole_N)
        type(lbe_site),intent(inout) :: lbe_N(0:,0:,0:)
        type(lbe_site),intent(inout) :: &
             &whole_N(1-halo_extent:,1-halo_extent:,1-halo_extent:)

        select case (boundary)
!!$        case ('periodic_inflow')
!!$           call bc_before_advection_periodic_inflow(whole_N)
        case ('periodic')
           ! nop
        case default
           call error('unknown boundary condition: boundary="'//boundary&
                &//'"')
        end select

        ! Halo exchange
        call start_timer(ti_halo)
        if (boundary_cond/=0.and.(inv_fluid==5.or.inv_fluid==6)) then
           call lbe_bc_le_halo_exchange(lbe_N,whole_N)
        else
           call lbe_bc_halo_exchange(whole_N)
        endif
        call stop_timer(ti_halo)
    end subroutine lbe_bc_before_advection

    !> calls the appropriate routines for halo exchange, depending on
    !> \c boundary_cond and on whether NOSURFACTANT is set
    !>
    !> \param[in,out] lbe_N local chunk of the lattice with halo
    !> extent 1 (old LB3D style)
    !>
    !> \param[in,out] whole_N local chunk of the lattice with full
    !> halo of depth \c halo_extent
    !>
    !> \param[in,out] d_adv at the moment just passed to \c
    !> lbe_adv_dipole_exchange() and \c le_adv_dipole_exchange()
    subroutine lbe_bc_complete_halo_exchange(lbe_N,whole_N,d_adv)
        type(lbe_site),intent(inout) :: lbe_N(0:,0:,0:)
        type(lbe_site),intent(inout) :: &
             &whole_N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
#ifndef NOSURFACTANT
        real(kind=rk),dimension(1:,1:,0:,0:,0:),intent(inout) :: d_adv
#else
        integer,intent(inout) :: d_adv
#endif

        call start_timer(ti_halo)
        if ((boundary_cond/=0).and.(5==inv_fluid.or.inv_fluid==6)) then
           call lbe_bc_le_halo_exchange(lbe_N,whole_N)
#ifndef NOSURFACTANT
           call le_adv_dipole_exchange(d_adv)
#endif
        else
           call lbe_bc_halo_exchange(whole_N)
#ifndef NOSURFACTANT
           call lbe_adv_dipole_exchange(d_adv)
#endif
        endif
        call stop_timer(ti_halo)
    end subroutine lbe_bc_complete_halo_exchange

    !> perform full halo exchange and possible adaption of the
    !> exchanged boundaries as required by \c lbe_invasion_module
    !>
    !> \param[in,out] whole_N local chunk of the lattice with full
    !> halo of depth \c halo_extent
    subroutine lbe_bc_halo_exchange(whole_N)
        type(lbe_site),intent(inout) :: &
             &whole_N(1-halo_extent:,1-halo_extent:,1-halo_extent:)

        call lbe_halo_exchange_new(whole_N)
        call lbe_invasion_adjust_halo(whole_N)
    end subroutine lbe_bc_halo_exchange

    !> supposed to contain all initialization code related to
    !> boundary-conditions
    !>
    !> \param[in,out] lbe_N local chunk of the lattice with halo
    !> extent 1 (old LB3D style)
    !>
    !> \param[in,out] whole_N local chunk of the lattice with full
    !> halo of depth \c halo_extent
    subroutine lbe_bc_init(lbe_N,whole_N)
        type(lbe_site),intent(inout) :: lbe_N(0:,0:,0:)
        type(lbe_site),intent(inout) :: &
             &whole_N(1-halo_extent:,1-halo_extent:,1-halo_extent:)

        ! Lees Edwards or skewed boundaries
        if (boundary_cond/=0.and.(inv_fluid==5.or.inv_fluid==6)) then
           call le_init()
           call lbe_bc_le_halo_exchange(lbe_N,whole_N)
        else
           select case (boundary)
!!$        case ('periodic_inflow')
!!$           call bc_init_periodic_inflow(whole_N)
           case ('periodic')
              call lbe_setup_halo_exchange_new(whole_N)
              call lbe_halo_exchange_new(whole_N)
!!$!!!
!!$           tsize_pi = tsize
!!$           last_periodic_z = tnz
!!$!!!
           case default
              call error('unknown boundary condition: boundary="'//boundary&
                   &//'"')
           end select
        end if
    end subroutine lbe_bc_init

    !> initialize boundary conditions for the current time step;
    !> intended to be called before all other calls in the time loop
    subroutine lbe_bc_init_step()
        if (boundary_cond/=0.and.(inv_fluid==5.or.inv_fluid==6)) &
             &call le_init_step()
    end subroutine lbe_bc_init_step

    !> read in boundary condition namelists and initialize things that
    !> are required before other initialization routines are called
    subroutine lbe_bc_input
        select case (boundary)
!!$        case ('periodic_inflow')
!!$           call bc_input_periodic_inflow
        case ('periodic')
           ! nop
        case default
           call error('unknown boundary condition: boundary="'//boundary&
                &//'"')
        end select
    end subroutine lbe_bc_input

    !> wraps everything needed for an halo update in the presence of
    !> Lees-Edwards planes (also for the non-LE halos)
    !>
    !> \param[in,out] lbe_N local chunk of the lattice with halo
    !> extent 1 (old LB3D style)
    !>
    !> \param[in] whole_N local lattice chunk with halo of extent \c
    !> halo_extent
    subroutine lbe_bc_le_halo_exchange(lbe_N,whole_N)
        type(lbe_site),intent(inout) :: lbe_N(0:,0:,0:)
        type(lbe_site),intent(in)&
             & :: whole_N(1-halo_extent:,1-halo_extent:,1-halo_extent:)

        call le_halo_exchange_x(whole_N)
#ifdef MD
        call md_before_le_adapt_halo_x(lbe_N)
#endif
        call le_adapt_halo_x(lbe_N)
        call le_halo_exchange_yz(lbe_N)
    end subroutine lbe_bc_le_halo_exchange

    !> clean-up code for LB boundary conditions
    subroutine lbe_bc_shutdown
        if ((boundary_cond.ne.0).and.(inv_fluid.eq.5.or.inv_fluid.eq.6)) then
           call le_cleanup()
        end if
    end subroutine lbe_bc_shutdown

    !> wraps a global integer lattice position around the periodic
    !> boundaries
    !>
    !> \param[in] p integer lattice position (global coordinates)
    !>
    !> \returns position \c p wrapped around into the real domain
    !> interval \c (1:tn[xyz])
    pure function periodically_wrapped(p)
        integer :: periodically_wrapped(3)
        integer,intent(in) :: p(3)
        integer :: ret(3),ts(3)

        ts = (/tnx,tny,tnz/)

        where (p>ts)
           ret = p-ts
        elsewhere (p<1)
           ret = p+ts
        elsewhere
           ret = p
        end where

        periodically_wrapped = ret
    end function periodically_wrapped

end module lbe_bc_module
