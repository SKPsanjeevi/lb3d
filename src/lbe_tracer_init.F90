#include "lbe.h"
!> initialization of \c TRACER part and initial condition for \c
!> TRACER particles
module lbe_tracer_init_module
#ifdef TRACER
  use lbe_globals_module, only: halo_extent,minpos,tsize
  use lbe_helper_module, only: is_fluid
  use lbe_log_module
  use lbe_parallel_module, only: calculate_displacements,comm_cart&
       &,get_send_partner,nprocs,start,tnx,tny
  use lbe_parms_module, only: nx,ny,nz
  use lbe_timer_module, only: register_timer
  use lbe_tracer_globals_module
  use lbe_tracer_helper_module, only: build_tracer_mpitype,fluid_velocity&
       &,log_msg_tracer
  use lbe_tracer_memory_module, only: boost_npmax
  use lbe_tracer_module, only: setup_list
  use lbe_types_module, only: lbe_site

  implicit none
  include 'mpif.h'
  private

  public place_tracers,setup_general,setup_parallel

  integer,save,public :: t_placement=0 !< time step to create tracers

  !> place a tracer only on every \c place_mod_pos'th lattice site
  integer,save,public :: place_mod_pos=1

contains

    !> initialize \c kind element of all tracers
    subroutine initialize_kinds
        integer i,ii

        i = atompnt
        do ii=1,nlocal
           ! all tracers in the first half of the system in x
           ! direction have kind 1, all others have kind 2
           if (T(i)%x(1)<minpos(1)+0.5_rk*tsize(1)) then
              T(i)%kind = 1
           else
              T(i)%kind = 2
           end if

           i = list(i)
        end do
    end subroutine initialize_kinds

    !> general initialization of \c TRACER part
    subroutine setup_general
        call register_timer('TRACER:Comm',ti_tracer_comm)
        call register_timer('TRACER:Dump',ti_tracer_dump)
        call register_timer('TRACER:Move',ti_tracer_move)
    end subroutine setup_general

    !> setup communication patterns
    subroutine setup_parallel
        integer i,k

        call setup_mpitypes

        do k=1,3                ! dimension
           do i=1,4             ! (pseudo-)direction
              sproc(k,i) = get_send_partner(k,i)
           end do
        end do
    end subroutine setup_parallel

    !> create mpi types that represent a single tracer in MPI
    !> communication
    subroutine setup_mpitypes
        call build_tracer_mpitype(tracer_exch_mpitype,x=.true.,v=.true.&
             &,akind=.true.,uid=.true.)
    end subroutine setup_mpitypes

    !> initialize tracer particles in the system
    !>
    !> \param[in] N local lattice chunk including full halo
    subroutine place_tracers(N)
        type(lbe_site),intent(in) :: N(&
             &1-halo_extent:,1-halo_extent:,1-halo_extent:)
        integer counts(0:nprocs-1),first_uids(0:nprocs-1),first_uid
        integer :: i,ii,ierror

        call setup_list()
        call place_tracers_on_lattice(N)

        ! initialize uids
        call MPI_Gather(nlocal,1,MPI_INTEGER,counts,1,MPI_INTEGER,0&
             &,comm_cart,ierror)
        call calculate_displacements(counts,first_uids)
        call MPI_Scatter(first_uids,1,MPI_INTEGER,first_uid,1,MPI_INTEGER,0&
             &,comm_cart,ierror)
        i = atompnt
        do ii=1,nlocal
           T(i)%uid = first_uid+ii
           i = list(i)
        end do

        call initialize_kinds()

        write (msgstr,"('Placed a total of ',I0,' tracers.')") sum(counts)
        call log_msg_tracer(msgstr)

        ! initialize velocities
        i = atompnt
        do ii = 1,nlocal
           call fluid_velocity(N,T(i)%x,T(i)%v)
           i = list(i)
        end do

        call log_msg_tracer("Finished tracer creation.")
    end subroutine place_tracers

    !> place a single tracer at position \c pos
    subroutine place_tracer(pos)
        real(kind=rk),intent(in) :: pos(3)
        integer itmp,k

        if (nlocal+1>npmax) call boost_npmax(nlocal+1)
        nlocal = nlocal + 1

        itmp = atompnt
        atompnt = freepnt
        freepnt = list(freepnt)
        list(atompnt) = itmp

        T(atompnt)%x = pos
    end subroutine place_tracer

    !> fills the local domain with tracers, one tracer exactly at each
    !> fluid site
    !>
    !> \param[in] N local lattice chunk (halo depth >=1 required)
    subroutine place_tracers_on_lattice(N)
        type(lbe_site),intent(in) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        integer :: x,y,z
        ! this should suffice to count all lattice sites
        integer,parameter :: ik=selected_int_kind(12)
        integer(kind=ik) :: gx(3),pos

        do x=1,nx
           do y=1,ny
              do z=1,nz
                 gx = (/x,y,z/)+start-2
                 pos = gx(1)+gx(2)*tnx+gx(3)*tnx*tny
                 if (mod(pos,int(place_mod_pos,kind=ik))==0&
                      &.and.is_fluid(N(x,y,z)%rock_state)) &
                      &call place_tracer(real((/x,y,z/)+start-1,kind=rk))
              end do
           end do
        end do
    end subroutine place_tracers_on_lattice

#endif
end module lbe_tracer_init_module
