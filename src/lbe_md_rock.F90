#include "lbe.h"

!> Delegates to different particle-wall potential implementations
module lbe_md_rock_module
#ifdef MD

    use lbe_globals_module, only: cx,cy,cz,halo_extent,nnonrest,myrankc,rock_value_surface
    use lbe_helper_module, only: is_restoring
    use lbe_md_globals_module
    use lbe_md_helper_module, only: count_particles_all,scatter_rock_state_halo&
         &,log_msg_md,error_md
    use lbe_md_rock_lj_module, only: input_rock_lj,local_rock_potential_lj&
         &,particle_rock_potential_lj,rock_f_interaction_lj,setup_rock_lj
    use lbe_md_rock_dlvo_module, only: input_rock_dlvo,local_rock_potential_dlvo&
         &,particle_rock_potential_dlvo,rock_f_interaction_dlvo,setup_rock_dlvo
    use lbe_md_rock_bprh_module, only: input_rock_bprh&
         &,local_rock_potential_bprh,particle_rock_potential_bprh&
         &,rock_ft_interaction_bprh,setup_rock_bprh
    use lbe_md_rock_friction_module, only: input_rock_friction&
         &,rock_f_interaction_friction,setup_rock_friction
    use lbe_parallel_module, only: check_allocate,comm_cart,gather_rock_state&
         &,tnx,tny,tnz
    use lbe_parms_module, only: nx,ny,nz
    use lbe_types_module, only: lbe_site

#ifdef DEBUG_MPI
    use lbe_parallel_module, only: log_mpi
#endif

    implicit none
    private

    include 'mpif.h'

    public input_rock,particle_rock_potential,rock_ft_interaction&
         &,rock_halo_extent,rock_potential,setup_rock

    character(len=32),save,public :: rock='none' !< particle-rock interaction
    character(len=32),save,public :: rock_friction='none' !< particle-rock interaction
contains

    !> read particle-rock coupling specific namelist from input file
    subroutine input_rock
        select case (rock)
        case ('bprh')
           call input_rock_bprh
        case ('lj')
           call input_rock_lj
        case ('dlvo')
           call input_rock_dlvo
        case ('none')
           ! nop
        case default
           call error_md('unknown value: rock="'//rock//'"')
        end select
        !> read md-rock-friction coupling
        select case (rock_friction)
        case ('stokes')
           call input_rock_friction
        case ('none')
           ! nop
        case default
           call error_md('unknown value: rock_friction="'//rock_friction//'"')
        end select
    end subroutine input_rock

    !> Returns the energy due to the rock potential according to the
    !> global rock state array \c rock_state for a single particle at
    !> position \c pos with orientation \c ori .
    real(kind=rk) function particle_rock_potential(rock_state,pos,ori)
        real(kind=rk),intent(in) :: rock_state&
             &(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        real(kind=rk),intent(in) :: pos(3),ori(0:3)

        select case (rock)
        case ('bprh')
           particle_rock_potential&
                & = particle_rock_potential_bprh(rock_state,pos,ori)
        case ('lj')
           particle_rock_potential&
                & = particle_rock_potential_lj(rock_state,pos)
        case ('dlvo')
           particle_rock_potential&
                & = particle_rock_potential_dlvo(rock_state,pos)
        case ('none')
           particle_rock_potential = 0.0_rk
        case default
           call error_md('unknown value: rock="'//rock//'"')
        end select
    end function particle_rock_potential

    subroutine rock_ft_interaction(N)
        type(lbe_site),intent(inout) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)

        DEBUG_MPI_MSG("Starting rock interaction")

        select case (rock)
        case ('bprh')
           call rock_ft_interaction_bprh(N)
        case ('lj')
           call rock_f_interaction_lj(N)
        case ('dlvo')
           call rock_f_interaction_dlvo(N)
        case ('none')
           ! nop
        case default
           call error_md('unknown value: rock="'//rock//'"')
        end select
  
        select case (rock_friction)
        case ('stokes')
           call rock_f_interaction_friction(N)
        case ('none')
           ! nop
        case default
           call error_md('unknown value: rock_friction="'//rock_friction//'"')
        end select
        DEBUG_MPI_MSG("Finished rock interaction")
    end subroutine rock_ft_interaction

    !> enlarge \c halo_extent if necessary for particle-rock interaction
    subroutine rock_halo_extent
        select case (rock)
        case ('bprh')
           ! nop
        case ('lj')
           ! nop
        case ('dlvo')
           ! nop
        case ('none')
           ! nop
        case default
           call error_md('unknown value: rock="'//rock//'"')
        end select
    end subroutine rock_halo_extent

    !> average particle-rock potential per particle
    subroutine rock_potential(N,rpot)
        type(lbe_site),intent(inout) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        real(kind=rk),intent(out) :: rpot
        real(kind=rk) rrpot
        integer ierror,n_global

        select case (rock)
        case ('bprh')
           call local_rock_potential_bprh(N,rpot)
        case ('lj')
           call local_rock_potential_lj(N,rpot)
        case ('dlvo')
           call local_rock_potential_dlvo(N,rpot)
        case ('none')
           ! nop
        case default
           call error_md('unknown value: rock="'//rock//'"')
        end select

        call mpi_allreduce(rpot,rrpot,1,MPI_REAL8,MPI_SUM,comm_cart,ierror)
        call count_particles_all(n_global)
        rpot = rrpot/n_global
    end subroutine rock_potential

    subroutine setup_rock(N)
        type(lbe_site),intent(inout),target :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
!!$        type(lbe_site),pointer :: lbe_N(:,:,:) ! the same as  N  in lbe.F90

        ! after checkpoint restoration everything is just fine already
        no_restore: if ( .not. is_restoring() ) then
           ! At this point, the lbe halo was already exchanged, which
           ! is enough to find rock sites at the surface of a wall.
           call mark_surface_rock(N)
!!$           ! The new rock_state values are exchanged into the md halo
!!$           ! regions later. The exchange into the lbe halos is done
!!$           ! here:
!!$           lbe_N => N(0:nx+1,0:ny+1,0:nz+1)
!!$           call lbe_halo_exchange(lbe_N)
        end if no_restore

        select case (rock)
        case ('bprh')
           call setup_rock_bprh
        case ('lj')
           call setup_rock_lj
        case ('dlvo')
           call setup_rock_dlvo
        case ('none')
           ! nop
        case default
           call error_md('unknown value: rock="'//rock//'"')
        end select
        select case (rock_friction)
        case ('stokes')
           call setup_rock_friction
        case ('none')
           ! nop
        case default
           call error_md('unknown value: rock="'//rock//'"')
        end select
    end subroutine setup_rock

    subroutine mark_surface_rock(N)
      type(lbe_site),intent(inout) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
      integer x,y,z

      call log_msg_md("Marking surface rock nodes.")
      do x=1,nx
        do y=1,ny
          do z=1,nz
            if (N(x,y,z)%rock_state/=0.0_rk.and.any(&
                 &N(x+cx(1:nnonrest),y+cy(1:nnonrest),z+cz(1:nnonrest))%&
                 &rock_state==0.0_rk))&
                 & N(x,y,z)%rock_state = rock_value_surface
          end do
        end do
      end do
    end subroutine mark_surface_rock

#endif
end module lbe_md_rock_module
