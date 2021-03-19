#include "lbe.h"

!> delegates to different md particle-lb fluid interaction
!> implementations
module lbe_md_fluid_module
#ifdef MD
    use lbe_globals_module, only: force_halo_extent,halo_extent
    use lbe_helper_module, only: density,massflow,velocity,request_halo_extent
    use lbe_md_fluid_friction_module, only: fluid_f_interaction_friction&
         &,input_fluid_friction,setup_fluid_friction
    use lbe_md_fluid_ladd_module, only: fluid_ft_interaction_ladd_pre&
         &,fluid_ft_interaction_ladd_post,fluid_ft_interaction_ladd_middle&
         &,init_fluid_ladd,input_fluid_ladd,local_particle_velocity&
         &,setup_fluid_ladd,shutdown_fluid_ladd,summary_fluid_ladd
    use lbe_md_fluid_tracer_module, only: fluid_v_interaction_tracer&
         &,input_fluid_tracer,setup_fluid_tracer,summary_fluid_tracer
    use lbe_md_globals_module
    use lbe_md_helper_module, only: error_md
    use lbe_types_module, only: lbe_site

#ifdef DEBUG_MPI
    use lbe_parallel_module, only: log_mpi
#endif

    implicit none
    private

    public fluid_f_interaction_pre,fluid_ft_interaction_pre&
         &,fluid_f_interaction_post,fluid_ft_interaction_post&
         &,fluid_force_halo_extent,fluid_halo_extent,fluid_post,fluid_middle&
         &,fluid_v_interaction,fluid_w_interaction,init_fluid,input_fluid&
         &,md_density,md_massflow,md_velocity,setup_fluid,shutdown_fluid&
         &,summary_fluid

    include 'mpif.h'

contains

    subroutine fluid_f_interaction_pre(N)
        type(lbe_site),intent(inout) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)

        DEBUG_MPI_MSG("Starting fluid_f_interaction_pre()")

        select case (interaction)
        case ('friction')
           call fluid_f_interaction_friction(N)
        case ('ladd')
           ! nop (everything is done in  fluid_ft_interaction_ladd_pre)
        case ('ladd_janus')
           ! nop (everything is done in  fluid_ft_interaction_ladd_pre)
        case ('none')
           ! nop
        case ('tracer')
           ! nop
        case default
           call error_md('unknown value: interaction="'//interaction//'"')
        end select

        DEBUG_MPI_MSG("Finished fluid_f_interaction_pre()")
    end subroutine fluid_f_interaction_pre

    !> branches off to particle-fluid coupling routines that may
    !> change particle forces and torques (and the fluid, maybe)
    !>
    !> \param[in,out] N local lattice chunk with halo of extent \c
    !> halo_extent
    !>
    !> \param[in] substep current md substep count (within \c
    !> 1..steps_per_lbe_step)
    subroutine fluid_ft_interaction_pre(N,substep)
        type(lbe_site),intent(inout) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        integer,intent(in) :: substep
        integer :: i,ii

        DEBUG_MPI_MSG("Starting fluid_ft_interaction_pre()")

        select case (interaction)
        case ('friction')
           ! nop
        case ('ladd')
           call fluid_ft_interaction_ladd_pre(N,substep)
        case ('ladd_janus')
           call fluid_ft_interaction_ladd_pre(N,substep)
        case ('none')
           ! nop
        case ('tracer')
           ! nop
        case default
           call error_md('unknown value: interaction="'//interaction//'"')
        end select

        DEBUG_MPI_MSG("Finished fluid_ft_interaction_pre()")
    end subroutine fluid_ft_interaction_pre

    subroutine fluid_f_interaction_post(N)
        type(lbe_site),intent(inout) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)

        select case (interaction)
        case ('friction')
           ! nop
        case ('ladd')
           ! nop (everything is done in  fluid_ft_interaction_ladd_post)
        case ('ladd_janus')
           ! nop (everything is done in  fluid_ft_interaction_ladd_post)
        case ('none')
           ! nop
        case ('tracer')
           ! nop
        case default
           call error_md('unknown value: interaction="'//interaction//'"')
        end select
    end subroutine fluid_f_interaction_post

    !> branches off to particle-fluid coupling routines that may
    !> change particle forces and torques (and the fluid, maybe)
    !>
    !> \param[in,out] N local lattice chunk with halo of extent \c
    !> halo_extent
    !>
    !> \param[in] substep current md substep count (within \c
    !> 1..steps_per_lbe_step)
    subroutine fluid_ft_interaction_post(N,substep)
        type(lbe_site),intent(inout) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        integer,intent(in) :: substep
        integer i,ii

        select case (interaction)
        case ('friction')
           ! nop
        case ('ladd')
           call fluid_ft_interaction_ladd_post(N,substep)
        case ('ladd_janus')
           call fluid_ft_interaction_ladd_post(N,substep)
        case ('none')
           ! nop
        case ('tracer')
           ! nop
        case default
           call error_md('unknown value: interaction="'//interaction//'"')
        end select

        ft_fluid_defined: if (use_ft_fluid) then
           which_substep: if (substep==1) then
              ! no substepping or 1st substep: apply average force/torque
              ! stored from the previous and the current LB step. The
              ! averaging helps to suppress possible lattice
              ! oscillations. Add current force/torque for owned and
              ! halo'ed particles, a collect() will follow that sums up
              ! f, t, f_fluid, and t_fluid onto the owning
              ! nodes. [ft]_fluid_prev is present only and complete
              ! already for the owner, so it is added only for local
              ! particles.
              i = atompnt
              if (average_ft_fluid) then
                 do ii = 1,nlocal
                    P(i)%f = P(i)%f + 0.5_rk*(P(i)%f_fluid+P(i)%f_fluid_prev)
                    P(i)%t = P(i)%t + 0.5_rk*(P(i)%t_fluid+P(i)%t_fluid_prev)
#ifdef FORCECOMPONENT                    
                    P(i)%f_n = P(i)%f_n + P(i)%f_normal
                    P(i)%f_t = P(i)%f_t + P(i)%f_tangent
#endif                    

                    i = list(i)
                 end do
                 do ii = 1,nother
                    P(i)%f = P(i)%f + 0.5_rk*P(i)%f_fluid
                    P(i)%t = P(i)%t + 0.5_rk*P(i)%t_fluid
#ifdef FORCECOMPONENT                    
                    P(i)%f_n = P(i)%f_n + P(i)%f_normal
                    P(i)%f_t = P(i)%f_t + P(i)%f_tangent
#endif                    

                    i = i+1
                 end do
              else
                 do ii = 1,nlocal
                    P(i)%f = P(i)%f + P(i)%f_fluid
                    P(i)%t = P(i)%t + P(i)%t_fluid
#ifdef FORCECOMPONENT                    
                    P(i)%f_n = P(i)%f_n + P(i)%f_normal
                    P(i)%f_t = P(i)%f_t + P(i)%f_tangent
#endif                    

                    !write(*, "('i : ' I5)") i
                    i = list(i)
                 end do
                 do ii = 1,nother
                    P(i)%f = P(i)%f + P(i)%f_fluid
                    P(i)%t = P(i)%t + P(i)%t_fluid
#ifdef FORCECOMPONENT                    
                    P(i)%f_n = P(i)%f_n + P(i)%f_normal
                    P(i)%f_t = P(i)%f_t + P(i)%f_tangent
#endif                    
                    !write(*, "('i : ' I5)") i
                    i = i+1
                 end do
              end if
           else which_substep
              ! later substep: apply average force/torque stored from the
              ! previous and the current LB step. Do this only for owned
              ! particles since f_fluid and t_fluid were collected
              ! already and will not be touched by the substep collect()
              ! operation. [ft]_fluid_prev is present and complete only
              ! for the owner, too. Owner also accumulates substep
              ! velocities which will---after averaging---be the v_fluid
              ! and ws_fluid of the next LB step and will then be
              ! averaged into v_fluid_avg and ws_fluid_avg and spread in
              ! communicate().
              i = atompnt
              do ii = 1,nlocal
                 if (average_ft_fluid) then
                    P(i)%f = P(i)%f + 0.5_rk*(P(i)%f_fluid+P(i)%f_fluid_prev)
                    P(i)%t = P(i)%t + 0.5_rk*(P(i)%t_fluid+P(i)%t_fluid_prev)
                 else
                    P(i)%f = P(i)%f + P(i)%f_fluid
                    P(i)%t = P(i)%t + P(i)%t_fluid
                 end if

                 P(i)%v_fluid_acc = P(i)%v_fluid_acc + P(i)%v
                 P(i)%ws_fluid_acc = P(i)%ws_fluid_acc + P(i)%ws

                 i = list(i)
              end do
           end if which_substep
	end if ft_fluid_defined
    end subroutine fluid_ft_interaction_post

    !> enlarge \c halo_extent if necessary for particle-fluid interaction
    subroutine fluid_halo_extent
        select case (interaction)
        case ('friction')
          call request_halo_extent(2,"MD friction")
        case ('ladd')
          call request_halo_extent(1,"MD Ladd")
        case ('ladd_janus')!keep?
          call request_halo_extent(1,"MD Ladd")
        case ('none')
           ! nop
        case ('tracer')
          call request_halo_extent(2,"MD tracer")
        case default
           call error_md('unknown value: interaction="'//interaction//'"')
        end select
    end subroutine fluid_halo_extent

    !> enlarge \c force_halo_extent if necessary for particle-fluid interaction
    subroutine fluid_force_halo_extent
        select case (interaction)
        case ('friction')
           force_halo_extent = max(force_halo_extent,1)
        case ('ladd')
           ! nop
        case ('ladd_janus')
           ! nop
        case ('none')
           ! nop
        case ('tracer')
           ! nop
        case default
           call error_md('unknown value: interaction="'//interaction//'"')
        end select
    end subroutine fluid_force_halo_extent

    !> called between lbe streaming and collision
    subroutine fluid_middle(N)
        type(lbe_site),intent(inout) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)

        select case (interaction)
        case ('friction')
           ! nop
        case ('ladd')
           call fluid_ft_interaction_ladd_middle(N)
        case ('ladd_janus')
           call fluid_ft_interaction_ladd_middle(N)
        case ('none')
           ! nop
        case ('tracer')
           ! nop
        case default
           call error_md('unknown value: interaction="'//interaction//'"')
        end select
    end subroutine fluid_middle

    !> called after lbe streaming and collision
    subroutine fluid_post(N)
        type(lbe_site),intent(inout) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)

        select case (interaction)
        case ('friction')
           ! nop
        case ('ladd')
           ! nop
        case ('ladd_janus')
           ! nop
        case ('none')
           ! nop
        case ('tracer')
           ! nop
        case default
           call error_md('unknown value: interaction="'//interaction//'"')
        end select
    end subroutine fluid_post

    subroutine fluid_v_interaction(N)
        type(lbe_site),intent(inout) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)

        select case (interaction)
        case ('friction')
           ! nop
        case ('ladd')
           ! nop
        case ('ladd_janus')
           ! nop
        case ('none')
           ! nop
        case ('tracer')
           call fluid_v_interaction_tracer(N)
        case default
           call error_md('unknown value: interaction="'//interaction//'"')
        end select
    end subroutine fluid_v_interaction

    subroutine fluid_w_interaction(N)
        type(lbe_site),intent(inout) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)

        select case (interaction)
        case ('friction')
           ! nop
        case ('ladd')
           ! nop
        case ('ladd_janus')
           ! nop
        case ('none')
           ! nop
        case ('tracer')
           ! nop
        case default
           call error_md('unknown value: interaction="'//interaction//'"')
        end select
    end subroutine fluid_w_interaction

    !> read particle-fluid coupling specific namelist from input file
    subroutine input_fluid
        select case (interaction)
        case ('friction')
           call input_fluid_friction
        case ('ladd')
           call input_fluid_ladd
        !case ('ladd_janus')
           !call input_fluid_ladd_janus
        case ('none')
           ! nop
        case ('tracer')
           call input_fluid_tracer
        case default
           call error_md('unknown value: interaction="'//interaction//'"')
        end select
    end subroutine input_fluid

    !> do things that have to be done after the particle positions are known
    subroutine init_fluid(N)
        type(lbe_site),intent(inout) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)

        select case (interaction)
        case ('friction')
           ! nop
        case ('ladd')
           call init_fluid_ladd(N)
        !case ('ladd_janus')
           !call init_fluid_ladd_janus(N)
        case ('none')
           ! nop
        case ('tracer')
           ! nop
        case default
           call error_md('unknown value: interaction="'//interaction//'"')
        end select
    end subroutine init_fluid

    function md_density(N,x,y,z,rhof)
        real(kind=rk) :: md_density
        type(lbe_site),intent(in) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        integer,intent(in) :: x,y,z
        real(kind=rk),intent(in) :: rhof
        integer :: rs

        rs = nint(N(x,y,z)%rock_state)
        if (rs==0) then
           md_density = density(N(x,y,z))
        else if (interaction=='ladd'.and.rs>0) then
           md_density = rhof
        else
           md_density = 0.0_rk
        end if
    end function md_density

    !> calculates the massflow at a given site while taking into
    !> account possible md particles there
    !>
    !> If the site belongs to a \c ladd particle, the mass flow is
    !> calculated based on the local particle velocity according on
    !> its rigid body motion and an assumed fluid mass density
    !> specified as argument \c rhof. This requires \c uid2i to be set
    !> up!
    !>
    !> \param[in] N local lattice chunk with full halo
    !>
    !> \param[in] x local lattice position in x direction
    !>
    !> \param[in] y local lattice position in x direction
    !>
    !> \param[in] z local lattice position in x direction
    !>
    !> \param[in] rhof assumed fluid mass density at \c ladd particle
    !> sites
    !>
    !> \returns mass flow
    function md_massflow(N,x,y,z,rhof)
        real(kind=rk) :: md_massflow(3)
        type(lbe_site),intent(in) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        integer,intent(in) :: x,y,z
        real(kind=rk),intent(in) :: rhof
        integer :: rs

        rs = nint(N(x,y,z)%rock_state)
        if (rs==0) then
           md_massflow = massflow(N(x,y,z))
        else if (interaction=='ladd'.and.rs>0) then
           md_massflow = rhof&
                &*local_particle_velocity(real((/x,y,z/),kind=rk),rs)
        else
           md_massflow = 0.0_rk
        end if
    end function md_massflow

    !> calculates the velocity at a given site while taking into
    !> account possible md particles there
    !>
    !> If the site belongs to a \c ladd particle, the velocity is
    !> calculated based on the local particle velocity according on
    !> its rigid body motion. This requires \c uid2i to be set up!
    !>
    !> \param[in] N local lattice chunk with full halo
    !>
    !> \param[in] x local lattice position in x direction
    !>
    !> \param[in] y local lattice position in x direction
    !>
    !> \param[in] z local lattice position in x direction
    !>
    !> \returns velocity
    function md_velocity(N,x,y,z)
        real(kind=rk) :: md_velocity(3)
        type(lbe_site),intent(in) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        integer,intent(in) :: x,y,z
        integer :: rs

        rs = nint(N(x,y,z)%rock_state)
        if (rs==0) then
           md_velocity = velocity(N(x,y,z))
        else if (interaction=='ladd'.and.rs>0) then
           md_velocity = local_particle_velocity(real((/x,y,z/),kind=rk),rs)
        else
           md_velocity = 0.0_rk
        end if
    end function md_velocity

    !> setup fluid coupling
    subroutine setup_fluid(N)
        type(lbe_site),intent(inout) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)

        select case (interaction)
        case ('friction')
           call setup_fluid_friction
        case ('ladd')
           call setup_fluid_ladd(N)
        !case ('ladd_janus')
           !call setup_fluid_ladd_janus(N)
        case ('none')
           ! nop
        case ('tracer')
           call setup_fluid_tracer
        case default
           call error_md('unknown value: interaction="'//interaction//'"')
        end select
    end subroutine setup_fluid

    !> do things that have to be done after the simulation
    subroutine shutdown_fluid
        select case (interaction)
        case ('friction')
           ! nop
        case ('ladd')
           call shutdown_fluid_ladd
        case ('ladd_janus')
           call shutdown_fluid_ladd
        case ('none')
           ! nop
        case ('tracer')
           ! nop
        case default
           call error_md('unknown value: interaction="'//interaction//'"')
        end select
    end subroutine shutdown_fluid

    !> writes conclusion concerning particle-fluid interaction to all units
    !> specified in \c units .
    subroutine summary_fluid(units)
        integer,intent(in) :: units(:)

        select case (interaction)
        case ('friction')
           ! nop
        case ('ladd')
           call summary_fluid_ladd(units)
        !case ('ladd_janus')
           !call summary_fluid_ladd_janus(units)
        case ('none')
           ! nop
        case ('tracer')
           call summary_fluid_tracer(units)
        case default
           call error_md('unknown value: interaction="'//interaction//'"')
        end select
    end subroutine summary_fluid
#endif
end module lbe_md_fluid_module
