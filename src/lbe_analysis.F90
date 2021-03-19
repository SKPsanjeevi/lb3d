#include "lbe.h"

!> routines to calculate physical properties
module lbe_analysis_module
    use lbe_globals_module
    use lbe_helper_module, only: density,is_fluid,massflow
    use lbe_parms_module, only: amass_b,amass_r,amass_s,nt,nx,ny,nz
    use lbe_types_module, only: lbe_site
#ifdef MD
    use lbe_md_analysis_module, only: md_total_mass,md_total_momentum
#endif

    implicit none
    private

    public avg_total_velocity,cached_avg_total_velocity,total_mass

    include 'mpif.h'

contains

    !> finds the average velocity of the total system
    !>
    !> \param[in] whole_N local chunk of the lattice with full halo of
    !> depth \c halo_extent
    !>
    !> \returns average velocity vector
    !>
    !> If compiled with \c MD, also the particle velocities are taken
    !> into account. In any case, weighting is done by the mass
    !> density.
    function avg_total_velocity(whole_N)
        real(kind=rk),dimension(3) :: avg_total_velocity
        type(lbe_site),intent(in) :: &
             &whole_N(1-halo_extent:,1-halo_extent:,1-halo_extent:)

        avg_total_velocity = total_momentum(whole_N)/total_mass(whole_N)
    end function avg_total_velocity

    !> finds the average velocity of the total system
    !>
    !> \param[in] whole_N local chunk of the lattice with full halo of
    !> depth \c halo_extent
    !>
    !> \returns average velocity vector
    !>
    !> If compiled with \c MD, also the particle velocities are taken
    !> into account. In any case, weighting is done by the mass
    !> density. The result is calculated only once per time step.
    function cached_avg_total_velocity(whole_N)
        real(kind=rk),dimension(3) :: cached_avg_total_velocity
        type(lbe_site),intent(in) :: &
             &whole_N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        integer,save :: last_calculation=-1
        real(kind=rk),save :: last_result(3)

        if (last_calculation/=nt) then
           last_result = avg_total_velocity(whole_N)
           last_calculation = nt
        end if
        cached_avg_total_velocity = last_result
    end function cached_avg_total_velocity

    !> finds the total fluid mass in the system
    !>
    !> \param[in] whole_N local chunk of the lattice with full halo of
    !> depth \c halo_extent
    !>
    !> \returns total mass of each species; elements 1, 2, and 3 are
    !> for species red, blue, and green, respectively, if existing
    function total_fluid_mass(whole_N)
        real(kind=rk),dimension(n_spec) :: total_fluid_mass
        type(lbe_site),intent(in) :: &
             &whole_N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        real(kind=rk) :: lsum(n_spec),ssum(n_spec)
        integer ierror,x,y,z

        lsum = 0.0_rk
        do x=1,nx
           do y=1,ny
              do z=1,nz
                 if (is_fluid(whole_N(x,y,z)%rock_state)) then
                    lsum(1) = lsum(1) + amass_r*density(whole_N(x,y,z)%n_r)
#ifndef SINGLEFLUID
                    lsum(2) = lsum(2) + amass_b*density(whole_N(x,y,z)%n_b)
#ifndef NOSURFACTANT
                    lsum(3) = lsum(3) + amass_s*density(whole_N(x,y,z)%n_s)
#endif
#endif
                 end if
              end do
           end do
        end do

        call MPI_Allreduce(lsum,ssum,n_spec,LBE_REAL,MPI_SUM,MPI_COMM_WORLD&
             &,ierror)
        total_fluid_mass = ssum
    end function total_fluid_mass

    !> finds the total fluid momentum in the system
    !>
    !> \param[in] whole_N local chunk of the lattice with full halo of
    !> depth \c halo_extent
    !>
    !> \returns total momentum vector
    function total_fluid_momentum(whole_N)
        real(kind=rk),dimension(3) :: total_fluid_momentum
        type(lbe_site),intent(in) :: &
             &whole_N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        real(kind=rk) :: lsum(3),ssum(3)
        integer ierror,x,y,z

        lsum = 0.0_rk
        do x=1,nx
           do y=1,ny
              do z=1,nz
                 if (is_fluid(whole_N(x,y,z)%rock_state)) &
                      &lsum = lsum + massflow(whole_N(x,y,z))
              end do
           end do
        end do

        call MPI_Allreduce(lsum,ssum,3,LBE_REAL,MPI_SUM,MPI_COMM_WORLD,ierror)
        total_fluid_momentum = ssum
    end function total_fluid_momentum

    !> finds the total mass in the system
    !>
    !> \param[in] whole_N local chunk of the lattice with full halo of
    !> depth \c halo_extent
    !>
    !> \returns total mass
    !>
    !> If compiled with \c MD, also the particle mass is taken
    !> into account.
    function total_mass(whole_N)
        real(kind=rk) :: total_mass
        type(lbe_site),intent(in) :: &
             &whole_N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        real(kind=rk) :: ret

        ret = sum(total_fluid_mass(whole_N))
#ifdef MD
        ret = ret + md_total_mass()
#endif
        total_mass = ret
    end function total_mass

    !> finds the total momentum in the system
    !>
    !> \param[in] whole_N local chunk of the lattice with full halo of
    !> depth \c halo_extent
    !>
    !> \returns total momentum vector
    !>
    !> If compiled with \c MD, also the particle momentum is taken
    !> into account.
    function total_momentum(whole_N)
        real(kind=rk),dimension(3) :: total_momentum
        type(lbe_site),intent(in) :: &
             &whole_N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        real(kind=rk) :: ret(3)

        ret = total_fluid_momentum(whole_N)
#ifdef MD
        ret = ret + md_total_momentum()
#endif
        total_momentum = ret
    end function total_momentum

end module lbe_analysis_module
