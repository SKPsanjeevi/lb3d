#include "lbe.h"

!> routines to calculate physical properties
module lbe_md_analysis_module
#ifdef MD
    use lbe_globals_module, only: rk
    use lbe_md_globals_module
    use lbe_md_helper_module, only: particle_mass

    implicit none
    private

    public md_total_mass,md_total_momentum

    include 'mpif.h'

contains

    !> finds the total mass of all MD particles
    !>
    !> \returns total mass
    function md_total_mass()
        real(kind=rk) :: md_total_mass
        integer :: i,ierror,ii
        real(kind=rk) :: lsum,ssum

        lsum = 0.0_rk
        i = atompnt
        do ii=1,nlocal
           lsum = lsum + particle_mass(P(i))

           i = list(i)
        end do

        call MPI_Allreduce(lsum,ssum,1,LBE_REAL,MPI_SUM,MPI_COMM_WORLD,ierror)
        md_total_mass = ssum
    end function md_total_mass

    !> finds the total momentum of all MD particles
    !>
    !> \returns total momentum vector
    function md_total_momentum()
        real(kind=rk),dimension(3) :: md_total_momentum
        integer :: i,ierror,ii
        real(kind=rk) :: lsum(3),ssum(3)

        lsum = 0.0_rk
        i = atompnt
        do ii=1,nlocal
           if (use_ft_fluid) then
              lsum = lsum + particle_mass(P(i)) * P(i)%v_fluid_avg
           else
              lsum = lsum + particle_mass(P(i)) * P(i)%v
           end if

           i = list(i)
        end do

        call MPI_Allreduce(lsum,ssum,3,LBE_REAL,MPI_SUM,MPI_COMM_WORLD,ierror)
        md_total_momentum = ssum
    end function md_total_momentum

#endif
end module lbe_md_analysis_module
