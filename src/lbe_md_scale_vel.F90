#include "lbe.h"

!> Scaling of velocities---Maybe David Sinz knows this?
module lbe_md_scale_vel_module
#ifdef MD
    use lbe_md_globals_module
    use lbe_md_helper_module, only: count_particles_all
    use lbe_parallel_module, only: comm_cart
    implicit none
    private
    include 'mpif.h'
    public scale_velocities
contains

    subroutine scale_velocities(temp1,molec_mass)
        integer i,ii,ierror,n_global
        real(kind=rk) factor,temp1,v_square,molec_mass,v_square_is
        real(kind=rk) vtot
        real(kind=rk),parameter :: univ_gas = 8.314472

        call count_particles_all(n_global)
!!$        if (natoms<2) call error_md('less than 2 particles'&
!!$             &//' - calculating a temperature makes absolutely no sense!')

        vtot = 0.0
        i = atompnt
        do ii = 1,nlocal
           vtot = vtot + dot_product(P(i)%v(:),P(i)%v(:))
           i = list(i)
        enddo

        !vtot is on local processor and after allreduce vtot is sum
        !over all processors of the square of the velocities
        call mpi_allreduce(vtot,vtot,3,MPI_REAL8,MPI_SUM,comm_cart,ierror)

        v_square =(temp1*univ_gas*3/(molec_mass*6.022E23)) ! v**2 SI-units

        v_square = v_square/(delta_x/delta_t)**2  ! v  L-units

        v_square_is =  vtot/n_global

        factor = sqrt(v_square/v_square_is)

        i = atompnt
        do ii = 1,nlocal
	   P(i)%v(:) = P(i)%v(:) * factor
           i = list(i)
        enddo
    end subroutine scale_velocities

#endif
end module lbe_md_scale_vel_module
