#include "lbe.h"

!> Call lubrication on neighbours
module lbe_md_potential_lubrication_mod
#ifdef MD
    use lbe_globals_module, only: tsize
    use lbe_md_fluid_ladd_module, only: particle_lubrication&
         &,set_growth_factor_particle_lubrication
    use lbe_md_globals_module

    implicit none
    private

    public force_lubrication,set_growth_factor_lubrication,setup_lubrication

contains

    !> rescales all length parameters involved in lubrication
    !> interactions
    !>
    !> \param[in] f linear length scale factor with respect to
    !> original lengths read from input file
    subroutine set_growth_factor_lubrication(f)
        real(kind=rk),intent(in) :: f

        call set_growth_factor_particle_lubrication(f)
    end subroutine set_growth_factor_lubrication

    subroutine setup_lubrication
    !NOP
    end subroutine setup_lubrication

    subroutine force_lubrication
        integer i,j,k,ii
        real(kind=rk) ptmp(3),del(3),rsq

        i = atompnt
        do ii = 1,nlocal
           ptmp(:) = P(i)%x(:)
           do k = nnlist(ii),nnlist(ii+1)-1
              j = nlist(k)
              del(:) = ptmp(:) - P(j)%x(:)
              ! ATTENTION: optimized minimum image criterion - result is
              ! useless if it is greater than cutsq1
              if (abs(del(1)).gt.tsize_mrc(1)) then
                 if (del(1).lt.0.0_rk) then
                    del(1) = del(1) + tsize(1)
                 else
                    del(1) = del(1) - tsize(1)
                 endif
              endif
              if (abs(del(2)).gt.tsize_mrc(2)) then
                 if (del(2).lt.0.0_rk) then
                    del(2) = del(2) + tsize(2)
                 else
                    del(2) = del(2) - tsize(2)
                 endif
              endif
              if (abs(del(3)).gt.tsize_mrc(3)) then
                 if (del(3).lt.0.0_rk) then
                    del(3) = del(3) + tsize(3)
                 else
                    del(3) = del(3) - tsize(3)
                 endif
              endif
              rsq = dot_product(del,del)

              if (rsq.lt.cutsq1) call particle_lubrication(rsq,del,i,j)
           enddo
           i = list(i)
        enddo
    end subroutine force_lubrication

#endif
end module lbe_md_potential_lubrication_mod
