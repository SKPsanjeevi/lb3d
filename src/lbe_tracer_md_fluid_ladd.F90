#include "lbe.h"

!> one-way coupling of \c TRACER particles to \c MD \c
!> interaction='ladd' particles
module lbe_tracer_md_fluid_ladd_module
#ifdef TRACER
#ifdef MD
    use lbe_helper_module, only: norm
    use lbe_md_boundary_condition_module, only: closest_image
    use lbe_md_fluid_ladd_module, only: in_particle
    use lbe_md_fluid_ladd_parms_module, only: R_orth,R_para
    use lbe_md_globals_module, only: p_atompnt=>atompnt,p_list=>list&
         &,p_nlocal=>nlocal,p_nother=>nother,P
    use lbe_tracer_globals_module

    implicit none
    private

    public tracer_md_fluid_ladd_interaction

contains

    !> calculates the distance vector from a point to the surface of
    !> an ellipsoid of revolution along the direction away from the
    !> ellipsoid center
    !>
    !> \param[in] o unit orientation vector parallel to the ellipsoids
    !> axis of rotational symmetry
    !>
    !> \param[in] px center position of the ellipsoid
    !>
    !> \param[in] tx point for which to calculate the distance to the surface
    !>
    !> \returns distance vector to the surface, parallel to \c tx-px
    pure function central_way_to_particle_surface(o,px,tx)
        real(kind=rk) :: central_way_to_particle_surface(3)
        real(kind=rk),intent(in) :: o(3),px(3),tx(3)
        real(kind=rk) :: rx(3),rxo(3),rxp(3),rx_scaled(3)

        rx = tx-px
        rxp = o*dot_product(rx,o)
        rxo = rx-rxp
        rx_scaled = rxp/R_para + rxo/R_orth
        central_way_to_particle_surface = rx*(1.0_rk/norm(rx_scaled)-1.0_rk)
    end function central_way_to_particle_surface

    !> interaction of a tracer and \c MD \c interaction=='ladd'
    !> particles with the purpose to avoid tracers entering the
    !> theoretical volume occupied by a finite size particle
    !>
    !> \param[in,out] t tracer particle
    subroutine tracer_md_fluid_ladd_interaction(t)
        type(tracer_particle_type),intent(inout) :: t
        integer :: j,jj,n_collisions
        real(kind=rk) :: offset(3),px(3)

        offset = 0.0_rk
        n_collisions = 0

        ! check for overlap of this tracer with any of the finite
        ! size particles
        j = p_atompnt
        particles: do jj=1,p_nlocal+p_nother
           ! nearest image of Ladd particle center
           px = closest_image(P(j)%x,t%x)

           within_particle: if (in_particle(P(j),px,t%x)) then
              n_collisions = n_collisions + 1

              ! suggest a way for the tracer to leave this finite
              ! size particle volume
              offset = offset&
                   &+central_way_to_particle_surface(P(j)%o,px,t%x)
           end if within_particle

           if (jj<=p_nlocal) then
              j = p_list(j)
           else
              j = j + 1
           endif
        end do particles

        if (n_collisions>0) then
           ! in case of multiple collisions, average the resulting
           ! offsets. This makes the result independent from the
           ! order of particles in P
           if (n_collisions>1) offset = offset/real(n_collisions,kind=rk)
           t%x = t%x + offset
        end if
    end subroutine tracer_md_fluid_ladd_interaction

#endif
#endif
end module lbe_tracer_md_fluid_ladd_module
