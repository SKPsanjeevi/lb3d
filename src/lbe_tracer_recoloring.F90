#include "lbe.h"

!> routines for changing the "kind" of tracers
module lbe_tracer_recoloring_module
#ifdef TRACER
    use lbe_tracer_globals_module
    use lbe_tracer_helper_module, only: passed_plane
    use lbe_tracer_helper_module, only: log_msg_tracer

    implicit none
    private

    public recolor_tracer_at_planes

    !> position of a x-plane through the system that gives passing
    !> tracers the kind "1"
    real(kind=rk),save,public :: x_recolor_plane_1=-999.0_rk
    !> position of a x-plane through the system that gives passing
    !> tracers the kind "2"
    real(kind=rk),save,public :: x_recolor_plane_2=-999.0_rk

contains

    !> recolors a tracer in case it just passed any of the recoloring
    !> planes
    !>
    !> \param[in,out] t tracer
    !>
    !> \param[in,out] px previous tracer position
    subroutine recolor_tracer_at_planes(t,px)
        type(tracer_particle_type),intent(inout) :: t
        real(kind=rk),intent(in) :: px(3)

        if (passed_plane(t%x(1),px(1),x_recolor_plane_1)) then
           t%kind = 1
        end if
        if (passed_plane(t%x(1),px(1),x_recolor_plane_2)) then
           t%kind = 2
        end if
    end subroutine recolor_tracer_at_planes

#endif
end module lbe_tracer_recoloring_module
