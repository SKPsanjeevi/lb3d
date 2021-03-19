#include "lbe.h"

!> adds random velocities to md particles---David Sinz knows this
module lbe_md_rand_walk_module
#ifdef MD
  use lbe_helper_module, only: local_coordinates
  use lbe_globals_module, only: halo_extent,pi
  use lbe_md_globals_module
  use lbe_md_helper_module, only: count_particles_all,sum_dx_all
  use lbe_types_module, only: lbe_site
  use luxury_rng_module, only: ranlux
  implicit none

  private
  public uni_dist,rand_walk_update_dx,rock_reflection

  !> system-wide avg dx per particle, used only during \c integrate_xq()
  real(kind=rk),save,public :: meandx

contains

  !> updates \c meandx and \c domaindx
  subroutine rand_walk_update_dx
    integer :: nparts
    real(kind=rk) :: systemdx

    call count_particles_all(nparts)
    call sum_dx_all(systemdx)
    meandx = systemdx / nparts
    if (meandx.ge.mean_free_path/delta_x) then
       domaindx = 0.0_rk
    end if
  end subroutine rand_walk_update_dx


  !> write normally distributed random number with expectancy 0 and
  !> variance 1 to  r_v . The normal distribution is obtained by the
  !> Box-Muller method.
  subroutine uni_dist(r_v)
    real*8,intent(out) :: r_v
    real*8,dimension(2) :: stor

    call RANLUX(stor,2)

    r_v = sqrt(dble(-2.0* log(1.0-stor(1)) ))*cos(2.0*pi*stor(2))
  end subroutine uni_dist

  subroutine rock_reflection(N,p_x,p_v,t1,p_v_new,p_x_new,rock)
    type(lbe_site),intent(in) :: &
         &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
    real*8,intent(in) :: p_x(3),t1,p_v(3)
    real*8 ::p_v_new(3),p_x_new(3)
    logical,intent(out) :: rock
    integer i,l
    real*8 :: r_v,k_boltz
    real*8 :: x(3),p_x_old(3),lx(3),signs(3)

    x(:)= NINT(p_x(:))

    k_boltz=1.3806504E-23

    call local_coordinates(x,lx)
    rock = .false.

    if (N(int(lx(1)),int(lx(2)),int(lx(3)))%rock_state.ne.0.0) then

      rock = .true.

      p_x_old(:)=p_x(:)-t1*p_v(:)

      x(:)= NINT(p_x_old(:))
      call local_coordinates(x,lx)
      if(N(int(lx(1)),int(lx(2)),int(lx(3)))%rock_state.ne.0.0) then

        print*,"old particle position was alreday rock probably --> stuck somehow put particle to",p_x_new(:)

      else
        do i=1,3
          p_x_new(i) = p_x_old(i)
        enddo
      endif

      do i=1,3
        call uni_dist(r_v)
        p_v_new(i) = abs(r_v)*sqrt(k_boltz*temperat/molec_mass)*(delta_t/delta_x)
      enddo

      do l=1,3
        do i=1,3
          x(i) = p_x_new(i)
        enddo
        x(l) = x(l) +t1*p_v(l)

        call local_coordinates(x,lx)
        lx(:) = NINT(lx(:))

        if(N(int(lx(1)),int(lx(2)),int(lx(3)))%rock_state.ne.0.0)then
          if(p_v(l).lt.0.0)then
            signs(l) = 1.
          else
            signs(l) = -1.
          endif
        else
          if(p_v(l).le.0.0)then
            signs(l) = -1.
          else
            signs(l) = 1.
          endif
        endif

      enddo
      do l=1,3
        p_v_new(l) = p_v_new(l)*signs(l)
      enddo

      if(.true.)then
!!!
        !Test if the new velocity causes the particle to fly in
        !rock againshouldnt happen normally but can in complex
        !geometries
!!!
        do i=1,3
          x(i)= NINT(p_x_new(i) + t1 * p_v_new(i) )
        enddo
        call local_coordinates(x,lx)
        if(N(int(lx(1)),int(lx(2)),int(lx(3)))%rock_state .ne. 0 )then
          print*,"particle will fly into rock again at",x(:),"local-coords",lx(:),"should only happen in complex geometries not at straight walls... particle position is",p_x_new(:),"new velocity is",p_v_new(:)
        endif
      endif
    endif
  end subroutine rock_reflection

#endif
end module lbe_md_rand_walk_module
