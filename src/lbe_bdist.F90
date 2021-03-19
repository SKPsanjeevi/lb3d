#include "lbe.h"

!>This contains the function to calculate the equilibrium distribution.
!>
!>The bdist parameter in the input file chooses the expression that is
!>used.
module lbe_bdist_module
  use lbe_globals_module
  use lbe_parallel_module
  use lbe_parms_module, only: bdist,kbT
  use luxury_rng_module, only : ranlux


  implicit none
  private

  public boltz_dist, boltz_dist_s, mrt_dist, mrt_init_dist, equilibrium_distribution, boltz_dist_constrho

contains
  subroutine boltz_dist_constrho(ux, uy, uz, dux, duy, duz, ddx, ddy, ddz, rho_eq, boltz)
    real(kind=rk) :: ux, uy, uz
    real(kind=rk) :: dux, duy, duz, ddx, ddy, ddz, rho
    real(kind=rk),intent(out) :: boltz(nvecs)

    real(kind=rk),parameter :: aN1=1.0_rk/36.0_rk,aN0=1.0_rk/3.0_rk
    real(kind=rk),parameter :: haN1 = 0.5_rk*aN1
    real(kind=rk),parameter :: aN1_6 = aN1/6.0_rk
    real(kind=rk),parameter :: holdit = 1.0_rk/T

    integer :: s, ifactor, nni

    real(kind=rk) :: const,const1,even,odd,tux,tuy,tuz,tu1x,tu1y,tu1z
    real(kind=rk) :: tuxt,tuyt,tuzt,uke,u1ke,dke
    real(kind=rk) :: cdotu,cdot1u,cdotd,cdotdu_s,udotdu_s
    real(kind=rk) :: ux_f,uy_f,uz_f,uke_eq,u_mag,delta
    real(kind=rk) :: expux,expuy,expuz,expek,a1eq
    real(kind=rk) :: n_iter
    real(kind=rk) :: rho1_tt,rho_tt,rho_eq
    real(kind=rk) :: px_tt,py_tt,pz_tt
    real(kind=rk) :: delt_ux,delt_uy,delt_uz
    real(kind=rk) :: ux_eq,uy_eq,uz_eq
    real(kind=rk) :: uholdit, choldit, cdholdit, ch2

    ux = ux+dux
    uy = uy+duy
    uz = uz+duz

!   write(*,"('const density scheme')") 
       uke = 0.5_rk*(ux*ux + uy*uy + uz*uz)
       uholdit = uke*holdit

       do s = 1, nnonrest
          cdotu = cx(s)*ux + cy(s)*uy + cz(s)*uz
          choldit = cdotu*holdit
          ch2 = choldit*choldit
!          boltz(s) = aN1*(rho_eq + 1.0_rk*(cdotu/T &
!               + 0.5_rk*((cdotu/T)**2) &
!               - uke/T ))
          boltz(s) = aN1*rho_eq*( 1.0_rk + cdotu/T &
               + 0.5_rk*((cdotu/T)**2) &
               - uke/T )

       end do
!       boltz(restvec) = aN0*(rho_eq - 1.0_rk*(uke)/T)
       boltz(restvec) = aN0*rho_eq*(1.0_rk - (uke)/T)
       return
end subroutine boltz_dist_constrho


    !> Equilibrium distribution function
    !>
    !> Takes a velocity specified by \c u, an increment \c du, and
    !> and array \c boltz of 19 elements; returns with \c boltz filled in with
    !> the equilibrium distribution of velocities.
    !>
    !> Uses either the Chen operator or the Luo operator, depending on the
    !> value of the input variable \c bdist.
  subroutine boltz_dist(ux, uy, uz, dux, duy, duz, ddx, ddy, ddz, boltz)
    real(kind=rk) :: ux, uy, uz
    real(kind=rk) :: dux, duy, duz, ddx, ddy, ddz, rho
    real(kind=rk),intent(out) :: boltz(nvecs)

    real(kind=rk),parameter :: aN1=1.0_rk/36.0_rk,aN0=1.0_rk/3.0_rk
    real(kind=rk),parameter :: haN1 = 0.5_rk*aN1
    real(kind=rk),parameter :: aN1_6 = aN1/6.0_rk
    real(kind=rk),parameter :: holdit = 1.0_rk/T

    integer :: s, ifactor, nni

    real(kind=rk) :: const,const1,even,odd,tux,tuy,tuz,tu1x,tu1y,tu1z
    real(kind=rk) :: tuxt,tuyt,tuzt,uke,u1ke,dke
    real(kind=rk) :: cdotu,cdot1u,cdotd,cdotdu_s,udotdu_s
    real(kind=rk) :: ux_f,uy_f,uz_f,uke_eq,u_mag,delta
    real(kind=rk) :: expux,expuy,expuz,expek,a1eq
    real(kind=rk) :: n_iter
    real(kind=rk) :: rho1_tt,rho_tt,rho_eq
    real(kind=rk) :: px_tt,py_tt,pz_tt
    real(kind=rk) :: delt_ux,delt_uy,delt_uz
    real(kind=rk) :: ux_eq,uy_eq,uz_eq
    real(kind=rk) :: uholdit, choldit, cdholdit, ch2

#ifndef FASTBDIST2
    select case (bdist)

       ! Luo operator
    case (0)
       uke = 0.5_rk*(ux*ux + uy*uy + uz*uz)
       udotdu_s = ux*dux + uy*duy + uz*duz
       uholdit = uke*holdit

       do s = 1, nnonrest
          cdotu = cx(s)*ux + cy(s)*uy + cz(s)*uz
          cdotdu_s = cx(s)*dux + cy(s)*duy + cz(s)*duz
          choldit = cdotu*holdit
          cdholdit = cdotdu_s*holdit
          ch2 = choldit*choldit

          !Elena: Replace division by T
          boltz(s) = aN1*( 1.0_rk + choldit + 0.5_rk*ch2 - uholdit &
               + (choldit*ch2)/6.0_rk - uholdit*choldit &
               + cdholdit + choldit*cdholdit - udotdu_s*holdit &
               + 0.5_rk - choldit*udotdu_s*holdit &
               - uholdit*cdholdit )
       end do
       boltz(restvec) = aN0*(1.0_rk - uke/T - udotdu_s/T)
       return

       ! Hudong operator.
    case(1)
       ! This is a third order expansion of f(u+delta_u) with no terms missing
       ! as described in Proc. R. Soc. Lond. A (2000) 456, 2043-2057
       ! PRL 81 1618 (1998) and PRE 58 6855 (1998) both claim that this is
       ! inconsistent with a systematic derivation of how force terms should be
       ! introduced. However, it is desirably more stable than the other
       ! offerings here.
       ! Phys. Rev. E 75 036712 use it as well...
       ux_f = ux + dux
       uy_f = uy + duy
       uz_f = uz + duz

       uke_eq = 0.5_rk*(ux_f*ux_f + uy_f*uy_f + uz_f*uz_f)

       u_mag = sqrt(2.0_rk*uke_eq)
       ifactor = 200.0_rk*u_mag

       n_iter = max(2,ifactor)
       delta = sqrt(1.0_rk/n_iter)

       ux_eq = ux_f
       uy_eq = uy_f
       uz_eq = uz_f

       rho_eq = 1.0_rk

       do nni = 1, int(n_iter)
          expux = exp(ux_eq/T)
          expuy = exp(uy_eq/T)
          expuz = exp(uz_eq/T)
          expek = 1.0_rk/exp(uke_eq/T)

          do s = 1, nnonrest
             a1eq = (expux**cx(s))*(expuy**cy(s))*(expuz**cz(s))
             boltz(s) = rho_eq*aN1*a1eq*expek
          end do
          boltz(restvec) = rho_eq*aN0*expek

          rho1_tt = 0.0_rk
          px_tt = 0.0_rk
          py_tt = 0.0_rk
          pz_tt = 0.0_rk

          do s = 1, nnonrest
             rho1_tt = rho1_tt + boltz(s)*g(s)
             px_tt = px_tt + boltz(s)*g(s)*cx(s)
             py_tt = py_tt + boltz(s)*g(s)*cy(s)
             pz_tt = pz_tt + boltz(s)*g(s)*cz(s)
          end do

          rho_tt = rho1_tt + boltz(nvecs)

          delt_ux = ux_f - px_tt
          delt_uy = uy_f - py_tt
          delt_uz = uz_f - pz_tt

          rho_eq = rho_eq/(1.0_rk + delta*(rho_tt - 1.0_rk )/rho_eq)
          ux_eq = ux_eq + delta*delt_ux
          uy_eq = uy_eq + delta*delt_uy
          uz_eq = uz_eq + delta*delt_uz
          uke_eq = 0.5_rk*(ux_eq*ux_eq + uy_eq*uy_eq + uz_eq*uz_eq)
       end do
       return

    case (2)
#endif
       tux = ux+dux
       tuy = uy+duy
       tuz = uz+duz

       tuxt = tux*holdit
       tuyt = tuy*holdit
       tuzt = tuz*holdit

       uke = 0.5_rk*(tux*tuxt + tuy*tuyt + tuz*tuzt)

       const = 1.0_rk-uke
       const1 = aN1*const

       cdotu = tuxt
       even = const1+haN1*cdotu**2
       odd = const1*cdotu+aN1_6*cdotu**3
       boltz(1) = even+odd
       boltz(2) = even-odd
       cdotu = tuyt
       even = const1+haN1*cdotu**2
       odd = const1*cdotu+aN1_6*cdotu**3
       boltz(3) = even+odd
       boltz(4) = even-odd
       cdotu = tuzt
       even = const1+haN1*cdotu**2
       odd = const1*cdotu+aN1_6*cdotu**3
       boltz(5) = even+odd
       boltz(6) = even-odd
       cdotu = tuxt+tuyt
       even = const1+haN1*cdotu**2
       odd = const1*cdotu+aN1_6*cdotu**3
       boltz(7) = even+odd
       boltz(12) = even-odd
       cdotu = tuxt-tuyt
       even = const1+haN1*cdotu**2
       odd = const1*cdotu+aN1_6*cdotu**3
       boltz(8) = even+odd
       boltz(11) = even-odd
       cdotu = tuxt+tuzt
       even = const1+haN1*cdotu**2
       odd = const1*cdotu+aN1_6*cdotu**3
       boltz(9) = even+odd
       boltz(14) = even-odd
       cdotu = tuxt-tuzt
       even = const1+haN1*cdotu**2
       odd = const1*cdotu+aN1_6*cdotu**3
       boltz(10) = even+odd
       boltz(13) = even-odd
       cdotu = tuyt+tuzt
       even = const1+haN1*cdotu**2
       odd = const1*cdotu+aN1_6*cdotu**3
       boltz(15) = even+odd
       boltz(18) = even-odd
       cdotu = tuyt-tuzt
       even = const1+haN1*cdotu**2
       odd = const1*cdotu+aN1_6*cdotu**3
       boltz(16) = even+odd
       boltz(17) = even-odd

       boltz(restvec) = aN0*const

#ifndef FASTBDIST2
    case (3)
       ! This is a second order expansion of f(u+delta_u) with no terms missing
       ! PRL 81 1618 (1998) and PRE 58 6855 (1998) both claim that this is
       ! inconsistent with a systematic derivation of how force terms should be
       ! introduced.
       uke = 0.5_rk*(ux*ux + uy*uy + uz*uz)
       udotdu_s = ux*dux + uy*duy + uz*duz

       do s = 1, nnonrest
          cdotu = cx(s)*ux + cy(s)*uy + cz(s)*uz
          cdotdu_s = cx(s)*dux + cy(s)*duy + cz(s)*duz

          boltz(s) = aN1*(1.0_rk + cdotu/T &
               + 0.5_rk*(cdotu/T)**2 - uke/T &
               + cdotdu_s/T + cdotu*cdotdu_s/T**2 - udotdu_s/T)
       end do
       boltz(restvec) = aN0*(1.0_rk - uke/T - udotdu_s/T)
       return


    case (4)
       ! Second order expansion for the equilibrium distribution
       ! and force implemmentation according to
       ! Phys. Rev. E 65,(2002): Guo

       ! u = u_tilde
       tux = ux + dux
       tuy = uy + duy
       tuz = uz + duz

       uke = 0.5_rk*(tux*tux + tuy*tuy + tuz*tuz)

       ! dd#  = (tau - 0.5)*g_accn_#
       dke = 0.5_rk*(ddx*ddx + ddy*ddy + ddz*ddz)

       do s = 1, nnonrest
          cdotu = cx(s)*tux + cy(s)*tuy + cz(s)*tuz
          cdotd = cx(s)*ddx + cy(s)*ddy + cz(s)*ddz

          boltz(s) = aN1*( 1.0_rk + cdotu/T &
               + 0.5_rk*((cdotu/T)**2 - (cdotd/T)**2) &
               - (uke-dke)/T )
       end do
       boltz(restvec) = aN0*(1.0_rk - (uke-dke)/T)
       return

    case (5)
       ! Third order expansion for the equilibrium distribution
       ! and force implemmentation according to
       ! Phys. Rev. E 65,(2002): Guo

       ! ux = u_tilde
       tux = ux + dux
       tuy = uy + duy
       tuz = uz + duz

       uke = 0.5_rk*(tux*tux + tuy*tuy + tuz*tuz)

       ! dd#  = (tau - 0.5)*g_accn
       dke = 0.5_rk*(ddx*ddx + ddy*ddy + ddz*ddz)

       ! u_tile + f_{Shan Chen}_r + 0.5*g_accn
       tu1x = ux + dux - ddx
       tu1y = uy + duy - ddy
       tu1z = uz + duz - ddz

       u1ke = 0.5_rk*(tu1x*tu1x + tu1y*tu1y + tu1z*tu1z)

       do s = 1, nnonrest
          cdotu = cx(s)*tux + cy(s)*tuy + cz(s)*tuz
          cdot1u = cx(s)*tu1x + cy(s)*tu1y + cz(s)*tu1z
          cdotd = cx(s)*ddx + cy(s)*ddy + cz(s)*ddz

          boltz(s) = aN1*( 1.0_rk + cdotu/T &
               + 0.5_rk*((cdotu/T)**2 - (cdotd/T)**2) &
               - (uke-dke)/T &
               - uke*cdot1u/T**2 + (cdot1u/T)**3/6.0_rk &
               + 0.5_rk*(u1ke/T)**2 - 0.5*(cdot1u**2*u1ke)/T**3 &
               + (cdot1u/T)**4/24.0_rk )
       end do
       boltz(restvec) = aN0*(1.0_rk - (uke-dke)/T + 0.5_rk*(uke/T)**2)
       return

       ! Corrected Luo 'operator'
       ! Extended to fourth order (as a test)
       ! Particle number and momentum are conserved - however
       ! I believe that other 'hydronamic moments' are not.
       ! Do not use/trust this code!
       ! However, it is interesting that this code agrees
       ! very closely with the Hudong code (bdist=0)
    case (6)
       tux = ux + dux
       tuy = uy + duy
       tuz = uz + duz

       uke = 0.5_rk*(tux*tux + tuy*tuy + tuz*tuz)

       do s = 1, nnonrest
          cdotu = cx(s)*tux + cy(s)*tuy + cz(s)*tuz

          boltz(s) = aN1*(1.0_rk + cdotu/T &
               + 0.5_rk*(cdotu/T)**2 - uke/T &
               - uke*cdotu/T**2 + (cdotu/T)**3/6.0_rk &
               + 0.5_rk*(uke/T)**2 - 0.5*(cdotu**2*uke)/T**3 &
               + (cdotu/T)**4/24.0_rk)
       end do
       boltz(restvec) = aN0*(1.0_rk - uke/T + 0.5_rk*(uke/T)**2)
       return

       ! This is a second order expansion of f(u+delta_u) with no
       ! terms missing PRL 81 1618 (1998) and PRE 58 6855 (1998) both
       ! claim that this is inconsistent with a systematic derivation
       ! of how force terms should be introduced. So this is the same
       ! as bdist==4, but optimized by manually unrolling the loop
       ! and exploiting zero-valued components of the lattice
       ! vectors.
    case (7)
       ! tux = ux+dux
       ! tuy = uy+duy
       ! tuz = uz+duz

       ! tuxt = tux*holdit
       ! tuyt = tuy*holdit
       ! tuzt = tuz*holdit

       ! tuxta = tuxt*aN1
       ! tuyta = tuyt*aN1
       ! tuzta = tuzt*aN1

       ! uke = 0.5d0*(tux*tuxt + tuy*tuyt + tuz*tuzt)

       ! const = 1.0_rk-uke
       ! const1 = aN1*const

       ! cdotu = tuxt
       ! even = const1+haN1*cdotu**2
       ! odd = tuxta
       ! boltz(1) = even+odd
       ! boltz(2) = even-odd

       ! cdotu = tuyt
       ! even = const1+haN1*cdotu**2
       ! odd = tuyta
       ! boltz(3) = even+odd
       ! boltz(4) = even-odd

       ! cdotu = tuzt
       ! even = const1+haN1*cdotu**2
       ! odd = tuzta
       ! boltz(5) = even+odd
       ! boltz(6) = even-odd

       ! cdotu = tuxt+tuyt
       ! even = const1+haN1*cdotu**2
       ! odd = tuxta+tuyta
       ! boltz(7) = even+odd
       ! boltz(12) = even-odd

       ! cdotu = tuxt-tuyt
       ! even = const1+haN1*cdotu**2
       ! odd = tuxta-tuyta
       ! boltz(8) = even+odd
       ! boltz(11) = even-odd

       ! cdotu = tuxt+tuzt
       ! even = const1+haN1*cdotu**2
       ! odd = tuxta+tuzta
       ! boltz(9) = even+odd
       ! boltz(14) = even-odd

       ! cdotu = tuxt-tuzt
       ! even = const1+haN1*cdotu**2
       ! odd = tuxta-tuzta
       ! boltz(10) = even+odd
       ! boltz(13) = even-odd

       ! cdotu = tuyt+tuzt
       ! even = const1+haN1*cdotu**2
       ! odd = tuyta-tuzta
       ! boltz(15) = even+odd
       ! boltz(18) = even-odd

       ! cdotu = tuyt-tuzt
       ! even = const1+haN1*cdotu**2
       ! odd = tuyta-tuzta
       ! boltz(16) = even+odd
       ! boltz(17) = even-odd

       ! boltz(restvec) = aN0*const
       return
    end select
#endif
end subroutine boltz_dist

!> returns only component \c s of \c boltz_dist() as defined above
subroutine boltz_dist_s(ux,uy,uz,dux,duy,duz,ddx,ddy,ddz,boltz_s,s)
    real(kind=rk),intent(in) :: ux,uy,uz,dux,duy,duz, ddx, ddy, ddz
    real(kind=rk),intent(out) :: boltz_s
    integer,intent(in) :: s
#ifdef FASTBDIST2
    real(kind=rk),parameter :: aN1=1.0_rk/36.0_rk,aN0=1.0_rk/3.0_rk
    real(kind=rk),parameter :: haN1 = 0.5_rk*aN1
    real(kind=rk),parameter :: aN1_6 = aN1/6.0_rk
    real(kind=rk),parameter :: holdit = 1.0_rk/T
    real(kind=rk) :: cdotu,const,tux,tuy,tuz
#else
    real(kind=rk) :: boltz(nvecs)
#endif

#ifdef FASTBDIST2
    tux = ux+dux
    tuy = uy+duy
    tuz = uz+duz
    const = 1.0_rk-0.5_rk*(tux*tux + tuy*tuy + tuz*tuz)*holdit

    select case (s)
       case (1)
           cdotu = tux*holdit
           boltz_s = aN1*const*(1.0_rk+cdotu)+haN1*cdotu**2+aN1_6*cdotu**3
       case (2)
           cdotu = tux*holdit
           boltz_s = aN1*const*(1.0_rk-cdotu)+haN1*cdotu**2-aN1_6*cdotu**3
       case (3)
           cdotu = tuy*holdit
           boltz_s = aN1*const*(1.0_rk+cdotu)+haN1*cdotu**2+aN1_6*cdotu**3
       case (4)
           cdotu = tuy*holdit
           boltz_s = aN1*const*(1.0_rk-cdotu)+haN1*cdotu**2-aN1_6*cdotu**3
       case (5)
           cdotu = tuz*holdit
           boltz_s = aN1*const*(1.0_rk+cdotu)+haN1*cdotu**2+aN1_6*cdotu**3
       case (6)
           cdotu = tuz*holdit
           boltz_s = aN1*const*(1.0_rk-cdotu)+haN1*cdotu**2-aN1_6*cdotu**3
       case (7)
           cdotu = (tux+tuy)*holdit
           boltz_s = aN1*const*(1.0_rk+cdotu)+haN1*cdotu**2+aN1_6*cdotu**3
       case (8)
           cdotu = (tux-tuy)*holdit
           boltz_s = aN1*const*(1.0_rk+cdotu)+haN1*cdotu**2+aN1_6*cdotu**3
       case (9)
           cdotu = (tux+tuz)*holdit
           boltz_s = aN1*const*(1.0_rk+cdotu)+haN1*cdotu**2+aN1_6*cdotu**3
       case (10)
           cdotu = (tux-tuz)*holdit
           boltz_s = aN1*const*(1.0_rk+cdotu)+haN1*cdotu**2+aN1_6*cdotu**3
       case (11)
           cdotu = (tux-tuy)*holdit
           boltz_s = aN1*const*(1.0_rk-cdotu)+haN1*cdotu**2-aN1_6*cdotu**3
       case (12)
           cdotu = (tux+tuy)*holdit
           boltz_s = aN1*const*(1.0_rk-cdotu)+haN1*cdotu**2-aN1_6*cdotu**3
       case (13)
           cdotu = (tux-tuz)*holdit
           boltz_s = aN1*const*(1.0_rk-cdotu)+haN1*cdotu**2-aN1_6*cdotu**3
       case (14)
           cdotu = (tux+tuz)*holdit
           boltz_s = aN1*const*(1.0_rk-cdotu)+haN1*cdotu**2-aN1_6*cdotu**3
       case (15)
           cdotu = (tuy+tuz)*holdit
           boltz_s = aN1*const*(1.0_rk+cdotu)+haN1*cdotu**2+aN1_6*cdotu**3
       case (16)
           cdotu = (tuy-tuz)*holdit
           boltz_s = aN1*const*(1.0_rk+cdotu)+haN1*cdotu**2+aN1_6*cdotu**3
       case (17)
           cdotu = (tuy-tuz)*holdit
           boltz_s = aN1*const*(1.0_rk-cdotu)+haN1*cdotu**2-aN1_6*cdotu**3
       case (18)
           cdotu = (tuy+tuz)*holdit
           boltz_s = aN1*const*(1.0_rk-cdotu)+haN1*cdotu**2-aN1_6*cdotu**3
       case (19)
           boltz_s = aN0*const
    end select
#else
    call boltz_dist(ux,uy,uz,dux,duy,duz,ddx,ddy,ddz,boltz)
    boltz_s = boltz(s)
#endif
end subroutine boltz_dist_s

! !> returns only 0th and 1st order terms of \c boltz_dist_s() with
! !> \c bdist=2 (3rd order bdist).
! !>
! !> \param[in] tux x velocity component
! !>
! !> \param[in] tux y velocity component
! !>
! !> \param[in] tux z velocity component
! !>
! !> \param[out] boltz_s component \c s of the distribution function
! !>
! !> \param[in] s specifying which component of the distribution
! !> function to return in \c boltz_s
! !>
! !> If necessary, \c (/dux,duy,duz/) should be added to \c (/ux,uy,uz/)
! !> and then passed to this subroutine as sum.
! subroutine boltz_dist1_s(tux,tuy,tuz,boltz_s,s)
!     real(kind=rk),intent(in) :: tux,tuy,tuz
!     real(kind=rk),intent(out) :: boltz_s
!     integer,intent(in) :: s
!     real(kind=rk),parameter :: aN1=1.0_rk/36.0_rk,aN0=1.0_rk/3.0_rk
!     real(kind=rk),parameter :: holdit = 1.0_rk/T

!     select case (s)
!        case (1)
!            boltz_s = aN1*(1.0_rk+tux*holdit)
!        case (2)
!            boltz_s = aN1*(1.0_rk-tux*holdit)
!        case (3)
!            boltz_s = aN1*(1.0_rk+tuy*holdit)
!        case (4)
!            boltz_s = aN1*(1.0_rk-tuy*holdit)
!        case (5)
!            boltz_s = aN1*(1.0_rk+tuz*holdit)
!        case (6)
!            boltz_s = aN1*(1.0_rk-tuz*holdit)
!        case (7)
!            boltz_s = aN1*(1.0_rk+(tux+tuy)*holdit)
!        case (8)
!            boltz_s = aN1*(1.0_rk+(tux-tuy)*holdit)
!        case (9)
!            boltz_s = aN1*(1.0_rk+(tux+tuz)*holdit)
!        case (10)
!            boltz_s = aN1*(1.0_rk+(tux-tuz)*holdit)
!        case (11)
!            boltz_s = aN1*(1.0_rk-(tux-tuy)*holdit)
!        case (12)
!            boltz_s = aN1*(1.0_rk-(tux+tuy)*holdit)
!        case (13)
!            boltz_s = aN1*(1.0_rk-(tux-tuz)*holdit)
!        case (14)
!            boltz_s = aN1*(1.0_rk-(tux+tuz)*holdit)
!        case (15)
!            boltz_s = aN1*(1.0_rk+(tuy+tuz)*holdit)
!        case (16)
!            boltz_s = aN1*(1.0_rk+(tuy-tuz)*holdit)
!        case (17)
!            boltz_s = aN1*(1.0_rk-(tuy-tuz)*holdit)
!        case (18)
!            boltz_s = aN1*(1.0_rk-(tuy+tuz)*holdit)
!        case (19)
!            boltz_s = aN0
!     end select
! end subroutine boltz_dist1_s

  subroutine mrt_dist(ueq, aaccn, locN, S)
    ! ueq = sum(ni*ci)/sum(ni) + f_{Shan Chen} + g_accn/2
    ! aaccn = g_accn

    !density
    real(kind=rk) :: rho, rho2
    !velocity
    real(kind=rk), dimension(3) :: u, ueq, aaccn
    real(kind=rk) :: jx, jy, jz
    real(kind=rk) :: jx2, jy2, jz2, jj

    !local occupation densities
    real(kind=rk), dimension(19) :: locN

    !collision "matrix" is diag so need a vec only
    real(kind=rk), dimension(19) :: S

    !moments
    real(kind=rk), dimension(19) :: m, meq, mnew, meqf

    !helper moments
    real(kind=rk), dimension(9)  :: mp, mn
    real(kind=rk), dimension(19) :: minv
    real(kind=rk) :: m0

    real(kind=rk),parameter :: holdit = 1.0_rk/T

    real(kind=rk) :: fluctuating_amplitude,gaussian1,gaussian2

    !local density
    rho = sum( locN(:)*g(:) )
    rho = max(rho,1.d-9)
    rho2 = rho**2

    !local velocities
    u(1) = sum( locN(:)*cx(:)*g(:) )
    u(2) = sum( locN(:)*cy(:)*g(:) )
    u(3) = sum( locN(:)*cz(:)*g(:) )

    !construct equilibrium density flux
    jx = ueq(1)*rho
    jy = ueq(2)*rho
    jz = ueq(3)*rho

    jx2 = jx**2
    jy2 = jy**2
    jz2 = jz**2
    jj  = jx2 + jy2 + jz2

    !calculate helper moments for transformation to moment space
    mp(1) = locN(1)  * g(1)  + locN(2)  * g(2) !x
    mp(2) = locN(3)  * g(3)  + locN(4)  * g(4) !y
    mp(3) = locN(5)  * g(5)  + locN(6)  * g(6) !z
    mp(4) = locN(7)  * g(7)  + locN(8)  * g(8)  !x xy
    mp(5) = locN(9)  * g(9)  + locN(10) * g(10) ! x xz
    mp(6) = locN(11) * g(11) + locN(12) * g(12) !-x xy
    mp(7) = locN(13) * g(13) + locN(14) * g(14) !-x xz
    mp(8) = locN(15) * g(15) + locN(16) * g(16) ! y yz
    mp(9) = locN(17) * g(17) + locN(18) * g(18) !-y yz

    mn(1) = locN(1)  * g(1)  - locN(2)  * g(2) !x
    mn(2) = locN(3)  * g(3)  - locN(4)  * g(4) !y
    mn(3) = locN(5)  * g(5)  - locN(6)  * g(6) !z
    mn(4) = locN(7)  * g(7)  - locN(8)  * g(8)  ! x xy
    mn(5) = locN(9)  * g(9)  - locN(10) * g(10) ! x xz
    mn(6) = locN(11) * g(11) - locN(12) * g(12) !-x xy
    mn(7) = locN(13) * g(13) - locN(14) * g(14) !-x xz
    mn(8) = locN(15) * g(15) - locN(16) * g(16) ! y yz
    mn(9) = locN(17) * g(17) - locN(18) * g(18) !-y yz

    m0= locN(19) * g(19)

    !calculate helper moments for transformation to moment space

    !transform current dist to moment space
    m(1) = rho
    m(2) = -11.0_rk*( mp(1) + mp(2) + mp(3) ) &
         + 8.0_rk*( mp(4) + mp(5) + mp(6) + mp(7) + mp(8) + mp(9) ) &
         - 30.0_rk*m0
    m(3) = -4.0_rk*( mp(1) + mp(2) + mp(3) ) &
         + ( mp(4) + mp(5) + mp(6) + mp(7) + mp(8) + mp(9) ) &
         + 12.0_rk*m0
    m(4) = u(1)
    m(5) = -4.0_rk*mn(1) + mp(4) + mp(5) - mp(6) - mp(7)
    m(6) = u(2)
    m(7) = -4.0_rk*mn(2) + mn(4) + mn(6) + mp(8) - mp(9)
    m(8) = u(3)
    m(9) = -4.0_rk*mn(3) + mn(5) + mn(7) + mn(8) + mn(9)
    m(10) = 2.0_rk*mp(1) - mp(2) - mp(3) +mp(4) + mp(5) + mp(6) + mp(7) &
         - 2.0_rk*( mp(8) + mp(9) )
    m(11) = -4.0_rk*mp(1) + 2.0_rk*( mp(2) + mp(3) ) + mp(4) + mp(5) + mp(6) + mp(7) &
         - 2.0_rk*( mp(8) + mp(9) )
    m(12) = mp(2) - mp(3) + mp(4) - mp(5) + mp(6) - mp(7)
    m(13) = 2.0_rk*( - mp(2) + mp(3) ) + mp(4) - mp(5) + mp(6) - mp(7)
    m(14) = mn(4) - mn(6)
    m(15) = mn(8) - mn(9)
    m(16) = mn(5) - mn(7)
    m(17) = mp(4) - mp(5) - mp(6) + mp(7)
    m(18) = - mn(4) - mn(6) + mp(8) - mp(9)
    m(19) = mn(5) + mn(7) - mn(8) - mn(9)

!    !calculate eq dist in moment space
!    meq(1) = rho
!    meq(2) = -11.0_rk*rho + 19.0_rk*jj/rho
!    meq(3) = 3.0_rk*rho - 5.50_rk*jj/rho  !!!!!!!!!!!!!!
!    meq(4) = jx
!    meq(5) = -2.0_rk/3.0_rk*jx + 2.50_rk*jx/rho2*( jy2 + jz2 )  !!!
!    meq(6) = jy
!    meq(7) = -2.0_rk/3.0_rk*jy + 2.50_rk*jy/rho2*( jz2 + jx2 )  !!!!
!    meq(8) = jz
!    meq(9) = -2.0_rk/3.0_rk*jz + 2.50_rk*jz/rho2*( jx2 + jy2 ) !!!!
!    meq(10) = ( 2.0_rk*jx2 - jy2 - jz2 )/rho
!    meq(11) = -meq(10)/2.0_rk !!!!!!!!!!!
!    meq(12) = ( jy2 - jz2 )/rho
!    meq(13) = -meq(12)/2.0_rk !!!!!!!!!!!
!    meq(14) = jx*jy/rho
!    meq(15) = jy*jz/rho
!    meq(16) = jx*jz/rho
!    meq(17) = 1.50_rk*jx/rho2*( jy2 - jz2 ) !!!!!!!!
!    meq(18) = 1.50_rk*jy/rho2*( jz2 - jx2 ) !!!!!!!!
!    meq(19) = 1.50_rk*jz/rho2*( jx2 - jy2 ) !!!!!!!!


    !calculate eq dist in moment space
    meq(1) = rho
    meq(2) = -11.0_rk*rho + 19.0_rk*jj
    meq(3) = 3.0_rk*rho - 5.50_rk*jj  !!!!!!!!!!!!!!
    meq(4) = jx
    meq(5) = -2.0_rk/3.0_rk*jx! + 2.50_rk*jx/rho2*( jy2 + jz2 )  !!!
    meq(6) = jy
    meq(7) = -2.0_rk/3.0_rk*jy! + 2.50_rk*jy/rho2*( jz2 + jx2 )  !!!!
    meq(8) = jz
    meq(9) = -2.0_rk/3.0_rk*jz! + 2.50_rk*jz/rho2*( jx2 + jy2 ) !!!!
    meq(10) = ( 2.0_rk*jx2 - jy2 - jz2 )!/rho
    meq(11) = -meq(10)/2.0_rk !!!!!!!!!!!
    meq(12) = ( jy2 - jz2 )!/rho
    meq(13) = -meq(12)/2.0_rk !!!!!!!!!!!
    meq(14) = jx*jy!/rho
    meq(15) = jy*jz!/rho
    meq(16) = jx*jz!/rho
    meq(17) = 0.0_rk !1.50_rk*jx/rho2*( jy2 - jz2 ) !!!!!!!!
    meq(18) = 0.0_rk !1.50_rk*jy/rho2*( jz2 - jx2 ) !!!!!!!!
    meq(19) = 0.0_rk !1.50_rk*jz/rho2*( jx2 - jy2 ) !!!!!!!!

#ifndef MRT_SINGLEPARTICLE
    !calculate eq dist to add external accelaration
    ! force implemmentation according to
    ! Phys. Rev. E 65,(2002): Guo
    ! S^(-1)*M*rho*( b*ci/cs**2 + b ueq*(cici - cs**2 I)/cs**4)
    ! ueq = u_tilde + f_{Shan Chen} + g_accn/2
    meqf(1) = -rho/3.0_rk*(3.0_rk/holdit - 1.0_rk)* &
         ( aaccn(1)*ueq(1) + aaccn(2)*ueq(2) + aaccn(3)*ueq(3) )*holdit**2/S(1)
    meqf(2) = rho/9.0_rk*(99.0_rk/holdit + 5.0_rk)* &
         ( aaccn(1)*ueq(1) + aaccn(2)*ueq(2) + aaccn(3)*ueq(3) )*holdit**2/S(2)
    meqf(3) = -rho/9.0_rk*(27.0_rk/holdit + 2.0_rk)* &
         ( aaccn(1)*ueq(1) + aaccn(2)*ueq(2) + aaccn(3)*ueq(3) )*holdit**2/S(3)
    meqf(4) = rho/3.0_rk*aaccn(1)*holdit/S(4)
    meqf(5) = -2.0_rk*rho/9.0_rk*aaccn(1)*holdit/S(5)
    meqf(6) = rho/3.0_rk*aaccn(2)*holdit/S(6)
    meqf(7) = -2.0_rk*rho/9.0_rk*aaccn(2)*holdit/S(7)
    meqf(8) = rho/3.0_rk*aaccn(3)*holdit/S(8)
    meqf(9) = -2.0_rk*rho/9.0_rk*aaccn(3)*holdit/S(9)
    meqf(10) = 2.0_rk*rho/9.0_rk*( 2.0_rk*aaccn(1)*ueq(1) - aaccn(2)*ueq(2) - aaccn(3)*ueq(3) )*holdit**2/S(10)
    meqf(11) = -rho/9.0_rk*( 2.0_rk*aaccn(1)*ueq(1) - aaccn(2)*ueq(2) - aaccn(3)*ueq(3) )*holdit**2/S(11)
    meqf(12) = 2.0_rk*rho/9.0_rk*( aaccn(2)*ueq(2) - aaccn(3)*ueq(3) )*holdit**2/S(12)
    meqf(13) = -rho/9.0_rk*( aaccn(2)*ueq(2) - aaccn(3)*ueq(3) )*holdit**2/S(13)
    meqf(14) = rho/9.0_rk*( aaccn(1)*ueq(2) + aaccn(2)*ueq(1) )*holdit**2/S(14)
    meqf(15) = rho/9.0_rk*( aaccn(2)*ueq(3) + aaccn(3)*ueq(2) )*holdit**2/S(15)
    meqf(16) = rho/9.0_rk*( aaccn(1)*ueq(3) + aaccn(3)*ueq(1) )*holdit**2/S(16)
    meqf(17) = 0.0_rk
    meqf(18) = 0.0_rk
    meqf(19) = 0.0_rk

    ! rescale in velocity field
    meqf(:) = (1.0_rk - 0.5_rk*S(:))*meqf(:)
#endif    

    !relax - don't want to come to it!
    !
#ifdef MRT_SINGLEPARTICLE    
    mnew(:) = S(:) * ( m(:) - meq(:))
#else    
    mnew(:) = S(:) * ( m(:) - meq(:) - meqf(:) )


    !!!REM PRINT HERE

    IF (kbT>0) THEN 
        fluctuating_amplitude=sqrt(rho*kbT*3.0_rk)
        call gaussianBoxMuller(gaussian1,gaussian2)
        mnew( 2)=mnew( 2)+fluctuating_amplitude*gaussian1*sqrt((402.0_rk)*			(1.0_rk-1.0_rk/S(2)**2.0_rk))
        mnew( 3)=mnew( 3)+fluctuating_amplitude*gaussian2*sqrt(59.0_rk)
        call gaussianBoxMuller(gaussian1,gaussian2)
        mnew( 5)=mnew( 5)+fluctuating_amplitude*gaussian1*sqrt(34.0_rk/9.0_rk)
        mnew( 7)=mnew( 7)+fluctuating_amplitude*gaussian2*sqrt(34.0_rk/9.0_rk)
        call gaussianBoxMuller(gaussian1,gaussian2)
        mnew( 9)=mnew( 9)+fluctuating_amplitude*gaussian1*sqrt(34.0_rk/9.0_rk)
        mnew(10)=mnew(10)+fluctuating_amplitude*gaussian2*sqrt((2.0_rk)*			(1.0_rk-1.0_rk/S(10)**2.0_rk))
        call gaussianBoxMuller(gaussian1,gaussian2)
        mnew(11)=mnew(11)+fluctuating_amplitude*gaussian1*sqrt(6.0_rk)
        mnew(12)=mnew(12)+fluctuating_amplitude*gaussian2*sqrt((2.0_rk/3.0_rk)*			(1.0_rk-1.0_rk/S(12)**2.0_rk))
        call gaussianBoxMuller(gaussian1,gaussian2)
        mnew(13)=mnew(13)+fluctuating_amplitude*gaussian1*sqrt(2.0_rk)
        mnew(14)=mnew(14)+fluctuating_amplitude*gaussian2*sqrt((  1.0_rk/ 9.0_rk)*		(1.0_rk-1.0_rk/S(14)**2.0_rk))
        call gaussianBoxMuller(gaussian1,gaussian2)
        mnew(15)=mnew(15)+fluctuating_amplitude*gaussian1*sqrt((  1.0_rk/ 9.0_rk)*		(1.0_rk-1.0_rk/S(15)**2.0_rk))
        mnew(16)=mnew(16)+fluctuating_amplitude*gaussian2*sqrt((  1.0_rk/ 9.0_rk)*		(1.0_rk-1.0_rk/S(16)**2.0_rk))
        call gaussianBoxMuller(gaussian1,gaussian2)
        mnew(17)=mnew(17)+fluctuating_amplitude*gaussian1*sqrt(  2.0_rk/ 9.0_rk)
        mnew(18)=mnew(18)+fluctuating_amplitude*gaussian2*sqrt(  2.0_rk/ 9.0_rk)
        call gaussianBoxMuller(gaussian1,gaussian2)
        mnew(19)=mnew(19)+fluctuating_amplitude*gaussian1*sqrt(  2.0_rk/ 9.0_rk)
    END IF
#endif    

! *sqrt(1085.0_rk)*	
! *sqrt(161.0_rk)*	
!                        
! *sqrt(2.0_rk)*		
! *sqrt(2.0_rk)*		
!                        
! *sqrt(2.0_rk)*		
! *sqrt(4.0_rk/3.0_rk)*	
!                        
! *sqrt(10.0_rk/9.0_rk)*	
! *sqrt(4.0_rk/9.0_rk)*	
!                        
! *sqrt(10.0_rk/9.0_rk)*	
! *sqrt(  1.0_rk/ 9.0_rk)
!                        
! *sqrt(  1.0_rk/ 9.0_rk)
! *sqrt(  1.0_rk/ 9.0_rk)
!                        
! *sqrt(  2.0_rk/ 9.0_rk)
! *sqrt(  2.0_rk/ 9.0_rk)
!                        
! *sqrt(  2.0_rk/ 9.0_rk)



! *sqrt((402.0_rk)*		(1.0_rk-1.0_rk/S(2)**2.0_rk))!sqrt(402.0_rk)
! *sqrt(59.0_rk)!sqrt(59.0_rk )
!                                                                                       
! *sqrt(34.0_rk/9.0_rk)!*sqrt( 34.0_rk/9.0_rk)    
! *sqrt(34.0_rk/9.0_rk)!*sqrt( 34.0_rk/9.0_rk)    
!                                                                                       
! *sqrt(34.0_rk/9.0_rk)!*sqrt(  34.0_rk/9.0_rk)    
! *sqrt((2.0_rk)*			(1.0_rk-1.0_rk/S(10)**2.0_rk))!*sqrt(  2.0_rk)
!                                                                                       
! *sqrt((6.0_rk)*			(1.0_rk-1.0_rk/S(10)**2.0_rk))!*sqrt( 6.0_rk)
! *sqrt((2.0_rk/3.0_rk)*		(1.0_rk-S(10)**2.0_rk))!*sqrt(  2.0_rk/ 3.0_rk)
!                                                                                       
! *sqrt((2.0_rk)*			(1.0_rk-1.0_rk/S(10)**2.0_rk))!*sqrt(  2.0_rk)
! *sqrt((  1.0_rk/ 9.0_rk)*		(1.0_rk-1.0_rk/S(10)**2.0_rk))
!                                                                                       
! *sqrt((  1.0_rk/ 9.0_rk)*		(1.0_rk-1.0_rk/S(10)**2.0_rk))
! *sqrt((  1.0_rk/ 9.0_rk)*		(1.0_rk-1.0_rk/S(10)**2.0_rk))
!                                                                                       
! *sqrt(  2.0_rk/ 9.0_rk)
! *sqrt(  2.0_rk/ 9.0_rk)
!                                                                                       
! *sqrt(  2.0_rk/ 9.0_rk)




    !calculate helpers for tranformation back to velocity space
    minv(1) = mnew(1)/19.0_rk - 11.0_rk*mnew(2)/2394.0_rk  - mnew(3)/63.0_rk
    minv(2) = mnew(1)/19.0_rk +  4.0_rk*mnew(2)/1197.0_rk  + mnew(3)/252.0_rk
    minv(3) = mnew(1)/19.0_rk -  5.0_rk*mnew(2)/399.0_rk   + mnew(3)/21.0_rk
    minv(4) = ( mnew(4) - mnew(5) )/10.0_rk
    minv(5) =   mnew(4)/10.0_rk + mnew(5)/40.0_rk
    minv(6) = ( mnew(6) - mnew(7) )/10.0_rk
    minv(7) =   mnew(6)/10.0_rk + mnew(7)/40.0_rk
    minv(8) = ( mnew(8) - mnew(9) )/10.0_rk
    minv(9) =   mnew(8)/10.0_rk + mnew(9)/40.0_rk
    minv(10) = ( mnew(10)  -  mnew(11) )/36.0_rk
    minv(11) =  mnew(10)/36.0_rk +  mnew(11)/72.0_rk
    minv(12) = ( mnew(12)  -  mnew(13) )/12.0_rk
    minv(13) =  mnew(12)/12.0_rk + mnew(13)/24.0_rk
    minv(14) = mnew(14)/4.0_rk
    minv(15) = mnew(15)/4.0_rk
    minv(16) = mnew(16)/4.0_rk
    minv(17) = mnew(17)/8.0_rk
    minv(18) = mnew(18)/8.0_rk
    minv(19) = mnew(19)/8.0_rk


    !do transformation to velocity space and apply to local occupation density
    locN(1)  = locN(1) - ( minv(1) + minv(4) + 2.0_rk*minv(10) )/g(1)
    locN(2)  = locN(2) - ( minv(1) - minv(4) + 2.0_rk*minv(10) )/g(2)
    locN(3)  = locN(3) - ( minv(1) + minv(6) - minv(10) + minv(12) )/g(3)
    locN(4)  = locN(4) - ( minv(1) - minv(6) - minv(10) + minv(12) )/g(4)
    locN(5)  = locN(5) - ( minv(1) + minv(8) - minv(10) - minv(12) )/g(5)
    locN(6)  = locN(6) - ( minv(1) - minv(8) - minv(10) - minv(12) )/g(6)
    locN(7)  = locN(7) - ( minv(2) + minv(5) + minv(7) + minv(11) + minv(13) + minv(14) + minv(17) - minv(18) )/g(7)
    locN(8)  = locN(8) - ( minv(2) + minv(5) - minv(7) + minv(11) + minv(13) - minv(14) + minv(17) + minv(18) )/g(8)
    locN(9)  = locN(9) - ( minv(2) + minv(5) + minv(9) + minv(11) - minv(13) + minv(16) - minv(17) + minv(19) )/g(9)
    locN(10) = locN(10) - ( minv(2) + minv(5) - minv(9) + minv(11) - minv(13) - minv(16) - minv(17) - minv(19) )/g(10)
    locN(11) = locN(11) - ( minv(2) - minv(5) + minv(7) + minv(11) + minv(13) - minv(14) - minv(17) - minv(18) )/g(11)
    locN(12) = locN(12) - ( minv(2) - minv(5) - minv(7) + minv(11) + minv(13) + minv(14) - minv(17) + minv(18) )/g(12)
    locN(13) = locN(13) - ( minv(2) - minv(5) + minv(9) + minv(11) - minv(13) - minv(16) + minv(17) + minv(19) )/g(13)
    locN(14) = locN(14) - ( minv(2) - minv(5) - minv(9) + minv(11) - minv(13) + minv(16) + minv(17) - minv(19) )/g(14)
    locN(15) = locN(15) - ( minv(2) + minv(7) + minv(9) - 2.0_rk*minv(11) + minv(15) - minv(19) + minv(18) )/g(15)
    locN(16) = locN(16) - ( minv(2) + minv(7) - minv(9) - 2.0_rk*minv(11) - minv(15) + minv(19) + minv(18) )/g(16)
    locN(17) = locN(17) - ( minv(2) - minv(7) + minv(9) - 2.0_rk*minv(11) - minv(15) - minv(19) - minv(18) )/g(17)
    locN(18) = locN(18) - ( minv(2) - minv(7) - minv(9) - 2.0_rk*minv(11) + minv(15) + minv(19) - minv(18) )/g(18)
    locN(19) = locN(19) - minv(3)/g(19)
    !...and done.
  end subroutine mrt_dist

  subroutine mrt_init_dist(ueq, mass, locN)
    !density
    real(kind=rk) :: rho, rho2, mass
    !velocity
    real(kind=rk), dimension(3) :: ueq
    real(kind=rk) :: jx,jy,jz
    real(kind=rk) :: jx2,jy2,jz2,jj

    !local occupation densities
    real(kind=rk), dimension(19) :: locN

    !moments
    real(kind=rk), dimension(19) :: meq,minv

    !local density
    rho = mass
    rho = max(rho,1.d-9)
    rho2 = rho**2

    !construct equilibrium density flux
    jx = ueq(1)*rho
    jy = ueq(2)*rho
    jz = ueq(3)*rho

    jx2 = jx**2
    jy2 = jy**2
    jz2 = jz**2
    jj  = jx2 + jy2 + jz2


!    !calculate eq dist in moment space
!    meq(1) = rho
!    meq(2) = -11.0_rk*rho + 19.0_rk*jj/rho
!    meq(3) = 3.0_rk*rho - 5.50_rk*jj/rho  !!!!!!!!!!!!!!
!    meq(4) = jx
!    meq(5) = -2.0_rk/3.0_rk*jx + 2.50_rk*jx/rho2*( jy2 + jz2 )  !!!
!    meq(6) = jy
!    meq(7) = -2.0_rk/3.0_rk*jy + 2.50_rk*jy/rho2*( jz2 + jx2 )  !!!!
!    meq(8) = jz
!    meq(9) = -2.0_rk/3.0_rk*jz + 2.50_rk*jz/rho2*( jx2 + jy2 ) !!!!
!    meq(10) = ( 2.0_rk*jx2 - jy2 - jz2 )/rho
!    meq(11) = -meq(10)/2.0_rk !!!!!!!!!!!
!    meq(12) = ( jy2 - jz2 )/rho
!    meq(13) = -meq(12)/2.0_rk !!!!!!!!!!!
!    meq(14) = jx*jy/rho
!    meq(15) = jy*jz/rho
!    meq(16) = jx*jz/rho
!    meq(17) = 1.50_rk*jx/rho2*( jy2 - jz2 ) !!!!!!!!
!    meq(18) = 1.50_rk*jy/rho2*( jz2 - jx2 ) !!!!!!!!
!    meq(19) = 1.50_rk*jz/rho2*( jx2 - jy2 ) !!!!!!!!


    !calculate eq dist in moment space
    meq(1) = rho
    meq(2) = -11.0_rk*rho + 19.0_rk*jj
    meq(3) = 3.0_rk*rho - 5.50_rk*jj  !!!!!!!!!!!!!!
    meq(4) = jx
    meq(5) = -2.0_rk/3.0_rk*jx! + 2.50_rk*jx/rho2*( jy2 + jz2 )  !!!
    meq(6) = jy
    meq(7) = -2.0_rk/3.0_rk*jy! + 2.50_rk*jy/rho2*( jz2 + jx2 )  !!!!
    meq(8) = jz
    meq(9) = -2.0_rk/3.0_rk*jz! + 2.50_rk*jz/rho2*( jx2 + jy2 ) !!!!
    meq(10) = ( 2.0_rk*jx2 - jy2 - jz2 )!/rho
    meq(11) = -meq(10)/2.0_rk !!!!!!!!!!!
    meq(12) = ( jy2 - jz2 )!/rho
    meq(13) = -meq(12)/2.0_rk !!!!!!!!!!!
    meq(14) = jx*jy!/rho
    meq(15) = jy*jz!/rho
    meq(16) = jx*jz!/rho
    meq(17) = 0.0_rk !1.50_rk*jx/rho2*( jy2 - jz2 ) !!!!!!!!
    meq(18) = 0.0_rk !1.50_rk*jy/rho2*( jz2 - jx2 ) !!!!!!!!
    meq(19) = 0.0_rk !1.50_rk*jz/rho2*( jx2 - jy2 ) !!!!!!!!

    !calculate helpers for tranformation back to velocity space
    minv(1) = meq(1)/19.0_rk - 11.0_rk*meq(2)/2394.0_rk  - meq(3)/63.0_rk
    minv(2) = meq(1)/19.0_rk +  4.0_rk*meq(2)/1197.0_rk  + meq(3)/252.0_rk
    minv(3) = meq(1)/19.0_rk -  5.0_rk*meq(2)/399.0_rk   + meq(3)/21.0_rk
    minv(4) = ( meq(4) - meq(5) )/10.0_rk
    minv(5) =   meq(4)/10.0_rk + meq(5)/40.0_rk
    minv(6) = ( meq(6) - meq(7) )/10.0_rk
    minv(7) =   meq(6)/10.0_rk + meq(7)/40.0_rk
    minv(8) = ( meq(8) - meq(9) )/10.0_rk
    minv(9) =   meq(8)/10.0_rk + meq(9)/40.0_rk
    minv(10) = ( meq(10)  -  meq(11) )/36.0_rk
    minv(11) =  meq(10)/36.0_rk +  meq(11)/72.0_rk
    minv(12) = ( meq(12)  -  meq(13) )/12.0_rk
    minv(13) =  meq(12)/12.0_rk + meq(13)/24.0_rk
    minv(14) = meq(14)/4.0_rk
    minv(15) = meq(15)/4.0_rk
    minv(16) = meq(16)/4.0_rk
    minv(17) = meq(17)/8.0_rk
    minv(18) = meq(18)/8.0_rk
    minv(19) = meq(19)/8.0_rk

    !do transformation to velocity space and apply to local occupation density
    locN(1)  =  ( minv(1) + minv(4) + 2.0_rk*minv(10) )/g(1)
    locN(2)  =  ( minv(1) - minv(4) + 2.0_rk*minv(10) )/g(2)
    locN(3)  =  ( minv(1) + minv(6) - minv(10) + minv(12) )/g(3)
    locN(4)  =  ( minv(1) - minv(6) - minv(10) + minv(12) )/g(4)
    locN(5)  =  ( minv(1) + minv(8) - minv(10) - minv(12) )/g(5)
    locN(6)  =  ( minv(1) - minv(8) - minv(10) - minv(12) )/g(6)
    locN(7)  =  ( minv(2) + minv(5) + minv(7) + minv(11) + minv(13) + minv(14) + minv(17) - minv(18) )/g(7)
    locN(8)  =  ( minv(2) + minv(5) - minv(7) + minv(11) + minv(13) - minv(14) + minv(17) + minv(18) )/g(8)
    locN(9)  =  ( minv(2) + minv(5) + minv(9) + minv(11) - minv(13) + minv(16) - minv(17) + minv(19) )/g(9)
    locN(10) =  ( minv(2) + minv(5) - minv(9) + minv(11) - minv(13) - minv(16) - minv(17) - minv(19) )/g(10)
    locN(11) =  ( minv(2) - minv(5) + minv(7) + minv(11) + minv(13) - minv(14) - minv(17) - minv(18) )/g(11)
    locN(12) =  ( minv(2) - minv(5) - minv(7) + minv(11) + minv(13) + minv(14) - minv(17) + minv(18) )/g(12)
    locN(13) =  ( minv(2) - minv(5) + minv(9) + minv(11) - minv(13) - minv(16) + minv(17) + minv(19) )/g(13)
    locN(14) =  ( minv(2) - minv(5) - minv(9) + minv(11) - minv(13) + minv(16) + minv(17) - minv(19) )/g(14)
    locN(15) =  ( minv(2) + minv(7) + minv(9) - 2.0_rk*minv(11) + minv(15) - minv(19) + minv(18) )/g(15)
    locN(16) =  ( minv(2) + minv(7) - minv(9) - 2.0_rk*minv(11) - minv(15) + minv(19) + minv(18) )/g(16)
    locN(17) =  ( minv(2) - minv(7) + minv(9) - 2.0_rk*minv(11) - minv(15) - minv(19) - minv(18) )/g(17)
    locN(18) =  ( minv(2) - minv(7) - minv(9) - 2.0_rk*minv(11) + minv(15) + minv(19) - minv(18) )/g(18)
    locN(19) =  minv(3)/g(19)

    !...and done.
  end subroutine mrt_init_dist

! subroutine mrt_dist_old(ueq, locN, S)
!         !density
! 	real*8 :: rho, rho2
!         !velocity
!         real*8, dimension(3) :: u,ueq
!         real*8 :: jx,jy,jz
!         real*8 :: jx2,jy2,jz2,jj

!         integer :: x,y,z

!         !local occupation densities
!         real*8, dimension(19) :: locN

!         !collision "matrix" is diag so need a vec only
!         real*8, dimension(19) :: S
!         !moments
!         real*8, dimension(19) :: m,meq,mnew
!         !helper moments
!         real*8, dimension(9)  :: mp,mn
!         real*8, dimension(19) :: minv
! 	real*8 :: m0

!         !local density
!         rho = sum( locN(:) *g(:) )
!         rho = max(rho,1.d-19)
!         rho2 = rho**2

!         !local velocities
!         u(1) = sum( locN(:) * cx(:) * g(:) )
!         u(2) = sum( locN(:) * cy(:) * g(:) )
!         u(3) = sum( locN(:) * cz(:) * g(:) )

!         !construct equilibrium density flux
!         jx = ueq(1) * rho
!         jy = ueq(2) * rho
!         jz = ueq(3) * rho

!         jx2 = jx**2
!         jy2 = jy**2
!         jz2 = jz**2
!         jj  = jx2 + jy2 + jz2

!         !calculate helper moments for transformation to moment space
!         mp(1) = locN(1)  * g(1)  + locN(2)  * g(2)
!         mp(2) = locN(3)  * g(3)  + locN(4)  * g(4)
!         mp(3) = locN(5)  * g(5)  + locN(6)  * g(6)
!         mp(4) = locN(7)  * g(7)  + locN(8)  * g(8)
!         mp(5) = locN(9)  * g(9)  + locN(10) * g(10)
!         mp(6) = locN(11) * g(11) + locN(12) * g(12)
!         mp(7) = locN(13) * g(13) + locN(14) * g(14)
!         mp(8) = locN(15) * g(15) + locN(16) * g(16)
!         mp(9) = locN(17) * g(17) + locN(18) * g(18)

!         mn(1) = locN(1)  * g(1)  - locN(2)  * g(2)
!         mn(2) = locN(3)  * g(3)  - locN(4)  * g(4)
!         mn(3) = locN(5)  * g(5)  - locN(6)  * g(6)
!         mn(4) = locN(7)  * g(7)  - locN(8)  * g(8)
!         mn(5) = locN(9)  * g(9)  - locN(10) * g(10)
!         mn(6) = locN(11) * g(11) - locN(12) * g(12)
!         mn(7) = locN(13) * g(13) - locN(14) * g(14)
!         mn(8) = locN(15) * g(15) - locN(16) * g(16)
!         mn(9) = locN(17) * g(17) - locN(18) * g(18)

!         m0= locN(19) * g(19)

!         !calculate helper moments for transformation to moment space

!         !transform current dist to moment space
!         m(1) = rho
!         m(2) = -11.d0 * ( mp(1) + mp(2) + mp(3) ) + 8.d0 * ( mp(4) + mp(5) + mp(6) + mp(7) + mp(8) + mp(9) ) - 30.d0 * m0
!         m(3) = -4.d0  * ( mp(1) + mp(2) + mp(3) ) + ( mp(4) + mp(5) + mp(6) + mp(7) + mp(8) + mp(9) ) + 12.d0 * m0
!         m(4) = u(1)
!         m(5) = -4.d0  * mn(1) + mp(4) + mp(5) - mp(6) - mp(7)
!         m(6) = u(2)
!         m(7) = -4.d0  * mn(2) + mn(4) + mn(6) + mp(8) - mp(9)
!         m(8) = u(3)
!         m(9) = -4.d0  * mn(3) + mn(5) + mn(7) + mn(8) + mn(9)
!         m(10) = 2.d0 * mp(1) - mp(2) - mp(3) +mp(4) + mp(5) + mp(6) + mp(7) - 2.d0 * ( mp(8) + mp(9) )
!         m(11) = -4.d0 * mp(1) + 2.d0 * ( mp(2) + mp(3) ) + mp(4) + mp(5) + mp(6) + mp(7) - 2.d0 * ( mp(8) + mp(9) )
!         m(12) = mp(2) - mp(3) + mp(4) - mp(5) + mp(6) - mp(7)
!         m(13) = 2.d0 * ( - mp(2) + mp(3) ) + mp(4) - mp(5) + mp(6) - mp(7)
!         m(14) = mn(4) - mn(6)
!         m(15) = mn(8) - mn(9)
!         m(16) = mn(5) - mn(7)
!         m(17) = mp(4) - mp(5) - mp(6) + mp(7)
!         m(18) = - mn(4) - mn(6) + mp(8) - mp(9)
!         m(19) = mn(5) + mn(7) - mn(8) - mn(9)

!         !calculate eq dist in moment space
!         meq(1) = rho
!         meq(2) = ( -11.d0 * rho ) + ( (19.d0 / rho)  * jj )
!         meq(3) = ( 3.d0 * rho ) - ( (11.d0 / ( 2.d0 * rho ) ) * jj )
!         meq(4) = jx
!         meq(5) = jx * ( -2.d0 / 3.d0 ) + jx * ( 15.d0 / ( 6.d0*rho2 ) ) * ( jy2 + jz2 )
!         meq(6) = jy
!         locN(12) = locN(12) - ( minv(2) - minv(5) - minv(7) + minv(11) + minv(13) + minv(14) - minv(17) + minv(18) ) / g(12)
!         locN(13) = locN(13) - ( minv(2) - minv(5) + minv(9) + minv(11) - minv(13) - minv(16) + minv(17) + minv(19) ) / g(13)
!         locN(14) = locN(14) - ( minv(2) - minv(5) - minv(9) + minv(11) - minv(13) + minv(16) + minv(17) - minv(19) ) / g(14)
!         locN(15) = locN(15) - ( minv(2) + minv(7) + minv(9) - ( 2.d0 * minv(11) ) + minv(15) - minv(19) + minv(18) ) / g(15)
!         locN(16) = locN(16) - ( minv(2) + minv(7) - minv(9) - ( 2.d0 * minv(11) ) - minv(15) + minv(19) + minv(18) ) / g(16)
!         locN(17) = locN(17) - ( minv(2) - minv(7) + minv(9) - ( 2.d0 * minv(11) ) - minv(15) - minv(19) - minv(18) ) / g(17)
!         locN(18) = locN(18) - ( minv(2) - minv(7) - minv(9) - ( 2.d0 * minv(11) ) + minv(15) + minv(19) - minv(18) ) / g(18)
!         locN(19) = locN(19) - ( minv(3) ) / g(19)

!         !...and done.
!         return

! end subroutine mrt_dist_old


!> Compute equilibrium populations derived from fluid density and velocity
!>
!> added by Timm, 07 May 2012
!> The equilibrium distributions are computed from the fluid density and velocity.
!> The distributions are written into the array \c eq_dist which has to be allocated before.
!> This subroutine has to be called once for each lattice node and for each fluid component.
!> The passed density \c den and velocity \c vel are the fluid (component) density
!> and the corrected fluid (component) velocity, respectively.
!> Note that the user has to make sure that the velocity already contains the force correction.
!> This is an experimental subroutine which has to be extended and made for efficient in the future.
!> Ideally, this subroutine should replace the existing ones eventually.

subroutine equilibrium_distribution(den, vel, eq_dist)
  real(kind=rk), intent(in) :: den ! density of fluid (component)
  real(kind=rk), dimension(3), intent(in) :: vel ! true equilibrium velocity of fluid (component)
  real(kind=rk), dimension(nvecs), intent(out) :: eq_dist ! computed equilibrium distribution for fluid (component)

  ! Declare variables.
  integer :: i ! index
  real(kind=rk) :: w0, w1, w2 ! lattice weights
  real(kind=rk), dimension(nvecs) :: weight ! lattice weights
  real(kind=rk) :: cs2, cs4 ! sound speed^2, sound speed^4

  ! Set variables.
  w0 = 1.0_rk / 3.00_rk
  w1 = 1.0_rk / 18.00_rk
  w2 = 1.0_rk / 36.00_rk
  weight = (/w1, w1, w1, w1, w1, w1, w2, w2, w2, w2, w2, w2, w2, w2, w2, w2, w2, w2, w0/)
  cs2 = 1.0_rk / 3.00_rk
  cs4 = 1.0_rk / 9.00_rk

  ! Loop over all directions and compute the equilibrium distributions.
  do i = 1, nvecs
    eq_dist(i) = den * weight(i) * &
      & (1.0_rk + dot_product(vel(:), c(i, :)) / cs2 &
      & + (dot_product(vel(:), c(i, :)))**2 / (2.0_rk * cs4) &
      & - dot_product(vel(:), vel(:)) / (2.0_rk * cs2))
  end do
end subroutine equilibrium_distribution

  !> uses the Box Mueller transform to generate to Gaussian random numbers with
  !> mean value 0 and variance 1
  subroutine gaussianBoxMuller(gaussian1,gaussian2)
    implicit none

    real(kind=rk) :: gaussian1, gaussian2
    real(kind=rk) :: u1,u2
    real(kind=rk),dimension(2) :: u

!!  With the better ranlux random numbers
    call ranlux(u,2)
    u1=u(1)
    u2=u(2)

!! Without the more expensive ranlux
!    call random_number(u1)
!    call random_number(u2)

    gaussian1 = sqrt(-2.0_rk*log(u1)) * cos(2.0_rk*pi*u2)
    gaussian2 = sqrt(-2.0_rk*log(u1)) * sin(2.0_rk*pi*u2)

!! Clipping to 6 sigma for stability purposes
   gaussian1=mod(gaussian1,6.0_rk)
   gaussian2=mod(gaussian2,6.0_rk)
  end subroutine gaussianBoxMuller

end module lbe_bdist_module
