#include "lbe.h"

!> Code for multiple relaxation time collision operator
module lbe_MRT

  use lbe_globals_module, only: rk,pi
  use lbe_parms_module, only: kBT
  use luxury_rng_module, only : ranlux

  implicit none
  private

  public lbe_MRTdist

contains

 subroutine lbe_MRTdist(nd,mass,tau,tau_bulk,f,dist)
  
  implicit none

  !type(lbe_site),dimension(0:,0:,0:)  :: N
  real(kind=rk), dimension(19), intent(in)  :: nd
  real(kind=rk),                intent(in)  :: mass
  real(kind=rk),                intent(in)  :: tau, tau_bulk
  real(kind=rk), dimension(3),  intent(in)  :: f
  real(kind=rk), dimension(19), intent(out) :: dist

  real(kind=rk), dimension(19) :: m, m_eq, delta
  real(kind=rk), dimension(3)  :: u
  real(kind=rk) :: rho
  real(kind=rk) :: fluctuating_amplitude
  real(kind=rk) :: gaussian1, gaussian2

  ! Calculating moments from number densities

  ! rho
  m(1) =  nd(19) +  nd(1) +  nd(2) +  nd(3) +  nd(4) + &
       nd(5) +  nd(6) +  nd(7) +  nd(11) +  nd(8) + &
       nd(12) +  nd(9) +  nd(13) +  nd(10) +  nd(14) +  nd(15) + &
       nd(17) +  nd(16) +  nd(18)
  ! e
  m(2) = -30 * nd(19) - 11 * nd(1) - 11 * nd(2) - 11 * nd(3) - 11 * nd(4) - &
       11 * nd(5) - 11 * nd(6) + 8 * nd(7) + 8 * nd(11) + 8 * nd(8) + 8 * nd(12) + 8 * &
       nd(9) + 8 * nd(13) + 8 * nd(10) + 8 * nd(14) + 8 * nd(15) + 8 * nd(17) + 8 * nd(16) + 8 * nd(18)
  ! epsilon
  m(3) = 12 * nd(19) - 4 * nd(1) - 4 * nd(2) - 4 * nd(3) - 4 * nd(4) - 4 * &
       nd(5) - 4 * nd(6) + nd(7) + nd(11) + nd(8) + nd(12) + nd(9) + nd(13) + nd(10) + nd(14) + &
       nd(15) + nd(17) + nd(16) + nd(18)
  ! j_x
  m(4) = nd(1) - nd(2) + nd(7) - nd(11) + nd(8) - nd(12) + nd(9) - nd(13) + nd(10) - nd(14)
  ! q_x
  m(5) = -4 * nd(1) + 4 * nd(2) + nd(7) - nd(11) + nd(8) - nd(12) + nd(9) - nd(13) + nd(10) - nd(14)
  ! j_y
  m(6) = nd(3) - nd(4) + nd(7) + nd(11) - nd(8) - nd(12) + nd(15) - nd(17) + nd(16) - nd(18)
  ! q_y
  m(7) = -4 * nd(3) + 4 * nd(4) + nd(7) + nd(11) - nd(8) - nd(12) + nd(15) - nd(17) + nd(16) - nd(18)
  ! j_z
  m(8) = nd(5) - nd(6) + nd(9) + nd(13) - nd(10) - nd(14) + nd(15) + nd(17) - nd(16) - nd(18)
  ! q_z
  m(9) = -4 * nd(5) + 4 * nd(6) + nd(9) + nd(13) - nd(10) - nd(14) + nd(15) + nd(17) - nd(16) - nd(18)
  ! 3p_xx
  m(10) = 2 * nd(1) + 2 * nd(2) - nd(3) - nd(4) - nd(5) - nd(6) + nd(7) + nd(11) + nd(8) + nd(12) + &
       nd(9) + nd(13) + nd(10) + nd(14) - 2 * nd(15) - 2 * nd(17) - 2 * nd(16) - 2 * nd(18)
  ! 3pi_xx
  m(11) = -4 * nd(1) - 4 * nd(2) + 2 * nd(3) + 2 * nd(4) + 2 * nd(5) + 2 * nd(6) + nd(7) + nd(11) + &
       nd(8) + nd(12) + nd(9) + nd(13) + nd(10) + nd(14) - 2 * nd(15) - 2 * nd(17) - 2 * nd(16) - 2 * nd(18)
  ! p_ww
  m(12) = nd(3) + nd(4) - nd(5) - nd(6) + nd(7) + nd(11) + nd(8) + nd(12) - nd(9) - nd(13) - nd(10) - nd(14)
  ! pi_ww
  m(13) = -2 * nd(3) - 2 * nd(4) + 2 * nd(5) + 2 * nd(6) + nd(7) + nd(11) + nd(8) + nd(12) - nd(9) - nd(13) - nd(10) - nd(14)
  ! p_xy
  m(14) = nd(7) - nd(11) - nd(8) + nd(12)
  ! p_yz
  m(15) = nd(15) - nd(17) - nd(16) + nd(18)
  ! p_xz
  m(16) = nd(9) - nd(13) - nd(10) + nd(14)
  ! m_x
  m(17) = nd(7) - nd(11) + nd(8) - nd(12) - nd(9) + nd(13) - nd(10) + nd(14)
  ! m_y
  m(18) = -nd(7) - nd(11) + nd(8) + nd(12) + nd(15) - nd(17) + nd(16) - nd(18)
  ! m_y
  m(19) = nd(9) + nd(13) - nd(10) - nd(14) - nd(15) - nd(17) + nd(16) + nd(18)


  rho = m(1)

  ! Calculating mean velocity and adding Shan Chen forces
  u(1) = (nd(1)-nd(2)+nd(7)-nd(12)+nd(8)-nd(11)+nd(9)-nd(14)+nd(10)-nd(13)) &
       &+ ((( 1.d0 / 10.d0 * tau_bulk) + (9.d0 / 10.d0 * tau)) * f(1) / mass )
  u(2) = (nd(3)-nd(4)+nd(7)+nd(11)-nd(8)-nd(12)+nd(15)-nd(17)+nd(16)-nd(18))+ ( tau * f(2) / mass )
  u(3) = (nd(5)-nd(6)+nd(9)+nd(13)-nd(10)-nd(14)+nd(15)+nd(17)-nd(16)-nd(18))+ ( tau * f(3) / mass )

  ! Calcualting equilibrium momenta
  m_eq(1) =  rho;
  m_eq(2) = -11.d0*rho+19.d0/rho*(u(1)**2+u(2)**2+u(3)**2)
  m_eq(3) = -475.d0/63.d0/rho*(u(1)**2+u(2)**2+u(3)**2)
  m_eq(4) = u(1)
  m_eq(5) = -2.d0/3.d0*u(1)
  m_eq(6) = u(2)
  m_eq(7) = -2.d0/3.d0*u(2)
  m_eq(8) = u(3)
  m_eq(9) = -2.d0/3.d0*u(3)
  m_eq(10) = 1.d0/rho*(2*u(1)**2-(u(2)**2+u(3)**2))
  m_eq(11) = 0.d0
  m_eq(12) = 1.d0/rho*(u(2)**2-u(3)**2)
  m_eq(13) = 0.d0
  m_eq(14) = 1.d0/rho*u(1)*u(2)
  m_eq(15) = 1.d0/rho*u(2)*u(3)
  m_eq(16) = 1.d0/rho*u(1)*u(3)
  m_eq(17) = 0.d0
  m_eq(18) = 0.d0
  m_eq(19) = 0.d0


  ! md_mdeq is the result of the expresion
  ! m-meq, where m are the momemtums of distribution 
  ! function n (ND(x,y,z)%n_r)
  ! meq is the equilibrium part of m
  ! m = M·n, M is the momentum matrix (19,19)
  ! meq = is a function of the macroscopic variables
  ! density,momentum

  !md_mdeq(1) =  0.0d0


  ! the advection step is
  ! nd(x+edt,t+dt)=n(x,t)-M^{-1}·SS·(m-meq)
  ! where SS is a diagonal matrix, and it has 
  ! the information about the relaxation times

  ! the variable boltz_r = -M^{-1}·SS·(m-meq) 

  delta = m - m_eq

  dist(19) = 17.0d0/1140.0d0*delta(2) - delta(3)/15.0d0

  dist(1) = 187.0d0/68400.0d0*delta(2) + delta(3)/90.0d0 +      &
       3.0d0/50.0d0*delta(5) - (1.d0/tau)*delta(10)/36.0d0 +           &
       7.0d0/180.0d0*delta(11)

  dist(2) = 187.0d0/68400.0d0*delta(2) + delta(3)/90.0d0 -      &
       3.0d0/50.0d0*delta(5) - (1.d0/tau)*delta(10)/36.0d0 +           &
       7.0d0/180.0d0*delta(11) 

  dist(3) = 187.0d0/68400.0d0*delta(2) + delta(3)/90.0d0 +      &
       3.0d0/50.0d0*delta(7) + (1.d0/tau)*delta(10)/72.0d0 -           &
       7.0d0/360.0d0*delta(11) - (1.d0/tau)*delta(12)/24.0d0 +         &
       7.0d0/120.0d0*delta(13) 

  dist(4) = 187.0d0/68400.0d0*delta(2) + delta(3)/90.0d0 -      &
       3.0d0/50.0d0*delta(7) + (1.d0/tau)*delta(10)/72.0d0 -           &
       7.0d0/360.0d0*delta(11) - (1.d0/tau)*delta(12)/24.0d0 +         &
       7.0d0/120.0d0*delta(13) 

  dist(5) = 187.0d0/68400.0d0*delta(2) + delta(3)/90.0d0 +      & 
       3.0d0/50.0d0*delta(9) + (1.d0/tau)*delta(10)/72.0d0 -           &
       7.0d0/360.0d0*delta(11) + (1.d0/tau)*delta(12)/24.0d0 -         &
       7.0d0/120.0d0*delta(13) 

  dist(6) = 187.0d0/68400.0d0*delta(2) + delta(3)/90.0d0 -      &
       3.0d0/50.0d0*delta(9) + (1.d0/tau)*delta(10)/72.0d0 -           &
       7.0d0/360.0d0*delta(11) + (1.d0/tau)*delta(12)/24.0d0 -         &
       7.0d0/120.0d0*delta(13)

  dist(7)  = -17.0d0/4275.0d0*delta(2) - delta(3)/180.0d0 -     &
       3.0d0/100.0d0*delta(5) - 3.0d0/100.0d0*delta(7) -            &
       (1.d0/tau)*delta(10)/36.0d0 - 7.0d0/360.0d0*delta(11) -         &
       (1.d0/tau)*delta(12)/12.0d0 - 7.0d0/120.0d0*delta(13) -         &
       (1.d0/tau)*delta(14)/4.0d0 - 99.0d0/400.0d0*delta(17) +         &
       99.0d0/400.0d0*delta(18)

  dist(11) = -17.0d0/4275.0d0*delta(2) - delta(3)/180.0d0 +     &
       3.0d0/100.0d0* delta(5) - 3.0d0/100.0d0*delta(7) -           &
       (1.d0/tau)*delta(10)/36.0d0 - 7.0d0/360.0d0*delta(11) -         &
       (1.d0/tau)*delta(12)/12.0d0 - 7.0d0/120.0d0*delta(13) +         &
       (1.d0/tau)*delta(14)/4.0d0 + 99.0d0/400.0d0*delta(17) +         &
       99.0d0/400.0d0*delta(18)

  dist(8) = -17.0d0/4275.0d0*delta(2) - delta(3)/180.0d0 -      &
       3.0d0/100.0d0*delta(5) + 3.0d0/100.0d0*delta(7) -            &
       (1.d0/tau)*delta(10)/36.0d0 - 7.0d0/360.0d0*delta(11) -         &
       (1.d0/tau)*delta(12)/12.0d0 - 7.0d0/120.0d0*delta(13) +         &
       (1.d0/tau)*delta(14)/4.0d0 - 99.0d0/400.0d0*delta(17) -         &
       99.0d0/400.0d0*delta(18)

  dist(12) = -17.0d0/4275.0d0*delta(2) - delta(3)/180.0d0 +     &
       3.0d0/100.0d0* delta(5) + 3.0d0/100.0d0*delta(7) -           &
       (1.d0/tau)*delta(10)/36.0d0 - 7.0d0/360.0d0*delta(11) -         &
       (1.d0/tau)*delta(12)/12.0d0 - 7.0d0/120.0d0*delta(13) -         &
       (1.d0/tau)*delta(14)/4.0d0 + 99.0d0/400.0d0*delta(17) -         &
       99.0d0/400.0d0*delta(18)

  dist(9) = -17.0d0/4275.0d0*delta(2) - delta(3)/180.0d0 -      &
       3.0d0/100.0d0* delta(5) - 3.0d0/100.0d0*delta(9) -           &
       (1.d0/tau)*delta(10)/36.0d0 - 7.0d0/360.0d0*delta(11) +         &
       (1.d0/tau)*delta(12)/12.0d0 + 7.0d0/120.0d0*delta(13) -         &
       (1.d0/tau)*delta(16)/4.0d0 + 99.0d0/400.0d0*delta(17) -         &
       99.0d0/400.0d0*delta(19)

  dist(13) = -17.0d0/4275.0d0*delta(2) - delta(3)/180.0d0 +     &
       3.0d0/100.0d0* delta(5) - 3.0d0/100.0d0*delta(9) -           &
       (1.d0/tau)*delta(10)/36.0d0 - 7.0d0/360.0d0*delta(11) +         &
       (1.d0/tau)*delta(12)/12.0d0 + 7.0d0/120.0d0*delta(13) +         &
       (1.d0/tau)*delta(16)/4.0d0 - 99.0d0/400.0d0*delta(17) -         &
       99.0d0/400.0d0*delta(19)

  dist(10) = -17.0d0/4275.0d0*delta(2) - delta(3)/180.0d0 -      &
       3.0d0/100.0d0* delta(5) + 3.0d0/100.0d0*delta(9) -           &
       (1.d0/tau)*delta(10)/36.0d0 - 7.0d0/360.0d0*delta(11) +         &
       (1.d0/tau)*delta(12)/12.0d0 + 7.0d0/120.0d0*delta(13) +         &
       (1.d0/tau)*delta(16)/4.0d0 + 99.0d0/400.0d0*delta(17) +         &
       99.0d0/400.0d0*delta(19)

  dist(14) = -17.0d0/4275.0d0*delta(2) - delta(3)/180.0d0 +     &
       3.0d0/100.0d0*delta(5) + 3.0d0/100.0d0*delta(9) -            &
       (1.d0/tau)*delta(10)/36.0d0 - 7.0d0/360.0d0*delta(11) +         &
       (1.d0/tau)*delta(12)/12.0d0 + 7.0d0/120.0d0*delta(13) -         &
       (1.d0/tau)*delta(16)/4.0d0 - 99.0d0/400.0d0*delta(17) +         &
       99.0d0/400.0d0*delta(19)

  dist(15) = -17.0d0/4275.0d0*delta(2) - delta(3)/180.0d0 -     &
       3.0d0/100.0d0*delta(7) - 3.0d0/100.0d0*delta(9) +            &
       (1.d0/tau)*delta(10)/18.0d0 + 7.0d0/180.0d0*delta(11) -         &
       (1.d0/tau)*delta(15)/4.0d0 - 99.0d0/400.0d0*delta(18) +         &
       99.0d0/400.0d0*delta(19)

  dist(17) = -17.0d0/4275.0d0*delta(2) - delta(3)/180.0d0 +     &
       3.0d0/100.0d0*delta(7) - 3.0d0/100.0d0*delta(9) +            &
       (1.d0/tau)*delta(10)/18.0d0 + 7.0d0/180.0d0*delta(11) +         &
       (1.d0/tau)*delta(15)/4.0d0 + 99.0d0/400.0d0*delta(18) +         &
       99.0d0/400.0d0*delta(19)

  dist(16) = -17.0d0/4275.0d0*delta(2) - delta(3)/180.0d0 -     &
       3.0d0/100.0d0*delta(7) + 3.0d0/100.0d0*delta(9) +            &
       (1.d0/tau)*delta(10)/18.0d0 + 7.0d0/180.0d0*delta(11) +         &
       (1.d0/tau)*delta(15)/4.0d0 - 99.0d0/400.0d0*delta(18) -         &
       99.0d0/400.0d0*delta(19)

  dist(18) = -17.0d0/4275.0d0*delta(2) - delta(3)/180.0d0 +     &
       3.0d0/100.0d0*delta(7) + 3.0d0/100.0d0*delta(9) +            &
       (1.d0/tau)*delta(10)/18.0d0 + 7.0d0/180.0d0*delta(11) -         &
       (1.d0/tau)*delta(15)/4.0d0 + 99.0d0/400.0d0*delta(18) -         &
       99.0d0/400.0d0*delta(19)

end subroutine lbe_MRTdist

end module lbe_MRT
