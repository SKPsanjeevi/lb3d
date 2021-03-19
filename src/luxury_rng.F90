#include "lbe.h"

!> Random number generator
!>
!> Subtract-and-borrow random number generator proposed by
!> Marsaglia and Zaman, implemented by F. James with the name
!> RCARRY in 1991, and later improved by Martin Luescher
!> in 1993 to produce "Luxury Pseudorandom Numbers".
!> Fortran 77 coded by F. James, 1993
!>
!> \verbatim
!>  References:
!>  M. Luscher, Computer Physics Communications  79 (1994) 100
!>  F. James, Computer Physics Communications 79 (1994) 111
!>
!>   LUXURY LEVELS.
!>   ------ ------      The available luxury levels are:
!>
!>  level 0  (p=24): equivalent to the original RCARRY of Marsaglia
!>           and Zaman, very long period, but fails many tests.
!>  level 1  (p=48): considerable improvement in quality over level 0,
!>           now passes the gap test, but still fails spectral test.
!>  level 2  (p=97): passes all known tests, but theoretically still
!>           defective.
!>  level 3  (p=223): DEFAULT VALUE.  Any theoretically possible
!>           correlations have very small chance of being observed.
!>  level 4  (p=389): highest possible luxury, all 24 bits chaotic.
!>
!> ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!>  Calling sequences for RANLUX:                                  ++
!>      CALL RANLUX (RVEC, LEN)   returns a vector RVEC of LEN     ++
!>                   32-bit random floating point numbers between  ++
!>                   zero (not included) and one (also not incl.). ++
!>      CALL RLUXGO(LUX,INT,K1,K2) initializes the generator from  ++
!>               one 32-bit integer INT and sets Luxury Level LUX  ++
!>               which is integer between zero and MAXLEV, or if   ++
!>               LUX .GT. 24, it sets p=LUX directly.  K1 and K2   ++
!>               should be set to zero unless restarting at a break++
!>               point given by output of RLUXAT (see RLUXAT).     ++
!>      CALL RLUXAT(LUX,INT,K1,K2) gets the values of four integers++
!>               which can be used to restart the RANLUX generator ++
!>               at the current point by calling RLUXGO.  K1 and K2++
!>               specify how many numbers were generated since the ++
!>               initialization with LUX and INT.  The restarting  ++
!>               skips over  K1+K2*E9   numbers, so it can be long.++
!>   A more efficient but less convenient way of restarting is by: ++
!>      CALL RLUXIN(ISVEC)    restarts the generator from vector   ++
!>                   ISVEC of 25 32-bit integers (see RLUXUT)      ++
!>      CALL RLUXUT(ISVEC)    outputs the current values of the 25 ++
!>                 32-bit integer seeds, to be used for restarting ++
!>      ISVEC must be dimensioned 25 in the calling program        ++
!> ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> \endverbatim
module luxury_rng_module

  use lbe_log_module

  implicit none
  private

  integer            :: iseeds(24), isdext(25)
  integer, parameter :: maxlev = 4, lxdflt = 3
  integer            :: ndskip(0:maxlev) = (/ 0, 24, 73, 199, 365 /)
  integer            :: next(24), igiga = 1000000000, jsdflt = 314159265
  real*8, parameter  :: twop12 = 4096.d0
  integer, parameter :: itwo24 = 2**24, icons = 2147483563
  integer            :: luxlev = lxdflt, nskip, inseed, jseed
  logical            :: notyet = .true.
  integer            :: in24 = 0, kount = 0, mkount = 0, i24 = 24, j24 = 10
  real*8             :: seeds(24), carry = 0.d0, twom24, twom12

  !                            default
  !  Luxury Level     0   1   2  *3*    4
  !    ndskip        /0, 24, 73, 199, 365/
  ! Corresponds to p=24  48  97  223  389
  !     time factor   1   2   3    6   10   on slow workstation
  !                   1 1.5   2    3    5   on fast mainframe

  ! public notyet, i24, j24, carry, seeds, twom24, twom12, luxlev
  ! public nskip, ndskip, in24, next, kount, mkount, inseed

  public ranlux, rluxgo

contains

subroutine ranlux(rvec, lenv)
  implicit none
  integer, intent(in) :: lenv
  real*8, intent(out) :: rvec(lenv)

  ! Local variables
  integer                :: i, k, lp, ivec, isk
  real*8                 :: uni

  ! NOTYET is .TRUE. if no initialization has been performed yet.
  ! Default Initialization by Multiplicative Congruential
  if (notyet) then
    notyet = .false.
    jseed = jsdflt
    inseed = jseed

    write(msgstr,"('RANLUX default initialization: ',I0)") jseed
    call log_msg(msgstr)

    luxlev = lxdflt
    nskip = ndskip(luxlev)
    lp = nskip + 24
    in24 = 0
    kount = 0
    mkount = 0

    write(msgstr,"('RANLUX default luxury level = ',I0,' , p = ',I0)") luxlev, lp
    call log_msg(msgstr)

    twom24 = 1.
    do i = 1, 24
      twom24 = twom24 * 0.5d0
      k = jseed / 53668
      jseed = 40014 * (jseed-k*53668) - k * 12211
      if (jseed.LT.0) jseed = jseed + icons
      iseeds(i) = MOD(jseed,itwo24)
    end do
    twom12 = twom24 * 4096.d0
    do i = 1, 24
      seeds(i) = dble(iseeds(i)) * twom24
      next(i) = i - 1
    end do
    next(1) = 24
    i24 = 24
    j24 = 10
    carry = 0.
    if (seeds(24).EQ.0.) carry = twom24
  end if

  ! The Generator proper: "Subtract-with-borrow",
  ! as proposed by Marsaglia and Zaman,
  ! Florida State University, March, 1989

  do ivec = 1, lenv
    uni = seeds(j24) - seeds(i24) - carry
    if (uni.LT.0.) then
      uni = uni + 1.0d0
      carry = twom24
    else
      carry = 0.
    end if
    seeds(i24) = uni
    i24 = next(i24)
    j24 = next(j24)
    rvec(ivec) = uni
    !  small numbers (with less than 12 "significant" bits) are "padded".
    if (uni.lt.twom12) then
      rvec(ivec) = rvec(ivec) + twom24 * seeds(j24)
      ! and zero is forbidden in case someone takes a logarithm
      if (rvec(ivec).eq.0.) rvec(ivec) = twom24 * twom24
    end if
    ! Skipping to luxury.  As proposed by Martin Luscher.
    in24 = in24 + 1
    if (in24.eq.24) then
      in24 = 0
      kount = kount + nskip
      do isk = 1, nskip
        uni = seeds(j24) - seeds(i24) - carry
        if (uni.lt.0.) then
          uni = uni + 1.0d0
          carry = twom24
        else
          carry = 0.
        end if
        seeds(i24) = uni
        i24 = next(i24)
        j24 = next(j24)
      end do
    end if
  end do
  kount = kount + lenv
  if (kount.ge.igiga) then
    mkount = mkount + 1
    kount = kount - igiga
  end if
  return

end subroutine ranlux

!> Subroutine to initialize from one or three integers
subroutine rluxgo(lux, ins, k1, k2)
  implicit none
  integer, intent(in) :: lux, ins, k1, k2

  !  Local variables
  integer                :: ilx, i, iouter, isk, k, inner, izip, izip2
  real*8                 :: uni

  call log_ws()

  if (lux.lt.0) then
    luxlev = lxdflt
  else if (lux.le.maxlev) then
    luxlev = lux
  else if (lux.lt.24.or.lux.gt.2000) then
    luxlev = maxlev
    write(msgstr,"('RANLUX illegal luxury RLUXGO: ',I0)") lux
    call log_msg(msgstr)
  else
    luxlev = lux
    do ilx = 0, maxlev
      if (lux.eq.ndskip(ilx)+24) luxlev = ilx
    end do
  end if

  if (luxlev.le.maxlev) then
    nskip = ndskip(luxlev)
    write(msgstr,"('RANLUX luxury level set by rluxgo to ',I0,', P = ',I0)") luxlev, nskip + 24
    call log_msg(msgstr)
  else
    nskip = luxlev - 24
    write(msgstr,"('RANLUX p-value set by RLUXGO to: ',I0)") luxlev
    call log_msg(msgstr)
  end if

  in24 = 0
  if (ins.lt.0) then
    call log_msg("RANLUX illegal initialization by RLUXGO, negative input seed")
  end if

  if (ins.gt.0) then
    jseed = ins
    write(msgstr,"('RANLUX initialized by rluxgo from seeds: ',I0,', ',I0,', ',I0)") jseed, k1, k2
  else
    jseed = jsdflt
    call log_msg("RANLUX initialized by RLUXGO from default seed")
  end if

  inseed = jseed
  notyet = .false.
  twom24 = 1.
  do i = 1, 24
    twom24 = twom24 * 0.5d0
    k = jseed / 53668
    jseed = 40014 * (jseed-k*53668) - k * 12211
    if (jseed.lt.0) jseed = jseed + icons
    iseeds(i) = mod(jseed,itwo24)
  end do

  twom12 = twom24 * 4096.d0

  do i = 1, 24
    seeds(i) = dble(iseeds(i)) * twom24
    next(i) = i - 1
  end do

  next(1) = 24
  i24 = 24
  j24 = 10
  carry = 0.d0

  if (seeds(24).eq.0.) carry = twom24
  ! If restarting at a break point, skip K1 + IGIGA*K2
  ! Note that this is the number of numbers delivered to
  ! the user PLUS the number skipped (if luxury .GT. 0).
  kount = k1
  mkount = k2

  if (k1+k2.ne.0) then
    do iouter = 1, k2 + 1
      inner = igiga
      if (iouter.eq.k2+1) inner = k1
      do isk = 1, inner
        uni = seeds(j24) - seeds(i24) - carry
        if (uni.lt.0.) then
          uni = uni + 1.0
          carry = twom24
        else
          carry = 0.
        end if
        seeds(i24) = uni
        i24 = next(i24)
        j24 = next(j24)
      end do
    end do
    ! Get the right value of IN24 by direct calculation
    in24 = mod(kount,nskip+24)
    if (mkount.gt.0) then
      izip = mod(igiga, nskip+24)
      izip2 = mkount * izip + in24
      in24 = mod(izip2, nskip+24)
    end if
    ! Now IN24 had better be between zero and 23 inclusive
    if (in24.gt.23) then
      write(msgstr,"('RANLUX error in restarting with RLUXGO: the values ',I0,', ',I0,', ',I0, ' cannot occur at luxury level ',I0)") &
        ins, k1, k2, luxlev
      call log_msg(msgstr)
      in24 = 0
    end if
  end if
  return

end subroutine rluxgo


! =====================================================
! | NOTE: rluxin, rluxut and rluxat are never called! |
! =====================================================

! !           Subroutine to input and float integer seeds from previous run
! SUBROUTINE rluxin
! !     the following IF BLOCK added by Phillip Helbig, based on conversation
! !     with Fred James; an equivalent correction has been published by James.
! 
! IMPLICIT NONE
! 
! !     Local variables
! 
! INTEGER             :: i, isd
! 
! IF (notyet) THEN
!   WRITE (6,'(A)') ' Proper results ONLY with initialisation from 25 ',  &
!   'integers obtained with RLUXUT'
!   notyet = .false.
! END IF
! 
! twom24 = 1.d0
! DO i = 1, 24
!   next(i) = i - 1
!   twom24 = twom24 * 0.5d0
! END DO
! next(1) = 24
! twom12 = twom24 * 4096.d0
! WRITE (6,'(A)') ' FULL INITIALIZATION OF RANLUX WITH 25 INTEGERS:'
! WRITE (6,'(5X,5I12)') isdext
! DO i = 1, 24
!   seeds(i) = dble(isdext(i)) * twom24
! END DO
! carry = 0.d0
! IF (isdext(25).LT.0) carry = twom24
! isd = IABS(isdext(25))
! i24 = MOD(isd,100)
! isd = isd / 100
! j24 = MOD(isd,100)
! isd = isd / 100
! in24 = MOD(isd,100)
! isd = isd / 100
! luxlev = isd
! IF (luxlev.LE.maxlev) THEN
!   nskip = ndskip(luxlev)
!   WRITE (6,'(A,I2)') ' RANLUX LUXURY LEVEL SET BY RLUXIN TO: ', luxlev
! ELSE IF (luxlev.GE.24) THEN
!   nskip = luxlev - 24
!   WRITE (6,'(A,I5)') ' RANLUX P-VALUE SET BY RLUXIN TO:', luxlev
! ELSE
!   nskip = ndskip(maxlev)
!   WRITE (6,'(A,I5)') ' RANLUX ILLEGAL LUXURY RLUXIN: ', luxlev
!   luxlev = maxlev
! END IF
! inseed = -1
! RETURN
! 
! END SUBROUTINE rluxin


! !                    Subroutine to ouput seeds as integers
! SUBROUTINE rluxut
! 
! IMPLICIT NONE
! 
! !     Local variables
! 
! INTEGER             :: i
! 
! DO i = 1, 24
!   isdext(i) = INT(seeds(i)*twop12*twop12)
! END DO
! isdext(25) = i24 + 100 * j24 + 10000 * in24 + 1000000 * luxlev
! IF (carry.GT.0.) isdext(25) = -isdext(25)
! RETURN
! 
! END SUBROUTINE rluxut


! !                    Subroutine to output the "convenient" restart point
! SUBROUTINE rluxat(lout, inout, k1, k2)
! 
! IMPLICIT NONE
! 
! INTEGER, INTENT(OUT) :: lout, inout, k1, k2
! 
! lout = luxlev
! inout = inseed
! k1 = kount
! k2 = mkount
! RETURN
! 
! END SUBROUTINE rluxat

end module luxury_rng_module
