#include "lbe.h"

!> general helper routines
module lbe_helper_module
  use lbe_globals_module
  use lbe_log_module
#ifdef MD
  use lbe_md_globals_module, only: interaction
#endif
  use lbe_parms_module
  use lbe_types_module

#ifndef NOIEEEARITHMETIC
  use ieee_arithmetic
#endif

  implicit none
  private
  include 'mpif.h'

  public check_dump_now,check_NaN,cross_product,density,every_n_time_steps&
       &,gaussianBoxMuller,iif,interpolate,lbe_count_sites,local_coordinates&
       &,makedatestr,massflow,norm,present_and_non_negative,present_and_true&
       &,report_check_NaN_function,unit_vector,unixtime,velocity,rotate_coordinates&
       &,is_restoring,is_wall,is_colloid,is_fluid,is_rock,is_surface&
       &,diag_symm_matrix,set_halo_extent,request_halo_extent, CS_EOS&
       &,quadratic_solver

  !> if-then-else-function, acts similarly the ternary operator in C
  interface iif
     module procedure iif_INTEGER
  end interface

  !> calculates local fluid densities
  interface density
     module procedure density_REAL_NVECS,density_LBE_SITE
  end interface

  !> interpolates scalar or vector fields
  interface interpolate
     module procedure interpolate_r,interpolate_r3
  end interface

  !> calculates local mass flow rates
  interface massflow
     module procedure massflow_REAL_NVECS,massflow_LBE_SITE
  end interface

  !> calculates local velocity vectors
  interface velocity
     module procedure velocity_REAL_NVECS,velocity_LBE_SITE
  end interface

  !> total number of LB lattice sites in the system
  integer(kind=8),save,public :: n_sites

  !> total number of fluid sites in the system
  integer(kind=8),save,public :: n_sites_fluid

  !> total number of particle-occupied sites in the system (\c
  !> interaction=='ladd' )
  integer(kind=8),save,public :: n_sites_particle

  !> total number of solid rock sites in the system (walls)
  integer(kind=8),save,public :: n_sites_rock

  !> total number of rock surface sites in the system
  integer(kind=8),save,public :: n_sites_surface

contains

!> Returns the cross (vector) product of two three-dimensional vectors
!>
!> \param[in] a vector \c a
!>
!> \param[in] b vector \c b
!>
!> \returns \f[\mathbf{a}\times\mathbf{b}\f]
pure function cross_product(a,b)
  real(kind=8),dimension(3) :: cross_product
  real(kind=8),intent(in) :: a(3),b(3)

  cross_product(1) = a(2)*b(3)-a(3)*b(2)
  cross_product(2) = a(3)*b(1)-a(1)*b(3)
  cross_product(3) = a(1)*b(2)-a(2)*b(1)
end function cross_product

!> function for calculating upper part of Carnahan-Starling EOS
!> a is density, b is g parameter in Shan-Chen force, Tcs is critical
!>temperature

pure function CS_EOS(a, b)
real(kind=8) :: CS_EOS
real(kind=8), intent(in) :: a, b

if (b==0) then
CS_EOS = 0.d0
else
CS_EOS = sqrt(2*((a*Tcs*(1+a+a**2-a**3)/(1-a)**3)-4/3*a**2)/(6*b))
end if
end function CS_EOS


!> uses the Box Mueller transform to generate to Gaussian random numbers with
!> mean value 0 and variance 1
subroutine gaussianBoxMuller(gaussian1,gaussian2)
    real(kind=rk) :: gaussian1, gaussian2
    real(kind=rk) :: u1,u2

    call random_number(u1)
    call random_number(u2)

    gaussian1 = sqrt(-2.0_rk*log(u1)) * cos(2.0_rk*pi*u2)
    gaussian2 = sqrt(-2.0_rk*log(u1)) * sin(2.0_rk*pi*u2)
end subroutine gaussianBoxMuller

!> if-then-else function, acts similarly as the ternary operator in C
!>
!> \param[in] cond condition to test
!>
!> \param[in] true return value if \c cond is true
!>
!> \param[in] false return value if \c cond is false
!>
!> \returns \c true if \c cond is \c .true., otherwise \c false
pure elemental function iif_INTEGER(cond,true,false)
    integer :: iif_INTEGER
    logical,intent(in) :: cond
    integer,intent(in) :: true,false

    if (cond) then
       iif_INTEGER = true
    else
       iif_INTEGER = false
    end if
end function iif_INTEGER

!> One should use this function to provide a request for enlarging the halo.
!> All requests are collected and then later executed by \c set_halo_extent .
!> Virtual fluids density, particle color, filled in the ourtermost layer of particle, [Fabian2011]. 
subroutine request_halo_extent(extent, name)
  implicit none

  integer, intent(in) :: extent
  character(len=*),intent(in) :: name
  integer :: stat, ti

  type(halo_request), allocatable, dimension(:) :: tmp
  
  if (.not.allocated(halo_requests)) then
    allocate (halo_requests(0),stat=stat)
    ! call check_allocate(stat,'request_halo_extent(): halo_requests')
  end if
  
  allocate (tmp(size(halo_requests)),stat=stat)
  ! call check_allocate(stat,'request_halo_extent(): tmp')

  tmp(:) = halo_requests(:)

  deallocate (halo_requests)
  allocate (halo_requests(size(tmp)+1),stat=stat)

  ! call check_allocate(stat,'request_halo_extent(): halo_requests')
  halo_requests(1:size(tmp)) = tmp(:)

  deallocate (tmp)
  
  ti = size(halo_requests)
  halo_requests(ti)%extent = extent
  halo_requests(ti)%name = name
end subroutine request_halo_extent

!> Report on all requested halo sizes and set \c halo_extent to the maximum value.
subroutine set_halo_extent()
  implicit none

  integer :: nr, he, i

  call log_msg_hdr("Setting halo_extent")

  write(msgstr,"('Initial halo_extent = ',I0,'.')") halo_extent
  call log_msg(msgstr)

  if ( allocated(halo_requests) ) then
    call log_msg("Requested halo extensions:")

    he = halo_extent

    nr = size(halo_requests)
    do i=1, nr
      write(msgstr,"('  ',A,': ',I0)") trim(halo_requests(i)%name), halo_requests(i)%extent
      call log_msg(msgstr)
      if ( halo_requests(i)%extent > he ) he = halo_requests(i)%extent
    end do

    halo_extent = he

  write(msgstr,"('Taking maximum value and setting halo_extent = ',I0,'.')") halo_extent
  call log_msg(msgstr)

  else
    call log_msg("No halo modification requested.")
  end if

end subroutine set_halo_extent

!> trilinear interpolation of a discrete scalar field
!>
!> \param[in] Nr floating point scalar field (position indices are
!> local as in N(x,y,z); there is a halo of depth 1)
!>
!> \param[in] x global position. Must be at some place that is covered
!> by the local \c N, however, this subroutine is smart enough to deal
!> with pbc.
!>
!> \param[out] r interpolated value
subroutine interpolate_r(Nr,x,r)
  real(kind=8),intent(in) :: Nr(0:,0:,0:)
  real(kind=8),intent(in) :: x(3)
  real(kind=8),intent(out) :: r
  integer l
  integer lx(3),lp(3),opp_lp(3) ! lattice point coordinates
  real(kind=8) :: xx(3),weight  

  CALL local_coordinates(x(:),xx(:))
  lx(:) = floor(xx(:))

  r = 0.0_8
  lattice_points: do l=1,n_lp
    ! weight each lattice point with the volume of the hypercube
    ! spanned by the particle position and the opposite lattice
    ! point
    lp(:) = lx(:) + lp_sur(:,l)
    opp_lp(:) = lx(:) + opp_lp_sur(:,l)
    weight = abs(product(real(opp_lp(:),kind=8)-xx(:)))
    r = r + weight*Nr(lp(1),lp(2),lp(3))
  end do lattice_points
end subroutine interpolate_r

!> trilinear interpolation of a discrete 3-vector field
!>
!> \param[in] Nr3 floating point 3-vector field (position indices are
!> local as in N(x,y,z); there is a halo of depth 1)
!>
!> \param[in] x global position. Must be at some place that is covered
!> by the local \c N, however, this subroutine is smart enough to deal
!> with pbc.
!>
!> \param[out] r3 interpolated value
subroutine interpolate_r3(Nr3,x,r3)
  real(kind=8),intent(in) :: Nr3(0:,0:,0:,1:)
  real(kind=8),intent(in) :: x(3)
  real(kind=8),intent(out) :: r3(3)
  integer l
  integer lx(3),lp(3),opp_lp(3) ! lattice point coordinates
  real(kind=8) :: xx(3),weight

  CALL local_coordinates(x(:),xx(:))
  lx(:) = floor(xx(:))

  r3(:) = 0.0_8
  lattice_points: do l=1,n_lp
    ! weight each lattice point with the volume of the hypercube
    ! spanned by the particle position and the opposite lattice
    ! point
    lp(:) = lx(:) + lp_sur(:,l)
    opp_lp(:) = lx(:) + opp_lp_sur(:,l)
    weight = abs(product(real(opp_lp(:),kind=8)-xx(:)))
    r3 = r3 + weight*Nr3(lp(1),lp(2),lp(3),:)
  end do lattice_points
end subroutine interpolate_r3

!> stores the frequencies of different site occupations into global
!> variables
!>
!> \param[in] N local lattice chunk (expects halo of extent 1 though
!> the halo is not needed)
!>
!> All sites are counted and the results are stored in \c n_sites, \c
!> n_sites_fluid, \c n_sites_particle, and \c n_sites_rock on all
!> tasks.
subroutine lbe_count_sites(N)
    type(lbe_site),dimension(0:,0:,0:),intent(in) :: N
    real(kind=rk) :: rbuf(5),sbuf(5)
    integer ierror,x,y,z
    integer(kind=8) :: ns,nfs,nps,nrs,nss

    ns = 0
    nfs = 0
    nps = 0
    nrs = 0
    nss = 0
    do x=1,nx
       do y=1,ny
          do z=1,nz
             ns = ns+1
             if (N(x,y,z)%rock_state==0.0_rk) then
                nfs = nfs+1
#ifdef MD
             else if (interaction=='ladd'.and.N(x,y,z)%rock_state>0.0_rk) then
                nps = nps+1
#endif
             else
                nrs = nrs+1
                if ( is_surface( N(x,y,z)%rock_state ) ) nss = nss+1
             end if
          end do
       end do
    end do

    sbuf = real((/ns,nfs,nps,nrs,nss/),kind=rk)
    call MPI_Allreduce(sbuf,rbuf,5,LBE_REAL,MPI_SUM,MPI_COMM_WORLD,ierror)
    n_sites = nint(rbuf(1),kind=8)
    n_sites_fluid = nint(rbuf(2),kind=8)
    n_sites_particle = nint(rbuf(3),kind=8)
    n_sites_rock = nint(rbuf(4),kind=8)
    n_sites_surface = nint(rbuf(5),kind=8)
end subroutine lbe_count_sites

!> convert global to local coordinates
!>
!> \param[in] global global position vector
!>
!> \param[out] local local position vector
!>
!> \warning Thereby check for pbc, nevertheless result will be
!> worthless if \c global has no local representation but is in the
!> domain of another process (even if it is still in the local halo!).
subroutine local_coordinates(global,local)
  real(kind=8),intent(in) :: global(3)
  real(kind=8),intent(out) :: local(3)
  integer k

  ! If  global  is wrapped around due to pbc, wrap it back to get
  ! into the valid index range of  N(:,:,:) . Notice that this does not
  ! help if  global  is out of range for some other reason/bug!
  do k=1,3
    ! This is very tricky. Why was there halo_extent and -1.5 ? These
    ! values are surely wrong since they lead to a wrap-around already
    ! for global(1)==tnx if halo_extent==1 . This problem does not
    ! occur for the new code but I am not completely sure that this
    ! change will not break something else, maybe in connection with
    ! the stuff implemented by David. The reason for the 0.5 is that
    ! it is likely that particles are that much outside border because
    ! exchange() is called only on demand (see list_still_valid() !).
    ! FJ, 2009-09-03
!!$           if (global(k)>=border(2,k)+(real(halo_extent,kind=8)-1.5_8)) &
    if (global(k)>=border(2,k)+0.5_rk) then
      local(k) = global(k) - tsize(k)
!!$           else if (global(k)<border(1,k)-(real(halo_extent,kind=8)-1.5_8)) &
    else if (global(k)<border(1,k)-0.5_rk) then
      local(k) = global(k) + tsize(k)
    else
      local(k) = global(k)
    end if
  end do
  ! convert  local(:)  to local lbe coordinates
  local = local - (border(1,:)-0.5_rk)
end subroutine local_coordinates

!> Returns the length of a vector
!>
!> \param[in] v vector (3-dimensional)
!>
!> \returns \f[|\mathbf{v}|\f]
pure function norm(v)
  real(kind=8) :: norm
  real(kind=8),intent(in) :: v(3)
  norm = sqrt(dot_product(v,v))
end function norm

!> Normalizes a vector to unit length
!>
!> \param[in] v vector (3-dimensional)
!>
!> \returns \f[\mathbf{v}/|\mathbf{v}|\f]
pure function unit_vector(v)
  real(kind=8),dimension(3) :: unit_vector
  real(kind=8),intent(in) :: v(3)
  unit_vector = v/sqrt(dot_product(v,v))
end function unit_vector

!> Finds solution for quadratic equation 
!>
!> \param[in] 3 coefficients of quadratic equation
function quadratic_solver(a,b,c,cond)
  logical :: cond
  real(kind=8),intent(in) :: a,b,c
  real(kind=8)            :: quadratic_solver(2),disc
  quadratic_solver(:) = 0.0_rk
  disc = b*b-4.0_rk*a*c
  if(disc.ge.0.0_rk) then
    disc =  sqrt(disc)
    quadratic_solver(1) = (-b+disc)/(2.0_rk*a)
    quadratic_solver(2) = (-b-disc)/(2.0_rk*a)
    cond = .true.
  else
    cond = .false.
!    write(msgstr,"('No real roots. Discriminant = ',F16.10)") disc
  end if

end function quadratic_solver


!> Calculates the density for a given set of distributions.
!>
!> \param[in] n local population vector
!>
!> \returns density (not scaled with \c amass_{rbs})
function density_REAL_NVECS(n)
  real(kind=rk) :: density_REAL_NVECS
  real(kind=rk),intent(in) :: n(nvecs)

  density_REAL_NVECS = 2.0_rk*sum(n(1:6))+sum(n(7:19))
end function density_REAL_NVECS

!> Calculates the total density for a given lattice site.
!>
!> \param[in] s lattice site
!>
!> \returns density (not scaled with \c amass_{rbs})
function density_LBE_SITE(s)
  real(kind=rk) :: density_LBE_SITE
  type(lbe_site),intent(in) :: s

  density_LBE_SITE = &
#ifndef SINGLEFLUID
       &amass_b*density(s%n_b)+&
#ifndef NOSURFACTANT
       &amass_s*density(s%n_s)+&
#endif
#endif
       &amass_r*density(s%n_r)
end function density_LBE_SITE

!> Calculates the mass flow (that is the momentum) for a given set of
!> distributions, a possible rescaling of masses using \c amass_[rbs]
!> is not taken into account.
!>
!> \param[in] n local population vector
!>
!> \returns mass flow vector
function massflow_REAL_NVECS(n)
  real(kind=rk) :: massflow_REAL_NVECS(3)
  real(kind=rk),intent(in) :: n(nvecs)

  massflow_REAL_NVECS(1) = &
       &2.0_rk*(n(1)-n(2))+n(7)+n(8)+n(9)+n(10)-n(11)-n(12)-n(13)-n(14)
  massflow_REAL_NVECS(2) = &
       &2.0_rk*(n(3)-n(4))+n(7)+n(11)+n(15)+n(16)-n(8)-n(12)-n(17)-n(18)
  massflow_REAL_NVECS(3) = &
       &2.0_rk*(n(5)-n(6))+n(9)+n(13)+n(15)+n(17)-n(10)-n(14)-n(16)-n(18)
end function massflow_REAL_NVECS

!> Calculates the total mass flow (that is the momentum) for a given
!> lattice site taking into account \c amass_[rbs].
!>
!> \param[in] s lattice site
!>
!> \returns mass flow vector
function massflow_LBE_SITE(s)
  real(kind=rk) :: massflow_LBE_SITE(3)
  type(lbe_site),intent(in) :: s


!massflow_LBE_SITE = &
!#ifdef COMMON_VEL_FIX
!#ifndef OLD_VEL
!#ifndef SINGLEFLUID
!       &amass_b*massflow(s%n_b_pre)+&
!#ifndef NOSURFACTANT
!       &amass_s*massflow(s%n_s_pre)+&
!#endif
!#endif
!       &amass_r*massflow(s%n_r_pre)
!#else 
!#ifndef SINGLEFLUID
!       &amass_b*massflow(s%n_b)+&
!#ifndef NOSURFACTANT
!       &amass_s*massflow(s%n_s)+&
!#endif
!#endif
!       &amass_r*massflow(s%n_r)
!#endif

  massflow_LBE_SITE = &
#ifndef SINGLEFLUID
       &amass_b*massflow(s%n_b)+&
#ifndef NOSURFACTANT
       &amass_s*massflow(s%n_s)+&
#endif
#endif
       &amass_r*massflow(s%n_r)
end function massflow_LBE_SITE

!> Calculates the velocity for a given lattice site.
!>
!> \param[in] s lattice site
!>
!> \returns velocity vector
function velocity_LBE_SITE(s)
  real(kind=rk) :: velocity_LBE_SITE(3)
  type(lbe_site),intent(in) :: s

  velocity_LBE_SITE = massflow_LBE_SITE(s)/density_LBE_SITE(s)
end function velocity_LBE_SITE

!> Calculates the velocity for a given set of distributions.
!>
!> \param[in] n local population vector
!>
!> \returns velocity vector
function velocity_REAL_NVECS(n)
  real(kind=rk) :: velocity_REAL_NVECS(3)
  real(kind=rk),intent(in) :: n(nvecs)

  velocity_REAL_NVECS = massflow_REAL_NVECS(n)/density_REAL_NVECS(n)
end function velocity_REAL_NVECS

subroutine makedatestr(datestr)
  integer,dimension(8) :: datevalues
  character(len=24)    :: datestr
  CALL date_and_time(values=datevalues)
  write(datestr,"(I0.4,'-',I0.2,'-',I0.2,X,I0.2,':',I0.2,':',I0.2)") &
    datevalues(1), datevalues(2), datevalues(3), datevalues(5), datevalues(6), datevalues(7)
end subroutine makedatestr

!> Some code to calculate UNIX time_ts
elemental integer function is_leap_year(year)
  integer, intent(in) :: year

  ! http://www.timeanddate.com/date/leapyear.html says:
  !
  ! In the Gregorian calendar, which is the calendar used by
  ! most modern countries, the following rules decides which
  ! years are leap years:
  !
  ! 1. Every year divisible by 4 is a leap year.
  !
  ! 2. But every year divisible by 100 is NOT a leap year
  !
  ! 3. Unless the year is also divisible by 400, then it is
  ! still a leap year
  !
  ! The original normative standard is the papal bull
  ! "Inter Gravissimas", of 1582.
  !

  if (modulo(year,400) == 0) then
    is_leap_year = 1
  else
    if (modulo(year,100) == 0) then
      is_leap_year = 0
    else
      if (modulo(year,4) == 0) then
        is_leap_year = 1
      else
        is_leap_year = 0
      endif
    endif
  endif

end function is_leap_year

elemental integer function days_in_year(year)
  integer, intent(in) :: year
  if (1 == is_leap_year(year)) then
    days_in_year = 366
  else
    days_in_year = 365
  end if
end function days_in_year

!> Return the number of days between January 1 1970 and January 1 of
!> given year
elemental integer function ydays_since_epoch(year)
  integer, intent(in) :: year
  integer :: i
  integer :: total

  total = 0
  do i=1970,(year-1)
    total = total + days_in_year(i)
  end do
  ydays_since_epoch = total
end function ydays_since_epoch

!> Return the time_t, i.e. the number of seconds elapsed
!> since 00:00:00 UTC, January 1, 1970.
!>
!> Does not correct for leap seconds, consistent with POSIX.1
function unixtime()
  implicit none
  integer :: unixtime
  character(len=10) :: datestr,timestr,zonestr
  integer, dimension(8) :: values
  integer, dimension(12) :: monthlengths = &
          (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
  integer :: days_since_epoch

  integer :: i
  integer :: delta_seconds ! Timezone offset

  call date_and_time(datestr,timestr,zonestr,values)
  ! values returns with:
  ! values(1) = year
  ! values(2) = month of year
  ! values(3) = day of month
  ! values(4) = time different in minutes wrt utc
  ! values(5) = hour of day
  ! values(6) = minutes of hour
  ! values(7) = seconds of minute
  ! values(8) = milliseconds of second

  ! Figure out if this is a leap year, and
  ! adjust accordingly.

  if (1 == is_leap_year(values(1))) then
    monthlengths(2) = 29
  end if

  ! Days between Jan 1 this year and Jan 1 1970
  days_since_epoch = ydays_since_epoch(values(1))

  ! Now add up days between start of this month and Jan 1
  do i=1,(values(2)-1)
    days_since_epoch = days_since_epoch + monthlengths(i)
  end do
  days_since_epoch = days_since_epoch + values(3)-1

  ! days_since_epoch is now number of complete days since Jan 1 1970
  ! now correct for timezone

  delta_seconds = values(4)*60

  ! Note that the Intel compiler *will* adjust the time when you
  ! set the TZ environment variable, but it won't set values(4),
  ! so you can't tell whether you're in UTC or not.
  !
  !  FORTRAN! Because your blood pressure is too low, and there
  !  aren't enough head-shaped dents in the wall!
  !

  unixtime = values(7)                    &
            + 60*( values(6)               &
            + 60*( values(5)               &
            + 24* days_since_epoch ))      &
            - delta_seconds

end function unixtime

!> This function makes sure that the boolean is checked before the mod.
!>
!> \param[in] sci disables or enables check of \c n_sci
!>
!> \param[in] n_sci if checked, returns \c .true. every \c n_sci time steps
!>
!> \returns \c .false. if \c .not.sci, \c every_n_time_steps(n_sci) otherwise
!>
!> 2010-05-06. Added by Stefan. Putting this in a function avoids
!> nasty crashes, and keeps the nested if-loops confined to one place.
logical function check_dump_now(sci, n_sci)
  logical,intent(in) :: sci
  integer,intent(in) :: n_sci

  if ( sci ) then
     check_dump_now = every_n_time_steps(n_sci)
  else
     check_dump_now = .false.
  endif
end function check_dump_now

!> checks whether the n-ths time step is reached while safely
!> returning \c .false. if \c n==0
!>
!> \param[in] n interval against which to check \c nt
!>
!> \returns \c .false. if \c n==0, \c mod(nt.n)==0 otherwise
!>
!> This can be used in querying md-style time interval parameters
!> where a value of 0 is valid and means that the respective feature
!> is disabled.
logical function every_n_time_steps(n)
  integer,intent(in) :: n

  if (n==0) then
    every_n_time_steps = .false.
  else
    every_n_time_steps = mod(nt,n)==0
  end if
end function every_n_time_steps

!> safely checks whether an optional real argument is present and
!> larger or equal to \c 0.0_rk
!>
!> \param[in] arg optional \c real(kind=rk) argument which might be
!> not present
!>
!> \returns \c .false. if \c arg is missing or \c arg<0.0_rk ,
!> otherwise \c .true.
elemental logical function present_and_non_negative(arg)
  real(kind=rk),intent(in),optional :: arg

  if (present(arg)) then
    if (arg>=0.0_rk) then
      present_and_non_negative = .true.
    else
      present_and_non_negative = .false.
    end if
  else
    present_and_non_negative = .false.
  end if
end function present_and_non_negative

!> safely checks whether an optional logical argument is present and
!> \c .true.
!>
!> \param[in] arg optional logical argument which might be not present
!>
!> \returns \c .false. if \c arg is missing or \c .not.arg ,
!> otherwise \c .true.
elemental logical function present_and_true(arg)
  logical,intent(in),optional :: arg

  if (present(arg)) then
    if (arg) then
      present_and_true = .true.
    else
      present_and_true = .false.
    end if
  else
    present_and_true = .false.
  end if
end function present_and_true

elemental logical function check_NaN(n)
  implicit none
  real(kind=rk),intent(in) :: n
#ifdef NOIEEEARITHMETIC
#ifndef NOISNAN
  check_NaN = isnan(n)
#else
  ! UNABLE TO CHECK FOR NANs!
  check_NaN = .false.
#endif
#else
  check_NaN = ieee_is_nan(n)
#endif
  return
end function check_NaN

elemental subroutine report_check_NaN_function(str)
  implicit none
  character(len=*), intent(out) :: str
#ifdef NOIEEEARITHMETIC
#ifndef NOISNAN
  str = "NOIEEEARITHMETIC is set, NaNs are checked through function 'isnan'."
#else
  str = "WARINING: NOIEEEARITHMETIC is set and NOISNAN is set, NaNs are not checked."
#endif
#else
  str = "NOIEEEARITHMETIC is not set, NaNs are checked through function 'ieee_is_nan'."
#endif
end subroutine report_check_NaN_function

! One can restore in three ways:
!   Set init_cond = 7 in the input file. Variable restore_string from the input file will be used.
!   Set restore = .true. in the input file. Variable restore_string from the input file will be used.
!   Set command line argument -r <restore_string>
logical function is_restoring()
  implicit none
  is_restoring = ( ( init_cond .eq. 7 ) .or. restore .or. arg_restore_string_set )
  return
end function is_restoring

elemental subroutine rotate_coordinates(rotstr, ox, oy, oz, rx, ry, rz)
  implicit none
  integer, intent(in)  :: ox, oy, oz ! Original coordinates
  integer, intent(out) :: rx, ry, rz ! Rotated coordinates
  character(len=*), intent(in) :: rotstr

  select case (rotstr)
    case ('xyz')
      rx = ox
      ry = oy
      rz = oz
    case ('xzy')
      rx = ox
      ry = oz
      rz = oy
    case ('yxz')
      rx = oy
      ry = ox
      rz = oz
    case ('yzx')
      rx = oy
      ry = oz
      rz = ox
    case ('zxy')
      rx = oz
      ry = ox
      rz = oy
    case ('zyx')
      rx = oz
      ry = oy
      rz = ox
    case default
      !call log_msg("WARNING: Unknown rotation string '"//rotstr//"', assuming 'xyz'.")
      rx = ox
      ry = oy
      rz = oz
    end select
end subroutine rotate_coordinates

elemental logical function is_colloid(rs)
  implicit none
  real(kind=rk),intent(in) :: rs
#ifdef MD
  if ( rs > 0.0_rk ) then
    is_colloid = .true.
  else
    is_colloid = .false.
  endif
#else
  is_colloid = .false.
#endif
end function is_colloid

elemental logical function is_wall(rs)
  implicit none
  real(kind=rk),intent(in) :: rs
#ifdef MD
  if ( rs < 0.0_rk) then
    is_wall = .true.
  else
    is_wall = .false.
  endif
#else
  if (rs > 0.0_rk) then
    is_wall = .true.
  else
    is_wall = .false.
  endif
#endif
end function is_wall

elemental logical function is_rock(rs)
  implicit none
  real(kind=rk),intent(in) :: rs
  is_rock = ( is_colloid(rs) .or. is_wall(rs) )
end function is_rock

elemental logical function is_fluid(rs)
  implicit none
  real(kind=rk),intent(in) :: rs
  is_fluid = .not. ( is_colloid(rs) .or. is_wall(rs) )
end function is_fluid

elemental logical function is_surface(rs)
  implicit none
  real(kind=rk), intent(in) :: rs
  is_surface = rs .eq. rock_value_surface
end function is_surface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Diagonalize symmetric 3x3-matrix
  !>
  !> A symmetric 3x3-matrix is diagonalized, and the eigenvalues and normalized eigenvectors are stored.

! ----------------------------------------------------------------------------
      SUBROUTINE diag_symm_matrix(A, Q, W)
! ----------------------------------------------------------------------------
! Original name of the subroutine: DSYEVJ3 (renamed for convenience),
! free code taken from http://www.mpi-hd.mpg.de/personalhomes/globes/3x3/index.html
! If you use one of the algorithms to prepare a scientific publication or talk, please cite:
! Joachim Kopp
! Efficient numerical diagonalization of hermitian 3x3 matrices
! Int. J. Mod. Phys. C 19 (2008) 523-548
! arXiv.org: physics/0610206
!
! Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
! matrix A using the Jacobi algorithm.
! The upper triangular part of A is destroyed during the calculation,
! the diagonal elements are read but not destroyed, and the lower
! triangular elements are not referenced at all.
! ----------------------------------------------------------------------------
! Parameters:
!   A: The symmetric input matrix
!   Q: Storage buffer for eigenvectors
!   W: Storage buffer for eigenvalues
! ----------------------------------------------------------------------------
!     .. Arguments ..
      DOUBLE PRECISION A(3,3)
      DOUBLE PRECISION Q(3,3)
      DOUBLE PRECISION W(3)

!     .. Parameters ..
      INTEGER          N
      PARAMETER        ( N = 3 )

!     .. Local Variables ..
      DOUBLE PRECISION SD, SO
      DOUBLE PRECISION S, C, T
      DOUBLE PRECISION G, H, Z, THETA
      DOUBLE PRECISION THRESH
      INTEGER          I, X, Y, R

!     Initialize Q to the identitity matrix
!     --- This loop can be omitted if only the eigenvalues are desired ---
      DO 10 X = 1, N
        Q(X,X) = 1.0D0
        DO 11, Y = 1, X-1
          Q(X, Y) = 0.0D0
          Q(Y, X) = 0.0D0
   11   CONTINUE
   10 CONTINUE

!     Initialize W to diag(A)
      DO 20 X = 1, N
        W(X) = A(X, X)
   20 CONTINUE

!     Calculate SQR(tr(A))
      SD = 0.0D0
      DO 30 X = 1, N
        SD = SD + ABS(W(X))
   30 CONTINUE
      SD = SD**2

!     Main iteration loop
      DO 40 I = 1, 50
!       Test for convergence
        SO = 0.0D0
        DO 50 X = 1, N
          DO 51 Y = X+1, N
            SO = SO + ABS(A(X, Y))
   51     CONTINUE
   50   CONTINUE
        IF (SO .EQ. 0.0D0) THEN
          RETURN
        END IF

        IF (I .LT. 4) THEN
          THRESH = 0.2D0 * SO / N**2
        ELSE
          THRESH = 0.0D0
        END IF

!       Do sweep
        DO 60 X = 1, N
          DO 61 Y = X+1, N
            G = 100.0D0 * ( ABS(A(X, Y)) )
            IF ( I .GT. 4 .AND. ABS(W(X)) + G .EQ. ABS(W(X)) .AND. ABS(W(Y)) + G .EQ. ABS(W(Y)) ) THEN
              A(X, Y) = 0.0D0
            ELSE IF (ABS(A(X, Y)) .GT. THRESH) THEN
!             Calculate Jacobi transformation
              H = W(Y) - W(X)
              IF ( ABS(H) + G .EQ. ABS(H) ) THEN
                T = A(X, Y) / H
              ELSE
                THETA = 0.5D0 * H / A(X, Y)
                IF (THETA .LT. 0.0D0) THEN
                  T = -1.0D0 / (SQRT(1.0D0 + THETA**2) - THETA)
                ELSE
                  T = 1.0D0 / (SQRT(1.0D0 + THETA**2) + THETA)
                END IF
              END IF

              C = 1.0D0 / SQRT( 1.0D0 + T**2 )
              S = T * C
              Z = T * A(X, Y)

!             Apply Jacobi transformation
              A(X, Y) = 0.0D0
              W(X)    = W(X) - Z
              W(Y)    = W(Y) + Z
              DO 70 R = 1, X-1
                T       = A(R, X)
                A(R, X) = C * T - S * A(R, Y)
                A(R, Y) = S * T + C * A(R, Y)
   70         CONTINUE
              DO 80, R = X+1, Y-1
                T       = A(X, R)
                A(X, R) = C * T - S * A(R, Y)
                A(R, Y) = S * T + C * A(R, Y)
   80         CONTINUE
              DO 90, R = Y+1, N
                T       = A(X, R)
                A(X, R) = C * T - S * A(Y, R)
                A(Y, R) = S * T + C * A(Y, R)
   90         CONTINUE

!             Update eigenvectors
!             --- This loop can be omitted if only the eigenvalues are desired ---
              DO 100, R = 1, N
                T       = Q(R, X)
                Q(R, X) = C * T - S * Q(R, Y)
                Q(R, Y) = S * T + C * Q(R, Y)
  100         CONTINUE
            END IF
   61     CONTINUE
   60   CONTINUE
   40 CONTINUE

      END SUBROUTINE

end module lbe_helper_module
