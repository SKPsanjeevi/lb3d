#include "lbe.h"

!> This module extends the LB3D timer module with some ELEC-specific functionality.
module lbe_elec_timer_module
#ifdef ELEC

  use lbe_elec_globals_module
  use lbe_elec_helper_module
  use lbe_parallel_module, only: stats_rk
  use lbe_timer_module
  use mpi

  implicit none

  private

  public elec_init_timers, elec_report_timers, elec_reset_timers
  public ti_elec_apply_fluxes, ti_elec_calc_E_fd, ti_elec_calc_fluxes_advection, &
       &ti_elec_calc_fluxes_diffusion, ti_elec_equilibration, ti_elec_init, &
       &ti_elec_io, ti_elec_p3m, ti_elec_p3m_setup, ti_elec_p3m_poisson, &
       &ti_elec_p3m_unpack, ti_elec_p3m_calc_E, &
       &ti_elec_sor, ti_elec_sor_allreduce, ti_elec_sor_phi_halo


  public register_timer, start_timer, stop_timer, timers ! From timer_module

  integer :: ti_elec_apply_fluxes, ti_elec_calc_E_fd, ti_elec_calc_fluxes_advection, &
       &ti_elec_calc_fluxes_diffusion, ti_elec_equilibration, ti_elec_init, &
       &ti_elec_io, ti_elec_p3m, ti_elec_p3m_setup, ti_elec_p3m_poisson, &
       &ti_elec_p3m_unpack, ti_elec_p3m_calc_E, &
       &ti_elec_sor, ti_elec_sor_allreduce, ti_elec_sor_phi_halo

  character(*), parameter :: elec_timer_prefix = "ELEC:"

contains

!> Register all ELEC timers.
subroutine elec_init_timers()
  implicit none

  call log_msg_elec_hdr("Registering ELEC timers")

  call elec_register_timer("ApplyFluxes"       , ti_elec_apply_fluxes)
  call elec_register_timer("CalcE_FD"          , ti_elec_calc_E_fd)
  call elec_register_timer("CalcFluxesAdv"     , ti_elec_calc_fluxes_advection)
  call elec_register_timer("CalcFluxesDiff"    , ti_elec_calc_fluxes_diffusion)
  call elec_register_timer("Equilibration"     , ti_elec_equilibration)
  call elec_register_timer("Init"              , ti_elec_init)
  call elec_register_timer("IO"                , ti_elec_io)
  call elec_register_timer("P3M"               , ti_elec_p3m)
  call elec_register_timer("P3Msetup"          , ti_elec_p3m_setup)
  call elec_register_timer("P3MPoisson"        , ti_elec_p3m_poisson)
  call elec_register_timer("P3Munpack"         , ti_elec_p3m_unpack)
  call elec_register_timer("P3MCalcE"          , ti_elec_p3m_calc_E)
  call elec_register_timer("Sor"               , ti_elec_sor)
  call elec_register_timer("SorAllReduce"      , ti_elec_sor_allreduce)
  call elec_register_timer("SorPhiHalo"        , ti_elec_sor_phi_halo)

  call log_msg_elec_ws("Finished registering ELEC timers.")

end subroutine elec_init_timers

!> Register a single ELEC timer, taking care of the prefix.
subroutine elec_register_timer(tname, ti)
  implicit none

  character(len=*), intent(in) :: tname !> Name of the timer (without prefix)
  integer, intent(out) :: ti !> Timer code

  write(msgstr,"('Registering timer ', A, A)") elec_timer_prefix, trim(tname)
  call log_msg_elec(msgstr)
  call register_timer(elec_timer_prefix // trim(tname), ti)
  
end subroutine elec_register_timer

!> Writes a summary of all timer data - to stdout if no argument is passed, to unit = iounit_in otherwise.
subroutine elec_report_timers(iounit_in)
  implicit none

  integer, optional, intent(in) :: iounit_in

  integer, parameter :: nhist = 10 ! Number of bins in the histogram.

  real(kind=rk), allocatable :: avgs(:), avgsl(:)
  real(kind=rk) :: tavg, tmax, tmin
  real(kind=rk) :: total, ttavg, ttmax, ttmin
  integer :: i, num_timers
  integer :: mpierror, ierror
  integer :: ihisto(nhist), ihistotmp(nhist)

  integer :: iounit

  if ( present(iounit_in) ) then
    iounit = iounit_in
  else
    iounit = stdout
  end if

  num_timers = size(timers)

  allocate(avgsl(num_timers), avgs(num_timers), stat = ierror)
  call check_allocate(ierror, 'Unable to allocate buffer for ELEC timers.')

  avgsl(:) = timers(:)%total / nprocs

  call MPI_Reduce(avgsl, avgs, num_timers, LBE_REAL, MPI_SUM, 0, comm_cart, mpierror)

  call log_msg_elec_hdr("Reporting ELEC timers (with histograms)")

  ! Write header
  if ( myrankc == 0 ) then
    write(iounit,"()")
    write(iounit,"('#')")
    write(iounit,"(A2,A30,X,A12,X,5(A13,X),10(A6,X))") "#?", "Name", "Count","TotRank0","AvgRank0", "TotCPUAvg", "TotCPUMax", "TotCPUMin",&
         "H01", "H02", "H03", "H04","H05","H06","H07","H08","H09","H10"
    write(iounit,"('#')")
  end if

  total = 0.0_rk
  ttavg = 0.0_rk
  ttmax = 0.0_rk
  ttmin = 0.0_rk

  ! Calculate and write timer data
  do i = 1, num_timers
    if ( index(trim(timers(i)%name), elec_timer_prefix) == 1 ) then
      call stats_rk(timers(i)%total, tavg, tmax, tmin, ihisto, ihistotmp, nhist)
      if ( myrankc == 0 ) then
        write(iounit,"(A32,X,I12,X,5(ES13.6,X),10(I6,:,X))") &
             trim(timers(i)%name(len(elec_timer_prefix)+1:)), timers(i)%count, &
             &timers(i)%total, timers(i)%total / real(timers(i)%count, kind=rk), tavg, tmax, tmin, ihisto
        ! Sum only top-level timers.
        if ( trim(timers(i)%name(len(elec_timer_prefix)+1:)) == "ApplyFluxes"    .or. &
             trim(timers(i)%name(len(elec_timer_prefix)+1:)) == "CalcE_FD"       .or. &
             trim(timers(i)%name(len(elec_timer_prefix)+1:)) == "CalcFluxesAdv"  .or. &
             trim(timers(i)%name(len(elec_timer_prefix)+1:)) == "CalcFluxesDiff" .or. &
             trim(timers(i)%name(len(elec_timer_prefix)+1:)) == "Init"           .or. &
             trim(timers(i)%name(len(elec_timer_prefix)+1:)) == "IO"             .or. &
             trim(timers(i)%name(len(elec_timer_prefix)+1:)) == "P3M"            .or. &
             trim(timers(i)%name(len(elec_timer_prefix)+1:)) == "Sor" ) then 
          total = total + timers(i)%total
          ttavg = ttavg + tavg
          ttmax = ttmax + tmax
          ttmin = ttmin + tmin
        end if
      end if
    end if
  end do

  ! Write sum totals
  if ( myrankc == 0 ) then
    write(iounit,"()")
    write(iounit,"(A32,X,12X,X,ES13.6,X,14X,3(ES13.6,X))") &
         "Total", total, ttavg, ttmax, ttmin
  end if

  ! Write footer
  if ( myrankc == 0 ) then
    write(iounit,"()")
  end if

  deallocate(avgsl, avgs, stat = ierror)
  call check_allocate(ierror, "Unable to deallocate buffer for ELEC timers.")

end subroutine elec_report_timers

!> Reset all ELEC timers to zero total time and zero count.
subroutine elec_reset_timers()
  implicit none

  integer :: i, num_timers

  num_timers = size(timers)

  call log_msg_elec("Resetting ELEC timers.")

  do i = 1, num_timers
    ! Reset ONLY ELEC timers
    if ( index(trim(timers(i)%name), elec_timer_prefix) == 1 ) then
      timers(i)%start = 0.0_rk
      timers(i)%total = 0.0_rk
      timers(i)%count = 0
    end if
  end do
end subroutine elec_reset_timers

#endif
end module lbe_elec_timer_module
