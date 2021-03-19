#include "lbe.h"

!> easy to use timer module for code profiling

!> The first timer created has a special meaning: It should be started
!> and stopped in a way that the total runtime is detected, including
!> the time measured by all other timers. It is used to determine the
!> count for the timer "Other" that is automatically added at the end
!> of the list of timers and contains the time spent in all unassessed
!> parts of the code.

!> \par
!> Please note that the plain time counts are averaged over all
!> processes but the avg/min/max and histogram output might give a
!> rough idea of load balance.
module lbe_timer_module
    use lbe_globals_module, only: rk, myrankc
    use lbe_log_module
    use lbe_parallel_module, only: check_allocate,comm_cart,nprocs,stats_rk

    implicit none
    include 'mpif.h'
    private

    public register_timer,start_timer,stop_timer,sum_timer
    public timers,lbe_md_timer ! For ELEC

    !> specifies one timer
    type lbe_md_timer
       real(kind=rk) start      !< start time of the current counting period
       real(kind=rk) total      !< cumulative timer count
       integer :: count         !< number of times the timer is started
       character(len=32) name   !< short human-readable description
    end type lbe_md_timer

    !> contains all timers
    type(lbe_md_timer),save,allocatable,dimension(:) :: timers
contains

    !> Create a new timer with name \c name. Its id is returned in \c ti.
    subroutine register_timer(name,ti)
        character(len=*),intent(in) :: name
        integer,intent(out) :: ti
        integer stat
        type(lbe_md_timer),allocatable,dimension(:) :: tmp

        if (.not.allocated(timers)) then
           allocate (timers(0),stat=stat)
           call check_allocate(stat,'register_timer(): timers')
        end if

        allocate (tmp(size(timers)),stat=stat)
        call check_allocate(stat,'register_timer(): tmp')
        tmp(:) = timers(:)
        deallocate (timers)
        allocate (timers(size(tmp)+1),stat=stat)
        call check_allocate(stat,'register_timer(): timers')
        timers(1:size(tmp)) = tmp(:)
        deallocate (tmp)

        ti = size(timers)
        timers(ti)%total = 0.0
        timers(ti)%count = 0
        timers(ti)%name = name
    end subroutine register_timer

    !> let timer \c ti start counting
    subroutine start_timer(ti)
        integer,intent(in) :: ti

        timers(ti)%start = mpi_wtime()
        timers(ti)%count = timers(ti)%count + 1
    end subroutine start_timer

    !> let timer \c ti stop counting
    subroutine stop_timer(ti)
        integer,intent(in) :: ti

        timers(ti)%total = timers(ti)%total + mpi_wtime() - timers(ti)%start
    end subroutine stop_timer

    !> Writes a summary of all timer data to all units specified in  units .
    subroutine sum_timer(units)
        integer,intent(in) :: units(:)
        integer i,j,nt,ierror,stat
        integer ihisto(10),ihistotmp(10)
        real(kind=rk),allocatable :: avg(:),aavg(:)
        real(kind=rk) tavg,ave,xmax,xmin

        ! add another timer that will hold the differenc between the "Total"
        ! timer (always the first) and the sum of all other timers
        call register_timer('Other',nt)

        timers(nt)%total = timers(1)%total - sum(timers(2:nt-1)%total)

        allocate(avg(nt),aavg(nt),stat=stat)
        call check_allocate(stat,'sum_timer(): avg,aavg')

        aavg(:) = timers(:)%total/nprocs
        call mpi_reduce(aavg,avg,size(timers),LBE_REAL,MPI_SUM,0,comm_cart,&
             &ierror)

        tavg = avg(1)           ! average of "Total" timer

        call log_msg_hdr("Reporting timers (absolute / relative)")

        if (myrankc==0) then
           do j=1,size(units)
              write (unit=units(j),fmt='()')
           end do
           do i=1,nt
              do j=1,size(units)
                 write (unit=units(j),fmt='(A," time/%",F15.6,F13.4)') &
                      &timers(i)%name,avg(i),100*avg(i)/tavg
              end do
           end do
           do j=1,size(units)
              write (unit=units(j),fmt='()')
           end do
        end if

        call log_msg_hdr("Reporting timers (with histograms)")

        if (myrankc==0) then
          do j=1,size(units)
            write (unit=units(j),fmt='()')
          end do
        end if
        do i=1,nt
           call stats_rk(timers(i)%total,ave,xmax,xmin,ihisto,ihistotmp,10)
           if (myrankc==0) then
              do j=1,size(units)
                 write (unit=units(j),fmt=&
                      &'(A," time: ",F13.4," ave",F13.4," max",F13.4," min")')&
                      &timers(i)%name,ave,xmax,xmin
                 write (unit=units(j),fmt='("  Histogram: ",10(I6,:,X))') ihisto
              end do
           end if
        end do
        if (myrankc==0) then
          do j=1,size(units)
            write (unit=units(j),fmt='()')
          end do
        end if
    end subroutine sum_timer
end module lbe_timer_module
