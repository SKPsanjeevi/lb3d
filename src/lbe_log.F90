#include "lbe.h"

!> Functions for logging output to stdout.

module lbe_log_module

  use lbe_globals_module, only: myrankc

  implicit none
  include 'mpif.h' ! This is for MPI_Abort only

  private

  public log_msg, log_msg_hdr, log_ws, log_msg_ws, error, msgstr, stderr, stdin, stdout

  integer, parameter :: stderr = 0
  integer, parameter :: stdin  = 5
  integer, parameter :: stdout = 6

  character(len=512) :: msgstr     !< Global string buffer for messages

contains

!> Writes a message msg with a timestamp to stdout.
!> If forall is true, every processor does this, otherwise only rank 0 will.
!> If the optional forall argument is not specified, only rank 0 will.
subroutine log_msg(msg, forall, error)
  implicit none

  character(len=*), intent(in)  :: msg
  logical, optional, intent(in) :: forall
  logical, optional, intent(in) :: error

  integer, dimension(8) :: datevalues
  character(len=24)     :: datestr

  integer :: stream

  call date_and_time(values = datevalues)

  write(datestr,"(I0.2,':',I0.2,':',I0.2)") datevalues(5), datevalues(6), datevalues(7)

  if ( present(error) ) then
    if (error) then
      stream = stderr
    else
      stream = stdout
    end if
  else
    stream = stdout
  end if

  if ( present(forall) ) then
    if (forall) then
      write(stdout, "(A, ' - (', I6.6, ') ', A)" )  trim(datestr), myrankc, trim(msg)
    else if (myrankc == 0) then
      write(stdout, "(A, ' - ', A)") trim(datestr), trim(msg)
    end if
  else
    if (myrankc == 0) then
      write(stdout, "(A, ' - ', A)") trim(datestr), trim(msg)
    end if
  end if

end subroutine log_msg

!> Wraps \c log_msg and writes whitespace and fancy header symbols around the 
!> message.
subroutine log_msg_hdr(msg, forall)
  implicit none

  character(len=*), intent(in)  :: msg
  logical, optional, intent(in) :: forall

  integer, parameter :: header_len = 80
  character(len=256) :: msg_with_header
  integer :: length, i

  length = len_trim(msg)
  length = header_len - length

  msg_with_header = "( " // trim(msg) // " )"

  do i = 1, length
    if ( mod(i,2) == 0 ) then
      msg_with_header = "=" // trim(msg_with_header)
    else
      msg_with_header = trim(msg_with_header) // "="
    end if
  end do

  if ( present(forall) ) then
    call log_msg_ws(msg_with_header, forall)
  else
    call log_msg_ws(msg_with_header)
  end if

end subroutine log_msg_hdr  
  
!> Wraps \c log_msg and writes whitespace around the message.
subroutine log_msg_ws(msg, forall)
  implicit none

  character(len=*), intent(in)  :: msg
  logical, optional, intent(in) :: forall

  if ( present(forall) ) then
    call log_msg("",  forall)
    call log_msg(msg, forall)
    call log_msg("",  forall)
  else
    call log_msg("")
    call log_msg(msg)
    call log_msg("")
  end if

end subroutine log_msg_ws

!> Wraps \c log_msg and writes an empty line.
subroutine log_ws(forall)
  implicit none

  logical, optional, intent(in) :: forall

  if ( present(forall) ) then
    call log_msg("", forall)
  else
    call log_msg("")
  end if

end subroutine log_ws

!> Wraps \c log_msg, writes the message as an error and aborts the program.
subroutine error(msg)
  implicit none

  character(len=*), intent(in) :: msg
  
  integer :: mpierror

  call log_msg("ERROR: " // trim(msg), forall = .true., error = .false.)
  call log_msg("ERROR: " // trim(msg), forall = .true., error = .true. )
  ! Calling Abend here would create a circular dependency, so use this nasty direct MPI call instead.
  call MPI_Abort(MPI_COMM_WORLD, -1, mpierror)
  stop

end subroutine error

end module lbe_log_module

