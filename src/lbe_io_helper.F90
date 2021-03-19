#include "lbe.h"

!> helper routines for IO
module lbe_io_helper_module
  use lbe_log_module
  use lbe_parallel_module
  use lbe_parms_module, only: chk_uid,dump_double,dump_format,folder,cpfolder&
       &,gr_out_file,nt,restore_string,srccpfolder

  implicit none
  private

  public lbe_delete_file
  public lbe_make_filename_output, lbe_make_filename_output_substep
  public lbe_make_filename_cp_rank, lbe_make_filename_rank
  public lbe_make_filename_restore_rank, lbe_make_filename_cp, lbe_make_filename_restore, lbe_make_filename_append
  public dump_avs_fld
  public die_unless_exists, die_unless_nonzero
  public restore_fmt_p, restore_fmt_t, restore_fmt_t_read

  character(len=32), parameter :: restore_fmt_t = "('t',I8.8,'-',I10.10)"
  character(len=32), parameter :: restore_fmt_t_read = "(1X,I8.8,1X,I10.10)"	
  character(len=32), parameter :: restore_fmt_p = "('p',I8.8,'-',I10.10)"

contains

!> Wrapper for file deletion.
subroutine lbe_delete_file(filename)
  implicit none

  character(len=*)   :: filename
  integer, parameter :: file_unit = 10
  logical :: file_exists

  inquire(file=trim(filename), exist=file_exists)
  if (file_exists) then
    open(unit=file_unit, file=trim(filename), status="OLD")
    close(unit=file_unit, status="DELETE")
  end if

end subroutine lbe_delete_file

!> Makes a filename by appending the given suffix and
!> timestep information, to the stem formed from the prefix and the value
!> of \c gr_out_file.
!> 
!> \c buffer will return this filename.
!>  This would be trivial in a sane language.
subroutine lbe_make_filename_output(buffer, prefix, suffix, timestep)
  implicit none

  character(len=*), intent(out) :: buffer
  character(len=*), intent(in) :: prefix, suffix
  integer, intent(in) :: timestep
  character(len=20)  :: chkuidstr

  write(chkuidstr, FMT=restore_fmt_t) timestep, chk_uid

  write(buffer,"('./',A,'/',A,'_',A,'_',A,A)") &
    trim(folder), trim(prefix), trim(gr_out_file), trim(chkuidstr), trim(suffix)

end subroutine lbe_make_filename_output

!> Makes a filename by appending the given suffix and (sub-)timestep
!> information, to the stem formed from the prefix and the value of \c
!> gr_out_file.
!>
!> \param[out] buffer will return this filename.
!>
!> \param[in] prefix file name prefix, for example "md-cfg"
!>
!> \param[in] suffix file name extension, for example ".xdr"
!>
!> \param[in] time number of the current time step
!>
!> \param[in] substep number of the current substep
subroutine lbe_make_filename_output_substep(buffer,prefix,suffix,timestep,substep)
  implicit none

  character(len=*),intent(out) :: buffer
  character(len=*),intent(in) :: prefix,suffix
  integer,intent(in) :: timestep,substep
  character(len=20)  :: chkuidstr

  write(chkuidstr, FMT=restore_fmt_t) timestep, chk_uid

  write (buffer,"('./',A,'/',A,'_',A,'_',A,'_s',I6.6,A)") &
       &trim(folder),trim(prefix),trim(gr_out_file),trim(chkuidstr),substep&
       &,trim(suffix)
end subroutine lbe_make_filename_output_substep

!> Makes a filename by appending the given suffix, plus CPU rank and
!> timestep information, to the stem formed from the prefix and the value
!> of \c gr_out_file.
!> 
!> \c buffer will return this filename.
!>  This would be trivial in a sane language.
subroutine lbe_make_filename_cp_rank(buffer, prefix, suffix, timestep, rank)
  implicit none
  character(len=*), intent(out) :: buffer
  character(len=*), intent(in) :: prefix, suffix
  integer, intent(in) :: timestep, rank
  character(len=20)  :: chkuidstr

  write(chkuidstr, FMT=restore_fmt_t) timestep, chk_uid

  write(buffer,"('./',A,'/',A,'/',A,'_',A,'_',A,'_p',I6.6,A)") & 
    trim(folder), trim(cpfolder), trim(prefix), trim(gr_out_file), trim(chkuidstr), rank, trim(suffix)
end subroutine lbe_make_filename_cp_rank

!> Makes a filename by appending the given suffix, plus CPU rank and
!> timestep information, to the stem formed from the prefix and the value
!> of \c gr_out_file.
!> 
!> \c buffer will return this filename.
!>  This would be trivial in a sane language.
subroutine lbe_make_filename_rank(buffer, prefix, suffix, timestep, rank)
  implicit none
  character(len=*), intent(out) :: buffer
  character(len=*), intent(in) :: prefix, suffix
  integer, intent(in) :: timestep, rank
  character(len=20)  :: chkuidstr

  write(chkuidstr, FMT=restore_fmt_t) timestep, chk_uid

  write(buffer,"('./',A,'/',A,'_',A,'_',A,'_p',I6.6,A)") & 
    trim(folder), trim(prefix), trim(gr_out_file), trim(chkuidstr), rank, trim(suffix)
end subroutine lbe_make_filename_rank

!> Makes a filename by appending the given suffix, plus CPU rank and
!> timestep information, to the stem formed from the prefix and the value
!> of \c gr_out_file.
!> 
!> \c buffer will return this filename.
!>  This would be trivial in a sane language.
subroutine lbe_make_filename_restore_rank(buffer, prefix, suffix, rank)
  implicit none

  character(len=*), intent(out) :: buffer
  character(len=*), intent(in) :: prefix, suffix
  integer, intent(in) :: rank

  if (trim(srccpfolder) .ne. "") then
    write(buffer,"(A,'/',A,'_',A,'_',A,'_p',I6.6,A)") & 
      trim(srccpfolder), trim(prefix), trim(gr_out_file), trim(restore_string), rank, trim(suffix)
  else
    write(buffer,"('./',A,'/',A,'/',A,'_',A,'_',A,'_p',I6.6,A)") & 
      trim(folder), trim(cpfolder), trim(prefix), trim(gr_out_file), trim(restore_string), rank, trim(suffix)
  endif

end subroutine lbe_make_filename_restore_rank

!> Makes a filename by appending the given suffix and
!> timestep information, to the stem formed from the prefix and the value
!> of \c gr_out_file.
!> 
!> \c buffer will return this filename.
!>  This would be trivial in a sane language.
subroutine lbe_make_filename_cp(buffer, prefix, suffix, timestep)
  character(len=*), intent(out) :: buffer
  character(len=*), intent(in) :: prefix, suffix
  integer, intent(in) :: timestep
  character(len=20)  :: chkuidstr

  write(chkuidstr, FMT=restore_fmt_t) timestep, chk_uid

  write(buffer,"('./',A,'/',A,'/',A,'_',A,'_',A,A)") & 
    trim(folder), trim(cpfolder), trim(prefix), trim(gr_out_file), trim(chkuidstr), trim(suffix)
end subroutine lbe_make_filename_cp


subroutine lbe_make_filename_restore(buffer, prefix, suffix)
  implicit none

  character(len=*), intent(out) :: buffer
  character(len=*), intent(in) :: prefix, suffix

  if (trim(srccpfolder) .ne. "") then
    write(buffer,"(A,'/',A,'_',A,'_',A,A)") & 
      trim(srccpfolder), trim(prefix), trim(gr_out_file), trim(restore_string), trim(suffix)
  else
    write(buffer,"('./',A,'/',A,'/',A,'_',A,'_',A,A)") & 
      trim(folder), trim(cpfolder), trim(prefix), trim(gr_out_file), trim(restore_string), trim(suffix)
  endif
end subroutine lbe_make_filename_restore

subroutine lbe_make_filename_append(buffer, prefix, suffix)
  implicit none

  character(len=*), intent(out) :: buffer
  character(len=*), intent(in) :: prefix, suffix
  character(len=10)  :: chkuidstr

  write(chkuidstr,"(I10.10)") chk_uid

  write(buffer,"('./',A,'/',A,'_',A,'-',A,A)") &
    trim(folder), trim(prefix), trim(gr_out_file), trim(chkuidstr), trim(suffix)

end subroutine lbe_make_filename_append

!>  This routine writes an AVS field file
!>  Variable veclen determines fiedltype:
!>  1=scalar,2=2scalar,3=3scalar,4=vector
!> 
!>  Added 21.06.02 by Jens
subroutine dump_avs_fld(prefix,nxi,nyi,nzi,vectmp)
  implicit none
  character(len=1024) :: filename,fldname
  character(len=1024) :: tmpstring
  character(len=*) :: prefix
  character(len=1)   :: vecstr,countstr,stridestr
  character(len=7)   :: skipstr
  integer :: nxi,nyi,nzi,veclen,i,i2,i3,vectmp

  ! Stupid compiler bug workaround on IRIX
  veclen = vectmp

  call lbe_make_filename_output(fldname,prefix,'.fld',nt)
  open(10,file=fldname)
  write(10,'(a)') '# AVS field file'
  write(10,'(a)') 'ndim=3'
  write(10,'(a,i3)') 'dim1=',nxi
  write(10,'(a,i3)') 'dim2=',nyi
  write(10,'(a,i3)') 'dim3=',nzi
  write(10,'(a)') 'nspace=3'
  write(10,'(a)') 'field=uniform'

  if (index(dump_format,'all').gt.0) then
    call lbe_make_filename_output(filename,prefix,'.all',nt)
    if (veclen.eq.4) veclen=3
    write(stridestr,'(i1.1)') veclen+3
    write(10,'(a)') 'data=float'
    write(vecstr,'(i1.1)') veclen
    write(10,'(a)') 'veclen='//vecstr
      do i=1,veclen
        write(countstr,'(i1.1)') i 
        write(skipstr,'(i3.3)') i+2
        tmpstring='variable '//trim(countstr)//' file='//      &
            trim(filename)//' filetype=ascii skip=0 offset='// &
            trim(skipstr)//' stride='//trim(stridestr)
        write(10,'(a,a,a)') trim(tmpstring)
      end do
  else
    if (index(dump_format,'xdr').gt.0) then
      if (veclen.eq.4) veclen=3
      call lbe_make_filename_output(filename,prefix,'.xdr',nt)
      write(vecstr,'(i1.1)') veclen
      write(10,'(a)') 'veclen='//vecstr
      if (dump_double) then
        write(10,'(a)') 'data=xdr_double'
      else
        write(10,'(a)') 'data=xdr_float'
      end if
      do i=1,veclen
        write(countstr,'(i1.1)') i
        write(skipstr,'(i7.7)') 8*(i-1)
        tmpstring='variable '//countstr//' file='//trim(filename)// &
          ' filetype=binary skip='//skipstr//' stride='//vecstr
        write(10,'(a,a,a)') trim(tmpstring)
      end do
    else
      call lbe_make_filename_output(filename,prefix,'.bin',nt)

      if (dump_double) then
        write(10,'(a)') 'data=double'
        i3 = 8
      else
        write(10,'(a)') 'data=float'
        i3 = 4
      end if

      i2 = veclen
      if (veclen.eq.4) veclen=3
      write(vecstr,'(i1.1)') veclen
      write(10,'(a)') 'veclen='//vecstr

      do i=1,veclen
        write(countstr,'(i1.1)') i
          if (i2.lt.4) then
          write(skipstr,'(i7.7)') i3*nxi*nyi*nzi*(i-1)
          stridestr='1'
          else
          write(skipstr,'(i7.7)') i3*(i-1)
          write(stridestr,'(i1.1)') veclen
          end if
        tmpstring='variable '//countstr//' file='//trim(filename)// &
        ' filetype=unformatted skip='//skipstr//&
        ' stride='//trim(stridestr)
        write(10,'(a,a,a)') trim(tmpstring)
      end do
    endif
  endif

  close(10)
end subroutine dump_avs_fld


!>Checks to see if the given file exists. If not,
!>prints its name and an error message, and aborts.
!>Otherwise, silently returns.
subroutine die_unless_exists(filename)
  character(len=*), intent(in) :: filename
  logical :: exists_p
  inquire(file=filename ,exist = exists_p)
  if (.not. exists_p) then
    write(msgstr,"('File <',A,'> not found. Aborting...')") trim(filename)
    call error(msgstr)
  end if
  write(msgstr,"('File <',A,'> found OK.')") trim(filename)
  call log_msg(msgstr)
end subroutine die_unless_exists

subroutine die_unless_nonzero(ierr,string)
  integer :: ierr
  character(len=*),intent(in) :: string

  if (ierr .eq. 0) then
    print*,string
    call abend
  end if
end subroutine die_unless_nonzero

end module lbe_io_helper_module
