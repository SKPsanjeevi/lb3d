#include "lbe.h"

!> Dumping of velocity field statistics
module lbe_io_arrstats_module
    use lbe_globals_module
    use lbe_helper_module, only: every_n_time_steps, velocity
    use lbe_log_module
#ifdef USEHDF
    use lbe_io_hdf5_module, only: read_iscalar_phdf5,read_vector_phdf5
#endif
    use lbe_io_helper_module, only: lbe_delete_file,lbe_make_filename_cp&
         &,lbe_make_filename_output,lbe_make_filename_restore
    use lbe_io_module, only: dump_iscalar,dump_vector
#ifdef USEXDRF
    use lbe_io_xdrf_module, only: check_xdrfopen
#endif
    use lbe_parallel_module, only: abend, check_allocate
    use lbe_parms_module, only: arrstats_intervals,arrstats_name,dump_format&
         &,n_sci_arrstats_dump,nt,nx,ny,nz,sci_arrstats
    use lbe_types_module, only: lbe_site
    implicit none
    private
    public arrstats_delete_checkpoint,dump_arrstats,dump_arrstats_checkpoint&
         &,lbe_io_setup_arrstats,restore_arrstats_checkpoint

    !> type representing data accumulated at one sampling interval
    type arrstats_data_type
       integer :: n_samples     !< number of sampled time steps
       !> number of samples per lattice site
       integer,allocatable,dimension(:,:,:) :: n
       !> accumulation buffer, contains 9 elements per lattice, the accumulated
       !> \c v(1:3), \c v(1:3)*v(1:3), and the quadratic cross terms \c v(1)*v(
       real(kind=rk),allocatable,dimension(:,:,:,:) :: m
    end type arrstats_data_type

    !> data accumulated at different sampling intervals
    type(arrstats_data_type),save :: arrstats_data(arrstats_intervals)

contains

    !> delete old arrstats checkpoint files.
    !>
    !> \param[in] last timestep of the checkpoint to delete
    subroutine arrstats_delete_checkpoint(last)
        integer, intent(in)  :: last
        character(len=1024) :: filename
        integer :: i

        if (.not.sci_arrstats) return

        call log_msg("  Deleting arrstats checkpoint files ...")

        do i=1,arrstats_intervals
           if (n_sci_arrstats_dump(i)==0) cycle

           if (myrankc==0) then
              ! this is always xdr format so far...
              call lbe_make_filename_output(filename&
                   &,'checkarrstats-desc-'//trim(arrstats_name(i)),'.xdr',last)
              call lbe_delete_file(filename)

              if (dump_format=='hdf') then
                 call lbe_make_filename_cp(filename&
                      &,'checkarrstats-n-'//trim(arrstats_name(i)),'.h5',last)
                 call lbe_delete_file(filename)

                 call lbe_make_filename_cp(filename&
                      &,'checkarrstats-m-'//trim(arrstats_name(i)),'.h5',last)
                 call lbe_delete_file(filename)
              else
                 call lbe_make_filename_cp(filename&
                      &,'checkarrstats-n-'//trim(arrstats_name(i)),'.xdr',last)
                 call lbe_delete_file(filename)

                 call lbe_make_filename_cp(filename&
                      &,'checkarrstats-m-'//trim(arrstats_name(i)),'.xdr',last)
                 call lbe_delete_file(filename)
              end if
           end if
        end do
    end subroutine arrstats_delete_checkpoint

    !> accumulates data at all sampling intervals, eventually dumps them
    !>
    !> \param[in] N local lattice chunk with halo of size 1
    subroutine dump_arrstats(N)
        type(lbe_site),intent(in) :: N(0:,0:,0:)
        integer i

        do i=1,arrstats_intervals
           if (n_sci_arrstats_dump(i)==0) cycle

           call sample_arrstats(N,i,arrstats_data(i))

           if (every_n_time_steps(n_sci_arrstats_dump(i))) then
              call log_msg(" -now dumping...")
              call dump_arrstats_data(i,arrstats_data(i))
              call reset_arrstats_data(arrstats_data(i))
           end if

        end do
    end subroutine dump_arrstats

    !> initialize \c arrstats output
    subroutine lbe_io_setup_arrstats
        integer i

        if (.not.sci_arrstats) return

        do i=1,arrstats_intervals
           if (n_sci_arrstats_dump(i)==0) cycle

           call init_arrstats_data(i,arrstats_data(i))
        end do
    end subroutine lbe_io_setup_arrstats

    !> initialize data objects related to a specific sampling interval
    !>
    !> \param[in] i index of the sampling interval (must be within \c
    !> (1:arrstats_intervals))
    !>
    !> \param[in,out] asd corresponding accumulation buffer
    subroutine init_arrstats_data(i,asd)
        integer,intent(in) :: i
        type(arrstats_data_type),intent(inout) :: asd
        integer :: stat

        if (trim(arrstats_name(i))=='') &
             &write (arrstats_name(i),fmt='("d",I0)') n_sci_arrstats_dump(i)

        allocate (asd%n(1:nx,1:ny,1:nz),asd%m(1:9,1:nx,1:ny,1:nz),stat=stat)
        call check_allocate(stat&
             &,'init_arrstats(): asd%n,asd%m (arrstats_name <'&
             &//trim(arrstats_name(i))//'>)')

        call reset_arrstats_data(asd)
    end subroutine init_arrstats_data

    !> samples the current velocity field into one accumulation buffer
    !>
    !> \param[in] N local lattice chunk (halo extent 1)
    !>
    !> \param[in] i index of the sampling interval (must be within \c
    !> (1:arrstats_intervals))
    !>
    !> \param[in,out] asd corresponding accumulation buffer
    subroutine sample_arrstats(N,i,asd)
        type(lbe_site),intent(in) :: N(0:,0:,0:)
        integer,intent(in) :: i
        type(arrstats_data_type),intent(inout) :: asd
        integer :: x,y,z
        real(kind=rk) :: v(3)

        call log_msg(' -"'//trim(arrstats_name(i))//'"')

        do x=1,nx
           do y=1,ny
              do z=1,nz
                 if (N(x,y,z)%rock_state==0.0_rk) then
                    asd%n(x,y,z) = asd%n(x,y,z)+1
                    v = velocity(N(x,y,z))
                    asd%m(1:3,x,y,z) = asd%m(1:3,x,y,z)+v ! velocity
                    asd%m(4:6,x,y,z) = asd%m(4:6,x,y,z)+v*v ! square terms
                    asd%m(7,x,y,z) = asd%m(7,x,y,z)+v(2)*v(1) ! cross terms...
                    asd%m(8,x,y,z) = asd%m(8,x,y,z)+v(3)*v(1)
                    asd%m(9,x,y,z) = asd%m(9,x,y,z)+v(3)*v(2)
                 end if
              end do
           end do
        end do
        asd%n_samples = asd%n_samples+1
    end subroutine sample_arrstats

    !> Writes the complete data related to one sampling interval to disk
    !>
    !> \param[in] i index of the sampling interval (must be within \c
    !> (1:arrstats_intervals))
    !>
    !> \param[in,out] asd corresponding accumulation buffer
    !>
    !> \param[in] pprefix prefix to add to the filenames (optional)
    subroutine dump_arrstats_data(i,asd,pprefix)
        integer,intent(in) :: i
        type(arrstats_data_type),intent(inout) :: asd
        character(len=*),intent(in),optional :: pprefix
        character(len=32) :: prefix

        if (present(pprefix)) then
           prefix = pprefix
        else
           prefix = ''
        end if

        call dump_arrstats_desc(i,asd,prefix)
        call dump_iscalar(asd%n&
             &,trim(prefix)//'arrstats-n-'//trim(arrstats_name(i)))
        call dump_vector(asd%m&
             &,trim(prefix)//'arrstats-m-'//trim(arrstats_name(i)))
    end subroutine dump_arrstats_data

    !> Writes additional descriptive data regarding a state of
    !> accumulated data to disk
    !>
    !> \param[in] i index of the sampling interval (must be within \c
    !> (1:arrstats_intervals))
    !>
    !> \param[in] asd corresponding accumulation buffer
    !>
    !> \param[in] pprefix prefix to add to the filenames (can be \c '')
    !>
    !> At the moment the file format is always XDR. The data written
    !> encompasses the number of sampled time steps and the sampling
    !> interval. From the first, the time interval covered by the data
    !> can be obtained if \c n_sci_arrstats is known.
    subroutine dump_arrstats_desc(i,asd,prefix)
        integer,intent(in) :: i
        type(arrstats_data_type),intent(in) :: asd
        character(len=32),intent(in) :: prefix
        character(len=1024) :: filename
        integer :: file_id,ierror

        if (myrankc/=0) return

        call lbe_make_filename_output(filename&
             &,trim(prefix)//'arrstats-desc-'//trim(arrstats_name(i)),'.xdr',nt)
#ifndef USEXDRF
        call error('dump_arrstats_desc() requires still XDRF')
#else
        call xdrfopen(file_id,filename,"w",ierror)
        call check_xdrfopen(ierror,filename)

        call xdrfint(file_id,asd%n_samples,ierror)
        call xdrfint(file_id,n_sci_arrstats_dump(i),ierror)

        call xdrfclose(file_id,ierror)
#endif
    end subroutine dump_arrstats_desc

    !> Restores the state accumulated data for one sampling interval
    !> from a checkpoint
    !>
    !> \param[in] i index of the sampling interval (must be within \c
    !> (1:arrstats_intervals))
    !>
    !> \param[in,out] asd corresponding accumulation buffer
    subroutine restore_arrstats_desc(i,asd)
        integer,intent(in) :: i
        type(arrstats_data_type),intent(inout) :: asd
        character(len=1024) :: filename
        integer :: dummy_n_sci_arrstats_dump,file_id,ierror

        if (myrankc/=0) return

        call lbe_make_filename_restore(filename&
             &,'checkarrstats-desc-'//trim(arrstats_name(i)),'.xdr')
#ifndef USEXDRF
        call error('restore_arrstats_desc() requires still XDRF')
        dummy_n_sci_arrstats_dump = -1 ! this just avoids Cray compiler warnings
#else
        call xdrfopen(file_id,filename,'r',ierror)
        call check_xdrfopen(ierror,filename)

        call xdrfint(file_id,asd%n_samples,ierror)
        call xdrfint(file_id,dummy_n_sci_arrstats_dump,ierror)

        call xdrfclose(file_id,ierror)
#endif

        if (dummy_n_sci_arrstats_dump/=n_sci_arrstats_dump(i)) then
           write (unit=msgstr,fmt='("WARNING: read n_sci_arrstats_dump==",I0,'&
                &//'" from <",A,"> (which is ignored) but the input files '&
                &//'set n_sci_arrstats_dump==",I0," already.")') &
                &dummy_n_sci_arrstats_dump,filename,n_sci_arrstats_dump(i)
           call log_msg(msgstr)
        end if
    end subroutine restore_arrstats_desc

    !> Resets an instance of \c arrstats_data_type for further accumulation
    !>
    !> \param[in,out] asd corresponding accumulation buffer
    subroutine reset_arrstats_data(asd)
        type(arrstats_data_type),intent(inout) :: asd

        asd%n_samples = 0
        asd%n = 0
        asd%m = 0.0_rk
    end subroutine reset_arrstats_data

    !> Writes checkpoints for all active sampling intervals
    subroutine dump_arrstats_checkpoint()
        integer i,n_wrote

        if (.not.sci_arrstats) return

        call log_msg("Writing arrstats checkpoint ...")

        n_wrote = 0
        do i=1,arrstats_intervals
           if (n_sci_arrstats_dump(i)==0) cycle
           call dump_arrstats_data(i,arrstats_data(i),'check')
           n_wrote = n_wrote+1
        end do

        write (unit=msgstr,fmt='("  Wrote ",I0," arrstats checkpoints.")') &
             &n_wrote
        call log_msg(msgstr)

        call log_msg("Finished writing arrstats checkpoint.")

    end subroutine dump_arrstats_checkpoint

    !> Restore all accumulated data from a checkpoint
    subroutine restore_arrstats_checkpoint
        character(len=1024) :: filename
        logical :: exists
        integer :: i

        if (.not.sci_arrstats) return

        do i=1,arrstats_intervals
           if (n_sci_arrstats_dump(i)==0) cycle

           call lbe_make_filename_restore(filename&
                &,'checkarrstats-desc-'//trim(arrstats_name(i)),'.xdr')
           inquire (file=filename,exist=exists)
           file_exists: if (exists) then
              call restore_arrstats_desc(i,arrstats_data(i))

              call lbe_make_filename_restore(filename&
                   &,'checkarrstats-n-'//trim(arrstats_name(i)),'.h5')
#ifdef USEHDF
              call read_iscalar_phdf5(arrstats_data(i)%n,filename)
#else
              call log_msg("Need HDF5 to restore arrstats")
              call Abend
#endif

              call lbe_make_filename_restore(filename&
                   &,'checkarrstats-m-'//trim(arrstats_name(i)),'.h5')
#ifdef USEHDF
              call read_vector_phdf5(arrstats_data(i)%m,filename&
                   &,size(arrstats_data(i)%m,1))
#else
              call log_msg("Need HDF5 to restore arrstats")
              call Abend
#endif

              write (unit=msgstr&
                   &,fmt='("Restored ",I0," samples for arrstats <",A,">.")') &
                   &arrstats_data(i)%n_samples,trim(arrstats_name(i))
              call log_msg(msgstr)
           else file_exists
              call log_msg('arrstats checkpoint file <'//trim(filename)&
                   &//'> not found---resetting arrstats <'&
                   &//trim(arrstats_name(i))//'> to zero!')
              call reset_arrstats_data(arrstats_data(i))
           end if file_exists
        end do
    end subroutine restore_arrstats_checkpoint

end module lbe_io_arrstats_module
