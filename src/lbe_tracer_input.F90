#include "lbe.h"

!> input for \c TRACER part of the code
module lbe_tracer_input_module
#ifdef TRACER
    use lbe_globals_module, only: input_dfile_unit,maxpos,minpos,tsize,myrankc
    use lbe_io_helper_module, only: lbe_make_filename_restore
    use lbe_log_module
    use lbe_parallel_module, only: check_allocate,comm_cart
    use lbe_parms_module, only: chk_uid,inp_file,nt,restore_string&
         &,arg_input_dfile_set,arg_input_dfile
    use lbe_tracer_globals_module
    use lbe_tracer_helper_module, only: log_msg_tracer,error_tracer&
         &,log_msg_tracer_hdr
    use lbe_tracer_init_module, only: place_mod_pos,t_placement
    use lbe_tracer_module, only: setup_list
    use lbe_tracer_parallel_module, only: scatter_tracers
    use lbe_tracer_output_module, only: dump_double,dump_format,dump_ids&
         &,dump_kinds,dump_mod_uid,dump_positions,dump_velocities&
	 &,msdx_name,msdx_subtract_drift&
	 &,n_msdx,n_msdx_max,n_msdx_sample&
	 &,msdx_len,msdx_minx,msdx_miny,msdx_minz,msdx_maxx,msdx_maxy,msdx_maxz&
	 &,msdy_name,msdy_subtract_drift&
	 &,n_msdy,n_msdy_max,n_msdy_sample&
         &,msdy_len,msdy_minx,msdy_miny,msdy_minz,msdy_maxx,msdy_maxy,msdy_maxz&
	 &,msdz_name,msdz_subtract_drift&
	 &,n_msdz,n_msdz_max,n_msdz_sample&
	 &,msdz_len,msdz_minx,msdz_miny,msdz_minz,msdz_maxx,msdz_maxy,msdz_maxz&
         &,n_dump,n_flow_rate&
         &,n_profile&
         &,tracer_flow_rate_restore_checkpoint_xdr
    use lbe_tracer_recoloring_module, only: x_recolor_plane_1,x_recolor_plane_2

    implicit none
    include 'mpif.h'
    private

    public read_tracer_input,tracer_restore_checkpoint

    namelist /tracer_input/ dump_double,dump_format,dump_ids,dump_kinds&
         &,dump_mod_uid,dump_positions,dump_velocities&
	 &,msdx_name&
	 &,n_msdx,n_msdx_sample&
	 &,msdx_len,msdx_minx,msdx_miny,msdx_minz,msdx_maxx,msdx_maxy,msdx_maxz&
         &,msdx_subtract_drift&
	 &,msdy_name&
         &,msdy_len,msdy_minx,msdy_miny,msdy_minz,msdy_maxx,msdy_maxy,msdy_maxz&
         &,msdy_subtract_drift&
	 &,n_msdy,n_msdy_sample&
	 &,msdz_name&
	 &,msdz_len,msdz_minx,msdz_miny,msdz_minz,msdz_maxx,msdz_maxy,msdz_maxz&
<<<<<<< HEAD
         &,msdz_subtract_drift&
=======
>>>>>>> c5db403... bug fix on tracer MSD calculation
	 &,n_msdz,n_msdz_sample&
         &,n_dump,n_flow_rate&
         &,n_profile,place_mod_pos,t_placement,x_recolor_plane_1&
         &,x_recolor_plane_2

contains

    !> read \c TRACER input file and broadcast its content
    subroutine read_tracer_input
        integer i,ierror

        call log_msg_tracer_hdr("Reading TRACER input")

        if (myrankc.eq.0) then
           open (unit=tracer_input_file_unit,file=trim(inp_file)//'.tracer'&
                &,err=100)
           read (unit=tracer_input_file_unit,nml=tracer_input,err=110)
           close (unit=tracer_input_file_unit,err=120)
        endif

        if ( arg_input_dfile_set ) then
           call log_msg_tracer("  Getting differential input...")
           open (unit=input_dfile_unit,file=arg_input_dfile,status='UNKNOWN')
           read (unit=input_dfile_unit,nml=tracer_input,iostat=ierror)
           if (ierror/=0) then
              call log_msg_tracer('    WARNING: Differential namelist not '&
                   &//'found or errors encountered.')
           end if
           close (unit=input_dfile_unit)
           call log_ws()
        end if

        write (msgstr,"('place_mod_pos      = ',I0)") place_mod_pos
        call log_msg_tracer(msgstr)
        write (msgstr,"('t_placement        = ',I0)") t_placement
        call log_msg_tracer(msgstr)
        write (msgstr,"('n_dump             = ',I0,', dump_format = <',A,'>')")&
             &n_dump, trim(dump_format)
        call log_msg_tracer(msgstr)
        write (msgstr,"('dump_double        = ',L1)") dump_double
        call log_msg_tracer(msgstr)
        write (msgstr,"('dump_mod_uid       = ',I0)") dump_mod_uid
        call log_msg_tracer(msgstr)
        write (msgstr&
             &,"('dump_ids           = ',L1,', dump_positions    = ',L1)") &
             &dump_ids,dump_positions
        call log_msg_tracer(msgstr)
        write (msgstr&
             &,"('dump_velocities    = ',L1,', dump_kinds       = ',L1)") &
             &dump_velocities,dump_kinds
        call log_msg_tracer(msgstr)
        write (msgstr,"('n_profile          = ',I0)") n_profile
        call log_msg_tracer(msgstr)
        write (msgstr,"('n_flow_rate        = ',I0)") n_flow_rate
        call log_msg_tracer(msgstr)
        write (msgstr,"('x_recolor_plane_1  = ',F16.10)") x_recolor_plane_1
        CALL log_msg_tracer(msgstr)
        write (msgstr,"('x_recolor_plane_2  = ',F16.10)") x_recolor_plane_2
        CALL log_msg_tracer(msgstr)

        write (msgstr,"('n_msdx             = ',I0)") n_msdx
        call log_msg_tracer(msgstr)
        if (n_msdx/=0) then
           write (msgstr,"('n_msdx_sample      = ',I0)") n_msdx_sample
           call log_msg_tracer(msgstr)
           write (msgstr,"('msdx_len           = ',I0)") msdx_len
           call log_msg_tracer(msgstr)
           if (all(len_trim(msdx_name)==0)) call log_msg_tracer(&
                &'All mean square x-displacement regions are disabled.')
           do i=1,n_msdx_max
              if (len_trim(msdx_name(i))==0) cycle

              write (msgstr,&
                   &"('msdx region ""',A,'"" referred to as <tracer-x-',A,'>')"&
                   &) trim(msdx_name(i)),trim(msdx_name(i))
              call log_msg_tracer(msgstr)
              write (msgstr,"('  msdx_min[xyz]      =',3(X,F16.10))") &
                   &msdx_minx(i),msdx_miny(i),msdx_minz(i)
              call log_msg_tracer(msgstr)
              write (msgstr,"('  msdx_max[xyz]      =',3(X,F16.10))") &
                   &msdx_maxx(i),msdx_maxy(i),msdx_maxz(i)
              call log_msg_tracer(msgstr)
              if (any((/msdx_minx(i),msdx_miny(i),msdx_minz(i)&
                   &,msdx_maxx(i),msdx_maxy(i),msdx_maxz(i)/)<0.0_rk)) then
                 write (msgstr,"(' -> defaulting to')")
                 call log_msg_tracer(msgstr)

                 if (msdx_minx(i)<0.0_rk) msdx_minx(i) = minpos(1)
                 if (msdx_maxx(i)<0.0_rk) msdx_maxx(i) = maxpos(1)
                 if (msdx_miny(i)<0.0_rk) msdx_miny(i) = minpos(2)
                 if (msdx_maxy(i)<0.0_rk) msdx_maxy(i) = maxpos(2)
                 if (msdx_minz(i)<0.0_rk) msdx_minz(i) = minpos(3)
                 if (msdx_maxz(i)<0.0_rk) msdx_maxz(i) = maxpos(3)

                 write (msgstr,"('  msdx_min[xyz]      =',3(X,F16.10))") &
                      &msdx_minx(i),msdx_miny(i),msdx_minz(i)
                 call log_msg_tracer(msgstr)
                 write (msgstr,"('  msdx_max[xyz]      =',3(X,F16.10))") &
                      &msdx_maxx(i),msdx_maxy(i),msdx_maxz(i)
                 call log_msg_tracer(msgstr)
              end if
              write (msgstr,"('  msdx_subtract_drift = ',L1)") &
                   &msdx_subtract_drift(i)
              call log_msg_tracer(msgstr)
           end do
        end if

        write (msgstr,"('n_msdy             = ',I0)") n_msdy
        call log_msg_tracer(msgstr)
        if (n_msdy/=0) then
           write (msgstr,"('n_msdy_sample      = ',I0)") n_msdy_sample
           call log_msg_tracer(msgstr)
           write (msgstr,"('msdy_len           = ',I0)") msdy_len
           call log_msg_tracer(msgstr)
           if (all(len_trim(msdy_name)==0)) call log_msg_tracer(&
                &'All mean square y-displacement regions are disabled.')
           do i=1,n_msdy_max
              if (len_trim(msdy_name(i))==0) cycle

              write (msgstr,&
                   &"('msdy region ""',A,'"" referred to as <tracer-y-',A,'>')"&
                   &) trim(msdy_name(i)),trim(msdy_name(i))
              call log_msg_tracer(msgstr)
              write (msgstr,"('  msdy_min[xyz]      =',3(X,F16.10))") &
                   &msdy_minx(i),msdy_miny(i),msdy_minz(i)
              call log_msg_tracer(msgstr)
              write (msgstr,"('  msdy_max[xyz]      =',3(X,F16.10))") &
                   &msdy_maxx(i),msdy_maxy(i),msdy_maxz(i)
              call log_msg_tracer(msgstr)
              if (any((/msdy_minx(i),msdy_miny(i),msdy_minz(i)&
                   &,msdy_maxx(i),msdy_maxy(i),msdy_maxz(i)/)<0.0_rk)) then
                 write (msgstr,"(' -> defaulting to')")
                 call log_msg_tracer(msgstr)

                 if (msdy_minx(i)<0.0_rk) msdy_minx(i) = minpos(1)
                 if (msdy_maxx(i)<0.0_rk) msdy_maxx(i) = maxpos(1)
                 if (msdy_miny(i)<0.0_rk) msdy_miny(i) = minpos(2)
                 if (msdy_maxy(i)<0.0_rk) msdy_maxy(i) = maxpos(2)
                 if (msdy_minz(i)<0.0_rk) msdy_minz(i) = minpos(3)
                 if (msdy_maxz(i)<0.0_rk) msdy_maxz(i) = maxpos(3)

                 write (msgstr,"('  msdy_min[xyz]      =',3(X,F16.10))") &
                      &msdy_minx(i),msdy_miny(i),msdy_minz(i)
                 call log_msg_tracer(msgstr)
                 write (msgstr,"('  msdy_max[xyz]      =',3(X,F16.10))") &
                      &msdy_maxx(i),msdy_maxy(i),msdy_maxz(i)
                 call log_msg_tracer(msgstr)
              end if
              write (msgstr,"('  msdy_subtract_drift = ',L1)") &
                   &msdy_subtract_drift(i)
              call log_msg_tracer(msgstr)
           end do
        end if

        write (msgstr,"('n_msdz             = ',I0)") n_msdz
        call log_msg_tracer(msgstr)
        if (n_msdz/=0) then
           write (msgstr,"('n_msdz_sample      = ',I0)") n_msdz_sample
           call log_msg_tracer(msgstr)
           write (msgstr,"('msdz_len           = ',I0)") msdz_len
           call log_msg_tracer(msgstr)
           if (all(len_trim(msdz_name)==0)) call log_msg_tracer(&
                &'All mean square z-displacement regions are disabled.')
           do i=1,n_msdz_max
              if (len_trim(msdz_name(i))==0) cycle

              write (msgstr,&
                   &"('msdz region ""',A,'"" referred to as <tracer-z-',A,'>')"&
                   &) trim(msdz_name(i)),trim(msdz_name(i))
              call log_msg_tracer(msgstr)
              write (msgstr,"('  msdz_min[xyz]      =',3(X,F16.10))") &
                   &msdz_minx(i),msdz_miny(i),msdz_minz(i)
              call log_msg_tracer(msgstr)
              write (msgstr,"('  msdz_max[xyz]      =',3(X,F16.10))") &
                   &msdz_maxx(i),msdz_maxy(i),msdz_maxz(i)
              call log_msg_tracer(msgstr)
              if (any((/msdz_minx(i),msdz_miny(i),msdz_minz(i)&
                   &,msdz_maxx(i),msdz_maxy(i),msdz_maxz(i)/)<0.0_rk)) then
                 write (msgstr,"(' -> defaulting to')")
                 call log_msg_tracer(msgstr)

                 if (msdz_minx(i)<0.0_rk) msdz_minx(i) = minpos(1)
                 if (msdz_maxx(i)<0.0_rk) msdz_maxx(i) = maxpos(1)
                 if (msdz_miny(i)<0.0_rk) msdz_miny(i) = minpos(2)
                 if (msdz_maxy(i)<0.0_rk) msdz_maxy(i) = maxpos(2)
                 if (msdz_minz(i)<0.0_rk) msdz_minz(i) = minpos(3)
                 if (msdz_maxz(i)<0.0_rk) msdz_maxz(i) = maxpos(3)

                 write (msgstr,"('  msdz_min[xyz]      =',3(X,F16.10))") &
                      &msdz_minx(i),msdz_miny(i),msdz_minz(i)
                 call log_msg_tracer(msgstr)
                 write (msgstr,"('  msdz_max[xyz]      =',3(X,F16.10))") &
                      &msdz_maxx(i),msdz_maxy(i),msdz_maxz(i)
                 call log_msg_tracer(msgstr)
              end if
              write (msgstr,"('  msdz_subtract_drift = ',L1)") &
                   &msdz_subtract_drift(i)
              call log_msg_tracer(msgstr)
           end do
        end if

        call log_ws()

        call MPI_Bcast(place_mod_pos,1,MPI_INTEGER,0,comm_cart,ierror)
        call MPI_Bcast(t_placement,1,MPI_INTEGER,0,comm_cart,ierror)
        call MPI_Bcast(n_dump,1,MPI_INTEGER,0,comm_cart,ierror)
        call MPI_Bcast(n_flow_rate,1,MPI_INTEGER,0,comm_cart,ierror)
        call MPI_Bcast(n_profile,1,MPI_INTEGER,0,comm_cart,ierror)
        call MPI_Bcast(dump_mod_uid,1,MPI_INTEGER,0,comm_cart,ierror)
        call MPI_Bcast(dump_format,32,MPI_CHARACTER,0,comm_cart,ierror)
        call MPI_Bcast(dump_double,1,MPI_LOGICAL,0,comm_cart,ierror)
        call MPI_Bcast(dump_ids,1,MPI_LOGICAL,0,comm_cart,ierror)
        call MPI_Bcast(dump_kinds,1,MPI_LOGICAL,0,comm_cart,ierror)
        call MPI_Bcast(dump_positions,1,MPI_LOGICAL,0,comm_cart,ierror)
        call MPI_Bcast(dump_velocities,1,MPI_LOGICAL,0,comm_cart,ierror)
        call MPI_Bcast(x_recolor_plane_1,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(x_recolor_plane_2,1,MPI_REAL8,0,comm_cart,ierror)

        call MPI_Bcast(n_msdx,1,MPI_INTEGER,0,comm_cart,ierror)
        call MPI_Bcast(n_msdx_sample,1,MPI_INTEGER,0,comm_cart,ierror)
        call MPI_Bcast(msdx_len,1,MPI_INTEGER,0,comm_cart,ierror)

        do i=1,n_msdx_max
           call MPI_Bcast(msdx_name(i),16,MPI_CHARACTER,0,comm_cart,ierror)
        end do
        call MPI_Bcast(msdx_minx,n_msdx_max,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(msdx_maxx,n_msdx_max,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(msdx_miny,n_msdx_max,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(msdx_maxy,n_msdx_max,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(msdx_minz,n_msdx_max,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(msdx_maxz,n_msdx_max,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(msdx_subtract_drift,n_msdx_max,MPI_LOGICAL,0,comm_cart&
             &,ierror)

        call MPI_Bcast(n_msdy,1,MPI_INTEGER,0,comm_cart,ierror)
        call MPI_Bcast(n_msdy_sample,1,MPI_INTEGER,0,comm_cart,ierror)
        call MPI_Bcast(msdy_len,1,MPI_INTEGER,0,comm_cart,ierror)

        do i=1,n_msdy_max
           call MPI_Bcast(msdy_name(i),16,MPI_CHARACTER,0,comm_cart,ierror)
        end do
        call MPI_Bcast(msdy_minx,n_msdy_max,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(msdy_maxx,n_msdy_max,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(msdy_miny,n_msdy_max,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(msdy_maxy,n_msdy_max,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(msdy_minz,n_msdy_max,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(msdy_maxz,n_msdy_max,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(msdy_subtract_drift,n_msdy_max,MPI_LOGICAL,0,comm_cart&
             &,ierror)

        call MPI_Bcast(n_msdz,1,MPI_INTEGER,0,comm_cart,ierror)
        call MPI_Bcast(n_msdz_sample,1,MPI_INTEGER,0,comm_cart,ierror)
        call MPI_Bcast(msdz_len,1,MPI_INTEGER,0,comm_cart,ierror)

        do i=1,n_msdz_max
           call MPI_Bcast(msdz_name(i),16,MPI_CHARACTER,0,comm_cart,ierror)
        end do
        call MPI_Bcast(msdz_minx,n_msdz_max,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(msdz_maxx,n_msdz_max,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(msdz_miny,n_msdz_max,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(msdz_maxy,n_msdz_max,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(msdz_minz,n_msdz_max,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(msdz_maxz,n_msdz_max,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(msdz_subtract_drift,n_msdz_max,MPI_LOGICAL,0,comm_cart&
             &,ierror)

        return
100     continue
        call error_tracer(&
             &'Error opening TRACER input file "'//trim(inp_file)//'.tracer"')
        return
110     continue
        call error_tracer(&
             &'Error reading TRACER input file "'//trim(inp_file)//'.tracer"')
        return
120     continue
        call error_tracer(&
             &'Error closing TRACER input file "'//trim(inp_file)//'.tracer"')
        return
    end subroutine read_tracer_input

    !> restores from an \c TRACER checkpoint file
    subroutine tracer_restore_checkpoint
        type(tracer_particle_type),allocatable,dimension(:) :: tp
        character(len=1024) chk_file_name
        integer i,chk_file_handle,restore_uid,ierror,n_global,stat
        real(kind=rk) :: pf(3)

        ! This would be a function pointer in C. In fact,  xdrfdouble  is
        ! an external function, but fortran does not seem to care about the
        ! type anyway...
        integer,external :: xdrfdouble

#ifdef USEXDRF
        rank_0: if (myrankc==0) then
           call lbe_make_filename_restore(chk_file_name,'tracer-checkpoint'&
                &,'.xdr')

           call xdrfopen(chk_file_handle,trim(chk_file_name),'r',ierror)
           if (ierror==0) then
              ! Checking for XDRSLOPPY probably makes no sense here,
              ! the checkpoint MUST be read.
              call error_tracer(&
                   &'tracer_restore_checkpoint(): xdrfopen() failed opening "'&
                   &//trim(chk_file_name)//'"')
           end if

           ! ...at last, read particle data:
           call xdrfint(chk_file_handle,n_global,ierror)

           allocate (tp(n_global),stat=stat)
           call check_allocate(stat,'tracer_restore_checkpoint(): tp')

           do i=1,n_global
              call xdrfvector(chk_file_handle,tp(i)%x(1:3),3,xdrfdouble,ierror)
              call xdrfvector(chk_file_handle,tp(i)%v(1:3),3,xdrfdouble,ierror)
              call xdrfint(chk_file_handle,tp(i)%uid,ierror)
              call xdrfint(chk_file_handle,tp(i)%kind,ierror)
           end do
        end if rank_0

        call tracer_flow_rate_restore_checkpoint_xdr(chk_file_handle)

        if (myrankc==0) then
           call xdrfclose(chk_file_handle,ierror)
        end if

        call setup_list
        call scatter_tracers(n_global,tp)

        if (myrankc==0) then
           deallocate (tp)
           write(msgstr&
                &,"('Placed a total of ',I0,' tracers from checkpoint.')") &
                &n_global
           call log_msg_tracer(msgstr)
        end if
#else
        call error_tracer('failed to restore from checkpoint file'&
             &//' - recompile with USEXDRF set!')
#endif
    end subroutine tracer_restore_checkpoint

#endif
end module lbe_tracer_input_module
