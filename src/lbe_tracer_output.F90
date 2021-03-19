#include "lbe.h"

!> output for \c TRACER part of the code
module lbe_tracer_output_module
#ifdef TRACER
    use lbe_globals_module, only: halo_extent,minpos,tsize,myrankc
    use lbe_helper_module, only: iif,is_fluid
    use lbe_log_module, only: msgstr
    use lbe_io_helper_module, only: lbe_delete_file,lbe_make_filename_cp&
         &,lbe_make_filename_output
    use lbe_mean_square_displacement_mod, only: close_trace,register_msd,sample
    use lbe_parallel_module, only: check_allocate,comm_cart,start,tnx&
         &,tny,tnz
    use lbe_parms_module, only: checkpoint_format,checkpoint_safe,gr_out_file&
         &,nt,num_chkp_files,nx,ny,nz,n_checkpoint
    use lbe_tracer_globals_module
    use lbe_tracer_helper_module, only: build_tracer_mpitype,count_tracers&
         &,error_tracer,gather_tracers,gather_tracers_provide_rbuf&
         &,log_msg_tracer,passed_plane
    use lbe_types_module, only: lbe_site

#ifdef USEXDRF
    use lbe_io_xdrf_module, only: check_xdrfsloppy
#endif

    implicit none
    private

    include 'mpif.h'

    public count_tracer_flow,dump_configuration,dump_tracer_profile,input_dump&
         &,report_tracer_flow_rate,setup_dump&
	 &,sample_msdx,sample_msdy,sample_msdz&
         &,tracer_dump_checkpoint,tracer_delete_checkpoint&
         &,tracer_flow_rate_restore_checkpoint_xdr

    integer,save,public :: n_dump=100 !< dump cfg how often? (0: never)
    !> report (and average) flow rate how often? (0: never)
    integer,save,public :: n_flow_rate=1000
    integer,save,public :: n_profile=1000 !< dump profile how often? (0: never)

    character(len=32),save,public :: dump_format='asc' !< dump file format
    logical,save,public :: dump_double=.false. !< dump doubles or floats?

    !> dump only tracers with \c mod(uid,dump_mod_uid)==0
    integer,save,public :: dump_mod_uid=1

    logical,save,public :: dump_positions=.true. !< include \c x into output?
    logical,save,public :: dump_velocities=.true. !< include \c v into output?
    logical,save,public :: dump_ids=.true. !< include \c uid into output?
    logical,save,public :: dump_kinds=.true. !< include \c kind into output?

    !> \{
    !> \name options for mean square x-displacement of tracers
    integer,save,public :: n_msdx=0 !< dump msdx how often? (0:never)
    integer,save,public :: n_msdx_sample=10 !< sample msdx how often
    !> number of data points in output (spacing between points is \c
    !> n_msdx_sample )
    integer,save,public :: msdx_len=100
    !> number of slots available to define sub-volumes for mean square
    !> x-displacement calculation
    integer,parameter,public :: n_msdx_max=10
    !> human readable name for a mean square x-displacement
    !> calculation region, a value of '' disables the slot
    character(len=16),save,public :: msdx_name(n_msdx_max)=''
    !> don't sample tracers at momentary x-positions below this value
    !> (negative value: automatically fill in \c 0.5 )
    real(kind=rk),save,public :: msdx_minx(n_msdx_max)=-1.0_rk
    !> don't sample tracers at momentary x-positions above this value
    !> (negative value: automatically fill in \c 0.5+tnx )
    real(kind=rk),save,public :: msdx_maxx(n_msdx_max)=-1.0_rk
    !> don't sample tracers at momentary y-positions below this value
    !> (negative value: automatically fill in \c 0.5 )
    real(kind=rk),save,public :: msdx_miny(n_msdx_max)=-1.0_rk
    !> don't sample tracers at momentary y-positions above this value
    !> (negative value: automatically fill in \c 0.5+tny )
    real(kind=rk),save,public :: msdx_maxy(n_msdx_max)=-1.0_rk
    !> don't sample tracers at momentary z-positions below this value
    !> (negative value: automatically fill in \c 0.5 )
    real(kind=rk),save,public :: msdx_minz(n_msdx_max)=-1.0_rk
    !> don't sample tracers at momentary z-positions above this value
    !> (negative value: automatically fill in \c 0.5+tnz )
    real(kind=rk),save,public :: msdx_maxz(n_msdx_max)=-1.0_rk
    !> correct for a possible drift of the whole system?
    logical,save,public :: msdx_subtract_drift(n_msdx_max)=.false.
    !> handles for msd objects for mean square tracer x displacement
    !> (\c -1 means: slot not in use)
    integer,save :: msd_tracer_x(n_msdx_max)=-1
    !> \}

    !> \{
    !> \name options for mean square y-displacement of tracers
    integer,save,public :: n_msdy=0 !< dump msdy how often? (0:never)
    integer,save,public :: n_msdy_sample=10 !< sample msdy how often
    !> number of data points in output (spacing between points is \c
    !> n_msdy_sample )
    integer,save,public :: msdy_len=100
    !> number of slots available to define sub-volumes for mean square
    !> y-displacement calculation
    integer,parameter,public :: n_msdy_max=10
    !> human readable name for a mean square y-displacement
    !> calculation region, a value of '' disables the slot
    character(len=16),save,public :: msdy_name(n_msdy_max)=''
    !> don't sample tracers at momentary x-positions below this value
    !> (negative value: automatically fill in \c 0.5 )
    real(kind=rk),save,public :: msdy_minx(n_msdy_max)=-1.0_rk
    !> don't sample tracers at momentary x-positions above this value
    !> (negative value: automatically fill in \c 0.5+tnx )
    real(kind=rk),save,public :: msdy_maxx(n_msdy_max)=-1.0_rk
    !> don't sample tracers at momentary y-positions below this value
    !> (negative value: automatically fill in \c 0.5 )
    real(kind=rk),save,public :: msdy_miny(n_msdy_max)=-1.0_rk
    !> don't sample tracers at momentary y-positions above this value
    !> (negative value: automatically fill in \c 0.5+tny )
    real(kind=rk),save,public :: msdy_maxy(n_msdy_max)=-1.0_rk
    !> don't sample tracers at momentary z-positions below this value
    !> (negative value: automatically fill in \c 0.5 )
    real(kind=rk),save,public :: msdy_minz(n_msdy_max)=-1.0_rk
    !> don't sample tracers at momentary z-positions above this value
    !> (negative value: automatically fill in \c 0.5+tnz )
    real(kind=rk),save,public :: msdy_maxz(n_msdy_max)=-1.0_rk
    !> correct for a possible drift of the whole system?
    logical,save,public :: msdy_subtract_drift(n_msdy_max)=.false.
    !> handles for msd objects for mean square tracer y displacement
    !> (\c -1 means: slot not in use)
    integer,save :: msd_tracer_y(n_msdy_max)=-1
    !> \}

    !> \{
    !> \name options for mean square z-displacement of tracers
    integer,save,public :: n_msdz=0 !< dump msdz how often? (0:never)
    integer,save,public :: n_msdz_sample=10 !< sample msdz how often
    !> number of data points in output (spacing between points is \c
    !> n_msdz_sample )
    integer,save,public :: msdz_len=100
    !> number of slots available to define sub-volumes for mean square
    !> x-displacement calculation
    integer,parameter,public :: n_msdz_max=10
    !> human readable name for a mean square z-displacement
    !> calculation region, a value of '' disables the slot
    character(len=16),save,public :: msdz_name(n_msdz_max)=''
    !> don't sample tracers at momentary x-positions below this value
    !> (negative value: automatically fill in \c 0.5 )
    real(kind=rk),save,public :: msdz_minx(n_msdz_max)=-1.0_rk
    !> don't sample tracers at momentary x-positions above this value
    !> (negative value: automatically fill in \c 0.5+tnx )
    real(kind=rk),save,public :: msdz_maxx(n_msdz_max)=-1.0_rk
    !> don't sample tracers at momentary y-positions below this value
    !> (negative value: automatically fill in \c 0.5 )
    real(kind=rk),save,public :: msdz_miny(n_msdz_max)=-1.0_rk
    !> don't sample tracers at momentary y-positions above this value
    !> (negative value: automatically fill in \c 0.5+tny )
    real(kind=rk),save,public :: msdz_maxy(n_msdz_max)=-1.0_rk
    !> don't sample tracers at momentary z-positions below this value
    !> (negative value: automatically fill in \c 0.5 )
    real(kind=rk),save,public :: msdz_minz(n_msdz_max)=-1.0_rk
    !> don't sample tracers at momentary z-positions above this value
    !> (negative value: automatically fill in \c 0.5+tnz )
    real(kind=rk),save,public :: msdz_maxz(n_msdz_max)=-1.0_rk
    !> correct for a possible drift of the whole system?
    logical,save,public :: msdz_subtract_drift(n_msdz_max)=.false.
    !> handles for msd objects for mean square tracer z displacement
    !> (\c -1 means: slot not in use)
    integer,save :: msd_tracer_z(n_msdz_max)=-1
    !> \}

    !> contains what is necessary for what was chosen to be dumped
    integer,save :: tracer_dump_mpitype

    !> recv buffer for dumping and checkpointing. Used only on rank
    !> zero, but with MPI_ALLGV_FASTER_THAN_GV allocation is required
    !> on all ranks
    type(tracer_particle_type),allocatable,dimension(:),save :: tp

    !> x-position of plane where flow rates are measured
    real(kind=rk),save :: x_flow_rate_plane
    !> last time step flow rates were calculated
    integer,save :: last_report_flow_rate=0
    !> number of particles of kind "1" crossing \c x_flow_rate_plane
    !> in positive direction since \c last_report_flow_rate
    integer,save :: flow_count_1=0
    !> number of particles of kind "2" crossing \c x_flow_rate_plane
    !> in positive direction since \c last_report_flow_rate
    integer,save :: flow_count_2=0

contains

    !> count tracer in case it just passed the plane at which the
    !> tracer flow is measured
    !>
    !> \param[in,out] t tracer
    !>
    !> \param[in,out] px previous tracer position
    subroutine count_tracer_flow(t,px)
        type(tracer_particle_type),intent(in) :: t
        real(kind=rk),intent(in) :: px(3)
        integer :: cnt

        if (passed_plane(t%x(1),px(1),x_flow_rate_plane)) then
           if (px(1)<x_flow_rate_plane) then
              cnt = 1         ! flow in positive x-direction
           else
              cnt = -1        ! flow in negative x-direction
           end if

           if (t%kind==1) then
              flow_count_1 = flow_count_1+cnt
           else
              flow_count_2 = flow_count_2+cnt ! assume this is kind "2"
           end if
        end if
    end subroutine count_tracer_flow

    !> initialization steps required to run during "input" phase
    !> before later initialization steps
    subroutine input_dump
        if (n_msdx/=0) call input_tracer_dump_msdx()
        if (n_msdy/=0) call input_tracer_dump_msdy()
	if (n_msdz/=0) call input_tracer_dump_msdz()
    end subroutine input_dump

    !> initialization of tracer mean square x displacement calculation
    !> and output
    !>
    !> This needs to be run at input stage so \c msd_init() can be
    !> called before \c tracer_init() where positions are sampled for
    !> the first time.
    subroutine input_tracer_dump_msdx()
        integer i

        do i=1,n_msdx_max
           if (len_trim(msdx_name(i))/=0) then
              call register_msd('tracer-x-'//trim(msdx_name(i))&
                   &,n_msdx,n_msdx_sample,msdx_len,msd_tracer_x(i)&
		   &,periodic=.true.,d_periodic=tsize(1)&
                   &,subtract_drift=iif(msdx_subtract_drift(i),1,-1))
           end if
        end do
    end subroutine input_tracer_dump_msdx

    !> initialization of tracer mean square y displacement calculation
    !> and output
    !>
    !> This needs to be run at input stage so \c msd_init() can be
    !> called before \c tracer_init() where positions are sampled for
    !> the first time.
    subroutine input_tracer_dump_msdy()
        integer i

        do i=1,n_msdy_max
           if (len_trim(msdy_name(i))/=0) then
              call register_msd('tracer-y-'//trim(msdy_name(i))&
                   &,n_msdy,n_msdy_sample,msdy_len,msd_tracer_y(i)&
                   &,periodic=.true.,d_periodic=tsize(2)&
                   &,subtract_drift=iif(msdy_subtract_drift(i),2,-1))
           end if
        end do
    end subroutine input_tracer_dump_msdy

    !> initialization of tracer mean square z displacement calculation
    !> and output
    !>
    !> This needs to be run at input stage so \c msd_init() can be
    !> called before \c tracer_init() where positions are sampled for
    !> the first time.
    subroutine input_tracer_dump_msdz()
        integer i

        do i=1,n_msdz_max
           if (len_trim(msdz_name(i))/=0) then
              call register_msd('tracer-z-'//trim(msdz_name(i))&
                   &,n_msdz,n_msdz_sample,msdz_len,msd_tracer_z(i)&
                   &,periodic=.true.,d_periodic=tsize(3)&
                   &,subtract_drift=iif(msdz_subtract_drift(i),3,-1))
           end if
        end do
    end subroutine input_tracer_dump_msdz

    !> sample tracer mean square x displacements
    subroutine sample_msdx()
        integer :: i,ii,j

        msdx: do j=1,n_msdx_max
           if (msd_tracer_x(j)<1) cycle

           call log_msg_tracer(&
                &'Sampling positions for tracer mean square x displacement '&
                &//'<tracer-x-'//trim(msdx_name(j))//'>')

           i = atompnt
           tracers: do ii = 1,nlocal
              if (all(T(i)%x>=(/msdx_minx(j),msdx_miny(j),msdx_minz(j)/)).and.&
                   &all(T(i)%x<(/msdx_maxx(j),msdx_maxy(j),msdx_maxz(j)/))) then
                 call sample(msd_tracer_x(j),T(i)%uid,T(i)%x(1))
              else
                 call close_trace(msd_tracer_x(j),T(i)%uid)
              end if

              i = list(i)
           end do tracers
        end do msdx
    end subroutine sample_msdx

    !> sample tracer mean square y displacements
    subroutine sample_msdy()
        integer :: i,ii,j

        msdy: do j=1,n_msdy_max
           if (msd_tracer_y(j)<1) cycle

           call log_msg_tracer(&
                &'Sampling p ositions for tracer mean square y displacement '&
                &//'<tracer-y-'//trim(msdy_name(j))//'>')

           i = atompnt
           tracers: do ii = 1,nlocal
              if (all(T(i)%x>=(/msdy_minx(j),msdy_miny(j),msdy_minz(j)/)).and.&
                   &all(T(i)%x<(/msdy_maxx(j),msdy_maxy(j),msdy_maxz(j)/))) then
                 call sample(msd_tracer_y(j),T(i)%uid,T(i)%x(2))
              else
                 call close_trace(msd_tracer_y(j),T(i)%uid)
              end if

              i = list(i)
           end do tracers
        end do msdy
    end subroutine sample_msdy

    !> sample tracer mean square z displacements
    subroutine sample_msdz()
        integer :: i,ii,j

        msdz: do j=1,n_msdz_max
           if (msd_tracer_x(j)<1) cycle

           call log_msg_tracer(&
                &'Sampling positions for tracer mean square z displacement '&
                &//'<tracer-z-'//trim(msdz_name(j))//'>')

           i = atompnt
           tracers: do ii = 1,nlocal
              if (all(T(i)%x>=(/msdz_minx(j),msdz_miny(j),msdz_minz(j)/)).and.&
                   &all(T(i)%x<(/msdz_maxx(j),msdz_maxy(j),msdz_maxz(j)/))) then
                 call sample(msd_tracer_x(j),T(i)%uid,T(i)%x(1))
              else
                 call close_trace(msd_tracer_x(j),T(i)%uid)
              end if

              i = list(i)
           end do tracers
        end do msdz
    end subroutine sample_msdz

    !> creates \c tracer_dump_mpitype and initializes tracer output
    !> routines
    subroutine setup_dump()
        call build_tracer_mpitype(tracer_dump_mpitype,x=dump_positions&
             &,v=dump_velocities,uid=dump_ids,akind=dump_kinds)

        call setup_tracer_flow_rate()
    end subroutine setup_dump

    !> initialization of tracer flow rate calculation and output
    subroutine setup_tracer_flow_rate
        ! obtain flow rate at center position in x direction
        x_flow_rate_plane = minpos(1)+0.5_rk*tsize(1)
    end subroutine setup_tracer_flow_rate

    !> writes tracer data to file
    subroutine dump_configuration
        select case (dump_format)
        case ('asc')
           call prepare_serial_dump
           call dump_asc
        case ('vtk')
           call prepare_serial_dump
           call dump_vtk
        case ('xdr')
           call prepare_serial_dump
           call dump_xdr
        case default
           call error_tracer('unknown value: dump_format="'//dump_format//'"')
        end select
    end subroutine dump_configuration

    !> calculate and report tracer flow rates
    subroutine report_tracer_flow_rate
        integer :: ierror,sample_steps
        integer :: fc(2),fcsum(2)
        real(kind=rk) :: fr(2)

        fc = (/flow_count_1,flow_count_2/)
        call MPI_Reduce(fc,fcsum,2,MPI_INTEGER,MPI_SUM,0,comm_cart,ierror)

        if (myrankc==0) then
           sample_steps = nt-last_report_flow_rate
           fr = real(fcsum,kind=rk)/real(sample_steps,kind=rk)

           write (msgstr,"('nt= ',I0,' flow rates during last ',I0"&
                &//",' time steps crossing plane x= ',F16.10,' : Q(kind=1)= '"&
                &//",F16.10,' Q(kind=2)= ',F16.10)") nt,sample_steps&
                &,x_flow_rate_plane,fr
           call log_msg_tracer(msgstr)

           last_report_flow_rate = nt
        end if

        flow_count_1 = 0
        flow_count_2 = 0
    end subroutine report_tracer_flow_rate

    !> calculate and dump tracer concentration profiles
    !>
    !> \param[in] N local lattice chunk including full halo
    subroutine dump_tracer_profile(N)
        type(lbe_site),intent(in) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)

        call log_msg_tracer('Dumping tracer profile...')

        ! for the moment dump only profile along x direction
        call dump_tracer_profile_dir(N,1)
    end subroutine dump_tracer_profile

    !> calculates and dumps a profile of tracer concentration, one kind in a row
    !>
    !> \param[in] N local lattice chunk including full halo
    !>
    !> \param[in] dir profile direction, 1,2,3 for x,y,z
    subroutine dump_tracer_profile_dir(N,dir)
        type(lbe_site),intent(in) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        integer,intent(in) :: dir
        integer,parameter :: fu=21,n_kinds=2
        character,parameter :: dirname(3)=(/'x','y','z'/)
        character(len=1024) :: fn
        integer,allocatable :: nft(:,:),nftsum(:,:)
        integer :: i,ierror,ii,tkind,pos,stat,ts(3),x,xx(3),y,z

        ts = (/tnx,tny,tnz/)

        ! allocate send and recv buffers
        allocate (nft(0:n_kinds,maxval(ts)),stat=stat)
        call check_allocate(stat,'dump_tracer_profile_dir(): nft')
        nft = 0
        if (myrankc==0) then
           allocate (nftsum(0:n_kinds,maxval(ts)),stat=stat)
           call check_allocate(stat,'dump_tracer_profile_dir(): nftsum')
        end if

        ! count tracers in each slice of the system sort them into nft
        ! according to their kind
        i = atompnt
        tracers: do ii = 1,nlocal
           tkind = T(i)%kind
           pos = int(T(i)%x(dir))
           if (pos<1) pos = pos + ts(dir)
           if (pos>ts(dir)) pos = pos - ts(dir)

           nft(tkind,pos) = nft(tkind,pos) + 1

           i = list(i)
        end do tracers

        ! count fluid sites in each slice of the system, put them into
        ! the 0th element of nft
        do x=1,nx
           do y=1,ny
              do z=1,nz
                 xx = (/x,y,z/)
                 pos = xx(dir)+start(dir)-1
                 if (is_fluid(N(x,y,z)%rock_state)) nft(0,pos) = nft(0,pos) + 1
              end do
           end do
        end do

        ! sum up everything on the root process
        call MPI_Reduce(nft,nftsum,(n_kinds+1)*ts(dir),MPI_INTEGER,MPI_SUM,0&
             &,comm_cart,ierror)

        ! write calculated profile to file
        if (myrankc==0) then
           call lbe_make_filename_output(fn,'tracer-profile-'//dirname(dir)&
                &,'.asc',nt)
           open (unit=fu,file=fn,status='REPLACE',action='WRITE',recl=80)
           write (unit=fu,advance='no',fmt='("# pos")')
           do ii=1,n_kinds
              write (unit=fu,advance='no',fmt='("         conc(",I1,")")') ii
           end do
           write (unit=fu,advance='yes',fmt='()')
           do i=1,ts(dir)
              write (unit=fu,advance='no',fmt='(I5)') i
              do ii=1,n_kinds
                 write (unit=fu,advance='no',fmt='(X,ES15.8)') &
                      &real(nftsum(ii,i),kind=rk)/real(nftsum(0,i),kind=rk)
              end do
              write (unit=fu,advance='yes',fmt='()')
           end do
           close (fu)
        end if

        ! clean up
        deallocate (nft)
        if (myrankc==0) deallocate (nftsum)
    end subroutine dump_tracer_profile_dir

    !> collects the required data on rank 0
    subroutine prepare_serial_dump
        call log_msg_tracer('Dumping tracer data...')
        call gather_tracers_provide_rbuf(tp,'prepare_serial_dump(): tp'&
             &,dump_mod_uid)
        call gather_tracers(tp,tracer_dump_mpitype,dump_mod_uid)
    end subroutine prepare_serial_dump

    !> writes prepared data into asc file on rank 0
    subroutine dump_asc
        integer,parameter :: cfg_file_unit=12
        character(len=1024) cfg_file_name
        integer i,n_global

        call dump_description_on_first_call

        call count_tracers(n_global,0,dump_mod_uid)
        rank0: if (myrankc==0) then
           call lbe_make_filename_output(cfg_file_name,'tracer-cfg','.asc',nt)
           open (unit=cfg_file_unit,file=cfg_file_name,status='REPLACE'&
                &,action='WRITE',recl=650)

           do i=1,n_global
              if (dump_positions) write (unit=cfg_file_unit&
                   &,fmt='(SP,3(ES15.8,X))',advance='no') tp(i)%x(:)
              if (dump_velocities) write (unit=cfg_file_unit&
                   &,fmt='(SP,3(ES15.8,X))',advance='no') tp(i)%v
              if (dump_kinds) write (unit=cfg_file_unit&
                   &,fmt='(SS,I10.10,X)',advance='no') tp(i)%kind
              if (dump_ids) write (unit=cfg_file_unit&
                   &,fmt='(SS,I10.10)',advance='no') tp(i)%uid
              write (unit=cfg_file_unit,fmt='()',advance='yes')
           end do
           close (cfg_file_unit)
        end if rank0
    end subroutine dump_asc

    !> writes prepared data into ascii vtk file on rank 0
    subroutine dump_vtk
        integer,parameter :: cfg_file_unit=12
        character(len=1024) cfg_file_name
        integer i,n_global

        call count_tracers(n_global,0,dump_mod_uid)
        rank0: if (myrankc==0) then
           if (.not.dump_positions) call log_msg_tracer(&
                &'ignoring dump_positions==.false.'&
                &//' - positions are mandatory for vtk output!')

           call lbe_make_filename_output(cfg_file_name,'tracer-cfg','.vtk',nt)
           open (unit=cfg_file_unit,file=cfg_file_name,status='REPLACE'&
                &,action='WRITE',recl=80)
           write (unit=cfg_file_unit,fmt='(A)') '# vtk DataFile Version 3.0'
           write (unit=cfg_file_unit,fmt='(A)') 'tracer cfg "'&
                &//trim(gr_out_file)//'"'
           write (unit=cfg_file_unit,fmt='(A)') 'ASCII'
           write (unit=cfg_file_unit,fmt='(A)') 'DATASET POLYDATA'
           write (unit=cfg_file_unit,fmt='(A,I0,A)') 'POINTS ',n_global,' FLOAT'
           do i=1,n_global
             write (unit=cfg_file_unit,fmt='(3(ES15.8,:,X))') tp(i)%x
           end do

           if (dump_velocities) then
              write (unit=cfg_file_unit,fmt='(A,I0)') 'POINT_DATA ',n_global
              write (unit=cfg_file_unit,fmt='(A)') 'VECTORS velocity FLOAT'
              do i=1,n_global
                 write (unit=cfg_file_unit,fmt='(SP,3(ES15.8,:,X))') tp(i)%v
              end do
           end if

           if (dump_kinds) then
              write (unit=cfg_file_unit,fmt='(A,I0)') 'POINT_DATA ',n_global
              write (unit=cfg_file_unit,fmt='(A)') 'SCALARS kind INT'
              do i=1,n_global
                 write (unit=cfg_file_unit,fmt='(SS,I10.10)') tp(i)%kind
              end do
           end if

           if (dump_ids) then
              write (unit=cfg_file_unit,fmt='(A,I0)') 'POINT_DATA ',n_global
              write (unit=cfg_file_unit,fmt='(A)') 'SCALARS uid INT'
              do i=1,n_global
                 write (unit=cfg_file_unit,fmt='(SS,I10.10)') tp(i)%uid
              end do
           end if

           close (cfg_file_unit)
        end if rank0
    end subroutine dump_vtk

    !> writes prepared data into xdr file on rank 0
    subroutine dump_xdr
        integer,parameter :: fk=4
        character(len=1024) cfg_file_name
        integer i,ierror,cfg_file_handle,n_global
        real(kind=rk) :: o(3)
        real(kind=fk) :: tmp(0:3)
        real(kind=rk) :: tmp2(0:3)

        ! These would be function pointers in C. In fact,  xdrf(double|float)
        ! ares external functions, but fortran does not seem to care about the
        ! type anyway...
        integer,external :: xdrffloat, xdrfdouble

        call dump_description_on_first_call

        call count_tracers(n_global,0,dump_mod_uid)
        rank0: if (myrankc==0) then
#ifdef USEXDRF
           call lbe_make_filename_output(cfg_file_name,'tracer-cfg','.xdr',nt)
           call xdrfopen(cfg_file_handle,cfg_file_name,'w',ierror)
           if (ierror==0) then
              call log_msg_tracer('dump_xdr(): xdrfopen() failed')
              call check_xdrfsloppy()
           end if

           double: if (dump_double) then
              do i=1,n_global
                 if (dump_positions) call xdrfvector&
                      &(cfg_file_handle,tp(i)%x(1:3),3,xdrfdouble,ierror)
                 if (dump_velocities) call xdrfvector&
                      &(cfg_file_handle,tp(i)%v(1:3),3,xdrfdouble,ierror)
                 if (dump_kinds) call xdrfint(cfg_file_handle,tp(i)%kind,ierror)
                 if (dump_ids) call xdrfint(cfg_file_handle,tp(i)%uid,ierror)
              end do
           else double
              do i=1,n_global
                 if (dump_positions) then
                    tmp(1:3) = real(tp(i)%x(1:3),kind=fk)
                    call xdrfvector(cfg_file_handle,tmp(1:3),3,xdrffloat,ierror)
                 end if
                 if (dump_velocities) then
                    tmp(1:3) = real(tp(i)%v(1:3),kind=fk)
                    call xdrfvector(cfg_file_handle,tmp(1:3),3,xdrffloat,ierror)
                 end if
                 if (dump_kinds) call xdrfint(cfg_file_handle,tp(i)%kind,ierror)
                 if (dump_ids) call xdrfint(cfg_file_handle,tp(i)%uid,ierror)
              end do
           end if double

           call xdrfclose(cfg_file_handle,ierror)
#else
           call log_msg_tracer('failed to dump '//cfg_file_name&
                &//' - recompile with USEXDRF set!')
#endif
        end if rank0
    end subroutine dump_xdr

    !> write a checkpoint of the current state
    subroutine tracer_dump_checkpoint
	character(len=1024) chk_file_name
        integer chk_file_handle,i,ierror,n_global

        ! This would be a function pointer in C. In fact,  xdrfdouble  is
        ! an external function, but fortran does not seem to care about the
        ! type anyway...
        integer,external :: xdrfdouble

        call log_msg_tracer("  Dumping TRACER checkpoint ...")

        if ( checkpoint_format .eq. 'xdr') then
#ifdef USEXDRF
           if (myrankc==0) then
              call lbe_make_filename_cp(chk_file_name,'tracer-checkpoint'&
                   &,'.xdr',nt)
              call xdrfopen(chk_file_handle,chk_file_name,'w',ierror)
              if (ierror==0) then
                 call log_msg_tracer(&
                      &'tracer_dump_checkpoint(): xdrfopen() failed')
                 call check_xdrfsloppy()
              end if
           end if

           ! ...at last, collect and dump the configuration of all particles:
           call gather_tracers_provide_rbuf(tp,'tracer_dump_checkpoint(): tp')
           call gather_tracers(tp)

           call count_tracers(n_global,0)
           if (myrankc==0) then
              call xdrfint(chk_file_handle,n_global,ierror)

              do i=1,n_global
                 call xdrfvector(chk_file_handle,tp(i)%x(1:3),3,xdrfdouble&
                      &,ierror)
                 call xdrfvector(chk_file_handle,tp(i)%v(1:3),3,xdrfdouble&
                      &,ierror)
                 call xdrfint(chk_file_handle,tp(i)%uid,ierror)
                 call xdrfint(chk_file_handle,tp(i)%kind,ierror)
              end do
           end if

           call tracer_flow_rate_dump_checkpoint_xdr(chk_file_handle)

           if (myrankc==0) call xdrfclose(chk_file_handle,ierror)
#else
           ! ifndef USEXDRF
           call log_msg_tracer('failed to dump '//chk_file_name&
                &//' - recompile with USEXDRF set!')
#endif

        else if ( checkpoint_format .eq. 'hdf' ) then
#ifdef USEHDF
           call log_msg_tracer(&
                &'WARNING: TRACER HDF5 checkpointing not yet implemented!')
#else
           ! ifndef USEHDF
           call log_msg_tracer('WARNING: failed to dump TRACER checkpoint '&
                &//'- recompile with USEHDF set!')
#endif
        endif

        call log_msg_tracer("  Finished dumping TRACER checkpoint.")

    end subroutine tracer_dump_checkpoint

    !> write tracer flow rate-related data to an open \c TRACER
    !> checkpoint file
    !>
    !> \param[in,out] h \c xdrf file handle
    subroutine tracer_flow_rate_dump_checkpoint_xdr(h)
        integer,intent(inout) :: h
        integer :: fc(2),fcsum(2),ierror

        fc = (/flow_count_1,flow_count_2/)
        call MPI_Reduce(fc,fcsum,2,MPI_INTEGER,MPI_SUM,0,comm_cart,ierror)

        if (myrankc==0) then
           call xdrfint(h,last_report_flow_rate,ierror)
           call xdrfint(h,fcsum(1),ierror)
           call xdrfint(h,fcsum(2),ierror)
        end if
    end subroutine tracer_flow_rate_dump_checkpoint_xdr

    !> read tracer flow rate-related data from an open \c TRACER
    !> checkpoint file
    !>
    !> \param[in,out] h \c xdrf file handle
    subroutine tracer_flow_rate_restore_checkpoint_xdr(h)
        integer,intent(inout) :: h
        integer :: ierror

        ! There is no need to broadcast or distribute flow counts
        ! since afterwards they will just be summed globally.
        if (myrankc==0) then
           call xdrfint(h,last_report_flow_rate,ierror)
           call xdrfint(h,flow_count_1,ierror)
           call xdrfint(h,flow_count_2,ierror)
        else
           flow_count_1 = 0
           flow_count_2 = 0
        end if
    end subroutine tracer_flow_rate_restore_checkpoint_xdr

    !> This routine deletes old TRACER checkpoint files.
    subroutine tracer_delete_checkpoint(last)
        integer, intent(in)  :: last
        character(len=1024) :: filename

        call log_msg_tracer("  Deleting TRACER checkpoint files ...")

        if (myrankc==0) then
           if ( checkpoint_format .eq. 'xdr') then
              call lbe_make_filename_cp(filename,'tracer-checkpoint','.xdr'&
                   &,last)
              call lbe_delete_file(filename)
           else if (checkpoint_format .eq. 'hdf' ) then
              call lbe_make_filename_cp(filename, 'cp_tracer_', '.h5', last)
              call lbe_delete_file(filename)
           end if
        end if
    end subroutine tracer_delete_checkpoint

    subroutine dump_description_on_first_call
        integer,parameter :: cfgdesc_file_unit=12
        logical,save :: first_call=.true.
        character(len=1024) cfgdesc_file_name

        if (.not.first_call) return
        first_call = .false.

        rank0: if (myrankc==0) then
           call lbe_make_filename_output(cfgdesc_file_name,'tracer-cfg-desc'&
                &,'.txt',nt)
           open (unit=cfgdesc_file_unit,file=cfgdesc_file_name&
                &,status='REPLACE',action='WRITE',recl=650)

           if (dump_positions) write (unit=cfgdesc_file_unit&
                &,fmt='(3(A15,X))',advance='no') 'x','y','z'
           if (dump_velocities) write (unit=cfgdesc_file_unit&
                &,fmt='(3(A15,X))',advance='no') 'vx','vy','vz'
           if (dump_kinds) write (unit=cfgdesc_file_unit&
                &,fmt='(A10,X)',advance='no') 'kind'
           if (dump_ids) write (unit=cfgdesc_file_unit&
                &,fmt='(A10)',advance='no') 'uid'

           write (unit=cfgdesc_file_unit,fmt='()',advance='yes')
           close (cfgdesc_file_unit)
        end if rank0
    end subroutine dump_description_on_first_call

#endif
end module lbe_tracer_output_module
