#include "lbe.h"

!> output for md part of the code
module lbe_md_output_module
#ifdef MD
    use lbe_globals_module, only: cz,g,halo_extent,maxpos,minpos,n_spec,tsize,myrankc
    use lbe_helper_module, only: iif,is_restoring
    use lbe_md_dynamic_module, only: first_inserted_puid
    use lbe_md_fluid_module, only: summary_fluid
#ifdef PARTICLESTRESS
    use lbe_md_fluid_ladd_module, only: A_cut
#endif
    use lbe_md_fluid_ladd_mc_module, only: pfr, pfb, pfg, global_mass_change, global_mass_target
    use lbe_md_globals_module
    use lbe_md_helper_module, only: build_particle_mpitype,count_particles&
         &,count_particles_all,error_md,log_msg_md,md_make_filename_particles&
         &,orientation,particle_radii,space_to_body_matrix,log_msg_md_hdr&
         &,warning
    use lbe_md_parallel_module, only: gather_particles&
         &,gather_particles_provide_rbuf
    use lbe_mean_square_displacement_mod, only: close_trace,register_msd,sample
    use lbe_timer_module, only: sum_timer
    use lbe_parallel_module, only: calculate_displacements,check_allocate,cdims&
         &,comm_cart,nprocs,start,stats_rk,tnx,tny,tnz
    use lbe_parms_module, only: &
#ifndef SINGLEFLUID
         &amass_b,amass_r,&
#ifndef NOSURFACTANT
         &amass_s,&
#endif
#endif
         &boundary,n_checkpoint,chk_uid,gr_out_file&
         &,n_iteration,nt,num_chkp_files,n_restore,checkpoint_format&
         &,checkpoint_safe,nx,ny,nz
    use lbe_types_module, only: lbe_site
    use lbe_io_helper_module, only: lbe_delete_file,lbe_make_filename_cp&
         &,lbe_make_filename_output,lbe_make_filename_output_substep

#ifdef USEXDRF
  use lbe_io_xdrf_module, only: check_xdrfsloppy
#endif

    implicit none
    private

    include 'mpif.h'

    public dump_configuration,dump_rvz_profile,input_dump,md_dump_checkpoint&
         &,md_delete_checkpoint,sample_msdx,sample_msdy,sample_msdz&
	 &,setup_dump,summary

    !> dump cfg interval in lbe steps (0: 'never')
    integer,save,public :: n_dump=100
    !> dump particle configurations each that many substeps (0: no
    !> substep dumping)
    integer,save,public :: n_dump_substep=0
    !> disable substep dumping before this LB time step
    integer,save,public :: n_dump_substep_start=0

    character(len=32),save,public :: dump_format='asc' !< dump file format
    logical,save,public :: dump_double=.false. !< dump doubles or floats?

    logical,save,public :: dump_time=.false. !< include \c t into output?
    logical,save,public :: dump_positions=.true. !< include \c x into output?
    logical,save,public :: dump_pnd=.false. !< include \c PND into output?
    logical,save,public :: dump_velocities=.true. !< include \c v into output?
    logical,save,public :: dump_forces=.false.     !< include \c f into output?
#ifdef FORCECOMPONENT
    logical,save,public :: dump_force_components=.true. !< include \c f_normal and f_tangent
#endif    
    logical,save,public :: dump_quaternions=.false. !< include \c q into output?

    !> include unit vector parallel to symmetry axis into output?
    logical,save,public :: dump_orientations=.true.

    logical,save,public :: dump_rotations=.true. !< include \c w into output?
    logical,save,public :: dump_torques=.false. !< include \c t into output?
    !> include \c R_orth and R_para into output?
    logical,save,public :: dump_radii=.false.
    !> include interaction between magnetic particles into output? 
    logical,save,public :: dump_mag=.false.
#ifdef PARTICLESTRESS
    !> include averaged body fixed stresses into output?
    logical,save,public :: dump_stresses_b=.false.
    !> include averaged space fixed stresses into output?
    logical,save,public :: dump_stresses=.false.
#endif
    !> include potential energy \c e_pot into output?
    logical,save,public :: dump_potentials=.false.
    logical,save,public :: dump_ids=.true. !< include \c uid into output?
#ifdef LADD_SURR_RHOF
    logical,save,public :: dump_surr_rhof=.false. !< dump \c rhof and \c n_acc?
#endif

    logical,save,public :: dump_summary=.true. !< dump a summary at the end?

    !> \{
    !> \name radial profile output
    !>
    !> sample \f$v_z(r)\f$ every \c n_rvz_profile_sample timesteps (0
    !> means never) and dump it every \c n_rvz_profile_dump timesteps.
    !> \c n_rvz_profile_dump should be an integer multiple of
    !> \c n_rvz_profile_sample!
    integer,save,public :: n_rvz_profile_dump=1
    integer,save,public :: n_rvz_profile_sample=0
    !> \}

    !> \{
    !> \name options for mean square x-displacement of md particles
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
    !> don't sample MD particles before this time step
    integer(kind=rk),save,public :: msdx_start(n_msdx_max)=-1
    !> don't sample MD particles at momentary x-positions below this value
    !> (negative value: automatically fill in \c 0.5 )
    real(kind=rk),save,public :: msdx_minx(n_msdx_max)=-1.0_rk
    !> don't sample MD particles at momentary x-positions above this value
    !> (negative value: automatically fill in \c 0.5+tnx )
    real(kind=rk),save,public :: msdx_maxx(n_msdx_max)=-1.0_rk
    !> don't sample MD particles at momentary y-positions below this value
    !> (negative value: automatically fill in \c 0.5 )
    real(kind=rk),save,public :: msdx_miny(n_msdx_max)=-1.0_rk
    !> don't sample MD particles at momentary y-positions above this value
    !> (negative value: automatically fill in \c 0.5+tny )
    real(kind=rk),save,public :: msdx_maxy(n_msdx_max)=-1.0_rk
    !> don't sample MD particles at momentary z-positions below this value
    !> (negative value: automatically fill in \c 0.5 )
    real(kind=rk),save,public :: msdx_minz(n_msdx_max)=-1.0_rk
    !> don't sample MD particles at momentary z-positions above this value
    !> (negative value: automatically fill in \c 0.5+tnz )
    real(kind=rk),save,public :: msdx_maxz(n_msdx_max)=-1.0_rk
    !> correct for a possible drift of the whole system?
    logical,save,public :: msdx_subtract_drift(n_msdx_max)=.false.
    !> handle for msd object for mean square particle x displacement (-1
    !> means: slot not in use)
    integer,save :: msd_md_x(n_msdx_max)=-1
    !> \}

    !> \{
    !> \name options for mean square y-displacement of md particles
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
    !> don't sample MD particles before this time step
    integer(kind=rk),save,public :: msdy_start(n_msdy_max)=-1
    !> don't sample MD particles at momentary x-positions below this value
    !> (negative value: automatically fill in \c 0.5 )
    real(kind=rk),save,public :: msdy_minx(n_msdy_max)=-1.0_rk
    !> don't sample MD particles at momentary x-positions above this value
    !> (negative value: automatically fill in \c 0.5+tnx )
    real(kind=rk),save,public :: msdy_maxx(n_msdy_max)=-1.0_rk
    !> don't sample MD particles at momentary y-positions below this
    !> value (negative value: automatically fill in \c 0.5
    !> ). Independently from the value, for the msd calculation, the
    !> system is always expected to be periodic in y-direction.
    real(kind=rk),save,public :: msdy_miny(n_msdy_max)=-1.0_rk
    !> don't sample MD particles at momentary y-positions above this
    !> value (negative value: automatically fill in \c 0.5+tny
    !> ). Independently from the value, for the msd calculation, the
    !> system is always expected to be periodic in y-direction.
    real(kind=rk),save,public :: msdy_maxy(n_msdy_max)=-1.0_rk
    !> don't sample MD particles at momentary z-positions below this value
    !> (negative value: automatically fill in \c 0.5 )
    real(kind=rk),save,public :: msdy_minz(n_msdy_max)=-1.0_rk
    !> don't sample MD particles at momentary z-positions above this value
    !> (negative value: automatically fill in \c 0.5+tnz )
    real(kind=rk),save,public :: msdy_maxz(n_msdy_max)=-1.0_rk
    !> correct for a possible drift of the whole system?
    logical,save,public :: msdy_subtract_drift(n_msdy_max)=.false.
    !> handle for msd object for mean square particle y displacement (-1
    !> means: slot not in use)
    integer,save :: msd_md_y(n_msdy_max)=-1
    !> \}

    !> \{
    !> \name options for mean square z-displacement of md particles
    integer,save,public :: n_msdz=0 !< dump msdz how often? (0:never)
    integer,save,public :: n_msdz_sample=10 !< sample msdz how often
    !> number of data points in output (spacing between points is \c
    !> n_msdz_sample )
    integer,save,public :: msdz_len=100
    !> number of slots available to define sub-volumes for mean square
    !> z-displacement calculation
    integer,parameter,public :: n_msdz_max=10
    !> human readable name for a mean square z-displacement
    !> calculation region, a value of '' disables the slot
    character(len=16),save,public :: msdz_name(n_msdz_max)=''
    !> don't sample MD particles before this time step
    integer(kind=rk),save,public :: msdz_start(n_msdz_max)=-1
    !> don't sample MD particles at momentary z-positions below this value
    !> (negative value: automatically fill in \c 0.5 )
    real(kind=rk),save,public :: msdz_minx(n_msdz_max)=-1.0_rk
    !> don't sample MD particles at momentary z-positions above this value
    !> (negative value: automatically fill in \c 0.5+tnx )
    real(kind=rk),save,public :: msdz_maxx(n_msdz_max)=-1.0_rk
    !> don't sample MD particles at momentary z-positions below this value
    !> (negative value: automatically fill in \c 0.5 )
    real(kind=rk),save,public :: msdz_miny(n_msdz_max)=-1.0_rk
    !> don't sample MD particles at momentary z-positions above this value
    !> (negative value: automatically fill in \c 0.5+tny )
    real(kind=rk),save,public :: msdz_maxy(n_msdz_max)=-1.0_rk
    !> don't sample MD particles at momentary z-positions below this value
    !> (negative value: automatically fill in \c 0.5 )
    real(kind=rk),save,public :: msdz_minz(n_msdz_max)=-1.0_rk
    !> don't sample MD particles at momentary z-positions above this value
    !> (negative value: automatically fill in \c 0.5+tnz )
    real(kind=rk),save,public :: msdz_maxz(n_msdz_max)=-1.0_rk
    !> correct for a possible drift of the whole system?
    logical,save,public :: msdz_subtract_drift(n_msdz_max)=.false.
    !> handle for msd object for mean square particle x displacement (-1
    !> means: slot not in use)
    integer,save :: msd_md_z(n_msdz_max)=-1
    !> \}

    !> contains what is necessary for what was chosen to be dumped
    integer,save :: particle_dump_mpitype

    !> recv buffer for dumping and checkpointing. Used only on rank
    !> zero, but with MPI_ALLGV_FASTER_THAN_GV allocation is required
    !> on all ranks
    type(md_particle_type),allocatable,dimension(:),save :: tp

contains

    !> initialization steps required to run during "input" phase
    !> before later initialization steps
    subroutine input_dump
        if (dump_potentials) collect_forces = .true.
        if (n_msdx/=0) call input_md_dump_msdx()
        if (n_msdy/=0) call input_md_dump_msdy()
	if (n_msdz/=0) call input_md_dump_msdz()
    end subroutine input_dump

    !> initialization of MD particle mean square x displacement calculation
    !> and output
    !>
    !> This needs to be run at input stage so \c msd_init() can be
    !> called before \c md_init() where positions are sampled for
    !> the first time.
    subroutine input_md_dump_msdx()
        integer i

        do i=1,n_msdx_max
           if (len_trim(msdx_name(i))/=0) then
              call register_msd('md-x-'//trim(msdx_name(i))&
                   &,n_msdx,n_msdx_sample,msdx_len,msd_md_x(i)&
		   &,periodic=.true.,d_periodic=tsize(1)&
                   &,subtract_drift=iif(msdx_subtract_drift(i),1,-1))
           end if
        end do
    end subroutine input_md_dump_msdx

    !> initialization of MD particle mean square y displacement calculation
    !> and output
    !>
    !> This needs to be run at input stage so \c msd_init() can be
    !> called before \c md_init() where positions are sampled for
    !> the first time.
    subroutine input_md_dump_msdy()
        integer i

        do i=1,n_msdy_max
           if (len_trim(msdy_name(i))/=0) then
              call register_msd('md-y-'//trim(msdy_name(i))&
                   &,n_msdy,n_msdy_sample,msdy_len,msd_md_y(i)&
                   &,periodic=.true.,d_periodic=tsize(2)&
                   &,subtract_drift=iif(msdy_subtract_drift(i),2,-1))
           end if
        end do
    end subroutine input_md_dump_msdy

    !> initialization of MD particle mean square z displacement calculation
    !> and output
    !>
    !> This needs to be run at input stage so \c msd_init() can be
    !> called before \c md_init() where positions are sampled for
    !> the first time.
    subroutine input_md_dump_msdz()
        integer i

        do i=1,n_msdz_max
           if (len_trim(msdz_name(i))/=0) then
              call register_msd('md-z-'//trim(msdz_name(i))&
                   &,n_msdz,n_msdz_sample,msdz_len,msd_md_z(i)&
                   &,periodic=.true.,d_periodic=tsize(3)&
                   &,subtract_drift=iif(msdz_subtract_drift(i),3,-1))
           end if
        end do
    end subroutine input_md_dump_msdz

    !> sample MD particle mean square x displacements
    subroutine sample_msdx()
        integer :: i,ii,j

        msdx: do j=1,n_msdx_max
           if (msd_md_x(j)<1) cycle
           if (msdx_start(j)>nt) cycle

           call log_msg_md(&
                &'Sampling positions for particle mean square x displacement '&
                &//'<md-x-'//trim(msdx_name(j))//'>')

           i = atompnt
           particles: do ii = 1,nlocal
              if (all(P(i)%x>=(/msdx_minx(j),msdx_miny(j),msdx_minz(j)/)).and.&
                   &all(P(i)%x<(/msdx_maxx(j),msdx_maxy(j),msdx_maxz(j)/))) then
                 call sample(msd_md_x(j),P(i)%uid,P(i)%x(1))
              else
                 call close_trace(msd_md_x(j),P(i)%uid)
              end if

              i = list(i)
           end do particles
        end do msdx
    end subroutine sample_msdx

    !> sample MD particle mean square y displacements
    subroutine sample_msdy()
        integer :: i,ii,j

        msdy: do j=1,n_msdy_max
           if (msd_md_y(j)<1) cycle
           if (msdy_start(j)>nt) cycle

           call log_msg_md(&
                &'Sampling positions for particle mean square y displacement '&
                &//'<md-y-'//trim(msdy_name(j))//'>')

           i = atompnt
           particles: do ii = 1,nlocal
              if (all(P(i)%x>=(/msdy_minx(j),msdy_miny(j),msdy_minz(j)/)).and.&
                   &all(P(i)%x<(/msdy_maxx(j),msdy_maxy(j),msdy_maxz(j)/))) then
                 call sample(msd_md_y(j),P(i)%uid,P(i)%x(2))
              else
                 call close_trace(msd_md_y(j),P(i)%uid)
              end if

              i = list(i)
           end do particles
        end do msdy
    end subroutine sample_msdy

    !> sample MD particle mean square z displacements
    subroutine sample_msdz()
        integer :: i,ii,j

        msdz: do j=1,n_msdz_max
           if (msd_md_z(j)<1) cycle
           if (msdz_start(j)>nt) cycle

           call log_msg_md(&
                &'Sampling positions for particle mean square z displacement '&
                &//'<md-z-'//trim(msdy_name(j))//'>')

           i = atompnt
           particles: do ii = 1,nlocal
              if (all(P(i)%x>=(/msdz_minx(j),msdz_miny(j),msdz_minz(j)/)).and.&
                   &all(P(i)%x<(/msdz_maxx(j),msdy_maxz(j),msdz_maxz(j)/))) then
                 call sample(msd_md_z(j),P(i)%uid,P(i)%x(3))
              else
                 call close_trace(msd_md_z(j),P(i)%uid)
              end if

              i = list(i)
           end do particles
        end do msdz
    end subroutine sample_msdz

    !> creates \c particle_dump_mpitype , check input parameters
    subroutine setup_dump
        ! with md_leesedwards, every type that is used with
        ! gather_particles() must contain x
        call build_particle_mpitype(particle_dump_mpitype&
             &,x=dump_positions.or.dump_pnd.or.md_leesedwards,v=dump_velocities&
             &,f=dump_forces,w=dump_rotations,t=dump_torques&
#ifdef FORCECOMPONENT	     
	     &,f_normal=dump_force_components&
	     &,f_tangent=dump_force_components&
	     &,f_n=dump_force_components&
	     &,f_t=dump_force_components&
#endif
             &,uid=dump_ids&
             &,e_pot=dump_potentials,R=polydispersity.and.dump_radii&
             &,q=dump_quaternions.or.dump_orientations& ! watch out, this one is
                                                        ! continued...
#ifdef PARTICLESTRESS
             &.or.dump_stresses&
             &,tau=dump_stresses.or.dump_stresses_b&
#endif
#ifdef RWALK
             &,v_r=dump_velocities&
#endif
#ifdef LADD_SURR_RHOF
             &,rhof=dump_surr_rhof,rhofn_acc=dump_surr_rhof&
#endif
!#ifdef   MD_MAG
             &,ftmi=dump_mag&
!#endif
             )
#ifdef PARTICLESTRESS
        if ((dump_stresses.or.dump_stresses_b).and.dump_format=='vtk') &
             &call error_md('stress tensor not available for vtk md-cfg output'&
             &//'---disable dump_stresses and dump_stresses_b or set '&
             &//'dump_format to ''asc'' or ''vtk''!')
#endif
    end subroutine setup_dump

    !> writes positions and velocities of all particles to file
    !>
    !> \param[in] substep indicates that the dump belongs to a substep
    !> with the given number (optional, defaults to main step dump)
    !>
    !> \param[in] prefix optional, if specified, prepend the file name
    !> with this string, followed by a '_'
    !>
    !> \param[in] t optional, if specified, \c t is used as time step
    !> to build the file name, otherwise \c nt is taken
    subroutine dump_configuration(substep,prefix,t)
        integer,intent(in),optional :: substep
        character(len=*),intent(in),optional :: prefix
        integer,intent(in),optional :: t
        character(len=64) :: pprefix
        integer :: tt

        if (present(prefix)) then
           pprefix = trim(prefix)//'_md-cfg'
        else
           pprefix = 'md-cfg'
        end if

        if (present(t)) then
           tt = t
        else
           tt = nt
        end if

        select case (dump_format)
        case ('asc')
           call prepare_serial_dump
           call dump_asc(pprefix,tt,substep)
        case ('vtk')
           call prepare_serial_dump
                 write(*,"('inside dump_vtk')")
           call dump_vtk(pprefix,tt,substep)
        case ('xdr')
           call prepare_serial_dump
                 write(*,"('inside dump_xdr')")
           call dump_xdr(pprefix,tt,substep)
        case default
           call error_md('unknown value: dump_format="'//dump_format//'"')
        end select
    end subroutine dump_configuration

    !> collects the required data on rank 0
    subroutine prepare_serial_dump
#ifdef PARTICLESTRESS
        integer i,k,n_global
#endif

        call log_msg_md('Dumping particle data...')

        call gather_particles_provide_rbuf(tp,'prepare_serial_dump(): tp')
        call gather_particles(tp,particle_dump_mpitype)
#ifdef PARTICLESTRESS
        call count_particles(n_global,0)
        if (myrankc==0) then
           if (dump_stresses.or.dump_stresses_b) then
              do i=1,n_global
                 do k=1,3
                    ! The 0.5 is because tau contains the summed up
                    ! force from BOTH sides of each particle but for
                    ! example pressure is defined as force per ONE
                    ! surface A_cut.
                    tp(i)%tau(k,:) = 0.5_rk*tp(i)%tau(k,:)/A_cut(k)
                 end do
              end do
           end if
        end if
#endif
    end subroutine prepare_serial_dump

    !> writes prepared data into asc file on rank 0
    !>
    !> \param[in] prefix prefix used to build the file name
    !>
    !> \param[in] t time step used to build the file name and to write
    !> in case of \c dump_time
    !>
    !> \param[in] substep indicates that the dump belongs to a substep
    !> with the given number (optional, defaults to main step dump)
    subroutine dump_asc(prefix,t,substep)
        character(len=*),intent(in) :: prefix
        integer :: t
        integer,intent(in),optional :: substep
        integer i,n_global,Nx,Ny,Nz,x,y,z,j,k,l,Indx
        character(len=1024) cfg_file_name
        integer,parameter :: cfg_file_unit=12
        integer, allocatable :: B_n(:,:,:) 
        real(kind=rk) :: o(3)
        logical :: fexist
!#ifdef PARTICLESTRESS
        real(kind=rk) :: A(3,3)
!#endif

        call dump_description_on_first_call

        if (dump_pnd) then
           Nx = tnx
           Ny = tny
           Nz = tnz
           allocate (B_n(Nx,Ny,Nz))
           B_n=0
        end if

        call count_particles(n_global,0)

        rank0: if (myrankc==0) then
           if (time_output) then
              if (present(substep)) then
                 call lbe_make_filename_output_substep(cfg_file_name,prefix&
                      &,'.asc',t,substep)
              else
                 call lbe_make_filename_output(cfg_file_name,prefix,'.asc',t)
              end if
              open (unit=cfg_file_unit,file=cfg_file_name,status='REPLACE'&
                   &,action='WRITE',recl=650)

              do i=1,n_global
                if (dump_time) write (unit=cfg_file_unit,&
                     &fmt='(I10.10,1x)',advance='no') t

                if (dump_positions) write (unit=cfg_file_unit&
                     &,fmt='(SP,3(ES15.8,X))',advance='no') tp(i)%x(:)
               if(dump_pnd) then
                      x = nint(tp(i)%x(1))
                      y = nint(tp(i)%x(2))
                      z = nint(tp(i)%x(3))
                  B_n(x,y,z) = B_n(x,y,z) + 1

               endif

#ifndef RWALK
                 if (dump_velocities) write (unit=cfg_file_unit&
                      &,fmt='(SP,3(ES15.8,X))',advance='no') tp(i)%v
#else
                 tp(i)%v =tp(i)%v+tp(i)%v_r
                 if (dump_velocities) write (unit=cfg_file_unit&
                      &,fmt='(SP,3(ES15.8,X))',advance='no') tp(i)%v
#endif
                 if (dump_forces) write (unit=cfg_file_unit&
                      &,fmt='(SP,3(ES15.8,X))',advance='no') tp(i)%f
                 if (dump_quaternions) write (unit=cfg_file_unit&
                      &,fmt='(SP,4(ES15.8,X))',advance='no') tp(i)%q
                 if (dump_orientations) then
                    o = orientation(tp(i)%q)
                    write (unit=cfg_file_unit,fmt='(SP,3(ES15.8,X))'&
                         &,advance='no') o
                 end if
!> dump space-fixed angular velocity. Do conversion from body-fixed to space fixed.
!>                 if (dump_rotations) write (unit=cfg_file_unit&
!>                      &,fmt='(SP,3(ES15.8,X))',advance='no') tp(i)%w
                 if (dump_rotations) then
                 A = space_to_body_matrix(tp(i)%q)
                  write (unit=cfg_file_unit&
                       &,fmt='(SP,3(ES15.8,X))',advance='no') &
                  &matmul(transpose(A),tp(i)%w)                 
                 end if
!>
                 if (dump_torques) write (unit=cfg_file_unit&
                      &,fmt='(SP,3(ES15.8,X))',advance='no') tp(i)%t
                   
                 if (dump_mag) then
                    write (unit=cfg_file_unit,fmt='(SP,3(ES15.8,X))',advance='no') tp(i)%fmi
                    write (unit=cfg_file_unit,fmt='(SP,3(ES15.8,X))',advance='no') tp(i)%tmi
                end if

#ifdef PARTICLESTRESS
                 if (dump_stresses) then
                    A = space_to_body_matrix(tp(i)%q)
                    write (unit=cfg_file_unit,fmt='(SP,9(ES15.8,X))'&
                         &,advance='no') &
                         &matmul(matmul(transpose(A),tp(i)%tau),A)
                 end if
                 if (dump_stresses_b) write (unit=cfg_file_unit&
                      &,fmt='(SP,9(ES15.8,X))',advance='no') tp(i)%tau
#endif
                 if (dump_potentials) write (unit=cfg_file_unit&
                      &,fmt='(SP,ES15.8,X)',advance='no') tp(i)%e_pot

#ifdef LADD_SURR_RHOF
                 if (dump_surr_rhof) then
                    write (unit=cfg_file_unit&
                         &,fmt='(SP,ES15.8,X)',advance='no') tp(i)%rhof
                    write (unit=cfg_file_unit&
                         &,fmt='(SS,I10.10,X)',advance='no') tp(i)%n_acc
                 end if

#endif
                 if (dump_radii) write (unit=cfg_file_unit&
                      &,fmt='(SP,2(ES15.8,X))',advance='no') &
                      &particle_radii(tp(i))
                 if (dump_ids) write (unit=cfg_file_unit&
                      &,fmt='(SS,I10.10)',advance='no') tp(i)%uid
                 write (unit=cfg_file_unit,fmt='()',advance='yes')
              end do
              close (cfg_file_unit)
           end if

           if (prtcl_output.and.present(substep)) then
              call warning(&
                   &'substep dumping not implemented yet for prtcl_output')
           else if (prtcl_output) then
              do i=1,n_global
                 call md_make_filename_particles(cfg_file_name,prefix,'.asc'&
                      &,tp(i)%uid)

                 inquire(file=cfg_file_name,exist=fexist)
                 if (fexist) then
                    open(unit=cfg_file_unit,file=cfg_file_name,status='OLD'&
                         &,position='APPEND',recl=650)
                 else
                    open(unit=cfg_file_unit,file=cfg_file_name,status='NEW'&
                         &,position='APPEND',recl=650)
                 endif

                 if (dump_time) write (unit=cfg_file_unit,&
                      &fmt='(I10.10,1x)',advance='no') t

                 if (dump_positions) write (unit=cfg_file_unit&
                      &,fmt='(SP,3(ES15.8,X))',advance='no') tp(i)%x
#ifndef RWALK
                 if (dump_velocities) write (unit=cfg_file_unit&
                      &,fmt='(SP,3(ES15.8,X))',advance='no') tp(i)%v
#else
                 tp(i)%v = tp(i)%v+tp(i)%v_r
                 if (dump_velocities) write (unit=cfg_file_unit&
                      &,fmt='(SP,3(ES15.8,X))',advance='no') tp(i)%v
#endif
                 if (dump_forces) write (unit=cfg_file_unit&
                      &,fmt='(SP,3(ES15.8,X))',advance='no') tp(i)%f

                 if (dump_quaternions) write (unit=cfg_file_unit&
                      &,fmt='(SP,4(ES15.8,X))',advance='no') tp(i)%q

                 if (dump_orientations) then
                    o = orientation(tp(i)%q)
                    write (unit=cfg_file_unit,fmt='(SP,3(ES15.8,X))'&
                         &,advance='no') o
                 end if

                 if (dump_rotations) then 
                 A = space_to_body_matrix(tp(i)%q)
                 write (unit=cfg_file_unit&
                      &,fmt='(SP,3(ES15.8,X))',advance='no') &
                 & matmul(transpose(A),tp(i)%w)
                 end if
                 if (dump_torques) write (unit=cfg_file_unit&
                      &,fmt='(SP,3(ES15.8,X))',advance='no') tp(i)%t
#ifdef FORCECOMPONENT		      
                 if (dump_force_components) write (unit=cfg_file_unit&
                      &,fmt='(SP,3(ES15.8,X))',advance='no') tp(i)%f_n
                 if (dump_force_components) write (unit=cfg_file_unit&
                      &,fmt='(SP,3(ES15.8,X))',advance='no') tp(i)%f_t
#endif		      
!#ifdef  MD_MAG
                 if (dump_mag) then 
                  write (unit=cfg_file_unit,fmt='(SP,3(ES15.8,X))',advance='no') tp(i)%fmi 
                  write (unit=cfg_file_unit,fmt='(SP,3(ES15.8,X))',advance='no') tp(i)%tmi
                 end if 
!#endif

#ifdef PARTICLESTRESS
                 if (dump_stresses) then
                    A = space_to_body_matrix(tp(i)%q)
                    write (unit=cfg_file_unit,fmt='(SP,9(ES15.8,X))'&
                         &,advance='no') &
                         &matmul(matmul(transpose(A),tp(i)%tau),A)
                 end if
                 if (dump_stresses_b) write (unit=cfg_file_unit&
                      &,fmt='(SP,9(ES15.8,X))',advance='no') tp(i)%tau
#endif
                 if (dump_potentials) write (unit=cfg_file_unit&
                      &,fmt='(SP,ES15.8,X)',advance='no') tp(i)%e_pot
#ifdef LADD_SURR_RHOF
                 if (dump_surr_rhof) then
                    write (unit=cfg_file_unit&
                         &,fmt='(SP,ES15.8,X)',advance='no') tp(i)%rhof
                    write (unit=cfg_file_unit&
                         &,fmt='(SS,I10.10,X)',advance='no') tp(i)%n_acc
                 end if
#endif
                 if (dump_radii) write (unit=cfg_file_unit&
                      &,fmt='(SP,2(ES15.8,X))',advance='no') &
                      &particle_radii(tp(i))
                 if (dump_ids) write (unit=cfg_file_unit&
                      &,fmt='(SS,I10.10)',advance='no') tp(i)%uid

                 write (unit=cfg_file_unit,fmt='()',advance='yes')
                 close (cfg_file_unit)
              end do
           end if

           if (dump_pnd) then
              call log_msg_md('Dumping particle number density...')
              call lbe_make_filename_output(cfg_file_name,'Particle_Nr_Density'&
                   &,'.asc',t)
              open (unit=cfg_file_unit,file=cfg_file_name,status='REPLACE'&
                   &,action='WRITE',recl=362)
              do j=1,Nz
               	 do k=1,Ny
                    do l=1,Nx
                       Indx = l + Nx*k + Nx*Ny*j
                       write (unit=cfg_file_unit,fmt='(5(I5,X))')&
                            & l,k,j,Indx,B_n(l,k,j)
                    enddo
                 enddo
              end do
              close(cfg_file_unit)
           endif
        end if rank0 
    end subroutine dump_asc

    !> writes prepared data into ascii vtk file on rank 0
    !>
    !> \param[in] prefix prefix used to build the file name
    !>
    !> \param[in] t time step used to build the file name
    !>
    !> \param[in] substep indicates that the dump belongs to a substep
    !> with the given number (optional, defaults to main step dump)
    subroutine dump_vtk(prefix,t,substep)
        character(len=*),intent(in) :: prefix
        integer :: t
        integer,intent(in),optional :: substep
        integer i,n_global
        character(len=1024) cfg_file_name
        integer,parameter :: cfg_file_unit=12
        real(kind=rk) :: o(3),r(2),A(3,3)

        call count_particles(n_global,0)
        rank0: if (myrankc==0) then
           if (.not.dump_positions) call log_msg_md(&
                &'ignoring dump_positions==.false.'&
                &//' - positions are mandatory for vtk output!')

           if (present(substep)) then
              call lbe_make_filename_output_substep(cfg_file_name,prefix&
                   &,'.vtk',t,substep)
           else
              call lbe_make_filename_output(cfg_file_name,prefix,'.vtk',t)
           end if
           open (unit=cfg_file_unit,file=cfg_file_name,status='REPLACE'&
                &,action='WRITE',recl=80)
           write (unit=cfg_file_unit,fmt='(A)') '# vtk DataFile Version 3.0'
           write (unit=cfg_file_unit,fmt='(A)') 'md cfg "'//trim(gr_out_file)&
                &//'"'
           write (unit=cfg_file_unit,fmt='(A)') 'ASCII'
           write (unit=cfg_file_unit,fmt='(A)') 'DATASET POLYDATA'
           write (unit=cfg_file_unit,fmt='(A,I0,A)') 'POINTS ',n_global,' FLOAT'
           do i=1,n_global
!             write (unit=cfg_file_unit,fmt='(SP,3(ES15.8,:,X))') tp(i)%x
             write (unit=cfg_file_unit,fmt='(3(ES15.8,:,X))') tp(i)%x(:)
           end do

#ifndef RWALK
           if (dump_velocities) then
              write (unit=cfg_file_unit,fmt='(A,I0)') 'POINT_DATA ',n_global
              write (unit=cfg_file_unit,fmt='(A)') 'VECTORS velocity FLOAT'
              do i=1,n_global
                 write (unit=cfg_file_unit,fmt='(SP,3(ES15.8,:,X))') tp(i)%v
              end do
           end if
#else
           if (dump_velocities) then
              write (unit=cfg_file_unit,fmt='(A,I0)') 'POINT_DATA ',n_global
              write (unit=cfg_file_unit,fmt='(A)') 'VECTORS velocity FLOAT'
              do i=1,n_global
                 write (unit=cfg_file_unit,fmt='(SP,3(ES15.8,:,X))')&
                      & tp(i)%v+tp(i)%v_r
              end do
           end if
#endif

           if (dump_forces) then
              write (unit=cfg_file_unit,fmt='(A,I0)') 'POINT_DATA ',n_global
              write (unit=cfg_file_unit,fmt='(A)') 'VECTORS force FLOAT'
              do i=1,n_global
                 write (unit=cfg_file_unit,fmt='(SP,3(ES15.8,:,X))') tp(i)%f
              end do
           end if

           if (dump_quaternions) then
              write (unit=cfg_file_unit,fmt='(A,I0)') 'POINT_DATA ',n_global
              write (unit=cfg_file_unit,fmt='(A)') 'SCALARS quaternion FLOAT 4'
              do i=1,n_global
                 write (unit=cfg_file_unit,fmt='(SP,4(ES15.8,:,X))') tp(i)%q
              end do
           end if

           if (dump_orientations) then
              write (unit=cfg_file_unit,fmt='(A,I0)') 'POINT_DATA ',n_global
              write (unit=cfg_file_unit,fmt='(A)') 'VECTORS orientation FLOAT'
              do i=1,n_global
                 o = orientation(tp(i)%q)
                 write (unit=cfg_file_unit,fmt='(SP,3(ES15.8,:,X))') o
              end do
           end if

           if (dump_rotations) then
              write (unit=cfg_file_unit,fmt='(A,I0)') 'POINT_DATA ',n_global
              write (unit=cfg_file_unit,fmt='(A)') 'VECTORS rotation FLOAT'
              do i=1,n_global
                  A = space_to_body_matrix(tp(i)%q)  
                  write (unit=cfg_file_unit,fmt='(SP,3(ES15.8,:,X))') matmul(transpose(A),tp(i)%w)
              end do
           end if

           if (dump_torques) then
              write (unit=cfg_file_unit,fmt='(A,I0)') 'POINT_DATA ',n_global
              write (unit=cfg_file_unit,fmt='(A)') 'VECTORS torque FLOAT'
              do i=1,n_global
                 write (unit=cfg_file_unit,fmt='(SP,3(ES15.8,:,X))') tp(i)%t
              end do
           end if

           if (dump_potentials) then
              write (unit=cfg_file_unit,fmt='(A,I0)') 'POINT_DATA ',n_global
              write (unit=cfg_file_unit,fmt='(A)') 'SCALARS e_pot FLOAT'
              do i=1,n_global
                 write (unit=cfg_file_unit,fmt='(SP,ES15.8)') tp(i)%e_pot
              end do
           end if

#ifdef LADD_SURR_RHOF
           if (dump_surr_rhof) then
              write (unit=cfg_file_unit,fmt='(A,I0)') 'POINT_DATA ',n_global
              write (unit=cfg_file_unit,fmt='(A)') 'SCALARS rhof FLOAT'
              do i=1,n_global
                 write (unit=cfg_file_unit,fmt='(SP,ES15.8)') tp(i)%rhof
              end do

              write (unit=cfg_file_unit,fmt='(A,I0)') 'POINT_DATA ',n_global
              write (unit=cfg_file_unit,fmt='(A)') 'SCALARS n_acc INT'
              do i=1,n_global
                 write (unit=cfg_file_unit,fmt='(SS,I10.10)') tp(i)%n_acc
              end do
           end if

#endif
           if (dump_radii) then
              write (unit=cfg_file_unit,fmt='(A,I0)') 'POINT_DATA ',n_global
              write (unit=cfg_file_unit,fmt='(A)') 'SCALARS R_orth FLOAT'
              do i=1,n_global
                 r = particle_radii(tp(i))
                 write (unit=cfg_file_unit,fmt='(SP,ES15.8)') r(1)
              end do

              write (unit=cfg_file_unit,fmt='(A,I0)') 'POINT_DATA ',n_global
              write (unit=cfg_file_unit,fmt='(A)') 'SCALARS R_para FLOAT'
              do i=1,n_global
                 r = particle_radii(tp(i))
                 write (unit=cfg_file_unit,fmt='(SP,ES15.8)') r(2)
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
    !>
    !> \param[in] prefix prefix used to build the file name
    !>
    !> \param[in] t time step used to build the file name
    !>
    !> \param[in] substep indicates that the dump belongs to a substep
    !> with the given number (optional, defaults to main step dump)
    subroutine dump_xdr(prefix,t,substep)
        character(len=*),intent(in) :: prefix
        integer :: t
        integer,intent(in),optional :: substep
        integer i,ierror,cfg_file_handle,n_global
        character(len=1024) cfg_file_name
        integer,parameter :: fk=4
        real(kind=rk) :: o(3)
        real(kind=fk) :: tmp(0:3)
        real(kind=rk) :: A(3,3),ws(3)
#ifdef RWALK
        real(kind=rk) :: tmp2(0:3)
#endif
#ifdef PARTICLESTRESS
        real(kind=fk) :: tmpt(3,3)
        real(kind=rk) :: T(3,3)
#endif

        ! These would be function pointers in C. In fact,  xdrf(double|float)
        ! ares external functions, but fortran does not seem to care about the
        ! type anyway...
        integer,external :: xdrffloat, xdrfdouble

        call dump_description_on_first_call

        call count_particles(n_global,0)
        rank0: if (myrankc==0) then
#ifdef USEXDRF
           if (present(substep)) then
              call lbe_make_filename_output_substep(cfg_file_name,prefix&
                   &,'.xdr',t,substep)
           else
              call lbe_make_filename_output(cfg_file_name,prefix,'.xdr',t)
           end if
           call xdrfopen(cfg_file_handle,cfg_file_name,'w',ierror)
           if (ierror==0) then
              call log_msg_md('dump_xdr(): xdrfopen() failed')
              call check_xdrfsloppy()
           end if

           double: if (dump_double) then
              do i=1,n_global
                 if (dump_positions) call xdrfvector&
                      &(cfg_file_handle,tp(i)%x(1:3),3,xdrfdouble,ierror)
#ifndef RWALK
                 if (dump_velocities) call xdrfvector&
                      &(cfg_file_handle,tp(i)%v(1:3),3,xdrfdouble,ierror)
#else
                 tmp2(1:3)= tp(i)%v(1:3) + tp(i)%v_r(1:3)
                 if (dump_velocities) call xdrfvector&
                      &(cfg_file_handle,tmp2(1:3),3,xdrfdouble,ierror)
#endif
                 if (dump_forces) call xdrfvector&
                      &(cfg_file_handle,tp(i)%f(1:3),3,xdrfdouble,ierror)
                 if (dump_quaternions) call xdrfvector&
                      &(cfg_file_handle,tp(i)%q(0:3),4,xdrfdouble,ierror)
                 if (dump_orientations) then
                    o = orientation(tp(i)%q)
                    call xdrfvector&
                         &(cfg_file_handle,o(1:3),3,xdrfdouble,ierror)
                 end if
                 if (dump_rotations) then 
                        A = space_to_body_matrix(tp(i)%q)
                        ws= matmul(transpose(A),tp(i)%w)  
                      call xdrfvector&
                      &(cfg_file_handle,ws,3,xdrfdouble,ierror)
                    end if 
                 if (dump_torques) call xdrfvector&
                      &(cfg_file_handle,tp(i)%t(1:3),3,xdrfdouble,ierror)
!#ifdef MD_MAG
                if (dump_mag) then 
                 call xdrfvector&
                 &(cfg_file_handle,tp(i)%fmi(1:3),3,xdrfdouble,ierror)
                call xdrfvector&
                 &(cfg_file_handle,tp(i)%tmi(1:3),3,xdrfdouble,ierror)
                end if 
!#endif
#ifdef PARTICLESTRESS
                 if (dump_stresses) then
                    A = space_to_body_matrix(tp(i)%q)
                    T = matmul(matmul(transpose(A),tp(i)%tau),A)
                    call xdrfvector(cfg_file_handle,T,9,xdrfdouble,ierror)
                 end if
                 if (dump_stresses_b) call xdrfvector&
                      &(cfg_file_handle,tp(i)%tau,9,xdrfdouble,ierror)
#endif
                 if (dump_potentials) call xdrfvector(&
                      &cfg_file_handle,tp(i)%e_pot,1,xdrfdouble,ierror)
#ifdef LADD_SURR_RHOF
                 if (dump_surr_rhof) then
                    call xdrfvector(cfg_file_handle,tp(i)%rhof,1,xdrfdouble&
                         &,ierror)
                    call xdrfint(cfg_file_handle,tp(i)%n_acc,ierror)
                 end if
#endif
                 if (dump_radii) call xdrfvector&
                      &(cfg_file_handle,particle_radii(tp(i)),2,xdrfdouble&
                      &,ierror)
                 if (dump_ids) call xdrfint(cfg_file_handle,tp(i)%uid,ierror)
              end do
           else double
              do i=1,n_global
                 if (dump_positions) then
                    tmp(1:3) = real(tp(i)%x(1:3),kind=fk)
                    call xdrfvector(cfg_file_handle,tmp(1:3),3,xdrffloat,ierror)
                 end if
                 if (dump_velocities) then
#ifndef RWALK
                    tmp(1:3) = real(tp(i)%v(1:3),kind=fk)
#else
                    tmp(1:3) = real(tp(i)%v(1:3),kind=fk)&
                         & + real(tp(i)%v_r(1:3),kind=fk)
#endif
                    call xdrfvector(cfg_file_handle,tmp(1:3),3,xdrffloat,ierror)
                 end if
                 if (dump_forces) then
                    tmp(1:3) = real(tp(i)%f(1:3),kind=fk)
                    call xdrfvector(cfg_file_handle,tmp(1:3),3,xdrffloat,ierror)
                 end if
                 if (dump_quaternions) then
                    tmp(0:3) = real(tp(i)%q(0:3),kind=fk)
                    call xdrfvector(cfg_file_handle,tmp(0:3),4,xdrffloat,ierror)
                 end if
                 if (dump_orientations) then
                    o = orientation(tp(i)%q)
                    tmp(1:3) = real(o(1:3),kind=fk)
                    call xdrfvector(cfg_file_handle,tmp(1:3),3,xdrffloat,ierror)
                 end if
                 if (dump_rotations) then
                    A = space_to_body_matrix(tp(i)%q)
                    tmp(1:3) = real(matmul(transpose(A),tp(i)%w),kind=fk)
                    call xdrfvector(cfg_file_handle,tmp(1:3),3,xdrffloat,ierror)
                 end if
                 if (dump_torques) then
                    tmp(1:3) = real(tp(i)%t(1:3),kind=fk)
                    call xdrfvector(cfg_file_handle,tmp(1:3),3,xdrffloat,ierror)
                 end if
!#ifdef MD_MAG
                 if (dump_mag) then 
                    tmp(1:3) = real(tp(i)%fmi(1:3),kind=fk)                                                                                                                                  
                    call xdrfvector(cfg_file_handle,tmp(1:3),3,xdrffloat,ierror)                                                                                                             
                    tmp(1:3) = real(tp(i)%tmi(1:3),kind=fk)                                                                                                                                 
                    call xdrfvector(cfg_file_handle,tmp(1:3),3,xdrffloat,ierror) 
                 end if
!#endif

#ifdef PARTICLESTRESS
                 if (dump_stresses) then
                    A = space_to_body_matrix(tp(i)%q)
                    tmpt = real(matmul(matmul(transpose(A),tp(i)%tau),A)&
                         &,kind=rk)
                    call xdrfvector(cfg_file_handle,tmpt,9,xdrffloat,ierror)
                 end if
                 if (dump_stresses_b) then
                    tmpt = real(tp(i)%tau,kind=fk)
                    call xdrfvector(cfg_file_handle,tmpt,9,xdrffloat,ierror)
                 end if
#endif
                 if (dump_potentials) then
                    tmp(0) = real(tp(i)%e_pot,kind=fk)
                    call xdrfvector(cfg_file_handle,tmp(0),1,xdrffloat,ierror)
                 end if
#ifdef LADD_SURR_RHOF
                 if (dump_surr_rhof) then
                    tmp(0) = real(tp(i)%rhof,kind=fk)
                    call xdrfvector(cfg_file_handle,tmp(0),1,xdrffloat,ierror)
                    call xdrfint(cfg_file_handle,tp(i)%n_acc,ierror)
                 end if
#endif
                 if (dump_radii) then
                    tmp(1:2) = real(particle_radii(tp(i)),kind=fk)
                    call xdrfvector(cfg_file_handle,tmp(1:2),2,xdrffloat,ierror)
                 end if
                 if (dump_ids) then
                    call xdrfint(cfg_file_handle,tp(i)%uid,ierror)
                 end if
              end do
           end if double

           call xdrfclose(cfg_file_handle,ierror)
#else
           call log_msg_md('failed to dump '//cfg_file_name&
                &//' - recompile with USEXDRF set!')
#endif
        end if rank0
    end subroutine dump_xdr

    !> dumps z-averaged v_z(r)-profile to file "md-rvz_<...>.asc".
    !>
    !> Each line contains six columns: r-coordinate, averaged fluid
    !> velocity, averaged particle velocity, average particle (center)
    !> density, average fluid node density, average moving rock node
    !> density, and average resting rock node density. Depending on
    !> n_rvz_profile_sample and n_rvz_profile_dump the dumped result
    !> is the time-average of several samples.
    subroutine dump_rvz_profile(N)
        type(lbe_site),intent(in) ::&
             & N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        integer,parameter :: rvz_file_unit=12
        character(len=1024) rvz_file_name
        real(kind=rk),allocatable,save :: vf(:),vp(:),vfsum(:),vpsum(:)
        integer,allocatable,save :: nf(:),np(:),nnr(:),nmr(:),nfr(:),nfsum(:)&
             &,npsum(:),nnrsum(:),nmrsum(:),nfrsum(:)
        integer,save :: maxr,n_samples
        real(kind=rk),save :: midpos(2)
        integer i,ii,x,y,z,r,ierror,stat
        real(kind=rk) :: dst(2),total_nodes

        initialization: if (.not.allocated(vf)) then
           ! first call of this subroutine: perform initialization
           midpos = minpos(1:2)+real((/tnx,tny/),kind=rk)/2.0_rk
           maxr = nint(sqrt(dot_product(1.0_rk-midpos(1:2),1.0_rk-midpos(1:2))))

           allocate(vf(0:maxr),vp(0:maxr),nf(0:maxr),np(0:maxr),nnr(0:maxr)&
                &,nmr(0:maxr),nfr(0:maxr),stat=stat)
           call check_allocate(stat,'dump_rvz_profile:vf,vp,nf,np,nnr,nmr,nfr')

           vf(:) = 0.0_rk
           vp(:) = 0.0_rk
           nf(:) = 0               ! number of fluid sites
           np(:) = 0               ! number of particles
           nnr(:) = 0              ! number of non-solid-rock sites
           nmr(:) = 0              ! number of moving rock sites
           nfr(:) = 0              ! number of fixed rock sites

           n_samples = 0

           if (myrankc==0) then
              allocate(vfsum(0:maxr),vpsum(0:maxr),nfsum(0:maxr),npsum(0:maxr)&
                   &,nnrsum(0:maxr),nmrsum(0:maxr),nfrsum(0:maxr),stat=stat)
              call check_allocate(stat&
                   &,'dump_rvz_profile:vfsum,vpsum,nfsum,npsum,nnrsum,nmrsum'&
                   &//',nfrsum')

              if (n_rvz_profile_dump==0) &
                   &call error_md('n_rvz_profile_dump  must not be zero '&
                   &//'if  n_rvz_profile_sample/=0 !')
           end if
        end if initialization

        ! fluid velocity
        do x=1,nx
           do y=1,ny
              dst = real((/x+start(1)-1,y+start(2)-1/),kind=rk)-midpos
              r = nint(sqrt(dot_product(dst,dst)))
              do z=1,nz
                 if (N(x,y,z)%rock_state==0.0_rk) then
                    nf(r) = nf(r) + 1
                    vf(r) = vf(r) + &
#ifdef SINGLEFLUID
                         &sum(g*N(x,y,z)%n_r*cz)/&
                         &sum(g*N(x,y,z)%n_r)
#else
                         &sum(g*cz*&
                              &(N(x,y,z)%n_r*amass_r&
                              &+N(x,y,z)%n_b*amass_b&
#ifndef NOSURFACTANT
                              &+N(x,y,z)%n_s*amass_s&
#endif
                              &))/&
                         &sum(g*&
                              &(N(x,y,z)%n_r*amass_r&
                              &+N(x,y,z)%n_b*amass_b&
#ifndef NOSURFACTANT
                              &+N(x,y,z)%n_s*amass_s&
#endif
                              &))
#endif
                    nnr(r) = nnr(r) + 1
                 else if (interaction=='ladd'.and.N(x,y,z)%rock_state>0.0_rk) &
                      &then
                    ! for Ladd-coupling solid rock has rock_state<0,
                    ! rock_state>0 means that it is part of a particle
                    ! and should be counted to the fluid for
                    ! calculating the particle-center-density
                    nnr(r) = nnr(r) + 1
                    nmr(r) = nmr(r) + 1
                 else
                    ! neither fluid nor ladd particle - it must be solid rock
                    nfr(r) = nfr(r) + 1
                 end if
              end do
           end do
        end do

        ! particle velocity
        i = atompnt
        do ii=1,nlocal
           dst = P(i)%x(1:2)-midpos
           r = nint(sqrt(sum(dst**2)))
           np(r) = np(r) + 1
           vp(r) = vp(r) + P(i)%v(3)
           i=list(i)
        end do

        n_samples = n_samples + 1

        dump: if (mod(nt,n_rvz_profile_dump)==0) then
           if (myrankc==0) then
              write (unit=6,fmt='(A,SS,I9,A)') &
                   &'v_z(r) was averaged over',n_samples,' samples'
              call log_msg_md('dumping v_z(r)...')
           end if

           call mpi_reduce(vf,vfsum,maxr+1,MPI_REAL8,MPI_SUM,0,comm_cart,ierror)
           call mpi_reduce(vp,vpsum,maxr+1,MPI_REAL8,MPI_SUM,0,comm_cart,ierror)
           call mpi_reduce(nf,nfsum,maxr+1,MPI_INTEGER,MPI_SUM,0,comm_cart&
                &,ierror)
           call mpi_reduce(np,npsum,maxr+1,MPI_INTEGER,MPI_SUM,0,comm_cart&
                &,ierror)
           call mpi_reduce(nnr,nnrsum,maxr+1,MPI_INTEGER,MPI_SUM,0,comm_cart&
                &,ierror)
           call mpi_reduce(nmr,nmrsum,maxr+1,MPI_INTEGER,MPI_SUM,0,comm_cart&
                &,ierror)
           call mpi_reduce(nfr,nfrsum,maxr+1,MPI_INTEGER,MPI_SUM,0,comm_cart&
                &,ierror)

           rank0: if (myrankc==0) then
              vfsum = vfsum/(real(nfsum,kind=rk))
              vpsum = vpsum/(real(npsum,kind=rk))

              call lbe_make_filename_output(rvz_file_name,'md-rvz','.asc',nt)
              open (unit=rvz_file_unit,file=rvz_file_name,status='REPLACE'&
                   &,action='WRITE',recl=160)

              do i=0,maxr
                 total_nodes = real(nfsum(i)+nmrsum(i)+nfrsum(i),kind=rk)
                 write (unit=rvz_file_unit,fmt='(SS,I5.5,X,SP,6(ES15.8,:,X))') &
                      &i,vfsum(i),vpsum(i)&
                      &,real(npsum(i),kind=rk)/real(nnrsum(i),kind=rk)&
                      &,real(nfsum(i),kind=rk)/total_nodes&
                      &,real(nmrsum(i),kind=rk)/total_nodes&
                      &,real(nfrsum(i),kind=rk)/total_nodes
              end do

              close (rvz_file_unit)
           end if rank0

           vf = 0.0_rk
           vp = 0.0_rk
           nf = 0               ! number of fluid sites
           np = 0               ! number of particles
           nnr = 0              ! number of non-solid-rock sites
           nmr = 0              ! number of moving rock sites
           nfr = 0              ! number of solid rock sites

           n_samples = 0
        end if dump
    end subroutine dump_rvz_profile

    !> write a checkpoint of the current state
    subroutine md_dump_checkpoint
	character(len=1024) chk_file_name
        real(kind=rk) :: tmp(3)
        integer chk_file_handle,i,ierror,n_global

        ! This would be a function pointer in C. In fact,  xdrfdouble  is
        ! an external function, but fortran does not seem to care about the
        ! type anyway...
        integer,external :: xdrfdouble

        call log_msg_md("  Dumping MD checkpoint ...")

        if ( checkpoint_format .eq. 'xdr') then
#ifdef USEXDRF
           if (myrankc==0) then
              call lbe_make_filename_cp(chk_file_name,'md-checkpoint','.xdr',nt)
              call xdrfopen(chk_file_handle,chk_file_name,'w',ierror)
              if (ierror==0) then
                 call log_msg_md('md_dump_checkpoint(): xdrfopen() failed')
                 call check_xdrfsloppy()
              end if
           end if

           ! ...at last, collect and dump the configuration of all particles:
           call gather_particles_provide_rbuf(tp,'md_dump_checkpoint(): tp')
           call gather_particles(tp,particle_cp_mpitype)

           call count_particles(n_global,0)
           if (myrankc==0) then
              call xdrfint(chk_file_handle,n_global,ierror)

              do i=1,n_global
                 call xdrfvector(chk_file_handle,tp(i)%x(1:3),3,xdrfdouble&
                      &,ierror)
                 call xdrfvector(chk_file_handle,tp(i)%v(1:3),3,xdrfdouble&
                      &,ierror)
                 if (use_rotation) then
                    call xdrfvector(chk_file_handle,tp(i)%q(0:3),4,xdrfdouble&
                         &,ierror)
                    call xdrfvector(chk_file_handle,tp(i)%w(1:3),3,xdrfdouble&
                         &,ierror)
                 end if
                 call xdrfint(chk_file_handle,tp(i)%uid,ierror)
                 call xdrfvector(chk_file_handle,tp(i)%vnew(1:3),3,xdrfdouble&
                      &,ierror)
                 if (use_rotation) then
                    call xdrfvector(chk_file_handle,tp(i)%qnew(0:3),4&
                         &,xdrfdouble,ierror)
                    call xdrfvector(chk_file_handle,tp(i)%wnew(1:3),3&
                         &,xdrfdouble,ierror)
                 end if
                 call xdrfint(chk_file_handle,tp(i)%master,ierror)

                 if (polydispersity) then
                    call xdrfvector(chk_file_handle,tp(i)%R_orth,1,xdrfdouble&
                         &,ierror)
                    call xdrfvector(chk_file_handle,tp(i)%R_para,1,xdrfdouble&
                         &,ierror)
                 end if
                 if (magdispersity) then
                      call xdrfvector(chk_file_handle,tp(i)%mag,1,xdrfdouble&
                           &,ierror)
                   end if
                

                 ! write averaged substep velocities in a way that
                 ! allows to change steps_per_lbe_step arbitrarily
                 ! when restoring the checkpoint
                 if (steps_per_lbe_step>1) then
                    ! store preliminary averages to make accumulators
                    ! independent from steps_per_lbe_step
                    tmp = tp(i)%v_fluid_acc/real(steps_per_lbe_step-1,kind=rk)
                    call xdrfvector(chk_file_handle,tmp,3,xdrfdouble,ierror)
                    if (use_rotation) then
                       tmp = tp(i)%ws_fluid_acc&
                            &/real(steps_per_lbe_step-1,kind=rk)
                       call xdrfvector(chk_file_handle,tmp,3,xdrfdouble,ierror)
                    end if
                 else
                    ! for steps_per_lbe_step==1, accumulators were not
                    ! set, store current velocities instead
                    call xdrfvector(chk_file_handle,tp(i)%v,3,xdrfdouble,ierror)
                    if (use_rotation) then
                       call xdrfvector(chk_file_handle,tp(i)%ws,3,xdrfdouble&
                            &,ierror)
                    end if
                 end if

                 call xdrfvector(chk_file_handle,tp(i)%v_fluid(1:3),3&
                      &,xdrfdouble,ierror)
                 if (use_rotation) then
                    call xdrfvector(chk_file_handle,tp(i)%ws_fluid(1:3),3&
                         &,xdrfdouble,ierror)
                 end if

                 call xdrfvector(chk_file_handle,tp(i)%f_fluid(1:3),3&
                      &,xdrfdouble,ierror)
                 if (use_rotation) then
                    call xdrfvector(chk_file_handle,tp(i)%t_fluid(1:3),3&
                         &,xdrfdouble,ierror)
                 end if

#ifdef LADD_SURR_RHOF
                 call xdrfvector(chk_file_handle,tp(i)%rhof,1,xdrfdouble,ierror)
#endif
              end do

              ! global variables used only by ladd code
              ladd: if (interaction=='ladd') then

                call xdrfvector(chk_file_handle,(/ pfr, pfb, pfg /),3,xdrfdouble&
                     &,ierror)
                call xdrfvector(chk_file_handle,global_mass_change,n_spec,xdrfdouble,ierror)
                call xdrfvector(chk_file_handle,global_mass_target,n_spec,xdrfdouble,ierror)
              end if ladd

              periodic_inflow: if (boundary=='periodic_inflow') then
                 call xdrfint(chk_file_handle,first_inserted_puid,ierror)
              end if periodic_inflow

              call xdrfclose(chk_file_handle,ierror)
           end if
#else
           ! ifndef USEXDRF
           call log_msg_md('failed to dump '//chk_file_name&
                &//' - recompile with USEXDRF set!')
#endif

        else if ( checkpoint_format .eq. 'hdf' ) then
#ifdef USEHDF
           call log_msg_md(&
                &'WARNING: MD HDF5 checkpointing not yet implemented!')
#else
           ! ifndef USEHDF
           call log_msg_md('WARNING: failed to dump MD checkpoint '&
                &//'- recompile with USEHDF set!')
#endif
        endif

        call log_msg_md("  Finished dumping MD checkpoint.")

    end subroutine md_dump_checkpoint

    !> This routine deletes old MD checkpoint files.
    subroutine md_delete_checkpoint(last)
      implicit none

      integer, intent(in) :: last

      character(len=1024) :: filename

      call log_msg_md("  Deleting MD checkpoint files ...")

      if (myrankc == 0) then
        if ( checkpoint_format .eq. 'xdr') then
          call lbe_make_filename_cp(filename, 'md-checkpoint', '.xdr', last)
          call lbe_delete_file(filename)
        else if (checkpoint_format .eq. 'hdf' ) then
          call lbe_make_filename_cp(filename, 'cp_md_', '.h5', last)
          call lbe_delete_file(filename)
        end if
      end if

    end subroutine md_delete_checkpoint

    !> output to screen and file, add long-range correction to energy and
    !> pressure
    subroutine summary
        character(len=1024) sum_file_name,stat_file_name
        integer ihisto(10),ihistotmp(10)
        integer i,ii,n_global,nlost,nnlost,tmpmax,tstat,ierror,mbinmax
        real(kind=rk) rtmp,ave,xmax,xmin
901     format(' ',a,f13.4,a,f13.4,a,f13.4,a)
902     format(' ',A,10(I6,:,X))

        call log_msg_md_hdr("Summary")

        call count_particles_all(n_global)

        ! check for lost atoms
        nnlost = 0
        i = atompnt
        do ii = 1,nlocal
           if (any(P(i)%x(:)<minpos(:)).or.any(P(i)%x(:)>=maxpos(:)))&
                &nnlost = nnlost + 1
           i = list(i)
        enddo
        call mpi_allreduce(nnlost,nlost,1,MPI_INTEGER,MPI_SUM,comm_cart,&
             &ierror)

        if (nlost.gt.0) then
           write (6,*) 'Particle counts: n_global=',n_global,', nlost=',nlost
           call error_md('Lost particles')
        endif

        ! output

        if (myrankc==0) then
           ! Summary still belongs to timestep  nt-1 , because this subroutine
           ! is called outside the time loop after the last iteration.
           call lbe_make_filename_output(sum_file_name,'md-summary','.txt',nt-1)
           open (unit=30,file=sum_file_name,recl=260)
           stats: if (mstat>0) then
              call lbe_make_filename_output(stat_file_name,'md-stats','.asc',nt-1)
              open (unit=31,file=stat_file_name,recl=260)
              write (30,*) 'Timestep, T*, T_trans, T_rot, U*, U_rock*, P*, '&
                   &//'Conservation, Momentum:'
              write (31,fmt='(A)') '# t        temp            '&
                   &//'T_trans         T_rot           U*              '&
                   &//'U_rock*         P*              conservation    '&
                   &//'momentum(x)     momentum(y)     momentum(z)'
              do i = 1,mstat
                 tstat = (i-1)*n_stat
                 if ( is_restoring() ) tstat = tstat + n_restore
                 write (30,fmt='(I10.10,X,SP,10(ES15.8,:,X))') tstat,tmparr(i)&
                      &,e_trans(i),e_rot(i),engarr(i),rpotarr(i),prsarr(i)&
                      &,conarr(i),momentumarr(:,i)
                 write (31,fmt='(I10.10,X,SP,10(ES15.8,:,X))') tstat,tmparr(i)&
                      &,e_trans(i),e_rot(i),engarr(i),rpotarr(i),prsarr(i)&
                      &,conarr(i),momentumarr(:,i)
              enddo
              if (mstat>1) then
                 write (30,*) 'Averages: (without timestep 0)'
                 write (30,fmt='(11X,SP,10(ES15.8,:,X))') tmpave,e_trans_ave&
                      &,e_rot_ave,engave,rpotave,prsave,conave,momentumave(:)
              endif

              print *
              write (30,*)
           end if stats
        endif

        if (myrankc==0) then
           print *,'Run on',nprocs,' procs for',n_global,' particles'
           write (30,*) 'Run on',nprocs,' procs for' ,n_global,' particles'
           print *
           write (30,*)
        end if

        rtmp = nnlist(nlocal+1) - 1
        call stats_rk(rtmp,ave,xmax,xmin,ihisto,ihistotmp,10)
        if (myrankc==0) then
           print 901,'Neighs:    ',ave,' ave',xmax,' max',xmin,' min'
           print 902,' Histogram:',(ihisto(i),i=1,10)
           write (30,901) 'Neighs:    ',ave,' ave',xmax,' max',xmin,' min'
           write (30,902) ' Histogram:',(ihisto(i),i=1,10)
        endif

        rtmp = nlocal
        call stats_rk(rtmp,ave,xmax,xmin,ihisto,ihistotmp,10)
        if (myrankc==0) then
           print 901,'Nlocal:    ',ave,' ave',xmax,' max',xmin,' min'
           print 902,' Histogram:',(ihisto(i),i=1,10)
           write (30,901) 'Nlocal:    ',ave,' ave',xmax,' max',xmin,' min'
           write (30,902) ' Histogram:',(ihisto(i),i=1,10)
        endif

        rtmp = nother
        call stats_rk(rtmp,ave,xmax,xmin,ihisto,ihistotmp,10)
        if (myrankc==0) then
           print 901,'Nother:    ',ave,' ave',xmax,' max',xmin,' min'
           print 902,' Histogram:',(ihisto(i),i=1,10)
           write (30,901) 'Nother:    ',ave,' ave',xmax,' max',xmin,' min'
           write (30,902) ' Histogram:',(ihisto(i),i=1,10)
        endif

        rtmp = nslistmax
        call stats_rk(rtmp,ave,xmax,xmin,ihisto,ihistotmp,10)
        if (myrankc==0) then
           print 901,'Nswaps:    ',ave,' ave',xmax,' max',xmin,' min'
           print 902,' Histogram:',(ihisto(i),i=1,10)
           write (30,901) 'Nswaps:    ',ave,' ave',xmax,' max',xmin,' min'
           write (30,902) ' Histogram:',(ihisto(i),i=1,10)
        endif

        if (myrankc==0) then
           print *
           write (30,*)
        endif

        tmpmax = nlocalmax
        call mpi_allreduce(tmpmax,nlocalmax,1,MPI_INTEGER,MPI_MAX,&
             &comm_cart,ierror)
        tmpmax = nothermax
        call mpi_allreduce(tmpmax,nothermax,1,MPI_INTEGER,MPI_MAX,&
             &comm_cart,ierror)
        tmpmax = neighmax
        call mpi_allreduce(tmpmax,neighmax,1,MPI_INTEGER,MPI_MAX,&
             &comm_cart,ierror)
        tmpmax = nslistmax
        call mpi_allreduce(tmpmax,nslistmax,1,MPI_INTEGER,MPI_MAX,&
             &comm_cart,ierror)
        tmpmax = nexcmax
        call mpi_allreduce(tmpmax,nexcmax,1,MPI_INTEGER,MPI_MAX,&
             &comm_cart,ierror)
        tmpmax = nswpmax
        call mpi_allreduce(tmpmax,nswpmax,1,MPI_INTEGER,MPI_MAX,&
             &comm_cart,ierror)
        tmpmax = mbinx*mbiny*mbinz
        call mpi_allreduce(tmpmax,mbinmax,1,MPI_INTEGER,MPI_MAX,&
             &comm_cart,ierror)

        if (myrankc==0) then
           write (30,*) 'Max # of local atoms =',nlocalmax,' out of',npmax
           write (30,*) 'Max # of other atoms =',nothermax,' out of',nomax
           write (30,*) 'Max # of neighbors =',neighmax,' out of',nnmax*npmax
           write (30,*) 'Max size of swap list =',nslistmax,' out of',nemax
           write (30,*) 'Max used in exchange buffer =',nexcmax,'&
                & out of',nfmax
           write (30,*) 'Max used in swap buffer =',nswpmax
           if (mbinmax.gt.0) write (30,*) 'Max # of bins =',mbinmax
           write (30,*)
           write (30,*) 'rs =',rs,' Cut/Box =',rs*cdims(:)/tsize(:)
           write (30,*)
           write (30,*) '# of reneighborings =',mneigh
        endif

        call summary_fluid((/6,30/))

        if (myrankc==0) then
           close (30)
           if (mstat>0) close (31)
        end if
    end subroutine summary

    subroutine dump_description_on_first_call
        logical,save :: first_call=.true.
        character(len=1024) cfgdesc_file_name
        integer,parameter :: cfgdesc_file_unit=12

        if (.not.first_call) return
        first_call = .false.

        rank0: if (myrankc==0) then
           call lbe_make_filename_output(cfgdesc_file_name,'md-cfg-desc','.txt'&
                &,nt)
           open (unit=cfgdesc_file_unit,file=cfgdesc_file_name&
                &,status='REPLACE',action='WRITE',recl=650)

           if (dump_time) write (unit=cfgdesc_file_unit&
                &,fmt='(A10,1X)',advance='no') 'time'
           if (dump_positions) write (unit=cfgdesc_file_unit&
                &,fmt='(3(A15,X))',advance='no') 'x','y','z'
           if (dump_velocities) write (unit=cfgdesc_file_unit&
                &,fmt='(3(A15,X))',advance='no') 'vx','vy','vz'
           if (dump_forces) write (unit=cfgdesc_file_unit&
                &,fmt='(3(A15,X))',advance='no') 'fx','fy','fz'
           if (dump_quaternions) write (unit=cfgdesc_file_unit&
                &,fmt='(4(A15,X))',advance='no') 'q0','q1','q2','q3'
           if (dump_orientations) write (unit=cfgdesc_file_unit&
                &,fmt='(3(A15,X))',advance='no') 'ox','oy','oz'
           if (dump_rotations) write (unit=cfgdesc_file_unit&
                &,fmt='(3(A15,X))',advance='no') 'wxb','wyb','wzb'
           if (dump_torques) write (unit=cfgdesc_file_unit&
                &,fmt='(3(A15,X))',advance='no') 'tx','ty','tz'
#ifdef FORCECOMPONENT		
           if (dump_force_components) write (unit=cfgdesc_file_unit&
                &,fmt='(3(A15,X))',advance='no') 'fnx','fny','fnz'
           if (dump_force_components) write (unit=cfgdesc_file_unit&
                &,fmt='(3(A15,X))',advance='no') 'ftx','fty','ftz'
#endif		
      !#ifdef MD_MAG

           if (dump_mag) then 
            write (unit=cfgdesc_file_unit&
            &,fmt='(6(A15,X))',advance='no') 'fmix','fmiy','fmiz'
            write (unit=cfgdesc_file_unit&
            &,fmt='(6(A15,X))',advance='no') 'tmix','tmiy','tmiz'
          end if
!#endif
#ifdef PARTICLESTRESS
           if (dump_stresses) write (unit=cfgdesc_file_unit&
                &,fmt='(9(A15,X))',advance='no')&
                & 'sxx','sxy','sxz'&
                &,'syx','syy','syz'&
                &,'szx','szy','szz'
           if (dump_stresses_b) write (unit=cfgdesc_file_unit&
                &,fmt='(9(A15,X))',advance='no')&
                & 'sxxb','sxyb','sxzb'&
                &,'syxb','syyb','syzb'&
                &,'szxb','szyb','szzb'
#endif
           if (dump_potentials) write (unit=cfgdesc_file_unit&
                &,fmt='(A15,X)',advance='no') 'e_pot'
#ifdef LADD_SURR_RHOF
           if (dump_surr_rhof) then
              write (unit=cfgdesc_file_unit,fmt='(A15,X)',advance='no') 'rhof'
              write (unit=cfgdesc_file_unit,fmt='(A10,X)',advance='no') 'n_acc'
           end if
#endif



           if (dump_radii) write (unit=cfgdesc_file_unit&
                &,fmt='(2(A15,X))',advance='no') 'R_orth','R_para'
           if (dump_ids) write (unit=cfgdesc_file_unit&
                &,fmt='(A10)',advance='no') 'uid'

           write (unit=cfgdesc_file_unit,fmt='()',advance='yes')
           close (cfgdesc_file_unit)
        end if rank0
    end subroutine dump_description_on_first_call

#endif
end module lbe_md_output_module
