#include "lbe.h"
#include "lbe_version.h"

!>Contains simulation parameters, which may or may not vary across PEs.
module lbe_parms_module

  use lbe_globals_module, only: rk
  use lbe_types_module, only : lbe_site

  implicit none

  public

  ! All use statements put at the top of the module should be made private here
  ! to preserve the cleanliness of the module structures.
  private rk, lbe_site

  character(len=256),   parameter :: lbeversion  = 'LB3D-version6 ('//trim(GIT_DESC)//trim(GIT_LOCALCHANGES)//', '//trim(GIT_BRANCH)//')'
  character(len=1024), parameter :: lbeflags    = 'LB3D flags: '//trim(LBE_FLAGS)
  character(len=1024), parameter :: lbeplatform = 'LB3D platform: '//trim(LBE_PLATFORM)

  integer, parameter :: max_flag_len = 256
  character(len=max_flag_len), allocatable, dimension(:) :: flag_list

  ! ======================================================================
  !                                FIXED INPUT
  ! ======================================================================

  integer, save :: nx = 16
  integer, save :: ny = 16
  integer, save :: nz = 16
  integer, save :: seed = 1111
  character(len=256) :: obs_file = 'empty.dat'
  logical, save :: gw_rock_colour = .false.
  logical, save :: rock_colour_double = .false. 
  logical, save :: g_value_double = .false.
  logical, save :: gw_double_wet = .false.
  logical, save :: gw_wet_global = .true.
  character(len=256) :: obs_file_r = 'empty.dat'
  character(len=256) :: obs_file_b = 'empty.dat'
  character(len=256) :: obs_folder = '../rocks'
  character(len=3)  :: obs_rotation = 'xyz'
  integer, save :: boundary_cond  = 0
  integer, save :: boundary_width = 2
  character(len=64) :: boundary = 'periodic'
  character(len=64) :: force = 'none'
  !> enable \c lbe_galilean_stabilizer_module
  logical, save :: galilean_stabilizer = .false.
  !> enable \c lbe_mass_scaler_module
  logical, save :: mass_scaler = .false.
  integer, save :: cdx = 0
  integer, save :: cdy = 0
  integer, save :: cdz = 0
  logical, save :: dbg_report_topology = .false.
  integer, save :: dbg_n_start_mpi_debug = 0

  namelist /FIXED_INPUT/ nx,ny,nz,seed,obs_file,obs_file_r, obs_file_b,obs_folder,obs_rotation&
       &,boundary_cond,boundary_width,boundary,force,galilean_stabilizer,cdx&
       &,cdy,cdz,dbg_report_topology,dbg_n_start_mpi_debug,mass_scaler&
       &,gw_rock_colour, rock_colour_double, g_value_double, gw_double_wet, gw_wet_global

  ! ======================================================================
  !                                VARIABLE INPUT
  ! ======================================================================

  logical, save :: sci_int = .false.
  logical, save :: sci_sur = .false.
  logical, save :: sci_od  = .false.
  logical, save :: sci_wd  = .false.
  logical, save :: sci_dir = .false.
  logical, save :: sci_vel = .false.
  logical, save :: sci_flo = .false.
  logical, save :: sci_arrows = .false.
  logical, save :: sci_velocities = .false.
  logical, save :: sci_velocities_od = .false.
  logical, save :: sci_velocities_wd = .false.
  logical, save :: sci_flux_od = .false.
  logical, save :: sci_flux_wd = .false.
  logical, save :: sci_rock = .false.
  logical, save :: sci_rock_colour = .false.
  logical, save :: sci_rock_rho_r = .false.
  logical, save :: sci_rock_rho_b = .false.
  logical, save :: sci_pressure = .false.
  logical, save :: sci_fluxz = .false.
  logical, save :: sci_massfluxz = .false.
  logical, save :: sci_profile = .false.
  logical, save :: sci_arrstats = .false.
  logical, save :: sci_stress = .false.
  logical, save :: sci_colour_clusters = .false.
  logical, save :: post = .false.
  logical, save :: sci_vel_correction = .true.  
  character(len=256), save :: gr_out_file = 'default'
  character(len=256), save :: folder = 'default'
  character(len=256), save :: cpfolder = '.'
  character(len=256), save :: srccpfolder = ''

  logical, save :: restore = .false.
  integer, save :: init_cond = 0
  real(kind=rk), save :: vel_poiseuille_max = 0.0_rk
  integer, save :: sci_pressure_init = 1

  integer, save :: n_iteration       = 100
  integer, save :: n_sci_start       = 0
  integer, save :: n_sci_int         = 100
  integer, save :: n_sci_sur         = 100
  integer, save :: n_sci_od          = 100
  integer, save :: n_sci_wd          = 100
  integer, save :: n_sci_dir         = 100
  integer, save :: n_sci_vel         = 100
  integer, save :: n_sci_flo         = 100
  integer, save :: n_sci_arrows      = 100
  integer, save :: n_sci_velocities  = 100
  integer, save :: n_sci_velocities_od  = 100
  integer, save :: n_sci_velocities_wd  = 100
  integer, save :: n_sci_flux_od  = 100
  integer, save :: n_sci_flux_wd  = 100
  integer, save :: n_sci_rock        = 100
  integer, save :: n_sci_rock_colour = 100
  integer, save :: n_sci_rock_rho_r  = 100
  integer, save :: n_sci_rock_rho_b  = 100
  integer, save :: n_sci_pressure    = 100
  integer, save :: n_sci_fluxz       = 100
  integer, save :: n_sci_massfluxz   = 100
  integer, save :: n_sci_profile     = 100
  !> \c n_sci_profile is the sampling interval for profile data. The
  !> data sampled since the last dump is averaged and dumped directly
  !> after sampling but only if \c mod(nt,n_sci_profile_dump)==0
  !> . This has two consequences: 1) As long as \c n_sci_profile_dump
  !> stays at its default value of \c 1, \c n_sci_profile controls the
  !> sampling and the dumping interval as well. 2) If averaged output
  !> is desired \c n_sci_profile_dump should be set to an integer
  !> multiple of \c n_sci_profile .
  integer, save :: n_sci_profile_dump = 1
  integer, save :: n_sci_arrstats = 100
  integer, save :: n_sci_stress = 100
  integer, save :: n_sci_colour_clusters = 100
  integer, save :: n_sci_colour_clusters_index = 100
  integer, save :: n_sci_start_colour_clusters = 0

  logical, save :: dbg_report_hk_timing = .false.

  integer, save :: stress_minx = -1
  integer, save :: stress_maxx = -1
  real(kind=rk), save :: stress_dd_cutoff = 0.0_rk  ! Droplet density cutoff

  !> maximum number of regions to calculate and dump fluxz data for
  integer, parameter :: fluxz_regions = 10

  !> Human-readable identifier for each fluxz region, default values
  !> are assigned later in case they are not overwritten in the input
  !> file.
  character(len=32), save :: fluxz_name(fluxz_regions) = ""

  !> \name lower boundaries for each fluxz region
  !> \{
  integer, save :: fluxz_xlo(fluxz_regions) = 1
  integer, save :: fluxz_ylo(fluxz_regions) = 1
  integer, save :: fluxz_zlo(fluxz_regions) = 1
  !> \}

  !> \name upper boundaries for each fluxz region
  !> \{
  !> value -1 for x-boundary indicates an unused slot, value 0 the
  !> maximum lattice position
  integer :: fxhidx                !< dummy index for in-line initialization
  integer, save :: fluxz_xhi(fluxz_regions)&
       & = (/0,(-1,fxhidx=2,fluxz_regions)/)
  integer, save :: fluxz_yhi(fluxz_regions) = 0
  integer, save :: fluxz_zhi(fluxz_regions) = 0
  !> \}

  !> maximum number of different intervals over which arrstats data is
  !> accumulated separately
  integer, parameter :: arrstats_intervals = 10
  integer :: asidx              !< dummy index for in-line initialization
  !> dump accumulated arrstats data each time step arrstats data is
  !> sampled and \c mod(nt,n_sci_arrstats_dump)==0. A value of \c 0
  !> disables the respective slot.
  integer, save :: n_sci_arrstats_dump(arrstats_intervals)&
       & = (/1,(0,asidx=2,arrstats_intervals)/)
  !> Human-readable identifier for each arrstats slot, default values
  !> are assigned later in case they are not overwritten in the input
  !> file.
  character(len=32), save :: arrstats_name(arrstats_intervals) = ""

  integer, save :: nr = 16
  real(kind=rk), save :: fr = 1.0_rk, fb = 0.5_rk, fg = 0.1_rk
  integer, save :: fd = 0
  real(kind=rk), save :: fr1 = 0.2_rk, fr2 = 0.3_rk
  real(kind=rk), save :: pr = 0.5_rk, pb = 1.0_rk, pg = 0.1_rk
  integer, save :: pd = 0
  real(kind=rk), save :: qr = 0.5_rk, qb = 0.5_rk, qg = 0.1_rk
  integer, save :: qd = 0
  real(kind=rk), save :: r1 = 10.0_rk !< WARNING Not used anywhere, but is in namelist
  real(kind=rk), save :: m_evp = 0.000_rk
  real(kind=rk), save :: m_evp_gr = 0.000_rk
  real(kind=rk), save :: m_evp_gb = 0.000_rk
  real(kind=rk), save :: m_evp_freq_f = 0.000_rk
  real(kind=rk), save :: m_evp_freq_a = 0.000_rk
  logical, save :: m_evp_set_density = .false.
  logical, save :: in_evp(3) = .false.
  logical, save :: out_evp(3) = .false.
  logical, save :: rock_colour_init = .false.
  real(kind=rk), save :: rock_colour = 0.0_rk
  real(kind=rk), save :: rock_colour_r = 0.0_rk
  real(kind=rk), save :: rock_colour_b = 0.0_rk


  integer, save :: inv_fluid = 0
  integer, save :: inv_type = 0
  real(kind=rk), save :: beta = 1.0_rk !< Dipole temperature

  !added for shifting/cutting droplets
  integer, save :: drop_xshift = 0
  integer, save :: drop_yshift = 0
  integer, save :: drop_zshift = 0
  integer, save :: drop_xcut = 0
  integer, save :: drop_ycut = 0
  integer, save :: drop_zcut = 0

  namelist /VARIABLE_INPUT/ &
       &arrstats_name,sci_arrstats,sci_int,sci_sur,sci_od,sci_wd,sci_dir&
       &,sci_vel,sci_flo,sci_arrows,sci_velocities,sci_rock,sci_rock_colour&
       &,sci_rock_rho_r,sci_rock_rho_b,sci_pressure&
       &,sci_fluxz,sci_massfluxz,sci_profile,sci_stress,post,gr_out_file&
       &,folder,cpfolder,srccpfolder,restore,init_cond,vel_poiseuille_max,sci_pressure_init&
       &,n_iteration,n_sci_arrstats,n_sci_arrstats_dump,n_sci_start,n_sci_int&
       &,n_sci_sur,n_sci_od,n_sci_wd,n_sci_dir,n_sci_vel,n_sci_flo,n_sci_arrows&
       &,n_sci_velocities,n_sci_rock,n_sci_rock_colour,n_sci_rock_rho_r,n_sci_rock_rho_b&
       &,n_sci_pressure,n_sci_fluxz,n_sci_massfluxz&
       &,n_sci_profile,n_sci_profile_dump,n_sci_stress,fluxz_name&
       &,sci_colour_clusters,n_sci_colour_clusters&
       &,n_sci_colour_clusters_index,n_sci_start_colour_clusters,dbg_report_hk_timing&
       &,fluxz_xlo,fluxz_ylo,fluxz_zlo,fluxz_xhi,fluxz_yhi,fluxz_zhi,fr,fb&
       &,fg,fd,fr1,fr2,pr,pb,pg,pd,qr,qb,qg,qd,rock_colour_init,rock_colour,inv_fluid,inv_type&
       &,beta,m_evp,m_evp_gr,m_evp_gb,m_evp_freq_f,m_evp_freq_a,in_evp,out_evp,m_evp_set_density,drop_xshift,drop_yshift,drop_zshift&
       &,drop_xcut,drop_ycut,drop_zcut,stress_minx,stress_maxx,stress_dd_cutoff&
       &,r1,rock_colour_r,rock_colour_b,sci_vel_correction&
       &,sci_velocities_od,sci_velocities_wd,n_sci_velocities_od,n_sci_velocities_wd& 
       &,sci_flux_od,sci_flux_wd,n_sci_flux_od,n_sci_flux_wd  !Deprecated

  ! ======================================================================
  !                                LBE INPUT
  ! ======================================================================

  ! These variables are intrinsic to a given simulation run.

  ! Interaction types
  ! Single Component Multi Phase - pos. intra-comp. SC force
  logical, save :: SCMP = .false. 
  ! Multi Component Multi Phase - pos. intra-comp.-, neg. inter-comp.-SC-force
  logical, save :: MCMP = .false. !< WARNING Not used anywhere, but is in namelist
  logical, save :: MRT = .false.
  ! Declare different collision model types:
  ! BGK - Single relaxation time model (=BGK)
  ! MRT - multiple relaxation time scheme (from Sebastian)
  ! FLUCTUATING_MRT - fluctuating mrt scheme from Duenweg et al (implemented by Philipp)
  integer, parameter :: BGK_id = 1, MRT_id = 2, fluctuating_MRT_id = 3
  integer, save :: collisiontype_id = BGK_id
  character(len=32), save :: collisiontype = 'BGK'
  ! For fluctuating MRT, we also need to know the value k_B*T
  ! where k_B is Boltzmann's constant and T is the temperature of the system. Default
  ! value is 1e-7
  ! This is also used for ELEC temperatures.
  real(kind=rk), save :: kbT = -0.0000001_rk
  ! Temperature used for Carnahan-Starling equation of state
  real(kind=rk), save :: Tcs = 0.3773_rk
  logical, save :: OXFORD = .false. !< WARNING Not used anywhere, but is in namelist
  logical, save :: zeroforceoffset = .true.
  logical, save :: ZFOSdiag = .true.

  ! Molecular masses
  real(kind=rk), save :: amass_r = 1.0_rk, amass_b = 1.0_rk, amass_s = 1.0_rk
  ! Relaxation times
  real(kind=rk), save :: tau_r = 1.0_rk, tau_b = 1.0_rk, tau_s = 1.0_rk, tau_d = 2.0_rk
  real(kind=rk), save :: taubulk_r = 0.84_rk, taubulk_b = 0.84_rk, taubulk_s = 0.84_rk, taubulk_d = 1.68_rk
  real(kind=rk), save :: s02_r = 1.19_rk, s03_r = 1.4_rk, s05_r = 1.2_rk, s11_r = 1.4_rk, s14_r = 1.0_rk, s17_r = 1.98_rk
  real(kind=rk), save :: s03_b = 1.0_rk, s05_b = 1.0_rk, s11_b = 1.0_rk, s14_b = 1.0_rk, s17_b = 1.0_rk
  real(kind=rk), save :: s03_s = 1.0_rk, s05_s = 1.0_rk, s11_s = 1.0_rk, s14_s = 1.0_rk, s17_s = 1.0_rk
  integer, save :: bcsel = 0
  integer, save :: interpolation_order = 1
  integer, save :: interpolation_scheme = 1
  real(kind=rk), save :: acccoef = 2.0_rk
  ! Inverse relaxation times
  real(kind=rk), save :: omega_b, omega_r, omega_s, omega_d ! Not in namelist
  real(kind=rk), save :: omegabulk_b, omegabulk_r, omegabulk_s, omegabulk_d  ! Not in namelist

  ! Interaction strengths
  real(kind=rk), save :: g_br = 0.00_rk , g_rb = 0.00_rk, g_bs = -0.00_rk , g_ss = 0.00_rk
  real(kind=rk), save :: g_rr = 0.00_rk, g_bb = 0.00_rk
  real(kind=rk), save :: g_wr = 0.0_rk, g_wb = 0.0_rk
  real(kind=rk), save :: tau_wr = 1.0_rk, tau_wb = 1.0_rk
  
  logical, save :: g_accn_fluid_r = .false.
  logical, save :: g_accn_fluid_b = .false.
  logical, save :: g_accn_fluid_s = .false.
  real(kind=rk), save :: g_accn = 0.0_rk    ! Acceleration due to gravity
  real(kind=rk), save :: g_accn_x = 0.0_rk  ! Acceleration due to gravity
  real(kind=rk), save :: g_accn_y = 0.0_rk  ! Acceleration due to gravity
  integer, save :: g_accn_min = 0    ! Acceleration due to gravity
  integer, save :: g_accn_max = 0    ! Acceleration due to gravity
  integer, save :: g_accn_min_x = 0    ! Acceleration due to gravity in x
  integer, save :: g_accn_max_x = 0    ! Acceleration due to gravity in x 
  integer, save :: g_accn_min_y = 0    ! Acceleration due to gravity in y
  integer, save :: g_accn_max_y = 0    ! Acceleration due to gravity in y
  integer, save :: nt = 0             ! No of timesteps since last output - not in namelist

  real(kind=rk), save :: perturbation = 0.0_rk ! Perturbation to initial state.

  ! Edge fluid's speed in shear system and shear rate frequency.
  real(kind=rk), save :: shear_u = 0.0_rk
  real(kind=rk), save :: shear_omega = 0.0_rk

  integer, save :: n_checkpoint = 300 ! Number of timesteps b/t checkpoints
  character(len=3), save :: checkpoint_format = 'xdr'
  logical, save :: checkpoint_safe = .true.
  logical, save :: checkpoint_shearsum = .true.
  integer, save :: num_chkp_files = 0

  integer, save :: chk_uid ! Unique identifier for chkpoint files - not in namelist
  character(len=32), save :: restore_string = "t00000000-0000000000"
  integer, save :: n_restore = 0 ! Timestep to restore from - not in namelist

  integer, save :: psifunc = 2    ! Which form of \psi to use in force terms
  integer, save :: bdist = 0
  real(kind=rk), save :: d_0 = 1.0_rk        ! Maximum modulus of dipoles

  character(len=3), save :: dump_format = 'bin'
  logical, save :: xdrfsloppy = .false.
  logical, save :: write_AVS_fld = .false.
  logical, save :: dump_double = .true.

  logical, save :: hdf_use_ibm_largeblock_io = .false.
  logical, save :: hdf_use_independent_io = .false.
  logical, save :: dbg_report_hdf5 = .false.
  logical, save :: dbg_report_hdf5_timing = .false.

  character(len=80), save :: inp_file ! Name of the input file - not in namelist

  ! These variables are set through command-line options
  ! arg_foo_p is set to 1 if arg_foo is defined.
  ! Not in namelist!
  character(len=256), save :: arg_restore_string
  logical, save            :: arg_restore_string_set = .false.
  character(len=256),save  :: arg_input_file
  logical, save            :: arg_input_file_set = .false.
  character(len=256),save  :: arg_input_dfile
  logical, save            :: arg_input_dfile_set = .false.

  ! interval after which sanity check is performed (0 means never)
  integer, save :: n_sanity_check = 1000

#ifdef VELCONTROL 
  ! VELCONTROL is a feedback control to achieve desired setpoint velocity
  ! using Proportional-Integral-Derivative (PID) control
  real(kind=rk), save :: k_p = 1.0_rk            ! Proportional constant
  real(kind=rk), save :: k_i = 0.01_rk           ! Integral constant
  real(kind=rk), save :: k_d = 0.01_rk           ! Derivative constant
  real(kind=rk), save :: u_refx = 0._rk          ! Desired setpoint velocity
  real(kind=rk), save :: u_refy = 0._rk          ! Desired setpoint velocity
  real(kind=rk), save :: u_refz = 0._rk          ! Desired setpoint velocity
  integer, save :: fluid_ramptime = 5000         ! Timesteps to ramp fluid
#endif

  namelist /LBE_INPUT/ &
    SCMP, MRT, collisiontype, kbT, Tcs, ZFOSdiag, &
    amass_r, amass_b, amass_s, &
    tau_r, tau_b, tau_s, tau_d, &
    taubulk_r, taubulk_b, taubulk_s, taubulk_d, s02_r, &
    s03_r, s05_r, s11_r, s14_r, s17_r, &
    s03_b, s05_b, s11_b, s14_b, s17_b, &
    s03_s, s05_s, s11_s, s14_s, s17_s, &
    bcsel, acccoef, &
    g_br, g_bs, g_ss, g_rr, g_bb, g_wr, g_wb, tau_wr, tau_wb, bdist, &
    g_accn_fluid_r, g_accn_fluid_b, g_accn_fluid_s, &
    g_accn, g_accn_x, g_accn_y, g_accn_min, g_accn_max, &
    g_accn_min_x, g_accn_max_x, g_accn_min_y, g_accn_max_y, &
    perturbation, shear_omega,shear_u, zeroforceoffset, &
    n_checkpoint, checkpoint_format, checkpoint_safe, checkpoint_shearsum, &
    restore_string, num_chkp_files, psifunc, d_0, &
    dump_format, xdrfsloppy, write_AVS_fld, dump_double, &
    hdf_use_ibm_largeblock_io, hdf_use_independent_io, dbg_report_hdf5, &
    dbg_report_hdf5_timing, n_sanity_check, &
    interpolation_order, interpolation_scheme, &
#ifdef VELCONTROL
     k_p, k_i, k_d, u_refx, u_refy, u_refz,&
#endif
    MCMP, OXFORD ! Deprecated
  contains

  ! ======================================================================
  !                        SUBROUTINES / FUNCTIONS
  ! ======================================================================

  !> Get relaxation times at a given lattice site
  !>
  !> If VARTAU is not defined, the return value is merely the global value tau_*.
  !> If VARTAU is defined, the function returns the relaxation time stored locally in N(X, Y, Z)%tau_*.
  !> For efficiency reasons, there is one wrapper function for each component (r, b, s) with repeated code.

  function get_tau_r(N, X, Y, Z)
    implicit none
    real(kind=rk) :: get_tau_r !< return value
    type(lbe_site), dimension(0:, 0:, 0:), intent(in) :: N !< lattice
    integer :: X, Y, Z !< lattice coordinates

#ifdef VARTAU
    get_tau_r = N(X, Y, Z)%taupos_r
#else
    get_tau_r = tau_r
#endif
  end function get_tau_r

#ifndef SINGLEFLUID
  function get_tau_b(N, X, Y, Z)
    implicit none
    real(kind=rk) :: get_tau_b !< return value
    type(lbe_site), dimension(0:, 0:, 0:), intent(in) :: N !< lattice
    integer :: X, Y, Z !< lattice coordinates

#ifdef VARTAU
    get_tau_b = N(X, Y, Z)%taupos_b
#else
    get_tau_b = tau_b
#endif
  end function get_tau_b
#endif

#ifndef NOSURFACTANT
  function get_tau_s(N, X, Y, Z)
    implicit none
    real(kind=rk) :: get_tau_s !< return value
    type(lbe_site), dimension(0:, 0:, 0:), intent(in) :: N !< lattice
    integer :: X, Y, Z !< lattice coordinates

#ifdef VARTAU
    get_tau_s = N(X, Y, Z)%taupos_s
#else
    get_tau_s = tau_s
#endif
  end function get_tau_s
#endif

  !> Get inverse relaxation times at a given lattice site
  !>
  !> If VARTAU is not defined, the return value is merely the global value omega_*.
  !> If VARTAU is defined, the function returns the inverse of the relaxation time stored locally in N(X, Y, Z)%tau_*.
  !> For efficiency reasons, there is one wrapper function for each component (r, b, s) with repeated code.

  function get_omega_r(N, X, Y, Z)
    implicit none
    real(kind=rk) :: get_omega_r !< return value
    type(lbe_site), dimension(0:, 0:, 0:), intent(in) :: N !< lattice
    integer :: X, Y, Z !< lattice coordinates

#ifdef VARTAU
    get_omega_r = 1.0_rk / N(X, Y, Z)%taupos_r
#else
    get_omega_r = omega_r
#endif
  end function get_omega_r

#ifndef SINGLEFLUID
  function get_omega_b(N, X, Y, Z)
    implicit none
    real(kind=rk) :: get_omega_b !< return value
    type(lbe_site), dimension(0:, 0:, 0:), intent(in) :: N !< lattice
    integer :: X, Y, Z !< lattice coordinates

#ifdef VARTAU
    get_omega_b = 1.0_rk / N(X, Y, Z)%taupos_b
#else
    get_omega_b = omega_b
#endif
  end function get_omega_b
#endif

#ifndef NOSURFACTANT
  function get_omega_s(N, X, Y, Z)
    implicit none
    real(kind=rk) :: get_omega_s !< return value
    type(lbe_site), dimension(0:, 0:, 0:), intent(in) :: N !< lattice
    integer :: X, Y, Z !< lattice coordinates

#ifdef VARTAU
    get_omega_s = 1.0_rk / N(X, Y, Z)%taupos_s
#else
    get_omega_s = omega_s
#endif
  end function get_omega_s
#endif

end module lbe_parms_module
