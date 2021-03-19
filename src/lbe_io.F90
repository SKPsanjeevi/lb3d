#include "lbe.h"

!> Loading and saving of files
!>
!> Most of it was stolen directly from the corresponding code
!> in ME3D, since the idea is to have compatible input and output files.
module lbe_io_module

  use hoshen_kopelman_module
  use lbe_analysis_module, only: cached_avg_total_velocity,total_mass
  use lbe_collision_module, only: lbe_calculate_sc_forces, md_calculate_sc_forces
  use lbe_globals_module
  use lbe_helper_module, only: check_dump_now,check_nan,cross_product&
       &,interpolate,lbe_count_sites,local_coordinates,massflow,n_sites&
       &,n_sites_fluid,n_sites_particle,n_sites_rock,norm,unit_vector&
       &,is_restoring, is_wall, is_fluid, is_colloid, n_sites_surface
  use lbe_io_helper_module
  use lbe_log_module
  use lbe_parms_module
  use lbe_parallel_module
  use lbe_types_module, only: lbe_site
  

#ifdef USEHDF
  use lbe_io_hdf5_module
#endif
#ifdef USEXDRF
  use lbe_io_xdrf_module
#endif

#ifdef MD
  use lbe_md_globals_module, only: atompnt,communicate_velocities,interaction&
       &,list,nlocal,nother,P,uid2i
  use map_module, only: Mii_map
use lbe_md_helper_module, only: error_md
#endif

  implicit none
  include 'mpif.h'
  private

  public lbe_get_lbe_input, lbe_parse_arguments, lbe_get_fixed_input, lbe_get_variable_input
  public lbe_report_flags, lbe_detect_flags
  public lbe_write_topology, lbe_print_site_counts, postprocess, dump_data, lbe_sanity_check
  public dump_scalar, dump_iscalar, dump_vector, dump_fluxz, dump_profile, check_dump_now
  public read_rock_all_par, corrected_velocity, corrected_velocity_new

  !> \name declarations for PROFILE output
  !> \{

  !> type representing one point of a profile containing data accumulated
  !> within the layer perpendicular to the profile direction
  type layer
     !> total fluid velocity
     real(kind=rk) :: vf(3)
     !> total mass flow
     real(kind=rk) ::  mf(3)
     !> total fluid density
     real(kind=rk) :: rhof
#ifndef SINGLEFLUID
     !> velocity of red fluid
     real(kind=rk) :: v_r(3)
     !> mass flow of red fluid
     real(kind=rk) :: m_r(3)
     !> density of red fluid
     real(kind=rk) :: rho_r
     !> velocity of blue fluid
     real(kind=rk) :: v_b(3)
     !> mass flow of blue fluid
     real(kind=rk) :: m_b(3)
     !> density of blue fluid
     real(kind=rk) :: rho_b
#endif
     !> number density of fluid sites
     real(kind=rk) :: nfs
     !> number density of solid rock sites
     real(kind=rk) :: nrs
#ifdef MD
     !> particle velocity
     real(kind=rk) :: vp(3)
     !> number densities of particle centers
     real(kind=rk) :: npc
     !> number density of particle sites
     real(kind=rk) :: nps
#endif
  end type layer

  !> type representing a complete profile along a specific direction
  !>
  !> \note maybe \c l should be pointer instead of allocatable, but it
  !> seems to work
  type profile
     character(len=5) :: name   !< distinct part of file names
     !> unit vectors in profile and averaging directions
     real(kind=rk) :: d(3,3)
     real(kind=rk) :: vd(3,3)    !< directions to project velocities at
     integer :: s(3)          !< starting point of profile
     integer :: size(3) !< dimensions in profile and averaging directions
     type(layer),allocatable :: l(:) !< data for each point of the profile
  end type profile

  !> profiles for PROFILE output
  type(profile),allocatable,save,private :: prf(:)
#ifdef MD
  !> particle center velocity at every site
  real(kind=rk),allocatable,save,private :: Nvp(:,:,:,:)
  !> particle number density at every site
  real(kind=rk),allocatable,save,private :: Nnpc(:,:,:)
#endif
  integer,save,private :: layer_mpitype !< custom mpi data type for layer
  integer,save,private :: sum_layer_mpiop !< custom mpi reduction operation

  !> summation buffer, used only by root process
  type(profile),allocatable,save,private :: prfsum(:)
  !> number of samples taken (also included in output; used only by
  !> root process)
  integer,save :: n_samples
  !> \}

  character(len=16), parameter :: input_source = '.input-file' !< Name of the .input-file

#ifdef BUGGYIFORT11
  !> bugfix for ifort 11
  !>
  !> There is a bug in Intel Fortran Compiler 11.0 and 11.1 (but
  !> probably not in versions 9.1 and 10.1) which leads to erroneous
  !> optimization of
  !> fluid_velocity_and_density_and_site_occupation(). Doing
  !> meaningless operations with the variable below seems to prevent at
  !> least ifort 11.1 from showing this behavior. Essential for this
  !> bug is that
  !> \li the code is compiled with the options -O3 and -ipo and that
  !> \li somewhere else in the code there is a call to a not
  !>     user-provided subroutine with another subroutine as
  !>     argument. It is not necessary that this subroutine call will
  !>     be executed. Here mpi_op_create() in lbe_io_module is the
  !>     call that causes the trouble.
  integer,save,public :: ifort_11_bug_dummy = 0
#endif

contains

!> This routine is called for all values of init_cond apart from
!> \c INIT_HUDONG; it reads the LBE-specific namelist from the input file and
!> broadcasts the data to all processors.
subroutine lbe_get_lbe_input()
  implicit none
  integer :: ierror, itmp
  
  if (myrankc == 0) then
    call log_msg_hdr("Reading LBE input")
    open(UNIT = input_file_unit, FILE = inp_file, STATUS = 'UNKNOWN')
    read(UNIT = input_file_unit, NML = lbe_input)
    close(unit=input_file_unit)
   
    if ( arg_input_dfile_set ) then
      call log_msg("  Getting differential input...")
      open(UNIT = input_dfile_unit, FILE = arg_input_dfile, STATUS = 'UNKNOWN')
      read(UNIT = input_dfile_unit, NML = lbe_input, IOSTAT = ierror)
      if (ierror .ne. 0) then
        call log_msg("    WARNING: Differential namelist not found or errors encountered.")
      endif
      close(UNIT = input_dfile_unit)
      call log_ws()
    end if

  end if
  
  call MPI_Bcast(amass_r,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(tau_r,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(taubulk_r,1,MPI_REAL8,0,comm_cart,ierror)

  call MPI_Bcast(amass_b,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(tau_b,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(taubulk_b,1,MPI_REAL8,0,comm_cart,ierror)

  call MPI_Bcast(amass_s,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(tau_s,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(taubulk_s,1,MPI_REAL8,0,comm_cart,ierror)

  call MPI_Bcast(tau_d,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(taubulk_d,1,MPI_REAL8,0,comm_cart,ierror)

  call MPI_Bcast(g_rr,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(g_br,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(g_bb,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(g_bs,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(g_ss,1,MPI_REAL8,0,comm_cart,ierror)

  call MPI_Bcast(g_wr,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(g_wb,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(tau_wr,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(tau_wb,1,MPI_REAL8,0,comm_cart,ierror)

  call MPI_Bcast(d_0,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(perturbation,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(acccoef,1,MPI_REAL8,0,comm_cart,ierror)

  call MPI_Bcast(g_accn  ,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(g_accn_x,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(g_accn_y,1,MPI_REAL8,0,comm_cart,ierror)

  call MPI_Bcast(shear_u,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(shear_omega,1,MPI_REAL8,0,comm_cart,ierror)

  call MPI_Bcast(s03_r,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(s05_r,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(s11_r,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(s14_r,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(s17_r,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(s03_b,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(s05_b,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(s11_b,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(s14_b,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(s17_b,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(s03_s,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(s05_s,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(s11_s,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(s14_s,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(s17_s,1,MPI_REAL8,0,comm_cart,ierror)

  ! Calculate some derived quantities
  omega_b = 1.0_8/tau_b
  omega_r = 1.0_8/tau_r
  omega_s = 1.0_8/tau_s
  omega_d = 1.0_8/tau_d
  omegabulk_b = 1.0_8/taubulk_b
  omegabulk_r = 1.0_8/taubulk_r
  omegabulk_s = 1.0_8/taubulk_s
  omegabulk_d = 1.0_8/taubulk_d

  ! Broadcast integer, logical and character data as well.
  call MPI_Bcast(n_sanity_check,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(interpolation_order,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(n_checkpoint,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(num_chkp_files,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(bcsel,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(bdist,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(psifunc,1,MPI_INTEGER,0,comm_cart,ierror)

  call MPI_Bcast(g_accn_min,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(g_accn_max,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(g_accn_min_x,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(g_accn_max_x,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(g_accn_min_y,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(g_accn_max_y,1,MPI_INTEGER,0,comm_cart,ierror)

#ifdef VELCONTROL
  call MPI_Bcast(u_refx,1,LBE_REAL,0,comm_cart,ierror)
  call MPI_Bcast(u_refy,1,LBE_REAL,0,comm_cart,ierror)
  call MPI_Bcast(u_refz,1,LBE_REAL,0,comm_cart,ierror)
  call MPI_Bcast(k_p,1,LBE_REAL,0,comm_cart,ierror)
  call MPI_Bcast(k_i,1,LBE_REAL,0,comm_cart,ierror)
  call MPI_Bcast(k_d,1,LBE_REAL,0,comm_cart,ierror)
  call MPI_Bcast(fluid_ramptime,1,MPI_INTEGER,0,comm_cart,ierror)
#endif

  call MPI_Bcast(g_accn_fluid_r,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(g_accn_fluid_b,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(g_accn_fluid_s,1,MPI_LOGICAL,0,comm_cart,ierror)
 
  call MPI_Bcast(MRT,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(collisiontype,32,MPI_CHARACTER,0,comm_cart,ierror)
  call MPI_Bcast(kbT,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(Tcs,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(zeroforceoffset,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(ZFOSdiag,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(write_AVS_fld,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(dump_double,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(hdf_use_ibm_largeblock_io,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(hdf_use_independent_io,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(dbg_report_hdf5,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(dbg_report_hdf5_timing,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(SCMP,1,MPI_LOGICAL,0,comm_cart,ierror)

  call MPI_Bcast(dump_format,3,MPI_CHARACTER,0,comm_cart,ierror)
  call MPI_Bcast(xdrfsloppy,1,MPI_LOGICAL,0,comm_cart,ierror)

  call MPI_Bcast(checkpoint_format,3,MPI_CHARACTER,0,comm_cart,ierror)
  call MPI_Bcast(checkpoint_safe,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(checkpoint_shearsum,1,MPI_LOGICAL,0,comm_cart,ierror)

  halo_extent = interpolation_order

  ! Calculate some derived quantities

  if ( g_accn_max == 0 ) g_accn_max = tnz
  if ( g_accn_max_x == 0 ) g_accn_max_x = tnx
  if ( g_accn_max_y == 0 ) g_accn_max_y = tny

  ! Hack to keep allowing the old-style collision selection.
  if (MRT .eqv. .true. ) then
    collisiontype = "MRT"
  end if
  
  call log_ws()
  write(msgstr,"('collisiontype   = <',A,'>')") trim(collisiontype)
  call log_msg(msgstr)

  select case(collisiontype)
  case("BGK")
    call log_msg("  -> Using BGK collision model.")
    collisiontype_id = BGK_id
  case("MRT")
    call log_msg("  -> Using MRT collision model.")
    collisiontype_id = MRT_id
  case("FLUCTUATING_MRT")
    call log_msg("  -> Using fluctuating MRT collision model.")
    collisiontype_id = fluctuating_MRT_id
  case default
    write(msgstr,"('FATAL ERROR: Unknown collisiontype <',A,'>. Aborting...')") trim(collisiontype)
    call log_msg(msgstr)
    call Abend
  end select
  call log_ws()

  ! Currently, MRT and VARTAU are not compatible.
#ifdef VARTAU
  if (collisiontype_id .eq. MRT_id) then
    call log_msg("FATAL ERROR: Compiler flag VARTAU and collision type MRT are not compatible. Aborting...")
    call Abend
  end if
#endif

  write(msgstr,"('kbT             = ',F16.10)") kbT
  call log_msg(msgstr)
   write(msgstr,"('Tcs             = ',F16.10)") Tcs
  call log_msg(msgstr)
  write(msgstr,"('zeroforceoffset = ',L1)") zeroforceoffset
  call log_msg(msgstr)
  write(msgstr,"('write_AVS_fld   = ',L1)") write_AVS_fld
  call log_msg(msgstr)
  write(msgstr,"('dump_double     = ',L1)") dump_double
  call log_msg(msgstr)
  write(msgstr,"('SCMP            = ',L1)") SCMP
  call log_msg(msgstr)
  write(msgstr,"('dump_format     = <',A,'>')") trim(dump_format)
  call log_msg(msgstr)
#ifdef USEHDF
  write(msgstr,"('hdf_use_ibm_largeblock_io = ',L1)") hdf_use_ibm_largeblock_io
  call log_msg(msgstr)
  write(msgstr,"('hdf_use_independent_io    = ',L1)") hdf_use_independent_io
  call log_msg(msgstr)
  if (dbg_report_hdf5) then
    write(msgstr,"('dbg_report_hdf5           = ',L1)") dbg_report_hdf5
    call log_msg(msgstr)
  end if
  if (dbg_report_hdf5_timing) then
    write(msgstr,"('dbg_report_hdf5_timing    = ',L1)") dbg_report_hdf5_timing
    call log_msg(msgstr)
  end if
#endif
  write(msgstr,"('xdrfsloppy      = ',L1)") xdrfsloppy
  call log_msg(msgstr)
  call log_ws()

  write(msgstr,"('n_sanity_check      = ',I0)") n_sanity_check
  call log_msg(msgstr)
  write(msgstr,"('n_checkpoint        = ',I0)") n_checkpoint
  call log_msg(msgstr)
  write(msgstr,"('checkpoint_format   = <',A,'>')") trim(checkpoint_format)
  call log_msg(msgstr)
  write(msgstr,"('checkpoint_safe     = ',L1)") checkpoint_safe
  call log_msg(msgstr)
  write(msgstr,"('checkpoint_shearsum = ',L1)") checkpoint_shearsum
  call log_msg(msgstr)
  write(msgstr,"('num_chkp_files      = ',I0)") num_chkp_files
  call log_msg(msgstr)
  if ( is_restoring() ) then
    write(msgstr,"('restore_string      = <',A,'>')") trim(restore_string)
    call log_msg(msgstr)
    write(msgstr,"('  -> n_restore      = ',I0)") n_restore
    call log_msg(msgstr)
  end if
  call log_ws()

  write(msgstr,"('interpolation_order        = ',I0)") interpolation_order
  call log_msg(msgstr)
  write(msgstr,"('bcsel        = ',I0)") bcsel
  call log_msg(msgstr)
  if ( bcsel .eq. 1 ) then
    write(msgstr,"('acccoef      = ',F16.10)") acccoef
    call log_msg(msgstr)
  end if

#if defined (MD) && !defined (INTERPOLATEDBB)
   if ( bcsel .eq. 2 ) then
   call error_md("bcsel='2' "&
               &//"does not yet support particle---use "&
               &//"bscel='0' for particles simulation!")
  end if
#endif
#ifdef VELCONTROL
  write(msgstr,"('Reference velocity      = ',3F16.10)") u_refx, u_refy, u_refz
  call log_msg(msgstr)
  write(msgstr,"('PID constants (k_p, k_i, k_d)  = ',3F16.10)") k_p, k_i, k_d
  call log_msg(msgstr)
#endif


  write(msgstr,"('bdist        = ',I0)") bdist
  call log_msg(msgstr)
  write(msgstr,"('psifunc      = ',I0)") psifunc
  call log_msg(msgstr)
  call log_ws()

  write(msgstr,"('amass_r      = ',F16.10)") amass_r
  call log_msg(msgstr)
  write(msgstr,"('tau_r        = ',F16.10)") tau_r
  call log_msg(msgstr)
#ifndef SINGLEFLUID
  write(msgstr,"('amass_b      = ',F16.10)") amass_b
  call log_msg(msgstr)
  write(msgstr,"('tau_b        = ',F16.10)") tau_b
  call log_msg(msgstr)
#endif
#ifndef NOSURFACTANT
  write(msgstr,"('amass_s      = ',F16.10)") amass_s
  call log_msg(msgstr)
  write(msgstr,"('tau_s        = ',F16.10)") tau_s
  call log_msg(msgstr)
  write(msgstr,"('tau_d        = ',F16.10)") tau_d
  call log_msg(msgstr)
  call log_ws()
#endif
  write(msgstr,"('g_rr         = ',F16.10)") g_rr
  call log_msg(msgstr)
#ifndef SINGLEFLUID
  write(msgstr,"('g_br         = ',F16.10)") g_br
  call log_msg(msgstr)
  write(msgstr,"('g_bb         = ',F16.10)") g_bb
  call log_msg(msgstr)
#ifndef NOSURFACTANT
  write(msgstr,"('g_bs         = ',F16.10)") g_bs
  call log_msg(msgstr)
  write(msgstr,"('g_ss         = ',F16.10)") g_ss
  call log_msg(msgstr)
#endif
#endif

  write(msgstr,"('g_wr         = ',F16.10)") g_wr
  call log_msg(msgstr)
  write(msgstr,"('g_wb         = ',F16.10)") g_wb
  call log_msg(msgstr)
  call log_ws()
  write(msgstr,"('tau_wr       = ',F16.10)") tau_wr
  call log_msg(msgstr)
  write(msgstr,"('tau_wb       = ',F16.10)") tau_wb
  call log_msg(msgstr)
  call log_ws()

#ifndef NOSURFACTANT
  write(msgstr,"('d_0          = ',F16.10)") d_0
  call log_msg(msgstr)
#endif
  write(msgstr,"('perturbation = ',F16.10)") perturbation
  call log_msg(msgstr)
  call log_ws()

  write(msgstr,"('shear_u      = ',F16.10)") shear_u
  call log_msg(msgstr)
  write(msgstr,"('shear_omega  = ',F16.10)") shear_omega
  call log_msg(msgstr)
  call log_ws()


  write(msgstr,"('g_accn_fluid_r = ',L1)") g_accn_fluid_r
  call log_msg(msgstr)
  write(msgstr,"('g_accn_fluid_b = ',L1)") g_accn_fluid_b
  call log_msg(msgstr)
  write(msgstr,"('g_accn_fluid_s = ',L1)") g_accn_fluid_s
  call log_msg(msgstr)
  write(msgstr,"('g_accn         = ',F16.10)") g_accn
  call log_msg(msgstr)
  write(msgstr,"('g_accn_x       = ',F16.10)") g_accn_x
  call log_msg(msgstr)
  write(msgstr,"('g_accn_y       = ',F16.10)") g_accn_y
  call log_msg(msgstr)
  write(msgstr,"('g_accn_min     = ',I0)") g_accn_min
  call log_msg(msgstr)
  write(msgstr,"('g_accn_min_x   = ',I0)") g_accn_min_x
  call log_msg(msgstr)
  write(msgstr,"('g_accn_min_y   = ',I0)") g_accn_min_y
  call log_msg(msgstr)
  write(msgstr,"('g_accn_max     = ',I0)") g_accn_max
  call log_msg(msgstr)
  write(msgstr,"('g_accn_max_x   = ',I0)") g_accn_max_x
  call log_msg(msgstr)
  write(msgstr,"('g_accn_max_y   = ',I0)") g_accn_max_y
  call log_msg(msgstr)

  if ( collisiontype_id .eq. MRT_id .or. collisiontype_id .eq. fluctuating_MRT_id ) then
    call log_ws()
    write(msgstr,"('taubulk_r    = ',F16.10)") taubulk_r
    call log_msg(msgstr)
#ifndef SINGLEFLUID
    write(msgstr,"('taubulk_b    = ',F16.10)") taubulk_b
    call log_msg(msgstr)
#endif
#ifndef NOSURFACTANT
    write(msgstr,"('taubulk_s    = ',F16.10)") taubulk_s
    call log_msg(msgstr)
    write(msgstr,"('taubulk_d    = ',F16.10)") taubulk_d
    call log_msg(msgstr)
#endif
    write(msgstr,"('s03_r        = ',F16.10)") s03_r
    call log_msg(msgstr)
    write(msgstr,"('s05_r        = ',F16.10)") s05_r
    call log_msg(msgstr)
    write(msgstr,"('s11_r        = ',F16.10)") s11_r
    call log_msg(msgstr)
    write(msgstr,"('s14_r        = ',F16.10)") s14_r
    call log_msg(msgstr)
    write(msgstr,"('s17_r        = ',F16.10)") s17_r
    call log_msg(msgstr)
#ifndef SINGLEFLUID
    write(msgstr,"('s03_b        = ',F16.10)") s03_b
    call log_msg(msgstr)
    write(msgstr,"('s05_b        = ',F16.10)") s05_b
    call log_msg(msgstr)
    write(msgstr,"('s11_b        = ',F16.10)") s11_b
    call log_msg(msgstr)
    write(msgstr,"('s14_b        = ',F16.10)") s14_b
    call log_msg(msgstr)
    write(msgstr,"('s17_b        = ',F16.10)") s17_b
    call log_msg(msgstr)
#endif
#ifndef NOSURFACTANT
    write(msgstr,"('s03_s        = ',F16.10)") s03_s
    call log_msg(msgstr)
    write(msgstr,"('s05_s        = ',F16.10)") s05_s
    call log_msg(msgstr)
    write(msgstr,"('s11_s        = ',F16.10)") s11_s
    call log_msg(msgstr)
    write(msgstr,"('s14_s        = ',F16.10)") s14_s
    call log_msg(msgstr)
    write(msgstr,"('s17_s        = ',F16.10)") s17_s
    call log_msg(msgstr)
#endif

  end if

#ifdef FASTBDIST2
  if (bdist/=2) call error('lb3d was compiled with FASTBDIST2 but bdist/=2'&
       &//'---set bdist=2 or recompile withou FASTBDIST2!')
#endif

  call log_msg_ws("Read and broadcast LBE input.")
end subroutine lbe_get_lbe_input

!> Sets the \c inp_file variable to contain the name of the input file.
!> If the -f option is passed on the command line, then
!> this is used. Otherwise, .input-file is checked -- if it contains the
!> word "INTERACTIVE", then the input file name is read from stdin.
!> Otherwise, the file named in .input-file is read.
subroutine lbe_define_inp_file()
  logical :: inp_files_p
  ! assume only called by rank 0
  if ( arg_input_file_set ) then
    ! Take it from command line
    inp_file = arg_input_file
  else
    inquire(FILE = input_source, EXIST = inp_files_p)
    if ( inp_files_p ) then
      open(UNIT = 10, FILE = input_source, STATUS = 'UNKNOWN')
      read(UNIT = 10, FMT = '(A)') inp_file
      close(UNIT = 10)
      if ( index(inp_file,'INTERACTIVE') .gt. 0 ) then
        print *,'Input file ?'
        read(*,*) inp_file
      end if
    else
      call log_msg("FATAL ERROR: Could not open .input-file. Aborting...")
      call Abend
    end if
  endif
  inp_file = trim(inp_file)
  call log_msg("Variable inp_file = <"//trim(inp_file)//">.")
end subroutine lbe_define_inp_file

!> Reads the input file, and parses the \c /FIXED_INPUT/ namelist - see
!> the User's Guide for a description of the variables.
!>
!> Rank zero does the reading, and broadcasts the values to all other CPUs.
subroutine lbe_get_fixed_input()
  implicit none
  integer                :: ierror

  ! Get fixed input:
  ! Read this data only on processor 0 and then broadcast it to the other processors.
  ! Processor reads in .input-file from where it expects to read the main file

  ! Use myrankw instead of myrankc here because the cartesian grid
  ! hasn't been set up yet
  if (myrankw == 0) then
    call log_msg_hdr("Reading fixed input")
    call lbe_define_inp_file()
    open (unit=input_file_unit,file=inp_file,status='UNKNOWN')
    read (unit=input_file_unit,nml=FIXED_INPUT)
    close (unit=input_file_unit)

    if ( arg_input_dfile_set ) then
      call log_msg("  Getting differential input...")
      open(UNIT = input_dfile_unit, FILE = arg_input_dfile, STATUS = 'UNKNOWN')
      read(UNIT = input_dfile_unit, NML = FIXED_INPUT, IOSTAT = ierror)
      if (ierror .ne. 0) then
        call log_msg("    WARNING: Differential namelist not found or errors encountered.")
      end if
      close(UNIT = input_dfile_unit)
      call log_ws()
    end if

    write(msgstr,"('nx = ',i0,', ny = ',i0,', nz = ', i0)") nx, ny, nz
    call log_msg(msgstr)
    write(msgstr,"('seed                = ',i0)") seed
    call log_msg(msgstr)
    write(msgstr,"('boundary_cond       = ',i0)") boundary_cond
    if (boundary_cond .gt. 0) then
      call log_msg(msgstr)
      write(msgstr,"('  boundary_width = ',i0)") boundary_width
    end if
    call log_msg(msgstr)
    write(msgstr,"('boundary            = <',A,'>')") trim(boundary)
    call log_msg(msgstr)
    write(msgstr,"('obs_file            = <',A,'>')") trim(obs_file)
    call log_msg(msgstr)
    if (trim(obs_file) .ne. 'empty.dat') then
      write(msgstr,"('obs_folder = <',A,'>')") trim(obs_folder)
      call log_msg(msgstr)
      write(msgstr,"('obs_rotation = <',A,'>')") trim(obs_rotation)
      call log_msg(msgstr)
    end if
    write(msgstr,"('gw_rock_colour   = ',L1)") rock_colour_double
    call log_msg(msgstr)
    write(msgstr,"('rock_colour_double   = ',L1)") rock_colour_double
    call log_msg(msgstr)
    write(msgstr,"('g_value_double   = ',L1)") g_value_double
    call log_msg(msgstr)
    write(msgstr,"('gw_double_wet   = ',L1)") gw_double_wet
    call log_msg(msgstr)
	 write(msgstr,"('gw_wet_global  = ',L1)") gw_wet_global
    call log_msg(msgstr)
    if (rock_colour_double .or. g_value_double .or. gw_double_wet) then
    ! red rock file
    write(msgstr,"('obs_file_r            = <',A,'>')") trim(obs_file_r)
      call log_msg(msgstr)
    if (trim(obs_file_r) .ne. 'empty.dat') then
        write(msgstr,"('obs_folder = <',A,'>')") trim(obs_folder)
        call log_msg(msgstr)
    end if
    ! blue rock file
    write(msgstr,"('obs_file_b            = <',A,'>')") trim(obs_file_b)
        call log_msg(msgstr)
    if (trim(obs_file_b) .ne. 'empty.dat') then
          write(msgstr,"('obs_folder = <',A,'>')") trim(obs_folder)
          call log_msg(msgstr)
    end if
    end if 
    write(msgstr,"('force               = <',A,'>')") trim(force)
    call log_msg(msgstr)
    write(msgstr,"('galilean_stabilizer = ',L1)") galilean_stabilizer
    call log_msg(msgstr)
    write(msgstr,"('mass_scaler         = ',L1)") mass_scaler
    call log_msg(msgstr)
    write(msgstr,"('cdx = ',i0,', cdy = ',i0,', cdz = ', i0)") cdx, cdy, cdz
    call log_msg(msgstr)
    if ( dbg_report_topology ) then
      write(msgstr,"('dbg_report_topology = ',L1)") dbg_report_topology
      call log_msg(msgstr)
    end if
#ifdef DEBUG_MPI
    write(msgstr,"('dbg_n_start_mpi_debug = ',I0)") dbg_n_start_mpi_debug
    call log_msg(msgstr)
#endif
  end if

  ! note: at this point, comm_cart is still set to MPI_COMM_WORLD
  call MPI_Bcast(inp_file,80,MPI_CHARACTER,0,comm_cart,ierror)

  call MPI_Bcast(nx,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(ny,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(nz,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(seed,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(obs_file,256,MPI_CHARACTER,0,comm_cart,ierror)
  call MPI_Bcast(gw_rock_colour,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(rock_colour_double,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(g_value_double,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(gw_double_wet,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(gw_wet_global,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(obs_file_r,256,MPI_CHARACTER,0,comm_cart,ierror)
  call MPI_Bcast(obs_file_b,256,MPI_CHARACTER,0,comm_cart,ierror)
  call MPI_Bcast(obs_folder,256,MPI_CHARACTER,0,comm_cart,ierror)
  call MPI_Bcast(obs_rotation,3,MPI_CHARACTER,0,comm_cart,ierror)
  call MPI_Bcast(boundary_cond,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(boundary_width,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(boundary,64,MPI_CHARACTER,0,comm_cart,ierror)
  call MPI_Bcast(force,64,MPI_CHARACTER,0,comm_cart,ierror)
  call MPI_Bcast(galilean_stabilizer,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(mass_scaler,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(cdx,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(cdy,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(cdz,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(dbg_report_topology,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(dbg_n_start_mpi_debug,1,MPI_INTEGER,0,comm_cart,ierror)

  ! Set the suggested MPI Cartesian dimensions
  cdims(1) = cdx
  cdims(2) = cdy
  cdims(3) = cdz

  call log_msg_ws("Read and broadcast fixed input.")

end subroutine lbe_get_fixed_input

!> Reads any arguments passed on the command-line, and prints them out.
subroutine lbe_parse_arguments()
  implicit none
  integer :: iargc
  integer :: i, argc
  character(len=1024) :: argbuf

  call log_msg_hdr("Parsing command line arguments")
  argc = iargc() ! Number of arguments on command line.
  if (argc == 0) then
    call log_msg("No command line arguments supplied")
  end if
  i = 1
  do
    if (i .gt. argc) exit
    call getarg(i, argbuf)

    ! -f <input-file name>
    if ("-f" == argbuf) then
      if (i == argc) then
        call log_msg("Need filename after -f, ignoring")
        call Abend
      end if
      i = i+1
      call getarg(i,argbuf)
      arg_input_file = trim(argbuf)
      arg_input_file_set = .true.
      write(msgstr,"('Found -f ',A)") arg_input_file
      call log_msg(msgstr)
    end if

    ! -d <diff input-file name>
    if ("-d" == argbuf) then
      if (i == argc) then
        call log_msg("Need filename after -d, ignoring")
        call Abend
      end if
      i = i+1
      call getarg(i,argbuf)
      arg_input_dfile = trim(argbuf)
      arg_input_dfile_set = .true.
      write(msgstr,"('Found -d ',A)") arg_input_dfile
      call log_msg(msgstr)
    end if

    ! -r <restore-string>
    if ("-r" == argbuf) then
      if (i == argc) then
        call log_msg("Need restore-string after -r, ignoring")
        call Abend
      end if
      i = i+1
      call getarg(i,argbuf)
      arg_restore_string = trim(argbuf)
      arg_restore_string_set = .true.
      write(msgstr,"('Found -r ',A)") arg_restore_string
      call log_msg(msgstr)
    end if
    i = i+1
  end do
end subroutine lbe_parse_arguments

!> Reads the input file, and parses the \c /VARIABLE_INPUT/ namelist -
!> see the User's Guide for a description of the variables.
!>
!> Rank zero does the reading, and broadcasts the values to all other
!> CPUs.
subroutine lbe_get_variable_input()
  implicit none
  integer :: i,ierror,eof_err

  if (myrankc == 0) then ! Let only processor 0 read from file.
    call log_msg_hdr("Reading variable input")
    open (unit=input_file_unit,file=inp_file,status='UNKNOWN')
    read (unit=input_file_unit,nml=VARIABLE_INPUT,iostat=eof_err)
    if(eof_err .lt. 0) then
      call log_msg("FATAL ERROR: End of input file encountered, aborting...")
      close (unit=input_file_unit)
      call Abend
    else if(eof_err .gt. 0) then
      call log_msg("FATAL ERROR: Error in input file encountered, aborting...")
      close (unit=input_file_unit)
      call Abend
    end if
    close(UNIT = input_file_unit)

    if ( arg_input_dfile_set ) then
      call log_msg("  Getting differential input...")
      open(UNIT = input_dfile_unit, FILE = arg_input_dfile, STATUS = 'UNKNOWN')
      read(UNIT = input_dfile_unit, NML = VARIABLE_INPUT, IOSTAT = ierror)
      if (ierror .ne. 0) then
        call log_msg("    WARNING: Differential namelist not found or errors encountered.")
      end if
      close(UNIT = input_dfile_unit)
      call log_ws()
    end if

    if ( arg_restore_string_set ) then
      restore = .true.
    end if

    write(msgstr,"('n_iteration    = ',I0)") n_iteration
    call log_msg(msgstr)
    write(msgstr,"('n_sci_start    = ',I0)") n_sci_start
    call log_msg(msgstr)

    call log_ws()
    write(msgstr,"('sci_od         = ',L1,', n_sci_od           = ',I0)") sci_od, n_sci_od
    call log_msg(msgstr)
#ifndef SINGLEFLUID
    write(msgstr,"('sci_wd         = ',L1,', n_sci_wd           = ',I0)") sci_wd, n_sci_wd
    call log_msg(msgstr)
    write(msgstr,"('sci_int        = ',L1,', n_sci_int          = ',I0)") sci_int, n_sci_int
    call log_msg(msgstr)
#ifndef NOSURFACTANT
    write(msgstr,"('sci_sur        = ',L1,', n_sci_sur          = ',I0)") sci_sur, n_sci_sur
    call log_msg(msgstr)
    write(msgstr,"('sci_dir        = ',L1,', n_sci_dir          = ',I0)") sci_dir, n_sci_dir
    call log_msg(msgstr)
#endif
#endif
    write(msgstr,"('sci_vel        = ',L1,', n_sci_vel          = ',I0)") sci_vel, n_sci_vel
    call log_msg(msgstr)
    write(msgstr,"('sci_flo        = ',L1,', n_sci_flo          = ',I0)") sci_flo, n_sci_flo
    call log_msg(msgstr)
    write(msgstr,"('sci_vel_correction        = ',L1)") sci_vel_correction
    call log_msg(msgstr)

    write(msgstr,"('sci_arrows     = ',L1,', n_sci_arrows       = ',I0)") sci_arrows, n_sci_arrows
    call log_msg(msgstr)
    write(msgstr,"('sci_velocities = ',L1,', n_sci_velocities   = ',I0)") sci_velocities, n_sci_velocities
    call log_msg(msgstr)
    write(msgstr,"('sci_velocities_od = ',L1,', n_sci_velocities_od   = ',I0)") sci_velocities_od, n_sci_velocities_od
    call log_msg(msgstr)
    write(msgstr,"('sci_velocities_wd = ',L1,', n_sci_velocities_wd   = ',I0)") sci_velocities_wd, n_sci_velocities_wd
    call log_msg(msgstr)
    write(msgstr,"('sci_flux_od = ',L1,', n_sci_flux_od   = ',I0)") sci_flux_od, n_sci_flux_od
    call log_msg(msgstr)
    write(msgstr,"('sci_flux_wd = ',L1,', n_sci_flux_wd   = ',I0)") sci_flux_wd, n_sci_flux_wd
    call log_msg(msgstr)

    write(msgstr,"('sci_rock       = ',L1,', n_sci_rock         = ',I0)") sci_rock, n_sci_rock
    call log_msg(msgstr)
    write(msgstr,"('sci_rock_colour= ',L1,', n_sci_rock_colour  = ',I0)") sci_rock_colour, n_sci_rock_colour
    call log_msg(msgstr)
    write(msgstr,"('sci_rock_rho_r = ',L1,', n_sci_rock_rho_r   = ',I0)") sci_rock_rho_r, n_sci_rock_rho_r
    call log_msg(msgstr)
    write(msgstr,"('sci_rock_rho_b = ',L1,', n_sci_rock_rho_b   = ',I0)") sci_rock_rho_b, n_sci_rock_rho_b
    call log_msg(msgstr)
    write(msgstr,"('sci_pressure   = ',L1,', n_sci_pressure     = ',I0)") sci_pressure, n_sci_pressure
    call log_msg(msgstr)
    write(msgstr,"('                    sci_pressure_init  = ',I0)") sci_pressure_init
    call log_msg(msgstr)
    write(msgstr,"('sci_profile    = ',L1,', n_sci_profile      = ',I0)") sci_profile, n_sci_profile
    call log_msg(msgstr)
    write(msgstr,"('                    n_sci_profile_dump = ',I0)") n_sci_profile_dump
    call log_msg(msgstr)
    write(msgstr,"('sci_arrstats   = ',L1,', n_sci_arrstats     = ',I0)") sci_arrstats, n_sci_arrstats
    call log_msg(msgstr)
    write(msgstr,"('sci_fluxz      = ',L1,', n_sci_fluxz        = ',I0)") sci_fluxz, n_sci_fluxz
    call log_msg(msgstr)
    write(msgstr,"('sci_massfluxz  = ',L1,', n_sci_massfluxz    = ',I0)") sci_massfluxz, n_sci_massfluxz
    call log_msg(msgstr)
    write(msgstr,"('sci_stress     = ',L1,', n_sci_stress       = ',I0)") sci_stress, n_sci_stress
    call log_msg(msgstr)
#ifndef SINGLEFLUID
    write(msgstr,"('sci_colour_clusters = ',L1)") sci_colour_clusters
    call log_msg(msgstr)
    if ( sci_colour_clusters ) then
      write(msgstr,"('n_sci_colour_clusters       = ',I0)") n_sci_colour_clusters
      call log_msg(msgstr)
      write(msgstr,"('n_sci_colour_clusters_index = ',I0)") n_sci_colour_clusters_index
      call log_msg(msgstr)
      write(msgstr,"('n_sci_start_colour_clusters = ',I0)") n_sci_start_colour_clusters
      call log_msg(msgstr)
    end if
    if (dbg_report_hk_timing) then
       write(msgstr,"('dbg_report_hk_timing      = ',L1)") dbg_report_hk_timing
       call log_msg(msgstr)
    end if
#endif
    call log_ws()

    if (any(fluxz_xhi /= -1)) then
      call log_msg("Regions for fluxz / massfluxz:")
      do i = 1, fluxz_regions
        if (fluxz_xhi(i) .ne. -1) then
          write(msgstr,"('  <',A,'>: x = [',I0,':',I0,'], y = [',I0,':',I0,'], z = [',I0,':',I0,']')") &
            trim(fluxz_name(i)), fluxz_xlo(i), fluxz_xhi(i), fluxz_ylo(i), fluxz_yhi(i), fluxz_zlo(i), fluxz_zhi(i)
          call log_msg(msgstr)
        end if
      end do
    else
      call log_msg("No active fluxz regions defined")
    end if

    if (any(n_sci_arrstats_dump/=0)) then
      call log_msg("Intervals for arrstats:")
      do i = 1,arrstats_intervals
        if (n_sci_arrstats_dump(i)/=0) then
          write(msgstr,"('  <',A,'>: n_sci_arrstats_dump = ',I0)") &
               &trim(arrstats_name(i)),n_sci_arrstats_dump(i)
          call log_msg(msgstr)
        end if
      end do
    else
      call log_msg("No active arrstats intervals defined")
    end if

    call log_ws()
    write(msgstr,"('post        = ',L1)") post
    call log_msg(msgstr)

    write(msgstr,"('folder      = <',A,'>')") trim(folder)
    call log_msg(msgstr)
    write(msgstr,"('cpfolder    = <',A,'>')") trim(cpfolder)
    call log_msg(msgstr)
    if (trim(srccpfolder) .ne. '') then
      write(msgstr,"('srccpfolder = <',A,'>')") trim(srccpfolder)
      call log_msg(msgstr)
    end if
    write(msgstr,"('gr_out_file = <',A,'>')") trim(gr_out_file)
    call log_msg(msgstr)

    write(msgstr,"('restore     = ',L1)") restore
    call log_msg(msgstr)
    write(msgstr,"('init_cond   = ',I0)") init_cond
    call log_msg(msgstr)
#ifdef INTERPOLATEDBB
  if (.not.(init_cond.eq.-4 .or. init_cond.eq.-5)) then
    call error("INTERPOLATED bounceback scheme works with init_cond = -4,-5! "&
              &// " Try copying the rho_0 line inside the cases -4 or -5 to " &
              &// "reqd init_cond. rho_0 is needed for computing modified f_eq")
  endif
#endif
    write(msgstr,"('vel_poiseuille_max = ',F16.10)") vel_poiseuille_max
    call log_msg(msgstr)
    write(msgstr,"('inv_fluid   = ',I0,', inv_type = ',I0)") inv_fluid, inv_type
    call log_msg(msgstr)

    call log_ws()

    write(msgstr,"('fr          = ',F16.10)") fr
    call log_msg(msgstr)
#ifndef SINGLEFLUID
    write(msgstr,"('fb          = ',F16.10)") fb
    call log_msg(msgstr)
#ifndef NOSURFACTANT
    write(msgstr,"('fg          = ',F16.10)") fg
    call log_msg(msgstr)
    write(msgstr,"('fd          = ',I0)") fd
    call log_msg(msgstr)
#endif
#endif
    write(msgstr,"('fr1         = ',F16.10)") fr1
    call log_msg(msgstr)
    write(msgstr,"('fr2         = ',F16.10)") fr2
    call log_msg(msgstr)
    write(msgstr,"('pr          = ',F16.10)") pr
    call log_msg(msgstr)
    if ( inv_fluid .eq. 17 ) then
      write(msgstr,"('pb          = ',F16.10)") pb
      call log_msg(msgstr)
      write(msgstr,"('pg          = ',F16.10)") pg
      call log_msg(msgstr)
    else
#ifndef SINGLEFLUID
      write(msgstr,"('pb          = ',F16.10)") pb
      call log_msg(msgstr)
#ifndef NOSURFACTANT
      write(msgstr,"('pg          = ',F16.10)") pg
      call log_msg(msgstr)
      write(msgstr,"('pd          = ',I0)") pd
      call log_msg(msgstr)
#endif
#endif
    end if
    write(msgstr,"('qr          = ',F16.10)") qr
    call log_msg(msgstr)
#ifndef SINGLEFLUID
    write(msgstr,"('qb          = ',F16.10)") qb
    call log_msg(msgstr)
#ifndef NOSURFACTANT
    write(msgstr,"('qg          = ',F16.10)") qg
    call log_msg(msgstr)
    write(msgstr,"('qd          = ',I0)") qd
    call log_msg(msgstr)
#endif
#endif

    write(msgstr,"('m_evp       = ',F16.10)") m_evp
    call log_msg(msgstr)
    write(msgstr,"('m_evp_gr       = ',F16.10)") m_evp_gr
    call log_msg(msgstr)
    write(msgstr,"('m_evp_gb       = ',F16.10)") m_evp_gb
    call log_msg(msgstr)
    write(msgstr,"('m_evp_freq_f       = ',F16.10)") m_evp_freq_f
    call log_msg(msgstr)
    write(msgstr,"('m_evp_freq_a       = ',F16.10)") m_evp_freq_a
    call log_msg(msgstr)
    write(msgstr,"('in_evp      = (',L1,',',L1,',',L1,')')") in_evp(1), in_evp(2), in_evp(3)
    call log_msg(msgstr)
    write(msgstr,"('out_evp     = (',L1,',',L1,',',L1,')')") out_evp(1), out_evp(2), out_evp(3)
    call log_msg(msgstr)
    write(msgstr,"('m_evp_set_density       = ',L1)") m_evp_set_density
    call log_msg(msgstr)

    write(msgstr,"('rock_colour_init   = ',L1)") rock_colour_init
    call log_msg(msgstr)
    if ( rock_colour_init ) then
      write(msgstr,"('rock_colour = ',F16.10)") rock_colour
      call log_msg(msgstr)
    end if
    !write(msgstr,"('rock_colour_double   = ',L1)") rock_colour_double
    !call log_msg(msgstr)

    write(msgstr,"('beta        = ',F16.10)") beta
    call log_msg(msgstr)

    if ( init_cond == 11 ) then
      call log_ws()
      write(msgstr,"('drop_xshift = ',I0,', drop_xcut = ',I0)") drop_xshift, drop_xcut
      call log_msg(msgstr)
      write(msgstr,"('drop_yshift = ',I0,', drop_ycut = ',I0)") drop_yshift, drop_ycut
      call log_msg(msgstr)
      write(msgstr,"('drop_zshift = ',I0,', drop_zcut = ',I0)") drop_zshift, drop_zcut
      call log_msg(msgstr)
    end if

    if ( sci_stress ) then
      call log_ws()
      write(msgstr,"('stress_minx      = ',I0,', stress_maxx = ',I0)") stress_minx,stress_maxx
      call log_msg(msgstr)
      write(msgstr,"('stress_dd_cutoff = ',F16.10)") stress_dd_cutoff
      call log_msg(msgstr)
    end if

  end if

  ! These variables match the order in lbe_params.F90, for clarity
 
  call MPI_Bcast(sci_int,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(sci_sur,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(sci_od,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(sci_wd,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(sci_dir,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(sci_vel,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(sci_vel_correction,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(sci_flo,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(sci_arrows,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(sci_velocities,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(sci_velocities_od,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(sci_velocities_wd,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(sci_flux_od,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(sci_flux_wd,1,MPI_LOGICAL,0,comm_cart,ierror)

  call MPI_Bcast(sci_rock,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(sci_rock_colour,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(sci_rock_rho_r,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(sci_rock_rho_b,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(sci_pressure,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(sci_fluxz,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(sci_massfluxz,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(sci_profile,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(sci_arrstats,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(sci_stress,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(sci_colour_clusters,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(post,1,MPI_LOGICAL,0,comm_cart,ierror)

  call MPI_Bcast(gr_out_file,256,MPI_CHARACTER,0,comm_cart,ierror)
  call MPI_Bcast(folder,256,MPI_CHARACTER,0,comm_cart,ierror)
  call MPI_Bcast(cpfolder,256,MPI_CHARACTER,0,comm_cart,ierror)
  call MPI_Bcast(srccpfolder,256,MPI_CHARACTER,0,comm_cart,ierror)

  call MPI_Bcast(restore,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(init_cond,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(vel_poiseuille_max,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(sci_pressure_init,1,MPI_INTEGER,0,comm_cart,ierror)

  call MPI_Bcast(n_iteration,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(n_sci_start,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(n_sci_int,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(n_sci_sur,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(n_sci_od,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(n_sci_wd,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(n_sci_dir,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(n_sci_vel,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(n_sci_flo,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(n_sci_arrows,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(n_sci_velocities,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(n_sci_velocities_od,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(n_sci_velocities_wd,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(n_sci_flux_od,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(n_sci_flux_wd,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(n_sci_rock,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(n_sci_rock_colour,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(n_sci_rock_rho_r,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(n_sci_rock_rho_b,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(n_sci_pressure,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(n_sci_fluxz,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(n_sci_massfluxz,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(n_sci_profile,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(n_sci_profile_dump,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(n_sci_arrstats,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(n_sci_stress,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(n_sci_colour_clusters,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(n_sci_colour_clusters_index,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(n_sci_start_colour_clusters,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(dbg_report_hk_timing,1,MPI_LOGICAL,0,comm_cart,ierror)

  do i = 1, fluxz_regions
    call MPI_Bcast(fluxz_name(i),32,MPI_CHARACTER,0,comm_cart,ierror)
  end do

  call MPI_Bcast(fluxz_xlo,fluxz_regions,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(fluxz_ylo,fluxz_regions,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(fluxz_zlo,fluxz_regions,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(fluxz_xhi,fluxz_regions,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(fluxz_yhi,fluxz_regions,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(fluxz_zhi,fluxz_regions,MPI_INTEGER,0,comm_cart,ierror)

  call MPI_Bcast(n_sci_arrstats_dump,arrstats_intervals,MPI_INTEGER,0,comm_cart,ierror)
  do i = 1,arrstats_intervals
    call MPI_Bcast(arrstats_name(i),32,MPI_CHARACTER,0,comm_cart,ierror)
  end do

  call MPI_Bcast(nr,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(fr,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(fb,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(fg,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(fd,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(fr1,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(fr2,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(pr,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(pb,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(pg,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(pd,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(qr,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(qb,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(qg,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(qd,1,MPI_INTEGER,0,comm_cart,ierror)

  call MPI_Bcast(m_evp,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(m_evp_gr,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(m_evp_gb,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(m_evp_freq_f,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(m_evp_freq_a,1,MPI_REAL8,0,comm_cart,ierror)
  call MPI_Bcast(in_evp,3,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(out_evp,3,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(m_evp_set_density,1,MPI_LOGICAL,0,comm_cart,ierror)

  call MPI_Bcast(rock_colour_init,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(rock_colour,1,MPI_REAL8,0,comm_cart,ierror)
  !call MPI_Bcast(rock_colour_double,1,MPI_LOGICAL,0,comm_cart,ierror)
  call MPI_Bcast(inv_fluid,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(inv_type,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(beta,1,MPI_REAL8,0,comm_cart,ierror)

  call MPI_Bcast(stress_minx,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(stress_maxx,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(stress_dd_cutoff,1,MPI_REAL8,0,comm_cart,ierror)

  call MPI_Bcast(drop_xshift,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(drop_yshift,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(drop_zshift,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(drop_xcut,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(drop_ycut,1,MPI_INTEGER,0,comm_cart,ierror)
  call MPI_Bcast(drop_zcut,1,MPI_INTEGER,0,comm_cart,ierror)

  ! Filter incompatibilities related to VARTAU.
#ifdef VARTAU
  if (inv_fluid .eq. 11) then
    call log_msg("FATAL ERROR: Compiler flag VARTAU and invading fluid 11 (Zou-He) are not compatible. Aborting...")
    call Abend
  end if

  if(sci_profile .eqv. .true.) then
    call log_msg("FATAL ERROR: Compiler flag VARTAU and sci_profiles = .true. are not compatible. Aborting...")
    call Abend
  end if
#endif

#ifdef MD
  ! Profile calculation requires the velocities of neighboring
  ! particles to be available. Of course, one could optimize
  ! this by introducing separate communication calls in
  !  dump_profile() .
  if (sci_profile) communicate_velocities = .true.
#endif

  call log_msg_ws("Read and broadcast variable input.")
end subroutine lbe_get_variable_input

!> update site occupation counters and print them nicely formatted to
!> standard output
!>
!> \param[in] N local lattice chunk (expects halo of extent 1)
subroutine lbe_print_site_counts(N)
    type(lbe_site),dimension(0:,0:,0:),intent(in) :: N

    call lbe_count_sites(N)

    call log_ws()
    call log_msg('Current site occupation      %      sites')
    write (unit=msgstr,fmt='("-fluid                  ",F6.2,X,I10)') &
         &100.0_rk*real(n_sites_fluid,kind=rk)/real(n_sites,kind=rk)&
         &,n_sites_fluid
    call log_msg(msgstr)
    write (unit=msgstr,fmt='("-particle               ",F6.2,X,I10)') &
         &100.0_rk*real(n_sites_particle,kind=rk)/real(n_sites,kind=rk)&
         &,n_sites_particle
    call log_msg(msgstr)
    write (unit=msgstr,fmt='("-rock                   ",F6.2,X,I10)') &
         &100.0_rk*real(n_sites_rock,kind=rk)/real(n_sites,kind=rk)&
         &,n_sites_rock
    call log_msg(msgstr)
    write (unit=msgstr,fmt='("-rock surface           ",F6.2,X,I10)') &
         &100.0_rk*real(n_sites_surface,kind=rk)/real(n_sites,kind=rk)&
         &,n_sites_surface
    call log_msg(msgstr)
    write (unit=msgstr,fmt='("-total                  ",F6.2,X,I10)') &
         &100.0_rk,n_sites
    call log_msg(msgstr)
    call log_ws()
end subroutine lbe_print_site_counts

!>This function collects the data from all processors and dumps the output
!>into a single file.
!>
!>Added 14.06.02 by Jens
subroutine postprocess(N)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N
  type(lbe_site), dimension(:,:,:),allocatable :: Nm
  integer :: ierror,x,y,z
  integer :: upstream(3), wakestream(3), shearstream(3)

!  if (check_dump_now(.true., 1) ) then
!    upstream    = (/144, 144,  70/)
!    wakestream  = (/144, 154, 204/)
!    shearstream = (/144, 154, 294/)
!!    upstream    = (/tnx/2,tny/2   ,(tnz*2/5- 60)/)
!!    wakestream  = (/tnx/2,tny/2+10,(tnz*2/5+ 30)/)
!!    shearstream = (/tnx/2,tny/2+20,(tnz*2/5+150)/)
!!    write(*, "('upstream    = ',3I3)") upstream
!!    write(*, "('wakestream  = ',3I3)") wakestream
!!    write(*, "('shearstream = ',3I3)") shearstream
!    call log_msg("  Dumping VEL_PROBE...")
!    call dump_vel_probe(N,upstream)
!    call dump_vel_probe(N,shearstream)
!    call dump_vel_probe(N,wakestream)
!  end if

   if (    ( check_dump_now(sci_int, n_sci_int) ) &
       .or.( check_dump_now(sci_sur, n_sci_sur) ) &
       .or.( check_dump_now(sci_od , n_sci_od ) ) &
       .or.( check_dump_now(sci_wd , n_sci_wd ) ) &
       .or.( check_dump_now(sci_dir, n_sci_dir) ) &
       .or.( check_dump_now(sci_vel, n_sci_vel) ) &
       .or.( check_dump_now(sci_flo, n_sci_flo) ) &
       .or.( check_dump_now(sci_arrows, n_sci_arrows) ) &
       .or.( check_dump_now(sci_velocities, n_sci_velocities) ) &
       .or.( check_dump_now(sci_velocities_od, n_sci_velocities_od) ) &
       .or.( check_dump_now(sci_velocities_wd, n_sci_velocities_wd) ) &
       .or.( check_dump_now(sci_flux_od, n_sci_flux_od) ) &
       .or.( check_dump_now(sci_flux_wd, n_sci_flux_wd) ) &
       .or.( check_dump_now(sci_rock, n_sci_rock) ) &
       .or.( check_dump_now(sci_pressure, n_sci_pressure) ) &
       .or.( check_dump_now(sci_fluxz, n_sci_fluxz) ) &
       .or.( check_dump_now(sci_massfluxz, n_sci_massfluxz) ) &
       .or.( check_dump_now(sci_colour_clusters, n_sci_colour_clusters) ) &
       .or.( check_dump_now(sci_colour_clusters, n_sci_colour_clusters_index) ) &
       ) then

    call log_msg_hdr("Dumping data")

    if(index(dump_format,'hdf').gt.0)then
      call log_msg("Dumping HDF5 data:")
#ifdef USEHDF
      !N is local array of size N(0:nx+1,0:ny+1,0:nz+1) and each bit will be
      !kept on its PE
      call dump_data(N)
#else
      call log_msg("FATAL ERROR: HDF5 support is switched off but HDF5 output is requested. Aborting...")
      call Abend
#endif
      call log_msg('Finished dumping HDF5 data.')
    else ! not phdf5
      call log_msg("Postprocessing, no HDF5")
      if (myrankc == 0) then
        ! Allocate the whole arena in the master CPU
        allocate(Nm(0:tnx+1,0:tny+1,0:tnz+1),stat=ierror)
        if (ierror .ne. 0) then
          call log_msg("FATAL ERROR: Unable to allocate memory for postprocessing data buffer. Aborting...",.true.)
          call Abend
        end if
        ! Collect the data
        call recv_final_lattice(N,Nm)
        ! Finally, dump data to disk.
        call log_msg("Dumping data.")
        call dump_data(Nm)

        deallocate(Nm)
      else    ! myrankc != 0
        ! Do nothing than sending my chunk to the master
        call send_final_lattice(N)
      end if  !myrankc
    end if !phdf5
    call log_msg_ws("Finished dumping data.")
  end if
end subroutine postprocess

!> write info on which compiler flags were used to stdout
!>
!> 2010-05-03: Added by Stefan, use the 'find_flags.sh' script to easily
!> check if this is still up to date!
subroutine lbe_detect_flags()
  implicit none
#ifdef BOUNCEBACK
  call register_flag("BOUNCEBACK")
#endif
#ifdef BUGGYIFORT11
  call register_flag("BUGGYIFORT11")
#endif
#ifdef BUGGYSENDINCOLLECT
  call register_flag("BUGGYSENDINCOLLECT")
#endif
#ifdef OLD_VEL
  call register_flag("OLD_VEL")
#endif
!#ifdef COMMON_VEL_FIX
!  call register_flag("COMMON_VEL_FIX")
!#endif
#ifdef DEBUG
  call register_flag("DEBUG")
#endif
#ifdef DEBUG_LE
  call register_flag("DEBUG_LE")
#endif
#ifdef DEBUG_MPI
  call register_flag("DEBUG_MPI")
#endif
#ifdef DEBUG_REPORTMDCOMM
  call register_flag("DEBUG_REPORTMDCOMM")
#endif
#ifdef DIST
  call register_flag("DIST")
#endif
#ifdef ELEC
  call register_flag("ELEC")
#endif
#ifdef ELEC_NNONREST
  call register_flag("ELEC_NNONREST")
#endif
#ifdef ELEC_NEWFORCE
  call register_flag("ELEC_NEWFORCE")
#endif
#ifdef FASTBDIST2
  call register_flag("FASTBDIST2")
#endif
#ifdef HDF5_FLIP
  call register_flag("HDF5_FLIP")
#endif
#ifdef IBM_BINARYIBM
  call register_flag("IBM_BINARYIBM")
#endif
#ifdef IBM_DEBUG
  call register_flag("IBM_DEBUG")
#endif
#ifdef IBM_EMULATETENSION
  call register_flag("IBM_EMULATETENSION")
#endif
#ifdef IBM_INDEXFIELD
  call register_flag("IBM_INDEXFIELD")
#endif
#ifdef IBM_FIXED
	call register_flag("IBM_FIXED")
#endif
#ifdef IBM_DRAG
	call register_flag("IBM_DRAG")
#endif
#ifdef IBM_PART
  call register_flag("IBM_PART")
#endif
#ifdef IBM_SWIMMER
  call register_flag("IBM_SWIMMER")
#endif
#ifdef LADD_SSD
  call register_flag("LADD_SSD")
#endif
#ifdef LADD_DLUB
  call register_flag("LADD_DLUB")
#endif
#ifdef LADD_DLUB_CENTRAL
  call register_flag("LADD_DLUB_CENTRAL")
#endif
#ifdef LADD_GALILEAN
  call register_flag("LADD_GALILEAN")
#endif
#ifdef LADD_SURR_RHOF
  call register_flag("LADD_SURR_RHOF")
#endif
#ifdef LBE_RS
  call register_flag("LBE_RS")
#endif
#ifdef LE_NOPVEL
  call register_flag("LE_NOPVEL")
#endif
#ifdef LOCALBC
  call register_flag("LOCALBC")
#endif
#ifdef MCMP
  call register_flag("MCMP")
#endif
#ifdef MD
  call register_flag("MD")
#endif
#ifdef MPI_ALLGV_FASTER_THAN_GV
  call register_flag("MPI_ALLGV_FASTER_THAN_GV")
#endif
#ifdef NEGGAL
  call register_flag("NEGGAL")
#endif
#ifdef NEWELLIPSOIDMETHOD
  call register_flag("NEWELLIPSOIDMETHOD")
#endif
#ifdef NOEDGESTEP
  call register_flag("NOEDGESTEP")
#endif
#ifdef NOFLUIDDYNAMICS
  call register_flag("NOFLUIDDYNAMICS")
#endif
#ifdef NOIEEEARITHMETIC
  call register_flag("NOIEEEARITHMETIC")
#endif
#ifdef NOISNAN
  call register_flag("NOISNAN")
#endif
#ifdef NOSURFACTANT
  call register_flag("NOSURFACTANT")
#endif
#ifdef OLDRRFORCE
  call register_flag("OLDRRFORCE")
#endif
#ifdef PARTICLESTRESS
  call register_flag("PARTICLESTRESS")
#endif
#ifdef P3M
  call register_flag("P3M")
#endif
#ifdef REDUCED_EXCHANGE
  call register_flag("REDUCED_EXCHANGE")
#endif
#ifdef RWALK
  call register_flag("RWALK")
#endif
#ifdef SINGLEFLUID
  call register_flag("SINGLEFLUID")
#endif
#ifdef TRACER
  call register_flag("TRACER")
#endif
#ifdef USEHDF
  call register_flag("USEHDF")
#endif
#ifdef USEXDRF
  call register_flag("USEXDRF")
#endif
#ifdef VARTAU
  call register_flag("VARTAU")
#endif
#ifdef AXISYM
  call register_flag("AXISYM")
#endif


end subroutine lbe_detect_flags

subroutine register_flag(flag)
  implicit none

  character(len=*),intent(in) :: flag
  integer :: ierror
  character(len=max_flag_len) ,allocatable,dimension(:) :: tmp

  if (.not.allocated(flag_list)) then
    allocate (flag_list(0),stat=ierror)
    call check_allocate(ierror,'register_flag(): flag_list')
  end if

  allocate (tmp(size(flag_list)),stat=ierror)
  call check_allocate(ierror,'register_flag(): tmp')
  tmp(:) = flag_list(:)
  deallocate (flag_list)
  allocate (flag_list(size(tmp)+1),stat=ierror)
  call check_allocate(ierror,'register_timer(): flag_list')
  flag_list(1:size(tmp)) = tmp(:)
  deallocate (tmp)

  flag_list(size(flag_list)) = flag
end subroutine register_flag

subroutine lbe_report_flags()
  implicit none
  integer :: i

  call log_msg_hdr("Reporting compiler flags")
  do i=1,size(flag_list)
    call log_msg("  "//trim(flag_list(i)))
  end do
  call log_ws()

end subroutine lbe_report_flags
!!!Attention: this function calculates momentum based on the distribution after collision.
!!!Please use sci_vel_correction = .true. if you want to calculate the
!physical momentum based on Eq.(12) in paper Shan Doolen PRE 81,379,1995 
!subroutine corrected_velocity(local_distribution, local_force, velocity)

!  type(lbe_site),intent(in) :: local_distribution
!  real(kind=rk),dimension(1:3,0:n_spec), intent(in) :: local_force
!
!  real(kind=rk),dimension(1:3,0:n_spec), intent(out) :: velocity

!  real(kind=rk) :: rho, rho_tot
!  real(kind=rk),dimension(1:3) :: momentum
!  integer ::  s 
  ! This is the most simple version of a force correction
  ! Redundant with the force correction in lsuperstruct_interface
  ! Still reimplemented here to keep things apart.
  
!  rho = 0.0_rk
!  rho_tot = 0.0_rk
!  momentum = 0.0_rk
!  velocity = 0.0_rk

 
!#ifdef COMMON_VEL_FIX
!#ifndef OLD_VEL

!if (sci_vel_correction) then
!    rho=sum(local_distribution%n_r_pre(:) * g(:))
!    momentum(1)=amass_r * sum(local_distribution%n_r_pre(:) * g(:) * c(:,1))
!    momentum(2)=amass_r * sum(local_distribution%n_r_pre(:) * g(:) * c(:,2))
!    momentum(3)=amass_r * sum(local_distribution%n_r_pre(:) * g(:) * c(:,3))
!else
!  rho=sum(local_distribution%n_r(:) * g(:))
!  momentum(1)=amass_r * sum(local_distribution%n_r(:) * g(:) * c(:,1))
!  momentum(2)=amass_r * sum(local_distribution%n_r(:) * g(:) * c(:,2))
!  momentum(3)=amass_r * sum(local_distribution%n_r(:) * g(:) * c(:,3))
!
!end if
!  rho_tot = rho_tot + rho * amass_r
!  velocity(:,0) = ( momentum(:) + 0.5_rk * (local_force(:,0) + local_force(:,1) ) )
!  !velocity(:,1) = ( momentum(:) - 0.5_rk * (local_force(:,0) + local_force(:,1) ) ) / max(esmall, rho)
!  velocity(:,1) = 0.0 !( momentum(:) - 0.5_rk * (local_force(:,0) + local_force(:,1) ) ) / max(esmall, rho)

!#ifndef SINGLEFLUID
! 
!if (sci_vel_correction) then
!		rho=sum(local_distribution%n_b_pre(:) * g(:))
!		momentum(1)=amass_b * sum(local_distribution%n_b_pre(:) * g(:) * c(:,1))
!      momentum(2)=amass_b * sum(local_distribution%n_b_pre(:) * g(:) * c(:,2))
!      momentum(3)=amass_b * sum(local_distribution%n_b_pre(:) * g(:) * c(:,3))
!else
!		rho=sum(local_distribution%n_b(:) * g(:))
!  		momentum(1)=amass_b * sum(local_distribution%n_b(:) * g(:) * c(:,1))
!  		momentum(2)=amass_b * sum(local_distribution%n_b(:) * g(:) * c(:,2))
!  		momentum(3)=amass_b * sum(local_distribution%n_b(:) * g(:) * c(:,3)) 
!end if 
!  rho_tot = rho_tot + rho * amass_b
!  velocity(:,0) = velocity(:,0) + ( momentum(:) + 0.5_rk * (local_force(:,2) ) )
!  velocity(:,2) = ( momentum(:) + 0.5_rk * (local_force(:,0) + local_force(:,2) ) ) / max(esmall, rho)
!#ifndef NOSURFACTANT 
!  
!if (sci_vel_correction) then
!        rho=sum(local_distribution%n_s_pre(:) * g(:))
!		  momentum(1)=amass_s * sum(local_distribution%n_s_pre(:) * g(:) * c(:,1))
!        momentum(2)=amass_s * sum(local_distribution%n_s_pre(:) * g(:) * c(:,2))
!        momentum(3)=amass_s * sum(local_distribution%n_s_pre(:) * g(:) * c(:,3))
!else 
!  rho=sum(local_distribution%n_s(:) * g(:))
!  momentum(1)=amass_s * sum(local_distribution%n_s(:) * g(:) * c(:,1))
!  momentum(2)=amass_s * sum(local_distribution%n_s(:) * g(:) * c(:,2))
!  momentum(3)=amass_s * sum(local_distribution%n_s(:) * g(:) * c(:,3))
!end if
!  rho_tot = rho_tot + rho * amass_s
!  velocity(:,0) = velocity(:,0) + ( momentum(:) + 0.5_rk * (local_force(:,3) ) )
!  velocity(:,3) = ( momentum(:) + 0.5_rk * (local_force(:,0) + local_force(:,3) ) ) / max(esmall, rho)
! endif ifndef NOSURFACTANT
!#endif
! endif ifndef SINGLEFLUID
!#endif

!  velocity(:,0) = velocity(:,0) / max(esmall, rho_tot)
!#endif
 ! endif COMMON_VEL_FIX
 !endif OLD_VEL

!end subroutine corrected_velocity
subroutine corrected_velocity(local_distribution, local_force, velocity)

  type(lbe_site),intent(in) :: local_distribution
  real(kind=rk),dimension(1:3,0:n_spec), intent(in) :: local_force
  real(kind=rk),dimension(1:3,0:n_spec), intent(out) :: velocity
  !real(kind=rk),dimension(1:3,0:n_spec), intent(out) :: vpre
  

  real(kind=rk) :: rho, rho_tot,rho_od, rho_wd,rho_s
  real(kind=rk),dimension(1:3) :: momentum
  real(kind=rk),dimension(1:3) :: momentum_od
  real(kind=rk),dimension(1:3) :: momentum_wd
  real(kind=rk),dimension(1:3) :: momentum_s
  real(kind=rk),dimension(1:3,0:n_spec):: flux
   integer :: x, y, z, s
  ! This is the most simple version of a force correction
  ! Redundant with the force correction in lsuperstruct_interface
  ! Still reimplemented here to keep things apart.
  
  rho = 0.0_rk
  rho_od = 0.0_rk
  rho_tot = 0.0_rk
  momentum = 0.0_rk
  velocity = 0.0_rk
  flux = 0.0_rk

!#ifdef COMMON_VEL_FIX 
#ifndef OLD_VEL 
 
if (sci_vel_correction) then
    rho=sum(local_distribution%n_r_pre(:) * g(:))
    momentum(1)=amass_r * sum(local_distribution%n_r_pre(:) * g(:) * c(:,1))
    momentum(2)=amass_r * sum(local_distribution%n_r_pre(:) * g(:) * c(:,2))
    momentum(3)=amass_r * sum(local_distribution%n_r_pre(:) * g(:) * c(:,3))
else
  rho=sum(local_distribution%n_r(:) * g(:))
  momentum(1)=amass_r * sum(local_distribution%n_r(:) * g(:) * c(:,1))
  momentum(2)=amass_r * sum(local_distribution%n_r(:) * g(:) * c(:,2))
  momentum(3)=amass_r * sum(local_distribution%n_r(:) * g(:) * c(:,3))

end if
  momentum_od = momentum
  !vpre(:,1) = momentum_od/rho
  rho_tot = rho_tot + rho * amass_r
  velocity(:,0) = ( momentum_od(:) + 0.5_rk * (local_force(:,0) + local_force(:,1) ) )
  velocity(:,1) = ( momentum(:) + 0.5_rk * (local_force(:,0) + local_force(:,1) ) ) / max(esmall, rho)
  rho_od = rho
  flux(:,1) = ( momentum_od(:) + (local_force(:,0) + local_force(:,1) ) )
#ifndef SINGLEFLUID
 
if (sci_vel_correction) then
		rho=sum(local_distribution%n_b_pre(:) * g(:))
		momentum(1)=amass_b * sum(local_distribution%n_b_pre(:) * g(:) * c(:,1))
      momentum(2)=amass_b * sum(local_distribution%n_b_pre(:) * g(:) * c(:,2))
      momentum(3)=amass_b * sum(local_distribution%n_b_pre(:) * g(:) * c(:,3))
else
		rho=sum(local_distribution%n_b(:) * g(:))
  		momentum(1)=amass_b * sum(local_distribution%n_b(:) * g(:) * c(:,1))
  		momentum(2)=amass_b * sum(local_distribution%n_b(:) * g(:) * c(:,2))
  		momentum(3)=amass_b * sum(local_distribution%n_b(:) * g(:) * c(:,3)) 
end if
  momentum_wd = momentum
  !vpre(:,2) = momentum_wd/rho 
  rho_tot = rho_tot + rho * amass_b
  velocity(:,0) = velocity(:,0) + ( momentum(:) + 0.5_rk * (local_force(:,2) ) )
  !velocity(:,2) = ( momentum(:) + 0.5_rk * (local_force(:,0) + local_force(:,2) ) ) / max(esmall, rho)
  rho_wd = rho
  !flux(:,2) = momentum(:) + 0.5_rk * (local_force(:,2) ) 
  flux(:,2) = momentum_wd(:) +  (local_force(:,0) + local_force(:,2) ) 
#ifdef NOSURFACTANT
  flux(:,1) = 0.5_rk*(flux(:,1) + rho_od/rho_tot*(momentum_wd+momentum_od)) 
  flux(:,2) = 0.5_rk*(flux(:,2) + rho_wd/rho_tot*(momentum_wd+momentum_od)) 
  velocity(:,1) = flux(:,1) / max(esmall, rho_od)
  velocity(:,2) = flux(:,2) / max(esmall, rho_wd) 
!#ifndef NOSURFACTANT 
#else  
if (sci_vel_correction) then
        rho=sum(local_distribution%n_s_pre(:) * g(:))
      	momentum(1)=amass_s * sum(local_distribution%n_s_pre(:) * g(:) * c(:,1))
        momentum(2)=amass_s * sum(local_distribution%n_s_pre(:) * g(:) * c(:,2))
        momentum(3)=amass_s * sum(local_distribution%n_s_pre(:) * g(:) * c(:,3))
else 
  rho=sum(local_distribution%n_s(:) * g(:))
  momentum(1)=amass_s * sum(local_distribution%n_s(:) * g(:) * c(:,1))
  momentum(2)=amass_s * sum(local_distribution%n_s(:) * g(:) * c(:,2))
  momentum(3)=amass_s * sum(local_distribution%n_s(:) * g(:) * c(:,3))
end if
  rho_tot = rho_tot + rho * amass_s
  rho_s  = rho
  momentum_s = momentum
  velocity(:,0) = velocity(:,0) + ( momentum(:) + 0.5_rk * (local_force(:,3) ) )
  !velocity(:,3) = ( momentum(:) + 0.5_rk * (local_force(:,0) + local_force(:,3) ) ) / max(esmall, rho)
  !flux(:,3) = momentum_wd(:) +  (local_force(:,0) + local_force(:,2) )
  !flux(:,1) = 0.5_rk*(flux(:,1) + rho_od/rho_tot*(momentum_wd+momentum_od+momentum_s))
  !flux(:,2) = 0.5_rk*(flux(:,2) + rho_wd/rho_tot*(momentum_wd+momentum_od+ momentum_s)) 
  !flux(:,3) = 0.5_rk*(flux(:,3) + rho_s/rho_tot*(momentum_wd+momentum_od+ momentum_s)) 
   velocity(:,1) = flux(:,1) / max(esmall, rho_od)
   velocity(:,2) = flux(:,2) / max(esmall, rho_wd)
   velocity(:,3) = flux(:,3) / max(esmall, rho_s)
  !flux(:,3) = flux(:,3) - rho_s*velocity(:,0) 
! endif ifndef NOSURFACTANT
#endif
! endif ifndef SINGLEFLUID
!velocity(:,0) = velocity(:,0) / max(esmall, rho_tot)
!flux(:,2) = flux(:,2) - rho_wd*velocity(:,0)
#endif

  velocity(:,0) = velocity(:,0) / max(esmall, rho_tot)
  !flux(:,1) = flux(:,1) - rho_od*velocity(:,0)
!#ifndef SINGLEFLUID
!  flux(:,2) = flux(:,2) - rho_wd*velocity(:,0)
!#ifndef NOSURFACTANT
!   flux(:,3) = flux(:,3) - rho_s*velocity(:,0)  
!#endif
!#endif

#endif 
! endif COMMON_VEL_FIX 
! endif OLD_VEL 
end subroutine corrected_velocity


subroutine corrected_velocity_new(local_distribution, local_force, velocity, flux,x,y,z)

  type(lbe_site),intent(in) :: local_distribution
  real(kind=rk),dimension(1:3,0:n_spec), intent(in) :: local_force

  real(kind=rk),dimension(1:3,0:n_spec), intent(out) :: velocity
  !real(kind=rk),dimension(1:3,0:n_spec), intent(out) :: vpre
  real(kind=rk),dimension(1:3,0:n_spec), intent(out) :: flux

  real(kind=rk) :: rho, rho_tot,rho_od, rho_wd,rho_s
  real(kind=rk),dimension(1:3) :: momentum
  real(kind=rk),dimension(1:3) :: momentum_od
  real(kind=rk),dimension(1:3) :: momentum_wd
  real(kind=rk),dimension(1:3) :: momentum_s
   integer :: x, y, z, s
  ! This is the most simple version of a force correction
  ! Redundant with the force correction in lsuperstruct_interface
  ! Still reimplemented here to keep things apart.
  
  rho = 0.0_rk
  rho_od = 0.0_rk
  rho_tot = 0.0_rk
  momentum = 0.0_rk
  velocity = 0.0_rk
  flux = 0.0_rk

!#ifdef COMMON_VEL_FIX 
#ifndef OLD_VEL 
 
if (sci_vel_correction) then
    rho=sum(local_distribution%n_r_pre(:) * g(:))
    momentum(1)=amass_r * sum(local_distribution%n_r_pre(:) * g(:) * c(:,1))
    momentum(2)=amass_r * sum(local_distribution%n_r_pre(:) * g(:) * c(:,2))
    momentum(3)=amass_r * sum(local_distribution%n_r_pre(:) * g(:) * c(:,3))
else
  rho=sum(local_distribution%n_r(:) * g(:))
  momentum(1)=amass_r * sum(local_distribution%n_r(:) * g(:) * c(:,1))
  momentum(2)=amass_r * sum(local_distribution%n_r(:) * g(:) * c(:,2))
  momentum(3)=amass_r * sum(local_distribution%n_r(:) * g(:) * c(:,3))

end if
  momentum_od = momentum
  !vpre(:,1) = momentum_od/rho
  rho_tot = rho_tot + rho * amass_r
  velocity(:,0) = ( momentum_od(:) + 0.5_rk * (local_force(:,0) + local_force(:,1) ) )
  velocity(:,1) = ( momentum(:) + 0.5_rk * (local_force(:,0) + local_force(:,1) ) ) / max(esmall, rho)
  rho_od = rho
  flux(:,1) = ( momentum_od(:) + (local_force(:,0) + local_force(:,1) ) )
!if (x==2.AND.y ==10 .AND. z ==10) then
! print*, "density \n", rho
! print*, "mom \n",   velocity(1,0), velocity(2,0),velocity(3,0)
! print*, "force: \n",  local_force(1,1), local_force(2,1), local_force(3,1)
 !print*, "momentum 1 2 x , y, z: \n",rho,  momentum(1), momentum(2),momentum(3), p_r(1), p_r(2),p_r(3)
!do s = 1, nvecs 
!print*, s, local_distribution%n_r(s)*g(s), local_distribution%n_r_pre(s)*g(s), cx(s), cy(s), cz(s)
!end do
! READ*
!  endif
#ifndef SINGLEFLUID
 
if (sci_vel_correction) then
		rho=sum(local_distribution%n_b_pre(:) * g(:))
		momentum(1)=amass_b * sum(local_distribution%n_b_pre(:) * g(:) * c(:,1))
      momentum(2)=amass_b * sum(local_distribution%n_b_pre(:) * g(:) * c(:,2))
      momentum(3)=amass_b * sum(local_distribution%n_b_pre(:) * g(:) * c(:,3))
else
		rho=sum(local_distribution%n_b(:) * g(:))
  		momentum(1)=amass_b * sum(local_distribution%n_b(:) * g(:) * c(:,1))
  		momentum(2)=amass_b * sum(local_distribution%n_b(:) * g(:) * c(:,2))
  		momentum(3)=amass_b * sum(local_distribution%n_b(:) * g(:) * c(:,3)) 
end if
  momentum_wd = momentum
  !vpre(:,2) = momentum_wd/rho 
  rho_tot = rho_tot + rho * amass_b
  velocity(:,0) = velocity(:,0) + ( momentum(:) + 0.5_rk * (local_force(:,2) ) )
  !velocity(:,2) = ( momentum(:) + 0.5_rk * (local_force(:,0) + local_force(:,2) ) ) / max(esmall, rho)
  rho_wd = rho
  !flux(:,2) = momentum(:) + 0.5_rk * (local_force(:,2) ) 
  flux(:,2) = momentum_wd(:) +  (local_force(:,0) + local_force(:,2) ) 
#ifdef NOSURFACTANT
  flux(:,1) = 0.5_rk*(flux(:,1) + rho_od/rho_tot*(momentum_wd+momentum_od)) 
  flux(:,2) = 0.5_rk*(flux(:,2) + rho_wd/rho_tot*(momentum_wd+momentum_od)) 
  velocity(:,1) = flux(:,1) / max(esmall, rho_od)
  velocity(:,2) = flux(:,2) / max(esmall, rho_wd) 
!#ifndef NOSURFACTANT 
#else  
if (sci_vel_correction) then
        rho=sum(local_distribution%n_s_pre(:) * g(:))
      	momentum(1)=amass_s * sum(local_distribution%n_s_pre(:) * g(:) * c(:,1))
        momentum(2)=amass_s * sum(local_distribution%n_s_pre(:) * g(:) * c(:,2))
        momentum(3)=amass_s * sum(local_distribution%n_s_pre(:) * g(:) * c(:,3))
else 
  rho=sum(local_distribution%n_s(:) * g(:))
  momentum(1)=amass_s * sum(local_distribution%n_s(:) * g(:) * c(:,1))
  momentum(2)=amass_s * sum(local_distribution%n_s(:) * g(:) * c(:,2))
  momentum(3)=amass_s * sum(local_distribution%n_s(:) * g(:) * c(:,3))
end if
  rho_tot = rho_tot + rho * amass_s
  rho_s  = rho
  momentum_s = momentum
  velocity(:,0) = velocity(:,0) + ( momentum(:) + 0.5_rk * (local_force(:,3) ) )
  !velocity(:,3) = ( momentum(:) + 0.5_rk * (local_force(:,0) + local_force(:,3) ) ) / max(esmall, rho)
  flux(:,3) = momentum_wd(:) +  (local_force(:,0) + local_force(:,2) )
  flux(:,1) = 0.5_rk*(flux(:,1) + rho_od/rho_tot*(momentum_wd+momentum_od+momentum_s))
  flux(:,2) = 0.5_rk*(flux(:,2) + rho_wd/rho_tot*(momentum_wd+momentum_od+ momentum_s)) 
  flux(:,3) = 0.5_rk*(flux(:,3) + rho_s/rho_tot*(momentum_wd+momentum_od+ momentum_s)) 
   velocity(:,1) = flux(:,1) / max(esmall, rho_od)
   velocity(:,2) = flux(:,2) / max(esmall, rho_wd)
   velocity(:,3) = flux(:,3) / max(esmall, rho_s)
  !flux(:,3) = flux(:,3) - rho_s*velocity(:,0) 
! endif ifndef NOSURFACTANT
#endif
! endif ifndef SINGLEFLUID
!velocity(:,0) = velocity(:,0) / max(esmall, rho_tot)
!flux(:,2) = flux(:,2) - rho_wd*velocity(:,0)
#endif

  velocity(:,0) = velocity(:,0) / max(esmall, rho_tot)
  flux(:,1) = flux(:,1) - rho_od*velocity(:,0)
#ifndef SINGLEFLUID
  flux(:,2) = flux(:,2) - rho_wd*velocity(:,0)
#ifndef NOSURFACTANT
   flux(:,3) = flux(:,3) - rho_s*velocity(:,0)  
#endif
#endif

#endif 
! endif COMMON_VEL_FIX 
! endif OLD_VEL 
end subroutine corrected_velocity_new
!!$subroutine mass_aware_corrected_velocity_stub(local_distribution, local_force, corrected_velocity)

!!$ 
!!$  type(lbe_site),intent(in) :: local_distribution
!!$  real(kind=rk),dimension(1:3,0:n_spec), intent(in) :: local_force
!!$  real(kind=rk),dimension(1:3), intent(out) :: corrected_velocity
!!$
!!$  real(rind=rk),dimension(1:3) :: local_density
!!$  real(rind=rk),dimension(1:3,0:n_spec) :: local_velocity
!!$
!!$
!!$  ! This is just a collection of formulas, Need to test effects of force-scaling
!!$
!!$  local_density(1)=sum(local_distribution%n_r(:)*g(:))
!!$  local_density(2)=sum(local_distribution%n_b(:)*g(:))
!!$  local_density(3)=sum(local_distribution%n_s(:)*g(:))
!!$
!!$  local_mass_density(1)=amass_r*local_density(1)
!!$  local_mass_density(2)=amass_b*local_density(2)
!!$  local_mass_density(3)=amass_s*local_density(3)
!!$
!!$  local_velocity(1,1)=sum(local_distribution%n_r(:)*g(:)*c(:,1))
!!$  local_velocity(2,1)=sum(local_distribution%n_r(:)*g(:)*c(:,2))
!!$  local_velocity(3,1)=sum(local_distribution%n_r(:)*g(:)*c(:,3))
!!$
!!$  local_velocity(1,2)=sum(local_distribution%n_b(:)*g(:)*c(:,1))
!!$  local_velocity(2,2)=sum(local_distribution%n_b(:)*g(:)*c(:,2))
!!$  local_velocity(3,2)=sum(local_distribution%n_b(:)*g(:)*c(:,3))
!!$
!!$  local_velocity(1,3)=sum(local_distribution%n_s(:)*g(:)*c(:,1))
!!$  local_velocity(2,3)=sum(local_distribution%n_s(:)*g(:)*c(:,2))
!!$  local_velocity(3,3)=sum(local_distribution%n_s(:)*g(:)*c(:,3))
!!$
!!$  local_momentum(:,1) = amass_r*local_velocity(:,1)
!!$  local_momentum(:,2) = amass_b*local_velocity(:,2)
!!$  local_momentum(:,3) = amass_s*local_velocity(:,3)
!!$
!!$  local_velocity(1,0)=sum(local_momentum(1,:))/sum(local_mass_density(:))
!!$  local_velocity(2,0)=sum(local_momentum(2,:))/sum(local_mass_density(:))
!!$  local_velocity(3,0)=sum(local_momentum(3,:))/sum(local_mass_density(:))
!!$
!!$  corrected_velocity(1,1) = (local_momentum(1,1) - delta_t*(local_force(1,0)+local_force(1,1))/2.0d0)/local_mass_density(1)
!!$  corrected_velocity(2,1) = (local_momentum(2,1) - delta_t*(local_force(2,0)+local_force(2,1))/2.0d0)/local_mass_density(1)
!!$  corrected_velocity(3,1) = (local_momentum(3,1) - delta_t*(local_force(3,0)+local_force(3,1))/2.0d0)/local_mass_density(1) 
!!$
!!$  corrected_velocity(1,2) = (local_momentum(1,2) - delta_t*(local_force(1,0)+local_force(1,2))/2.0d0)/local_mass_density(2)
!!$  corrected_velocity(2,2) = (local_momentum(2,2) - delta_t*(local_force(2,0)+local_force(2,2))/2.0d0)/local_mass_density(2)
!!$  corrected_velocity(3,2) = (local_momentum(3,2) - delta_t*(local_force(3,0)+local_force(3,2))/2.0d0)/local_mass_density(2)
!!$
!!$  corrected_velocity(1,3) = (local_momentum(1,3) - delta_t*(local_force(1,0)+local_force(1,3))/2.0d0)/local_mass_density(3)
!!$  corrected_velocity(2,3) = (local_momentum(2,3) - delta_t*(local_force(2,0)+local_force(2,3))/2.0d0)/local_mass_density(3)
!!$  corrected_velocity(3,3) = (local_momentum(3,3) - delta_t*(local_force(3,0)+local_force(3,3))/2.0d0)/local_mass_density(3)
!!$
!!$! Version A not so sure about this, since it is summing up forces unscaled
!!$  corrected_velocity(1,0) = (local_velocity(1,0) - (0.5_rk * sum(local_force(1,:)))
!!$  corrected_velocity(2,0) = (local_velocity(2,0) - (0.5_rk * sum(local_force(2,:)))
!!$  corrected_velocity(3,0) = (local_velocity(3,0) - (0.5_rk * sum(local_force(3,:)))
!!$
!!$! Version B
!!$  corrected_velocity(1,0) = (&
!!$       (local_momentum(1,1) - delta_t*(local_force(1,0)+local_force(1,1))/2.0d0) +&
!!$       (local_momentum(1,2) - delta_t*(local_force(1,0)+local_force(1,2))/2.0d0) +&
!!$       (local_momentum(1,3) - delta_t*(local_force(1,0)+local_force(1,3))/2.0d0) &
!!$       ) / (&
!!$       local_mass_density(1) + local_mass_density(2) + local_mass_density(3)&
!!$       )
!!$
!!$  corrected_velocity(2,0) = (&
!!$       (local_momentum(2,1) - delta_t*(local_force(2,0)+local_force(2,1))/2.0d0) +&
!!$       (local_momentum(2,2) - delta_t*(local_force(2,0)+local_force(2,2))/2.0d0) +&
!!$       (local_momentum(2,3) - delta_t*(local_force(2,0)+local_force(2,3))/2.0d0) &
!!$       ) / (&
!!$       local_mass_density(1) + local_mass_density(2) + local_mass_density(3)&
!!$       )
!!$
!!$  corrected_velocity(3,0) = (&
!!$       (local_momentum(3,1) - delta_t*(local_force(3,0)+local_force(3,1))/2.0d0) +&
!!$       (local_momentum(3,2) - delta_t*(local_force(3,0)+local_force(3,2))/2.0d0) +&
!!$       (local_momentum(3,3) - delta_t*(local_force(3,0)+local_force(3,3))/2.0d0) &
!!$       ) / (&
!!$       local_mass_density(1) + local_mass_density(2) + local_mass_density(3)&
!!$       )

!!$  corrected_velocity
!!$
!!$  v_r(1,x,y,z)=sum(N(x,y,z)%n_r(:)*g*cx)*amass_r/max(10.e-9,real( sum(N(x,y,z)%n_r(:)*g)*amass_r))
!!$  v_r(2,x,y,z)=sum(N(x,y,z)%n_r(:)*g*cy)*amass_r/max(10.e-9,real( sum(N(x,y,z)%n_r(:)*g)*amass_r))
!!$  v_r(3,x,y,z)=sum(N(x,y,z)%n_r(:)*g*cz)*amass_r/max(10.e-9,real( sum(N(x,y,z)%n_r(:)*g)*amass_r))
!!$
!!$  v_r(1,x,y,z)=v_r(1,x,y,z)-(lbe_force(1,0.x,y,z)+lbe_force(1,1.x,y,z))/2.0d0
!!$  v_r(2,x,y,z)=v_r(2,x,y,z)-(lbe_force(2,0.x,y,z)+lbe_force(2,1.x,y,z))/2.0d0
!!$  v_r(3,x,y,z)=v_r(3,x,y,z)-(lbe_force(3,0.x,y,z)+lbe_force(3,1.x,y,z))/2.0d0
!!$
!!$
!!$  v(1,x,y,z)=sum(N(x,y,z)%n_r(:)*g*cx)*amass_r&
!!$       &-(f_r(x,y,z,1)/2.0d0*sum(N(x,y,z)%n_r(:)*g)*amass_r)&
!!$       &+sum(N(x,y,z)%n_b(:)*g*cx)*amass_b&
!!$       &-(f_b(x,y,z,1)/2.0d0*sum(N(x,y,z)%n_b(:)*g)*amass_b)&
!!$       &+sum(N(x,y,z)%n_s(:)*g*cx)*amass_s&
!!$       &-(f_s(x,y,z,1)/2.0d0*sum(N(x,y,z)%n_s(:)*g)*amass_s)
!!$
!!$
!!$        v(2,x,y,z)=sum(N(x,y,z)%n_r(:)*g*cy)*amass_r&
!!$          &-(f_r(x,y,z,2)/2.0d0*sum(N(x,y,z)%n_r(:)*g)*amass_r)&
!!$          &+sum(N(x,y,z)%n_b(:)*g*cy)*amass_b&
!!$          &-(f_b(x,y,z,2)/2.0d0*sum(N(x,y,z)%n_b(:)*g)*amass_b)&
!!$          &+sum(N(x,y,z)%n_s(:)*g*cy)*amass_s&
!!$          &-(f_s(x,y,z,2)/2.0d0*sum(N(x,y,z)%n_s(:)*g)*amass_s)
!!$        v(3,x,y,z)=sum(N(x,y,z)%n_r(:)*g*cz)*amass_r&
!!$          &-(f_r(x,y,z,3)/2.0d0*sum(N(x,y,z)%n_r(:)*g)*amass_r)&
!!$          &+sum(N(x,y,z)%n_b(:)*g*cz)*amass_b&
!!$          &-(f_b(x,y,z,3)/2.0d0*sum(N(x,y,z)%n_b(:)*g)*amass_b)&
!!$          &+sum(N(x,y,z)%n_s(:)*g*cz)*amass_s&
!!$          &-(f_s(x,y,z,3)/2.0d0*sum(N(x,y,z)%n_s(:)*g)*amass_s)
!!$
!!$        v(1,x,y,z)=v(1,x,y,z)/max(10.e-9, real(sum(N(x,y,z)%n_r(:)*g)*amass_r &
!!$                                            & +sum(N(x,y,z)%n_b(:)*g)*amass_b &
!!$                                            & +sum(N(x,y,z)%n_s(:)*g)*amass_s) )
!!$        v(2,x,y,z)=v(2,x,y,z)/max(10.e-9, real(sum(N(x,y,z)%n_r(:)*g)*amass_r &
!!$                                            & +sum(N(x,y,z)%n_b(:)*g)*amass_b &
!!$                                            & +sum(N(x,y,z)%n_s(:)*g)*amass_s) )
!!$        v(3,x,y,z)=v(3,x,y,z)/max(10.e-9, real(sum(N(x,y,z)%n_r(:)*g)*amass_r &
!!$                                            & +sum(N(x,y,z)%n_b(:)*g)*amass_b &
!!$                                            & +sum(N(x,y,z)%n_s(:)*g)*amass_s) )


!!$end subroutine mass_aware_corrected_velocity_stub

!> performs checks on the physical sanity of the LB system and dumps a
!> short summary to standard output
!>
!> The checks currently encompass looking for NaNs, looking for
!> negative distribution functions or too large velocities (and the
!> frequency of both), maximum velocity (and its location in the
!> lattice), averages of velocities and densities, and total momentum
!> and densities.
!>
!> \param[in] whole_N local chunk of the lattice with full halo of
!> depth \c halo_extent
!>
!> \param[in,out] insanity set to \c .true. in case of issues, never
!> set to \c .false. here
!>
!> \note All checks should be put into the same routine in order to
!> have a single loop over the whole lattice only.
subroutine lbe_sanity_check(N,insanity)
  type(lbe_site),intent(in) ::& 
       &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
  logical,intent(inout) :: insanity
  real(kind=rk),parameter :: max_vel=0.1_rk,small=1.0e-12_rk
  logical :: found_NaN,site_too_fast,tmps
  integer :: i,j,k,s,l,ierror
  real(kind=rk) :: tot
  real(kind=rk) :: nred,nblue,nsurf,rhored,rhoblue,rhosurf,sumred,sumblue&
       &,sumsurf,sumdens,rhotot
  real(kind=rk) :: p(3),sum_p(3)
  real(kind=rk), dimension(3) :: velred,velredavg,velredavgsum
  real(kind=rk), dimension(3) :: velblue,velblueavg,velblueavgsum
  real(kind=rk), dimension(3) :: velsurf,velsurfavg,velsurfavgsum
  real(kind=rk) :: maxvel(n_spec,3) ! maximum velocity for all species
  ! array dimensions: species, component of vel...max, component of
  ! lattice position where vel...max has its maximum value
  integer :: maxx(n_spec,3,3),buf_pos(3)
  ! counters might exceed 32bit integer range, MPI_INTEGER8 seems to
  ! be not available everywhere (?) so use double to be safe
  real(kind=rk) :: n_fs,n_nn,n_tf,sum_fs,sum_nn,sum_tf
  real(kind=rk) :: pair(2),rpair(2) ! for MPI_MAXLOC
  character(len=4) :: spec_name(3)=(/' red','blue','surf'/)
!#ifdef COMMON_VEL_FIX
#ifndef OLD_VEL
    real(kind=rk), dimension(1:3,0:n_spec) :: corr_vel
#endif
  call log_msg_hdr("Start LBE sanity check")

  ! initialize
  nred=0.0_rk
  nblue=0.0_rk
  nsurf=0.0_rk
  rhored=0.0_rk
  rhoblue=0.0_rk
  rhosurf=0.0_rk
  rhotot = 0.0_rk
  sumred=0.0_rk
  sumblue=0.0_rk
  sumsurf=0.0_rk
  sumdens=0.0_rk
  velredavg=0.0_rk
  velredavgsum=0.0_rk
  velblueavg=0.0_rk
  velblueavgsum=0.0_rk
  velsurfavg=0.0_rk
  velsurfavgsum=0.0_rk
  p = 0.0_rk
  sum_p = 0.0_rk
  found_NaN = .false.
  n_fs = 0.0_rk
  n_nn = 0.0_rk
  n_tf = 0.0_rk
  maxvel=0.0_rk
  maxx = 0

  do i = 1, nx
    do j = 1, ny
      do k = 1, nz
        no_rock: if ( is_fluid(N(i,j,k)%rock_state) ) then
          n_fs = n_fs+1.0_rk

          ! look for NaNs
          directions: do s = 1, nvecs
            if (check_NaN(N(i,j,k)%n_r(s))) found_NaN = .true.
#ifndef SINGLEFLUID
            if (check_NaN(N(i,j,k)%n_b(s))) found_NaN = .true.
#ifndef NOSURFACTANT
            if (check_NaN(N(i,j,k)%n_s(s))) found_NaN = .true.
#endif
#endif
          end do directions

          ! look for negative distribution functions
          if (any(N(i,j,k)%n_r<0.0_rk)&
#ifndef SINGLEFLUID
               &.or.any(N(i,j,k)%n_b<0.0_rk)&
#ifndef NOSURFACTANT
               &.or.any(N(i,j,k)%n_s<0.0_rk)&
#endif
#endif
               &) n_nn = n_nn+1.0_rk

          ! accumulate total mass
           rhotot = 0.0_rk
          rhored  = sum(N(i,j,k)%n_r*g)
          nred  = nred  + rhored
          rhotot = rhotot+rhored
#ifndef SINGLEFLUID
          rhoblue = sum(N(i,j,k)%n_b*g)
          nblue = nblue + rhoblue
          rhotot = rhotot+rhoblue
#endif
#ifndef NOSURFACTANT
          rhosurf = sum(N(i,j,k)%n_s*g)
          nsurf = nsurf + rhosurf
          rhotot = rhotot + rhosurf
#endif
          ! accumulate avg and max velocities
          site_too_fast = .false.
          if (abs(rhored)>small) then
            velred(1) = sum(N(i,j,k)%n_r*g*cx)*amass_r/rhored
            velred(2) = sum(N(i,j,k)%n_r*g*cy)*amass_r/rhored
            velred(3) = sum(N(i,j,k)%n_r*g*cz)*amass_r/rhored
            velredavg = velredavg+velred
            do l=1,3
              if (abs(velred(l))>=maxvel(1,l)) then
                maxvel(1,l) = abs(velred(l))
                maxx(1,l,:) = (/i,j,k/)+start-1
              end if
            end do
            if (norm(velred)>max_vel) site_too_fast = .true.
          end if
#ifndef SINGLEFLUID
          if (abs(rhoblue)>small) then
            velblue(1) = sum(N(i,j,k)%n_b*g*cx)*amass_b/rhoblue
            velblue(2) = sum(N(i,j,k)%n_b*g*cy)*amass_b/rhoblue
            velblue(3) = sum(N(i,j,k)%n_b*g*cz)*amass_b/rhoblue
            velblueavg = velblueavg+velblue
            do l=1,3
              if (abs(velblue(l))>=maxvel(2,l)) then
                maxvel(2,l) = abs(velblue(l))
                maxx(2,l,:) = (/i,j,k/)+start-1
              end if
            end do
            if (norm(velblue)>max_vel) site_too_fast = .true.
          end if
#endif
#ifndef NOSURFACTANT
          if (abs(rhosurf)>small) then
            velsurf(1) = sum(N(i,j,k)%n_s*g*cx)*amass_s/rhosurf
            velsurf(2) = sum(N(i,j,k)%n_s*g*cy)*amass_s/rhosurf
            velsurf(3) = sum(N(i,j,k)%n_s*g*cz)*amass_s/rhosurf
            velsurfavg = velsurfavg+velsurf
            do l=1,3
              if (abs(velsurf(l))>=maxvel(3,l)) then
                maxvel(3,l) = abs(velsurf(l))
                maxx(3,l,:) = (/i,j,k/)+start-1
              end if
            end do
            if (norm(velsurf)>max_vel) site_too_fast = .true.
          end if
#endif
          if (site_too_fast) n_tf = n_tf+1.0_rk

          ! accumulate total momentum
!#ifdef COMMON_VEL_FIX
#ifndef OLD_VEL

call corrected_velocity(N(i,j,k),lbe_force(:,:,i,j,k),corr_vel)
p = p + corr_vel(:,0)*rhotot
!if (i==1.AND.j ==1 .AND. k ==64) then
!    !print*, "density \n", rho,rho_test
!    print*, "mom \n",   p(1), p(2),p(3)
!READ*
!endif
#else
          p = p+massflow(N(i,j,k))
#endif
        end if no_rock
      end do
    end do
  end do

  if (found_NaN) then
    call log_msg("Found NaN",.true.)
    insanity = .true.
  end if

  ! check mass conservation:
  ! Now, add up all the subdomains' values of nred,,nblue,nsurf and place
  ! the sum in rank 0's bigsum.
  call MPI_Reduce(nred,sumred,1,LBE_REAL,MPI_SUM,0,Comm_Cart,ierror)
  msgstr = 'MPI_Reduce() of red failed'
  call checkmpi(ierror,msgstr)

#ifndef SINGLEFLUID
  call MPI_Reduce(nblue,sumblue,1,LBE_REAL,MPI_SUM,0,Comm_Cart,ierror)
  msgstr = 'MPI_Reduce() of blue failed'
  call checkmpi(ierror,msgstr)
#endif

#ifndef NOSURFACTANT
  call MPI_Reduce(nsurf,sumsurf,1,LBE_REAL,MPI_SUM,0,Comm_Cart,ierror)
  msgstr = 'MPI_Reduce() of surf failed'
  call checkmpi(ierror,msgstr)
#endif
  if (myrankc == 0) then
    tot=product(tsize)
    write(msgstr,"('Total densities (r,b,g)      :',3(2X,ES15.8))") &
      max(0.d-30,sumred), max(0.d-15,sumblue), max(0.d-30,sumsurf)
    call log_msg(msgstr)
    write(msgstr,"('Average densities (r,b,g)    :',3(2X,ES15.8))") &
      max(0.d-30,sumred/tot), max(0.d-30,sumblue/tot), max(0.d-30,sumsurf/tot)
    call log_msg(msgstr)
    sumdens=sumred+sumblue+sumsurf
    write(msgstr,"('Total/average density        :',2(2X,ES15.8))") &
      sumdens, max(0.d-30,sumdens/tot)
    call log_msg(msgstr)
    call log_ws()
  end if

  call MPI_Reduce(velredavg,velredavgsum,3,LBE_REAL,MPI_SUM,0,Comm_Cart,ierror)
  msgstr = 'MPI_Reduce() of avg red velocity failed'
  call checkmpi(ierror,msgstr)
#ifndef SINGLEFLUID
  call MPI_Reduce(velblueavg,velblueavgsum,3,LBE_REAL,MPI_SUM,0,Comm_Cart,ierror)
  msgstr = 'MPI_Reduce() of avg blue velocity failed'
  call checkmpi(ierror,msgstr)
#endif
#ifndef NOSURFACTANT
  call MPI_Reduce(velsurfavg,velsurfavgsum,3,LBE_REAL,MPI_SUM,0,Comm_Cart,ierror)
  msgstr = 'MPI_Reduce() of avg red velocity failed'
  call checkmpi(ierror,msgstr)
#endif

  ! determine position of maximum values of velocity components for
  ! all species present
  species: do i=1,n_spec
    vel_comp: do j=1,3
      ! If we want to use the predefined MPI_MAXLOC operation in
      ! Fortran we have to pack the rank into a floating point
      ! number here.
      pair = (/maxvel(i,j),real(myrankc,kind=rk)/)
      call MPI_Allreduce(pair,rpair,1,MPI_2DOUBLE_PRECISION,MPI_MAXLOC&
           &,comm_cart,ierror)
      ! rank holding the highest velocity sends it to rank 0
      if (myrankc==0) then
        maxvel(i,j) = rpair(1)
        if (nint(rpair(2))/=0) then ! else maxx(i,j,:) is on rank 0 already
          call MPI_Recv(buf_pos,3,LBE_REAL,nint(rpair(2)),0,comm_cart&
               &,MPI_STATUS_IGNORE,ierror)
          maxx(i,j,:) = buf_pos
        end if
      else if (myrankc==nint(rpair(2))) then
        buf_pos = maxx(i,j,:)
        call MPI_Send(buf_pos,3,LBE_REAL,0,0,comm_cart,ierror)
      end if
    end do vel_comp
  end do species

  myrankc0: if (myrankc == 0) then
    call log_msg("Maximum absolute velocity    : ")
    do i=1,n_spec
      write (msgstr,"('     ',A4,' (x,y,z)            :',3(2X,ES15.8))") &
           &spec_name(i),max(0.0e-30_rk,maxvel(i,:))
      call log_msg(msgstr)
      write (msgstr&
           &,"('      at position',12X,': ',3('(',2(I4,','),I4,')',:,X))") &
           &maxx(i,1,:),maxx(i,2,:),maxx(i,3,:)
      call log_msg(msgstr)
    end do

    call log_msg("Average velocity             : ")
    write(msgstr,"('      red (x,y,z)            :',3(2X,ES15.8))") &
         velredavgsum(1)/tot, velredavgsum(2)/tot, velredavgsum(3)/tot
    call log_msg(msgstr)
#ifndef SINGLEFLUID
    write(msgstr,"('     blue (x,y,z)            :',3(2X,ES15.8))") &
         velblueavgsum(1)/tot, velblueavgsum(2)/tot, velblueavgsum(3)/tot
    call log_msg(msgstr)
#endif
#ifndef NOSURFACTANT
    write(msgstr,"('     surf (x,y,z)            :',3(2X,ES15.8))") &
         velsurfavgsum(1)/tot, velsurfavgsum(2)/tot, velsurfavgsum(3)/tot
    call log_msg(msgstr)
#endif
    call log_ws()

    do i=1,n_spec
      if (any(maxvel(i,:)>0.1_rk)) then
        call log_msg("********  WARNING: MAXIMUM "//trim(spec_name(i))&
             &//" VELOCITY IS TOO HIGH!  ********")
!!$             insanity = .true.
      end if

    end do
  end if myrankc0

  ! total number of fluid sites
  call MPI_Reduce(n_fs,sum_fs,1,LBE_REAL,MPI_SUM,0,comm_cart,ierror)
  call checkmpi(ierror,'lbe_sanity_check(): MPI_Reduce() of n_fs failed')

  call MPI_Reduce(n_tf,sum_tf,1,LBE_REAL,MPI_SUM,0,comm_cart,ierror)
  call checkmpi(ierror,'lbe_sanity_check(): MPI_Reduce() of n_tf failed')
  if (myrankc==0 .and. sum_tf/=0.0_rk) then
    write (msgstr&
         &,"('********  WARNING: absolute species velocity higher than ',"&
         &//"ES15.8,' on ',I0,' lattice sites ( ',ES15.8,"&
         &//"' % of total fluid sites)  ********')") &
         &max_vel,int(sum_tf),100.0_rk*sum_tf/sum_fs
    call log_msg(msgstr)
  end if

  call MPI_Reduce(n_nn,sum_nn,1,LBE_REAL,MPI_SUM,0,comm_cart,ierror)
  call checkmpi(ierror,'lbe_sanity_check(): MPI_Reduce() of n_nn failed')
  if (myrankc==0 .and. sum_nn/=0.0_rk) then
    write (msgstr&
         &,"('********  WARNING: negative distribution functions on '"&
         &//",I0,' lattice sites ( ',ES15.8,"&
         &//"' % of total fluid sites)  ********')") &
         &int(sum_nn),100.0_rk*sum_nn/sum_fs
    call log_msg(msgstr)
  end if

  ! total momentum
  call MPI_Reduce(p,sum_p,3,LBE_REAL,MPI_SUM,0,comm_cart,ierror)
  call checkmpi(ierror,'lbe_sanity_check(): MPI_Reduce() of momentum failed')
  if (myrankc==0) then
    write (msgstr,"('Total fluid momentum (x,y,z) :',3(2X,ES15.8))") sum_p
    call log_msg(msgstr)
  end if

  write (msgstr,"('Total mass                   :',1(2X,ES15.8))") total_mass(N)
  call log_msg(msgstr)
  write (msgstr,"('Total avg velocity (x,y,z)   :',3(2X,ES15.8))") &
       &cached_avg_total_velocity(N)
  call log_msg(msgstr)

  tmps = insanity
  call MPI_Allreduce(tmps,insanity,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ierror)
  call checkmpi(ierror,'MPI_Allreduce() of insanity check failed')

  call log_msg_hdr("End LBE sanity check")
end subroutine lbe_sanity_check

!>Writes out an ASCII file detailing the position of each of
!>the processors in the virtual topology.
subroutine lbe_write_topology()
  integer, dimension(:,:),allocatable :: pcoords
  character(len=1024) :: filename
  integer :: i

  call lbe_make_filename_output(filename,'coords','.txt',nt)

  write(msgstr,"('wtfile=<',A,'>')") trim(filename)
  call log_msg(msgstr)

  open(unit=10,file=filename)
  allocate(pcoords(3,0:nprocs-1))
  call find_topology(pcoords)
  do i=nprocs-1,0,-1
    write(10,'(I6.6,a,I6.6,a,I6.6,a,I6.6)') &
         &i,' ',pcoords(1,i),' ',pcoords(2,i),' ',pcoords(3,i)
  end do
  deallocate(pcoords)
  close(unit=10)
end subroutine lbe_write_topology

!> This subroutine just calls the specific routines to dump the data to disk
subroutine dump_data(N)
  implicit none
  type(lbe_site), dimension(0:,0:,0:), intent(inout) :: N

  if (check_dump_now(sci_od, n_sci_od) ) then
    call log_msg("  Dumping OD...")
    call dump_od(N)
  end if
#ifndef SINGLEFLUID
  if (check_dump_now(sci_int, n_sci_int) ) then
    call log_msg("  Dumping COLOUR...")
    call dump_colour(N)
  end if
  if (check_dump_now(sci_colour_clusters, n_sci_colour_clusters) .or. &
      check_dump_now(sci_colour_clusters, n_sci_colour_clusters_index) ) then
    if (nt .ge. n_sci_start_colour_clusters) then
      call log_msg("  Dumping COLOUR CLUSTERS...")
      call dump_colour_clusters(N)
    end if
  end if
  if (check_dump_now(sci_wd, n_sci_wd) ) then
    call log_msg("  Dumping WD...")
    call dump_wd(N)
  end if
#endif
#ifndef NOSURFACTANT
  if (check_dump_now(sci_sur, n_sci_sur) ) then
    call log_msg("  Dumping SUR...")
    call dump_sur(N)
  end if
  if (check_dump_now(sci_dir, n_sci_dir) ) then
    call log_msg("  Dumping DIR...")
    call dump_dir(N)
  end if
#endif
  if (check_dump_now(sci_vel, n_sci_vel) ) then
    call log_msg("  Dumping VEL...")
    call dump_vel(N)
  end if
  if (check_dump_now(sci_flo, n_sci_flo) ) then
    call log_msg("  Dumping FLO...")
    call dump_flo(N)
  end if
  if (check_dump_now(sci_arrows, n_sci_arrows) .or. &
     check_dump_now(sci_velocities, n_sci_velocities) .or. &
     check_dump_now(sci_velocities_od, n_sci_velocities_od) .or. &
     check_dump_now(sci_velocities_wd, n_sci_velocities_wd) .or. &
     check_dump_now(sci_flux_od, n_sci_flux_od) .or. &
     check_dump_now(sci_flux_wd, n_sci_flux_wd) ) then
    call log_msg("  Dumping ARROWS / VELOCITIES...")
    call dump_arrows(N)
  end if
  if (check_dump_now(sci_rock, n_sci_rock) ) then
    call log_msg("  Dumping ROCK STATE...")
    call dump_rock(N)
  end if
#ifdef SINGLEFLUID
  if (SCMP) then
     if (check_dump_now(sci_rock_colour, n_sci_rock_colour) ) then
        call log_msg("  Dumping ROCK COLOUR...")
        call dump_rock_colour(N)
     end if
     if (check_dump_now(sci_rock_rho_r, n_sci_rock_rho_r) ) then
        call log_msg("  Dumping ROCK DENSITY R...")
        call dump_rock_rho_r(N)
     end if
  end if
#else
  if (check_dump_now(sci_rock_colour, n_sci_rock_colour) ) then
     call log_msg("  Dumping ROCK COLOUR...")
     call dump_rock_colour(N)
  end if
  if (check_dump_now(sci_rock_rho_r, n_sci_rock_rho_r) ) then
    call log_msg("  Dumping ROCK DENSITY R...")
    call dump_rock_rho_r(N)
  end if
  if (check_dump_now(sci_rock_rho_b, n_sci_rock_rho_b) ) then
    call log_msg("  Dumping ROCK DENSITY B...")
    call dump_rock_rho_b(N)
  end if
#endif

  ! Note: dump_pressure is not called if the flag NOSURFACTANT is defined
  ! This is automatically excluded by the init conditions (see lbe_init).
  if (check_dump_now(sci_pressure, n_sci_pressure) ) then
    call log_msg("  Dumping PRESSURE...")
    if(sci_pressure_init == 1) then
      call dump_pressure(N,droplet)
      call dump_pressure(N,nondroplet)
    else
      ! if((sci_pressure_init == 8) .or.                     &
      !   (sci_pressure_init == 3) .or.                     &
      !   (sci_pressure_init == -2)) then
      call dump_pressure(N,nondroplet)
      ! CALL dump_popul_zx(N)
    end if
  end if
end subroutine dump_data

  !> dumps output of  fluxz()  (see there) for every region
  !> specified by  fluxz_[xyz](lo|hi)  to a file named according
  !> to  fluxz_name
subroutine dump_fluxz(N, mass_flux)
  type(lbe_site),intent(in) :: N(0:,0:,0:)
  logical, intent(in) :: mass_flux
  integer :: i,dz,nf,hipos(3)
  real(kind=8) :: sf
#ifdef MD
  integer :: np
  real(kind=8) :: sp
#endif
  integer,parameter :: fz_file_unit=12
  character(len=1024) fz_file_name

  do i=1,fluxz_regions
    if (fluxz_xhi(i)==-1) cycle ! value -1 indicates unused slot

    ! value 0 means maximum position in respective direction
    hipos = (/fluxz_xhi(i),fluxz_yhi(i),fluxz_zhi(i)/)
    where (hipos==0) hipos = (/tnx,tny,tnz/)

    CALL fluxz(N, mass_flux, fluxz_xlo(i), fluxz_ylo(i), fluxz_zlo(i), &
      hipos(1), hipos(2), hipos(3), dz, nf, sf&
#ifdef MD
      &, np, sp&
#endif
      &)

    if (myrankc==0) then

      if (mass_flux) then
        CALL lbe_make_filename_output(fz_file_name,'massfluxz_'//trim(fluxz_name(i)),'.asc',nt)
      else
        CALL lbe_make_filename_output(fz_file_name,'fluxz_'//trim(fluxz_name(i)),'.asc',nt)
      end if
      open (unit=fz_file_unit,file=fz_file_name,status='REPLACE',action='WRITE',recl=80)
#ifdef MD
      write (unit=fz_file_unit,fmt='(SS,I9,X,SS,I4,X,I7,X,SP,ES15.8,X,SS,I7,X,SP,ES15.8)') nt, dz, nf, sf, np, sp
#else
      write (unit=fz_file_unit,fmt='(SS,I9,X,SS,I4,X,I7,X,SP,ES15.8)') nt, dz, nf, sf
#endif
      close (fz_file_unit)
    end if
  end do
end subroutine dump_fluxz

!> Returns in  dz  the dimension in z-direction, in  n[fp]  the number of
!> lattice nodes occupied by fluid/ladd particles, and in s[fp] the sum of
!> the z-velocities at all fluid/ladd particle sites for the cuboidal chunk
!> parallel to the coordinate axis spanned between the points  (xl,yl,zl)
!> and  (xh,yh,zh) .
!>
!> See documentation on how to use these values. If compiled
!> without  MD  set, only the fluid related data is calculated and returned.
!>
!> \warning possible particle rotations are not taken into account yet
subroutine fluxz(N,mass_flux, xl,yl,zl,xh,yh,zh,dz,nf,sf&
#ifdef MD
  &,np,sp&
#endif
  &)

  implicit none

  type(lbe_site),intent(in) :: N(0:,0:,0:)
  logical,intent(in) :: mass_flux
  integer,intent(in) :: xl,yl,zl,xh,yh,zh
  integer,intent(out) :: dz,nf
  real(kind=rk),intent(out) :: sf
#ifdef MD
  integer,intent(out) :: np
  real(kind=rk),intent(out) :: sp

  integer,parameter :: fzc = 2 ! fluid and ladd particles
#else
  integer,parameter :: fzc = 1 ! fluid only
#endif
  integer :: m(fzc),nsum(fzc)
  real(kind=rk) :: s(fzc),ssum(fzc)
  integer :: x,y,z,ierror
  integer :: minx(3),maxx(3)

!#ifdef COMMON_VEL_FIX
#ifndef OLD_VEL
  real(kind=rk),dimension(1:3,0:n_spec) :: corr_vel

! #else COMMON_VEL_FIX
#else 
  real(kind=rk), dimension(nx,ny,nz,3) :: f_r
#ifndef SINGLEFLUID
  real(kind=rk), dimension(nx,ny,nz,3) :: f_b
#ifndef NOSURFACTANT
  real(kind=rk), dimension(nx,ny,nz,3) :: f_s
#endif
#endif
! #endif COMMON_VEL_FIX
#endif

  m = 0
  s = 0.0_rk

  minx = max(1,(/xl,yl,zl/)+1-start)
  maxx = min((/nx,ny,nz/),(/xh,yh,zh/)+1-start)

  do x=minx(1),maxx(1)
    do y=minx(2),maxx(2)
      do z=minx(3),maxx(3)
        fluid_or_particle: if ( is_fluid(N(x,y,z)%rock_state) ) then

!#ifdef COMMON_VEL_FIX
#ifndef OLD_VEL
             call corrected_velocity(N(x,y,z),lbe_force(:,:,x,y,z),corr_vel)

             m(1) = m(1) + 1

             if (mass_flux) then
                s(1) = s(1) + amass_r * sum(g(:)*N(x,y,z)%n_r(:))* corr_vel(3,1)
#ifndef SINGLEFLUID
                s(1) = s(1) + amass_b * sum(g(:)*N(x,y,z)%n_b(:))* corr_vel(3,2)
#else
#ifndef NOSURFACTANT
                s(1) = s(1) + amass_s * sum(g(:)*N(x,y,z)%n_s(:))* corr_vel(3,3)
#endif
#endif
             else
                s(1) = s(1) + corr_vel(3,1)
#ifndef SINGLEFLUID
                s(1) = s(1) + corr_vel(3,2)
#else
#ifndef NOSURFACTANT
                s(1) = s(1) + corr_vel(3,3)
#endif
#endif
            end if
! #else COMMON_VEL_FIX
#else

#ifdef SINGLEFLUID
          call lbe_calculate_sc_forces(N,x,y,z,f_r)
#else
#ifdef NOSURFACTANT
          call lbe_calculate_sc_forces(N,x,y,z,f_b,f_r)
#else
          call lbe_calculate_sc_forces(N,x,y,z,f_b,f_r,f_s)
#endif
#endif

            m(1) = m(1) + 1

            if (mass_flux) then
              s(1) = s(1) + (&
                &amass_r*(sum(g(:)*N(x,y,z)%n_r(:)*cz(:))-&
                &sum(g(:)*N(x,y,z)%n_r(:))*&
                &((f_r(x,y,z,3)+g_accn)/2.0_rk))&
#ifndef SINGLEFLUID
                &+amass_b*(sum(g(:)*N(x,y,z)%n_b(:)*cz(:))-&
                &sum(g(:)*N(x,y,z)%n_b(:))*&
                &((f_b(x,y,z,3)+g_accn)/2.0_rk))&
#ifndef NOSURFACTANT
                &+amass_s*(sum(g(:)*N(x,y,z)%n_s(:)*cz(:))-&
                &sum(g(:)*N(x,y,z)%n_s(:))*&
                &((f_s(x,y,z,3)+g_accn)/2.0_rk))&
#endif
#endif
                 &)
          else

            s(1) = s(1) + ( ( &
                 &amass_r*(sum(g(:)*N(x,y,z)%n_r(:)*cz(:))-&
                 &sum(g(:)*N(x,y,z)%n_r(:))*&
                 &((f_r(x,y,z,3)+g_accn)/2.0_rk))&
#ifndef SINGLEFLUID
                 &+amass_b*(sum(g(:)*N(x,y,z)%n_b(:)*cz(:))-&
                 &sum(g(:)*N(x,y,z)%n_b(:))*&
                 &((f_b(x,y,z,3)+g_accn)/2.0_rk))&
#ifndef NOSURFACTANT
                 &+amass_s*(sum(g(:)*N(x,y,z)%n_s(:)*cz(:))-&
                 &sum(g(:)*N(x,y,z)%n_s(:))*&
                 &((f_s(x,y,z,3)+g_accn)/2.0_rk))&
#endif
#endif
                 & )/( &
                 &+sum(g(:)*N(x,y,z)%n_r(:))*amass_r&
#ifndef SINGLEFLUID
                 &+sum(g(:)*N(x,y,z)%n_b(:))*amass_b&
#ifndef NOSURFACTANT
                 &+sum(g(:)*N(x,y,z)%n_s(:))*amass_s&
#endif
#endif
                 &) )
          end if

! #endif COMMON_VEL_FIX
#endif

#ifdef MD
        else if ( is_colloid((N(x,y,z)%rock_state) ) ) then
!          call log_msg("In no fluid_or_particle")
          m(2) = m(2)+1
!          write(msgstr,"('rock_state: ',I0)") nint(N(x,y,z)%rock_state)
!          call log_msg(msgstr)
          s(2) = s(2)+P(Mii_map(uid2i,nint(N(x,y,z)%rock_state)))%v(3)
#endif
        end if fluid_or_particle
      end do
    end do
  end do

  CALL MPI_Reduce(m,nsum,fzc,MPI_INTEGER,MPI_SUM,0,comm_cart,ierror)
  CALL MPI_Reduce(s,ssum,fzc,LBE_REAL,MPI_SUM,0,comm_cart,ierror)
  dz = 1+zh-zl
  nf = nsum(1)
  sf = ssum(1)
#ifdef MD
  np = nsum(2)
  sp = ssum(2)
#endif
end subroutine fluxz

!> \{
!> \name Rock input routines
!>
!> These routines should only be called by rank zero, and are used to
!> read
!>in the entire rock description file.
!>
!>Each takes an array (\c rockery) of integers, the size of the global lattice,
!>which is filled with rock data, and the \c filename to read.
!>
!>There is a routine called \c read_rock_foo for eadch possible rock file
!>format \c foo.

!> Reads from a file in "all" format. Note that sparse rocks are
!> permitted: a file consisting of the line
!> 
!> \c 3 \c 4 \c 5 \c 1
!> 
!> would put a single stone at coordinates (3,5,1), and leave the rest of
!> the lattice sites free.
!>
!> After reading a line from the file the entry is passed to the according
!> CPU.
!!!!Do not use this function if compiled with -DMD, because it will
!!!!set rock_state of rock lattice to positive, equal to rock_colou 
subroutine read_rock_all_par(filename,N)
  type(lbe_site),dimension(0:,0:,0:) :: N
  integer :: x, y, z, rock
  integer :: funit = 10
  integer :: ios
  character(len=*), intent(in) :: filename
  ! integer :: np,pcoords(nd),buffer(4)
  integer, dimension(:,:),allocatable :: cpcoords

  ! Create Array containing a cpu's SX SY SZ
  allocate(cpcoords(0:nprocs-1,nd))
  call set_cpcoords(cpcoords)

  ! Open file and loop
  open(unit = funit, file = filename)
  do
    read(funit, *, iostat = ios) x, y, z, rock
    if (ios .ne. 0) exit
    ! Send rock state to appropriate CPU.
    if ( (x .le. tnx ) .and. ( y .le. tny ) .and. ( z .le. tnz ) ) then
      !print*,'Set ',x,y,z,' to ',stone,'.'
      call send_rock_bit(N, x, y, z, rock, cpcoords)
    end if
  end do
  close(unit = funit)
  deallocate(cpcoords)

  ! SEND SOMETHING TO RELEASE WAITING PROCESSORS
  call send_rock_release()

end subroutine read_rock_all_par
!> \}

!> \{
!> \name Debugging output routines
!>
!> These dump some part of the system state to a file, for debugging
!>purposes. Not intended for use in production code, but very useful to
!>see what's going on inside.

!> This routine is unique, and intended for debugging. It writes a file in
!> "all" format, ie containing one line per lattice site. Each line
!> begins with three integers giving the coordinates of that site in the
!> local lattice, and then six floating-point numbers for the red, blue
!> and surfactant densities and dipoles, followed by an integer containing
!> the rock state at that site.
! subroutine dump_sites_all(N)
! 	implicit none
! 	type(lbe_site), dimension(0:,0:,0:) :: N
! 	character(len=1024) :: filename
! 	real*8 :: nred,nblue,nsurf
! 	integer :: x,y,z,s,nxi,nyi,nzi
! 	real*8 :: red,surf,blue
! 
!         nxi = size(N,1)-2
!         nyi = size(N,2)-2
!         nzi = size(N,3)-2
! 
! 	call lbe_make_filename(filename,'sites','.all',nt)
! 	
! 	open(10,file=filename)
! 
! 	do z=1,nzi
! 	 do y=1,nyi
! 	  do x=1,nxi
! 		red=sum(N(x,y,z)%n_r(:)*g)
! #ifndef SINGLEFLUID
! 		blue=sum(N(x,y,z)%n_b(:)*g)
! #endif
! #ifndef NOSURFACTANT
! 		surf=sum(N(x,y,z)%n_s(:)*g)
! 		write(10,'(3i3,6f12.8,i3)')				&
! 			ccoords(1)*nxi+x,				&
! 			ccoords(2)*nyi+y,				&
! 			ccoords(3)*nzi+z,				&
! 			red,blue,surf,					&
! 			N(x,y,z)%d(:),					&
! 			int(N(x,y,z)%rock_state)
! #else
! #ifndef SINGLEFLUID
! 		write(10,'(3i3,2f12.8,i3)')				&
! 			ccoords(1)*nxi+x,				&
! 			ccoords(2)*nyi+y,				&
! 			ccoords(3)*nzi+z,				&
! 			red,blue,					&
! 			int(N(x,y,z)%rock_state)
! #else
! 		write(10,'(3i3,1f12.8,i3)')				&
! 			ccoords(1)*nxi+x,				&
! 			ccoords(2)*nyi+y,				&
! 			ccoords(3)*nzi+z,				&
! 			red,					&
! 			int(N(x,y,z)%rock_state)
! #endif
! #endif
! 	  end do
! 	 end do
! 	end do
! 
! 	close(10)
! 
! end subroutine dump_sites_all

!> This routine dumps a subdomain-plus-haloes worth of reals.
!> It's intended for debugging purposes - useful to examine the state of a
!> single CPU's subdomain.
!> 
!> It also writes an AVS field file for easy visualization.
! subroutine dump_debug_real_all(prefix,dump)
! 	implicit none
! 	real*8, dimension(0:,0:,0:) :: dump
! 	character(len=1024) :: filename
! 	character(len=*) :: prefix
! 	integer :: x,y,z,nxi,nyi,nzi
!         nxi = size(dump,1)-2
!         nyi = size(dump,2)-2
!         nzi = size(dump,3)-2
! 
! 	call lbe_make_filename(filename,prefix,'.all',nt)
! 	
! 	open(10,file=filename)
! 
! 	do z=0,nzi+1
! 	 do y=0,nyi+1
! 	  do x=0,nxi+1
! 		write(10,'(3i3,f24.12)')	x,y,z, dump(x,y,z)
! 	  end do
! 	 end do
! 	end do
! 
! 	close(10)
! 	!write(*,'(a,a,a)') 'Wrote <',trim(filename),'>'
! 
! 	! Write a field file, too.
!         call dump_avs_fld(prefix,nxi,nyi,nzi,int(1))
! 
! end subroutine dump_debug_real_all
!> \}

!> \{
!> \name Output routines - Science-specific routines
!>
!>These routines are responsible for writing the science outputs when
!>required, such as the colour at each site.
!>
!>If a given science output, \c foo, is required, then \p lbe.f90 makes
!>a call to the subroutine \c dump_foo.
!>
!>\c dump_foo then compiles all the data required into one or more arrays,
!>and passes these to \c dump_scalar, \c dump_nscalar, \c dump_vector, etc,
!>depending on how many arrays are to be dumped and of what type.
!>
!>These routines in turn call \c dump_scalar_all, \c dump_scalar_bin,
!>etc depending on whether binary or ASCII output has been specified.
!>
!>Each of these routines is responsible for a given sort of science
!>output. It calls the appropriate array routine, which in turn calls the
!>appropriate file-format routine.

#ifndef SINGLEFLUID
!>Calculates the colour field at each site, and then calls \c dump_scalar.
subroutine dump_colour(N)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N
  real*8, dimension(:,:,:), allocatable :: scalar
  integer :: ierror,x,y,z,nxi,nyi,nzi

  nxi = size(N,1)-2
  nyi = size(N,2)-2
  nzi = size(N,3)-2

  allocate(scalar(1:nxi,1:nyi,1:nzi),stat=ierror)

  if (ierror .ne. 0) then
    call log_msg("WARNING: unable to allocate scalar buffer for colour output",.true.)
    return ! Inability to write output is not quite fatal.
  end if

  do z=1,nzi
    do y=1,nyi
      do x=1,nxi
        if ( is_fluid( N(x,y,z)%rock_state) ) then
          scalar(x,y,z)=sum(N(x,y,z)%n_r(:)*g)-sum(N(x,y,z)%n_b(:)*g)
        else
          scalar(x,y,z)=0.0
        endif
      end do
    end do
  end do
  CALL dump_scalar(scalar,'colour')
  deallocate(scalar)
end subroutine dump_colour

!>Calculates the colour field at each site, and then calls \c dump_scalar.
subroutine dump_colour_clusters(N)
  implicit none
  type(lbe_site), dimension(0:,0:,0:), intent(inout) :: N

  real(kind=rk), dimension(:,:,:), allocatable :: colour
  integer, dimension(:,:,:), allocatable :: index_field

  ! The hoshen_kopelman subroutine requires an upper and a lower limit to determine if a site
  ! is part of a cluster. We choose large_value such that a colour field value will never reach it.
  ! Similarly, non-fluid sites, which should never be part of a cluster, get this as a negative
  ! value as a choice for a value < 0.
  real(kind=rk), parameter :: large_value = 42.0_rk
  integer :: nxi,nyi,nzi
  integer :: x, y, z
  integer :: ierror
  integer :: n_clusters

  nxi = nx + 2
  nyi = ny + 2
  nzi = nz + 2

  allocate(colour(1:nxi,1:nyi,1:nzi),stat=ierror)
  if (ierror .ne. 0) then
    call log_msg("WARNING: unable to allocate scalar buffer for colour clusters.",.true.)
    return ! Inability to write output is not quite fatal.
  end if

  ! Loop over everything, including width 1 halo
  do z=0,nz+1
    do y=0,ny+1
      do x=0,nx+1
        if ( is_fluid( N(x,y,z)%rock_state) ) then
          colour(x+1,y+1,z+1) = sum(N(x,y,z)%n_r(:)*g) - sum(N(x,y,z)%n_b(:)*g)
        else
          colour(x+1,y+1,z+1) = -large_value
        endif
      end do
    end do
  end do

  allocate(index_field(1:nxi,1:nyi,1:nzi),stat=ierror)
  if (ierror .ne. 0) then
    call log_msg("WARNING: unable to allocate scalar buffer for index field.",.true.)
    return ! Inability to write output is not quite fatal.
  end if

  call log_msg("    Calling Hoshen Kopelman algorithm for colour field ...")
  call hoshen_kopelman(colour, index_field, 0.0_rk, large_value, parallelize = .true., n_clusters = n_clusters)

  if ( check_dump_now(sci_colour_clusters, n_sci_colour_clusters_index) ) then
    call log_msg("    Dumping cluster index scalar field ...")
    call dump_iscalar(index_field(2:nx+1,2:ny+1,2:nz+1),'cluster-index')
  end if

  if ( check_dump_now(sci_colour_clusters, n_sci_colour_clusters) ) then
    call log_msg("    Dumping cluster sizes ...")
    call dump_colour_cluster_sizes(index_field(2:nx+1,2:ny+1,2:nz+1), n_clusters)
  end if

  deallocate(index_field, stat=ierror)
  call check_allocate(ierror, "Failed deallocation of index_field")
  deallocate(colour, stat=ierror)
  call check_allocate(ierror, "Failed deallocation of colour")

end subroutine dump_colour_clusters

subroutine dump_colour_cluster_sizes(p_index_field, n_clusters)
  implicit none
  integer, dimension(1:,1:,1:), intent(inout) :: p_index_field
  integer, intent(in) :: n_clusters

  integer, dimension(:), allocatable :: local_csizes, global_csizes

  integer, parameter :: funit = 10

  character(len=1024) :: filename

  integer :: x, y, z ! lattice coordinates
  integer :: i
  integer :: idx
  integer :: dim_x, dim_y, dim_z ! array dimensions
  integer :: mpierror
  integer :: ierror

  ! Find array dimensions.
  dim_x = size(p_index_field, 1)
  dim_y = size(p_index_field, 2)
  dim_z = size(p_index_field, 3)

  allocate( local_csizes(0:n_clusters), stat=ierror )
  if (ierror .ne. 0) then
    call log_msg("WARNING: unable to allocate local_csizes buffer for colour clusters.",.true.)
    return ! Inability to write output is not quite fatal.
  end if

  allocate( global_csizes(0:n_clusters), stat=ierror )
  if (ierror .ne. 0) then
    call log_msg("WARNING: unable to allocate global_csizes buffer for colour clusters.",.true.)
    return ! Inability to write output is not quite fatal.
  end if

  local_csizes(:) = 0
  global_csizes(:) = 0

  !> Count cluster sites
  do x = 1, dim_x
    do y = 1, dim_y
      do z = 1, dim_z
        idx = p_index_field(x,y,z)
        if ( idx .lt. 0 .or. idx .gt. n_clusters ) then
          write(msgstr,"('Over- / underflow in n_clusters = ',I0,' , index ', I0, ' found at ',3(I0,X))") n_clusters, idx, x, y, z
          call error(msgstr)
        end if
        local_csizes( idx ) = local_csizes( idx ) + 1
      end do
    end do
  end do

  !> Parallel sum
  call MPI_Allreduce( local_csizes, global_csizes, n_clusters + 1, MPI_INTEGER, MPI_SUM, comm_cart, mpierror)
  ! call MPI_Barrier(comm_cart, mpierror)

  if (myrankc == 0) then
    !> Prepare to write file
    call lbe_make_filename_output(filename, 'cluster-sizes', '.asc', nt)
    ! write(msgstr,"('Writing cluster data to <',A,'>')") trim(filename) ; call log_msg(msgstr)

    open(unit = funit, file = filename)

    !> Write header
    write(funit, "(A6,X,A12)") "#?  id", "size"
    !> Dump data
    do i = 0, n_clusters
      write(funit,"(I6,X,I12)") i, global_csizes(i)
    end do

    close(unit = funit)
  end if

  deallocate( global_csizes, stat=ierror )
  call check_allocate(ierror, "Failed deallocation of global_csizes")
  deallocate( local_csizes, stat=ierror )
  call check_allocate(ierror, "Failed deallocation of local_csizes")

end subroutine dump_colour_cluster_sizes
! endif nSINGLEFLUID
#endif

#ifndef NOSURFACTANT
!> Calculates the surfactant density at each site, then calls \c dump_scalar.
subroutine dump_sur(N)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N
  real*8, dimension(:,:,:), allocatable :: scalar
  integer :: ierror,x,y,z,nxi,nyi,nzi

  nxi = size(N,1)-2
  nyi = size(N,2)-2
  nzi = size(N,3)-2

  allocate(scalar(1:nxi,1:nyi,1:nzi),stat=ierror)

  if (ierror .ne. 0) then
    call log_msg("WARNING: unable to allocate scalar buffer for surfactant output",.true.)
    return ! Inability to write output is not quite fatal.
  end if

  do z=1,nzi
    do y=1,nyi
      do x=1,nxi
        scalar(x,y,z)=amass_s*sum(N(x,y,z)%n_s(:)*g)
      end do
    end do
  end do
  call dump_scalar(scalar,'sur')
  deallocate(scalar)
end subroutine dump_sur
! endif NOSURFACTANT
#endif

#ifndef NOSURFACTANT
!> 
!> Calculates the dipole moment at each site, then calls \c dump_vector.
subroutine dump_dir(N)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N
  real*8, dimension(:,:,:,:), allocatable :: d
  integer :: ierror,x,y,z,nxi,nyi,nzi

  nxi = size(N,1)-2
  nyi = size(N,2)-2
  nzi = size(N,3)-2

  allocate(d(3,1:nxi,1:nyi,1:nzi),stat=ierror)
  if (ierror .ne. 0) then
    call log_msg("WARNING: unable to allocate scalar buffer for dir output",.true.)
    return
  end if

  do z=1,nzi
    do y=1,nyi
      do x=1,nxi
        d(:,x,y,z)=N(x,y,z)%d(:)
      end do
    end do
  end do
  CALL dump_vector(d,'dir')
  deallocate(d)
end subroutine dump_dir

#endif

!> Calculates the averaged Z velocity at each site, then calls
!> \c dump_scalar.
subroutine dump_vel(N)
  type(lbe_site),dimension(0:,0:,0:),intent(in) :: N
  real*8, dimension(:,:,:), allocatable :: scalar
  integer :: ierror,x,y,z,nxi,nyi,nzi,tz,tx,ty
!#ifdef COMMON_VEL_FIX
#ifndef OLD_VEL
  real(kind=rk),dimension(1:3,0:n_spec) :: corr_vel
! #else COMMON_VEL_FIX
#else
  real*8,dimension(:,:,:,:),allocatable :: f_r
#ifndef SINGLEFLUID
  real*8,dimension(:,:,:,:),allocatable :: f_b
#endif
#ifndef NOSURFACTANT
  real*8,dimension(:,:,:,:),allocatable :: f_s
#endif
! #endif COMMON_VEL_FIX
#endif
  nxi = size(N,1)-2
  nyi = size(N,2)-2
  nzi = size(N,3)-2
  allocate(scalar(1:nxi,1:nyi,1:nzi),stat=ierror)
  if (ierror .ne. 0) then
     call log_msg("WARNING: unable to allocate scalar buffer for vel output"&
          &,.true.)
     return	! Inability to write output is not quite fatal.
  end if
!#ifndef COMMON_VEL_FIX
#ifdef OLD_VEL
  allocate(f_r(nxi,nyi,nzi,3)&
#ifndef SINGLEFLUID
       &,f_b(nxi,nyi,nzi,3)&
#endif
#ifndef NOSURFACTANT
       &,f_s(nxi,nyi,nzi,3)&
#endif
       &,stat=ierror)
  if (ierror/=0) then
     call log_msg("WARNING: unable to allocate force buffer for vel output"&
          &,.true.)
     return
  end if
! #endifn COMMON_VEL_FIX
#endif

  do z=1,nzi
    do y=1,nyi
      do x=1,nxi
        ! Note this is u not u_tilde (no funny averaging with
        ! \tau If all your \tau's are the same then u = u_tilde,
        ! if not and you want u_tilde, strip the code out of
        ! lbe_collision.F90   (before 22.05.2015)
         !!!!Above statement is wrong. u != u_tilde, because u is after
         !collision step and u_tilde is before collision step. (after 22.05.2015)
        ! Get Shan Chen forces
        ! only calculate forces on fluid at non-rock sites
        if ( is_fluid(N(x,y,z)%rock_state) ) then

!#ifdef COMMON_VEL_FIX
#ifndef OLD_VEL
           call corrected_velocity(N(x,y,z),lbe_force(:,:,x,y,z),corr_vel)
           scalar(x,y,z)=corr_vel(3,0)
           
! #else COMMON_VEL_FIX
#else
#ifdef SINGLEFLUID
          CALL lbe_calculate_sc_forces(N,x,y,z,f_r)
#else
#ifdef NOSURFACTANT
          CALL lbe_calculate_sc_forces(N,x,y,z,f_b,f_r)
#else
          CALL lbe_calculate_sc_forces(N,x,y,z,f_b,f_r,f_s)
#endif
#endif
        else
          f_r(x,y,z,:) = 0.d0
#ifndef SINGLEFLUID
          f_b(x,y,z,:) = 0.d0
#endif
#ifndef NOSURFACTANT
          f_s(x,y,z,:) = 0.d0
#endif
        end if

        if ( is_fluid(N(x,y,z)%rock_state) ) then

           scalar(x,y,z) = sum(N(x,y,z)%n_r(:)*g*cz)*amass_r*get_omega_r(N, x, y, z)
           if (SCMP) then
             scalar(x,y,z) = sum(N(x,y,z)%n_r(:)*g*cz)*amass_r*get_omega_r(N, x, y, z) - tau_r*f_r(x,y,z,3)/2.d0 
           end if
#ifndef SINGLEFLUID
           scalar(x,y,z) = scalar(x,y,z) + sum(N(x,y,z)%n_b(:)*g*cz)*amass_b*get_omega_b(N, x, y, z) !- tau_b*f_b(x,y,z,3)/2.d0
#endif
#ifndef NOSURFACTANT
           scalar(x,y,z) = scalar(x,y,z) + sum(N(x,y,z)%n_s(:)*g*cz)*amass_s*get_omega_s(N, x, y, z) !- tau_s*f_s(x,y,z,3)/2.d0
           scalar(x,y,z) = scalar(x,y,z) / &
                & max(10.D-9,dble( sum(N(x,y,z)%n_r(:)*g)*amass_r*get_omega_r(N, x, y, z) &
                + sum(N(x,y,z)%n_b(:)*g)*amass_b*get_omega_b(N, x, y, z) + sum(N(x,y,z)%n_s(:)*g)*amass_s*get_omega_s(N, x, y, z)) )
#else 
#ifndef SINGLEFLUID
           scalar(x,y,z) = scalar(x,y,z) / &
                & max(10.D-9,dble( sum(N(x,y,z)%n_r(:)*g)*amass_r*get_omega_r(N, x, y, z)  + sum(N(x,y,z)%n_b(:)*g)*amass_b*get_omega_b(N, x, y, z) ))
#else
           scalar(x,y,z) = scalar(x,y,z) / max(10.D-9,dble( sum(N(x,y,z)%n_r(:)*g)*amass_r*get_omega_r(N, x, y, z) ))
#endif
#endif
           
           ! force extra-term
           tx = x + ccoords(1)*nx
           ty = y + ccoords(2)*ny
           tz = z + ccoords(3)*nz
           if( (tz.ge.g_accn_min).and.(tz.le.g_accn_max) .and. &
                (tx.ge.g_accn_min_x).and.(tx.le.g_accn_max_x) .and. &
                (ty.ge.g_accn_min_y).and.(ty.le.g_accn_max_y)) then
              if ( is_fluid(N(x,y,z)%rock_state) ) then
                 ! external force term, in the zone where the
                 ! external force works
                 scalar(x,y,z) = scalar(x,y,z) - g_accn/2.0d0
              end if
           endif
! #endif COMMON_VEL_FIX
#endif
        else
           ! rock_site
           scalar(x,y,z) = 0.0d0
        end if

      end do
    end do
  end do
  call dump_scalar(scalar,'vel')
  deallocate(scalar)

!#ifndef COMMON_VEL_FIX
#ifdef OLD_VEL
  deallocate(f_r&
#ifndef SINGLEFLUID
       &,f_b&
#endif
#ifndef NOSURFACTANT
       &,f_s&
#endif
       &,stat=ierror)
#endif

end subroutine dump_vel

!> Calculates the oil density at each site, then calls
!> \c dump_scalar.
subroutine dump_od(N)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N
  real*8, dimension(:,:,:), allocatable :: oil
  integer :: ierror,x,y,z,nxi,nyi,nzi

  nxi = size(N,1)-2
  nyi = size(N,2)-2
  nzi = size(N,3)-2

  allocate(oil(1:nxi,1:nyi,1:nzi),stat=ierror)

  if (ierror .ne. 0) then
    call log_msg("WARNING: unable to allocate scalar buffer for od output",.true.)
    return	! Inability to write output is not quite fatal.
  end if

  do z=1,nzi
    do y=1,nyi
      do x=1,nxi
        if ( is_fluid(N(x,y,z)%rock_state) ) then
          oil(x,y,z)=amass_r*sum(N(x,y,z)%n_r(:)*g)
        else
          oil(x,y,z)=0.0
        endif
      end do
    end do
  end do
  call dump_scalar(oil,'od')
  deallocate(oil)
end subroutine dump_od


!> Calculates the oil density at each site, then calls
!> \c dump_scalar.
subroutine dump_vel_probe(N, pos)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N
  integer, dimension(3), intent(in)   :: pos
!  integer, intent(in)                 :: nprobe
!  character(len=*), intent(in) :: filename
  real(kind=rk),dimension(1:3,0:n_spec) :: corr_vel
  integer             :: lpos(3),minpos(3),maxpos(3),ierror
!  integer             :: status(MPI_STATUS_SIZE)
!  character(len=512)  :: buffer,suffix
!  integer,parameter   :: fz_file_unit=121
#ifndef SINGLEFLUID
  call error('DUMP_VEL_PROBE implemented only for SINGLE_FLUID')
#endif
!  suffix = ".asc"
!  write(buffer,"('./',A,'/'A,A)") &
!    trim(folder), trim(filename), trim(suffix) 
  
  if(pos(1).lt.1.or.pos(2).lt.1.or.pos(3).lt.1.or.&
    &pos(1).gt.tnx.or.pos(2).gt.tny.or.pos(3).gt.tnz) then
     write(msgstr,"('Probe location = ',3I5,' is out of domain : ', 3I5)") pos,&
     &tnx,tny,tnz
     call error(msgstr)
  end if
!  minpos = ccoords*(/nx,ny,nz/)
!  maxpos = (ccoords+1)*(/nx,ny,nz/)
!  call MPI_FILE_OPEN(MPI_COMM_WORLD, 'testfile', &  
!                     MPI_MODE_WRONLY + MPI_MODE_CREATE, &  
!                     MPI_INFO_NULL, thefile, ierr)  
    lpos(1) = pos(1) - nx*ccoords(1)
    lpos(2) = pos(2) - ny*ccoords(2)
    lpos(3) = pos(3) - nz*ccoords(3)

  if(lpos(1).gt.0 .and.lpos(2).gt.0 .and.lpos(3).gt.0 .and.&
     lpos(1).le.nx.and.lpos(2).le.ny.and.lpos(3).le.nz)then
    ! compute local position from global
!    lpos(1) = pos(1) - nx*ccoords(1)
!    lpos(2) = pos(2) - ny*ccoords(2)
!    lpos(3) = pos(3) - nz*ccoords(3)
    call corrected_velocity(N(lpos(1),lpos(2),lpos(3)),&
                          &lbe_force(:,:,lpos(1),lpos(2),lpos(3)),corr_vel)

    write(*,"('Probe = ',3I5,(1X,I10.10),(SP,3(X,ES15.8)))") pos,nt,corr_vel(:,0)
!    call MPI_FILE_WRITE(thefile, buf, BUFSIZE, MPI_INTEGER, & 
!    MPI_STATUS_IGNORE, ierr)
!    call MPI_Sendrecv&
!    &(corr_vel,3,LBE_REAL,      0,1&
!    &,corr_vel,3,LBE_REAL,myrankc,1&
!    &,comm_cart,status,ierror)

!    call checkmpi(ierror,'failed to complete sendrecv')

!    open (unit=fz_file_unit,file=buffer, position="append")
!    write (unit=fz_file_unit,fmt='(I9,3F16.10)') nt,corr_vel(:,0)
!    close (fz_file_unit)
  endif     
!  if(myrankc==0) then
!    write(*,"('Probe = ',3I5,' velocity : ', 3F16.10)") pos,corr_vel
!  end if
!  call MPI_FILE_CLOSE(thefile, ierr) 


end subroutine dump_vel_probe

#ifndef SINGLEFLUID
!> 
!> Calculates water density at each site, then calls
!> \c dump_scalar.
subroutine dump_wd(N)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N
  real*8, dimension(:,:,:), allocatable :: water
  integer :: ierror,x,y,z,nxi,nyi,nzi

  nxi = size(N,1)-2
  nyi = size(N,2)-2
  nzi = size(N,3)-2

  allocate(water(1:nxi,1:nyi,1:nzi),stat=ierror)

  if (ierror .ne. 0) then
    call log_msg("WARNING: unable to allocate scalar buffer for wd output",.true.)
    return	! Inability to write output is not quite fatal.
  end if

  do z=1,nzi
    do y=1,nyi
      do x=1,nxi
        if ( is_fluid(N(x,y,z)%rock_state) ) then
          water(x,y,z)=amass_b*sum(N(x,y,z)%n_b(:)*g)
        else
          water(x,y,z)=0.0
        endif
      end do
    end do
  end do
  call dump_scalar(water,'wd')
  deallocate(water)
end subroutine dump_wd
#endif

!> Calculates the Z velocity of each species at each site, then calls
!> \c dump_3scalar.
subroutine dump_flo(N)
  type(lbe_site),dimension(0:,0:,0:),intent(in) :: N
  real*8, dimension(:,:,:), allocatable :: oil,water,surf
  integer :: ierror,x,y,z,nxi,nyi,nzi,tz,tx,ty
!#ifdef COMMON_VEL_FIX
#ifndef OLD_VEL
  real(kind=rk), dimension(1:3,0:n_spec) :: corr_vel
! #else COMMON_VEL_FIX
#else

  real*8,dimension(:,:,:,:),allocatable :: f_r
#ifndef SINGLEFLUID
  real*8,dimension(:,:,:,:),allocatable :: f_b
#endif
#ifndef NOSURFACTANT
  real*8,dimension(:,:,:,:),allocatable :: f_s
#endif

!#endif COMMON_VEL_FIX
#endif

  nxi = size(N,1)-2
  nyi = size(N,2)-2
  nzi = size(N,3)-2
  
!#ifdef COMMON_VEL_FIX
#ifndef OLD_VEL

  allocate(oil(1:nxi,1:nyi,1:nzi),stat=ierror)
  if (ierror .ne. 0) then
     call log_msg("WARNING: unable to allocate scalar buffer for flo oil output",.true.)
     return  ! Inability to write output is not quite fatal.
  end if
#ifndef SINGLEFLUID
  allocate(water(1:nxi,1:nyi,1:nzi),stat=ierror)
  if (ierror .ne. 0) then
     call log_msg("WARNING: unable to allocate scalar buffer for flo water output",.true.)
     return  ! Inability to write output is not quite fatal.
  end if
#endif
#ifndef NOSURFACTANT
  allocate(surf(1:nxi,1:nyi,1:nzi),stat=ierror)
  if (ierror .ne. 0) then
     call log_msg("WARNING: unable to allocate scalar buffer for flo surf output",.true.)
     return  ! Inability to write output is not quite fatal.
  end if
#endif  

  do z=1,nzi
     do y=1,nyi
        do x=1,nxi
           call corrected_velocity(N(x,y,z),lbe_force(:,:,x,y,z),corr_vel)
           oil(x,y,z) = corr_vel(3,1)
#ifndef SINGLEFLUID
           water(x,y,z) = corr_vel(3,2)
#endif
#ifndef NOSURFACTANT
           surf(x,y,z) = corr_vel(3,3)
#endif           
        end do
     end do
  end do

  call dump_scalar(oil,'flooil')
  deallocate(oil)
#ifndef SINGLEFLUID
  call dump_scalar(water,'flowater')
  deallocate(water)
#endif
#ifndef NOSURFACTANT
  call dump_scalar(surf,'flosurf')
  deallocate(surf)
#endif

! #else COMMON_VEL_FIX
#else

  allocate(f_r(nxi,nyi,nzi,3)&
#ifndef SINGLEFLUID
       &,f_b(nxi,nyi,nzi,3)&
#endif
#ifndef NOSURFACTANT
       &,f_s(nxi,nyi,nzi,3)&
#endif
       &,stat=ierror)
  if (ierror/=0) then
     call log_msg("WARNING: unable to allocate force buffer for flo output"&
          &,.true.)
     return
  end if

  allocate(oil(1:nxi,1:nyi,1:nzi),stat=ierror)
  if (ierror .ne. 0) then
    call log_msg("WARNING: unable to allocate scalar buffer for flo oil output",.true.)
    return  ! Inability to write output is not quite fatal.
  end if
  do z=1,nzi
    do y=1,nyi
      do x=1,nxi
        ! Get Shan Chen forces
        ! only calculate forces on fluid at non-rock sites
        if ( is_fluid(N(x,y,z)%rock_state) ) then
#ifdef SINGLEFLUID
          CALL lbe_calculate_sc_forces(N,x,y,z,f_r)
#else
#ifdef NOSURFACTANT
          CALL lbe_calculate_sc_forces(N,x,y,z,f_b,f_r)
#else
          CALL lbe_calculate_sc_forces(N,x,y,z,f_b,f_r,f_s)
#endif
#endif
        else
          f_r(x,y,z,:) = 0.d0
#ifndef SINGLEFLUID
          f_b(x,y,z,:) = 0.d0
#endif
#ifndef NOSURFACTANT
          f_s(x,y,z,:) = 0.d0
#endif
        end if
        oil(x,y,z)=sum(N(x,y,z)%n_r(:)*g*cz)/max(10.e-5,real(sum(N(x,y,z)%n_r(:)*g)))-f_r(x,y,z,3)/2.0d0
        tz=z+ccoords(3)*nz
        tx=x+ccoords(1)*nx
        ty=y+ccoords(2)*ny
        if( (tz.ge.g_accn_min).and.(tz.le.g_accn_max) .and. &
            (tx.ge.g_accn_min_x).and.(tx.le.g_accn_max_x) .and. &
            (ty.ge.g_accn_min_y).and.(ty.le.g_accn_max_y)) then
          if ( is_fluid( N(x,y,z)%rock_state ) ) then
            ! external force term, in the zone where the external force works
            oil(x,y,z) = oil(x,y,z)-g_accn/2.0d0
          end if
        endif
      end do
    end do
  end do
  call dump_scalar(oil,'flooil')
  deallocate(oil)

#ifndef SINGLEFLUID
  allocate(water(1:nxi,1:nyi,1:nzi),stat=ierror)
  if (ierror .ne. 0) then
    call log_msg("WARNING: unable to allocate scalar buffer for flo water output",.true.)
    return  ! Inability to write output is not quite fatal.
  end if
  do z=1,nzi
    do y=1,nyi
      do x=1,nxi
        water(x,y,z)=sum(N(x,y,z)%n_b(:)*g*cz)/max(10.e-5,real(sum(N(x,y,z)%n_b(:)*g)))-f_b(x,y,z,3)/2.0d0
        tz=z+ccoords(3)*nz
        tx=x+ccoords(1)*nx
        ty=y+ccoords(2)*ny
        if( (tz.ge.g_accn_min).and.(tz.le.g_accn_max) .and. &
            (tx.ge.g_accn_min_x).and.(tx.le.g_accn_max_x) .and. &
            (ty.ge.g_accn_min_y).and.(ty.le.g_accn_max_y)) then
          if ( is_fluid(N(x,y,z)%rock_state) ) then
            ! external force term, in the zone where the external force works
            water(x,y,z) = water(x,y,z)-g_accn/2.0d0
          end if
        endif
      end do
    end do
  end do
  call dump_scalar(water,'flowater')
  deallocate(water)
#endif
#ifndef NOSURFACTANT
  allocate(surf(1:nxi,1:nyi,1:nzi),stat=ierror)
  if (ierror .ne. 0) then
    call log_msg("WARNING: unable to allocate scalar buffer for flo surf output",.true.)
    return  ! Inability to write output is not quite fatal.
  end if
  do z=1,nzi
    do y=1,nyi
      do x=1,nxi
        surf(x,y,z)=sum(N(x,y,z)%n_s(:)*g*cz)/max(10.e-5,real(sum(N(x,y,z)%n_s(:)*g)))-f_s(x,y,z,3)/2.0d0
        tz=z+ccoords(3)*nz
        tx=x+ccoords(1)*nx
        ty=y+ccoords(2)*ny
        if( (tz.ge.g_accn_min).and.(tz.le.g_accn_max) .and. &
            (tx.ge.g_accn_min_x).and.(tx.le.g_accn_max_x) .and. &
            (ty.ge.g_accn_min_y).and.(ty.le.g_accn_max_y)) then
          if ( is_fluid( N(x,y,z)%rock_state )) then
            ! external force term, in the zone where the external force works
            surf(x,y,z) = surf(x,y,z)-g_accn/2.0d0
          end if
        endif
      end do
    end do
  end do
  CALL dump_scalar(surf,'flosurf')
  deallocate(surf)

  deallocate(f_r&
#ifndef SINGLEFLUID
       &,f_b&
#endif
#ifndef NOSURFACTANT
       &,f_s&
#endif
       &,stat=ierror)

!#endif COMMON_VEL_FIX
#endif
#endif

end subroutine dump_flo

!>Calculates the nett flow direction at each site, then calls \c dump_vector.
!>
!>Added by Nelido, edited by Jens 26.06.02
!>Edited by FrankR 19.07.07 to calculate the averaged velocity, not the momentum
subroutine dump_arrows(N)
  type(lbe_site),dimension(0:,0:,0:),intent(in) :: N
  real*8, dimension(:,:,:,:), allocatable :: v
  real*8, dimension(:,:,:,:), allocatable :: vod
  real*8, dimension(:,:,:,:), allocatable :: vwd
  real*8, dimension(:,:,:,:), allocatable :: fluxod
  real*8, dimension(:,:,:,:), allocatable :: fluxwd
  !real*8, dimension(:,:,:,:), allocatable :: vpre
  integer :: ierror,x,y,z,nxi,nyi,nzi,tx,ty,tz, ierror_wd,ierror_od,ierror_fod,ierror_fwd
!#ifdef COMMON_VEL_FIX
#ifndef OLD_VEL
  real(kind=rk), dimension(1:3,0:n_spec) :: corr_vel
  real(kind=rk), dimension(1:3,0:n_spec) :: corr_flux
  !real(kind=rk), dimension(1:3,0:n_spec) :: corr_vpre
! #else COMMON_VEL_FIX
#else
  real*8,dimension(:,:,:,:),allocatable :: f_r
#ifndef SINGLEFLUID
  real*8,dimension(:,:,:,:),allocatable :: f_b
#endif
#ifndef NOSURFACTANT
  real*8,dimension(:,:,:,:),allocatable :: f_s
#endif
! #endif COMMON_VEL_FIX
#endif

  nxi = size(N,1)-2
  nyi = size(N,2)-2
  nzi = size(N,3)-2
  allocate(v(3,1:nxi,1:nyi,1:nzi),stat=ierror)
  if (ierror .ne. 0) then
     call log_msg("WARNING: unable to allocate scalar buffer for arrows output"&
          &,.true.)
     return
  end if
 allocate(vod(3,1:nxi,1:nyi,1:nzi),stat=ierror_od)
    if (ierror_od .ne. 0) then
       call log_msg("WARNING: unable to allocate scalar buffer for velocities_od output"&
            &,.true.)
       return
    end if
 allocate(vwd(3,1:nxi,1:nyi,1:nzi),stat=ierror_wd)
    if (ierror_wd .ne. 0) then
       call log_msg("WARNING: unable to allocate scalar buffer for velocities_wd output"&
            &,.true.)
       return
    end if
   allocate(fluxod(3,1:nxi,1:nyi,1:nzi),stat=ierror_fod)
    if (ierror_fod .ne. 0) then
       call log_msg("WARNING: unable to allocate scalar buffer for flux_od output"&
            &,.true.)
       return
    end if
   allocate(fluxwd(3,1:nxi,1:nyi,1:nzi),stat=ierror_fwd)
      if (ierror_fwd .ne. 0) then
         call log_msg("WARNING: unable to allocate scalar buffer for flux_wd output"&
              &,.true.)
         return
      end if




!#ifndef COMMON_VEL_FIX
#ifdef OLD_VEL

  allocate(f_r(nxi,nyi,nzi,3)&
#ifndef SINGLEFLUID
       &,f_b(nxi,nyi,nzi,3)&
#endif
#ifndef NOSURFACTANT
       &,f_s(nxi,nyi,nzi,3)&
#endif
       &,stat=ierror)
  if (ierror/=0) then
     call log_msg("WARNING: unable to allocate force buffer for arrows output"&
          &,.true.)
     return
  end if
! #endifn COMMON_VEL_FIX
#endif



  do z=1,nzi
    do y=1,nyi
      do x=1,nxi
      ! Get Shan Chen forces
      ! only calculate forces on fluid at non-rock sites
        if ( is_fluid(N(x,y,z)%rock_state) ) then

!#ifdef COMMON_VEL_FIX
#ifndef OLD_VEL
           call corrected_velocity_new(N(x,y,z),lbe_force(:,:,x,y,z),corr_vel,corr_flux,x,y,z)
           !call corrected_velocity(N(x,y,z),lbe_force(:,:,x,y,z),corr_vel)
           v(:,x,y,z) = corr_vel(:,0)
           vod(:,x,y,z) = corr_vel(:,1)
           vwd(:,x,y,z) = corr_vel(:,2)
           fluxod(:,x,y,z) = corr_flux(:,1)
           fluxwd(:,x,y,z) = corr_flux(:,2)
            
! #else COMMON_VEL_FIX
#else
#ifdef SINGLEFLUID
          CALL lbe_calculate_sc_forces(N,x,y,z,f_r)
#else
#ifdef NOSURFACTANT
          CALL lbe_calculate_sc_forces(N,x,y,z,f_b,f_r)
#else
          CALL lbe_calculate_sc_forces(N,x,y,z,f_b,f_r,f_s)
#endif
#endif
        else
          f_r(x,y,z,:) = 0.d0
#ifndef SINGLEFLUID
          f_b(x,y,z,:) = 0.d0
#endif
#ifndef NOSURFACTANT
          f_s(x,y,z,:) = 0.d0
#endif
        end if

#ifndef SINGLEFLUID
#ifndef NOSURFACTANT
        v(1,x,y,z)=sum(N(x,y,z)%n_r(:)*g*cx)*amass_r&
          &-(f_r(x,y,z,1)/2.0d0*sum(N(x,y,z)%n_r(:)*g)*amass_r)&
          &+sum(N(x,y,z)%n_b(:)*g*cx)*amass_b&
          &-(f_b(x,y,z,1)/2.0d0*sum(N(x,y,z)%n_b(:)*g)*amass_b)&
          &+sum(N(x,y,z)%n_s(:)*g*cx)*amass_s&
          &-(f_s(x,y,z,1)/2.0d0*sum(N(x,y,z)%n_s(:)*g)*amass_s)
        v(2,x,y,z)=sum(N(x,y,z)%n_r(:)*g*cy)*amass_r&
          &-(f_r(x,y,z,2)/2.0d0*sum(N(x,y,z)%n_r(:)*g)*amass_r)&
          &+sum(N(x,y,z)%n_b(:)*g*cy)*amass_b&
          &-(f_b(x,y,z,2)/2.0d0*sum(N(x,y,z)%n_b(:)*g)*amass_b)&
          &+sum(N(x,y,z)%n_s(:)*g*cy)*amass_s&
          &-(f_s(x,y,z,2)/2.0d0*sum(N(x,y,z)%n_s(:)*g)*amass_s)
        v(3,x,y,z)=sum(N(x,y,z)%n_r(:)*g*cz)*amass_r&
          &-(f_r(x,y,z,3)/2.0d0*sum(N(x,y,z)%n_r(:)*g)*amass_r)&
          &+sum(N(x,y,z)%n_b(:)*g*cz)*amass_b&
          &-(f_b(x,y,z,3)/2.0d0*sum(N(x,y,z)%n_b(:)*g)*amass_b)&
          &+sum(N(x,y,z)%n_s(:)*g*cz)*amass_s&
          &-(f_s(x,y,z,3)/2.0d0*sum(N(x,y,z)%n_s(:)*g)*amass_s)

        v(1,x,y,z)=v(1,x,y,z)/max(10.e-9, real(sum(N(x,y,z)%n_r(:)*g)*amass_r &
                                            & +sum(N(x,y,z)%n_b(:)*g)*amass_b &
                                            & +sum(N(x,y,z)%n_s(:)*g)*amass_s) )
        v(2,x,y,z)=v(2,x,y,z)/max(10.e-9, real(sum(N(x,y,z)%n_r(:)*g)*amass_r &
                                            & +sum(N(x,y,z)%n_b(:)*g)*amass_b &
                                            & +sum(N(x,y,z)%n_s(:)*g)*amass_s) )
        v(3,x,y,z)=v(3,x,y,z)/max(10.e-9, real(sum(N(x,y,z)%n_r(:)*g)*amass_r &
                                            & +sum(N(x,y,z)%n_b(:)*g)*amass_b +sum(N(x,y,z)%n_s(:)*g)*amass_s) )
#else
        v(1,x,y,z)=sum(N(x,y,z)%n_r(:)*g*cx)*amass_r&
          &-(f_r(x,y,z,1)/2.0d0*sum(N(x,y,z)%n_r(:)*g)*amass_r)&
          &+sum(N(x,y,z)%n_b(:)*g*cx)*amass_b&
          &-(f_b(x,y,z,1)/2.0d0*sum(N(x,y,z)%n_b(:)*g)*amass_b)
        v(2,x,y,z)=sum(N(x,y,z)%n_r(:)*g*cy)*amass_r&
          &-(f_r(x,y,z,2)/2.0d0*sum(N(x,y,z)%n_r(:)*g)*amass_r)&
          &+sum(N(x,y,z)%n_b(:)*g*cy)*amass_b&
          &-(f_b(x,y,z,2)/2.0d0*sum(N(x,y,z)%n_b(:)*g)*amass_b)
        v(3,x,y,z)=sum(N(x,y,z)%n_r(:)*g*cz)*amass_r&
          &-(f_r(x,y,z,3)/2.0d0*sum(N(x,y,z)%n_r(:)*g)*amass_r)&
          &+sum(N(x,y,z)%n_b(:)*g*cz)*amass_b&
          &-(f_b(x,y,z,3)/2.0d0*sum(N(x,y,z)%n_b(:)*g)*amass_b)

        v(1,x,y,z)=v(1,x,y,z)/max(10.e-9,real( sum(N(x,y,z)%n_r(:)*g)*amass_r+ &
          &sum(N(x,y,z)%n_b(:)*g)*amass_b))
        v(2,x,y,z)=v(2,x,y,z)/max(10.e-9,real( sum(N(x,y,z)%n_r(:)*g)*amass_r+ &
          &sum(N(x,y,z)%n_b(:)*g)*amass_b))
        v(3,x,y,z)=v(3,x,y,z)/max(10.e-9,real( sum(N(x,y,z)%n_r(:)*g)*amass_r+ &
          &sum(N(x,y,z)%n_b(:)*g)*amass_b))
#endif
#else
        v(1,x,y,z)=sum(N(x,y,z)%n_r(:)*g*cx)*amass_r
        v(2,x,y,z)=sum(N(x,y,z)%n_r(:)*g*cy)*amass_r
        v(3,x,y,z)=sum(N(x,y,z)%n_r(:)*g*cz)*amass_r

        v(1,x,y,z)=(v(1,x,y,z)/max(10.e-9,real( sum(N(x,y,z)%n_r(:)*g)*amass_r)))-f_r(x,y,z,1)/2.0d0
        v(2,x,y,z)=(v(2,x,y,z)/max(10.e-9,real( sum(N(x,y,z)%n_r(:)*g)*amass_r)))-f_r(x,y,z,2)/2.0d0
        v(3,x,y,z)=(v(3,x,y,z)/max(10.e-9,real( sum(N(x,y,z)%n_r(:)*g)*amass_r)))-f_r(x,y,z,3)/2.0d0
#endif

        ! force extra-term

        tx=x+ccoords(1)*nx
        ty=y+ccoords(2)*ny
        tz=z+ccoords(3)*nz
        if( (tz.ge.g_accn_min).and.(tz.le.g_accn_max) .and. &
            (tx.ge.g_accn_min_x).and.(tx.le.g_accn_max_x) .and. &
            (ty.ge.g_accn_min_y).and.(ty.le.g_accn_max_y)) then
          if ( is_fluid(N(x,y,z)%rock_state) ) then
            ! external force term, in the zone where the
            ! external force works
            v(1,x,y,z) = v(1,x,y,z)-g_accn_x/2.0d0
            v(2,x,y,z) = v(2,x,y,z)-g_accn_y/2.0d0
            v(3,x,y,z) = v(3,x,y,z)-g_accn/2.0d0
          end if
        endif

! #endif COMMON_VEL_FIX
#endif
!#ifdef COMMON_VEL_FIX
#ifndef OLD_VEL
        !no rock
     else
        v(:,x,y,z) = 0.0_rk
        vod(:,x,y,z) = 0.0_rk
        vwd(:,x,y,z) = 0.0_rk
        fluxod(:,x,y,z) = 0.0_rk
        fluxwd(:,x,y,z) = 0.0_rk
     endif
#endif
      end do
    end do
  end do
  if ( check_dump_now(sci_arrows, n_sci_arrows) ) then
    call dump_vector(v,'arr')
  endif

  if ( check_dump_now(sci_velocities, n_sci_velocities) ) then
    call dump_scalar(v(1,:,:,:),'velx')
    call dump_scalar(v(2,:,:,:),'vely')
    call dump_scalar(v(3,:,:,:),'velz')
  endif
  if ( check_dump_now(sci_velocities_od, n_sci_velocities_od) ) then
      call dump_scalar(vod(1,:,:,:),'velx_od')
      call dump_scalar(vod(2,:,:,:),'vely_od')
      call dump_scalar(vod(3,:,:,:),'velz_od')
    endif
 if ( check_dump_now(sci_flux_od, n_sci_flux_od) ) then
 call dump_scalar(fluxod(1,:,:,:),'fluxx_od')
 call dump_scalar(fluxod(2,:,:,:),'fluxy_od')
 call dump_scalar(fluxod(3,:,:,:),'fluxz_od')
 endif 
 if ( check_dump_now(sci_velocities_wd, n_sci_velocities_wd) ) then
        call dump_scalar(vwd(1,:,:,:),'velx_wd')
        call dump_scalar(vwd(2,:,:,:),'vely_wd')
        call dump_scalar(vwd(3,:,:,:),'velz_wd')
       !call dump_scalar(fluxwd(3,:,:,:),'fluxz_wd')
      endif
if ( check_dump_now(sci_flux_wd, n_sci_flux_wd) ) then
   call dump_scalar(fluxwd(1,:,:,:),'fluxx_wd')
   call dump_scalar(fluxwd(2,:,:,:),'fluxy_wd')
   call dump_scalar(fluxwd(3,:,:,:),'fluxz_wd')
   endif


  deallocate(v)
  deallocate(vod)
  deallocate(vwd)
  deallocate(fluxod)
  deallocate(fluxwd)

!#ifndef COMMON_VEL_FIX
#ifdef OLD_VEL
  deallocate(f_r&
#ifndef SINGLEFLUID
       &,f_b&
#endif
#ifndef NOSURFACTANT
       &,f_s&
#endif
       &,stat=ierror)
#endif

end subroutine dump_arrows

!>  dumps rock_state of whole system
subroutine dump_rock(N)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N
  real*8, dimension(:,:,:), allocatable :: scalar
  integer :: ierror,x,y,z,nxi,nyi,nzi

  nxi = size(N,1)-2
  nyi = size(N,2)-2
  nzi = size(N,3)-2

  allocate(scalar(1:nxi,1:nyi,1:nzi),stat=ierror)

  if (ierror .ne. 0) then
    call log_msg("WARNING: unable to allocate scalar buffer for rock output",.true.)
    return	! Inability to write output is not quite fatal.
  end if

  do z=1,nzi
     do y=1,nyi
        do x=1,nxi
           scalar(x,y,z) = N(x,y,z)%rock_state
        end do
     end do
  end do

  call dump_scalar(scalar,'rock_state')
  deallocate(scalar)
end subroutine dump_rock


!>  dumps rock_colour of whole system
subroutine dump_rock_colour(N)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N
  real*8, dimension(:,:,:), allocatable :: scalar
  integer :: ierror,x,y,z,nxi,nyi,nzi

  nxi = size(N,1)-2
  nyi = size(N,2)-2
  nzi = size(N,3)-2

  allocate(scalar(1:nxi,1:nyi,1:nzi),stat=ierror)

  if (ierror .ne. 0) then
    call log_msg("WARNING: unable to allocate scalar buffer for rock colour output",.true.)
    return	! Inability to write output is not quite fatal.
  end if

  do z=1,nzi
     do y=1,nyi
        do x=1,nxi
           scalar(x,y,z) = N(x,y,z)%rock_colour 
        end do
     end do
  end do

  call dump_scalar(scalar,'rock_colour')
  deallocate(scalar)
end subroutine dump_rock_colour

!>  dumps rock_rho_r of whole system
subroutine dump_rock_rho_r(N)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N
  real*8, dimension(:,:,:), allocatable :: scalar
  integer :: ierror,x,y,z,nxi,nyi,nzi
  real*8 :: rockfloat

  nxi = size(N,1)-2
  nyi = size(N,2)-2
  nzi = size(N,3)-2

  allocate(scalar(1:nxi,1:nyi,1:nzi),stat=ierror)

  if (ierror .ne. 0) then
    call log_msg("WARNING: unable to allocate scalar buffer for rock rho r output",.true.)
    return	! Inability to write output is not quite fatal.
  end if

  do z=1,nzi
     do y=1,nyi
        do x=1,nxi
           if ( is_wall(N(x,y,z)%rock_state) ) then
              scalar(x,y,z) = N(x,y,z)%n_r(restvec)
           else 
              scalar(x,y,z) = 0.0
           endif
        end do
     end do
  end do

  call dump_scalar(scalar,'rock_rho_r')
  deallocate(scalar)
end subroutine dump_rock_rho_r

#ifndef SINGLEFLUID
!>  dumps rock_rho_b of whole system
subroutine dump_rock_rho_b(N)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N
  real*8, dimension(:,:,:), allocatable :: scalar
  integer :: ierror,x,y,z,nxi,nyi,nzi

  nxi = size(N,1)-2
  nyi = size(N,2)-2
  nzi = size(N,3)-2

  allocate(scalar(1:nxi,1:nyi,1:nzi),stat=ierror)

  if (ierror .ne. 0) then
    call log_msg("WARNING: unable to allocate scalar buffer for rock rho b output",.true.)
    return	! Inability to write output is not quite fatal.
  end if

  do z=1,nzi
     do y=1,nyi
        do x=1,nxi
           if ( is_wall(N(x,y,z)%rock_state) ) then
              scalar(x,y,z) = N(x,y,z)%n_b(restvec)
           else 
              scalar(x,y,z) = 0.0
           endif
        end do
     end do
  end do

  call dump_scalar(scalar,'rock_rho_b')
  deallocate(scalar)
end subroutine dump_rock_rho_b
#endif

!>Calculates total fluid velocity at each site, then calls \c dump_vector.
!>
!>Based on the kinetic theory definition for gas mixtures
!>(cf. Chapman & Cowling, ``The mathematical theory of non-uniform
!>gases." (CUP: Cambridge, 1970, 3rd ed.)), which for a discrete
!>velocity set numbered by \f$k\f$ and species masses \f$m^{\alpha}\f$, reads:
!>\f[\mathbf{u}=\frac{\rho\mathbf{u}}{\rho}
!>\equiv
!>\frac{\sum_{\alpha}m^{\alpha}\sum_k\mathbf{c}_k n_k^{\alpha}}
!>{\sum_{\alpha}m^{\alpha}\sum_k n_k^{\alpha}}\f]
!>where \f$n_k^{alpha}\f$ is the number density, of the real var
!>\c N(x,y,z)%n_\f$\alpha\f$\c (k), \f$u\f$=total fluid velocity (the unknown),
!>\f$f\f$=single-particle d.f.
!>
!>Added by N. Gonzalez-Segredo 26.09.03.
subroutine tot_vel(N,u,switch,verb)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N
  real*8, dimension(:,:,:,:) :: u
  integer :: ierror,x,y,z,nxi,nyi,nzi
  real*8  :: dens
  real*8, parameter :: eps = 1.e-20
  logical :: switch, verb
  logical :: rep = .false.

  nxi = size(N,1)-2
  nyi = size(N,2)-2
  nzi = size(N,3)-2

  do z=1,nzi
    do y=1,nyi
      do x=1,nxi
        if (.not. is_fluid(N(x,y,z)%rock_state) ) then
          u(:,x,y,z) = 0.0_8
        else
#ifndef NOSURFACTANT
          dens = amass_r*sum(N(x,y,z)%n_r(:)*g) + &
                 amass_b*sum(N(x,y,z)%n_b(:)*g) + &
                 amass_s*sum(N(x,y,z)%n_s(:)*g)
          if( (verb .eqv. .true.).and.(dens < eps).and.(rep .eqv. .false.) ) then
            write(msgstr,"('WARNING lbe_io.F90: tot_vel(): dens < ',F16.10)") eps
            call log_msg(msgstr)
            rep = .true.
          endif
          u(1,x,y,z)= amass_r*sum(N(x,y,z)%n_r(:)*g*cx) + &
                      amass_b*sum(N(x,y,z)%n_b(:)*g*cx) + &
                      amass_s*sum(N(x,y,z)%n_s(:)*g*cx)
          u(1,x,y,z)= u(1,x,y,z)/dens
          u(2,x,y,z)= amass_r*sum(N(x,y,z)%n_r(:)*g*cy) + &
                      amass_b*sum(N(x,y,z)%n_b(:)*g*cy) + &
                      amass_s*sum(N(x,y,z)%n_s(:)*g*cy)
          u(2,x,y,z)= u(2,x,y,z)/dens
          u(3,x,y,z)= amass_r*sum(N(x,y,z)%n_r(:)*g*cz) + &
                      amass_b*sum(N(x,y,z)%n_b(:)*g*cz) + &
                      amass_s*sum(N(x,y,z)%n_s(:)*g*cz)
          u(3,x,y,z)=	u(3,x,y,z)/dens
#else
#ifndef SINGLEFLUID
          dens = amass_r*sum(N(x,y,z)%n_r(:)*g) + &
                 amass_b*sum(N(x,y,z)%n_b(:)*g)
          if( (verb .eqv. .true.).and.(dens < eps).and.(rep .eqv. .false.) ) then
            write(msgstr,"('WARNING lbe_io.F90: tot_vel(): dens < ',F16.10)") eps
            call log_msg(msgstr)
            rep = .true.
          endif

          u(1,x,y,z)= amass_r*sum(N(x,y,z)%n_r(:)*g*cx) + &
                      amass_b*sum(N(x,y,z)%n_b(:)*g*cx)
          u(2,x,y,z)= amass_r*sum(N(x,y,z)%n_r(:)*g*cy) + &
                      amass_b*sum(N(x,y,z)%n_b(:)*g*cy)
          u(3,x,y,z)= amass_r*sum(N(x,y,z)%n_r(:)*g*cz) + &
                      amass_b*sum(N(x,y,z)%n_b(:)*g*cz)
#else
          dens = amass_r*sum(N(x,y,z)%n_r(:)*g)
          if( (verb .eqv. .true.).and.(dens < eps).and.(rep .eqv. .false.) ) then
            write(msgstr,"('WARNING lbe_io.F90: tot_vel(): dens < ',F16.10)") eps
            call log_msg(msgstr)
            rep = .true.
          endif
          u(1,x,y,z) = amass_r*sum(N(x,y,z)%n_r(:)*g*cx)
          u(2,x,y,z) = amass_r*sum(N(x,y,z)%n_r(:)*g*cy)
          u(3,x,y,z) = amass_r*sum(N(x,y,z)%n_r(:)*g*cz)
#endif
          if (dens.lt.eps) then
            u(1,x,y,z) = 0.
            u(2,x,y,z) = 0.
            u(3,x,y,z) = 0.
          else
            u(1,x,y,z) = u(1,x,y,z)/dens
            u(2,x,y,z) = u(2,x,y,z)/dens
            u(3,x,y,z) = u(3,x,y,z)/dens
          endif
#endif
        end if
      end do
    end do
  end do
  if(switch .eqv. .true.) then
    CALL dump_vector(u,'tvel')
  endif
end subroutine tot_vel
!> \}

!> \{
!> \name Array-specific stubroutines - Format-specific routines
!>
!> For an array of type \c foo, there will be a routine called \c
!> dump_foo(), which will in turn call \c dump_foo_all(), \c
!> dump_foo_bin(), or some other routine \c dump_foo_bar(), for file
!> format \c bar.

!> Dumps a floating point scalar field in a format depending on \c
!> dump_format.
!>
!> \param[in,out] scalar local portion of the scalar field
!>
!> \param[in] name string to be used to generate the filename from
subroutine dump_scalar(scalar,name)
    real*8, dimension(1:,1:,1:),intent(inout) :: scalar
    character(len=*)            :: name

    select case(trim(dump_format))
    case ('hdf')
#ifdef USEHDF
       call dump_scalar_phdf5(scalar, name)
#else
       call log_msg("HDF5 support switched off")
#endif
    case ('mpi')
#ifdef USEXDRF
       call dump_scalar_xdr_parallel(scalar, name)
#else
       call log_msg("XDRF support switched off")
#endif
    case ('all')
       call dump_scalar_all(scalar, name)
    case ('xdr')
#ifdef USEXDRF
       call dump_scalar_xdr(scalar, name)
#else
       call log_msg("XDRF support switched off")
#endif
    case ('vtk')
#ifdef USEXDRF
       call dump_scalar_vtk(scalar, name)
#else
       call log_msg("XDRF/VTK support switched off")
#endif
    case ('bin')
       call dump_scalar_bin(scalar, name)
    case default
       call error('unknown value: dump_format="'//trim(dump_format)//'"')
    end select
end subroutine dump_scalar

!> Dumps an integer scalar field in a format depending on \c
!> dump_format.
!>
!> \param[in,out] iscalar local portion of the scalar field
!>
!> \param[in] name string to be used to generate the filename from
subroutine dump_iscalar(iscalar,name)
    integer,dimension(1:,1:,1:),intent(inout) :: iscalar
    character(len=*)            :: name

    select case (trim(dump_format))
    case ('hdf')
#ifdef USEHDF
       call dump_iscalar_phdf5(iscalar,name)
#else
       call log_msg("HDF5 support switched off")
#endif
    case ('mpi')
#ifdef USEXDRF
       call error('dump_iscalar_xdr_parallel() not implemented yet')
!!$          call dump_iscalar_xdr_parallel(iscalar, name)
#else
       call log_msg("XDRF support switched off")
#endif
    case ('all')
       call error('dump_iscalar_all() not implemented yet')
!!$          call dump_iscalar_all(iscalar, name)
    case ('xdr')
#ifdef USEXDRF
       call dump_iscalar_xdr(iscalar,name)
#else
       call log_msg("XDRF support switched off")
#endif
    case ('vtk')
#ifdef USEXDRF
       call error('dump_iscalar_vtk() not implemented yet')
!!$          call dump_iscalar_vtk(iscalar, name)
#else
       call log_msg("XDRF/VTK support switched off")
#endif
    case ('bin')
       call error('dump_iscalar_bin() not implemented yet')
!!$          call dump_iscalar_bin(iscalar, name)
    case default
       call error('unknown value: dump_format="'//trim(dump_format)//'"')
    end select
end subroutine dump_iscalar

!> Dumps a scalar value in ASCII "all" form.
subroutine dump_scalar_all(scalar,name)
  implicit none
  real*8, dimension(1:,1:,1:),intent(inout) :: scalar
  character(len=*), intent(in) :: name
  character(len=1024) :: filename
  integer :: x,y,z,nxi,nyi,nzi

  nxi = size(scalar,1)
  nyi = size(scalar,2)
  nzi = size(scalar,3)

  call lbe_make_filename_output(filename,trim(name),'.all',nt)
  
  open(10,file=filename)

  do z=1,nzi
   do y=1,nyi
    do x=1,nxi
             if (dump_double) then
    write(10,'(3i3,f12.8)')         &
      ccoords(1)*nxi+x,       &
      ccoords(2)*nyi+y,       &
      ccoords(3)*nzi+z,       &
        scalar(x,y,z)
             else
    write(10,'(3i3,f8.4)')          &
      ccoords(1)*nxi+x,       &
      ccoords(2)*nyi+y,       &
      ccoords(3)*nzi+z,       &
        scalar(x,y,z)
             end if
    end do
   end do
  end do

  close(10)
  !write(*,'(a,a,a)') 'Wrote <',trim(filename),'>'

        if (write_AVS_fld) call dump_avs_fld(trim(name),nxi,nyi,nzi,1)
end subroutine dump_scalar_all

!> Dumps a scalar in binary format.
subroutine dump_scalar_bin(scalar,name)
  implicit none
  real*8, dimension(1:,1:,1:),intent(inout) :: scalar
  character(len=*), intent(in) :: name
  character(len=1024) :: filename
  integer :: x,y,z,nxi,nyi,nzi,ierror
  real*4, dimension(:,:,:),allocatable :: scalar2

  nxi = size(scalar,1)
  nyi = size(scalar,2)
  nzi = size(scalar,3)
 
  call lbe_make_filename_output(filename,trim(name),'.bin',nt)
  
  open(10,file=filename,form="unformatted")
  if (dump_double) then
    write(10) scalar(:,:,:)
  else

    ! FIXME: Very unelegegant way of doing this, 
    ! but conversion functions did not work
    allocate(scalar2(nxi,nyi,nzi),stat=ierror)
    scalar2 = scalar
    write(10) scalar2(:,:,:)
    deallocate(scalar2)
  end if
  close(10)
        
  !write(*,'(a,a,a)') 'Wrote <',trim(filename),'>'

        if (write_AVS_fld) call dump_avs_fld(trim(name),nxi,nyi,nzi,1)
end subroutine dump_scalar_bin

!> Calls \c dump_vector_all or \c dump_vector_bin depending on the value
!> of the variable \c dump_format.
subroutine dump_vector(vector,name)
    real*8, dimension(1:,1:,1:,1:),intent(inout) :: vector
    character(len=*), intent(in) :: name
    character(len=1024) :: filename

    if (index(dump_format,'all').gt.0) then
       call dump_vector_all(vector,name)

    elseif(index(dump_format,'hdf').gt.0)then
#ifdef USEHDF
       call lbe_make_filename_output(filename, name, '.h5', nt)
       call dump_vector_phdf5(vector,filename)
#else
       call log_msg("HDF5 support switched off")
#endif

    elseif (index(dump_format,'mpi').gt.0) then
#ifdef USEXDRF
       call dump_vector_xdr_parallel(vector,name)
#else
       call log_msg("XDRF support switched off")
#endif

    elseif (index(dump_format,'xdr').gt.0) then
#ifdef USEXDRF
       call dump_vector_xdr(vector,name)
#else
       call log_msg("XDRF support switched off")
#endif

    else
       call dump_vector_bin(vector,name)
    endif
end subroutine dump_vector

!> Dumps a vector field in ASCII "all" form.
subroutine dump_vector_all(vector,name)
	implicit none
	real*8, dimension(1:,1:,1:,1:),intent(inout) :: vector
	character(len=*), intent(in) :: name
	character(len=1024) :: filename
	integer :: x,y,z,nxi,nyi,nzi

        nxi = size(vector,2)
        nyi = size(vector,3)
        nzi = size(vector,4)

	call lbe_make_filename_output(filename,trim(name),'.all',nt)
	open(10,file=filename)

	do z=1,nzi
	 do y=1,nyi
	  do x=1,nxi
             if (dump_double) then
		write(10,'(3i3,3f12.8)')				&
			ccoords(1)*nxi+x,				&
			ccoords(2)*nyi+y,				&
			ccoords(3)*nzi+z,				&
	  		vector(:,x,y,z)
             else
		write(10,'(3i3,3f8.4)')				&
			ccoords(1)*nxi+x,				&
			ccoords(2)*nyi+y,				&
			ccoords(3)*nzi+z,				&
	  		vector(:,x,y,z)
             end if
	  end do
	 end do
	end do

	close(10)
	!write(*,'(a,a,a)') 'Wrote <',trim(filename),'>'
        if (write_AVS_fld) call dump_avs_fld(trim(name),nxi,nyi,nzi,int(4))
end subroutine dump_vector_all

!> Dump an array of vectors in binary format.
subroutine dump_vector_bin(vector,name)
	implicit none
	real*8, dimension(1:,1:,1:,1:),intent(inout) :: vector
        real*4, dimension(:,:,:,:),allocatable :: vector2 
	character(len=*), intent(in) :: name
	character(len=1024) :: filename
	integer :: x,y,z,nxi,nyi,nzi,ierror

        nxi = size(vector,2)
        nyi = size(vector,3)
        nzi = size(vector,4)

	call lbe_make_filename_output(filename,trim(name),'.bin',nt)
	
	open(10,file=filename,form="unformatted")

          if (dump_double) then
	     write(10) vector(:,:,:,:)
          else
             allocate(vector2(size(vector,1),nxi,nyi,nzi),stat=ierror)
             vector2 = vector
	     write(10) vector2(:,:,:,:)
             deallocate(vector2)
          end if

	close(10)
	!write(*,'(a,a,a)') 'Wrote <',trim(filename),'>'
        if (write_AVS_fld) call dump_avs_fld(trim(name),nxi,nyi,nzi,int(4))
end subroutine dump_vector_bin
!> \}

!> \{
!> \name lbe_pressure routines
!>
!>T hese routines:
!>     [1] compute the components of the pressure tensor
!>     [2] compute Pzz-Pxx along a line perpendicular to
!>     a planar interface. (Not used.)
!>
!>Added by NGS. Modified by NGS on 26.09.03.

!> Calculates all the components of the pressure tensor.
!> 's' decides whether or not to include scalar pressure calculation.
subroutine dump_pressure(N,s)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N
  real*8, dimension(:,:,:), allocatable :: pxx, pyy, pzz,    &
                                           pxy, pyz, pxz,    &
                                           scalarpressure
  real*8, dimension(3,3) :: sum_psixc_r
  real*8, parameter :: omega = 0.25
  integer :: ierror,x,y,z,s,xc,yc,zc,i,j,nxi,nyi,nzi
  real*8 :: psixc_r
  real*8 :: psix_r
  real*8, dimension(:,:,:,:), allocatable :: u !Modif. NGS
#ifndef SINGLEFLUID
  real*8, dimension(3,3) :: sum_psixc_b, p_br
  real*8 :: psixc_b, psix_b
#ifndef NOSURFACTANT
  real*8, dimension(3,3) :: p_bs, p_ss, p_rs, sum_psixc_s
  real*8 :: psixc_s, psix_s
#endif
#endif

  ! !DEBUG
  !   integer, dimension(5) :: sample_point
  !   sample_point(1)=size(N,1)/8
  !   sample_point(2)=size(N,1)/4
  !   sample_point(3)=size(N,1)/2
  !   sample_point(4)=size(N,1)*(1.0-1.0/4.0)
  !   sample_point(5)=size(N,1)*(1.0-1.0/8.0)
  ! !END_DEBUG

  nxi = size(N,1)-2
  nyi = size(N,2)-2
  nzi = size(N,3)-2

  allocate(u(3,1:nxi,1:nyi,1:nzi),stat=ierror) !Modif. NGS
  if (ierror .ne. 0) then
    call log_msg("WARNING: unable to allocate total velocity real array u(:,:,:,:)",.true.)
    return  ! Inability to write output is not quite fatal.
  end if
  allocate(pxx(1:nxi,1:nyi,1:nzi),stat=ierror)
  if (ierror .ne. 0) then
    call log_msg("WARNING: unable to allocate scalar buffer p_xx(:,:,:)",.true.)
    return  ! Inability to write output is not quite fatal.
  end if
  allocate(pyy(1:nxi,1:nyi,1:nzi),stat=ierror)
  if (ierror .ne. 0) then
    call log_msg("WARNING: unable to allocate scalar buffer p_yy(:,:,:)",.true.)
    return  ! Inability to write output is not quite fatal.
  end if
  allocate(pzz(1:nxi,1:nyi,1:nzi),stat=ierror)
  if (ierror .ne. 0) then
    call log_msg("WARNING: unable to allocate scalar buffer p_zz(:,:,:)",.true.)
    return  ! Inability to write output is not quite fatal.
  end if
  allocate(pxy(1:nxi,1:nyi,1:nzi),stat=ierror)
  if (ierror .ne. 0) then
    call log_msg("WARNING: unable to allocate scalar buffer p_xy(:,:,:)",.true.)
    return  ! Inability to write output is not quite fatal.
  end if
  allocate(pyz(1:nxi,1:nyi,1:nzi),stat=ierror)
  if (ierror .ne. 0) then
    call log_msg("WARNING: unable to allocate scalar buffer p_yz(:,:,:)",.true.)
    return  ! Inability to write output is not quite fatal.
  end if
  allocate(pxz(1:nxi,1:nyi,1:nzi),stat=ierror)
  if (ierror .ne. 0) then
    call log_msg("WARNING: unable to allocate scalar buffer p_xz(:,:,:)",.true.)
    return  ! Inability to write output is not quite fatal.
  end if
  if(s == droplet) then
    allocate(scalarpressure(1:nxi,1:nyi,1:nzi),stat=ierror)
    if (ierror .ne. 0) then
      call log_msg("WARNING: unable to allocate scalar buffer scalarpressure(:,:,:)",.true.)
      return  ! Inability to write output is not quite fatal.
    end if
  endif

  ! Health measure:
  pxx = 0.
  pyy = 0.
  pzz = 0.
  pxy = 0.
  pyz = 0.
  pxz = 0.
  if (s == droplet) then
    scalarpressure = 0.
  end if

  call tot_vel(N,u,.false.,.false.)
  do z=1,nzi
    do y=1,nyi
      do x=1,nxi

        ! INTER-COMPONENT-INDEPENDENT PART
        ! OF THE PRESSURE TENSOR.
        ! DEFINED IN TERMS OF THE BARICENTRIC
        ! MICROSCOPIC VELOCITY    --NGS
        pxx(x,y,z) = amass_r*&
            &sum(N(x,y,z)%n_r(:)*g*(cx-u(1,x,y,z))*(cx-u(1,x,y,z)))
        pyy(x,y,z) = amass_r*&
            &sum(N(x,y,z)%n_r(:)*g*(cy-u(2,x,y,z))*(cy-u(2,x,y,z)))
        pzz(x,y,z) = amass_r*&
            &sum(N(x,y,z)%n_r(:)*g*(cz-u(3,x,y,z))*(cz-u(3,x,y,z)))
        pxy(x,y,z) = amass_r*&
            &sum(N(x,y,z)%n_r(:)*g*(cx-u(1,x,y,z))*(cy-u(2,x,y,z)))
        pxz(x,y,z) = amass_r*&
            &sum(N(x,y,z)%n_r(:)*g*(cx-u(1,x,y,z))*(cz-u(3,x,y,z)))
        pyz(x,y,z) = amass_r*&
            &sum(N(x,y,z)%n_r(:)*g*(cy-u(2,x,y,z))*(cz-u(3,x,y,z)))
#ifndef SINGLEFLUID
        pxx(x,y,z) = pxx(x,y,z) + amass_b*&
            &sum(N(x,y,z)%n_b(:)*g*(cx-u(1,x,y,z))*(cx-u(1,x,y,z)))
        pyy(x,y,z) = pyy(x,y,z) + amass_b*&
            &sum(N(x,y,z)%n_b(:)*g*(cy-u(2,x,y,z))*(cy-u(2,x,y,z)))
        pzz(x,y,z) = pzz(x,y,z) + amass_b*&
            &sum(N(x,y,z)%n_b(:)*g*(cz-u(3,x,y,z))*(cz-u(3,x,y,z)))
        pxy(x,y,z) = pxy(x,y,z) + amass_b*&
            &sum(N(x,y,z)%n_b(:)*g*(cx-u(1,x,y,z))*(cy-u(2,x,y,z)))
        pxz(x,y,z) = pxz(x,y,z) + amass_b*&
            &sum(N(x,y,z)%n_b(:)*g*(cx-u(1,x,y,z))*(cz-u(3,x,y,z)))
        pyz(x,y,z) = pyz(x,y,z) + amass_b*&
            &sum(N(x,y,z)%n_b(:)*g*(cy-u(2,x,y,z))*(cz-u(3,x,y,z)))
#ifndef NOSURFACTANT
        pxx(x,y,z) = pxx(x,y,z) + amass_s*&
            &sum(N(x,y,z)%n_s(:)*g*(cx-u(1,x,y,z))*(cx-u(1,x,y,z)))
        pyy(x,y,z) = pyy(x,y,z) + amass_s*&
            &sum(N(x,y,z)%n_s(:)*g*(cy-u(2,x,y,z))*(cy-u(2,x,y,z)))
        pzz(x,y,z) = pzz(x,y,z) + amass_s*&
            &sum(N(x,y,z)%n_s(:)*g*(cz-u(3,x,y,z))*(cz-u(3,x,y,z)))
        pxy(x,y,z) = pxy(x,y,z) + amass_s*&
            &sum(N(x,y,z)%n_s(:)*g*(cx-u(1,x,y,z))*(cy-u(2,x,y,z)))
        pxz(x,y,z) = pxz(x,y,z) + amass_s*&
            &sum(N(x,y,z)%n_s(:)*g*(cx-u(1,x,y,z))*(cz-u(3,x,y,z)))
        pyz(x,y,z) = pyz(x,y,z) + amass_s*&
            &sum(N(x,y,z)%n_s(:)*g*(cy-u(2,x,y,z))*(cz-u(3,x,y,z)))
#endif
#endif

        !DEBUG
        !do i=1,5
        !   if(x==sample_point(i) .and.				  &
        !   y==nyi/2 .and. z==nzi/2) then
        !      print*,'At (',x,',',y,',',z,'):'
        !      print*,'  pxz_kinetic = ', pxz(x,y,z)
        !      print*,'  u_x =', u(1,x,y,z)
        !      print*,'  u_z =', u(3,x,y,z)
        !      print*,'  amass_r * sum(g * (cx-ux) * (cz-uz) * n) =',&
        !      amass_r * sum(N(x,y,z)%n_r(:) *			  &
        !      g * (cx-u(1,x,y,z)) * (cz-u(3,x,y,z)))
        !      print*,'  amass_b * sum(g * (cx-ux) * (cz-uz) * n) =',&
        !      amass_b * sum(N(x,y,z)%n_b(:) *			  &
        !      g * (cx-u(1,x,y,z)) * (cz-u(3,x,y,z)))
        !      print*,'  amass_s * sum(g * (cx-ux) * (cz-uz) * n) =',&
        !      amass_s * sum(N(x,y,z)%n_s(:) *			  &
        !      g * (cx-u(1,x,y,z)) * (cz-u(3,x,y,z)))
        !      print*,'0 =? sum_alpha(mass * sum((c_k - u) * n)) = (',&
        !      amass_r * sum(N(x,y,z)%n_r(:) * g * (cx-u(1,x,y,z))) +&
        !      amass_b * sum(N(x,y,z)%n_b(:) * g * (cx-u(1,x,y,z))) +&
        !      amass_s * sum(N(x,y,z)%n_s(:) * g * (cx-u(1,x,y,z))), &
        !      ',',						    &
        !      amass_r * sum(N(x,y,z)%n_r(:) * g * (cy-u(2,x,y,z))) +&
        !      amass_b * sum(N(x,y,z)%n_b(:) * g * (cy-u(2,x,y,z))) +&
        !      amass_s * sum(N(x,y,z)%n_s(:) * g * (cy-u(2,x,y,z))), &
        !      ',',						    &
        !      amass_r * sum(N(x,y,z)%n_r(:) * g * (cz-u(3,x,y,z))) +&
        !      amass_b * sum(N(x,y,z)%n_b(:) * g * (cz-u(3,x,y,z))) +&
        !      amass_s * sum(N(x,y,z)%n_s(:) * g * (cz-u(3,x,y,z))), &
        !      ')'
        !   endif
        !enddo
        !END

        ! INTER-COMPONENT-DEPENDENT PART
        ! OF THE PRESSURE TENSOR
        ! DOES NOT INCLUDE VELOCITY SINCE IT'S NOT A 
        ! KINETIC CONTRIBUTION    --NGS

        ! Set the contributions to zero initially
        !
        ! They are array(3x3):
        sum_psixc_r = 0.
        psix_r = sum( N(x,y,z)%n_r(:)*g(:) )
#ifndef SINGLEFLUID
        sum_psixc_b = 0.
        psix_b = sum( N(x,y,z)%n_b(:)*g(:) )
#ifndef NOSURFACTANT
        sum_psixc_s = 0.
        psix_s = sum( N(x,y,z)%n_s(:)*g(:) )
#endif
#endif

        select case(psifunc)
          case (0)
            ! Original code: psi=number density
            ! plus funny clipping routine
            psix_r = min(1., real(psix_r))
#ifndef SINGLEFLUID
            psix_b = min(1., real(psix_b))
#ifndef NOSURFACTANT
            psix_s = min(1., real(psix_s))
#endif
#endif
          case (1)
            ! psi = n, with no clipping
            ! nothing further to do.
          case (2)
            ! psi = 1 - exp(-n)
            ! no clipping
            psix_r = 1 - exp(-psix_r)
#ifndef SINGLEFLUID
            psix_b = 1 - exp(-psix_b)
#ifndef NOSURFACTANT
            psix_s = 1 - exp(-psix_s)
#endif
#endif
          case(3)
  ! Carnahan-Starling EOS
   ! psi =
   ! sqrt(((2*psi*Tcs*(1+psi+psi**2-psi**3)/(1-psi)**3)-4/3*psi**2)/(6*g_br))
            psix_r = sqrt(2*((psix_r*Tcs*(1+psix_r+psix_r**2-psix_r**3)/(1-psix_r)**3)-4/3*psix_r**2)/(6*g_br))
#ifndef SINGLEFLUID
            psix_b = sqrt(2*((psix_b*Tcs*(1+psix_b+psix_b**2-psix_b**3)/(1-psix_b)**3)-4/3*psix_b**2)/(6*g_br))
#ifndef NOSURFACTANT
            psix_s = sqrt(2*((psix_s*Tcs*(1+psix_s+psix_s**2-psix_s**3)/(1-psix_s)**3)-4/3*psix_s**2)/(6*g_br))
#endif
#endif

          case default
            call log_msg("ERROR: Unknown psi functional, aborting...")
            call Abend
        end select

        do i=1,nnonrest
          xc = x + cx(i)
          yc = y + cy(i)
          zc = z + cz(i)
          ! psi_b = no density of blue particles at the site (ixa,iya,iza),
          ! calculated by summing all those advecting from that site.

          psixc_r = sum( N(xc,yc,zc)%n_r(:)*g(:) )
#ifndef SINGLEFLUID
          psixc_b = sum( N(xc,yc,zc)%n_b(:)*g(:) )
#ifndef NOSURFACTANT
          psixc_s = sum( N(xc,yc,zc)%n_s(:)*g(:) )
#endif
#endif

          select case(psifunc)
            case (0)
              ! Original code: psi=number density
              ! plus funny clipping routine
              psixc_r = min(1., real(psixc_r))
#ifndef SINGLEFLUID
              psixc_b = min(1., real(psixc_b))
#ifndef NOSURFACTANT
              psixc_s = min(1., real(psixc_s))
#endif
#endif
            case (1)
              ! psi = n, with no clipping
              ! nothing further to do.
            case (2)
              ! psi = 1 - exp(-n)
              ! no clipping
              psixc_r = 1 - exp(-psixc_r)
#ifndef SINGLEFLUID
              psixc_b = 1 - exp(-psixc_b)
#ifndef NOSURFACTANT
              psixc_s = 1 - exp(-psixc_s)
#endif
#endif

 				case(3)
    ! Carnahan-Starling EOS
     ! psi = sqrt(((2*psi*Tcs*(1+psi+psi**2-psi**3)/(1-psi)**3)-4/3*psi**2)/(6*g_br))
              psixc_r = sqrt(2*((psixc_r*Tcs*(1+psixc_r+psixc_r**2-psixc_r**3)/(1-psixc_r)**3)-4/3*psixc_r**2)/(6*g_br))
#ifndef SINGLEFLUID
              psixc_b = sqrt(2*((psixc_b*Tcs*(1+psixc_b+psixc_b**2-psixc_b**3)/(1-psixc_b)**3)-4/3*psixc_b**2)/(6*g_br))
#ifndef NOSURFACTANT
              psixc_s = sqrt(2*((psixc_s*Tcs*(1+psixc_s+psixc_s**2-psixc_s**3)/(1-psixc_s)**3)-4/3*psixc_s**2)/(6*g_br))
#endif
#endif
            case default
              call log_msg("ERROR: Unknown psi functional, aborting...")
              call Abend
          end select

          sum_psixc_r(1,1) = sum_psixc_r(1,1) + psixc_r * g(i) * cx(i) * cx(i)
          sum_psixc_r(2,2) = sum_psixc_r(2,2) + psixc_r * g(i) * cy(i) * cy(i)
          sum_psixc_r(3,3) = sum_psixc_r(3,3) + psixc_r * g(i) * cz(i) * cz(i)
          sum_psixc_r(1,2) = sum_psixc_r(1,2) + psixc_r * g(i) * cx(i) * cy(i)
          sum_psixc_r(1,3) = sum_psixc_r(1,3) + psixc_r * g(i) * cx(i) * cz(i)
          sum_psixc_r(2,3) = sum_psixc_r(2,3) + psixc_r * g(i) * cy(i) * cz(i)
#ifndef SINGLEFLUID
          sum_psixc_b(1,1) = sum_psixc_b(1,1) + psixc_b * g(i) * cx(i) * cx(i)
          sum_psixc_b(2,2) = sum_psixc_b(2,2) + psixc_b * g(i) * cy(i) * cy(i)
          sum_psixc_b(3,3) = sum_psixc_b(3,3) + psixc_b * g(i) * cz(i) * cz(i)
          sum_psixc_b(1,2) = sum_psixc_b(1,2) + psixc_b * g(i) * cx(i) * cy(i)
          sum_psixc_b(1,3) = sum_psixc_b(1,3) + psixc_b * g(i) * cx(i) * cz(i)
          sum_psixc_b(2,3) = sum_psixc_b(2,3) + psixc_b * g(i) * cy(i) * cz(i)
#ifndef NOSURFACTANT
          sum_psixc_s(1,1) = sum_psixc_s(1,1) + psixc_s * g(i) * cx(i) * cx(i)
          sum_psixc_s(2,2) = sum_psixc_s(2,2) + psixc_s * g(i) * cy(i) * cy(i)
          sum_psixc_s(3,3) = sum_psixc_s(3,3) + psixc_s * g(i) * cz(i) * cz(i)
          sum_psixc_s(1,2) = sum_psixc_s(1,2) + psixc_s * g(i) * cx(i) * cy(i)
          sum_psixc_s(1,3) = sum_psixc_s(1,3) + psixc_s * g(i) * cx(i) * cz(i)
          sum_psixc_s(2,3) = sum_psixc_s(2,3) + psixc_s * g(i) * cy(i) * cz(i)
#endif
#endif
        end do !Next velocity i

        ! PRESSURE TENSOR IS DECOMPOSED INTO INTERACTIONS:
#ifndef SINGLEFLUID
        p_br = g_br * (psix_r * sum_psixc_b + psix_b * sum_psixc_r)
#ifndef NOSURFACTANT
        p_bs = g_bs * (psix_b * sum_psixc_s + psix_s * sum_psixc_b)
        p_ss = g_ss * (psix_s * sum_psixc_s) * 2.0
        p_rs = g_bs * (psix_r * sum_psixc_s + psix_s * sum_psixc_r)
#endif
#endif
        ! g_rs = g_bs

        ! NOW, FULL PRESSURE TENSOR
#ifndef SINGLEFLUID
        pxx(x,y,z) = pxx(x,y,z) + omega * (p_br(1,1))
        pyy(x,y,z) = pyy(x,y,z) + omega * (p_br(2,2))
        pzz(x,y,z) = pzz(x,y,z) + omega * (p_br(3,3))
        pxy(x,y,z) = pxy(x,y,z) + omega * (p_br(1,2))
        pxz(x,y,z) = pxz(x,y,z) + omega * (p_br(1,3))
        pyz(x,y,z) = pyz(x,y,z) + omega * (p_br(2,3))
#ifndef NOSURFACTANT
        pxx(x,y,z) = pxx(x,y,z) +                                   &
                    omega * (p_bs(1,1) + p_ss(1,1) + p_rs(1,1))
        pyy(x,y,z) = pyy(x,y,z) +                                   &
                    omega * (p_bs(2,2) + p_ss(2,2) + p_rs(2,2))
        pzz(x,y,z) = pzz(x,y,z) +                                   &
                    omega * (p_bs(3,3) + p_ss(3,3) + p_rs(3,3))
        pxy(x,y,z) = pxy(x,y,z) +                                   &
                    omega * (p_bs(1,2) + p_ss(1,2) + p_rs(1,2))
        pxz(x,y,z) = pxz(x,y,z) +                                   &
                    omega * (p_bs(1,3) + p_ss(1,3) + p_rs(1,3))
        pyz(x,y,z) = pyz(x,y,z) +                                   &
                    omega * (p_bs(2,3) + p_ss(2,3) + p_rs(2,3))
#endif
#endif
        if(s == droplet) then
          scalarpressure(x,y,z) = (pxx(x,y,z) + pyy(x,y,z) + pzz(x,y,z))/3.
        endif

      end do
    end do
  end do

  if(s == droplet) then
    CALL dump_scalar(scalarpressure,'scp')
  else if(s == nondroplet) then
    CALL dump_scalar(pxx,'pxx')
    CALL dump_scalar(pyy,'pyy')
    CALL dump_scalar(pzz,'pzz')
    CALL dump_scalar(pxy,'pxy')
    CALL dump_scalar(pyz,'pyz')
    CALL dump_scalar(pxz,'pxz')
  end if

  !DEBUG
  !do i=1,5
  !   print*,'At (',sample_point(i),',',nyi/2,',',nzi/2, &
  !   '):'
  !   print*,'  pxz_virial = ', &
  !   omega * (p_br(1,3) + p_bs(1,3) + p_ss(1,3) + p_rs(1,3))
  !   print*,'  pxz_kinetic+virial =', &
  !   pxz(sample_point(i),nyi/2,nzi/2)
  !enddo
  !END

  deallocate(pxx)
  deallocate(pyy)
  deallocate(pzz)
  deallocate(pxy)
  deallocate(pyz)
  deallocate(pxz)
  if(s == droplet) then
      deallocate(scalarpressure)
  endif

  deallocate(u)
end subroutine dump_pressure

!> Calculates \f$ p(tz)=(P_zz - P_xx)(tz) \f$ and then calls
!> \c dump_scalar to write it onto file in order to compute:
!> surface tension = \f$ \int pressure(z) dz \f$
!> at postprocessing time.
!>
!> s decides whether including scalar pressure calculation
subroutine dump_popul_zx(N)
  implicit none
  type(lbe_site), dimension(0:,0:,0:) :: N
  real*8, dimension(:,:,:), allocatable :: populz, populx
  integer :: ierror,x,y,z,nxi,nyi,nzi

  nxi = size(N,1)-2
  nyi = size(N,2)-2
  nzi = size(N,3)-2

  allocate(populz(1:nxi,1:nyi,1:nzi),stat=ierror)
  if (ierror .ne. 0) then
    call log_msg("WARNING: unable to allocate scalar buffer populz(:,:,:)",.true.)
    return  ! Inability to write output is not quite fatal.
  end if

  allocate(populx(1:nxi,1:nyi,1:nzi),stat=ierror)
  if (ierror .ne. 0) then
    call log_msg("WARNING: unable to allocate scalar buffer populx(:,:,:)",.true.)
    return  ! Inability to write output is not quite fatal.
  end if


  ! Following is unnecessary if not using accumulators.
  do z=1,nzi
    do y=1,nyi
    do x=1,nxi
        populz(x,y,z) = 0.
        populx(x,y,z) = 0.
    enddo
    enddo
  enddo

  do z=1,nzi
    do y=1,nyi
      do x=1,nxi

        populz(x,y,z) =   N(x,y,z)%n_r(13) + N(x,y,z)%n_r(15) + &
                          N(x,y,z)%n_r(9)  + N(x,y,z)%n_r(17) + &
                          N(x,y,z)%n_r(5)  + &
                          N(x,y,z)%n_r(14) + N(x,y,z)%n_r(16) + &
                          N(x,y,z)%n_r(10) + N(x,y,z)%n_r(18) + &
                          N(x,y,z)%n_r(6)
#ifndef SINGLEFLUID
        populz(x,y,z) =   populz(x,y,z)    + &
                          N(x,y,z)%n_b(13) + N(x,y,z)%n_b(15) + &
                          N(x,y,z)%n_b(9)  + N(x,y,z)%n_b(17) + &
                          N(x,y,z)%n_b(5)  + &
                          N(x,y,z)%n_b(14) + N(x,y,z)%n_b(16) + &
                          N(x,y,z)%n_b(10) + N(x,y,z)%n_b(18) + &
                          N(x,y,z)%n_b(6)
#endif
        populx(x,y,z) =   N(x,y,z)%n_r(7)  + N(x,y,z)%n_r(10) + &
                          N(x,y,z)%n_r(8)  + N(x,y,z)%n_r(9)  + &
                          N(x,y,z)%n_r(1)  + &
                          N(x,y,z)%n_r(11) + N(x,y,z)%n_r(14) + &
                          N(x,y,z)%n_r(12) + N(x,y,z)%n_r(13) + &
                          N(x,y,z)%n_r(2)
#ifndef SINGLEFLUID
        populx(x,y,z) = populx(x,y,z)      + &
                          N(x,y,z)%n_b(7)  + N(x,y,z)%n_b(10) + &
                          N(x,y,z)%n_b(8)  + N(x,y,z)%n_b(9)  + &
                          N(x,y,z)%n_b(1)  + &
                          N(x,y,z)%n_b(11) + N(x,y,z)%n_b(14) + &
                          N(x,y,z)%n_b(12) + N(x,y,z)%n_b(13) + &
                          N(x,y,z)%n_b(2)
#endif
      end do
    end do
  end do

  CALL dump_scalar(populz,'popz')
  CALL dump_scalar(populx,'popx')

  deallocate(populz)
  deallocate(populx)
end subroutine dump_popul_zx
!> \}

!> \{
!> \name PROFILE output

!> dump profiles
!>
!> \param[in] N local lattice chunk with halo of depth \c halo_extent
!>
!> Dumps profiles along the coordinate axis and along and
!> perpendicular to the body force vector (only if this itself is
!> not parallel to one of the coordinate axis). Beside the
!> self-explaining information in the profile header each file
!> consist of one line per profile position containing profile
!> position, averaged fluid velocity vector, averaged fluid
!> density, averaged fluid site concentration, and rock site
!> concentration. If compiled with MD additional columns containing
!> averaged particle velocity, averaged particle center
!> concentration, and averaged particle site concentration are
!> printed. Depending on n_sci_profile and n_sci_profile_dump the
!> dumped result is the time-average of several samples.
subroutine dump_profile(N)
  implicit none

  type(lbe_site),intent(in) :: N(1-halo_extent:,1-halo_extent:,1-halo_extent:)

  integer,parameter :: profile_file_unit=12

  character(len=1024) :: profile_file_name
  integer i, ii, k, l, lp(3), lx(3), opp_lp(3), u, v, w, ierror
  real(kind=rk) :: rp(3), total_sites, weight
  real(kind=rk) :: ifs, irs, nfs, nrs
  real(kind=rk) :: ivf(3), vf(3), imf(3), mf(3), irf, rhof
  real(kind=rk) :: perm
#ifndef SINGLEFLUID
  real(kind=rk) :: iv_r(3), v_r(3), iv_b(3), v_b(3)
  real(kind=rk) :: im_r(3), m_r(3), im_b(3), m_b(3)
  real(kind=rk) :: ir_r, ir_b, rho_r, rho_b
#endif
#ifdef MD
  real(kind=rk) :: ipc, ips, ivp(3), vp(3), npc, nps
#endif

  ! first call of this subroutine: perform initialization
  if ( .not. allocated(prf) ) then
    call log_msg("  Setting up profile (one time only).")
    call setup_dump_profile()
  endif

#ifdef MD
  ! calculate on-site particle center density and velocity from
  ! off-site particle positions
  call log_msg("  Calculating on-site particle centres and velocities.")
  Nvp = 0.0_rk
  Nnpc = 0.0_rk
  i = atompnt
  particles: do ii = 1,nlocal+nother
    ! Iterate through all possible periodic images; this is
    ! surely not efficient but most safe.
    pbc_x: do u=-1,1
      pbc_y: do v=-1,1
        pbc_z: do w=-1,1
          lx = floor(P(i)%x) + (/ u*tnx, v*tny, w*tnz /)
          lattice_points: do l=1,n_lp
            ! lattice point to look at now, these are still
            ! global coordinates
            lp(:) = lx(:) + lp_sur(:,l)

            ! transform to local coordinates, clipping to
            ! local chunk follows
            lp = lp-(start-1)

            ! only calculate data for real nodes and halo
            ! nodes (halo extent 1)
            if ( any(lp<0) .or. any( lp > (/ nx, ny, nz /) + 1 ) ) cycle

            ! weight for each lattice point is the volume of the
            ! hypercube spanned by the particle position and the
            ! opposite lattice point
            opp_lp = lx + opp_lp_sur(:,l)
            weight = abs( product( real( opp_lp - (/ u*tnx, v*tny, w*tnz /), kind=rk) - P(i)%x ) )
            Nvp( lp(1), lp(2), lp(3), : ) = Nvp( lp(1), lp(2), lp(3), : ) + weight*P(i)%v
            Nnpc( lp(1), lp(2), lp(3) ) = Nnpc( lp(1), lp(2), lp(3) ) + weight
          end do lattice_points
        end do pbc_z
      end do pbc_y
    end do pbc_x

    if ( ii <= nlocal ) then
      i = list(i)
    else
      i = i + 1
    endif
  enddo particles
#endif

  profiles: do k = 1,size(prf)
    write(msgstr,"('  Calculating profile <', A , '> ...')") trim(prf(k)%name)
    call log_msg(msgstr)

    prf_direction: do u = 1,prf(k)%size(1)
      avg1_direction: do v = 1,prf(k)%size(2)
        avg2_direction: do w = 1,prf(k)%size(3)
          ! position along profile directions
           rp = prf(k)%s + real(u-1,kind=rk)*prf(k)%d(1,:) &
             + real(v-1,kind=rk)*prf(k)%d(2,:) + real(w-1,kind=rk)*prf(k)%d(3,:)

          ! normalize according to periodic boundaries
          where ( rp < 0.5_rk ) rp = rp + real( (/ tnx, tny, tnz /), kind=rk )
          where ( rp >= real( (/ tnx, tny, tnz /), kind=rk ) + 0.5_rk ) rp = rp - real( (/ tnx, tny, tnz /), kind=rk )

          ! accumulate data only on that process that owns the
          ! node which is closest to  rp
          if ( any( rp < real(start,kind=rk)-0.5_rk) .or. any( rp >= real( start + (/ nx, ny, nz /), kind=rk ) - 0.5_rk ) ) cycle

          ! interpolate mass density, velocity, and site occupation
          call fluid_velocity_and_density_and_site_occupation( N, rp, ivf, imf, irf &
#ifndef SINGLEFLUID
            &, iv_r, im_r, ir_r, iv_b, im_b, ir_b &
#endif
            &, ifs, irs &
#ifdef MD
            &, ips &
#endif
            &)
!          print *, u,v,w,k, ivf

          ! accumulate interpolated quantities
          prf(k)%l(u)%rhof = prf(k)%l(u)%rhof + irf
          prf(k)%l(u)%vf = prf(k)%l(u)%vf + matmul(prf(k)%vd,ivf)
          prf(k)%l(u)%mf = prf(k)%l(u)%mf + matmul(prf(k)%vd,imf)
#ifndef SINGLEFLUID
          prf(k)%l(u)%rho_r = prf(k)%l(u)%rho_r + ir_r
          prf(k)%l(u)%v_r = prf(k)%l(u)%v_r + matmul(prf(k)%vd,iv_r)
          prf(k)%l(u)%m_r = prf(k)%l(u)%m_r + matmul(prf(k)%vd,im_r)
          prf(k)%l(u)%rho_b = prf(k)%l(u)%rho_b + ir_b
          prf(k)%l(u)%v_b = prf(k)%l(u)%v_b + matmul(prf(k)%vd,iv_b)
          prf(k)%l(u)%m_b = prf(k)%l(u)%m_b + matmul(prf(k)%vd,im_b)
#endif
          prf(k)%l(u)%nfs = prf(k)%l(u)%nfs + ifs
          prf(k)%l(u)%nrs = prf(k)%l(u)%nrs + irs
#ifdef MD
          prf(k)%l(u)%nps = prf(k)%l(u)%nps + ips

          ! interpolate particle center density and velocity
          ! at the continuous profile position  rp.
          call interpolate( Nnpc, rp, ipc)
          call interpolate( Nvp, rp, ivp)

          ! accumulate additional interpolated quantities
          prf(k)%l(u)%npc = prf(k)%l(u)%npc + ipc
          prf(k)%l(u)%vp = prf(k)%l(u)%vp + matmul( prf(k)%vd, ivp )
#endif
        end do avg2_direction
      end do avg1_direction
    end do prf_direction
  end do profiles

  n_samples = n_samples + 1

  dump: if (mod(nt,n_sci_profile_dump)==0) then

    do k=1,size(prf)

      write(msgstr,"('  Writing profile <', A , '>...')") trim(prf(k)%name)
      call log_msg(msgstr)

      call MPI_Reduce(prf(k)%l,prfsum(k)%l,size(prf(k)%l),layer_mpitype,sum_layer_mpiop,0,comm_cart,ierror)

      rank0: if (myrankc==0) then
        call lbe_make_filename_output(profile_file_name,'profile-'//trim(prf(k)%name),'.asc',nt)
        open (unit=profile_file_unit,file=profile_file_name,status='REPLACE',action='WRITE',recl=500)
        write (unit=profile_file_unit,fmt='("# direction=(/",2(ES15.8,","),ES15.8,"/)")') prf(k)%d(1,:)
        write (unit=profile_file_unit,fmt='("# position 1 corresponds to (/",2(I5,","),I5,"/)")') prf(k)%s
        write (unit=profile_file_unit,fmt='("# profile length: ",I5)') prf(k)%size(1)
        write (unit=profile_file_unit,fmt='("# averaging lengths: ",I5,"x",I5)') prf(k)%size(2:3)
        write (unit=profile_file_unit,fmt='("# base for velocity components: '//'(/",2(ES15.8,","),ES15.8,"/), '&
            &//'(/",2(ES15.8,","),ES15.8,"/), '//'(/",2(ES15.8,","),ES15.8,"/)")') prf(k)%vd(1,:),prf(k)%vd(2,:),prf(k)%vd(3,:)
        write (unit=profile_file_unit,fmt='("# average of ",I9," samples")') n_samples
        write (unit=profile_file_unit,fmt='("# ")')
#ifndef SINGLEFLUID
#ifdef MD
        write (unit=profile_file_unit,fmt='("#? pos'&
             &//' vf1             vf2             vf3             mf1             mf2             mf3             rhof'&
             &//'            vr1             vr2             vr3             mr1             mr2             mr3             rho_r'&
             &//'           vb1             vb2             vb3             mb1             mb2             mb3             rho_b'&
             &//'           vp1             vp2             vp3             npc             nps'&
             &//'              nfs             nrs'&
             &//'")')
#else
        write (unit=profile_file_unit,fmt='("#? pos'&
             &//' vf1             vf2             vf3             mf1             mf2             mf3             rhof'&
             &//'            vr1             vr2             vr3             mr1             mr2             mr3             rho_r'&
             &//'           vb1             vb2             vb3             mb1             mb2             mb3             rho_b'&
             &//'           nfs             nrs'&
             &//'")')
#endif
#else
#ifdef MD
        write (unit=profile_file_unit,fmt='("#? pos'&
             &//' vf1             vf2             vf3             mf1             mf2             mf3             rhof'&
             &//'            vp1             vp2             vp3             npc             nps'&
             &//'              nfs             nrs'&
             &//'")')
#else
        write (unit=profile_file_unit,fmt='("#? pos'&
             &//' vf1             vf2             vf3             mf1             mf2             mf3             rhof'&
             &//'            nfs             nrs             perm'&
             &//'")')
#endif
#endif

        do i = lbound(prfsum(k)%l,1),ubound(prfsum(k)%l,1)
          ! normalize and write cumulated numbers
          total_sites = prfsum(k)%l(i)%nfs + prfsum(k)%l(i)%nrs
#ifdef MD
          total_sites = total_sites + prfsum(k)%l(i)%nps
#endif

          ! avoid output of NaN if site or particle counts
          ! are zero, instead output zero, this is still
          ! understandable from the output since also the
          ! counts are dumped (as fractions, however).
          if (prfsum(k)%l(i)%nfs==0.0_rk) then
            vf = 0.0_rk
            mf = 0.0_rk
            rhof = 0.0_rk
#ifndef SINGLEFLUID
            v_r = 0.0_rk
            m_r = 0.0_rk
            rho_r = 0.0_rk
            v_b = 0.0_rk
            m_b = 0.0_rk
            rho_b = 0.0_rk
#endif
         else
            vf = prfsum(k)%l(i)%vf/prfsum(k)%l(i)%nfs
            mf = prfsum(k)%l(i)%mf
            rhof = prfsum(k)%l(i)%rhof/prfsum(k)%l(i)%nfs
#ifndef SINGLEFLUID
            v_r = prfsum(k)%l(i)%v_r/prfsum(k)%l(i)%nfs
            m_r = prfsum(k)%l(i)%m_r
            rho_r = prfsum(k)%l(i)%rho_r/prfsum(k)%l(i)%nfs
            v_b = prfsum(k)%l(i)%v_b/prfsum(k)%l(i)%nfs
            m_b = prfsum(k)%l(i)%m_b
            rho_b = prfsum(k)%l(i)%rho_b/prfsum(k)%l(i)%nfs
#endif
          end if

          if (total_sites==0.0_rk) then
            nfs = 0.0_rk
            nrs = 0.0_rk
#ifdef MD
            npc = 0.0_rk
            nps = 0.0_rk
#endif
          else
            nfs = prfsum(k)%l(i)%nfs/total_sites
            nrs = prfsum(k)%l(i)%nrs/total_sites
#ifdef MD
            npc = prfsum(k)%l(i)%npc/total_sites
            nps = prfsum(k)%l(i)%nps/total_sites
#endif
          end if
#ifdef MD
          if (prfsum(k)%l(i)%npc==0.0_rk) then
            vp = 0.0_rk
          else
            vp = prfsum(k)%l(i)%vp/prfsum(k)%l(i)%npc
          end if
#endif

#ifndef MD
#ifndef SINGLEFLUID
          write (unit=profile_file_unit,fmt='(SS,I5.5,X,SP,23(ES15.8,:,X))') i,vf,mf,rhof,v_r,m_r,rho_r,v_b,m_b,rho_b,nfs,nrs
#else
          perm = ( (2.0_rk*tau_r - 1.0_rk) / 6.0_rk ) * ( 3.0_rk * mf(3) * &
               ( real(tnz, kind=rk) - 2.0_rk*real(boundary_width, kind=rk) + 1.0_rk ) ) / ( ( fr - pr ) * &
               ( real(tnx, kind=rk) - 2.0_rk*real(boundary_width, kind=rk) ) * &
               ( real(tny, kind=rk) - 2.0_rk*real(boundary_width, kind=rk) ) )
          write (unit=profile_file_unit,fmt='(SS,I5.5,X,SP,10(ES15.8,:,X))') i,vf,mf,rhof,nfs,nrs, perm
#endif
#else
#ifndef SINGLEFLUID
          write (unit=profile_file_unit,fmt='(SS,I5.5,X,SP,28(ES15.8,:,X))') i,vf,mf,rhof,v_r,m_r,rho_r,v_b,m_b,rho_b,vp,npc,nps,nfs,nrs
#else
          write (unit=profile_file_unit,fmt='(SS,I5.5,X,SP,14(ES15.8,:,X))') i,vf,mf,rhof,vp,npc,nps,nfs,nrs
#endif
#endif
        end do

        close (profile_file_unit)
      end if rank0
    end do
    call log_msg("  Resetting profile.")
    call reset_profile(prf)

    n_samples = 0
  end if dump
end subroutine dump_profile

!> reset accumulation buffers to zero
elemental subroutine reset_profile(p)
  ! Be careful:  intent(out)  would mean to assign a new object and thus
  ! lose everything not explicitly defined here, namely the allocation
  ! status of the array  p%l .
  type(profile),intent(inout) :: p

  p%l(:)%vf(1) = 0.0_rk
  p%l(:)%vf(2) = 0.0_rk
  p%l(:)%vf(3) = 0.0_rk
  p%l(:)%mf(1) = 0.0_rk
  p%l(:)%mf(2) = 0.0_rk
  p%l(:)%mf(3) = 0.0_rk
  p%l(:)%rhof  = 0.0_rk
#ifndef SINGLEFLUID
  p%l(:)%v_r(1) = 0.0_rk
  p%l(:)%v_r(2) = 0.0_rk
  p%l(:)%v_r(3) = 0.0_rk
  p%l(:)%m_r(1) = 0.0_rk
  p%l(:)%m_r(2) = 0.0_rk
  p%l(:)%m_r(3) = 0.0_rk
  p%l(:)%rho_r  = 0.0_rk
  p%l(:)%v_b(1) = 0.0_rk
  p%l(:)%v_b(2) = 0.0_rk
  p%l(:)%v_b(3) = 0.0_rk
  p%l(:)%m_b(1) = 0.0_rk
  p%l(:)%m_b(2) = 0.0_rk
  p%l(:)%m_b(3) = 0.0_rk
  p%l(:)%rho_b  = 0.0_rk
#endif
  p%l(:)%nfs = 0.0_rk
  p%l(:)%nrs = 0.0_rk
#ifdef MD
  p%l(:)%vp(1) = 0.0_rk
  p%l(:)%vp(2) = 0.0_rk
  p%l(:)%vp(3) = 0.0_rk
  p%l(:)%npc = 0.0_rk
  p%l(:)%nps = 0.0_rk
#endif
end subroutine reset_profile

!> initialization of profiles and associated mpi stuff
subroutine setup_dump_profile()
  implicit none

  integer ierror,k,stat
  character,parameter :: axischar(3)=(/'x','y','z'/) ! for file names
  ! unit vectors parallel to coordinate axes
  real(kind=rk), parameter :: e(3,3) = reshape( (/1.0_rk,0.0_rk,0.0_rk,0.0_rk,1.0_rk,0.0_rk,0.0_rk,0.0_rk,1.0_rk/),(/3,3/) )
  real(kind=rk) :: norms(3) ! for choice of base containing body force
  integer :: maxnorm(1)    ! index of direction to build 5th vector

  ! creation of custom mpi data type for type(layer)
#ifndef SINGLEFLUID
#ifndef MD
  integer,parameter :: lm_n=13
#else
  integer,parameter :: lm_n=16
#endif
#else
#ifndef MD
  integer,parameter :: lm_n=7
#else
  integer,parameter :: lm_n=10
#endif
#endif
  integer lm_lengths(lm_n),lm_types(lm_n)
  integer(kind=MPI_ADDRESS_KIND) :: lm_base,lm_addrs(lm_n),lm_displs(lm_n)
  ! external sum_layer

  ! calculate additional profiles for each vector of the
  ! orthogonal base containing the body force vector which is
  ! not parallel to one of the coordinate axes
  select case ( count( (/g_accn_x, g_accn_y, g_accn /) == 0.0_rk ) )
  case (3,2) ! body force switched off or parallel to one axis
    allocate(prf(3),stat=stat)
    call check_allocate(stat,'setup_dump_profile(): prf(3)')
  case (1)             ! body force in plane with two axes
    allocate(prf(5),stat=stat)
    call check_allocate(stat,'setup_dump_profile(): prf(5)')
  case default         ! general case
    allocate(prf(6),stat=stat)
    call check_allocate(stat,'setup_dump_profile(): prf(6)')
  end select

  ! Setup profile names, directions, starting points, and
  ! dimensions. There are 3, 5, or 6 profiles. Choose the same
  ! starting point  prf%s  for all profiles even if that means
  ! that the associated direction might point out of the real
  ! system, but we assume pbc.
  do k=1,3
      prf(k)%name = axischar(k)
      prf(k)%d = cshift(e,k-1)
      prf(k)%size = cshift((/tnx,tny,tnz/),k-1)
      prf(k)%vd(:,:) = e(:,:)
      prf(k)%s = (/1,1,1/)
  end do

  if (size(prf)>=4) then
    ! 4th profile
    prf(4)%name = 'g'
    prf(4)%d(1,:) = unit_vector((/g_accn_x,g_accn_y,g_accn/))
    prf(4)%s = (/1,1,1/)

    ! choose system diagonal as off-axis profile or averaging length
    prf(4)%size(1) = nint(norm(real((/tnx,tny,tnz/), kind=rk)))

    ! fifth profile: direction is the cross product of the body
    ! force vector with that unit vector which results in the
    ! largest absolute value
    do k=1,3
      norms(k) = norm(cross_product((/g_accn_x,g_accn_y,g_accn/),e(k,:)))
    end do

    maxnorm = maxloc(norms)
    prf(5)%name = 'gx'//axischar(maxnorm(1))
    prf(5)%d(1,:) = unit_vector(cross_product((/g_accn_x,g_accn_y,g_accn/),e(maxnorm(1),:)))
    prf(5)%s = (/1,1,1/)

    ! maybe sixth profile
    if (size(prf)==6) then
      prf(6)%name = 'gx'//trim(prf(5)%name)
      prf(6)%d(1,:) = unit_vector(cross_product((/g_accn_x,g_accn_y,g_accn/),prf(5)%d(1,:)))
      prf(6)%s = (/1,1,1/)

      ! profile and averaging directions
      prf(4)%d(2,:) = prf(5)%d(1,:)
      prf(4)%d(3,:) = prf(6)%d(1,:)
      prf(5)%d = cshift(prf(4)%d,1)
      prf(6)%d = cshift(prf(4)%d,2)

      ! profile and averaging dimensions
      prf(4)%size(2:3) = prf(4)%size(1)
      prf(5)%size = prf(4)%size
      prf(6)%size = prf(4)%size

      ! directions for velocity projection
      prf(4)%vd = prf(4)%d
      prf(5)%vd = prf(4)%vd
      prf(6)%vd = prf(4)%vd
    else
      ! make sure  d  and  vd  contain a right-hand-rule-oriented
      ! system also for the case of 5 profiles
      prf(4)%d(2,:) = e(maxnorm(1),:)
      prf(4)%d(3,:) = prf(5)%d(1,:)
      prf(5)%d = cshift(prf(4)%d,2)

      prf(4)%vd = prf(4)%d
      prf(5)%vd = prf(4)%vd

      ! profile and averaging dimensions
      prf(4)%size(2) = prf(1)%size(maxnorm(1))
      prf(4)%size(3) = prf(4)%size(1)
      prf(5)%size = cshift(prf(4)%size,2)
    end if
  end if

  ! allocate accumulation buffers
  do k=1,size(prf)
    allocate(prf(k)%l(prf(k)%size(1)),stat=stat)
    call check_allocate(stat,'setup_dump_profile(): prf(k)%l(prf(k)%size(1))')
  end do

  ! reset buffers
  call reset_profile(prf)

  ! build custom mpi data type for  layer
  lm_lengths(1) = 1    ! start of layer in memory
  lm_types(1) = MPI_LB
  call MPI_Get_address(prf(1)%l(1),lm_addrs(1),ierror)
  lm_lengths(2) = 3    ! fluid velocity
  lm_types(2) = LBE_REAL
  call MPI_Get_address(prf(1)%l(1)%vf(1),lm_addrs(2),ierror)
  lm_lengths(3) = 3    ! mass flow
  lm_types(3) = LBE_REAL
  call MPI_Get_address(prf(1)%l(1)%mf(1),lm_addrs(3),ierror)
  lm_lengths(4) = 1    ! fluid mass density
  lm_types(4) = LBE_REAL
  call MPI_Get_address(prf(1)%l(1)%rhof,lm_addrs(4),ierror)
  lm_lengths(5) = 1    ! number of fluid sites
  lm_types(5) = LBE_REAL
  call MPI_Get_address(prf(1)%l(1)%nfs,lm_addrs(5),ierror)
  lm_lengths(6) = 1    ! number of solid rock sites
  lm_types(6) = LBE_REAL
  call MPI_Get_address(prf(1)%l(1)%nrs,lm_addrs(6),ierror)
#ifndef SINGLEFLUID
  lm_lengths(7) = 3    ! fluid velocity r
  lm_types(7) = LBE_REAL
  call MPI_Get_address(prf(1)%l(1)%v_r(1),lm_addrs(7),ierror)
  lm_lengths(8) = 3    ! mass flow r
  lm_types(8) = LBE_REAL
  call MPI_Get_address(prf(1)%l(1)%m_r(1),lm_addrs(8),ierror)
  lm_lengths(9) = 1    ! fluid mass density r
  lm_types(9) = LBE_REAL
  call MPI_Get_address(prf(1)%l(1)%rho_r,lm_addrs(9),ierror)
  lm_lengths(10) = 3    ! fluid velocity b
  lm_types(10) = LBE_REAL
  call MPI_Get_address(prf(1)%l(1)%v_b(1),lm_addrs(10),ierror)
  lm_lengths(11) = 3    ! mass flow b
  lm_types(11) = LBE_REAL
  call MPI_Get_address(prf(1)%l(1)%m_b(1),lm_addrs(11),ierror)
  lm_lengths(12) = 1    ! fluid mass density b
  lm_types(12) = LBE_REAL
  call MPI_Get_address(prf(1)%l(1)%rho_b,lm_addrs(12),ierror)
#ifdef MD
  lm_lengths(13) = 3    ! particle velocity
  lm_types(13) = LBE_REAL
  call MPI_Get_address(prf(1)%l(1)%vp(1),lm_addrs(13),ierror)
  lm_lengths(14) = 1    ! number of particle centers
  lm_types(14) = LBE_REAL
  call MPI_Get_address(prf(1)%l(1)%npc,lm_addrs(14),ierror)
  lm_lengths(15) = 1    ! number of moving rock sites
  lm_types(15) = LBE_REAL
  call MPI_Get_address(prf(1)%l(1)%nps,lm_addrs(15),ierror)
#endif
#else
#ifdef MD
  lm_lengths(7) = 3    ! particle velocity
  lm_types(7) = LBE_REAL
  call MPI_Get_address(prf(1)%l(1)%vp(1),lm_addrs(7),ierror)
  lm_lengths(8) = 1    ! number of particle centers
  lm_types(8) = LBE_REAL
  call MPI_Get_address(prf(1)%l(1)%npc,lm_addrs(8),ierror)
  lm_lengths(9) = 1    ! number of moving rock sites
  lm_types(9) = LBE_REAL
  call MPI_Get_address(prf(1)%l(1)%nps,lm_addrs(9),ierror)
#endif
#endif
  lm_lengths(lm_n) = 1
  lm_types(lm_n) = MPI_UB
  call MPI_Get_address(prf(1)%l(2),lm_addrs(lm_n),ierror)

  call MPI_Get_address(prf(1)%l(1),lm_base,ierror) ! base address
  lm_displs(1:lm_n) = lm_addrs(1:lm_n) - lm_base

  call MPI_Type_create_struct(lm_n,lm_lengths,lm_displs,lm_types,layer_mpitype,ierror)
  call MPI_Type_commit(layer_mpitype,ierror)

  ! register custom mpi summation operation for  layer
  call MPI_Op_create(sum_layer,.true.,sum_layer_mpiop,ierror)

  ! allocate more profiles just to have the summation buffers
  allocate (prfsum(size(prf)),stat=stat)
  call check_allocate(stat,'setup_dump_profile(): prfsum(size(prf))')
  do k=1,size(prf)
    allocate (prfsum(k)%l(lbound(prf(k)%l,1):ubound(prf(k)%l,1)),stat=stat)
    call check_allocate(stat,'setup_dump_profile(): '//'prfsum(k)%l(lbound(prf(k)%l,1):ubound(prf(k)%l,1))')
  end do

#ifdef MD
  ! allocate arrays to store particle data on-site
  allocate (Nvp(0:nx+1,0:ny+1,0:nz+1,1:3),Nnpc(0:nx+1,0:ny+1,0:nz+1),stat=stat)
  call check_allocate(stat,'Nvp(0:nx+1,0:ny+1,0:nz+1,1:3),Nnpc(0:nx+1,0:ny+1,0:nz+1)')
#endif

  n_samples = 0

  if (myrankc==0) then
    if (n_sci_profile_dump==0) then 
      call error("FATAL ERROR: n_sci_profile_dump == 0 but sci_profile == .true. . Aborting...")
    end if
  end if
end subroutine setup_dump_profile

!> custom mpi reduction operation to sum objects of  type(layer)
subroutine sum_layer(invec,inoutvec,length,type)
  !    use lbe_io_module, only: layer
  implicit none
  type(layer),intent(in) :: invec(length)
  type(layer),intent(inout) :: inoutvec(length)
  integer,intent(in) :: length,type
  integer i

  do i=1,length
    inoutvec(i)%vf(:) = invec(i)%vf(:) + inoutvec(i)%vf(:)
    inoutvec(i)%mf(:) = invec(i)%mf(:) + inoutvec(i)%mf(:)
    inoutvec(i)%rhof = invec(i)%rhof + inoutvec(i)%rhof
#ifndef SINGLEFLUID
    inoutvec(i)%v_r(:) = invec(i)%v_r(:) + inoutvec(i)%v_r(:)
    inoutvec(i)%m_r(:) = invec(i)%m_r(:) + inoutvec(i)%m_r(:)
    inoutvec(i)%rho_r = invec(i)%rho_r + inoutvec(i)%rho_r
    inoutvec(i)%v_b(:) = invec(i)%v_b(:) + inoutvec(i)%v_b(:)
    inoutvec(i)%m_b(:) = invec(i)%m_b(:) + inoutvec(i)%m_b(:)
    inoutvec(i)%rho_b = invec(i)%rho_b + inoutvec(i)%rho_b
#endif
    inoutvec(i)%nfs = invec(i)%nfs + inoutvec(i)%nfs
    inoutvec(i)%nrs = invec(i)%nrs + inoutvec(i)%nrs
#ifdef MD
    inoutvec(i)%vp(:) = invec(i)%vp(:) + inoutvec(i)%vp(:)
    inoutvec(i)%npc = invec(i)%npc + inoutvec(i)%npc
    inoutvec(i)%nps = invec(i)%nps + inoutvec(i)%nps
#endif
  end do
end subroutine sum_layer
!> \}

!> return the interpolated fluid velocity and mass density at global
!> postion  x  in  v  and  rho  and the interpolated fluid site, (solid)
!> rock site, and--for MD--particle site number density in  nrs ,  nfs ,
!> and  nps .  x  must be at some place that is covered by the local  N ,
!> however, this subroutine is smart enough to deal with pbc.
subroutine fluid_velocity_and_density_and_site_occupation(N, x, v, m, rho &
#ifndef SINGLEFLUID
     &, v_r, m_r, rho_r, v_b, m_b, rho_b &
#endif 
     &, nfs, nrs &
#ifdef MD
     &, nps &
#endif
     &)
  type(lbe_site),intent(in) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
  real(kind=rk),intent(in) :: x(3)
  real(kind=rk),intent(out) :: v(3), m(3), rho, nfs, nrs
#ifdef MD
  real(kind=rk),intent(out) :: nps
#endif
  integer l,s, tx, ty, tz
  integer lx(3), lp(3), opp_lp(3) ! lattice point coordinates
  real(kind=rk) :: xx(3), rho_lp, v_lp(3), weight

#ifndef SINGLEFLUID
  real(kind=rk) :: rho_lp_r, v_lp_r(3)
  real(kind=rk) :: rho_lp_b, v_lp_b(3)
  real(kind=rk) :: v_r(3), v_b(3), m_r(3), m_b(3)
  real(kind=rk) :: rho_r, rho_b
#endif
  real(kind=rk) :: FF_r(19), FF_b(19), ar, ab, rho_r2, rho_b2 

  logical Fminz,  Fmaxz

!#ifdef COMMON_VEL_FIX
#ifndef OLD_VEL
  real(kind=rk), dimension(1:3,0:n_spec) :: corr_vel
#else

  real(kind=rk), dimension(0:nx+1,0:ny+1,0:nz+1,3) :: f_r
#ifndef SINGLEFLUID
  real(kind=rk), dimension(0:nx+1,0:ny+1,0:nz+1,3) :: f_b
#ifndef NOSURFACTANT
  real(kind=rk), dimension(0:nx+1,0:ny+1,0:nz+1,3) :: f_s
#endif
#endif
! endif COMMON_VEL_FIX
#endif

!#ifdef COMMON_VEL_FIX
#ifndef OLD_VEL
  corr_vel(:,:) = 0.0_rk
#endif

  call local_coordinates(x(:),xx(:))
  lx(:) = floor(xx(:))

  v(:) = 0.0_rk
  m(:) = 0.0_rk
  rho = 0.0_rk
#ifndef SINGLEFLUID
  v_r(:) = 0.0_rk
  m_r(:) = 0.0_rk
  rho_r = 0.0_rk
  v_b(:) = 0.0_rk
  m_b(:) = 0.0_rk
  rho_b = 0.0_rk
#endif
  nfs = 0.0_rk
  nrs = 0.0_rk
#ifdef MD
  nps = 0.0_rk
#endif

  lattice_points: do l=1,n_lp
    ! weight each lattice point with the volume of the hypercube
    ! spanned by the particle position and the opposite lattice
    ! point
    lp(:) = lx(:) + lp_sur(:,l)
    opp_lp(:) = lx(:) + opp_lp_sur(:,l)
    weight = abs(product(real(opp_lp(:),kind=rk)-xx(:)))

!#ifndef COMMON_VEL_FIX
#ifdef OLD_VEL

    f_r(lp(1),lp(2),lp(3),:) = 0.0_rk
#ifndef SINGLEFLUID
    f_b(lp(1),lp(2),lp(3),:) = 0.0_rk
#ifndef NOSURFACTANT
    f_s(lp(1),lp(2),lp(3),:) = 0.0_rk
#endif
#endif

#endif
    lp_no_rock: if ( is_fluid(N(lp(1),lp(2),lp(3))%rock_state) ) then
      nfs = nfs + weight

!#ifdef COMMON_VEL_FIX
#ifndef OLD_VEL
      call corrected_velocity(N(lp(1),lp(2),lp(3)),lbe_force(:,:,lp(1),lp(2),lp(3)),corr_vel)
      v_lp = corr_vel(:,1)
      
#ifndef SINGLEFLUID
      v_lp_r = corr_vel(:,1)
      v_lp_b = corr_vel(:,2)
#ifndef NOSURFACTANT
      ! this is meaningless
      !v_lp_s = corr_vel(:,3)
#endif
#endif
      
#ifdef BUGGYIFORT11
      ifort_11_bug_dummy = ifort_11_bug_dummy + 1
#endif
      
      ! else COMMON_VEL_FIX
#else     
      
      v_lp(:) = 0.0_rk
#ifndef SINGLEFLUID
      v_lp_r(:) = 0.0_rk
      v_lp_b(:) = 0.0_rk
#endif
      
#ifdef SINGLEFLUID
      call md_calculate_sc_forces(N,lp(1),lp(2),lp(3),f_r)
#else
#ifdef NOSURFACTANT
      call md_calculate_sc_forces(N,lp(1),lp(2),lp(3),f_b,f_r)
#else
      call md_calculate_sc_forces(N,lp(1),lp(2),lp(3),f_b,f_r,f_s)
#endif
#endif

      directions: do s=1,nnonrest ! ignore resting populations
        v_lp(:) = v_lp(:) + c(s,:)*g(s)*&
             &( N(lp(1),lp(2),lp(3))%n_r(s)*amass_r&
#ifndef SINGLEFLUID
             & + N(lp(1),lp(2),lp(3))%n_b(s)*amass_b&
#endif
#ifndef NOSURFACTANT
             & + N(lp(1),lp(2),lp(3))%n_s(s)*amass_s&
#endif
             &)

#ifdef BUGGYIFORT11
        ifort_11_bug_dummy = ifort_11_bug_dummy + 1
#endif

#ifndef SINGLEFLUID
        v_lp_r(:) = v_lp_r(:) + c(s,:)*g(s)*N(lp(1),lp(2),lp(3))%n_r(s)*amass_r
        v_lp_b(:) = v_lp_b(:) + c(s,:)*g(s)*N(lp(1),lp(2),lp(3))%n_b(s)*amass_b
#endif
      end do directions

      ! force correction because Shan Chen forces
      ! Always call this after streaming because the correction is added with a positive sign!
      v_lp(:) = v_lp(:)  + &
           &( + tau_r*f_r(lp(1),lp(2),lp(3),:)/2.0_rk&
#ifndef SINGLEFLUID
           & + tau_b*f_b(lp(1),lp(2),lp(3),:)/2.0_rk&
#endif
#ifndef NOSURFACTANT
           & + tau_s*f_s(lp(1),lp(2),lp(3),:)/2.0_rk&
#endif
           &)

#ifndef SINGLEFLUID
      v_lp_r(:) = v_lp_r(:) + tau_r*f_r(lp(1),lp(2),lp(3),:)/2.0_rk
      v_lp_b(:) = v_lp_b(:) + tau_b*f_b(lp(1),lp(2),lp(3),:)/2.0_rk
#endif
      
! endif COMMON_VEL_FIX
#endif

      rho_lp = sum(g(:)*& ! sum loops through rest vector, too!
           &( N(lp(1),lp(2),lp(3))%n_r(:)*amass_r &
#ifndef SINGLEFLUID
           & + N(lp(1),lp(2),lp(3))%n_b(:)*amass_b &
#endif
#ifndef NOSURFACTANT
           & + N(lp(1),lp(2),lp(3))%n_s(:)*amass_s &
#endif
           &) )

#ifndef SINGLEFLUID
      rho_lp_r = sum(g(:)*N(lp(1),lp(2),lp(3))%n_r(:))*amass_r
      rho_lp_b = sum(g(:)*N(lp(1),lp(2),lp(3))%n_b(:))*amass_b
#endif

      rho = rho + weight*rho_lp
#ifndef SINGLEFLUID
      rho_r = rho_r + weight*rho_lp_r
      rho_b = rho_b + weight*rho_lp_b
#endif

! #ifdef COMMON_VEL_FIX
!       v(:) = v(:) + weight*v_lp(:)
!       m(:) = m(:) + weight*v_lp(:) * rho_lp
! #ifndef SINGLEFLUID
!       v_r(:) = v_r(:) + weight*v_lp_r(:)
!       m_r(:) = m_r(:) + weight*v_lp_r(:) * rho_lp_r
!       v_b(:) = v_b(:) + weight*v_lp_b(:)
!       m_b(:) = m_b(:) + weight*v_lp_b(:) * rho_lp_b
! #endif


! else COMMON_VEL_FIX
! #else
      v(:) = v(:) + weight*v_lp(:)/max(esmall,rho_lp)
      m(:) = m(:) + weight*v_lp(:)
#ifndef SINGLEFLUID
      v_r(:) = v_r(:) + weight*v_lp_r(:)/max(esmall,rho_lp_r)
      m_r(:) = m_r(:) + weight*v_lp_r(:)
      v_b(:) = v_b(:) + weight*v_lp_b(:)/max(esmall,rho_lp_b)
      m_b(:) = m_b(:) + weight*v_lp_b(:)
#endif



!#ifndef COMMON_VEL_FIX
#ifdef OLD_VEL
      ! force extra-term
      tx = lp(1) + ccoords(1)*nx
      ty = lp(2) + ccoords(2)*ny
      tz = lp(3) + ccoords(3)*nz
      if( (tz.ge.g_accn_min)    .and. (tz.le.g_accn_max)   .and. &
          (tx.ge.g_accn_min_x)  .and. (tx.le.g_accn_max_x) .and. &
          (ty.ge.g_accn_min_y)  .and. (ty.le.g_accn_max_y) ) then
          ! external force term, in the zone where the external force works
          v(:) = v(:) - weight*(/g_accn_x,g_accn_y,g_accn/)/2.0_rk
          m(:) = m(:) - weight*rho_lp*(/g_accn_x,g_accn_y,g_accn/)/2.0_rk
#ifndef SINGLEFLUID
          v_r(:) = v_r(:) - weight*(/g_accn_x,g_accn_y,g_accn/)/2.0_rk
          m_r(:) = m_r(:) - weight*rho_lp_r*(/g_accn_x,g_accn_y,g_accn/)/2.0_rk
          v_b(:) = v_b(:) - weight*(/g_accn_x,g_accn_y,g_accn/)/2.0_rk
          m_b(:) = m_b(:) - weight*rho_lp_b*(/g_accn_x,g_accn_y,g_accn/)/2.0_rk
#endif
      endif
#endif

#ifdef MD
    else if ( interaction=='ladd' .and. is_colloid(N(lp(1),lp(2),lp(3))%rock_state) ) then
      nps = nps + weight
#endif
    else lp_no_rock
      nrs = nrs + weight
    end if lp_no_rock
  end do lattice_points
  
#ifdef BUGGYIFORT11
  ifort_11_bug_dummy = ifort_11_bug_dummy / 2
#endif
end subroutine fluid_velocity_and_density_and_site_occupation

end module lbe_io_module
