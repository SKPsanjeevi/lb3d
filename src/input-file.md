!!-*-f90-*-
&md_input
! steps_per_lbe_step =1   ! # md steps performed after each lbe step
! mass = 0.3 !1.423713094                     ! mass of an md particle
rho = 1.0
! principal moments of inertia orthogonal and parallel to axis of rotational
! symmetry
! inertia_orth = 0.2
! inertia_para = 0.2

ineigh = 1                      ! Neighboring: (0) N^2 (1) Binned
rc = 5.7                       ! inner cutoff
rs = 5.9                      ! outer cutoff

prtcl_output = .true.     ! when true, the program writes one file per particle
time_output = .true.     ! when true, the program writes one file per timestep

n_dump = 1 ! output configuration every this many steps (0 means 'never')
dump_format = 'asc'             ! dump file format (asc|vtk|xdr)
dump_double = .true.            ! dump doubles or floats?
dump_time = .true.            ! dump time?
dump_positions = .true.         ! include positions into output?
dump_pnd = .false.              ! dump PND as output ? 
dump_velocities = .true.        ! include velocities into output?
dump_forces = .true.            ! include forces into output?
!dump_quaternions = .false.      ! include quaternions into output?
! include unit vector parallel to symmetry axis into output?
dump_orientations = .true.
dump_rotations = .true.         ! include angular velocities into output?
dump_torques = .true.           ! include torques into output?
dump_ids = .true.               ! include particle id into output?

n_stat = 0                 ! status calculation every this many steps

!n_msdx = 0
!n_msdx_sample = 1
!msdx_len = 0
!
!msdx_name = 'wide'
!msdx_minx = 0.5
!msdx_maxx = 32.5
!
!n_msdy = 0
!n_msdy_sample = 1
!msdy_len = 0
!
!msdy_name = 'widey'
!msdy_minx = 0.5
!msdy_maxx = 32.5
!msdy_subtract_drift = .true.
!
!n_msdz = 0
!n_msdz_sample = 1
!msdz_len = 0
!
!msdz_name = 'widez'
!msdz_minx = .5
!msdz_maxx = 32.5

! specifies initial particle positions
! (sc|fcc|face[xyz]|file:asc/xvow|file:asc/xvqw|lyapunov|random)
! initial_placing = 'file:asc/xvow'
initial_placing = sc
! path to file to read initial configuration from (when required)
!init_file = 'init.cfg'
alat = 2
! lattice constant, interpretation depends on  initial_placing !
! min/max position for particle placement, clipped to  minpos(:)/maxpos(:) .
! Every particle consumes a space of  alat  in every direction, so
! for  all(x0hi(:)==x0lo(:)+alat)  only 1 particle is placed instead of 8.

! specifies kind of velocity initialization (none|temp0|v0)
! initial_velocities = 'temp0'
! v0 = 0.0 0.0 0.1          ! initial velocity vector for all particles
! temp0 = 0.00            ! initial temperature (reduced units)

!lattice units in SI-units
! delta_x = 7.0E-4	!in meters
! delta_t = 3.15E-7	!in seconds
! delta_m = 8.38978E-15	!in kg
! molec_mass = 1.660538782E-25 ! in KG
! molec_mass_lbm = 3.321077564E-27 ! in KG
! temperat = 295.0
! reflective_rocks=.false.
! mean_free_path = 1.5732E-3

! rdm_seed = 84692
! reenter_rdm_min = .false. .false. .false.
! reenter_rdm_max = .false. .false. .false.

! semiper_max = .false. .false. .false.
! semiper_min = .false. .false. .false.

! diffuse_x = .true.
! diffuse_y = .true.
! diffuse_z = .true.
! count_periodic = .false.

! kind of orientation initialization (none|q0)
! initial_orientations = 'q0'
! ! initial orientation for all particles (rotation defined by this quaternion is
! ! applied to (0 0 1) vector as the particles axis of symmetry)
! q0 = 1.0 0.0 0.0 0.0

! specifies kind of angular velocity initialization (none|w0)
! initial_rotations = 'w0'
! initial angular velocity vector for all particles
! w0 = 0.0 0.0 0.0

! hold (angular) velocity fixed?
! fix_v = .false.
! fix_w = .false.

 constant_f = 0.0 0.0 0.00 ! constant force on all particles (like gravitation)
! constant_t = 0.0 0.0 0.0 ! constant torque on all particles

potential = 'none' ! particle-particle potential (none|bprh|bprhag|gb|lj|dlvo|hertz)
interaction = 'ladd' ! particle-fluid interaction (none|friction|ladd|tracer)
rock = 'none'                     ! particle-rock interaction (none|bprh|lj|dlvo)

! lyapunov = .false.         ! switch lyapunov code on/off
/

&md_fluid_ladd
R_orth = 3.0
R_para = 3.0
particle_colour = 0.0
/
