!!-*-f90-*-
&md_input

!rho = 10.0
!mass = -1.4                     ! mass of an md particle
!mass = 110587.746982998          ! mass of an md particle
!inertia_orth = 8449790.38511598
!inertia_para = 8449790.38511598
fix_v = .true.
fix_w = .true.

!potential = 'hertz'
rc = 69.0 
rs = 69.5
ineigh=0

n_dump = 10 
dump_format = 'asc'
prtcl_output = .true.
time_output = .false.
dump_time = .true.
dump_double = .false.
dump_positions = .false. !.true.
dump_velocities = .false. !.true.
dump_forces = .true.
dump_quaternions = .false.
dump_orientations = .false. !.true.
dump_rotations = .false.
dump_torques = .true.
dump_ids = .false.
dump_mag = .false.
n_stat = 0

initial_placing = 'file:asc/xvow'
!initial_placing = 'fcc'
!alat = 41.5692193817 
!alat = 32 
!alat = 55.4256258422 
!alat = 64 
!min_dist = 10.0
init_file = 'init.cfg'

ramping = .false.
!ramptime = 50 !5000
vmax_ramping = 0.0 0.0 0.0 ! 1.75430071e-2  ! 1.75430071e-003
forcebalance = .false.
initial_velocities = 'none'
initial_orientations = 'none'
average_ft_fluid = .false.


potential = 'none'   ! particle-particle potential (none|bprh|bprhag|gb|lj|hertz)
interaction = 'ladd' ! particle-fluid interaction (none|friction|ladd|tracer)
rock = 'none'                     ! particle-rock interaction (none|bprh|lj)

constant_f = 0.0 0.0 0.0
!constant_f = 0.0 0.0 0.0003657744
!constant_t = 0.0 0.0 0.0
!cargo_par = .true.
!cargo_par_id = 1
!constant_fc = 0.0 0.0 0.0
!constant_tc = 0.0 0.0 0.0

&md_potential_hertz
K_hertz = 100.0
d_hertz =10
/

&md_fluid_ladd
R_orth = 15.0 !18.0 !4.5  !14.25
R_para = 60.0 !18.0 !4.5  !14.25
!R_orth = 7.0
!R_para = 17.5
!lubrication = .false.
particle_colour = 0.0
!mass_correction = .true.
!C0=10000
/

