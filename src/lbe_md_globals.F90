#include "lbe.h"

!> Molecular dynamics global variables module
module lbe_md_globals_module
#ifdef MD

    use lbe_globals_module, only: rk
    use map_module, only: Mii_type

    implicit none
    public

    !> io unit for md input file, must not be the number used for the
    !> lbe input file (currently 17) because that file stays open for
    !> several function calls
    integer,parameter :: md_input_file_unit=18

    integer,parameter :: nnmax=400 !< max # of neighbors per owned atom
    integer,parameter :: nsmax=250 !< max # of swaps to do at each timestep

    !> a type for md particles:
    type md_particle_type
       real(kind=rk),dimension(3) :: x   !< position
       !> orientation (quaternion)
       !>
       !> rotation defined by this quaternion is applied to (0 0 1)
       !> vector as the particle's axis of symmetry
       real(kind=rk),dimension(0:3) :: q
       real(kind=rk),dimension(0:3) :: qnew !< updated \c q
       real(kind=rk),dimension(3) :: v      !< velocity
       real(kind=rk),dimension(3) :: w !< body fixed angular velocity
       real(kind=rk),dimension(3) :: wnew !< updated \c w
       real(kind=rk),dimension(3) :: ws !< space fixed angular velocity
#ifdef RWALK
       real(kind=rk),dimension(3) :: v_r !< for \c RWALK (ask David Sinz)
       real(kind=rk)              :: sdx !< for \c RWALK (ask David Sinz)
#endif
       real(kind=rk),dimension(3) :: f !< space fixed force
       real(kind=rk),dimension(3) :: t !< space fixed torque
       
!#ifdef MD_MAG
       real(kind=rk), dimension(2) :: mag        ! magnitude and orientation angle of magnetic dipole
       real(kind=rk),dimension(3) :: fmi ! < force by magnetic dispole interaction
       real(kind=rk),dimension(3) :: tmi ! < torque by magnetic dispole interaction
!#endif
       real(kind=rk),dimension(3) :: freeze ! fix particle due to friction force. 
       real(kind=rk),dimension(3) :: cf !< constant force on particles
       real(kind=rk),dimension(3) :: ct !< constant torque on particles

       !> \name elements required if \c decouple_fluid is set
       !> \{
       !> space fixed force from fluid interaction
       real(kind=rk),dimension(3) :: f_fluid
       !> space fixed torque from fluid interaction
       real(kind=rk),dimension(3) :: t_fluid
       !> space fixed force from fluid interaction in previous LB step
       real(kind=rk),dimension(3) :: f_fluid_prev
       !> space fixed torque from fluid interaction in previous LB step
       real(kind=rk),dimension(3) :: t_fluid_prev
       !> velocity averaged during last LB step
       real(kind=rk),dimension(3) :: v_fluid
       !> space fixed angular velocity averaged during last LB step
       real(kind=rk),dimension(3) :: ws_fluid
       !> velocity used for fluid interaction in LB steps, that is,
       !> average of current and previous \c v_fluid (unless \c
       !> average_vw_fluid is \c .false.)
       real(kind=rk),dimension(3) :: v_fluid_avg
       !> space fixed angular velocity used for fluid interaction in
       !> LB steps, that is, average of current and previous \c
       !> ws_fluid (unless \c average_vw_fluid is \c .false.)
       real(kind=rk),dimension(3) :: ws_fluid_avg
       !> accumulator to compute next \c v_fluid from \c v during substeps
       real(kind=rk),dimension(3) :: v_fluid_acc
       !> accumulator to compute next \c ws_fluid from \c ws during substeps
       real(kind=rk),dimension(3) :: ws_fluid_acc
       !> \}

       !> \name elements required if \c polydispersity is set
       !> \{
       !> particle radius orthogonal to rotational symmetry axis,
       !> replaces the global \c R_orth in \c
       !> lbe_md_fluid_ladd_parms_module
       real(kind=rk) :: R_orth
       !> particle radius parallel to rotational symmetry axis,
       !> replaces the global \c R_para in \c
       !> lbe_md_fluid_ladd_parms_module
       real(kind=rk) :: R_para
       !> \}
#ifdef FORCECOMPONENT
       real(kind=rk), dimension(3) :: f_normal
       real(kind=rk), dimension(3) :: f_tangent
       real(kind=rk), dimension(3) :: f_n
       real(kind=rk), dimension(3) :: f_t
#endif
       

       real(kind=rk),dimension(3) :: dx !< displacement since neighbor update
       real(kind=rk),dimension(3) :: o !< orientation (unit vector)
       real(kind=rk),dimension(3) :: vnew !< new \c v
       real(kind=rk) :: e_pot             !< potential energy
#ifdef PARTICLESTRESS
       real(kind=rk),dimension(3,3) :: tau !< body fixed cumulative stress
#endif
       integer :: uid           !< unique particle id
       integer :: master        !< uid of possible master particle
#ifdef LADD_SURR_RHOF
       real(kind=rk) :: rhof_acc !< accumulated surrounding rhof
       real(kind=rk) :: rhof     !< surrounding rhof (previous time step)
       integer :: n_acc          !< number of samples in \c rhof_acc
#endif
#ifdef LADD_SSD
       !> \name friction (inverse mobility) matrices for \c 'ladd'
       !> particle sub-step damping
       !> \{
       !> force due to \c v-v_fluid_avg
       real(kind=rk),dimension(3,3) :: FvF(3,3)
       !> force due to \c ws-ws_fluid_avg
       real(kind=rk),dimension(3,3) :: FwF(3,3)
       !> torque due to \c v-v_fluid_avg
       real(kind=rk),dimension(3,3) :: FvT(3,3)
       !> torque due to \c ws-ws_fluid_avg
       real(kind=rk),dimension(3,3) :: FwT(3,3)
       !> \}
#endif
    end type md_particle_type

    !> type used for communcation of partial forces and torques
    type ft_buffer_type
       real(kind=rk),dimension(3) :: f !< force
       real(kind=rk),dimension(3) :: t !< space fixed torque
#ifdef FORCECOMPONENT
       real(kind=rk), dimension(3) :: f_normal
       real(kind=rk), dimension(3) :: f_tangent
       real(kind=rk), dimension(3) :: f_n
       real(kind=rk), dimension(3) :: f_t
#endif
       !> space fixed force from fluid interaction (with \c decouple_fluid )
       real(kind=rk),dimension(3) :: f_fluid
       !> space fixed torque from fluid interaction (with \c decouple_fluid )
       real(kind=rk),dimension(3) :: t_fluid
       real(kind=rk) :: e_pot          !< potential energy
#ifdef PARTICLESTRESS
       real(kind=rk),dimension(3,3) :: tau !< body fixed cumulative stress
#endif
#ifdef LADD_SURR_RHOF
       real(kind=rk) :: rhof_acc !< accumulated surrounding rhof
       integer :: n_acc          !< number of samples in \c rhof_acc
#endif
    end type ft_buffer_type

    !> contains what is needed for \c exchange()
    integer,save :: particle_exch_mpitype
    !> contains what is needed for \c communicate()
    integer,save :: particle_comm_mpitype
    !> contains what is needed for \c collect()
    integer,save :: particle_coll_mpitype
    !> contains what is needed for \c collect() in timesteps where \c
    !> md-cfg data is dumped
    integer,save :: particle_coll_dump_mpitype
    !> contains what is needed for \c collect() in MD substeps
    integer,save :: particle_coll_substep_mpitype
    !> contains what is necessary to dump or restore an MD checkpoint
    integer,save :: particle_cp_mpitype

    !> representation of \c ft_buffer_type for \c collect()
    integer,save :: ft_buffer_coll_mpitype
    !> representation of \c ft_buffer_type for \c collect() in
    !> timesteps where \c md-cfg data is dumped
    integer,save :: ft_buffer_coll_dump_mpitype
    !> representation of \c ft_buffer_type for \c collect() in MD substeps
    integer,save :: ft_buffer_coll_substep_mpitype

    !> number of (pseudo-)directions
    integer,parameter :: n_dir_max=8

    !> \name definition of (pseudo-)direction indices
    !> \{
    !> send in negative coordinate direction
    integer,parameter :: DIR_DOWN=1
    !> send in positive coordinate direction
    integer,parameter :: DIR_UP=2
    !> negative direction, non-periodic particle (for \c periodic_inflow)
    integer,parameter :: DIR_NONPERIODIC_DOWN=3
    !> positive direction, non-periodic particle (for \c periodic_inflow)
    integer,parameter :: DIR_NONPERIODIC_UP=4
    !> negative direction across Lees-Edwards plane into neighbor with lower z
    integer,parameter :: DIR_LE_LO_Z_DOWN=5
    !> positive direction across Lees-Edwards plane into neighbor with lower z
    integer,parameter :: DIR_LE_LO_Z_UP=6
    !> negative direction across Lees-Edwards plane into neighbor with higher z
    integer,parameter :: DIR_LE_HI_Z_DOWN=7
    !> positive direction across Lees-Edwards plane into neighbor with higher z
    integer,parameter :: DIR_LE_HI_Z_UP=8
    !> \}

#ifdef DEBUG_REPORTMDCOMM
    !> human-readable names for (pseudo-)directions
    character(len=*),parameter :: DIR_NAMES(n_dir_max)=(/'DOWN','UP','NP_DOWN'&
         &,'NP_UP','LE_LO_Z_DOWN','LE_LO_Z_UP','LE_HI_Z_DOWN','LE_HI_Z_UP'/)
#endif

    !> type for data related to a specific swap
    type swap_type
       integer :: s_mpitype !< specifies all elements of \c P to be sent

       real(kind=rk) :: sbndlo !< low slab boundary on atom positions to send
       real(kind=rk) :: sbndhi !< high slab boundary on atom positions to send

       integer :: sicnt      !< number of indexes in \c sidx
       integer :: sipos      !< position of indexes in \c sidx
       integer :: rpcnt      !< number of particles in \c P received
       integer :: rppos      !< position of particles in \c P received
    end type swap_type

    !> type for data related to a specific swap direction
    type swap_direction_type
       integer :: dim           !< coordinate dimension index
       integer :: dir           !< (pseudo-)direction index
       integer :: sproc         !< rankc of cpu to send to
       integer :: rproc         !< rankc of cpu to recv from
       integer :: scnt          !< number of sends
       integer :: rcnt          !< number of recvs

       real(kind=rk) :: blo_z   !< special low z slab boundary (where needed)
       real(kind=rk) :: bhi_z   !< special high z slab boundary (where needed)

       type(swap_type) :: s(nsmax) !< data specific for each swap
    end type swap_direction_type

    !> data specific for swaps in each coordinate dimension and
    !> (pseudo-)direction
    type(swap_direction_type),save :: sdirs(3,n_dir_max)

    !> for building \c sdirs(:,:)%swaps(:)%s_mpitype and for adding
    !> the forces to the right particles in \c collect()
    integer,allocatable,save :: sidxs(:)

    !> for \c mpi_type_indexed() in \c exchange() and \c borders()
    integer,allocatable,save :: slens(:)

    !> particles local on each processor:
    type(md_particle_type),allocatable,save :: P(:)

    !> particle recv buffer for  exchange() :
    type(md_particle_type),allocatable,save :: pbuf(:)

    !> force/torque recv buffer for  collect() :
    type(ft_buffer_type),allocatable,save :: ftbuf(:)

    !> indicates whether \c md_init() was run already
    logical,save :: md_initialized=.false.

    !> changed by  setup_fluid_<couplingname>  if necessary
    logical,save :: halo_exchange_before_md_run=.false.
    !> changed by  setup_fluid_<couplingname>  if necessary
    logical,save :: halo_exchange_after_md_run_middle=.false.

    !> interval to calculate z-forces acting on chunk dfz_minx<=x<=dfz_maxx
    integer,save :: n_dfz=0
    !> lowest x coordinate of chunk for which z-forces are calculated
    integer,save :: dfz_minx=-1
    !> highest x coordinate of chunk for which z-forces are calculated
    integer,save :: dfz_maxx=-1
    !> \name z-momentum transfer due to various effects
    !> \{
    !> z-forces/transferred momentum acting through  dfz_minx  (1) and
    !> dfz_maxx  (2) due to fluid advection (f), particle advection (p),
    !> fluid-particle coupling (fp), and particle-particle interaction
    !> (pp). All variables hold only the local contributions of each
    !> process.
    real(kind=rk),save :: dfz_f(1:2),dfz_p(1:2),dfz_fp(1:2),dfz_pp(1:2)
    !> \}

    !> are operations necessary for non-point particles performed?
    logical,save :: use_rotation=.false.

    !> single particle: gradient matrix to compute surface normal
    !> of ellipsoid.
    !> CAUTION: used only for single particle. For multiparticles
    !> will return garbage values.
    real(kind=rk), dimension(3,3) :: gradient_matrix

    !> indicates whether particle's own R_orth and R_para are used,
    !> supposed to be enabled automatically if polydisperse systems
    !> are initialized
    logical,save :: polydispersity=.false.

     !> indicates whether the different magnitude of magnetic dipole is used or not.
      logical,save :: magdispersity=.false.

    !> are \c md_particle_type elements \c f_fluid[_prev] and \c
    !> t_fluid[_prev] in use?
    logical,save :: use_ft_fluid=.false.

    !> decouple particle-fluid forces from other particle forces
    !> assuming that the first ones do not change during one LB time
    !> step while the latter ones can change during the MD substeps
    !> (only useful in case \c steps_per_lbe_step>1 )
    logical,save :: decouple_fluid=.false.

    !> \name options for communicate()
    !> \{
    !> specify whether velocities, angular velocities, angular
    !> velocities in space fixed frame, and uids respectively are
    !> communicated along with the positions in communicate()
    logical,save :: communicate_velocities=.false.
    logical,save :: communicate_rotations=.false. !< body fixed \c w
    logical,save :: communicate_rotations_s=.false. !< space fixed \c w
    logical,save :: communicate_uids=.false.
#ifdef RWALK
    logical,save :: communicate_velocities_r=.true.
    logical,save :: communicate_dsx=.true.
#endif
    !> \}

    real(kind=rk),save :: domaindx !< for \c RWALK (ask David Sinz)

    !> specify whether forces (and torques, if use_rotation is set) of
    !> halo'ed particles are sent back to the owning processor
    logical,save :: collect_forces=.false.

    !> specify whether orientation vectors  P(:)%o(1:3)  are calculated from the
    !> quaternions  P(:)%q(0:3)  each timestep (they are calculated for own and
    !> halo'ed particles all the same)
    logical,save :: calculate_orientations=.false.

    !> specify whether uid2i is provided
    logical,save :: provide_uid2i=.false.

    !> local lookup table converting from particle uid to local index
    !> position within P
    type(Mii_type),save :: uid2i

    character(len=32),save :: interaction='none' !< fluid-particle interaction

    integer,save :: steps_per_lbe_step=1 !< # md steps per lbe step

    !> average particle force and torque exerted by the fluid over two
    !> consecutive time steps?
    logical,save :: average_ft_fluid=.true.

    integer,save :: nt_substep=0 !< substep counter (\c 1..steps_per_lbe_step)

    !> \name global sim space boundaries
    !> \{
    real(kind=rk),save :: minx,miny,minz,maxx,maxy,maxz
    !> \}

    integer,save :: nemax !< max # of atoms to send to all neighbors (swap list)
    integer,save :: nfmax       !< size of  pbuf(:)
    integer,save :: nomax       !< max # of nearby atoms
    !> maximum number of owned particles, halo'ed particles start at P(npmax+1)
    integer,save :: npmax
    integer,save :: ntmax       !< max # of thermodynamic calls
    !> swap list of atoms to send out for each swap
    integer,allocatable,save :: slist(:)

    !> linked list of local atoms (last one -> npmax+1)
    integer,allocatable,save :: list(:)
    integer,allocatable,save :: nlist(:) !< neighbor lists of my atoms
    !> pointers to start of neighbor lists for my atoms
    integer,allocatable,save :: nnlist(:)
    !> linked list pointers from one atom to next in bin
    integer,allocatable,save :: bin(:)
    integer,allocatable,save :: binpnt(:) !< pointer to 1st atom in each bin

    !> if \c .true. , enforces a neighbor list update as soon as possible
    logical,save :: list_update_required=.false.

    !> particle displacement that is effective for deciding about list
    !> updates in addition to the displacement obtained in the
    !> integrator
    real(kind=rk),save :: list_update_displ_add=0.0_rk

    real(kind=rk),allocatable,save :: tmparr(:) !< temperature (reduced units)
    real(kind=rk),allocatable,save :: engarr(:) !< particle-particle potential
    real(kind=rk),allocatable,save :: rpotarr(:) !< rock-particle potential
    real(kind=rk),allocatable,save :: prsarr(:) !< presssure (in reduced units)
    real(kind=rk),allocatable,save :: conarr(:) !< energy conservation
    real(kind=rk),allocatable,save :: momentumarr(:,:) !< total md momentum
    real(kind=rk),save :: tmpave       !< average temperature
    real(kind=rk),save :: engave       !< average particle-particle potential
    real(kind=rk),save :: rpotave      !< average rock-particle potential
    real(kind=rk),save :: prsave       !< average presssure
    real(kind=rk),save :: conave       !< average energy conservation
    real(kind=rk),save :: momentumave(3) !< average total md momentum

    !> evolution of translational kinetic energy
    real(kind=rk),allocatable,save :: e_trans(:)
    !> evolution of rotational kinetic energy
    real(kind=rk),allocatable,save :: e_rot(:)
    !> average translational kinetic energy
    real(kind=rk),save :: e_trans_ave
    !> average rotational kinetic energy
    real(kind=rk),save :: e_rot_ave

    integer,save :: atompnt !< pointer to 1st atom in my list
    integer,save :: freepnt !< pointer to 1st free space in list (last one -> 0)
    !> \name indices referencing timers
    !> \{
    integer,save :: ti_md_comm,ti_md_fluid,ti_md_force,ti_md_neigh,ti_md_dump&
         &,ti_md_rock
    !> \}
    real(kind=rk),save :: rc=-1.0_rk    !< inner (force) cutoff
    real(kind=rk),save :: rs=-2.0_rk    !< outer (neighbor) cutoff
    integer,save :: n_stat=100 !< status calculation interval in lbe steps
    integer,save :: ineigh=1    !< neighbor flag - N^2 or binning
    integer,save :: nbinx !< # of global neighbor bins in each dimension
    integer,save :: nbiny !< # of global neighbor bins in each dimension
    integer,save :: nbinz !< # of global neighbor bins in each dimension

    !> specifies initial particle positions
    character(len=32),save,public :: initial_placing='sc'
    !> lattice constant, interpretation depends on  initial_placing
    real(kind=rk),save,public :: alat=8.0_rk

    !> mass density of an md particle (only used to calculate  mass  with
    !> R_para and R_orth (interaction=='ladd') if mass was set to -1)
    real(kind=rk),save :: rho=1.0_rk

    !> mass of an md particle (negative value means initialization
    !> based on rho, R_para, and R_orth (interaction=='ladd'))
    real(kind=rk),save :: mass=-1.0_rk

    !> \name principal moments of inertia
    !> \{
    !> negative value means initialization based on mass,
    !> R_para, R_orth (interaction=='ladd'), and the assumption of a
    !> homogenous mass density
    real(kind=rk),save :: inertia_orth=-1.0_rk,inertia_para=-1.0_rk
    !> \}

    real(kind=rk),save :: cutsq1       !< inner cutoff squared
    real(kind=rk),save :: cutsq2       !< outer cutoff squared

    !> for optimized min. image crit.: tsize(:)-rc
    real(kind=rk),save :: tsize_mrc(3)
    !> for optimized min. image crit.: tsize(:)-rs
    real(kind=rk),save :: tsize_mrs(3)

    integer,save :: nlocal      !< # of atoms I currently own
    integer,save :: nother      !< # of nearby atoms I currently store
    real(kind=rk),save :: binsizex !< size of bins in each dimension
    real(kind=rk),save :: binsizey !< size of bins in each dimension
    real(kind=rk),save :: binsizez !< size of bins in each dimension
    integer,save :: mbinx !< # of bins in my box (with nearby atoms included)
    integer,save :: mbiny !< # of bins in my box (with nearby atoms included)
    integer,save :: mbinz !< # of bins in my box (with nearby atoms included)
    !> 1st global bin indices (offset) at lower left of my box
    integer,save :: mbinxlo
    !> 1st global bin indices (offset) at lower left of my box
    integer,save :: mbinylo
    !> 1st global bin indices (offset) at lower left of my box
    integer,save :: mbinzlo
    real(kind=rk),save :: enginit !< total energy at time 0
    integer,save :: mstat       !< # of times status routine called
    integer,save :: mneigh      !< # of times neighbor routine called
    integer,save :: nlocalmax   !< most atoms I ever owned
    integer,save :: nothermax   !< most nearby atoms I ever stored
    integer,save :: neighmax !< most neighbors every stored in neighbor list
    integer,save :: nslistmax   !< biggest size swap list ever reached
    integer,save :: nexcmax !< most atoms ever leaving my box (in one dimension)
    integer,save :: nswpmax     !< most atoms ever sent in one swap

    !> \name SI units
    !> Probably David knows about this.
    !> \{
    real(kind=rk),save :: delta_x = 1.0E-10
    real(kind=rk),save :: delta_t = 1.0E-10
    real(kind=rk),save :: delta_m = 1.0E-10
    real(kind=rk),save :: molec_mass = 1.66025E-27     	!< molecular mass Kg
    real(kind=rk),save :: temperat = 273.0      !< Temperature K
    real(kind=rk),save :: mean_free_path = 1.0	!< mean free path

    !> molecular mass of the lbm particles  Kg
    real(kind=rk),save :: molec_mass_lbm = 1.66025E-27
    !> \}

    logical,save :: prtcl_output=.false.
    logical,save :: time_output=.true.

    logical,save :: semiper_max(3)=.false.
    logical,save :: semiper_min(3)=.false.
    logical,save :: reenter_rdm_max(3)=.false.
    logical,save :: reenter_rdm_min(3)=.false.

    real,save :: cutoff_v=1000.0

    logical,save :: reflective_rocks = .false.
    logical,save :: diffuse_x=.true.
    logical,save :: diffuse_y=.true.
    logical,save :: diffuse_z=.true.
    logical,save :: count_periodic = .false.

    integer,save :: np_sphere = 0
    integer,save :: np_sphere_cap = 0

    !> set to \c .true. automatically if MD Lees-Edwards code needs to
    !> be run
    logical,save :: md_leesedwards
#endif
end module lbe_md_globals_module
