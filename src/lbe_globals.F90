#include "lbe.h"

!> This module defines constants which all parts of the program will
!> need to refer to such as the lattice vectors.
module lbe_globals_module
  implicit none
  public

  integer, save :: myrankc            !< rank of the CPU (Comm_Cart)

  integer, parameter :: input_file_unit = 17
  integer, parameter :: input_dfile_unit = 42

  integer, parameter :: rk = 8 !< type parameter (kind) of most of the reals

  !> number of lbe particle species
#ifdef SINGLEFLUID
  integer, parameter :: n_spec = 1
#else
#ifdef NOSURFACTANT
  integer, parameter :: n_spec = 2
#else
  integer, parameter :: n_spec = 3
#endif
#endif

  real(kind=rk), parameter :: pi = 3.1415926535897932384626433_rk
  ! > 10^9
  real(kind=rk), parameter :: esmall = 0.000000001_rk

  real(kind=rk), parameter :: no_rock_value = 0.0_rk
#ifdef MD
  !> Sum of forces acting on all particles
  real(kind=rk),dimension(3) :: forcesum = (/0.0_rk,0.0_rk,0.0_rk/)
  real(kind=rk),dimension(3) :: totalforce = (/0.0_rk,0.0_rk,0.0_rk/)
  real(kind=rk), parameter :: rock_value = -1.0_rk
  real(kind=rk), parameter :: rock_value_surface = -2.0_rk
#else
  real(kind=rk), parameter :: rock_value = 1.0_rk
  real(kind=rk), parameter :: rock_value_surface = 2.0_rk
#endif

#ifdef VELCONTROL
  ! VELCONTROL is a feedback control to achieve desired setpoint velocity
  ! using Proportional-Integral-Derivative (PID) control

  !> Difference between setpoint velocity and current average velocity
  real*8, dimension(3)  :: err_0 = 0._rk       ! difference at current time t
  real*8, dimension(3)  :: err_1 = 0._rk       ! difference at t-1
  real*8, dimension(3)  :: err_2 = 0._rk       ! difference at t-2
  real*8, dimension(3)  :: gravity_new = 0._rk ! new gravity computed at time t
  real*8, dimension(3)  :: gravity_old = 0._rk ! gravity computed from t-1 
#endif

  !> \{
  !> \name These parameters are a result solely of the lattice chosen.

  integer, parameter :: D = 4         !< No dimensions of lattice
  !> Number of different lattice vectors.
  !>
  !> Contains the number of independent vectors (19), which
  !> corresponds to the size of the \c n_r, \c n_b, and
  !> \c n_s arrays of each site. When operating on all elements
  !> of such an array, it is strongly advised that the loop run
  !> from \c 1 to \c nvecs, so that the code will not need
  !> changing should another lattice be used.
  integer, parameter :: nvecs = 19    !< No of vectors
  integer, parameter :: ngvecs = 25   !< No of vectors (inc. degen.)
  integer, parameter :: nnonrest = 18 !< No of nonrest vectors
  real(kind=rk) , parameter :: T = 1.0_rk/3.0_rk !< Lattice temp.
  integer, parameter :: nd = 3    !< Number of dimensions of the system.
  !> The index of the rest vector.
  !>
  !> When referring explicitly to the rest particles, it is strongly
  !> advised that an index of \c restvec rather than 19 is used, to
  !> ensure code readability and minimise problems should the lattice
  !> be changed.
  integer, parameter :: restvec = 19
  !> \}
  !>
  !> \{
  !> \name Lattice vectors and associated degeneracies.
  !>
  !> \verbatim
  !> the vectors are the columns of the following matrix
  !>
  !> 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19             (index)
  !>
  !> x,-x, y,-y, z,-z, x+-y, x+-z,-x+-y,-x+-z, y+-z,-y+-z, 0    (short notation)
  !>
  !> 1,-1, 0, 0, 0, 0, 1, 1, 1, 1,-1,-1,-1,-1, 0, 0, 0, 0, 0           (vectors)
  !> 0, 0, 1,-1, 0, 0, 1,-1, 0, 0, 1,-1, 0, 0, 1, 1,-1,-1, 0
  !> 0, 0, 0, 0, 1,-1, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 1,-1, 0
  !>
  !> 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1(projection weights)
  !>
  !> the weights here are not the lattice weights, but the weights
  !> used for the downprojection from a D3Q25 to a D3Q19 lattice
  !>
  !> 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 12  (lattice weights)
  !> \endverbatim

  !> the X components of the lattice vectors
  integer, parameter :: cx(nvecs) = &
    (/ 1,-1, 0, 0, 0, 0, &
       1, 1, 1, 1,-1,-1,-1,-1, 0, 0, 0, 0, &
       0 &
    /)

  !> the Y components of the lattice vectors
  integer, parameter :: cy(nvecs) = &
    (/ 0, 0, 1,-1, 0, 0, &
       1,-1, 0, 0, 1,-1, 0, 0, 1, 1,-1,-1, &
       0 &
    /)

  !> the Z components of the lattice vectors
  integer, parameter :: cz(nvecs) = &
    (/ 0, 0, 0, 0, 1,-1, &
       0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 1,-1, &
       0 &
    /)

  !> degeneracies of the lattice vectors
  integer, parameter :: g(nvecs) = &
    (/ 2, 2, 2, 2, 2, 2, &
       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
       1 &
    /)

  !> lattice weights
  real(kind=rk), parameter :: w0 = 1.0_rk/3.0_rk
  real(kind=rk), parameter :: w1 = 1.0_rk/18.0_rk
  real(kind=rk), parameter :: w2 = 1.0_rk/36.0_rk
  real(kind=rk), parameter :: w(nvecs) = &
    (/ w1, w1, w1, w1, w1, w1, &
       w2, w2, w2, w2, w2, w2, w2, w2, w2, w2, w2, w2, &
       w0 &
    /)

#ifdef LADD_DLUB
  !> absolute lengths of all non-rest lattice directions
  !> (D3Q19-specific)
  real(kind=rk),parameter &
       &:: c_length(1:nnonrest)=(2.0_rk/g(1:nnonrest))**0.5_rk

  !> unit vectors in all non-rest lattice directions
  real(kind=rk),parameter :: c_dir(3,1:nnonrest)=reshape((/&
       &(/cx(1),cy(1),cz(1)/)/c_length(1),&
       &(/cx(2),cy(2),cz(2)/)/c_length(2),&
       &(/cx(3),cy(3),cz(3)/)/c_length(3),&
       &(/cx(4),cy(4),cz(4)/)/c_length(4),&
       &(/cx(5),cy(5),cz(5)/)/c_length(5),&
       &(/cx(6),cy(6),cz(6)/)/c_length(6),&
       &(/cx(7),cy(7),cz(7)/)/c_length(7),&
       &(/cx(8),cy(8),cz(8)/)/c_length(8),&
       &(/cx(9),cy(9),cz(9)/)/c_length(9),&
       &(/cx(10),cy(10),cz(10)/)/c_length(10),&
       &(/cx(11),cy(11),cz(11)/)/c_length(11),&
       &(/cx(12),cy(12),cz(12)/)/c_length(12),&
       &(/cx(13),cy(13),cz(13)/)/c_length(13),&
       &(/cx(14),cy(14),cz(14)/)/c_length(14),&
       &(/cx(15),cy(15),cz(15)/)/c_length(15),&
       &(/cx(16),cy(16),cz(16)/)/c_length(16),&
       &(/cx(17),cy(17),cz(17)/)/c_length(17),&
       &(/cx(18),cy(18),cz(18)/)/c_length(18)&
       &/),(/3,nnonrest/))
#endif

  !> contains the index of the vector pointing in the opposite
  !> direction to c(i).
  integer, parameter :: bounce(nvecs) = &
    (/  2, 1, 4, 3, 6, 5, &
       12,11,14,13, 8, 7, 10, 9, 18,17,16,15, &
       19 &
    /)

  !> contains the index of the vector pointing in the specular-reflection
  !> direction to c(i).

  integer, parameter :: reflect_y(nvecs) = &
    (/  2,  1, 4 , 3,  6 , 5 , &
        8,  7, 14, 13, 12, 11, &
        10, 9, 17, 18, 15, 16, &
        19 &
    /)

!!$  integer, parameter :: reflect_z(nvecs) = &
!!$    (/  2, 1, 4 , 3, 6 , 5 , &
!!$        8, 7, 10, 9, 12, 11, &
!!$        14, 13, 16,15,18,17, &
!!$        19 &
!!$    /)


  integer, parameter :: nnp = 5 !< size of the \c negx,\c negy, etc arrays

  !> \{
  !> \name Arrays used in invasive flow.
  !>
  !> \c negx contains a list of the indices of each vector which
  !> has a negative X componentl \c posy contains a list of the
  !> indices of each vector which has a positive Y component, etc.
  integer, parameter :: negx(nnp) = (/ 1 , 7 , 8 , 9 , 10 /)
  integer, parameter :: negy(nnp) = (/ 3 , 7 , 11 , 15 , 16 /)
  integer, parameter :: negz(nnp) = (/ 5 , 9 , 13 , 15 , 17 /)
  integer, parameter :: posx(nnp) = (/ 2 , 11 , 12 , 13 , 14 /)
  integer, parameter :: posy(nnp) = (/ 4 , 8 , 12 , 17 , 18 /)
  integer, parameter :: posz(nnp) = (/ 6 , 10 , 14 , 16 , 18 /)
  !> \}

  !> Array containing the lattice vectors.
  !>
  !> contains a vector corresponding to the ith lattice vector, but has
  !> to be initialised in code, since you can't initialise an array of
  !> rank greater than 1.  Fortran sucks.
  integer, save :: c(nvecs,3)
  !> \}

  !> depth of halo region in each direction
  !>
  !> meant to be enlarged on demand. Currently, only in the MD branch
  !> the value is enlarged and currently only in the MD branch the
  !> enlargement has some effect.
  !>
  !> \note However, also there currently only the rock state in the
  !> additional halo is communicated and this is done only once during
  !> initialization. This is for performance reason. The code to
  !> completely replace lbe_halo_exchange() was implemented already
  !> and can be found in lbe_md_fluid.F90.
  integer, save :: halo_extent = 1

  type halo_request
    integer :: extent           !< halo width requested
    character(len=64) :: name   !< short human-readable description
  end type halo_request

  type(halo_request), allocatable, dimension(:) :: halo_requests

  !> # of surrounding lattice nodes for an arbitrary off-lattice position
  integer, parameter :: n_lp=2**3
  !> relative coordinates of these nodes
  integer, save :: lp_sur(3,n_lp)
  !> relative coordinates of the opposing nodes
  integer, save :: opp_lp_sur(3,n_lp)

  !> sizes of total system
  real(kind=rk), save :: tsize(3)

  !> global maximum simulation space boundaries in all directions
  real(kind=rk), save :: maxpos(3)

  !> global minimum simulation space boundaries in all directions
  real(kind=rk), save :: minpos(3)

  !> sizes of my part
  real(kind=rk), save :: chunksize(3)

  !> lo/hi boundaries of my box in each dimension
  real(kind=rk), save :: border(2,3)

  !> \{
  !> \name indices into the radial model array for lbe_init_radial
  integer, parameter :: rad_inner = 1, rad_middle = 2, rad_outer = 3
  !> \}

  !>NEW: For lbe_io.F90:dump_pressure()
  integer, parameter :: droplet = 1, nondroplet = 0

  integer, save :: timesteps_count = 0 !< Global count of timesteps performed

  !> \name common external forcing
  !> \{
  !> has to be set to \c .true. if \c lbe_force should be taken into
  !> account as additional forcing during collision
  logical, save :: use_lbe_force = .false.

  !> array holding the force in the local domain resulting from
  !> different sources, indices are: component, species (0: common
  !> force, 1-3: red, blue, green), x, y, z
  real(kind=rk), save, allocatable, dimension(:,:,:,:,:) :: lbe_force

  !> depth of halo region for \c lbe_force in each direction, enlarged
  !> on demand
  integer, save :: force_halo_extent = 1
  !> \}

  !> \name indices referencing timers
  !> \{
  integer, save :: ti_total,ti_adv,ti_intf,ti_intc,ti_halo,ti_inv,ti_dump
  integer, save :: ti_IBM_init, ti_IBM_intspread, ti_IBM_forces, ti_IBM_MPI, ti_IBM_update, ti_IBM_dump
  !> \}
end module lbe_globals_module
