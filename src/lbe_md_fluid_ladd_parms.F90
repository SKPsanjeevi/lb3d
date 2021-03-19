#include "lbe.h"

!> input file parameters for \c lbe_md_fluid_ladd_module, moved here
!> to avoid circular dependencies and keep other modules clean
module lbe_md_fluid_ladd_parms_module
#ifdef MD
    use lbe_globals_module, only: rk

    implicit none
    private

    !> radius of a rock site assumed for lubrication
    real(kind=rk),parameter,public :: R_rock=0.5_rk

    !> cutoff surface distance for normal translational lubrication
    !> (Ladd (2001) suggests 2/3 for \f$\tau=1\f$)
    real(kind=rk),save,public :: delta_c=2.0_rk/3.0_rk

    !> cutoff surface distance for tangential translational
    !> lubrication (Nguyen&Ladd (2002) suggest 1/2 for all \f$\tau\f$)
    real(kind=rk),save,public :: delta_c_T=0.5_rk

    !> cutoff surface distance for tangential rotational lubrication
    !> (Nguyen&Ladd (2002) suggest 0.43 for \f$\tau=1\f$)
    real(kind=rk),save,public :: delta_c_R=0.43_rk

    !> switch particle-particle normal lubrication force on or off
    logical,save,public :: lubrication=.false.

    !> lubrication force will not grow beyond the value at this
    !> surface separation, measured in lu, even in case of closer
    !> contact. This affects only the analytical correction as in
    !> Ladd, 2001 for spherical particles. Default is 1% of the value
    !> suggested by Ladd for \c delta_c .
    real(kind=rk),save,public :: lubrication_clip_dist=(2.0_rk/3.0_rk)/100.0_rk

    !> If with \c lubrication==.true. for spherical particles, the
    !> surface separation is smaller than \c lubrication_clip_dist , a
    !> repulsive spring force with a stiffness of \c
    !> lubrication_stiffness starts to act. For stability reasons it
    !> is limited to the value it reaches at zero separation.
    real(kind=rk),save,public :: lubrication_stiffness=0.0_rk

    !> enables analytical per-particle lubrication corrections in
    !> tangential direction for force and torque
    logical,save,public :: lubrication_tangential=.false.

    !> print counters every how many time steps to standard output (0
    !> means never)
    integer,save,public :: n_print_overlap = 10

    !> colour of particle sites
    real(kind=rk),save,public :: particle_colour=0.0_rk

    !> colour of particle sites (target in case of changing pc)
    real(kind=rk),save,public :: particle_colour_target=0.0_rk

    !> when to start changing particle colour
    integer,save,public :: n_particle_colour_change_start = -1

    !> how long to change particle colour
    integer,save,public :: n_particle_colour_change = -1

    !> particle radii orthogonal to rotational symmetry axis
    real(kind=rk),save,public :: R_orth=2.5_rk

    !> particle radius parallel to rotational symmetry axis
    real(kind=rk),save,public :: R_para=2.5_rk

    !> specifies optional polydisperse particle initialization
    character(len=32),save,public :: initial_radii='none'

    !> maximum particle radius, used for polydisperse initialization
    !> via \c initial_radii
    real(kind=rk),save,public :: R_max=-1.0_rk

    !> minimum particle radius, used for polydisperse initialization
    !> via \c initial_radii
    real(kind=rk),save,public :: R_min=-1.0_rk

    !> standard deviation of particle half-axes, used for polydisperse
    !> initialization via \c initial_radii=='normal'
    real(kind=rk),save,public :: stdev_R=0.0_rk

    real(kind=rk),save,public :: inv_delta_c   !< \f$1/\delta_c\f$
    real(kind=rk),save,public :: log_delta_c_T !< \f$\ln\delta_{c_T}\f$
    real(kind=rk),save,public :: log_delta_c_R !< \f$\ln\delta_{c_R}\f$

    !> parameters for Janus particles
    !> switch on Janus particle for the parallel direction
    logical,save,public :: Janus_para=.false.
    real(kind=rk),save,public :: particle_colour_top=0.0_rk
    real(kind=rk),save,public :: particle_colour_bottom=0.0_rk
    real(kind=rk),save,public :: alpha=0.0_rk
    real(kind=rk),save,public :: n_janus=0.0_rk

#endif
end module lbe_md_fluid_ladd_parms_module
