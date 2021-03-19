#include "lbe.h"

!> Electrolyte plugin global variables module. All these variables should be public for easy access.
module lbe_elec_globals_module
#ifdef ELEC

  use lbe_globals_module, only: halo_extent, pi, rk, nvecs, nnonrest, myrankc, &
       & g, cx, cy, cz, c, w, force_halo_extent, n_spec
  use lbe_parallel_module, only: nprocs, comm_cart, ccoords, halo_exchange, tnx, tny, tnz
  use lbe_parms_module, only: nt, nx, ny, nz
  use lbe_parms_module, only: inv_fluid, pb, pg ! For electro-osmosis + slip
  use lbe_types_module, only: lbe_site

  implicit none
  public

  ! System
  integer, parameter       :: elec_input_file_unit = 19   !< IO unit for elec input file (must not be the number used for the lbe input file (currently 17) because that file stays open for several function calls).
#ifndef ELEC_NNONREST
  integer, parameter       :: nnn             = 6         !< The first 6 direction vectors are direct nearest neighbours.
#endif
  integer, parameter       :: rk_eps          = rk / 2    !< Real kind to calculate machine epsilon from.

  real(kind=rk), parameter :: meps            = epsilon(1.0_rk_eps) !< Machine epsilon.

  ! Charges
  real(kind=rk), parameter :: ec              = 1.0_rk    !< Elementary unit of charge.
  character(len=64), save  :: rho_init        = "none"    !< Initial condition for the charge densities.
  real(kind=rk), save      :: Q_colloid       = 0.0_rk    !< Charge per colloid (requires MD Ladd particles).
  real(kind=rk), save      :: Q_wall          = 0.0_rk    !< Total charge to be distributed over wall sites.
  real(kind=rk), save      :: Q_slip          = 0.0_rk    !< Charge density to be placed on slip stripes (matching inv_fluid = 17).
  real(kind=rk), save      :: Q_noslip        = 0.0_rk    !< Charge density to be placed on no-slip stripes (matching inv_fluid = 17).
  real(kind=rk), save      :: debye_length    = 0.0_rk    !< Debye screening length.
  real(kind=rk), save      :: acc_neutrality  = meps      !< Requested accuracy of neutrality check for uniform rho.
  real(kind=rk), save      :: rho_lj          = 0.0_rk    !< Charge density when using rho_init = 'liquidjunction'.
  real(kind=rk), save      :: delta_rho_lj    = 0.0_rk    !< Difference between charge densities for positive and negative species when using rho_init = 'liquidjunction'.

  ! Diffusivity
  real(kind=rk), save      :: beta                        !< Derived inverse kbT.
  real(kind=rk), save      :: D_therm         = 0.0_rk    !< Thermal diffusion coefficient.
  real(kind=rk), save      :: delta_D_therm   = 0.0_rk    !< Difference between thermal diffusion 
  real(kind=rk), save      :: D_plus                      !< Derived thermal diffusion coefficient for positive ion species.
  real(kind=rk), save      :: D_minus                     !< Derived thermal diffusion coefficient for negative ion species.

  ! Poisson solver
  character(len=64), save  :: poisson_solver  = "sor"     !< Algorithm used to solve the Poisson equation.
  integer, save            :: poisson_solver_id           !< Derived integer value to store the Poisson solver type.
  integer, parameter       :: poisson_solver_SOR = 1      !< SOR solver id.
  integer, parameter       :: poisson_solver_p3m = 2      !< P3M solver id.
  integer, save            :: maxits_SOR      = 10000     !< Maximum number of iterations of the SOR.
  integer, save            :: n_check_SOR     = 1         !< Iteration interval between MPI_Allreduce calls to determine convergence in SOR solver.
  real(kind=rk), save      :: radius_SOR                  !< Derived spectral radius of Jacobi iteration for SOR solver.
  real(kind=rk), save      :: acc_SOR         = meps      !< Requested accuracy of the SOR solver.
  real(kind=rk), save      :: acc_fluxes      = meps      !< Requested accuracy of the link fluxes.
  character(len=64), save  :: E_solver        = "fd"      !< Method used to calculate the electric field.
  integer, save            :: E_solver_id                 !< Derived integer value to store the Poisson solver type.
  integer, parameter       :: E_solver_fd  = 1            !< Finite difference solver id.
  integer, parameter       :: E_solver_p3m = 2            !< P3M solver id.

  ! Permittivity
  real(kind=rk), save      :: bjerrum_length  = 0.0_rk    !< Bjerrum length.
  logical, save            :: local_eps       = .false.   !< Use local dielectric constant in calculations.
  character(len=64), save  :: eps_init        = "uniform" !< Initial condition for the local dielectric constant.
  real(kind=rk), save      :: eps_global      = 0.0_rk    !< Dielectric constant to be set for uniform epsilon.
  real(kind=rk), save      :: eps_uniform                 !< Derived epsilon, calculated from bjerrum_length if eps_uniform < 0.
  real(kind=rk), save      :: eps_r           = 0.0_rk    !< Dielectric constant for red fluid.
  real(kind=rk), save      :: eps_b           = 0.0_rk    !< Dielectric constant for blue fluid.
  real(kind=rk), save      :: eps_wall        = 0.0_rk    !< Dielectric constant for wall sites.
  real(kind=rk), save      :: eps_avg                     !< Derived average epsilon of red and blue fluid.
  real(kind=rk), save      :: eps_gamma                   !< Derived dielectric contrast.

  ! Geometry
  character(len=64), save  :: rock_init       = "none"    !< Initial condition for the rock.

  ! Boundary conditions
  character(len=64), save  :: boundary_phi    = "periodic"!< Boundary condition for the phi field.
  integer, save            :: boundary_phi_id             !< Derived integer value to store the boundary condition of the phi field.
  integer, parameter       :: boundary_phi_periodic = 1   !< Periodic boundary for phi field.
  integer, parameter       :: boundary_phi_neumannx = 2   !< Neumann boundary condition in x-direction for phi field.
  integer, parameter       :: boundary_phi_neumannz = 3   !< Neumann boundary condition in z-direction for phi field.
  integer, parameter       :: boundary_phi_dropz    = 4   !< Potential drop in z-direction.
  real(kind=rk), save      :: phi_dropz       = 0.0_rk    !< Value of the potential drop.

  ! Couplings
  logical, save            :: fluid_on_elec   = .true.    !< Enable charge advection.
  logical, save            :: elec_on_fluid   = .true.    !< Enable fluid forcing.
  real(kind=rk), save      :: delta_mu_plus   = 0.0_rk    !< Solvation free energy difference between the red and blue fluids for the positive ion species.
  real(kind=rk), save      :: delta_mu_minus  = 0.0_rk    !< Solvation free energy difference between the red and blue fluids for the negative ion species.
  integer, save            :: noforce_offset_ux = 0       !< Increased region of no force in upper region - used with inv_fluid = 17.

  ! External electric field
  real(kind=rk), save      :: Ex              = 0.0_rk    !< x-component of the external electric field.
  real(kind=rk), save      :: Ey              = 0.0_rk    !< y-component of the external electric field.
  real(kind=rk), save      :: Ez              = 0.0_rk    !< z-component of the external electric field.

  ! Dumping
  integer, save :: n_sci_rho_p = 0          !< Length of the interval between dumps of the positive charge scalar field.
  integer, save :: n_sci_rho_m = 0          !< Length of the interval between dumps of the negative charge scalar field.
  integer, save :: n_sci_phi   = 0          !< Length of the interval between dumps of the electric potential scalar field.
  integer, save :: n_sci_eps   = 0          !< Length of the interval between dumps of the local permittivity scalar field.
  integer, save :: n_sci_E     = 0          !< Length of the interval between dumps of the electric field. This will be dumped as three scalar fields Ex, Ey, and Ez.
  integer, save :: n_sci_elec  = 0          !< Write ASCII data of rho_p, rho_m, phi, eps, E(:) and rock_state at all positions. This writes one file per process.

  ! Debugging
  logical, save :: dump_elec_halo = .false. !< Include halo sites in ascii dump.
  integer, save :: n_show_eq   = 0          !< Length of the interval (in iteration loops) during charge equilibration between showing the number of performed iterations.
  integer, save :: n_dump_eq   = 0          !< Length of the interval (in iteration loops) during charge equilibration between dumping a full ASCII file.
  integer, save :: n_show_SOR  = 0          !< Length of the interval (in iteration loops) during SOR between showing the number of performed iterations.
  ! If these are set to >= 0, use this many iterations during equilibration, instead of using convergence.
  integer, save :: n_eq_noE    = -1         !< Number of iterations to perform during equilibration with the external electric field disabled. If this is set to a negative number, the convergence criterion is used instead.
  integer, save :: n_eq_E      = -1         !< Number of iterations to perform during equilibration with the external electric field enabled. If this is set to a negative number, the convergence criterion is used instead.
  logical, save :: dbg_check_zero_flux = .false. !< Enable check of zero total flux.

  ! Miscellaneous
  integer, save :: n_colloids = 0           !< Will hold the number of MD particles present in the system.
  integer, save :: ti_elec_calc_fluxes, ti_elec_move_charges, ti_elec_calc_E

  real(kind=rk), dimension(nvecs), parameter :: flux_link_weight = (/ 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, &
       sqrt(2.0_rk), sqrt(2.0_rk), sqrt(2.0_rk), sqrt(2.0_rk), sqrt(2.0_rk), sqrt(2.0_rk), &
       sqrt(2.0_rk), sqrt(2.0_rk), sqrt(2.0_rk), sqrt(2.0_rk), sqrt(2.0_rk), sqrt(2.0_rk), &
       0.0_rk /)                            !< Weights of the flux links.

  namelist /elec_input/ debye_length, bjerrum_length, D_therm, delta_D_therm, acc_SOR, &
       acc_fluxes, Ex, Ey, Ez, Q_colloid, Q_wall, Q_slip, Q_noslip, rho_lj, delta_rho_lj, maxits_SOR, acc_neutrality,&
       local_eps, eps_init, eps_global, eps_r, eps_b, eps_wall, rho_init, rock_init, boundary_phi,&
       phi_dropz, poisson_solver, n_check_SOR, E_solver,&
       fluid_on_elec, elec_on_fluid, delta_mu_plus, delta_mu_minus, noforce_offset_ux,&
       n_sci_rho_p, n_sci_rho_m, n_sci_phi, n_sci_eps, n_sci_E,&
       n_sci_elec, dump_elec_halo, n_show_eq, n_dump_eq, n_show_SOR, n_eq_noE, n_eq_E, dbg_check_zero_flux

#endif
end module lbe_elec_globals_module
