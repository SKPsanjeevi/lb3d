#include "lbe.h"

!> Process input for ELEC part of the code.
module lbe_elec_input_module
#ifdef ELEC
  use lbe_elec_globals_module
  use lbe_elec_helper_module
  use lbe_globals_module, only: input_dfile_unit
  use lbe_log_module, only: log_msg, log_ws ! To display input variables without the [EL] tag.
  use lbe_parms_module, only: inp_file, arg_input_dfile_set, arg_input_dfile, kbT
  use mpi

  implicit none

  private

  public elec_read_input

contains

!> Read ELEC input file, calculate derived quantities, display the variables, and broadcast them.
subroutine elec_read_input()
  implicit none

  character(len=256) :: elec_inp_file
  integer :: ioerror, mpierror

  elec_inp_file = trim(inp_file)//'.elec'

  call log_msg_elec_hdr("Reading ELEC input")

  ! Read input file(s).
  if (myrankc .eq. 0) then
    open(unit = elec_input_file_unit, file = trim(elec_inp_file), status = 'UNKNOWN')
    read(unit = elec_input_file_unit, nml = elec_input, iostat = ioerror)
    if (ioerror .ne. 0) then
      call error_elec("Namelist not found or errors encountered.")
    end if
    close(unit = elec_input_file_unit)

    if ( arg_input_dfile_set ) then
      call log_msg_elec("  Getting differential input...")
      open(unit = input_dfile_unit, file = arg_input_dfile, status = 'UNKNOWN')
      read(unit = input_dfile_unit, nml = elec_input, iostat = ioerror)
      if (ioerror .ne. 0) then
        call log_msg_elec("    WARNING: Differential namelist not found or errors encountered.")
      end if
      close(unit = input_dfile_unit)
      call log_ws()
    end if
  end if

  ! Compute derived quantities and report on values.

  ! Charges
  call log_msg("Charges:")
  write(msgstr,"('  ec (fixed)         = ',F16.10)") ec                ; call log_msg(msgstr)
  write(msgstr,"('  rho_init           = <',A,'>')") trim(rho_init)    ; call log_msg(msgstr)
  if ( trim(rho_init) == "uniform" ) then
    write(msgstr,"('   Q_colloid         = ',F16.10)") Q_colloid       ; call log_msg(msgstr)
    write(msgstr,"('   Q_wall            = ',F16.10)") Q_wall          ; call log_msg(msgstr)
    write(msgstr,"('   bjerrum_length    = ',F16.10)") bjerrum_length  ; call log_msg(msgstr)
    write(msgstr,"('   debye_length      = ',F16.10)") debye_length    ; call log_msg(msgstr)
    write(msgstr,"('   acc_neutrality    = ',F16.10)") acc_neutrality  ; call log_msg(msgstr)
  else if ( trim(rho_init) == "capacitor" ) then
    write(msgstr,"('   Q_wall            = ',F16.10)") Q_wall          ; call log_msg(msgstr)
  else if ( trim(rho_init) == "liquidjunction" ) then
    write(msgstr,"('   rho_lj            = ',F16.10)") rho_lj          ; call log_msg(msgstr)
    write(msgstr,"('   delta_rho_lj      = ',F16.10)") delta_rho_lj    ; call log_msg(msgstr)
  else if ( trim(rho_init) == "slipbc" ) then
    write(msgstr,"('   Q_slip            = ',F16.10)") Q_slip          ; call log_msg(msgstr)
    write(msgstr,"('   Q_noslip          = ',F16.10)") Q_noslip        ; call log_msg(msgstr)
  end if
  call log_ws()

  ! Diffusivity
  call log_msg("Diffusivity:")
  ! Calculate derived inverse temperature.
  if ( kbT .le. 0.0_rk ) then
    call error_elec("Temperature kbT needs to be greater than 0.")
  end if
  beta = 1.0_rk / kbT
  D_plus  = D_therm + delta_D_therm
  D_minus = D_therm - delta_D_therm
  write(msgstr,"('  beta (derived)     = ',F16.10)") beta              ; call log_msg(msgstr)
  write(msgstr,"('  D_therm            = ',F16.10)") D_therm           ; call log_msg(msgstr)
  write(msgstr,"('  delta_D_therm      = ',F16.10)") delta_D_therm     ; call log_msg(msgstr)
  write(msgstr,"('  D_plus (derived)   = ',F16.10)") D_plus            ; call log_msg(msgstr)
  write(msgstr,"('  D_minus (derived)  = ',F16.10)") D_minus           ; call log_msg(msgstr)
  call log_ws()

  ! Poisson solver
  call log_msg("Poisson solver:")
  ! Set Poisson solver ID from string value.
  select case( trim(poisson_solver) )
    case ( "sor" )
      poisson_solver_id = poisson_solver_SOR
    case ( "p3m" )
#ifdef P3M
      poisson_solver_id = poisson_solver_p3m
#else
      call error_elec("When using P3M as a Poisson solver, please set -DP3M.")
#endif
    case default
      write(msgstr,"('  Invalid poisson_solver value <',A,'> .')") trim(poisson_solver)
      call error_elec(msgstr)
  end select
  ! Set E solver ID from string value.
  select case( trim(E_solver) )
    case ( "fd" )
      E_solver_id = E_solver_fd
    case ( "p3m" )
#ifdef P3M
      E_solver_id = E_solver_p3m
#else
      call error_elec("When using P3M as a solver for the electric field, please set -DP3M.")
#endif
    case default
      write(msgstr,"('  Invalid E_solver value <',A,'> .')") trim(E_solver)
      call error_elec(msgstr)
  end select
  ! Check for invalid solver combinations.
  if ( ( E_solver_id .eq. E_solver_p3m ) .and. ( poisson_solver_id .ne. poisson_solver_p3m ) ) then
    call error_elec("When using P3M as E solver, one also needs to set P3M as Poisson solver.")
  end if
  ! Calculate the spectral radius of Jacobi iteration for SOR solver.
  radius_SOR = 1.0_rk - ( 0.5_rk*( ( pi / max( real(tnx,kind=rk), max( real(tny,kind=rk), real(tnz,kind=rk) ) ) )**2 ) )
  write(msgstr,"('  poisson_solver     = <',A,'>')") trim(poisson_solver) ; call log_msg(msgstr)
  if ( poisson_solver_id == poisson_solver_SOR ) then
    if ( n_check_SOR .lt. 1 ) then
      call error_elec("When using SOR as solver, n_check_SOR has to be at least 1.")
    end if
    write(msgstr,"('   maxits_SOR        = ',I0)") maxits_SOR          ; call log_msg(msgstr)
    write(msgstr,"('   n_check_SOR       = ',I0)") n_check_SOR         ; call log_msg(msgstr)
    write(msgstr,"('   radius_SOR (der)  = ',F16.10)") radius_SOR      ; call log_msg(msgstr)
    write(msgstr,"('   acc_SOR           = ',F16.10)") acc_SOR         ; call log_msg(msgstr)
  end if
  write(msgstr,"('  acc_fluxes         = ',F16.10)") acc_fluxes        ; call log_msg(msgstr)
  write(msgstr,"('  E_solver           = <',A,'>')") trim(E_solver)    ; call log_msg(msgstr)
  call log_ws()

  ! Permittivity
  call log_msg("Permittivity:")
  ! Calculate derived permittivity.
  if ( eps_global .lt. 0.0_rk ) then
    if ( bjerrum_length .le. 0.0_rk ) then
      call error_elec("When eps_uniform is calculated from the Bjerrum length, bjerrum_length needs to be strictly positive.")
    else
      eps_uniform = beta * ec * ec / ( 4.0_rk * pi * bjerrum_length ) ! Inverting definition of Bjerrum length.
    end if
  else
    eps_uniform = eps_global
  end if
  write(msgstr,"('  local_eps          = ',L1)"    ) local_eps         ; call log_msg(msgstr)
  if ( local_eps ) then
    if ( fluid_on_elec ) then
      if ( ( eps_r .ge. 0.0_rk  ) .and. ( eps_b .ge. 0.0_rk ) .and. ( eps_wall .ge. 0.0_rk) ) then
        eps_avg = 0.5_rk * ( eps_r + eps_b )
        eps_gamma = ( eps_b - eps_r ) / ( eps_b + eps_r )
      else
        call error_elec("Using local_eps combined with fluid_on_elec requires positive eps_r, eps_b and eps_wall.")
      end if
      write(msgstr,"('   eps_r             = ',F16.10)") eps_r           ; call log_msg(msgstr)
#ifndef SINGLEFLUID
      write(msgstr,"('   eps_b             = ',F16.10)") eps_b           ; call log_msg(msgstr)
      write(msgstr,"('   eps_gamma (der)   = ',F16.10)") eps_gamma       ; call log_msg(msgstr)
      write(msgstr,"('   eps_avg (der)     = ',F16.10)") eps_avg         ; call log_msg(msgstr)
#endif
      write(msgstr,"('   eps_wall          = ',F16.10)") eps_wall        ; call log_msg(msgstr)
    else ! no fluid_on_elec
      write(msgstr,"('   eps_init          = <',A,'>')") trim(eps_init)  ; call log_msg(msgstr)
      if ( ( trim(eps_init) == "uniform") .or. ( trim(eps_init) == "capacitor" ) ) then
        write(msgstr,"('   bjerrum_length    = ',F16.10)") bjerrum_length; call log_msg(msgstr)
        write(msgstr,"('   eps_global        = ',F16.10)") eps_global    ; call log_msg(msgstr)
        write(msgstr,"('   eps_uniform (der) = ',F16.10)") eps_uniform   ; call log_msg(msgstr)
      end if
    end if
  else ! no local eps
    write(msgstr,"('   bjerrum_length    = ',F16.10)") bjerrum_length    ; call log_msg(msgstr)
    write(msgstr,"('   eps_global        = ',F16.10)") eps_global        ; call log_msg(msgstr)
    write(msgstr,"('   eps_uniform (der) = ',F16.10)") eps_uniform       ; call log_msg(msgstr)
  end if
  call log_ws()

  ! Geometry
  call log_msg("Geometry:")
  write(msgstr,"('  rock_init          = <',A,'>')") trim(rock_init)   ; call log_msg(msgstr)
  call log_ws()

  ! Boundary conditions
  call log_msg("Boundary conditions:")
  ! Set boundary condition ID from string value.
  select case( trim(boundary_phi) )
    case( "neumann_x" )
      boundary_phi_id = boundary_phi_neumannx
    case( "neumann_z" )
      boundary_phi_id = boundary_phi_neumannz
    case( "periodic" )
      boundary_phi_id = boundary_phi_periodic
    case( "dropz" )
      boundary_phi_id = boundary_phi_dropz
    case default
      write(msgstr,"('  Invalid boundary_phi value <',A,'> .')") trim(boundary_phi)
      call error_elec(msgstr)
  end select
  if ( ( boundary_phi_id == boundary_phi_neumannz ) .or. ( boundary_phi_id == boundary_phi_dropz ) ) then
    if ( poisson_solver_id == poisson_solver_p3m ) then
      call error_elec("Using Neumann or potential drop boundary conditions are not yet supported for the P3M solver.")
    end if
  end if
  write(msgstr,"('  boundary_phi       = <',A,'>')") trim(boundary_phi); call log_msg(msgstr)
  if ( boundary_phi_id == boundary_phi_dropz ) then
    write(msgstr,"('   phi_dropz         = ',F16.10)") phi_dropz; call log_msg(msgstr)
  end if
  call log_ws()

  ! Couplings
  call log_msg("Couplings:")
  write(msgstr,"('  fluid_on_elec      = ',L1)"    ) fluid_on_elec     ; call log_msg(msgstr)
  write(msgstr,"('  elec_on_fluid      = ',L1)"    ) elec_on_fluid     ; call log_msg(msgstr)
#ifndef SINGLEFLUID
  if ( fluid_on_elec ) then
    write(msgstr,"('  delta_mu_plus      = ',F16.10)") delta_mu_plus     ; call log_msg(msgstr)
    write(msgstr,"('  delta_mu_minus     = ',F16.10)") delta_mu_minus    ; call log_msg(msgstr)
  end if
#endif
  if ( elec_on_fluid .and. ( inv_fluid .eq. 17 ) ) then
    write(msgstr,"('  noforce_offset_ux  = ',I0)") noforce_offset_ux ; call log_msg(msgstr)
  end if
  call log_ws()

  ! External electric field
  call log_msg("External electric field:")
  write(msgstr,"('  Ex                 = ',F16.10)") Ex                ; call log_msg(msgstr)
  write(msgstr,"('  Ey                 = ',F16.10)") Ey                ; call log_msg(msgstr)
  write(msgstr,"('  Ez                 = ',F16.10)") Ez                ; call log_msg(msgstr)
  call log_ws()

  ! Dumping
  call log_msg("Dumping:")
  write(msgstr,"('  n_sci_rho_p        = ',I0)") n_sci_rho_p    ; call log_msg(msgstr)
  write(msgstr,"('  n_sci_rho_m        = ',I0)") n_sci_rho_m    ; call log_msg(msgstr)
  write(msgstr,"('  n_sci_phi          = ',I0)") n_sci_phi      ; call log_msg(msgstr)
  write(msgstr,"('  n_sci_eps          = ',I0)") n_sci_eps      ; call log_msg(msgstr)
  write(msgstr,"('  n_sci_E            = ',I0)") n_sci_E        ; call log_msg(msgstr)
  write(msgstr,"('  n_sci_elec         = ',I0)") n_sci_elec     ; call log_msg(msgstr)
  call log_ws()

  ! Debugging
  call log_msg("Debug:")
  write(msgstr,"('  dump_elec_halo     = ',L1)") dump_elec_halo ; call log_msg(msgstr)
  write(msgstr,"('  n_show_eq          = ',I0)") n_show_eq      ; call log_msg(msgstr)
  if ( poisson_solver_id == poisson_solver_SOR ) then
    write(msgstr,"('  n_show_SOR         = ',I0)") n_show_SOR   ; call log_msg(msgstr)
  end if
  write(msgstr,"('  n_dump_eq          = ',I0)") n_dump_eq      ; call log_msg(msgstr)
  write(msgstr,"('  n_eq_noE           = ',I0)") n_eq_noE       ; call log_msg(msgstr)
  write(msgstr,"('  n_eq_E             = ',I0)") n_eq_E         ; call log_msg(msgstr)
  if ( dbg_check_zero_flux ) then
    write(msgstr,"('  dbg_check_zero_flux= ',L1)") dbg_check_zero_flux ; call log_msg(msgstr)
  end if

  ! Broadcast variables. List to be kept in sync with the lbe_elec_globals variable list.

  ! Charges
  call MPI_Bcast(rho_init      , 64, MPI_CHARACTER, 0, comm_cart, mpierror)
  call MPI_Bcast(Q_colloid     , 1 , LBE_REAL,      0, comm_cart, mpierror)
  call MPI_Bcast(Q_wall        , 1 , LBE_REAL,      0, comm_cart, mpierror)
  call MPI_Bcast(Q_slip        , 1 , LBE_REAL,      0, comm_cart, mpierror)
  call MPI_Bcast(Q_noslip      , 1 , LBE_REAL,      0, comm_cart, mpierror)
  call MPI_Bcast(debye_length  , 1 , LBE_REAL,      0, comm_cart, mpierror)
  call MPI_Bcast(acc_neutrality, 1 , LBE_REAL,      0, comm_cart, mpierror)
  call MPI_Bcast(rho_lj        , 1 , LBE_REAL,      0, comm_cart, mpierror)
  call MPI_Bcast(delta_rho_lj  , 1 , LBE_REAL,      0, comm_cart, mpierror)

  ! Diffusivity
  call MPI_Bcast(beta          , 1 , LBE_REAL,      0, comm_cart, mpierror)
  call MPI_Bcast(D_therm       , 1 , LBE_REAL,      0, comm_cart, mpierror)
  call MPI_Bcast(delta_D_therm , 1 , LBE_REAL,      0, comm_cart, mpierror)
  call MPI_Bcast(D_plus        , 1 , LBE_REAL,      0, comm_cart, mpierror)
  call MPI_Bcast(D_minus       , 1 , LBE_REAL,      0, comm_cart, mpierror)

  ! Permittivity
  call MPI_Bcast(bjerrum_length, 1 , LBE_REAL,      0, comm_cart, mpierror)
  call MPI_Bcast(local_eps     , 1 , MPI_LOGICAL,   0, comm_cart, mpierror)
  call MPI_Bcast(eps_init      , 64, MPI_CHARACTER, 0, comm_cart, mpierror)
  call MPI_Bcast(eps_global    , 1 , LBE_REAL,      0, comm_cart, mpierror)
  call MPI_Bcast(eps_uniform   , 1 , LBE_REAL,      0, comm_cart, mpierror)
  call MPI_Bcast(eps_r         , 1 , LBE_REAL,      0, comm_cart, mpierror)
  call MPI_Bcast(eps_b         , 1 , LBE_REAL,      0, comm_cart, mpierror)
  call MPI_Bcast(eps_wall      , 1 , LBE_REAL,      0, comm_cart, mpierror)
  call MPI_Bcast(eps_avg       , 1 , LBE_REAL,      0, comm_cart, mpierror)
  call MPI_Bcast(eps_gamma     , 1 , LBE_REAL,      0, comm_cart, mpierror)

  ! Geometry
  call MPI_Bcast(rock_init     , 64, MPI_CHARACTER, 0, comm_cart, mpierror)

  ! Boundary conditions
  call MPI_Bcast(boundary_phi  , 64, MPI_CHARACTER, 0, comm_cart, mpierror)
  call MPI_Bcast(boundary_phi_id, 1, MPI_INTEGER,   0, comm_cart, mpierror)
  call MPI_Bcast(phi_dropz     , 1 , LBE_REAL,      0, comm_cart, mpierror)

  ! Poisson solver
  call MPI_Bcast(poisson_solver, 64, MPI_CHARACTER, 0, comm_cart, mpierror)
  call MPI_Bcast(poisson_solver_id,1,MPI_INTEGER,   0, comm_cart, mpierror)
  call MPI_Bcast(maxits_SOR    , 1 , MPI_INTEGER,   0, comm_cart, mpierror)
  call MPI_Bcast(n_check_SOR   , 1 , MPI_INTEGER,   0, comm_cart, mpierror)
  call MPI_Bcast(radius_SOR    , 1 , LBE_REAL,      0, comm_cart, mpierror)
  call MPI_Bcast(acc_SOR       , 1 , LBE_REAL,      0, comm_cart, mpierror)
  call MPI_Bcast(acc_fluxes    , 1 , LBE_REAL,      0, comm_cart, mpierror)
  call MPI_Bcast(E_solver      , 64, MPI_CHARACTER, 0, comm_cart, mpierror)
  call MPI_Bcast(E_solver_id   , 1 , MPI_INTEGER,   0, comm_cart, mpierror)

  ! Couplings
  call MPI_Bcast(fluid_on_elec , 1 , MPI_LOGICAL,   0, comm_cart, mpierror)
  call MPI_Bcast(elec_on_fluid , 1 , MPI_LOGICAL,   0, comm_cart, mpierror)
  call MPI_Bcast(delta_mu_plus , 1 , LBE_REAL,      0, comm_cart, mpierror)
  call MPI_Bcast(delta_mu_minus, 1 , LBE_REAL,      0, comm_cart, mpierror)
  call MPI_Bcast(noforce_offset_ux,1,MPI_INTEGER,   0, comm_cart, mpierror)

  ! External electric field
  call MPI_Bcast(Ex            , 1 , LBE_REAL,      0, comm_cart, mpierror)
  call MPI_Bcast(Ey            , 1 , LBE_REAL,      0, comm_cart, mpierror)
  call MPI_Bcast(Ez            , 1 , LBE_REAL,      0, comm_cart, mpierror)

  ! Dumping
  call MPI_Bcast(n_sci_rho_p   , 1 , MPI_INTEGER,   0, comm_cart, mpierror)
  call MPI_Bcast(n_sci_rho_m   , 1 , MPI_INTEGER,   0, comm_cart, mpierror)
  call MPI_Bcast(n_sci_phi     , 1 , MPI_INTEGER,   0, comm_cart, mpierror)
  call MPI_Bcast(n_sci_eps     , 1 , MPI_INTEGER,   0, comm_cart, mpierror)
  call MPI_Bcast(n_sci_E       , 1 , MPI_INTEGER,   0, comm_cart, mpierror)
  call MPI_Bcast(n_sci_elec    , 1 , MPI_INTEGER,   0, comm_cart, mpierror)

  ! Debug
  call MPI_Bcast(dump_elec_halo, 1 , MPI_LOGICAL,   0, comm_cart, mpierror)
  call MPI_Bcast(n_show_eq     , 1 , MPI_INTEGER,   0, comm_cart, mpierror)
  call MPI_Bcast(n_show_sor    , 1 , MPI_INTEGER,   0, comm_cart, mpierror)
  call MPI_Bcast(n_dump_eq     , 1 , MPI_INTEGER,   0, comm_cart, mpierror)
  call MPI_Bcast(n_eq_noE      , 1 , MPI_INTEGER,   0, comm_cart, mpierror)
  call MPI_Bcast(n_eq_E        , 1 , MPI_INTEGER,   0, comm_cart, mpierror)
  call MPI_Bcast(dbg_check_zero_flux, 1 , MPI_INTEGER,   0, comm_cart, mpierror)

  call log_msg_elec_ws("Read and broadcast ELEC parameters.")

  end subroutine elec_read_input

#endif
end module lbe_elec_input_module
