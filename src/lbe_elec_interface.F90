#include "lbe.h"

!> This module exposes all ELEC functions used by external code (i.e. basic LB3D).
module lbe_elec_interface_module
#ifdef ELEC

  use lbe_elec_timestep_module, only: elec_timestep
  use lbe_elec_forces_module,   only: elec_add_forces, elec_force_reset
  use lbe_elec_input_module,    only: elec_read_input
  use lbe_elec_init_module,     only: elec_init_system, elec_restore_init_system, &
       & elec_initial_equilibration, elec_shutdown
  use lbe_elec_output_module,   only: elec_postprocess
  use lbe_elec_parallel_module, only: elec_setup_parallel, elec_halo
  use lbe_elec_timer_module,    only: elec_init_timers, elec_report_timers, elec_reset_timers

  implicit none
  public

#endif
end module lbe_elec_interface_module

