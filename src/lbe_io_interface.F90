#include "lbe.h"

module lbe_io_interface_module
    use lbe_bc_module, only: lbe_bc_halo_exchange,lbe_bc_le_halo_exchange
    use lbe_globals_module, only: halo_extent, myrankc
    use lbe_helper_module, only: check_dump_now, is_restoring, request_halo_extent
    use lbe_io_arrstats_module, only: dump_arrstats,lbe_io_setup_arrstats
#ifdef USEHDF
    use lbe_io_hdf5_module, only: lbe_init_metadata_hdf5
#endif
    use lbe_io_module, only: dump_fluxz,dump_profile
    use lbe_io_stress_module, only: setup_stress
    use lbe_log_module
    use lbe_parms_module, only: boundary_cond,inv_fluid,n_restore&
         &,n_sci_arrstats,n_sci_fluxz,n_sci_massfluxz,n_sci_profile&
         &,nt,sci_arrstats,sci_fluxz,sci_massfluxz,sci_profile&
         &,sci_stress,n_sci_stress
    use lbe_types_module, only: lbe_site
    implicit none
    private
    public dump_parallel,lbe_io_halo_extent,lbe_io_init
contains

    !> Output functions that are parallel by themselves should be
    !> placed here and not in \c dump_data().
    !>
    !> \note The additional halo exchange is probably unnecessary in most
    !> cases.
    subroutine dump_parallel(lbe_N,whole_N)
        type(lbe_site),intent(inout) :: lbe_N(0:,0:,0:)
        type(lbe_site),intent(inout) :: &
             &whole_N(1-halo_extent:,1-halo_extent:,1-halo_extent:)

        if (check_dump_now(sci_arrstats,n_sci_arrstats)) then
           ! It does not make sense to do this outside the time loop.
           if ( .not. ( nt==0 .or. (is_restoring() .and. nt==n_restore) ) ) then
              call log_msg("Sampling ARRSTATS...")
              call dump_arrstats(lbe_N)
           end if
        end if
        if (check_dump_now(sci_fluxz, n_sci_fluxz) ) then
           call log_msg("Dumping FLUXZ...")
           call dump_fluxz(lbe_N,.false.)
        endif
        if (check_dump_now(sci_massfluxz, n_sci_massfluxz) ) then
           call log_msg("Dumping MASSFLUXZ...")
           call dump_fluxz(lbe_N,.true.)
        endif
        if (check_dump_now(sci_profile, n_sci_profile) ) then
           if ((boundary_cond/=0).and.(inv_fluid==5.or.inv_fluid==6)) then
              call lbe_bc_le_halo_exchange(lbe_N,whole_N)
           else
              call lbe_bc_halo_exchange(whole_N)
           endif
           call log_msg("Dumping PROFILE...")
           call dump_profile(whole_N)
        endif
    end subroutine dump_parallel

    !> increase \c halo_extent if required by some of the IO routines
    subroutine lbe_io_halo_extent
        if (sci_profile) then
           ! this is necessary because of
           ! fluid_velocity_and_density_and_site_occupation() which is
           ! called by dump_profile() and in turn calls
           ! md_calculate_sc_forces()
          call request_halo_extent(2,"sci_profile")
        end if
    end subroutine lbe_io_halo_extent

    subroutine lbe_io_init
#ifdef USEHDF
      if (myrankc==0) call lbe_init_metadata_hdf5()
#endif
      call lbe_io_setup_arrstats
      call setup_stress()
    end subroutine lbe_io_init

end module lbe_io_interface_module

