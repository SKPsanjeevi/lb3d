#include "lbe.h"

!> forcing that is sinusoidally modulated in x-direction, producing
!> Kolmogorov flow
module lbe_force_kolmogorov_module
    use lbe_force_interface_module, only: add_force_to_all
    use lbe_globals_module
    use lbe_log_module
    use lbe_parallel_module, only: check_allocate,comm_cart,start
    use lbe_parms_module, only: arg_input_dfile,arg_input_dfile_set,inp_file,nx&
         &,ny,nz

    implicit none
    include 'mpif.h'
    private

    public force_init_kolmogorov,force_input_kolmogorov,force_apply_kolmogorov

    !> buffer to avoid force re-calculation at each time step and lattice site
    real(kind=rk),save,allocatable :: kolmogorov_flow_forces(:)

    real(kind=rk),save :: amplitude=0.0_rk !< maximum absolute value of force

    namelist /force_kolmogorov/ amplitude

contains

    !> apply sinusoidally modulated forces for Kolmogorov flow
    subroutine force_apply_kolmogorov
        integer x,y,z

        do x=1,nx
           do y=1,ny
              do z=1,nz
                 call add_force_to_all(0.0_rk,0.0_rk,kolmogorov_flow_forces(x)&
                      &,x,y,z)
              end do
           end do
        end do
    end subroutine force_apply_kolmogorov

    !> mainly pre-calculate position-dependent force
    subroutine force_init_kolmogorov
        integer stat,x

        use_lbe_force = .true.

        allocate (kolmogorov_flow_forces(1:nx),stat=stat)
        call check_allocate(stat&
             &,'force_init_kolmogorov():kolmogorov_flow_forces')

        do x=1,nx
           kolmogorov_flow_forces(x) &
                &= sin(2.0_rk*pi*(real(start(1)+x-1,kind=rk)-0.5_rk)/tsize(1))&
                &*amplitude
        end do
    end subroutine force_init_kolmogorov

    !> read namelist
    subroutine force_input_kolmogorov
        integer ierror

        call log_msg_hdr("Reading FORCE KOLMOGOROV input")

        if (myrankc==0) then
           open (unit=input_file_unit,file=trim(inp_file),err=100)
           read (unit=input_file_unit,nml=force_kolmogorov,err=100)
           close (unit=input_file_unit,err=100)
        end if

        if ( arg_input_dfile_set ) then
           call log_msg("  Getting differential input...")
           open (unit=input_dfile_unit,file=arg_input_dfile,status='UNKNOWN')
           read (unit=input_dfile_unit,nml=force_kolmogorov,iostat=ierror)
           if (ierror/=0) then
              call log_msg("    WARNING: Differential namelist not found "&
                   &//'or errors encountered.')
           end if
           close (unit=input_dfile_unit)
           call log_ws()
        end if

        write (msgstr,"('amplitude = ',F16.10)") amplitude
        call log_msg(msgstr)
        call log_ws()

        call MPI_Bcast(amplitude,1,LBE_REAL,0,comm_cart,ierror)

        return
100     continue
        call error('Error reading input file "'//trim(inp_file))
    end subroutine force_input_kolmogorov

end module lbe_force_kolmogorov_module
