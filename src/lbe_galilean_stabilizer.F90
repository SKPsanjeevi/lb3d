#include "lbe.h"

!> performs a Galilean transformation of the whole system to remove
!> its center of mass motion
module lbe_galilean_stabilizer_module
    use lbe_analysis_module, only: cached_avg_total_velocity
    use lbe_bdist_module, only: boltz_dist
    use lbe_globals_module
    use lbe_helper_module, only: density,every_n_time_steps,velocity
    use lbe_log_module
#ifdef MD
    use lbe_md_interface_module, only: md_shift_velocities
#endif
    use lbe_parallel_module, only: comm_cart
    use lbe_parms_module, only: arg_input_dfile,arg_input_dfile_set&
         &,galilean_stabilizer,inp_file,nt,nx,ny,nz,tau_r
#ifndef SINGLEFLUID
    use lbe_parms_module, only: tau_b
#endif
#ifdef TRACER
    use lbe_tracer_interface_module, only: tracer_shift_velocities
#endif
    use lbe_types_module, only: lbe_site

    implicit none
    private
    include 'mpif.h'

    public lbe_galilean_stabilizer_init,lbe_galilean_stabilizer_input&
         &,lbe_galilean_stabilizer_run

    !> remove center of mass motion every how many time steps? (0: never)
    integer,save :: n_reset=0

    namelist /lbe_galilean_stabilizer/ n_reset

contains

    !> remove the whole system's center of mass motion
    !>
    !> \param[in,out] N local lattice chunk with full halo
    subroutine galilean_stabilizer_apply(N)
        type(lbe_site),intent(inout) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        real(kind=rk) :: dist(nvecs),rho_r,u_r(3),v_com(3)
        integer x,y,z
#ifndef SINGLEFLUID
        real(kind=rk) :: rho_b,u_b(3)
#endif

        v_com = cached_avg_total_velocity(N)

        do x=1,nx
           do y=1,ny
              do z=1,nz
                 rho_r = density(N(x,y,z)%n_r)
                 u_r = velocity(N(x,y,z)%n_r)-v_com
                 call boltz_dist(u_r(1),u_r(2),u_r(3)&
                      &,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,dist)
                 N(x,y,z)%n_r = dist*rho_r
#ifndef SINGLEFLUID
                 rho_b = density(N(x,y,z)%n_b)
                 u_b = velocity(N(x,y,z)%n_b)-v_com
                 call boltz_dist(u_b(1),u_b(2),u_b(3)&
                      &,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,dist)
                 N(x,y,z)%n_b = dist*rho_b
#endif
              end do
           end do
        end do

#ifdef MD
        call md_shift_velocities(-v_com)
#endif
#ifdef TRACER
        call tracer_shift_velocities(-v_com)
#endif

        write (msgstr,"('Galilean stabilizer at nt= ',I0,' : removed v_com = ( ',3(ES15.8,X),')')") nt,v_com
        call log_msg(msgstr)
    end subroutine galilean_stabilizer_apply

    !> set up Galilean stabilizer
    subroutine lbe_galilean_stabilizer_init
        if (.not.galilean_stabilizer) return

#ifndef NOSURFACTANT
        ! a generalization for more than one fluid component should be
        ! straightforward but is was not necessary so far
        call error('Galilean stabilizer not implemented for surfactant '&
             &//'species---set galilean_stabilizer=.false. or compile '&
             &//'with NOSURFACTANT')
#endif
        ! with tau_[rb]/=1 there would be no complete relaxation after
        ! the collision step, so the non-equilibrium part of the
        ! distributions would need to be conserved somehow while the
        ! current implementation just creates new equilibrium
        ! distributions everywhere. A solution is probably present
        ! already in the Lees-Edwards code and the related literature.
        if (tau_r/=1.0_rk) call error('Galilean stabilizer not implemented for '&
             &//'tau_r/=1---set tau_r=1 or galilean_stabilizer=.false.')
#ifndef SINGLEFLUID
        if (tau_b/=1.0_rk) call error('Galilean stabilizer not implemented for '&
             &//'tau_b/=1---set tau_b=1 or galilean_stabilizer=.false.')
#endif
    end subroutine lbe_galilean_stabilizer_init

    !> read namelist \c /galilean_stabilizer/ and bcast parameters
    subroutine lbe_galilean_stabilizer_input
        integer ierror

        if (.not.galilean_stabilizer) return

        call log_msg_hdr("Reading Galilean stabilizer input")

        if (myrankc==0) then
           open (unit=input_file_unit,file=trim(inp_file),err=100)
           read (unit=input_file_unit,nml=lbe_galilean_stabilizer,err=100)
           close (unit=input_file_unit,err=100)
        end if

        if ( arg_input_dfile_set ) then
          call log_msg("  Getting differential input...")
          open (unit=input_dfile_unit,file=arg_input_dfile,status='UNKNOWN')
          read (unit=input_dfile_unit,nml=lbe_galilean_stabilizer,iostat=ierror)
          if (ierror/=0) call log_msg('    WARNING: Differential namelist '&
               &//'not found or errors encountered.')
          close (unit=input_dfile_unit)
          call log_ws()
        end if

        write (msgstr,"('n_reset = ',I8)") n_reset
        call log_msg(msgstr)
        call log_ws()

        call MPI_Bcast(n_reset,1,MPI_INTEGER,0,comm_cart,ierror)

        return
100     continue
        call error('Error reading input file "'//trim(inp_file)//'"')
    end subroutine lbe_galilean_stabilizer_input

    !> Keeps Galilean stabilizer running. Intended to be called once
    !> per LB time step, after the collision stage.
    !>
    !> \param[in,out] N local lattice chunk with full halo
    subroutine lbe_galilean_stabilizer_run(N)
        type(lbe_site),intent(inout) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)

        if (.not.galilean_stabilizer) return

        if (every_n_time_steps(n_reset)) call galilean_stabilizer_apply(N)
    end subroutine lbe_galilean_stabilizer_run

end module lbe_galilean_stabilizer_module
