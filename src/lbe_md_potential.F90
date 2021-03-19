#include "lbe.h"

!> delegates to different particle-particle potentials
module lbe_md_potential_module
#ifdef MD

    use lbe_md_globals_module
    use lbe_md_helper_module, only: count_particles_all,error_md,log_msg_md&
         &,number_density
    use lbe_md_potential_bprh_module, only: force_and_torque_bprh&
         &,input_bprh,local_energy_bprh,local_virial_bprh&
         &,set_growth_factor_bprh,setup_bprh
    use lbe_md_potential_bprhag_module, only: force_and_torque_bprhag&
         &,input_bprhag,local_energy_bprhag,local_virial_bprhag&
         &,setup_bprhag
    use lbe_md_potential_gb_module, only: force_and_torque_gb,input_gb&
         &,local_energy_gb,local_virial_gb,setup_gb
    use lbe_md_potential_hertz_module, only: force_hertz&
         &,input_hertz,input_hertz_softrep,local_energy_hertz,local_energy_hertz_softrep,local_virial_hertz,setup_hertz&
         &,decide_hertz,setup_hertz_softrep,force_hertz_softrep
    use lbe_md_potential_lj_module, only: corrections_lj,force_lj&
         &,input_lj,local_energy_lj,local_virial_lj,setup_lj
    use lbe_md_potential_dlvo_module, only: force_dlvo&
         &,input_dlvo,local_energy_dlvo,local_virial_dlvo,setup_dlvo
     use lbe_md_potential_spring_module, only: force_spring,input_spring,setup_spring, local_energy_spring,local_virial_spring
    use lbe_md_potential_lubrication_mod, only: force_lubrication&
         &,set_growth_factor_lubrication,setup_lubrication
    use lbe_md_fluid_ladd_parms_module, only: lubrication
    use lbe_parallel_module, only: comm_cart

    implicit none
    private

    include 'mpif.h'

    public corrections,energy,force_and_torque,pressure,input_potential&
         &,set_potential_growth_factor,setup_potential

    character(len=32),save,public :: potential='none' !< potential type

contains

    !> read potential section from md input file
    subroutine input_potential
        select case (potential)
        case ('bprh')
           call input_bprh
        case ('bprhag')
           call input_bprhag
        case ('gb')
           call input_gb
        case ('hertz')
           call input_hertz
        case ('hertzandsoftrep')
           call input_hertz_softrep
        case ('lj')
           call input_lj
        case ('dlvo')
           call input_dlvo
         case ('spring')
             call input_spring
        case ('none')
           ! nop
        case default
           call error_md('unknown value: potential="'//potential//'"')
        end select
    end subroutine input_potential

    !> rescales all potential length parameters
    !>
    !> \param[in] f linear length scale factor with respect to
    !> original lengths read from input file
    subroutine set_potential_growth_factor(f)
        real(kind=rk),intent(in) :: f

        select case (potential)
        case ('bprh')
           call set_growth_factor_bprh(f)
        case ('bprhag')
           call error_md('set_growth_factor_bprhag() not implemented yet')
        case ('gb')
           call error_md('set_growth_factor_gb() not implemented yet')
        case ('hertz')
           call error_md('set_growth_factor_hertz() not implemented yet')
        case ('hertzandsoftrep')
           call error_md('set_growth_factor_hertz_softrep() not implemented yet')
        case ('lj')
           call error_md('set_growth_factor_lj() not implemented yet')
        case ('dlvo')
           call error_md('set_growth_factor_dlvo() not implemented yet')
        case ('spring')
             call error_md('set_growth_factor_spring() not implemented yet')
        case ('none')
           if (lubrication) call set_growth_factor_lubrication(f)
        case default
           call error_md('unknown value: potential="'//potential//'"')
        end select
    end subroutine set_potential_growth_factor

    !> initialization
    subroutine setup_potential
        select case (potential)
        case ('bprh')
           call setup_bprh
        case ('bprhag')
           call setup_bprhag
        case ('gb')
           call setup_gb
        case ('hertz')
           call setup_hertz
        case ('hertzandsoftrep')
           call setup_hertz_softrep
        case ('lj')
           call setup_lj
        case ('dlvo')
           call setup_dlvo
        case ('spring')
             call setup_spring
        case ('none')
           if (lubrication) call setup_lubrication
        case default
           call error_md('unknown value: potential="'//potential//'"')
        end select
    end subroutine setup_potential

    !> branch to potentials for force and possibly also torque calculation
    !>
    !> Newton's 3rd law, force/torque components stored for
    !> atom \c i always and for atom \c j IF I own it
    subroutine force_and_torque
        select case (potential)
        case ('bprh')
           call force_and_torque_bprh
        case ('bprhag')
           call force_and_torque_bprhag
        case ('gb')
           call force_and_torque_gb
        case ('hertz')
           call decide_hertz
        case ('hertzandsoftrep')
           call force_hertz_softrep
        case ('lj')
           call force_lj
        case ('dlvo')
           call force_dlvo
        case ('spring')
             call force_spring
        case ('none')
           if (lubrication) call force_lubrication
        case default
           call error_md('unknown value: potential="'//potential//'"')
        end select
    end subroutine force_and_torque

    !> averaged particle-particle potential per particle
    !> \param[out] eng potential energy
    subroutine energy(eng)
        real(kind=rk),intent(out) :: eng
        real(kind=rk) eeng
        integer ierror,n_global

        select case (potential)
        case ('bprh')
           call local_energy_bprh(eng)
        case ('bprhag')
           call local_energy_bprhag(eng)
        case ('gb')
           call local_energy_gb(eng)
        case ('hertz')
           call local_energy_hertz(eng)
        case ('hertzandsoftrep')
           call local_energy_hertz_softrep(eng)
        case ('lj')
           call local_energy_lj(eng)
        case ('dlvo')
           call local_energy_dlvo(eng)
       case ('spring')
               call local_energy_spring(eng) 
       case ('none')
           eng = 0.0_rk
        case default
           call error_md('unknown value: potential="'//potential//'"')
        end select

        call mpi_allreduce(eng,eeng,1,MPI_REAL8,MPI_SUM,comm_cart,ierror)
        call count_particles_all(n_global)
        eng = eeng/n_global
    end subroutine energy

    !> reduced pressure from virial
    !> \param[in] t temperature
    !> \param[out] p pressure
    subroutine pressure(p,t)
        real(kind=rk),intent(in) :: t
        real(kind=rk),intent(out) :: p
        real(kind=rk) vir,vvir,density
        integer ierror,n_global

        call number_density(density)

        select case (potential)
        case ('bprh')
           call local_virial_bprh(vir)
        case ('bprhag')
           call local_virial_bprhag(vir)
        case ('gb')
           call local_virial_gb(vir)
        case ('hertz')
           call local_virial_hertz(vir)
        case ('hertzandsoftrep')
           call local_virial_hertz(vir)!to replace
        case ('lj')
           call local_virial_lj(vir)
        case ('dlvo')
           call local_virial_dlvo(vir)
        case ('spring')
               call local_virial_spring(vir)
        case ('none')
           vir = 0.0
        case default
           call error_md('unknown value: potential="'//potential//'"')
        end select

        call mpi_allreduce(vir,vvir,1,MPI_REAL8,MPI_SUM,comm_cart,ierror)
        call count_particles_all(n_global)
        p = t*density + density/3.0/n_global*vvir
    end subroutine pressure

    !> calculate possible corrections to energy and pressure
    !>
    !> \param[out] engcorr energy correction
    !> \param[out] prscorr pressure correction
    !>
    !> \todo Corrections need to be implemented in
    !> \c lbe_md_potential_gb_module and in
    !> \c lbe_md_potential_bprhag_module. Both potentials reach zero only
    !> asymptotically for large distances.
    subroutine corrections(engcorr,prscorr)
        real(kind=rk),intent(out) :: engcorr
        real(kind=rk),intent(out) :: prscorr

        select case (potential)
        case ('bprh')
           ! This potential is an analytically exact short range
           ! potential, no correction is needed.
           engcorr = 0.0
           prscorr = 0.0
        case ('bprhag')
           ! implement this
           engcorr = 0.0
           prscorr = 0.0
        case ('gb')
           ! implement this
           engcorr = 0.0
           prscorr = 0.0
        case ('hertz')
           engcorr = 0.0
           prscorr = 0.0
        case ('hertzandsoftrep')
           engcorr = 0.0
           prscorr = 0.0
        case ('lj')
           call corrections_lj(engcorr,prscorr)
        case ('dlvo')
           ! implement this
           engcorr = 0.0
           prscorr = 0.0
        case ('spring')
               engcorr = 0.0
              prscorr = 0.0
        case ('none')
           engcorr = 0.0
           prscorr = 0.0
        case default
           call error_md('unknown value: potential="'//potential//'"')
        end select
    end subroutine corrections

#endif
end module lbe_md_potential_module
