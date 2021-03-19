#include "lbe.h"

!> simple automatic rescaling of the total mass in the system
!>
!> Intended to be run at specific time intervals to prevent a slow
!> drift in the total mass as it might occur in setups that conserve
!> mass only approximately. If there is a massive drift, probably
!> something is going wrong in the simulations or the drift is desired
!> for some reasons. In such cases, using this module is probably no
!> good idea. This is really meant for cases where the changes would
!> affect physics only on long time scales, so they can be safely
!> extinguished before they do harm.
!>
!> Simple SCSP systems should not depend at all on the actual mass
!> densities since \c N%n_r could be rescaled with a pre-factor without
!> changing physics. This is not the case for special driving
!> mechanisms such as lbe_force_module because of F=m*a, so again, the
!> changes done here should stay small.
!>
!> Even though it is possible to take only a sub-volume of the whole
!> system into account to compute the proper rescaling factor, this
!> factor will be applied to the whole system. Otherwise shock waves
!> could occur at the edges of this sub-volume.
module lbe_mass_scaler_module
    use lbe_globals_module
    use lbe_helper_module, only: density,every_n_time_steps,is_fluid
    use lbe_log_module
#ifdef MD
    use lbe_md_interface_module, only: md_lbe_mass_scaler_apply
#endif
    use lbe_parallel_module, only: comm_cart,start,tnx,tny,tnz
    use lbe_parms_module, only: amass_r,arg_input_dfile,arg_input_dfile_set&
         &,inp_file,mass_scaler,nt,nx,ny,nz
    use lbe_types_module, only: lbe_site

    implicit none
    private
    include 'mpif.h'

    public lbe_mass_scaler_init,lbe_mass_scaler_input,lbe_mass_scaler_run



    !> target value for avg mass of red fluid within \c minx and \c maxx
    real(kind=rk),save :: m_r=1.0_rk
    !> max lattice position where \c m_r should hold (<0: upper system boundary)
    integer,save :: maxx(3)=-1
    !> min lattice position where \c m_r should hold (<0: lower system boundary)
    integer,save :: minx(3)=-1
    integer,save :: n_scale=0  !< time step interval to rescale (0: never)

    integer,save :: lmaxx(3) !< local version of \c maxx
    integer,save :: lminx(3) !< local version of \c minx

    namelist /lbe_mass_scaler/ m_r,maxx,minx,n_scale	
contains

    !> rescale the mass in the system
    !>
    !> \param[in,out] N local lattice chunk including full halo
    subroutine mass_scaler_apply(N)
        type(lbe_site),intent(inout) &
             &:: N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        integer ierror,x,y,z
        real(kind=rk) :: factor,m(2),tm(2)

        ! find the global avg mass density on all fluid sites, use
        ! double also for the number of sites to avoid 32bit integer
        ! wrap around for large systems and to allow for a joint
        ! reduction operation
        m = 0.0_rk
        do x=lminx(1),lmaxx(1)
           do y=lminx(2),lmaxx(2)
              do z=lminx(3),lmaxx(3)
                 if (is_fluid(N(x,y,z)%rock_state)) then
                    m(1) = m(1)+amass_r*density(N(x,y,z)%n_r)
                    m(2) = m(2)+1.0_rk
                 end if
              end do
           end do
        end do
        call MPI_Allreduce(m,tm,2,LBE_REAL,MPI_SUM,MPI_COMM_WORLD,ierror)

        ! rescale densities in the whole system, so that m_r is
        ! achieved as avg mass density within the sub-volume specified
        ! by minx and maxx
        factor = m_r/(tm(1)/tm(2))
        do x=1,nx
           do y=1,ny
              do z=1,nz
                 ! if n_r at non-fluid sites is abused for something
                 ! that should not be rescaled, this functionality
                 ! will be broken here...
                 N(x,y,z)%n_r = factor*N(x,y,z)%n_r
              end do
           end do
        end do

        ! report on what just happened
        write (msgstr,"('mass_scaler: nt= ',I0,' : old avg mass density '"&
             &//",ES15.8,' on ',ES15.8,' fluid sites, rescaled by factor '"&
             &//",ES15.8)") nt,tm(1)/tm(2),tm(2),factor
        call log_msg(msgstr)

#ifdef MD
        ! give also the md modules a chance to rescale things which
        ! should be rescaled in sync with the masses in N
        call md_lbe_mass_scaler_apply(factor)
#endif
    end subroutine mass_scaler_apply

    !> set up everything after the parameters were read
    subroutine lbe_mass_scaler_init()
        if (.not.mass_scaler) return

#ifndef SINGLEFLUID
        call error('mass_scaler_module() not yet implemented '&
             &//'for multi-component case---set mass_scaler=.false. '&
             &//'or compile with SINGLEFLUID')
#endif
        ! clip to reasonable boundaries
        lminx = max(1,minx)
        lmaxx = min((/tnx,tny,tnz/),maxx)
        where (lmaxx<0) lmaxx=(/tnx,tny,tnz/)

        ! convert to local coordinates
        lminx = lminx+1-start
        lmaxx = lmaxx+1-start

        ! clip to local chunk, no mass calculation in the halo
        lminx = max(1,lminx)
        lmaxx = min((/nx,ny,nz/),lmaxx)
    end subroutine lbe_mass_scaler_init

    !> read \c /mass_scaler/ from input file
    subroutine lbe_mass_scaler_input()
        integer ierror

        if (.not.mass_scaler) return

        call log_msg_hdr("Reading mass scaler input")

        if (myrankc==0) then
           open (unit=input_file_unit,file=trim(inp_file),err=100)
           read (unit=input_file_unit,nml=lbe_mass_scaler,err=100)
           close (unit=input_file_unit,err=100)
        end if

        if (arg_input_dfile_set) then
           call log_msg("  Getting differential input...")
           open (unit=input_dfile_unit,file=arg_input_dfile,status='UNKNOWN')
           read (unit=input_dfile_unit,nml=lbe_mass_scaler,iostat=ierror)
           if (ierror/=0) then
              call log_msg("    WARNING: Differential namelist not found "&
                   &//'or errors encountered.')
           end if
           close (unit=input_dfile_unit)
           call log_ws()
        end if

        disabled: if (n_scale==0) then
           write (msgstr&
                &,"('n_scale         = ',I0,' (mass scaler disabled)')") n_scale
           call log_msg(msgstr)
        else disabled
           write (msgstr,"('m_r             = ',F16.10)") m_r
           call log_msg(msgstr)
           write (msgstr,"('minx            = (',2(I0,','),I0,')')") minx
           call log_msg(msgstr)
           write (msgstr,"('maxx            = (',2(I0,','),I0,')')") maxx
           call log_msg(msgstr)
           write (msgstr,"('n_scale         = ',I0)") n_scale
           call log_msg(msgstr)
           call log_ws()
        end if disabled

        call MPI_Bcast(m_r,1,LBE_REAL,0,comm_cart,ierror)
        call MPI_Bcast(minx,3,MPI_INTEGER,0,comm_cart,ierror)
        call MPI_Bcast(maxx,3,MPI_INTEGER,0,comm_cart,ierror)
        call MPI_Bcast(n_scale,1,MPI_INTEGER,0,comm_cart,ierror)

        return
100     continue
        call error('Error reading input file "'//trim(inp_file))
    end subroutine lbe_mass_scaler_input

    !> depending on \c n_scale and the current time step, rescale the
    !> mass in the system
    !>
    !> \param[in,out] N local lattice chunk including full halo
    subroutine lbe_mass_scaler_run(N)
        type(lbe_site),intent(inout) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)

        if (.not.mass_scaler) return

        if (every_n_time_steps(n_scale)) call mass_scaler_apply(N)
    end subroutine lbe_mass_scaler_run

end module lbe_mass_scaler_module
