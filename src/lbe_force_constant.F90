#include "lbe.h"

!> constant forcing in a predefined, rectangular chunk of the system
!>
!> This is similar to what can be done with g_accn etc. but note that
!> here we have a force instead of an acceleration. Also, this module
!> supports automatic steering towards a desired rate of massflow.
module lbe_force_constant_module
    use lbe_force_interface_module, only: add_force_to_all
    use lbe_globals_module
    use lbe_helper_module, only: every_n_time_steps,massflow
    use lbe_log_module
#ifdef MD
    use lbe_md_fluid_module, only: md_massflow
    use lbe_md_globals_module, only: provide_uid2i
#endif
    use lbe_parallel_module, only: avg_fluid_density_facez,comm_cart&
         &,start,tnx,tny,tnz
    use lbe_parms_module, only: arg_input_dfile,arg_input_dfile_set,inp_file,nt&
         &,nx,ny,nz
    use lbe_types_module, only: lbe_site

    implicit none
    private
    include 'mpif.h'

    public force_init_constant,force_input_constant,force_apply_constant

    real(kind=rk),save :: force(3)=0.0_rk !< force
    integer,save :: maxx(3)=-1 !< maximum lattice position to apply force
    integer,save :: minx(3)=-1 !< minimum lattice position to apply force

    !> \{
    !> \name basic auto-steering to attain a predefined massflow in the system
    !>
    !> The idea is to update \c force*=auto_massfluxz/actual_massflow
    !> each \c n_auto time steps. The user is responsible to set \c
    !> n_auto to a reasonable value or else this will not work at
    !> all. The new value for \c force is printed to stdout, in case
    !> of restoring from a checkpoint, it needs to be copied manually
    !> into the input file replacing the initial value for \c force.

    !> time step interval after which force is updated (0 means never)
    integer,save :: n_auto=0
    !> targeted total massflow in whole plane z==1
    real(kind=rk),save :: auto_massfluxz=0.0_rk
    !> \}

    namelist /force_constant/ auto_massfluxz,force,maxx,minx,n_auto

    integer,save :: lmaxx(3) !< local maximum position for forcing
    integer,save :: lminx(3) !< local minimum position for forcing

contains

    !> update \c force based on actual massflow at global z==1 plane
    !> and \c auto_massfluxz
    !>
    !> \param[in] lbe_N local lattice chunk with halo 1
    !>
    !> \param[in] whole_N local lattice chunk with full halo
    subroutine auto_steer(lbe_N,whole_N)
        type(lbe_site),intent(in) :: lbe_N(0:,0:,0:)
        type(lbe_site),intent(in) :: &
             &whole_N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        real(kind=rk) :: mf(3),mfluxz,old_forcez,rhof_avg
        integer :: x,y,ierror

#ifdef MD
        rhof_avg = avg_fluid_density_facez(lbe_N,1)
#endif
        mf = 0.0_rk
        if (start(3)==1) then
           do x=1,nx
              do y=1,ny
                 mf = mf+&
#ifdef MD
                      &md_massflow(whole_N,x,y,1,rhof_avg)
#else
                      &massflow(whole_N(x,y,1))
#endif
              end do
           end do
        end if
        call MPI_Allreduce(mf(3),mfluxz,1,LBE_REAL,MPI_SUM,MPI_COMM_WORLD&
             &,ierror)

        !> calculate new force assuming a linear relation between
        !> force and mass flow
        old_forcez = force(3)
        force(3) = force(3)*auto_massfluxz/mfluxz

        write (msgstr,"('auto_steer: nt= ',I0,' old_forcez= ',ES15.8,"&
             &//"' mfluxz= ',ES15.8,' force(3)= ',ES15.8)") &
             &nt,old_forcez,mfluxz,force(3)
        call log_msg(msgstr)
    end subroutine auto_steer

    !> apply constant force within the chosen volume
    !>
    !> \param[in] lbe_N local lattice chunk with halo 1
    !>
    !> \param[in] whole_N local lattice chunk with full halo
    ! Comment by Timm (2012-02-28):
    ! Potential bug!
    ! The force should only be applied to physical sites.
    ! If it is also added to the halo, there will be multiple forcings for the sites
    ! which are physical and in the halo of a neighboring process.
    subroutine force_apply_constant(lbe_N,whole_N)
        type(lbe_site),intent(in) :: lbe_N(0:,0:,0:)
        type(lbe_site),intent(in) :: &
             &whole_N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        integer x,y,z

        if (every_n_time_steps(n_auto)) call auto_steer(lbe_N,whole_N)

        do x=lminx(1),lmaxx(1)
           do y=lminx(2),lmaxx(2)
              do z=lminx(3),lmaxx(3)
                 call add_force_to_all(force(1),force(2),force(3),x,y,z)
              end do
           end do
        end do
    end subroutine force_apply_constant

    !> mainly pre-calculate spatial ranges to apply force
    subroutine force_init_constant
        integer stat,x

        use_lbe_force = .true.

        ! clip to reasonable boundaries
        lminx = max(1,minx)
        lmaxx = min((/tnx,tny,tnz/),maxx)
        where (lmaxx<0) lmaxx=(/tnx,tny,tnz/)

        ! convert to local coordinates
        lminx = lminx+1-start
        lmaxx = lmaxx+1-start

        ! clip to local chunk, no forcing in the halo
        lminx = max(1,lminx)
        lmaxx = min((/nx,ny,nz/),lmaxx)
    end subroutine force_init_constant

    !> read namelist and do early steps of initialization
    subroutine force_input_constant
        integer ierror

        call log_msg_hdr("Reading FORCE CONSTANT input")

        if (myrankc==0) then
           open (unit=input_file_unit,file=trim(inp_file),err=100)
           read (unit=input_file_unit,nml=force_constant,err=100)
           close (unit=input_file_unit,err=100)
        end if

        if ( arg_input_dfile_set ) then
           call log_msg("  Getting differential input...")
           open (unit=input_dfile_unit,file=arg_input_dfile,status='UNKNOWN')
           read (unit=input_dfile_unit,nml=force_constant,iostat=ierror)
           if (ierror/=0) then
              call log_msg("    WARNING: Differential namelist not found "&
                   &//'or errors encountered.')
           end if
           close (unit=input_dfile_unit)
           call log_ws()
        end if

        write (msgstr,"('force           = (',2(F16.10,','),F16.10,')')") force
        call log_msg(msgstr)
        write (msgstr,"('minx            = (',2(I0,','),I0,')')") minx
        call log_msg(msgstr)
        write (msgstr,"('maxx            = (',2(I0,','),I0,')')") maxx
        call log_msg(msgstr)
        write (msgstr,"('n_auto          = ',I0)") n_auto
        call log_msg(msgstr)
        if (n_auto>0) then
           write (msgstr,"('auto_massfluxz  = ',F16.10)") auto_massfluxz
        else
           write (msgstr,"('auto_massfluxz  = <ignored>')")
        end if
        call log_msg(msgstr)
        call log_ws()

        call MPI_Bcast(force,3,LBE_REAL,0,comm_cart,ierror)
        call MPI_Bcast(minx,3,MPI_INTEGER,0,comm_cart,ierror)
        call MPI_Bcast(maxx,3,MPI_INTEGER,0,comm_cart,ierror)
        call MPI_Bcast(n_auto,1,MPI_INTEGER,0,comm_cart,ierror)
        call MPI_Bcast(auto_massfluxz,1,LBE_REAL,0,comm_cart,ierror)

#ifdef MD
        ! auto steering relies on md_massflow() which relies on uid2i
        if (n_auto/=0) provide_uid2i = .true.
#endif

        return
100     continue
        call error('Error reading input file "'//trim(inp_file))
    end subroutine force_input_constant

end module lbe_force_constant_module
