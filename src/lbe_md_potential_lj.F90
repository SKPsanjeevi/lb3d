#include "lbe.h"

!> Lennard-Jones potential
module lbe_md_potential_lj_module
#ifdef MD

    use lbe_globals_module, only: pi,tsize, input_dfile_unit,myrankc
    use lbe_log_module
    use lbe_md_boundary_condition_module, only: rdstsq
    use lbe_md_fluid_ladd_module, only: particle_lubrication&
         &,rc_lubrication_particle
    use lbe_md_fluid_ladd_parms_module, only: lubrication
    use lbe_md_globals_module
    use lbe_md_helper_module, only: number_density,log_msg_md,error_md,log_msg_md_hdr
    use lbe_parallel_module, only: comm_cart
    use lbe_parms_module, only: inp_file,arg_input_dfile_set,arg_input_dfile

    implicit none
    include 'mpif.h'
    private

    public corrections_lj,force_lj,input_lj,local_energy_lj,local_virial_lj&
         &,setup_lj

    !> \name Lennard-Jones parameters
    !> \{
    real(kind=rk),save :: epsilon=0.01_rk !< energy
    real(kind=rk),save :: sigma=7.0_rk !< range
    !> \}

    real(kind=rk),save :: sigsq         !< sigma squared
    real(kind=rk),save :: epsilon4      !< epsilon*4
    real(kind=rk),save :: epsilon48     !< epsilon*48

    namelist /md_potential_lj/ epsilon,sigma

contains

    !> read lj section from md input file
    subroutine input_lj
        integer ierror

        call log_msg_md_hdr("Reading MD P LJ input")
        if (myrankc.eq.0) then
           open (unit=md_input_file_unit,file=trim(inp_file)//'.md',err=100)
           read (unit=md_input_file_unit,nml=md_potential_lj,err=100)
           close (unit=md_input_file_unit,err=100)
           !call log_msg_md('read /md_potential_lj/ from file "'//trim(inp_file)&
           !     &//'.md"',.false.)
           !write (6,nml=md_potential_lj)
        end if

        if ( arg_input_dfile_set ) then
          call log_msg_md("  Getting differential input...")
          open(UNIT = input_dfile_unit, FILE = arg_input_dfile, STATUS = 'UNKNOWN')
          read(UNIT = input_dfile_unit, NML = md_potential_lj, IOSTAT = ierror)
          if (ierror .ne. 0) then
            call log_msg_md("    WARNING: Differential namelist not found or errors encountered.")
          endif
          close(UNIT = input_dfile_unit)
          call log_ws()
        end if


        write(msgstr,"('epsilon            = ',F16.10)") epsilon
        call log_msg(msgstr)
        write(msgstr,"('sigma              = ',F16.10)") sigma
        call log_msg(msgstr)
        call log_ws()

        call MPI_Bcast(epsilon,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(sigma,1,MPI_REAL8,0,comm_cart,ierror)

        return
100     continue
        call error_md('Error reading md input file "'//trim(inp_file)//'.md"')
    end subroutine input_lj

    subroutine setup_lj
        if (rc<2.0_rk*sigma.and.myrankc==0) call error_md(&
             &'rc is significantly smaller than 2.5*sigma.')
        if (lubrication.and.rc<rc_lubrication_particle) call error_md(&
             &'rc too small! '&
             &//'(rc<rc_lubrication_particle but lubrication==.true.')
        sigsq = sigma*sigma
        epsilon4 = epsilon*4.0_rk
        epsilon48 = epsilon*48.0_rk
    end subroutine setup_lj

    subroutine force_lj
        integer i,j,k,ii
        real(kind=rk) ptmp(3),del(3),rsq,sr2,sr6,tmp

        i = atompnt
        do ii = 1,nlocal
           ptmp(:) = P(i)%x(:)
           do k = nnlist(ii),nnlist(ii+1)-1
              j = nlist(k)
              del(:) = ptmp(:) - P(j)%x(:)
              ! ATTENTION: optimized minimum image criterion - result is
              ! useless if it is greater than cutsq1
              if (abs(del(1)).gt.tsize_mrc(1)) then
                 if (del(1).lt.0.0_rk) then
                    del(1) = del(1) + tsize(1)
                 else
                    del(1) = del(1) - tsize(1)
                 endif
              endif
              if (abs(del(2)).gt.tsize_mrc(2)) then
                 if (del(2).lt.0.0_rk) then
                    del(2) = del(2) + tsize(2)
                 else
                    del(2) = del(2) - tsize(2)
                 endif
              endif
              if (abs(del(3)).gt.tsize_mrc(3)) then
                 if (del(3).lt.0.0_rk) then
                    del(3) = del(3) + tsize(3)
                 else
                    del(3) = del(3) - tsize(3)
                 endif
              endif
              rsq = dot_product(del,del)
              if (rsq.lt.cutsq1) then
                 sr2 = sigsq/rsq
                 sr6 = sr2*sr2*sr2
                 tmp = epsilon48*sr6*(sr6-0.5_rk)/rsq
                 P(i)%f(:) = P(i)%f(:) + del(:)*tmp
                 if (j.le.npmax) then
                    P(j)%f(:) = P(j)%f(:) - del(:)*tmp
                 endif
                 if (lubrication) call particle_lubrication(rsq,del,i,j)
              endif
           enddo
           i = list(i)
        enddo
    end subroutine force_lj

    subroutine local_energy_lj(eng)
        real(kind=rk),intent(out) :: eng
        integer i,j,ii
        real(kind=rk) factor,rsq

        eng = 0.0_rk
        i = atompnt
        do ii = 1,nlocal
           do j = nnlist(ii),nnlist(ii+1)-1
              rsq = rdstsq(P(i),P(nlist(j)))
              if (rsq.lt.cutsq1) then
                 factor = 1.0_rk
                 if (nlist(j).gt.npmax) factor = 0.5_rk
                 eng = eng + phi_lj(rsq)*factor
              endif
           enddo
           i = list(i)
        enddo
    end subroutine local_energy_lj

    subroutine local_virial_lj(virial)
        real(kind=rk),intent(out) :: virial
        integer i,j,ii
        real(kind=rk) factor,rsq

        virial = 0.0_rk
        i = atompnt
        do ii = 1,nlocal
           do j = nnlist(ii),nnlist(ii+1)-1
              rsq = rdstsq(P(i),P(nlist(j)))
              if (rsq.lt.cutsq1) then
                 factor = 1.0_rk
                 if (nlist(j).gt.npmax) factor = 0.5_rk
                 virial = virial + fphi_lj(rsq)*factor
              endif
           enddo
           i = list(i)
        enddo
    end subroutine local_virial_lj

    !> derivative of LJ energy (force)  ...TIMES THE DISTANCE! (fj, 2008-02-05)
    real(kind=rk) function fphi_lj(rsq)
        real(kind=rk),intent(in) :: rsq
        real(kind=rk) sr2,sr6

        sr2 = sigsq/rsq
        sr6 = sr2*sr2*sr2
        fphi_lj = epsilon48*sr6*(sr6-0.5_rk)
    end function fphi_lj

    !> LJ energy for two particles
    !> \param[in] rsq squared distance
    !> \returns LJ energy
    real(kind=rk) function phi_lj(rsq)
        real(kind=rk),intent(in) :: rsq
        real(kind=rk) sr2,sr6,sr12

        sr2 = sigsq/rsq
        sr6 = sr2*sr2*sr2
        sr12 = sr6*sr6
        phi_lj = epsilon4*(sr12-sr6)
    end function phi_lj

    !> long range corrections
    !>
    !> mandatory due to the cut-off potential,
    !> see J.M. Haile: Molecular Dynamics Simulation, ch. 6, (pp. 226-230)
    !>
    !> \param[out] engcorr energy correction
    !> \param[out] prscorr pressure correction
    !>
    !> \todo This was not used/edited for a long time, it might be wrong for
    !> specific simulation setups.
    subroutine corrections_lj(engcorr,prscorr)
        real(kind=rk),intent(out) :: engcorr,prscorr
        real(kind=rk) density

        call number_density(density)

        engcorr = epsilon*8.0_rk*pi*density&
             & * (1.0_rk/(9.0_rk*(rc/sigma)**9) - 1.0_rk/(3.0_rk*(rc/sigma)**3))
        prscorr = epsilon*8.0_rk*pi*density*density&
             & * (4.0_rk/(9.0_rk*(rc/sigma)**9) - 2.0_rk/(3.0_rk*(rc/sigma)**3))
    end subroutine corrections_lj

#endif
end module lbe_md_potential_lj_module
