#include "lbe.h"

!> DLVO potential
module lbe_md_potential_dlvo_module
#ifdef MD

    use lbe_globals_module, only: pi,tsize, input_dfile_unit,myrankc
    use lbe_log_module
    use lbe_md_boundary_condition_module, only: rdstsq
    use lbe_md_fluid_ladd_module, only: particle_lubrication&
         &,rc_lubrication_particle
    use lbe_md_fluid_ladd_parms_module, only: lubrication,R_orth,R_para
    use lbe_md_globals_module
    use lbe_md_helper_module, only: number_density,log_msg_md,error_md,log_msg_md_hdr
    use lbe_parallel_module, only: comm_cart
    use lbe_parms_module, only: inp_file,arg_input_dfile_set,arg_input_dfile

    implicit none
    include 'mpif.h'
    private

    public force_dlvo,input_dlvo,local_energy_dlvo,local_virial_dlvo&
         &,setup_dlvo

    !> \name DLVO parameters
    !> \{
    real(kind=rk),save :: ionic_strength=0.0000001_rk !< ionic strength or salt concentration for NaCl
    real(kind=rk),save :: kbT=0.0000001_rk            !< kbT, unsure if it scales 1:1 with thermal fluxes kbT or elec kbT, hence the name
    real(kind=rk),save :: zeta=0.000000001_rk         !< Zeta potential
    real(kind=rk),save :: AH=10.0_rk                  !< Hamacker constant
    real(kind=rk),save :: hertz=10.0_rk               !< Hertz force
    !> \}

    !> \name DLVO constants, not from namelist but hardcoded
    real(kind=rk),save :: er                          !< relative dielectric constant of the solvent
    real(kind=rk),save :: eo                          !< permittivity of the vacuum
    real(kind=rk),save :: echarge                     !< Elementary charge
    real(kind=rk),save :: IonicStrength               !< Ionic strength
    real(kind=rk),save :: zvalency                    !< valency electrons
    real(kind=rk),save :: Khertz
    !Internal variables
    real(kind=rk),save :: bjerrumLength
    real(kind=rk),save :: debyeKappa
    !more variables
    real(kind=rk),save :: rn
    real(kind=rk),save :: vcoul
    real(kind=rk),save :: vwaals
    real(kind=rk),save :: vhertz
!    real(kind=rk),save :: 
    !> \}

    namelist /md_potential_dlvo/ ionic_strength, kbT, zeta, AH, hertz

contains

    !> read dlvo section from md input file
    subroutine input_dlvo
        integer ierror

        call log_msg_md_hdr("Reading MD P DLVO input")
        if (myrankc.eq.0) then
           open (unit=md_input_file_unit,file=trim(inp_file)//'.md',err=100)
           read (unit=md_input_file_unit,nml=md_potential_dlvo,err=100)
           close (unit=md_input_file_unit,err=100)
           !call log_msg_md('read /md_potential_dlvo/ from file "'//trim(inp_file)&
           !     &//'.md"',.false.)
           !write (6,nml=md_potential_dlvo)
        end if

        if ( arg_input_dfile_set ) then
          call log_msg_md("  Getting differential input...")
          open(UNIT = input_dfile_unit, FILE = arg_input_dfile, STATUS = 'UNKNOWN')
          read(UNIT = input_dfile_unit, NML = md_potential_dlvo, IOSTAT = ierror)
          if (ierror .ne. 0) then
            call log_msg_md("    WARNING: Differential namelist not found or errors encountered.")
          endif
          close(UNIT = input_dfile_unit)
          call log_ws()
        end if

        write(msgstr,"('ionic_strength            = ',F16.10)") ionic_strength
        call log_msg(msgstr)
        write(msgstr,"('kbT                       = ',F16.10)") kbT
        call log_msg(msgstr)
        write(msgstr,"('zeta                      = ',F16.10)") zeta
        call log_msg(msgstr)
        write(msgstr,"('AH                        = ',F16.10)") AH
        call log_msg(msgstr)
        write(msgstr,"('hertz                     = ',F16.10)") hertz
        call log_msg(msgstr)
        call log_ws()

        call MPI_Bcast(ionic_strength,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(kbT           ,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(zeta          ,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(AH            ,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(hertz         ,1,MPI_REAL8,0,comm_cart,ierror)

        return
100     continue
        call error_md('Error reading md input file "'//trim(inp_file)//'.md"')
    end subroutine input_dlvo

    subroutine setup_dlvo
        if (rc<4.0_rk*R_orth.and.myrankc==0) call error_md(&
             &'rc is significantly smaller than 4*R_orth.')
        er=81.0_rk                 
        !eo=8.854187817e-12        
        eo=8.854187817            !No idea about the rescaling for now 
        echarge=1.60217662e-19     
        echarge=1.0_rk             
        IonicStrength=0.001_rk     
        zvalency=1.0_rk            
        Khertz=10.0_rk

    end subroutine setup_dlvo

    subroutine force_dlvo
        integer i,j,k,ii
        real(kind=rk) ptmp(3),del(3),rsq,sr2,sr6,tmp

!Continue to edit from here
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
                  tmp=0.0_rk
                  P(j)%f(:) = P(j)%f(:) - del(:)*tmp
                  
!                 sr2 = sigsq/rsq
!                 sr6 = sr2*sr2*sr2
!                 tmp = epsilon48*sr6*(sr6-0.5_rk)/rsq
!                 P(i)%f(:) = P(i)%f(:) + del(:)*tmp
!                 if (j.le.npmax) then
!                    P(j)%f(:) = P(j)%f(:) - del(:)*tmp
!                 endif
                  tmp = fphi_dlvo(rsq)
                  P(i)%f(:) = P(i)%f(:) + del(:)*tmp
                 if (lubrication) call particle_lubrication(rsq,del,i,j)
              endif
           enddo
           i = list(i)
        enddo
    end subroutine force_dlvo

    subroutine local_energy_dlvo(eng)
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
                 eng = eng + phi_dlvo(rsq)*factor
              endif
           enddo
           i = list(i)
        enddo
    end subroutine local_energy_dlvo

    subroutine local_virial_dlvo(virial)
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
                 virial = virial + fphi_dlvo(rsq)*factor
              endif
           enddo
           i = list(i)
        enddo
    end subroutine local_virial_dlvo

    !> derivative of DLVO energy (force)  ...TIMES THE DISTANCE! (dh copy pasting with changes fj, 2008-02-05)
    real(kind=rk) function fphi_dlvo(rsq)
        real(kind=rk),intent(in) :: rsq
        real(kind=rk) sr2,sr6

        !Coulomb calc
        rn=sqrt(rsq)
        bjerrumLength=echarge**2.0_rk/(4.0_rk*pi*eo*er*kbT)
        debyeKappa=sqrt(8.0_rk*pi*bjerrumLength*IonicStrength)
        vcoul=pi*er*eo*( (2+2*R_orth*debyeKappa)/(1+2*R_orth*debyeKappa) * 4*kbT/(zvalency*echarge)*tanh( (zvalency*echarge*Zeta)/(4.0*kbT)  )  )**2  *4 *rn   *exp( -debyeKappa*(rn-2.0*R_orth))*( -debyeKappa*(rn-2.0*R_orth)+1.0)

        !Van der Waals calc
        vwaals=AH*( 2.0/3.0*R_orth**2.0/rn**3.0 - 4.0/3.0*R_orth**2/(rn*(rn**2.0-4.0*R_orth**2.0))  + 2.0/3.0*rn*R_orth**2.0/( (rn**2-4.0*R_orth**2)**2.0 )    )
        !Hertz calc
        if (rn.lt.(2.0_rk*R_orth)) then
            vhertz=-5.0/2.0*Khertz* (2.0*R_orth-rn)**(3.0/2.0)
        else
            vhertz=0.0_rk
        endif


        fphi_dlvo = vcoul+vwaals+vhertz
        fphi_dlvo = fphi_dlvo*rn
    end function fphi_dlvo

    !> DLVO energy for two particles
    !> \param[in] rsq squared distance
    !> \returns DLVO energy
    real(kind=rk) function phi_dlvo(rsq)
        real(kind=rk),intent(in) :: rsq
        real(kind=rk) sr2,sr6,sr12

        !Coulomb calc
        rn=sqrt(rsq)
        bjerrumLength=echarge**2.0_rk/(4.0_rk*pi*eo*er*kbT)
        debyeKappa=sqrt(8.0_rk*pi*bjerrumLength*IonicStrength)
        vcoul=pi*er*eo*( (2+2*R_orth*debyeKappa)/(1+2*R_orth*debyeKappa) * 4*kbT/(zvalency*echarge)*tanh( (zvalency*echarge*Zeta)/(4.0*kbT)  )  )**2  *(2*rn)**2/rn   *exp( -debyeKappa*(rn-2.0*R_orth))

        !Van der Waals calc
        vwaals=-AH/12.0* (  (4.0*R_orth**2.0)/(rn**2.0-4.0*R_orth**2.0) + (4.0*R_orth**2)/(rn**2.0) + 2.0*log( (rn**2.0-4.0*R_orth**2)/rn**2.0 )  ) 
        !Hertz calc
        if (rn.lt.(2.0_rk*R_orth)) then
            vhertz=(Khertz*(2.0*R_orth-rn)**(5.0/2.0))
        else
            vhertz=0.0_rk
        endif

        !sr2 = sigsq/rsq
        !sr6 = sr2*sr2*sr2
        !sr12 = sr6*sr6
        !phi_dlvo = epsilon4*(sr12-sr6)
        phi_dlvo = vcoul+vwaals+vhertz
    end function phi_dlvo

#endif
end module lbe_md_potential_dlvo_module
