#include "lbe.h"

!> implementation of Hertz potential
module lbe_md_potential_hertz_module
#ifdef MD

    use lbe_globals_module, only: pi,tsize, input_dfile_unit,myrankc
    use lbe_helper_module, only: cross_product
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

    public force_hertz,force_hertz_softrep,input_hertz,input_hertz_softrep,local_energy_hertz,local_energy_hertz_softrep,local_virial_hertz&
         &,setup_hertz,setup_hertz_softrep,decide_hertz

    ! hertz parameters
    real(kind=rk),save :: K_hertz=0.01_rk
    real(kind=rk),save :: d_hertz=10.0_rk

    ! hertz_soft parameters
    real(kind=rk),save :: K_hsr=0.01_rk
    real(kind=rk),save :: d_hsr=10.0_rk
    real(kind=rk),save :: epsilon_hsr=1.0_rk
    real(kind=rk),save :: sigma_hsr=10.0_rk
    real(kind=rk),save :: r12cut_hsr=20.0_rk

    real(kind=rk),save :: dsq         ! d squared
    real(kind=rk),save :: dsq_hsr         ! d squared
    real(kind=rk),save :: sigsq_hsr         !< sigma squared
    real(kind=rk),save :: rcsq_hsr         ! squared cutoff of r-12 potential, defined via rc
    real(kind=rk),save :: epsilon48_hsr     ! epsilon*48
    real(kind=rk),save :: epsilon96_hsr     ! epsilon*96
    real(kind=rk),save :: epsilon4_hsr     ! epsilon*4
    real(kind=rk),save :: sr2_hsr,sr6_hsr

#ifdef NEWELLIPSOIDMETHOD
    real(kind=rk),dimension(3,3) :: EMM1
    real(kind=rk),dimension(3,3) :: EMM2
#endif

    namelist /md_potential_hertz/ K_hertz,d_hertz
    namelist /md_potential_hertz_softrep/ K_hsr,d_hsr,epsilon_hsr,sigma_hsr,r12cut_hsr

contains

    !> read hertz section from md input file
    subroutine input_hertz
        integer ierror

        call log_msg_md_hdr("Reading MD P Hertz input")
        use_rotation = .true.
        calculate_orientations = .true.

        if (myrankc.eq.0) then
           open (unit=md_input_file_unit,file=trim(inp_file)//'.md',err=100)
           read (unit=md_input_file_unit,nml=md_potential_hertz,err=100)
           close (unit=md_input_file_unit,err=100)
           !call log_msg_md('read /md_potential_hertz/ from file "'//trim(inp_file)//'.md"')
           !write (6,nml=md_potential_hertz)
        end if

        if ( arg_input_dfile_set ) then
          call log_msg_md("  Getting differential input...")
          open(UNIT = input_dfile_unit, FILE = arg_input_dfile, STATUS = 'UNKNOWN')
          read(UNIT = input_dfile_unit, NML = md_potential_hertz, IOSTAT = ierror)
          if (ierror .ne. 0) then
            call log_msg_md("    WARNING: Differential namelist not found or errors encountered.")
          endif
          close(UNIT = input_dfile_unit)
          call log_ws()
        end if

        write(msgstr,"('K_hertz            = ',F16.10)") K_hertz
        call log_msg(msgstr)
        write(msgstr,"('d_hertz            = ',F16.10)") d_hertz
        call log_msg(msgstr)
        call log_ws()

        call MPI_Bcast(K_hertz,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(d_hertz,1,MPI_REAL8,0,comm_cart,ierror)

        return
100     continue
        call error_md('Error reading md input file "'//trim(inp_file)//'.md"')
    end subroutine input_hertz

    !> read hertz_softrep section from md input file
    subroutine input_hertz_softrep
        integer ierror

        call log_msg_md_hdr("Reading MD P Hertz and soft rep input")
        use_rotation = .true.
        calculate_orientations = .true.

        if (myrankc.eq.0) then
           open (unit=md_input_file_unit,file=trim(inp_file)//'.md',err=100)
           read (unit=md_input_file_unit,nml=md_potential_hertz_softrep,err=100)
           close (unit=md_input_file_unit,err=100)
        end if

        if ( arg_input_dfile_set ) then
          call log_msg_md("  Getting differential input...")
          open(UNIT = input_dfile_unit, FILE = arg_input_dfile, STATUS = 'UNKNOWN')
          read(UNIT = input_dfile_unit, NML = md_potential_hertz_softrep, IOSTAT = ierror)
          if (ierror .ne. 0) then
            call log_msg_md("    WARNING: Differential namelist not found or errors encountered.")
          endif
          close(UNIT = input_dfile_unit)
          call log_ws()
        end if

        write(msgstr,"('K_hsr            = ',F16.10)") K_hsr
        call log_msg(msgstr)
        write(msgstr,"('d_hsr            = ',F16.10)") d_hsr
        call log_msg(msgstr)
        write(msgstr,"('epsilon_hsr            = ',F16.10)") epsilon_hsr
        call log_msg(msgstr)
        write(msgstr,"('sigma_hsr              = ',F16.10)") sigma_hsr
        call log_msg(msgstr)
        write(msgstr,"('r12cut_hsr              = ',F16.10)") r12cut_hsr
        call log_msg(msgstr)
        call log_ws()

        call MPI_Bcast(K_hsr,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(d_hsr,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(epsilon_hsr,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(sigma_hsr,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(r12cut_hsr,1,MPI_REAL8,0,comm_cart,ierror)

        return
100     continue
        call error_md('Error reading md input file "'//trim(inp_file)//'.md"')
    end subroutine input_hertz_softrep

    subroutine setup_hertz
        !if (polydispersity) call error_md("potential='hertz' does not yet "&
        !     &//"support polydisperse particles---disable polydispersity "&
        !     &//"or select potential/='hertz'!")
        if (rc<d_hertz.and.myrankc==0) then
             print *,'rc',rc,'d_hertz',d_hertz
             call error_md(&
             &'rc is smaller than d_hertz.')
        end if
        dsq=d_hertz**2
    end subroutine setup_hertz

    subroutine setup_hertz_softrep
        if (polydispersity) call error_md("potential='hertzandsoftrep' "&
             &//"does not yet support polydisperse particles---disable "&
             &//"polydispersity or select potential/='hertzandsoftrep'!")
        if (rc<d_hsr.and.myrankc==0) then
             print *,'rc',rc,'d_hsr',d_hsr
             call error_md(&
             &'rc is smaller than d_hsr.')
        end if
        dsq_hsr=d_hsr**2
        sigsq_hsr = sigma_hsr*sigma_hsr
	rcsq_hsr = r12cut_hsr**2
        epsilon4_hsr = epsilon_hsr*4.0_rk
        epsilon48_hsr = epsilon_hsr*48.0_rk
        epsilon96_hsr = epsilon_hsr*96.0_rk
    end subroutine setup_hertz_softrep

    subroutine force_hertz
        integer i,j,k,ii
        real(kind=rk) ptmp(3),del(3),rsq,tmp,r

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
                if (rsq.lt.dsq) then
                  r=rsq**0.5_rk
                  tmp = 2.5_rk*K_hertz*(d_hertz-r)**1.5_rk/r
                  P(i)%f(:) = P(i)%f(:) + del(:)*tmp
                  if (j.le.npmax) then
                    P(j)%f(:) = P(j)%f(:) - del(:)*tmp
                  endif
                endif
                if (lubrication) call particle_lubrication(rsq,del,i,j)
              endif
           enddo
           i = list(i)
        enddo
    end subroutine force_hertz

    subroutine force_hertz_softrep
        integer i,j,k,ii
        real(kind=rk) ptmp(3),del(3),rsq,tmp,r

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
		if (rsq.lt.rcsq_hsr) then 
		  sr2_hsr = sigsq_hsr/rsq
                  sr6_hsr = sr2_hsr*sr2_hsr*sr2_hsr
                  tmp = tmp +epsilon48_hsr*sr6_hsr*sr6_hsr/rsq
		endif
                if (rsq.lt.dsq_hsr) then
                  r=rsq**0.5_rk
                  tmp = tmp + 2.5_rk*K_hsr*(d_hsr-r)**1.5_rk/r
                  endif
                  P(i)%f(:) = P(i)%f(:) + del(:)*tmp
                  if (j.le.npmax) then
                    P(j)%f(:) = P(j)%f(:) - del(:)*tmp !why a subtraction of what was just added?
                endif
                if (lubrication) call particle_lubrication(rsq,del,i,j)
              endif
           enddo
           i = list(i)
        enddo
    end subroutine force_hertz_softrep

    subroutine force_and_torque_hertz_ellipsoid
        !> evtl noch input_hertz_ellipsoid setup_hertz_ellipsoid
        !> und schreiben und entsprechend eintragen
        !> evtl analog zu force_hertz radius unabhaenig von
        !> anderen ww waehlen und so das ganze unabhaenhig von zB
        !> lubrication oder ladd zu machen
        integer i,j,k,ii
        real(kind=rk) ptmp(3),del(3),rsq,r
        real(kind=rk) :: e,omxoioj,opxoioj,xoioj,rijoiprijoj,rijoimrijoj,sqrt1&
                             &,sqrt2,fracp,fracm,chi,sigma,sigmasq,rijoi&
                             &,chi_half,rijoj,bracket(3),inv_s,x_2s2&
                             &,r_X2,grad_ri_r(3),s2,eps,arij,urij(3),r2&
                             &,fracm_arij,fracp_arij,phit,dphitdr,f(3)&
                             &,ti(3),tj(3)

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

              ! r_{ij}^2
              rsq = dot_product(del,del)

              if (rsq.lt.cutsq1) then
                  ! noch abklaeren, ob diese Zeile ersetzt wird
                  ! ODER DOCH WIEDER ACTIVIERT, DA rc > R_o bzw R_p
                ! Berechnung von rijoiprijoj schon hier
                !d_sqrt2 = 1.0_rk&
                       !&-(rijoiprijoj*fracp+rijoimrijoj*fracm)*chi_half/rij2
                !d_hertz_sigma_sq = ????/d_sqrt2
                !!d_hertz_sigma_sq=d_hertz_sigma**2

                chi = (R_para**2-R_orth**2)&
                    &/(R_para**2+R_orth**2)
                chi_half = 0.5_rk*chi

                ! dot products of rij(del) and o(i|j)
                rijoi = dot_product(del,P(i)%o)
                rijoj = dot_product(del,P(j)%o)

                ! \chi(\mathbf{o}_i\mathbf{o}_j)
                xoioj = chi*dot_product(P(i)%o,P(j)%o)

                ! \mathbf{r}_{ij}\mathbf{o}_i\pm\mathbf{r}_{ij}\mathbf{o}_j)
                ! ----------------------------------------------------------
                !              1\pm\chi\mathbf{o}_i\mathbf{o}_j
                rijoiprijoj = rijoi+rijoj
                rijoimrijoj = rijoi-rijoj
                omxoioj = 1.0_rk-xoioj
                opxoioj = 1.0_rk+xoioj
                fracp = rijoiprijoj/opxoioj
                fracm = rijoimrijoj/omxoioj

                !sigma = 2.0_rk*sigma_orth
                sigma = 2.0_rk*R_orth
                sigmasq = sigma**2

                !x_2s2 = 0.5_rk*chi/sigmasq
                !inv_s = 1.0_rk/sigma
 
                ! \sqrt{1-\frac{\chi}{2}\left[\ldots\right]}^2
                sqrt2 = 1.0_rk&
                     &-(rijoiprijoj*fracp+rijoimrijoj*fracm)*chi_half/rsq

                ! (angle dependent range parameter \sigma)**2
                s2 = sigmasq/sqrt2

                ! (normalized input for radial symmetric potential \Phi(r))**2
                r2 = rsq/s2

                !if (rsq.lt.d_hertz_sigma_sq) then
                    ! alternativ auch r(=rij/sigma).lt.1 ????
                    ! dann sigma Berechnung vorziehen
                    ! altes "r" umbenennen !!!!
                cutoff: if (r2<1.0_rk) then
                    ! r2=r*r, neues r!!!!!!!!!!!!
                    ! rsq neues r2


                    ! \sqrt{1-\frac{\chi}{2}\left[\ldots\right]}
                    sqrt1 = sqrt(sqrt2)

                    ! eps=(2*R_o)^{\frac{5}{2}}K_H
                    eps = (sigma**2.5)*K_hertz

                    ! angle dependent strength parameter \epsilon
                    e = eps/sqrt(omxoioj*opxoioj)

                    ! |\mathbf{r}_{ij}|
                    arij = sqrt(rsq)

                    ! unit vector in direction of rij
                    urij(:) = del(:)/arij

                    ! frac(p|m) with rij replaced by urij
                    fracm_arij = fracm/arij
                    fracp_arij = fracp/arij

                    ! normalized input for radial symmetric potential \Phi(r)
                    r = arij*sqrt1/sigma
                    !r = sqrt(r2)
                    ! aequivalente gleichungen???? wennja, welche braucht
                    ! weniger rechenzeit????

                    ! \frac{\partial\tilde{\phi}_H}{\partial r}
                    dphitdr = -2.5_rk*(1-r)**1.5_rk
                    ! Achtung: neues r noch zu definieren !!!!

                    ! \tilde{\phi}_H
                    phit = (1-r)**2.5_rk

                    x_2s2 = 0.5_rk*chi/sigmasq
                    inv_s = 1.0_rk/sigma

                    r_X2 = r/sqrt2

                    ! ...         ...
                    ! --- (...) + --- (...)
                    ! ...         ...
                    bracket(:) = fracp*(P(i)%o+P(j)%o)+fracm*(P(i)%o-P(j)%o)



                    grad_ri_r(:) = &
                         &(dot_product(bracket,urij)*urij(:)-bracket(:))&
                         &*x_2s2/r&
                         &+inv_s*sqrt1*urij(:)

                    ! force on particle i
                    f(:) = -e*dphitdr*grad_ri_r(:)

                    ! torque on particle i
                    ti(:) = e*chi&
                        &*cross_product(&
                        &r_X2*(-1.0_rk)*dphitdr*(fracm_arij+fracp_arij)*urij(:)&
                        &+(r_X2*chi_half*(-1.0_rk)*dphitdr*(fracm_arij**2-fracp_arij**2)&
                        &+xoioj*phit/(omxoioj*opxoioj))&
                        &*P(j)%o(:)&
                        &,P(i)%o)
                  
                    !P(i)%f(:) = P(i)%f(:) + del(:)*tmp
                    P(i)%f = P(i)%f + f
                    P(i)%t = P(i)%t + ti
                    if (j.le.npmax) then
                       ! torque on particle j
                       tj(:) = e*chi&
                           &*cross_product(&
                           &r_X2*(-1.0_rk)*dphitdr*(fracp_arij-fracm_arij)*urij(:)&
                           &+(r_X2*chi_half*(-1.0_rk)*dphitdr*&
                           &(fracm_arij**2-fracp_arij**2)&
                           &+xoioj*phit/(omxoioj*opxoioj))&
                           &*P(i)%o(:)&
                           &,P(j)%o)
                        !P(j)%f(:) = P(j)%f(:) - del(:)*tmp
                       ! Newton's 3rd law
                       P(j)%f = P(j)%f - f
                       P(j)%t = P(j)%t + tj
                    endif
                end if cutoff
                !endif
                if (lubrication) call particle_lubrication(rsq,del,i,j)
              endif
           enddo
           i = list(i)
        enddo
    end subroutine force_and_torque_hertz_ellipsoid

#ifdef NEWELLIPSOIDMETHOD
    subroutine force_and_torque_hertz_ellipsoid_new_method
        integer i,j,k,ii
        real(kind=rk) ptmp(3),del(3),rsq,r2
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
                cutoff: if (r2<1.0_rk) then
                end if cutoff


           enddo
           i = list(i)
        enddo
    end subroutine force_and_torque_hertz_ellipsoid_new_method
#endif

    subroutine decide_hertz
        if (R_orth.eq.R_para) then
            call force_hertz
        else
               if (polydispersity) call error_md("potential='hertz' does not yet "&
               &//"support ellipsoid polydisperse particles---disable polydispersity "&
               &//"or select potential/='hertz'!")
#ifndef NEWELLIPSOIDMETHOD
            call force_and_torque_hertz_ellipsoid
#else
            call force_and_torque_hertz_ellipsoid_new_method
#endif
        endif
    end subroutine decide_hertz

    subroutine local_energy_hertz(eng)
        real(kind=rk),intent(out) :: eng
        integer i,j,ii
        real(kind=rk) factor,rsq

        eng = 0.0_rk
        i = atompnt
        do ii = 1,nlocal
           do j = nnlist(ii),nnlist(ii+1)-1
              rsq = rdstsq(P(i),P(nlist(j)))
              if (rsq.lt.dsq) then
                 factor = 1.0_rk
                 if (nlist(j).gt.npmax) factor = 0.5_rk
                 eng = eng + phi_hertz(rsq)*factor
              endif
           enddo
           i = list(i)
        enddo
    end subroutine local_energy_hertz

    subroutine local_energy_hertz_softrep(eng)
        real(kind=rk),intent(out) :: eng
        integer i,j,ii
        real(kind=rk) factor,rsq,sr2_hsr,sr6_hsr
	real(kind=rk) r12correctionfactor !energy at the cutoff-distance for shifting U
		  sr2_hsr = sigsq_hsr/rcsq_hsr 
                  sr6_hsr = sr2_hsr*sr2_hsr*sr2_hsr
r12correctionfactor = epsilon4_hsr*sr6_hsr*sr6_hsr
        eng = 0.0_rk
        i = atompnt
        do ii = 1,nlocal
           do j = nnlist(ii),nnlist(ii+1)-1
              rsq = rdstsq(P(i),P(nlist(j)))
              if (rsq.lt.dsq_hsr) then
                 factor = 1.0_rk
                 if (nlist(j).gt.npmax) factor = 0.5_rk 
                 eng = eng + phi_hertz(rsq)*factor
              endif
              if (rsq.lt.rcsq_hsr) then
                 factor = 1.0_rk
                 if (nlist(j).gt.npmax) factor = 0.5_rk
		  sr2_hsr = sigsq_hsr/rsq
                  sr6_hsr = sr2_hsr*sr2_hsr*sr2_hsr
                 eng = eng + (epsilon4_hsr*sr6_hsr*sr6_hsr-r12correctionfactor)*factor
              endif
           enddo
           i = list(i)
        enddo
    end subroutine local_energy_hertz_softrep

    subroutine local_virial_hertz(virial)
        real(kind=rk),intent(out) :: virial
        integer i,j,ii
        real(kind=rk) factor,rsq

        virial = 0.0_rk
        i = atompnt
        do ii = 1,nlocal
           do j = nnlist(ii),nnlist(ii+1)-1
              rsq = rdstsq(P(i),P(nlist(j)))
              if (rsq.lt.dsq) then
                 factor = 1.0_rk
                 if (nlist(j).gt.npmax) factor = 0.5_rk
                 virial = virial + fphi_hertz(rsq)*factor
              endif
           enddo
           i = list(i)
        enddo
    end subroutine local_virial_hertz

    !> derivative of hertz energy (force) ...TIMES THE DISTANCE! (fj,
    !> 2008-02-05)
    real(kind=rk) function fphi_hertz(rsq)
        real(kind=rk),intent(in) :: rsq
        real(kind=rk) r
        r=rsq**0.5_rk
        fphi_hertz = 2.5_rk*K_hertz*(d_hertz-r)**1.5_rk
    end function fphi_hertz

    !> Hertz energy
    real(kind=rk) function phi_hertz(rsq)
        real(kind=rk),intent(in) :: rsq
        real(kind=rk) r

        r=rsq**0.5_rk
        phi_hertz = K_hertz*(d_hertz-r)**2.5_rk
    end function phi_hertz

#endif
end module lbe_md_potential_hertz_module
