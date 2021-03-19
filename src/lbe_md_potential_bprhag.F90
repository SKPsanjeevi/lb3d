#include "lbe.h"

!> Berne-Pechukas-scaled repulsive Hooke potential with an additional
!> attractive Gaussian
!>
!> Implementation of a repulsive Hooke potential \f$(1-r)^2-\delta, 0<r<1\f$
!> and an attractive Gaussian \f$-\delta*\exp[-4(1-r)^2], 1<r\f$ fitted to
!> pairs of two-axial ellipsoids as proposed by Berne and Pechukas in
!> J. Chem. Phys. 56 (8), 4213 (1972)
module lbe_md_potential_bprhag_module
#ifdef MD
    use lbe_helper_module, only: cross_product
    use lbe_globals_module, only: pi,tsize, input_dfile_unit,myrankc
    use lbe_log_module
    use lbe_md_boundary_condition_module, only: rdst
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

    public force_and_torque_bprhag,input_bprhag,local_energy_bprhag&
         &,local_virial_bprhag,setup_bprhag

    !> inverse gaussian width parameter
    !>
    !> \note This might change in the future.
    real(kind=rk),parameter :: inv_2sigma2=4.0_rk

    !> hard-coded cutoff distance in reduced units
    !>
    !> \note This might change in the future.
    real(kind=rk),parameter :: r_cut=3.0_rk

    !> square of \c r_cut
    !>
    !> \note This might change in the future.
    real(kind=rk),parameter :: r2_cut=r_cut**2

    !> \{
    !> \name Potential parameters
    !>
    !> strength/energy and range parameters orthogonal and parallel to
    !> axis of rotational symmetry. Be aware that also the Gaussian
    !> part of the potential is scaled by epsilon, so delta determines
    !> not the absolute potential dip but its deepness as a fraction
    !> of epsilon.
    real(kind=rk),save :: epsilon=1.0_rk
    real(kind=rk),save :: sigma_orth=4.0_rk
    real(kind=rk),save :: sigma_para=1.0_rk
    real(kind=rk),save :: delta=0.1_rk
    !> \}

    !> anisotropy parameter \f$ \chi \f$
    real(kind=rk),save :: chi
    !> \f$ \chi/2 \f$
    real(kind=rk),save :: chi_half
    !> range parameter \f$ \sigma \f$
    real(kind=rk),save :: sigma
    !>  \f$ \sigma^2 \f$
    real(kind=rk),save :: sigmasq

    !> largest distance between two particles for which--depending on
    !> their orientations--their potential might still be greater than
    !> zero
    real(kind=rk),save :: maxr
    real(kind=rk),save :: maxrsq !< square of \c maxr

    !> \f$ \frac{\chi}{2\sigma^2} \f$
    real(kind=rk),save :: x_2s2

    !> \f$ \frac{1}{\sigma} \f$
    real(kind=rk),save :: inv_s

    namelist /md_potential_bprhag/ epsilon,sigma_orth,sigma_para,delta

contains

    !> read bprhag section from md input file
    subroutine input_bprhag
        integer ierror

        call log_msg_md_hdr("Reading MD P BPRHAG input")
        ! These features are essential for this coupling module. They are
        ! enabled here because in  setup_bprhag()  it would be too late.
        use_rotation = .true.
        calculate_orientations = .true.

        if (myrankc.eq.0) then
           open (unit=md_input_file_unit,file=trim(inp_file)//'.md',err=100)
           read (unit=md_input_file_unit,nml=md_potential_bprhag,err=100)
           close (unit=md_input_file_unit,err=100)
           !call log_msg_md('read /md_potential_bprhag/ from file "'&
           !     &//trim(inp_file)//'.md"',.false.)
           !write (6,nml=md_potential_bprhag)
        end if

        if ( arg_input_dfile_set ) then
          call log_msg_md("  Getting differential input...")
          open(UNIT = input_dfile_unit, FILE = arg_input_dfile, STATUS = 'UNKNOWN')
          read(UNIT = input_dfile_unit, NML = md_potential_bprhag, IOSTAT = ierror)
          if (ierror .ne. 0) then
            call log_msg_md("    WARNING: Differential namelist not found or errors encountered.")
          end if
          close(UNIT = input_dfile_unit)
          call log_ws()
        end if

        write(msgstr,"('epsilon            = ',F16.10)") epsilon
        call log_msg(msgstr)
        write(msgstr,"('sigma_orth         = ',F16.10)") sigma_orth
        call log_msg(msgstr)
        write(msgstr,"('sigma_para         = ',F16.10)") sigma_para
        call log_msg(msgstr)
        write(msgstr,"('delta              = ',F16.10)") delta
        call log_msg(msgstr)
        call log_ws()

        call MPI_Bcast(epsilon,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(sigma_orth,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(sigma_para,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(delta,1,MPI_REAL8,0,comm_cart,ierror)

        return
100     continue
        call error_md('Error reading md input file "'//trim(inp_file)&
             &//'.md" /md_potential_bprhag/')
    end subroutine input_bprhag

    !> calculate constants
    subroutine setup_bprhag
        maxr = 2.0_rk*r_cut*max(sigma_orth,sigma_para)
        maxrsq = maxr**2

        if (rc<maxr.and.myrankc==0) call error_md(&
             &'rc is smaller than 2*r_cut*max(sigma_orth,sigma_para).')

        if (lubrication.and.maxr<rc_lubrication_particle.and.myrankc==0) &
             &call error_md('2*r_cut*max(sigma_orth,sigma_para) too small! '&
             &//'(2*r_cut*max(sigma_orth,sigma_para)<rc_lubrication_particle '&
             &//'but lubrication==.true.')

        chi = (sigma_para**2-sigma_orth**2)/(sigma_para**2+sigma_orth**2)
        chi_half = 0.5_rk*chi

        ! The factor must differ by sqrt(2) from the one in the paper
        ! to achieve zero potential when two ellipsoids are just
        ! touching each other.
        sigma = 2.0_rk*sigma_orth
        sigmasq = sigma**2

        x_2s2 = 0.5_rk*chi/sigmasq

        inv_s = 1.0_rk/sigma
    end subroutine setup_bprhag

    subroutine force_and_torque_bprhag
        integer i,j,k,ii
        real(kind=rk) :: arij,bracket(3),dphi_dr,e,f(3),fracm,fracm_arij,fracp&
             &,fracp_arij,grad_ri_r(3),omr,omxoioj,opxoioj,phi,r,r2,rdphi_X22dr&
             &,rij(3),rij2,rijoi,rijoimrijoj,rijoiprijoj,rijoj,s2,sqrt1,sqrt2&
             &,ti(3),tj(3),urij(3),xoioj

        i = atompnt
        do ii = 1,nlocal
           do k = nnlist(ii),nnlist(ii+1)-1
              j = nlist(k)
              rij = P(i)%x - P(j)%x
              ! ATTENTION: optimized minimum image criterion - result is
              ! useless if its square is greater than cutsq1
              if (abs(rij(1)).gt.tsize_mrc(1)) then
                 if (rij(1).lt.0.0_rk) then
                    rij(1) = rij(1) + tsize(1)
                 else
                    rij(1) = rij(1) - tsize(1)
                 endif
              endif
              if (abs(rij(2)).gt.tsize_mrc(2)) then
                 if (rij(2).lt.0.0_rk) then
                    rij(2) = rij(2) + tsize(2)
                 else
                    rij(2) = rij(2) - tsize(2)
                 endif
              endif
              if (abs(rij(3)).gt.tsize_mrc(3)) then
                 if (rij(3).lt.0.0_rk) then
                    rij(3) = rij(3) + tsize(3)
                 else
                    rij(3) = rij(3) - tsize(3)
                 endif
              endif

              ! (absolute value of rij)**2
              rij2 = dot_product(rij,rij)

              lt_maxrsq: if (rij2<maxrsq) then
                 ! dot products of rij and o(i|j)
                 rijoi = dot_product(rij,P(i)%o)
                 rijoj = dot_product(rij,P(j)%o)

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

                 ! \sqrt{1-\frac{\chi}{2}\left[\ldots\right]}^2
                 sqrt2 = 1.0_rk&
                      &-(rijoiprijoj*fracp+rijoimrijoj*fracm)*chi_half/rij2

                 ! (angle dependent range parameter \sigma)**2
                 s2 = sigmasq/sqrt2

                 ! (normalized input for radial symmetric potential \Phi(r))**2
                 r2 = rij2/s2

                 cutoff: if (r2<r2_cut) then
                    ! angle dependent strength parameter \epsilon
                    e = epsilon/sqrt(omxoioj*opxoioj)

                    ! \sqrt{1-\frac{\chi}{2}\left[\ldots\right]}
                    sqrt1 = sqrt(sqrt2)

                    ! |\mathbf{r}_{ij}|
                    arij = sqrt(rij2)

                    ! unit vector in direction of rij
                    urij(:) = rij(:)/arij

                    ! frac(p|m) with rij replaced by urij
                    fracm_arij = fracm/arij
                    fracp_arij = fracp/arij

                    ! normalized input for radial symmetric potential \Phi(r)
                    r = arij*sqrt1/sigma

                    ! ...         ...
                    ! --- (...) + --- (...)
                    ! ...         ...
                    bracket(:) = fracp*(P(i)%o+P(j)%o)+fracm*(P(i)%o-P(j)%o)

                    ! \frac{\partial r}{\partial\mathbf{r}_i}
                    grad_ri_r(:) = &
                         &(dot_product(bracket,urij)*urij(:)-bracket(:))&
                         &*x_2s2/r&
                         &+inv_s*sqrt1*urij(:)

                    ! 1-r
                    omr = 1.0_rk-r

                    ! \Phi(r) and \frac{\partial\Phi(r)}{\partial r}
                    rh_ag: if (r<1.0_rk) then
                       phi = omr**2-delta
                       dphi_dr = -2.0_rk*omr
                    else rh_ag
                       phi = -delta*exp(-inv_2sigma2*omr**2)
                       dphi_dr = -2.0_rk*inv_2sigma2*delta&
                            &*exp(-inv_2sigma2*omr**2)*omr
                    end if rh_ag

                    ! force on particle i
                    f(:) = -e*dphi_dr*grad_ri_r(:)

                    rdphi_X22dr = 0.5_rk*dphi_dr*r/sqrt2

                    ! torque on particle i
                    ti(:) = e*chi&
                         &*cross_product(&
                         &P(i)%o&
                         &,rdphi_X22dr*(fracm_arij+fracp_arij)*urij(:)&
                         &+(rdphi_X22dr*chi_half&
                         &*(fracm_arij**2-fracp_arij**2)&
                         &-phi*xoioj/(omxoioj*opxoioj))&
                         &*P(j)%o(:))

                    P(i)%f = P(i)%f + f
                    P(i)%t = P(i)%t + ti
                    if (j.le.npmax) then
                       ! Newton's 3rd law
                       P(j)%f = P(j)%f - f

                       ! torque on particle j
                       tj(:) = e*chi&
                            &*cross_product(&
                            &P(j)%o&
                            &,rdphi_X22dr*(fracp_arij-fracm_arij)*urij(:)&
                            &+(rdphi_X22dr*chi_half&
                            &*(fracm_arij**2-fracp_arij**2)&
                            &-phi*xoioj/(omxoioj*opxoioj))&
                            &*P(i)%o(:))

                       P(j)%t = P(j)%t + tj
                    end if

                    ! The conditions do not work if a periodic
                    ! boundary lies between both particles.  The
                    ! force is always counted on that process that
                    ! owns the particle the force is acting on.
                    if (P(i)%x(1)<real(dfz_minx,kind=rk)-0.5_rk&
                         &.and.P(j)%x(1)>=real(dfz_minx,kind=rk)-0.5_rk&
                         &.and.j.le.npmax)&
                         &dfz_pp(1) = dfz_pp(1) - f(3)
                    if (P(j)%x(1)<real(dfz_minx,kind=rk)-0.5_rk&
                         &.and.P(i)%x(1)>=real(dfz_minx,kind=rk)-0.5_rk)&
                         &dfz_pp(1) = dfz_pp(1) + f(3)
                    if (P(i)%x(1)<real(dfz_maxx,kind=rk)+0.5_rk&
                         &.and.P(j)%x(1)>=real(dfz_maxx,kind=rk)+0.5_rk)&
                         &dfz_pp(2) = dfz_pp(2) + f(3)
                    if (P(j)%x(1)<real(dfz_maxx,kind=rk)+0.5_rk&
                         &.and.P(i)%x(1)>=real(dfz_maxx,kind=rk)+0.5_rk&
                         &.and.j.le.npmax)&
                         &dfz_pp(2) = dfz_pp(2) - f(3)
                 end if cutoff
                 if (lubrication) call particle_lubrication(rij2,rij,i,j)
              end if lt_maxrsq
           enddo
           i = list(i)
        enddo
    end subroutine force_and_torque_bprhag

    subroutine local_energy_bprhag(eng)
        real(kind=rk),intent(out) :: eng
        integer i,j,ii
        real(kind=rk) :: factor,r(3)

        eng = 0.0_rk
        i = atompnt
        do ii = 1,nlocal
           do j = nnlist(ii),nnlist(ii+1)-1
              r = rdst(P(i),P(nlist(j)))
              if (dot_product(r,r).lt.maxrsq) then
                 factor = 1.0_rk
                 if (nlist(j).gt.npmax) factor = 0.5_rk
                 eng = eng + factor*pair_potential(P(i)%o,P(nlist(j))%o,r)
              endif
           enddo
           i = list(i)
        enddo
    end subroutine local_energy_bprhag

    subroutine local_virial_bprhag(virial)
        real(kind=rk),intent(out) :: virial
        integer i,j,ii
        real(kind=rk) factor,r(3)

        virial = 0.0_rk
        i = atompnt
        do ii = 1,nlocal
           do j = nnlist(ii),nnlist(ii+1)-1
              r = rdst(P(i),P(nlist(j)))
              if (dot_product(r,r).lt.maxrsq) then
                 factor = 1.0_rk
                 if (nlist(j).gt.npmax) factor = 0.5_rk
                 virial = virial&
                      &+factor*dot_product(pair_force(P(i)%o,P(nlist(j))%o,r),r)
              endif
           enddo
           i = list(i)
        enddo
    end subroutine local_virial_bprhag

    !> returns the force on particle i exerted by particle j for
    !> orientations  oi(:)  and  oj(:)  and distance vector  rij(:)
    function pair_force(oi,oj,rij)
        real(kind=rk),dimension(3) :: pair_force
        real(kind=rk),intent(in) :: oi(3),oj(3),rij(3)
        real(kind=rk) :: arij,bracket(3),e,fracm,fracp,grad_ri_r(3),omr,omxoioj&
             &,opxoioj,r,r2,rij2,rijoi,rijoimrijoj,rijoiprijoj,rijoj,s2,sqrt1&
             &,sqrt2,urij(3),xoioj

        ! dot products of rij and o(i|j)
        rijoi = dot_product(rij,oi)
        rijoj = dot_product(rij,oj)

        ! \chi(\mathbf{o}_i\mathbf{o}_j)
        xoioj = chi*dot_product(oi,oj)

        omxoioj = 1.0_rk-xoioj
        opxoioj = 1.0_rk+xoioj

        ! (absolute value of rij)**2
        rij2 = dot_product(rij,rij)

        ! \mathbf{r}_{ij}\mathbf{o}_i\pm\mathbf{r}_{ij}\mathbf{o}_j)
        ! ----------------------------------------------------------
        !              1\pm\chi\mathbf{o}_i\mathbf{o}_j
        rijoiprijoj = rijoi+rijoj
        rijoimrijoj = rijoi-rijoj
        fracp = rijoiprijoj/opxoioj
        fracm = rijoimrijoj/omxoioj

        ! \sqrt{1-\frac{\chi}{2}\left[\ldots\right]}^2
        sqrt2 = 1.0_rk-(rijoiprijoj*fracp+rijoimrijoj*fracm)*chi_half/rij2

        ! (angle dependent range parameter \sigma)**2
        s2 = sigmasq/sqrt2

        ! (normalized input for radial symmetric potential \Phi(r))**2
        r2 = rij2/s2

        cutoff: if (r2<r2_cut) then
           ! angle dependent strength parameter \epsilon
           e = epsilon/sqrt(omxoioj*opxoioj)

           ! \sqrt{1-\frac{\chi}{2}\left[\ldots\right]}
           sqrt1 = sqrt(sqrt2)

           ! |\mathbf{r}_{ij}|
           arij = sqrt(rij2)

           ! unit vector in direction of rij
           urij(:) = rij(:)/arij

           ! normalized input for radial symmetric potential \Phi(r)
           r = arij*sqrt1/sigma

           ! ...         ...
           ! --- (...) + --- (...)
           ! ...         ...
           bracket(:) = fracp*(oi(:)+oj(:))+fracm*(oi(:)-oj(:))

           ! \frac{\partial r}{\partial\mathbf{r}_i}
           grad_ri_r(:) = &
                &(dot_product(bracket,urij)*urij(:)-bracket(:))*x_2s2/r&
                &+inv_s*sqrt1*urij(:)

           omr = 1.0_rk-r

           ! repulsive Hooke or attractive Gaussian?
           rh_ag: if (r<=1.0_rk) then
              pair_force(:) = e*2.0_rk*omr*grad_ri_r(:)
           else rh_ag
              pair_force(:) = e*2.0_rk*inv_2sigma2*delta&
                   &*exp(-inv_2sigma2*omr**2)*omr*grad_ri_r(:)
           end if rh_ag
        else cutoff
           pair_force(:) = 0.0_rk
        end if cutoff
    end function pair_force

    !> returns the potential energy of two particles with orientations  oi(:)
    !> and  oj(:)  and distance vector  rij(:)
    real(kind=rk) function pair_potential(oi,oj,rij)
        real(kind=rk),intent(in) :: oi(3),oj(3),rij(3)
        real(kind=rk) :: e,omxoioj,opxoioj,Phi,r,r2,rij2,rijoi,rijoj,s2,xoioj

        rijoi = dot_product(rij,oi)
        rijoj = dot_product(rij,oj)

        ! \chi(\mathbf{o}_i\mathbf{o}_j)
        xoioj = chi*dot_product(oi,oj)

        omxoioj = 1.0_rk-xoioj
        opxoioj = 1.0_rk+xoioj

        ! (absolute value of rij)**2
        rij2 = dot_product(rij,rij)

        ! (angle dependent range parameter \sigma)**2
        s2 = sigmasq&
             &/(1.0_rk&
             &-((rijoi+rijoj)**2/(opxoioj)+(rijoi-rijoj)**2/(omxoioj))&
             &*chi_half/rij2)

        ! (normalized input for radial symmetric potential \Phi(r))**2
        r2 = rij2/s2

        cutoff: if (r2<r2_cut) then
           ! angle dependent strength parameter \epsilon
           e = epsilon/sqrt(omxoioj*opxoioj)

           ! normalized input for radial symmetric potential \Phi(r)
           r = sqrt(r2)

           ! radial symmetric repulsive Hooke/attractive Gaussian
           ! potential \Phi(r)
           rh_ag: if (r<=1.0_rk) then
              Phi = (1.0_rk-r)**2-delta
           else rh_ag
              Phi = -delta*exp(-inv_2sigma2*(1.0_rk-r)**2)
           end if rh_ag

           pair_potential = e*Phi
        else cutoff
           pair_potential = 0.0_rk
        end if cutoff
    end function pair_potential

#endif
end module lbe_md_potential_bprhag_module
