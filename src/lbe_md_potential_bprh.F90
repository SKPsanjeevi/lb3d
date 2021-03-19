#include "lbe.h"

!> Berne-Pechukas-scaled repulsive Hooke potential
!>
!> Implementation of a repulsive Hooke potential \f$(1-r)^2, 0<r<1\f$ fitted
!> to pairs of two-axial ellipsoids as proposed by Berne and Pechukas
!> in J. Chem. Phys. 56 (8), 4213 (1972)
module lbe_md_potential_bprh_module
#ifdef MD

    use lbe_helper_module, only: cross_product,every_n_time_steps
    use lbe_globals_module, only: pi,tsize, input_dfile_unit,myrankc
    use lbe_log_module
    use lbe_md_boundary_condition_module, only: rdst
    use lbe_md_fluid_ladd_module, only: particle_lubrication&
         &,rc_lubrication_particle
    use lbe_md_fluid_ladd_parms_module, only: lubrication
    use lbe_md_globals_module
    use lbe_md_helper_module, only: number_density,log_msg_md,error_md,log_msg_md_hdr
    use lbe_md_output_module, only: dump_potentials,n_dump
    use lbe_parallel_module, only: comm_cart
    use lbe_parms_module, only: inp_file,arg_input_dfile_set,arg_input_dfile

    implicit none
    include 'mpif.h'
    private

    public force_and_torque_bprh,input_bprh&
         &,local_energy_bprh,local_virial_bprh,set_growth_factor_bprh&
         &,setup_bprh

    !> \name input file parameters
    !> \{
    !> strength/energy parameter
    real(kind=rk),save :: epsilon=1.0_rk
    !> range parameter orthogonal to axis of rotational symmetry
    real(kind=rk),save :: sigma_orth=4.0_rk
    !> range parameter parallel to axis of rotational symmetry
    real(kind=rk),save :: sigma_para=1.0_rk
    !> \}

    ! anisotropy parameter \chi and its half
    real(kind=rk),save :: chi
    real(kind=rk),save :: chi_half
    ! range parameter \sigma and its square
    real(kind=rk),save :: sigma
    real(kind=rk),save :: sigmasq

    ! largest distance between two particles for which--depending on
    ! their orientations--their potential might still be greater than
    ! zero; and its square
    real(kind=rk),save :: maxr
    real(kind=rk),save :: maxrsq

    ! \frac{\chi}{2\sigma^2}
    real(kind=rk),save :: x_2s2

    ! \frac{1}{\sigma}
    real(kind=rk),save :: inv_s

    namelist /md_potential_bprh/ epsilon,sigma_orth,sigma_para

contains

    !> read bprh section from md input file
    subroutine input_bprh
        integer ierror

        call log_msg_md_hdr("Reading MD P BPRH input")

        ! These features are essential for this coupling module. They are
        ! enabled here because in  setup_bprh()  it would be too late.
        use_rotation = .true.
        calculate_orientations = .true.

        if (myrankc.eq.0) then
           open (unit=md_input_file_unit,file=trim(inp_file)//'.md',err=100)
           read (unit=md_input_file_unit,nml=md_potential_bprh,err=100)
           close (unit=md_input_file_unit,err=100)
           !call log_msg_md('read /md_potential_bprh/ from file "'&
           !     &//trim(inp_file)//'.md"',.false.)
           !write (6,nml=md_potential_bprh)
        endif

        if ( arg_input_dfile_set ) then
          call log_msg_md("  Getting differential input...")
          open(UNIT = input_dfile_unit, FILE = arg_input_dfile, STATUS = 'UNKNOWN')
          read(UNIT = input_dfile_unit, NML = md_potential_bprh, IOSTAT = ierror)
          if (ierror .ne. 0) then
            call log_msg_md("    WARNING: Differential namelist not found or errors encountered.")
          endif
          close(UNIT = input_dfile_unit)
          call log_ws()
        end if

        write(msgstr,"('epsilon            = ',F16.10)") epsilon
        call log_msg(msgstr)
        write(msgstr,"('sigma_orth         = ',F16.10)") sigma_orth
        call log_msg(msgstr)
        write(msgstr,"('sigma_para         = ',F16.10)") sigma_para
        call log_msg(msgstr)
        call log_ws()

        call MPI_Bcast(epsilon,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(sigma_orth,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(sigma_para,1,MPI_REAL8,0,comm_cart,ierror)

        return
100     continue
        call error_md('Error reading md input file "'//trim(inp_file)//'.md"')
    end subroutine input_bprh

    !> precompute actual length parameters from \c sigma_para and \c
    !> sigma_orth
    !>
    !> \param[in] f (optional, defaults to \c 1) linear length scale
    !> factor with respect to original lengths read from input file,
    !> needed for \c growth_stage
    !>
    !> Note that except for initialization, \c sigma_para and \c
    !> sigma_orth are not used but instead the parameters computed
    !> here.
    subroutine precompute_length_parameters(f)
        real(kind=rk),intent(in),optional :: f
        real(kind=rk) :: factor

        if (present(f)) then
           factor = f
        else
           factor = 1.0_rk
        end if

        ! The factor must differ by sqrt(2) from the one in the paper
        ! to achieve zero potential when two ellipsoids are just
        ! touching each other.
        sigma = 2.0_rk*sigma_orth*factor
        sigmasq = sigma**2

        x_2s2 = 0.5_rk*chi/sigmasq

        inv_s = 1.0_rk/sigma

        maxr = 2.0_rk*max(sigma_orth,sigma_para)*factor
        maxrsq = maxr**2

        if (rc<maxr.and.myrankc==0) call error_md('potential_bprh: rc is '&
             &//'smaller than 2*max(sigma_orth,sigma_para). It must be at '&
             &//'least equally large. If final_length_factor>1 in '&
             &//'&md_growing_stage, rc must be enlarged by this factor.')
    end subroutine precompute_length_parameters

    !> rescales all length parameters of the BPrH potential
    !>
    !> \param[in] f linear length scale factor with respect to
    !> original lengths read from input file
    subroutine set_growth_factor_bprh(f)
        real(kind=rk),intent(in) :: f

        call precompute_length_parameters(f)
    end subroutine set_growth_factor_bprh

    !> calculate constants and check input parameters
    subroutine setup_bprh
        ! these parameters depend on the aspect ratio only, not on the
        ! absolute lengths
        chi = (sigma_para**2-sigma_orth**2)/(sigma_para**2+sigma_orth**2)
        chi_half = 0.5_rk*chi

        ! compute actual length parameters from sigma_para and sigma_orth
        call precompute_length_parameters

        if (lubrication.and.maxr<rc_lubrication_particle.and.myrankc==0) &
             &call error_md('2*max(sigma_o,sigma_p) too small! '&
             &//'(2*max(sigma_o,sigma_p)<rc_lubrication_particle '&
             &//'but lubrication==.true.')
    end subroutine setup_bprh

    subroutine force_and_torque_bprh
        integer i,j,k,ii
        real(kind=rk) :: arij,bracket(3),e,e_pot_2,eomr,f(3),fracm,fracm_arij&
             &,fracp,fracp_arij,grad_ri_r(3),omr,omxoioj,opxoioj,r,r_X2,r2&
             &,rij(3),rij2,rijoi,rijoimrijoj,rijoiprijoj,rijoj,s2,sqrt1,sqrt2&
             &,ti(3),tj(3),urij(3),xoioj
        logical :: energy_dump_follows

        energy_dump_follows = every_n_time_steps(n_dump).and.dump_potentials

        i = atompnt
        do ii = 1,nlocal
           do k = nnlist(ii),nnlist(ii+1)-1
              j = nlist(k)
              rij = rdst(P(i),P(j))

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

                 cutoff: if (r2<1.0_rk) then
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

                    ! epsilon*(1-r)
                    eomr = e*omr

                    ! force on particle i
                    f(:) = 2.0_rk*eomr*grad_ri_r(:)

                    r_X2 = r/sqrt2

                    ! torque on particle i
                    ti(:) = eomr*chi&
                         &*cross_product(&
                         &r_X2*(fracm_arij+fracp_arij)*urij(:)&
                         &+(r_X2*chi_half*(fracm_arij**2-fracp_arij**2)&
                         &+xoioj*omr/(omxoioj*opxoioj))&
                         &*P(j)%o(:)&
                         &,P(i)%o)

                    P(i)%f = P(i)%f + f
                    P(i)%t = P(i)%t + ti
                    if (j.le.npmax) then
                       ! Newton's 3rd law
                       P(j)%f = P(j)%f - f

                       ! torque on particle j
                       tj(:) = eomr*chi&
                            &*cross_product(&
                            &r_X2*(fracp_arij-fracm_arij)*urij(:)&
                            &+(r_X2*chi_half*(fracm_arij**2-fracp_arij**2)&
                            &+xoioj*omr/(omxoioj*opxoioj))&
                            &*P(i)%o(:)&
                            &,P(j)%o)

                       P(j)%t = P(j)%t + tj
                    end if

                    ! The conditions do not work if a periodic
                    ! boundary lies between both particles.  The force
                    ! is always counted on that process that owns the
                    ! particle the force is acting on.
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

                 if (energy_dump_follows) then
                    ! share potential energy between involved
                    ! particles, every particles gets potential energy
                    ! from its owning process only
                    e_pot_2 = 0.5_rk*pair_potential(P(i)%o,P(j)%o,rij)
                    P(i)%e_pot = P(i)%e_pot + e_pot_2
                    if (j<=npmax) P(j)%e_pot = P(j)%e_pot + e_pot_2
                 end if
              end if lt_maxrsq
           enddo
           i = list(i)
        enddo
    end subroutine force_and_torque_bprh

    subroutine local_energy_bprh(eng)
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
    end subroutine local_energy_bprh

    subroutine local_virial_bprh(virial)
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
    end subroutine local_virial_bprh

    !> returns the force on particle i exerted by particle j for
    !> orientations  oi(:)  and  oj(:)  and distance vector  rij(:)
    function pair_force(oi,oj,rij)
        real(kind=rk),dimension(3) :: pair_force
        real(kind=rk),intent(in) :: oi(3),oj(3),rij(3)
        real(kind=rk) :: arij,bracket(3),e,fracm,fracp,grad_ri_r(3),omxoioj&
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

        cutoff: if (r2<1.0_rk) then
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

           pair_force(:) = 2.0_rk*e*(1.0_rk-r)*grad_ri_r(:)
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

        cutoff: if (r2<1.0_rk) then
           ! angle dependent strength parameter \epsilon
           e = epsilon/sqrt(omxoioj*opxoioj)

           ! normalized input for radial symmetric potential \Phi(r)
           r = sqrt(r2)

           ! radial symmetric repulsive Hooke potential \Phi(r)
           Phi = (1.0_rk-r)**2

           pair_potential = e*Phi
        else cutoff
           pair_potential = 0.0_rk
        end if cutoff
    end function pair_potential

#endif
end module lbe_md_potential_bprh_module
