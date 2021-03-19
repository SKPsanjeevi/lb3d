#include "lbe.h"

!> Gay-Berne potential (that is Gay-Berne-scaled Lennard-Jones)
!>
!> Implementation of classic Gay-Berne potential as in gay81 but with
!> additional range scaling parameter \c sigma_star which is independent
!> on the range-shifting inherent in the GB formalism. Later
!> publications (for example bates96) usually apply a fixed parameter
!> that is equal to one of the two half-axes of the hard-core volume
!> instead. Our approach is more flexible. As in gay81, \c mu=2 and \c nu=1.
!>
!> The potential was developed having oblate particles
!> \c sigma_orth>sigma_para  in mind. Most likely it will work for prolate
!> particles as well without modifications.
module lbe_md_potential_gb_module
#ifdef MD

    use lbe_helper_module, only: cross_product
    use lbe_globals_module, only: pi,tsize, input_dfile_unit, myrankc
    use lbe_log_module
    use lbe_md_boundary_condition_module, only: rdst
    use lbe_md_fluid_ladd_module, only: particle_lubrication&
         &,rc_lubrication_particle
    use lbe_md_fluid_ladd_parms_module, only: lubrication
    use lbe_md_globals_module
    use lbe_md_helper_module, only: log_msg_md,error_md,log_msg_md_hdr
    use lbe_parallel_module, only: comm_cart
    use lbe_parms_module, only: inp_file,arg_input_dfile_set,arg_input_dfile

    implicit none
    include 'mpif.h'
    private

    public force_and_torque_gb,input_gb,local_energy_gb&
         &,local_virial_gb,setup_gb

    ! Lennard-Jones-equivalent dimensionless cutoff distance - the
    ! actual cutoff is calculated from this number taking into account
    ! shifting and scaling of the range dependency
    real(kind=rk),parameter :: rLJ=2.5_rk

    ! exponents for energy scaling functions - may be changed here,
    ! maybe even to non-integer numbers...
    integer,parameter :: mu=2
    integer,parameter :: nu=1

    ! range-shifting parameters orthogonal and parallel to axis of
    ! rotational symmetry.
    real(kind=rk),save :: sigma_orth=4.0_rk
    real(kind=rk),save :: sigma_para=1.0_rk

    ! universal range-scaling parameter
    real(kind=rk),save :: sigma_star=1.0_rk

    ! energy scaling parameters orthogonal and parallel to axis of
    ! rotational symmetry.
    real(kind=rk),save :: epsilon_orth=1.0_rk
    real(kind=rk),save :: epsilon_para=4.0_rk

    ! range anisotropy parameter
    real(kind=rk),save :: chi

    ! energy anisotropy parameter
    real(kind=rk),save :: chi_prime

    ! scaling constants in GB-scaling functions
    real(kind=rk),save :: sigma0,epsilon0

    ! precalculated constants: 1/sigma_star, 2*mu*chi', sigma0*chi, nu*chi
    real(kind=rk),save :: inv_sigma_star,mu2chi_prime,sigma0chi,nuchi

    ! largest distance between two particles for which--depending on
    ! their orientations--their potential might still be
    ! non-negligibl; and its square
    real(kind=rk),save :: maxr
    real(kind=rk),save :: maxrsq

    namelist /md_potential_gb/ epsilon_orth,epsilon_para,sigma_orth,sigma_para&
         &,sigma_star

contains

    !> read gb section from md input file
    subroutine input_gb
        integer ierror

        call log_msg_md_hdr("Reading MD P GB input")

        ! These features are essential for this coupling module. They are
        ! enabled here because in  setup_gb()  it would be too late.
        use_rotation = .true.
        calculate_orientations = .true.

        if (myrankc.eq.0) then
           open (unit=md_input_file_unit,file=trim(inp_file)//'.md',err=100)
           read (unit=md_input_file_unit,nml=md_potential_gb,err=100)
           close (unit=md_input_file_unit,err=100)
           !call log_msg_md('read /md_potential_gb/ from file "'&
           !     &//trim(inp_file)//'.md"',.false.)
           !write (6,nml=md_potential_gb)
        end if

        if ( arg_input_dfile_set ) then
          call log_msg_md("  Getting differential input...")
          open(UNIT = input_dfile_unit, FILE = arg_input_dfile, STATUS = 'UNKNOWN')
          read(UNIT = input_dfile_unit, NML = md_potential_gb, IOSTAT = ierror)
          if (ierror .ne. 0) then
            call log_msg_md("    WARNING: Differential namelist not found or errors encountered.")
          endif
          close(UNIT = input_dfile_unit)
          call log_ws()
        end if

        write(msgstr,"('epsilon_orth       = ',F16.10)") epsilon_orth
        call log_msg(msgstr)
        write(msgstr,"('epsilon_para       = ',F16.10)") epsilon_para
        call log_msg(msgstr)
        write(msgstr,"('sigma_orth         = ',F16.10)") sigma_orth
        call log_msg(msgstr)
        write(msgstr,"('sigma_para         = ',F16.10)") sigma_para
        call log_msg(msgstr)
        write(msgstr,"('sigma_star         = ',F16.10)") sigma_star
        call log_msg(msgstr)
        call log_ws()

        call MPI_Bcast(epsilon_orth,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(epsilon_para,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(sigma_orth,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(sigma_para,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(sigma_star,1,MPI_REAL8,0,comm_cart,ierror)

        return
100     continue
        call error_md('Error reading md input file "'//trim(inp_file)&
             &//'.md" /md_potential_gb/')
    end subroutine input_gb

    !> calculate constants
    subroutine setup_gb
        ! would be more nice to have this as parameter but then it
        ! would not compile with pathscale anymore...
        real(kind=rk) :: imu
        imu = 1.0_rk/real(mu,kind=rk)

        chi = (sigma_para**2-sigma_orth**2)/(sigma_para**2+sigma_orth**2)
        chi_prime = (epsilon_orth**imu-epsilon_para**imu)&
             &/(epsilon_orth**imu+epsilon_para**imu)

        sigma0 = 2.0_rk*sigma_orth

        epsilon0 = epsilon_orth**(1.0_rk/real(nu,kind=rk))*sqrt(1.0_rk-chi**2)

        inv_sigma_star = 1.0_rk/sigma_star
        mu2chi_prime = 2.0_rk*real(mu,kind=rk)*chi_prime
        sigma0chi = sigma0*chi
        nuchi = real(nu,kind=rk)*chi

        ! take into account both shifting and scaling of distances
        maxr = 1.5_rk*sigma_star+2.0_rk*max(sigma_orth,sigma_para)
        maxrsq = maxr**2

        if (rc<maxr.and.myrankc==0) then
           write (unit=6,fmt='(A,ES15.8)') 'maxr=',maxr
           call error_md('setup_gb(): rc is smaller than maxr.')
        end if

        if (lubrication.and.maxr<rc_lubrication_particle.and.myrankc==0) then
           write (unit=6,fmt='(A,ES15.8,A,ES15.8)') 'maxr=',maxr&
                &,', rc_lubrication_particle=',rc_lubrication_particle
           call error_md('setup_gb(): maxr too small! '&
                &//'maxr<rc_lubrication_particle but lubrication==.true.')
        end if
    end subroutine setup_gb

    subroutine force_and_torque_gb
        integer i,j,k,ii
        real(kind=rk) :: f(3),oi(3),oj(3),rij(3),ti(3),tj(3),urij(3)
        real(kind=rk) :: arij,DD,dPhi,EE,inv_arij,inv_NN,inv_RR,inv_SS,IO,IQ,JO&
             &,JQ,oioj,OO,Phi,PP,QQ,r,RI,rij2,RJ,RR,sigma,SS,X0,X1,XO,XQ,YO,YQ

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
                 oi = P(i)%o
                 oj = P(j)%o

                 oioj = dot_product(oi,oj)

                 XO = chi*oioj

                 arij = sqrt(rij2)
                 inv_arij = 1.0_rk/arij

                 RI = dot_product(rij,oi)*inv_arij
                 RJ = dot_product(rij,oj)*inv_arij

                 RR = 1.0_rk-XO**2
                 IO = RI-XO*RJ
                 JO = RJ-XO*RI

                 YO = -RI*IO-RJ*JO
                 inv_NN = 1.0_rk/(RR+chi*YO)

                 ! orientation dependent range parameter
                 sigma = sigma0*sqrt(RR*inv_NN)

                 cutoff: if (arij<1.5_rk*sigma_star+sigma) then
                    urij = inv_arij*rij

                    inv_RR = 1.0_rk/RR

                    ! normalized input for radial symmetric potential
                    r = inv_sigma_star*(arij-sigma)+1.0_rk

                    ! radial symmetric Lennard-Jones potential and its
                    ! derivative
                    Phi = 4.0_rk*(r**(-12)-r**(-6))
                    dPhi = 24.0_rk*r**(-7)-48.0_rk*r**(-13)

                    XQ = chi_prime*oioj

                    IQ = RI-XQ*RJ
                    JQ = RJ-XQ*RI

                    SS = 1.0_rk-XQ**2
                    inv_SS = 1.0_rk/SS
                    YQ = -RI*IQ-RJ*JQ

                    PP = SS+chi_prime*YQ

                    OO = chi*sigma*inv_NN
                    QQ = mu2chi_prime/PP
                    EE = (epsilon0*sqrt(inv_RR))**nu*(PP*inv_SS)**mu

                    ! scaling of energy - force on particle i
                    f(:) = EE*inv_arij*(Phi*QQ*(IQ*oi+JQ*oj+YQ*urij)&
                         &+inv_sigma_star*dPhi*(OO*(IO*oi+JO*oj+YO*urij)-rij))

                    DD = nuchi*XO*inv_RR
                    X0 = -chi*IO*JO*inv_RR
                    X1 = -chi_prime*IQ*JQ*inv_SS

                    ! torque on particle i
                    ti(:) = EE*cross_product(oi&
                         &,inv_sigma_star*dPhi*OO*(IO*urij+X0*oj)&
                         &+Phi*(QQ*(IQ*urij+X1*oj)-DD*oj))

                    P(i)%f = P(i)%f + f
                    P(i)%t = P(i)%t + ti
                    if (j.le.npmax) then
                       ! Newton's 3rd law
                       P(j)%f = P(j)%f - f

                       ! torque on particle j
                       tj(:) = EE*cross_product(oj&
                            &,inv_sigma_star*dPhi*OO*(JO*urij+X0*oi)&
                            &+Phi*(QQ*(JQ*urij+X1*oi)-DD*oi))

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
              end if lt_maxrsq
           enddo
           i = list(i)
        enddo
    end subroutine force_and_torque_gb

    subroutine local_energy_gb(eng)
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
    end subroutine local_energy_gb

    subroutine local_virial_gb(virial)
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
    end subroutine local_virial_gb

    !> returns the force on particle i exerted by particle j for
    !> orientations  oi(:)  and  oj(:)  and distance vector  rij(:)
    function pair_force(oi,oj,rij)
        real(kind=rk),dimension(3) :: pair_force
        real(kind=rk),intent(in) :: oi(3),oj(3),rij(3)
        real(kind=rk) :: arij,dPhi,EE,inv_arij,inv_NN,IO,IQ,JO,JQ,oioj,OO,Phi&
             &,PP,QQ,r,RI,rij2,RJ,RR,sigma,SS,urij(3),XO,XQ,YO,YQ

        ! (absolute value of rij)**2
        rij2 = dot_product(rij,rij)

        lt_maxrsq: if (rij2<maxrsq) then
           oioj = dot_product(oi,oj)

           XO = chi*oioj

           arij = sqrt(rij2)
           inv_arij = 1.0_rk/arij

           RI = dot_product(rij,oi)*inv_arij
           RJ = dot_product(rij,oj)*inv_arij

           RR = 1.0_rk-XO**2
           IO = RI-XO*RJ
           JO = RJ-XO*RI

           YO = -RI*IO-RJ*JO
           inv_NN = 1.0_rk*(RR+chi*YO)

           ! orientation dependent range parameter
           sigma = sigma0*sqrt(RR*inv_NN)

           cutoff: if (arij<1.5_rk*sigma_star+sigma) then
              urij = inv_arij*rij

              ! normalized input for radial symmetric potential
              r = inv_sigma_star*(arij-sigma)+1.0_rk

              ! radial symmetric Lennard-Jones potential and its derivative
              Phi = 4.0_rk*(r**(-12)-r**(-6))
              dPhi = 24.0_rk*r**(-7)-48.0_rk*r**(-13)

              XQ = chi_prime*oioj

              IQ = RI-XQ*RJ
              JQ = RJ-XQ*RI

              SS = 1.0_rk-XQ**2
              YQ = -RI*IQ-RJ*JQ

              PP = SS+chi_prime*YQ

              OO = chi*sigma*inv_NN
              QQ = mu2chi_prime/PP
              EE = (epsilon0/sqrt(RR))**nu*(PP/SS)**mu

              ! scaling of energy
              pair_force = EE*inv_arij*(Phi*QQ*(IQ*oi+JQ*oj+YQ*urij)&
                   &+inv_sigma_star*dPhi*(OO*(IO*oi+JO*oj+YO*urij)-rij))
           else cutoff
              pair_force = 0.0_rk
           end if cutoff
        end if lt_maxrsq
    end function pair_force

    !> returns the potential energy of two particles with orientations  oi(:)
    !> and  oj(:)  and distance vector  rij(:)
    real(kind=rk) function pair_potential(oi,oj,rij)
        real(kind=rk),intent(in) :: oi(3),oj(3),rij(3)
        real(kind=rk) :: arij,EE,inv_arij,inv_NN,IO,IQ,JO,JQ,oioj,Phi,PP,r,RI&
             &,rij2,RJ,RR,sigma,SS,XO,XQ,YO,YQ

        ! (absolute value of rij)**2
        rij2 = dot_product(rij,rij)

        lt_maxrsq: if (rij2<maxrsq) then
           oioj = dot_product(oi,oj)

           XO = chi*oioj

           arij = sqrt(rij2)
           inv_arij = 1.0_rk/arij

           RI = dot_product(rij,oi)*inv_arij
           RJ = dot_product(rij,oj)*inv_arij

           RR = 1.0_rk-XO**2
           IO = RI-XO*RJ
           JO = RJ-XO*RI

           YO = -RI*IO-RJ*JO
           inv_NN = 1.0_rk/(RR+chi*YO)

           ! orientation dependent range parameter
           sigma = sigma0*sqrt(RR*inv_NN)

           cutoff: if (arij<1.5_rk*sigma_star+sigma) then
              ! normalized input for radial symmetric potential
              r = inv_sigma_star*(arij-sigma)+1.0_rk

              ! radial symmetric Lennard-Jones potential
              Phi = 4.0_rk*(r**(-12)-r**(-6))

              XQ = chi_prime*oioj

              IQ = RI-XQ*RJ
              JQ = RJ-XQ*RI

              SS = 1.0_rk-XQ**2
              YQ = -RI*IQ-RJ*JQ

              PP = SS+chi_prime*YQ

              EE = (epsilon0/sqrt(RR))**nu*(PP/SS)**mu

              ! scaling of energy
              pair_potential = EE*Phi
           else cutoff
              pair_potential = 0.0_rk
           end if cutoff
        end if lt_maxrsq
    end function pair_potential

#endif
end module lbe_md_potential_gb_module
