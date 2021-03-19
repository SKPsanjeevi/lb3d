#include "lbe.h"

!> Berne-Pechukas-scaled repulsive Hooke potential for particel-wall
!> interaction
!>
!> Repulsive Hooke potential \f$(1-r)^2\f$, \f$0<r<1\f$ fitted to
!> pairs of a two-axial ellipsoid (md particle) and a sphere (surface
!> rock node) as proposed by Berne and Pechukas in J. Chem. Phys. 56
!> (8), 4213 (1972)
module lbe_md_rock_bprh_module
#ifdef MD
    use lbe_bc_module, only: periodically_wrapped
    use lbe_globals_module, only: halo_extent, input_dfile_unit,myrankc
    use lbe_helper_module, only: cross_product,local_coordinates
    use lbe_log_module
    use lbe_md_boundary_condition_module, only: local_chunks,local_chunk_type&
         &,n_max_local_chunks
    use lbe_md_fluid_ladd_module, only: lubrication_force_rock&
         &,rc_lubrication_rock
    use lbe_md_fluid_ladd_parms_module, only: lubrication
    use lbe_md_globals_module
    use lbe_md_helper_module, only: orientation,log_msg_md,error_md,log_msg_md_hdr
    use lbe_parallel_module, only: comm_cart
    use lbe_parms_module, only: inp_file,arg_input_dfile_set,arg_input_dfile
    use lbe_types_module, only: lbe_site

    implicit none
    include 'mpif.h'
    private

    public input_rock_bprh,local_rock_potential_bprh&
         &,particle_rock_potential_bprh,rock_ft_interaction_bprh,setup_rock_bprh

    ! strength/energy and range parameters orthogonal and parallel to
    ! axis of rotational symmetry of the ellipsoid and of the sphere assumed
    ! at each surface rock site
    real(kind=rk),save :: epsilon=1.0_rk
    real(kind=rk),save :: sigma_orth=4.0_rk
    real(kind=rk),save :: sigma_para=1.0_rk
    real(kind=rk),save :: sigma=0.5_rk

    ! largest distance between a particle and a surface rock site for
    ! which--depending on the particle orientation--their potential might
    ! still be greater than zero; and its square
    real(kind=rk),save :: r_cut
    real(kind=rk),save :: r_cut_sq

    ! anisotropy parameter \chi
    real(kind=rk),save :: chi
    ! square and inverse of the range parameter \overline{\sigma}
    real(kind=rk),save :: sigma_barsq
    real(kind=rk),save :: inv_sigma_bar

    ! 2*epsilon/sigma_bar
    real(kind=rk),save :: epsilon_2_by_sigma_bar

    namelist /md_rock_bprh/ epsilon,sigma,sigma_orth,sigma_para

contains

    !> read namelist \c /md_rock_bprh/  from input file
    subroutine input_rock_bprh
        integer ierror

        call log_msg_md_hdr("Reading MD R BPRH input")

        ! These features are essential for this rock module. They are
        ! enabled here because in  setup_rock_bprh()  it would be too late.
        use_rotation = .true.
        calculate_orientations = .true.
        collect_forces = .true.

        if (myrankc.eq.0) then
           open (unit=md_input_file_unit,file=trim(inp_file)//'.md',err=100)
           read (unit=md_input_file_unit,nml=md_rock_bprh,err=100)
           close (unit=md_input_file_unit,err=100)
           !call log_msg_md('read /md_rock_bprh/ from file "'//trim(inp_file)&
           !     &//'.md"',.false.)
           !write (6,nml=md_rock_bprh)
        end if

        if ( arg_input_dfile_set ) then
          call log_msg_md("  Getting differential input...")
          open(UNIT = input_dfile_unit, FILE = arg_input_dfile, STATUS = 'UNKNOWN')
          read(UNIT = input_dfile_unit, NML = md_rock_bprh, IOSTAT = ierror)
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
        write(msgstr,"('sigma              = ',F16.10)") sigma
        call log_msg(msgstr)
        call log_ws()

        call MPI_Bcast(epsilon,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(sigma_orth,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(sigma_para,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(sigma,1,MPI_REAL8,0,comm_cart,ierror)

        r_cut = sigma+max(sigma_orth,sigma_para)

        return
100     continue
        call error_md('Error reading md input file "'//trim(inp_file)//'.md"')
    end subroutine input_rock_bprh

    subroutine local_rock_potential_bprh(N,pot)
        type(lbe_site),intent(in) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        real(kind=rk),intent(out) :: pot
        integer i,ii,j,x,y,z
        real(kind=rk) :: rij(3)
        integer n_c
        type(local_chunk_type) :: c(n_max_local_chunks)

        pot = 0.0_rk
        i = atompnt
        particles: do ii = 1,nlocal+nother
           call local_chunks(P(i),r_cut,0,n_c,c)
           chunks: do j=1,n_c
              do x=c(j)%minx(1),c(j)%maxx(1)
                 do y=c(j)%minx(2),c(j)%maxx(2)
                    do z=c(j)%minx(3),c(j)%maxx(3)
                       surface: if (N(x,y,z)%rock_state==-2.0_rk) then
                          rij = c(j)%xc-real((/x,y,z/),kind=rk)
                          lt_r_cut: if (dot_product(rij,rij)<r_cut_sq) then
                             pot = pot + pair_potential(P(i)%o,rij)
                          end if lt_r_cut
                       end if surface
                    end do
                 end do
              end do
           end do chunks
           if (ii<=nlocal) then
              i = list(i)
           else
              i = i+1
           endif
        end do particles
    end subroutine local_rock_potential_bprh

    !> Returns the energy due to the rock potential according to the
    !> global rock state array \c rock_state for a single particle at
    !> position \c pos with orientation \c ori .
    real(kind=rk) function particle_rock_potential_bprh(rock_state,pos,ori)
        real(kind=rk),intent(in) :: rock_state&
             &(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        real(kind=rk),intent(in) :: pos(3),ori(0:3)
        integer lx(3),lxlo(3),lxhi(3),x,y,z
        real(kind=rk) :: pot,o(3),rij(3)

        o = orientation(ori)

        lxlo = ceiling(pos - r_cut)
        lxhi = floor(pos + r_cut)

        pot = 0.0_rk
        do x=lxlo(1),lxhi(1)
           do y=lxlo(2),lxhi(2)
              do z=lxlo(3),lxhi(3)
                 lx = periodically_wrapped((/x,y,z/))
                 surface: if (rock_state(lx(1),lx(2),lx(3))==-2.0_rk) then
                    rij = pos-real((/x,y,z/),kind=rk)
                    lt_r_cut: if (dot_product(rij,rij)<r_cut_sq) then
                       pot = pot + pair_potential(o,rij)
                    end if lt_r_cut
                 end if surface
              end do
           end do
        end do
        particle_rock_potential_bprh = pot
    end function particle_rock_potential_bprh

    subroutine rock_ft_interaction_bprh(N)
        type(lbe_site),intent(in) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        integer i,j,ii,x,y,z
        real(kind=rk) :: rij(3),rij2sqrt2,rij2,rijo,rijsqrt1,sc,xrijo
        integer n_c
        type(local_chunk_type) :: c(n_max_local_chunks)

        i = atompnt
        particles: do ii = 1,nlocal+nother
           call local_chunks(P(i),r_cut,0,n_c,c)
           chunks: do j=1,n_c
              x_loop: do x=c(j)%minx(1),c(j)%maxx(1)
                 y_loop: do y=c(j)%minx(2),c(j)%maxx(2)
                    z_loop: do z=c(j)%minx(3),c(j)%maxx(3)
                       surface: if (N(x,y,z)%rock_state==-2.0_rk) then
                          rij = c(j)%xc-real((/x,y,z/),kind=rk)
                          rij2 = dot_product(rij,rij)

                          lt_r_cut: if (rij2<r_cut_sq) then
                             rijo = dot_product(rij,P(i)%o)

                             ! \chi \mathbf{r}_{ij}\mathbf{o}
                             xrijo = chi*rijo

                             ! r_{ij}^2\sqrt{\ldots}^2
                             rij2sqrt2 = rij2-rijo*xrijo

                             cutoff: if (rij2sqrt2<sigma_barsq) then
                                ! r_{ij}\sqrt{\ldots}
                                rijsqrt1 = sqrt(rij2sqrt2)

                                ! common scalar expression in force and torque
                                sc = epsilon_2_by_sigma_bar&
                                     &*(1.0_rk/rijsqrt1-inv_sigma_bar)

                                ! force on particle
                                P(i)%f = P(i)%f + sc*(rij-xrijo*P(i)%o)

                                ! torque on particle
                                P(i)%t = P(i)%t&
                                     & + sc*xrijo*cross_product(rij,P(i)%o)
                             end if cutoff
                             if (lubrication)&
                                  & call lubrication_force_rock(rij2,rij,i)
                          end if lt_r_cut
                       end if surface
                    end do z_loop
                 end do y_loop
              end do x_loop
           end do chunks
           if (ii<=nlocal) then
              i = list(i)
           else
              i = i+1
           endif
        end do particles
    end subroutine rock_ft_interaction_bprh

    subroutine setup_rock_bprh
        ! save expensive instructions in the force loop
        r_cut_sq = r_cut*r_cut
        sigma_barsq = (sigma+sigma_orth)**2
        inv_sigma_bar = 1.0_rk/(sigma+sigma_orth)
        chi = 1.0_rk-(sigma+sigma_orth)**2/(sigma+sigma_para)**2
        epsilon_2_by_sigma_bar = 2.0_rk*epsilon*inv_sigma_bar

        ! this check was added due to the dependence on
        ! rc_lubrication_rock below and because
        ! lubrication_force_rock() is not compatible with
        ! polydispersity yet
        if (polydispersity) call error_md("rock='bprh' does not yet "&
             &//"support polydisperse particles---disable polydispersity "&
             &//"or select rock/='bprh'!")

        if (lubrication.and.r_cut<rc_lubrication_rock.and.myrankc==0) &
             &call error_md('lbe_md_rock_bprh_module: r_cut too small! '&
             &//'(sigma+max(sigma_orth,sigma_para)<max(R_orth,R_para)+1/2+2/3 '&
             &//'but lubrication==.true.)')
    end subroutine setup_rock_bprh

    !> returns the potential energy of one ellipsoid with orientations
    !> \c o and one sphere, their distance vector being \c rij
    real(kind=rk) function pair_potential(o,rij)
        real(kind=rk),intent(in) :: o(3),rij(3)
        real(kind=rk) :: Phi,r,r2,rij2,s2

        rij2 = dot_product(rij,rij)

        ! (angle dependent range parameter \sigma)**2
        s2 = sigma_barsq/(1.0_rk-chi*dot_product(rij,o)**2/rij2)

        ! (normalized input for radial symmetric potential \Phi(r))**2
        r2 = rij2/s2

        cutoff: if (r2<1.0_rk) then
           ! normalized input for radial symmetric potential \Phi(r)
           r = sqrt(r2)

           ! radial symmetric repulsive Hooke potential \Phi(r)
           Phi = (1.0_rk-r)**2

           pair_potential = epsilon*Phi
        else cutoff
           pair_potential = 0.0_rk
        end if cutoff
    end function pair_potential

#endif
end module lbe_md_rock_bprh_module
