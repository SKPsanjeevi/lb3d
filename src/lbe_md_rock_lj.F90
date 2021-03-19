#include "lbe.h"

!> Lennard-Jones potential between particles and rock sites at the
!> surface of walls
module lbe_md_rock_lj_module
#ifdef MD
    use lbe_bc_module, only: periodically_wrapped
    use lbe_globals_module, only: halo_extent,tsize, input_dfile_unit,myrankc
    use lbe_helper_module, only: local_coordinates
    use lbe_log_module
    use lbe_md_boundary_condition_module, only: local_chunks,local_chunk_type&
         &,n_max_local_chunks
    use lbe_md_fluid_ladd_module, only: lubrication_force_rock&
         &,rc_lubrication_rock
    use lbe_md_fluid_ladd_parms_module, only: lubrication
    use lbe_md_globals_module
    use lbe_md_helper_module, only: log_msg_md,error_md,log_msg_md_hdr
    use lbe_parallel_module, only: comm_cart
    use lbe_parms_module, only: inp_file,arg_input_dfile_set,arg_input_dfile
    use lbe_types_module, only: lbe_site

    implicit none
    include 'mpif.h'
    private

    public input_rock_lj,local_rock_potential_lj,particle_rock_potential_lj&
         &,rock_f_interaction_lj,setup_rock_lj

    ! lennard jones parameters
    real(kind=rk),save :: epsilon=0.001_rk
    real(kind=rk),save :: sigma=5.0_rk

    real(kind=rk),save :: r_cut         ! lennard jones cutoff radius
    real(kind=rk),save :: r_cut_sq      ! r_cut**2
    real(kind=rk),save :: sigma_sq      ! sigma**2
    real(kind=rk),save :: epsilon4      ! epsilon*4
    real(kind=rk),save :: epsilon48     ! epsilon*48

    namelist /md_rock_lj/ epsilon,sigma

contains

    !> read namelist  /md_rock_lj/  from input file
    subroutine input_rock_lj
        integer ierror

        call log_msg_md_hdr("Reading MD R LJ input")

        ! These features are essential for this rock module. They are
        ! enabled here because in  setup_rock_bprh()  it would be too late.
        collect_forces = .true.

        if (myrankc.eq.0) then
           open (unit=md_input_file_unit,file=trim(inp_file)//'.md',err=100)
           read (unit=md_input_file_unit,nml=md_rock_lj,err=100)
           close (unit=md_input_file_unit,err=100)
           !call log_msg_md('read /md_rock_lj/ from file "'//trim(inp_file)&
           !     &//'.md"',.false.)
           !write (6,nml=md_rock_lj)
        end if

        if ( arg_input_dfile_set ) then
          call log_msg_md("  Getting differential input...")
          open(UNIT = input_dfile_unit, FILE = arg_input_dfile, STATUS = 'UNKNOWN')
          read(UNIT = input_dfile_unit, NML = md_rock_lj, IOSTAT = ierror)
          if (ierror .ne. 0) then
              call log_msg_md("    WARNING: Differential namelist not found or errors encountered.")
          end if
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

        r_cut = 2.5_rk*sigma

        return
100     continue
        call error_md('Error reading md input file "'//trim(inp_file)//'.md"')
    end subroutine input_rock_lj

    subroutine local_rock_potential_lj(N,pot)
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
                             pot = pot + pair_potential(rij)
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
    end subroutine local_rock_potential_lj

    !> Returns the energy due to the rock potential according to the
    !> global rock state array \c rock_state for a single particle at
    !> position \c pos.
    real(kind=rk) function particle_rock_potential_lj(rock_state,pos)
        real(kind=rk),intent(in) :: rock_state&
             &(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        real(kind=rk),intent(in) :: pos(3)
        integer lx(3),lxlo(3),lxhi(3),x,y,z
        real(kind=rk) :: pot,rij(3)

        lxlo(:) = ceiling(pos(:) - r_cut)
        lxhi(:) = floor(pos(:) + r_cut)

        pot = 0.0_rk
        do x=lxlo(1),lxhi(1)
           do y=lxlo(2),lxhi(2)
              do z=lxlo(3),lxhi(3)
                 lx = periodically_wrapped((/x,y,z/))
                 surface: if (rock_state(lx(1),lx(2),lx(3))==-2.0_rk) then
                    rij = pos-real((/x,y,z/),kind=rk)
                    lt_r_cut: if (dot_product(rij,rij)<r_cut_sq) then
                       pot = pot + pair_potential(rij)
                    end if lt_r_cut
                 end if surface
              end do
           end do
        end do
        particle_rock_potential_lj = pot
    end function particle_rock_potential_lj

    subroutine rock_f_interaction_lj(N)
        type(lbe_site),intent(in) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        integer i,ii,j,x,y,z
        real(kind=rk) :: dist(3),dist_sq,sr2,sr6
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
                          dist = c(j)%xc-real((/x,y,z/),kind=rk)
                          dist_sq = dot_product(dist,dist)
                          cutoff: if (dist_sq<r_cut_sq) then
                             sr2 = sigma_sq/dist_sq
                             sr6 = sr2*sr2*sr2
                             P(i)%f(:) = P(i)%f(:) + &
                                  &dist(:)*epsilon48*sr6*(sr6-0.5_rk)/dist_sq

                             if (lubrication)&
                                  & call lubrication_force_rock(dist_sq,dist,i)
                          end if cutoff
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
    end subroutine rock_f_interaction_lj

    subroutine setup_rock_lj
        ! save expensive instructions in the force loop
        r_cut_sq = r_cut*r_cut
        sigma_sq = sigma*sigma
        epsilon4 = epsilon*4.0_rk
        epsilon48 = epsilon*48.0_rk

        ! this check was added due to the dependence on
        ! rc_lubrication_rock below and because
        ! lubrication_force_rock() is not compatible with
        ! polydispersity yet
        if (polydispersity) call error_md("rock='lj' does not yet "&
             &//"support polydisperse particles---disable polydispersity "&
             &//"or select rock/='lj'!")

        if (lubrication.and.r_cut<rc_lubrication_rock) call error_md(&
             &'r_cut too small! '&
             &//'(r_cut<rc_lubrication_rock but lubrication==.true.')
    end subroutine setup_rock_lj

    !> returns the potential energy of two particles at a distance \c rij
    real(kind=rk) function pair_potential(rij)
        real(kind=rk),intent(in) :: rij(3)
        real(kind=rk) sr2,sr6,sr12

        sr2 = sigma_sq/dot_product(rij,rij)
        sr6 = sr2*sr2*sr2
        sr12 = sr6*sr6
        pair_potential = epsilon4*(sr12-sr6)
    end function pair_potential

#endif
end module lbe_md_rock_lj_module
