#include "lbe.h"
!> routines implementing Lees-Edwards boundary conditions also for MD
!> \c interaction=='ladd' particles
module lbe_md_bc_leesedwards_module
#ifdef MD
    use lbe_bdist_module, only: boltz_dist
    use lbe_globals_module, only: border,chunksize,halo_extent,maxpos,minpos&
         &,myrankc,nvecs,tsize
    use lbe_helper_module, only: is_colloid,is_fluid
    use lbe_leesedwards_module, only: le_fractional_offset&
         &,le_find_bottom_neighbours,le_find_top_neighbours,le_mrbuf,le_prbuf&
         &,le_shear_offset,le_shear_omega,le_shear_velocity
    use lbe_md_fluid_ladd_module, only: average_density_r&
         &,local_particle_velocity
    use lbe_md_globals_module
    use lbe_md_helper_module, only: error_md,list_update_add_displacement&
         &,log_msg_md,trigger_list_update
    use lbe_parallel_module, only: ccoords,cdims,comm_cart,nprocs
    use lbe_parms_module, only: boundary_cond,inv_fluid,nx,ny,nz
    use lbe_types_module, only: lbe_site
#ifdef LADD_SURR_RHOF
    use map_module
#endif

    implicit none
    include 'mpif.h'
    private

    !> current (MD substep-precise) z-offset between lower and upper
    !> global x system boundary
    real(kind=rk),save :: shear_offset

    !> allows to enforce having MD Lees-Edwards turned off
    logical,save,public :: no_md_leesedwards=.true.

    !> if \c .true. checkpoints are exactly reproducible at the cost
    !> of somewhat reduced performance and the risk of treating the
    !> particle \c rock_state somewhat different in LE halos than in
    !> the rest of the domain.
    !>
    !> The LE halo \c rock_state at the begin of a time step is
    !> defined by its state in the previous time step plus the changed
    !> done by MD fluid Ladd \c update_rock_state(). When restarting
    !> from a checkpoint, update_rock_state() is applied to a
    !> rock-free halo, so wherever a site is claimed by two particles,
    !> the one with the lower \c uid wins while before also the
    !> history of the site would matter. With this option, the LE halo
    !> \c rock_state is cleared in every time step, so in the LE halo,
    !> the particle with the lower \c uid always wins.
    logical,save,public :: md_leesedwards_reproducible_cp=.false.

    public input_md_leesedwards,md_before_le_adapt_halo_x&
         &,md_leesedwards_adapt_after_communicate&
         &,md_leesedwards_adapt_after_exchange&
         &,md_leesedwards_adapt_after_gather,md_leesedwards_adapt_fluid_rbufs&
#ifdef DEBUG_REPORTMDCOMM
         &,md_leesedwards_debug_reportmdcomm&
#endif
         &,md_leesedwards_init_step&
         &,md_leesedwards_provisional_communication_setup&
         &,md_leesedwards_z_split_down,md_leesedwards_z_split_up&
         &,md_shear_offset,setup_md_leesedwards

contains

    !> adapt a single Lees-Edwards fluid recv buffer according to
    !> possible MD Ladd particles for the following LE interpolation
    !>
    !> The subroutine assumes that the particle rock state is set
    !> correctly in the LE halo planes and that the fluid populations
    !> matter only at fluid sites. It also assumes that the
    !> corresponding MD particles can be accessed via \c uid2i to
    !> obtain their local velocity and possibly also the surrounding
    !> fluid density (with LADD_SURR_RHOF).
    !>
    !> \param[in] lbe_N local chunk of the lattice with halo extent 1
    !> (old LB3D style)
    !>
    !> \param[in,out] rbuf Lees-Edwards recv buffer (one of \c
    !> le_mrbuf and \c le_prbuf ). Dimensions are reals communicated
    !> per site, y-, and z-direction.
    !>
    !> \param[in] x x-position of the corresponding halo layer in \c
    !> lbe_N
    !>
    !> \param[in] dv z-velocity across the LE planes; the sign is
    !> opposite as the correction applied to everything entering at
    !> layer \c x (see above)
    !>
    !> \param[in] dz_frac non-integer part of the z-offset across both
    !> LE planes; always measured in negative z-direction starting at
    !> an integer site position, therefore different for both halos
    subroutine adapt_fluid_rbuf(lbe_N,rbuf,x,dv,dz_frac)
        type(lbe_site),intent(in) :: lbe_N(0:,0:,0:)
        real(kind=rk),intent(inout) :: rbuf(1:,1:,1:)
        integer,intent(in) :: x
        real(kind=rk),intent(in) :: dv,dz_frac
        ! begin and end indexes for red fluid, and the position of the
        ! rock_state in the first dimension of le_mrbuf and le_prbuf
        integer,parameter :: rb=1,re=19,rs=20
        real(kind=rk) :: dist(nvecs),rhof_r,v(3),vpos(3)
        integer :: y,z,uid

        ! don't care about the y- and z-halo here since a halo
        ! exchange in these directions will follow
        do y=1,ny
           do z=1,nz+1
              ! only do something where necessary, that is, where a
              ! particle site in the recv buffer will be used for
              ! interpolation of a fluid site in the halo
              if (is_colloid(rbuf(rs,y+1,z+1))&
                   &.and.any(is_fluid(lbe_N(x,y,z-1:z)%rock_state))) then
                 uid = nint(rbuf(rs,y+1,z+1))
#ifdef LADD_SURR_RHOF
                 rhof_r = P(Mii_map(uid2i,uid))%rhof
#else
                 rhof_r = average_density_r()
#endif
                 vpos = real((/x,y,z/),kind=rk)
                 vpos(3) = vpos(3)-dz_frac
!!!
#ifdef LE_NOPVEL
                 v(1:2) = 0.0_rk
                 v(3) = -0.5_rk*dv
#else
!!!
                 v = local_particle_velocity(vpos,uid)
!!!
#endif
!!!

                 ! add z-velocity offset here to compensate for the
                 ! subtraction to come in le_adapt_halo_x()
                 call boltz_dist(v(1),v(2),v(3)+dv,0.0_rk,0.0_rk,0.0_rk&
                      &,0.0_rk,0.0_rk,0.0_rk,dist)

                 rbuf(rb:re,y+1,z+1) = rhof_r*dist
              end if
           end do
        end do
    end subroutine adapt_fluid_rbuf

    !> determine whether MD Lees-Edwards code is active, if yes check
    !> for conflicting options
    subroutine input_md_leesedwards()
        if (no_md_leesedwards) then
           md_leesedwards = .false.
        else
           md_leesedwards = boundary_cond/=0.and.inv_fluid==5
        end if

        if (.not.md_leesedwards) return

        ! ensure that particles across a LE plane come always from
        ! another process
        if (cdims(1)<2) call error_md('MD Lees-Edwards requires the '&
             &//'x-direction to be decomposed into two or more processes'&
             &//'---set cdx>1 or no_md_leesedwards=.true.!')

        ! ensure that no particle needs to be present in the halo of
        ! the same process twice (once via a normal swap and once via
        ! a LE swap) !!! THIS NEEDS TO BE CHECKED ONCE MORE !!!
        if (rs>=tsize(1)/3.0_rk) call error_md('MD Lees-Edwards requires '&
             &//'rs<nx/3---increase nx or decrease rs!')

        ! Otherwise there would be more than one x-swap per time step
        ! which would require normal swaps and LE swaps to be
        ! interleaved. This would actually be far easier to implement
        ! than something to get rid of the two requirements above.
        if (rs>=chunksize(1)) call error_md('MD Lees-Edwards requires rs to '&
             &//'be smaller than the x-domain size---reduce rs or/and cdx!')

        if (le_shear_omega()/=0.0_rk) call error_md('MD Lees-Edwards was not '&
             &//'tested for shear_omega/=0---check the code carefully, start '&
             &//'with md_leesedwards_adapt_after_exchange()!')

#ifndef NOSURFACTANT
        call error_md('MD Lees-Edwards probably does not take sufficient care '&
             &//'of d_adv. Compile with NOSURFACTANT or run with '&
             &//'no_md_leesedwards==.true.!')
#endif
#ifndef SINGLEFLUID
        call error_md('MD Lees-Edwards does not account for more than the red '&
             &//'fluid component. Compile with SINGLEFLUID or run with '&
             &//'no_md_leesedwards==.true.!')
#endif
    end subroutine input_md_leesedwards

    !> adapt Lees-Edwards fluid halo according to possible MD
    !> particles at the LE planes
    !>
    !> \param[in] lbe_N local chunk of the lattice with halo extent 1
    !> (old LB3D style)
    subroutine md_before_le_adapt_halo_x(lbe_N)
        type(lbe_site),intent(in) :: lbe_N(0:,0:,0:)

        if (md_leesedwards.and.md_initialized) &
             &call md_leesedwards_adapt_fluid_rbufs(lbe_N)
    end subroutine md_before_le_adapt_halo_x

    !> adapt particle positions and velocities assuming they were just
    !> sent by \c borders() or \c communicate() across an LE plane
    !>
    !> \param[in,out] p array containing the particles to adapt
    subroutine md_leesedwards_adapt_after_communicate(p)
        type(md_particle_type),intent(inout) :: p(1:)
        real(kind=rk) :: dvz,dz
        integer :: i

        if (ccoords(1)==0) then
           ! particles came upward
           dz = -shear_offset
           dvz = -le_shear_velocity()
        else
           ! particles came downward
           dz = shear_offset
           dvz = le_shear_velocity()
        end if

        do i=1,size(p)
           p(i)%x(3) = modulo(p(i)%x(3)+dz-0.5_rk,tsize(3))+0.5_rk

           p(i)%v(3) = p(i)%v(3) + dvz
           p(i)%v_fluid_avg(3) = p(i)%v_fluid_avg(3) + dvz
        end do
    end subroutine md_leesedwards_adapt_after_communicate

    !> adapt particle positions and velocities assuming they were just
    !> \c exchange() 'd via an LE plane
    !>
    !> \param[in,out] p array containing the particles to adapt
    subroutine md_leesedwards_adapt_after_exchange(p)
        type(md_particle_type),intent(inout) :: p(1:)
        real(kind=rk) :: dvz,dvz_acc,dz
        integer :: i

        if (ccoords(1)==0) then
           ! particles came upward
           dz = -shear_offset
           dvz = -le_shear_velocity()
        else
           ! particles came downward
           dz = shear_offset
           dvz = le_shear_velocity()
        end if

        ! correct v_fluid_acc as often as v was added to it. If
        ! shear_u should change over time, a small error is made here
        ! as dvz comes from the current timestep and v_fluid_acc maybe
        ! from the previous one. (!!!)
        dvz_acc = dvz*modulo(nt_substep-2,steps_per_lbe_step)

        do i=1,size(p)
           p(i)%x(3) = modulo(p(i)%x(3)+dz-0.5_rk,tsize(3))+0.5_rk

           p(i)%v(3) = p(i)%v(3) + dvz
           p(i)%v_fluid(3) = p(i)%v_fluid(3) + dvz
           p(i)%v_fluid_acc(3) = p(i)%v_fluid_acc(3) + dvz_acc
           p(i)%v_fluid_avg(3) = p(i)%v_fluid_avg(3) + dvz
           p(i)%vnew(3) = p(i)%vnew(3) + dvz
        end do
    end subroutine md_leesedwards_adapt_after_exchange

    !> adapt particle positions and velocities assuming they were just
    !> gathered by \c gather_particles()
    !>
    !> \param[in,out] p array containing the particles to adapt
    !>
    !> \param[in] counts particle counts received from each rank in \c
    !> comm_cart as in \c gather_particles()
    !>
    !> \param[in] counts displacements of particles in \c p as
    !> received from each rank in \c comm_cart in \c
    !> gather_particles()
    !>
    !> The adaption is necessary because \c exchange() is not called
    !> every time step. This allows particles that travelled across
    !> the periodic x-boundary to stay on the previously owning node
    !> from the other side of the system for a few time steps. If they
    !> are gathered naively, they still lack the adaption of z
    !> position and velocity they would normally undergo with the next
    !> \c exchange(). This routine is supposed to perform this
    !> adaption for particles that were gathered from one process but
    !> have a position in the opposite x-boundary layer of the system.
    !>
    !> \warning This routine requires the particle x position to be
    !> checked. Therefore, the mpi type used with \c
    !> gather_particles() must contain \c p(:)%x . This is always the
    !> case for all current used of \c gather_particles() except for
    !> the one used for \c md-cfg output where a check was implemented
    !> that gets activated with \c md_leesedwards.
    subroutine md_leesedwards_adapt_after_gather(p,counts,displs)
        type(md_particle_type),intent(inout) :: p(1:)
        integer,intent(in) :: counts(0:nprocs-1),displs(0:nprocs-1)
        real(kind=rk) :: dvz,dvz_acc,dz
        integer :: i,ierror,r,y,z

        if (myrankc/=0) return

        do y=0,cdims(2)-1
           do z=0,cdims(3)-1
              ! r first loops through ranks at bottom x layer
              call MPI_Cart_Rank(comm_cart,(/0,y,z/),r,ierror)

              ! particle needs to be adapted as if it was exchange()'d
              ! to the top x layer
              dz = shear_offset
              dvz = le_shear_velocity()

              ! correct v_fluid_acc as often as v was added to it. If
              ! shear_u should change over time, a small error is made
              ! here as dvz comes from the current timestep and
              ! v_fluid_acc maybe from the previous one. (!!!)
              dvz_acc = dvz*modulo(nt_substep-2,steps_per_lbe_step)

              ! loop through particles sent from rank r
              do i=displs(r)+1,displs(r)+counts(r)
                 if (p(i)%x(1)>=maxpos(1)-1.0_rk) then
                    p(i)%x(3) = modulo(p(i)%x(3)+dz-0.5_rk,tsize(3))+0.5_rk

                    p(i)%v(3) = p(i)%v(3) + dvz
                    p(i)%v_fluid(3) = p(i)%v_fluid(3) + dvz
                    p(i)%v_fluid_acc(3) = p(i)%v_fluid_acc(3) + dvz_acc
                    p(i)%v_fluid_avg(3) = p(i)%v_fluid_avg(3) + dvz
                    p(i)%vnew(3) = p(i)%vnew(3) + dvz
                 end if
              end do

              ! r then loops through ranks at top x layer
              call MPI_Cart_Rank(comm_cart,(/cdims(1)-1,y,z/),r,ierror)

              ! particle needs to be adapted as if it was exchange()'d
              ! to the bottom x layer
              dz = -shear_offset
              dvz = -le_shear_velocity()

              ! correct v_fluid_acc as often as v was added to it. If
              ! shear_u should change over time, a small error is made
              ! here as dvz comes from the current timestep and
              ! v_fluid_acc maybe from the previous one. (!!!)
              dvz_acc = dvz*modulo(nt_substep-2,steps_per_lbe_step)

              ! loop through particles sent from rank r
              do i=displs(r)+1,displs(r)+counts(r)
                 if (p(i)%x(1)<minpos(1)+1.0_rk) then
                    p(i)%x(3) = modulo(p(i)%x(3)+dz-0.5_rk,tsize(3))+0.5_rk

                    p(i)%v(3) = p(i)%v(3) + dvz
                    p(i)%v_fluid(3) = p(i)%v_fluid(3) + dvz
                    p(i)%v_fluid_acc(3) = p(i)%v_fluid_acc(3) + dvz_acc
                    p(i)%v_fluid_avg(3) = p(i)%v_fluid_avg(3) + dvz
                    p(i)%vnew(3) = p(i)%vnew(3) + dvz
                 end if
              end do
           end do
        end do
    end subroutine md_leesedwards_adapt_after_gather

    !> prepare the fluid recv buffers \c le_mrbuf and \c le_prbuf from
    !> \c lbe_leesedwards_module at the LE planes according to
    !> possible MD Ladd particles at the interface in a way that the
    !> following LE interpolation leads to reasonable fluid
    !> populations in the LE halos
    !>
    !> \param[in] lbe_N local chunk of the lattice with halo extent 1
    !> (old LB3D style)
    subroutine md_leesedwards_adapt_fluid_rbufs(lbe_N)
        type(lbe_site),intent(in) :: lbe_N(0:,0:,0:)

        if (interaction/='ladd') return

        if (ccoords(1)==0) then
           call adapt_fluid_rbuf(lbe_N,le_prbuf,0,le_shear_velocity()&
                &,le_fractional_offset())
        else if (ccoords(1)==cdims(1)-1) then
           call adapt_fluid_rbuf(lbe_N,le_mrbuf,nx+1,-le_shear_velocity()&
                &,1.0_rk-le_fractional_offset())
        end if
    end subroutine md_leesedwards_adapt_fluid_rbufs

#ifdef DEBUG_REPORTMDCOMM
    !> reports the communication setup and z send boundaries for a
    !> simulated period of shear
    subroutine md_leesedwards_debug_reportmdcomm()
        integer,parameter :: n_steps=1000
        real(kind=rk) :: prev_shear_offset
        logical :: lox,hix
        integer :: i,ierror,p

        prev_shear_offset = shear_offset

        do i=1,n_steps
           call update_communication_topology(shear_offset)

           lox = ccoords(1)==0
           hix = ccoords(1)==cdims(1)-1

           do p=0,nprocs
              call MPI_Barrier(MPI_COMM_WORLD,ierror)

              if (myrankc==p) then
                 if (lox.or.hix) write (unit=6,advance='no',fmt=&
                      &'("so= ",ES15.8," r= ",I0," (cx cz)=( ",2(I0,X),") ")') &
                      &shear_offset,myrankc,ccoords((/1,3/))
                 if (lox) write (unit=6,fmt='("DN (slo shi)=( ",2(I0,X)'&
                      &//',") (rlo rhi)=( ",2(I0,X)'&
                      &//',") z_lo=[ ",2(ES15.8,X)," ] z_hi=[ ",2(ES15.8,X)'&
                      &//'," ]")')&
                      & sdirs(1,(/DIR_LE_LO_Z_DOWN,DIR_LE_HI_Z_DOWN/))%sproc&
                      &,sdirs(1,(/DIR_LE_LO_Z_UP,DIR_LE_HI_Z_UP/))%rproc&
                      &,sdirs(1,DIR_LE_LO_Z_DOWN)%blo_z&
                      &,sdirs(1,DIR_LE_LO_Z_DOWN)%bhi_z&
                      &,sdirs(1,DIR_LE_HI_Z_DOWN)%blo_z&
                      &,sdirs(1,DIR_LE_HI_Z_DOWN)%bhi_z
                 if (hix) write (unit=6,fmt='("UP (slo shi)=( ",2(I0,X)'&
                      &//',") (rlo rhi)=( ",2(I0,X)'&
                      &//',") z_lo=[ ",2(ES15.8,X)," ] z_hi=[ ",2(ES15.8,X)'&
                      &//'," ]")')&
                      & sdirs(1,(/DIR_LE_LO_Z_UP,DIR_LE_HI_Z_UP/))%sproc&
                      &,sdirs(1,(/DIR_LE_LO_Z_DOWN,DIR_LE_HI_Z_DOWN/))%rproc&
                      &,sdirs(1,DIR_LE_LO_Z_UP)%blo_z&
                      &,sdirs(1,DIR_LE_LO_Z_UP)%bhi_z&
                      &,sdirs(1,DIR_LE_HI_Z_UP)%blo_z&
                      &,sdirs(1,DIR_LE_HI_Z_UP)%bhi_z
              end if
           end do

           call MPI_Barrier(MPI_COMM_WORLD,ierror)

           shear_offset = shear_offset + le_shear_velocity()
        end do

        shear_offset = prev_shear_offset
    end subroutine md_leesedwards_debug_reportmdcomm
#endif

    !> initialization at the beginning of each new MD (sub-)step
    !>
    !> \param[in,out] whole_N local lattice chunk with halo of extent
    !> \c halo_extent
    subroutine md_leesedwards_init_step(whole_N)
        type(lbe_site),intent(inout),optional :: &
             &whole_N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        real(kind=rk) :: displ_per_substep

        displ_per_substep = le_shear_velocity()/real(steps_per_lbe_step,kind=rk)

        ! extrapolate shear_offset for current MD substep
        shear_offset = le_shear_offset() &
             &+ displ_per_substep*real(nt_substep-1,kind=rk)

        call update_communication_topology(shear_offset)
        call trigger_list_update_on_topology_change()
        call list_update_add_displacement(displ_per_substep)

        if (md_leesedwards_reproducible_cp.and.nt_substep==1) then
           ! enforce checkpoint reproducibility by clearing the LE
           ! halo rockstate each time step because when restoring,
           ! there also is no LE halo rock state initially
           if (ccoords(1)==0) whole_N(0,:,:)%rock_state = 0.0_rk
           if (ccoords(1)==cdims(1)-1) whole_N(nx+1,:,:)%rock_state = 0.0_rk
        end if
    end subroutine md_leesedwards_init_step

    !> setup LE communication in a consistent way without caring about
    !> the actual \c shear_offset
    subroutine md_leesedwards_provisional_communication_setup()
        call update_communication_topology(0.0_rk)
    end subroutine md_leesedwards_provisional_communication_setup

    !> finds the z position where the two process domains touch that
    !> lie downward across the LE plane
    !>
    !> \returns obtained z position
    pure function md_leesedwards_z_split_down()
        real(kind=rk) :: md_leesedwards_z_split_down

        md_leesedwards_z_split_down = border(2,3)-modulo(shear_offset,chunksize(3))
    end function md_leesedwards_z_split_down

    !> finds the z position where the two process domains touch that
    !> lie upward across the LE plane
    !>
    !> \returns obtained z position
    pure function md_leesedwards_z_split_up()
        real(kind=rk) :: md_leesedwards_z_split_up

        md_leesedwards_z_split_up = border(1,3)+modulo(shear_offset,chunksize(3))
    end function md_leesedwards_z_split_up

    !> provides access to Lees-Edwards \c shear offset including
    !> possible extrapolation for MD substeps
    !>
    !> \returns spatial z-offset between lower and upper global x
    !> boundary of the system
    pure function md_shear_offset()
        real(kind=rk) :: md_shear_offset

        md_shear_offset = shear_offset
    end function md_shear_offset

    !> setup of MD Lees-Edwards code
    !>
    !> \param[in,out] whole_N local lattice chunk with halo of extent
    !> \c halo_extent
    !>
    !> The basic communication scheme must have been set up by \c
    !> setup_parallel() before, what comes here is the linking the two
    !> planes according to the current \c shear_offset.
    subroutine setup_md_leesedwards(whole_N)
        type(lbe_site),intent(inout) :: &
             &whole_N(1-halo_extent:,1-halo_extent:,1-halo_extent:)

        if (.not.md_leesedwards) then
           call log_msg_md(&
                &'Lees-Edwards boundaries for MD particles are inactive.')
           return
        end if

        call log_msg_md('Lees-Edwards boundaries for MD particles are active.')

        call md_leesedwards_init_step(whole_N)
    end subroutine setup_md_leesedwards

    !> sets up send and recv partners and z send boundaries according
    !> to a given shear offset
    !>
    !> \param[in] so assumed z-offset across the Lees-Edwards planes
    subroutine update_communication_topology(so)
        real(kind=rk),intent(in) :: so
        integer hiz_proc,loz_proc
        real(kind=rk) :: h,s
        logical hix,lox

        h = huge(h)                  ! just for convencienc
        lox = ccoords(1)==0          ! are we at lower...
        hix = ccoords(1)==cdims(1)-1 ! or upper Lees-Edwards plane?

        if (lox) then
           s = md_leesedwards_z_split_down()
           call le_find_bottom_neighbours(so,loz_proc,hiz_proc)

           sdirs(1,DIR_LE_LO_Z_DOWN)%sproc = loz_proc
           sdirs(1,DIR_LE_HI_Z_UP)%rproc = loz_proc
           if (loz_proc/=hiz_proc) then
              sdirs(1,DIR_LE_HI_Z_DOWN)%sproc = hiz_proc
              sdirs(1,DIR_LE_LO_Z_UP)%rproc = hiz_proc

              sdirs(1,DIR_LE_LO_Z_DOWN)%blo_z = -h
              sdirs(1,DIR_LE_LO_Z_DOWN)%bhi_z = s
              sdirs(1,DIR_LE_HI_Z_DOWN)%blo_z = s
              sdirs(1,DIR_LE_HI_Z_DOWN)%bhi_z = h
           else
              ! this can happen only for cdims(3)==1, so one send is
              ! enough but it needs to cover the full z range
              sdirs(1,DIR_LE_HI_Z_DOWN)%sproc = MPI_PROC_NULL
              sdirs(1,DIR_LE_LO_Z_UP)%rproc = MPI_PROC_NULL

              sdirs(1,DIR_LE_LO_Z_DOWN)%blo_z = -h
              sdirs(1,DIR_LE_LO_Z_DOWN)%bhi_z = h
           end if
        end if

        if (hix) then
           s = md_leesedwards_z_split_up()
           call le_find_top_neighbours(so,loz_proc,hiz_proc)

           sdirs(1,DIR_LE_HI_Z_UP)%sproc = hiz_proc
           sdirs(1,DIR_LE_LO_Z_DOWN)%rproc = hiz_proc
           if (loz_proc/=hiz_proc) then
              sdirs(1,DIR_LE_LO_Z_UP)%sproc = loz_proc
              sdirs(1,DIR_LE_HI_Z_DOWN)%rproc = loz_proc

              sdirs(1,DIR_LE_LO_Z_UP)%blo_z = -h
              sdirs(1,DIR_LE_LO_Z_UP)%bhi_z = s
              sdirs(1,DIR_LE_HI_Z_UP)%blo_z = s
              sdirs(1,DIR_LE_HI_Z_UP)%bhi_z = h
           else
              ! this can happen only for cdims(3)==1, so one send is
              ! enough but it needs to cover the full z range
              sdirs(1,DIR_LE_LO_Z_UP)%sproc = MPI_PROC_NULL
              sdirs(1,DIR_LE_HI_Z_DOWN)%rproc = MPI_PROC_NULL

              sdirs(1,DIR_LE_HI_Z_UP)%blo_z = -h
              sdirs(1,DIR_LE_HI_Z_UP)%bhi_z = h
           end if
        end if
    end subroutine update_communication_topology

    !> trigger a neighbor list update if the LE communication topology
    !> should has changed since the last call
    subroutine trigger_list_update_on_topology_change()
        integer,save :: prev_topo(4)=MPI_PROC_NULL
        integer :: cur_topo(4)

        cur_topo = (/&
             &sdirs(1,DIR_LE_LO_Z_DOWN)%sproc&
             &,sdirs(1,DIR_LE_LO_Z_UP)%sproc&
             &,sdirs(1,DIR_LE_HI_Z_DOWN)%sproc&
             &,sdirs(1,DIR_LE_HI_Z_UP)%sproc&
             &/)
        if (any(cur_topo/=prev_topo)) call trigger_list_update()
        prev_topo = cur_topo
    end subroutine trigger_list_update_on_topology_change

#endif
end module lbe_md_bc_leesedwards_module
