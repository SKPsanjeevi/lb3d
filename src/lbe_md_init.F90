#include "lbe.h"

!> initialization of molecular dynamics part
module lbe_md_init_module
#ifdef MD
  use lbe_globals_module, only: border,chunksize,halo_extent,maxpos,minpos,nd&
       &,pi,tsize,myrankc
  use lbe_helper_module, only: local_coordinates
  use lbe_init_rock_module, only: inside_rock_global,rock_is_present
  use lbe_log_module
  use lbe_md_bc_leesedwards_module, only:&
       & md_leesedwards_provisional_communication_setup
  use lbe_md_boundary_condition_module, only: md_clip_periodic,rdstsq
  use lbe_md_debug_module, only: debug_check_send_bounds
  use lbe_md_fluid_ladd_module, only: avg_particle_volume_fluid_ladd&
       &,init_radii_normal,init_radii_random,init_radii_volume
  use lbe_md_fluid_ladd_parms_module, only: initial_radii,R_orth,R_para
  use lbe_md_globals_module
  use lbe_md_helper_module, only: build_collect_mpitypes,build_particle_mpitype&
       &,count_particles_all,inside_rock,quaternion,space_to_body_matrix&
       &,log_msg_md,maximum_particle_radius,error_md
  use lbe_md_memory_module, only: boost_npmax
  use lbe_md_module, only: borders,calculate_rotations_s,setup_list,temperature, calculate_all_orientations
  use lbe_md_output_module, only: dump_potentials
  use lbe_md_parallel_module, only: scatter_particles
  use lbe_md_rand_walk_module, only: uni_dist
  use lbe_md_rock_module, only: particle_rock_potential
  use lbe_parallel_module, only: calculate_displacements,check_allocate,ccoords&
       &,cdims,comm_cart,gather_rock_state,get_send_partner,nnprocs,nprocs&
       &,owning_rankc,start,tnx,tny,tnz
  use lbe_parms_module, only: boundary,inv_fluid,n_iteration&
       &,drop_xshift,drop_yshift,drop_zshift,drop_xcut,drop_ycut,drop_zcut,fr1&
       &,fr2,nx,ny,nz
  use lbe_timer_module, only: register_timer
  use lbe_types_module, only: lbe_site
  use luxury_rng_module, only: ranlux
  use map_module, only: Mii_init
#ifdef DEBUG_REPORTMDCOMM
  use lbe_md_bc_leesedwards_module, only: md_leesedwards_debug_reportmdcomm
  use lbe_md_debug_module, only: debug_reportmdcomm
#endif
#ifdef FORCECOMPONENT
 use lbe_md_output_module, only: dump_force_components
#endif 
  
  implicit none
  include 'mpif.h'
  private
  public setup_general,setup_parallel,setup_particles,setup_uid2i

  !> path to file to read initial configuration from (when required)
  character(len=1024),save,public :: init_file='init.cfg'
  !> \{
  !> \name min/max position for particle placement
  !>
  !> This will be clipped to \c minpos/maxpos. Every particle
  !> consumes a space of \c alat in every direction, so for \c
  !> all(x0hi==x0lo+alat) only 1 particle is placed instead of 8.
  real(kind=rk),save,public :: x0lo(3)=(/0.0_rk,0.0_rk,0.0_rk/)
  real(kind=rk),save,public :: x0hi(3)=(/-1.0_rk,-1.0_rk,-1.0_rk/)
  !> \}

  !> minimum particle center distance for \c initial_placing=='random'
  real(kind=rk),save,public :: min_dist=6.0_rk

  !> particle volume concentration for \c initial_placing=='random',
  !> negative \c phi means that as many particles are placed as
  !> would fit into a \c sc lattice with \c alat.
  real(kind=rk),save,public :: phi=-1.0_rk

  !> type of orientation initialization
  character(len=32),save,public :: initial_orientations='q0'
  !> initial orientation for all particles
  real(kind=rk),save,public :: q0(0:3)=(/0.0_rk,0.0_rk,0.0_rk,1.0_rk/)

  !> type of velocity initialization
  character(len=32),save,public :: initial_velocities='v0'
  real(kind=rk),save,public :: v0(3)=(/0.0_rk,0.0_rk,0.0_rk/) !< initial v

  !> type of angular velocity initialization
  character(len=32),save,public :: initial_rotations='w0'
  real(kind=rk),save,public :: w0(3)=(/0.0_rk,0.0_rk,0.0_rk/) !< initial w

  !> initial temperature (reduced units, this means: temp0=k_B*T/mass)
  real(kind=rk),save,public :: temp0=1.0e-4_rk

contains

    subroutine setup_general()
        integer ::d,k

        nt_substep = steps_per_lbe_step

        minx = minpos(1)
        miny = minpos(2)
        minz = minpos(3)

        maxx = maxpos(1)
        maxy = maxpos(2)
        maxz = maxpos(3)

        if (rs<rc) call error_md('rs  must not be smaller than  rc .')

        cutsq1 = rc*rc
        cutsq2 = rs*rs
        tsize_mrc = tsize - rc
        tsize_mrs = tsize - rs

        call register_timer('MD:Force',ti_md_force)
        call register_timer('MD:Rock',ti_md_rock)
        call register_timer('MD:Neigh',ti_md_neigh)
        call register_timer('MD:Comm',ti_md_comm)
        call register_timer('MD:Dump',ti_md_dump)
        call register_timer('MD:Fluid:Other',ti_md_fluid)

        mstat = 0
        mneigh = 0

        nlocalmax = 0
        nothermax = 0
        neighmax = 0
        nslistmax = 0
        nexcmax = 0
        nswpmax = 0

        ! The initialization value MPI_DATATYPE_NULL serves as flag to
        ! avoid one of sdirs(:,:)%s(:)%s_mpitype being freed
        ! before its first commitment.
        do k=1,3
           do d=1,n_dir_max
              sdirs(k,d)%dim = k
              sdirs(k,d)%dir = d
              sdirs(k,d)%s%s_mpitype = MPI_DATATYPE_NULL

              ! just initialize z send bounds with reasonable values,
              ! this is not really necessary. If these are requried,
              ! they are initialized separately anyway.
              sdirs(k,d)%blo_z = -huge(sdirs(k,d)%blo_z)
              sdirs(k,d)%bhi_z = huge(sdirs(k,d)%bhi_z)
           end do
        end do

        ! simplified initialization of ladd particle parameters
        if (0==myrankc.and.any((/mass,inertia_orth,inertia_para/)<0.0_rk)&
             &.and.interaction/='ladd') call error_md&
             &('Simplified initialization of particle parameters requires '&
             &//'interaction==''ladd''. Set mass, inertia_orth, and '&
             &//'inertia_para manually!')

        polydisperse: if (polydispersity) then
           call log_msg_md('Polydispersity is enabled. Masses and moments of '&
                &//'inertia are computed for every particle based on rho. '&
                &//'mass, inertia_orth, and inertia_para from the input file '&
                &//'are ignored.')
        else polydisperse
           if (mass<0.0_rk) then
              mass = 4.0_rk*pi*R_orth**2*R_para*rho/3.0_rk
              write(msgstr,"('Initialized mass = ',ES15.8)") mass
              call log_msg_md(msgstr)
           end if
           if (inertia_orth<0.0_rk) then
              inertia_orth = (R_orth**2+R_para**2)*mass/5.0_rk
              write(msgstr,"('Initialized inertia_orth = ',ES15.8)") &
                   &inertia_orth
              call log_msg_md(msgstr)
           end if
           if (inertia_para<0.0_rk) then
              inertia_para = 2.0_rk*R_orth**2*mass/5.0_rk
              write(msgstr,"('Initialized inertia_para = ',ES15.8)") &
                   &inertia_para
              call log_msg_md(msgstr)
           end if
        end if polydisperse

        if (inv_fluid==12.and.interaction=='ladd'.and&
             &.any(reenter_rdm_max.or.reenter_rdm_min)) &
             &call error_md('for inv_fluid=12, we did not find an easy way to '&
             &//'set the density in the node occupied by the particles before '&
             &//'the rdm boundary conditions so pleaso dont use together ladd '&
             &//'and inv_fluid=12')

        if (polydispersity.and.any(reenter_rdm_min.or.reenter_rdm_max)) &
             &call error_md("random reentering does not "&
             &//"support polydisperse particles yet---disable polydispersity "&
             &//"or set reenter_rdm_min=.false. and reenter_rdm_max=.false.!")

!#ifndef SINGLEFLUID
        ! this is because near_particle_edge() was not adjusted to
        ! work with different radii for each particle yet
        ! Adjusted by Xie, so it should work now
!        if (polydispersity) call error_md("multi-component particle-fluid "&
!             &//"coupling does not work with particle polydispersity yet"&
!             &//"---disable polydispersity or compile with SINGLEFLUID!")
!#endif

#ifdef PARTICLESTRESS
        if (polydispersity) call error_md("PARTICLESTRESS does not "&
             &//"support polydisperse particles yet---disable polydispersity "&
             &//"or compile without PARTICLESTRESS!")
#endif

#ifdef VARTAU
#error "MD particle-fluid interactions do not support VARTAU yet---compile without MD or without VARTAU!"
#endif
    end subroutine setup_general

    !> calculates the number of mutual communication steps between
    !> neighboring pairs of cpus that are necessary for md
    !> communication in a given dimension and (pseudo-)direction
    !>
    !> \param[in] curdim Cartesian dimension (1-3 representing x,y,z)
    !>
    !> \param[in] dir (pseudo-)direction
    !>
    !> \returns number of swaps
    !>
    !> In the case of periodic_inflow, directions 1 and 2 are still
    !> used for swapping within the periodic sub-volume.
    integer function get_sends_required(curdim,dir)
        integer,intent(in) :: curdim,dir
        integer cdims_periodic_z,fnc,n
        real(kind=rk) need,rfnz,rlpz
        logical hix,lox

        if (boundary/='periodic_inflow'.or.curdim<3) then
           if (dir==DIR_NONPERIODIC_DOWN.or.dir==DIR_NONPERIODIC_UP) then
              ! There are no non-periodic sends in this case.
              get_sends_required = 0
           else if (cdims(curdim)==1) then
              ! don't exchange if only 1 box in a dimension
              get_sends_required = 0
           else if (md_leesedwards.and.curdim==1) then
              ! special treatment for Lees-Edwards x-communication
              lox = ccoords(1)==0 ! are we at lower...
              hix = ccoords(1)==cdims(1)-1 ! or upper Lees-Edwards plane?
              n = rs/chunksize(curdim)+1 ! natural amount of sends due to rs
              select case (dir)
              case (DIR_DOWN)
                 if (lox) then
                    get_sends_required = 0
                 else
                    get_sends_required = n
                 end if
              case (DIR_UP)
                 if (hix) then
                    get_sends_required = 0
                 else
                    get_sends_required = n
                 end if
              case (DIR_LE_LO_Z_DOWN)
                 if (lox) then
                    get_sends_required = n
                 else
                    get_sends_required = 0
                 end if
              case (DIR_LE_LO_Z_UP)
                 ! skip one of the LE sends if there is just one
                 ! process in z direction
                 if (hix.and.cdims(3)>1) then
                    get_sends_required = n
                 else
                    get_sends_required = 0
                 end if
              case (DIR_LE_HI_Z_DOWN)
                 ! skip one of the LE sends if there is just one
                 ! process in z direction
                 if (lox.and.cdims(3)>1) then
                    get_sends_required = n
                 else
                    get_sends_required = 0
                 end if
              case (DIR_LE_HI_Z_UP)
                 if (hix) then
                    get_sends_required = n
                 else
                    get_sends_required = 0
                 end if
              end select
           else if (any(dir==(/DIR_LE_LO_Z_DOWN,DIR_LE_LO_Z_UP&
                &,DIR_LE_HI_Z_DOWN,DIR_LE_HI_Z_UP/))) then
              ! ...but no Lees-Edwards sends without md_leesedwards
              ! and in dimensions different from x
              get_sends_required = 0
           else
              ! this is the default case: normal communication in the
              ! bulk and for directions with simply periodic
              ! boundaries
              n = rs/chunksize(curdim)+1
              if (2*n>cdims(curdim)) then
                 ! don't exchange more than 1/2 way over (e.g. 3 boxes
                 ! away when cdims = 5). This case cannot happen for
                 ! md_leesedwards x-communication as there we always
                 ! require rs<tsize(3)/3, therefore it is omitted
                 ! above.
                 get_sends_required = n-1
              else
                 get_sends_required = n
              end if
           end if

!!$        else if (dir<3) then
!!$           rlpz = real(last_periodic_z,kind=rk)
!!$           ! dir is one of {1,2} --- z-send with BC periodic_inflow
!!$           ! within periodic sub-volume
!!$
!!$           if (start(3)>last_periodic_z) then
!!$              get_sends_required = 0
!!$           else if (last_periodic_z<=nz) then
!!$              ! don't exchange if only 1 box in a dimension
!!$              get_sends_required = 0
!!$           else
!!$              cdims_periodic_z = ceiling(real(last_periodic_z,kind=rk)/nz)
!!$
!!$              if (mod(last_periodic_z,nz)==0) then
!!$                 ! periodic sub-volume ends with a domain boundary -
!!$                 ! this is easy
!!$
!!$                 get_sends_required = ceiling(rs/chunksize(3))
!!$
!!$              else ! sub-volume boundary cuts a domain
!!$
!!$                 dir_if: if (dir==1) then ! negative z
!!$
!!$                    ! never send further than 50% of the system
!!$                    ! without the target node
!!$                    need = min(rs&
!!$                         &,0.5_rk*real(last_periodic_z-&
!!$                         &calculate_chunksize_pi_z(ccoords(3)-1),kind=rk))
!!$
!!$                    n=0
!!$                    do
!!$                       n = n+1
!!$                       if (real(n_boxes_pi_upward_z(n),kind=rk)>=need) exit
!!$                    end do
!!$                    get_sends_required = n
!!$
!!$                 else dir_if    ! positive z (dir==2)
!!$
!!$                    ! never send further than 50% of the system
!!$                    ! without the target node
!!$                    need = min(rs&
!!$                         &,0.5_rk*real(last_periodic_z-&
!!$                         &calculate_chunksize_pi_z(ccoords(3)+1),kind=rk))
!!$
!!$                    n=0
!!$                    do
!!$                       n = n+1
!!$                       if (real(n_boxes_pi_downward_z(n),kind=rk)>=need) exit
!!$                    end do
!!$                    get_sends_required = n
!!$                 end if dir_if
!!$              end if
!!$           end if
!!$        else
!!$           ! dir is one of {3,4} --- z-swap with BC periodic_inflow
!!$           ! within non-periodic sub-volume
!!$           if (border(2,3)<enslavement_threshold) then
!!$              get_sends_required = 0
!!$           else
!!$              ! first cpu that might own non-periodic particles
!!$              fnc = floor((enslavement_threshold-0.5_rk)/real(nz,kind=rk))
!!$
!!$              ! Never send more often in one direction than there are
!!$              ! cpus in the other one plus the local cpu since no new
!!$              ! data would be available for these sends.
!!$
!!$              if (dir==3) then  ! negative z
!!$                 if (ccoords(3)==fnc) then
!!$                    get_sends_required = 0
!!$                 else
!!$                    get_sends_required = &
!!$                         &min(ceiling(rs/chunksize(3)),cdims(3)-ccoords(3))
!!$                 end if
!!$              else              ! positive z (dir==4)
!!$                 if (ccoords(3)==cdims(3)-1) then
!!$                    get_sends_required = 0
!!$                 else
!!$                    get_sends_required = &
!!$                         &min(ceiling(rs/chunksize(3)),1+ccoords(3)-fnc)
!!$                 end if
!!$              end if
!!$           end if
        end if
    end function get_sends_required

    !> finds the rank of the communication partner to receive from for a given
    !> dimension and (pseudo-)direction
    !>
    !> \param[in] curdim Cartesian dimension (1-3 representing x,y,z)
    !>
    !> \param[in] dir (pseudo-)direction
    !>
    !> \returns rank in \c comm_cart
    !>
    !> Note that \c dir refers to the direction along that the data is
    !> transfered, it thus is the same as on the sending site.
    !>
    !> Even in the case of \c periodic_inflow, directions 1 and 2 are
    !> still used for swapping within the periodic sub-volume.
    integer function get_recv_partner(curdim,dir)
        integer,intent(in) :: curdim,dir
        integer fnc,n
        real(kind=rk) rlpz

        if (boundary/='periodic_inflow'.or.curdim<3) then
           if (dir>2.or.cdims(curdim)==1) then
              ! Don't send with dir>2 in this case; also don't
              ! exchange if only 1 box in a dimension, MD Lees-Edwards
              ! directions are set up anew in each (sub-)step
              get_recv_partner = MPI_PROC_NULL
           else if (md_leesedwards.and.curdim==1.and.&
                &((ccoords(1)==0.and.dir==DIR_UP)&
                &.or.(ccoords(1)==cdims(1)-1.and.dir==DIR_DOWN))) then
              ! no conventional DIR_DOWN/DIR_UP recv across
              ! Lees-Edwards plane
              get_recv_partner = MPI_PROC_NULL
           else
              select case (dir)
              case (DIR_DOWN)
                 get_recv_partner = nnprocs(curdim,2)
              case (DIR_UP)
                 get_recv_partner = nnprocs(curdim,1)
              end select
           end if

!!$        else if (dir<3) then
!!$           rlpz = real(last_periodic_z,kind=rk)
!!$           ! dir is one of {1,2} --- z-swap with BC periodic_inflow
!!$           ! within periodic sub-volume
!!$
!!$           if (start(3)>last_periodic_z) then
!!$              get_recv_partner = MPI_PROC_NULL
!!$           else if (last_periodic_z<=nz) then
!!$              ! don't exchange if only 1 box in a dimension
!!$              get_recv_partner = MPI_PROC_NULL
!!$           else
!!$              if (dir==1) then
!!$                 if (border(2,3)>rlpz) then
!!$                    get_recv_partner = owning_rankc((/start(1),start(2),1/))
!!$                 else
!!$                    get_recv_partner = nnprocs(3,2)
!!$                 end if
!!$              else              ! dir==2 (positive z)
!!$                 if (ccoords(3)==0) then
!!$                    get_recv_partner = &
!!$                         &owning_rankc((/start(1),start(2),last_periodic_z/))
!!$                 else
!!$                    get_recv_partner = nnprocs(3,1)
!!$                 end if
!!$              end if
!!$           end if
!!$        else
!!$           ! dir is one of {3,4} --- z-swap with BC periodic_inflow
!!$           ! within non-periodic sub-volume
!!$           if (border(2,3)<enslavement_threshold) then
!!$              get_recv_partner = MPI_PROC_NULL
!!$           else
!!$              if (dir==3) then  ! negative z
!!$                 if (border(2,3)>real(tnz,kind=rk)) then
!!$                    get_recv_partner = MPI_PROC_NULL
!!$                 else
!!$                    get_recv_partner = nnprocs(3,2)
!!$                 end if
!!$              else              ! positive z (dir==4)
!!$                 if (border(1,3)<enslavement_threshold) then
!!$                    get_recv_partner = MPI_PROC_NULL
!!$                 else
!!$                    get_recv_partner = nnprocs(3,1)
!!$                 end if
!!$              end if
!!$           end if
        end if
    end function get_recv_partner

    !> calculate swap boundaries for all swaps in a given direction
    !>
    !> \param[in,out] sd swap direction
    subroutine calculate_send_bounds(sd)
        type(swap_direction_type),intent(inout) :: sd
        integer n
!!$        integer cdims_pi_z

        if (sd%scnt==0) return

        if (boundary/='periodic_inflow'.or.sd%dim<3) then
           select case (sd%dir)
           case (DIR_DOWN,DIR_LE_LO_Z_DOWN,DIR_LE_HI_Z_DOWN)
              ! simple periodic downward, in the case of Lees-Edwards,
              ! x send boundaries are just the same

              do n=1,sd%scnt-1
                 sd%s(n)%sbndlo = border(1,sd%dim)+(n-1)*chunksize(sd%dim)
                 call md_clip_periodic(sd%dim,sd%s(n)%sbndlo)
                 sd%s(n)%sbndhi = sd%s(n)%sbndlo+chunksize(sd%dim)
              end do

              ! last send
              sd%s(sd%scnt)%sbndlo &
                   &= border(1,sd%dim)+(sd%scnt-1)*chunksize(sd%dim)
              call md_clip_periodic(sd%dim,sd%s(sd%scnt)%sbndlo)
              sd%s(sd%scnt)%sbndhi = border(1,sd%dim)+rs
              call md_clip_periodic(sd%dim,sd%s(sd%scnt)%sbndhi)

              ! never send particles around both directions. This can
              ! apply to odd and even cdims(sd%dim) when rs close to
              ! tsize(sd%dim)/2, but only to the last send.
              if (sd%scnt==cdims(sd%dim)/2) sd%s(sd%scnt)%sbndhi = &
                   &min(sd%s(sd%scnt)%sbndhi&
                   &,sd%s(sd%scnt)%sbndlo+chunksize(sd%dim)/2.0_rk)

           case (DIR_UP,DIR_LE_LO_Z_UP,DIR_LE_HI_Z_UP)
              ! simple periodic upward, in the case of Lees-Edwards, x
              ! send boundaries are just the same

              do n=1,sd%scnt-1
                 sd%s(n)%sbndhi = border(2,sd%dim)-(n-1)*chunksize(sd%dim)
                 if (sd%s(n)%sbndhi<=minpos(sd%dim)) &
                      &sd%s(n)%sbndhi = sd%s(n)%sbndhi+tsize(sd%dim)
                 sd%s(n)%sbndlo = sd%s(n)%sbndhi-chunksize(sd%dim)
              end do

              ! last send
              sd%s(sd%scnt)%sbndhi &
                   &= border(2,sd%dim)-(sd%scnt-1)*chunksize(sd%dim)
              ! md_clip_periodic() would not clip at equality here...
              if (sd%s(sd%scnt)%sbndhi<=minpos(sd%dim)) &
                   &sd%s(sd%scnt)%sbndhi &
                   &= sd%s(sd%scnt)%sbndhi+tsize(sd%dim)
              sd%s(sd%scnt)%sbndlo = border(2,sd%dim)-rs
              call md_clip_periodic(sd%dim,sd%s(sd%scnt)%sbndlo)

              ! never send particles around both directions. This can
              ! apply to odd and even cdims(sd%dim) when rs close to
              ! tsize(sd%dim)/2, but only to the last send.
              if (sd%scnt==cdims(sd%dim)/2) sd%s(sd%scnt)%sbndlo = &
                   &max(sd%s(sd%scnt)%sbndlo&
                   &,sd%s(sd%scnt)%sbndhi-chunksize(sd%dim)/2.0_rk)
           end select

!!$        else                  ! now we have sd%dim==3 and periodic_inflow
!!$           select case(sd%dir)
!!$           case (1)             ! periodic downward
!!$              cdims_pi_z = ceiling(real(last_periodic_z,kind=rk)/nz)
!!$
!!$              do n=1,sd%scnt-1
!!$                 sd%s(n)%sbndlo = border(1,sd%dim)+n_boxes_pi_upward_z(n-1)
!!$                 call md_clip_periodic(3,sd%s(n)%sbndlo)
!!$                 ! don't clip any but the last send. If already an
!!$                 ! earlier send would require clipping, the next
!!$                 ! interval would definitely be negative. Trust
!!$                 ! get_sends_required() that such sends do not exist.
!!$                 sd%s(n)%sbndhi = border(1,3)+real(n_boxes_pi_upward_z(n-1)&
!!$                      &+calculate_chunksize_pi_z(ccoords(3)+n-1),kind=rk)
!!$                 if (sd%s(n)%sbndhi>maxpos_pi(3)) &
!!$                      &sd%s(n)%sbndhi = sd%s(n)%sbndhi-tsize_pi(3)
!!$              end do
!!$
!!$              ! last send
!!$              sd%s(sd%scnt)%sbndlo = border(1,3)+n_boxes_pi_upward_z(sd%scnt-1)
!!$              call md_clip_periodic(sd%dim,sd%s(sd%scnt)%sbndlo)
!!$              ! never send particles around both directions. Limit
!!$              ! boundary to 1/2 of the periodic range (without the
!!$              ! current target chunk length)
!!$              sd%s(n)%sbndhi = border(1,3)+min(rs&
!!$                   &,+0.5_rk*(last_periodic_z&
!!$                   &-calculate_chunksize_pi_z(ccoords(3)-1)))
!!$              call md_clip_periodic(sd%dim,sd%s(sd%scnt)%sbndhi)
!!$
!!$           case (2)             ! periodic upward
!!$              cdims_pi_z = ceiling(real(last_periodic_z,kind=rk)/nz)
!!$
!!$              do n=1,sd%scnt-1
!!$                 sd%s(n)%sbndhi = border_pi(2,sd%dim)-&
!!$                      &n_boxes_pi_downward_z(n-1)
!!$                 if (sd%s(n)%sbndhi<=minpos(sd%dim)) &
!!$                      &sd%s(n)%sbndhi = sd%s(n)%sbndhi+tsize_pi(sd%dim)
!!$                 ! don't clip any but the last send. If already an
!!$                 ! earlier send would require clipping, the next
!!$                 ! interval would definitely be negative. Trust
!!$                 ! get_sends_required() that such sends do not exist.
!!$                 sd%s(n)%sbndlo = border_pi(2,3)&
!!$                      &-real(n_boxes_pi_downward_z(n),kind=rk)
!!$                 call md_clip_periodic(3,sd%s(n)%sbndlo)
!!$              end do
!!$
!!$              ! last send
!!$              sd%s(sd%scnt)%sbndhi = border_pi(2,sd%dim)-n_boxes_pi_downward_z(sd%scnt-1)
!!$              ! md_clip_periodic() would not clip at equality here...
!!$              if (sd%s(sd%scnt)%sbndhi<=minpos(sd%dim)) &
!!$                   &sd%s(sd%scnt)%sbndhi = sd%s(sd%scnt)%sbndhi+tsize_pi(sd%dim)
!!$              ! never send particles around both directions. Limit
!!$              ! boundary to 1/2 of the periodic range (without the
!!$              ! current target chunk length)
!!$              sd%s(sd%scnt)%sbndlo = border_pi(2,3)+max(-rs&
!!$                   &,-0.5_rk*(last_periodic_z&
!!$                   &-calculate_chunksize_pi_z(ccoords(3)+1)))
!!$              call md_clip_periodic(sd%dim,sd%s(sd%scnt)%sbndlo)
!!$
!!$           case (3)             ! non-periodic downward
!!$
!!$              ! This is simple: just do the same as for sd%dir==1 in a
!!$              ! solely periodic dimension, no need to worry about
!!$              ! periodic boundaries, if a boundary gets too high,
!!$              ! there will be no particles anyway. Correctly set sd%scnt
!!$              ! makes the interval [sbndlo:sbndhi] automatically larger than
!!$              ! one chunk for domains close to the lower non-periodic
!!$              ! boundary (if this is required), so everything can be
!!$              ! transferred in a single send provided that it is done
!!$              ! as part of the last swap.
!!$
!!$              do n=1,sd%scnt-1
!!$                 sd%s(n)%sbndlo = border(1,sd%dim)+(n-1)*chunksize(sd%dim)
!!$                 sd%s(n)%sbndhi = sd%s(n)%sbndlo+chunksize(sd%dim)
!!$              end do
!!$
!!$              ! last send
!!$              sd%s(sd%scnt)%sbndlo = border(1,sd%dim)+(sd%scnt-1)*chunksize(sd%dim)
!!$
!!$              ! clip to the upper non-periodic system boundary
!!$              sd%s(sd%scnt)%sbndhi = min(border(1,sd%dim)+rs,maxpos(3))
!!$
!!$           case (4)             ! non-periodic upward
!!$
!!$              ! Also here, don't worry about periodic boundaries. If
!!$              ! everything is set up correctly, minz will never be
!!$              ! reached since enslavement_threshold is sufficiently
!!$              ! close to last_periodic_z+0.5 with always
!!$              ! last_periodic_z>2rs and no communication beyond rs is
!!$              ! necessary.
!!$
!!$              do n=1,sd%scnt-1
!!$                 sd%s(n)%sbndhi = border(2,sd%dim)-(n-1)*chunksize(sd%dim)
!!$                 sd%s(n)%sbndlo = sd%s(n)%sbndhi-chunksize(sd%dim)
!!$              end do
!!$
!!$              ! last send
!!$              sd%s(sd%scnt)%sbndhi = border(2,sd%dim)-(sd%scnt-1)*chunksize(sd%dim)
!!$
!!$              ! clip to the minimum non-periodic position
!!$              sd%s(sd%scnt)%sbndlo = max(border(2,sd%dim)-rs,enslavement_threshold)
!!$           end select
        end if
    end subroutine calculate_send_bounds

    !> setup neighbor bins, spatial-decomposition communication patterns,
    !> error check
    subroutine setup_parallel
        real(kind=rk),parameter :: small=1.0E-12_rk
        integer i,ierror,k,ix,iy,iz,stat,status(MPI_STATUS_SIZE)

        call setup_mpitypes

        do k=1,3                ! dimension
           do i=1,n_dir_max     ! (pseudo-)direction
              sdirs(k,i)%sproc = &
                   &get_send_partner(k,i,leesedwards=md_leesedwards)
              sdirs(k,i)%scnt = get_sends_required(k,i)
              if (sdirs(k,i)%scnt>nsmax) then
                 write (unit=6,fmt=&
                      &'("dim=",I1,", dir=",I1,", scnt=",I6,", nsmax=",I6)') &
                      &k,i,sdirs(k,i)%scnt,nsmax
                 call error_md('nsmax<scnt - boost nsmax!')
              end if

              call calculate_send_bounds(sdirs(k,i))

              sdirs(k,i)%rproc = get_recv_partner(k,i)
           end do
        end do

        ! preliminary MD Lees-Edwards setup so recvs_required can be
        ! sent around
        if (md_leesedwards) &
             &call md_leesedwards_provisional_communication_setup()

        ! obtain recvs_required just from the respective sender, this
        ! also is a test that things were set up not completely
        ! incosistently
        do k=1,3                ! dimension
           do i=1,n_dir_max             ! all non-Lees-Edwards (pseudo-)directions
              call MPI_Sendrecv&
                   &(sdirs(k,i)%scnt,1,MPI_INTEGER,sdirs(k,i)%sproc,i&
                   &,sdirs(k,i)%rcnt,1,MPI_INTEGER,sdirs(k,i)%rproc,i&
                   &,comm_cart,status,ierror)
           end do
        end do

#ifdef DEBUG_REPORTMDCOMM
        call debug_reportmdcomm()
        if (md_leesedwards) call md_leesedwards_debug_reportmdcomm()
#endif
        call debug_check_send_bounds ! this cannot cost much but may help a lot!

        ! setup neighbor binning parameters in box owned by each processor
        !  addding small -> bins slightly larger
        !  prevents round-off error (bin # ix = nbinx) when atoms are binned
        if (ineigh.eq.1) then
           nbinx = tnx/rs
           nbiny = tny/rs
           nbinz = tnz/rs
           binsizex = (tnx+small*tnx) / nbinx
           binsizey = (tny+small*tny) / nbiny
           binsizez = (tnz+small*tnz) / nbinz
           ! now binsize is basically rs, but a bit larger

           ! coordinates of bin just below the own processing domain
           mbinxlo = int((border(1,1) - minx) / binsizex) - 1
           mbinylo = int((border(1,2) - miny) / binsizey) - 1
           mbinzlo = int((border(1,3) - minz) / binsizez) - 1

           ! coordinates of bin just above the own processing domain...
           ix = int((border(2,1) - minx) / binsizex) + 1
           iy = int((border(2,2) - miny) / binsizey) + 1
           iz = int((border(2,3) - minz) / binsizez) + 1

           ! ...but not further than the last bin
           if (ix.gt.nbinx) ix = nbinx
           if (iy.gt.nbiny) iy = nbiny
           if (iz.gt.nbinz) iz = nbinz

           ! number of bins on the local processor plus halo of 1 on both ends
           mbinx = ix - mbinxlo + 1
           mbiny = iy - mbinylo + 1
           mbinz = iz - mbinzlo + 1

           ! clip to total number of bins in case of only one cpu per direction
           mbinx = min(mbinx,nbinx)
           mbiny = min(mbiny,nbiny)
           mbinz = min(mbinz,nbinz)

           ! clip  mbin[xyz]lo  to non-negative values
           if (mbinxlo.lt.0) mbinxlo = nbinx - 1
           if (mbinylo.lt.0) mbinylo = nbiny - 1
           if (mbinzlo.lt.0) mbinzlo = nbinz - 1

           ! RESULT
           ! - mbin[xyz]
           !   number of bins that are of interest in each direction on the
           !   local process
           ! - mbin[xyz]lo
           !   coordinates of first bin of interest for the local process,
           !   maybe wrapped around periodic boundary conditions
        endif

        if (any(rs>=tsize(:)/2.0_rk))&
             & call error_md('Outer cutoff rs >= 1/2 system size---'&
             &//'reduce rs or increase min((/nx,ny,nz/))!')
        if (ineigh.eq.1.and.(nbinx.le.2.or.nbiny.le.2.or.nbinz.le.2))&
             &call error_md('Two or less bins in a dimension---'&
             &//'reduce rs, set ineigh==0, or increase min((/nx,ny,nz/))!')

        if (ineigh==1) then
           allocate(binpnt(mbinx*mbiny*mbinz),stat=stat)
           call check_allocate(stat,'setup_parallel(): binpnt')
        end if
    end subroutine setup_parallel

    !> Create mpi types that represent a single particle in different
    !> situations where communication of particles is necessary.
    subroutine setup_mpitypes
        call build_particle_mpitype(particle_comm_mpitype,x=.true.&
             &,v=communicate_velocities,q=use_rotation,w=communicate_rotations&
             &,vws_fluid_avg=use_ft_fluid,uid=communicate_uids&
             &,ws=communicate_rotations_s,master=.true.,R=polydispersity,mag=magdispersity&
#ifdef RWALK
             &,v_r=communicate_velocities_r,sdx=communicate_dsx&
#endif
#ifdef LADD_SURR_RHOF
             &,rhof=.true.&
#endif
             &)

        call build_particle_mpitype(particle_exch_mpitype,x=.true.,v=.true.&
             &,q=use_rotation,w=use_rotation,ft_fluid=use_ft_fluid&
             &,ft_fluid_prev=decouple_fluid,vws_fluid=use_ft_fluid&
             &,vws_fluid_acc=decouple_fluid,vws_fluid_avg=decouple_fluid&
             &,uid=.true.,vqw_new=.true.,master=.true.,R=polydispersity, mag=magdispersity&
#ifdef RWALK
             &,v_r=.true.,sdx=.true.&
#endif
#ifdef LADD_SURR_RHOF
             &,rhof=.true.&
#endif
             &)

        call build_particle_mpitype(particle_cp_mpitype,x=.true.,v=.true.&
             &,q=use_rotation,w=use_rotation,ft_fluid=use_ft_fluid&
             &,vws_fluid=use_ft_fluid,vws_fluid_acc=decouple_fluid,uid=.true.&
             &,vqw_new=.true.,ws=steps_per_lbe_step==1.and.use_rotation&
             &,master=.true.,R=polydispersity, mag=magdispersity&
#ifdef RWALK
             &,v_r=.true.,sdx=.true.&
#endif
#ifdef LADD_SURR_RHOF
             &,rhof=.true.&
#endif
             &)

        if (collect_forces) then
           call build_collect_mpitypes(particle_coll_mpitype&
                &,ft_buffer_coll_mpitype,f=.true.,t=use_rotation&
                &,ft_fluid=use_ft_fluid&
#ifdef FORCECOMPONENT
                &,f_normal=dump_force_components&
                &,f_tangent=dump_force_components&
                &,f_n=dump_force_components&
                &,f_t=dump_force_components&
#endif
#ifdef PARTICLESTRESS
                &,tau=.true.&
#endif
#ifdef LADD_SURR_RHOF
                &,rhofn_acc=.true.&
#endif
                &)

           call build_collect_mpitypes(particle_coll_dump_mpitype&
                &,ft_buffer_coll_dump_mpitype,f=.true.,t=use_rotation&
                &,ft_fluid=use_ft_fluid,e_pot=dump_potentials&
#ifdef FORCECOMPONENT
                &,f_normal=dump_force_components&
                &,f_tangent=dump_force_components&
                &,f_n=dump_force_components&
                &,f_t=dump_force_components&
#endif
#ifdef PARTICLESTRESS
                &,tau=.true.&
#endif
#ifdef LADD_SURR_RHOF
                &,rhofn_acc=.true.&
#endif
                &)

           call build_collect_mpitypes(particle_coll_substep_mpitype&
                &,ft_buffer_coll_substep_mpitype,f=.true.,t=use_rotation&
                &,e_pot=dump_potentials&
#ifdef FORCECOMPONENT
                &,f_normal=dump_force_components&
                &,f_tangent=dump_force_components&
                &,f_n=dump_force_components&
                &,f_t=dump_force_components&
#endif
#ifdef PARTICLESTRESS
                &,tau=.true.&
#endif
                &)
        end if
    end subroutine setup_mpitypes

    !> initialize particles according to initial_placing
    !>
    !> \todo Several things...
    !> \li Finish parallel random initialization and add \c 'ndom'
    !>     to \c 'ra' in the \c if statement.
    !> \li Parallelize all initialization routines.
    !> \li One could check for overlap with walls in the future
    !>     and maybe delete misplaced particles again.
    subroutine setup_particles(N)
      type(lbe_site),intent(in) :: N(1-halo_extent:,1-halo_extent:&
           &,1-halo_extent:)
      real(kind=rk),allocatable,dimension(:,:,:) :: rock_state
      type(md_particle_type),allocatable,dimension(:) :: tp
      integer counts(0:nprocs-1),first_uids(0:nprocs-1),first_uid
      integer :: i,ii,ierror,n_global,stat
      logical :: assigned_uids

      ! Parallel initialization routines should be used in all cases
      ! in the future, at the moment they are able to handle only
      ! special cases:

      ! Changed the nice and compact array statement to an ugly set of
      ! 'or' clauses because FORTRAN doesn't always like different
      ! string lengths in array initialization (WHY???)...
      parallel: if (initial_placing=='file_parallel:asc/xvow'&
           &.or.initial_placing=='file_parallel:asc/xvoW'&
           &.or.initial_placing=='file_parallel:asc/xvowr'&
           &.or.initial_placing=='file_parallel:asc/xvowrm'&
           &.or.initial_placing=='file_parallel:asc/xvowm'&
           &.or.initial_placing=='sc'.or.initial_placing=='random'&
           &.or.initial_placing=='sphere') then

        ! This is parallel particle placement.
        call setup_list()

        ! this is the default but some initial_placings could assign
        ! uids themselves
        assigned_uids = .false.

        select case (initial_placing)
          case ('file_parallel:asc/xvow')
            call log_msg_md(&
                 &"Initializing particles using <file_parallel:asc/xvow>...")
            call place_particles_file_parallel_asc_xvow()
            assigned_uids = .true.
          case ('file_parallel:asc/xvoW')
            call log_msg_md(&
                 &"Initializing particles using <file_parallel:asc/xvoW>...")
            call place_particles_file_parallel_asc_xvo_ws()
            assigned_uids = .true.
          case ('file_parallel:asc/xvowr')
            call log_msg_md(&
                 &"Initializing particles using <file_parallel:asc/xvowr>...")
            call place_particles_file_parallel_asc_xvowr()
           assigned_uids = .true.
          case ('file_parallel:asc/xvowrm')
             call log_msg_md(&
                   &"Initializing particles using <file_parallel:asc/xvowrm>...")
             call place_particles_file_parallel_asc_xvowrm()
            assigned_uids = .true.
           case ('file_parallel:asc/xvowm')
             call log_msg_md(&
                   &"Initializing particles using <file_parallel:asc/xvowrm>...")
             call place_particles_file_parallel_asc_xvowm()
            assigned_uids = .true.
          case ('random')
            call log_msg_md(&
                 &"Initializing particles using <random> (parallel)...")
            call place_particles_random(N)
          case ('sc')
            call log_msg_md("Initializing particles using <sc> (parallel)...")
            call place_particles_sc(N)
          case ('sphere')
            call log_msg_md("Initializing particles using <sphere>...")
            call place_particles_sphere(N)
          case default
            call error_md('Unknown value: initial_placing = <'//initial_placing&
                 &//'>')
        end select

        select case (initial_orientations)
          case ('none')
            ! NOOP
            call log_msg_md("Skipped explicit initialization of orientations.")
          case ('q0')
            call log_msg_md("Initializing orientations from q0.")
            call init_orientations_q0()
          case ('random')
            call log_msg_md("Initializing random orientations.")
            call init_orientations_random()
          case ('sphere')
            call log_msg_md("Initializing orientations normal to sphere.")
            call init_orientations_sphere()
          case default
            call error_md('Unknown value: initial_orientations = <'&
                 &//trim(initial_orientations)//'>')
        end select

        call MPI_Gather(nlocal,1,MPI_INTEGER,counts,1,MPI_INTEGER,0&
             &,comm_cart,ierror)

        ! assign unique ids, if not done already
        if (.not.assigned_uids) then
           call calculate_displacements(counts,first_uids)
           call MPI_Scatter(first_uids,1,MPI_INTEGER,first_uid,1,MPI_INTEGER,0&
                &,comm_cart,ierror)
           i = atompnt
           do ii=1,nlocal
              P(i)%uid = first_uid+ii
              i = list(i)
           end do
        end if

        write(msgstr,"('Placed a total of ',I0,' particles (parallel).')")&
             & sum(counts)
        call log_msg_md(msgstr)
      else parallel
        ! These lines are for serial placement on the root process.

        ! That much memory is only required if there is actual
        ! rock. Still it would be better not to rely on global
        ! arrays also if there is rock.
        if (rock_is_present()) then
          if (myrankc==0) then
            allocate (rock_state&
                 &(1-halo_extent:tnx+halo_extent&
                 &,1-halo_extent:tny+halo_extent&
                 &,1-halo_extent:tnz+halo_extent)&
                 &,stat=stat)
            call check_allocate(stat,'setup_particles(): rock_state')
          end if
          call gather_rock_state(N,rock_state)
        end if

        rank0: if (myrankc==0) then
          n_global = 0
          select case (initial_placing)
            case ('facex','facey','facez')
              call log_msg_md("Initializing particles using <faceX>...")
              call place_particles_face(rock_state,tp,n_global)
            case ('fcc')
              call log_msg_md("Initializing particles using <fcc>...")
              call place_particles_fcc(rock_state,tp,n_global)
            case ('file:asc/xvow')
              call log_msg_md(&
                   &"Initializing particles using <file:asc/xvow>...")
              call place_particles_file_asc_xvow(rock_state,tp,n_global)
            case ('file:asc/xvqw')
              call log_msg_md(&
                   &"Initializing particles using <file:asc/xvqw>...")
              call place_particles_file_asc_xvqw(rock_state,tp,n_global)
             case ('file:asc/xvowrm')
                call log_msg_md(&
                     &"Initializing particles using <file:asc/xvowrm>...")
              call place_particles_file_asc_xvowrm(rock_state,tp,n_global)
            case ('lyapunov')
              call log_msg_md("Initializing particles using <lyapunov>...")
              call place_particles_lyapunov(tp,n_global)
            case ('random_serial')
              call log_msg_md("Initializing particles using <random_serial>...")
              call place_particles_random_serial(rock_state,tp,n_global)
            case ('sc_serial')
              call log_msg_md(&
                   &"Initializing particles using <sc_serial>...")
              call place_particles_sc_serial(rock_state,tp,n_global)
            case default
              call error_md('Unknown value: initial_placing = <'&
                   &//trim(initial_placing)//'>')
            end select
            write(msgstr,"('Placed a total of ',I0,' particles (serial).')")&
                 & n_global
            call log_msg_md(msgstr)
        end if rank0

        call setup_list()
        call scatter_particles(n_global,tp)
        if (myrankc==0) deallocate (tp)
        if (rock_is_present().and.myrankc==0) deallocate (rock_state)
      end if parallel

      select case (initial_velocities)
        case ('none')
          ! NOOP
          call log_msg_md("Skipped explicit initialization of velocities.")
        case ('v0')
          call log_msg_md("Initializing velocities from v0.")
          call init_velocities_v0()
        case ('temp0')
          call log_msg_md("Initializing velocities from temp0.")
          call init_velocities_temp0()
#ifdef RWALK
          call log_msg_md("Initializing random velocities.")
          call init_velocities_tempr()
#endif
        case default
          call error_md('Unknown value: initial_velocities = <'&
               &//trim(initial_velocities)//'>')
      end select

      select case (initial_rotations)
        case ('none')
          ! NOOP
          call log_msg_md("Skipped explicit initialization of rotations.")
        case ('w0')
          call log_msg_md("Initializing rotations from w0.")
          call init_rotations_w0()
        case default
          call error_md('Unknown value: initial_rotations = <'&
                &//trim(initial_rotations)//'>')
      end select

      select case (initial_radii)
      case ('none')
         ! NOOP
      case ('normal')
         call log_msg_md("Initializing normally distributed random particle radii.")
         call init_radii_normal()
      case ('random')
         call log_msg_md("Initializing random particle radii.")
         call init_radii_random()
      case ('volume')
         call log_msg_md("Initializing random particle radii leading to constant volume.")
         call init_radii_volume()
      case default
         call error_md('unknown value initial_radii=<'//trim(initial_radii)&
              &//'>---choose either <none>, <normal>, <random>, or <volume>')
      end select

      i = atompnt
      do ii = 1,nlocal
        ! initialize forces and torques to zero (not required for
        ! results, just avoid confusing output of
        ! dump_configuration() called in md_init() which is before
        ! integrate() )
        P(i)%f(:) = 0.0_rk
        P(i)%t(:) = 0.0_rk
#ifdef FORCECOMPONENT
        P(i)%f_normal(:)  = 0.0_rk
        P(i)%f_tangent(:) = 0.0_rk
        P(i)%f_n(:)  = 0.0_rk
        P(i)%f_t(:) = 0.0_rk
#endif
        i = list(i)
      end do

      ! needed to initialize ws_fluid and ws_fluid_acc below
      if (decouple_fluid.or.use_ft_fluid) call calculate_rotations_s()

      ! This is not exactly correct but should not matter for the first
      ! timestep.
      i = atompnt
      do ii = 1,nlocal
        P(i)%vnew = P(i)%v
        P(i)%wnew = P(i)%w
        P(i)%qnew = P(i)%q

        if (use_ft_fluid) then
          P(i)%v_fluid = P(i)%v
          P(i)%ws_fluid = P(i)%ws

          P(i)%v_fluid_avg = P(i)%v_fluid
          P(i)%ws_fluid_avg = P(i)%ws_fluid
        end if

        if (decouple_fluid) then
          ! this results in v_fluid==v and ws_fluid==ws during the
          ! first LB step
          P(i)%v_fluid_acc = real(steps_per_lbe_step-1,kind=rk)*P(i)%v
          P(i)%ws_fluid_acc = real(steps_per_lbe_step-1,kind=rk)*P(i)%ws
        end if

        ! necessary since [ft]_fluid_prev will be initialized from
        ! [ft]_fluid later on
        P(i)%f_fluid = 0.0_rk
        P(i)%t_fluid = 0.0_rk
        if (calculate_orientations) call calculate_all_orientations
        i = list(i)
      end do
      call log_msg_md("Finished particle setup.")
    end subroutine setup_particles

    !> allocate and initialize \c uid2i
    subroutine setup_uid2i
        call Mii_init(uid2i)
    end subroutine setup_uid2i

    !> place a single particle on the local process for initialization
    !>
    !> \param[in] pos particle position
    !>
    !> \param[in] v particle velocity (optional, defaults to undefined)
    !>
    !> \param[in] q particle orientation (quaternion; optional,
    !> defaults to undefined))
    !>
    !> \param[in] w particle body-fixed angular velocity (optional,
    !> defaults to undefined)
    !>
    !> \param[in] uid particle uid (optional, defaults to undefined)
    subroutine place_particle(pos,v,q,w,R_orth,R_para,mag,uid)
        real(kind=rk),intent(in) :: pos(3)
        real(kind=rk),intent(in),optional :: v(3),q(0:3),w(3),R_orth,R_para,mag(2)
        integer,intent(in),optional :: uid
        integer itmp,k

        if (nlocal+1>npmax) call boost_npmax(nlocal+1)
        nlocal = nlocal + 1

        itmp = atompnt
        atompnt = freepnt
        freepnt = list(freepnt)
        list(atompnt) = itmp

        P(atompnt)%x = pos
        if (present(v)) P(atompnt)%v = v
        if (present(q)) P(atompnt)%q = q
        if (present(w)) P(atompnt)%w = w
        if (present(R_orth)) P(atompnt)%R_orth = R_orth
        if (present(R_para)) P(atompnt)%R_para = R_para
        if (present(mag)) P(atompnt)%mag = mag
        if (present(uid)) P(atompnt)%uid = uid
    end subroutine place_particle




    !> place a single particle, but check whether \c pos  is valid and do not
    !> set particles straight into rock. The subroutine (re)allocates memory
    !> for \c tp if needed, assuming \c n_global to be the number of particles
    !> that is in \c tp already.
    subroutine place_particle_serial(rock_state,pos,tp,n_global,success)
        real(kind=rk),intent(in) :: rock_state&
             &(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        real(kind=rk),intent(in) :: pos(3)
        type(md_particle_type),allocatable,intent(inout) :: tp(:)
        integer,intent(inout) :: n_global
        logical,intent(out),optional :: success
        type(md_particle_type),allocatable :: tmp(:)
        real(kind=rk) :: ori(0:3)
        integer stat,k
        logical problem_with_rock

        if (present(success)) success = .false.

        if (any(pos(:)<minpos(:)).or.any(pos(:)>=maxpos(:))) then
           print *,'pos(:)=',pos
           call error_md(&
                &'tried to place particle outside of the global sim space')
        end if

        select case (initial_orientations)
        case ('none')
           ! nop
        case ('q0')
           ori = q0
        case ('random')
           do k=0,3
              call uni_dist(ori(k))
           end do
           ori = ori/sqrt(dot_product(ori,ori))
        case default
           call error_md('unknown value: initial_orientations="'&
                &//initial_orientations//'"')
        end select

        ! only check for too close rock if  rock_state  is available
        if (rock_is_present()) then
           problem_with_rock = too_close_to_rock(rock_state,pos,ori)
        else
           problem_with_rock = .false.
        end if

        if (.not.problem_with_rock) then
           if (.not.allocated(tp)) then
              allocate (tp(1),stat=stat)
              call check_allocate(stat,'place_particle_serial(): tp')
           end if
           if (size(tp)==n_global) then
              allocate (tmp(n_global),stat=stat)
              call check_allocate(stat,'place_particle_serial(): tmp')
              tmp(:) = tp(:)
              deallocate (tp)
              allocate (tp(2*n_global),stat=stat)
              call check_allocate(stat,'place_particle_serial(): tp')
              tp(1:n_global) = tmp(:)
              deallocate(tmp)
           end if
           if (present(success)) success = .true.
           n_global = n_global+1
           tp(n_global)%x(:) = pos(:)
           tp(n_global)%q(:) = ori(:)
           tp(n_global)%uid = n_global
        endif
    end subroutine place_particle_serial

    !> places particles on a single layer
    subroutine place_particles_face(rock_state,tp,n_global)
        real(kind=rk),intent(in) :: rock_state&
             &(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        type(md_particle_type),allocatable,intent(inout) :: tp(:)
        integer,intent(inout) :: n_global
        integer ndir,fdirs(2) ! face directions and direction normal to them
        integer ncells(2)       ! number of "cells" in face directions
        integer i,j
        real(kind=rk) :: tmp(3),xlo(3),xhi(3)
        real(kind=rk),parameter :: small=1.0E-12_rk ! prevents roundoff errors

        xlo(:) = max(x0lo(:),minpos(:))
        call find_highest_lattice_position(xhi)

        dirs: select case (initial_placing(5:5))
        case ('x')
           ndir = 1
           fdirs = (/ 2,3 /)
        case ('y')
           ndir = 2
           fdirs = (/ 1,3 /)
        case ('z')
           ndir = 3
           fdirs = (/ 1,2 /)
        end select dirs

        tmp(ndir) = xlo(ndir)

        ncells(:) = (small+xhi(fdirs(:))-xlo(fdirs(:)))/alat
        do i=0,ncells(1)-1
           do j=0,ncells(2)-1
              tmp(fdirs(:)) = (/i,j/)*alat + xlo(fdirs(:))
              call place_particle_serial(rock_state,tmp(:),tp,n_global)
           end do
        end do
    end subroutine place_particles_face

    !> fills the global domain with particles on an fcc lattice
    subroutine place_particles_fcc(rock_state,tp,n_global)
        real(kind=rk),intent(in) :: rock_state&
             &(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        type(md_particle_type),allocatable,intent(inout) :: tp(:)
        integer,intent(inout) :: n_global
        integer i,j,k,ncells(3)
        real(kind=rk) :: xlo(3),xhi(3)
        real(kind=rk),parameter :: small=1.0E-12_rk ! prevents roundoff errors

        xlo(:) = max(x0lo(:),minpos(:))
        call find_highest_lattice_position(xhi)

        ncells(:) = (small+xhi(:)-xlo(:))/alat

        do k = 1,ncells(3)*2
           do j = 1,ncells(2)*2
              do i = 1,ncells(1)*2
                 if (mod(i+j+k,2).eq.1) then
                    call place_particle_serial(rock_state&
                         &,xlo(:)+(/i-1,j-1,k-1/)*alat/2.0_rk,tp,n_global)
                 endif
              enddo
           enddo
        enddo
    end subroutine place_particles_fcc

    !> Initialize particles ( x ,  v ,  q , and  w ) from the file specified
    !> by  init_file .  uid  is set according to the order of the particles
    !> read.
    subroutine place_particles_file_asc_xvow(rock_state,tp,n_global)
        real(kind=rk),intent(in) :: rock_state&
             &(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        type(md_particle_type),allocatable,intent(inout) :: tp(:)
        integer,intent(inout) :: n_global
        integer,parameter :: init_file_unit=12
        real(kind=rk) :: x(3),v(3),o(3),w(3)
        logical :: success

        write(msgstr,"('Reading initial particle configuration (x,v,o,w) from file <',A,'>...')") trim(init_file)
        call log_msg_md(msgstr)

        open (unit=init_file_unit,file=init_file,status='OLD',action='READ'&
             &,err=200)

        do
           read (unit=init_file_unit,fmt=*,end=100) x(:),v(:),o(:),w(:)
           call place_particle_serial(rock_state,x,tp,n_global,success)
           if (success) then
              tp(n_global)%v = v
              tp(n_global)%q = quaternion(o)
              tp(n_global)%w = w
           end if
        end do
100     continue

        close (init_file_unit,err=200)

        return

200     continue
        call error_md('Error reading md initial configuration file "'&
             &//trim(init_file)//'"')
    end subroutine place_particles_file_asc_xvow

    !> Initialize particles ( x ,  v ,  q , and  w ) from the file specified
    !> by  init_file .  uid  is set according to the order of the particles
    !> read.
    subroutine place_particles_file_asc_xvqw(rock_state,tp,n_global)
        real(kind=rk),intent(in) :: rock_state&
             &(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        type(md_particle_type),allocatable,intent(inout) :: tp(:)
        integer,intent(inout) :: n_global
        integer,parameter :: init_file_unit=12
        real(kind=rk) :: x(3),v(3),q(0:3),w(3)
        logical :: success

        write(msgstr,"('Reading initial particle configuration (x,v,q,w) from file <',A,'>...')") trim(init_file)
        call log_msg_md(msgstr)

        open (unit=init_file_unit,file=init_file,status='OLD',action='READ'&
             &,err=200)

        do
           read (unit=init_file_unit,fmt=*,end=100) x(:),v(:),q(:),w(:)
           call place_particle_serial(rock_state,x,tp,n_global,success)
           if (success) then
              tp(n_global)%v(:) = v(:)
              tp(n_global)%q(:) = q(:)
              tp(n_global)%w(:) = w(:)
           end if
        end do
100     continue

        close (init_file_unit,err=200)

        return
200     continue
        call error_md('Error reading md initial configuration file "'&
             &//trim(init_file)//'"')
    end subroutine place_particles_file_asc_xvqw
   !> Initialize particles ( x ,  v ,  q , and  w ) from the file specified
      !> by  init_file .  uid  is set according to the order of the particles
      !> read.
      subroutine place_particles_file_asc_xvowrm(rock_state,tp,n_global)
          real(kind=rk),intent(in) :: rock_state&
               &(1-halo_extent:,1-halo_extent:,1-halo_extent:)
          type(md_particle_type),allocatable,intent(inout) :: tp(:)
          integer,intent(inout) :: n_global
          integer,parameter :: init_file_unit=12
          real(kind=rk) :: x(3),v(3),o(3),w(3),ro,rp,ma(2)
          logical :: success
  
          write(msgstr,"('Reading initial particle configuration (x,v,o,w,R,mag) from file <',A,'>...')") trim(init_file)
          call log_msg_md(msgstr)
  
          open (unit=init_file_unit,file=init_file,status='OLD',action='READ'&
               &,err=200)
  
          do
             read (unit=init_file_unit,fmt=*,end=100) x(:),v(:),o(:),w(:),ro,rp,ma(:)
             call place_particle_serial(rock_state,x,tp,n_global,success)
             if (success) then
                tp(n_global)%v = v
                tp(n_global)%q = quaternion(o)
                tp(n_global)%w = w
               tp(n_global)%R_orth = ro
               tp(n_global)%R_para = rp
					tp(n_global)%mag = ma
             end if
          end do
  100     continue
  
          close (init_file_unit,err=200)
  
          return

 200       continue
          call error_md('Error reading md initial configuration file "'&
               &//trim(init_file)//'"')
      end subroutine place_particles_file_asc_xvowrm


    !> Initialize particles \c x , \c v , \c q , and \c w from the
    !> file specified by \c init_file . \c uid is set according to the
    !> order of the particles read. Every process reads all particles
    !> and picks out the ones that sit in the local domain. No check
    !> for walls is performed.
    subroutine place_particles_file_parallel_asc_xvow()
        integer,parameter :: init_file_unit=12
        real(kind=rk) :: x(3),v(3),o(3),w(3)
        integer :: uid

        write(msgstr,"('All processes reading initial particle configuration "&
             &//"(x,v,o,w) from file <',A,'>...')") trim(init_file)
        call log_msg_md(msgstr)

        open (unit=init_file_unit,file=init_file,status='OLD',action='READ'&
             &,err=200)

        uid = 0
        do
           read (unit=init_file_unit,fmt=*,end=100) x,v,o,w
           uid = uid + 1        ! use record (particle) number as uid

           if (all(x>=border(1,:)).and.all(x<border(2,:))) then
              call place_particle(x,v=v,q=quaternion(o),w=w,uid=uid)
           end if
        end do
100     continue

        close (init_file_unit,err=200)

        return

200     continue
        call error_md('Error reading md initial configuration file "'&
             &//trim(init_file)//'"')
    end subroutine place_particles_file_parallel_asc_xvow

    !> Initialize particles \c x , \c v , \c q , and \c w from the
    !> file specified by \c init_file . \c uid is set according to the
    !> order of the particles read. Every process reads all particles
    !> and picks out the ones that sit in the local domain. No check
    !> for walls is performed.
    subroutine place_particles_file_parallel_asc_xvo_ws()
        integer,parameter :: init_file_unit=12
        real(kind=rk) :: x(3),v(3),o(3),ws(3)
        integer :: uid

        write(msgstr,"('All processes reading initial particle configuration "&
             &//"(x,v,o,ws) from file <',A,'>...')") trim(init_file)
        call log_msg_md(msgstr)

        open (unit=init_file_unit,file=init_file,status='OLD',action='READ'&
             &,err=200)

        uid = 0
        do
           read (unit=init_file_unit,fmt=*,end=100) x,v,o,ws
           uid = uid + 1        ! use record (particle) number as uid

           if (all(x>=border(1,:)).and.all(x<border(2,:))) then
              call place_particle(x,v=v,q=quaternion(o)&
                   &,w=matmul(space_to_body_matrix(quaternion(o)),ws),uid=uid)
           end if
        end do
100     continue

        close (init_file_unit,err=200)

        return

200     continue
        call error_md('Error reading md initial configuration file "'&
             &//trim(init_file)//'"')
    end subroutine place_particles_file_parallel_asc_xvo_ws

    !> Initialize particles \c x , \c v , \c q , \c w , \c R_orth ,
    !> and \c R_para from the file specified by \c init_file . \c uid
    !> is set according to the order of the particles read. Every
    !> process reads all particles and picks out the ones that sit in
    !> the local domain. No check for walls is performed.
    subroutine place_particles_file_parallel_asc_xvowr()
        integer,parameter :: init_file_unit=12
        real(kind=rk) :: x(3),v(3),o(3),w(3),ro,rp
        integer :: uid

        write(msgstr,"('All processes reading initial particle configuration "&
             &//"(x,v,o,w,R_orth,R_para) from file <',A,'>...')") &
             &trim(init_file)
        call log_msg_md(msgstr)

        open (unit=init_file_unit,file=init_file,status='OLD',action='READ'&
             &,err=200)

        uid = 0
        do
           read (unit=init_file_unit,fmt=*,end=100) x,v,o,w,ro,rp
           uid = uid + 1        ! use record (particle) number as uid

           if (max(ro,rp)>maximum_particle_radius()) then
              write (msgstr,fmt='("particle ",I0," has maximum radius "'&
                   &//',ES15.8," which is larger than the maximum expected '&
                   &//'radius ",ES15.8,"---define smaller particles or '&
                   &//'increase R_orth or R_para!")') &
                   &uid,max(ro,rp),maximum_particle_radius()
              call error_md(msgstr)
           end if

           if (all(x>=border(1,:)).and.all(x<border(2,:))) then
              call place_particle(x,v=v,q=quaternion(o),w=w,R_orth=ro,R_para=rp&
                   &,uid=uid)
           end if
        end do
100     continue

        close (init_file_unit,err=200)

        return

200     continue
        call error_md('Error reading md initial configuration file "'&
             &//trim(init_file)//'"')
    end subroutine place_particles_file_parallel_asc_xvowr

  !> Initialize particles \c x , \c v , \c q , \c w , \c R_orth ,
      !> ,\c R_para, \c mag from the file specified by \c init_file . \c uid
      !> is set according to the order of the particles read. Every
      !> process reads all particles and picks out the ones that sit in
      !> the local domain. No check for walls is performed.
      subroutine place_particles_file_parallel_asc_xvowrm()
          integer,parameter :: init_file_unit=12
          real(kind=rk) :: x(3),v(3),o(3),w(3),ro,rp, ma(2)
          integer :: uid
  
          write(msgstr,"('All processes reading initial particle configuration "&
               &//"(x,v,o,w,R_orth,R_para, mag) from file <',A,'>...')") &
               &trim(init_file)
          call log_msg_md(msgstr)
  
          open (unit=init_file_unit,file=init_file,status='OLD',action='READ'&
               &,err=200)
  
          uid = 0
          do
             read (unit=init_file_unit,fmt=*,end=100) x,v,o,w,ro,rp,ma
             uid = uid + 1        ! use record (particle) number as uid
  
             if (max(ro,rp)>maximum_particle_radius()) then
                write (msgstr,fmt='("particle ",I0," has maximum radius "'&
                     &//',ES15.8," which is larger than the maximum expected '&
                     &//'radius ",ES15.8,"---define smaller particles or '&
                     &//'increase R_orth or R_para!")') &
                     &uid,max(ro,rp),maximum_particle_radius()
                call error_md(msgstr)
             end if
  
             if (all(x>=border(1,:)).and.all(x<border(2,:))) then
                call place_particle(x,v=v,q=quaternion(o),w=w,R_orth=ro,R_para=rp&
                     &,mag=ma, uid=uid)
             end if
          end do
  100     continue
  
          close (init_file_unit,err=200)
  
          return
  
  200     continue
          call error_md('Error reading md initial configuration file "'&
               &//trim(init_file)//'"')
      end subroutine place_particles_file_parallel_asc_xvowrm
  !> Initialize particles \c x , \c v , \c q , \c w , \c mag from the file specified by \c init_file . \c uid
      !> is set according to the order of the particles read. Every
      !> process reads all particles and picks out the ones that sit in
      !> the local domain. No check for walls is performed.
      subroutine place_particles_file_parallel_asc_xvowm()
          integer,parameter :: init_file_unit=12
          real(kind=rk) :: x(3),v(3),o(3),w(3), ma(2)
          integer :: uid
  
          write(msgstr,"('All processes reading initial particle configuration "&
               &//"(x,v,o,w,mag) from file <',A,'>...')") &
               &trim(init_file)
          call log_msg_md(msgstr)
  
          open (unit=init_file_unit,file=init_file,status='OLD',action='READ'&
               &,err=200)
  
          uid = 0
          do
             read (unit=init_file_unit,fmt=*,end=100) x,v,o,w,ma
             uid = uid + 1        ! use record (particle) number as uid
  
      
  
             if (all(x>=border(1,:)).and.all(x<border(2,:))) then
                call place_particle(x,v=v,q=quaternion(o),w=w,mag=ma,uid=uid)
             end if
          end do
  100     continue
  
          close (init_file_unit,err=200)
  
          return
  
  200     continue
          call error_md('Error reading md initial configuration file "'&
               &//trim(init_file)//'"')
      end subroutine place_particles_file_parallel_asc_xvowm
    !> Initial condition for lbe_md_lyapunov_module
    !>
    !> Fills the system with particles on a lattice which consists of
    !> points similar to that of  'facez'  and 4 additional points arranged
    !> in a cross around each of the first points. The  uid  for the center
    !> particles is k*5+1 (k integer), the uids for the 4 particles around a
    !> center particle with uid==u are u+1,u+2,u+3,u+4. No check concerning
    !> rock_state  is performed, but there is some safety interval to the
    !> edges of the z=0 plane. The distance between two crosses is 3 lu,
    !> the distance between the center particles and the ones around them is
    !> 1 lu.
    !>
    !> \verbatim
    !> example (uid 'X' means 10, but '10' would break the layout):
    !>  2  7  .
    !> 513X68...
    !>  4  9  .
    !>  .  .  .
    !> .........
    !>  .  .  .
    !> \endverbatim
    !>
    !> Unlike the other particle placement routines, this one initializes not
    !> only the positions and the uids but also the velocities (to zero) of the
    !> particles.
    subroutine place_particles_lyapunov(tp,n_global)
        type(md_particle_type),allocatable,intent(inout) :: tp(:)
        integer,intent(inout) :: n_global
        integer i,j,k,p,stat
        integer sur(2,4) ! offset of surrounding points to their center point
        ! outermost lattice points that could become center points
        integer :: begin(2),end(2)
        integer :: extent(2)    ! # of centers in x and y direction

        begin(:) = (/3,3/)
        end(:) = (/tnx-2,tny-2/)
        extent(:) = (end(:)-begin(:))/3+1
        n_global = product(extent(:))*5

        sur(:,1) = (/0,-1/)
        sur(:,2) = (/-1,0/)
        sur(:,3) = (/0,1/)
        sur(:,4) = (/1,0/)

        allocate (tp(n_global),stat=stat)
        call check_allocate(stat,'place_particles_lyapunov(): tp')

        p = 1
        x: do i=begin(1),end(1),+3
           y: do j=begin(2),end(2),+3
              tp(p)%x(:) = real((/i,j,1/),kind=rk)
              tp(p)%v(:) = 0.0_rk
              tp(p)%uid = p
              sur_points: do k=1,4
                 tp(p+k)%x(:) = real((/i+sur(1,k),j+sur(2,k),1/),kind=rk)
                 tp(p+k)%v(:) = 0.0_rk
                 tp(p+k)%uid = p+k
              end do sur_points
              p = p + 5
           end do y
        end do x
    end subroutine place_particles_lyapunov

    !> Fills the local domain with particles at random positions in parallel.
    !>
    !> \param[in] N local lattice chunk (halo depth >=1 required)
    !>
    !> A minimum particle distance of \c min_dist is maintained. Note
    !> that though this is better than \c
    !> place_particles_random_serial(), it is still slow/does not work
    !> for high particle densities.
    subroutine place_particles_random(N)
        type(lbe_site),intent(in) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        integer i,ierror,j,k,pp,stat,x,y,z
        integer rc(3),rp(3),sp(3)
        integer lnp,np ! local and estimated total number of particles
        integer rlnp   ! remaining local number of particles
        integer send_more       ! indicator for sending remaining particles
        integer,allocatable :: npp(:) ! number of particles per process
        real(kind=rk) :: fnp,rnp ! actual and remaining number of particles
        real(kind=rk) :: r(1),sumrnp,V_chunk,V_fill

        if (min_dist>rs) call error_md(&
             &'Parallel random particle placement requires that rs>=min_dist'&
             &//'---increase rs or try initial_placing=''random_serial''.')

        if (myrankc==0) then
           ! start with the number of particles in the whole system
           ! (neglecting  x0lo  and  x0hi )
           if (phi<0.0_rk) then
              ! number of particles is determined as for an sc-lattice
              np = nint(product(tsize/alat))
              if (np<1) call error_md&
                   &('Volume concentration derived from alat would lead to '&
                   &//'less than 1 particle!')
           else
              ! number of particles is determined via volume concentration phi
              if (interaction/='ladd') call error_md&
                   &('Random particle placement based on volume concentration '&
                   &//'phi works only with interaction=''ladd''. Set phi<0 '&
                   &//'and use alat instead!')

              np = nint(phi*product(tsize)/avg_particle_volume_fluid_ladd())
              if (np<1) call error_md&
                   &('Desired volume concentration phi would lead to less '&
                   &//'than 1 particle!')
           end if

           ! distribute these particles to all processes, give each
           ! process the same number of particles and distribute the
           ! remainder randomly.
           allocate (npp(0:nprocs-1),stat=stat)
           call check_allocate(stat,'place_particles_random(): npp')
           npp = np/nprocs
           do i=(np/nprocs)*nprocs+1,np
              call ranlux(r,1)
              pp = floor(r(1)*nprocs)
              npp(pp) = npp(pp)+1
           end do
        end if

        call mpi_scatter(npp,1,MPI_INTEGER,lnp,1,MPI_INTEGER,0,comm_cart,ierror)

        ! reduce the local number of particles according to the actual
        ! volume to fill
        V_fill = fillable_volume(N)
        V_chunk = product(chunksize)
        fnp = real(lnp,kind=rk)*V_fill/V_chunk

        lnp = floor(fnp)

        ! redistribute remainder of fractional particles but not to
        ! nodes whose fillable volumes are less than 50% of their
        ! total volumes - that should reduce the probability of
        ! pathological overcrowding
        rnp = fnp-real(lnp,kind=rk)
        call mpi_reduce(rnp,sumrnp,1,MPI_REAL8,MPI_SUM,0,comm_cart,ierror)

        if (V_fill>0.5_rk*V_chunk) then
           send_more = 0
        else
           send_more = -1
        end if
        call mpi_gather(send_more,1,MPI_INTEGER,npp,1,MPI_INTEGER,0,comm_cart&
             &,ierror)

        if (myrankc==0) then
           ! prevent endless loop if  V_fill<=V_chunk/2  everywhere
           if (all(npp<0)) npp = 0

           do i=1,nint(sumrnp)
              do
                 call ranlux(r,1)
                 pp = floor(r(1)*nprocs)
                 if (npp(pp)>=0) exit
              end do

              npp(pp) = npp(pp)+1
           end do
        end if
        call mpi_scatter(npp,1,MPI_INTEGER,rlnp,1,MPI_INTEGER,0,comm_cart&
             &,ierror)
        if (rlnp>0) lnp = lnp+rlnp

        ! Now lnp particles have to be placed locally on each process,
        ! To reduce communication, this is done only at a time for
        ! domains which are further away from each other than
        ! min_dist. Because of the periodic boundary conditions, in
        ! each direction, first fill the lower domains up to rc in
        ! strides of size sp, then fill the remaining rp domains at
        ! the end.

        ! stride in each direction for that part of the decomposition
        ! that can be filled completely in parallel
        sp = 1+ceiling(min_dist/real((/nx,ny,nz/),kind=rk))

        ! number of remaining domains in each direction which need to
        ! be filled separately at the end of each direction
        rp = mod(cdims,sp)

        ! first process coordinate of the remainder in each direction
        rc = (cdims/sp)*sp

        do i=1,sp(1)+rp(1)
           do j=1,sp(2)+rp(2)
              do k=1,sp(3)+rp(3)
                 call borders
                 if (all(&
                      &((/i,j,k/)<=sp&
                      &.and.ccoords<rc.and.mod(ccoords,sp)==(/i,j,k/)-1)&
                      &.or.((/i,j,k/)>sp&
                      &.and.ccoords==rc+(/i,j,k/)-(sp+1)))) &
                      &call place_particles_random_local(N,lnp)
              end do
           end do
        end do

        if (myrankc==0) deallocate (npp)
    end subroutine place_particles_random

    !> calculates the volume in the local domain which could be filled
    !> by particles
    !>
    !> Counted as available volume are all fluid sites as far as they
    !> are within the range \c x0lo..x0hi .
    real(kind=rk) function fillable_volume(N)
        type(lbe_site),intent(in) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        real(kind=rk) :: V_fill,xlo(3),xhi(3),gxlo(3),gxhi(3),lxlo(3),lxhi(3)
        integer :: x,y,z

        ! global extremal positions to place particles
        xlo = max(x0lo,minpos)
        call find_highest_lattice_position(xhi)

        ! local extremal positions (still in global coordinates)
        gxlo = max(border(1,:),xlo)
        gxhi = min(border(2,:),xhi)

        ! go to local coordinates
        call local_coordinates(gxlo,lxlo)
        call local_coordinates(gxhi,lxhi)

        V_fill = 0.0_rk
        ! (/x,y,z/) could be in a halo of depth 1, so rock_state must
        ! be present there
        do x=nint(lxlo(1)),nint(lxhi(1))
           do y=nint(lxlo(2)),nint(lxhi(2))
              do z=nint(lxlo(3)),nint(lxhi(3))
                 if (N(x,y,z)%rock_state==0.0_rk) then
                    ! sum up a volume of 1x1x1 centered around each
                    ! fluid lattice site---but only as far as it is
                    ! within the bounds given by x0lo and x0hi
                    V_fill = V_fill&
                         &+product(min(lxhi,real((/x,y,z/),kind=rk)+0.5_rk)&
                         &-max(lxlo,real((/x,y,z/),kind=rk)-0.5_rk))
                 end if
              end do
           end do
        end do
        fillable_volume = V_fill
    end function fillable_volume

    !> place particles randomly in the local process domain
    !>
    !> \param[in] N local lattice chunk (halo depth >=1 required)
    !>
    !> \param[in] lnp number of particles to place
    !>
    !> The minimum distance between the new particles themselves to
    !> possibly existing particles that are either local or existing
    !> in \c P as a neighbor particle will be \c min_dist.  Particles
    !> will not be placed directly into rock or outside the interval
    !> \c \c x0lo..x0hi .
    subroutine place_particles_random_local(N,lnp)
        type(lbe_site),intent(in) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        integer,intent(in) :: lnp
        real(kind=rk) :: lxlo(3),lxhi(3),min_distsq,x(3),xhi(3),xlo(3)
        integer :: i,ii,np
        type(md_particle_type) :: new_p

        xlo(:) = max(x0lo(:),minpos(:))
        call find_highest_lattice_position(xhi)

        ! local extremal positions (in global coordinates)
        lxlo = max(border(1,:),xlo)
        lxhi = min(border(2,:),xhi)

        min_distsq = min_dist**2

        particles: do np=1,lnp
           position: do
              call ranlux(x,3)
              new_p%x = lxlo+x*(lxhi-lxlo)
              new_p%uid = -1 ! to let rdstsq work for sure

              ! don't place particles straight into rock
              if (inside_rock(N,new_p%x)) cycle position

              i = atompnt
              others: do ii=1,nlocal+nother
                 if (rdstsq(new_p,P(i))<min_distsq) cycle position
                 if (ii<=nlocal) then
                    i = list(i)
                 else
                    i = i+1
                 end if
              end do others
              exit position
           end do position
           call place_particle(new_p%x)
        end do particles
    end subroutine place_particles_random_local

    !> fills the system with particles at random positions (serial
    !> implementation)
    !>
    !> \param[in] rock_state global rock_state array with halo depth
    !> \c halo_extent, only read if a rock file is actually used
    !>
    !> \param[inout] tp will be filled with the particles created
    !> here, may be reallocated
    !>
    !> \param[inout] n_global is increased according to the number of
    !> particles created
    !>
    !> The routine takes into account \c x0lo and \x0hi. The number of
    !> particles to create is determined from the volume concentration
    !> given by \c phi in the input file. This works only with \c
    !> interaction=='ladd'. If \c phi<0, the number of sc volume cells
    !> with a lattice constant of \c alat that would fit into the
    !> volume determines the particle count.
    !>
    !> For simplicity usher algorithm is not used here, instead
    !> particles are replaced until their distance to other particles
    !> is larger than or equal to \c min_dist.
    !>
    !> \todo This is very slow/will hang for high particle densities!
    subroutine place_particles_random_serial(rock_state,tp,n_global)
        real(kind=rk),intent(in) :: rock_state&
             &(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        type(md_particle_type),allocatable,intent(inout) :: tp(:)
        integer,intent(inout) :: n_global
        integer j,np,x,y,z
        real(kind=rk) :: free_volume,xlo(3),xhi(3),mrsq,min_distsq,pos(3)
        real(kind=rk),parameter :: small=1.0E-12_rk ! prevents roundoff errors
        type(md_particle_type) :: new_p
        logical rock_site

        min_distsq = min_dist**2

        xlo(:) = max(x0lo(:),minpos(:))
        call find_highest_lattice_position(xhi)

        if (phi<0.0_rk) then
           ! number of particles is determined as for an sc-lattice
           np = nint(product((small+xhi(:)-xlo(:))/alat))
           if (np<1) call error_md&
                &('Volume concentration derived from alat would lead to less '&
                &//'than 1 particle!')
        else
           if (interaction/='ladd') call error_md&
                &('Random particle placement based on volume concentration phi'&
                &//' works only with interaction=''ladd''. Set phi<0 and use '&
                &//'alat instead!')

           free_volume = 0.0_rk
           do x=nint(xlo(1)),nint(xhi(1))
              do y=nint(xlo(2)),nint(xhi(2))
                 do z=nint(xlo(3)),nint(xhi(3))
                    ! could be that no rock file was read---too bad
                    ! Fortran has no built-in short-cut evaluation...
                    rock_site = .false.
                    if (rock_is_present()) then
                       if (rock_state(x,y,z)/=0.0_rk) rock_site = .true.
                    end if

                    if (.not.rock_site) then
                       ! sum up a volume of 1x1x1 centered around each
                       ! fluid lattice site---but only as far as it is
                       ! within the bounds given by x0lo and x0hi
                       free_volume = free_volume&
                            &+product(min(xhi,real((/x,y,z/),kind=rk)+0.5_rk)&
                            &-max(xlo,real((/x,y,z/),kind=rk)-0.5_rk))
                    end if
                 end do
              end do
           end do
           np = nint(phi*free_volume/avg_particle_volume_fluid_ladd())
           if (np<1) call error_md&
                &('Desired volume concentration phi would lead to less than 1'&
                &//' particle!')
        end if

        particles: do
           position: do
              call RANLUX(pos,3)
              new_p%x = pos*(xhi-xlo)+xlo
              new_p%uid = -1 ! to let rdstsq work for sure

              ! find minimum distance from  pos  to other particle
              if (n_global>=1) then
                 mrsq = rdstsq(new_p,tp(1))
              else
                 mrsq = min_distsq
              end if
              existing: do j=2,n_global
                 mrsq = min(mrsq,rdstsq(new_p,tp(j)))
              end do existing

              ! accept  pos  as position for new particle if all other
              ! particles are separated from it by more than  min_dist .
              if (mrsq>=min_distsq) exit
           end do position
           call place_particle_serial(rock_state,new_p%x,tp,n_global)

           ! the particle was placed only if there was no problem,
           ! e.g. with rock---exit only if enough particles are placed
           if (n_global==np) exit particles
        end do particles
    end subroutine place_particles_random_serial

!> Matches the parameters for init_cond = 11 and puts particles on the surface
!> of the droplet on the path of a spiral. Places np_sphere particles if
!> drop_Xcut and drop_Xshift are zero, but particles which would be on the part
!> of the droplet that is cut / shifted are not placed.
!> This uses the principle of a Bauer spiral.
!> 2010-06-28: Added by Stefan
!>
!> \param[in] N local lattice chunk (halo depth >=1 required)
!>
!> Particles that would be set straight into rock are not placed.
subroutine place_particles_sphere(N)
  type(lbe_site),intent(in) :: N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
  integer :: k
  real(kind=rk) :: phik, thetak, x, y, z, L, zk
  real*8 :: if1, if2 ! Real radii of the sphere
  integer, dimension(3) :: centre ! Global coords of centre of system.
  integer, dimension(3) :: base ! Global coords of corner of subdomain
  integer, dimension(3) :: offset
  real(kind=rk) :: ppos(3)

  if (np_sphere .le. 0) then
    return
    ! CALL error_md("initial_placing = <sphere> requires the parameter np_sphere to be set to a positive value.")
  endif

  ! From lbe_init_new_drop in lbe_init.F90:
  ! Define the actual radii of the spheres.
  if1 = fr1*min(tnx,tny,tnz)/2.d0
  if2 = fr2*min(tnx,tny,tnz)/2.d0

  centre = (/ tnx/2+drop_xshift, tny/2+drop_yshift, tnz/2+drop_zshift /)
  base = ccoords*(/ nx, ny, nz /)
  offset = base - centre

  L = sqrt( pi*np_sphere )
  do k=1,np_sphere

    zk = 1.0 - ( (2.0*k - 1.0) / np_sphere )
    phik = acos (zk)
    thetak = L*phik

    ! x,y,z on the unit sphere
    x = sin ( phik ) * cos ( thetak );
    y = sin ( phik ) * sin ( thetak );
    ! /* z = cos ( theta ); But z==h, so: */
    z = zk;

    ppos(1) = if1 * x + 1.0*centre(1)
    ppos(2) = if1 * y + 1.0*centre(2)
    ppos(3) = if1 * z + 1.0*centre(3)

    ! write(msgstr,"('Reporting particle ',I0,' at (',F16.10,', ',F16.10,', ',F16.10,').')") k, ppos(1), ppos(2), ppos(3)
    ! call log_msg_md(msgstr)

    ! Slight modification of the code in lbe_init_new_drop
    if ((((drop_xcut.lt.0).and.(ppos(1)+offset(1).gt.drop_xcut)).or.((drop_xcut.gt.0).and.(ppos(1)+offset(1).lt.drop_xcut)).or.(drop_xcut.eq.0)).and.&
        (((drop_ycut.lt.0).and.(ppos(2)+offset(2).gt.drop_ycut)).or.((drop_ycut.gt.0).and.(ppos(2)+offset(2).lt.drop_ycut)).or.(drop_ycut.eq.0)).and.&
        (((drop_zcut.lt.0).and.(ppos(3)+offset(3).gt.drop_zcut)).or.((drop_zcut.gt.0).and.(ppos(3)+offset(3).lt.drop_zcut)).or.(drop_zcut.eq.0))) then

      ! Only actually place the particle if it's on the domain of the core
      if ( ( ppos(1) .ge. base(1) .and. ppos(1) .lt. base(1) + nx ) .and. &
          ( ppos(2) .ge. base(2) .and. ppos(2) .lt. base(2) + ny ) .and. &
          ( ppos(3) .ge. base(3) .and. ppos(3) .lt. base(3) + nz ) ) then
        !write(msgstr,"('Placing particle ',I0,' at (',F16.10,', ',F16.10,', ',F16.10,').')") k, ppos(1), ppos(2), ppos(3)
        ! call log_msg_md(msgstr,.true.)
        ! Particles near the poles can still cause problems, so one can use the np_sphere_cap parameter to not place a number of them.
        if (k .ge. 1 + np_sphere_cap .and. k .le. np_sphere - np_sphere_cap) then
           if (.not.inside_rock(N,ppos)) call place_particle(ppos)
        end if
      end if
    end if
  end do
end subroutine place_particles_sphere

    !> fills the local domain with particles on an sc lattice
    !>
    !> \param[in] N local lattice chunk (halo depth >=1 required)
    !>
    !> Particles that would be set straight into rock are not placed.
    subroutine place_particles_sc(N)
        type(lbe_site),intent(in) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        integer i,j,k,ncells(3)
        real(kind=rk) :: p(3),xlo(3),xhi(3)
        real(kind=rk),parameter :: small=1.0E-12_rk ! prevents roundoff errors

        xlo = max(x0lo,minpos) ! global starting position
        call find_highest_lattice_position(xhi) ! global upper limit

        ! full unit cell has to fit below xhi even for the last particle
        xhi = small+xhi-alat

        ! local upper limit
        xhi = min(xhi,border(2,:))

        do k=1,3
           do
              if (xlo(k)>=border(1,k)) exit
              xlo(k) = xlo(k)+alat
           end do

           ncells(k) = 0
           do
              if (xlo(k)+alat*real(ncells(k),kind=rk)>=xhi(k)) exit
              ncells(k) = ncells(k)+1
           end do
        end do

        do k = 0,ncells(3)-1
           do j = 0,ncells(2)-1
              do i = 0,ncells(1)-1
                 p = xlo+alat*real((/i,j,k/),kind=rk)
                 if (.not.inside_rock(N,p)) call place_particle(p)
              end do
           end do
        end do
    end subroutine place_particles_sc

    !> fills the whole domain with particles on an sc lattice
    subroutine place_particles_sc_serial(rock_state,tp,n_global)
        real(kind=rk),intent(in) :: rock_state&
             &(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        type(md_particle_type),allocatable,intent(inout) :: tp(:)
        integer,intent(inout) :: n_global
        integer i,j,k,ncells(3)
        real(kind=rk) :: xlo(3),xhi(3)
        real(kind=rk),parameter :: small=1.0E-12_rk ! prevents roundoff errors

        xlo(:) = max(x0lo(:),minpos(:))
        call find_highest_lattice_position(xhi)

        ncells(:) = (small+xhi(:)-xlo(:))/alat

        do k = 0,ncells(3)-1
           do j = 0,ncells(2)-1
              do i = 0,ncells(1)-1
                 call place_particle_serial(rock_state,xlo(:)+alat*(/i,j,k/),tp&
                      &,n_global)
              enddo
           enddo
        enddo
    end subroutine place_particles_sc_serial

    !> initialize the particles with random velocities that result in the
    !> temperature  temp0
    subroutine init_velocities_temp0
        integer i,ii,ierror,n
        real(kind=rk) factor,t
        real(kind=rk) vtmp(nd),vtot(nd)

        if (polydispersity) call error_md("initial_velocities=='temp0' not "&
             &//"implemented yet for polidisperse particles---disable "&
             &//"polydispersity or set initial_velocities/='temp0'")

        call count_particles_all(n)
        if (n<2) call error_md('less than 2 particles'&
             &//' - calculating a temperature makes absolutely no sense!')

        vtot(:) = 0.0
        i = atompnt
        do ii = 1,nlocal
           call ranlux(vtmp, 3)
           P(i)%v(:) = vtmp(:)
           vtot(:) = vtot(:) + P(i)%v(:)
           i = list(i)
        enddo

        vtmp(:) = vtot(:)
        call mpi_allreduce(vtmp,vtot,3,MPI_REAL8,MPI_SUM,comm_cart,ierror)

        ! remove center of mass motion
        i = atompnt
        do ii = 1,nlocal
           P(i)%v(:) = P(i)%v(:) - vtot(:)/n
           i = list(i)
        enddo

        ! scale to requested temperature
        call temperature(t)
        factor = sqrt(temp0/t)
        i = atompnt
        do ii = 1,nlocal
           P(i)%v(:) = P(i)%v(:) * factor
           i = list(i)
        enddo
    end subroutine init_velocities_temp0



#ifdef RWALK
    !> initialize the particles with random velocities that result in the
    !> temperature  termerat
    subroutine init_velocities_tempr
        integer i,ii,ierror,l
        real(kind=rk) factor,t,k_boltz,r_v
        real(kind=rk) vtmp(nd),vtot(nd)

        k_boltz = 1.3806504E-23
        print*,"initializing random velocities"
!        if (natoms<2) call error_md('less than 2 particles'&
!             &//' - calculating a temperature makes absolutely no sense!')

        i = atompnt
        do ii = 1,nlocal
           do l=1,3
           call uni_dist(r_v)
           P(i)%v_r(l) = r_v* sqrt(k_boltz*temperat/molec_mass)*(delta_t/delta_x)
           enddo
           i = list(i)
        enddo
    end subroutine init_velocities_tempr

#endif

    !> assign  v0  as initial velocity for every particle
    subroutine init_velocities_v0
        integer i,ii

        i = atompnt
        do ii = 1,nlocal
           P(i)%v(:) = v0(:)
           i = list(i)
        enddo
    end subroutine init_velocities_v0

    !> assign  w0  as initial angular velocity for every particle
    subroutine init_rotations_w0
        integer i,ii

        i = atompnt
        do ii = 1,nlocal
           P(i)%w(:) = w0(:)
           i = list(i)
        enddo
    end subroutine init_rotations_w0

    !> assign  q0  as initial orientation for every particle
    subroutine init_orientations_q0
        integer i,ii

        i = atompnt
        do ii = 1,nlocal
           P(i)%q = q0
           i = list(i)
        enddo
    end subroutine init_orientations_q0

    !> assign random initial orientations to all particles
    subroutine init_orientations_random
        integer i,ii,k

        i = atompnt
        do ii = 1,nlocal
           do k=0,3
              call uni_dist(P(i)%q(k))
           end do
           P(i)%q = P(i)%q/sqrt(dot_product(P(i)%q,P(i)%q))

           i = list(i)
        enddo
    end subroutine init_orientations_random

    !> assign sphere normal as initial orientation for every particle
    subroutine init_orientations_sphere
        integer i,ii
        real(kind=rk), dimension(3) :: centre ! Global coords of centre of system.
        real(kind=rk), dimension(3) :: offset, dir
        real(kind=rk) :: offset_norm

        centre = (/ tnx+drop_xshift, tny+drop_yshift, tnz+drop_zshift /) /2

        i = atompnt
        do ii = 1,nlocal
          offset = P(i)%x - centre
          offset_norm = sqrt(dot_product(offset,offset))
          dir = offset / offset_norm
          P(i)%q = quaternion(dir)
          i = list(i)
        enddo
    end subroutine init_orientations_sphere

    !> Returns  .true.  if a particle with position  pos(1:3)  and orientation
    !>  ori(0:3)  would be too close to rock sites. Otherwise  .false.  is
    !> returned.
    logical function too_close_to_rock(rock_state,pos,ori)
        real(kind=rk),intent(in) :: rock_state&
             &(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        real(kind=rk),intent(in) :: pos(3),ori(0:3)

        too_close_to_rock = inside_rock_global(rock_state,pos)&
             &.or.particle_rock_potential(rock_state,pos,ori)>100.0_rk*temp0
    end function too_close_to_rock

    !> Returns in  xhi(1:3)  the highest position in every direction at which
    !> particles are placed by the initialization routines.
    subroutine find_highest_lattice_position(xhi)
        real(kind=rk),intent(out) :: xhi(3)
        integer k

        do k=1,3
           if (x0hi(k)<0.0_rk) then
              xhi(k) = maxpos(k)
           else
              xhi(k) = x0hi(k)
           end if
        end do
    end subroutine find_highest_lattice_position

#endif
end module lbe_md_init_module
