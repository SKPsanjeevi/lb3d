
#include "lbe.h"
!> main molecular dynamics module
module lbe_md_module
#ifdef MD
    use lbe_parms_module, only: nt,kbT 
    use lbe_bdist_module, only: boltz_dist
    use lbe_helper_module, only: every_n_time_steps,local_coordinates,check_NaN
    use lbe_globals_module, only: border,chunksize,halo_extent,maxpos,minpos&
         &,tsize,myrankc,pi,forcesum
    use lbe_log_module
    use lbe_md_dynamic_module, only: dissolve,insert,remove,insert_particle,delete_particle
    use lbe_md_bc_leesedwards_module, only:&
         & md_leesedwards_adapt_after_communicate&
         &,md_leesedwards_adapt_after_exchange,md_leesedwards_init_step&
         &,md_leesedwards_z_split_down,md_leesedwards_z_split_up
    use lbe_md_boundary_condition_module, only: closest_image
    use lbe_md_fluid_module, only: fluid_f_interaction_pre&
         &,fluid_ft_interaction_pre,fluid_f_interaction_post&
         &,fluid_ft_interaction_post,fluid_v_interaction,fluid_w_interaction
    use lbe_md_globals_module
    use lbe_md_fluid_ladd_module, only: in_particle
    use lbe_md_fluid_ladd_parms_module, only: R_orth,R_para
    use lbe_md_helper_module, only: count_particles_all, count_particles, error_md&
         &,fluid_velocity_and_viscosity,log_msg_md,orientation,particle_inertia&
         &,particle_mass, check_quaternion
    use lbe_md_memory_module, only: boost_npmax
    use lbe_md_output_module, only: dump_potentials,n_dump
    use lbe_md_parallel_module, only: scatter_particles,gather_particles
    use lbe_md_potential_module, only: corrections,energy,force_and_torque&
         &,pressure
    use lbe_md_magnetic_module, only: md_magnetic,force_and_torque_magnetic,force_and_torque_magnetic_inter
    use lbe_md_rand_walk_module, only: meandx,rand_walk_update_dx&
         &,rock_reflection,uni_dist
    use lbe_md_rock_module, only: rock_ft_interaction,rock_potential
    use lbe_timer_module, only: start_timer,stop_timer
    use lbe_parallel_module, only: check_allocate,cdims,comm_cart&
         &,nnprocs,nprocs,tnx,tny,tnz,owning_rankc
    use lbe_parms_module, only: boundary,inv_fluid,nt,pr,shear_u,nx,ny,nz
    use lbe_types_module, only: lbe_site
    use luxury_rng_module, only: ranlux
    use map_module, only: Mii_clear,Mii_commit,Mii_map,Mii_preinsert,Mii_exist&
         &,Mii_provide_capacity,Mii_rewind
    use lbe_md_boundary_condition_module, only: rdstsq

#ifdef DEBUG_MPI
    use lbe_parallel_module, only: log_mpi
#endif

    implicit none
    private

    public calculate_all_orientations,calculate_rotations_s&
         &,check_shear_boundary,clear_uid2i,borders,communicate,create_uid2i&
         &,exchange,fix_shear_boundary,integrate_q,integrate_v,integrate_vw&
         &,integrate_w,integrate_x,integrate_xq,list_still_valid,neighbor&
         &,reset_forces_and_torques,setup_computations,setup_list,statistics,force_gaussian_plane,random_forcekon&
         &,status,temperature

    include 'mpif.h'

    !> average particle (angular) velocities seen by the fluid over
    !> two consecutive time steps?
    logical,save,public :: average_vw_fluid=.true.

    logical,save,public :: fix_v=.false. !< hold velocity fixed?
    logical,save,public :: fix_w=.false. !< hold angular velocity fixed?
    logical,save,public :: fix_vx=.false. !< hold position in x direction fixed?
    logical,save,public :: fix_vz=.false. !< hold position in z direction fixed?
    integer,save,public :: fix_vx_t= 0
    integer,save,public :: fix_vz_t= 0
    !> \name reenter_rdm
    !> \{
    logical,save,public :: check_rdm=.false.
    logical,save,public :: check2_rdm=.false.

    integer,save,public :: pp = 0
    integer,save,public :: guid
    real(kind=rk),save,public :: gp(3)
    real(kind=rk),save,public :: go(3)
    real(kind=rk),save,public :: gv(3)
    real(kind=rk),save,public :: rntr_rdm(1)

    integer,save,public :: n_reenter = 0
    integer,save,public :: g_reenter = 0
    integer,save,public :: total_reenter = 0
    integer,save,public,allocatable :: reenter_list(:)
    integer,save,public,allocatable :: reenter_pos_list(:)
    type(md_particle_type),allocatable,dimension(:),save :: tp
    integer,save :: stat = 0
    !> \}

    !> enable ramping i.e., linear increase in particle velocity upto vmax_ramping
    logical, save,public :: ramping = .false.
  
    !> max velocity of particle prescribed to be at the end of linear ramping
    real(kind=rk), save,public :: vmax_ramping(3)=(/0.0_rk,0.0_rk,0.0_rk/)
  
    !> time to achieve vmax_ramping
    integer, save,public :: ramptime = 0 

    !> enables momentum conservation. If particle is prescribed with velocity
    !> then, force on particle is smeared through the fluid
    logical, save,public :: forcebalance = .false.

    !> constant force like gravitation
    real(kind=rk),save,public :: constant_f(3)=(/0.0_rk,0.0_rk,0.0_rk/)
    !> consant force for cargo particle
    logical, save, public :: cargo_par = .false.
    integer,save,public :: cargo_par_id  = 4
    real(kind=rk),save,public :: constant_fc(3)=(/0.0_rk,0.0_rk,0.0_rk/)
     real(kind=rk),save,public :: constant_tc(3)=(/0.0_rk,0.0_rk,0.0_rk/)
    !> constant torque
    real(kind=rk),save,public :: constant_t(3)=(/0.0_rk,0.0_rk,0.0_rk/)
    real(kind=rk),save,public :: constant_f_y_min=0.0_rk 
!gaussian potential around a plane
    logical,save,public :: gaussian_plane_active = .false.
    real(kind=rk),save,public :: gaussian_plane_pos=0.0_rk    
	real(kind=rk),save,public :: gaussian_plane_E=0.0_rk  
	real(kind=rk),save,public :: gaussian_plane_sqwidth=0.0_rk  !will lead to undefined number when this force is called and the parameter is on default
	real(kind=rk),save,public :: gaussian_plane_c=0.0_rk
  
logical,save,public :: randompartforce = .false.
!insertion / deletion in an area of the box
logical, save, public :: runtime_particle_mod = .false.
real(kind=rk),save,public :: rpm_miny=0.0_rk 
integer,save,public :: rpm_numpart=0
real(kind=rk),save,public :: rpm_minpartdist=0.0_rk 
integer,save,public :: rpm_addintervall=0
integer,save,public :: rpm_deleteintervall=0


    !> thickness of the zones at both x-boundaries of the system in
    !> which v_x is reverted
    real(kind=rk),save,public :: boundary_width=0.0_rk

    real(kind=rk),save :: dt !< 1/steps_per_lbe_step
    real(kind=rk),save :: half_dt    !< dt/2

contains

#ifdef LADD_SURR_RHOF
    !> calculate the surrounding fluid densities \c rhof for all owned
    !> particles from the accumulated values
    subroutine average_surr_rhof
        integer i,ii

        i=atompnt
        do ii=1,nlocal
           ! keep rhof of last time step if there was no single
           ! neighboring fluid site to calculate the average
           if (P(i)%n_acc>0) P(i)%rhof = P(i)%rhof_acc/real(P(i)%n_acc,rk)

           i = list(i)
        end do
    end subroutine average_surr_rhof
#endif

  subroutine gaussianBoxMuller(gaussian1,gaussian2)
    implicit none

    real(kind=rk) :: gaussian1, gaussian2
    real(kind=rk) :: u1,u2
    real(kind=rk),dimension(2) :: u

!    call ranlux(u,2)
!    u1=u(1)
!    u2=u(2)
    call random_number(u1)
    call random_number(u2)

    gaussian1 = sqrt(-2.0_rk*log(u1)) * cos(2.0_rk*pi*u2)
    gaussian2 = sqrt(-2.0_rk*log(u1)) * sin(2.0_rk*pi*u2)

  end subroutine gaussianBoxMuller

    real(kind=rk) function random_forcekon()
    real(kind=rk) :: gaussian1,gaussian2
call gaussianBoxMuller(gaussian1,gaussian2)
random_forcekon=gaussian1*kbT
end function random_forcekon

!calculate force from an attractive potential in y-plane
    real(kind=rk) function force_gaussian_plane(particle_id)
        integer,intent(in) :: particle_id
real(kind=rk) dist, distabs

dist = gaussian_plane_pos - P(particle_id)%x(2)
distabs = ABS(dist)
if(distabs < gaussian_plane_c) then

force_gaussian_plane = gaussian_plane_E * EXP(-0.5_rk*dist*dist/gaussian_plane_sqwidth )*(dist/gaussian_plane_sqwidth)
else
force_gaussian_plane = 0
endif

!if(dist < 0) then
!sign change
!force_gaussian_plane = -1.0_rk*force_gaussian_plane
!endif
 
end function force_gaussian_plane

    !> first part of leapfrog integrator
    !>
    !> \param[in,out] N local lattice chunk with halo of extent \c halo_extent
    !>
    !> \param[in] substep current md substep count (within \c
    !> 1..steps_per_lbe_step)
    !>
    !> integrate \c (v|w) to \c (x|q) and do everything else but not
    !> the integration from \c (f|t) to \c (v|w) and other things that
    !> affect \c (v|w)
    subroutine integrate_xq(N,substep)
        type(lbe_site),intent(inout) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        integer,intent(in) :: substep
        integer i,ii,k,ierror
        integer :: cp(2)
        real(kind=rk),parameter :: small = 1.0E-9_rk

        if (md_leesedwards) call md_leesedwards_init_step(N)

        call setup_random_reentering_on_demand
        cp = 0

#ifdef RWALK
        call rand_walk_update_dx
#endif

        i = atompnt
        do ii = 1,nlocal
           ! Q(T,W(T-1/2))
           if (use_rotation) call integrate_q(P(i),half_dt)

           ! X(T,V(T-1/2))
           call integrate_x(&
#ifdef RWALK
                &N,&
#endif
                &P(i),dt)

           call prepare_random_reentering_ladd_inv_fluid_11_z(i)

           do k=1,3
              if (P(i)%x(k) >= maxpos(k)) then
                 if (semiper_max(k)) then
                    call semiper(i,k,.true.)
                 else
                    ! periodic bcs
                    P(i)%x(k) = P(i)%x(k) - tsize(k)

                    if (reenter_rdm_max(k)) call mark_for_random_reentering(k,i)
                    if (count_periodic.and.k==3) cp(2) = cp(2) + 1
                 endif
              endif
              if (P(i)%x(k) < minpos(k)) then
                 if (semiper_min(k)) then
                    call semiper(i,k,.false.)
                 else
                    ! periodic bcs
                    P(i)%x(k) = P(i)%x(k) + tsize(k)

                    if (reenter_rdm_min(k)) call mark_for_random_reentering(k,i)
                    if (count_periodic.and.k==3) cp(1) = cp(1) + 1
                 endif
              end if
              ! enforce that all coordinates are inside intervals with
              ! boundaries of type [minpos|maxpos[ . Particles exactly
              ! at maxpos cause problems!
              if (P(i)%x(k)==maxpos(k)) P(i)%x(k) = P(i)%x(k) - small
           end do

           i = list(i)
        enddo

        DEBUG_MPI_MSG("After atom loop")
if(runtime_particle_mod) then 
if(every_n_time_steps(rpm_addintervall).and.substep==1) call modify_particle_add 
endif
        call insert

        if (any(reenter_rdm_min.or.reenter_rdm_max)) then
           call mpi_allreduce(n_reenter,total_reenter,1,MPI_INTEGER,MPI_SUM&
                &,comm_cart,ierror)
           n_reenter = 0
        else
           total_reenter = 0
        end if

        if (total_reenter/=0) then
          call do_random_reentering(N)
        else
           if (list_still_valid()) then
              if (communicate_rotations_s) call calculate_rotations_s
              if (substep==1) call update_vws_fluid()
              call start_timer(ti_md_comm)
              call communicate
              call stop_timer(ti_md_comm)
           else
              if (provide_uid2i) call clear_uid2i
              call remove
              call start_timer(ti_md_comm)
              call exchange
              call stop_timer(ti_md_comm)
              if (communicate_rotations_s) call calculate_rotations_s
              if (substep==1) call update_vws_fluid()
              call start_timer(ti_md_comm)
              call borders
              call stop_timer(ti_md_comm)
              if (provide_uid2i) call create_uid2i
              call start_timer(ti_md_neigh)
              call neighbor
              call stop_timer(ti_md_neigh)
           end if
        end if
if(runtime_particle_mod) then
if(every_n_time_steps(rpm_deleteintervall).and.substep==1)call modify_particle_delete
end if
!!$#ifdef PERIODIC_INFLOW
!!$        call dissolve
!!$#endif

        if (calculate_orientations) call calculate_all_orientations

        call start_timer(ti_md_force)
        call reset_forces_and_torques(substep)

        ! particle-particle interaction
        call force_and_torque
        call stop_timer(ti_md_force)

        ! magnetic force
        if (md_magnetic) then 
        call force_and_torque_magnetic
        call force_and_torque_magnetic_inter
        end if 
        !call check_quaternion
        ! give fluid the opportunity to change  P(:)%f(:)  and  P(:)%t(:)
        call start_timer(ti_md_fluid)
        if (use_rotation) then
           call fluid_ft_interaction_pre(N,substep)
        else
           call fluid_f_interaction_pre(N)
        end if
        call stop_timer(ti_md_fluid)

        ! give rock the opportunity to change  P(:)%f(:)  and  P(:)%t(:)
        call start_timer(ti_md_rock)
        call rock_ft_interaction(N)
        call stop_timer(ti_md_rock)

        if (count_periodic) &
             &call count_periodic_summary(cp,substep==steps_per_lbe_step)
    end subroutine integrate_xq

subroutine modify_particle_delete 
integer uid,ilocal,ierror,ii,partcounter
integer natoms,chosen_part
real (kind=rk)::random(1)
logical partnotfound
partnotfound=.true.
call count_particles(natoms,0)
        if (myrankc==0) then
           ! allocate buffers only once on first invocation
           if (.not.allocated(tp)) then
              allocate (tp(natoms),stat=stat)
              call check_allocate(stat,'tp')
           end if

           ! allow  natoms  to alter during the simulation
           if (natoms>size(tp)) then
             call log_msg_md(&
                   &'random reentering buffers are too small,'&
                   &//' allocating more memory...')
              deallocate (tp)
              allocate (tp(2*natoms),stat=stat)
              call check_allocate(stat,'tp')
           end if
        end if
call gather_particles(tp)

if (myrankc==0) then
partcounter=0
!grab random particle id, check if in space
do ii = 1,natoms
if(tp(ii)%x(2)>rpm_miny) then
GO TO 111
endif
end do
!print*,"found no particles to delete"
GO TO 222

111 do while(partnotfound)
call ranlux(random,1)
chosen_part=nint(natoms*random(1))
if(tp(chosen_part)%x(2)>rpm_miny) partnotfound=.false.
end do

uid = tp(chosen_part)%uid 

endif
222    call MPI_Bcast(partnotfound,1,MPI_LOGICAL,0,comm_cart,ierror)
if(.not. partnotfound) then
    call MPI_Bcast(uid,1,MPI_INTEGER,0,comm_cart,ierror)

    if(Mii_exist(uid2i,uid)) then
    ilocal = Mii_map(uid2i,uid)
    call delete_particle(P(ilocal))
    endif
endif
end subroutine modify_particle_delete

subroutine modify_particle_add 
integer uid,partcounter,ierror,loopcounter
integer:: natoms,ii
integer::newposint(3)
type(md_particle_type):: newparticle
real(kind=rk)::newpos(3)
real(kind=rk)::sqrpm_minpartdist,rpmmaxy
logical::addparticle
addparticle=.false.
sqrpm_minpartdist=rpm_minpartdist*rpm_minpartdist
partcounter=0
loopcounter=0
rpmmaxy=tsize(2)+minpos(2)-rpm_minpartdist/2.0
call count_particles(natoms,0)
        if (myrankc==0) then
           ! allocate buffers only once on first invocation
           if (.not.allocated(tp)) then
              allocate (tp(natoms),stat=stat)
              call check_allocate(stat,'tp')
           end if

           ! allow  natoms  to alter during the simulation
           if (natoms>size(tp)) then
             call log_msg_md(&
                   &'random reentering buffers are too small,'&
                   &//' allocating more memory...')
              deallocate (tp)
              allocate (tp(2*natoms),stat=stat)
              call check_allocate(stat,'tp')
           end if
        end if
call gather_particles(tp)
if (myrankc==0) then
!count particles in area of interest
	do ii = 1,natoms
		if(tp(ii)%x(2)>rpm_miny) partcounter=partcounter+1
	end do
	if(partcounter<rpm_numpart) then
!decide where to add particle
		144 call ranlux(newpos,3)
		newpos=minpos+newpos*tsize
		if(newpos(2)<rpm_miny) GO TO 144 !in allowed region
		if(newpos(2)>rpmmaxy) GO TO 144 !not touching the wall
		newparticle%x(:)=newpos

		do ii=1,natoms
			if (rdstsq(newparticle,tp(ii))<sqrpm_minpartdist)then
			loopcounter=loopcounter+1
				if(loopcounter>499) then
				print*,"found no space to insert particle, giving up"
				GO TO 222 
				else
 				GO TO 144
				end if
			end if
		end do
	addparticle=.true.
	end if
222 continue
end if !leave process on proc0

call MPI_Bcast(addparticle,1,MPI_LOGICAL,0,comm_cart,ierror)
call MPI_Bcast(newpos,3,LBE_REAL,0,comm_cart,ierror)
if(addparticle) then
newposint = nint(newpos)
if(myrankc==owning_rankc(newposint)) then
newparticle%x(:)=newpos
newparticle%v=(/0,0,0/)
newparticle%q=(/0,0,0,1/)
newparticle%qnew=(/0,0,0,1/)
newparticle%w=(/0,0,0/)
newparticle%wnew=(/0,0,0/)
newparticle%ws=(/0,0,0/)
newparticle%f_fluid=(/0,0,0/)
newparticle%t_fluid=(/0,0,0/)
newparticle%f_fluid_prev=(/0,0,0/)
newparticle%t_fluid_prev=(/0,0,0/)
newparticle%v_fluid=(/0,0,0/)
    newparticle%ws_fluid=(/0,0,0/)
newparticle%v_fluid_avg=(/0,0,0/)
newparticle%ws_fluid_avg=(/0,0,0/)
newparticle%v_fluid_acc=(/0,0,0/)
newparticle%ws_fluid_acc=(/0,0,0/)
newparticle%dx=(/0,0,0/)
    newparticle%o=(/0,0,0/)
newparticle%vnew=(/0,0,0/)
#ifdef LADD_SURR_RHOF
newparticle%rhof_acc=0
newparticle%rhof=1
     newparticle%n_acc=0
#endif
!print*,"add particle at ", newpos, " on proc ",myrankc
call insert_particle(newparticle)
endif
endif
end subroutine modify_particle_add
    !> second part of leapfrog integrator
    !>
    !> \param[in,out] N local lattice chunk with halo of extent \c halo_extent
    !>
    !> \param[in] substep current md substep count (within \c
    !> 1..steps_per_lbe_step)
    !>
    !> integrate from (f|t) to (v|w) and do other things that
    !> affect (v|w)
    subroutine integrate_vw(N,substep)
        type(lbe_site),intent(inout) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        integer,intent(in) :: substep
        integer i,ii

	! give fluid the opportunity to change  P(:)%f(:)  and  P(:)%t(:) again
        call start_timer(ti_md_fluid)
        if (use_rotation) then
           call fluid_ft_interaction_post(N,substep)
        else
           call fluid_f_interaction_post(N)
        end if
        call stop_timer(ti_md_fluid)

        ! send forces (and torques) back to owning processor if required
        call start_timer(ti_md_comm)
        if (collect_forces) call collect()
        call stop_timer(ti_md_comm)


        i = atompnt
        do ii = 1,nlocal
           ! V(T+1/2,F(T))
           if(fix_v.and.ramping) call error_md(&
             &'Both fix_v and ramping cant be true')
           if (forcebalance) then
              forcesum = forcesum+P(i)%f   
           endif
           if (.not.fix_v) call integrate_v(P(i),dt)
           ! W(T+1/2,T(T))
           if (use_rotation) call integrate_w(P(i),dt,fix_w)

           i = list(i)
        enddo

        ! give fluid the opportunity to change  P(:)%v(:)  and  P(:)%w(:)
        call start_timer(ti_md_fluid)
        if (.not.fix_v) call fluid_v_interaction(N)
        if (use_rotation.and..not.fix_w) call fluid_w_interaction(N)
        call stop_timer(ti_md_fluid)

        ! This is supposed to change the x-velocity only.
        if (boundary_width/=0.0_rk) call fix_shear_boundary
    end subroutine integrate_vw

    !> integrate the quaternion defining the orientation of a particle
    !>
    !> \param[in,out] p the particle
    !>
    !> \param[in] dt_by_2 half of the integration time step
    subroutine integrate_q(p,dt_by_2)
        type(md_particle_type),intent(inout) :: p
        real(kind=rk),intent(in) :: dt_by_2
        real(kind=rk) :: S_4x3_new(4,3)

        ! first row of matrix S is ignored because the first component of
        ! 4-vector w is constantly zero (see Allen/Tildesley, eq. (3.37)).
        S_4x3_new = reshape(&
             &(/&
             &-p%qnew(1),+p%qnew(0),+p%qnew(3),-p%qnew(2),&
             &-p%qnew(2),-p%qnew(3),+p%qnew(0),+p%qnew(1),&
             &-p%qnew(3),+p%qnew(2),-p%qnew(1),+p%qnew(0)&
             &/)&
             &,(/4,3/))
        p%q = p%q + dt_by_2*matmul(S_4x3_new,p%wnew)

        ! normalize quaternion to 1
        p%q = p%q / sqrt(dot_product(p%q,p%q))
    end subroutine integrate_q

    !> integrate the vector defining the translatory position of a particle
#ifdef RWALK
    !>
    !> \param[in] N local lattice chunk with halo of extent \c halo_extent
#endif
    !>
    !> \param[in,out] p the particle
    !>
    !> \param[in] ddt the integration time step
    subroutine integrate_x(&
#ifdef RWALK
         &N,&
#endif
         &p,ddt)
#ifdef RWALK
        type(lbe_site),intent(in) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
#endif
        type(md_particle_type),intent(inout) :: p
        real(kind=rk),intent(in) :: ddt
        real(kind=rk) :: ddx(3)

        ddx = ddt*p%vnew

#ifdef RWALK
        if (reflective_rocks) call handle_reflective_rocks(N,p)
!!$        if (p%sdx.ge.mean_free_path/delta_x) then
        if (meandx.ge.mean_free_path/delta_x) then
           call pseudo_particle_collision(N,p)
        end if
        ddx = ddx + ddt*p%v_r
        p%sdx = p%sdx + sqrt(sum(ddx**2))
        domaindx = domaindx + sqrt(sum(ddx**2))
#endif

        p%dx = p%dx + ddx
        p%x = p%x + ddx
    end subroutine integrate_x

    !> integrate the translatory velocity of a particle
    !>
    !> \param[in,out] p the particle
    !>
    !> \param[in] ddt the integration time step
    subroutine integrate_v(p,ddt)
        type(md_particle_type),intent(inout) :: p
        real(kind=rk),intent(in) :: ddt
        real(kind=rk) :: vold(3)


        vold = p%vnew
        if(.not.ramping) then
           p%vnew = vold + p%f * ddt / particle_mass(p)
          
           if (fix_vx .AND. nt < fix_vx_t) p%vnew(1)=0.0_rk 
           if (fix_vz .AND. nt < fix_vz_t) p%vnew(3)=0.0_rk
           if (p%freeze(1) == 1.0_rk) then 
           vold(1)=0.0_rk
           p%vnew(1)=0.0_rk
           end if
           if (p%freeze(2) == 1.0_rk) then
           vold(2)=0.0_rk
           p%vnew(2)=0.0_rk
           end if
           if (p%freeze(3) == 1.0_rk) then 
           vold(3) = 0.0_rk
           p%vnew(3)=0.0_rk
           end if
        else
           if(nt.le.ramptime) then
             p%vnew = vmax_ramping * (nt-1)/ramptime
           else
             p%vnew = vmax_ramping
           endif
        endif
        p%v = 0.5_rk*(vold+p%vnew)

    end subroutine integrate_v

    !> integrate the rotational velocity of a particle
    !>
    !> \param[in,out] p the particle
    !>
    !> \param[in] ddt the integration time step
    !>
    !> \param[in] ffix_w (optional, defaults to \c .false.) keep
    !> angular velocity fixed?
    subroutine integrate_w(p,ddt,ffix_w)
        type(md_particle_type),intent(inout) :: p
        real(kind=rk),intent(in) :: ddt
        logical,intent(in),optional :: ffix_w
        logical :: fffix_w
        real(kind=rk) :: A(3,3),A_new(3,3),inertia(3),l_b(3),l_b_new(3)&
             &,l_b_now(3),l_s(3),l_s_new(3),l_s_now(3),S_4X3(4,3),w_b_now(3)&
             &,wold(3)

        if (present(ffix_w)) then
           fffix_w = ffix_w
        else
           fffix_w = .false.
        end if

        inertia = particle_inertia(p)

        ! matrix to transform from space fixed to body fixed coordinates
        ! (see Allen/Tildesley, eq. (3.36)):
        A = reshape((/&
             &p%q(0)**2+p%q(1)**2-p%q(2)**2-p%q(3)**2,&
             &2.0_rk*(p%q(1)*p%q(2)-p%q(0)*p%q(3)),&
             &2.0_rk*(p%q(1)*p%q(3)+p%q(0)*p%q(2)),&

             &2.0_rk*(p%q(1)*p%q(2)+p%q(0)*p%q(3)),&
             &p%q(0)**2-p%q(1)**2+p%q(2)**2-p%q(3)**2,&
             &2.0_rk*(p%q(2)*p%q(3)-p%q(0)*p%q(1)),&

             &2.0_rk*(p%q(1)*p%q(3)-p%q(0)*p%q(2)),&
             &2.0_rk*(p%q(2)*p%q(3)+p%q(0)*p%q(1)),&
             &p%q(0)**2-p%q(1)**2-p%q(2)**2+p%q(3)**2&
             &/),&
             &(/3,3/))

        ! body-fixed angular momentum
        l_b = p%wnew * inertia

        ! space fixed angular momentum
        l_s = matmul(transpose(A),l_b)

        ! space fixed angular momentum now
        l_s_now = l_s
        if (.not.fffix_w) l_s_now = l_s_now + 0.5_rk*ddt*p%t

        ! body fixed angular momentum now
        l_b_now = matmul(A,l_s_now)

        ! body fixed angular velocity now
        w_b_now = l_b_now / inertia

        ! first row of matrix S is ignored because the first component of
        ! 4-vector w is constantly zero (see Allen/Tildesley, eq. (3.37)).
        S_4x3 = reshape(&
             &(/&
             &-p%q(1),+p%q(0),+p%q(3),-p%q(2),&
             &-p%q(2),-p%q(3),+p%q(0),+p%q(1),&
             &-p%q(3),+p%q(2),-p%q(1),+p%q(0)&
             &/)&
             &,(/4,3/))

        ! quaternion at time t+1/2dt
        p%qnew = p%q + 0.5_rk*ddt*0.5_rk*matmul(S_4x3,w_b_now)

        A_new = reshape((/&
             &p%qnew(0)**2+p%qnew(1)**2-p%qnew(2)**2-p%qnew(3)**2,&
             &2.0_rk*(p%qnew(1)*p%qnew(2)-p%qnew(0)*p%qnew(3)),&
             &2.0_rk*(p%qnew(1)*p%qnew(3)+p%qnew(0)*p%qnew(2)),&

             &2.0_rk*(p%qnew(1)*p%qnew(2)+p%qnew(0)*p%qnew(3)),&
             &p%qnew(0)**2-p%qnew(1)**2+p%qnew(2)**2-p%qnew(3)**2,&
             &2.0_rk*(p%qnew(2)*p%qnew(3)-p%qnew(0)*p%qnew(1)),&

             &2.0_rk*(p%qnew(1)*p%qnew(3)-p%qnew(0)*p%qnew(2)),&
             &2.0_rk*(p%qnew(2)*p%qnew(3)+p%qnew(0)*p%qnew(1)),&
             &p%qnew(0)**2-p%qnew(1)**2-p%qnew(2)**2+p%qnew(3)**2&
             &/),&
             &(/3,3/))

        ! space fixed angular momentum at time t+1/2dt
        l_s_new = l_s
        if (.not.fffix_w) l_s_new = l_s_new + ddt*p%t

        ! body fixed angular momentum at time t+1/2dt
        l_b_new = matmul(A_new,l_s_new)

        wold = p%wnew
        ! body fixed angular velocity at time t+1/2dt
        p%wnew = l_b_new / inertia
        p%w = 0.5_rk*(wold+p%wnew)
    end subroutine integrate_w

    !> reenter_rdm : the particle leaving the domain enter once again
    !> as periodic BC but the position is chossen randomly
    !>
    !> \param[in] k Cartesian dimension where reentering happens (\c 1..3 )
    !>
    !> \param[in] i index of the particle in \c P
    subroutine mark_for_random_reentering(k,i)
        integer,intent(in) :: k,i

        n_reenter = n_reenter + 1
        reenter_list(n_reenter) = i
        reenter_pos_list(n_reenter) = k
    end subroutine mark_for_random_reentering

    !> prepares random reentering for a particle. This is a special
    !> case for \c inv_fluid==11, \c interaction=='ladd' and the
    !> z-direction
    !>
    !> \param[in] i index of the particle in \c P
    subroutine prepare_random_reentering_ladd_inv_fluid_11_z(i)
        integer,intent(in) :: i

        if (reenter_rdm_max(3)) then
           if (interaction=='ladd'.and.inv_fluid==11) then
              if (P(i)%x(3) >= maxpos(3) - 0.5_rk - max(R_orth,R_para)) then
                 call mark_for_random_reentering(3,i)
                 P(i)%x(3) = P(i)%x(3) - tsize(3) + 2*max(R_orth,R_para) + 1.5_rk
              end if
              if (P(i)%x(3) < minpos(3) + 0.5_rk + max(R_orth,R_para)) then
                 call mark_for_random_reentering(3,i)
                 P(i)%x(3) = P(i)%x(3) + tsize(3) - 2*max(R_orth,R_para) - 1.5_rk
              end if
           end if
        end if
    end subroutine prepare_random_reentering_ladd_inv_fluid_11_z

    !> accumulate and possibly dump counters for periodic wrap-around
    !> of particles
    !>
    !> \param[in] counters number of particles that wrapped around in
    !> negative (1st element) or positive (2nd element) z direction
    !> during the current integration step
    !>
    !> \param[in] dump if \c .true., counters are dumped to standard output
    subroutine count_periodic_summary(counters,dump)
        integer,intent(in) :: counters(2)
        logical,intent(in) :: dump
        integer :: sum_counters(2),ierror
        integer,save :: count_zmax=0
        integer,save :: count_zmin=0

        call MPI_Reduce(counters,sum_counters,2,MPI_INTEGER,MPI_SUM,0,comm_cart&
             &,ierror)
        if (myrankc == 0) then
           count_zmin = count_zmin + sum_counters(1)
           count_zmax = count_zmax + sum_counters(2)
           if (dump) write (unit=6,fmt='(A,SS,I9,SP,2(A,ES15.8))') &
                &'nt=',nt,' count_zmin=',count_zmin,' count_zmax=',count_zmax
        end if
    end subroutine count_periodic_summary

    !> swap slabs of atoms with other processors in all 3 dimensions
    !> and each direction
    !>
    !> This (or equivalently \c borders()) is supposed to be done
    !> every timestep.
    subroutine communicate
        integer k,d

        DEBUG_MPI_MSG("Entered communicate")

        do k=1,3
           do d=1,n_dir_max
              call communicate_dir(sdirs(k,d))
           end do
        end do

        DEBUG_MPI_MSG("Returning from communicate")
    end subroutine communicate

    !> actually swap slabs of atoms within a given (pseudo-)direction
    !>
    !> \param[in] sd swap direction
    !>
    !> We loop until the larger one of send and recv counts while
    !> sending and recving simultanously in order to keep things
    !> synchronized. Usually, data recved will be sent further in the
    !> following iteration.
    subroutine communicate_dir(sd)
        type(swap_direction_type),intent(in) :: sd
        integer ierror,nr,reqs(2),s
        logical :: adapt

        select case (sd%dir)
        case (DIR_LE_LO_Z_DOWN,DIR_LE_LO_Z_UP,DIR_LE_HI_Z_DOWN,DIR_LE_HI_Z_UP)
           adapt = .true.
        case default
           adapt = .false.
        end select

        do s = 1,max(sd%scnt,sd%rcnt)
           nr = 0

           if (s<=sd%rcnt) then
              nr = nr+1
              call MPI_Irecv(P(sd%s(s)%rppos),sd%s(s)%rpcnt&
                   &,particle_comm_mpitype,sd%rproc,sd%dir,comm_cart,reqs(nr)&
                   &,ierror)
           end if

           if (s<=sd%scnt) then
              nr = nr+1
              call MPI_Isend(P(0),1,sd%s(s)%s_mpitype,sd%sproc,sd%dir,comm_cart&
                   &,reqs(nr),ierror)
           end if

           call MPI_Waitall(nr,reqs,MPI_STATUSES_IGNORE,ierror)

           if (s<=sd%rcnt) then
              if (md_leesedwards.and.adapt) &
                   &call md_leesedwards_adapt_after_communicate(&
                   &P(sd%s(s)%rppos:sd%s(s)%rppos+sd%s(s)%rpcnt-1))
           end if
        end do
    end subroutine communicate_dir

    !> send forces back to owning process
    !>
    !> swap forces and torques back for particles swapped by
    !> \c communicate() and add them to the forces and torques in \c P
    !>
    !> This is also used to collect the potential particle energy on
    !> the respective owning process if required by \c dump_potentials.
    subroutine collect()
        integer k,d
#ifdef BUGGYSENDINCOLLECT
        integer stat
        real(kind=rk),allocatable :: sendbuf(:)

#ifdef PARTICLESTRESS
        call error_md(&
             &'PARTICLESTRESS not yet implemented for BUGGYSENDINCOLLECT')
#endif
#ifdef LADD_SURR_RHOF
        call error_md(&
             &'LADD_SURR_RHOF not yet implemented for BUGGYSENDINCOLLECT')
#endif
        if (decouple_fluid) call error_md(&
             &'decouple_fluid not yet implemented for BUGGYSENDINCOLLECT')
        if (dump_potentials) call error_md(&
             &'dump_potentials not yet implemented for BUGGYSENDINCOLLECT')
        allocate (sendbuf(1:6*maxval(rprt_cnt)),stat=stat)
        call check_allocate(stat,'collect(): sendbuf')
#endif

        do k=3,1,-1
           do d=n_dir_max,1,-1
#ifdef BUGGYSENDINCOLLECT
              call collect_dir(sdirs(k,d),sendbuf)
#else
              call collect_dir(sdirs(k,d))
#endif
           end do
        end do
#ifdef BUGGYSENDINCOLLECT
        deallocate (sendbuf)
#endif
#ifdef LADD_SURR_RHOF
        if (nt_substep==1) call average_surr_rhof()
#endif
    end subroutine collect

    !> actually perform force[, torque][, potential energy] collection
    !> in different directions and swaps
    !>
    !> \param[in] sd swap directions; be aware that a 'collect' works
    !> in reverse order with respect to \c communicate() and how \c sd%dir
    !> is defined
#ifdef BUGGYSENDINCOLLECT
    !>
    !> \param[in] sbuf allows allocating buffer just once by caller
    subroutine collect_dir(sd,sbuf)
        real(kind=rk),allocatable,intent(inout) :: sbuf(:)
#else
    subroutine collect_dir(sd)
#endif
        type(swap_direction_type),intent(in) :: sd
        integer ft_buffer_mpitype,i,ierror,ii,nr,particle_mpitype,reqs(2),s
        logical :: energy_dump_follows

        energy_dump_follows = every_n_time_steps(n_dump).and.dump_potentials

        if (nt_substep>1) then
           ft_buffer_mpitype = ft_buffer_coll_substep_mpitype
           particle_mpitype = particle_coll_substep_mpitype
        else
           if (energy_dump_follows) then
              ft_buffer_mpitype = ft_buffer_coll_dump_mpitype
              particle_mpitype = particle_coll_dump_mpitype
           else
              ft_buffer_mpitype = ft_buffer_coll_mpitype
              particle_mpitype = particle_coll_mpitype
           end if
        end if

        do s = max(sd%scnt,sd%rcnt),1,-1
           nr = 0

           initiate_recv: if (s<=sd%scnt) then
              nr = nr+1
              call MPI_Irecv(ftbuf,sd%s(s)%sicnt,ft_buffer_mpitype,sd%sproc&
                   &,sd%dir,comm_cart,reqs(nr),ierror)
           end if initiate_recv

           initiate_send: if (s<=sd%rcnt) then ! send in reverse direction
              nr = nr+1
#ifdef BUGGYSENDINCOLLECT
              do i=sd%s(s)%rppos,sd%s(s)%rppos+sd%s(s)%rpcnt-1
                 ii = (i-sd%s(s)%rppos)*6
                 sbuf(ii+1:ii+3) = P(i)%f
                 sbuf(ii+4:ii+6) = P(i)%t
              end do
              call MPI_Isend(sbuf,6*sd%s(s)%rpcnt,MPI_REAL8,sd%rproc,sd%dir&
                   &,comm_cart,reqs(nr),ierror)
#else
              call MPI_Isend(P(sd%s(s)%rppos:sd%s(s)%rppos+sd%s(s)%rpcnt-1)&
                   &,sd%s(s)%rpcnt,particle_mpitype,sd%rproc,sd%dir,comm_cart&
                   &,reqs(nr),ierror)
#endif
           end if initiate_send

           call MPI_Waitall(nr,reqs,MPI_STATUSES_IGNORE,ierror)

           finish_recv: if (s<=sd%scnt) then
              do ii=1,sd%s(s)%sicnt
                 i = sidxs(sd%s(s)%sipos+ii-1)
                 P(i)%f = P(i)%f + ftbuf(ii)%f
#ifdef FORCECOMPONENT                 
                 P(i)%f_n = P(i)%f_n + ftbuf(ii)%f_n
                 P(i)%f_t = P(i)%f_t + ftbuf(ii)%f_t
#endif                 
                 if (use_rotation) P(i)%t = P(i)%t + ftbuf(ii)%t
                 if (nt_substep==1) then
                    if (energy_dump_follows) &
                         &P(i)%e_pot = P(i)%e_pot + ftbuf(ii)%e_pot
                    if (use_ft_fluid) then
                       P(i)%f_fluid = P(i)%f_fluid + ftbuf(ii)%f_fluid
                       P(i)%t_fluid = P(i)%t_fluid + ftbuf(ii)%t_fluid
                    end if
#ifdef LADD_SURR_RHOF
                    P(i)%rhof_acc = P(i)%rhof_acc + ftbuf(ii)%rhof_acc
                    P(i)%n_acc = P(i)%n_acc + ftbuf(ii)%n_acc
#endif
                 end if
#ifdef PARTICLESTRESS
                 P(i)%tau = P(i)%tau + ftbuf(ii)%tau
#endif
              end do
           end if finish_recv
        enddo
    end subroutine collect_dir

    !> send out atoms that have left my box, receive ones entering my
    !> box since last reneighboring for all 3 dimensions
    !>
    !> supposed to be done before every reneighboring
    subroutine exchange()
        integer k

        DEBUG_MPI_MSG("Entered exchange")

        do k=1,3
           ! conventional periodic exchange
           call exchange_dir(sdirs(k,DIR_DOWN),sdirs(k,DIR_UP))

           ! non-periodic particles (with periodic_inflow)
           call exchange_dir(&
                &sdirs(k,DIR_NONPERIODIC_DOWN),sdirs(k,DIR_NONPERIODIC_UP))

           ! possibly diagonal exchange across Lees-Edwards plane
           call exchange_dir(sdirs(k,DIR_LE_LO_Z_DOWN),sdirs(k,DIR_LE_HI_Z_UP))
           call exchange_dir(sdirs(k,DIR_LE_HI_Z_DOWN),sdirs(k,DIR_LE_LO_Z_UP))
        end do

        DEBUG_MPI_MSG("Returning from exchange")
    end subroutine exchange

    !> actually send and recv atoms for a pair of opposite directions
    !>
    !> \param[in] sdlo downward swap direction
    !>
    !> \param[in] sdhi upward swap direction
    !>
    !> \note The caller has to ensure that \c sdlo%dim==sdhi%dim and
    !> that \c sdlo%dir and \c sdhi%dir are defined as opposite
    !> directions. The check for \c freepnt!=0 should not be
    !> necessary...
    subroutine exchange_dir(sdlo,sdhi)
        type(swap_direction_type),intent(in) :: sdlo,sdhi
        real(kind=rk) :: blo,blo_z,bhi,bhi_z,pos(3)
        integer :: actrcount,curdim,i,ierror,ii,iprev,itmp,j,maxrcount,ndelete&
             &,totrcount,status(MPI_STATUS_SIZE),stype
        integer,parameter :: tag=0
        logical adapt,ignore_per,ignore_z,per

        if (sdlo%sproc==MPI_PROC_NULL.and.sdhi%sproc==MPI_PROC_NULL) return

        curdim = sdlo%dim          ! assume this is the same as sdhi%dim

        ! don't sort out particles that no cpu would take
        if (sdlo%sproc==MPI_PROC_NULL) then
           blo = -huge(blo)
        else
!!$           if (boundary=='periodic_inflow') then
!!$              if (sdlo%dir==3) then ! this implies that sdhi%dir==4
!!$                 blo = border_npi(1,curdim)
!!$              else
!!$                 blo = border_pi(1,curdim)
!!$             end if
!!$           else
              blo = border(1,curdim)
!!$           end if
        end if
        if (sdhi%sproc==MPI_PROC_NULL) then
           bhi = huge(bhi)
        else
!!$           if (boundary=='periodic_inflow') then
!!$              if (sdhi%dir==4) then
!!$                 bhi = border_npi(2,curdim)
!!$              else
!!$                 bhi = border_pi(2,curdim)
!!$              end if
!!$           else
              bhi = border(2,curdim)
!!$           end if
        end if

        if (boundary=='periodic_inflow'.and.curdim==3) then
           ignore_per = .false.
           if (sdlo%dir==3) then
              per = .false.
           else
              per = .true.
           end if
        else
           ignore_per = .true.
           per = .true.         ! set this just to please debugging tools
        end if

        ! assume that if sdlo%dir is a LE direction, also sdhi%dir is
        if (any(sdlo%dir==(/&
             &DIR_LE_LO_Z_DOWN,DIR_LE_HI_Z_DOWN,DIR_LE_LO_Z_UP,DIR_LE_HI_Z_UP/)&
             &)) then
           adapt = .true.
           ignore_z = .false.
           if (sdlo%sproc/=MPI_PROC_NULL) then
              ! this means that the 1st sendrecv will send
              blo_z = sdlo%blo_z
              bhi_z = sdlo%bhi_z
           else
              ! this means that the 2nd sendrecv will send
              blo_z = sdhi%blo_z
              bhi_z = sdhi%bhi_z
           end if
        else
           adapt = .false.
           ! non LE-direction: don't care about blo_z and bhi_z
           ignore_z = .true.
           ! just don leave these undefined
           blo_z = huge(blo_z)
           bhi_z = -huge(bhi_z)
        end if

        ndelete = 0
        iprev = 0
        j = 0
        i = atompnt

        ! index particles leaving my box, update local list
        do ii = 1,nlocal
           ! if P(i)%x was wrapped around a periodic boundary, wrap it
           ! back here because the question is which way a particle
           ! left. This could be optimized by wrapping back only the
           ! component(s) that are actually needed for the test.
           pos = closest_image(P(i)%x,border(1,:)+0.5_rk*chunksize)
!!$           if ((ignore_per.or.(is_periodic(P(i)).eqv.per))&
           if ((ignore_per)&
                &.and.(pos(curdim)<blo.or.pos(curdim)>=bhi)&
                &.and.(ignore_z.or.(pos(3)>=blo_z.and.pos(3)<bhi_z))) then
              j = j + 1
              if (j<=nemax) sidxs(j) = i

              ndelete = ndelete + 1
              if (iprev.eq.0) then
                 atompnt = list(i)
              else
                 list(iprev) = list(i)
              endif
              itmp = list(i)
              list(i) = freepnt
              freepnt = i
              i = itmp
           else
              iprev = i
              i = list(i)
           endif
        enddo
        nlocal = nlocal-ndelete
        nexcmax = max(nexcmax,j)

        if (j>nfmax.or.j>nemax) then
           write (6,*) myrankc,curdim,blo,bhi,'j=',j,',nfmax=',nfmax&
                &,',nemax=',nemax
           call error_md('Sending too many exchange atoms. '&
                &//'(Boost nfmax or nemax to values>j !)')
        endif

        call mpi_type_indexed(j,slens,sidxs,particle_exch_mpitype,stype,ierror)
        call mpi_type_commit(stype,ierror)

        ! send them out in both directions (if neighboring nodes are
        ! different)
        totrcount = 0
        maxrcount = nfmax
        call mpi_sendrecv(P(0),1,stype,sdlo%sproc,tag&
             &,pbuf(1),maxrcount,particle_exch_mpitype,sdhi%sproc,tag&
             &,comm_cart,status,ierror)
        call mpi_get_count(status,particle_exch_mpitype,actrcount,ierror)
        totrcount = totrcount + actrcount

        if (sdlo%sproc/=sdhi%sproc) then
           maxrcount = nfmax - totrcount
           call mpi_sendrecv(P(0),1,stype,sdhi%sproc,tag&
                &,pbuf(totrcount+1),maxrcount,particle_exch_mpitype,sdlo%sproc&
                &,tag&
                &,comm_cart,status,ierror)
           call mpi_get_count(status,particle_exch_mpitype,actrcount,ierror)
           totrcount = totrcount + actrcount
        end if

        call mpi_type_free(stype,ierror)

        if (md_leesedwards.and.adapt) &
             &call md_leesedwards_adapt_after_exchange(pbuf(1:totrcount))
        sort_incoming: do j = 1,totrcount
           ! check incoming atoms to see if they are in my box (could
           ! be in node's box on other side of the sender)
           if (pbuf(j)%x(curdim)>=blo.and.pbuf(j)%x(curdim)<bhi) then
              if (nlocal+1>npmax) call boost_npmax(nlocal+1)
              nlocal = nlocal + 1
              ! The case freepnt==0 should never happen thanks
              ! to boost_npmax()
              if (freepnt.ne.0) then
                 itmp = atompnt
                 atompnt = freepnt
                 freepnt = list(freepnt)
                 list(atompnt) = itmp
                 ! assign only elements from  particle_exch_mpitype ???
                 P(atompnt) = pbuf(j)
              endif
           endif
        enddo sort_incoming
    end subroutine exchange_dir

    !> Make lists of nearby atoms to send to neighboring nodes and
    !> actually send them.
    !>
    !> For each dimensions and direction, first, one list is made for
    !> every send that will be made as list is made, then actually do
    !> swaps. This does equivalent of a \c communicate() (so don't
    !> need to explicitly call \c communicate() routine on
    !> reneighboring timestep). This routine is called before every
    !> reneighboring.
    subroutine borders
        integer k,d,nsend

        DEBUG_MPI_MSG("Entered borders")

        nother = 0              ! recalculate this in the following
        nsend = 0

        do k=1,3
           do d=1,n_dir_max
              call borders_dir(sdirs(k,d),nsend)
           end do
        end do

        if (collect_forces.and.size(ftbuf)<nswpmax) then
           deallocate (ftbuf)
           allocate (ftbuf(2*nswpmax),stat=stat)
           call check_allocate(stat,'borders(): ftbuf')
        end if

        nslistmax = max(nslistmax,nsend) ! just for statistics

        DEBUG_MPI_MSG("Returning from borders")
    end subroutine borders

    !> Actually make lists of nearby atoms to send to neighboring
    !> nodes for a given dimension and direction and then send them.
    !>
    !> \param[in] sd swap direction
    !>
    !> \param[in,out] nsend particles marked for sending so far
    !>
    !> We loop until the larger one of send and recv counts while
    !> sending and recving simultanously in order to keep things
    !> synchronized. Usually, data recved will be sent further in the
    !> following iteration according to \c sbndlo and \c sbndhi. \c
    !> s_mpitype is set up for usage by \c communicate() in th
    !> following time steps.
    subroutine borders_dir(sd,nsend)
        type(swap_direction_type),intent(inout) :: sd
        integer,intent(inout) :: nsend
        integer i,ierror,ii,maxrcount,nr,reqs(2),s
        integer statuses(MPI_STATUS_SIZE,2)
        logical :: adapt,ignore_per,ignore_z,per

        if (boundary=='periodic_inflow'.and.sd%dim==3) then
           ignore_per = .false.
           select case (sd%dir)
           case (3,4)
              per = .false.
           case default
              per = .true.
           end select
        else
           ignore_per = .true.
           per = .true.         ! set this just to please debugging tools
        end if

        select case(sd%dir)
        case (DIR_LE_LO_Z_DOWN,DIR_LE_LO_Z_UP,DIR_LE_HI_Z_DOWN,DIR_LE_HI_Z_UP)
           adapt = .true.
           ignore_z = .false.
        case default
           adapt = .false.
           ignore_z = .true.
        end select

        do s = 1,max(sd%scnt,sd%rcnt)
           nr = 0 ! number of MPI requests for this send (finally either 1 or 2)

           initiate_recv: if (s<=sd%rcnt) then
              nr = nr+1

              ! put incoming particles at the end of my P

              ! starting position for this recv
              sd%s(s)%rppos = npmax+nother+1

              maxrcount = nomax-nother
              call MPI_Irecv(P(sd%s(s)%rppos),maxrcount,particle_comm_mpitype&
                   &,sd%rproc,sd%dir,comm_cart,reqs(nr),ierror)
           end if initiate_recv

           initiate_send: if (s<=sd%scnt) then
              nr = nr+1

              ! starting position for this send
              sd%s(s)%sipos = nsend + 1

              i = atompnt
              do ii = 1,nlocal+nother
!!$                 if ((ignore_per.or.(is_periodic(P(i)).eqv.per)).and.&
                 if (ignore_per.and.&
                      &P(i)%x(sd%dim)>=sd%s(s)%sbndlo.and.&
                      &P(i)%x(sd%dim)<sd%s(s)%sbndhi.and.&
                      &(ignore_z.or.&
                      &(P(i)%x(3)>=sd%blo_z.and.P(i)%x(3)<sd%bhi_z))) &
                      &then
                    nsend = nsend+1
                    if (nsend<=nemax) sidxs(nsend) = i
                 end if
                 if (ii<=nlocal) then
                    i = list(i)
                 else
                    i = i + 1
                 end if
              end do

              if (nsend>nemax) then
                 write(msgstr,&
                      &"('dim=',I0,', dir=',I0,', send=',I0,', nsend=',I0,"&
                      &//"', nemax=',I0)") sd%dim,sd%dir,s,nsend,nemax
                 call log_msg_md(msgstr,.true.)
                 call error_md('Too many atoms in border list. Boost nemax!')
              end if

              sd%s(s)%sicnt = nsend+1-sd%s(s)%sipos

              if (sd%s(s)%s_mpitype/=MPI_DATATYPE_NULL)&
                   & call mpi_type_free(sd%s(s)%s_mpitype,ierror)
              call mpi_type_indexed(sd%s(s)%sicnt,slens(1:sd%s(s)%sicnt)&
                   &,sidxs(sd%s(s)%sipos:sd%s(s)%sipos+sd%s(s)%sicnt-1)&
                   &,particle_comm_mpitype,sd%s(s)%s_mpitype,ierror)
              call mpi_type_commit(sd%s(s)%s_mpitype,ierror)

              nswpmax = max(nswpmax,sd%s(s)%sicnt)

              call MPI_Isend(P(0),1,sd%s(s)%s_mpitype,sd%sproc,sd%dir,comm_cart&
                   &,reqs(nr),ierror)
           end if initiate_send

           call MPI_Waitall(nr,reqs,statuses,ierror)

           finish_recv: if (s<=sd%rcnt) then
              call MPI_Get_count(statuses(:,1),particle_comm_mpitype&
                   &,sd%s(s)%rpcnt,ierror)
              nother = nother + sd%s(s)%rpcnt

              if (md_leesedwards.and.adapt) &
                   &call md_leesedwards_adapt_after_communicate(&
                   &P(sd%s(s)%rppos:sd%s(s)%rppos+sd%s(s)%rpcnt-1))

              if (nother>=nomax) then
                 ! not sure whether this is reached at all... probably
                 ! a "message truncated" error is issued before
                 ! already due to maxrcount being smaller than the
                 ! data recv'ed.
                 write(msgstr,&
                      &"('dim=',I0,', dir=',I0,', send=',I0,', nother=',I0,"&
                      &//"', nomax=',I0)") sd%dim,sd%dir,s,nother,nomax
                 call log_msg_md(msgstr,.true.)
                 call error_md('Received too many border atoms. Boost nomax!')
              end if
           end if finish_recv
        enddo
    end subroutine borders_dir

    !> driver for neighbor-list creation
    subroutine neighbor
        integer i,ii

        DEBUG_MPI_MSG("Entered neighbor")

        if (ineigh.eq.0) then
           call neighbor0
        else
           call neighbor1
        endif

        i = atompnt
        do ii = 1,nlocal
           P(i)%dx(:) = 0.0
           i = list(i)
        end do

        mneigh = mneigh + 1
        nlocalmax = max(nlocalmax,nlocal)
        nothermax = max(nothermax,nother)
        neighmax = max(neighmax,nnlist(nlocal+1)-1)

        DEBUG_MPI_MSG("Returning from neighbor")
    end subroutine neighbor

    !> neighbor list creation without binning
    !>
    !> Newton's 3rd law N^2 / 2 search for neighbor pairs in
    !> my box pair stored once in list if atoms i AND j are in my box
    !> (and i < j) pair stored by me if j is NOT in my box (also
    !> stored by node owning j)
    subroutine neighbor0
        real(kind=rk) xtmp,ytmp,ztmp,delx,dely,delz,rsq
        integer npnt,i,j,ii,jj

        npnt = 1
        i = atompnt
        do ii = 1,nlocal
           nnlist(ii) = npnt
           xtmp = P(i)%x(1)
           ytmp = P(i)%x(2)
           ztmp = P(i)%x(3)
           j = list(i)
           do jj = ii+1,nlocal+nother
              delx = xtmp - P(j)%x(1)
              dely = ytmp - P(j)%x(2)
              delz = ztmp - P(j)%x(3)
              ! ATTENTION: optimized minimum image criterion - result is
              ! useless if it is greater than cutsq2
              if (abs(delx).gt.tsize_mrs(1)) then
                 if (delx.lt.0.0) then
                    delx = delx + tnx
                 else
                    delx = delx - tnx
                 endif
              endif
              if (abs(dely).gt.tsize_mrs(2)) then
                 if (dely.lt.0.0) then
                    dely = dely + tny
                 else
                    dely = dely - tny
                 endif
              endif
              if (abs(delz).gt.tsize_mrs(3)) then
                 if (delz.lt.0.0) then
                    delz = delz + tnz
                 else
                    delz = delz - tnz
                 endif
              endif
              rsq = delx*delx + dely*dely + delz*delz
              if (rsq.le.cutsq2) then
                 nlist(npnt) = j
                 npnt = npnt + 1
              endif
              if (jj.le.nlocal) then
                 j = list(j)
              else
                 j = j + 1
              endif
           enddo
           if (npnt>size(nlist)-nnmax) then
              write (6,*) myrankc,npnt,ii,nlocal
              call error_md('Neighbor list too big.')
           endif
           i = list(i)
        enddo
        nnlist(nlocal+1) = npnt
    end subroutine neighbor0

    !> neighbor list creation with binning
    !>
    !> Newton's 3rd law all of mine and nearby atoms binned
    !> once each owned atom i checks 27 surrounding boxes pair stored
    !> once in list if atoms i AND j are in my box (and i < j) pair
    !> stored by me if j is NOT in my box (also stored by node owning
    !> j)
    subroutine neighbor1
        integer i,j,k,ii,ix,iy,iz,npnt,ib,ixx,iyy,izz
        real(kind=rk) xtmp,ytmp,ztmp,delx,dely,delz,rsq

        if (boundary=='periodic_inflow') &
             &call error_md('Currently, BC periodic_inflow requires ineigh=0')

        do i = 1,mbinx*mbiny*mbinz
           binpnt(i) = 0
        enddo

        i = atompnt
        do ii = 1,nlocal+nother
           ix = (P(i)%x(1) - minx) / binsizex
           iy = (P(i)%x(2) - miny) / binsizey
           iz = (P(i)%x(3) - minz) / binsizez
           ix = ix - mbinxlo
           if (ix.lt.0) ix = ix + nbinx
           iy = iy - mbinylo
           if (iy.lt.0) iy = iy + nbiny
           iz = iz - mbinzlo
           if (iz.lt.0) iz = iz + nbinz
           if (ix.le.mbinx.and.iy.le.mbiny.and.iz.le.mbinz) then
              ib = iz*mbiny*mbinx + iy*mbinx + ix + 1
              bin(i) = binpnt(ib)
              binpnt(ib) = i
           endif
           if (ii.le.nlocal) then
              i = list(i)
           else
              i = i + 1
           endif
        enddo

        npnt = 1
        i = atompnt
        do ii = 1,nlocal
           nnlist(ii) = npnt
           xtmp = P(i)%x(1)
           ytmp = P(i)%x(2)
           ztmp = P(i)%x(3)
           ixx = (xtmp - minx) / binsizex
           iyy = (ytmp - miny) / binsizey
           izz = (ztmp - minz) / binsizez
           do k = 0,26
              ix = ixx + mod(k,3) - 1
              iy = iyy + mod(k/3,3) - 1
              iz = izz + k/9 - 1
              ix = ix - mbinxlo
              if (ix.lt.0) ix = ix + nbinx
              if (ix.eq.nbinx) ix = 0
              iy = iy - mbinylo
              if (iy.lt.0) iy = iy + nbiny
              if (iy.eq.nbiny) iy = 0
              iz = iz - mbinzlo
              if (iz.lt.0) iz = iz + nbinz
              if (iz.eq.nbinz) iz = 0
              ib = iz*mbiny*mbinx + iy*mbinx + ix + 1
              j = binpnt(ib)
30            if (j.ne.0) then
                 if (j.le.i) goto 40
                 delx = xtmp - P(j)%x(1)
                 dely = ytmp - P(j)%x(2)
                 delz = ztmp - P(j)%x(3)
                 ! ATTENTION: optimized minimum image criterion - result is
                 ! useless if it is greater than cutsq2
                 if (abs(delx).gt.tsize_mrs(1)) then
                    if (delx.lt.0.0) then
                       delx = delx + tnx
                    else
                       delx = delx - tnx
                    endif
                 endif
                 if (abs(dely).gt.tsize_mrs(2)) then
                    if (dely.lt.0.0) then
                       dely = dely + tny
                    else
                       dely = dely - tny
                    endif
                 endif
                 if (abs(delz).gt.tsize_mrs(3)) then
                    if (delz.lt.0.0) then
                       delz = delz + tnz
                    else
                       delz = delz - tnz
                    endif
                 endif
                 rsq = delx*delx + dely*dely + delz*delz
                 if (rsq.le.cutsq2) then
                    nlist(npnt) = j
                    npnt = npnt + 1
                 endif
40               j = bin(j)
                 goto 30
              endif
           enddo
           if (npnt>size(nlist)-nnmax) then
              write (6,*) myrankc,npnt,ii,nlocal
              call error_md('Neighbor list too big.')
           endif
           i = list(i)
        enddo
        nnlist(nlocal+1) = npnt
    end subroutine neighbor1

    !> calculate orientations for all own and halo'ed particles
    subroutine calculate_all_orientations
        integer i,ii

        i = atompnt
        do ii = 1,nlocal
           P(i)%o = orientation(P(i)%q)
           i = list(i)
        end do
        do i=npmax+1,npmax+nother
           P(i)%o = orientation(P(i)%q)
        end do
    end subroutine calculate_all_orientations

    !> calculate angular velocity in space fixed coordinates for own
    !> particles
    subroutine calculate_rotations_s
        integer :: i,ii
        real(kind=rk) :: At(3,3)

        i = atompnt
        do ii = 1,nlocal
           At = transpose(reshape((/&
                &P(i)%q(0)**2+P(i)%q(1)**2-P(i)%q(2)**2-P(i)%q(3)**2,&
                &2.0_rk*(P(i)%q(1)*P(i)%q(2)-P(i)%q(0)*P(i)%q(3)),&
                &2.0_rk*(P(i)%q(1)*P(i)%q(3)+P(i)%q(0)*P(i)%q(2)),&

                &2.0_rk*(P(i)%q(1)*P(i)%q(2)+P(i)%q(0)*P(i)%q(3)),&
                &P(i)%q(0)**2-P(i)%q(1)**2+P(i)%q(2)**2-P(i)%q(3)**2,&
                &2.0_rk*(P(i)%q(2)*P(i)%q(3)-P(i)%q(0)*P(i)%q(1)),&

                &2.0_rk*(P(i)%q(1)*P(i)%q(3)-P(i)%q(0)*P(i)%q(2)),&
                &2.0_rk*(P(i)%q(2)*P(i)%q(3)+P(i)%q(0)*P(i)%q(1)),&
                &P(i)%q(0)**2-P(i)%q(1)**2-P(i)%q(2)**2+P(i)%q(3)**2&
                &/),&
                &(/3,3/)))
           P(i)%ws = matmul(At,P(i)%w)

           i = list(i)
        end do
    end subroutine calculate_rotations_s

    !> calculate every particle's new \c v_fluid, \c ws_fluid, \c
    !> v_fluid_avg, and \c ws_fluid_avg on the respective owning node
    !>
    !> Calculate owned particles' \c v_fluid and \c ws_fluid, at the
    !> beginning of a new LB step, that are, the averaged values of
    !> the current and the last \c steps_per_lbe_step-1 sub-steps and
    !> reset accumulators.
    !>
    !> Calculate \c v_fluid_avg and \c ws_fluid_avg as the average of
    !> the previous and the new \c v_fluid and \c ws_fluid,
    !> respectively.
    subroutine update_vws_fluid()
        integer i,ii
        real(kind=rk) :: v_fluid_new(3),ws_fluid_new(3)

        i = atompnt
        if (decouple_fluid) then
           if (average_vw_fluid) then
              do ii=1,nlocal
                 v_fluid_new = &
                      &(P(i)%v+P(i)%v_fluid_acc)/real(steps_per_lbe_step,kind=rk)
                 P(i)%v_fluid_avg = 0.5_rk*(P(i)%v_fluid+v_fluid_new)
                 P(i)%v_fluid = v_fluid_new
                 P(i)%v_fluid_acc = 0.0_rk

                 ws_fluid_new = &
                      &(P(i)%ws+P(i)%ws_fluid_acc)/real(steps_per_lbe_step,kind=rk)
                 P(i)%ws_fluid_avg = 0.5_rk*(P(i)%ws_fluid+ws_fluid_new)
                 P(i)%ws_fluid = ws_fluid_new
                 P(i)%ws_fluid_acc = 0.0_rk

                 i = list(i)
              end do
           else
              do ii=1,nlocal
                 P(i)%v_fluid = &
                      &(P(i)%v+P(i)%v_fluid_acc)/real(steps_per_lbe_step,kind=rk)
                 P(i)%v_fluid_avg = P(i)%v_fluid
                 P(i)%v_fluid_acc = 0.0_rk

                 P(i)%ws_fluid = &
                      &(P(i)%ws+P(i)%ws_fluid_acc)/real(steps_per_lbe_step,kind=rk)
                 P(i)%ws_fluid_avg = P(i)%ws_fluid
                 P(i)%ws_fluid_acc = 0.0_rk

                 i = list(i)
              end do
           end if
        else
           ! simplified version of the above for steps_per_lbe_step==1
           if (average_vw_fluid) then
              do ii=1,nlocal
                 v_fluid_new = P(i)%v
                 P(i)%v_fluid_avg = 0.5_rk*(P(i)%v_fluid+v_fluid_new)
                 P(i)%v_fluid = v_fluid_new

                 ws_fluid_new = P(i)%ws
                 P(i)%ws_fluid_avg = 0.5_rk*(P(i)%ws_fluid+ws_fluid_new)
                 P(i)%ws_fluid = ws_fluid_new

                 i = list(i)
              end do
           else
              do ii=1,nlocal
                 P(i)%v_fluid = P(i)%v
                 P(i)%v_fluid_avg = P(i)%v_fluid

                 P(i)%ws_fluid = P(i)%ws
                 P(i)%ws_fluid_avg = P(i)%ws_fluid

                 i = list(i)
              end do
           end if
        end if
    end subroutine update_vws_fluid

    !> Checks for neighbor list validity
    !>
    !> checks whether two particles together have moved farther than
    !> the skin radius since the last list update
    !>
    !> \note The second condition \c displ1<0.5_rk is probably not
    !> enough for particles with finite size and it is placed here
    !> somewhat unsystematically and it is unnecessary for particles
    !> that do not interact with the lbe sites
    logical function list_still_valid()
      integer i,ii,ierror
      real(kind=rk) displ,displtmp
      real(kind=rk) displ1,displ2    ! largest and 2nd-largest displacement
      logical manual_trigger
      logical :: displ1_NaN = .false.
      logical :: displ2_NaN = .false.

      call start_timer(ti_md_comm)
      call mpi_allreduce(list_update_required,manual_trigger,1,MPI_LOGICAL&
           &,MPI_LOR,MPI_COMM_WORLD,ierror)
      call stop_timer(ti_md_comm)
      list_update_required = .false.

      if (.not.manual_trigger) then
         displ1 = 0.0
         displ2 = 0.0

         i = atompnt
         do ii = 1,nlocal
            displ = sqrt(dot_product(P(i)%dx(:),P(i)%dx(:)))
            if (displ >= displ1) then
               displ2 = displ1
               displ1 = displ
            else if (displ >= displ2) then
               displ2 = displ
            end if

            i = list(i)
         end do

         ! Having NaNs here is very dangerous, because it could make
         ! different CPUs take different code paths, resulting in MPI
         ! blocking.  So check for them, and commit seppuku if there's
         ! a NaN.

         if (check_NaN(displ1)) displ1_NaN = .true.
         if (check_NaN(displ2)) displ2_NaN = .true.

         if ( displ1_NaN .or. displ2_NaN ) then
            call error_md("NaN detected in list_still_valid")
         endif

         !   search the two largest displacements on all processors
         call start_timer(ti_md_comm)
         displtmp = displ1
         call MPI_Allreduce(displtmp,displ1,1,MPI_REAL8,MPI_MAX,comm_cart,&
              &ierror)
         if (displtmp == displ1) then
            displtmp = displ2
         end if
         call MPI_Allreduce(displtmp,displ2,1,MPI_REAL8,MPI_MAX,comm_cart,&
              &ierror)
         call stop_timer(ti_md_comm)
      end if

      ! second condition is for tracer particles: If a particle
      ! moved further than 0.5 lattice units out of the local box
      ! then the interpolation of the velocity field at its position
      ! would require information outside the halo.
      list_still_valid = (displ1 + displ2 + list_update_displ_add <= rs - rc)&
           &.and.(displ1 + list_update_displ_add < 0.5_rk)&
           &.and..not.manual_trigger

      ! reset accumulator for additional displacement if the lists
      ! were invalidated (so they will be rebuilt)
      if (.not.list_still_valid) list_update_displ_add = 0.0_rk
    end function list_still_valid

    !> reset local forces and torques, apply constant force/torque if set
    !>
    !> \param[in] substep current md substep count (within \c
    !> 1..steps_per_lbe_step)
    !>
    !> Also \c e_pot is reset here when this is necessary.
#ifdef LADD_SURR_RHOF
    !> Also \c rhof_acc and \c n_acc are reset here when this is
    !> necessary.
#endif
    subroutine reset_forces_and_torques(substep)
        integer,intent(in) :: substep
        integer i,ii
        logical :: energy_dump_follows

        energy_dump_follows = every_n_time_steps(n_dump).and.dump_potentials

        if (substep==1.and.use_ft_fluid) call reset_fluid_forces_and_torques

        i = atompnt
        if (use_rotation) then
           do ii = 1,nlocal
              if(P(i)%x(2)>constant_f_y_min) then
             !  print*,"particle id", P(i)%uid    
            if (cargo_par .AND. P(i)%uid==cargo_par_id) then
                P(i)%f = constant_fc
        			 else   
                P(i)%f = constant_f
					 end if            
   else
                  P(i)%f = 0.0_rk
#ifdef FORCECOMPONENT
                  P(i)%f_normal  = 0.0_rk
                  P(i)%f_tangent = 0.0_rk
                  P(i)%f_n = 0.0_rk
                  P(i)%f_t = 0.0_rk
#endif
              endif
              if(P(i)%x(2)>constant_f_y_min) then
                 if (cargo_par .AND. P(i)%uid==cargo_par_id) then
                P(i)%t = constant_tc
            else   
               P(i)%t = constant_t
            end if
              else
                  P(i)%t = 0.0_rk
              endif
	      if(gaussian_plane_active) then
		P(i)%f(2)=P(i)%f(2)+force_gaussian_plane(i)
	      endif
if(randompartforce) then
	P(i)%f(1)=P(i)%f(1)+random_forcekon()
	P(i)%f(2)=P(i)%f(2)+random_forcekon()
	P(i)%f(3)=P(i)%f(3)+random_forcekon()
endif
              if (energy_dump_follows) P(i)%e_pot = 0.0_rk
#ifdef PARTICLESTRESS
              P(i)%tau = 0.0_rk
#endif
#ifdef LADD_SURR_RHOF
              P(i)%rhof_acc = 0.0_rk
              P(i)%n_acc = 0
#endif
!#ifdef MD_MAG
              P(i)%fmi = 0.0_rk
              P(i)%tmi = 0.0_rk
!#endif
              P(i)%freeze = 0.0_rk
              i = list(i)
           end do
        else
           do ii = 1,nlocal
              if(P(i)%x(2)>constant_f_y_min) then
                  		 if (cargo_par .AND. P(i)%uid== cargo_par_id) then
               		 P(i)%f = constant_fc
							else
               		   P(i)%f = constant_f
								end if   
            else
                  P(i)%f = 0.0_rk
              endif
	      if(gaussian_plane_active) then
		P(i)%f(2)=P(i)%f(2)+force_gaussian_plane(i)
	      endif
if(randompartforce) then
	P(i)%f(1)=P(i)%f(1)+random_forcekon()
	P(i)%f(2)=P(i)%f(2)+random_forcekon()
	P(i)%f(3)=P(i)%f(3)+random_forcekon()
endif
              if (energy_dump_follows) P(i)%e_pot = 0.0_rk
#ifdef PARTICLESTRESS
              P(i)%tau = 0.0_rk
#endif
#ifdef LADD_SURR_RHOF
              P(i)%rhof_acc = 0.0_rk
              P(i)%n_acc = 0
#endif
!#ifdef MD_MAG
              P(i)%fmi = 0.0_rk
              P(i)%tmi = 0.0_rk
!#endif
              P(i)%freeze = 0.0_rk
              i = list(i)
           end do
        end if

        if (collect_forces) then
           do i=npmax+1,npmax+nother
              P(i)%f = 0.0_rk
              if (use_rotation) P(i)%t = 0.0_rk
              if (energy_dump_follows) P(i)%e_pot = 0.0_rk
#ifdef PARTICLESTRESS
              P(i)%tau = 0.0_rk
#endif
#ifdef LADD_SURR_RHOF
              P(i)%rhof_acc = 0.0_rk
              P(i)%n_acc = 0
#endif
!#ifdef MD_MAG
              P(i)%fmi = 0.0_rk
              P(i)%tmi = 0.0_rk
!#endif
           end do
        end if
    end subroutine reset_forces_and_torques

    !> resets \c f_fluid and \c t_fluid to \c 0 for owned and halo'ed
    !> particles
    subroutine reset_fluid_forces_and_torques
        integer i,ii

        i = atompnt
        do ii = 1,nlocal+nother
           if (ii<=nlocal) then
              P(i)%f_fluid_prev = P(i)%f_fluid
              P(i)%t_fluid_prev = P(i)%t_fluid
           end if

           P(i)%f_fluid = 0.0_rk
           P(i)%t_fluid = 0.0_rk
#ifdef FORCECOMPONENT
           P(i)%f_normal  = 0.0_rk
           P(i)%f_tangent = 0.0_rk
           P(i)%f_n = 0.0_rk
           P(i)%f_t = 0.0_rk
#endif
#ifdef LADD_SSD
           P(i)%FvF = 0.0_rk
           P(i)%FvT = 0.0_rk
           P(i)%FwF = 0.0_rk
           P(i)%FwT = 0.0_rk
#endif

           if (ii<=nlocal) then
              i = list(i)
           else
              i = i+1
           end if
        end do
    end subroutine reset_fluid_forces_and_torques

    !> initialize constants that simplify the computations
    subroutine setup_computations
        dt = 1.0_rk/real(steps_per_lbe_step,kind=rk)
        if (use_rotation) then
           half_dt = 0.5_rk*dt
        end if
    end subroutine setup_computations

    !> thermodynamic computations
    subroutine status(N)
        type(lbe_site),intent(inout) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        real(kind=rk) :: eng,m(3),p,rpot,t

        mstat = mstat + 1
        if (mstat > ntmax) then
           write (6,fmt='(a,i4,a,i6,a,i6)')&
                & 'myrankc=',myrankc,', ntmax=',ntmax,', mstat=',mstat
           call error_md('ntmax too small')
        end if

        ! this should be ok...
        call translational_kinetic_energy(e_trans(mstat))
        if (use_rotation) then
           call rotational_kinetic_energy(e_rot(mstat))
        else
           e_rot(mstat) = 0.0_rk
        end if

        ! this is probably not correct...
        call temperature(t)
        call energy(eng)
        call rock_potential(N,rpot)
        call pressure(p,t)
        call momentum(m)
        tmparr(mstat) = t
        engarr(mstat) = eng
        rpotarr(mstat) = rpot
        prsarr(mstat) = p
        momentumarr(:,mstat) = m(:)
        if (mstat.eq.1) then
           enginit = 1.5*t + eng
           conarr(mstat) = 1.0
        else
           conarr(mstat) = (1.5*t+eng)/enginit
        endif
    end subroutine status

    !> total md particle momentum
    subroutine momentum(m)
        real(kind=rk),intent(out) :: m(3)
        integer i,ii,ierror
        real(kind=rk) :: mm(3)

        m = 0.0_rk
        i = atompnt
        do ii = 1,nlocal
           m = m + particle_mass(P(i)) * P(i)%v
           i = list(i)
        enddo

        call mpi_allreduce(mm,m,3,MPI_REAL8,MPI_SUM,comm_cart,ierror)
    end subroutine momentum

    !> reduced temperature
    subroutine temperature(t)
        real(kind=rk),intent(out) :: t
        integer i,ii,ierror,n
        real(kind=rk) vv(3),tt

        t = 0.0
        i = atompnt
        do ii = 1,nlocal
           vv = P(i)%v
           t = t + dot_product(vv,vv)
           i = list(i)
        enddo
        call mpi_allreduce(t,tt,1,MPI_REAL8,MPI_SUM,comm_cart,ierror)
        call count_particles_all(n)
        t = tt / (3.0_rk*n)
    end subroutine temperature

    !> calculates the total translational kinetic energy of all particles
    subroutine rotational_kinetic_energy(er)
        real(kind=rk),intent(out) :: er
        real(kind=rk) :: tmp
        integer i,ii,ierror

        tmp = 0.0_rk
        i = atompnt
        do ii = 1,nlocal
           tmp = tmp + sum(particle_inertia(P(i)) * P(i)%w**2)
           i = list(i)
        enddo
        tmp = 0.5_rk * tmp
        call mpi_allreduce(tmp,er,1,MPI_REAL8,MPI_SUM,comm_cart,ierror)
    end subroutine rotational_kinetic_energy

    !> calculates the total translational kinetic energy of all particles
    subroutine translational_kinetic_energy(et)
        real(kind=rk),intent(out) :: et
        real(kind=rk) tmp
        integer i,ii,ierror

        tmp = 0.0_rk
        i = atompnt
        do ii = 1,nlocal
           tmp = tmp + particle_mass(P(i)) * dot_product(P(i)%v,P(i)%v)
           i = list(i)
        enddo
        tmp = 0.5_rk*tmp
        call mpi_allreduce(tmp,et,1,MPI_REAL8,MPI_SUM,comm_cart,ierror)
    end subroutine translational_kinetic_energy

    subroutine statistics
        integer k
        real(kind=rk) engcorr,prscorr

        call corrections(engcorr,prscorr)

        engarr(1:mstat) = engarr(1:mstat) + engcorr
        prsarr(1:mstat) = prsarr(1:mstat) + prscorr

        ! exclude timestep 0 from averaging
        if (mstat>1) then
           tmpave = sum(tmparr(2:mstat))/(mstat-1)
           e_trans_ave = sum(e_trans(2:mstat))/(mstat-1)
           e_rot_ave = sum(e_rot(2:mstat))/(mstat-1)
           engave = sum(engarr(2:mstat))/(mstat-1)
           rpotave = sum(rpotarr(2:mstat))/(mstat-1)
           prsave = sum(prsarr(2:mstat))/(mstat-1)
           do k=1,3
              momentumave(k) = sum(momentumarr(k,2:mstat))/(mstat-1)
           end do
           conave = sum(conarr(2:mstat))/(mstat-1)
        end if
    end subroutine statistics

    !> prepare list structure
    subroutine setup_list
        integer i

        freepnt = 1
        do i = 1,npmax
           list(i) = i + 1
        enddo
        list(npmax) = 0
        nlocal = 0
        atompnt = npmax + 1
    end subroutine setup_list

    !> sets  P(:)%v(3)  according to Lees-Edwards boundary conditions for all
    !> particles whose distance to the planes x=0.5 or x=tnx+0.5 is smaller than
    !>  boundary_width .
    subroutine fix_shear_boundary
        integer i,ii

        i = atompnt
        do ii = 1,nlocal
           if ((P(i)%x(1)-minx<boundary_width.and.P(i)%vnew(1)<0.0_rk).or.&
                &(P(i)%x(1)>maxx-boundary_width.and.P(i)%vnew(1)>0.0_rk)) then
              P(i)%vnew(1) = -P(i)%vnew(1)
              ! vold is lost, therefore correct v instead of recalculating it
              P(i)%v(1) = P(i)%v(1)+P(i)%vnew(1)
           end if
           i = list(i)
        end do
    end subroutine fix_shear_boundary

    !> issues an error if particles located within the boundary layer
    !> defined by \c boundary_width
    !>
    !> This is intended to check the initial configuration after
    !> particle placement.
    subroutine check_shear_boundary
        integer i,ii

        i = atompnt
        do ii = 1,nlocal
           if (P(i)%x(1)<minx+boundary_width.or.P(i)%x(1)>=maxx-boundary_width)&
                & call error_md('found particle in boundary layer---'&
                &//'disable enforced shear boundary layer (set '&
                &//'boundary_width=0.0) or change particle placement, '&
                &//'maybe with x0lo and x0hi!')

           i = list(i)
        end do
    end subroutine check_shear_boundary

    subroutine semiper(ident,direction,maxborder)
        integer ident,direction,counter
        integer not_direction(2)
        logical diffuse_direction,positive_v,maxborder
        real(kind=rk) ampl_v,ampl_v_old,k_b,r_v,ratio

#ifdef RWALK
        real(kind=rk) ampl_v_r_old, ampl_v_r
        logical positive_v_r
#endif

        k_b = 1.3806504E-23

        if(direction == 1)then
           not_direction(1)=2
           not_direction(2)=3
           diffuse_direction = diffuse_x
        endif
        if(direction == 2)then
           not_direction(1)=1
           not_direction(2)=3
           diffuse_direction = diffuse_y
        endif
        if(direction == 3)then
           not_direction(1)=1
           not_direction(2)=2
           diffuse_direction = diffuse_z
        endif

        if(maxborder)then
           P(ident)%x(direction) = maxpos(direction) - (P(ident)%x(direction)- maxpos(direction))
        else
           P(ident)%x(direction) = minpos(direction) + (minpos(direction)-P(ident)%x(direction))
        endif
        if (diffuse_direction) then
           ampl_v_old = P(ident)%v(not_direction(1))**2+P(ident)%v(not_direction(2))**2+ P(ident)%v(direction)**2
#ifdef RWALK
           ampl_v_r_old = P(ident)%v_r(not_direction(1))**2+P(ident)%v_r(not_direction(2))**2+ P(ident)%v_r(direction)**2
#endif

           if(P(ident)%v(direction).le.0.0)then
              positive_v=.true.
           else
              positive_v=.false.
           endif
#ifdef RWALK
           if(P(ident)%v_r(direction).le.0.0)then
              positive_v_r=.true.
           else
              positive_v_r=.false.
           endif
#endif
           call uni_dist(r_v)
           if(r_v==0.)print*,"rv returned 0 membrane on rank", myrankc
           P(ident)%vnew(not_direction(1)) = r_v*sqrt(k_b*temperat/molec_mass)*(delta_t/delta_x)
           call uni_dist(r_v)
           if(r_v==0.)print*,"rv returned 0 membrane on rank", myrankc
           P(ident)%vnew(not_direction(2)) = r_v*sqrt(k_b*temperat/molec_mass)*(delta_t/delta_x)
           call uni_dist(r_v)
           if(r_v==0.)print*,"rv returned 0 membrane on rank", myrankc
           P(ident)%vnew(direction) = r_v*sqrt(k_b*temperat/molec_mass)*(delta_t/delta_x)
           ampl_v = P(ident)%vnew(not_direction(1))**2+P(ident)%v(not_direction(2))**2+ P(ident)%v(direction)**2

           ratio = ampl_v_old/ampl_v
           do counter=1,3
              P(ident)%vnew(counter)= P(ident)%v(counter)*sqrt(ratio)
           end do

           if(positive_v)then
              P(ident)%vnew(direction) = abs(P(ident)%v(direction))
           else
              P(ident)%vnew(direction) = - abs(P(ident)%v(direction))
           endif
#ifdef RWALK
           call uni_dist(r_v)
           if(r_v==0.)print*,"rv returned 0 membrane on rank", myrankc
           P(ident)%v_r(not_direction(1)) = r_v*sqrt(k_b*temperat/molec_mass)*(delta_t/delta_x)
           call uni_dist(r_v)
           if(r_v==0.)print*,"rv returned 0 membrane on rank", myrankc
           P(ident)%v_r(not_direction(2)) = r_v*sqrt(k_b*temperat/molec_mass)*(delta_t/delta_x)
           call uni_dist(r_v)
           if(r_v==0.)print*,"rv returned 0 membrane on rank", myrankc
           P(ident)%v_r(direction) = r_v*sqrt(k_b*temperat/molec_mass)*(delta_t/delta_x)
           ampl_v_r = P(ident)%v_r(not_direction(1))**2+P(ident)%v_r(not_direction(2))**2+ P(ident)%v_r(direction)**2
           ratio = ampl_v_r_old/ampl_v_r
           do counter=1,3
              P(ident)%v_r(counter)= P(ident)%v_r(counter)*sqrt(ratio)
           end do
           if(positive_v_r)then
              P(ident)%v_r(direction) = abs(P(ident)%v_r(direction))
           else
              P(ident)%v_r(direction) = - abs(P(ident)%v_r(direction))
           endif
#endif

        else

#ifdef RWALK
           P(ident)%v_r(direction) = -P(ident)%v_r(direction)
#endif
           P(ident)%vnew(direction) = -P(ident)%vnew(direction)
        endif
    end subroutine semiper

#ifdef RWALK
    !> handle reflective rocks for a given particle
    !>
    !> \param[in] N local lattice chunk with halo of extent \c halo_extent
    !>
    !> \param[in,out] p the particle
    subroutine handle_reflective_rocks(N,p)
        type(lbe_site),intent(in) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        type(md_particle_type),intent(inout) :: p
        real(kind=rk) :: v_new(3),x_new(3)
        logical rock

        v_new(:) =0.0
        x_new(:) =0.0
        call rock_reflection(N,p%x,p%v_r,dt,v_new,x_new,rock)
        if(rock)then
           if (any(v_new==0.0_rk)) &
                &print*,"returned 0 vel at rock, on proc",myrankc
           p%x = x_new
           p%v_r = v_new
        endif
    end subroutine handle_reflective_rocks

    !> performs a pseudo-particle collision for a given particle
    !>
    !> \param[in] N local lattice chunk with halo of extent \c halo_extent
    !>
    !> \param[in,out] p the particle
    subroutine pseudo_particle_collision(N,p)
        type(lbe_site),intent(in) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        type(md_particle_type),intent(inout) :: p
	real(kind=rk),parameter :: k_b = 1.3806504E-23
	real(kind=rk) :: mu,p_abs,r_squaresum,r_v,v_abs
        real(kind=rk) :: lbm_cap(3),lbm_vel(3),p_cap(3),p_ort(3),sphere(3)
        real(kind=rk) :: c1(3),c2(3),cm(3),cr(3),cr_val,cr_new(3)
        real(kind=rk),dimension(2) :: r
        integer :: l
        logical :: not_yet,stuck

        p%vnew = 0.
        p%v = 0.
        lbm_vel(:) = 0.

        ! p%v=0
        ! set normal velocities zero if only r walk is used...

        call fluid_velocity_and_viscosity&
             &(N,(/1.0_rk,1.0_rk,1.0_rk/),p%x,lbm_vel,mu,stuck)
        if (stuck) then
           print '(A,I10,A,3ES10.2,A,3ES10.2,A)'&
                &,'WARNING: lost a particle in rock: uid='&
                &,p%uid,',x=(',p%x,'),v_r=(',p%v_r,')'
        else
!!$           lbm_vel = (/ 0.0_rk , 0.0_rk , 0.0_rk /)

           do l=1,3
              ! Get normal distributed number r_v by Box-Muller
              call uni_dist(r_v)
              ! Calculate first pre-collision velocity c1
              ! Maxwell velocity component representing the LBM fluid
              ! This is expectation value + normal-distributed value * standard deviation
              !         lbm_vel + r_v + sqrt( k_b * T  / m )
              ! Scaled by the lattice units
              !         v_SI = v_LBM / (dx/dt)
              c1(l) = lbm_vel(l)+&
                   & r_v*sqrt((k_b*temperat)/(molec_mass_lbm))&
                   &*(delta_t/delta_x)
              ! Second pre-collision velocity is contaminant particle velocity from MD
              c2(l) = p%v_r(l)
           end do

           ! Calculate velocity of center of mass (conserved)
           cm(:) = ((molec_mass_lbm * c1(:)) + (molec_mass * c2(:))) / (molec_mass_lbm + molec_mass)
           ! Calculate relative velocity of particle (absolute value is conserved s.b.)
           cr(:) = (c1(:) - c2(:))
           ! Calculate absolute value of relative velocity (conserved)
           cr_val = sqrt(cr(1)**2+cr(2)**2+cr(3)**2)
          
           ! For hard spheres the orientation of post-collision relative velocity is uniformly distributed
           ! Calculate uniformly distributed point on the unit-sphere by Marsaglia
           not_yet = .true.
           do
              call RANLUX(r,2)
              r(:) = 1.0_rk-2.0_rk*r(:)
              r_squaresum = r(1)**2 + r(2)**2
              if(r_squaresum .le. 1)then
                 sphere(1) = 2.0*r(1)*(sqrt(1.0-r_squaresum))
                 sphere(2) = 2.0*r(2)*(sqrt(1.0-r_squaresum))
                 sphere(3) = 1.0-2.*r_squaresum
                 not_yet = .false.
              endif
              if(not_yet .eqv. .false.)exit
           enddo

           ! scale random unit-vector by conserved absolute relative velocity
           cr_new(:) = cr_val * sphere(:)

           ! Calculate post-collision velocity of the contaminant particle
           ! LBM representing particle is 'forgotten'
           p%v_r = cm(:) - (molec_mass_lbm / (molec_mass_lbm + molec_mass) ) * cr_new(:)

           p%sdx=0.0
        endif
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!!! Random Speed for Pseudo Particle, First option creates a randomly
!!$!!! directed vector of unit length and then scales it Second option
!!$!!! directly creates 3 gauss distributed velocities in the three
!!$!!! dimensions
!!$!!!
!!$!!! Both are equal but Mr. Hilfer thinks not, so the first is used
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!!! First option
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$           not_yet = .true.
!!$           do
!!$              call RANLUX(r,2)
!!$              r(:) = 1.0_rk-2.0_rk*r(:)
!!$              r_squaresum = r(1)**2 + r(2)**2
!!$              if(r_squaresum .le. 1)then
!!$                 sphere(1) = 2.0*r(1)*(sqrt(1.0-r_squaresum))
!!$                 sphere(2) = 2.0*r(2)*(sqrt(1.0-r_squaresum))
!!$                 sphere(3) = 1.0-2.*r_squaresum
!!$                 not_yet = .false.
!!$              endif
!!$              if(not_yet == .false.)exit
!!$           enddo
!!$           v_abs= 0.0
!!$
!!$           do l=1,3
!!$              call uni_dist(r_v)
!!$              v_abs = v_abs + (r_v*sqrt((k_b*temperat)/(molec_mass_lbm))&
!!$                   &*(delta_t/delta_x))**2
!!$           end do
!!$           v_abs = sqrt(v_abs)
!!$           lbm_vel(:) = lbm_vel(:) + sphere(:)*v_abs
!!$#ifdef RWALKTEST
!!$          p%v_r=lbm_vel
!!$#else
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!!! Second option
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           do l=1,3
!              call uni_dist(r_v)
!              lbm_vel(l) = lbm_vel(l)+&
!                   & r_v*sqrt((k_b*temperat)/(molec_mass_lbm))&
!                   &*(delta_t/delta_x)
!           end do
!!$
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!!!
!!$!!! Splitting the velocity of the MD - particle in an orthogonal part
!!$!!! and a part along the velocity axis of the pseudo particle
!!$!!!
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$           lbm_cap = lbm_vel/sqrt(dot_product(lbm_vel,lbm_vel))
!!$           p_abs = dot_product(lbm_cap,p%v_r)
!!$           p_cap(:) = p_abs*lbm_cap(:)
!!$           p_ort = p%v_r-p_cap
!!$
!!$           do l=1,3
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!!! Perform the collision and add the old orthogonal part
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$              p_cap(l)=((molec_mass-&
!!$                   &molec_mass_lbm)*p_cap(l)+2*molec_mass_lbm*lbm_vel(l))&
!!$                   &/(molec_mass+molec_mass_lbm)
!!$              p%v_r(l) = p_cap(l) + p_ort(l)
!!$           end do
!!$#endif
!!$           p%sdx=0.0
!!$        endif
    end subroutine pseudo_particle_collision
#endif

!!$    ! random number
!!$    function rdm(msd)
!!$        integer msd
!!$        real rdm
!!$        INTEGER, DIMENSION(1) :: OLD, SEED ! THIS PROGRAM ASSUMES K = 1
!!$	INTEGER :: I, K
!!$        REAL, DIMENSION(1) :: HARVEST
!!$        SEED(1) = msd
!!$	CALL RANDOM_SEED
!!$        CALL RANDOM_SEED(SIZE=K)
!!$	CALL RANDOM_SEED(GET=OLD(1:K))
!!$	CALL RANDOM_NUMBER(HARVEST)
!!$	CALL RANDOM_SEED(GET=OLD(1:K))
!!$	CALL RANDOM_SEED(PUT=SEED(1:K))
!!$	CALL RANDOM_SEED(GET=OLD(1:K))
!!$	CALL RANDOM_NUMBER(HARVEST)
!!$        rdm = HARVEST(1)
!!$    end function rdm

    !> allocte memory needed for random_reentering, but only on the
    !> first call
    subroutine setup_random_reentering_on_demand
!!$        integer stat
!!$
!!$        if (.not.allocated(reenter_list)) then
!!$           allocate(reenter_list(natoms),stat=stat)
!!$           call check_allocate(stat&
!!$                &,'setup_random_reentering_on_demand(): reenter_list')
!!$        end if
!!$
!!$        if (.not.allocated(reenter_pos_list)) then
!!$           allocate(reenter_pos_list(natoms),stat=stat)
!!$           call check_allocate(stat&
!!$                &,'setup_random_reentering_on_demand(): reenter_pos_list')
!!$        end if
    end subroutine setup_random_reentering_on_demand

    !> perform random reentering
    !>
    !> Done for particles that were marked to be repositioned, the
    !> communication steps at the end of this subroutine replace the
    !> usual exchange/borders/neighbor sequence.
    subroutine do_random_reentering(N)
        type(lbe_site),intent(inout) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        integer i,ierror,x,y,z

        DEBUG_MPI_MSG("Entered do_random_reentering")

!!! reenable this one time... !!!
        call error_md('do_random_reentering() currently not working in '&
             &//'version6. The reason is mostly that particle uids are not '&
             &//'contiguous anymore.')
!!$        ! all procesors
!!$        do pp=0,nprocs-1
!!$           if (pp==myrankc) then
!!$
!!$              ! sending the particle leaving the domain in this processor
!!$              call mpi_bcast(n_reenter,1,MPI_INTEGER,pp,comm_cart,ierror)
!!$
!!$              do i=1,n_reenter
!!$                 call mpi_bcast(P(reenter_list(i))%uid,1,MPI_INTEGER,pp,comm_cart,ierror)
!!$                 call mpi_bcast(P(reenter_list(i))%v(:),3,MPI_REAL8,pp,comm_cart,ierror)
!!$                 do x=0,nx+1
!!$                    do y=0,ny+1
!!$                       do z=0,nz+1
!!$                          if (N(x,y,z)%rock_state==real(P(reenter_list(i))%uid,kind=rk)) then
!!$                             N(x,y,z)%rock_state=0
!!$                             call boltz_dist(P(reenter_list(i))%v(1), &
!!$                                  & P(reenter_list(i))%v(2), &
!!$                                  & P(reenter_list(i))%v(3), &
!!$                                  0.0_rk,0.0_rk,0.0_rk, &
!!$                                  0.0_rk,0.0_rk,0.0_rk,N(x,y,z)%n_r)
!!$
!!$                             N(x,y,z)%n_r=N(x,y,z)%n_r*pr
!!$                          end if
!!$                       end do
!!$                    end do
!!$                 end do
!!$
!!$                 check_rdm=.false.
!!$                 do while (.not.check_rdm) 
!!$                    ! calculating new position for this particles
!!$                    if (reenter_pos_list(i)==1) then
!!$                       call  RANLUX(rntr_rdm,1)
!!$                       P(reenter_list(i))%x(2) = tsize(2)*rntr_rdm(1) + 0.5
!!$                       call  RANLUX(rntr_rdm,1)
!!$                       P(reenter_list(i))%x(3) = tsize(3)*rntr_rdm(1) + 0.5
!!$                    endif
!!$
!!$                    if (reenter_pos_list(i)==2) then
!!$                       call  RANLUX(rntr_rdm,1)
!!$                       P(reenter_list(i))%x(1) = tsize(1)*rntr_rdm(1) + 0.5
!!$                       call  RANLUX(rntr_rdm,1)
!!$                       P(reenter_list(i))%x(3) = tsize(3)*rntr_rdm(1) + 0.5
!!$                    endif
!!$
!!$                    if (reenter_pos_list(i)==3) then
!!$                       call  RANLUX(rntr_rdm,1)
!!$                       P(reenter_list(i))%x(1) = tsize(1)*rntr_rdm(1) + 0.5
!!$                       call  RANLUX(rntr_rdm,1)
!!$                       P(reenter_list(i))%x(2) = tsize(2)*rntr_rdm(1) + 0.5
!!$                    end if
!!$
!!$                    ! checking in own procesor
!!$                    call checking_rdm(P(reenter_list(i))%x(:), &
!!$                         & P(reenter_list(i))%o, P(reenter_list(i))%uid, N, check_rdm)
!!$
!!$                    ! sending the positions and orientation to check out
!!$                    call mpi_bcast(P(reenter_list(i))%x(:),3,MPI_REAL8,pp,comm_cart,ierror)
!!$                    call mpi_bcast(P(reenter_list(i))%o(:),3,MPI_REAL8,pp,comm_cart,ierror)
!!$                    call mpi_bcast(P(reenter_list(i))%uid,1,MPI_INTEGER,pp,comm_cart,ierror)
!!$
!!$                    ! getting the answer of the check out
!!$                    call mpi_allreduce(check_rdm,check2_rdm,1,MPI_LOGICAL,MPI_LAND,comm_cart,ierror)
!!$                    check_rdm = check2_rdm
!!$                 end do
!!$              end do
!!$           else
!!$              ! getting the information of the particle leaving in the
!!$              ! other processors
!!$              call mpi_bcast(g_reenter,1,MPI_INTEGER,pp,comm_cart,ierror)
!!$
!!$              do i=1,g_reenter
!!$                 call mpi_bcast(guid,1,MPI_INTEGER,pp,comm_cart,ierror)
!!$                 call mpi_bcast(gv,3,MPI_REAL8,pp,comm_cart,ierror)
!!$                 do x=0,nx+1
!!$                    do y=0,ny+1
!!$                       do z=0,nz+1
!!$                          if (N(x,y,z)%rock_state==real(guid,kind=rk)) then
!!$                             N(x,y,z)%rock_state=0
!!$                             call boltz_dist(gv(1),gv(2),gv(3),0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,N(x,y,z)%n_r)
!!$                             N(x,y,z)%n_r=N(x,y,z)%n_r*pr
!!$                          end if
!!$                       end do
!!$                    end do
!!$                 end do
!!$
!!$                 check_rdm=.false.
!!$
!!$                 do while (.not.check_rdm)
!!$                    ! getting the positions and orientation to check
!!$                    ! out
!!$                    call mpi_bcast(gp(:),3,MPI_REAL8,pp,comm_cart,ierror)
!!$                    call mpi_bcast(go(:),3,MPI_REAL8,pp,comm_cart,ierror)
!!$                    call mpi_bcast(guid,1,MPI_INTEGER,pp,comm_cart,ierror)
!!$
!!$                    ! checking
!!$                    call checking_rdm(gp(:),go(:),guid,N,check_rdm)
!!$
!!$                    ! sending answer
!!$                    call mpi_allreduce(check_rdm,check2_rdm,1,MPI_LOGICAL,MPI_LAND,comm_cart,ierror)
!!$                    check_rdm = check2_rdm
!!$                 end do
!!$              end do
!!$           end if
!!$        end do
!!$
!!$        if (myrankc==0) then
!!$           ! allocate buffers only once on first invocation
!!$           if (.not.allocated(tp)) then
!!$              allocate (tp(natoms),stat=stat)
!!$              call check_allocate(stat,'tp')
!!$           end if
!!$
!!$           ! allow  natoms  to alter during the simulation
!!$           if (natoms>size(tp)) then
!!$              call log_msg_md(&
!!$                   &'random reentering buffers are too small,'&
!!$                   &//' allocating more memory...')
!!$              deallocate (tp)
!!$              allocate (tp(2*natoms),stat=stat)
!!$              call check_allocate(stat,'tp')
!!$           end if
!!$        end if
!!$
!!$        ! taking care, thid step could be really slow if many
!!$        ! processor are involve because all the particle are
!!$        ! redistributed in all processor
!!$        if (provide_uid2i) call clear_uid2i
!!$        call start_timer(ti_md_comm)
!!$        call gather_particles(tp)
!!$        call setup_list
!!$        call scatter_particles(tp(:))
!!$        call stop_timer(ti_md_comm)
!!$        ! when at least one particle trajectory was reset then update
!!$        ! neighbor lists (the equivalent of exchange() is already done
!!$        ! in scatter_particles() )
!!$        if (communicate_rotations_s) call calculate_rotations_s
!!$        call start_timer(ti_md_comm)
!!$        call borders
!!$        call stop_timer(ti_md_comm)
!!$        if (provide_uid2i) call create_uid2i
!!$        call start_timer(ti_md_neigh)
!!$        call neighbor
!!$        call stop_timer(ti_md_neigh)

          DEBUG_MPI_MSG("Returning from do_random_reentering")
    end subroutine do_random_reentering

    subroutine checking_rdm(px,po,puid,Nr,check_rdm)
       real(kind=rk) :: px(3), po(3), local(3),dist(9),dist2(27),gdist
       integer,intent(in) :: puid
       logical :: check_rdm
       integer :: x,y,z,ii,i,j
       type(lbe_site),intent(in) :: &
            & Nr(1-halo_extent:,1-halo_extent:,1-halo_extent:)
       type(md_particle_type) :: dummy

       if (interaction=='tracer') then
          call local_coordinates(px(:),local)
          if (local(1)+0.5>=1.and.local(1)-0.5<=nx .and. &
               & local(2)+0.5>=1.and.local(2)-0.5<=ny .and. &
               & local(3)+0.5>=1.and.local(3)-0.5<=nz ) then
             if (Nr(int(local(1)+0.5),int(local(2)+0.5),int(local(3)+0.5))&
                  &%rock_state/=0.0) then
                check_rdm = .false.
             else
                check_rdm = .true.
             end if
          else
             check_rdm = .true.
          end if
       else if (interaction=='friction') then
          call local_coordinates(px(:),local)
          if(local(1)+0.5>=1.and.local(1)-0.5<=nx .and. &
               & local(2)+0.5>=1.and.local(2)-0.5<=ny .and. &
               & local(3)+0.5>=1.and.local(3)-0.5<=nz) then
             if (Nr(int(local(1)+0.5),int(local(2)+0.5),int(local(3)+0.5))&
                  &%rock_state/=0.0) then
                check_rdm = .false.
             else
                check_rdm = .true.
             end if
          else
             check_rdm = .true.
          end if
       else if (interaction=='ladd') then
          call local_coordinates(px(:),local)
          check_rdm = .true.
          do x=1,nx
             do y=1,ny
                do z=1,nz
                   ! dummy particle to pass to in_particle(). Only
                   ! it's orientation will be read, polydispersity
                   ! does not support random reentering, anyway
                   dummy%o = po

                   ! checking particle-solid overlaping
                   if (in_particle(dummy,local,real((/x,y,z/),kind=rk)).and. &
                        & Nr(x,y,z)%rock_state/=0.0) then
                      check_rdm = .false.
                   end if

                end do
             end do
          end do


          ! checking particle-particle overlaping
          i = atompnt
          do ii=1,nlocal
             !       do_something_with_atom_i(P(i)%x, P(i)%o)
             if(puid/=P(i)%uid) then
                dist2(1)  = (px(1)-P(i)%x(1))**2+(px(2)-P(i)%x(2))**2+(px(3)-P(i)%x(3))**2
                dist2(2)  = (px(1)-P(i)%x(1))**2+(px(2)-P(i)%x(2))**2+(px(3)-P(i)%x(3)+tsize(3))**2
                dist2(3)  = (px(1)-P(i)%x(1))**2+(px(2)-P(i)%x(2))**2+(px(3)-P(i)%x(3)-tsize(3))**2
                dist2(4)  = (px(1)-P(i)%x(1))**2+(px(2)-P(i)%x(2)+tsize(2))**2+(px(3)-P(i)%x(3))**2
                dist2(5)  = (px(1)-P(i)%x(1))**2+(px(2)-P(i)%x(2)+tsize(2))**2+(px(3)-P(i)%x(3)+tsize(3))**2
                dist2(6)  = (px(1)-P(i)%x(1))**2+(px(2)-P(i)%x(2)+tsize(2))**2+(px(3)-P(i)%x(3)-tsize(3))**2
                dist2(7)  = (px(1)-P(i)%x(1))**2+(px(2)-P(i)%x(2)-tsize(2))**2+(px(3)-P(i)%x(3))**2
                dist2(8)  = (px(1)-P(i)%x(1))**2+(px(2)-P(i)%x(2)-tsize(2))**2+(px(3)-P(i)%x(3)+tsize(3))**2
                dist2(9)  = (px(1)-P(i)%x(1))**2+(px(2)-P(i)%x(2)-tsize(2))**2+(px(3)-P(i)%x(3)-tsize(3))**2
                dist2(10) = (px(1)-P(i)%x(1)+tsize(1))**2+(px(2)-P(i)%x(2))**2+(px(3)-P(i)%x(3))**2
                dist2(11) = (px(1)-P(i)%x(1)+tsize(1))**2+(px(2)-P(i)%x(2))**2+(px(3)-P(i)%x(3)+tsize(3))**2
                dist2(12) = (px(1)-P(i)%x(1)+tsize(1))**2+(px(2)-P(i)%x(2))**2+(px(3)-P(i)%x(3)-tsize(3))**2
                dist2(13) = (px(1)-P(i)%x(1)+tsize(1))**2+(px(2)-P(i)%x(2)+tsize(2))**2+(px(3)-P(i)%x(3))**2
                dist2(14) = (px(1)-P(i)%x(1)+tsize(1))**2+(px(2)-P(i)%x(2)+tsize(2))**2+(px(3)-P(i)%x(3)+tsize(3))**2
                dist2(15) = (px(1)-P(i)%x(1)+tsize(1))**2+(px(2)-P(i)%x(2)+tsize(2))**2+(px(3)-P(i)%x(3)-tsize(3))**2
                dist2(16) = (px(1)-P(i)%x(1)+tsize(1))**2+(px(2)-P(i)%x(2)-tsize(2))**2+(px(3)-P(i)%x(3))**2
                dist2(17) = (px(1)-P(i)%x(1)+tsize(1))**2+(px(2)-P(i)%x(2)-tsize(2))**2+(px(3)-P(i)%x(3)+tsize(3))**2
                dist2(18) = (px(1)-P(i)%x(1)+tsize(1))**2+(px(2)-P(i)%x(2)-tsize(2))**2+(px(3)-P(i)%x(3)-tsize(3))**2
                dist2(19) = (px(1)-P(i)%x(1)-tsize(1))**2+(px(2)-P(i)%x(2))**2+(px(3)-P(i)%x(3))**2
                dist2(20) = (px(1)-P(i)%x(1)-tsize(1))**2+(px(2)-P(i)%x(2))**2+(px(3)-P(i)%x(3)+tsize(3))**2
                dist2(21) = (px(1)-P(i)%x(1)-tsize(1))**2+(px(2)-P(i)%x(2))**2+(px(3)-P(i)%x(3)-tsize(3))**2
                dist2(22) = (px(1)-P(i)%x(1)-tsize(1))**2+(px(2)-P(i)%x(2)+tsize(2))**2+(px(3)-P(i)%x(3))**2
                dist2(23) = (px(1)-P(i)%x(1)-tsize(1))**2+(px(2)-P(i)%x(2)+tsize(2))**2+(px(3)-P(i)%x(3)+tsize(3))**2
                dist2(24) = (px(1)-P(i)%x(1)-tsize(1))**2+(px(2)-P(i)%x(2)+tsize(2))**2+(px(3)-P(i)%x(3)-tsize(3))**2
                dist2(25) = (px(1)-P(i)%x(1)-tsize(1))**2+(px(2)-P(i)%x(2)-tsize(2))**2+(px(3)-P(i)%x(3))**2
                dist2(26) = (px(1)-P(i)%x(1)-tsize(1))**2+(px(2)-P(i)%x(2)-tsize(2))**2+(px(3)-P(i)%x(3)+tsize(3))**2
                dist2(27) = (px(1)-P(i)%x(1)-tsize(1))**2+(px(2)-P(i)%x(2)-tsize(2))**2+(px(3)-P(i)%x(3)-tsize(3))**2

                gdist=minval(dist2(:))
                if (gdist.lt.(2.0_rk*max(r_orth,r_para))**2) then
                   check_rdm=.false.
                end if
             end if

             i=list(i)
          end do


! checking moving particle-particle overlaping
          do j=1,n_reenter
             if (puid/=P(reenter_list(j))%uid) then
                if (reenter_pos_list(i)==1) then
                   dist(1) = (px(2)-P(reenter_list(j))%x(2)-tsize(2))**2+(px(3)-P(reenter_list(j))%x(3))**2
                   dist(2) = (px(2)-P(reenter_list(j))%x(2)+tsize(2))**2+(px(3)-P(reenter_list(j))%x(3))**2
                   dist(3) = (px(2)-P(reenter_list(j))%x(2))**2+(px(3)-P(reenter_list(j))%x(3)-tsize(3))**2
                   dist(4) = (px(2)-P(reenter_list(j))%x(2))**2+(px(3)-P(reenter_list(j))%x(3)+tsize(3))**2
                   dist(5) = (px(2)-P(reenter_list(j))%x(2)-tsize(2))**2+(px(3)-P(reenter_list(j))%x(3)-tsize(3))**2
                   dist(6) = (px(2)-P(reenter_list(j))%x(2)-tsize(2))**2+(px(3)-P(reenter_list(j))%x(3)+tsize(3))**2
                   dist(7) = (px(2)-P(reenter_list(j))%x(2)+tsize(2))**2+(px(3)-P(reenter_list(j))%x(3)-tsize(3))**2
                   dist(8) = (px(2)-P(reenter_list(j))%x(2)+tsize(2))**2+(px(3)-P(reenter_list(j))%x(3)+tsize(3))**2
                   dist(9) = (px(2)-P(reenter_list(j))%x(2))**2+(px(3)-P(reenter_list(j))%x(3))**2

                   gdist=minval(dist(:))

                   if (gdist.lt.(2.0_rk*max(r_orth,r_para))**2) then
                      check_rdm=.false.
                   end if
                end if

                if (reenter_pos_list(i)==2) then
                   dist(1) = (px(1)-P(reenter_list(j))%x(1)-tsize(1))**2+(px(3)-P(reenter_list(j))%x(3))**2
                   dist(2) = (px(1)-P(reenter_list(j))%x(1)+tsize(1))**2+(px(3)-P(reenter_list(j))%x(3))**2
                   dist(3) = (px(1)-P(reenter_list(j))%x(1))**2+(px(3)-P(reenter_list(j))%x(3)-tsize(3))**2
                   dist(4) = (px(1)-P(reenter_list(j))%x(1))**2+(px(3)-P(reenter_list(j))%x(3)+tsize(3))**2
                   dist(5) = (px(1)-P(reenter_list(j))%x(1)-tsize(1))**2+(px(3)-P(reenter_list(j))%x(3)-tsize(3))**2
                   dist(6) = (px(1)-P(reenter_list(j))%x(1)-tsize(1))**2+(px(3)-P(reenter_list(j))%x(3)+tsize(3))**2
                   dist(7) = (px(1)-P(reenter_list(j))%x(1)+tsize(1))**2+(px(3)-P(reenter_list(j))%x(3)-tsize(3))**2
                   dist(8) = (px(1)-P(reenter_list(j))%x(1)+tsize(1))**2+(px(3)-P(reenter_list(j))%x(3)+tsize(3))**2
                   dist(9) = (px(1)-P(reenter_list(j))%x(1))**2+(px(3)-P(reenter_list(j))%x(3))**2
                   gdist=minval(dist(:))

                   if (gdist.lt.(2.0_rk*max(r_orth,r_para))**2) then
                      check_rdm=.false.
                   end if
                end if

                if (reenter_pos_list(i)==3) then
                   dist(1) = (px(2)-P(reenter_list(j))%x(2)-tsize(2))**2+(px(1)-P(reenter_list(j))%x(1))**2
                   dist(2) = (px(2)-P(reenter_list(j))%x(2)+tsize(2))**2+(px(1)-P(reenter_list(j))%x(1))**2
                   dist(3) = (px(2)-P(reenter_list(j))%x(2))**2+(px(1)-P(reenter_list(j))%x(1)-tsize(1))**2
                   dist(4) = (px(2)-P(reenter_list(j))%x(2))**2+(px(1)-P(reenter_list(j))%x(1)+tsize(1))**2
                   dist(5) = (px(2)-P(reenter_list(j))%x(2)-tsize(2))**2+(px(1)-P(reenter_list(j))%x(1)-tsize(1))**2
                   dist(6) = (px(2)-P(reenter_list(j))%x(2)-tsize(2))**2+(px(1)-P(reenter_list(j))%x(1)+tsize(1))**2
                   dist(7) = (px(2)-P(reenter_list(j))%x(2)+tsize(2))**2+(px(1)-P(reenter_list(j))%x(1)-tsize(1))**2
                   dist(8) = (px(2)-P(reenter_list(j))%x(2)+tsize(2))**2+(px(1)-P(reenter_list(j))%x(1)+tsize(1))**2
                   dist(9) = (px(2)-P(reenter_list(j))%x(2))**2+(px(1)-P(reenter_list(j))%x(1))**2

                   gdist=minval(dist(:))

                   if (gdist.lt.(2.0_rk*max(r_orth,r_para))**2) then
                      check_rdm=.false.
                   end if
                end if
             end if
          end do

          else if (interaction=='none') then
          call local_coordinates(px(:),local)
          if(local(1)+0.5>=1.and.local(1)-0.5<=nx.and.local(2)+0.5>=1 .and. &
               & local(2)-0.5<=ny.and.local(3)+0.5>=1 .and. &
               & local(3)-0.5<=nz) then
             if (Nr(int(local(1)+0.5),int(local(2)+0.5),int(local(3)+0.5))&
                  &%rock_state/=0.0) then
                check_rdm = .false.
             else
                check_rdm = .true.
             end if
           else
              check_rdm = .true.
           end if
        endif
      end subroutine checking_rdm

!     subroutine checking_rdm(px,po,Nr,check_rdm)
!     use lbe_md_globals_module
!	real(kind=rk) :: px(3), po(3), local(3)
!        logical :: check_rdm
!        integer :: x,y,z,ii,i
!        type(lbe_site),intent(in) :: &
!             &Nr(1-halo_extent:,1-halo_extent:,1-halo_extent:)

!        if (interaction=='tracer') then                    
!           call local_coordinates(px(:),local)
!           if (local(1)+0.5>=1.and.local(1)-0.5<=nx.and.local(2)+0.5>=1.and.local(2)-0.5<=ny.and.local(3)+0.5>=1.and.local(3)-0.5<=nz) then
!              if (Nr(int(local(1)+0.5),int(local(2)+0.5),int(local(3)+0.5))%rock_state/=0.0) then
!                 check_rdm = .false.
!              else
!                 check_rdm = .true.
!              end if
!           else
!              check_rdm = .true.
!           end if
!        else if (interaction=='friction') then  
!           call local_coordinates(px(:),local)
!           if(local(1)+0.5>=1.and.local(1)-0.5<=nx.and.local(2)+0.5>=1.and.local(2)-0.5<=ny.and.local(3)+0.5>=1.and.local(3)-0.5<=nz) then
!              if (Nr(int(local(1)+0.5),int(local(2)+0.5),int(local(3)+0.5))%rock_state/=0.0) then
!                 check_rdm = .false.
!              else
!                 check_rdm = .true.
!              end if
!           else
!              check_rdm = .true.
!           end if
!        else if (interaction=='ladd') then
!           call local_coordinates(px(:),local)
!           check_rdm = .true.
!           do x=1,nx       
!           do y=1,ny       
!           do z=1,nz       

      ! checking particle-solid overlaping
!              if (in_particle(po,px,real((/x,y,z/),kind=rk)).and.&
!                     &Nr(x,y,z)%rock_state/=0.0) then
!                 check_rdm = .false.             
!              end if

      ! checking particle-particle overlaping   

!              i = atompnt
!              do ii=1,nlocal
       !       do_something_with_atom_i(P(i)%x, P(i)%o)


!                 if( po(1)/=P(i)%o(1) .and. po(2)/=P(i)%o(2) .and. po(3)/=P(i)%o(3) &
!                      & .and. px(1)/=P(i)%x(1).and. po(2)/=P(i)%x(2).and. po(3)/=P(i)%x(3)) then              
!                   if (in_particle(po,px,real((/x,y,z/),kind=rk)).and.&
!                        & in_particle(P(i),P(i)%x,real((/x,y,z/),kind=rk))) then
!                      check_rdm = .false.             
!                    end if
!                 end if
!                 i=list(i)       
!              end do
              
!           end do
!           end do
!           end do

!        else if (interaction=='none') then
!           call local_coordinates(px(:),local)
!           if(local(1)+0.5>=1.and.local(1)-0.5<=nx.and.local(2)+0.5>=1.and.local(2)-0.5<=ny.and.local(3)+0.5>=1.and.local(3)-0.5<=nz) then
!              if (Nr(int(local(1)+0.5),int(local(2)+0.5),int(local(3)+0.5))%rock_state/=0.0) then
!                 check_rdm = .false.
!              else
!                 check_rdm = .true.
!              end if
!           else
!!              check_rdm = .true.
!           end if
!        endif
!    end subroutine checking_rdm

      !> remove all entries from uid2i
      subroutine clear_uid2i
          call Mii_clear(uid2i)
      end subroutine clear_uid2i

      !> fill owned and halo'ed particles into uid2i
      subroutine create_uid2i
          integer i,ii

          call check_allocate(Mii_provide_capacity(uid2i,nlocal+nother),'uid2i')

          call Mii_rewind(uid2i)

          i = atompnt
          fill_particles_own: do ii = 1,nlocal
             call Mii_preinsert(uid2i,P(i)%uid,i)
             i = list(i)
          enddo fill_particles_own

          fill_particles_haloed: do ii = 1,nother
             call Mii_preinsert(uid2i,P(i)%uid,i)
             i = i+1
          enddo fill_particles_haloed

          call Mii_commit(uid2i)
      end subroutine create_uid2i
#endif
end module lbe_md_module
