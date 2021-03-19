#include "lbe.h"

!> per-particle pair lubrication correction for (in principle)
!> arbitrary convex bodies inspired by Cox, 1974 and Claeys&Brady,
!> 1989.
module lbe_md_fluid_ladd_club_module
#ifdef MD
    use lbe_globals_module, only: esmall,myrankc,pi,rk
    use lbe_helper_module, only: cross_product,norm,unit_vector
    use lbe_log_module
    use lbe_md_fluid_ladd_mc_module, only: pfr
    use lbe_md_fluid_ladd_parms_module, only: delta_c,delta_c_R,delta_c_T&
         &,inv_delta_c,log_delta_c_T,log_delta_c_R,lubrication_clip_dist&
         &,lubrication_stiffness,lubrication_tangential
    use lbe_md_globals_module
    use lbe_md_helper_module, only: error_md,log_msg_md,particle_radii&
         &,space_to_body_matrix
    use lbe_parallel_module, only: check_allocate,comm_cart
    use lbe_parms_module, only: nt,tau_r
    use map_module, only: Mi6r_clear,Mi6r_count,Mi6r_exist,Mi6r_init,Mi6r_map&
         &,Mi6r_provide_capacity,Mi6r_safe_assign,Mi6r_type

    implicit none
    include 'mpif.h'
    private

    public club_force_and_torque,club_set_growth_factor,club_setup&
         &,club_summary_fluid

    !> enable per-particle pair lubrication correction for (in
    !> principle) arbitrary convex bodies inspired by Cox, 1974 and
    !> Claeys&Brady, 1989
    logical,save,public :: lubrication_cox=.false.

    !> don't store more than this many pairs of surface points. If new
    !> pairs should be added, first clear all existing pairs. This
    !> way, avoid accumulating unused pairs in \c surface_points and
    !> spending too much time on the binary search for looking them
    !> up.
    integer,save,public :: lubrication_cox_max_store_points=100

    !> maximum upper gap cutoff
    real(kind=rk),save :: delta_c_max

    !> number of times that the closest distance between two particles
    !> needed to be computed
    real(kind=rk),save :: closest_surface_points_count=0.0_rk

    !> total number of iterations (including the unsuccessful ones)
    !> performed to find the closest distance between particles
    !> (double precision real number should prevent integer overflow)
    real(kind=rk),save :: closest_surface_points_iterations=0.0_rk

    !> holds surface points of minimum distance for particle pairs
    type(Mi6r_type),save :: surface_points

    !> local minimum number of time steps between two resets of \c
    !> surface_points
    integer,save :: min_surface_points_lifetime=&
         &huge(min_surface_points_lifetime)

    !> last time step where the local \c surface_points was cleaned
    integer,save :: surface_points_reset_time

    !> maximum expected particle uid. Set at the begin of a
    !> simulation. This assumption is made to ease the implementation
    !> in \c pair_uid(). \warning It is incompatible with inserting
    !> particles as in \c lbe_md_dynamic_module !!!
    integer,save :: max_particle_uid

    !> prefactor to all particle dimensions, used during \c growing_stage
    real(kind=rk),save :: growth_factor=1.0_rk

contains

    !> apply a force and a torque acting on a specific point on a particle
    !>
    !> \param[in] k index in \c P of particle
    !>
    !> \param[in] x point where force attacks in global space-fixed coordinates
    !>
    !> \param[in] f force
    !>
    !> \param[in] t torque
    subroutine apply_force_and_torque(k,x,f,t)
        integer,intent(in) :: k
        real(kind=rk),intent(in) :: x(3),f(3),t(3)

        P(k)%f = P(k)%f + f
        P(k)%t = P(k)%t + (t + cross_product(x-P(k)%x,f))
    end subroutine apply_force_and_torque

    !> try different strategies to find the two points on the surfaces
    !> of two particles that minimize the distance between
    !>
    !> \param[in] i first particle index in \c P
    !>
    !> \param[in] j second particle index in \c P
    !>
    !> \param[in] rij minimum image vector from particle \c j to \c i
    !>
    !> \param[in] ropi half axes \c R_orth and \c R_para for particle \c i
    !>
    !> \param[in] ropj half axes \c R_orth and \c R_para for particle \c j
    !>
    !> \param[out] xi point on particle \c i
    !>
    !> \param[out] xj point on particle \c j
    !>
    !> \param[out] g unit gap vector from \c i to \c j
    !>
    !> \param[out] gap gap width
    subroutine compute_closest_surface_points(i,j,rij,ropi,ropj,xi,xj,g,gap)
        integer,intent(in) :: i,j
        real(kind=rk),intent(in) :: rij(3),ropi(2),ropj(2)
        real(kind=rk),intent(out) :: xi(3),xj(3),g(3),gap
        logical :: restart,success
        integer :: n,n_tot
        real(kind=rk) :: rrxi(3),rrxj(3),rxi(3),rxj(3)

        n_tot = 0

        ! start iteration from previous positions, if available
        call retrieve_surface_points(i,j,rrxi,rrxj,restart)

        call closest_surface_points(i,j,rij,ropi,ropj,10,.true.,restart&
             &,rrxi,rrxj,rxi,rxj,g,gap,success,n)
        n_tot = n_tot + n
        if (.not.success) then
           call closest_surface_points(i,j,rij,ropi,ropj,1000,.false.,restart&
                &,rrxi,rrxj,rxi,rxj,g,gap,success,n)
           n_tot = n_tot + n
           if (restart.and..not.success) then
              ! last try: constant (small) step size, start from centers
              call closest_surface_points(i,j,rij,ropi,ropj&
                   &,1000,.false.,.false.,rrxi,rrxj,rxi,rxj,g,gap,success,n)
              n_tot = n_tot + n
           end if
        end if

        if (.not.success) then
           write (msgstr,'("compute_closest_surface_points(): '&
                &//'no convergence achieved. '&
                &//'rij=(",2(ES15.8,X),ES15.8,"), '&
                &//'oi=(",2(ES15.8,X),ES15.8,"), '&
                &//'oj=(",2(ES15.8,X),ES15.8,"), '&
                &//'qi=(",3(ES15.8,X),ES15.8,"), '&
                &//'qj=(",3(ES15.8,X),ES15.8,"), '&
                &//'uidi=",I0,", uidj=",I0," restart=",L1,", '&
                &//'rrxi=(",2(ES15.8,X),ES15.8,"), '&
                &//'rrxj=(",2(ES15.8,X),ES15.8,")")') &
                &rij,P(i)%o,P(j)%o,P(i)%q,P(j)%q,P(i)%uid,P(j)%uid,restart&
                &,rrxi,rrxj
           call error_md(msgstr)
        end if

        ! remember surface points to speed up convergence next time
        call store_surface_points(i,j,rxi,rxj)

        xi = rxi + P(i)%x
        xj = rxj + P(j)%x

        closest_surface_points_iterations = &
             &closest_surface_points_iterations + real(n_tot,kind=rk)
        closest_surface_points_count = closest_surface_points_count + 1.0_rk
    end subroutine compute_closest_surface_points

    !> find the respective two points in between of which the distance
    !> between the surfaces of two particles is smallest
    !>
    !> The algorithm closely follows the ideas in Lin&Han. SIAM
    !> J. Optim. 13, 298 (2002). It is implemented, however, based on
    !> direct geometric considerations for spheroids, thus avoiding
    !> the general description of ellipsoids as quadrics (which might
    !> or might not have some disadvantages with respect to
    !> performance). The additional exit condition on spheroid overlap
    !> is omitted which means that a meaningful (then negative) gap
    !> width and the corresponding direction is returned also in case
    !> of moderate overlap. The algorithm still fails if the overlap
    !> is large with respect to the radius of curvature.
    !>
    !> \param[in] i first particle index in \c P
    !>
    !> \param[in] j second particle index in \c P
    !>
    !> \param[in] rij minimum image vector from particle \c j to \c i
    !>
    !> \param[in] ropi half axes \c R_orth and \c R_para for particle \c i
    !>
    !> \param[in] ropj half axes \c R_orth and \c R_para for particle \c j
    !>
    !> \param[in] n_it_max if no convergence is achieved after this
    !> many iterations, this is considered a failure and the algorithm
    !> stops
    !>
    !> \param[in] local_radius compute iteration step based on local
    !> radius of curvature instead of the particle's global minimum
    !> radius of curvature. In some cases, this was found to speed up
    !> convergence tremendously but it is less stable. According to
    !> Lin&Han, the sphere around \c c1/2 touching the surface at \c
    !> x1/2 being fully inside is a sufficient requirement for
    !> convergence of the algorithm.
    !>
    !> \param[in] restart restart iteration with \c rrxi and \c rrxj
    !> as initial starting points instead of starting from the
    !> particle centers
    !>
    !> \param[in] rrxi with \c restart, initial point for iteration
    !> touching particle \c i, measured from its center
    !>
    !> \param[in] rrxj with \c restart, initial point for iteration
    !> touching particle \c j, measured from its center
    !>
    !> \param[out] rxi resulting point on particle \c i, measured from
    !> its center
    !>
    !> \param[out] rxj resulting point on particle \c j, measured from
    !> its center
    !>
    !> \param[out] g unit gap vector from \c i to \c j
    !>
    !> \param[out] gap gap width
    !>
    !> \param[out] success \c .true. if the computation converged
    !> without errors within \c n_it_max iterations
    !>
    !> \param[out] n_it actual number of iterations performed
    subroutine closest_surface_points(i,j,rij,ropi,ropj&
         &,n_it_max,local_radius,restart,rrxi,rrxj,rxi,rxj,g,gap,success,n_it)
        integer,intent(in) :: i,j,n_it_max
        real(kind=rk),intent(in) :: rij(3),ropi(2),ropj(2),rrxi(3),rrxj(3)
        logical,intent(in) :: local_radius,restart
        real(kind=rk),intent(out) :: rxi(3),rxj(3),g(3),gap
        logical,intent(out) :: success
        integer,intent(out) :: n_it

        ! This defines the divergence criterion. If it is chosen too
        ! weak, the planes spanned by the principal curvature
        ! directions on both surfaces will not be parallel which will
        ! make the result of gap_eigenproblem() differ for both
        ! particles. In consequence, actio/=reactio for the normal
        ! direction (probably more bad things might follow...).
        real(kind=rk),parameter :: esmall_csp=1.0e-6_rk
        real(kind=rk),parameter :: eso=1.0_rk-esmall_csp

        real(kind=rk) :: c1(3),p1(3),un1(3),x1(3)
        real(kind=rk) :: c2(3),p2(3),un2(3),x2(3)
        real(kind=rk) :: ac12,c12(3),ls1,ls2,R1,R2
        integer :: n

        ! assume particle i to be centered at the origin
        p1 = 0.0_rk
        p2 = -rij

        if (restart) then
           ! restart iteration from previous points on the surface
           un1 = spheroid_unit_normal(P(i), rrxi, ropi)
           un2 = spheroid_unit_normal(P(j), rrxj, ropj)

           if (local_radius) then
              R1 = minval(principal_curvature_radii(P(i),rrxi,ropi))
              R2 = minval(principal_curvature_radii(P(j),rrxj,ropj))
           else
              R1 = minimal_radius(ropi)
              R2 = minimal_radius(ropj)
           end if

           c1 = p1+rrxi-un1*R1
           c2 = p2+rrxj-un2*R2
        else
           ! start iteration from scratch, that means, from particle centers
           R1 = minimal_radius(ropi)
           R2 = minimal_radius(ropj)

           ! start iteration with spheres centered in each spheroid
           c1 = p1
           c2 = p2
        end if

        success = .false.
        n = 0
        do
           if (n>n_it_max) exit
           n = n+1

           ! line between the centers of both spheres
           c12 = c1-c2
           ac12 = norm(c12)

           ! its intersection points with the spheroids
           ls1 = line_spheroid_intersection(P(i), c1-p1, -c12, ropi)
           ls2 = line_spheroid_intersection(P(j), c2-p2, c12, ropj)
           ! for such amounts of overlap, c1 and c2 move beyond each
           ! other, so the direction of c12 flips and the present
           ! algorithm would fail
           if (ac12*(1.0_rk-(ls1+ls2)) <= -(R1+R2)) exit

           x1 = c1 - c12*ls1
           x2 = c2 + c12*ls2

           ! exit if normals in both points are close to parallel to c12
           un1 = spheroid_unit_normal(P(i), x1-p1, ropi)
           un2 = spheroid_unit_normal(P(j), x2-p2, ropj)
           if (eso*ac12 < dot_product(un2,c12).and.&
                &eso*ac12 < -dot_product(un1,c12)) then
              ! overlap and/or large R1/2 combined with unfortunate
              ! starting position could cause that an extremum is
              ! found but not the minimum. In this case, un1/2 are
              ! pointing away from the other particle but this is not
              ! noticed since also c12 points in the other
              ! direction. Do not accept such results.
              success = dot_product(un2,rij) > 0.0_rk
              exit
           end if

           if (local_radius) then
              R1 = minval(principal_curvature_radii(P(i),x1-p1,ropi))
              R2 = minval(principal_curvature_radii(P(j),x2-p2,ropj))
           end if

           ! otherwise move tangent spheres to now touch at x1 and x2
           c1 = x1-un1*R1
           c2 = x2-un2*R2
        end do
        n_it = n

        rxi = x1-p1
        rxj = x2-p2
        g = -c12/ac12
        gap = dot_product(x2-x1,g)

        ! if particles are (moderately) overlapping the surface points
        ! on each of the spheroids must mutually lie within the other
        ! spheroid's volume. Otherwise x1 and x2 are not the points of
        ! minimum distance but another extremum was found where the
        ! normals point in the wrong directions and the spheroids
        ! actually are behind each other.
        if (gap<0.0_rk.and.success) then
           success &
                &= in_spheroid(P(i),x2-p1,ropi).and.in_spheroid(P(j),x1-p2,ropj)
        end if
    end subroutine closest_surface_points

    !> calculate sub-lattice lubrication interactions for a pair of particles
    !>
    !> \param[in] rij2 squared minimum image distance
    !>
    !> \param[in] rij minimum image vector from particle \c j to \c i
    !>
    !> \param[in] i first particle index in \c P
    !>
    !> \param[in] j second particle index in \c P
    subroutine club_force_and_torque(rij2,rij,i,j)
        real(kind=rk),intent(in) :: rij2,rij(3)
        integer,intent(in) :: i,j
        real(kind=rk) :: g(3),gap,mu
        real(kind=rk) :: xi(3),pi(2,3),ci(2),cci(0:3),vi(3),wi(3),ropi(2)
        real(kind=rk) :: xj(3),pj(2,3),cj(2),ccj(0:3),vj(3),wj(3),ropj(2)
        real(kind=rk) :: fi(3),ti(3)
        real(kind=rk) :: fj(3),tj(3)
        real(kind=rk),parameter :: accepted_rel_asym=1.0e-2_rk
        real(kind=rk),parameter :: accepted_abs_asym=1.0e-12_rk
        real(kind=rk) :: ferr,favg,terr,tavg

        ! if this is called at time step 0 it means it's during
        ! growing stage, switch to a stripped-down routine then that
        ! computes only the conservative repulsion forces and no
        ! lubrication interactions
        if (nt==0) then
           call club_force_and_torque_growth(rij2,rij,i,j)
           return
        end if

        ! half axes R_orth and R_para for both particles
        ropi = particle_radii(P(i))
        ropj = particle_radii(P(j))

        ! particle i is always owned locally. If j is not owned, the
        ! same pair of particles is treated on j's owning process as
        ! well. In this case, compute interactions for both particles
        ! only on the process that owns the particle with the smaller
        ! uid. Since lubrication_cox makes sense only with MD Fluid
        ! Ladd, collect_forces can be relied on to send the forces
        ! back to the owner.
        if (j>npmax.and.P(i)%uid>P(j)%uid) return

        ! avoid expensive iterations if particles are too far apart
        if (contact_impossible(rij2,rij,i,j,ropi,ropj)) return

        ! find the two closest points on the surfaces of particle i
        ! and j and the resulting unit gap vector and width applying
        ! short-range clipping
        call compute_closest_surface_points(i,j,rij,ropi,ropj,xi,xj,g,gap)

        if (gap>=delta_c_max) return

        ! find unit directions of and principal curvatures themselves
        call principal_curvatures(i,xi,ropi,pi,ci)
        call principal_curvatures(j,xj,ropj,pj,cj)

        ! compute coefficients of 3rd order correction to local
        ! surface shape
        call cubic_coefficients(i,xi,ropi,cci)
        call cubic_coefficients(j,xj,ropj,ccj)

        ! compute local speed of surface translation and rotation
        call surface_velocity(i,xi,vi,wi)
        call surface_velocity(j,xj,vj,wj)

        call dynamic_viscosity(i,j,mu)

        ! in principle fi=-fj and ti=-tj. One could save computations
        ! and enforce perfect momentum conservation using this. At the
        ! moment, this is not done but the 2nd computation is used as
        ! consistency check instead.
        call force_and_torque(mu,gap,ci,cj,pi,pj,cci,ccj,g,vj-vi,wj-wi,fi,ti)
        call force_and_torque(mu,gap,cj,ci,pj,pi,ccj,cci,-g,vi-vj,wi-wj,fj,tj)

        call apply_force_and_torque(i,xi,fi,ti)
        call apply_force_and_torque(j,xj,fj,tj)

        ferr = norm(fi+fj)
        favg = 0.5_rk*norm(fi-fj)
        if (ferr>accepted_rel_asym*favg.and.ferr>accepted_abs_asym) then
           write (msgstr,'("club_force_and_torque(): poor anti-symmetry of'&
                &//' force; relative/absolute error=",ES15.8,", ",ES15.8,"; '&
                &//'uidi=",I0,", uidj=",I0)') ferr/favg,ferr,P(i)%uid,P(j)%uid
           call log_msg_md(msgstr,forall=.true.)
        end if
        terr = norm(ti+tj)
        tavg = 0.5_rk*norm(ti-tj)
        if (terr>accepted_rel_asym*tavg.and.terr>accepted_abs_asym) then
           write (msgstr,'("club_force_and_torque(): poor anti-symmetry of'&
                &//' torque; relative/absolute error=",ES15.8,", ",ES15.8,"; '&
                &//'uidi=",I0,", uidj=",I0)') terr/tavg,terr,P(i)%uid,P(j)%uid
           call log_msg_md(msgstr,forall=.true.)
        end if

        if (gap<0.0_rk) then
           write (msgstr,'("club_force_and_torque(): '&
                &//'particle overlap detected. '&
                &//'gap=",ES15.8,", '&
                &//'uidi=",I0,", uidj=",I0)') &
                &gap,P(i)%uid,P(j)%uid
           call log_msg_md(msgstr,forall=.true.)
        end if
    end subroutine club_force_and_torque

    !> stripped-down version of \c club_force_and_torque() that
    !> computes only the conservative repulsion force which is used
    !> during particle growth
    !>
    !> \param[in] rij2 squared minimum image distance
    !>
    !> \param[in] rij minimum image vector from particle \c j to \c i
    !>
    !> \param[in] i first particle index in \c P
    !>
    !> \param[in] j second particle index in \c P
    subroutine club_force_and_torque_growth(rij2,rij,i,j)
        real(kind=rk),intent(in) :: rij2,rij(3)
        integer,intent(in) :: i,j

        real(kind=rk) :: f(3),g(3),gap
        real(kind=rk) :: xi(3),ropi(2)
        real(kind=rk) :: xj(3),ropj(2)

        ! half axes R_orth and R_para for both particles
        ropi = particle_radii(P(i)) * growth_factor
        ropj = particle_radii(P(j)) * growth_factor

        ! particle i is always owned locally. If j is not owned, the
        ! same pair of particles is treated on j's owning process as
        ! well. In this case, compute interactions for both particles
        ! only on the process that owns the particle with the smaller
        ! uid. Since lubrication_cox makes sense only with MD Fluid
        ! Ladd, collect_forces can be relied on to send the forces
        ! back to the owner.
        if (j>npmax.and.P(i)%uid>P(j)%uid) return

        ! avoid expensive iterations if particles are too far apart
        if (contact_impossible(rij2,rij,i,j,ropi,ropj)) return

        ! find the two closest points on the surfaces of particle i
        ! and j and the resulting unit gap vector and width applying
        ! short-range clipping
        call compute_closest_surface_points(i,j,rij,ropi,ropj,xi,xj,g,gap)

        ! compute only repulsive short-range force
        if (gap < lubrication_clip_dist) then
           f = g * lubrication_stiffness*(lubrication_clip_dist-gap)

           P(i)%f = P(i)%f - f
           P(i)%t = P(i)%t - cross_product(xi-P(i)%x,f)

           P(j)%f = P(j)%f + f
           P(j)%t = P(j)%t + cross_product(xj-P(j)%x,f)
        end if
    end subroutine club_force_and_torque_growth

    !> sets \c growth_factor to a new value
    !>
    !> \param[in] f new \c growth_factor
    subroutine club_set_growth_factor(f)
        real(kind=rk),intent(in) :: f

        growth_factor = f
    end subroutine club_set_growth_factor

    !> setup stage for \c lbe_md_fluid_ladd_club_module
    subroutine club_setup
        ! maximum upper gap cutoff
        delta_c_max = max(delta_c,delta_c_T,delta_c_R)

        call Mi6r_init(surface_points)
        call check_allocate(Mi6r_provide_capacity(surface_points&
             &,lubrication_cox_max_store_points),'club_setup(): surface_points')
        surface_points_reset_time = nt

        ! This assumes that there will never be a higer particle uid
        ! during the simulation!!!
        max_particle_uid = floor(sqrt(real(huge(max_particle_uid),kind=rk)))
        write (msgstr,fmt='("club_setup(): assuming maximum particle uid ever '&
             &//'reached to be ",I0," (or smaller)")') max_particle_uid
        call log_msg_md(msgstr)

#ifndef SINGLEFLUID
        call error_md('At the moment, lubrication_cox requires SINGLEFLUID---'&
             &//'set lubrication_cox=.false. or recompile with SINGLEFLUID')
#endif
    end subroutine club_setup

    !> prints conclusive information concerning \c lubrication_cox
    !>
    !> \param[in] units vector of output units to which information is
    !> printed
    subroutine club_summary_fluid(units)
        integer,intent(in) :: units(:)
        integer :: ierror,gmin,u
        real(kind=rk) :: rbuf(2),sbuf(2)

        sbuf = &
             &(/closest_surface_points_count,closest_surface_points_iterations/)
        call MPI_Reduce(sbuf,rbuf,2,LBE_REAL,MPI_SUM,0,comm_cart,ierror)

        call MPI_Reduce(min_surface_points_lifetime,gmin,1,MPI_INTEGER,MPI_MIN&
             &,0,comm_cart,ierror)

        rank0: if (myrankc==0) then
           do u = 1,size(units)
              write (unit=units(u),fmt='("MD Fluid Ladd lubrication_cox '&
                   &//'needed to determine the minimum particle distance '&
                   &//'",ES15.8," times.")') rbuf(1)
              write (unit=units(u),fmt='("On average, this required '&
                   &//'",ES15.8," iterations.")') rbuf(2)/rbuf(1)
              write (unit=units(u),fmt='("Minimum lifetime of stored surface '&
                   &//'points was ",I0," time steps. '&
                   &//'lubrication_cox_max_store_points=",I0)') &
                   &gmin,lubrication_cox_max_store_points
              write (unit=units(u),fmt='()')
           end do
        end if rank0
    end subroutine club_summary_fluid

    !> try to exclude small gaps within two particles in an
    !> inexpensive way
    !>
    !> \param[in] rij2 squared minimum image distance
    !>
    !> \param[in] rij minimum image vector from particle \c j to \c i
    !>
    !> \param[in] i first particle index in \c P
    !>
    !> \param[in] j second particle index in \c P
    !>
    !> \param[in] ropi half axes \c R_orth and \c R_para for particle \c i
    !>
    !> \param[in] ropj half axes \c R_orth and \c R_para for particle \c j
    !>
    !> \returns \c .true. if a gap small enough to make lubrication
    !> corrections necessary can be excluded, \c .false. otherwise
    pure logical function contact_impossible(rij2,rij,i,j,ropi,ropj)
        real(kind=rk),intent(in) :: rij2,rij(3),ropi(2),ropj(2)
        integer,intent(in) :: i,j
        real(kind=rk) :: arij,n_i(3),n_j(3),R,R_plus_delta_i,R_plus_delta_j&
             &,urij(3)
        logical intersection_test_trustworthy_i,intersection_test_trustworthy_j
        logical ret

        R_plus_delta_i = maxval(ropi) + delta_c_max
        R_plus_delta_j = maxval(ropj) + delta_c_max

        arij = sqrt(rij2)
        urij = rij/arij

        ! if particles are closer, the
        ! plane_spheroid_intersection_exists() below might fail
        ! erroneously (checking for intersections with a plane on the
        ! other side of the particle)
        intersection_test_trustworthy_i = arij >= R_plus_delta_i
        intersection_test_trustworthy_j = arij >= R_plus_delta_j

        ret = .false.
        ! if just one of the particles does not intersect with a plane
        ! normal to rij in a distance of the other particle's
        ! R_plus_delta away from the other particle, it is certainly
        ! impossible that the gap is smaller than delta_c[_T/R]
        if (intersection_test_trustworthy_j) then
           n_i = -rij + R_plus_delta_j*urij
           ret = .not.plane_spheroid_intersection_exists(P(i),n_i,ropi)
        end if
        if (intersection_test_trustworthy_i.and..not.ret) then
           ! do this  only if a contact was not excluded before already
           n_j = rij - R_plus_delta_i*urij
           ret = .not.plane_spheroid_intersection_exists(P(j),n_j,ropj)
        end if
        contact_impossible = ret
    end function contact_impossible

    !> compute 3rd order coefficients in the expansion of the surface
    !> of a spheroidal particle
    !>
    !> \param[in] k particle index in \c P
    !>
    !> \param[in] x point on the surface of particle \c k
    !>
    !> \param[in] rop half axes \c R_orth and \c R_para for \c P(k)
    !>
    !> \param[out] cc 3rd order coefficients (called
    !> \f$\Gamma_{0--3}\f$ in Claeys&Brady 1989)
    !>
    !> The sagittal direction is identified with \f$x_1\f$ in
    !> Claeys&Brady 1989 and it is assumed that its angle with \c
    !> P(k)%o is 90 degree or less as defined in \c
    !> principal_curvatures() . Such an assumption is important to
    !> define the sign of the cubic coefficients properly.
    subroutine cubic_coefficients(k,x,rop,cc)
        integer,intent(in) :: k
        real(kind=rk),intent(in) :: x(3),rop(2)
        real(kind=rk),intent(out) :: cc(0:3)
        real(kind=rk) :: aq,it0sq,q(3),r(3),ro,rq,ruq,r2diff,r_o,r_p

        r_o = rop(1)
        r_p = rop(2)

        r = x-P(k)%x
        ro = dot_product(r,P(k)%o)
        q = r - ro*P(k)%o
        rq = dot_product(r,q)
        aq = norm(q)

        if (aq>esmall) then
           ruq = rq/aq
        else
           ! cure singularity for small q in the same way as in
           ! principal_curvatures(). At such small q also the cubic
           ! coefficients will be zero for symmetry reasons anyway.
           ruq = 0.0_rk
        end if

        r2diff = r_p**2 - r_o**2
        it0sq = 1.0_rk/(r_o**2 + r2diff*(ruq/r_o)**2)

        ! these coefficients were obtained from a 3rd order expansion
        ! of the spheroid surface. Note that while ruq is always
        ! positive, ro changes sign at the equatorial plane. Different
        ! from spheroid_unit_normal(), the sign change is desired
        ! here.
        cc(0) = 0.5_rk*r2diff*ro*ruq*it0sq**3
        cc(1) = 0.0_rk
        cc(2) = 0.5_rk*r2diff*ro*ruq*it0sq**2/r_o**2
        cc(3) = 0.0_rk
    end subroutine cubic_coefficients

    !> obtains the dynamic fluid viscosity assumed in the gap between
    !> two particles
    !>
    !> \param[in] i 1st particle index in \c P
    !>
    !> \param[in] i 2nd particle index in \c P
    !>
    !> \param[out] mu dynamic viscosity
    subroutine dynamic_viscosity(i,j,mu)
        integer,intent(in) :: i,j
        real(kind=rk),intent(out) :: mu
        real(kind=rk) :: nu_r,rhof_r

        ! kinematic viscosity
        nu_r = (2.0_rk*tau_r-1.0_rk)/6.0_rk

        ! local fluid density
#ifdef LADD_SURR_RHOF
        rhof_r = 0.5_rk*(P(i)%rhof+P(j)%rhof)
#else
        rhof_r = pfr
#endif

        mu = rhof_r*nu_r
    end subroutine dynamic_viscosity

    !> calculate effect of lubrication interactions on one particle based on
    !> the knowledge of the facing surface points
    !>
    !> \param[in] mu local dynamic viscosity
    !>
    !> \param[in] gap gap width
    !>
    !> \param[in] c principal curvatures
    !>
    !> \param[in] co principal curvatures of other particle
    !>
    !> \param[in] pp unit vectors in direction of principal curvatures
    !>
    !> \param[in] ppo unit vectors in direction of principal
    !> curvatures of other particle
    !>
    !> \param[in] cc 3rd order coefficients (called
    !> \f$\Gamma_{0--3}\f$ in Claeys&Brady 1989)
    !>
    !> \param[in] cc 3rd order coefficients of other particle
    !>
    !> \param[in] g unit surface normal, forms right-handed
    !> orthonormal base together with \c pp
    !>
    !> \param[in] dv surface velocity difference defined as in Claeys&Brady
    !> 1989: \f$\mathbf{v}_j-\mathbf{v}_i\f$
    !>
    !> \param[in] dw surface angular velocity difference defined as in
    !> Claeys&Brady 1989: \f$\mathbf{\omega}_j-\mathbf{\omega}_i\f$
    !>
    !> \param[out] f_ret resulting force in point of closest approach
    !>
    !> \param[out] t_ret resulting torque
    subroutine force_and_torque(mu,gap,c,co,pp,ppo,cc,cco,g,dv,dw,f_ret,t_ret)
        real(kind=rk),intent(in) :: mu,gap,c(2),co(2),pp(2,3),ppo(2,3)&
             &,cc(0:3),cco(0:3),g(3),dv(3),dw(3)
        real(kind=rk),intent(out) :: f_ret(3),t_ret(3)
        real(kind=rk) :: f(3),t(3),l1,l2,V(6),Chi(6,6),cgap,f_rep,sinp,cosp&
             &,sinc,cosc
        real(kind=rk) :: gammap(0:3),beta(0:3),K(0:3),kappa(0:3)
        integer :: i

        ! for small enough gaps, clip gap and apply repulsive force
        if (gap<lubrication_clip_dist) then
           f_rep = lubrication_stiffness*(lubrication_clip_dist-max(0.0_rk,gap))
           cgap = lubrication_clip_dist
        else
           f_rep = 0.0_rk
           cgap = gap
        end if

        ! find eigenvalues (l1, l2) of the 2nd order-expanded gap
        ! height, a quadratic form of the position in the plane
        ! defined by pp, and the angle to the related eigenvectors
        ! (sinc, cosc) and to the other particle's directions of
        ! principal curvature (sinp, cosp)
        call gap_eigenproblem(c,co,pp,ppo,sinp,cosp,l1,l2,sinc,cosc)
        ! swap order of cubic coefficients of other surface because
        ! also the principal directions ppo are interpreted in
        ! gap_eigenproblem() as being swapped. A positive cubic
        ! coefficient means that the gap grows stronger than
        ! quadratically in the corresponding direction. Since both
        ! cubic coefficients and directions are reordered, this
        ! relation is kept and no additional sign changes are
        ! required.
        do i=0,3
           gammap(i) = cco(3-i)
        end do
        ! cubic coefficients of other surface in own coordinate frame
        call transform_cubic_coefficients(gammap,sinp,cosp,beta)
        ! combined cubic coefficients in eigenframe
        call transform_cubic_coefficients(cc+beta,sinc,cosc,K)
        do i=0,3
           kappa(i) = K(i) / (sqrt(l1)**(3-i) * sqrt(l2)**i)
        end do

        ! translation and rotation velocity vector in gap coordinate
        ! system
        V(1) = dot_product(dv,pp(1,:))
        V(2) = dot_product(dv,pp(2,:))
        V(3) = dot_product(dv,g)
        V(4) = dot_product(dw,pp(1,:))
        V(5) = dot_product(dw,pp(2,:))
        V(6) = dot_product(dw,g)

        ! diverging components of the resistance tensor Chi. The sign
        ! is consistent with Cox 1974 but opposite to what one might
        ! expect from Claeys&Brady 1989 because there the stress on
        ! the fluid is given, here it's the stress on the particle.

        ! The three cutoffs delta_c{_[TR]} are used in the same way as
        ! in Ladd's Susp3D and described (less precisely, however) in
        ! Nguyen&Ladd 2002: The "rotational" cutoff is employed only
        ! for the coupling of rotations to torques, the "tangential"
        ! one for everything else which is not a normal force due to a
        ! normal translation.

        f = 0.0_rk
        t = 0.0_rk

        if (cgap<delta_c) then
           ! from claeys89, consistent with cox74
           Chi(3,3) = (1.0_rk/cgap-inv_delta_c) &
                &* 3.0_rk*pi/(sqrt(l1*l2)*(l1+l2))

           f = f + mu*g*Chi(3,3)*V(3) - g*f_rep
        end if

        if (lubrication_tangential.and.cgap<delta_c_T) then
           ! these were computed only in claeys89
           Chi(1,3) = -(log(cgap)-log_delta_c_T) * 1.5_rk*pi * (&
                &(2.0_rk*sqrt(l1)*(3.0_rk*kappa(0)*l1+kappa(2)*l2)&
                &-((7.0_rk*l1+2.0_rk*l2)*kappa(0)+(l1+2.0_rk*l2)*kappa(2))&
                &*1.5_rk*c(1)/sqrt(l1)) * cosc/(3.0_rk*l1+2.0_rk*l2)&
                &+(2.0_rk*sqrt(l2)*(kappa(1)*l1+3.0_rk*kappa(3)*l2)&
                &-((2.0_rk*l1+l2)*kappa(1)+(2.0_rk*l1+7.0_rk*l2)*kappa(3))&
                &*1.5_rk*c(1)/sqrt(l2)) * sinc/(2.0_rk*l1+3.0_rk*l2)&
                &+3.0_rk*cc(0)*(cosc**2/l1+sinc**2/l2)&
                &+2.0_rk*cc(1)*sinc*cosc*(1.0_rk/l2-1.0_rk/l1)&
                &+cc(2)*(sinc**2/l1+cosc**2/l2)&
                &) / (sqrt(l1*l2)*(l1+l2))
           Chi(2,3) = -(log(cgap)-log_delta_c_T) * 1.5_rk*pi * (&
                &(-2.0_rk*sqrt(l1)*(3.0_rk*kappa(0)*l1+kappa(2)*l2)&
                &+((7.0_rk*l1+2.0_rk*l2)*kappa(0)+(l1+2.0_rk*l2)*kappa(2))&
                &*1.5_rk*c(2)/sqrt(l1)) * sinc/(3.0_rk*l1+2.0_rk*l2)&
                &+(2.0_rk*sqrt(l2)*(kappa(1)*l1+3.0_rk*kappa(3)*l2)&
                &-((2.0_rk*l1+l2)*kappa(1)+(2.0_rk*l1+7.0_rk*l2)*kappa(3))&
                &*1.5_rk*c(2)/sqrt(l2)) * cosc/(2.0_rk*l1+3.0_rk*l2)&
                &+3.0_rk*cc(3)*(sinc**2/l1+cosc**2/l2)&
                &+2.0_rk*cc(2)*sinc*cosc*(1.0_rk/l2-1.0_rk/l1)&
                &+cc(1)*(cosc**2/l1+sinc**2/l2)&
                &) / (sqrt(l1*l2)*(l1+l2))
           Chi(4,3) = (log(cgap)-log_delta_c_T) * 2.25_rk*pi * (&
                &((2.0_rk*l1+l2)*kappa(1) + (2.0_rk*l1+7.0_rk*l2)*kappa(3))&
                &*cosc / (sqrt(l2)*(2.0_rk*l1+3.0_rk*l2))&
                &-((l1+2.0_rk*l2)*kappa(2) + (7.0_rk*l1+2.0_rk*l2)*kappa(0))&
                &*sinc / (sqrt(l1)*(3.0_rk*l1+2.0_rk*l2))&
                &) / (sqrt(l1*l2)*(l1+l2))
           Chi(5,3) = -(log(cgap)-log_delta_c_T) * 2.25_rk*pi * (&
                &((l1+2.0_rk*l2)*kappa(2) + (7.0_rk*l1+2.0_rk*l2)*kappa(0))&
                &*cosc / (sqrt(l1)*(3.0_rk*l1+2.0_rk*l2))&
                &+((2.0_rk*l1+l2)*kappa(1) + (2.0_rk*l1+7.0_rk*l2)*kappa(3))&
                &*sinc / (sqrt(l2)*(2.0_rk*l1+3.0_rk*l2))&
                &) / (sqrt(l1*l2)*(l1+l2))

           ! from claeys89, consistent with cox74 which is identical
           ! to what can be obtained later from claeys89 as Chi(3,6)
           Chi(6,3) = -(log(cgap)-log_delta_c_T) * 1.5_rk*pi * sinc*cosc * &
                &(c(1)-c(2))*(1.0_rk/l1-1.0_rk/l2) / (sqrt(l1*l2)*(l1+l2))

           ! these are from cox74, check with notation in claeys89 pending...
           Chi(1,1) = -(log(cgap)-log_delta_c_T) * pi*c(1)**2 * (&
                &3.0_rk*cosc**2*(1.0_rk-l1/c(1))**2 &
                &/ ((3.0_rk*l1+2.0_rk*l2)*l1)&
                &+3.0_rk*sinc**2*(1.0_rk-l2/c(1))**2 &
                &/ ((2.0_rk*l1+3.0_rk*l2)*l2)&
                &+1.0_rk/c(1)**2&
                &) / sqrt(l1*l2)
           Chi(1,2) = (log(cgap)-log_delta_c_T) &
                &* 3.0_rk*pi*sinc*cosc*c(1)*c(2) *(&
                &(1.0_rk-l1/c(1))*(1.0_rk-l1/c(2)) &
                &/ ((3.0_rk*l1+2.0_rk*l2)*l1)&
                &-(1.0_rk-l2/c(1))*(1.0_rk-l2/c(2)) &
                &/ ((2.0_rk*l1+3.0_rk*l2)*l2)&
                &) / sqrt(l1*l2)
           Chi(1,4) = (log(cgap)-log_delta_c_T) * 3.0_rk*pi*sinc*cosc*c(1) *(&
                &(1.0_rk-l1/c(1)) / ((3.0_rk*l1+2.0_rk*l2)*l1)&
                &-(1.0_rk-l2/c(1)) / ((2.0_rk*l1+3.0_rk*l2)*l2)&
                &) / sqrt(l1*l2)
           Chi(1,5) = (log(cgap)-log_delta_c_T) * 3.0_rk*pi*c(1) * (&
                &cosc**2*(1.0_rk-l1/c(1)) / ((3.0_rk*l1+2.0_rk*l2)*l1)&
                &+sinc**2*(1.0_rk-l2/c(1)) / ((2.0_rk*l1+3.0_rk*l2)*l2)&
                &) / sqrt(l1*l2)
           Chi(2,2) = -(log(cgap)-log_delta_c_T) * pi*c(2)**2 * (&
                &3.0_rk*sinc**2*(1.0_rk-l1/c(2))**2 &
                &/ ((3.0_rk*l1+2.0_rk*l2)*l1)&
                &+3.0_rk*cosc**2*(1.0_rk-l2/c(2))**2 &
                &/ ((2.0_rk*l1+3.0_rk*l2)*l2)&
                &+1.0_rk/c(2)**2&
                &) / sqrt(l1*l2)
           Chi(2,4) = -(log(cgap)-log_delta_c_T) * 3.0_rk*pi*c(2) * (&
                &sinc**2*(1.0_rk-l1/c(2)) / ((3.0_rk*l1+2.0_rk*l2)*l1)&
                &+cosc**2*(1.0_rk-l2/c(2)) / ((2.0_rk*l1+3.0_rk*l2)*l2)&
                &) / sqrt(l1*l2)
           Chi(2,5) = (log(cgap)-log_delta_c_T) * 3.0_rk*pi*sinc*cosc*c(2) * (&
                &-(1.0_rk-l1/c(2)) / ((3.0_rk*l1+2.0_rk*l2)*l1)&
                &+(1.0_rk-l2/c(2)) / ((2.0_rk*l1+3.0_rk*l2)*l2)&
                &) / sqrt(l1*l2)

           ! exploit symmetry of Chi
           Chi(2,1) = Chi(1,2)
           Chi(3,1) = Chi(1,3)
           Chi(3,2) = Chi(2,3)
           Chi(3,4) = Chi(4,3)
           Chi(3,5) = Chi(5,3)
           Chi(3,6) = Chi(6,3)
           Chi(4,1) = Chi(1,4)
           Chi(4,2) = Chi(2,4)
           Chi(5,1) = Chi(1,5)
           Chi(5,2) = Chi(2,5)

           f = f + mu*pp(1,:)*Chi(1,1)*V(1)
           f = f + mu*pp(1,:)*Chi(1,2)*V(2)
           f = f + mu*pp(1,:)*Chi(1,3)*V(3)
           f = f + mu*pp(1,:)*Chi(1,4)*V(4)
           f = f + mu*pp(1,:)*Chi(1,5)*V(5)

           f = f + mu*pp(2,:)*Chi(2,1)*V(1)
           f = f + mu*pp(2,:)*Chi(2,2)*V(2)
           f = f + mu*pp(2,:)*Chi(2,3)*V(3)
           f = f + mu*pp(2,:)*Chi(2,4)*V(4)
           f = f + mu*pp(2,:)*Chi(2,5)*V(5)

           f = f + mu*g*Chi(3,1)*V(1)
           f = f + mu*g*Chi(3,2)*V(2)
           f = f + mu*g*Chi(3,4)*V(4)
           f = f + mu*g*Chi(3,5)*V(5)
           f = f + mu*g*Chi(3,6)*V(6)

           t = t + mu*pp(1,:)*Chi(4,1)*V(1)
           t = t + mu*pp(1,:)*Chi(4,2)*V(2)
           t = t + mu*pp(1,:)*Chi(4,3)*V(3)

           t = t + mu*pp(2,:)*Chi(5,1)*V(1)
           t = t + mu*pp(2,:)*Chi(5,2)*V(2)
           t = t + mu*pp(2,:)*Chi(5,3)*V(3)

           t = t + mu*g*Chi(6,3)*V(3)
        end if

        if (lubrication_tangential.and.cgap<delta_c_R) then
           ! these are from cox74, check with notation in claeys89 pending...
           Chi(4,4) = -(log(cgap)-log_delta_c_R) * 3.0_rk*pi * (&
                &sinc**2 / ((3.0_rk*l1+2.0_rk*l2)*l1)&
                &+cosc**2 / ((2.0_rk*l1+3.0_rk*l2)*l2)&
                &) / sqrt(l1*l2)
           Chi(4,5) = (log(cgap)-log_delta_c_R) * 3.0_rk*pi*sinc*cosc * (&
                &-1.0_rk / ((3.0_rk*l1+2.0_rk*l2)*l1)&
                &+1.0_rk / ((2.0_rk*l1+3.0_rk*l2)*l2)&
                &) / sqrt(l1*l2)
           Chi(5,5) = -(log(cgap)-log_delta_c_R) * 3.0_rk*pi * (&
                &cosc**2 / ((3.0_rk*l1+2.0_rk*l2)*l1)&
                &+sinc**2 / ((2.0_rk*l1+3.0_rk*l2)*l2)&
                &) / sqrt(l1*l2)

           Chi(5,4) = Chi(4,5)

           t = t + mu*pp(1,:)*Chi(4,4)*V(4)
           t = t + mu*pp(1,:)*Chi(4,5)*V(5)
           t = t + mu*pp(2,:)*Chi(5,4)*V(4)
           t = t + mu*pp(2,:)*Chi(5,5)*V(5)
        end if

        f_ret = f
        t_ret = t
    end subroutine force_and_torque

    !> finds the curvature eigenvalues of the gap height expressed as
    !> quadratic form and the angle between the eigensystem and the
    !> principal curvatures
    !>
    !> \param[in] c principal curvatures
    !>
    !> \param[in] co principal curvatures of other particle
    !>
    !> \param[in] pp unit vectors in direction of principal curvatures
    !>
    !> \param[in] ppo unit vectors in direction of principal
    !> curvatures of other particle
    !>
    !> \param[out] sinp sine of the angle between both particles'
    !> directions of principal curvature
    !>
    !> \param[out] cosp cosine of the angle between both particles'
    !> directions of principal curvature
    !>
    !> \param[out] l1 1st curvature eigenvalue
    !>
    !> \param[out] l2 2nd curvature eigenvalue
    !>
    !> \param[out] sinc sine of the angle of the rotation matrix that
    !> transforms to the eigensystem
    !>
    !> \param[out] cosc cosine of the angle of the rotation matrix
    !> that transforms to the eigensystem
    subroutine gap_eigenproblem(c,co,pp,ppo,sinp,cosp,l1,l2,sinc,cosc)
        real(kind=rk),intent(in) :: c(2),co(2),pp(2,3),ppo(2,3)
        real(kind=rk),intent(out) :: sinp,cosp,l1,l2,sinc,cosc

        real(kind=rk) :: B0,B1,B2,C11,C12,C21,C22,ll1,ll2,inorm

        ! In the following, swap the two elements of each co and ppo
        ! with respect to the indices found for the other surface's
        ! curvature in the literature. This is because ppo and the
        ! inverse gap vector form a right-handed base, so the swapped
        ! ppo and the gap vector itself (which is used when looking at
        ! this particle) also is one.

        ! cosine and sine of the angle between this particle's and the
        ! other particle's directions of principal curvature, that is,
        ! with the above comment, the angle between pp(1,:) and
        ! ppo(2,:) or between pp(2,:) and ppo(1,:). cosine is trivial,
        ! sine can be deducted easily.
        cosp = dot_product(pp(1,:),ppo(2,:))
        sinp = -dot_product(pp(1,:),ppo(1,:))

        ! B as in Claeys&Brady 1989: 2nd order coefficients of other
        ! particle's surface expressed in system of this particle's
        ! principal curvature directions (note comment above on
        ! swapped co)
        B0 = 0.5_rk*(cosp**2*co(2) + sinp**2*co(1))
        B1 = (co(2) - co(1))*sinp*cosp
        B2 = 0.5_rk*(sinp**2*co(2) + cosp**2*co(1))

        ! C matrix as in Claeys&Brady 1989: describing the gap height
        ! as a quadratic form
        C11 = 0.5_rk*c(1) + B0
        C12 = 0.5_rk*B1
        C21 = C12
        C22 = 0.5_rk*c(2) + B2

        ! eigenvalues of a 2x2 matrix
        ll1 = 0.5_rk*(C11+C22)
        ll2 = 0.5_rk*sqrt((C11-C22)**2+4.0_rk*C12*C21)
        l1 = ll1-ll2
        l2 = ll1+ll2

        if (abs(ll2)<esmall) then
           ! (numerically) zero discriminant ll2 indicates that the
           ! eigenvectors are indeterminate. In this case, arbitrarily
           ! assume angle chi to be zero.
           sinc = 0.0_rk
           cosc = 1.0_rk
        else
           ! normalized eigenvector to l2, it's the second column of
           ! the transformation matrix in (2.6) of Claeys&Brady 1989,
           ! so its components are directly sin and cos of the angle
           ! chi.
           inorm = 1.0_rk/sqrt(C12**2 + (l2 - C11)**2)
           sinc = C12 * inorm
           cosc = (l2 - C11) * inorm
        end if
    end subroutine gap_eigenproblem

    !> checks whether a position is inside a spheroidal particle
    !>
    !> \param[in] p particle
    !>
    !> \param[in] rx position to check, measured from the spheroid's
    !> center
    !>
    !> \param[in] rop half axes \c R_orth and \c R_para for \c p
    !>
    !> \returns \c .true. if \c rx is located within the spheroid
    !> (excluding positions exactly on its surface), \c
    !> .false. otherwise
    pure logical function in_spheroid(p,rx,rop)
        type(md_particle_type),intent(in) :: p
        real(kind=rk),intent(in) :: rx(3),rop(2)
        real(kind=rk) :: r(3),rho(3),ro,rp,u,Q

        ro = rop(1)
        rp = rop(2)

        u = dot_product(p%o,rx)
        rho = rx-u*p%o
        Q = u*ro/rp
        in_spheroid = Q**2 + dot_product(rho,rho) < ro**2
    end function in_spheroid

    !> returns the distance from a point inside a spheroidal particle
    !> to its surface, measured along and in units of a given
    !> direction
    !>
    !> \param[in] p particle
    !>
    !> \param[in] rx position of the point measured from the
    !> particle's center
    !>
    !> \param[in] v direction in which to look for an intersection,
    !> starting from \c rx
    !>
    !> \param[in] rop half axes \c R_orth and \c R_para for \c p
    !>
    !> \returns \c w with \c rx+w*v hitting the spheroid surface
    pure real(kind=rk) function line_spheroid_intersection(p,rx,v,rop)
        type(md_particle_type),intent(in) :: p
        real(kind=rk),intent(in) :: rx(3),v(3),rop(2)
        real(kind=rk) :: A(3,3),vb(3),vbs(3),rb(3),rbs(3),inv_R(3),iQQ,PP

        ! v and rx in body-fixed coordinates
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
        vb = matmul(A,v)
        rb = matmul(A,rx)

        ! rescale all lengths to the case of a unit sphere for the
        ! sake of convenience and performance
        ! inverse particle half-axes
        inv_R = (/1.0_rk/rop(1),1.0_rk/rop(1),1.0_rk/rop(2)/)

        vbs = vb*inv_R
        rbs = rb*inv_R

        ! this one can be deducted setting the equations for a sphere
        ! and a line equal
        iQQ = 1.0_rk/dot_product(vbs,vbs)
        PP = dot_product(vbs,rbs)*iQQ
        line_spheroid_intersection &
             &= -PP+sqrt(PP*PP-(dot_product(rbs,rbs)-1.0_rk)*iQQ)
    end function line_spheroid_intersection

    !> returns the minimal radius of curvature for a particle
    !>
    !> \param[in] rop half axes \c R_orth and \c R_para of the
    !> particle
    !>
    !> \returns minimal radius of curvature
    pure real(kind=rk) function minimal_radius(rop)
        real(kind=rk),intent(in) :: rop(2)
        real(kind=rk) :: ro,rp

        ro = rop(1)
        rp = rop(2)

        ! for a spheroid, the minimal radius of curvature is realized
        ! either at the tip (prolate) or at the equatorial
        ! circumference in tangential direction (oblate); also see \c
        ! principal_curvatures() for reference
        if (rp<=ro) then
           ! oblate (or sphere)
           minimal_radius = rp**2/ro
        else
           ! prolate
           minimal_radius = ro**2/rp
        end if
    end function minimal_radius

    !> creates a unique identifier for a pair of particles
    !>
    !> \param[in] i index of 1st particle in \c P
    !>
    !> \param[in] j index of 2nd particle in \c P
    pure integer function pair_uid(i,j)
        integer,intent(in) :: i,j
        integer :: ui,uj

        ui = P(i)%uid
        uj = P(j)%uid
        pair_uid = min(ui,uj)*max_particle_uid + max(ui,uj)
    end function pair_uid

    !> checks whether an intersection of a spheroidal particle with a
    !> given plane exists
    !>
    !> \param[in] p particle
    !>
    !> \param[in] n distance (and therfore normal) vector to the
    !> plane, measured from the spheroid's center
    !>
    !> \param[in] rop half axes \c R_orth and \c R_para for \c p
    !>
    !> \returns \c .true. if an intersection exists, \c
    !> .false. otherwise
    pure logical function plane_spheroid_intersection_exists(p,n,rop)
        type(md_particle_type),intent(in) :: p
        real(kind=rk),intent(in) :: n(3),rop(2)
        real(kind=rk) :: A(3,3),nb(3)

        ! n in body-fixed coordinates
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
        nb = matmul(A,n)

        plane_spheroid_intersection_exists = &
             &sum((nb*(/rop(1),rop(1),rop(2)/))**2) >= dot_product(nb,nb)**2
    end function plane_spheroid_intersection_exists

    !> compute radii of principal curvatures at a given point of a
    !> particle surface
    !>
    !> \param[in] p particle
    !>
    !> \param[in] rx point on particle surface measured from its
    !> center
    !>
    !> \param[in] rop half axes \c R_orth and \c R_para for \c p
    !>
    !> radius of curvature and terminology follow Harris,
    !> Ophthal. Physiol. Opt. 26, 497 (2006) but could be computed
    !> from a local 2nd order expansion of the surface relatively
    !> easily as well.
    pure function principal_curvature_radii(p,rx,rop)
        real(kind=rk) :: principal_curvature_radii(2)
        type(md_particle_type),intent(in) :: p
        real(kind=rk),intent(in) :: rx(3),rop(2)

        real(kind=rk) :: aq,q(3),ro,rp,rq,ruq

        ro = rop(1)
        rp = rop(2)

        q = rx - dot_product(rx,p%o)*p%o
        rq = dot_product(rx,q)
        aq = norm(q)

        ! direction of tangential principal curvature
        if (aq>esmall) then
           ruq = rq/aq
        else
           ! cure the singularity for r almost parallel to o:
           ! choose an arbitrary direction normal to o
           ruq = 0.0_rk
        end if

        ! sagittal radius of curvature
        principal_curvature_radii(2) = &
             &sqrt((ro**2/rp)**2+(1.0_rk-ro**2/rp**2)*ruq**2)

        ! tangential radius of curvature
        principal_curvature_radii(1) = &
             &principal_curvature_radii(2)**3 * rp**2/ro**4
    end function principal_curvature_radii

    !> compute principal curvatures and their directions of a particle
    !> surface at a given point
    !>
    !> \param[in] k particle index in \c P
    !>
    !> \param[in] x point on the surface of particle \c k
    !>
    !> \param[in] rop half axes \c R_orth and \c R_para for \c P(p)
    !>
    !> \param[out] pp unit directions of principal curvatures arranged
    !> to form a right-handed coordinate system together with the
    !> outward pointing surface normal; \c pp(1,:) is the direction of
    !> tangential principal curvature, \c pp(2,:) the sagittal
    !> one. The angle between \c pp(1,:) and \c P(k)%o is 90 degrees
    !> or less.
    !>
    !> \param[out] c curvatures (inverse radii) corresponding to \c pp
    !>
    !> The computation of the tangent directions follow simple
    !> geometric considerations for an ellipsoid of revolution. The
    !> related curvatures and terminology follow Harris,
    !> Ophthal. Physiol. Opt. 26, 497 (2006) but could be computed
    !> from a local 2nd order expansion of the surface relatively
    !> easily as well.
    subroutine principal_curvatures(k,x,rop,pp,c)
        integer,intent(in) :: k
        real(kind=rk),intent(in) :: x(3),rop(2)
        real(kind=rk),intent(out) :: pp(2,3),c(2)

        real(kind=rk) :: aq,q(3),r(3),ro,rp,rq,ruq
        real(kind=rk),parameter :: ex(3)=(/1.0_rk,0.0_rk,0.0_rk/)
        real(kind=rk),parameter :: ey(3)=(/0.0_rk,1.0_rk,0.0_rk/)

        ro = rop(1)
        rp = rop(2)

        r = x-P(k)%x
        q = r - dot_product(r,P(k)%o)*P(k)%o
        rq = dot_product(r,q)
        aq = norm(q)

        ! direction of tangential principal curvature
        if (aq>esmall) then
           pp(1,:) = unit_vector(-q*dot_product(r,P(k)%o)*ro/rp&
                & + P(k)%o*rq*rp/ro)
           ruq = rq/aq
        else
           ! cure the singularity for r almost parallel to P(k)%o:
           ! choose an arbitrary direction normal to P(k)%o
           if (abs(dot_product(P(k)%o,ex)) < abs(dot_product(P(k)%o,ey))) then
              pp(1,:) = unit_vector(cross_product(P(k)%o,ex))
           else
              pp(1,:) = unit_vector(cross_product(P(k)%o,ey))
           end if
           ruq = 0.0_rk
        end if

        ! direction of sagittal principal curvature
        pp(2,:) = unit_vector(cross_product(r,pp(1,:)))

        ! sagittal curvature
        c(2) = 1.0_rk&
             &/sqrt((ro**2/rp)**2+(1.0_rk-ro**2/rp**2)*ruq**2)

        ! tangential curvature
        c(1) = c(2)**3*ro**4/rp**2
    end subroutine principal_curvatures

    !> try and retrieve a pair of surface points on two particles that
    !> was possibly stored by \c store_surface_points() before
    !>
    !> \param[in] i index of 1st particle in \c P
    !>
    !> \param[in] j index of 2nd particle in \c P
    !>
    !> \param[out] rxi point on particle \c i, measured from its center
    !>
    !> \param[out] rxj point on particle \c j, measured from its center
    !>
    !> \param[out] success \c .true. if a pair of points could be
    !> found, \c .false. otherwise
    subroutine retrieve_surface_points(i,j,rxi,rxj,success)
        integer,intent(in) :: i,j
        real(kind=rk),intent(out) :: rxi(3),rxj(3)
        logical,intent(out) :: success
        integer :: key
        real(kind=rk) :: values(6),rxib(3),rxjb(3)

        key = pair_uid(i,j)
        if (Mi6r_exist(surface_points,key)) then
           values = Mi6r_map(surface_points,key)

           ! position on particle with lower uid was stored first
           if (P(i)%uid<P(j)%uid) then
              rxib = values(1:3)
              rxjb = values(4:6)
           else
              rxjb = values(1:3)
              rxib = values(4:6)
           end if

           ! transform back to space-fixed frame
           rxi = matmul(transpose(space_to_body_matrix(P(i)%q)),rxib)
           rxj = matmul(transpose(space_to_body_matrix(P(j)%q)),rxjb)

           success = .true.
        else
           success = .false.
        end if
    end subroutine retrieve_surface_points

    !> store a pair of points on two particle surfaces for later use
    !>
    !> \param[in] i index of 1st particle in \c P
    !>
    !> \param[in] j index of 2nd particle in \c P
    !>
    !> \param[in] rxi point on particle \c i, measured from its center
    !>
    !> \param[in] rxj point on particle \c j, measured from its center
    subroutine store_surface_points(i,j,rxi,rxj)
        integer,intent(in) :: i,j
        real(kind=rk),intent(in) :: rxi(3),rxj(3)
        integer :: key,surface_points_lifetime
        real(kind=rk) :: values(6),rxib(3),rxjb(3)

        ! don't let surface_points get too big or retrieving the data
        ! will cost too much time
        if (Mi6r_count(surface_points)>=lubrication_cox_max_store_points) then
           call Mi6r_clear(surface_points)

           ! record minimum lifetime of data in surface_points
           surface_points_lifetime = nt-surface_points_reset_time
           if (surface_points_lifetime<5) then
              write (msgstr,fmt='("store_surface_points(): lifetime of stored '&
                   &//'surface points was only ",I0," time steps or less. If '&
                   &//'this happens regularly, peformance might be increased '&
                   &//'by increasing lubrication_cox_max_store_points in '&
                   &//'/md_fluid_ladd/")') surface_points_lifetime
              call log_msg_md(msgstr,forall=.true.)
           end if
           min_surface_points_lifetime = &
                &min(min_surface_points_lifetime, surface_points_lifetime)
           surface_points_reset_time = nt
        end if

        ! particles can rotate, so store positions in body-fixed frames
        rxib = matmul(space_to_body_matrix(P(i)%q),rxi)
        rxjb = matmul(space_to_body_matrix(P(j)%q),rxj)

        ! store position on particle with lower uid first
        if (P(i)%uid<P(j)%uid) then
           values = (/rxib(1),rxib(2),rxib(3),rxjb(1),rxjb(2),rxjb(3)/)
        else
           values = (/rxjb(1),rxjb(2),rxjb(3),rxib(1),rxib(2),rxib(3)/)
        end if

        key = pair_uid(i,j)
        call Mi6r_safe_assign(surface_points,key,values)
    end subroutine store_surface_points

    !> returns the outward directed unit normal vector for a point
    !> located on the surface of a spheroidal particle
    !>
    !> \param[in] p particle
    !>
    !> \param[in] rx position of the point measured from the spheroid
    !> center
    !>
    !> \param[in] rop half axes \c R_orth and \c R_para for \c p
    !>
    !> \returns outward directed unit normal vector
    pure function spheroid_unit_normal(p,rx,rop)
        real(kind=rk) :: spheroid_unit_normal(3)
        type(md_particle_type),intent(in) :: p
        real(kind=rk),intent(in) :: rx(3),rop(2)
        real(kind=rk) :: q(3),QQ,ro,rp,rxo

        ro = rop(1)
        rp = rop(2)

        ! computation is based on a projection of the plane containing
        ! p%o and rx, where the spheroid's axes are p%o and q
        rxo = dot_product(rx,p%o)
        q = rx - rxo*p%o

        ! if rx is on the spheroid's surface, in theory, QQ cannot
        ! exceed ro**2, to be sure, limit the argument of sqrt() below
        ! to positive numbers in order to avoid NaNs
        QQ = dot_product(q,q)

        spheroid_unit_normal = &
             &(rp*q &
             &+ sign(1.0_rk,rxo)*ro*sqrt(max(0.0_rk,ro**2-QQ))*p%o)&
             &/sqrt(QQ*(rp**2-ro**2) + ro**4)
    end function spheroid_unit_normal

    !> computes the translation and rotation velocity at the surface
    !> of a particle
    !>
    !> \param[in] k particle index in P
    !>
    !> \param[in] x position on the surface
    !>
    !> \param[out] v translation velocity
    !>
    !> \param[out] w rotation velocity
    subroutine surface_velocity(k,x,v,w)
        integer,intent(in) :: k
        real(kind=rk),intent(in) :: x(3)
        real(kind=rk),intent(out) :: v(3),w(3)

        v = P(k)%v + cross_product(P(k)%ws,x-P(k)%x)
        w = P(k)%ws
    end subroutine surface_velocity

    !> compute cubic coefficients in a rotated coordinate system
    !>
    !> \param[in] in input set of cubic coefficients
    !>
    !> \param[in] s sine of the rotation angle
    !>
    !> \param[in] c cosine of the rotation angle
    !>
    !> \param[out] out cubic coefficients in rotated system
    subroutine transform_cubic_coefficients(in,s,c,out)
        real(kind=rk),intent(in) :: in(0:3),s,c
        real(kind=rk),intent(out) :: out(0:3)

        ! this can be derived easily by substituting for the original
        ! coordinates in the cubic term expressions of the rotated
        ! ones and sorting the result for powers of the rotated
        ! coordinates
        out(0) = in(0)*c**3 - in(1)*s*c**2 + in(2)*s**2*c - in(3)*s**3
        out(1) = in(0)*3.0_rk*s*c**2 + in(1)*(c**3-2.0_rk*s**2*c) &
             &+ in(2)*(s**3-2.0_rk*s*c**2) + in(3)*3.0_rk*s**2*c
        out(2) = in(0)*3.0_rk*c*s**2 + in(1)*(2.0_rk*c**2*s-s**3) &
             &+ in(2)*(c**3-2.0_rk*c*s**2) - in(3)*(3.0_rk*c**2*s)
        out(3) = in(0)*s**3 + in(1)*c*s**2 + in(2)*c**2*s + in(3)*c**3
    end subroutine transform_cubic_coefficients

#endif
end module lbe_md_fluid_ladd_club_module
