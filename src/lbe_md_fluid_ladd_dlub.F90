#include "lbe.h"

!> per-link lubrication correction for \c interaction=='ladd'
!> particles, very similar to Ding&Aidun. J. Stat. Phys. 112, 685-708
!> (2003) or also MacMeccan et al. JFM 618, 13-39 (2009).
module lbe_md_fluid_ladd_dlub_module
#ifdef MD
#ifdef LADD_DLUB
    use lbe_globals_module, only: bounce,c_dir,c_length,g,nnonrest,rk
    use lbe_helper_module, only: cross_product
    use lbe_md_boundary_condition_module, only: min_distance
    use lbe_md_fluid_ladd_mc_module, only: pfr
    use lbe_md_fluid_ladd_parms_module, only: R_orth,R_para
    use lbe_md_globals_module
    use lbe_md_helper_module, only: error_md
    use lbe_parallel_module, only: check_allocate,start
    use lbe_parms_module, only: tau_r
    use lbe_timer_module, only: register_timer
    use map_module, only: Mii_map
#ifdef LADD_DLUB_CENTRAL
    use lbe_helper_module, only: unit_vector
#endif

    implicit none
    private

    public dlub_clear_links,dlub_force_and_torque,dlub_setup,dlub_store_link

    !> enable link-wise lubrication correction similar to Ding&Aidun, 2003
    logical,save,public :: ding_lubrication=.false.

    !> for the sake of stability, linkwise lubrication corrections are
    !> clipped at the value they would achieve at a gap width of \c
    !> clip_frac times the length of the respective lattice
    !> link. Forces still can diverge if the velocities diverge.
    real(kind=rk),save,public :: clip_frac=0.01_rk

    !> For numerical reasons, particle overlap has to be avoided
    !> manually. Therefore, for every link-projected gap width that
    !> (relatively) is smaller than \c clip_frac , a repulsive spring
    !> force with a stiffness of \c contact_stiffness starts to act.
    real(kind=rk),save,public :: contact_stiffness=1000.0_rk

    !> friction coefficient for per-link particle-particle friction
    real(kind=rk),save,public :: contact_friction=0.0_rk

    !> timer for \c dlub_force_and_torque()
    integer,save,public :: ti_md_fluid_dlubft

    !> information about a particle-particle link
    type link_type
       integer :: uid1          !< uid of particle 1
       integer :: uid2          !< uid of particle 2
       integer :: x1(3) !< global lattice position of site 1 (halo is allowed)
       integer :: x2(3) !< global lattice position of site 2 (halo is allowed)
       integer :: dir !< index of direction from 1 to 2 (\c 1..nnonrest)
       integer :: length !< number of links between \c x1 and \c x2 (\c 1..2)
    end type link_type

    !> collection of all particle-particle links which might require
    !> the calculation of lubrication forces
    type(link_type),save,allocatable :: links(:)

    !> number of links currently stored in \c links
    integer,save :: n_links

    !> kinematic viscosity of red component
    real(kind=rk),save :: nu_r

    !> inverse particle half-axes
    real(kind=rk),save :: inv_R(3)

contains

    !> clear all entries from \c links
    subroutine dlub_clear_links
        n_links = 0
    end subroutine dlub_clear_links

    !> add forces and torques due to linkwise lubrication correction
    !> according to the content of \c links to all local and halo'd
    !> particles
    subroutine dlub_force_and_torque
        integer :: i

        do i=1,n_links
           call link_force_and_torque(links(i))
        end do
    end subroutine dlub_force_and_torque

    !> setup stage for \c lbe_md_fluid_ladd_dlub_module
    subroutine dlub_setup
        integer stat

        ! kinematic viscosity
        nu_r = (2.0_rk*tau_r-1.0_rk)/6.0_rk

        ! inverse particle half-axes
        inv_R = 1.0_rk/(/R_orth,R_orth,R_para/)

        ! this is increased later on demand
        allocate (links(1),stat=stat)
        call check_allocate(stat,'dlub_setup: links')

        call register_timer('MD:Fluid:Dlub_ft',ti_md_fluid_dlubft)

        if (polydispersity) call error_md("LADD_DLUB does not "&
             &//"support polydisperse particles yet---disable polydispersity "&
             &//"or compile without LADD_DLUB!")

#ifdef LADD_DLUB_CENTRAL
        if (R_orth/=R_para) call error_md('LADD_DLUB_CENTRAL does not make '&
             &//'sense for non-spherical particles---set R_orth==R_para or '&
             &//'compile without LADD_DLUB_CENTRAL')
#endif
    end subroutine dlub_setup

    !> stores a particle-particle link for later calculation of
    !> lubrication forces
    !>
    !> \param[in] uid1 uid of particle 1
    !>
    !> \param[in] uid2 uid of particle 2
    !>
    !> \param[in] lx1 local lattice position of site 1
    !>
    !> \param[in] lx2 local lattice position of site 2
    !>
    !> \param[in] dir index of direction from 1 to 2 (\c 1..nnonrest)
    !>
    !> \param[in] length number of links between \c x1 and \c x2 (\c 1..2)
    subroutine dlub_store_link(uid1,uid2,lx1,lx2,dir,length)
        integer,intent(in) :: uid1,uid2,lx1(3),lx2(3),dir,length

        if (n_links>=size(links)) call boost_links(2*n_links)
        n_links = n_links+1

        links(n_links)%uid1 = uid1
        links(n_links)%uid2 = uid2
        links(n_links)%x1 = lx1+start-1
        links(n_links)%x2 = lx2+start-1
        links(n_links)%dir = dir
        links(n_links)%length = length
    end subroutine dlub_store_link

    !> enlarges the capacity of \c links
    !>
    !> \param[in] new_size desired new capacity
    subroutine boost_links(new_size)
        integer,intent(in) :: new_size
        type(link_type),allocatable :: tmp(:)
        integer :: stat

        allocate (tmp(n_links),stat=stat)
        call check_allocate(stat&
             &,'lbe_md_fluid_ladd_dlub_module: boost_links(): tmp')

        tmp = links(1:n_links)
        deallocate (links)

        allocate (links(new_size),stat=stat)
        call check_allocate(stat&
             &,'lbe_md_fluid_ladd_dlub_module: boost_links(): links')

        links(1:n_links) = tmp
        deallocate (tmp)
    end subroutine boost_links

    !> calculates per-link particle-particle lubrication correction force
    !>
    !> \param[in] l link between two particles
    !>
    !> The resulting per-link contribution to the lubrication
    !> correction in case of actual particle surface displacements of
    !> less than one lattice vector is applied as force and torque
    !> increment for both particles.
    subroutine link_force_and_torque(l)
        type(link_type),intent(in) :: l
        ! suggested by Ding&Aidun (2003): qbar=0.6. This is 3*qbar.
        real(kind=rk),parameter :: qbar_3=0.6_rk*3.0_rk
        real(kind=rk) :: clip_length,f(3),f_lub, f_rep,gap,lambda,lsi,lsj&
             &,rhof_r,rsi(3),rsj(3),rxi(3),rxj(3),U,U_t(3),vi(3),vj(3)
        integer :: i,j
#ifdef LADD_DLUB_CENTRAL
        real(kind=rk) :: rij_norm(3)
#endif

        ! find particles in P array
        i = Mii_map(uid2i,l%uid1)
        j = Mii_map(uid2i,l%uid2)

!!$        if (P(i)%uid/=l%uid1.or.P(j)%uid/=l%uid2) then
!!$           print *,'l%uid1=',l%uid1
!!$           print *,'l%uid2=',l%uid2
!!$           print *,'P(i)%uid=',P(i)%uid
!!$           print *,'P(j)%uid=',P(j)%uid
!!$           call error_md('link_force_and_torque(): particle not present')
!!$        end if

        ! vector from particle center to respective boundary node
        ! (with minimum image criterion)
        rxi=min_distance(real(l%x1,kind=rk),P(i)%x)
        rxj=min_distance(real(l%x2,kind=rk),P(j)%x)

        ! gap width projected on link direction
        lsi = link_surface_intersection(rxi,P(i)%q,bounce(l%dir))
        lsj = link_surface_intersection(rxj,P(j)%q,l%dir)
        gap = real(l%length,kind=rk)*c_length(l%dir)-(lsi+lsj)

        cutoff: if (gap<c_length(l%dir)) then
           ! never apply a force larger than that which would result
           ! from a gap of size clip_frac times the respective lattice
           ! direction. If the gap will be clipped, apply repulsive
           ! force first.
           clip_length = clip_frac*c_length(l%dir)

           ! letting the repulsive force set in already at
           ! c*clip_length with c>1 would enhance stability but also
           ! introduce one more parameter and is therefore tried to be
           ! avoided.
           if (gap<clip_length) then
              ! The repulsive force grows linearly when going down to
              ! zero gap width where it stays constant, so the minimum
              ! gap width effective for f_rep calculation is 0.
              f_rep = contact_stiffness*(clip_length-max(0.0_rk,gap))
           else
              f_rep = 0.0_rk
           end if
           gap = max(gap,clip_length)

           ! local fluid density
#ifdef LADD_SURR_RHOF
           rhof_r = 0.5_rk*(P(i)%rhof+P(j)%rhof)
#else
           rhof_r = pfr
#endif

           ! local velocity difference between both
           ! link-surface-intersections projected on the link
           rsi = rxi-lsi*c_dir(:,l%dir)
           rsj = rxj+lsj*c_dir(:,l%dir)
           vi = P(i)%v + cross_product(P(i)%ws,rsi)
           vj = P(j)%v + cross_product(P(j)%ws,rsj)
           U = dot_product(c_dir(:,l%dir),(vj-vi))

           ! tangential part of the velocity differenc, needed for friction
           U_t = (vj-vi) - U*c_dir(:,l%dir)

           ! local curvature. Not so sure actually whether mean or
           ! Gaussian curvature should be taken---assuming mean curvature
           ! to be correct because the Gaussian one would be zero for a
           ! cylinder which looks weird (cylinders might be weird,
           ! though)
           lambda = 0.5_rk*&
                &(local_mean_curvature(P(i)%o,rsi)&
                &+local_mean_curvature(P(j)%o,rsj))

           ! compute lubrication force
           f_lub = qbar_3*nu_r*rhof_r*U&
                &*(1.0_rk/gap**2-1.0_rk/c_length(l%dir)**2)&
                &/(2.0_rk*c_length(l%dir)**2*lambda)

           ! sum forces on this link
           f = c_dir(:,l%dir)*(f_lub+f_rep) + contact_friction*U_t*f_rep
#ifdef LADD_DLUB_CENTRAL
           rij_norm = unit_vector(min_distance(P(i)%x,P(j)%x))
           f = dot_product(f,rij_norm)*rij_norm
#endif

           ! apply force and torques on both particles
           P(i)%f = P(i)%f + f
#ifndef LADD_DLUB_CENTRAL
           P(i)%t = P(i)%t + cross_product(rsi,f)
#endif

           P(j)%f = P(j)%f - f
#ifndef LADD_DLUB_CENTRAL
           P(j)%t = P(j)%t - cross_product(rsj,f)
#endif
        end if cutoff
    end subroutine link_force_and_torque

    !> returns the distance from a boundary node to the particle's
    !> surface, measured along a link
    !>
    !> \param[in] rx lattice position of the boundary node measured
    !> from the particle's center position
    !>
    !> \param[in] q quaternion describing the particle's orientation
    !>
    !> \param[in] s distribution index of the outward directed link
    !> starting from \c rx
    !>
    !> \returns \c w with \c rx+w*u(s) hitting exactly the theoretical
    !> particle surface where \c u(s) is the unit vector along lattice
    !> direction \c s
    function link_surface_intersection(rx,q,s)
        real(kind=rk) :: link_surface_intersection
        real(kind=rk),intent(in) :: rx(3),q(0:3)
        integer,intent(in) :: s
        real(kind=rk) :: A(3,3),lb(3),lbs(3),PP,iQQ,rb(3),rbs(3)

        ! lattice vector and boundary node in body-fixed coordinates
        A = reshape((/&
             &q(0)**2+q(1)**2-q(2)**2-q(3)**2,&
             &2.0_rk*(q(1)*q(2)-q(0)*q(3)),&
             &2.0_rk*(q(1)*q(3)+q(0)*q(2)),&

             &2.0_rk*(q(1)*q(2)+q(0)*q(3)),&
             &q(0)**2-q(1)**2+q(2)**2-q(3)**2,&
             &2.0_rk*(q(2)*q(3)-q(0)*q(1)),&

             &2.0_rk*(q(1)*q(3)-q(0)*q(2)),&
             &2.0_rk*(q(2)*q(3)+q(0)*q(1)),&
             &q(0)**2-q(1)**2-q(2)**2+q(3)**2&
             &/),&
             &(/3,3/))
        lb = matmul(A,c_dir(:,s))
        rb = matmul(A,rx)

        ! rescale all lengths to the case of a unit sphere for the
        ! sake of convenience and performance
        lbs = lb*inv_R
        rbs = rb*inv_R

        ! avoid duplicate calculations and code
        iQQ = 1.0_rk/dot_product(lbs,lbs)
        PP = dot_product(lbs,rbs)*iQQ

        link_surface_intersection &
             &= -PP+sqrt(PP*PP-(dot_product(rbs,rbs)-1.0_rk)*iQQ)
!!$    ! This would be valid for spheres only...
!!$    rb = rx
!!$    R = R_orth
!!$    link_surface_intersection = -dot_product(c_dir(:,s),rb)&
!!$         &+sqrt(dot_product(c_dir(:,s),rb)**2-dot_product(rb,rb)+R**2)
    end function link_surface_intersection

    !> calculates the local mean curvature for a point on the theoretical
    !> surface of a particle
    !>
    !> \param[in] o orientation vector (directed along the axis of
    !> rotational symmetry of the ellipsoid)
    !>
    !> \param[in] rs position on the particle surface to calculate the
    !> curvature for
    !>
    !> \returns mean curvature (source: Wikipedia)
    !>
    !> The center of the ellipsoid is assumed to be at the origin.
    function local_mean_curvature(o,rs)
        real(kind=rk) :: local_mean_curvature
        real(kind=rk),intent(in) :: o(3),rs(3)
        real(kind=rk) :: B,rho(3),u

        u = dot_product(o,rs)
        rho = rs-u*o
        B = R_orth**2&
             &+(R_para**2-R_orth**2)*dot_product(rho,rho)/dot_product(rs,rs)

        local_mean_curvature = R_para*(R_orth**2+B)/(2.0_rk*R_orth*B*sqrt(B))
    end function local_mean_curvature
#endif
#endif
end module lbe_md_fluid_ladd_dlub_module
