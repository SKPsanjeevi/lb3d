#include "lbe.h"

!> Calculation of electric potential and fields.
module lbe_elec_poisson_solver_module

#ifdef ELEC
  use lbe_elec_globals_module
  use lbe_elec_helper_module
  use lbe_elec_parallel_module, only: phi_halo
  use lbe_elec_timer_module
  use lbe_parallel_module
  use mpi

#ifdef P3M
  use lammps_global_module
#endif

  implicit none

  private

  public solve_poisson, calc_E_fd
#ifdef P3M
  public E_prev
#endif

#ifdef P3M
  real(kind=rk), dimension(:,:,:,:), allocatable :: E_prev
#endif

contains

!> Depending on the value of poisson_solver in the input file, the
!> (generalized) Poisson equation will be solved as decided by this
!> wrapper function.
subroutine solve_poisson(N)
  implicit none

  type(lbe_site), intent(inout) :: N(1-halo_extent:, 1-halo_extent:, 1-halo_extent:)

  select case( poisson_solver_id )
    case( poisson_solver_SOR )
      call solve_poisson_sor(N)
    case( poisson_solver_p3m )
#ifdef P3M
      call solve_poisson_p3m(N)
      ! P3M needs an explicit halo exchange done at the end.
      call halo_exchange(N, phi_halo)
#endif
  end select

end subroutine solve_poisson

!> The SOR scheme can be used to solve the Poisson equation.
!> Input: charge distribution (rho_p, rho_m (and eps) in N).
!> Output: potential (phi in N).
subroutine solve_poisson_sor(N)
  implicit none

  type(lbe_site), intent(inout) :: N(1-halo_extent:, 1-halo_extent:, 1-halo_extent:)

  real(kind=rk) :: omega, depsphi, epstot, epsh, residual, absres
  real(kind=rk) :: tol_abs, tol_rel
  integer :: it, ipass
  integer :: isw, jsw, ksw
  integer :: i, j, k, ip, jp, kp
  integer :: veldir
  integer :: mpierror
  integer :: oddity, moddity

  real(kind=rk) :: cur_eps, cur_phi, cur_rho, nb_phi

  real(kind=rk), dimension(2) :: rnorm, rnorm_local
  real(kind=rk) :: maxres, maxres_local = 0.0

  call start_timer(ti_elec_sor)

  tol_rel = acc_SOR
  tol_abs = acc_SOR

  rnorm_local(1) = 0.0_rk

  ! We might need phi three sites away from the physical region, but we can just calculate the physical region and
  ! exchange the halo later.
  do i = 1, nx
    do j = 1, ny
      do k = 1, nz
        depsphi = 0.0_rk
        cur_rho = N(i,j,k)%rho_p - N(i,j,k)%rho_m
        cur_phi = N(i,j,k)%phi
#ifdef ELEC_NNONREST
        if ( local_eps ) then
          cur_eps = get_local_eps(N,i,j,k)
          do veldir = 1,nnonrest
            nb_phi = get_nb_phi(N, veldir, i, j, k, ip, jp, kp)
            depsphi = depsphi + g(veldir) * ( get_local_eps(N,ip,jp,kp) + cur_eps ) * ( nb_phi - cur_phi )
          end do
          depsphi = depsphi / 12.0_rk
        else
          do veldir = 1,nnonrest
            nb_phi = get_nb_phi(N, veldir, i, j, k, ip, jp, kp)
            depsphi = depsphi + g(veldir) * nb_phi
          end do
          depsphi = eps_uniform * ( depsphi - 24.0_rk * cur_phi ) / 6.0_rk
        end if
        rnorm_local(1) = rnorm_local(1) + abs(depsphi + cur_rho )
#else
        if ( local_eps ) then
          cur_eps = get_local_eps(N,i,j,k)
          do veldir = 1,nnn
            nb_phi = get_nb_phi(N,veldir, i, j, k, ip, jp, kp)
            depsphi = depsphi + 0.5_rk * ( get_local_eps(N,ip,jp,kp) + cur_eps ) * ( nb_phi - cur_phi )
          end do
        else
          do veldir = 1,nnn
            nb_phi = get_nb_phi(N,veldir, i, j, k, ip, jp, kp)
            depsphi = depsphi + nb_phi
          end do
          depsphi = eps_uniform * ( depsphi - 6.0_rk * cur_phi )
        end if
        rnorm_local(1) = rnorm_local(1) + abs(depsphi + cur_rho )
#endif
      end do
    end do
  end do

  ! See if we need to invert our checkerboard pattern due to parallelization

  oddity = ccoords(1) * nx + ccoords(2) * ny + ccoords(3) * nz
  moddity = mod(oddity,2)

  ! Iterate to solution
  omega = 1.0_rk

  do it = 1, maxits_SOR

    ! Compute current normal of the residual

    rnorm_local(2) = 0.0_rk
    maxres_local   = 0.0_rk

    ksw = 1

    do ipass = 1, 2
      jsw = ksw
      do k = 1, nz
        isw = jsw
        ! Checkerboard inversion?
        if (moddity .eq. 1) then
          isw = 3 - isw
        end if
        do j = 1, ny
          do i = isw, nx, 2
            depsphi = 0.0_rk
            cur_phi = N(i,j,k)%phi
            cur_rho = N(i,j,k)%rho_p - N(i,j,k)%rho_m
#ifdef ELEC_NNONREST
            if ( local_eps ) then
              epstot = 0.0_rk
              cur_eps = get_local_eps(N,i,j,k)
              do veldir = 1,nnonrest ! FIX
                nb_phi = get_nb_phi(N, veldir, i, j, k, ip, jp, kp)
                epsh = g(veldir) * ( get_local_eps(N,ip,jp,kp) + cur_eps ) / 6.0_rk
                epstot = epstot + epsh ! sum of weighted averages
                depsphi = depsphi + epsh * 0.5_rk * ( nb_phi - cur_phi )
              end do
              residual = depsphi + cur_rho
              N(i,j,k)%phi = N(i,j,k)%phi + ( omega * residual / ( 0.75_rk * epstot ) )
            else
              do veldir = 1,nnonrest
                nb_phi = get_nb_phi(N, veldir, i, j, k, ip, jp, kp)
                depsphi = depsphi + g(veldir) * ( nb_phi - cur_phi )
              end do
              residual = eps_uniform * depsphi / 6.0_rk + cur_rho
              N(i,j,k)%phi = N(i,j,k)%phi + ( omega * residual / ( 6.0_rk * eps_uniform ) )
            end if
            absres = abs(residual)
            rnorm_local(2) = rnorm_local(2) + absres
            if (absres > maxres_local) then
               maxres_local = absres
            end if
#else
            if ( local_eps ) then
              epstot = 0.0_rk
              cur_eps = get_local_eps(N,i,j,k)
              do veldir = 1,nnn
                nb_phi = get_nb_phi(N, veldir, i, j, k, ip, jp, kp)
                epsh = 0.5_rk * ( get_local_eps(N,ip,jp,kp) + cur_eps )
                epstot = epstot + epsh
                depsphi  = depsphi  + epsh * ( nb_phi - cur_phi )
              end do
              residual = depsphi + cur_rho
              N(i,j,k)%phi = N(i,j,k)%phi - ( omega * residual / ( -1.0_rk * epstot ) )
            else
              do veldir = 1,nnn
                nb_phi = get_nb_phi(N, veldir, i, j, k, ip, jp, kp)
                depsphi = depsphi + nb_phi
              end do
              residual = eps_uniform * ( depsphi - 6.0_rk * cur_phi ) + cur_rho
              N(i,j,k)%phi = N(i,j,k)%phi - ( omega * residual / ( -6.0_rk * eps_uniform ) )
            end if
            absres = abs(residual)
            rnorm_local(2) = rnorm_local(2) + absres
            if (absres > maxres_local) then
               maxres_local = absres
            end if
#endif
          end do ! end i loop
          isw = 3 - isw
        end do ! end j loop
        jsw = 3 - jsw
      end do ! end k loop
      ksw = 3 - ksw

      omega = 1.0_rk / ( 1.0_rk - 0.25_rk * radius_SOR * radius_SOR * omega )

      call start_timer(ti_elec_sor_phi_halo)
      call halo_exchange(N, phi_halo)
      call stop_timer(ti_elec_sor_phi_halo)

    end do ! end ipass loop

    if ( mod(it-1, n_check_sor) == 0 ) then
      call start_timer(ti_elec_sor_allreduce)
      call MPI_Allreduce(rnorm_local, rnorm, 2, LBE_REAL, MPI_SUM, Comm_cart, mpierror)
      call MPI_Allreduce(maxres_local, maxres, 1, LBE_REAL, MPI_MAX, Comm_cart, mpierror)
      call stop_timer(ti_elec_sor_allreduce)

      if (rnorm(2) .lt. acc_SOR .or. maxres .lt. meps) then
        if ( n_show_SOR > 0) then
          write(msgstr,"('  Finished  SOR iteration = ',I6,' at t = ',I0,' , maxres = ',(E16.8,X),', rnorm = ',2(E16.8,X),' , tol = ',2(E16.8,X))") it, nt, maxres, rnorm, tol_abs, tol_rel
          call log_msg_elec(msgstr)
        end if

        call stop_timer(ti_elec_sor)
        return
      else
        if ( n_show_SOR > 0 ) then
          if ( mod(it-1, n_show_SOR) == 0 ) then
            write(msgstr,"('  Performed SOR iteration = ',I6,' at t = ',I0,' , maxres = ',(E16.8,X),' , rnorm = ',2(E16.8,X),' , tol = ',2(E16.8,X))") it, nt, maxres, rnorm, tol_abs, tol_rel
            call log_msg_elec(msgstr)
          end if
        end if
      end if
    end if

  end do ! end iteration loop

  call stop_timer(ti_elec_sor)

  call error_elec("  SOR won't converge!")

end subroutine solve_poisson_sor

#ifdef P3M

!> The Poisson solver from the implementation of the P3M method from LAMMPS can 
!> be used to solve the Poisson equation.
!> Input: (same as SOR) charge distribution (rho_p, rho_m (and eps) in N).
!> Output: (same as SOR) potential (phi in N).
subroutine solve_poisson_p3m(N)
  implicit none

  type(lbe_site), intent(inout) :: N(1-halo_extent:, 1-halo_extent:, 1-halo_extent:)

  logical, save :: first_call = .true.

  real(kind=rk) :: E_grad_eps, grad_eps(3)
  integer :: veldir, ip, jp, kp, iflag, icount, i, j, k, ipos, jpos, kpos, mpierror
  logical :: phi_and_E

  call start_timer(ti_elec_p3m)

  call start_timer(ti_elec_p3m_setup)
  !gjp call to map invariant data
  if ( first_call ) then 
    first_call = .false.

    !gjp maps lb3dv6 parameters and packs into lammps parameters
    !gjp LAMMPS counts from zeros, in terms of both real space and arrays, etc.
    !gjp This version employs only the poisson solver from the
    !gjp P3M algorithm, thus we do not require rho/phi mesh ghost cells
    !gjp and certain subroutines.  These routines are retained in this
    !gjp release for completness.

    call log_msg_elec("  Setting up P3M solver (first time only) ...")

    !gjp map basic MPI stuff
    node = myrankc ! rank
    ncores = nprocs ! number of cores
    !gjp extent of virtual topolgy 
    nx_pppm_input = tnx ! x extent on virtual topology
    ny_pppm_input = tny ! y extent on virtual topology
    nz_pppm_input = tnz ! z extent on virtual topology

    !gjp proc's coordinate in virtual topology
    me(1) = ccoords(1)
    me(2) = ccoords(2)
    me(3) = ccoords(3)
    !gjp extent of virtual topogology
    pgrid(1) = cdims(1)
    pgrid(2) = cdims(2)
    pgrid(3) = cdims(3)
    !gjp use lb3d neighbours and not employ lammps grid_3d routine as this breaks code
    mpart(1,3)=nnprocs(3,1)
    mpart(2,3)=nnprocs(3,2)
    mpart(1,2)=nnprocs(2,1)
    mpart(2,2)=nnprocs(2,2)
    mpart(1,1)=nnprocs(1,1)
    mpart(2,1)=nnprocs(1,2)

    !gjp box is lammps' computational box, with corner at origin
    !gjp when using poisson solver alone (and not full P3M code),
    !gjp box encompoases the FFT mesh.
    !gjp NB box is expanded by 0.5 in all directions so that the
    !gjp extent of the compuational box matches extent of FFT mesh
    box(1,1)=0.0-0.5
    box(2,1)=real(tnx)-1.0+0.5
    box(1,2)=0.0-0.5
    box(2,2)=real(tny)-1.0+0.5
    box(1,3)=0.0-0.5
    box(2,3)=real(tnz)-1.0+0.5
    xprd = box(2,1) - box(1,1) ! extend of box in x-direction
    yprd = box(2,2) - box(1,2) ! extend of box in y-direction
    zprd = box(2,3) - box(1,3) ! extend of box in z-direction

    !gjp sets 3D distribution, rather than 2D
    slabflag = 0
    slab_volfactor = 1.0
    idimension=3

    !gjp sets lammps border variable
    border(1,1) = float(me(1))/pgrid(1) * xprd + box(1,1)
    border(2,1) = float(me(1)+1)/pgrid(1) * xprd + box(1,1)
    border(1,2) = float(me(2))/pgrid(2) * yprd + box(1,2)
    border(2,2) = float(me(2)+1)/pgrid(2) * yprd + box(1,2)
    border(1,3) = float(me(3))/pgrid(3) * zprd + box(1,3)
    border(2,3) = float(me(3)+1)/pgrid(3) * zprd + box(1,3)

    !gjp orderflag = order of PPPM, how far into grid the charge overlaps
    !gjp orderflag = 0 breaks pppm_coeff
    orderflag = 1
    long_prec = .0001
    cutcoul = 10.0
    natoms = 1  ! makes no difference to computation of poisson solver alone
    qsqsum = natoms
    skin=0.01
    !gjp meshflag = 1 if user sets PPPM mesh, 0 otherwise code sets PPPM mesh
    meshflag = 1

    !gjp set ghost_cells=.false. as we are using poisson solver only
    !gjp   which does not require ghost cells
    ghost_cells=.false.

    !gjp set up pppm routine
    !gjp 0 arg sets iflag=0 inside pppm_coeff routine
    call pppm_coeff(0)
  end if

  if ( ghost_cells ) then
    !gjp pack charge into rho/phi-mesh
    ! --------------------------------------
    icount=0
    do k=1,nz
      do j=1,ny
        do i=1,nx
          icount = icount + 1
          if ( local_eps ) then
            call error_elec("Local epsilon not implemented for P3M with ghost cells as Poisson solver.")
          else
            density_brick(icount) = ( N(i,j,k)%rho_p - N(i,j,k)%rho_m )
          end if
        end do
      end do
    end do

    if ( .not. local_eps ) then
      density_brick = density_brick / eps_uniform
    end if
    !gjp brick2fft first updates density P3M mesh ghost cells then packs FFT mesh, 
    !gjp      and then remaps remap fft mesh from rho/phi-mesh distribution to 
    !gjp      optimial fft-mesh distribution
    call brick2fft(density_brick,density_fft)
  else ! no ghost cells means we can avoid brick2fft
    !gjp packing charge direclty into density_fft mesh bypasses some comp and comms
    !gjp only possible when solving poisson equ alone and not as part of p3m
    icount = 0
    do k=1,nz
      do j=1,ny
        do i=1,nx
          icount = icount + 1
          if ( local_eps ) then
            ! Calculate contribution from local dieletric properties.
#ifdef ELEC_NNONREST
            call error_elec("ELEC_NNONREST not implemented for grad eps calculation.")
#else
            grad_eps(:) = 0.0_rk
            do veldir = 1, nnn
              ip = i + cx(veldir)
              jp = j + cy(veldir)
              kp = k + cz(veldir)
              grad_eps(1) = grad_eps(1) + ( cx(veldir) * get_local_eps(N,ip,jp,kp) )
              grad_eps(2) = grad_eps(2) + ( cy(veldir) * get_local_eps(N,ip,jp,kp) )
              grad_eps(3) = grad_eps(3) + ( cz(veldir) * get_local_eps(N,ip,jp,kp) )
            end do
            grad_eps(:) = 0.5_rk * grad_eps(:)
#endif
            E_grad_eps = grad_eps(1) * ( 2.0_rk * N(i,j,k)%E(1) - E_prev(i,j,k,1) ) &
                       + grad_eps(2) * ( 2.0_rk * N(i,j,k)%E(2) - E_prev(i,j,k,2) ) &
                       + grad_eps(3) * ( 2.0_rk * N(i,j,k)%E(3) - E_prev(i,j,k,3) )

            ! if ( i .eq. 1 .and. j .eq. 1 .and. k .eq. 1 ) then
            !   write(msgstr,"('E = ',3(E16.8,X),' E_prev = ',3(E16.8,X), ' grad_eps = ',3(E16.8,X), ' E_grad_eps = ',E16.8)") N(i,j,k)%E(:), E_prev(i,j,k,:), grad_eps(:), E_grad_eps
            !   call log_msg_elec(msgstr)
            ! end if

            density_fft(icount) = ( ( N(i,j,k)%rho_p - N(i,j,k)%rho_m ) / get_local_eps(N, i, j, k) ) - E_grad_eps
          else
            ! Division by epsilon will be done outside the loop.
            density_fft(icount) = ( N(i,j,k)%rho_p - N(i,j,k)%rho_m )
          end if
        end do
      end do
    end do

    ! If we don't have local epsilon, do this general division outside the loop.
    if ( .not. local_eps ) then
      density_fft = density_fft / eps_uniform
    end if

  end if

  !gjp remap fft mesh from rho/phi-mesh distribution to optimial fft-mesh distribution
  call remap_3d(density_fft,density_fft,workvec1,plan_remap)
  call stop_timer(ti_elec_p3m_setup)

  !gjp iflag determins if lammps' 'new' or 'old' method employed
  !gjp iflag set to zero for safety
  iflag = 0

  !gjp original LAMMPS poisson code calculates E and not phi
  if ( E_solver_id == E_solver_p3m ) then
    if ( local_eps ) then
      call store_E_prev(N)
    end if

    phi_and_E = .true. ! compute both phi and E field using FFTs
  else
    phi_and_E = .false. ! compute only phi using FFTs
  end if

  call start_timer(ti_elec_p3m_poisson)
  call poisson(phi_brick,vdx_brick,vdy_brick,vdz_brick,iflag,phi_and_E)
  call stop_timer(ti_elec_p3m_poisson)

  call start_timer(ti_elec_p3m_unpack)
  !gjp map fft-mesh back to density/phi-mesh
  call fillbrick_phi(phi_brick)
  !gjp  assumes no ghost cells for phi_brick
  icount=0
  do k=1,nz
    do j=1,ny
      do i=1,nx
        icount = icount + 1
        N(i,j,k)%phi = phi_brick(icount)
      end do
    end do
  end do
  if ( E_solver_id == E_solver_p3m ) then
    call start_timer(ti_elec_p3m_calc_E)
    !gjp map fft-mesh back to density/phi-mesh
    call fillbrick(vdx_brick,vdy_brick,vdz_brick)
    !gjp assumes no ghost cells for vdx_brick
    icount=0
    do k=1,nz
      do j=1,ny
        do i=1,nx
          icount=icount+1
          N(i,j,k)%E(1) = -vdx_brick(icount)
          N(i,j,k)%E(2) = -vdy_brick(icount)
          N(i,j,k)%E(3) = -vdz_brick(icount)
        end do
      end do
    end do
    call stop_timer(ti_elec_p3m_calc_E)
  end if
  call stop_timer(ti_elec_p3m_unpack)

  call stop_timer(ti_elec_p3m)

end subroutine solve_poisson_p3m

!> When using P3M with local_eps = .true. , store the old electric field in a buffer.
subroutine store_E_prev(N)
  implicit none

  type(lbe_site), intent(inout) :: N(1-halo_extent:, 1-halo_extent:, 1-halo_extent:)

  integer :: i, j, k

  ! For charge advection, the required force at halo depth one in turn requires the electric field
  ! one more site into the halo. As the electric field E is not exchanged as part of the elec_halo
  ! we explicitly loop over the physical domain + 2.
  do k=-1,nz+2
    do j=-1,ny+2
      do i=-1,nx+2
        E_prev(i,j,k,1) = N(i,j,k)%E(1)
        E_prev(i,j,k,2) = N(i,j,k)%E(2)
        E_prev(i,j,k,3) = N(i,j,k)%E(3)
      end do    
    end do
  end do

end subroutine store_E_prev

#endif

!> Calculate the electric field E(:) using the potential phi as stored on the lattice sites, using finite difference.
subroutine calc_E_fd(N)
  implicit none

  type(lbe_site), intent(inout) :: N(1-halo_extent:, 1-halo_extent:, 1-halo_extent:)

  integer :: i, j, k
  integer :: veldir, ip, jp, kp
  real(kind=rk) :: E_ext(3)
  real(kind=rk) :: nb_phi

#ifdef P3M
  if ( poisson_solver_id == poisson_solver_p3m .and. local_eps ) then
    call store_E_prev(N)
  end if
#endif

  call start_timer(ti_elec_calc_E_fd)

  E_ext = (/ Ex, Ey, Ez /)

  ! For charge advection, the required force at halo depth one in turn requires the electric field
  ! one more site into the halo. As the electric field E is not exchanged as part of the elec_halo
  ! we explicitly loop over the physical domain + 2.
  do i = -1, nx+2
    do j = -1, ny+2
      do k = -1, nz+2
        ! Contribution of the charge density
        ! Calculate E = - \grad \phi
        N(i,j,k)%E(:) = 0.0_rk
#ifdef ELEC_NNONREST
        do veldir = 1, nnonrest
          nb_phi = get_nb_phi(N,veldir, i, j, k, ip, jp, kp)
          N(i,j,k)%E(1) = N(i,j,k)%E(1) + cx(veldir) * g(veldir) * nb_phi
          N(i,j,k)%E(2) = N(i,j,k)%E(2) + cy(veldir) * g(veldir) * nb_phi
          N(i,j,k)%E(3) = N(i,j,k)%E(3) + cz(veldir) * g(veldir) * nb_phi
        end do
        N(i,j,k)%E(:) = -N(i,j,k)%E(:) / 12.0_rk
#else
        do veldir = 1, nnn
          nb_phi = get_nb_phi(N,veldir, i, j, k, ip, jp, kp)
          N(i,j,k)%E(1) = N(i,j,k)%E(1) + cx(veldir) * nb_phi
          N(i,j,k)%E(2) = N(i,j,k)%E(2) + cy(veldir) * nb_phi
          N(i,j,k)%E(3) = N(i,j,k)%E(3) + cz(veldir) * nb_phi
        end do
        N(i,j,k)%E(:) = -0.5_rk * N(i,j,k)%E(:)
#endif
        ! Add external field
        N(i,j,k)%E(:) = N(i,j,k)%E(:) + E_ext(:)

        ! write(msgstr,"('E_tot = ',6(ES15.8,X))") N(i,j,k)%E(:), E_ext(:)
        ! call log_msg_elec(msgstr)

      end do
    end do
  end do

  call stop_timer(ti_elec_calc_E_fd)

end subroutine calc_E_fd
#endif

end module lbe_elec_poisson_solver_module

