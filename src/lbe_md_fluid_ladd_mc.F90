#include "lbe.h"

module lbe_md_fluid_ladd_mc_module
#ifdef MD

  use lbe_globals_module, only: g, halo_extent, n_spec, myrankc, tsize, nvecs
  use lbe_log_module
  use lbe_md_globals_module
  use lbe_md_helper_module, only: log_msg_md, md_make_filename_particles
  use lbe_helper_module, only: is_fluid, every_n_time_steps
  use lbe_parallel_module, only: comm_cart
  use lbe_parms_module, only: n_sanity_check, nt, nx, ny, nz, kbT, collisiontype_id, MRT_id,omegabulk_b,omegabulk_s,  omegabulk_r, s03_r, s05_r, omega_r, s11_r, s17_r, s03_b, s05_b, omega_b, s11_b, s17_b, s03_s, s05_s, omega_s, s11_s, s17_s
  use lbe_types_module, only: lbe_site
  use lbe_bdist_module


  implicit none
  include 'mpif.h'
  private

  public set_initial_average_density, calculate_total_mass, update_global_mass_change, reset_local_mass_change, create_mc_fluid_site, delete_mc_fluid_site, correct_mc_fluid_site
  public mass_correction, C0, n_recalculate_dmg, dump_masschange, dbg_report_mass_correction, max_mc
  public global_mass_target, global_mass_change, pfr, pfb, pfg, n_renew_rho

  !> \name Mass correction
  !> initial average densities and variables that hold the global
  !> mass error (global_mass_change) and the local mass error in the current step
  !> (local_mass_change)
  !> \{
  real(kind=rk),save :: pfr, pfb, pfg
  real(kind=rk),save :: local_mass_change(n_spec), global_mass_change(n_spec), global_mass_target(n_spec)
  !> \}

  !> \{
  !> \name mass correction
  !>
  !> enable mass correction for better mass conservation
  logical, save :: mass_correction = .false.
  real(kind=rk), save :: C0 = 100.0
  integer, save :: n_recalculate_dmg = 100
  integer, save :: n_renew_rho = 100
  logical, save :: dump_masschange = .false.
  logical, save :: dbg_report_mass_correction = .false.
  !> \}

  real(kind=rk),save :: max_mc = 0.5_rk

  integer :: n_ladd_removed = 0
  integer :: n_ladd_placed = 0
  integer :: n_ladd_lim(n_spec) = 0

  contains

subroutine reset_local_mass_change()
  implicit none

  if (dbg_report_mass_correction) call log_msg_md("Resetting local_mass_change to zero ...")
  local_mass_change(:) = 0.0_rk
end subroutine reset_local_mass_change

#ifdef SINGLEFLUID
subroutine create_mc_fluid_site(N, x, y, z, density, dist)
  implicit none
  type(lbe_site), intent(inout), target &
       &:: N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
  integer, intent(in) :: x, y, z
  real(kind=rk), intent(inout) :: density
  real(kind=rk), intent(in) :: dist(nvecs)

  real(kind=rk) :: corrected_density

  ! Calculate corrected density
  corrected_density = density * (1.0_rk - C0*global_mass_change(1)/product(tsize))
  ! Bookkeeping
  local_mass_change = local_mass_change + corrected_density
  n_ladd_placed = n_ladd_placed + 1
  ! Place fluid
  N(x,y,z)%n_r(:) = corrected_density * dist

! Thermalize newly created fluid 
if (collisiontype_id .eq. MRT_id) then
 if (kbT>0) then
  call mrt_dist( (/0.d0,0.d0,0.d0/) , (/0.d0,0.d0,0.d0/) , N(x, y, z)%n_r, (/ 1.d0, omegabulk_r, s03_r, 1.d0, s05_r, 1.d0, s05_r, 1.d0, s05_r, omega_r, s11_r, omega_r, s11_r, omega_r, omega_r, omega_r, s17_r, s17_r, s17_r /) )
 end if
end if

end subroutine create_mc_fluid_site

#else
! else BINARY

subroutine create_mc_fluid_site(N, x, y, z, density, dist)
  implicit none
  type(lbe_site), intent(inout), target &
       &:: N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
  integer, intent(in) :: x, y, z
  real(kind=rk), intent(inout) :: density(n_spec)
  real(kind=rk), intent(in) :: dist(nvecs)

  real(kind=rk) :: corrected_density(n_spec), C1(n_spec)

#ifdef NOSURFACTANT
  C1 = C0 * ( (pfr+pfb) / (/pfr,pfb/) ) * global_mass_change / product(tsize)
  C1 = max(min(C1,(/max_mc,max_mc/)),(/-max_mc,-max_mc/))
  if (dbg_report_mass_correction) then
    write(msgstr,"('density: ',2(E16.8,X))") density
    call log_msg_md(msgstr,.true.)
    write(msgstr,"('pfr / pfb:    ',2(E16.8,X), ' ; C1:     ',2(E16.8,X))") pfr, pfb, C1
    call log_msg_md(msgstr,.true.)
  end if
#else
  C1 = C0 * ( (pfr+pfb+pfg) / (/pfr,pfb,pfg/) ) * global_mass_change / product(tsize)
  C1 = max(min(C1,(/max_mc,max_mc,max_mc/)),(/-max_mc,-max_mc,-max_mc/))
  if (dbg_report_mass_correction) then
    write(msgstr,"('density: ',3(E16.8,X))") density
    call log_msg_md(msgstr,.true.)
    write(msgstr,"('pfr / pfb:    ',3(E16.8,X), ' ; C1:     ',3(E16.8,X))") pfr, pfb, pfg, C1
    call log_msg_md(msgstr,.true.)
  end if
#endif

  ! Calculate corrected density
  corrected_density = max(density * (1.0_rk - C1), 0.0_rk)

  ! Bookkeeping
  local_mass_change = local_mass_change + corrected_density
  n_ladd_placed = n_ladd_placed + 1

  ! Place fluid
  N(x,y,z)%n_r(:) = corrected_density(1) * dist
  N(x,y,z)%n_b(:) = corrected_density(2) * dist
#ifndef NOSURFACTANT
  N(x,y,z)%n_s(:) = corrected_density(3) * dist
#endif

  ! Set return value
  density = corrected_density

  if (dbg_report_mass_correction) then
#ifdef NOSURFACTANT
    write(msgstr,"('dist: ',19(E16.8,X), ' ; len: ',(E16.8,X))") dist, sum(dist*g)
    call log_msg_md(msgstr,.true.)
    write(msgstr,"('new fluid:    ',2(E16.8,X),' at ',3(I6,X))") &
         & sum(g*N(x,y,z)%n_r), sum(g*N(x,y,z)%n_b), x, y, z
    call log_msg_md(msgstr,.true.)
#else
    write(msgstr,"('Creating fluid ',3(E16.8,X),' at ',3(I6,X))") &
         & sum(g*N(x,y,z)%n_r),sum(g*N(x,y,z)%n_b),sum(g*N(x,y,z)%n_s), x, y, z
    call log_msg_md(msgstr,.true.)
#endif
! endif NOSURFACTANT
  end if

! Thermalize newly created fluid 
if (collisiontype_id .eq. MRT_id) then
 if (kbT>0) then
  call mrt_dist( (/0.d0,0.d0,0.d0/) , (/0.d0,0.d0,0.d0/) , N(x, y, z)%n_r, (/ 1.d0, omegabulk_r, s03_r, 1.d0, s05_r, 1.d0, s05_r, 1.d0, s05_r, omega_r, s11_r, omega_r, s11_r, omega_r, omega_r, omega_r, s17_r, s17_r, s17_r /) )
  call mrt_dist( (/0.d0,0.d0,0.d0/) , (/0.d0,0.d0,0.d0/) , N(x, y, z)%n_b, (/ 1.d0, omegabulk_b, s03_b, 1.d0, s05_b, 1.d0, s05_b, 1.d0, s05_b, omega_b, s11_b, omega_b, s11_b, omega_b, omega_b, omega_b, s17_b, s17_b, s17_b /) )
#ifndef NOSURFACTANT
  call mrt_dist( (/0.d0,0.d0,0.d0/) , (/0.d0,0.d0,0.d0/) , N(x, y, z)%n_s, (/ 1.d0, omegabulk_s, s03_s, 1.d0, s05_s, 1.d0, s05_s, 1.d0, s05_s, omega_s, s11_s, omega_s, s11_s, omega_s, omega_s, omega_s, s17_s, s17_s, s17_s /) )
#endif
 end if
end if

end subroutine create_mc_fluid_site

#endif
! endif SINGLEFLUID

subroutine delete_mc_fluid_site(N, x, y, z)
  implicit none
  type(lbe_site), intent(inout), target &
       &:: N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
  integer, intent(in) :: x, y, z

  local_mass_change(1) = local_mass_change(1) - sum(g*N(x,y,z)%n_r)
  N(x,y,z)%n_r(:) = 0.0_rk
#ifndef SINGLEFLUID
  local_mass_change(2) = local_mass_change(2) - sum(g*N(x,y,z)%n_b)
  N(x,y,z)%n_b(:) = 0.0_rk
#ifndef NOSURFACTANT
  local_mass_change(3) = local_mass_change(3) - sum(g*N(x,y,z)%n_s)
  N(x,y,z)%n_s(:) = 0.0_rk
#endif
#endif

  if (dbg_report_mass_correction) then
#ifdef SINGLEFLUID
    write(msgstr,"('Deleting fluid ',1(E16.8,X),' at ',3I6)") sum(g*N(x,y,z)%n_r),x,y,z
#else
#ifdef NOSURFACTANT
    write(msgstr,"('Deleting fluid ',2(E16.8,X),' at ',3I6)") sum(g*N(x,y,z)%n_r),sum(g*N(x,y,z)%n_b),x,y,z
#else
    write(msgstr,"('Deleting fluid ',3(E16.8,X),' at ',3I6)") sum(g*N(x,y,z)%n_r),sum(g*N(x,y,z)%n_b),sum(g*N(x,y,z)%n_s),x,y,z
#endif
#endif
    call log_msg_md(msgstr)
  end if

  ! keep track of number of removed fluid sites, but ignore timestep zero, where fluid will be removed during initial placement of particles
  if (nt > 1) n_ladd_removed = n_ladd_removed + 1

end subroutine delete_mc_fluid_site

subroutine correct_mc_fluid_site(N, x, y, z, s, correction)
  implicit none
  type(lbe_site), intent(inout), target &
       &:: N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
  integer, intent(in) :: x, y, z, s
  real(kind=rk), intent(in) :: correction(n_spec)

  local_mass_change(1) = local_mass_change(1) - correction(1)*g(s)
  N(x,y,z)%n_r(s) = N(x,y,z)%n_r(s) - correction(1)
#ifndef SINGLEFLUID
  local_mass_change(2) = local_mass_change(2) - correction(2)*g(s)
  N(x,y,z)%n_b(s) = N(x,y,z)%n_b(s) - correction(2)
#ifndef NOSURFACTANT
  local_mass_change(3) = local_mass_change(3) - correction(3)*g(s)
  N(x,y,z)%n_s(s) = N(x,y,z)%n_s(s) - correction(3)
#endif
#endif

  if (dbg_report_mass_correction) then
    if (any(.not.(correction*g(s) == 0.0_rk)) ) then
#ifdef SINGLEFLUID
      write(msgstr,"('Correcting fluid ',1(E16.8,X),'at ',3I6)") correction(:)*g(s),x,y,z
#else
#ifdef NOSURFACTANT
      write(msgstr,"('Correcting fluid ',2(E16.8,X),'at ',3I6)") correction(:)*g(s),x,y,z
#else
      write(msgstr,"('Correcting fluid ',3(E16.8,X),'at ',3I6)") correction(:)*g(s),x,y,z
#endif
#endif
      call log_msg_md(msgstr,.true.)
    end if
  end if

end subroutine correct_mc_fluid_site

!> This computes the initial average density of the different fluid species
!> and sets the global variables \c pfr[, \c pfb, \c pfg]
subroutine set_initial_average_density(N)
  implicit none
  type(lbe_site), intent(in), target &
       &:: N(1-halo_extent:,1-halo_extent:,1-halo_extent:)

  integer :: i,j,k
  real(kind=rk) :: lfr, lfb, lfg !< local fr, fb, fg
  integer :: mpierror
  integer :: sitecount
  real(kind=rk) :: ssitecount,totalsites

  call log_msg_md("  Setting initial average densities...")

  sitecount = 0
  pfr = 0.0_rk
  pfb = 0.0_rk
  pfg = 0.0_rk
  lfr = 0.0_rk
  lfb = 0.0_rk
  lfg = 0.0_rk

  do i=1,nx
    do j=1,ny
      do k=1,nz
        if ( is_fluid( N(i,j,k)%rock_state ) ) then
          sitecount = sitecount + 1
          lfr = lfr + sum(N(i,j,k)%n_r*g)
#ifndef SINGLEFLUID
          lfb = lfb + sum(N(i,j,k)%n_b*g)
#ifndef NOSURFACTANT
          lfg = lfg + sum(N(i,j,k)%n_s*g)
#endif
#endif
        end if
      end do
    end do
  end do

  ssitecount = real(sitecount,kind=rk) ! avoid overflow for huge systems
  call MPI_Allreduce(ssitecount, totalsites, 1,LBE_REAL, MPI_SUM, comm_cart, mpierror)
  call MPI_Allreduce(lfr, pfr, 1, LBE_REAL, MPI_SUM, Comm_Cart, mpierror)
  pfr = pfr / totalsites
  write(msgstr,"('    pfr = ',E16.8)") pfr
  call log_msg_md(msgstr)
#ifndef SINGLEFLUID
  call MPI_Allreduce(lfb, pfb, 1, LBE_REAL, MPI_SUM, Comm_Cart, mpierror)
  pfb = pfb / totalsites
  write(msgstr,"('    pfb = ',E16.8)") pfb
  call log_msg_md(msgstr)
#ifndef NOSURFACTANT
  call MPI_Allreduce(lfg, pfg, 1, LBE_REAL, MPI_SUM, Comm_Cart, mpierror)
  pfg = pfg / totalsites
  write(msgstr,"('    pfg = ',E16.8)") pfg
  call log_msg_md(msgstr)
#endif
#endif

end subroutine set_initial_average_density

!> This computes the initial average density of the different fluid species
!> and stores it in \c tm
subroutine calculate_total_mass(N, total_mass)
  implicit none
  type(lbe_site),intent(in),target &
       &:: N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
  real(kind=rk),intent(out) :: total_mass(n_spec)

  integer :: i, j, k !< Loop dummy variables
  real(kind=rk) :: lm(n_spec) !< Local masses
  integer :: mpierror

  if (dbg_report_mass_correction) call log_msg_md("Calculating total masses through calculate_total_mass() ...")

  lm(:) = 0.0_rk

  ! Sum local masses
  do i=1,nx
    do j=1,ny
      do k=1,nz
        if ( is_fluid( N(i,j,k)%rock_state ) ) then
          lm(1) = lm(1) + sum(N(i,j,k)%n_r*g)
#ifndef SINGLEFLUID
          lm(2) = lm(2) + sum(N(i,j,k)%n_b*g)
#ifndef NOSURFACTANT
          lm(3) = lm(3) + sum(N(i,j,k)%n_s*g)
#endif
#endif
        end if
      end do
    end do
  end do

  call MPI_Allreduce(lm,total_mass,n_spec,LBE_REAL,MPI_SUM,Comm_Cart,mpierror)

  if (dbg_report_mass_correction) then
    do i=1,n_spec
      write(msgstr,"('  total_mass(',I0,') = ',E16.8)") i, total_mass(i)
      call log_msg_md(msgstr)
    end do
  end if

end subroutine calculate_total_mass

!> This routine collects the local mass errors calculated during the current
!> timestep and broadcasts the total to all CPUs.
subroutine update_global_mass_change(N)
  implicit none
  type(lbe_site),intent(in) :: N(1-halo_extent:,1-halo_extent:,1-halo_extent:)

  real(kind=rk) :: gmctmp(n_spec) !< Global Mass Change temporary variable
  integer :: mpierror

  integer :: i

  ! Sum all local mass changes into gmctmp -- these are mass changes for this timestep
  call MPI_Allreduce(local_mass_change, gmctmp, n_spec, LBE_REAL, MPI_SUM, Comm_Cart, mpierror)

  if (dbg_report_mass_correction) then
    call log_msg_md("Updating global mass change through update_global_mass_change() ...")
    do i=1,n_spec
      write(msgstr,"('  local_mass_change(',I0,') = ',E16.8)") i, local_mass_change(i)
      call log_msg_md(msgstr)
      write(msgstr,"('  gmctmp(',I0,')            = ',E16.8)") i, gmctmp(i)
      call log_msg_md(msgstr)
    end do
  end if

  ! Add the mass change of this timestep to the global mass change
  global_mass_change = global_mass_change + gmctmp

  ! Report values
  if ( every_n_time_steps(n_sanity_check) ) then
    call log_msg_md("Updating global mass change through update_global_mass_change() ...")
    do i=1,n_spec
      write(msgstr,"('  global_mass_target(',I0,') = ',E16.8)") i, global_mass_target(i)
      call log_msg_md(msgstr)
      write(msgstr,"('  global_mass_change(',I0,') = ',E16.8, F16.10)") i, global_mass_change(i), 1.0 + global_mass_change(i)/global_mass_target(i)
      call log_msg_md(msgstr)
    end do
  end if

  if ( every_n_time_steps(n_recalculate_dmg) ) then
    call calculate_total_mass(N, gmctmp)
    global_mass_change = gmctmp - global_mass_target
    call log_msg_md("Updating global mass change through synchronization ...")
    do i=1,n_spec
      write(msgstr,"('  global_mass_target(',I0,') = ',E16.8)") i, global_mass_target(i)
      call log_msg_md(msgstr)
      write(msgstr,"('  global_mass_change(',I0,') = ',E16.8, F16.10)") i, global_mass_change(i), 1.0 + global_mass_change(i)/global_mass_target(i)
      call log_msg_md(msgstr)
    end do
    if (dump_masschange) call write_masschange()
  end if

end subroutine update_global_mass_change

!> writes prepared data into asc file on rank 0
subroutine write_masschange()
  implicit none

  character(len=1024) cfg_file_name
  integer,parameter :: cfg_file_unit=12
  logical :: fexist
  integer :: n_ladd_removed_g, n_ladd_placed_g
  integer :: n_ladd_lim_g(n_spec)
  integer :: mpierror
  integer :: i

  call MPI_Reduce(n_ladd_removed,n_ladd_removed_g,1,MPI_INTEGER,MPI_SUM,0,Comm_Cart,mpierror)
  call MPI_Reduce(n_ladd_placed,n_ladd_placed_g,1,MPI_INTEGER,MPI_SUM,0,Comm_Cart,mpierror)
  call MPI_Reduce(n_ladd_lim, n_ladd_lim_g, n_spec ,MPI_INTEGER,MPI_SUM,0,Comm_Cart,mpierror)

  rank0: if (myrankc==0) then

    call log_msg_md("Dumping mass data ...")

    call md_make_filename_particles(cfg_file_name,'massgain','.asc',0)

    inquire(file=cfg_file_name,exist=fexist)
    file_exists: if (fexist) then
      open(unit=cfg_file_unit,file=cfg_file_name,status='OLD'&
           &,position='APPEND',recl=650)
    else
      open(unit=cfg_file_unit,file=cfg_file_name,status='NEW'&
           &,position='APPEND',recl=650)
      ! Write header data (only do this for a new file)
      write(unit=cfg_file_unit, fmt="('# ',A8,X)",advance='no') "nt"

      write(unit=cfg_file_unit, fmt="(A16,X,A16,X)",advance='no') "gmc(1)", "gmcrel(1)"
#ifndef SINGLEFLUID
      write(unit=cfg_file_unit, fmt="(A16,X,A16,X)",advance='no') "gmc(2)", "gmcrel(2)"
#ifndef NOSURFACTANT
      write(unit=cfg_file_unit, fmt="(A16,X,A16,X)",advance='no') "gmc(3)", "gmcrel(3)"
#endif
#endif

      write(unit=cfg_file_unit, fmt="(2(A10,X))",advance='no') "removed", "placed"
      write(unit=cfg_file_unit, fmt="(1(A10,X))",advance='no') "lim_r"
#ifndef SINGLEFLUID
      write(unit=cfg_file_unit, fmt="(1(A10,X))",advance='no') "lim_b"
#ifndef NOSURFACTANT
      write(unit=cfg_file_unit, fmt="(1(A10,X))",advance='no') "lim_s"
#endif
#endif
      write(unit=cfg_file_unit, fmt='()',advance='yes')
    end if file_exists ! End new file / header writing

    write (unit=cfg_file_unit, fmt='(I10.10,1X)',advance='no') nt

    do i=1,n_spec
      write (unit=cfg_file_unit, fmt='(E16.8,X,F16.10,X)',advance='no') global_mass_change(i), 1.0+global_mass_change(i)/global_mass_target(i)
    end do

    write (unit=cfg_file_unit, fmt='(2(I10,1X))',advance='no') n_ladd_removed_g, n_ladd_placed_g
    do i=1,n_spec
      write (unit=cfg_file_unit, fmt='(1(I10,1X))',advance='no') n_ladd_lim_g(i)
    end do

    write (unit=cfg_file_unit, fmt='()',advance='yes')
    close (cfg_file_unit)
  end if rank0
end subroutine write_masschange

#endif
end module lbe_md_fluid_ladd_mc_module
