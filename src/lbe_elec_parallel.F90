#include "lbe.h"

module lbe_elec_parallel_module
#ifdef ELEC

  use lbe_elec_globals_module
  use lbe_elec_helper_module
  use lbe_parallel_module, only: build_all_chunk_mpitypes, halo, checkmpi
  use mpi

  implicit none

  private

  public elec_setup_parallel, elec_halo, phi_halo

  type(halo) :: elec_halo
  type(halo) :: phi_halo

  ! Worst-case estimate: charge advection requires forces in the first halo layer,
  ! which requires the electric field in the second halo layer (dielectrophoretic term),
  ! which requires the electric potential in the third halo layer.
  ! A first attempt to change this parameter to be more flexible resulted in code
  ! too complicated for the associated performance gain.
  integer, parameter :: elec_halo_extent = 3
  integer, parameter :: phi_halo_extent = 3

contains

!> Create MPI types: basic types for lattice-site elec-only halo exchange and
!> phi-only halo exchange as well as types for the halo chunks.
subroutine elec_setup_parallel(N)
  implicit none

  type(lbe_site), intent(inout) :: N(1-halo_extent:, 1-halo_extent:, 1-halo_extent:)

  integer :: elec_mpitype, phi_mpitype

  call log_msg_elec_hdr("Setting up ELEC parallelization")

  call log_msg_elec("Creating MPI datatype for ELEC fields ...")
  call build_elec_site_mpitype(elec_mpitype)
  call log_msg_elec("Creating MPI datatype for ELEC phi ...")
  call build_phi_site_mpitype(phi_mpitype)
  call log_msg_elec("Creating MPI datatypes for ELEC field chunks ...")
  call build_all_chunk_mpitypes(N, elec_halo, elec_mpitype, elec_halo_extent)
  call log_msg_elec("Creating MPI datatypes for ELEC phi chunks ...")
  call build_all_chunk_mpitypes(N, phi_halo, phi_mpitype, phi_halo_extent)

  call log_msg_elec_ws("Finished setting up ELEC parallelization.")

end subroutine elec_setup_parallel

!> Set up the MPI type for the rho_p, rho_m, phi and eps components of the lattice.
subroutine build_elec_site_mpitype(mpitype)
  implicit none

  integer, intent(out) :: mpitype ! type to be built

  integer, parameter   :: cnt = 4 ! rho_p, rho_m, phi, eps

  type(lbe_site) :: sample
  integer(kind = MPI_ADDRESS_KIND) :: addrs(cnt), base, displs(cnt)
  integer :: blocklengths(cnt), types(cnt)
  integer :: mpierror

  ! Get base address of an lbe_site
  call MPI_Get_address(sample, base, mpierror)
  DEBUG_CHECKMPI(mpierror, "build_elec_site_mpitype, MPI_Get_address : base")

  ! rho_p
  blocklengths(1) = 1
  types(1) = LBE_REAL
  call MPI_Get_address(sample%rho_p,addrs(1), mpierror)
  DEBUG_CHECKMPI(mpierror, "build_elec_site_mpitype, MPI_Get_address : 1")

  ! rho_m
  blocklengths(2) = 1
  types(2) = LBE_REAL
  call MPI_Get_address(sample%rho_m,addrs(2), mpierror)
  DEBUG_CHECKMPI(mpierror, "build_elec_site_mpitype, MPI_Get_address : 2")

  ! phi
  blocklengths(3) = 1
  types(3) = LBE_REAL
  call MPI_Get_address(sample%phi  ,addrs(3), mpierror)
  DEBUG_CHECKMPI(mpierror, "build_elec_site_mpitype, MPI_Get_address : 3")

  ! eps
  blocklengths(4) = 1
  types(4) = LBE_REAL
  call MPI_Get_address(sample%eps  ,addrs(4), mpierror)
  DEBUG_CHECKMPI(mpierror, "build_elec_site_mpitype, MPI_Get_address : 4")

  displs(:) = addrs(:) - base

  call MPI_Type_create_struct(cnt, blocklengths, displs, types, mpitype, mpierror)
  DEBUG_CHECKMPI(mpierror, "build_elec_site_mpitype, MPI_Type_create_struct")
  call MPI_Type_commit(mpitype, mpierror)
  DEBUG_CHECKMPI(mpierror, "build_elec_site_mpitype, MPI_Type_commit")

end subroutine build_elec_site_mpitype

!> Set up the MPI type for the phi component of the lattice only.
subroutine build_phi_site_mpitype(mpitype)
  implicit none

  integer, intent(out) :: mpitype ! type to be built

  integer, parameter :: cnt = 1 ! only phi

  type(lbe_site) :: sample
  integer(kind = MPI_ADDRESS_KIND) :: addrs(cnt), base, displs(cnt)
  integer :: blocklengths(cnt), types(cnt)
  integer :: mpierror

  ! Get base address of an lbe_site
  call MPI_Get_address(sample, base, mpierror)
  DEBUG_CHECKMPI(mpierror, "build_phi_site_mpitype, MPI_Get_address : base")

  ! phi
  blocklengths(1) = 1
  types(1) = LBE_REAL
  call MPI_Get_address(sample%phi,addrs(1), mpierror)
  DEBUG_CHECKMPI(mpierror, "build_phi_site_mpitype, MPI_Get_address : 1")

  ! Calculate displacements
  displs(:) = addrs(:) - base

  call MPI_Type_create_struct(cnt, blocklengths, displs, types, mpitype, mpierror)
  DEBUG_CHECKMPI(mpierror, "build_phi_site_mpitype, MPI_Type_create_struct")
  call MPI_Type_commit(mpitype, mpierror)
  DEBUG_CHECKMPI(mpierror, "build_phi_site_mpitype, MPI_Type_commit")

end subroutine build_phi_site_mpitype

#endif
end module lbe_elec_parallel_module

