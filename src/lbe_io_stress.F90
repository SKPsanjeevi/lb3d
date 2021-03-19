#include "lbe.h"

module lbe_io_stress_module

  use lbe_parms_module
  use lbe_globals_module
  use lbe_log_module
  use lbe_parallel_module
  use lbe_helper_module
  use lbe_io_helper_module
  use lbe_types_module, only: lbe_site

  implicit none
  include 'mpif.h'

  private
  public dump_stress, local_fluid_momentum_transfer, setup_stress 

contains

  subroutine dump_stress(N)
    implicit none
    type(lbe_site),intent(in) :: N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
    real(kind=rk) :: stress(4),stress_sum(4)
    character(len=1024) :: cfg_file_name
    integer,parameter :: cfg_file_unit=12
    logical :: fexist
    integer :: mpierror
    integer :: grange(2)

    stress(:) = 0.0d0
    
    call local_fluid_momentum_transfer(N,stress,grange)

    ! The momentum transfer through each of both planes is caused
    ! by forces from both sides. So the values for both planes are
    ! multiplied by one half here in order to let their sum be the
    ! force from both sides times one timestep.
    stress(:) = stress(:)*0.5_rk

    call MPI_Reduce(stress,stress_sum,4,LBE_REAL,MPI_SUM,0,Comm_Cart,mpierror)

    rank0: if (myrankc==0) then

      call lbe_make_filename_append(cfg_file_name,'stress','.asc')

      inquire(file=cfg_file_name,exist=fexist)
      if (fexist) then
        open(unit=cfg_file_unit,file=cfg_file_name,status='OLD',position='APPEND',recl=650)
      else
        open(unit=cfg_file_unit,file=cfg_file_name,status='NEW',position='APPEND',recl=650)
        ! If this is a new file, first write a header
        write (unit=cfg_file_unit, fmt='(A,X,A,X,F16.10)',advance='yes') "#","stress_dd_cutoff =", stress_dd_cutoff
        write (unit=cfg_file_unit, fmt='(A,X,A6,1X,4(A16,X))',advance='no') "#","t", "avg", "total", "stress1", "stress2"
        if (init_cond .eq. 14) then
          write (unit=cfg_file_unit, fmt='(4(A16,X),2(A6,X))',advance='no') "avg_drop", "total_drop", "stress1_drop", "stress2_drop", "minz", "maxz"
        endif
        write (unit=cfg_file_unit, fmt='()',advance='yes')
      endif

      write (unit=cfg_file_unit, fmt='(I8.8,1X)',advance='no') nt
      write (unit=cfg_file_unit, fmt='(4(E16.8,X))',advance='no') (-stress_sum(1)+stress_sum(2))/(real(tnz)*real(tny)),-stress_sum(1)+stress_sum(2),stress_sum(1), stress_sum(2)

      if (init_cond .eq. 14) then
        write (unit=cfg_file_unit, fmt='(4(E16.8,X))',advance='no') (-stress_sum(3)+stress_sum(4))/(real(grange(2)-grange(1)+1)*real(tny)),-stress_sum(3)+stress_sum(4),stress_sum(3), stress_sum(4)
        write (unit=cfg_file_unit, fmt='(2(I6,X))',advance='no') grange
      endif

      write (unit=cfg_file_unit, fmt='()',advance='yes')
      close (cfg_file_unit)
    end if rank0

  end subroutine dump_stress

  !> Adds to \c lfmt(1:2) the momentum transfer through the layers
  !> \c x==dfz_minx and \c x==dfz_maxx (as far as part of the local
  !> chunk) caused by the fluid itself.
  !>
  !> \todo This doesn't work for \c amass_(r|b|s)/=1.0
  subroutine local_fluid_momentum_transfer(N,lfmt,grange)
    implicit none
    type(lbe_site),intent(in) :: N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
    real(kind=rk),intent(out) :: lfmt(4)
    ! vectors that are advected in x direction AND carry momentum in z
    ! direction (p:positive x-direction, m:negative x-direction)
    integer,parameter,dimension(2) :: pvecs=(/9,10/),mvecs=(/13,14/)
    integer :: local_stress_minx,local_stress_maxx
    integer :: r,s,t,y,z
    integer, intent(out) :: grange(:)

    if (init_cond .eq. 14) then ! Only do this for single droplets
      call droplet_extent(N,stress_dd_cutoff,grange)
    else
      lfmt(3) = 0.0_rk
      lfmt(4) = 0.0_rk
      grange(:) = 0
    endif

    local_stress_minx = stress_minx+1-start(1)
    if (1<=local_stress_minx.and.local_stress_minx<=nx) then
      ! write(msgstr,"('prev stresses: ',2(F16.10,X))") lfmt
      ! call log_msg(msgstr,.true.)
      ! loop through lattice nodes on plane x=dfz_minx
      ! call log_msg("Looping through plane 1",.true.)
      do y=1,ny
        do z=1,nz
          do t=1,size(pvecs)
            s = pvecs(t)
            r = bounce(s)
            ! both source and target node have to be fluid nodes
            if (N(local_stress_minx-1,y,z-cz(s))%rock_state==0.0_rk&
                 &.and.N(local_stress_minx,y,z)%rock_state==0.0_rk) then
              ! sum incoming and outgoing momentum
              lfmt(1) = lfmt(1)&
                   &+N(local_stress_minx-1,y,z-cz(s))%n_r(s)*cz(s)*g(s)&
                   &-N(local_stress_minx,y,z)%n_r(r)*cz(r)*g(r)
#ifndef SINGLEFLUID
              lfmt(1) = lfmt(1)&
                   &+N(local_stress_minx-1,y,z-cz(s))%n_b(s)*cz(s)*g(s)&
                   &-N(local_stress_minx,y,z)%n_b(r)*cz(r)*g(r)
#ifndef NOSURFACTANT
              lfmt(1) = lfmt(1)&
                   &+N(local_stress_minx-1,y,z-cz(s))%n_s(s)*cz(s)*g(s)&
                   &-N(local_stress_minx,y,z)%n_s(r)*cz(r)*g(r)
#endif
#endif
              if (init_cond .eq. 14) then
                if ( (z + ccoords(3)*nz .ge. grange(1) ) .and. &
                     &(z + ccoords(3)*nz .le. grange(2) ) ) then
                  lfmt(3) = lfmt(3)&
                       &+N(local_stress_minx-1,y,z-cz(s))%n_r(s)*cz(s)*g(s)&
                       &-N(local_stress_minx,y,z)%n_r(r)*cz(r)*g(r)
#ifndef SINGLEFLUID
                  lfmt(3) = lfmt(3)&
                       &+N(local_stress_minx-1,y,z-cz(s))%n_b(s)*cz(s)*g(s)&
                       &-N(local_stress_minx,y,z)%n_b(r)*cz(r)*g(r)
#ifndef NOSURFACTANT
                  lfmt(3) = lfmt(3)&
                       &+N(local_stress_minx-1,y,z-cz(s))%n_s(s)*cz(s)*g(s)&
                       &-N(local_stress_minx,y,z)%n_s(r)*cz(r)*g(r)
#endif
#endif
                endif
              endif
            end if
          end do
        end do
      end do
      ! write(msgstr,"('summed stresses: ',4(F16.10,X))") lfmt
      ! call log_msg(msgstr,.true.)
    end if

    local_stress_maxx = stress_maxx+1-start(1)
    if (1<=local_stress_maxx.and.local_stress_maxx<=nx) then
      ! write(msgstr,"('prev stresses: ',2(F16.10,X))") lfmt
      ! call log_msg(msgstr,.true.)
      ! loop through lattice nodes on plane x=dfz_maxx
      ! call log_msg("Looping through plane 2",.true.)
      do y=1,ny
        do z=1,nz
          do t=1,size(pvecs)
            s = mvecs(t)
            r = bounce(s)
            ! both source and target node have to be fluid nodes
            if (N(local_stress_maxx+1,y,z-cz(s))%rock_state==0.0_rk&
                 &.and.N(local_stress_maxx,y,z)%rock_state==0.0_rk) then
              ! sum incoming and outgoing momentum
              lfmt(2) = lfmt(2)&
                   &+N(local_stress_maxx+1,y,z-cz(s))%n_r(s)*cz(s)*g(s)&
                   &-N(local_stress_maxx,y,z)%n_r(r)*cz(r)*g(r)
#ifndef SINGLEFLUID
              lfmt(2) = lfmt(2)&
                   &+N(local_stress_maxx+1,y,z-cz(s))%n_b(s)*cz(s)*g(s)&
                   &-N(local_stress_maxx,y,z)%n_b(r)*cz(r)*g(r)
#ifndef NOSURFACTANT
              lfmt(2) = lfmt(2)&
                   &+N(local_stress_maxx+1,y,z-cz(s))%n_s(s)*cz(s)*g(s)&
                   &-N(local_stress_maxx,y,z)%n_s(r)*cz(r)*g(r)
#endif
#endif
              if (init_cond .eq. 14) then
                if ( (z + ccoords(3)*nz .ge. grange(1) ) .and. &
                     &(z + ccoords(3)*nz .le. grange(2) ) ) then
                  lfmt(4) = lfmt(4)&
                       &+N(local_stress_maxx+1,y,z-cz(s))%n_r(s)*cz(s)*g(s)&
                       &-N(local_stress_maxx,y,z)%n_r(r)*cz(r)*g(r)
#ifndef SINGLEFLUID
                  lfmt(4) = lfmt(4)&
                       &+N(local_stress_maxx+1,y,z-cz(s))%n_b(s)*cz(s)*g(s)&
                       &-N(local_stress_maxx,y,z)%n_b(r)*cz(r)*g(r)
#ifndef NOSURFACTANT
                  lfmt(4) = lfmt(4)&
                       &+N(local_stress_maxx+1,y,z-cz(s))%n_s(s)*cz(s)*g(s)&
                       &-N(local_stress_maxx,y,z)%n_s(r)*cz(r)*g(r)
#endif
#endif
                endif
              endif
            end if
          end do
        end do
      end do
      ! write(msgstr,"('summed stresses: ',4(F16.10,X))") lfmt
      ! call log_msg(msgstr,.true.)
    end if

  end subroutine local_fluid_momentum_transfer

  !> setup stress boundary coordinates
  subroutine setup_stress()
    implicit none
    if (stress_minx<0) stress_minx = 1
    if (stress_maxx<0) stress_maxx = tnx
  end subroutine setup_stress

  subroutine droplet_extent(N,cutoff,grange)
    implicit none
    type(lbe_site),intent(in) :: N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
    real(kind=rk), intent(in) :: cutoff
    integer, intent(out) :: grange(2)
    integer :: lrange(2)
    integer :: mpierror

    lrange(1) = droplet_minz(N,cutoff)
    lrange(2) = droplet_maxz(N,cutoff)

    ! write(msgstr,"('minz = ',I0,' maxz = ',I0)") lrange(1), lrange(2)
    ! call log_msg(msgstr,.true.)

    call MPI_Allreduce(lrange(1),grange(1),1,MPI_INTEGER,MPI_MIN,Comm_cart,mpierror)
    call MPI_Allreduce(lrange(2),grange(2),1,MPI_INTEGER,MPI_MAX,Comm_cart,mpierror)

    ! write(msgstr,"('minz = ',I0,' maxz = ',I0)") grange(1), grange(2)
    ! call log_msg(msgstr)

  end subroutine droplet_extent

  pure integer function droplet_minz(N,cutoff)
    implicit none
    type(lbe_site),intent(in) :: N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
    real(kind=rk), intent(in) :: cutoff
    integer :: i,j,k

    droplet_minz = tnz

    do k=1,nz
      do j=1,ny
        do i=1,nx
          if (amass_r*sum(N(i,j,k)%n_r(:)*g) .gt. cutoff) then
            droplet_minz = k + ccoords(3)*nz
            return
          endif
        enddo
      enddo
    enddo
    return
  end function droplet_minz

  pure integer function droplet_maxz(N,cutoff)
    implicit none
    type(lbe_site),intent(in) :: N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
    real(kind=rk), intent(in) :: cutoff
    integer :: i,j,k

    droplet_maxz = 0

    do k=1,nz
      do j=1,ny
        do i=1,nx
          if (amass_r*sum(N(i,j,nz+1-k)%n_r(:)*g) .gt. cutoff) then
            droplet_maxz = (nz+1-k) + ccoords(3)*nz
            return
          endif
        enddo
      enddo
    enddo
    return
  end function droplet_maxz

end module lbe_io_stress_module
