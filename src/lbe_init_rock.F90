#include "lbe.h"

module lbe_init_rock_module

  use lbe_globals_module
  use lbe_io_module, only: read_rock_all_par
#ifdef USEHDF
  use lbe_io_hdf5_module, only: read_scalar_phdf5
#endif
#ifdef USEXDRF
  use lbe_io_xdrf_module, only: read_bit_rock_xdrf_par, read_rock_xdrf_par

! Intermediate Knudsen initialisation requires additional information from files
#ifdef LOCALBC
  use lbe_io_xdrf_module, only: read_bc_xdrf_par
#endif
#ifdef VARTAU 
  use lbe_io_xdrf_module, only: read_rel_xdrf_par
#endif
#endif
  use lbe_log_module
  use lbe_parallel_module, only: checkmpi,recv_rock_par, tnx, tny, tnz, ccoords
  use lbe_parallel_module, only: build_all_chunk_mpitypes, halo_exchange, halo, rock_halo
  use lbe_parms_module, only: bcsel, boundary_cond, boundary_width, obs_file,&
       obs_folder, nx, ny, nz, rock_colour_init, rock_colour,& 
       gw_rock_colour, rock_colour_double, g_value_double, gw_double_wet,& 
       obs_rotation, obs_file_r, obs_file_b
  use lbe_types_module
  use lbe_helper_module

#ifdef DEBUG_MPI
    use lbe_parallel_module, only: checkmpi
#endif

  implicit none
  include 'mpif.h'
  private

  public inside_rock_global,lbe_init_rock,lbe_init_rock_parallel,rock_is_present

contains

    !> Same as \c inside_rock() but this one works on global \c
    !> rock_state(:,:,:) instead of the local chunk \c N(:,:,:) , so it
    !> can examine every possible position.
    logical function inside_rock_global(rock_state,pos)
        real(kind=rk),intent(in) :: rock_state&
             &(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        real(kind=rk),intent(in) :: pos(3)
        integer fp(3)

        fp = floor(pos)
        inside_rock_global = all(rock_state&
             &(fp(1)+lp_sur(1,:),fp(2)+lp_sur(2,:),fp(3)+lp_sur(3,:))&
             &/=0.0_rk)
    end function inside_rock_global

!> This routine will handle initialization of the rock geometry.
!> It is the only routine that is exposed to other modules.
subroutine lbe_init_rock(whole_N)
  implicit none
  type(lbe_site),dimension(1-halo_extent:,1-halo_extent:,1-halo_extent:),intent(inout),target :: whole_N
  type(lbe_site),pointer :: N(:,:,:)

  N => whole_N(0:nx+1,0:ny+1,0:nz+1)

  call log_msg("Setting up rock geometry.")

  ! First, read rock from file, if any...
  if ( trim(obs_file) .ne. 'empty.dat') then
    call lbe_read_rock(N)
  else
    call log_msg("Special value obs_file = <empty.dat> : no rock file will be read.")
  end if
  if (rock_colour_double .or. g_value_double .or. gw_double_wet) then
  if ( (trim(obs_file_r) .ne. 'empty.dat') .AND. trim(obs_file_b) .ne. 'empty.dat') then
    call lbe_read_rock_rb(N)
  else
    call log_msg("Special value obs_file_rb = <empty.dat> : no rock file will be read.")
   call error("ERROR:,double rock are used, but no red or blue rock file supplied - can'tinitialize rock. Aborting...") 
  end if
end if 
  ! Then, modify, if boundary_cond > 0.
  if ( boundary_cond .gt. 0 ) then
    call log_msg("  boundary_cond > 0: Adding additional rock geometry...")
    call lbe_modify_rock(N)
  else if ( boundary_cond .eq. 0 ) then
    call log_msg("  boundary_cond == 0: No additional rock geometry will be added (and no special BCs will be applied).")
  else
    call log_msg("  boundary_cond < 0: No additional rock geometry will be added (but special BCs might be applied).")
  end if

  ! We only need this if we do not get wettability information from rock
  ! files and do not use read_rock_xdrf_par.
  call lbe_init_rock_colour(N)

  ! Exchange rock information in halo (needed to mark the surface nodes)
  call lbe_init_rock_parallel(whole_N)
  call halo_exchange(whole_N, rock_halo)

  ! Mark surface nodes, i.e., rock nodes with at least one fluid neighbour
  call lbe_init_rock_surface(whole_N)

end subroutine lbe_init_rock

!> Wraps calls to various rock reading functions.
subroutine lbe_read_rock(N)
  implicit none

  type(lbe_site),dimension(0:,0:,0:),intent(inout) :: N
  character(len=1024) :: full_obs_fname
  
  full_obs_fname = trim(obs_folder)//'/'//trim(obs_file)
  write(msgstr,"('Reading rock data file <',A,'>')") trim(full_obs_fname)
  call log_msg(msgstr)

  if (index(full_obs_fname,'.xdr') .gt. 0) then
#ifdef USEXDRF
    ! Call XDR-read code
    call read_rock_xdrf_par(full_obs_fname,N,0)

#ifdef DIST
    call read_dist_xdrf_par(full_obs_fname,N)
#endif

#ifdef VARTAU
! Timm: This subroutine call has been deactivated for now as the feature does not seem to be required right now.
! If necessary, it should be checked later if this subroutine still does what it was supposed to do a few years ago.
     call read_rel_xdrf_par(full_obs_fname,N)
#endif

#ifdef LOCALBC
    call read_bc_xdrf_par(full_obs_fname,N)
#endif

#else
    ! If USEXDRF is not set:
    call error("FATAL ERROR: XDRF disabled, but xdr rock file supplied - can't initialize rock. Aborting...")
#endif

  else if (index(full_obs_fname,'.h5') .gt. 0) then
#ifdef USEHDF
    call lbe_read_rock_hdf5(full_obs_fname,N)
#else
    call error("FATAL ERROR: HDF5 disabled, but hdf rock file supplied - can't initialize rock. Aborting...")
#endif

  else if (index(full_obs_fname,'.bit') .gt. 0) then
#ifdef USEXDRF
    ! Call XDR-read code
    call read_bit_rock_xdrf_par(full_obs_fname,N)
#else
    ! If USEXDRF is not set:
    call error("FATAL ERROR: XDRF disabled, but xdr rock file supplied - can't initialize rock. Aborting...")
#endif

  else
    if (myrankc == 0) then
      ! Call parallel ASCII read code
      ! This will not give correct rock_state value to rock lattice if compiled
      ! with -DMD.
      call read_rock_all_par(full_obs_fname,N)
    else
      ! Receive rock from rank zero
      call recv_rock_par(N)
    endif
  endif
end subroutine lbe_read_rock

!> red rock reading functions.
subroutine lbe_read_rock_rb(N)
  implicit none

  type(lbe_site),dimension(0:,0:,0:),intent(inout) :: N
  character(len=1024) :: full_obs_fname_r
  character(len=1024) :: full_obs_fname_b
	full_obs_fname_r = trim(obs_folder)//'/'//trim(obs_file_r)	
   full_obs_fname_b = trim(obs_folder)//'/'//trim(obs_file_b)

  if ((index(full_obs_fname_r,'.xdr') .gt. 0) .AND.(index(full_obs_fname_b,'.xdr') .gt. 0) ) then
#ifdef USEXDRF
    ! Call XDR-read code
     write(msgstr,"('Reading red rock data file <',A,'>')")trim(full_obs_fname_r)
     call log_msg(msgstr)
     call read_rock_xdrf_par(full_obs_fname_r,N,1)

     write(msgstr,"('Reading blue rock data file <',A,'>')") trim(full_obs_fname_b)
     call log_msg(msgstr)
     call read_rock_xdrf_par(full_obs_fname_b,N,2)
#else
    ! If USEXDRF is not set:
    call error("FATAL ERROR: XDRF disabled, but xdr rock file supplied - can't initialize rock. Aborting...")
#endif
else if ((index(full_obs_fname_r,'.h5') .gt. 0) .AND.(index(full_obs_fname_b,'.h5') .gt. 0) ) then
#ifdef USEHDF
write(msgstr,"('Reading red rock data file <',A,'>')")trim(full_obs_fname_r)
       call log_msg(msgstr)    
call lbe_read_rock_hdf5_rb(full_obs_fname_r,N,1)
    
write(msgstr,"('Reading blue rock data file <',A,'>')")trim(full_obs_fname_b)
       call log_msg(msgstr)
call lbe_read_rock_hdf5_rb(full_obs_fname_b,N,2)
#else
    call error("FATAL ERROR: HDF5 disabled, but hdf rock file supplied - can't initialize rock. Aborting...")
#endif

 else 
    call error("FATAL ERROR: no xdr red or blue rock file supplied - can't initialize rock. Aborting...")
  endif
end subroutine lbe_read_rock_rb

#ifdef USEHDF
subroutine lbe_read_rock_hdf5(filename, N)
  implicit none
  character(len=*), intent(in) :: filename
  type(lbe_site),dimension(0:,0:,0:), intent(inout) :: N

  real(kind=rk), dimension(:,:,:), allocatable :: rock

  real(kind=rk) :: rockfloat  ! Temporary variable to store a value of scalar in
  integer :: rnx, rny, rnz, roti, rotj, rotk
  integer :: i,j,k            ! Dummy loop variables

  integer :: ierror ! For allocation errors

  call rotate_coordinates(obs_rotation, nx, ny, nz, rnx, rny, rnz)
  allocate(rock(1:rnx,1:rny,1:rnz),stat=ierror)

  if (ierror .ne. 0) then
    call error("WARNING: unable to allocate scalar buffer for rock reading")
  end if

  call read_scalar_phdf5(rock, filename, obs_rotation)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        call rotate_coordinates(obs_rotation,i,j,k,roti,rotj,rotk)
        rockfloat = rock(roti,rotj,rotk)
        ! TODO: Fix this to match new rock writing routine?
        ! Code below copied from read_rock_xdrf_par function

        if ( is_wall(rockfloat) ) then
          N(i,j,k)%rock_state = rock_value
          N(i,j,k)%rock_colour = real(rockfloat, kind=rk) - 5.0_rk
       else
          N(i,j,k)%rock_state = no_rock_value
          N(i,j,k)%rock_colour = 0.d0
        endif

      enddo
    enddo
  enddo
  deallocate(rock)

end subroutine lbe_read_rock_hdf5

subroutine lbe_read_rock_hdf5_rb(filename, N,rb)
  implicit none
  character(len=*), intent(in) :: filename
  type(lbe_site),dimension(0:,0:,0:), intent(inout) :: N

  real(kind=rk), dimension(:,:,:), allocatable :: rock

  real(kind=rk) :: rockfloat  ! Temporary variable to store a value of scalar in
  integer :: rnx, rny, rnz, roti, rotj, rotk
  integer :: i,j,k            ! Dummy loop variables
  integer :: rb
  integer :: ierror ! For allocation errors

  call rotate_coordinates(obs_rotation, nx, ny, nz, rnx, rny, rnz)
  allocate(rock(1:rnx,1:rny,1:rnz),stat=ierror)

  if (ierror .ne. 0) then
    call error("WARNING: unable to allocate scalar buffer for rock reading")
  end if

  call read_scalar_phdf5(rock, filename, obs_rotation)
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        call rotate_coordinates(obs_rotation,i,j,k,roti,rotj,rotk)
        rockfloat = rock(roti,rotj,rotk)
        ! TODO: Fix this to match new rock writing routine?
        ! Code below copied from read_rock_xdrf_par function

       ! if ( is_wall(rockfloat) ) then
      !    N(i,j,k)%rock_state = rock_value
     !     N(i,j,k)%rock_colour = real(rockfloat, kind=rk) - 5.0_rk
       !else
      !    N(i,j,k)%rock_state = no_rock_value
       !   N(i,j,k)%rock_colour = 0.d0
     !   endif

if (rb==0) then   ! read rock_colour for all
    if (rockfloat /= 0.0) then
       N(i,j,k)%rock_state = rock_value
             N(i,j,k)%rock_colour = real(rockfloat, kind=rk) - 5.0_rk
 				 !rcount = rcount+1     
   else 
       N(i,j,k)%rock_state = 0.d0
       N(i,j,k)%rock_colour = 0.d0
      end if
   end if
   
   if (rb == 1) then ! read rock_colour_r for red rock
      
      if (rockfloat /= 0.0) then
         N(i,j,k)%rock_state = rock_value
         if (g_value_double) then
            N(i,j,k)%rock_colour = 1.d0 !dble(rockfloat) - 5.0d0
            N(i,j,k)%rock_colour_r = real(rockfloat, kind=rk) ! - 5.0d0
         else
            N(i,j,k)%rock_colour = real(rockfloat, kind=rk) - 5.0_rk
            N(i,j,k)%rock_colour_r = real(rockfloat, kind=rk) - 5.0_rk
         end if
         !rcount = rcount+1
         else
            N(i,j,k)%rock_state = 0.d0
            N(i,j,k)%rock_colour = 0.d0
            N(i,j,k)%rock_colour_r = 0.d0
         end if
      end if

      if (rb == 2) then ! read rock_colour_b for blue rock, do not rewrite rock_state of red part
         if (rockfloat /= 0.0) then
            N(i,j,k)%rock_state = rock_value
			 if (g_value_double) then
            N(i,j,k)%rock_colour = 1.d0 ! N(x,y,z)%rock_colour- (dble(rockfloat) - 5.0d0) ! rock_colour_r - rock_colour_b
            N(i,j,k)%rock_colour_b = real(rockfloat, kind=rk) ! - 5.0d0
			 else
				  N(i,j,k)%rock_colour = N(i,j,k)%rock_colour- (real(rockfloat, kind=rk) - 5.0_rk) !rock_colour_r - rock_colour_b
              N(i,j,k)%rock_colour_b = real(rockfloat, kind=rk) - 5.0_rk
			 end if
          ! rcount = rcount+1
        else
            N(i,j,k)%rock_colour_b = 0.d0
         end if
      end if


      enddo
    enddo
  enddo
  deallocate(rock)

end subroutine lbe_read_rock_hdf5_rb


#endif

subroutine lbe_init_rock_colour(N)
  implicit none
  type(lbe_site),dimension(0:,0:,0:),intent(inout) :: N
  integer :: x,y,z

  ! Clear the nonrest vectors of each rock site.
  ! FIXME - this could be cleaned up, although it only
  ! gets called once, so it's not that important.

  do z = 1, nz
    do y = 1, ny
      do x = 1, nx
       ! if ( .not. is_fluid(N(x,y,z)%rock_state) ) then
        if ( is_wall(N(x,y,z)%rock_state) ) then
         
          ! N(x,y,z)%n_r(:nnonrest)=g_inv*fr/nvecs

 if (rock_colour_double .or. g_value_double .or. gw_rock_colour) then

!!Code before restoring single rock_colour function
            N(x,y,z)%n_r(restvec)=0.d0
#ifndef SINGLEFLUID
            N(x,y,z)%n_b(restvec)=0.d0
#ifndef NOSURFACTANT
            N(x,y,z)%n_s(restvec)=0.d0
#endif
#endif

         
else 

 if (N(x,y,z)%rock_colour > 0.d0) then
            N(x,y,z)%n_r(restvec)= N(x,y,z)%rock_colour
          else 
            N(x,y,z)%n_r(restvec)= 0.d0
          end if
#ifndef SINGLEFLUID
          if (N(x,y,z)%rock_colour < 0.d0) then
            N(x,y,z)%n_b(restvec)= -N(x,y,z)%rock_colour
          else
            N(x,y,z)%n_b(restvec)= 0.d0
          end if
#ifndef NOSURFACTANT
          N(x,y,z)%n_s(restvec)=0.d0
#endif
#endif
end if

 if (bcsel.ne.2) then
            N(x,y,z)%n_r(:nnonrest)=0.d0
#ifndef SINGLEFLUID
            N(x,y,z)%n_b(:nnonrest)=0.d0
#ifndef NOSURFACTANT
            N(x,y,z)%n_s(:nnonrest)=0.d0
#endif
#endif
          end if

          if (rock_colour_init) then
            ! setting the rock colour
            N(x,y,z)%rock_colour = rock_colour
          ! this "else" code is a bit funny, since it resets the rock_colour to zero which is
          ! read from rock file previously, thus making dump rock_colour not
          ! consistent with input. So I comment it. (Qingguang)
          !else
          !  N(x,y,z)%rock_colour = 0.d0
          endif

        endif
      end do
    end do
  end do
end subroutine lbe_init_rock_colour

!> Mark surface nodes, i.e., rock nodes with at least one fluid neighbour
subroutine lbe_init_rock_surface(N)
  implicit none

  type(lbe_site),intent(inout) :: N(1-halo_extent:,1-halo_extent:,1-halo_extent:)

  integer :: x, y, z

  ! loop over all nodes (halo excluded)
  do x = 1, nx
     do y = 1, ny
        do z = 1, nz
           ! check if rock node is surface site by looking at all neighbours
           if ( is_rock( N(x,y,z)%rock_state ) &
                .and. &
                ! check if any neighbour is a fluid node
                any( is_fluid( N(x+cx(1:nnonrest),y+cy(1:nnonrest),z+cz(1:nnonrest))%rock_state ) ) ) &
                ! mark rock node as a surface site
                N(x,y,z)%rock_state = rock_value_surface
        end do
     end do
  end do

end subroutine lbe_init_rock_surface

!> Modify geometries to add simple walls and tunnels.
subroutine lbe_modify_rock(N)
  implicit none
  type(lbe_site),dimension(0:,0:,0:),intent(inout) :: N

  ! Values to set the 3 sides to
  integer :: xval, yval, zval, mval

  ! Dummy indices
  integer :: i,j,k
  integer :: ti, tj, tk

  ! Values to set on the edges:
  !  1: rock
  !  0: no rock
  ! -1: do not modify
  ! The mapping to boundary_cond values are chosen to keep backwards compatibility
  select case(boundary_cond)
  case(1)
    call log_msg("  Adding rock cage.")
    xval = 1
    yval = 1
    zval = 1
  case(2)
    call log_msg("  Adding square duct in z-direction with cleared sites.")
    xval = 1
    yval = 1
    zval = 0
  case(3)
    call log_msg("  Clearing x and y sides.")
    xval = 0
    yval = 0
    zval = -1
  case(4)
    call log_msg("  Adding square duct in z-direction.")
    xval = 1
    yval = 1
    zval = -1
  case(5)
    call log_msg("  Clearing z sides.")
    xval = -1
    yval = -1
    zval = 0
  case(6)
    call log_msg("  Adding walls on x sides.")
    xval = 1
    yval = -1
    zval = -1
  case default
    write(msgstr,"('Unknown value <',I0,'> for additional rock geometry.')") boundary_cond
    call error(msgstr)
  end select

  ! Loop over the local domain, and set relevant rock states
  ! Rock wins over no rock wins over NOOP
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        ! Positions on the global lattice
        ti = i + ccoords(1) * nx
        tj = j + ccoords(2) * ny
        tk = k + ccoords(3) * nz
        
        if ( ( ti .le. boundary_width ) .or. ( ti .gt. tnx - boundary_width ) ) then

          if ( ( tj .le. boundary_width ) .or. ( tj .gt. tny - boundary_width ) ) then
            
            if ( ( tk .le. boundary_width ) .or. ( tk .gt. tnz - boundary_width ) ) then
              mval = max(xval, max(yval, zval) )
              if (mval .ge. 0) N(i,j,k)%rock_state = rock_value
            else
              mval = max(xval, yval)
              if (mval .ge. 0) N(i,j,k)%rock_state = rock_value
            endif

          else
            if ( ( tk .le. boundary_width ) .or. ( tk .gt. tnz - boundary_width ) ) then
              mval = max(xval, zval)
              if (mval .ge. 0) N(i,j,k)%rock_state = rock_value
            else
              mval = xval
              if (mval .ge. 0) N(i,j,k)%rock_state = rock_value
            endif
          endif
          
        else
          if ( ( tj .le. boundary_width ) .or. ( tj .gt. tny - boundary_width ) ) then

            if ( ( tk .le. boundary_width ) .or. ( tk .gt. tnz - boundary_width ) ) then
              mval= max(yval, zval)
              if (mval .ge. 0) N(i,j,k)%rock_state = rock_value
            else
              mval = yval
              if (mval .ge. 0) N(i,j,k)%rock_state = rock_value
            endif

          else
            if ( ( tk .le. boundary_width ) .or. ( tk .gt. tnz - boundary_width ) ) then
              mval = zval
              if (mval .ge. 0) N(i,j,k)%rock_state = rock_value
            else
              ! NOOP
            endif
          endif

        endif

      enddo
    enddo
  enddo

end subroutine lbe_modify_rock

!> reports whether stationary rock sites are present in the simulation
!>
!> \returns \c .true. if there is a geometry and \c .false. otherwise
logical function rock_is_present()
  select case (boundary_cond)
  case(1,2,4,6)
    rock_is_present = .true.
  case default
    rock_is_present = .false.
  end select

  rock_is_present = rock_is_present.or.trim(obs_file)/='empty.dat'
end function rock_is_present

subroutine lbe_init_rock_parallel(N)
  implicit none

  type(lbe_site),intent(inout) :: N(1-halo_extent:,1-halo_extent:,1-halo_extent:)

  integer :: rock_mpitype, rock_halo_extent

  rock_halo_extent = halo_extent

  write(msgstr,"('Creating MPI datatype for rock fields with halo_extent = ',I0)") rock_halo_extent
  call log_msg(msgstr)
  call build_rock_site_mpitype(rock_mpitype)
  call log_msg("Creating MPI datatypes for rock field chunks ...")
  call build_all_chunk_mpitypes(N, rock_halo, rock_mpitype, rock_halo_extent)
end subroutine lbe_init_rock_parallel

subroutine build_rock_site_mpitype(mpitype)
  implicit none

  integer, intent(out) :: mpitype ! type to be built

  integer, parameter :: cnt = 2 ! only rock_state & rock_colour

  type(lbe_site) :: sample
  integer(kind = MPI_ADDRESS_KIND) :: addrs(cnt), base, displs(cnt)
  integer :: blocklengths(cnt), types(cnt)
  integer :: mpierror

  ! Get base address of an lbe_site
  call MPI_Get_address(sample, base, mpierror)
  DEBUG_CHECKMPI(mpierror, "build_rock_site_mpitype, MPI_Get_address : base")

  ! rock_state
  blocklengths(1) = 1
  types(1) = LBE_REAL
  call MPI_Get_address(sample%rock_state,addrs(1), mpierror)
  DEBUG_CHECKMPI(mpierror, "build_rock_site_mpitype, MPI_Get_address : 1")

  ! rock_state
  blocklengths(2) = 1
  types(2) = LBE_REAL
  call MPI_Get_address(sample%rock_colour,addrs(2), mpierror)
  DEBUG_CHECKMPI(mpierror, "build_rock_site_mpitype, MPI_Get_address : 2")

  ! Calculate displacements
  displs(:) = addrs(:) - base

  call MPI_Type_create_struct(cnt, blocklengths, displs, types, mpitype, mpierror)
  DEBUG_CHECKMPI(mpierror, "build_rock_site_mpitype, MPI_Type_create_struct")
  call MPI_Type_commit(mpitype, mpierror)
  DEBUG_CHECKMPI(mpierror, "build_rock_site_mpitype, MPI_Type_commit")

end subroutine build_rock_site_mpitype

end module lbe_init_rock_module

