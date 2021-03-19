#include "lbe.h"

!> Lees-Edwards boundary conditions, or other boundary conditions
!> which involve changing the velocity and/or position of particles
!>
!> In the case of LE to produce shear.
!>
!>When special boundary conditions are in force, code here replaces
!>the halo exchange code in lbe_parallel_module.
module lbe_leesedwards_module
    use lbe_bdist_module
    use lbe_globals_module
    use lbe_log_module
    use lbe_parallel_module, only: ccoords,cdims,check_allocate,checkmpi&
         &,comm_cart,nnprocs,tag_mx,tag_my,tag_mz,tag_px,tag_py,tag_pz,tnx
    use lbe_parms_module, only: amass_b,amass_r,amass_s,omega_b,omega_r,omega_s&
         &,shear_omega,shear_u,inv_fluid,nt, nx, ny, nz, get_omega_r
#ifndef SINGLEFLUID
    use lbe_parms_module, only: get_omega_b
#endif
#ifndef NOSURFACTANT
    use lbe_parms_module, only: get_omega_s
#endif
    use lbe_types_module, only: lbe_site

    implicit none
    include 'mpif.h'

    private

    public le_adapt_halo_x,le_cleanup,le_fill_x_sbufs,le_find_bottom_neighbours&
         &,le_find_top_neighbours,le_fractional_offset,le_halo_exchange_x&
         &,le_halo_exchange_yz,le_init,le_init_step,le_shear_offset&
         &,le_shear_omega,le_shear_velocity
    public shear_sum
#ifndef NOSURFACTANT
    public le_adv_dipole_exchange,le_all_dipole_exchange
#endif
    !> \name x-recv buffers
    !>
    !> Between the calls to \c le_halo_exchange_x() and \c
    !> le_adapt_halo_x(), these buffers are available for adaption
    !> through plugins, such as MD. The three indexes are for all the
    !> reals communicated per lattice site (red, blue, surfactant
    !> populations plus dipole direction and rock state; as
    !> applicable), local y-position (1..ny+2), and local z-position
    !> (1..nz+3). The buffer-z-positions \c z+1 and \c z+2 might be
    !> used to interpolate the state at each local z-position \c z.
    !>
    !> \{
    real(kind=rk),allocatable,save,public :: le_prbuf(:,:,:) !< at position 0
    real(kind=rk),allocatable,save,public :: le_mrbuf(:,:,:) !< at position nx+1
    !> \}

    !> \name x-send buffers
    !>
    !> In case the populations should change at the edge of the system
    !> before \c le_halo_exchange_x() is executed (as they do in MD
    !> Ladd \c moving_wall_prebounce() ), it is possible to fill these
    !> buffers earlier calling \c le_fill_x_sbufs() . The three
    !> indexes are for all the reals communicated per lattice site
    !> (red, blue, surfactant populations plus dipole direction and
    !> rock state; as applicable), local y-position (1..ny+2), and
    !> local z-position (1..nz+3)
    !>
    !> \{
    real(kind=rk),allocatable,save :: le_psbuf(:,:,:) !< for position nx
    real(kind=rk),allocatable,save :: le_msbuf(:,:,:) !< for position 1
    !> \}

    !> Are \c le_psbuf and \c le_msbuf filled already?
    logical,save :: x_sbufs_filled=.false.

    integer :: lelo            ! Lattice offset
    real*8  :: lefrac          ! Fraction used in interpolation
    real*8  :: shear_sum = 0.d0

    !> holds the maximum amplitude of shear_u in case \c shear_omega/=0
    real(kind=rk),save :: shear_tmp

    !Array of ranks of processors will need to write to in LE scheme
    !Two arrays are needed as a processor can be both on the top and
    !bottom.
    integer, dimension(:), allocatable :: topprocs
    integer, dimension(:), allocatable :: botprocs

    ! Ranks of Lees Edwards neighbours (p|m,chunk1|chunk2)
    integer :: leprocs(2,3)=-1

    ! Unique tags to identify which part of a halo transfer
    ! data belongs to.
    integer, parameter :: tag_le_px1 = 9
    integer, parameter :: tag_le_px2 = 10
    integer, parameter :: tag_le_px3 = 11
    integer, parameter :: tag_le_px4 = 12
    integer, parameter :: tag_le_mx1 = 13
    integer, parameter :: tag_le_mx2 = 14
    integer, parameter :: tag_le_mx3 = 15
    integer, parameter :: tag_le_mx4 = 16

    ! number of real numbers communicated per lattice site
    integer,parameter :: ibuf&
#ifdef NOSURFACTANT
#ifdef SINGLEFLUID
         &=20
#else
         &=39
#endif
#else
         &=61
#endif

contains

    !> determines the two ranks located downward across the LE plane
    !> touching the local process domain for a given spatial z-offset
    !> across the LE plane
    !>
    !> \param[in] so spatial offset across the LE plane
    !>
    !> \param[out] loz neighbor rank touching the z-lower part of the
    !> local domain
    !>
    !> \param[out] hiz neighbor rank touching the z-upper part of the
    !> local domain (can be the same as \c loz)
    !>
    !> This routine must be called by a process that actually is located
    !> just above the LE plane.
    subroutine le_find_bottom_neighbours(so,loz,hiz)
        real(kind=rk),intent(in) :: so
        integer,intent(out) :: loz,hiz
        integer :: lower

        lower = ccoords(3)+floor(so/chunksize(3))

        ! different from the mod() intrinsic, the result of
        ! modulo(a,b) will always be in the range [0,b], even if a<0
        ! which is what we want here.
        loz = topprocs(modulo(lower,cdims(3))+1)
        hiz = topprocs(modulo(lower+1,cdims(3))+1)
    end subroutine le_find_bottom_neighbours

    !> determines the two ranks located upward across the LE plane
    !> touching the local process domain for a given spatial z-offset
    !> across the LE plane
    !>
    !> \param[in] so spatial offset across the LE plane
    !>
    !> \param[out] loz neighbor rank touching the z-lower part of the
    !> local domain
    !>
    !> \param[out] hiz neighbor rank touching the z-upper part of the
    !> local domain (can be the same as \c loz)
    !>
    !> This routine must be called by a process that actually is located
    !> just below the LE plane.
    subroutine le_find_top_neighbours(so,loz,hiz)
        real(kind=rk),intent(in) :: so
        integer,intent(out) :: loz,hiz
        integer :: lower

        lower = ccoords(3)-floor(so/chunksize(3))

        ! different from the mod() intrinsic, the result of
        ! modulo(a,b) will always be in the range [0,b], even if a<0
        ! which is what we want here.
        loz = botprocs(modulo(lower-1,cdims(3))+1)
        hiz = botprocs(modulo(lower,cdims(3))+1)
    end subroutine le_find_top_neighbours

    !> perform initial tasks for Lees Edwards code
    !>
    !> - setup shear velocity
    !> - allocate and fill topprocs and botprocs arrays holding ranks
    !>   of processors this rank will need to exchange information with.
    !> - find current communication partners
    subroutine le_init
        integer :: ierror
        integer :: i              ! looping variable

        ! Correct shear_u so that velocity on side planes in
        ! simple fluid tends to value in input file.
        shear_u = (tnx/(tnx-1.d0))*shear_u

        ! store this shear_u as the maximum amplitude of shear
        shear_tmp = shear_u

        ! If this processor is on minimum x-plane, remember the
        ! choice of top processors.
        if (ccoords(1) .eq. 0) then
           allocate(topprocs(cdims(3)))
           do i = 1, cdims(3)
              call MPI_Cart_Rank(Comm_Cart, &
                   (/ cdims(1)-1, ccoords(2), i-1 /), &
                   topprocs(i), ierror )
           enddo
        endif

        ! and/or if processor is on max x-plane remember bottom processor ranks
        if ((ccoords(1)+1) .eq. cdims(1)) then
           allocate(botprocs(cdims(3)))
           do i = 1, cdims(3)
              call MPI_Cart_Rank(Comm_Cart, &
                   (/ 0, ccoords(2), i-1 /), &
                   botprocs(i),ierror )
           enddo
        endif

        ! Work out how processors lie
        call le_neighbours()
    end subroutine le_init

    !> initialize LE boundary condition for the current time step
    subroutine le_init_step()
        shear_u = shear_tmp*cos(shear_omega*nt)
        if (myrankc==0.and.shear_omega>0.0_rk) then
           write(msgstr,"('SHEAR_U is set to ',ES15.8)") shear_u
           call log_msg(msgstr)
        endif

        shear_sum = shear_sum + shear_u

        call le_neighbours()
    end subroutine le_init_step

    !>Set up Lees Edwards for current timestep. leprocs stores the
    !>ranks of this processor's neighbours for interfaces where
    !>boundary conditions are strange or changing. It replaces nnprocs
    !>when not equal to -1.
    subroutine le_neighbours()
        integer :: proco
        real(kind=rk) :: offset

        select case (inv_fluid)
        case (5) ! Lees-Edwards
           offset = 2.0_rk*shear_sum

           lelo = mod(int(offset),nz)
           if (lelo .lt. 0) then
              write(msgstr,"('Lelo = ',I0,' , shear_sum = ',F16.10,"&
                   &//"' , shear_u = ',F16.10,' , shear_omega = ',F16.10)") &
                   &lelo, shear_sum, shear_u, shear_omega
              call log_msg(msgstr)
           end if

           lefrac = abs(offset - aint(offset))
           proco = int(offset/nz)

        case (6) ! rotated boundary conditions
           shear_u = 0.d0
           lelo = nz/2 + 1
           lefrac = 0.d0
           proco = cdims(3)/2
        end select

        ! If this processor is on minimum x-plane, find the
        ! two processors to write to.
        if (ccoords(1) .eq. 0) then
           !leprocs changed order, one added, as suggested by Jens Oct 2003
           leprocs(1,1) = topprocs(modulo(ccoords(3)+proco-1,cdims(3))+1)
           leprocs(1,2) = topprocs(modulo(ccoords(3)+proco  ,cdims(3))+1)
           leprocs(1,3) = topprocs(modulo(ccoords(3)+proco+1,cdims(3))+1)
        endif

        ! and/or if processor is on max x-plane...
        if (ccoords(1) .eq. cdims(1) - 1) then
           !leproc(2,3) added as suggested by Jens, Oct 2003
           leprocs(2,1) = botprocs(modulo(ccoords(3)-proco-1,cdims(3))+1)
           leprocs(2,2) = botprocs(modulo(ccoords(3)-proco  ,cdims(3))+1)
           leprocs(2,3) = botprocs(modulo(ccoords(3)-proco+1,cdims(3))+1)
        endif

#ifdef DEBUG_LE
        ! Who I am and where I'm going:
	print *, 'Rank ',myrankc,' writing to'&
            &,leprocs(1,1:2),leprocs(2,1:2),nnprocs(1,1:2)
#endif
    end subroutine le_neighbours

    !> fill the send buffers for the x halo exchange both for the LE
    !> planes and for the bulk and mark them as filled
    !>
    !> \param[in] whole_N local lattice chunk with halo of extent \c
    !> halo_extent
    subroutine le_fill_x_sbufs(whole_N)
        type(lbe_site),intent(in)&
             & :: whole_N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        integer :: j,k,stat

        ! prepare send from top buffer destined for bottom halo
	if (leprocs(2,1) .eq. -1) then
	   allocate (le_psbuf(ibuf,ny+2,nz+2),stat=stat)
           call check_allocate(stat,'le_fill_x_sbufs(): le_psbuf')
           do k=1,nz+2
              do j=1,ny+2
                 le_psbuf(1:19,j,k) = whole_N(nx,j-1,k-1)%n_r(:)
#ifndef SINGLEFLUID
                 le_psbuf(20:38,j,k) = whole_N(nx,j-1,k-1)%n_b(:)
#endif
#ifndef NOSURFACTANT
                 le_psbuf(39:57,j,k) = whole_N(nx,j-1,k-1)%n_s(:)
                 le_psbuf(58:60,j,k) = whole_N(nx,j-1,k-1)%d(:)
#endif
                 le_psbuf(ibuf,j,k) = whole_N(nx,j-1,k-1)%rock_state
              enddo
           enddo
	else
           allocate (le_psbuf(ibuf,ny+2,nz+3),stat=stat)
           call check_allocate(stat,'le_fill_x_sbufs(): le_psbuf')
           lelo1: if (lelo.lt.0) then
              do k=1,nz+1+lelo
                 do j=1,ny+2
                    le_psbuf(1:19,j,k)  = whole_N(nx,j-1,k)%n_r(:)
#ifndef SINGLEFLUID
                    le_psbuf(20:38,j,k) = whole_N(nx,j-1,k)%n_b(:)
#endif
#ifndef NOSURFACTANT
                    le_psbuf(39:57,j,k) = whole_N(nx,j-1,k)%n_s(:)
                    le_psbuf(58:60,j,k) = whole_N(nx,j-1,k)%d(:)
#endif
                    le_psbuf(ibuf,j,k) = whole_N(nx,j-1,k)%rock_state
                 enddo
              enddo

!        integer :: x,y,z
!        integer :: i,j,k,s          ! looping variables
              do k=nz+2+lelo,nz+3
                 do j=1,ny+2
                    le_psbuf(1:19,j,k)  = whole_N(nx,j-1,k-3)%n_r(:)
#ifndef SINGLEFLUID
                    le_psbuf(20:38,j,k) = whole_N(nx,j-1,k-3)%n_b(:)
#endif
#ifndef NOSURFACTANT
                    le_psbuf(39:57,j,k) = whole_N(nx,j-1,k-3)%n_s(:)
                    le_psbuf(58:60,j,k) = whole_N(nx,j-1,k-3)%d(:)
#endif
                    le_psbuf(ibuf,j,k) = whole_N(nx,j-1,k-3)%rock_state
                 enddo
              enddo
           else lelo1 ! lelo >= 0
              do k=1,lelo+2
                 do j=1,ny+2
                    le_psbuf(1:19,j,k)  = whole_N(nx,j-1,k)%n_r(:)
#ifndef SINGLEFLUID
                    le_psbuf(20:38,j,k) = whole_N(nx,j-1,k)%n_b(:)
#endif
#ifndef NOSURFACTANT
                    le_psbuf(39:57,j,k) = whole_N(nx,j-1,k)%n_s(:)
                    le_psbuf(58:60,j,k) = whole_N(nx,j-1,k)%d(:)
#endif
                    le_psbuf(ibuf,j,k) = whole_N(nx,j-1,k)%rock_state
                 enddo
              enddo

              do k=lelo+3,nz+3
                 do j=1,ny+2
                    le_psbuf(1:19,j,k)  = whole_N(nx,j-1,k-3)%n_r(:)
#ifndef SINGLEFLUID
                    le_psbuf(20:38,j,k) = whole_N(nx,j-1,k-3)%n_b(:)
#endif
#ifndef NOSURFACTANT
                    le_psbuf(39:57,j,k) = whole_N(nx,j-1,k-3)%n_s(:)
                    le_psbuf(58:60,j,k) = whole_N(nx,j-1,k-3)%d(:)
#endif
                    le_psbuf(ibuf,j,k) = whole_N(nx,j-1,k-3)%rock_state
                 enddo
              enddo
           end if lelo1
	endif

        ! prepare send from bottom buffer destined for top halo
	if (leprocs(1,2) .eq. -1) then
	   allocate (le_msbuf(ibuf,ny+2,nz+2),stat=stat)
           call check_allocate(stat,'le_fill_x_sbufs(): le_msbuf')
           do k=1,nz+2
              do j=1,ny+2
                 le_msbuf(1:19,j,k) = whole_N(1,j-1,k-1)%n_r(:)
#ifndef SINGLEFLUID
                 le_msbuf(20:38,j,k) = whole_N(1,j-1,k-1)%n_b(:)
#endif
#ifndef NOSURFACTANT
                 le_msbuf(39:57,j,k) = whole_N(1,j-1,k-1)%n_s(:)
                 le_msbuf(58:60,j,k) = whole_N(1,j-1,k-1)%d(:)
#endif
                 le_msbuf(ibuf,j,k) = whole_N(1,j-1,k-1)%rock_state
              enddo
           enddo
	else
           allocate (le_msbuf(ibuf,ny+2,nz+3),stat=stat)
           call check_allocate(stat,'le_fill_x_sbufs(): le_msbuf')
           lelo2: if (lelo.lt.0) then
              do k=1,-lelo+2
                 do j=1,ny+2
                    le_msbuf(1:19,j,k) = whole_N(1,j-1,k)%n_r(:)
#ifndef SINGLEFLUID
                    le_msbuf(20:38,j,k) = whole_N(1,j-1,k)%n_b(:)
#endif
#ifndef NOSURFACTANT
                    le_msbuf(39:57,j,k) = whole_N(1,j-1,k)%n_s(:)
                    le_msbuf(58:60,j,k) = whole_N(1,j-1,k)%d(:)
#endif
                    le_msbuf(ibuf,j,k) = whole_N(1,j-1,k)%rock_state
                 enddo
              enddo

              do k=3-lelo,nz+3
                 do j=1,ny+2
                    le_msbuf(1:19,j,k) = whole_N(1,j-1,k-3)%n_r(:)
#ifndef SINGLEFLUID
                    le_msbuf(20:38,j,k) = whole_N(1,j-1,k-3)%n_b(:)
#endif
#ifndef NOSURFACTANT
                    le_msbuf(39:57,j,k) = whole_N(1,j-1,k-3)%n_s(:)
                    le_msbuf(58:60,j,k) = whole_N(1,j-1,k-3)%d(:)
#endif
                    le_msbuf(ibuf,j,k) = whole_N(1,j-1,k-3)%rock_state
                 enddo
              enddo
           else lelo2 !lelo >= 0
              do k=1,nz+1-lelo
                 do j=1,ny+2
                    le_msbuf(1:19,j,k) = whole_N(1,j-1,k)%n_r(:)
#ifndef SINGLEFLUID
                    le_msbuf(20:38,j,k) = whole_N(1,j-1,k)%n_b(:)
#endif
#ifndef NOSURFACTANT
                    le_msbuf(39:57,j,k) = whole_N(1,j-1,k)%n_s(:)
                    le_msbuf(58:60,j,k) = whole_N(1,j-1,k)%d(:)
#endif
                    le_msbuf(ibuf,j,k) = whole_N(1,j-1,k)%rock_state
                 enddo
              enddo

              do k=nz+2-lelo,nz+3
                 do j=1,ny+2
                    le_msbuf(1:19,j,k) = whole_N(1,j-1,k-3)%n_r(:)
#ifndef SINGLEFLUID
                    le_msbuf(20:38,j,k) = whole_N(1,j-1,k-3)%n_b(:)
#endif
#ifndef NOSURFACTANT
                    le_msbuf(39:57,j,k) = whole_N(1,j-1,k-3)%n_s(:)
                    le_msbuf(58:60,j,k) = whole_N(1,j-1,k-3)%d(:)
#endif
                    le_msbuf(ibuf,j,k) = whole_N(1,j-1,k-3)%rock_state
                 enddo
              enddo
           end if lelo2
	endif

        x_sbufs_filled = .true.
    end subroutine le_fill_x_sbufs

    !> performs a halo exchange in x direction which on the maximum
    !> and minimum x-planes is subject to Lees-Edwards Boundary
    !> conditions
    !>
    !> The send buffers are filled only if this has not happened
    !> already.
    !>
    !> After this call, the z-shifted halo site data is available in
    !> \c le_mrbuf and \c le_prbuf for interpolation by \c
    !> le_adapt_halo_x(). So to complete the LE boundary update, \c
    !> le_adapt_halo_x() needs to be called, to be on the safe side,
    !> this should be done before the \c le_halo_exchange_yz() that
    !> follows.
    !>
    !> \param[in] whole_N local lattice chunk with halo of extent \c
    !> halo_extent
    subroutine le_halo_exchange_x(whole_N)
        type(lbe_site),intent(in)&
             & :: whole_N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        ! number of MPI requests. Away from the LE plane, a silly
        ! dummy integer is sent after the normal halo exchange, so
        ! there are always 2x2x2=8 requests for every process
        integer,parameter :: n_req=8
        integer :: requests(n_req)
        integer :: dummy,ierror,stat

!!!!!!
        ! Begin asynchronous receives
!!!!!!

	! Receive from top into bottom halo
	if (leprocs(1,2) .eq. -1) then
	   allocate (le_prbuf(ibuf,ny+2,nz+2),stat=stat)
           call check_allocate(stat,'le_halo_exchange_x(): le_prbuf')
	   call MPI_Irecv(	le_prbuf,		& ! buf
                ibuf*(ny+2)*(nz+2),	& ! count
                LBE_REAL,	& ! datatype
                nnprocs(1,1),	& ! source
                tag_px,		& ! tag
                Comm_Cart,	& ! communicator
                requests(1),	& ! request
                ierror)
	   call MPI_Irecv(	dummy,		& ! buf
                1,	& ! count
                MPI_INTEGER,	& ! datatype
                nnprocs(1,1),	& ! source
                1,		& ! tag
                Comm_Cart,	& ! communicator
                requests(2),	& ! request
                ierror)
	   call checkmpi(ierror,'Bad +x async receive')
	else
           allocate (le_prbuf(ibuf,ny+2,nz+3),stat=stat)
           call check_allocate(stat,'le_halo_exchange_x(): le_prbuf')
           if(lelo.lt.0)then
              call MPI_Irecv(	le_prbuf(1,1,3-lelo), &
                   ibuf*(ny+2)*(nz+lelo+1),          &
                   LBE_REAL,	& ! datatype
                   leprocs(1,2),	& ! source
                   tag_le_px3,	& ! tag
                   Comm_Cart,	& ! communicator
                   requests(1),	& ! request
                   ierror)
              call checkmpi(ierror,'Bad +x_le1 async receive')
              call MPI_Irecv(	le_prbuf(1,1,1), &
                   ibuf*(ny+2)*(-lelo+2),	   &
                   LBE_REAL,	& ! datatype
                   leprocs(1,1),	& ! source
                   tag_le_px4,	& ! tag
                   Comm_Cart,	& ! communicator
                   requests(2),	& ! request
                   ierror)

           else !lelo >=0
              call MPI_Irecv(	le_prbuf(1,1,nz-lelo+2), &
                   ibuf*(ny+2)*(lelo+2),          &
                   LBE_REAL,	& ! datatype
                   leprocs(1,3),	& ! source
                   tag_le_px1,	& ! tag
                   Comm_Cart,	& ! communicator
                   requests(1),	& ! request
                   ierror)
              call checkmpi(ierror,'Bad +x_le1 async receive')
              call MPI_Irecv(	le_prbuf(1,1,1), &
                   ibuf*(ny+2)*(nz-lelo+1),	   &
                   LBE_REAL,	& ! datatype
                   leprocs(1,2),	& ! source
                   tag_le_px2,	& ! tag
                   Comm_Cart,	& ! communicator
                   requests(2),	& ! request
                   ierror)

              call checkmpi(ierror,'Bad +x_le2 async receive')
              !! check this line
              !xreq = xreq + 1
              call checkmpi(ierror,'Bad +x_le2 async receive')
           endif !lelo >= or < 0
	endif

 ! Receive from bottom into top halo
	if (leprocs(2,1) .eq. -1) then
	   allocate (le_mrbuf(ibuf,ny+2,nz+2),stat=stat)
           call check_allocate(stat,'le_halo_exchange_x(): le_mrbuf')
	   call MPI_Irecv(	le_mrbuf,		& ! buf
                ibuf*(ny+2)*(nz+2),	& ! count
                LBE_REAL,	& ! datatype
                nnprocs(1,2),	& ! source
                tag_mx,		& ! tag
                Comm_Cart,	& ! communicator
                requests(3),	& ! request
                ierror)
	   call MPI_Irecv(	dummy,		& ! buf
                1,	& ! count
                MPI_INTEGER,	& ! datatype
                nnprocs(1,2),	& ! source
                2,		& ! tag
                Comm_Cart,	& ! communicator
                requests(4),	& ! request
                ierror)
	   call checkmpi(ierror,'Bad -x async receive')
	else
	   allocate (le_mrbuf(ibuf,ny+2,nz+3),stat=stat)
           call check_allocate(stat,'le_halo_exchange_x(): le_mrbuf')
           if(lelo.lt.0)then

              call MPI_Irecv(	le_mrbuf(1,1,nz+lelo+2),  &
                   ibuf*(ny+2)*(2-lelo),     &
                   LBE_REAL,	& ! datatype
                   leprocs(2,3),	& ! source
                   tag_le_mx3,	& ! tag
                   Comm_Cart,	& ! communicator
                   requests(3),	& ! request
                   ierror)
              call checkmpi(ierror,'Bad -x_le1 async receive')
              ! check this line
              call MPI_Irecv(	le_mrbuf(1,1,1), &
                   ibuf*(ny+2)*(nz+1+lelo),    &
                   LBE_REAL,	& ! datatype
                   leprocs(2,2),	& ! source
                   tag_le_mx4,	& ! tag
                   Comm_Cart,	& ! communicator
                   requests(4),	& ! request
                   ierror)

           else !lelo >= 0

              call MPI_Irecv(	le_mrbuf(1,1,lelo+3),  &
                   ibuf*(ny+2)*(nz+1-lelo),     &
                   LBE_REAL,	& ! datatype
                   leprocs(2,2),	& ! source
                   tag_le_mx1,	& ! tag
                   Comm_Cart,	& ! communicator
                   requests(3),	& ! request
                   ierror)
              call checkmpi(ierror,'Bad -x_le1 async receive')
              ! check this line
              !call MPI_Irecv(	minbuf(1,1,1), &
              call MPI_Irecv(	le_mrbuf(1,1,1), &
                   ibuf*(ny+2)*(lelo+2),    &
                   LBE_REAL,	& ! datatype
                   leprocs(2,1),	& ! source
                   tag_le_mx2,	& ! tag
                   Comm_Cart,	& ! communicator
                   requests(4),	& ! request
                   ierror)

              call checkmpi(ierror,'Bad -x_le2 async receive')
              ! check this line
              !xreq = xreq + 1
           endif !lelo < or >= 0
	endif

!!!!!!
 ! Begin asynchronous sends
!!!!!!

        if (.not.x_sbufs_filled) call le_fill_x_sbufs(whole_N)

        ! Send from top buffer destined for bottom halo
	if (leprocs(2,1) .eq. -1) then
	   call MPI_ISend(le_psbuf,	& ! buf
                ibuf*(ny+2)*(nz+2),	& ! count
                LBE_REAL,	& ! datatype
                nnprocs(1,2),	& ! dest
                tag_px,		& ! tag
                Comm_Cart,	& ! communicator
                requests(5),	& ! request
                ierror)
	   call MPI_ISend(	2,	& ! buf
                1,	& ! count
                MPI_INTEGER,	& ! datatype
                nnprocs(1,2),	& ! dest
                1,		& ! tag
                Comm_Cart,	& ! communicator
                requests(6),	& ! request
                ierror)
	   call checkmpi(ierror,'Bad +x async send')
	else
           lelo1: if (lelo.lt.0) then
              call MPI_ISend(le_psbuf(1,1,1), &
                   ibuf*(ny+2)*(nz+lelo+1),	   &
                   LBE_REAL,	& ! datatype
                   leprocs(2,2),	& ! dest
                   tag_le_px3,	& ! tag
                   Comm_Cart,	& ! communicator
                   requests(5),	& ! request
                   ierror)
              call checkmpi(ierror,'Bad +x async send')
              ! check this line
              !call MPI_ISend(	poutbuf(1,1,nz+2+lelo), &
              call MPI_ISend(le_psbuf(1,1,nz+2+lelo), &
                   ibuf*(ny+2)*(-lelo+2),     &
                   LBE_REAL,	& ! datatype
                   leprocs(2,3),	& ! dest
                   tag_le_px4,	& ! tag
                   Comm_Cart,	& ! communicator
                   requests(6),	& ! request
                   ierror)
           else lelo1 ! lelo >= 0
              call MPI_ISend(le_psbuf(1,1,1), &
                   ibuf*(ny+2)*(lelo+2),	   &
                   LBE_REAL,	& ! datatype
                   leprocs(2,1),	& ! dest
                   tag_le_px1,	& ! tag
                   Comm_Cart,	& ! communicator
                   requests(5),	& ! request
                   ierror)
              call checkmpi(ierror,'Bad +x async send')
              ! check this line
              !call MPI_ISend(	poutbuf(1,1,lelo+3), &
              call MPI_ISend(le_psbuf(1,1,lelo+3), &
                   ibuf*(ny+2)*(nz-lelo+1),     &
                   LBE_REAL,	& ! datatype
                   leprocs(2,2),	& ! dest
                   tag_le_px2,	& ! tag
                   Comm_Cart,	& ! communicator
                   requests(6),	& ! request
                   ierror)
              call checkmpi(ierror,'Bad +x async send')
!check this line
 !             xreq = xreq + 1
  !         endif ! lelo>= or < 0
           end if lelo1
	endif

        ! Send from bottom buffer destined for top halo
	if (leprocs(1,2) .eq. -1) then
	   call MPI_ISend(le_msbuf,	& ! buf
                ibuf*(ny+2)*(nz+2),	& ! count
                LBE_REAL,	& ! datatype
                nnprocs(1,1),	& ! dest
                tag_mx,		& ! tag
                Comm_Cart,	& ! communicator
                requests(7),	& ! request
                ierror)
	   call MPI_ISend(	2,	& ! buf
                1,	& ! count
                MPI_INTEGER,	& ! datatype
                nnprocs(1,1),	& ! dest
                2,		& ! tag
                Comm_Cart,	& ! communicator
                requests(8),	& ! request
                ierror)
	   call checkmpi(ierror,'Bad -x async send')
	else
           lelo2: if(lelo.lt.0) then
              call MPI_ISend(le_msbuf(1,1,1), &
                   ibuf*(ny+2)*(2-lelo),     &
                   LBE_REAL,	& ! datatype
                   leprocs(1,1),	& ! dest
                   tag_le_mx3,	& ! tag
                   Comm_Cart,	& ! communicator
                   requests(7),	& ! request
                   ierror)
              call checkmpi(ierror,'Bad +x async send')
              ! check this line
              !call MPI_ISend(	moutbuf(1,1,3-lelo), &
              call MPI_ISend(le_msbuf(1,1,3-lelo), &
                   ibuf*(ny+2)*(nz+1+lelo), &
                   LBE_REAL,	& ! datatype
                   leprocs(1,2),	& ! dest
                   tag_le_mx4,	& ! tag
                   Comm_Cart,	& ! communicator
                   requests(8),	& ! request
                   ierror)
           else lelo2 ! lelo >= 0
              call MPI_ISend(le_msbuf(1,1,1), &
                   ibuf*(ny+2)*(nz+1-lelo),     &
                   LBE_REAL,	& ! datatype
                   leprocs(1,2),	& ! dest
                   tag_le_mx1,	& ! tag
                   Comm_Cart,	& ! communicator
                   requests(7),	& ! request
                   ierror)
              call checkmpi(ierror,'Bad +x async send')
              !check this line
              !call MPI_ISend(	moutbuf(1,1,nz+2-lelo), &
              !=======
              call MPI_ISend(le_msbuf(1,1,nz+2-lelo), &
                   ibuf*(ny+2)*(lelo+2), &
                   LBE_REAL,	& ! datatype
                   leprocs(1,3),	& ! dest
                   tag_le_mx2,	& ! tag
                   Comm_Cart,	& ! communicator
                   requests(8),	& ! request
                   ierror)
           end if lelo2
	endif

        ! wait for all I/O to complete
	call MPI_Waitall(n_req,requests,MPI_STATUSES_IGNORE,ierror)
	call checkmpi(ierror,'MPI_Waitall() failed in X direction')

        x_sbufs_filled = .false.
	if (allocated(le_psbuf)) deallocate(le_psbuf)
	if (allocated(le_msbuf)) deallocate(le_msbuf)
    end subroutine le_halo_exchange_x

    !check this line
    !	call checkmpi(ierror,'MPI_Waitall() failed in X direction')
    !> where applicable, LE interpolate halo sites and apply LE
    !> velocity offset
    !>
    !> For interpolation, \c le_mrbuf and \c le_prbuf are taken as
    !> input.
    !>
    !> \param[in,out] N local chunk of the lattice with halo extent 1
    !> (old LB3D style)
    subroutine le_adapt_halo_x(N)
        type(lbe_site),dimension(0:,0:,0:),intent(inout) :: N
        integer :: y,z,j,k,s    ! looping variables
        real*8,dimension(nvecs) :: boltz_du, boltz_u, boltz_u_du
        real*8,dimension(3) :: p_tilde, u_tilde
        real*8 :: rho_tilde, rho_tmp

        real*8 :: rho_r,rN_s,p_r(3),df_r(nvecs),fnew_r(nvecs)
#ifndef SINGLEFLUID
        real*8 :: rho_b,bN_s,p_b(3),df_b(nvecs),fnew_b(nvecs)
#ifndef NOSURFACTANT
        real*8 :: rho_s,sN_s,p_s(3),df_s(nvecs),fnew_s(nvecs)
#endif
#endif
        real*8 :: df_scale

 ! Write bottom halo

	if (leprocs(1,2) .ne. -1) then
    ! FIXME: So far I've failed to get a terse array version of this
    ! with an overloaded * operator working
	   do y=0,ny+1
	      do z=0,nz+1
		 N(0,y,z)%n_r(:) = (1-lefrac)*le_prbuf(1:19,y+1,z+1) + &
                      lefrac*le_prbuf(1:19,y+1,z+2)
#ifndef SINGLEFLUID
		 N(0,y,z)%n_b(:) = (1-lefrac)*le_prbuf(20:38,y+1,z+1) + &
                      lefrac*le_prbuf(20:38,y+1,z+2)
#endif
#ifndef NOSURFACTANT
		 N(0,y,z)%n_s(:) = (1-lefrac)*le_prbuf(39:57,y+1,z+1) + &
                      lefrac*le_prbuf(39:57,y+1,z+2)
		 N(0,y,z)%d(:) = (1-lefrac)*le_prbuf(58:60,y+1,z+1) +  &
                      lefrac*le_prbuf(58:60,y+1,z+2)
#endif

   !d ----- under construction - change velocity at boundary ------------

   ! bN_s = density of blue particles at this site, etc.
   ! p_b = velocity of blue particles at this site, etc.

		 p_r = 0.d0
		 rN_s = 0.d0
#ifndef SINGLEFLUID
		 p_b = 0.d0
		 bN_s = 0.d0
#ifndef NOSURFACTANT
		 p_s = 0.d0
		 sN_s = 0.d0
#endif
#endif

		 do s=1,nvecs
		    rN_s = rN_s + N(0,y,z)%n_r(s) * g(s)
		    p_r = p_r + N(0,y,z)%n_r(s)*g(s)*c(s,:)
#ifndef SINGLEFLUID
		    p_b = p_b + N(0,y,z)%n_b(s)*g(s)*c(s,:)
		    bN_s = bN_s + N(0,y,z)%n_b(s) * g(s)
#endif
#ifndef NOSURFACTANT
		    sN_s = sN_s + N(0,y,z)%n_s(s) * g(s)
		    p_s = p_s + N(0,y,z)%n_s(s)*g(s)*c(s,:)
#endif
		 end do

   ! Calculate weighted total momentum  p_tilde ,
   ! rho_i = mass of species i at this site, and
   ! rho_tilde = sum of masses over relaxation times.

                 p_tilde = amass_r*get_omega_r(N, 0, y, z)*p_r
		 rho_r = amass_r*rN_s
                 rho_tilde = rho_r*get_omega_r(N, 0, y, z)
#ifndef SINGLEFLUID
                 p_tilde = p_tilde + amass_b*get_omega_b(N, 0, y, z)*p_b
		 rho_b = amass_b*bN_s
                 rho_tilde = rho_tilde + rho_b*get_omega_b(N, 0, y, z)
#ifndef NOSURFACTANT
                 p_tilde = p_tilde + amass_s*get_omega_s(N, 0, y, z)*p_s
		 rho_s = amass_s*sN_s
                 rho_tilde = rho_tilde + rho_s*get_omega_s(N, 0, y, z)
#endif
#endif

		 ! rho_tmp = clipped version of this.

		 rho_tmp = max(rho_tilde,dble(10.e-9))

   ! Calculate averaged velocity

		 u_tilde = p_tilde / rho_tmp


   ! FIXME This should be optimized.

		 call boltz_dist(u_tilde(1),u_tilde(2),&
                      &u_tilde(3)-shear_u-shear_u,&
                      &0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,boltz_u_du)

		 call boltz_dist(u_tilde(1),u_tilde(2),u_tilde(3),&
                      &0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,boltz_u)

		 boltz_du = boltz_u_du - boltz_u

#ifdef DEBUG_LE
                 !Should be zero (or very little) change in flux in x
                 !direction
                 print *,sum(boltz_du*cx*g)
#endif
		 df_r = rN_s * boltz_du
		 fnew_r = N(0,y,z)%n_r + df_r
#ifndef SINGLEFLUID
		 df_b = bN_s * boltz_du
		 fnew_b = N(0,y,z)%n_b + df_b
#ifndef NOSURFACTANT
		 df_s = sN_s * boltz_du
		 fnew_s = N(0,y,z)%n_s + df_s
#endif
#endif

   ! If the velocity adjustment would make the site density
   ! negative anywhere, then scale down the velocity adjustment
   ! so that it will not
		 if (minval(fnew_r)<0.0_8&
#ifndef SINGLEFLUID
                    &.or.minval(fnew_b)<0.0_8&
#ifndef NOSURFACTANT
                    &.or.minval(fnew_s)<0.0_8&
#endif
#endif
                    &) then

      ! Avoid dividing by zero
                    if (minval((/N(0,y,z)%n_r&
#ifndef SINGLEFLUID
                       &,N(0,y,z)%n_b&
#ifndef NOSURFACTANT
                       &,N(0,y,z)%n_s&
#endif
#endif
                       &/))>0.0_8) then
                       df_scale = maxval(-(/df_r&
#ifndef SINGLEFLUID
                       &,df_b&
#ifndef NOSURFACTANT
                       &,df_s&
#endif
#endif
                       &/)/(/N(0,y,z)%n_r&
#ifndef SINGLEFLUID
                       &,N(0,y,z)%n_b&
#ifndef NOSURFACTANT
                       &,N(0,y,z)%n_s&
#endif
#endif
                       &/))
                    else
                       df_scale = 0.0_8
                    end if
                    if (abs(df_scale).gt.1.0_8) then
                       fnew_r = N(0,y,z)%n_r + df_r/df_scale
#ifndef SINGLEFLUID
                       fnew_b = N(0,y,z)%n_b + df_b/df_scale
#ifndef NOSURFACTANT
                       fnew_s = N(0,y,z)%n_s + df_s/df_scale
#endif
#endif
                    else
                       fnew_r = N(0,y,z)%n_r
#ifndef SINGLEFLUID
                       fnew_b = N(0,y,z)%n_b
#ifndef NOSURFACTANT
                       fnew_s = N(0,y,z)%n_s
#endif
#endif
                    endif
                 endif
                 N(0,y,z)%n_r = fnew_r
#ifndef SINGLEFLUID
                 N(0,y,z)%n_b = fnew_b
#ifndef NOSURFACTANT
                 N(0,y,z)%n_s = fnew_s
#endif
#endif
              enddo
           enddo
	else
           do k=1,nz+2
              do j=1,ny+2
                 N(0,j-1,k-1)%n_r(:) = le_prbuf(1:19,j,k)
#ifndef SINGLEFLUID
                 N(0,j-1,k-1)%n_b(:) = le_prbuf(20:38,j,k)
#endif
#ifndef NOSURFACTANT
                 N(0,j-1,k-1)%n_s(:) = le_prbuf(39:57,j,k)
                 N(0,j-1,k-1)%d(:) = le_prbuf(58:60,j,k)
#endif
                 N(0,j-1,k-1)%rock_state = le_prbuf(ibuf,j,k)
              enddo
           enddo
	endif

 ! Write top halo
	if (leprocs(2,1) .ne. -1) then
    ! FIXME: So far I've failed to get a terse array version of this
    ! with an overloaded * operator working
	   do y=0,ny+1
	      do z=0,nz+1
		 N(nx+1,y,z)%n_r(:) = lefrac*le_mrbuf(1:19,y+1,z+1) &
                      &+ (1-lefrac)*le_mrbuf(1:19,y+1,z+2)
#ifndef SINGLEFLUID
		 N(nx+1,y,z)%n_b(:) = lefrac*le_mrbuf(20:38,y+1,z+1) &
                      &+ (1-lefrac)*le_mrbuf(20:38,y+1,z+2)
#ifndef NOSURFACTANT
		 N(nx+1,y,z)%n_s(:) = lefrac*le_mrbuf(39:57,y+1,z+1) &
                      &+ (1-lefrac)*le_mrbuf(39:57,y+1,z+2)
		 N(nx+1,y,z)%d(:) = lefrac*le_mrbuf(58:60,y+1,z+1) &
                      &+ (1-lefrac)*le_mrbuf(58:60,y+1,z+2)
#endif
#endif

   ! bN_s = density of blue particles at this site, etc.
   ! p_b = velocity of blue particles at this site, etc.

                 p_r = 0.d0
                 rN_s = 0.d0
#ifndef SINGLEFLUID
                 p_b = 0.d0
                 bN_s = 0.d0
#ifndef NOSURFACTANT
                 p_s = 0.d0
                 sN_s = 0.d0
#endif
#endif
                 do s=1,nvecs
                    rN_s = rN_s + N(nx+1,y,z)%n_r(s) * g(s)
                    p_r = p_r + N(nx+1,y,z)%n_r(s)*g(s)*c(s,:)
#ifndef SINGLEFLUID
                    bN_s = bN_s + N(nx+1,y,z)%n_b(s) * g(s)
                    p_b = p_b + N(nx+1,y,z)%n_b(s)*g(s)*c(s,:)
#ifndef NOSURFACTANT
                    sN_s = sN_s + N(nx+1,y,z)%n_s(s) * g(s)
                    p_s = p_s + N(nx+1,y,z)%n_s(s)*g(s)*c(s,:)
#endif
#endif
                 end do

                 ! Calculate weighted total momentum
                 ! rho_i = mass of species i at this site.
                 ! rho_tilde = sum of masses over relaxation times.
                 p_tilde = amass_r*get_omega_r(N, nx+1, y, z)*p_r
                 rho_r = amass_r*rN_s
                 rho_tilde = rho_r*get_omega_r(N, nx+1, y, z)
#ifndef SINGLEFLUID
                 p_tilde = p_tilde + amass_b*get_omega_b(N, nx+1, y, z)*p_b
                 rho_b = amass_b*bN_s
                 rho_tilde = rho_tilde + rho_b*get_omega_b(N, nx+1, y, z)
#ifndef NOSURFACTANT
                 p_tilde = p_tilde + amass_s*get_omega_s(N, nx+1, y, z)*p_s
                 rho_s = amass_s*sN_s
                 rho_tilde = rho_tilde + rho_s*get_omega_s(N, nx+1, y, z)
#endif
#endif

                 ! rho_tmp = clipped version of this.
                 rho_tmp = max(rho_tilde,dble(10.e-9))

                 ! Calculate averaged velocity
                 u_tilde = p_tilde / rho_tmp

		 ! FIXME This should be optimized.

		 call boltz_dist(u_tilde(1),u_tilde(2)&
                      &,u_tilde(3)+shear_u+shear_u&
                      &,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,boltz_u_du)
		 call boltz_dist(u_tilde(1),u_tilde(2),u_tilde(3)&
                      &,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,boltz_u)

		 boltz_du = boltz_u_du - boltz_u

		 df_r = rN_s * boltz_du
		 fnew_r = N(nx+1,y,z)%n_r + df_r
#ifndef SINGLEFLUID
		 df_b = bN_s * boltz_du
		 fnew_b = N(nx+1,y,z)%n_b + df_b
#ifndef NOSURFACTANT
		 df_s = sN_s * boltz_du
		 fnew_s = N(nx+1,y,z)%n_s + df_s
#endif
#endif

   ! If the velocity adjustment would make the site density
   ! negative anywhere, then scale down the velocity adjustment
   ! so that it will not
		 if (minval(fnew_r)<0.0_8&
#ifndef SINGLEFLUID
                    &.or.minval(fnew_b)<0.0_8&
#ifndef NOSURFACTANT
                    &.or.minval(fnew_s)<0.0_8&
#endif
#endif
                    &) then

      ! Avoid dividing by zero
                    if (minval((/N(nx+1,y,z)%n_r&
#ifndef SINGLEFLUID
                       &,N(nx+1,y,z)%n_b&
#ifndef NOSURFACTANT
                       &,N(nx+1,y,z)%n_s&
#endif
#endif
                       &/))>0.0_8) then
                       df_scale = maxval(-(/df_r&
#ifndef SINGLEFLUID
                       &,df_b&
#ifndef NOSURFACTANT
                       &,df_s&
#endif
#endif
                       &/)/(/N(nx+1,y,z)%n_r&
#ifndef SINGLEFLUID
                       &,N(nx+1,y,z)%n_b&
#ifndef NOSURFACTANT
                       &,N(nx+1,y,z)%n_s&
#endif
#endif
                       &/))
                    else
                       df_scale = 0.0_8
                    end if
                    if (abs(df_scale).gt.1.0_8) then
                       fnew_r = N(nx+1,y,z)%n_r + df_r/df_scale
#ifndef SINGLEFLUID
                       fnew_b = N(nx+1,y,z)%n_b + df_b/df_scale
#ifndef NOSURFACTANT
                       fnew_s = N(nx+1,y,z)%n_s + df_s/df_scale
#endif
#endif
                    else
                       fnew_r = N(nx+1,y,z)%n_r
#ifndef SINGLEFLUID
                       fnew_b = N(nx+1,y,z)%n_b
#ifndef NOSURFACTANT
                       fnew_s = N(nx+1,y,z)%n_s
#endif
#endif
                    endif
                 endif

                 N(nx+1,y,z)%n_r = fnew_r
#ifndef SINGLEFLUID
                 N(nx+1,y,z)%n_b = fnew_b
#ifndef NOSURFACTANT
                 N(nx+1,y,z)%n_s = fnew_s
#endif
#endif
              enddo
           enddo
        else
           do k=1,nz+2
              do j=1,ny+2
                 N(nx+1,j-1,k-1)%n_r = le_mrbuf(1:19,j,k)
#ifndef SINGLEFLUID
                 N(nx+1,j-1,k-1)%n_b = le_mrbuf(20:38,j,k)
#endif
#ifndef NOSURFACTANT
                 N(nx+1,j-1,k-1)%n_s = le_mrbuf(39:57,j,k)
                 N(nx+1,j-1,k-1)%d = le_mrbuf(58:60,j,k)
#endif
                 N(nx+1,j-1,k-1)%rock_state = le_mrbuf(ibuf,j,k)
              enddo
           enddo
        endif

#ifdef DEBUG_LE
	do s=1,nvecs
	   if (minval( N(nx+1,:,:)%n_r(s) ) .lt. 0 ) then
	      print *,'Ack! Negative density alert stage'
	   endif
	   if (minval( N(0,:,:)%n_r(s) ) .lt. 0 ) then
	      print *,'Ack! Negative density alert stage'
	   endif
	enddo
#endif

	if (allocated(le_prbuf)) deallocate(le_prbuf)
	if (allocated(le_mrbuf)) deallocate(le_mrbuf)
    end subroutine le_adapt_halo_x

    !> sub-lattice part of the z-offset between both Lees-Edwards planes
    !>
    !> \warning This is valid only after \c le_init_step() was called
    !> until the end of the LB step.
    pure function le_fractional_offset()
        real(kind=rk) :: le_fractional_offset

        le_fractional_offset = lefrac
    end function le_fractional_offset

    !> Performs a halo exchange in y and z direction after a
    !> Lees-Edwards halo exchange was already performed in the
    !> x-direction
    !>
    !> \param[in,out] N local chunk of the lattice with halo extent 1
    !> (old LB3D style)
    subroutine le_halo_exchange_yz(N)
        type(lbe_site),dimension(0:,0:,0:),intent(inout) :: N
        real*8,dimension(:,:,:),allocatable :: pinbuf,poutbuf,minbuf,moutbuf
        integer,parameter :: pin=1,mmin=2,pout=3,mout=4 !defined request numbers
        integer,parameter :: n_req=4
        integer :: requests(n_req)
        integer :: ierror,stat
        integer :: i,j,k          ! looping variables

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Swap in the Y direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	allocate(pinbuf(ibuf,nx+2,nz+2),poutbuf(ibuf,nx+2,nz+2)&
             &,moutbuf(ibuf,nx+2,nz+2),minbuf(ibuf,nx+2,nz+2),stat=stat)
        call check_allocate(stat&
             &,'le_halo_exchange_yz(): pinbuf,poutbuf,moutbuf,minbuf')

!!!!!!
 ! Begin asynchronous receives
!!!!!!

	call MPI_Irecv(		pinbuf,		& ! buf
             ibuf*(nx+2)*(nz+2),	& ! count
             LBE_REAL,	& ! datatype
             nnprocs(2,1),	& ! source
             tag_py,		& ! tag
             Comm_Cart,	& ! communicator
             requests(pin),	& ! request
             ierror)
	call checkmpi(ierror,'Bad +y async receive')

	call MPI_Irecv(		minbuf,		& ! buf
             ibuf*(nx+2)*(nz+2),	& ! count
             LBE_REAL,	& ! datatype
             nnprocs(2,2),	& ! source
             tag_my,		& ! tag
             Comm_Cart,	& ! communicator
             requests(mmin),	& ! request
             ierror)
	call checkmpi(ierror,'Bad -y async receive')

        do k=1,nz+2
           do i=1,nx+2
              poutbuf(1:19,i,k) = N(i-1,ny,k-1)%n_r
#ifndef SINGLEFLUID
              poutbuf(20:38,i,k) = N(i-1,ny,k-1)%n_b
#endif
#ifndef NOSURFACTANT
              poutbuf(39:57,i,k) = N(i-1,ny,k-1)%n_s
              poutbuf(58:60,i,k) = N(i-1,ny,k-1)%d
#endif
              poutbuf(ibuf,i,k) = N(i-1,ny,k-1)%rock_state

              moutbuf(1:19,i,k) = N(i-1,1,k-1)%n_r
#ifndef SINGLEFLUID
              moutbuf(20:38,i,k) = N(i-1,1,k-1)%n_b
#endif
#ifndef NOSURFACTANT
              moutbuf(39:57,i,k) = N(i-1,1,k-1)%n_s
              moutbuf(58:60,i,k) = N(i-1,1,k-1)%d
#endif
              moutbuf(ibuf,i,k) = N(i-1,1,k-1)%rock_state
           enddo
        enddo

!!!!!!
	! Begin asynchronous sends
!!!!!!

	call MPI_ISend(		poutbuf,	& ! buf
             ibuf*(nx+2)*(nz+2),	& ! count
             LBE_REAL,	& ! datatype
             nnprocs(2,2),	& ! dest
             tag_py,		& ! tag
             Comm_Cart,	& ! communicator
             requests(pout),	& ! request
             ierror)
	call checkmpi(ierror,'Bad +y async send')

	call MPI_ISend(		moutbuf,	& ! buf
             ibuf*(nx+2)*(nz+2),	& ! count
             LBE_REAL,	& ! datatype
             nnprocs(2,1),	& ! dest
             tag_my,		& ! tag
             Comm_Cart,	& ! communicator
             requests(mout),	& ! request
             ierror)
	call checkmpi(ierror,'Bad -y async send')

!!!!!!
 ! Now wait for all I/O to complete
!!!!!!

        ! check this line
        ! call MPI_Waitall(4,requests,statuses,ierror)
	call MPI_Waitall(4,requests,MPI_STATUSES_IGNORE,ierror)
	call checkmpi(ierror,'MPI_Waitall() failed in Y direction')

        do k=1,nz+2
           do i=1,nx+2
              N(i-1,ny+1,k-1)%n_r = minbuf(1:19,i,k)
#ifndef SINGLEFLUID
              N(i-1,ny+1,k-1)%n_b = minbuf(20:38,i,k)
#endif
#ifndef NOSURFACTANT
              N(i-1,ny+1,k-1)%n_s = minbuf(39:57,i,k)
              N(i-1,ny+1,k-1)%d = minbuf(58:60,i,k)
#endif
              N(i-1,ny+1,k-1)%rock_state = minbuf(ibuf,i,k)

              N(i-1,0,k-1)%n_r = pinbuf(1:19,i,k)
#ifndef SINGLEFLUID
              N(i-1,0,k-1)%n_b = pinbuf(20:38,i,k)
#endif
#ifndef NOSURFACTANT
              N(i-1,0,k-1)%n_s = pinbuf(39:57,i,k)
              N(i-1,0,k-1)%d = pinbuf(58:60,i,k)
#endif
              N(i-1,0,k-1)%rock_state = pinbuf(ibuf,i,k)
           enddo
        enddo

	deallocate(pinbuf,poutbuf,minbuf,moutbuf)

!!!!
 ! Y swap done.
!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Swap in the Z direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        allocate (pinbuf(ibuf,nx+2,ny+2),poutbuf(ibuf,nx+2,ny+2)&
             &,moutbuf(ibuf,nx+2,ny+2),minbuf(ibuf,nx+2,ny+2),stat=stat)
        call check_allocate(stat&
             &,'le_halo_exchange_yz(): pinbuf,poutbuf,moutbuf,minbuf')

!!!!!!
 ! Begin asynchronous receives
!!!!!!

	call MPI_Irecv(		pinbuf,		& ! buf
             ibuf*(nx+2)*(ny+2),	& ! count
             LBE_REAL,	& ! datatype
             nnprocs(3,1),	& ! source
             tag_pz,		& ! tag
             Comm_Cart,	& ! communicator
             requests(pin),	& ! request
             ierror)
	call checkmpi(ierror,'Bad +z async receive')

	call MPI_Irecv(		minbuf,		& ! buf
             ibuf*(nx+2)*(ny+2),	& ! count
             LBE_REAL,	& ! datatype
             nnprocs(3,2),	& ! source
             tag_mz,		& ! tag
             Comm_Cart,	& ! communicator
             requests(mmin),	& ! request
             ierror)
	call checkmpi(ierror,'Bad -z async receive')

        do j=1,ny+2
           do i=1,nx+2
              poutbuf(1:19,i,j) = N(i-1,j-1,nz)%n_r
#ifndef SINGLEFLUID
              poutbuf(20:38,i,j) = N(i-1,j-1,nz)%n_b
#endif
#ifndef NOSURFACTANT
              poutbuf(39:57,i,j) = N(i-1,j-1,nz)%n_s
              poutbuf(58:60,i,j) = N(i-1,j-1,nz)%d
#endif
              poutbuf(ibuf,i,j) = N(i-1,j-1,nz)%rock_state

              moutbuf(1:19,i,j) = N(i-1,j-1,1)%n_r
#ifndef SINGLEFLUID
              moutbuf(20:38,i,j) = N(i-1,j-1,1)%n_b
#endif
#ifndef NOSURFACTANT
              moutbuf(39:57,i,j) = N(i-1,j-1,1)%n_s
              moutbuf(58:60,i,j) = N(i-1,j-1,1)%d
#endif
              moutbuf(ibuf,i,j) = N(i-1,j-1,1)%rock_state
           enddo
        enddo

!!!!!!
	! Begin asynchronous sends
!!!!!!

	call MPI_ISend(		poutbuf,	& ! buf
             ibuf*(nx+2)*(ny+2),	& ! count
             LBE_REAL,	& ! datatype
             nnprocs(3,2),	& ! dest
             tag_pz,		& ! tag
             Comm_Cart,	& ! communicator
             requests(pout),	& ! request
             ierror)
	call checkmpi(ierror,'Bad +z async send')

	call MPI_ISend(		moutbuf,	& ! buf
             ibuf*(nx+2)*(ny+2),	& ! count
             LBE_REAL,	& ! datatype
             nnprocs(3,1),	& ! dest
             tag_mz,		& ! tag
             Comm_Cart,	& ! communicator
             requests(mout),	& ! request
             ierror)
	call checkmpi(ierror,'Bad -z async send')

!!!!!!
 ! Now wait for all I/O to complete
!!!!!!

        !check this line
        !	call MPI_Waitall(4,requests,statuses,ierror)
        !=======
	call MPI_Waitall(4,requests,MPI_STATUSES_IGNORE,ierror)
	call checkmpi(ierror,'MPI_Waitall() failed in Z direction')

        do j=1,ny+2
           do i=1,nx+2
              N(i-1,j-1,nz+1)%n_r = minbuf(1:19,i,j)
#ifndef SINGLEFLUID
              N(i-1,j-1,nz+1)%n_b = minbuf(20:38,i,j)
#endif
#ifndef NOSURFACTANT
              N(i-1,j-1,nz+1)%n_s = minbuf(39:57,i,j)
              N(i-1,j-1,nz+1)%d   = minbuf(58:60,i,j)
#endif
              N(i-1,j-1,nz+1)%rock_state = minbuf(ibuf,i,j)

              N(i-1,j-1,0)%n_r = pinbuf(1:19,i,j)
#ifndef SINGLEFLUID
              N(i-1,j-1,0)%n_b = pinbuf(20:38,i,j)
#endif
#ifndef NOSURFACTANT
              N(i-1,j-1,0)%n_s = pinbuf(39:57,i,j)
              N(i-1,j-1,0)%d   = pinbuf(58:60,i,j)
#endif
              N(i-1,j-1,0)%rock_state = pinbuf(ibuf,i,j)
           enddo
        enddo

	deallocate(pinbuf,poutbuf,minbuf,moutbuf)
!!!!
 ! Z swap done.
!!!!
    end subroutine le_halo_exchange_yz

#ifndef NOSURFACTANT
!> Like a halo exchange, but just swap the advected dipole moment vectors
!> according to the diplacements of Lees-Edwards boundary conditions.
subroutine le_adv_dipole_exchange(d_adv)

	implicit none

	real*8,dimension(:,:,:,:),allocatable :: pinbuf
	real*8,dimension(:,:,:,:),allocatable :: poutbuf
	real*8,dimension(:,:,:,:),allocatable :: minbuf
	real*8,dimension(:,:,:,:),allocatable :: moutbuf
	real*8,dimension(1:,1:,0:,0:,0:) :: d_adv

	integer, parameter :: pin=1,mmin=2,pout=3,mout=4
	! When more complicated boundary conditions became possible
	! more (and an uncertain number of) request numbers were needed.
	! This is coped by xreq, the number of the next undefined request
	! All request numbers could be replaced with this, but it
	! may hamper performance or upset someone's feelings?
	integer :: xreq=5

!	integer, dimension(4) :: requests
	! FIXME: Should allocate this properly
	! - as it is we'll not get more than 8 requests
	integer, dimension(10) :: requests
	integer, dimension(MPI_STATUS_SIZE,10) :: statuses
	integer :: ierror

	integer :: x,y,z,s,i
	real*8 :: t1,t2,t3

	!print*,'Processor ',myrankc,'beginning exchange'

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Swap in the X direction
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	xreq = 5

	!!!!!!
	! Begin asynchronous receives
	!!!!!!

	! FIXME if the lattice changes
	! below 57 = 3*19 = nd*nvecs


	! Receive from top into bottom halo
	if (leprocs(1,2) .eq. -1) then
	   ! Why change to zero based now when one based in halo_exchange? Grrr
	   allocate(pinbuf(nd,nvecs,0:ny+1,0:nz+1))
	   call MPI_Irecv(	pinbuf,		& ! buf
				(ny+2)*(nz+2)*57,&! count
				LBE_REAL,	& ! datatype
				nnprocs(1,1),	& ! source
				tag_px,		& ! tag
				Comm_Cart,	& ! communicator
				requests(1),	& ! request
				ierror)
	   call MPI_Irecv(	xreq,		& ! buf
				1,	& ! count
				MPI_INTEGER,	& ! datatype
				nnprocs(1,1),	& ! source
				1,		& ! tag
				Comm_Cart,	& ! communicator
				requests(2),	& ! request
				ierror)
	   call checkmpi(ierror,'(d)Bad +x async receive')
	else       
	   ! I'm sticking with one based for LE bits throughout,
	   ! or I'll get very confused
	   allocate(pinbuf(nd,nvecs,ny+2,nz+3))
           if(lelo.lt.0)then
	   call MPI_Irecv(	pinbuf(1,1,1,3-lelo), &
				57*(ny+2)*(nz+lelo+1),          &
				LBE_REAL,	& ! datatype
				leprocs(1,2),	& ! source
				tag_le_px3,	& ! tag
				Comm_Cart,	& ! communicator
				requests(1),	& ! request
				ierror)
	   call checkmpi(ierror,'Bad +x_le1 async receive')
	   call MPI_Irecv(	pinbuf(1,1,1,1), &
				57*(ny+2)*(-lelo+2),	   &
				LBE_REAL,	& ! datatype
				leprocs(1,1),	& ! source
				tag_le_px4,	& ! tag
				Comm_Cart,	& ! communicator
				requests(2),	& ! request
				ierror)

	   else !lelo >=0
         call MPI_Irecv(	pinbuf(1,1,1,nz-lelo+2), &
				(ny+2)*(lelo+2)*57,           &
				LBE_REAL,	& ! datatype
				leprocs(1,3),	& ! source
				tag_le_px1,	& ! tag
				Comm_Cart,	& ! communicator
				requests(1),	& ! request
				ierror)
	   call checkmpi(ierror,'(d)Bad +x_le1 async receive')
	   call MPI_Irecv(	pinbuf(1,1,1,1), &
				(ny+2)*(nz-lelo+1)*57,	   &
				LBE_REAL,	& ! datatype
				leprocs(1,2),	& ! source
				tag_le_px2,	& ! tag
				Comm_Cart,	& ! communicator
				requests(2),	& ! request
				ierror)
           endif !lelo >= or < 0
	   call checkmpi(ierror,'(d)Bad +x_le2 async receive')
	   xreq = xreq + 1
	endif

	! Receive from bottom into top halo
	if (leprocs(2,1) .eq. -1) then
	   allocate(minbuf(nd,nvecs,0:ny+1,0:nz+1))
	   call MPI_Irecv(	minbuf,		& ! buf
				(ny+2)*(nz+2)*57,	& ! count
				LBE_REAL,	& ! datatype
				nnprocs(1,2),	& ! source
				tag_mx,		& ! tag
				Comm_Cart,	& ! communicator
				requests(3),	& ! request
				ierror)
	   call MPI_Irecv(	xreq,		& ! buf
				1,	& ! count
				MPI_INTEGER,	& ! datatype
				nnprocs(1,2),	& ! source
				2,		& ! tag
				Comm_Cart,	& ! communicator
				requests(4),	& ! request
				ierror)
	   call checkmpi(ierror,'(d)Bad -x async receive')
	else
	   allocate(minbuf(nd,nvecs,ny+2,nz+3))
         if(lelo.lt.0)then
	   call MPI_Irecv(	minbuf(1,1,1,nz+lelo+2),  &
				57*(ny+2)*(2-lelo),     &
				LBE_REAL,	& ! datatype
				leprocs(2,3),	& ! source
				tag_le_mx3,	& ! tag
				Comm_Cart,	& ! communicator
				requests(3),	& ! request
				ierror)
	   call checkmpi(ierror,'Bad -x_le1 async receive')
	   call MPI_Irecv(	minbuf(1,1,1,1), &
				57*(ny+2)*(nz+1+lelo),    &
				LBE_REAL,	& ! datatype
				leprocs(2,2),	& ! source
				tag_le_mx4,	& ! tag
				Comm_Cart,	& ! communicator
				requests(4),	& ! request
				ierror)
         else !lelo >= 0
	   call MPI_Irecv(	minbuf(1,1,1,lelo+3),  &
				(ny+2)*(nz+1-lelo)*57,      &
				LBE_REAL,	& ! datatype
				leprocs(2,2),	& ! source
				tag_le_mx1,	& ! tag
				Comm_Cart,	& ! communicator
				requests(3),	& ! request
				ierror)
	   call checkmpi(ierror,'(d)Bad -x_le1 async receive')
	   call MPI_Irecv(	minbuf(1,1,1,1), &
				(ny+2)*(lelo+2)*57,     &
				LBE_REAL,	& ! datatype
				leprocs(2,1),	& ! source
				tag_le_mx2,	& ! tag
				Comm_Cart,	& ! communicator
				requests(4),	& ! request
				ierror)

	   call checkmpi(ierror,'(d)Bad -x_le2 async receive')
	   xreq = xreq + 1
         endif !lelo < or >= 0
 	endif


	!!!!!!
	! Begin asynchronous sends
	!!!!!!


	! Send from top buffer destined for bottom halo
	if (leprocs(2,1) .eq. -1) then
	   allocate(poutbuf(nd,nvecs,0:ny+1,0:nz+1))
!	   poutbuf = d_adv(:,:,nx,:,:)
	   do z=0,nz+1
	      do y=0,ny+1
		 poutbuf(:,:,y,z)=d_adv(:,:,nx,y,z)
	      end do
	   end do
	   call MPI_ISend(	poutbuf,	& ! buf
				(ny+2)*(nz+2)*57,	& ! count
				LBE_REAL,	& ! datatype
				nnprocs(1,2),	& ! dest
				tag_px,		& ! tag
				Comm_Cart,	& ! communicator
				requests(5),	& ! request
				ierror)
	   call MPI_ISend(	2,	& ! buf
				1,	& ! count
				MPI_INTEGER,	& ! datatype
				nnprocs(1,2),	& ! dest
				1,		& ! tag
				Comm_Cart,	& ! communicator
				requests(6),	& ! request
				ierror)
	   call checkmpi(ierror,'(d)Bad +x async send')
	else
	   allocate(poutbuf(nd,nvecs,ny+2,nz+3))
         if(lelo.lt.0)then
	   poutbuf(:,:,:,1:nz+lelo+1) = d_adv(:,:,nx,:,1:nz+lelo+1)
	   poutbuf(:,:,:,nz+2+lelo:nz+3) = d_adv(:,:,nx,:,nz+lelo-1:nz)
       	   call MPI_ISend(	poutbuf(1,1,1,1), &
				57*(ny+2)*(nz+lelo+1),	   &
				LBE_REAL,	& ! datatype
				leprocs(2,2),	& ! dest
				tag_le_px3,	& ! tag
				Comm_Cart,	& ! communicator
				requests(5),	& ! request
				ierror)
	   call checkmpi(ierror,'Bad +x async send')
	   call MPI_ISend(	poutbuf(1,1,1,nz+2+lelo), &
				57*(ny+2)*(-lelo+2),     &
				LBE_REAL,	& ! datatype
				leprocs(2,3),	& ! dest
				tag_le_px4,	& ! tag
				Comm_Cart,	& ! communicator
				requests(6),	& ! request
				ierror)           

         else  !lelo >= 0
	   ! array notation is terse:
	   poutbuf(:,:,:,1:lelo+2) = d_adv(:,:,nx,:,1:lelo+2)
	   poutbuf(:,:,:,lelo+3:nz+3) = d_adv(:,:,nx,:,lelo:nz)
	   call MPI_ISend(	poutbuf(1,1,1,1), &
				(ny+2)*(lelo+2)*57,  &
				LBE_REAL,	& ! datatype
				leprocs(2,1),	& ! dest
				tag_le_px1,	& ! tag
				Comm_Cart,	& ! communicator
				requests(5),	& ! request
				ierror)
	   call checkmpi(ierror,'(d)Bad +x async send')
	   call MPI_ISend(	poutbuf(1,1,1,lelo+3), &
				(ny+2)*(nz-lelo+1)*57,  &
				LBE_REAL,	& ! datatype
				leprocs(2,2),	& ! dest
				tag_le_px2,	& ! tag
				Comm_Cart,	& ! communicator
				requests(6),	& ! request
				ierror)
	   call checkmpi(ierror,'(d)Bad +x async send')
	   xreq = xreq + 1
         endif ! lelo>= or < 0
	   ! This works, so it is the MPI sending that has gone yucky
	   ! pinbuf(:,:,:,1:nz-lelo+1) = poutbuf(:,:,:,lelo+3:nz+3)
	endif

	! Send from bottom buffer destined for top halo
	if (leprocs(1,2) .eq. -1) then
	   allocate(moutbuf(nd,nvecs,0:ny+1,0:nz+1))
	   !moutbuf = d_adv(:,:,1,:,:)
	   do z=0,nz+1
	      do y=0,ny+1
		 moutbuf(:,:,y,z)=d_adv(:,:, 1,y,z)
	      end do
	   end do
	   call MPI_ISend(	moutbuf,	& ! buf
				(ny+2)*(nz+2)*57,	& ! count
				LBE_REAL,	& ! datatype
				nnprocs(1,1),	& ! dest
				tag_mx,		& ! tag
				Comm_Cart,	& ! communicator
				requests(7),	& ! request
				ierror)
	   call MPI_ISend(	2,	& ! buf
				1,	& ! count
				MPI_INTEGER,	& ! datatype
				nnprocs(1,1),	& ! dest
				2,		& ! tag
				Comm_Cart,	& ! communicator
				requests(8),	& ! request
				ierror)
	   call checkmpi(ierror,'(d)Bad -x async send')
	else
	   allocate(moutbuf(nd,nvecs,ny+2,nz+3))
         if(lelo.lt.0)then
	   moutbuf(:,:,:,1:2-lelo) = d_adv(:,:,1,:,1:2-lelo)
	   moutbuf(:,:,:,3-lelo:nz+3) = d_adv(:,:,1,:,-lelo:nz)
	   call MPI_ISend(	moutbuf(1,1,1,1), &
				57*(ny+2)*(2-lelo),     &
				LBE_REAL,	& ! datatype
				leprocs(1,1),	& ! dest
				tag_le_mx3,	& ! tag
				Comm_Cart,	& ! communicator
				requests(7),	& ! request
				ierror)
	   call checkmpi(ierror,'Bad +x async send')
	   call MPI_ISend(	moutbuf(1,1,1,3-lelo), &
				57*(ny+2)*(nz+1+lelo), &
				LBE_REAL,	& ! datatype
				leprocs(1,2),	& ! dest
				tag_le_mx4,	& ! tag
				Comm_Cart,	& ! communicator
				requests(8),	& ! request
				ierror)
         else !lelo >= 0
	   moutbuf(:,:,:,1:nz+1-lelo) = d_adv(:,:,1,:,1:nz+1-lelo)
	   moutbuf(:,:,:,nz+2-lelo:nz+3) = d_adv(:,:,1,:,nz-1-lelo:nz)
	   call MPI_ISend(	moutbuf(1,1,1,1), &
				(ny+2)*(nz+1-lelo)*57,     &
				LBE_REAL,	& ! datatype
				leprocs(1,2),	& ! dest
				tag_le_mx1,	& ! tag
				Comm_Cart,	& ! communicator
				requests(7),	& ! request
				ierror)
	   call checkmpi(ierror,'Bad +x async send')
	   call MPI_ISend(	moutbuf(1,1,1,nz+2-lelo), &
				(ny+2)*(lelo+2)*57, &
				LBE_REAL,	& ! datatype
				leprocs(1,3),	& ! dest
				tag_le_mx2,	& ! tag
				Comm_Cart,	& ! communicator
				requests(8),	& ! request
				ierror)
	   call checkmpi(ierror,'Bad +x async send')
	   xreq = xreq + 1
         endif !lwlo < or >= 0
	endif

	!!!!!!
	! Now wait for all I/O to complete
	!!!!!!

	call MPI_Waitall(8,requests,statuses,ierror)
!	print*,'X statuses:',statuses
	call checkmpi(ierror,'(d)MPI_Waitall() failed in X direction')


	!!!!!!
	! Interpolate points and write halos:
	!!!!!!

	! Write bottom halo
	if (leprocs(1,2) .ne. -1) then
	   ! FIXME A lot of this could be terser if array versions were done
	   do y=0,ny+1
	      do z=0,nz+1
		 d_adv(:,:,0,y,z) = (1-lefrac)*pinbuf(:,:,y+1,z+1) + lefrac*pinbuf(:,:,y+1,z+2)
	      end do
	   end do
	else
	   !d_adv(:,:,nx+1,:,:) = minbuf
	   !d_adv(:,:,0,:,:) = pinbuf
	   do z=0,nz+1
	      do y=0,ny+1
		 d_adv(:,:,0,y,z)=pinbuf(:,:,y,z)
	      end do
	   end do
	endif


	! Write top halo
	if (leprocs(2,1) .ne. -1) then
	   ! FIXME A lot of this could be terser if array versions were done
	   do y=0,ny+1
	      do z=0,nz+1
		 d_adv(:,:,nx+1,y,z) = lefrac*minbuf(:,:,y+1,z+1) + (1-lefrac)*minbuf(:,:,y+1,z+2)
	      end do
	   end do
	else
	   !d_adv(:,:,nx+1,:,:) = minbuf
	   !d_adv(:,:,0,:,:) = pinbuf
	   do z=0,nz+1
	      do y=0,ny+1
		 d_adv(:,:,nx+1,y,z)=minbuf(:,:,y,z)
	      end do
	   end do
	endif


	deallocate(pinbuf)
	deallocate(poutbuf)
	deallocate(minbuf)
	deallocate(moutbuf)

	!print*,'Rank ',myrankc,'completed X exchange.'

	!!!!
	! X swap done.
	!!!!
!#endif
!eyb
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Swap in the Y direction
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	allocate(pinbuf(3,19,0:nx+1,0:nz+1))
	allocate(poutbuf(3,19,0:nx+1,0:nz+1))
	allocate(moutbuf(3,19,0:nx+1,0:nz+1))
	allocate(minbuf(3,19,0:nx+1,0:nz+1))

	!poutbuf = d_adv(:,:,:,ny,:)
	!moutbuf = d_adv(:,:,:,1,:)



	!!!!!!
	! Begin asynchronous receives
	!!!!!!

	call MPI_Irecv(		pinbuf,		& ! buf
				(nx+2)*(nz+2)*57,	& ! count
				LBE_REAL,	& ! datatype
				nnprocs(2,1),	& ! source
				tag_py,		& ! tag
				Comm_Cart,	& ! communicator
				requests(pin),	& ! request
				ierror)
	call checkmpi(ierror,'(d)Bad +y async receive')

	call MPI_Irecv(		minbuf,		& ! buf
				(nx+2)*(nz+2)*57,	& ! count
				LBE_REAL,	& ! datatype
				nnprocs(2,2),	& ! source
				tag_my,		& ! tag
				Comm_Cart,	& ! communicator
				requests(mmin),	& ! request
				ierror)
	call checkmpi(ierror,'(d)Bad -y async receive')

	do z=0,nz+1
	 do x=0,nx+1
		poutbuf(:,:,x,z)=d_adv(:,:,x,ny,z)
		moutbuf(:,:,x,z)=d_adv(:,:,x, 1,z)
	 end do
	end do

	!!!!!!
	! Begin asynchronous sends
	!!!!!!

	call MPI_ISend(		poutbuf,	& ! buf
				(nx+2)*(nz+2)*57,	& ! count
				LBE_REAL,	& ! datatype
				nnprocs(2,2),	& ! dest
				tag_py,		& ! tag
				Comm_Cart,	& ! communicator
				requests(pout),	& ! request
				ierror)
	call checkmpi(ierror,'(d)Bad +y async send')

	call MPI_ISend(		moutbuf,	& ! buf
				(nx+2)*(nz+2)*57,	& ! count
				LBE_REAL,	& ! datatype
				nnprocs(2,1),	& ! dest
				tag_my,		& ! tag
				Comm_Cart,	& ! communicator
				requests(mout),	& ! request
				ierror)
	call checkmpi(ierror,'(d)Bad -y async send')

	!!!!!!
	! Now wait for all I/O to complete
	!!!!!!

	call MPI_Waitall(4,requests,statuses,ierror)
!	print*,'Y statuses:',statuses
	call checkmpi(ierror,'(d)MPI_Waitall() failed in Y direction')


	!d_adv(:,:,:,ny+1,:) = minbuf
	!d_adv(:,:,:,0,:) = pinbuf
	do z=0,nz+1
	 do x=0,nx+1
		d_adv(:,:,x,0,z)=pinbuf(:,:,x,z)
		d_adv(:,:,x,ny+1,z)=minbuf(:,:,x,z)
	 end do
	end do

	deallocate(pinbuf)
	deallocate(poutbuf)
	deallocate(minbuf)
	deallocate(moutbuf)

	!print*,'Rank ',myrankc,'completed Y exchange.'

	!!!!
	! Y swap done.
	!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Swap in the Z direction
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	allocate(pinbuf(3,19,0:nx+1,0:ny+1))
	allocate(poutbuf(3,19,0:nx+1,0:ny+1))
	allocate(moutbuf(3,19,0:nx+1,0:ny+1))
	allocate(minbuf(3,19,0:nx+1,0:ny+1))

	!poutbuf = d_adv(:,:,:,:,nz)
	!moutbuf = d_adv(:,:,:,:,1)


	!print*,'Rank ',myrankc,'Starting Z exchange.'

	!!!!!!
	! Begin asynchronous receives
	!!!!!!

	call MPI_Irecv(		pinbuf,		& ! buf
				(nx+2)*(ny+2)*57,	& ! count
				LBE_REAL,	& ! datatype
				nnprocs(3,1),	& ! source
				tag_pz,		& ! tag
				Comm_Cart,	& ! communicator
				requests(pin),	& ! request
				ierror)
	call checkmpi(ierror,'(d)Bad +z async receive')

	call MPI_Irecv(		minbuf,		& ! buf
				(nx+2)*(ny+2)*57,	& ! count
				LBE_REAL,	& ! datatype
				nnprocs(3,2),	& ! source
				tag_mz,		& ! tag
				Comm_Cart,	& ! communicator
				requests(mmin),	& ! request
				ierror)
	call checkmpi(ierror,'(d)Bad -z async receive')

	do y=0,ny+1
	 do x=0,nx+1
		poutbuf(:,:,x,y)=d_adv(:,:,x,y,nz)
		moutbuf(:,:,x,y)=d_adv(:,:,x,y, 1)
	 end do
	end do
	!print*,'send'

	!!!!!!
	! Begin asynchronous sends
	!!!!!!

	call MPI_ISend(		poutbuf,	& ! buf
				(nx+2)*(ny+2)*57,	& ! count
				LBE_REAL,	& ! datatype
				nnprocs(3,2),	& ! dest
				tag_pz,		& ! tag
				Comm_Cart,	& ! communicator
				requests(pout),	& ! request
				ierror)
	call checkmpi(ierror,'(d)Bad +z async send')

	call MPI_ISend(		moutbuf,	& ! buf
				(nx+2)*(ny+2)*57,	& ! count
				LBE_REAL,	& ! datatype
				nnprocs(3,1),	& ! dest
				tag_mz,		& ! tag
				Comm_Cart,	& ! communicator
				requests(mout),	& ! request
				ierror)
	call checkmpi(ierror,'(d)Bad -z async send')

	!!!!!!
	! Now wait for all I/O to complete
	!!!!!!

	call MPI_Waitall(4,requests,statuses,ierror)
	call checkmpi(ierror,'(d)MPI_Waitall() failed in Z direction')

	do y=0,ny+1
	 do x=0,nx+1
		d_adv(:,:,x,y,   0)=pinbuf(:,:,x,y)
		d_adv(:,:,x,y,nz+1)=minbuf(:,:,x,y)
	 end do
	end do

	deallocate(pinbuf)
	deallocate(poutbuf)
	deallocate(minbuf)
	deallocate(moutbuf)

	!print*,'Rank ',myrankc,'completed Z exchange.'

	!!!!
	! Z swap done.
	!!!!
end subroutine le_adv_dipole_exchange

!> Exchange both advected and real dipoles.
!>
!> Unlike \c lbe_adv_dipole_exchange(), this is optimised but simply
!> calls \c le_halo_exchange_*(), \c le_adapt_halo_x(), and \c
!> le_adv_dipole_exchange()
!>
!> \param d_adv just passed to \c lbe_adv_dipole_exchange()
!>
!> \param[in,out] lbe_N local chunk of the lattice with halo extent 1
!> (old LB3D style)
!>
!> \param[in] whole_N local lattice chunk with halo of extent \c
!> halo_extent
subroutine le_all_dipole_exchange(d_adv,lbe_N,whole_N)
	type(lbe_site),intent(inout) :: lbe_N(0:,0:,0:)
        type(lbe_site),intent(in)&
             & :: whole_N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
	real*8,dimension(1:,1:,0:,0:,0:) :: d_adv

        call le_halo_exchange_x(whole_N)
        call le_adapt_halo_x(lbe_N)
        call le_halo_exchange_yz(lbe_N)

	call le_adv_dipole_exchange(d_adv)
end subroutine le_all_dipole_exchange

#endif

!>Holds any code relating to Lees Edwards that needs running
!>at the end of the simulation.
subroutine le_cleanup

      implicit none

      if (ccoords(1) .eq. 0) then
         deallocate(topprocs)
      endif

      if ((ccoords(1)+1) .eq. cdims(1)) then
         deallocate(botprocs)
      endif

end subroutine le_cleanup

!> time-frequency of oscillatory shear
!>
!> \returns oscillation frequency
pure function le_shear_omega()
    real(kind=rk) :: le_shear_omega

    le_shear_omega = shear_omega
end function le_shear_omega

!> spatial offset between both Lees-Edwards planes
!>
!> \returns z-position offset
!>
!> \warning This is valid only after \c le_init_step() was called
!> until the end of the LB step.
pure function le_shear_offset()
    real(kind=rk) :: le_shear_offset

    le_shear_offset = 2.0_rk*shear_sum
end function le_shear_offset

!> velocity offset between both Lees-Edwards planes
!>
!> \returns z-velocity offset
!>
!> \warning This is valid only after \c le_init_step() was called
!> until the end of the LB step.
pure function le_shear_velocity()
    real(kind=rk) :: le_shear_velocity

    le_shear_velocity = 2.0_rk*shear_u
end function le_shear_velocity

end module lbe_leesedwards_module
