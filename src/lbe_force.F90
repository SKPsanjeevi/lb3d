#include "lbe.h"

!> general code for a common forcing of both fluid and particles,
!> branches off to other modules that actually implement different
!> types of forcing
module lbe_force_module
    use lbe_force_constant_module, only: force_apply_constant&
         &,force_init_constant,force_input_constant
    use lbe_force_kolmogorov_module, only: force_apply_kolmogorov&
         &,force_init_kolmogorov,force_input_kolmogorov
    use lbe_globals_module
    use lbe_log_module
    use lbe_parallel_module, only: check_allocate, comm_cart, nnprocs
    use lbe_parms_module, only: force,nx,ny,nz
    use lbe_types_module, only: lbe_site

    implicit none
    include 'mpif.h'
    private

    public lbe_add_force_halo,lbe_force_apply,lbe_force_init,lbe_force_input&
         &,lbe_force_reset

    !> \name buffers to recv halo forces into
    !> \{
    !> indices: component, species (0: common force, 1-3: red, blue,
    !> green), x, y, z
    real*8,save,allocatable,dimension(:,:,:,:,:) :: x_rbuf,y_rbuf,z_rbuf
    !> \}

    !> mpi datatype representing the whole \c lbe_force data at a lattice
    !> point (in fortran something like \c real*8,dimension(3,0:n_spec) )
    integer,save :: forcedata_mpitype

    !> mpi datatypes representing \c [xyz]_rbuf .
    integer,save :: r_mpitype(3)

    !> \name more custom MPI data types
    !> \{
    !> MPI datatypes representing the lower (\c l)/upper (\c u) part
    !> of the force halo that is sent (\c s) during the exchange for
    !> each direction (indices 1/2/3 represent x/y/z). Datatype
    !> definitions are relative to the memory address of \c lbe_force
    !> as a whole.
    integer,save :: ls_mpitype(3) !< for lower halos
    integer,save :: us_mpitype(3) !< for upper halos
    !> \}

contains

    !> Builds a custom mpi datatype that represents the whole force
    !> data on a lattice point (in fortran: \c
    !> real*8,dimension(3,0:n_spec) )
    !>
    !> \param[out] fdmt mpi datatype created
    subroutine build_forcedata_mpitype(fdmt)
        integer,intent(out) :: fdmt ! type to be built
        real*8 :: sample(3,0:n_spec) ! sample of a force lattice point
        integer s_mt ! temporary mpi data type
        integer ierror
        integer(kind=MPI_ADDRESS_KIND) adr1,adr2,stride

        ! mpi datatype for all vector components of a species ( sample(:,s) )
        call mpi_get_address(sample(1,0),adr1,ierror)
        call mpi_get_address(sample(2,0),adr2,ierror)
        stride = adr2 - adr1
        call mpi_type_create_hvector(3,1,stride,LBE_REAL,s_mt,ierror)

        ! mpi datatype for whole point data ( sample(:,:) ), depending
        ! on n_spec we need (2 to 4) times s_mt :
        call mpi_get_address(sample(1,0),adr1,ierror)
        call mpi_get_address(sample(1,1),adr2,ierror)
        stride = adr2 - adr1
        call mpi_type_create_hvector(1+n_spec,1,stride,s_mt,fdmt,ierror)

        call mpi_type_commit(fdmt,ierror)
    end subroutine build_forcedata_mpitype

    !> Builds a custom mpi datatype that represents the part of \c f
    !> specified by the coordinate intervals \c xr(2),yr(2),zr(2) and
    !> stores it in \c fcmt . The starting indices of f are \c si(:).
    !>
    !> \param[in] si starting indices of \c f (differ depending on \c
    !> force_halo_extent )
    !>
    !> \param[in] f sample force lattice
    !>
    !> \param[in] xr chunk range in x-direction
    !>
    !> \param[in] yr chunk range in y-direction
    !>
    !> \param[in] zr chunk range in z-direction
    !>
    !> \param[out] fmct type to be built
    !>
    !> \warning The resulting MPI type will always relate to the
    !> spatial position \c (/0,0,0/) of \c f, thus if it is passed in
    !> an MPI call, the 1st-vector-component, 0th-species element at
    !> position \c (/0,0,0/) of \c f needs to be passed. This
    !> implementation assumes that this works even if position \c
    !> (/0,0,0/) does not exist as it is always the case for the
    !> lbe_force halo recv buffers.
    subroutine build_force_chunk_mpitype(si,f,xr,yr,zr,fcmt)
        integer,intent(in) :: si(3)
        real*8,intent(in) :: f(1:,0:,si(1):,si(2):,si(3):)
        integer,intent(in) :: xr(2),yr(2),zr(2)
        integer,intent(out) :: fcmt
        integer xrow_mt,xyplane_mt,xyzchunk_mt ! temporary mpi data types
        integer blocks,ierror
        integer(kind=MPI_ADDRESS_KIND) :: adr1,adr2,base,offset,stride
        integer(kind=MPI_ADDRESS_KIND) :: displs(1) 
        integer lengths(1)

        ! mpi datatype for slices of  f  like  f(:,:,xr(1):xr(2),y,z)
        blocks = 1 + xr(2) - xr(1)
        call mpi_get_address(f(1,0,1,1,1),adr1,ierror)
        call mpi_get_address(f(1,0,2,1,1),adr2,ierror)
        stride = adr2 - adr1
        call mpi_type_create_hvector(blocks,1,stride,forcedata_mpitype,xrow_mt,ierror)

        ! mpi datatype for slices of  f  like  f(:,:,xr(1):xr(2),yr(1):yr(2),z)
        blocks = 1 + yr(2) - yr(1)
        call mpi_get_address(f(1,0,1,1,1),adr1,ierror)
        call mpi_get_address(f(1,0,1,2,1),adr2,ierror)
        stride = adr2 - adr1
        call mpi_type_create_hvector(blocks,1,stride,xrow_mt,xyplane_mt,ierror)

        ! mpi datatype for whole chunk
        ! f(:,:,xr(1):xr(2),yr(1):yr(2),zr(1):zr(2))
        blocks = 1 + zr(2) - zr(1)
        call mpi_get_address(f(1,0,1,1,1),adr1,ierror)
        call mpi_get_address(f(1,0,1,1,2),adr2,ierror)
        stride = adr2 - adr1
        call mpi_type_create_hvector(blocks,1,stride,xyplane_mt,xyzchunk_mt,ierror)

        ! position of the beginning of the chunk relative to the beginning of  f
        call mpi_get_address(f(1,0,0,0,0),base,ierror)
        call mpi_get_address(f(1,0,xr(1),yr(1),zr(1)),offset,ierror)
        offset = offset - base

        ! fcmt becomes a datatype like  xyzchunk_mt  but relative to  base
        lengths = (/ 1 /)
        displs  = (/ offset /)
        call mpi_type_create_hindexed(1,lengths,displs,xyzchunk_mt,fcmt,ierror)
        call mpi_type_commit(fcmt,ierror)
    end subroutine build_force_chunk_mpitype

    !> initializes buffers and custom mpi types concerning \c lbe_force
    subroutine init_buffers_and_types
        integer stat
        integer :: h1,he

        h1=force_halo_extent-1
        he=force_halo_extent

        ! allocate and initialize arrays
        allocate (lbe_force(3,0:n_spec,-h1:nx+he,-h1:ny+he,-h1:nz+he)&
             &,x_rbuf(3,0:n_spec,     1:he,-h1:ny+he,-h1:nz+he)&
             &,y_rbuf(3,0:n_spec,-h1:nx+he,     1:he,-h1:nz+he)&
             &,z_rbuf(3,0:n_spec,-h1:nx+he,-h1:ny+he,     1:he),stat=stat)
        call check_allocate(stat&
             &,'init_buffers_and_types(): lbe_force,[xyz]_rbuf')
        lbe_force(:,:,:,:,:) = 0.0_8

        ! create custom mpi types for communication of force halo

        call build_forcedata_mpitype(forcedata_mpitype)

        ! types used for sending
        call build_force_chunk_mpitype((/-h1,-h1,-h1/),lbe_force&
             &,(/ -h1,    0/),(/ -h1,ny+he/),(/ -h1,nz+he/),ls_mpitype(1))
        call build_force_chunk_mpitype((/-h1,-h1,-h1/),lbe_force&
             &,(/nx+1,nx+he/),(/ -h1,ny+he/),(/ -h1,nz+he/),us_mpitype(1))
        call build_force_chunk_mpitype((/-h1,-h1,-h1/),lbe_force&
             &,(/ -h1,nx+he/),(/ -h1,    0/),(/ -h1,nz+he/),ls_mpitype(2))
        call build_force_chunk_mpitype((/-h1,-h1,-h1/),lbe_force&
             &,(/ -h1,nx+he/),(/ny+1,ny+he/),(/ -h1,nz+he/),us_mpitype(2))
        call build_force_chunk_mpitype((/-h1,-h1,-h1/),lbe_force&
             &,(/ -h1,nx+he/),(/ -h1,ny+he/),(/ -h1,    0/),ls_mpitype(3))
        call build_force_chunk_mpitype((/-h1,-h1,-h1/),lbe_force&
             &,(/ -h1,nx+he/),(/ -h1,ny+he/),(/nz+1,nz+he/),us_mpitype(3))

        ! types for receiving
        call build_force_chunk_mpitype((/  1,-h1,-h1/),x_rbuf&
             &,(/  1,   he/),(/-h1,ny+he/),(/-h1,nz+he/),r_mpitype(1))
        call build_force_chunk_mpitype((/-h1,  1,-h1/),y_rbuf&
             &,(/-h1,nx+he/),(/  1,   he/),(/-h1,nz+he/),r_mpitype(2))
        call build_force_chunk_mpitype((/-h1,-h1,  1/),z_rbuf&
             &,(/-h1,nx+he/),(/-h1,ny+he/),(/  1,   he/),r_mpitype(3))
    end subroutine init_buffers_and_types

    !> Add halo forces from neighbor processes to own \c lbe_force.
    !> Send own force halo to neighbors.
    !> Reset force halo to zero.
    !> Return directly if lbe_force is deactivated or if there is no force halo.
    subroutine lbe_add_force_halo
        integer ierror,status(MPI_STATUS_SIZE)
        integer :: h1,he

        if ( force_halo_extent == 0 ) return

        h1=force_halo_extent-1
        he=force_halo_extent

        ! send "downward" in x direction
        call mpi_sendrecv&
             &(lbe_force(1,0,0,0,0),1,ls_mpitype(1),nnprocs(1,1),0&
             &,x_rbuf(1,0,0,0,0),1,r_mpitype(1),nnprocs(1,2),0&
             &,comm_cart,status,ierror)
        lbe_force         (:,:,nx-h1:nx,-h1:ny+he,-h1:nz+he)&
             & = lbe_force(:,:,nx-h1:nx,-h1:ny+he,-h1:nz+he)&
             & +    x_rbuf(:,:,    1:he,-h1:ny+he,-h1:nz+he)

        ! send "upward" in x direction
        call mpi_sendrecv&
             &(lbe_force(1,0,0,0,0),1,us_mpitype(1),nnprocs(1,2),0&
             &,x_rbuf(1,0,0,0,0),1,r_mpitype(1),nnprocs(1,1),0&
             &,comm_cart,status,ierror)
        lbe_force         (:,:,    1:he,-h1:ny+he,-h1:nz+he)&
             & = lbe_force(:,:,    1:he,-h1:ny+he,-h1:nz+he)&
             & +    x_rbuf(:,:,    1:he,-h1:ny+he,-h1:nz+he)

        ! send "downward" in y direction
        call mpi_sendrecv&
             &(lbe_force(1,0,0,0,0),1,ls_mpitype(2),nnprocs(2,1),0&
             &,y_rbuf(1,0,0,0,0),1,r_mpitype(2),nnprocs(2,2),0&
             &,comm_cart,status,ierror)
        lbe_force         (:,:,-h1:nx+he,ny-h1:ny,-h1:nz+he)&
             & = lbe_force(:,:,-h1:nx+he,ny-h1:ny,-h1:nz+he)&
             & +    y_rbuf(:,:,-h1:nx+he,    1:he,-h1:nz+he)

        ! send "upward" in y direction
        call mpi_sendrecv&
             &(lbe_force(1,0,0,0,0),1,us_mpitype(2),nnprocs(2,2),0&
             &,y_rbuf(1,0,0,0,0),1,r_mpitype(2),nnprocs(2,1),0&
             &,comm_cart,status,ierror)
        lbe_force         (:,:,-h1:nx+he,    1:he,-h1:nz+he)&
             & = lbe_force(:,:,-h1:nx+he,    1:he,-h1:nz+he)&
             & +    y_rbuf(:,:,-h1:nx+he,    1:he,-h1:nz+he)

        ! send "downward" in z direction
        call mpi_sendrecv&
             &(lbe_force(1,0,0,0,0),1,ls_mpitype(3),nnprocs(3,1),0&
             &,z_rbuf(1,0,0,0,0),1,r_mpitype(3),nnprocs(3,2),0&
             &,comm_cart,status,ierror)
        lbe_force         (:,:,-h1:nx+he,-h1:ny+he,nz-h1:nz)&
             & = lbe_force(:,:,-h1:nx+he,-h1:ny+he,nz-h1:nz)&
             & +    z_rbuf(:,:,-h1:nx+he,-h1:ny+he,    1:he)

        ! send "upward" in z direction
        call mpi_sendrecv&
             &(lbe_force(1,0,0,0,0),1,us_mpitype(3),nnprocs(3,2),0&
             &,z_rbuf(1,0,0,0,0),1,r_mpitype(3),nnprocs(3,1),0&
             &,comm_cart,status,ierror)
        lbe_force         (:,:,-h1:nx+he,-h1:ny+he,    1:he)&
             & = lbe_force(:,:,-h1:nx+he,-h1:ny+he,    1:he)&
             & +    z_rbuf(:,:,-h1:nx+he,-h1:ny+he,    1:he)

        ! Clear halo.
        ! Comment by Timm (2012-02-28):
        ! I think that it is not necessary to clear the halo here.
        ! This could be done together with the reset of lbe_force in lbe_force_reset.
        lbe_force(:,:,     -h1:0, -h1:ny+he, -h1:nz+he) = 0.0_8
        lbe_force(:,:,nx+1:nx+he, -h1:ny+he, -h1:nz+he) = 0.0_8
        lbe_force(:,:, -h1:nx+he,     -h1:0, -h1:nz+he) = 0.0_8
        lbe_force(:,:, -h1:nx+he,ny+1:ny+he, -h1:nz+he) = 0.0_8
        lbe_force(:,:, -h1:nx+he, -h1:ny+he,     -h1:0) = 0.0_8
        lbe_force(:,:, -h1:nx+he, -h1:ny+he,nz+1:nz+he) = 0.0_8
    end subroutine lbe_add_force_halo

    !> apply force depending on the chosen implementation
    !>
    !> \param[in] lbe_N local lattice chunk with halo 1
    !>
    !> \param[in] whole_N local lattice chunk with full halo
    subroutine lbe_force_apply(lbe_N,whole_N)
        type(lbe_site),intent(in) :: lbe_N(0:,0:,0:)
        type(lbe_site),intent(in) :: &
             &whole_N(1-halo_extent:,1-halo_extent:,1-halo_extent:)

        if (.not.use_lbe_force) return

        select case (force)
        case ('constant')
           call force_apply_constant(lbe_N,whole_N)
        case ('kolmogorov')
           call force_apply_kolmogorov
        case ('none')
        end select
    end subroutine lbe_force_apply

    !> initialize common force data structures and branch off into
    !> actual force implementations
    subroutine lbe_force_init
        select case (force)
        case ('constant')
           call force_init_constant
        case ('kolmogorov')
           call force_init_kolmogorov
        case ('none')
        end select

        call init_buffers_and_types
    end subroutine lbe_force_init

    !> read in force namelists and initialize things that are required
    !> before other initialization routines are called
    subroutine lbe_force_input
        select case (force)
        case ('constant')
           call force_input_constant
        case ('kolmogorov')
           call force_input_kolmogorov
        case ('none')
        case default
           call error('unknown type of force: force="'//force//'"')
        end select
    end subroutine lbe_force_input

    !> Reset \c lbe_force to zero (as far as it is used).
    !>
    !> A possible force halo is reset already after the force halo
    !> summation.
    ! Comment by Timm (2012-02-28):
    ! It would be better to clear the total lbe_force here, instead of clearing the physical region here
    ! and the halo region only in lbe_add_force_halo.
    subroutine lbe_force_reset
      lbe_force(:, 0:n_spec, 1:nx, 1:ny, 1:nz) = 0.0_rk
    end subroutine lbe_force_reset

end module lbe_force_module
