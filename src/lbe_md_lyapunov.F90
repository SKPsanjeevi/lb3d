#include "lbe.h"

!> Calculation of Lyapunov exponents for tracer particles---Aniruddha
!> knows this!
module lbe_md_lyapunov_module
#ifdef MD

    use lbe_globals_module, only: halo_extent,input_dfile_unit,rk,tsize,myrankc
    use lbe_init_rock_module, only: inside_rock_global
    use lbe_io_helper_module, only: lbe_make_filename_output
    use lbe_log_module
    use lbe_md_globals_module, only: communicate_rotations_s,md_input_file_unit&
         &,md_particle_type,ti_md_comm,ti_md_dump
    use lbe_md_helper_module, only: count_particles,log_msg_md,error_md&
         &,log_msg_md_hdr
    use lbe_md_module, only: borders,calculate_rotations_s,neighbor,setup_list
    use lbe_md_parallel_module, only: gather_particles,scatter_particles
    use lbe_timer_module, only: start_timer,stop_timer
    use lbe_parallel_module, only: check_allocate,comm_cart,gather_rock_state&
         &,tnx,tny,tnz
    use lbe_parms_module, only: chk_uid,inp_file,nt,arg_input_dfile_set&
         &,arg_input_dfile
    use lbe_types_module, only: lbe_site

    implicit none
    include 'mpif.h'
    private

    public input_lyapunov,lyapunov_run,setup_lyapunov,shutdown_lyapunov

    !> switches lyapunov code on/off
    logical,save,public :: lyapunov=.false.

    !> \{
    !> \name Resetting of pair trajectories
    !>
    !> If a pair distance is greater than \c d_max it is reset to \c d_aim.
    real(kind=rk),save :: d_aim=5.0_rk
    real(kind=rk),save :: d_max=10.0_rk
    !> \}

    !> perform lyapunov check every this many steps (0 means 'never')
    integer,save :: n_check=10
    !> perform lyapunov output every this many steps (0 means 'never')
    integer,save :: n_dump=1

    !> \{
    !> \name Buffers for different quantities
    !>
    !> All of them are (re-)allocated on demand by the root process.
    type(md_particle_type),allocatable,dimension(:),save :: tp !< particle cfg
    real(kind=rk),allocatable,dimension(:),save :: td !< pair distances
    real(kind=rk),allocatable,dimension(:),save :: trst !< traj. reset statuses
    integer,allocatable,dimension(:),save :: idx        !< uid indexes
    !> \}

    !> In this buffer, root holds a copy of the rock_state of the
    !> whole lattice plus a halo of extent \c halo_extent.
    !>
    !> It is needed for \c normalize() . The shape of the rock must
    !> not change after \c setup_lyapunov() was called!
    real(kind=rk),allocatable,dimension(:,:,:),save :: rock_state

    namelist /md_lyapunov/ d_aim,d_max,n_check,n_dump

contains

    !> read  \c /md_lyapunov/  from the md input file
    subroutine input_lyapunov
        integer ierror

        call log_msg_md_hdr("Reading MD Lyapunov input")

        if (myrankc.eq.0) then
           open (unit=md_input_file_unit,file=trim(inp_file)//'.md',err=100)
           read (unit=md_input_file_unit,nml=md_lyapunov,err=100)
           close (unit=md_input_file_unit,err=100)
           !call log_msg_md('read /md_lyapunov/ from file "'//trim(inp_file)&
           !     &//'.md"')
           !write (6,nml=md_lyapunov)
        end if

        if ( arg_input_dfile_set ) then
          call log_msg_md("  Getting differential input...")
          open(UNIT = input_dfile_unit, FILE = arg_input_dfile, STATUS = 'UNKNOWN')
          read(UNIT = input_dfile_unit, NML = md_lyapunov, IOSTAT = ierror)
          if (ierror .ne. 0) then
            call log_msg_md("    WARNING: Differential namelist not found or errors encountered.")
          endif
          close(UNIT = input_dfile_unit)
          call log_ws()
        end if

        write(msgstr,"('d_aim              = ',F16.10)") d_aim
        call log_msg(msgstr)
        write(msgstr,"('d_max              = ',F16.10)") d_max
        call log_msg(msgstr)
        write(msgstr,"('n_check            = ',I0)") n_check
        call log_msg(msgstr)
        write(msgstr,"('n_dump             = ',I0)") n_dump
        call log_msg(msgstr)
        call log_ws()

        call MPI_Bcast(d_aim,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(d_max,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(n_check,1,MPI_INTEGER,0,comm_cart,ierror)
        call MPI_Bcast(n_dump,1,MPI_INTEGER,0,comm_cart,ierror)

        return
100     continue
        call error_md('Error reading md input file "'//trim(inp_file)//'.md"')
    end subroutine input_lyapunov

    !> Creates an index in  index  that maps the unique id of each particle
    !> to the array index, at which the id is found in  uid . Both arguments
    !> need to be arrays of size  natoms  and must be allocated/deallocated
    !> by the caller.
    !>
    !> \warning This assumes that particle uids are contiguous and in
    !> the range \c (1:n_global) which the do not need to be.
    subroutine create_uid_index(n_global,uid,idx)
        integer,intent(in) :: n_global
        integer,intent(in) :: uid(:)
        integer,intent(out) :: idx(:)
        integer i

        do i=1,n_global
           idx(uid(i)) = i
        end do
    end subroutine create_uid_index

    !> This should be called every timestep, it branches to all other lyapunov
    !> routines needed during simulation after collecting all particles on
    !> root. At the end, the particles are scattered to all processors again,
    !> if their positions have been changed.
    subroutine lyapunov_run(N)
        type(lbe_site),intent(in) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        integer i,ii,j,ij,n_global,stat

        ! save the communication costs if there is nothing to do... (This
        ! could be written much easier if one could rely on shortcut
        ! evaluation as in C)
        if ((n_check==0.or.(n_check/=0.and.mod(nt,n_check)/=0)).and.&
             &(n_dump==0.or.(n_dump/=0.and.mod(nt,n_dump)/=0))) return

        if (myrankc==0) then
           call count_particles(n_global,0)

           ! allocate buffers only once on first invocation
           if (.not.allocated(tp)) then
              allocate (tp(n_global),td(n_global),trst(n_global),idx(n_global)&
                   &,stat=stat)
              call check_allocate(stat,'lyapunov_run(): tp,td,trst,idx')
           end if

           ! allow n_global to alter during the simulation
           if (n_global>size(tp)) then
              call log_msg_md('lyapunov buffers are too small,'&
                   &//' allocating more memory...')
              deallocate (tp,td,trst,idx)
              allocate (tp(2*n_global),td(2*n_global),trst(2*n_global)&
                   &,idx(2*n_global),stat=stat)
              call check_allocate(stat,'lyapunov_run(): tp,td,trst,idx')
           end if
        end if

        call start_timer(ti_md_comm)
        call gather_particles(tp)
        call stop_timer(ti_md_comm)

        if (myrankc==0) then
           call create_uid_index(n_global,tp(:)%uid,idx)

           ! calculate distances in  td(:)  - they are needed for checking
           ! as well as for dumping
           do i=1,n_global,+5
              ii = idx(i)
              do j=i+1,i+4
                 ij = idx(j)
                 td(ij) = distance(tp(ii)%x(:),tp(ij)%x(:))
                 ! this is overwritten by  check()  where appropriate:
                 trst(ij) = 0.0
              end do
              trst(ii) = 0.0       ! fixed trajectories are never reset
              td(ii) = 0.0 ! dumped distance is zero for fixed trajectories
           end do
        end if

        if (n_check/=0.and.mod(nt,n_check)==0) call check(n_global,N)
        if (n_dump/=0.and.mod(nt,n_dump)==0) then
           call start_timer(ti_md_dump)
           call dump_configuration_lyapunov(n_global)
           call stop_timer(ti_md_dump)
        end if
    end subroutine lyapunov_run

    !> root fills \c rock_state(:,:,:)  with data gathered from all processes
    subroutine setup_lyapunov(N)
        type(lbe_site),intent(in) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        integer stat

        if (myrankc==0) then
           allocate (rock_state&
                &(1-halo_extent:tnx+halo_extent&
                &,1-halo_extent:tny+halo_extent&
                &,1-halo_extent:tnz+halo_extent)&
                &,stat=stat)
           call check_allocate(stat,'setup_lyapunov(): rock_state')
        end if
        call gather_rock_state(N,rock_state)
    end subroutine setup_lyapunov

    !> this is called at the end of every simulation...
    subroutine shutdown_lyapunov
        if (myrankc==0) then
           deallocate (rock_state)
           if (allocated(tp)) deallocate(tp,td,trst,idx)
        end if
    end subroutine shutdown_lyapunov

    !> performs normalization if necessary
    subroutine check(n_global,N)
        integer,intent(in) :: n_global
        type(lbe_site),intent(in) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        integer i,j,ii,ij,ierror
        logical reset_something ! true if at least one particle has been reset

        ! root process does all work
        rank0: if (myrankc==0) then
           reset_something = .false.
           do i=1,n_global,+5
              ii = idx(i)
              do j=i+1,i+4
                 ij = idx(j)
                 if (td(ij)>=d_max) then
                    call normalize(N,ii,ij)
                    reset_something = .true.
                    ! reset status /= 0.0 indicates that normalization
                    ! occurred - abuse this to print the normalized distance...
                    trst(ij) = distance(tp(ii)%x(:),tp(ij)%x(:))
                 end if
              end do
           end do
        end if rank0

        call mpi_bcast(reset_something,1,MPI_LOGICAL,0,comm_cart,ierror)
        if (reset_something) then
           call setup_list
           call scatter_particles(n_global,tp)
           ! when at least one particle trajectory was reset then update
           ! neighbor lists (the equivalent of  exchange()  is already done
           ! in  scatter_particles() )
           if (communicate_rotations_s) call calculate_rotations_s
           call borders
           call neighbor
        end if
    end subroutine check

    !> calculates the distance between the two positions \c x1 and \c
    !> x2 without caring about pbc except in the z-direction
    real(kind=rk) function distance(x1,x2)
        real(kind=rk),intent(in) :: x1(3),x2(3)
        real(kind=rk) vec(3)

        vec(:) = x2(:)-x1(:)

        ! minimum image criterion applied to z component
        if (abs(vec(3))>tsize(3)/2.0_rk) then
           if (vec(3)<0.0_rk) then
              vec(3) = vec(3)+tsize(3)
           else
              vec(3) = vec(3)-tsize(3)
           end if
        end if

        distance = sqrt(dot_product(vec,vec))
    end function distance

    !> writes an ascii dump of  tx , tv , tuid , td , and  trst  to file
    subroutine dump_configuration_lyapunov(n_global)
        integer,intent(in) :: n_global
        integer i,ii
        character(len=1024) cfg_file_name
        integer,parameter :: cfg_file_unit=12

        rank0: if (myrankc==0) then
           call log_msg_md('dumping lyapunov cfg...')

           call lbe_make_filename_output(cfg_file_name,'lyapunov','.asc',nt)
           open (unit=cfg_file_unit,file=cfg_file_name,status='REPLACE'&
                &,action='WRITE',recl=160)

           do i=1,n_global
              ii = idx(i)
              write (unit=cfg_file_unit,fmt=&
                   &'(SP,3(ES15.8,X),3(ES15.8,X),SS,I10.10,ES15.8,ES15.8)'&
                   &) tp(ii)%x(:),tp(ii)%v(:),tp(ii)%uid,td(ii),trst(ii)
           end do
           close (cfg_file_unit)
        end if rank0
    end subroutine dump_configuration_lyapunov

    !> moves  tx(:,ij)  to some place on the line between  tx(:,ii)
    !> and  tx(:,ij) , for which the distance between  tx(:,ii)
    !> and  tx(:,ij)  is just greater than or equal to  d_aim  and
    !> which is not completely surrounded by rock sites. The subroutine
    !> reads the pair distance from  td(ij) .
    subroutine normalize(N,ii,ij)
        type(lbe_site),intent(in) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        integer,intent(in) :: ii,ij
        real(kind=rk) dir(3)

        ! direction from  tx(:,ii)  to  tx(:,ij)  normalized to 1
        dir(:) = (tp(ij)%x(:)-tp(ii)%x(:))/td(ij)

        ! this is the ideal position...
        tp(ij)%x(:) = tp(ii)%x(:) + dir(:)*d_aim

        ! ...but if it is inside rock, choose the next rock-free position
        ! in the direction of  tx(:,ij) :
        do while (inside_rock_global(rock_state,tp(ij)%x(:)))
           tp(ij)%x(:) = tp(ij)%x(:) + dir(:)*0.01
        end do
    end subroutine normalize

#endif
end module lbe_md_lyapunov_module
