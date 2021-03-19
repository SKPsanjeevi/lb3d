#include "lbe.h"

!> implementation of massless and pointlike tracer particles
module lbe_md_fluid_tracer_module
#ifdef MD

    use lbe_globals_module, only: halo_extent, input_dfile_unit, myrankc
    use lbe_md_globals_module
    use lbe_md_helper_module, only: count_particles,fluid_velocity,log_msg_md&
         &,error_md,log_msg_md_hdr
    use lbe_log_module
    use lbe_parallel_module, only: calculate_displacements,check_allocate&
         &,comm_cart,nprocs,start
    use lbe_parms_module, only: inp_file,arg_input_dfile_set,arg_input_dfile
    use lbe_types_module, only: lbe_site

    implicit none
    include 'mpif.h'
    private

    public fluid_v_interaction_tracer,input_fluid_tracer,setup_fluid_tracer&
         &,summary_fluid_tracer

    !> \{
    !> \name Trace mask
    !>
    !> additional weight for each lbe particle species (will be mostly
    !> 1 or 0 to just switch tracing of a particular species on or
    !> off)
    real(kind=rk),save :: trace_mask_r=1.0_rk
    real(kind=rk),save :: trace_mask_b=1.0_rk
    real(kind=rk),save :: trace_mask_s=1.0_rk
    !> \}

    !> list of tracer particles lost in rock
    integer,dimension(:),allocatable,save :: stuck_particles

    namelist /md_fluid_tracer/ trace_mask_r,trace_mask_b,trace_mask_s

contains

    !> read  \c /md_fluid_tracer/
    subroutine input_fluid_tracer
        integer ierror

        call log_msg_md_hdr("Reading MD F Tracer input")

        if (myrankc.eq.0) then
           open (unit=md_input_file_unit,file=trim(inp_file)//'.md',err=100)
           read (unit=md_input_file_unit,nml=md_fluid_tracer,err=100)
           close (unit=md_input_file_unit,err=100)
           !call log_msg_md('read /md_fluid_tracer/ from file "'//trim(inp_file)&
           !     &//'.md"',.false.)
           !write (6,nml=md_fluid_tracer)
        end if

        if ( arg_input_dfile_set ) then
          call log_msg_md("  Getting differential input...")
          open(UNIT = input_dfile_unit, FILE = arg_input_dfile, STATUS = 'UNKNOWN')
          read(UNIT = input_dfile_unit, NML = md_fluid_tracer, IOSTAT = ierror)
          if (ierror .ne. 0) then
            call log_msg_md("    WARNING: Differential namelist not found or errors encountered.")
          endif
          close(UNIT = input_dfile_unit)
          call log_ws()
        end if

        write(msgstr,"('trace_mask_r       = ',F16.10)") trace_mask_r
        call log_msg(msgstr)
        write(msgstr,"('trace_mask_b       = ',F16.10)") trace_mask_b
        call log_msg(msgstr)
        write(msgstr,"('trace_mask_s       = ',F16.10)") trace_mask_s
        call log_msg(msgstr)
        call log_ws()

        call MPI_Bcast(trace_mask_r,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(trace_mask_b,1,MPI_REAL8,0,comm_cart,ierror)
        call MPI_Bcast(trace_mask_s,1,MPI_REAL8,0,comm_cart,ierror)

        return
100     continue
        call error_md('Error reading md input file "'//trim(inp_file)//'.md"')
    end subroutine input_fluid_tracer

    !> allocate arrays
    subroutine setup_fluid_tracer
        integer stat
        halo_exchange_before_md_run = .true.

        allocate (stuck_particles(0),stat=stat)
        call check_allocate(stat,'stuck_particles')
    end subroutine setup_fluid_tracer

    !> assign to every particle the velocity of the fluid field at its position,
    !> keep track of particles that have been trapped by rock
    subroutine fluid_v_interaction_tracer(N)
        type(lbe_site),intent(in) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        integer i,ii,j,stat
        logical stuck        ! true if a particle sees only rock sites
        logical stuck_before ! true if a particle was already stuck earlier
        integer,dimension(:),allocatable :: tmp ! tmp list of stuck particles

        i = atompnt
        particles: do ii = 1,nlocal
           ! set  P(i)%v(:)  just to the fluid velocity at point  P(i)%x(:)
           call fluid_velocity(N,(/trace_mask_r,trace_mask_b,trace_mask_s/)&
                &,P(i)%x,P(i)%vnew,stuck)

           ! clipping the interfacial velocity to a FIXED value 
           
           if(P(i)%vnew(1).gt.cutoff_v) then 
              P(i)%vnew(1) = cutoff_v
	      print *, 'cutoff_v just clipped a velocity'
           endif
           
           if(P(i)%vnew(1).lt.-cutoff_v) then
              P(i)%vnew(1) = -cutoff_v
	      print *, 'cutoff_v just clipped a velocity'
           endif
           
           
           if(P(i)%vnew(2).gt.cutoff_v) then
              P(i)%vnew(2) = cutoff_v
	      print *, 'cutoff_v just clipped a velocity'
           endif
           
           if(P(i)%vnew(2).lt.-cutoff_v) then
              P(i)%vnew(2) = -cutoff_v
	      print *, 'cutoff_v just clipped a velocity'
           endif
           
           if(P(i)%vnew(3).gt.cutoff_v) then
              P(i)%vnew(3) = cutoff_v
	      print *, 'cutoff_v just clipped a velocity'
           endif
           
           if(P(i)%vnew(3).lt.-cutoff_v) then
              P(i)%vnew(3) = -cutoff_v
	      print *, 'cutoff_v just clipped a velocity'
           endif

           ! make a list of stuck tracers
           if (stuck) then
              stuck_before = .false.
              do j=1,size(stuck_particles)
                 if (P(i)%uid==stuck_particles(j)) then
                    stuck_before = .true.
                    exit
                 end if
              end do
              if (.not.stuck_before) then
                 print '(A,I10,A,3ES10.2,A,3ES10.2,A)'&
                      &,'WARNING: lost a tracer in rock: uid='&
                      &,P(i)%uid,',x=(',P(i)%x(:),'),v=(',P(i)%v(:),')'
                 allocate (tmp(size(stuck_particles)),stat=stat)
                 call check_allocate(stat,'tmp')
                 tmp(:) = stuck_particles(:)
                 deallocate (stuck_particles)
                 allocate (stuck_particles(size(tmp)+1),stat=stat)
                 call check_allocate(stat,'stuck_particles')
                 stuck_particles(1:size(tmp)) = tmp(:)
                 deallocate (tmp)
                 stuck_particles(size(stuck_particles)) = P(i)%uid
              end if
           end if
           i = list(i)
        enddo particles
    end subroutine fluid_v_interaction_tracer

    !> writes conclusion concerning tracers to all units specified in
    !> \c units .
    subroutine summary_fluid_tracer(units)
        integer,intent(in) :: units(:)
        integer i,j,n,n_global,n_real,u,ierror,stat
        integer counts(0:nprocs-1) ! # of lost particles on each processor
        integer displs(0:nprocs-1) ! buffer displacements for  mpi_gatherv()
        integer,dimension(:),allocatable :: recvbuf

        !  myrankc==0  collects stuck particles from all processors and dumps
        ! their  uid .

        call mpi_gather(size(stuck_particles),1,MPI_INTEGER&
             &,counts,1,MPI_INTEGER,0,comm_cart,ierror)
        if (myrankc==0) n = sum(counts)
        call mpi_bcast(n,1,MPI_INTEGER,0,comm_cart,ierror)

        n_stuck: if (n==0) then          ! n==0, we are done
           if (myrankc==0) then
              do u=1,size(units)
                 write (unit=units(u),fmt='(A)') 'No single tracer was lost!'
                 write (unit=units(u),fmt='()')
              end do
           end if

        else n_stuck            ! n/=0, dump ids on root

           if (myrankc==0) then
              allocate (recvbuf(n),stat=stat)
              call check_allocate(stat,'recvbuf')
              call calculate_displacements(counts,displs)
           end if

           call mpi_gatherv(stuck_particles,size(stuck_particles),MPI_INTEGER&
                &,recvbuf,counts,displs,MPI_INTEGER,0,comm_cart,ierror)

           call count_particles(n_global,0)
           rank0: if (myrankc==0) then
              ! it is possible (yet very unlikely) that some tracers were
              ! registered as stuck on two processors, delete (set to -1) the
              ! duplicates...
              n_real = n
              do i=1,n
                 if (recvbuf(i)<0) cycle
                 do j=i+1,n
                    if (recvbuf(j)==recvbuf(i)) then
                       recvbuf(j) = -1
                       n_real = n_real-1
                    end if
                 end do
              end do

              do u=1,size(units)
                 write (unit=units(u),fmt='(A,I10,A,I10,A,F12.3,A)') &
                      &'Of ',n_global,' tracers, ',n_real,' were lost ('&
                      &,(n_real*100.0)/n_global,'%). Dumping ids...'
              end do

              do i=1,n
                 do u=1,size(units)
                    if (recvbuf(i)>0) &
                         &write (unit=units(u),fmt='(I10)') recvbuf(i)
                 end do
              end do

              do u=1,size(units)
                 write (unit=units(u),fmt='()')
              end do

              deallocate (recvbuf)
           end if rank0
        end if n_stuck
    end subroutine summary_fluid_tracer

#endif
end module lbe_md_fluid_tracer_module
