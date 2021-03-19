#include "lbe.h"

!> routines that supports debugging
module lbe_md_debug_module
#ifdef MD
    use lbe_globals_module, only: halo_extent,lbe_force,tsize,myrankc
#ifdef USEXDRF
    use lbe_io_xdrf_module, only: dump_scalar_xdr
#endif
    use lbe_md_globals_module
    use lbe_parallel_module, only: check_allocate,comm_cart,find_topology&
         &,nprocs,start,tnx,tny,tnz
    use lbe_parms_module, only: nt,nx,ny,nz
    use lbe_types_module, only: lbe_site
    use lbe_md_helper_module, only: log_msg_md,error_md
    implicit none
    private
    public debug_check_send_bounds,debug_dump_lbe_force&
         &,debug_dump_rock_state_halo,debug_dump_total_rock_state_halo&
         &,debug_print_site,debug_print_nx_halo_rs,debug_dump_nx_halo_site&
         &,debug_print_nz_halo_rs&
#ifdef DEBUG_REPORTMDCOMM
         &,debug_reportmdcomm&
#endif
         &,debug_print_pz_halo_rs

    include 'mpif.h'
contains

    !> consistency check of position intervals in communication swaps
    subroutine debug_check_send_bounds()
        integer d,ierror,k,s

        do k=1,3
           do d=1,n_dir_max
              do s=1,sdirs(k,d)%scnt
                 if (sdirs(k,d)%s(s)%sbndhi<=sdirs(k,d)%s(s)%sbndlo) &
                      &then
                    write (unit=6,fmt='("myrankc=",I0,",dim=",I1,",dir=",I1,'&
                         &//'",send=",I0,",interval [",F6.2,":",F6.2,"[")') &
                         &myrankc,k,d,s&
                         &,sdirs(k,d)%s(s)%sbndlo,sdirs(k,d)%s(s)%sbndhi
                    call error_md&
                         &('debug_check_send_bounds(): empty send interval')
                 end if

                 if (s>1) then
                    if (mod(d,2)==1) then ! uneven direction means downward
                       if (sdirs(k,d)%s(s)%sbndlo/=sdirs(k,d)%s(s-1)%sbndhi&
!!$                            &.and.sdirs(k,d)%s(s)%sbndlo&
!!$                            &/=sdirs(k,d)%s(s-1)%sbndhi-tsize_pi(k)) &
                            &.and.sdirs(k,d)%s(s)%sbndlo&
                            &/=sdirs(k,d)%s(s-1)%sbndhi-tsize(k)) &
                            &then
                          write (unit=6,fmt='("myrankc=",I0,",dim=",I1,'&
                               &//'",dir=",I1,",sends=",I0,"/",I0,'&
                               &//'",intervals [",F6.2,":",F6.2,"[/[",F6.2,'&
                               &//'":",F6.2,"[")') &
                               &myrankc,k,d,s-1,s&
                               &,sdirs(k,d)%s(s-1)%sbndlo,sdirs(k,d)%s(s-1)%sbndhi&
                               &,sdirs(k,d)%s(s)%sbndlo,sdirs(k,d)%s(s)%sbndhi
                          call error_md('debug_check_send_bounds(): '&
                               &//'discontinuous send intervals')
                       end if
                    else if (mod(d,2)==0) then ! even direction means upward
                       if (sdirs(k,d)%s(s)%sbndhi/=sdirs(k,d)%s(s-1)%sbndlo&
!!$                            &.and.sdirs(k,d)%s(s)%sbndhi&
!!$                            &/=sdirs(k,d)%s(s-1)%sbndlo+tsize_pi(k)) &
                            &.and.sdirs(k,d)%s(s)%sbndhi&
                            &/=sdirs(k,d)%s(s-1)%sbndlo+tsize(k)) &
                            &then
                          write (unit=6,fmt='("myrankc=",I0,",dim=",I1,'&
                               &//'",dir=",I1,",sends=",I0,"/",I0,'&
                               &//'",intervals [",F6.2,":",F6.2,"[/[",F6.2,'&
                               &//'":",F6.2,"[")') &
                               &myrankc,k,d,s-1,s&
                               &,sdirs(k,d)%s(s-1)%sbndlo,sdirs(k,d)%s(s-1)%sbndhi&
                               &,sdirs(k,d)%s(s)%sbndlo,sdirs(k,d)%s(s)%sbndhi
                          call error_md('debug_check_send_bounds(): '&
                               &//'discontinuous send intervals')
                       end if
                    end if
                 end if
              end do
           end do
        end do
        call MPI_Barrier(MPI_COMM_WORLD,ierror)
        call log_msg_md('Successfully passed debug_check_send_bounds().')
    end subroutine debug_check_send_bounds

    !> Dump  sum(lbe_force(:,1,1:nx,1:ny,1:nz))  for everey lattice point to
    !> an xdr file.
    subroutine debug_dump_lbe_force
        real(kind=rk),allocatable :: sbuf(:),rbuf(:),lf(:,:,:)
        integer pcoords(3,0:nprocs-1) ! cartesian topology coordinates
        integer os(3)           ! lattice offset for different processors
        integer i,j,k,p,r,ierror,stat

        allocate (sbuf(nx*ny*nz),stat=stat)
        call check_allocate(stat,'debug_dump_lbe_force(): sbuf')
        if (myrankc==0) then
           allocate (rbuf(tnx*tny*tnz),stat=stat)
           call check_allocate(stat,'debug_dump_lbe_force(): rbuf')
        end if

        ! fill  sbuf
        p = 0
        do i=1,nx
           do j=1,ny
              do k=1,nz
                 p = p + 1
                 sbuf(p) = sum(lbe_force(:,1,i,j,k))
              end do
           end do
        end do

        call mpi_gather(sbuf,nx*ny*nz,MPI_REAL8&
             &,rbuf,nx*ny*nz,MPI_REAL8,0,comm_cart,ierror)
        deallocate (sbuf)

        rank0: if (myrankc==0) then
           ! copy data from  rbuf  to  lf
           allocate (lf(1:tnx,1:tny,1:tnz),stat=stat)
           call check_allocate(stat,'debug_dump_lbe_force(): lf')
           call find_topology(pcoords)
           p = 0
           do r=0,nprocs-1
              os(:) = pcoords(:,r)*(/nx,ny,nz/)
              do i=os(1)+1,os(1)+nx
                 do j=os(2)+1,os(2)+ny
                    do k=os(3)+1,os(3)+nz
                       p = p + 1
                       lf(i,j,k) = rbuf(p)
                    end do
                 end do
              end do
           end do
           deallocate (rbuf)
#ifdef USEXDRF
           call dump_scalar_xdr(lf,'sum-lbe-force')
#else
            call log_msg_md("XDRF support switched off")
#endif
           deallocate (lf)
        end if rank0
    end subroutine debug_dump_lbe_force

    !> dumps \c N(:,:,:)\%rock_state  of rank  \c i  to file  \p rs\<i\>.vtk .
    subroutine debug_dump_rock_state_halo(N,i)
        type(lbe_site),intent(in) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        integer, intent(in) :: i
        integer x,y,z
        character(3) :: counter
        character(80) :: f_name

        write (counter,'(I3.1)') i
        f_name = "rs-proc"//trim(counter)//".vtk"

        print *,"trying to put out vtk file",trim(f_name)

        open (unit=70,file=f_name,status='REPLACE',action='WRITE')
        write (unit=70,fmt='(A)') '# vtk DataFile Version 2.0'
        write (unit=70,fmt='(A)') 'debug_dump_rock_state_halo'
        write (unit=70,fmt='(A)') 'ASCII'
        write (unit=70,fmt='(A)') 'DATASET STRUCTURED_POINTS'
        write (unit=70,fmt='(A,3I4)') 'DIMENSIONS '&
             &,nx+2*halo_extent,ny+2*halo_extent,nz+2*halo_extent
        write (unit=70,fmt='(A,3I4)') 'ORIGIN '&
             &,1-halo_extent,1-halo_extent,1-halo_extent
        write (unit=70,fmt='(A)') 'SPACING 1 1 1'
        write (unit=70,fmt='(A,I12)') 'POINT_DATA '&
             &,product((/nx,ny,nz/)+2*halo_extent)
        write (unit=70,fmt='(A)') 'SCALARS OutArray double 1'
        write (unit=70,fmt='(A)') 'LOOKUP_TABLE default'
        do z=1-halo_extent,nz+halo_extent
           do y=1-halo_extent,ny+halo_extent
              do x=1-halo_extent,nx+halo_extent
                 write (unit=70,fmt='(ES15.8)') N(x,y,z)%rock_state
              end do
           end do
        end do
        close (unit=70)
    end subroutine debug_dump_rock_state_halo

    !> dumps  rock_state(:,:,:)  to file  'rs.vtk' .
    subroutine debug_dump_total_rock_state_halo(rock_state)
        real(kind=rk),intent(in) :: &
             &rock_state(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        integer x,y,z

        open (unit=70,file='rs.vtk',status='REPLACE',action='WRITE')
        write (unit=70,fmt='(A)') '# vtk DataFile Version 2.0'
        write (unit=70,fmt='(A)') 'debug_dump_total_rock_state_halo'
        write (unit=70,fmt='(A)') 'ASCII'
        write (unit=70,fmt='(A)') 'DATASET STRUCTURED_POINTS'
        write (unit=70,fmt='(A,3I4)') 'DIMENSIONS '&
             &,tnx+2*halo_extent,tny+2*halo_extent,tnz+2*halo_extent
        write (unit=70,fmt='(A,3I4)') 'ORIGIN '&
             &,1-halo_extent,1-halo_extent,1-halo_extent
        write (unit=70,fmt='(A)') 'SPACING 1 1 1'
        write (unit=70,fmt='(A,I12)') 'POINT_DATA '&
             &,product((/tnx,tny,tnz/)+2*halo_extent)
        write (unit=70,fmt='(A)') 'SCALARS OutArray double 1'
        write (unit=70,fmt='(A)') 'LOOKUP_TABLE default'
        do z=1-halo_extent,tnz+halo_extent
           do y=1-halo_extent,tny+halo_extent
              do x=1-halo_extent,tnx+halo_extent
                 write (unit=70,fmt='(ES15.8)') rock_state(x,y,z)
              end do
           end do
        end do
        close (unit=70)
    end subroutine debug_dump_total_rock_state_halo

    subroutine debug_print_nx_halo_rs(N,x,id_str,t_min)
        type(lbe_site),intent(in) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        integer,intent(in) :: x,t_min
        character(len=3),intent(in) :: id_str
        integer ierror,lx,p,y,z

        if (nt<t_min) return

        do p=0,nprocs-1
           call MPI_Barrier(MPI_COMM_WORLD,ierror)

           if (myrankc==p) then
              lx = x+1-start(1)
              if (lx==0) then
                 do y=1,ny
                    do z=1,nz
                       write (unit=6&
                            &,fmt='("DBG: ",A3,X,I6,X,I8,X,2(I4,X),F16.10)') &
                            &id_str,myrankc,nt,start(2)+y-1,start(3)+z-1&
                            &,N(lx,y,z)%rock_state
                    end do
                 end do
              end if
           end if

           call MPI_Barrier(MPI_COMM_WORLD,ierror)
        end do
    end subroutine debug_print_nx_halo_rs

    subroutine debug_dump_nx_halo_site(N,x,id_str,t_min)
        type(lbe_site),intent(in) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        integer,intent(in) :: x,t_min
        character(len=*),intent(in) :: id_str
        character(len=128) :: buf
        integer ierror,lx,p,y,z

        if (nt<t_min) return

        lx = x+1-start(1)
        if (lx==0) then
           write (unit=buf,fmt='(A,"_t",I8.8,"-p",I6.6,".asc")') trim(id_str),nt,myrankc
           open (unit=70,file=trim(buf),status='REPLACE',action='WRITE')
           do y=1,ny
              do z=1,nz
                 write (unit=70,fmt='(2(I4,X),20(F16.10,:,X))') &
                      &start(2:3)+(/y,z/)-1,N(lx,y,z)%rock_state,N(lx,y,z)%n_r
              end do
           end do
           close (unit=70)
        end if
    end subroutine debug_dump_nx_halo_site

    subroutine debug_print_nz_halo_rs(N,z,id_str,t_min)
        type(lbe_site),intent(in) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        integer,intent(in) :: z,t_min
        character(len=3),intent(in) :: id_str
        integer ierror,lz,p,x,y

        if (nt<t_min) return

        do p=0,nprocs-1
           call MPI_Barrier(MPI_COMM_WORLD,ierror)

           if (myrankc==p) then
              lz = z+1-start(3)
              if (lz==0) then
                 do x=1,nx
                    do y=1,ny
                       write (unit=6&
                            &,fmt='("DBG: ",A3,X,I6,X,I8,X,2(I4,X),F16.10)') &
                            &id_str,myrankc,nt,start(1)+x-1,start(2)+y-1&
                            &,N(x,y,lz)%rock_state
                    end do
                 end do
              end if
           end if

           call MPI_Barrier(MPI_COMM_WORLD,ierror)
        end do
    end subroutine debug_print_nz_halo_rs

    subroutine debug_print_pz_halo_rs(N,z,id_str,t_min)
        type(lbe_site),intent(in) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        integer,intent(in) :: z,t_min
        character(len=3),intent(in) :: id_str
        integer ierror,lz,p,x,y

        if (nt<t_min) return

        do p=0,nprocs-1
           call MPI_Barrier(MPI_COMM_WORLD,ierror)

           if (myrankc==p) then
              lz = z+1-start(3)
              if (lz==nz+1) then
                 do x=1,nx
                    do y=1,ny
                       write (unit=6&
                            &,fmt='("DBG: ",A3,X,I6,X,I8,X,2(I4,X),F16.10)') &
                            &id_str,myrankc,nt,start(1)+x-1,start(2)+y-1&
                            &,N(x,y,lz)%rock_state
                    end do
                 end do
              end if
           end if

           call MPI_Barrier(MPI_COMM_WORLD,ierror)
        end do
    end subroutine debug_print_pz_halo_rs

    subroutine debug_print_site(N,x,y,z,id_str)
        type(lbe_site),intent(in) :: &
             &N(1-halo_extent:,1-halo_extent:,1-halo_extent:)
        integer,intent(in) :: x,y,z
        character(len=3),intent(in) :: id_str
        integer :: lp(3)

        lp = (/x,y,z/)-(start-1)

        if (all(lp>=1).and.all(lp<=(/nx,ny,nz/))) then
           write (unit=6,fmt='("DBG: ",A3,X,I6,X,I8,X,3(I4,X),20(F16.10,:,X))') &
             &id_str,myrankc,nt,x,y,z&
             &,N(lp(1),lp(2),lp(3))%rock_state&
             &,N(lp(1),lp(2),lp(3))%n_r
        end if
    end subroutine debug_print_site

#ifdef DEBUG_REPORTMDCOMM
    subroutine debug_reportmdcomm()
        integer i,ierror,k,j,p

        do p=0,nprocs-1
           if (myrankc==p) then
              write (unit=6,fmt='(I6,": x=[",I3,":",I3,"], y=[",'&
                   &//'I3,":",I3,"], z=[",I3,":",I3,"]...")') &
                   &p,start(1),start(1)+nx-1,start(2),start(2)+ny-1&
                   &,start(3),start(3)+nz-1
              do k=1,3
                 write (unit=6,fmt='(" dim=",I1)') k
                 do i=1,n_dir_max
                    write (unit=6,advance='no',fmt=&
                         &'("  +(",I1,":",A12")->",I6,": ",I2,"r, ",I2,"s:")') &
                         &i,DIR_NAMES(i)&
                         &,sdirs(k,i)%sproc,sdirs(k,i)%rcnt,sdirs(k,i)%scnt
                    do j=1,sdirs(k,i)%scnt
                       write (unit=6,advance='no'&
                            &,fmt='(" [",F6.2,":",F6.2,"]")')&
                            &sdirs(k,i)%s(j)%sbndlo&
                            &,sdirs(k,i)%s(j)%sbndhi
                    end do
                    write (unit=6,fmt='()')
                 end do
              end do
           end if

           call MPI_Barrier(MPI_COMM_WORLD,ierror)
        end do
    end subroutine debug_reportmdcomm
#endif

#endif
end module lbe_md_debug_module
