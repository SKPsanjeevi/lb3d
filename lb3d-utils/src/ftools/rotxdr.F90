! This program has been written for debuugging purposes and does nothing
! else than reading a single precision scalar xdr file and writing it with x
! and z swapped.
!
! Written by Jens, 11.08.03


      program rotxdr 
      implicit none
      
      character(LEN=80) OUTFILE,INFILE,prefix,ASCIIFILE
      integer IA,I,J,K,L,IL,IB,num
      integer NX,NY,NZ,ixdrs,ixdrr,ierr,ieer,type
      real*4,allocatable :: scalar(:,:,:)
	real :: bla
      
      PRINT *, 'DOUBLE PRECISION XDR FILE ?'
      READ(*,*) INFILE
      PRINT *, 'NX,NY,NZ ?'
      READ(*,*) NX,NY,NZ

	allocate(scalar(nz,ny,nx))
      IL = INDEX(INFILE,' ')-1
      PREFIX = INFILE(1:IL-4)//'.2'
      OUTFILE = INFILE(1:IL-4)//'.2.xdr'

        print *,'Output file: ',OUTFILE
        call xdrfopen(ixdrr,INFILE(1:IL),"r",ierr)
        call xdrfopen(ixdrs,OUTFILE,"w",ierr)

          do i=1,nz
	   do j=1,ny
	    do k=1,nx
               call xdrffloat(ixdrr, scalar(i,j,k), ieer)
            end do
           end do
          end do

     call xdrfclose(ixdrr,ierr)
	print *,'Read infile',nx,ny,nz

	    do k=1,nx
	   do j=1,ny
          do i=1,nz
               call xdrffloat(ixdrs, scalar(i,j,k), ierr)
            end do
           end do
          end do
     call xdrfclose(ixdrs,ierr)
     call dump_avs_fld(trim(prefix),nz,ny,nx,1)

      END

!
!=head2 C<dump_avs_fld(prefix,nxi,nyi,nzi,1)>
!
! This routine writes an AVS field file
! Variable veclen determines fiedltype:
! 1=scalar,2=2scalar,3=3scalar,4=vector
!
! Added 21.06.02 by Jens
!
subroutine dump_avs_fld(prefix,nxi,nyi,nzi,vectmp)
        implicit none
        character(LEN=128) :: filename,fldname
        character(LEN=256) :: tmpstring
        character(LEN=*) :: prefix
        character(LEN=1)   :: vecstr,countstr,stridestr
        character(LEN=7)   :: skipstr
        integer :: nxi,nyi,nzi,veclen,i,i2,i3,vectmp,il

        ! Stupid compiler bug workaround on IRIX
        veclen = vectmp

        fldname=trim(prefix)//'.fld'
        open(10,file=fldname)
        write(10,'(a)') '# AVS field file'
        write(10,'(a)') 'ndim=3'
        write(10,'(a,i3)') 'dim1=',nxi
        write(10,'(a,i3)') 'dim2=',nyi
        write(10,'(a,i3)') 'dim3=',nzi
        write(10,'(a)') 'nspace=3'
        write(10,'(a)') 'field=uniform'

           if (veclen.eq.4) veclen=3
           filename=trim(prefix)//'.xdr'
           write(vecstr,'(i1.1)') veclen
           write(10,'(a)') 'veclen='//vecstr
           write(10,'(a)') 'data=xdr_float'
           do i=1,veclen
             write(countstr,'(i1.1)') i
             write(skipstr,'(i7.7)') 8*(i-1)
             tmpstring='variable '//countstr//' file='//trim(filename)// &
                ' filetype=binary skip='//skipstr//' stride='//vecstr
             write(10,'(a,a,a)') trim(tmpstring)
           end do
        close(10)
end subroutine dump_avs_fld

