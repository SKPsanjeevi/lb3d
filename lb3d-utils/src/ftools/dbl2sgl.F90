! This program converts a double precision XDR file to single precision.
! In addition, an AVS field file is generated.
! Jens, 09.10.02

      program dbl2sgl
      implicit none
      
      character*80 OUTFILE,INFILE,prefix,ASCIIFILE
      integer IA,I,J,L,IL,IB,num
      integer NX,NY,NZ,ixdrs,ixdrr,ierr,ieer,type
      real*8 scalar
      real*4 scalartmp 
      
      PRINT *, 'DOUBLE PRECISION XDR FILE ?'
      READ(*,*) INFILE
      PRINT *, 'NX,NY,NZ ?'
      READ(*,*) NX,NY,NZ
      PRINT *,'FIELD TYPE (1=scalar, 2=2scalar, 3=3scalar, 4=vector) ?'
      READ(*,*) TYPE

      IL = INDEX(INFILE,' ')-1
      PREFIX = INFILE(1:IL-4)//'.2'
      OUTFILE = INFILE(1:IL-4)//'.2.xdr'

        print *,'Output file: ',OUTFILE
        call xdrfopen(ixdrr,INFILE(1:IL),"r",ierr)
        call xdrfopen(ixdrs,OUTFILE,"w",ierr)

          do i=1,100000000
               call xdrfdouble(ixdrr, scalar, ieer)
               if (ieer.eq.0) goto 100
               scalartmp = real(scalar)
!               print *,scalartmp
               call xdrffloat(ixdrs, scalartmp, ierr)
          end do

100     call xdrfclose(ixdrs,ierr)
        call xdrfclose(ixdrr,ierr)
        call dump_avs_fld(trim(prefix),nx,ny,nz,type)

      END

!
!=head2 C<dump_avs_fld(prefix,nxi,nyi,nzi,veclen)>
!
! This routine writes an AVS field file
! Variable veclen determines fiedltype:
! 1=scalar,2=2scalar,3=3scalar,4=vector
!
! Added 21.06.02 by Jens
!
subroutine dump_avs_fld(prefix,nxi,nyi,nzi,vectmp)
        implicit none
        character*128 :: filename,fldname
        character*256 :: tmpstring
        character*(*) :: prefix
        character*1   :: vecstr,countstr,stridestr
        character*7   :: skipstr
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

