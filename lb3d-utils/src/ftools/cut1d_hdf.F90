! This program 'cuts' a 1D line out of a 3D data file. 
! Jens, 29.10.02

      program cut1d 

USE HDF5

      implicit none

      INTEGER(HID_T) :: file_id      
      INTEGER(HID_T) :: dset_id
      INTEGER(HID_T) :: type_id
      INTEGER(HID_T) :: space_id
      INTEGER(HSIZE_T) :: dimensions(3),maxdims(3)
      INTEGER(HSIZE_T), DIMENSION(7) :: dims
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE:: dset_double
      REAL*4, DIMENSION(:,:,:), ALLOCATABLE:: dset_float

      character*80 OUTFILE,INFILE,prefix,ASCIIFILE
      integer A,I,J,K,IL,B,x,y,z
      integer NX,NY,NZ,ixdrs,ixdrr,ierr,ieer,type,ftype
      character*1 dir
      character*3 ACHR,BCHR 
      
      character(LEN=8), parameter :: dsetname='OutArray'
      integer :: error

      PRINT *, 'HDF FILE ?'
      READ(*,*) INFILE
      PRINT *,'FILE TYPE (1=double, 2=float) ?'
      READ(*,*) FTYPE
       IF ((FTYPE.gt.2).or.(FTYPE.lt.1)) STOP
      PRINT *,'FIELD TYPE (1=scalar, 2=2scalar) ?'
      READ(*,*) TYPE
       IF ((TYPE.gt.2).or.(TYPE.lt.1)) STOP
      PRINT *,'VARIABLE DIRECTION (X,Y,Z) ?'
      READ(*,*) DIR

      IF ((INDEX(DIR,'X').EQ.1).OR.(INDEX(DIR,'x').EQ.1)) THEN
        DIR='X' 
	PRINT *,'Y,Z ?'
      ELSE
       IF ((INDEX(DIR,'Y').EQ.1).OR.(INDEX(DIR,'y').EQ.1)) THEN
         DIR='Y' 
 	 PRINT *,'X,Z ?'
       ELSE	
        IF ((INDEX(DIR,'Z').EQ.1).OR.(INDEX(DIR,'z').EQ.1)) THEN
         DIR='Z' 
	 PRINT *,'X,Y ?'
        ELSE
         PRINT *,'UNKNOWN DIRECTION.'
         STOP
        ENDIF
       ENDIF
      ENDIF
      
      READ(*,*) A,B

      IL = INDEX(INFILE,' ')-1

!      call int2str(A,ACHR,3)
!      call int2str(B,BCHR,3)
      write(ACHR,'(i3.3)') A
      write(BCHR,'(i3.3)') B
      OUTFILE = INFILE(1:IL-3)//trim(ACHR)//'.'//trim(BCHR)//'.'//DIR
      open(10,file=OUTFILE)


      call h5open_f(error)
      call h5fopen_f(INFILE,H5F_ACC_RDONLY_F,file_id,error)
      call h5dopen_f(file_id,dsetname,dset_id,error)

      call h5dget_space_f(dset_id,space_id,error)
      call h5sget_simple_extent_dims_f(space_id,dimensions,maxdims,error)
      write(*,*)'Dims: ',dimensions

      nx=dimensions(1)
      ny=dimensions(2)
      nz=dimensions(3)

      dims=(/nx,ny,nz,0,0,0,0/)
      
      if(FTYPE.eq.1)then
      allocate(dset_double(nx,ny,nz))
       type_id=H5T_NATIVE_DOUBLE
       call h5dread_f(dset_id,type_id,dset_double,dims,error)
      else
      allocate(dset_float(nx,ny,nz))
       type_id=H5T_NATIVE_REAL
       call h5dread_f(dset_id,type_id,dset_float,dims,error)
      endif

        print *,'Output file: ',OUTFILE
        print *,'Input file:  ',INFILE

          do z=1,NZ
           do y=1,NY
            do x=1,NX

            if(FTYPE.eq.1)then
              write(*,*)'Doubles: ',dset_double(x,y,z),x,y,z

              IF ((INDEX(DIR,'X').EQ.1).and.(y.eq.a).and.(z.eq.b)) THEN
                 write (10,*) x,dset_double(x,y,z)
              ENDIF
              IF ((INDEX(DIR,'Y').EQ.1).and.(x.eq.a).and.(z.eq.b)) THEN
                 write (10,*) y,dset_double(x,y,z)
              ENDIF
              IF ((INDEX(DIR,'Z').EQ.1).and.(x.eq.a).and.(y.eq.b)) THEN
                 write (10,*) z,dset_double(x,y,z)
              ENDIF

            else
              write(*,*)'Floats: ',dset_float(x,y,z),x,y,z
              IF ((INDEX(DIR,'X').EQ.1).and.(y.eq.a).and.(z.eq.b)) THEN
                 write (10,*) x,dset_float(x,y,z)
              ENDIF
              IF ((INDEX(DIR,'Y').EQ.1).and.(x.eq.a).and.(z.eq.b)) THEN
                 write (10,*) y,dset_float(x,y,z)
              ENDIF
              IF ((INDEX(DIR,'Z').EQ.1).and.(x.eq.a).and.(y.eq.b)) THEN
                 write (10,*) z,dset_float(x,y,z)
              ENDIF
            endif

            end do
           end do
          end do

      close(10)

call h5dclose_f(dset_id,error)
call h5fclose_f(file_id,error)
CALL h5close_f(error)

      END

!
!===============================================
!
      SUBROUTINE INT2STR(IDIGIT,STRING,ISTR)
!
!.....Convert an integer number to a string
!
!.....INPUT / OUTPUT PARAMETERS
!.....IDIGIT      |  INTEGER TO BE CONVERTED         (IN)
!.....STRING      |  STRING RESULT OF CONVERSION     (OUT)
!.....ISTR        |  LENGTH OF STRING                (OUT)
!
      CHARACTER STRING*(*)
      INTEGER IDIGIT, ISTR, IASC(0:9)
      REAL*8 DIGIT
      IASC(0) = ICHAR('0')
      IASC(1) = ICHAR('1')
      IASC(2) = ICHAR('2')
      IASC(3) = ICHAR('3')
      IASC(4) = ICHAR('4')
      IASC(5) = ICHAR('5')
      IASC(6) = ICHAR('6')
      IASC(7) = ICHAR('7')
      IASC(8) = ICHAR('8')
      IASC(9) = ICHAR('9')

      INOW = 0
      IF ( IDIGIT .LT. 0 ) THEN
         INOW = INOW + 1
         STRING(1:1) = '-'
      ENDIF
      IF ( IDIGIT .EQ. 0 ) THEN
         INOW = INOW + 1
         STRING(INOW:INOW) = '0'
         NDIG = INOW
         GOTO 99
      ENDIF
      ITMP = IABS(IDIGIT)
      DIGIT = DBLE(ITMP)
      NDIG  = INT( LOG10(DIGIT) ) + 1
      ILEN = LEN(STRING)
      DO I=1, NDIG
         IAUX = ITMP / 10**(NDIG-I)
        INOW = INOW + 1
         STRING(INOW:INOW) = CHAR( IASC(IAUX) )
         ITMP = ITMP - IAUX*10**(NDIG-I)
      ENDDO
99    DO I=INOW+1, ILEN
         STRING(I:I) = ' '
      ENDDO
      ISTR = INOW
      RETURN
      END

