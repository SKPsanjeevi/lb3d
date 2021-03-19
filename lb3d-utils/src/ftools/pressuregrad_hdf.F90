! This program 'cuts' a 1D line out of a 3D data file. 
! Jens, 29.10.02

      program cut1d 

USE HDF5

      implicit none

      INTEGER(HID_T) :: file_id,file_id2
      INTEGER(HID_T) :: dset_id,dset_id2
      INTEGER(HID_T) :: type_id,type_id2
      INTEGER(HID_T) :: space_id,space_id2
      INTEGER(HSIZE_T) :: dimensions(3),maxdims(3)
      INTEGER(HSIZE_T), DIMENSION(7) :: dims
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE:: dset_double,dset_double2
      REAL*4, DIMENSION(:,:,:), ALLOCATABLE:: dset_float,dset_float2
      REAL*8 AVG,AVGVEL,SUMXY,DP,SUMDP,PIN,POUT,DPINOUT,TMP,PERM,PERM2,PERM3,BETAPERM,BETAPERM2,BETAPERM3,CORR
      REAL*8, DIMENSION(:), ALLOCATABLE:: DPTMP,DPTMP2
      INTEGER NUMVOX

      REAL*8 n,m,b,r,sumz,sump,sumzp,sumz2,sump2
      character*80 OUTFILE,OUTFILE2,INFILE,INFILE2,prefix,ASCIIFILE
      integer A,I,J,K,IL,x,y,z,xmin,xmax,ymin,ymax,zmin,zmax
      integer NX,NY,NZ,ixdrs,ixdrr,ierr,ieer,type,ftype
      character*1 dir
      character*3 ACHR,BCHR 
      
      character(LEN=8), parameter :: dsetname='OutArray'
      integer :: error

      n=0.d0
      m=0.d0
      b=0.d0
      sumz=0.d0
      sumz2=0.d0
      sump=0.d0
      sump2=0.d0
      PRINT *, 'OD HDF FILE ?'
      READ(*,*) INFILE
      PRINT *, 'VEL HDF FILE ?'
      READ(*,*) INFILE2
      PRINT *,'FILE TYPE (1=double, 2=float) ?'
      READ(*,*) FTYPE
       IF ((FTYPE.gt.2).or.(FTYPE.lt.1)) STOP
      PRINT *,'XRANGE (X_MIN,X_MAX (1..NX))?'
      READ(*,*) XMIN,XMAX
      PRINT *,'YRANGE (Y_MIN,Y_MAX (1..NY))?'
      READ(*,*) YMIN,YMAX
      PRINT *,'ZRANGE (Z_MIN,Z_MAX (1..NZ))?'
      READ(*,*) ZMIN,ZMAX
      AVG=0.0
      AVGVEL=0.0
      NUMVOX=0

      IL = INDEX(INFILE,' ')-1

!      call int2str(A,ACHR,3)
!      call int2str(B,BCHR,3)
      write(ACHR,'(i3.3)') A
      write(BCHR,'(i3.3)') B
      OUTFILE = 'perm_'//INFILE(4:IL-3)//'.Dp'
      OUTFILE2 = INFILE(1:IL-3)//'.Dpz'
      open(10,file=OUTFILE)
      open(11,file=OUTFILE2)


      call h5open_f(error)
      call h5fopen_f(INFILE,H5F_ACC_RDONLY_F,file_id,error)
      call h5dopen_f(file_id,dsetname,dset_id,error)

      call h5fopen_f(INFILE2,H5F_ACC_RDONLY_F,file_id2,error)
      call h5dopen_f(file_id2,dsetname,dset_id2,error)

      call h5dget_space_f(dset_id,space_id,error)
      call h5sget_simple_extent_dims_f(space_id,dimensions,maxdims,error)
      call h5dget_space_f(dset_id2,space_id2,error)
      call h5sget_simple_extent_dims_f(space_id2,dimensions,maxdims,error)
      write(*,*)'Dims: ',dimensions

      nx=dimensions(1)
      ny=dimensions(2)
      nz=dimensions(3)
      if (nx.lt.xmax) STOP
      if (ny.lt.ymax) STOP
      if (nz.lt.zmax) STOP

      dims=(/nx,ny,nz,0,0,0,0/)
      
      if(FTYPE.eq.1)then
      allocate(dset_double(nx,ny,nz))
      allocate(dset_double2(nx,ny,nz))
       type_id=H5T_NATIVE_DOUBLE
       call h5dread_f(dset_id,type_id,dset_double,dims,error)
       call h5dread_f(dset_id2,type_id,dset_double2,dims,error)
      else
      allocate(dset_float(nx,ny,nz))
      allocate(dset_float2(nx,ny,nz))
       type_id=H5T_NATIVE_REAL
       call h5dread_f(dset_id,type_id,dset_float,dims,error)
       call h5dread_f(dset_id2,type_id,dset_float2,dims,error)
      endif

        print *,'Output file: ',OUTFILE
        print *,'Input file:  ',INFILE

          allocate(DPTMP(zmax))
          allocate(DPTMP2(zmax))
          do z=zmin,zmax
           SUMXY=0.d0
           DPTMP(z)=0.d0
           do y=ymin,ymax
            do x=xmin,xmax

              NUMVOX=NUMVOX+1
            if(FTYPE.eq.1)then
              AVG=AVG+dset_double(x,y,z)
              AVGVEL=AVGVEL+dset_double2(x,y,z)
            else
              AVG=AVG+dble(dset_float(x,y,z))
              AVGVEL=AVGVEL+dble(dset_float2(x,y,z))

              DPTMP(z)=DPTMP(z)+dble(dset_float(x,y,z))/3.d0
              !DPTMP(z)=DPTMP(z)+dble(dset_float(x,y,z))
            endif

            end do
           end do
             DPTMP(z)=DPTMP(z)/dble(((xmax-xmin+1)*(ymax-ymin+1)))
          end do
          AVGVEL=AVGVEL/dble(NUMVOX)
          AVG=AVG/dble(NUMVOX)
          ! avg DP
          DP=0.d0
          PIN=DPTMP(zmin)
          POUT=DPTMP(zmax)
          do z=zmin,zmax-1
           ! regression code:
           n = n + 1.0d0
           sumz  = sumz + z    ! compute sum of x
           sumz2 = sumz2 + z * z   ! compute sum of x**2
           sumzp = sumzp + z * DPTMP(z)   ! compute sum of x * y
           sump  = sump + DPTMP(z)        ! compute sum of y
           sump2 = sump2 + DPTMP(z)*DPTMP(z)

           ! back to pressure:
           DPTMP2(z)=DPTMP(z)-DPTMP(z+1)
           DP=DP+DPTMP2(z) 
          end do

          m = (n * sumzp  -  sumz * sump) / (n * sumz2 - sumz**2)            ! compute slope
          b = (sump * sump2  -  sumz * sumzp) / (n * sumz2  -  sumz**2)      ! compute y-intercept
          r = (sumzp - sumz * sump / n) /                              &     ! compute correlation coefficient
                     sqrt((sumz2 - sumz**2/n) * (sump2 - sump**2/n))

!   write (unit=*, fmt="(/a,es15.6)") " Slope        m = ", m                        ! print results
!   write (unit=*, fmt="(a, es15.6)") " y-intercept  b = ", b
!   write (unit=*, fmt="(a, es15.6)") " Correlation  r = ", r

          do z=zmin,zmax-1
           write (11,*) z,DPTMP(z),DPTMP2(z)
          end do

          DP=DP/dble(zmax-zmin)
          print *,XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX,NUMVOX
          write (10,*) '#XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX'
          write (10,*) XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX
          write (10,*) '#AVG density, velocity'
          write (10,*) AVG, AVGVEL
          write (10,*) '#pressure gradient: <DP>, (Pin-Pout)/L, lin reg <Pxy>'
          write (10,*) DP, (PIN-POUT)/dble(zmax-zmin+1), -m
          ! For tau=1, exact boundary position is slightly different => k(theocorr)
          CORR=0.1d0
          write (10,*) '#PERMEABILITY, k(<DP>), k(DPIN-DPOUT), k(linreg), k(theocorr)'
          PERM=AVGVEL/(DP)*(1.d0/6.d0*AVG)
          PERM2=AVGVEL/((PIN-POUT)/dble(zmax-zmin+1))*(1.d0/6.d0*AVG)
          PERM3=AVGVEL/(-m)*(1.d0/6.d0*AVG)
          write (10,*) PERM,PERM2,PERM3,dble((xmax-xmin+1-CORR)**2)/12.d0
          write (10,*) '#beta, (<DP>), (DPIN-DPOUT), (lin_reg)'
          betaperm=2*PERM/(xmax-xmin+1-CORR)**2-1.d0/6.d0
          betaperm2=2*PERM2/(xmax-xmin+1-CORR)**2-1.d0/6.d0
          betaperm3=2*PERM3/(xmax-xmin+1-CORR)**2-1.d0/6.d0
          write (10,*) betaperm,betaperm2,betaperm3

      close(10)
      close(11)

call h5dclose_f(dset_id,error)
call h5fclose_f(file_id,error)
call h5dclose_f(dset_id2,error)
call h5fclose_f(file_id2,error)
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

