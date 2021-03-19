! This program 'cuts' a 1D line out of a 3D data file. 
! Jens, 29.10.02

      program cut1d 
      implicit none
      
      character*80 OUTFILE,OUTFILE2,INFILE,prefix,ASCIIFILE
      integer A,I,J,K,IL,B,x,y,z
      integer NX,NY,NZ,ixdrs,ixdrr,ierr,ieer,type,ftype
      integer count
      real*8 scalar,scalar2
      real*4 scalartmp,scalartmp2,lobo,hibo 
      character*1 dir
      character*3 ACHR,BCHR 
! NGS
      real*8 min, max, avg
      real*8 out(3), sc
!.....out	vector out(1)==max
!		       out(2)==min
!		       out(3)==accum for avg
!
      max = -1.0e+30
      min = 1.0e+30
      out(3) = 0.0
      count = 0
!
      
      PRINT *, 'XDR FILE ?'
      READ(*,*) INFILE
      PRINT *,'FILE TYPE (1=double, 2=float) ?'
      READ(*,*) FTYPE
       IF ((FTYPE.gt.2).or.(FTYPE.lt.1)) STOP
      PRINT *, 'NX,NY,NZ ?'
      READ(*,*) NX,NY,NZ
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

      call int2str(A,ACHR,3)
      call int2str(B,BCHR,3)
      OUTFILE = INFILE(1:IL-4)//'_'//trim(ACHR)//'.'//trim(BCHR)//'.'//DIR
!
! NGS: Stats
      OUTFILE2 = INFILE(1:IL-4)//'_'//trim(ACHR)//'.'//trim(BCHR)//'.'//DIR//'.Stats'
      open(10,file=OUTFILE)

        print *,'Output file: ',OUTFILE
        call xdrfopen(ixdrr,INFILE(1:IL),"r",ierr)
	lobo = 0.20*NZ    !low boundary
        hibo = 0.80*NZ    !high boundary
          do z=1,NZ
           do y=1,NY
            do x=1,NX
               if (FTYPE.eq.1) then
                 call xdrfdouble(ixdrr, scalar, ieer)
                 IF (TYPE.gt.1) call xdrfdouble(ixdrr, scalar2, ieer)
               else
                 call xdrffloat(ixdrr, scalartmp, ierr)
                 scalar = real(scalartmp)
                 IF (TYPE.gt.1) THEN
                  call xdrffloat(ixdrr, scalartmp2, ierr)
                  scalar2 = real(scalartmp2)
                 ENDIF
               endif

              IF ((INDEX(DIR,'X').EQ.1).and.(y.eq.a).and.(z.eq.b)) THEN
                IF (TYPE.eq.1) THEN
                 write (10,*) x,scalar
! NGS: Stats
		 sc = abs(scalar)
		 if (sc .gt. max) max = sc
		 if (sc .lt. min) min = sc
		 out(3) = out(3) + sc
		 count = count + 1
!
                ELSE
                 write (10,*) x,scalar,scalar2
                ENDIF
              ENDIF
              IF ((INDEX(DIR,'Y').EQ.1).and.(x.eq.a).and.(z.eq.b)) THEN
                IF (TYPE.eq.1) THEN
                 write (10,*) y,scalar
! NGS: Stats
		 sc = abs(scalar)
		 if (sc .gt. max) max = sc
		 if (sc .lt. min) min = sc
		 out(3) = out(3) + sc
		 count = count + 1
!
                ELSE
                 write (10,*) y,scalar,scalar2
                ENDIF
              ENDIF
              IF ((INDEX(DIR,'Z').EQ.1).and.(x.eq.a).and.(y.eq.b)) THEN
                IF (TYPE.eq.1) THEN
                 write (10,*) z,scalar
! NGS: Stats
! not needed:	 sc = abs(scalar)
	         if((z.ge.lobo).and.(z.le.hibo)) then
		   if (scalar .gt. max) max = scalar
		   if (scalar .lt. min) min = scalar
		   out(3) = out(3) + scalar
		   count = count + 1
		 endif
!
                ELSE
                 write (10,*) z,scalar,scalar2
                ENDIF
              ENDIF

            end do
           end do
          end do
      call xdrfclose(ixdrr,ierr)
      close(10)
! ................timestep is INFILE(IL-9:IL-4)
      out(1) = min
      out(2) = max
      out(3) = out(3) / count
      open(11,file=OUTFILE2)
      write (11,*) INFILE(IL-9:IL-4),out(1),out(2),out(3)
      close(11)
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

