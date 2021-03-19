! This program converts an ASCII rock data file to XDR. 
! Input format has to be x y z stone, output format is x y z only
! Jens, 28.04.03

      program rockall2xdr
      implicit none
      
      character*80 OUTFILE,INFILE,prefix,ASCIIFILE
      integer I,IL,IB,num,x,y,z,stone
      integer ixdrs,ixdrr,ierr,ieer,type,ios
      
      PRINT *, 'File containing rock data ?'
      READ(*,*) INFILE

      IB = 0
      IL = INDEX(INFILE,' ')-1
      OUTFILE = INFILE(1:IL-4)//'.xdr'

        print *,'Output file: ',OUTFILE
	open(unit=20,file=INFILE(1:IL))
        call xdrfopen(ixdrs,OUTFILE,"w",ierr)

        do
           read(20,*,iostat=ios) x,y,z,stone
           if (ios .ne. 0) exit
	   ib = ib +1
	   if (stone.eq.1) then
	   call xdrfshort(ixdrs,x,ieer)
	   call xdrfshort(ixdrs,y,ieer)
	   call xdrfshort(ixdrs,z,ieer)
	   end if
        end do
        close(unit=20)
        call xdrfclose(ixdrs,ierr)
	PRINT *,'READ ',IB,' lines.'

      END

