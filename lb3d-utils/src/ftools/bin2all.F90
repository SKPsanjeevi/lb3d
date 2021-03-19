!
!=head1 C<bin2all.F90>
!
!This stand-alone program is part of the postprocessor.
!
!It reads in a Fortran unformatted data file, converts it to the ASCII 
!"all" format, and writes this to standard output.
!
!Takes eight command-line arguments:
!
!C<bin2ascii> I<nx> I<ny> I<nz> I<ox> I<oy> I<oz> I<type> I<name>
!
!=over 4
!
!=item *
!
!I<nx>, I<ny>, I<nz> : the size of the subdomain being processed.
!
!=item *
!
!I<ox>, I<oy>, I<oz> : The offsets to use
!
!=item *
!
!I<name> : the filename to read
!
!=item *
!
!I<type> : can be one of C<scalar>, C<2scalar>, C<3scalar>, C<vector>
!depending on the type of data.
!
!=back
!
!Some platforms (notably UNICOS/mk) do not support the C<getarg> call
!which returns the arguments passed on the command-line. In this case,
!the value C<PXFGETARG> should be defined when preprocessing this file,
!which will cause calls to C<getarg> to be replaced by a similar
!function.
!
!No arguments may be more than 64 chars long or buffers will overflow.
!
!=cut
!
! Yes, I know static buffers suck. Doing decent dynamic string handling
! in Fortran is like kicking dead whales down a beach.


program bin2all

	integer :: nx, ny, nz,ox,oy,oz
	integer :: x,y,z
	character(LEN=128) :: filename
	character(LEN=64) :: filetype
	character(LEN=64) :: buffer

	real, dimension(:,:,:), allocatable :: a,b,c
	real, dimension(:,:,:,:), allocatable :: vector

	logical	:: existp
	logical :: donep = .false.
	integer	:: ierror = 0

	if (iargc() .ne. 8) then
		print*,'Error: badly invoked.'
		print*,'Syntax: bin2ascii nx ny nz ox oy oz filetype filename'
		print*,'nx, ny, nz = subdomain dimensions'
		print*,'ox, oy, oz = subdomain offset'
		print*,'filename = name of unformatted file to read'
		print*,'filetype = scalar, 2scalar, 3scalar, vector, etc'
		stop
	endif

#ifdef PXFGETARG
	call pxfgetarg(1,buffer,ilen,ierror)
#else
	call getarg(1,buffer)
#endif
	read (buffer,*) nx
	if (nx <= 0 ) then
		print*,'Error: bad nx'
		stop
	endif

#ifdef PXFGETARG
	call pxfgetarg(2,buffer,ilen,ierror)
#else
	call getarg(2,buffer)
#endif
	read (buffer,*) ny
	if (ny <= 0 ) then
		print*,'Error: bad ny'
		stop
	endif

#ifdef PXFGETARG
	call pxfgetarg(3,buffer,ilen,ierror)
#else
	call getarg(3,buffer)
#endif
	read (buffer,*) nz
	if (nz <= 0 ) then
		print*,'Error: bad nz'
		stop
	endif

#ifdef PXFGETARG
	call pxfgetarg(4,buffer,ilen,ierror)
#else
	call getarg(4,buffer)
#endif
	read (buffer,*) ox
	if (ox < 0 ) then
		print*,'Error: bad ox'
		stop
	endif

#ifdef PXFGETARG
	call pxfgetarg(5,buffer,ilen,ierror)
#else
	call getarg(5,buffer)
#endif
	read (buffer,*) oy
	if (oy < 0 ) then
		print*,'Error: bad oy'
		stop
	endif

#ifdef PXFGETARG
	call pxfgetarg(6,buffer,ilen,ierror)
#else
	call getarg(6,buffer)
#endif
	read (buffer,*) oz
	if (oz < 0 ) then
		print*,'Error: bad oz'
		stop
	endif

#ifdef PXFGETARG
	call pxfgetarg(8,filename,ilen,ierror)
#else
	call getarg(8,filename)
#endif

#ifdef PXFGETARG
	call pxfgetarg(7,filetype,ilen,ierror)
#else
	call getarg(7,filetype)
#endif

	inquire(file=filename,exist=existp)
	if (.not. existp) then
		write(*,'(a,a,a)') 'Error: file <',filename,		&
					'> does not exist.'
		stop
	endif

	open(unit=10,file=filename,form="unformatted")

	if (filetype == 'scalar') then
		allocate(a(1:nx,1:ny,1:nz),stat=ierror)
		if (ierror .ne. 0) then
			print*,'Error: could not allocate array'
			stop
		endif

		read(unit=10,iostat=ierror) a(:,:,:)
		if (ierror .ne. 0) then
			print*,'Error ',ierror,' reading data'
			stop
		endif

		do z=1,nz
		 do y=1,ny
		  do x=1,nx
			write(*,'(3i3,f25.20)')		&
				x+ox,y+oy,z+oz,a(x,y,z)
		  end do
		 end do
		end do

		deallocate(a)

		donep = .true.
	end if
	if (filetype == '2scalar') then
		allocate(a(1:nx,1:ny,1:nz),stat=ierror)
		if (ierror .ne. 0) then
			print*,'Error: could not allocate array'
			stop
		endif
		allocate(b(1:nx,1:ny,1:nz),stat=ierror)
		if (ierror .ne. 0) then
			print*,'Error: could not allocate array'
			stop
		endif

		read(unit=10,iostat=ierror) a(:,:,:),b(:,:,:)
		if (ierror .ne. 0) then
			print*,'Error ',ierror,' reading data'
			stop
		endif

		do z=1,nz
		 do y=1,ny
		  do x=1,nx
			write(*,'(3i3,2f25.20)')		&
				x+ox,y+oy,z+oz,a(x,y,z),b(x,y,z)
		  end do
		 end do
		end do

		deallocate(a)
		deallocate(b)

		donep = .true.
	end if
	if (filetype == '3scalar') then
		allocate(a(1:nx,1:ny,1:nz),stat=ierror)
		if (ierror .ne. 0) then
			print*,'Error: could not allocate array'
			stop
		endif
		allocate(b(1:nx,1:ny,1:nz),stat=ierror)
		if (ierror .ne. 0) then
			print*,'Error: could not allocate array'
			stop
		endif
		allocate(c(1:nx,1:ny,1:nz),stat=ierror)
		if (ierror .ne. 0) then
			print*,'Error: could not allocate array'
			stop
		endif

		read(unit=10,iostat=ierror) a(:,:,:),b(:,:,:),c(:,:,:)
		if (ierror .ne. 0) then
			print*,'Error ',ierror,' reading data'
			stop
		endif

		do z=1,nz
		 do y=1,ny
		  do x=1,nx
			write(*,'(3i3,3f25.20)')		&
				x+ox,y+oy,z+oz,a(x,y,z),b(x,y,z),c(x,y,z)
		  end do
		 end do
		end do

		deallocate(a)
		deallocate(b)
		deallocate(c)

		donep = .true.
	end if

	if (filetype == 'vector') then
		allocate(vector(3,1:nx,1:ny,1:nz),stat=ierror)
		if (ierror .ne. 0) then
			print*,'Error: could not allocate array'
			stop
		endif

		read(unit=10,iostat=ierror) vector(:,:,:,:)
		if (ierror .ne. 0) then
			print*,'Error ',ierror,' reading data'
			stop
		endif

		do z=1,nz
		 do y=1,ny
		  do x=1,nx
			write(*,'(3i3,3f25.20)')		&
				x+ox,y+oy,z+oz,vector(:,x,y,z)
		  end do
		 end do
		end do

		deallocate(vector)

		donep = .true.
	end if

	close(unit=10)

	if (.not. donep) then
		print*,'Error: bad filetype'
		stop
	endif

end program bin2all
