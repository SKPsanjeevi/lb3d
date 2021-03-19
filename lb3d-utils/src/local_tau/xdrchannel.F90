!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!creates a pseudo 2-D channel with walls at x=1 and x=tnx
!also creates a file for local BC sets all rock to BB
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


program xdrchannel
implicit none

character(LEN=80) :: filename,bcfilename
integer :: x,y,z,tnx,tny,tnz,ierr,ixdrs,printer=10,localbc
real*8 :: ran, ran2


type	relax_list
		sequence	
		real*4	::	rock_state
end type relax_list

type	relax_list2
		sequence	
		real*8	::	bc_local
end type relax_list2



type(relax_list),dimension(:,:,:), allocatable	:: 	rocks
type(relax_list2),dimension(:,:,:), allocatable	:: 	rocks2

Print*,"Enter rockfile name without extension .xdr"
Read*, filename
Print*,"Size in x-dimension"
Read*,tnx
Print*,"Size in y-dimension"
Read*,tny
Print*,"Size in z-dimension"
Read*,tnz
!tny=4
!tnz=64


filename =trim(filename)//".xdr"
bcfilename = "localbc"//trim(filename)

allocate(rocks(0:tnx,0:tny,0:tnz))
allocate(rocks2(0:tnx,0:tny,0:tnz))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Put Rocks in the rocks structure
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do z=1,tnz
	do y=1,tny
	rocks(1,y,z)%rock_state = 1.0
	rocks(tnx,y,z)%rock_state = 1.0
	end do
end do
do x=2,tnx-1
	do z=1,tnz
		do y=1,tny
		rocks(x,y,z)%rock_state = 0.0
		end do
	end do
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!write into the xdr file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call xdrfopen(ixdrs,filename,"w",ierr)
do z=1,tnz
	do y=1,tny
		do x=1,tnx
		call xdrffloat(ixdrs,rocks(x,y,z)%rock_state,ierr)
		end do
	end do
end do
call xdrfclose(ixdrs,ierr)

call xdrfopen(ixdrs,bcfilename,"w",ierr)
do z=1,tnz
	do y=1,tny
		do x=1,tnx
		if(x==1)then
		rocks2(x,y,z)%bc_local = dble(rocks(x,y,z)%rock_state )+1.0
		endif
		if(x==tnx)then
		rocks2(x,y,z)%bc_local = dble(rocks(x,y,z)%rock_state )+1.0
		endif
		call xdrfdouble(ixdrs,rocks2(x,y,z)%bc_local,ierr)
		end do
	end do
end do
call xdrfclose(ixdrs,ierr)







end program xdrchannel
