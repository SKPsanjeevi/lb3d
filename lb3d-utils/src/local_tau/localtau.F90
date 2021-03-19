!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!-reads in a xdr "rockfile" written in floates
!-callculates the distance from every lattice point to the nearest rock site
!-using this a local relaxation time according to the paper:
! "Capturing Knudsen layer phenomena using a lattice Boltzman model"
! by Zhang, Gu, Barber and Emerson
!
!the output filenames are taupos"rockfile"compi.xdr with i being the number of the component and their local !relaxation times
!and rockdist"rockfile".xdr is a file containing the distance of each lattice site to the next wall
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program localtau
implicit none

character(LEN=80) :: filename
character(LEN=80) :: filenamefull
character(LEN=80) :: filename3

character(LEN=80),dimension(:), allocatable :: filenametaupos
character(LEN=80),dimension(:), allocatable :: filenamedistance

integer :: x,y,z,tnx,tny,tnz
integer :: stonenumber,i
integer :: ierr,ixdrs,ixdrs2,ixdrs3,IO
integer :: components,simtype
real*4	:: stone
real*8, dimension(:), allocatable	:: kn

type	rock_list
		sequence	
		integer	::	x
		integer	::	y
		integer	::	z
end type rock_list
type	relax_list
		sequence	
		real*8	::	rock_state
		real*8	:: 	dist
		real*8	:: 	vartau1
		real*8	:: 	vartau2
		real*8	:: 	vartau3

end type relax_list

type physical_values
		sequence	
		real*8	::	pressure		!Pa
		real*8	::	molecmass		!kg
		real*8	::	dynvisc			!Pa*s
		real*8	::	lmean			!m
		real*8	::	dx			!m
		real*8	::	tau
end type physical_values

type phys_constants
		sequence	
		real*8	::	boltzmann_konst		!J/K---> J=Kg*m^2/s^2
		real*8	::	R			!J/(mol*K)
end type phys_constants

type(physical_values),dimension(:), allocatable	::	physvalues

type(rock_list),dimension(:), allocatable	:: 	rocklist
type(relax_list),dimension(:,:,:), allocatable	:: 	relax



real*8	::	temp_dist
real*8	::	maxl
real*8	::	temperature=292				!!K

real*8,parameter	::	boltzmann=1.3816504E-23			!J/K---> J=Kg*m^2/s^2
real*8,parameter	::	pi=3.141593
real*8,parameter	::	avo=6.02214179E23

!real*8	::	molecmass=1.00794*6.022E-26	!!Kg
!real*8	::	dynvisc=8.76E-6			!!Pas

real*8,	dimension(:), allocatable	::	lmean




!real*8	::	R=8.314472				!J/(mol*K)
real*8,	dimension(:), allocatable	::	tau
real*8	::	dx

Print*,"Enter rockfile name (without .xdr)"
Read*, filename
Print*,"Size in x-dimension"
Read*,tnx
Print*,"Size in y-dimension"
Read*,tny
Print*,"Size in z-dimension"
Read*,tnz


do
print*,"number of simulated components? (1,2 or 3)"
read*, components
if(components ==1 .or. components ==2 .or. components == 3)exit
end do

allocate(filenametaupos(1:components))
allocate(filenamedistance(1:components))



do 
print*," simulation of a channel with known Kn (press 1) or general simulation (press 2) "
read*,simtype
if(simtype ==1 .or. simtype ==2)exit
end do


if(simtype==1)then

allocate(kn(1:components))
allocate(lmean(1:components))
allocate(tau(1:components))
Print*,"Kn-number component 1 (channel flow)"
Read*,kn(1)

if(components ==2 .or. components == 3)then
Print*,"Kn-number component 2 (channel flow)"
Read*,kn(2)
end if

if(components == 3)then
Print*,"Kn-number component 3 (channel flow)"
Read*,kn(3)
end if

do i=1,components
lmean(i)=kn(i)*real(tnx-2)
print*,"mean free path is calculated as:", lmean(i), "simulation type 1first time"
tau(i)=kn(i)*sqrt((3.0*3.14159265358979323846)/8.0)*real(tnx-2)+0.5
end do

end if







if(simtype==2)then
allocate(physvalues(1:components))


Print*,"avg pressure in the system (Pa)"
Read*, physvalues(1)%pressure
Print*,"avg pressure in the system (Pa)", physvalues(1)%pressure

Print*,"molecular mass component1"
Read*,physvalues(1)%molecmass
Print*,"molecular mass component1",physvalues(1)%molecmass

Print*,"dynamic viscosity component1"
Read*,physvalues(1)%dynvisc
Print*,"dynamic viscosity component1",physvalues(1)%dynvisc

Print*,"dx in meters"
Read*,physvalues(1)%dx
Print*,"dx in meters",physvalues(1)%dx

Print*,"Relaxation time component 1 (red)"
Read*,physvalues(1)%tau
Print*,"Relaxation time component 1 (red)",physvalues(1)%tau


if(components ==2 .or. components == 3)then
physvalues(2)%pressure = physvalues(1)%pressure
Print*,"molecular mass component2"
Read*,physvalues(2)%molecmass
Print*,"dynamic viscosity component2"
Read*,physvalues(2)%dynvisc
physvalues(2)%dx=physvalues(1)%dx
Print*,"Relaxation time component 2 (blue)"
Read*,physvalues(2)%tau
end if

if(components == 3)then
physvalues(3)%pressure=physvalues(1)%pressure
Print*,"molecular mass component3"
Read*,physvalues(3)%molecmass
Print*,"dynamic viscosity component3"
Read*,physvalues(3)%dynvisc
physvalues(3)%dx=physvalues(1)%dx
Print*,"Relaxation time component 3 (surfactand)"
Read*,physvalues(3)%tau
end if

do i=1,components
physvalues(i)%lmean= ( physvalues(i)%dynvisc / physvalues(i)%pressure ) * &
     & ( (pi*boltzmann*temperature) / (2*physvalues(i)%molecmass) )**0.5 
print*,"lmean component:",i,"is",physvalues(i)%lmean,"in meters"
physvalues(i)%lmean = physvalues(i)%lmean/physvalues(i)%dx
print*,"lmean component:",i,"is",physvalues(i)%lmean
end do

end if





filenamefull = (trim(filename) // ".xdr")

allocate(relax(0:tnx,0:tny,0:tnz))
stonenumber=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read Rock-file, count rocksites
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call xdrfopen(ixdrs,filenamefull,"r",ierr)
do z=1,tnz
	do y=1,tny
		do x=1,tnx
			call xdrffloat(ixdrs,stone, ierr)
			relax(x,y,z)%rock_state = dble(stone)
			if(stone .ne. 0.)stonenumber=stonenumber+1
			
		end do
	end do
end do
call xdrfclose(ixdrs,ierr)




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Put rocksites in the rocklist
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate(rocklist(0:stonenumber))
i=0
do z=1,tnz
	do y=1,tny
		do x=1,tnx
			if(relax(x,y,z)%rock_state .ne. 0.) then
			i=i+1
			rocklist(i)%x=x
			rocklist(i)%y=y
			rocklist(i)%z=z
			end if
		end do
	end do
end do



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!calculate distance to next rock-site
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

relax(0:tnx,0:tny,0:tnz)%dist = ((tnx)**2+(tny)**2+(tnz)**2)**0.5
do z=1,tnz
	do y=1,tny
		do x=1,tnx
		if(relax(x,y,z)%rock_state == 0.) then
do i=1,stonenumber
temp_dist =((x-rocklist(i)%x)**2+(y-rocklist(i)%y)**2+(z-rocklist(i)%z)**2)**0.5

if(temp_dist .le. relax(x,y,z)%dist) then
relax(x,y,z)%dist = temp_dist
end if
!if(z==1)print*,temp_dist,x,y,z,relax(x,y,z)%dist
end do
		end if
		end do
	end do
end do


do i=1,stonenumber
!print*,rocklist(i)%x,rocklist(i)%y,rocklist(i)%z
end do




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! write distance and local relaxation-time files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




filenametaupos(1)=("taupos"//trim(filename) // "comp1.xdr")
if(components == 2 .or. components == 3)then
filenametaupos(2)=("taupos"//trim(filename) // "comp2.xdr")
if(components == 3)then
filenametaupos(3)=("taupos"//trim(filename) // "comp3.xdr")
end if 
end if


filename3 =("rockdist"//trim(filename) // ".xdr")



if(simtype == 1)then

call xdrfopen(ixdrs2,filenametaupos(1),"w",ierr)
do z=1,tnz
	do y=1,tny
		do x=1,tnx
		relax(x,y,z)%vartau1=(tau(1)-0.5)/(1.+0.7*(EXP(-1.*relax(x,y,z)%dist/lmean(1))))+0.5

		call xdrfdouble(ixdrs2,relax(x,y,z)%vartau1,ierr)
		!print*,x,y,z,relax(x,y,z)%vartau1
		end do
	end do
end do
call xdrfclose(ixdrs2,ierr)





if (components == 2 .or. components == 3)then


call xdrfopen(ixdrs2,filenametaupos(2),"w",ierr)
do z=1,tnz
	do y=1,tny
		do x=1,tnx
		relax(x,y,z)%vartau2=(tau(2)-0.5)/(1.+0.7*(EXP(-1.*relax(x,y,z)%dist/lmean(2))))+0.5

		call xdrfdouble(ixdrs2,relax(x,y,z)%vartau2,ierr)
		!print*,x,y,z,relax(x,y,z)%vartau2
		end do
	end do
end do
call xdrfclose(ixdrs2,ierr)






if (components == 3)then

call xdrfopen(ixdrs2,filenametaupos(3),"w",ierr)
do z=1,tnz
	do y=1,tny
		do x=1,tnx
		relax(x,y,z)%vartau3=(tau(3)-0.5)/(1.+0.7*(EXP(-1.*relax(x,y,z)%dist/lmean(3))))+0.5

		call xdrfdouble(ixdrs2,relax(x,y,z)%vartau3,ierr)
		!print*,x,y,z,relax(x,y,z)%vartau3
		end do
	end do
end do
call xdrfclose(ixdrs2,ierr)


end if
end if



else if(simtype == 2)then


call xdrfopen(ixdrs2,filenametaupos(1),"w",ierr)
do z=1,tnz
	do y=1,tny
		do x=1,tnx
		relax(x,y,z)%vartau1=(physvalues(1)%tau-0.5)/(1.+0.7*(EXP(-1.*relax(x,y,z)%dist/physvalues(1)%lmean)))+0.5

		call xdrfdouble(ixdrs2,relax(x,y,z)%vartau1,ierr)
		!print*,x,y,z,relax(x,y,z)%vartau1
		end do
	end do
end do
call xdrfclose(ixdrs2,ierr)

if (components == 2 .or. components == 3)then


call xdrfopen(ixdrs2,filenametaupos(2),"w",ierr)
do z=1,tnz
	do y=1,tny
		do x=1,tnx
		relax(x,y,z)%vartau2=(physvalues(2)%tau-0.5)/(1.+0.7*(EXP(-1.*relax(x,y,z)%dist/physvalues(2)%lmean)))+0.5

		call xdrfdouble(ixdrs2,relax(x,y,z)%vartau2,ierr)
		!print*,x,y,z,relax(x,y,z)%vartau2
		end do
	end do
end do
call xdrfclose(ixdrs2,ierr)


if (components == 3)then

call xdrfopen(ixdrs2,filenametaupos(3),"w",ierr)
do z=1,tnz
	do y=1,tny
		do x=1,tnx
		relax(x,y,z)%vartau3=(physvalues(3)%tau-0.5)/(1.+0.7*(EXP(-1.*relax(x,y,z)%dist/physvalues(3)%lmean)))+0.5

		call xdrfdouble(ixdrs2,relax(x,y,z)%vartau3,ierr)
		!print*,x,y,z,relax(x,y,z)%vartau3
		end do
	end do
end do
call xdrfclose(ixdrs2,ierr)


end if
end if


end if





call xdrfopen(ixdrs3,filename3,"w",ierr)
do z=1,tnz
	do y=1,tny
		do x=1,tnx
		call xdrfdouble(ixdrs3,relax(x,y,z)%dist,ierr)
		end do
	end do
end do
call xdrfclose(ixdrs3,ierr)

if(simtype==1)then

do i=1,components
print*,"mean free path for component",i,"is calculated as:", lmean(i),";simtype",simtype
end do
else 
do i=1,components
print*,"mean free path for component",i,"is calculated as:", physvalues(i)%lmean,";simtype",simtype
end do

end if
end program localtau















