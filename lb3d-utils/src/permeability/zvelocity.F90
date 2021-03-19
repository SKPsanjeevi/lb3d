! z averaging and permeability calculation 
! Frank, 11/2007
!based on :

! Attention !!!!!!!
! a) This program assumes that fr was set to 0.7 in lb3d's input file.
! b) Although performing the correct physical calculations, they 
!    are carried out in a numerical weak fashion.
! c) There exists a c-code permcalc.c which is ding the same calculations.
!    Both programms should return apprx. the same results.    
   

      program zvelocity

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
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE:: pdset_double
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE:: vdset_double
      REAL*4, DIMENSION(:,:,:), ALLOCATABLE:: pdset_float 
      REAL*4, DIMENSION(:,:,:), ALLOCATABLE:: vdset_float 
      REAL*8, DIMENSION(:), ALLOCATABLE:: zvel_double 
      REAL*8, DIMENSION(:), ALLOCATABLE:: density_double
      REAL*8, DIMENSION(:), ALLOCATABLE:: pressure_double
      REAL*8, DIMENSION(:), ALLOCATABLE:: massflow_double
      REAL*8, DIMENSION(:), ALLOCATABLE:: massflow_tot
      REAL*8, DIMENSION(:), ALLOCATABLE:: filling_double
      REAL*8 ::vvalue,pvalue, vsample, porsample,perm, tau, lattice,dx , c , factor
      REAL*8 ::psample,chanx,chany,c1,c2,vp,rhop,dynvis,deltap_p,rhosample_p,perm_L,vz_P
      REAL*8 ::pmin,pmax,perm_mm,perm_mu,mflow_mean,rcount, delta_rho_rel, perm_au
      REAL*8 ::perm_tot,vsample_tot,vz_P_tot
      REAL*8 ::vmin,vmax
      REAL*8, DIMENSION(:), ALLOCATABLE:: zvel_double_tot
      character*80 OUTFILE,VINFILE,PINFILE,prefix,ASCIIFILE
      integer A,I,J,K,IL,B,x,y,z
      integer NX,NY,NZ,ixdrs,ixdrr,ierr,ieer,type,ftype
      character*1 dir
      character*3 ACHR,BCHR 
      
      character(LEN=8), parameter :: dsetname='OutArray'
      integer :: error

      ! read in this section with mflowscript.sh
      PRINT *, 'HDF vel FILE ?' ! z_velocity file
      READ(*,*) VINFILE
      PRINT *, 'HDF od FILE ?'  ! oil density file
      READ(*,*) PINFILE
      PRINT *,'FILE TYPE (1=double, 2=float) ?'
      READ(*,*) FTYPE
       IF ((FTYPE.gt.2).or.(FTYPE.lt.1)) STOP
      PRINT *,'FIELD TYPE (1=scalar, 2=2scalar) ?'
      READ(*,*) TYPE
      IF ((TYPE.gt.2).or.(TYPE.lt.1)) STOP
      PRINT *,'channel width in x,y, direction (int,int) ?' !channel half width, i.e. 
      READ(*,*) chanx,chany   ! porous medium extends from 1+chanx ... nx-chanx
      PRINT *, 'Z-range of sample  ?' 
      READ(*,*) a,b
      write(ACHR,'(i3.3)') A
      write(BCHR,'(i3.3)') B
      print *,  'tau (~viscosity) of liquid (double)?' 
      read(*,*) tau
      print *, 'lattice constant of rock in micron? (double)' 
      read(*,*) lattice
      OUTFILE='flow.sum'
      open(10,file=OUTFILE)

      
      ! open vel file and read data to vdset_float(double)
      call h5open_f(error)
      call h5fopen_f(VINFILE,H5F_ACC_RDONLY_F,file_id,error)
      call h5dopen_f(file_id,dsetname,dset_id,error)

      call h5dget_space_f(dset_id,space_id,error)
      call h5sget_simple_extent_dims_f(space_id,dimensions,maxdims,error)
      write(*,*)'Dims: ',dimensions

      nx=dimensions(1)
      ny=dimensions(2)
      nz=dimensions(3)

      dims=(/nx,ny,nz,0,0,0,0/)
      
      if(FTYPE.eq.1)then
      allocate(vdset_double(nx,ny,nz))
       type_id=H5T_NATIVE_DOUBLE
       call h5dread_f(dset_id,type_id,vdset_double,dims,error)
       !DEBUG write(*,*) 'some values from double  ', vdset_double(3,3,:)
      else
      allocate(vdset_float(nx,ny,nz))
       type_id=H5T_NATIVE_REAL
       call h5dread_f(dset_id,type_id,vdset_float,dims,error)
       !DEBUG write(*,*) 'some values from float ', vdset_float(2,2,1), vdset_float(2,2,2), vdset_float(2,2,3)
      endif

        print *,'Output file: ',OUTFILE
        print *,'V Input file:  ',VINFILE


      close(10)

      call h5dclose_f(dset_id,error)
      call h5fclose_f(file_id,error)
      CALL h5close_f(error)
      
      ! convert lattice from mu -> m 
      dx=lattice	
      lattice=lattice/1.0E6
      ! speed of sound of water in SI:
      c1=1470*sqrt(3.)
      ! other constant
      c2=1/(lattice*lattice*lattice)

      ! memory allocation for z-cut values:
      allocate(zvel_double(nz))
      allocate(zvel_double_tot(nz))
      allocate(density_double(nz))
      allocate(pressure_double(nz))
      allocate(massflow_double(nz))      
      allocate(massflow_tot(nz))      
      allocate(filling_double(nz))

      
      open(10,file=OUTFILE)

      ! read in density data:
      call h5open_f(error)
      call h5fopen_f(PINFILE,H5F_ACC_RDONLY_F,file_id,error)
      call h5dopen_f(file_id,dsetname,dset_id,error)
      call h5dget_space_f(dset_id,space_id,error)
      call h5sget_simple_extent_dims_f(space_id,dimensions,maxdims,error)
      write(*,*)'Dims: ',dimensions
      nx=dimensions(1)
      ny=dimensions(2)
      nz=dimensions(3)
      dims=(/nx,ny,nz,0,0,0,0/)
      if(FTYPE.eq.1)then
         allocate(pdset_double(nx,ny,nz))
         type_id=H5T_NATIVE_DOUBLE
         call h5dread_f(dset_id,type_id,pdset_double,dims,error)
         !DEBUG       write(*,*) 'some values from double  ', pdset_double(3,3,:)
      else
         allocate(pdset_float(nx,ny,nz))
         type_id=H5T_NATIVE_REAL
         call h5dread_f(dset_id,type_id,pdset_float,dims,error)
      endif
      call h5dclose_f(dset_id,error)
      call h5fclose_f(file_id,error)
      call  h5close_f(error)
      
      print *,'Output file: ',OUTFILE
      print *,'OD Input file:  ',PINFILE
      !DEBUG    write(10, *) 'Test'
      !DEBUG    write(*,*) 'zvel at z=87: ', zvel_double(87)
      
      ! initialize average values:
       porsample=0.0
       psample=0.0
       vsample=0.0
       vsample_tot=0.0
       mflow_mean=0.0
       pmin=1.
       pmax=0.
       vmin=1e8
       vmax=-1e8
       rcount=0.    
 
       ! average over system in, creating z-cuts:
       do z=1,nz
          filling_double(z)=0.0
          zvel_double(z)=0.0
          zvel_double_tot(z)=0.0
          density_double(z)=0.0
          pressure_double(z)=0.0
          massflow_double(z)=0.0
          massflow_tot(z)=0.0
          do y=1+chany,ny-chany
             do x=1+chanx,nx-chanx
                !DEBUG             write(10, *) vvalue
                IF(TYPE.eq.1)then
                   vvalue= vdset_double(x,y,z)
                   pvalue= pdset_double(x,y,z)
                ELSE
                   vvalue= dble(vdset_float(x,y,z))
                   pvalue= dble(pdset_float(x,y,z))
                ENDIF

                IF(pvalue.eq.0.0)then ! i.e. rock
                   zvel_double_tot(z)=zvel_double_tot(z)+vvalue
                   IF(((z.ge.A).and.(z.le.B)))then !inside porous medium z-range
                      rcount=rcount+1.
                      vsample_tot=vsample_tot+vvalue
                      massflow_tot(z)=massflow_tot(z)+(pvalue*vvalue) ! = +0 , but write it like this for consistency
                   ENDIF
                ELSE ! i.e. pore
                   filling_double(z)=filling_double(z)+1
                   zvel_double_tot(z)=zvel_double_tot(z)+vvalue
                   zvel_double(z)=zvel_double(z)+vvalue ! add summatiomn also for rock sites in corresp. IF loop, i.e. zvel_double_tot
                   density_double(z)=density_double(z)+pvalue
                   pressure_double(z)=pressure_double(z)+pvalue
                   massflow_double(z)=massflow_double(z)+(pvalue*vvalue)
                   massflow_tot(z)=massflow_tot(z)+(pvalue*vvalue) ! = +0 , but write it like this for consistency
                   IF(((z.ge.A).and.(z.le.B)))then
                      porsample=porsample+1.0
                      vsample=vsample+vvalue
                      vsample_tot=vsample_tot+vvalue
                      psample=psample+pvalue
                      mflow_mean=mflow_mean+(pvalue*vvalue)
                   ENDIF
                   ! look for maxima, minima of rho, v_z:
                   IF(pmax.lt.pvalue)then
                      pmax=pvalue
                   ENDIF
                   IF(pmin.gt.pvalue)then
                      pmin=pvalue
                   ENDIF
                   IF(vmax.lt.vvalue)then
                      vmax=vvalue
                   ENDIF
                   IF(vmin.gt.vvalue)then
                      vmin=vvalue
                   ENDIF
                ENDIF
             end do
          end do
       end do
       !DEBUG    write(*,*) 'vvalue at z=87: ', vdset_double(:,:,87)
       !DEBUG    write(*,*) 'zvel at z=87: ', zvel_double(87)

       vsample_tot=vsample_tot/(nx-(2*chanx))/(ny-(2*chany))/(B-A+1)  ! average velocity in p.m., including rock sites
       vsample=vsample/porsample                                      ! "", exlcuding rock sites:
       psample=psample/porsample
       mflow_mean=mflow_mean/porsample

       ! calculate average from sums:
       do z=1,nz
          IF(filling_double(z).ne.0.0)then
             zvel_double(z)=zvel_double(z)/filling_double(z)
             zvel_double_tot(z)=zvel_double_tot(z)/(nx-(2*chanx))/(ny-(2*chany))!/(nz)
             massflow_tot(z)=massflow_tot(z)/(nx-(2*chanx))/(ny-(2*chany))! /(nz)
             density_double(z)=density_double(z)/filling_double(z)
             pressure_double(z)=pressure_double(z)/filling_double(z)
             massflow_double(z)=massflow_double(z)/filling_double(z)
             filling_double(z)=filling_double(z)/(nx-(2*chanx))/(ny-(2*chany))
             ! write out z-cuts
             ! cave: new (11/2007) massflow (total instead of pores-only)
             write (10,FMT='(1X, I3.3, E14.6, E14.6, E14.6, E14.6, E14.6, E14.6)'), z, zvel_double(z), zvel_double_tot(z), &
             density_double(z), (1000*(density_double(z)/0.7)) , massflow_tot(z), filling_double(z)
          !DEBUG          write(10,*), zvel_double(z)
       ENDIF
    end do
    
    write(10,*)'!Sample sites visited= ', (rcount+porsample)
    write(10,*)'!Rock sites visited=   ', rcount
    write(10,*)'!Pore sites visited=   ', porsample
    porsample=porsample/(b-a+1)/(nx-(2*chanx))/(ny-(2*chany))

    vz_P=c1*vsample
    vz_P_tot=vsample_tot*c1
    dynvis=(1000*psample/0.7)*c1*lattice*(2.*tau -1)/6.
    deltap_p=(density_double(a)-density_double(b))*1000*1470*1470/0.7
    c = (2.*tau -1)*(b-a+1)*psample/2.
    rhosample_p=1000*psample/0.7
    delta_rho_rel = (density_double(a)-density_double(b))/psample
    perm=vz_P*dynvis*(b-a+1)*lattice/deltap_p
    !Output permeability in micrometres:
    perm_tot=1e12*vz_P_tot*dynvis*(b-a+1)*lattice/deltap_p
    !perm_tot=vz_P_tot*dynvis*(b-a+1)*lattice/deltap_p
    perm_L=vsample*((2.*tau -1)/6.)*(b-a+1)*lattice*psample/(density_double(a)-density_double(b))
    perm_au=vsample*(b-a+1)/(delta_rho_rel*2.)
    perm_mm=vz_P*dynvis*(b-a+1)*lattice/((pmax-pmin)*1000*1470*1470/psample)

    perm_mu= (vsample_tot*psample*dx*dx*(2.0*tau-1)*(b-a+1) / (density_double(a)-density_double(b)) ) * (1.0/2.0)

    !perm_mu= 0.5 / (density_double(a)/(vsample_tot*psample*dx*dx*(2.0*tau-1)*(b-a+1)) - density_double(b)/(vsample_tot*psample*dx*dx*(2.0*tau-1)*(b-a+1)) )

    
    write(10,*)'!tau= ', tau
    write(10,*)'!dim= ', nx,ny,nz
    write(10,*)'!chanx ,y= ', chanx, chany
    write(10,*)'!range_z= ', A,B
    write(10,*)'!lattice= ', lattice
    write(10,*)'!dynvis= ', dynvis
    write(10,*)'!--------------------------------'
    write(10,*)'!sample_por=', porsample
    write(10,*)'!sample_rhomin,max= ', pmin,pmax !i.e., in pores
    write(10,*)'!sample_vmin,vmax= ', vmin,vmax !i.e., in pores
    write(10,*)'!sample_vz_L=', vsample
    write(10,*)'!vsample_tot=', vsample_tot
    write(10,*)'!sample_rho_L=', psample
!    write(10,*)'!sample_perm_L [l.u.]=',perm_L 
    write(10,*)'!sample_vz_P=', vz_P
    write(10,*)'!sample_perm_tot=', perm_tot
    write(10,*)'!rhosample_p=', rhosample_p
    write(10,*)'!deltap_p=', deltap_p
!    write(10,*)'!sample_perm_P [D]=', perm*1e12
!    write(10,*)'!sample_perm_mm [D]=', perm_mm*1e12
    write(10,*)'!sample_mflow_L =', mflow_mean
    write(10,*)'!--------------------------------'
    write(10,*)'!perm_au = ', perm_au
    write(10,*)'!perm[mu^2] = ', perm_mu

 	close(10)

END program zvelocity

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

