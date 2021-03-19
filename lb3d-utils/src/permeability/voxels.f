! voxels.f: discretizes a file containing sphere positions and radii.
! Originally by B. Biswal
! Modified by Jens, 03/2008

      MODULE parameters
 
      ! Number of voxels
      INTEGER    :: voxels
      ! Number of spheres
      INTEGER    :: inc
      ! Box dimension
      REAL(8)    :: Lx,Ly,Lz

      ! input file "parameters.dat"
      ! Compile and run "ifort voxels.f -o voxels;./voxels"
      ! Output stored in "voxels.dat"

      REAL(8)    :: xb, yb, zb
      REAL(8)    :: Lxb,Lyb,Lzb

      INTEGER,PARAMETER    :: maxpts=216
      REAL(8),DIMENSION(maxpts)       :: xi,yi,zi

      REAL(8)       :: ar
      REAL(8)       :: ahalf
      REAL(8)       :: xbg, ybg, zbg
      REAL(8)       :: xnd, ynd, znd
      INTEGER    :: xmax
      INTEGER    :: ymax
      INTEGER    :: zmax

      REAL(8),DIMENSION(:),ALLOCATABLE    :: xc,yc,zc,rc
      INTEGER,DIMENSION(:),ALLOCATABLE    :: grx,gry,grz
      REAL(8),DIMENSION(:),ALLOCATABLE    :: rxb,ryb,rzb,rxe,rye,rze
      INTEGER    :: jj=10

      INTEGER                             :: in,id
      REAL(8)                             :: xx,yy,zz
      
      CONTAINS
      SUBROUTINE CollocationPoints
      INTEGER :: i,j,k,n,ind,step,m
      REAL(8) :: aq
      m=(maxpts)**(1./3)
      n = 0 ; ind = m-1 ;aq=1.*ar/2/m
      DO k=-ind,ind,2
         DO j=-ind,ind,2
            DO i=-ind,ind,2
               n=n+1
               xi(n)=i*aq
               yi(n)=j*aq
               zi(n)=k*aq
            END DO
         END DO
      END DO
      END SUBROUTINE Collocationpoints
      
      END MODULE parameters

!--------------------------------------------------------------------

      PROGRAM BinaryImage
      USE parameters
      IMPLICIT NONE
      INTEGER :: i,j,AllocStat
      REAL(8)    :: ax,ay,az,dia
      CHARACTER*80 :: INFILE,OUTFILE
      
      PRINT *, 'Input file?'
      READ(*,*) INFILE
      PRINT *, 'Output file?'
      READ(*,*) OUTFILE
      PRINT *, 'Number of spheres?'
      READ(*,*) inc
      PRINT *, 'Number of voxels?'
      READ(*,*) voxels
      PRINT *, 'Box dimension?'
      READ(*,*) Lx
      !Lx=4.44E-6
      Ly=Lx
      Lz=Lx
      xb=0.
      yb=0.
      zb=0.
      Lxb=Lx+xb
      Lyb=Ly+yb
      Lzb=Lz+zb
      ar=Lx/voxels
      ahalf=ar/2
      xbg=xb
      ybg=yb
      zbg=zb
      xnd=Lxb
      ynd=Lyb
      znd=Lzb
      xmax=(xnd-xbg)/ar
      ymax=(ynd-ybg)/ar
      zmax=(znd-zbg)/ar


      PRINT *, 'OUTFILE', OUTFILE
      PRINT *, 'INFILE', INFILE
      OPEN (11,FILE=INFILE,STATUS='OLD')
      OPEN (15,FILE=OUTFILE)

      CALL CollocationPoints

      ALLOCATE(xc(inc),yc(inc),zc(inc),rc(inc),
     &     grx(inc),gry(inc),grz(inc),rxb(inc),ryb(inc),rzb(inc),
     &     rxe(inc),rye(inc),rze(inc), STAT = AllocStat)

      DO i=1,inc
         READ(11,*)xc(i),yc(i),zc(i),rc(i)
         dia=rc(i)
         rxb(i)=xc(i)+dia; rxe(i)=xc(i)-dia 
         ryb(i)=yc(i)+dia; rye(i)=yc(i)-dia 
         rzb(i)=zc(i)+dia; rze(i)=zc(i)-dia
      ENDDO
      CLOSE(11)
      CALL Discretization
      WRITE(*,*)' '

      END PROGRAM BinaryImage

!--------------------------------------------------------------------

      SUBROUTINE Discretization
      USE parameters
      IMPLICIT none
      INTEGER :: i,j,k,ic,nfac,count,inx,iny,inz,ip
      REAL(8) :: x1,x2,y1,y2,z1,z2,x,y,z

      WRITE(*,*) 'Discretizing'
      nfac = zmax/100
      IF(nfac==0)nfac=1
      WRITE(*,*)xmax,ymax,zmax
      DO k=0,zmax-1
        z=zbg+ahalf+(k*ar)
        IF (MOD(k,nfac)==0)WRITE(*,*)'DONE',(100.*k)/zmax,'%'
        inz=0; z1=z-ahalf; z2=z+ahalf
        DO ic=1,inc
          IF ( rzb(ic) >= z1 .AND. rze(ic) <= z2)THEN
            inz=inz+1
            grz(inz)=ic
          ENDIF
        ENDDO
        DO j=0,ymax-1
          y=ybg+ahalf+(j*ar)
          iny=0; y1=y-ahalf; y2=y+ahalf
          DO ic=1,inz
            IF ( ryb(grz(ic)) >= y1 .AND. rye(grz(ic)) <= y2)THEN
              iny=iny+1
              gry(iny)=grz(ic)
            ENDIF
          ENDDO
          DO i=0,xmax-1
            count=0
            x=xbg+ahalf+(i*ar)
            inx=0; x1=x-ahalf; x2=x+ahalf
            DO ic=1,iny
              IF ( rxb(gry(ic)) >= x1 .AND. rxe(gry(ic)) <= x2)THEN
                inx=inx+1
                grx(inx)=gry(ic)
              ENDIF
            ENDDO
            points: DO ip=1,maxpts
              xx=x+xi(ip)
              yy=y+yi(ip)
              zz=z+zi(ip)
              DO ic=1,inx
                id=grx(ic)
                CALL InsideSphere
                IF(in==1) THEN
                  count=count+1
                  cycle points
                ENDIF
              ENDDO
            ENDDO points
            WRITE(15,303) count
          ENDDO
        ENDDO
      ENDDO
      CLOSE(15)
      
 301  FORMAT(I1)
 302  FORMAT(I2)
 303  FORMAT(I3)
      END SUBROUTINE Discretization

!--------------------------------------------------------------------

      SUBROUTINE InsideSphere
      USE parameters, only: xx,yy,zz,xc,yc,zc,rc,id,in
      IMPLICIT none
      REAL(8)              :: dist
      REAL(8),DIMENSION(3) :: di
      in=0
      di(1) = xx - xc(id)
      di(2) = yy - yc(id)
      di(3) = zz - zc(id)
      dist=SQRT(di(1)**2+di(2)**2+di(3)**2)
      IF (dist <= rc(id))in=1
      END SUBROUTINE InsideSphere

!--------------------------------------------------------------------
