      PROGRAM threshold_greyscale
      
      IMPLICIT NONE
      INTEGER,PARAMETER   :: th=127     ! threshold value
      INTEGER,PARAMETER   :: M=100      ! size of the input grid MxMxM
      INTEGER,PARAMETER   :: N=M        ! size of the output grid NxNxN
      
      INTEGER,DIMENSION(M,M,M)  :: g
      INTEGER                   :: i,j,k,ind,p

      ! Reading the greyscale data
      OPEN(10,FILE ='zeta0.raw.ascii',STATUS='OLD')
      OPEN(20,FILE ='zeta0.raw.ascii.thresh')
    
      p=0
      DO k=1,M
        DO j=1,M
          DO i=1,M
            READ(10,*)ind
            IF (ind <= th)THEN
              g(i,j,k)=0
              p=p+1
            ELSE
              g(i,j,k)=1
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      WRITE(*,*)'Porosity: ',FLOAT(p)/M/M/M

      DO k=1,N
        DO j=1,N
          DO i=1,N
            WRITE(20,101)g(i,j,k)
          ENDDO
        ENDDO
      ENDDO       

 101  FORMAT(I1)

      END PROGRAM threshold_greyscale
