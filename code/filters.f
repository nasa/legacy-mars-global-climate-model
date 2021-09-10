C*********************************************************************
C********************* SUBROUTINE FILTML *****************************
C*********************************************************************


!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center


      SUBROUTINE FILTML(A,B,ULATS,IM,JM, IE, IW, C,D, ORDER, TAU)
      
      IMPLICIT NONE

      INTEGER     ORDER, PASSES
      REAL TAU, SHAPFAC, ZERO, ONE, TWO
      
      PARAMETER   (ZERO=0.0, ONE=1.0, TWO=2.)

C   ARGUMENTS

      REAL    A(IM,JM), B(IM,JM), C(IM,JM), D(IM,JM)
      INTEGER    IE(IM*JM), IW(IM*JM)
      LOGICAL      ULATS
      INTEGER IM, JM

      INTEGER I, K, LEN, JJ

      DO 8090 I=1,IM*JM

      IF(MOD(I,IM).NE.0) THEN
       IE(I) = I + 1
      ELSE
       IE(I) = I + 1 - IM
      ENDIF

      IF(MOD(I,IM).NE.1) THEN
       IW(I) = I - 1
      ELSE
       IW(I) = I - 1 + IM
      ENDIF

8090  CONTINUE

      PASSES=ORDER/2
      SHAPFAC=(-.25)**PASSES

      IF(ULATS) THEN
        LEN = IM*(JM-1)
        JJ = 1
      ELSE
        LEN = IM*(JM-2)
        JJ = 2
      ENDIF

C =====>  DO ZONAL SMOOTHING

      DO I =1,LEN
       C(I,JJ) = A(I,JJ)
      ENDDO

      DO 100 K=1,PASSES/2
        DO I =1,LEN
         B(I,JJ) = C(IE(I),JJ) - TWO*C(I,JJ) + C(IW(I),JJ)
        ENDDO
        DO I =1,LEN
         C(I,JJ) = B(IE(I),JJ) - TWO*B(I,JJ) + B(IW(I),JJ)
        ENDDO
100   CONTINUE

C   SAVE ZONALLY SMOOTH FIELD IN D

      DO I =1,LEN
       D(I,JJ) = -C(I,JJ) * SHAPFAC
       C(I,JJ) =  A(I,JJ) + D(I,JJ)
      ENDDO

C =====>   DO MERIDIONAL SMOOTHING

      IF(ULATS) THEN
        DO I =1,IM
          B(I,1 ) = ZERO
          B(I,JM) = ZERO
        ENDDO
      ELSE
        DO I =1,IM
          B(I,1 ) = ZERO
          B(I,JM-1) = ZERO
        ENDDO
      ENDIF

      DO 101 K=1,PASSES

      IF (ULATS) THEN
        DO I =1,LEN-IM
          B(I,2) = C(I,2) - C(I,1)
        ENDDO
        DO I =1,LEN
          C(I,1) = B(I,2) - B(I,1)
        ENDDO
      ELSE
        DO I =1,LEN-IM
          B(I,2) = C(I,3) - C(I,2)
        ENDDO
        DO I =1,LEN
          C(I,2) = B(I,2) - B(I,1)
        ENDDO
      ENDIF

101   CONTINUE


C    TAU IS THE DAMPING TIME IN SECONDS.

      DO I =1,LEN
C Zonal filtering only
C        B(I,JJ) = D(I,JJ) * (ONE/TAU)
C Zonal and meridional
        B(I,JJ) = (D(I,JJ) - C(I,JJ)*SHAPFAC) * (ONE/TAU)
      ENDDO

      RETURN
      END



C*********************************************************************
C********************* SUBROUTINE FILTPQ *****************************
C*********************************************************************


      SUBROUTINE FILTPQ(A, B, IM,JM, IE,IW, C,D, ORDER, TAU)
      
      IMPLICIT NONE

      INTEGER     ORDER, PASSES
      REAL        SHAPFAC, TAU

      REAL FOUR_THIRDS
      PARAMETER   (FOUR_THIRDS = 4./3.)

      REAL ONE
      PARAMETER   (ONE = 1.)

      REAL ONE_THIRD
      PARAMETER   (ONE_THIRD = 1./3.)

C   ARGUMENTS

      REAL         A(IM,JM)
      INTEGER     IE(IM*JM), IW(IM*JM)
      INTEGER      IM, JM

C   WORK SPACE

      REAL         B(IM,JM), C(IM,JM), D(IM,JM)
      INTEGER      I, K, J
      REAL         CSP, CNP
      REAL P8SSUM

CFPP$ EXPAND (P8SSUM)
      DO 8090 I=1,IM*JM

      IF(MOD(I,IM).NE.0) THEN
       IE(I) = I + 1
      ELSE
       IE(I) = I + 1 - IM
      ENDIF

      IF(MOD(I,IM).NE.1) THEN
       IW(I) = I - 1
      ELSE
       IW(I) = I - 1 + IM
      ENDIF

8090  CONTINUE

      PASSES=ORDER/2
      SHAPFAC=(-.25)**PASSES

C =====>  DO ZONAL SMOOTHING

      DO  I=1,IM*(JM-1)
       B(I,1) = A(IE(I),1) - 2.*A(I,1) + A(IW(I),1)
      ENDDO
       
      DO I=1,IM*(JM-1)
       C(I,1) = B(IE(I),1) - 2.*B(I,1) + B(IW(I),1)
      ENDDO
      
      DO K=1,PASSES/2 - 1
      
       DO  I=1,IM*(JM-1)
        B(I,1) = C(IE(I),1) - 2.*C(I,1) + C(IW(I),1)
       ENDDO
       
       DO I=1,IM*(JM-1)
        C(I,1) = B(IE(I),1) - 2.*B(I,1) + B(IW(I),1)
       ENDDO
       
      ENDDO

C   SAVE ZONAL SMOOTHING TENDENCY IN D

      DO I=1,IM*(JM-1)
       D(I,1) = -C(I,1) * SHAPFAC
       C(I,1) =  A(I,1) + D(I,1)
      ENDDO

C =====>   DO MERIDIONAL SMOOTHING


      DO K=1,PASSES
      
       CSP  = FOUR_THIRDS*P8SSUM(C(1,   1),IM)*(1./FLOAT(IM))
     *      - ONE_THIRD  *P8SSUM(C(1,   2),IM)*(1./FLOAT(IM))
       CNP  = FOUR_THIRDS*P8SSUM(C(1,JM-1),IM)*(1./FLOAT(IM))
     *      - ONE_THIRD  *P8SSUM(C(1,JM-2),IM)*(1./FLOAT(IM))

       DO I=1,IM
        B(I, 1) = C(I,1) - CSP
        B(I,JM) = CNP - C(I,JM-1)
       ENDDO
       
       DO I=1,IM*(JM-2)
        B(I,2) = C(I,2) - C(I,1)
       ENDDO
       
       DO I=1,IM*(JM-1)
        C(I,1) = B(I,2) - B(I,1)
       ENDDO
       
      ENDDO

C    TAU IS THE DAMPING TIME IN SECONDS.

      DO I=1,IM*(JM-1)
C Zonal filtering only
C       B(I,1) = D(I,1) * (ONE/TAU)
C Zonal and meridional
       B(I,1) = (D(I,1) - C(I,1) * SHAPFAC) * (ONE/TAU)
      ENDDO


      RETURN
      END

c...  begin p8ssum code
      FUNCTION P8SSUM(A,N)
      REAL A(N)
      P8SSUM = SSUM(N,A,1)
      RETURN
      END
c...  end p8ssum code
