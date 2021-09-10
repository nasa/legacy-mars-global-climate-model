      SUBROUTINE GRIDVEL(JCMN,ICMN,U,V,QTRACE,QCOND,UPI,VPI,QPI,QPIG)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

C  PURPOSE
C      GRDVEL CALCULATES THE U AND V VALUES (WIND VELOCITY) AT THE
C      MIDPOINT OF EACH LAYER FOR A GIVEN 'PI' POINT BY AVERAGING THE
C      U AND V VALUES FOR THE NEIGHBORING U,V GRID POINTS.
C
C  AUTHOR
C      STEVE POHORSKY    INFORMATICS     TASK 605    OCT 81
C
C  FOR
C      JIM POLLACK
C
C  ENVIRONMENT
C      Cray-2            UNICOS 3.0      FORTRAN
C
C  REVISION HISTORY
C      Re-written 9/21/93  
C      removed references to all commons and created an argument list
C      to pass all variables.
C
C  INPUT PARAMETERS
C      JCMN & ICMN          - THE J AND I COORDINATES OF THE 'PI' POINT.
C      U(J,I,K) ARRAY       - THE LOCAL EAST-WEST WIND VELOCITY
C                             COMPONENT AT EACH OF THE U,V-GRID POINTS
C                             AT THE MIDPOINTS OF EACH LAYER.
C      V(J,I,K) ARRAY       - THE LOCAL NORTH-SOUTH WIND VELOCITY
C                             COMPONENT AT EACH OF THE U,V-GRID POINTS
C                             AT THE MIDPOINTS OF EACH LAYER.
C
C  OUTPUT PARAMETERS:     (K = 4 TO NLEVSM1 BY 2S.)
C      UPI(K) ARRAY          - U VALUES AT THE 'PI' POINT FOR EACH LAYER
C      VPI(K) ARRAY          - V VALUES AT THE 'PI' POINT FOR EACH LAYER
C
C  CALLED BY
C      COMP3
C
      use grid_h
      use defines_h
      use standard_h, only: microphysics

      implicit none

C######################################################################

      real*8 :: U(L_JSIZE,L_ISIZE,L_LAYERS)
      real*8 :: V(L_JSIZE,L_ISIZE,L_LAYERS)
      real*8 :: UPI(L_LEVELS), VPI(L_LEVELS)
      real*8 :: QTRACE(L_JSIZE,L_ISIZE,L_LAYERS,NTRACE)
      real*8 :: QCOND(L_JSIZE,L_ISIZE,NTRACE)
      real*8 :: QPI(L_LEVELS,NTRACE)
      real*8 :: QPIG(NTRACE) 

!     implicit none

      integer :: JCMN,ICMN, L, K, M

C#=====================================================================

C     CALCULATE WIND VELOCITY AT THE MIDPOINT OF EACH LAYER FOR A
C     GIVEN 'PI' POINT BY AVERAGING THE U AND V VALUES FOR THE
C     NEIGHBORING U,V GRID POINTS.

C     FOR 'PI' POINTS OTHER THAN THE NORTH OR SOUTH POLE

C     L is the level index of the middle of layer K.

      DO 100  K = 1, L_LAYERS
        L      = 2*K+2
        if(ICMN .EQ. 1)THEN
          UPI(L) = 0.5*(U(JCMN,ICMN,K)+U(JCMN,L_ISIZE,K))
        ELSE
          UPI(L) = 0.5*(U(JCMN,ICMN,K)+U(JCMN,ICMN-1,K))
        ENDIF
        VPI(L) = 0.5*(V(JCMN,ICMN,K)+V(JCMN+1,ICMN,K))
        DO M = 1,NTRACE
          QPI(L,M) = QTRACE(JCMN,ICMN,K,M)
        END DO
  100 CONTINUE

!     if micro-physics is on, qcond has valid values.  If we are not
!     doing micro-physics, tracers are not updated, and so set to zero.

      if(MICROPHYSICS .eqv. .TRUE.) then
        DO M = 1,NTRACE
          QPIG(M) = QCOND(JCMN,ICMN,M)
        ENDDO 
      else
        DO M = 1,NTRACE
          QPIG(M) = 0.0D0
        ENDDO
      end if

      RETURN
      END
