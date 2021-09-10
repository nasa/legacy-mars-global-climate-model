      SUBROUTINE POTEMP2(JCMN,ICMN,PTROP,PSL,P,SIGMA,TSTRAT,
     *                   AADJ,BADJ,PLOGADJ,PL,OM,TL,TETA)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

!  PURPOSE
!      POTEMP CALCULATES THE POTENTIAL TEMPERATURE AND SEVERAL OTHER
!      QUANTITIES AT EACH LEVEL FOR A GIVEN 'PI' GRID POINT.
!
!  AUTHOR
!      STEVE POHORSKY    INFORMATICS     TASK 605    OCT 81
!
!  FOR
!      JIM POLLACK
!
!  ENVIRONMENT
!      Cray-2            UNICOS 3.0      FORTRAN
!
!  REVISION HISTORY
!     1-82  SP   SCHEME FOR POTENTIAL TEMPERATURE INTERPOLATION REVISED
!                ACCORDING TO NOTE OF 1-18-82.
!     6-82  SP   TL(1), AND TL(2) ADDED AS OUTPUT PARAMETERS.
!
!      Re-written September 1993.
!      Removed references to all commons and created an argument list
!      to pass all variables.
!  INPUT PARAMETERS
!      JCMN & ICMN          - THE J AND I COORDINATES OF THE 'PI' POINT.
!      P(JCMN,ICMN)         - THE CURRENT 'PI' VALUE (SURFACE PRESSURE
!                             MINUS TROPOPAUSE PRESSURE) AT THE 'PI'
!                             POINT.
!      SIGMA(K) ARRAY       - THE SIGMA (VERTICAL COORDINATE) FOR EACH
!                             LEVEL OF THE ATMOSPHERE. (K = 3 TO
!                             NLEVELS)
!      AADJ(K) AND  -  THE FIRST TWO COEFFICIENTS IN THE TAYLOR'S SERIES
!      BADJ(K) ARRAYS  FOR THE DIFFERENCE IN LOG OF PRESSURE BETWEEN
!                      ADJACENT LAYERS. AADJ AND BADJ ARE, RESPECTIVELY,
!                      THE FUNCTION AND ITS DERIVATIVE EVALUATED AT PSL.
!                      ( K = 3, NLEVSM1. )
!                      (THE DIFFERENCE BETWEEN LOG OF PRESSURE AT LEVEL
!                      K+1 AND LEVEL K IS A FUNCTION OF THE SURFACE
!                      PRESSURE. THIS FUNCTION IS EVALUATED AT A GIVEN
!                      SURFACE PRESSURE IN THIS ROUTINE BY MEANS OF ITS
!                      TAYLOR SERIES INSTEAD OF DIRECTLY TO MAKE
!                      EXECUTION FASTER. THE TAYLOR SERIES APPROXIMATION
!                      USED HERE IS
!                      THE DIFFERENCE BETWEEN LOG OF PRESSURE AT LEVEL
!                      K+1 AND LEVEL K =  AADJ(K) + BADJ(K) * DELTA P
!                      WHERE DELTA P IS THE GIVEN SURFACE PRESSURE
!                      MINUS PSL.)
!      T(JCMN,ICMN,K)       - THE TEMPERATURE AT THE MIDPOINT OF EACH
!        ARRAY ELEMENTS       LAYER FOR THE 'PI' POINT. (K = 1 TO NLAY)
!      PTROP                - PRESSURE AT THE TROPOPAUSE.
!
!  OUTPUT PARAMETERS:     (K = 3 TO NLEVELS UNLESS STATED OTHERWISE)
!      PL(K) ARRAY          - PRESSURE AT EACH LEVEL FOR THE 'PI' POINT.
!      TL(K) ARRAY          - TEMPERATURE AT EACH LEVEL FOR THE 'PI'
!                             POINT. (FOR K = 1 TO NLEVELS.)
!      TETA(K) ARRAY        - POTENTIAL TEMPERATURE AT EACH LEVEL.
!      OM(K) ARRAY          - CONVERSION FACTOR FROM TETA TO TL AT EACH
!                             LEVEL FOR THE 'PI' POINT.
!      PLOGADJ(K) ARRAY     - DIFFERENCE IN LOG OF PRESSURE BETWEEN
!                             ADJACENT LEVELS. PLOGADJ(K) = LN(PL(K+1))-
!                             LN(PL(K)). (K = 3 TO NLEVELS - 1)
!
!  CALLED BY
!      COMP3
!
      use grid_h
      use defines_h
      use constants_h, only: KAPA
      implicit none

!######################################################################

      REAL*8 :: P(L_JSIZE,L_ISIZE), SIGMA(L_LEVELS), AADJ(L_LEVELS)
      real*8 :: BADJ(L_LEVELS), TSTRAT(L_JSIZE,L_ISIZE)
      real*8 :: T(L_JSIZE,L_ISIZE,L_LAYERS), PL(L_LEVELS), OM(L_LEVELS)
      real*8 :: TL(L_LEVELS), TETA(L_LEVELS), PLOGADJ(L_LEVELS)

!     implicit none

      integer :: K, L, JCMN,ICMN
      real*8  :: ptrop, psl, dpsurf, slope

!#=====================================================================

!     CALCULATE PRESSURE AT ALL LEVELS IN TROPOSPHERE FROM
!     PI & SIGMA VALUES.

!     At the surface

      PL(L_LEVELS) = PTROP+P(JCMN,ICMN)
      OM(L_LEVELS) = 1.0

!     OM is the conversion factor from TETA to T.

      DO 100  K = 3, L_LEVELM1
        PL(K) = PTROP+SIGMA(K)*P(JCMN,ICMN)
        OM(K) = (PL(K)/PL(L_LEVELS))**KAPA
  100 CONTINUE

      PL(2) = PTROP/2.0
      OM(2) = (PL(2)/PL(L_LEVELS))**KAPA

!     Potential temperature at the layer midpoints.

      DO 200 K=2,L_LEVELM1,2
        TETA(K) = TL(K)/OM(K)
  200 CONTINUE

!     CALCULATE DIFFERENCES IN LOG OF PRESSURE BETWEEN ADJACENT
!     LEVELS FOR USE IN EXTRAPOLATION.

!     Delta in surface pressure

      DPSURF = PL(L_LEVELS)-PSL

!     Use a Taylor series expression for LOG(PL(K+1))-LOG(PL(K))

      DO 300  K = 3, L_LEVELM1
        PLOGADJ(K) = AADJ(K)+BADJ(K)*DPSURF
  300 CONTINUE

!     INTERPOLATE TO FIND POTENTIAL TEMPERATURE AT LAYER BOUNDARIES.

      DO 400 K=3,L_LEVELM4,2
        SLOPE   = (TETA(K+1)-TETA(K-1))/(PLOGADJ(K)+PLOGADJ(K-1))
        TETA(K) = TETA(K-1)+SLOPE*PLOGADJ(K-1)
  400 CONTINUE

!     Extrapolate from above two midpoints.

      K       = L_LEVELM2
      SLOPE   = (TETA(K-1)-TETA(K-3))/(PLOGADJ(K-2)+PLOGADJ(K-3))
      TETA(K) = TETA(K-1)+SLOPE*PLOGADJ(K-1)

!     Extrapolate for the surface from the 2 levels immediately above
!     for a good approximation to the physics.

      K       = L_LEVELS
      SLOPE   = (TETA(K-1)-TETA(K-2))/PLOGADJ(K-2)
      TETA(K) = TETA(K-1)+SLOPE*PLOGADJ(K-1)

!     Temperatures at the layer boundaries.

      DO 500  K = 3, L_LEVELS, 2
        TL(K) = OM(K)*TETA(K)
  500 CONTINUE

      RETURN
      END
