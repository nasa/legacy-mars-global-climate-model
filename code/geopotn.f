      SUBROUTINE GEOPOTN(PTROP,SIGMA,P,T,DSIG,NLAY,TOPOG,GEOP)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

C
C  PURPOSE
C      GEOPOTN CALCULATES THE GEOPOTENTIAL (GRAV. TIMES HEIGHT) AT THE
C      MIDPOINTS OF EACH LAYER FOR EACH 'PI' POINT.
C  AUTHOR
C      STEVE POHORSKY    INFORMATICS     TASK 605    NOV 81
C  FOR
C      JIM POLLACK
C  ENVIRONMENT
C      Cray-2            UNICOS 3.0      FORTRAN
C
C  REVISION HISTORY
C      Re-written September 1993.
C      Removed references to all commons and created an argument list
C      to pass all variables.
C
C  INPUT PARAMETERS
C      P(J,I) ARRAY   - THE CURRENT 'PI' VALUE (SURFACE PRESSURE
C                             MINUS TROPOPAUSE PRESSURE) AT THE 'PI'
C                             POINTS.
C      T(J,I,K) ARRAY - THE TEMPERATURE AT THE MIDPOINT OF EACH
C                             LAYER FOR THE 'PI' POINTS. (K = 1 TO NLAY)
C  OUTPUT PARAMETERS
C      GEOP(J,I,K) ARRAY   - GEOPOTENTIAL (GRAV. TIMES HEIGHT) AT THE
C                            MIDPOINTS OF EACH LAYER FOR EACH 'PI'
C                            POINT. THE HEIGHT IS MEASURED FROM SOME
C                            REFERENCE ELEVATION.
C  .........NOTE GEOP IS LEVEL 2..............
C  CALLED BY
C      COMP3
C
      use grid_h
      use defines_h
      use constants_h, only: RGAS, CP, KAPA

      implicit none
C
C#######################################################################

      REAL*8 P(L_JSIZE,L_ISIZE),     T(L_JSIZE,L_ISIZE,L_LAYERS)
      REAL*8 SIGMA(L_LEVELS),        DSIG(L_LAYERS)
      REAL*8 TOPOG(L_JSIZE,L_ISIZE), GEOP(L_JSIZE,L_ISIZE,L_LAYERS)
      REAL*8 DPHIL(L_LEVELS),        TL(L_LEVELS)
      REAL*8 PL(L_LEVELS)

      integer :: i, j, k, l, nlay
      real*8  :: sum1, sum2, ptrop, rkdn, ptem, altrm

C#======================================================================

C     Loop over all grid points.

c     DO 1050 J = 1,L_JSIZE
C Only go from 1 to JM-1 in the C grid
      DO 1050 J = 1,L_JSIZE-1

        DO 1000 I = 1,L_ISIZE

C         Pressure and temperature at the layer midpoints

          DO 100 K = 4, L_LEVELM1, 2
            PL(K) = PTROP+SIGMA(K)*P(J,I)
  100     CONTINUE

          DO 200 K = 1, NLAY
            TL(2*K+2) = T(J,I,K)
  200     CONTINUE

C         GEOPOTENTIAL IS CALCULATED FROM AN ENERGY-CONSERVING FORM OF
C         THE HYDROSTATIC EQUATION.  A PRESENTATION OF THIS METHOD IS
C         FOUND IN THE SET OF WORKSHOP NOTES 'THE UCLA ATMOSPHERIC
C         GENERAL CIRCULATION MODEL' BY A. ARAKAWA AND Y. MINTZ OF
C         MARCH 25, 1974. (SEE ESPECIALLY EQUATIONS III.48 AND III.49.)

          DO 300 K = 4, L_LEVELM3, 2
            RKDN  = (PL(K+2)/PL(K))**KAPA
            PTEM  = RKDN*TL(K)/TL(K+2)
            ALTRM = TL(K)

!  PGI compiler would (sometimes) barf on "IF(PTEM .NE. 1.0)  THEN"
!  and end up dividing by zero.  The IF test modified to get around the
!  problem.   10/18/01

c           IF(PTEM .NE. 1.0)  THEN
            IF(PTEM.lt.0.99999 .or. PTEM.gt.1.0001) THEN
              ALTRM = ALTRM*LOG(PTEM)/(PTEM-1.0)
            ENDIF

C           Differential in geopotential (PHI) between adjacent midpoint

            DPHIL(K) = CP*ALTRM*(RKDN-1.0)
  300     CONTINUE

          SUM1 = 0.0
          SUM2 = 0.0

          DO 400 K = 1, NLAY
            L    = 2*K+2
            SUM1 = SUM1+SIGMA(L)*TL(L)*DSIG(K)/PL(L)

            IF (K.EQ.NLAY)  GOTO 400

            SUM2 = SUM2+SIGMA(L+1)*DPHIL(L)
  400     CONTINUE

          GEOP(J,I,NLAY)=-TOPOG(J,I)+P(J,I)*RGAS*SUM1-SUM2

          DO 500 L = 1, L_LAYERM1
            K = L_LAYERM1 - L + 1
            GEOP(J,I,K) = GEOP(J,I,K+1) + DPHIL( 2*K+2 )
  500     CONTINUE

 1000   CONTINUE
 1050 CONTINUE

      RETURN
      END
