      SUBROUTINE CMP3OUT(JCMN,IM,UTDELTA,VTDELTA,YM,NDT,FRY)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

C
C  PURPOSE:
C      CMP3OUT ACCUMULATES GLOBAL AND ZONAL AVERAGES, CALCULATES A
C      VARIETY OF OUTPUT QUANTITIES, AND PRODUCES SELECTED PRINTER
C      OUTPUT.
C  AUTHOR
C      STEVE POHORSKY    INFORMATICS     TASK 605    OCT 81
C  FOR
C      JIM POLLACK
C  ENVIRONMENT
C      Cray-2            UNICOS 3.0      FORTRAN
C  REVISION HISTORY
C     4-82  SP   MODIFICATIONS MADE AS IN NOTE OF 3-8-82.
C      Re-written September 1993.
C      Removed references to all commons and created an argument list
C      to pass all variables.
C  OUTPUT PARAMETERS
C      FRY(J,L) - The rate at which heat is being gained by the surface
C                 in exchange with the atmosphere through thermal
C                 radiation and convection.
C  CALLED BY
C      COMP3
C
      use grid_h
      use defines_h
      implicit none
C
C######################################################################

      REAL*8 NDT

      REAL*8 FKSX(L_LEVELS), FKSY(L_LEVELS)
      REAL*8 UTDELTA(L_LAYERS), VTDELTA(L_LAYERS), YM(L_LEVELS)
      REAL*8 FRY(L_JSIZE,L_LAYERS)

      integer :: k, l, jcmn, im

C#=====================================================================

      FKSX(3) = 0.0
      FKSY(3) = 0.0
      FKSX(5) = UTDELTA(1)*YM(4)/NDT
      FKSY(5) = VTDELTA(1)*YM(4)/NDT

C     Velocity change and mass

      DO L = 2, L_LAYERS
        K       = 2*L+3
        FKSX(K) = UTDELTA(L)*YM(K-1)/NDT+FKSX(K-2)
        FKSY(K) = VTDELTA(L)*YM(K-1)/NDT+FKSY(K-2)
      END DO


      DO L=1,L_LAYERS
        FRY(JCMN,L) = FRY(JCMN,L)+(FKSY(2*L+3)-FKSY(2*L+1))/FLOAT(IM)
      END DO

      RETURN
      END
