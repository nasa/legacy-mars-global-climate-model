      SUBROUTINE COLDAIR(PTROP,TSTRAT,NDT,SIGMA,PL,TL,DSIG,YM,IM,GT,
     *                   CO2ICE,SQRDY,ZIN,CO2LATST,CO2LAT,
     *                   DMADT,ATMCOND,JCMN,ICMN)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center
!
!  PURPOSE
!      COLDAIR CHECKS TO SEE IF THE ATMOSPHERE IS FREEZING. FOR EACH
!      GRID POINT COLDAIR DETERMINES IF THE ATMOSPHERIC TEMPERATURE
!      HAS DROPPED BELOW THE CO2 FROST POINT. IF SO, THE AMOUNT OF
!      CO2 CONDENSATION IS CALCULATED AND APPROPRIATE ADJUSTMENTS
!      ARE MADE.
!  AUTHOR
!      STEVE POHORSKY    INFORMATICS     TASK 605    FEB 82
!
!  FOR
!      JIM POLLACK      (SEE NOTE OF 2-1-82)
!
!  ENVIRONMENT
!        Cray-2         UNICOS 3.0         FORTRAN
!
!  REVISION HISTORY
!      JB WHITE         TASK 904           NOV 85
!      STRUCTURE CODE
!      INCORPORATE STRATOSPHERE (SEE NOTE OF 15 SEP 85)
!
!      JR SCHAEFFER     TASK 904           JAN 86
!      STRATOSPHERE INCORPORATED
!      JR SCHAEFFER     TASK 904           DEC 86
!      THERMAL INERTIA DATA/CODE CHANGES INCORPORATED
!                              (SEE NOTE OF 10/6/86).
!
!      Re-written September 1993.
!      Removed references to all commons and created an argument list
!      to pass all variables.
!  INPUT PARAMETERS
!      GT(J,I) ARRAY  - GROUND TEMPERATURE OR (-1) TIMES CO2 ICE
!                             MASS AT THE 'PI' POINTS.
!      P(J,I) ARRAY   - THE CURRENT 'PI' VALUE (SURFACE PRESSURE
!                             MINUS TROPOPAUSE PRESSURE) AT THE 'PI'
!                             POINTS.
!      T(J,I,K) ARRAY - THE TEMPERATURE AT THE MIDPOINT OF EACH
!                             LAYER FOR THE 'PI' POINTS. (K = 1 TO NLAY)
!
!  OUTPUT PARAMETERS
!      T(J,I,K) ARRAY - WHERE THE AIR TEMPERATURE IS FOUND TO BE BELOW
!                       THE CO2 FROST POINT, IT IS SET TO THE CO2 FROST
!                       POINT.
!      GT(J,I) ARRAY  - AS IN INPUT, BUT ADJUSTED ACCORDING TO THE
!                       EFFECT OF ANY CONDENSATION IN THE AIR.
!      DMADT(J,I)     - THE RATE OF CONDENSATION PER UNIT SURFACE AREA
!       ARRAY           IN THE WHOLE ATMOSPHERIC COLUMN AT A 'PI' GRID
!                       POINT.
!  CALLED BY
!      COMP3
!
      use grid_h
      use defines_h
      use constants_h, only: GRAV, CP, SCALEP, XLHTC

      implicit none
!
!######################################################################

      REAL*8  :: CO2LATST(L_JSIZE),      CO2LAT(L_JSIZE,L_LAYERS)
      REAL*8  :: PL(L_LEVELS),           TL(L_LEVELS)
      REAL*8  :: DSIG(L_LAYERS),         SIGMA(L_LEVELS)
      REAL*8  :: GT(L_JSIZE,L_ISIZE),    ZIN(L_JSIZE,L_ISIZE,NL)
      REAL*8  :: CO2ICE(L_JSIZE,L_ISIZE)
      REAL*8  :: DMADT(L_JSIZE,L_ISIZE), TSTRAT(L_JSIZE,L_ISIZE)
      REAL*8  :: YM(L_LEVELS)
      REAL*8  :: NDT
 
      REAL*4  :: ATMCOND(L_JSIZE,L_ISIZE,L_LAYERS)

!     implicit none

      integer :: I, J, K, L, ICMN, JCMN, IM
      real*8  :: TGP, TINP, SQRDY, PIFACT, PTROP, TSAT, YMSTRAT
      real*8  :: ACONDNS, CONDENS, PSAT

!#=====================================================================

      I = ICMN
      J = JCMN

!  This used to loop over all I,J

      CO2LATST(J) = 0.0

      DO 100 L=1,L_LAYERS
         CO2LAT(J,L) = 0.0
  100 CONTINUE

      CONDENS=0.0D0

!     STRATOSPHERE

      PSAT = PTROP/2.0
      TSAT = 3182.48/(23.3494-LOG(PSAT))

!     COMPUTE STRATOSPHERIC CONDENSATION

      IF(TSTRAT(J,I).LT.TSAT) THEN
        YMSTRAT     = PTROP*SCALEP/GRAV
        ACONDNS     = CP*(TSAT-TSTRAT(J,I))*YMSTRAT/XLHTC
        CO2LATST(J) = CO2LATST(J)+ACONDNS*XLHTC/(NDT*FLOAT(IM))
        TSTRAT(J,I) = TSAT
        CONDENS     = CONDENS+ACONDNS/NDT
      END IF

!     TROPOSPHERE

      DO 110 L=1,L_LAYERS
        K    = 2*L+2
        PSAT = PL(K)
        TSAT = 3182.48/(23.3494-LOG(PSAT))

        IF(TL(K).LT.TSAT) THEN
          PIFACT      = PL(L_LEVELS)-PTROP
          YM(K)       = DSIG(L)*PIFACT*SCALEP/GRAV
          ACONDNS     = CP* (TSAT-TL(K)) * YM(K) /XLHTC
          CO2LAT(J,L) = CO2LAT(J,L)+ACONDNS*XLHTC/(NDT*FLOAT(IM))
          TL(K)       = TSAT
          CONDENS     = CONDENS+ACONDNS/NDT
          ATMCOND(J,I,L) = ATMCOND(J,I,L) + ACONDNS
        END IF

  110 CONTINUE

        IF(CONDENS.GT.0.0) THEN

!         CO2 FROST POINT AT THIS SURFACE PRESSURE
          PSAT = PL(L_LEVELS)
          TSAT = 3182.48/(23.3494-LOG(PSAT))

!         CASE 1:  CO2 ICE ALREADY ON GROUND
!         ADD CONDENSATION TO EXISTING CO2 ICE MASS

          IF(CO2ICE(J,I) .GT. 0.0 ) THEN
            CO2ICE(J,I) = CO2ICE(J,I) + NDT*CONDENS
            GT(J,I)     = TSAT

!           CASE 2:  NO CO2 ICE ON GROUND; GROUND IS WARMER

          ELSE

!           GROUND TEMPERATURE DROPS WHEN CO2 ICE SUBLIMES
!           ON WARMER SURFACE

            TINP = SQRDY/ZIN(J,I,1)
            TGP  = GT(J,I) - TINP*NDT*CONDENS*XLHTC

!           COMPUTE HOW MUCH CO2 ICE SUBLIMES ON HITTING GROUND

            IF(TGP.GE.TSAT) THEN

!             CASE 2A:  ALL CO2 ICE SUBLIMES
!             NO NET CONDENSATION

              CONDENS     = 0.0
              GT(J,I)     = TGP
              CO2ICE(J,I) = 0.0D0

            ELSE

!             CASE 2B:  GROUND COOLED TO CO2 FROST POINT
!             AND SOME CO2 ICE REMAINS
!             COMPUTE AMOUNT OF ICE REMAINING

              CONDENS     = CONDENS*(TSAT-TGP)/(GT(J,I)-TGP)
              CO2ICE(J,I) = CONDENS*NDT
              GT(J,I)     = TSAT
            END IF
          END IF

        END IF

!       OUTPUT VALUES BASED ON NET CONDENSATION

        DMADT(J,I)  = CONDENS

      RETURN
      END
