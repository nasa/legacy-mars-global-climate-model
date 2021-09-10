      SUBROUTINE SDET

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

C  PURPOSE:
C      THIS SUBROUTINE CALCULATES QUANTITIES RELATED TO THE
C      POSITION OF MARS IN ITS ORBIT.
C
C  MODIFIED FROM ORIGINAL BY
C      STEVE POHORSKY    INFORMATICS     TASK 605    JUL 82
C  FOR
C      JIM POLLACK
C  ENVIRONMENT
C      Cray-2            UNICOS 3.0      FORTRAN
C  REVISION HISTORY
C
C  INPUT PARAMETERS
C      SDEYR, SDEDY, & ?
C  OUTPUT PARAMETERS
C      COSD, SIND, RDIST, SDEYR, SDEDY, & ?
C  CALLED BY
C      MAIN,INSDET
C
      use grid_h
      use defines_h
      use constants_h, only: PI
      use standard_h

      implicit none

C######################################################################

      INTEGER DY

      integer :: maxday
      real*8  :: quad, si, co, x, arsin, arcos, raddeg, hellon, eta
      real*8  :: eta2, eta1, beta, anomtp, cose, sine, v1, v2, anomt
      real*8  :: dec

      QUAD(SI,CO) = PI-SIGN(PI-CO,SI)
      ARSIN(X)    = ASIN(X)
      ARCOS(X)    = ACOS(X)

C#=====================================================================

C     CALCULATES ORBITAL INFORMATION:  DISTANCE TO SUN, SOLAR
C     DECLINATION, HELIOCENTRIC LONGITUDE, AND TRUE ANOMOLY FROM
C     SPRING EQUINOX IN NORTHERN HEMISPHERE

      IGROW = IGROW+1

      IF(IGROW.GT.IGMAX) THEN
        DCLK =.TRUE.
      ENDIF

      IF(IGROW.EQ.1) GOTO 210

      IF(IGROW.LE.IGMAX) THEN
        SDEDY = SDEDY+1
      ENDIF

      IF(IGROW.GT.IGMAX) THEN
        SDEDY = SDEDY+SUNSTP
      ENDIF

  210 CONTINUE

      DY     = SDEDY+1
      MAXDAY = DAYPYR+1.0E-2

      IF(SDEDY.LE.MAXDAY) GOTO 211

      SDEDY = SDEDY-MAXDAY-1
      SDEYR = SDEYR+1

  211 CONTINUE

      COSE   = COS(ANOME(DY))
      SINE   = SIN(ANOME(DY))
      RSDIST = (AORB*(1.0- ECCN*COSE))**2
      V1     = ARCOS((COSE-ECCN)/(1.-ECCN*COSE))
      V2     = ARSIN(ESQ*SINE/(1.-ECCN*COSE))
      ANOMT  = QUAD(V2,V1)
      SIND   = SIN(DECMAX)*COS(ANOMT-VINC)
      DEC    = ARSIN(SIND)
      COSD   = COS(DEC)
      ANOMTP = ANOMT-VINC+PI*0.5

      IF(ANOMTP.LT.0.) THEN
        ANOMTP = ANOMTP+2.0*PI
      ENDIF

      BETA   = ARSIN(SIN(ORBINC)*SIN(ANOMT-ANMTSP))
      ETA1   = ARCOS(COS(ANOMT-ANMTSP)/COS(BETA))
      ETA2   = ARSIN(COS(ORBINC)*SIN(ANOMT-ANMTSP)/COS(BETA))
      ETA    = QUAD(ETA2,ETA1)
      HELLON = ETA+ETAS

      IF(HELLON.GT.2.*PI) THEN
        HELLON = HELLON-2.0*PI
      ENDIF

      RADDEG = 180./PI
      VOUT   = ANOMT*RADDEG
      VPOUT  = ANOMTP*RADDEG
      ETAOUT = HELLON*RADDEG

      IF(VOUT.LT.0.) THEN
        VOUT = VOUT+360.0
      ENDIF

      IF(VPOUT.LT.0.) THEN
        VPOUT = VPOUT+360.0
      ENDIF

      IF(ETAOUT.LT.0.) THEN
        ETAOUT = ETAOUT+360.0
      ENDIF

      VOUT   = MOD(VOUT,360.0D0)
      VPOUT  = MOD(VPOUT,360.0D0)
      ETAOUT = MOD(ETAOUT,360.0D0)

      WRITE(MTP,50) VOUT, VPOUT, ETAOUT
      WRITE(MTP,51) SDEDY, SDEYR, IGROW
      WRITE(MTP,52) DEC, RSDIST

   50 FORMAT('TRUE ANAMOLY=',F15.7,3X,'LS=',F15.7,3X,
     *           'HELIOCENTRIC LONGITUDE =',F15.7)
   51 FORMAT(1X, 'DAY=',I5,5X, 'YEAR=',I3,5X, 'IGROW=',I5)
   52 FORMAT(5X, 'DEC=',E13.6,5X, 'RSDIST=',E12.6)

      RETURN
      END
