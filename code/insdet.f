      SUBROUTINE INSDET(rsetsw)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

C
C  PURPOSE:
C      THIS INITIALIZATION SECTION SUBROUTINE CALLS SUBROUTINE SDET,
C      AND IN THE CASE OF A COLD START, CALCULATES THE POLAR
C      CAP AND ICE CLOUD BOUNDARIES.
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
C
C  OUTPUT PARAMETERS
C      JBPS, JBPN, JBCS, JBCN,  & ?
C  CALLED BY
C
C  SUBROUTINES CALLED:
C
      use grid_h
      use defines_h
      use standard_h

      implicit none

C######################################################################

      INTEGER CLKSW, RSETSW

      integer :: RESET = 1
      integer :: OFF   = 0

      integer lday
      namelist / insdetnl / LDAY

!     implicit none
   
      integer :: i, j, inu, lyr, immpn, ibound, ibp, jmdel, immcN
      integer :: immps, immcs
      real*8  :: dellat, xlatbd

C#=====================================================================

C New C grid
      DO 11 J=1,JM-1
        SINL(J)=SIN(LAT(J))
        COSL(J)=COS(LAT(J))
   11 CONTINUE

      INU = 5

      READ(INU,insdetnl)
      CLKSW = 0
      LYR   = 0

   31 CONTINUE

      IF(RSETSW.NE.RESET) GO TO 14
      SDEDY = LDAY
      SDEYR = LYR

   14 CONTINUE

      DCLK  = .FALSE.
      IGROW = 0
      CALL SDET

C     write out the GCM log file

      CALL GCMLOG(CLKSW,RSETSW,LDAY,LYR)

      IF(TAU.GT.0.) GO TO 2000
C      DELLAT = 180.0/FLOAT(JM-1)
C C grid
      DELLAT = 180.0/FLOAT(JM)

C     NORTH POLAR CAP

      IMMPN = IMAXPN-1

      DO 1110 I=1,IMMPN
        IBOUND = I
        IBP    = I+1
        IF(VPOUT.LT.PNLS(IBP)) GO TO 1120
 1110 CONTINUE

      JBPN = JM+1
      GO TO 1130

 1120 CONTINUE

      XLATBD = PNLAT(IBOUND)+(PNLAT(IBP)-PNLAT(IBOUND))*(VPOUT-
     *           PNLS(IBOUND))/(PNLS(IBP)-PNLS(IBOUND))
      JMDEL  = (90.0-XLATBD)/DELLAT
C C grid: if JBPN = JMDEL this point isn't used in C grid
      JBPN   = JM-JMDEL

      IF(XLATBD.GT.(90.-DELLAT/4.)) THEN
        JBPN = JM + 1
      ENDIF

C     NORTH CLOUD

 1130 CONTINUE

      IMMCN = IMAXCN-1

      DO 1140 I=1,IMMCN
        IBOUND = I
        IBP    = I+1
        IF(VPOUT.LT.CNLS(IBP)) GO TO 1150
 1140 CONTINUE

 1150 CONTINUE

      XLATBD = CNLAT(IBOUND)+(CNLAT(IBP)-CNLAT(IBOUND))*(VPOUT-
     *           CNLS(IBOUND))/(CNLS(IBP)-CNLS(IBOUND))
      JMDEL  = (90.-XLATBD)/DELLAT
C C grid: if JBPN = JMDEL this point isn't used in C grid
      JBCN   = JM-JMDEL

      IF(XLATBD.GT.(90.-DELLAT/4.)) THEN
        JBCN = JM + 1
      ENDIF

C     SOUTH POLAR CAP

      IMMPS = IMAXPS-1
      IF(VPOUT.LT.PSLS(1)) GO TO 1165

      DO 1160 I=1,IMMPS
        IBOUND = I
        IBP    = I+1
        IF(VPOUT.LT.PSLS(IBP)) GO TO 1170
 1160 CONTINUE

 1165 CONTINUE

      JBPS = 0
      GO TO 1180

 1170 CONTINUE

      XLATBD = PSLAT(IBOUND)+(PSLAT(IBP)-PSLAT(IBOUND))*(VPOUT-
     *           PSLS(IBOUND))/(PSLS(IBP)-PSLS(IBOUND))
      IF(XLATBD.GE.90.) GO TO 1165
      JMDEL  = (90.-XLATBD)/DELLAT
C      JBPS   = JMDEL+1
C C PI grid does not start at pole
      JBPS   = JMDEL

C     SOUTH CLOUD

 1180 CONTINUE

      IMMCS = IMAXCS-1

      DO 1190 I=1,IMMCS
        IBOUND = I
        IBP    = I+1
        IF(VPOUT.LT.CSLS(IBP)) GO TO 1192
 1190 CONTINUE

 1192 CONTINUE

      XLATBD = CSLAT(IBOUND)+(CSLAT(IBP)-CSLAT(IBOUND))*(VPOUT-
     *           CSLS(IBOUND))/(CSLS(IBP)-CSLS(IBOUND))
      JMDEL  = (90.-XLATBD)/DELLAT

      IF(XLATBD.GE.90.) THEN
        JMDEL=-1
      ENDIF

C      JBCS = JMDEL+1
C C PI grid does not start at pole
      JBCS = JMDEL

 2000 CONTINUE

      IF(CLKSW.NE.OFF) DCLK = .TRUE.

      RETURN
      END
