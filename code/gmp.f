      SUBROUTINE GMP

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

C
C  PURPOSE:
C      GMP CALCULATES THE POLAR CAP AND ICE CLOUD BOUNDARIES AND
C      MAKES SURFACE PRESSURE ADJUSTMENTS TO INSURE CONSERVATION OF
C      MASS.
C  MODIFIED FROM ORIGINAL BY
C      STEVE POHORSKY    INFORMATICS     TASK 605    SEP 82
C  FOR
C      JIM POLLACK
C  ENVIRONMENT
C      Cray-2            UNICOS 3.0      FORTRAN
C  REVISION HISTORY
C      JR SCHAEFFER     TASK 904           DEC 86
C      MODIFICATIONS TO INSURE TREATMENT OF THE SOUTH CAP AND CLOUDS
C      IS THE SAME AS FOR THE NORTH CAP AND CLOUDS.
C                                              (SEE NOTE OF 10/24/86).
C
C  INPUT PARAMETERS
C      VPOUT          - THE CURRENT SEASONAL DATE INDEX, LS,
C                             COMPUTED IN SDET ONCE PER DAY.
C      P(J,I) ARRAY   - THE CURRENT 'PI' VALUE (SURFACE PRESSURE
C                             MINUS TROPOPAUSE PRESSURE) AT THE 'PI'
C                             POINTS.
C      GT(J,I) ARRAY  - GROUND TEMPERATURE OR (-1) TIMES CO2 ICE
C                             MASS AT THE 'PI' POINTS.
C
C  OUTPUT PARAMETERS
C      JBPS           - APPROXIMATE LATITUDE INDEX OF EDGE OF SOUTHERN
C                       POLAR CAP.
C      JBPN           - APPROXIMATE LATITUDE INDEX OF EDGE OF NORTHERN
C                       POLAR CAP.
C      JBCS           - APPROXIMATE LATITUDE INDEX OF EDGE OF SOUTHERN
C                       POLAR WATER ICE CLOUD.
C      JBCN           - APPROXIMATE LATITUDE INDEX OF EDGE OF NORTHERN
C                       POLAR WATER ICE CLOUD.
C      P(J,I) ARRAY   - AS IN INPUT, BUT ADJUSTED TO PROVIDE
C                       CONSERVATION OF MASS.
C      GT(J,I) ARRAY  - AS IN INPUT, BUT ADJUSTED IN POLAR REGIONS
C                       ACCORDING TO NEW POLAR CAP BOUNDARIES AND
C                       PRESSURE VALUES.
C  CALLED BY
C      MAIN
C
      use grid_h
      use defines_h
      use constants_h, only: GRAV
      use standard_h

      implicit none

C######################################################################

      real*8  :: ZM(L_JSIZE), ZMG(L_JSIZE)

      integer :: i, j, immpn, ibound, ibp, jmdel, immcn, immps, immcs
      real*8  :: xls, ratio, zmtot, wtm, zmmg, xlatbd, dellat

C#=====================================================================

C     ALLOWING FOR CO2 FROST ON GROUND WHEN SDEDY IS ADVANCED BY SUNSTP.
C     WE CORRECT CO2 AMOUNT BY THIS FACTOR FIRST.  LETTER G DENOTES
C     GROUND.  GMP MUST BE CALLED ONCE A DAY.

c     DELLAT = 180./FLOAT(JM-1)
C  C-grid core change
      DELLAT = 180./FLOAT(JM)

C     NORTH POLAR CAP

      IMMPN = IMAXPN-1

      DO 1110 I=1,IMMPN
        IBOUND = I
        IBP    = I+1
        IF(VPOUT.LT.PNLS(IBP)) GO TO 1120
 1110 CONTINUE

      JBPN = JM+1
      GO TO 1130

 1120 XLATBD = PNLAT(IBOUND)+(PNLAT(IBP)-PNLAT(IBOUND))*(VPOUT-
     *         PNLS(IBOUND))/(PNLS(IBP)-PNLS(IBOUND))
      JMDEL  = (90.-XLATBD)/DELLAT

C  c-grid:  if JBPN = JMDEL this point isn't used in the c-grid
      JBPN   = JM-JMDEL

      IF(XLATBD.GT.(90.-DELLAT/4.0)) THEN
        JBPN = JM+1
      ENDIF

C     NORTH CLOUD

 1130 IMMCN = IMAXCN-1

      DO 1140 I=1,IMMCN
        IBOUND = I
        IBP    = I+1
        IF(VPOUT.LT.CNLS(IBP)) GO TO 1150
 1140 CONTINUE

 1150 XLATBD = CNLAT(IBOUND)+(CNLAT(IBP)-CNLAT(IBOUND))*(VPOUT-
     *         CNLS(IBOUND))/(CNLS(IBP)-CNLS(IBOUND))
      JMDEL  = (90.0-XLATBD)/DELLAT
C  c-grid:  if JBCN = JMDEL this point isn't used in the c-grid
      JBCN   = JM - JMDEL

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

 1165 JBPS = 0
      GO TO 1180

 1170 XLATBD = PSLAT(IBOUND)+(PSLAT(IBP)-PSLAT(IBOUND))*(VPOUT-
     *         PSLS(IBOUND))/(PSLS(IBP)-PSLS(IBOUND))
      JMDEL  = (90.0-XLATBD)/DELLAT

C     JBPS   = JMDEL+1
C  C-grid PI grid does not start at the pole
      JBPS   = JMDEL

      IF(XLATBD.GT.(90.0-DELLAT/4.0)) THEN
        JBPS = 0
      ENDIF

C     SOUTH CLOUD

 1180 IMMCS = IMAXCS-1

      DO 1190 I=1,IMMCS
        IBOUND = I
        IBP    = I+1
        IF(VPOUT.LT.CSLS(IBP)) GO TO 1200
 1190 CONTINUE

 1200 XLATBD = CSLAT(IBOUND)+(CSLAT(IBP)-CSLAT(IBOUND))*(VPOUT-
     *         CSLS(IBOUND))/(CSLS(IBP)-CSLS(IBOUND))
      JMDEL  = (90.0-XLATBD)/DELLAT
c     JBCS   = JMDEL+1
C  C-grid PI grid does not start at the pole
      JBCS   = JMDEL

      IF(XLATBD.GT.(90.0-DELLAT/4.0)) THEN
        JBCS = 0
      ENDIF

C     MAKE SURFACE PRESSURE ADJUSTMENTS TO INSURE CONSERVATION OF MASS.

      FIM = IM

c     DO 1360 J=1,JM
C  c-grid change
      DO 1360 J=1,JM-1
        ZM(J) = 0.0

        DO 1350 I=1,IM
          ZM(J) = ZM(J)+P(J,I)
 1350   CONTINUE

        ZM(J) = ZM(J)/FIM
 1360 CONTINUE

c     DO 1460 J=1,JM
C  c-grid change
      DO 1460 J=1,JM-1
        ZMG(J) = 0.0

        DO 1450 I=1,IM
          ZMG(J) = ZMG(J) + CO2ICE(J,I)
 1450   CONTINUE

        ZMG(J) = ZMG(J)/FIM*0.01*GRAV
 1460 CONTINUE

      ZMMG = 0.0
      WTM  = 0.0
      ZMM  = 0.0

c     DO 1465 J=1,JM
C  c-grid change
      DO 1465 J=1,JM-1
        WTM = WTM+ABS(DXYP(J))
        ZMM = ZMM+ZM(J)*ABS(DXYP(J))
 1465 CONTINUE

      ZMM = ZMM/WTM+PTROP

c     DO 1470 J=1,JM
C  c-grid change
      DO 1470 J=1,JM-1
        ZMMG = ZMMG+ZMG(J)*ABS(DXYP(J))
 1470 CONTINUE

      ZMMG  = ZMMG/WTM
      ZMTOT = ZMM+ZMMG

C     Add pressure of CO2 ice on the ground.  Set in BKDATA.

      RATIO = (PSF+PSFGND-ZMMG)/ZMM

C     MAKE ADJUSTMENTS TO P.

      DO 3015 I=1,IM
c       DO 3010 J=1,JM
C  c-grid change
        DO 3010 J=1,JM-1

          P(J,I) = RATIO*(P(J,I)+PTROP)-PTROP

          IF(P(J,I).LE.0.) THEN
            WRITE(MTP,3002) J,I,P(J,I)
 3002       FORMAT(5X,'TROUBLE IN GMP.  P.LE.0.',5X,'J=',I4,3X,
     *                'I=', I4, 3X, 'P=', E12.5 )
            STOP
          ENDIF

 3010   CONTINUE
 3015 CONTINUE

      ZMM   = RATIO*ZMM
      ZMTOT = ZMM+ZMMG

C     Scale the TAUDST values to take into account varying dust
C     optical depth.  The next sol's value of dust optical depth is
C     found by calling nextls to get the Ls of the new day, and
C     then calling interpdust to interpolate the dust values (at
C     integer values of Ls, i.e., at Ls=0,1,2,3,4...) to the real
C     value of the next sol's Ls.  The value is returned in TAUTOT.

      if(ACTIVE_DUST .eqv. .TRUE.) then

        call nextls(igrow,igmax,sdedy,sunstp,anome,eccn,vinc,xls)
        call getvdust(XLS,NVDUST,VDUSTLS,TAUTOTJI)

      else

        call nextls(igrow,igmax,sdedy,sunstp,anome,eccn,vinc,xls)

C  VDUST method of opacity - read opacity maps if VDUST = .TRUE.,
C  else use old GCM "TAUTOT" method  /7/19/01

        if(VDUST) then
          call getvdust(XLS,NVDUST,VDUSTLS,TAUTOTJI)
        else
          call interpdust(xls,dustod,TAUTOT)
          do J=1,L_J
            do I=1,L_I
              TAUTOTJI(J,I) = TAUTOT*TAUPAT(J,I)
            end do
          end do
        endif

        call dustprofile
      end if

C     Polar retreating section of the code.
 
      RETURN
      END
