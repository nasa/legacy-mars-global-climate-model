      SUBROUTINE EMISS

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

C
C  PURPOSE
C      EMISS COMPUTES THE CO2 EMISSIVITY INTERPOLATION TABLE
C      AND VARIOUS OTHER VALUES USED IN SUBSEQUENT ROUTINES,
C      PARTICULARLY COMP3.
C  AUTHOR
C      STEVE POHORSKY    INFORMATICS     TASK 605    JAN 82
C
C  FOR
C      JIM POLLACK      (SEE NOTES OF 9-16-81,3-30-82,5-17-82)
C
C  ENVIRONMENT
C      Cray-2            UNICOS 3.0      FORTRAN
C
C  REVISION HISTORY
C     6-82  SP   EMISSIVITY CALCULATIONS ADDED AS IN NOTE OF 3-30-82.
C     7-82  SP   NAME OF THIS SUBROUTINE CHANGED FROM EQWID TO EMISS.
C     7-82  SP   EMISSIVITY CALCULATIONS ADDED AS IN NOTE OF 5-17-82.
C     2-83  SP   DUST AND ICE EMISSIVITIES PUT IN EMISSES BLOCK DATA.
C      OCT 84  JB WHITE  TASK 904  FINE-MESH REVISIONS INCORPORATED
C                                  (SEE NOTE OF SEPTEMBER 10, 1984).
C      JR SCHAEFFER        TASK 904          JAN 87
C      TAYLOR SERIES APPROXIMATION TO LOG(PRESSURE) INCLUDED
C      USED IN IN15IR                   (SEE NOTE OF 1/13/87)
C      J. SCHAEFFER        TASK 904          JAN 87
C      INCORPORATED VARIABLE DUST  ( SEE NOTE OF 1/19/87 ).
C  INPUT PARAMETERS
C      NLAY            - THE NUMBER OF LAYERS (PARALLEL TO THE SURFACE)
C                        THAT THE ATMOSPHERE IS DIVIDED INTO IN THE
C                        3-DIMENSIONAL GRID REPRESENTATION.
C      DSIG(K) ARRAY   - THE DIFFERENCE IN SIGMA (VERTICAL COORDINATE)
C                        BETWEEN THE TOP AND BOTTOM OF THE KTH LAYER
C                        OF THE ATMOSPHERE. ( K = 1 TO NLAY. )
C      PSL             - A CONSTANT EQUAL TO THE PRESSURE AT A MEAN
C                        SURFACE ELEVATION (CALLED SEA LEVEL) FOR MARS.
C      DT              - ROUGH APPROXIMATION OF NUMBER OF SECONDS IN
C                        ONE TIME STEP.
C  OUTPUT PARAMETERS
C      SIGMA(K) ARRAY  - THE SIGMA (VERTICAL COORDINATE) FOR EACH LEVEL
C                        OF THE ATMOSPHERE. ( K = 3 TO NLEVELS. )
C      AADJ(K) AND  -  THE FIRST TWO COEFFICIENTS IN THE TAYLOR'S SERIES
C      BADJ(K) ARRAYS  FOR THE DIFFERENCE IN LOG OF PRESSURE BETWEEN
C                      ADJACENT LAYERS. AADJ AND BADJ ARE, RESPECTIVELY,
C                      THE FUNCTION AND ITS DERIVATIVE EVALUATED AT PSL.
C                      ( K = 3, NLEVSM1. )
C      RICHC1          - A COEFFICIENT USED IN SUBROUTINE RICHNUM.
C      RICHC2          - A COEFFICIENT USED IN SUBROUTINE RICHNUM.
C      DAY             - THE NUMBER OF (EARTH) SECONDS IN A MARTIAN DAY
C                        = 88775.
C      DT              - EXACT NUMBER OF SECONDS IN ONE TIME STEP. THE
C                        INPUT DT VALUE IS DECREASED TO AN EVEN FRACTION
C                        OF DAY.
C      SINA(I) AND     - SINE & COSINE OF LONGITUDES OF U,V GRID POINTS.
C      COSA(I) ARRAYS
C      ANOME(N) ARRAY  - ECCENTRIC ANOMOLIES AT MIDPOINT OF EACH DAY
C      FMEM(K,L,2,MP) ARRAY  EMISSIVITY BETWEEN LAYER K AND FINE-MESH
C                            LEVEL L FOR DOWN AND UP FLUXES AND WITH
C                            RESPECT TO SURFACE PRESSURE DEFINED BY MP;
C                            PASSED THROUGH COMMON /FMESH0/ TO IN15IR.
C
C  CALLED BY
C      INPUT
C
      use grid_h
      use defines_h
      use constants_h, only: PI, GRAV, CP, KAPA
      use standard_h
      use comp3cmn_h

      implicit none

C----------------------------------------------------------------------C

      real*8  :: OUTCOEF(3), QRDST(L_NPDST)
      integer :: i, k, l, n, ntd1, ntd, ntwo
      real*8  :: degrad, e1, e0, anomm, time, alf, dvar

C======================================================================C

C     SET UP SOME CONSTANTS FOR USE IN OTHER ROUTINES.

C     Number of vertical levels in the model atmosphere

      NLEVELS = 2 * NLAY + 3
      NLEVSM1 = NLEVELS - 1
      NLEVSM2 = NLEVELS - 2
      NLEVSM3 = NLEVELS - 3
      NLEVSM4 = NLEVELS - 4
      NLAYM1  = NLAY-1
      TWOPI   = 2.0*PI
      PLAPR   = FUDG*GRAV/CP
      PKAPA   = KAPA*FUDG
      FIM     = IM

C     SIGMA at layer boundaries

      SIGMA(3) = 0.0
      DYSIG(1) = 0.0
      DO  50  K = 1, NLAY
        L = 2*K+3
        SIGMA(L)   = SIGMA(L-2)+DSIG(K)
        DYSIG(K+1) = SIGMA(L)
   50 CONTINUE

C     SIGMA at the layer mid-points.

      DO  55  K = 4, NLEVSM1, 2
        SIGMA(K) = 0.5*(SIGMA(K+1)+SIGMA(K-1))
   55 CONTINUE

C     COMPUTE THE FIRST TWO COEFFICIENTS IN THE TAYLOR'S SERIES FOR
C     THE DIFFERENCE IN LOG OF PRESSURE BETWEEN ADJACENT LAYERS.

      PL(2) = PTROP/2.0

      DO 1100 K = 3, NLEVELS
        PL(K) = PTROP + SIGMA(K) * ( PSL - PTROP )
 1100 CONTINUE

C     Evaluate the function and it's derivative at 'sea level' pressure

      DO 1200 K = 3, NLEVSM1
        AADJ(K) = LOG(PL(K+1)/PL(K))
        BADJ(K) = SIGMA(K+1)/PL(K+1)-SIGMA(K)/PL(K)
 1200 CONTINUE

      PLOGADJ(2)=LOG(2.0)

C     COMPUTE SINE & COSINE OF LONGITUDES OF U,V GRID POINTS.

      DO 2100 I=1,IM
        ALF     = TWOPI*(FLOAT(I)+0.5)/FLOAT(IM)
        SINA(I) = SIN(ALF)
        COSA(I) = COS(ALF)
 2100 CONTINUE

C     TIME CONVERSIONS

C     PREV is the number of Earth seconds in a Martian year.

      PREV   = PREV*8.64E4
      DAY    = TWOPI/ROT
      PROTCR = DAY

C     DAY is the number of (Earth) seconds in a Martian day.

      DAY    = DAY*PREV/(PREV-DAY)
      SQRDY  = SQRT(4.0*PI/DAY)
      TSTART = TSTART*DAY

C     DAYPYR is the number of Martian days in a Martian year.

      DAYPYR = PREV/DAY

C     MAKING DT INTO A SIMPLE FRACTION OF A MARTIAN DAY
C     DVAR is an even fraction of DAY.
C     DT is an even fraction of DVAR.

      DVAR = DAY*MIN(TAUO2,TAUH)/ROTPER
      DVAR = DVAR/NC3
      NTWO = DVAR/DT+1.0
      DT   = DVAR/NTWO

      WRITE(MTP,2110) DT,DAY,PROTCR
 2110 FORMAT(10X,'ADJUSTED DT=',E12.6,5X,'DAY=',E12.6,2X,'SEC',5X,
     *           'COR.DAY=',E12.6)

      CNIN = DT/(NSMTH*DAY)

C     CALCULATION OF ECCENTRIC ANOMOLIES AT THE MIDPOINT OF EACH DAY

      ESQ = SQRT(1.0-ECCN**2)
      NTD = INT(IGMAX+(TAURUN-IGMAX+1)*SUNSTP)+1

      IF(NTD.LT.IGMAX) NTD=IGMAX

      NTD1= NTD +1

      DO 3085 N=1,NTD1
        TIME  = TSTART+(FLOAT(N)-0.5)*DAY
        ANOMM = TWOPI*TIME/PREV
        E0    = ANOMM

        DO 3081 I=1,100
          E1 = ANOMM+ECCN*SIN(E0)
          IF(ABS(1.-E0/E1).LT.1.E-4) GO TO 3082
          E0 = E1
 3081   CONTINUE

 3082   CONTINUE
        ANOME(N) = E1
 3085 CONTINUE

      DO 4000 L=1,L_LAYERS
        RAYK(L) = 0.0
 4000 CONTINUE

      DO 4010 L=1,LRAY
        RAYK(L) = ((PL(K1)/PL(2*L+2))**ALFRAY)/(TREFR*DAY)
 4010 CONTINUE

C     CONVERSION OF ORBITAL PARAMETERS FROM DEGREES TO RADIANS

      DEGRAD = PI/180.0

      DECMAX = DECMAX*DEGRAD
      VINC   = VINC*DEGRAD
      ETPER  = ETPER*DEGRAD
      ETAS   = ETAS*DEGRAD
      ORBINC = ORBINC*DEGRAD
      ANMTSP = ANMTSP*DEGRAD

C     CALCULATION OF DUST MIXING RATIO AND CUMULATIVE DUST OPTICAL
C     DEPTH AT A SEQUENCE OF REFERENCE PRESSURE LEVELS.

C Bob's updates 9/17/99 
C Reference the dust optical depth to RPTAU, and modify the way
C the dust mixing ratio is calculated to more accurately reflect the
C pressure-optical depth relationship.

      CALL DUSTPROFILE

      RETURN
      END
