       SUBROUTINE DYCORE ( IM,JM,JDIM,LM,SIG,PTOP,KM,DT,
     .                    OMEGA, CP, RGAS, AE, P00,
     .                    PHS  ,PKH,
     .                    PIB  ,UOB  ,VOB  ,POB  ,QOB,
     .                    PIM  ,UOM  ,VOM  ,POM  ,QOM,
     .                    PII,  UOI,  VOI,  POI,  QOI,
     .                    OMG,  VOR,  DIAG, JY1)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

!  No pointers version: 06-26-09
!  Mel's dycore
C **********************************************************************
C  PURPOSE
C     Update Time-Tendencies of Prognostic Fields 
C     due to Hydrodynamical Processes
C
C  INPUT ARGUMENT DESCRIPTION
C
C     IM .... Number of Grid Intervals in Zonal      Direction
C     JM .... Number of Grid Intervals in Meridional Direction
C     JDIM .. Meridional (Second) Dimension of Input  Fields 
C     LM .... Number of Vertical Levels
C     SIG ... (LM+1): Sigma at Interfaces. SIG(1)=0; SIG(LM+1)=1
C     PTOP .. Model Top Pressure at Sigma = 0.0 
C     KM .... Number of Scalars, Including H2O, but not Theta.
C     DT .... Time-Step fron n-1 to n+1 (Seconds)
C
C     OMEGA.. Rotation rate (rad/sec)
C     CP..... Specific heat at constant pressure (J/(kg K))
C     RGAS... Gas contant (J/(kg K))
C     AE..... 'Earth' radius (meters)
C
C     PHS ... (IM,JDIM): Surface Geopotential (m * m/sec**2)
C     PKH ... (IM,JDIM,LM+1): (P/P00)**KAPPA (Not used if PTOP=0)
C
C     PIB ... (IM,JDIM):    Mass (Psurf-Ptop) mb  at Current  Time-Level
C     UOB ... (IM,JDIM,LM): Zonal      Wind   m/s at Current  Time-Level
C     VOB ... (IM,JDIM,LM): Meridional Wind   m/s at Current  Time-Level
C     POB ... (IM,JDIM,LM): Potential Temp.   K   at Current  Time-Level
C     QOB ... (IM,JDIM,LM,KM): Scalar Fields      at Current  Time-Level
C
C     PIM ... (IM,JDIM):    Mass (Psurf-Ptop)     at Previous Time-Level
C     UOM ... (IM,JDIM,LM): Zonal      Wind       at Previous Time-Level
C     VOM ... (IM,JDIM,LM): Meridional Wind       at Previous Time-Level
C     POM ... (IM,JDIM,LM): Potential Temperature at Previous Time-Level
C     QOM ... (IM,JDIM,LM,KM): Scalar Fields      at Previous Time-Level
C
C  OUTPUT ARGUMENT DESCRIPTION-- Tendencies are in per second.
c
C     PII ... (IM,JDIM):       Updated Surface Pressure   Time-Tendency
C     UOI ... (IM,JDIM,LM):    Updated Zonal Wind         Time-Tendency
C     VOI ... (IM,JDIM,LM):    Updated Meridional Wind    Time-Tendency
C     POI ... (IM,JDIM,LM):    Updated PI-Weighted Theta  Time-Tendency
C     QOI ... (IM,JDIM,LM,KM): Updated PI-Weighted Scalar Time-Tendency
C
C     OMG ... (IM,JDIM,LM): Omega Diagnostic (mb/sec)
C     VOR ... (IM,JDIM,LM): Vorticity Diagnostic  (1/sec)
C     DIAG .. Logical (On/Off) Flag for Diagnostics
C
C  NOTES:
C     (1) JDIM is to be used for DIMENSION Purposes ONLY
C         Vectorization will be performed over IM*JM.
C     (2) The Vertical Layers are numbered from TOP(1) to BOTTOM(LM).
C     (3) All Time-Tendencies are INCREMENTED (bumped).
C         The Momentum Time-Tendencies ARE NOT mass-weighted.
C         The Potential Temperature and Scalar Time-Tendencies ARE 
C         mass-weighted (by PI).
C     (4) JM is 180 degrees divided by the meridional grid size.
C     (5) UXX(I,J) are located half a grid interval EAST of PXX(I,J).
C         VXX(I,J) are located half a grid interval SOUTH of PXX(I,J).
C     (6) If PTOP>0, the PKH MUST be defined. 
C     (7) The previous time level fields (PIM,UOM,etc) are used for the
C         economical explicit calculation done in conjunction with
C         leap-frog steps. If you are not doing leap-frog or do not
C         wish to have economical explicit tendencies, pass the current
C         time-level fields twice (i.e., in PIB,UOB,etc and again in
C         PIM,UOM,etc.). 
C
C  SPACE REQUIREMENTS:
C     (1) Takes IM*JM*19+4*JM+2*IM words from the heap for STATIC storag
C         these are kept throughout the run.
C     (2) Takes IM*JM*(LM+25) + 3*LM + 1 words from the heap for DYNAMIC
C         storage when PTOP=0; for PTOP
C         All of this storage is freed before returning.
C
C **********************************************************************


      IMPLICIT NONE

      integer :: ierror

      REAL ZERO, ONE, TWO, THREE, FOUR, FIVE, SIX
      REAL HALF, FOURTH, THIRD

      PARAMETER (ZERO=0.0)
      PARAMETER (ONE=1.0)
      PARAMETER (TWO=2.0)
      PARAMETER (THREE=3.0)
      PARAMETER (FOUR=4.0)
      PARAMETER (FIVE=5.0)
      PARAMETER (SIX =6.0)

      PARAMETER (FOURTH=ONE/FOUR)
      PARAMETER (HALF=ONE/TWO)
      PARAMETER (THIRD=ONE/THREE)

      INTEGER TMPVRS
      PARAMETER (TMPVRS=25)

      REAL ALPHA
      PARAMETER (ALPHA=.2475)

C ARGUMENTS

      INTEGER IM
      INTEGER JM
      INTEGER JDIM
      INTEGER LM
      INTEGER KM

      LOGICAL DIAG

      REAL OMEGA
      REAL AE
      REAL P00
      REAL CP
      REAL RGAS

      REAL PTOP
      REAL DT

      REAL SIG(LM+1)

      REAL UOB(IM,JDIM,LM)
      REAL VOB(IM,JDIM,LM)
      REAL POB(IM,JDIM,LM)
      REAL QOB(IM,JDIM,LM,KM)
      REAL PIB(IM,JM   )

      REAL UOM(IM,JDIM,LM)
      REAL VOM(IM,JDIM,LM)
      REAL POM(IM,JDIM,LM)
      REAL QOM(IM,JDIM,LM,KM)
      REAL PIM(IM,JM   )

      REAL UOI(IM,JDIM,LM)
      REAL VOI(IM,JDIM,LM)
      REAL POI(IM,JDIM,LM)
      REAL QOI(IM,JDIM,LM,KM)
      REAL PII(IM,JM   )

c     Franck's
      REAL QOI2(IM,JDIM,LM,KM)  ! Debug
      REAL USBSAV(IM,JM,LM)
      REAL VSBSAV(IM,JM,LM)
      INTEGER JY1
c     Franck's

      REAL PHS(IM,JM     )
      REAL PKH(IM,JDIM,LM+1)

      REAL OMG(IM,JDIM,LM)
      REAL VOR(IM,JDIM,LM)

C DYNAMIC LOCALS

C   THESE ARE DYNAMIC (RESET EACH TIME) LOCALS FOR NOW

      REAL DSG(LM)
      REAL PRJ(LM+1)
      REAL PRH(LM)

C   VECTOR TEMPORARIES

      REAL PSD(IM,JM,LM)
      REAL PKL(IM,JM,LM)
      REAL VT1(IM,JM)
      REAL VT2(IM,JM)
      REAL VT3(IM,JM)
      REAL VT4(IM,JM)
      REAL BET(IM,JM)
      REAL GAM(IM,JM)
      REAL DEL(IM,JM)
      REAL ALF(IM,JM)
      REAL PBI(IM,JM)
      REAL PBJ(IM,JM)
      REAL BIP(IM,JM)
      REAL BJP(IM,JM)
      REAL ACH(IM,JM)
      REAL ZOB(IM,JM)
      REAL DDX(IM,JM)
      REAL DDY(IM,JM)
      REAL USB(IM,JM)
      REAL VSB(IM,JM)
      REAL DPX(IM,JM)
      REAL DPY(IM,JM)
      REAL LAM(IM,JM)
      REAL MUU(IM,JM)
      REAL EPS(IM,JM)
      REAL PIV(IM,JM)
      REAL PIK(IM,JM)

      REAL PHI(IM,JM,LM)

C  SCALAR TEMPORARIES

      REAL PSUMS
      REAL PSUMN

      REAL SUMSO
      REAL SUMNO

      REAL ST1
      REAL ST2

      LOGICAL LEAP

C STATIC LOCALS

      REAL, allocatable :: DXPIJ(:,:)
      REAL, allocatable :: DYPIJ(:,:)
      REAL, allocatable :: DXUIJ(:,:)
      REAL, allocatable :: DYUIJ(:,:)
      REAL, allocatable :: DXVIJ(:,:)
      REAL, allocatable :: DYVIJ(:,:)

      REAL, allocatable :: D2PIJ(:,:)
      REAL, allocatable :: D2UIJ(:,:)
      REAL, allocatable :: D2VIJ(:,:)
      REAL, allocatable :: D2ZIJ(:,:)

      REAL, allocatable :: DYVIN(:,:)
      REAL, allocatable :: DXUIN(:,:)
      REAL, allocatable :: D2PIN(:,:)
      REAL, allocatable :: D2ZIN(:,:)

      REAL, allocatable :: FFFIJ(:,:)
      REAL, allocatable :: SU(:,:)
      REAL, allocatable :: SV(:,:)

      INTEGER, allocatable :: IE(:)
      INTEGER, allocatable :: IW(:)

      SAVE DXPIJ
      SAVE DYPIJ
      SAVE DXUIJ
      SAVE DYUIJ
      SAVE DXVIJ
      SAVE DYVIJ

      SAVE D2PIJ
      SAVE D2UIJ
      SAVE D2VIJ
      SAVE D2ZIJ

      SAVE DYVIN
      SAVE DXUIN
      SAVE D2PIN
      SAVE D2ZIN

      SAVE FFFIJ
      SAVE SU
      SAVE SV
      SAVE IE
      SAVE IW

C  THESE DETERMINE THE REINITIALIZATION OF THE HORIZONTAL GRID

      INTEGER IM0, JM0
      DATA IM0/0/, JM0/0/
      REAL AE0
      DATA AE0/0./

      INTEGER I, L, K, LL

      REAL P2K
      REAL SDOT

CFPP$ NOVECTOR R

      IF( IM.NE.IM0 .OR. JM.NE.JM0 .OR. AE.NE.AE0 ) THEN

        allocate(DXPIJ(IM,JM),stat=ierror)
        write(6,'("   ====== ""No Pointers"" dycore.f ======")')
        if(ierror.ne.0) then
          write(6,'("Could not allocate DXPIJ array in avrx.")')
          stop
        end if

        allocate(DYPIJ(IM,JM),stat=ierror)
        if(ierror.ne.0) then
          write(6,'("Could not allocate DYPIJ array in avrx.")')
          stop
        end if

        allocate(DXUIJ(IM,JM),stat=ierror)
        if(ierror.ne.0) then
          write(6,'("Could not allocate DXUIJ array in avrx.")')
          stop
        end if

        allocate(DYUIJ(IM,JM),stat=ierror)
        if(ierror.ne.0) then
          write(6,'("Could not allocate DYUIJ array in avrx.")')
          stop
        end if

        allocate(DXVIJ(IM,JM),stat=ierror)
        if(ierror.ne.0) then
          write(6,'("Could not allocate DXVIJ array in avrx.")')
          stop
        end if

        allocate(DYVIJ(IM,JM),stat=ierror)
        if(ierror.ne.0) then
          write(6,'("Could not allocate DYVIJ array in avrx.")')
          stop
        end if

        allocate(D2PIJ(IM,JM),stat=ierror)
        if(ierror.ne.0) then
          write(6,'("Could not allocate D2PIJ array in avrx.")')
          stop
        end if

        allocate(D2UIJ(IM,JM),stat=ierror)
        if(ierror.ne.0) then
          write(6,'("Could not allocate D2UIJ array in avrx.")')
          stop
        end if

        allocate(D2VIJ(IM,JM),stat=ierror)
        if(ierror.ne.0) then
          write(6,'("Could not allocate D2VIJ array in avrx.")')
          stop
        end if

        allocate(D2ZIJ(IM,JM),stat=ierror)
        if(ierror.ne.0) then
          write(6,'("Could not allocate D2ZIJ array in avrx.")')
          stop
        end if

        allocate(DYVIN(IM,JM),stat=ierror)
        if(ierror.ne.0) then
          write(6,'("Could not allocate DYVIN array in avrx.")')
          stop
        end if

        allocate(DXUIN(IM,JM),stat=ierror)
        if(ierror.ne.0) then
          write(6,'("Could not allocate DXUIN array in avrx.")')
          stop
        end if

        allocate(D2PIN(IM,JM),stat=ierror)
        if(ierror.ne.0) then
          write(6,'("Could not allocate D2PIN array in avrx.")')
          stop
        end if

        allocate(D2ZIN(IM,JM),stat=ierror)
        if(ierror.ne.0) then
          write(6,'("Could not allocate D2ZIN array in avrx.")')
          stop
        end if

        allocate(FFFIJ(IM,JM),stat=ierror)
        if(ierror.ne.0) then
          write(6,'("Could not allocate FFFIJ array in avrx.")')
          stop
        end if

        allocate(SU(IM+2,JM),stat=ierror)
        if(ierror.ne.0) then
          write(6,'("Could not allocate SU array in avrx.")')
          stop
        end if

        allocate(SV(IM+2,JM),stat=ierror)
        if(ierror.ne.0) then
          write(6,'("Could not allocate SV array in avrx.")')
          stop
        end if

        allocate(IE(IM*JM),stat=ierror)
        if(ierror.ne.0) then
          write(6,'("Could not allocate IE array in avrx.")')
          stop
        end if

        allocate(IW(IM*JM),stat=ierror)
        if(ierror.ne.0) then
          write(6,'("Could not allocate IW array in avrx.")')
          stop
        end if

         CALL GRIDH(
     *              IM,JM, AE, OMEGA
     *,             DXPIJ, DYPIJ, DXUIJ, DYUIJ, DXVIJ, DYVIJ
     *,             D2PIJ, D2UIJ, D2VIJ, D2ZIJ
     *,             DYVIN, DXUIN, D2PIN, D2ZIN, FFFIJ, SU, SV
     *,             IE, IW
     *           )
         IM0 = IM
         JM0 = JM
         AE0 = AE
      ENDIF

      CALL GRIDV(
     *           LM, RGAS/CP, SIG, DSG, PRJ, PRH, P00
     *          )

CFPP$ EXPAND(P2K,HADVECT)

C   L-INDEPENDENT QUANTITIES

      LEAP = DT .GT. ZERO


CMIC$ DO ALL AUTOSCOPE VECTOR
      DO I =1,IM*(JM-1)
       PBI(I,1) = (PIB(I,1) + PIB(IE(I),1)) * HALF
       BIP(I,1) = ONE / PBI(I,1)
       PIV(I,1) = ONE / PIB(I,1)
       PIK(I,1) = P2K(PIB(I,1))
       VT2(I,1) = D2PIJ(I,1) * PIB(I,1)
      ENDDO

CMIC$ DO ALL AUTOSCOPE VECTOR
      DO I=1,IM*(JM-1)
       VT1(I,1) = (VT2(I,1) + VT2(IE(I),1)) * HALF
      ENDDO

CMIC$ DO ALL AUTOSCOPE VECTOR
      DO I =1,IM*(JM-2)
       PBJ(I,2) = (PIB(I,2) + PIB(I,1)) * HALF
       BJP(I,2) = ONE / PBJ(I,2 )
       ACH(I,2) = D2ZIJ(I,2) / (HALF* (VT1(I,2) + VT1(I,1)))
      ENDDO

      PSUMS = (HALF/FLOAT(IM)) * SDOT(IM,PIB(1,   1),1,D2PIJ(1,   1),1)
      PSUMN = (HALF/FLOAT(IM)) * SDOT(IM,PIB(1,JM-1),1,D2PIJ(1,JM-1),1)

      DO I =1,IM
       ACH(I, 1) = D2ZIJ(I, 1) * (ONE / PSUMS)
       ACH(I,JM) = D2ZIJ(I,JM) * (ONE / PSUMN)
      ENDDO

C***************** BEGIN FIRST L LOOP ***********************
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

CMIC$  DO ALL AUTOSCOPE
CMIC$* PRIVATE(USB,VSB,VT1,VT2,VT3,ZOB,ALF,BET,GAM,DEL,LAM,MUU,EPS,ST1)

      DO 1000 L=1,LM

C  P TO THE KAPPA 

      IF(PTOP.GT.ZERO) THEN
       DO I=1,IM*(JM-1)
        ST1 = PIB(I,1)*SIG(L  )+PTOP
        ST2 = PIB(I,1)*SIG(L+1)+PTOP
        PKL(I,1,L) = (ST2*PKH(I,1,L+1)-ST1*PKH(I,1,L))
     *             / ( (ST2-ST1)*((RGAS/CP)+ONE) )
       ENDDO
      ENDIF

C   MASS FLUXES

      DO I =1,IM*(JM-1)
       USB(I,1) = DYUIJ(I,1) * PBI(I,1) * UOB(I,1,L)
       USBSAV(I,1,L) = USB(I,1)
      ENDDO

      DO I =1,IM*(JM-2)
       VSB(I,2) = DXVIJ(I,2) * PBJ(I,2) * VOB(I,2,L)
       VSBSAV(I,2,L) = VSB(I,2)
      ENDDO

      DO I =1,IM
       VSB(I,1 ) = ZERO
       VSB(I,JM) = ZERO
       VSBSAV(I,1,L ) = ZERO
       VSBSAV(I,JM,L) = ZERO
      ENDDO

C  HORIZONTAL ADVECTION OF TRACERS

cc      IF(KM.GT.0) THEN
cc        DO K =1,KM
cc          CALL HADVECT(USB(1,1),VSB(1,1),QOB(1,1,L,K)
cc     *,                IM,JM,IE(1),IW(1)
cc     *,                VT3(1,1),VT1(1,1),VT2(1,1))
cc
cc          DO I =1,IM*(JM-1)
cc           QOI(I,1,L,K) =  QOI(I,1,L,K) - VT3(I,1)*D2PIN(I,1)
cc          ENDDO
cc        ENDDO
cc      ENDIF
c      IF(KM.GT.0) THEN
c        IF (JY1.EQ.1) THEN
c        CALL TRACZONALADV(TWO,DT,IM,JM,LM,KM,L,QOM,QOI,
c     *                    PIM,PII,USB,IW,IE,D2PIN)
c        CALL TRACMERIDADV(TWO,DT,IM,JM,LM,KM,L,QOM,QOI,
c     *                    PIM,PII,VSB,D2PIN)
c        ELSE
c        CALL TRACMERIDADV(TWO,DT,IM,JM,LM,KM,L,QOM,QOI,
c     *                    PIM,PII,VSB,D2PIN)
c        CALL TRACZONALADV(TWO,DT,IM,JM,LM,KM,L,QOM,QOI,
c     *                    PIM,PII,USB,IW,IE,D2PIN)
c        ENDIF
c      ENDIF

C  HORIZONTAL ADVECTION OF POTENTIAL TEMPERATURE

      CALL HADVECT(USB(1,1),VSB(1,1),POB(1,1,L)
     *,           IM,JM,IE(1),IW(1)
     *,           VT3(1,1),VT1(1,1),VT2(1,1))

      DO I =1,IM*(JM-1)
       POI(I,1,L) =  POI(I,1,L) - VT3(I,1)*D2PIN(I,1)
      ENDDO

C  COMPUTE CONVERGENCE

      DO I =1,IM*(JM-1)
       PSD(I,1,L) = ( (USB(IW(I),1)-USB(I,1))
     *              + (VSB(I    ,1)-VSB(I,2)) ) * (D2PIN(I,1)*DSG(L))
      ENDDO

C   COMPUTE VORTICITY

      DO I =1,IM*(JM-1)
       VT2(I,1) = DXUIJ(I,1) * UOB(I,1,L)
      ENDDO

      DO I =1,IM*(JM-2)
       ZOB(I,2) = ( (VOB(IE(I),2,L) - VOB(I,2,L)) * DYVIJ(I,2)
     *          +   (VT2(I    ,1  ) - VT2(I,2  ))   ) * D2ZIN(I,2)
      ENDDO

      SUMSO = SDOT(IM,UOB(1,1   ,L),1,DXUIJ(1,1   ),1) * (-ONE/FLOAT(IM)
     >)
      SUMNO = SDOT(IM,UOB(1,JM-1,L),1,DXUIJ(1,JM-1),1) * ( ONE/FLOAT(IM)
     >)

      DO I =1,IM
       ZOB(I, 1) = SUMSO * D2ZIN(1, 1)
       ZOB(I,JM) = SUMNO * D2ZIN(1,JM)
      ENDDO

      IF(DIAG) THEN
       DO I=1,IM*JM
        VOR(I,1,L) = ZOB(I,1)
        OMG(I,1,L) = ZERO
       ENDDO
      ENDIF

      DO I =1,IM*JM
       ZOB(I,1) = (ZOB(I,1)+FFFIJ(I,1)) * ACH(I,1)
      ENDDO

C  COMPUTE VORTICITY COEFFICIENTS

      DO I =1,IM*(JM-1)

       EPS(I,1) = (ZOB(I    ,1) + ZOB(I    ,2))

       LAM(I,1) = (ZOB(IW(I),1) + ZOB(IE(I),2))
       MUU(I,1) = (ZOB(IW(I),2) + ZOB(IE(I),1))

      ENDDO
       
      DO I =1,IM*(JM-2)

       BET(I,1) = (EPS(I,1)+ZOB(IW(I),2))
       GAM(I,2) = (EPS(I,2)+ZOB(IW(I),2))
       DEL(I,2) = (EPS(I,2)+ZOB(IE(I),2))
       ALF(I,1) = (EPS(I,1)+ZOB(IE(I),2))

       BET(I,1) = (THREE*HALF)*BET(I,1)-HALF*(LAM(I,1)+ZOB(I,3))
       GAM(I,2) = (THREE*HALF)*GAM(I,2)-HALF*(MUU(I,2)+ZOB(I,1))
       DEL(I,2) = (THREE*HALF)*DEL(I,2)-HALF*(LAM(I,2)+ZOB(I,1))
       ALF(I,1) = (THREE*HALF)*ALF(I,1)-HALF*(MUU(I,1)+ZOB(I,3))

      ENDDO


      DO I =1,IM*(JM-2)
       MUU(I,2) = HALF*(ZOB(I    ,1) - ZOB(I    ,3))
      ENDDO

      DO I =1,IM*(JM-4)
       LAM(I,3) = HALF*(ZOB(IE(I),3) - ZOB(IW(I),3))
      ENDDO

      ST1 = DYPIJ(1,   1)/(DYPIJ(1,   1)+DYPIJ(1,   2))
      ST2 = DYPIJ(1,JM-1)/(DYPIJ(1,JM-1)+DYPIJ(1,JM-2))

      DO I =1,IM
C       LAM(I,   2) = (THREE/FIVE)*LAM(I,   3)
C       LAM(I,JM-1) = (THREE/FIVE)*LAM(I,JM-2)
       LAM(I,   2) = ST1*LAM(I,   3)
       LAM(I,JM-1) = ST2*LAM(I,JM-2)
      ENDDO


C  U INCREMENT

      DO I =1,IM*(JM-2)
       EPS(I,1) = BET(I,1) * VSB(I    ,2)
     *          + ALF(I,1) * VSB(IE(I),2)

     *          - LAM(I,2) * USB(I    ,2)

      ENDDO

      DO I =1,IM
       EPS(I,JM-1) = ZERO
      ENDDO

      DO I =1,IM*(JM-2)
       EPS(I,2) = EPS(I,2)
     *          + GAM(I,2) * VSB(I    ,2)
     *          + DEL(I,2) * VSB(IE(I),2)

     *          + LAM(I,2) * USB(I    ,1)

      ENDDO

      DO I =1,IM*(JM-1)
       UOI(I,1,L) = UOI(I,1,L) + (THIRD*FOURTH)*EPS(I,1)*DXUIN(I,1)
      ENDDO

C   V INCREMENT

      DO I =1,IM*(JM-2)
       ST1 = BET(I    ,1) * USB(I    ,1)
     *     + GAM(I    ,2) * USB(I    ,2)
     *     + DEL(IW(I),2) * USB(IW(I),2)
     *     + ALF(IW(I),1) * USB(IW(I),1)

     *     + MUU(I    ,2) * VSB(IE(I),2)
     *     - MUU(IW(I),2) * VSB(IW(I),2)


       VOI(I,2,L) = VOI(I,2,L) - (THIRD*FOURTH)*ST1*DYVIN(I,2)
      ENDDO

1000  CONTINUE


C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C***************** END FIRST L LOOP *************************

C***************** BEGIN L-CRITICAL SECTION *****************
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

C  VERTICAL INTERGRAL OF CONTINUITY EQUATION

      DO L =2,LM
CMIC$ DO ALL AUTOSCOPE VECTOR
       DO I =1,IM*(JM-1)
        PSD(I,1,L) = PSD(I,1,L-1) + PSD(I,1,L)
       ENDDO
      ENDDO

C  SURFACE PRESSURE TENDENCY

CMIC$ DO ALL AUTOSCOPE VECTOR
      DO I =1,IM*(JM-1)
        PII(I,1) = PII(I,1) + PSD(I,1,LM)
      ENDDO

C /* DDX AND DDY */ ARE USED IN OMEGA CALCULATION

      IF(DIAG) THEN
CMIC$ DO ALL AUTOSCOPE VECTOR
       DO I =1,IM*(JM-1)
        DDX(I,1) = (PIB(IE(I),1)-PIB(I,1))*DYUIJ(I,1)*PBI(I,1)
       ENDDO

CMIC$ DO ALL AUTOSCOPE VECTOR
       DO I =1,IM*(JM-2)
        DDY(I,2) = (PIB(I    ,2)-PIB(I,1))*DXVIJ(I,2)*PBJ(I,2)
       ENDDO
      ENDIF

C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C***************** END L-CRITICAL SECTION *******************

C***************** BEGIN SECOND L LOOP **********************
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

      DO 2000 LL=1,2
CMIC$ DO ALL AUTOSCOPE PRIVATE(VT1,VT2,VT3,VT4,ALF,BET)
      DO 2000 L=LL,LM-1,2

C  VERTICAL ADVECTION

      DO I =1,IM*(JM-1)
       PSD(I,1,L) = PSD(I,1,L) - PSD(I,1,LM)*SIG(L+1)
      ENDDO

      IF(PTOP.EQ.ZERO) THEN
       DO I =1,IM*(JM-1)
        ALF(I,1) = (PRJ(L+1)-PRH(L))/(PRH(L+1)-PRH(L))
       ENDDO
      ELSE
       DO I =1,IM*(JM-1)
        ALF(I,1) = (PKH(I,1,L+1)-PKL(I,1,L))
     *           / (PKL(I,1,L+1)-PKL(I,1,L))
       ENDDO
      ENDIF

      DO I =1,IM*(JM-1)
       ST1 = ( (POB(I,1,L  ) - POB(I,1,L+1)) * ALF(I,1)
     *       +  POB(I,1,L+1)                            )
     *       *  PSD(I,1,L  )

       POI(I,1,L  ) = POI(I,1,L  ) - ST1 * (ONE/DSG(L  ))
       POI(I,1,L+1) = POI(I,1,L+1) + ST1 * (ONE/DSG(L+1))
      ENDDO

c      OLD VERSION OF TRACER VERTICAL ADVECTION
c      IF(KM.GT.0) THEN
c       DO K =1,KM
c        DO I =1,IM*(JM-1)
c         ST1 = ( (QOB(I,1,L,K) + QOB(I,1,L+1,K)) * HALF )
c     *         *  PSD(I,1,L  )

c         QOI(I,1,L  ,K) = QOI(I,1,L  ,K) - ST1 * (ONE/DSG(L  ))
c         QOI(I,1,L+1,K) = QOI(I,1,L+1,K) + ST1 * (ONE/DSG(L+1))
c        ENDDO
c       ENDDO
c      ENDIF

      DO I =1,IM*(JM-1)
       ST1 = (UOB(I    ,1,L+1) - UOB(I,1,L) )
     *     * (PSD(IE(I),1,L  ) + PSD(I,1,L) )
     *     *  BIP(I    ,1    )

       UOI(I,1,L  ) = UOI(I,1,L  ) - ST1 * (FOURTH/DSG(L  ))
       UOI(I,1,L+1) = UOI(I,1,L+1) - ST1 * (FOURTH/DSG(L+1))
      ENDDO

      DO I =1,IM*(JM-2)
       ST1 = (VOB(I,2,L+1) - VOB(I,2,L) )
     *     * (PSD(I,2,L  ) + PSD(I,1,L  ) )
     *     *  BJP(I,2    )

       VOI(I,2,L  ) = VOI(I,2,L  ) - ST1 * (FOURTH/DSG(L  ))
       VOI(I,2,L+1) = VOI(I,2,L+1) - ST1 * (FOURTH/DSG(L+1))
      ENDDO

      IF(DIAG) THEN

       DO I=1,IM*(JM-1)
        VT1(I,1) = UOB(I,1,L  )*DDX(I,1)
        VT3(I,1) = UOB(I,1,L+1)*DDX(I,1)
       ENDDO

       DO I=1,IM*(JM-2)
        VT2(I,2) = VOB(I,2,L  )*DDY(I,2)
        VT4(I,2) = VOB(I,2,L+1)*DDY(I,2)
       ENDDO

       DO I=1,IM
        VT2(I, 1) = ZERO
        VT2(I,JM) = ZERO
        VT4(I, 1) = ZERO
        VT4(I,JM) = ZERO
       ENDDO

      IF(PTOP.EQ.ZERO) THEN
       DO I =1,IM*(JM-1)
        ALF(I,1) = (PRJ(L+1)-PRH(L  ))/(PRJ(L+1)-PRJ(L  ))
        BET(I,1) = (PRH(L+1)-PRJ(L+1))/(PRJ(L+2)-PRJ(L+1))
       ENDDO
      ELSE
       DO I =1,IM*(JM-1)
        ALF(I,1) = (PKH(I,1,L+1)-PKL(I,1,L  ))
     *           / (PKH(I,1,L+1)-PKH(I,1,L  ))
        BET(I,1) = (PKL(I,1,L+1)-PKH(I,1,L+1))
     *           / (PKH(I,1,L+2)-PKH(I,1,L+1))
       ENDDO
      ENDIF

       DO I=1,IM*(JM-1)
        ST1 = SIG(L+1)*PSD(I,1,LM) + PSD(I,1,L)
        ST2 = HALF*(VT1(IW(I),1)+VT1(I,1)+VT2(I,2)+VT2(I,1))
     *      * SIG(L+1)*D2PIN(I,1)*PIV(I,1)
        OMG(I,1,L  ) = OMG(I,1,L  ) + ALF(I,1) * (ST1 + ST2)
        ST2 = HALF*(VT3(IW(I),1)+VT3(I,1)+VT4(I,2)+VT4(I,1))
     *      * SIG(L+1)*D2PIN(I,1)*PIV(I,1)
        OMG(I,1,L+1) = OMG(I,1,L+1) + BET(I,1) * (ST1 + ST2)
       ENDDO

       IF(L.EQ.LM-1) THEN
        IF(PTOP.EQ.ZERO) THEN
          DO I =1,IM*(JM-1)
           BET(I,1) = (PRJ(L+2)-PRH(L+1))/(PRJ(L+2)-PRJ(L+1))
          ENDDO
         ELSE
          DO I =1,IM*(JM-1)
           BET(I,1) = (PKH(I,1,L+2)-PKL(I,1,L+1))
     *              / (PKH(I,1,L+2)-PKH(I,1,L+1))
          ENDDO
        ENDIF
        DO I=1,IM*(JM-1)
         ST1 =  PSD(I,1,LM)
         ST2 = HALF*(VT3(IW(I),1)+VT3(I,1)+VT4(I,2)+VT4(I,1))
     *       * D2PIN(I,1)*PIV(I,1)
         OMG(I,1,L+1) = OMG(I,1,L+1) + BET(I,1) * (ST1 + ST2)
        ENDDO
       ENDIF


      ENDIF

2000  CONTINUE

C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C***************** END SECOND L LOOP ************************

      IF (KM.GT.0) THEN
        IF (JY1.EQ.1) THEN
        CALL TRACZONALADV(TWO,DT,IM,JM,LM,KM,QOM,QOI,
     *                    PIB,PII,USBSAV,IW,IE,D2PIN)
        CALL TRACMERIDADV(TWO,DT,IM,JM,LM,KM,QOM,QOI,
     *                    PIB,PII,VSBSAV,D2PIN)
        ELSE
        CALL TRACMERIDADV(TWO,DT,IM,JM,LM,KM,QOM,QOI,
     *                    PIB,PII,VSBSAV,D2PIN)
        CALL TRACZONALADV(TWO,DT,IM,JM,LM,KM,QOM,QOI,
     *                    PIB,PII,USBSAV,IW,IE,D2PIN)
        ENDIF
        CALL TRACVERTADV(TWO,DT,IM,JM,LM,KM,DSG,PSD,QOM,QOI,QOI2,PIM,
     *                   PIB)
      ENDIF

C***************** BEGIN L-CRITICAL SECTION *****************
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

C   AVRX OF SCALARS

      CALL AVRX(POI(1,1,1  ),IM,JM,JDIM,LM,SU(1,1))

      IF(KM.GT.0) THEN
       DO K =1,KM
        CALL AVRX(QOI(1,1,1,K),IM,JM,JDIM,LM,SU(1,1))
       ENDDO
      ENDIF

      IF(DIAG .AND. LM.GT.1) THEN
       CALL AVRX(OMG(1,1,1  ),IM,JM,JDIM,LM,SU(1,1))
      ENDIF

      CALL AVRX(PII(1,1),IM,JM,JM,1,SU(1,1))

C  AT THIS POINT THETA, Q, AND PS INCREMENTS ARE COMPLETE
C  AND CAN BE USED FOR COMPUTING THE PRESSURE GRADIENT
C  IN THE ECONOMICAL EXPLICIT SCHEME.

      IF(LEAP) THEN
CMIC$ DO ALL AUTOSCOPE VECTOR
       DO I=1,IM*(JM-1)
        ST1 = PIM(I,1) + DT*PII(I,1)
        VT3(I,1) = ONE / ST1
        VT1(I,1) = ALPHA*(PIM(I,1)+ST1) + (ONE-TWO*ALPHA)*PIB(I,1)
       ENDDO
      ELSE
CMIC$ DO ALL AUTOSCOPE VECTOR
       DO I=1,IM*(JM-1)
        VT1(I,1) = PIB(I,1)
       ENDDO
      ENDIF

C  GRADIENT OF SURFACE PRESSURE

CMIC$ DO ALL AUTOSCOPE VECTOR
      DO I=1,IM*(JM-1)
       PHI(I,1,LM) = PHS(I,1)
      ENDDO


CMIC$ DO ALL AUTOSCOPE VECTOR
      DO I=1,IM*(JM-1)
       DPX(I,1) = (VT1(IE(I),1) - VT1(I,1))
      ENDDO

CMIC$ DO ALL AUTOSCOPE VECTOR
      DO I =1,IM*(JM-2)
       DPY(I,2) = (VT1(I    ,2) - VT1(I,1))
      ENDDO

C   HYDROSTATIC EQUATION


      DO L=LM,2,-1

       IF(LEAP) THEN
CMIC$ DO ALL AUTOSCOPE VECTOR
        DO I =1,IM*(JM-1)
         VT1(I,1) = ( POM(I,1,L)*PIM(I,1)
     *            + DT*POI(I,1,L) )*VT3(I,1)
         VT1(I,1) = ALPHA*(POM(I,1,L) + VT1(I,1))
     *            + (ONE-TWO*ALPHA)*POB(I,1,L)
        ENDDO
       ELSE
CMIC$ DO ALL AUTOSCOPE VECTOR
        DO I =1,IM*(JM-1)
         VT1(I,1) = POB(I,1,L)
        ENDDO
       ENDIF

       IF(PTOP.EQ.ZERO) THEN
CMIC$ DO ALL AUTOSCOPE VECTOR
        DO I =1,IM*(JM-1)
         ST1 = VT1(I,1)*PIK(I,1)
         PHI(I,1,L-1) = PHI(I,1,L) + ST1*(CP*(PRJ(L+1)-PRJ(L)))
         PHI(I,1,L  ) = PHI(I,1,L) + ST1*(CP*(PRJ(L+1)-PRH(L)))
        ENDDO
       ELSE
CMIC$ DO ALL AUTOSCOPE VECTOR
        DO I =1,IM*(JM-1)
         ST1 = VT1(I,1)
         PHI(I,1,L-1) = PHI(I,1,L)
     *                + ST1*(CP*(PKH(I,1,L+1)-PKH(I,1,L)))
         PHI(I,1,L  ) = PHI(I,1,L)
     *                + ST1*(CP*(PKH(I,1,L+1)-PKL(I,1,L)))
        ENDDO
       ENDIF

      ENDDO

       IF(LEAP) THEN
CMIC$ DO ALL AUTOSCOPE VECTOR
        DO I =1,IM*(JM-1)
         VT1(I,1) = ( POM(I,1,1)*PIM(I,1)
     *            + DT*POI(I,1,1) )*VT3(I,1)
         VT1(I,1) = ALPHA*(POM(I,1,1)
     *            + VT1(I,1)) + (ONE-TWO*ALPHA)*POB(I,1,1)
        ENDDO
       ELSE
CMIC$ DO ALL AUTOSCOPE VECTOR
        DO I =1,IM*(JM-1)
         VT1(I,1) = POB(I,1,1)
        ENDDO
       ENDIF

       IF(PTOP.EQ.ZERO) THEN
CMIC$ DO ALL AUTOSCOPE VECTOR
        DO I =1,IM*(JM-1)
         ST1 = VT1(I,1)*PIK(I,1)
         PHI(I,1,1  ) = PHI(I,1,1) + ST1*(CP*(PRJ(2)-PRH(1)))
        ENDDO
       ELSE
CMIC$ DO ALL AUTOSCOPE VECTOR
        DO I =1,IM*(JM-1)
         ST1 = VT1(I,1)
         PHI(I,1,1  ) = PHI(I,1,1)
     *                + ST1*(CP*(PKH(I,1,2)-PKL(I,1,1)))
        ENDDO
       ENDIF

C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C***************** END L-CRITICAL SECTION *******************

C***************** BEGIN THIRD L LOOP ***********************
Cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

CMIC$ DO ALL AUTOSCOPE PRIVATE(VT1,VT2,GAM)
      DO 3000 L=1,LM

C   KINETIC ENERGY  (SADOURNEY SICKPROOF)
C   THERE IS A 1/2 IN DXDY[UV]


      ST2 = FIVE/SIX

      DO I =1,IM*(JM-3)
       ST1      = HALF*(UOB(I    ,1,L)+UOB(I    ,3,L))
       VT1(I,2) = ( (    ST2)*UOB(I,2,L)*UOB(I,2,L)
     *            + (ONE-ST2)*ST1       *ST1        )*D2UIJ(I,2)
      ENDDO

      DO I =1,IM*(JM-2)
       ST1      = HALF*(VOB(IE(I),2,L)+VOB(IW(I),2,L))
       VT2(I,2) = ( (    ST2)*VOB(I,2,L)*VOB(I,2,L)
     *            + (ONE-ST2)*ST1       *ST1        )*D2VIJ(I,2)
      ENDDO

      DO I=1,IM
       VT1(I,1   ) = UOB(I,   1,L)*UOB(I,   1,L)*D2UIJ(I,   1)
       VT1(I,JM-1) = UOB(I,JM-1,L)*UOB(I,JM-1,L)*D2UIJ(I,JM-1)
      ENDDO

      DO I =1,IM
       VT2(I,1 ) = ZERO
       VT2(I,JM) = ZERO
      ENDDO

      DO  I =1,IM*(JM-1)
       PHI(I,1,L) = PHI(I,1,L)
     *            + HALF*( HALF*(VT1(IW(I),1 ) + VT1(I,1 ))
     *            +        HALF*(VT2(I    ,2 ) + VT2(I,1 ))
     *                           )*D2PIN(I,1 )
      ENDDO

C  PRESSURE GRADIENT FORCE, INCLUDING GRAD OF KE
c  GAM IS (CP THETA DP/DPI)

      IF(PTOP.EQ.ZERO) THEN
       DO I =1,IM*(JM-1)
        GAM(I,1) = ( SIG(L  )*(PRH(L  )-PRJ(L))
     *             + SIG(L+1)*(PRJ(L+1)-PRH(L)) ) * PIK(I,1)
     *           * (CP/DSG(L)) * PIV(I,1) * POB(I,1,L)
       ENDDO
      ELSE
       DO I =1,IM*(JM-1)
        GAM(I,1) = ( SIG(L  )*(PKL(I,1,L  )-PKH(I,1,L))
     *             + SIG(L+1)*(PKH(I,1,L+1)-PKL(I,1,L)) )
     *           * (CP/DSG(L)) * PIV(I,1) * POB(I,1,L)
       ENDDO
      ENDIF

      DO I =1,IM*(JM-1)
       ST1        = (PHI(IE(I),1,L) - PHI(I,1,L)) 
     *            + DPX(I,1) * HALF*(GAM(IE(I),1)+GAM(I,1))
       UOI(I,1,L) = UOI(I,1,L) - ST1*DXUIN(I,1)
      ENDDO

      DO I =1,IM*(JM-2)
       ST1        = (PHI(I    ,2,L) - PHI(I,1,L)) 
     *            + DPY(I,2) * HALF*(GAM(I    ,2)+GAM(I,1))
       VOI(I,2,L) = VOI(I,2,L) - ST1*DYVIN(I,2)
      ENDDO

3000  CONTINUE

C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C***************** END THIRD L LOOP *************************

C  AVRX U INCREMENTS

      CALL AVRX(UOI(1,1,1),IM,JM,JDIM,LM,SU(1,1))
      CALL AVRX(VOI(1,1,1),IM,JM,JDIM,LM,SV(1,1))

C******************************************************
C***************** END DYCORE *************************
C******************************************************

      RETURN
      END

      SUBROUTINE GRIDH(
     *              IM,JM, AE, OMEGA
     *,             DXPIJ, DYPIJ, DXUIJ, DYUIJ, DXVIJ, DYVIJ
     *,             D2PIJ, D2UIJ, D2VIJ, D2ZIJ
     *,             DYVIN, DXUIN, D2PIN, D2ZIN, FFFIJ, SU, SV
     *,             IE, IW
     *                )


C*********************************************************************
C*********************** ARIES   MODEL *******************************
C********************** SUBROUTINE GRIDH *****************************
C********************** 9 FEBRUARY 1987 ******************************
C*********************************************************************

      IMPLICIT NONE

      REAL ONE
      REAL TWO
      REAL HALF
      REAL FOURTH
      REAL ZERO
      REAL PI

      PARAMETER (ONE=1.0)
      PARAMETER (TWO=2.0)
      PARAMETER (ZERO=0.0)
      PARAMETER (HALF=0.5)
      PARAMETER (FOURTH=0.25)
      PARAMETER (PI=3.1415926535898)

      INTEGER NSYM
      REAL SINPNP
      REAL LAMNP

      PARAMETER (NSYM=1)
      PARAMETER (SINPNP=1.)
      PARAMETER (LAMNP=0.)

      INTEGER IM, JM

      REAL AE
      REAL OMEGA


      REAL DXPIJ(IM,JM)
      REAL DYPIJ(IM,JM)
      REAL DXUIJ(IM,JM)
      REAL DYUIJ(IM,JM)
      REAL DYVIJ(IM,JM)
      REAL DXVIJ(IM,JM)

      REAL DYVIN(IM,JM)
      REAL DXUIN(IM,JM)
      REAL D2ZIN(IM,JM)
      REAL D2PIN(IM,JM)

      REAL D2PIJ(IM,JM)
      REAL D2UIJ(IM,JM)
      REAL D2VIJ(IM,JM)
      REAL D2ZIJ(IM,JM)

      REAL FFFIJ(IM,JM)

      REAL SU(IM+2,JM)
      REAL SV(IM+2,JM)

      INTEGER IE(IM*JM)
      INTEGER IW(IM*JM)


C   LOCALS

      REAL COSPNP, FS, FC

      REAL DLAM (IM)
      REAL DPHI (JM)

      REAL SL (IM)
      REAL CL (IM)

      REAL LAMBDA(IM)

      REAL CC  (0:JM)
      REAL SC  (0:JM)

      REAL CV  (JM)
      REAL PHI (0:JM)

      INTEGER I,J, LATS

      REAL SSUM, SDOT

CFPP$ NOVECTOR R

      COSPNP = SQRT(ONE - SINPNP**2)

C   INDEXING

      DO I=1,IM*JM
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
      ENDDO

C  HORIZONTAL GRID

      DO I=1,IM
       DLAM(I) = (TWO*PI)/FLOAT(IM*NSYM)
      ENDDO

      LAMBDA(1) = -PI/FLOAT(NSYM)
      DO I=2,IM
       LAMBDA(I) = LAMBDA(I-1) + DLAM(I-1)
      ENDDO

      DO I=1,IM
       SL(I) = SIN(LAMBDA(I)-LAMNP*(PI/180.))
       CL(I) = COS(LAMBDA(I)-LAMNP*(PI/180.))
      ENDDO

      DO J=1,JM
       DPHI(J) = PI/FLOAT(JM)
      ENDDO

      PHI(0) = -PI*HALF
      CC(0)  = ZERO
      SC(0)  = -ONE

      DO J=1,JM
       PHI(J) =  PHI(J-1) + DPHI(J)
       CC(J)  = COS(PHI(J))
       SC(J)  = SIN(PHI(J))
      ENDDO

      CV( 1) = ZERO
      CV(JM) = ZERO

      DO J=2,JM-1
       CV(J)  = HALF*(CC(J)+CC(J-1))
      ENDDO

C   DXU & DYV ARE DEFINED FROM THE GRID;
C   ALL OTHER FACTORS ARE DEFINED IN TERMS OF THEM.

      DO J=1,JM-1
      DO I=1,IM
       DXUIJ(I,J)   = AE*DLAM(I)*CC(J)
      ENDDO
      ENDDO

      DO J=1,JM
      DO I=1,IM
       DYVIJ(I,J)   = AE*DPHI(J)
      ENDDO
      ENDDO


C   DXP AND DXV

      DO I=1,IM*(JM-1)
       DXPIJ(I,1)   = HALF*(DXUIJ(I,1) + DXUIJ(IW(I),1))
      ENDDO

      DO I=1,IM*(JM-2)
       DXVIJ(I,2)   = HALF*(DXPIJ(I,2) + DXPIJ(I    ,1))
      ENDDO

      DO I=1,IM
       DXVIJ(I, 1) = ZERO
       DXVIJ(I,JM) = ZERO
       DXUIJ(I,JM) = ZERO
       DXPIJ(I,JM) = ZERO
      ENDDO


C  DYU & DYP

      DO I=1,IM*(JM-3)
       DYPIJ(I,2)   = HALF*(DYVIJ(I,2) + DYVIJ(I,3))
      ENDDO

      DO I=1,IM
       DYPIJ(I,   1) = AE*(DPHI( 1) + HALF*DPHI(   2))
       DYPIJ(I,JM-1) = AE*(DPHI(JM) + HALF*DPHI(JM-1))
       DYPIJ(I,  JM) = ZERO
       DYUIJ(I,  JM) = ZERO
      ENDDO

      DO I=1,IM*(JM-1)
       DYUIJ(I,1)   = DYPIJ(I,1)
      ENDDO


C  AREA FACTORS
 
      DO I=1,IM*(JM-3)
       D2PIJ(I,2) = HALF*(DXVIJ(I,2)+DXVIJ(I,3)) * DYPIJ(I,2)
      ENDDO

      DO I=1,IM
       D2PIJ(I,   1) = DXPIJ(I,   1) * DYVIJ(I, 1)
       D2PIJ(I,JM-1) = DXPIJ(I,JM-1) * DYVIJ(I,JM)
       D2UIJ(I,   1) = D2PIJ(I,   1)
       D2UIJ(I,JM-1) = D2PIJ(I,JM-1)
       D2UIJ(I,  JM) = ZERO
       D2PIJ(I,  JM) = ZERO
      ENDDO
 
      DO I=1,IM*(JM-3)
       D2UIJ(I,2) = DXUIJ(I,2) * DYUIJ(I,2)
      ENDDO

      DO I=1,IM*(JM-2)
       D2VIJ(I,2) = DXVIJ(I,2) * DYVIJ(I,2)
      ENDDO

      DO I=1,IM*(JM-2)
       D2ZIJ(I,2) = FOURTH*(  D2PIJ(I    ,1)+D2PIJ(I    ,2)
     *                      + D2PIJ(IE(I),1)+D2PIJ(IE(I),2) )
      ENDDO

      DO I=1,IM
       D2VIJ(I, 1) = ZERO
       D2VIJ(I,JM) = ZERO
       D2ZIJ(I, 1) = (HALF/FLOAT(IM))*SSUM(IM,D2PIJ(1,   1),1)
       D2ZIJ(I,JM) = (HALF/FLOAT(IM))*SSUM(IM,D2PIJ(1,JM-1),1)
      ENDDO


C  CORIOLIS PARAMETER

      DO J=2,JM-1
      DO I=1,IM
       FS = (CC(J)*DXUIJ(I,J)-CC(J-1)*DXUIJ(I,J-1))
       FC = ( DYVIJ(I,J)*(SL(IE(I))-SL(I)) + HALF*(CL(IE(I))+CL(I))
     *       * (DXUIJ(I,J)*SC(J)-DXUIJ(I,J-1)*SC(J-1)) )
       FFFIJ(I,J) = -(AE*OMEGA/D2ZIJ(I,J)) * (SINPNP*FS - COSPNP*FC)
      ENDDO
      ENDDO

      FS = CC(1)*(ONE/FLOAT(IM))*SSUM(IM,DXUIJ(1,1),1)
      FC = SC(1)*(ONE/FLOAT(IM))*SDOT(IM,DXUIJ(1,1),1,CL(1),1)

      DO I=1,IM
       FFFIJ(I, 1) = -(AE*OMEGA/D2ZIJ(I, 1))*(SINPNP*FS - COSPNP*FC)
      ENDDO

      FS = CC(JM-1)*(ONE/FLOAT(IM))*SSUM(IM,DXUIJ(1,JM-1),1)
      FC = SC(JM-1)*(ONE/FLOAT(IM))*SDOT(IM,DXUIJ(1,JM-1),1,CL(1),1)

      DO I=1,IM
       FFFIJ(I,JM) =  (AE*OMEGA/D2ZIJ(I,JM))*(SINPNP*FS - COSPNP*FC)
      ENDDO

C  PRINT ARRAYS
c
c      PRINT *, '  DYCORE: Version @(#)dycore.F	1.3 6/1/95 '
c
c      PRINT *, "  DYCORE: COMPILED ", "Wed Nov 20 17:46:40 EST 1996"
c
c      PRINT *
c
c      PRINT 400, (PHI(J)*180./PI, CC(J), SC(J), FFFIJ(1,J)*1.E4
c     *,       D2PIJ(1,J)*1.E-6, D2UIJ(1,J)*1.E-6, D2VIJ(1,J)*1.E-6
c     *,       J=1,JM)
c      PRINT *
c      PRINT 500, (
c     *        DXVIJ(1,J)*1.E-3, DXPIJ(1,J)*1.E-3, DXUIJ(1,J)*1.E-3
c     *,       DYPIJ(1,J)*1.E-3, DYVIJ(1,J)*1.E-3
c     *,       J=1,JM)
c400   FORMAT(7F12.4)
c500   FORMAT(5F12.4)
c
C  PRE-COMPUTE INVERSES

      DO J=1,JM-1
      DO I=1,IM
       DXUIN(I,J) = ONE  / DXUIJ(I,J)
       D2PIN(I,J) = ONE  / D2PIJ(I,J)
      ENDDO
      ENDDO

      DO J=2,JM-1
      DO I=1,IM
       DYVIN(I,J) = ONE / DYVIJ(I,J)
      ENDDO
      ENDDO

      DO J=1,JM
      DO I=1,IM
       D2ZIN(I,J) = ONE / D2ZIJ(I,J)
      ENDDO
      ENDDO

C  AVRX FILTER ARRAYS. THESE ASSUME EQUALLY SPACED LAT-LON GRID.

      DO I=1,(IM+2)*JM
       SU(I,1) = ONE
       SV(I,1) = ONE
      ENDDO

      LATS  = (JM/4)

      DO J=1,JM-1
       DO I=3,IM+2
        SU(I,J) = AMIN1(
     *             (CC(J)/CC(LATS))/SIN(INT((I-1)/2)*PI/FLOAT(IM))
     *             ,ONE  ) ** 2
       ENDDO
      ENDDO

      DO J=2,JM-1
       DO I=3,IM+2
        SV(I,J) = AMIN1(
     *             (CV(J)/CV(LATS))/SIN(INT((I-1)/2)*PI/FLOAT(IM))
     *             ,ONE  ) ** 2
       ENDDO
      ENDDO

C   WRITE FILTER ARRAYS
C
C      WRITE(6,201)  '  THETA   ', '   SU(1)  ','   SU(2)  '
C     *,             '   SU(3)  ', '   SU(4)  '
C     *,             '   SU(5)  ', '   SU(6)  '
C     *,             ' SU(MM-5) ', ' SU(MM-4) '
C     *,             ' SU(MM-3) ', ' SU(MM-2)   '
C     *,             ' SU(MM-1) ', ' SU(MM) '
C      WRITE(6,200) (PHI(J)*(180./PI),SU(4,J),SU(6,J),SU(8,J)
C     *,             SU(10,J),SU(12,J),SU(14,J),SU(IM-8,J)
C     *,             SU(IM-6,J),SU(IM-4,J),SU(IM-2,J),SU(IM,J)
C     *,             SU(IM+2,J),J=1,JM-1)
C
C      WRITE(6,201)  '  THETA   ', '   SV(1)  ','   SV(2)  '
C     *,             '   SV(3)  ', '   SV(4)  '
C     *,             '   SV(5)  ', '   SV(6)  '
C     *,             ' SV(MM-5) ', ' SV(MM-4) '
C     *,             ' SV(MM-3) ', ' SV(MM-2)   '
C     *,             ' SV(MM-1) ', ' SV(MM) '
C      WRITE(6,200) ((PHI(J)-DPHI(J)/TWO)*(180./PI),SV(4,J),SV(6,J)
C     *,             SV(8,J),SV(10,J),SV(12,J),SV(14,J),SV(IM-8,J)
C     *,             SV(IM-6,J),SV(IM-4,J),SV(IM-2,J),SV(IM,J)
C     *,             SV(IM+2,J),J=2,JM-1)
C
C201   FORMAT(/,13(1X,A10),/)
C200   FORMAT(13(1X,F10.6))


      RETURN
      END

      SUBROUTINE GRIDV(LM, RKAP
     *,                SIG, DSG, PRJ, PRH, P00
     *                )

      IMPLICIT NONE

      REAL P00
c ...      PARAMETER (P00=1000.)
C      PARAMETER (P00=7.6)

      INTEGER LM

      REAL RKAP

      REAL SIG(LM+1)
      REAL PRJ(LM+1)
      REAL DSG(LM)
      REAL PRH(LM)

      INTEGER L

      PRJ(1) = 0.0

      DO L=1,LM
       DSG(L  ) = SIG(L+1) - SIG(L)
       PRJ(L+1) = (SIG(L+1)/P00)**RKAP
       PRH(L  ) = (SIG(L+1)*PRJ(L+1)-SIG(L)*PRJ(L))/(DSG(L)*(RKAP+1.0))
      ENDDO

      RETURN
      END


      SUBROUTINE AVRX(U,IM,JM,JDIM,LM,S)

      IMPLICIT NONE

      INTEGER   FORWARD, BACKWARD
      PARAMETER (FORWARD=-1, BACKWARD=1)


      REAL ONE
      REAL ZERO

      PARAMETER (ONE=1.0)
      PARAMETER (ZERO=0.0)

C  ARGUMENTS

      INTEGER IM, JM, JDIM, LM

C Changed to real*8 (input variable)
      REAL U(IM,JDIM,LM)
      REAL S(IM+2,JM)


C   SCRATCH SPACE

      REAL Z(IM+2,64  )
      REAL Y(IM+2,64  )
      REAL B(IM+2,128 )

      INTEGER NLAT, I, J, L, J1(JM*LM), L1(JM*LM) 
      INTEGER NFFTS,LOTS,NLOT,OLOT,LL,II
      INTEGER IM0
      DATA IM0/0/

      INTEGER         IX(19)
      REAL, allocatable :: TR(:)
      integer :: ierror

      SAVE TR
      SAVE IX

CFPP$ NOVECTOR R

      IF (IM.NE.IM0) THEN
        allocate(TR(IM*2),stat=ierror)
        if(ierror.ne.0) then
          write(6,'("Could not allocate TR array in avrx.")')
          stop
        end if

        CALL FFTFAX(IM,IX(1),TR(1))
        IM0=IM
      ENDIF


C   COPY LATITUDES TO BE SMOOTHED INTO CONTIGUOUS BUFFER

      NFFTS = 0
      DO J=1,JM-1
       DO L=1,LM
        IF(S(IM+2,J).LT..9999) THEN
          NFFTS = NFFTS + 1
          J1(NFFTS) = J
          L1(NFFTS) = L
        ENDIF
       ENDDO
      ENDDO

      IF(NFFTS.EQ.0) RETURN

      NLOT = (NFFTS-1)/64 + 1
      OLOT = MOD(NFFTS-1,64)+1


CMIC$ DO ALL AUTOSCOPE PRIVATE(Z,Y,B,LL,II) SHARED(TR,IX,S,U,J1,L1)
      DO LOTS=1,NLOT

        IF(LOTS.EQ.NLOT) THEN
          LL=OLOT
        ELSE
          LL=64
        ENDIF

        DO L=1,LL
          II = (LOTS-1)*64+L
          DO I=1,IM
           Z(I,L) = U(I,J1(II),L1(II))
          ENDDO
          Z(IM+1,L) = ZERO
          Z(IM+2,L) = ZERO
          DO I=1,IM+2
           Y(I,L) = S(I,J1(II)       )
          ENDDO
        ENDDO


C   TRANSFORM AND FILTER

        CALL RFFTMLT(Z(1,1),B(1,1)
     *,            TR(1),IX(1),1,IM+2,IM,LL,FORWARD)

        DO I=1,(IM+2)*LL
         Z(I,1) = Z(I,1)*Y(I,1)
        ENDDO

        CALL RFFTMLT(Z(1,1),B(1,1)
     *,            TR(1),IX(1),1,IM+2,IM,LL,BACKWARD)


C   COPY BACK FROM CONTIGUOUS BUFFER

        DO L=1,LL
          II = (LOTS-1)*64+L
          DO I=1,IM
           U(I,J1(II),L1(II)) = Z(I,L)
          ENDDO
        ENDDO

      ENDDO


      RETURN
      END

C FFT991 RENAMED RFFTMLT
C
C PURPOSE      PERFORMS MULTIPLE FAST FOURIER TRANSFORMS.  THIS PACKAGE
C              WILL PERFORM A NUMBER OF SIMULTANEOUS REAL/HALF-COMPLEX
C              PERIODIC FOURIER TRANSFORMS OR CORRESPONDING INVERSE
C              TRANSFORMS, I.E.  GIVEN A SET OF REAL DATA VECTORS, THE
C              PACKAGE RETURNS A SET OF 'HALF-COMPLEX' FOURIER
C              COEFFICIENT VECTORS, OR VICE VERSA.  THE LENGTH OF THE
C              TRANSFORMS MUST BE AN EVEN NUMBER GREATER THAN 4 THAT HAS
C              NO OTHER FACTORS EXCEPT POSSIBLY POWERS OF 2, 3, AND 5.
C              THIS IS AN ALL FORTRAN VERSION OF THE CRAYLIB PACKAGE
C              THAT IS MOSTLY WRITTEN IN CAL.
C
C              THE PACKAGE FFT99F CONTAINS SEVERAL USER-LEVEL ROUTINES:
C
C            SUBROUTINE FFTFAX
C                AN INITIALIZATION ROUTINE THAT MUST BE CALLED ONCE
C                BEFORE A SEQUENCE OF CALLS TO THE FFT ROUTINES
C                (PROVIDED THAT N IS NOT CHANGED).
C
C            SUBROUTINES FFT99 AND FFT991
C                TWO FFT ROUTINES THAT RETURN SLIGHTLY DIFFERENT
C                ARRANGEMENTS OF THE DATA IN GRIDPOINT SPACE.
C
C
C ACCESS       THIS FORTRAN VERSION MAY BE ACCESSED WITH
C
C                   *FORTRAN,P=XLIB,SN=FFT99F
C
C              TO ACCESS THE CRAY OBJECT CODE, CALLING THE USER ENTRY
C              POINTS FROM A CRAY PROGRAM IS SUFFICIENT.  THE SOURCE
C              FORTRAN AND CAL CODE FOR THE CRAYLIB VERSION MAY BE
C              ACCESSED USING
C
C                   FETCH P=CRAYLIB,SN=FFT99
C                   FETCH P=CRAYLIB,SN=CAL99
C
C USAGE        LET N BE OF THE FORM 2**P * 3**Q * 5**R, WHERE P .GE. 1,
C              Q .GE. 0, AND R .GE. 0.  THEN A TYPICAL SEQUENCE OF
C              CALLS TO TRANSFORM A GIVEN SET OF REAL VECTORS OF LENGTH
C              N TO A SET OF 'HALF-COMPLEX' FOURIER COEFFICIENT VECTORS
C              OF LENGTH N IS
C
C                   DIMENSION IFAX(13),TRIGS(3*N/2+1),A(M*(N+2)),
C                  +          WORK(M*(N+1))
C
C                   CALL FFTFAX (N, IFAX, TRIGS)
C                   CALL FFT99 (A,WORK,TRIGS,IFAX,INC,JUMP,N,M,ISIGN)
C
C              SEE THE INDIVIDUAL WRITE-UPS FOR FFTFAX, FFT99, AND
C              FFT991 BELOW, FOR A DETAILED DESCRIPTION OF THE
C              ARGUMENTS.
C
C HISTORY      THE PACKAGE WAS WRITTEN BY CLIVE TEMPERTON AT ECMWF IN
C              NOVEMBER, 1978.  IT WAS MODIFIED, DOCUMENTED, AND TESTED
C              FOR NCAR BY RUSS REW IN SEPTEMBER, 1980.
C
C-----------------------------------------------------------------------
C
C SUBROUTINE FFTFAX (N,IFAX,TRIGS)
C
C PURPOSE      A SET-UP ROUTINE FOR FFT99 AND FFT991.  IT NEED ONLY BE
C              CALLED ONCE BEFORE A SEQUENCE OF CALLS TO THE FFT
C              ROUTINES (PROVIDED THAT N IS NOT CHANGED).
C
C ARGUMENT     IFAX(13),TRIGS(3*N/2+1)
C DIMENSIONS
C
C ARGUMENTS
C
C ON INPUT     N
C               AN EVEN NUMBER GREATER THAN 4 THAT HAS NO PRIME FACTOR
C               GREATER THAN 5.  N IS THE LENGTH OF THE TRANSFORMS (SEE
C               THE DOCUMENTATION FOR FFT99 AND FFT991 FOR THE
C               DEFINITIONS OF THE TRANSFORMS).
C
C              IFAX
C               AN INTEGER ARRAY.  THE NUMBER OF ELEMENTS ACTUALLY USED
C               WILL DEPEND ON THE FACTORIZATION OF N.  DIMENSIONING
C               IFAX FOR 13 SUFFICES FOR ALL N LESS THAN A MILLION.
C
C              TRIGS
C               A FLOATING POINT ARRAY OF DIMENSION 3*N/2 IF N/2 IS
C               EVEN, OR 3*N/2+1 IF N/2 IS ODD.
C
C ON OUTPUT    IFAX
C               CONTAINS THE FACTORIZATION OF N/2.  IFAX(1) IS THE
C               NUMBER OF FACTORS, AND THE FACTORS THEMSELVES ARE STORED
C               IN IFAX(2),IFAX(3),...  IF FFTFAX IS CALLED WITH N ODD,
C               OR IF N HAS ANY PRIME FACTORS GREATER THAN 5, IFAX(1)
C               IS SET TO -99.
C
C              TRIGS
C               AN ARRAY OF TRIGNOMENTRIC FUNCTION VALUES SUBSEQUENTLY
C               USED BY THE FFT ROUTINES.
C
C-----------------------------------------------------------------------
C
C SUBROUTINE FFT991 (A,WORK,TRIGS,IFAX,INC,JUMP,N,M,ISIGN)
C                       AND
C SUBROUTINE FFT99 (A,WORK,TRIGS,IFAX,INC,JUMP,N,M,ISIGN)
C
C PURPOSE      PERFORM A NUMBER OF SIMULTANEOUS REAL/HALF-COMPLEX
C              PERIODIC FOURIER TRANSFORMS OR CORRESPONDING INVERSE
C              TRANSFORMS, USING ORDINARY SPATIAL ORDER OF GRIDPOINT
C              VALUES (FFT991) OR EXPLICIT CYCLIC CONTINUITY IN THE
C              GRIDPOINT VALUES (FFT99).  GIVEN A SET
C              OF REAL DATA VECTORS, THE PACKAGE RETURNS A SET OF
C              'HALF-COMPLEX' FOURIER COEFFICIENT VECTORS, OR VICE
C              VERSA.  THE LENGTH OF THE TRANSFORMS MUST BE AN EVEN
C              NUMBER THAT HAS NO OTHER FACTORS EXCEPT POSSIBLY POWERS
C              OF 2, 3, AND 5.  THESE VERSION OF FFT991 AND FFT99 ARE
C              OPTIMIZED FOR USE ON THE CRAY-1.
C
C ARGUMENT     A(M*(N+2)), WORK(M*(N+1)), TRIGS(3*N/2+1), IFAX(13)
C DIMENSIONS
C
C ARGUMENTS
C
C ON INPUT     A
C               AN ARRAY OF LENGTH M*(N+2) CONTAINING THE INPUT DATA
C               OR COEFFICIENT VECTORS.  THIS ARRAY IS OVERWRITTEN BY
C               THE RESULTS.
C
C              WORK
C               A WORK ARRAY OF DIMENSION M*(N+1)
C
C              TRIGS
C               AN ARRAY SET UP BY FFTFAX, WHICH MUST BE CALLED FIRST.
C
C              IFAX
C               AN ARRAY SET UP BY FFTFAX, WHICH MUST BE CALLED FIRST.
C
C              INC
C               THE INCREMENT (IN WORDS) BETWEEN SUCCESSIVE ELEMENTS OF
C               EACH DATA OR COEFFICIENT VECTOR (E.G.  INC=1 FOR
C               CONSECUTIVELY STORED DATA).
C
C              JUMP
C               THE INCREMENT (IN WORDS) BETWEEN THE FIRST ELEMENTS OF
C               SUCCESSIVE DATA OR COEFFICIENT VECTORS.  ON THE CRAY-1,
C               TRY TO ARRANGE DATA SO THAT JUMP IS NOT A MULTIPLE OF 8
C               (TO AVOID MEMORY BANK CONFLICTS).  FOR CLARIFICATION OF
C               INC AND JUMP, SEE THE EXAMPLES BELOW.
C
C              N
C               THE LENGTH OF EACH TRANSFORM (SEE DEFINITION OF
C               TRANSFORMS, BELOW).
C
C              M
C               THE NUMBER OF TRANSFORMS TO BE DONE SIMULTANEOUSLY.
C
C              ISIGN
C               = +1 FOR A TRANSFORM FROM FOURIER COEFFICIENTS TO
C                    GRIDPOINT VALUES.
C               = -1 FOR A TRANSFORM FROM GRIDPOINT VALUES TO FOURIER
C                    COEFFICIENTS.
C
C ON OUTPUT    A
C               IF ISIGN = +1, AND M COEFFICIENT VECTORS ARE SUPPLIED
C               EACH CONTAINING THE SEQUENCE:
C
C               A(0),B(0),A(1),B(1),...,A(N/2),B(N/2)  (N+2 VALUES)
C
C               THEN THE RESULT CONSISTS OF M DATA VECTORS EACH
C               CONTAINING THE CORRESPONDING N+2 GRIDPOINT VALUES:
C
C               FOR FFT991, X(0), X(1), X(2),...,X(N-1),0,0.
C               FOR FFT99, X(N-1),X(0),X(1),X(2),...,X(N-1),X(0).
C                   (EXPLICIT CYCLIC CONTINUITY)
C
C               WHEN ISIGN = +1, THE TRANSFORM IS DEFINED BY:
C                 X(J)=SUM(K=0,...,N-1)(C(K)*EXP(2*I*J*K*PI/N))
C                 WHERE C(K)=A(K)+I*B(K) AND C(N-K)=A(K)-I*B(K)
C                 AND I=SQRT (-1)
C
C               IF ISIGN = -1, AND M DATA VECTORS ARE SUPPLIED EACH
C               CONTAINING A SEQUENCE OF GRIDPOINT VALUES X(J) AS
C               DEFINED ABOVE, THEN THE RESULT CONSISTS OF M VECTORS
C               EACH CONTAINING THE CORRESPONDING FOURIER COFFICIENTS
C               A(K), B(K), 0 .LE. K .LE N/2.
C
C               WHEN ISIGN = -1, THE INVERSE TRANSFORM IS DEFINED BY:
C                 C(K)=(1/N)*SUM(J=0,...,N-1)(X(J)*EXP(-2*I*J*K*PI/N))
C                 WHERE C(K)=A(K)+I*B(K) AND I=SQRT(-1)
C
C               A CALL WITH ISIGN=+1 FOLLOWED BY A CALL WITH ISIGN=-1
C               (OR VICE VERSA) RETURNS THE ORIGINAL DATA.
C
C               NOTE: THE FACT THAT THE GRIDPOINT VALUES X(J) ARE REAL
C               IMPLIES THAT B(0)=B(N/2)=0.  FOR A CALL WITH ISIGN=+1,
C               IT IS NOT ACTUALLY NECESSARY TO SUPPLY THESE ZEROS.
C
C EXAMPLES      GIVEN 19 DATA VECTORS EACH OF LENGTH 64 (+2 FOR EXPLICIT
C               CYCLIC CONTINUITY), COMPUTE THE CORRESPONDING VECTORS OF
C               FOURIER COEFFICIENTS.  THE DATA MAY, FOR EXAMPLE, BE
C               ARRANGED LIKE THIS:
C
C FIRST DATA   A(1)=    . . .                A(66)=             A(70)
C VECTOR       X(63) X(0) X(1) X(2) ... X(63) X(0)  (4 EMPTY LOCATIONS)
C
C SECOND DATA  A(71)=   . . .                                  A(140)
C VECTOR       X(63) X(0) X(1) X(2) ... X(63) X(0)  (4 EMPTY LOCATIONS)
C
C               AND SO ON.  HERE INC=1, JUMP=70, N=64, M=19, ISIGN=-1,
C               AND FFT99 SHOULD BE USED (BECAUSE OF THE EXPLICIT CYCLIC
C               CONTINUITY).
C
C               ALTERNATIVELY THE DATA MAY BE ARRANGED LIKE THIS:
C
C                FIRST         SECOND                          LAST
C                DATA          DATA                            DATA
C                VECTOR        VECTOR                          VECTOR
C
C                 A(1)=         A(2)=                           A(19)=
C
C                 X(63)         X(63)       . . .               X(63)
C        A(20)=   X(0)          X(0)        . . .               X(0)
C        A(39)=   X(1)          X(1)        . . .               X(1)
C                  .             .                               .
C                  .             .                               .
C                  .             .                               .
C
C               IN WHICH CASE WE HAVE INC=19, JUMP=1, AND THE REMAINING
C               PARAMETERS ARE THE SAME AS BEFORE.  IN EITHER CASE, EACH
C               COEFFICIENT VECTOR OVERWRITES THE CORRESPONDING INPUT
C               DATA VECTOR.
C
C-----------------------------------------------------------------------
      SUBROUTINE FFT99A(A,WORK,TRIGS,INC,JUMP,N,LOT)
      DIMENSION A(N),WORK(N)
      REAL      TRIGS(N)
C
C     SUBROUTINE FFT99A - PREPROCESSING STEP FOR FFT99, ISIGN=+1
C     (SPECTRAL TO GRIDPOINT TRANSFORM)
C
      NH=N/2
      NX=N+1
      INK=INC+INC
C
C     A(0) AND A(N/2)
      IA=1
      IB=N*INC+1
      JA=1
      JB=2
CDIR$ IVDEP
      DO 10 L=1,LOT
      WORK(JA)=A(IA)+A(IB)
      WORK(JB)=A(IA)-A(IB)
      IA=IA+JUMP
      IB=IB+JUMP
      JA=JA+NX
      JB=JB+NX
   10 CONTINUE
C
C     REMAINING WAVENUMBERS
      IABASE=2*INC+1
      IBBASE=(N-2)*INC+1
      JABASE=3
      JBBASE=N-1
C
      DO 30 K=3,NH,2
      IA=IABASE
      IB=IBBASE
      JA=JABASE
      JB=JBBASE
      C=TRIGS(N+K)
      S=TRIGS(N+K+1)
CDIR$ IVDEP
      DO 20 L=1,LOT
      WORK(JA)=(A(IA)+A(IB))-
     *    (S*(A(IA)-A(IB))+C*(A(IA+INC)+A(IB+INC)))
      WORK(JB)=(A(IA)+A(IB))+
     *    (S*(A(IA)-A(IB))+C*(A(IA+INC)+A(IB+INC)))
      WORK(JA+1)=(C*(A(IA)-A(IB))-S*(A(IA+INC)+A(IB+INC)))+
     *    (A(IA+INC)-A(IB+INC))
      WORK(JB+1)=(C*(A(IA)-A(IB))-S*(A(IA+INC)+A(IB+INC)))-
     *    (A(IA+INC)-A(IB+INC))
      IA=IA+JUMP
      IB=IB+JUMP
      JA=JA+NX
      JB=JB+NX
   20 CONTINUE
      IABASE=IABASE+INK
      IBBASE=IBBASE-INK
      JABASE=JABASE+2
      JBBASE=JBBASE-2
   30 CONTINUE
C
      IF (IABASE.NE.IBBASE) GO TO 50
C     WAVENUMBER N/4 (IF IT EXISTS)
      IA=IABASE
      JA=JABASE
CDIR$ IVDEP
      DO 40 L=1,LOT
      WORK(JA)=2.0*A(IA)
      WORK(JA+1)=-2.0*A(IA+INC)
      IA=IA+JUMP
      JA=JA+NX
   40 CONTINUE
C
   50 CONTINUE
      RETURN
      END
      SUBROUTINE FFT99B(WORK,A,TRIGS,INC,JUMP,N,LOT)
      DIMENSION WORK(N),A(N)
      REAL      TRIGS(N)
C
C     SUBROUTINE FFT99B - POSTPROCESSING STEP FOR FFT99, ISIGN=-1
C     (GRIDPOINT TO SPECTRAL TRANSFORM)
C
      NH=N/2
      NX=N+1
      INK=INC+INC
C
C     A(0) AND A(N/2)
      SCALE=1.0/FLOAT(N)
      IA=1
      IB=2
      JA=1
      JB=N*INC+1
CDIR$ IVDEP
      DO 10 L=1,LOT
      A(JA)=SCALE*(WORK(IA)+WORK(IB))
      A(JB)=SCALE*(WORK(IA)-WORK(IB))
      A(JA+INC)=0.0
      A(JB+INC)=0.0
      IA=IA+NX
      IB=IB+NX
      JA=JA+JUMP
      JB=JB+JUMP
   10 CONTINUE
C
C     REMAINING WAVENUMBERS
      SCALE=0.5*SCALE
      IABASE=3
      IBBASE=N-1
      JABASE=2*INC+1
      JBBASE=(N-2)*INC+1
C
      DO 30 K=3,NH,2
      IA=IABASE
      IB=IBBASE
      JA=JABASE
      JB=JBBASE
      C=TRIGS(N+K)
      S=TRIGS(N+K+1)
CDIR$ IVDEP
      DO 20 L=1,LOT
      A(JA)=SCALE*((WORK(IA)+WORK(IB))
     *   +(C*(WORK(IA+1)+WORK(IB+1))+S*(WORK(IA)-WORK(IB))))
      A(JB)=SCALE*((WORK(IA)+WORK(IB))
     *   -(C*(WORK(IA+1)+WORK(IB+1))+S*(WORK(IA)-WORK(IB))))
      A(JA+INC)=SCALE*((C*(WORK(IA)-WORK(IB))-S*(WORK(IA+1)+WORK(IB+1)))
     *    +(WORK(IB+1)-WORK(IA+1)))
      A(JB+INC)=SCALE*((C*(WORK(IA)-WORK(IB))-S*(WORK(IA+1)+WORK(IB+1)))
     *    -(WORK(IB+1)-WORK(IA+1)))
      IA=IA+NX
      IB=IB+NX
      JA=JA+JUMP
      JB=JB+JUMP
   20 CONTINUE
      IABASE=IABASE+2
      IBBASE=IBBASE-2
      JABASE=JABASE+INK
      JBBASE=JBBASE-INK
   30 CONTINUE
C
      IF (IABASE.NE.IBBASE) GO TO 50
C     WAVENUMBER N/4 (IF IT EXISTS)
      IA=IABASE
      JA=JABASE
      SCALE=2.0*SCALE
CDIR$ IVDEP
      DO 40 L=1,LOT
      A(JA)=SCALE*WORK(IA)
      A(JA+INC)=-SCALE*WORK(IA+1)
      IA=IA+NX
      JA=JA+JUMP
   40 CONTINUE
C
   50 CONTINUE
      RETURN
      END

      SUBROUTINE RFFTMLT(A,WORK,TRIGS,IFAX,INC,JUMP,N,LOT,ISIGN)
      DIMENSION A(N),WORK(N),TRIGS(N),IFAX(1)
C
C     SUBROUTINE "FFT991" - MULTIPLE REAL/HALF-COMPLEX PERIODIC
C     FAST FOURIER TRANSFORM
C
C     SAME AS FFT99 EXCEPT THAT ORDERING OF DATA CORRESPONDS TO
C     THAT IN MRFFT2
C
C     PROCEDURE USED TO CONVERT TO HALF-LENGTH COMPLEX TRANSFORM
C     IS GIVEN BY COOLEY, LEWIS AND WELCH (J. SOUND VIB., VOL. 12
C     (1970), 315-337)
C
C     A IS THE ARRAY CONTAINING INPUT AND OUTPUT DATA
C     WORK IS AN AREA OF SIZE (N+1)*LOT
C     TRIGS IS A PREVIOUSLY PREPARED LIST OF TRIG FUNCTION VALUES
C     IFAX IS A PREVIOUSLY PREPARED LIST OF FACTORS OF N/2
C     INC IS THE INCREMENT WITHIN EACH DATA 'VECTOR'
C         (E.G. INC=1 FOR CONSECUTIVELY STORED DATA)
C     JUMP IS THE INCREMENT BETWEEN THE START OF EACH DATA VECTOR
C     N IS THE LENGTH OF THE DATA VECTORS
C     LOT IS THE NUMBER OF DATA VECTORS
C     ISIGN = +1 FOR TRANSFORM FROM SPECTRAL TO GRIDPOINT
C           = -1 FOR TRANSFORM FROM GRIDPOINT TO SPECTRAL
C
C     ORDERING OF COEFFICIENTS:
C         A(0),B(0),A(1),B(1),A(2),B(2),...,A(N/2),B(N/2)
C         WHERE B(0)=B(N/2)=0; (N+2) LOCATIONS REQUIRED
C
C     ORDERING OF DATA:
C         X(0),X(1),X(2),...,X(N-1)
C
C     VECTORIZATION IS ACHIEVED ON CRAY BY DOING THE TRANSFORMS IN
C     PARALLEL
C
C     *** N.B. N IS ASSUMED TO BE AN EVEN NUMBER
C
C     DEFINITION OF TRANSFORMS:
C     -------------------------
C
C     ISIGN=+1: X(J)=SUM(K=0,...,N-1)(C(K)*EXP(2*I*J*K*PI/N))
C         WHERE C(K)=A(K)+I*B(K) AND C(N-K)=A(K)-I*B(K)
C
C     ISIGN=-1: A(K)=(1/N)*SUM(J=0,...,N-1)(X(J)*COS(2*J*K*PI/N))
C               B(K)=-(1/N)*SUM(J=0,...,N-1)(X(J)*SIN(2*J*K*PI/N))
C
C THE FOLLOWING CALL IS FOR MONITORING LIBRARY USE AT NCAR
C     CALL Q8QST4 ( 4HXLIB, 6HFFT99F, 6HFFT991, 10HVERSION 01)
      NFAX=IFAX(1)
      NX=N+1
      NH=N/2
      INK=INC+INC
      IF (ISIGN.EQ.+1) GO TO 30
C
C     IF NECESSARY, TRANSFER DATA TO WORK AREA
      IGO=50
      IF (MOD(NFAX,2).EQ.1) GOTO 40
      IBASE=1
      JBASE=1
      DO 20 L=1,LOT
      I=IBASE
      J=JBASE
CDIR$ IVDEP
      DO 10 M=1,N
      WORK(J)=A(I)
      I=I+INC
      J=J+1
   10 CONTINUE
      IBASE=IBASE+JUMP
      JBASE=JBASE+NX
   20 CONTINUE
C
      IGO=60
      GO TO 40
C
C     PREPROCESSING (ISIGN=+1)
C     ------------------------
C
   30 CONTINUE
      CALL FFT99A(A,WORK,TRIGS,INC,JUMP,N,LOT)
      IGO=60
C
C     COMPLEX TRANSFORM
C     -----------------
C
   40 CONTINUE
      IA=1
      LA=1
      DO 80 K=1,NFAX
      IF (IGO.EQ.60) GO TO 60
   50 CONTINUE
      CALL VPASSM(A(IA),A(IA+INC),WORK(1),WORK(2),TRIGS,
     *   INK,2,JUMP,NX,LOT,NH,IFAX(K+1),LA)
      IGO=60
      GO TO 70
   60 CONTINUE
      CALL VPASSM(WORK(1),WORK(2),A(IA),A(IA+INC),TRIGS,
     *    2,INK,NX,JUMP,LOT,NH,IFAX(K+1),LA)
      IGO=50
   70 CONTINUE
      LA=LA*IFAX(K+1)
   80 CONTINUE
C
      IF (ISIGN.EQ.-1) GO TO 130
C
C     IF NECESSARY, TRANSFER DATA FROM WORK AREA
      IF (MOD(NFAX,2).EQ.1) GO TO 110
      IBASE=1
      JBASE=1
      DO 100 L=1,LOT
      I=IBASE
      J=JBASE
CDIR$ IVDEP
      DO 90 M=1,N
      A(J)=WORK(I)
      I=I+1
      J=J+INC
   90 CONTINUE
      IBASE=IBASE+NX
      JBASE=JBASE+JUMP
  100 CONTINUE
C
C     FILL IN ZEROS AT END
  110 CONTINUE
      IB=N*INC+1
CDIR$ IVDEP
      DO 120 L=1,LOT
      A(IB)=0.0
      A(IB+INC)=0.0
      IB=IB+JUMP
  120 CONTINUE
      GO TO 140
C
C     POSTPROCESSING (ISIGN=-1):
C     --------------------------
C
  130 CONTINUE
      CALL FFT99B(WORK,A,TRIGS,INC,JUMP,N,LOT)
C
  140 CONTINUE
      RETURN
      END
      SUBROUTINE FFTFAX(N,IFAX,TRIGS)
      DIMENSION IFAX(13)
      REAL      TRIGS(1)
C
C MODE 3 IS USED FOR REAL/HALF-COMPLEX TRANSFORMS.  IT IS POSSIBLE
C TO DO COMPLEX/COMPLEX TRANSFORMS WITH OTHER VALUES OF MODE, BUT
C DOCUMENTATION OF THE DETAILS WERE NOT AVAILABLE WHEN THIS ROUTINE
C WAS WRITTEN.
C
      DATA MODE /3/
      CALL FAX (IFAX, N, MODE)
      I = IFAX(1)
      IF (IFAX(I+1) .GT. 5 .OR. N .LE. 4) IFAX(1) = -99
      IF (IFAX(1) .LE. 0 ) WRITE(6,1) N
      IF (IFAX(1) .LE. 0 ) STOP 999
      CALL FFTRIG (TRIGS, N, MODE)
    1 FORMAT(//5X, ' FFTFAX -- INVALID N =', I5,/)
      RETURN
      END
      SUBROUTINE FAX(IFAX,N,MODE)
      DIMENSION IFAX(10)
      NN=N
      IF (IABS(MODE).EQ.1) GO TO 10
      IF (IABS(MODE).EQ.8) GO TO 10
      NN=N/2
      IF ((NN+NN).EQ.N) GO TO 10
      IFAX(1)=-99
      RETURN
   10 K=1
C     TEST FOR FACTORS OF 4
   20 IF (MOD(NN,4).NE.0) GO TO 30
      K=K+1
      IFAX(K)=4
      NN=NN/4
      IF (NN.EQ.1) GO TO 80
      GO TO 20
C     TEST FOR EXTRA FACTOR OF 2
   30 IF (MOD(NN,2).NE.0) GO TO 40
      K=K+1
      IFAX(K)=2
      NN=NN/2
      IF (NN.EQ.1) GO TO 80
C     TEST FOR FACTORS OF 3
   40 IF (MOD(NN,3).NE.0) GO TO 50
      K=K+1
      IFAX(K)=3
      NN=NN/3
      IF (NN.EQ.1) GO TO 80
      GO TO 40
C     NOW FIND REMAINING FACTORS
   50 L=5
      INC=2
C     INC ALTERNATELY TAKES ON VALUES 2 AND 4
   60 IF (MOD(NN,L).NE.0) GO TO 70
      K=K+1
      IFAX(K)=L
      NN=NN/L
      IF (NN.EQ.1) GO TO 80
      GO TO 60
   70 L=L+INC
      INC=6-INC
      GO TO 60
   80 IFAX(1)=K-1
C     IFAX(1) CONTAINS NUMBER OF FACTORS
      NFAX=IFAX(1)
C     SORT FACTORS INTO ASCENDING ORDER
      IF (NFAX.EQ.1) GO TO 110
      DO 100 II=2,NFAX
      ISTOP=NFAX+2-II
      DO 90 I=2,ISTOP
      IF (IFAX(I+1).GE.IFAX(I)) GO TO 90
      ITEM=IFAX(I)
      IFAX(I)=IFAX(I+1)
      IFAX(I+1)=ITEM
   90 CONTINUE
  100 CONTINUE
  110 CONTINUE
      RETURN
      END
      SUBROUTINE FFTRIG(TRIGS,N,MODE)
      REAL      TRIGS(1)
      PI=2.0*ASIN(1.0)
      IMODE=IABS(MODE)
      NN=N
      IF (IMODE.GT.1.AND.IMODE.LT.6) NN=N/2
      DEL=(PI+PI)/FLOAT(NN)
      L=NN+NN
      DO 10 I=1,L,2
      ANGLE=0.5*FLOAT(I-1)*DEL
      TRIGS(I)=COS(ANGLE)
      TRIGS(I+1)=SIN(ANGLE)
   10 CONTINUE
      IF (IMODE.EQ.1) RETURN
      IF (IMODE.EQ.8) RETURN
      DEL=0.5*DEL
      NH=(NN+1)/2
      L=NH+NH
      LA=NN+NN
      DO 20 I=1,L,2
      ANGLE=0.5*FLOAT(I-1)*DEL
      TRIGS(LA+I)=COS(ANGLE)
      TRIGS(LA+I+1)=SIN(ANGLE)
   20 CONTINUE
      IF (IMODE.LE.3) RETURN
      DEL=0.5*DEL
      LA=LA+NN
      IF (MODE.EQ.5) GO TO 40
      DO 30 I=2,NN
      ANGLE=FLOAT(I-1)*DEL
      TRIGS(LA+I)=2.0*SIN(ANGLE)
   30 CONTINUE
      RETURN
   40 CONTINUE
      DEL=0.5*DEL
      DO 50 I=2,N
      ANGLE=FLOAT(I-1)*DEL
      TRIGS(LA+I)=SIN(ANGLE)
   50 CONTINUE
      RETURN
      END
      SUBROUTINE VPASSM(A,B,C,D,TRIGS,INC1,INC2,INC3,INC4,LOT,N,IFAC,LA)
      DIMENSION A(N),B(N),C(N),D(N)
      REAL      TRIGS(N)
C
C     SUBROUTINE "VPASSM" - MULTIPLE VERSION OF "VPASSA"
C     PERFORMS ONE PASS THROUGH DATA
C     AS PART OF MULTIPLE COMPLEX FFT ROUTINE
C     A IS FIRST REAL INPUT VECTOR
C     B IS FIRST IMAGINARY INPUT VECTOR
C     C IS FIRST REAL OUTPUT VECTOR
C     D IS FIRST IMAGINARY OUTPUT VECTOR
C     TRIGS IS PRECALCULATED TABLE OF SINES " COSINES
C     INC1 IS ADDRESSING INCREMENT FOR A AND B
C     INC2 IS ADDRESSING INCREMENT FOR C AND D
C     INC3 IS ADDRESSING INCREMENT BETWEEN A"S & B"S
C     INC4 IS ADDRESSING INCREMENT BETWEEN C"S & D"S
C     LOT IS THE NUMBER OF VECTORS
C     N IS LENGTH OF VECTORS
C     IFAC IS CURRENT FACTOR OF N
C     LA IS PRODUCT OF PREVIOUS FACTORS
C
      DATA SIN36/0.587785252292473/,COS36/0.809016994374947/,
     *     SIN72/0.951056516295154/,COS72/0.309016994374947/,
     *     SIN60/0.866025403784437/
C
      M=N/IFAC
      IINK=M*INC1
      JINK=LA*INC2
      JUMP=(IFAC-1)*JINK
      IBASE=0
      JBASE=0
      IGO=IFAC-1
      IF (IGO.GT.4) RETURN
      GO TO (10,50,90,130),IGO
C
C     CODING FOR FACTOR 2
C
   10 IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      DO 20 L=1,LA
      I=IBASE
      J=JBASE
CDIR$ IVDEP
      DO 15 IJK=1,LOT
      C(JA+J)=A(IA+I)+A(IB+I)
      D(JA+J)=B(IA+I)+B(IB+I)
      C(JB+J)=A(IA+I)-A(IB+I)
      D(JB+J)=B(IA+I)-B(IB+I)
      I=I+INC3
      J=J+INC4
   15 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
   20 CONTINUE
      IF (LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO 40 K=LA1,M,LA
      KB=K+K-2
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      DO 30 L=1,LA
      I=IBASE
      J=JBASE
CDIR$ IVDEP
      DO 25 IJK=1,LOT
      C(JA+J)=A(IA+I)+A(IB+I)
      D(JA+J)=B(IA+I)+B(IB+I)
      C(JB+J)=C1*(A(IA+I)-A(IB+I))-S1*(B(IA+I)-B(IB+I))
      D(JB+J)=S1*(A(IA+I)-A(IB+I))+C1*(B(IA+I)-B(IB+I))
      I=I+INC3
      J=J+INC4
   25 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
   30 CONTINUE
      JBASE=JBASE+JUMP
   40 CONTINUE
      RETURN
C
C     CODING FOR FACTOR 3
C
   50 IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      IC=IB+IINK
      JC=JB+JINK
      DO 60 L=1,LA
      I=IBASE
      J=JBASE
CDIR$ IVDEP
      DO 55 IJK=1,LOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IC+I))
      C(JB+J)=(A(IA+I)-0.5*(A(IB+I)+A(IC+I)))-(SIN60*(B(IB+I)-B(IC+I)))
      C(JC+J)=(A(IA+I)-0.5*(A(IB+I)+A(IC+I)))+(SIN60*(B(IB+I)-B(IC+I)))
      D(JB+J)=(B(IA+I)-0.5*(B(IB+I)+B(IC+I)))+(SIN60*(A(IB+I)-A(IC+I)))
      D(JC+J)=(B(IA+I)-0.5*(B(IB+I)+B(IC+I)))-(SIN60*(A(IB+I)-A(IC+I)))
      I=I+INC3
      J=J+INC4
   55 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
   60 CONTINUE
      IF (LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO 80 K=LA1,M,LA
      KB=K+K-2
      KC=KB+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      DO 70 L=1,LA
      I=IBASE
      J=JBASE
CDIR$ IVDEP
      DO 65 IJK=1,LOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IC+I))
      C(JB+J)=
     *    C1*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))-(SIN60*(B(IB+I)-B(IC+I))))
     *   -S1*((B(IA+I)-0.5*(B(IB+I)+B(IC+I)))+(SIN60*(A(IB+I)-A(IC+I))))
      D(JB+J)=
     *    S1*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))-(SIN60*(B(IB+I)-B(IC+I))))
     *   +C1*((B(IA+I)-0.5*(B(IB+I)+B(IC+I)))+(SIN60*(A(IB+I)-A(IC+I))))
      C(JC+J)=
     *    C2*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))+(SIN60*(B(IB+I)-B(IC+I))))
     *   -S2*((B(IA+I)-0.5*(B(IB+I)+B(IC+I)))-(SIN60*(A(IB+I)-A(IC+I))))
      D(JC+J)=
     *    S2*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))+(SIN60*(B(IB+I)-B(IC+I))))
     *   +C2*((B(IA+I)-0.5*(B(IB+I)+B(IC+I)))-(SIN60*(A(IB+I)-A(IC+I))))
      I=I+INC3
      J=J+INC4
   65 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
   70 CONTINUE
      JBASE=JBASE+JUMP
   80 CONTINUE
      RETURN
C
C     CODING FOR FACTOR 4
C
   90 IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      IC=IB+IINK
      JC=JB+JINK
      ID=IC+IINK
      JD=JC+JINK
      DO 100 L=1,LA
      I=IBASE
      J=JBASE
CDIR$ IVDEP
      DO 95 IJK=1,LOT
      C(JA+J)=(A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I))
      C(JC+J)=(A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I))
      D(JA+J)=(B(IA+I)+B(IC+I))+(B(IB+I)+B(ID+I))
      D(JC+J)=(B(IA+I)+B(IC+I))-(B(IB+I)+B(ID+I))
      C(JB+J)=(A(IA+I)-A(IC+I))-(B(IB+I)-B(ID+I))
      C(JD+J)=(A(IA+I)-A(IC+I))+(B(IB+I)-B(ID+I))
      D(JB+J)=(B(IA+I)-B(IC+I))+(A(IB+I)-A(ID+I))
      D(JD+J)=(B(IA+I)-B(IC+I))-(A(IB+I)-A(ID+I))
      I=I+INC3
      J=J+INC4
   95 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  100 CONTINUE
      IF (LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO 120 K=LA1,M,LA
      KB=K+K-2
      KC=KB+KB
      KD=KC+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      C3=TRIGS(KD+1)
      S3=TRIGS(KD+2)
      DO 110 L=1,LA
      I=IBASE
      J=JBASE
CDIR$ IVDEP
      DO 105 IJK=1,LOT
      C(JA+J)=(A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I))
      D(JA+J)=(B(IA+I)+B(IC+I))+(B(IB+I)+B(ID+I))
      C(JC+J)=
     *    C2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I)))
     *   -S2*((B(IA+I)+B(IC+I))-(B(IB+I)+B(ID+I)))
      D(JC+J)=
     *    S2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I)))
     *   +C2*((B(IA+I)+B(IC+I))-(B(IB+I)+B(ID+I)))
      C(JB+J)=
     *    C1*((A(IA+I)-A(IC+I))-(B(IB+I)-B(ID+I)))
     *   -S1*((B(IA+I)-B(IC+I))+(A(IB+I)-A(ID+I)))
      D(JB+J)=
     *    S1*((A(IA+I)-A(IC+I))-(B(IB+I)-B(ID+I)))
     *   +C1*((B(IA+I)-B(IC+I))+(A(IB+I)-A(ID+I)))
      C(JD+J)=
     *    C3*((A(IA+I)-A(IC+I))+(B(IB+I)-B(ID+I)))
     *   -S3*((B(IA+I)-B(IC+I))-(A(IB+I)-A(ID+I)))
      D(JD+J)=
     *    S3*((A(IA+I)-A(IC+I))+(B(IB+I)-B(ID+I)))
     *   +C3*((B(IA+I)-B(IC+I))-(A(IB+I)-A(ID+I)))
      I=I+INC3
      J=J+INC4
  105 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  110 CONTINUE
      JBASE=JBASE+JUMP
  120 CONTINUE
      RETURN
C
C     CODING FOR FACTOR 5
C
  130 IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      IC=IB+IINK
      JC=JB+JINK
      ID=IC+IINK
      JD=JC+JINK
      IE=ID+IINK
      JE=JD+JINK
      DO 140 L=1,LA
      I=IBASE
      J=JBASE
CDIR$ IVDEP
      DO 135 IJK=1,LOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I))
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IE+I))+(B(IC+I)+B(ID+I))
      C(JB+J)=(A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))
     *  -(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I)))
      C(JE+J)=(A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))
     *  +(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I)))
      D(JB+J)=(B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))
     *  +(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I)))
      D(JE+J)=(B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))
     *  -(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I)))
      C(JC+J)=(A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))
     *  -(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I)))
      C(JD+J)=(A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))
     *  +(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I)))
      D(JC+J)=(B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))
     *  +(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I)))
      D(JD+J)=(B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))
     *  -(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I)))
      I=I+INC3
      J=J+INC4
  135 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  140 CONTINUE
      IF (LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO 160 K=LA1,M,LA
      KB=K+K-2
      KC=KB+KB
      KD=KC+KB
      KE=KD+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      C3=TRIGS(KD+1)
      S3=TRIGS(KD+2)
      C4=TRIGS(KE+1)
      S4=TRIGS(KE+2)
      DO 150 L=1,LA
      I=IBASE
      J=JBASE
CDIR$ IVDEP
      DO 145 IJK=1,LOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I))
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IE+I))+(B(IC+I)+B(ID+I))
      C(JB+J)=
     *    C1*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))
     *      -(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))
     *   -S1*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))
     *      +(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
      D(JB+J)=
     *    S1*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))
     *      -(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))
     *   +C1*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))
     *      +(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
      C(JE+J)=
     *    C4*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))
     *      +(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))
     *   -S4*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))
     *      -(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
      D(JE+J)=
     *    S4*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))
     *      +(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))
     *   +C4*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))
     *      -(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
      C(JC+J)=
     *    C2*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))
     *      -(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))
     *   -S2*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))
     *      +(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
      D(JC+J)=
     *    S2*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))
     *      -(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))
     *   +C2*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))
     *      +(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
      C(JD+J)=
     *    C3*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))
     *      +(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))
     *   -S3*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))
     *      -(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
      D(JD+J)=
     *    S3*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))
     *      +(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))
     *   +C3*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))
     *      -(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
      I=I+INC3
      J=J+INC4
  145 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  150 CONTINUE
      JBASE=JBASE+JUMP
  160 CONTINUE
      RETURN
      END
      SUBROUTINE HADVECT(U,V,P,IM,JM,IE,IW,DIV,VT1,VT2)
      IMPLICIT NONE

      INTEGER IM, JM

      INTEGER IE(IM*JM)
      INTEGER IW(IM*JM)

      REAL U(IM,JM)
      REAL V(IM,JM)
C Equivalenced to real*8 input variable
      REAL P(IM,JM)

      REAL DIV(IM,JM)
      REAL VT1(IM,JM)
      REAL VT2(IM,JM)

      INTEGER I

      REAL ST1
      REAL ST2



C  CORNER FLUXES

      DO I=1,IM*(JM-1)
       DIV(I,1) = (U(I,1)-U(IW(I),1)) - (V(I,2)-V(I,1))
      ENDDO

      DO I=1,IM*(JM-2)
       VT1(I,2) = 0.5*(P  (I    ,1)+P  (IE(I),2))
     *          *     (DIV(I    ,2)-DIV(IE(I),1))

       VT2(I,2) = 0.5*(P  (I    ,2)+P  (IE(I),1))
     *          *     (DIV(IE(I),2)-DIV(I    ,1))
      ENDDO

      DO I =1,IM
       VT1(I,1 ) = 0.
       VT1(I,JM) = 0.
       VT2(I,1 ) = 0.
       VT2(I,JM) = 0.
      ENDDO

      DO I =1,IM*(JM-1)
       DIV(I,1) = (1./48.)*( (VT1(I    ,2) - VT1(IW(I),1))
     *                     + (VT2(IW(I),2) - VT2(I    ,1)) )
      ENDDO

c NEAR FLUXES

      DO I =1,IM*(JM-3)
       ST1 = 28.*U(I,2) + 3.*(U(IE(I),2)+U(IW(I),2))
     *                  -    (U(I    ,1)+U(I    ,3))
       VT1(I,2) = 0.5*(P(IE(I),2) + P(I,2))*ST1
      ENDDO

      DO I =1,IM
       ST1 = 0.5*(U(IE(I),   1)+U(IW(I),   1)+V(I,   2)-V(IE(I),   2))
       ST2 = 0.5*(U(IE(I),JM-1)+U(IW(I),JM-1)-V(I,JM-1)+V(IE(I),JM-1))

       ST1 = 28.*U(I,   1) + 3.*(U(IE(I),   1)+U(IW(I),   1))
     *                     -    (ST1          +U(I    ,   2))
       ST2 = 28.*U(I,JM-1) + 3.*(U(IE(I),JM-1)+U(IW(I),JM-1))
     *                     -    (ST2          +U(I    ,JM-2))

       VT1(I,   1) = 0.5*(P(IE(I),   1) + P(I,   1))*ST1
       VT1(I,JM-1) = 0.5*(P(IE(I),JM-1) + P(I,JM-1))*ST2
      ENDDO


      DO I =1,IM*(JM-4)
       ST1 = 28.*V(I,3) + 3.*(V(I    ,4)+V(I    ,2))
     *                  -    (V(IE(I),3)+V(IW(I),3))
       VT2(I,3) = 0.5*(P(I,3) + P(I,2))*ST1
      ENDDO

      DO I =1,IM
       ST1 = 26.*V(I,   2) + 3.*(V(I    ,   3)              )
     *                     -    (V(IE(I),   2)+V(IW(I),   2))
       ST2 = 26.*V(I,JM-1) + 3.*(              V(I    ,JM-2))
     *                     -    (V(IE(I),JM-1)+V(IW(I),JM-1))

       VT2(I,   2) = 0.5*(P(I,   2) + P(I,   1))*ST1
       VT2(I,JM-1) = 0.5*(P(I,JM-1) + P(I,JM-2))*ST2

       VT2(I,1 ) = 0.
       VT2(I,JM) = 0.
      ENDDO


      DO I =1,IM*(JM-1)
       DIV(I,1) = DIV(I,1)
     *          + (1./24.)*( (VT1(I,1) - VT1(IW(I),1))
     *                     + (VT2(I,2) - VT2(I    ,1)) )
      ENDDO

C FAR FLUXES

      DO I =1,IM*(JM-1)
       ST1 = U(I,1) + U(IW(I),1)
       VT1(I,1) = 0.5*(P(IE(I),1) + P(IW(I),1))*ST1
      ENDDO

      DO I =1,IM*(JM-3)
       ST1 = V(I,3) + V(I,2)
       VT2(I,2) = 0.5*(P(I,3) + P(I,1))*ST1
      ENDDO

      DO I =1,IM
       VT2(I,JM-1) = 0.
       VT2(I,   1) = 0.
      ENDDO

      DO I =1,IM*(JM-3)
       DIV(I,2) = DIV(I,2) 
     *          - (1./12.)*( (VT1(IE(I),2) - VT1(IW(I),2))
     *                     + (VT2(I    ,3) - VT2(I    ,1)) )
      ENDDO

      DO I =1,IM
       DIV(I,   1) = DIV(I,   1) 
     *          - (1./12.)*( (VT1(IE(I),   1) - VT1(IW(I),   1))
     *                     + (VT2(I    ,   2)                  ) )
       DIV(I,JM-1) = DIV(I,JM-1) 
     *          - (1./12.)*( (VT1(IE(I),JM-1) - VT1(IW(I),JM-1))
     *                     + (                - VT2(I    ,JM-2)) )
      ENDDO



      RETURN
      END
c...  begin p2k code
C*********************************************************************
C******************** FUNCTION    P2K     ***************************
C*********************************************************************
 
      FUNCTION P2K(PL)
      REAL PREF, A1, A2, A3, A4, A5, A6, A7
c ...      PARAMETER (PREF=900.,N=7)
      PARAMETER (PREF=6.3,N=7)
      PARAMETER (A1= 6.99702390562      )
      PARAMETER (A2= 2.2234987077859E-3 )
      PARAMETER (A3=-8.8198782075508E-7 )
      PARAMETER (A4= 5.5989893510156E-10)
      PARAMETER (A5=-4.2210158607379E-13)
      PARAMETER (A6= 3.4837450903957E-16)
      PARAMETER (A7=-3.0411804363195E-19)
 
      REAL T, P2K, PL
 
      T = PL - PREF
      P2K = (T*(T*(T*(T*(T*(T*A7+A6)+A5)+A4)+A3)+A2)+A1)
 
      RETURN
      END
      REAL FUNCTION SSUM(LEN,A,ISTR)
      IMPLICIT NONE
      INTEGER LEN
      INTEGER ISTR
      REAL A(ISTR,LEN)
      INTEGER I
C
      SSUM = 0.
      DO I=1,LEN
      SSUM=SSUM+A(1,I)
      ENDDO
      RETURN
      END
 
      REAL FUNCTION SDOT(LEN,A,ISTRA,B,ISTRB)
      IMPLICIT NONE
      INTEGER LEN
      INTEGER ISTRA,ISTRB
      REAL A(ISTRA,LEN),B(ISTRB,LEN)
      INTEGER I
C
      SDOT = 0.
      DO I=1,LEN
      SDOT=SDOT+A(1,I)*B(1,I)
      ENDDO
      RETURN
      END

*********************************************************************
      SUBROUTINE TRACVERTADV(MAXSLOPE,DT,IM,JM,LM,KM,DSG,PSD,QOM,QOI
     * ,QOI2,PIM,PII)
*     Replaces previous formulation of tracer vertical advection.   *
*     Van-Leer scheme I with slope limitations.                     *
*     Change made by F. Montmessin - Nov. 2004                      *
*                 From Hourdin and Armengaud                        *
*     ( Monthly Wheather Review, vol. 127, 822-837,  1999)          *
*********************************************************************
      IMPLICIT NONE

*     ARGUMENTS
      INTEGER IM,JM,LM,KM
      REAL MAXSLOPE
      REAL DT
      REAL DSG(LM)
      REAL PSD(IM,JM,LM)
      REAL QOM(IM,JM,LM,KM)
      REAL QOI(IM,JM,LM,KM)
      REAL QOI2(IM,JM,LM,KM)
      REAL PIM(IM,JM)
      REAL PII(IM,JM)

*     LOCAL VARIABLES
      INTEGER I,K,L,LL

      REAL P(IM,JM)
      REAL PINV(IM,JM)

      REAL Q(IM,JM,LM,KM)

      REAL DZQW(IM*JM,LM)
      REAL ADZQW(IM*JM,LM)
      REAL DZQ(IM*JM,LM)
      REAL WQ(IM*JM,LM)
      REAL MASSDT(IM*JM,LM)

      REAL DZQMAX
      REAL SIGW
      REAL VERTVEL

      REAL DTINV
      REAL DSGINV(LM)
      REAL DELTAQ

      DTINV = 1. / DT

      DO L=1,LM
        DSGINV(L) = 1. / DSG(L)
      ENDDO

      DO I=1,IM*(JM-1)
        P(I,1)    = PIM(I,1) + PII(I,1) * DT
        PINV(I,1) = 1. / P(I,1)
      ENDDO

      DO L=1,LM
        DO I=1,IM*(JM-1)
          MASSDT(I,L) = DT * DSGINV(L) * PINV(I,1)
        ENDDO
      ENDDO

      DO K=1,KM

      DO L=1,LM
        DO I=1,IM*(JM-1)
          Q(I,1,L,K) = QOM(I,1,L,K) + QOI(I,1,L,K) * DT * PINV(I,1)
        ENDDO
      ENDDO

      DO L=1,LM-1
        DO I=1,IM*(JM-1)
          DZQW(I,L)  = Q(I,1,L+1,K) - Q(I,1,L,K)
          ADZQW(I,L) = ABS(DZQW(I,L))
        ENDDO
      ENDDO

      DO L=2,LM-1
        DO I=1,IM*(JM-1)
          IF (DZQW(I,L)*DZQW(I,L-1).GT.0) THEN
            DZQ(I,L) = 0.5 * (DZQW(I,L) + DZQW(I,L-1))
          ELSE
            DZQ(I,L) = 0.
          ENDIF
          DZQMAX   = MAXSLOPE * MIN(ADZQW(I,L),ADZQW(I,L-1))
          DZQ(I,L) = SIGN( MIN(ABS(DZQ(I,L)),DZQMAX),DZQ(I,L) )
        ENDDO
      ENDDO

      DO I=1,IM*(JM-1)
        DZQ(I,1)  = 0.
        DZQ(I,LM) = 0.
      ENDDO

      DO L=2,LM
        DO I=1,IM*(JM-1)
          VERTVEL = PSD(I,1,L-1)
          IF (VERTVEL.GT.0) THEN
            VERTVEL   = MIN(VERTVEL,DSG(L-1)*P(I,1)*DTINV)
            SIGW      = VERTVEL * MASSDT(I,L-1)
            WQ(I,L-1) = VERTVEL * (Q(I,1,L-1,K) +
     *                             0.5 * (1.-SIGW) * DZQ(I,L-1))
          ELSE
            VERTVEL   = MAX(VERTVEL,-DSG(L)*P(I,1)*DTINV)
            SIGW      = VERTVEL * MASSDT(I,L)
            WQ(I,L-1) = VERTVEL * (Q(I,1,L,K) -
     *                             0.5 * (1.+SIGW) * DZQ(I,L))
          ENDIF
        ENDDO
      ENDDO

      DO I=1,IM*(JM-1)
        WQ(I,LM) = 0.
      ENDDO

      DO L=2,LM
        DO I=1,IM*(JM-1)
          DELTAQ       = + WQ(I,L-1) * DSGINV(L)
     *                   - WQ(I,L)   * DSGINV(L)
          QOI(I,1,L,K) = QOI(I,1,L,K) + DELTAQ
!         QOI2(I,1,L,K)=   WQ(I,L-1) - WQ(I,L)
        ENDDO
      ENDDO

      L = 1
      DO I=1,IM*(JM-1)
        DELTAQ       = - WQ(I,L) * DSGINV(L)
        QOI(I,1,L,K) = QOI(I,1,L,K) + DELTAQ 
!       QOI2(I,1,L,K)=   WQ(I,L-1) - WQ(I,L)
      ENDDO

      ENDDO

      RETURN
      END
*********************************************************************
      SUBROUTINE TRACZONALADV(MAXSLOPE,DT,IM,JM,LM,KM,QOM,QOI,
     *                        PIM,PII,USB,IW,IE,D2PIN)
*     Replaces previous formulation of tracer zonal advection.      *
*     Van-Leer scheme I with slope limitations.                     *
*     Change made by F. Montmessin - Nov. 2004                      *
*                 From Hourdin and Armengaud                        *
*     ( Monthly Weather Review, vol. 127, 822-837,  1999)           *
*								    *
*     Accounts for advection over more than one gridbox in 1        *
*     timestep.  						    *
*     Additional changes made by M. Kahre - Jan. 2007		    *
*		  Also from Hourdin and Armengaud		    *
*								    *
*********************************************************************
      IMPLICIT NONE

*     ARGUMENTS
      INTEGER IM,JM,LM,KM
      REAL MAXSLOPE
      REAL DT
      REAL QOM(IM,JM,LM,KM)
      REAL QOI(IM,JM,LM,KM)
      REAL PIM(IM,JM)
      REAL PII(IM,JM)
      REAL USB(IM,JM,LM)
      REAL D2PIN(IM,JM)

      INTEGER IE(IM*JM)
      INTEGER IW(IM*JM)

*     LOCAL VARIABLES
      INTEGER I,K,L

      REAL P(IM,JM)
      REAL PINV(IM,JM)

      REAL Q(IM,JM,LM,KM)

      REAL DXQU(IM*JM)
      REAL ADXQU(IM*JM)
      REAL DXQ(IM*JM)
      REAL UQ(IM*JM)
      REAL MASSDT(IM*JM)
      REAL MAXMASS(IM*JM)

      REAL DXQMAX
      REAL SIGU
      REAL ZONVEL

      REAL DTINV
      REAL DELTAQ

c!  Mel's Variables
      
      LOGICAL MGRID
      REAL MTOTLEFT
      REAL MTOTRIGHT
      REAL MQTOT , QHAT
      REAL AZONVEL
      REAL UQNEW, QNEW
      INTEGER NCOUNT, JJ
      INTEGER INEW, JNEW, N, NC
      INTEGER IEALL(IM*JM,IM+1)
      INTEGER IWALL(IM*JM,IM+1)
      
      MGRID = .TRUE.		! flag allowing for advection over more
c!      MGRID = .FALSE.         ! than one gridbox in one timestep

      DTINV = 1. / DT

      DO I=1,IM*(JM-1)
        P(I,1)    = PIM(I,1) + PII(I,1) * DT
        PINV(I,1) = 1. / P(I,1)
        MASSDT(I)  = DT    * PINV(I,1) * D2PIN(I,1)
        MAXMASS(I) = DTINV * P(I,1)    / D2PIN(I,1)
c!	if(maxmass(i) .lt. 0.) print*,maxmass(i),p(i,1)
	IEALL(I,1) = I
	IWALL(I,1) = I
        IF(MOD(I,IM) .NE. 0) THEN
         IEALL(I,2) = I + 1
        ELSE
         IEALL(I,2) = I + 1 - IM
        ENDIF
        IF(MOD(I,IM) .NE. 1) THEN
         IWALL(I,2) = I - 1
        ELSE
         IWALL(I,2) = I - 1 + IM
        ENDIF
	DO NC = 3 , IM + 1
	   IF(MOD(IEALL(I,NC-1),IM) .NE. 0) THEN
		IEALL(I,NC) = IEALL(I,NC-1) + 1
	   ELSE
		IEALL(I,NC) = IEALL(I,NC-1) + 1 - IM
	   ENDIF
	   IF(MOD(IWALL(I,NC-1),IM) .NE. 1) THEN
		IWALL(I,NC) = IWALL(I,NC-1) - 1
	   ELSE
		IWALL(I,NC) = IWALL(I,NC-1) - 1 + IM
	   ENDIF
	ENDDO
      ENDDO

      DO K=1,KM
      DO L=1,LM

        DO I=1,IM*(JM-1)
          Q(I,1,L,K) = QOM(I,1,L,K) + QOI(I,1,L,K)*DT*PINV(I,1)
	  Q(I,1,L,K) = MAX(Q(I,1,L,K),0.)
	  UQ(I) = 0.
        ENDDO

        DO I=1,IM*(JM-1)
          DXQU(I)  = Q(IE(I),1,L,K) - Q(I,1,L,K)
          ADXQU(I) = ABS(DXQU(I))
        ENDDO

        DO I=1,IM*(JM-1)
          IF (DXQU(IW(I))*DXQU(I).GT.0) THEN
            DXQ(I) = 0.5 * (DXQU(IW(I)) + DXQU(I))
          ELSE
            DXQ(I) = 0.
          ENDIF
          DXQMAX   = MAXSLOPE * MIN(ADXQU(IW(I)),ADXQU(I))
          DXQ(I)   = SIGN( MIN(ABS(DXQ(I)),DXQMAX),DXQ(I) )
        ENDDO

        DO I=1,IM*(JM-1)
          ZONVEL = USB(I,1,L)
          IF (MGRID .EQV. .FALSE.) THEN         !the old way
           IF (ZONVEL.GT.0.) THEN
            ZONVEL= MIN(ZONVEL,MAXMASS(I))
            SIGU  = ZONVEL * MASSDT(I)
            UQ(I) = ZONVEL * (Q(I,1,L,K) + 0.5 * (1.-SIGU) * DXQ(I))
           ELSE
            ZONVEL= MAX(ZONVEL,-MAXMASS(IE(I)))
            SIGU  = ZONVEL * MASSDT(IE(I))
            UQ(I) = ZONVEL * (Q(IE(I),1,L,K) -
     *                        0.5 * (1.+SIGU) * DXQ(IE(I)))
           ENDIF
          ELSE					! the new way
           IF (ZONVEL.GT.0.) THEN
	      NCOUNT = 0
              IF (ZONVEL.LE.MAXMASS(I)) THEN
                 SIGU = ZONVEL * MASSDT(I)
                 UQ(I) = ZONVEL * (Q(I,1,L,K) + 0.5 *
     *			 (1.-SIGU) * DXQ(I))
              ELSE
                 NCOUNT = 1
1002   CONTINUE
                 MTOTLEFT  = 0.
                 MTOTRIGHT = 0.
                 MQTOT     = 0.
		 NC = NCOUNT + 1
                 DO JJ = I - NCOUNT, I
		    JNEW      = IWALL(I,NC)
                    MTOTRIGHT = MTOTRIGHT + MAXMASS(JNEW)
		    NC        = NC - 1
                 ENDDO
                 IF(ZONVEL.LE.MTOTRIGHT)THEN
		    NC = NCOUNT     
                    DO JJ = I - NCOUNT + 1, I
                       JNEW     = IWALL(I,NC)
                       MTOTLEFT = MTOTLEFT + MAXMASS(JNEW)
                       MQTOT    = MQTOT + MAXMASS(JNEW) * Q(JNEW,1,L,K)
		       NC       = NC - 1
                    ENDDO
                 ELSE
                    NCOUNT = NCOUNT + 1
		    GO TO	 1002
                 ENDIF
	      	 INEW = IWALL(I,NCOUNT+1)
	      	 QHAT = Q(INEW,1,L,K)+.5*(1-(ZONVEL-MTOTLEFT)
     *                    /MAXMASS(INEW))*DXQ(INEW)
              	 UQ(I) = MQTOT + (ZONVEL-MTOTLEFT) * QHAT
              ENDIF	      
           ELSE
              AZONVEL = ABS(ZONVEL)
	      NCOUNT = 0
              IF (AZONVEL.LE.MAXMASS(IE(I))) THEN
                 SIGU  = ZONVEL * MASSDT(IE(I))
                 UQ(I) = ZONVEL * (Q(IE(I),1,L,K) -
     *                        0.5 * (1.+SIGU) * DXQ(IE(I)))
              ELSE
                 NCOUNT = 1
1003   CONTINUE
                 MTOTLEFT  = 0.
                 MTOTRIGHT = 0.
                 MQTOT     = 0.
		 NC = 2
                 DO JJ = I + 1, I + NCOUNT + 1
		    JNEW = IEALL(I,NC)
                    MTOTRIGHT = MTOTRIGHT + MAXMASS(JNEW)
		    NC = NC + 1
                 ENDDO
                 IF(AZONVEL.LE.MTOTRIGHT)THEN
		    NC = 2
                    DO JJ = I + 1, I + NCOUNT
                       JNEW = IEALL(I,NC)
                       MTOTLEFT = MTOTLEFT + MAXMASS(JNEW)
                       MQTOT = MQTOT + MAXMASS(JNEW) * Q(JNEW,1,L,K)
		       NC = NC + 1
                    ENDDO
                 ELSE
                    NCOUNT = NCOUNT + 1
                    GO TO 	1003
                 ENDIF
	         INEW = IEALL(I,NCOUNT + 2)
	         QHAT = Q(INEW,1,L,K) +.5*(1-(AZONVEL-MTOTLEFT)
     *			/MAXMASS(INEW))*DXQ(INEW)
                 UQ(I) = (-1.) * MQTOT - (AZONVEL - MTOTLEFT) * QHAT
              ENDIF
           ENDIF
          ENDIF
        ENDDO

        DO I=1,IM*(JM-1)

c	the following can be added (uncommented out) to ensure that
c       we don't try to remove too much from a gridbox.
c	  IF((MAXMASS(I)*Q(I,1,L,K)+UQ(IW(I))-UQ(I)).LT. 0.) THEN	
c		UQ(I)=MAXMASS(I)*Q(I,1,L,K)+UQ(IW(I))
c	  ENDIF			
		 
          DELTAQ       = ( UQ(IW(I)) - UQ(I) ) * D2PIN(I,1)
          QOI(I,1,L,K) = QOI(I,1,L,K) + DELTAQ
        ENDDO

      ENDDO
      ENDDO

      RETURN
      END
*********************************************************************
      SUBROUTINE TRACMERIDADV(MAXSLOPE,DT,IM,JM,LM,KM,QOM,QOI,
     *                        PIM,PII,VSB,D2PIN)
*     Replaces previous formulation of tracer meridional advection. *
*     Van-Leer scheme I with slope limitations.                     *
*     Change made by F. Montmessin - Nov. 2004                      *
*                 From Hourdin and Armengaud                        *
*     ( Monthly Wheather Review, vol. 127, 822-837,  1999)          *
*********************************************************************
      IMPLICIT NONE

*     ARGUMENTS
      INTEGER IM,JM,LM,KM
      REAL MAXSLOPE
      REAL DT
      REAL QOM(IM,JM,LM,KM)
      REAL QOI(IM,JM,LM,KM)
      REAL PIM(IM,JM)
      REAL PII(IM,JM)
      REAL VSB(IM,JM,LM)
      REAL D2PIN(IM,JM)

*     LOCAL VARIABLES
      INTEGER I,K,L

      REAL P(IM,JM)
      REAL PINV(IM,JM)

      REAL Q(IM,JM,LM,KM)

      REAL DYQV(IM,JM)
      REAL ADYQV(IM,JM)
      REAL DYQ(IM,JM)
      REAL VQ(IM,JM)
      REAL MASSDT(IM,JM)

      REAL DYQMAX
      REAL SIGV
      REAL MERVEL

      REAL DTINV
      REAL DELTAQ

      DTINV = 1. / DT

      DO I=1,IM*(JM-1)
        P(I,1)      = PIM(I,1) + PII(I,1) * DT
        PINV(I,1)   = 1. / P(I,1)
        MASSDT(I,1) = DT * PINV(I,1) * D2PIN(I,1)
      ENDDO

      DO K=1,KM
      DO L=1,LM

        DO I=1,IM*(JM-1)
          Q(I,1,L,K) = QOM(I,1,L,K) + QOI(I,1,L,K) * DT * PINV(I,1)
        ENDDO

        DO I=1,IM*(JM-2)
          DYQV(I,2)  = Q(I,2,L,K) - Q(I,1,L,K)
          ADYQV(I,2) = ABS(DYQV(I,2))
        ENDDO

        DO I=IM+1,IM*(JM-2)
          IF (DYQV(I,1)*DYQV(I,2).GT.0) THEN
            DYQ(I,1) = 0.5 * (DYQV(I,1) + DYQV(I,2))
          ELSE
            DYQ(I,1) = 0.
          ENDIF
          DYQMAX   = MAXSLOPE * MIN(ADYQV(I,1),ADYQV(I,2))
          DYQ(I,1) = SIGN( MIN(ABS(DYQ(I,1)),DYQMAX),DYQ(I,1) )
        ENDDO

        DO I=1,IM
          DYQ(I,JM-1)  = 0.
          DYQ(I,1)     = 0.
        ENDDO

        DO I=1,IM*(JM-2)
          MERVEL = VSB(I,2,L)
          IF (MERVEL.GT.0) THEN
            SIGV    = MERVEL * MASSDT(I,1)
            VQ(I,2) = MERVEL * (Q(I,1,L,K) + 0.5 * (1.-SIGV) * DYQ(I,1))
          ELSE
            SIGV    = MERVEL * MASSDT(I,2)
            VQ(I,2) = MERVEL * (Q(I,2,L,K) - 0.5 * (1.+SIGV) * DYQ(I,2))
          ENDIF
        ENDDO

        DO I=1,IM
          VQ(I,1) = 0. 
        ENDDO

        DO I=1,IM*(JM-2)
          DELTAQ       = ( VQ(I,1) - VQ(I,2) ) * D2PIN(I,1)
          QOI(I,1,L,K) = QOI(I,1,L,K) + DELTAQ
        ENDDO

        DO I=1,IM
          DELTAQ          = VQ(I,JM-1) * D2PIN(I,JM-1)
          QOI(I,JM-1,L,K) = QOI(I,JM-1,L,K) + DELTAQ
        ENDDO

      ENDDO
      ENDDO

      RETURN
      END
