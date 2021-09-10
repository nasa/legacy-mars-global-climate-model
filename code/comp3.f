      SUBROUTINE COMP3

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

!
!  PURPOSE:
!      COMP3 COMPUTES THE CHANGES IN GROUND AND ATMOSPHERIC TEMPERATURES
!      CAUSED BY RADIATIVE (SOLAR AND INFRARED) AND NON-RADIATIVE
!      (CONVECTION, TURBULENCE, ETC.) EFFECTS. CHANGES IN WIND VELOCITY
!      CAUSED BY FRICTION ARE ALSO COMPUTED. ( SEE UNDER OUTPUT
!      PARAMETERS FOR OTHER QUANTITIES COMPUTED. )
!  AUTHOR
!      STEVE POHORSKY    INFORMATICS     TASK 605    OCT 81
!  FOR
!      JIM POLLACK
!  ENVIRONMENT
!      Cray-2            UNICOS 3.0      FORTRAN
!  REVISION HISTORY
!     6-82  SP   CALL DUSTIR & CALL CLOUDIR REPLACED WITH CALL OUT15IR.
!     7-82  SP   CO2IR SUBROUTINE NAME CHANGED TO IN15IR.
!     7-82  SP   IRTOTAL CALCULATION REVISED FOR ABOVE CHANGES.
!     7-82  SP   EMISSIVITY FACTORS INCLUDED IN FA CALCULATION AS IN
!                NOTE OF 4-28-82.
!     7-82  SP   TAUICET CALCULATION ADDED.
!     9-82  SP   IRNET CALCULATION ADDED.
!    11-82  SP   TETA(K) SET UP FOR ODD K AS WELL BEFORE CALL TO CONVECT
!
!    JR SCHAEFFER     TASK 904           MAR 86
!    STRUCTURE THE CODE AND INCORPORATE THE TIME-AVERAGED
!    SOLAR ZENITH ANGLE (SEE NOTE OF 16-SEP-85).
!
!    JR SCHAEFFER     TASK 904           DEC 86
!    MODIFY THE TIME-AVERAGED SOLAR ZENITH ANGLE (SEE NOTE OF 9/15/86).
!    INCLUDE RAYLEIGH FRICTION (SEE NOTE OF 9/15/86).
!    MISCELLANEOUS MODIFICATIONS TO THE CODE STRUCTURE (9/20/86).
!    J SCHAEFFER      TASK 904          FEB 87
!    ADDED NEW VARIABLES TO THE HISTORY TAPE (BOB'S NOTE OF 2/18/87).
!    J. SCHAEFFER     TASK 904          JULY/AUGUST 87
!    HISTORY TAPE STRUCTURE COMPLETELY REDONE, WITH NEW VARIABLES
!    COMPUTED IN COMP3.  SEE NOTES OF 6/29/87 AND 7/31/87.
!  INPUT PARAMETERS
!      NSTEP          - NUMBER OF THE CURRENT TIME STEP. THIS SERVES AS
!                       A CONTROL PARAMETER IN THAT COMP3 PERFORMS THE
!                       ATMOSPHERIC TEMPERATURE CALCULATIONS ONLY IF
!                       NSTEP IS A MULTIPLE OF NC3 (= A CONSTANT READ IN
!                       AS A RUN PARAMETER).
!      GT(J,I) ARRAY  - GROUND TEMPERATURE OR (-1) TIMES CO2 ICE
!                             MASS AT THE 'PI' POINTS.
!      P(J,I) ARRAY   - THE CURRENT 'PI' VALUE (SURFACE PRESSURE
!                             MINUS TROPOPAUSE PRESSURE) AT THE 'PI'
!                             POINTS.
!      T(J,I,K) ARRAY - THE TEMPERATURE AT THE MIDPOINT OF EACH
!                             LAYER FOR THE 'PI' POINTS. (K = 1 TO NLAY)
!      U(J,I,K) ARRAY       - THE LOCAL EAST-WEST WIND VELOCITY
!                             COMPONENT AT EACH OF THE U,V-GRID POINTS
!                             AT THE MIDPOINTS OF EACH LAYER.
!      V(J,I,K) ARRAY       - THE LOCAL NORTH-SOUTH WIND VELOCITY
!                             COMPONENT AT EACH OF THE U,V-GRID POINTS
!                             AT THE MIDPOINTS OF EACH LAYER.
!      ALSP(J,I) ARRAY- THE LOCAL SURFACE ALBEDO FOR THE 'PI' POINTS
!                             ASSUMING NO CO2 ICE IS PRESENT.
!      QDUST(J,I,K)   - THE MIXING RATIO OF DUST TO AIR IN THE KTH
!        ARRAY                LAYER. (K = 1 TO NLAY.)
!      TINF(J,I) ARRAY- THE DEEP SUBSURFACE TEMPERATURES FOR THE 'PI'
!                       POINTS.
!      DMADT(J,I)     - THE RATE OF CONDENSATION PER UNIT SURFACE AREA
!       ARRAY           IN THE WHOLE ATMOSPHERIC COLUMN AT A 'PI' GRID
!                       POINT.
!      FA(J,I) ARRAY  - THE RATE AT WHICH HEAT IS BEING GAINED BY THE
!                       SURFACE IN EXCHANGE WITH THE ATMOSPHERE
!                       (THROUGH THERMAL RADIATION AND CONVECTION)
!                       FOR EACH 'PI' POINT. (THE FA VALUES COMPUTED IN
!                       THE ATMOSPHERIC TEMPERATURE CALCULATIONS
!                       ARE USED AS INPUT IN THE NEXT TIME STEPS UNTIL
!                       FA IS CALCULATED AGAIN.)
!      JBPS           - APPROXIMATE LATITUDE INDEX OF EDGE OF SOUTHERN
!                       POLAR CAP.
!      JBPN           - APPROXIMATE LATITUDE INDEX OF EDGE OF NORTHERN
!                       POLAR CAP.
!  OUTPUT PARAMETERS
!      GT(J,I) ARRAY  - AS IN INPUT, BUT ADJUSTED ACCORDING TO
!                       HEATING AT THE GROUND.
!      TINF(J,I) ARRAY- AS IN INPUT, BUT ADJUSTED IN PROPORTION TO
!                       TO THE NEW GROUND TEMPERATURE.
!      SDGR(J,I) ARRAY- THE VERTICAL (SIGMA) COMPONENT OF WIND VELOCITY
!                       AT THE GROUND CALCULATED FROM THE CHANGE IN THE
!                       AMOUNT OF CO2 ICE ON THE GROUND.
!    THE REST OF THE OUTPUT PARAMETERS ARE PRODUCED ONLY WHEN THE
!    ATMOSPHERIC TEMPERATURE CALCULATIONS ARE PERFORMED.
!      T(J,I,K) ARRAY - AS IN INPUT, BUT ADJUSTED ACCORDING TO
!                       RADIATIVE & NON-RADIATIVE EFFECTS.
!      U(J,I,K) ARRAY - AS IN INPUT, BUT ADJUSTED ACCORDING TO FRICTION.
!      V(J,I,K) ARRAY - AS IN INPUT, BUT ADJUSTED ACCORDING TO FRICTION.
!      FA(J,I) ARRAY  - AS IN INPUT.
!  CALLED BY
!      MAIN
!  SUBROUTINES CALLED
!      TEMPGR,CMP3SET,GRIDVEL,POTEMP,CO2SUN,DUSTSUN,IN15IR,OUT15IR,
!      CONVECT,TURBO,RICHNUM,CMP3OUT,C3OUT2,GEOPOTN,EFFLUX,COLDAIR
!
      use grid_h
      use defines_h
      use constants_h, only: PI, GRAV, CP, SCALEP, XLHTC
      use radinc_h 
      use radcommon_h
      use fccsave_h
      use standard_h
      use comp3cmn_h
      use dtcommon_h
      use cldcommon_h

      implicit none

!######################################################################

      real*8  :: COSPTG(L_JSIZE),THETPTG(L_JSIZE)
      real*8  :: TETAORG(L_LEVELS)
      real*8  :: TETACON(L_LEVELS), UPICON(L_LEVELS)
      real*8  :: VPICON(L_LEVELS)

      REAL*8  :: Z(2*L_LAYERS+1)
      REAL*8  :: KD(2*L_LAYERS+1)

! Tracer mixing in vertical (only in CONVECT at the moment)

      REAL*8  :: test1(ntrace),test2(ntrace)

      REAL*8  :: QPI(L_LEVELS,NTRACE),QPISAV(L_LEVELS,NTRACE)
      real*8  :: tetanew(l_levels),deltateta(l_levels)
      real*8  :: qnew(l_levels,ntrace),deltawv(l_levels)

! Amount of Tracer Deposit on the ground (kg/m2): QPIG

      REAL*8  :: QPIG(NTRACE)

!      LOGICAL :: CLOUDON
!      LOGICAL :: CO2SCAV

      REAL*8  :: RHOUCHpbl

      REAL*8  :: QRAD(L_LEVELS)

!  Flag to turn on of off surface source of dust (Surface Source ON)

      LOGICAL :: SSON

!  Radiation code variables

      REAL*8 :: PLEV(L_LEVELS), TLEV(L_LEVELS)
      REAL*8 :: PMID(L_LEVELS), TMID(L_LEVELS)
      REAL*8 :: SOL(L_NSPECTV), SFMARS

!  VISIBLE

      real*8 :: DTAUV(L_NLAYRAD,L_NSPECTV,L_NGAUSS)
      real*8 :: TAUV(L_NLEVRAD,L_NSPECTV,L_NGAUSS)
      real*8 :: TAUCUMV(L_LEVELS,L_NSPECTV,L_NGAUSS)
      real*8 :: COSBV(L_NLAYRAD,L_NSPECTV,L_NGAUSS)
      real*8 :: WBARV(L_NLAYRAD,L_NSPECTV,L_NGAUSS)
      real*8 :: FMNETV(L_NLAYRAD)
      real*8 :: fluxupv(L_NLAYRAD), fluxdnv(L_NLAYRAD), NFLUXTOPV
      real*8 :: taugsurf(L_NSPECTV,L_NGAUSS-1)
      real*8 :: dnvflux(L_J,L_I)

!  IR

      real*8 :: DTAUI(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
      real*8 :: TAUCUMI(L_LEVELS,L_NSPECTI,L_NGAUSS)
      real*8 :: COSBI(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
      real*8 :: WBARI(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
      real*8 :: FMNETI(L_NLAYRAD)
      real*8 :: fluxupi(L_NLAYRAD), fluxdni(L_NLAYRAD), NFLUXTOPI
      real*8 :: taugsurfi(L_NSPECTI,L_NGAUSS-1)

!  Water mixing 

      real*8  :: QH2O(L_LEVELS)

      real*8  :: longit,loctime
      real*8  :: ecart

      integer :: univtime

!     implicit none

      integer :: I, J, K, L, M, MM, NN, ILON, NG, NW
      integer :: L_SCAVUP, L_SCAVDN, NS
      real*8  :: x3, c6, c7, t1, t2, albi, diffvt
      real*8  :: directsol, rinc3, xpsat, xtsat, slope
      real*8  :: atmass, rkef, rkei, pie, sunlte, pnlte, dumz0, ht_pbl
      real*8  :: sx_pbl, sy_pbl, dustref
      real*8  :: taucrt = 0.0
      real*8  :: h2oflux

!#=====================================================================

!  Total water mass: atm, cld, surface

        if(mod(nstep,nc3*5).eq.0) then
          write(6,'("------------------- TIME: ",f10.2,5x,"Ls: ",f6.2)')
     *                                     tau, vpout
          call watermass
        endif


      if(TIMESPLIT .eqv. .FALSE.) then
        nsplit = 1
        sdt = dt*float(nc3)/float(nsplit)
      else
        NSPLIT = INT(DT*FLOAT(NC3)/DTSPLIT)
        SDT = DT*FLOAT(NC3)/FLOAT(NSPLIT)
      endif

!      icealb=0.4

!      icethresh_depth=5.0         !microns of ice required to change albedo

!      icethresh_kgm2=icethresh_depth*dpden_ice*1.e-6

      if (firstcall) then
!******************************************************************************
!    Find the universal time corresponding to 2PM (TES obs) for each longitude  
!******************************************************************************

       print*,sdt,nsplit

        if (mod(88775.0D0,nc3*dt).lt.0.5) then
          NIT = int(88775./float(nc3)/dt)
        else
          NIT = int(88775./float(nc3)/dt) + 1
        endif
        DO ILON=1,L_I
           ecart=24.
           longit = -180.+360./L_I*(ILON-1)
           do i=1,NIT
             loctime = i/float(NIT)*24. + longit*12./180.
             if (loctime.gt.24.) loctime = loctime - 24.
             if (loctime.lt.0)   loctime = loctime + 24.
             if (abs(14.-loctime).lt.ecart) then
               ecart      = abs(14.-loctime)
               LON2PM(ilon) = i
             endif
           enddo
      
        ENDDO
      
        call settozero4(L_I*L_J,LATHEAT)
        call settozero4(L_I*L_J*L_LAYERS,ATMCOND)
        call ini_optdst(QEXTV,QSCATV,GV,QEXTI,QSCATI,GI,  
     *                  QXVDST,QXIDST,QSVDST,QSIDST,GVDST,GIDST,
     *                  QEXTREFDST)
        call ini_optcld(QXVCLD,QXICLD,QSVCLD,QSICLD,GVCLD,GICLD,
     *                  QEXTREFCLD,TAUREFCLD)
        call firstcomp3(L_J,L_I,L_NSPECTV,L_NGAUSS,detau,dnirflux,
     *                  dndiffv)
        firstcall = .FALSE.
      endif

!  Spatial dust pattern is fixed through spin-up, then allowed to
!  change via the dust tracer scheme.

!      CLOUDON      = .TRUE.
!      ACTIVE_DUST  = .TRUE.
!      CO2SCAV      = .TRUE.
!      ACTIVE_WATER = .FALSE.
!      MICROPHYSICS = .TRUE.
!      TIMESPLIT    = .TRUE.
!      ALBFEED      = .FALSE.


!     If microphysics is turned on, the above flags can have any set of
!     values. If microphysics is turned off, the flags must be set as 
!     follows:

      if(MICROPHYSICS .eqv. .FALSE.) then
        CLOUDON      = .FALSE.
        ACTIVE_DUST  = .FALSE.
        CO2SCAV      = .FALSE.
        ACTIVE_WATER = .FALSE.
        TIMESPLIT    = .FALSE.
        ALBFEED      = .FALSE.
      end if

      if(TIMESPLIT .eqv. .FALSE.) then
        nsplit = 1
        sdt = dt*float(nc3)/float(nsplit)
      endif
 
!      print*,albfeed,cloudon,active_dust,co2scav,active_water,
!     *       microphysics,timesplit,albfeed,latent_heat,vdust
!      print*,'dtsplit ',dtsplit,nsplit,sdt

!  Calculate solar flux at the current mars distance
!  SFMARS: SOLAR FLUX AT MARS

!  Early mars 1-D  SOL(NW)*0.24986*0.75

      sfmars = 0.0
      do NW=1,L_NSPECTV
        SOL(nw) = SOLARF(NW)/RSDIST
        SFMARS  = SFMARS + SOL(NW)
      end do

!     write(6,'("Solar flux at Mars: ",f7.2)') sfmars

!     GET RADIATION FOR TEMPGR:
!     LOOP OVER ALL I AND J, THAT IS, OVER ALL 'PI' POINTS.

      DO I=1,IM
        DO J=1,JM-1

!         SET UP VARIOUS VARIABLES FOR THIS TIME THROUGH THE LOOP.

          ICMN = I
          JCMN = J

          X3    = SIND*SINL(J)
          ACOSZ = SIND*SINL(J)+COSD*COSL(J)*
     *              COS(2.*PI*(TOFDAY/ROTPER+FLOAT(I-1)/FLOAT(IM)))

          if(ACOSZ.lt.1.0D-5) ACOSZ = 0.0D0

!         AMOUNT OF SOLAR ENERGY ENTERING AT THE TOP OF THE
!         ATMOSPHERE.  EQUAL TO 0 IF THE SUN IS DOWN, OR VERY LOW.
!         SCOSZ = 1.353E3*ACOSZ/RSDIST

!          RADIATIVE CALCULATIONS.

!  Set up, and solve for, the solar (visible) fluxes, if the sun
!  is up


          if(acosz.ge.1.0e-5) then

            call dsolflux(SOL,ACOSZ,GWEIGHT,FZEROV,DETAU,JCMN,ICMN,
     *                    DIRECTSOL)

          else

!  If the sun is down, no solar flux, nor downward flux. . .

            DIRECTSOL = 0.0

          end if    ! condition on sun elevation 

          DNVFLUX(J,I)  = DNDIFFV(J,I) + DIRECTSOL

        END DO         ! J-loop
      END DO           ! I-loop

!     WANTIT IS TRUE IF PRINTER OUTPUT DESIRED FOR THIS TIME STEP.

      WANTIT = (MPRINT .NE. 0)

      IF(WANTIT) THEN
        WANTIT = NSTEP.EQ.NC3 .OR. MOD(NSTEP,NC3*MPRINT).EQ.0
      END IF

!     GROUND TEMPERATURE CALCULATION.

      CALL TEMPGR(GT,CO2ICE,IM,JM,PTROP,P,T,ASYM,TAUCRT,ALSP,FA,EG15GND,
     *            EGOGND,SQRDY,ZIN,DT,ROT,DMADT,SDGR,TINF,STEMP,
     *            STHICK,SCOND,DNIRFLUX,DNVFLUX,RHOUCH,QCOND,
     *            NPCFLAG,LATHEAT,CPSOIL,RHOSOIL,QTRACE)

!     THE FULL COMP3 IS DONE ONLY EVERY NC3 TIME STEPS.

      IF((NSTEP.eq.1) .or. (MOD(NSTEP,NC3) .EQ. 0)) THEN

        NSTEPC3 = NSTEPC3 + 1

!       ATMOSPHERIC TEMPERATURE CALCULATIONS PLUS SOME FRICTION
!       CALCULATIONS.  THIS IS THE MAIN LOOP OF COMP3.

!       LOOP INITIALIZATION.

        CALL CMP3SET(IM,NC3,DT,IMOVER2,IMOVER4,NDT,UT,VT,TT,
     *               QTDELTA,FRY, HSOLCO2, HSOLDST, HTH15, HTHOUT,
     *               HCONADJ, HTURBO, HRINUM, FCONADJ, FTURBO,
     *               FRINUM, FRAYFR, HRAYFR, DISRAY,CO2LAT,QCDEL)

        call settozero(L_I*L_J*NTRACE,SRFUPFLX)
        call settozero(L_I*L_J*NTRACE,SRFDNFLX)
        call settozero(l_i*l_j*l_levels,tauref3d)

!       THE FOLLOWING (THETPTG & COSPTG) PERTAIN TO GRAPHICS.

        C6 = 2.0*PI*RAD
        C7 = (XLHTC/GRAV)*C6

!  Only go from 1 to JM-1 in the c-grid
!  180/24 = 7.5, L_JSIZE = 24, so we want that, instead of L_JSIZEM1 (=23)

        DO 100 J=1,L_JSIZE-1
          THETPTG(J) = -90.0+180.0*(FLOAT(J)-1.0)/FLOAT(L_JSIZE)
          COSPTG(J)  = COS(PI*THETPTG(J)/180.0)
  100   CONTINUE

        if (mod(nstep,nc3*5).eq.0) then
          DO M=1,NTRACE
            test1(M)=0.
            test2(M)=0.
          ENDDO
        endif

!       MAIN COMPUTATIONAL LOOP.
!       LOOP OVER ALL I AND J, THAT IS, OVER ALL 'PI' POINTS.

        DO 500 I=1,IM

!  Only go from 1 to JM-1 in the c-grid

            DO 490 J=1,JM-1
      
!              SET UP VARIOUS VARIABLES FOR THIS TIME THROUGH THE LOOP.

               ICMN = I
               JCMN = J
               IM1  = MOD(I+IM-2,IM) + 1
               IP1  = MOD(I,IM) + 1
               JP1  = J + 1
               JM1  = J - 1

!              Calculate the average cosine of the solar zenith angle -
!              ACOSZ - for this grid point, at this local time.

               call solarza(I,IM,sind,cosd,sinl(j),cosl(j),tofday,day,
     *                      rotper,ndt,dt,acosz)

!              CALCULATE CO2 CONDENSATION TEMPERATURE AT THE SURFACE
!              FOR EACH GRID POINT.  IF THE GROUND TEMPERATURE IS LOWER
!              THAN THE CO2 FROST POINT, SET THE GROUND TEMPERATURE
!              EQUAL TO THE FROST POINT.

               PSAT = P(J,I) + PTROP
               TSAT = 3182.48/(23.3494-LOG(PSAT))

               IF(GT(J,I).LT.TSAT .or. CO2ICE(J,I).gt.0.0) THEN
                 GT(J,I) = TSAT
               END IF

!              CALCULATE THE ICE CLOUD OPTICAL DEPTH, WHERE PRESENT.

               DO 120 K=1,NLEVELS,2
                  TAUICET(K) = 0.0
                  TAUICEL(K) = 0.0
  120          CONTINUE

!              SET UP WIND VELOCITY VALUES AT THIS GRID POINT.

               CALL GRIDVEL(JCMN,ICMN,U,V,QTRACE,QCOND,UPI,VPI,QPI,QPIG)

!              SET UP POTENTIAL TEMPERATURE AND RELATED VALUES AT
!              THIS GRID POINT.

               CALL POTEMP1(JCMN,ICMN,PTROP,PSL,P,SIGMA,TSTRAT,T,
     *                      AADJ,BADJ,PLOGADJ,PL,OM,TL,TETA)

!              STORE INITIAL VALUES FOR TETA, UPI, AND VPI.  ALL AT
!              'PI' POINTS.

               TETAORG(2) = TETA(2)

               DO 140 K=4,NLEVSM1,2
                 TETAORG(K) = TETA(K)
! tracer mixing
                 DO M = 1,NTRACE
                   QPISAV(K,M)  = QPI(K,M)
                 END DO

                 UPISAV(K)  = UPI(K)
                 VPISAV(K)  = VPI(K)
  140          CONTINUE

!              SET UP MASS PER UNIT AREA VALUES FOR EACH LAYER.
!              PRESSURES CONVERTED TO MKS BY SCALEP( = 100).

               YM(2) = PTROP*SCALEP/GRAV

               DO 150 K = 4, NLEVSM1, 2
                 YM(K) = (PL(K+1)-PL(K-1))*SCALEP/GRAV
  150          CONTINUE

!
!  call coldair here - gcm1.7.3_b version
!
               CALL COLDAIR(PTROP,TSTRAT,NDT,SIGMA,PL,TL,DSIG,YM,IM,GT,
     *                      CO2ICE,SQRDY,ZIN,CO2LATST,
     *                      CO2LAT,DMADT,ATMCOND,JCMN,ICMN) 

               CALL POTEMP2(JCMN,ICMN,PTROP,PSL,P,SIGMA,TSTRAT,
     *                      AADJ,BADJ,PLOGADJ,PL,OM,TL,TETA)


!              RADIATIVE CALCULATIONS.

!  Fill the new radiation code variables.
!  PLEV and TLEV are the pressure and tempertures on a vertical grid
!  that the new radiation code uses.

               CALL FILLPT(PL,P(J,I),PTROP,GT(J,I),TSTRAT(J,I),TL,
     *                     PLEV,TLEV,PMID,TMID)

! AEROSOL SCAVENGING BY CO2 SNOW FALL
               IF (CO2SCAV) THEN

               L_SCAVUP = 0
               L_SCAVDN = 0
c              Determine the portion of atmosphere affected by CO2 condensation
               DO L=1,L_LAYERS
                 IF (ATMCOND(JCMN,ICMN,L).GT.0) THEN
                   L_SCAVUP = L
                   EXIT
                 ENDIF
               ENDDO
               DO L=L_LAYERS,1,-1
                 IF (ATMCOND(JCMN,ICMN,L).GT.0) THEN
                   L_SCAVDN = L
                   EXIT
                 ENDIF
               ENDDO

               IF (L_SCAVUP.NE.0.AND.L_SCAVDN.NE.0) THEN
               DO L=L_SCAVUP,L_SCAVDN
!              If condensation occurs down to the surface, put all
!              aerosols on the surface (in fact, only the cloud mass
!              matters...qpig(iMa_vap) changes accordingly)).
                 IF (L_SCAVDN.EQ.L_LAYERS) THEN
                   QPIG(iMa_vap) = QPIG(iMa_vap) + SCAVEFF *
     *                             QPI(2*L+2,iMa_cld) * SCALEP *
     *                             (PL(2*L+3)-PL(2*L+1)) / GRAV
                   qpig(iMa_dt) = qpig(iMa_dt) + SCAVEFF *
     *                             qpi(2*L+2,iMa_dt) * SCALEP *
     *                             (PL(2*L+3)-PL(2*L+1)) / GRAV
                   qpig(iMa_cor) = qpig(iMa_cor) + SCAVEFF *
     *                              qpi(2*L+2,iMa_cor) * SCALEP *
     *                              (PL(2*L+3)-PL(2*L+1)) / GRAV
                   srfdnflx(jcmn,icmn,ima_dt) =
     *                          srfdnflx(jcmn,icmn,ima_dt) +
     *                             (1./ndt) * SCAVEFF *
     *                             qpi(2*L+2,iMa_dt) * SCALEP *
     *                             (PL(2*L+3)-PL(2*L+1)) / GRAV
                   srfdnflx(jcmn,icmn,ima_cor) =
     *                          srfdnflx(jcmn,icmn,ima_cor) +
     *                             (1./ndt) * SCAVEFF *
     *                             qpi(2*L+2,iMa_dt) * SCALEP *
     *                             (PL(2*L+3)-PL(2*L+1)) / GRAV
                   srfdnflx(jcmn,icmn,ima_cld) =
     *                          srfdnflx(jcmn,icmn,ima_cld) +
     *                             (1./ndt) * SCAVEFF *
     *                             qpi(2*L+2,iMa_dt) * SCALEP *
     *                             (PL(2*L+3)-PL(2*L+1)) / GRAV

                 ELSE
!              If condensation occurs in a restricted portion, put
!              aerosols in the highest layer unaffected by CO2 condensation.
                   DO M = 1,NAER
                     QPI(2*L_SCAVDN+4,M) = QPI(2*L_SCAVDN+4,M) + SCAVEFF
     *                                    *QPI(2*L+2,M) 
     *                                    *(PL(2*L+3)-PL(2*L+1))
     *                             /(PL(2*L_SCAVDN+5)-PL(2*L_SCAVDN+3))
                   ENDDO
                 ENDIF
                 DO M = 1,NAER
                   QPI(2*L+2,M) = QPI(2*L+2,M) * (1.-SCAVEFF)
                 ENDDO
                 ATMCOND(JCMN,ICMN,L) = 0.
               ENDDO
               ENDIF

               ENDIF   ! End condition on CO2SCAV

!  Fill QPI with water information
               if(ACTIVE_WATER .eqv. .TRUE.) then
                 M = iMa_vap
                 DO  L = 1, L_LAYERS
                   K = 2*L+2
                   QH2O(K)   = mwratio*QTRACE(JCMN,ICMN,L,M)
                   QH2O(K+1) = QH2O(K)
                 END DO
               else
                 DO  L = 1, L_LAYERS
                   K = 2*L+2
                   QH2O(K)   = 1.0E-7
                   QH2O(K+1) = QH2O(K)
                 END DO
               end if

!  Fill the TAUREF array, the dust column density for each GCM sub-layer
               IF(ACTIVE_DUST .EQV. .FALSE.) then
                 call filltaucum(J,I)
               ELSE

!  No dust in the stratosphere

              CALL opt_dst(JCMN,ICMN,QTRACE,PL,
     *                     QXVDST,QXIDST,QSVDST,QSIDST,GVDST,GIDST,
     *                     QEXTREFDST,TAUREF)

               DO K = 1,3
                 TAUREF(K) = 0.0
                 TAUCUM(K) = 0.0
               END DO

               DO L=1,L_LAYERS
                 DO NN=1,2
                   K=2*L+1+NN
                   TAUCUM(K) = TAUCUM(K-1) + TAUREF(K)
                 END DO
               END DO

               END IF

!  Fill special bottom radiation level to zero.

                TAUREF(L_LEVELS+1) = 0.0
                TAUSURF(J,I)       = TAUCUM(L_LEVELS)

                dustref = TAUTOTJI(J,I)
                IF (CLOUDON) THEN
                  call opt_cld(JCMN,ICMN,QTRACE,PL,
     *                         QXVCLD, QXICLD,QSVCLD,QSICLD,GVCLD,GICLD,
     *                         QEXTREFCLD,TAUREFCLD)
                ELSE
                  do K=1,L_LEVELS
                    TAUREFCLD(K) = 0.0D0
                  end do

                ENDIF

!  End of the "fill T, P arrays to pass to the new radiation code"
!  Set up, and solve for, the solar (visible) fluxes, if the sun
!  is up

               if(acosz.ge.1.0e-5) then


!  Check for ground ice.  Change albedo if there is any ice.

                 ALS   = ALSP(JCMN,ICMN)
                 IF(CO2ICE(J,I) .gt. 0.0) THEN
                   IF(JCMN.LT.JEQUATOR) THEN
                     ALS = ALICES
                   ELSE
                     ALS = ALICEN
                   END IF
                 ELSEIF (ALBFEED .and. QPIG(iMa_vap).gt. 
     *                   icethresh_kgm2 .and. 
     *                   .not. NPCFLAG(J,I) ) THEN
                   ALS = icealb
                 ENDIF
                 surfalb(J,I) = ALS

!  Get the optical depth (due to all sources) in the optical.

                 call optcv(DTAUV,TAUV,TAUCUMV,CO2V,PLEV,PFGASREF,
     *                      TGASREF,QXVDST,QSVDST,GVDST,WBARV,COSBV,
     *                      TAURAY,TAUREF,TMID,PMID,TAUGSURF,QH2O,
     *                      WREFH2O,
     *                      QEXTREFCLD,TAUREFCLD,QXVCLD,QSVCLD,GVCLD)

                 DO NW=1,L_NSPECTV
                   DO NG=1,L_NGAUSS
                     CUMTAUV(JCMN,ICMN,NW,NG) = TAUCUMV(L_LEVELS,NW,NG)
                   ENDDO
                 ENDDO

!  Calculate the fluxes in the visible

                 call sfluxv(DTAUV,TAUV,TAUCUMV,ALS,DWNV,WBARV,COSBV,
     *                       ACOSZ,SOL,GWEIGHT,NFLUXTOPV,FMNETV,
     *                       fluxupv,fluxdnv,diffvt,fzerov,taugsurf,
     *                       detau,jcmn,icmn)

                 SUNTOT(3) = FMNETV(1) - NFLUXTOPV

                 DO L=2,L_NLAYRAD
                   SUNTOT(2*L+1) = FMNETV(L) - FMNETV(L-1)
                 END DO

               else

!  If the sun is down, no solar flux, nor downward flux. . .

                 DO L=1,L_NLAYRAD
                   SUNTOT(2*L+1) = 0.0
                   FLUXDNV(L)    = 0.0
                 END DO
                 DIFFVT             = 0.0
                 fluxupv(1)         = 0.0
                 fluxdnv(1)         = 0.0
                 fluxupv(L_NLAYRAD) = 0.0
                 fluxdnv(L_NLAYRAD) = 0.0
               end if

               DNVFLUX(J,I) = FLUXDNV(L_NLAYRAD)
               DNDIFFV(J,I) = DIFFVT

!  History file output:  comp3cmn.h variables

               fuptopv(J,I)  = fluxupv(1)
               fdntopv(J,I)  = fluxdnv(1)
               fupsurfv(J,I) = fluxupv(L_NLAYRAD)
               fdnsurfv(J,I) = fluxdnv(L_NLAYRAD)

!  Set up, and solve for, the infrared fluxes

!  Check for ground ice.  Change the IR albedo if there is any ice.

               ALBI   = 1.0D0-EGOGND

               IF(CO2ICE(J,I) .gt. 0.0) THEN
                 IF(JCMN.LT.JEQUATOR) THEN
                   ALBI = 1.0D0-EGOCO2S
                 ELSE
                   ALBI = 1.0D0-EGOCO2N
                 END IF
               ENDIF

!  Get the optical depth (due to all sources) in the infrared.

               call optci(DTAUI,TAUCUMI,CO2I,PLEV,PFGASREF,TGASREF,
     *                    QEXTREFDST,QXIDST,QSIDST,GIDST,COSBI,WBARI,
     *                    TAUREF,TMID,PMID,TAUGSURFI,QH2O,WREFH2O,
     *                    QEXTREFCLD,TAUREFCLD,QXICLD,QSICLD,GICLD)

!  Calculate the fluxes in the IR.

               call sfluxi(PLEV,TLEV,DTAUI,TAUCUMI,UBARI,ALBI,WNOI,DWNI,
     *                     COSBI,WBARI,GWEIGHT,NFLUXTOPI,FMNETI,
     *                     FLUXUPI,FLUXDNI,FZEROI,TAUGSURFI)

               IRTOTAL(3) = FMNETI(1) - NFLUXTOPI
               DO L=2,L_NLAYRAD
                 IRTOTAL(2*L+1) = FMNETI(L) - FMNETI(L-1)
               END DO

               DNIRFLUX(J,I) = FLUXDNI(L_NLAYRAD)
               FLUXSURF(J,I) = (1.0-ALS)*FLUXDNV(L_NLAYRAD)

!  History file output:  comp3cmn_h variables

               fuptopir(J,I)  = fluxupi(1)
               fupsurfir(J,I) = fluxupi(L_NLAYRAD)
               fdnsurfir(J,I) = fluxdni(L_NLAYRAD)

!              THE FOLLOWING DEALING WITH HTH PERTAIN TO GRAPHICS.

               HTHST(J,I) = IRTOTAL(3)

               DO 180 K=1,L_LAYERS
                  HTH15(J,K)  = HTH15(J,K)+IR15(2*K+3)/FLOAT(IM)
                  HTHOUT(J,K) = HTHOUT(J,K)+IROUT(2*K+3)/FLOAT(IM)
  180          CONTINUE

!              CHANGE ATMOSPHERIC TEMPERATURES FOR SOLAR AND
!              INFRARED HEATING.
 
! Jan 2002  - Include heating rates in boundary layer scheme
! Store total radiative heating rates in QRAD and later pass them into
! the new boundary layer scheme.  Also, add in the non-LTE correction.

              DO K=2,NLEVSM1,2
                PNLTE   = SIGMA(K)*(P(J,I)-PTROP)+PTROP
                SUNLTE  = SUNTOT(K+1)*2.2E4*PNLTE/
     *                     (1.0+2.2E4*PNLTE)
                QRAD(K) = (SUNLTE+IRTOTAL(K+1))/
     *                    (CP*YM(K)*OM(K))
              END DO

! Non LTE fudge

!  Now, only update the stratosphere 

               K = 2
               PNLTE   = SIGMA(K)*(P(J,I)-PTROP)+PTROP
               SUNLTE  = SUNTOT(K+1)*2.2E4*PNLTE/
     *                     (1.0+2.2E4*PNLTE)
               TETA(K) = TETA(K)+NDT*(SUNLTE+IRTOTAL(K+1))/
     *                                  (CP*YM(K)*OM(K))

! End Jan 2002 update

!              UPDATE STRATOSPHERIC TEMPERATURE CHANGE.

               DELTAT(J,I,0) = OM(2)*(TETA(2)-TETAORG(2))/FLOAT(NC3)

!              SAVE POTENTIAL TEMPERATURE AND WIND VELOCITY BEFORE
!              DOING THE NON-RADIATIVE CALCULATIONS.

! Jan 2002  - Include heating rates in boundary layer scheme

               K = 2
               TETASAV(K) = TETA(K)

! End Jan 2002 update

!              END OF THE RADIATION CALCULATIONS.

! New boundary layer code

! Jan 2002  - Include heating rates in boundary layer scheme

!              NON-RADIATIVE CALCULATIONS.  CHANGES IN TEMPERATURE
!              AND WIND VELOCITY.
 
               PIE  = P(J,I)
               I1D = I
               J1D = J
 
!  Modify surface roughness when CO2 ice is on the ground.
!  Arya, pg 164.  Update added 7/21/99

               dumz0 = z0_pbl
               IF(CO2ICE(J,I) .gt. 0.0) dumz0 = 1.0E-4

!  Do the same for water ice if thick enough. 
!  Slab ice roughness  has been measured to be 0.01 mm.
               if(qpig(iMa_vap) .GT. 100
     *         .or. NPCFLAG(J,I) )  dumz0 = 1.0E-4

!  QRAD added Jan 2002 - pass radiation heating rates to the boundary
!  layer scheme
!  QPI added to argument list for dust tracer scheme AUG 2002.

               h2oflux = subflux(J,I)

               CALL NEWPBL(dumz0, epsl0_pbl, vk_pbl, alphl0_pbl,
     *                     ric_pbl, dtmu_pbl,
     *                     rmu1mu_pbl,dxi_pbl,du_pbl,I1D,J1D,
     *                     NSTEP,PTROP,GT(J1D,I1D),
     *                     PIE,NDT,SIGMA,PL,OM,UPI,VPI,TETA,
     *                     HT_PBL,SX_PBL,SY_PBL,RHOUCHpbl,KD,QRAD,
     *                     QPI,QPIG,LATHEAT(J1D,I1D),
     *                     NPCFLAG(J1D,I1D),sup_pbl,sdn_pbl,h2oflux)

               RHOUCH(J,I)  = RHOUCHpbl
               STRESSX(J,I) = SX_PBL
               STRESSY(J,I) = SY_PBL

               srfupflx(j,i,ima_vap)=sup_pbl
               srfdnflx(j,i,ima_vap)=sdn_pbl

! Calculate cloud tendencies

       if(MICROPHYSICS .eqv. .true.) then

         do l=2,2*l_layers+2,2
           deltateta(l)=teta(l)-tsave(jcmn,icmn,l)
           deltawv(l)=qpi(l,ima_vap)-qsave(jcmn,icmn,l,ima_vap)
           tetanew(l)=tsave(jcmn,icmn,l)
           qnew(l,ima_vap)=qsave(jcmn,icmn,l,ima_vap)
           do m=1,ntrace-1
             qnew(l,m)=qpi(l,m)
           enddo
         enddo

         do ns=1,nsplit

             do l = 4, 2*l_layers+2, 2
              tetanew(l)=tetanew(l)+deltateta(l)/float(nsplit)
              qnew(l,ima_vap)=qnew(l,ima_vap)+deltawv(l)/float(nsplit)
              tl(l) = om(l) * tetanew(l)
!              TL(L) = OM(L) * TETA(L)
             enddo
               CALL MICROPHYS(SDT,PL,TL,KD,SX_PBL,SY_PBL,HT_PBL,
     *              PCONSAV(J,I),GT(J,I),CO2ICE(J,I),dustref,QNEW,QPIG,
     *              QEXTREFDST)

             IF (LATENT_HEAT) THEN
               DO L = 4, 2*L_LAYERS+2, 2
                 TETA(L) = TL(L) / OM(L)
               ENDDO
             ENDIF

         enddo
   
         do k=4,nlevsm1,2
           do m=1,ntrace
              qpi(k,m)=qnew(k,m)
              qsave(jcmn,icmn,k,m)=qpi(k,m)
           enddo
           tsave(jcmn,icmn,k)=teta(k)
!           print*,jcmn,icmn,k,qsave(jcmn,icmn,k,ima_vap)
         enddo

         end if

!            COMPUTE TEMPERATURE ADJUSTMENTS DUE TO CONVECTION.

             CALL CONVECT(OM,YM,PLOGADJ,TETA,UPI,VPI,PL,TECON,PCON,
     *                    QPI)

!  Calculate heating and acceleration due to convective adjustment,
!  and store the convective boundary layer potential temperature
!  and zonal wind.

               DO 220 L=1,L_LAYERS
                  HCONADJ(J,L)   = HCONADJ(J,L)+CP*YM(2*L+2)*OM(2*L+2)*
     *                             (TETA(2*L+2)-TETASAV(2*L+2))/
     *                             (NDT*FLOAT(IM))
                  FCONADJ(J,L)   = FCONADJ(J,L)+YM(2*L+2)*(UPI(2*L+2)-
     *                             UPISAV(2*L+2))/(NDT*FLOAT(IM))
                  TETACON(2*L+2) = TETA(2*L+2)
                  UPICON(2*L+2)  = UPI(2*L+2)
  220          CONTINUE

!              STORE PRESSURE AT THE TOP OF THE SURFACE BOUNDARY LAYER.

               PCONSAV(J,I) = PCON
          
!              COMPUTE RAYLEIGH FRICTION SLOWING OF THE TOP 'LRAY'
!              LAYERS OF THE TROPOSPHERE
          
               DO 250 L=1,LRAY
                  RKEI        = 0.5*UPI(2*L+2)*UPI(2*L+2)+
     *                          0.5*VPI(2*L+2)*VPI(2*L+2)
                  UPI(2*L+2)  = UPI(2*L+2)-RAYK(L)*UPI(2*L+2)*NDT
                  VPI(2*L+2)  = VPI(2*L+2)-RAYK(L)*VPI(2*L+2)*NDT
                  RKEF        = 0.5*UPI(2*L+2)*UPI(2*L+2)+
     *                          0.5*VPI(2*L+2)*VPI(2*L+2)
                  TETA(2*L+2) = TETA(2*L+2)+(RKEI-RKEF)/(CP*OM(2*L+2))
  250          CONTINUE

!  Calculate dissipation and heating due to Rayleigh friction.
 
               DO 260 L=1,L_LAYERS
                  DISRAY(J,L) = DISRAY(J,L)-UPI(2*L+2)*YM(2*L+2)*
     *                          (UPI(2*L+2)-UPICON(2*L+2))/NDT-
     *                          VPI(2*L+2)*YM(2*L+2)*(VPI(2*L+2)-
     *                          VPICON(2*L+2))/NDT

                  HRAYFR(J,L) = HRAYFR(J,L)+CP*YM(2*L+2)*OM(2*L+2)*
     *                          (TETA(2*L+2)-TETACON(2*L+2))/
     *                          (NDT*FLOAT(IM))

                  FRAYFR(J,L) = FRAYFR(J,L)+YM(2*L+2)*(UPI(2*L+2)-
     *                          UPICON(2*L+2))/(NDT*FLOAT(IM))
  260          CONTINUE
           
!              COMPLETION STEPS.

!              Check local time, fill the 2pm opacity array if needed.
               univtime= NSTEPC3  - int(NSTEPC3 / NIT) * NIT
               if (univtime.eq.0) univtime = NIT
               if (univtime.eq.LON2PM(I)) then
                 taucld2pm(J,I) = taucloud(J,I,2)
                 taudst2pm(J,I) = taudust(J,I,2)
               endif

               if (mod(nstep,nc3*5).eq.0) then
               DO M = 1,NTRACE
                 test1(M)=test1(M)+qpig(M)*dxyp(J)
                 test2(M)=test2(M)+qcond(J,I,M)*dxyp(J)
                 DO K=1,L_LAYERS
                   L          = 2*K+2
                   atmass  = 100. * (PL(L+1)-PL(L-1)) / GRAV
                   test1(M)=test1(M)+qpi(L,M)*atmass*dxyp(J)
                   test2(M)=test2(M)+qtrace(j,i,K,M)*atmass*dxyp(J)
                  ENDDO
                ENDDO
              endif

!              STORE RESULTS FOR MAJOR VARIABLES IN PREPARATION FOR
!              SMOOTHING OVER THE GRID.
               
               DO 270 K=1,NLAY
                  L          = 2*K+2
                  DTETA(L)   = TETA(L)-TETASAV(L)
                  TTDELTA(K) = OM(L)*(TETA(L)-TETAORG(L))
                  DO M = 1,NTRACE
                    QTDELTA(J,I,K,M) = QPI(L,M)-QPISAV(L,M)
                  END DO
                  UTDELTA(K) = UPI(L)-UPISAV(L)
                  VTDELTA(K) = VPI(L)-VPISAV(L)
                  TT(J,I,K)  = TTDELTA(K)
                  UT(J,I,K)  = UTDELTA(K)
                  VT(J,I,K)  = VTDELTA(K)
  270          CONTINUE

               DO M = 1,NTRACE
                 QCDEL(J,I,M) = QPIG(M) - QCOND(J,I,M)
               END DO
       
!              TEMPERATURE OF THE AIR AT THE SURFACE
!              EXTRAPOLATE FROM ABOVE TWO MIDPOINTS.

               K       = NLEVSM2
               SLOPE   = (TETA(K-1)-TETA(K-3))/
     *                              (PLOGADJ(K-2)+PLOGADJ(K-3))
               TETA(K) = TETA(K-1)+SLOPE*PLOGADJ(K-1)

!              EXTRAPOLATE FOR THE SURFACE FROM THE TWO LEVELS
!              IMMEDIATELY ABOVE FOR A GOOD APPROXIMATION TO THE
!              PHYSICS.

               K       = NLEVELS
               SLOPE   = (TETA(K-1)-TETA(K-2))/PLOGADJ(K-2)
               TETA(K) = TETA(K-1)+SLOPE*PLOGADJ(K-1)

               TS(J,I) = TETA(K)*OM(K)
            
!              Don't allow TS(J,I) to fall below the saturation 
!              temperature.

               XPSAT = P(J,I)+PTROP
               XTSAT = 3182.48/(23.3494-LOG(XPSAT))

               IF(TS(J,I).lt.XTSAT) TS(J,I) = XTSAT

!              EXCHANGE OF HEAT BETWEEN SURFACE AND AIR.

               FA(J,I) = EG15GND*DNIRFLUX(J,I) - HT_PBL
         
!              TOTAL ABSORPTION OF SOLAR ENERGY BY THE ATMOSPHERE.

               SSUN(J,I) = 0.0

               DO 340 K=3,NLEVELS,2
                 SSUN(J,I) = SSUN(J,I)+SUNTOT(K)
  340          CONTINUE

!              THE FOLLOWING DEALING WITH HSOL, HRAD, AND HCOND
!              PERTAIN TO GRAPHICS.

               HSOLST(J,I) = SUNTOT(3)

!              HRAD(J,I)     = C6*COSPTG(J)*(FLUXES(NLEVELS)+
!    *                         SSUN(J,I)-IRNET(1))
!              HCOND(J,I)    = C7*COSPTG(J)*SDGR(J,I)

!               FLUXSURF(J,I) = FLUXES(NLEVELS)
               RIROUT(J,I)   = IRNET(1)

!              WRITE OUTPUT, IF DESIRED.

               IF(WANTIT) THEN
                 CALL CMP3OUT(J,IM,UTDELTA,VTDELTA,YM,NDT,FRY)
               END IF
         
  490       CONTINUE
  500    CONTINUE

!      if (mod(nstep,nc3*5).eq.0 .and. microphysics .eqv. .true.) then
!       write(*,'("At hour #",f10.1)')TAU
!       write(*,'("LS =",f10.1)')VPOUT
!       print*,'Check tot before and after microphys'
!       print*,test2(3)+test2(6),test1(3)+test1(6),
!    *         test2(3)+test2(6)-test1(3)-test1(6)
!       print*,'Dust'
!       print*,test2(1)+test2(5),test1(1)+test1(5),
!    *          test2(1)+test2(5)-test1(1)-test1(5)
!       print*,'--------------------------------------------'
!     endif

!        THE END OF THE MAIN COMPUTATIONAL LOOPS (I,J).

!        PUT DELTAS INTO FINAL FORM.  DIVIDE DELTAS BY NC3 TO GIVE
!        THE CHANGE PER TIME STEP.  THE U AND V DELTAS ARE CHANGED
!        FROM THE 'PI' GRID TO THE U,V-GRID.

         RINC3 = 1.0/FLOAT(NC3)
        
         do i=1,im
           do j=1,jm-1
             do m=1,ntrace
               qcdel(j,i,m) = qcdel(j,i,m)*rinc3
             enddo
           enddo
         enddo
 
         DO 570 L=1,NLAY

! *** BEGIN C GRID CORE INTERPOLATION CHANGE
! Do U, T and V separately for the moment as they lie on different grids
! QTDELTA - the mixing ratio tendency - is passed into newstep
! We do not change mixing ratios here.  But we do change U, V, & T.

           DO 560 I=1,IM
             DO 550 J=1,JM-1
               DO M = 1,NTRACE
                 QTDELTA(J,I,L,M) = QTDELTA(J,I,L,M)*RINC3/DT
               END DO
               DELTAT(J,I,L) = TT(J,I,L)*RINC3
               if(i.eq.im)then
                 DELTAU(J,I,L) = 0.5*RINC3*(UT(J,IM,L)+UT(J,1,L))
               else
                 DELTAU(J,I,L) = 0.5*RINC3*(UT(J,I,L)+UT(J,I+1,L))
               endif
  550        CONTINUE
  560      CONTINUE

! Now do V away from the poles

           DO 562 J=2,JM-1
             DO 552 I=1,IM
               DELTAV(J,I,L) = 0.5*RINC3*(VT(J,I,L)+VT(J-1,I,L))
  552        CONTINUE
  562      CONTINUE

! dV/dt at the poles is zero

           DO I=1,IM
             DELTAV(1,I,L)  = 0.0
             DELTAV(JM,I,L) = 0.0
           end do

! *** END C GRID CORE INTERPOLATION CHANGE

  570    CONTINUE

!  End of full comp3 
!  Used in TEMPGR to know when to set subflux(j,i) = 0

        fullcomp3 = .TRUE.
      
      END IF

!################################-#-#-#-################################

!     END THE MAJOR IF LOOP (ENTERED EVERY NC3 TIME STEPS).

!     THIS SECTION OF CODE IS DONE EACH TIME THROUGH COMP3.

!     ADD DELTAS TO T, U, V, AND TSTRAT EACH TIME STEP.
!     U,V AND DELTAU, DELTAV ARE AT THE U,V GRID POINTS
!     T, TSTRAT AND DELTAT ARE AT THE 'PI' POINTS.

      DO 640 I=1,IM

         DO 630 J=1,JM-1
            TSTRAT(J,I) = TSTRAT(J,I)+DELTAT(J,I,0)
  630    CONTINUE

  640 CONTINUE

!     CALCULATE GEOPOTENTIALS AT LAYER MIDPOINTS FOR EVERY 'PI'
!     GRID POINT.

      CALL GEOPOTN(PTROP,SIGMA,P,T,DSIG,NLAY,TOPOG,GEOT)

      RETURN
      END
      subroutine watermass
      use grid_h
      use defines_h
      use standard_h
      use constants_h

      implicit none

      integer, parameter :: iMacld = 3
      integer, parameter :: iMavap = 6

      real*8  :: atm, cld, ice, total, icenpc
      real*8  :: atmc, cldc, atmos
      real*8  :: pl(L_LEVELS)
      integer :: I, J, K, L

      atm    = 0.0D0
      cld    = 0.0D0
      ice    = 0.0D0
      icenpc = 0.0D0

      do J=1,L_J-1
        do I = 1,L_I
          PL(3) = PTROP
          DO K = 4, L_LEVELS
            PL(K)  = PTROP+SIGMA(K)*P(J,I)
          end do

          atmc = 0.0D0
          cldc = 0.0D0
          do L=1,L_LAYERS
            k = 2*L+2
            atmos = 100.0D0*(PL(K+1)-PL(K-1)) / GRAV
            atmc  = atmc + qtrace(J,I,L,iMavap)*atmos
            cldc  = cldc + qtrace(J,I,L,iMacld)*atmos
          end do
          atm = atm + atmc*dxyp(j)
          cld = cld + cldc*dxyp(j)
          if(npcflag(j,i) .and. qcond(j,i,iMavap).lt.0.0) then
            icenpc = icenpc + qcond(j,i,iMavap)*dxyp(j)
          elseif(qcond(j,i,iMavap).gt.0.0) then
            ice = ice + qcond(j,i,iMavap)*dxyp(j)
          end if
        end do
      end do
      total = atm + cld + ice + icenpc

      write(6,'("Water:  atm",8x,"cld",11x,"ice",10x,"icenpc",8x,
     *           "Total")')
      write(6,'(1pd12.5,4(2x,1pd12.5))'), atm, cld, ice, icenpc, total
      write(6,'(69("="))')

      return
      end
      subroutine firstcomp3(jm,im,lspv,ngauss,detau,dnirflux,
     *                      dndiffv)

!  Subroutine created 2/6/13
!  First time step in comp3, set variables to well defined values
!  but not the correct values.  After this initial pass through
!  TEMPGR, the code goes through a full comp3, and these values are
!  replaced by correct ones.  The first time step uses these values,
!  all other time steps use the recently computed (correct) values.

      implicit none
      integer  :: JM, IM, LSPV, ngauss, j, i, nw, ng
      real*8   :: detau(JM,IM,LSPV,NGAUSS)
      real*8   :: dnirflux(JM,IM), dndiffv(JM,IM)

      do j = 1,JM
        do i=1,IM
          dnirflux(J,I) = 1.0D0
          dndiffv(J,I)  = 0.0D0
          do nw=1,LSPV
            do ng=1,ngauss
              detau(j,i,nw,ng) = 0.1D0
            end do
          end do
        end do
      end do

      return
      end
