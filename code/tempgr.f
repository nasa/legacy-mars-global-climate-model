      SUBROUTINE TEMPGR(GT,CO2ICE,IM,JM,PTROP,P,T,ASYM,TAUCRT,ALSP,FA,
     *                  EG15GND,EGOGND,SQRDY,ZIN,DT,ROT,DMADT,SDGR,
     *                  TINF,STEMP,STHICK,SCOND,DNIRFLUX,DNVFLUX,RHOUCH,
     *                  QCOND,NPCFLAG,LATHEAT,CPSOIL,RHOSOIL,QTRACE)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

!  PURPOSE
!      TEMPGR CALCULATES THE GROUND TEMPERATURE FOR EACH 'PI' GRID POINT
!      ON THE SURFACE OF MARS. WHERE CO2 ICE IS PRESENT ON THE GROUND,
!      THIS ALSO INVOLVES CALCULATION OF THE MASS OF CO2 ICE ON THE
!      GROUND AND SIGMA DOT AT THE GROUND.
!  INPUT PARAMETERS
!      GT(J,I) ARRAY  - GROUND TEMPERATURE OR (-1) TIMES CO2 ICE
!                             MASS AT THE 'PI' POINTS.
!      P(J,I) ARRAY   - THE CURRENT 'PI' VALUE (SURFACE PRESSURE
!                             MINUS TROPOPAUSE PRESSURE) AT THE 'PI'
!                             POINTS.
!      ALSP(J,I) ARRAY- THE LOCAL SURFACE ALBEDO FOR THE 'PI' POINTS
!                             ASSUMING NO CO2 ICE IS PRESENT.
!      T(J,I,K) ARRAY - THE TEMPERATURE AT THE MIDPOINT OF EACH
!                             LAYER FOR THE 'PI' POINTS. (K = 1 TO NLAY)
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
!                       FOR EACH 'PI' POINT.
!      JBPS           - APPROXIMATE LATITUDE INDEX OF EDGE OF SOUTHERN
!                       POLAR CAP.
!      JBPN           - APPROXIMATE LATITUDE INDEX OF EDGE OF NORTHERN
!                       POLAR CAP.
!
!  OUTPUT PARAMETERS
!      GT(J,I) ARRAY  - AS IN INPUT, BUT ADJUSTED ACCORDING TO
!                       HEATING AT THE GROUND.
!      TINF(J,I) ARRAY- AS IN INPUT, BUT ADJUSTED IN PROPORTION TO
!                       TO THE NEW GROUND TEMPERATURE.
!      SDGR(J,I) ARRAY- THE VERTICAL (SIGMA) COMPONENT OF WIND VELOCITY
!                       AT THE GROUND CALCULATED FROM THE CHANGE IN THE
!                       AMOUNT OF CO2 ON THE GROUND.
!  CALLED BY
!      COMP3
!
!  SUBROUTINES CALLED
!      GETFLUX
!

      use grid_h
      use defines_h
      use constants_h, only: GRAV, STBO, XLHTC, CP
      use standard_h, only: ALICEN, ALICES, EGOCO2N, EGOCO2S, EG15CO2N,
     *                      EG15CO2S, JEQUATOR, SURFALB,ALBFEED,ICEALB,
     *                      ICETHRESH_KGM2
      use cldcommon_h

      implicit none

!######################################################################

      REAL*8  :: GT(L_JSIZE,L_ISIZE), CO2ICE(L_JSIZE,L_ISIZE)
      REAL*8  :: P(L_JSIZE,L_ISIZE)
      REAL*8  :: T(L_JSIZE,L_ISIZE,L_LAYERS)
      REAL*8  :: ALSP(L_JSIZE,L_ISIZE)
      REAL*8  :: TINF(L_JSIZE,L_ISIZE), FA(L_JSIZE,L_ISIZE)
      REAL*8  :: ZIN(L_JSIZE,L_ISIZE,NL), DMADT(L_JSIZE,L_ISIZE)
      REAL*8  :: SDGR(L_JSIZE,L_ISIZE)
      REAL*8  :: DTSDT(NL)

      REAL*8  :: RHOUCH(L_JSIZE,L_ISIZE), RHOUCHT
      REAL*8  :: DOWNIR
      real*8  :: DNIRFLUX(L_J,L_I), DNVFLUX(L_J,L_I)

      REAL*8  :: QCOND(L_JSIZE,L_ISIZE,NTRACE)
      LOGICAL :: NPCFLAG(L_JSIZE,L_ISIZE)

      REAL*4  :: LATHEAT(L_J,L_I)

!  10-07-11 latent heat updates

      real*8  :: qtrace(l_jsize,l_isize,l_layers,ntrace)
      real*8  :: h2oflux
      real*8  :: qgnd, wflux

!  FCC updates

      REAL*8  :: FLUX(2*NL+1)
      REAL*8  :: RHOSOIL(L_J,L_I,NL), CPSOIL(L_J,L_I,NL)
      REAL*8  :: STEMP(L_JSIZE,L_ISIZE,2*NL+1)
      REAL*8  :: STHICK(2*NL+1)
      REAL*8  :: SCOND(L_JSIZE,L_ISIZE,2*NL+1)

!     Implicit none

      integer :: i, j, k, l, JBCS, JBCN, IM, JM
      real*8  :: EMGOUT, TINP, SQRDY,  FCDN, TGP, DT, DMGDT, XA
      real*8  :: EMG15, EG15GND, EGOGND, X1, ASYM, TAUCRT, ROT, GTEMP
      real*8  :: PSAT, TSAT, PTROP, RGTOT, TGINF, ALS, ALCDS

!  Latent heating updates
!  Passed into TEMPGR vi cldcommon_h module

!     logical :: fullcomp3
!     real*8  :: subflux(L_JSIZE,L_ISIZE)
!     real*8  :: gndice(L_JSIZE,L_ISIZE)

!#=====================================================================

!  First time in tempgr after a full trip through comp3?

      if(FULLCOMP3) then
        FULLCOMP3 = .FALSE.
        do J=1,JM-1
          do I=1,IM
            subflux(j,i) = 0.0D0
            gndice(J,I)  = QCOND(J,I,iMa_vap)
          end do
        end do
      end if
 
!     GROUND TEMPERATURE CALCULATION DONE EVERY TIME STEP

      DO 4150 I=1,IM
         DO 4151 J=1,JM-1

            DMGDT     = 0.0D0
            SDGR(J,I) = 0.0D0

!           FIRST CALCULATE 'SURFACE' ALBEDO.

            ALS = ALSP(J,I)

!           Check for CO2 ice on the ground.

            IF(CO2ICE(J,I) .gt. 0.0) THEN

              IF(J.LT.JEQUATOR) THEN
                ALS = ALICES
              ELSE
                ALS = ALICEN
              END IF
            ELSEIF (ALBFEED .and. QCOND(J,I,iMa_vap).gt.
     *              icethresh_kgm2
     *              .and. (NPCFLAG(J,I) .eqv. .FALSE.)) THEN
                 ALS = icealb
            END IF

            surfalb(J,I) = ALS

!           TINF is the deep subsurface temperature, FA the heat
!           exchange from convection and thermal radiation from
!           the atmosphere.

            RGTOT = FA(J,I)

!           Get TSAT, the CO2 frost point at this surface pressure.

            PSAT = P(J,I)+PTROP
            TSAT = 3182.48D0/(23.3494D0-LOG(PSAT))

!  Begin major modifications

            IF(CO2ICE(J,I).le.0.0) THEN  ! no CO2 on the ground

              CO2ICE(J,I) = 0.0D0

!             Emissivities for bare ground

              EMG15  = EG15GND
              EMGOUT = EGOGND

!             CHANGE TG FOR HEATING & COOLING EFFECTS.

!             New soil scheme: surface boundary condition
!             Changes to the new ground temperature computation

              DOWNIR  = EMG15*DNIRFLUX(J,I)
              RHOUCHT = RHOUCH(J,I)*T(J,I,L_LAYERS)

!             Newton-Raphson method of solving for Tg, with new 
!             Tg-dependent terms (RHOUCH, RHOUCHT replace RGTOT) added
!             to the mix

              CALL NEWTG(ALS,DNVFLUX(J,I),DOWNIR,RHOUCH(J,I),RHOUCHT,
     *                   SCOND(J,I,2),STEMP(J,I,2),STHICK(2),
     *                   GT(J,I),J,I,PSAT,QTRACE(J,I,L_LAYERS,iMa_vap),
     *                   gndice(J,I),dt,subflux(j,i),npcflag(j,i))

              IF(GT(J,I).lt.TSAT) THEN

                GT(J,I) = TSAT

!               Emissivities of bare ground since CO2 ice just starting

                EMG15  = EG15GND
                EMGOUT = EGOGND

                FCDN = -2.0D0*SCOND(J,I,2)*(STEMP(J,I,2)-TSAT)/
     *                        STHICK(2)

                TGP  = DT*((1.0D0-ALS)*DNVFLUX(J,I)+RGTOT-
     *                 EMG15*(STBO*TSAT**4) - FCDN)/XLHTC

!               Check to see if there is any CO2 ice accumulation.

                IF(TGP.LT.0.0) THEN

                  DMGDT       = -TGP/DT
                  CO2ICE(J,I) = -TGP

                ELSE

                  DMGDT     = 0.0D0

!                 This term represents the last amounts of ice
!                 evaporating resulting in an increase in Tg. It still 
!                 depends on TINP until I can figure out what to do 
!                 with it.

!                 added 5/13/04,modified 5/13/08
                  TINP        = SQRDY/zin(j,i,1)
                  GT(J,I)     = TSAT+(TGP*XLHTC*TINP)
                  CO2ICE(J,I) = 0.0D0

                ENDIF

              ENDIF

            ELSE    ! CO2 ice on the ground

              GT(J,I) = TSAT

!  Only modify if we have water condensation
!  Latent heating update Oct. 2011

              wflux = 0.0

              qgnd  = (18.0/44.0)*
     *                6.11*exp(22.5*(1.0-(273.16/GT(j,i))))/psat
              wflux = -rhouch(j,i)*
     *                 (QTRACE(J,I,L_LAYERS,iMa_vap)-qgnd)/cp
              if(wflux.lt.0.0) then
                gndice(j,i) = gndice(j,i) - wflux*dt

!  Note:  wflux > 0 is sublimation
!         wflux < 0 is condensation

                SUBFLUX(J,I) = SUBFLUX(J,I) + wflux*dt
              end if

!             Emissivities of CO2 ice on the ground.

              IF(J.LT.JEQUATOR) THEN
                EMG15  = EG15CO2S
                EMGOUT = EGOCO2S
               ELSE
                EMG15  = EG15CO2N
                EMGOUT = EGOCO2N
              END IF

!             New soil scheme: surface boundary condition with ice on 
!             ground

              FCDN = -2.0D0*SCOND(J,I,2)*(STEMP(J,I,2)-TSAT)/
     *                      STHICK(2)

              TGP  = -CO2ICE(J,I) + DT*((1.0D0-ALS)*DNVFLUX(J,I)+RGTOT-
     *                EMG15*(STBO*TSAT**4) - FCDN)/XLHTC

!             Check to see if there is any CO2 ice still left.

              IF(TGP.LT.0.0) THEN

                DMGDT       = -(CO2ICE(J,I)+TGP)/DT
                CO2ICE(J,I) = -TGP

              ELSE

                DMGDT = -CO2ICE(J,I)/DT

!               This term represents the last amounts of ice evaporating
!               resulting in an increase in Tg. It still depends on TINP
!               until I can figure out what to do with it.

                TINP        = SQRDY/zin(j,i,1)
                GT(J,I)     = TSAT+(TGP*XLHTC*TINP)
                CO2ICE(J,I) = 0.0D0

              ENDIF

            ENDIF

!  End major modification

            X1 = GT(J,I)
!
! Begin new soil scheme
!

!  Calculate fluxes at layer boundaries
!  Fluxes defined as positive downward

           do K=3,2*NL-1,2
             flux(K) = -scond(J,I,k)*(stemp(J,I,K+1) - stemp(J,I,K-1))/
     *                  sthick(K)
           end do

!  Calculate flux at top and bottom

           flux(1)      = -2.0D0*scond(J,I,2)*(stemp(J,I,2) - X1)/
     *                         sthick(2)
           flux(2*NL+1) = 0.0

!  Update soil temperatures

           do L=1,NL
             K = 2*L
             stemp(J,I,K) = stemp(J,I,K) - dt*(flux(k+1) - flux(K-1))/
     *                      (rhosoil(J,I,L)*cpsoil(J,I,L)*sthick(K))
           end do

!
! End new soil scheme

      SDGR(J,I) = GRAV * (DMGDT + DMADT(J,I))

 4151 CONTINUE
 4150 CONTINUE

      RETURN
      END
