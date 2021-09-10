      SUBROUTINE NEWSTEP

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

C  PURPOSE:
C  Step time 
C  CALLED BY
C      MAIN
C  SUBROUTINES CALLED:
C     DYCORE
C
      use grid_h
      use defines_h
      use constants_h, only: PI, GRAV, RGAS, CP, KAPA
      use standard_h
      implicit none

C Intensity of Shapiro filters (order and timescale in hours)
      integer FILTML_ORDER
      parameter(FILTML_ORDER=8)
      real*4 FILTML_TAU
c      parameter(FILTML_TAU=2.4*3600.0)
      parameter(FILTML_TAU=12.0*3600.0)  ! works with 16x24, 40x24, 60x36
c      parameter(FILTML_TAU=2.4*3600.0)  ! initial setup
      integer FILTPQ_ORDER
      parameter(FILTPQ_ORDER=8)
      real*4 FILTPQ_TAU
      parameter(FILTPQ_TAU=12.0*3600.0) ! initial setup, works with 16x24, 40x24
c      parameter(FILTPQ_TAU=12.0*3600.0) ! works with 60x36
C Exner function multiplier
      real*4 pkh(l_isize,l_jsize,l_layers+1)
C Tendencies
      real*4 pii(l_isize,l_jsize)
      real*4 uoi(l_isize,l_jsize,l_layers)
      real*4 voi(l_isize,l_jsize,l_layers)
      real*4 poi(l_isize,l_jsize,l_layers)
      real*4 qoi(l_isize,l_jsize,l_layers,ntrace)
C Converting Real*8 common block variables to real*4 for input to DYCORE
      real*4 RFPTROP,RFDT,RFOM,RFCP,RFRGAS,RFRAD,RFPSF
C Used in Shapiro filters
      real*4 VT1(l_isize,l_jsize)
      real*4 VT2(l_isize,l_jsize)
      real*4 VT3(l_isize,l_jsize)
      real*4 VT4(l_isize,l_jsize)
      real*4 VT5(l_isize,l_jsize)
      real*4 VT6(l_isize,l_jsize)
      real*4 VT7(l_isize,l_jsize)
      real*4 qci(l_isize,l_jsize,ntrace)

      REAL*8 NEGMASS(NTRACE) 
      REAL*8 TOTMASS(NTRACE)
      REAL*8 NEGFAC(NTRACE)
      REAL*8 MASSINT

!     implicit none

      integer :: i, j, k, l, lm
      real*8  :: TIMFILT, PRESB, PRESM, PTERM
      REAL*4  :: OMG(L_I,L_J,L_LAYERS)
      REAL*4  :: VOR(L_I,L_J,L_LAYERS)

!======================================================================!

      IM = L_ISIZE
      JM = L_JSIZE
      LM = L_LAYERS

C Set up array positions for leapfrog timestepping

      JY2 = JY1
      JY1 = MOD(JY2,2) + 1

C Diagnostics from DYCORE (mainly Omega and Vorticity)

      DIAG = .FALSE.

C Robert Filter

      TIMFILT = 0.05

C Calculate Exner function multiplier AT LAYER BOUNDARIES

      DO L = 1,LM+1
        DO I = 1,IM
          DO J = 1,JM-1
            PRESB      = (DYSIG(L)*P(J,I))+PTROP
            PKH(I,J,L) = (PRESB/PSF)**KAPA
          END DO
        END DO
      END DO

C Set up real*4 variables

      RFPTROP = PTROP
      RFDT    = 2.0*DT
C     RFDT    = DT
      RFOM    = 2.0*PI/DAY
      RFCP    = CP
      RFRGAS  = RGAS
      RFRAD   = RAD
      RFPSF   = PSF

C Update X(t) from GCM physics fields ready for DYCORE
C And preset tendencies to zero

      DO I = 1,IM
        DO J = 1,JM-1
          PIB(I,J,JY2) = P(J,I)
c          PII(I,J) = 0.0
          PII(I,J) = - SDGR(J,I) * 0.01

          DO K =1, NTRACE
            QOC(I,J,K,JY2) = QCOND(J,I,K)
            QCI(I,J,K)     = QCDEL(J,I,K) / DT
          END DO

          DO L = 1,LM
            UOB(I,J,L,JY2) = U(J,I,L)
            VOB(I,J,L,JY2) = V(J,I,L)
            POB(I,J,L,JY2) = T(J,I,L)

            UOI(I,J,L) = DELTAU(J,I,L)/DT
            VOI(I,J,L) = DELTAV(J,I,L)/DT

            PRESM = ( SIGMA(2*L+2) * P(J,I) ) + PTROP
            PTERM = ( (PSF/PRESM)**KAPA )

!            POI(I,J,L) = DELTAT(J,I,L) * P(J,I) * PTERM / DT

            DHEAT(J,I,L)=DELTAT(J,I,L)/DT*88775.


            POI(I,J,L) = POB(I,J,L,JY2) * PII(I,J) + P(J,I) *
     *               (DELTAT(J,I,L) * PTERM / DT - (POB(I,J,L,JY2)
     *               * KAPA / PRESM * PII(I,J) * SIGMA(2*L+2)))


c            UOI(I,J,L) = 0.0
c            VOI(I,J,L) = 0.0
c            POI(I,J,L) = 0.0

            DO K = 1,NTRACE
              QOB(I,J,L,K,JY2) = QTRACE(J,I,L,K)
              QOI(I,J,L,K)     = P(J,I) * QTDELTA(J,I,L,K)
            END DO
          END DO
        END DO
      END DO

      CALL DYCORE ( IM,JM,JM,LM,DYSIG,
     *                    RFPTROP,NTRACE,RFDT,
     *                    RFOM, RFCP, RFRGAS, RFRAD, RFPSF,
     *                    PHS  ,PKH,
     *           PIB(1,1,JY2),UOB(1,1,1,JY2),VOB(1,1,1,JY2),
     *           POB(1,1,1,JY2),QOB(1,1,1,1,JY2),
     *           PIB(1,1,JY1),UOB(1,1,1,JY1),VOB(1,1,1,JY1),
     *           POB(1,1,1,JY1),QOB(1,1,1,1,JY1),
     .           PII,  UOI,  VOI,  POI,  QOI,
     .                    OMG,  VOR,  DIAG, JY1)
 
C Add tendencies onto fields using leapfrog method. 
C Set X(t-1) to average X(t) using Robert filter. Honest. No really.
C It is optimal with ALPHA = 0.2475 (The Brown Campana Scheme) in DYCORE
C (See Aries documentation Volume 5 page 37)
C Functional form taken from Haltiner and Williams p146 (1979)

C Firstly do PI

      DO I=1,IM*(JM-1)
        PII(I,1) = PIB(I,1,JY1) + PII(I,1)*RFDT
        VT4(I,1) = 1.0 / PII(I,1)
        VT1(I,1   ) = PIB(I,1,JY2)
        VT3(I,1   ) = PIB(I,1,JY1)
        PIB(I,1,JY2) = (PIB(I,1,JY1) + PII(I,1)) * (0.5*TIMFILT)
     *            +  PIB(I,1,JY2)             * (1.0 - TIMFILT)
        VT2(I,1   ) = 1.0 / PIB(I,1,JY2)
        PIB(I,1,JY1) = PII(I,1)
      ENDDO

      DO I =1, IM*(JM-1)
        DO K = 1,NTRACE
          QCI(I,1,K) = QOC(I,1,K,JY1) + QCI(I,1,K)*RFDT
        END DO
      END DO

C *** Start BIG L LOOP

      DO L = 1,LM
 
C Add tendencies onto fields

      DO I =1,IM*(JM-1)
        POB(I,1,L,JY1)= POB(I,1,L,JY1) * VT3(I,1   )
        POI(I,1,L   )= POB(I,1,L,JY1) + POI(I,1,L )*RFDT
        DO K = 1,NTRACE
          QOB(I,1,L,K,JY1)= QOB(I,1,L,K,JY1) * VT3(I,1   )
          QOI(I,1,L   ,K) = QOB(I,1,L,K,JY1) + QOI(I,1,L,K)*RFDT
        END DO
        UOI(I,1,L   )= UOB(I,1,L,JY1) + UOI(I,1,L )*RFDT
        VOI(I,2,L   )= VOB(I,2,L,JY1) + VOI(I,2,L )*RFDT
      ENDDO

C V at poles is zero

      DO I=1,IM
        VOI(I,JM,L) = 0.0
        VOI(I, 1,L) = 0.0
      ENDDO

C Call Shapiro filters: Off for tracers (including theta) for now

      CALL FILTPQ(POI(1,1,L),VT5, IM,JM,IYE,IYW, VT6, VT7,
     &             FILTPQ_ORDER, FILTPQ_TAU)

      do I = 1,IM*(JM-1)
        POI(I,1,L) = POI(I,1,L) + VT5(I,1)*RFDT
      end do

      CALL FILTML(UOI(1,1,L),VT5,.TRUE.,IM,JM,IYE,IYW,VT6,VT7,
     &             FILTML_ORDER, FILTML_TAU)
      DO I=1,IM*(JM-1)
        UOI(I,1,L) = UOI(I,1,L) + VT5(I,1)*RFDT
      ENDDO

      CALL FILTML(VOI(1,1,L),VT5,.FALSE.,IM,JM,IYE,IYW,VT6,VT7,
     &             FILTML_ORDER, FILTML_TAU)
      DO I =1,IM*(JM-2)
        VOI(I,2,L) = VOI(I,2,L) + VT5(I,2)*RFDT
      ENDDO

C Use Robert filter on all tendencies apart from PI

      DO I =1,IM*(JM-1)

        POI(I,1,L )= POI(I,1,L) * VT4(I,1)
        DO K = 1,NTRACE
          QOI(I,1,L,K)= QOI(I,1,L,K) * VT4(I,1)
        END DO
        POB(I,1,L,JY2) = (POI(I,1,L ) * PII(I,1) + POB(I,1,L,JY1))
     *                                            * (0.5*TIMFILT)
     *              + POB(I,1,L,JY2) * VT1(I,1) * (1.0 - TIMFILT)
        POB(I,1,L,JY2) = POB(I,1,L,JY2) * VT2(I,1)
        POB(I,1,L,JY1) = POI(I,1,L  )

        DO K = 1,NTRACE
          QOB(I,1,L,K,JY2) = (QOI(I,1,L,K)*PII(I,1) + QOB(I,1,L,K,JY1))
     *                                            * (0.5*TIMFILT)
     *                    + QOB(I,1,L,K,JY2)*VT1(I,1)*(1.0 - TIMFILT)
          QOB(I,1,L,K,JY2) = QOB(I,1,L,K,JY2)*VT2(I,1)
          QOB(I,1,L,K,JY1) = QOI(I,1,L,K)
        END DO
        UOB(I,1,L,JY2) = (UOB(I,1,L,JY1) + UOI(I,1,L)) * (0.5*TIMFILT)
     *              +  UOB(I,1,L,JY2)               * (1.0 - TIMFILT)
        UOB(I,1,L,JY1) = UOI(I,1,L   )
      ENDDO

      DO  I =1,IM*(JM-2)
        VOB(I,2,L,JY2) = (VOB(I,2,L,JY1) + VOI(I,2,L)) * (0.5*TIMFILT)
     *              +  VOB(I,2,L,JY2)               * (1.0 - TIMFILT)
        VOB(I,2,L,JY1) = VOI(I,2,L   )
      ENDDO

C *** End BIG L LOOP

      ENDDO

      DO I =1,IM*(JM-1)
        DO K = 1,NTRACE
          QOC(I,1,K,JY2) = (QCI(I,1,K) + QOC(I,1,K,JY1))
     *                                            * (0.5*TIMFILT)
     *                    + QOC(I,1,K,JY2)*(1.0 - TIMFILT)
          QOC(I,1,K,JY1) = QCI(I,1,K)
        END DO
      ENDDO

C Initialize total neg mass and tot mass arrays

      DO K = 1,NTRACE
        NEGMASS(K) = 0.
        TOTMASS(K) = 1.d-50
      ENDDO

C Now write values to physics fields

      DO I = 1,IM
        DO J = 1,JM-1
          P(J,I) = PIB(I,J,JY1)

          DO K = 1, NTRACE
            QCOND(J,I,K) = QOC(I,J,K,JY1)
          END DO

          DO L = 1,LM
            U(J,I,L) = UOB(I,J,L,JY1)
            V(J,I,L) = VOB(I,J,L,JY1)
            T(J,I,L) = POB(I,J,L,JY1)
            MASSINT  = 100. * DSIG(L)* (P(J,I) + PTROP) / GRAV * DXYP(J)
            DO K = 1,NTRACE
              QTRACE(J,I,L,K) = QOB(I,J,L,K,JY1) 
C Eliminate negative mixing ratios

              IF (QTRACE(J,I,L,K).LT.0.) THEN
                NEGMASS(K) = NEGMASS(K) + QTRACE(J,I,L,K) * MASSINT
                QTRACE(J,I,L,K) = 1.d-50
              ENDIF
              TOTMASS(K) = TOTMASS(K) + QTRACE(J,I,L,K) * MASSINT

            END DO
          END DO
        END DO    ! END loop on longitude
      END DO      ! END loop on latitude

c Scale tracer fileds by total negative mass ratio

      DO K = 1,NTRACE
        NEGFAC(K) = max((TOTMASS(K)+NEGMASS(K))/TOTMASS(K),1.d-50)
      ENDDO

      DO K = 1,NTRACE
        DO L = 1, LM
          DO I = 1,IM
            DO J = 1,JM-1
              QTRACE(J,I,L,K) = QTRACE(J,I,L,K) * NEGFAC(K)
            ENDDO
          ENDDO
        ENDDO
      ENDDO


      RETURN
      END
