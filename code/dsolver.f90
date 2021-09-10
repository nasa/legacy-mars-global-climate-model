      SUBROUTINE DSOLVER(NL,GAMA,CP,CM,CPM1,CMM1,E1,E2,E3,E4,BTOP,     &
                         BSURF,RSF,XK1,XK2)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

! DOUBLE PRECISION VERSION OF SOLVER

      implicit none

      integer, parameter :: NMAX=201

      integer :: L, NL, I, LM2, LM1, N
      real*8  :: GAMA(NL),CP(NL),CM(NL),CPM1(NL),CMM1(NL),XK1(NL)
      real*8  :: XK2(NL),E1(NL),E2(NL),E3(NL),E4(NL)
      real*8  :: AF(NMAX),BF(NMAX),CF(NMAX),DF(NMAX),XK(NMAX)
      real*8  :: BTOP, BSURF, RSF

!*********************************************************
!* THIS SUBROUTINE SOLVES FOR THE COEFFICIENTS OF THE    *
!* TWO STREAM SOLUTION FOR GENERAL BOUNDARY CONDITIONS   *
!* NO ASSUMPTION OF THE DEPENDENCE ON OPTICAL DEPTH OF   *
!* C-PLUS OR C-MINUS HAS BEEN MADE.                      *
!* NL     = NUMBER OF LAYERS IN THE MODEL                *
!* CP     = C-PLUS EVALUATED AT TAO=0 (TOP)              *
!* CM     = C-MINUS EVALUATED AT TAO=0 (TOP)             *
!* CPM1   = C-PLUS  EVALUATED AT TAOSTAR (BOTTOM)        *
!* CMM1   = C-MINUS EVALUATED AT TAOSTAR (BOTTOM)        *
!* EP     = EXP(LAMDA*DTAU)                              *
!* EM     = 1/EP                                         *
!* E1     = EP + GAMA *EM                                *
!* E2     = EP - GAMA *EM                                *
!* E3     = GAMA*EP + EM                                 *
!* E4     = GAMA*EP - EM                                 *
!* BTOP   = THE DIFFUSE RADIATION INTO THE MODEL AT TOP  *
!* BSURF  = THE DIFFUSE RADIATION INTO THE MODEL AT      *
!*          THE BOTTOM: INCLUDES EMMISION AND REFLECTION *
!*          OF THE UNATTENUATED PORTION OF THE DIRECT    *
!*          BEAM. BSTAR+RSF*FO*EXP(-TAOSTAR/U0)          *
!* RSF    = REFLECTIVITY OF THE SURFACE                  *
!* XK1    = COEFFICIENT OF THE POSITIVE EXP TERM         *
!* XK2    = COEFFICIENT OF THE NEGATIVE EXP TERM         *
!*********************************************************

!======================================================================C

      L=2*NL
 
!     ************MIXED COEFFICENTS**********
!     THIS VERSION AVOIDS SINGULARITIES ASSOC.
!     WITH W0=0 BY SOLVING FOR XK1+XK2, AND XK1-XK2.

      AF(1) = 0.0
      BF(1) = GAMA(1)+1.
      CF(1) = GAMA(1)-1.
      DF(1) = BTOP-CMM1(1)
      N     = 0
      LM2   = L-2

!     EVEN TERMS
 
      DO I=2,LM2,2
        N     = N+1
        AF(I) = (E1(N)+E3(N))*(GAMA(N+1)-1.)       
        BF(I) = (E2(N)+E4(N))*(GAMA(N+1)-1.)
        CF(I) = 2.0*(1.-GAMA(N+1)**2)
        DF(I) = (GAMA(N+1)-1.) * (CPM1(N+1) - CP(N)) +                 &
                  (1.-GAMA(N+1))* (CM(N)-CMM1(N+1))
      END DO
 
      N   = 0
      LM1 = L-1
      DO I=3,LM1,2
        N     = N+1
        AF(I) = 2.0*(1.-GAMA(N)**2)
        BF(I) = (E1(N)-E3(N))*(1.+GAMA(N+1))
        CF(I) = (E1(N)+E3(N))*(GAMA(N+1)-1.)
        DF(I) = E3(N)*(CPM1(N+1) - CP(N)) + E1(N)*(CM(N) - CMM1(N+1))
      END DO
 
      AF(L) = E1(NL)-RSF*E3(NL)
      BF(L) = E2(NL)-RSF*E4(NL)
      CF(L) = 0.0
      DF(L) = BSURF-CP(NL)+RSF*CM(NL)
 
      CALL DTRIDGL(L,AF,BF,CF,DF,XK)
 
!     ***UNMIX THE COEFFICIENTS****

      DO 28 N=1,NL
        XK1(N) = XK(2*N-1)+XK(2*N)
        XK2(N) = XK(2*N-1)-XK(2*N)

!       NOW TEST TO SEE IF XK2 IS REALLY ZERO TO THE LIMIT OF THE
!       MACHINE ACCURACY  = 1 .E -30
!       XK2 IS THE COEFFICEINT OF THE GROWING EXPONENTIAL AND MUST
!       BE TREATED CAREFULLY

        IF(XK2(N) .EQ. 0.0) GO TO 28
        IF (ABS (XK2(N)/(XK(2*N-1)+1.e-20)) .LT. 1.E-30) XK2(N)=0.0 

   28 CONTINUE
 
      RETURN
      end subroutine dsolver
