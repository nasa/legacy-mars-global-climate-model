      SUBROUTINE OPTCI(DTAUI,TAUCUMI,CO2I,PLEV,PFGASREF,TGASREF,       &
                       QREFV,QXIDST,QSIDST,GIDST,COSBI,WBARI,TAUREF,   &
                       TMID,PMID,TAUGSURF,QH2O,WREFH2O,                &
                       Qextrefcld,TAUREFCLD,Qxicld,Qsicld,gicld)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

! THIS SUBROUTINE SETS THE OPTICAL CONSTANTS IN THE INFRARED
! IT CALCUALTES FOR EACH LAYER, FOR EACH SPECRAL INTERVAL IN THE IR
! LAYER: WBAR, DTAU, COSBAR
! LEVEL: TAU
!
! Qrefv is the extinction coefficient at the reference (visible) 
! wavelength - 0.67 microns.
!
! TAUI(L,LW) is the cumulative optical depth at level L (or alternatively
! at the *bottom* of layer L), LW is the spectral wavelength interval.
!
!     TLEV(L) - Temperature at the layer boundary (i.e. level)
!     PLEV(L) - Pressure at the layer boundary (i.e. level)
!     CO2_KI(NT,NP,NW,NG) - IR CO2 k-coefficients 
!                           CO2_K(temp,Pres,Waveln,gauss)
!                           currently: CO2_K(7,11,5,17)
!
!----------------------------------------------------------------------C

      use grid_h
      use constants_h, only: Cmk
      use radinc_h

      implicit none

      real*8  :: DTAUI(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
      real*8  :: DTAUKI(L_LEVELS+1,L_NSPECTI,L_NGAUSS)
      real*8  :: TAUI(L_NLEVRAD,L_NSPECTI,L_NGAUSS)
      real*8  :: TAUCUMI(L_LEVELS,L_NSPECTI,L_NGAUSS)
      real*8  :: TAUGAS
      real*8  :: PLEV(L_LEVELS)
      real*8  :: TMID(L_LEVELS), PMID(L_LEVELS)
      real*8  :: CO2I(L_NTREF,L_PINT,L_REFH2O,L_NSPECTI,L_NGAUSS)
      real*8  :: TGASREF(L_NTREF)
      real*8  :: PFGASREF(L_PINT)
      real*8  :: COSBI(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
      real*8  :: WBARI(L_NLAYRAD,L_NSPECTI,L_NGAUSS)

      integer :: L, NW, NG, K, LK
      integer :: MT(L_LEVELS), MP(L_LEVELS), NP(L_LEVELS)
      real*8  :: ANS, TAUREFL
      real*8  :: TAEROS(L_LEVELS,L_NSPECTI)
      real*8  :: DPR(L_LEVELS), U(L_LEVELS), TAUAC
      real*8  :: LCOEF(4), LKCOEF(L_LEVELS,4)

!     For dust
      real*8  :: QXIDST(L_LEVELS+1,L_NSPECTI)
      real*8  :: QSIDST(L_LEVELS+1,L_NSPECTI)
      real*8  :: GIDST(L_LEVELS+1,L_NSPECTI)
      real*8  :: QREFV(L_LEVELS+1)
      real*8  :: TAUREF(L_LEVELS+1)
      real*8  :: TAUREF_save(L_LEVELS+1)

!     For clouds
      real*8  :: Qxicld(L_LEVELS+1,L_NSPECTI)
      real*8  :: Qsicld(L_LEVELS+1,L_NSPECTI)
      real*8  :: gicld(L_LEVELS+1,L_NSPECTI)
      real*8  :: Qextrefcld(L_LEVELS+1)
      real*8  :: TAUREFCLD(L_LEVELS+1)

      real*8  :: TCLOUD(L_LEVELS,L_NSPECTI),TAUREFLCLD

      real*8  :: TAUREFLK(L_LEVELS+1,L_NSPECTI)
      real*8  :: TAUCLDK(L_LEVELS+1,L_NSPECTI)

! fraction of zeros in each spectral interval, as a function of T, P

      real*8  :: dt, tt
      real*8  :: taugsurf(L_NSPECTI,L_NGAUSS-1)

!  Water mixing ratio variables

      real*8  :: QH2O(L_LEVELS), WREFH2O(L_REFH2O), WRATIO(L_LEVELS)
      real*8  :: KCOEF(4)
      integer :: NH2O(L_LEVELS)

!======================================================================C

      do nw=1,L_NSPECTI
        do ng=1,L_NGAUSS
          DTAUKI(L_LEVELS+1,nw,ng)   = 0.0D0
        end do
        TAUREFLK(L_LEVELS+1,nw) = 0.0D0
        TAUCLDK(L_LEVELS+1,nw)  = 0.0D0
      end do

!  Save old tauref values

      do k=1,L_LEVELS+1
        tauref_save(k) = tauref(k)
      end do

!  Determine the total gas opacity throughout the column, for each
!  spectral interval, NW, and each Gauss point, NG.

      DO NG=1,L_NGAUSS-1
        do NW=1,L_NSPECTI
          TAUGSURF(NW,NG) = 0.0D0
        end do
      end do

      do K=2,L_LEVELS
        DPR(k) = PLEV(K)-PLEV(K-1)
        U(k)   = Cmk*DPR(k)

        call tpindex(PMID(K),TMID(K),QH2O(K),pfgasref,tgasref,WREFH2O, &
                     LCOEF,MT(K),MP(K),NH2O(K),WRATIO(K))

        do LK=1,4
          LKCOEF(K,LK) = LCOEF(LK)
        end do

        TAUREF(K)    = TAUREF(K)    / Qrefv(K)
        TAUREFCLD(K) = TAUREFCLD(K) / Qextrefcld(K)

        DO NW=1,L_NSPECTI
          TAEROS(K,NW) = TAUREF(K)    * Qxidst(K,NW)
          TCLOUD(K,NW) = TAUREFCLD(K) * Qxicld(K,NW)
        END DO
      end do

      do K=2,L_LEVELS
        do nw=1,L_NSPECTI
          do ng=1,L_NGAUSS-1

!           NOW COMPUTE TAUGAS

!  Interpolate between water mixing ratios
!  WRATIO = 0.0 if the requested water amount is equal to, or outside the
!  the range of water amount data.

            KCOEF(1) = CO2I(MT(K),MP(K),NH2O(K),NW,NG) + WRATIO(K)*    &
                      (CO2I(MT(K),MP(K),NH2O(K)+1,NW,NG) -             &
                       CO2I(MT(K),MP(K),NH2O(K),NW,NG))

            KCOEF(2) = CO2I(MT(K),MP(K)+1,NH2O(K),NW,NG) + WRATIO(K)*  &
                      (CO2I(MT(K),MP(K)+1,NH2O(K)+1,NW,NG) -           &
                       CO2I(MT(K),MP(K)+1,NH2O(K),NW,NG))

            KCOEF(3) = CO2I(MT(K)+1,MP(K)+1,NH2O(K),NW,NG) + WRATIO(K)*&
                      (CO2I(MT(K)+1,MP(K)+1,NH2O(K)+1,NW,NG) -         &
                       CO2I(MT(K)+1,MP(K)+1,NH2O(K),NW,NG))

            KCOEF(4) = CO2I(MT(K)+1,MP(K),NH2O(K),NW,NG) + WRATIO(K)*  &
                      (CO2I(MT(K)+1,MP(K),NH2O(K)+1,NW,NG) -           &
                       CO2I(MT(K)+1,MP(K),NH2O(K),NW,NG))

!  Interpolate the CO2 k-coefficients to the requested T,P


            ANS = LKCOEF(K,1)*KCOEF(1) + LKCOEF(K,2)*KCOEF(2) +        &
                  LKCOEF(K,3)*KCOEF(3) + LKCOEF(K,4)*KCOEF(4)

            TAUGAS          = U(k)*ANS
            TAUGSURF(NW,NG) = TAUGSURF(NW,NG) + TAUGAS
            DTAUKI(K,nw,ng) = TAUGAS+TAEROS(K,NW)+TCLOUD(K,NW)
          end do

!  Now fill in the "clear" part of the spectrum (NG = L_NGAUSS)
!  Which holds continuum opacity only

          NG              = L_NGAUSS
          DTAUKI(K,nw,ng) = TAEROS(K,NW)+TCLOUD(K,NW)
        end do
      end do

!  Now the full treatment for the layers, where besides the opacity
!  we need to calculate the scattering albedo and asymmetry factors
!  for each layer

      DO NW=1,L_NSPECTI
        DO K=2,L_LEVELS+1
          TAUREFLK(K,NW) = TAUREF(K)    * QSIDST(K,NW)
          TAUCLDK(K,NW)  = TAUREFCLD(K) * QSICLD(K,NW)
        ENDDO
      ENDDO

      DO NW=1,L_NSPECTI

!  First, the special "clear" channel

        NG = L_NGAUSS

        DO L=1,L_NLAYRAD
          K              = 2*L+1
          DTAUI(L,nw,ng) = DTAUKI(K,NW,NG) + DTAUKI(K+1,NW,NG) + 1.d-50
          if(DTAUI(L,NW,NG) .GT. 1.0E-9) then
            WBARI(L,nw,ng) = (TAUREFLK(K,NW)+ TAUREFLK(K+1,NW) +       &
                              TAUCLDK(K,NW) + TAUCLDK(K+1,NW)   ) /    &
                              DTAUI(L,NW,NG)
          else
            WBARI(L,nw,ng) = 0.0D0
            DTAUI(L,NW,NG) = 1.0E-9
          endif

          TAUAC = TAUREFLK(K,NW)+ TAUREFLK(K+1,NW) + TAUCLDK(K,NW) +   &
                  TAUCLDK(K+1,NW)

          if(TAUAC .GT. 0.0) then
            cosbi(L,NW,NG) = ( GIDST(K,NW)   * TAUREFLK(K,NW) +        &
                               GIDST(K+1,NW) * TAUREFLK(K+1,NW) +      &
                               GICLD(K,NW)   * TAUCLDK(K,NW) +         &
                               GICLD(K+1,NW) * TAUCLDK(K+1,NW) ) /     &
                              (TAUREFLK(K,NW)+ TAUREFLK(K+1,NW) +      &
                               TAUCLDK(K,NW) + TAUCLDK(K+1,NW)   )
          else
            cosbi(L,NW,NG) = 0.0D0
          end if
        END DO

!  . . .Now the other Gauss points, if needed.

        DO NG=1,L_NGAUSS-1

          DO L=1,L_NLAYRAD
            K              = 2*L+1
            DTAUI(L,nw,ng) = DTAUKI(K,NW,NG)+DTAUKI(K+1,NW,NG)+1.d-50
            if(DTAUI(L,NW,NG) .GT. 1.0E-9) then
              WBARI(L,nw,ng) = (TAUREFLK(K,NW)+ TAUREFLK(K+1,NW) +     &
                                TAUCLDK(K,NW) + TAUCLDK(K+1,NW)   ) /  &
                                DTAUI(L,NW,NG)
            else
              WBARI(L,nw,ng) = 0.0D0
              DTAUI(L,NW,NG) = 1.0E-9
            endif

            cosbi(L,NW,NG) = cosbi(L,NW,L_NGAUSS)
          END DO
        END DO

      END DO     ! NW spectral loop

!     TOTAL EXTINCTION OPTICAL DEPTHS

      DO NW=1,L_NSPECTI
        NG = L_NGAUSS
        TAUI(1,NW,NG) = 0.0D0
        DO L=1,L_NLAYRAD
          TAUI(L+1,NW,NG) = TAUI(L,NW,NG)+DTAUI(L,NW,NG)
        END DO

        TAUCUMI(1,NW,NG)=0.0D0
        DO K=2,L_LEVELS
          TAUCUMI(K,NW,NG)=TAUCUMI(K-1,NW,NG)+DTAUKI(K,NW,NG)
        END DO

        DO NG=1,L_NGAUSS-1
          TAUI(1,NW,NG)=0.0D0
          DO L=1,L_NLAYRAD
            TAUI(L+1,NW,NG)=TAUI(L,NW,NG)+DTAUI(L,NW,NG)
          END DO

          TAUCUMI(1,NW,NG)=0.0D0
          DO K=2,L_LEVELS
            TAUCUMI(K,NW,NG)=TAUCUMI(K-1,NW,NG)+DTAUKI(K,NW,NG)
          END DO
        END DO
      END DO

!  Restore old tauref values

      do k=1,L_LEVELS+1
        tauref(k) = tauref_save(k)
      end do

      RETURN
      end subroutine optci
