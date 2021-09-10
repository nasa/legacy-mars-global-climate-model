      SUBROUTINE GETDETAU(DTDEL,TDEL,TAUCUMIN,WDEL,CDEL,UBAR0,DETAU)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

!======================================================================!

      use grid_h
      use radinc_h

      implicit none

      REAL*8  :: W0(L_NLAYRAD), COSBAR(L_NLAYRAD), DTAU(L_NLAYRAD)
      REAL*8  :: TAU(L_NLEVRAD), WDEL(L_NLAYRAD), CDEL(L_NLAYRAD)
      REAL*8  :: DTDEL(L_NLAYRAD), TDEL(L_NLEVRAD)
      REAL*8  :: FACTOR, TAUCUMIN(L_LEVELS), TAUCUM(L_LEVELS)
      real*8  :: detau

      integer :: NAYER, L, K
      real*8  :: ubar0

!======================================================================C
 
      NAYER  = L_NLAYRAD
      
!  Delta-Eddington Scaling

      FACTOR    = 1.0D0 - WDEL(1)*CDEL(1)**2

      TAU(1)    = TDEL(1)*FACTOR
      TAUCUM(1) = 0.0D0
      TAUCUM(2) = TAUCUMIN(2)*FACTOR
      TAUCUM(3) = TAUCUM(2) +(TAUCUMIN(3)-TAUCUMIN(2))*FACTOR

      DO L=1,L_NLAYRAD-1
        FACTOR      = 1.0D0 - WDEL(L)*CDEL(L)**2
        W0(L)       = WDEL(L)*(1.0D0-CDEL(L)**2)/FACTOR
        COSBAR(L)   = CDEL(L)/(1.0D0+CDEL(L))
        DTAU(L)     = DTDEL(L)*FACTOR
        TAU(L+1)    = TAU(L)+DTAU(L)
        K           = 2*(L+1)
        TAUCUM(K)   = TAU(L+1)
        TAUCUM(K+1) = TAUCUM(K) + (TAUCUMIN(K+1)-TAUCUMIN(K))*FACTOR
      END DO

!  Bottom layer

      L             = L_NLAYRAD
      FACTOR        = 1.0D0 - WDEL(L)*CDEL(L)**2
      W0(L)         = WDEL(L)*(1.0D0-CDEL(L)**2)/FACTOR
      COSBAR(L)     = CDEL(L)/(1.0D0+CDEL(L))
      DTAU(L)       = DTDEL(L)*FACTOR
      TAU(L+1)      = TAU(L)+DTAU(L)
      TAUCUM(2*L+1) = TAU(L+1)
      detau         = TAUCUM(2*L+1)

      return
      end subroutine getdetau
