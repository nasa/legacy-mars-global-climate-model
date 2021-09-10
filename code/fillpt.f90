      subroutine fillpt(pl,p,ptrop,tg,tstrat,tl,plev,tlev,pmid,   &
                        tmid)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

!  Put the T & P GCM arrays onto the NRC grid:  PLEV, PMID, TLEV, TMID
!  Sept 2002
!
!  PMID and TMID are the pressure and temperature at the GCM layer 
!  mid-points.  PLEV and TLEV are the pressures and temperatures at
!  the GCM layer boundaries, i.e. at GCM levels.

      use grid_h
 
      implicit none

      integer :: K, L
      real*8  :: PL(L_LEVELS), PLEV(L_LEVELS), PMID(L_LEVELS)
      real*8  :: TSTRAT, TL(L_LEVELS), TLEV(L_LEVELS), TMID(L_LEVELS)
      real*8  :: P, PTROP, TG

!======================================================================C

!  Fill the new radiation code variables.
!  PLEV and TLEV are the pressure and tempertures on a vertical grid
!  that the new radiation code uses.

      PLEV(1) = ptrop/2.0D0
      PLEV(2) = ptrop/2.0D0

      DO K=3,L_LEVELS
        PLEV(K) = PL(K)
      END DO

      DO K=1,3
        TLEV(K) = TSTRAT
      END DO

      DO K=4,L_LEVELS-1,2
        TLEV(K) = TL(K)
      END DO

      DO K=5,L_LEVELS-2,2
        TLEV(K) = TLEV(K+1) + (TLEV(K-1)-TLEV(K+1))*                   &
                  DLOG(PLEV(K)/PLEV(K+1))/                             &
                  DLOG(PLEV(K-1)/PLEV(K+1))
      END DO

!  Temperature of the bottom level is the ground temperature.

      TLEV(L_LEVELS) = TG

!  Fill the PMID & TMID arrays used by OPTCI and OPTCV subroutines.
!  TMID and PMID used to get the index for CO2 k-coefficient 
!  interpolation.

      TMID(1) = TLEV(2)
      TMID(2) = TLEV(2)
      PMID(1) = PLEV(1)
      PMID(2) = PLEV(2)

      DO L=1,L_LAYERS
        TMID(2*L+1) = TLEV(2*L+1)
        TMID(2*L+2) = TLEV(2*L+1)
        PMID(2*L+1) = PLEV(2*L+1)
        PMID(2*L+2) = PLEV(2*L+1)
      END DO

      TMID(L_LEVELS) = TLEV(L_LEVELS)
      PMID(L_LEVELS) = PLEV(L_LEVELS)

      RETURN
      end subroutine fillpt
