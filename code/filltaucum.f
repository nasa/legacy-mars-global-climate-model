      SUBROUTINE FILLTAUCUM(JFF,IFF)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

C
C  PURPOSE
C      FILLTAUCUM fills the TAUCUM array - the dust column density of each
C      GCM level.  NDUSK(K) is the column density of the sublayer whose
C      bottom is level K, and whose top is at K-1.
C
      use grid_h
      use defines_h
      use radinc_h
      use radcommon_h
      use standard_h
      use comp3cmn_h
      use dtcommon_h

      implicit none

      external JSRCHGT

      INTEGER Nstar, N

!     implicit none

      integer :: JFF, J, IFF, I, K, JSRCHGT
      real*8  :: PSTAR, PSTAR1
 
C#=====================================================================

      J = JFF
      I = IFF

C  Dust optical depth at the bottom of each sub-layer.

      DO N=1,3
        TAUCUM(N) = 0.0
        TAUREF(N) = 0.0
      END DO

      DO K=4,L_LEVELS
        PSTAR     = P(J,I)*SIGMA(K)+PTROP
        PSTAR1    = MAX(PSTAR,PRDST(1))
        NSTAR     = MIN0(JSRCHGT(NPDST-1,PRDST,1,PSTAR1)-1,NPDST-1)
        TAUCUM(K) = TAUDST(J,I,NSTAR)+(PSTAR1-PRDST(NSTAR))*
     *              (TAUDST(J,I,NSTAR+1) - TAUDST(J,I,NSTAR))/
     *              (PRDST(NSTAR+1)-PRDST(NSTAR))
        TAUREF(K) = TAUCUM(K) - TAUCUM(K-1)
        TAUREF3D(J,I,K) = TAUREF(K)
      END DO

      TAUREF(L_LEVELS+1) = 0.0

      RETURN
      END
