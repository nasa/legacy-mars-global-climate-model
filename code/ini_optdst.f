      subroutine ini_optdst(QEXTV,QSCATV,GrefV,QEXTI,QSCATI,GrefI,
     *                  Qxv,Qxi,Qsv,Qsi,gv,gi,Qextref)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

      use grid_h
      use radinc_h
      use dtcommon_h
      use cldcommon_h

      implicit none

c  Arguments
c  ---------

      real*8  QEXTV(L_NSPECTV)
      real*8  QSCATV(L_NSPECTV)
      real*8  GrefV(L_NSPECTV)

      real*8  QEXTI(L_NSPECTI)
      real*8  QSCATI(L_NSPECTI)
      real*8  GrefI(L_NSPECTI)

      real*8  Qxv(L_LEVELS+1,L_NSPECTV)
      real*8  Qsv(L_LEVELS+1,L_NSPECTV)
      real*8  gv(L_LEVELS+1,L_NSPECTV)

      real*8  Qxi(L_LEVELS+1,L_NSPECTI)
      real*8  Qsi(L_LEVELS+1,L_NSPECTI)
      real*8  gi(L_LEVELS+1,L_NSPECTI)

      real*8  Qextref(L_LEVELS+1)

c  Local variables
c  ---------------

      integer i,k

c Initialyze various variables
c ----------------------------

      DO K = 1, L_LEVELS+1
        Qextref(K) = QEXTV(L_NREFV)
      ENDDO

      do i = 1, nlonv
        DO K = 1, L_LEVELS+1
          Qxv(K,i) = QEXTV(i)
          Qsv(K,i) = QSCATV(i)
          gv(K,i)  = GrefV(i)
        ENDDO
      enddo
      do i = 1, nloni
        DO K = 1, L_LEVELS+1
          Qxi(K,i) = QEXTI(i)
          Qsi(K,i) = QSCATI(i)
          gi(K,i)  = GrefI(i)
        ENDDO
      enddo

      return

      end
