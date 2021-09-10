      subroutine ini_optcld(Qxv,Qxi,Qsv,Qsi,gv,gi,Qextref,TAUREF_CLD)

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

      real*8  Qxv(L_LEVELS+1,L_NSPECTV)
      real*8  Qsv(L_LEVELS+1,L_NSPECTV)
      real*8  gv(L_LEVELS+1,L_NSPECTV)

      real*8  Qxi(L_LEVELS+1,L_NSPECTI)
      real*8  Qsi(L_LEVELS+1,L_NSPECTI)
      real*8  gi(L_LEVELS+1,L_NSPECTI)

      real*8  Qextref(L_LEVELS+1)
      real*8  TAUREF_CLD(L_LEVELS+1)

c  Local variables
c  ---------------

      integer i,k

c Initialyze various variables
c ----------------------------

      DO K = 1, L_LEVELS+1
        Qextref(K)    = 1.
      ENDDO

      call settozero(L_LEVELS+1,TAUREF_CLD)

      call settozero(nlonv*(L_LEVELS+1),Qxv)
      call settozero(nlonv*(L_LEVELS+1),Qsv)
      call settozero(nlonv*(L_LEVELS+1),gv)

      call settozero(nloni*(L_LEVELS+1),Qxi)
      call settozero(nloni*(L_LEVELS+1),Qsi)
      call settozero(nloni*(L_LEVELS+1),gi)

      return

      end
