      subroutine radsetup

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

!     PURPOSE:
!        Bundle the new radiation code setup subroutines and call
!     this one subroutine from main, where the three included files
!     are also listed.  Quantities are passed between this driver
!     and the radiation code via modules (eg radcommon_h).
!
!----------------------------------------------------------------------C

      use grid_h
      use radinc_h
      use radcommon_h

      implicit none

      REAL*8  :: VIS2IR, FACTOR
      INTEGER :: NW

!======================================================================C

      call setspv(WNOV,DWNV,WAVEV,SOLARF,TAURAY)
      call setspi(WNOI,DWNI,WAVEI)
      call setrad(TGASREF,PFGASREF,CO2V,CO2I,QEXTV,QSCATV,WV,GV,       &
                        QEXTI,QSCATI,WI,GI,FZEROI,FZEROV)

!  Scale IR opacities (Qexti and Qscati) such that 
!  TAU(0.67 micron)/TAU(9 micron) = VIS2IR, which nominally is 2.
!
!      VIS2IR  = 2.75D0
!
!      factor  = Qextv(L_NREFV)/(VIS2IR*Qexti(L_NREFI))
!
!      DO NW=1,L_NSPECTI
!        Qexti(NW)  = Qexti(NW)*factor
!        Qscati(NW) = Qscati(NW)*factor
!      END DO

      PTOP = 10.0**PFGASREF(1)
      
      return
      end subroutine radsetup
