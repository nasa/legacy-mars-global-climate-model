      subroutine funcd(astar,downir,rhouch,rhoucht,scond,stemp,sthick,
     *                 tg,f,df,psfg,qmh2o,qice,polarcap)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

      use constants_h, only: stbo, cp
      use cldcommon_h, only: latent_heat

      implicit none 

C  TG is the "X" that we are solving for in the grand scheme of things.
C  It is, of course, the ground temperature variable.

      real*8 tg, f, df

      real*8 astar, scond, stemp, sthick
      real*8 downir, rhoucht, rhouch

!  10-07-11 Latent heating variables

!  psfg    - PSAT
!  qmh2o   - QTRACE(J,I,NLAY,iMa_vap)
!  qice    - QCOND(J,I,iMa_vap)
!  polarcap- true if inside north polar water ice cap; always water
!            ice on the surface

      real*8  :: psfg, qmh2o, qice
      real*8  :: qg
      logical :: polarcap

C======================================================================C

C     F is the function value, df the derivitave

      f  = astar + downir + rhoucht - rhouch*TG +
     *                      2.0*scond*(stemp-TG)/sthick - stbo*tg**4
      df = -rhouch -2.0*scond/sthick - 4.0*stbo*TG**3

      if(LATENT_HEAT) then
        if(qice.gt.0.0 .or. polarcap) then
          qg = (18.0/44.0)*6.11*exp(22.5*(1.0-(273.16/tg)))/psfg
          f  = f + rhouch*(2.8e6/cp)*(qmh2o - qg)
          df = df - 6146.1*rhouch*(2.8e6/cp)*qg/tg**2
        end if
      end if

      return
      end
