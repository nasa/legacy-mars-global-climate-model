      subroutine isource(sfcden,threshold,sx,sy,co2ice,fs)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

C  Dust tracer
C  Compute surface dust flux, in a fully interactive manner
C  interactive source
C  April 2002
C
C  SX         - surface stress x-direction  (N/m^2)
C  SY         - surface stress y-direction  (N/m^2)
C  CO2ICE     - Amount of CO2 ice on the ground
C  FS         - mass flux of dust lifted from the surface (kg m^-2 s^-1)
C  THRESHOLD  - stress threshold 
C----------------------------------------------------------------------C

      use constants_h, only: GRAV
      implicit none

      real*8 sx, sy, co2ice, fs
      real*8 stress, sfcden ,threshold, hflux, alfa

!======================================================================!

c      alfa = 1.0e-5  ! Global dust storm
c      alfa = 3.0e-6  ! 24*16 dust interactive, no clouds
      alfa = 5.e-7  ! 24*16 dust interactive, interactive clouds

      stress = sqrt(sx**2 + sy**2)
 
C  No lifting if stress is below the threshold value (THRESHOLD), or if
C  there is CO2 ice on the ground (TG < 0).

      if(stress.le.threshold .or. co2ice.gt.0.0) then
        fs = 0.0
      else
C First the horizontal saltation flux
        hflux = 2.61 * ( stress**1.5/grav/sqrt(sfcden) ) *
     &           ( 1. - sqrt(threshold/stress) ) *
     &           ( 1. + sqrt(threshold/stress) ) ** 2
C Then the saltation flux
        fs = alfa*hflux
      end if

      return
      end
