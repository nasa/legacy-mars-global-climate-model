      subroutine fallvel(scale,dpden,r,t,rho,vf)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

C  Dust tracer 
C  Calculate the fall velocity of dust particles in the Martian
C  atmosphere at the model levels.  Part of the dust tracer scheme for 
C  the c-grid model.
C
C  GCM1.7.2    4/24/02
C  GCM2.0  Sept 2002
C
C  VARIABLES:
C
C  NLEV      - number of levels in the model (GCM: L_LEVELS)
C  T(NLEV)   - temperature at the levels (layer boundaries)
C  RHO(NLEV) - atmosphere density at the layer boundaries
C  VF(NLEV)  - gravitational settling velocity at the layer boundaries
C  Kb        - Boltzmann constant:  Kb  = 1.38066E-23 J K^-1
C  AMU       - atomic mass unit:    Amu = 1.66054E-27 Kg
C  KBAMU     - Kb/AMU
C  SCALE     - 3.0*KBAMU/44.0 (44.0 is the mean molecular weight of
C              the atmosphere)
C  GRAV      - acceleration due to gravity
C  DPDEN     - dust particle density
C  R         - dust particle radius
C  WT        - mean thermal velocity
C  MFP       - mean free path
C  DV        - dynamic viscosity (kg m^-1 s^-1)
C  KN        - Knudsen number
C  ALPHA     - ALPHA: from (1 + ALPHA*Kn) which is the Cunningham
C              slip-flow correction
C  CONST     - the level-independent part of the gravitational settling
C              velocity:  (2*dpden*GRAV*r^2)/9
C
C----------------------------------------------------------------------C

      use constants_h, only: GRAV

      implicit none

      real*8  t, rho, vf, r
      real*8  dpden, wt, mfp, dv, kn, alpha
      real*8  scale, const, kbamu

C======================================================================C

      const = 2.0*dpden*r*r*grav/9.0

      wt    = sqrt(scale*t)
      dv    = (1.59E-6*(t**1.5))/(t+244.4)
      mfp   = 2.0*dv/(rho*wt)
      kn    = mfp/r
      alpha = 1.246 + 0.42*exp(-0.87/kn)
      vf = const*(1.0+alpha*kn)/dv

      return
      end
