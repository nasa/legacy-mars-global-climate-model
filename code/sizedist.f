      subroutine sizedist(taustorm,Td,ndp,r,dr,qext,atot,as,
     *                    prho,flux)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

C  Calculate the mass flux from the surface assuming a given source
C  region and particle size distribution.  From the note of 6/11/02.

C  VARIABLES:
C  H            - near surface scale height (R*T/g) (m)
C  ALPHA, GAMMA - Modified GAMMA distribution parameters
C  RM           - Particle size distribution mode radius (m)
C  NDP          - Number of dust particle size bins
C  R(NDP)       - Particle radius of size bin N (m)
C  QEXT(NDP)    - Extinction coefficient of each particle bin
C  TAUSTORM     - Total optical depth produced by the specified surface
C                 source (See note of 6/11/02 for details).
C  NZERO(NDP)   - Surface number density of dust particles in each bin
C                 (particles/m**3/radii-interval)
C  DR(NDP)      - Particle distribution size bin width (m)
C  TAU(NDP)     - Dust opacity of particle bin N
C  MASS(NDP)    - Dust mass loading of the atmosphere for each bin (kg)
C  PRHO         - Dust particle mass density (kg/m**3)
C  Atot         - Total surface area of the planet (m**2)
C  As           - Area of specified surface dust source (m**2)
C  Td           - Time (duration) of dust lifting (See note of 6/11/02
C                 for details of how TAUSTORM, Td and As combine to specify
C                 the amount of dust lifted)
C  FLUX(N)      - Dust flux into the atmosphere (kg/m**2/s) for each
C                 particle bin.

      use constants_h, only: PI

      implicit none

      integer ndp, n
      real*8  c, h
      real*8  taustorm, r(ndp), dr(ndp), prho, fluz(ndp), atot, as
      real*8  nzero(ndp), qext(ndp), sum, ratio, tau(ndp), mass(ndp)
      real*8  flux(ndp), Td

      real*8 :: alpha = 2.0
      real*8 :: gamma = 0.5
      real*8 :: rm    = 4.0e-7

C======================================================================C

      H     = 1.0e+4
      ratio = alpha/gamma

      sum = 0.0
      do n=1,ndp
        sum = sum + (r(n)**alpha)*exp(-ratio*((r(n)/rm)**gamma)) *
     *              pi*r(n)*r(n)*qext(n)*h*dr(n)
      end do
 
      c = taustorm/sum

      do n=1,ndp
        nzero(n) = c*(r(n)**alpha)*exp(-ratio*((r(n)/rm)**gamma))
        tau(n)   = nzero(n)*pi*r(n)*r(n)*qext(n)*h*dr(n)
        mass(n)  = h*nzero(n)*(4.0*pi*r(n)*r(n)*r(n)/3.0)*prho*dr(n)
        flux(n)  = (mass(n)/Td)*(Atot/As)
      end do

      return
      end
