      subroutine eddycoef(nstep,ric,vk,alphl0,epsl0,dxi,q0,
     *                    du,I1D,J1D,
     *                    xinum,z,rkm,rkh,psi,dt)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

c
c A new scheme for determining the eddy mixing coefficients
c 2 (count 'em) fiddles...errr...parameterizations are present:
c dudvdz smoothed in time to remove oscillations in lowest model level
c Kh and Km have minimum values: From Arya, Into to Micromet, p164-166
c
c Apart from all this the scheme is completely kosher
c
      use grid_h
      use constants_h, only: GRAV
      use pbl_defines_h

      implicit none

c     global variables:

      real*8 dxi(2*n+1),q0(2*n+1),du(l_jsize,l_isize,2*n+1),
     *          xinum(2*n+1),z(2*n+1),rkm(2*n+1),rkh(2*n+1)
      real*8 psi(2*n+1,nvar)

!     implicit none

      integer :: i1d, j1d, nstep, k
      real*8  :: ric, vk, alphl0, epsl0, rl0, dudvdz, ri, dt, rim, rkh0
      real*8  :: rkm0, rkmin, rl, beta, dudz, dvdz, dthdz

C======================================================================C

C Specify constants

      rl0    = 150.0

      do 100 k = 3, 2*n-1, 2

c Calculate mixing length and beta

        rl        = rl0 * vk * z(k) / (rl0 + vk * z(k))
        beta      = (dxi(k-1) + dxi(k+1)) / (dxi(k-1) *
     *               psi(k-1,3) + dxi(k+1) * psi(k+1,3))

c Calculate shears and Richardson number

        dudz      = (psi(k+1,1) - psi(k-1,1)) * xinum(k)
        dvdz      = (psi(k+1,2) - psi(k-1,2)) * xinum(k)
        dthdz     = (psi(k+1,3) - psi(k-1,3)) * xinum(k)
        dudvdz    = dudz*dudz + dvdz*dvdz
        ri        = beta*grav*dthdz / (dudvdz+1.0E-9)

        du(j1d,i1d,k) = du(j1d,i1d,k)-(du(j1d,i1d,k)-dudvdz)*dt/1.0E4
        rim       = beta*grav*dthdz / (du(j1d,i1d,k)+1.0E-9)

        q0(k)     = rim

c Neutrally stable mixing coefficients = q(Ri=0)*l*l*S(Ri=0)

        rkh0      = sqrt(du(j1d,i1d,k)/0.153) * rl * rl * 0.493
        rkm0      = sqrt(du(j1d,i1d,k)/0.153) * rl * rl * 0.393

c Variation of Km and Kh with Ri (from Arya p.164)

        if(rim .le. 0.0)then

          rkh(k)  = rkh0*((1.0-15.0*rim)**0.50)
          rkm(k)  = rkm0*((1.0-15.0*rim)**0.25)

        else

          rkh(k)  = rkh0*(1.0-rim/ric)
          rkm(k)  = rkm0*(1.0-rim/ric)

        endif

C set limiting value to RKM

        rkmin = 0.001
        if(z(k) .lt. 300.0)rkmin=0.1

        rkh(k)  = max(rkh(k),rkmin)
        rkm(k)  = max(rkm(k),rkmin)

  100 continue

      return
      end
