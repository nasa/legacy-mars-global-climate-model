       subroutine ddevil(p_pbl,ps,Fs,Ts,Tair,FLUX)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

c------------------------------------------------c
c     Subroutine generating dust devils.         c
c Two methods : threshold sensitive(DTH) and not c
c        Details of theory are given in          c
c       Newman et al. (JGR,V107,E12,2002)        c 
c                                                c
c       F. Montmessin (october, 2003)            c
c------------------------------------------------c

      use grid_h
      use constants_h, only: GRAV, RGAS
      use dtcommon_h

      implicit none

c     Arguments
c     ~~~~~~~~~

      real*8 p_pbl              ! Pressure at boundary layer top (pconsav in comp3) in Pascals (mbar/100)
      real*8 ps                 ! Surface pressure in Pascals (mbar/100)
      real*8 Fs                 ! Sensible heat flux
      real*8 Ts                 ! Surface temperature
      real*8 Tair               ! Temperature in the first layer
      real*8 FLUX               ! Lifted mass flux  (kg.m-2.s-1)

c     Local variables
c     ~~~~~~~~~~~~~~~

!     Dust devil lifting efficiency 
      real*8 :: alpha_d = 5.e-11

!     Threshold sensitive formulation
      logical :: DTH = .false.

      real*8 rho                ! near surface atmospheric density
      real*8 b       
      real*8 delta_p            ! Pressure difference across the vortex
      real*8 vtan_t             ! Threshold tangential velocity
      real*8 vtan               ! Actual tangential velocity
      real*8 eta                ! = (1-b), Thermodynamic efficiency of the dust devil
      real*8 eta_H              ! Horizontal thermodynamic efficiency
     
      real*8 :: cpp = 744.5

!     Diameter of lifted dust particles
      real*8 :: Dlift = 6.e-6

      real*8 :: gama = 0.5
  
      real*8 xsi                ! r / cp

c     Treatment
c     ~~~~~~~~~

      vtan = 0.0    ! added 6/26/08  Tim's bug list
      xsi = rgas / cpp + 1.

      FLUX = 0.

      IF (DTH) THEN
 
        rho    = ps / (rgas*Ts)

        vtan_t =  sqrt( 1. + 15. / (dpden_dt*grav*Dlift) ) 
     *           *sqrt( dpden_dt*grav*Dlift/rho )

        b     = (ps**xsi - p_pbl**xsi) / ( (ps-p_pbl) * xsi * ps**xsi )
        eta   = 1. - b
c        eta_H = (Tair-Ts) / Ts
        eta_H = (Ts-Tair) / Tair
        delta_p = ps * ( 1. - 
     *            exp( (gama*eta)/(gama*eta-1.)*(eta_H/(xsi-1.)) ) )
        if (vtan.ge.vtan_t) then
          FLUX = alpha_d * (rho*vtan**2. - 15.) / grav
        endif

      ELSE
     
        IF (Fs.GT.0) then
        b = (ps**xsi - p_pbl**xsi) / ( (ps-p_pbl) * xsi * ps**xsi )
        FLUX = alpha_d * (1.-b) * Fs
        ENDIF
        
      ENDIF

      return
      end
