      subroutine microphys(dt,pl,tl,kd,sx,sy,Fs,p_pbl,tsurf,co2ice,
     *                     dustref,qpi,qpig,qextrefdst)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

      use grid_h
      use constants_h, only: PI, GRAV, RGAS
      use dtcommon_h
      use comp3cmn_h, only: jcmn,icmn
      use cldcommon_h
      use standard_h, only: srfupflx,srfdnflx,npcflag

      implicit none

c  Arguments
c  ---------

      REAL*8 dt              ! physical time step (s)
      REAL*8 pl(l_levels)    ! pressure at each level (mbar)
      REAL*8 tl(l_levels)    ! temperature at each level (K)
      REAL*8 kd(2*l_layers+1)! eddy mixing coefficient (m2/s)
      REAL*8 sx,sy           ! Friction velocity
      REAL*8 Fs              ! Sensible heat flux exchanged with the surface
      REAL*8 p_pbl           ! Pressure at pbl top
      REAL*8 tsurf           ! ground temperature (equivalent to GT)
      REAL*8 co2ice          ! Amount of CO2 ice on the surface     
      REAL*8 dustref         ! cumulated dust opacity used in radiative transfer
      real*8 qextrefdst(L_LEVELS+1)

c    Tracers :
      REAL*8 qpi(l_levels,ntrace)  ! tracer (kg/kg)
      REAL*8 qpig(ntrace)          ! tracer on surface (kg/m2)


c  Local variables
c  ---------------

      logical add_dt             ! true if some dust has to be be added in the first layer
 
      real*8 m_aer(nz,nbin,ntype)
      real*8 n_aer(nz,nbin,ntype-1)
      real*8 h2ovap(nz)
      real*8 rho(nz)

      real*8 dev,dev2,cst
      real*8 Mo,No,Rn,Rm
      real*8 deposit(ntrace)
      real*8 Rcor(nz,naer)
      real*8 dens(nz,naer)
      real*8   g

      integer i,j,l,n   ! Loop integers
      integer ilay,iq


c debug variables
c ---------------
      real*8 qsat,test(10)
      integer k

c  Treatment 
c  ---------

      g  = grav

c     This is where dust is injected in the atmsophere

      call dust_update(dt,pl,tl,tsurf,co2ice,sx,sy,Fs,p_pbl,
     *                 dustref,qpi,qpig,add_dt,qextrefdst)

c     Compute atmospheric density
      do l = 1, nz
        ilay   = l * 2 + 2
        rho(l) = pl(ilay) * 100.0 / ( rgas*tl(ilay) )
      enddo

      do l = 1, nz
        ilay   = l * 2 + 2
        do n = 1, naer
        qpi(ilay,n) = max(qpi(ilay,n),0.)
        select case(n)
        case (1)
        Mo = qpi(ilay,iMa_dt)
        No = qpi(ilay,iNb_dt) + 1.d-50
        dens(l,n) = aerdens(iMa_dt)
        Rcor(l,n) = ( Mo / No * 0.75 / pi / dens(l,n) )**(athird)
     *              * exp( 3.* stdv(n)**2. )
        if (Mo.lt.1.e-20) Rcor(l,n)=1.e-8
        case (2)
        dens(l,n) = aerdens(iMa_dt)
        Rcor(l,n) = Rcor(l,iMa_dt) * exp( -3.* stdv(n)**2. )
        case (3)
        Mo = qpi(ilay,iMa_cld) + qpi(ilay,iMa_cor) + 1.d-50
        No = qpi(ilay,iNb_cld) + 1.d-50
        dens(l,n) = qpi(ilay,iMa_cld) / Mo * aerdens(iMa_cld)
     *             +qpi(ilay,iMa_cor) / Mo * aerdens(iMa_dt)
        dens(l,n) = min(max(dens(l,n),dpden_ice),dpden_dt)
        Rcor(l,n) = ( Mo / No * 0.75 / pi / dens(l,n) )**(athird)
     *              * exp( 3.* stdv(n)**2. )
        if (Mo.lt.1.e-20) Rcor(l,n)=1.e-8
        case (4)
        dens(l,n) = dens(l,iMa_cld)
        Rcor(l,n) = Rcor(l,iMa_cld) * exp( -3.* stdv(n)**2. )
        case (5)
        dens(l,n) = dens(l,iMa_cld)
        Rcor(l,n) = Rcor(l,iMa_cld)
        end select
        enddo
      enddo

c  Now compute the microphysical processes.
c  Always check the order of sedim and nucleacond.
c  newpbl called before microphys, so do sedim first !

      call sedim(dt,pl,tl,rho,kd,add_dt,qpi,Rcor,dens,deposit)

      if (h2ocloudform .eqv. .true.) then
       call nucleacond(dt,pl,tl,rho,qpi)
      endif

c     Update water ice on the surface
      do i = 1, ntrace
        qpig(i) = qpig(i) + deposit(i)
      enddo

      return
      END

********************************************************************
*							           *
      subroutine dust_update(dt,pl,tl,tsurf,co2ice,sx,sy,Fs,p_pbl,
     *                       dustref,qpi,qpig,add_dt,qextrefdst)
*      							           *
*                Computing dust source                             *
*							           *
********************************************************************

      use grid_h
      use constants_h, only: PI, GRAV, RGAS
      use dtcommon_h
      use comp3cmn_h, only: jcmn,icmn
      use cldcommon_h
      use standard_h, only: srfupflx,srfdnflx,npcflag

      implicit none

c Arguments
c ---------

      real*8   dt
      real*8 pl(l_levels)
      real*8 tl(l_levels)
      real*8 tsurf
      real*8 co2ice
      real*8 dustref
      real*8 sx,sy           ! Friction velocity
      real*8 Fs,p_pbl
      real*8 qpi(l_levels,ntrace)
      real*8 qpig(ntrace)
      real*8 qextrefdst(l_levels+1)

      logical add_dt

c Local
c -----

      integer l,ilay

      real*8 Madd,Nadd
      real*8 dev,dens
      real*8 cst,cst2
      real*8 tau_ref,tau_add,tau_current
      real*8 Reff,Rn,Rs
      real*8 ym

      logical :: Newman     = .false.
      logical :: devil      = .false.
      logical :: background = .true.
      logical :: kmh        = .false.
      logical :: MGS        = .false.
      save MGS,Newman,background
      real*8 taueq,tauN,tauS,latpoint

      real*8 Reff_lift,Rm,mass_part
      real*8 FLUX,threshold,stmag
      real*8 rho,alfa

!     corresponds to injecting an opacity of 1 per day
      real*8 :: fluxmax = 2.e-8
    
      integer l_pbl

      real*8 defla ! local deflation in mm 

c Treatment
c ---------

      add_dt  = .false.
      Nadd    = 0. 
      Madd    = 0.
      ym   = 100. * ( pl(l_levels)-pl(l_levels-2) ) / grav
      if(ym.lt.0.0) ym = 1.0D-200

        if (Newman) then
c         Dust Flux at the base generated by saltation
          rho  = 100.* pl(l_levels) / (rgas*tsurf)
!         ym   = 100. * ( pl(l_levels)-pl(l_levels-2) ) / grav
c          defla = -(qpig(iMa_dt)+qpig(iMa_cor)) / 2.5
          FLUX  = 0.
          call isource(rho,0.0225,sx,sy,co2ice,FLUX)
          FLUX = max(FLUX*dt,0.)
          if (FLUX.gt.0) then
            Reff_lift  = 2.e-6
            Rm  = Reff_lift * exp( -0.5 * dev_dt**2. )
            mass_part = 4./3. * pi * 2500. * Rm**3.
            Madd  = FLUX / ym
            Nadd  = Madd / mass_part
            qpi(2*nz+2,iMa_dt) = qpi(2*nz+2,iMa_dt) + Madd
            qpi(2*nz+2,iNb_dt) = qpi(2*nz+2,iNb_dt) + Nadd
            qpig(iMa_dt) = qpig(iMa_dt) - FLUX
          endif
        endif  ! end Newman
        if (kmh) then
           FLUX = 0.
           alfa = 0.03
           threshold=0.0225
           stmag=sqrt(sx**2+sy**2)
           if (stmag .le. threshold .or. co2ice .gt. 0.0) then
            flux =0.
           else
            flux=alfa*2.43E-3*(stmag**2)*(stmag-threshold)/threshold
           endif
           flux = max(flux*dt,0.)
           if (flux .gt. 0.) then
              Reff_lift = 2.5E-6
              Rm = Reff_lift * exp(-0.5*dev_dt**2.)
              mass_part = 4./3.* pi * 2500. * Rm**3.
              Madd = flux / ym
              Nadd = Madd / mass_part
              qpi(2*nz+2,ima_dt)=qpi(2*nz+2,ima_dt) + Madd
              qpi(2*nz+2,inb_dt)=qpi(2*nz+2,inb_dt) + Nadd
              qpig(ima_dt) = qpig(ima_dt)-flux
           endif
        srfupflx(jcmn,icmn,ima_dt) = flux/dt
!        if(flux .gt. 0.) print*,srfupflx(jcmn,icmn,ima_dt)
        endif ! end kmh


c       Dust Flux put in pbl by dust devil
        if (devil.and.p_pbl.ne.0) then
          call ddevil(p_pbl*100.,pl(l_levels)*100.,Fs,tsurf,tl(2*nz+2)
     *                ,FLUX)
          FLUX = max(FLUX*dt,0.)
          if (FLUX.gt.0) then
            l_pbl = 2*nz+2  ! Let dust be injected only in the 1st layer
            do l = 2*nz+2, 4, -2
              if (pl(l).ge.p_pbl) l_pbl = l
            enddo
!           ym = 100. * ( pl(l_levels)-pl(l_pbl-1) ) / grav
            Reff_lift = 2.e-6
            Rm  = Reff_lift * exp( -0.5 * dev_dt**2. )
            mass_part = 4./3. * pi * 2500. * Rm**3.
            Madd = FLUX / ym
            Nadd = Madd / mass_part
            do l = 2*nz+2, l_pbl, -2
              qpi(l,iMa_dt) = qpi(l,iMa_dt) + Madd
              qpi(l,iNb_dt) = qpi(l,iNb_dt) + Nadd
            enddo
            qpig(iMa_dt) = qpig(iMa_dt) - FLUX
          endif
        endif  ! end Devil

        if (Background) then

        dens    = dpden_dt
        dev     = dev_dt
        cst     = .75 / (pi*dens) * exp( -4.5*dev**2. )
c        tau_ref = 0.3
        tau_ref = dustref
        tau_current = 0.

        do l = 1, nz
          ilay = 2 * l + 2
          cst2 = 100. * (pl(2*l+3) - pl(2*l+1)) / grav
     *           *  qextrefdst(ilay)
c         Compute opacity of dust alone
          Rn   = ( qpi(ilay,iMa_dt)/( qpi(ilay,iNb_dt)+1.e-30) 
     *           *cst )**(athird)
          Rs   = Rn * dexp( dev**2. )
          tau_current = tau_current + pi * Rs**2. * qpi(ilay,iNb_dt)
     *                   * cst2

        enddo
c       Scale current opacity to a reference pressure of 6.1 mbar
        tau_current = tau_current / pl(l_levels) * 6.1
        tau_add = tau_ref - tau_current

        if ( tau_add .gt. 0. ) then
!         if (.not. npcflag(jcmn,icmn)) then
          add_dt = .true.
!         ym     = 100. * ( pl(l_levels)-pl(l_levels-2) ) / grav
!          Reff   = 2.5e-6
          Reff   = 2.0e-6
          Rn     = Reff * dexp( -2.5 * dev**2. )
          Rs     = Rn * dexp( dev**2. )

          Nadd   = 0.2 / (Qext_dt(1) * pi * Rs**2. ) / ym 
          Madd   = 4./3.*pi * Nadd * dens * Rn**3. * dexp( 4.5*dev**2. )
          Nadd   = Nadd / 88775. * dt
          Madd   = Madd / 88775. * dt
!         endif
        endif
    
        qpi(2*nz+2,iMa_dt) = qpi(2*nz+2,iMa_dt) + Madd
        qpi(2*nz+2,iNb_dt) = qpi(2*nz+2,iNb_dt) + Nadd

        FLUX         = Madd * ym
        qpig(iMa_dt) = qpig(iMa_dt) - FLUX
        srfupflx(jcmn,icmn,ima_dt) = flux/dt
        endif ! End of background test

      return
      end
   

****************************************************************
*							       *
      subroutine sedim(dt,pl,tl,rho,kd,add_dt,qpi,Rcor,dens,deposit)
*							       *
*                Computing aerosol sedimentation               *
*							       *
****************************************************************

      use grid_h
      use constants_h, only: GRAV, RGAS
      use dtcommon_h
      use comp3cmn_h, only: jcmn,icmn
      use cldcommon_h
      use standard_h, only: srfupflx,srfdnflx

      implicit none

c Arguments
c ---------
 
      integer ndim           ! dimension of the c array

      real*8 dt              ! comp3 time step
      real*8 pl(2*nz+3),tl(2*nz+3)
      real*8 rho(nz)         ! atmospheric density
      real*8 kd(2*nz+1)        ! Eddy mixing coefficient (m2/s)
      real*8 qpi(l_levels,ntrace)! Integrated concentrations of particle/each bin (#/m3)
      real*8 deposit(ntrace)         ! amount of water ice (kg/m2) felt on the ground during dt
      real*8 Rcor(nz,naer)
      real*8 dens(nz,naer)      ! density of aerosol (kg/m3)

      logical add_dt

c Local variables
c ---------------

      real*8 rhob(nz)       ! atmospheric density at layer boundaries
      real*8 dz(nz)         ! layer thickness (m)
      real*8 vf             ! fall velocity of particle (m/s)
      real*8 cour,w1        ! Courant number & corrected velocity (vf+kd/H with h scale height)
      real*8 dzbx
      real*8 sigma,theta,hc,lg,rap,cmp,w,wp 
      real*8 fs(nz+1),ft(nz+1)
      real*8 as(nz),bs(nz),cs(nz),ds(nz)
      real*8 asi(nz),bsi(nz),csi(nz),dsi(nz),xsol(nz)
      real*8 cold(2)
      real*8 c(nz)

      integer i,l,n
      integer ilay,ilev

c debug variables
c ---------------
      real*8 test(10)

      theta = 0.0    ! added 6/26/08  Tim's bug list

c     Water ice deposit reset to zero
      do i=1,ntrace
        deposit(i) = 0.
      enddo

c     Layer thickness: dz
      do l = 1, nz
        dz(l) = 100. * ( pl(2*l+3)-pl(2*l+1) ) / grav / rho(l)
      enddo

c     Compute density at the layer boundaries
      do l = 1, nz
        ilev = l * 2 + 3
        if (l.lt.nz) then
          rhob(l) = pl(ilev)  * 100.0 / (rgas*tl(ilev))
        else
          rhob(l) = pl(ilev-1)* 100.0 / (rgas*tl(ilev-1))
        endif
      enddo

c Loop over type of aerosols (1 for dust, 2 for water ice, and 3 for dust cores 
c                             if the c array is the mass distribution)
      DO 1 n = 1, naer

c Loop over layers

      do 20 l = 1, nz

      ilay = 2 * l + 2

      c(l) = qpi(ilay,n) * rho(l)

      if (l.eq.1) goto 20
 
c     Compute fall velocity
      ilev = 2 * l + 3
      call fallvel(scale_dt,dens(l,n),
     &             Rcor(l,n),tl(ilev),rhob(l),vf)

      dzbX = ( dz(l)+dz(l-1) ) / 2.

      w  = -1. * vf * exp(-stdv(n)**2.)

c     Get the corrected fall velocity (virtual speed accounting for mixing)

      if (kd(2*l-1) .ne. 0.) then
        theta = 0.5 * ( w*dzbX/kd(2*l-1) + log(rho(l-1)/rho(l)) )
        if (theta.ne.0) then
          sigma = 1./dtanh(theta) - 1./theta
        else
          sigma = 1.
        endif
      else
        sigma = 1.
      endif

      if (c(l).eq.0.) then
        rap=10.
        if (c(l-1).eq.0.) then
          rap=1.
        endif
      else
        rap = min( max(c(l-1)/c(l),0.1), 10.)
      endif

      cour=abs(w*dt)

      if (rap.gt.0.9 .and. rap.lt.1.1 .or. cour.gt.dz(l)) then
        w1 = w
      else
        if (w.lt.0) then
          hc = dzbX / dlog(rap)
          lg = dzbX / (w*dt) * (dexp(-w*dt/hc)-1.) / (1.-rap)
          wp = w * 1.d0
          cmp= dlog(-wp) + abs(sigma) * dlog(lg)
          w1 = -dexp(cmp)
        else
          w1 = 0.
        endif
      endif

c  Fluxes at layer boundaries

      if (kd(2*l-1).ne.0.) then
        if (theta.ne.0.) then
          ft(l)=( w1 + log(rho(l-1)/rho(l))*kd(2*l-1)/dzbX ) 
     &          / ( dexp(2.*theta) - 1. )
          fs(l) = ft(l) * dexp(2.*theta)
        else
          ft(l) = kd(2*l-1) / dzbX
          fs(l) = kd(2*l-1) / dzbX
        endif
      else
        if (w1.lt.0.)then
          ft(l) = -w1
          fs(l) = 0.
        else
          ft(l) = 0.
          fs(l) = w1
        endif
      endif

20    continue

c Boundary conditions for the fluxes

      fs(1)    =  0.
      ft(1)    =  0.
      fs(nz+1) =  0.
      ft(nz+1) = -w1

      if (add_dt .and. n.eq.iMa_dt) ft(nz+1) = 0.
      if (add_dt .and. n.eq.iNb_dt) ft(nz+1) = 0.

c Compute the coefficient of the continuity equation

      do l=1,nz
        cs(l) =  ft(l+1) + fs(l) - dz(l) / dt
        if ( cs(l) .gt. 0. ) goto 1000 
        as(l) = -dz(l) / dt
        bs(l) = -ft(l)
        ds(l) = -fs(l+1)
      enddo

c Depending on the cs value, switch to an explicit or an implicit scheme

c Explicit case 

      cold(1)  = c(1)
      c(1) = ( cs(1)*c(1) + ds(1)*c(2) ) / as(1)

      do l = 2, nz-1
        cold(2)  = c(l)
        c(l) = ( bs(l)*cold(1) + cs(l)*c(l)
     &             + ds(l)*c(l+1) ) / as(l)
        cold(1)  = cold(2)
      enddo 

c Compute the mass of water ice falling on the ground
      if (n.eq.iMa_cld) deposit(iMa_vap) = deposit(iMa_vap)
     *                                    + c(nz) * ft(nz+1) * dt
      if (n.eq.iMa_dt)  deposit(iMa_dt)   = deposit(iMa_dt) 
     *                                    + c(nz) * ft(nz+1) * dt
      if (n.eq.iMa_cor) deposit(iMa_cor)  = deposit(iMa_cor)
     *                                    + c(nz) * ft(nz+1) * dt

      if (n.eq.iMa_dt)  then

        srfdnflx(jcmn,icmn,iMa_dt) =
     *                  srfdnflx(jcmn,icmn,ima_dt) +
     *                  c(nz) * ft(nz+1)
      endif
      if (n.eq.iMa_cor)  srfdnflx(jcmn,icmn,iMa_cor) =
     *                  srfdnflx(jcmn,icmn,ima_cor) +
     *                                  c(nz) * ft(nz+1)
      if (n.eq.iMa_cld)  srfdnflx(jcmn,icmn,iMa_cld) =
     *                  srfdnflx(jcmn,icmn,ima_cld) +
     *                                  c(nz) * ft(nz+1)

      c(nz) = ( bs(nz)*cold(1) + cs(nz)*c(nz) ) / as(nz)

      do l = 1, nz
        qpi(2*l+2,n) = c(l) / rho(l)
      enddo

      GOTO 1

1000  continue

c Implicit case 

      do l = 1, nz
        asi(l) =  ft(l)
        bsi(l) = -( ft(l+1) + fs(l) + dz(l)/dt )
        csi(l) =  fs(l+1)
        dsi(l) = -dz(l) / dt * c(l)
      enddo

c Matrix inversion

      call dtridgl(nz,asi,bsi,csi,dsi,xsol)

      do l = 1, nz
        c(l) = xsol(l)
        qpi(2*l+2,n) = c(l) / rho(l)
      enddo

c Compute the mass of water ice falling on the ground
      if (n.eq.iMa_cld) deposit(iMa_vap) = deposit(iMa_vap) 
     *                                    + c(nz) * ft(nz+1) * dt
      if (n.eq.iMa_dt)  deposit(iMa_dt)   = deposit(iMa_dt)  
     *                                    + c(nz) * ft(nz+1) * dt
      if (n.eq.iMa_cor) deposit(iMa_cor)  = deposit(iMa_cor) 
     *                                    + c(nz) * ft(nz+1) * dt

      if (n.eq.iMa_dt) then
        srfdnflx(jcmn,icmn,iMa_dt) =
     *                  srfdnflx(jcmn,icmn,ima_dt) +
     *                                  c(nz) * ft(nz+1)

      endif
      if (n.eq.iMa_cor)  srfdnflx(jcmn,icmn,iMa_cor) =
     *                  srfdnflx(jcmn,icmn,ima_cor) +
     *                                  c(nz) * ft(nz+1)
      if (n.eq.iMa_cld)  srfdnflx(jcmn,icmn,iMa_cld) =
     *                  srfdnflx(jcmn,icmn,ima_cld) +
     *                                  c(nz) * ft(nz+1)

1     CONTINUE

      RETURN
      END

****************************************************************
      subroutine nucleacond(dt,pl,tl,rho,qpi)
*                                                              *
*     This routine updates species concentrations due          *
*     to both nucleation and condensation-induced variations.  *
*     Gain and loss rates associated to each one of these      *
*     processes are computed separately in other routines.     *
*                                                              *
****************************************************************

      use grid_h
      use constants_h, only: PI
      use dtcommon_h
      use comp3cmn_h, only: jcmn,icmn
      use cldcommon_h
 
      implicit none

c  Arguments
c  ---------

      real*8 dt                    ! comp3 time step
      real*8 pl(2*nz+3),tl(2*nz+3)
      real*8 rho(nz)               ! Atmospheric density (kg/m3)
      real*8 qpi(2*nz+3,ntrace)

c  Local
c  -----

      integer i,l

      real*8 n_aer(nbin)   ! number concentrations of particle/each size bin
      real*8 m_aer(nbin)   ! number concentrations of particle/each size bin
      real*8 sat_ratio     ! Water vapor saturation ratio over ice
      real*8 ph2o          ! Water vapor partial pressure (Pa) 
      real*8 qsat          ! Water vapor mass mixing ratio at saturation (kg/kg of air)
      real*8 rate(nbin)    ! Nucleation rate (s-1)
      real*8 h2ovap        ! Water vapor mass mixing ratio (kg/m3)
       
      real*8 qpisav(naer)
      real*8   Cste
      real*8 p,temp
      real*8 Mo,No,dens
      real*8 up,dwn,Ctot,gr,rad,seq
      real*8 newvap,dN,dM
      real*8 Rn,Rm,dev2
      real*8   sig

      real*8 newT,newS
      real*8 Qcond
      real*8 :: lw  = 2.8e+6
      real*8 :: cpp = 744.5

      integer ilay

      real*8   derf

      real*8 dqpi

!     added jrs

      real*8 :: tjs

c  Treatment
c  ---------

      Cste = dt * 4 * pi * dpden_ice

c     Start loop over heights
      DO 100 l = 1, nz

        ilay = 2 * l + 2

        p      = pl(ilay)
        temp   = tl(ilay)
        h2ovap = qpi(ilay,iMa_vap)

c       Save the values of the c arrays before condensation
        do i = 1, naer
          qpisav(i) = qpi(ilay,i)
        enddo

        call watsat(temp,p,qsat)

c       Get the partial presure of water vapor and its saturation ratio
        ph2o      = h2ovap * (44./18.) * p * 100.
        sat_ratio = h2ovap / qsat

c       Expand the dust moments into a binned distribution
        Mo = qpi(ilay,iMa_dt)
        No = qpi(ilay,iNb_dt)+ 1.d-50
        Rn = ( Mo / No * 0.75 / pi / aerdens(iMa_dt) )**(athird)
     *       * exp( -1.5 * stdv(iMa_dt)**2.)
        Rn = min( max(Rn,aerad(1)) , aerad(nbin) )
        Rm = Rn * exp( 3 * stdv(iMa_dt)**2. )
        Rn = 1. / Rn
        Rm = 1. / Rm
        dev2 = 1 / ( sqrt(2.) * stdv(iMa_dt) )
        do i = 1, nbin
          n_aer(i) = 0.5 * No * ( derf( dlog(rb(i+1)*Rn) * dev2 )
     &                           -derf( dlog(rb(i) * Rn) * dev2 ) )
          m_aer(i) = 0.5 * Mo * ( derf( dlog(rb(i+1)*Rm) * dev2 )
     &                           -derf( dlog(rb(i) * Rm) * dev2 ) )
        enddo

c       Get the rates of nucleation
        call nuclea(ph2o,temp,sat_ratio,n_aer,rate)

        dN = 0.
        dM = 0.
        do i = 1, nbin
          n_aer(i) = n_aer(i) / ( 1. + rate(i)*dt )
          m_aer(i) = m_aer(i) / ( 1. + rate(i)*dt )
          dN       = dN + n_aer(i) * rate(i) * dt
          dM       = dM + m_aer(i) * rate(i) * dt
        enddo

        Qcond = 0.

        IF (qpi(ilay,iNb_cld).ge.1.e-20.and.sat_ratio.ne.1) THEN

        Mo   = qpisav(iMa_cld) + qpisav(iMa_cor) + 1.d-50
        No   = qpisav(iNb_cld)
        dens = qpisav(iMa_cld) / Mo * aerdens(iMa_cld)
     *        +qpisav(iMa_cor) / Mo * aerdens(iMa_dt)
        dens = max(dens,dpden_ice)
        rad  = ( Mo / No * 0.75 / pi / dens ) **(athird)
        if (Mo.lt.1.e-20) rad = 1.e-8
        seq  = exp( 2.*sig(temp)*mh2o / (dpden_ice*rgp*temp*rad) )

        tjs = ph2o/sat_ratio
        call growthrate(dt,temp,p,ph2o,tjs,seq,rad,gr)
!       call growthrate(dt,temp,p,ph2o,ph2o/sat_ratio,seq,rad,gr)

        up  = Cste * gr * rad * No * seq + h2ovap
        dwn = Cste * gr * rad * No / qsat+ 1.

        Ctot = qpisav(iMa_cld) + h2ovap

        newvap = min(up/dwn,Ctot)

        gr = gr * ( newvap/qsat - seq )
        Qcond = max( Cste * No * rad * gr , -qpisav(iMa_cld) )
        if (latent_heat) then
          newT = temp + Qcond * lw / cpp
          call watsat(newT,p,qsat)
          newS = ( h2ovap - Qcond ) / qsat
          if (sat_ratio.lt.1.and.newS.gt.1.01) then
            call findQ(h2ovap,p,temp,Qcond)
          elseif(sat_ratio.gt.1.and.newS.lt.0.99) then
            call findQ(h2ovap,p,temp,Qcond)
          endif
        endif

        qpi(ilay,iMa_cld) = qpisav(iMa_cld) + Qcond

        if (qpi(ilay,iMa_cld).le.0) then
          qpi(ilay,iMa_cld) = 0.
          qpi(ilay,iMa_dt ) = qpisav(iMa_dt) + qpisav(iMa_cor)
          qpi(ilay,iNb_dt ) = qpisav(iNb_dt) + qpisav(iNb_cld)
          qpi(ilay,iMa_cor) = 0.
          qpi(ilay,iNb_cld) = 0.
        endif

        ENDIF

        qpi(ilay,iMa_dt ) = qpi(ilay,iMa_dt ) - dM
        qpi(ilay,iNb_dt ) = qpi(ilay,iNb_dt ) - dN
        qpi(ilay,iMa_cor) = qpi(ilay,iMa_cor) + dM
        qpi(ilay,iNb_cld) = qpi(ilay,iNb_cld) + dN

        qpi(ilay,iMa_vap) = h2ovap - (qpi(ilay,iMa_cld)-qpisav(iMa_cld))
        if (latent_heat) tl(ilay) = temp + Qcond * lw / cpp

100   CONTINUE

      return
      end

*******************************************************
*                                                     *
      subroutine nuclea(ph2o,temp,sat,n,nucrate)
*                                                     *
*   This subroutine computes the nucleation rate      *
*   as given in Pruppacher & Klett (1978) in the      *
*   case of water ice forming on a solid substrate.   *
*     Definition refined by Keese (jgr,1989)          *
*                                                     *
*******************************************************

      use grid_h
      use constants_h, only: PI
      use cldcommon_h

      implicit none

      integer i

      real*8 n(nbin)

      real*8 nucrate(nbin)
      real*8 ph2o,temp,sat

      real*8 nh2o
      real*8 sig          ! Water-ice/air surface tension  (N.m)
      real*8 rstar        ! Radius of the critical germ (m)
      real*8 gstar        ! # of molecules forming a critical embryo
      real*8 x            ! Ratio rstar/radius of the nucleating dust particle
      real*8 fistar       ! Activation energy required to form a critical embryo (J)
      real*8 zeldov       ! Zeldovitch factor (no dim)
      real*8 fshape       ! function defined at the end of the file
      real*8 sinteta      ! sine of the contact angle
      real*8 deltaf


      if (sat .gt. 1.) then    ! minimum condition to activate nucleation

        nh2o    = ph2o / kbz / temp
        rstar  = 2. * sig(temp) * vo1 / (rgp*temp*dlog(sat))
        gstar  = 4. * nav * pi * (rstar**3) / (3.*vo1)

c       Loop over size bins
        do 200 i=1,nbin

       if ( n(i) .eq. 0. ) then  ! no dust, no need to compute nucleation
          nucrate(i)=0.
          goto 200
        endif

        x      = aerad(i) / rstar
        fistar = (4./3.*pi) * sig(temp) * (rstar**2.) * fshape(mteta,x)
        deltaf = min( max((2.*desorp-surfdif-fistar)/(kbz*temp)
     &           , -100.), 100.)

        if (deltaf.eq.-100.) then
          nucrate(i) = 0.
        else
          zeldov = sqrt ( fistar / (3.*pi*kbz*temp*(gstar**2.)) )
          nucrate(i)= zeldov * kbz* temp * rstar
     &                * rstar * 4. * pi * ( nh2o*aerad(i) )**2.
     &                / ( fshape(mteta,x) * nus * m0 )
     &                * dexp (deltaf)
        endif

200   continue

      else

        do i=1,nbin
          nucrate(i) = 0.
        enddo

      endif

      return
      end

*********************************************************
      real*8 function sig(t)
*    this function computes the surface tension (N.m)   *
*   between water ice and air as a function of temp.    *
*********************************************************

      real*8 t

      sig = (141. - 0.15 * t) * 1.e-3

      return
      end

*********************************************************
      real*8 function fshape(cost,rap)
*        function computing the f(m,x) factor           *
* related to energy required to form a critical embryo  *
*********************************************************

      real*8 cost,rap
      real*8 phi

      phi = sqrt( 1. - 2.*cost*rap + rap**2 )
      a = 1. + ( (1.-cost*rap)/phi )**3
      b = (rap**3) * (2.-3.*(rap-cost)/phi+((rap-cost)/phi)**3)
      c = 3. * cost * (rap**2) * ((rap-cost)/phi-1.)

      fshape = 0.5*(a+b+c)

      if (rap.gt.3000.) fshape = ((2.+cost)*(1.-cost)**2)/4.

      return
      end

****************************************************************
      subroutine findQ(Qvap,p,temp,Qcond)
*     Start moddfied Newton-Raphson method to find             *
*     the amount of condensed/sublimed ice required to         *
*     reach saturation. Iterations are made to                 *
*     solve the F(Qc)=S-1=0 equation (Qc=condensed mass,       *
*     Saturation ratio) in the case of latent heat             *
*     release associated to Qc..                               *
****************************************************************
      implicit none

      real*8 Qvap,p,temp,Qcond

      real*8 x1,x2,xl,xh,f,df,fl,fh
      real*8 dx,dxold,tempo
      real*8 rtsafe

      integer i

      x1 = 0.
      x2 = Qcond
      call newton(x1,Qvap,p,temp,fl,df)
      call newton(x2,Qvap,p,temp,fh,df)

      if (fl*fh.ge.0.) then
        print*,'root not bracketed'
        stop
      endif
      if (fl.lt.0.) then
        xl = x1
        xh = x2
      else
        xh = x1
        xl = x2
      endif
      rtsafe = 0.5 * (x1+x2)
      dxold  = abs(x2-x1)
      dx     = dxold
      call newton(rtsafe,Qvap,p,temp,f,df)
      do i = 1, 500
        if ( ((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f).ge.0.
     *      .or. abs(2.*f).gt.abs(dxold*df) ) then
          dxold = dx
          dx    = 0.5 * (xh-xl)
          rtsafe= xl+dx
          if (xl.eq.rtsafe) goto 667
        else
          dxold = dx
          dx    = f / df
          tempo  = rtsafe
          rtsafe= rtsafe - dx
          if (tempo.eq.rtsafe) goto 667
        endif
        if (abs(dx).lt.1.e-8) goto 667
          call newton(rtsafe,Qvap,p,temp,f,df)
          if (f.lt.0) then
            xl = rtsafe
          else
            xh = rtsafe
          endif
        enddo
        print*,'500 reached'
667     continue

        Qcond = rtsafe

        return
        end

*********************************************************
      subroutine newton(dq,q,press,temp,f,df)
*     subroutine called during the Newton-Raphson loop  *
*     to get the function and its derivative for the    *
*              condensed mass determination             *
*********************************************************

      implicit none

      real*8 dq,q,press,temp
      real*8 newT
      real*8 qsat

      real*8 :: lw  = 2.8e+6
      real*8 :: cpp = 744.5
      save cpp,lw

      real*8 f,df

      newT = temp + dq * lw / cpp
      newT = max(newT,60.)

      call watsat(newT,press,qsat)
      qsat = max(qsat,1.d-50)

      f    = (q-dq) / qsat - 1.
      df   = - 2. / qsat - (f+1.) * (6146.1*lw/cpp) / newT**2.

      return
      end
