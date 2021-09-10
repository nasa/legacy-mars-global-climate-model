      subroutine newtg(als,s6star,downir,rhouch,rhoucht,scond,stemp,
     *                 sthick,tg,jin,iin,psfg,qmh2o,qice,dt,
     *                 subflx,polarcap)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

C  Use modified Newton-Raphson method for finding the ground 
C  temperature.

      use constants_h, only: cp

      implicit none
  
      integer J, N, jin, iin
      real*8  df, dx, dxold, f, fh, fl, temp, TH, TL
      real*8  T1, T2, xacc
      real*8  astar, s6star, scond, stemp, sthick, als, tg
      real*8  downir, rhouch, rhoucht

!  10-07-11 Latent heating variables

!  psfg     - PSAT
!  qmh2o    - QTRACE(J,I,NLAY,iMa_vap)
!  qice     - GNDICE(J,I)
!  dt       - dynamical time step
!  polarcap - true if grid point is part of the north polar water ice
!             cap

      real*8  :: psfg, qmh2o, qice, dt
      real*8  :: qgnd, wflux
      real*8  :: subflx
      logical :: polarcap

C  MAXIT is the maximum number of iterations allowed for convergence.

      integer MAXIT
      parameter(MAXIT = 30)

C  T1 is the initial low temperature
C  T2 is the initial high temperature
C  XACC is the accuracy (in Kelvins)

      parameter(T1   = 50.0D0)
      parameter(T2   = 350.0D0)
      parameter(XACC = 1.0D-3)

C======================================================================C

      astar = (1.0D0-ALS)*S6STAR
 
      call funcd(astar,downir,rhouch,rhoucht,scond,stemp,sthick,T1,FL,
     *           df,psfg,qmh2o,qice,polarcap)
      call funcd(astar,downir,rhouch,rhoucht,scond,stemp,sthick,T2,FH,
     *           df,psfg,qmh2o,qice,polarcap)

      if(FL.eq.0.0) then
        TG = T1
        goto 100
      elseif(FH.eq.0.0) then
        TG = T2
        goto 100
      elseif(FL.LT.0.0) then
        TL = T1
        TH = T2
      else
        TL = T2
        TH = T1
      end if

      TG    = 0.5*(T1+T2)
      dxold = abs(T2-T1)
      dx    = dxold

      call funcd(astar,downir,rhouch,rhoucht,scond,stemp,sthick,TG,f,
     *           df,psfg,qmh2o,qice,polarcap)

      do J=1,MAXIT
        if(((TG-TH)*df-f)*((TG-TL)*df-f) .ge. 0.0   .or.
     *       abs(2.0*f) .gt. abs(dxold*df)               ) then

          dxold = dx
          dx    = 0.5*(TH-TL)
          TG    = TL+dx
          if(TL.eq.TG) then
            goto 100
          end if
        else

          dxold = dx
          dx    = f/df
          temp  = TG
          TG    = TG - dx
          if(TEMP.eq.TG) then
            goto 100
          end if
        end if

        if(abs(dx).lt.xacc) then
          goto 100
        end if

        call funcd(astar,downir,rhouch,rhoucht,scond,stemp,sthick,TG,f,
     *             df,psfg,qmh2o,qice,polarcap)

        if(F.lt.0.0) then
          TL = TG
        else
          TH = TG
        end if

      end do

C  If we reach this statement, we've done MAXIT number of iterations 
C  without convergence.

      write(6,'("------- SUBROUTINE NEWTG")')
      write(6,'("Maximum number of iterations ",i3," exceeded.")') MAXIT
      write(6,'("J = ",i3,5x,"I = ",I3)') JIN, Iin
      write(6,'("ALS     = ",1pe12.5)') ALS
      write(6,'("s6star  = ",1pe12.5)') s6star
      write(6,'("DOWNIR  = ",1pe12.5)') DOWNIR
      write(6,'("RHOUCH  = ",1pe12.5)') RHOUCH
      write(6,'("RHOUCHT = ",1pe12.5)') RHOUCHT
      write(6,'("SCOND   = ",1pe12.5)') SCOND
      write(6,'("STEMP   = ",1pe12.5)') STEMP
      write(6,'("STHICK  = ",1pe12.5)') STHICK
      write(6,'("TG      = ",1pe12.5)') TG
      stop

  100 continue

      wflux = 0.0
 
      if(qice.gt.0.0 .and. .not.polarcap) then
        qgnd  = (18.0/44.0)*6.11*exp(22.5*(1.0-(273.16/tg)))/psfg
        wflux = -rhouch*(qmh2o-qgnd)/cp
        if(wflux*dt.ge.qice) then
          wflux = qice/dt
          qice  = 0.0
        else
          qice  = qice - wflux*dt
        end if
      elseif(polarcap) then
        qgnd  = (18.0/44.0)*6.11*exp(22.5*(1.0-(273.16/tg)))/psfg
        wflux = -rhouch*(qmh2o-qgnd)/cp
        qice  = qice - wflux*dt
      end if

!  Note:  wflux > 0 is sublimation
!         wflux < 0 is condensation

      subflx = subflx + wflux*dt

      return
      end
