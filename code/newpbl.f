      subroutine newpbl(z0, epsl0, vk, alphl0, ric, dtmu, 
     *                  rmu1mu,dxi,du,I1D,J1D,nstep,ptrop,tg,
     *                  pi,dt,sigma,pl,om,upi,vpi,teta,
     *                  htflux,strx,stry,rhouch,rkh,qrad,qpi,qpig,
     *                  latent,polarcap,sup,sdn,h2oflux)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

      use grid_h
      use constants_h, only: GRAV, RGAS, CP
      use pbl_defines_h
      use cldcommon_h

      implicit none

c     First dimension variables passed from GCM

      real*8 om(2*n+3), rhouch
      real*8 upi(2*n+3),vpi(2*n+3),teta(2*n+3),
     *          sigma(2*n+3),pl(2*n+3),z(2*n+1),qrad(2*n+3)

C  May 2002  -  dust tracer variables (QPI)

      real*8 qpi(2*n+3,ntrace)
      real*8 qpig(ntrace)

c     Now dimension variables local to newpbl

      real*8 r1roxi(n),r2roxi(n),r1roro(n),r2roro(n),vect(n)
      real*8 mat(n,nvar)
      real*8 dxi(2*n+1),q0(2*n+1),du(l_jsize,l_isize,2*n+1),
     *          ro(2*n+1),rnum(2*n+1),
     *          vg(2*n+1),xinum(2*n+1),rnumnum(2*n+1),ronum(2*n+1),
     *          rkm(2*n+1),rkh(2*n+1),
     *          dudz
      real*8 psi(2*n+1,nvar),
     *          pc(2*n+1,nvar),qc(2*n+1,nvar),rc(2*n+1,nvar)


c     Local variables
      real*8 qsat    ! mass mixing ratio of water vapor near the surface at temp=tg
      real*8 h2oflux, qtmp, qold
      real*8 checkm(n,nvar),checkv(n)

      real*4 latent
      real*8   coef

      logical polarcap

!     implicit none

      integer :: i, j, k, l, i1d, j1d, nstep
      real*8  :: cdh, w, dphi, thstar, ustar, tg, dt, ptrop, rmu1mu
      real*8  :: uold, htflux, stry, strx, cdm, alphl0, epsl0, pi
      real*8  :: dtmu, ric, vk, alpha, tbar, hscale, rlnzz, z0
      real*8  :: sup, sdn

C======================================================================C
     
C     Constants for boundary layer routines are initialized in initpbl.f
C     These variables are:

C     z0, epsl0, vk, alphl0, ric, rf1, rf2, rf3, rf4, sh1, sh2, gh1, 
C     gh2, sm1, sm2, sm3, sm4, sm5, gm1, gm2, gm3, gm4, gm5, dtmu, 
C     rmu1mu, dxi 

c     These steps must be done each time through:

      do k = 2,2*n,2
        psi(k,1) = upi(k+2)
        psi(k,2) = vpi(k+2)
        psi(k,3) = teta(k+2)
        psi(k,4) = qpi(k+2,iMa_vap)
      end do

c     Calculate heights

      z(2*n+1) = 0.0

      do k = 2*n, 2, -2
        tbar = 0.0

        do l = k,2*n,2
          tbar = tbar + teta(l+2)*om(l+2)* (sigma(l+3) - sigma(l+1))
        end do

        tbar   = tbar / (sigma(2*n+3) - sigma(k+1))
        hscale = rgas*tbar/grav
        z(k)   = hscale * log( pl(2*n+3) / pl(k+2) )
        z(k-1) = hscale * log( pl(2*n+3) / pl(k+1) )
      end do

      rlnzz = log(z(2*n) / z0)

      do k = 1,2*n+1
        ro(k)   = + 100.0 *pl(k+2) / rgas / (om(k+2)*teta(k+2))
        rnum(k) = - 100.0 *pi / ( ro(k)*grav )
      end do

      do k = 2, 2*n
        ronum(k) = ro(k) * rnum(k)
      end do

      do j = 1, n-1
        r1roxi(j)   = -dtmu * ronum(2*j+1) / ronum(2*j) /
     *                    dxi(2*j+2)   / dxi(2*j+1)
        r2roxi(j)   = -dtmu * ronum(2*j+1) / ronum(2*j+2) /
     *                    dxi(2*j)     / dxi(2*j+1)
        r1roro(j+1) = ronum(2*j)   / ronum(2*j+2)
        r2roro(j)   = ronum(2*j+2) / ronum(2*j)
      end do

      do k = 3, 2*n-1, 2
        xinum(k)   = 1. / dxi(k) / rnum(k)
        rnumnum(k) = rnum(k) * rnum(k)
      end do

      call eddycoef(nstep,ric,vk,alphl0,epsl0,dxi,q0,du,I1D,J1D,
     *              xinum,z,rkm,rkh,psi,dt)

      call bndcond(nstep,tg,vk,rlnzz,z0,ustar,
     *             thstar,ric,dphi,w,z,psi,cdh,cdm)

      alpha = atan2( psi(2*n,2), psi(2*n,1) )
      strx = + ro(2*n+1)*ustar*ustar*cos(alpha)
      stry = + ro(2*n+1)*ustar*ustar*sin(alpha)

      htflux = - ro(2*n+1)*cp*ustar*thstar

C  Compute RHOUCH for the new method of determing Tg in TEMPGR:

      RHOUCH = RO(2*N+1)*CP*CDH*USTAR

      uold = psi(2*n,1)
      qold = qpi(2*n+2,iMa_vap)

      call scldef(ro,rnum,rkm,rkh,rnumnum,pc,qc,rc,psi,qrad)
      call watsat(tg,pl(2*n+3),qsat)

c     Reduce the sublimation flux by a coefficient coef.
c     Tunable parameter to avoid the formation of low-lying
c     clouds in summer above the north permanent cap.
      if (qsat.gt.qpi(2*n+2,iMa_vap)) then
        coef = 1.0D0
       else
        coef = 1.0D0
      endif

c     loop over variables

      do i = 1, nvar
        call matrix(i,dt,dtmu,rmu1mu,ustar,vk,rlnzz,tg,
     *              r1roro,r2roro,r1roxi,r2roxi,rnum,dxi,ro,vect,
     *              pc,qc,rc,psi,mat,cdh,cdm,qsat,coef)
c       Check the amount of subliming water ice.  - Feb 2003
c       If it exceeds the available surface water, change it.
        if (i .eq. 4) then
c       Make a preliminary computation of vect(n) (the true vect(n) is computed in solve.f) 
c       After scaling, vect(n) becomes qpi of water in the first layer.
          checkm(1,2) = mat(1,2) 
          checkv(1)   = vect(1)
          do j = 1, n-1
            checkm(j+1,2) = mat(j+1,2) - mat(j+1,1)*mat(j,3)/checkm(j,2)
            checkv(j+1)   = vect(j+1) - mat(j+1,1)*checkv(j)/checkm(j,2)
          end do
          qtmp = checkv(n) / checkm(n,2) / ro(2*n) 
c         Use a temporary value of qpi (qtmp) to compute the near surface water flux*dt    
c         (the flux accounts through dtmu for the type of scheme used; i.e. expl. or impl.) 
!         h2oflux = coef * dt * ro(2*n) / rnum(2*n) * cdh * ustar
!    *              * (qsat - dtmu/dt*qtmp - (1.0D0-dtmu/dt)*qold)
c         If subliming more than available (qpig), change vect(n).
c         This time, no boundary condition on the ground but water vapor flux is passed like 
c         a mass mixing ratio tendency in the first layer. 
            
          if((.not. polarcap) .and. (h2oflux.gt.qpig(iMa_vap))) then
             h2oflux  = qpig(iMa_vap)
          end if

          mat(n,2) = 1. - r1roro(n) * mat(n,1) - dtmu * qc(2*n,i)
          vect(n)  = - rmu1mu * mat(n,1) * psi(2*n-2,i) +
     *                 (1. + rmu1mu * (1. - mat(n,2))) * psi(2*n,i) +
     *                 rc(2*n,i) * dt + ro(2*n) * h2oflux /
     *                 ( (pl(2*n+3)-pl(2*n+1))*100./grav)
        endif

        call solve(i,vect,mat,psi)
      end do

      call descale (ro,rnum,rkm,rkh,rnumnum,psi)

      do k = 2,2*n,2
        upi(k+2)  = psi(k,1)
        vpi(k+2)  = psi(k,2)
        teta(k+2) = psi(k,3)
        qpi(k+2,iMa_vap) = psi(k,4)
      end do

c  updating water ice budget on the surface - Feb 2003
  
      qpig(iMa_vap) = qpig(iMa_vap) - h2oflux 
      if(.not.polarcap .and. qpig(iMa_vap).lt.0.0) then
         qpig(iMa_vap) = 0.0D0
      end if

      if (h2oflux .gt. 0.) then
        sup = h2oflux/dt
        sdn = 0.
      else
        sdn = (-1.) * h2oflux/dt
        sup = 0.
      endif

      if(latent_heat)  latent = 2.8e+6*(-h2oflux)/dt

      return
      end

      subroutine watsat(temp,press,qsat)
      
      real*8 pvs,temp,press,qsat

c Bob's vapor pressure
      pvs  = 6.11*exp(22.5*(1.0-(273.16/temp)))

      qsat = pvs * 18.0 / (44.0*press)
      
      return
      end
