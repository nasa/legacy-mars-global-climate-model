      subroutine matrix(i,dt,dtmu,rmu1mu,ustar,vk,rlnzz,tg,
     *                  r1roro,r2roro,r1roxi,r2roxi,rnum,dxi,ro,vect,
     *                  pc,qc,rc,psi,mat,cdh,cdm,qsat,coef)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

      use grid_h
      use pbl_defines_h

      implicit none

c     local areas: none
c     local variables:
      real*8 qsat   ! mass mixing ratio of water vapor near the surface at temp=tg
      real*8 coef

c     global arrays:

      real*8 mat(n,3)
      real*8 r1roro(n),r2roro(n),r1roxi(n),r2roxi(n),vect(n),
     *       rnum(2*n+1),dxi(2*n+1),ro(2*n+1),
     *       pc(2*n+1,nvar),qc(2*n+1,nvar),rc(2*n+1,nvar),
     *       psi(2*n+1,nvar)

      integer :: i, j
      real*8  :: dtmu, ustar, cdm, cdh, tg, vk, rlnzz, rmu1mu, dt

C======================================================================C

      do j = 1, n-1
         mat(j+1, 1) = r1roxi(j)  *  pc(2*j+1,i)      
         mat(j, 3)   = r2roxi(j)  *  pc(2*j+1,i)
      end do

      mat(1, 1) = 0.0
      mat(n, 3) = 0.0

      do j = 2, n-1
         mat(j, 2) = 1. - r1roro(j) * mat(j,1) -
     *                    r2roro(j) * mat(j,3)      -
     *                    dtmu * qc(2*j,i) 
      end do

      mat(1, 2) = 1. -  r2roro(1) * mat(1,3) - dtmu  *  qc(2,i) 

      do j = 2, n-1
         vect(j) = -rmu1mu * mat(j,1) * psi(2*j-2,i) +
     *              (1. + rmu1mu * (1. - mat(j,2)))   *
     *              psi(2*j,i) - rmu1mu * mat(j,3)    *
     *              psi(2*j+2,i) + rc(2*j,i) * dt 
      end do

      vect(1) = (1. + rmu1mu * (1. - mat(1,2)))   *
     *           psi(2,i) - rmu1mu   * mat(1,3)    *
     *           psi(4,i) + rc(2,i)   * dt        


      if(i.eq.1 .or. i.eq.2) then

c     Surface stress boundary condition in u and v

         mat(n,2) = 1. - r1roro(n) * mat(n,1) - dtmu * qc(2*n,i) -
     *                   ( ustar / rnum(2*n) ) *
     *                   ( sqrt(cdm) ) * ( dtmu / dxi(2*n) )

         vect(n)  = - rmu1mu * mat(n,1) * psi(2*n-2,i) +
     *                (1. + rmu1mu * (1. - mat(n,2))) * psi(2*n,i) +
     *                rc(2*n,i) * dt

c     Surface heat flux boundary condition in theta

      elseif (i .eq. 3) then

         mat(n,2) = 1. - r1roro(n) * mat(n,1) - dtmu * qc(2*n,i) -
     *              ( ustar / rnum(2*n) ) * ( cdh ) *
     *              ( dtmu / dxi(2*n) )

         vect(n)  = - rmu1mu * mat(n,1) * psi(2*n-2,i) +
     *               (1. + rmu1mu * (1. - mat(n,2))) * psi(2*n,i) +
     *               rc(2*n,i) * dt - ( ustar / rnum(2*n) ) *
     *               ( cdh ) * 
     *               ( dt * ro(2*n)* tg/dxi(2*n)) 

c     Water vapor flux condition (sublimation from surface)

      elseif (i .eq. 4) then

         mat(n,2) = 1. - r1roro(n) * mat(n,1) - dtmu * qc(2*n,i) -
     *              ( ustar / rnum(2*n) ) * ( cdh ) * coef *
     *              ( dtmu / dxi(2*n) )

         vect(n)  = - rmu1mu * mat(n,1) * psi(2*n-2,i) +
     *               (1. + rmu1mu * (1. - mat(n,2))) * psi(2*n,i) +
     *               rc(2*n,i) * dt - ( ustar / rnum(2*n) ) *
     *               ( cdh ) * coef *
     *               ( dt * ro(2*n)* qsat/dxi(2*n)) 
c
c      else
c
c     Surface conditions for tracer quantities 
c
c         mat(n,2) = 1. - r1roro(n) * mat(n,1) - dtmu * qc(2*n,i)
c
c         vect(n)  = - rmu1mu * mat(n,1) * psi(2*n-2,i) +
c     *                (1. + rmu1mu * (1. - mat(n,2))) * psi(2*n,i) +
c     *                rc(2*n,i) * dt
c
      end if

      return
      end
