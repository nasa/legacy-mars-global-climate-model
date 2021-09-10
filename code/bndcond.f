      subroutine bndcond(nstep,tg,vk,rlnzz,z0,ustar,
     *                   thstar,ric,dphi,w,z,psi,cdh,cdm)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

! New routine calculates ustar and thstar from rib, 
! based on Savijarvi, Icarus Vol 117 p121, section 2.1 (rib < 0)
! and Hourdin et al, JGR Vol 100 p5505, equations 5 and 6 (rib > 0)

      use grid_h
      use constants_h, only: grav
      use pbl_defines_h 

      implicit none

!     global variables:

      real*8 z(2*n+1), psi(2*n+1,nvar)

      integer :: nstep
      real*8  :: tg, vk, rlnzz, z0, ustar, thstar, ric, dphi, w, cdh 
      real*8  :: cdm, rib, frih, frim

!======================================================================C

      w  = sqrt( psi(2*n,1)*psi(2*n,1) + psi(2*n,2)*psi(2*n,2) )

!     For now use psi at midpoint of bottom layer.
!     Later interpolate to surface

      dphi = (psi(2*n,3) - tg)
      rib  = ((grav * z(2*n)) / ( psi(2*n,3)*w*w + 1.0e-09))*dphi

      if ( rib .ge. 0.0 ) then

        frih   = 1.0/(1.0+(15.0*rib/sqrt(1.0+5.0*rib)))
        frim   = 1.0/(1.0+(10.0*rib/sqrt(1.0+5.0*rib)))

      else if ( rib .lt. 0.0 ) then

! LMD formulation
 
!           frih   = 1.0 - 15.0*rib/(1+75.0*((vk/rlnzz)**2)*
!     &              sqrt(-rib*exp(rlnzz)))
!           frim   = 1.0 - 15.0*rib/(1+75.0*((vk/rlnzz)**2)*
!     &              sqrt(-rib*(1.0+exp(rlnzz))))

! Original, incorrect version

!       frih   = sqrt(1.0-16.0*rib)
!       frim   = sqrt(1.0-64.0*rib)

! Corrected 5-30-06

        frih = sqrt(1.0-64.0*rib)
        frim = sqrt(1.0-16.0*rib)

      end if

      cdh    = sqrt(frih)*vk/rlnzz
      cdm    = frim*((vk/rlnzz)**2)

      ustar  = sqrt(cdm)*w 
      thstar = cdh*dphi

      return
      end
