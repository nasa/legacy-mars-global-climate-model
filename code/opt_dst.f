      subroutine opt_dst(jcmn,icmn,qtrace,pl,
     *                   Qxv,Qxi,Qsv,Qsi,gv,gi,
     *                   Qextref,TAUREF)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

      use grid_h
      use constants_h, only: PI, GRAV
      use radinc_h
      use dtcommon_h
      use cldcommon_h
      use standard_h, only: tauref3d
 
      implicit none

c  Arguments
c  ---------

      integer jcmn,icmn

      real*8  qtrace(l_jsize,l_isize,l_layers,ntrace)
      real*8  pl(l_levels)

      real*8  Qxv(L_LEVELS+1,L_NSPECTV)
      real*8  Qsv(L_LEVELS+1,L_NSPECTV)
      real*8  gv(L_LEVELS+1,L_NSPECTV)

      real*8  Qxi(L_LEVELS+1,L_NSPECTI)
      real*8  Qsi(L_LEVELS+1,L_NSPECTI)
      real*8  gi(L_LEVELS+1,L_NSPECTI)

      real*8  Qextref(L_LEVELS+1)
      real*8  TAUREF(L_LEVELS+1)

c  Local variables
c  ---------------

      integer i,iwav,j,k,l
      integer NN

      logical do_it

      real*8 dens,cst,dev,dev2
      real*8 Rn,Rs,Ao,Mo,No
      real*8 surf(nbin_rt)
      real*8 MASS

      real*8 derf

c Initialyze various variables
c ----------------------------

      DO K = 1, L_LEVELS+1
        Qextref(K) = 1.
        TAUREF(K)  = 0.
      ENDDO

      do i = 1, nlonv
        DO K = 1, L_LEVELS+1
          Qxv(K,i) = Qextref(K)
          Qsv(K,i) = Qextref(K) * 0.99
          gv(K,i)  = 0.
        ENDDO
      enddo
      do i = 1, nloni
        DO K = 1, L_LEVELS+1
          Qxi(K,i) = Qextref(K)
          Qsi(K,i) = Qextref(K) * 0.99
          gi(K,i)  = 0.
        ENDDO
      enddo

      dev  = dev_dt
      dens = dpden_dt 
      cst  = 0.75 / (pi*dens)
      dev2 = 1. / ( sqrt(2.)*dev )

c  Treatment 
c  ---------

      taudust(jcmn,icmn,1) = 0.
      taudust(jcmn,icmn,2) = 0.

      DO L = 1, L_LAYERS

      do_it = qtrace(jcmn,icmn,l,iMa_dt) .gt. 1.e-8
     *                .and.
     *        qtrace(jcmn,icmn,l,iNb_dt) .gt. 1.

      if (do_it) then

c     Get the cross-section mean radius (Rs) of the log-normal distribution
      Mo = qtrace(jcmn,icmn,l,iMa_dt)          ! Mass mixing ratio
      No = qtrace(jcmn,icmn,l,iNb_dt) + 1.     ! Number mixing ratio

      Rs = ( Mo/No*cst )**(athird) * dexp( -0.5*dev**2. )

c     Get the total cross sectional area Ao of water ice particles
      Ao = No * pi * Rs**2.

c *********************************************************************

c     Define the cross-section weighted distribution, i.e. surface/size 
c	 bin.  Change Rs to Reff.  MAK 6 May 2008.

      Rs = Rs * dexp ( 1.5 * dev**2. )

c *********************************************************************
 
      Rs = min( max(Rs,1.e-7) , 50.e-6 )
      Rs = 1. / Rs

      do i = 1, nbin_rt
        surf(i) = 0.5 * ( derf( dlog(radb_rt(i+1)*Rs) * dev2 )
     *                   -derf( dlog(radb_rt(i)  *Rs) * dev2 ) )
      enddo

c     Get the average values of <Qext>, <Qscat>, and <g> for the whole distribution.
      DO NN = 1, 2
        K = 2*L + 1 + NN

        do iwav = 1, nlonv
          Qxv(K,iwav) = 0.
          Qsv(K,iwav) = 0.
          do i = 1, nbin_rt
            Qxv(K,iwav) = Qxv(K,iwav)+ surf(i) * qextv_dst(i,iwav)
            Qsv(K,iwav) = Qsv(K,iwav)+ surf(i) * qscatv_dst(i,iwav)
            gv(K,iwav)  = gv(K,iwav) + surf(i) * gv_dst(i,iwav)
          enddo
          Qsv(K,iwav) = min( Qsv(K,iwav) , 0.99999*Qxv(K,iwav) )
        enddo

        do iwav = 1, nloni
          Qxi(K,iwav) = 0.
          Qsi(K,iwav) = 0.
          do i = 1, nbin_rt
            Qxi(K,iwav) = Qxi(K,iwav)+ surf(i) * qexti_dst(i,iwav)
            Qsi(K,iwav) = Qsi(K,iwav)+ surf(i) * qscati_dst(i,iwav)
            gi(K,iwav)  = gi(K,iwav) + surf(i) * gi_dst(i,iwav)
          enddo
          Qsi(K,iwav) = min( Qsi(K,iwav) , 0.99999*Qxi(K,iwav) )
        enddo

        Qextref(K) = Qxv(K,L_NREFV)
        MASS       = 100. * (pl(K) - pl(K-1)) / grav
        TAUREF(K)  = Ao * Qextref(K) * MASS
        TAUREF3D(JCMN,ICMN,K)=TAUREF(K)

c       For diagnostics: cloud opacity at ref wavelengths (in the vis and ir)
        taudust(jcmn,icmn,1) = taudust(jcmn,icmn,1) + TAUREF(K)
        taudust(jcmn,icmn,2) = taudust(jcmn,icmn,2)
     &                        +TAUREF(K) / Qextref(K) * 
     &                         ( Qxi(K,4)-Qsi(K,4) )

      ENDDO

      endif    ! End condition on do_it

      ENDDO    ! End loop over layers

      return

      end
