      subroutine opt_cld(jcmn,icmn,qtrace,pl,
     *                   Qxv,Qxi,Qsv,Qsi,gv,gi,
     *                   Qextrefcld,TAUREFCLD)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

      use grid_h
      use constants_h, only: PI, GRAV
      use radinc_h
      use dtcommon_h
      use cldcommon_h

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

      real*8  Qextrefcld(L_LEVELS+1)
      real*8  TAUREFCLD(L_LEVELS+1)

c  Local variables
c  ---------------

      integer i,iwav,j,k,l,irap
      integer NN

      real*8 dens,cst,dev,dev2
      real*8 Rn,Rs,Ao,Mo,No
      real*8 surf(nbin_rt)
      real*8 MASS

      real*4 mantletocore
      logical do_it

      real*8 derf

c Initialyze various variables
c ----------------------------

      DO K = 1, L_LEVELS+1
        Qextrefcld(K) = 1.
        TAUREFCLD(K) = 0.
      ENDDO

      do i = 1, nlonv
        DO K = 1, L_LEVELS+1
          Qxv(K,i) = Qextrefcld(K)
          Qsv(K,i) = Qextrefcld(K) * 0.99
          gv(K,i)  = 0.
        ENDDO
      enddo

      do i = 1, nloni
        DO K = 1, L_LEVELS+1
          Qxi(K,i) = Qextrefcld(K)
          Qsi(K,i) = Qextrefcld(K) * 0.99
          gi(K,i)  = 0.
        ENDDO
      enddo

c  Treatment 
c  ---------

      dev= dev_ice
      dev2 = 1. / ( sqrt(2.)*dev )

      taucloud(jcmn,icmn,1) = 0.
      taucloud(jcmn,icmn,2) = 0.

      DO L = 1, L_LAYERS

c     Do not do the computations if the amount of cloud is too low
      do_it = qtrace(jcmn,icmn,l,iMa_cld)+qtrace(jcmn,icmn,l,iMa_cor) 
     *       .gt. 1.e-7
     *       .and.
     *        qtrace(jcmn,icmn,l,iNb_cld) .gt. 1.

      if (do_it) then

c     Determine the ratio of the dust core radius over that of the ice mantle
      mantletocore =   qtrace(jcmn,icmn,l,iMa_cor)/dpden_dt 
     *              /( qtrace(jcmn,icmn,l,iMa_cld)/dpden_ice
     *                +qtrace(jcmn,icmn,l,iMa_cor)/dpden_dt )
      mantletocore = mantletocore**(athird)

c     Find the index to which corresponds the optical properties of the 
c     core to mantle radius ratio. Those properties were determined off-line
c     using the Toon and Ackerman coated spheres code.
      irap = nratio
      do i = 1, nratio
        if (mantletocore.lt.cor_ratio(i) .and. i.ne.1) then 
          irap = i - 1
          exit
        elseif (mantletocore.eq.cor_ratio(i)) then 
          irap = i
          exit
        elseif (mantletocore.lt.cor_ratio(i) .and. i.eq.1) then
          irap = 1
          exit
        endif
      enddo

c     Get the cross-section mean radius (Rs) of the log-normal distribution
      Mo = qtrace(jcmn,icmn,l,iMa_cld) 
     *    +qtrace(jcmn,icmn,l,iMa_cor)       ! Mass mixing ratio
      No = qtrace(jcmn,icmn,l,iNb_cld)       ! Number mixing ratio

      dens =  qtrace(jcmn,icmn,l,iMa_cld) / Mo * dpden_ice
     *       +qtrace(jcmn,icmn,l,iMa_cor) / Mo * dpden_dt

      cst  = 0.75 / (pi*dens)
      Rs = ( Mo/No*cst )**(athird) * dexp( -0.5*dev**2. )

c     Get the total cross sectional area Ao of the water ice particles
      Ao = No * pi * Rs**2.

c *********************************************************************

c     Define the cross-section weighted distribution, i.e. surface/size 
c	bin.   Change Rs to Reff.  MAK 6 May 2008.

      Rs = Rs * dexp ( 1.5 * dev**2. )

c *********************************************************************

      Rs = min( max(Rs,1.e-7) , 100.e-6 )
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
            Qxv(K,iwav) = Qxv(K,iwav)+ surf(i) * qextv_cld(irap,i,iwav)
            Qsv(K,iwav) = Qsv(K,iwav)+ surf(i) * qscatv_cld(irap,i,iwav)
            gv(K,iwav)  = gv(K,iwav) + surf(i) * gv_cld(irap,i,iwav)
          enddo
          Qsv(K,iwav) = min( Qsv(K,iwav) , 0.99999*Qxv(K,iwav) )
        enddo

        do iwav = 1, nloni
          Qxi(K,iwav) = 0.
          Qsi(K,iwav) = 0.
          do i = 1, nbin_rt
            Qxi(K,iwav) = Qxi(K,iwav)+ surf(i) * qexti_cld(irap,i,iwav)
            Qsi(K,iwav) = Qsi(K,iwav)+ surf(i) * qscati_cld(irap,i,iwav)
            gi(K,iwav)  = gi(K,iwav) + surf(i) * gi_cld(irap,i,iwav)
          enddo
          Qsi(K,iwav) = min( Qsi(K,iwav) , 0.99999*Qxi(K,iwav) )
        enddo

        Qextrefcld(K) = Qxv(K,L_NREFV)
        MASS          = 100. * (pl(K) - pl(K-1)) / grav
        TAUREFCLD(K)  = Ao * Qextrefcld(K) * MASS

c       For diagnostics: cloud opacity at ref wavelengths (in the vis and ir)
        taucloud(jcmn,icmn,1) = taucloud(jcmn,icmn,1) + TAUREFCLD(K)
        taucloud(jcmn,icmn,2) = taucloud(jcmn,icmn,2) 
     &                         +TAUREFCLD(K) / Qextrefcld(K) * 
     &                          ( Qxi(K,4) - Qsi(K,4) )
      ENDDO

      endif    ! Condition on cloud amount 

      ENDDO    ! End loop over layers

      return

      end
