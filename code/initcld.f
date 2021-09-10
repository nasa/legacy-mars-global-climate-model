      subroutine initcld

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

C  Water ice tracer
C  GCM2.0  Feb 2003
C
C  Initialize the variables used in the cloud scheme

      use grid_h
      use constants_h, only: PI
      use dtcommon_h
      use cldcommon_h

      implicit none

      integer i,j,l

      real*8 :: rmin  = 0.1e-6
      real*8 :: rmax  = 10.e-6
      real*8 :: rbmin = 0.0001e-6
      real*8 :: rbmax = 1.e-2
 
      real*8 vrat_rt

      real*8 factor

      real*8 :: vistoir = 2.75

      real*8 a0

C======================================================================C

C  Set the various tracer index

c  For dust: Mass and Number 
      iMa_dt  = 1
      iNb_dt  = 2 
c  For water ice: Mass and Number 
      iMa_cld = 3 
      iNb_cld = 4 
c  For dust core: only Mass
      iMa_cor = 5
c  For water vapor: only Mass
      iMa_vap = 6 

C  Some general constants - fill variables that will be passed to the 
C  cloud scheme via cldcommon.h.

c  1/3 is computed here once for all

      athird  = (1. / 3.)

C  Water ice particle density

      dpden_ice = 917.

c  Avogadro number
      nav    = 6.02e23
c  Perfect gas constant
      rgp    = 8.3143
c  Boltzman constant
      kbz    = rgp / nav
c  Molecular weight of water  (kg)
      mh2o   = 18.01e-3
c  Weight of a water molecule (kg)
      m0     = mh2o / nav
c  Volume of a water molecule (m3)
      vo1    = mh2o / dpden_ice      
c  Activation energy for desorption of water on a dust-like substrate (J/molecule)
      desorp = 0.288e-19
c  Estimated activation energy for surface diffusion of water molecules (J/molecule)
      surfdif= desorp / 10.
c  Jump frequency of a water molecule (s-1)
      nus    = 1.e+13
c  
      a0     = 4.52E-10
c
      d_nuc  = a0 / sqrt(2.)

c  Contact parameter of water ice on dust ( m=cos(theta) )
!     mteta  = 0.95        !Franck's number
      mteta  = 0.975       !Franck recommends - newer value

c  nu_eff = exp( dev^2) - 1
c  Standard deviation of the dust distribution
      dev_dt = 0.63676    ! gives an effective variance of 0.5  for dust
c  Standard deviation of the water ice distribution
      dev_ice= 0.3087     ! gives an effective variance of 0.1  for water ice
!      dev_ice= 0.05     ! gives an effective variance of 0.02  for water ice

      aerdens(iMa_dt) = dpden_dt
      aerdens(iNb_dt) = dpden_dt
      aerdens(iMa_cld)= dpden_ice
      aerdens(iNb_cld)= dpden_ice
      aerdens(iMa_cor)= dpden_ice

      stdv(iMa_dt) = dev_dt
      stdv(iNb_dt) = dev_dt
      stdv(iMa_cld)= dev_ice
      stdv(iNb_cld)= dev_ice
      stdv(iMa_cor)= dev_ice

c     Definition of the size grid
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     aerad is the primary radius grid used for microphysics computation.
c     the grid spacing is based on a volume ratio between two
c     consecutive bins; i.e. vrat.
c     rb defines the boundary values for each aerad bin.

c     Volume ratio between two adjacent bins
      vrat = log(rmax/rmin) / float(nbin-1) *3.
      vrat = exp(vrat)

      rb(1)      = rbmin
      aerad(1)   = rmin
      vol(1)     = 4./3. * pi * rmin**3.

      do i=1,nbin-1
        aerad(i+1)  = aerad(i) * vrat**(athird)
        vol(i+1)    = vol(i) * vrat
      enddo

      do i=1,nbin
        rb(i+1)= ( (2.*vrat) / (vrat+1.) )**(athird) * aerad(i)
        dr(i)  = rb(i+1) - rb(i)
      enddo
      rb(nbin+1) = rbmax
      dr(nbin)   = rb(nbin+1) - rb(nbin)

      print*, ' '
      print*, ' Here are the size bins (lower bnd, center values & dr):'
      print*, ' -------------------------------------------------------'
      do i=1,nbin
         write(*,'(i2,3x,3(e12.6,4x))') i,rb(i), aerad(i),dr(i)
      enddo
       write(*,'(i2,3x,e12.6)') nbin+1,rb(nbin+1)
      print*, ' ----------------------------------------------------'

c     Now, initialize the water ice cloud radiative properties.
c     Here, we use a different size grid (more refined) called rad_rt.

      rad_rt(1)      = 1.e-7
      rad_rt(nbin_rt)= 50.e-6
      radb_rt(1)     = rbmin
      radb_rt(nbin_rt+1)= rbmax

      vrat_rt = log(rad_rt(nbin_rt)/rad_rt(1)) / float(nbin_rt-1) *3.
      vrat_rt = exp(vrat_rt)

      do i = 1, nbin_rt-1
        rad_rt(i+1)  = rad_rt(i) * vrat_rt**(athird)
        radb_rt(i+1)=((2.*vrat_rt) / (vrat_rt+1.))**(athird) * rad_rt(i)
      enddo

c     Read in the data files the Qext,Qscat and g values for each
c     size bins and for each spectral interval.

      open(60,file='data/waterCoated_vis_JD_12bands.dat')
      open(61,file='data/waterCoated_ir_JD_12bands.dat')
      open(62,file='data/Dust_vis_wolff2010_JD_12bands.dat')
      open(63,file='data/Dust_ir_wolff2010_JD_12bands.dat')

      do j = 1, nratio
      do i = 1, nbin_rt
        read(60,'(7(e12.7,x))') (qextv_cld(j,i,l), l=1,nlonv)
        read(60,'(7(e12.7,x))') (qscatv_cld(j,i,l), l=1,nlonv)
        read(60,'(7(e12.7,x))') (gv_cld(j,i,l), l=1,nlonv)
        read(61,'(5(e12.7,x))') (qexti_cld(j,i,l), l=1,nloni)
        read(61,'(5(e12.7,x))') (qscati_cld(j,i,l), l=1,nloni)
        read(61,'(5(e12.7,x))') (gi_cld(j,i,l), l=1,nloni)
      enddo
      enddo

      do i = 1, nbin_rt
        read(62,'(7(e11.5,x))') (qextv_dst(i,l), l=1,nlonv)
        read(62,'(7(e11.5,x))') (qscatv_dst(i,l), l=1,nlonv)
        read(62,'(7(e11.5,x))') (gv_dst(i,l), l=1,nlonv)
        read(63,'(5(e11.5,x))') (qexti_dst(i,l), l=1,nloni)
        read(63,'(5(e11.5,x))') (qscati_dst(i,l), l=1,nloni)
        read(63,'(5(e11.5,x))') (gi_dst(i,l), l=1,nloni)
!        factor  = qextv_dst(i,6) / (vistoir*qexti_dst(i,4))
        factor  = 1.0 
        do l = 1, nloni
          qexti_dst(i,l)  = qexti_dst(i,l)  * factor 
          qscati_dst(i,l) = qscati_dst(i,l) * factor 
        enddo
      enddo

      close(60)
      close(61)
      close(62)
      close(63)

      return
      end
