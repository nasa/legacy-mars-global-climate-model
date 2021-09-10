      subroutine initpbl

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

!     Initialize variables/constants for newpbl.f subroutine.  Variables
!     whoes names end in "_pbl" are passed, via common in standard.h,
!     into COMP3, and there are passed to newpbl in the argument list.
!     These variables are used in NEWPBL without the "_pbl" appended to
!     their names.

      use grid_h
      use defines_h
      use standard_h

      implicit none

!----------------------------------------------------------------------!

!  c-grid 25x40 configuration

!     real*4 :: zavgtg(L_J-1) = [
!    *                   164.964, 168.532, 172.986, 180.763,
!    *                   188.178, 195.825, 203.212, 209.604, 213.451,
!    *                   214.594, 213.797, 212.902, 212.568, 212.073,
!    *                   210.583, 207.491, 203.669, 198.750, 192.880,
!    *                   185.400, 177.236, 170.440, 162.364         ]

!  C-grid 36 latitude points

      real*4 :: zavgtg(L_J-1) = [
     *            163.766, 166.153, 168.532, 171.501, 175.578,
     *   180.763, 185.706, 190.727, 195.825, 200.750, 205.343,
     *   209.604, 212.169, 213.832, 214.594, 214.063, 213.499,
     *   212.902, 212.679, 212.403, 212.073, 211.080, 209.552,
     *   207.491, 204.943, 202.029, 198.750, 194.837, 190.387,
     *   185.400, 179.957, 174.971, 170.440, 165.056, 162.546 ]

!  North and south boundaries (35x60 c-grid) of GRS ground ice 

      integer :: grsn(L_I) =
     *            [28, 27, 27, 27, 27, 26, 26, 26, 27, 27, 27, 27, 27,
     *             27, 27, 27, 28, 28, 28, 29, 29, 29, 30, 30, 30, 30,
     *             29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29,
     *             29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 28, 28,
     *             28, 27, 27, 27, 27, 27, 27, 28 ]

      integer :: grss(L_I) =
     *            [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 6, 6, 7, 7, 7,
     *             7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 7,
     *             7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 8, 8, 8,
     *             8, 8, 7, 7, 7, 7 ]

!     implicit none

      integer :: i, j, k, l
      real*8  :: rmu, factl, factm, rhoc, skind

!======================================================================C

!     Numerical Constants that need to be defined for newpbl

      rmu        = 1.0
      z0_pbl     = 0.01
      epsl0_pbl  = 0.1 
      vk_pbl     = 0.4
      alphl0_pbl = 0.1

!     Calculated Constants

      ric_pbl = 0.195

!     speed constants

      dtmu_pbl   = rmu * dt * real(NC3)
      rmu1mu_pbl = (1. - rmu) / rmu

!
! SOIL MODEL INITIALISATION
!

!     Now set in input.f
!     rho = 1.5e3        ! Soil Density
!     soilcp = 627.9     ! Soil Specific Heat
 
      factl  = 0.25    ! these have to do with defining layer thickness 
      factm  = 1.2     !                    ditto
 
! Soil properties

      do i=1,l_isize
       do j=1,l_jsize-1

!  at layer mid-points

          do k=1,NL
            rhoc           = rhosoil(J,I,K)*cpsoil(J,I,K)
            scond(j,i,2*k) = (zin(j,i,k)**2)/rhoc
          end do

!  at layer boundaries

          scond(J,I,1) = scond(J,I,2)

          do k=3,2*NL-1,2
            scond(J,I,K) = 0.5*(scond(J,I,k+1)+scond(J,I,k-1))
          end do

        end do
      end do

!    Set constant depth with location

      skind = 0.06

!  Set up layer thickness (even K)

      sthick(2) = factl*skind
      do k=4,2*NL,2
        STHICK(K) = STHICK(K-2)*factm
      end do

!  Now the depth to layer boundaries (odd K)

      SDEPTH(1) = 0.0
      do K=3,2*NL+1,2
        sdepth(K) = sdepth(K-2)+sthick(k-1)
      end do

!  Now depth to layer mid-points (even K)

      do K=2,2*NL,2
        sdepth(K) = sdepth(K-1) + sthick(K)/2.0
      end do

!  Finally dz (STHICK) at layer boundaries (odd K)
!  Don't need dz at bottom boundary

      do K=3,2*NL-1,2
        sthick(K) = 0.5*(sthick(k-1) + sthick(k+1))
      end do

!  Soil model initialized

      do k = 2,2*L_LAYERS
         dxi_pbl(k) = sigma(k+3)-sigma(k+1)
      end do

!  To put water ice into the soil model, uncomment the following lines.
!  Newsoil updates - make soil properties lat/depth dependent to account
!  for ice.  Here assumed as poleward of +/- 60 degrees.  Modify values
!  only below 10 cm.

      do I=1,L_ISIZE

!  Southern hemisphere
!  GIDS - Ground Ice Depth in the Southern hemisphere

        do j=1,grss(i)
          do L=1,NL
            if(sdepth(2*L-1).gt.gids) then
              zin(j,i,L)     = 2236.995
              rhosoil(j,i,L) = 1781.99
              cpsoil(j,i,L)  = 1404.09
              rhoc           = rhosoil(J,I,L)*cpsoil(J,I,L)
              scond(j,i,2*L) = (zin(j,i,L)**2)/rhoc
            end if
          end do
 
!         at layer boundaries
 
          SCOND(J,I,1) = SCOND(J,I,2)
 
          do k=3,2*NL-1,2
            SCOND(J,I,K) = 0.5*(SCOND(J,I,K+1)+SCOND(J,I,K-1))
          end do
 
        end do

!  Northern hemisphere
!  GIDN - Ground Ice Depth in the Northern hemisphere

        do j=grsn(i),L_JSIZE-1
          do L=1,NL
            if(sdepth(2*L-1).gt.gidn) then
!              zin(j,i,L)     = 2236.995
              zin(j,i,L)     = 1100.00
              rhosoil(j,i,L) = 1781.99
              cpsoil(j,i,L)  = 1404.09
              rhoc           = rhosoil(J,I,L)*cpsoil(J,I,L)
              scond(j,i,2*L) = (zin(j,i,L)**2)/rhoc
            end if
          end do

!  at layer boundaries

          SCOND(J,I,1) = SCOND(J,I,2)

          do k=3,2*NL-1,2
            SCOND(J,I,K) = 0.5*(SCOND(J,I,K+1)+SCOND(J,I,K-1))
          end do

        end do

      end do  ! I-loop

!     Variables initialised only at COLD start

      IF(TAUI.NE.0.0) RETURN

!     independent variable (Windshear smoothed)

      do k = 2,2*L_LAYERS
         do i=1,l_isize
           do j=1,l_jsize-1
             du_pbl(j,i,k+1) = 0.0
           end do
         end do
      end do

! Temperature of soil varies from 170K to 200K 
! tanh function to accelerate soil spin up time

      do i=1,l_isize
        do j=1,l_jsize-1
          do l=1,2*nl+1
            stemp(J,I,L) = zavgtg(J)
          end do
        end do
      end do

      return
      end
