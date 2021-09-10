      SUBROUTINE SETSPV(WNOV,DWNV,WAVEV,SOLARF,TAURAY)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

!     PURPOSE:
!        Set up the spectral intervals in the visible (solar).  Based
!     on Chris McKay's SETSPV code.
!
!     INPUT PARAMETERS
!     L_NSPECTV  - Number of spectral intervals in the visible
!
!     OUTPUT PARAMETERS
!     WNOV       - Array of wavenumbers at the spectral interval
!                  center for the visible.  Array is NSPECTV
!                  elements long.
!     DWNV       - Array of "delta wavenumber", i.e., the width,
!                  in wavenumbers (cm^-1) of each visible spectral
!                  interval.  NSPECTV elements long.
!     WAVEV      - Array (NSPECTV elements long) of the wavelenght
!                  (in microns) at the center of each visible spectral
!                  interval.
!     SOLARF     - Array (NSPECTV elements) of solar flux (W/M^2) in
!                  each spectral interval.  Values are for 1 AU, and
!                  are scaled to the Mars distance elsewhere.
!     TAURAY     - Array (NSPECTV elements) of the wavelength independent
!                  part of Rayleigh Scattering.  The pressure dependent 
!                  part is computed elsewhere (OPTCV).
!
!**********************************************************************C

      use grid_h
      use constants_h, only: GRAV, SCALEP
      use radinc_h

      implicit none

!     BWNV - Bin wavenumber of the edges of the visible spectral bins
!     units are inverse centimeters.  Dimension needs to be changed
!     if the number of visible bins changes.

      REAL*8  :: WNOV(L_NSPECTV), DWNV(L_NSPECTV), WAVEV(L_NSPECTV)
      REAL*8  :: SOLARF(L_NSPECTV), TAURAY(L_NSPECTV)

      REAL*8  :: SUM, WL
      INTEGER :: N, M

!     P0      - Rayleigh scattering reference pressure in pascals.
!     GRAV    - Acceleration due to gravity (g) - MKS

      real*8  :: P0 = 9.423D+6

!     Bin wavenumber - wavenumber [cm^(-1)] at the edges of the visible
!     spectral bins.  Go from smaller to larger wavenumbers, the same as
!     in the IR.

!     2222.22D0    ->   4.50 microns
!     3087.37D0    ->   3.24 microns
!     4030.63D0    ->   2.48 microns
!     5370.57D0    ->   1.86 microns
!     7651.11D0    ->   1.31 microns
!     12500.00D0   ->   0.80 microns
!     25000.00D0   ->   0.40 microns
!     41666.67D0   ->   0.24 microns

      REAL*8 :: BWNV(L_NSPECTV+1) = [ 2222.22D0, 3087.37D0, 4030.63D0, &
                                     5370.57D0, 7651.11D0, 12500.00D0, &
                                     25000.00D0, 41666.67D0 ]

!     Solar flux within each spectral interval, at 1AU (W/M^2)
!     Sum equals 1356 W/m^2 (values from Wehrli, 1985)

      real*8 :: SOLAR(L_NSPECTV) = [ 12.7, 24.2, 54.6, 145.9, 354.9,   &
                                     657.5, 106.3 ]

!======================================================================C

!     Set up mean wavenumbers and wavenumber deltas.  Units of 
!     wavenumbers is cm^(-1); units of wavelengths is microns.

      do M=1,L_NSPECTV
        WNOV(M)  = 0.5*(BWNV(M+1)+BWNV(M))
        DWNV(M)  = BWNV(M+1)-BWNV(M)
        WAVEV(M) = 1.0E+4/WNOV(M)
      end do

!     Sum the solar flux, and write out the result.  

      sum = 0.0
      do N=1,L_NSPECTV
        SOLARF(N) = SOLAR(N)
        sum       = sum+SOLARF(N)
      end do
      write(6,'("Solar flux at 1AU = ",f7.2," W/M^2")') sum

!     Set up the wavelength independent part of Rayleigh Scattering.
!     The pressure dependent part will be computed elsewhere (OPTCV).
!     WAVEV is in microns.  There is no Rayleigh scattering in the IR.

      do N=1,L_NSPECTV
        WL        = WAVEV(N)
        TAURAY(N) = (8.7/grav)*(1.527*(1.0+0.013/wl**2)/wl**4)*        &
                     scalep/P0
      end do

      RETURN
      end subroutine setspv
