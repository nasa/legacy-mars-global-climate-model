      subroutine setspi(WNOI,DWNI,WAVEI)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

!     PURPOSE:
!        Set up the spectral intervals in the infrared.  Based on
!     Chris McKay's SETSPI code.
!
!     INPUT PARAMETERS
!     L_NSPECTI  - Number of spectral intervals in the INFRARED
!
!     OUTPUT PARAMETERS
!     WNOI       - Array of wavenumbers at the spectral interval
!                  centers for the infrared.  Array is NSPECTI
!                  elements long.
!     DWNI       - Array of "delta wavenumber", i.e., the width,
!                  in wavenumbers (cm^-1) of each IR spectral
!                  interval.  NSPECTI elements long.
!     WAVEI      - Array (NSPECTI elements long) of the wavelenght
!                  (in microns) at the center of each IR spectral
!                  interval.
!     
!**********************************************************************C

      use grid_h
      use constants_h, only: PI
      use radinc_h
      use radcommon_h, only: planckir

      implicit none

!     BWNI - Bin wavenumber of the edges of the IR spectral bins
!     units are inverse centimeters.  Dimension needs to be changed
!     if the number of IR bins changes.

      REAL*8  :: BWNI(L_NSPECTI+1)
      REAL*8  :: WNOI(L_NSPECTI), DWNI(L_NSPECTI), WAVEI(L_NSPECTI)

      real*8  :: a, b, ans, y, bpa, bma, T
      real*8  :: wn1, wn2
      integer :: n, nw, nt, m

!  C1 and C2 values from Goody and Yung (2nd edition)  MKS units
!  These values lead to a "sigma" (sigma*T^4) of 5.67032E-8 W m^-2 K^-4

      real*8 :: c1 = 3.741832D-16      ! W m^-2
      real*8 :: c2 = 1.438786D-2       ! m K
      
      real*8 :: x(12) = [ -0.981560634246719D0,  -0.904117256370475D0, &
                          -0.769902674194305D0,  -0.587317954286617D0, &
                          -0.367831498998180D0,  -0.125233408511469D0, &
                           0.125233408511469D0,   0.367831498998180D0, &
                           0.587317954286617D0,   0.769902674194305D0, &
                           0.904117256370475D0,   0.981560634246719D0 ]

      real*8 :: w(12) = [  0.047175336386512D0,   0.106939325995318D0, &
                           0.160078328543346D0,   0.203167426723066D0, &
                           0.233492536538355D0,   0.249147045813403D0, &
                           0.249147045813403D0,   0.233492536538355D0, &
                           0.203167426723066D0,   0.160078328543346D0, &
                           0.106939325995318D0,   0.047175336386512D0 ]

!======================================================================C

!     Bin wavenumber - wavenumber [cm^(-1)] at the edges of the IR
!     spectral bins.

      BWNI( 1) =   10.000D0       ! 1000.0 microns
      BWNI( 2) =  166.667D0       !   60.0 microns
      BWNI( 3) =  416.667D0       !   24.0 microns
      BWNI( 4) =  833.333D0       !   12.0 microns
      BWNI( 5) = 1250.000D0       !    8.0 microns
      BWNI( 6) = 2222.222D0       !    4.5 microns

!     Set up mean wavenumbers and wavenumber deltas.  Units of 
!     wavenumbers is cm^(-1); units of wavelengths is microns.

      do M=1,L_NSPECTI
        WNOI(M)  = 0.5*(BWNI(M+1)+BWNI(M))
        DWNI(M)  = BWNI(M+1)-BWNI(M)
        WAVEI(M) = 1.0E+4/WNOI(M)
      end do

!  For each IR wavelength interval, compute the integral of B(T), the
!  Planck function, divided by the wavelength interval, in cm-1.  The
!  integration is in MKS units, the final answer is the same as the
!  original planck.f; W m^-2 wavenumber^-1, where wavenumber is in CM^-1.

      DO NW=1,L_NSPECTI
        a = 1.0D-2/BWNI(NW+1)
        b = 1.0D-2/BWNI(NW)
        bpa = (b+a)/2.0
        bma = (b-a)/2.0
        do nt=500,9000
          T   = dble(NT)/1.0D+1
          ans = 0.0D0
          do m=1,12
            y    = bma*x(m)+bpa
            ans  = ans + w(m)*c1/(y**5*(exp(c2/(y*T))-1.0D0))
          end do
          planckir(NW,nt-499) = ans*bma/(PI*DWNI(NW))
        end do
      END DO

      return
      end subroutine setspi
