      subroutine setrad(TGASREF,PFGASREF,CO2V,CO2I,QEXTV,QSCATV,WV,GV, &
                        QEXTI,QSCATI,WI,GI,FZEROI,FZEROV)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

!     PURPOSE:
!        Set up values used by the radiation code, such as the CO2 gas
!     absorption coefficients.  True constants are defined, and the 
!     time-independent quantities used by the radiation code are 
!     calculated. 
!
!     INPUT PARAMETERS
!     DTAU(L,M)      - Dust optical depth of layer L, and for aerosol 
!                      species M.
!     ptrop          - Pressure of the tropopause (mb)
!
!     OUTPUT PARAMETERS
!
!     AEROSOL RADIATIVE OPTICAL CONSTANTS
!     Values are at the wavelenght interval center
!
!     MIE SCATTERING - Size distribution weighted
!     Qextv    - Extinction efficiency - in the visible.
!     Qscatv   - Scattering efficiency - in the visible.
!     WV       - Single scattering albedo - in the visible.
!     GV       - Asymmetry parameter - in the visible.
!
!     Qexti    - Extinction efficiency - in the infrared.
!     Qscati   - Scattering efficiency - in the infrared.
!     WI       - Single scattering albedo - in the infrared.
!     GI       - Asymmetry parameter - in the infrared.
!     
!----------------------------------------------------------------------C

      use grid_h
      use radinc_h

      implicit none

      integer :: N, NS

      real*8  :: CO2I(L_NTREF,L_PINT,L_REFH2O,L_NSPECTI,L_NGAUSS)
      real*8  :: CO2V(L_NTREF,L_PINT,L_REFH2O,L_NSPECTV,L_NGAUSS)
      real*8  :: PFGASREF(L_PINT)
      real*8  :: PGASREF(L_NPREF), TGASREF(L_NTREF)

      real*8  :: qextv(L_NSPECTV)
      real*8  :: qscatv(L_NSPECTV)
      real*8  :: wv(L_NSPECTV)
      real*8  :: gv(L_NSPECTV)

      real*8  :: qexti(L_NSPECTI)
      real*8  :: qscati(L_NSPECTI)
      real*8  :: wi(L_NSPECTI)
      real*8  :: gi(L_NSPECTI)

!     real*8 kvis, kir
      integer :: nt, np, nw, ng

      real*8  :: fzeroi(L_NSPECTI)
      real*8  :: fzerov(L_NSPECTV)

!----------------------------------------------------------------------C

!  Visible dust properties:  M. Wolff Planck-weighted values (T=6000K)
!  Log-normal size distribution:  Reff = 1.5 microns, Veff = 0.5

!     Qext - M. Wolff values (order is increasing waveNUMBER)
!     VISULAL WAVELENGTHS.

      real*8 :: qev1(L_NSPECTV) =  [ 1.834D0, 2.296D0, 2.672D0,        &
                            2.829D0, 2.698D0, 2.452D0, 2.261D0   ]

!     Qscat - M. Wolff values
!     VISIBLE wavelengths

      real*8 :: qsv1(L_NSPECTV) =  [ 1.695D0, 2.031D0, 2.583D0,        &
                            2.744D0, 2.626D0, 2.225D0, 1.525D0   ]

!     G - M. Wolff values
!     VISIBLE wavelengths

      real*8 :: gv1(L_NSPECTV)  =  [ 0.551D0, 0.640D0, 0.661D0,        &
                            0.678D0, 0.690D0, 0.743D0, 0.868D0   ]

!     And now the INFRARED

!     M. Wolff Planck-weighted values (T=215K)
!     INFRARED wavelengths.  (The order is increasing waveNUMBER.)
!     Qext for the IR

      real*8 :: qei1(L_NSPECTI) =  [ 0.008D0, 0.262D0, 0.491D0,        &
                                              1.017D0, 0.444D0   ]

!     Qsca for M. Wolff      INFRARED wavelengths

      real*8 :: qsi1(L_NSPECTI) =  [ 0.001D0, 0.037D0, 0.122D0,        &
                                              0.351D0, 0.336D0   ]

!     g for M. Wolff values  INFRARED wavelengths

      real*8 :: gi1(L_NSPECTI)  =  [ 0.004D0, 0.030D0, 0.095D0,        &
                                              0.214D0, 0.316D0   ]

!=======================================================================

!     Set the reference pressure and temperature arrays.  These are
!     the pressures and temperatures at which we have k-coefficients.

      pgasref( 1) = 1.0E-6
      pgasref( 2) = 1.0E-5
      pgasref( 3) = 1.0E-4
      pgasref( 4) = 1.0E-3
      pgasref( 5) = 1.0E-2
      pgasref( 6) = 1.0E-1
      pgasref( 7) = 1.0
      pgasref( 8) = 1.0E+1
      pgasref( 9) = 1.0E+2
      pgasref(10) = 1.0E+3
      pgasref(11) = 1.0E+4

      tgasref(1)  =  50.0
      tgasref(2)  = 100.0
      tgasref(3)  = 150.0
      tgasref(4)  = 200.0
      tgasref(5)  = 250.0
      tgasref(6)  = 300.0
      tgasref(7)  = 350.0
 
!     Fill the (VISIBLE) arrays Qextv, Qscatv, WV, GV

      DO N=1,L_NSPECTV
        Qextv(n)  = qev1(n)
        Qscatv(n) = qsv1(n)
        IF(Qscatv(n).GE.Qextv(n)) then
          Qscatv(n) = 0.99999*Qextv(n)
        END IF
        WV(n)     = Qscatv(n)/Qextv(n)
        GV(n)     = gv1(n)
      END DO

!     Fill the (INFRARED) arrays Qexti, Qscati, WI, GI

      DO N=1,L_NSPECTI
        Qexti(n)  = qei1(n)
        Qscati(n) = qsi1(n)
        IF(Qscati(n).GE.Qexti(n)) then
          Qscati(n) = 0.99999*Qexti(n)
        END IF
        WI(n)     = Qscati(n)/Qexti(n)
        GI(n)     = gi1(n)
      END DO

!     Interpolate CO2 k coefficients to the finer pressure grid.

      call laginterp(PGASREF,PFGASREF,CO2I,CO2V,FZEROI,FZEROV)

      return
      end subroutine setrad
