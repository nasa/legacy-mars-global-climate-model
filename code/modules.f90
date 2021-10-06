!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                            modules.f90
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                            version_h   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module version_h

      character(len=7), parameter :: version = "V24-001"

      end module version_h
      module grid_h

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                             grid_h
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

!     X dimension of grid arrays
      integer, PARAMETER :: L_ISIZE   = 60

!     Y dimension of grid arrays
      integer, PARAMETER :: L_JSIZE   = 36

!     Number of atmospheric layers
      integer, PARAMETER :: L_LAYERS  = 24

!     Number of atmospheric levels:   2 * L_LAYERS + 3
      integer, PARAMETER :: L_LEVELS  = 2*L_LAYERS+3

      integer, PARAMETER :: L_J = L_JSIZE
      integer, PARAMETER :: L_I = L_ISIZE

!     Number of dust particle sizes
      integer, parameter :: NDP = 2

!     Number of other tracers (water, clouds. . .)

      integer, parameter :: L_NOT = 4

!     Number of tracers in GCM
!     The first NDP tracers are dust, the next L_NOT are the other 
!     tracers, water, clouds, . . .
! 
!     This MUST be greater than zero or errors appear in the loops

      integer, PARAMETER :: NTRACE = NDP + L_NOT

!     From defines.h

!     Number of soil layers
      integer, PARAMETER :: NL  = 40

      end module grid_h
      module defines_h

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                            defines_h
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use grid_h
      implicit none

!     L_LAYERS + 1
      integer, PARAMETER :: L_LAYERP1 = L_LAYERS+1

!     L_LAYERS - 1
      integer, PARAMETER :: L_LAYERM1 = L_LAYERS-1

!     L_ISIZE / 2
      integer, PARAMETER :: L_ISIZEO2 = L_ISIZE/2

!     L_ISIZE / 2
      integer, PARAMETER :: L_ISIZED2 = L_ISIZE/2

!     L_ISIZE - 1
      integer, PARAMETER :: L_ISIZEM1 = L_ISIZE-1

!     L_ISIZE - 2
      integer, PARAMETER :: L_ISIZEM2 = L_ISIZE-2

!     L_JSIZE - 1
      integer, PARAMETER :: L_JSIZEM1 = L_JSIZE-1

!     L_JSIZE - 2
      integer, PARAMETER :: L_JSIZEM2 = L_JSIZE-2

!     L_LEVELS - 1
      integer, PARAMETER :: L_LEVELM1 = L_LEVELS-1

!     L_LEVELS - 2
      integer, PARAMETER :: L_LEVELM2 = L_LEVELS-2

!     L_LEVELS - 3
      integer, PARAMETER :: L_LEVELM3 = L_LEVELS-3

!     L_LEVELS - 4
      integer, PARAMETER :: L_LEVELM4 = L_LEVELS-4

      integer, PARAMETER :: L_NPDST   = 100

      end module defines_h
      module constants_h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                            constants_h
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

!     The mathematical constant PI.
      real*8, parameter :: PI = 3.1415926535897932D0

!     Acceleration due to gravity (mks)
      real*8, parameter :: GRAV = 3.72D0

!     Gas constant for mars.
      real*8, parameter :: RGAS = 1.8902D+2

!     Stefan-Boltzmann constant
      real*8, parameter :: STBO = 5.67051D-8

!     Heat capacity (or specific heat) of CO2 gas.
!     ( units of joules per ( kg * degrees kelvin))
      real*8, parameter :: Cp = 7.3594D+2

!     Factor to convert pressures from millibars to Pascals
      real*8, parameter :: SCALEP = 1.00D+2

!     CO2 latent heat.
      real*8, parameter :: XLHTC = 5.902D+5

!     A thermodynamic constant.
      real*8, parameter :: KAPA = 0.25684D0

!     A radiation code conversion factor.
      real*8, parameter :: Cmk = 3.51D+22 

      end module constants_h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                            radinc_h
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module radinc_h

      use grid_h, only: L_LAYERS
      implicit none

!======================================================================C
!
!     RADINC.H    RADiation INCludes
!
!     Includes for the radiation code; RADIATION LAYERS, LEVELS,
!     number of spectral intervals. . .
!
!     GCM2.0  Feb 2003
! 
!======================================================================C

!     RADIATION parameters

!     In radiation code, layer 1 corresponds to the stratosphere.  Level
!     1 is the top of the stratosphere.  The dummy layer is at the same
!     temperature as the (vertically isothermal) stratosphere, and
!     any time it is explicitly needed, the appropriate quantities will
!     be dealt with (aka "top". . .)

!     L_NLEVRAD corresponds to the surface - i.e., the GCM Level that
!     is at the surface.  PLEV(L_NLEVRAD) = P(J,I)+PTROP, 
!     PLEV(2) = PTROP/2.0, PLEV(1) = ptrop/2.0

!     L_NLAYRAD is the number of radiation code layers
!     L_NLEVRAD is the number of radiation code levels.  Level N is the
!               top of layer N. 
!
!     L_NSPECTI is the number of IR spectral intervals
!     L_NSPECTV is the number of visible (or Solar) spectral intervals
!     L_NGAUSS  is the number of Gauss points for K-coefficients
!               GAUSS POINT 9 (aka the last one) is the special case
!     L_NWNGI   is L_NSPECTI*L_NGAUSS;  the total number of "intervals"
!               in the IR
!     L_NWNGV   is L_NSPECTV*L_NGAUSS;  the total number of "intervals"
!               in the VISIBLE
!
!     L_NPREF   is the number of reference pressures that the 
!               k-coefficients are calculated on
!     L_PINT    is the number of Lagrange interpolated reference
!               pressures for the CO2 k-coefficients.
!     L_NTREF   is the number of refernce temperatures for the
!               k-coefficients
!     L_TAUMAX  is the largest optical depth - larger ones are set
!               to this value.
!
!     L_REFH2O  The number of different water-mixing ratio values for
!               the k-coefficients that are now CO2+H2O. 
!
!     L_NREFI   The spectral interval number of the IR reference
!               wavelength (i.e. the 9 micron band) (8-12 microns)
!
!     L_NREFV   The spectral interval number of the visible reference
!               wavelength (i.e. the 0.67 micron band) 
!
!----------------------------------------------------------------------C

      integer, parameter :: L_NLAYRAD  = L_LAYERS+1
      integer, parameter :: L_NLEVRAD  = L_LAYERS+2
      
      integer, parameter :: L_NSPECTI =  5
      integer, parameter :: L_NSPECTV =  7
      integer, parameter :: L_NGAUSS  = 17

      integer, parameter :: L_NPREF   = 11
      integer, parameter :: L_NTREF   =  7
      integer, parameter :: L_TAUMAX  = 35

      real*8, parameter  :: MAXEXP    = 35.0D0

      integer, parameter :: L_PINT    = 51

      integer, parameter :: L_REFH2O  = 10

      integer, parameter :: L_NREFV   = 6
      integer, parameter :: L_NREFI   = 4

      end module radinc_h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                            radcommon_h
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module radcommon_h

      use grid_h, only: L_LEVELS, L_ISIZE, L_JSIZE
      use radinc_h
      implicit none

!----------------------------------------------------------------------C
!
!                             radcommon.h
!                         FORTRAN PARAMETERS
!                          GCM2.0  Feb 2003
!
!----------------------------------------------------------------------C
!
!  "Include" grid.h and radinc.h before this file in code that uses
!  some or all of this common data set
!
!     WNOI       - Array of wavenumbers at the spectral interval
!                  centers for the infrared.  Array is NSPECTI
!                  elements long.
!     DWNI       - Array of "delta wavenumber", i.e., the width,
!                  in wavenumbers (cm^-1) of each IR spectral
!                  interval.  NSPECTI elements long.
!     WAVEI      - Array (NSPECTI elements long) of the wavelenght
!                  (in microns) at the center of each IR spectral
!                  interval.
!     WNOV       - Array of wavenumbers at the spectral interval
!                  center for the VISIBLE.  Array is NSPECTV
!                  elements long.
!     DWNV       - Array of "delta wavenumber", i.e., the width,
!                  in wavenumbers (cm^-1) of each VISIBLE spectral
!                  interval.  NSPECTV elements long.
!     WAVEV      - Array (NSPECTV elements long) of the wavelenght
!                  (in microns) at the center of each VISIBLE spectral
!                  interval.
!     SOLARF     - Array (NSPECTV elements) of solar flux (W/M^2) in
!                  each spectral interval.  Values are for 1 AU, and
!                  are scaled to the Mars distance elsewhere.
!     TAURAY     - Array (NSPECTV elements) of the pressure-independent
!                  part of Rayleigh scattering optical depth.
!     PTOP       - Pressure at the top of the radiation code coordinate;
!                  = smallest k-coefficient pressure (1.0E-6 mbar)
!     FZEROI     - Fraction of zeros in the IR CO2 k-coefficients, for
!                  each temperature, pressure, and spectral interval
!     FZEROV     - Fraction of zeros in the VISIBLE CO2 k-coefficients, for
!                  each temperature, pressure, and spectral interval
!
!     AEROSOL RADIATIVE OPTICAL CONSTANTS
!     Values are at the wavelenght interval center
!
!     MIE SCATTERING - Size distribution weighted
!     Qextv    - Extinction efficiency - in the visible.
!     QextREF  - Reference visible wavelength (.67 micron band)
!     Qscatv   - Scattering efficiency - in the visible.
!     WV       - Single scattering albedo - in the visible.
!     GV       - Asymmetry parameter - in the visible.
!
!     Qexti    - Extinction efficiency - in the infrared.
!     Qscati   - Scattering efficiency - in the infrared.
!     WI       - Single scattering albedo - in the infrared.
!     GI       - Asymmetry parameter - in the infrared.
!     

      REAL*8 WNOI(L_NSPECTI), DWNI(L_NSPECTI), WAVEI(L_NSPECTI)
      REAL*8 WNOV(L_NSPECTV), DWNV(L_NSPECTV), WAVEV(L_NSPECTV)
      REAL*8 SOLARF(L_NSPECTV), TAURAY(L_NSPECTV)

      real*8 CO2I(L_NTREF,L_PINT,L_REFH2O,L_NSPECTI,L_NGAUSS)
      real*8 CO2V(L_NTREF,L_PINT,L_REFH2O,L_NSPECTV,L_NGAUSS)
      real*8 FZEROI(L_NSPECTI)
      real*8 FZEROV(L_NSPECTV)
      real*8 PGASREF(L_NPREF), TGASREF(L_NTREF)
      real*8 DETAU(L_JSIZE,L_ISIZE,L_NSPECTV,L_NGAUSS)

      real*8 qextv(L_NSPECTV), qscatv(L_NSPECTV), wv(L_NSPECTV)
      real*8 gv(L_NSPECTV)
      real*8 QextREF

      real*8 qexti(L_NSPECTI), qscati(L_NSPECTI), wi(L_NSPECTI)
      real*8 gi(L_NSPECTI)

      real*8 planckir(L_NSPECTI,8501)

      real*8 PTOP, TAUREF(L_LEVELS+1)
      real*8 PFGASREF(L_PINT)

      real*8, parameter :: UBARI = 0.5D0

!     These are for the Gauss-split 0.95 case

      real*8, parameter :: GWEIGHT(L_NGAUSS) = [                       &
                      4.8083554740D-02, 1.0563099137D-01,              &
                      1.4901065679D-01, 1.7227479710D-01,              &
                      1.7227479710D-01, 1.4901065679D-01,              &
                      1.0563099137D-01, 4.8083554740D-02,              &
                      2.5307134073D-03, 5.5595258613D-03,              &
                      7.8426661469D-03, 9.0670945845D-03,              &
                      9.0670945845D-03, 7.8426661469D-03,              &
                      5.5595258613D-03, 2.5307134073D-03,  0.0D0 ]  

!     These are for the CO2+H2O k-coefficients

      real*8, parameter :: WREFCO2(L_REFH2O) = [                       &
                     9.999999D-1, 9.99999D-1, 9.9999D-1, 9.999D-1,     &
                     9.99D-1, 9.9D-1, 9.0D-1, 8.0D-1, 7.0D-1, 6.0D-1 ]

      real*8, parameter :: WREFH2O(L_REFH2O) = [                       &
                     1.0D-7, 1.0D-6, 1.0D-5, 1.0D-4, 1.0D-3, 1.0D-2,   &
                     1.0D-1, 2.0D-1, 3.0D-1, 4.0D-1                  ]

!     If the CO2 optical depth (top to the surface) is less than
!     this value, we place that Gauss-point into the "zeros"
!     channel.  NRC parameter.

!  TLIMITS - TLIMIT for solar part of the spectrum
!  TLIMITI - TLIMIT for the IR

      real*8, parameter :: TLIMITS = 1.0D-3
      real*8, parameter :: TLIMITI = 5.0D-3

      end module radcommon_h
      module fccsave_h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                              fccsave_h
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Cloud Radiative properties

      use grid_h,   only: L_LEVELS, L_LAYERS, L_I, L_J
      use radinc_h, only: L_NSPECTV, L_NSPECTI, L_NGAUSS
      implicit none

      real*8  TAUREFCLD(L_LEVELS+1)
      real*8  QEXTREFCLD(L_LEVELS+1)

      real*8  QXVCLD(L_LEVELS+1,L_NSPECTV)
      real*8  QSVCLD(L_LEVELS+1,L_NSPECTV)
      real*8  GVCLD(L_LEVELS+1,L_NSPECTV)

      real*8  QXICLD(L_LEVELS+1,L_NSPECTI)
      real*8  QSICLD(L_LEVELS+1,L_NSPECTI)
      real*8  GICLD(L_LEVELS+1,L_NSPECTI)

!  Dust Radiative properties

      real*8  QEXTREFDST(L_LEVELS+1)

      real*8  QXVDST(L_LEVELS+1,L_NSPECTV)
      real*8  QSVDST(L_LEVELS+1,L_NSPECTV)
      real*8  GVDST(L_LEVELS+1,L_NSPECTV)

      real*8  QXIDST(L_LEVELS+1,L_NSPECTI)
      real*8  QSIDST(L_LEVELS+1,L_NSPECTI)
      real*8  GIDST(L_LEVELS+1,L_NSPECTI)

      logical :: firstcall = .true.

      real*4 ATMCOND(L_J,L_I,L_LAYERS)
      real*4 :: SCAVEFF = 0.6
      real*4 LATHEAT(L_J,L_I)

      real*8 CUMTAUV(L_J,L_I,L_NSPECTV,L_NGAUSS)

      integer LON2PM(L_I)
      integer NIT
      integer :: NSTEPC3 = 0

      end module fccsave_h
      module pbl_defines_h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                            pbl_defines_h
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use grid_h, only: L_LAYERS

      integer, parameter :: n = L_LAYERS

!      Altered for tracer scheme
!      parameter (nvar=3+ntrace)

      integer, parameter :: nvar = 4

      end module pbl_defines_h
      module standard_h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                            standard_h
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use grid_h
      use defines_h

      implicit none

!     Variables that once made up the HIST array

      real*8 ::        FA(L_JSIZE,L_ISIZE)
      real*8 ::      SDGR(L_JSIZE,L_ISIZE)
      real*8 ::    DELTAT(L_JSIZE,L_ISIZE,0:L_LAYERS)
      real*8 ::    DHEAT(L_JSIZE,L_ISIZE,L_LAYERS)
      real*8 ::    DELTAU(L_JSIZE,L_ISIZE,L_LAYERS)
      real*8 ::    DELTAV(L_JSIZE,L_ISIZE,L_LAYERS)
      real*8 ::     HTHST(L_JSIZE,L_ISIZE)
      real*8 ::    HSOLST(L_JSIZE,L_ISIZE)
      real*8 ::      HRAD(L_JSIZE,L_ISIZE)
      real*8 ::     HCOND(L_JSIZE,L_ISIZE)
      real*8 ::     DMADT(L_JSIZE,L_ISIZE)
      real*8 ::   PCONSAV(L_JSIZE,L_ISIZE)
      real*8 ::   STRESSX(L_JSIZE,L_ISIZE)
      real*8 ::   STRESSY(L_JSIZE,L_ISIZE)
      real*8 ::        SD(L_JSIZE,L_ISIZE,L_LAYERM1)

      real*8 ::      GEOT(L_JSIZE,L_ISIZE,L_LAYERS)
      real*8 ::        TS(L_JSIZE,L_ISIZE)
      real*8 ::         P(L_JSIZE,L_ISIZE)
      real*8 ::         U(L_JSIZE,L_ISIZE,L_LAYERS)
      real*8 ::         V(L_JSIZE,L_ISIZE,L_LAYERS)
      real*8 ::         T(L_JSIZE,L_ISIZE,L_LAYERS)
      real*8 ::    TSTRAT(L_JSIZE,L_ISIZE)
      real*8 ::        GT(L_JSIZE,L_ISIZE)
      real*8 ::    CO2ICE(L_JSIZE,L_ISIZE)
      real*8 ::      TINF(L_JSIZE,L_ISIZE)
      real*8 ::      SSUN(L_JSIZE,L_ISIZE)
      real*8 ::  FLUXSURF(L_JSIZE,L_ISIZE)
      real*8 ::    RIROUT(L_JSIZE,L_ISIZE)
      real*8 ::    DISGRN(L_JSIZE,L_ISIZE)

!     Variables that were collapsed from 3-d  -->  2-d are now in
!     common block HISTJL, which has one big array, HIST2D.

      real*8 ::          FRY(L_JSIZE,L_LAYERS)
      real*8 ::      HSOLCO2(L_JSIZE,L_LAYERS)
      real*8 ::      HSOLDST(L_JSIZE,L_LAYERS)
      real*8 ::        HTH15(L_JSIZE,L_LAYERS)
      real*8 ::       HTHOUT(L_JSIZE,L_LAYERS)
      real*8 ::      HCONADJ(L_JSIZE,L_LAYERS)
      real*8 ::       HTURBO(L_JSIZE,L_LAYERS)
      real*8 ::       HRINUM(L_JSIZE,L_LAYERS)
      real*8 ::      FCONADJ(L_JSIZE,L_LAYERS)
      real*8 ::       FTURBO(L_JSIZE,L_LAYERS)
      real*8 ::       FRINUM(L_JSIZE,L_LAYERS)
      real*8 ::       FRAYFR(L_JSIZE,L_LAYERS)
      real*8 ::       HRAYFR(L_JSIZE,L_LAYERS)
      real*8 ::       DISRAY(L_JSIZE,L_LAYERS)
      real*8 ::       CO2LAT(L_JSIZE,L_LAYERS)

      real*8 :: PT(L_JSIZE,L_ISIZE)
      real*8 :: UT(L_JSIZE,L_ISIZE,L_LAYERS)
      real*8 :: VT(L_JSIZE,L_ISIZE,L_LAYERS)
      real*8 :: TT(L_JSIZE,L_ISIZE,L_LAYERS)

!     This common block holds a number of other related arrays:

      real*8 :: TOPOG(L_JSIZE,L_ISIZE)
      real*8 :: PU(L_JSIZE,L_ISIZE)
      real*8 :: PV(L_JSIZE,L_ISIZE)
      real*8 :: ALSP(L_JSIZE,L_ISIZE)

!  Surfalb is the surface albedo used in the run (ALS in comp3.f)
!  and output on the history tape

      real*4  :: surfalb(L_JSIZE,L_ISIZE)

!  icealb is the albedo of surface water ice when albfeed is true.

      real*4  :: icealb

!  icethresh is the amount of surface ice required to change surface albedo

      real*4  :: icethresh_depth
      real*4  :: icethresh_kgm2

      real*8 :: SDT = 0.0
      real*8 :: DTSPLIT = 0.0

      integer :: NSPLIT
      integer :: i1d, j1d
      real*8  :: tloc1d, ylat1d, ustar, fksx, fksy, wmag

!     PBL commons - in standard.h and initpbl.f
 
      real*8  :: z0_pbl, epsl0_pbl, vk_pbl, alphl0_pbl, ric_pbl
      real*8  :: dtmu_pbl, rmu1mu_pbl,  sup_pbl, sdn_pbl
      real*8  :: dxi_pbl(2*L_LAYERS+1)
      real*8  :: du_pbl(L_JSIZE,L_ISIZE,2*L_LAYERS+1)
      real*8  :: STEMP(L_JSIZE,L_ISIZE,2*NL+1), SDEPTH(2*NL+1)
      real*8  :: STHICK(2*NL+1)
      real*8  :: SCOND(L_JSIZE,L_ISIZE,2*NL+1)

      REAL*8  :: RHOSOIL(L_J,L_I,NL), CPSOIL(L_J,L_I,NL)
      real*8  :: TAUCUM(L_LEVELS)
      real*8  :: CO2LATST(L_JSIZE)
      INTEGER :: SDEDY, SDEYR
      logical :: RESTRT
      real*8  :: LAT(L_JSIZE),     DXU(L_JSIZE)
      real*8  :: DYU(L_JSIZE),     SINL(L_JSIZE)
      real*8  :: COSL(L_JSIZE),    AXU(L_JSIZE),    AXV(L_JSIZE)
      real*8  :: AYU(L_JSIZE),     AYV(L_JSIZE),    DXV(L_JSIZE)
      real*8  :: DYV(L_JSIZE),                      F(L_JSIZE)

      integer :: NLAYM1, NLEVELS, NLEVSM1, NLEVSM2, NLEVSM3, NLEVSM4
      LOGICAL :: KEY1, KEY2, KEY8, KEY9, KEY13, KEY15, GOODINP
      real*8  :: SIGMA(L_LEVELS), TWOPI

      real*8  :: CNIN, COSA(L_ISIZE), SINA(L_ISIZE), FIM, PLAPR, PKAPA
      real*8  :: PROTCR
      real*8  :: ESQ, ANOME(1000), DAYPYR
      real*8  :: PRDST(L_NPDST),TAUDST(L_J,L_I,L_NPDST)
      real*8  :: SQRDY, ZIN(L_JSIZE,L_ISIZE,NL)

!     For the new way of calculating the ground temperature, we need 
!     IR fluxes at the surface - these are done in IN15IR and OUT15IR, 
!     and need to be saved to pass to TEMPGR.

!     DNDIFFV  - downward diffuse solar flux at the surface - added 
!                Dec 2002
!     DNIRFLUX - downward diffuse IR flux at the surface - added to 
!                common Dec 2002

      REAL*8  :: IR15DN_S(L_JSIZE,L_ISIZE), IROUTDN_S(L_JSIZE,L_ISIZE)
      REAL*8  :: RHOUCH(L_JSIZE,L_ISIZE)
      REAL*8  :: DNDIFFV(L_JSIZE,L_ISIZE), DNIRFLUX(L_JSIZE,L_ISIZE)

!     Spatially varying dust - TAUSURF(J,I) is the (visible 0.67 micron)
!     dust opacity at the surface for each grid point
!     TAUTOTJI, TAUPAT added 6/28/01  GCM1.7
!     TAUTOTJI(J,I) is the dust opacity at the RPTAU reference pressure
!     for each grid point.  TAUPAT(J,I) is the spatially varying scale
!     factor for each grid point that relates the global opacity 
!     (TAUTOT) to TAUTOTJI:  TAUTOTJI(J,I) = TAUTOT * TAUPAT(J,I)
!     TAUTOT can vary as a function of season, TAUPAT(J,I) is fixed for
!     now, but can vary.

      REAL*8  :: TAUSURF(L_J,L_I), TAUTOTJI(L_J,L_I), TAUPAT(L_J,L_I)

      real*8  :: DUSTOD(0:360)

!     Variable dust updates - allow code to read in opacity maps and use
!     those to fill TAUTOTJI 7/19/01
!     VDUST    - .TRUE. if this new dust prescription is to be used
!                .FALSE. if old/standard GCM dust method used
!                Currently set in BKDATA - might change to input file 
!                later
!     NVDUST   - number of dust opacity maps (i.e. Ls entries vdustls 
!                array)
!     VDUSTLS  - array of Ls values, one for each opacity map

      logical :: VDUST  
      INTEGER :: NVDUST
      REAL*8  :: VDUSTLS(1000)

!     Dust tracer updates for fixed dust - fixed dust means fixed 
!     spatially, no transport.  ACTIVE_DUST means radiatively active
!     Check comp3.f where these variables are dynamically set.

      logical :: ACTIVE_DUST

!     Water tracer flag - ACTIVE_WATER = .TRUE. means transported 
!     water is radiatively active.
!     Check comp3.f where this variable is dynamically set.

      logical :: ACTIVE_WATER

!     Flag to turn on/off call to microphys.f
!     Check comp3.f where this variable is dynamically set.

      logical :: MICROPHYSICS 

      logical :: CLOUDON 

      logical :: CO2SCAV 

      logical :: TIMESPLIT

      logical :: ALBFEED

      logical :: VARY_CONR 

!     Surface dust lifting: kg m^-2 s^-1
!     Surface dust downward flux (Gravitational settling flux, at
!     the surface): kg m^-2 s^-1

!     REAL*8  :: SLIFTING(L_J,L_I), SFLUXG(L_J,L_I)
      REAL*8  :: SFLUXG(L_J,L_I)

!     variables that were in bkdata.f
!
!   ALICEN  :  SURFACE ALBEDO FOR CO2 ICE ON GROUND - Northern cap
!   ALICES  :  SURFACE ALBEDO FOR CO2 ICE ON GROUND - Southern cap
!
!   CP      :  HEAT CAPACITY (OR SPECIFIC HEAT) OF CO2 GAS.
!                  ( UNITS OF JOULES PER ( KG * DEGREES KELVIN))
!
!   DAY     :  DEFAULT NUMBER OF HOURS PER DAY;  CHANGED IN EQWID
!                  TO NUMBER OF EARTH SECONDS PER MARTIAN DAY.
!   DECMAX  :  MAXIMUM SOLAR DECLINATION.
!   DLAT    :  SPACING OF GRID POINTS (DEGREES LAT.)
!   DSIG    :  COORDINATE FOR ATM. LAYERS.
!   DTM     :  UNROUNDED DEFAULT NUMBER OF MINUTES PER TIME STEP.
!
!   ECCN    :  ORBITAL ECCENTRICITY FOR MARS.
!   EGOGND  :  EMISSIVITY OF BARE GROUND OUTSIDE 15 MICRON BAND WIDTH.
!   EG15GND :  EMISSIVITY OF BARE GROUND INSIDE 15 MICRON BAND WIDTH.
!   EG15CO2 :  EMISSIVITY OF GROUND WITH CO2 ICE INSIDE 15 MICRON BAND
!                   WIDTH.
!   EPS     :  SHEAR < EPS  MEANS  NEAR ZERO SHEAR.
!
!   GRAV    :  ACCELERATION DUE TO GRAVITY (MKS)
!
!   IM      :  DEFAULT LONGITUDINAL DIMENSION OF GRID
!   ISIZE   :  X DIMENSION OF GRID ARRAYS
!
!   JM      :  DEFAULT LATITUDINAL DIMENSION OF GRID
!   JSIZE   :  Y DIMENSION OF GRID ARRAYS
!
!   KAPA    :  A THERMODYNAMIC CONSTANT.
!
!   LAYERS  :  Z DIMENSION OF GRID ARRAYS
!
!   NCYCLE  :  DEFAULT NUMBER OF TIME STEPS PER TIME-STEPPING CYCLE.
!                   (USED IN STEP.)
!   NC3     :  THE FULL COMP3 ROUTINE IS EXECUTED ONLY ONCE PER NC3
!                   TIME STEPS.
!   NLAY    :  DEFAULT NUMBER OF ATMOSPHERIC LAYERS
!   NSIZE2D :  SIZE OF 2-D GRID ARRAYS =ISIZE*JSIZE.
!   NTAPE   :  Number of tape volume (ex. fort.11_032   NTAPE = 32)
!              GCM1.7 update
!
!   PI      :  THE MATHEMATICAL CONSTANT PI.
!   PMIN    :  MIN. PRESSURE VALUE FOR EQSTD TABLE
!   PSL     :  DEFAULT PRESSURE AT A MEAN 'SEA LEVEL' SURFACE ELEVATION.
!   PTROP   :  DEFAULT PRESSURE AT THE TROPOPAUSE.
!
!   RAD     :  RADIUS OF MARS IN KM. CONVERTED TO METERS IN
!                   SUBROUTINE INPUT.
!   RGAS    :  GAS CONSTANT FOR MARS.
!   ROT     :  ANGULAR (ROTATIONAL) VELOCITY OF MARS IN RADIANS
!                   PER SECOND.
!   ROTPER  :  ROTATIONAL PERIOD IN MARTIAN HOURS.
!   RPTAU   :  Reference Pressure optical depth;  6.1 mbar for now
!
!   SCALEP  :  MULTIPLY BY 100 TO CONVERT PRESSURE FROM MILLIBARS
!                   TO PASCALS.
!   STBO    :  STEFAN-BOLTZMAN CONSTANT.
!
!   TAUCRT  :  ICE-CLOUD OPTICAL DEPTH
!   TAUH    :  DEFAULT NUMBER OF HOURS BETWEEN WRITING OUTPUT TO HISTORY
!                     TAPE.
!   TAUID   :  DEFAULT STARTING DAY FOR A GCM RUN.
!   TAUIH   :  DEFAULT STARTING HOUR FOR A GCM RUN.
!   TAUTOT  :  TOTAL GLOBAL OPTICAL DEPTH
!   TDINCR  :  MINIMUM AND INCREMENT VALUES
!   TDMIN   :  FOR SETS OF DUST OPT. DEPTHS.
!   THD     :  THERMAL DIFFUSIVITY
!   TIINCR  :  MINIMUM AND INCREMENT VALUES
!   TIMIN   :  FOR SET OF ICE OPT. DEPTHS.
!
!   VDUST   :  LOGICAL FLAG - TRUE IF WE USE OPACITY MAPS TO DEFINE
!              DUST OPACITY, FALSE IF WE USE THE STANDARD OLD GCM
!              "TAUTOT" OPACITY METHOD.
!
!   XLHTC   :  CO2 LATENT HEAT.
!
!   end original bkdata.f comments

!     THE C ARRAY IS A SMALLER RECORD-TYPE STRUCTURE USED TO GROUP
!     TOGETHER A NUMBER OF GLOBALLY USED CONSTANTS.

      logical :: DCLK
      integer :: KTP, MTP, ID
      integer :: IGROW, JBPN, JBCN
      integer :: JBPS, JBCS, IDAY, NSTEP, IBLKCT
      real*8  :: TAU, TAUI, DT, DLON, RSDIST, SIND, COSD, TOFDAY
      real*8  :: VOUT, ETAOUT, ZMM, VPOUT

      real*8  :: DAY         = 24.
      real*8  :: DLAT        = 5.0
      real*8  :: DLIC        = 90.0

      real*8, parameter :: DSIG(L_LAYERS) = [                          &
                  0.0001237D0, 0.0002186D0, 0.0003877D0, 0.0006877D0,  &
                  0.0012199D0, 0.0021622D0, 0.0038572D0, 0.0057287D0,  &
                  0.0094287D0, 0.0155443D0, 0.0256262D0, 0.0422609D0,  &
                  0.0696418D0, 0.1148966D0, 0.1530633D0, 0.1572269D0,  &
                  0.1352191D0, 0.1212412D0, 0.0718975D0, 0.0435678D0,  &
                  0.014D0,     0.008D0,     0.003D0,     0.001D0      ]

      real*8  :: dxyp(L_J)

!     Annual run parameters

      integer, parameter :: JEQUATOR = 18          ! For 36 lat points

      real*8  :: PSFGND   = 0.0000

!  Best-fit values when ground ice is used

      real*8  :: ALICEN   = 0.600
!      real*8  :: ALICEN   = 0.500
      real*8  :: ALICES   = 0.500
!      real*8  :: EGOCO2N  = 1.000
      real*8  :: EGOCO2N  = 0.800
      real*8  :: EGOCO2S  = 1.000
!      real*8  :: EG15CO2N = 1.000
      real*8  :: EG15CO2N = 0.800
      real*8  :: EG15CO2S = 1.000

!  Best-fit Values when no ground ice is used

!     real*8  :: ALICEN   = 0.700
!     real*8  :: ALICES   = 0.500
!     real*8  :: EGOCO2N  = 0.500
!     real*8  :: EGOCO2S  = 0.700
!     real*8  :: EG15CO2N = 0.500
!     real*8  :: EG15CO2S = 0.700

!     End of annual run parameters

      real*8  :: MWCO2 = 4.41D-2
      real*8  :: MWH2O = 1.80153D-2
      real*8  :: MWRATIO

      real*8  :: ANMTSP    = 73.809
      real*8  :: AORB      = 1.52369
      real*8  :: ASYM      = 0.16374
      real*8  :: DECMAX    = 25.2193
      real*8  :: ECCN      = 0.093379
      real*8  :: EGOGND    = 1.0
      real*8  :: EG15GND   = 1.0
      real*8  :: ETAS      = 49.3391
      real*8  :: ETPER     = 335.538
      real*8  :: NSMTH     = 0.2
      real*8  :: ORBINC    = 1.84987
      real*8  :: PMIN      = 0.12
      real*8  :: ROT       = 7.0882E-5
      real*8  :: RUNNUM
      real*8  :: TINP
      real*8  :: TSTART    = 10.0
      real*8  :: VINC      = 199.075      ! nominal
!     real*8  :: VINC      = 19.075       ! at Ls=180

      real*8  :: PNLS(6)   = [ 0., 40., 67., 77., 90., 100.]
      real*8  :: PNLAT(6)  = [60., 65., 75., 80., 85., 90.]
      real*8  :: CNLS(20)  = [0.,8.,17.,25.,36.,42.,60.,135.,152.,     &
                              156.,160.,163.,172.,187.,222.,280.,320., &
                              345.,355.,360. ]

      real*8  :: CNLAT(20) = [60.,65.,70.,75.,80.,85.,90.,90.,80.,75., &
                              70.,65.,60.,55.,50.,45.,45.,50.,55.,60.  ]


      real*8  :: PSLS(9)   = [163.,180.,197.,213.,225.,237.,251.,270., &
                              280.]
      real*8  :: PSLAT(9)  = [50.,55.,60.,65.,70.,75.,80.,85.,90.]
      real*8  :: CSLS(24)  = [0.,80.,135.,152.,162.,170.,176.,180.,    &
                              183.,186.,188.,190.,240.,310.,315.,320., &
                              325.,328.,330.,332.,335.,338.,340.,360.  ]

      real*8  :: CSLAT(24) = [44.,40.,40.,45.,50.,55.,60.,65.,70.,75., &
                              80.,85.,90.,90.,85.,80.,75.,70.,65.,60., &
                              55.,50.,45.,44.   ]

      integer :: IMAXCN    = 20
      integer :: IMAXCS    = 24
      integer :: IMAXPN    = 6
      integer :: IMAXPS    = 9

      integer :: NPDST     = L_NPDST
      integer :: ISIZE     = L_ISIZE
      integer :: JSIZE     = L_JSIZE
      integer :: LAYERS    = L_LAYERS

!     EMISS DATA

      real*8  :: FUDG      = 0.5
      real*8  :: PREV      = 686.980
      real*8  :: TREFR     = 0.5
      real*8  :: ALFRAY    = 2.0
      integer :: K1        = 4
      integer :: LRAY      = 3

!     remainder of bkdata

      real*8  :: DTM     = 4.8
      real*8  :: ED      = 10.0

      integer :: IDUSTSW = 1
      integer :: IGMAX   = 250
      integer :: IM      = L_ISIZE
      integer :: JM      = L_JSIZE
      integer :: MPRINT  = 4
      integer :: NCYCLE  = 15
      integer :: NC3     = 5
      integer :: NICETOP = 29
      integer :: NLAY    = L_LAYERS
      integer :: NTAPE   =  1
      real*8  :: PSF     = 6.10
      real*8  :: PSL     = 6.10
      real*8  :: PTROP   = 0.00050
      real*8  :: RAD     = 3393.0
      real*8  :: RIBMIN  = -3.0E+5
      real*8  :: ROTPER  = 24.0
      real*8  :: RPTAU   = 6.1
      real*8  :: SUNSTP  = 1.0
      real*8  :: TAUC    = 5.0
      real*8  :: TAUD    = 24.0
      real*8  :: TAUE    = 600.0
      real*8  :: TAUH    = 1.5
      real*8  :: TAUID   = 0.0
      real*8  :: TAUIH   = 0.0
      real*8  :: TAUO    = 36.0
      real*8  :: TAUO2   = 6.0
      real*8  :: TAURUN  = 25.0
      real*8  :: TAUTOT  = 0.0

!  rename on warm-start
!  rename = 1:  set tau =0 and ntape=1 on warm start

      integer :: rename = 0

!  Standard value = 0.03

      real*8  :: CONRNU = 0.03 
!     DATA  CONRNU / 0.03  /   ! Standard value  ~25km half-height
!     DATA  CONRNU / 0.003 /   ! ~50 km half-height
!     DATA  CONRNU / 0.5   /   ! ~10 km half-height

!  c-grid commons/ variables
! Topography
      real*4  :: phs(l_isize,l_jsize)
! Prognostics
      real*4  :: pib(l_isize,l_jsize,2)
      real*4  :: uob(l_isize,l_jsize,l_layers,2)
      real*4  :: vob(l_isize,l_jsize,l_layers,2)
      real*4  :: pob(l_isize,l_jsize,l_layers,2)
      real*4  :: qob(l_isize,l_jsize,l_layers,NTRACE,2)
      real*4  :: qoc(l_isize,l_jsize,NTRACE,2)

! Q array that can be passed to physics routines
! (NOTE Y,X,Z indexing for consistency with GCM physics)

      real*8  :: qtrace(l_jsize,l_isize,l_layers,ntrace)

! Q tendency calculated in COMP3 (does similar job to DELTAT, DELTAU etc)

      real*8  :: QTDELTA(L_JSIZE,L_ISIZE,L_LAYERS,NTRACE)
      real*8  :: QCDEL(L_JSIZE,L_ISIZE,NTRACE)
      real*8  :: DSTFALL(L_JSIZE,L_ISIZE,NTRACE)
      real*8  :: DSTLIFT(L_JSIZE,L_ISIZE,NTRACE)
      real*8  :: SRFUPFLX(L_JSIZE,L_ISIZE,NTRACE)
      real*8  :: SRFDNFLX(L_JSIZE,L_ISIZE,NTRACE)
      real*8  :: TAUREF3D(L_JSIZE,L_ISIZE,L_LEVELS)
      real*8  :: TSAVE(L_JSIZE,L_ISIZE,L_LEVELS)
      real*8  :: QSAVE(L_JSIZE,L_ISIZE,L_LEVELS,NTRACE)

! QCOND is boundary value of QTRACE: eg: amount of water on ground
! Also keep track of globally averaged amounts

      real*8  :: QCOND(L_JSIZE,L_ISIZE,NTRACE)

! 10-07-11 latent heat update

      real*8  ::qcondsave(L_JSIZE,L_ISIZE,NTRACE)
! end

      real*4  :: dysig(l_layers+1)
      logical :: diag
      integer :: jy1, jy2, iye(l_isize*l_jsize), iyw(l_isize*l_jsize)

      REAL*4 taucld2pm(L_JSIZE,L_ISIZE),taudst2pm(L_JSIZE,L_ISIZE)

      LOGICAL NPCFLAG(L_JSIZE,L_ISIZE)
      real*8, parameter :: npcwikg = 0.0D0  !north polar cap water ice in kg

!  GRS ground ice
!  GIDN (Ground Ice Depth North) - the depth below the surface where the 
!  water ice ground ice begins in the northern hemisphere; in meters.
!  GIDS (Ground Ice Depth South) - the depth below the surface where the 
!  water ice ground ice begins in the southern hemisphere; in meters.

      real*8, parameter  :: gidn = 0.0545
      real*8, parameter  :: gids = 0.0805

      end module standard_h
      module comp3cmn_h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                            comp3cmn_h
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use grid_h, only: L_J, L_I, L_LAYERS, L_LEVELS
      implicit none

      logical :: WANTIT

      integer :: ICMN, IM1, IP1, JCMN, JM1, JP1, IMOVER2, IMOVER4
      integer :: JEQ, IT1, JML1, JML2

      real*8  :: NDT
      real*8  :: ACOSZ, ALS, BBG, BBP(L_LEVELS), COSZ
      real*8  :: DTETA(L_LEVELS), EMGOUT, EMG15
      real*8  :: FLUXES(L_LEVELS), H(L_LEVELS)
      real*8  :: OM(L_LEVELS), PCON, PDX, PDY
      real*8  :: PL(L_LEVELS), PLOGADJ(L_LEVELS), PSAT, SCOSZ
      real*8  :: SDUST(L_LEVELS)
      real*8  :: SUNTOT(L_LEVELS), TAUICEL(L_LEVELS)
      real*8  :: TAUICET(L_LEVELS), TECON, TETA(L_LEVELS)
      real*8  :: TETASAV(L_LEVELS), TG, TL(L_LEVELS), TSAT
      real*8  :: TTDELTA(L_LAYERS), UPI(L_LEVELS)
      real*8  :: UPISAV(L_LEVELS), UTDELTA(L_LAYERS)
      real*8  :: VPI(L_LEVELS), VPISAV(L_LEVELS)
      real*8  :: VTDELTA(L_LAYERS), YM(L_LEVELS)

      real*8  :: IR15DN(L_LEVELS), IR15UP(L_LEVELS)
      real*8  :: IR15NET(L_LEVELS), IR15(L_LEVELS)
      real*8  :: IROUTDN(L_LEVELS), IROUTUP(L_LEVELS)
      real*8  :: IRONET(L_LEVELS), IROUT(L_LEVELS)
      real*8  :: IRTOTAL(L_LEVELS), IRNET(L_LEVELS)

      real*8  :: SPUPI(L_LEVELS),  SPVPI(L_LEVELS)
      real*8  :: NPUPI(L_LEVELS),  NPVPI(L_LEVELS)

      real*8  :: AADJ(L_LEVELS), BADJ(L_LEVELS)
      real*8  :: RICHN(L_LEVELS)
      real*8  :: RAYK(L_LAYERS)

      real*8  :: WTM, WPTM

!  For history tape output - radiation variables
!  History file output:  comp3cmn_h variables

      real*4  :: fuptopv(L_J,L_I), fdntopv(L_J,L_I), fupsurfv(L_J,L_I)
      real*4  :: fdnsurfv(L_J,L_I), fuptopir(L_J,L_I)
      real*4  :: fupsurfir(L_J,L_I), fdnsurfir(L_J,L_I)
!
      real*8  :: TAUTS, TAUTE
      integer :: VOLNO

      end module comp3cmn_h
      module dtcommon_h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                            dtcommon_h
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  VARIABLE_DT     - the "_DT" is used on all the "dust tracer" variables
!
!  VARIABLES:
!
!  NLAY_DT         - number of layers in the model (GCM: L_LAYERS)
!  NLEV_DT         - number of levels in the model (GCM: L_LEVELS)
!  JM_DT           - number of GCM latitude grid points
!  IM_DT           - number of GCM longitude grid points
!  J1_DT, J2_DT    - latitude indicies of specified dust lifting source 
!                    region
!  I1_DT, I2_DT    - longitude indicies of specified dust lifting source 
!                    region
!
!                       J2_DT  +-------------+
!                              |             |
!                       J1_DT  +-------------+
!                              I1_DT         I2_DT
!
!  AREA_DT         - surface area of the planet (m^2)
!  SSAREA_DT       - area of dust lifting source region (m^2)
!  DPDEN_DT        - dust particle density  (Kg/m^3)
!  DRAD_DT(NDP_DT) - dust particle bin radius
!  DSD_DT          - Dust storm duration: the total duration (seconds) 
!                    of source lifting
!  DXYP_DT(JM)     - area of GCM grid box at each latitude (m^2)
!  FLUX_DT(NDP_DT) - dust flux into the atmosphere (kg/m^2/s) from the 
!                    surface, from a specified source region
!  GRAV_DT         - acceleration due to gravity  (3.72 m/s^2)
!  GTAU_DT         - global dust optical due to dust lifting
!  NDP_DT          - number of dust particle sizes
!  NOT_DT          - Number of Other Tracers (such as water, clouds. . .)
!  RGAS_DT         - gas constant for mars
!  SCALE_DT        - 3.0*(1.38066E-23/1.66054E-27)/44.0, that is
!                    3 * Kb / (AMU * Mco2) where Kb is Boltzmann's
!                    constant, AMU the atomic mass unit, Mco2 the mass
!                    of a CO2 molecule.
!  SPECIFIED_DT    - logical variable:  true if specified source is to 
!                    be used for dust lifting, false for fully 
!                    interactive dust lifting
!  THRESHOLD_DT    - interactive source threshold value = 22.5E-3 N/m^2
!
!----------------------------------------------------------------------!

      use grid_h, only: L_J, ndp
      implicit none

      integer :: nlay_dt, nlev_dt, JM_dt, IM_dt 
      integer :: ndp_dt, not_dt
      integer :: J1_dt, J2_dt, I1_dt, I2_dt

      real*8  :: rgas_dt, dpden_dt, drad_dt(NDP), dr_dt(NDP)
      real*8  :: dxyp_dt(L_J), area_dt, ssarea_dt, dsd_dt, threshold_dt
      real*8  :: scale_dt, gtau_dt, flux_dt(NDP), qext_dt(NDP)

      logical :: specified_dt

      end module dtcommon_h
      module cldcommon_h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                            cldcommon_h
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use grid_h, only: L_LAYERS, L_JSIZE, L_ISIZE
      use radinc_h, only: L_NSPECTI, L_NSPECTV
      implicit none

      integer, parameter :: nz    = L_LAYERS
      integer, parameter :: nbin  = 4
      integer, parameter :: ntype = 3
      integer, parameter :: naer  = 5
 
      integer iMa_dt,iNb_dt,iMa_cor
      integer iMa_vap,iMa_cld,iNb_cld

      real*8 vrat
      real*8 aerad(nbin),rb(nbin+1),dr(nbin)
      real*8 vol(nbin)

      real*8 athird
      real*8 stdv(naer),aerdens(naer)
      real*8 dpden_ice,dev_dt,dev_ice 

      real*8 nav,rgp,kbz,mh2o,vo1,m0
      real*8 desorp,surfdif,nus,d_nuc,mteta

      logical :: latent_heat 
      logical :: fullcomp3
      logical :: h2ocloudform
      real*8  :: subflux(L_JSIZE,L_ISIZE)
      real*8  :: gndice(L_JSIZE,L_ISIZE)

      integer, parameter :: nbin_rt = 20
      integer, parameter :: nratio  = 15
      integer, parameter :: nlonv   = L_NSPECTV
      integer, parameter :: nloni   = L_NSPECTI

      real*8 rad_rt(nbin_rt),radb_rt(nbin_rt+1)
      real*4 qextv_cld(nratio,nbin_rt,nlonv)
      real*4 qscatv_cld(nratio,nbin_rt,nlonv)
      real*4 gv_cld(nratio,nbin_rt,nlonv)
      real*4 qexti_cld(nratio,nbin_rt,nloni)
      real*4 qscati_cld(nratio,nbin_rt,nloni)
      real*4 gi_cld(nratio,nbin_rt,nloni)

      real*4 qextv_dst(nbin_rt,nlonv),qscatv_dst(nbin_rt,nlonv)
      real*4 qexti_dst(nbin_rt,nloni),qscati_dst(nbin_rt,nloni)
      real*4 gv_dst(nbin_rt,nlonv),   gi_dst(nbin_rt,nloni)

      real*4 taucloud(l_jsize,l_isize,2)
      real*4 taudust(l_jsize,l_isize,2)

      real*4, parameter :: cor_ratio(nratio) = [ 0.1,0.2,0.25,0.3,     &
                              0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,     &
                              0.8,0.9,0.99 ]

      end module cldcommon_h
