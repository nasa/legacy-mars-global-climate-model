      subroutine tpindex(pw,tw,qh2o,pref,tref,wrefh2o,LCOEF,MT,MP,     &
                         NH2O,wratio)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center
 
!  PURPOSE
!    Get the TI, UI values for a 2-dimensional interpolation
!    based on the following (The interpolation is done in interpco2):
!    Interpolate the CO2 K-coefficients to the current P,T values.
!    The CO2 coefficients are given on a P,T grid:
!    P = {1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1, 1, 1E+1, 1E+2, 1E+3, 1E+4},
!    T = {50, 100, 150, 200, 250, 300, 350}.
!
!    The interpolation is the usual interpolation in 2-dimensions given
!    in "Numerical Recipes", where the "X" are P, the "Y" are
!    T, and the F(X,Y) are the CO2 K-coefficients.
!
!     The interpolating box is designated as follows:
!
!           (PL,TU)                        (PRR,TU)
!
!                          (TW,PW)
!
!           
!           (PL,TL)                        (PRR,TL)
!
!     PL  - Pressure left
!     PRR - Pressure right
!     TL  - Temperature lower
!     TU  - Temperature upper
!     PW  - Pressure wanted
!     TW  - Temperature wanted
!
!
!  INPUT PARAMETERS
!    PW                 - The pressure to interpolate to
!    TW                 - The temperature to interpolate to
!    Pref(NP)           - The pressure grid array.
!    Tref(NT)           - The temperature grid array.
!    
!  OUTPUT PARAMETERS
!    TI                 - Interpolation term (pressure)
!    UI                 - Interpolation term (temperature)
!    MT                 - Temperature index (bottom left Temperature)
!                         of bounding box
!    MP                 - Pressure index (bottom left pressure)
!                         of bounding box
!
!----------------------------------------------------------------------C

      use grid_h
      use radinc_h

      implicit none

      real*8  :: Tref(L_NTREF)
      real*8  :: pref(L_PINT)
      real*8  :: wrefh2o(L_REFH2O)

      integer :: MT, MP, N, M, NP, NH2O
      real*8  :: PW, TW, Qh2o, wratio
      real*8  :: PWL, LCOEF(4), T, U

!======================================================================C
 
!     Get the upper and lower Temperature-grid indicies that bound the
!     requested temperature.  If the requested temperature is outside
!     the T-grid, set up to extrapolate from the appropriate end.

      IF(TW.LE.TREF(1)) THEN
        MT = 1
      ELSE
        do n=1,L_NTREF-1
          if(tw.gt.Tref(n) .and. TW.LE.TREF(N+1)) then
            MT = n
            goto 10
          end if
        end do

        MT = L_NTREF-1
      
   10   continue
      END IF

      U = (TW-TREF(MT))/(TREF(MT+1)-TREF(MT))

!     Get the upper and lower Pressure-grid indicies that bound the
!     requested pressure.  If the requested pressure is outside
!     the P-grid, set up to extrapolate from the appropiate end.

      pwl = log10(pw)

      do n=2,L_PINT-1
        if(pwl.le.Pref(n)) then
          MP = n-1
          goto 20
        end if
      end do

      MP = L_PINT-1

   20 continue

      T = (PWL-PREF(MP))/(PREF(MP+1)-PREF(MP))

!  Fill the interpolation coeficients:

      LCOEF(1) = (1.0-T)*(1.0-U)
      LCOEF(2) = T*(1.0-U)
      LCOEF(3) = T*U
      LCOEF(4) = (1.0-T)*U

!  Get the indicies for water abundance.  There are 10 sets of 
!  k-coefficients with differing amounts of water vs. CO2.

      IF(Qh2o.le.WREFH2O(1)) then
        NH2O   = 1
        WRATIO = 0.0D0
      ELSEIF(Qh2o.ge.WREFH2O(L_REFH2O)) then
        NH2O   = L_REFH2O
        WRATIO = 0.0D0
      ELSE
        DO N=2,L_REFH2O
          IF(QH2O.GE.WREFH2O(N-1) .and. QH2O.lt.WREFH2O(N)) then
            NH2O   = N-1
            WRATIO = (QH2O - WREFH2O(N-1))/(WREFH2O(N) - WREFH2O(N-1))
            GOTO 30
          END IF
        END DO
      END IF

   30 CONTINUE
  
      return
      end subroutine tpindex
