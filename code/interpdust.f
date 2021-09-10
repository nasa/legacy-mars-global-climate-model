      subroutine interpdust(xls,dustod,opacity)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

C     PURPOSE:
C        Interpolate the dust optical depth to the current Ls.  The
C     dust opacities are input at 1-degree values of Ls (0 -> 360)
C     and this subroutine does a linear interpolation to get the values
C     at the current GCM Ls.
C
C     AUTHOR
C        Jim Schaeffer
C
C     UPDATES FOR
C        Bob Haberle
C
C     ENVIRONMENT
C        CANALI - Solaris 2.5.1         FORTRAN
C
C     REVISION HISTORY
C        Original 10/20/97
C
C     INPUT PARAMETERS
C     xls     - Current GCM Ls
C     dustod  - Dust opacity history used for this run.  At one degree
C               intervals in Ls, 0 -> 360.
C
C     OUTPUT PARAMETERS
C     opacity - interpolated dust opacity
C
C     CALLED BY
C        INPUT, GMP
C
C     SUBROUTINES CALLED
C        None
C
C**********************************************************************C

      implicit none

      real*8  :: xls, dustod(0:360), opacity
      integer :: ls

C=======================================================================

      if(xls.lt.0.0) then
        ls = 0
      elseif(xls.ge.360.0) then
        ls = xls - 360.0
      else
        ls = xls
      end if

      opacity = (DUSTOD(ls+1)-DUSTOD(ls))*(XLS-LS)+DUSTOD(LS)

      return
      end
