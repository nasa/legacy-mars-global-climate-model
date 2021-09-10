      SUBROUTINE DTRIDGL(L,AF,BF,CF,DF,XK)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

      implicit none

!     DOUBLE PRESCISION VERSION OF TRIDGL

      integer, parameter :: NMAX = 201
      integer :: L, I

      real*8  :: AF(L),BF(L),CF(L),DF(L),XK(L)
      real*8  :: AS(NMAX),DS(NMAX), X, XKB

!*    THIS SUBROUTINE SOLVES A SYSTEM OF TRIDIAGIONAL MATRIX
!*    EQUATIONS. THE FORM OF THE EQUATIONS ARE:
!*    A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = D(I)
!*    WHERE I=1,L  LESS THAN 103.
!* ..............REVIEWED -CP........

!======================================================================C

      AS(L) = AF(L)/BF(L)
      DS(L) = DF(L)/BF(L)

      DO I=2,L
        X         = 1./(BF(L+1-I) - CF(L+1-I)*AS(L+2-I))
        AS(L+1-I) = AF(L+1-I)*X
        DS(L+1-I) = (DF(L+1-I)-CF(L+1-I)*DS(L+2-I))*X
      END DO
 
      XK(1)=DS(1)
      DO I=2,L
        XKB   = XK(I-1)
        XK(I) = DS(I)-AS(I)*XKB
      END DO

      RETURN
      end subroutine dtridgl
