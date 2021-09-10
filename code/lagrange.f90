      subroutine lagrange(x, xi, yi, ans)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

!  Lagrange interpolation - Polynomial interpolation at point x 
!  xi(1) <= x <= xi(4).  Yi(n) is the functional value at XI(n).

      implicit none

      real*8 :: x, xi(4), yi(4), ans
      real*8 :: fm1, fm2, fm3, fm4

!======================================================================!

      fm1 = x - XI(1)
      fm2 = x - XI(2)
      fm3 = x - XI(3)
      fm4 = x - XI(4)

!  Get the "answer" at the requested X
 
      ans = fm2*fm3*fm4*YI(1)/                                         &
                      ((XI(1)-XI(2))*(XI(1)-XI(3))*(XI(1)-XI(4)))  +   &
            fm1*fm3*fm4*YI(2)/                                         &
                      ((XI(2)-XI(1))*(XI(2)-XI(3))*(XI(2)-XI(4)))  +   &
            fm1*fm2*fm4*YI(3)/                                         &
                      ((XI(3)-XI(1))*(XI(3)-XI(2))*(XI(3)-XI(4)))  +   &
            fm1*fm2*fm3*YI(4)/                                         &
                      ((XI(4)-XI(1))*(XI(4)-XI(2))*(XI(4)-XI(3))) 

      return
      end subroutine lagrange
