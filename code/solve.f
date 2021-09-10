      subroutine solve(i,vect,mat,psi)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

      use grid_h
      use pbl_defines_h

      implicit none

c     local arrays:

c     global arrays:

      real*8  :: mat(n,nvar)
      real*8  :: vect(n),psi(2*n+1,nvar)
      integer :: i, j

C======================================================================C

      do j = 1, n-1
         mat(j+1, 2) = mat(j+1,2) - mat(j+1,1) * mat(j,3) / mat(j,2)
         vect(j+1)   = vect(j+1)  - mat(j+1,1) * vect(j)  / mat(j,2)
      end do

      vect(n) = vect(n) / mat(n,2)

      do j = n-1, 1, -1
         vect(j) = (vect(j) - mat(j,3) * vect(j+1)) / mat(j,2)
      end do

      do j = 1, n
         psi(2*j, i) = vect (j)
      end do

      return
      end
