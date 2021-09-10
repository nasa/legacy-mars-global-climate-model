      subroutine descale(ro,rnum,rkm,rkh,rnumnum,psi)  

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

      use grid_h
      use pbl_defines_h

      implicit none

c     local arrays: none

c     global arrays:

      real*8 ro(2*n+1),rnum(2*n+1),rkm(2*n+1),rkh(2*n+1),
     *          rnumnum(2*n+1)
      real*8 psi(2*n+1,nvar)

      integer :: k, m

C======================================================================C

C  m-loop added 5/16/02  - dust tracer updates

      do m=1,nvar
        do k = 2, 2*n, 2
          psi(k,m) = psi(k,m) / ro(k)
        end do
      end do

      do k = 2, 2*n
         ro(k) = ro(k) / rnum(k)
      end do

      do k = 3, 2*n-1, 2
         rkm(k) = rkm(k) * rnumnum(k)
         rkh(k) = rkh(k) * rnumnum(k)
      end do

      return
      end
