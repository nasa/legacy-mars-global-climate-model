      subroutine scldef(ro,rnum,rkm,rkh,rnumnum,pc,qc,rc,psi,qrad)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

      use grid_h
      use pbl_defines_h

      implicit none

      real*8 ro(2*n+1), rnum(2*n+1), rkm(2*n+1), rkh(2*n+1)
      real*8 rnumnum(2*n+1)

      real*8 psi(2*n+1,nvar), pc(2*n+1,nvar), qc(2*n+1,nvar)
      real*8 rc(2*n+1,nvar), qrad(2*n+3)

      integer :: k, m
 
C======================================================================C
 
      do k = 2, 2*n
         ro(k) = ro(k) * rnum(k)
      end do

C  dust tracer modification - May 2002

      do m=1,nvar
        do k=2,2*n,2
          psi(k,m) = ro(k)*psi(k,m)
        end do
      end do

      do k = 3, 2*n-1, 2
         rkm(k) = rkm(k) / rnumnum(k) 
         rkh(k) = rkh(k) / rnumnum(k) 
      end do

C  dust tracer modification - May 2002

      do m=1,nvar
        do k=2,2*n,2
          pc(k,m) = 0.0
          qc(k,m) = 0.0
          rc(k,m) = 0.0
        end do
      end do

      do k=2,2*n,2
        rc(k,3) = ro(k)*qrad(k+2)
        rc(k,4) = 0.
      end do

      do k = 3, 2*n-1, 2
        pc(k, 1) = rkm(k)
        pc(k, 2) = rkm(k)
        pc(k, 3) = rkh(k)
c        do m=1,ntrace
c          pc(k,m+3) = rkh(k)
c        end do
         pc(k,4) = rkh(k)
      end do

      return
      end
