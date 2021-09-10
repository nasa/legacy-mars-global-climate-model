      subroutine dsolflux(SOL,UBAR0,GWEIGHT,FZEROV,DETAU,J1D,I1D,
     *                    DIRECTSOL)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

C  Calculate the direct surface solar flux (solar flux at the surface due
C  to the direct solar beam.
C
C  Dec 2002

      use grid_h
      use radinc_h

      implicit none
   
      real*8 FZEROV(L_NSPECTV)
      real*8 SOL(L_NSPECTV), UBAR0, GWEIGHT(L_NGAUSS), DIRECTSOL
      real*8 FACTOR
      real*8  :: detau(L_J,L_I,L_NSPECTV,L_NGAUSS)
      integer :: j1d,i1d

      integer NW, NG

C======================================================================C

      DIRECTSOL = 0.0D0

      do NW=1,L_NSPECTV
        FACTOR = UBAR0*SOL(NW)
      
        do NG=1,L_NGAUSS-1
          if(DETAU(J1D,I1D,NW,NG) .LE. 5.0) then
            DIRECTSOL = DIRECTSOL + FACTOR*
     *                  EXP(-DETAU(J1D,I1D,NW,NG)/UBAR0)*
     *                  GWEIGHT(NG)*(1.0D0-FZEROV(NW))
          end if
        end do
        NG        = L_NGAUSS
        if(DETAU(J1D,I1D,NW,NG) .LE. 5.0) then
          DIRECTSOL = DIRECTSOL + FACTOR*
     *                EXP(-DETAU(J1D,I1D,NW,NG)/UBAR0)*FZEROV(NW)
        endif
      end do

      return
      end
