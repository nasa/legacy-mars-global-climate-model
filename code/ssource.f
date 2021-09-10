      subroutine ssource(jcmn,icmn,j1,j2,i1,i2,ndp,n,flux_dt,fs)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

C  Dust tracer
C  Calculates the surface dust flux (FS) for a specified source
C  April 2002
C  MODIFIED    Sept. 2002
C  GCM2.0  Sept 2002
C
C  JCMN, ICMN - GCM (J,I) grid point
C  J1, J2     - latitude indicies of specified dust lifting source region
C  I1, I2     - longitude indicies of specified dust lifting source region
C
C                       J2  +-------------+
C                           |             |
C                       J1  +-------------+
C                           I1            I2
C
C  NDP          - number of dust particle bins
C  FLUX_DT(NDP) - Mass flux of dust lifted from the surface in each of
C                 the NDP dust bins
C  FS           - Total mass flux of dust lifted from the surface 
C                 for all NDP particle bin sizes (kg m^-2 s^-1)
C
C----------------------------------------------------------------------C

      implicit none

      integer n, j1, j2, i1, i2, j, jcmn, icmn, ndp
      real*8  flux_dt(ndp), fs

C======================================================================C

C  If outside the lifting area, Fs = 0.0

      if((JCMN.ge.J1 .and. JCMN.le.J2)   .and. 
     *   (ICMN.ge.I1 .and. ICMN.le.I2))        then
 
        fs = flux_dt(N)

      else

        fs = 0.0
 
      end if

      return
      end
