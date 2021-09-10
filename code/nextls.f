      subroutine nextls(igrow,igmax,sdedy,sunstp,anome,eccn,vinc,xls)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

      use constants_h, only: PI

      implicit none

      real*8 anome(1000)
      integer dy, sdedy

      real*8  :: quad, si, co, x, arsin, arcos

      quad(si,co) = PI-SIGN(PI-CO,si)
      arsin(X)    = ASIN(X)
      arcos(X)    = ACOS(X)

      integer :: ig, igmax, igrow
      real*8  :: esq, xls, raddeg, vinc, anomtp, anomt, v2, v1
      real*8  :: cose, sine, sd, eccn, sunstp

C======================================================================C

      esq = SQRT(1.0-ECCN**2)
 
      ig = igrow+1
      sd = sdedy

      if(ig.gt.1) then
        if(ig.le.igmax) then
          sd = sdedy+1
        end if

        if(ig.gt.igmax) then
          sd = sdedy+sunstp
        end if
      end if
 
      dy = sd+1

      cose   = cos(anome(dy))
      sine   = sin(anome(dy))
      v1     = arcos((cose-eccn)/(1.0-eccn*cose))
      v2     = arsin(esq*sine/(1.0-eccn*cose))
      anomt  = quad(v2,v1)
      anomtp = anomt-vinc+PI*0.5
 
      if(anomtp.lt.0.0) then
        anomtp = anomtp+2.0*PI
      end if

      raddeg = 180.0/PI
      xls    = anomtp*raddeg

      if(xls.lt.0.0) then
        xls = xls+360.0
      end if
      xls = mod(xls,360.0)

      return
      end
