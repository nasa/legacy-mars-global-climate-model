      SUBROUTINE GCMLOG(CLKSW,RSETSW,LDAY,LYR)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

      use version_h
      use grid_h
      use defines_h
      use constants_h, only: PI
      use standard_h
      use cldcommon_h, only: latent_heat

      implicit none

C     Machine name - for GCM log purposes.

      CHARACTER(len=72) :: MNAME
      integer :: hostnm, hostst

      character(len=10):: gcmdate
      integer CLKSW, RSETSW, values(8), lday, lyr

!     implicit none

      integer :: i, jsday, irem, isig, is, l
      real*8  :: xjlat

C======================================================================C

      gcmdate = "  -  -    "

      if(rsetsw.eq.1) then
        open(89,file='GCM_LOG')
        write(89,'("    Cold start")')
      else
        open(89,file='GCM_LOG_WARMSTART')
        write(89,'("    Warm start")')
      end if

      call date_and_time(values=values)
      write(gcmdate(1:2),'(I2.2)') values(2)
      write(gcmdate(4:5),'(I2.2)') values(3)
      write(gcmdate(7:10),'(I4)') values(1)

      hostst = hostnm(mname)

      write(89,10) version, RUNNUM, gcmdate, mname
   10 format(' ',3x,'MARS GCM LOG',25x,'GCM2.1',
     *       ' ',3x,'Version: ',a7,//,
     *       ' ',3x,'RUN:      ',f7.2,/,
     *       ' ',3x,'Date:     ',a10,/,
     *       ' ',3x,'Computer: ',a21,/)

      write(89,12) NLAY, TAUTOT, VPOUT, JM, IM, NLAY
   12 format(' ',3x,'NLAY:',i3,3x,'TAUTOT:',f8.5,3x,'Ls:',f7.2,3x,
     *       'Configuration: [',i2,',',i2,',',i2,']')

      jsday = NINT(TAUE/24.0)

      write(89,14) PTROP, PSF, DTM, jsday
   14 format(' ',3x,'PTROP:',f9.6,3x,'PSF:',f7.3,3x,'Time step:',
     *       f5.2,3x,'Length of run:',i4,' sols')
      write(89,15) rptau, conrnu
   15 format(' ',3x,'RPTAU: ',f6.2,3x,'CONRNU: ',1pe10.3,/)

      write(89,96) ALICEN, EGOCO2N, EG15CO2N
   96 format(' ',3x,'ALICEN:',f7.4,6x,'EGOCO2N:',f7.4,6x,'EG15CO2N:',
     *           f7.4)
      write(89,98) ALICES, EGOCO2S, EG15CO2S
   98 format(' ',3x,'ALICES:',f7.4,6x,'EGOCO2S:',f7.4,6x,'EG15CO2S:',
     *           f7.4,/)

      write(89,90) (DSIG(I),I=1,5)
   90 format(' ',3x,'DSIG:  ',f9.7,4(4x,f9.7))

      irem = mod(NLAY,5)
      isig = NLAY/5

      do I=2,isig
        is = (i-1)*5+1
        write(89,18) (DSIG(L),L=IS,IS+4)
   18   format(' ',10x,f9.7,4(4x,f9.7))
      end do

      if(ISIG*5+1.lt.NLAY) then
        write(89,18) (DSIG(I),I=isig*5+1,NLAY)
      endif

      write(89,20)
   20 format(' ',3x,72('_'))
      write(89,21)
   21 format(' ',3x,72('_'),/)
      write(89,22)
   22 format(' ',3x,'Purpose:',///)
      write(89,21)
      write(89,24)
   24 format(' ',3x,'Remarks:',///)
      write(89,21)
      write(89,26)
   26 format(' ',3x,'Results:',//)
      write(89,21)
      write(89,28)
   28 format(' ',3x,'Run values:',/)
      write(89,30) RUNNUM
   30 format(' ',3x,'RUNNUM: ',f7.2)

      xjlat = DLAT*180.0/PI
      write(89,32) XJLAT, JM, IM, NLAY
   32 format(' ',3x,'  DLAT:',f7.2,7x,'JM:',i3,12x,'IM:',i3,
     *           8x,'NLAY:',i3)
      write(89,34) PSF, PSL, PTROP
   34 format(' ',3x,'   PSF:',f8.3,5x,'PSL:',f7.3,5x,'PTROP:',
     *              f9.6)
      write(89,36) DTM, TAUTOT
   36 format(' ',3x,'   DTM:',f8.3,2x,'TAUTOT:',f8.4)
      write(89,38) TAUE, TAUH, TAUID, TAUIH
   38 format(' ',3x,'  TAUE:',f8.1,4x,'TAUH:',f7.3,5x,
     *           'TAUID:',f7.3,3x,'TAUIH:',f8.1)
      write(89,40) KEY2, KEY13, KEY15
   40 format(' ',3x,'  KEY2:',L3,8x,'KEY13:',L3,9x,'KEY15:',
     *           L3)
      write(89,42) IGMAX, MPRINT, NC3, NICETOP
   42 format(' ',3x,' IGMAX:',i4,6x,'MPRINT:',i3,11x,'NC3:',
     *           I3,5x,'NICETOP:',i3)
      write(89,44) CLKSW, RSETSW, LDAY, LYR
   44 format(' ',3x,' CLKSW:',i2,8x,'RSETSW:',i2,11x,'LDAY:',
     *           i4,8x,'LYR:',i3)
      write(89,'("     VDUST: ",L1,2x,"MICROPHYSICS: ",L1,4x,
     *     "LATENT_HEAT: ",L1 )') VDUST, MICROPHYSICS, LATENT_HEAT
      write(89,21)
      write(89,46)
   46 format(' ',3x,'Baseline file:',/,' ',3x,'  Update file:')  
      write(89,21)
      write(89,50)
   50 format(' ',3x,'Run time:')


      close(89)
      return
      end
