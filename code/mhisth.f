      SUBROUTINE MHISTH(TOPOG,ALSP,ZIN,sdepth,dsig,dxyp,
     *                  runnum,decmax,eccn,orbinc,vinc)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

C     Mintz HISTory Header - WRITE out the history header information
C     for the start of each new tape, in real*4 format.  The HISTORY
C     variables are written out in HISTV.f.
C
C     Convert GCM history tapes from ieee real*8 to ieee real*4
C     for output.  A special startup file with real*8 will also be
C     written in HISTORY every time a new history "tape" is started.
C
C     Original starting code was sgi2sun.f
C
C     Version 0.1
C     Started Feb. 3, 1994
C

      use version_h
      use grid_h
      use defines_h
      use constants_h, only: PI, GRAV, RGAS, STBO, CP, XLHTC, KAPA, 
     *                       Cmk
      use standard_h, only: ALICEN, ALICES, EGOCO2N, EGOCO2S, EG15CO2N,
     *                      EG15CO2S, JEQUATOR, NPCWIKG, NPCFLAG,
     *                      GIDN, GIDS

      implicit none

C----------------------------------------------------------------------#

      REAL*8 TOPOG(L_JSIZE,L_ISIZE), ALSP(L_JSIZE,L_ISIZE),
     *       ZIN(L_JSIZE,L_ISIZE,NL)

      REAL*4 TOPOG4(L_JSIZE,L_ISIZE), ALSP4(L_JSIZE,L_ISIZE),
     *       ZIN4(L_JSIZE,L_ISIZE,NL)

      real*8 dsig(L_LAYERS), dxyp(L_J)
      real*4 dsig4(L_LAYERS), dxyp4(L_J), runnum4, grav4, rgas4

      real*4 cp4, stbo4, xlhtc4, kapa4, cmk4

      real*8 sdepth(2*NL+1)
      real*4 sdepth4(2*NL+1)

      real*4 alicen4, alices4, egoco2n4, egoco2s4
!  Grid size

      integer jm4, im4, layers4

!  Orbital parameters
      real*8 runnum, decmax, eccn, orbinc, vinc
      REAL*4 DECMAX4, ECCN4, ORBINC4, VINC4

!  Implicit none
  
      integer i, j, l, n

!  Water cyle

      real*4  :: npcwikg4

!  Sub-surface water ice depth

      real*4  :: GIDN4, GIDS4

C======================================================================C

      JM4     = L_JSIZE
      IM4     = L_ISIZE
      LAYERS4 = L_LAYERS

      runnum4 = runnum 
      grav4   = Grav 
      rgas4   = Rgas 
      cp4     = cp
      stbo4   = stbo
      xlhtc4  = xlhtc
      kapa4   = kapa
      cmk4    = cmk

      decmax4 = 180.0*decmax/PI
      eccn4   = eccn
      orbinc4 = 180.0*orbinc/PI
      vinc4   = 180.0*vinc/PI

      alicen4  = alicen
      alices4  = alices
      egoco2n4 = egoco2n
      egoco2s4 = egoco2s

      npcwikg4 = npcwikg
      gidn4    = gidn
      gids4    = gids

      do J=1,L_J-1
        dxyp4(J) = dxyp(J)
      end do

      do L=1,L_LAYERS
        dsig4(L) = dsig(L)
      end do

      do L=1,2*NL+1
        sdepth4(L) = sdepth(L)
      end do

C     Fill in the "real*4" arrays.
!  What was the CH-array is now a set of variables, real*8 and integer

      do I=1,L_ISIZE
        do J=1,L_JSIZE-1
          topog4(J,I) = topog(J,I)
          alsp4(J,I)  = alsp(J,I)
          do N=1,NL
            zin4(J,I,N)   = zin(J,I,N)
          end do
        end do
      end do

      write(11) RUNNUM4, JM4, IM4, LAYERS4, NL, NTRACE, version
      write(11) DSIG4, DXYP4, GRAV4, RGAS4, cp4, stbo4, xlhtc4, kapa4, 
     *          cmk4, decmax4, eccn4, orbinc4, vinc4, sdepth4, alicen4,
     *          alices4, egoco2n4, egoco2s4, npcwikg4, gidn4, gids4
      write(11) TOPOG4, ALSP4, ZIN4, NPCFLAG

      RETURN
      END
