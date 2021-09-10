      SUBROUTINE M4READ(KTP)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

C
C     Read the real*4 history tape first record, to position it
C     for subsequent output.
C
C     Version 0.1
C     Started Feb. 4, 1994
C
C----------------------------------------------------------------------#
      implicit none
      integer :: ktp
C======================================================================C

C     Read the header part of the file
      read(KTP)! RUNNUM4, JM4, IM4, LAYERS4, NL, NTRACE, version
      read(KTP)! DSIG4, DXYP4, GRAV4, RGAS4, cp4, stbo4, xlhtc4, kapa4,
      read(KTP)! TOPOG4, ALSP4, ZIN4

C     And now, the first record

      read(KTP)! TAU4, VPOUT4, RSDIST4, TOFDAY4, PSF4, PTROP4, TAUTOT4,
      read(KTP)! NC3, NCYCLE

      read(KTP)! P4
      read(KTP)! T4
      read(KTP)! U4
      read(KTP)! V4
      read(KTP)! GT4
      read(KTP)! CO2ICE4
      read(KTP)! STRESSX4
      read(KTP)! STRESSY4
      read(KTP)! TSTRAT4
      read(KTP)! TAUSURF4
      read(KTP)! SSUN4
      read(KTP)! QTRACE4
      read(KTP)! QCOND4
      read(KTP)! STEMP4
      read(KTP)! fuptopv, fdntopv, fupsurfv, fdnsurfv
      read(KTP)! fuptopir, fupsurfir, fdnsurfir
      read(KTP)! surfalb
      read(KTP)! dheat
      read(KTP)! geot

      read(45) ! SRFUPFLX,SRFDNFLX

      RETURN
      END
