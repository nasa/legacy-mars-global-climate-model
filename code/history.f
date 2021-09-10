      SUBROUTINE HISTORY

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

C
C  PURPOSE:
C      HISTORY WRITES THE ATMOSPHERIC VARIABLES FOR EACH GRID POINT TO
C      A HISTORY TAPE FILE. HISTORY IS CALLED TO STORE THE DATA FOR THE
C      ATMOSPHERE FOR THE CURRENT TIME STEP. THUS IT IS SAID: HISTORY
C      WRITES A TIME STEP TO THE HISTORY TAPE.
C  MODIFIED FROM ORIGINAL BY
C      STEVE POHORSKY    INFORMATICS     TASK 605    SEP 82
C  FOR
C      JIM POLLACK
C  ENVIRONMENT
C      SGI Indigo        IRIX 4.0.5.h    FORTRAN
C  REVISION HISTORY
C      Feburary 2, 1994 - modified to write out history tapes as
C                         real*4 values, instead of real*8.  Special
C                         startup files using real*8 values are written
C                         with the beginning of each new tape.
C                         Variables XXX4 are the real*4 versions of
C                         variable XXX.
C  INPUT PARAMETERS
C      KTP        - THE UNIT NUMBER OF THE CURRENT HISTORY TAPE FILE.
C      IBLKCT     - THE NUMBER OF TIME STEPS ALREADY WRITTEN TO THE
C                   CURRENT HISTORY TAPE FILE.
C  OUTPUT PARAMETERS
C      KTP        - UNCHANGED; OR INCREMENTED BY 1 IF STARTING A
C                   NEW HISTORY TAPE.
C      IBLKCT     - INCREMENTED BY 1; OR SET TO 1 IF STARTING A NEW
C                   HISTORY TAPE.
C  HISTORY TAPE OUTPUT
C      THE DATA FOR A GIVEN TIME STEP THAT IS WRITTEN TO THE HISTORY
C      TAPE CONSISTS OF TAU (THE CURRENT SIMULATED TIME IN HOURS) AND
C      THE ARRAYS  C, P, U, V, T, GT, TS?, TINF, FA, AND SDGR.
C      IN THE FUTURE THE QDUST ARRAY MAY ALSO BE WRITTEN. THE
C      HISTORY TAPE CONTAINS THE VARIABLE DATA FOR A TIME STEP THAT
C      WOULD BE REQUIRED BY A SUBSEQUENT RUN STARTING AT THAT TIME STEP.
C  CALLED BY
C      INPUT, MAIN
C  SUBROUTINES CALLED:
C      RETURNF FROM AMESLIB
C

      use grid_h
      use defines_h
      use fccsave_h
      use standard_h
      use comp3cmn_h
      use cldcommon_h

      implicit none

C######################################################################

      integer :: IBLKMX, IBLKMX1, KTPX, KTP8

      PARAMETER (IBLKMX=161,IBLKMX1=160)
c     PARAMETER (IBLKMX=224,IBLKMX1=223)
c     PARAMETER (IBLKMX=241,IBLKMX1=240)

C  GCM1.7 updates
C  Output files are now all written to unit 11, 51, and 91.  The 
C  output files are opened with names such as fort.11_002, fort.51_002,
C  and fort.91_002 (say), the fort.11_XXX file is always acessed as
C  unit=11, and so on.
C  unit 45 is for QTRACE

      character*12 fort11, fort51, fort91, fort45

C#=====================================================================

      fort11 = 'fort.11_0000'
      fort45 = 'fort.45_0000'
      fort51 = 'fort.51_0000'
      fort91 = 'fort.91_0000'

      IBLKCT = IBLKCT+1

      IF(IBLKCT.GT.IBLKMX) GO TO 100

   10 CONTINUE

C     Write out the variables to the history tape for this time step.
C     This is the REAL*4 history write.

C     Write out the variables to the history tape for this time step.
C     This is the REAL*4 history write.

      CALL MHISTV(TAU,P,T,U,V,GT,CO2ICE,TAUSURF,STRESSX,STRESSY,DXYP,
     *            TSTRAT,
     *            SSUN,QTRACE,QCOND,STEMP,NC3,NCYCLE,PTROP,vpout,rsdist,
     *            tofday,psf,tautot,rptau,sind,fuptopv,fdntopv,fupsurfv,
     *            fdnsurfv,fuptopir,fupsurfir,fdnsurfir,
     *            srfupflx,srfdnflx,tauref3d,DHEAT,GEOT)

      IF(.NOT. KEY15) GO TO 50
      
 50   CONTINUE

      IF(IBLKCT .NE. IBLKMX) GO TO 7000

C**********************************************
C
C     THIS TAPE IS DONE...    END IT

  100 CONTINUE

C     SET UP NEW TAPE

      close(11)

      NTAPE = NTAPE + 1

      write(FORT11(9:12),'(I4.4)') NTAPE
      write(FORT45(9:12),'(I4.4)') NTAPE
      write(FORT51(9:12),'(I4.4)') NTAPE
      write(FORT91(9:12),'(I4.4)') NTAPE

      open(11,FILE=FORT11,FORM='UNFORMATTED')
      open(45,FILE=FORT45,FORM='UNFORMATTED')
      open(51,FILE=FORT51,FORM='UNFORMATTED')
      open(91,FILE=FORT91,FORM='UNFORMATTED')

  200 CONTINUE

      IBLKCT = 1
      KTPX   = KTP
      KTP    = KTP+1
      VOLNO  = VOLNO+1
      TAUTS  = TAU
      TAUTE  = TAU+IBLKMX1*TAUH
      KTP8   = KTP+40

C     Write out the header information to the history file.

      CALL MHISTH(TOPOG,ALSP,ZIN,sdepth,dsig,dxyp,
     *            runnum,decmax,eccn,orbinc,vinc)

C     Write out the initial tape header - written on each tape
C     This is the warm start file.

      WRITE(51) runnum, tau, tauts, taute
      write(51) PNLS, PNLAT, CNLS, CNLAT, PSLS, PSLAT, CSLS, CSLAT,
     *          FUDG, PREV, ANMTSP, AORB, ASYM, DECMAX, ECCN,
     *          EGOGND, EG15GND, ETAS, ETPER, NSMTH, ORBINC, PMIN,
     *          ROT, RUNNUM, TINP, TSTART, VINC, TREFR, ALFRAY
      write(51) IMAXCN, IMAXCS, IMAXPN, IMAXPS, ISIZE, JSIZE, LAYERS,
     *          K1, LRAY, NPDST
      write(51) TOPOG, ALSP, ZIN

C     Write out the data for this time step, for warm start 
C     purposes.  The write to the history tape is done during the
C     "GO TO 10" branch.

      write(51) TAU, dclk, sind, cosd, tofday, sunstp, rotper, sdedy,
     *          sdeyr, vout, day, psf, psl, ptrop, etaout, ntape, ed,
     *          vpout, igrow, igmax, lat, dxu, dyu, dxyp, f,
     *          sinl, cosl, axu, axv, ayu, ayv, jbpn, jbcn, jbps, jbcs,
     *          conrnu, nstep, dlic, iblkct, co2latst

      write(51) FA, SDGR, DELTAT, DELTAU,
     *          DELTAV, HTHST, HSOLST, HRAD, HCOND, DMADT,
     *          PCONSAV, STRESSX, STRESSY, SD, GEOT, TS, P,
     *          U, V, T, TSTRAT, GT, CO2ICE, TINF, SSUN, FLUXSURF,
     *          RIROUT, DISGRN, TAUSURF
      write(51) FRY, HSOLCO2, HSOLDST, HTH15, HTHOUT, HCONADJ, HTURBO,
     *          HRINUM, FCONADJ, FTURBO, FRINUM, FRAYFR, HRAYFR, 
     *          DISRAY, CO2LAT
      write(51) DNDIFFV, DNIRFLUX, RHOUCH

!  FCC Save

      WRITE(51) TAUREFCLD, QEXTREFCLD, QXVCLD, QSVCLD, GVCLD,
     *          QXICLD, QSICLD, GICLD, QEXTREFDST, QXVDST,
     *          QSVDST, GVDST, QXIDST, QSIDST, GIDST, ATMCOND,
     *          LATHEAT, CUMTAUV, SCAVEFF, LON2PM, NIT,
     *          NSTEPC3, firstcall
 
      close(51)

C WRITE TSOIL AND SMOOTHED WINDSHEAR TO FORT.87:
C THIS IS NEEDED FOR WARM STARTING A RUN

      write(91) DU_PBL
      write(91) STEMP
      write(91) PIB,UOB,VOB,POB,QOB,QOC,JY1,JY2
      write(91) QTRACE
      write(91) QCOND
      write(91) QTDELTA
      write(91) NPCFLAG
      close(91)

      write(MTP,9810) KTPX,KTP,TAU
      GO TO 10

 7000 CONTINUE

      RETURN

C     FORMAT STATEMENTS

 9810 FORMAT(2X,'CHANGING FROM TAPE ',I2,' TO TAPE ',
     *       I2,' AFTER WRITING  TAU=',F15.4)
 9820 FORMAT(2X,'HISTORY RECORD WRITTEN ON TAPE ',I4,' AT TAU =',F15.5)

      END
