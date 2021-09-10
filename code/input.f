      SUBROUTINE INPUT

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

C
C  PURPOSE:
C      INPUT READS INPUT VALUES AND RUN PARAMETERS AND CONTROLS THE
C      INITIALIZATION SEQUENCE IN THIS PROGRAM.
C  MODIFIED FROM ORIGINAL BY
C      STEVE POHORSKY    INFORMATICS     TASK 605    JUL 82
C  FOR
C      JIM POLLACK
C  ENVIRONMENT
C      Cray-2            UNICOS 3.0      FORTRAN
C  REVISION HISTORY
C
C  INPUT PARAMETERS
C      SUBROUTINE INPUT READS JOB INPUT INCLUDING RUN PARAMETERS.
C      INCLUDED AMONG THESE ARE THE STARTING DAY AND HOUR OF SIMULATED
C      TIME AT WHICH THIS RUN OF THE PROGRAM IS TO START. IF THESE ARE
C      BOTH ZERO, THEN SUBROUTINE INIT1 IS CALLED TO
C      PERFORM A COLD START OF THIS MODEL FROM A SOMEWHAT ARTIFICIAL
C      (ISOTHERMAL WITH NO WINDS) SET OF DATA. OTHERWISE THE STARTING
C      DAY AND HOUR REFER TO A TIME STEP FROM A PREVIOUS RUN OF THIS
C      MODEL FROM WHICH THIS RUN IS TO CONTINUE. FOR THIS CASE THE
C      HISTORY TAPE OUTPUT FROM THE PREVIOUS RUN IS READ, AND THE MODEL
C      IS INITIALIZED TO THE CONDITIONS IT HAD IN THE SPECIFIED
C      TIME STEP.
C  OUTPUT PARAMETERS
C      INPUT INITIALIZES THE MAJOR ARRAYS (P,U,V,T,GT,TOPOG,ALSP,...)
C      THAT SPECIFY THE STATE OF THE ATMOSPHERE AND THE SURFACE FOR THE
C      STARTING TIME STEP. (SEE SUBROUTINE INIT1.) ALSO, A NUMBER OF
C      CONSTANTS ARE CALCULATED FOR USE BY THE ROUTINES THAT PERFORM
C      THE TIME INTEGRATION. THE HISTORY TAPE FILE IS ALSO SET UP AND/OR
C      POSITIONED FOR SUBSEQUENT OUTPUT FROM THIS RUN.
C  CALLED BY
C      MAIN
C  SUBROUTINES CALLED:
C      EMISS, MAGFAC, INSDET (WHICH CALLS SDET); INIT1 CALLED
C      ONLY WHEN DOING A COLD START.
C

      use grid_h
      use defines_h
      use constants_h, only: PI, GRAV, RGAS, KAPA
      use fccsave_h
      use standard_h
      use dtcommon_h
      use cldcommon_h
      use comp3cmn_h, only: tauts, taute, volno, imover2, imover4, ndt,
     *                      wantit

      implicit none

C#######################################################################

      INTEGER :: TID

      real*8  :: BOUNDUM(L_JSIZE+1,L_ISIZE)

      CHARACTER(len=64) :: DUMMY

      integer, parameter :: coldstart = 1
      integer, parameter :: warmstart = 0
      integer, parameter :: renameyes = 1

!  Namelist input

      integer :: RSETSW
      real*8  :: runnumx
      namelist / inputnl / RUNNUMX, DLAT, JM, IM, NLAY, PSF, PTROP,
     *                     DTM, TAUTOT, RPTAU, CONRNU, TAUE, TAUH,
     *                     TAUID, TAUIH, NC3, RSETSW, RENAME,
     *                     CLOUDON, ACTIVE_DUST, ACTIVE_WATER,
     *                     MICROPHYSICS, CO2SCAV, TIMESPLIT, ALBFEED,
     *                     LATENT_HEAT,VDUST, ICEALB, ICETHRESH_DEPTH,
     *                     DTSPLIT,H2OCLOUDFORM

      integer :: i, j, k, n, ios, inu, jkl, ktpr,l,m
      real*8  :: taur, xls, taux, tauxm

!     NPCFLAG input variable

      logical :: npcflagin(L_J+1,L_I)


      !! Zero argument function that returns the number of bytes
      !! in a "recl=" unit
      integer, external :: recl_unit

      !! record length of the direct-access file
      integer, parameter :: record_length_in_bytes = 17768


C#================================================================

C     HISTORY FILE

      KTP  = 11

C     DATA CARD IMAGE FILE

      INU = 5

C     OUTPUT  STREAM

      MTP = 6

C     Read in the run values.

      READ(INU,inputnl)

      NCYCLE  = NC3
      NICETOP = 2*L_LAYERS + 3
      PSL     = PSF
      KEY2    = .TRUE.
      KEY13   = .FALSE.
      KEY15   = .TRUE.
      IGMAX   = 999
      MPRINT  = 48

C     Write out a copy of the run values to the output file.

      WRITE(MTP,7005)
      WRITE(MTP,7006) RUNNUMX
      WRITE(MTP,7007) TAUID, TAUIH, TAUD, TAUE, TAUO, TAUH, TAUC
      WRITE(MTP,7008) DTM, NCYCLE, NC3, JM, IM, DLAT
      WRITE(MTP,7009) (DSIG(JKL),JKL=1,NLAY)
      WRITE(MTP,7010) ED, RAD, GRAV, DAY, RGAS, KAPA, PSL, PSF
      WRITE(MTP,7011) PTROP, DLIC, TAUO2, TAURUN, KEY2, KEY13, KEY15
      WRITE(MTP,7012) TSTART, MPRINT, NLAY, NICETOP, IDUSTSW, IGMAX
      WRITE(MTP,7014) TAUTOT, TREFR, CONRNU, TAUR
      WRITE(MTP,7015) RPTAU, VDUST, MICROPHYSICS, RENAME
      WRITE(MTP,7016)

 7005 FORMAT(' ',//,25X,'************* RUN VALUES *************',/)
 7006 FORMAT(' ','Run number: ',f7.2,/)
 7007 FORMAT(' ','TAUID =',F5.1,3X,'TAUIH =',F8.1,3X,'TAUD =',
     *       F6.1,2X,'TAUE =',F7.1,3X,'TAUO =',F5.1,/,' TAUH =',
     *       F5.1,3X,'TAUC =',F5.1)
 7008 FORMAT(' ','DTM =',F5.1,3X,'NCYCLE =',I3,3X,'NC3 =',I3,3X,
     *       'JM =',I3,3X,'IM =',I3,3X,'DLAT =',F6.2)
 7009 FORMAT(' ','DSIG =',5(2X,F9.6))
 7010 FORMAT(' ','ED =',F5.1,3X,'RAD =',F6.0,3X,'GRAV =',F6.2,3X,
     *       'DAY =',F6.2,3X,'RGAS =',F7.2,/,' KAPA =',F8.5,3x,
     *       'PSL =',F6.2,3X,'PSF =',F6.2)
 7011 FORMAT(' ','PTROP =',F8.5,3X,'DLIC =',F6.2,3X,'TAUO2 =',F5.1,
     *       3X,'TAURUN =',F6.1,/,' KEY2 =',L3,3X,'KEY13 =',L3,3X,
     *       'KEY15 =',L3)
 7012 FORMAT(' ','TSTART =',F6.2,3X,'MPRINT =',I3,3X,'NLAY =',I3,3x,
     *       'NICETOP =',I3,/,' IDUSTSW =',I2,3X,'IGMAX =',I4)
 7014 FORMAT(' ','TAUTOT =',F8.5,3X,'TREFR =',F6.3,3X,'CONRNU =',
     *       F7.4,3X,'TAUR = ',1PE8.2)
 7015 FORMAT(' ','RPTAU = ',f6.2,3x,'VDUST = ',L1,3x,'MICROPHYS = ',
     *           L1,3x,'RENAME = ',i1)
 7016 FORMAT(' ',/,25X,'**************************************',///)

      IF(MOD(NICETOP,2).EQ.1.AND.NICETOP.LE.2*NLAY+3) GOTO 30

      WRITE (MTP,9010) NICETOP
      GOODINP = .FALSE.
      GOTO 8000

   30 CONTINUE

      TAUI    = TAUID*24.0+TAUIH
      DT      = DTM*60.0
      FIM     = IM
      DLAT    = DLAT*PI/180.0
      DLON    = 2.0*PI/FIM

!     QH2O is the mixing ratio.  The GCM computes it as mass mixing 
!     ratio.  The radiation code wants number mixing ratio.  The TOTAL 
!     mass in each layer is just the mass of CO2, not CO2+H2O.  If we 
!     change the total mass, then the expresion for MWRATIO would 
!     change to read MWRATIO = (MWCO2+MWH2O)/MWH2O

      MWRATIO = MWCO2/MWH2O

C     INITIALIZE CONSTANT ELEMENTS OF VECTOR OF CUMULATIVE
C     OPTICAL DEPTHS FOR INDIVIDUAL GRID POINTS

      TAUCUM(1) = 0.0
      TAUCUM(3) = 0.0

      do J=1,L_J-1
        do I=1,L_I
          DMADT(J,I) = 0.0
        end do
      end do

C     READ IN THERMAL INERTIA, ALTITUDE, AND ALBEDO FROM DATA FILES

      OPEN(UNIT=9,
     *  FILE='data/osu_ti_5x6_2011',
     *  STATUS='OLD',FORM='UNFORMATTED')
      READ(9) BOUNDUM
      CLOSE(9)

      DO J = 1,L_JSIZE-1
        DO I = 1,L_ISIZE
          DO K=1,NL
            ZIN(J,I,K) = BOUNDUM(J+1,I)
          END DO
        END DO
      END DO

C SOIL model 3-D arrays, RHOSOIL & CPSOIL
C We'll eventually read these values in from a data file

      DO J=1,L_JSIZE-1
        DO I=1,L_ISIZE
          DO K=1,NL
            RHOSOIL(J,I,K) = 1481.39
            CPSOIL(J,I,K)  = 840.0
          END DO
        END DO
      END DO

      OPEN(UNIT=9,
     *  FILE='data/topog37x60.mola_intel',
     *  STATUS='OLD',FORM='UNFORMATTED')
      READ(9) BOUNDUM
      CLOSE(9)

      DO J = 1,L_JSIZE
        DO I = 1,L_ISIZE
          TOPOG(J,I) = BOUNDUM(J+1,I)
        END DO
      END DO

C Initialise dynamical core geopotential WHICH IS NEGATIVE OF TOPOG FILES

      DO I = 1,L_ISIZE
        DO J = 1,L_JSIZE
          PHS(I,J) = -TOPOG(J,I)
        END DO
      END DO

      OPEN(UNIT=9,
     *  FILE='data/osu_albedo_5x6_2011',
     *  STATUS='OLD',FORM='UNFORMATTED')
      READ(9) BOUNDUM
      CLOSE(9)

      DO J = 1,L_JSIZE-1
        DO I = 1,L_ISIZE
          ALSP(J,I) = BOUNDUM(J+1,I)
        END DO
      END DO

!  NPCFLAG - template that defines where the north polar ice cap is.

      open(20,file='data/npcflag_osu_550_intel',
     *        status='old',form='unformatted')
      read(20) npcflagin
      close(20)

      do j=1,L_J-1
        do i=1,L_I
          npcflag(j,i) = npcflagin(j+1,i)
        end do
      end do


C  Read in the opacity maps if VDUST = .TRUE., otherwise read in
C  the opacity ascii file (original GCM method).  VDUST modification
C  made 7/19/01.  UNIT 85 is to remain open throughout the run.

!  opacitymap_0.3 has 3 Ls values:  359.7, 359.8, 359.9  & all
!  opacity values are 0.30000000000000

      IF(VDUST) THEN
      open(85,file=
     * 'data/TES_my24_dustscenario_zvary_37x60_6ls_intel',
!     * 'data/TES_my24_dustscenario_zmean_37x60_6ls_intel',
!     * 'data/TES_yr26_dust_opacity_37x60_intel',
     *          form='unformatted',access='direct',
     *          recl=(record_length_in_bytes / recl_unit()),
     *          status='old',iostat=ios)

        if(ios.ne.0) then
          write(6,'("Could not open data file: opacity.map")')
          stop
        end if

        do n=1,1000
          read(85,rec=n,iostat=ios) xls
!         print*,xls
          if(ios.ne.0) exit
          nvdust     = n
          vdustls(n) = xls
        end do
        
        if(nvdust.lt.1) then
          write(6,'("Problem with opacity.map data.")')
          stop
        endif

C  Get initial TAUTOTJI map (Ls < 0 returns first opacity record).
 
        xls = 0.0
        call getvdust(xls,nvdust,vdustls,TAUTOTJI)

      ELSE

C  Read in dust optical depths for each day of the run.

c       open(UNIT=4,FILE='/home/jschaef/gcm/gcmdata/opacity.ascii',
c    *       status='OLD')
c       read(4,36) dummy
c  36   format(a64)

        do n=0,360
c         read(4,*) xls, dustod(n)
          dustod(n) = TAUTOT
        end do
        close(4)

C  GCM1.7   6/28/01
C  This is where spatially varying dust (at the RPTAU reference
C  pressure) is set.  In this picture, dustod is a global scale
C  factor acting on TAUTOTJI.  TAUTOTJI(J,I) is the dust optical
C  depth at the RPTAU reference pressure (6.1 mbar is the current
C  standard reference pressure).  TAUPAT(J,I) is the spatially varying
C  scale factor, TAUTOT is the global factor (aka the viking dust
C  history and so on) and can change with time.  TAUPAT(J,I) can
C  vary in time, but for now is constant.  When we start to model
C  the TES results using the TES opacities, that will enter through
C  TAUPAT.

C  Remember, the c-grid J=1 is NOT the south pole, adjust possible
C  patterns accordingly.

        do J=1,L_J-1
          do I = 1,L_I
            TAUPAT(J,I) = 1.0D0
          end do
        end do
 
        do J=1,L_J-1
          do I = 1,L_I
            TAUTOTJI(J,I) = TAUTOT*TAUPAT(J,I)
          end do
        end do

      END IF

      CALL EMISS
      RAD = RAD*1000.0
      CALL MAGFAC

C     Initialize new pbl code
 
      CALL INITPBL

C     Initialize the dust tracer code
 
      CALL INITDT

      CALL INITCLD

      icethresh_kgm2=icethresh_depth*dpden_ice*1.e-6

!     IF(TAUI.NE.0.0)  GOTO 200
      if(rsetsw.eq.warmstart) goto 200

C     START UP FROM DATA FILES AND CONSTANTS...   COLD START
C     ISOTHERMAL ATMOSPHERE WITH NO WINDS.

      RUNNUM = RUNNUMX
      TAU    = 0.0
      CALL INIT1

      CALL INSDET(rsetsw)
      TOFDAY = 0.0

C     Because the initial TAUTOT value is anywhere in the dustod array,
C     not necessarily the first value, we need to get the Ls now, and
C     the appropiate TAUTOT value.
C  VDUST update 7/19/01:  if VDUST = .TRUE. use opacity maps, otherwise
C  use the old GCM "TAUTOT" method.

      call nextls(igrow,igmax,sdedy,sunstp,anome,eccn,vinc,xls)

      if(VDUST) then
        call getvdust(xls,nvdust,vdustls,TAUTOTJI)
      else
        call interpdust(xls,dustod,TAUTOT)

C  GCM1.7   6/28/01

        do J=1,L_J-1
          do I = 1,L_I
            TAUTOTJI(J,I) = TAUTOT*TAUPAT(J,I)
          end do
        end do
      endif

C     SET UP HISTORY TAPE

      VOLNO = 1
      TAUTS = 0.0
      TAUTE = TAUE+0.001
      KTP   = 11

C     Write out the initial tape header - written on each tape
C     Used for standard 50-day history stuff

      CALL MHISTH(TOPOG,ALSP,ZIN,sdepth,dsig,dxyp,
     *            runnum,decmax,eccn,orbinc,vinc)

      IBLKCT = 0
      IF(KEY13)GO TO 8000
      CALL HISTORY

      GOTO 8000

C     HISTORY TAPE RESTART

  200 CONTINUE

      KTPR = 51
      KTP  = 11

!     Modified warmstart logic

      READ(KTPR) runnum, taux, tauts, taute

      if(abs(taux-tauih).gt.0.01) then
        write(6,'(//)')
        write(6,'("Problem with warm start:  check requested time.")')
        stop
      end if

!     We've got the correct file - proceed with the warm start

      CALL M4READ(KTP)

      READ(KTPR) PNLS, PNLAT, CNLS, CNLAT, PSLS, PSLAT, CSLS, CSLAT,
     *           FUDG, PREV, ANMTSP, AORB, ASYM, DECMAX, ECCN,
     *           EGOGND, EG15GND, ETAS, ETPER, NSMTH, ORBINC, PMIN,
     *           ROT, RUNNUM, TINP, TSTART, VINC, TREFR, ALFRAY
      READ(ktpr) IMAXCN, IMAXCS, IMAXPN, IMAXPS, ISIZE, JSIZE, LAYERS, 
     *           K1, LRAY, NPDST
 
      READ(KTPR)

      read(ktpr) TAUx, dclk, sind, cosd, tofday, sunstp, rotper, sdedy,
     *           sdeyr, vout, day, psf, psl, ptrop, etaout, ntape, ed,
     *           vpout, igrow, igmax, lat, dxu, dyu, dxyp, f,
     *           sinl, cosl, axu, axv, ayu, ayv, jbpn, jbcn, jbps, jbcs,
     *           conrnu, nstep, dlic, iblkct, co2latst

      TAUXM = TAUX
      read(ktpr) FA, SDGR, DELTAT, DELTAU,
     *           DELTAV, HTHST, HSOLST, HRAD, HCOND, DMADT,
     *           PCONSAV, STRESSX, STRESSY, SD, GEOT, TS, P,
     *           U, V, T, TSTRAT, GT, CO2ICE, TINF, SSUN, FLUXSURF,
     *           RIROUT, DISGRN, TAUSURF
      READ(KTPR) FRY, HSOLCO2, HSOLDST, HTH15, HTHOUT, HCONADJ, HTURBO,
     *           HRINUM, FCONADJ, FTURBO, FRINUM, FRAYFR, HRAYFR,
     *           DISRAY, CO2LAT
      READ(KTPR) DNDIFFV, DNIRFLUX, RHOUCH

!  FCC Save

      READ(KTPR) TAUREFCLD, QEXTREFCLD, QXVCLD, QSVCLD, GVCLD,
     *           QXICLD, QSICLD, GICLD, QEXTREFDST, QXVDST,
     *           QSVDST, GVDST, QXIDST, QSIDST, GIDST, ATMCOND,
     *           LATHEAT, CUMTAUV, SCAVEFF, LON2PM, NIT,
     *           NSTEPC3, firstcall

C WARM START: READ TSOIL FROM FORT.91: THE RECORD CONTAINING TSOIL
C AND SMOOTHED WINDSHEAR FROM TIME T-1

      READ(91) DU_PBL
      READ(91) STEMP
C Read C grid dynamics- a bit of a waste but it's only done every
C 160 timesteps
      READ(91) PIB,UOB,VOB,POB,QOB,QOC,JY1,JY2
      READ(91) QTRACE
      READ(91) QCOND
      READ(91) QTDELTA
      READ(91) NPCFLAG
    

      do j=1,l_j
       do i=1,l_i
        do l=1,l_layers
         k=2*l+2
         tsave(j,i,k)=t(j,i,l)
         do m=1,ntrace
          qsave(j,i,k,m)=qtrace(j,i,l,m)
         enddo
        enddo
       enddo
      enddo
 
!  end warm start updates

      if(rename.eq.renameyes) then
        TAU   = 0.0
        NTAPE = 0
      else
        TAU    = TAUX
      end if

      RUNNUM = RUNNUMX
      CALL INSDET(rsetsw)

C     Get the current dust-loading for this Ls.
C     VDUST update added 7/19/01

      if(VDUST) then
        call getvdust(vpout,nvdust,vdustls,TAUTOTJI)
      else
        call interpdust(vpout,dustod,TAUTOT)

C  GCM1.7   6/28/01

        do J=1,L_J-1
          do I = 1,L_I
            TAUTOTJI(J,I) = TAUTOT*TAUPAT(J,I)
          end do
        end do
      endif

      call dustprofile  !  Added 08-07-06

      TAUE   = TAUI+TAUE
      TAUI   = TAU
      TOFDAY = MOD(TAU,ROTPER)  ! generic MOD 
      WRITE(MTP,9839) TAU
      GOTO 8000

C     TAPE ID DOESNT MATCH EXPERIMENT ID

 1000 CONTINUE

      WRITE(MTP,9835) TID,ID
      GOODINP = .FALSE.

 8000 CONTINUE

c     Check that the albedo and thermal inertia of the north polar cap
c     is not too low.

!     DO J=1,L_J
!       DO I=1,L_I
!         IF (NPCFLAG(J,I)) THEN
!            ALSP(J,I)  = max(ALSP(J,I),0.325)
!            ZIN(J,I,1) = max(ZIN(J,I,1),500.)
!            DO N=1,NL
!              ZIN(J,I,N) = ZIN(J,I,1)
!            END DO
!         ENDIF
!       ENDDO
!     ENDDO

C     Put input value of runnum into the CH array.

      RUNNUM = RUNNUMX

      RETURN

 9002 FORMAT (1X,'IGMAX=',I5,5X,'TSTART=',F10.2,5X,'TAURUN=',F7.2,
     *        5X,' MPRINT=',I4)
 9010 FORMAT( ///,1X, '#########  ERROR: ', I3,
     *                ' IS AN INVALID VALUE FOR NICETOP.  #########')
 9050 FORMAT (10A4)
 9190 FORMAT('1*****  DEFAULT INPUT VALUES  **********')
 9193 FORMAT('1********  RUNNING VALUES  *************')
 9835 FORMAT(2X, 'TAPE ID (',A4,') DOESNT MATCH EXPERIMENT ID (',A4,')')
 9837 FORMAT(2X, 'END OF TAPES FOR EXPERIMENT (', A4, ') ENCOUNTERED',/,
     *       2X, ' STARTING AT END OF VOLUME ',I2)
 9839 FORMAT(2X, 'EXPERIMENT RESTARTING FROM TAPE AT TAU=',
     *            F15.5)
 9841 FORMAT(2X, 'TAU=',F15.5,'   NOT FOUND ON TAPE (',A4,')  VOLUME ',
     *       I2, / ,
     *       2X, 'CHANGING FROM TAPE ',I2,' TO TAPE ', I2)
 9843 FORMAT(1X, '==>  HISTORY RECORD READ AT TAU=', F15.5)

      END


      integer function recl_unit()
      implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Figure out the units that the "recl=" specifier wants when opening
!! a file for direct access.  Some compilers specify recl in bytes,
!! some specify it in "words", which is usually 4 bytes, but not always.
!!
!! The return value of this function is the number of bytes that
!! "recl=1" would use.  It can thus be used to create a system
!! independent "open" specification: have the code internally specify
!! the record length in bytes, and then open the file with, e.g.:
!!      open(unit_number, filename, access='direct',
!!     *     recl= record_length_in_bytes / recl_unit() )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !! Find out how many "units" a real*8 value would take.
      !! We use real*8 because we know by definition that it is 8 bytes
      !! long, and no compiler I'm aware of counts recl in units that
      !! are bigger than that.
      real*8  :: fake
      integer :: num_units

      inquire(iolength=num_units) fake

      if (num_units .eq. 1) then
          recl_unit = 8
      else if (num_units .eq. 2) then
          recl_unit = 4
      else if (num_units .eq. 4) then
          recl_unit = 2
      else if (num_units .eq. 8) then
          recl_unit = 1
      else
          stop 'Impossible "recl=" specification'
      endif

      return
      end



