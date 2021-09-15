      PROGRAM htest_sa

!     GCM history tape analysis program.  Modelled after HTEST.f
!     used on the Crays to access the history tapes and report
!     the results foar ANY grid point/layer combination for
!     ANY HIST or HIST2D variable.
!
!
!     Version 0.1
!     Started April 23, 1992
!
!
!     compile:  f90 -o htest htest.f90 historymod.o
!
      use historymod

!----------------------------------------------------------------------#

      CHARACTER(len=72) fname
      real*4 tlocal
      integer L_J, L_I, L_LAYERS
      integer :: nr

!----------------------------------------------------------------------#

!     Get history file name and store it in fname.
 
      write(6,'("History file name:  ")',advance='no')
      read(5,'(A)') fname

!     ok, try and OPEN the history file, stop if can't open file.

      open(unit=20,file=fname,status='old',form='unformatted',   &
           IOSTAT=ios)

      if(ios.ne.0) then
        write(6,'("Could not open data file.")')
        stop
      end if

!     input parameters for wanted history record.

      write(6,'("Record number?  ")',advance='no')
      read(5,*) nrec
      write(6,'("J, I, L (Which are: Lat, Lon, Layer)  ")',advance='no')
      read(5,*) js1, is1, ls1

!     Make a few simple checks to make sure that pure garbage is
!     not asked for.

      if(js1.lt.1.or.is1.lt.1.or.ls1.lt.1) then
        js1  = 1
        is1  = 1
        ls1  = 1
      end if

      if(nrec.lt.1) nrec=1
      nskip = nrec-1

!     Read header information.  These variables are written only once
!     at the beginning of each history tape.

      read(20) RUNNUM, JM, IM, LAYERS, NL, ntrace

      L_J      = JM
      L_I      = IM
      L_LAYERS = LAYERS

!  Allocate history arrays

      call myalloc(IM,JM,LAYERS,NL,NTRACE)
      call readheader

      tp2 = -topog(JS1,IS1)/3720.0
      write(6,'(" ")')
      write(6,34) runnum
   34 format(' Run number: ',f7.2)

!     Now read the part of the history tape that is written each time
!     HISTRY is called (currently once every 1.5 simulated hours).

      call skiprecords(nskip)
      call readvariables
      call localtime(IM,tofday,IS1,TLOCAL)
     
      write(6,'(" ")')
      write(6,'(" ")')
      write(6,'(3x,"History file name:  ",A)') trim(fname)
      write(6,'(10x,"Run number:  ",f7.2)') runnum
      write(6,'(7x,"Record number:    ",i3)') nrec
      write(6,                                                         &
          '(10x,"      Grid:  J = ",i2,"    I = ",i2,"    L = ",i2)')  &
                 JS1, IS1, LS1
      write(6,'(" ")')
      write(6,'(10x,"       Ls = ",f11.2)') VPOUT
      write(6,'(10x,"   RSDIST = ",f11.4)') RSDIST
      write(6,'(10x,"   DECMAX = ",f11.4)') DECMAX
      write(6,'(10x,"      TAU = ",f11.2)') TAU
      write(6,'(10x,"   TOFDAY = ",f11.2)') TOFDAY
      write(6,'(" Time at Grid Point = ",f11.2)') TLOCAL
      write(6,'(" ")')
      write(6,'(10x,"   TAUTOT = ",F11.4)') TAUTOT
!     write(6,'(10x,"   CONRNU = ",1pe11.4)') CONRNU
      write(6,'(10x,"    RPTAU = ",F11.2)') RPTAU
      write(6,'(" ")')
      write(6,'(9x,"TOPOG(J,I) = ",f11.4,"  -----> ",f8.4," km")')     &
                    topog(js1,is1), tp2
      write(6,'(10x,"ALSP(J,I) = ",f11.4)') alsp(js1,is1)
      write(6,'(7x,"SURFALB(J,I) = ",f11.4)') surfalb(js1,is1)
      write(6,'(9x,"ZIN(J,I,1) = ",f11.4)') zin(js1,is1,1)
      write(6,'(10x,"     GIDN = ",f11.4)',advance='no') gidn
      write(6,'(7x,"GIDS = ",f11.4)') gids
      write(6,'(" ")')

      write(6,'(10x,"      PSF = ",F11.4)') PSF
      write(6,'(10x,"     GASP = ",F11.4)') GASP
      write(6,'(" GASP: Global Average Surface Pressure")')
      write(6,'(" ")')

!     Print individual elements of the HIST array

      write(6,'(10x,"    PTROP = ",1pe11.4)') PTROP
      write(6,'(10x,"   P(J,I) = ",F11.4)') P(JS1,IS1)
      write(6,'(8x,"TSTRAT(J,I) = ",F11.4)') TSTRAT(JS1,IS1)
      write(6,'(10x," T(J,I,L) = ",F11.4)') T(JS1,IS1,LS1)
      write(6,'(10x," U(J,I,L) = ",F11.4)') U(JS1,IS1,LS1)
      write(6,'(10x," V(J,I,L) = ",F11.4)') V(JS1,IS1,LS1)
      write(6,'(" ")')
      write(6,'(10x,"  GT(J,I) = ",f11.4)') GT(JS1,IS1)
      write(6,'(7x,"STEMP(J,I,1) = ",f11.4)',advance='no')       &
                    STEMP(JS1,IS1,1)
      write(6,'(8x,"SDEPTH( 2) = ",f9.4," m")') SDEPTH(2)
      write(6,'(7x,"STEMP(J,I,5) = ",f11.4)',advance='no')       &
                    STEMP(JS1,IS1,5)
      write(6,'(8x,"SDEPTH(10) = ",f9.4," m")') SDEPTH(10)
      write(6,'(" ")')

      write(6,'(8x,"CO2ICE(J,I) = ",1pe11.4)') CO2ICE(JS1,IS1)
      write(6,'(8x,"ALICEN      = ",f7.4,12x,"ALICES       = ",f7.4)') &
                    ALICEN, ALICES
      write(6,'(8x,"EGOCO2N     = ",f7.4,12x,"EGOCO2S      = ",f7.4)') &
                    EGOCO2N, EGOCO2S
      write(6,'(7x,"STRESSX(J,I) = ",1pe11.4)',advance='no')      &
                    STRESSX(JS1,IS1)
      write(6,'(8x,"STRESSY(J,I) = ",1pe11.4)') STRESSY(JS1,IS1)
!     write(6,'(10x,"SSUN(J,I) = ",f11.5)') SSUN(JS1,IS1)
      write(6,'(7x,"TAUSURF(J,I) = ",f11.5)') tausurf(JS1,IS1)

      write(6,'(" ")')
      write(6,'(7x,"fuptopv(J,I) = ",f11.5)',advance='no')  &
                     fuptopv(JS1,IS1)
      write(6,'(8x,"fuptopir(J,I)  = ",f11.5)') fuptopir(JS1,IS1)
      write(6,'(7x,"fdntopv(J,I) = ",f11.5)') fdntopv(JS1,IS1)
      write(6,'(6x,"fupsurfv(J,I) = ",f11.5)',advance='no') &
                     fupsurfv(JS1,IS1)
      write(6,'(8x,"fupsurfir(J,I) = ",f11.5)') fupsurfir(JS1,IS1)
      write(6,'(6x,"fdnsurfv(J,I) = ",f11.5)',advance='no') &
                     fdnsurfv(JS1,IS1)
      write(6,'(8x,"fdnsurfir(J,I) = ",f11.5)') fdnsurfir(JS1,IS1)

      write(6,'(" ")')
      write(6,'(8x,"NPCFLAG =  ",L1)') NPCFLAG(js1,is1)
      write(6,'(4x,"Water vapor = ",1pe11.4)')               &
                             qtrace(js1,is1,ls1,iMa_vap)
      end
      subroutine localtime(IM,tofday,Igrid,tlocal)
      implicit none
      integer, parameter :: rk4 = selected_real_kind(6)
      real(rk4) tofday, tlocal, dtime
      integer  igrid, IM, IMo2

      !================================================!

      IMO2  = IM/2 + 1
      DTIME = 24.0/float(IM)
 
      if(Igrid.eq.IMO2) then
        tlocal = tofday
      elseif(Igrid.lt.IMO2) then
        tlocal = tofday-DTIME*float(IMO2-Igrid)
        if(tlocal.lt.0.0) tlocal = tlocal+24.0
      else
        tlocal = tofday+DTIME*float(Igrid-IMO2)
        if(tlocal.gt.24.0) tlocal = tlocal-24.0
      end if

      return

      end
