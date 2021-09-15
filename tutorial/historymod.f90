      module historymod

!  Cloud code parameters

!  Set the various tracer index

!  For dust: Mass and Number
      integer, parameter :: iMa_dt  = 1
      integer, parameter :: iNb_dt  = 2

!  For water ice: Mass and Number
      integer, parameter :: iMa_cld = 3
      integer, parameter :: iNb_cld = 4

!  For dust core: only Mass
      integer, parameter :: iMa_cor = 5

!  For water vapor: only Mass
      integer, parameter :: iMa_vap = 6

!  History scalars

      real*4 runnum, grav, rgas

      INTEGER KTP4, NC3, NCYCLE, JM, IM, LAYERS, NL, NTRACE

      real*4 TAU, VPOUT, RSDIST, TOFDAY, PSF, PTROP, TAUTOT, RPTAU
      real*4 SIND, decmax, eccn, orbinc, vinc, gasp
      real*4 cp, stbo, xlhtc, kapa, cmk
      real*4 alicen, alices, egoco2n, egoco2s
      real*4 npcwikg, gidn, gids

      character(len=7) version

!  Surface data sets

      real*4, allocatable :: topog(:,:), alsp(:,:), zin(:,:,:)
      logical, allocatable :: npcflag(:,:)
      real*4, allocatable :: dsig(:), dxyp(:)

      real*4, allocatable :: STRESSX(:,:), STRESSY(:,:)
      real*4, allocatable :: P(:,:), T(:,:,:), U(:,:,:), V(:,:,:)
      real*4, allocatable :: TSTRAT(:,:), GT(:,:), CO2ICE(:,:)
      real*4, allocatable :: QTRACE(:,:,:,:)
      real*4, allocatable :: QCOND(:,:,:)
      real*4, allocatable :: surfalb(:,:)
      real*4, allocatable :: srfdnflx(:,:,:)
      real*4, allocatable :: srfupflx(:,:,:)
      real*4, allocatable :: srfupflx_dd(:,:,:)
      real*4, allocatable :: tau3d(:,:,:)
      real*4, allocatable :: taucld3d(:,:,:)
      real*4, allocatable :: SSUN(:,:),TAUSURF(:,:),tautotji(:,:)
      real*4, allocatable :: STEMP(:,:,:), SDEPTH(:)
      real*4, allocatable :: DHEAT(:,:,:),GEOP(:,:,:)

!  Radiation output

      real*4, allocatable :: fuptopv(:,:), fdntopv(:,:)
      real*4, allocatable :: fupsurfv(:,:), fdnsurfv(:,:)
      real*4, allocatable :: fuptopir(:,:)
      real*4, allocatable :: fupsurfir(:,:), fdnsurfir(:,:)

      contains
      subroutine readheader
      implicit none

      read(20) DSIG, DXYP, GRAV, RGAS, cp, stbo, xlhtc, kapa, cmk, &
               decmax, eccn, orbinc, vinc, sdepth, alicen, alices, &
               egoco2n, egoco2s, npcwikg, gidn, gids

      read(20) TOPOG, ALSP, ZIN, NPCFLAG
      return
      end subroutine readheader
      subroutine readvariables
      implicit none
      read(20) TAU, VPOUT, RSDIST, TOFDAY, PSF, PTROP, TAUTOT, RPTAU, &
               SIND, GASP

      read(20) NC3, NCYCLE

      read(20) P
      read(20) T
      read(20) U
      read(20) V
      read(20) GT
      read(20) CO2ICE
      read(20) STRESSX
      read(20) STRESSY
      read(20) TSTRAT
      read(20) TAUSURF
      read(20) SSUN
      read(20) QTRACE
      read(20) QCOND
      read(20) STEMP
      read(20) fuptopv, fdntopv, fupsurfv, fdnsurfv
      read(20) fuptopir, fupsurfir, fdnsurfir
      read(20) surfalb
      read(20) dheat
      read(20) geop

      return
      end subroutine readvariables
      subroutine skiprecords(nskip)
      implicit none
      integer :: nskip, i, ios, nr
      nr = 0
      write(6,'(" ")')
      do 40 i=1,nskip
        read(20,IOSTAT=ios) !
        read(20,IOSTAT=ios) !
        if(ios.ne.0) then
          write(6,'("File shorter than expected.")')
          write(6,'("Last record: ",i3)') nr
          stop
        end if
        nr = nr + 1
        read(20) !P
        read(20) !T
        read(20) !U
        read(20) !V
        read(20) !GT
        read(20) !CO2ICE
        read(20) !STRESSX
        read(20) !STRESSY
        read(20) !TSTRAT
        read(20) !TAUSURF
        read(20) !SSUN
        read(20) !QTRACE
        read(20) !QCOND
        read(20) !STEMP
        read(20) !fuptopv, fdntopv, fupsurfv, fdnsurfv
        read(20) !fuptopir, fupsurfir, fdnsurfir
        read(20) !surfalb
        read(20) !dheat
        read(20) !geop
   40 continue
      return
      end subroutine skiprecords
      subroutine myalloc(NI,NJ,NLAY,NL,NT)
      integer, intent(IN) :: NI, NJ, NLAY, NL, NT

      allocate(dsig(NLAY),stat=ierror)
      if(ierror.ne.0) call errorcode("dsig")

      allocate(dxyp(NJ),stat=ierror)
      if(ierror.ne.0) call errorcode("dxyp")

      allocate(topog(NJ,NI),stat=ierror)
      if(ierror.ne.0) call errorcode("topog")

      allocate(alsp(NJ,NI),stat=ierror)
      if(ierror.ne.0) call errorcode("alsp")

      allocate(zin(NJ,NI,NL),stat=ierror)
      if(ierror.ne.0) call errorcode("zin")

      allocate(npcflag(NJ,NI),stat=ierror)
      if(ierror.ne.0) call errorcode("npcflag")

      allocate(P(NJ,NI),stat=ierror)
      if(ierror.ne.0) call errorcode("P")

      allocate(T(NJ,NI,NLAY),stat=ierror)
      if(ierror.ne.0) call errorcode("T")

      allocate(U(NJ,NI,NLAY),stat=ierror)
      if(ierror.ne.0) call errorcode("U")

      allocate(V(NJ,NI,NLAY),stat=ierror)
      if(ierror.ne.0) call errorcode("V")

      allocate(surfalb(NJ,NI),stat=ierror)
      if(ierror.ne.0) call errorcode("surfalb")

      allocate(STRESSX(NJ,NI),stat=ierror)
      if(ierror.ne.0) call errorcode("STRESSX")

      allocate(STRESSY(NJ,NI),stat=ierror)
      if(ierror.ne.0) call errorcode("STRESSY")

      allocate(TSTRAT(NJ,NI),stat=ierror)
      if(ierror.ne.0) call errorcode("TSTRAT")

      allocate(GT(NJ,NI),stat=ierror)
      if(ierror.ne.0) call errorcode("GT")

      allocate(CO2ICE(NJ,NI),stat=ierror)
      if(ierror.ne.0) call errorcode("CO2ICE")

      allocate(SSUN(NJ,NI),stat=ierror)
      if(ierror.ne.0) call errorcode("SSUN")

      allocate(TAUSURF(NJ,NI),stat=ierror)
      if(ierror.ne.0) call errorcode("TAUSURF")

      allocate(TAUTOTJI(NJ,NI),stat=ierror)
      if(ierror.ne.0) call errorcode("TAUTOTJI")

      allocate(QTRACE(NJ,NI,NLAY,NT),stat=ierror)
      if(ierror.ne.0) call errorcode("QTRACE")

      allocate(QCOND(NJ,NI,NT),stat=ierror)
      if(ierror.ne.0) call errorcode("QCOND")

      allocate(srfdnflx(NJ,NI,NT),stat=ierror)
      if(ierror.ne.0) call errorcode("srfdnflx")

      allocate(srfupflx(NJ,NI,NT),stat=ierror)
      if(ierror.ne.0) call errorcode("srfupflx")

      allocate(srfupflx_dd(NJ,NI,NT),stat=ierror)
      if(ierror.ne.0) call errorcode("srfupflx_dd")

      allocate(tau3d(NJ,NI,2*NLAY+3),stat=ierror)
      if(ierror.ne.0) call errorcode("tau3d")

      allocate(taucld3d(NJ,NI,2*NLAY+3),stat=ierror)
      if(ierror.ne.0) call errorcode("taucld3d")

      allocate(STEMP(NJ,NI,NL),stat=ierror)
      if(ierror.ne.0) call errorcode("STEMP")

      allocate(DHEAT(NJ,NI,NLAY),stat=ierror)
      if(ierror.ne.0) call errorcode("DHEAT")

      allocate(GEOP(NJ,NI,NLAY),stat=ierror)
      if(ierror.ne.0) call errorcode("GEOP")

      allocate(SDEPTH(2*NL+1),stat=ierror)
      if(ierror.ne.0) call errorcode("SDEPTH")

      allocate(FUPTOPV(NJ,NI),stat=ierror)
      if(ierror.ne.0) call errorcode("FUPTOPV")

      allocate(FDNTOPV(NJ,NI),stat=ierror)
      if(ierror.ne.0) call errorcode("FDNTOPV")

      allocate(FUPSURFV(NJ,NI),stat=ierror)
      if(ierror.ne.0) call errorcode("FUPSURFV")

      allocate(FDNSURFV(NJ,NI),stat=ierror)
      if(ierror.ne.0) call errorcode("FDNSURFV")

      allocate(FUPTOPIR(NJ,NI),stat=ierror)
      if(ierror.ne.0) call errorcode("FUPTOPIR")

      allocate(FUPSURFIR(NJ,NI),stat=ierror)
      if(ierror.ne.0) call errorcode("FUPSURFIR")

      allocate(FDNSURFIR(NJ,NI),stat=ierror)
      if(ierror.ne.0) call errorcode("FDNSURFIR")

      return
      end subroutine myalloc

      subroutine errorcode(name)
      character(len=*) name
      write(6,'("Error allocating array ",(A))') name
      stop
      return
      end subroutine errorcode
      end module historymod
