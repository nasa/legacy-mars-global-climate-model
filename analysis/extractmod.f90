      module extractmod

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
      real*4 SIND, decmax, eccn, orbinc, vinc, gasp2
      real*4 cp, stbo, xlhtc, kapa, cmk
      real*4 alicen, alices, egoco2n, egoco2s
      real*4 npcwikg, gidn, gids

      character(len=7) version

!  Surface data sets

      real*4, allocatable :: topog(:,:),alsp(:,:),zin(:,:,:)
      logical, allocatable :: npcflag(:,:)
      real*4, allocatable :: dsig(:), dxyp(:)

      real*4, allocatable :: STRESSX(:,:), STRESSY(:,:)
      real*4, allocatable :: P(:,:), T(:,:,:), U(:,:,:), V(:,:,:)
      real*4, allocatable :: TSTRAT(:,:), GT(:,:), CO2ICE(:,:)
      real*4, allocatable :: QTRACE(:,:,:,:)
      real*4, allocatable :: DHEAT(:,:,:)
      real*4, allocatable :: GEOP(:,:,:)
      real*4, allocatable :: QCOND(:,:,:)
      real*4, allocatable :: surfalb(:,:)
      real*4, allocatable :: srfdnflx(:,:,:)
      real*4, allocatable :: srfupflx(:,:,:)
      real*4, allocatable :: srfupflx_dd(:,:,:)
      real*4, allocatable :: tau3d(:,:,:)
      real*4, allocatable :: taucld3d(:,:,:)
      real*4, allocatable :: SSUN(:,:),TAUSURF(:,:),tautotji(:,:)
      real*4, allocatable :: STEMP(:,:,:), SDEPTH(:)

!  Radiation output

      real*4, allocatable :: fuptopv(:,:), fdntopv(:,:)
      real*4, allocatable :: fupsurfv(:,:), fdnsurfv(:,:)
      real*4, allocatable :: fuptopir(:,:)
      real*4, allocatable :: fupsurfir(:,:), fdnsurfir(:,:)

!  Calculated Arrays

      real*4, allocatable :: lat(:), lon(:), lsseas(:), nlsbin(:)
      real*4, allocatable :: gasp(:), gagt(:), gd(:), gwv(:), plx(:)
      real*4, allocatable :: gwi(:), hco2_n(:), hco2_s(:), lssol(:)
      real*4, allocatable :: vl1(:), vl2(:)
      real*4, allocatable :: zd(:,:), zwv(:,:), zwi(:,:)
      real*4, allocatable :: zws(:,:), zss(:,:), zt(:,:)
      real*4, allocatable :: zt2pm(:,:), zt2am(:,:), zco2i(:,:)
      real*4, allocatable :: znumt2pm(:,:), znumt2am(:,:)
      real*4, allocatable :: z0pt5t(:,:), z0pt5t2pm(:,:), z0pt5t2am(:,:)
      real*4, allocatable :: uwnd(:,:,:), vwnd(:,:,:), ws(:,:,:)
      real*4, allocatable :: ustr(:,:,:), vstr(:,:,:), mstr(:,:,:)
      real*4, allocatable :: co2i(:,:,:), tsrf(:,:,:), tsrf2pm(:,:,:)
      real*4, allocatable :: tsrf2am(:,:,:), t0pt5(:,:,:)
      real*4, allocatable :: numtsrf2pm(:,:,:),numtsrf2am(:,:,:)
      real*4, allocatable :: t0pt52pm(:,:,:), t0pt52am(:,:,:)
      real*4, allocatable :: zpu(:,:,:), zpv(:,:,:), zpt(:,:,:)
      real*4, allocatable :: zpt2pm(:,:,:), zpt2am(:,:,:), zpmsf(:,:,:)
      real*4, allocatable :: numt(:,:,:),numt2pm(:,:,:),numt2am(:,:,:)
      real*4, allocatable :: numw(:,:,:), sigma(:),zpsrf(:,:)


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
               SIND, GASP2

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
      read(20) DHEAT
      read(20) GEOP
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
        read(20) !DHEAT
        read(20) !GEOP
   40 continue
      return
      end subroutine skiprecords
      subroutine myalloc1(NI,NJ,NLAY,NL,NT)
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

      allocate(STEMP(NJ,NI,NL),stat=ierror)
      if(ierror.ne.0) call errorcode("STEMP")

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

      allocate(DHEAT(NJ,NI,NLAY),stat=ierror)
      if(ierror.ne.0) call errorcode("DHEAT")

      allocate(GEOP(NJ,NI,NLAY),stat=ierror)
      if(ierror.ne.0) call errorcode("GEOP")

      allocate(QCOND(NJ,NI,NT),stat=ierror)
      if(ierror.ne.0) call errorcode("QCOND")

      allocate(srfdnflx(NJ,NI,NT),stat=ierror)
      if(ierror.ne.0) call errorcode("srfdnflx")

      allocate(srfupflx(NJ,NI,NT),stat=ierror)
      if(ierror.ne.0) call errorcode("srfupflx")

      allocate(srfupflx_dd(NJ,NI,NT),stat=ierror)
      if(ierror.ne.0) call errorcode("srfupflx_dd")

      return
      end subroutine myalloc1

      subroutine myalloc2(NI,NJ,NLAY,NL,NT,NS,ND,NP)
      integer, intent(IN) :: NI, NJ, NLAY, NL, NT, NS, ND, NP

      allocate(tau3d(NJ,NI,2*NLAY+3),stat=ierror)
      if(ierror.ne.0) call errorcode("tau3d")

      allocate(taucld3d(NJ,NI,2*NLAY+3),stat=ierror)
      if(ierror.ne.0) call errorcode("taucld3d")

      allocate(lat(nj-1),stat=ierror)
      if(ierror.ne.0) call errorcode("lat")

      allocate(lon(ni),stat=ierror)
      if(ierror.ne.0) call errorcode("lon")

      allocate(plx(nlay),stat=ierror)
      if(ierror.ne.0) call errorcode("plx")

      allocate(sigma(2*nlay+3),stat=ierror)
      if(ierror.ne.0) call errorcode("sigma")
 
      allocate(lsseas(ns),stat=ierror)
      if(ierror.ne.0) call errorcode("lsseas")

      allocate(nlsbin(ns),stat=ierror)
      if(ierror.ne.0) call errorcode("nlsbin")

      allocate(lssol(nd),stat=ierror)
      if(ierror.ne.0) call errorcode("lssol")

      allocate(gasp(nd),stat=ierror)
      if(ierror.ne.0) call errorcode("gasp")

      allocate(gagt(nd),stat=ierror)
      if(ierror.ne.0) call errorcode("gagt")

      allocate(gd(nd),stat=ierror)
      if(ierror.ne.0) call errorcode("gd")

      allocate(gwv(nd),stat=ierror)
      if(ierror.ne.0) call errorcode("gwv")

      allocate(gwi(nd),stat=ierror)
      if(ierror.ne.0) call errorcode("gwi")

      allocate(hco2_n(nd),stat=ierror)
      if(ierror.ne.0) call errorcode("hco2_n")

      allocate(hco2_s(nd),stat=ierror)
      if(ierror.ne.0) call errorcode("hco2_s")

      allocate(vl1(nd),stat=ierror)
      if(ierror.ne.0) call errorcode("vl1")

      allocate(vl2(nd),stat=ierror)
      if(ierror.ne.0) call errorcode("vl2")

      allocate(zd(nd,nj-1),stat=ierror)
      if(ierror.ne.0) call errorcode("zd")

      allocate(zwv(nd,nj-1),stat=ierror)
      if(ierror.ne.0) call errorcode("znumt2amwv")

      allocate(zwi(nd,nj-1),stat=ierror)
      if(ierror.ne.0) call errorcode("zwi")

      allocate(zws(nd,nj-1),stat=ierror)
      if(ierror.ne.0) call errorcode("zws")

      allocate(zss(nd,nj-1),stat=ierror)
      if(ierror.ne.0) call errorcode("zss")

      allocate(zt(nd,nj-1),stat=ierror)
      if(ierror.ne.0) call errorcode("zt")

      allocate(zt2pm(nd,nj-1),stat=ierror)
      if(ierror.ne.0) call errorcode("zt2pm")

      allocate(zt2am(nd,nj-1),stat=ierror)
      if(ierror.ne.0) call errorcode("zt2am")

      allocate(zco2i(nd,nj-1),stat=ierror)
      if(ierror.ne.0) call errorcode("zco2i")

      allocate(z0pt5t(nd,nj-1),stat=ierror)
      if(ierror.ne.0) call errorcode("z0pt5t")

      allocate(z0pt5t2pm(nd,nj-1),stat=ierror)
      if(ierror.ne.0) call errorcode("z0pt5t2pm")

      allocate(z0pt5t2am(nd,nj-1),stat=ierror)
      if(ierror.ne.0) call errorcode("z0pt5t2am")

      allocate(znumt2pm(nd,nj-1),stat=ierror)
      if(ierror.ne.0) call errorcode("znumt2pm")

      allocate(znumt2am(nd,nj-1),stat=ierror)
      if(ierror.ne.0) call errorcode("znumt2am")

      allocate(uwnd(ni,nj-1,ns),stat=ierror)
      if(ierror.ne.0) call errorcode("uwnd")

      allocate(vwnd(ni,nj-1,ns),stat=ierror)
      if(ierror.ne.0) call errorcode("vwnd")

      allocate(ws(ni,nj-1,ns),stat=ierror)
      if(ierror.ne.0) call errorcode("ws")

      allocate(ustr(ni,nj-1,ns),stat=ierror)
      if(ierror.ne.0) call errorcode("ustr")

      allocate(vstr(ni,nj-1,ns),stat=ierror)
      if(ierror.ne.0) call errorcode("vstr")

      allocate(mstr(ni,nj-1,ns),stat=ierror)
      if(ierror.ne.0) call errorcode("mstr")

      allocate(co2i(ni,nj-1,ns),stat=ierror)
      if(ierror.ne.0) call errorcode("co2i")

      allocate(tsrf(ni,nj-1,ns),stat=ierror)
      if(ierror.ne.0) call errorcode("tsrf")

      allocate(tsrf2pm(ni,nj-1,ns),stat=ierror)
      if(ierror.ne.0) call errorcode("tsrf2pm")

      allocate(tsrf2am(ni,nj-1,ns),stat=ierror)
      if(ierror.ne.0) call errorcode("tsrf2am")

      allocate(t0pt5(ni,nj-1,ns),stat=ierror)
      if(ierror.ne.0) call errorcode("t0pt5")

      allocate(t0pt52pm(ni,nj-1,ns),stat=ierror)
      if(ierror.ne.0) call errorcode("t0pt52pm")

      allocate(t0pt52am(ni,nj-1,ns),stat=ierror)
      if(ierror.ne.0) call errorcode("t0pt52am")

      allocate(zpu(nj-1,np,ns),stat=ierror)
      if(ierror.ne.0) call errorcode("zpu")

      allocate(zpv(nj-1,np,ns),stat=ierror)
      if(ierror.ne.0) call errorcode("zpv")

      allocate(zpt(nj-1,np,ns),stat=ierror)
      if(ierror.ne.0) call errorcode("zpt")

      allocate(zpsrf(nj-1,ns),stat=ierror)
      if(ierror.ne.0) call errorcode("zpsrf")

      allocate(zpt2pm(nj-1,np,ns),stat=ierror)
      if(ierror.ne.0) call errorcode("zpt2pm")

      allocate(zpt2am(nj-1,np,ns),stat=ierror)
      if(ierror.ne.0) call errorcode("zpt2am")

      allocate(zpmsf(nj-1,np,ns),stat=ierror)
      if(ierror.ne.0) call errorcode("zpmsf")

      allocate(numt(nj-1,np,ns),stat=ierror)
      if(ierror.ne.0) call errorcode("numt")

      allocate(numt2pm(nj-1,np,ns),stat=ierror)
      if(ierror.ne.0) call errorcode("numt2pm")

      allocate(numt2am(nj-1,np,ns),stat=ierror)
      if(ierror.ne.0) call errorcode("numt2am")

      allocate(numw(nj-1,np,ns),stat=ierror)
      if(ierror.ne.0) call errorcode("numw")

      allocate(numtsrf2pm(ni,nj-1,ns),stat=ierror)
      if(ierror.ne.0) call errorcode("numtsrf2pm")

      allocate(numtsrf2am(ni,nj-1,ns),stat=ierror)
      if(ierror.ne.0) call errorcode("numtsrf2am")

          lssol(:)=0.
          gasp(:)=0.
          gagt(:)=0.
          gd(:)=0.
          gwv(:)=0.
          gwi(:)=0.
          hco2_n(:)=0.
          hco2_s(:)=0.
          vl1(:)=0.
          vl2(:)=0.

          zd(:,:)=0.
          zwv(:,:)=0.
          zwi(:,:)=0.
          zws(:,:)=0.
          zss(:,:)=0.
          zt(:,:)=0.
          zt2pm(:,:)=0.
          zt2am(:,:)=0.
          zco2i(:,:)=0.
          z0pt5t(:,:)=0.
          z0pt5t2pm(:,:)=0.
          z0pt5t2am(:,:)=0.
          znumt2pm(:,:)=0.
          znumt2am(:,:)=0.

          uwnd(:,:,:)=0.
          vwnd(:,:,:)=0.
          ws(:,:,:)=0.
          ustr(:,:,:)=0.
          vstr(:,:,:)=0.
          mstr(:,:,:)=0.
          co2i(:,:,:)=0.
          tsrf(:,:,:)=0.
          tsrf2pm(:,:,:)=0.
          tsrf2am(:,:,:)=0.
          t0pt5(:,:,:)=0.
          t0pt52pm(:,:,:)=0.
          t0pt52am(:,:,:)=0.
          numtsrf2pm(:,:,:)=0.
          numtsrf2am(:,:,:)=0.
          zpu(:,:,:)=0.
          zpv(:,:,:)=0.
          zpt(:,:,:)=0.
          zpsrf(:,:)=0.
          zpt2pm(:,:,:)=0.
          zpt2am(:,:,:)=0.
          zpmsf(:,:,:)=0.
          numt(:,:,:)=0.
          numt2pm(:,:,:)=0.
          numt2am(:,:,:)=0.
          numw(:,:,:)=0.


      return
      end subroutine myalloc2

      subroutine errorcode(name)
      character(len=*) name
      write(6,'("Error allocating array ",(A))') name
      stop
      return
      end subroutine errorcode
      end module extractmod
