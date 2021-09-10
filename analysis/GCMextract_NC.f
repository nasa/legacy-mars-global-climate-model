      PROGRAM GCMextract_NC
!  This is an analysis package designed to compute basic fields from
!  the GCM.  Files written out from this routine can be read in and
!  plotted in IDL, Matlab, FORTRAN, etc.
!
!  M. Kahre
!
!  Array size definitions:
!
!	nsoltot: 		number of simulated days to be analyzed
!	nlat (jm-1): 		number of latitude grid points
!	nlon (im): 		number of longitude grid points
!	npres (np_interp): 	number of pressure levels in the vertical
!	nseas:  		number of seasons (taken to be 12)
!
!  Computed Fields:
!
!	1.  gasp(nsoltot): sol mean global mean surface pressure
!	2.  gagt(nsoltot): sol mean global mean surface temperature
!	3.  gd(nsoltot):  sol mean globally integrated dust (kg)
!	4.  gwv(nsoltot): sol mean globally integrated water vapor (kg)
!       5.  gwi(nsoltot): sol mean globally integrated water ice (kg)
!       6.  hco2_n(nsoltot): sol mean n-hemispherically integrated surface CO2 ice (kg)
!	7.  hco2_s(nsoltot): sol mean s-hemispherically integrated surface CO2 ice (kg)
!	8.  vl1(nsoltot): sol mean hydrostatically adjusted surface pressure at VL1 (mb)
!	9.  vl2(nsoltot): sol mean hydrostatically adjusted surface pressure at VL2 (mb)
!	
!	10. zd(nsoltot,nlat): sol mean zonal mean column integrated dust opacity at 9 um
!       11. zwv(nsoltot,nlat): sol mean zonal mean column integrated water vapor (pr-um)
!	12. zwi(nsoltot,nlat): sol mean zonal mean column integrated water ice (pr-um)
!	13. zws(nsoltot,nlat): sol mean zonal mean surface wind speed (m/s)
!	14. zss(nsoltot,nlat): sol mean zonal mean surface stress magnitude (mN/m^2)
!	15. zt(nsoltot,nlat): sol mean zonal mean surface temperature (K)
!	16. zt2pm(nsoltot,nlat): 1-3 pm mean zonal mean surface temperature (K)
!  	17. zt2am(nsoltot,nlat): 1-3 am mean zonal mean surface temperature (K)
!	18. zco2i(nsoltot,nlat): sol mean zonal mean surface CO2 ice (kg/m^2)
!	19. z0pt5t(nsoltot,nlat): sol mean zonal mean 0.5 mb temperature (K)
!	20. z0pt5t2pm(nsoltot,nlat): 1-3 pm mean zonal mean 0.5 mb temperature (K)
!	21. z0pt5t2am(nsoltot,nlat): 1-3 am mean zonal mean 0.5 mb temperature (K)
!
!	22. zpu(nlat,npres,nseas): 30 Ls mean zonal mean cross-section of U wind
!	23. zpv(nlat,npres,nseas): 30 Ls mean zonal mean cross-section of V wind
! 	24. zpt(nlat,npres,nseas): 30 Ls mean zonal mean cross-section of temperature (K)
!	25. zpt2pm(nlat,npres,nseas): 30 Ls 1-3 pm mean zonal mean cross-section of temperature (K)
! 	26. zpt2am(nlat,npres,nseas): 30 Ls 1-3 am mean zonal mean cross-section of temperature (K)
!	27. zpmsf(nlat,npres,nseas): 30 Ls mean mass stream function (10^8 kg/s)
!	
!	28. uwnd(nlon,nlat,nseas): 30 Ls mean U wind at surface (m/s)
!	29. vwnd(nlon,nlat,nseas): 30 Ls mean V wind at surface (m/s)
!	30. ws(nlon,nlat,nseas): 30 Ls mean wind speed at surface (m/s)
!	31. ustr(nlon,nlat,nseas): 30 Ls mean zonal surface stress (mN/m^2)
!	32. vstr(nlon,nlat,nseas): 30 Ls mean meridional surface stress (mN/m^2)
!	33. mstr(nlon,nlat,nseas): 30 Ls mean surface stress magnitude (mN/m^2)
!	34. co2i(nlon,nlat,nseas): 30 Ls mean surface CO2 ice budget (kg/m^2)
!	35. tsrf(nlon,nlat,nseas): 30 Ls mean surface temperature (K)
!	36. tsrf2pm(nlon,nlat,nseas): 30 Ls 1-3 pm mean surface temperature (K)
! 	37. tsrf2am(nlon,nlat,nseas): 30 Ls 1-3 am mean surface temperature (K)
!	38: t0pt5(nlon,nlat,nseas): 30 Ls mean 0.5 mb temperature (K)
!	39: t0pt52pm(nlon,nlat,nseas): 30 Ls 1-3 pm mean 0.5 mb temperature (K)
! 	40: t0pt52am(nlon,nlat,nseas): 30 Ls 1-3 am mean 0.5 mb temperature (K)
!
!       41: zpsrf(nlat,nseas): 30 Ls mean zonal mean surface pressure (mb)
!
!----------------------------------------------------------------------!

      use extractmod
      implicit none
      include 'netcdf.inc'
 
      integer, parameter :: np_interp = 41
      integer :: nsol  = 10
      integer :: nseas = 12
      integer :: wsol  = 16
      integer :: nsolcnt = 0
      integer :: nsoltot

      integer :: ifile1, ifile2, m, lm, k, lp, jj, ii, lsbin,nyr
      real*4  :: pi, wvcol,wicol, wu, wv, ss, wuv, tlocal, psx
      real*4  :: uavgl, uavglp1, vavgl, vavglp1, cosy, fl, flp1
      real*4  :: pavesfc, radius, zl, zlp1, sum
      character(len=72) :: fname
      integer :: ll,l1,ii1,ii2,iu1,iu2
      real*4  :: area = 0.0


! Surface Pressure Stuff

      real*4  :: psvl1(4),psvl2(4),gtvl1(4),gtvl2(4)
      real*4  :: psfgm
      integer :: j1(4),i1(4),j2(4),i2(4)
      integer :: nn,nk
      real*4  :: vl1z(4),vl2z(4)
      real*4  :: vlp(3),fact1,fact2,fact3,term1,term2,term3
      real*4  :: vlt(3),h1,h2,d1,d2
      real*4  :: vl1gpp(4),vl1gpt(4)
      real*4  :: vl2gpp(4),vl2gpt(4)
      real*4  :: rad

      integer :: num = 0

      real*4,  allocatable :: latitude(:), longitude(:)

! zonal tums stuff

      real*4  :: pinterp(np_interp)
      real*4  :: zoverh(np_interp)

      real*4  :: mls, tot

      integer :: i, j, l, n, nd, ierror,nf,ios
      real*4  :: als, ati

      data vl1z/-3.155,-3.658,-3.443,-3.797/
      data vl2z/-4.559,-4.280,-4.500,-4.281/

      data pinterp / 0.00056, 0.00075, 0.001, 0.0013, 0.0018,
     *               0.0024, 0.0032, 0.0042, 0.0056, 0.0075, 0.010,
     *               0.013, 0.018, 0.024, 0.032, 0.042, 0.056, 0.075,
     *               0.100, 0.133, 0.178, 0.237, 0.316, 0.422, 0.562,
     *               0.750, 1.000, 1.333, 1.780, 2.370, 3.160, 4.217,
     *               5.0,
     *               5.620, 6.0, 6.5, 7.0, 7.500, 8.0, 10.00, 13.00 / 
C======================================================================
C      Defititions relative to Netcdf files
C======================================================================
C     This is the name of the data file we will create.
      character*(*) FILE_NAME
      parameter (FILE_NAME='GCM_data.nc')
      integer ncid !file ID    
C     Defining netCDF dimensions.
      integer NDIMS
      parameter (NDIMS=3) 
      character*(*) LAT_NAME, LON_NAME,NSOLTOT_NAME,
     * NSEAS_NAME, NPRES_NAME
      parameter(LAT_NAME='lat', LON_NAME='lon',NSOLTOT_NAME='nsoltot',
     *  NSEAS_NAME='nseas',NPRES_NAME='npres') !npress=np_interp
      integer lon_dimid, lat_dimid, nsoltot_dimid,
     * nseas_dimid,npres_dimid

C     Definine the "coordinate variables"
      integer lat_varid, lon_varid,nsoltot_varid,
     * nseas_varid,npres_varid
C     Define "longname" attributes. 
      character*(*) LONGNAMES    
      parameter (LONGNAMES = 'longname')
      character*(*)LAT_LN, LON_LN, NPRES_LN
      parameter (LAT_LN = 'latitude, degrees_north')
      parameter (LON_LN = ' longitude, degrees_east')
      parameter (NPRES_LN = 'pressure interpolation levels')
C======================================================================
C     1D var (nsoltot)  
C======================================================================
      integer num_1d
      parameter(num_1d=10)
      integer ::varid_1d(num_1d)
      integer dimids1d(NDIMS-2)
      character(len=20) :: NAME_1D(num_1d)
      data NAME_1D /'lssol',!1
     +'gasp',               !2
     +'gagt',               !3  
     +'gd',                 !4 
     +'gwv',                !5
     +'gwi',                !6
     +'hco2_n',             !7
     +'hco2_s',             !8
     +'vl1',                !9
     +'vl2'/                !10 
      character(len=75) :: LN_1D(num_1d) !LONGNAME 1D variables           
      data LN_1D /'solar longitude' ,                                   !1
     +'sol mean global mean surface pressure',                          !2
     +'sol mean global mean surface temperature',                       !3
     +'sol mean globally integrated dust (kg)',                         !4
     +'sol mean globally integrated water vapor (kg)',                  !5
     +'sol mean globally integrated water ice (kg)',                    !6
     +'sol mean n-hemispherically integrated surface CO2 ice (kg)',     !7
     +'sol mean s-hemispherically integrated surface CO2 ice (kg)',     !8
     +'sol mean hydrostatically adjusted surface pressure at VL1 (mb)', !9
     +'sol mean hydrostatically adjusted surface pressure at VL2 (mb)'/ !10

C======================================================================
C     1D var (npres)  
C======================================================================
      integer num_1d_p
      parameter(num_1d_p=1)
      integer ::varid_1d_p(num_1d_p)
      integer dimids1d_p(NDIMS-2)
      character(len=20) :: NAME_1D_P(num_1d_p)
      data NAME_1D_P /'pinterp'/ !1
      character(len=75) :: LN_1D_P(num_1d_p) !LONGNAME 1D variables           
      data LN_1D_P /'pressure grid (mb)'/        !1

C======================================================================
C     2D var (nsoltot,nlat)
C======================================================================
      integer num_2d_SL ! SL--> Sol and Latitude
      parameter(num_2d_SL=12)
      integer ::varid_2d_SL(num_2d_SL)
      integer dimids2d_SL(NDIMS-1) 
      character(len=20) :: NAME_2D_SL(num_2d_SL)
      data NAME_2D_SL /'zd',!1
     +'zwv',                !2 
     +'zwi',                !3
     +'zws',                !4
     +'zss',                !5
     +'zt',                 !6
     +'zt2pm',              !7
     +'zt2am',              !8
     +'zco2i',              !9
     +'z0pt5t',             !10
     +'z0pt5t2pm',          !11
     +'z0pt5t2am'/          !12
      character(len=75) :: LN_2D_SL(num_2d_SL)! LONGNAME
      data LN_2D_SL/
     +'sol mean zonal mean column integrated dust opacity at 9um',      !1 
     +'sol mean zonal mean column integrated water vapor (pr-um)',      !2
     +'sol mean zonal mean column integrated water ice (pr-um)',        !3
     +'sol mean zonal mean surface wind speed (m/s)',                   !4
     +'sol mean zonal mean surface stress magnitude (mN/m^2)',          !5
     +'sol mean zonal mean surface temperature (K)',                    !6
     +'1-3 pm mean zonal mean surface temperature (K)',                 !7
     +'1-3 am mean zonal mean surface temperature (K)',                 !8
     +'sol mean zonal mean surface CO2 ice (kg/m^2)',                   !9
     +'sol mean zonal mean 0.5 mb temperature (K)',                     !10
     +'1-3 pm mean zonal mean 0.5 mb temperature (K)',                  !11
     +'1-3 am mean zonal mean 0.5 mb temperature (K)'/                  !12


C======================================================================
C     2D var (nlat,nseas)
C======================================================================
      integer num_2d_LS ! SL--> Latitude and Season
      parameter(num_2d_LS=1)
      integer ::varid_2d_LS(num_2d_LS)
      integer dimids2d_LS(NDIMS-1)
      character(len=20) :: NAME_2D_LS(num_2d_LS)
      data NAME_2D_LS /'zpsrf'/ !1
      character(len=75) :: LN_2D_LS(num_2d_LS)! LONGNAME
      data LN_2D_LS/
     +'30 Ls mean zonal mean surface pressure (mb)'/      !1 

C======================================================================
C     3D var (nlat,npres,nseas)
C======================================================================
      integer num_3d_LPS ! LPS--> Lat, Press and Season
      parameter(num_3d_LPS=6)
      integer ::varid_3d_LPS(num_3d_LPS)
      integer dimids3d_LPS(NDIMS) 
      character(len=20) :: NAME_3D_LPS(num_3d_LPS)
      data NAME_3D_LPS /'zpu',   !1
     +'zpv',                     !2 
     +'zpt',                     !3
     +'zpt2pm',                  !4
     +'zpt2am',                  !5
     +'zpmsf'/                   !6
      character(len=75) :: LN_3D_LPS(num_3d_LPS)! LONGNAME
      data LN_3D_LPS/
     +'30 Ls mean zonal mean cross-section of U wind',                 !1 
     +'30 Ls mean zonal mean cross-section of V wind',                 !2
     +'30 Ls mean zonal mean cross-section of temperature (K)',        !3
     +'30 Ls 1-3 pm mean zonal mean cross-section of temperature (K)', !4
     +'30 Ls 1-3 am mean zonal mean cross-section of temperature (K)', !5
     +'30 Ls mean mass stream function (10^8 kg/s)'/                   !6
C======================================================================
C     3D var (nlon,nlat,nseas)
C======================================================================
      integer num_3d_LLS ! LLS--> Lon, Lat and Season
      parameter(num_3d_LLS=13)
      integer ::varid_3d_LLS(num_3d_LLS)
      integer dimids3d_LLS(NDIMS) 
      character(len=20) :: NAME_3D_LLS(num_3d_LLS)
      data NAME_3D_LLS /'uwnd',   !1
     +'vwnd',                     !2 
     +'ws',                       !3
     +'ustr',                     !4
     +'vstr',                     !5
     +'mstr',                     !6
     +'co2i',                     !7
     +'tsrf',                     !8
     +'tsrf2pm',                  !9
     +'tsrf2am',                  !10
     +'t0pt5',                    !11
     +'t0pt52pm',                 !12
     +'t0pt52am'/                 !13
      character(len=75) :: LN_3D_LLS(num_3d_LLS)! LONGNAME
      data LN_3D_LLS/
     +'30 Ls mean U wind at surface (m/s)',           !1
     +'30 Ls mean V wind at surface (m/s)',           !2
     +'30 Ls mean wind speed at surface (m/s)',       !3
     +'30 Ls mean zonal surface stress (mN/m^2)',     !4
     +'30 Ls mean meridional surface stress (mN/m^2)',!5
     +'30 Ls mean surface stress magnitude (mN/m^2)', !6
     +'30 Ls mean surface CO2 ice budget (kg/m^2)',   !7
     +'30 Ls mean surface temperature (K)',           !8
     +'30 Ls 1-3 pm mean surface temperature (K)',    !9
     +'30 Ls 1-3 am mean surface temperature (K)',    !10
     +'30 Ls mean 0.5 mb temperature (K)',            !11
     +'30 Ls 1-3 pm mean 0.5 mb temperature (K)',     !12
     +'30 Ls 1-3 am mean 0.5 mb temperature (K)'/     !13

C     Error handling.
      integer retval
C======================================================================
C     END Defititions relative to Netcdf files
C======================================================================

!  initialize some stuff here

      radius  = 3.393E+6
      pavesfc = 7.0

      pi = 4.*atan(1.)
      rad = pi/180.


c for 36x60 grid

      j1(1) = 22 ! VL-1 J. Corresponds to 20N
      i1(1) = 23 ! VL-1 I. Corresponds to 48W / -48 E

      j2(1) = 27 ! VL-2 J   Corresponds to 45N
      i2(1) = 53 ! VL-2 I.  Corresponds to 228 W / 132E


c now set the other three j,i pairs

      j1(2) = j1(1)
      j1(3) = j1(1)+1
      j1(4) = j1(1)+1
      i1(2) = i1(1)+1
      i1(3) = i1(1)
      i1(4) = i1(1)+1

      j2(2) = j2(1)
      j2(3) = j2(1)+1
      j2(4) = j2(1)+1
      i2(2) = i2(1)+1
      i2(3) = i2(1)
      i2(4) = i2(1)+1


C     get history file name and store it in fname.

      write(6,'(" ")')
      write(6,'("First history file number:  ")',advance='no')
      read(5,*) ifile1
      write(6,'(" Last history file number:  ")',advance='no')
      read(5,*) ifile2

      nsoltot=(ifile2-ifile1+1)*nsol

      do l=1,np_interp
        zoverh(l)=-log(pinterp(l)/pavesfc)
      enddo

      do nf = ifile1, ifile2
        write(fname,'("fort.11_",i4.4)') nf
        write(6,'(A)') fname

C     try to open the history file, stop if can't open file.

        open(unit=20,file=fname,status='old',form='unformatted',
     *       IOSTAT=ios)


        if(ios.ne.0) then
          write(6,'("Could not open data file ",A)') fname
          stop
        end if

        read(20) RUNNUM, JM, IM, LM, NL, ntrace, version

        if(nf.eq.ifile1) then

!  allocate arrays and initialize to zero (or 1.e20)

          call myalloc1(im,jm,lm,nl,ntrace)
          call myalloc2(im,jm,lm,nl,ntrace,nseas,nsoltot,np_interp)

          do j=1,jm-1
            lat(j) = -90.0 + (180.0/float(jm))*j
          end do
          do i=1,im
            lon(i)=-180.+(360./float(im))*i
          enddo
          do n=1,nseas
            lsseas(n) = (n-1)*360./float(nseas)
          end do

        end if

C     read header information.  These variables are written only once per file

        call readheader

        if (nf .eq. ifile1) then
         do j=1,jm
           area = area + dxyp(j) * float(im)
         enddo
         sigma = 0.0
         do l=1,lm
           k=2*l+3
           sigma(k) = sigma(k-2)+dsig(l)
         end do
         do k=4,2*lm+3,2
           sigma(k) = 0.5*(sigma(k+1)+sigma(k-1))
         end do
        endif

! loop over number of sols per file

        do nd=1,nsol

          nsolcnt=nsolcnt+1

! initialize arrays for sol means 

          do m=1,4
           psvl1(m) = 0.
           psvl2(m) = 0.
           gtvl1(m) = 0.
           gtvl2(m) = 0.
          end do

! loop over number of writeouts per sol

          do n=1,wsol

            call readvariables

            lssol(nsolcnt)=lssol(nsolcnt)+vpout/float(wsol)

            als = als + vpout/float(wsol)
            ati = ati + tau/float(wsol)

            lsbin=int((vpout+15.)/30.)+1
            if (lsbin .gt. 12) lsbin=1
            nlsbin(lsbin)=nlsbin(lsbin)+1.

            do j=1,jm-1
              do i=1,im

                 do l=1,lm
                   k=2*l+2
                   plx(l)=sigma(k)*p(j,i)+ptrop
                 enddo 

                 call localtime(im,tofday,i,tlocal)

                 gasp(nsolcnt)=gasp(nsolcnt)+p(j,i)*dxyp(j)/
     *                                       (area*float(wsol))
                 gagt(nsolcnt)=gagt(nsolcnt)+gt(j,i)*dxyp(j)/
     *                                       (area*float(wsol))
                 gd(nsolcnt)=gd(nsolcnt)+tausurf(j,i)*dxyp(j)/
     *                                       (area*float(wsol))

                 if (j .gt. int(jm/2.))
     *             hco2_n(nsolcnt)=hco2_n(nsolcnt)+co2ice(j,i)
                 if (j .le. int(jm/2.))
     *             hco2_s(nsolcnt)=hco2_s(nsolcnt)+co2ice(j,i)
                 if (j .eq. int(jm/2.)) then
                   hco2_n(nsolcnt)=hco2_n(nsolcnt)+co2ice(j,i)/2.
                   hco2_s(nsolcnt)=hco2_s(nsolcnt)+co2ice(j,i)/2.
                 endif

                 wvcol=0.
                 wicol=0.
                 do l=1,lm
                   wvcol=wvcol+qtrace(j,i,l,ima_vap)*1.e5*p(j,i)*
     *                         dsig(l)/grav
                   wicol=wicol+qtrace(j,i,l,ima_cld)*1.e5*p(j,i)*
     *                         dsig(l)/grav
                 enddo

                 gwv(nsolcnt)=gwv(nsolcnt)+wvcol*dxyp(j)/
     *                                      (area*float(wsol))
                 gwi(nsolcnt)=gwi(nsolcnt)+wicol*dxyp(j)/
     *                                      (area*float(wsol))
                 zd(nsolcnt,j)=zd(nsolcnt,j)+tausurf(j,i)/
     *                                      (float(im*wsol))
                 zwv(nsolcnt,j)=zwv(nsolcnt,j)+wvcol/
     *                                      (float(im*wsol))
                 zwi(nsolcnt,j)=zwi(nsolcnt,j)+wicol/
     *                                      (float(im*wsol))

                 if (i .eq. 1) then        ! put U winds on PI points
                        wu=0.5*(u(j,i,lm)+u(j,im,lm))
                 else
                        wu=0.5*(u(j,i,lm)+u(j,i-1,lm))
                 endif
                 if(j.gt.1 .and. j.lt.jm-1) then   ! put V winds on PI points
                      wv   = 0.5*(v(j+1,i,lm)+v(j,i,lm))
                 elseif(j.eq.1) then
                      wv   = 2.0*v(j+1,i,lm)/3.0
                 elseif(j.eq.jm-1) then
                      wv   = 2.0*v(j,i,lm)/3.0
                 endif

                 wuv = (wu**2. + wv**2.)**0.5
                 ss = (stressx(j,i)**2. + stressy(j,i)**2.)**0.5

                 zws(nsolcnt,j)=zws(nsolcnt,j)+wuv/
     *                                      (float(im*wsol))
                 zss(nsolcnt,j)=zss(nsolcnt,j)+ss/
     *                                      (float(im*wsol))
                 zt(nsolcnt,j)=zt(nsolcnt,j)+gt(j,i)/
     *                                      (float(im*wsol))
                 zco2i(nsolcnt,j)=zco2i(nsolcnt,j)+co2ice(j,i)/
     *                                      (float(im*wsol))

                 do l=1,lm-1
                   if(plx(l+1) .gt. 0.5 .and. plx(l) .le. 0.5) then
                     z0pt5t(nsolcnt,j) = z0pt5t(nsolcnt,j)+(t(j,i,l+1)+ 
     *                              (t(j,i,l)-t(j,i,l+1)) *  
     *                              alog(0.5/plx(l+1))/ 
     *                              alog(plx(l)/plx(l+1)))/
     *                                      (float(im*wsol))
                     exit
                   end if
                 end do

                 if(tlocal .gt. 13.0 .and. tlocal .lt. 15.0) then
                    zt2pm(nsolcnt,j)=zt2pm(nsolcnt,j)+gt(j,i)
                    znumt2pm(nsolcnt,j)=znumt2pm(nsolcnt,j)+1.
                    do l=1,lm-1
                      if(plx(l+1) .gt. 0.5 .and. plx(l) .le. 0.5) then
                        z0pt5t2pm(nsolcnt,j) = z0pt5t2pm(nsolcnt,j)+
     *                              (t(j,i,l+1)+ 
     *                              (t(j,i,l)-t(j,i,l+1)) *  
     *                              alog(0.5/plx(l+1))/ 
     *                              alog(plx(l)/plx(l+1)))
                        exit
                       end if
                    end do
                 elseif(tlocal .gt. 1.0 .and. tlocal .lt. 3.0) then
                    zt2am(nsolcnt,j)=zt2am(nsolcnt,j)+gt(j,i)
                    znumt2am(nsolcnt,j)=znumt2am(nsolcnt,j)+1.
                    do l=1,lm-1
                      if(plx(l+1) .gt. 0.5 .and. plx(l) .le. 0.5) then
                        z0pt5t2am(nsolcnt,j) = z0pt5t2am(nsolcnt,j)+
     *                              (t(j,i,l+1)+ 
     *                              (t(j,i,l)-t(j,i,l+1)) *  
     *                              alog(0.5/plx(l+1))/ 
     *                              alog(plx(l)/plx(l+1)))
                        exit
                       end if
                    end do
                 endif

! binning 30 degrees of Ls (centered at 0, 30, 60, etc)


                 zpsrf(j,lsbin)=zpsrf(j,lsbin)+p(j,i)/float(im)
                 uwnd(i,j,lsbin)=uwnd(i,j,lsbin)+wu
                 vwnd(i,j,lsbin)=vwnd(i,j,lsbin)+wv
                 ws(i,j,lsbin)=ws(i,j,lsbin)+wuv
                 ustr(i,j,lsbin)=ustr(i,j,lsbin)+stressx(j,i)
                 vstr(i,j,lsbin)=vstr(i,j,lsbin)+stressy(j,i)
                 mstr(i,j,lsbin)=mstr(i,j,lsbin)+ss
                 co2i(i,j,lsbin)=co2i(i,j,lsbin)+co2ice(j,i)
                 tsrf(i,j,lsbin)=tsrf(i,j,lsbin)+gt(j,i)

                 do l=1,lm-1
                   if(plx(l+1) .gt. 0.5 .and. plx(l) .le. 0.5) then
                     t0pt5(i,j,lsbin) = t0pt5(i,j,lsbin)+t(j,i,l+1)+ 
     *                              (t(j,i,l)-t(j,i,l+1)) *  
     *                              alog(0.5/plx(l+1))/ 
     *                              alog(plx(l)/plx(l+1))
                     exit
                   end if
                 end do

                 if(tlocal .gt. 13.0 .and. tlocal .lt. 15.0) then
                    tsrf2pm(i,j,lsbin)=tsrf2pm(i,j,lsbin)+gt(j,i)
                    numtsrf2pm(i,j,lsbin)=numtsrf2pm(i,j,lsbin)+1.
                    do l=1,lm-1
                      if(plx(l+1) .gt. 0.5 .and. plx(l) .le. 0.5) then
                        t0pt52pm(i,j,lsbin) = t0pt52pm(i,j,lsbin)+
     *                              t(j,i,l+1)+ 
     *                              (t(j,i,l)-t(j,i,l+1)) *  
     *                              alog(0.5/plx(l+1))/ 
     *                              alog(plx(l)/plx(l+1))
                        exit
                      end if
                    end do
                 elseif(tlocal .gt. 1.0 .and. tlocal .lt. 3.0) then
                    tsrf2am(i,j,lsbin)=tsrf2am(i,j,lsbin)+gt(j,i)
                    numtsrf2am(i,j,lsbin)=numtsrf2am(i,j,lsbin)+1.
                    do l=1,lm-1
                      if(plx(l+1) .gt. 0.5 .and. plx(l) .le. 0.5) then
                        t0pt52am(i,j,lsbin) = t0pt52am(i,j,lsbin)+
     *                              t(j,i,l+1)+ 
     *                              (t(j,i,l)-t(j,i,l+1)) *  
     *                              alog(0.5/plx(l+1))/ 
     *                              alog(plx(l)/plx(l+1))
                        exit
                      end if
                    end do
                 endif

! pressure interpolation here (zonal mean fields)

                 psx          = P(J,I)+PTROP
 
                 if(i.eq.1) then
                   iu1=im
                   iu2=1
                 else
                   iu1=i-1
                   iu2=i
                 endif

                 do lp=1,np_interp
                   if(pinterp(lp).lt.ptrop) then
                     zpt(j,lp,lsbin)=zpt(j,lp,lsbin)+tstrat(j,i)
                     numt(j,lp,lsbin)=numt(j,lp,lsbin)+1
                     if(tlocal .gt. 13.0 .and. tlocal .lt. 15.0) then
                       zpt2pm(j,lp,lsbin)=zpt2pm(j,lp,lsbin)+
     *                                    tstrat(j,i)
                       numt2pm(j,lp,lsbin)=numt2pm(j,lp,lsbin)+1
                     elseif(tlocal .gt. 1.0 .and. tlocal .lt. 3.0) then
                       zpt2am(j,lp,lsbin)=zpt2am(j,lp,lsbin)+
     *                                    tstrat(j,i)
                       numt2am(j,lp,lsbin)=numt2am(j,lp,lsbin)+1
                     endif
                   elseif(pinterp(lp).lt.psx .and. 
     *                    pinterp(lp).ge.plx(lm)) then
                     zpt(j,lp,lsbin)=zpt(j,lp,lsbin)+t(j,i,lm)+
     *                         (t(j,i,lm-1)-t(j,i,lm))*
     *                         alog(pinterp(lp)/plx(lm))/
     *                         alog(plx(lm-1)/plx(lm))
                     numt(j,lp,lsbin)=numt(j,lp,lsbin)+1
                     if(tlocal .gt. 13.0 .and. tlocal .lt. 15.0) then 
                       zpt2pm(j,lp,lsbin)=zpt2pm(j,lp,lsbin)+
     *                          t(j,i,lm) + (t(j,i,lm-1)-t(j,i,lm))* 
     *                          alog(pinterp(lp)/plx(lm))/
     *                          alog(plx(lm-1)/plx(lm))
                       numt2pm(j,lp,lsbin)=numt2pm(j,lp,lsbin)+1
                     elseif(tlocal .gt. 1.0 .and. tlocal .lt. 3.0) then
                       zpt2am(j,lp,lsbin)=zpt2am(j,lp,lsbin)+
     *                          t(j,i,lm) + (t(j,i,lm-1)-t(j,i,lm))* 
     *                          alog(pinterp(lp)/plx(lm))/
     *                          alog(plx(lm-1)/plx(lm))
                       numt2am(j,lp,lsbin)=numt2am(j,lp,lsbin)+1
                     endif

                     uavgl   = 0.5*(u(j,iu1,lm-1)+u(j,iu2,lm-1))
                     uavglp1 = 0.5*(u(j,iu1,lm)+u(j,iu2,lm))
                     if(j.gt.1 .and. j.lt.jm-1) then
                      vavgl   = 0.5*(v(j+1,i,lm-1)+v(j,i,lm-1))
                      vavglp1 = 0.5*(v(j+1,i,lm)+v(j,i,lm))
                     elseif(j.eq.1) then
                      vavgl   = 2.0*v(j+1,i,lm-1)/3.0
                      vavglp1 = 2.0*v(j+1,i,lm)/3.0
                     elseif(j.eq.jm-1) then
                      vavgl   = 2.0*v(j,i,lm-1)/3.0
                      vavglp1 = 2.0*v(j,i,lm)/3.0
                     endif
                     zpu(j,lp,lsbin)=zpu(j,lp,lsbin)+uavglp1+
     *                          (uavgl-uavglp1)*
     *                          (pinterp(lp)-plx(lm))/
     *                          (plx(lm-1)-plx(lm))
                     zpv(j,lp,lsbin)=zpv(j,lp,lsbin)+vavglp1+
     *                         (vavgl-vavglp1)*
     *                         (pinterp(lp)-plx(lm))/
     *                         (plx(lm-1)-plx(lm))
                     numw(j,lp,lsbin)=numw(j,lp,lsbin)+1
                    else
                     do l=1,lm-1
                       if(pinterp(lp).lt.plx(l+1) .and. 
     *                    pinterp(lp).ge.plx(l)) then
                         zpt(j,lp,lsbin) = zpt(j,lp,lsbin)+t(j,i,l+1)+ 
     *                              (t(j,i,l)-t(j,i,l+1)) *  
     *                              alog(pinterp(lp)/plx(l+1))/ 
     *                              alog(plx(l)/plx(l+1))
                         numt(j,lp,lsbin) = numt(j,lp,lsbin)+1
                         if(tlocal.gt.13.0 .and. tlocal.lt.15.0) then 
                           zpt2pm(j,lp,lsbin) = zpt2pm(j,lp,lsbin) +
     *                              t(j,i,l+1)+(t(j,i,l)-t(j,i,l+1))*  
     *                              alog(pinterp(lp)/plx(l+1))/ 
     *                              alog(plx(l)/plx(l+1))
                           numt2pm(j,lp,lsbin) = numt2pm(j,lp,lsbin)+1
                         elseif(tlocal.gt.1.0 .and. tlocal.lt.3.0) then
                           zpt2am(j,lp,lsbin) = zpt2am(j,lp,lsbin) +
     *                              t(j,i,l+1)+(t(j,i,l)-t(j,i,l+1))*  
     *                              alog(pinterp(lp)/plx(l+1))/ 
     *                              alog(plx(l)/plx(l+1))
                           numt2am(j,lp,lsbin) = numt2am(j,lp,lsbin)+1

                         endif
                         uavgl   = 0.5*(u(j,iu1,l)+u(j,iu2,l))
                         uavglP1 = 0.5*(u(j,iu1,l+1)+u(j,iu2,l+1))
                         if(j.gt.1 .and. j.lt.jm-1) then
                           vavgl   = 0.5*(v(j+1,i,l)+v(j,i,l))
                           vavglP1 = 0.5*(v(j+1,i,l+1)+v(j,i,l+1))
                         elseif(j.eq.1) then
                           vavgl   = 2.0*v(j+1,i,l)/3.0
                           vavgLP1 = 2.0*v(j+1,i,l+1)/3.0
                         elseif(j.eq.jm-1) then
                           vavgl   = 2.0*v(j,i,l)/3.0
                           vavgLP1 = 2.0*v(j,i,l+1)/3.0
                         endif
                         zpu(j,lp,lsbin)=zpu(j,lp,lsbin)+uavglp1+
     *                              (uavgl-uavglp1)*
     *                              (pinterp(lp)-plx(l+1))/
     *                              (plx(L)-plx(l+1))
                         zpv(j,lp,lsbin)=zpv(j,lp,lsbin)+vavglp1+
     *                              (vavgl-vavglp1)*
     *                              (pinterp(Lp)-Plx(L+1))/
     *                              (Plx(L)-Plx(L+1))
                         numw(j,lp,lsbin)=numw(j,lp,lsbin)+1
                         exit
                       end if
                     end do
                   endif
                end do

              end do  ! end i loop
            end do    ! end j loop

! VL 1 and VL 2 Surface Pressures

            do m=1,4
              psvl1(m) = psvl1(m) + p(j1(m),i1(m))/float(wsol)
              psvl2(m) = psvl2(m) + p(j2(m),i2(m))/float(wsol)
              gtvl1(m) = gtvl1(m) + gt(j1(m),i1(m))/float(wsol)
              gtvl2(m) = gtvl2(m) + gt(j2(m),i2(m))/float(wsol)
            end do

          end do   ! end n=1,wday loop
            

C Interpolate to get Pressures at Viking Sites
 
          do l=1,4
            vl1gpp(l)   = psvl1(l)
            vl2gpp(l)   = psvl2(l)
            vl1gpt(l)   = gtvl1(l)
            vl2gpt(l)   = gtvl2(l)
          enddo

C First Correct GCM grid point pressures for Elevation

          do l = 1,4
            h1 = 189.*vl1gpt(l)/3.72/1.0e3
            h2 = 189.*vl2gpt(l)/3.72/1.0e3
            d1 = -3.627-vl1z(l)
            d2 = -4.505-vl2z(l)
            vl1gpp(l) = vl1gpp(l)*exp(-d1/h1)
            vl2gpp(l) = vl2gpp(l)*exp(-d2/h2)
          end do

C Now interpoate in latitude and longitude

          fact1 = (0.03/6.0)*cos(20.0*rad)
          fact2 = (0.03/6.0)*cos(25.0*rad)
          fact3 = (2.48/5.)
          term1 = vl1gpp(1) + fact1*(vl1gpp(2)-vl1gpp(1))
          term2 = vl1gpp(3) + fact2*(vl1gpp(4)-vl1gpp(3))
          vlp(1) = term1 + fact3*(term2-term1)
          term1 = vl1gpt(1) + fact1*(vl1gpt(2)-vl1gpt(1))
          term2 = vl1gpt(3) + fact2*(vl1gpt(4)-vl1gpt(3))
          vlt(1) = term1 + fact3*(term2-term1)

          fact1 = (2.26/6.0)*cos(45.0*rad)
          fact2 = (2.26/6.0)*cos(50.0*rad)
          fact3 = (2.97/5.0)
          term1 = vl2gpp(1) + fact1*(vl2gpp(2)-vl2gpp(1))
          term2 = vl2gpp(3) + fact2*(vl2gpp(4)-vl2gpp(3))
          vlp(2) = term1 + fact3*(term2-term1)
          term1 = vl2gpt(1) + fact1*(vl2gpt(2)-vl2gpt(1))
          term2 = vl2gpt(3) + fact2*(vl2gpt(4)-vl2gpt(3))
          vlt(2) = term1 + fact3*(term2-term1)
      
          vl1(nsolcnt)=vlp(1)
          vl2(nsolcnt)=vlp(2)

        end do   ! end nd=1,nsol loop

      end do  ! end nf=ifile1,ifile2 loop

!  complete averaging

      do n=1,nseas
        do j=1,jm-1
          zpsrf(j,n)=zpsrf(j,n)/nlsbin(n)
          do i = 1,im
            uwnd(i,j,n)=uwnd(i,j,n)/nlsbin(n)
            vwnd(i,j,n)=vwnd(i,j,n)/nlsbin(n)
            ws(i,j,n)=ws(i,j,n)/nlsbin(n)
            ustr(i,j,n)=ustr(i,j,n)/nlsbin(n)
            vstr(i,j,n)=vstr(i,j,n)/nlsbin(n)
            mstr(i,j,n)=mstr(i,j,n)/nlsbin(n)
            co2i(i,j,n)=co2i(i,j,n)/nlsbin(n)
            tsrf(i,j,n)=tsrf(i,j,n)/nlsbin(n)
            tsrf2pm(i,j,n)=tsrf2pm(i,j,n)/numtsrf2pm(i,j,n)
            tsrf2am(i,j,n)=tsrf2am(i,j,n)/numtsrf2am(i,j,n)
            t0pt5(i,j,n)=t0pt5(i,j,n)/nlsbin(n)
            t0pt52pm(i,j,n)=t0pt52pm(i,j,n)/numtsrf2pm(i,j,n)
            t0pt52am(i,j,n)=t0pt52am(i,j,n)/numtsrf2am(i,j,n)
          enddo
          cosy=cos(lat(j)*rad)
          do l=1,np_interp
            if (numt(j,l,n) .gt. 0.) then
               zpt(j,l,n)=zpt(j,l,n)/numt(j,l,n)
            else
               zpt(j,l,n) = 1.e20
            endif
            if (numt2pm(j,l,n) .gt. 0.) then
               zpt2pm(j,l,n)=zpt2pm(j,l,n)/numt2pm(j,l,n)
            else
               zpt2pm(j,l,n) = 1.e20
            endif
            if (numt2am(j,l,n) .gt. 0.) then
               zpt2am(j,l,n)=zpt2am(j,l,n)/numt2am(j,l,n)
            else
               zpt2am(j,l,n) = 1.e20
            endif
            if (numw(j,l,n) .gt. 0.) then
               zpu(j,l,n) = zpu(j,l,n)/numw(j,l,n)
               zpv(j,l,n) = zpv(j,l,n)/numw(j,l,n)
            else
               zpu(j,l,n) = 1.e20
               zpv(j,l,n) = 1.e20
            endif
          enddo
          do k=2,np_interp
            sum = 0.0
            do l=2,k
              zl   = zoverh(l)
              zlp1 = zoverh(l-1)
              if(numw(j,l,n).gt.0 .and. numw(j,l-1,n).gt.0) then
                fl   = zpv(j,l,n)*exp(-zl)
                flp1 = zpv(j,l-1,n)*exp(-zlp1)
                sum  = sum + 0.5*(zl-zlp1)*(fl+flp1)*pavesfc*100.0*
     *                       1.0E-8
              end if
            end do
            if(sum.ne.0.0) then
              zpmsf(j,k,n) = (2.0*pi*radius*cosy/grav)*sum
            else
              zpmsf(j,k,n) = 1.e20
            end if
          end do
        enddo
      enddo

      nyr = -1

      do n=1,nsolcnt
        do j=1,jm-1
          zt2pm(n,j)=zt2pm(n,j)/znumt2pm(n,j)
          zt2am(n,j)=zt2am(n,j)/znumt2am(n,j)
          z0pt5t2pm(n,j)=z0pt5t2pm(n,j)/znumt2pm(n,j)
          z0pt5t2am(n,j)=z0pt5t2am(n,j)/znumt2am(n,j)
        enddo
        if (n .eq.1) then
           if (lssol(n) .lt. 300.) then
             nyr=nyr+1
           endif
        else
           if((lssol(n-1).gt.(360*(nyr+1)-60)) .and. (lssol(n).lt.60.0))
     *       nyr=nyr+1
        endif
        lssol(n)=lssol(n)+360.*float(nyr)
      enddo 

! write to files

!      write(40) nsolcnt
!      write(40) lssol  
!      write(40) gasp
!      write(40) gagt
!      write(40) gd
!      write(40) gwv
!      write(40) gwi
!      write(40) hco2_n
!      write(40) hco2_s
!      write(40) vl1
!      write(40) vl2

!      write(41) nsolcnt
!      write(41) jm
!      write(41) lssol
!      write(41) lat
!      write(41) zd
!      write(41) zwv
!      write(41) zwi
!      write(41) zws
!      write(41) zss
!      write(41) zt
!      write(41) zt2pm
!      write(41) zt2am
!      write(41) zco2i
!      write(41) z0pt5t
!      write(41) z0pt5t2pm
!      write(41) z0pt5t2am

!      write(42) nseas
!      write(42) jm
!      write(42) im
!      write(42) lsseas
!      write(42) nlsbin
!      write(42) lat
!      write(42) lon
!      write(42) uwnd
!      write(42) vwnd
!      write(42) ws
!      write(42) ustr
!      write(42) vstr
!      write(42) mstr
!      write(42) co2i
!      write(42) tsrf
!      write(42) tsrf2pm
!      write(42) tsrf2am
!      write(42) t0pt5
!      write(42) t0pt52pm
!      write(42) t0pt52am

!      write(43) nseas
!      write(43) jm
!      write(43) np_interp
!      write(43) pinterp
!      write(43) zpt
!      write(43) zpt2pm
!      write(43) zpt2am
!      write(43) zpu
!      write(43) zpv
!      write(43) zpmsf
!      write(43) zprsf

C======================================================================      
C     Create the file.
C======================================================================
      retval = nf_create(FILE_NAME, nf_clobber, ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

C======================================================================      
C     Define the dimensions.
C======================================================================
      retval = nf_def_dim(ncid, NSOLTOT_NAME, nsoltot, nsoltot_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid, LAT_NAME, JM-1, lat_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid, LON_NAME, IM, lon_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid, NSEAS_NAME, nseas, nseas_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid, NPRES_NAME, np_interp, npres_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
C======================================================================
C     Define the coordinate variables. They will hold the coordinate
C     information: the latitudes and longitudes. A varid is
C     returned for each.
C======================================================================
      write(*,*)'Defining coordinate variables...'
C     Latitude
      retval = nf_def_var(ncid, LAT_NAME, NF_REAL, 1, lat_dimid, 
     +     lat_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, lat_varid, LONGNAMES, len(LAT_LN), 
     +     LAT_LN)
      if (retval .ne. nf_noerr) call handle_err(retval)
      write(*,*)'i=1','NAME: ',LAT_NAME,', varid=',
     +   lat_varid,'|', LAT_LN
C     Longitude
      retval = nf_def_var(ncid, LON_NAME, NF_REAL, 1, lon_dimid, 
     +     lon_varid)     
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, lon_varid, LONGNAMES, len(LON_LN), 
     +     LON_LN)
      if (retval .ne. nf_noerr) call handle_err(retval)
      write(*,*)'i=2','NAME: ',LON_NAME,', varid=',
     +   lon_varid, '|',LON_LN
C     Pressure interpolation levels
      retval = nf_def_var(ncid, NPRES_NAME, NF_REAL, 1, npres_dimid, 
     +     npres_varid)     
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, npres_varid, LONGNAMES,
     + len(NPRES_LN),    NPRES_LN)
      if (retval .ne. nf_noerr) call handle_err(retval)
      write(*,*)'i=3','NAME: ',NPRES_NAME,', varid=',
     +   npres_varid, '|',NPRES_LN
C======================================================================
C     Define the netCDF variables. The dimids array is used to pass the
C     dimids of the dimensions of the netCDF variables.
C======================================================================    
C     1D variables (nsoltot)
      write(*,*)'Defining 1D variables...(nsoltot)'
      dimids1d = nsoltot_dimid
C     Define the netCDF variables.
      do i = 1, num_1d
         retval = nf_def_var(ncid, NAME_1D(i), NF_REAL,
     +   NDIMS-2, dimids1d, varid_1d(i))
         if (retval .ne. nf_noerr) call handle_err(retval)
C     Assign longname attributes 
         retval = nf_put_att_text(ncid, varid_1d(i), LONGNAMES, 
     +	 len(trim(LN_1D(i))), trim(LN_1D(i)))
         if (retval .ne. nf_noerr) call handle_err(retval)     
         write(*,*)'i=',i,'NAME:',trim(NAME_1D(i)),', varid=',
     +   varid_1d(i),'|',LN_1D(i)     
      end  do
C     1D variables (npres)
      write(*,*)'Defining 1D variables...(npres)'
      dimids1d_p = npres_dimid
C     Define the netCDF variables.
      do i = 1, num_1d_p
         retval = nf_def_var(ncid, NAME_1D_P(i), NF_REAL,
     +   NDIMS-2, dimids1d_p, varid_1d_p(i))
         if (retval .ne. nf_noerr) call handle_err(retval)
C     Assign longname attributes 
         retval = nf_put_att_text(ncid, varid_1d_p(i), LONGNAMES,
     +   len(trim(LN_1D_P(i))), trim(LN_1D_P(i)))
         if (retval .ne. nf_noerr) call handle_err(retval)
         write(*,*)'i=',i,'NAME:',trim(NAME_1D_P(i)),', varid=',
     +   varid_1d_p(i),'|',LN_1D_P(i)
      end  do

C-----------------------------------------------
C     2D variables (nsoltot,nlat)
      write(*,*)'Defining 2D variables... (nsoltot,nlat)'
      dimids2d_SL(1) = nsoltot_dimid 
      dimids2d_SL(2) = lat_dimid 
C     Define the netCDF variables. 
      do i = 1, num_2d_SL
         retval = nf_def_var(ncid, NAME_2D_SL(i), NF_REAL,
     +   NDIMS-1, dimids2d_SL, varid_2d_SL(i))
         if (retval .ne. nf_noerr) call handle_err(retval)
C     Assign longname attributes 
         retval = nf_put_att_text(ncid, varid_2d_SL(i), LONGNAMES, 
     +	 len(trim(LN_2D_SL(i))), trim(LN_2D_SL(i)))
         if (retval .ne. nf_noerr) call handle_err(retval)     
         write(*,*)'i=',i,'NAME: ',trim(NAME_2D_SL(i)),', varid=',
     +   varid_2d_SL(i),'|', LN_2D_SL(i)     
      end  do

      write(*,*)'Defining 2D variables... (nlat,nseas)'
      dimids2d_LS(1) = lat_dimid
      dimids2d_LS(2) = nseas_dimid
C     Define the netCDF variables. 
      do i = 1, num_2d_LS
         retval = nf_def_var(ncid, NAME_2D_LS(i), NF_REAL,
     +   NDIMS-1, dimids2d_LS, varid_2d_LS(i))
         if (retval .ne. nf_noerr) call handle_err(retval)
C     Assign longname attributes 
         retval = nf_put_att_text(ncid, varid_2d_LS(i), LONGNAMES,
     +   len(trim(LN_2D_LS(i))), trim(LN_2D_LS(i)))
         if (retval .ne. nf_noerr) call handle_err(retval)
         write(*,*)'i=',i,'NAME: ',trim(NAME_2D_LS(i)),', varid=',
     +   varid_2d_LS(i),'|', LN_2D_LS(i)
      end  do

C----------------------------------------------
C     3D variables (nlat,npres,nseas)
      write(*,*)'Defining 3D variables (nlat,npres,nseas)'
      dimids3d_LPS(1) = lat_dimid 
      dimids3d_LPS(2) = npres_dimid 
      dimids3d_LPS(3) = nseas_dimid 
C     Define the netCDF variables. 
      do i = 1, num_3d_LPS
         retval = nf_def_var(ncid, NAME_3D_LPS(i), NF_REAL,
     +   NDIMS, dimids3d_LPS, varid_3d_LPS(i))
         if (retval .ne. nf_noerr) call handle_err(retval)
C     Assign longname attributes 
         retval = nf_put_att_text(ncid, varid_3d_LPS(i), LONGNAMES, 
     +	 len(trim(LN_3D_LPS(i))), trim(LN_3D_LPS(i)))
         if (retval .ne. nf_noerr) call handle_err(retval)     
         write(*,*)'i=',i,'NAME: ',trim(NAME_3D_LPS(i)),', varid=',
     +   varid_3d_LPS(i),'|', LN_3D_LPS(i)     
      end  do
C----------------------------------------------
C     3D variables (nlon,nlat,nseas)
      write(*,*)'Defining 3D variables (nlon,nlat,nseas)'
      dimids3d_LLS(1) = lon_dimid 
      dimids3d_LLS(2) = lat_dimid 
      dimids3d_LLS(3) = nseas_dimid 
C     Define the netCDF variables. 
      do i = 1, num_3d_LLS
         retval = nf_def_var(ncid, NAME_3D_LLS(i), NF_REAL,
     +   NDIMS, dimids3d_LLS, varid_3d_LLS(i))
         if (retval .ne. nf_noerr) call handle_err(retval)
C     Assign longname attributes 
         retval = nf_put_att_text(ncid, varid_3d_LLS(i), LONGNAMES, 
     +	 len(trim(LN_3D_LLS(i))), trim(LN_3D_LLS(i)))
         if (retval .ne. nf_noerr) call handle_err(retval)     
         write(*,*)'i=',i,'NAME: ',trim(NAME_3D_LLS(i)),', varid=',
     +   varid_3d_LLS(i),'|',LN_3D_LLS(i)     
      end  do
C======================================================================
C     End define mode.
C======================================================================
      retval = nf_enddef(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
C======================================================================
C     Write the coordinate variable data. 
C======================================================================
C     Latitude
      retval = nf_put_var_real(ncid, lat_varid, lat)
      if (retval .ne. nf_noerr) call handle_err(retval)
C     Longitude
      retval = nf_put_var_real(ncid, lon_varid, lon)
      if (retval .ne. nf_noerr) call handle_err(retval)
C     Pressure levels
      retval = nf_put_var_real(ncid, npres_varid, pinterp)
      if (retval .ne. nf_noerr) call handle_err(retval)
C     note: nseas does not have coordinate variable data
C======================================================================
C     Write the data.
C======================================================================
C     1D variables (nsoltot)      
      write(*,*)'Writing 1D variables (ntotsol)...'
      retval = nf_put_var_real(ncid, varid_1d(1),lssol)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, varid_1d(2),gasp)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, varid_1d(3),gagt)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, varid_1d(4),gd)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, varid_1d(5),gwv)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, varid_1d(6),gwi)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, varid_1d(7),hco2_n)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, varid_1d(8),hco2_s)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, varid_1d(9),vl1)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, varid_1d(10),vl2)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, varid_1d_p(1),pinterp)
      if (retval .ne. nf_noerr) call handle_err(retval)
      write(*,*)'Writing 2D variables (nsoltot,nlat)...'
C     2D variables (nsoltot,nlat)
      retval = nf_put_var_real(ncid, varid_2d_SL(1),zd)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, varid_2d_SL(2),zwv)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, varid_2d_SL(3),zwi)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, varid_2d_SL(4),zws)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, varid_2d_SL(5),zss)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, varid_2d_SL(6),zt)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, varid_2d_SL(7),zt2pm)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, varid_2d_SL(8),zt2am)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, varid_2d_SL(9),zco2i)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, varid_2d_SL(10),z0pt5t)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, varid_2d_SL(11),z0pt5t2pm)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, varid_2d_SL(12),z0pt5t2pm)
      if (retval .ne. nf_noerr) call handle_err(retval)
      write(*,*)'Writing 2D variables (nlat,nseas)...'
C     2D variables (nsoltot,nlat)
      retval = nf_put_var_real(ncid, varid_2d_LS(1),zpsrf)
      if (retval .ne. nf_noerr) call handle_err(retval)
C     3D variables (nlat,npres,nseas):
      write(*,*)'Writing 3D variables (nlat,npres,nseas)...'
      retval = nf_put_var_real(ncid, varid_3d_LPS(1),zpu)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, varid_3d_LPS(2),zpv)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, varid_3d_LPS(3),zpt)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, varid_3d_LPS(4),zpt2pm)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, varid_3d_LPS(5),zpt2am)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, varid_3d_LPS(6),zpmsf)
      if (retval .ne. nf_noerr) call handle_err(retval)
C     3D variables (nlat,npres,nseas):
      write(*,*)'Writing 3D variables (nlon,nlat,nseas)...'
      retval = nf_put_var_real(ncid, varid_3d_LLS(1),uwnd)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, varid_3d_LLS(2),vwnd)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, varid_3d_LLS(3),ws)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, varid_3d_LLS(4),ustr)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, varid_3d_LLS(5),vstr)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, varid_3d_LLS(6),mstr)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, varid_3d_LLS(7),co2i)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, varid_3d_LLS(8),tsrf)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, varid_3d_LLS(9),tsrf2pm)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, varid_3d_LLS(10),tsrf2am)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, varid_3d_LLS(11),t0pt5)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, varid_3d_LLS(12),t0pt52pm)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, varid_3d_LLS(13),t0pt52am)
      if (retval .ne. nf_noerr) call handle_err(retval)
C======================================================================
C     Close the file.
C======================================================================
      retval = nf_close(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
   
C     If we got this far, everything worked as expected. Yipee!
      print *,'*** SUCCESS writing example file ',FILE_NAME
C======================================================================
C     END Netcdf 
C======================================================================
      end

      subroutine localtime(IM,tofday,Igrid,tlocal)
      implicit none
      real*4  :: tofday, tlocal, dtime
      integer :: igrid, IM, IMo2
               
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
C======================================================================
C     Netcdf subroutine defintion
C======================================================================
      subroutine handle_err(errcode)
      implicit none
      include 'netcdf.inc'
      integer errcode
      print *, 'Error: ', nf_strerror(errcode)
      stop 2
      end
