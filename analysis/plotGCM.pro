pro plotGCM

readfort,lssol,gasp,gagt,gd,gwv,gwi,hco2_n,hco2_s,vl1,vl2,lat,zd,zwv,$
zwi,zws,zss,zt,zt2pm,zt2am,zco2i,z0pt5t,z0pt5t2pm,z0pt5t2am,lon,lsseas,$
nlsbin,uwnd,vwnd,ws,ustr,vstr,mstr,co2i,tsrf,tsrf2pm,tsrf2am,t0pt5,$
t0pt52pm,t0pt52am,plevels,zpt,zpt2pm,zpt2am,zpu,zpv,zpmsf,zpsrf,$
nlat,nlon,nsolcnt,nseas,npres


; readnc,lssol,gasp,gagt,gd,gwv,gwi,hco2_n,hco2_s,vl1,vl2,lat,zd,zwv,$
; zwi,zws,zss,zt,zt2pm,zt2am,zco2i,z0pt5t,z0pt5t2pm,z0pt5t2am,lon,lsseas,$
; nlsbin,uwnd,vwnd,ws,ustr,vstr,mstr,co2i,tsrf,tsrf2pm,tsrf2am,t0pt5,$
; t0pt52pm,t0pt52am,plevels,zpt,zpt2pm,zpt2am,zpu,zpv,zpmsf,zpsrf,$
; nlat,nlon,nsolcnt,nseas,npres

print,' '
print,'Assumes that there is directory called "plots" in your current directory'
print,' '

path='plots/'

set_plot,'ps'
!p.font=0.

!p.multi=[0,1,2]
device,file=path+'gavgPsTs.eps',/inch,xsize=7.5,ysize=9.5,/encaps,/helvet
plot,lssol,gasp,yrange=[5.,9.],/ystyle,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,xtitle='Ls',title='Global Average Surface Pressure'
plot,lssol,gagt,yrange=[150.,250.],/ystyle,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,xtitle='Ls',title='Global Average Surface Temperature'
device,/close

!p.multi=[0,1,3]
device,file=path+'gavgDWvWi.eps',/inch,xsize=7.5,ysize=9.5,/encaps,/helvet
plot,lssol,gd,xthick=2.5,ythick=2.5,thick=2.5,charsize=1.5,$
     charthick=2.5,xtitle='Ls',title='Global Average Dust Visible Opacity'
plot,lssol,gwv,xthick=2.5,ythick=2.5,thick=2.5,charsize=1.5,$
     charthick=2.5,xtitle='Ls',title='Global Average Water Vapor'
plot,lssol,gwi,xthick=2.5,ythick=2.5,thick=2.5,charsize=1.5,$
     charthick=2.5,xtitle='Ls',title='Global Average Water Ice'
device,/close

!p.multi=[0,1,2]
device,file=path+'htotCO2i.eps',/inch,xsize=7.5,ysize=9.5,/encaps,/helvet
plot,lssol,hco2_n,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,xtitle='Ls',title='Hemispheric Total Surface CO2 Ice'
oplot,lssol,hco2_s,line=2,thick=2.5
plot,lssol,hco2_n+hco2_s,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,xtitle='Ls',title='Global Total Surface CO2 Ice'
device,/close

!p.multi=[0,1,2]
device,file=path+'PsVL.eps',/inch,xsize=7.5,ysize=9.5,/encaps,/helvet
plot,lssol,vl1,psym=1,yrange=[6.5,9.5],/ystyle,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,xtitle='Ls',title='Adjusted Surface Pressure at VL1'
plot,lssol,vl2,psym=2,yrange=[7.0,10.],/ystyle,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,xtitle='Ls',title='Adjusted Surface Pressure at VL2'
device,/close

levdst=findgen(10)*2/10.
levwv=findgen(30)*5
levwi=levwv/10.

!p.multi=[0,1,3]
device,file=path+'zavgDWvWi.eps',/inch,xsize=7.5,ysize=9.5,/encaps,/helvet
contour,zd,lssol,lat,lev=levdst,yrange=[-90.,90],/ystyle,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,xtitle='Ls',ytitle='Latitude',title='Zonal Average Dust Visible Opacity',$
     charsize=1.5,c_lab=levdst*0.+1
contour,zwv,lssol,lat,lev=levwv,yrange=[-90.,90],/ystyle,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,xtitle='Ls',ytitle='Latitude',title='Zonal Average Water Vapor',$
     charsize=1.5,c_lab=levwv*0.+1
contour,zwi,lssol,lat,lev=levwi,yrange=[-90.,90],/ystyle,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,xtitle='Ls',ytitle='Latitude',title='Zonal Average Water Ice',$
     charsize=1.5,c_lab=levwi*0.+1
device,/close

levws=findgen(20)*2.
levss=findgen(20)*4.
levco2i=findgen(20)*100.

!p.multi=[0,1,3]
device,file=path+'zavgWsSsCO2i.eps',/inch,xsize=7.5,ysize=9.5,/encaps,/helvet
contour,zws,lssol,lat,lev=levws,yrange=[-90.,90],/ystyle,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,xtitle='Ls',ytitle='Latitude',title='Zonal Average Surface Wind Speed',$
     charsize=1.5,c_lab=levws*0.+1
contour,zss*1000.,lssol,lat,lev=levss,yrange=[-90.,90],/ystyle,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,xtitle='Ls',ytitle='Latitude',title='Zonal Average Surface Stress Magnitude',$
     charsize=1.5,c_lab=levss*0.+1
contour,zco2i,lssol,lat,lev=levco2i,yrange=[-90.,90],/ystyle,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,xtitle='Ls',ytitle='Latitude',title='Zonal Average Surface CO2 Ice',$
     charsize=1.5,c_lab=levco2i*0.+1
device,/close

levts=findgen(20)*10+120.

!p.multi=[0,1,3]
device,file=path+'zavgTs.eps',/inch,xsize=7.5,ysize=9.5,/encaps,/helvet
contour,zt,lssol,lat,lev=levts,yrange=[-90.,90],/ystyle,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,xtitle='Ls',ytitle='Latitude',title='Zonal Average Surface Temperature',$
     charsize=1.5,c_lab=levts*0.+1
contour,zt2pm,lssol,lat,lev=levts,yrange=[-90.,90],/ystyle,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,xtitle='Ls',ytitle='Latitude',title='Zonal Average 1-3 PM Surface Temperature',$
     charsize=1.5,c_lab=levts*0.+1
contour,zt2am,lssol,lat,lev=levts,yrange=[-90.,90],/ystyle,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,xtitle='Ls',ytitle='Latitude',title='Zonal Average 1-3 AM Surface Temperature',$
     charsize=1.5,c_lab=levts*0.+1
device,/close

!p.multi=[0,1,3]
device,file=path+'zavgT0pt5mb.eps',/inch,xsize=7.5,ysize=9.5,/encaps,/helvet
contour,z0pt5t,lssol,lat,lev=levts,yrange=[-90.,90],/ystyle,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,xtitle='Ls',ytitle='Latitude',title='Zonal Average 0.5 hPa Temperature',$
     charsize=1.5,c_lab=levts*0.+1
contour,z0pt5t2pm,lssol,lat,lev=levts,yrange=[-90.,90],/ystyle,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,xtitle='Ls',ytitle='Latitude',title='Zonal Average 1-3 PM 0.5 hPa Temperature',$
     charsize=1.5,c_lab=levts*0.+1
contour,z0pt5t2am,lssol,lat,lev=levts,yrange=[-90.,90],/ystyle,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,xtitle='Ls',ytitle='Latitude',title='Zonal Average 1-3 AM 0.5 hPa Temperature',$
     charsize=1.5,c_lab=levts*0.+1
device,/close

levu=findgen(40)*10.-200.
levv=findgen(50)*2.-50.
levmsf=findgen(50)*4.-100.

pxval = [ lat(0), lat, lat(nlat-2)]

if (min(zpu(*,*,0) gt 1000.)) then goto, skip0

!p.multi=[0,2,3]
device,file=path+'zavgXsect_ls0.eps',/inch,xsize=7.5,ysize=9.5,/encaps,/helvet
contour,zpt(*,*,0),lat,plevels,lev=levts,yrange=[10.,0.01],/ystyle,$
        /ylog,xrange=[-90,90],/xstyle,max_val=1000.,$
        c_lab=levts*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,xtitle='Latitude',ytitle='Pressure (hPa)',$
     title='Zonal Average Temperature at Ls 0',$
     charsize=1.5
pyval = [ 10.0, zpsrf(*,0), 10.0]
polyfill, pxval, pyval, color=40
oplot, lat, zpsrf(*,0), thick=6
contour,zpu(*,*,0),lat,plevels,lev=levu,yrange=[10.,0.01],/ystyle,$
        /ylog,xrange=[-90,90],/xstyle,c_linestyle=(levu lt 0.),max_val=1000.,$
        c_lab=levu*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,xtitle='Latitude',ytitle='Pressure (hPa)',$
     title='Zonal Average U Wind at Ls 0',$
     charsize=1.5
pyval = [ 10., zpsrf(*,0), 10.]
polyfill, pxval, pyval, color=40
oplot, lat, zpsrf(*,0), thick=6
contour,zpv(*,*,0),lat,plevels,lev=levv,yrange=[10.,0.01],/ystyle,$
        /ylog,xrange=[-90,90],/xstyle,c_linestyle=(levv lt 0.),max_val=1000.,$
        c_lab=levv*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,xtitle='Latitude',ytitle='Pressure (hPa)',$
     title='Zonal Average V Wind at Ls 0',$
     charsize=1.5
pyval = [ 10., zpsrf(*,0), 10.]
polyfill, pxval, pyval, color=40
oplot, lat, zpsrf(*,0), thick=6
contour,-1.*zpmsf(*,*,0),lat,plevels,lev=levmsf,yrange=[10.,0.01],/ystyle,$
        /ylog,xrange=[-90,90],/xstyle,c_linestyle=(levmsf lt 0.),max_val=1000.,$
        c_lab=levmsf*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,xtitle='Latitude',ytitle='Pressure (hPa)',$
     title='Mass Stream Function at Ls 0',$
     charsize=1.5
pyval = [ 10., zpsrf(*,0), 10.]
polyfill, pxval, pyval, color=40
oplot, lat, zpsrf(*,0), thick=6
contour,zpt2pm(*,*,0),lat,plevels,lev=levts,yrange=[10.,0.01],/ystyle,$
        /ylog,xrange=[-90,90],/xstyle,max_val=1000.,$
        c_lab=levts*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,xtitle='Latitude',ytitle='Pressure (hPa)',$
     title='Zonal Average 1-3 PM Temperature at Ls 0',$
     charsize=1.5
pyval = [ 10., zpsrf(*,0), 10.]
polyfill, pxval, pyval, color=40
oplot, lat, zpsrf(*,0), thick=6
contour,zpt2am(*,*,0),lat,plevels,lev=levts,yrange=[10.,0.01],/ystyle,$
        /ylog,xrange=[-90,90],/xstyle,max_val=1000.,$
        c_lab=levts*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,xtitle='Latitude',ytitle='Pressure (hPa)',$
     title='Zonal Average 1-3 AM Temperature at Ls 0',$
     charsize=1.5
pyval = [ 10., zpsrf(*,0), 10.]
polyfill, pxval, pyval, color=40
oplot, lat, zpsrf(*,0), thick=6
device,/close

skip0:

if (min(zpu(*,*,3) gt 1000.)) then goto, skip90
!p.multi=[0,2,3]
device,file=path+'zavgXsect_ls90.eps',/inch,xsize=7.5,ysize=9.5,/encaps,/helvet
contour,zpt(*,*,3),lat,plevels,lev=levts,yrange=[10.,0.01],/ystyle,$
        /ylog,xrange=[-90,90],/xstyle,max_val=1000.,$
        c_lab=levts*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,xtitle='Latitude',ytitle='Pressure (hPa)',$
     title='Zonal Average Temperature at Ls 90',$
     charsize=1.5
pyval = [ 10.0, zpsrf(*,3), 10.0]
polyfill, pxval, pyval, color=40
oplot, lat, zpsrf(*,3), thick=6
contour,zpu(*,*,3),lat,plevels,lev=levu,yrange=[10.,0.01],/ystyle,$
        /ylog,xrange=[-90,90],/xstyle,c_linestyle=(levu lt 0.),max_val=1000.,$
        c_lab=levu*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,xtitle='Latitude',ytitle='Pressure (hPa)',$
     title='Zonal Average U Wind at Ls 90',$
     charsize=1.5
pyval = [ 10.0, zpsrf(*,3), 10.0]
polyfill, pxval, pyval, color=40
oplot, lat, zpsrf(*,3), thick=6
contour,zpv(*,*,3),lat,plevels,lev=levv,yrange=[10.,0.01],/ystyle,$
        /ylog,xrange=[-90,90],/xstyle,c_linestyle=(levv lt 0.),max_val=1000.,$
        c_lab=levv*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,xtitle='Latitude',ytitle='Pressure (hPa)',$
     title='Zonal Average V Wind at Ls 90',$
     charsize=1.5
pyval = [ 10.0, zpsrf(*,3), 10.0]
polyfill, pxval, pyval, color=40
oplot, lat, zpsrf(*,3), thick=6
contour,-1.*zpmsf(*,*,3),lat,plevels,lev=levmsf,yrange=[10.,0.01],/ystyle,$
        /ylog,xrange=[-90,90],/xstyle,c_linestyle=(levmsf lt 0.),max_val=1000.,$
        c_lab=levmsf*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,xtitle='Latitude',ytitle='Pressure (hPa)',$
     title='Mass Stream Function at Ls 90',$
     charsize=1.5
pyval = [ 10.0, zpsrf(*,3), 10.0]
polyfill, pxval, pyval, color=40
oplot, lat, zpsrf(*,3), thick=6
contour,zpt2pm(*,*,3),lat,plevels,lev=levts,yrange=[10.,0.01],/ystyle,$
        /ylog,xrange=[-90,90],/xstyle,max_val=1000.,$
        c_lab=levts*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,xtitle='Latitude',ytitle='Pressure (hPa)',$
     title='Zonal Average 1-3 PM Temperature at Ls 90',$
     charsize=1.5
pyval = [ 10.0, zpsrf(*,3), 10.0]
polyfill, pxval, pyval, color=40
oplot, lat, zpsrf(*,3), thick=6
contour,zpt2am(*,*,3),lat,plevels,lev=levts,yrange=[10.,0.01],/ystyle,$
        /ylog,xrange=[-90,90],/xstyle,max_val=1000.,$
        c_lab=levts*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,xtitle='Latitude',ytitle='Pressure (hPa)',$
     title='Zonal Average 1-3 AM Temperature at Ls 90',$
     charsize=1.5
pyval = [ 10.0, zpsrf(*,3), 10.0]
polyfill, pxval, pyval, color=40
oplot, lat, zpsrf(*,3), thick=6
device,/close

skip90:

if (min(zpu(*,*,6) gt 1000.)) then goto, skip180
!p.multi=[0,2,3]
device,file=path+'zavgXsect_ls180.eps',/inch,xsize=7.5,ysize=9.5,/encaps,/helvet
contour,zpt(*,*,6),lat,plevels,lev=levts,yrange=[10.,0.01],/ystyle,$
        /ylog,xrange=[-90,90],/xstyle,max_val=1000.,$
        c_lab=levts*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,xtitle='Latitude',ytitle='Pressure (hPa)',$
     title='Zonal Average Temperature at Ls 180',$
     charsize=1.5
pyval = [ 10.0, zpsrf(*,6), 10.0]
polyfill, pxval, pyval, color=40
oplot, lat, zpsrf(*,6), thick=6
contour,zpu(*,*,6),lat,plevels,lev=levu,yrange=[10.,0.01],/ystyle,$
        /ylog,xrange=[-90,90],/xstyle,c_linestyle=(levu lt 0.),max_val=1000.,$
        c_lab=levu*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,xtitle='Latitude',ytitle='Pressure (hPa)',$
     title='Zonal Average U Wind at Ls 180',$
     charsize=1.5
pyval = [ 10.0, zpsrf(*,6), 10.0]
polyfill, pxval, pyval, color=40
oplot, lat, zpsrf(*,6), thick=6
contour,zpv(*,*,6),lat,plevels,lev=levv,yrange=[10.,0.01],/ystyle,$
        /ylog,xrange=[-90,90],/xstyle,c_linestyle=(levv lt 0.),max_val=1000.,$
        c_lab=levv*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,xtitle='Latitude',ytitle='Pressure (hPa)',$
     title='Zonal Average U Wind at Ls 180',$
     charsize=1.5
pyval = [ 10.0, zpsrf(*,6), 10.0]
polyfill, pxval, pyval, color=40
oplot, lat, zpsrf(*,6), thick=6
contour,-1.*zpmsf(*,*,6),lat,plevels,lev=levmsf,yrange=[10.,0.01],/ystyle,$
        /ylog,xrange=[-90,90],/xstyle,c_linestyle=(levmsf lt 0.),max_val=1000.,$
        c_lab=levmsf*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,xtitle='Latitude',ytitle='Pressure (hPa)',$
     title='Mass Stream Function at Ls 180',$
     charsize=1.5
pyval = [ 10.0, zpsrf(*,6), 10.0]
polyfill, pxval, pyval, color=40
oplot, lat, zpsrf(*,6), thick=6
contour,zpt2pm(*,*,6),lat,plevels,lev=levts,yrange=[10.,0.01],/ystyle,$
        /ylog,xrange=[-90,90],/xstyle,max_val=1000.,$
        c_lab=levts*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,xtitle='Latitude',ytitle='Pressure (hPa)',$
     title='Zonal Average 1-3 PM Temperature at Ls 180',$
     charsize=1.5
pyval = [ 10.0, zpsrf(*,6), 10.0]
polyfill, pxval, pyval, color=40
oplot, lat, zpsrf(*,6), thick=6
contour,zpt2am(*,*,6),lat,plevels,lev=levts,yrange=[10.,0.01],/ystyle,$
        /ylog,xrange=[-90,90],/xstyle,max_val=1000.,$
        c_lab=levts*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,xtitle='Latitude',ytitle='Pressure (hPa)',$
     title='Zonal Average 1-3 AM Temperature at Ls 180',$
     charsize=1.5
pyval = [ 10.0, zpsrf(*,6), 10.0]
polyfill, pxval, pyval, color=40
oplot, lat, zpsrf(*,6), thick=6
device,/close

skip180:

if (min(zpu(*,*,9) gt 1000.)) then goto, skip270
!p.multi=[0,2,3]
device,file=path+'zavgXsect_ls270.eps',/inch,xsize=7.5,ysize=9.5,/encaps,/helvet
contour,zpt(*,*,9),lat,plevels,lev=levts,yrange=[10.,0.01],/ystyle,$
        /ylog,xrange=[-90,90],/xstyle,max_val=1000.,$
        c_lab=levts*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,xtitle='Latitude',ytitle='Pressure (hPa)',$
     title='Zonal Average Temperature at Ls 270',$
     charsize=1.5
pyval = [ 10.0, zpsrf(*,9), 10.0]
polyfill, pxval, pyval, color=40
oplot, lat, zpsrf(*,9), thick=6
contour,zpu(*,*,9),lat,plevels,lev=levu,yrange=[10.,0.01],/ystyle,$
        /ylog,xrange=[-90,90],/xstyle,c_linestyle=(levu lt 0.),max_val=1000.,$
        c_lab=levu*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,xtitle='Latitude',ytitle='Pressure (hPa)',$
     title='Zonal Average U Wind at Ls 270',$
     charsize=1.5
pyval = [ 10.0, zpsrf(*,9), 10.0]
polyfill, pxval, pyval, color=40
oplot, lat, zpsrf(*,9), thick=6
contour,zpv(*,*,9),lat,plevels,lev=levv,yrange=[10.,0.01],/ystyle,$
        /ylog,xrange=[-90,90],/xstyle,c_linestyle=(levv lt 0.),max_val=1000.,$
        c_lab=levv*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,xtitle='Latitude',ytitle='Pressure (hPa)',$
     title='Zonal Average V Wind at Ls 270',$
     charsize=1.5
pyval = [ 10.0, zpsrf(*,9), 10.0]
polyfill, pxval, pyval, color=40
oplot, lat, zpsrf(*,9), thick=6
contour,-1.*zpmsf(*,*,9),lat,plevels,lev=levmsf,yrange=[10.,0.01],/ystyle,$
        /ylog,xrange=[-90,90],/xstyle,c_linestyle=(levmsf lt 0.),max_val=1000.,$
        c_lab=levmsf*0.+1,xtitle='Latitude',ytitle='Pressure (hPa)',$
     title='Mass Stream Function at Ls 270',$
     charsize=1.5
pyval = [ 10.0, zpsrf(*,9), 10.0]
polyfill, pxval, pyval, color=40
oplot, lat, zpsrf(*,9), thick=6
contour,zpt2pm(*,*,9),lat,plevels,lev=levts,yrange=[10.,0.01],/ystyle,$
        /ylog,xrange=[-90,90],/xstyle,max_val=1000.,$
        c_lab=levts*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,xtitle='Latitude',ytitle='Pressure (hPa)',$
     title='Zonal Average 1-3 PM Temperature at Ls 270',$
     charsize=1.5
pyval = [ 10.0, zpsrf(*,9), 10.0]
polyfill, pxval, pyval, color=40
oplot, lat, zpsrf(*,9), thick=6
contour,zpt2am(*,*,9),lat,plevels,lev=levts,yrange=[10.,0.01],/ystyle,$
        /ylog,xrange=[-90,90],/xstyle,max_val=1000.,$
        c_lab=levts*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,xtitle='Latitude',ytitle='Pressure (hPa)',$
     title='Zonal Average 1-3 AM Temperature at Ls 270',$
     charsize=1.5
pyval = [ 10.0, zpsrf(*,9), 10.0]
polyfill, pxval, pyval, color=40
oplot, lat, zpsrf(*,9), thick=6
device,/close

skip270:

levco2i=findgen(20)*100.
!p.multi=[0,1,3]
device,file=path+'wssrf_ls0.eps',/inch,xsize=7.5,ysize=9.5,/encaps,/helvet
contour,ws(*,*,0),lon,lat,xrange=[-180,180],/xstyle,yrange=[-90,90],$
        /ystyle,lev=levws,c_lab=levws*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,ytitle='Latitude',xtitle='E Longitude',$
     title='Surface Wind at Ls 0',$
     charsize=1.5
velovect,uwnd(*,*,0),vwnd(*,*,0),lon,lat,/overplot
contour,mstr(*,*,0)*1000.,lon,lat,xrange=[-180,180],/xstyle,yrange=[-90,90],$
        /ystyle,lev=levss,c_lab=levss*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,ytitle='Latitude',xtitle='E Longitude',$
     title='Surface Stress at Ls 0',$
     charsize=1.5
velovect,ustr(*,*,0),vstr(*,*,0),lon,lat,/overplot
contour,co2i(*,*,0),lon,lat,xrange=[-180,180],/xstyle,yrange=[-90,90],$
        /ystyle,lev=levco2i,c_lab=levco2i*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,ytitle='Latitude',xtitle='E Longitude',$
     title='Surface CO2 Ice at Ls 0',$
     charsize=1.5
device,/close

!p.multi=[0,1,3]
device,file=path+'tsrf_ls0.eps',/inch,xsize=7.5,ysize=9.5,/encaps,/helvet
contour,tsrf(*,*,0),lon,lat,xrange=[-180,180],/xstyle,yrange=[-90,90],$
        /ystyle,lev=levts,c_lab=levts*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,ytitle='Latitude',xtitle='E Longitude',$
     title='Surface Temperature at Ls 0',$
     charsize=1.5
contour,tsrf2pm(*,*,0),lon,lat,xrange=[-180,180],/xstyle,yrange=[-90,90],$
        /ystyle,lev=levts,c_lab=levts*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,ytitle='Latitude',xtitle='E Longitude',$
     title='1-3 PM Surface Temperature at Ls 0',$
     charsize=1.5
contour,tsrf2am(*,*,0),lon,lat,xrange=[-180,180],/xstyle,yrange=[-90,90],$
        /ystyle,lev=levts,c_lab=levts*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,ytitle='Latitude',xtitle='E Longitude',$
     title='1-3 AM Surface Temperature at Ls 0',$
     charsize=1.5
device,/close

!p.multi=[0,1,3]
device,file=path+'t0pt5mb_ls0.eps',/inch,xsize=7.5,ysize=9.5,/encaps,/helvet
contour,t0pt5(*,*,0),lon,lat,xrange=[-180,180],/xstyle,yrange=[-90,90],$
        /ystyle,lev=levts,c_lab=levts*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,ytitle='Latitude',xtitle='E Longitude',$
     title='0.5 hPa Temperature at Ls 0',$
     charsize=1.5
contour,t0pt52pm(*,*,0),lon,lat,xrange=[-180,180],/xstyle,yrange=[-90,90],$
        /ystyle,lev=levts,c_lab=levts*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,ytitle='Latitude',xtitle='E Longitude',$
     title='1-3 PM 0.5 hPa Temperature at Ls 0',$
     charsize=1.5
contour,t0pt52am(*,*,0),lon,lat,xrange=[-180,180],/xstyle,yrange=[-90,90],$
        /ystyle,lev=levts,c_lab=levts*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,ytitle='Latitude',xtitle='E Longitude',$
     title='1-3 AM 0.5 hPa Temperature at Ls 0',$
     charsize=1.5
device,/close

!p.multi=[0,1,3]
device,file=path+'wssrf_ls90.eps',/inch,xsize=7.5,ysize=9.5,/encaps,/helvet
contour,ws(*,*,3),lon,lat,xrange=[-180,180],/xstyle,yrange=[-90,90],$
        /ystyle,lev=levws,c_lab=levws*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,ytitle='Latitude',xtitle='E Longitude',$
     title='Surface Wind at Ls 90',$
     charsize=1.5
velovect,uwnd(*,*,3),vwnd(*,*,3),lon,lat,/overplot
contour,mstr(*,*,3)*1000.,lon,lat,xrange=[-180,180],/xstyle,yrange=[-90,90],$
        /ystyle,lev=levss,c_lab=levss*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,ytitle='Latitude',xtitle='E Longitude',$
     title='Surface Stress at Ls 90',$
     charsize=1.5
velovect,ustr(*,*,3),vstr(*,*,3),lon,lat,/overplot
contour,co2i(*,*,3),lon,lat,xrange=[-180,180],/xstyle,yrange=[-90,90],$
        /ystyle,lev=levco2i,c_lab=levco2i*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,ytitle='Latitude',xtitle='E Longitude',$
     title='Surface CO2 Ice at Ls 90',$
     charsize=1.5
device,/close

!p.multi=[0,1,3]
device,file=path+'tsrf_ls90.eps',/inch,xsize=7.5,ysize=9.5,/encaps,/helvet
contour,tsrf(*,*,3),lon,lat,xrange=[-180,180],/xstyle,yrange=[-90,90],$
        /ystyle,lev=levts,c_lab=levts*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,ytitle='Latitude',xtitle='E Longitude',$
     title='Surface Temperature at Ls 90',$
     charsize=1.5
contour,tsrf2pm(*,*,3),lon,lat,xrange=[-180,180],/xstyle,yrange=[-90,90],$
        /ystyle,lev=levts,c_lab=levts*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,ytitle='Latitude',xtitle='E Longitude',$
     title='1-3 PM Surface Temperature at Ls 90',$
     charsize=1.5
contour,tsrf2am(*,*,3),lon,lat,xrange=[-180,180],/xstyle,yrange=[-90,90],$
        /ystyle,lev=levts,c_lab=levts*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,ytitle='Latitude',xtitle='E Longitude',$
     title='1-3 AM Surface Temperature at Ls 90',$
     charsize=1.5
device,/close

!p.multi=[0,1,3]
device,file=path+'t0pt5mb_ls90.eps',/inch,xsize=7.5,ysize=9.5,/encaps,/helvet
contour,t0pt5(*,*,3),lon,lat,xrange=[-180,180],/xstyle,yrange=[-90,90],$
        /ystyle,lev=levts,c_lab=levts*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,ytitle='Latitude',xtitle='E Longitude',$
     title='0.5 hPa Temperature at Ls 90',$
     charsize=1.5
contour,t0pt52pm(*,*,3),lon,lat,xrange=[-180,180],/xstyle,yrange=[-90,90],$
        /ystyle,lev=levts,c_lab=levts*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,ytitle='Latitude',xtitle='E Longitude',$
     title='1-3 PM 0.5 hPa Temperature at Ls 90',$
     charsize=1.5
contour,t0pt52am(*,*,3),lon,lat,xrange=[-180,180],/xstyle,yrange=[-90,90],$
        /ystyle,lev=levts,c_lab=levts*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,ytitle='Latitude',xtitle='E Longitude',$
     title='1-3 AM 0.5 hPa Temperature at Ls 90',$
     charsize=1.5
device,/close

!p.multi=[0,1,3]
device,file=path+'wssrf_ls180.eps',/inch,xsize=7.5,ysize=9.5,/encaps,/helvet
contour,ws(*,*,6),lon,lat,xrange=[-180,180],/xstyle,yrange=[-90,90],$
        /ystyle,lev=levws,c_lab=levws*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,ytitle='Latitude',xtitle='E Longitude',$
     title='Surface Wind at Ls 180',$
     charsize=1.5
velovect,uwnd(*,*,6),vwnd(*,*,3),lon,lat,/overplot
contour,mstr(*,*,6)*1000.,lon,lat,xrange=[-180,180],/xstyle,yrange=[-90,90],$
        /ystyle,lev=levss,c_lab=levss*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,ytitle='Latitude',xtitle='E Longitude',$
     title='Surface Stress at Ls 180',$
     charsize=1.5
velovect,ustr(*,*,6),vstr(*,*,6),lon,lat,/overplot
contour,co2i(*,*,6),lon,lat,xrange=[-180,180],/xstyle,yrange=[-90,90],$
        /ystyle,lev=levco2i,c_lab=levco2i*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,ytitle='Latitude',xtitle='E Longitude',$
     title='Surface CO2 Ice at Ls 180',$
     charsize=1.5
device,/close

!p.multi=[0,1,3]
device,file=path+'tsrf_ls180.eps',/inch,xsize=7.5,ysize=9.5,/encaps,/helvet
contour,tsrf(*,*,6),lon,lat,xrange=[-180,180],/xstyle,yrange=[-90,90],$
        /ystyle,lev=levts,c_lab=levts*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,ytitle='Latitude',xtitle='E Longitude',$
     title='Surface Temperature at Ls 180',$
     charsize=1.5
contour,tsrf2pm(*,*,6),lon,lat,xrange=[-180,180],/xstyle,yrange=[-90,90],$
        /ystyle,lev=levts,c_lab=levts*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,ytitle='Latitude',xtitle='E Longitude',$
     title='1-3 PM Surface Temperature at Ls 180',$
     charsize=1.5
contour,tsrf2am(*,*,6),lon,lat,xrange=[-180,180],/xstyle,yrange=[-90,90],$
        /ystyle,lev=levts,c_lab=levts*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,ytitle='Latitude',xtitle='E Longitude',$
     title='1-3 AM Surface Temperature at Ls 180',$
     charsize=1.5
device,/close

!p.multi=[0,1,3]
device,file=path+'t0pt5mb_ls180.eps',/inch,xsize=7.5,ysize=9.5,/encaps,/helvet
contour,t0pt5(*,*,6),lon,lat,xrange=[-180,180],/xstyle,yrange=[-90,90],$
        /ystyle,lev=levts,c_lab=levts*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,ytitle='Latitude',xtitle='E Longitude',$
     title='0.5 hPa Temperature at Ls 180',$
     charsize=1.5
contour,t0pt52pm(*,*,6),lon,lat,xrange=[-180,180],/xstyle,yrange=[-90,90],$
        /ystyle,lev=levts,c_lab=levts*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,ytitle='Latitude',xtitle='E Longitude',$
     title='1-3 PM 0.5 hPa Temperature at Ls 180',$
     charsize=1.5
contour,t0pt52am(*,*,6),lon,lat,xrange=[-180,180],/xstyle,yrange=[-90,90],$
        /ystyle,lev=levts,c_lab=levts*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,ytitle='Latitude',xtitle='E Longitude',$
     title='1-3 AM 0.5 hPa Temperature at Ls 180',$
     charsize=1.5
device,/close

!p.multi=[0,1,3]
device,file=path+'wssrf_ls270.eps',/inch,xsize=7.5,ysize=9.5,/encaps,/helvet
contour,ws(*,*,9),lon,lat,xrange=[-180,180],/xstyle,yrange=[-90,90],$
        /ystyle,lev=levws,c_lab=levws*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,ytitle='Latitude',xtitle='E Longitude',$
     title='Surface Wind at Ls 270',$
     charsize=1.5
velovect,uwnd(*,*,9),vwnd(*,*,9),lon,lat,/overplot
contour,mstr(*,*,9)*1000.,lon,lat,xrange=[-180,180],/xstyle,yrange=[-90,90],$
        /ystyle,lev=levss,c_lab=levss*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,ytitle='Latitude',xtitle='E Longitude',$
     title='Surface Stress at Ls 270',$
     charsize=1.5
velovect,ustr(*,*,9),vstr(*,*,0),lon,lat,/overplot
contour,co2i(*,*,9),lon,lat,xrange=[-180,180],/xstyle,yrange=[-90,90],$
        /ystyle,lev=levco2i,c_lab=levco2i*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,ytitle='Latitude',xtitle='E Longitude',$
     title='Surface CO2 Ice at Ls 270',$
     charsize=1.5
device,/close

!p.multi=[0,1,3]
device,file=path+'tsrf_ls270.eps',/inch,xsize=7.5,ysize=9.5,/encaps,/helvet
contour,tsrf(*,*,9),lon,lat,xrange=[-180,180],/xstyle,yrange=[-90,90],$
        /ystyle,lev=levts,c_lab=levts*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,ytitle='Latitude',xtitle='E Longitude',$
     title='Surface Temperature at Ls 270',$
     charsize=1.5
contour,tsrf2pm(*,*,9),lon,lat,xrange=[-180,180],/xstyle,yrange=[-90,90],$
        /ystyle,lev=levts,c_lab=levts*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,ytitle='Latitude',xtitle='E Longitude',$
     title='1-3 PM Surface Temperature at Ls 270',$
     charsize=1.5
contour,tsrf2am(*,*,9),lon,lat,xrange=[-180,180],/xstyle,yrange=[-90,90],$
        /ystyle,lev=levts,c_lab=levts*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,ytitle='Latitude',xtitle='E Longitude',$
     title='1-3 AM Surface Temperature at Ls 270',$
     charsize=1.5
device,/close

!p.multi=[0,1,3]
device,file=path+'t0pt5mb_ls270.eps',/inch,xsize=7.5,ysize=9.5,/encaps,/helvet
contour,t0pt5(*,*,9),lon,lat,xrange=[-180,180],/xstyle,yrange=[-90,90],$
        /ystyle,lev=levts,c_lab=levts*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,ytitle='Latitude',xtitle='E Longitude',$
     title='0.5 hPa Temperature at Ls 270',$
     charsize=1.5
contour,t0pt52pm(*,*,9),lon,lat,xrange=[-180,180],/xstyle,yrange=[-90,90],$
        /ystyle,lev=levts,c_lab=levts*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,ytitle='Latitude',xtitle='E Longitude',$
     title='1-3 PM 0.5 hPa Temperature at Ls 270',$
     charsize=1.5
contour,t0pt52am(*,*,9),lon,lat,xrange=[-180,180],/xstyle,yrange=[-90,90],$
        /ystyle,lev=levts,c_lab=levts*0.+1,xthick=2.5,ythick=2.5,thick=2.5,$
     charthick=2.5,ytitle='Latitude',xtitle='E Longitude',$
     title='1-3 AM 0.5 hPa Temperature at Ls 270',$
     charsize=1.5
device,/close


set_plot,'x'


end
