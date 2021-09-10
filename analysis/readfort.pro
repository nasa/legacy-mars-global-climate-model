pro readfort,lssol,gasp,gagt,gd,gwv,gwi,hco2_n,hco2_s,vl1,vl2,lat,zd,zwv,$
zwi,zws,zss,zt,zt2pm,zt2am,zco2i,z0pt5t,z0pt5t2pm,z0pt5t2am,lon,lsseas,$
nlsbin,uwnd,vwnd,ws,ustr,vstr,mstr,co2i,tsrf,tsrf2pm,tsrf2am,t0pt5,$
t0pt52pm,t0pt52am,plevels,zpt,zpt2pm,zpt2am,zpu,zpv,zpmsf,zpsrf,$
nlat,nlon,nsolcnt,nseas,npres

nsolcnt=intarr(1)
nseas=intarr(1)
nlat=intarr(1)
nlon=intarr(1)
npres=intarr(1)

openr,10,'fort.40',/f77_unformatted
readu,10,nsolcnt

lssol=fltarr(nsolcnt)
gasp=fltarr(nsolcnt)
gagt=fltarr(nsolcnt)
gd=fltarr(nsolcnt)
gwv=fltarr(nsolcnt)
gwi=fltarr(nsolcnt)
hco2_n=fltarr(nsolcnt)
hco2_s=fltarr(nsolcnt)
vl1=fltarr(nsolcnt)
vl2=fltarr(nsolcnt)

readu,10,lssol
readu,10,gasp
readu,10,gagt
readu,10,gd
readu,10,gwv
readu,10,gwi
readu,10,hco2_n
readu,10,hco2_s
readu,10,vl1
readu,10,vl2
close,10

openr,11,'fort.41',/f77_unformatted
readu,11,nsolcnt
readu,11,nlat

lat=fltarr(nlat-1)
zd=fltarr(nsolcnt,nlat-1)
zwv=fltarr(nsolcnt,nlat-1)
zwi=fltarr(nsolcnt,nlat-1)
zws=fltarr(nsolcnt,nlat-1)
zss=fltarr(nsolcnt,nlat-1)
zt=fltarr(nsolcnt,nlat-1)
zt2pm=fltarr(nsolcnt,nlat-1)
zt2am=fltarr(nsolcnt,nlat-1)
zco2i=fltarr(nsolcnt,nlat-1)
z0pt5t=fltarr(nsolcnt,nlat-1)
z0pt5t2pm=fltarr(nsolcnt,nlat-1)
z0pt5t2am=fltarr(nsolcnt,nlat-1)

readu,11,lssol
readu,11,lat
readu,11,zd
readu,11,zwv
readu,11,zwi
readu,11,zws
readu,11,zss
readu,11,zt
readu,11,zt2pm
readu,11,zt2am
readu,11,zco2i
readu,11,z0pt5t
readu,11,z0pt5t2pm
readu,11,z0pt5t2am
close,11

openr,12,'fort.42',/f77_unformatted
readu,12,nseas
readu,12,nlat
readu,12,nlon

lat=fltarr(nlat-1)
lon=fltarr(nlon)
lsseas=fltarr(nseas)
nlsbin=fltarr(nseas)
uwnd=fltarr(nlon,nlat-1,nseas)
vwnd=fltarr(nlon,nlat-1,nseas)
ws=fltarr(nlon,nlat-1,nseas)
ustr=fltarr(nlon,nlat-1,nseas)
vstr=fltarr(nlon,nlat-1,nseas)
mstr=fltarr(nlon,nlat-1,nseas)
co2i=fltarr(nlon,nlat-1,nseas)
tsrf=fltarr(nlon,nlat-1,nseas)
tsrf2pm=fltarr(nlon,nlat-1,nseas)
tsrf2am=fltarr(nlon,nlat-1,nseas)
t0pt5=fltarr(nlon,nlat-1,nseas)
t0pt52pm=fltarr(nlon,nlat-1,nseas)
t0pt52am=fltarr(nlon,nlat-1,nseas)

readu,12,lsseas
readu,12,nlsbin
readu,12,lat
readu,12,lon
readu,12,uwnd
readu,12,vwnd
readu,12,ws
readu,12,ustr
readu,12,vstr
readu,12,mstr
readu,12,co2i
readu,12,tsrf
readu,12,tsrf2pm
readu,12,tsrf2am
readu,12,t0pt5
readu,12,t0pt52pm
readu,12,t0pt52am

close,12

openr,13,'fort.43',/f77_unformatted
readu,13,nseas
readu,13,nlat
readu,13,npres

plevels=fltarr(npres)
zpt=fltarr(nlat-1,npres,nseas)
zpt2pm=fltarr(nlat-1,npres,nseas)
zpt2am=fltarr(nlat-1,npres,nseas)
zpu=fltarr(nlat-1,npres,nseas)
zpv=fltarr(nlat-1,npres,nseas)
zpmsf=fltarr(nlat-1,npres,nseas)
zpsrf=fltarr(nlat-1,nseas)

readu,13,plevels
readu,13,zpt
readu,13,zpt2pm
readu,13,zpt2am
readu,13,zpu
readu,13,zpv
readu,13,zpmsf
readu,13,zpsrf

close,13



return
end
