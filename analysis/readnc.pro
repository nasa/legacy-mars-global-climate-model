pro readnc,lssol,gasp,gagt,gd,gwv,gwi,hco2_n,hco2_s,vl1,vl2,lat,zd,zwv,$
zwi,zws,zss,zt,zt2pm,zt2am,zco2i,z0pt5t,z0pt5t2pm,z0pt5t2am,lon,lsseas,$
nlsbin,uwnd,vwnd,ws,ustr,vstr,mstr,co2i,tsrf,tsrf2pm,tsrf2am,t0pt5,$
t0pt52pm,t0pt52am,plevels,zpt,zpt2pm,zpt2am,zpu,zpv,zpmsf,zpsrf,$
nlat,nlon,nsolcnt,nseas,npres

read_ncdf,'GCM_data.nc',data

lat=data.lat
lon=data.lon
npres=data.npres
lssol=data.lssol
gasp=data.gasp
gagt=data.gagt
gd=data.gd
gwv=data.gd
gwi=data.gwi
hco2_n=data.hco2_n
hco2_s=data.hco2_s
vl1=data.vl1
vl2=data.vl2
zd=data.zd
zwv=data.zwv
zwi=data.zwi
zws=data.zws
zss=data.zss
zt=data.zt
zt2pm=data.zt2pm
zt2am=data.zt2am
zco2i=data.zco2i
z0pt5t=data.z0pt5t
z0pt5t2pm=data.z0pt5t2pm
z0pt5t2am=data.z0pt5t2am
zpt=data.zpt
zpt2pm=data.zpt2pm
zpt2am=data.zpt2am
zpu=data.zpu
zpv=data.zpv
zpmsf=data.zpmsf
;plevels=data.pinterp
uwnd=data.uwnd
vwnd=data.vwnd
ws=data.ws
ustr=data.ustr
vstr=data.vstr
mstr=data.mstr
co2i=data.co2i
tsrf=data.tsrf
tsrf2pm=data.tsrf2pm
tsrf2am=data.tsrf2am
t0pt5=data.t0pt5
t0pt52pm=data.t0pt52pm
t0pt52am=data.t0pt52am
zpsrf=data.zpsrf
plevels=data.pinterp

nsolcnt=n_elements(lssol)
nlat=n_elements(lat)
nlon=n_elements(lon)
npres=n_elements(plevels)

stop

return
end
