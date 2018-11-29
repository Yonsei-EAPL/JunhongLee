;from
;http://www2.mmm.ucar.edu/wrf/OnLineTutorial/Graphics/NCL/Examples/BASIC/wrf_wps_ter1.ncl

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
begin

;##############################################################################

Varn = "UV"
xdim = 174
ydim = 126

;##############################################################################

file1 = systemfunc("ls ../chem_ctl_v1.8.3/wrfout_d01_2016-12-*:00:00")
file2 = systemfunc("ls ../chem_rsl_v1.9.2/wrfout_d01_2016-12-*:00:00")
file3 = systemfunc("ls ../chem_rsl_v1.9.2/wrfout_d01_2016-12-01_01:00:00")
c = addfile(file3,"r")
msk = c->XLAND

pvar11 = new((/241,1,30,ydim,xdim/),"float")
pvar12 = new((/241,1,30,ydim,xdim/),"float")
pvar21 = new((/241,1,30,ydim,xdim/),"float")
pvar22 = new((/241,1,30,ydim,xdim/),"float")


do i = 1, 97
;do i = 1, 240
a = addfile(file1(i-1),"r")
b = addfile(file2(i-1),"r")
pvar11(i-1,0,:,:,:) = wrf_user_getvar(a,"ua",0)
pvar12(i-1,0,:,:,:) = wrf_user_getvar(a,"va",0)
pvar21(i-1,0,:,:,:) = wrf_user_getvar(b,"ua",0)
pvar22(i-1,0,:,:,:) = wrf_user_getvar(b,"va",0)
end do

count = 0
do i = 1, 97
;do i = 1, 240
  mods = mod(i-1,24)
    if ( count .eq. 0 ) then
      var1 = sqrt(pvar11(i-1,0,0,:,:)^2+pvar12(i-1,0,0,:,:)^2)
      var2 = sqrt(pvar21(i-1,0,0,:,:)^2+pvar22(i-1,0,0,:,:)^2)
    else
      var1(:,:) = var1(:,:) + sqrt(pvar11(i-1,0,0,:,:)^2+pvar12(i-1,0,0,:,:)^2)
      var2(:,:) = var2(:,:) + sqrt(pvar21(i-1,0,0,:,:)^2+pvar22(i-1,0,0,:,:)^2)
    end if
    count = count + 1
end do
do i = 0, ydim-1
do j = 0, xdim-1
;  var1(i,j) = mask(var1(i,j)/count,msk(0,i,j),1)
;  var2(i,j) = mask(var2(i,j)/count,msk(0,i,j),1)
  var1(i,j) = var1(i,j)/count
  var2(i,j) = var2(i,j)/count
end do
end do

type = "x11"
wks = gsn_open_wks("x11","win_dif")
wks1 = gsn_open_wks("png","win_dif")

;printVarSummary(var1)
;print(var1(:,:))

res = True
res@cnFillOn = True
res@cnFillMode = "RasterFill"
res@cnLinesOn = False

res@cnLineDrawOrder = "Predraw"
res@cnFillDrawOrder = "Predraw"

res@cnFillPalette = "MPL_bwr"
res@cnLevelSelectionMode = "ExplicitLevels"
res@cnLevels = fspan(-3.0,3.0,101)
res@cnFillColors = ispan(13,114,1)
print(ispan(15,116,1))

res@cnInfoLabelOn = False
res@lbLabelBarOn = True
res@lbBoxLinesOn = False
res@lbLabelStride = 25
res@Footer = False
res@NoHeaderFooter = False
res@InitTime = False
res@tmBorderThicknessF = 3.0

res@lbLabelOffsetF = 0.045
res@lbTitleOn = False


;plot = gsn_contour(wks,var2(:,:)-var1(:,:),True)
;plot = gsn_csm_contour(wks,var2(:,:)-var1(:,:),res)
plot = wrf_contour(c,wks,var2(:,:)-var1(:,:),res)
plot1 = wrf_contour(c,wks1,var2(:,:)-var1(:,:),res)
alot = wrf_contour(c,wks,var2(20:120,35:125)-var1(20:120,35:125),res)
alot1 = wrf_contour(c,wks1,var2(20:120,35:125)-var1(20:120,35:125),res)


pltres = True
pltres@NoTitles = True
pltres@CommonTitle = False
;pltres@FramePlot = False
pltres@PanelPlat = False
;pltres@InitTime = False
;pltres@ValidTime = False
;Pltres@Footer = False
;Pltres@NoHeaderFooter = False

mpres = True
mpres@mpGeophysicalLineColor = "gray75"
mpres@mpGeophysicalLineThicknessF = 2
mpres@ZoomIn = True        ; Tell wrf_map_resources we want to zoom in.
mpres@Xstart = 35     ; Set these four special WRF resources
mpres@Xend   = 125       ; required for zooming.
mpres@Ystart = 20
mpres@Yend   = 120
plot2 = wrf_map_overlays(c,wks,(/alot/),pltres,mpres)
plot22 = wrf_map_overlays(c,wks1,(/alot1/),pltres,mpres)

end

