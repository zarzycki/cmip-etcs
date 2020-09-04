filename = "composite_CESM1-LENS.nc"
varname="PRECT"
f     = addfile(filename,"r")
var = f->$varname$
var!0="lat"
var!1="lon"

lat2d=conform(var,var&lat,0)
lon2d=conform(var,var&lon,1)
copy_VarMeta(var,lat2d)
copy_VarMeta(var,lon2d)
gcdist = gc_latlon(0.,0.,lat2d,lon2d,2,2)
copy_VarMeta(var,gcdist)

var = where(gcdist .lt. 24.0,var,var@_FillValue)

if (varname .eq. "PRECT") then
  var = 8.64e7*var
end if
if (varname .eq. "QREFHT") then
  var = var*1000.
end if
;******************************************************
; create plot
;******************************************************
wks = gsn_open_wks("pdf","fig-composite-"+varname)               ; open workstation
gsn_define_colormap(wks,"BlAqGrYeOrReVi200")   ; choose colormap

res                     = True         ; plot modifications desired
;  res@gsnMaximize         = True         ; Maximize size of plot in frame

res@cnFillOn            = True         ; Turn on contour fill
;  res@cnFillMode          = "AreaFill"   ; Style of fill. You can also

res@cnLinesOn           = True        ; Turn off contour lines
res@cnLineLabelsOn      = False        ; Turn off contour line labels

res@mpOutlineOn =False
res@mpFillOn = False

res@gsnAddCyclic = False



;res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
;res@cnLevels    = (/ 0.5,1.5,2.5,3.5/)   ; set levels
;res@cnFillColors = (/ "lightblue","paleturquoise","palegreen","yellow","red"/)  ; set the colors to be used

if (varname .eq. "PRECT") then
  ncontlev=9
  res@cnLevelSelectionMode  = "ManualLevels"
  res@cnMinLevelValF        = 2
  res@cnMaxLevelValF        = 20
  res@cnLevelSpacingF       = (res@cnMaxLevelValF - res@cnMinLevelValF) / tofloat(ncontlev)
  res@cnFillPalette="precip3_16lev"
res@gsnLeftString="Precipitation rate"
res@gsnRightString="mm/day"
end if
if (varname .eq. "TREFHT") then
  ncontlev=10
  res@cnLevelSelectionMode  = "ManualLevels"
  res@cnMinLevelValF        = 235
  res@cnMaxLevelValF        = 285
  res@cnLevelSpacingF       = (res@cnMaxLevelValF - res@cnMinLevelValF) / tofloat(ncontlev)
  res@cnFillPalette="MPL_rainbow"
res@gsnLeftString="2-m temperature"
res@gsnRightString="K"
end if
if (varname .eq. "QREFHT") then
  ncontlev=9
  res@cnLevelSelectionMode  = "ManualLevels"
  res@cnMinLevelValF        = 1.
  res@cnMaxLevelValF        = 10
  res@cnLevelSpacingF       = (res@cnMaxLevelValF - res@cnMinLevelValF) / tofloat(ncontlev)
  res@cnFillPalette="MPL_PuBuGn"
  res@gsnLeftString="2-m specific humidity"
  res@gsnRightString="g/kg"
end if

res@mpMinLatF             = 0. - 24.
res@mpMaxLatF             = 0. + 24.
res@mpMinLonF             = 0. - 24.
res@mpMaxLonF             = 0. + 24.

res@lbLabelAutoStride   = True         ; Clean up labelbar labels.

map = gsn_csm_contour_map_ce(wks,var,res)