filename = "composite_CESM1-LENS.nc"
varname="IVT"   ; T850, IVT, PRECT
f     = addfile(filename,"r")

if (varname .eq. "IVT") then
  tmpq = f->Q
  tmpu = f->U
  tmpv = f->V
  tmpps = f->PS
  P0=100000.
  tmpdp = dpres_hybrid_ccm(tmpps,P0,f->hyai,f->hybi)
  tmpivt = tmpdp*tmpq*tmpu
  ivt_u = dim_sum_n_Wrap(tmpivt,0)
  tmpivt = tmpdp*tmpq*tmpv
  ivt_v = dim_sum_n_Wrap(tmpivt,0)
  var = sqrt(ivt_u^2. + ivt_v^2.)
  copy_VarMeta(tmpps,var)
else if (varname .eq. "T850") then
  var = f->T({850},:,:)
else
  var = f->$varname$
end if
end if



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

if (varname .eq. "T850") then
  U = f->U({850},:,:)
  V = f->V({850},:,:)
  U = where(gcdist .lt. 24.0,U,U@_FillValue)
  V = where(gcdist .lt. 24.0,V,V@_FillValue)
end if

if (varname .eq. "IVT") then
  U = f->U(:,:,:)
  V = f->V(:,:,:)
  PS = f->PS(:,:)
  P0=100000.
  U!0="lev"
  U!1="lat"
  U!2="lon"
  copy_VarMeta(U,V)
  copy_VarMeta(var,PS)
  OMEGA = omega_ccm_driver(P0,PS,U,V,f->hyam,f->hybm,f->hyai,f->hybi)
  OMEGA!0="lev"
  OMEGA!1="lat"
  OMEGA!2="lon"
  nlev = dimsizes(f->hybm)
  do ii = 0,nlev-1
    do jj = 0,1
      OMEGA(ii,:,:) = smth9_Wrap(OMEGA(ii,:,:), 0.50,  0.25, False)
    end do
    OMEGA(ii,:,:) = where(gcdist .lt. 21.0,OMEGA(ii,:,:),OMEGA@_FillValue)
  end do
  OMEGA = OMEGA * 3600. / 100.
end if

;******************************************************
; create plot
;******************************************************
wks = gsn_open_wks("pdf","fig-composite-"+varname)               ; open workstation
gsn_define_colormap(wks,"BlAqGrYeOrReVi200")   ; choose colormap

res                     = True         ; plot modifications desired
;  res@gsnMaximize         = True         ; Maximize size of plot in frame

res@gsnDraw             = False           ; don't draw
res@gsnFrame            = False           ; don't advance frame

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
  markerColor="black"
  markerColor2="white"
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
  markerColor="white"
  ncontlev=10
  res@cnLevelSelectionMode  = "ManualLevels"
  res@cnMinLevelValF        = 235
  res@cnMaxLevelValF        = 285
  res@cnLevelSpacingF       = (res@cnMaxLevelValF - res@cnMinLevelValF) / tofloat(ncontlev)
  res@cnFillPalette="MPL_rainbow"
  res@gsnLeftString="2-m temperature"
  res@gsnRightString="K"
end if
if (varname .eq. "T850") then
  markerColor="black"
  markerColor2="white"
  ncontlev=10
  res@cnLevelSelectionMode  = "ManualLevels"
  res@cnMinLevelValF        = 255
  res@cnMaxLevelValF        = 280
  res@cnLevelSpacingF       = (res@cnMaxLevelValF - res@cnMinLevelValF) / tofloat(ncontlev)
  res@cnFillPalette="MPL_rainbow"
  res@gsnLeftString="850 hPa temperature"
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
if (varname .eq. "IVT") then
  markerColor="black"
  markerColor2="green"
  ncontlev=9
  res@cnLevelSelectionMode  = "ManualLevels"
  res@cnMinLevelValF        = 280.
  res@cnMaxLevelValF        = 2800.
  res@cnLevelSpacingF       = (res@cnMaxLevelValF - res@cnMinLevelValF) / tofloat(ncontlev)
  res@gsnLeftString="Integrated Vapor Transport"
  res@gsnRightString="kg/m/s"
  cmap = read_colormap_file("MPL_gnuplot")
  tmpdims = dimsizes(cmap)
  newcmap = new((/tmpdims(0)+1,tmpdims(1)/),typeof(cmap))
  newcmap(1:tmpdims(0),:) =  cmap(::-1,:)
  newcmap(0,:) = (/1,1,1,1/)
  res@cnFillPalette = newcmap
end if

res@mpMinLatF             = 0. - 24.
res@mpMaxLatF             = 0. + 24.
res@mpMinLonF             = 0. - 24.
res@mpMaxLonF             = 0. + 24.

res@lbLabelAutoStride   = True         ; Clean up labelbar labels.

map = gsn_csm_contour_map_ce(wks,var,res)

if (varname .eq. "TREFHT" .or. varname .eq. "T850") then
  vecres                  = True            ; vector only resources
  vecres@gsnDraw          = False           ; don't draw
  vecres@gsnFrame         = False           ; don't advance frame
  ;vecres@vcGlyphStyle     = "CurlyVector"   ; curly vectors
  vecres@vcRefMagnitudeF  = 18              ; define vector ref mag
  vecres@vcRefLengthF     = 0.045           ; define length of vec ref
  vecres@gsnRightString   = " "             ; turn off right string
  vecres@gsnLeftString    = " "             ; turn off left string
  vecres@tiXAxisString    = " "             ; turn off axis label 
  vecres@vcRefAnnoOrthogonalPosF = -.44    ; move ref vector into plot
  vecres@vcLineArrowThicknessF =6.0
  vecstride=4
  mapV   = gsn_csm_vector(wks,U(::vecstride,::vecstride),V(::vecstride,::vecstride),vecres)
  overlay(map,mapV)                      ; result will be plotA
end if

if (varname .eq. "IVT") then
  sres                      = True              ; set up a second resource list
  sres@gsnDraw              = False             ; do not draw the plot
  sres@gsnFrame             = False             ; do not advance the frame
  sres@cnLineLabelsOn  = False
  sres@cnLevelSelectionMode = "ExplicitLevels"  ; use explicit levels
  sres@cnLevels             = (/-2,-3,-4,-5,-6/)   ; set the contour levels
  sres@cnLineThicknessF  = 5.0
  sres@cnLineColor ="white"
  sres@cnInfoLabelOn = False
  sres@gsnLeftString=""
  sres@gsnRightString=""
  plot_ov = gsn_csm_contour(wks,OMEGA({600.},::1,::1),sres)      ; create the U-wind plot
  overlay(map,plot_ov)                      ; result will be plotA
end if

pmres=True
pmres@gsMarkerColor = markerColor
pmres@gsMarkerIndex = 5
pmres@gsMarkerThicknessF = 15.
pmres@gsMarkerSizeF=20.
dum = gsn_add_polymarker(wks, map, 0., 0., pmres)

pmres@gsMarkerColor = markerColor2
pmres@gsMarkerThicknessF = pmres@gsMarkerThicknessF*0.45
;pmres@gsMarkerSizeF=pmres@gsMarkerSizeF*0.6
dum2 = gsn_add_polymarker(wks, map, 0., 0., pmres)

draw(map)
frame(wks)

