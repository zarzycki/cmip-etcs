load "~/Software/cymep/ncl/functions/getTrajectoriesGeneric.ncl"
load "~/Software/cymep/ncl/functions/trackDensity.ncl"

thefile="etc-all-traj.txt"
ncols = -1 ; number of columns in non-header rows of trajectory file
maxLines = -1  ; set to -1 for auto-detection
traj_data = getTrajectoriesGeneric(thefile,ncols,maxLines,"start","tab",False,False)

; load in trajectory data (need to be modified if using different ASCII output)
; USER_MODIFY
xlon  = traj_data(:,2,:)
xlat  = traj_data(:,3,:)
xpres  = tofloat(traj_data(:,4,:))/100.
xyear  = toint(traj_data(:,5,:))
xmonth  = toint(traj_data(:,6,:))
xday  = toint(traj_data(:,7,:))
xhour  = toint(traj_data(:,8,:))

nyears=16.
clat=ndtooned(xlat)
clon=ndtooned(xlon)
gridsize=5.0
trackdens = track_density(gridsize,0.0,clat,clon,True)
trackdens = trackdens/nyears

trackdens = where(trackdens.le.1.0,trackdens@_FillValue,trackdens)


wks = gsn_open_wks("pdf","fig-spatial-"+thefile)               ; open workstation

  res                       = True     ; plot mods desired
  res@gsnDraw = False
  res@gsnFrame = False
  res@gsnAddCyclic          = True    
  res@cnFillOn              = True     ; turn on color fill
  res@cnFillMode            = "RasterFill"       ; Raster Mode
  res@cnLinesOn             = False    ; turn of contour lines
  res@cnLineLabelsOn  = False

if (thefile .eq. "etc-strong-traj.txt") then
  res@gsnCenterString       = "Annual ETC track density (<=990hPa)"
else
  res@gsnCenterString       = "Annual ETC track density (all)"
end if

  ncontlev=10
  res@cnLevelSelectionMode  = "ManualLevels"
  res@cnMinLevelValF        = 0
  res@cnMaxLevelValF        = 25
  res@cnLevelSpacingF       = (res@cnMaxLevelValF - res@cnMinLevelValF) / tofloat(ncontlev)
  res@cnFillPalette="precip_11lev"

  res@mpProjection      = "Orthographic"
  res@mpOutlineBoundarySets = "National"
  res@mpCenterLatF      =  35.
  res@mpCenterLonF      =  290.
  res@mpLandFillColor        = "tan"
  res@mpOceanFillColor       = "slategray3"
  res@mpInlandWaterFillColor = "slategray3"
  res@mpGeophysicalLineThicknessF = 2.0
  
  plot = gsn_csm_contour_map_ce(wks,trackdens,res)

  

;---Create some text strings to attach outside of plot.
  text_str = "6-hr ETC fix per 5~S~o~N~x5~S~o~N~ box per yr."

  txres               = True
  txres@txFontHeightF = 0.011
  txres@txAngleF      = 0     ; Rotate 90 degrees

  txid1 = gsn_create_text(wks,text_str,txres)

;----------------------------------------------------------------------
; Use gsn_add_annotation to attach text string to outside of plot. 
; This is useful if you need to resize the plot later, like in
; a panel plot.
;----------------------------------------------------------------------
  amres                  = True
  amres@amJust           = "CenterCenter"  ; Center string
  amres@amParallelPosF   =  0.0          ; -0.5 - 0 - 0.5 from L -> R
  amres@amOrthogonalPosF =  0.60          ; Move to center (in Y direction)
  annoid1 = gsn_add_annotation(plot, txid1, amres)

;drawNDCGrid(wks)

draw(plot)
frame(wks)


