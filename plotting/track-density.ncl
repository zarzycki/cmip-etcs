;=================== Description ==================================
; This script reads in Tempeset
; tracking code and plots a track density
;==================================================================

load "~/ncl/functions/getTrajectoriesGeneric.ncl"
load "~/ncl//functions/trackDensity.ncl"

begin

;=================== User Settings ================================
gridsize = 8.0
out_type = "png"
nyrs=3.

;======================= Constants ================================
ms_to_kts = 1.94384449

;=================== Get traj data ================================
thefile="./trajectories.txt.JRA"  ; trajectory file full path
ncols = -1 ; number of columns in non-header rows of trajectory file
maxLines = -1  ; set to -1 for auto-detection
traj_data = getTrajectoriesGeneric(thefile,ncols,maxLines,"start","tab",False,False)

xlon  = traj_data(:,2,:)
xlat  = traj_data(:,3,:)
xpres  = traj_data(:,4,:)/100.
xwind  = traj_data(:,5,:)

clat = ndtooned(xlat)
clon = ndtooned(xlon)
cpres = ndtooned(xpres)
cwind = ndtooned(xwind)

count1 = track_density(gridsize,0.0,clat,clon,False)

count1=count1/nyrs

delete([/xlon,xlat,xpres,xwind,clon,clat,cpres,cwind,nyrs,traj_data,maxLines,ncols,thefile/])

;=================== change zeros to missing values ================================

count1 = where(count1.eq.0, count1@_FillValue,count1)


;************************************************
; create plot
plotname = "anntrackdens_all"

wks = gsn_open_wks(out_type,plotname)
plot = new(1,"graphic")

colorMap1="WhiteBlueGreenYellowRed"

res                       = True     ; plot mods desired
res@gsnDraw = False
res@gsnFrame = False
res@gsnAddCyclic          = False    

res@cnFillOn              = True     ; turn on color fill
res@cnFillMode            = "RasterFill"       ; Raster Mode
res@cnLinesOn             = False    ; turn of contour lines
res@cnLineLabelsOn  = False
res@lbLabelStride = 2
res@gsnCenterString       = ""

res@mpMinLatF             = -85.
res@mpMaxLatF             = 85.
res@mpMinLonF             = 0.
res@mpMaxLonF             = 359.
res@mpCenterLonF          = (res@mpMinLonF + res@mpMaxLonF)/2


print("min: "+min(count1({res@mpMinLatF:res@mpMaxLatF},{res@mpMinLonF:res@mpMaxLonF}))+"   max: "+max(count1({res@mpMinLatF:res@mpMaxLatF},{res@mpMinLonF:res@mpMaxLonF})))


; plot absolute data
res@cnLevelSelectionMode  = "ManualLevels"
; 0-120
res@cnMinLevelValF        = 0.0
res@cnMaxLevelValF        = max(count1({res@mpMinLatF:res@mpMaxLatF},{res@mpMinLonF:res@mpMaxLonF}))
res@cnLevelSpacingF       = res@cnMaxLevelValF/20.
res@cnFillPalette=colorMap1



plot = gsn_csm_contour_map(wks,count1(:,:),res)

draw(plot)
frame(wks)

end
