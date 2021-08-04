;========================================================
;
;========================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;========================================================
begin

write_output=True
isglobal=True
isStoN=False
keepStoN=False

;  doLonShift=True
doLonShift=False ; True for MERRA, False for JRA

topoDir = getenv("TOPOORIGDIR")
print("TopoDir is: "+topoDir)
topoFileName=topoDir+"/2deg_cesm_topo_latlon.nc"
;outFileName = whatData+".topo.nc"

;trackerFileName="/glade/collections/cmip/CMIP6/CMIP/MRI/MRI-ESM2-0/historical/r1i1p1f1/6hrPlevPt/psl/gn/v20190722/psl/psl_6hrPlevPt_MRI-ESM2-0_historical_r1i1p1f1_gn_195001010000-195912311800.nc"
;trackerFileName="/glade/collections/cmip/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/historical/r1i1p1f1/6hrPlevPt/psl/gn/v20191115/psl/psl_6hrPlevPt_ACCESS-ESM1-5_historical_r1i1p1f1_gn_195001010600-196001010000.nc"
;trackerFileName="/glade/collections/cmip/CMIP6/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r1i1p1f1/6hrPlevPt/psl/gn/v20190815/psl/psl_6hrPlevPt_MPI-ESM1-2-LR_historical_r1i1p1f1_gn_185001010600-187001010000.nc"
;trackerFileName="/glade/collections/cmip/CMIP6/CMIP/MIROC/MIROC6/historical/r1i1p1f1/6hrPlevPt/psl/gn/v20191204/psl/psl_6hrPlevPt_MIROC6_historical_r1i1p1f1_gn_201401010600-201501010000.nc"
;trackerFileName="/glade/u/home/zarzycki/scratch/ERA5/CAT.ERA5.U850.201701.nc"

topoFile=addfile(topoFileName,"r")
stride=64
topoLat=topoFile->lat(::stride)
topoLon=topoFile->lon(::stride)
topoTopo=topoFile->PHIS(::stride,::stride)
;topoTopo=topoTopo/9.81
print(topoLat)

if(doLonShift) then
  topoTopoShift=lonFlip(topoTopo)
  delete(topoTopo)
  topoTopo = topoTopoShift
  delete(topoTopoShift)
end if

trackerFile=addfile(trackerFileName,"r")

; first check if 2D
doCurvilinear=False
if (isfilevar(trackerFile, "XLAT")) then
  print("found XLAT, doing curvilinear")
  trackerLat=trackerFile->XLAT
  trackerLon=trackerFile->XLONG
  doCurvilinear=True
else
  print("doing RLL")
  trackerLat=trackerFile->lat
  if (isStoN) then
    trackerLat = trackerLat(::-1)
  end if
  trackerLon=trackerFile->lon
end if

if (doCurvilinear) then
  trackerTopo = rgrid2rcm_Wrap(topoTopo&lat,topoTopo&lon,topoTopo,trackerLat,trackerLon,1)
else
  trackerTopo = linint2_Wrap(topoTopo&lon,topoTopo&lat,topoTopo,isglobal,trackerLon,trackerLat,0)
  if (isStoN .and. keepStoN) then
    trackerTopo = trackerTopo(::-1,:)
  end if
end if

if (write_output) then
	print("Writing output file!")
	system("/bin/rm -f "+outFileName)   ; remove any pre-existing file
	ncdf     = addfile(outFileName ,"c")  ; open output netCDF file

	; make time and UNLIMITED dimension      ; recommended  for most applications
	filedimdef(ncdf,"time",-1,True) 

	ncdf->PHIS = trackerTopo         
end if

end
