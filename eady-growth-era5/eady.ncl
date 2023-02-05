begin

;YYYYMMDD="20051231"
outbase="/glade/scratch/zarzycki/k1files/ERA5/"
dtime_map = (/4,2,2/)
tst_dtime = str_split_by_length(YYYYMMDD, dtime_map)
yyyy=tst_dtime(0)
mm=tst_dtime(1)
dd=tst_dtime(2)
print("Processing: "+yyyy+" "+mm+" "+dd)

times=(/0,6,12,18/)
plevels = (/ 100, 200, 300, 500, 700, 850, 925 /)

afile="/glade/u/home/zarzycki/rda/ds633.0/e5.oper.an.pl/"+yyyy+mm+"/e5.oper.an.pl.128_131_u.ll025uv."+yyyy+mm+dd+"00_"+yyyy+mm+dd+"23.nc"
bfile="/glade/u/home/zarzycki/rda/ds633.0/e5.oper.an.pl/"+yyyy+mm+"/e5.oper.an.pl.128_129_z.ll025sc."+yyyy+mm+dd+"00_"+yyyy+mm+dd+"23.nc"
cfile="/glade/u/home/zarzycki/rda/ds633.0/e5.oper.an.pl/"+yyyy+mm+"/e5.oper.an.pl.128_130_t.ll025sc."+yyyy+mm+dd+"00_"+yyyy+mm+dd+"23.nc"

a  = addfile(afile,"r")
b  = addfile(bfile,"r")
c  = addfile(cfile,"r")

; get variables
U = a->U(times,{plevels},:,:)
Z = b->Z(times,{plevels},:,:)
T = c->T(times,{plevels},:,:)

lev = T&level
; convert lev from hPa to Pa
lev = lev * 100.

; calculate theta
th = pot_temp(lev, T, 1, False)

; conform latitude
XLAT = conform(T, T&latitude, 2)

; egr function
egr = eady_growth_rate(th, U, Z, XLAT, 0,  1)
printMinMax(egr, 0)

; write file

outdir = outbase+"/"+yyyy
outfile = "ERA5.k1."+yyyy+mm+dd+".nc"

system("mkdir -p -v "+outdir)   ; remove any pre-existing file

system("/bin/rm -f -v "+outdir+"/"+outfile)   ; remove any pre-existing file

;setfileoption ("nc","Format","NetCDF4")
ncdf = addfile(outdir+"/"+outfile ,"c")  ; open output netCDF file
;setfileoption (ncdf,"CompressionLevel",1)

fAtt               = True            ; assign file attributes
fAtt@source_file   = afile
fAtt@creation_date = systemfunc ("date")
fileattdef( ncdf, fAtt )            ; copy file attributes

filedimdef(ncdf,"time",-1,True) 

ncdf->egr = egr                ; name on file different from name in code

print("Wrote NetCDF at: "+outdir+"/"+outfile)

end