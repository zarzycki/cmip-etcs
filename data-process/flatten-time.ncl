begin

;print("X"+filename+"X")
;filename="/glade/u/home/zarzycki/scratch/CMIPTMP/BCC-CSM2-MR/pr_FILT_3hr_BCC-CSM2-MR_historical_r1i1p1f1_gn_195001010130-195212312230.nc"
;print("X"+filename+"X")

f = addfile(filename,"w")

time = f->time
;time@calendar="standard"
print("BEFORE "+cd_calendar(time(0),3))

new_time_units = "days since 1900-01-01 00:00"
time = time_to_newtime_fix(time, new_time_units)

print("AFTER "+cd_calendar(time(0),3))

f->time = time

end