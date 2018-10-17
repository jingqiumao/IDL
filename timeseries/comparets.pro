pro comparets
input1='f0'
input3='f1Spin'

ModelInfo = ctm_type( 'GEOS5_47L', res = 2 )
GridInfo  = ctm_grid(ModelInfo)
ztop=12
; Find the index of the vertical gridbox nearest to the top of the
; domain, specified by ztop.
near_z = Min(abs(gridinfo.zmid-ztop),iztop)
zmid = GridInfo.zmid(0:iztop)

nlon=144
nlat=47
nlev=47

;dates=['20040701','20040702','20040703','20040704','20040705','20040706','20040707','20040708','20040709']
;tracers=[ 1, 2, 3, 4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25]

dates=['20040710']
tracers=1
tracername=['NOx']
;tracername=['NOx','Ox','PAN','CO','ALK4','ISOP','HNO3','H2O2','ACET','MEK','ALD2','RCHO','MVK','MACR',$
;            'PMN','PPN','R4N2','PRPE','C3H8','CH2O','C2H6','N2O5','HNO4','MP','DMS']
;DiagN=

dir='/as/scratch/2011-04/jmao/ts_compare/'
for idate=0,n_elements(dates)-1  do begin
infile1=dir+'ts'+strtrim(string(dates(idate)),2)+'.'+input3+'.bpch'

tmp=fltarr(nlon,nlat,nlev,24)&tmp(*,*,*)=!values.f_NaN
O3a=fltarr(nlon,nlat,nlev,24)&O3a(*,*,*)=!values.f_NaN

ctm_diaginfo,/all,/force_reading,filename=dir+'diaginfo.dat'
ctm_tracerinfo,/all,/force_reading, filename=dir+'tracerinfo.dat'


 CTM_Get_Data, DataInfo1,File=InFile1,'IJ-AVG-$',tracer=tracernumber
;CTM_Get_Data, DataInfo2,File=InFile1,'TIME-SER',tracer=2

  Tau0a = DataInfo1[*].Tau0
  Tau0a = Tau0a[Uniq( Tau0a, Sort( Tau0a ) ) ]
;  undefine,datainfo 
; get time in the format of YYYYMMDDHH
time1=lonarr(24)
for n_tau=0,n_elements(Tau0a)-1 do begin
  t=tau2yymmdd(Tau0a[n_tau],/nformat)&time1[n_tau]=t[1]/10000;time[n_tau]=t[0]*100+t[1]/10000
  ;end of the loop for Tau0
endfor


for j=0,n_elements(tracername)-1 do begin
;  CTM_Get_Data, DataInfo1,File=InFile1,'IJ-AVG-$'
  ; THISDATAINFO is an array of data blocks for time TAU0
  ;O3 is ppbv
  Ind = Where( DataInfo1[*].Tracername eq tracername(j) )
  for n=0,n_elements(Ind)-1 do begin
    i=Ind[n] 
    if Datainfo1[i].tau0 eq Tau0a[n] then begin
    tmp[*,*,*,n]=(*(DataInfo1[i].data))*1e-9&O3a[*,*,*,n]=tmp[*,*,*,n]
    endif else begin
    print,'problems in tau'&stop
    endelse
  endfor
tmp=mean(o3a,4)
;endfor for the tracername
endfor
;endfor for the dates
endfor
stop

end

