function get_curtain,infile1=infile1,infile2=infile2,tracenumber=tracenumber,tracername=tracername
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

tmp=fltarr(nlon,nlat,nlev,24)&tmp(*,*,*)=!values.f_NaN
O3a=fltarr(nlon,nlat,nlev,24)&O3a(*,*,*)=!values.f_NaN
O3b=fltarr(nlon,nlat,nlev,24)&O3b(*,*,*)=!values.f_NaN


 CTM_Get_Data, DataInfo1,File=InFile1,'IJ-AVG-$',tracer=tracernumber
 CTM_Get_Data, DataInfo2,File=InFile2,'IJ-AVG-$',tracer=tracernumber

  Tau0a = DataInfo1[*].Tau0
  Tau0a = Tau0a[Uniq( Tau0a, Sort( Tau0a ) ) ]
;  undefine,datainfo 
; get time in the format of YYYYMMDDHH
time1=lonarr(24)
for n_tau=0,n_elements(Tau0a)-1 do begin
  t=tau2yymmdd(Tau0a[n_tau],/nformat)&time1[n_tau]=t[1]/10000;time[n_tau]=t[0]*100+t[1]/10000
  ;end of the loop for Tau0
endfor

  Tau0b = DataInfo2[*].Tau0
  Tau0b = Tau0b[ Uniq( Tau0b, Sort( Tau0b ) ) ]
time2=lonarr(24)
for n_tau=0,n_elements(Tau0b)-1 do begin
  t=tau2yymmdd(Tau0b[n_tau],/nformat)&time2[n_tau]=t[1]/10000;time[n_tau]=t[0]*100+t[1]/10000
  ;end of the loop for Tau0
endfor


;  CTM_Get_Data, DataInfo1,File=InFile1,'IJ-AVG-$'
  ; THISDATAINFO is an array of data blocks for time TAU0
  ;O3 is ppbv
  Ind = Where( DataInfo1[*].Tracername eq tracername )
  for n=0,n_elements(Ind)-1 do begin
    i=Ind[n] 
    if Datainfo1[i].tau0 eq Tau0a[n] then begin
    tmp[*,*,*,n]=(*(DataInfo1[i].data))*1e-9&O3a[*,*,*,n]=tmp[*,*,*,n]
    endif else begin
    print,'problems in tau'&stop
    endelse
  endfor

;  CTM_Get_Data, DataInfo2,File=InFile2,'IJ-AVG-$'
  ; THISDATAINFO is an array of data blocks for time TAU0
  ;O3 is ppbv
  Ind = Where( DataInfo2[*].Tracername eq tracername )
  for n=0,n_elements(Ind)-1 do begin
    i=Ind[n] 
    if Datainfo2[i].tau0 eq Tau0b[n] then begin
    tmp[*,*,*,n]=(*(DataInfo2[i].data))*1e-9&O3b[*,*,*,n]=tmp[*,*,*,n]
    endif else begin
    print,'problems in tau'&stop
    endelse
  endfor

;90w is -85, and index is 39? 35N is 16 (since I only output 45 to 91) 
curtain1=O3a[39,16,0:iztop,*]*1e9
curtain2=O3b[39,16,0:iztop,*]*1e9
timearray=[time1,time2]
curtain=[transpose(reform(curtain1)),transpose(reform(curtain2))]

undefine, DataInfo1
undefine, DataInfo2
return, curtain
end
