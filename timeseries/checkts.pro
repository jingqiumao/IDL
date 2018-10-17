pro checkts
;MACR 14 MVK 13 HNO3 7 HCHO 20
input1='f0'
input2='f1'
input3='f1all'
tracenumber=20
tracername='CH2O'
date=20040707

dir='/as/scratch/2011-04/jmao/timeseries/'
infile1a=dir+'ts'+strtrim(string(date),2)+'.'+input1+'.bpch'
infile1b=dir+'ts'+strtrim(string(date+1),2)+'.'+input1+'.bpch'
infile1c=dir+'ts'+strtrim(string(date+2),2)+'.'+input1+'.bpch'
infile2a=dir+'ts'+strtrim(string(date),2)+'.'+input2+'.bpch'
infile2b=dir+'ts'+strtrim(string(date+1),2)+'.'+input2+'.bpch'
infile2c=dir+'ts'+strtrim(string(date+2),2)+'.'+input2+'.bpch'
infile3a=dir+'ts'+strtrim(string(date),2)+'.'+input3+'.bpch'
infile3b=dir+'ts'+strtrim(string(date+1),2)+'.'+input3+'.bpch'
infile3c=dir+'ts'+strtrim(string(date+2),2)+'.'+input3+'.bpch'


ModelInfo = ctm_type( 'GEOS5_47L', res = 2 )
GridInfo  = ctm_grid(ModelInfo)
ztop=12
; Find the index of the vertical gridbox nearest to the top of the
; domain, specified by ztop.
near_z = Min(abs(gridinfo.zmid-ztop),iztop)
zmid = GridInfo.zmid(0:iztop)

window,1
curtain1=get_curtain(infile1=infile1a,infile2=infile1b,tracername=tracername,tracenumber=tracenumber)
timearray=indgen(48)+1
tvcurtain_mao, curtain1, timearray, zmid, /ystyle, color = 1, $
           mindata=mmindata, maxdata=mmaxdata, ytitle=ytitle, title=title,$
           /CBar,/CBVertical,/NoAdvance,/NoZInterp
curtain=curtain1
save,curtain,filename=input1+'.'+tracername+'.'+strtrim(string(date),2)+'.sav'
ctm_cleanup

window,2
curtain2=get_curtain(infile1=infile2a,infile2=infile2b,tracername=tracername,tracenumber=tracenumber)
tvcurtain_mao, curtain2, timearray, zmid, /ystyle, color = 1, $
           mindata=mmindata, maxdata=mmaxdata, ytitle=ytitle, title=title,$
           /CBar,/CBVertical,/NoAdvance,/NoZInterp
curtain=curtain2
save,curtain,filename=input2+'.'+tracername+'.'+strtrim(string(date),2)+'.sav'
ctm_cleanup

window,3
curtain3=get_curtain(infile1=infile3a,infile2=infile3b,tracername=tracername,tracenumber=tracenumber)
tvcurtain_mao, curtain3, timearray, zmid, /ystyle, color = 1, $
           mindata=mmindata, maxdata=mmaxdata, ytitle=ytitle, title=title,$
           /CBar,/CBVertical,/NoAdvance,/NoZInterp
curtain=curtain3
save,curtain,filename=input3+'.'+tracername+'.'+strtrim(string(date),2)+'.sav'
ctm_cleanup
end

pro plotcurtain
;!X.OMARGIN = [0, 0]
;!Y.OMARGIN = [0, 0]
;!X.MARGIN = [0.5, 0.5]
;!Y.MARGIN = [1.3, 2.0]
;!P.CHARTHICK = 1.5
;!X.THICK = 2
;!Y.THICK = 2
;!P.font = 1
;!P.Charsize = 2.2
;
;nrow = 3 & ncol = 1
;!P.Multi = [0, ncol, nrow, 0, 1]
;multipanel,rows=nrow,cols=ncol

;1 NOx 3 PAN 7 HNO3 14 MACR 15 PMN 20 CH2O
input1='f0'
input2='f1'
input3='f1all'
tracenumber=3
tracername='PAN'
timearray=indgen(48)+1
date=20040707

ModelInfo = ctm_type( 'GEOS5_47L', res = 2 )
GridInfo  = ctm_grid(ModelInfo)
ztop=12
; Find the index of the vertical gridbox nearest to the top of the
; domain, specified by ztop.
near_z = Min(abs(gridinfo.zmid-ztop),iztop)
zmid = GridInfo.zmid(0:iztop)

window,1
file1=input1+'.'+tracername+'.'+strtrim(string(date),2)+'.sav'
restore,file1
print,'loading',file1
curtain1=curtain
timearray=indgen(48)+1
tvcurtain_mao, curtain1, timearray, zmid, /ystyle, color = 1, $
           mindata=mmindata, maxdata=mmaxdata, ytitle=ytitle, title='f0=0',$
           /CBar,/CBVertical,/NoAdvance,/NoZInterp
ctm_cleanup

;window,2
;file2=input2+'.'+tracername+'.'+strtrim(string(date),2)+'.sav'
;restore,file2
;print,'loading',file2
;curtain2=curtain
;tvcurtain_mao, curtain2, timearray, zmid, /ystyle, color = 1, $
;           mindata=mmindata, maxdata=mmaxdata, ytitle=ytitle, title=title,$
;           /CBar,/CBVertical,/NoAdvance,/NoZInterp
;ctm_cleanup

window,3
file3=input3+'.'+tracername+'.'+strtrim(string(date),2)+'.sav'
restore,file3
print,'loading',file3
curtain3=curtain
tvcurtain_mao, curtain3, timearray, zmid, /ystyle, color = 1, $
           mindata=mmindata, maxdata=mmaxdata, ytitle=ytitle, title='f0=1',$
           /CBar,/CBVertical,/NoAdvance,/NoZInterp
window,4
tvcurtain_mao, curtain3-curtain1, timearray, zmid, /ystyle, color = 1, $
           mindata=-0.05, maxdata=0.05, ytitle=ytitle, title='(f0=1)-(f0=0)',$
           /CBar,/CBVertical,/NoAdvance,/NoZInterp
stop
end

pro comparets
date=20040705
input1='f0'
input2='f1all'
tracenumber=1
tracername='NOx'
timearray=indgen(48)+1

ModelInfo = ctm_type( 'GEOS5_47L', res = 2 )
GridInfo  = ctm_grid(ModelInfo)
ztop=12
; Find the index of the vertical gridbox nearest to the top of the
; domain, specified by ztop.
near_z = Min(abs(gridinfo.zmid-ztop),iztop)
zmid = GridInfo.zmid(0:iztop)

window,4
restore,input1+'.'+tracername+'.'+strtrim(string(date),2)+'.sav'
curtainA=curtain
restore,input2+'.'+tracername+'.'+strtrim(string(date),2)+'.sav'
curtainB=curtain
tvcurtain_mao, curtainB-curtainA, timearray, zmid, /ystyle, color = 1, $
           mindata=mmindata, maxdata=mmaxdata, ytitle=ytitle, title=title,$
           /CBar,/CBVertical,/NoAdvance,/NoZInterp


ctm_cleanup
end
