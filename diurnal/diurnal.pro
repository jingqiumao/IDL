pro diurnal

;restore,'~/INTEXA/run.std2011.cents.fast/midlatfast.sav'
restore,'/mnt/as/home2/j/jmao/INTEXA/IDL/barkley/OVOCoutput/fast_NO3uptk_snowscavHNO3/midlatfast.sav'
JD         = tau2yymmdd( newGC.tau )
maxPBL=fltarr(n_elements(JD.day))&maxPBL(*)=!values.f_NaN
UTC  = JULDAY( JD.month , JD.day , JD.year , JD.hour , JD.minute )
Longitude = -86
LocalTime=newGC.tau+(Longitude/15)
JD2         = tau2yymmdd( LocalTime ) 

index=uniq(JD2.day)

 for i=0,n_elements(index)-1 do begin
    index1=where(JD2.day eq JD2.day(index(i)))
    maxPBL(index1)=max(newGC.PBL(index1))
 endfor

highpbl=where (maxPBL ge 100)
;highpbl=where (maxPBL ge 1500)

level=0

NOx_median=tapply(newGC.NOx(highpbl,level),JD2.hour(highpbl),'median')
O3_median=tapply(newGC.O3(highpbl,level),JD2.hour(highpbl),'median')
ISOPemis_median=tapply(newGC.ISOPemis(highpbl,0),JD2.hour(highpbl),'median')
ISOP_median=tapply(newGC.ISOP(highpbl,level),JD2.hour(highpbl),'median')
MACR_median=tapply(newGC.MACR(highpbl,level),JD2.hour(highpbl),'median')
MVK_median=tapply(newGC.MVK(highpbl,level),JD2.hour(highpbl),'median')
PAN_median=tapply(newGC.PAN(highpbl,level),JD2.hour(highpbl),'mean')
NOy_median=tapply(newGC.NOy(highpbl,level),JD2.hour(highpbl),'median')
HNO3_median=tapply(newGC.HNO3(highpbl,level),JD2.hour(highpbl),'mean')
NO3_median=tapply(newGC.NO3(highpbl,level),JD2.hour(highpbl),'mean')
OH_median=tapply(newGC.OH(highpbl,level),JD2.hour(highpbl),'median')
CH2O_median=tapply(newGC.CH2O(highpbl,level),JD2.hour(highpbl),'median')
ALD2_median=tapply(newGC.ALD2(highpbl,level),JD2.hour(highpbl),'median')
N2O5_median=tapply(newGC.N2O5(highpbl,level),JD2.hour(highpbl),'median')
MPAN_median=tapply(newGC.PMN(highpbl,level),JD2.hour(highpbl),'median')
pbl_median=tapply(newGC.pbl(highpbl),JD2.hour(highpbl),'median')
ratio_median=tapply(newGC.MVK(highpbl,level)/newGC.MACR(highpbl,level),JD2.hour(highpbl),'median')
ratio2_median=tapply(newGC.MVK(highpbl,level)/newGC.ISOP(highpbl,level),JD2.hour(highpbl),'median')
ratio3_median=tapply(newGC.MACR(highpbl,level)/newGC.ISOP(highpbl,level),JD2.hour(highpbl),'median')
HNO3dvel_median=tapply(newGC.HNO3dvel(highpbl),JD2.hour(highpbl),'mean')

hour_median=findgen(24)+0.5
stop
window,0
plot,hour_median,NOx_median,title='NOx',xtitle='hour',ytitle='ppb'
window,1
plot,hour_median,O3_median,color=1,title='O3',xtitle='hour',ytitle='ppb'
window,2
plot,hour_median,ISOPemis_median,color=1,title='ISOPemis'
window,3
plot,hour_median,ISOP_median/5,color=1,title='ISOP',xtitle='hour',ytitle='ppb'
window,4
plot,hour_median,MACR_median,color=1,title='MACR',xtitle='hour',ytitle='ppb'
window,5
plot,hour_median,ratio_median,color=1,title='MVK/MACR',xtitle='hour'
window,6
plot,hour_median,ratio2_median,color=1,title='MVK/ISOP'
window,7
plot,hour_median,ratio3_median,color=1,title='MACR/ISOP'
window,8
plot,hour_median,HNO3_median,color=1,title='HNO3&PAN',xtitle='hour',ytitle='ppb'
oplot,hour_median,PAN_median,color=2
legend, lcolor=[1,2],line=[0,0], lthick = [1, 1], $
        label=['HNO3', 'PAN'],$
        halign=0.9, valign=0.9, charsize=1.2, /colo
;oplot,hour_median,HNO3_median,color=3
window,9
plot,hour_median, N2O5_median,color=1,title='N2O5',xtitle='hour',ytitle='ppb'
window,10
plot,hour_median, CH2O_median,color=1,title='CH2O'
window,11
plot,hour_median, MPAN_median,color=1,title='MPAN',xtitle='hour',ytitle='ppb'
window,12
plot,hour_median, pbl_median,color=1,title='PBLheight',xtitle='Hour',ytitle='m'
window,13
plot,hour_median, NO3_median,color=1,title='NO3',xtitle='Hour',ytitle='ppb'
window,14
plot,hour_median, HNO3dvel_median,color=1,title='HNO3dvel',xtitle='Hour',ytitle='cm/s'
window,15
plot,hour_median, Oxdvel_median,color=1,title='O3dvel',xtitle='Hour',ytitle='cm/s'

end

pro comp2diurnal
restore,'/mnt/as/home2/j/jmao/INTEXA/IDL/barkley/OVOCoutput/fast/midlatfast.sav'
;restore,'/mnt/as/home2/j/jmao/INTEXA/IDL/barkley/OVOCoutput/fast_N2O5uptake_NO3chem/midlatfast.sav'
;restore,'/mnt/as/home2/j/jmao/INTEXA/IDL/barkley/OVOCoutput/fast_NO3uptk/midlatfast.sav'
;restore,'/mnt/as/home2/j/jmao/INTEXA/IDL/barkley/OVOCoutput/slow/midlatslow.sav'
fast1=newGC

restore,'/mnt/as/home2/j/jmao/INTEXA/IDL/barkley/OVOCoutput/fast_NO3uptk_snowscavHNO3/midlatfast.sav'
;restore,'/mnt/as/home2/j/jmao/INTEXA/IDL/barkley/OVOCoutput/fast_N2O5uptake_NO3chem/midlatfast.sav'
;restore,'/mnt/as/home2/j/jmao/INTEXA/IDL/barkley/OVOCoutput/slow/midlatslow.sav'
fast2=newGC

JD         = tau2yymmdd( newGC.tau )
maxPBL=fltarr(n_elements(JD.day))&maxPBL(*)=!values.f_NaN
UTC  = JULDAY( JD.month , JD.day , JD.year , JD.hour , JD.minute )
Longitude = -86
LocalTime=newGC.tau+(Longitude/15)
JD2         = tau2yymmdd( LocalTime ) 

index=uniq(JD2.day)

 for i=0,n_elements(index)-1 do begin
    index1=where(JD2.day eq JD2.day(index(i)))
    maxPBL(index1)=max(newGC.PBL(index1))
 endfor

highpbl=where (maxPBL ge 100)
;highpbl=where (maxPBL ge 1500)

level=0

NOx_median1=tapply(fast1.NOx(highpbl,level),JD2.hour(highpbl),'mean')
NOx_median2=tapply(fast2.NOx(highpbl,level),JD2.hour(highpbl),'mean')
O3_median1=tapply(fast1.O3(highpbl,level),JD2.hour(highpbl),'mean')
O3_median2=tapply(fast2.O3(highpbl,level),JD2.hour(highpbl),'mean')
OH_median1=tapply(fast1.OH(highpbl,level),JD2.hour(highpbl),'mean')
OH_median2=tapply(fast2.OH(highpbl,level),JD2.hour(highpbl),'mean')
HNO3_median1=tapply(fast1.HNO3(highpbl,level),JD2.hour(highpbl),'mean')
HNO3_median2=tapply(fast2.HNO3(highpbl,level),JD2.hour(highpbl),'mean')
MVK_median1=tapply(fast1.MVK(highpbl,level),JD2.hour(highpbl),'mean')
MVK_median2=tapply(fast2.MVK(highpbl,level),JD2.hour(highpbl),'mean')
MACR_median1=tapply(fast1.MACR(highpbl,level),JD2.hour(highpbl),'mean')
MACR_median2=tapply(fast2.MACR(highpbl,level),JD2.hour(highpbl),'mean')
N2O5_median1=tapply(fast1.N2O5(highpbl,level),JD2.hour(highpbl),'mean')
N2O5_median2=tapply(fast2.N2O5(highpbl,level),JD2.hour(highpbl),'mean')
PAN_median1=tapply(fast1.PAN(highpbl,level),JD2.hour(highpbl),'mean')
PAN_median2=tapply(fast2.PAN(highpbl,level),JD2.hour(highpbl),'mean')
Ox_median1=tapply(fast1.Ox(highpbl,level),JD2.hour(highpbl),'mean')
Ox_median2=tapply(fast2.Ox(highpbl,level),JD2.hour(highpbl),'mean')
NOy_median1=tapply(fast1.NOy(highpbl,level),JD2.hour(highpbl),'mean')
NOy_median2=tapply(fast2.NOy(highpbl,level),JD2.hour(highpbl),'mean')
ISOP_median1=tapply(fast1.ISOP(highpbl,level),JD2.hour(highpbl),'mean')
ISOP_median2=tapply(fast2.ISOP(highpbl,level),JD2.hour(highpbl),'mean')
CH2O_median1=tapply(fast1.CH2O(highpbl,level),JD2.hour(highpbl),'mean')
CH2O_median2=tapply(fast2.CH2O(highpbl,level),JD2.hour(highpbl),'mean')
pbl_median1=tapply(fast1.pbl(highpbl),JD2.hour(highpbl),'mean')
pbl_median2=tapply(fast2.pbl(highpbl),JD2.hour(highpbl),'mean')
HNO3dvel_median1=tapply(fast1.HNO3dvel(highpbl),JD2.hour(highpbl),'mean')
HNO3dvel_median2=tapply(fast2.HNO3dvel(highpbl),JD2.hour(highpbl),'mean')
ISOPemis_median1=tapply(fast1.ISOPemis(highpbl,0),JD2.hour(highpbl),'mean')
ISOPemis_median2=tapply(fast2.ISOPemis(highpbl,0),JD2.hour(highpbl),'mean')

hour_median=findgen(24)+0.5
print,'O3_1:',O3_median1
print,'O3_2:',O3_median2
stop
window,0
plot,hour_median,O3_median1,title='O3',xtitle='hour',ytitle='ppb'
oplot,hour_median,O3_median2,color=2
legend, lcolor=[1,2],line=[0,0], lthick = [1, 1], $
        label=['GEOS-4', 'GEOS-5'],$
        halign=0.9, valign=0.9, charsize=1.2, /colo
window,1
plot,hour_median,NOx_median1,color=1,title='NOx',xtitle='hour',ytitle='ppb',yrange=[0,7]
oplot,hour_median,NOx_median2,color=2
window,2
plot,hour_median,HNO3_median1,color=1,title='HNO3',xtitle='hour',ytitle='ppb',yrange=[0,2]
oplot,hour_median,HNO3_median2,color=2
window,3
plot,hour_median,OH_median1,color=1,title='OH',xtitle='hour',ytitle='molec/cm3/s',yrange=[0,1e7]
oplot,hour_median,OH_median2,color=2
window,4
plot,hour_median,PAN_median1,color=1,title='PAN',xtitle='hour',ytitle='ppb'
oplot,hour_median,PAN_median2,color=2
window,5
plot,hour_median,MVK_median1,color=1,title='MVK',xtitle='hour',ytitle='ppb'
oplot,hour_median,MVK_median2,color=2
window,6
plot,hour_median,MACR_median1,color=1,title='MACR',xtitle='hour',ytitle='ppb'
oplot,hour_median,MACR_median2,color=2
window,7
plot,hour_median,N2O5_median1,color=1,title='N2O5',xtitle='hour',ytitle='ppb'
oplot,hour_median,N2O5_median2,color=2
window,8
plot,hour_median,Ox_median1,color=1,title='Ox',xtitle='hour',ytitle='ppb'
oplot,hour_median,Ox_median2,color=2
window,9
plot,hour_median,NOy_median1,color=1,title='NOy',xtitle='hour',ytitle='ppb'
oplot,hour_median,NOy_median2,color=2
window,10
plot,hour_median,ISOP_median1/5,color=1,title='ISOP',xtitle='hour',ytitle='ppb',yrange=[0,3]
oplot,hour_median,ISOP_median2/5,color=2
window,11
plot,hour_median,pbl_median1,color=1,title='pbl',xtitle='hour',ytitle='m',yrange=[0,2000]
oplot,hour_median,pbl_median2,color=2
window,12
plot,hour_median,HNO3dvel_median1,color=1,title='HNO3dvel',xtitle='hour',ytitle='m/s',yrange=[0,3.5]
oplot,hour_median,HNO3dvel_median2,color=2
window,13
plot,hour_median,ISOPemis_median1,color=1,title='ISOPemis',xtitle='hour',ytitle='m/s'
oplot,hour_median,ISOPemis_median2,color=2
window,14
plot,hour_median,CH2O_median1,color=1,title='CH2O',xtitle='hour',ytitle='ppb'
oplot,hour_median,CH2O_median2,color=2
end

pro diurnal4proposal
!X.OMARGIN = [3, 3]
!Y.OMARGIN = [5, 5]
!X.MARGIN = [3, 3]
!Y.MARGIN = [2.0, 2.0]
!P.CHARTHICK = 1.5
!X.THICK = 2
!Y.THICK = 2
!P.font = 1
!P.Charsize = 1.8


nrow = 1 & ncol = 2
!P.Multi = [0, ncol, nrow, 0, 0]
multipanel,rows=nrow,cols=ncol

restore,'/mnt/as/home2/j/jmao/INTEXA/IDL/barkley/OVOCoutput/fast/midlatfast.sav'
fast1=newGC

restore,'/mnt/as/home2/j/jmao/INTEXA/IDL/barkley/OVOCoutput/fast_NO3uptk_snowscavHNO3/midlatfast.sav'
fast2=newGC

JD         = tau2yymmdd( newGC.tau )
maxPBL=fltarr(n_elements(JD.day))&maxPBL(*)=!values.f_NaN
UTC  = JULDAY( JD.month , JD.day , JD.year , JD.hour , JD.minute )
Longitude = -86
LocalTime=newGC.tau+(Longitude/15)
JD2         = tau2yymmdd( LocalTime ) 

index=uniq(JD2.day)

 for i=0,n_elements(index)-1 do begin
    index1=where(JD2.day eq JD2.day(index(i)))
    maxPBL(index1)=max(newGC.PBL(index1))
 endfor

highpbl=where (maxPBL ge 100)
;highpbl=where (maxPBL ge 1500)

level=0

O3_median1=tapply(fast1.O3(highpbl,level),JD2.hour(highpbl),'mean')
O3_median2=tapply(fast2.O3(highpbl,level),JD2.hour(highpbl),'mean')
HNO3_median1=tapply(fast1.HNO3(highpbl,level),JD2.hour(highpbl),'mean')
HNO3_median2=tapply(fast2.HNO3(highpbl,level),JD2.hour(highpbl),'mean')
hour_median=findgen(24)+0.5
print,'O3_1:',O3_median1
print,'O3_2:',O3_median2
print,'diff',O3_median1-O3_median2
;window,0
plot,hour_median,O3_median1,title='O3',xtitle='hour',ytitle='ppb',yrange=[30,70]
oplot,hour_median,O3_median2,color=2
;window,2
plot,hour_median,HNO3_median1,color=1,title='HNO3',xtitle='hour',ytitle='ppb',yrange=[0,3]
oplot,hour_median,HNO3_median2,color=2

legend, lcolor=[1,2],line=[0,0],lthick=[5,5],$
        label=['std chem','improved'], $
        position=[0.70,0.6,0.75,0.7], charsize=1.8, /color
end
