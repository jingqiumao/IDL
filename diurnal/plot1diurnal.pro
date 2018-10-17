pro plot1diurnal
;open_device,  /ps, filename='diurnal.ps', $
;     Bits=8, /color
;open_device, winparam=[2,800,600] 
!X.OMARGIN = [2, 2]
!Y.OMARGIN = [2, 0]
!X.MARGIN = [2, 2]
!Y.MARGIN = [2, 2]
!P.CHARTHICK = 1.5
!P.CHARSIZE = 3
!X.THICK = 1
!Y.THICK = 1
!X.CHARSIZE = 1
!Y.CHARSIZE = 1
!P.font = 5

nrow=3
ncol=3
!P.Multi = [0, ncol, nrow, 0, 1]
multipanel,rows=nrow,cols=ncol

level=0

restore,'midlat.20040702.terp.fastPANdep.RO2NO3.sav'
ISOP2012=GC
restore,'midlat.20040702.dblterp.fastPANdep.RO2NO3.sav'
NO3sink=GC
plot,ISOP2012.Tau,ISOP2012.OH(*,level),title='OH',color=1
oplot,NO3sink.Tau,NO3sink.OH(*,level),color=2

plot,ISOP2012.Tau,ISOP2012.MACR(*,level),title='MACR',color=1,yrange=[0,1]
oplot,NO3sink.Tau,NO3sink.MACR(*,level),color=2

plot,ISOP2012.Tau,ISOP2012.MVK(*,level),title='MVK',color=1
oplot,NO3sink.Tau,NO3sink.MVK(*,level),color=2

plot,ISOP2012.Tau,ISOP2012.NO3(*,level),title='NO3',color=1
oplot,NO3sink.Tau,NO3sink.NO3(*,level),color=2

plot,ISOP2012.Tau,ISOP2012.HAC(*,level),title='HAC',color=1
oplot,NO3sink.Tau,NO3sink.HAC(*,level),color=2

plot,ISOP2012.Tau,ISOP2012.ALD2(*,level),title='ALD2',color=1,yrange=[0,5]
oplot,NO3sink.Tau,NO3sink.ALD2(*,level),color=2

plot,ISOP2012.Tau,ISOP2012.GLYC(*,level),title='GLYC',color=1
oplot,NO3sink.Tau,NO3sink.GLYC(*,level),color=2

plot,ISOP2012.Tau,ISOP2012.MOBA(*,level),title='MOBA',color=1,yrange=[0,2]
oplot,NO3sink.Tau,NO3sink.MOBA(*,level),color=2

plot,ISOP2012.Tau,ISOP2012.O3(*,level),title='O3',color=1
oplot,NO3sink.Tau,NO3sink.O3(*,level),color=2
print,'O3(*,0)',ISOP2012.O3(*,0)
print,'O3(*,0)',NO3sink.O3(*,0)
close,/all
;device, /close
end
pro plotnitrogen
;open_device,  /ps, filename='diurnal.ps', $
;     Bits=8, /color
;open_device, winparam=[2,800,600] 
!X.OMARGIN = [2, 2]
!Y.OMARGIN = [2, 0]
!X.MARGIN = [2, 2]
!Y.MARGIN = [2, 2]
!P.CHARTHICK = 1.5
!P.CHARSIZE = 3
!X.THICK = 1
!Y.THICK = 1
!X.CHARSIZE = 1
!Y.CHARSIZE = 1
!P.font = 5

nrow=3
ncol=3
!P.Multi = [0, ncol, nrow, 0, 1]
multipanel,rows=nrow,cols=ncolb

level=0

restore,'midlat.20040702.terp.fastPANdep.RO2NO3.sav'
ISOP2012=GC
restore,'midlat.20040702.dblterp.fastPANdep.RO2NO3.sav'
NO3sink=GC
plot,ISOP2012.Tau,ISOP2012.PAN(*,level),title='PAN',color=1,yrange=[0,2]
oplot,NO3sink.Tau,NO3sink.PAN(*,level),color=2

plot,ISOP2012.Tau,ISOP2012.PMN(*,level),title='PMN',color=1,yrange=[0,0.1]
oplot,NO3sink.Tau,NO3sink.PMN(*,level),color=2

plot,ISOP2012.Tau,ISOP2012.PPN(*,level),title='PPN',color=1,yrange=[0,0.3]
oplot,NO3sink.Tau,NO3sink.PPN(*,level),color=2

plot,ISOP2012.Tau,ISOP2012.R4N2(*,level),title='R4N2',color=1
oplot,NO3sink.Tau,NO3sink.R4N2(*,level),color=2

plot,ISOP2012.Tau,ISOP2012.ISOPN(*,level),title='ISOPN',color=1
oplot,NO3sink.Tau,NO3sink.ISOPN(*,level),color=2

plot,ISOP2012.Tau,ISOP2012.MMN(*,level),title='MMN',color=1
oplot,NO3sink.Tau,NO3sink.MMN(*,level),color=2

plot,ISOP2012.Tau,ISOP2012.HNO3(*,level),title='HNO3',yrange=[0,5],color=1
oplot,NO3sink.Tau,NO3sink.HNO3(*,level),color=2

plot,ISOP2012.Tau,ISOP2012.PROPNN(*,level),title='PROPNN',color=1
oplot,NO3sink.Tau,NO3sink.PROPNN(*,level),color=2

plot,ISOP2012.Tau,ISOP2012.NO2(*,level),title='NO2',color=1
oplot,NO3sink.Tau,NO3sink.NO2(*,level),color=2

close,/all
;device, /close
end
