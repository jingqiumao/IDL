;my idea here is to make curtain for each day in makecurtain by saving
;to GC structure
;then merge them in mergecurtain
;finally plot them in mycurtain


pro makecurtain
;input3='slow'
  ;90w is -85, and index is 39? 35N is 16 (since I only output 45 to 91)
  ;curtain1=O3a[39,16,0:iztop,*]*1e9
  ;-60 index is 49, -2s index is 0.
  ;input='Amazon'
  ;ilon=49&ilat=0
  input='midlat'
  ilon=39&ilat=16

ModelInfo = ctm_type( 'GEOS5_47L', res = 2 )
GridInfo  = ctm_grid(ModelInfo)
ztop=10
; Find the index of the vertical gridbox nearest to the top of the
; domain, specified by ztop.
near_z = Min(abs(gridinfo.zmid-ztop),iztop)
zmid = GridInfo.zmid(0:iztop)

nlon=144
nlat=47
nlev=30
;
;dates=['20040701','20040702','20040703','20040704','20040705','20040706','20040707',$
;      '20040708','20040709','20040710','20040711','20040712','20040713','20040714',$
;      '20040715','20040716','20040717','20040718','20040719','20040720','20040721',$
;      '20040722','20040723','20040724','20040725','20040726','20040727','20040728',$
;      '20040728','20040729','20040730']
dates=['20040701','20040702']
tracers=[ 1, 2, 3, 4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,44,45,46,47,48,49,50,51,52,53,54,55,56]
;for 52-tracer isoprene
;tracername=['NOx','Ox','PAN','CO','ALK4','ISOP','HNO3','H2O2','ACET','MEK','ALD2','RCHO','MVK','MACR',$
;            'PMN','PPN','R4N2','PRPE','C3H8','CH2O','C2H6','N2O5','HNO4','MP','DMS','ISOPN','MOBA','PROPNN',$
;            'HAC','GLYC','MMN','RIP','IEPOX','MAP','API','LIM','ONIT','O3']
;for 55-tracer terpene
tracername=['NOx','Ox','PAN','CO','ALK4','ISOP','HNO3','H2O2','ACET','MEK','ALD2','RCHO','MVK','MACR',$
            'PMN','PPN','R4N2','PRPE','C3H8','CH2O','C2H6','N2O5','HNO4','MP','DMS','ISOPN','MOBA','PROPNN',$
            'HAC','GLYC','MMN','RIP','IEPOX','MAP','API','LIM','ONIT','O3']
;DiagN=

;dir='/as/scratch/2011-04/jmao/ts_compare/'
dir='~/INTEXA/test.terp/'
ctm_diaginfo,/all,/force_reading,filename=dir+'diaginfo.dat'
ctm_tracerinfo,/all,/force_reading, filename=dir+'tracerinfo.dat'

for idate=0,n_elements(dates)-1  do begin
infile1=dir+'ts'+strtrim(string(dates(idate)),2)+'.bpch'

tmp=fltarr(nlon,nlat,nlev,24)&tmp(*,*,*)=!values.f_NaN
O3a=fltarr(nlon,nlat,nlev,24)&O3a(*,*,*)=!values.f_NaN
O3b=fltarr(nlon,nlat,24)&O3b(*,*,*)=!values.f_NaN
tmpb=fltarr(nlon,nlat,24)&tmpb(*,*,*)=!values.f_NaN


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


  curtain1=O3a[ilon,ilat,0:iztop,*]*1e9

  timearray=[time1]
  curtain=[transpose(reform(curtain1))]
  s=tracername(j)+'=curtain'
  result=execute(s)

;for tracers
endfor

;OH(molec/cm3) is 2, NOy (ppbv) is 3, NO2(ppbv) is 25,NO(ppbv) is 9
CTM_Get_Data, DataInfo2,File=InFile1,'TIME-SER',tracer=2;unit is m

;here I use a differen strategy, which seesm simpler.
  for n=0,n_elements(DataInfo2)-1 do begin
    if Datainfo2[n].tau0 eq Tau0a[n] then begin
    tmp[*,*,*,n]=(*(DataInfo2[n].data))&O3a[*,*,*,n]=tmp[*,*,*,n]
    endif else begin
    print,'problems in tau'&stop
    endelse
  endfor

  curtain1=O3a[ilon,ilat,0:iztop,*]

  OH=[transpose(reform(curtain1))]

CTM_Get_Data, DataInfo2,File=InFile1,'TIME-SER',tracer=3;unit is m

;here I use a differen strategy, which seesm simpler.
  for n=0,n_elements(DataInfo2)-1 do begin
    if Datainfo2[n].tau0 eq Tau0a[n] then begin
    tmp[*,*,*,n]=(*(DataInfo2[n].data))&O3a[*,*,*,n]=tmp[*,*,*,n]
    endif else begin
    print,'problems in tau'&stop
    endelse
  endfor

  curtain1=O3a[ilon,ilat,0:iztop,*]

  NOy=[transpose(reform(curtain1))]

CTM_Get_Data, DataInfo2,File=InFile1,'TIME-SER',tracer=9;unit is m

;here I use a differen strategy, which seesm simpler.
  for n=0,n_elements(DataInfo2)-1 do begin
    if Datainfo2[n].tau0 eq Tau0a[n] then begin
    tmp[*,*,*,n]=(*(DataInfo2[n].data))&O3a[*,*,*,n]=tmp[*,*,*,n]
    endif else begin
    print,'problems in tau'&stop
    endelse
  endfor

  curtain1=O3a[ilon,ilat,0:iztop,*]

  NO=[transpose(reform(curtain1))]

CTM_Get_Data, DataInfo2,File=InFile1,'TIME-SER',tracer=25;unit is m

;here I use a differen strategy, which seesm simpler.
  for n=0,n_elements(DataInfo2)-1 do begin
    if Datainfo2[n].tau0 eq Tau0a[n] then begin
    tmp[*,*,*,n]=(*(DataInfo2[n].data))&O3a[*,*,*,n]=tmp[*,*,*,n]
    endif else begin
    print,'problems in tau'&stop
    endelse
  endfor

  curtain1=O3a[ilon,ilat,0:iztop,*]

  NO2=[transpose(reform(curtain1))]

CTM_Get_Data, DataInfo2,File=InFile1,'TIME-SER',tracer=14;unit is m

;here I use a differen strategy, which seesm simpler.
  for n=0,n_elements(DataInfo2)-1 do begin
    if Datainfo2[n].tau0 eq Tau0a[n] then begin
    tmp[*,*,*,n]=(*(DataInfo2[n].data))&O3a[*,*,*,n]=tmp[*,*,*,n]
    endif else begin
    print,'problems in tau'&stop
    endelse
  endfor

  curtain1=O3a[ilon,ilat,0:iztop,*]

  NO3=[transpose(reform(curtain1))]

CTM_Get_Data, DataInfo2,File=InFile1,'TIME-SER',tracer=34;unit is m

;here I use a differen strategy, which seesm simpler.
  for n=0,n_elements(DataInfo2)-1 do begin
    if Datainfo2[n].tau0 eq Tau0a[n] then begin
    tmp[*,*,*,n]=(*(DataInfo2[n].data))&O3a[*,*,*,n]=tmp[*,*,*,n]
    endif else begin
    print,'problems in tau'&stop
    endelse
  endfor

  curtain1=O3a[ilon,ilat,0:iztop,*]

  HO2=[transpose(reform(curtain1))]
CTM_Get_Data, DataInfo2,File=InFile1,'BIOGSRCE',tracer=7;unit is m

;here I use a differen strategy, which seesm simpler.
  for n=0,n_elements(DataInfo2)-1 do begin
    if Datainfo2[n].tau0 eq Tau0a[n] then begin
    tmp[*,*,*,n]=(*(DataInfo2[n].data))&O3a[*,*,*,n]=tmp[*,*,*,n]
    endif else begin
    print,'problems in tau'&stop
    endelse
  endfor

  curtain1=O3a[ilon,ilat,0:iztop,*]

  APIE=[transpose(reform(curtain1))]

CTM_Get_Data, DataInfo2,File=InFile1,'BIOGSRCE',tracer=8;unit is m

;here I use a differen strategy, which seesm simpler.
  for n=0,n_elements(DataInfo2)-1 do begin
    if Datainfo2[n].tau0 eq Tau0a[n] then begin
    tmp[*,*,*,n]=(*(DataInfo2[n].data))&O3a[*,*,*,n]=tmp[*,*,*,n]
    endif else begin
    print,'problems in tau'&stop
    endelse
  endfor

  curtain1=O3a[ilon,ilat,0:iztop,*]

  BPIE=[transpose(reform(curtain1))]
CTM_Get_Data, DataInfo2,File=InFile1,'BIOGSRCE',tracer=9;unit is m

;here I use a differen strategy, which seesm simpler.
  for n=0,n_elements(DataInfo2)-1 do begin
    if Datainfo2[n].tau0 eq Tau0a[n] then begin
    tmp[*,*,*,n]=(*(DataInfo2[n].data))&O3a[*,*,*,n]=tmp[*,*,*,n]
    endif else begin
    print,'problems in tau'&stop
    endelse
  endfor

  curtain1=O3a[ilon,ilat,0:iztop,*]

  LIME=[transpose(reform(curtain1))]

CTM_Get_Data, DataInfo2,File=InFile1,'BXHGHT-$',tracer=1;unit is m

;here I use a differen strategy, which seesm simpler.
  for n=0,n_elements(DataInfo2)-1 do begin
    if Datainfo2[n].tau0 eq Tau0a[n] then begin
    tmp[*,*,*,n]=(*(DataInfo2[n].data))&O3a[*,*,*,n]=tmp[*,*,*,n]
    endif else begin
    print,'problems in tau'&stop
    endelse
  endfor

  curtain1=O3a[ilon,ilat,0:iztop,*]

  boxheight=[transpose(reform(curtain1))]

CTM_Get_Data, DataInfo3,File=InFile1,'PBLDEPTH',tracer=1;unit is m

;here I use a differen strategy, which seesm simpler.
  for n=0,n_elements(DataInfo3)-1 do begin
    if Datainfo3[n].tau0 eq Tau0a[n] then begin
    tmpb[*,*,n]=(*(DataInfo3[n].data))&O3b[*,*,n]=tmpb[*,*,n]
    endif else begin
    print,'problems in tau'&stop
    endelse
  endfor

  ;90w is -85, and index is 39? 35N is 16 (since I only output 45 to 91) 
  curtain2=O3b[ilon,ilat,*]

  PBLheight=[reform(curtain2)]

CTM_Get_Data, DataInfo3,File=InFile1,'PBLDEPTH',tracer=2;unit is m

;here I use a differen strategy, which seesm simpler.
  for n=0,n_elements(DataInfo3)-1 do begin
    if Datainfo3[n].tau0 eq Tau0a[n] then begin
    tmpb[*,*,n]=(*(DataInfo3[n].data))&O3b[*,*,n]=tmpb[*,*,n]
    endif else begin
    print,'problems in tau'&stop
    endelse
  endfor

  ;90w is -85, and index is 39? 35N is 16 (since I only output 45 to 91) 
  curtain2=O3b[ilon,ilat,*]

  PBLlevel=[reform(curtain2)]



      ; Create structure with fields common to all simulations
      GC = { Tau:Tau0a, OH:OH, NO:NO, NO2:NO2, HO2:HO2, NO3:NO3,NOy:NOy,NOx:NOx,Ox:Ox,O3:O3,PAN:PAN,CO:CO,ALK4:ALK4,$
            ISOP:ISOP,HNO3:HNO3,H2O2:H2O2,ACET:ACET,MEK:MEK,ALD2:ALD2,RCHO:RCHO,MVK:MVK,MACR:MACR,$
            PMN:PMN,PPN:PPN,R4N2:R4N2,PRPE:PRPE,C3H8:C3H8,CH2O:CH2O,C2H6:C2H6,N2O5:N2O5,HNO4:HNO4,MP:MP,DMS:DMS,ISOPN:ISOPN,$
            MOBA:MOBA,PROPNN:PROPNN,HAC:HAC,GLYC:GLYC,MMN:MMN,RIP:RIP,IEPOX:IEPOX,MAP:MAP,API:API,LIM:LIM,ONIT:ONIT,APIE:APIE,BPIE:BPIE,LIME:LIME,$            
            PBL:PBLheight,PBLlevel:PBLlevel,box:boxheight}

save,GC,filename=dir+input+'.'+strtrim(string(dates(idate)),2)+'.sav'
print,input+'.'+strtrim(string(dates(idate)),2)+'.sav'
undefine, DataInfo1
ctm_cleanup
;for dates
endfor
end

pro mergecurtain
;I copy this code mainly from merge_cat_file.pro
input1='fast'
;input2='Amazon'
input2='midlat'
StructName='GC'
dir='~/INTEXA/run.std2011.cents.'+input1+'/'
;dates=['20040706','20040707',$
;       '20040708','20040709','20040710','20040711','20040712','20040713','20040714']
;dates=['20040701','20040702','20040703','20040704','20040705','20040706','20040707','20040708','20040709',$
;'20040710','20040711','20040712','20040713','20040714']
dates=['20040701','20040702','20040703','20040704','20040705','20040706','20040707',$
      '20040708','20040709','20040710','20040711','20040712','20040713','20040714',$
      '20040715','20040716','20040717','20040718','20040719','20040720','20040721',$
      '20040722','20040723','20040724','20040725','20040726','20040727','20040728',$
      '20040728','20040729','20040730']
   StructName = StrTrim( StructName, 2 )

   FIRST = 1L

   ; Make list of files for Houston flights
   FOR idate=0L, n_elements( dates )-1L DO BEGIN

;       restore,dir+'Amazon'+input1+'.'+strtrim(string(dates(idate)),2)+'.sav'
      restore,dir+input2+input1+'.'+strtrim(string(dates(idate)),2)+'.sav' 
      cmd = 'tagNames = tag_names( ' + StructName +' )'
       
       status = Execute( cmd )

       ; Number of tags
       nTags = n_elements( tagNames ) 

       ; Add next this structure to existing structure
       IF (FIRST) THEN BEGIN

           FIRST = 0L

           FOR T=0L, nTags-1L DO BEGIN

              ; Construct command to put all data into an array
               cmd = tagNames[T] + ' = '+StructName+'.'+tagNames[T]

               status = Execute( cmd )

           ENDFOR

       ENDIF ELSE BEGIN

           FOR T=0L, n_elements(tagNames)-1L DO BEGIN

               ; Construct command to concatenate arrays
               cmd = tagNames[T] + ' = [ ' + tagNames[T] + ', '+ $
                 StructName+'.'+tagNames[T] +' ]'

               ; Concatenate
               status = Execute( cmd )
 
           ENDFOR

       ENDELSE

   ENDFOR

   FIRST = 1L
   
   cmd = 'outStruct  = {'

   ; Reconstruct structure with all data
   FOR T=0L, nTags-1L DO BEGIN

       IF (FIRST) THEN BEGIN

           FIRST = 0L 

          cmd = cmd + ''+ tagNames[T] +':'+ tagNames[T]
       
       ENDIF ELSE BEGIN

           cmd = cmd + ',' + tagNames[T] +':'+ tagNames[T]

       ENDELSE       

 
   ENDFOR

   cmd = cmd + '}'
   status= Execute( cmd )
;save output structure 
newGC=outStruct

save,newGC,filename=dir+input2+input1+'.sav'
CLOSE,/ALL
end

pro  shortcurtain
;input1='slow'
;dir='~/INTEXA/run.std2011.cents.'+input1+'/'
;print,dir+input1+'.sav'
;restore,dir+input1+'.sav'
restore,'~/INTEXA/run.std2011.cents.slow/slow.sav'&slow=newGC
restore,'~/INTEXA/run.std2011.cents.fast/fast.sav'&fast=newGC

MaxMin = FLTARR(20,2)
OverRide = 1
; Handle the display....
psf          = 0

IF psf THEN BEGIN & psf_open, 'test.ps' , PORTRAIT=PORT & END
display, psf=psf
scols , 33 , top , bottom , ncols , black , white
IF NOT psf THEN BEGIN & dp & END
IF psf THEN !P.FONT = 0 ELSE !P.FONT = 0
IF psf THEN BEGIN & DEVICE, /HELVETICA & END
Ltk = 1.25 & !X.THICK = Ltk & !Y.THICK = Ltk


; Default character size
Csz         = 1.0
!P.CHARSIZE = 1.2


; Include the colour bar settings...
NDivs       =  4 ; 6
cb_charsize =  1.2 ; 1.175 
yp          = 0.051 
dp          = 0.008 
dx          = 0.015


; Set plot positions...
MULTIPANEL, ROWS=5, COLS=4 , OMARGIN=[0.001,0.01,0.0,0.002] 

; ...default individual margins for this set up... 
imargin = [0.0135,0.06,0.012,0.015]

; margins =  [left,bottom,right,top] 

; ******** IMPORTANT **********
; NOTE FOR 16 - 20 we set :: nTimeSteps   = 4 AND LastHour   = JULDAY( 09 , 20 , 2004 , 00 , 00 )
; NOTE FOR 18 - 21 we set :: nTimeSteps   = 3 AND LastHour   = JULDAY( 09 , 21 , 2004 , 00 , 00 )
; AND RESET FirstHour   = JULDAY( 09 , 18 , 2004 , 00 , 00 )

; Set-up the axes how we want
nTimeSteps   = 3              ; # of DAYs for the time-axis
TimeTag      = 'Time [day/month of 2004]'
DateFmt      = '%d/%n'           ; Format of time-axis %m = jul $n = 07


;PlotPosition = [0.1,0.2,0.9,0.95]  ; Plot position
xInfo        = ''                  ; x-axis title
yInfo        = 'Alt (km)'     ; y-axis title
pInfo        = 'Example plot'      ; plot title
dHt          = 0.5                 ; height resolution for y axis
nhr          = 2                   ; # of minor xtick marks (hourly res')
DryFlxCol    = 3                   ; Colour of dry deposition fluxes
PBLCol       = black               ; Colour of PBL height


; Convert boxheights & PBL to to km
Box = newGC.Box / 1000.0 & PBL = newGC.PBL / 1000.0

; Convert the tau format to julian times...
JD         = tau2yymmdd( newGC.Tau ) & jdIN = JD
TimeStamp  = JULDAY( JD.month , JD.day , JD.year , JD.hour , JD.minute )

; Create extra hour at the start to TIMESTAMP ...
FirstHour  = JULDAY( JD.month[0] , JD.day[0] , JD.year[0] , JD.hour[0] - 1 , JD.minute[0] )

; For 18 to 21 plots
;FirstHour   = JULDAY( 07 , 10 , 2004 , 00 , 00 )

; Create final date if we want
LastHour   = JULDAY( 07 , 13 , 2004 , 00 , 00 )

; Add in a first extra hour so we can plot boxes from 't' to 't + dt'...
TimeStamp  = [ FirstHour , TimeStamp ]    ; Extended time values
nTime      = N_ELEMENTS( newGC.Tau )           ; Original # of timesteps (w/o extra hour)
dt         = TimeStamp[1] - TimeStamp[0]  ; Delta time step
eTime      = N_ELEMENTS( TimeStamp )      ; New # of timesteps (with extra hour)

; Save these times to hold the axis in 'GMT' time (even though times are local)...
SafeTime = TimeStamp

; Next convert to local time for this location's longitude...
Longitude = -86 & UTC = TimeStamp
local_time, Longitude, UTC , LocalTime

; & Rename them back...
TimeStamp = LocalTime

; Caculate the pressure/box levels and time..
nlevs = 17 
BoxEdge   = FLTARR( N_ELEMENTS( newGC.Tau ) , nlevs + 1 )

FOR t = 0L , N_ELEMENTS( newGC.Tau ) - 1L DO BEGIN
    FOR l = 1, nlevs DO BEGIN
        BoxEdge[t,l] = TOTAL( Box[t,0:l-1] )
    ENDFOR
ENDFOR

; Trim the top level off - looks nicer for plotting...
CutOff           =  FLOAT( ROUND( BoxEdge[0, nlevs] ) )
BoxEdge[*,nlevs] = CutOff
YR               = [ 0.0 , CutOff ]

; Define the (adjusted) vertical grid boxes...
CellEdge = BoxEdge

; Most important bit really for handling the local times:
; Clip the data from midnight to the endtime of the data
Clipped = [ FirstHour , 0.0 , TimeStamp[eTime-1] , CutOff ] 
Clipped = [ FirstHour , 0.0 , LastHour , CutOff ] 
;Clipped = [ StartHour , 0.0 , LastHour , CutOff ] 

; Order here will dictate order of plots as they appear in array...
;Species = [ 'NO'   , 'NO2' , 'NOy'  , 'N2O5' , $
;            'OH'   , 'POH' , 'O3'   , 'HNO3' , $
;            'ISOP' , 'MVK' , 'MACR' , 'HNO4' , $
;            'HCHO' , 'PAN' , 'PMN'  , 'H2O2' , $
;            'CO'   , 'MP'  , 'ALD'  , 'RCHO' ] 

Species=['Ox','PAN','MVK','MACR','NOx','HCHO','HNO4','H2O2','HNO3','ISOP']
PlotCounter = 0L

FOR i = 0 , N_ELEMENTS(Species) - 1  DO BEGIN

    pInfo = ''
    scols , 33 , top , bottom , ncols , black , white

    CASE Species[i]  OF 
    'NO'  : BEGIN 
         cb_title = 'NO [ppt]'                    & GData   = newgc.no    
         zeroed   = 2                             & Scaling = 1.0E03  
         scols , 33 , top , bottom , ncols , black , white, LOG=1
         ;IF OverRide THEN BEGIN & Maxv = 200.0 & Minv = 0.0 & NDivs = 4& END
         IF OverRide THEN BEGIN & Maxv = 100.0 & Minv = 0.0 & NDivs = 4& END

         END
    'NO2' : BEGIN 
         cb_title = 'NO!D2!N [ppb]'               & GData   = newgc.no2   
         zeroed   = 2                             & Scaling = 1.0 
         scols , 33 , top , bottom , ncols , black , white, LOG=1
         ;IF OverRide THEN BEGIN & Maxv = 1.5 & Minv = 0.0 & NDivs = 3 & END
         IF OverRide THEN BEGIN & Maxv = 0.6 & Minv = 0.0 & NDivs = 3 & END

         END
    'NOx' : BEGIN 
         cb_title = 'NO!Dx!N [ppb]'               & GData   = slow.nox-fast.nox   
         zeroed   = 2                             & Scaling = 1.0 
         scols , 33 , top , bottom , ncols , black , white, LOG=1
         ;IF OverRide THEN BEGIN & Maxv = 1.5 & Minv = 0.0 & NDivs = 3 & END
         IF OverRide THEN BEGIN & Maxv = 0.1 & Minv = 0.0 & NDivs = 4 & END

         END
    'NOy' : BEGIN 
         cb_title = 'NO!Dy!N [ppb]'               & GData   = newgc.noy 
         zeroed   = 1                             & Scaling = 1.0 
         scols , 33 , top , bottom , ncols , black , white, LOG=1
         IF OverRide THEN BEGIN & Maxv = 0.5 & Minv = 0.0 & NDivs = 3 & END
         END
    'N2O5'  : BEGIN 
         cb_title = 'N!D2!NO!D5!N [ppt]'          & GData   = newgc.n2o5   
         zeroed   = 1                             & Scaling = 1.0E03 
         scols , 33 , top , bottom , ncols , black , white, LOG=1
         IF OverRide THEN BEGIN & Maxv = 3.0 & Minv = 0.0 & NDivs = 3 & END
         END
    'HNO3' : BEGIN 
         cb_title = 'HNO!D3!N [ppb]'              & GData   = newgc.hno3 
         zeroed   = 1                             & Scaling = 1.0
         IF OverRide THEN BEGIN & Maxv = 3.0 & Minv = 0.0 & NDivs = 3 & END
         END
    'HNO4'  : BEGIN 
         cb_title = 'HNO!D4!N [ppt]'              & GData   = newgc.hno4 
         zeroed   = 1                             & Scaling = 1.0E3
         IF OverRide THEN BEGIN & Maxv = 4.0 & Minv = 0.0 & NDivs = 4 & END
         END
    'PAN' : BEGIN 
         cb_title = 'PAN [ppt]'                   & GData   = newgc.pan
         zeroed   = 1                             & Scaling = 1.0E3
         IF OverRide THEN BEGIN & Maxv = 100.0 & Minv = 0.0 & NDivs = 3 & END
         END
    'PMN'  : BEGIN 
         cb_title = 'PMN [ppt]'                   & GData   = newgc.pmn
         zeroed   = 1                             & Scaling = 1.0E3
         scols , 33 , top , bottom , ncols , black , white, LOG=1
         IF OverRide THEN BEGIN & Maxv = 100.0 & Minv = 0.0 & NDivs = 4 & END
         END
    'OH' : BEGIN 
         cb_title = 'OH [x10!U6!N molecules cm!U-3!N]'  & GData   = newgc.oh 
         zeroed   = 2                                   & Scaling = 1.0E-6
         scols , 33 , top , bottom , ncols , black , white, LOG=1
         IF OverRide THEN BEGIN & Maxv = 12.0 & Minv = 0.0 & NDivs = 4 & END
         END
    'POH' : BEGIN 
         cb_title = 'P(OH) [x10!U-6!N s!U-1!N]'          & GData   = newgc.poh
         zeroed   = 1                                    & Scaling = 1.0E6
         scols , 33 , top , bottom , ncols , black , white, LOG=1
         IF OverRide THEN BEGIN & Maxv = 8 & Minv = 0.0 & NDivs = 4 & END
         END
    'Ox' : BEGIN 
         cb_title = 'Ozone [ppb]'                 & GData   = newgc.ox
         zeroed   = 1                             & Scaling = 1.0
         IF OverRide THEN BEGIN & Maxv = 90.0 & Minv = 0.0 & NDivs = 3 & END
         END
    'H2O2'  : BEGIN 
         cb_title = 'H!D2!NO!D2!N [ppb]'          & GData   = newgc.h2o2 
         zeroed   = 1                             & Scaling = 1.0
         IF OverRide THEN BEGIN & Maxv = 15.0 & Minv = 0.0 & NDivs = 3 & END
         END
    'ISOP' : BEGIN 
         cb_title = 'Isoprene [ppb]'               & GData   = newgc.isop 
         zeroed   = 1                              & Scaling = 1.0
         IF OverRide THEN BEGIN & Maxv = 2.0 & Minv = 0.0 & NDivs = 4 & END
         END
    'MVK' : BEGIN 
         cb_title = 'MVK [ppb]'                    & GData   = newgc.mvk
         zeroed   = 1                              & Scaling = 1.0
         ;IF OverRide THEN BEGIN & Maxv = 2.0 & Minv = 0.0 & NDivs = 4 & END
         IF OverRide THEN BEGIN & Maxv = 1.0 & Minv = 0.0 & NDivs = 5 & END
         END
    'MACR' : BEGIN 
         cb_title = 'MACR [ppb]'                   & GData   = newgc.macr
         zeroed   = 1                              & Scaling = 1.0
         ;IF OverRide THEN BEGIN & Maxv = 2.0 & Minv = 0.0 & NDivs = 4 & END
         IF OverRide THEN BEGIN & Maxv = 1.0 & Minv = 0.0 & NDivs = 5 & END
         END
    'HCHO' : BEGIN 
         cb_title = 'HCHO [ppb]'                   & GData   = newgc.ch2o
         zeroed   = 1                              & Scaling = 1.0
         IF OverRide THEN BEGIN & Maxv = 8.0 & Minv = 0.0 & NDivs = 4 & END
         END
    'CO' : BEGIN 
         cb_title = 'CO [ppb]'                    & GData   = newgc.co 
         zeroed   = 0                             & Scaling = 1.0
         IF OverRide THEN BEGIN & Maxv = 240.0 & Minv = 80.0 & NDivs = 4 & END
         END
    'ALD' : BEGIN 
         cb_title = 'Acetaldehyde [ppb]'          & GData   = newgc.ald2
         zeroed   = 1                             & Scaling = 1.0
         IF OverRide THEN BEGIN & Maxv = 1.5 & Minv = 0.0 & NDivs = 3 & END
         END
    'ALK4' : BEGIN 
         cb_title = 'ALK4 (C4,5 alkanes) [ppb]'   & GData   = newgc.alk4
         zeroed   = 1                             & Scaling = 1.0
         END
    'C2H6' : BEGIN 
         cb_title = 'Ethane (C!D2!NH!d6!N) [ppb]' & GData   = newgc.c2h6
         zeroed   = 1                             & Scaling = 1.0
         END
    'C3H8' : BEGIN 
         cb_title = 'Propane (C!D3!NH!D8!N) [ppb]' & GData   = newgc.c3h8
         zeroed   = 1                              & Scaling = 1.0
         END
    'MEK' : BEGIN 
         cb_title = 'MEK (>C3 ketones) [ppb]'      & GData   = newgc.mek
         zeroed   = 1                              & Scaling = 1.0
         END
    'MP' : BEGIN 
         cb_title = 'MP (Methyl hydroperoxide) [ppb]'  & GData   = newgc.mp
         zeroed   = 1                                  & Scaling = 1.0
         IF OverRide THEN BEGIN & Maxv = 4.0 & Minv = 0.0 & NDivs = 4 & END
         END
    'RCHO' : BEGIN 
         cb_title = 'RCHO (>C2 aldehydes) [ppb]'   & GData   = newgc.rcho
         zeroed   = 1                              & Scaling = 1.0
         IF OverRide THEN BEGIN & Maxv = 0.2 & Minv = 0.0 & NDivs = 2 & END
         END

        ; Default case:
        ELSE :  BEGIN & zeroed = 0 & Scaling = 1.0 & END

    ENDCASE

    ; 
    MULTIPANEL, POSITION=PlotPosition, MARGIN=imargin
    
    ; ==============================================================

    
    ; Scale data accordingly...
    GData = GData * Scaling

    IF NOT (OverRide) THEN BEGIN
        CASE zeroed OF
            -1 : BEGIN 
                Maxv = MAX( ABS(GData) ,/NAN ) 
                Maxv = ROUND( Maxv )        & Minv = Maxv * (-1.0)     
            END
            0 : BEGIN 
                Maxv = MAX(GData,/NAN)      & Minv = MIN(GData,/NAN)
            END
            1 : BEGIN 
                Maxv = MAX(GData,/NAN)*0.85 & Minv = 0.0 
            END
            2 : BEGIN 
                Maxv = MAX(GData,/NAN)*0.6 & Minv = 0.0 
            END
            ELSE:
        ENDCASE 

        ;PRINT, Species[i] , Maxv , Minv, FORMAT='(A10,2X,F10.2,2X,F10.2)'

    ENDIF

    MaxMin[i,0] = Minv & MaxMin[i,1] = Maxv


    ; Handle the style of the Y-axes
    IF ( PlotCounter EQ 0 OR PlotCounter EQ 4  OR $
         PlotCounter EQ 8 OR PlotCounter EQ 12 OR PlotCounter EQ 16 ) THEN Y_STY=9 ELSE Y_STY=5

    ; Set up the plot axes (use SafeTime to hold axis)...
    plot_setup , nTimeSteps , YR , DateFmt , SafeTime , $
      PlotPosition , xInfo , yInfo , pInfo, Ltk , dHt, black, nHr, Y_STY=Y_STY

    ; Assign colours to the data...
    gcols =  get_gcols( nTime , nlevs , GData , TOP=top , $
                    BOTTOM=Bottom, MAXV=Maxv , MINV=Minv )

    ; Plot the data...
    plot_evol , CellEdge , nTime , nlevs , TimeStamp , gcols, Clipped=Clipped

    ; Overplot the planetary boundary layer...
    OPLOT, TimeStamp , PBL, COLOR=PBLCol, THICK=PBL_Thick

    ; Redraw 'axes' around plot....
    Rectangle, PlotPosition, XVec, YVec
    PLOTS, xvec,yvec, COLOR=black, /NORMAL, THICK=Ltk

    ; For a null color bar
    IF ( Maxv EQ 0.0 AND MinV EQ 0.0 ) THEN BEGIN
        Maxv = 1.0 & Minv = 0.0
    ENDIF 

    ; Add the appropriate colour bar...    
    cb_minv     = Minv & cb_maxv     = Maxv ;* 0.4
    cb_position = [ PlotPosition[0]+dx , PlotPosition[1] - yp ,  $
                    PlotPosition[2]-dx , PlotPosition[1] - yp+dp ]
    
    MY_COLORBAR, BOTTOM=4, CHARSIZE=cb_charsize, COLOR=black,            $
             DIVISIONS=NDivs, FORMAT='(F6.1)', POSITION=cb_position, $
             NCOLORS=ncols, TITLE=cb_title,  MINOR=10,               $
	     RANGE=[cb_minv,cb_maxv], TICKLEN=ticklen, _EXTRA=extra, $
             CHARTHICK=1.0,FONT=font,XTHICK=1,YTHICK=1

    IF ( i EQ 12 OR i EQ 24 OR i EQ 36 OR i EQ 48 OR i EQ 49 ) THEN BEGIN
        xwp = 0.5 & ywp = 0.998 ; 0.99
        XYOUTS, xwp, ywp, MainTitle , COLOR=black, /NORMAL, $
                ALIGNMENT=0.5, CHARSIZE=1.75, CHARTHICK=8
    ENDIF 

    PlotCounter = PlotCounter + 1L
    MULTIPANEL, /ADVANCE, /NOERASE
ENDFOR

IF psf THEN BEGIN & psf_close & END

end

pro  longcurtain
;for this code,the color is controlled by tropics and midlat.

;restore,'~/INTEXA/run.std2011.cents.slow/slow.sav'&slow=newGC
restore,'~/INTEXA/run.std2011.cents.fast/midlatfast.sav'&fast=newGC

MaxMin = FLTARR(20,2)
OverRide = 1
tropics=0
midlat=1
; Handle the display....
psf          = 0

IF psf THEN BEGIN & psf_open, 'test.ps' , PORTRAIT=PORT & END
display, psf=psf
scols , 33 , top , bottom , ncols , black , white
IF NOT psf THEN BEGIN & dp & END
IF psf THEN !P.FONT = 0 ELSE !P.FONT = 0
IF psf THEN BEGIN & DEVICE, /HELVETICA & END
Ltk = 1.25 & !X.THICK = Ltk & !Y.THICK = Ltk


; Default character size
Csz         = 1.0
!P.CHARSIZE = 1.2


; Include the colour bar settings...
NDivs       =  4 ; 6
cb_charsize =  1.2 ; 1.175 
yp          = 0.051 
dp          = 0.008 
dx          = 0.015


; Set plot positions...
MULTIPANEL, ROWS=8, COLS=1 , OMARGIN=[0.1,0.01,0.0,0.002] 

; ...default individual margins for this set up... 
imargin = [0.0135,0.06,0.012,0.015]

; Set-up the axes how we want
nTimeSteps   = 9              ; # of DAYs for the time-axis
TimeTag      = 'Time [day/month of 2004]'
DateFmt      = '%d/%n'           ; Format of time-axis %m = jul $n = 07


;PlotPosition = [0.1,0.2,0.9,0.95]  ; Plot position
xInfo        = ''                  ; x-axis title
yInfo        = 'Alt (km)'     ; y-axis title
pInfo        = 'Example plot'      ; plot title
dHt          = 0.5                 ; height resolution for y axis
nhr          = 2                   ; # of minor xtick marks (hourly res')
DryFlxCol    = 3                   ; Colour of dry deposition fluxes
PBLCol       = black               ; Colour of PBL height


; Convert boxheights & PBL to to km
Box = newGC.Box / 1000.0 & PBL = newGC.PBL / 1000.0&ISOPemis=newGC.ISOPemis/max(newGC.ISOPemis)*3

; Convert the tau format to julian times...
JD         = tau2yymmdd( newGC.Tau ) & jdIN = JD
TimeStamp  = JULDAY( JD.month , JD.day , JD.year , JD.hour , JD.minute )

; Create extra hour at the start to TIMESTAMP ...
FirstHour  = JULDAY( JD.month[0] , JD.day[0] , JD.year[0] , JD.hour[0] - 1 , JD.minute[0] )
nTime      = N_ELEMENTS( newGC.Tau )           ; Original # of timesteps (w/o extra hour)
; Create final date if we want
LastHour   = JULDAY( JD.month[ntime-1], JD.day[ntime-1], JD.year[ntime-1], JD.hour[ntime-1], JD.minute[ntime-1] )

; Add in a first extra hour so we can plot boxes from 't' to 't + dt'...
TimeStamp  = [ FirstHour , TimeStamp ]    ; Extended time values
dt         = TimeStamp[1] - TimeStamp[0]  ; Delta time step
eTime      = N_ELEMENTS( TimeStamp )      ; New # of timesteps (with extra hour)

; Save these times to hold the axis in 'GMT' time (even though times are local)...
SafeTime = TimeStamp

; Next convert to local time for this location's longitude...
Longitude = -86 & UTC = TimeStamp
local_time, Longitude, UTC , LocalTime

; & Rename them back...
TimeStamp = LocalTime

; Caculate the pressure/box levels and time..
nlevs = 17 
BoxEdge   = FLTARR( N_ELEMENTS( newGC.Tau ) , nlevs + 1 )

FOR t = 0L , N_ELEMENTS( newGC.Tau ) - 1L DO BEGIN
    FOR l = 1, nlevs DO BEGIN
        BoxEdge[t,l] = TOTAL( Box[t,0:l-1] )
    ENDFOR
ENDFOR

; Trim the top level off - looks nicer for plotting...
CutOff           =  FLOAT( ROUND( BoxEdge[0, nlevs] ) )
BoxEdge[*,nlevs] = CutOff
YR               = [ 0.0 , CutOff ]

; Define the (adjusted) vertical grid boxes...
CellEdge = BoxEdge

; Most important bit really for handling the local times:
; Clip the data from midnight to the endtime of the data
Clipped = [ FirstHour , 0.0 , TimeStamp[eTime-1] , CutOff ] 
Clipped = [ FirstHour , 0.0 , LastHour , CutOff ] 
;Clipped = [ StartHour , 0.0 , LastHour , CutOff ] 

; Order here will dictate order of plots as they appear in array...
;Species = [ 'NO'   , 'NO2' , 'NOy'  , 'N2O5' , $
;            'OH'   , 'POH' , 'O3'   , 'HNO3' , $
;            'ISOP' , 'MVK' , 'MACR' , 'HNO4' , $
;            'HCHO' , 'PAN' , 'PMN'  , 'H2O2' , $
;            'CO'   , 'MP'  , 'ALD'  , 'RCHO' ] 

Species=['O3','ISOP','MVK','N2O5','NOx','HNO3','PAN','CH2O']
PlotCounter = 0L

FOR i = 0 , N_ELEMENTS(Species) - 1  DO BEGIN

    pInfo = ''
    scols , 33 , top , bottom , ncols , black , white

    CASE Species[i]  OF 
    'NO'  : BEGIN 
         cb_title = 'NO [ppt]'                    & GData   = newgc.no    
         zeroed   = 2                             & Scaling = 1.0E03  
         scols , 33 , top , bottom , ncols , black , white, LOG=1
         ;IF OverRide THEN BEGIN & Maxv = 200.0 & Minv = 0.0 & NDivs = 4& END
         IF OverRide THEN BEGIN & Maxv = 100.0 & Minv = 0.0 & NDivs = 4& END

         END
    'NO2' : BEGIN 
         cb_title = 'NO!D2!N [ppb]'               & GData   = newgc.no2   
         zeroed   = 2                             & Scaling = 1.0 
         scols , 33 , top , bottom , ncols , black , white, LOG=1
         ;IF OverRide THEN BEGIN & Maxv = 1.5 & Minv = 0.0 & NDivs = 3 & END
         IF OverRide THEN BEGIN & Maxv = 0.6 & Minv = 0.0 & NDivs = 3 & END

         END
    'NOx' : BEGIN 
         cb_title = 'NO!Dx!N [ppb]'               & GData   = newgc.nox   
         zeroed   = 2                             & Scaling = 1.0 
         scols , 33 , top , bottom , ncols , black , white, LOG=1
         ;IF OverRide THEN BEGIN & Maxv = 1.5 & Minv = 0.0 & NDivs = 3 & END
         IF tropics THEN BEGIN & Maxv =  1 & Minv = 0.0 & NDivs = 4 & END
         IF midlat THEN BEGIN & Maxv =   5 & Minv = 0.0 & NDivs = 5 & END
         END
    'NOy' : BEGIN 
         cb_title = 'NO!Dy!N [ppb]'               & GData   = newgc.noy 
         zeroed   = 1                             & Scaling = 1.0 
         scols , 33 , top , bottom , ncols , black , white, LOG=1
         IF tropics THEN BEGIN & Maxv = 1 & Minv = 0.0 & NDivs = 4 & END
         IF midlat THEN BEGIN & Maxv = 5 & Minv = 0.0 & NDivs = 4 & END
         END
    'N2O5'  : BEGIN 
         cb_title = 'N!D2!NO!D5!N [ppt]'          & GData   = newgc.n2o5   
         zeroed   = 1                             & Scaling = 1.0E03 
         scols , 33 , top , bottom , ncols , black , white, LOG=1
         IF OverRide THEN BEGIN & Maxv = 100.0 & Minv = 0.0 & NDivs = 5 & END
         END
    'HNO3' : BEGIN 
         cb_title = 'HNO!D3!N [ppb]'              & GData   = newgc.hno3 
         zeroed   = 1                             & Scaling = 1.0
         IF tropics THEN BEGIN & Maxv = 0.2 & Minv = 0.0 & NDivs = 4 & END
         IF midlat THEN BEGIN & Maxv = 3.0 & Minv = 0.0 & NDivs = 3 & END
         END
    'HNO4'  : BEGIN 
         cb_title = 'HNO!D4!N [ppt]'              & GData   = newgc.hno4 
         zeroed   = 1                             & Scaling = 1.0E3
         IF OverRide THEN BEGIN & Maxv = 4.0 & Minv = 0.0 & NDivs = 4 & END
         END
    'PAN' : BEGIN 
         cb_title = 'PAN [ppt]'                   & GData   = newgc.pan
         zeroed   = 1                             & Scaling = 1.0E3
         IF tropics THEN BEGIN & Maxv = 100.0 & Minv = 0.0 & NDivs = 5 & END
         IF midlat THEN BEGIN & Maxv = 1000.0 & Minv = 0.0 & NDivs = 5 & END
         END
    'PMN'  : BEGIN 
         cb_title = 'PMN [ppt]'                   & GData   = newgc.pmn
         zeroed   = 1                             & Scaling = 1.0E3
         scols , 33 , top , bottom , ncols , black , white, LOG=1
         IF OverRide THEN BEGIN & Maxv = 100.0 & Minv = 0.0 & NDivs = 4 & END
         END
    'OH' : BEGIN 
         cb_title = 'OH [x10!U6!N molecules cm!U-3!N]'  & GData   = newgc.oh 
         zeroed   = 2                                   & Scaling = 1.0E-6
         scols , 33 , top , bottom , ncols , black , white, LOG=1
         IF OverRide THEN BEGIN & Maxv = 12.0 & Minv = 0.0 & NDivs = 4 & END
         END
    'POH' : BEGIN 
         cb_title = 'P(OH) [x10!U-6!N s!U-1!N]'          & GData   = newgc.poh
         zeroed   = 1                                    & Scaling = 1.0E6
         scols , 33 , top , bottom , ncols , black , white, LOG=1
         IF OverRide THEN BEGIN & Maxv = 8 & Minv = 0.0 & NDivs = 4 & END
         END
    'O3' : BEGIN 
         cb_title = 'Ozone [ppb]'                 & GData   = newgc.o3
         zeroed   = 1                             & Scaling = 1.0
         IF midlat THEN BEGIN & Maxv = 90.0 & Minv = 30.0 & NDivs = 6 & END
         IF tropics THEN BEGIN & Maxv = 40.0 & Minv = 0.0 & NDivs = 3 & END
         END
    'H2O2'  : BEGIN 
         cb_title = 'H!D2!NO!D2!N [ppb]'          & GData   = newgc.h2o2 
         zeroed   = 1                             & Scaling = 1.0
         IF OverRide THEN BEGIN & Maxv = 15.0 & Minv = 0.0 & NDivs = 3 & END
         END
    'ISOP' : BEGIN 
         cb_title = 'Isoprene [ppb]'               & GData   = newgc.isop/5
         zeroed   = 1                              & Scaling = 1.0
         IF midlat THEN BEGIN & Maxv = 2.0 & Minv = 0.0 & NDivs = 4 & END
         IF tropics THEN BEGIN & Maxv = 10.0 & Minv = 0.0 & NDivs = 5 & END
         END
    'MVK' : BEGIN 
         cb_title = 'MVK [ppb]'                    & GData   = newgc.mvk
         zeroed   = 1                              & Scaling = 1.0
         IF tropics THEN BEGIN & Maxv = 2.0 & Minv = 0.0 & NDivs = 4 & END
         IF midlat THEN BEGIN & Maxv = 1.0 & Minv = 0.0 & NDivs = 5 & END
         END
    'MACR' : BEGIN 
         cb_title = 'MACR [ppb]'                   & GData   = newgc.macr
         zeroed   = 1                              & Scaling = 1.0
         IF tropics THEN BEGIN & Maxv = 2.0 & Minv = 0.0 & NDivs = 4 & END
         IF midlat THEN BEGIN & Maxv = 1.0 & Minv = 0.0 & NDivs = 5 & END
         END
    'CH2O' : BEGIN 
         cb_title = 'HCHO [ppb]'                   & GData   = newgc.ch2o
         zeroed   = 1                              & Scaling = 1.0
         IF OverRide THEN BEGIN & Maxv = 6.0 & Minv = 0.0 & NDivs = 6 & END
         END
    'CO' : BEGIN 
         cb_title = 'CO [ppb]'                    & GData   = newgc.co 
         zeroed   = 0                             & Scaling = 1.0
         IF OverRide THEN BEGIN & Maxv = 240.0 & Minv = 80.0 & NDivs = 4 & END
         END
    'ALD' : BEGIN 
         cb_title = 'Acetaldehyde [ppb]'          & GData   = newgc.ald2
         zeroed   = 1                             & Scaling = 1.0
         IF OverRide THEN BEGIN & Maxv = 1.5 & Minv = 0.0 & NDivs = 3 & END
         END
    'ALK4' : BEGIN 
         cb_title = 'ALK4 (C4,5 alkanes) [ppb]'   & GData   = newgc.alk4
         zeroed   = 1                             & Scaling = 1.0
         END
    'C2H6' : BEGIN 
         cb_title = 'Ethane (C!D2!NH!d6!N) [ppb]' & GData   = newgc.c2h6
         zeroed   = 1                             & Scaling = 1.0
         END
    'C3H8' : BEGIN 
         cb_title = 'Propane (C!D3!NH!D8!N) [ppb]' & GData   = newgc.c3h8
         zeroed   = 1                              & Scaling = 1.0
         END
    'MEK' : BEGIN 
         cb_title = 'MEK (>C3 ketones) [ppb]'      & GData   = newgc.mek
         zeroed   = 1                              & Scaling = 1.0
         END
    'MP' : BEGIN 
         cb_title = 'MP (Methyl hydroperoxide) [ppb]'  & GData   = newgc.mp
         zeroed   = 1                                  & Scaling = 1.0
         IF OverRide THEN BEGIN & Maxv = 4.0 & Minv = 0.0 & NDivs = 4 & END
         END
    'RCHO' : BEGIN 
         cb_title = 'RCHO (>C2 aldehydes) [ppb]'   & GData   = newgc.rcho
         zeroed   = 1                              & Scaling = 1.0
         IF OverRide THEN BEGIN & Maxv = 0.2 & Minv = 0.0 & NDivs = 2 & END
         END

        ; Default case:
        ELSE :  BEGIN & zeroed = 0 & Scaling = 1.0 & END

    ENDCASE

    ; 
    MULTIPANEL, POSITION=PlotPosition, MARGIN=imargin
    
    ; ==============================================================

    
    ; Scale data accordingly...
    GData = GData * Scaling

    IF NOT (OverRide) THEN BEGIN
        CASE zeroed OF
            -1 : BEGIN 
                Maxv = MAX( ABS(GData) ,/NAN ) 
                Maxv = Maxv         & Minv = Maxv * (-1.0)     
            END
            0 : BEGIN 
                Maxv = MAX(GData,/NAN)      & Minv = MIN(GData,/NAN)
            END
            1 : BEGIN 
                Maxv = MAX(GData,/NAN)*0.85 & Minv = 0.0 
            END
            2 : BEGIN 
                Maxv = MAX(GData,/NAN)*0.6 & Minv = 0.0 
            END
            ELSE:
        ENDCASE 

        ;PRINT, Species[i] , Maxv , Minv, FORMAT='(A10,2X,F10.2,2X,F10.2)'

    ENDIF

    MaxMin[i,0] = Minv & MaxMin[i,1] = Maxv


    ; Handle the style of the Y-axes
;    IF ( PlotCounter EQ 0 OR PlotCounter EQ 4  OR $
;         PlotCounter EQ 8 OR PlotCounter EQ 12 OR PlotCounter EQ 16 ) THEN Y_STY=9 ELSE Y_STY=5

    Y_STY=9
    ; Set up the plot axes (use SafeTime to hold axis)...
    plot_setup , nTimeSteps , YR , DateFmt , SafeTime , $
      PlotPosition , xInfo , yInfo , pInfo, Ltk , dHt, black, nHr, Y_STY=Y_STY

    ; Assign colours to the data...
    gcols =  get_gcols( nTime , nlevs , GData , TOP=top , $
                    BOTTOM=Bottom, MAXV=Maxv , MINV=Minv )

    ; Plot the data...
    plot_evol , CellEdge , nTime , nlevs , TimeStamp , gcols, Clipped=Clipped

    ; Overplot the planetary boundary layer...
    OPLOT, TimeStamp , PBL, COLOR=PBLCol, THICK=PBL_Thick
   ; OPLOT, TimeStamp , ISOPemis, COLOR=1, THICK=PBL_Thick

    ; Redraw 'axes' around plot....
    Rectangle, PlotPosition, XVec, YVec
    PLOTS, xvec,yvec, COLOR=black, /NORMAL, THICK=Ltk

    ; For a null color bar
    IF ( Maxv EQ 0.0 AND MinV EQ 0.0 ) THEN BEGIN
        Maxv = 1.0 & Minv = 0.0
    ENDIF 

    ; Add the appropriate colour bar...    
    cb_minv     = Minv & cb_maxv     = Maxv ;* 0.4
    cb_position = [ PlotPosition[0]+dx , PlotPosition[1] - yp ,  $
                    PlotPosition[2]-dx , PlotPosition[1] - yp+dp ]
    
    MY_COLORBAR, BOTTOM=4, CHARSIZE=cb_charsize, COLOR=black,            $
             DIVISIONS=NDivs, FORMAT='(F7.2)', POSITION=cb_position, $
             NCOLORS=ncols, TITLE=cb_title,  MINOR=10,               $
	     RANGE=[cb_minv,cb_maxv], TICKLEN=ticklen, _EXTRA=extra, $
             CHARTHICK=1.0,FONT=font,XTHICK=1,YTHICK=1

    IF ( i EQ 12 OR i EQ 24 OR i EQ 36 OR i EQ 48 OR i EQ 49 ) THEN BEGIN
        xwp = 0.5 & ywp = 0.998 ; 0.99
        XYOUTS, xwp, ywp, MainTitle , COLOR=black, /NORMAL, $
                ALIGNMENT=0.5, CHARSIZE=1.75, CHARTHICK=8
    ENDIF 

    PlotCounter = PlotCounter + 1L
    MULTIPANEL, /ADVANCE, /NOERASE
ENDFOR

IF psf THEN BEGIN & psf_close & END

MULTIPANEL,/off
end

pro  compcurtain
restore,'~/INTEXA/run.std2011.cents.slow/Amazonslow.sav'&slow=newGC
restore,'~/INTEXA/run.std2011.cents.fast/Amazonfast.sav'&fast=newGC

MaxMin = FLTARR(20,2)
OverRide = 0
; Handle the display....
psf          = 0

IF psf THEN BEGIN & psf_open, 'test.ps' , PORTRAIT=PORT & END
display, psf=psf
;scols , 33 , top , bottom , ncols , black , white
IF NOT psf THEN BEGIN & dp & END
IF psf THEN !P.FONT = 0 ELSE !P.FONT = 0
IF psf THEN BEGIN & DEVICE, /HELVETICA & END
Ltk = 1.25 & !X.THICK = Ltk & !Y.THICK = Ltk


; Default character size
Csz         = 1.0
!P.CHARSIZE = 1.2


; Include the colour bar settings...
NDivs       =  4 ; 6
cb_charsize =  1.2 ; 1.175 
yp          = 0.051 
dp          = 0.008 
dx          = 0.015


; Set plot positions...
MULTIPANEL, ROWS=7, COLS=1 , OMARGIN=[0.1,0.01,0.0,0.002] 

; ...default individual margins for this set up... 
imargin = [0.0135,0.06,0.012,0.015]

; Set-up the axes how we want
nTimeSteps   = 9              ; # of DAYs for the time-axis
TimeTag      = 'Time [day/month of 2004]'
DateFmt      = '%d/%n'           ; Format of time-axis %m = jul $n = 07


;PlotPosition = [0.1,0.2,0.9,0.95]  ; Plot position
xInfo        = ''                  ; x-axis title
yInfo        = 'Alt (km)'     ; y-axis title
pInfo        = 'Example plot'      ; plot title
dHt          = 0.5                 ; height resolution for y axis
nhr          = 2                   ; # of minor xtick marks (hourly res')
DryFlxCol    = 3                   ; Colour of dry deposition fluxes
;PBLCol       = black               ; Colour of PBL height


; Convert boxheights & PBL to to km
Box = newGC.Box / 1000.0 & PBL = newGC.PBL / 1000.0

; Convert the tau format to julian times...
JD         = tau2yymmdd( newGC.Tau ) & jdIN = JD
TimeStamp  = JULDAY( JD.month , JD.day , JD.year , JD.hour , JD.minute )

; Create extra hour at the start to TIMESTAMP ...
FirstHour  = JULDAY( JD.month[0] , JD.day[0] , JD.year[0] , JD.hour[0] - 1 , JD.minute[0] )
nTime      = N_ELEMENTS( newGC.Tau )           ; Original # of timesteps (w/o extra hour)
; Create final date if we want
LastHour   = JULDAY( JD.month[ntime-1], JD.day[ntime-1], JD.year[ntime-1], JD.hour[ntime-1], JD.minute[ntime-1] )

; Add in a first extra hour so we can plot boxes from 't' to 't + dt'...
TimeStamp  = [ FirstHour , TimeStamp ]    ; Extended time values
dt         = TimeStamp[1] - TimeStamp[0]  ; Delta time step
eTime      = N_ELEMENTS( TimeStamp )      ; New # of timesteps (with extra hour)

; Save these times to hold the axis in 'GMT' time (even though times are local)...
SafeTime = TimeStamp

; Next convert to local time for this location's longitude...
Longitude = -86 & UTC = TimeStamp
local_time, Longitude, UTC , LocalTime

; & Rename them back...
TimeStamp = LocalTime

; Caculate the pressure/box levels and time..
nlevs = 17 
BoxEdge   = FLTARR( N_ELEMENTS( newGC.Tau ) , nlevs + 1 )

FOR t = 0L , N_ELEMENTS( newGC.Tau ) - 1L DO BEGIN
    FOR l = 1, nlevs DO BEGIN
        BoxEdge[t,l] = TOTAL( Box[t,0:l-1] )
    ENDFOR
ENDFOR

; Trim the top level off - looks nicer for plotting...
CutOff           =  FLOAT( ROUND( BoxEdge[0, nlevs] ) )
BoxEdge[*,nlevs] = CutOff
YR               = [ 0.0 , CutOff ]

; Define the (adjusted) vertical grid boxes...
CellEdge = BoxEdge

; Most important bit really for handling the local times:
; Clip the data from midnight to the endtime of the data
Clipped = [ FirstHour , 0.0 , TimeStamp[eTime-1] , CutOff ] 
Clipped = [ FirstHour , 0.0 , LastHour , CutOff ] 
;Clipped = [ StartHour , 0.0 , LastHour , CutOff ] 

Species=['Ox','ISOP','MVK','MACR','NOx','PAN','CH2O']
PlotCounter = 0L

FOR i = 0 , N_ELEMENTS(Species) - 1  DO BEGIN

    pInfo = ''
    scols , -1 , top , bottom , ncols , black , white
;   MyCt, 63, NColors=14, /MidCol, /White,/Reverse
    CASE Species[i]  OF 
    'NOx' : BEGIN 
         cb_title = 'NO!Dx!N [ppt]'               & GData   = slow.nox-fast.nox   
         zeroed   = -1                             & Scaling = 1e3 
         ;scols , 33 , top , bottom , ncols , black , white, LOG=1
         ;IF OverRide THEN BEGIN & Maxv = 1.5 & Minv = 0.0 & NDivs = 3 & END
         END
    'PAN' : BEGIN 
         cb_title = 'PAN [ppt]'                   & GData   = slow.pan-fast.pan
         zeroed   = -1                             & Scaling = 1.0E3
         IF OverRide THEN BEGIN & Maxv = 100.0 & Minv = 0.0 & NDivs = 3 & END
         END
    'Ox' : BEGIN 
         cb_title = 'Ozone [ppb]'                 & GData   = slow.ox-fast.ox
         zeroed   = -1                             & Scaling = 1.0
         IF OverRide THEN BEGIN & Maxv = 1.0 & Minv = 0.0 & NDivs = 4 & END
         END
    'ISOP' : BEGIN 
         cb_title = 'Isoprene [ppb]'               & GData   =slow.isop-fast.isop 
         zeroed   = -1                              & Scaling = 1.0
         IF OverRide THEN BEGIN & Maxv = 1.0 & Minv = 0.0 & NDivs = 4 & END
         END
    'MVK' : BEGIN 
         cb_title = 'MVK [ppt]'                    & GData   = slow.mvk-fast.mvk
         zeroed   = -1                              & Scaling = 1e3
         ;IF OverRide THEN BEGIN & Maxv = 2.0 & Minv = 0.0 & NDivs = 4 & END
         IF OverRide THEN BEGIN & Maxv = 1.0 & Minv = 0.0 & NDivs = 5 & END
         END
    'MACR' : BEGIN 
         cb_title = 'MACR [ppt]'                   & GData   = slow.macr-fast.macr
         zeroed   = -1                              & Scaling = 1e3
         ;IF OverRide THEN BEGIN & Maxv = 2.0 & Minv = 0.0 & NDivs = 4 & END
         END
    'CH2O' : BEGIN 
         cb_title = 'HCHO [ppb]'                   & GData   = slow.ch2o-fast.ch2o
         zeroed   = -1                              & Scaling = 1.0
         IF OverRide THEN BEGIN & Maxv = 1.0 & Minv = 0.0 & NDivs = 4 & END
         END
    ENDCASE

    ; 
    MULTIPANEL, POSITION=PlotPosition, MARGIN=imargin
    
    ; ==============================================================

    
    ; Scale data accordingly...
    GData = GData * Scaling

    IF NOT (OverRide) THEN BEGIN
        CASE zeroed OF
            -1 : BEGIN 
                Maxv = MAX( ABS(GData) ,/NAN ) 
                Maxv = Maxv         & Minv = Maxv * (-1.0)     
            END
            0 : BEGIN 
                Maxv = MAX(GData,/NAN)      & Minv = MIN(GData,/NAN)
            END
            1 : BEGIN 
                Maxv = MAX(GData,/NAN)*0.85 & Minv = 0.0 
            END
            2 : BEGIN 
                Maxv = MAX(GData,/NAN)*0.6 & Minv = 0.0 
            END
            ELSE:
        ENDCASE 

        PRINT, Species[i] , Maxv , Minv, FORMAT='(A10,2X,F10.3,2X,F10.3)'

    ENDIF

    MaxMin[i,0] = Minv & MaxMin[i,1] = Maxv


    ; Handle the style of the Y-axes
;    IF ( PlotCounter EQ 0 OR PlotCounter EQ 4  OR $
;         PlotCounter EQ 8 OR PlotCounter EQ 12 OR PlotCounter EQ 16 ) THEN Y_STY=9 ELSE Y_STY=5

    Y_STY=9
    ; Set up the plot axes (use SafeTime to hold axis)...
    plot_setup , nTimeSteps , YR , DateFmt , SafeTime , $
      PlotPosition , xInfo , yInfo , pInfo, Ltk , dHt, black, nHr, Y_STY=Y_STY

    ; Assign colours to the data...
    gcols =  get_gcols( nTime , nlevs , GData , TOP=top , $
                    BOTTOM=Bottom, MAXV=Maxv , MINV=Minv )

    ; Plot the data...
    plot_evol , CellEdge , nTime , nlevs , TimeStamp , gcols, Clipped=Clipped

    ; Overplot the planetary boundary layer...
;    OPLOT, TimeStamp , PBL, COLOR=0, THICK=PBL_Thick

    ; Redraw 'axes' around plot....
    Rectangle, PlotPosition, XVec, YVec
    PLOTS, xvec,yvec, COLOR=0, /NORMAL, THICK=Ltk

    ; For a null color bar
    IF ( Maxv EQ 0.0 AND MinV EQ 0.0 ) THEN BEGIN
        Maxv = 1.0 & Minv = 0.0
    ENDIF 

    ; Add the appropriate colour bar...    
    cb_minv     = Minv & cb_maxv     = Maxv ;* 0.4
    cb_position = [ PlotPosition[0]+dx , PlotPosition[1] - yp ,  $
                    PlotPosition[2]-dx , PlotPosition[1] - yp+dp ]
    
;    MY_COLORBAR, BOTTOM=4, CHARSIZE=cb_charsize, COLOR=black,            $
;             DIVISIONS=NDivs, FORMAT='(F8.3)', POSITION=cb_position, $
;             NCOLORS=ncols, TITLE=cb_title,  MINOR=10,               $
;	     RANGE=[cb_minv,cb_maxv], TICKLEN=ticklen, _EXTRA=extra, $
;             CHARTHICK=1.0,FONT=font,XTHICK=1,YTHICK=1

   ColorBar,                        $
   Max=Maxv,       Min=Minv,       NColors=ncols,     $
   Bottom=Bottom,       Color=black,       Position=CB_Position, $
   Unit=CBUnit,         Divisions=nDivs, Log=Log,             $
   Format='(F8.3)',     Charsize=cb_charsize,      TickLen=tickLen,   $
   C_Colors=CC_Colors,   C_Levels=C_Levels2, Vertical=CBVertical, $
   _EXTRA=e
   stop
   xyouts, plotposition[0]-0.15,plotposition[1]+0.1,cb_title
    IF ( i EQ 12 OR i EQ 24 OR i EQ 36 OR i EQ 48 OR i EQ 49 ) THEN BEGIN
        xwp = 0.5 & ywp = 0.998 ; 0.99
        XYOUTS, xwp, ywp, MainTitle , COLOR=black, /NORMAL, $
                ALIGNMENT=0.5, CHARSIZE=1.75, CHARTHICK=8
    ENDIF 

    PlotCounter = PlotCounter + 1L
    MULTIPANEL, /ADVANCE, /NOERASE
ENDFOR

IF psf THEN BEGIN & psf_close & END

end
