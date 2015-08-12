function read_accumek,filename
;+
; Name
;      read_accumek (function)
; Purpose
;      read the kinetic energy information of the accumulation
;      in NGC 5461
;      
; Written by sw, Feb 05, 2015
;-

  rec={rawval:0.0D,subval:0.0D,netval:0.0D,area:0L}
  list=replicate(rec,50)
  openr,lun,filename,/get_lun
  ;-Read record until en-of-file reached
  nrecords=0L & recnum  =1L
  while(eof(lun) ne 1) do begin
    ;-Read this record(jumps to bad_rec:on error)
    on_ioerror,bad_rec
    error=1
    readf,lun,rec,format='(3e10.3,i3)'
    error=0
    ;-Store data for this record
    list[nrecords]=rec
    nrecords=nrecords+1L
    ;- Check for bad imput record
    bad_rec:
    if(error eq 1) then $
      print,'Bad data at record',recnum
    recnum=recnum+1
  endwhile
  ;-Close input file
  free_lun,lun
  ;-Trim data array and return it to caller
  list=list[0:nrecords-1]
  return,list
end

pro accumek_plot
;+
; Name
;      accumek_plot
; Purpose
;      Plot the information of accumulative kinetic energy
; 
; Written by sw, Feb 05, 2015
;-

  genpath='../accum_ek/' & figpath=genpath+'figures/'
  sublist=['enlarge','integr'] & subdir =['wide','narrow']
  labtext=['1105','1098'] & divfac=double(1e30)*1e22
  divstr=strtrim(string(round(alog10(divfac))),2)
  
  ;Plot parameters
  set_plot,'ps'
  device,/encapsulated,preview=2,bits=8,/color
  tvlct,[0,255,  0,  0,  0,255,255,  0,  0,127,255], $
        [0,  0,255,  0,255,  0,127,127,255,  0,  0], $
        [0,  0,  0,255,255,255,  0,255,127,255,127]
    ;Black,  R,  G,B,cyan,purp,orag,...
  charsize=1.5 & thick=3 & lchsize=1.2
  xind=indgen(10)+1 & xrange=[0.5,10.5]
  legtext=['raw','sub','net']
  
  for i=0,1 do begin
    for j=0,1 do begin
      filepath=genpath+sublist[i]+'/'+subdir[j]+'/'
      for k=0,1 do begin
        data=read_accumek(filepath+'ek'+labtext[k]+'_'+ $
          sublist[i]+'_'+subdir[j]+'.info')
        mearaw=mean(data.rawval) & measub=mean(data.subval)
        meanet=mean(data.netval)
        print,'For '+sublist[i]+', '+subdir[j]+', '+'H'+labtext[k]+ $
          ', raw:'+string(mearaw,format='(e10.3)')+ $
          ', sub:'+string(measub,format='(e10.3)')+ $
          ', net:'+string(meanet,format='(e10.3)')
        for m=0,2 do data.(m)/=divfac
        ymin=min([data.rawval,data.subval,data.netval],max=ymax)
        ydel=ymax-ymin & yrange=[ymin-0.1*ydel,ymax+0.2*ydel]
        device,filename=figpath+'ek'+labtext[k]+'_'+ $
          sublist[i]+'_'+subdir[j]+'.eps'
        plot,xrange,yrange,/nodata,xrange=xrange,yrange=yrange, $
          xstyle=1,ystyle=1,charsize=charsize,xtitle='binsize (pixel)', $
          ytitle='Kinetic Energy (10!U'+divstr+'!Nerg s!U-1!N)'
        oplot,xind,data.rawval,thick=thick,color=1,psym=-4
        oplot,xind,data.subval,thick=thick,color=2,psym=-4
        oplot,xind,data.netval,thick=thick,color=3,psym=-4
        legend,legtext,color=[1,2,3],thick=thick+[0,0,0],psym=-4+[0,0,0]
        for m=0,9 do xyouts,xind[m],data[m].subval, $
          strtrim(string(data[m].area),2),charsize=lchsize
        device,/close
      endfor
    endfor
  endfor
  
  set_plot,'x'

END