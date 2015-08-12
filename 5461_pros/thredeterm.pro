pro thredeterm,velbegin=velbegin,velend=velend
;+
; Name
;      thredeterm
; Purpose
;      Determine the threshold of real data in the whole field of M101   
; Parameters
;      velbegin - the beginning velocity of the folding; def: 224.779 km/s
;        velend - the ending velocity of the folding; def: 302.257 km/s
; 
; Inherited from "accum_ek_5461.pro"    
; Modified by sw, Jan 20, 2015
;-

  ;Procedure keyword setting
  if keyword_set(velbegin) eq 0 then velbegin=80000.D;224000.D
  if keyword_set(velend)   eq 0 then velend=500000.D;303000.D
  ;Physical paramters
  dis=7.4           ;in Mpc
  bsize=7.49*6.07*!dpi/(4*alog(2)) ;in arcsec^2
  pixelsize=1.0     ;in arcsec
  delv=5.165198242D ;in km/s
  data2flux=pixelsize^2/bsize ;Conversion factor of data (in Jy/beam)
  ; to flux (in Jy/pixel)
  flux2mass=2.36e5*dis^2*delv ;Conversion factor of flux (in Jy/pixel) 
  ; to HI mass (in M_solar)
  data2mass=data2flux*flux2mass*!msolar ;Conversion factor of data 
  ; (in Jy/beam) to HI mass (in gram)

  ;System paramters
  genpath='../thredeterm/' & figpath=genpath+'figures/'
  cubename='m101_ro_cube.fits'
  intename='m101_ro_mom0.fits'

  ;Plot parameters
  set_plot,'ps'
  device,/encapsulated,preview=2,bits=8,/color
  tvlct,[0,255,  0,  0,  0,255,255,  0,  0,127,255], $
        [0,  0,255,  0,255,  0,127,127,255,  0,  0], $
        [0,  0,  0,255,255,255,  0,255,127,255,127]
    ;Black,  R,  G,B,cyan,purp,orag,...
  
  ;Read data
  intdata=mrdfits(genpath+intename,0,hdr,/silent)
  rawdata=mrdfits(genpath+cubename,0,hdr,/silent)
  rawsz=size(rawdata)
  ; backavg_5461,cubename,meaval,stdval,genpath=genpath
  ; for i=0,rawsz[1]-1 do for j=0,rawsz[2]-1 do $
  ;   rawdata[i,j,*]-=meaval
  crval3=sxpar(hdr,'crval3') & crpix3=sxpar(hdr,'crpix3')
  cdelt3=sxpar(hdr,'cdelt3') & delv=-cdelt3
  ;Pre-select data
  velarr=(dindgen(rawsz(3))-crpix3)*cdelt3+crval3
  velselind=where(velarr gt velbegin and velarr lt velend,nvelselind)
  if nvelselind eq 0 then begin
    print,'The velocity range seems not right!'
    return
  endif
  velarr=velarr[velselind]/1e3 & xrange=[min(velarr),max(velarr)]
  data=rawdata[*,*,velselind] & tmpdata=data
  meaint=mean(intdata) & velcur=dblarr(nvelselind)
  ndata=n_elements(intdata)
  
  threlist=[-100.,2.0,2.3,2.4,2.5,2.6,2.7,3.0]*0.48e-3
  nthrelist=n_elements(threlist)
  threstr=strtrim(string(threlist,format='(E9.2)'),2)
  meadev=dblarr(nthrelist) & meastd=dblarr(nthrelist)
  condev=dblarr(nthrelist) & constd=dblarr(nthrelist)
  meaval=dblarr(nthrelist) & conval=dblarr(nthrelist)
  for i=0,nthrelist-1 do begin
    selind=where(tmpdata lt threlist[i],nselind,comp=cselind)
    if nselind ne 0 then tmpdata[selind]=0.0
    tmpimg=total(tmpdata,3)*delv & meaimg=mean(tmpimg)
    meadev[i]=meaimg-meaint & meaval[i]=meaimg
    meastd[i]=mean(abs(tmpimg-intdata))
    for j=0,nvelselind-1 do velcur[j]=total(tmpdata[*,*,j])
    velcur*=data2flux
    ; ;Continous judgement
    ; condata=data &
    ;
    
    ymin=min(velcur,max=ymax) & ydel=ymax-ymin
    yrange=[ymin-ydel*0.1,ymax+ydel*0.1]
    device,filename=figpath+'velcur_'+string(i+1,format='(i1)')+'_thr'+ $
      threstr[i]+'.eps'
    plot,velarr,velcur,xstyle=1,ystyle=1,xrange=xrange,yrange=yrange, $
      xtitle='!6Velocity (km s!U-1!N)',ytitle='Intensity (Jy)', $
      charsize=1.5,thick=3
    device,/close
  endfor
  
  print,'N-sigma: ',threlist/0.48e-3,format='(a9,'+ $
    strtrim(string(nthrelist),2)+'f10.3)'
  print,'Target:  ',1102.
  print,'Values:  ',meaval*ndata*data2flux/1e3,format='(a9,'+ $
    strtrim(string(nthrelist),2)+'f10.3)'
  device,filename=figpath+'MeanValue.eps'
  plot,threlist[1:*],meadev[1:*],psym=-6,thick=3,xtitle='Threshold', $
    ytitle='Mean Values',charsize=1.5
  oplot,threlist,intarr(nthrelist),linestyle=2,thick=3
  device,/close
  device,filename=figpath+'StdValue.eps'
  plot,threlist[1:*],meastd[1:*],psym=-6,thick=3,xtitle='Threshold', $
    ytitle='Standard deviation',charsize=1.5
  device,/close
  
  set_plot,'x'
  
END