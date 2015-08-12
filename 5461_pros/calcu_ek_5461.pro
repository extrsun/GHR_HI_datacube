pro calcu_ek_5461,velbegin=velbegin,velend=velend;,subdir=subdir;,binsize=binsize
;+
; Name
;      calcu_ek_5461
; Purpose
;      Calculate the accumulative kinematic energy at different positions 
;      along the spiral arm around NGC5461.     
; Parameters
;      velbegin - the beginning velocity of the folding; def: 224.779 km/s
;        velend - the ending velocity of the folding; def: 302.257 km/s
;        subdir - the name of the sub-directory
;      
; Inherited from "accum_ek_5461.pro", and
; Modified by sw, Feb 06, 2015
;-

  ;Procedure keyword setting
  ; if keyword_set(velbegin) eq 0 then velbegin=80000.D
  ; if keyword_set(velend) eq 0 then velend=500000.D
  if keyword_set(velbegin) eq 0 then velbegin=170000.D;224000.D
  if keyword_set(velend)   eq 0 then velend=320000.D;303000.D
  ; if keyword_set(binsize) eq 0 then binsize=8
  
  ;Data threshold
  threshold=1.152e-3 ;2.4-sigma
  ; threshold=1.44e-3  ;3-sigma
  
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
  divfac_ek=double(1e29)*1e20 ;Dividing factor of kinetic energy
  divfac_ekd=1.0e12 ;Dividing factor of density of kinetic energy
  ; divfac_mass_g=1e40 ;Dividing factor of HI mass (in gram)
  divfac_mass=1.0e4 ;Dividing factor of HI mass (in M_solar)

  ;System paramters
  genpath='../'
  cubename=genpath+'5461_RO_CUBE_rot.fits' ;Use the enlarged HI CUBE
  ;The CUBE file of the other weighting, 'NA', is alternative.

  ; ;Plot parameters
  ; set_plot,'ps'
  ; device,/encapsulated,preview=2,bits=8,/color
  ; tvlct,[0,255,  0,  0,  0,255,255,  0,  0,127,255], $
  ;       [0,  0,255,  0,255,  0,127,127,255,  0,  0], $
  ;       [0,  0,  0,255,255,255,  0,255,127,255,127]
  ;   ;Black,  R,  G,B,cyan,purp,orag,...
  ; charsize=1.5 & thick=5 & ancthick=3 & lchsize=1.2 & ekchsize=1.8
  
  ;Read data
  rawdata=mrdfits(cubename,0,hdr,/silent)
  crval3=sxpar(hdr,'crval3') & crpix3=sxpar(hdr,'crpix3')
  cdelt3=sxpar(hdr,'cdelt3')  
  ;Pre-select data
  rawsz=size(rawdata) & velarr=(dindgen(rawsz(3))-crpix3)*cdelt3+crval3
  velselind=where(velarr gt velbegin and velarr lt velend,nvelselind)
  if nvelselind eq 0 then begin
    print,'The velocity range seems not right!'
    return
  endif
  velvalarr=velarr[velselind]
  
  
  velvalarr*=1e2 ;Convert the velocity from m/s to cm/s (cgs unit)
  data=rawdata[*,*,velselind]
  ;Convert the data into HI mass (in gram)
  data*=data2mass ;in gram/pixel
  
  sublist=['enlarge','rot'] & nsublist=n_elements(sublist)
  ;Use two kinds of velocity map: original ('enlarge') or integrated ('rot')
  for i=0,nsublist-1 do begin    
    filepath=genpath+'accum_ek/'+sublist[i]+'/'
    subpath=filepath+subdir+'/' & figpath=subpath+'figures/'
    ;Open the ascii files used to record the kinetic energy
    openu,lun1,subpath+'ek1105_'+sublist[i]+'_'+subdir+'.info',/get_lun,/append
    openu,lun2,subpath+'ek1098_'+sublist[i]+'_'+subdir+'.info',/get_lun,/append
    ;Calculate the mean velocity at different X-positions
    intdata=mrdfits(filepath+'5461_MOM0_'+sublist[i]+'.fits',0,hdr1,/silent)
    veldata=mrdfits(filepath+'5461_MOM1_'+sublist[i]+'.fits',0,hdr2,/silent)
    selind=where(finite(veldata,/nan),nselind)
    if nselind ne 0 then veldata[selind]=0.0
    intensity=intdata/(delv*1e3) & velocity=veldata
    fdat=intensity[*,ymin:ymax] & vdat=velocity[*,ymin:ymax]
    for k=0,xnum-1 do grdat[k,*]=total(fdat[xindl[k]:xindr[k],*],1)
    grdat=grdat>0 & selint=total(grdat,2)
    for k=0,xnum-1 do begin
      selvdat=vdat[xindl[k]:xindr[k],*]
      selfdat=fdat[xindl[k]:xindr[k],*]
      for l=0,ynum-1 do begin
        sind=where(finite(selvdat[*,l]))
        selgrfdat=total(selfdat[sind,l])
        if selgrfdat gt 0 then begin
          grvdat[k,l]=total(selvdat[sind,l]*selfdat[sind,l])/selgrfdat
        endif else grvdat[k,l]=0.0
        if total(selfdat[sind,l]) eq 0.0 then $
          print,sublist[i],xindl[k],xindr[k],l+ymin;,grvdat[k,l]
      endfor
    endfor
    selvel=total(grvdat*grdat,2)/selint
    ;Convert the velocity from m/s to cm/s (cgs unit)
    grvdat*=1e2 & selvel*=1e2
    ;Calculate the integrated mass profile (in M_solar/pixel)
    grmass=grdat*data2mass/!msolar/binsize
    selmass=selint*data2mass/!msolar/pixarea
  
    ;Calculate the accumulative kinematic energy
    for k=0,xnum-1 do begin
      fdata=data[xindl[k]:xindr[k],*,*]
      selind=where(fdata lt threshold*data2mass,nselind)
      if nselind ne 0 then fdata[selind]=0.0
      for l=0,ynum-1 do cenvel[*,l,*]=selvel[k]
      ekval=total(fdata*(velvalarr-cenvel)^2/2,3)
      grekv1[k,*]=total(ekval,1)
      for l=0,ynum-1 do cenvel[*,l,*]=grvdat[k,l]
      ekval=total(fdata*(velvalarr-cenvel)^2/2,3)
      nulind=where(grvdat[k,*] eq 0,nnulind)
      if nnulind ne 0 then begin
        tmpval=total(fdata[*,nulind,*]*velvalarr[*,nulind,*]^2,3)
        if total(tmpval) ne 0 then print,nnulind,total(tmpval,1)
      endif
      grekv2[k,*]=total(ekval,1)
    endfor
    accumek1=total(grekv1,2)/pixarea
    accumek2=total(grekv2,2)/pixarea
    grekv1/=binsize & grekv2/=binsize
        
  endfor

END