pro integr_5461,velbegin=velbegin,velend=velend
;+
; Name
;      integr_5461
; Purpose
;      Calculate the integral results of the HI data cube toward NGC 5461, 
;      which includes intensity (MOM0), velocity field (MOM1), and dispersion 
;      distribution (MOM2).
; Parameters
;      velbegin - the beginning velocity of the folding; def: 224.779 km/s
;        velend - the ending velocity of the folding; def: 302.257 km/s
;      
; Written by sw, Jan 06, 2015
;-

  ;Procedure keyword setting
  ; if keyword_set(velbegin) eq 0 then velbegin=204000.D
  ; if keyword_set(velend)   eq 0 then velend=308000.D
  if keyword_set(velbegin) eq 0 then velbegin=80000.D
  if keyword_set(velend) eq 0 then velend=500000.D

  ;Data threshold
  threshold=1.152e-3 ;2.4-sigma
  ; threshold=1.44e-3  ;3-sigma
  
  ;System parameters
  filepath='../' & fileapp='RO'
  figpath=filepath+'figure/'
  cubeapp=['enlarge','small' ,'rot'   ,'contisel']
  outapp =['enlint' ,'smlint','rotint','conint'  ]
  nfile=n_elements(cubeapp)
  for k=0,nfile-1 do begin
    cubename='5461_'+fileapp+'_CUBE_'+cubeapp[k]+'.fits' 
    ;Use the enlarged HI CUBE; the CUBE file of the other 
    ;weighting, 'NA', is alternative.
  
    ;Read data
    rawdata=mrdfits(filepath+cubename,0,hdr)
    crval3=sxpar(hdr,'crval3') & crpix3=sxpar(hdr,'crpix3')
    cdelt3=sxpar(hdr,'cdelt3')
    ;Pre-select data
    rawsz=size(rawdata) & velarr=(dindgen(rawsz(3))-crpix3)*cdelt3+crval3
    velselind=where(velarr gt velbegin and velarr lt velend,nvelselind)
    if nvelselind gt 0 then begin
      data=rawdata[*,*,velselind]
    endif else begin
      print,'The velocity range seems not right!'
      return
    endelse
    velo=ulindgen(nvelselind,rawsz(1),rawsz(2)) mod (nvelselind)
    velo=double(transpose(velo))*cdelt3+velarr[velselind[0]]
    routind=where(data lt threshold)
    data[routind]=0.0 & grdat=total(data,3)
    selind=where(grdat eq 0,nselind)
    if nselind ne 0 then grdat[selind]=1.0
    velocity=total(data*velo,3)/grdat
    velc=dblarr(rawsz(1),rawsz(2),nvelselind)
    for i=0,nvelselind-1 do velc[*,*,i]=velocity
    dispersion=sqrt(total(data*(velo-velc)^2,3)/grdat)
    if nselind ne 0 then begin
      velocity[selind]=0.0 & dispersion[selind]=0.0
      grdat[selind]=0.0
    endif
  
    ;Modify the header of the fits cube to be used for the fits image
    sxdelpar,hdr,'NAXIS3' & sxdelpar,hdr,'NAXIS4'
    sxdelpar,hdr,'CRVAL3' & sxdelpar,hdr,'CRPIX3' & sxdelpar,hdr,'CDELT3'
    sxdelpar,hdr,'CRVAL4' & sxdelpar,hdr,'CRPIX4' & sxdelpar,hdr,'CDELT4'
    sxaddpar,hdr,'NAXIS',2,'number of data axes'
    writefits,filepath+'5461_'+fileapp+'_MOM0_'+outapp[k]+'.fits', $
      grdat*(-cdelt3),hdr
    sxaddpar,hdr,'BITPIX',-64,'number of bits per data pixel'
    writefits,filepath+'5461_'+fileapp+'_MOM1_'+outapp[k]+'.fits',velocity,hdr
    writefits,filepath+'5461_'+fileapp+'_MOM2_'+outapp[k]+'.fits',dispersion,hdr
    
  endfor
  
END