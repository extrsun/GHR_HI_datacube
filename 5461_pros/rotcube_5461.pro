pro rotcube_5461,velbegin=velbegin,velend=velend,rotang=rotang
;+
; Name
;      rotcube_5461
; Purpose
;      Rotate the data cube to let the spiral arm follows the direction of 
;      the X-axis     
; Parameters
;      velbegin - the beginning velocity of the folding; def: 224.779 km/s
;        velend - the ending velocity of the folding; def: 302.257 km/s
;        rotang - the rotation angle of the data cube; def: 51 deg from 
;                 north to east
;      
; Inherited from "accum_ek_5461.pro"
; Modified by sw, Jan 11, 2015
;-

  ;Procedure keyword setting
  ; if keyword_set(velbegin) eq 0 then velbegin=204000.D
  ; if keyword_set(velend)   eq 0 then velend=308000.D
  if keyword_set(velbegin) eq 0 then velbegin=80000.D
  if keyword_set(velend) eq 0 then velend=500000.D
  if keyword_set(rotang) eq 0 then rotang=-51

  ;System paramters
  genpath='../'
  cubename=genpath+'5461_RO_CUBE_enlarge.fits' ;Use the enlarged HI CUBE
  ;The CUBE file of the other weighting, 'NA', is alternative.
  
  ;Read data
  rawdata=mrdfits(cubename,0,hdr)
  crval1=sxpar(hdr,'crval1') & crval2=sxpar(hdr,'crval2')
  crpix1=sxpar(hdr,'crpix1') & crpix2=sxpar(hdr,'crpix2')
  cdelt1=sxpar(hdr,'cdelt1') & cdelt2=sxpar(hdr,'cdelt2')
  crval3=sxpar(hdr,'crval3') & crpix3=sxpar(hdr,'crpix3')
  cdelt3=sxpar(hdr,'cdelt3')  
  ;Pre-select data
  rawsz=size(rawdata) & velarr=(dindgen(rawsz(3))-crpix3)*cdelt3+crval3
  velselind=where(velarr gt velbegin and velarr lt velend,nvelselind)
  if nvelselind eq 0 then begin
    print,'The velocity range seems not right!'
    return
  endif
  data=fltarr(rawsz[1],rawsz[2],nvelselind)
  orihdr=hdr
  sxaddpar,hdr,'NAXIS',2 & sxdelpar,hdr,'NAXIS3'
  sxdelpar,hdr,'NAXIS4'
  for i=0,nvelselind-1 do begin
    tmpdata=rawdata[*,*,velselind[i]]
    hrot,tmpdata,hdr,seldata,outhdr,rotang,123,121,1,mis=0
    data[*,*,i]=seldata
  endfor
  ajdlist=['CRVAL1','CRVAL2','CRPIX1','CRPIX2','CDELT1','CDELT2', $
    'CROTA1','CROTA2'] & najdlist=n_elements(ajdlist)
  ; ajdlist=['CROTA1','CROTA2'] & najdlist=n_elements(ajdlist)
  for i=0,najdlist-1 do begin
    tmpval=sxpar(outhdr,ajdlist[i])
    sxaddpar,orihdr,ajdlist[i],tmpval
  endfor
  sxaddpar,orihdr,'NAXIS3',nvelselind
  sxaddpar,orihdr,'CRPIX3',crpix3-velselind[0]+1
  writefits,genpath+'5461_RO_CUBE_rot.fits',data,orihdr
  
  for i=0,2 do begin
    data=mrdfits(genpath+'5461_RO_MOM'+string(i,format='(i1)')+ $
      '_enlarge.fits',0,hdr)
    hrot,data,hdr,-1,-1,rotang,123,121,1,mis=0
    writefits,genpath+'5461_RO_MOM'+string(i,format='(i1)')+ $
      '_rot.fits',data,hdr
  endfor

END