pro accum_ek_5461,velbegin=velbegin,velend=velend,ymin=ymin,ymax=ymax, $
  subdir=subdir,binsize=binsize
;+
; Name
;      accum_ek_5461
; Purpose
;      Calculate the accumulative kinematic energy at different positions 
;      along the spiral arm around NGC5461.     
; Parameters
;      velbegin - the beginning velocity of the folding; def: 224.779 km/s
;        velend - the ending velocity of the folding; def: 302.257 km/s
;          ymin - the minimum of the y-coordinate, def: 80
;          ymax - the maximum of the y-coordinate, def:160
;        subdir - the name of the sub-directory
;       binsize - the bin size along the X-axis
;      
; Written by sw, Jan 04, 2015
;-

  ;Procedure keyword setting
  ; if keyword_set(velbegin) eq 0 then velbegin=80000.D
  ; if keyword_set(velend) eq 0 then velend=500000.D
  if keyword_set(velbegin) eq 0 then velbegin=170000.D;224000.D
  if keyword_set(velend)   eq 0 then velend=320000.D;303000.D
  if keyword_set(ymin) eq 0 then ymin=80
  if keyword_set(ymax) eq 0 then ymax=159
  if keyword_set(subdir) eq 0 then subdir='wide'
  if keyword_set(rotang) eq 0 then rotang=-51.
  if keyword_set(binsize) eq 0 then binsize=8
  ynum=ymax-ymin+1 & ycenter=123 & xcenter=124.5
  xnum=250/binsize & illnum=3
  labtext=['H1105','H1098'] & xnamelab=[[115,132],[138,144]]-1;-xcenter
  pixarea=binsize*ynum
  areastr='Binsize: '+strtrim(string(pixarea,format='(i4)'),2)+' pixels'
  
  ;Data threshold
  threshold=1.152e-3 ;2.4-sigma
  ; threshold=1.44e-3  ;3-sigma
  
  ;Net fluctuation kinetic energy ((cm/s)^2)
  if subdir eq 'wide' then flucdek=[0.5,0.15] else flucdek=[0.5,0.2]
  flucdek*=1e12
  
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
  cubename=genpath+'5461_RO_CUBE_contisel.fits' ;Use the enlarged HI CUBE
  ;The CUBE file of the other weighting, 'NA', is alternative.

  ;Plot parameters
  set_plot,'ps'
  device,/encapsulated,preview=2,bits=8,/color
  tvlct,[0,255,  0,  0,  0,255,255,  0,  0,127,255], $
        [0,  0,255,  0,255,  0,127,127,255,  0,  0], $
        [0,  0,  0,255,255,255,  0,255,127,255,127]
    ;Black,  R,  G,B,cyan,purp,orag,...
  charsize=1.5 & thick=5 & ancthick=3 & lchsize=1.2 & ekchsize=1.8
  
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
  velvalarr=dblarr(binsize,ynum,nvelselind)
  for k=0,binsize-1 do for l=0,ynum-1 do $
    velvalarr[k,l,*]=velarr[velselind]
  velvalarr*=1e2 ;Convert the velocity from m/s to cm/s (cgs unit)
  presdata=fltarr(rawsz[1],rawsz[2],nvelselind)
  sxaddpar,hdr,'NAXIS',2 & sxdelpar,hdr,'NAXIS3'
  sxdelpar,hdr,'NAXIS4'
  data=rawdata[*,ymin:ymax,velselind]
  ;Convert the data into HI mass (in gram)
  data*=data2mass ;in gram/pixel
  
  sublist=['enlarge','integr'] & nsublist=n_elements(sublist)
  xindl=indgen(xnum)*binsize & xindr=(indgen(xnum)+1)*binsize-1
  grint=dblarr(xnum) & grvel=dblarr(xnum)
  grdat=dblarr(xnum,ynum) & grvdat=dblarr(xnum,ynum)
  grekv1=dblarr(xnum,ynum) & grekv2=dblarr(xnum,ynum)
  cenvel=dblarr(binsize,ynum,nvelselind)
  ;Use two kinds of velocity map: original ('enlarge') or integrated ('integr')
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
    
    ;Plot the kinematic energy profile along y-axis at variant positions
    device,filename=figpath+'ekprof_'+sublist[i]+'_'+subdir+'_'+ $
      string(binsize,format='(i02)')+'_r51.eps'
    !p.multi=[0,3,2] & colind=indgen(illnum)+1 & thkind=intarr(illnum)+3
    lstind=intarr(illnum) & psmind=intarr(illnum)+6 & xrange=[ymin,ymax]
    yind=indgen(ynum)+ymin
    illran=lindgen(illnum*6)+xnum/2-3*illnum
    plymin1=min([grekv1[illran,*],grekv2[illran,*]],max=plymax1)
    plydel1=plymax1-plymin1
    yrange1=[plymin1-plydel1*1e-4,plymax1+plydel1*0.2]/divfac_ek
    plymin2=min(grmass[illran,*],max=plymax2) & plydel2=plymax2-plymin2
    yrange2=[(plymin2-plydel2*1e-4)>0,plymax2+plydel2*0.2]/divfac_mass
    for q=0,1 do begin
      for p=0,2 do begin
        if q eq 0 then begin
          xtickname=replicate(' ',7) & xtitle=' '
        endif else xtitle='!6Y-Coor (pix)'
        ytickname=replicate(' ',9) & ytitle=' '
        if p eq 0 then ytitle='!6E!Ik!N (10!U'+ $
            strtrim(string(round(alog10(divfac_ek))),2)+ $
              '!Nerg s!U-1!Npixel!U-1!N)'
        position=[0.08+p*0.28,0.53-q*0.45,0.08+(p+1)*0.28,0.98-q*0.45]
        illind=indgen(illnum)+3*illnum*q+illnum*p+xnum/2-3*illnum
        pt=string(fix(alog10((illind[illnum-1]+1)*binsize)+1),format='(i1)')
        pixbeg=string(illind*binsize+1,format='(i'+pt+')')
        pixend=string((illind+1)*binsize,format='(i'+pt+')')
        text='!6X: '+pixbeg+'-'+pixend
        ;Plot the frame
        if p eq 0 then begin
          if q eq 0 then $
            plot,xrange,yrange1,xrange=xrange,/nodata,yrange=yrange1, $
              xstyle=1,ystyle=9,xtitle=xtitle,ytitle=ytitle, $
              charsize=ekchsize,xtickname=xtickname,position=position $
          else $
            plot,xrange,yrange1,xrange=xrange,/nodata,yrange=yrange1, $
              xstyle=1,ystyle=9,xtitle=xtitle,ytitle=ytitle, $
              charsize=ekchsize,position=position
        endif else begin
          if q eq 0 then $
            plot,xrange,yrange1,xrange=xrange,/nodata,yrange=yrange1, $
              xstyle=1,ystyle=9,xtitle=xtitle,ytitle=ytitle, $
              position=position,charsize=ekchsize,xtickname=xtickname, $
              ytickname=ytickname $
          else $
            plot,xrange,yrange1,xrange=xrange,/nodata,yrange=yrange1, $
              xstyle=1,ystyle=9,xtitle=xtitle,ytitle=ytitle, $
              charsize=ekchsize,ytickname=ytickname,position=position
        endelse
        for m=0,illnum-1 do begin
          pldat=grekv1[illind[m],*]/divfac_ek & selind=where(pldat gt 0.0)
          oplot,yind[selind],pldat[selind],color=colind[m], $
            thick=thkind[m],linestyle=2
          pldat=grekv2[illind[m],*]/divfac_ek & selind=where(pldat gt 0.0)
          oplot,yind[selind],pldat[selind],color=colind[m], $
            thick=thkind[m],linestyle=1
        endfor
        if p eq 2 then $
          axis,/yaxis,yrange=yrange2,/save,ystyle=1,charsize=ekchsize, $
            ytitle='!6HI Mass (10!U'+strtrim(string(alog10(divfac_mass), $
              format='(i1)'),2)+'!N!8M!6!I!9n!6!N pixel!U-1!N)' $
        else $
          axis,/yaxis,yrange=yrange2,/save,ystyle=1,charsize=ekchsize, $
            ytitle=' ',ytickname=replicate(' ',7)
        for m=0,illnum-1 do begin
          pldat=grmass[illind[m],*]/divfac_mass & selind=where(pldat gt 0.0)
          oplot,yind[selind],pldat[selind],color=colind[m], $
            thick=thkind[m],linestyle=0
        endfor
        legend,text,color=colind,thick=thkind,linestyle=lstind, $
          charsize=charsize/2.5,bottom=bottom,right=right
      endfor
    endfor
    !p.multi=0
    device,/close
    
    ;Plot the accumulative kinematic energy profile along x-axis
    device,filename=figpath+'acekpf_'+sublist[i]+'_'+subdir+'_'+ $
      string(binsize,format='(i02)')+'_r51.eps'
    xind=(xindl+xindr)/2.;-xcenter
    xrange=[xind[0],xind[xnum-1]]
    plymin1=min([accumek1,accumek2],max=plymax1)
    plymin1/=divfac_ek & plymax1/=divfac_ek & plydel1=plymax1-plymin1
    yrange1=[plymin1-plydel1*0.1,plymax1+plydel1*0.1]
    plymin2=min(selmass,max=plymax2)
    plymin2/=divfac_mass & plymax2/=divfac_mass & plydel2=plymax2-plymin2
    yrange2=[(plymin2-plydel2*0.1)>0,plymax2+plydel2*0.1]
    plot,xrange,yrange1,/nodata,xrange=xrange1,yrange=yrange,xstyle=1, $
      charsize=charsize,xmargin=[8,7],xtitle='!6X-Coordinate (pixel)', $
      ystyle=9,ytitle='Kinetic Energy (10!U'+ $
        strtrim(string(round(alog10(divfac_ek))),2)+ $
        '!Nerg s!U-1!Npixel!U-1!N)';,thick=thick
    oplot,xind,accumek1/divfac_ek,thick=thick,color=1
    oplot,xind,accumek2/divfac_ek,thick=thick,color=2
    oplot,[0,0]+xcenter,yrange1,thick=ancthick,linestyle=2
    xyouts,xcenter,yrange1[0]+0.05*plydel1,'H1105 Center',charsize=lchsize
    ylab=(plymax1+plydel1*[0.01,-0.23])
    for m=0,1 do begin
      oplot,xnamelab[*,m],[ylab[m],ylab[m]],thick=ancthick
      for n=0,1 do $
        oplot,xnamelab[n,m]*[1,1],ylab[m]+[-1,1]*2e-2*plydel1,thick=ancthick
      xyouts,xnamelab[m,m],ylab[m],labtext[m],alignment=1-m,charsize=lchsize
    endfor
    ylab=(plymax1-plydel1*0.2) & xlab=35
    xyouts,xlab-25,ylab+plydel1*0.17,subdir,charsize=lchsize*2
    xyouts,xlab-25,ylab+plydel1*0.06,areastr,charsize=lchsize
    xyouts,xlab,ylab-plydel1*0.06,'Beam Size',charsize=lchsize, $
      alignment=0.5
    oplot,xlab+[-3.75,3.75],ylab*[1,1],thick=ancthick
    oplot,-3.75*[1,1]+xlab,ylab+[-1,1]*1e-2*plydel1,thick=ancthick
    oplot, 3.75*[1,1]+xlab,ylab+[-1,1]*1e-2*plydel1,thick=ancthick
    axis,/yaxis,yrange=yrange2,ystyle=1,ytitle='HI Mass (10!U'+ $
      strtrim(string(alog10(divfac_mass),format='(i2)'),2)+ $
      '!N!8M!6!I!9n!6!N pixel!U-1!N)',charsize=charsize,/save
    oplot,xind,selmass/divfac_mass,thick=thick,color=0
    legend,['HI Mass','E!Ik!N (MeanVel)','E!Ik!N (YvarVel)'],thick=[3,3,3], $
      color=[0,1,2],linestyle=[0,0,0],/right
    device,/close
    
    ;Plot the density of the accumulative kinematic energy profile along x-axis
    device,filename=figpath+'acdekpf_'+sublist[i]+'_'+subdir+'_'+ $
      string(binsize,format='(i02)')+'_r51.eps'
    xind=(xindl+xindr)/2.;-xcenter
    xrange=[xind[0],xind[xnum-1]]
    pldat1=accumek1/(selmass*!msolar)/divfac_ekd
    pldat2=accumek2/(selmass*!msolar)/divfac_ekd
    plymin1=min([pldat1,pldat2],max=plymax1)
    plymax1=plymax1+0.15*(plymax1-plymin1) & plydel1=plymax1-plymin1
    yrange1=[plymin1-plydel1*0.1,plymax1+plydel1*0.1]
    pldat3=selmass/divfac_mass
    plymin2=min(pldat3,max=plymax2) & plydel2=plymax2-plymin2
    yrange2=[(plymin2-plydel2*0.1)>0,plymax2+plydel2*0.1]
    plot,xrange,yrange1,/nodata,xrange=xrange1,yrange=yrange,xstyle=1, $
      charsize=charsize,xmargin=[8,7],xtitle='!6X-Coordinate (pixel)', $
      ystyle=9,ytitle='Kinetic Energy Density (10!U'+ $
        strtrim(string(round(alog10(divfac_ekd))),2)+ $
        '!Ncm!U2!Ns!U-2!N)';,thick=thick
    oplot,xind,pldat1,thick=thick,color=1
    oplot,xind,pldat2,thick=thick,color=2
    oplot,[0,0]+xcenter,yrange1,thick=ancthick,linestyle=2
    xyouts,xcenter,yrange1[0]+0.05*plydel1,'H1105 Center',charsize=lchsize
    ylab=(plymax1+plydel1*[0.01,-0.23])
    for m=0,1 do begin
      oplot,xnamelab[*,m],[ylab[m],ylab[m]],thick=ancthick
      for n=0,1 do $
        oplot,xnamelab[n,m]*[1,1],ylab[m]+[-1,1]*2e-2*plydel1,thick=ancthick
      xyouts,xnamelab[m,m],ylab[m],labtext[m],alignment=1-m,charsize=lchsize
    endfor
    ylab=(plymax1-plydel1*0.2) & xlab=35
    xyouts,xlab-25,ylab+plydel1*0.17,subdir,charsize=lchsize*2
    xyouts,xlab-25,ylab+plydel1*0.06,areastr,charsize=lchsize
    xyouts,xlab,ylab-plydel1*0.06,'Beam Size',charsize=lchsize, $
      alignment=0.5
    oplot,xlab+[-3.75,3.75],ylab*[1,1],thick=ancthick
    oplot,-3.75*[1,1]+xlab,ylab+[-1,1]*1e-2*plydel1,thick=ancthick
    oplot, 3.75*[1,1]+xlab,ylab+[-1,1]*1e-2*plydel1,thick=ancthick
    axis,/yaxis,yrange=yrange2,ystyle=1,ytitle='HI Mass (10!U'+ $
      strtrim(string(alog10(divfac_mass),format='(i2)'),2)+ $
      '!N!8M!6!I!9n!6!N pixel!U-1!N)',charsize=charsize,/save
    oplot,xind,pldat3,thick=thick,color=0
    legend,['HI Mass','MeanVel','YvarVel'], $
      thick=thick+[0,0,0],color=[0,1,2],linestyle=[0,0,0],/right
    device,/close
    
    ;Plot the net (fluctuation subtracted) kinematic energy profile
    ;along the X-axis
    ; pldata3=(pldat2*divfac_ekd-flucdek[i])*selmass*!msolar/divfac_ek
    pldata1=(accumek2-selmass*!msolar*flucdek[i])/divfac_ek
    pldata2=selmass/divfac_mass
    device,filename=figpath+'acnetekpf_'+sublist[i]+'_'+subdir+'_'+ $
      string(binsize,format='(i02)')+'.eps'
    xind=(xindl+xindr)/2.;-xcenter
    xrange=[xind[0],xind[xnum-1]]
    plymin1=min(pldata1,max=plymax1) & plydel1=plymax1-plymin1
    yrange1=[plymin1-plydel1*0.1,plymax1+plydel1*0.1]
    plymin2=min(pldata2,max=plymax2) & plydel2=plymax2-plymin2
    yrange2=[(plymin2-plydel2*0.1)>0,plymax2+plydel2*0.1]
    plot,xrange,yrange1,/nodata,xrange=xrange1,yrange=yrange,xstyle=1, $
      charsize=charsize,xmargin=[8,7],xtitle='!6X-Coordinate (pixel)', $
      ystyle=9,thick=thick,ytitle='Net Kinetic Energy (10!U'+ $
        strtrim(string(round(alog10(divfac_ek))),2)+ $
          '!Nerg s!U-1!Npixel!U-1!N)'
    oplot,xind,pldata1,thick=thick,color=1
    ; oplot,xind,pldata3,thick=3,color=2
    oplot,[0,0]+xcenter,yrange1,thick=ancthick,linestyle=2
    xyouts,xcenter,yrange1[0]+0.05*plydel1,'H1105 Center',charsize=lchsize
    ylab=(plymax1+plydel1*[0.01,-0.23])
    for m=0,1 do begin
      oplot,xnamelab[*,m],[ylab[m],ylab[m]],thick=ancthick
      for n=0,1 do $
        oplot,xnamelab[n,m]*[1,1],ylab[m]+[-1,1]*2e-2*plydel1,thick=ancthick
      xyouts,xnamelab[m,m],ylab[m],labtext[m],alignment=1-m,charsize=lchsize
    endfor
    ylab=(plymax1-plydel1*0.2) & xlab=35
    xyouts,xlab-25,ylab+plydel1*0.17,subdir,charsize=lchsize*2
    xyouts,xlab-25,ylab+plydel1*0.06,areastr,charsize=lchsize
    xyouts,xlab,ylab-plydel1*0.06,'Beam Size',charsize=lchsize, $
      alignment=0.5
    oplot,xlab+[-3.75,3.75],ylab*[1,1],thick=ancthick
    oplot,-3.75*[1,1]+xlab,ylab+[-1,1]*1e-2*plydel1,thick=ancthick
    oplot, 3.75*[1,1]+xlab,ylab+[-1,1]*1e-2*plydel1,thick=ancthick
    axis,/yaxis,yrange=yrange2,ystyle=1,ytitle='HI Mass (10!U'+ $
      strtrim(string(alog10(divfac_mass),format='(i2)'),2)+ $
      '!N!8M!6!I!9n!6!N pixel!U-1!N)',charsize=charsize,/save
    oplot,xind,pldata2,thick=thick,color=0
    legend,['HI Mass','Net E!Ik!N (YvarVel)'],thick=thick+[0,0], $
      color=[0,1],linestyle=[0,0],/right
    device,/close
    
    ;Output the total kinetic energy
    print,'Use '+sublist[i]+' data in '+subdir+' range under x-bin of '+ $
      string(binsize,format='(i2)')+' pixels'
    for m=0,1 do begin
      ; print,xindl
      ; print,xindr
      ; print,xnamelab[*,m]
      selind=where(xindl+binsize ge xnamelab[0,m] and $
        xindr-binsize lt xnamelab[1,m],nselind)
      netval=total(pldata1[selind])*divfac_ek*pixarea
      rawval=total(accumek2[selind])*pixarea
      subval=total(selmass[selind]*!msolar*flucdek[i])*pixarea
      if m eq 0 then lun=lun1 else lun=lun2
      printf,lun,rawval,subval,netval,nselind*binsize,format='(3e10.3,i3)'
      print,'For '+labtext[m]+', raw:'+string(rawval,format='(e10.3)')+ $
        ', sub: '+string(subval,format='(e10.3)')+', net: '+ $
        string(netval,format='(e10.3)'),nselind*binsize
    endfor
    free_lun,lun1 & free_lun,lun2
    
  endfor

END