pro multiangle,ymin=ymin,ymax=ymax,subdir=subdir,binsize=binsize
;+
; Name
;      multiangle
; Purpose
;      Rotate the HI map by variant angles to check which one is optimal.
; Parameters
;      ymin - the minimum of the y-coordinate, def: 80
;      ymax - the maximum of the y-coordinate, def:160
;    subdir - the name of the sub-directory
;   binsize - the bin size along the X-axis
;
; Written by sw, Jan 07, 2015
;-

  ;Keyword setting
  if keyword_set(ymin) eq 0 then ymin=80
  if keyword_set(ymax) eq 0 then ymax=159
  if keyword_set(subdir) eq 0 then subdir='wide/' else subdir+='/'
  if keyword_set(binsize) eq 0 then binsize=8
  ynum=ymax-ymin+1 & ycenter=123 & xcenter=124.5
  xnum=250/binsize & illnum=4

  ;System paramters
  genpath='../rothi/'
  ;Plot parameters
  set_plot,'ps'
  device,/encapsulated,preview=2,bits=8,/color
  usersym,[-1,-1,1,1,-1],[-1,1,1,-1,-1],/fill
  tvlct,[0,255,  0,  0,  0,255,255,  0,  0,127,255], $
        [0,  0,255,  0,255,  0,127,127,255,  0,  0], $
        [0,  0,  0,255,255,255,  0,255,127,255,127]
    ;Black,  R,  G,B,cyan,purp,orag,...

  sublist=['enlarge','integr'] & nsublist=n_elements(sublist)
  ; anglist=-dindgen(7)-48 & nanglist=n_elements(anglist)
  anglist=-dindgen(1)-51 & nanglist=n_elements(anglist)
  xindl=indgen(xnum)*binsize & xindr=(indgen(xnum)+1)*binsize-1
  grint=dblarr(xnum) & grvel=dblarr(xnum)
  grdat=dblarr(xnum,ynum) & grvdat=dblarr(xnum,ynum)
  for i=0,nsublist-1 do begin
    filepath=genpath+sublist[i]+'/'
    subpath=filepath+subdir & figpath=subpath+'figures/'
    intdata=mrdfits(filepath+'5461_MOM0_'+sublist[i]+'.fits',0,hdr1)
    veldata=mrdfits(filepath+'5461_MOM1_'+sublist[i]+'.fits',0,hdr2)
    for j=0,nanglist-1 do begin
      angstr=string(-anglist[j],format='(i2)')
      hrot,intdata,hdr1,intensity,inthdr,anglist[j],123,121,1,mis=0;,/pivot
      hrot,veldata,hdr2,velocity,inthdr,anglist[j],123,121,1,mis=0;,/pivot
      writefits,subpath+'5461_MOM0_'+sublist[i]+'_r'+angstr+'.fits', $
        intensity,inthdr
      writefits,subpath+'5461_MOM1_'+sublist[i]+'_r'+angstr+'.fits', $
        velocity,inthdr
      ;Calculate the intensity-weighted y-coordinate at variant X-positions
      fdat=intensity[*,ymin:ymax]
      for k=0,xnum-1 do $
        grdat[k,*]=total(fdat[xindl[k]:xindr[k],*],1)
      grdat=grdat>0
      selint=total(grdat,2)
      yind=(dindgen(ynum,xnum) mod ynum)
      yind=transpose(yind)+ymin
      selpos=total(grdat*yind,2)/selint
      sposind=dblarr(xnum,ynum)
      for k=0,xnum-1 do sposind[k,*]=selpos[k]
      selposerr=sqrt(total(grdat*(yind-sposind)^2,2)/selint)
      ;Plot the intensity profile along y-axis at positions close to
      ;H1105 center
      ; if i eq 1 then print,!p.multi
      device,filename=figpath+'intprof_r'+angstr+'.eps'
      yind=indgen(ynum)+ymin & xrange=[ymin,ymax]
      maxval=max(grdat[13:17,*]) & yrange=[0,maxval*1.1]
      plot,xrange,yrange,/nodata,xrange=xrange,yrange=yrange,xstyle=1, $
        ystyle=1,xtitle='!6Y-Coordinate (pixel)',ytitle='Intensity', $
        charsize=1.5
      for k=13,17 do begin
        if k eq 15 then thick=5 else thick=2
        oplot,yind,grdat[k,*],thick=thick,color=k-12;,linestyle=linestyle
      endfor
      oplot,[ycenter,ycenter],yrange,linestyle=2
      xyouts,ycenter,yrange[0],'!6Y-center',charsize=1.5
      device,/close
      ;Plot the intensity weighted y-coordinates and corresponding error
      device,filename=figpath+'grpos_r'+angstr+'.eps'
      xind=(xindl+xindr)/2. & xrange=[xind[0],xind[xnum-1]]
      ; yrange=[(ycenter+ymin)/2.,(ycenter+ymax)/2.]
      yrange=[ymin,ymax]
      plot,xind,selpos,psym=6,xrange=xrange,yrange=yrange,xstyle=2, $
        ystyle=1,charsize=1.5,xtitle='!6X-Coordinate (pixel)', $
        ytitle='Weighted Y-coor (pixel)',thick=3
      errplot,xind,selpos-selposerr,selpos+selposerr,thick=3
      oplot,xrange,[ycenter,ycenter],thick=2,linestyle=2
      oplot,[xcenter,xcenter],yrange,thick=2,linestyle=2
      xyouts,xrange[0],ycenter,'!6Y-center',charsize=1.5;,alignment=1.0
      xyouts,xcenter,yrange[0],'!6X-center',charsize=1.5
      res=linfit(xind,selpos,measure_err=selposerr)
      xyouts,xcenter,ymax*0.9+ymin*0.1,'angle: '+ $
        strtrim(string(atan(res[1])*180/!dpi,format='(f7.3)'),2), $
        charsize=1.8
      device,/close

      ;Calculate the intensity-weighted velocity at variant X-positions
      vdat=velocity[*,ymin:ymax]
      for k=0,xnum-1 do begin
        selvdat=vdat[xindl[k]:xindr[k],*]
        selfdat=fdat[xindl[k]:xindr[k],*]>0.0
        for l=0,ynum-1 do begin
          sind=where(finite(selvdat[*,l]))
          selgrfdat=total(selfdat[sind,l])
          if selgrfdat gt 0 then begin
            grvdat[k,l]=total(selvdat[sind,l]*selfdat[sind,l])/selgrfdat
            ; grverr[k,l]=sqrt(total(selfdat[sind,l]* $
            ;   (selvdat[sind,l]-grvdat[k,l])^2)/selgrfdat)
          endif else grvdat[k,l]=0.0
          if total(selfdat[sind,l]) eq 0.0 then $
            print,sublist[i],angstr,xindl[k],xindr[k],l+ymin;,grvdat[k,l]
        endfor
      endfor
      grvdat/=1e3
      selvel=total(grvdat*grdat,2)/selint
      svelind=dblarr(xnum,ynum)
      for k=0,xnum-1 do svelind[k,*]=selvel[k]
      selvelerr=sqrt(total(grdat*(grvdat-svelind)^2,2)/selint)
      ;Plot the intensity-weighted velocity profile along the X-axis
      device,filename=figpath+'grvel_r'+angstr+'_'+ $
        string(binsize,format='(i02)')+'.eps'
      xind=(xindl+xindr)/2. & xrange=[xind[0],xind[xnum-1]]
      plymin=min(selvel-selvelerr) & plymax=max(selvel+selvelerr)
      plydel=plymax-plymin
      yrange=[plymin-plydel*0.1,plymax+plydel*0.1]
      plot,xind,selvel,psym=6,xrange=xrange,yrange=yrange,xstyle=2, $
        ystyle=1,charsize=1.5,xtitle='!6X-Coordinate (pixel)', $
        ytitle='Weighted Velocity (km/s)',thick=3
      errplot,xind,selvel-selvelerr,selvel+selvelerr,thick=3
      oplot,[xcenter,xcenter],yrange,thick=2,linestyle=2
      xyouts,xcenter,yrange[0],'!6X-center',charsize=1.5
      device,/close
      ;Plot the velocity profile along y-axis at variant positions
      device,filename=figpath+'velprof_r'+angstr+'_'+ $
        string(binsize,format='(i02)')+'.eps'
      !p.multi=[0,3,2] & colind=indgen(illnum)+1
      thkind=intarr(illnum)+3 & lstind=intarr(illnum)
      symind=intarr(illnum)-8 & xrange=[ymin,ymax]
      charsize=1.8
      selind=where(grvdat gt 0.0) 
      plymin=min(grvdat[selind],max=plymax) & plydel=plymax-plymin
      yrange=[plymin-plydel*1e-4,plymax+plydel*1e-4]
      for q=0,1 do begin
        for p=0,2 do begin
          if q eq 0 then begin
            xtickname=replicate(' ',7) & xtitle=' '
          endif else xtitle='!6Y-Coor (pix)'
          ytickname=replicate(' ',9) & ytitle=' '
          if p eq 0 then ytitle='!6Velocity (km/s)'
          position=[0.08+p*0.30,0.53-q*0.45,0.38+p*0.30,0.98-q*0.45]
          illind=indgen(illnum)+3*illnum*q+illnum*p+xnum/2-3*illnum
          pt=string(fix(alog10((illind[illnum-1]+1)*binsize)+1), $
            format='(i1)')
          pixbeg=string(illind*binsize+1,format='(i'+pt+')')
          pixend=string((illind+1)*binsize,format='(i'+pt+')')
          text='!6X-Pix:'+pixbeg+'-'+pixend
          ;Plot the frame
          if p eq 0 then begin
            if q eq 0 then $
              plot,xrange,yrange,xrange=xrange,/nodata,yrange=yrange, $
                xstyle=1,ystyle=1,xtitle=xtitle,ytitle=ytitle, $
                charsize=charsize,xtickname=xtickname,position=position $
            else $
              plot,xrange,yrange,xrange=xrange,/nodata,yrange=yrange, $
                xstyle=1,ystyle=1,xtitle=xtitle,ytitle=ytitle, $
                charsize=charsize,position=position
          endif else begin
            if q eq 0 then $
              plot,xrange,yrange,xrange=xrange,/nodata,yrange=yrange, $
                xstyle=1,ystyle=1,xtitle=xtitle,ytitle=ytitle, $
                position=position,charsize=charsize,xtickname=xtickname, $
                ytickname=ytickname $
            else $
              plot,xrange,yrange,xrange=xrange,/nodata,yrange=yrange, $
                xstyle=1,ystyle=1,xtitle=xtitle,ytitle=ytitle, $
                charsize=charsize,ytickname=ytickname,position=position
          endelse
          for m=0,illnum-1 do begin
            pldat=grvdat[illind[m],*] & selind=where(pldat gt 0.0)
            oplot,yind[selind],pldat[selind],color=colind[m],thick=thkind[m]
            oplot,selpos[illind[m]]*[1,1],selvel[illind[m]]*[1,1], $
              psym=-symind[m],color=colind[m]
            ; if m eq 0 and p eq 0 then print,s
          endfor
          if q eq 0 then begin
            bottom=1 & right=1
          endif else begin
            bottom=0 & right=0
          endelse
          legend,text,color=colind,thick=thkind,linestyle=lstind, $
            charsize=charsize/2.5,bottom=bottom,right=right
        endfor
      endfor
      !p.multi=0
      device,/close
    endfor
  endfor

  set_plot,'x'

END
