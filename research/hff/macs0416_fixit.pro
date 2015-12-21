pro macs0416_fixit
; convert the published magnitudes back to nJy    

    cc = rsex('z9_candidates.cat')
    cc = cc[where(strmatch(cc.name,'*0416*'))]    
    filters = hff_filterlist(short_filter=filt,/useirac,usehawki=1)  
    tags = strupcase(filt+'_flux')
    errtags = strupcase(filt+'_fluxerr')

    out = cc

    factor = 10D^(0.4D*31.4)

    nf = n_elements(filt)
    ngal = 6
    for jj = 0, ngal-1 do begin
       for ii = 0, nf-1 do begin
          mag = [cc[jj].(tag_indx(cc[0],tags[ii]))]
          magerr = [cc[jj].(tag_indx(cc[0],errtags[ii]))]
          maggies1 = mag2maggies(mag,magerr=magerr,ivarmaggies=ivar1)

          if ivar1 gt 0 then begin
             out[jj].(tag_indx(out[0],tags[ii])) = maggies1*factor
             out[jj].(tag_indx(out[0],errtags[ii])) = 1.0/sqrt(ivar1)*factor
          endif else begin
;         if magerr gt 0 and mag eq 0 then $
            out[jj].(tag_indx(out[0],errtags[ii])) = mag2maggies(magerr)*factor
         endelse
       endfor
    endfor

    wsex, out, outfile='z9_candidates.fix'

stop

return
end
