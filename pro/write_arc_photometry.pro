pro write_arc_photometry

    jj = read_maskpops()
    maskpops_to_maggies, jj, maggies, ivarmaggies

    ff = clash_filterlist(weff=weff)
    nff = n_elements(ff)

    out = {weff: 0.0, $
      flux1: 0.0, ivar1: 0.0,$
      flux2: 0.0, ivar2: 0.0,$
      flux3: 0.0, ivar3: 0.0}
    out = replicate(out,nff)
    out.weff = weff
    out.flux1 = maggies[*,0]    ; total
    out.ivar1 = ivarmaggies[*,0]    
    out.flux2 = abs(maggies[*,3]) ; redside
    out.ivar2 = ivarmaggies[*,3]    
    out.flux3 = maggies[*,5]    ; blueside
    out.ivar3 = ivarmaggies[*,5]    

    outpath = getenv('IM_GITREPOS')+'/siena-astrophysics/summer2013/clash/'
    wsex, out, outfile=outpath+'arcphot.txt'

return
end
    
