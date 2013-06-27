pro write_arc_photometry

    jj = read_maskpops()
    maskpops_to_maggies, jj, maggies, ivarmaggies

    ff = clash_filterlist(weff=weff)
    nff = n_elements(ff)

    out = {weff: 0.0, flux1: 0.0, ivar1: 0.0}
    out = replicate(out,nff)
    out.weff = weff
    out.flux1 = maggies[*,0]    
    out.ivar1 = ivarmaggies[*,0]    

    outpath = getenv('IM_GITREPOS')+'/siena-astrophysics/summer2013/clash/'
    wsex, out, outfile=outpath+'arcphot.txt'

return
end
    
