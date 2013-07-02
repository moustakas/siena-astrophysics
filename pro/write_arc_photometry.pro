pro write_arc_photometry

    outpath = getenv('IM_GITREPOS')+'/siena-astrophysics/summer2013/clash/'

; restore the best-fit models
    prefix = 'maskpops'
    isedfit_dir = maskpops_path(/isedfit)
    montegrids_dir = maskpops_path(/montegrids)
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'
    sfhgrid_paramfile = isedfit_dir+prefix+'_sfhgrid.par'
    supergrid_paramfile = isedfit_dir+prefix+'_supergrid.par'

    model = isedfit_restore(isedfit_paramfile,isedfit,params=params,$
      isedfit_dir=isedfit_dir,supergrid_paramfile=supergrid_paramfile,$
      thissupergrid=thissupergrid,montegrids_dir=montegrids_dir,index=index,$
      outprefix=outprefix,silent=silent)

    npix = n_elements(model[0].wave)
    out = replicate({wave1: 0.0, flux1: 0.0, wave2: 0.0, flux2: 0.0, wave3: 0.0, flux3: 0.0},npix)
    out.wave1 = model[0].wave
    out.flux1 = model[0].flux
    out.wave2 = model[3].wave
    out.flux2 = model[3].flux
    out.wave3 = model[5].wave
    out.flux3 = model[5].flux

    wsex, out, outfile=outpath+'arcphot_bestfit.txt'
    
stop

; photometry    
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

    wsex, out, outfile=outpath+'arcphot.txt'

return
end
    
