;+
; NAME:
;   GADGET2_MOVIE
;
; PURPOSE:
;   Build a movie from a sequence of Gadget2 snapshots.
;
; INPUTS: 
;
; OPTIONAL INPUTS: 
;
; KEYWORD PARAMETERS: `
;
; OUTPUTS: 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2014 Feb 07, Siena
;-

pro gadget2_movie, moviename, outdir=outdir

    if n_elements(outdir) eq 0 then outdir = 'output/'

    if n_elements(moviename) eq 0 then begin
       print, 'Please specify the output name of the movie (.mp4 file)'
       return
    endif

    set_plot, 'Z'
    cgdisplay, xsize=1000, ysize=333
    !p.multi=[0,3,1]
    
; movie stuff
    video = IDLffVideoWrite(moviename, Format='mp4')
    print, 'Creating Gadget2 movie '+moviename
    framerate = 5
    framedims = [600,200] ; [1804,600]
    stream = video.AddVideoStream(framedims[0],framedims[1],framerate)
    
    xlen = 300 ; physical diameter of the simulation [kpc?]

; gather the snapshot files
    snapfiles = file_search(outdir+'snapshot_*',count=nsnap)
    
    for ii = 0, 5 do begin
;   for ii = 0, nsnap-1 do begin
       print, format='("Processing Snapshot ",I0,"/",I0, A10,$)', ii+1, nsnap, string(13b)
       
; set up some variables       
       npart = lonarr(6) ; number of particles per type
       massarr = dblarr(6)
       time = 0D
       redshift = 0D
       flag_sfr = 0L
       flag_feedback = 0L
       npartTotal = lonarr(6)   
       bytesleft = 256-6*4 - 6*8 - 8 - 8 - 2*4-6*4
       la = intarr(bytesleft/2)

; read the snapshot       
       fname = snapfiles[ii]
       openr, 1, fname, /f77_unformatted
       readu, 1, npart, massarr, time, redshift, flag_sfr, flag_feedback, npartTotal, la
;      print, 'Time= ', time
       
       N = total(npart)         ; total number of particles
       pos = fltarr(3,N)        ; position
       vel = fltarr(3,N)        ; velocity
       id = lonarr(N)           ; type of particle
       
       ind = where((npart gt 0) and (massarr eq 0)) 
       if ind(0) ne -1 then begin
          Nwithmass = total(npart[ind])
          mass = fltarr(Nwithmass)
       endif else begin 
          Nwithmass= 0
       endelse
       
       readu, 1, pos
       readu, 1, vel
       readu, 1, id
       if Nwithmass gt 0 then begin
          readu, 1, mass
       endif
       
       NGas = npart[0]
       NHalo = npart[1]
       NDisk = npart[2]
       NBulge = npart[3]
       NStars = npart[4]

; read the gas 
       if Ngas gt 0 then begin
          u=fltarr(Ngas)
          readu,1,u
          
          rho=fltarr(Ngas)
          readu,1,rho
          
          if flag_sfr gt 0 then begin
             sfr=fltarr(Ngas)
             mfs=fltarr(Ngas)
             readu,1,sfr
             readu,1,mfs
          endif
       endif
       close,1

       if Ngas gt 0 then begin
          xgas = fltarr(Ngas) &  ygas=fltarr(Ngas)  & zgas=fltarr(Ngas) & mgas=fltarr(Ngas)
          xgas(*)=pos(0,0:Ngas-1)
          ygas(*)=pos(1,0:Ngas-1)
          zgas(*)=pos(2,0:Ngas-1)
          if massarr(0) eq 0 then begin
             mgas(*)=mass(0:Ngas-1)     
          endif else begin
             mgas(*)= massarr(0)
          endelse
       endif
       
; read the halo 
       if Nhalo gt 0 then begin
          xhalo=fltarr(NHalo) &  yhalo=fltarr(Nhalo) & zhalo=fltarr(Nhalo) & mhalo=fltarr(Nhalo)
          xhalo(*)=pos(0,0+Ngas:Nhalo+Ngas-1)
          yhalo(*)=pos(1,0+Ngas:Nhalo+Ngas-1)
          zhalo(*)=pos(2,0+Ngas:Nhalo+Ngas-1)
          if massarr(1) eq 0 then begin
             skip=0L
             for t=0,0 do begin 
                if (npart(t) gt 0) and (massarr(t) eq 0) then begin
                   skip=skip + npart(t)
                endif
             endfor
             mhalo(*)=mass(0+skip:Nhalo-1+skip) 
          endif else begin
             mhalo(*)= massarr(1)
          endelse
       endif
       
; read the disk 
       if Ndisk gt 0 then begin
          xdisk=fltarr(NDisk) &  ydisk=fltarr(NDisk) &  zdisk=fltarr(NDisk) & mdisk=fltarr(NDisk)
          xdisk(*)=pos(0,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)
          ydisk(*)=pos(1,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)
          zdisk(*)=pos(2,Nhalo+Ngas:Nhalo+Ndisk+Ngas-1)
          if massarr(2) eq 0 then begin
             skip=0L
             for t=0,1 do begin 
                if (npart(t) gt 0) and (massarr(t) eq 0) then begin
                   skip=skip + npart(t)
                endif
             endfor
             mdisk(*)=mass(0+skip:Ndisk-1+skip) 
          endif else begin
             mdisk(*)= massarr(2)
          endelse
       endif
       
; read the bulge
       if Nbulge gt 0 then begin
          xbulge=fltarr(NBulge) &  ybulge=fltarr(NBulge) &  zbulge=fltarr(NBulge) & mbulge=fltarr(NBulge)
          xbulge(*)=pos(0,Nhalo+Ndisk+Ngas:Nhalo+Ndisk+Ngas+Nbulge-1)
          ybulge(*)=pos(1,Nhalo+Ndisk+Ngas:Nhalo+Ndisk+Ngas+Nbulge-1)
          zbulge(*)=pos(2,Nhalo+Ndisk+Ngas:Nhalo+Ndisk+Ngas+Nbulge-1)
          if massarr(3) eq 0 then begin
             skip=0L
             for t=0,2 do begin 
                if (npart(t) gt 0) and (massarr(t) eq 0) then begin
                   skip=skip + npart(t)
                endif
             endfor
             mbulge(*)=mass(0+skip:Nbulge-1+skip)       
          endif else begin
             mbulge(*)= massarr(3)
        endelse
       endif
       
; read the stars
       if Nstars gt 0 then begin
          xstars = fltarr(Nstars)
          ystars = fltarr(Nstars)
          zstars = fltarr(Nstars)
          mstars = fltarr(Nstars)
          
          xstars = pos[0,Nhalo+Ngas+Ndisk+Nbulge:Nhalo+Ndisk+Ngas+NBulge+Nstars-1]
          ystars = pos[1,Nhalo+Ngas+Ndisk+Nbulge:Nhalo+Ndisk+Ngas+Nbulge+Nstars-1]
          zstars = pos[2,Nhalo+Ngas+Ndisk+Nbulge:Nhalo+Ndisk+Ngas+NBulge+Nstars-1]
          if massarr[4] eq 0 then begin
             skip = 0L
             for t=0,2 do begin 
                if (npart(t) gt 0) and (massarr(t) eq 0) then begin
                   skip = skip + npart(t)
                endif
             endfor
             mstars = mass[0+skip:Nstars-1+skip]
          endif else begin
             mstars = massarr[4]
          endelse
       endif

; plot the simulation points
       cgps_open, 'movie.ps', /quiet
       cgplot, xdisk, ydisk, psym=3, xrange= [-xlen,xlen], yrange= [-xlen,xlen], $
         xstyle=1, ystyle=1, xthick=4, ythick=4
       cgplot, xdisk, zdisk, psym=3, xrange= [-xlen,xlen], yrange= [-xlen,xlen], $
         xstyle=1, ystyle=1, xthick=4, ythick=4
       cgplot, ydisk, zdisk, psym=3, xrange= [-xlen,xlen], yrange= [-xlen,xlen], $
         xstyle=1, ystyle=1, xthick=4, ythick=4
       cgps_close, /png, /delete_ps, /nomessage, width=600
      
       image = read_png('movie.png')
       File_Delete, 'movie.png'

; dump into the video       
       void = video -> Put(stream, image)
    endfor                      ; close the snapshot loop
    print

; close the video    
    video -> Cleanup
    set_plot, 'X'
        
return
end
