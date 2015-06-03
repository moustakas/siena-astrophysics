pro show_galaxy

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Example showing how the file format of Gadget can    ;;
;; be read-in in IDL                                    ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

frun='./output/'

window,xsize=1000,ysize=333
!P.multi=[0,3,1]

xlen=300

 video = IDLffVideoWrite('testmerger.mp4', Format='mp4')
 framerate = 5
 framedims = [1804,600]
 stream = video.AddVideoStream(framedims[0], framedims[1], framerate)

;for num=0,10 do begin
for num=0,801 do begin

    exts='000'
    exts=exts+strcompress(string(num),/remove_all)
    exts=strmid(exts,strlen(exts)-3,3)

    fname=frun+"/snapshot_"+exts
    fname=strcompress(fname,/remove_all)

    
    npart=lonarr(6)	
    massarr=dblarr(6)
    time=0.0D
    redshift=0.0D
    flag_sfr=0L
    flag_feedback=0L
    npartTotal=lonarr(6)	
    bytesleft=256-6*4 - 6*8 - 8 - 8 - 2*4-6*4
    la=intarr(bytesleft/2)


    print,fname
    openr,1,fname,/f77_unformatted
    readu,1,npart,massarr,time,redshift,flag_sfr,flag_feedback,npartTotal,la

    print,"Time= ", time

    N=total(npart)
    pos=fltarr(3,N)
    vel=fltarr(3,N)
    id=lonarr(N)
 
    ind=where((npart gt 0) and (massarr eq 0)) 
    if ind(0) ne -1 then begin
	Nwithmass= total(npart(ind))
        mass=fltarr(Nwithmass)
    endif else begin	
        Nwithmass= 0
    endelse

    readu,1,pos
    readu,1,vel
    readu,1,id
    if Nwithmass gt 0 then begin
      readu,1,mass
    endif

    NGas=npart(0)
    NHalo=npart(1)
    NDisk=npart(2)
    NBulge=npart(3)
    NStars=npart(4)

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
        xgas=fltarr(Ngas) &  ygas=fltarr(Ngas)  & zgas=fltarr(Ngas) & mgas=fltarr(Ngas)
        xgas(*)=pos(0,0:Ngas-1)
        ygas(*)=pos(1,0:Ngas-1)
        zgas(*)=pos(2,0:Ngas-1)
        if massarr(0) eq 0 then begin
            mgas(*)=mass(0:Ngas-1)	
        endif else begin
            mgas(*)= massarr(0)
	endelse
    endif

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


    if Nstars gt 0 then begin
        xstars=fltarr(Nstars) &  ystars=fltarr(Nstars)  & zstars=fltarr(Nstars)  & mstars=fltarr(Nstars)
        xstars(*)=pos(0,Nhalo+Ngas+Ndisk+Nbulge:Nhalo+Ndisk+Ngas+NBulge+Nstars-1)
        ystars(*)=pos(1,Nhalo+Ngas+Ndisk+Nbulge:Nhalo+Ndisk+Ngas+Nbulge+Nstars-1)
        zstars(*)=pos(2,Nhalo+Ngas+Ndisk+Nbulge:Nhalo+Ndisk+Ngas+NBulge+Nstars-1)
        if massarr(4) eq 0 then begin
	    skip=0L
            for t=0,2 do begin	
                if (npart(t) gt 0) and (massarr(t) eq 0) then begin
			skip=skip + npart(t)
                endif
            endfor
            mstars(*)=mass(0+skip:Nstars-1+skip)	
        endif else begin
            mstars(*)= massarr(4)
	endelse
    endif



    cgPS_Open, 'movie.ps', /Quiet
   

    cgplot, xdisk, ydisk, psym=3, xrange= [-xlen,xlen], yrange= [-xlen,xlen], xstyle=1,ystyle=1
    cgplot, xdisk, zdisk, psym=3, xrange= [-xlen,xlen], yrange= [-xlen,xlen], xstyle=1,ystyle=1
    cgplot, ydisk, zdisk, psym=3, xrange= [-xlen,xlen], yrange= [-xlen,xlen], xstyle=1,ystyle=1

;   wait, 0.5	

      cgPS_Close, /PNG, Width=600, /Delete_PS ; Convert to PNG file.

      image = Read_PNG('movie.png')
      File_Delete, 'movie.png'
    
      void = video -> Put(stream, image)
      
   endfor

video -> Cleanup


return
end
