#!/usr/bin/env python
import sys
import os 
import asciidata
import numpy as np
import pyregion
import pyfits
run=os.system

if __name__ == "__main__":
    roidata = sys.argv[1]
else:
    roidata='roiphot_startup.dat'
if not (os.path.isfile(roidata)):
    print "missing startup file"
    sys.exit()

phothead = \
'# 1 filter\n'+\
'# 2 ABmag_zpt\n'+\
'# 3 threesigmaAB\n'+\
'# 4 threesigmafnu\n'+\
'# 5 area\n'+\
'# 6 fnu\n'+\
'# 7 fnuerr\n'+\
'# 8 fnuerr_shotnoise\n'+\
'# 9 fnuerr_skyrms\n'+\
'# 10 ABmag\n'+\
'# 11 ABmagerr\n'
    
topdir='/moustakas-archive/clash-archive/'
#topdir='/Volumes/Archive/CLASH/archive.stsci.edu/pub/clash/outgoing/'
if not (os.path.isdir(topdir)):
    print "CLASH data archive may not be mounted"
    sys.exit()
subdir='/HST/images/mosaicdrizzle_image_pipeline/'
filters=['wfc3uvis_f225w',\
         'wfc3uvis_f275w',\
         'wfc3uvis_f336w',\
         'wfc3uvis_f390w',\
         'acs_f435w',\
         'acs_f475w',\
         'acs_f606w',\
         'acs_f625w',\
         'acs_f775w',\
         'acs_f814w',\
         'acs_f850lp',\
         'wfc3ir_f105w',\
         'wfc3ir_f110w',\
         'wfc3ir_f125w',\
         'wfc3ir_f140w',\
         'wfc3ir_f160w']

# jphotfilename = roidata+'.cat'
# jphotfile=open(jphotfilename,'w')
# jphothead=\
# '# 1 rootname\n'+\
# '# 2 redshift\n'+\
# '# 3 cluster\n'+\
# '# 4 datestamp\n'+\
# '# 5 scaledir\n'+\
# 
# '# 1 filter\n'+\
# '# 2 ABmag_zpt\n'+\
# '# 3 threesigmaAB\n'+\
# '# 4 threesigmafnu\n'+\
# '# 5 area\n'+\
# '# 6 fnu\n'+\
# '# 7 fnuerr\n'+\
# '# 8 fnuerr_shotnoise\n'+\
# '# 9 fnuerr_skyrms\n'+\
# '# 10 ABmag\n'+\
# '# 11 ABmagerr\n'

roiinfo   = asciidata.open(roidata)
for index in range(roiinfo.nrows): 
    rootname  = str(roiinfo['rootname'][index])
    cluster   = str(roiinfo['cluster'][index])
    datestamp = str(roiinfo['datestamp'][index])
    scaledir  = str(roiinfo['scaledir'][index])
    redshift  = roiinfo['redshift'][index]

    dir=topdir+cluster+subdir+scaledir+'/'

    regionfile=rootname+'.reg'
    if not (os.path.isfile(regionfile)):
        print "object region file not found"
        sys.exit()
    print regionfile, dir

    photometryfile=rootname+'_'+scaledir+'_'+datestamp+'_phot.cat'
    photfile = open(photometryfile,'w')
    photfile.write(phothead)

    for band in (filters):
        if (scaledir=='scale_65mas'):   mosaicstring = '_mosaic_065mas_'
        elif (scaledir=='scale_30mas'): mosaicstring = '_mosaic_030mas_'
        elif (scaledir=='scale_native'):
            if 'acs' in band:           mosaicstring='_mosaic_050mas_'
            elif 'wfc3uvis' in band:    mosaicstring='_mosaic_040mas_'
            elif 'wfc3if' in band:      mosaicstring='_mosaic_130mas_'

        drzname=dir+cluster+mosaicstring+band+'_drz_'+datestamp+'.fits'
        if not (os.path.isfile(drzname)):
            drzname=drzname+'.gz'
            if not os.path.isfile(drzname):
                print "problem finding the drz file for "+band
                continue
        whtname=dir+cluster+mosaicstring+band+'_wht_'+datestamp+'.fits'
        if not (os.path.isfile(whtname)):
            whtname=whtname+'.gz'
            if not os.path.isfile(whtname):
                print "problem finding the wht file for "+band
                continue

        #        print drzname
        #        print whtname

        drz       = pyfits.open(drzname)
        wht       = pyfits.open(whtname)
        exptime   = drz[0].header['EXPTIME']
        gain      = drz[0].header['CCDGAIN']
        photflam  = drz[0].header['PHOTFLAM']
        photplam  = drz[0].header['PHOTPLAM']
        ABMAG_ZPT = -2.5*np.log10(photflam)-21.10-5*np.log10(photplam)+18.6921

        rwcs=pyregion.open(regionfile)
        rcoo=pyregion.open(regionfile).as_imagecoord(drz[0].header)
        mask=rwcs.get_mask(hdu=drz[0])
        data=drz[0].data
        weight=wht[0].data
        maskdata = data[(mask)]
        maskweight = weight[(mask)]
        area = maskdata.size

        # We calculate the rms from the weight map, which accounts for
        # everything except shot noise from the object.
        wrms = np.sqrt(np.sum(1.0/maskweight[np.where(maskweight>0.0)]))

        # units are electrons per second 
        countflux = np.sum(data[(mask)])

        # poisson error = sqrt(electrons), converted back to native
        # electrons/s units. Only assign this uncertainty if the
        # object is detected at an arbitrary 5-sigma or above.
        if (countflux<5.0*wrms): counterror=0.0
        else:                    counterror = np.sqrt(countflux*exptime)/exptime

        fnucofactor = pow((ABMAG_ZPT+23.9),-0.4)
        fnu    = fnucofactor*countflux
        fnuerr = fnucofactor * np.sqrt(counterror**2 + wrms**2)

        fnuerr_shotnoise = fnucofactor*counterror
        fnuerr_skyrms    = fnucofactor*wrms

        if (countflux<=0.0):
            ABmag=-99.0
            ABmagerr=-99.0
        else:
            ABmag = -2.5*np.log10(countflux)+ABMAG_ZPT
            ABmagerr = 2.5/np.log(10.0) * fnuerr/fnu
            
        threesigmafnu = 3.0*wrms
        threesigmaAB = -2.5*np.log10(3.0*wrms)+ABMAG_ZPT

        datastring = "%s %.3f %.4f %.4f %.1f %.4f %.4f %.4f %.4f %.3f %.3f" % \
            (band, ABMAG_ZPT, threesigmaAB, threesigmafnu, area, fnu, fnuerr, fnuerr_shotnoise, fnuerr_skyrms, ABmag, ABmagerr)
        print datastring
        photfile.write(datastring+'\n')
        drz.close()
        wht.close()

    photfile.close()

# phothead = \
#    '# 1 filter\n'+\
#    '# 2 area\n'+\
#    '# 3 countflux\n'+\
#    '# 4 ABmag_zpt\n'+\
#    '# 5 threesigma\n'+\
#    '# 6 ABmag\n'+\
#    '# 7 ABmag_unc\n'
    
"""
Drizzled data are in "counts" per second, which could be DN or electrons.
In our instrument and filter set, all counts are electrons.
Multiply DRZ images (which are in counts per second) by exposure time. 
"""
# Will want to loop over all entries in the setup / startup file
#for index in range(roiinfo.nrows):

        #        area = mask[(mask)].size

# g.write('#!/bin/csh -f\n\n')
# g.write('date\n')
# g.write('if (-e '+logfile+') then\n')
# g.write('  /bin/rm '+logfile+'\n')
# g.write('endif\n')
# g.write('touch '+logfile+'\n')
# g.close()

# +  wsigma) / countflux
#                        (2.5/np.log(10.0)*counterror)**2 +\
#                        (2.5/np.log(10.0)*wsigma/)**2 )

# logstart = "%04d" % (ScriptStart)
# logfile  = outputfilename+'_'+logstart+'.log'
# shfile   = outputfilename+'_'+logstart+'.sh'
# 


"""
The sky measurements should not be needed if the weight (inverse rms) maps are used. 
        swcs=pyregion.open(skyfile)
        skyarea=mask[(mask)].size
        skymask=swcs.get_mask(hdu=drz[0])
        skydata=data[(skymask)]
        skymean=np.mean(skydata)
        skyrms = np.sqrt(np.sum((skydata-skymean)**2)/skyarea)
"""

    # wsigma = sqrt(wrms) = sqrt(weight rms) = sqrt(quadrature sum
    # of equivalent single-pixel noise), Casertano+ 2000, Eq1&2
    #        wsigma = np.sqrt(np.sum(1.0/weight[(mask)]))
    #        photfnu = drz[0].header['PHOTFNU']

    #        \
    #            np.sqrt(area * skyrms**2 + \
    #                    countflux/gain)/countflux
#        Fnu = 10**(0.4*(23.9-ABmag))


    #    skyname   = str(roiinfo['skyname'][index])
    #    skyfile=skyname+'.reg'
    #    if not (os.path.isfile(regionfile) or os.path.isfile(skyfile)):
    #        print "object and corresponding sky region file missing"
    #        sys.exit()
