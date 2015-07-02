#!/usr/bin/env python

"""Insert artificial (fake) galaxies into the DECaLS imaging and reprocess
through The Tractor.

Questions/issues:
* Should the flux be drawn from a uniform or power-law distribution? 

TODO:
* Model a mixture of object types.

python decals_simulations.py -n 10 --zoom 1750 1850 1750 1850
runbrick --stage image_coadds --no-write --threads 16 --gpsf --radec 242.845 11.75 --width 500 --height 500
runbrick --no-sdss --no-wise --threads 8 --zoom 1800 1900 1800 1900

decals_simulations -b 2428p117 -n 10 -o STAR --zoom 1800 1900 1800 1900
runbrick -b 2428p117 --decals-dir /Users/ioannis/research/projects/decals/fake_decals_dir --stage image_coadds --no-write --threads 8 --zoom 1800 1900 1800 1900

decals_simulations -b 2428p117 -n 500 -o STAR --rmag-range 18 25 --zoom 1600 2000 1600 2000 

setenv DECALS_DIR /Users/ioannis/research/projects/decals/fake_decals_dir
runbrick -b 2428p117 --stage image_coadds --no-write --threads 8 --zoom 1600 2000 1600 2000 

"""
from __future__ import division, print_function

import os
import sys
import shutil
import logging
import argparse
import subprocess
import numpy as np

import galsim
from astropy.io import fits

from projects.desi.common import Decals, wcs_for_brick, ccds_touching_wcs

# Set up logging and our global directories.
logging.basicConfig(format='%(message)s',level=logging.INFO,stream=sys.stdout)
log = logging.getLogger('decals_simulations')

def get_simdir(brickname=None,objtype=None):
    """Get the simulation directory."""

    # Check for the environment variables we need.
    tractor_dir = 'TRACTOR_DIR'
    if tractor_dir not in os.environ:
        log.error('Missing ${} environment variable'.format(tractor_dir))
        sys.exit(1)
    else:
        tractor_dir = os.getenv('TRACTOR_DIR')

    decals_dir = 'DECALS_DIR'
    if decals_dir not in os.environ:
        log.error('Missing ${} environment variable'.format(decals_dir))
        sys.exit(1)
    else:
        decals_dir = os.getenv('DECALS_DIR')

    if 'DECALS_SIM_DIR' in os.environ:
        decals_sim_dir = os.getenv('DECALS_SIM_DIR')
    else:
        log.error('Missing $DECALS_SIM_DIR environment variable')
        sys.exit(1)

    log.info('DECALS_DIR {}'.format(decals_dir))
    log.info('DECALS_SIM_DIR {}'.format(decals_sim_dir))
    simdir = os.path.join(decals_sim_dir,objtype.lower())

    log.info('Creating directories...')
    if not os.path.isdir(simdir):
        log.info('  {}'.format(simdir))
        os.makedirs(simdir)
        
    simdir = os.path.join(simdir,brickname)
    if not os.path.isdir(simdir):
        log.info('  {}'.format(simdir))
        os.makedirs(simdir)
    
    #qadir = os.path.join(simdir,'qaplots')
    #if not os.path.isdir(qadir):
    #    log.info('  {}'.format(qadir))
    #    os.makedirs(qadir)
    
    return tractor_dir, decals_dir, simdir

def get_brickinfo(brickname=None,decals_dir=None):
    """Get info on this brick.

    """
    allbrickinfo = fits.getdata(os.path.join(decals_dir,'decals-bricks.fits'),1)
    brickinfo = allbrickinfo[np.where((allbrickinfo['brickname']==brickname)*1==1)]

    decals = Decals()
    wcs = wcs_for_brick(decals.get_brick_by_name(brickname))

    return brickinfo, wcs

def get_ccdinfo(brickwcs=None,decals_dir=None):
    """Get info on this brick and on the CCDs touching it.

    """
    allccdinfo =  fits.getdata(os.path.join(decals_dir,'decals-ccds-zeropoints.fits'))

    # Get all the CCDs that touch this brick
    decals = Decals()
    these = ccds_touching_wcs(brickwcs,decals.get_ccds())
    #ccdinfo = decals.ccds_touching_wcs(targetwcs)
    ccdinfo = allccdinfo[these]

    log.info('Got {} CCDs'.format(len(ccdinfo)))
    return ccdinfo

def build_simcat(nobj=None,brickname=None,brickwcs=None,objtype=None,
                 raminmax=None,decminmax=None,rmag_range=None,
                 decals_sim_dir=None,seed=None):
    """Build the simulated object catalog, which depends on the type of object.
       Will eventually generalize this so that a mixture of object types can be
       simulated.

    TODO (@moustakas): Remove simulated sources which are too near to one
    another.

    """
    from astropy.table import Table, Column
    rand = np.random.RandomState(seed=seed)

    # Assign central coordinates uniformly
    ra = rand.uniform(raminmax[0],raminmax[1],nobj)
    dec = rand.uniform(decminmax[0],decminmax[1],nobj)
    xxyy = brickwcs.radec2pixelxy(ra,dec)

    cat = Table()
    cat['ID'] = Column(np.arange(nobj,dtype='i4'))
    cat['OBJTYPE'] = Column(np.repeat(objtype,nobj),dtype='S10')
    cat['X'] = Column(xxyy[1][:],dtype='f4')
    cat['Y'] = Column(xxyy[2][:],dtype='f4')
    cat['RA'] = Column(ra,dtype='f8')
    cat['DEC'] = Column(dec,dtype='f8')

    if objtype.upper()=='ELG':
        sersicn_1_range = [1.0,1.0]
        r50_1_range = [0.5,2.5]
        ba_1_range = [0.2,1.0]

        sersicn_1 = rand.uniform(sersicn_1_range[0],sersicn_1_range[1],nobj)
        r50_1 = rand.uniform(r50_1_range[0],r50_1_range[1],nobj)
        ba_1 = rand.uniform(ba_1_range[0],ba_1_range[1],nobj)
        phi_1 = rand.uniform(0,180,nobj)

        cat['SERSICN_1'] = Column(sersicn_1,dtype='f4')
        cat['R50_1'] = Column(r50_1,dtype='f4')
        cat['BA_1'] = Column(ba_1,dtype='f4')
        cat['PHI_1'] = Column(phi_1,dtype='f4')

        ## Bulge parameters
        #bulge_r50_range = [0.1,1.0]
        #bulge_n_range = [3.0,5.0]
        #bdratio_range = [0.0,1.0] # bulge-to-disk ratio

        # Magnitudes and colors
        gr_range = [-0.3,0.5]
        rz_range = [0.0,1.5]

    if objtype.upper()=='STAR':
        gr_range = [0.0,0.5]
        rz_range = [0.0,1.5]

    # For convenience, also store the grz fluxes in nanomaggies.
    rmag = rand.uniform(rmag_range[0],rmag_range[1],nobj)
    gr = rand.uniform(gr_range[0],gr_range[1],nobj)
    rz = rand.uniform(rz_range[0],rz_range[1],nobj)

    gflux = 1E9*10**(-0.4*(rmag+gr))
    rflux = 1E9*10**(-0.4*rmag)
    zflux = 1E9*10**(-0.4*(rmag-rz))

    # Pack into a Table.
    cat['R'] = Column(rmag,dtype='f4')
    cat['GR'] = Column(gr,dtype='f4')
    cat['RZ'] = Column(rz,dtype='f4')
    cat['GFLUX'] = Column(gflux,dtype='f4')
    cat['RFLUX'] = Column(rflux,dtype='f4')
    cat['ZFLUX'] = Column(zflux,dtype='f4')

    # Write out.
    outfile = os.path.join(decals_sim_dir,'simcat-'+brickname+'-'+objtype.lower()+'.fits')
    log.info('Writing {}'.format(outfile))
    if os.path.isfile(outfile):
        os.remove(outfile)
    cat.write(outfile)

    return cat

def copy_cpdata(ccdinfo=None,decals_dir=None,decals_sim_dir=None):
    """Copy the CP-processed images, inverse variance maps, and bad-pixel masks we
    need from DECALS_DIR to DECALS_SIM_DIR, creating directories as necessary.
    Also construct the requisite soft links pointing back to DECALS_DIR.

    """
    from distutils.file_util import copy_file

    allcpimage = ccdinfo['CPIMAGE']
    allcpdir = set([cpim.split('/')[1] for cpim in allcpimage])

    log.info('Creating symbolic links...')
    linkfiles = ['calib','decals-ccds.fits','decals-bricks.fits',
                 'decals-ccds-zeropoints.fits']
    for lfile in linkfiles:
        log.info('  {}'.format(os.path.join(decals_sim_dir,lfile)))
        os.symlink(os.path.join(decals_dir,lfile),os.path.join(decals_sim_dir,lfile))
    
    log.info('Creating directories...')
    outdir = os.path.join(decals_sim_dir,'images')
    if not os.path.isdir(outdir):
        log.info('  {}'.format(outdir))
        os.makedirs(outdir)
    outdir = os.path.join(decals_sim_dir,'images','decam')
    if not os.path.isdir(outdir):
        log.info('  {}'.format(outdir))
        os.makedirs(outdir)

    for cpdir in list(allcpdir):
        outdir = os.path.join(decals_sim_dir,'images','decam',cpdir)
        if not os.path.isdir(outdir):
            log.info('  {}'.format(outdir))
            os.makedirs(outdir)
    
    log.info('Copying files...')
    for cpimage in list(set(allcpimage)):
        cpdir = cpimage.split('/')[1]
        indir = os.path.join(decals_dir,'images','decam',cpdir)
        outdir = os.path.join(decals_sim_dir,'images','decam',cpdir)

        imfile = cpimage.split('/')[2].split()[0]
        log.info('  {}'.format(imfile))
        copy_file(os.path.join(indir,imfile),os.path.join(outdir,imfile),update=0)

        imfile = imfile.replace('ooi','oow')
        #log.info('{}-->{}'.format(os.path.join(indir,imfile),os.path.join(outdir,imfile)))
        copy_file(os.path.join(indir,imfile),os.path.join(outdir,imfile),update=0)

        imfile = imfile.replace('oow','ood')
        #log.info('{}-->{}'.format(os.path.join(indir,imfile),os.path.join(outdir,imfile)))
        copy_file(os.path.join(indir,imfile),os.path.join(outdir,imfile),update=0)

class build_stamp():
    def __init__(self,objtype):
        """Build stamps of different types of objects.
        
        """        
        self.objtype = objtype

    def getlocal(self,objinfo,siminfo):
        self.pos = siminfo.wcs.toImage(galsim.CelestialCoord(objinfo['RA']*galsim.degrees,
                                                             objinfo['DEC']*galsim.degrees))
        self.xpos = int(self.pos.x)
        self.ypos = int(self.pos.y)
        self.offset = galsim.PositionD(self.pos.x-self.xpos,self.pos.y-self.ypos)

        # Get the WCS and PSF at the center of the stamp and the integrated
        # flux of the object (in the appropriate band).
        localwcs, pixscale = siminfo.getlocalwcs(image_pos=self.pos)
        self.localwcs = localwcs
        self.pixscale = pixscale
        self.localpsf = siminfo.getlocalpsf(image_pos=self.pos,pixscale=self.pixscale)
        self.objflux = siminfo.getobjflux(objinfo)

    def convolve_and_draw(self,obj):
        """Convolve the object with the PSF and then draw it."""
        obj = galsim.Convolve([obj,self.localpsf])
        stamp = obj.drawImage(offset=self.offset,wcs=self.localwcs,method='no_pixel')
        stamp.setCenter(self.xpos,self.ypos)
        return stamp

    def addnoise(self,stamp,varstamp,siminfo):
        varstamp.invertSelf()            # [ADU^2]
        medvar = np.median(varstamp.array[varstamp.array>0])
        varstamp.array[varstamp.array<(0.2*medvar)] = medvar

        # Convert to electrons
        stamp *= siminfo.gain         # [electron]
        varstamp *= (siminfo.gain**2) # [electron^2]

        stamp.addNoise(galsim.VariableGaussianNoise(galsim.BaseDeviate(),varstamp))
        varstamp += stamp
                
        stamp /= siminfo.gain          # [ADU]
        varstamp /= (siminfo.gain**2)  # [ADU^2]
        varstamp.invertSelf()          # [1/ADU^2]

        return stamp, varstamp

    def star(self):
        """Create a PSF source."""
        psf = self.localpsf
        psf = psf.withFlux(self.objflux)
        stamp = psf.drawImage(offset=self.offset,wcs=self.localwcs,method='no_pixel')
        stamp.setCenter(self.xpos,self.ypos)
        return stamp
    
    def elg(self,objinfo,siminfo):
        """Create an ELG (disk-like) galaxy."""
        obj = galsim.Sersic(float(objinfo['SERSICN_1']),half_light_radius=
                            float(objinfo['R50_1']),
                            flux=self.objflux,gsparams=siminfo.gsparams)
        obj = obj.shear(q=float(objinfo['BA_1']),beta=
                        float(objinfo['PHI_1'])*galsim.degrees)
        stamp = self.convolve_and_draw(obj)
        return stamp

    def lrg(self,objinfo,siminfo):
        """Create an LRG (spheroidal) galaxy."""
        obj = galsim.Sersic(float(objinfo['SERSICN_1']),half_light_radius=
                            float(objinfo['R50_1']),
                            flux=self.objflux,gsparams=siminfo.gsparams)
        obj = obj.shear(q=float(objinfo['BA_1']),beta=
                        float(objinfo['PHI_1'])*galsim.degrees)
        stamp = self.convolve_and_draw(obj)
        return stamp

class simobj_info():
    def __init__(self,ccdinfo,gsparams,decals_sim_dir):
        """Access everything we need about an individual CCD.
        
        """        
        from tractor.psfex import PsfEx
        from tractor.basics import GaussianMixtureEllipsePSF

        self.PsfEx = PsfEx
        self.GaussianMixtureEllipsePSF = GaussianMixtureEllipsePSF
        self.gsparams = gsparams
        self.cpimage = ccdinfo['CPIMAGE']
        self.cpimage_hdu = ccdinfo['CPIMAGE_HDU']
        self.calname = ccdinfo['CALNAME']
        self.filter = ccdinfo['FILTER']
        self.imfile = os.path.join(decals_sim_dir,'images',self.cpimage)
        self.imfile_root = self.cpimage.replace('decam/','')
        self.ivarfile = self.imfile.replace('ooi','oow')
        self.wcsfile = os.path.join(decals_sim_dir,'calib','decam',
                                    'astrom-pv',self.calname+'.wcs.fits')
        self.psffile = os.path.join(decals_sim_dir,'calib','decam',
                                    'psfex',self.calname+'.fits')
        self.magzpt = float(ccdinfo['CCDZPT'] + 2.5*np.log10(ccdinfo['EXPTIME']))
        self.gain = float(ccdinfo['ARAWGAIN']) # [electron/ADU]

    def getdata(self):
        """Read the CCD image and inverse variance data, and the corresponding headers. 
        
        """
        #log.info('Reading extension {} of image {}'.format(self.ccdnum,self.imfile))
        image = galsim.fits.read(self.imfile,hdu=self.cpimage_hdu)       # [ADU]
        invvar = galsim.fits.read(self.ivarfile,hdu=self.cpimage_hdu) # [1/ADU^2]
        
        imhdr = galsim.fits.FitsHeader(self.imfile,hdu=self.cpimage_hdu)
        ivarhdr = galsim.fits.FitsHeader(self.ivarfile,hdu=self.cpimage_hdu)
        
        self.width = image.xmax
        self.height = image.ymax
        
        return image, invvar, imhdr, ivarhdr

    def getwcs(self):
        """Read the global WCS for this CCD."""
        wcs, origin = galsim.wcs.readFromFitsHeader(
            galsim.fits.FitsHeader(self.wcsfile))
        self.wcs = wcs

        return wcs
    
    def getpsf(self):
        """Read the PSF for this CCD."""
        psf = self.PsfEx(self.psffile,self.width,self.height,ny=13,nx=7,
                          psfClass=self.GaussianMixtureEllipsePSF,K=2)
        self.psf = psf

        return psf

    def getlocalwcs(self,image_pos=None):
        """Get the local WCS, given a position."""
        localwcs = self.wcs.local(image_pos=image_pos)
        pixscale, shear, theta, flip = localwcs.getDecomposition() # get the pixel scale

        return localwcs, pixscale

    def getlocalpsf(self,image_pos=None,pixscale=0.262):
        """Get the local PSF, given a position.  Need to recentroid because this is a
        PSFeX PSF."""
        xpos = int(image_pos.x)
        ypos = int(image_pos.y)
        psfim = self.PsfEx.instantiateAt(self.psf,xpos,ypos)[5:-5,5:-5] # trim
        psf = galsim.InterpolatedImage(galsim.Image(psfim),scale=pixscale,flux=1.0)
        psf_centroid = psf.centroid()
        psf = psf.shift(-psf_centroid.x,-psf_centroid.y)
        return psf
    
    def getobjflux(self,objinfo):
        """Calculate the flux of a given object in ADU."""
        flux = objinfo[self.filter.upper()+'FLUX']
        flux *= 10**(0.4*(self.magzpt-22.5)) # [ADU]
        return float(flux)
    
def insert_simobj(simcat,ccdinfo,decals_sim_dir):
    """Simulate objects and place them into individual CCDs."""

    gsparams = galsim.GSParams(maximum_fft_size=2L**30L)

    width = int(ccdinfo['WIDTH'][0])
    height = int(ccdinfo['HEIGHT'][0])

    stampwidth = 45 # postage stamp width [pixels, roughly 14 arcsec]
    stampbounds = galsim.BoundsI(-stampwidth,stampwidth,-stampwidth,stampwidth)
    imagebounds = galsim.BoundsI(0,width,0,height)

    objtype = simcat['OBJTYPE'][0].upper()
    objstamp = build_stamp(objtype)
    
    for ccd in ccdinfo:
        # Gather some basic info on this CCD and then read the data, the WCS
        # info, and initialize the PSF.
        siminfo = simobj_info(ccd,gsparams,decals_sim_dir)
        wcs = siminfo.getwcs()

        # Loop on each object and figure out which, if any, objects will be
        # placed on this CCD.
        onccd = []
        for iobj, objinfo in enumerate(simcat):
            pos = wcs.toImage(galsim.CelestialCoord(objinfo['RA']*galsim.degrees,
                objinfo['DEC']*galsim.degrees))
            stampbounds1 = stampbounds.shift(galsim.PositionI(int(pos.x),int(pos.y)))
        
            overlap = stampbounds1 & imagebounds
            #if iobj<5:
            #    print(iobj, pos, stampbounds1, imagebounds, overlap)
            if (overlap.xmax>=0 and overlap.ymax>=0 and overlap.xmin<=imagebounds.xmax and
                overlap.ymin<=imagebounds.ymax and overlap.area()>0):
                onccd.append(iobj)

        nobj = len(onccd)
        if nobj>0:
            log.info('Adding {} objects to HDU {}'.format(nobj,siminfo.cpimage_hdu))

            image, invvar, imhdr, ivarhdr = siminfo.getdata()
            initpsf = siminfo.getpsf()
            
            for iobj in range(nobj):
                #print(iobj)
                objinfo = simcat[onccd[iobj]]

                # get the local coordinate, WCS, and PSF and then build the stamp
                objstamp.getlocal(objinfo,siminfo)
                if objtype=='STAR':
                    stamp = objstamp.star()
                if objtype=='ELG':
                    stamp = objstamp.elg(objinfo,siminfo)

                overlap = stamp.bounds & image.bounds
                if (overlap.xmax>=0 and overlap.ymax>=0 and
                    overlap.xmin<=image.bounds.xmax and
                    overlap.ymin<=image.bounds.ymax and overlap.area()>0):

                    # Add Poisson noise
                    stamp = stamp[overlap]            # [ADU]
                    varstamp = invvar[overlap].copy() # [1/ADU^2]

                    stamp, varstamp = objstamp.addnoise(stamp,varstamp,siminfo)
                    image[overlap] += stamp
                    invvar[overlap] = varstamp

            log.info('Writing {}[{}]'.format(siminfo.imfile_root,siminfo.cpimage_hdu))
            fits.update(siminfo.imfile,image.array,ext=siminfo.cpimage_hdu,
                        header=fits.Header(imhdr.items()))
            fits.update(siminfo.ivarfile,invvar.array,ext=siminfo.cpimage_hdu,
                        header=fits.Header(ivarhdr.items()))

def qaplots(brickinfo,ccdinfo,simcat,decals_sim_dir=None):
    """Build some simple QAplots of the simulation inputs."""
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    from matplotlib.patches import Rectangle

    brickname = brickinfo['BRICKNAME'][0]
    color = iter(cm.rainbow(np.linspace(0,1,len(ccdinfo))))

    fig = plt.figure()
    ax = fig.gca()
    ax.get_xaxis().get_major_formatter().set_useOffset(False) 
    ax.plot(simcat['RA'],simcat['DEC'],'gs',markersize=3)
    for ii, ccd in enumerate(ccdinfo):
        dy = ccd['WIDTH']*0.262/3600.0
        dx = ccd['HEIGHT']*0.262/3600.0
        rect = plt.Rectangle((ccd['RA']-dx/2,ccd['DEC']-dy/2),
                             dx,dy,fill=False,lw=1,color=next(color),
                             ls='solid')
        ax.add_patch(rect)
        rect = plt.Rectangle((brickinfo['RA1'],brickinfo['DEC1']),
                             brickinfo['RA2']-brickinfo['RA1'],
                             brickinfo['DEC2']-brickinfo['DEC1'],fill=False,lw=3,
                             color='b')
        ax.add_patch(rect)
        ax.set_xlim(np.array([brickinfo['RA2'][0],brickinfo['RA1'][0]])*[1.0002,0.9998])
        ax.set_ylim(np.array([brickinfo['DEC1'][0],brickinfo['DEC2'][0]])*[0.985,1.015])
        ax.set_xlabel('$RA\ (deg)$',fontsize=18)
        ax.set_ylabel('$Dec\ (deg)$',fontsize=18)

    qafile = os.path.join(decals_sim_dir,'qa_'+brickname+'_ccds.png')
    #qafile = os.path.join(decals_sim_dir,'qaplots','qa_'+brickname+'_ccds.png')
    log.info('Writing QAplot {}'.format(qafile))
    fig.savefig(qafile)

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.
                                     ArgumentDefaultsHelpFormatter,
                                     description='DECaLS simulations.')
    parser.add_argument('-n', '--nobj', type=long, default=None, metavar='', 
                        help='number of objects to simulate (required input)')
    parser.add_argument('-b', '--brick', type=str, default='2428p117', metavar='', 
                        help='simulate objects in this brick')
    parser.add_argument('-o', '--objtype', type=str, default='ELG', metavar='', 
                        help='object type (STAR, ELG, LRG, BGS)') 
    parser.add_argument('-t', '--threads', type=int, default=8, metavar='', 
                        help='number of threads to use when calling The Tractor')
    parser.add_argument('-s', '--seed', type=long, default=None, metavar='', 
                        help='random number seed')
    parser.add_argument('--zoom', nargs=4, type=int, metavar='', 
                        help='see runbrick.py (default is to populate the full brick)')
    parser.add_argument('--rmag-range', nargs=2, type=float, default=(18,25), metavar='', 
                        help='r-band magnitude range')
    parser.add_argument('--do-tractor', action='store_true',
                        help='run the Tractor')
    parser.add_argument('--no-qaplots', action='store_true',
                        help='do not generate QAplots')

    args = parser.parse_args()
    if args.nobj is None:
        parser.print_help()
        sys.exit(1)

    nobj = args.nobj
    brickname = args.brick
    objtype = args.objtype.upper()

    if objtype=='LRG':
        log.warning('{} objtype not yet supported!'.format(objtype))
        sys.exit(1)
    elif objtype=='BGS':
        log.warning('{} objtype not yet supported!'.format(objtype))
        sys.exit(1)

    log.info('Simulating {} objects of objtype {} in brick {}'.
             format(nobj,objtype,brickname))

    # Build the output directory names.
    tractor_dir, decals_dir, decals_sim_dir = get_simdir(brickname,objtype)

    # Get the brick info and corresponding WCS object.
    brickinfo, brickwcs = get_brickinfo(brickname,decals_dir)

    if args.zoom is None:
        raminmax = [brickinfo['ra1'],brickinfo['ra2']]
        decminmax = [brickinfo['dec1'],brickinfo['dec2']]
    else:
        pixscale = 0.262/3600.0 # average pixel scale [deg/pixel]
        zoom = args.zoom
        dx = zoom[1]-zoom[0]
        dy = zoom[3]-zoom[2]

        ra, dec = brickwcs.pixelxy2radec(zoom[0]+dx/2,zoom[2]+dy/2)
        raminmax = [ra-dx*pixscale/2,ra+dx*pixscale/2]
        decminmax = [dec-dy*pixscale/2,dec+dy*pixscale/2]

        brickwcs = brickwcs.get_subimage(zoom[0],zoom[2],dx,dy)

    # Identify the CCDs in the region of interest.
    ccdinfo = get_ccdinfo(brickwcs,decals_dir)

    log.info('RA range: {:.6f} to {:.6f}'.format(float(raminmax[0]),
                                                 float(raminmax[1])))
    log.info('DEC range: {:.6f} to {:.6f}'.format(float(decminmax[0]),
                                                  float(decminmax[1])))
        
    # Build the simulated object catalog and optionally make some QAplots.
    log.info('Building the simulated object catalog')
    simcat = build_simcat(nobj,brickname,brickwcs,objtype,raminmax,
                          decminmax,rmag_range=args.rmag_range,
                          decals_sim_dir=decals_sim_dir,seed=args.seed)
    if args.no_qaplots is False:
        qaplots(brickinfo,ccdinfo,simcat,decals_sim_dir)
        
    # Copy the CP-processed data we need to DECALS_SIM_DIR.
    copy_cpdata(ccdinfo,decals_dir,decals_sim_dir)

    # Do the simulation!
    insert_simobj(simcat,ccdinfo,decals_sim_dir)

    # Call The Tractor
    if args.do_tractor:
        env = os.environ.copy()
        env['DECALS_DIR'] = decals_sim_dir
        runbrick = 'python {}/projects/desi/runbrick.py'.format(tractor_dir)
        if args.zoom is None:
            cmd = '{} -b {} --threads {} --no-sdss --no-wise'.format(
                runbrick,brickname,args.threads).split()
        else:
            cmd = '{} -b {} --threads {} --no-sdss --no-wise --zoom {} {} {} {}'.format(
                runbrick,brickname,args.threads,zoom[0],zoom[1],zoom[2],zoom[3]).split()
        log.info(cmd)
        ok = subprocess.check_call(cmd,cwd=decals_sim_dir,env=env)

        # Copy the tractor catalogs and coadded images
        if ok==0:
            linkfiles = ['calib','decals-ccds.fits','decals-bricks.fits',
                         'decals-ccds-zeropoints.fits']
            for lfile in linkfiles:
                os.unlink(os.path.join(decals_sim_dir,lfile))

            shutil.move(os.path.join(decals_sim_dir,'tractor',brickname[:3],'tractor-'+brickname+'.fits'),
                        os.path.join(decals_sim_dir,'tractor-'+brickname+'.fits'))
            shutil.move(os.path.join(decals_sim_dir,'coadd',brickname[:3],brickname,
                                     'decals-'+brickname+'-image.jpg'),
                        os.path.join(decals_sim_dir,'decals-'+brickname+'-image.jpg'))
            shutil.move(os.path.join(decals_sim_dir,'coadd',brickname[:3],brickname,
                                     'decals-'+brickname+'-resid.jpg'),
                        os.path.join(decals_sim_dir,'decals-'+brickname+'-resid.jpg'))
            shutil.move(os.path.join(decals_sim_dir,'coadd',brickname[:3],brickname,
                                     'decals-'+brickname+'-model.jpg'),
                        os.path.join(decals_sim_dir,'decals-'+brickname+'-model.jpg'))

            # Clean up files, symbolic links, and directories
            shutil.rmtree(os.path.join(decals_sim_dir,'images'))
            shutil.rmtree(os.path.join(decals_sim_dir,'coadd'))
            shutil.rmtree(os.path.join(decals_sim_dir,'metrics'))
            shutil.rmtree(os.path.join(decals_sim_dir,'pickles'))
            shutil.rmtree(os.path.join(decals_sim_dir,'tractor'))

            # Write a log file!

            # Analyze the outputs.


if __name__ == "__main__":
    main()
