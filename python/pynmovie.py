#!/usr/bin/python2.7
import sys
from optparse import OptionParser
import pynbody.analysis.angmom as angmom
import pynbody as pyn
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
#import pynbody.plot as p_plt

if __name__ == "__main__":
	parser = OptionParser()
	parser.add_option("-r", "--range", action="store", dest="valrange", 
                          help="The range of values to plot (in the format 'min max')")
	parser.add_option("-f", "--fps", action="store", dest="fps", 
                          help="Frame Rate", default=5)
	parser.add_option("-t", "--type", action="store", dest="ptype", 
                          help="Type of particles to show (gas, dm, or stars)", default="stars")
	parser.add_option("-o", "--outname", action="store", dest="outname", 
                          help="Filename of the video to be produced", default="pynmovie.mp4")
	parser.add_option("-n", action="store_false", dest="video", 
                          help="Don't actually produce a video, just make pngs.", default=True)
	(opts, args) = parser.parse_args()
	imgcount = 0
	try:
                (vmin, vmax) = [float(i) for i in opts.valrange.split()]
        except AttributeError:
                (vmin, vmax) = (None, None)
                for i in args:
                        # load the snapshot and retrieve the xyz coordinates
                        sim = pyn.load(i)
                        xyz_dm = sim.dm['pos']
                        xyz_stars = sim.stars['pos']
                        xyz_gas = sim.gas['pos']
                        
                        # build the multi-panel plot                
                        
                        plt.figure(1)
                        # xy, yz, xz for the stars
                        plt.subplot(2,3,1)
                        plt.hist2d(xyz_stars[:,0],xyz_stars[:,1],bins=50,cmap='Greys',norm=LogNorm())
#                       plt.plot(xyz_stars[:,0],xyz_stars[:,1],'yo',markersize=3)
                        plt.axis([-75,75,-75,75])
                        plt.ylabel('Separation (kpc)')
                        plt.subplot(2,3,2)
                        plt.hist2d(xyz_stars[:,1],xyz_stars[:,2],bins=50,cmap='Greys',norm=LogNorm())
#                       plt.plot(xyz_stars[:,1],xyz_stars[:,2],'yo',markersize=3)
                        plt.axis([-75,75,-75,75])
                        plt.subplot(2,3,3)
                        plt.hist2d(xyz_stars[:,0],xyz_stars[:,2],bins=50,cmap='Greys',norm=LogNorm())
#                       plt.plot(xyz_stars[:,0],xyz_stars[:,2],'yo',markersize=3)
                        plt.axis([-75,75,-75,75])

                        # xy, yz, xz for the gas
                        plt.subplot(2,3,4)
                        plt.hist2d(xyz_gas[:,0],xyz_gas[:,1],bins=50,cmap='Greys',norm=LogNorm())
#                       pyn.plot.sph.image(sim.gas,qty='pos',width=75,cmap='Greys',norm=LogNorm())
#                       plt.plot(xyz_gas[:,0],xyz_gas[:,1],'ro',markersize=2)
#                       plt.axis([-75,75,-75,75])
                        plt.ylabel('Separation (kpc)')
                        plt.xlabel('Separation (kpc)')
                        plt.subplot(2,3,5)
                        plt.hist2d(xyz_gas[:,1],xyz_gas[:,2],bins=50,cmap='Greys',norm=LogNorm())
#                       plt.plot(xyz_gas[:,1],xyz_gas[:,2],'ro',markersize=2)
                        plt.axis([-75,75,-75,75])
                        plt.xlabel('Separation (kpc)')
                        plt.subplot(2,3,6)
#                       plt.plot(xyz_gas[:,0],xyz_gas[:,2],'ro',markersize=2)
                        plt.hist2d(xyz_gas[:,0],xyz_gas[:,2],bins=50,cmap='Greys',norm=LogNorm())
                        plt.axis([-75,75,-75,75])
                        plt.xlabel('Separation (kpc)')

#                        # xy, yz, xz for the DM
#                        plt.subplot(3,3,7)
#                        plt.plot(xyz_dm[:,0],xyz_dm[:,1],'bo',markersize=1)
#                        plt.axis([-75,75,-75,75])
#                        plt.ylabel('Separation (kpc)')
#                        plt.xlabel('Separation (kpc)')
#                        plt.subplot(3,3,8)
#                        plt.plot(xyz_dm[:,1],xyz_dm[:,2],'bo',markersize=1)
#                        plt.axis([-75,75,-75,75])
#                        plt.xlabel('Separation (kpc)')
#                        plt.subplot(3,3,9)
#                        plt.plot(xyz_dm[:,0],xyz_dm[:,2],'bo',markersize=1)
#                        plt.axis([-75,75,-75,75])
#                        plt.xlabel('Separation (kpc)')

# close the figure
                        plt.savefig("%09d.png" % (imgcount), dpi=150)
                        plt.close(1)
                        imgcount += 1
                        
                if opts.video:
                        import envoy
                        print 'ffmpeg -y -qscale 1 -r %d -i %%09d.png %s' % (int(opts.fps), opts.outname)
                        vid = envoy.run('ffmpeg  -qscale 1 -r %d -i %%09d.png %s' % (int(opts.fps), opts.outname))
#                       for i in range(imgcount):
#                               envoy.run("rm %09d.png" % (i))
