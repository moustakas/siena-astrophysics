#!/usr/bin/python2.7
import sys
import pynbody
from optparse import OptionParser
import pynbody.analysis.angmom as angmom
import pynbody.plot as p_plt
import pynbody as pyn
import matplotlib.pyplot as plt

# test code

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
	parser.add_option("-e", action="store", dest="script", 
	help="preprocessing script to run on each simulation.")
	parser.add_option("-n", action="store_false", dest="video", 
	help="Don't actually produce a video, just make pngs.", default=True)
	(opts, args) = parser.parse_args()
	imgcount = 0
	try:
		(vmin, vmax) = [float(i) for i in opts.valrange.split()]
	except AttributeError:
		(vmin, vmax) = (None, None)
	if opts.script != None:
		sys.path.append('.')
		scriptfile = __import__(opts.script.split('.')[0])
	for i in args:
		if opts.script != None:
			scriptfile.process(i)
		else:
			sim = pyn.load(i)
			ptypes = {"gas":sim.gas, "dm":sim.dm, "stars":sim.stars}
			plt.figure(1)
                        xyz_dm = sim.dm['pos']
                        xyz_stars = sim.stars['pos']
                        xyz_gas = sim.gas['pos']
#                       xyz = (ptypes[opts.ptype])['pos']
                        plt.plot(xyz_stars[:,0],xyz_stars[:,1],'yo',markersize=3)
                        plt.plot(xyz_gas[:,0],xyz_gas[:,1],'ro',markersize=2)
#                       plt.plot(xyz_dm[:,0],xyz_dm[:,1],'bo',markersize=1)
			plt.axis([-250,250,-250,250])
#                       plt.plot(xyz[:,0],xyz[:,1],'bo')
#			p_plt.image(ptypes[opts.ptype]) #, units='m_p cm**-3', cmap="jet", vmin=vmin, vmax=vmax)
		plt.savefig("%09d.png" % (imgcount), dpi=150)
		plt.close(1)
		imgcount += 1
	if opts.video:
		import envoy
                print 'ffmpeg -y -qscale 1 -r %d -i %%09d.png %s' % (int(opts.fps), opts.outname)
		vid = envoy.run('ffmpeg  -qscale 1 -r %d -i %%09d.png %s' % (int(opts.fps), opts.outname))
#		for i in range(imgcount):
#			envoy.run("rm %09d.png" % (i))
