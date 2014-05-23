#Script by Graziano Vernizzi (Siena College - 2014)
#Based on Brian R. Kent paper "Visualizing Astronomical Data with Blender"
#
# load the blender module first, then  
# use the command:
# blender -b -P BlenderGadgetGV.py

import sys
import pynbody as pyn
from optparse import OptionParser
import bpy
import bmesh

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
            sim = pyn.load(i)
            ptypes = {"gas": sim.gas, "dm": sim.dm, "stars": sim.stars}
            xyz_dm = sim.dm['pos']
            xyz_stars = sim.stars['pos']
            xyz_gas = sim.gas['pos']

#           xyz = (ptypes[opts.ptype])['pos']
#           plt.plot(xyz[:,0],xyz[:,1],'bo')
#           p_plt.image(ptypes[opts.ptype]) #, units='m_p cm**-3', cmap="jet", vmin=vmin, vmax=vmax)

# initialize blender
            bpy.ops.object.select_all()
            bpy.ops.object.select_all()
            bpy.ops.object.delete(use_global=False)   
            bpy.ops.mesh.primitive_circle_add(radius=1, location=(0,0,0))
            obj= bpy.data.objects['Circle']
            bpy.ops.object.mode_set(mode='EDIT')
            bpy.ops.mesh.merge(type='CENTER')
            bm=bmesh.from_edit_mesh(obj.data)

# stars            
            bm.verts.new(xyz_stars)

            bmesh.update_edit_mesh(obj.data)
            bpy.ops.object.mode_set(mode='OBJECT')

            mat = bpy.data.materials.new('HALO')
            mat.type= 'HALO'
            mat.halo.size=0.7
            txt=bpy.data.textures.new('yel','BLEND')
            mtex=mat.texture_slots.add()    
            mtex.color = (1, 0.788799, 0.193239)
            mtex.texture=txt
            
            bpy.context.object.data.materials.append(mat)
            bpy.context.scene.world.horizon_color=(0,0,0)

            bpy.ops.object.camera_add(view_align=True, enter_editmode=False, location=(200,0,400), rotation=(12.0/180*3.14,24.0/180*3.14,2
9.0/180*3.14))
            bpy.context.object.data.clip_end = 1000
            bpy.data.objects['Camera'].select=True

            bpy.data.scenes['Scene'].camera=bpy.data.objects[0]
            bpy.data.scenes['Scene'].render.filepath = "%09d.png" % (imgcount)
            bpy.ops.render.render( write_still=True ) 
            
            imgcount += 1
            if opts.video:
                import envoy
                print 'ffmpeg  -qscale 1 -r %d -i %%09d.png %s' % (int(opts.fps), opts.outname)
                vid = envoy.run('ffmpeg  -qscale 1 -r %d -i %%09d.png %s' % (int(opts.fps), opts.outname))
                #		for i in range(imgcount):
                #			envoy.run("rm %09d.png" % (i))
                
