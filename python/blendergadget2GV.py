#Script by Graziano Vernizzi (Siena College - 2014)
#Based on Brian R. Kent paper "Visualizing Astronomical Data with Blender"
#
# load the blender module first, then  
# use the command:
# blender -b -P BlenderGadgetGV.py


import bpy
import bmesh

Smin=0
Smax=30+1

#modify this with your own path
path='//home//vernizzi1//scenes//'

for i in range(Smin, Smax):
    print(i)
    
    bpy.ops.object.select_all()
    bpy.ops.object.select_all()
    bpy.ops.object.delete(use_global=False)   
    bpy.ops.mesh.primitive_circle_add(radius=1, location=(0,0,0))
    obj= bpy.data.objects['Circle']
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.merge(type='CENTER')
    bm=bmesh.from_edit_mesh(obj.data)


    filename=path+'scene_'+str(i)+'.txt'
    fin=open(filename,"r")
    a=fin.readlines()
    fin.close()
    L=len(a)
    print(L)


    count=0
    scale=1
    for l in a:
        l1=l.strip();        
        l2=l1.split(',');   
        x=float(l2[0])      
        y=float(l2[1])
        z=float(l2[2])
        bm.verts.new((x,y,z))


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

    bpy.ops.object.camera_add(view_align=True, enter_editmode=False, location=(200,0,400), rotation=(12.0/180*3.14,24.0/180*3.14,29.0/180*3.14))
    bpy.context.object.data.clip_end = 1000
    bpy.data.objects['Camera'].select=True
    



    bpy.data.scenes['Scene'].camera=bpy.data.objects[0]
    bpy.data.scenes['Scene'].render.filepath = 'scenes//scene_'+str(i)+'.jpg'
    
    bpy.ops.render.render( write_still=True ) 
