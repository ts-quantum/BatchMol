# created with BatchMol 4.3
import bpy, os, re

path_to_glb = "/Users/user/python/BatchMol/examples/ex1/bld-homo-mulit"
protected = ["Camera", "Plane", "Cylinder", "Sun", "World", "TRAJECTORY_CONTROL", "MASTER", "DUMMY"]

# 1. Setup Controller
if "TRAJECTORY_CONTROL" not in bpy.data.objects:
    cntrl = bpy.data.objects.new("TRAJECTORY_CONTROL", None)
    bpy.context.collection.objects.link(cntrl)
else:
    cntrl = bpy.data.objects["TRAJECTORY_CONTROL"]

# 2. File Discovery
files = sorted([f for f in os.listdir(path_to_glb) if f.endswith(".glb")])
frame_files = [f for f in files if "scalebar" not in f]

# 3. Import Loop
for i, filename in enumerate(frame_files):
    filepath = os.path.join(path_to_glb, filename)
    bpy.ops.import_scene.gltf(filepath=filepath)
    
    # All objects from THIS file
    new_objs = [o for o in bpy.context.selected_objects]
    # Sort them by their internal Blender name (mesh_0, mesh_1...)
    new_objs.sort(key=lambda o: [int(c) if c.isdigit() else c.lower() for c in re.split('(\\d+)', o.name)])
    
    current_frame = i + 1 

    # We only care about Meshes for materials and animation
    meshes_in_file = [o for o in new_objs if o.type == 'MESH']
    
    for slot_idx, obj in enumerate(meshes_in_file):
        # 1. Material Assignment (Slot Logic)
        # Based on your observation: 0=Mol, 1=Orb_Pos, 2=Orb_Neg, 3=Surface
        m_names = ["MASTER_Molecule", "MASTER_Orb_Pos", "MASTER_Orb_Neg", "MASTER_Surface"]
        if slot_idx < len(m_names):
            mat = bpy.data.materials.get(m_names[slot_idx])
            if mat:
                obj.data.materials.clear()
                obj.data.materials.append(mat)
                obj.color = mat.diffuse_color

        # 2. Parenting & Animation
        obj.parent = cntrl
        obj.scale = (0,0,0)
        obj.keyframe_insert(data_path="scale", frame=current_frame - 1)
        
        # Hide dummies (1 vertex)
        if len(obj.data.vertices) > 1:
            obj.scale = (1,1,1)
        obj.keyframe_insert(data_path="scale", frame=current_frame)
        
        obj.scale = (0,0,0)
        obj.keyframe_insert(data_path="scale", frame=current_frame + 1)
        
        # Force CONSTANT interpolation (using a more robust access method)
        anim_data = obj.animation_data
        if anim_data and anim_data.action:
            action = anim_data.action
            if hasattr(action, "fcurves"):
                for fc in action.fcurves:
                    for kp in fc.keyframe_points:
                        kp.interpolation = 'CONSTANT'

    # 4. Cleanup Hierarchy: Remove the empty "Camera/Renderer Nodes" from this import
    for o in new_objs:
        if o.type == 'EMPTY' or not o.data:
            bpy.data.objects.remove(o, do_unlink=True)

bpy.context.scene.frame_end = len(frame_files)
bpy.context.scene.frame_set(1)
print("Multi-File Import & Master-Link complete.")
