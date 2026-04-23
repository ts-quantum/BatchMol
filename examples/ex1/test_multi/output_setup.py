# created with BatchMol 4.3 (Multi-File Orbitals Frozen)
# ==============================================================================
# USER GUIDE: 
# 1. RUN THIS SCRIPT
# 2. CLEANUP: Search 'Renderer Node' in Outliner -> Select all (A) 
#    -> Right Click -> 'Delete Hierarchy'. This removes Cameras but keeps Meshes.
# ==============================================================================

import bpy, os, re

# --- Settings ---
path_to_glb = "/Users/user/python/BatchMol/examples/ex1/test_multi"
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
    
    # Track imported objects
    new_objs = [o for o in bpy.context.selected_objects]
    # Sort them by their internal Blender name (mesh_0, mesh_1...)
    new_objs.sort(key=lambda o: [int(c) if c.isdigit() else c.lower() for c in re.split('(\\d+)', o.name)])
    
    current_frame = i + 1 

    # We only care about Meshes for materials and animation
    meshes_in_file = [o for o in new_objs if o.type == 'MESH']
    
    for slot_idx, obj in enumerate(meshes_in_file):
        # A) Material Assignment (Slot Logic)
        m_names = ["MASTER_Molecule", "MASTER_Orb_Pos", "MASTER_Orb_Neg", "MASTER_Surface"]
        if slot_idx < len(m_names):
            mat = bpy.data.materials.get(m_names[slot_idx])
            if mat:
                obj.data.materials.clear()
                obj.data.materials.append(mat)
                obj.color = mat.diffuse_color

        # B) Parenting & Animation
        obj.parent = cntrl
        
        # Keyframe Sequence
        obj.scale = (0,0,0)
        obj.keyframe_insert(data_path="scale", frame=current_frame - 1)
        
        # Show only real geometry (Dummies stay hidden)
        if len(obj.data.vertices) > 1:
            obj.scale = (1,1,1)
        obj.keyframe_insert(data_path="scale", frame=current_frame)
        
        obj.scale = (0,0,0)
        obj.keyframe_insert(data_path="scale", frame=current_frame + 1)
        
        # C) Robust Constant Interpolation
        if obj.animation_data and obj.animation_data.action:
            action = obj.animation_data.action
            if hasattr(action, "fcurves"):
                for fc in action.fcurves:
                    for kp in fc.keyframe_points:
                        kp.interpolation = 'CONSTANT'

    # 4. Cleanup Hierarchy: Remove empty containers (keeps Meshes due to parenting)
    for o in new_objs:
        if o.name in bpy.data.objects:
            if o.type == 'EMPTY' or not o.data:
                bpy.data.objects.remove(o, do_unlink=True)

# Finalize
bpy.context.scene.frame_start = 1
bpy.context.scene.frame_end = len(frame_files)
bpy.context.scene.frame_set(1)

print(f"Multi-File Sync finished: {len(frame_files)} frames ready.")
