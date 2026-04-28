# created with BatchMol 4.3 (C) 2026 by Dr. Tobias Schulz
# ==============================================================================
# USER GUIDE for BatchMol Blender Animation (One File Mode)
# ==============================================================================
# (A) Run this Script within Blender
#
# 1. GLOBAL VISUAL CONTROL: 
#    This script links all imported meshes to "MASTER" materials in your template.
#    Edit these materials in the 'Material Properties' tab to update ALL frames:
#    - 'MASTER_Molecule'  -> Controls atoms and bonds (mol_***)
#    - 'MASTER_Orb_Pos'   -> Controls positive lobes (orb_pos_***, spin_pos_***)
#    - 'MASTER_Orb_Neg'   -> Controls negative lobes (orb_neg_***, spin_neg_***)
#
# 2. RETAINING COLORS (CPK):
#    'MASTER_Molecule' uses 'Vertex Colors' (Color Attributes).
#    In the Shader Editor, ensure a 'Color Attribute' node is connected to the 
#    'Base Color' and 'Emission Color' of the Principled BSDF.
#
# 3. POSITIONING:
#    Select the 'TRAJECTORY_CONTROL' (Empty) to move, rotate, or scale the 
#    entire animation sequence simultaneously over your scene.
# ==============================================================================

import bpy, re, os

# --- Helper: Natural Sort (mesh_2 comes before mesh_10) ---
def natural_key(text):
    return [int(c) if c.isdigit() else c.lower() for c in re.split('(\\d+)', text)]

# 1. Controller Setup
if "TRAJECTORY_CONTROL" not in bpy.data.objects:
    cntrl = bpy.data.objects.new("TRAJECTORY_CONTROL", None)
    bpy.context.collection.objects.link(cntrl)
else:
    cntrl = bpy.data.objects["TRAJECTORY_CONTROL"]

# 2. Discovery and FORCE SORT
container = bpy.data.objects.get("Renderer Node")
if not container:
    all_objs = [o for o in bpy.data.objects if o.type == 'MESH' and "MASTER" not in o.name]
else:
    all_objs = [child for child in container.children if child.type == 'MESH']

# force correct order
all_objs.sort(key=lambda o: natural_key(o.name))

# 3. Processing Slots
slots_per_frame = 4
for i, obj in enumerate(all_objs):
    frame_idx = (i // slots_per_frame) + 1
    slot_idx = i % slots_per_frame 
    
    obj.parent = cntrl
    
    # Material Assignment
    m_names = ["MASTER_Molecule", "MASTER_Orb_Pos", "MASTER_Orb_Neg", "MASTER_Surface"]
    mat = bpy.data.materials.get(m_names[slot_idx])
    if mat:
        obj.data.materials.clear()
        obj.data.materials.append(mat)
        obj.color = mat.diffuse_color

    # Animation
    is_dummy = len(obj.data.polygons) == 0
    obj.scale = (0, 0, 0)
    obj.keyframe_insert(data_path="scale", frame=frame_idx - 1)
    obj.scale = (1, 1, 1) if not is_dummy else (0, 0, 0)
    obj.keyframe_insert(data_path="scale", frame=frame_idx)
    obj.scale = (0, 0, 0)
    obj.keyframe_insert(data_path="scale", frame=frame_idx + 1)

    # --- Constant Interpolation Fix ---
    if obj.animation_data and obj.animation_data.action:
        act = obj.animation_data.action
        if hasattr(act, "fcurves"):
            for fc in act.fcurves:
                for kp in fc.keyframe_points:
                    kp.interpolation = 'CONSTANT'

# 4. Viewport Fix
for area in bpy.context.screen.areas:
    if area.type == 'VIEW_3D':
        for space in area.spaces:
            if space.type == 'VIEW_3D': space.shading.color_type = 'OBJECT'

bpy.context.scene.frame_end = len(all_objs) // slots_per_frame
bpy.context.scene.frame_set(1)
print(f"Sorted {len(all_objs)} objects. Timeline ready.")
