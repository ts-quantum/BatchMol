# created with BatchMol 4.3 (One-File ESP Hybrid Pro)
import bpy, re, os

def natural_key(text):
    return [int(c) if c.isdigit() else c.lower() for c in re.split('(\\d+)', text)]

# 1. Controller Setup
if "TRAJECTORY_CONTROL" not in bpy.data.objects:
    cntrl = bpy.data.objects.new("TRAJECTORY_CONTROL", None)
    bpy.context.scene.collection.objects.link(cntrl)
else:
    cntrl = bpy.data.objects["TRAJECTORY_CONTROL"]

# 2. Discovery & Sorting
protected = ["Camera", "Plane", "Sun", "World", "TRAJECTORY_CONTROL", "MASTER", "DUMMY", "STATIC_CB"]
container = bpy.data.objects.get("Renderer Node")

if container:
    all_objs = [child for child in container.children if child.type == 'MESH']
else:
    all_objs = [o for o in bpy.data.objects if o.type == 'MESH' and not any(p in o.name for p in protected)]

all_objs.sort(key=lambda o: natural_key(o.name))

# 3. Processing Slots (0=Mol, 1=D, 2=D, 3=Surf)
slots_per_frame = 2
for i, obj in enumerate(all_objs):
    frame_idx = (i // slots_per_frame) + 1
    slot_idx = i % slots_per_frame 
    
    obj.parent = cntrl
    
    # --- Material-Logik (Hybrid Pro) ---
    if slot_idx == 0:
        # Molecule -> Link to MASTER_Molecule
        mat = bpy.data.materials.get("MASTER_Molecule")
        if mat:
            obj.data.materials.clear()
            obj.data.materials.append(mat)
            obj.color = mat.diffuse_color
            
    elif slot_idx == 1 and len(obj.data.polygons) > 10:
        # ESP-Surface: Preserve Colors + Refine Properties
        if obj.data.materials:
            esp_mat = obj.data.materials[0]
            esp_mat.use_nodes = True
            esp_mat.blend_method = 'BLEND'  # Enable Transparency
            
            nodes = esp_mat.node_tree.nodes
            links = esp_mat.node_tree.links
            bsdf = nodes.get("Principled BSDF")
            
            if bsdf:
                # Adjust properties (Alpha, Metallic, Roughness)
                bsdf.inputs['Alpha'].default_value = 0.5
                bsdf.inputs['Metallic'].default_value = 0.3
                bsdf.inputs['Roughness'].default_value = 0.2
                
                # Link Base Color to Emission to make the rainbow glow
                if bsdf.inputs['Base Color'].is_linked:
                    source_socket = bsdf.inputs['Base Color'].links[0].from_socket
                    links.new(source_socket, bsdf.inputs['Emission Color'])
                
                if 'Emission Strength' in bsdf.inputs:
                    bsdf.inputs['Emission Strength'].default_value = 1.5

    # 4. Animation (Visibility via Scale)
    is_dummy = len(obj.data.polygons) == 0
    obj.scale = (0, 0, 0)
    obj.keyframe_insert(data_path="scale", frame=frame_idx - 1)
    obj.scale = (1, 1, 1) if not is_dummy else (0, 0, 0)
    obj.keyframe_insert(data_path="scale", frame=frame_idx)
    obj.scale = (0, 0, 0)
    obj.keyframe_insert(data_path="scale", frame=frame_idx + 1)

    # Constant Interpolation
    if obj.animation_data and obj.animation_data.action:
        act = obj.animation_data.action
        if hasattr(act, "fcurves"):
            for fc in act.fcurves:
                for kp in fc.keyframe_points: kp.interpolation = 'CONSTANT'

# 5. Timeline & Viewport
bpy.context.scene.frame_end = len(all_objs) // slots_per_frame
bpy.context.scene.frame_set(1)

for area in bpy.context.screen.areas:
    if area.type == 'VIEW_3D':
        for space in area.spaces:
            if space.type == 'VIEW_3D': space.shading.color_type = 'OBJECT'

print(f"Hybrid ESP Setup complete. Alpha and Emission applied to {len(all_objs)//4} frames.")
