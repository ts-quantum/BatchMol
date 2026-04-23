# created with BatchMol 4.3 (ESP Multi-File Fix)
import bpy, os, re

path_to_glb = "/Users/user/python/BatchMol/examples/ex1/bld-esp-multi"
protected = ["Camera", "Plane", "Cylinder", "Sun", "World", "TRAJECTORY_CONTROL", "MASTER", "DUMMY"]

# 1. Controller Setup
if "TRAJECTORY_CONTROL" not in bpy.data.objects:
    cntrl = bpy.data.objects.new("TRAJECTORY_CONTROL", None)
    bpy.context.scene.collection.objects.link(cntrl)
else:
    cntrl = bpy.data.objects["TRAJECTORY_CONTROL"]

# 2. File Discovery
files = sorted([f for f in os.listdir(path_to_glb) if f.endswith(".glb") and "scalebar" not in f])

# 3. Import Loop
for i, filename in enumerate(files):
    filepath = os.path.join(path_to_glb, filename)
    bpy.ops.import_scene.gltf(filepath=filepath)
    
    # WICHTIG: Wir sortieren die Objekte INNERHALB der Datei (0=Mol, 1=Surf)
    new_objs = [o for o in bpy.context.selected_objects if o.type == 'MESH']
    new_objs.sort(key=lambda o: [int(c) if c.isdigit() else c.lower() for c in re.split('(\\d+)', o.name)])
    
    current_frame = i + 1 

    for slot_idx, obj in enumerate(new_objs):
        # --- A) POSITIONING ---
        old_matrix = obj.matrix_world.copy()
        obj.parent = cntrl
        obj.matrix_world = old_matrix

        # --- B) MATERIAL LOGIC (Slot-Zuweisung innerhalb des Frames) ---
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
                esp_mat.blend_method = 'BLEND'
                
                nodes = esp_mat.node_tree.nodes
                links = esp_mat.node_tree.links
                bsdf = nodes.get("Principled BSDF")
                
                if bsdf:
                    bsdf.inputs['Alpha'].default_value = 0.5
                    bsdf.inputs['Metallic'].default_value = 0.3
                    bsdf.inputs['Roughness'].default_value = 0.2
                    
                    # Link Base Color to Emission (Rainbow Glow)
                    if bsdf.inputs['Base Color'].is_linked:
                        source_socket = bsdf.inputs['Base Color'].links[0].from_socket
                        links.new(source_socket, bsdf.inputs['Emission Color'])
                    
                    if 'Emission Strength' in bsdf.inputs:
                        bsdf.inputs['Emission Strength'].default_value = 1.5

        # --- C) ANIMATION ---
        obj.scale = (0, 0, 0)
        obj.keyframe_insert(data_path="scale", frame=current_frame - 1)
        obj.scale = (1, 1, 1) if len(obj.data.vertices) > 1 else (0, 0, 0)
        obj.keyframe_insert(data_path="scale", frame=current_frame)
        obj.scale = (0, 0, 0)
        obj.keyframe_insert(data_path="scale", frame=current_frame + 1)
        
        # Constant Interpolation
        if obj.animation_data and obj.animation_data.action:
            if hasattr(obj.animation_data.action, "fcurves"):
                for fc in obj.animation_data.action.fcurves:
                    for kp in fc.keyframe_points: kp.interpolation = 'CONSTANT'
    
    # Cleanup Junk Nodes (nur Empties/Nodes löschen)
    for o in bpy.context.selected_objects:
        if o.type != 'MESH':
            bpy.data.objects.remove(o, do_unlink=True)

# 4. Final Scene Sync
bpy.context.scene.frame_start = 1
bpy.context.scene.frame_end = len(files)
bpy.context.scene.frame_set(1)
