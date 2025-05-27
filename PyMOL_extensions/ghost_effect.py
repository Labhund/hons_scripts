from pymol import cmd

def create_ghosting_animation(object_prefix, max_transparency=0.8, num_offsets=None):
    """
    PyMOL command to create a ghosting animation effect for objects with a common prefix.

    Modifies the states of objects with a common prefix to create a ghosting delay effect.
    For objects with offsets, the first N states (where N is the offset value) will be
    empty, followed by the states from the reference object (offset0).
    
    Also applies a transparency gradient from 0.0 (reference object) to max_transparency
    (highest offset object).

    If offset objects don't exist, they can be created automatically.

    Args:
        object_prefix (str): The common prefix for the object names.
                             Example: "g2_rep1_strided_offset" for objects like
                             "g2_rep1_strided_offset0", "g2_rep1_strided_offset1", etc.
        max_transparency (float, optional): Maximum transparency value (0.0-1.0) to apply
                                           to the highest offset object. Default is 0.8.
        num_offsets (int or str, optional): Number of offset objects to create if they don't exist.
                                           If None, only works with existing objects.
    """
    # Convert parameters to the correct types
    try:
        max_transparency = float(max_transparency)
    except (ValueError, TypeError):
        print(f"Warning: Invalid max_transparency value '{max_transparency}'. Using default 0.8.")
        max_transparency = 0.8
        
    if num_offsets is not None:
        try:
            num_offsets = int(num_offsets)
        except (ValueError, TypeError):
            print(f"Warning: Invalid num_offsets value '{num_offsets}'. Only using existing objects.")
            num_offsets = None

    all_objects = cmd.get_names("objects")
    matching_objects_temp = [obj for obj in all_objects if obj.startswith(object_prefix)]

    # Check if we need to create offset objects
    if (not matching_objects_temp or len(matching_objects_temp) == 1) and num_offsets is not None and num_offsets > 0:
        print(f"Creating {num_offsets} offset objects from existing object...")
        
        # Find a source object to duplicate
        source_obj = None
        if matching_objects_temp:
            source_obj = matching_objects_temp[0]
        else:
            # If no object with the prefix exists, check if the prefix itself is an object
            if object_prefix in all_objects:
                source_obj = object_prefix
                
        if not source_obj:
            print(f"Error: No source object found to create offsets. Need either:")
            print(f"  - An object named '{object_prefix}'")
            print(f"  - At least one object starting with '{object_prefix}'")
            return
            
        # Create offset0 if it doesn't exist
        if not any(obj.startswith(f"{object_prefix}0") for obj in all_objects):
            offset0_name = f"{object_prefix}0"
            print(f"Creating reference object {offset0_name} from {source_obj}")
            cmd.create(offset0_name, source_obj)
            
        # Create additional offset objects
        for i in range(1, num_offsets + 1):
            offset_name = f"{object_prefix}{i}"
            if not any(obj == offset_name for obj in all_objects):
                print(f"Creating offset object {offset_name}")
                cmd.create(offset_name, f"{object_prefix}0")
                
        # Refresh the list of objects after creation
        all_objects = cmd.get_names("objects")
        matching_objects_temp = [obj for obj in all_objects if obj.startswith(object_prefix)]

    if not matching_objects_temp:
        print(f"No objects found with prefix '{object_prefix}'.")
        return

    parsed_objects = []
    for obj_name in matching_objects_temp:
        try:
            offset_str = obj_name[len(object_prefix):]
            if offset_str.isdigit() or (offset_str.startswith('-') and offset_str[1:].isdigit()):
                offset = int(offset_str)
                num_states = cmd.count_states(obj_name)
                if num_states > 0:
                    parsed_objects.append({"name": obj_name, "offset": offset, "num_states": num_states})
                else:
                    print(f"Warning: Object {obj_name} has 0 states. Skipping.")
        except ValueError:
            print(f"Warning: Could not parse numeric offset from suffix of {obj_name}. Skipping.")
            
    if not parsed_objects:
        print(f"No objects with a valid numeric suffix found for prefix '{object_prefix}'.")
        return

    # Sort by offset to process them in order
    parsed_objects.sort(key=lambda x: x["offset"])
    
    # Find the reference object (offset0)
    reference_obj = next((p_obj for p_obj in parsed_objects if p_obj["offset"] == 0), None)
    
    if not reference_obj:
        print(f"Error: Reference object with offset 0 (e.g., '{object_prefix}0') not found.")
        return
        
    reference_name = reference_obj["name"]
    reference_states = reference_obj["num_states"]
    
    print(f"Using {reference_name} with {reference_states} states as reference object.")
    
    # Find the maximum offset value for transparency calculation
    max_offset = max(p_obj["offset"] for p_obj in parsed_objects)
    
    # Process each object with offset > 0
    for p_obj in parsed_objects:
        obj_name = p_obj["name"]
        obj_offset = p_obj["offset"]
        
        # Calculate transparency based on offset
        if max_offset > 0:  # Avoid division by zero
            # For offset0, transparency will be 0.0
            # For highest offset, transparency will be max_transparency
            transparency = (obj_offset / max_offset) * max_transparency
        else:
            transparency = 0.0
            
        # Apply transparency to the object - using only valid PyMOL transparency settings
        cmd.set("cartoon_transparency", transparency, obj_name)
        cmd.set("stick_transparency", transparency, obj_name)
        cmd.set("sphere_transparency", transparency, obj_name)
        cmd.set("ribbon_transparency", transparency, obj_name)
        # For surfaces, use transparency setting (not surface_transparency)
        cmd.set("transparency", transparency, obj_name)
        
        print(f"Set transparency for {obj_name} to {transparency:.2f}")
        
        # Skip state modification for reference object or negative offsets
        if obj_offset <= 0:
            continue
            
        print(f"Processing {obj_name} (offset: {obj_offset})...")
        
        # Create temporary object to rebuild with proper state offsets
        tmp_obj_name = f"__tmp_{obj_name}"
        cmd.create(tmp_obj_name, "none")  # Create empty object
        
        # First, add empty states up to the offset value
        for i in range(1, obj_offset + 1):
            # Create an empty state by copying with a selection that matches nothing
            cmd.create(f"{tmp_obj_name}", "none", 0, i)
            print(f"  Created empty state {i}/{obj_offset}")
            
        # Then add states from the reference object
        for i in range(1, reference_states + 1):
            target_state = obj_offset + i
            cmd.create(f"{tmp_obj_name}", f"{reference_name}", i, target_state)
            if i % 50 == 0 or i == reference_states:
                print(f"  Copied reference state {i} to target state {target_state}")
                
        # Delete the original object and rename the temporary one
        cmd.delete(obj_name)
        cmd.set_name(tmp_obj_name, obj_name)
        
        # Re-apply transparency to the renamed object
        cmd.set("cartoon_transparency", transparency, obj_name)
        cmd.set("stick_transparency", transparency, obj_name)
        cmd.set("sphere_transparency", transparency, obj_name)
        cmd.set("ribbon_transparency", transparency, obj_name)
        cmd.set("transparency", transparency, obj_name)
        
        print(f"Rebuilt {obj_name} with {obj_offset} empty states followed by states from {reference_name}.")
    
    print("Ghosting effect applied to all objects.")
    print("Objects modified with transparency gradient:")
    for p_obj in parsed_objects:
        transparency = (p_obj["offset"] / max_offset) * max_transparency if max_offset > 0 else 0.0
        print(f"  - {p_obj['name']} (offset: {p_obj['offset']}, transparency: {transparency:.2f})")

# Make this function available as a command in PyMOL
cmd.extend("ghost_animate", create_ghosting_animation)
# Usage examples:
# ghost_animate g2_rep1_strided_offset
# ghost_animate g2_rep1_strided_offset, 0.7
# ghost_animate g2_rep1_strided_offset, 0.8, 5  # Create 5 offset objects if they don't exist