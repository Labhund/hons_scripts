# PyMOL Extensions

This directory contains custom Python scripts that extend the functionality of PyMOL.

## Table of Contents

- [PyMOL Extensions](#pymol-extensions)
  - [Table of Contents](#table-of-contents)
  - [1. `ghost_effect.py` (`ghost_animate` command)](#1-ghost_effectpy-ghost_animate-command)
    - [Purpose](#purpose)
    - [Key Features](#key-features)
    - [Installation / Making the Command Available](#installation--making-the-command-available)
    - [Usage (PyMOL Command)](#usage-pymol-command)
    - [Dependencies](#dependencies)
  - [License](#license)

## 1. `ghost_effect.py` (`ghost_animate` command)

This script defines a new PyMOL command `ghost_animate`.

### Purpose

The `ghost_animate` command creates a "ghosting" or "trailing" effect for molecular animations in PyMOL. It manipulates the states and transparency of a series of related molecular objects (e.g., different time points or replicas from a simulation) to make previous states linger with increasing transparency. This is particularly useful for visualizing dynamic processes, conformational changes, or ensemble movements, providing a clearer view of the trajectory or history of motion.

This effect can be used to complement other visualization techniques, such as PCA trajectory plots, by providing a corresponding molecular view of the system's dynamics.

### Key Features

*   **Object Targeting**: Works on PyMOL objects that share a common prefix and have a numeric suffix indicating their offset (e.g., `protein_offset0`, `protein_offset1`, `protein_offset2`). The script processes these objects to create the ghosting sequence.
*   **State Shifting**: For objects with an offset `N > 0`, the script modifies them so their first `N` states are empty. Subsequent states are then populated by copying states from a reference object (typically the one with offset 0). This creates the illusion of a delayed appearance.
*   **Transparency Gradient**: Applies a configurable transparency gradient across the offset objects. The reference object (offset 0) remains fully opaque, while objects representing earlier time points (higher offsets) become progressively more transparent, up to a user-defined maximum.
*   **Automatic Object Creation**: If specified, the script can automatically create the necessary offset objects by duplicating a source object. This simplifies setup if you only have a single trajectory or reference structure loaded.
*   **Flexible Control**: The command accepts parameters to customize the object prefix, the maximum transparency for the "oldest" ghost, and the number of offset objects to create if they don't already exist.
*   **Comprehensive Transparency**: Applies transparency settings to various molecular representations (cartoon, stick, sphere, ribbon, surface) for a consistent effect.

### Installation / Making the Command Available

To use the `ghost_animate` command within PyMOL:

1.  Place the `ghost_effect.py` script in a known directory.
2.  In PyMOL, source the script using the `run` command. For example, if `ghost_effect.py` is in the current working directory or a directory in PyMOL's Python path:
    ```python
    # In the PyMOL command line or a .pml script:
    run ghost_effect.py
    ```
    Alternatively, you can provide the full path to the script:
    ```python
    run /path/to/your/PyMOL_extensions/ghost_effect.py
    ```
    Once run, the `ghost_animate` command will be available for use in your PyMOL session.

### Usage (PyMOL Command)

After making the command available (see Installation), you can use `ghost_animate` as follows:

```python
# In the PyMOL command line:

# Basic usage:
# Assumes objects like 'g2_rep1_strided_offset0', 'g2_rep1_strided_offset1', etc., already exist.
# The script will identify them based on the 'g2_rep1_strided_offset' prefix.
ghost_animate g2_rep1_strided_offset

# Specify maximum transparency for the furthest ghost (e.g., 70% transparent):
ghost_animate g2_rep1_strided_offset, 0.7

# Automatically create 5 offset objects (e.g., g2_rep1_strided_offset0 through g2_rep1_strided_offset5)
# if they don't already exist. It will use an existing object named 'g2_rep1_strided_offset'
# or 'g2_rep1_strided_offset0' as a template for creating these new objects.
# The maximum transparency for g2_rep1_strided_offset5 will be 0.8.
ghost_animate g2_rep1_strided_offset, 0.8, 5
```

**Command Arguments (for `ghost_animate`):**

1.  `object_prefix` (str): The common prefix for the object names that the script will operate on.
2.  `max_transparency` (float, optional): The maximum transparency value (ranging from 0.0 for fully opaque to 1.0 for fully transparent) that will be applied to the object with the highest offset (the "oldest" ghost). Default is 0.8.
3.  `num_offsets` (int or str, optional): The number of offset objects (e.g., `prefix0`, `prefix1`, ..., `prefixN`) to manage or create. If these objects (or some of them) do not exist, and `num_offsets` is provided, the script will attempt to create them by duplicating a source object (either an object matching `object_prefix` or `object_prefix` + `0`). If `num_offsets` is not provided or set to `None`, the script will only operate on existing objects matching the prefix.

### Dependencies

*   PyMOL (with Python scripting enabled)

## License

This script is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.