# PCA Cloud Generator

This repository contains a Python script to generate 3D plot frames from Principal Component Analysis (PCA) data, typically derived from GROMACS trajectory projections. These frames can then be compiled into an animation to visualize the conformational landscape explored during a simulation.

## Table of Contents

- [PCA Cloud Generator](#pca-cloud-generator)
  - [Table of Contents](#table-of-contents)
  - [1. `generate_3d_plot_frames.py`](#1-generate_3d_plot_framespy)
    - [Purpose](#purpose)
    - [Author & Version](#author--version)
    - [Key Features](#key-features)
    - [Dependencies](#dependencies)
    - [Usage](#usage)
    - [Output](#output)
  - [Workflow](#workflow)
  - [License](#license)

## 1. `generate_3d_plot_frames.py`

### Purpose

This Python script generates a series of 3D plot frames from PCA (Principal Component Analysis) data, typically derived from GROMACS trajectory projections. These frames can then be compiled into an animation to visualize the conformational landscape explored during a simulation. It highlights the trajectory's movement against a backdrop of all sampled conformations.

### Author & Version

*   **Author**: Markus Williams
*   **Date**: 27/05/2025
*   **Version**: 1.0

### Key Features

*   **Input**: Reads GROMACS XVG files containing PCA projections. It can handle files where multiple projections (e.g., PC1, PC2, PC3) are concatenated and separated by an '&' symbol.
*   **3D Visualization**: Plots trajectory frames in 3D space using the first three principal components.
*   **Background Cloud**: Displays a transparent "cloud" of all data points from a provided XVG file, representing the overall conformational space.
*   **Frame Skipping**: Allows processing of every Nth frame to manage large trajectory datasets and reduce animation length.
*   **Ghost Frames**: Optionally shows a specified number of preceding frames (ghosts) with decreasing opacity to visualize recent history of the trajectory.
*   **Trailing Path**: Can draw a line showing the path taken by the trajectory up to the current frame.
*   **Customization**: Offers command-line arguments to control output directory, image DPI, cloud transparency, and ghost/trail effects.
*   **Output**: Generates individual PNG frames for the animation.

### Dependencies

*   Python 3
*   NumPy
*   Matplotlib

### Usage

```bash
# Example 1: Basic usage with frame skipping
python generate_3d_plot_frames.py --cloud proj_g2_background_cloud.xvg --traj proj_g2rep1_anim.xvg --skip 5 --outdir pca_frames

# Example 2: With ghost frames
python generate_3d_plot_frames.py --cloud proj_g2_background_cloud.xvg --traj proj_g2rep1_anim.xvg --ghost 10 --outdir pca_frames
```

**Command-line Arguments:**

*   `--cloud`: (Required) XVG file for the background conformational cloud.
*   `--traj`: (Required) XVG file for the trajectory to be animated.
*   `--skip`: (Optional) Process every Nth frame. Default: 5.
*   `--outdir`: (Optional) Directory to save output frames. Default: `pca_frames`.
*   `--dpi`: (Optional) DPI for output images. Default: 150.
*   `--cloud-alpha`: (Optional) Transparency for cloud points. Default: 0.05.
*   `--trail`: (Optional) Enable trailing path of the trajectory. Enabled by default.
*   `--no-trail`: (Optional) Disable trailing path.
*   `--ghost`: (Optional) Number of ghost frames to show. Default: 0 (disabled).

### Output

The script outputs:
1.  A series of PNG image files (e.g., `pca_frame_0000.png`, `pca_frame_0001.png`, etc.) in the specified output directory.
2.  A suggestion for an `ffmpeg` command to compile these frames into an MP4 video:
    ```bash
    ffmpeg -framerate 10 -i <outdir>/pca_frame_%04d.png -c:v libx264 -pix_fmt yuv420p pca_animation.mp4
    ```

## Workflow

1.  Perform PCA on your molecular dynamics trajectory (e.g., using GROMACS tools like `gmx covar` and `gmx anaeig`). Obtain XVG files for:
    *   Projections of all trajectory frames onto the principal components (for the background cloud).
    *   Projections of a specific trajectory (or part of it) onto the principal components (for the animated trajectory).
2.  Use `generate_3d_plot_frames.py` to create PNG frames visualizing the PCA trajectory.
    ```bash
    python generate_3d_plot_frames.py --cloud <cloud_data.xvg> --traj <trajectory_data.xvg> --outdir pca_frames [options]
    ```
3.  Compile the generated frames into a video using `ffmpeg` or a similar tool.
    ```bash
    ffmpeg -framerate 10 -i pca_frames/pca_frame_%04d.png -c:v libx264 -pix_fmt yuv420p pca_animation.mp4
    ```

## License

This project is licensed under the MIT License. 
See the [LICENSE](LICENSE) file for details.