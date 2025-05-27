#!/usr/bin/env python3
"""
Script to generate 3D PCA plot frames for animation.
Takes every Nth frame of a trajectory and plots it against the first three eigenvalues,
with a transparent cloud of all points in the background.
Takes a GROMACS XVG file with multiple projections separated by '&' and combines them.
The script supports ghost frames to show previous states in the animation.
This script is designed to visualize PCA trajectories in 3D, allowing for the creation of animations
of molecular dynamics simulations or other time-series data.
It uses matplotlib for plotting and can generate frames for use in video creation.
It can also handle large datasets by skipping frames and supports customization of output parameters.
Author: Markus Williams
Date: 27/05/2025
Version: 1.0
License: MIT License

Usage:
  python generate_3d_plot_frames.py --cloud proj_g2_background_cloud.xvg --traj proj_g2rep1_anim.xvg --skip 5 --outdir pca_frames
  python generate_3d_plot_frames.py --cloud proj_g2_background_cloud.xvg --traj proj_g2rep1_anim.xvg --ghost 10 --outdir pca_frames
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
import argparse
import sys

def read_xvg_multi_proj(filename, num_expected_pcs=3):
    """
    Read a GROMACS XVG file containing multiple projection data sets separated by '&'
    and combine them into a single array with time and PCs.
    
    Args:
        filename (str): Path to the XVG file.
        num_expected_pcs (int): Number of expected principal components (default: 3)
        
    Returns:
        numpy.ndarray: Array with columns [time, PC1, PC2, PC3, ...]
    """
    print(f"Reading {filename} for multiple projections...")
    
    # Initialize lists to store data from each section
    all_sections = []
    current_section = []
    current_pc = 0
    
    # Read the file and split into sections
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            
            # Handle section delimiter
            if line == '&':
                if current_section:
                    all_sections.append(np.array(current_section))
                    print(f"  Found PC{current_pc+1} section with {len(current_section)} data points")
                    current_pc += 1
                current_section = []
                continue
                
            # Skip comment and header lines
            if line.startswith(('#', '@')) or not line:
                continue
                
            # Try to parse data lines
            try:
                values = [float(x) for x in line.split()]
                if len(values) >= 2:  # At least time and one value
                    current_section.append(values)
            except ValueError:
                print(f"  Warning: Skipping unparseable line: {line}")
    
    # Don't forget the last section
    if current_section:
        all_sections.append(np.array(current_section))
        print(f"  Found PC{current_pc+1} section with {len(current_section)} data points")
    
    # Check if we have any data
    if not all_sections:
        raise ValueError(f"No valid data found in {filename}")
    
    # Extract times from the first section (should be the same for all sections)
    time_column = all_sections[0][:, 0]
    
    # Create the combined dataset
    # First column is time, followed by one column per PC
    combined_data = np.zeros((len(time_column), 1 + num_expected_pcs))
    combined_data[:, 0] = time_column  # Time column
    
    # Add data from each section
    for i, section in enumerate(all_sections):
        if i < num_expected_pcs:
            if len(section) != len(time_column):
                print(f"  Warning: PC{i+1} has {len(section)} points but expected {len(time_column)}. Using available data.")
                min_len = min(len(section), len(time_column))
                combined_data[:min_len, i+1] = section[:min_len, 1]
            else:
                combined_data[:, i+1] = section[:, 1]
    
    num_pcs_found = len(all_sections)
    if num_pcs_found < num_expected_pcs:
        print(f"  Warning: Found {num_pcs_found} PCs but expected {num_expected_pcs}. Missing PCs will be zero.")
    elif num_pcs_found > num_expected_pcs:
        print(f"  Warning: Found {num_pcs_found} PCs but only using the first {num_expected_pcs}.")
    
    return combined_data

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Generate 3D PCA plot frames for animation')
    
    parser.add_argument('--cloud', type=str, required=True,
                        help='XVG file containing the background cloud data')
    parser.add_argument('--traj', type=str, required=True,
                        help='XVG file containing the trajectory data for animation')
    parser.add_argument('--skip', type=int, default=5,
                        help='Use every Nth frame from the trajectory (default: 5)')
    parser.add_argument('--outdir', type=str, default='pca_frames',
                        help='Output directory for frames (default: pca_frames)')
    parser.add_argument('--dpi', type=int, default=150,
                        help='DPI for output images (default: 150)')
    parser.add_argument('--cloud-alpha', type=float, default=0.05,
                        help='Alpha (transparency) for cloud points (default: 0.05)')
    parser.add_argument('--trail', action='store_true',
                        help='Show trailing path of trajectory (default: enabled)')
    parser.add_argument('--no-trail', dest='trail', action='store_false',
                        help='Disable trailing path')
    parser.add_argument('--ghost', type=int, default=0,
                        help='Number of ghost frames to show (default: 0, disabled)')
    parser.set_defaults(trail=True)
    
    return parser.parse_args()

def main():
    """Main function to generate PCA plot frames."""
    args = parse_arguments()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.outdir, exist_ok=True)
    
    # Read the cloud data (all frames)
    try:
        cloud_data = read_xvg_multi_proj(args.cloud, num_expected_pcs=3)
        # First column is time, next 3 columns are projections on the first 3 eigenvectors
        cloud_pc1 = cloud_data[:, 1]  # First principal component
        cloud_pc2 = cloud_data[:, 2]  # Second principal component
        cloud_pc3 = cloud_data[:, 3]  # Third principal component
    except Exception as e:
        print(f"Error reading cloud data: {e}")
        sys.exit(1)
    
    # Read the trajectory data for animation
    try:
        traj_data = read_xvg_multi_proj(args.traj, num_expected_pcs=3)
        # First column is time, next 3 columns are projections on the first 3 eigenvectors
        traj_pc1 = traj_data[:, 1]  # First principal component
        traj_pc2 = traj_data[:, 2]  # Second principal component
        traj_pc3 = traj_data[:, 3]  # Third principal component
    except Exception as e:
        print(f"Error reading trajectory data: {e}")
        sys.exit(1)
    
    # Determine the total number of frames and how many we'll process
    total_frames = len(traj_data)
    selected_frames = list(range(0, total_frames, args.skip))  # Every Nth frame
    print(f"Total frames in trajectory: {total_frames}")
    print(f"Selecting every {args.skip}th frame ({len(selected_frames)} frames total)")
    
    # Create a figure with specific size for better resolution
    plt.rcParams['figure.figsize'] = [12, 10]
    plt.rcParams['figure.dpi'] = args.dpi
    
    # Find min/max values for consistent axes
    pc1_min, pc1_max = min(np.min(cloud_pc1), np.min(traj_pc1)), max(np.max(cloud_pc1), np.max(traj_pc1))
    pc2_min, pc2_max = min(np.min(cloud_pc2), np.min(traj_pc2)), max(np.max(cloud_pc2), np.max(traj_pc2))
    pc3_min, pc3_max = min(np.min(cloud_pc3), np.min(traj_pc3)), max(np.max(cloud_pc3), np.max(traj_pc3))
    
    # Add a buffer for aesthetics
    buffer = 0.1
    pc1_range = pc1_max - pc1_min
    pc2_range = pc2_max - pc2_min
    pc3_range = pc3_max - pc3_min
    pc1_min -= buffer * pc1_range
    pc1_max += buffer * pc1_range
    pc2_min -= buffer * pc2_range
    pc2_max += buffer * pc2_range
    pc3_min -= buffer * pc3_range
    pc3_max += buffer * pc3_range
    
    print("Generating frames...")
    # Generate a frame for each selected frame
    for i, frame_idx in enumerate(selected_frames):
        print(f"Processing frame {frame_idx} ({i+1}/{len(selected_frames)})")
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        
        # Plot elements in order of background to foreground
        
        # 1. Plot the cloud with high transparency (background)
        cloud = ax.scatter(cloud_pc1, cloud_pc2, cloud_pc3, 
                   alpha=args.cloud_alpha,  # Reduced transparency 
                   color='lightgray',  # Lighter color for better contrast
                   s=5,       # Smaller size for less visual clutter
                   label='All conformations')
        
        # 2. Add a trailing line showing the trajectory up to the current point
        if args.trail and frame_idx > 0:
            # Plot trajectory up to current frame
            step_size = args.skip  # Match our frame selection
            history_indices = range(0, frame_idx+1, step_size)
            trail = ax.plot(traj_pc1[history_indices], 
                    traj_pc2[history_indices], 
                    traj_pc3[history_indices], 
                    color='blue', 
                    linewidth=1.5,
                    alpha=0.7)
        
        # 3. Add ghost effect if enabled (intermediate layer)
        if args.ghost > 0 and i > 0:
            # Determine how many ghost frames to show (limited by available frames)
            num_ghost_frames = min(args.ghost, i)
            
            # Calculate indices of frames to show as ghosts
            ghost_indices = []
            for j in range(1, num_ghost_frames + 1):
                ghost_frame_index = i - j
                if ghost_frame_index >= 0:
                    ghost_indices.append(selected_frames[ghost_frame_index])
            
            # Calculate alpha values for ghost frames
            # Most recent ghost: alpha=0.8, Oldest ghost: alpha=0.1
            if num_ghost_frames > 1:
                alpha_values = np.linspace(0.8, 0.1, num_ghost_frames)
            else:
                alpha_values = [0.8]
            
            # Plot ghost frames with decreasing opacity
            for j, (ghost_idx, alpha) in enumerate(zip(ghost_indices, alpha_values)):
                lag = j + 1  # How many frames behind
                ghost = ax.scatter(traj_pc1[ghost_idx], traj_pc2[ghost_idx], traj_pc3[ghost_idx], 
                          color='orange',  # Different color for ghost frames
                          s=80 - j * (30 / max(num_ghost_frames, 1)),  # Decreasing size
                          alpha=alpha,
                          edgecolors='darkorange',  # Add edge for better visibility
                          label=f'Ghost (-{lag})' if j == 0 else None)  # Only label the first ghost
        
        # 4. Plot the current frame as a fully opaque dot (foreground)
        current = ax.scatter(traj_pc1[frame_idx], traj_pc2[frame_idx], traj_pc3[frame_idx], 
                   color='red', 
                   s=120,      # Larger size for better visibility
                   alpha=1.0,  # Fully opaque
                   edgecolors='darkred',  # Add edge for better visibility
                   linewidth=1.5,
                   label=f'Frame {frame_idx}',
                   zorder=10)  # Ensure it's drawn on top
        
        # Set consistent axis limits
        ax.set_xlim(pc1_min, pc1_max)
        ax.set_ylim(pc2_min, pc2_max)
        ax.set_zlim(pc3_min, pc3_max)
        
        # Label the axes
        ax.set_xlabel('PC1', fontsize=12)
        ax.set_ylabel('PC2', fontsize=12)
        ax.set_zlabel('PC3', fontsize=12)
        ax.set_title(f'PCA Trajectory Visualization - Frame {frame_idx}', fontsize=14)
        
        # Add a legend with better positioning
        ax.legend(loc='upper right', fontsize=10)
        
        # Save the figure
        frame_filename = os.path.join(args.outdir, f'pca_frame_{i:04d}.png')
        plt.tight_layout()
        plt.savefig(frame_filename, dpi=args.dpi)
        plt.close()
    
    print(f"Done! Generated {len(selected_frames)} frames in {args.outdir}/")
    print("To create a video from these frames, you can use ffmpeg:")
    print(f"ffmpeg -framerate 10 -i {args.outdir}/pca_frame_%04d.png -c:v libx264 -pix_fmt yuv420p pca_animation.mp4")

if __name__ == "__main__":
    main()