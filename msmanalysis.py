import mdtraj as md
import pyemma
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import os
import sys
import shutil
# Needed for the multiprocessing fix, especially on Windows/WSL with 'spawn'
from multiprocess import freeze_support
import traceback
import time
import argparse

# --- Configuration ---
topology_file = "reference.pdb"
trajectory_file = "protein_only_stride10_truncated.xtc"
base_output_dir = "msm_analysis_output_k100" # Base name for output
# stride = 1 # Stride used for analysis (adjust if different from traj loading)
tica_lag = 50 # TICA lag time in frames (needs tuning based on data)
tica_dim = 10 # Number of TICA dimensions
# n_clusters = 100 # Number of clusters
its_lags = range(1, 101, 5) # Lags for implied timescale calculation (frames)
# n_metastable = 4 # Number of metastable states (needs tuning based on ITS plot)
# msm_lag = 50 # MSM lag time in frames (needs tuning based on ITS plot)
n_jobs = -1 # Use all available CPU cores for parallelizable steps (-1 for all)
n_samples_per_state = 5 # Number of representative structures per state

# --- Helper Class for Logging ---
class Logger:
    """Redirects print statements to both console and a log file."""
    def __init__(self, filepath):
        self.terminal = sys.stdout
        # Ensure directory exists before opening file
        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        self.logfile = open(filepath, 'a') # Append mode
        self.closed = False

    def write(self, message):
        self.terminal.write(message)
        self.logfile.write(message)

    def flush(self):
        self.terminal.flush()
        self.logfile.flush()

    def close(self):
        if not self.closed:
            try:
                self.logfile.close()
            except Exception as e:
                self.terminal.write(f"Warning: Error closing log file: {e}\n")
            self.closed = True

    def __enter__(self):
        sys.stdout = self
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout = self.terminal # Restore original stdout
        self.close()

# --- Helper Function for Output Directory Backup ---
def prepare_output_directory(base_dir):
    """Checks if base_dir exists. If so, renames it with a _bak suffix."""
    output_dir = base_dir
    if os.path.exists(output_dir):
        i = 0
        backup_dir = f"{output_dir}_bak"
        while os.path.exists(backup_dir):
            i += 1
            backup_dir = f"{output_dir}_bak{i}"
        print(f"Output directory '{output_dir}' already exists. Renaming to '{backup_dir}'.")
        try:
            shutil.move(output_dir, backup_dir)
            print(f"Successfully renamed '{output_dir}' to '{backup_dir}'.")
        except Exception as e:
            print(f"Error renaming existing directory: {e}. Attempting to continue...")
            # Decide if this is fatal or not. For now, we continue.
            # exit() # Uncomment to make renaming failure fatal

    # Create the new directory
    try:
        os.makedirs(output_dir, exist_ok=False) # exist_ok=False ensures it's newly created
        print(f"Created output directory: {output_dir}")
    except FileExistsError:
         print(f"Warning: Output directory '{output_dir}' still exists after backup attempt. Contents might be overwritten.")
    except Exception as e:
        print(f"Fatal error creating output directory '{output_dir}': {e}")
        traceback.print_exc()
        exit() # Exit if creation fails
    return output_dir

# --- Main Execution Block ---
if __name__ == '__main__':
    freeze_support() # Essential for multiprocessing esp. on Windows/spawn systems

 # --- Argument Parsing ---
    parser = argparse.ArgumentParser(description='Perform MSM analysis on molecular dynamics trajectory data.')
    parser.add_argument('-k', '--clusters', type=int, default=100,
                        help='Number of clusters (microstates) for KMeans. Default: 100')
    parser.add_argument('-l', '--lag', type=int, default=50,
                        help='MSM lag time in frames (should be chosen based on ITS plot). Default: 50')
    parser.add_argument('-s', '--stride', type=int, default=1,
                        help='Stride to use when loading trajectory features for analysis. Default: 1')
    parser.add_argument('-n', '--meta', type=int, default=4, dest='n_metastable', # Use dest to match internal variable name
                        help='Number of metastable states to identify with PCCA++. Default: 4')
    # Add other arguments here if needed in the future (e.g., -t for TICA lag, -d for TICA dim)
    # parser.add_argument('-t', '--tica_lag', type=int, default=50, help='TICA lag time. Default: 50')
    # parser.add_argument('-d', '--tica_dim', type=int, default=10, help='TICA dimensions. Default: 10')

    args = parser.parse_args()

    # --- Assign Parsed Arguments to Variables ---
    n_clusters = args.clusters
    msm_lag = args.lag
    stride = args.stride             # Assign stride from args
    n_metastable = args.n_metastable # Assign n_metastable from args
    # If you add args for tica_lag, tica_dim, assign them here too
    # tica_lag = args.tica_lag
    # tica_dim = args.tica_dim

    # Define base_output_dir dynamically based on n_clusters
    base_output_dir = f"msm_analysis_output_k{n_clusters}_lag{msm_lag}_stride{stride}" # Make output more descriptive

    # --- Prepare Output Directory (with backup) ---
    output_dir = prepare_output_directory(base_output_dir)
    log_file_path = os.path.join(output_dir, 'msmanalysis.log')

    # --- Start Logging ---
    with Logger(log_file_path) as logger:
        start_time = time.time()
        print(f"--- MSM Analysis Script Started: {time.strftime('%Y-%m-%d %H:%M:%S')} ---")
        print(f"--- Configuration ---")
        print(f"  Topology: {topology_file}")
        print(f"  Trajectory: {trajectory_file}")
        print(f"  Output Directory: {output_dir}")
        print(f"  Analysis Stride: {stride}")
        print(f"  TICA Lag: {tica_lag}, Dimensions: {tica_dim}")
        print(f"  Clusters (k): {n_clusters}")
        print(f"  ITS Lags: {list(its_lags)}")
        print(f"  MSM Lag: {msm_lag}")
        print(f"  Metastable States: {n_metastable}")
        print(f"  Representative Samples: {n_samples_per_state}")
        print(f"  Parallel Jobs: {n_jobs if n_jobs > 0 else 'All'}")
        print(f"---------------------\n")

        current_step = "Initialization" # For error tracking
        msm = None
        pcca_successful = False
        metastable_assignments = None # Maps active state index -> metastable index
        all_dtrajs_concat = None
        tica_concat = None
        clustering = None
        dtrajs = None

        try:
            # --- Step 1: Initialize Featurizer ---
            current_step = "Step 1: Initialize Featurizer"
            print(f"--- {current_step} ---")
            try:
                # Featurizer needs the topology file path
                if not os.path.exists(topology_file):
                    raise FileNotFoundError(f"Topology file '{topology_file}' not found.")
                print(f"Using topology file: '{topology_file}'")
            except FileNotFoundError as e:
                print(f"Error in {current_step}: {e}")
                exit()

            print("Creating featurizer...")
            features = pyemma.coordinates.featurizer(topology_file)
            features.add_backbone_torsions(cossin=True, periodic=False)
            print(f"Using {features.dimension()} features.")

        except Exception as e:
            print(f"An unexpected error occurred during {current_step}: {e}")
            traceback.print_exc()
            exit()

        # --- Step 2: Feature Extraction (Load trajectory and apply features) ---
        current_step = "Step 2: Feature Extraction"
        print(f"\n--- {current_step} ---")
        feat_traj = None # Initialize
        try:
            print(f"Loading trajectory '{trajectory_file}' and extracting features with stride {stride} using {n_jobs if n_jobs > 0 else 'all'} cores...")
            # Pass the featurizer object to pyemma.coordinates.load
            feat_traj = pyemma.coordinates.load(
                trajectory_file,
                features=features, # Apply features during loading
                stride=stride,
                n_jobs=n_jobs
                # chunksize=1000 # Consider adding if memory is an issue
            )
            n_frames_loaded = sum(len(chunk) for chunk in feat_traj) if isinstance(feat_traj, list) else len(feat_traj)
            print(f"Feature extraction complete. Processed {n_frames_loaded} frames.")

            if n_frames_loaded == 0:
                 raise RuntimeError("No frames were processed during feature extraction. Check trajectory file and stride.")

        except FileNotFoundError:
            print(f"Error in {current_step}: Trajectory file '{trajectory_file}' not found.")
            exit()
        except RuntimeError as e:
            print(f"Error in {current_step} (check file integrity/stride?): {e}")
            exit()
        except Exception as e:
            print(f"An unexpected error occurred during {current_step}: {e}")
            traceback.print_exc()
            exit()

        # --- Step 3: Dimensionality Reduction (TICA) ---
        current_step = "Step 3: Dimensionality Reduction (TICA)"
        print(f"\n--- {current_step} ---")
        tica_output = None # Initialize
        tica = None
        try:
            print(f"Performing TICA with lag={tica_lag} frames, target dim={tica_dim}, using {n_jobs if n_jobs > 0 else 'all'} cores...")
            tica = pyemma.coordinates.tica(
                feat_traj, # Use the feature trajectory from Step 2
                lag=tica_lag,
                dim=tica_dim,
                kinetic_map=True,
                n_jobs=n_jobs
            )
            tica_output = tica.get_output() # List of numpy arrays
            print("TICA calculation complete.")
            del feat_traj # Free memory - raw features no longer needed

            print("Visualizing TICA projection...")
            plt.figure(figsize=(10, 8))
            tica_concat = np.concatenate(tica_output) # Needed later
            plt.hexbin(tica_concat[:, 0], tica_concat[:, 1], bins='log', cmap='viridis', mincnt=1)
            plt.xlabel('TIC 1')
            plt.ylabel('TIC 2')
            plt.colorbar(label='log(Counts)')
            plt.title(f'TICA Projection (Lag={tica_lag} frames)')
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, 'tica_projection.png'))
            plt.close()
            print(f"Saved TICA projection to {os.path.join(output_dir, 'tica_projection.png')}")

            print("Visualizing TICA timescales...")
            plt.figure(figsize=(8, 6))
            tica_timescales = tica.timescales
            plt.plot(range(1, len(tica_timescales) + 1), tica_timescales, 'o-')
            plt.xlabel('Index')
            plt.ylabel(f'TICA Timescale (frames, lag={tica_lag})')
            plt.title('TICA Timescales')
            plt.grid(True, which='both', linestyle='--', linewidth=0.5)
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, 'tica_timescales.png'))
            plt.close()
            print(f"Saved TICA timescales to {os.path.join(output_dir, 'tica_timescales.png')}")

        except Exception as e:
            print(f"An unexpected error occurred during {current_step}: {e}")
            traceback.print_exc()
            exit()

        # --- Step 4: Clustering ---
        current_step = "Step 4: Clustering"
        print(f"\n--- {current_step} ---")
        try:
            print(f"Performing KMeans clustering with k={n_clusters} clusters using {n_jobs if n_jobs > 0 else 'all'} cores...")
            clustering = pyemma.coordinates.cluster_kmeans(
                tica_output, # Use the list of arrays from TICA
                k=n_clusters,
                max_iter=100,
                stride=1, # Use all TICA frames for clustering
                n_jobs=n_jobs
            )
            dtrajs = clustering.dtrajs # Get list of discrete trajectory arrays
            print("Clustering complete.")
            del tica_output # Free memory if TICA output list not needed further

            dtrajs_save_path = os.path.join(output_dir, 'dtrajs.npy')
            np.save(dtrajs_save_path, dtrajs, allow_pickle=True) # Must allow pickle for list of arrays
            print(f"Saved discretized trajectories to {dtrajs_save_path}")

            print("Visualizing cluster populations...")
            plt.figure(figsize=(12, 6))
            all_dtrajs_concat = np.concatenate(dtrajs) # Needed later
            counts = np.bincount(all_dtrajs_concat, minlength=n_clusters)
            population = counts / counts.sum()
            plt.bar(range(n_clusters), population)
            plt.xlabel('Cluster Index')
            plt.ylabel('Fractional Population')
            plt.title(f'Cluster Populations (k={n_clusters})')
            plt.grid(True, axis='y', linestyle='--', linewidth=0.5)
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, 'cluster_populations.png'))
            plt.close()
            print(f"Saved cluster populations to {os.path.join(output_dir, 'cluster_populations.png')}")

        except Exception as e:
            print(f"An unexpected error occurred during {current_step}: {e}")
            traceback.print_exc()
            exit()

        # --- Step 5: Implied Timescales Test ---
        current_step = "Step 5: Implied Timescales Test"
        print(f"\n--- {current_step} ---")
        try:
            print(f"Calculating implied timescales for lags {list(its_lags)} using {n_jobs if n_jobs > 0 else 'all'} cores...")
            its = pyemma.msm.its(
                dtrajs, # Pass the list of discretized trajectories
                lags=its_lags,
                nits=min(10, n_clusters - 1), # Request fewer timescales than states
                errors='bayes', # Use Bayesian error estimation
                n_jobs=n_jobs
            )
            print("Implied timescales calculation complete.")

            print("Plotting implied timescales...")
            plt.figure(figsize=(10, 8))
            pyemma.plots.plot_implied_timescales(its, units='frames', dt=stride) # dt=stride as lags are in strided units
            plt.xlabel(f'Lag Time ({stride} frame units)')
            plt.ylabel(f'Implied Timescale ({stride} frame units)')
            plt.title(f'Implied Timescales Analysis (k={n_clusters})')
            plt.grid(True, which='both', linestyle='--', linewidth=0.5)
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, 'implied_timescales.png'))
            plt.close()
            print(f"Saved implied timescales plot to {os.path.join(output_dir, 'implied_timescales.png')}")
            print(f"*** Action Required: Examine '{os.path.join(output_dir, 'implied_timescales.png')}' to choose an appropriate MSM lag time (msm_lag) where timescales converge. Current guess: {msm_lag} ***")

        except Exception as e:
            print(f"An error occurred during {current_step}: {e}")
            traceback.print_exc()
            exit() # Exit if ITS fails, as MSM lag choice depends on it

        # --- Step 6: Build MSM ---
        current_step = "Step 6: Build Markov State Model"
        print(f"\n--- {current_step} ---")
        try:
            print(f"Estimating MSM with lag time = {msm_lag} frames (using {stride} frame units)...")
            # Use an estimate for score_k, will be validated/adjusted later
            score_k_estimate = min(10, n_clusters - 1) if n_clusters > 1 else 1

            msm = pyemma.msm.estimate_markov_model(
                dtrajs,
                lag=msm_lag, # Lag time in units of stride
                score_method='VAMP2',
                score_k=score_k_estimate
            )
            print("MSM estimation complete.")

            # --- Print MSM Properties ---
            print(f"  Number of states (in largest connected set): {msm.nstates}")
            print(f"  Number of clusters (microstates): {n_clusters}")
            print(f"  Fraction of states used (connected states / total clusters): {msm.active_state_fraction:.3f}")

            # Validate and print VAMP score using actual number of connected states
            score_k_actual = min(10, msm.nstates - 1) if msm.nstates > 1 else 1
            if msm.nstates > 1 and score_k_actual > 0 :
                 current_score = msm.score(dtrajs, score_k=score_k_actual, score_method='VAMP2')
                 print(f"  VAMP-2 score (k={score_k_actual}): {current_score:.3f}")
            else:
                 print(f"  VAMP-2 score: N/A (connected states = {msm.nstates}, cannot score)")

            msm_timescales = msm.timescales()
            print(f"  Timescales ({stride} frame units): {[f'{t:.2f}' for t in msm_timescales]}")
            if len(msm_timescales) > 0:
                print(f"  Longest timescale: {msm_timescales[0]:.2f} ({stride} frame units)")
            else:
                print("  No timescales calculated (disconnected MSM or only 1 state?).")

        except Exception as e:
            print(f"An error occurred during {current_step}: {e}")
            traceback.print_exc()
            exit() # Exit if MSM estimation fails
      
        # --- Step 7: Identify Metastable States (PCCA++) ---
        current_step = "Step 7: Identify Metastable States (PCCA++)"
        print(f"\n--- {current_step} ---")
        pcca_successful = False # Reset flag
        metastable_assignments = None # Reset assignments
        try:
            print(f"Attempting PCCA++ clustering to identify {n_metastable} metastable states...")
            # Ensure enough connected states for PCCA++
            if msm is None:
                 raise RuntimeError("MSM object not available for PCCA++.")
            if msm.nstates >= n_metastable and msm.nstates > 1:
                msm.pcca(n_metastable)
                metastable_assignments = msm.metastable_assignments # Maps active state index (0..nstates-1) -> metastable index (0..n_metastable-1)
                print("PCCA++ complete.")
                pcca_successful = True
  
                # --- Print Metastable State Properties ---
                print("Metastable state properties:")
                for i in range(n_metastable):
                    active_state_indices_for_meta_i = np.where(metastable_assignments == i)[0]
                    if len(active_state_indices_for_meta_i) > 0:
                        prob = msm.stationary_distribution[active_state_indices_for_meta_i].sum()
                        # Map active indices back to original cluster indices for clarity
                        original_cluster_indices = msm.active_set[active_state_indices_for_meta_i]
                        print(f"  Metastable State {i}: Probability = {prob:.4f}, Contains {len(active_state_indices_for_meta_i)} active microstates (Orig. Indices e.g.: {original_cluster_indices[:5]}...)")
                    else:
                        print(f"  Metastable State {i}: Probability = 0.0000, Contains 0 active microstates")
  
                # --- Visualize Metastable States on Free Energy Landscape ---
                print("Visualizing metastable states on free energy surface...")
                plt.figure(figsize=(10, 8))
                ax_fe = plt.gca() # Get current axes
  
                print("  Calculating weights for free energy plot from stationary distribution...")
                # Create a map from original cluster index to stationary probability
                full_stationary_prob_map = np.zeros(n_clusters)
                active_set_indices = msm.active_set
                active_set_probs = msm.stationary_distribution
                for active_idx, original_cluster_idx in enumerate(active_set_indices):
                     if active_idx < len(active_set_probs): # Ensure index is valid
                          full_stationary_prob_map[original_cluster_idx] = active_set_probs[active_idx]
  
                weights = None # Default to None (uniform weights)
                try:
                    # Map the concatenated dtraj (containing original cluster indices) through the probability map
                    # Ensure all_dtrajs_concat indices are within bounds
                    if np.any(all_dtrajs_concat >= n_clusters) or np.any(all_dtrajs_concat < 0):
                         print(f"  Warning: Indices in all_dtrajs_concat out of bounds for n_clusters={n_clusters}. Using uniform weights.")
                    else:
                         weights_calc = full_stationary_prob_map[all_dtrajs_concat]
                         if np.any(np.isnan(weights_calc)) or np.all(np.abs(weights_calc) < 1e-12):
                             print("  Warning: Calculated weights are NaN or effectively zero. Using uniform weights for FE plot.")
                         elif np.sum(weights_calc) < 1e-9:
                             print(f"  Warning: Sum of calculated weights ({np.sum(weights_calc):.4f}) is very small. Using uniform weights for FE plot.")
                         else:
                             weights = weights_calc # Use calculated weights
                             print(f"  Successfully calculated weights (Sum: {np.sum(weights):.4f}).")
                except IndexError as e:
                     print(f"  Warning: IndexError calculating weights from stationary dist ({e}). Using uniform weights for FE plot.")
                except Exception as e:
                     print(f"  Warning: Error calculating weights from stationary dist ({e}). Using uniform weights for FE plot.")
  
                # Plot the free energy surface
                pyemma.plots.plot_free_energy(
                    tica_concat[:, 0], tica_concat[:, 1],
                    weights=weights, # Use calculated weights or None (uniform)
                    ax=ax_fe, cmap='viridis',
                    nbins=100, # Adjust binning as needed
                    legacy=False, minener_zero=True # Set minimum energy to zero for reference
                )
  
                print("  Calculating and plotting metastable state centers...")
                metastable_centers_tica = np.zeros((n_metastable, tica_dim))
                plotted_centers = 0
                for i in range(n_metastable):
                     active_microstate_indices_in_active_set = np.where(metastable_assignments == i)[0] # Indices within the *active set*
                     if len(active_microstate_indices_in_active_set) > 0:
                         # Map back to original cluster indices
                         original_cluster_indices = msm.active_set[active_microstate_indices_in_active_set]
  
                         if clustering is None:
                             print("  Error: Clustering object not available for calculating centers.")
                             continue
                         if np.any(original_cluster_indices >= len(clustering.clustercenters)):
                             print(f"  Error: Original cluster indices {original_cluster_indices} out of bounds for cluster centers (size {len(clustering.clustercenters)}). Skipping state {i}.")
                             continue
  
                         # Get TICA coordinates of the cluster centers for these microstates
                         centers_for_state_tica = clustering.clustercenters[original_cluster_indices]
  
                         # Get stationary probabilities for these *active* microstates
                         weights_for_state = msm.stationary_distribution[active_microstate_indices_in_active_set]
  
                         # Calculate weighted average TICA center
                         weights_sum = weights_for_state.sum()
                         if weights_sum > 1e-9: # Avoid division by zero
                             metastable_centers_tica[i] = np.average(centers_for_state_tica, axis=0, weights=weights_for_state / weights_sum)
                         else: # Fallback if weights are zero/tiny (e.g., use simple average)
                             print(f"  Warning: Near-zero weights for state {i}, using unweighted average for center.")
                             metastable_centers_tica[i] = np.average(centers_for_state_tica, axis=0)
  
                         # Plot the calculated center
                         ax_fe.scatter(
                             metastable_centers_tica[i, 0], metastable_centers_tica[i, 1],
                             s=200, marker='*', c=f'C{plotted_centers % 10}', # Cycle through default colors
                             label=f'State {i}', edgecolors='black', linewidth=1.5, zorder=3
                         )
                         plotted_centers += 1
                     # else: state i has no microstates assigned, center remains [0,0,...] and is not plotted
  
                if plotted_centers > 0:
                    ax_fe.legend(title="Metastable States")
                ax_fe.set_xlabel('TIC 1')
                ax_fe.set_ylabel('TIC 2')
                ax_fe.set_title(f'Free Energy Landscape & Metastable States (MSM lag={msm_lag}, k={n_clusters})')
                plt.tight_layout()
                plt.savefig(os.path.join(output_dir, 'metastable_states_free_energy.png'))
                plt.close()
                print(f"Saved metastable states visualization to {os.path.join(output_dir, 'metastable_states_free_energy.png')}")
  
            else:
                print(f"Skipping PCCA++: Number of connected states ({msm.nstates}) is less than requested metastable states ({n_metastable}) or is <= 1.")
  
        except RuntimeError as e:
             print(f"RuntimeError during {current_step}: {e}")
             traceback.print_exc()
        except Exception as e:
            print(f"Error during {current_step} or visualization: {e}. Check MSM connectivity and n_metastable.")
            traceback.print_exc()
            # Continue if possible, but note PCCA failed
          
 # --- Step 8: Validate MSM (Chapman-Kolmogorov Test) ---
        current_step = "Step 8: Validate MSM (Chapman-Kolmogorov Test)"
        print(f"\n--- {current_step} ---")
        # <<< Keep the conditional check you added >>>
        if pcca_successful and msm is not None and n_metastable > 0:
            try:
                print(f"Performing Chapman-Kolmogorov test for {n_metastable} states using {n_jobs if n_jobs > 0 else 'all'} cores...")
                cktest = msm.cktest(
                    n_metastable,
                    mlags=10,       # Test up to 10*msm.lag
                    err_est=True,   # Request error estimation (might not store in .confidences)
                    n_jobs=n_jobs
                )
                print("Chapman-Kolmogorov test complete.")

                # --- MODIFIED Manual Plotting (No Error Bars) ---
                print("Plotting Chapman-Kolmogorov test results (manual workaround, no error bars)...")
                # Determine grid layout
                n_rows = int(np.ceil(np.sqrt(n_metastable)))
                n_cols = int(np.ceil(n_metastable / n_rows))
                fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 4, n_rows * 3), sharex=True, sharey=True, squeeze=False)
                axes = axes.flatten() # Flatten the 2D array of axes for easy iteration

                # Check if cktest object has the necessary attributes before proceeding
                if not hasattr(cktest, 'lagtimes') or not hasattr(cktest, 'estimates') or not hasattr(cktest, 'predictions'):
                     raise AttributeError("CK test object missing required attributes (lagtimes, estimates, or predictions). Check cktest object structure.")

                lagtimes = cktest.lagtimes # Lag times used in the test

                # --- Dimension Checks (Optional but Recommended) ---
                if cktest.estimates.shape[0] != len(lagtimes) or cktest.predictions.shape[0] != len(lagtimes):
                     raise ValueError(f"Dimension mismatch: lagtimes ({len(lagtimes)}) vs estimates ({cktest.estimates.shape[0]}) or predictions ({cktest.predictions.shape[0]}).")
                if cktest.estimates.shape[1] < n_metastable or cktest.predictions.shape[1] < n_metastable:
                     # This might happen if some states are disconnected at certain lags in the test
                     print(f"Warning: CK test results have fewer states ({cktest.estimates.shape[1]}) than requested ({n_metastable}). Plotting available states.")
                     actual_n_metastable_in_ck = cktest.estimates.shape[1]
                else:
                     actual_n_metastable_in_ck = n_metastable
                # --- End Dimension Checks ---

                # Iterate up to the number of states actually present in the CK results
                for i in range(actual_n_metastable_in_ck):
                    ax = axes[i]
                    # Directly access estimates and predictions for state i
                    estimates = cktest.estimates[:, i]
                    predictions = cktest.predictions[:, i]

                    # Plot estimates as points (without error bars)
                    ax.plot(
                        lagtimes,
                        estimates,
                        marker='o',
                        linestyle='', # No connecting line for estimates
                        markersize=5,
                        label='Estimation', # Simplified label
                        color=f'C{i % 10}' # Cycle through default colors
                    )

                    # Plot predictions as lines
                    ax.plot(
                        lagtimes,
                        predictions,
                        linestyle='-',
                        linewidth=2,
                        label='Prediction',
                        color=f'C{i % 10}', # Use same color
                        alpha=0.7,
                        marker=''
                    )

                    ax.set_title(f'State {i}')
                    ax.legend()
                    ax.grid(True, alpha=0.5)

                # Add overall labels and title
                fig.suptitle(f'Chapman-Kolmogorov Test (MSM Lag: {msm_lag}, k={n_clusters}, {n_metastable} req. states)', fontsize=14)
                fig.text(0.5, 0.02, f'Lag Time ({stride} frame units)', ha='center', va='center', fontsize=12) # Use stride units
                fig.text(0.02, 0.5, 'Conditional Probability P(lag)', ha='center', va='center', rotation='vertical', fontsize=12)

                # Hide unused axes if the grid is larger than needed
                for j in range(actual_n_metastable_in_ck, n_rows * n_cols):
                    if j < len(axes): # Ensure index is valid
                         fig.delaxes(axes[j])

                plt.tight_layout(rect=[0.03, 0.03, 0.97, 0.95]) # Adjust layout
                ck_plot_path = os.path.join(output_dir, 'chapman_kolmogorov_test_manual.png')
                plt.savefig(ck_plot_path)
                plt.close(fig)
                print(f"Saved manually generated Chapman-Kolmogorov test plot (no error bars) to {ck_plot_path}")
                print("*** Action Required: Examine CK test plot. Predictions (lines) should match estimations (dots) for a valid MSM. ***")

            except AttributeError as e:
                 print(f"\nAttributeError during CK plotting: {e}. The cktest object structure might have changed or lacks expected data.")
                 print("Traceback:")
                 traceback.print_exc()
                 print("\nSkipping CK test visualization due to error.")
            except ValueError as e:
                 print(f"\nValueError during CK plotting: {e}. Check dimensions of cktest results.")
                 print("Traceback:")
                 traceback.print_exc()
                 print("\nSkipping CK test visualization due to error.")
            except Exception as e:
                # Catch other potential errors during plotting
                print(f"\nError during Chapman-Kolmogorov test plotting: {e}")
                print("Traceback:")
                traceback.print_exc()
                print("\nSkipping CK test visualization due to error.")
        # <<< Keep the else block you added >>>
        else:
            # Provide reasons for skipping
            skip_reason = []
            if not pcca_successful:
                skip_reason.append("PCCA++ did not complete successfully or was skipped")
            if msm is None:
                skip_reason.append("MSM object is missing")
            if n_metastable <= 0:
                 skip_reason.append(f"Number of metastable states requested ({n_metastable}) is not positive")
            print(f"Skipping Chapman-Kolmogorov test because: {', '.join(skip_reason)}.")
 
        # --- Step 9: Extract Representative Structures ---
        current_step = "Step 9: Extract Representative Structures"
        print(f"\n--- {current_step} ---")
        # Check PCCA success and that assignments were actually made
        if pcca_successful and metastable_assignments is not None and msm is not None:
            print(f"Extracting up to {n_samples_per_state} representative structures per metastable state...")
            traj_full = None # Initialize traj_full outside try block for potential cleanup
            try:
                # Map microstate assignments (from clustering, stored in all_dtrajs_concat)
                # to metastable assignments (from PCCA). This requires mapping original cluster indices
                # to active set indices, then using metastable_assignments.

                # Create a mapping: original cluster index -> metastable index (or -1 if not active)
                cluster_to_metastable_map = -1 * np.ones(n_clusters, dtype=int)
                for active_idx, original_cluster_idx in enumerate(msm.active_set):
                    if active_idx < len(metastable_assignments): # Ensure index is valid
                        cluster_to_metastable_map[original_cluster_idx] = metastable_assignments[active_idx]

                # Apply this map to the concatenated discrete trajectory
                metastable_traj_assignments = cluster_to_metastable_map[all_dtrajs_concat]

                output_pdb_dir = os.path.join(output_dir, 'representative_structures')
                os.makedirs(output_pdb_dir, exist_ok=True)

                # Load trajectory once for efficient slicing using MDTraj directly
                print("Loading full trajectory for structure extraction (this might take time)...")
                # FIX: Use mdtraj.load to get an MDTraj object
                traj_full = md.load(trajectory_file, top=topology_file, stride=1)

                # Check if loading was successful and it's an MDTraj object before proceeding
                if not isinstance(traj_full, md.Trajectory) or traj_full.n_frames == 0:
                    raise RuntimeError(f"Failed to load trajectory using mdtraj or trajectory is empty. Check file: {trajectory_file}")

                print(f"Full trajectory loaded with {traj_full.n_frames} frames.") # This line should now work

                for i in range(n_metastable):
                    # Find indices in the *concatenated dtraj* that belong to metastable state i
                    # These indices correspond to frames in the *strided* trajectory used for analysis
                    state_frames_in_strided_traj = np.where(metastable_traj_assignments == i)[0]

                    if len(state_frames_in_strided_traj) == 0:
                        print(f"  Metastable State {i}: No frames assigned. Skipping structure extraction.")
                        continue

                    # Map these indices back to the original trajectory frame numbers (accounting for analysis stride)
                    # original_frame_index = strided_frame_index * stride + initial_offset (assume offset=0)
                    original_unstrided_indices = state_frames_in_strided_traj * stride

                    # Select a subset of these original frame indices
                    n_available = len(original_unstrided_indices)
                    n_to_sample = min(n_samples_per_state, n_available)
                    # Use linspace to get evenly spaced indices, convert to int
                    sample_indices_in_original_traj = original_unstrided_indices[np.linspace(0, n_available - 1, n_to_sample, dtype=int)]

                    # Ensure indices are within the bounds of the fully loaded trajectory
                    valid_indices = sample_indices_in_original_traj[sample_indices_in_original_traj < traj_full.n_frames]

                    if len(valid_indices) == 0:
                         print(f"  Warning: Metastable state {i} - No valid frame indices found after mapping/bounds check (Original indices: {sample_indices_in_original_traj[:5]}...). Skipping.")
                         continue

                    # --- ADDED: Print representative frame indices ---
                    # These indices are 0-based and refer to the original, unstrided trajectory file
                    print(f"  Metastable State {i}: Representative frame indices (0-based): {valid_indices.tolist()}")
                    # --- End ADDED section ---

                    outfile = os.path.join(output_pdb_dir, f'metastable_state_{i}_samples.pdb')
                    try:
                        # Slice the pre-loaded full trajectory using the valid indices
                        sampled_structures = traj_full.slice(valid_indices)

                        if sampled_structures and sampled_structures.n_frames > 0:
                            # Save as a multi-frame PDB
                            sampled_structures.save_pdb(outfile)
                            print(f"  Saved {sampled_structures.n_frames} representative structures for metastable state {i} to {outfile}")
                        else:
                             print(f"  Warning: Could not extract frames for metastable state {i} (Indices: {valid_indices}). Slice returned empty.")

                    except IndexError as e:
                         print(f"  Error slicing frames for state {i} (IndexError: {e}).")
                         print(f"  Attempted indices: {valid_indices[:10]}... Max frame: {traj_full.n_frames - 1}")
                    except AttributeError as ae:
                         # Catch if traj_full is somehow not an MDTraj object
                         print(f"  Error: 'traj_full' object does not have expected slicing method. Type: {type(traj_full)}. Error: {ae}")
                         traceback.print_exc()
                         break # Stop extraction if traj object is broken
                    except Exception as e:
                        print(f"  Error processing or saving structures for state {i}: {e}")
                        traceback.print_exc()

            except FileNotFoundError:
                print(f"Error in {current_step}: Trajectory file '{trajectory_file}' or topology '{topology_file}' not found for structure extraction.")
            except RuntimeError as e:
                 print(f"RuntimeError during {current_step}: {e}")
                 traceback.print_exc()
            except Exception as e:
                print(f"An unexpected error occurred during {current_step} setup or execution: {e}")
                traceback.print_exc()
            finally:
                 # Ensure trajectory object is deleted to free memory, even if errors occurred
                 if traj_full is not None:
                     del traj_full
                     print("Cleaned up full trajectory object from memory.")

        else:
            print("Skipping structure extraction because PCCA++ did not complete successfully, was skipped, or assignments/MSM are missing.")

 
        # --- Step 10: Save MSM Object ---
        current_step = "Step 10: Save MSM Object"
        print(f"\n--- {current_step} ---")
        if msm is not None:
            try:
                msm_save_path = os.path.join(output_dir, 'msm_model.pyemma')
                msm.save(msm_save_path, overwrite=True)
                print(f"Saved MSM object to {msm_save_path}")
            except Exception as e:
                print(f"An error occurred during {current_step}: {e}")
                traceback.print_exc()
        else:
            print("Skipping MSM save because the MSM object was not successfully created.")
  
  
        # --- Final Summary ---
        end_time = time.time()
        print("\n--- MSM Analysis Workflow Complete ---")
        print(f"Total execution time: {end_time - start_time:.2f} seconds")
        print(f"All outputs saved in: {output_dir}")
        print(f"Log file saved to: {log_file_path}")
        print(f"Run summary: k={n_clusters}, tica_lag={tica_lag}, msm_lag={msm_lag}, n_metastable={n_metastable}")
        if msm:
           print(f"MSM details: {msm.nstates} connected states from {n_clusters} clusters.")
        if pcca_successful:
            print("PCCA++ completed successfully.")
        else:
            print("PCCA++ was skipped or did not complete successfully.")
        print("--- End of Log ---")
  
  # --- End of Script ---
