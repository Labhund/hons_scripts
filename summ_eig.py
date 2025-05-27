### Add Eigenvalues from gmx covar

import numpy as np
import os

def read_eigenvalues(filename):
    """
    Read eigenvalues from a file.

    Parameters:
    filename (str): Path to the file containing eigenvalues.

    Returns:
    np.ndarray: Array of eigenvalues.
    """
    try:
        with open(filename, 'r') as file:
            lines = file.readlines()
        
        # Extract eigenvalues from the lines
        eigenvalues = []
        for line in lines:
            if line.strip() and not line.startswith('#','@'): # Skip empty lines and comments
                parts = line.split()
                if len(parts) > 1 and parts[0].isdigit():
                    eigenvalues.append(float(parts[1]))
        
        return np.array(eigenvalues)
    
    except Exception as e:
        raise IOError(f"Error reading file {filename}: {e}")

def summ_eig(eigenvalues, n_modes):
    """
    Sum the eigenvalues from a covariance matrix.

    Parameters:
    eigenvalues (np.ndarray): Array of eigenvalues.
    n_modes (int): Number of modes to consider.

    Returns:
    float: Sum of the first n_modes eigenvalues.
    """
    if len(eigenvalues) < n_modes:
        raise ValueError("Number of modes exceeds the number of eigenvalues.")
    
    return np.sum(eigenvalues[:n_modes])

def eighty_percent(eigenvalues):
    """
    Calculate the number of modes required to reach 80% of the total eigenvalues.
    Parameters:
    eigenvalues (np.ndarray): Array of eigenvalues.
    Returns:
    int: Number of modes required to reach 80% of the total eigenvalues.
    """
    total_eigenvalues = summ_eig(eigenvalues, len(eigenvalues))
    target_eigenvalues = 0.8 * total_eigenvalues
    cumulative_eigenvalues = np.cumsum(eigenvalues)
    modes = np.where(cumulative_eigenvalues >= target_eigenvalues)[0][0] + 1
    return modes

def main():
    # Read eigenvalues from the file
    eigenvalues = read_eigenvalues('eigenvalues.txt')
    # Calculate the sum of the first n_modes eigenvalues
    n_modes = len(eigenvalues)  # You can set this to a specific number if needed
    sum_eig = summ_eig(eigenvalues, n_modes)
    # Calculate the number of modes required to reach 80% of the total eigenvalues
    modes_80_percent = eighty_percent(eigenvalues)
    # Print the results
    print(f"Sum of the first {n_modes} eigenvalues: {sum_eig}")
    print(f"Number of modes required to reach 80% of the total eigenvalues: {modes_80_percent}")

if __name__ == "__main__":
    main()