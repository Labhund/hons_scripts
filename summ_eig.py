### Add Eigenvalues from gmx covar

import numpy as np
import os
import sys

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
            if line.strip() and not line.startswith(('#','@')): # Skip empty lines and comments
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

def percent_reached(eigenvalues, n_modes=0, percent=0.8):
    """
    Calculate the number of modes required to reach a given percentage of the total eigenvalues,
    or calculate the percentage reached by a given number of modes.
    
    Parameters:
    eigenvalues (np.ndarray): Array of eigenvalues.
    n_modes (int): Number of modes to consider. If 0, it will calculate the number of modes 
                  needed to reach the percentage.
    percent (float): Percentage of the total eigenvalues to reach (default is 0.8).
    
    Returns:
    int or float: If n_modes=0, returns the number of modes required to reach the given percentage.
                 If n_modes>0, returns the percentage of total eigenvalues reached by n_modes.
    """
    total_eigenvalues = summ_eig(eigenvalues, len(eigenvalues))
    
    if n_modes == 0:
        # Find the number of modes that reach the target percentage
        target_eigenvalues = total_eigenvalues * percent
        cumulative_eigenvalues = np.cumsum(eigenvalues)
        modes = np.where(cumulative_eigenvalues >= target_eigenvalues)[0][0] + 1
        return modes
    else:
        # Calculate the percentage reached by n_modes
        if n_modes > len(eigenvalues):
            raise ValueError("n_modes exceeds the number of eigenvalues.")
        sum_n_modes = summ_eig(eigenvalues, n_modes)
        return sum_n_modes / total_eigenvalues

def main():
    # Check if a filename is provided as a command line argument
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        filename = 'eigenvalues.txt'  # Default filename if none provided
    
    try:
        # Read eigenvalues from the file
        eigenvalues = read_eigenvalues(filename)
        
        # Calculate the sum of all eigenvalues
        total_eigenvalues = summ_eig(eigenvalues, len(eigenvalues))
        
        # Calculate the number of modes required to reach 80% of the total eigenvalues
        modes_80_percent = percent_reached(eigenvalues, percent=0.8)

        # Print the results
        print(f"Input file: {filename}")
        print(f"Total number of eigenvalues: {len(eigenvalues)}")
        print(f"Sum of all eigenvalues: {total_eigenvalues}")
        
        # Percent reached with specific n_modes values
        print("-" * 40)
        print("Percent reached for different n_modes:")
        print("n_modes\t% reached\tsum_eigenvalues")
        n_modes_list = [2, 3, 5, 10, 15, 20]  # Example n_modes values
        
        for n in n_modes_list:
            if n <= len(eigenvalues):
                percent_val = percent_reached(eigenvalues, n_modes=n) * 100
                sum_val = summ_eig(eigenvalues, n)
                print(f"{n}\t{percent_val:.2f}%\t\t{sum_val:.6f}")
        
        print("-" * 40)
        print(f"Number of modes required to reach 80% of the total eigenvalues: {modes_80_percent}")
    
    except IOError as e:
        print(f"Error: {e}")
        sys.exit(1)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()