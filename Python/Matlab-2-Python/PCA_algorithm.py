# ---------------------------------------------------------------------------- #
#                             CODE DESCRIPTION                                 #
# ---------------------------------------------------------------------------- #

'''
PCA_dim_reduction(data , variance)

   Computes the eigenvalues, eigenvectors and weights of a data set using 
   Singular Value Decomposition (SVD), then project the original data set into
   a Principal Components (PCs) space. Based on the defined amount of variance
   explained, the projected data will be cropped by only keeping the PCs that
   explain the desired amount of variance.

Author: Benjamin Dourthe

Created: September 15th, 2018
Last update: September 27th, 2018

Input: 
    - data: m x n matrix containing data, where m = trials/subjects/observations
    and n = number of variables
    - variance: percentage of variance desired (between 1 and 100)

Output:
    - eigen_values: arranged as a column vector
    - eigen_vectors: arranged as rows
    - weights: projections of the trials onto the eigenvectors (each column
    corresponds to one trial/subject/observation)
    - variance_explained_acc: vector showing the accumulative contribution of 
    each PC to the total variance (last row should be 100)
    - projected_data: original data set projected into PC space (only includes
    the PCs explaining the desired amount of variance)

Dependencies:
    Libraries:
        - numpy as np
    Functions:
        - None
    Files:
        - None
'''  

# ---------------------------------------------------------------------------- #
#                                 FUNCTION                                     #
# ---------------------------------------------------------------------------- #

# LIBRARIES
import numpy as np

def PCA_algorithm(data , variance):

# STEP 1 - Singular Value Decomposition (SVD) -------------------------------- #

    # Centralize data around zero by subtracting the mean 
    # (leads to a data set with zero mean)
    zero_mean_data = data - np.mean(data , 0)

    # Calculate number of observations and variables
    [r , c] = np.shape(zero_mean_data)
    num_observations = r
    num_variables = c    
    
    # SVD: computes eigenvalues and eigenvectors
    if r >= c:
        sim_data = zip(*zero_mean_data) * zero_mean_data / num_observations
        U , eigen_values , eigen_vectors = np.linalg.svd(sim_data , full_matrices = False)
        del U
        eigen_values = np.diagonal(eigen_values)
        eigen_vectors = zip(*eigen_vectors)
    else:
        sim_data = zero_mean_data * zip(*zero_mean_data) / num_observations
        U , eigen_values , eigen_vectors = np.linalg.svd(sim_data)
        del U
        eigen_values = np.diagonal(eigen_values)
        eigen_vectors = zip(*eigen_vectors)
        neg_sqrt_eigen_values = np.diagonal(-(np.sqrt(eigen_values))**-1)
        eigen_vectors = zip(*zero_mean_data) * eigen_vectors * neg_sqrt_eigen_values / np.sqrt(num_observations)
        eigen_vectors = zip(*eigen_vectors)

# STEP 2 - Residual Calculation ---------------------------------------------- #

    # Project mean on eigenvectors to obtain the weights that describe the mean
    mean_data_proj = eigen_vectors * zip(*zero_mean_data)

    # Reconstruct mean from reduced number of eigenvectors
    recon_data_mean = zip(*mean_data_proj) * eigen_vectors

    # Calculation of the residual between real mean and reconstruction
    data_residual_mean = zero_mean_data - recon_data_mean
    
    # Normalize the residual
    data_residual_mean = data_residual_mean / np.linalg.norm(data_residual_mean)

    # Combine residual with original eigenvectors
    eigen_vectors = [data_residual_mean , eigen_vectors]

# STEP 3 - Weights -Calculation ---------------------------------------------- #
        
    # Calculate the weights
    #   Note: Each row represents the weights of the observations of the corresponding eigenvector
    weights = eigen_vectors * zip(*data)

    # Note: To reconstruct the original data: data_reconstructed = zip(*weights) * eigen_vectors
    
# STEP 4 - Positivity -------------------------------------------------------- #
    
    # Change all negative weights and eigenvectors to positive

    mean_weights = np.mean(weights , 1)
    for i in range(1 , np.len(mean_weights)):
        if mean_weights(i) < 0:
            eigen_vectors[i , :] = -eigen_vectors[i , :]
            weights[i , :] = -weights[i , :]

# STEP 5 - Variance ---------------------------------------------------------- #            
                                    
    # Convert the eigenvalues into percentages to explain the contribution
    # of each PC into the total variance
    explained = eigen_values * 100 / np.sum(eigen_values)
    
    # Generate a matrix of accumulative contributions: used to define how many
    # PCs will be kept to explain the desired amount of variance
    variance_explained_acc = []
    for i in range(1 , np.len(explained)):
        variance_explained_acc[i] = np.sum(explained[1:i])
    
    # Identify how many PCs are needed to reach the corresponding threshold
    for i in range(1, np.len(variance_explained_acc)):
        if variance_explained_acc[i] >= variance:
            num_PCs = i
    
    # Generate a new data set projected into the PC space and cropped based on
    # the amount of PCs required to explain the defined amount of variance
    projected_data = zip(*weights[2:num_PCs(1)+1 , :])

# STEP 6 - RETURN OUTPUTS ---------------------------------------------------- #

    return eigen_values
    return eigen_vectors
    return weights
    return variance_explained_acc
    return projected_data
