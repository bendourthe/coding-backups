# ---------------------------------------------------------------------------- #
#                             CODE DESCRIPTION                                 #
# ---------------------------------------------------------------------------- #
'''
[eigen_vectors eigen_values weights] = V_Matrix_SVD_res(data,eval_max)

   Computes the eigenvectors and eigenvalues using a Singular Value
   Decomposition (SVD) method on the data.
   Compared with V_Matrix_SVD, this function only keep the rows of
   weights according the number of dimensions.

Author: VVT (edited by Martin Ullrich and Benjamin Dourthe)

Created: October 10th, 2006
Last update: April 7th, 2018

Input: 
   - data: m x n matrix containing data (m = trials/subjects/obeservations)
   - eval_max: if set to 0 then the maximal number of eigenvectors will be
     computed

Output:
  - eigenVectors: arranged as rows
  - eigenValues: arranged as a column vector
  - weights: projections of the trials onto the eigenvectors (each column
    corresponds to one trials/subjects/obeservations)

Dependencies:
  Functions:
      - None
   Files:
       - None
'''  

# ---------------------------------------------------------------------------- #
#                                 FUNCTION                                     #
# ---------------------------------------------------------------------------- #

def SVD_residual(data,num_eigenvectors):
    
    # Import libraries
    import numpy as np

    # Normalize data around the mean
    norm_data = data - np.mean(data,0)

    # Initialization
    [r,c] = np.shape(norm_data)
    num_observations = r;
    
    # Singular Value Decomposition (SVD):
    #     computation of eigenvalues and eigenvectors
    if r >= c:
        sim_data = zip(*norm_data) * norm_data / num_observations
        U,eigen_values,eigen_vectors = np.linalg.svd(sim_data,full_matrices=False)
        del U

        if num_eigenvectors in ['max']:
            eigen_values = np.diagonal(eigen_values)
            eigen_vectors = zip(*eigen_vectors)
            
        else:
            eigen_values = np.diagonal(eigen_values)
            if num_eigenvectors > len(eigen_values):
                num_eigenvectors = np.len(eigen_values)
                print('The requested number of eigenvectors is too large')
            eigen_values = eigen_values[1:num_eigenvectors]
            eigen_vectors = zip(*eigen_vectors[:,1:num_eigenvectors])            

    else:
        sim_data = norm_data * zip(*norm_data) / num_observations
        U,eigen_values,eigen_vectors = np.linalg.svd(sim_data)
        del U

        if num_eigenvectors in ['max']:
            eigen_values = np.diagonal(eigen_values)
            eigen_vectors = zip(*eigen_vectors)
            neg_sqrt_eigen_values = np.diagonal(-(np.sqrt(eigen_values))**-1)
            eigen_vectors = zip(*norm_data)*eigen_vectors*neg_sqrt_eigen_values / np.sqrt(num_observations)
            eigen_vectors = zip(*eigen_vectors)
        else:
            eigen_values = np.diagonal(eigen_values)
            eigen_values = eigen_values[1:num_eigenvectors]
            if num_eigenvectors > np.len(eigen_values):
                num_eigenvectors = np.len(eigen_values)
                print('The requested number of eigenvectors was too large')
            eigen_vectors = eigen_vectors[:,1:num_eigenvectors]
            neg_sqrt_eigen_values = np.diagonal(-(np.sqrt(eigen_values))**-1)
            eigen_vectors = zip(*norm_data)*eigen_vectors*neg_sqrt_eigen_values / np.sqrt(num_observations)
            eigen_vectors = zip(*eigen_vectors)

    # Project mean on eigenvectors to obtain the weights that describe the mean
    mean_data_proj = eigen_vectors * zip(*norm_data)

    # Reconstruct mean from reduced number of eigenvectors
    recon_data_mean = zip(*mean_data_proj) * eigen_vectors

    # Calculation of the residual between real mean and reconstruction
    data_residual_mean = norm_data - recon_data_mean
    
    # Normalize the residual
    data_residual_mean = data_residual_mean / np.linalg.norm(data_residual_mean)

    # Combine residual with original eigenvectors
    eigen_vectors = [data_residual_mean,eigen_vectors]
    
    # Calculate the weights
    #   Note: Each row represents the weights of the observations of the corresponding eigenvector
    weights = eigen_vectors *zip(*data)

    # Note: To reconstruct the original data: data_reconstructed = zip(*weights) * eigen_vectors

    
    # Change all negative weights and eigenvectors to positive

    mean_weights = np.mean(weights,1)
    for i in range(1 , len(mean_weights)):
        if mean_weights(i) < 0:
            eigen_vectors[i,:] = -eigen_vectors[i,:]
            weights[i,:] = -weights[i,:]

    
