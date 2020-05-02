function [eigenVectors eigenValues weights] = V_Matrix_SVD_res(data,eval_max)
%% DESCRIPTION:
%
% [eigenVectors eigenValues weights] = V_Matrix_SVD_res(data,eval_max)
%
%   Computes the eigenvectors and eigenvalues using a Singular Value
%       Decomposition (SVD) method on the data.
%       Compared with V_Matrix_SVD, this function only keep the rows of
%       weights according the number of dimensions.
%
%   Author: VVT (edited by Martin Ullrich and Benjamin Dourthe)
%
%   Created: October 10th, 2006
%   Last update: April 7th, 2018
%
%% Input: 
%   - data: m x n matrix containing data (m = trials/subjects/obeservations)
%   - eval_max: if set to 0 then the maximal number of eigenvectors will be
%       computed
%
%% Output:
%   - eigenVectors: arranged as rows
%   - eigenValues: arranged as a column vector
%   - weights: projections of the trials onto the eigenvectors (each column
%       corresponds to one trials/subjects/obeservations)
%
%% Dependencies:
%   Functions:
%       - None
%   Files:
%       - None
  
%% Function

    % Normalize data around the mean
    M_mean = mean(data,1);
    M = data - M_mean;

    % Initialization
    [r c] = size(M);
    trials = r;
    
    % Singular Value Decomposition (SVD) for computation of eigenvalues and
    % eigenvectors
    if r >= c
        C = M' * M / trials;
        [U,eigenValues,eigenVectors] = svd(C,'econ');
        clear U

        if ~(eval_max == 0)
            eigenValues = diag(eigenValues);
            if eval_max > length(eigenValues)
                eval_max = length(eigenValues);
                disp('M_Matrix_SVD_res requested number of eigenvectors was too large');
            end
            eigenValues = eigenValues(1:eval_max);
            eigenVectors = eigenVectors(:,1:eval_max)';
        else
            eigenValues = diag(eigenValues);
            eigenVectors = eigenVectors';
        end    

    else

        C = M * M' / trials;
        [U,eigenValues,eigenVectors] = svd(C);
        clear U

        if ~(eval_max == 0)
            eigenValues = diag(eigenValues);
            eigenValues = eigenValues(1:eval_max);
            if eval_max > length(eigenValues)
                eval_max = length(eigenValues);
                disp('V_Matrix_SVD_res requested number of eigenvectors was too large');
            end
            eigenVectors = eigenVectors(:,1:eval_max);
            neg_sqrt_S = diag(-(sqrt(eigenValues)).^-1);
            eigenVectors = M'*eigenVectors*neg_sqrt_S / sqrt(trials);
            eigenVectors = eigenVectors';
        else
            eigenValues = diag(eigenValues);
            neg_sqrt_S = diag(-(sqrt(eigenValues)).^-1);
            eigenVectors = M'*eigenVectors*neg_sqrt_S / sqrt(trials);
            eigenVectors = eigenVectors';
        end
    end


    % Project mean on eigenvectors -> weights that describe mean
    mean_M_proj = eigenVectors * M_mean';

    % Reconstruct mean from reduced number of eigenvectors
    recon_M_mean = mean_M_proj' * eigenVectors;

    % Calculation of the residual between real mean and reconstruction
    % + normalization
    M_residual_mean = M_mean - recon_M_mean;
    M_residual_mean = M_residual_mean / norm(M_residual_mean);

    eigenVectors = [M_residual_mean;eigenVectors];

    weights = eigenVectors * data';
    % Each row represents the weights for the trials on the corresponding eigenvector

    % The reconstruction is dataReconstructed = W' * V

%% Change all negative weights and vectors to positve

    % meanW = mean(W,2);
    % for j = 1 : length(meanW)
    %     if meanW(j) < 0
    %         V(j,:) = -V(j,:);
    %         W(j,:) = -W(j,:);
    %     end
    % end
    
