function [eigenVectors eigenValues weights] = V_Matrix_SVD(data,eval_max)
%% DESCRIPTION:
%
% [eigenVectors eigenValues weights] = V_Matrix_SVD(data,eval_max)
%
%   Computes the eigenvectors and eigenvalues using a Singular Value
%       Decomposition (SVD) method on the data.
%
%   Author: VVT (edited by Benjamin Dourthe)
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

    % Initialization
    [r c] = size(data);
    trials = r;

    % Singular Value Decomposition (SVD) for computation of eigenvalues and
    % eigenvectors
    if r >= c
        C = data' * data / trials;
        [U,eigenValues,eigenVectors] = svd(C,'econ');
        clear U

        if ~(eval_max == 0)
            eigenValues = diag(eigenValues);
            if eval_max > length(eigenValues)
                eval_max = length(eigenValues);
                disp('W_Matrix_SVD requested number of eigenvectors was too large');
            end
            eigenValues = eigenValues(1:eval_max);
            eigenVectors = eigenVectors(:,1:eval_max)';
        else
            eigenValues = diag(eigenValues);
            eigenVectors = eigenVectors';
        end

    else

        C = data * data' / trials;
        [U,eigenValues,eigenVectors] = svd(C);
        clear U
        if ~(eval_max == 0)
            eigenValues = diag(eigenValues);
            eigenValues = eigenValues(1:eval_max);
            if eval_max > length(eigenValues)
                eval_max = length(eigenValues);
                disp('W_Matrix_SVD requested number of eigenvectors was too large');
            end
            eigenVectors = eigenVectors(:,1:eval_max);
            neg_sqrt_S = diag(-(sqrt(eigenValues)).^-1);
            eigenVectors = data'*eigenVectors*neg_sqrt_S / sqrt(trials);
            eigenVectors = eigenVectors';
        else
            eigenValues = diag(eigenValues);
            neg_sqrt_S = diag(-(sqrt(eigenValues)).^-1);
            eigenVectors = data'*eigenVectors*neg_sqrt_S / sqrt(trials);
            eigenVectors = eigenVectors';
        end

    end

    weights = eigenVectors * data';
    % Each row represents the weights for the trials on the corresponding eigenvector

    % The reconstruction is dataReconstructed = W' * V
    
