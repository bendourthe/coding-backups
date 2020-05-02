function [dataSetPCA] = PCA_red(dataSet,dataState,variance)
%% DESCRIPTION:
%
% [dataSetPCA] = PCA_red(dataSet,dataState)
%
%   Performs a a Principal Component Analysis and a data set a generate the
%   corresponding project using the specific amount of variance explained
%
%   Author: Benjamin Dourthe
%
%   Created: April 10th, 2018
%   Last update: April 10th, 2018
%
%% Input: 
%   - dataSet: m x n matrix containing data to apply PCA on
%   (m = trials/subjects/obeservations)
%   - dataState: m x 1 matrix containing the states of each
%   trials/subjects/obeservations
%   - variance: amount of variance that should be explained by the
%   Principal Components (will determine how many PCs will be kept in the
%   projected dataSet)
%
%% Output:
%   - dataSetPCA: projection of the dataSet in PC space (reduced
%   dimensionality)
%
%% Dependencies:
%   Functions:
%       - PCA.m
%       - V_Matrix_SVD.m
%       - V_Matrix_SVD_res.m
%   Files:
%       - None
  
%% Function

    % Apply PCA on Data Set
    [whitenedWeights,PCAvars] = PCAnew(dataSet, 0, dataState);

    % Convert the eigenvalues into percentage to explain the contribution
    % of each PC into the total variance
    explained = PCAvars.eigenvalues*100/sum(PCAvars.eigenvalues);
    
    % Generate a matrix of added contribution (useful to define threshold,
    % i.e. how much of the variance do we want to explain with our
    % dimensionality reduction)
    for i=1:length(explained)
        sumExplained(i) = sum(explained(1:i));
    end
    
    % Find how many PC are needed to reach the corresponding threshold
    idx = find(sumExplained >= variance);
    
    % Generate a new training set based on the PCA (lower dimension)
    dataSetPCA = PCAvars.weights(2:idx(1)+1,:)';