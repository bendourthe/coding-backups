function [whiteWeights,PCAvars] = PCAnew(data, variance, state)
%% DESCRIPTION:
%
% [whiteWeights,PCAvars] = PCAnew(data, variance, state)
%
%   Performs two principle component analyses:
%       (1) used to determine the number of PCs that explain the specified
%       ammount of variance
%       (2) computes the actual Principle Components (PCs) and the residual
%
%   The resulting weights are withened by dividing them by the standard
%   deviation of the first group.
%   Returned are the whitened weights and a structure containing the
%   original PCA results. w(1) is always the residual.
%
%   Author: Fabian Hoitz (edited by Benjamin Dourthe)
%
%   Created: November 1st, 2017
%   Last update: April 7th, 2018
%
%% Input: 
%   - data: m x n matrix containing data for PCA (m = trials/subjects/obeservations)
%   - variance: amount of variance that should be explained by the PCs
%       (value between 1 and 100 - if 0 all PCs will be computed)
%   - state: m x 1 vector describing what state each row relates to (should
%       have same amount of rows as data)
%
%% Output:
%   - withWeights: weights after whitening process (usefull for SVM)
%   - PCAvars: structure with original PCA weights, eigenvalues and
%       eigenvectors
%
%% Dependencies:
%   Functions:
%       - V_Matrix_SVD.m
%       - V_Matrix_SVD_res.m
%   Files:
%       - None
  
%% Function

    % Find number of PC's that explain specified variance
    PCAvars.originalMean = mean(data);

    % Subtract mean from original input
    data = data - mean(data);

    % Run PCA without residuals to determine x% of variance
    [~,s,~] = V_Matrix_SVD(data, 0);

    % Find x% of variance
    if variance == 0
        indx_var = 0;
    elseif variance >= 1 && variance <= 100
        indx_var = sum(cumsum(s/sum(s)) <= variance/100);
        % if ammount of variance explained is smaller than 1st pc only
        % calculate pc 1
        if indx_var == 0
            indx_var = 1;
        end
    end

    % run PCA with residuals and selected ammount of varianze
    [v,s,w] = V_Matrix_SVD_res(data, indx_var);

    % subtract mean of principle components from weights
    w_mean_subtracted = w - mean(w,2);

    % Whitening
        % Take std of pc's of group 1 for whitening
        groupLabels = unique(state);
        std_w_grp1 = std(w_mean_subtracted(:,state == groupLabels(1)),0,2);

        % Divide weights by std of first group
        w_whitened = w_mean_subtracted ./ std_w_grp1;

    % Organize output variables
    PCAvars.weights = w;
    PCAvars.eigenvectors = v;
    PCAvars.eigenvalues = s;

    whiteWeights = w_whitened;
end

