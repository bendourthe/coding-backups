function [theta, cost] = logic_reg(trainingSet,trainingState)
%% DESCRIPTION:
%
% [theta, cost] = logic_reg(trainingSet,trainingState)
%
%   Performs a logistic regression algorithm on a training set to extract
%   the corresponding optimized Theta and Cost
%
%   Author: Benjamin Dourthe
%
%   Created: April 10th, 2018
%   Last update: April 10th, 2018
%
%% Input: 
%   - trainingSet: m x n matrix containing data to feed the logistic
%   regression algorithm (m = trials/subjects/obeservations)
%   - trainingState: m x 1 matrix containing the states of each
%   trials/subjects/obeservations
%
%% Output:
%   - theta: optimal theta obtained after runing the fminunc function
%   - cost: optimal cost obtained after runing the fminunc function
%
%% Dependencies:
%   Functions:
%       - sigmoid.m
%       - costFunctionLG.m
%   Files:
%       - None
  
%% Function

% Compute Cost and Gradient for Logistic Regression (classification)

    %  Setup the data matrix appropriately, and add ones for the intercept term
    [m, n] = size(trainingSet);

    % Add intercept term to x and X_test
    trainingSet = [ones(m, 1) trainingSet];

    % Initialize fitting parameters
    initial_theta = zeros(n + 1, 1);

    % Compute and display initial cost and gradient
    [cost, grad] = costFunctionLG(initial_theta, trainingSet, trainingState);
    
% Optimization of Theta
    
    % Set options for fminunc
    options = optimoptions('fminunc','Algorithm','quasi-newton');

    % Run fminunc to obtain the optimal theta
    % This function will return theta and the cost 
    [theta, cost] = fminunc(@(t)(costFunctionLG(t, trainingSet, trainingState)),...
        initial_theta, options);
