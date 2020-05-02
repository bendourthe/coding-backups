function p = predictLG(theta, X)
%% DESCRIPTION:
%
% p = predict(theta, X)
%
% Predict whether the label is 0 or 1 using learned logistic regression
%   parameters theta. Computes the predictions for X using a threshold at
%   0.5 (i.e., if sigmoid(theta'*x) >= 0.5, predict 1)
%
%   Author: Benjamin Dourthe
%
%   Created: October 10th, 2016
%   Last update: April 8th, 2016
%
%% Input: 
%   - signal: original sinusoidal signal (size: m x 1)
%   - samplingRate: frequency of the signal in Hz
%
%% Output:
%   - p: predicted state
%
%% Dependencies:
%   Functions:
%       - sigmoid.m
%   Files:
%       - None

%% Function

% Number of training examples

    m = size(X, 1); 

% Prediction

    p = sigmoid(theta' * X')';
    p(p >= 0.5) = 1;
    p(p < 0.5) = 0;

end
