function [J, grad] = costFunctionLG(theta, X, y)
%% DESCRIPTION:
%
% [J, grad] = costFunction(theta, X, y)
%
% Computes the cost of using theta as the parameter for logistic regression 
%   and the gradient of the cost w.r.t. to the parameters.
%
%   Author: Benjamin Dourthe
%
%   Created: October 10th, 2016
%   Last update: April 8th, 2018
%
%% Input: 
%   - theta
%   - X: training set
%   - y: corresponding state
%
%% Output:
%   - J: cost
%   - grad: gradient of the cost (same length as theta)
%
%% Dependencies:
%   Functions:
%       - sigmoid.m
%   Files:
%       - None
  
%% Function

    % Initialize some useful values
    m = length(y); % number of training examples

    % Compute J and gradient
    J = 0;
    grad = zeros(size(theta));

    % Compute the cost of a particular choice of theta
    
    for i=1:m
        J(i) = 1/m * (-y(i)*log(sigmoid(X(i,:)*theta)) - (1 - y(i))*log(1 - sigmoid(X(i,:)*theta)));
    end
    J = sum(J);

    % Compute the partial derivatives and set grad to the partial
    % derivatives of the cost w.r.t. each parameter in theta

    grad = 1/m .* (sigmoid(X*theta)-y')' * X;

end
