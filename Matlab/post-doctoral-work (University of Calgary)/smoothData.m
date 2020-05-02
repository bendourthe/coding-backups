function [smoothedX, smoothedY] = smoothData(x,y,samplingRateIncrease)
%% DESCRIPTION:
%
% [smoothedX, smoothedY] = smoothData(x,y,samplingRateIncrease)
%
% Artificially increases the sampling rate to increase the number of data
% points and generate smoother line plots
%
%   Author: Benjamin Dourthe
%
%   Last update: August 22nd, 2018
%
%% Input: 
%   - x: original x variable
%   - y: original y variable
%   - samplingRateIncrease: factor by which the sampling rate will be
%   increased (the higher the smoother, start with 10)
%
%% Output:
%   - smoothedX: Smoothed X variable
%   - smoothedY: Smoothed Y variable
%
%% Dependencies:
%   Functions:
%       - None
%   Files:
%       - None


%% FUNCTION  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

smoothedX = linspace(min(x), max(x), length(x) * samplingRateIncrease);
smoothedY = spline(x, y, smoothedX);

