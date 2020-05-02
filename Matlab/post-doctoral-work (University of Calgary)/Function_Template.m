function [smoothedX, smoothedY] = smoothData(x,y,samplingRateIncrease)
%% DESCRIPTION:
%
% [smoothedX] = smoothX(x,samplingRateIncrease)
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
%   - Smooth plot
%
%% Dependencies:
%   Functions:
%       - None
%   Files:
%       - None


%% FUNCTION  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = cfs;
    y = dataGL1(1,:);
    samplingRateIncrease = 10;
    newXSamplePoints = linspace(min(x), max(x), length(x) * samplingRateIncrease);
    smoothedY = spline(x, y, newXSamplePoints);
    % Now flip back
    xSmooth = newXSamplePoints;
    ySmooth = smoothedY;
    
