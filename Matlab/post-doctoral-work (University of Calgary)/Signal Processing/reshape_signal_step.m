function [reshapedSignal] = reshape_signal_step(signal, numSteps, hsLoc, w1, w2, samplingFrequency)
%% DESCRIPTION:
%
%[reshapedSignal] = reshape_signal_step(signal, numSteps, hsLoc, w1, w2, samplingFrequency)
%
% Divide the original signal in n frames (where n = number of steps)
%
%   Author: Benjamin Dourthe
%
%   Last update: Jan. 09th, 2018
%
%% Input: 
%   - signal: raw signal (size: m x 1)
%   - numSteps: number of steps recorded by the accelerometer
%   - hsLoc: vector including the indexes of each detected heel strike
%   - w1: how much time (ms) before heel strike for data reshape
%   - w2: how much time (ms) after heel strike for data reshape
%   - samplingFrequency: sampling frequency (in Hz)
%
%% Output:
%   - reshapedSignal: new matrix (size m x n) where the ith line corresponds 
%   to a frame of the full original signal centralized around the ith step
%       m = number of data points included (defined by w1 and w2)
%       n = number of steps (calculated by heel_strike_loc.m)

%% Function

    % Defines the window
    window1 = w1 * samplingFrequency / 1000;
    window2 = w2 * samplingFrequency / 1000;

    for i=1:numSteps
        reshapedSignal(:,i) = signal(hsLoc(i)-window1:hsLoc(i)+window2);
    end
 
