function normSignal = norm_sig(signal,dataPoints)
%% DESCRIPTION:
%
% normalized_signal = norm_sig(signal,dataPoints)
%
% Normalizes a signal to lower dimension specified by the amount of data
% points desired and the corresponding sampling rate of the original signal
%
%   Author: Benjamin Dourthe
%
%   Last update: April 7th, 2018
%
%% Input: 
%   - signal: original signal to normalize
%   - samplingRate: sampling frequency of the original signal (how many
%   data points per second captured)
%   - dataPoints: how many data points should compose the normalized signal
%
%% Output:
%   - normSignal: normalized signal
%
%% Dependencies:
%   Functions:
%       - None
%   Files:
%       - None

signalLength = length(signal(~isnan(signal)));
signal = signal(1:signalLength);

x0 = 0:1:signalLength-1;
% Determines size of step array
x1 = 0:signalLength/dataPoints:signalLength-1;
% Interpolate the original signal to generate a normalized one
normSignal = interp1(x0,signal,x1,'pchip');

end