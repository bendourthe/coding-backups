function [rampLoc] = ramp_det(signal, samplingFrequency, numRamps, pMax)
%% DESCRIPTION:
%
% [rampLoc] = ramp_det(signal, samplingFrequency, percMax)
%
% Detects the starting and ending points of a ramp function. Can be used on
% the positioning signal collected with the biodex to detect when motion
% happens and calculate the corresponding angular velocity.
%
%   Author: Benjamin Dourthe
%
%   Last update: Jan. 30th, 2018
%
%% Input: 
%   - signal: ramp signal (size: m x 1)
%   - samplingFrequency: frequency of the signal (in Hz)
%   - percMax: threshold used to determine which data points are part of
%   the ramp
%
%% Output:
%   - rampLoc: 2xm matrix where m is the number of detected ramps, with the
%   index of the starting point on the 1st row, and the index of the
%   ending point on the 2nd row


%% Function

%% Derivative

der = diff(signal(100:end-100))*samplingFrequency;

%% Filtering (removes noise)

    % Parameters
    dt = 1/samplingFrequency;
    mode = 1;
    type = 2;
    
    % Detects high frequency cut by calculating the mean distance between
    % peaks
    [pks,locs] = findpeaks(der);
    mD = mean(diff(locs));                  
    chf = samplingFrequency/mD;
    
    % Filters signal
    [filtDer,filtDer_test] = WaveletFilter(der',dt,chf,mode,type);


%% Ramps detection

for i=1:numRamps
    % Finds all data points that are above 75% of the max of the signal
    % (should include all data points that belong to the ramps)
    if i==1
        x1 = find(filtDer > pMax*max(filtDer));
        starting(i) = x1(1) + 100;
        x2 = find(filtDer(starting(i):end) < pMax*max(filtDer));
        ending(i) = starting(i) + x2(1) + 100;
    else
        x1 = find(filtDer(ending(i-1):end) > pMax*max(filtDer));
        starting(i) = ending(i-1) + x1(1) + 100;
        x2 = find(filtDer(starting(i):end) < pMax*max(filtDer));
        ending(i) = starting(i) + x2(1) + 100;
    end
end
    % Collects indexes of each ramp starting and ending points
    rampLoc = [starting;ending];
    
    