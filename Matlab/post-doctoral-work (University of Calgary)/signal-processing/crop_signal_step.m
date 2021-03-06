function [croppedSignal, reshapedSignal] = crop_signal_step(signal, numSteps, frameLength, firstPeakIdx, heelStrikeLoc, w1, w2, samplingFrequency)
%% DESCRIPTION:
%
%[croppedSignal] = signap_crop_step(signal, numSteps, frameLength, firstPeakIdx)
%
% Generate a new cropped signal that will fit the filtered accelerometer
% signal generated by the heel_strike_loc.m function
%
%   Author: Benjamin Dourthe
%
%   Last update: Dec. 6th, 2017
%
%% Input: 
%   - signal: raw signal (size: m x 1)
%   - numSteps: number of steps recorded by the accelerometer
%   - frameLength: number of data points within one frame (between 2 heel
%   strikes)
%   - firstPeakIdx: indexes of the first high peak of the signal
%   - heelStrikeLoc: vector including the indexes of each detected heel
%   strike
%   - w1: how much time (ms) before heel strike for data reshape
%   - w2: how much time (ms) after heel strike for data reshape
%   - samplingFrequency: sampling frequency (in Hz)
%
%% Output:
%   - croppedSignal: cropped signal that will fit the filtered accelerometer
%   signal generated by the heel_strike_loc.m function

%% Function

    % Defines the window
    window1 = w1 * samplingFrequency / 1000;
    window2 = w2 * samplingFrequency / 1000;

    % Defines the starting point for the cropped signal (= 500 datapoints
    % before first high peak) (first step removed)
    startPoint = firstPeakIdx+frameLength-500;
    
    % Crops the signal to make sure the length of the signal is equal to:
    %       frameLength x numSteps 
    fittedSignal = signal(startPoint:startPoint-1+frameLength*(numSteps-2));
    
    % Calculates the number of steps included in the cropped signal:
    includedSteps = length(fittedSignal)/frameLength;

for i=1:includedSteps
    % Isolate a portion of the accelerometer data (just 1 step)
    if i == 1
        data = fittedSignal(1:frameLength);    
    else
        data = fittedSignal(frameLength*(i-1):frameLength*i);
    end
    
    % Normalizes data around 0
    data = data-mean(data);    
    
    % Test whether the length of the data is a power of 2 (needed to use the
    % fft_real_vvt function)
    if mod(frameLength,2) == 1
        data = data(1:frameLength-1);
        frameLength = frameLength-1;
    else
    end
    
    % Saves the cropped signal of each step and reconstruct the signal over
    % the whole time period
    if i == 1
        croppedSignal(1:frameLength) = data;
    else
        croppedSignal(frameLength*(i-1):frameLength*i) = data;
    end
end


for i = 2:includedSteps-1
    reshapedSignal(:,i) = croppedSignal(heelStrikeLoc(i)-window1:heelStrikeLoc(i)+window2);

end
 
