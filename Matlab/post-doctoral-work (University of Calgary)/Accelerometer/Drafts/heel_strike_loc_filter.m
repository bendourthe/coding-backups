function [numSteps, frameLength, firstPeakIdx, filteredSignal, hsLoc] = heel_strike_loc(signal, samplingFrequency)
%% DESCRIPTION:
%
% [numSteps, frameLength, firstPeakIdx, filteredSignal, heelStrikeLoc] = heel_strike_loc(signal, samplingFrequency)
%
% STEP 1: Calculates the number of steps by applying a strong lowpass
% filter to the raw accelerometer signal
% STEP 2: Calculates main peaks location (following heel strike)
% STEP 3: Calculates heel strike location (crops first and last step)
%
%   Author: Benjamin Dourthe
%
%   Last update: Dec. 11th, 2017
%
%% Input: 
%   - signal: raw accelerometer signal (size: m x 1)
%   - samplingFrequency: frequency of the signal in Hz
%
%% Output:
%   - numSteps: number of steps recorded by the accelerometer
%   - frameLength: number of data points within one frame (between 2 heel
%   strikes)
%   - hellStrikeIdx: indexes of each heel strike within the raw
%   accelerometer signal
%   - firstPeakIdx: indexes of the first high peak of the signal
%   - filteredSignal: filtered accelerometer signal
%   - heelStrikeLoc: vector including the indexes of each detected heel
%   strike
%   - figure: plot of the raw (blue) and filtered (green) signals with red
%   star indicating the calculated location of each heel strike

%% Function

%% Parameters
    
    % Windows used to centralized each frame on the highest peak detected
    window1 = 500*samplingFrequency/1000;
    window2 = 500*samplingFrequency/1000;
    

%% STEP 1: Number of steps

    % Roughly filters the signal to transform it into sinusoidal waves

    Fn = samplingFrequency/2;
    Ts = 1/samplingFrequency;
    L = length(signal);
    t = linspace(0, L, L)*Ts;

    Wp = 1/Fn;                               % Passband Frequencies (Normalised)
    Ws = 1.1/Fn;                             % Stopband Frequencies (Normalised)
    Rp = 10;                                 % Passband Ripple (dB)
    Rs = 50;                                 % Stopband Ripple (dB)
    [n,Ws] = cheb2ord(Wp,Ws,Rp,Rs);          % Filter Order
    [z,p,k] = cheby2(n,Rs,Ws);               % Filter Design
    [sosbp,gbp] = zp2sos(z,p,k);             % Convert To Second-Order-Section For Stability

    % Filters the signal
    filteredSignal = filtfilt(sosbp,gbp, signal);

% Calculates the number of steps based on the peaks of the filtered signal
[steps,locSteps] = findpeaks(filteredSignal);
numSteps = length(steps);

% Calculates frame length (~how many data points for each steps)
frameLength = round(mean(diff(locSteps)));


%% STEP 2: Main peaks detection

for i=1:numSteps-1
    
    % Isolate a portion of the accelerometer data (just 1 step)
    if i == 1
        frame = signal(1:frameLength);    
    else
        frame = signal(frameLength*(i-1):frameLength*i);
    end
    
    frame = frame-mean(frame);    % Moving the signal to be aligned with the y = 0 axis

    % Finds location of the highest peak (around where heel strike happens)
    [M,I] = max(frame);
    
    % Finds the locations and values of peak
    if i == 1
        peakLoc(i) = I;                                
        peakVal(i) = signal(peakLoc(i));                      
    else
        peakLoc(i) = I+frameLength*(i-1)-1;              
        peakVal(i) = signal(peakLoc(i));        
    end
end


%% STEP 3: Heel strike(s) detection

% Creates new frames centralized around the detected peaks
for i=2:numSteps-1
    
    clear frame
    
    frame = signal(peakLoc(i)-window1:peakLoc(i)+window2);
    
    % Finds location of the highest peak on the corresponding frame
    [M1,I1] = max(frame);
    
    % Find the location of the max before the highest peak
    [M2,I2] = findpeaks(frame(I1-100:I1));
    
    % Find heel strike location on the corresponding frame
    [M3,I3] = max(M2);
    I3 = I2(I3);
    I4 = I1-100 + I3-1;
    
    % Find heel strike location and value on the original signal
    hsLoc(i) = peakLoc(i)-window1+I4-1;
    hsVal(i) = signal(hsLoc(i));

end

hsLoc = hsLoc(2:length(hsLoc));     % Heel strike(s) location
hsVal = hsVal(2:length(hsVal));     % Signal amplitude at heel strike

