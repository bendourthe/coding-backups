function [numSteps, frameLength, heelStrikeIdx] = heel_strike_loc(signal, sampling_frequency)
%% DESCRIPTION:
%
%[numSteps, frameLength, heelStrikeIdx] = heel_strike_loc(signal, sampling_frequency)
%
% STEP 1: Filters the raw accelerometer data to obtain a less noisy signal
%         used to calculate the total number of steps recorded
% STEP 2: Calculates the location (index) of each heel strike within an
%         accelerometer signal (crops 1st and last steps)
%
%   Author: Benjamin Dourthe
%
%   Last update: Nov. 27th, 2017
%
%% Input: 
%   - signal: raw accelerometer signal (size: m x 1)
%   - sampling_frequency: frequency of the signal in Hz
%
%% Output:
%   - numSteps: number of steps recorded by the accelerometer
%   - frameLength: number of data points within one frame (between 2 heel
%   strikes)
%   - hellStrikeIdx: indexes of each heel strike within the raw
%   accelerometer signal

%% Function

%% STEP 1: Number of steps

% Roughly filters the signal to transform it into sinusoidal waves

Fn = sampling_frequency/2;
Ts = 1/sampling_frequency;
L = length(signal);
t = linspace(0, L, L)*Ts;

% Calculates the Fast Fourrier Transform of the signal
Facc = fft(signal)/L;          
Fv = linspace(0, 1, fix(L/2)+1)*Fn;
Iv = 1:length(Fv);

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
steps = findpeaks(filteredSignal);
numSteps = length(steps);

%% STEP 2: Reshape signal into frame (number of frames = number of steps)

% Calculates frame length (~how many data points for each steps)
frameLength = round(length(signal)/numSteps);

% Crops the signal to make sure the length of the signal is equal to:
%       frameLength x numSteps
croppedSignal = signal(1:frameLength*(numSteps-1));

% Divides and reshapes the signal into frames
reshapedSignal = reshape(croppedSignal,frameLength,numSteps-1);

% Finds the location of each peak of the reshaped signal
[peakFull,peakFullIdx] = findpeaks(reshapedSignal);

% Finds the location of each max peak of the reshaped signal
[peak,peakIdx] = max(reshapedSignal);

% Finds the location of each heel strike peak of the reshaped signal (heel
% strike peak = peak before max peak)

[maxIdx] = find(peakIdx==peakFullIDX);
[heelstrikeidx] = peakFullIdx(maxIdx-1);

% Calculates the index of each heel strike within the whole signal
for i=1:numSteps-2
    heelStrikeIdx(i) = peakIdx(i+1)+(i)*frameLength;
end