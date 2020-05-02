function [numSteps, hsLoc] = heel_strike_loc(signal, samplingFrequency, window1, window2)
%% DESCRIPTION:
%
% [numSteps, hsLoc] = heel_strike_loc(signal, samplingFrequency, window1, window2)
%
% STEP 1: Calculates the number of steps by applying a strong lowpass
% filter to the raw accelerometer signal
% STEP 2: Calculates main peaks location (following heel strike)
% STEP 3: Removes outliers (detected steps that are not steps)
% STEP 4: Calculates heel strike location (crops first and last step)
% STEP 5: Visualization
%
%   Author: Benjamin Dourthe
%
%   Last update: May 8th, 2018
%
%% Input: 
%   - signal: raw accelerometer signal (size: m x 1)
%   - samplingFrequency: frequency of the signal (in Hz)
%   - window1: number of data points included in the centralized frames
%   (before and after main peaks, default: 500; might need to be smaller
%   for faster signals, e.g. running)
%   - window2: : number of data points included in the window before each
%   main peak (default: 50; might need to be smaller for faster signals, 
%   e.g. running)
%
%% Output:
%   - numSteps: number of steps recorded by the accelerometer (- first and
%   last step)
%   - hsLoc: vector including the indexes of each detected heel
%   strike on the raw accelerometer signal
%   - figure: plot showing the raw accelerometer signal with green stars
%   indicating peaks, and red stars indicating heel strikes


%% Function

%% Parameters
    
    % Windows used to centralized each frame on the highest peak detected
    window1 = window1*samplingFrequency/1000;
    

%% STEP 1: Number of steps

    % Normalize signal around zero mean
    signal = signal - mean(signal);

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
numSteps = length(steps)-1;

% Calculates frame length (~how many data points for each steps)
frameLength = round(mean(diff(locSteps)));


%% STEP 2: Main peaks detection

    % First peak
        % Isolate a portion of the accelerometer data (first step)
        frame = signal(1:frameLength);    
        % Align signal with x-axis
        frame = frame-mean(frame);    
        % Finds location of the first peak
        [M,I] = max(frame);
        
        clear frame

for i=2:numSteps-1
    
    % Finds the locations and values of every other peak
        % Isolate a portion of the accelerometer data (first step)
        frame = signal(I+50+frameLength*(i-1):I+frameLength*i);    
        % Align signal with x-axis
        frame = frame-mean(frame);    
        % Find peaks
        [Mp,Ip] = max(frame);
        peakLoc(i) = Ip + I + 50 + frameLength*(i-1) -1;
        peakVal(i) = signal(peakLoc(i));
        
        clear frame
        
end

%% STEP 3: Remove bad steps

idxOut1 = find(peakVal<0.5*mean(peakVal));
peakVal(idxOut1) = [];
peakLoc(idxOut1) = [];

%% STEP 4: Heel strike(s) detection

% Creates new frames centralized around the detected peaks

for i=1:length(peakLoc)
    
    % Find the location of the max before the highest peak
    [M1,I1] = findpeaks(signal(peakLoc(i)-window2:peakLoc(i)));
    
    % Find heel strike location on the corresponding frame
    [M2,I2] = max(M1);
    
    % Find heel strike location and value on the original signal
    hsLoc(i) = I1(I2)+peakLoc(i)-window2 - 1;
    hsVal(i) = signal(hsLoc(i));

end

hsLoc = hsLoc(2:length(hsLoc));     % Heel strike(s) location
hsVal = hsVal(2:length(hsVal));     % Signal amplitude at heel strike

%% STEP 5: Remove outliers

for i = 1:length(hsVal)
    if (peakVal(i+1)-hsVal(i)) < 0.2*peakVal(i+1)
        idxOut2(i) = i;
    else
        idxOut2(i) = 0;
    end
end
idxOut2(idxOut2==0) = [];

    hsVal(idxOut2) = [];
    hsLoc(idxOut2) = [];
    peakVal(idxOut2+1) = [];
    peakLoc(idxOut2+1) = [];
    numSteps = length(hsLoc);            % Number of steps corresponding to the number of heel strikes detected


%% STEP 6: Visualization

    % Plots the raw signal (blue), as well as the calculated location of
    % each peak(green stars) and heel strike (red stars)
    
        % Defines the time scale according to the frequency rate
        time = 0+1/samplingFrequency:1/samplingFrequency:length(signal)/samplingFrequency;
        
        fig = figure;
        maxfig(fig,1);
        plot(time,signal)
        hold on
        plot(hsLoc/samplingFrequency,hsVal,'r*','MarkerSize',10);
        xlim([0 max(time)])
        hold on
        plot(peakLoc/samplingFrequency,peakVal,'g*','MarkerSize',10);

