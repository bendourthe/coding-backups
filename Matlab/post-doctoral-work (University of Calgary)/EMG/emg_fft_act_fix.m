function [power, intensity, tot_intensity, freq_spectrum, mnf] = emg_fft_act_fix(signal, samplingFrequency, pStart, pEnd)
%% DESCRIPTION:
%
% [power, intensity, tot_intensity, freq_spectrum, mnf] = 
%    emg_fft_act_fix(signal, samplingFrequency, pStart, pEnd)
%
% Applies a Fast Fourrier Transform (FFT) on the active portion of an EMG
%   signal (using a fixed window) and calculates the corresponding power, 
%   intensity, total intensity, frequency spectrum and mean power frequency
%   (MNF)
%
%   Author: Benjamin Dourthe
%
%   Last update: June 14th, 2018
%
%% Input: 
%   - signal: EMG signal (size: m x 1)
%   - samplingFrequency: in Hz
%   - pStart: percentage of the signal when the muscle is expected to be
%   active
%   - pEnd: percentage of the signal when the muscle is expected to stop
%   being active
%
%% Output:
%   - power: power spectrum of the FFT of the signal
%   - intensity: = sqrt(power)
%   - total intensity: mean of the intensity
%   - frequency spectrum: spectrum of the different frequencies
%   - mnf: mean power frequency
%
%% Dependencies:
%
%   Function:
%       - fft_real_vvt.m

%% Function

    % Detection of active portions (i.e. when power is above threshold)
    startAct = round(pStart/100*length(signal));
    endAct = round(pEnd/100*length(signal));
    
    % Calculation of power, intensity, total intensity, and frequency 
    % spectrum during active portions of the signal
    signalAct = signal(startAct:endAct,:);
    if mod(length(signalAct),2) == 1
        signalAct = signalAct(2:end,:);
    else
    end    
    power = fft_real_vvt(signalAct,2);
    intensity = sqrt(power');
    tot_intensity = mean(intensity);
    freq_spectrum = (1:1:length(power));
    mnf = freq_spectrum*power/sum(power);


