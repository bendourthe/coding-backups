function [power, intensity, tot_intensity, freq_spectrum, pa] = emg_wavelet_act_fix(signal, ...
    sampling_rate, scale, control, nr_wavelets, pStart, pEnd)
%% DESCRIPTION:
%
% [power, intensity, tot_intensity, freq_spectrum,pa] = 
%    emg_wavelet_act_fix(signal,sampling_rate, new_sampling_rate, scale, control, nr_wavelets, pStart, pEnd)
%
% Applies a wavelet filter on an EMG signal and then isolates the active
%   portions of the signal using a fixed window and calculates the 
%   corresponding power, intensity, total intensity, and frequency spectrum
%
%   Author: Benjamin Dourthe
%
%   Last update: June 25th, 2018
%
%% Input: 
%   - signal: EMG signal (size: m x 1)
%   - sampling_rate: in Hz, between 1000 and 10 000 Hz
%   - scale: defines the frequency range (default: 0.3, allows max
%   frequency ~500Hz; 0.15: allows max frequency ~250Hz; etc.)
%   - control: keep it at 3 for no effect
%   - nr_wavelets: default 13 (more increase the frequency resolution)
%   - pStart: percentage of the signal when the muscle is expected to be
%   active
%   - pEnd: percentage of the signal when the muscle is expected to stop
%   being active
%
%% Output:
%   - power: Gauss filtered power of the complex wavelet transformed signal
%   - intensity: = sqrt(power) (size: number of wavelets x m)
%   - total intensity: sum of the intensities of each wavelet over time = sum(intensity)
%   - frequency spectrum: sum of each data point that represents the
%   evolution of the intensity in frequency space = sum(intensity.')
%
%% Dependencies:
%
%   Function:
%       - wavelet_transform_V70419.m

%% Function

    % Wavelet transform
    [power, cwt, pa] = wavelet_transform_V70419(signal, sampling_rate, sampling_rate,...
        scale, control, nr_wavelets);
    
    % Detection of active portions (i.e. when power is above threshold)
    startAct = round(pStart/100*length(signal));
    endAct = round(pEnd/100*length(signal));
    
    % Calculation of power, intensity, total intensity, and frequency 
    % spectrum during active portions of the signal
    power = power(startAct:endAct,:);
    intensity = sqrt(power');
    tot_intensity = sum(intensity(2:end,:));
    freq_spectrum = sum(intensity.');


