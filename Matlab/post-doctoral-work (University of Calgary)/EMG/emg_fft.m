function [power, intensity, tot_intensity, freq_spectrum, mnf] = emg_fft(signal, sampling_rate)
%% DESCRIPTION:
%
% [power, intensity, tot_intensity, freq_spectrum, mnfa] = emg_fft(signal, sampling_rate)
%
% Applies a Fast Fourier Transform on an EMG signal and calculates the 
%   corresponding power, intensity, total intensity, frequency spectrum and
%   mean power frequency (MNF)
%
%   Author: Benjamin Dourthe
%
%   Last update: June 14th, 2018
%
%% Input: 
%   - signal: EMG signal (size: m x 1)
%   - sampling_rate: in Hz, between 1000 and 10 000 Hz
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

    % FFT
    if mod(length(signal),2) == 1
        signal = signal(2:end,:);
    else
    end    
    power = fft_real_vvt(signal,2);
    intensity = sqrt(power');
    tot_intensity = mean(intensity);
    freq_spectrum = (1:1:length(power));
    mnf = freq_spectrum*power/sum(power);


