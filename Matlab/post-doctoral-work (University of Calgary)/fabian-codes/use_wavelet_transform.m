function [power, intensity, tot_intensity] = use_wavelet_transform(signal, sampling_rate, new_sampling_rate, scale, control, nr_wavelets)
% UTILIZE VINZENZ WAVELET FUNCTION
% inputs:
%   emg signal for one muscle
%   sampling rate
%   new sampling rate
%   scale
%   control
%   number of wavelets
% outputs:
%   power
%   intensity
%   total intensity
if mod(size(signal,1),2)
    signal = signal(1:end - 1,:);
end

power = wavelet_transform_V70419(signal, sampling_rate, new_sampling_rate, scale, control, nr_wavelets);
intensity = sqrt(power');
tot_intensity = sum(intensity(2:end,:));

end

