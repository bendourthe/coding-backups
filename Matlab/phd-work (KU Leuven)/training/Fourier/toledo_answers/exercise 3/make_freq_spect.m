function [data_ft,frequency_out,magnitude_out] = make_freq_spect( data,fs,bool_plot)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% opmerking Maarten: op lijn 13: waarom delen door sze en vermenigvuldigen
% met 2?

% opmerking waarom absolute waarde nemen van complex getal (amplitude)
%abs(X) returns an array Y such that each element of Y is the absolute value of the corresponding element of X.
%If X is complex, abs(X) returns the complex modulus (magnitude), which is the same as
%sqrt(real(X).^2 + imag(X).^2)

data_ft=fft(data);
sze = length(data);
ff= fix(sze/2) + 1;
f = [0:ff-1]*fs/sze;
if bool_plot
    figure
    plot(f(1:ff), abs(data_ft(1:ff)/sze*2));
    xlabel('Frequency in Hz');
    ylabel('Magnitude');
    title('Frequency spectrum');
    axis tight;
end
frequency_out=f(1:ff)';
magnitude_out=abs(data_ft(1:ff)/sze*2);
end

