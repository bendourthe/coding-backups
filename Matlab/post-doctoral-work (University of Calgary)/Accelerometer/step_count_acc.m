function numSteps = step_count_acc(signal, sampling_frequency)
%% DESCRIPTION:
%
% numSteps = step_count_acc(signal, sampling_frequency)
%
% Filters the raw accelerometer data to obtain a less noisy signal used to
% calculate the total number of steps recorded
%
%   Author: Benjamin Dourthe
%
%   Last update: Nov. 27th, 2017
%
%% Input: 
%   - signal: accelerometer signal (size: m x 1)
%   - sampling_rate: in Hz, between 1000 and 10 000 Hz
%
%% Output:
%   - numSteps: number of steps during acquisition

%% Function
acc = signal;
Fs = sampling_frequency;
Fn = Fs/2;
Ts = 1/Fs;
L = length(acc);
t = linspace(0, L, L)*Ts;

% Calculates the Fast Fourrier Transform of the signal
Facc = fft(acc)/L;          
Fv = linspace(0, 1, fix(L/2)+1)*Fn;
Iv = 1:length(Fv);

Wp = 1/Fn;                               % Passband Frequencies (Normalised)
Ws = 1.1/Fn;                             % Stopband Frequencies (Normalised)
Rp = 10;                                 % Passband Ripple (dB)
Rs = 50;                                 % Stopband Ripple (dB)
[n,Ws] = cheb2ord(Wp,Ws,Rp,Rs);          % Filter Order
[z,p,k] = cheby2(n,Rs,Ws);               % Filter Design
[sosbp,gbp] = zp2sos(z,p,k);             % Convert To Second-Order-Section For Stability

% Filtered signal
filtSignal = filtfilt(sosbp,gbp, acc);

% Calculate the number of steps based on the filtered signal
steps = findpeaks(filtSignal);
numSteps = length(steps);

