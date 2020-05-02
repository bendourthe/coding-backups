function [filteredSignal, figFilteredSignal] = acc_filter(signal, sampling_frequency)
%% DESCRIPTION:
%
% [filteredSignal, figFilteredSignal] = acc_filter(signal, sampling_frequency)
%
% Filters the raw accelerometer data to obtain a less noisy signal (used to
% calculate the total number of steps recorded)
%
%   Author: Star Strider (edited by Benjamin Dourthe)
%
%   Last update: Jul. 16th, 2017
%
%% Input: 
%   - signal: raw accelerometer signal (size: m x 1)
%   - sampling_frequency: frequency of the signal in Hz
%
%% Output:
%   - filteredSignal: lowpass filtered signal
%   - figFilteredSignal: figure showing the original signal and the
%   filtered signal

%% Function

Fn = sampling_frequency/2;
Ts = 1/sampling_frequency;
L = length(signal);
t = linspace(0, L, L)*Ts;

figure(1)
plot(t, signal);
grid

Facc = fft(signal)/L;
Fv = linspace(0, 1, fix(L/2)+1)*Fn;
Iv = 1:length(Fv);

figure(2)
semilogx(Fv, abs(Facc(Iv))*2);
grid

Wp = 1/Fn;                                              % Passband Frequencies (Normalised)
Ws = 1.1/Fn;                                            % Stopband Frequencies (Normalised)
Rp = 10;                                                % Passband Ripple (dB)
Rs = 50;                                                % Stopband Ripple (dB)
[n,Ws] = cheb2ord(Wp,Ws,Rp,Rs);                         % Filter Order
[z,p,k] = cheby2(n,Rs,Ws);                              % Filter Design
[sosbp,gbp] = zp2sos(z,p,k);                            % Convert To Second-Order-Section For Stability

figure(3)
freqz(sosbp, 2^16, sampling_frequency);                  % Filter Bode Plot

s_filt = filtfilt(sosbp,gbp, signal);                   % Filter Signal
filteredSignal = s_filt;

figFilteredSignal = figure(4);
maxfig(figFilteredSignal,1);
plot(t, signal, '-b');
hold on
plot(t, s_filt, '-r', 'LineWidth',1.5);
hold off
xlim([0, max(t)]);
xlabel('Time');
ylabel('Amplitude');
legend('Original', 'Lowpass Filtered');
fileName = 'figFilteredSignal.png';
saveas(gcf,fileName);