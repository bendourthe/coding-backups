function fSignal = lp_filter(signal, samplingRate)
%% DESCRIPTION:
%
% fSignal = lp_filter(signal, samplingRate)
%
% Uses a lowpass filter to remove noise from sinusoidal signal
%
%   Author: Star Strider (edited by Benjamin Dourthe)
%
%   Last update: Jul. 16th, 2017
%
%% Input: 
%   - signal: original signal (size: m x 1)
%   - samplingRate: frequency of the signal in Hz
%
%% Output:
%   - fSignal: lowpass filtered signal
%
%% Dependencies:
%   Functions:
%       - None
%   Files:
%       - None

%% Function

    % Signal properties
    Fn = samplingRate/2;

    % Filter parameters
    Wp = 1/Fn;                          % Passband Frequencies (Normalised)
    Ws = 1.1/Fn;                        % Stopband Frequencies (Normalised)
    Rp = 10;                            % Passband Ripple (dB)
    Rs = 50;                            % Stopband Ripple (dB)
    [n,Ws] = cheb2ord(Wp,Ws,Rp,Rs);     % Filter Order
    [z,p,k] = cheby2(n,Rs,Ws);          % Filter Design
    [sosbp,gbp] = zp2sos(z,p,k);        % Convert To Second-Order-Section For Stability

    % Filtered signal
    fSignal = filtfilt(sosbp,gbp, signal);