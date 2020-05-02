function [numSteps, frameLength, firstPeakIdx, filteredSignal, heelStrikeLoc] = heel_strike_loc(signal, sampling_frequency)
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
%   Last update: Dec. 5th, 2017
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
%   - firstPeakIdx: indexes of the first high peak of the signal
%   - filteredSignal: filtered accelerometer signal
%   - heelStrikeLoc: vector including the indexes of each detected heel
%   strike
%   - figure: plot of the raw (blue) and filtered (green) signals with red
%   star indicating the calculated location of each heel strike

%% Function

%% Filtering parameters
    
    dt = 1/ sampling_frequency;
    
    % High frequency filtering parameters
    hff = 1;   % High frequency factor: decreasing it will filter lower frequency data        
    hf_mode = 1;
    hf_type = 2;
    
    % Low frequency filtering parameters
    lff = 1;      % Low frequency factor: increasing it will filter higher frequency data
    lf_mode = 1;
    lf_type = 4;
    

%% STEP 1: Number of steps

    % Roughly filters the signal to transform it into sinusoidal waves

    Fn = sampling_frequency/2;

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
numSteps = length(steps);

% Calculates frame length (~how many data points for each steps)
frameLength = round(mean(diff(locSteps)));

% Finds the location of the first peak (~first heel strike)
[firstPeak,firstPeakIdx] = max(signal(1:frameLength));


%% STEP 2: Reshape signal into frames (number of frames = number of steps)

    % Defines the starting point for the cropped signal (= 500 datapoints
    % before first high peak) (first step removed)
    startPoint = firstPeakIdx+frameLength-500;
    
    % Crops the signal to make sure the length of the signal is equal to:
    %       frameLength x numSteps 
    croppedSignal = signal(startPoint:startPoint-1+frameLength*(numSteps-2));
    
    % Calculates the number of steps included in the cropped signal:
    includedSteps = length(croppedSignal)/frameLength;


%% STEP 3: Wavelet filtering of one frame (step) for heel strike detection

filteredSignal = zeros(1,length(croppedSignal));
for i=1:includedSteps
    % Isolate a portion of the accelerometer data (just 1 step)
    if i == 1
        data = croppedSignal(1:frameLength);    
    else
        data = croppedSignal(frameLength*(i-1):frameLength*i);
    end
    
    data = data-mean(data);    % Moving the signal to be aligned with the y = 0 axis
    ld = length(data);

    % Test whether the length of the data is a power of 2 (needed to use the
    % fft_real_vvt function)
    if mod(ld,2) == 1
        data = data(1:ld-1);
        ld = ld - 1;
        frameLength = frameLength-1;
    else
    end

    % Find cutting high frequency
    [pks,locs] = findpeaks(data);
    mD = mean(diff(locs));                  % Calculates the mean distance between peaks
    chf = hff*sampling_frequency/mD;        % Calculates the cutting high frequency

    % High frequency filter
    [hf_signal,hf_test] = WaveletFilter(data',dt,chf,hf_mode,hf_type);

    % Find cutting low frequency
    [M1,I1] = max(data);                                % Finds location of the highest peak (around where heel strike happens)
    data_crop = data(I1+1000:end);                       % Crops the raw accelerometer signal to remove the heel strike section
    [M2,I2] = max(data_crop);                           % Finds the location of the highest peak on the cropped signal
    [m1,i1] = min(data_crop);                           % Finds the location of the lowest peak on the cropped signal
    clf = lff*sampling_frequency/(abs(I2-i1));          % Calculates the cutting low frequency

    % Low frequency filter
    [lf_signal,lf_test] = WaveletFilter(hf_signal,dt,clf,lf_mode,lf_type);
    lf_signal = lf_signal-mean(lf_signal);

    lf_signal_crop = lf_signal(I1-50:I1+10);
    range = mean(lf_signal)+1*std(lf_signal);                      % Defines a threshold to detect actual peaks
    [A,idx]=find(lf_signal_crop>range);                            % Finds the datapoints that are out of the threshold
    [B,idx2]=findpeaks(lf_signal_crop(idx));
    [C,idx3]=find(lf_signal==B(1));

    % Saves the filtered signal of each step
    if i == 1
        filteredSignal(1:frameLength) = lf_signal;
        rawSignal(1:frameLength) = data;
    else
        filteredSignal(frameLength*(i-1):frameLength*i) = lf_signal;
        rawSignal(frameLength*(i-1):frameLength*i) = data;
    end
    
    % Finds the locations and values of each heel strike
    if i == 1
        heelStrikeLoc(i) = idx3(1);                                
        heelStrikeVal(i) = lf_signal(idx3(1));                      
    else
        heelStrikeLoc(i) = idx3(1)+frameLength*(i-1)-1;              
        heelStrikeVal(i) = filteredSignal(heelStrikeLoc(i));        
    end
end