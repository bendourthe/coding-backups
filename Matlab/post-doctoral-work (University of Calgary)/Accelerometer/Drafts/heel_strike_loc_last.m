function [numSteps, frameLength, firstPeakIdx, filteredSignal, heelStrikeLoc] = heel_strike_loc(signal, samplingFrequency)
%% DESCRIPTION:
%
% [numSteps, frameLength, firstPeakIdx, filteredSignal, heelStrikeLoc] = heel_strike_loc(signal, samplingFrequency)
%
% STEP 1: Calculates the number of steps by applying a strong lowpass
% filter to the raw accelerometer signal
% STEP 2: Reshapes signal into frames (number of frames = number of steps)
% STEP 3: Apply two layers of Wavelet filtering to the original signal
%           1) Cuts out high frequencies
%           2) Cuts out low frequencies
% STEP 4: Reconstructs the raw and filtered signals
% STEP 5: Calculates the position of every heel strike and plots the raw
% and filtered signal, along with the corresponding heel strikes locations
%
%   Author: Benjamin Dourthe
%
%   Last update: Dec. 11th, 2017
%
%% Input: 
%   - signal: raw accelerometer signal (size: m x 1)
%   - samplingFrequency: frequency of the signal in Hz
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

    dt = 1/samplingFrequency;

    % High frequency filtering parameters
    hff = 1;   % High frequency factor: decreasing it will filter lower frequency data        
    hf_mode = 1;
    hf_type = 2;
    
    % Low frequency filtering parameters
    lff = 3;      % Low frequency factor: increasing it will filter higher frequency data
    lf_mode = 1;
    lf_type = 4;
    
    
%% STEP 1: Number of steps

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
    i
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
    chf = hff*samplingFrequency/mD;        % Calculates the cutting high frequency

    % High frequency filter
    [hf_signal,hf_test] = WaveletFilter(data',dt,chf,hf_mode,hf_type);

    % Find cutting low frequency
    [M1,I1] = max(data);                                % Finds location of the highest peak (around where heel strike happens)
    data_crop = data(I1+1000:end);                       % Crops the raw accelerometer signal to remove the heel strike section
    [M2,I2] = max(data_crop);                           % Finds the location of the highest peak on the cropped signal
    [m1,i1] = min(data_crop);                           % Finds the location of the lowest peak on the cropped signal
    clf = lff*samplingFrequency/(abs(I2-i1));          % Calculates the cutting low frequency

    % Low frequency filter
    [lf_signal,lf_test] = WaveletFilter(hf_signal,dt,clf,lf_mode,lf_type);
    lf_signal = lf_signal-mean(lf_signal);

    lf_signal_crop = lf_signal(I1-50:I1+10);
    range = mean(lf_signal)+1*std(lf_signal);                      % Defines a threshold to detect actual peaks
    [A,idx]=find(lf_signal_crop>range);                            % Finds the datapoints that are out of the threshold
    [B,idx2]=findpeaks(lf_signal_crop(idx));
    [C,idx3]=find(lf_signal==B(1));

%% STEP 4: Reconstructed raw and signals

    % Saves the filtered signal of each step
    if i == 1
        rawSignal(1:frameLength) = data;
        filteredSignal(1:frameLength) = lf_signal;
    else
        rawSignal(frameLength*(i-1):frameLength*i) = data;
        filteredSignal(frameLength*(i-1):frameLength*i) = lf_signal;
    end
    
%% STEP 5: Heel strikes locations

    % Finds the locations and values of each heel strike
    if i == 1
        heelStrikeLoc(i) = idx3(1);                                
        heelStrikeVal(i) = lf_signal(idx3(1));                      
    else
        heelStrikeLoc(i) = idx3(1)+frameLength*(i-1)-1;              
        heelStrikeVal(i) = filteredSignal(heelStrikeLoc(i));        
    end
end

    % Plots the raw (blue) and filtered (green) signals, as well as the 
    % calculated location of each heel strikes(red star)
    
        % Defines the time scale according to the frequency rate
        timeRaw = 0+1/samplingFrequency:1/samplingFrequency:length(rawSignal)/samplingFrequency;
        timeFiltered = 0+1/samplingFrequency:1/samplingFrequency:length(filteredSignal)/samplingFrequency;
        
    fig1 = figure;
    maxfig(fig1,1);
    plot(timeRaw,rawSignal)
    hold on
    plot(timeFiltered,filteredSignal,'g')
    hold on
    plot(heelStrikeLoc/samplingFrequency,heelStrikeVal,'r*','MarkerSize',10);
    xlim([0 max(timeRaw)])