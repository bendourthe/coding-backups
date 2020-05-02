clear all; close all;
% Current directory (where the data are stored for one participant)
cd 'C:\Users\bdour\OneDrive\Work\Calgary\Projects\Dr. Scholls\Data\101\';
pathName = 'C:\Users\bdour\OneDrive\Work\Calgary\Projects\Dr. Scholls\Data\101';


%% Signal information
dataFull = read_data(pathName,'\5_101_emg_1.emg',4);
signal = dataFull(:,1);   % Retrieve accelerometer data from the full dataset
    samplingFrequency = 2400;  % Sampling frequency
    dt = 1/samplingFrequency;
    window1 = 500*samplingFrequency/1000;
    window2 = 500*samplingFrequency/1000;
    
    
%% Filtering parameters

    % High frequency filtering parameters
    hff = 1;   % High frequency factor: decreasing it will filter lower frequency data        
    hf_mode = 1;
    hf_type = 2;
    

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

    % Plots the raw and filtered signal for validity assessment
    fig1 = figure(1);
    maxfig(fig1,1);
    plot(t, signal, '-b');                           % Raw signal
    hold on
    plot(t, filteredSignal, '-r', 'LineWidth',1.5);  % Filtered signal
    hold off
    xlim([0, max(t)]);
    xlabel('Time');
    ylabel('Amplitude');
    legend('Original', 'Lowpass Filtered');

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
    chf = hff*samplingFrequency/mD;        % Calculates the cutting high frequency

    % High frequency filter
    [hf_signal,hf_test] = WaveletFilter(data',dt,chf,hf_mode,hf_type);

    % Finds location of the highest peak (around where heel strike happens)
    [M1,I1] = max(hf_signal);

    % Saves the filtered signal of each step
    if i == 1
        filteredSignal(1:frameLength) = hf_signal;
        rawSignal(1:frameLength) = data;
    else
        filteredSignal(frameLength*(i-1):frameLength*i) = hf_signal;
        rawSignal(frameLength*(i-1):frameLength*i) = data;
    end
    
    % Finds the locations and values of each heel strike
    if i == 1
        peakLoc(i) = I1;                                
        peakVal(i) = filteredSignal(peakLoc(i));                      
    else
        peakLoc(i) = I1+frameLength*(i-1)-1;              
        peakVal(i) = filteredSignal(peakLoc(i));        
    end
end

% Creates new frames centralized around the detected peaks
for i=2:includedSteps-1
    
    frame = filteredSignal(peakLoc(i)-window1:peakLoc(i)+window2);
    
    % Finds location of the highest peak on the corresponding frame
    [M1,I1] = max(frame);
    
    % Find the location of the max before the highest peak
    [M2,I2] = findpeaks(frame(I1-100:I1));
    
    % Find heel strike location and value on the corresponding frame
    [M3,I3] = find(frame==max(M2));
    
    % Find heel strike location and value on the cropped filtered signal
    hsfLoc(i) = peakLoc(i)-window1+I3-1;
    hsfVal(i) = filteredSignal(hsfLoc(i));
    
    % Find heel strike value on the cropped raw signal
    [M4, I4] = max(rawSignal(hsfLoc(i)-5:hsfLoc(i)+5));
    hsLoc(i) = hsfLoc(i)-5+I4;
    hsVal(i) = rawSignal(hsLoc(i));
    
end

%% STEP 5: Heel strike detection on original signal

% for i=1:length(hsrLoc)
%     [hsVal(i), hsLoc(i)] = find(signal==hsrVal(i));
% end



    % Plots the raw (blue) and filtered (green) signals, as well as the 
    % calculated location of each heel strikes(red star)
    
        % Defines the time scale according to the frequency rate
        timeRaw = 0+1/samplingFrequency:1/samplingFrequency:length(rawSignal)/samplingFrequency;
        timeFiltered = 0+1/samplingFrequency:1/samplingFrequency:length(filteredSignal)/samplingFrequency;
        
    fig2 = figure(2);
    maxfig(fig2,1);
    plot(timeRaw,rawSignal)
    hold on
    plot(timeFiltered,filteredSignal,'g')
    hold on
    plot(hsLoc/samplingFrequency,hsVal,'r*','MarkerSize',10);
    xlim([0 max(timeRaw)])


