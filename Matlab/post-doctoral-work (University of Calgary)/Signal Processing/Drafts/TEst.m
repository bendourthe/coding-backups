clear all; close all;
% Current directory (where the data are stored for one participant)
cd 'C:\Users\bdour\OneDrive\Work\Calgary\Projects\Dr. Scholls\Data\Pilot\';
pathName = 'C:\Users\bdour\OneDrive\Work\Calgary\Projects\Dr. Scholls\Data\Pilot';

%% Signal information
dataFull = read_data(pathName,'\1_00900_emg_1.emg',4);
signal = dataFull(:,1);   % Retrieve accelerometer data from the full dataset
    sampling_frequency = 2400;  % Sampling frequency
    dt = 1/sampling_frequency;
    
%% Filtering parameters
    % Studied step
    stepIdx = 2;
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
    Ts = 1/sampling_frequency;
    L = length(signal);
    t = linspace(0, L, L)*Ts;

    % Calculates the Fast Fourrier Transform of the signal
    Facc = fft(signal)/L;          
    Fv = linspace(0, 1, fix(L/2)+1)*Fn;
    Iv = 1:length(Fv);

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
    % before firt high peak) (first step removed)
    startPoint = firstPeakIdx+frameLength-500;
    
    % Crops the signal to make sure the length of the signal is equal to:
    %       frameLength x numSteps 
    croppedSignal = signal(startPoint:startPoint-1+frameLength*(numSteps-2));
    
    % Calculates the number of steps included in the cropped signal:
    includedSteps = length(croppedSignal)/frameLength;

% Divides and reshapes the signal into frames
reshapedSignal = reshape(croppedSignal,frameLength,includedSteps);


%% STEP 3: Wavelet filtering of one frame (step) for heel strike detection

    % Isolate a portion of the accelerometer data (just 1 step)
    if stepIdx ==1
        data = croppedSignal(1:frameLength);    
    else
        data = croppedSignal(frameLength*(stepIdx-1):frameLength*stepIdx);
    end
    
    data = data-mean(data);    % Moving the signal to be aligned with the y = 0 axis
    ld = length(data);

    % Test whether the length of the data is a power of 2 (needed to use the
    % fft_real_vvt function)
    if mod(ld,2) == 1
        data = data(1:ld-1);
        ld = ld - 1;
    else
    end

    % Defines the time scale according to the frequency rate
    time = 0+1/sampling_frequency:1/sampling_frequency:ld/sampling_frequency-1/sampling_frequency;

    df = sampling_frequency/ld;
    ind = (1:ld/2+1);
    freq = (ind - 1)*df;
    p_data = fft_real_vvt(data,2);
    
    % Plots the power spectrum of the raw signal
    figure(2);
    plot(freq,p_data)

    % Find cutting high frequency
    [pks,locs] = findpeaks(data);
    mD = mean(diff(locs));                  % Calculates the mean distance between peaks
    chf = hff*sampling_frequency/mD;        % Calculates the cutting high frequency

    % High frequency filter
    [hf_signal,hf_test] = WaveletFilter(data',dt,chf,hf_mode,hf_type);

    % Find cutting low frequency
    [M1,I1] = max(data);            % Finds location of the highest peak (around where heel strike happens)
    data_crop = data(I1+200:end);  % Crops the raw accelerometer signal to remove the heel strike section
    [M2,I2] = max(data_crop);       % Finds the location of the highest peak on the cropped signal
    [m,i] = min(data_crop);       % Finds the location of the lowest peak on the cropped signal
    clf = lff*sampling_frequency/(abs(I2-i));        % Calculates the cutting low frequency

    % Low frequency filter
    [lf_signal,lf_test] = WaveletFilter(hf_signal,dt,clf,lf_mode,lf_type);
    lf_signal = lf_signal;
    lf_signal = lf_signal-mean(lf_signal);

    % Plots the filter in the frequency space (shown in green)
    fig2 = figure(2);
    maxfig(fig2,1);
    hold on
    plot(hf_test.frequency,hf_test.value,'r')
    plot(lf_test.frequency,lf_test.value,'g')

    lf_signal_crop = lf_signal(I1-100:I1+10);
    range = mean(lf_signal)+1*std(lf_signal);   % Defines a threshold to detect actual peaks
    [A,idx]=find(lf_signal_crop>range);              % Finds the datapoints that are out of the threshold
    [B,idx2]=findpeaks(lf_signal_crop(idx));
    [C,idx3]=find(lf_signal==B(1));

    hs = lf_signal(idx3(1));                     % Finds the value of the heel strike

    % Plots the raw (blue) and filtered (green) signals, as well as the 
    % calculated location of the heel strike (red star)
    fig3 = figure(3);
    maxfig(fig3,1);
    plot(data)
    hold on
    plot(lf_signal,'g')
    hold on
    plot(idx3(1),hs,'r*','MarkerSize',10);



