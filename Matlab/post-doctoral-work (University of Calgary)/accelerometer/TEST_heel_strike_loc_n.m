clear all; close all;
% Current directory (where the data are stored for one participant)
cd 'D:\Dr. Scholls\00\';
pathName = 'D:\Dr. Scholls\00';


%% Signal import

dataFull = read_data(pathName,'\1_00_emg_1.emg',4);
signal = dataFull(:,1);   % Retrieve accelerometer data from the full dataset
signal = signal - mean(signal);
    
%% Parameters
    
    samplingFrequency = 2400;  % Sampling frequency
    window1 = 500;
    window2 = 100;
    
    % Windows used to centralized each frame on the highest peak detected
    window1 = window1*samplingFrequency/1000;
    

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
numSteps = length(steps) - 1;

% Calculates frame length (~how many data points for each steps)
frameLength = round(mean(diff(locSteps)));


%% STEP 2: Main peaks detection

    % First peak
        % Isolate a portion of the accelerometer data (first step)
        frame = signal(1:frameLength);    
        % Align signal with x-axis
        frame = frame-mean(frame);    
        % Finds location of the first peak
        [M,I] = max(frame);
        
        clear frame

for i=2:numSteps-1
    
    % Finds the locations and values of every other peak
        % Isolate a portion of the accelerometer data (first step)
        frame = signal(I+50+frameLength*(i-1):I+frameLength*i);    
        % Align signal with x-axis
        frame = frame-mean(frame);    
        % Find peaks
        [Mp,Ip] = max(frame);
        peakLoc(i) = Ip + I + 50 + frameLength*(i-1) -1;
        peakVal(i) = signal(peakLoc(i));
        
        clear frame
        
end

%% STEP 3: Remove bad steps

idxOut1 = find(peakVal<0.5*mean(peakVal));
peakVal(idxOut1) = [];
peakLoc(idxOut1) = [];

%% STEP 4: Heel strike(s) detection

% Creates new frames centralized around the detected peaks

for i=1:length(peakLoc)
    
    % Find the location of the max before the highest peak
    [M1,I1] = findpeaks(signal(peakLoc(i)-window2:peakLoc(i)));
    
    % Find heel strike location on the corresponding frame
    [M2,I2] = max(M1);
    
    % Find heel strike location and value on the original signal
    hsLoc(i) = I1(I2)+peakLoc(i)-window2 - 1;
    hsVal(i) = signal(hsLoc(i));

end

numSteps = length(hsLoc);;            % Number of steps corresponding to the number of heel strikes detected
hsLoc = hsLoc(2:length(hsLoc));     % Heel strike(s) location
hsVal = hsVal(2:length(hsVal));     % Signal amplitude at heel strike

%% STEP 5: Remove outliers

for i = 1:length(hsVal)
    if (peakVal(i+1)-hsVal(i)) < 0.2*peakVal(i+1)
        idxOut2(i) = i;
    else
        idxOut2(i) = 0;
    end
end
idxOut2(idxOut2==0) = [];

    hsVal(idxOut2) = [];
    hsLoc(idxOut2) = [];
    peakVal(idxOut2+1) = [];
    peakLoc(idxOut2+1) = [];
    numSteps = numSteps - length(idxOut2);

    % Plots the raw (blue) and filtered (green) signals, as well as the 
    % calculated location of each heel strikes(red star)
    
        % Defines the time scale according to the frequency rate
        time = 0+1/samplingFrequency:1/samplingFrequency:length(signal)/samplingFrequency;
        
    fig2 = figure(2);
    maxfig(fig2,1);
    plot(time,signal)
    hold on
    plot(hsLoc/samplingFrequency,hsVal,'r*','MarkerSize',10);
    xlim([0 max(time)])
    hold on
    plot(peakLoc/samplingFrequency,peakVal,'g*','MarkerSize',10);
