clear all; close all;
% Current directory (where the data are stored for one participant)
cd 'D:\Dr. Scholls\201\';
pathName = 'D:\Dr. Scholls\201';


%% Signal import

dataFull = read_data(pathName,'\4_201_emg_1.emg',4);
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
numSteps = length(steps);

% Calculates frame length (~how many data points for each steps)
frameLength = round(mean(diff(locSteps)));


%% STEP 2: Main peaks detection

for i=1:numSteps-1
    
    % Isolate a portion of the accelerometer data (just 1 step)
    if i == 1
        frame = signal(1:frameLength);    
    else
        frame = signal(frameLength*(i-1):frameLength*i);
    end
    
    frame = frame-mean(frame);    % Moving the signal to be aligned with the y = 0 axis

    % Finds location of the highest peak (around where heel strike happens)
    [M,I] = max(frame);
    
    % Finds the locations and values of peak
    if i == 1
        peakLoc(i) = I;                                
        peakVal(i) = signal(peakLoc(i));                      
    else
        peakLoc(i) = I+frameLength*(i-1)-1;              
        peakVal(i) = signal(peakLoc(i));        
    end
end

%% STEP 3: Remove outliers

idxOut = find(peakVal<0.5*mean(peakVal))
peakVal(idxOut) = [];
peakLoc(idxOut) = [];
numSteps = numSteps - length(idxOut);


%% STEP 4: Heel strike(s) detection

% Creates new frames centralized around the detected peaks

for i=2:numSteps-1
    
    clear frame
    
    frame = signal(peakLoc(i)-window1:peakLoc(i)+window1);
    
    % Finds location of the highest peak on the corresponding frame
    [M1,I1] = max(frame);
    
    % Find the location of the max before the highest peak
    [M2,I2] = findpeaks(frame(I1-window2:I1));
    
    % Find heel strike location on the corresponding frame
    [M3,I3] = max(M2);
    I3 = I2(I3);
    I4 = I1-window2 + I3-1;
    
    % Find heel strike location and value on the original signal
    hsLoc(i) = peakLoc(i)-window1 + I4-1;
    hsVal(i) = signal(hsLoc(i));

end

numSteps = numSteps - 2;            % Number of steps corresponding to the number of heel strikes detected
hsLoc = hsLoc(2:length(hsLoc));     % Heel strike(s) location
hsVal = hsVal(2:length(hsVal));     % Signal amplitude at heel strike

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
