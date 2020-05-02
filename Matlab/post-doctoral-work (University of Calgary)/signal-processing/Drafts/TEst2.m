clear all; close all;
% Current directory (where the data are stored for one participant)
cd 'C:\Users\bdour\OneDrive\Work\Calgary\Projects\Dr. Scholls\Data\Pilot\';
pathName = 'C:\Users\bdour\OneDrive\Work\Calgary\Projects\Dr. Scholls\Data\Pilot';

dataFull = read_data(pathName,'\1_00900_emg_1.emg',4);
signal = dataFull(:,1);   % Retrieve accelerometer data from the full dataset
sampling_frequency = 2400;  % Sampling frequency

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
    plot(t, signal, '-b');                   % Raw signal
    hold on
    plot(t, filteredSignal, '-r', 'LineWidth',1.5);  % Filtered signal
    hold off
    xlim([0, max(t)]);
    xlabel('Time');
    ylabel('Amplitude');
    legend('Original', 'Lowpass Filtered');

    % Plots 2 buttons for validation:
        % 'Continue' if filtered signal looks okay
        b1 = uicontrol('Position',[20 20 200 40],'String','Continue',...
              'Callback','uiresume(gcbf)');
        % 'Cancel' if filtered signal looks wrong (stops code)
        b2 = uicontrol('Position',[20 60 200 40],'String','Cancel',...
              'Callback','close(gcbf)');
        uiwait(gcf);
        close(fig1);

% Calculates the number of steps based on the peaks of the filtered signal
[steps,locSteps] = findpeaks(filteredSignal);
numSteps = length(steps);


%% STEP 2: Reshape signal into frame (number of frames = number of steps)

    % Calculates frame length (~how many data points for each steps)
    frameLength = round(mean(diff(locSteps)));

    % Crops the signal to make sure the length of the signal is equal to:
    %       frameLength x numSteps
    croppedSignal = signal(1:frameLength*(numSteps-1));

% Divides and reshapes the signal into frames
reshapedSignal = reshape(croppedSignal,frameLength,numSteps-1);

%% STEP 3: Wavelet filtering of each frame for heel strike detection

    dt = 1/sampling_frequency;
    % High frequency filtering parameters
    hff = 1;   % High frequency factor: decreasing it will filter lower frequency data        
    hf_mode = 20;
    hf_type = 2;
    % Low frequency filtering parameters
    lff = 2;      % Low frequency factor: increasing it will filter higher frequency data
    lf_mode = 10;
    lf_type = 4;

    % Random number generator: defines 5 random step number that will be
    % assessed for validity
    for i=1:5
        randomStep(i) = round(rand*numSteps);
    end

    % Starts filtering the signal
    for i=1:numSteps-1
        data = reshapedSignal(:,i)-mean(reshapedSignal(:,i));    % Aligns signal with the y = 0 axis
        ld = length(data);

        % Tests whether the length of the data is a power of 2 (needed to use the
        % fft_real_vvt function)
        if mod(ld,2) == 1
            data = data(1:ld-1);
            ld = ld - 1;
        else
        end

        df = sampling_frequency/ld;
        ind = (1:ld/2+1);
        freq = (ind - 1)*df;
        p_data = fft_real_vvt(data,2);

        % Finds cutting high frequency
        [pks,locs] = findpeaks(data);
        mD = mean(diff(locs));                      % Calculates the mean distance between peaks
        chf = hff*sampling_frequency/mD;            % Calculates the cutting high frequency

            % High frequency filter
            [hf_signal,hf_test] = WaveletFilter(data',dt,chf,hf_mode,hf_type);

        % Find cutting low frequency
        [M1,I1] = max(data);                        % Finds location of the highest peak (usually right after heel strike)
        data_crop = data(I1+200:end);               % Crops the raw step signal to remove the heel strike section
        [M2,I2] = max(data_crop);                   % Finds the location of the highest peak on the cropped signal
        [m,i] = min(data_crop);                     % Finds the location of the lowest peak on the cropped signal
        clf = lff*sampling_frequency/(abs(I2-i));   % Calculates the cutting low frequency

            % Low frequency filter
            [lf_signal,lf_test] = WaveletFilter(hf_signal,dt,clf,lf_mode,lf_type);
            lf_signal = lf_signal-mean(lf_signal);

        lf_signal_crop = lf_signal(I1-40:I1);
        [A,idx1]=findpeaks(lf_signal_crop);         % Finds the indexes of the peaks that are above the threshold

        [B,idx2]=find(lf_signal==A(1));             % Finds the index of the heel strike peak (first peak)
        hs = lf_signal(idx2(1));                    % Finds the value of the heel strike
    
        % Plots the raw and filtered signal and indicate the calculated heel
        % strike location with a red star for validity
        if ismember(i,randomStep) == 1
            fig2 = figure(2)
            maxfig(fig2,1);
            plot(data)
            hold on
            plot(lf_signal,'g')
            hold on
            plot(idx2(1),hs,'r*','MarkerSize',10);
            
            % Plots 2 buttons for validation:
                % 'Continue' if filtered signal looks okay
                b1 = uicontrol('Position',[20 20 200 40],'String','Continue',...
                    'Callback','uiresume(gcbf)');
                % 'Cancel' if filtered signal looks wrong (stops code)
                b2 = uicontrol('Position',[20 60 200 40],'String','Cancel',...
                    'Callback','close(gcbf)');
                uiwait(gcf);
                close(fig2);
        else
        end
    
        % Heel strike indexes within each frame
        heelStrikeFrameIdx(i) = idx2;
    end

% Calculates the index of each heel strike within the whole signal
for i=1:numSteps-2
    heelStrikeIdx(i) = heelStrikeFrameIdx(i+1)+(i)*frameLength;
end