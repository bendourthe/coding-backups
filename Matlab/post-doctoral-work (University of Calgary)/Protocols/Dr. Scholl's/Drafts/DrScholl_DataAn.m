clear all; close all;
%% DESCRIPTION:
%
% This code was designed to enable to visualization and analysis of the
% data provided via the Dr. Scholl's Research Project:
%
%   Project Title: Quantifying the effects of different insole
%   configurations on fatigue
%
%   Author: Benjamin Dourthe
%
%   Last update: Nov. 22nd, 2017
%
%% Input: 
%   - Currex data (3 static videos, 6 dynamic means)
%   - Biodex data (2 MVC trial, 2 50% trials)
%   - Pedar data (2 walking trials)
%   - EMG data (2 walking trials)
%
%% Output:
% 
%
%% Dependencies:
%
%   Files:
%       - Centroid Calculation_All Insoles_B.xlsx
%       - Pedar_Insole_area_B.xlsx
%
%   Functions:
%       - read_data.m
%       - use_wavelet_transform.m
%       - wavelet_transform_V70419.m
%       - pedarDiffPlot.m
%       - heel_strike_loc.m


%% ID and DIRECTORY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Subject ID
    subjectID = '00900';
    
% Current directory (where the data are stored for one participant)
cd 'C:\Users\bdour\OneDrive\Work\Calgary\Projects\Dr. Scholls\Data\Pilot\';
pathName = 'C:\Users\bdour\OneDrive\Work\Calgary\Projects\Dr. Scholls\Data\Pilot';


%% PROTOCOL SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Wavelet parameters
    freqRate = 2400;    % frenquence used during protocol (in Hz)
    scale = 0.3;
    control = 3;
    numWavelets = 13;
    
    % Pedar Insole
    insoleType = 'XS Insole';
    
    % Sessions
    numSessions = 5;
    session = num2str([1:numSessions].','%01d');
    
    % Condition titles
    Conditions = ['Condition 1'; 'Condition 1'; 'Condition 2'; 'Condition 2'; 'Condition 3';...
        'Condition 3'; 'Condition 4'; 'Condition 4'; 'Condition 5'; 'Condition 5'];
    Trials = ['Trial 1'; 'Trial 2'; 'Trial 1'; 'Trial 2'; 'Trial 1'; 'Trial 2';...
        'Trial 1'; 'Trial 2'; 'Trial 1'; 'Trial 2'];
    
    
%% FILENAMES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Currex
        % Static (videos)
    cs1 = [strcat('\',session,'_',subjectID,'_cs_v_1.csv')];
    cs2 = [strcat('\',session,'_',subjectID,'_cs_v_2.csv')];
    cs3 = [strcat('\',session,'_',subjectID,'_cs_v_3.csv')];
    currexStat = [cs1; cs2; cs3];
        % Dynamic (mean)
    cd0 = [strcat('\',session,'_',subjectID,'_cd_m_0.csv')];
    cd1 = [strcat('\',session,'_',subjectID,'_cd_m_1.csv')];
    cd2 = [strcat('\',session,'_',subjectID,'_cd_m_2.csv')];
    cd3 = [strcat('\',session,'_',subjectID,'_cd_m_3.csv')];
    cd4 = [strcat('\',session,'_',subjectID,'_cd_m_4.csv')];
    cd5 = [strcat('\',session,'_',subjectID,'_cd_m_5.csv')];
    currexDyn = [cd0; cd1; cd2; cd3; cd4; cd5];
    
    % Biodex
        % MVC Trials
    mvc1 = [strcat('\',session,'_',subjectID,'_mvc_1.csv')];
    mvc2 = [strcat('\',session,'_',subjectID,'_mvc_2.csv')];
    mvc = [mvc1; mvc2];
        % 50% Trials
    half1 = [strcat('\',session,'_',subjectID,'_50p_1.csv')];
    half2 = [strcat('\',session,'_',subjectID,'_50p_2.csv')];
    halfp = [half1; half2];

    % Pedar
    cellSheetPath = 'C:\Users\bdour\OneDrive\Work\Calgary\Documentation\Pedar\Centroid Calculation_All Insoles_B.xlsx';
    areaSheetPath = 'C:\Users\bdour\OneDrive\Work\Calgary\Documentation\Pedar\Pedar_Insole_area_B.xlsx';
    ped1 = [strcat('\',session,'_',subjectID,'_ped_1.asc')];
    ped2 = [strcat('\',session,'_',subjectID,'_ped_2.asc')];
    ped = [ped1; ped2];

    % EMG (Treadmill)
    emg1 = [strcat('\',session,'_',subjectID,'_emg_1.emg')];
    emg2 = [strcat('\',session,'_',subjectID,'_emg_2.emg')];
    emg = [emg1; emg2];
    

%% ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PRIMARY VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Peak TORQUE during MVC
    
    % Calculates the maximal torque for every trial (2 trials per session,
    % 5 sessions -> 10 values)
    
    for i=1:size(mvc,1)        
        
        % Load data
        data = csvread([pathName mvc(i,:)],2,1);  
        
        % Calculation
        mvcTorqueMax(i) = max(abs(data(:,1)));
        
    end
    
    % Reorganizes the values in a table (ROW = session, COLUMN = trial)
    mvcTorqueMax = [mvcTorqueMax(1:5); mvcTorqueMax(6:10)].';
    
    
%% Peak-to-peak MUSCLE ACTIVITY during MVC

    % Calculates the maximal amplitude of the EMG signal, then calculates
    % the peak-to-peak amplitude using 500 data points around the maximal
    % (100 data points before, 400 data points after).
    
    % Calculation done for every trial (2 trials per session, 5 sessions
    % -> 10 values)
    
    for i=1:size(mvc,1)
        
        % Load data
        data = csvread([pathName mvc(i,:)],2,1);
        emgData = data(:,2);
        
        % Normalizes data around 0
        emgData = emgData - mean(emgData);
        
        % Calculation
        [M, I] = max(abs(emgData));
        maxStim(i) = M;
        mvcPeak2Peak(i) = peak2peak(emgData(I-100:I+400));
        
    end
    
    % Reorganizes the values in a table (ROW = session, COLUMN = trial)
    mvcPeak2Peak = [mvcPeak2Peak(1:5); mvcPeak2Peak(6:10)].';
    
    % Plots the EMG data around the stimulation peak during MVC
    
    figPeakStim = figure;
    maxfig(figPeakStim,1);
    
    for i=1:size(mvc,1)
        
        % Defines the time scale according to the frequency rate
        time = 0:1/freqRate:length(data)/freqRate-1/freqRate;
        
        % Plots data
        subplot(numSessions,2,i)
            plot(time(I-100:I+400),emgData(I-100:I+400))
            xlim([time(I-100), time(I+400)]);
            ylim([0, max(maxStim)]);
            xlabel('Time (s)')
            ylabel('EMG Amplitude')
            title(strcat(Conditions(i,:),' - ',Trials(i,:)));
            
    end
    
    % Defines figure caption
    set(figPeakStim,'numbertitle','off','name','EMG activity during stimulation during MVC');
        fileNamePeakStim = strcat('EMG_Stimulation_MVC');
    
    % Save figure as Matlab figure
    savefig(figPeakStim,fileNamePeakStim);
    % Save figure as PNG
    saveas(gcf,'Total_Intensity_50p_contration.png');
    
    
%% TOTAL INTENSITY of GM during SUSTAINED CONTRACTION (50% trials)

    % Process the EMG signals recorded during the sustained contraction
    % with a wavelet filter and then transformed with a wavelet function to
    % be represented in time-frequency space.
    
    figTInt50p = figure;
    maxfig(figTInt50p,1);
    
    for i=1:size(halfp,1)
        
        % Load data
        data = csvread([pathName halfp(i,:)],2,1);
        emgData = data(:,2);
        
        % Normalizes data around 0
        emgData = emgData - mean(emgData);
        
        % Defines the time scale according to the frequency rate
        time = 0:1/freqRate:length(data)/freqRate-1/freqRate;
                
        % Uses Vinzenz wavelet transform to calculate power, intensity,
        % total intensity and frequency spectrum of the EMG signal
        [p50p, int50p, tInt50p, f50p] = emg_wavelet(emgData, freqRate, freqRate, 0.3, 3, 13);
        
        % Calculates the mean, max and min of the total intensity from
        % time = 5sec to time = 10sec and plots the evolution of the total
        % intensity over time
        meanTotInt50p(i) = mean(tInt50p(12001:24001));
        maxTotInt50p(i) = max(tInt50p(12001:24001));
        minTotInt50p(i) = min(tInt50p(12001:24001));
        
        % Plots data
        subplot(numSessions,2,i)
            plot(time,tInt50p)
            xlim([0, 15]);
            ylim([0, 1]);
            xlabel('Time (s)')
            ylabel('Total Intensity')
            title(strcat(Conditions(i,:),' - ',Trials(i,:)));
            
    end  
    
    % Defines figure caption
    set(figTInt50p,'numbertitle','off','name','Total Intensity of the EMG signals during the 50% contration trials');
        fileName50p = strcat('Total_Intensity_50p_contration');
        
    % Save figure as Matlab figure
    savefig(figTInt50p,fileName50p);
    % Save figure as PNG
    saveas(gcf,'Total_Intensity_50p_contration.png');
    
    % Reorganizes the values in a table (ROW = session, COLUMN = trial)
    meanTotInt_50p = [meanTotInt50p(1:5); meanTotInt50p(6:10)].';
    maxTotInt_50p = [maxTotInt50p(1:5); maxTotInt50p(6:10)].';
    minTotInt_50p = [minTotInt50p(1:5); minTotInt50p(6:10)].';


    %% TOTAL INTENSITY of GM, VL and GM during WALKING TRIALS

    % Process the EMG signals recorded during the walking trials with a
    % wavelet filter and then transformed with a wavelet function to be
    % represented in time-frequency space.
    
    % Calculates the power, intensity and total intensity for every trial,
    % then plots the evolution of the total intensity over time (5
    % conditions, 2 trials each = 10 plots), and calculate the mean, max
    % and min total intensity over a 60secs period
    
    figTIntWalk = figure;
    maxfig(figTIntWalk,1);
    
    for i=1:size(emg,1)
        
        % Load data
        data = read_data(pathName,emg(i,:),4);
        dataAcc = data(:,1);
        emgDataGM = data(:,2);
        emgDataVL = data(:,3);
        emgDataVM = data(:,4);
        
        % Normalizes data around 0
        dataAcc = dataAcc - mean(dataAcc);
        emgDataGM = emgDataGM - mean(emgDataGM);
        emgDataVL = emgDataVL - mean(emgDataVL);
        emgDataVM = emgDataVM - mean(emgDataVM);
        
        % Defines the time scale according to the frequency rate
        time = 0+1/freqRate:1/freqRate:length(data)/freqRate-1/freqRate;        
        
        % Calculates the number of steps recorded (= number of peaks
        % recorded by the accelerometer) and their corresponding location
        [numSteps, frameLength, firstPeakIdx, filteredSignal, heelStrikeLoc] = heel_strike_loc_n(dataAcc, freqRate);
        
        % Crops each signal to fit the filtered accelerometer signal
        cropDataAcc(:,j) = crop_signal_step(dataAcc, numSteps, frameLength, firstPeakIdx);
        cropEmgDataGM(:,j) = crop_signal_step(emgDataGM, numSteps, frameLength, firstPeakIdx);
        cropEmgDataVL(:,j) = crop_signal_step(emgDataVL, numSteps, frameLength, firstPeakIdx);
        cropEmgDataVM(:,j) = crop_signal_step(emgDataVM, numSteps, frameLength, firstPeakIdx);
        
        % Reshapes each signal in N frames (N = numSteps), starting from
        % the beginning of each heel strike
        for j=1:numSteps-2
            reshapedDataAcc(:,j) = dataAcc(heelStrikeIdx(j):heelStrikeIdx(j)+frameLength);
            reshapedEmgDataGM(:,j) = emgDataGM(heelStrikeIdx(j):heelStrikeIdx(j)+frameLength);
            reshapedEmgDataVL(:,j) = emgDataVL(heelStrikeIdx(j):heelStrikeIdx(j)+frameLength);
            reshapedEmgDataVM(:,j) = emgDataVM(heelStrikeIdx(j):heelStrikeIdx(j)+frameLength);
        end
        
        % Calculates the mean step for each signal
        meanDataAcc = mean(reshapedDataAcc.');
        meanEmgDataGM = mean(reshapedEmgDataGM.');
        meanEmgDataVL = mean(reshapedEmgDataVL.');
        meanEmgDataVM = mean(reshapedEmgDataVM.');

        % Uses Vinzenz wavelet transform to calculate power, intensity,
        % total intensity and frequency spectrum of each mean EMG signal
        [pGM, intGM, tIntGM, fGM] = emg_wavelet(meanEmgDataGM, freqRate, freqRate, scale, control, numWavelets);
        [pVL, intVL, tIntVL, fVL] = emg_wavelet(meanEmgDataVL, freqRate, freqRate, scale, control, numWavelets);
        [pVM, intVM, tIntVM, fVM] = emg_wavelet(meanEmgDataVM, freqRate, freqRate, scale, control, numWavelets);
        maxF(i) = max([fGM fVL fVM]);
        
        % Plots the intensity calculated for a mean step
        subplot(numSessions,2,i)
            contourf(intGM,20,'Linestyle','None')
            % Setup X ticks according to time (ms)
            xticks([0 600 1200 1800 2400])
            xticklabels({'0','250','500','750','1000'})
            % Removes Y axis label
            set(gca,'ytick',[])
        
        % Clears re-used variables
        clear newDataAcc
        clear newEmgDataGM
        clear newEmgDataVL
        clear newEmgDataVM
        
    end
    
    % Defines figure caption
    set(figTIntWalk,'numbertitle','off','name','Mean intensity of the GM for each step during walking trials');
        fileNameWalk = strcat('Mean_Intensity_Walking_GM');
        
    % Save figure as Matlab figure
    savefig(figTIntWalk,fileNameWalk);
    % Save figure as PNG
    saveas(gcf,'Mean_Intensity_Walking_GM.png');
    
    
    %% FREQUENCY SPECTRUMS of GM, VL and GM during WALKING TRIALS

    % Process the EMG signals recorded during the walking trials with a
    % wavelet filter and then transformed with a wavelet function to be
    % represented in time-frequency space.
    
    figFreqSpectrum = figure;
    for i=1:size(emg,1)
        data = read_data(pathName,emg(i,:),4);
        % Define the frequency scale according to the frequency rate and
        % number of wavelets
        frequency = [1:13]; % Temporary
        emgDataGM = data(:,2);
        emgDataVL = data(:,3);
        emgDataVM = data(:,4);
        [pGM, intGM, tIntGM, fGM] = emg_wavelet(emgDataGM, freqRate, freqRate, scale, control, numWavelets);
        [pVL, intVL, tIntVL, fVL] = emg_wavelet(emgDataVL, freqRate, freqRate, scale, control, numWavelets);
        [pVM, intVM, tIntVM, fVM] = emg_wavelet(emgDataVM, freqRate, freqRate, scale, control, numWavelets);
        % Plots the frequency spectrum of each muscle for each trial and
        % condition
            subplot(numSessions,2,i)
            plot(frequency,fGM,'r',frequency,fVL,'g',frequency,fVM,'b')
            xlim([1, 13]);
            ylim([0, max(maxF*1.10)]);
            xlabel('Wavelets')
            ylabel('Amplitude')
            title(strcat(Conditions(i,:),' - ',Trials(i,:)));
    end
    set(figFreqSpectrum,'numbertitle','off','name','Frequency spectrum of GM, VL and VM during walking trials');
        fileNameFreqSpectrum = strcat('Frequency_Spectrum');
        savefig(figFreqSpectrum,fileNameFreqSpectrum);
        
    
    %% SUMMARY MATRIX
    
    % Copy/Paste this matrix directly in the Excel sheet
    Results = [mvcTorqueMax mvcPeak2Peak minTotInt_50p meanTotInt_50p maxTotInt_50p];
    csvwrite([subjectID '_Full_Results.csv'],Results);
    
    
%% SECONDARY VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PLANTAR PRESSURE distribution DIFFERENCE (Pedar)

    % Plots the pressure distrubution measured using Pedar during the
    % Control trial (left) and one of the other conditions (middle) and
    % plots the difference between the 2 plots (right)

    % Trial 1
    for i=2:4
        map = pedarDiffPlot(cellSheetPath,areaSheetPath,insoleType,[pathName ped(1,:)],[pathName ped(i,:)]);
        set(map,'numbertitle','off','name',sprintf('Trial 1: Control VS Condition %d',i-1));
        fileName = strcat('Trial_1_',subjectID,'_Control_VS_Condition_');
        savefig(map,sprintf([fileName '%d.fig'],i-1)); % Save figures
    end
    
    % Trial 2
    for i=7:9
        map = pedarDiffPlot(cellSheetPath,areaSheetPath,insoleType,[pathName ped(6,:)],[pathName ped(i,:)]);
        set(map,'numbertitle','off','name',sprintf('Trial 2: Control VS Condition %d',i-6));
        fileName = strcat('Trial_2_',subjectID,'_Control_VS_Condition_');
        savefig(map,sprintf([fileName '%d.fig'],i-6)); % Save figures
    end
