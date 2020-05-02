clear all; close all;
%% DESCRIPTION:
%
% This code was designed to enable the visualization and analysis of the
% data provided via the Dr. Scholl's Research Project:
%
%   Project Title: Quantifying the effects of different insole
%   configurations on fatigue
%
%   Author: Benjamin Dourthe
%
%   Last update: May 4th, 2018
%
%% Input: 
%   - Currex data (3 static videos, 6 dynamic means)
%   - Biodex data (2 MVC trial, 2 50% trials)
%   - EMG data (2 walking trials)
%   - Pedar data (2 walking trials)
%
%% Output:
%   - Results_ArchIndexStat: Static arch indexes for each trial
%   - Results_ArchIndexStatMean: Mean static arch indexes
%   - Results_MVC: Max torques, percentages of voluntary activation and
%                   peak to peak muscle activation during all MVC trials
%   - Results_MVCMean: Means of max torques, percentages of voluntary activation (GM) 
%                   and peak to peak muscle activation during MVC (mean between trials)
%   - Results_50p: Mean torques, muscle activity (GM) and muscle frequency (GM)
%                   during all 50% trials
%   - Results_50p: Mean torques, muscle activity (GM) and muscle frequency (GM)
%                   during all 50% trials (mean between trials)
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
    subjectID = '101';
    
% Current directory (where the data are stored for one participant)
cd 'D:\Dr. Scholls\Phase 1\101\';
pathName = 'D:\Dr. Scholls\Phase 1\101';


%% PROTOCOL SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Sessions
    numSessions = 5;
    session = num2str([1:numSessions].','%01d');
    
    % Code execution: determines which parts of the code will be run (0:
    % won't be run; 1: will be run)
    currexS = 0;
    currexD = 0;
    mvc = 1;
    half = 1;
    walking = 1;
    pedar = 0;
    
    % Window definition
    w1 = 300;   % how much time (ms) before heel strike
    w2 = 600;   % how much time (ms) after heel strike
    
    % Wavelet parameters
    samplingFrequency = 2400;    % frenquence used during protocol (in Hz)
    scale = 0.3;
    control = 3;
    numWavelets = 13;
    
    % Power threshold
    pth = 0.2;  % percentage of the max mean power (allows selection of active portion of the EMG signal)
    
    % Pedar Insole
    cellSheetPath = 'C:\Users\bdour\OneDrive\Work\Calgary\Documentation\Pedar\Centroid Calculation_All Insoles_B.xlsx';
    areaSheetPath = 'C:\Users\bdour\OneDrive\Work\Calgary\Documentation\Pedar\Pedar_Insole_area_B.xlsx';
    insoleType = 'YS Insole';
    
    
%% FILENAMES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Currex
        % Static (videos)
    cs1 = [strcat('\',session,'_',subjectID,'_csv_1.csv')];
    cs2 = [strcat('\',session,'_',subjectID,'_csv_2.csv')];
    cs3 = [strcat('\',session,'_',subjectID,'_csv_3.csv')];
    currexStat_name = [cs1; cs2; cs3];
        % Dynamic (mean)
    cd0 = [strcat('\',session,'_',subjectID,'_cdm_0.csv')];
    cd1 = [strcat('\',session,'_',subjectID,'_cdm_1.csv')];
    cd2 = [strcat('\',session,'_',subjectID,'_cdm_2.csv')];
    cd3 = [strcat('\',session,'_',subjectID,'_cdm_3.csv')];
    cd4 = [strcat('\',session,'_',subjectID,'_cdm_4.csv')];
    cd5 = [strcat('\',session,'_',subjectID,'_cdm_5.csv')];
    currexDyn_name = [cd0; cd1; cd2; cd3; cd4; cd5];
    
    % Biodex
        % MVC Trials
    mvc1 = [strcat('\',session,'_',subjectID,'_mvc_1.csv')];
    mvc2 = [strcat('\',session,'_',subjectID,'_mvc_2.csv')];
    mvc_name = [mvc1; mvc2];
        % 50p Trials
    half1 = [strcat('\',session,'_',subjectID,'_50p_1.csv')];
    half2 = [strcat('\',session,'_',subjectID,'_50p_2.csv')];
    half_names = [half1; half2];

    % Pedar
    ped1 = [strcat('\',session,'_',subjectID,'_ped_1.asc')];
    ped2 = [strcat('\',session,'_',subjectID,'_ped_2.asc')];
    pedar_name = [ped1; ped2];

    % Walking (EMG data)
    emg1 = [strcat('\',session,'_',subjectID,'_emg_1.emg')];
    emg2 = [strcat('\',session,'_',subjectID,'_emg_2.emg')];
    emg_name = [emg1; emg2];
    

%% ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CURREX Static %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Execution
if currexS == 0
else

    % Calculates the Arch Index of both feet during standing static trials
    % on the Currex plate (3 trials per session, 5 sessions -> 30 values)

    for i=1:size(currexStat_name,1)
        current_trial = i
        % Load data
            data = csvread([pathName currexStat_name(i,:)]);  
        % Calculation
            [AILeftStat(i),AIRightStat(i)] = archindexStat_function(data);
    end

    % Reorganizes the values in a table (ROW = session, COLUMN = trial)
        AILeftStat = reshape(AILeftStat,numSessions,3).';
        AIRightStat = reshape(AIRightStat,numSessions,3).';

    % Calculates the mean Arch Index between trials
        AILeftMeanStat = mean(AILeftStat);
        AIRightMeanStat = mean(AIRightStat);

    % Exports Results
        ArchIndexStat = [AILeftStat AIRightStat;AILeftMeanStat AIRightMeanStat];
        csvwrite([subjectID '_AI_Stat.csv'],ArchIndexStat);
    
    disp('     *****     DONE WITH CURREX STATIC TRIALS     *****     ')

end

%% CURREX Dynamic %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Execution
if currexD == 0
else

    % Calculates the Arch Index of both feet during standing static trials
    % on the Currex plate (3 trials per session, 5 sessions -> 30 values)

    for i=1:size(currexDyn_name,1)
        current_trial = i
        % Load data
            data = csvread([pathName currexDyn_name(i,:)]);  
        % Calculation
            AIDyn(i) = archindexDyn_function(data);
    end

    % Reorganizes the values in a table (ROW = trial, COLUMN = session)
        AIDyn = reshape(AIDyn,numSessions,6).';
    
    % Isoles left from right foot (right = rows 1, 3, 5; left = rows 2, 4, 6)
        AIRightDyn = [AIDyn(1,:); AIDyn(3,:); AIDyn(5,:)];
        AILeftDyn = [AIDyn(2,:); AIDyn(4,:); AIDyn(6,:)];

    % Calculates the mean Arch Index between trials
        AILeftMeanDyn = mean(AILeftDyn);
        AIRightMeanDyn = mean(AIRightDyn);

    % Exports Results
        ArchIndexDyn = [AILeftDyn AIRightDyn;AILeftMeanDyn AIRightMeanDyn];
        csvwrite([subjectID '_AI_Dyn.csv'],ArchIndexDyn);
    
    disp('     *****     DONE WITH CURREX DYNAMIC TRIALS     *****     ')

end
    
%% MVC TRIALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Execution
if mvc == 0
else
    
    for i=1:size(mvc_name,1)
        current_trial = i
        % Load data
            data = csvread([pathName mvc_name(i,:)],2,1);
            dataTorque = abs(data(:,1))*144.72+2.3368;
            dataEMG = data(:,2);
            
        % Peak torque during MVC
            mvcTorqueMax(i) = max(dataTorque);

        % Percentage of voluntary activation
            PVA = figure;
            maxfig(PVA,1);
            plot(dataTorque)
            hold on
            plot(dataEMG)
            disp('Manually select:')
            disp('     1) Torque right before twitch during MVC')
            disp('     2) Peak torque resulting from twitch during MVC')
            disp('     3) Min torque after MVC')
            disp('     4) Peak torque resulting from twitch after MVC')
            mvcP = ginput(4);
            mvcP = round(mvcP(:,1));
            mvcT = [min(dataTorque(mvcP(1)-50:mvcP(1)+50));...
                max(dataTorque(mvcP(2)-50:mvcP(2)+50));...
                min(dataTorque(mvcP(3)-50:mvcP(3)+50));...
                max(dataTorque(mvcP(4)-50:mvcP(4)+50))];
            if mvcT(2) < mvcTorqueMax
                mvcRatio(i) = 100 - 100*(mvcT(2) - mvcT(1))/(mvcT(4) - mvcT(3));
                mvcVAR(i) = (mvcTorqueMax(i)*mvcRatio(i))/mvcT(2);
            elseif (mvcT(2) - mvcT(1)) > (mvcT(4) - mvcT(3))
                mvcRatio(i) = 0;
                mvcVAR(i) = mvcRatio(i);
            else
                mvcRatio(i) = 100 - 100*(mvcT(2) - mvcT(1))/(mvcT(4) - mvcT(3));
                mvcVAR(i) = mvcRatio(i);
            end
            if mvcRatio(i) > 100
                mvcRatio(i) = 100;
            end
            if mvcVAR(i) >100
                mvcVAR(i) = 100;
            end
            
            % Plots torque data with the selected point for verification
                figure
                clf
                hold on
                plot(dataTorque)
                plot(mvcP,mvcT,'r*','MarkerSize',10)
                disp('PRESS ENTER if points are correctly selected')
                disp('PRESS Ctr+C in the command window if not')
                pause
                
        % Peak-to-peak muscle activation
            % Select dataPoints around twitch (e.g. 2400 data points
            % corresponds to 1 sec before and 1 sec after the twitch)
                dataPoints = 2400;
            % Calculates the index of the first stimulation (happens within the
            % first 5 secs of the signal -> 12000 data points)
                [M, I1] = max(dataEMG(1:12000));
            % Calculates the index and amplitude of the second stimulation (during MVC,
            % happens within 10 secs following the first stimulation -> 24000
            % data points)
                [mvcStimMax(i), I2M(i)] = max(dataEMG(I1+2000:I1+20000));
                [mvcStimMin(i), I2m(i)] = min(dataEMG(I1+2000:I1+20000));
                if mvcStimMax(i)>abs(mvcStimMin(i))
                    I2(i) = I2M(i);
                else
                    I2(i) = I2m(i);
                end
                I2(i) = I2(i) + I1 + 2000;
            % Calculates the peak to peak amplitude of the stimulation during
            % the MVC
                mvcP2P(i) = peak2peak(dataEMG(I2(i)-dataPoints:I2(i)+dataPoints));
            % Plots torque data with the selected point for verification
                figure
                clf
                hold on
                plot(dataEMG)
                plot([I2(i)-dataPoints:I2(i)+dataPoints],dataEMG(I2(i)-dataPoints:I2(i)+dataPoints),'r')
                disp('PRESS ENTER if points are correctly selected')
                disp('PRESS Ctr+C in the command window if not')
                pause
            
            close all

    end
    
    % Reorganizes the values in a table (ROW = session, COLUMN = trial)
        mvcTorqueMax = reshape(mvcTorqueMax,numSessions,2).';
        mvcRatio = reshape(mvcRatio,numSessions,2).';
        mvcVAR = reshape(mvcVAR,numSessions,2).';
        mvcP2P = reshape(mvcP2P,numSessions,2).';
    
    % Calculates the means between trials
        mvcTorqueMaxMean = mean(mvcTorqueMax);
        mvcRatioMean = mean(mvcRatio);
        mvcVARMean = mean(mvcVAR);
        mvcP2PMean = mean(mvcP2P);
    
    % Export Results
        mvcResults = [mvcTorqueMax mvcRatio mvcVAR  mvcP2P; ...
            mvcTorqueMaxMean mvcRatioMean mvcVARMean mvcP2PMean];
        csvwrite([subjectID '_MVC.csv'],mvcResults);

    disp('     *****     DONE WITH MVC TRIALS     *****     ')
        
end
    
    
%% 50% TRIALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Execution
if half == 0
else

    % Apply a Wavelet filter to the EMG signal and then transformed with a wavelet function to
    % be represented in time-frequency space.
    
    for i=1:size(half_names,1)
        
        % Load data
            data = csvread([pathName half_names(i,:)],2,1);
            dataTorque = abs(data(:,1))*144.72+2.3368;
            dataEMG = data(:,2);
        % Normalizes data around 0
            dataEMG = dataEMG - mean(dataEMG);        
        % Uses Vinzenz Wavelet transform to calculate power, intensity,
        % total intensity and frequency spectrum of the EMG signal from 
        % time = 5sec to time = 10sec
            [p50p, int50p, tInt50p, f50p, pa50p] = emg_wavelet(dataEMG,...
                samplingFrequency, scale, control, numWavelets);        
        % Calculates the mean torque and total intensity  
            meanTorque50p(i) = mean(dataTorque);    
            meanTotInt50p(i) = mean(tInt50p);            
        % Calculates the mean frequency
            MNF50p(i) = f50p*pa50p.cfs'/sum(f50p);
         
    end  
    
    % Reorganizes the values in a table (ROW = session, COLUMN = trial)
        meanTorque50p = reshape(meanTorque50p,numSessions,2).';    
        meanTotInt50p = reshape(meanTotInt50p,numSessions,2).';        
        MNF50p = reshape(MNF50p,numSessions,2).';
        
    % Calculates the means between trials
        meanTorque50pMean = mean(meanTorque50p);
        meanTotInt50pMean = mean(meanTotInt50p);        
        MNF50pMean = mean(MNF50p);
    
    % Export Results
        halfResults = [meanTorque50p meanTotInt50p MNF50p;meanTorque50pMean meanTotInt50pMean MNF50pMean];
        csvwrite([subjectID '_50p.csv'],halfResults);
    
    disp('     *****     DONE WITH 50% TRIALS     *****     ')

end

%% WALKING TRIALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% Execution
if walking == 0
else

    % Apply a Wavelet filter to the EMG signal
    
    for i=1:size(emg_name,1)
        emg_name(i,:)
        % Load data
            data = read_data(pathName,emg_name(i,:),4);
            acc = data(:,1);
            emgGM = data(:,2);
            emgVL = data(:,3);
            emgVM = data(:,4);
            
        % Heel strike detection
            [numSteps, hsLoc] = heel_strike_loc_n(acc, samplingFrequency, 500, 100);
            
        % Validation: check if steps are correctly detected
            disp('PRESS ENTER if heel strikes are correctly location')
            disp('PRESS Ctr+C in the command window if not')
            pause
            close all            
            
        % Reshape signals (ROW = data point; COLUMN = step)
            reshapedEmgGM = reshape_signal_step(emgGM, numSteps, hsLoc, w1, w2, samplingFrequency);
            reshapedEmgVL = reshape_signal_step(emgVL, numSteps, hsLoc, w1, w2, samplingFrequency);
            reshapedEmgVM = reshape_signal_step(emgVM, numSteps, hsLoc, w1, w2, samplingFrequency);
            
        % Calculate the mean step for each signal 
            meanEmgGM = mean(reshapedEmgGM.');
            meanEmgVL = mean(reshapedEmgVL.');
            meanEmgVM = mean(reshapedEmgVM.');
            
        % Wavelet transform
            [pGM, intGM, tIntGM, fGM, paGM] = emg_wavelet_act(meanEmgGM,...
                samplingFrequency, scale, control, numWavelets, pth);
            [pVL, intVL, tIntVL, fVL, paVL] = emg_wavelet_act(meanEmgVL,...
                samplingFrequency, scale, control, numWavelets, pth);   
            [pVM, intVM, tIntVM, fVM, paVM] = emg_wavelet_act(meanEmgVM,...
                samplingFrequency, scale, control, numWavelets, pth);
            
        % Validation: check that the active portion of each signal was
        % correctly isolated
            ACT = figure;
            maxfig(ACT,1);
            subplot(1,3,1)
            plot(mean(pGM'))
            subplot(1,3,2)
            plot(mean(pVL'))
            subplot(1,3,3)
            plot(mean(pVM'))
            disp('PRESS ENTER if heel strikes are correctly location')
            disp('PRESS Ctr+C in the command window if not')
            pause
            close all  
        
        % Calculates the mean total intensity of the active portion of each signal   
            meanTotIntGM(i) = mean(tIntGM);
            meanTotIntVL(i) = mean(tIntVL); 
            meanTotIntVM(i) = mean(tIntVM); 
            
        % Calculates the mean frequency of the active portion of each signal 
            mnfGM(i) = fGM*paGM.cfs'/sum(fGM);
            mnfVL(i) = fVL*paVL.cfs'/sum(fVL);
            mnfVM(i) = fVM*paVM.cfs'/sum(fVM);
    
    end  
    
    % Reorganizes the values in a table (ROW = session, COLUMN = trial)
        meanTotIntGM = reshape(meanTotIntGM,numSessions,2).';    
        meanTotIntVL = reshape(meanTotIntVL,numSessions,2).';        
        meanTotIntVM = reshape(meanTotIntVM,numSessions,2).';
        mnfGM = reshape(mnfGM,numSessions,2).';    
        mnfVL = reshape(mnfVL,numSessions,2).';        
        mnfVM = reshape(mnfVM,numSessions,2).';
        
    % Calculates the means between trials
        meanTotIntGMMean = mean(meanTotIntGM);
        meanTotIntVLMean = mean(meanTotIntVL);        
        meanTotIntVMMean = mean(meanTotIntVM);
        mnfGMMean = mean(mnfGM);
        mnfVLMean = mean(mnfVL);        
        mnfVMMean = mean(mnfVM);
    
    % Export Results
        walkResults = [meanTotIntGM meanTotIntVL meanTotIntVM;mnfGM mnfVL mnfVM;...
            meanTotIntGMMean meanTotIntVLMean meanTotIntVMMean;mnfGMMean mnfVLMean mnfVMMean];
        csvwrite([subjectID '_Walk.csv'],walkResults);
     
        

%     % Process the EMG signals recorded during the walking trials with a
%     % wavelet filter and then transform with a wavelet function to be
%     % represented in time-frequency space.
%     
%     % Calculate the power, intensity and total intensity for every trial,
%     % then plots the evolution of the total intensity over time (5
%     % conditions, 2 trials each = 10 plots), and calculate the mean, max
%     % and min total intensity over a 60secs period
%     
%     for i=1:size(emg_name,1)
%         
%         % Load data
%         data = read_data(pathName,emg_name(i,:),4);
%         dataAcc = data(:,1);
%         emgDataGM = data(:,2);
%         emgDataVL = data(:,3);
%         emgDataVM = data(:,4);
%         
%         % Normalize data around 0
%         dataAcc = dataAcc - mean(dataAcc);
%         emgDataGM = emgDataGM - mean(emgDataGM);
%         emgDataVL = emgDataVL - mean(emgDataVL);
%         emgDataVM = emgDataVM - mean(emgDataVM);
%         
%         % Define the time scale according to the frequency rate
%         time = 0+1/samplingFrequency:1/samplingFrequency:length(data)/samplingFrequency-1/samplingFrequency;        
%         
%         % Calculate the number of steps recorded (= number of peaks
%         % recorded by the accelerometer) and their corresponding location
%         [numSteps, hsLoc] = heel_strike_loc(dataAcc, samplingFrequency, 500, 100);
%         
%         % Crop each signal to fit the filtered accelerometer signal
%         reshapedDataAcc = reshape_signal_step(dataAcc, numSteps, hsLoc, w1, w2, samplingFrequency);
%         reshapedEmgDataGM = reshape_signal_step(emgDataGM, numSteps, hsLoc, w1, w2, samplingFrequency);
%         reshapedEmgDataVL = reshape_signal_step(emgDataVL, numSteps, hsLoc, w1, w2, samplingFrequency);
%         reshapedEmgDataVM = reshape_signal_step(emgDataVM, numSteps, hsLoc, w1, w2, samplingFrequency);
%         
%         % Calculate the mean step for each signal
%         meanDataAcc = mean(reshapedDataAcc.');
%         meanEmgDataGM = mean(reshapedEmgDataGM.');
%         meanEmgDataVL = mean(reshapedEmgDataVL.');
%         meanEmgDataVM = mean(reshapedEmgDataVM.');
% 
%         % Use Vinzenz wavelet transform to calculate power, intensity,
%         % total intensity and frequency spectrum of each mean EMG signal
%         [pGM, intGM, tIntGM, fGM, paGM] = emg_wavelet(meanEmgDataGM, samplingFrequency, samplingFrequency, scale, control, numWavelets);
%         [pVL, intVL, tIntVL, fVL, paVL] = emg_wavelet(meanEmgDataVL, samplingFrequency, samplingFrequency, scale, control, numWavelets);
%         [pVM, intVM, tIntVM, fVM, paVM] = emg_wavelet(meanEmgDataVM, samplingFrequency, samplingFrequency, scale, control, numWavelets);
%         maxF(i) = max([fGM fVL fVM]);
%         
%         % Plots the intensity calculated for a mean step
%         
%             % GM
%             figIntGMWalk = figure(3);
%             subplot(2,numSessions,i)
%                 contourf(intGM,20,'Linestyle','None')
%                 % Setup X ticks according to time (ms)
%                 xticks([1 (w1/2)*samplingFrequency/1000 (w1)*samplingFrequency/1000 ...
%                     (w1+w2/2)*samplingFrequency/1000 (w1+w2)*samplingFrequency/1000])
%                 xticklabels({num2str(-w1),num2str(-w1/2),'0',num2str(w2/2),num2str(w2)})
%                 % Setup Y ticks according to time (ms)
%                 yticks([1 2 3 4 5 6 7 8 9 10 11 12 13])
%                 yticklabels({num2str(round(paGM.cfs(1))),num2str(round(paGM.cfs(2))),num2str(round(paGM.cfs(3))),...
%                     num2str(round(paGM.cfs(4))),num2str(round(paGM.cfs(5))),num2str(round(paGM.cfs(6))),...
%                     num2str(round(paGM.cfs(7))),num2str(round(paGM.cfs(8))),num2str(round(paGM.cfs(9))),...
%                     num2str(round(paGM.cfs(10))),num2str(round(paGM.cfs(11))),num2str(round(paGM.cfs(12))),...
%                     num2str(round(paGM.cfs(13)))})
%                 % Removes Y axis label
%                 set(gca,'ytick',[])
%                 
%             % VL
%             figIntVLWalk = figure(4);
%             subplot(2,numSessions,i)
%                 contourf(intVL,20,'Linestyle','None')
%                 % Setup X ticks according to time (ms)
%                 xticks([1 (w1/2)*samplingFrequency/1000 (w1)*samplingFrequency/1000 ...
%                     (w1+w2/2)*samplingFrequency/1000 (w1+w2)*samplingFrequency/1000])
%                 xticklabels({num2str(-w1),num2str(-w1/2),'0',num2str(w2/2),num2str(w2)})
%                 % Removes Y axis label
%                 set(gca,'ytick',[])
%                 
%             % VM
%             figIntVMWalk = figure(5);
%             subplot(2,numSessions,i)
%                 contourf(intVM,20,'Linestyle','None')
%                 % Setup X ticks according to time (ms)
%                 xticks([1 (w1/2)*samplingFrequency/1000 (w1)*samplingFrequency/1000 ...
%                     (w1+w2/2)*samplingFrequency/1000 (w1+w2)*samplingFrequency/1000])
%                 xticklabels({num2str(-w1),num2str(-w1/2),'0',num2str(w2/2),num2str(w2)})
%                 % Removes Y axis label
%                 set(gca,'ytick',[])
%         
%         % Clears re-used variables
%         clear newDataAcc
%         clear newEmgDataGM
%         clear newEmgDataVL
%         clear newEmgDataVM
%         
%     end
%     
%     % Defines figure caption
%     set(figIntGMWalk,'numbertitle','off','name','Mean intensity of the GM for each step during walking trials');
%         fileNameGMWalk = strcat('Mean_Intensity_Walking_GM');
%     set(figIntVLWalk,'numbertitle','off','name','Mean intensity of the VL for each step during walking trials');
%         fileNameVLWalk = strcat('Mean_Intensity_Walking_VL');
%     set(figIntVMWalk,'numbertitle','off','name','Mean intensity of the VM for each step during walking trials');
%         fileNameVMWalk = strcat('Mean_Intensity_Walking_VM');
%         
%     % Maximizes figures
%     maxfig(figIntGMWalk,1);
%     maxfig(figIntVLWalk,1);
%     maxfig(figIntVMWalk,1);
%         
%     % Save figure as Matlab figure
%     savefig(figIntGMWalk,fileNameGMWalk);
%     savefig(figIntVLWalk,fileNameVLWalk);
%     savefig(figIntVMWalk,fileNameVMWalk);
%     % Save figure as PNG
%     saveas(figIntGMWalk,'Mean_Intensity_Walking_GM.png');
%     saveas(figIntVLWalk,'Mean_Intensity_Walking_VL.png');
%     saveas(figIntVMWalk,'Mean_Intensity_Walking_VM.png');
%     
%     
%     % FREQUENCY SPECTRUMS of GM, VL and GM during WALKING TRIALS
% 
%     % Process the EMG signals recorded during the walking trials with a
%     % wavelet filter and then transformed with a wavelet function to be
%     % represented in time-frequency space.
%     
%     figFreqSpectrum = figure;
%     for i=1:size(emg_name,1)
%         data = read_data(pathName,emg_name(i,:),4);
%         % Define the frequency scale according to the frequency rate and
%         % number of wavelets
%         frequency = [1:13]; % Temporary
%         emgDataGM = data(:,2);
%         emgDataVL = data(:,3);
%         emgDataVM = data(:,4);
%         
%         % Normalizes data around 0
%         dataAcc = dataAcc - mean(dataAcc);
%         emgDataGM = emgDataGM - mean(emgDataGM);
%         emgDataVL = emgDataVL - mean(emgDataVL);
%         emgDataVM = emgDataVM - mean(emgDataVM);
%         
%         % Defines the time scale according to the frequency rate
%         time = 0+1/samplingFrequency:1/samplingFrequency:length(data)/samplingFrequency-1/samplingFrequency;        
%         
%         % Calculates the number of steps recorded (= number of peaks
%         % recorded by the accelerometer) and their corresponding location
%         [numSteps, hsLoc] = heel_strike_loc(dataAcc, samplingFrequency, 500, 100);
%         
%         % Crops each signal to fit the filtered accelerometer signal
%         reshapedDataAcc = reshape_signal_step(dataAcc, numSteps, hsLoc, w1, w2, samplingFrequency);
%         reshapedEmgDataGM = reshape_signal_step(emgDataGM, numSteps, hsLoc, w1, w2, samplingFrequency);
%         reshapedEmgDataVL = reshape_signal_step(emgDataVL, numSteps, hsLoc, w1, w2, samplingFrequency);
%         reshapedEmgDataVM = reshape_signal_step(emgDataVM, numSteps, hsLoc, w1, w2, samplingFrequency);
%         
%         % Calculates the mean step for each signal
%         meanDataAcc = mean(reshapedDataAcc.');
%         meanEmgDataGM = mean(reshapedEmgDataGM.');
%         meanEmgDataVL = mean(reshapedEmgDataVL.');
%         meanEmgDataVM = mean(reshapedEmgDataVM.');
% 
%         % Uses Vinzenz wavelet transform to calculate power, intensity,
%         % total intensity and frequency spectrum of each mean EMG signal
%         [pGM, intGM, tIntGM, fGM, paGM] = emg_wavelet(meanEmgDataGM, samplingFrequency, samplingFrequency, scale, control, numWavelets);
%         [pVL, intVL, tIntVL, fVL, paVL] = emg_wavelet(meanEmgDataVL, samplingFrequency, samplingFrequency, scale, control, numWavelets);
%         [pVM, intVM, tIntVM, fVM, paVM] = emg_wavelet(meanEmgDataVM, samplingFrequency, samplingFrequency, scale, control, numWavelets);
%         
%         % Calculates the Mean Power Frequency
% 
%         
%         % Plots the frequency spectrum of each muscle for each trial and
%         % condition
%             subplot(2,numSessions,i)
%             plot(frequency,fGM,'r',frequency,fVL,'g',frequency,fVM,'b')
%             ylim([0, max(maxF*1.10)]);
%             % Setup X ticks according to time (ms)
%                 xticks([1 2 3 4 5 6 7 8 9 10 11 12 13])
%                 xticklabels({num2str(round(paGM.cfs(1))),num2str(round(paGM.cfs(2))),num2str(round(paGM.cfs(3))),...
%                     num2str(round(paGM.cfs(4))),num2str(round(paGM.cfs(5))),num2str(round(paGM.cfs(6))),...
%                     num2str(round(paGM.cfs(7))),num2str(round(paGM.cfs(8))),num2str(round(paGM.cfs(9))),...
%                     num2str(round(paGM.cfs(10))),num2str(round(paGM.cfs(11))),num2str(round(paGM.cfs(12))),...
%                     num2str(round(paGM.cfs(13)))})
%                 % Removes Y axis label
%                 set(gca,'ytick',[])
% 
%     end
%     set(figFreqSpectrum,'numbertitle','off','name','Frequency spectrum of GM, VL and VM during walking trials');
%         fileNameFreqSpectrum = strcat('Frequency_Spectrum');
%         savefig(figFreqSpectrum,fileNameFreqSpectrum);
%         
    disp('     *****     DONE WITH WALKING TRIALS     *****     ')
        
end
        
    
%% PEDAR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Execution
if pedar == 0
else

    % Plots the pressure distrubution measured using Pedar during the
    % Control trial (left) and one of the other conditions (middle) and
    % plots the difference between the 2 plots (right)

    % Trial 1
    for i=2:numSessions
        map = pedarDiffPlot(cellSheetPath,areaSheetPath,insoleType,[pathName pedar_name(1,:)],[pathName pedar_name(i,:)]);
        set(map,'numbertitle','off','name',sprintf('Trial 1: Crl VS Cond %d',i-1));
        fileName = strcat(subjectID,'_T1__Ctrl_VS_Cond_');
        % Maximizes figures
        maxfig(map,1);
        % Save figure as PNG
        saveas(map,sprintf([fileName '%d.png'],i-1));
    end
    
    % Trial 2
    for i=7:numSessions+5
        map = pedarDiffPlot(cellSheetPath,areaSheetPath,insoleType,[pathName pedar_name(6,:)],[pathName pedar_name(i,:)]);
        set(map,'numbertitle','off','name',sprintf('Trial 2: Ctrl VS Cond %d',i-6));
        fileName = strcat(subjectID,'_T2_Ctrl_VS_Cond_');
        % Maximizes figures
        maxfig(map,1);
        % Save figure as PNG
        saveas(map,sprintf([fileName '%d.png'],i-6));
    end
    
    disp('     *****     DONE WITH PEDAR TRIALS     *****     ')
    
end

%% END

finish = strcat('_______________DONE WITH SUBJECT__',subjectID,'_______________');
disp(finish)

