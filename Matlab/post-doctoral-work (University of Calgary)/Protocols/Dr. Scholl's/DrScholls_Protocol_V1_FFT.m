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
%   Last update: June 14th, 2018
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
%       - heel_strike_loc_n.m
%       - reshape_signal_step
%       - fft_real_vvt.m
%       - emg_fft.m
%       - pedarDiffPlot.m
%       - heel_strike_loc.m


%% ID and DIRECTORY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Subject ID
    subjectID = '101';
    
% Current directory (where the data are stored for one participant)
cd 'D:\Dr. Scholls\Protocol V1\101\';
pathName = 'D:\Dr. Scholls\Protocol V1\101';


%% PROTOCOL SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Sessions
    numSessions = 5;
    session = num2str([1:numSessions].','%01d');
    
    % Code execution: determines which parts of the code will be run (0:
    % won't be run; 1: will be run)
    currexS = 0;
    currexD = 0;
    mvc = 0;
    half = 1;
    walking = 1;
    pedar = 0;
    
    % Signal parameters
    samplingFrequency = 2400;    % frenquence used during protocol (in Hz)
    
    % Window definition
    w1 = 300;   % how much time (ms) before heel strike
    w2 = 600;   % how much time (ms) after heel strike
    
    % Activity window
    pStartGL = 55;                  % percentage of the step when muscle is expected to start being active
    pEndGL = 95;                    % percentage of the step when muscle is expected to stop being active
    pStartGM = pStartGL;            % percentage of the step when muscle is expected to start being active
    pEndGM = pEndGL;                % percentage of the step when muscle is expected to stop being active
    pStartVL = 30;                  % percentage of the step when muscle is expected to start being active
    pEndVL = 50;                    % percentage of the step when muscle is expected to stop being active
    pStartVM = pStartVL;            % percentage of the step when muscle is expected to start being active
    pEndVM = pEndVL;                % percentage of the step when muscle is expected to stop being active
    
    % Validation
    valStepDetection = 0;           % 0: WON'T plot visual feedback NOR ask for validation; 1: WILL plot visual feedback AND ask for validation
    valPowerDetection = 0;          % 0: WON'T plot visual feedback NOR ask for validation; 1: WILL plot visual feedback AND ask for validation
    
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
            dataEMG = data(12000:24000,2);     
        % FFT
            [p50p, int50p, meanTotInt50p(i), f50p, MNF50p(i)] = emg_fft(dataEMG,samplingFrequency);        
        % Calculates the mean torque
            meanTorque50p(i) = mean(dataTorque);    
         
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

    % Processing
    
    for i=1:size(emg_name,1)
        emg_name(i,:)
        % Load data
            data = read_data(pathName,emg_name(i,:),4);
            acc = data(:,1);
            emgGM = data(:,2);
            emgVL = data(:,3);
            emgVM = data(:,4);
            
        %% Heel strike detection
            [numSteps, hsLoc] = heel_strike_loc_n(acc, samplingFrequency, 500, 100);
            
            % Validation: check if steps are correctly detected
            if valStepDetection == 1
                disp('PRESS ENTER if heel strikes are correctly location')
                disp('PRESS Ctr+C in the command window if not')
                pause
            else
            end
            close all            
            
            % Reshape signals (ROW = data point; COLUMN = step)
            reshapedEmgGM = reshape_signal_step(emgGM, numSteps, hsLoc, w1, w2, samplingFrequency);
            reshapedEmgVL = reshape_signal_step(emgVL, numSteps, hsLoc, w1, w2, samplingFrequency);
            reshapedEmgVM = reshape_signal_step(emgVM, numSteps, hsLoc, w1, w2, samplingFrequency);
            
            for k=1:numSteps
                % FFT
                [pGM, intGM, stepMtiGM(k), fGM, stepMnfGM(k)] = emg_fft_act_fix(reshapedEmgGM(:,k),...
                    samplingFrequency, pStartGM, pEndGM);
                [pVL, intVL, stepMtiVL(k), fVL, stepMnfVL(k)] = emg_fft_act_fix(reshapedEmgVL(:,k),...
                    samplingFrequency, pStartVL, pEndVL);   
                [pVM, intVM, stepMtiVM(k), fVM, stepMnfVM(k)] = emg_fft_act_fix(reshapedEmgVM(:,k),...
                    samplingFrequency, pStartVM, pEndVM);            
                    % Validation: check that the active portion of each signal was
                    % correctly isolated
                    if valPowerDetection == 1
                        ACT = figure;
                        maxfig(ACT,1);
                        subplot(1,4,1)
                        plot(mean(pGL'))
                        subplot(1,4,2)
                        plot(mean(pGM'))
                        subplot(1,4,3)
                        plot(mean(pVL'))
                        subplot(1,4,4)
                        plot(mean(pVM'))
                        disp('PRESS ENTER if heel strikes are correctly location')
                        disp('PRESS Ctr+C in the command window if not')
                        pause
                        close all
                    else
                    end                
                
            end
            
            % Mean of all the steps of the corresponding frame
                % Mean frequency
                mnfGM(i) = mean(stepMnfGM);
                mnfVL(i) = mean(stepMnfVL);
                mnfVM(i) = mean(stepMnfVM);
                % Mean total intensity
                mtiGM(i) = mean(stepMtiGM);
                mtiVL(i) = mean(stepMtiVL);
                mtiVM(i) = mean(stepMtiVM);
            
            % Standard deviation of all the steps of the corresponding frame
                % Mean frequency
                mnfGM_sd(i) = std(stepMnfGM);
                mnfVL_sd(i) = std(stepMnfVL);
                mnfVM_sd(i) = std(stepMnfVM);
                % Mean total intensity
                mtiGM_sd(i) = std(stepMtiGM);
                mtiVL_sd(i) = std(stepMtiVL);
                mtiVM_sd(i) = std(stepMtiVM);
                                        
    end
    
        % Reorganize the values in a table (ROW = trial, COLUMN = session)
        mnfGM = reshape(mnfGM,numSessions,2).';   
        mnfVL = reshape(mnfVL,numSessions,2).';       
        mnfVM = reshape(mnfVM,numSessions,2).';
        mtiGM = reshape(mtiGM,numSessions,2).';    
        mtiVL = reshape(mtiVL,numSessions,2).';        
        mtiVM = reshape(mtiVM,numSessions,2).';
        
        % Means between trials
        mnfGMMean = mean(mnfGM);
        mnfVLMean = mean(mnfVL);        
        mnfVMMean = mean(mnfVM);
        mtiGMMean = mean(mtiGM);
        mtiVLMean = mean(mtiVL);        
        mtiVMMean = mean(mtiVM);
        
        % Export Results
        walkResults = [mnfGM mnfVL mnfVM mtiGM mtiVL mtiVM;mnfGMMean mnfVLMean mnfVMMean mtiGMMean mtiVLMean mtiVMMean];
        csvwrite([subjectID '_Walk.csv'],walkResults);
        
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

