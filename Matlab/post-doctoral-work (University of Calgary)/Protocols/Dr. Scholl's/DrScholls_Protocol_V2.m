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
%   Last update: June 4th, 2018
%
%% Input: 
%   - Currex data (static videos, dynamic means)
%   - Biodex data (MVC trials, 50% trials)
%   - EMG data (walking trials)
%   - Pedar data (standing trials)
%
%% Output:
%   - AI_Stat: Static arch indexes
%   - AI_Dyn: Dynamic arch indexes
%   - MVC: Max torques, raw VAR and VAR
%   - 50p: Mean Torque, mean GM EMG total intensity, mean GM EMG frequency
% Note: each file has the same structure:
%   1. Baseline trials (first few rows)
%   2. Post-intervention trials (next few rows)
%   3. Baseline mean (3rd last row)
%   4. Post-intervention mean (2nd last row)
%   5. Change from Baseline to Post-intervention (in % where baseline
%   represents 100%)
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


%% Initialization

    clear ; close all; clc
    
%% Directory

    cd 'D:\Dr. Scholls\Phase 2\03\';
    pathName = 'D:\Dr. Scholls\Phase 2\03';

%% Subject ID   

    subjectID = '03';
    
%% Execution  

    % Determine which parts of the code will be run (0: won't be run; 1: will be run)
    currexS = 0;
    currexD = 0;
    mvc = 0;
    half = 0;
    walking = 1;
    pedar = 0;
    
%% Settings
    
    % Sessions
    numSessions = 1;            % how many sessions to process?
    startSession = 3;           % what is the first session number to process?
    endSession = 3;             % what is the last session number to process?
    sessions = num2str([startSession:endSession].','%01d');
    
    % Recording time during walking trials (all recordings together)
    recTime = 40;               % in mins
    
    % Number of channels during walking trials
    chanNum = 6;
    
    % Filtering parameters
    samplingFrequency = 2400;   % sampling frequency used during protocol (Hz)
    scale = 1.5;                % determines the frequency range of the wavelets
    control = 3;
    numWavelets = 20;           % determines the number of wavelets (more wavelets = more frequency resolution, but more processing time)
    
    % Heel strike window (to define which portion of each step we look at)
    w1 = 300;                   % how much time (ms) before heel strike
    w2 = 600;                   % how much time (ms) after heel strike
    
    % Activity window (when each muscle is expected to be active)
    pStartGL = 55;              % percentage of the step when muscle is expected to start being active
    pEndGL = 95;                % percentage of the step when muscle is expected to stop being active
    pStartGM = pStartGL;        % percentage of the step when muscle is expected to start being active
    pEndGM = pEndGL;            % percentage of the step when muscle is expected to stop being active
    pStartVL = 25;              % percentage of the step when muscle is expected to start being active
    pEndVL = 50;                % percentage of the step when muscle is expected to stop being active
    pStartVM = pStartVL;        % percentage of the step when muscle is expected to start being active
    pEndVM = pEndVL;            % percentage of the step when muscle is expected to stop being active
    
    % Validation
    valStepDetection = 0;       % 0: WON'T plot visual feedback NOR ask for validation; 1: WILL plot visual feedback AND ask for validation
    valPowerDetection = 0;      % 0: WON'T plot visual feedback NOR ask for validation; 1: WILL plot visual feedback AND ask for validation
    
    % Pedar insole documents (needed to determine pedar coordinate system)
    cellSheetPath = 'C:\Users\bdour\OneDrive\Work\Calgary\Documentation\Pedar\Centroid Calculation_All Insoles_B.xlsx';
    areaSheetPath = 'C:\Users\bdour\OneDrive\Work\Calgary\Documentation\Pedar\Pedar_Insole_area_B.xlsx';
    insoleType = 'XS Insole';
    
    % Figures properties
    set(groot,'defaultFigureColor','w')
    set(groot,'defaultLineLineWidth',1);
    set(groot,'defaultAxesLineWidth',2);
    set(groot,'defaultAxesTickDir','out')
    set(groot,'defaultAxesFontSize',15)
    set(groot,'defaultAxesFontWeight','bold')    
    
%% Filenames

    % Currex
        % Static (videos)
        cs1 = [strcat('\',sessions,'_1_',subjectID,'_csv_1.csv')];
        cs2 = [strcat('\',sessions,'_1_',subjectID,'_csv_2.csv')];
        cs3 = [strcat('\',sessions,'_1_',subjectID,'_csv_3.csv')];
        cs4 = [strcat('\',sessions,'_2_',subjectID,'_csv_1.csv')];
        cs5 = [strcat('\',sessions,'_2_',subjectID,'_csv_2.csv')];
        cs6 = [strcat('\',sessions,'_2_',subjectID,'_csv_3.csv')];
        currexStat_name = [cs1; cs2; cs3; cs4; cs5; cs6];
        % Dynamic (mean)
        cd1 = [strcat('\',sessions,'_1_',subjectID,'_cdm_l.csv')];
        cd2 = [strcat('\',sessions,'_1_',subjectID,'_cdm_r.csv')];
        cd3 = [strcat('\',sessions,'_2_',subjectID,'_cdm_l.csv')];
        cd4 = [strcat('\',sessions,'_2_',subjectID,'_cdm_r.csv')];
        currexDyn_name = [cd1; cd2; cd3; cd4];
    
    % Biodex
        % MVC Trials
        mvc1 = [strcat('\',sessions,'_1_',subjectID,'_mvc_1.csv')];
        mvc2 = [strcat('\',sessions,'_1_',subjectID,'_mvc_2.csv')];
        mvc3 = [strcat('\',sessions,'_2_',subjectID,'_mvc_1.csv')];
        mvc_name = [mvc1; mvc2; mvc3];
        % 50p Trials
        half1 = [strcat('\',sessions,'_1_',subjectID,'_50p_1.csv')];
        half2 = [strcat('\',sessions,'_1_',subjectID,'_50p_2.csv')];
        half3 = [strcat('\',sessions,'_2_',subjectID,'_50p_1.csv')];
        half_names = [half1; half2; half3];

    % Walking
    if walking == 1
        % Pedar
        ped1 = [strcat('\',sessions,'_',subjectID,'_ped_1.asc')];
        ped2 = [strcat('\',sessions,'_',subjectID,'_ped_2.asc')];
        pedar_name = [ped1; ped2];
        % EMG
            % Data reconstruction
            for i=1:numSessions
                % Condition 1
                if i==1
                    [emg_1] = import_emg(pathName,subjectID,numSessions,sessions,recTime,chanNum);
                % Condition 2
                elseif i==2
                    [emg_1, emg_2] = import_emg(pathName,subjectID,numSessions,sessions,recTime,chanNum);              
                % Condition 3
                elseif i==3
                    [emg_1, emg_2, emg_3] = import_emg(pathName,subjectID,numSessions,sessions,recTime,chanNum);
                % Condition 4
                elseif i==4
                    [emg_1, emg_2, emg_3, emg_4] = import_emg(pathName,subjectID,numSessions,sessions,recTime,chanNum);
                % Condition 5
                elseif i==5
                    [emg_1, emg_2, emg_3, emg_4, emg_5] = import_emg(pathName,subjectID,numSessions,sessions,recTime,chanNum);
                end      
            end
    else
    end

%% Analysis

%% Currex Static

    % Execution
    if currexS == 0
    else

        % Arch Index
        for i=1:size(currexStat_name,1)
            current_trial = i
            % Load data
            data = csvread([pathName currexStat_name(i,:)]);  
            % Calculation
            [AILeftStat(i),AIRightStat(i)] = archindexStat_function(data);
        end

        % Reorganize the values in a table (ROW = trial, COLUMN = session)
        AILeftStat = reshape(AILeftStat,numSessions,6).';
        AIRightStat = reshape(AIRightStat,numSessions,6).';

        % Mean Arch Index between trials
            % Baseline
            AILeftMeanStat1 = mean(AILeftStat(1:3,:));
            AIRightMeanStat1 = mean(AIRightStat(1:3,:));
            % Post-Intervention
            AILeftMeanStat2 = mean(AILeftStat(4:6,:));
            AIRightMeanStat2 = mean(AIRightStat(4:6,:));
        
        % Mean Arch Index increase/decrease
        for i=1:numSessions
            AILeftDiffStat(i) = AILeftMeanStat2(i)*100/AILeftMeanStat1(i) - 100;
            AIRightDiffStat(i) = AIRightMeanStat2(i)*100/AIRightMeanStat1(i) - 100;
        end

        % Exports Results
        ArchIndexStat = [AILeftStat AIRightStat; AILeftMeanStat1 AIRightMeanStat1; ...
            AILeftMeanStat2 AIRightMeanStat2; AILeftDiffStat AIRightDiffStat];
        csvwrite([subjectID '_AI_Stat.csv'],ArchIndexStat);

    disp('     *****     DONE WITH CURREX STATIC TRIALS     *****     ')

    end

    %% Currex Dynamic

    % Execution
    if currexD == 0
    else

        % Arch Index
        for i=1:size(currexDyn_name,1)
            current_trial = i
            % Load data
            data = csvread([pathName currexDyn_name(i,:)]);  
            % Calculation
            AIDyn(i) = archindexDyn_function(data);
        end

        % Reorganizes the values in a table (ROW = trial, COLUMN = session)
        AIDyn = reshape(AIDyn,numSessions,4).';
    
        % Isoles left from right foot (left = rows 1, 3; right = rows 2, 4)
        % (note: rows 1 and 2 = baseline; rows 3 and 4 = Post-Intervention)
        AILeftDyn = [AIDyn(1,:); AIDyn(3,:)];
        AIRightDyn = [AIDyn(2,:); AIDyn(4,:)];

        % Mean Arch Index increase/decrease
        for i=1:numSessions
            AILeftDiffDyn(i) = AILeftDyn(2,i)*100/AILeftDyn(1,i) - 100;
            AIRightDiffDyn(i) = AIRightDyn(2,i)*100/AIRightDyn(1,i) - 100;
        end

        % Exports Results
        ArchIndexDyn = [AILeftDyn AIRightDyn; AILeftDiffDyn AIRightDiffDyn];
        csvwrite([subjectID '_AI_Dyn.csv'],ArchIndexDyn);

    disp('     *****     DONE WITH CURREX DYNAMIC TRIALS     *****     ')

end
    
%% MVC Trials

    % Execution
    if mvc == 0
    else

    for i=1:size(mvc_name,1)
        current_trial = i
        
        % Load data
        data = csvread([pathName mvc_name(i,:)],2,1);
        dataTorque = abs(data(:,1))*144.72+2.3368;
        dataEMG = data(:,2);

        % Peak torque (PT)
        mvcPT(i) = max(dataTorque);

        % Percentage of voluntary activation (PVA)
        PVA = figure;
        maxfig(PVA,1);
        plot(dataTorque)
        hold on
        plot(dataEMG)
        disp('Manually select:')
        disp('     1) Torque right before twitch during MVC')
        disp('     2) Peak torque resulting from twitch during MVC')
        disp('     3) Mean torque when foot is inactive after MVC')
        disp('     4) Peak torque resulting from twitch after MVC')
        mvcP = ginput(4);
        mvcP = round(mvcP(:,1));
        mvcT = [min(dataTorque(mvcP(1)-100:mvcP(1)+100));...
            max(dataTorque(mvcP(2)-100:mvcP(2)+100));...
            min(dataTorque(mvcP(3)-200:mvcP(3)+200));...
            max(dataTorque(mvcP(4)-200:mvcP(4)+200))];
        if mvcT(2) < mvcPT
            mvcRatio(i) = 100 - 100*(mvcT(2) - mvcT(1))/(mvcT(4) - mvcT(3));
            mvcVAR(i) = (mvcPT(i)*mvcRatio(i))/mvcT(2);
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
            % Plot torque data with the selected point for verification
            figure
            clf
            hold on
            plot(dataTorque)
            plot(mvcP,mvcT,'r*','MarkerSize',10)
            disp('PRESS ENTER if points are correctly selected - PRESS Ctr+C in the command window if not')
            pause

            close all

    end

    % Reorganize the values in a table (ROW = trial, COLUMN = session)
    mvcPT = reshape(mvcPT,numSessions,3).';
    mvcRatio = reshape(mvcRatio,numSessions,3).';
    mvcVAR = reshape(mvcVAR,numSessions,3).';

    % Means between trials
        % Baseline
        mvcPTMean1 = mean(mvcPT(1:2,:));
        mvcRatioMean1 = mean(mvcRatio(1:2,:));
        mvcVARMean1 = mean(mvcVAR(1:2,:));
        % Post-Intervention
        mvcPTMean2 = mvcPT(3,:);
        mvcRatioMean2 = mvcRatio(3,:);
        mvcVARMean2 = mvcVAR(3,:);
    
    % Mean increase/decrease
    for i=1:numSessions
        mvcPTMeanDiff(i) = mvcPTMean2(i)*100/mvcPTMean1(i) - 100;
        mvcRatioMeanDiff(i) = mvcRatioMean2(i)*100/mvcRatioMean1(i) - 100;
        mvcVARMeanDiff(i) = mvcVARMean2(i)*100/mvcVARMean1(i) - 100;
    end

    % Export Results
    mvcResults = [mvcPT mvcRatio mvcVAR  ; mvcPTMean1 mvcRatioMean1 mvcVARMean1;...
        mvcPTMean2 mvcRatioMean2 mvcVARMean2; mvcPTMeanDiff mvcRatioMeanDiff mvcVARMeanDiff];
    csvwrite([subjectID '_MVC.csv'],mvcResults);

    disp('     *****     DONE WITH MVC TRIALS     *****     ')

end
    
    
%% 50% Trials

    % Execution
    if half == 0
    else

    % Wavelet filtering    
    for i=1:size(half_names,1)        
                % Load data
        data = csvread([pathName half_names(i,:)],2,1);
        dataTorque = abs(data(:,1))*144.72+2.3368;
        dataEMG = data(:,2);
        % Normalize data around 0
        dataEMG = dataEMG - mean(dataEMG);        
        % Define the time scale according to the frequency rate
        time = 0:1/samplingFrequency:length(data)/samplingFrequency-1/samplingFrequency;
        time = time*1000;   % in ms                
        % Apply Vinzenz Wavelet transform
        [p50p, i50p, tI50p, f50p, pa50p] = emg_wavelet(dataEMG,...
            samplingFrequency, scale, control, numWavelets);        
        % Mean torque and total intensity  
        meanT50p(i) = mean(dataTorque);    
        meanTI50p(i) = mean(tI50p);            
        % Mean frequency
        MNF50p(i) = f50p*pa50p.cfs'/sum(f50p);         
    end  
    
    % Reorganize the values in a table (ROW = trial, COLUMN = session)
    meanT50p = reshape(meanT50p,numSessions,3).';    
    meanTI50p = reshape(meanTI50p,numSessions,3).';        
    MNF50p = reshape(MNF50p,numSessions,3).';
        
    % Means between trials
        % Baseline
        meanT50pMean1 = mean(meanT50p(1:2,:));
        meanTI50pMean1 = mean(meanTI50p(1:2,:));        
        MNF50pMean1 = mean(MNF50p(1:2,:));
        % Post-Intervention
        meanT50pMean2 = meanT50p(3,:);
        meanTI50pMean2 = meanTI50p(3,:);        
        MNF50pMean2 = MNF50p(3,:);
    
    % Mean increase/decrease
    for i=1:numSessions
        meanT50pMeanDiff(i) = meanT50pMean2(i)*100/meanT50pMean1(i) - 100;
        meanTI50pMeanDiff(i) = meanTI50pMean2(i)*100/meanTI50pMean1(i) - 100;
        MNF50pMeanDiff(i) = MNF50pMean2(i)*100/MNF50pMean1(i) - 100;
    end
    
    % Export Results
    halfResults = [meanT50p meanTI50p MNF50p; meanT50pMean1 meanTI50pMean1 MNF50pMean1;...
        meanT50pMean2 meanTI50pMean2 MNF50pMean2; meanT50pMeanDiff meanTI50pMeanDiff MNF50pMeanDiff];
    csvwrite([subjectID '_50p.csv'],halfResults);
    
    disp('     *****     DONE WITH 50% TRIALS     *****     ')

end

%% Walking Trials
    
    % Execution
    if walking == 0
    else

    % Processing
    for i=1:numSessions
        
        % Load data
            % Condition 1
            if i==1             
                data = emg_1;
            % Condition 2
            elseif i==2         
                data = emg_2;
            % Condition 3
            elseif i==3         
                data = emg_3;
            % Condition 4
            elseif i==4         
                data = emg_4;
            % Condition 5
            elseif i==5         
                data = emg_5;
            end
            
        % Identify data
        acc = data(:,1);
        emgGM = data(:,2);
        emgVL = data(:,3);
        emgVM = data(:,4);
        
        % Reshape data (every ROW = 1 min)
        accMin = reshape(acc,samplingFrequency*60,recTime)';
        emgGMMin = reshape(emgGM,samplingFrequency*60,recTime)';
        emgVLMin = reshape(emgVL,samplingFrequency*60,recTime)';
        emgVMMin = reshape(emgVM,samplingFrequency*60,recTime)';
      
        for j=1:recTime
            
            % Heel strike detection
            [numSteps, hsLoc] = heel_strike_loc_n(accMin(j,:), samplingFrequency, 500, 100);
            randStep = round(rand(1)*numSteps);         % Select one random step to plot for validation
            
            % Validation: check if steps are correctly detected
            if valStepDetection == 1
                disp('PRESS ENTER if heel strikes are correctly location')
                disp('PRESS Ctr+C in the command window if not')
                pause
            else
            end
            close all            
            
            % Reshape signals (ROW = data point; COLUMN = step)
            reshapedEmgGM = reshape_signal_step(emgGMMin(j,:), numSteps, hsLoc, w1, w2, samplingFrequency);
            reshapedEmgVL = reshape_signal_step(emgVLMin(j,:), numSteps, hsLoc, w1, w2, samplingFrequency);
            reshapedEmgVM = reshape_signal_step(emgVMMin(j,:), numSteps, hsLoc, w1, w2, samplingFrequency);
            
            for k=1:numSteps
                disp('    Step |Out of |Frame |Condition')
                disp([k  numSteps  j  i])
                % Wavelet transform
                [pGM, intGM, tIntGM, fGM, paGM] = emg_wavelet_act_fix(reshapedEmgGM(:,k),...
                    samplingFrequency, scale, control, numWavelets, pStartGM, pEndGM);
                [pVL, intVL, tIntVL, fVL, paVL] = emg_wavelet_act_fix(reshapedEmgVL(:,k),...
                    samplingFrequency, scale, control, numWavelets, pStartVL, pEndVL);   
                [pVM, intVM, tIntVM, fVM, paVM] = emg_wavelet_act_fix(reshapedEmgVM(:,k),...
                    samplingFrequency, scale, control, numWavelets, pStartVM, pEndVM);            
                    % Validation: check that the active portion of each signal was
                    % correctly isolated
                    if valPowerDetection == 1                        
                        if k==randStep                          
                            test = figure;
                            maxfig(test,1);
                            subplot(1,3,1)
                            plot(reshapedEmgGM(:,k))
                            hold on
                            plot(round(pStartGM/100*length(reshapedEmgGM(:,k))):round(pEndGM/100*length(reshapedEmgGM(:,k))),...
                                reshapedEmgGM(round(pStartGM/100*length(reshapedEmgGM(:,k))):round(pEndGM/100*length(reshapedEmgGM(:,k))),k),'r')
                            subplot(1,3,2)
                            plot(reshapedEmgVL(:,k))
                            hold on
                            plot(round(pStartVL/100*length(reshapedEmgVL(:,k))):round(pEndVL/100*length(reshapedEmgVL(:,k))),...
                                reshapedEmgVL(round(pStartVL/100*length(reshapedEmgVL(:,k))):round(pEndVL/100*length(reshapedEmgVL(:,k))),k),'r')
                            subplot(1,3,3)
                            plot(reshapedEmgVM(:,k))
                            hold on
                            plot(round(pStartVM/100*length(reshapedEmgVM(:,k))):round(pEndVM/100*length(reshapedEmgVM(:,k))),...
                                reshapedEmgVM(round(pStartVM/100*length(reshapedEmgVM(:,k))):round(pEndVM/100*length(reshapedEmgVM(:,k))),k),'r')
                            set(test,'numbertitle','off','name',...
                                sprintf('Validation plot: is the active portion of each EMG signal (in red) correctly identified by the fixed window? STEP #%d',k));
        
                            disp('PRESS ENTER if active portion of each EMG signal (in red) is correctly identified by the fixed window')
                            disp('PRESS Ctr+C in the command window if not')
                            pause
                            close all
                        else
                        end
                    else
                    end                
                    % Mean frequency of the active portion of each step
                    stepMnfGM(k) = fGM*paGM.cfs'/sum(fGM);
                    stepMnfVL(k) = fVL*paVL.cfs'/sum(fVL);
                    stepMnfVM(k) = fVM*paVM.cfs'/sum(fVM);
                    
                    % Mean total intensity of the active portion of each step
                    stepMtiGM(k) = mean(tIntGM);
                    stepMtiVL(k) = mean(tIntVL);
                    stepMtiVM(k) = mean(tIntVM);
                
            end
            
            % Mean of all the steps of the corresponding frame
                % Mean frequency
                mnfGM(i,j) = mean(stepMnfGM);
                mnfVL(i,j) = mean(stepMnfVL);
                mnfVM(i,j) = mean(stepMnfVM);
                % Mean total intensity
                mtiGM(i,j) = mean(stepMtiGM);
                mtiVL(i,j) = mean(stepMtiVL);
                mtiVM(i,j) = mean(stepMtiVM);
            
            % Standard deviation of all the steps of the corresponding frame
                % Mean frequency
                mnfGM_sd(i,j) = std(stepMnfGM);
                mnfVL_sd(i,j) = std(stepMnfVL);
                mnfVM_sd(i,j) = std(stepMnfVM);
                % Mean total intensity
                mtiGM_sd(i,j) = std(stepMtiGM);
                mtiVL_sd(i,j) = std(stepMtiVL);
                mtiVM_sd(i,j) = std(stepMtiVM);
                
        end
        
        walkFig = figure;
        
        subplot(2,2,1)
        hold on
        plot(mnfGM(i,:),'y')
        plot(mnfVL(i,:),'b')
        plot(mnfVM(i,:),'c')
            
        subplot(2,2,2)
        hold on
        plot(mnfGM_sd(i,:),'y')
        plot(mnfVL_sd(i,:),'b')
        plot(mnfVM_sd(i,:),'c')
       
        subplot(2,2,3)
        hold on
        plot(mtiGM(i,:),'y')
        plot(mtiVL(i,:),'b')
        plot(mtiVM(i,:),'c')
        
        subplot(2,2,4)
        hold on
        plot(mtiGM_sd(i,:),'y')
        plot(mtiVL_sd(i,:),'b')
        plot(mtiVM_sd(i,:),'c')
        
        set(walkFig,'numbertitle','off','name',...
            sprintf('Walking trial: top left (MNF), top right (MNFvar), bottom left (MTI), bottom right (MTIvar) - GM (yellow), VL (blue) and VM (cyan) during Condition %d',i));
        
        % Maximize figures
        maxfig(walkFig,1);
        
        % Save figure as JPG
        saveas(walkFig,sprintf([subjectID '_Walk_%d.png'],i));
        
    end
    
    % Reorganize the values in a table (ROW = trial, COLUMN = session)
    mnfGM = reshape(mnfGM,numSessions,recTime).';    
    mnfVL = reshape(mnfVL,numSessions,recTime).';        
    mnfVM = reshape(mnfVM,numSessions,recTime).';
    mtiGM = reshape(mtiGM,numSessions,recTime).';    
    mtiVL = reshape(mtiVL,numSessions,recTime).';        
    mtiVM = reshape(mtiVM,numSessions,recTime).';

    % Means between trials
    mnfGMMean = mean(mnfGM);
    mnfVLMean = mean(mnfVL);        
    mnfVMMean = mean(mnfVM);
    mtiGMMean = mean(mtiGM);
    mtiVLMean = mean(mtiVL);        
    mtiVMMean = mean(mtiVM);

    % Export Results
    walkResults = [mnfGMMean mnfVLMean mnfVMMean mtiGMMean mtiVLMean mtiVMMean; ...
        mnfGM mnfVL mnfVM mtiGM mtiVL mtiVM];
    csvwrite([subjectID '_Walk.csv'],walkResults);
        
    disp('     *****     DONE WITH WALKING TRIALS     *****     ')
        
end
        
    
%% Pedar

    % Execution
    if pedar == 0
    else

    % Plot the pressure distrubution measured using Pedar during the
    % Control trial (left) and one of the other conditions (middle) and
    % plots the difference between the 2 plots (right)

    % Trial 1
    for i=2:4
        map = pedarDiffPlot(cellSheetPath,areaSheetPath,insoleType,[pathName pedar_name(1,:)],[pathName pedar_name(i,:)]);
        set(map,'numbertitle','off','name',sprintf('Trial 1: Control VS Condition %d',i-1));
        fileName = strcat('Trial_1_',subjectID,'_Control_VS_Condition_');
        % Maximizes figures
        maxfig(map,1);
        % Save figure as Matlab figure
        savefig(map,sprintf([fileName '%d.fig'],i-1));
        % Save figure as PNG
        saveas(map,sprintf([fileName '%d.png'],i-1));
    end
    
    % Trial 2
    for i=7:9
        map = pedarDiffPlot(cellSheetPath,areaSheetPath,insoleType,[pathName pedar_name(6,:)],[pathName pedar_name(i,:)]);
        set(map,'numbertitle','off','name',sprintf('Trial 2: Control VS Condition %d',i-6));
        fileName = strcat('Trial_2_',subjectID,'_Control_VS_Condition_');
        % Maximizes figures
        maxfig(map,1);
        % Save figure as Matlab figure
        savefig(map,sprintf([fileName '%d.fig'],i-6));
        % Save figure as PNG
        saveas(map,sprintf([fileName '%d.png'],i-6));
    end
    
    disp('     *****     DONE WITH PEDAR TRIALS     *****     ')
    
end
