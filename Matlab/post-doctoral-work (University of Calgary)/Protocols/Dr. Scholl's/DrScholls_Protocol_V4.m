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
%   Last update: October 10th, 2018
%
%% Input: 
%   - Currex data (static videos, dynamic means)
%   - Biodex data (EMG GL and GM during 50% trials)
%   - EMG data (GL, GM, VL, VM during walking trials)
%   - Pedar data (standing trials)
%
%% Output:
%   - AI_Stat.csv: Static arch index data
%   - AI_Dyn.csv: Dynamic arch index data
%   - 50p.csv: 50% isometric contraction data
%       The rows of each file are as followed:
%           1. Absolute change from Baseline to Post-intervention
%           2. Relative change from Baseline to Post-intervention
%           3. Baseline results
%           4. Post-intervention results
%       The columns of each file are as followed:
%           Each series of five columns represent the results for the
%           different conditions in this order: CTRL, CF, BPR, MG, REP
%           (zeros are added when the corresponding condition wasn't
%           tested)
%   - Walk1, Walk2 and Walk3.csv: Processed EMG data from the intervention
%       The rows of each file are as followed:
%           1. Absolute change from First to Last minute of the intervention
%           2. Mean for the whole intervention
%           3 to end. Mean of each minute recorded during the intervention
%       The columns of each file are as followed:
%           1 to 4. EMG frequency data for the GL, GM, VL and VM
%           5 to 8. EMG activity data for the GL, GM, VL and VM
%   - Freq1, Freq2 and Freq3.csv: Detailed power spectrum of each muscle
%       The rows of each file are as followed:
%           1. Central frequencies of the wavelets
%           2. Detailed power spectrum of each muscle (GL, GM, VL, VM) for 
%           each minute frame (one line = one power spectrum)
%   - pedar_CTRL, pedar_CF, pedar_C3: pedar difference plots comparing
%           first vs. last minute among each condition
%   - pedar_First_CTRL_vs_CF, pedar_Last_CTRL_vs_CF: pedar difference plots
%           comparing first vs. first and last vs. last minutes of the CTRL
%           and CF conditions
%
%% Dependencies:
%
%   Files:
%       - Centroid Calculation_All Insoles_B.xlsx
%       - Pedar_Insole_area_B.xlsx
%
%   Functions:
%       - import_emg.m
%       - read_data.m
%       - archindexStat_function.m
%       - archindexDyn_function.m
%       - emg_wavelet.m
%       - wavelet_transform_V70419.m
%       - heel_strike_loc_n.m
%       - reshape_signal_step.m
%       - emg_wavelet_act_fix.m
%       - pedarDiffPlot.m


%% Initialization

    clear ; close all; clc
    
%% Execution  

    % Determine which parts of the code will be run (0: won't be run; 1: will be run)
    currexS = 0;
    currexD = 0;
    half = 0;
    walking = 0;
    pedar = 1;
    
%% ID and Directory

    subjectID = '028';
    pathName = ['C:\Users\bdour\OneDrive\Work\Calgary\Projects\Dr. Scholls\Pedar\' subjectID];
    cd(pathName);
    
%% Settings
    
    % Sessions
    numSessions = 2;            % how many sessions to process?
    startSession = 1;           % what is the first session number to process?
    endSession = 3;             % what is the last session number to process?
    sessions = num2str([startSession:endSession].','%01d'); 
    condition3 = 3;             % define third session: type 3 for bpr, 4 for mg, 5 for repeated
    insole = 'R';               % define pedar insole used (type B for black and R for red)
    
    % Filtering parameters
    samplingFrequency = 2400;   % sampling frequency used during protocol (Hz)
    scale = 1.5;                % determines the frequency range of the wavelets
    control = 3;
    numWavelets = 20;           % determines the number of wavelets (more wavelets = more frequency resolution, but more processing time)
   
    % Walking parameters
    if walking == 1 
        % Recording time during walking trials (all recordings together)
        recTime = 60;               % in mins

        % Number of channels during walking trials
        chanNum = 5;

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
    else
    end
    
    % Pedar insole documents (needed to determine pedar coordinate system)
    if pedar == 1
        cellSheetPath = 'C:\Users\bdour\OneDrive\Work\Calgary\Documentation\Pedar\Centroid Calculation_All Insoles_B.xlsx';
        areaSheetPath = 'C:\Users\bdour\OneDrive\Work\Calgary\Documentation\Pedar\Pedar_Insole_area_B.xlsx';
        if insole == 'B'
            insoleType = 'XS Insole';
        elseif insole == 'R'
            insoleType = 'YS Insole';
        end
    else
    end
    
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
        if currexS == 1
            cs1 = [strcat('\',sessions,'_1_',subjectID,'_csv_1.csv')];
            cs2 = [strcat('\',sessions,'_1_',subjectID,'_csv_2.csv')];
            cs3 = [strcat('\',sessions,'_1_',subjectID,'_csv_3.csv')];
            cs4 = [strcat('\',sessions,'_2_',subjectID,'_csv_1.csv')];
            cs5 = [strcat('\',sessions,'_2_',subjectID,'_csv_2.csv')];
            cs6 = [strcat('\',sessions,'_2_',subjectID,'_csv_3.csv')];
            currexStat_name = [cs1; cs2; cs3; cs4; cs5; cs6];
        else
        end
        % Dynamic (mean)
        if currexD == 1
            cd1 = [strcat('\',sessions,'_1_',subjectID,'_cdm_l.csv')];
            cd2 = [strcat('\',sessions,'_1_',subjectID,'_cdm_r.csv')];
            cd3 = [strcat('\',sessions,'_2_',subjectID,'_cdm_l.csv')];
            cd4 = [strcat('\',sessions,'_2_',subjectID,'_cdm_r.csv')];
            currexDyn_name = [cd1; cd2; cd3; cd4];
        else
        end
    
    % Biodex 50% Trials
    if half == 1
        half1 = [strcat('\',sessions,'_1_',subjectID,'_50p_1.csv')];
        half2 = [strcat('\',sessions,'_1_',subjectID,'_50p_2.csv')];
        half3 = [strcat('\',sessions,'_2_',subjectID,'_50p_1.csv')];
        half_names = [half1; half2; half3];
    else
    end

    % Walking
    if walking == 1
        % Data reconstruction
        for i = 1:numSessions
            % Condition 1
            if i == 1
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
    
    % Pedar
    if pedar == 1
        ped1 = [strcat('\',sessions,'_',subjectID,'_pedSD_1.asc')];
        ped2 = [strcat('\',sessions,'_',subjectID,'_pedSD_2.asc')];
        pedar_name = [ped1; ped2];
    else
    end

%% Analysis

%% Currex Static

    % Execution
    if currexS == 0
    else

        % Arch Index
        for i = 1:size(currexStat_name,1)
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
        
        % Absolute difference from Baseline to Post-Intervention
        for i = 1:numSessions
            AILeftDiffStatAbs(i) = AILeftMeanStat2(i) - AILeftMeanStat1(i);
            AIRightDiffStatAbs(i) = AIRightMeanStat2(i) - AIRightMeanStat1(i);
        end
        
        % Relative difference from Baseline to Post-Intervention (i.e.
        % Baseline = 100%)
        for i = 1:numSessions
            AILeftDiffStatRel(i) = (AILeftMeanStat2(i) - AILeftMeanStat1(i))*100 / AILeftMeanStat1(i);
            AIRightDiffStatRel(i) = (AIRightMeanStat2(i) - AIRightMeanStat1(i))*100 / AIRightMeanStat1(i);
        end
        
        % Attribute 0 value to the conditions that are not being tested
        % (to match Excel format)
        if condition3 == 3
            AILeftDiffStatAbs = [AILeftDiffStatAbs 0 0];
            AIRightDiffStatAbs = [AIRightDiffStatAbs 0 0];
            AILeftDiffStatRel = [AILeftDiffStatRel 0 0];
            AIRightDiffStatRel = [AIRightDiffStatRel 0 0];
            AILeftMeanStat1 = [AILeftMeanStat1 0 0];
            AIRightMeanStat1 = [AIRightMeanStat1 0 0];
            AILeftMeanStat2 = [AILeftMeanStat2 0 0];
            AIRightMeanStat2 = [AIRightMeanStat2 0 0];
        elseif condition3 == 4
            AILeftDiffStatAbs = [AILeftDiffStatAbs(1) AILeftDiffStatAbs(2) 0 AILeftDiffStatAbs(3) 0];
            AIRightDiffStatAbs = [AIRightDiffStatAbs(1) AIRightDiffStatAbs(2) 0 AIRightDiffStatAbs(3) 0];
            AILeftDiffStatRel = [AILeftDiffStatRel(1) AILeftDiffStatRel(2) 0 AILeftDiffStatRel(3) 0];
            AIRightDiffStatRel = [AIRightDiffStatRel(1) AIRightDiffStatRel(2) 0 AIRightDiffStatRel(3) 0];
            AILeftMeanStat1 = [AILeftMeanStat1(1) AILeftMeanStat1(2) 0 AILeftMeanStat1(3) 0];
            AIRightMeanStat1 = [AIRightMeanStat1(1) AIRightMeanStat1(2) 0 AIRightMeanStat1(3) 0];
            AILeftMeanStat2 = [AILeftMeanStat2(1) AILeftMeanStat2(2) 0 AILeftMeanStat2(3) 0];
            AIRightMeanStat2 = [AIRightMeanStat2(1) AIRightMeanStat2(2) 0 AIRightMeanStat2(3) 0];
        elseif condition3 == 5
            AILeftDiffStatAbs = [AILeftDiffStatAbs(1) AILeftDiffStatAbs(2) 0 0 AILeftDiffStatAbs(3)];
            AIRightDiffStatAbs = [AIRightDiffStatAbs(1) AIRightDiffStatAbs(2) 0 0 AIRightDiffStatAbs(3)];
            AILeftDiffStatRel = [AILeftDiffStatRel(1) AILeftDiffStatRel(2) 0 0 AILeftDiffStatRel(3)];
            AIRightDiffStatRel = [AIRightDiffStatRel(1) AIRightDiffStatRel(2) 0 0 AIRightDiffStatRel(3)];
            AILeftMeanStat1 = [AILeftMeanStat1(1) AILeftMeanStat1(2) 0 0 AILeftMeanStat1(3)];
            AIRightMeanStat1 = [AIRightMeanStat1(1) AIRightMeanStat1(2) 0 0 AIRightMeanStat1(3)];
            AILeftMeanStat2 = [AILeftMeanStat2(1) AILeftMeanStat2(2) 0 0 AILeftMeanStat2(3)];
            AIRightMeanStat2 = [AIRightMeanStat2(1) AIRightMeanStat2(2) 0 0 AIRightMeanStat2(3)];
        end
        
        % Exports Results
        ArchIndexStat = [AILeftDiffStatAbs AIRightDiffStatAbs; AILeftDiffStatRel AIRightDiffStatRel; ...
            AILeftMeanStat1 AIRightMeanStat1; AILeftMeanStat2 AIRightMeanStat2];
        xlswrite([subjectID '_AI_Stat.xlsx'],ArchIndexStat);

    disp('     *****     DONE WITH CURREX STATIC DATA ANALYSIS     *****     ')

    end

%% Currex Dynamic

    % Execution
    if currexD == 0
    else

        % Arch Index
        for i = 1:size(currexDyn_name,1)
            current_trial = i
            % Load data
            data = csvread([pathName currexDyn_name(i,:)]);  
            % Calculation
            AIDyn(i) = archindexDyn_function(data);
        end

        % Reorganize the values in a table (ROW = trial, COLUMN = session)
        AIDyn = reshape(AIDyn,numSessions,4).';
    
        % Isolate left from right foot (left = rows 1, 3; right = rows 2, 4)
        % (note: rows 1 and 2 = baseline; rows 3 and 4 = Post-Intervention)
        AILeftDyn = [AIDyn(1,:); AIDyn(3,:)];
        AIRightDyn = [AIDyn(2,:); AIDyn(4,:)];

        % Absolute difference from Baseline to Post-Intervention
        for i = 1:numSessions
            AILeftDiffDynAbs(i) = AILeftDyn(2,i) - AILeftDyn(1,i);
            AIRightDiffDynAbs(i) = AIRightDyn(2,i) - AIRightDyn(1,i);
        end
        
        % Relative difference from Baseline to Post-Intervention (i.e.
        % Baseline = 100%)
        for i = 1:numSessions
            AILeftDiffDynRel(i) = (AILeftDyn(2,i) - AILeftDyn(1,i))*100 / AILeftDyn(1,i);
            AIRightDiffDynRel(i) = (AIRightDyn(2,i) - AIRightDyn(1,i))*100 / AIRightDyn(1,i);
        end
        
        % Attribute 0 value to the conditions that are not being tested
        % (to match Excel format)
        if condition3 == 3
            AILeftDiffDynAbs = [AILeftDiffDynAbs 0 0];
            AIRightDiffDynAbs = [AIRightDiffDynAbs 0 0];
            AILeftDiffDynRel = [AILeftDiffDynRel 0 0];
            AIRightDiffDynRel = [AIRightDiffDynRel 0 0];
            AILeftDyn1 = [AILeftDyn(1,:) 0 0];
            AIRightDyn1 = [AIRightDyn(1,:) 0 0];
            AILeftDyn2 = [AILeftDyn(2,:) 0 0];
            AIRightDyn2 = [AIRightDyn(2,:) 0 0];
        elseif condition3 == 4
            AILeftDiffDynAbs = [AILeftDiffDynAbs(1) AILeftDiffDynAbs(2) 0 AILeftDiffDynAbs(3) 0];
            AIRightDiffDynAbs = [AIRightDiffDynAbs(1) AIRightDiffDynAbs(2) 0 AIRightDiffDynAbs(3) 0];
            AILeftDiffDynRel = [AILeftDiffDynRel(1) AILeftDiffDynRel(2) 0 AILeftDiffDynRel(3) 0];
            AIRightDiffDynRel = [AIRightDiffDynRel(1) AIRightDiffDynRel(2) 0 AIRightDiffDynRel(3) 0];
            AILeftDyn1 = [AILeftDyn(1,1) AILeftDyn(1,2) 0 AILeftDyn(1,3) 0];
            AIRightDyn1 = [AIRightDyn(1,1) AIRightDyn(1,2) 0 AIRightDyn(1,3) 0];
            AILeftDyn2 = [AILeftDyn(2,1) AILeftDyn(2,2) 0 AILeftDyn(2,3) 0];
            AIRightDyn2 = [AIRightDyn(2,1) AIRightDyn(2,2) 0 AIRightDyn(2,3) 0];
        elseif condition3 == 5
            AILeftDiffDynAbs = [AILeftDiffDynAbs(1) AILeftDiffDynAbs(2) 0 0 AILeftDiffDynAbs(3)];
            AIRightDiffDynAbs = [AIRightDiffDynAbs(1) AIRightDiffDynAbs(2) 0 0 AIRightDiffDynAbs(3)];
            AILeftDiffDynRel = [AILeftDiffDynRel(1) AILeftDiffDynRel(2) 0 0 AILeftDiffDynRel(3)];
            AIRightDiffDynRel = [AIRightDiffDynRel(1) AIRightDiffDynRel(2) 0 0 AIRightDiffDynRel(3)];
            AILeftDyn1 = [AILeftDyn(1,1) AILeftDyn(1,2) 0 0 AILeftDyn(1,3)];
            AIRightDyn1 = [AIRightDyn(1,1) AIRightDyn(1,2) 0 0 AIRightDyn(1,3)];
            AILeftDyn2 = [AILeftDyn(2,1) AILeftDyn(2,2) 0 0 AILeftDyn(2,3)];
            AIRightDyn2 = [AIRightDyn(2,1) AIRightDyn(2,2) 0 0 AIRightDyn(2,3)];
        end

        % Exports Results
        ArchIndexDyn = [AILeftDiffDynAbs AIRightDiffDynAbs; AILeftDiffDynRel AIRightDiffDynRel; ...
            AILeftDyn1 AIRightDyn1; AILeftDyn2 AIRightDyn2];
        xlswrite([subjectID '_AI_Dyn.xlsx'],ArchIndexDyn);

    disp('     *****     DONE WITH CURREX DYNAMIC DATA ANALYSIS     *****     ')

end
    
%% 50% Trials

    % Execution
    if half == 0
    else

        % Wavelet filtering    
        for i = 1:size(half_names,1)        
            % Load data
            data = csvread([pathName half_names(i,:)],2,1);
            dataTorque = abs(data(:,1))*144.72+2.3368;
            emgGL = data(:,2);
            emgGM = data(:,3);
            % Normalize data around 0
            emgGL = emgGL - mean(emgGL); 
            emgGM = emgGM - mean(emgGM); 
            % Apply Vinzenz Wavelet transform
            [p50pGL, i50pGL, tI50pGL, f50pGL, pa50pGL] = emg_wavelet(emgGL,...
                samplingFrequency, scale, control, numWavelets);
            [p50pGM, i50pGM, tI50pGM, f50pGM, pa50pGM] = emg_wavelet(emgGM,...
                samplingFrequency, scale, control, numWavelets); 
            % Mean torque and total intensity   
            meanTI50pGL(i) = mean(tI50pGL); 
            meanTI50pGM(i) = mean(tI50pGM);
            % Mean frequency
            MNF50pGL(i) = f50pGL*pa50pGL.cfs'/sum(f50pGL);
            MNF50pGM(i) = f50pGM*pa50pGM.cfs'/sum(f50pGM);
        end  

        % Reorganize the values in a table (ROW = trial, COLUMN = session)  
        meanTI50pGL = reshape(meanTI50pGL,numSessions,3).';  
        meanTI50pGM = reshape(meanTI50pGM,numSessions,3).';
        MNF50pGL = reshape(MNF50pGL,numSessions,3).';
        MNF50pGM = reshape(MNF50pGM,numSessions,3).';

        % Means between trials
            % Baseline
            meanTI50pGLMean1 = mean(meanTI50pGL(1:2,:));  
            meanTI50pGMMean1 = mean(meanTI50pGM(1:2,:));
            MNF50pGLMean1 = mean(MNF50pGL(1:2,:));
            MNF50pGMMean1 = mean(MNF50pGM(1:2,:));
            % Post-Intervention
            meanTI50pGLMean2 = meanTI50pGL(3,:); 
            meanTI50pGMMean2 = meanTI50pGM(3,:);
            MNF50pGLMean2 = MNF50pGL(3,:);
            MNF50pGMMean2 = MNF50pGM(3,:);

        % Absolute difference from Baseline to Post-Intervention
        for i = 1:numSessions
            meanTI50pGLMeanDiffAbs(i) = meanTI50pGLMean2(i) - meanTI50pGLMean1(i);
            meanTI50pGMMeanDiffAbs(i) = meanTI50pGMMean2(i) - meanTI50pGMMean1(i);
            MNF50pGLMeanDiffAbs(i) = MNF50pGLMean2(i) - MNF50pGLMean1(i);
            MNF50pGMMeanDiffAbs(i) = MNF50pGMMean2(i) - MNF50pGMMean1(i);
        end
        
        % Relative difference from Baseline to Post-Intervention
        for i = 1:numSessions
            meanTI50pGLMeanDiffRel(i) = (meanTI50pGLMean2(i) - meanTI50pGLMean1(i))*100 / meanTI50pGLMean1(i);
            meanTI50pGMMeanDiffRel(i) = (meanTI50pGMMean2(i) - meanTI50pGMMean1(i))*100 / meanTI50pGMMean1(i);
            MNF50pGLMeanDiffRel(i) = (MNF50pGLMean2(i) - MNF50pGLMean1(i))*100 / MNF50pGLMean1(i);
            MNF50pGMMeanDiffRel(i) = (MNF50pGMMean2(i) - MNF50pGMMean1(i))*100 / MNF50pGMMean1(i);
        end
        
        % Attribute 0 value to the conditions that are not being tested
        % (to match Excel format)
        if condition3 == 3
            meanTI50pGLMeanDiffAbs = [meanTI50pGLMeanDiffAbs 0 0];
            meanTI50pGMMeanDiffAbs = [meanTI50pGMMeanDiffAbs 0 0];
            MNF50pGLMeanDiffAbs = [MNF50pGLMeanDiffAbs 0 0];
            MNF50pGMMeanDiffAbs = [MNF50pGMMeanDiffAbs 0 0];
            meanTI50pGLMeanDiffRel = [meanTI50pGLMeanDiffRel 0 0];
            meanTI50pGMMeanDiffRel = [meanTI50pGMMeanDiffRel 0 0];
            MNF50pGLMeanDiffRel = [MNF50pGLMeanDiffRel 0 0];
            MNF50pGMMeanDiffRel = [MNF50pGMMeanDiffRel 0 0];
        elseif condition3 == 4
            meanTI50pGLMeanDiffAbs = [meanTI50pGLMeanDiffAbs(1) meanTI50pGLMeanDiffAbs(2) 0 meanTI50pGLMeanDiffAbs(3) 0];
            meanTI50pGMMeanDiffAbs = [meanTI50pGMMeanDiffAbs(1) meanTI50pGMMeanDiffAbs(2) 0 meanTI50pGMMeanDiffAbs(3) 0];
            MNF50pGLMeanDiffAbs = [MNF50pGLMeanDiffAbs(1) MNF50pGLMeanDiffAbs(2) 0 MNF50pGLMeanDiffAbs(3) 0];
            MNF50pGMMeanDiffAbs = [MNF50pGMMeanDiffAbs(1) MNF50pGMMeanDiffAbs(2) 0 MNF50pGMMeanDiffAbs(3) 0];
            meanTI50pGLMeanDiffRel = [meanTI50pGLMeanDiffRel(1) meanTI50pGLMeanDiffRel(2) 0 meanTI50pGLMeanDiffRel(3) 0];
            meanTI50pGMMeanDiffRel = [meanTI50pGMMeanDiffRel(1) meanTI50pGMMeanDiffRel(2) 0 meanTI50pGMMeanDiffRel(3) 0];
            MNF50pGLMeanDiffRel = [MNF50pGLMeanDiffRel(1) MNF50pGLMeanDiffRel(2) 0 MNF50pGLMeanDiffRel(3) 0];
            MNF50pGMMeanDiffRel = [MNF50pGMMeanDiffRel(1) MNF50pGMMeanDiffRel(2) 0 MNF50pGMMeanDiffRel(3) 0];
        elseif condition3 == 5
            meanTI50pGLMeanDiffAbs = [meanTI50pGLMeanDiffAbs(1) meanTI50pGLMeanDiffAbs(2) 0 0 meanTI50pGLMeanDiffAbs(3)];
            meanTI50pGMMeanDiffAbs = [meanTI50pGMMeanDiffAbs(1) meanTI50pGMMeanDiffAbs(2) 0 0 meanTI50pGMMeanDiffAbs(3)];
            MNF50pGLMeanDiffAbs = [MNF50pGLMeanDiffAbs(1) MNF50pGLMeanDiffAbs(2) 0 0 MNF50pGLMeanDiffAbs(3)];
            MNF50pGMMeanDiffAbs = [MNF50pGMMeanDiffAbs(1) MNF50pGMMeanDiffAbs(2) 0 0 MNF50pGMMeanDiffAbs(3)];
            meanTI50pGLMeanDiffRel = [meanTI50pGLMeanDiffRel(1) meanTI50pGLMeanDiffRel(2) 0 0 meanTI50pGLMeanDiffRel(3)];
            meanTI50pGMMeanDiffRel = [meanTI50pGMMeanDiffRel(1) meanTI50pGMMeanDiffRel(2) 0 0 meanTI50pGMMeanDiffRel(3)];
            MNF50pGLMeanDiffRel = [MNF50pGLMeanDiffRel(1) MNF50pGLMeanDiffRel(2) 0 0 MNF50pGLMeanDiffRel(3)];
            MNF50pGMMeanDiffRel = [MNF50pGMMeanDiffRel(1) MNF50pGMMeanDiffRel(2) 0 0 MNF50pGMMeanDiffRel(3)];
        end

        % Export Results
        halfResults = [meanTI50pGLMeanDiffAbs meanTI50pGMMeanDiffAbs MNF50pGLMeanDiffAbs MNF50pGMMeanDiffAbs; ...
            meanTI50pGLMeanDiffRel meanTI50pGMMeanDiffRel MNF50pGLMeanDiffRel MNF50pGMMeanDiffRel];
        xlswrite([subjectID '_50p.xlsx'],halfResults);

        disp('     *****     DONE WITH 50% EMG DATA ANALYSIS     *****     ')

    end

%% Walking Trials
    
    % Execution
    if walking == 0
    else

        % Processing
        for i = 1:numSessions

            % Load data
                % Condition 1
                if i == 1             
                    data = emg_1;
                % Condition 2
                elseif i == 2         
                    data = emg_2;
                % Condition 3
                elseif i == 3         
                    data = emg_3;
                end

            % Identify data
            acc = data(:,1);
            emgGL = data(:,2);
            emgGM = data(:,3);
            emgVL = data(:,4);
            emgVM = data(:,5);

            % Reshape data (every ROW = 1 min)
            accMin = reshape(acc,samplingFrequency*60,recTime)';
            emgGLMin = reshape(emgGL,samplingFrequency*60,recTime)';
            emgGMMin = reshape(emgGM,samplingFrequency*60,recTime)';
            emgVLMin = reshape(emgVL,samplingFrequency*60,recTime)';
            emgVMMin = reshape(emgVM,samplingFrequency*60,recTime)';

            for j = 1:recTime

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
                reshapedEmgGL = reshape_signal_step(emgGLMin(j,:), numSteps, hsLoc, w1, w2, samplingFrequency);
                reshapedEmgGM = reshape_signal_step(emgGMMin(j,:), numSteps, hsLoc, w1, w2, samplingFrequency);
                reshapedEmgVL = reshape_signal_step(emgVLMin(j,:), numSteps, hsLoc, w1, w2, samplingFrequency);
                reshapedEmgVM = reshape_signal_step(emgVMMin(j,:), numSteps, hsLoc, w1, w2, samplingFrequency);

                for k = 1:numSteps
                    disp('    Step |Out of |Min |Out of| Condition')
                    disp([k  numSteps  j  recTime i])
                    % Wavelet transform
                    [pGL, intGL, tIntGL, fGL(k,:), paGL] = emg_wavelet_act_fix(reshapedEmgGL(:,k),...
                        samplingFrequency, scale, control, numWavelets, pStartGL, pEndGL);
                    [pGM, intGM, tIntGM, fGM(k,:), paGM] = emg_wavelet_act_fix(reshapedEmgGM(:,k),...
                        samplingFrequency, scale, control, numWavelets, pStartGM, pEndGM);
                    [pVL, intVL, tIntVL, fVL(k,:), paVL] = emg_wavelet_act_fix(reshapedEmgVL(:,k),...
                        samplingFrequency, scale, control, numWavelets, pStartVL, pEndVL);   
                    [pVM, intVM, tIntVM, fVM(k,:), paVM] = emg_wavelet_act_fix(reshapedEmgVM(:,k),...
                        samplingFrequency, scale, control, numWavelets, pStartVM, pEndVM);            
                        % Validation: check that the active portion of each signal was
                        % correctly isolated
                        if valPowerDetection == 1                        
                            if k == randStep                          
                                test = figure;
                                maxfig(test,1);
                                subplot(2,2,1)
                                hold on
                                plot(reshapedEmgGL(:,k))
                                plot(round(pStartGL/100*length(reshapedEmgGL(:,k))):round(pEndGL/100*length(reshapedEmgGL(:,k))),...
                                    reshapedEmgGL(round(pStartGL/100*length(reshapedEmgGL(:,k))):round(pEndGL/100*length(reshapedEmgGL(:,k))),k),'r')
                                subplot(2,2,2)
                                hold on
                                plot(reshapedEmgGM(:,k))
                                plot(round(pStartGM/100*length(reshapedEmgGM(:,k))):round(pEndGM/100*length(reshapedEmgGM(:,k))),...
                                    reshapedEmgGM(round(pStartGM/100*length(reshapedEmgGM(:,k))):round(pEndGM/100*length(reshapedEmgGM(:,k))),k),'r')
                                subplot(2,2,3)
                                hold on
                                plot(reshapedEmgVL(:,k))
                                plot(round(pStartVL/100*length(reshapedEmgVL(:,k))):round(pEndVL/100*length(reshapedEmgVL(:,k))),...
                                    reshapedEmgVL(round(pStartVL/100*length(reshapedEmgVL(:,k))):round(pEndVL/100*length(reshapedEmgVL(:,k))),k),'r')
                                subplot(2,2,4)
                                hold on
                                plot(reshapedEmgVM(:,k))
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
                        stepMnfGL(k) = fGL(k,:)*paGL.cfs'/sum(fGL(k,:));
                        stepMnfGM(k) = fGM(k,:)*paGM.cfs'/sum(fGM(k,:));
                        stepMnfVL(k) = fVL(k,:)*paVL.cfs'/sum(fVL(k,:));
                        stepMnfVM(k) = fVM(k,:)*paVM.cfs'/sum(fVM(k,:));

                        % Mean total intensity of the active portion of each step
                        stepMtiGL(k) = mean(tIntGL);
                        stepMtiGM(k) = mean(tIntGM);
                        stepMtiVL(k) = mean(tIntVL);
                        stepMtiVM(k) = mean(tIntVM);

                end

                % Mean of all the steps of the corresponding frame
                    % Mean frequency
                    mnfGL(j,1) = mean(stepMnfGL);
                    mnfGM(j,1) = mean(stepMnfGM);
                    mnfVL(j,1) = mean(stepMnfVL);
                    mnfVM(j,1) = mean(stepMnfVM);
                    % Mean total intensity
                    mtiGL(j,1) = mean(stepMtiGL);
                    mtiGM(j,1) = mean(stepMtiGM);
                    mtiVL(j,1) = mean(stepMtiVL);
                    mtiVM(j,1) = mean(stepMtiVM);
                    % Mean frequency spectrum
                    mfsGL(j,:) = mean(fGL);
                    mfsGM(j,:) = mean(fGM);
                    mfsVL(j,:) = mean(fVL);
                    mfsVM(j,:) = mean(fVM);       

                % Standard deviation of all the steps of the corresponding frame
                    % Mean frequency
                    mnfGL_sd(j,1) = std(stepMnfGL);
                    mnfGM_sd(j,1) = std(stepMnfGM);
                    mnfVL_sd(j,1) = std(stepMnfVL);
                    mnfVM_sd(j,1) = std(stepMnfVM);
                    % Mean total intensity
                    mtiGL_sd(j,1) = std(stepMtiGL);
                    mtiGM_sd(j,1) = std(stepMtiGM);
                    mtiVL_sd(j,1) = std(stepMtiVL);
                    mtiVM_sd(j,1) = std(stepMtiVM);

            end

                % Means between trials
                mnfGLMean = mean(mnfGL);
                mnfGMMean = mean(mnfGM);
                mnfVLMean = mean(mnfVL);        
                mnfVMMean = mean(mnfVM);
                mtiGLMean = mean(mtiGL);
                mtiGMMean = mean(mtiGM);
                mtiVLMean = mean(mtiVL);        
                mtiVMMean = mean(mtiVM);
                
                % Absolute difference between first and last minute of the
                % intervention
                mnfGLDiff = mnfGL(end) - mnfGL(1);
                mnfGMDiff = mnfGM(end) - mnfGM(1);
                mnfVLDiff = mnfVL(end) - mnfVL(1);
                mnfVMDiff = mnfVM(end) - mnfVM(1);
                mtiGLDiff = mtiGL(end) - mtiGL(1);
                mtiGMDiff = mtiGM(end) - mtiGM(1);
                mtiVLDiff = mtiVL(end) - mtiVL(1);
                mtiVMDiff = mtiVM(end) - mtiVM(1);

                % Export Results
                frequencyResults = [paGL.cfs paGM.cfs paVL.cfs paVM.cfs; mfsGL mfsGM mfsVL mfsVM];
                csvwrite(sprintf([subjectID '_Freq%d.csv'],i),frequencyResults);
                walkResults = [mnfGLDiff mnfGMDiff mnfVLDiff mnfVMDiff mtiGLDiff mtiGMDiff mtiVLDiff mtiVMDiff; ...
                    mnfGLMean mnfGMMean mnfVLMean mnfVMMean mtiGLMean mtiGMMean mtiVLMean mtiVMMean; ...
                    mnfGL mnfGM mnfVL mnfVM mtiGL mtiGM mtiVL mtiVM];
                csvwrite(sprintf([subjectID '_Walk%d.csv'],i),walkResults);

        end
        
        disp('     *****     DONE WITH EMG WALKING DATA ANALYSIS     *****     ')
        
    end
            
%% Pedar

    % Execution
    if pedar == 0
    else

        % Plot the pressure distrubution measured using Pedar during the
        % first (left) and last (middle) minute of the intervention and
        % plots the difference between the 2 plots (right)
        for i = 1:numSessions
            pedarDiffPlot(cellSheetPath,areaSheetPath,insoleType,[pathName pedar_name(i,:)],[pathName pedar_name(i+3,:)]);
            if i == 1
                set(gcf,'numbertitle','off','name','Pressure distribution - First vs. Last minute - CTRL');
                fileName = strcat(subjectID,'_pedar_CTRL');
            elseif i == 2
                set(gcf,'numbertitle','off','name','Pressure distribution - First vs. Last minute - CF');
                fileName = strcat(subjectID,'_pedar_CF');
            elseif i == 3
                set(gcf,'numbertitle','off','name','Pressure distribution - First vs. Last minute - C3');
                fileName = strcat(subjectID,'_pedar_C3');
            end
            % Maximizes figures
            maxfig(gcf,1);
            % Save figure as JPG
            saveas(gcf,[fileName '.jpg']);
        end
        
        % Plot the pressure distrubution measured using Pedar during the
        % first minute of the CTRL (left) and the first minute of the CF
        % (middle) and plots the difference between the 2 plots (right)
        pedarDiffPlot(cellSheetPath,areaSheetPath,insoleType,[pathName pedar_name(1,:)],[pathName pedar_name(2,:)]);
        set(gcf,'numbertitle','off','name','Pressure distribution - First vs. First minute - CTRL vs. CF');
        fileName = strcat(subjectID,'_pedar_First_CTRL_vs_CF');
        % Maximizes figures
        maxfig(gcf,1);
        % Save figure as JPG
        saveas(gcf,[fileName '.jpg']);
        
        % Plot the pressure distrubution measured using Pedar during the
        % last minute of the CTRL (left) and the last minute of the CF
        % (middle) and plots the difference between the 2 plots (right)
        pedarDiffPlot(cellSheetPath,areaSheetPath,insoleType,[pathName pedar_name(4,:)],[pathName pedar_name(5,:)]);
        set(gcf,'numbertitle','off','name','Pressure distribution - Last vs. Last minute - CTRL vs. CF');
        fileName = strcat(subjectID,'_pedar_Last_CTRL_vs_CF');
        % Maximizes figures
        maxfig(gcf,1);
        % Save figure as JPG
        saveas(gcf,[fileName '.jpg']);

        disp('     *****     DONE WITH PEDAR DATA ANALYSIS     *****     ')
    
    end
