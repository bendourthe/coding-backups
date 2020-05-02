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
%   Last update: June 6th, 2018
%
%% Input: 
%   - Currex data (static videos, dynamic means)
%   - Biodex data (50% trials)
%   - EMG data (walking trials)
%   - Pedar data (standing trials)
%
%% Output:
%   - AI_Stat: Static arch indexes
%   - AI_Dyn: Dynamic arch indexes
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

    cd 'D:\Dr. Scholls\Pilot\Protocol V3\Ben_incline\10\';
    pathName = 'D:\Dr. Scholls\Pilot\Protocol V3\Ben_incline\10';

%% Subject ID   

    subjectID = '10';
    
%% Execution  

    % Determine which parts of the code will be run (0: won't be run; 1: will be run)
    currexS = 1;
    currexD = 1;
    half = 1;
    walking = 0;
    pedar = 0;
    
%% Settings
    
    % Sessions
    numSessions = 1;
    startSession = 1;
    endSession = 1;
    session = num2str([startSession:endSession].','%01d');
    
    % Filtering parameters
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
    
    % Recording time during walking (all recordings together)
    recTime = 40;               % in mins
    
    % Pedar Insole
    cellSheetPath = 'C:\Users\bdour\OneDrive\Work\Calgary\Documentation\Pedar\Centroid Calculation_All Insoles_B.xlsx';
    areaSheetPath = 'C:\Users\bdour\OneDrive\Work\Calgary\Documentation\Pedar\Pedar_Insole_area_B.xlsx';
    insoleType = 'XS Insole';
    
    
%% Filenames

    % Currex
        % Static (videos)
        cs1 = [strcat('\',session,'_1_',subjectID,'_csv_1.csv')];
        cs2 = [strcat('\',session,'_1_',subjectID,'_csv_2.csv')];
        cs3 = [strcat('\',session,'_1_',subjectID,'_csv_3.csv')];
        cs4 = [strcat('\',session,'_2_',subjectID,'_csv_1.csv')];
        cs5 = [strcat('\',session,'_2_',subjectID,'_csv_2.csv')];
        cs6 = [strcat('\',session,'_2_',subjectID,'_csv_3.csv')];
        currexStat_name = [cs1; cs2; cs3; cs4; cs5; cs6];
        % Dynamic (mean)
        cd1 = [strcat('\',session,'_1_',subjectID,'_cdm_l.csv')];
        cd2 = [strcat('\',session,'_1_',subjectID,'_cdm_r.csv')];
        cd3 = [strcat('\',session,'_2_',subjectID,'_cdm_l.csv')];
        cd4 = [strcat('\',session,'_2_',subjectID,'_cdm_r.csv')];
        currexDyn_name = [cd1; cd2; cd3; cd4];
    
    % Biodex 50% Trials
        half1 = [strcat('\',session,'_1_',subjectID,'_50p_1.csv')];
        half2 = [strcat('\',session,'_1_',subjectID,'_50p_2.csv')];
        half3 = [strcat('\',session,'_2_',subjectID,'_50p_1.csv')];
        half_names = [half1; half2; half3];

    % Walking
    if walking == 1
        % Pedar
        ped1 = [strcat('\',session,'_',subjectID,'_ped_1.asc')];
        ped2 = [strcat('\',session,'_',subjectID,'_ped_2.asc')];
        pedar_name = [ped1; ped2];
        % EMG
        emg1 = [strcat('\',session,'_',subjectID,'_emg_1.emg')];
        emg2 = [strcat('\',session,'_',subjectID,'_emg_2.emg')];
        emg3 = [strcat('\',session,'_',subjectID,'_emg_3.emg')];
        emg4 = [strcat('\',session,'_',subjectID,'_emg_4.emg')];
        emg5 = [strcat('\',session,'_',subjectID,'_emg_5.emg')];
        emg6 = [strcat('\',session,'_',subjectID,'_emg_6.emg')];
        emg7 = [strcat('\',session,'_',subjectID,'_emg_7.emg')];
        emg8 = [strcat('\',session,'_',subjectID,'_emg_8.emg')];
            % Data reconstruction
            for i=1:numSessions
                data1 = read_data(pathName,emg1(i,:),4);
                data2 = read_data(pathName,emg2(i,:),4);
                data3 = read_data(pathName,emg3(i,:),4);
                data4 = read_data(pathName,emg4(i,:),4);
                data5 = read_data(pathName,emg5(i,:),4);
                data6 = read_data(pathName,emg6(i,:),4);
                data7 = read_data(pathName,emg7(i,:),4);
                data8 = read_data(pathName,emg8(i,:),4);
                % Condition 1
                if i==1
                    emg_1 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:);...
                        data5(2:end,:);data6(2:end,:);data7(2:end,:);data8(2:end,:)];
                % Condition 2
                elseif i==2
                    emg_2 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:);...
                        data5(2:end,:);data6(2:end,:);data7(2:end,:);data8(2:end,:)];              
                % Condition 3
                elseif i==3
                    emg_3 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:);...
                        data5(2:end,:);data6(2:end,:);data7(2:end,:);data8(2:end,:)];
                % Condition 4
                elseif i==4
                    emg_4 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:);...
                        data5(2:end,:);data6(2:end,:);data7(2:end,:);data8(2:end,:)];
                % Condition 5
                elseif i==5
                    emg_5 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:);...
                        data5(2:end,:);data6(2:end,:);data7(2:end,:);data8(2:end,:)];
                end      
                clear data1 data2 data3 data4 data5 data6 data7 data8
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
        xlswrite([subjectID '_AI_Stat.xlsx'],ArchIndexStat);

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
        xlswrite([subjectID '_AI_Dyn.xlsx'],ArchIndexDyn);

    disp('     *****     DONE WITH CURREX DYNAMIC TRIALS     *****     ')

end
    
%% 50% Trials

    % Execution
    if half == 0
    else

    % Processing   
    for i=1:size(half_names,1)        
        % Load data
        data = csvread([pathName half_names(i,:)],2,1);
        dataTorque = abs(data(:,1))*144.72+2.3368;
        emgGL = data(:,2);
        emgGM = data(:,3);
        % Normalize data around 0
        emgGL = emgGL - mean(emgGL); 
        emgGM = emgGM - mean(emgGM); 
        % Calculate power, intensity, and mean frequency using FFT
        [p50pGL, i50pGL, meanTI50pGL(i), f50pGL, MNF50pGL(i)] = emg_fft_act_fix(emgGL,...
                    samplingFrequency, pStartGL, pEndGL);
        [p50pGM, i50pGM, meanTI50pGM(i), f50pGM, MNF50pGM(i)] = emg_fft_act_fix(emgGM,...
                    samplingFrequency, pStartGM, pEndML);
        % Mean torque  
        meanT50p(i) = mean(dataTorque);    

    end  
    
    % Reorganize the values in a table (ROW = trial, COLUMN = session)
    meanT50p = reshape(meanT50p,numSessions,3).';    
    meanTI50pGL = reshape(meanTI50pGL,numSessions,3).';  
    meanTI50pGM = reshape(meanTI50pGM,numSessions,3).';
    MNF50pGL = reshape(MNF50pGL,numSessions,3).';
    MNF50pGM = reshape(MNF50pGM,numSessions,3).';
        
    % Means between trials
        % Baseline
        meanT50pMean1 = mean(meanT50p(1:2,:));
        meanTI50pGLMean1 = mean(meanTI50pGL(1:2,:));  
        meanTI50pGMMean1 = mean(meanTI50pGM(1:2,:));
        MNF50pGLMean1 = mean(MNF50pGL(1:2,:));
        MNF50pGMMean1 = mean(MNF50pGM(1:2,:));
        % Post-Intervention
        meanT50pMean2 = meanT50p(3,:);
        meanTI50pGLMean2 = meanTI50pGL(3,:); 
        meanTI50pGMMean2 = meanTI50pGM(3,:);
        MNF50pGLMean2 = MNF50pGL(3,:);
        MNF50pGMMean2 = MNF50pGM(3,:);
    
    % Mean increase/decrease
    for i=1:numSessions
        meanT50pMeanDiff(i) = meanT50pMean2(i)*100/meanT50pMean1(i) - 100;
        meanTI50pGLMeanDiff(i) = meanTI50pGLMean2(i)*100/meanTI50pGLMean1(i) - 100;
        meanTI50pGMMeanDiff(i) = meanTI50pGMMean2(i)*100/meanTI50pGMMean1(i) - 100;
        MNF50pGLMeanDiff(i) = MNF50pGLMean2(i)*100/MNF50pGLMean1(i) - 100;
        MNF50pGMMeanDiff(i) = MNF50pGMMean2(i)*100/MNF50pGMMean1(i) - 100;
    end
    
    % Export Results
    halfResults = [meanT50p meanTI50pGL meanTI50pGM MNF50pGL MNF50pGM;...
        meanT50pMean1 meanTI50pGLMean1 meanTI50pGMMean1 MNF50pGLMean1 MNF50pGMMean1;...
        meanT50pMean2 meanTI50pGLMean2 meanTI50pGMMean2 MNF50pGLMean2 MNF50pGMMean2;...
        meanT50pMeanDiff meanTI50pGLMeanDiff meanTI50pGMMeanDiff MNF50pGLMeanDiff MNF50pGMMeanDiff];
    xlswrite([subjectID '_50p.xlsx'],halfResults);
    
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
            iteration_number = [i j]
            
            % Heel strike detection
            [numSteps, hsLoc] = heel_strike_loc_n(accMin(j,:), samplingFrequency, 500, 100);
            
            % Validation: check if steps are correctly detected
%             disp('PRESS ENTER if heel strikes are correctly location')
%             disp('PRESS Ctr+C in the command window if not')
%             pause
            close all            
            
            % Reshape signals (ROW = data point; COLUMN = step)
            reshapedEmgGM = reshape_signal_step(emgGMMin(j,:), numSteps, hsLoc, w1, w2, samplingFrequency);
            reshapedEmgVL = reshape_signal_step(emgVLMin(j,:), numSteps, hsLoc, w1, w2, samplingFrequency);
            reshapedEmgVM = reshape_signal_step(emgVMMin(j,:), numSteps, hsLoc, w1, w2, samplingFrequency);
            
            % Wavelet transform
            for k=1:numSteps
                [pGM, intGM, tIntGM, fGM, paGM] = emg_wavelet_act(reshapedEmgGM(:,k),...
                    samplingFrequency, scale, control, numWavelets, pth);
                [pVL, intVL, tIntVL, fVL, paVL] = emg_wavelet_act(reshapedEmgVL(:,k),...
                    samplingFrequency, scale, control, numWavelets, pth);   
                [pVM, intVM, tIntVM, fVM, paVM] = emg_wavelet_act(reshapedEmgVM(:,k),...
                    samplingFrequency, scale, control, numWavelets, pth);
            
                % Validation: check that the active portion of each signal was
                % correctly isolated
%                 ACT = figure;
%                 maxfig(ACT,1);
%                 subplot(1,3,1)
%                 plot(mean(pGM'))
%                 subplot(1,3,2)
%                 plot(mean(pVL'))
%                 subplot(1,3,3)
%                 plot(mean(pVM'))
%                 disp('PRESS ENTER if heel strikes are correctly location')
%                 disp('PRESS Ctr+C in the command window if not')
%                 pause
%                 close all
                
                % Calculates the mean frequency of the active portion of each step 
                stepMnfGM(k) = fGM*paGM.cfs'/sum(fGM);
                stepMnfVL(k) = fVL*paVL.cfs'/sum(fVL);
                stepMnfVM(k) = fVM*paVM.cfs'/sum(fVM);
                
            end
            
            % Calculates the mean of all the steps of the corresponding frame 
            mnfGM(i,j) = mean(stepMnfGM);
            mnfVL(i,j) = mean(stepMnfVL);
            mnfVM(i,j) = mean(stepMnfVM);
            
        end
        
        mnf = figure;
        plot(mnfGM(i,:),'r')
        hold on
        plot(mnfVL(i,:),'g')
        plot(mnfVM(i,:),'b')
        ylim([0 300])
        set(mnf,'numbertitle','off','name',...
            sprintf('Mean Frequency evolution for the GM (red), VL (greeen) and VM (VM) during Trial %d',i));
        
        % Maximize figures
        maxfig(mnf,1);
        
        % Save figure as PNG
        saveas(mnf,sprintf([subjectID '_Walk_%d.png'],i));
                
    end
    
        % Reorganizes the values in a table (ROW = trial, COLUMN = session)
        mnfGM = reshape(mnfGM,numSessions,recTime).';    
        mnfVL = reshape(mnfVL,numSessions,recTime).';        
        mnfVM = reshape(mnfVM,numSessions,recTime).';
        
        % Calculates the means between trials
        mnfGMMean = mean(mnfGM);
        mnfVLMean = mean(mnfVL);        
        mnfVMMean = mean(mnfVM);
        
        % Export Results
        walkResults = [mnfGMMean mnfVLMean mnfVMMean; mnfGM mnfVL mnfVM];
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
