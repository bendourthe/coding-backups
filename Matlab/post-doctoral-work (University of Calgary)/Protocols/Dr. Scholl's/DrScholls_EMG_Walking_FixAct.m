%% DESCRIPTION:
%
% This code was designed to enable the visualization and analysis of the
% EMG data provided via the Dr. Scholl's Research Project:
%
%   Project Title: Quantifying the effects of different insole
%   configurations on fatigue
%
%   Author: Benjamin Dourthe
%
%   Last update: June 7th, 2018
%
%% Input: 
%   - EMG data (walking trials)
%
%% Output:
%   - 
%
%% Dependencies:
%
%   Files:
%       - None
%
%   Functions:
%       - read_data.m
%       - use_wavelet_transform.m
%       - wavelet_transform_V70419.m
%       - heel_strike_loc.m


%% Initialization

    clear ; close all; clc
    
%% Directory

    cd 'D:\Dr. Scholls\Pilot\Protocol V3\Kevin_incline\';
    pathName = 'D:\Dr. Scholls\Pilot\Protocol V3\Kevin_incline';

%% Subject ID   

    subjectID = '00';
    
%% Settings
    
    % Sessions
    numSessions = 1;
    startSession = 1;
    endSession = 1;
    sessions = num2str([startSession:endSession].','%01d');
    
    % Recording time during walking (all recordings together)
    recTime = 25;                   % in mins
    
    % Number of channels
    chanNum = 5;
    
    % Window definition pre- and post-heel strike
    w1 = 300;                       % time before heel strike (ms)
    w2 = 600;                       % time after heel strike (ms)
    
    % Wavelet filtering parameters
    samplingFrequency = 2400;       % in Hz
    scale = 0.15;                    % defines the frequency range covered by the wavelets (e.g. scale of 0.3 covers from 0 to ~500Hz, scale of 0.15 covers from 0 to ~250Hz)
    control = 3;
    numWavelets = 20;               % defines the frequency resolution (more wavelets = more resolution)
        
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
    
%% EMG data import

    % Note: for more than one condition, change function output to 
    % [emg_1, emg_2, emg_3, etc.] (up to emg_5)

    [emg_1] = import_emg(pathName,subjectID,numSessions,sessions,recTime,chanNum);
    
%% Data processing

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
      
        for j=1:recTime
            iteration_number = [i j]
            
            % Heel strike detection
            [numSteps, hsLoc] = heel_strike_loc_n(accMin(j,:), samplingFrequency, 500, 100);
            
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
            
            for k=1:numSteps
                % Wavelet transform
                [pGL, intGL, tIntGL, fGL, paGL] = emg_wavelet_act_fix(reshapedEmgGL(:,k),...
                    samplingFrequency, scale, control, numWavelets, pStartGL, pEndGL);
                [pGM, intGM, tIntGM, fGM, paGM] = emg_wavelet_act_fix(reshapedEmgGM(:,k),...
                    samplingFrequency, scale, control, numWavelets, pStartGM, pEndGM);
                [pVL, intVL, tIntVL, fVL, paVL] = emg_wavelet_act_fix(reshapedEmgVL(:,k),...
                    samplingFrequency, scale, control, numWavelets, pStartVL, pEndVL);   
                [pVM, intVM, tIntVM, fVM, paVM] = emg_wavelet_act_fix(reshapedEmgVM(:,k),...
                    samplingFrequency, scale, control, numWavelets, pStartVM, pEndVM);            
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
                    % Mean frequency of the active portion of each step
                    stepMnfGL(k) = fGL*paGL.cfs'/sum(fGL);
                    stepMnfGM(k) = fGM*paGM.cfs'/sum(fGM);
                    stepMnfVL(k) = fVL*paVL.cfs'/sum(fVL);
                    stepMnfVM(k) = fVM*paVM.cfs'/sum(fVM);
                    
                    % Mean total intensity of the active portion of each step
                    stepMtiGL(k) = mean(tIntGL);
                    stepMtiGM(k) = mean(tIntGM);
                    stepMtiVL(k) = mean(tIntVL);
                    stepMtiVM(k) = mean(tIntVM);
                    
                % Root Mean Square (RMS)
                yGL = (reshapedEmgGL(:,k) - mean(reshapedEmgGL(:,k))).^2;
                stepRmsGL(k) = sqrt(sum(yGL));
                yGM = (reshapedEmgGM(:,k) - mean(reshapedEmgGM(:,k))).^2;
                stepRmsGM(k) = sqrt(sum(yGM));
                yVL = (reshapedEmgVL(:,k) - mean(reshapedEmgVL(:,k))).^2;
                stepRmsVL(k) = sqrt(sum(yVL));
                yVM = (reshapedEmgVM(:,k) - mean(reshapedEmgVM(:,k))).^2;
                stepRmsVM(k) = sqrt(sum(yVM));
                
            end
            
            % Mean of all the steps of the corresponding frame
                % Mean frequency
                mnfGL(i,j) = mean(stepMnfGL);
                mnfGM(i,j) = mean(stepMnfGM);
                mnfVL(i,j) = mean(stepMnfVL);
                mnfVM(i,j) = mean(stepMnfVM);
                % Mean total intensity
                mtiGL(i,j) = mean(stepMtiGL);
                mtiGM(i,j) = mean(stepMtiGM);
                mtiVL(i,j) = mean(stepMtiVL);
                mtiVM(i,j) = mean(stepMtiVM);
                % Mean RMS
                rmsGL(i,j) = mean(stepRmsGL);
                rmsGM(i,j) = mean(stepRmsGM);
                rmsVL(i,j) = mean(stepRmsVL);
                rmsVM(i,j) = mean(stepRmsVM);
            
            % Standard deviation of all the steps of the corresponding frame
                % Mean frequency
                mnfGL_sd(i,j) = std(stepMnfGL);
                mnfGM_sd(i,j) = std(stepMnfGM);
                mnfVL_sd(i,j) = std(stepMnfVL);
                mnfVM_sd(i,j) = std(stepMnfVM);
                % Mean total intensity
                mtiGL_sd(i,j) = std(stepMtiGL);
                mtiGM_sd(i,j) = std(stepMtiGM);
                mtiVL_sd(i,j) = std(stepMtiVL);
                mtiVM_sd(i,j) = std(stepMtiVM);
                % Mean RMS
                rmsGL_sd(i,j) = std(stepRmsGL);
                rmsGM_sd(i,j) = std(stepRmsGM);
                rmsVL_sd(i,j) = std(stepRmsVL);
                rmsVM_sd(i,j) = std(stepRmsVM);  
                        
        end
        
        meanMnf = figure;
        hold on
        plot(mnfGL(i,:),'r')
        plot(mnfGM(i,:),'y')
        plot(mnfVL(i,:),'b')
        plot(mnfVM(i,:),'c')
        set(meanMnf,'numbertitle','off','name',...
            sprintf('Mean Frequency evolution for the GL (red) GM (yellow), VL (blue) and VM (cyan) during Condition %d',i));
        
        varMnf = figure;
        hold on
        plot(mnfGL_sd(i,:),'r')
        plot(mnfGM_sd(i,:),'y')
        plot(mnfVL_sd(i,:),'b')
        plot(mnfVM_sd(i,:),'c')
        set(varMnf,'numbertitle','off','name',...
            sprintf('Mean Frequency variablity (SD) evolution for the GL (red) GM (yellow), VL (blue) and VM (cyan) during Condition %d',i));
        
        meanMti = figure;
        hold on
        plot(mtiGL(i,:),'r')
        plot(mtiGM(i,:),'y')
        plot(mtiVL(i,:),'b')
        plot(mtiVM(i,:),'c')
        set(meanMti,'numbertitle','off','name',...
            sprintf('Mean total intensity evolution for the GL (red) GM (yellow), VL (blue) and VM (cyan) during Condition %d',i));
        
        varMti = figure;
        hold on
        plot(mtiGL_sd(i,:),'r')
        plot(mtiGM_sd(i,:),'y')
        plot(mtiVL_sd(i,:),'b')
        plot(mtiVM_sd(i,:),'c')
        set(varMti,'numbertitle','off','name',...
            sprintf('Mean total intensity variability (SD) evolution for the GL (red) GM (yellow), VL (blue) and VM (cyan) during Condition %d',i));
                
        meanRms = figure;
        hold on
        plot(rmsGL(i,:),'r')
        plot(rmsGM(i,:),'y')
        plot(rmsVL(i,:),'b')
        plot(rmsVM(i,:),'c')
        ylim([0 25])
        set(meanRms,'numbertitle','off','name',...
            sprintf('Mean RMS evolution for the GL (red) GM (yellow), VL (blue) and VM (cyan) during Condition %d',i));
        
        varRms = figure;
        hold on
        plot(rmsGL_sd(i,:),'r')
        plot(rmsGM_sd(i,:),'y')
        plot(rmsVL_sd(i,:),'b')
        plot(rmsVM_sd(i,:),'c')
        set(varRms,'numbertitle','off','name',...
            sprintf('Mean RMS variablity (SD) evolution for the GL (red) GM (yellow), VL (blue) and VM (cyan) during Condition %d',i));
                
        % Maximize figures
        maxfig(meanMnf,1);
        maxfig(varMnf,1);
        maxfig(meanMti,1);
        maxfig(varMti,1);
        maxfig(meanRms,1);
        maxfig(varRms,1);
                
    end
    
        % Reorganize the values in a table (ROW = trial, COLUMN = session)
        mnfGL = reshape(mnfGL,numSessions,recTime).';
        mnfGM = reshape(mnfGM,numSessions,recTime).';    
        mnfVL = reshape(mnfVL,numSessions,recTime).';        
        mnfVM = reshape(mnfVM,numSessions,recTime).';
        mtiGL = reshape(mtiGL,numSessions,recTime).';
        mtiGM = reshape(mtiGM,numSessions,recTime).';    
        mtiVL = reshape(mtiVL,numSessions,recTime).';        
        mtiVM = reshape(mtiVM,numSessions,recTime).';
        
        % Means between trials
        mnfGLMean = mean(mnfGL);
        mnfGMMean = mean(mnfGM);
        mnfVLMean = mean(mnfVL);        
        mnfVMMean = mean(mnfVM);
        mtiGLMean = mean(mtiGL);
        mtiGMMean = mean(mtiGM);
        mtiVLMean = mean(mtiVL);        
        mtiVMMean = mean(mtiVM);
        
        % Export Results
        walkResults = [mnfGLMean mnfGMMean mnfVLMean mnfVMMean mtiGLMean mtiGMMean mtiVLMean mtiVMMean; ...
            mnfGL mnfGM mnfVL mnfVM mtiGL mtiGM mtiVL mtiVM];
        csvwrite([subjectID '_Walk.csv'],walkResults);
        
    disp('     *****     DONE WITH WALKING TRIALS     *****     ')
        

