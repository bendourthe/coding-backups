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
%   Last update: March 12th, 2018
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


%% Initialization

    clear ; close all; clc
    
%% Directory

    cd 'D:\Dr. Scholls\00\';
    pathName = 'D:\Dr. Scholls\00';

%% Subject ID   

    subjectID = '00';
    
%% Execution  

    % Determine which parts of the code will be run (0: won't be run; 1: will be run)
    currexS = 0;
    currexD = 0;
    mvc = 0;
    half = 1;
    walking = 0;
    pedar = 0;
    
%% Settings
    
    % Sessions
    numSessions = 1;
    session = num2str([1:numSessions].','%01d');
    
    % Filtering parameters
    samplingFrequency = 2400;    % frenquence used during protocol (in Hz)
    scale = 0.3;
    control = 3;
    numWavelets = 13;
    
    % Walking trial
    walkTime = 45;               % in mins
    
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
    
    % Biodex
        % MVC Trials
        mvc1 = [strcat('\',session,'_1_',subjectID,'_mvc_1.csv')];
        mvc2 = [strcat('\',session,'_1_',subjectID,'_mvc_2.csv')];
        mvc3 = [strcat('\',session,'_2_',subjectID,'_mvc_1.csv')];
        mvc_name = [mvc1; mvc2; mvc3];
        % 50p Trials
        half1 = [strcat('\',session,'_1_',subjectID,'_50p_1.csv')];
        half2 = [strcat('\',session,'_1_',subjectID,'_50p_2.csv')];
        half3 = [strcat('\',session,'_2_',subjectID,'_50p_1.csv')];
        half_names = [half1; half2; half3];

    % Walking
        % Pedar
        ped1 = [strcat('\',session,'_',subjectID,'_ped_1.asc')];
        ped2 = [strcat('\',session,'_',subjectID,'_ped_2.asc')];
        pedar_name = [ped1; ped2];
        % EMG
        emg_name = [strcat('\',session,'_',subjectID,'_emg_1.emg')];
    

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
            mean(dataTorque(mvcP(3)-100:mvcP(3)+100));...
            max(dataTorque(mvcP(4)-100:mvcP(4)+100))];
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

        % Peak-to-peak muscle activation (P2P)
            % Select dataPoints around twitch (e.g. 2400 data points
            % corresponds to 1 sec before and 1 sec after the twitch)
            dataPoints = 2400;
            % Calculate the index of the first stimulation (happens within the
            % first 5 secs of the signal -> 12000 data points)
            [M, I1] = max(dataEMG(1:12000));
            % Calculate the index and amplitude of the second stimulation (during MVC,
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
            % Calculate the peak to peak amplitude of the stimulation during
            % the MVC
            mvcP2P(i) = peak2peak(dataEMG(I2(i)-dataPoints:I2(i)+dataPoints));
            % Plot torque data with the selected point for verification
            figure
            clf
            hold on
            plot(dataEMG)
            plot([I2(i)-dataPoints:I2(i)+dataPoints],dataEMG(I2(i)-dataPoints:I2(i)+dataPoints),'r')
            disp('PRESS ENTER if points are correctly selected - PRESS Ctr+C in the command window if not')
            pause

            close all

    end

    % Reorganize the values in a table (ROW = trial, COLUMN = session)
    mvcPT = reshape(mvcPT,numSessions,3).';
    mvcRatio = reshape(mvcRatio,numSessions,3).';
    mvcVAR = reshape(mvcVAR,numSessions,3).';
    mvcP2P = reshape(mvcP2P,numSessions,3).';

    % Means between trials
        % Baseline
        mvcPTMean1 = mean(mvcPT(1:2,:));
        mvcRatioMean1 = mean(mvcRatio(1:2,:));
        mvcVARMean1 = mean(mvcVAR(1:2,:));
        mvcP2PMean1 = mean(mvcP2P(1:2,:));
        % Post-Intervention
        mvcPTMean2 = mean(mvcPT(3,:));
        mvcRatioMean2 = mean(mvcRatio(3,:));
        mvcVARMean2 = mean(mvcVAR(3,:));
        mvcP2PMean2 = mean(mvcP2P(3,:));
    
    % Mean increase/decrease
    for i=1:numSessions
        mvcPTMeanDiff(i) = mvcPTMean2(i)*100/mvcPTMean1(i) - 100;
        mvcRatioMeanDiff(i) = mvcRatioMean2(i)*100/mvcRatioMean1(i) - 100;
        mvcVARMeanDiff(i) = mvcVARMean2(i)*100/mvcVARMean1(i) - 100;
        mvcP2PMeanDiff(i) = mvcP2PMean2(i)*100/mvcP2PMean1(i) - 100;
    end

    % Export Results
    mvcResults = [mvcPT mvcRatio mvcVAR  mvcP2P; mvcPTMean1 mvcRatioMean1 mvcVARMean1 mvcP2PMean1;...
        mvcPTMean2 mvcRatioMean2 mvcVARMean2 mvcP2PMean2; mvcPTMeanDiff mvcRatioMeanDiff mvcVARMeanDiff mvcP2PMeanDiff];
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
        meanT50pMean2 = mean(meanT50p(3,:));
        meanTI50pMean2 = mean(meanTI50p(3,:));        
        MNF50pMean2 = mean(MNF50p(3,:));
    
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
    for i=1:size(emg_name,1)

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
    meanTotIntGM = reshape(meanTotIntGM,numSessions,1).';    
    meanTotIntVL = reshape(meanTotIntVL,numSessions,1).';        
    meanTotIntVM = reshape(meanTotIntVM,numSessions,1).';
    mnfGM = reshape(mnfGM,numSessions,1).';    
    mnfVL = reshape(mnfVL,numSessions,1).';        
    mnfVM = reshape(mnfVM,numSessions,1).';
        
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
        
end
        
    disp('     *****     DONE WITH WALKING TRIALS     *****     ')
       
    
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
