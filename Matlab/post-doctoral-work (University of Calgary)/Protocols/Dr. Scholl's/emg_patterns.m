%% DESCRIPTION:
%
% This code was designed to enable the analysis of the frequency spectrum
% obtained via the Dr. Scholl's Research Project:
%
%   Project Title: Quantifying the effects of different insole
%   configurations on fatigue
%
%   Author: Benjamin Dourthe
%
%   Last update: August 21st, 2018
%
%% Input: 
%   - Frequency spectrum obtained using the DrScholls_Protocol_V4 code
%
%% Output:
%   - 
%
%% Dependencies:
%
%   Files:
%       - Frequency spectrum obtained using the DrScholls_Protocol_V4 code
%
%   Functions:
%       - None


%% Initialization

    clear ; close all; clc
    
%% Directory

    cd 'F:\Dr. Scholls\Phase 4\005\Results\';
    pathName = 'F:\Dr. Scholls\Phase 4\005\Results';

%% Subject ID   

    subjectID = '005';
    
%% Settings

    % Sessions
    numSessions = 3;            % how many sessions to process?
    startSession = 1;           % what is the first session number to process?
    endSession = 3;             % what is the last session number to process?
    sessions = num2str([1:3].','%01d');
    
    % Execution: determine which parts of the code will be run
    % (0: won't be run; 1: will be run)
    patterns = 1;
    bands = 1;

    % Wavelets
    numWavelets = 20;           % how many wavelets were used to process the raw EMG signals

    % Figures properties
    set(groot,'defaultFigureColor','w')
    set(groot,'defaultLineLineWidth',1.5);
    set(groot,'defaultAxesLineWidth',1.5);
    set(groot,'defaultAxesTickDir','out')
    set(groot,'defaultAxesFontSize',15)
    set(groot,'defaultAxesFontWeight','bold')
    
    % Colormap path
    addpath('C:\Users\bdour\Documents\Matlab\Calgary\bluewhitered\');
    
    % Colors for plots
    ctrl = [153/255 153/255 153/255];
    cf = [0/255 204/255 102/255];
    bpr = [0/255 115/255 230/255];
    mg = [0/255 102/255 204/255];
    
%% Filenames

    freq = [strcat('\',subjectID,'_Freq',sessions,'.csv')];
    
%% Execution

% Difference in EMG patterns between conditions
if patterns == 1
    % Load data
    data1 = csvread([pathName freq(1,:)]);
    data2 = csvread([pathName freq(2,:)]);
    data3 = csvread([pathName freq(3,:)]);
    
    % Remove first line (central frequencies)
    dataEMG1 = data1(2:end,:);
    dataEMG2 = data2(2:end,:);
    dataEMG3 = data3(2:end,:);
    
    % Calculates the difference between conditions 2 and 3 and the control
    % (condition 1)
    dataEMG21 = dataEMG2 - dataEMG1;
    dataEMG31 = dataEMG3 - dataEMG1;
    
    % Isolate muscles
    dataGL21 = dataEMG21(1:end,1:numWavelets);
    dataGM21 = dataEMG21(1:end,numWavelets+1:2*numWavelets);
    dataVL21 = dataEMG21(1:end,2*numWavelets+1:3*numWavelets);
    dataVM21 = dataEMG21(1:end,3*numWavelets+1:4*numWavelets);
    dataGL31 = dataEMG31(1:end,1:numWavelets);
    dataGM31 = dataEMG31(1:end,numWavelets+1:2*numWavelets);
    dataVL31 = dataEMG31(1:end,2*numWavelets+1:3*numWavelets);
    dataVM31 = dataEMG31(1:end,3*numWavelets+1:4*numWavelets);
    
    % Reshape full EMG matrices
    dataEMG21 = [dataGL21.' dataGM21.' dataVL21.' dataVM21.'];
    dataEMG31 = [dataGL31.' dataGM31.' dataVL31.' dataVM31.'];
    
    % Time
    time = [4/size(dataEMG21,2):4/size(dataEMG21,2):4];
    
    % Central frequencies
    cfs = data1(1,1:numWavelets);
    
    % Plot corresponding difference patterns
        
        % Condition 2 - Condition 1
        figure
        maxfig(gcf,1)
        contourf(time,cfs,dataEMG21,100,'LineStyle','None')        
        colormap(bluewhitered(256))
            % Add straight lines to seperate muscles
            hold on
            for i=1:3
                plot([i i],[0 max(cfs)],'k')
            end
            % Title and axes
            title('Difference pattern Condition 2 - Condition 1','FontSize',30);
            xticks([0.5 1.5 2.5 3.5])
            xticklabels({'GL','GM','VL','VM'})
            colorbar
            set(gcf,'numbertitle','off','name','Difference pattern Condition 2 - Condition 1');
            % Save figure as PNG
            saveas(gcf,[subjectID '_FullSpectrum_c2-c1.fig']);
            
        % Condition 3 - Condition 1
        figure
        maxfig(gcf,1)
        contourf(time,cfs,dataEMG31,100,'LineStyle','None')        
        colormap(bluewhitered(256))
            % Add straight lines to seperate muscles
            hold on
            for i=1:3
                plot([i i],[0 max(cfs)],'k')
            end
            % Title and axes
            title('Difference pattern Condition 3 - Condition 1','FontSize',30);
            xticks([0.5 1.5 2.5 3.5])
            xticklabels({'GL','GM','VL','VM'})
            colorbar
            set(gcf,'numbertitle','off','name','Difference pattern Condition 3 - Condition 1');
            % Save figure as PNG
            saveas(gcf,[subjectID '_FullSpectrum_c3-c1.fig']);
            
else
end
    
% Time series of each frequency bands
if bands == 1
    
    for i=startSession:endSession

        % Import data
        data = csvread([pathName freq(i,:)]);

        % Central frequencies
        cfs = data(1,1:numWavelets);

        % Isolate muscles
        dataGL = data(2:end,1:numWavelets);
        dataGM = data(2:end,numWavelets+1:2*numWavelets);
        dataVL = data(2:end,2*numWavelets+1:3*numWavelets);
        dataVM = data(2:end,3*numWavelets+1:4*numWavelets);

        % Plot power spectrum corresponding to first and last minute of the
        % intervention
            % GL
            figure
            maxfig(gcf,1)
            plot(cfs,dataGL(1,:),'b')
            hold on
            plot(cfs,dataGL(end,:),'r')
                % Title and axes
                title('GL','FontSize',30);
                xlabel('Frequency (Hz)');   
                set(gcf,'numbertitle','off','name',sprintf('Full spectrum GL - Condition %d',i));
                % Save figure as PNG
                saveas(gcf,sprintf([subjectID '_SpectrumFirstLastMin_GL_%d.png'],i));

            % GM
            figure
            maxfig(gcf,1)
            plot(cfs,dataGM(1,:),'b')
            hold on
            plot(cfs,dataGM(end,:),'r')
                % Title and axes
                title('GM','FontSize',30);
                xlabel('Frequency (Hz)');   
                set(gcf,'numbertitle','off','name',sprintf('Full spectrum GM - Condition %d',i));
                % Save figure as PNG
                saveas(gcf,sprintf([subjectID '_SpectrumFirstLastMin_GM_%d.png'],i));

            % VL
            figure
            maxfig(gcf,1)
            plot(cfs,dataVL(1,:),'b')
            hold on
            plot(cfs,dataVL(end,:),'r')
                % Title and axes
                title('VL','FontSize',30);
                xlabel('Frequency (Hz)');   
                set(gcf,'numbertitle','off','name',sprintf('Full spectrum VL - Condition %d',i));
                % Save figure as PNG
                saveas(gcf,sprintf([subjectID '_SpectrumFirstLastMin_VL_%d.png'],i));

            % VM
            figure
            maxfig(gcf,1)
            plot(cfs,dataVM(1,:),'b')
            hold on
            plot(cfs,dataVM(end,:),'r')
                % Title and axes
                title('VM','FontSize',30);
                xlabel('Frequency (Hz)');   
                set(gcf,'numbertitle','off','name',sprintf('Full spectrum VM - Condition %d',i));
                % Save figure as PNG
                saveas(gcf,sprintf([subjectID '_SpectrumFirstLastMin_VM_%d.png'],i));

    end
    
else
end

