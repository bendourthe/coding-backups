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

%% Subject ID   

    subjectID = '028';
    
%% Directory

    pathName = ['F:\Dr. Scholls\Phase 3\' subjectID '\Results'];
    cd(pathName);
    
%% Settings

    % Sessions
    sessions = num2str([1:3].','%01d');
    
    % Define third session: type 1 for ctrl, 2 for cf, 3 for bpr, 4 for mg
    condition3 = 4;

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
    mg = [255/255 128/255 0/255];
    if condition3 == 1
        c3 = ctrl;
    elseif condition3 == 2
        c3 = cf;
    elseif condition3 == 3
        c3 = bpr;
    elseif condition3 == 4
        c3 = mg;
    end
    
%% Filenames

    freq = [strcat('\',subjectID,'_Freq',sessions,'.csv')];
    
%% Execution

    % Load data
    data1 = csvread([pathName freq(1,:)]);
    data2 = csvread([pathName freq(2,:)]);
    data3 = csvread([pathName freq(3,:)]);
    
    % Central frequencies
    cfs = data1(1,1:numWavelets);
    
    % Remove first line (central frequencies)
    dataEMG1 = data1(2:end,:);
    dataEMG2 = data2(2:end,:);
    dataEMG3 = data3(2:end,:);
    
%% Difference in EMG patterns between conditions

    % Isolate muscles
        % Conditions
        dataGL1 = dataEMG1(1:end,1:numWavelets)/mean(mean(dataEMG1(1:end,1:numWavelets)));
        dataGM1 = dataEMG1(1:end,numWavelets+1:2*numWavelets)/mean(mean(dataEMG1(1:end,numWavelets+1:2*numWavelets)));
        dataVL1 = dataEMG1(1:end,2*numWavelets+1:3*numWavelets)/mean(mean(dataEMG1(1:end,2*numWavelets+1:3*numWavelets)));
        dataVM1 = dataEMG1(1:end,3*numWavelets+1:4*numWavelets)/mean(mean(dataEMG1(1:end,3*numWavelets+1:4*numWavelets)));
        dataGL2 = dataEMG2(1:end,1:numWavelets)/mean(mean(dataEMG2(1:end,1:numWavelets)));
        dataGM2 = dataEMG2(1:end,numWavelets+1:2*numWavelets)/mean(mean(dataEMG2(1:end,numWavelets+1:2*numWavelets)));
        dataVL2 = dataEMG2(1:end,2*numWavelets+1:3*numWavelets)/mean(mean(dataEMG2(1:end,2*numWavelets+1:3*numWavelets)));
        dataVM2 = dataEMG2(1:end,3*numWavelets+1:4*numWavelets)/mean(mean(dataEMG2(1:end,3*numWavelets+1:4*numWavelets)));
        dataGL3 = dataEMG3(1:end,1:numWavelets)/mean(mean(dataEMG3(1:end,1:numWavelets)));
        dataGM3 = dataEMG3(1:end,numWavelets+1:2*numWavelets)/mean(mean(dataEMG3(1:end,numWavelets+1:2*numWavelets)));
        dataVL3 = dataEMG3(1:end,2*numWavelets+1:3*numWavelets)/mean(mean(dataEMG3(1:end,2*numWavelets+1:3*numWavelets)));
        dataVM3 = dataEMG3(1:end,3*numWavelets+1:4*numWavelets)/mean(mean(dataEMG3(1:end,3*numWavelets+1:4*numWavelets)));
        % Differences
        dataGL21 = dataGL2 - dataGL1;
        dataGM21 = dataGM2 - dataGM1;
        dataVL21 = dataVL2 - dataVL1;
        dataVM21 = dataVM2 - dataVM1;        
        dataGL31 = dataGL3 - dataGL1;
        dataGM31 = dataGM3 - dataGM1;
        dataVL31 = dataVL3 - dataVL1;
        dataVM31 = dataVM3 - dataVM1;
        
    % Reshape full EMG matrices
        % Conditions
        dataEMG1 = [dataGL1.' dataGM1.' dataVL1.' dataVM1.'];
        dataEMG2 = [dataGL2.' dataGM2.' dataVL2.' dataVM2.'];
        dataEMG3 = [dataGL3.' dataGM3.' dataVL3.' dataVM3.'];
        % Differences
        dataEMG21 = [dataGL21.' dataGM21.' dataVL21.' dataVM21.'];
        dataEMG31 = [dataGL31.' dataGM31.' dataVL31.' dataVM31.'];
        
    % Define limits for graphs
        % Pattern plots
        cmin = 0;
        cmax = max([max(max(dataGL1)) max(max(dataGL2)) max(max(dataGL3)) ...
            max(max(dataGM1)) max(max(dataGM2)) max(max(dataGM3)) ...
            max(max(dataVL1)) max(max(dataVL2)) max(max(dataVL3)) ...
            max(max(dataVM1)) max(max(dataVM2)) max(max(dataVM3))]);
        % Difference pattern plots
        diff_min = min([min(min(dataEMG21)) min(min(dataEMG31))]);
        diff_max = max([max(max(dataEMG21)) max(max(dataEMG31))]);
        cmin_diff = -max([abs(diff_min) abs(diff_max)]);
        cmax_diff = max([abs(diff_min) abs(diff_max)]);
        % First vs last plots
        lim_y_first = max([max(max(dataGL1(1,:))) max(max(dataGL2(1,:))) max(max(dataGL3(1,:))) ...
            max(max(dataGM1(1,:))) max(max(dataGM2(1,:))) max(max(dataGM3(1,:))) ...
            max(max(dataVL1(1,:))) max(max(dataVL2(1,:))) max(max(dataVL3(1,:))) ...
            max(max(dataVM1(1,:))) max(max(dataVM2(1,:))) max(max(dataVM3(1,:)))]);
        lim_y_last = max([max(max(dataGL1(end,:))) max(max(dataGL2(end,:))) max(max(dataGL3(end,:))) ...
            max(max(dataGM1(end,:))) max(max(dataGM2(end,:))) max(max(dataGM3(end,:))) ...
            max(max(dataVL1(end,:))) max(max(dataVL2(end,:))) max(max(dataVL3(end,:))) ...
            max(max(dataVM1(end,:))) max(max(dataVM2(end,:))) max(max(dataVM3(end,:)))]);
        lim_y = 1.2*max([lim_y_first lim_y_last]);
    
    % Time
    time = [4/size(dataEMG21,2):4/size(dataEMG21,2):4];
    
    % Locate peaks of each power spectrum (1 peak per minute)
    [mGL1, iGL1] = max(dataGL1.');
    [mGM1, iGM1] = max(dataGM1.');
    [mVL1, iVL1] = max(dataVL1.');
    [mVM1, iVM1] = max(dataVM1.');
    [mGL2, iGL2] = max(dataGL2.');
    [mGM2, iGM2] = max(dataGM2.');
    [mVL2, iVL2] = max(dataVL2.');
    [mVM2, iVM2] = max(dataVM2.');
    [mGL3, iGL3] = max(dataGL3.');
    [mGM3, iGM3] = max(dataGM3.');
    [mVL3, iVL3] = max(dataVL3.');
    [mVM3, iVM3] = max(dataVM3.');
        
    % Plot corresponding patterns (1 per condition)
        
        % Condition 1
        figure
        maxfig(gcf,1)
        contourf(time,cfs,dataEMG1,100,'LineStyle','None')        
        colormap(redblue)
            % Add straight lines to seperate muscles
            hold on
            for i=1:3
                plot([i i],[0 max(cfs)],'k')
            end
            % Add lines to show where the peak of each spectrum is located
            plot([4/size(dataEMG21,2):4/size(dataEMG21,2):4],cfs([iGL1 iGM1 iVL1 iVM1]),'k')
            % Title and axes
            title('Intensity pattern Condition 1','FontSize',30);
            xticks([0.5 1.5 2.5 3.5])
            xticklabels({'GL','GM','VL','VM'})
            colorbar
            caxis([cmin cmax])
            set(gcf,'numbertitle','off','name','Intensity pattern Condition 1');
            % Save figure as JPG
            saveas(gcf,[subjectID '_FullSpectrum_c1.jpg']);
            
        % Condition 2
        figure
        maxfig(gcf,1)
        contourf(time,cfs,dataEMG2,100,'LineStyle','None')        
        colormap(redblue)
            % Add straight lines to seperate muscles
            hold on
            for i=1:3
                plot([i i],[0 max(cfs)],'k')
            end
            % Add lines to show where the peak of each spectrum is located
            plot([4/size(dataEMG21,2):4/size(dataEMG21,2):4],cfs([iGL2 iGM2 iVL2 iVM2]),'k')
            % Title and axes
            title('Intensity pattern Condition 2','FontSize',30);
            xticks([0.5 1.5 2.5 3.5])
            xticklabels({'GL','GM','VL','VM'})
            colorbar
            caxis([cmin cmax])
            set(gcf,'numbertitle','off','name','Intensity pattern Condition 2');
            % Save figure as JPG
            saveas(gcf,[subjectID '_FullSpectrum_c2.jpg']);
            
        % Condition 3
        figure
        maxfig(gcf,1)
        contourf(time,cfs,dataEMG3,100,'LineStyle','None')        
        colormap(redblue)
            % Add straight lines to seperate muscles
            hold on
            for i=1:3
                plot([i i],[0 max(cfs)],'k')
            end
            % Add lines to show where the peak of each spectrum is located
            plot([4/size(dataEMG21,2):4/size(dataEMG21,2):4],cfs([iGL3 iGM3 iVL3 iVM3]),'k')
            % Title and axes
            title('Intensity pattern Condition 3','FontSize',30);
            xticks([0.5 1.5 2.5 3.5])
            xticklabels({'GL','GM','VL','VM'})
            colorbar
            caxis([cmin cmax])
            set(gcf,'numbertitle','off','name','Intensity pattern Condition 3');
            % Save figure as JPG
            saveas(gcf,[subjectID '_FullSpectrum_c3.jpg']);
    
    % Plot corresponding difference patterns
        
        % Condition 2 - Condition 1
        figure
        maxfig(gcf,1)
        contourf(time,cfs,dataEMG21,100,'LineStyle','None')        
        colormap(redblue)
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
            caxis([cmin_diff cmax_diff])
            set(gcf,'numbertitle','off','name','Difference pattern Condition 2 - Condition 1');
            % Save figure as JPG
            saveas(gcf,[subjectID '_FullSpectrum_c2-c1.jpg']);
            
        % Condition 3 - Condition 1
        figure
        maxfig(gcf,1)
        contourf(time,cfs,dataEMG31,100,'LineStyle','None')        
        colormap(redblue)
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
            caxis([cmin_diff cmax_diff])
            set(gcf,'numbertitle','off','name','Difference pattern Condition 3 - Condition 1');
            % Save figure as JPG
            saveas(gcf,[subjectID '_FullSpectrum_c3-c1.jpg']);
            
%% Time series of each frequency bands

    % Smoothing data
    [scfs, sGL1] = smoothData(cfs,dataGL1,10);
    [scfs, sGM1] = smoothData(cfs,dataGM1,10);
    [scfs, sVL1] = smoothData(cfs,dataVL1,10);
    [scfs, sVM1] = smoothData(cfs,dataVM1,10);
    [scfs, sGL2] = smoothData(cfs,dataGL2,10);
    [scfs, sGM2] = smoothData(cfs,dataGM2,10);
    [scfs, sVL2] = smoothData(cfs,dataVL2,10);
    [scfs, sVM2] = smoothData(cfs,dataVM2,10);
    [scfs, sGL3] = smoothData(cfs,dataGL3,10);
    [scfs, sGM3] = smoothData(cfs,dataGM3,10);
    [scfs, sVL3] = smoothData(cfs,dataVL3,10);
    [scfs, sVM3] = smoothData(cfs,dataVM3,10);

    % Plot power spectrum corresponding to first and last minute of the
    % intervention
        % GL
        figure
        maxfig(gcf,1)
        plot(scfs,sGL1(1,:),'--','color', ctrl)
        hold on
        plot(scfs,sGL1(end,:),'color', ctrl)
        plot(scfs,sGL2(1,:),'--','color', cf)
        plot(scfs,sGL2(end,:),'color', cf)
        plot(scfs,sGL3(1,:),'--','color', c3)
        plot(scfs,sGL3(end,:),'color', c3)
            % Title and axes
            title('GL','FontSize',30);
            xlabel('Frequency (Hz)');
            ylim([0 lim_y])
            set(gcf,'numbertitle','off','name','Full spectrum GL');
            % Save figure as PNG
            saveas(gcf,[subjectID '_SpectrumFirstLastMin_GL.png']);

        % GM
        figure
        maxfig(gcf,1)
        plot(scfs,sGM1(1,:),'--','color', ctrl)
        hold on
        plot(scfs,sGM1(end,:),'color', ctrl)
        plot(scfs,sGM2(1,:),'--','color', cf)
        plot(scfs,sGM2(end,:),'color', cf)
        plot(scfs,sGM3(1,:),'--','color', c3)
        plot(scfs,sGM3(end,:),'color', c3)
            % Title and axes
            title('GM','FontSize',30);
            xlabel('Frequency (Hz)');
            ylim([0 lim_y])
            set(gcf,'numbertitle','off','name','Full spectrum GM');
            % Save figure as PNG
            saveas(gcf,[subjectID '_SpectrumFirstLastMin_GM.png']);

        % VL
        figure
        maxfig(gcf,1)
        plot(scfs,sVL1(1,:),'--','color', ctrl)
        hold on
        plot(scfs,sVL1(end,:),'color', ctrl)
        plot(scfs,sVL2(1,:),'--','color', cf)
        plot(scfs,sVL2(end,:),'color', cf)
        plot(scfs,sVL3(1,:),'--','color', c3)
        plot(scfs,sVL3(end,:),'color', c3)
            % Title and axes
            title('VL','FontSize',30);
            xlabel('Frequency (Hz)');  
            ylim([0 lim_y])
            set(gcf,'numbertitle','off','name','Full spectrum VL');
            % Save figure as PNG
            saveas(gcf,[subjectID '_SpectrumFirstLastMin_VL.png']);

        % VM
        figure
        maxfig(gcf,1)
        plot(scfs,sVM1(1,:),'--','color', ctrl)
        hold on
        plot(scfs,sVM1(end,:),'color', ctrl)
        plot(scfs,sVM2(1,:),'--','color', cf)
        plot(scfs,sVM2(end,:),'color', cf)
        plot(scfs,sVM3(1,:),'--','color', c3)
        plot(scfs,sVM3(end,:),'color', c3)
            % Title and axes
            title('VM','FontSize',30);
            xlabel('Frequency (Hz)');
            ylim([0 lim_y])
            set(gcf,'numbertitle','off','name','Full spectrum VM');
            % Save figure as PNG
            saveas(gcf,[subjectID '_SpectrumFirstLastMin_VM.png']);
    
    disp('     *****     DONE WITH EMG SPECTRUM ANALYSIS     *****     ')
