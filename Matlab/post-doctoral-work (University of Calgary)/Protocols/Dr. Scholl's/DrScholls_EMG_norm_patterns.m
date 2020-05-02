%% Initialization

    clear ; close all; clc

%% Subject ID   

    subjectID = '002';
    
%% Directory

    pathName = ['F:\Dr. Scholls\Phase 3\' subjectID '\Results'];
    cd(pathName);
    
%% Settings

    % Sessions
    sessions = num2str([1:3].','%01d');

    % Number of data points for normalization
    dataPoints = 50;
    
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
        % Differences
        dataEMG21 = [dataGL21.' dataGM21.' dataVL21.' dataVM21.'];
        dataEMG31 = [dataGL31.' dataGM31.' dataVL31.' dataVM31.'];
        
    % Define limits for graphs
        % Difference pattern plots
        diff_min = min([min(min(dataEMG21)) min(min(dataEMG31))]);
        diff_max = max([max(max(dataEMG21)) max(max(dataEMG31))]);
        cmin_diff = -max([abs(diff_min) abs(diff_max)]);
        cmax_diff = max([abs(diff_min) abs(diff_max)]);
    
    % Time
    time = [4/size(dataEMG21,2):4/size(dataEMG21,2):4];
    
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
    
%% Difference in EMG patterns between conditions (normalized to less data points)

for i = 1:numWavelets
    dataGL1norm(:,i) = norm_sig(dataGL1(:,i),dataPoints);
    dataGM1norm(:,i) = norm_sig(dataGM1(:,i),dataPoints);
    dataVL1norm(:,i) = norm_sig(dataVL1(:,i),dataPoints);
    dataVM1norm(:,i) = norm_sig(dataVM1(:,i),dataPoints);
    dataGL2norm(:,i) = norm_sig(dataGL2(:,i),dataPoints);
    dataGM2norm(:,i) = norm_sig(dataGM2(:,i),dataPoints);
    dataVL2norm(:,i) = norm_sig(dataVL2(:,i),dataPoints);
    dataVM2norm(:,i) = norm_sig(dataVM2(:,i),dataPoints);
end

dataGL21norm = dataGL2norm - dataGL1norm;
dataGM21norm = dataGM2norm - dataGM1norm;
dataVL21norm = dataVL2norm - dataVL1norm;
dataVM21norm = dataVM2norm - dataVM1norm;  

dataEMG21norm = [dataGL21norm.' dataGM21norm.' dataVL21norm.' dataVM21norm.'];

timeNorm = [4/size(dataEMG21norm,2):4/size(dataEMG21norm,2):4];

figure
maxfig(gcf,1)
contourf(timeNorm,cfs,dataEMG21norm,100,'LineStyle','None')        
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