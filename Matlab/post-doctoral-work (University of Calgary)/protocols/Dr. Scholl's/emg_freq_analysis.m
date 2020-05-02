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
%   Last update: August 28th, 2018
%
%% Input: 
%   - Frequency spectrum obtained using the DrScholls_Protocol_V4 code
%   - Times when subject transitioned toward each new fatigue stage
%
%% Output:
%   - Graphs showing the evolution of each frequency band with the trends
%   for each corresponding fatigue stage
%   - Matrix including the slopes of each trend
%   - Matrix including the R2 of each trend
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
    numSessions = 3;            % how many sessions to process?
    startSession = 1;           % what is the first session number to process?
    endSession = 3;             % what is the last session number to process?
    sessions = num2str([1:3].','%01d');
    
    % Perceived fatigue
        % Each line is a condition, each column corresponds to the time at
        % which each fatigue stage was reached (column 1 = stage 1, column
        % 2 = stage 2, etc.). Last column should be the recording time.
    fatigue = [1 15 25 25 25; 1 5 10 15 25; 1 5 10 25 25];
    
    % Plots
    raw = 1;                    % type 1 if you want to see the raw frequency time series, type 0 otherwise
    trend = 1;                  % type 1 if you want to see the trend lines calculated from the raw frequency time series,
                                % type 0 otherwise
    
    % Wavelets
    numWavelets = 20;           % how many wavelets were used to process the raw EMG signals

    % Figures properties
    set(groot,'defaultFigureColor','w')
    set(groot,'defaultLineLineWidth',1.5);
    set(groot,'defaultAxesLineWidth',1.5);
    set(groot,'defaultAxesTickDir','out')
    set(groot,'defaultAxesFontSize',15)
    set(groot,'defaultAxesFontWeight','bold')
    
    % Colors
    band1 = [204/255 153/255 0/255];
    band2 = [0/255 204/255 102/255];
    band3 = [0/255 204/255 204/255];
    band4 = [0/255 102/255 204/255];
    
    % Plot colors and legend
    figure
    maxfig(gcf,1)
    X = [7 35 55 140; 25 45 95 246];
    Y = [1 1 1 1; 1 1 1 1];
    bands = area(X,Y);
    bands(1).FaceColor = band1;
    bands(2).FaceColor = band2;
    bands(3).FaceColor = band3;
    bands(4).FaceColor = band4;
    xlabel('Frequency (Hz)')
    xticklabels({'7','25','35','45','55','95','140','246'})
    xticks([7 25 35 45 55 95 140 246]);
    yticklabels({'0','1'})
    yticks([0 1]);
    legend('7 to 25 Hz','35 to 45 Hz','55 to 95 Hz','140 to 246 Hz','Location','bestoutside')
    % Save figure as PNG
    saveas(gcf,sprintf([subjectID '_Bands.png']));
    
%% Filenames

    freq = [strcat('\',subjectID,'_Freq',sessions,'.csv')];
    
%% Execution

    % Define lim for graphs
    for i=startSession:endSession
        
        % Import data
        data = csvread([pathName freq(i,:)]);
        
        % Isolate muscles
        dataGL = data(2:end,1:numWavelets)/mean(mean(data(2:end,1:numWavelets)));
        dataGM = data(2:end,numWavelets+1:2*numWavelets)/mean(mean(data(2:end,numWavelets+1:2*numWavelets)));
        dataVL = data(2:end,2*numWavelets+1:3*numWavelets)/mean(mean(data(2:end,2*numWavelets+1:3*numWavelets)));
        dataVM = data(2:end,3*numWavelets+1:4*numWavelets)/mean(mean(data(2:end,3*numWavelets+1:4*numWavelets)));

        % Isolate bands
            % GL
            bandGL_7_25(:,i) = mean(dataGL(:,3:6).').';
            bandGL_35_45(:,i) = mean(dataGL(:,7:8).').';
            bandGL_55_95(:,i) = mean(dataGL(:,9:12).').';
            bandGL_140_246(:,i) = mean(dataGL(:,15:20).').';
            maxGL = max([bandGL_7_25(:,i); bandGL_35_45(:,i); bandGL_55_95(:,i); bandGL_140_246(:,i)]);
            % GM
            bandGM_7_25(:,i) = mean(dataGM(:,3:6).').';
            bandGM_35_45(:,i) = mean(dataGM(:,7:8).').';
            bandGM_55_95(:,i) = mean(dataGM(:,9:12).').';
            bandGM_140_246(:,i) = mean(dataGM(:,15:20).').';
            maxGM = max([bandGM_7_25(:,i); bandGM_35_45(:,i); bandGM_55_95(:,i); bandGM_140_246(:,i)]);
            % GL
            bandVL_7_25(:,i) = mean(dataVL(:,3:6).').';
            bandVL_35_45(:,i) = mean(dataVL(:,7:8).').';
            bandVL_55_95(:,i) = mean(dataVL(:,9:12).').';
            bandVL_140_246(:,i) = mean(dataVL(:,15:20).').';
            maxVL = max([bandVL_7_25(:,i); bandVL_35_45(:,i); bandVL_55_95(:,i); bandVL_140_246(:,i)]);
            % GL
            bandVM_7_25(:,i) = mean(dataVM(:,3:6).').';
            bandVM_35_45(:,i) = mean(dataVM(:,7:8).').';
            bandVM_55_95(:,i) = mean(dataVM(:,9:12).').';
            bandVM_140_246(:,i) = mean(dataVM(:,15:20).').';
            maxVM = max([bandVM_7_25(:,i); bandVM_35_45(:,i); bandVM_55_95(:,i); bandVM_140_246(:,i)]);

        % Find max intensity among all bands
        lim_y(i) = 1.2*max([maxGL maxGM maxVL maxVM]);
        
    end
    lim_x = size(data,1)-1;
    lim_y = max(lim_y);
    clear data

    % Frequency bands analysis
    for i=startSession:endSession

        % Import data
        data = csvread([pathName freq(i,:)]);

        % Central frequencies
        cfs = data(1,1:numWavelets);

        % Isolate muscles
        dataGL = data(2:end,1:numWavelets)/mean(mean(data(2:end,1:numWavelets)));
        dataGM = data(2:end,numWavelets+1:2*numWavelets)/mean(mean(data(2:end,numWavelets+1:2*numWavelets)));
        dataVL = data(2:end,2*numWavelets+1:3*numWavelets)/mean(mean(data(2:end,2*numWavelets+1:3*numWavelets)));
        dataVM = data(2:end,3*numWavelets+1:4*numWavelets)/mean(mean(data(2:end,3*numWavelets+1:4*numWavelets)));

        % Frequency bands analysis
            % GL
            bandGL_7_25(:,i) = mean(dataGL(:,3:6).').';
            bandGL_35_45(:,i) = mean(dataGL(:,7:8).').';
            bandGL_55_95(:,i) = mean(dataGL(:,9:12).').';
            bandGL_140_246(:,i) = mean(dataGL(:,15:20).').';
            maxGL(i) = max([bandGL_7_25(:,i); bandGL_35_45(:,i); bandGL_55_95(:,i); bandGL_140_246(:,i)]);
            % GM
            bandGM_7_25(:,i) = mean(dataGM(:,3:6).').';
            bandGM_35_45(:,i) = mean(dataGM(:,7:8).').';
            bandGM_55_95(:,i) = mean(dataGM(:,9:12).').';
            bandGM_140_246(:,i) = mean(dataGM(:,15:20).').';
            maxGM(i) = max([bandGM_7_25(:,i); bandGM_35_45(:,i); bandGM_55_95(:,i); bandGM_140_246(:,i)]);
            % GL
            bandVL_7_25(:,i) = mean(dataVL(:,3:6).').';
            bandVL_35_45(:,i) = mean(dataVL(:,7:8).').';
            bandVL_55_95(:,i) = mean(dataVL(:,9:12).').';
            bandVL_140_246(:,i) = mean(dataVL(:,15:20).').';
            maxVL(i) = max([bandVL_7_25(:,i); bandVL_35_45(:,i); bandVL_55_95(:,i); bandVL_140_246(:,i)]);
            % GL
            bandVM_7_25(:,i) = mean(dataVM(:,3:6).').';
            bandVM_35_45(:,i) = mean(dataVM(:,7:8).').';
            bandVM_55_95(:,i) = mean(dataVM(:,9:12).').';
            bandVM_140_246(:,i) = mean(dataVM(:,15:20).').';
            maxVM(i) = max([bandVM_7_25(:,i); bandVM_35_45(:,i); bandVM_55_95(:,i); bandVM_140_246(:,i)]);

            % Plot time series of each frequency band, muscle and condition
            %   (4 bands per figure, 1 figure per muscle and condition)

                % GL
                figure
                maxfig(gcf,1)
                    % Plot the intensity of each frequency band over time
                    if raw == 1
                        plot(bandGL_7_25(:,i),'color',band1)
                        hold on
                        plot(bandGL_35_45(:,i),'color',band2)
                        plot(bandGL_55_95(:,i),'color',band3)
                        plot(bandGL_140_246(:,i),'color',band4)
                    end
                    % Title and axes
                    title('GL','FontSize',30);
                    xlabel('Time (mins)');
                    xlim([0 lim_x])
                    ylim([0 lim_y])     
                    set(gcf,'numbertitle','off','name',sprintf('GL - Condition %d',i));
                    % Calculate trend lines for each fatigue section and each
                    % frequency band
                    if trend == 1
                        % 7 to 25 Hz
                        for j = 1:size(fatigue,2)-1
                            coef = polyfit(fatigue(i,j):fatigue(i,j+1), ...
                                bandGL_7_25(fatigue(i,j):fatigue(i,j+1),i).', 1);
                            xFitting = fatigue(i,j):fatigue(i,j+1); 
                            yFitted = polyval(coef, xFitting);
                            coeffs_GL_7_25(i,j) = yFitted(end)/yFitted(1)*100-100;
                            R = corr2(bandGL_7_25(fatigue(i,j):fatigue(i,j+1),i),yFitted.');                        
                            R2_GL_7_25(i,j) = R^2;
                            plot(xFitting, yFitted, '--','color', band1);
                        end
                        % 35 to 45 Hz
                        for j = 1:size(fatigue,2)-1
                            coef = polyfit(fatigue(i,j):fatigue(i,j+1), ...
                                bandGL_35_45(fatigue(i,j):fatigue(i,j+1),i).', 1);
                            xFitting = fatigue(i,j):fatigue(i,j+1); 
                            yFitted = polyval(coef, xFitting);
                            coeffs_GL_35_45(i,j) = yFitted(end)/yFitted(1)*100-100;
                            R = corr2(bandGL_35_45(fatigue(i,j):fatigue(i,j+1),i),yFitted.');
                            R2_GL_35_45(i,j) = R^2;
                            plot(xFitting, yFitted, '--','color', band2);
                        end
                        % 55 to 95 Hz
                        for j = 1:size(fatigue,2)-1
                            coef = polyfit(fatigue(i,j):fatigue(i,j+1), ...
                                bandGL_55_95(fatigue(i,j):fatigue(i,j+1),i).', 1);
                            xFitting = fatigue(i,j):fatigue(i,j+1); 
                            yFitted = polyval(coef, xFitting);
                            coeffs_GL_55_95(i,j) = yFitted(end)/yFitted(1)*100-100;
                            R = corr2(bandGL_55_95(fatigue(i,j):fatigue(i,j+1),i),yFitted.');
                            R2_GL_55_95(i,j) = R^2;
                            plot(xFitting, yFitted, '--','color', band3);
                        end
                        % 140 to 246 Hz
                        for j = 1:size(fatigue,2)-1
                            coef = polyfit(fatigue(i,j):fatigue(i,j+1), ...
                                bandGL_140_246(fatigue(i,j):fatigue(i,j+1),i).', 1);
                            xFitting = fatigue(i,j):fatigue(i,j+1); 
                            yFitted = polyval(coef, xFitting);
                            coeffs_GL_140_246(i,j) = yFitted(end)/yFitted(1)*100-100;
                            R = corr2(bandGL_140_246(fatigue(i,j):fatigue(i,j+1),i),yFitted.');
                            R2_GL_140_246(i,j) = R^2;
                            plot(xFitting, yFitted, '--','color', band4);
                        end            
                    end
                    % Plot straight lines when perceived fatigue changes
                    for j = 2:size(fatigue,2)
                        hold on
                        plot([fatigue(i,j) fatigue(i,j)],[0 lim_y],'k')
                    end
                    % Save figure as JPG
                    saveas(gcf,sprintf([subjectID '_Bands_GL_%d.jpg'],i));

                % GM
                figure
                maxfig(gcf,1)
                    % Plot the intensity of each frequency band over time
                    if raw == 1
                        plot(bandGM_7_25(:,i),'color',band1)
                        hold on
                        plot(bandGM_35_45(:,i),'color',band2)
                        plot(bandGM_55_95(:,i),'color',band3)
                        plot(bandGM_140_246(:,i),'color',band4)
                    end
                    % Title and axes
                    title('GM','FontSize',30);
                    xlabel('Time (mins)');
                    xlim([0 lim_x])
                    ylim([0 lim_y])     
                    set(gcf,'numbertitle','off','name',sprintf('GM - Condition %d',i));
                    % Calculate trend lines for each fatigue section and each
                    % frequency band
                    if trend == 1
                        % 7 to 25 Hz
                        for j = 1:size(fatigue,2)-1
                            coef = polyfit(fatigue(i,j):fatigue(i,j+1), ...
                                bandGM_7_25(fatigue(i,j):fatigue(i,j+1),i).', 1);
                            xFitting = fatigue(i,j):fatigue(i,j+1); 
                            yFitted = polyval(coef, xFitting);
                            coeffs_GM_7_25(i,j) = yFitted(end)/yFitted(1)*100-100;
                            R = corr2(bandGM_7_25(fatigue(i,j):fatigue(i,j+1),i),yFitted.');
                            R2_GM_7_25(i,j) = R^2;
                            plot(xFitting, yFitted, '--','color', band1);
                        end
                        % 35 to 45 Hz
                        for j = 1:size(fatigue,2)-1
                            coef = polyfit(fatigue(i,j):fatigue(i,j+1), ...
                                bandGM_35_45(fatigue(i,j):fatigue(i,j+1),i).', 1);
                            xFitting = fatigue(i,j):fatigue(i,j+1); 
                            yFitted = polyval(coef, xFitting);
                            coeffs_GM_35_45(i,j) = yFitted(end)/yFitted(1)*100-100;
                            R = corr2(bandGM_35_45(fatigue(i,j):fatigue(i,j+1),i),yFitted.');
                            R2_GM_35_45(i,j) = R^2;
                            plot(xFitting, yFitted, '--','color', band2);
                        end
                        % 55 to 95 Hz
                        for j = 1:size(fatigue,2)-1
                            coef = polyfit(fatigue(i,j):fatigue(i,j+1), ...
                                bandGM_55_95(fatigue(i,j):fatigue(i,j+1),i).', 1);
                            xFitting = fatigue(i,j):fatigue(i,j+1); 
                            yFitted = polyval(coef, xFitting);
                            coeffs_GM_55_95(i,j) = yFitted(end)/yFitted(1)*100-100;
                            R = corr2(bandGM_55_95(fatigue(i,j):fatigue(i,j+1),i),yFitted.');
                            R2_GM_55_95(i,j) = R^2;
                            plot(xFitting, yFitted, '--','color', band3);
                        end
                        % 140 to 246 Hz
                        for j = 1:size(fatigue,2)-1
                            coef = polyfit(fatigue(i,j):fatigue(i,j+1), ...
                                bandGM_140_246(fatigue(i,j):fatigue(i,j+1),i).', 1);
                            xFitting = fatigue(i,j):fatigue(i,j+1); 
                            yFitted = polyval(coef, xFitting);
                            coeffs_GM_140_246(i,j) = yFitted(end)/yFitted(1)*100-100;
                            R = corr2(bandGM_140_246(fatigue(i,j):fatigue(i,j+1),i),yFitted.');
                            R2_GM_140_246(i,j) = R^2;
                            plot(xFitting, yFitted, '--','color', band4);
                        end            
                    end
                    % Plot straight lines when perceived fatigue changes
                    for j = 2:size(fatigue,2)
                        hold on
                        plot([fatigue(i,j) fatigue(i,j)],[0 lim_y],'k')
                    end    
                    % Save figure as JPG
                    saveas(gcf,sprintf([subjectID '_Bands_GM_%d.jpg'],i));

                % VL
                figure
                maxfig(gcf,1)
                    % Plot the intensity of each frequency band over time
                    if raw == 1
                        plot(bandVL_7_25(:,i),'color',band1)
                        hold on
                        plot(bandVL_35_45(:,i),'color',band2)
                        plot(bandVL_55_95(:,i),'color',band3)
                        plot(bandVL_140_246(:,i),'color',band4)
                    end
                    % Title and axes
                    title('VL','FontSize',30);
                    xlabel('Time (mins)');
                    xlim([0 lim_x])
                    ylim([0 lim_y])     
                    set(gcf,'numbertitle','off','name',sprintf('VL - Condition %d',i));
                    % Calculate trend lines for each fatigue section and each
                    % frequency band
                    if trend == 1
                        % 7 to 25 Hz
                        for j = 1:size(fatigue,2)-1
                            coef = polyfit(fatigue(i,j):fatigue(i,j+1), ...
                                bandVL_7_25(fatigue(i,j):fatigue(i,j+1),i).', 1);
                            xFitting = fatigue(i,j):fatigue(i,j+1); 
                            yFitted = polyval(coef, xFitting);
                            coeffs_VL_7_25(i,j) = yFitted(end)/yFitted(1)*100-100;
                            R = corr2(bandVL_7_25(fatigue(i,j):fatigue(i,j+1),i),yFitted.');
                            R2_VL_7_25(i,j) = R^2;
                            plot(xFitting, yFitted, '--','color', band1);
                        end
                        % 35 to 45 Hz
                        for j = 1:size(fatigue,2)-1
                            coef = polyfit(fatigue(i,j):fatigue(i,j+1), ...
                                bandVL_35_45(fatigue(i,j):fatigue(i,j+1),i).', 1);
                            xFitting = fatigue(i,j):fatigue(i,j+1); 
                            yFitted = polyval(coef, xFitting);
                            coeffs_VL_35_45(i,j) = yFitted(end)/yFitted(1)*100-100;
                            R = corr2(bandVL_35_45(fatigue(i,j):fatigue(i,j+1),i),yFitted.');
                            R2_VL_35_45(i,j) = R^2;
                            plot(xFitting, yFitted, '--','color', band2);
                        end
                        % 55 to 95 Hz
                        for j = 1:size(fatigue,2)-1
                            coef = polyfit(fatigue(i,j):fatigue(i,j+1), ...
                                bandVL_55_95(fatigue(i,j):fatigue(i,j+1),i).', 1);
                            xFitting = fatigue(i,j):fatigue(i,j+1); 
                            yFitted = polyval(coef, xFitting);
                            coeffs_VL_55_95(i,j) = yFitted(end)/yFitted(1)*100-100;
                            R = corr2(bandVL_55_95(fatigue(i,j):fatigue(i,j+1),i),yFitted.');
                            R2_VL_55_95(i,j) = R^2;
                            plot(xFitting, yFitted, '--','color', band3);
                        end
                        % 140 to 246 Hz
                        for j = 1:size(fatigue,2)-1
                            coef = polyfit(fatigue(i,j):fatigue(i,j+1), ...
                                bandVL_140_246(fatigue(i,j):fatigue(i,j+1),i).', 1);
                            xFitting = fatigue(i,j):fatigue(i,j+1); 
                            yFitted = polyval(coef, xFitting);
                            coeffs_VL_140_246(i,j) = yFitted(end)/yFitted(1)*100-100;
                            R = corr2(bandVL_140_246(fatigue(i,j):fatigue(i,j+1),i),yFitted.');
                            R2_VL_140_246(i,j) = R^2;
                            plot(xFitting, yFitted, '--','color', band4);
                        end            
                    end
                    % Plot straight lines when perceived fatigue changes
                    for j = 2:size(fatigue,2)
                        hold on
                        plot([fatigue(i,j) fatigue(i,j)],[0 lim_y],'k')
                    end    
                    % Save figure as JPG
                    saveas(gcf,sprintf([subjectID '_Bands_VL_%d.jpg'],i));

                % VM
                figure
                maxfig(gcf,1)
                    % Plot the intensity of each frequency band over time
                    if raw == 1
                        plot(bandVM_7_25(:,i),'color',band1)
                        hold on
                        plot(bandVM_35_45(:,i),'color',band2)
                        plot(bandVM_55_95(:,i),'color',band3)
                        plot(bandVM_140_246(:,i),'color',band4)
                    end
                    % Title and axes
                    title('VM','FontSize',30);
                    xlabel('Time (mins)');
                    xlim([0 lim_x])
                    ylim([0 lim_y])     
                    set(gcf,'numbertitle','off','name',sprintf('VM - Condition %d',i));
                    % Calculate trend lines for each fatigue section and each
                    % frequency band
                    if trend == 1
                        % 7 to 25 Hz
                        for j = 1:size(fatigue,2)-1
                            coef = polyfit(fatigue(i,j):fatigue(i,j+1), ...
                                bandVM_7_25(fatigue(i,j):fatigue(i,j+1),i).', 1);
                            xFitting = fatigue(i,j):fatigue(i,j+1); 
                            yFitted = polyval(coef, xFitting);
                            coeffs_VM_7_25(i,j) = yFitted(end)/yFitted(1)*100-100;
                            R = corr2(bandVM_7_25(fatigue(i,j):fatigue(i,j+1),i),yFitted.');
                            R2_VM_7_25(i,j) = R^2;
                            plot(xFitting, yFitted, '--','color', band1);
                        end
                        % 35 to 45 Hz
                        for j = 1:size(fatigue,2)-1
                            coef = polyfit(fatigue(i,j):fatigue(i,j+1), ...
                                bandVM_35_45(fatigue(i,j):fatigue(i,j+1),i).', 1);
                            xFitting = fatigue(i,j):fatigue(i,j+1); 
                            yFitted = polyval(coef, xFitting);
                            coeffs_VM_35_45(i,j) = yFitted(end)/yFitted(1)*100-100;
                            R = corr2(bandVM_35_45(fatigue(i,j):fatigue(i,j+1),i),yFitted.');
                            R2_VM_35_45(i,j) = R^2;
                            plot(xFitting, yFitted, '--','color', band2);
                        end
                        % 55 to 95 Hz
                        for j = 1:size(fatigue,2)-1
                            coef = polyfit(fatigue(i,j):fatigue(i,j+1), ...
                                bandVM_55_95(fatigue(i,j):fatigue(i,j+1),i).', 1);
                            xFitting = fatigue(i,j):fatigue(i,j+1); 
                            yFitted = polyval(coef, xFitting);
                            coeffs_VM_55_95(i,j) = yFitted(end)/yFitted(1)*100-100;
                            R = corr2(bandVM_55_95(fatigue(i,j):fatigue(i,j+1),i),yFitted.');
                            R2_VM_55_95(i,j) = R^2;
                            plot(xFitting, yFitted, '--','color', band3);
                        end
                        % 140 to 246 Hz
                        for j = 1:size(fatigue,2)-1
                            coef = polyfit(fatigue(i,j):fatigue(i,j+1), ...
                                bandVM_140_246(fatigue(i,j):fatigue(i,j+1),i).', 1);
                            xFitting = fatigue(i,j):fatigue(i,j+1); 
                            yFitted = polyval(coef, xFitting);
                            coeffs_VM_140_246(i,j) = yFitted(end)/yFitted(1)*100-100;
                            R = corr2(bandVM_140_246(fatigue(i,j):fatigue(i,j+1),i),yFitted.');
                            R2_VM_140_246(i,j) = R^2;
                            plot(xFitting, yFitted, '--','color', band4);
                        end            
                    end                
                    % Plot straight lines when perceived fatigue changes
                    for j = 2:size(fatigue,2)
                        hold on
                        plot([fatigue(i,j) fatigue(i,j)],[0 lim_y],'k')
                    end    
                    % Save figure as JPG
                    saveas(gcf,sprintf([subjectID '_Bands_VM_%d.jpg'],i));

    end

    % Combine coeffs in one matrix
    coeffs_GL = [coeffs_GL_7_25 coeffs_GL_35_45 coeffs_GL_55_95 coeffs_GL_140_246];
    coeffs_GM = [coeffs_GM_7_25 coeffs_GM_35_45 coeffs_GM_55_95 coeffs_GM_140_246];
    coeffs_VL = [coeffs_VL_7_25 coeffs_VL_35_45 coeffs_VL_55_95 coeffs_VL_140_246];
    coeffs_VM = [coeffs_VM_7_25 coeffs_VM_35_45 coeffs_VM_55_95 coeffs_VM_140_246];
    coeffs = [coeffs_GL; coeffs_GM; coeffs_VL; coeffs_VM];
    % Combine R-squared in one matrix
    R2_GL = [R2_GL_7_25 R2_GL_35_45 R2_GL_55_95 R2_GL_140_246];
    R2_GM = [R2_GM_7_25 R2_GM_35_45 R2_GM_55_95 R2_GM_140_246];
    R2_VL = [R2_VL_7_25 R2_VL_35_45 R2_VL_55_95 R2_VL_140_246];
    R2_VM = [R2_VM_7_25 R2_VM_35_45 R2_VM_55_95 R2_VM_140_246];
    R2 = [R2_GL; R2_GM; R2_VL; R2_VM];
    % Export results
    csvwrite([subjectID '_coeffs.csv'],coeffs);
    csvwrite([subjectID '_R2.csv'],R2);
        
    disp('     *****     DONE WITH FREQUENCY BANDS ANALYSIS     *****     ')
    