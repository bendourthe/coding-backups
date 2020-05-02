clear all; close all;
%% DESCRIPTION:
%
% This code was designed to enable the processing of the ACTIVE data
% collected for the Per4ma project (manually)
%
%   Project Title: The effects of X-Act Compression garments technology on
%   muscle fatigue, active and passive joint torque and performance.
%
%   Author: Benjamin Dourthe
%
%   Last update: Feb. 9th, 2018
%
%% Input: 
%   - Biodex active data
%
%% Output:
%   - Maximal torque for each trial
%   - Angular velocity for each trial
%
%% Dependencies:
%   Functions:
%       - WaveletFilter.m


%% CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % SHORTS or SHIRT
    cond = 2;       % Type 1 for SHORTS, 2 for SHIRT
    
    % 0deg, YES or NO
    deg = 0;        % Type 0 if you want to look at the 0deg condition, Type 1 otherwise
    
    % Signal sign
    s = 1;          % Type -1 if signal appears upside down and an error shows up, Type 1 otherwise
    

%% ID and DIRECTORY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Subject ID
    ID = 42;     % Type subject ID number (2 digits after PX)
    if ID < 10
        subjectID = sprintf('PX0%d',(ID));
    else
        subjectID = sprintf('PX%d',(ID));
    end

    % Current directory (where the data are stored)
    if cond == 1        
        cd 'C:\Users\bdour\OneDrive\Work\Calgary\Projects\Per4ma - X-ACT\Data Collection\Data\SHORTS\Valid\'
        directory = 'C:\Users\bdour\OneDrive\Work\Calgary\Projects\Per4ma - X-ACT\Data Collection\Data\SHORTS\Valid\';
    else
        cd 'C:\Users\bdour\OneDrive\Work\Calgary\Projects\Per4ma - X-ACT\Data Collection\Data\SHIRT\Valid\'
        directory = 'C:\Users\bdour\OneDrive\Work\Calgary\Projects\Per4ma - X-ACT\Data Collection\Data\SHIRT\Valid\';
    end

    % Data folder
    datafolder = 'Active\';

    % Labels
    act = '_act_';
    label = 'XA';
    deg0 = '_0';
    degX = '_30';


%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Signal sampling frequency (in Hz)
    samplingFrequency = 1000;
    
    % Set up Filter
    cf = 20; %center frequency
    mode = 6; % essentially the scale (determines the frequency bandwidth
    type = 1; %type of filter 1(low pass) 4(high pass)

    scale=1;
    

%% CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Files directories
    filedir = [directory subjectID '\' datafolder label '\' subjectID act label];
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if deg == 0
        % 0 deg/s
        deg0_file = [filedir deg0 '.csv'];
        deg0_data4 = csvread([deg0_file],2,1);
        
        le=size(deg0_data4,1);
        dt=1/samplingFrequency;
        T=le*dt;
        df=1/T;
        frequency=(0:le/2)*df;
        [deg0_data2,bb]=WaveletFilter(deg0_data4',dt,cf,mode,type); % low pass
        deg0_data3 = deg0_data2';
        
        if cond == 1
            deg0_data(:,1) = ((deg0_data3(:,1)*144.72)+2.3368);
            deg0_data(:,2) = ((deg0_data3(:,2)*34.634)+90);
            diffangdeg0 = diff(deg0_data(:,2))/0.001; % Superfluous - do not need
        else
            deg0_data(:,1) = s*(((deg0_data3(:,1)*144.72)+2.3368)*-1);
            deg0_data(:,2) = s*(((deg0_data3(:,2)*34.634)+90)*-1);
            diffangdeg0 = diff(deg0_data(:,2))/0.001; % Superfluous - do not need
        end
        torque0 = deg0_data(:,1);
        position0 = deg0_data(:,2);
        
        % Pick out the ramp and determine peak torque for each trial
        figure
        clf
        hold on
        plot(deg0_data)
        disp('CLICK around the beginning and ending of each ramp (from left to right)')
        x = round(ginput(4));
        starting0 = [x(1) x(2)];
        ending0 = [x(3) x(4)];
        r = [starting0;ending0];
        clf
        figure
        clf
        hold on
        plot(deg0_data)
        plot(r(1,:),torque0(r(1,:)),'b*','MarkerSize',10);
        plot(r(2,:),torque0(r(2,:)),'b*','MarkerSize',10);
        disp('PRESS ENTER if ramps are correctly located - PRESS Ctr+C in the command window if not')
        pause
        clf
        tor0(1) = max(torque0(r(1,1):r(2,1)));
        vel0(1) = 0;
        
        tor0(2) = max(torque0(r(1,2):r(2,2)));
        vel0(2) = 0;
        clear r
        clear deg0_data
        
        disp('Completed 0deg/s')
        tor0

else        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % X deg/s
        degX_file = [filedir degX '.csv'];        
        degX_data4 = csvread([degX_file],2,1);
        
        le=size(degX_data4,1);
        dt=1/samplingFrequency;
        T=le*dt;
        df=1/T;
        frequency=(0:le/2)*df;
        [degX_data2,bb]=WaveletFilter(degX_data4',dt,cf,mode,type);
        degX_data3 = degX_data2';
        
        if cond == 1
            degX_data(:,1) = ((degX_data3(:,1)*144.72)+2.3368);
            degX_data(:,2) = ((degX_data3(:,2)*34.634)+90);
            diffangdegX = diff(degX_data(:,2))./.001;
        else
            degX_data(:,1) = s*(((degX_data3(:,1)*144.72)+2.3368)*-1);
            degX_data(:,2) = s*(((degX_data3(:,2)*34.634)+90)*-1);
            diffangdegX = diff(degX_data(:,2))./.001;
        end
        torX = degX_data(:,1);
        positionX = degX_data(:,2);
        
        % Detects ramps and determine peak torque for each trial
        figure
        clf
        hold on
        plot(degX_data)
        plot(diffangdegX,'g')
        disp('CLICK around the beginning and ending of each ramp (from left to right)')
        x = round(ginput(4));
        starting0 = [x(1) x(2)];
        ending0 = [x(3) x(4)];
        clf
        figure
        clf
        hold on
        plot(degX_data)
        plot(diffangdegX,'g')
        r = [starting0;ending0];
        plot(r(1,:),positionX(r(1,:)),'b*','MarkerSize',10);
        plot(r(2,:),positionX(r(2,:)),'b*','MarkerSize',10);
        disp('PRESS ENTER if ramps are correctly located - PRESS Ctr+C in the command window if not')
        pause
        clf
        torX(1) = max(torX(r(1,1):r(2,1)));
        velX(1) = mean(diffangdegX(r(1,1):r(2,1)));
        
        torX(2) = max(torX(r(1,2):r(2,2)));
        velX(2) = mean(diffangdegX(r(1,2):r(2,2)));
        clear r
        clear degX_data
        
        disp('Completed Xdeg/s')
        torX
end
        