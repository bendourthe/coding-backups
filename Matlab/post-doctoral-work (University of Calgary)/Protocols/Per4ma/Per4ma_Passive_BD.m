clear all; close all;
%% DESCRIPTION:
%
% This code was designed to enable the processing of the PASSIVE data
% collected for the Per4ma project
%
%   Project Title: The effects of X-Act Compression garments technology on
%   muscle fatigue, active and passive joint torque and performance.
%
%   Author: Benjamin Dourthe
%
%   Last update: Jan. 26th, 2018
%
%% Input: 
%   - Biodex passive data (4 CTRL, 2 UA, 2 XA)
%
%% Output:
%   - Passive Range of Motion (ROM)
%   - Peak and baseline torques for each garment (CTRL, UA, XA)
%
%% Dependencies:
%   Functions:
%       - WaveletFilter.m


%% SHORTS or SHIRT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Condition
    cond = 1;       % Type 1 for SHORTS, 2 for SHIRT
    

%% ID and DIRECTORY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Subject ID
    ID = 34;     % Type subject ID number (2 digits after PX)
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
    datafolder = 'Passive\';

    % Label
    pass = '_pass_';


%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Crop factor
    cropf = 0.5;

    % Passive velocity
    vel = 10;   % degrees/s

    % Set up Filter
    cf = 30; % center frequency
    mode = 6; % essentially the scale (determines the frequency bandwidth)
    type = 1; % type of filter 1(low pass) 4(high pass)
    sampling_rate = 1000;
    scale = 1;
    

%% CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for s = ID

    % CONTROL TRIALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for a = 1:4

        % Control
        CTRL = sprintf('CTRL%d',(a));
        ctrlfile = [directory subjectID '\' datafolder subjectID pass CTRL '.csv'];
        CTRL_data2 = csvread([ctrlfile],2,1); % channel 1 is time, 2 is torque, 3 is position and 4 is EMG

        % Filter
        le=length(CTRL_data2);
        dt=1/sampling_rate;
        T=le*dt;
        df=1/T;
        frequency=(0:le/2)*df;
        [CTRL_data3,bb]=WaveletFilter(CTRL_data2',dt,cf,mode,type); %low pass
        CTRL_data = CTRL_data3';

        % Convert from voltage to Nm
        CTRL_data(:,1) = ((CTRL_data(:,1)*144.72)+2.3368);  % Torque
        CTRL_data(:,2) = ((CTRL_data(:,2)*34.634)+90);      % Position

        % Isolate torque from position
        if cond == 1 % SHORTS condition requires to crop the end of the signal to have only one full cylcle
            crop = round(length(CTRL_data)*cropf);
            torqueCTRL = CTRL_data((1:crop),1); 
            posCTRL = CTRL_data((1:crop),2);
        else
            torqueCTRL = CTRL_data(:,1); 
            posCTRL = CTRL_data(:,2);
        end

        % Passive ROM
        [MROM, IROM] = max(posCTRL);
        [mROM, iROM] = min(posCTRL);
        passROM(a) = vel*(iROM-IROM)/sampling_rate;

        % Identify beginning and end of motion cycle
        if iROM < IROM
            limit1 = iROM+500;
            limit2 = IROM;
        else
            limit1 = IROM+500;
            limit2 = iROM;
        end

        % Peak torque and baseline
        [MP, IP] = max(torqueCTRL(limit1:limit2));
        [mP, iP] = min(torqueCTRL(limit1:limit2));
        IP = limit1 + IP;
        iP = limit1 + iP;
        CTRL_peak(a) = MP;
        CTRL_base(a) = mP;

        % Plot figure to assess validity
        fig = figure;
        hold on
        plot(torqueCTRL)
        plot(posCTRL,'r')
        plot(IROM,MROM,'r*','MarkerSize',10)
        plot(iROM,mROM,'r*','MarkerSize',10)
        plot(IP,MP,'b*','MarkerSize',10)
        plot(iP,mP,'b*','MarkerSize',10)
        set(fig,'numbertitle','off','name','CTRL-PASSIVE: Evolution of position (in red) and torque (in blue)');
        disp('Check if detected min/max are correct, PRESS ENTER if valid, PRESS Ctrl+C in the Command Window if not')
        pause
        clf

    end
    disp('Completed CTRL trials')


    % COMPRESSIVE GARMENTS

    for b = 1:2

        % Under Armor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        UA = sprintf('UA%d',(b));
        uafile = [directory subjectID '\' datafolder subjectID pass UA '.csv'];
        UA_data2 = csvread([uafile],2,1); %channel 1 is time, 2 is torque, 3 is position and 4 is EMG

        % Filter
        le=length(UA_data2);
        dt=1/sampling_rate;
        T=le*dt;
        df=1/T;
        frequency=(0:le/2)*df;
        [UA_data3,bb]=WaveletFilter(UA_data2',dt,cf,mode,type); %low pass
        UA_data = UA_data3';

        % Convert from voltage to Nm
        UA_data(:,1) = ((UA_data(:,1)*144.72)+2.3368);
        UA_data(:,2) = ((UA_data(:,2)*34.634)+90);

        % Isolate torque from position
        if cond == 1 % SHORTS condition requires to crop the end of the signal to have only one full cylcle
            crop = round(length(UA_data)*cropf);
            torqueUA = UA_data((1:crop),1); 
            posUA = UA_data((1:crop),2);
        else
            torqueUA = UA_data(:,1); 
            posUA = UA_data(:,2);
        end

        % Passive ROM
        [MROM, IROM] = max(posUA);
        [mROM, iROM] = min(posUA);

        % Identify beginning and end of motion cycle
        if iROM < IROM
            limit1 = iROM+500;
            limit2 = IROM;
        else
            limit1 = IROM+500;
            limit2 = iROM;
        end

        % Peak torque and baseline
        [MP, IP] = max(torqueUA(limit1:limit2));
        [mP, iP] = min(torqueUA(limit1:limit2));
        IP = limit1 + IP;
        iP = limit1 + iP;
        UA_peak(b) = MP;
        UA_base(b) = mP;

        % Plot figure to assess validity
        fig = figure;
        hold on
        plot(torqueUA)
        plot(posUA,'r')
        plot(IROM,MROM,'r*','MarkerSize',10)
        plot(iROM,mROM,'r*','MarkerSize',10)
        plot(IP,MP,'b*','MarkerSize',10)
        plot(iP,mP,'b*','MarkerSize',10)
        set(fig,'numbertitle','off','name','UA-PASSIVE: Evolution of position (in red) and torque (in blue)');
        disp('Check if detected min/max are correct, PRESS ENTER if valid, PRESS Ctrl+C in the Command Window if not')
        pause
        clf


        % X-Act %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        XA = sprintf('XA%d',(b));
        xafile = [directory subjectID '\' datafolder subjectID pass XA '.csv'];
        XA_data2 = csvread([xafile],2,1); %channel 1 is time, 2 is torque, 3 is position and 4 is EMG

        % Filter
        le=length(XA_data2);
        dt=1/sampling_rate;
        T=le*dt;
        df=1/T;
        frequency=(0:le/2)*df;
        [XA_data3,bb]=WaveletFilter(XA_data2',dt,cf,mode,type); %low pass
        XA_data = XA_data3';

        % Convert from voltage to Nm
        XA_data(:,1) = ((XA_data(:,1)*144.72)+2.3368);
        XA_data(:,2) = ((XA_data(:,2)*34.634)+90);

        % Isolate torque from position
        if cond == 1 % SHORTS condition requires to crop the end of the signal to have only one full cylcle
            crop = round(length(XA_data)*cropf);
            torqueXA = XA_data((1:crop),1); 
            posXA = XA_data((1:crop),2);
        else
            torqueXA = XA_data(:,1); 
            posXA = XA_data(:,2);
        end

        % Passive ROM
        [MROM, IROM] = max(posXA);
        [mROM, iROM] = min(posXA);

        % Identify beginning and end of motion cycle
        if iROM < IROM
            limit1 = iROM+500;
            limit2 = IROM;
        else
            limit1 = IROM+500;
            limit2 = iROM;
        end

        % Peak torque and baseline
        [MP, IP] = max(torqueXA(limit1:limit2));
        [mP, iP] = min(torqueXA(limit1:limit2));
        IP = limit1 + IP;
        iP = limit1 + iP;
        XA_peak(b) = MP;
        XA_base(b) = mP;

        % Plot figure to assess validity
        fig = figure;
        hold on
        plot(torqueXA)
        plot(posXA,'r')
        plot(IROM,MROM,'r*','MarkerSize',10)
        plot(iROM,mROM,'r*','MarkerSize',10)
        plot(IP,MP,'b*','MarkerSize',10)
        plot(iP,mP,'b*','MarkerSize',10)
        set(fig,'numbertitle','off','name','XA-PASSIVE: Evolution of position (in red) and torque (in blue)');
        disp('Check if detected min/max are correct, PRESS ENTER if valid, PRESS Ctrl+C in the Command Window if not')
        pause
        clf

    end    
        disp('Completed UA and XA trials')     
end
close all

%% RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Copy/Paste this matrix directly in the Excel sheet
    Results = [abs(passROM) CTRL_base CTRL_peak UA_base UA_peak XA_base XA_peak];
    csvwrite(['Results_' subjectID '_pass' '.csv'], Results);
        