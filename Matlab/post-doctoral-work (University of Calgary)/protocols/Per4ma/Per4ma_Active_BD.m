clear all; close all;
%% DESCRIPTION:
%
% This code was designed to enable the processing of the ACTIVE data
% collected for the Per4ma project
%
%   Project Title: The effects of X-Act Compression garments technology on
%   muscle fatigue, active and passive joint torque and performance.
%
%   Author: Benjamin Dourthe
%
%   Last update: Feb. 8th, 2018
%
%% Input: 
%   - Biodex active data (7 CTRL, 7 UA + 1 CTRL30, 7 XA + 1 CTRL30)
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
    
    % Signal sign
    s = 1;          % Type -1 if signal appears upside down and an error shows up, Type 1 otherwise
    
    % Control 30 trial
    ctrl30 = 1;     % Type 0 if ctrl30 trials were not performed for the UA and XA conditions, Type 1 otherwise
    

%% ID and DIRECTORY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Subject ID
    ID = 41;     % Type subject ID number (2 digits after PX)
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
    ctrl = 'CTRL';
    ua = 'UA';
    xa = 'XA';
    CTRL30deg = '_CTRL30';
    deg0 = '_0';
    deg30 = '_30';
    deg60 = '_60';
    deg90 = '_90';
    deg120 = '_120';
    deg300 = '_300';
    deg500 = '_500';


%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Signal sampling frequency (in Hz)
    samplingFrequency = 1000;
    
    % Ramp detection
    numRamps = 2;
    pMax = 0.65;
    pMin = 0.30;
    ylim1 = -300;
    ylim2 = 200;

    % Set up Filter
    cf = 20; %center frequency
    mode = 6; % essentially the scale (determines the frequency bandwidth
    type = 1; %type of filter 1(low pass) 4(high pass)

    scale=1;
    

%% CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for d = 1:3

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Files directories

            % CTRL
            if d == 1
                filedir = [directory subjectID '\' datafolder ctrl '\' subjectID act ctrl];

            % UA    
            elseif d == 2
                filedir = [directory subjectID '\' datafolder ua '\' subjectID act ua];
                
            % XA
            elseif d == 3
                filedir = [directory subjectID '\' datafolder xa '\' subjectID act xa];
                
            end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        figure(1)
        clf
        hold on
        plot(deg0_data)
        
        for i=1:numRamps
            % Finds all data points that are above 75% of the max of the signal
            % (should include all data points that belong to the ramps)
            if i==1
                x1 = find(torque0 > pMin*max(torque0));
                starting0(i) = x1(1);
                x2 = find(torque0(starting0(i):end) < pMin*max(torque0));
                ending0(i) = starting0(i) + x2(1);
            else
                x1 = find(torque0(ending0(i-1):end) > pMin*max(torque0));
                starting0(i) = ending0(i-1) + x1(1);
                x2 = find(torque0(starting0(i):end) < pMin*max(torque0));
                ending0(i) = starting0(i) + x2(1);
            end
        end
        r = [starting0;ending0];
        plot(r(1,:),torque0(r(1,:)),'b*','MarkerSize',10);
        plot(r(2,:),torque0(r(2,:)),'b*','MarkerSize',10);
        ylim([ylim1 ylim2]);
        disp('PRESS ENTER if ramps are correctly located - PRESS Ctr+C in the command window if not')
        pause
        
        tor(d,1) = max(torque0(r(1,1):r(2,1)));
        vel(d,1) = 0;
        
        tor(d,2) = max(torque0(r(1,2):r(2,2)));
        vel(d,2) = 0;
        clear r
        clear deg0_data
        
        disp('Completed 0deg/s')

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 30 deg/s
        deg30_file = [filedir deg30 '.csv'];        
        deg30_data4 = csvread([deg30_file],2,1);
        
        le=size(deg30_data4,1);
        dt=1/samplingFrequency;
        T=le*dt;
        df=1/T;
        frequency=(0:le/2)*df;
        [deg30_data2,bb]=WaveletFilter(deg30_data4',dt,cf,mode,type);
        deg30_data3 = deg30_data2';
        
        if cond == 1
            deg30_data(:,1) = ((deg30_data3(:,1)*144.72)+2.3368);
            deg30_data(:,2) = ((deg30_data3(:,2)*34.634)+90);
            diffangdeg30 = diff(deg30_data(:,2))./.001;
        else
            deg30_data(:,1) = s*(((deg30_data3(:,1)*144.72)+2.3368)*-1);
            deg30_data(:,2) = s*(((deg30_data3(:,2)*34.634)+90)*-1);
            diffangdeg30 = diff(deg30_data(:,2))./.001;
        end
        torque30 = deg30_data(:,1);
        position30 = deg30_data(:,2);
        
        % Detects ramps and determine peak torque for each trial
        figure(1)
        clf
        hold on
        plot(deg30_data)
        plot(diffangdeg30,'g')
        r = ramp_det(position30,samplingFrequency,numRamps,pMax);
        plot(r(1,:),position30(r(1,:)),'b*','MarkerSize',10);
        plot(r(2,:),position30(r(2,:)),'b*','MarkerSize',10);
        ylim([ylim1 ylim2]);
        disp('PRESS ENTER if ramps are correctly located - PRESS Ctr+C in the command window if not')
        pause
        
        tor(d,3) = max(torque30(r(1,1):r(2,1)));
        vel(d,3) = mean(diffangdeg30(r(1,1):r(2,1)));
        
        tor(d,4) = max(torque30(r(1,2):r(2,2)));
        vel(d,4) = mean(diffangdeg30(r(1,2):r(2,2)));
        clear r
        clear deg30_data
        
        disp('Completed 30deg/s')
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 60 deg/s
        deg60_file = [filedir deg60 '.csv'];
        deg60_data4 = csvread([deg60_file],2,1);
        
        le=size(deg60_data4,1);
        dt=1/samplingFrequency;
        T=le*dt;
        df=1/T;
        frequency=(0:le/2)*df;
        [deg60_data2,bb]=WaveletFilter(deg60_data4',dt,cf,mode,type);
        deg60_data3 = deg60_data2';
        
        if cond == 1
        deg60_data(:,1) = ((deg60_data3(:,1)*144.72)+2.3368);
        deg60_data(:,2) = ((deg60_data3(:,2)*34.634)+90);
        diffangdeg60 = diff(deg60_data(:,2))/.001;
        else
        deg60_data(:,1) = s*(((deg60_data3(:,1)*144.72)+2.3368)*-1);
        deg60_data(:,2) = s*(((deg60_data3(:,2)*34.634)+90)*-1);
        diffangdeg60 = diff(deg60_data(:,2))/.001;
        end
        torque60 = deg60_data(:,1);
        position60 = deg60_data(:,2);
        
        % Pick out the ramp and determine peak torque for each trial
        figure(1)
        clf
        hold on
        plot(deg60_data)
        plot(diffangdeg60,'g')
        r = ramp_det(position60,samplingFrequency,numRamps,pMax);
        plot(r(1,:),position60(r(1,:)),'b*','MarkerSize',10);
        plot(r(2,:),position60(r(2,:)),'b*','MarkerSize',10);
        ylim([ylim1 ylim2]);
        disp('PRESS ENTER if ramps are correctly located - PRESS Ctr+C in the command window if not')
        pause
        
        tor(d,5) = max(torque60(r(1,1):r(2,1)));
        vel(d,5) = mean(diffangdeg60(r(1,1):r(2,1)));
        
        tor(d,6) = max(torque60(r(1,2):r(2,2)));
        vel(d,6) = mean(diffangdeg60(r(1,2):r(2,2)));
        clear r
        clear deg60_data
        
        disp('Completed 60deg/s')
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 90 deg/s
        deg90_file = [filedir deg90 '.csv'];
        deg90_data4 = csvread([deg90_file],2,1);
        
        le=size(deg90_data4,1);
        dt=1/samplingFrequency;
        T=le*dt;
        df=1/T;
        frequency=(0:le/2)*df;
        [deg90_data2,bb]=WaveletFilter(deg90_data4',dt,cf,mode,type);
        deg90_data3 = deg90_data2';
        
        if cond == 1
            deg90_data(:,1) = ((deg90_data3(:,1)*144.72)+2.3368);
            deg90_data(:,2) = ((deg90_data3(:,2)*34.634)+90);
            diffangdeg90 = diff(deg90_data(:,2))/.001;
        else
            deg90_data(:,1) = s*(((deg90_data3(:,1)*144.72)+2.3368)*-1);
            deg90_data(:,2) = s*(((deg90_data3(:,2)*34.634)+90)*-1);
            diffangdeg90 = diff(deg90_data(:,2))/.001;
        end
        torque90 = deg90_data(:,1);
        position90 = deg90_data(:,2);
        
        % Pick out the ramp and determine peak torque for each trial
        figure(1)
        clf
        hold on
        plot(deg90_data)
        plot(diffangdeg90,'g')
        r = ramp_det(position90,samplingFrequency,numRamps,pMax);
        plot(r(1,:),position90(r(1,:)),'b*','MarkerSize',10);
        plot(r(2,:),position90(r(2,:)),'b*','MarkerSize',10);
        ylim([ylim1 ylim2]);
        disp('PRESS ENTER if ramps are correctly located - PRESS Ctr+C in the command window if not')
        pause
        
        tor(d,7) = max(torque90(r(1,1):r(2,1)));
        vel(d,7) = mean(diffangdeg90(r(1,1):r(2,1)));
        
        tor(d,8) = max(torque90(r(1,2):r(2,2)));
        vel(d,8) = mean(diffangdeg90(r(1,2):r(2,2)));
        clear r
        clear deg90_data
        
        disp('Completed 90deg/s')
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 120 deg/s
        deg120_file = [filedir deg120 '.csv'];
        deg120_data4 = csvread([deg120_file],2,1);
        
        le=size(deg120_data4,1);
        dt=1/samplingFrequency;
        T=le*dt;
        df=1/T;
        frequency=(0:le/2)*df;
        [deg120_data2,bb]=WaveletFilter(deg120_data4',dt,cf,mode,type);
        deg120_data3 = deg120_data2';
        
        if cond == 1
            deg120_data(:,1) = ((deg120_data3(:,1)*144.72)+2.3368);
            deg120_data(:,2) = ((deg120_data3(:,2)*34.634)+90);
            diffangdeg120 = diff(deg120_data(:,2))/.001;
        else
            deg120_data(:,1) = s*(((deg120_data3(:,1)*144.72)+2.3368)*-1);
            deg120_data(:,2) = s*(((deg120_data3(:,2)*34.634)+90)*-1);
            diffangdeg120 = diff(deg120_data(:,2))/.001;
        end
        torque120 = deg120_data(:,1);
        position120 = deg120_data(:,2);
        
        % Pick out the ramp and determine peak torque for each trial
        figure(1)
        clf
        hold on
        plot(deg120_data)
        plot(diffangdeg120,'g')
        r = ramp_det(position120,samplingFrequency,numRamps,pMax);
        plot(r(1,:),position120(r(1,:)),'b*','MarkerSize',10);
        plot(r(2,:),position120(r(2,:)),'b*','MarkerSize',10);
        ylim([ylim1 ylim2]);
        disp('PRESS ENTER if ramps are correctly located - PRESS Ctr+C in the command window if not')
        pause
        
        tor(d,9) = max(torque120(r(1,1):r(2,1)));
        vel(d,9) = mean(diffangdeg120(r(1,1):r(2,1)));
        
        tor(d,10) = max(torque120(r(1,2):r(2,2)));
        vel(d,10) = mean(diffangdeg120(r(1,2):r(2,2)));
        clear r
        clear deg120_data
        
        disp('Completed 120deg/s')
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %300 deg/s
        deg300_file = [filedir deg300 '.csv'];
        deg300_data4 = csvread([deg300_file],2,1);
        
        le=size(deg300_data4,1);
        dt=1/samplingFrequency;
        T=le*dt;
        df=1/T;
        frequency=(0:le/2)*df;
        cf = 45; %center frequency
        [deg300_data2,bb]=WaveletFilter(deg300_data4',dt,cf,mode,type);
        deg300_data3 = deg300_data2';
        
        if cond == 1
            deg300_data(:,1) = ((deg300_data3(:,1)*144.72)+2.3368);
            deg300_data(:,2) = ((deg300_data3(:,2)*34.634)+90);
            diffangdeg300 = diff(deg300_data(:,2))/.001;
        else
            deg300_data(:,1) = s*(((deg300_data3(:,1)*144.72)+2.3368)*-1);
            deg300_data(:,2) = s*(((deg300_data3(:,2)*34.634)+90)*-1);
            diffangdeg300 = diff(deg300_data(:,2))/.001;
        end
        torque300 = deg300_data(:,1);
        position300 = deg300_data(:,2);
        
        % Pick out the ramp and determine peak torque for each trial
        figure(1)
        clf
        hold on
        plot(deg300_data)
        plot(diffangdeg300,'g')
        r = ramp_det(position300,samplingFrequency,numRamps,pMax);
        plot(r(1,:),position300(r(1,:)),'b*','MarkerSize',10);
        plot(r(2,:),position300(r(2,:)),'b*','MarkerSize',10);
        ylim([ylim1 ylim2]);
        disp('PRESS ENTER if ramps are correctly located - PRESS Ctr+C in the command window if not')
        pause
        
        tor(d,11) = max(torque300(r(1,1):r(2,1)));
        vel(d,11) = mean(diffangdeg300(r(1,1):r(2,1)));
        
        tor(d,12) = max(torque300(r(1,2):r(2,2)));
        vel(d,12) = mean(diffangdeg300(r(1,2):r(2,2)));
        clear r
        clear deg300_data
        
        disp('Completed 300deg/s')
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 500 deg/s
        deg500_file = [filedir deg500 '.csv'];
        deg500_data4 = csvread([deg500_file],2,1);
        
        le=size(deg500_data4,1);
        dt=1/samplingFrequency;
        T=le*dt;
        df=1/T;
        frequency=(0:le/2)*df;
        cf = 70; %center frequency
        [deg500_data2,bb]=WaveletFilter(deg500_data4',dt,cf,mode,type);
        deg500_data3 = deg500_data2';
        
        if cond == 1
            deg500_data(:,1) = ((deg500_data3(:,1)*144.72)+2.3368);
            deg500_data(:,2) = ((deg500_data3(:,2)*34.634)+90);
            diffangdeg500 = diff(deg500_data(:,2))/.001;
        else
            deg500_data(:,1) = s*(((deg500_data3(:,1)*144.72)+2.3368)*-1);
            deg500_data(:,2) = s*(((deg500_data3(:,2)*34.634)+90)*-1);
            diffangdeg500 = diff(deg500_data(:,2))/.001;
        end
        torque500 = deg500_data(:,1);
        position500 = deg500_data(:,2);
        
        % Pick out the ramp and determine peak torque for each trial
        figure(1)
        clf
        hold on
        plot(deg500_data)
        plot(diffangdeg500,'g')
        r = ramp_det(position500,samplingFrequency,numRamps,pMax);
        plot(r(1,:),position500(r(1,:)),'b*','MarkerSize',10);
        plot(r(2,:),position500(r(2,:)),'b*','MarkerSize',10);
        ylim([ylim1 ylim2]);
        disp('PRESS ENTER if ramps are correctly located - PRESS Ctr+C in the command window if not')
        pause
        
        tor(d,13) = max(torque500(r(1,1):r(2,1)));
        vel(d,13) = mean(diffangdeg500(r(1,1):r(2,1)));
        
        tor(d,14) = max(torque500(r(1,2):r(2,2)));
        vel(d,14) = mean(diffangdeg500(r(1,2):r(2,2)));
        clear r
        clear deg500_data
        
        disp('Completed 500deg/s')
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % CTRL 30 deg/s       
    
        if d == 2
            
            if ctrl30 == 1
            
            % CTRL30 deg/s
            CTRL30deg_file = [filedir CTRL30deg '.csv'];
            CTRL30deg_data4 = csvread([CTRL30deg_file],2,1);

            le=size(CTRL30deg_data4,1);
            dt=1/samplingFrequency;
            T=le*dt;
            df=1/T;
            frequency=(0:le/2)*df;
            [CTRL30deg_data2,bb]=WaveletFilter(CTRL30deg_data4',dt,cf,mode,type);
            CTRL30deg_data3 = CTRL30deg_data2';

            if cond == 1
            CTRL30deg_data(:,1) = ((CTRL30deg_data3(:,1)*144.72)+2.3368);
            CTRL30deg_data(:,2) = ((CTRL30deg_data3(:,2)*34.634)+90);
            diffangUACTRL30 = diff(CTRL30deg_data(:,2))/0.001;
            else
            CTRL30deg_data(:,1) = s*(((CTRL30deg_data3(:,1)*144.72)+2.3368)*-1);
            CTRL30deg_data(:,2) = s*(((CTRL30deg_data3(:,2)*34.634)+90)*-1);
            diffangUACTRL30 = diff(CTRL30deg_data(:,2))/0.001;
            end
            torqueUACTRL30 = CTRL30deg_data(:,1);
            positionUACTRL30 = CTRL30deg_data(:,2);

            % Pick out the ramp and determine peak torque for each trial
            figure(1)
            clf
            hold on
            plot(CTRL30deg_data)
            plot(diffangUACTRL30,'g')
            r = ramp_det(positionUACTRL30,samplingFrequency,numRamps,pMax);
            plot(r(1,:),positionUACTRL30(r(1,:)),'b*','MarkerSize',10);
            plot(r(2,:),positionUACTRL30(r(2,:)),'b*','MarkerSize',10);
            ylim([ylim1 ylim2]);
            disp('PRESS ENTER if ramps are correctly located - PRESS Ctr+C in the command window if not')
            pause

            torCTRL30(1,1) = max(torqueUACTRL30(r(1,1):r(2,1)));
            velCTRL30(1,1) = mean(diffangUACTRL30(r(1,1):r(2,1)));

            torCTRL30(1,2) = max(torqueUACTRL30(r(1,2):r(2,2)));
            velCTRL30(1,2) = mean(diffangUACTRL30(r(1,2):r(2,2)));
            
            else
            torCTRL30(1,1) = 0;
            velCTRL30(1,1) = 0;

            torCTRL30(1,2) = 0;
            velCTRL30(1,2) = 0;
            end
            
            clear r
            clear CTRL30deg_data

            disp('Completed CTRL30')

        elseif d == 3
            
            if ctrl30 == 1

            % CTRL30 deg/s
            CTRL30deg_file = [filedir CTRL30deg '.csv'];
            CTRL30deg_data4 = csvread([CTRL30deg_file],2,1);

            le=size(CTRL30deg_data4,1);
            dt=1/samplingFrequency;
            T=le*dt;
            df=1/T;
            frequency=(0:le/2)*df;
            [CTRL30deg_data2,bb]=WaveletFilter(CTRL30deg_data4',dt,cf,mode,type);
            CTRL30deg_data3 = CTRL30deg_data2';

            if cond == 1
            CTRL30deg_data(:,1) = ((CTRL30deg_data3(:,1)*144.72)+2.3368);
            CTRL30deg_data(:,2) = ((CTRL30deg_data3(:,2)*34.634)+90);
            diffangXACTRL30 = diff(CTRL30deg_data(:,2))/0.001;
            else
            CTRL30deg_data(:,1) = s*(((CTRL30deg_data3(:,1)*144.72)+2.3368)*-1);
            CTRL30deg_data(:,2) = s*(((CTRL30deg_data3(:,2)*34.634)+90)*-1);
            diffangXACTRL30 = diff(CTRL30deg_data(:,2))/0.001;
            end
            torqueXACTRL30 = CTRL30deg_data(:,1);
            positionXACTRL30 = CTRL30deg_data(:,2);            

            % Pick out the ramp and determine peak torque for each trial
            figure(1)
            clf
            hold on
            plot(CTRL30deg_data)
            plot(diffangXACTRL30,'g')
            r = ramp_det(positionXACTRL30,samplingFrequency,numRamps,pMax);
            plot(r(1,:),positionXACTRL30(r(1,:)),'b*','MarkerSize',10);
            plot(r(2,:),positionXACTRL30(r(2,:)),'b*','MarkerSize',10);
            ylim([ylim1 ylim2]);
            disp('PRESS ENTER if ramps are correctly located - PRESS Ctr+C in the command window if not')
            pause

            torCTRL30(2,1) = max(torqueXACTRL30(r(1,1):r(2,1)));
            velCTRL30(2,1) = mean(diffangXACTRL30(r(1,1):r(2,1)));

            torCTRL30(2,2) = max(torqueXACTRL30(r(1,2):r(2,2)));
            velCTRL30(2,2) = mean(diffangXACTRL30(r(1,2):r(2,2)));
            
            else
            torCTRL30(2,1) = 0;
            velCTRL30(2,1) = 0;

            torCTRL30(2,2) = 0;
            velCTRL30(2,2) = 0;
            end
            
            clear r
            clear CTRL30deg_data
            
            disp('Completed CTRL30')

        end
        
end

Results = [tor(1,:) tor(2,:) torCTRL30(1,:) tor(3,:) torCTRL30(2,:);...
            vel(1,:) vel(2,:) velCTRL30(1,:) vel(3,:) velCTRL30(2,:)];
csvwrite(['Results_' subjectID '_act' '.csv'], Results);
    
close all

        