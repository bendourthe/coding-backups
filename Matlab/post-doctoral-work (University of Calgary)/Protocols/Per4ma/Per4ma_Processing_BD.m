%%%%%%%%%%%%%    PER4MA ANALYSIS    %%%%%%%%%%%%%%%%%%%%%%%%%%%  May , 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This program will load raw(.asc) and force(.fgt) PEDAR files to anaylize
% COP (ML & AP), peak force, and partitians the insole to different regions

% Created by Chris Lam

% Start Fresh2death
close all;clc;clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Allocate directory
directory = 'C:\Users\bdour\OneDrive\Work\Calgary\Projects\Per4ma - X-ACT\Data Collection\Data\SHORTS\';

% Save folder
savefolder = 'Results\';

% Data folder
datafolder = 'Valid\';

% Load Per4ma Matrix to add additional trials to it
% load('Per4madata.mat')


%%   PASSIVE ROM   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load('Per4madata.mat')

% Passive label
pass = '_pass_';

% Set up Filter
cf = 30; % center frequency
mode = 6; % essentially the scale (determines the frequency bandwidth)
type = 1; % type of filter 1(low pass) 4(high pass)
sampling_rate = 1000;
scale = 1;


%%% Subject Loop %%%

for s = 1
    
    subjectID = sprintf('PX0%d',(s));
        
    for a = 1:4
        
        % Control label
        CTRL = sprintf('CTRL%d',(a));
        ctrlfile = [directory datafolder subjectID '\' subjectID pass CTRL '.csv'];
        CTRL_data2 = csvread([ctrlfile],2,1); % channel 1 is time, 2 is torque, 3 is position and 4 is EMG
        
        % Filter
        le=size(CTRL_data2,1);
        dt=1/sampling_rate;
        T=le*dt;
        df=1/T;
        frequency=(0:le/2)*df;
        [CTRL_data3,bb]=WaveletFilter(CTRL_data2',dt,cf,mode,type); %low pass
        CTRL_data = CTRL_data3';
        
        % Convert from voltage to Nm
        CTRL_data(:,1) = ((CTRL_data(:,1)*144.72)+2.3368);
        CTRL_data(:,2) = ((CTRL_data(:,2)*34.634)+90);
        
        % Display figure to determine area to measure torque
        figure(2);
        clf
        hold on
        plot(CTRL_data)
        disp('Click graph around just before lowest valley (baseline), then around peak flexion')
        f = ginput(4);
        peakctrl = CTRL_data(f(3,1):f(4,1),1:2); 
        baseline2 = CTRL_data(f(1,1):f(2,1),1);
        
        
        CTRL_torque(s,a) = max(peakctrl(:,1));  %%FIND MIN FOR LEFTIES
        CTRL_baseline(s,a) = mean(baseline2);
        
        temp = peakctrl(:,2)+2.3;
        clf
        hold on
        plot(peakctrl(:,1))
        plot(temp,'r')
        disp('Just to make sure the value makes sense')
        disp(CTRL_torque(a))
        pause
        
    end
    
    disp('Completed CTRL trials')
    
    for b = 1:2
        % Under Armor
        UA = sprintf('UA%d',(b));
        uafile = [directory datafolder subjectID '\' subjectID pass UA '.csv'];
        UA_data2 = csvread([uafile],2,1); %channel 1 is time, 2 is torque, 3 is position and 4 is EMG
        
        % Filter
        le=size(UA_data2,1);
        dt=1/sampling_rate;
        T=le*dt;
        df=1/T;
        frequency=(0:le/2)*df;
        [UA_data3,bb]=WaveletFilter(UA_data2',dt,cf,mode,type); %low pass
        UA_data = UA_data3';
        
        % Convert from voltage to Nm
        UA_data(:,1) = ((UA_data(:,1)*144.72)+2.3368);
        UA_data(:,2) = ((UA_data(:,2)*34.634)+90);
        
        % Display figure to determine area to measure torque
        figure(1);
        clf
        hold on
        plot(UA_data)
        disp('Click graph around just before lowest valley (baseline), then around peak flexion')
        f = ginput(4);
        peakUA = UA_data(f(3,1):f(4,1),1:2); 
        baseline2 = UA_data(f(1,1):f(2,1),1);
        
        UA_torque(s,b) = max(peakUA(:,1));  %% FIND MIN FOR LEFTIES
        UA_baseline(s,b) = mean(baseline2);
        
        temp = peakUA(:,2)+2.3;
        clf
        hold on
        plot(peakUA(:,1))
        plot(temp,'r')
        disp('Just to make sure the value makes sense')
        disp(UA_torque(b))
        pause

        % X-Act
        XA = sprintf('XA%d',(b));
        xafile = [directory datafolder subjectID '\' subjectID pass XA '.csv'];
        XA_data2 = csvread([xafile],2,1); %channel 1 is time, 2 is torque, 3 is position and 4 is EMG
        
        % Filter
        le=size(XA_data2,1);
        dt=1/sampling_rate;
        T=le*dt;
        df=1/T;
        frequency=(0:le/2)*df;
        [XA_data3,bb]=WaveletFilter(XA_data2',dt,cf,mode,type); %low pass
        XA_data = XA_data3';
        
        % Convert from voltage to Nm
        XA_data(:,1) = ((XA_data(:,1)*144.72)+2.3368);
        XA_data(:,2) = ((XA_data(:,2)*34.634)+90);
        
        % Display figure to determine area to measure torque
        figure(1);
        clf
        hold on
        plot(XA_data)
        disp('Click graph around just before lowest valley (baseline), then around peak flexion')
        f = ginput(4);
        peakXA = XA_data(f(3,1):f(4,1),1:2); 
        baseline2 = XA_data(f(1,1):f(2,1),1);
        
        XA_torque(s,b) = max(peakXA(:,1)); %%FIND MIN FOR LEFTIES
        XA_baseline(s,b) = mean(baseline2);
        
        temp = peakXA(:,2)+2.3;
        clf
        hold on
        plot(peakXA(:,1))
        plot(temp,'r')
        disp('Just to make sure the value makes sense')
        disp(XA_torque(b))
        pause
        
    end    
        disp('Completed UA and XA trials')

end

% eval(['Per4madata.passive.CTRL.peaktorq = CTRL_torque;']);
% eval(['Per4madata.passive.CTRL.baseline = CTRL_baseline;']);
% eval(['Per4madata.passive.UA.peaktorq = UA_torque;']);
% eval(['Per4madata.passive.UA.baseline = UA_baseline;']);
% eval(['Per4madata.passive.XA.peaktorq = XA_torque;']);
% eval(['Per4madata.passive.XA.baseline = XA_baseline;']);
% 
% save('Per4madata.mat','Per4madata')


%%   TORQUE-VELOCITY   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load('Per4madata.mat')  %NOTE Need to clear all before running this for shirts after running shorts and vice versa

% Active labels
act = '_Act_';
ctrl = 'CTRL';
ua = 'UA';
xa = 'XA';
% CTRL30deg = '_CTRL30';
deg0 = '0';
deg30 = '30';
deg60 = '60';
deg90 = '90';
deg120 = '120';
deg300 = '300';
deg500 = '500';

% Set up Filter
cf = 20; %center frequency
mode = 6; % essentially the scale (determines the frequency bandwidth
type = 1; %type of filter 1(low pass) 4(high pass)
sampling_rate = 1000;
scale=1;

close all;clc;clear all;

%%% Subject Loop %%%

for sss = 10:23
    
    subjectID = sprintf('PX%d',(sss));
    
    prompt2 = 'Press 0 for RIGHT ARM/LEFT LEG\nPress 1 for LEFT ARM/RIGHT LEG\n';
    z = input(prompt2)
    
    for d = 1:3
            
%         load('Per4madata.mat')  %NOTE Need to clear all before running this for shirts after running shorts and vice versa
        
        directory = 'C:\Users\bdour\OneDrive\Work\Calgary\Projects\Per4ma - X-ACT\Data Collection\Data\SHORTS\';
        savefolder = 'Results\';
        datafolder = 'Valid\';
%         load('Per4madata.mat')
        
        % Active labels
        act = '_act_';
        ctrl = 'CTRL';
        ua = 'UA';
        xa = 'XA';
        CTRL30deg = 'CTRL30';
        deg0 = '0';
        deg30 = '30';
        deg60 = '60';
        deg90 = '90';
        deg120 = '120';
        deg300 = '300';
        deg500 = '500';
        
        % Set up Filter
        cf = 20; % center frequency
        mode = 6; % essentially the scale (determines the frequency bandwidth
        type = 1; % type of filter 1(low pass) 4(high pass)
        sampling_rate = 1000;
        scale=1;        
        
        if d == 1
            filedir = [directory subjectID '\' datafolder ctrl '\' subjectID act ctrl];
            e = 1;
            cf = 20; % center frequency
            
        elseif d == 2
            filedir = [directory datafolder subjectID '\' ua '\' subjectID act ua];
            e = 3;
            cf = 20; % center frequency

            % CTRL30 deg/s
            CTRL30deg_file = [filedir CTRL30deg '.csv'];
            CTRL30deg_data4 = csvread([CTRL30deg_file],2,1);
            
            le=size(CTRL30deg_data4,1);
            dt=1/sampling_rate;
            T=le*dt;
            df=1/T;
            frequency=(0:le/2)*df;
            [CTRL30deg_data2,bb]=WaveletFilter(CTRL30deg_data4',dt,cf,mode,type);
            CTRL30deg_data3 = CTRL30deg_data2';
            
            if z == 1
            CTRL30deg_data(:,1) = ((CTRL30deg_data3(:,1)*144.72)+2.3368);
            CTRL30deg_data(:,2) = ((CTRL30deg_data3(:,2)*34.634)+90);
            diffangCTRL30deg = diff(CTRL30deg_data(:,2))/0.001;
            else
            CTRL30deg_data(:,1) = ((CTRL30deg_data3(:,1)*144.72)+2.3368)*-1;
            CTRL30deg_data(:,2) = ((CTRL30deg_data3(:,2)*34.634)+90)*-1;
            diffangCTRL30deg = diff(CTRL30deg_data(:,2))/0.001;
            end
            
            % Pick out the ramp and determine peak torque for each trial
            figure(1)
            clf
            hold on
            plot(CTRL30deg_data)
            plot(diffangCTRL30deg,'g')
            disp('Click around onset of FIRST ramp')
            pause
            r = ginput(2);
            tor(e,8) = max(CTRL30deg_data(r(1,1):r(2,1),1));
            vel(e,8) = mean(diffangCTRL30deg(r(1,1):r(2,1)));
            disp('Click around onset of SECOND ramp')
            pause
            q = ginput(2);
            tor(e+1,8) = max(CTRL30deg_data(q(1,1):q(2,1),1));
            vel(e+1,8) = mean(diffangCTRL30deg(q(1,1):q(2,1)));
            clear r q
            
            disp('Completed CTRL30')
 
        elseif d == 3
            filedir = [directory datafolder subjectID '\XA\' subjectID act xa];
            e = 5;
            cf = 20; % center frequency

            % CTRL30 deg/s
            CTRL30deg_file = [filedir CTRL30deg '.csv'];
            CTRL30deg_data4 = csvread([CTRL30deg_file],2,1);
            
            le=size(CTRL30deg_data4,1);
            dt=1/sampling_rate;
            T=le*dt;
            df=1/T;
            frequency=(0:le/2)*df;
            [CTRL30deg_data2,bb]=WaveletFilter(CTRL30deg_data4',dt,cf,mode,type);
            CTRL30deg_data3 = CTRL30deg_data2';
            
            if z == 1
            CTRL30deg_data(:,1) = ((CTRL30deg_data3(:,1)*144.72)+2.3368);
            CTRL30deg_data(:,2) = ((CTRL30deg_data3(:,2)*34.634)+90);
            diffangCTRL30deg = diff(CTRL30deg_data(:,2))/0.001;
            else
            CTRL30deg_data(:,1) = ((CTRL30deg_data3(:,1)*144.72)+2.3368)*-1;
            CTRL30deg_data(:,2) = ((CTRL30deg_data3(:,2)*34.634)+90)*-1;
            diffangCTRL30deg = diff(CTRL30deg_data(:,2))/0.001;
            end
            
            % Pick out the ramp and determine peak torque for each trial
            figure(1)
            clf
            hold on
            plot(CTRL30deg_data)
            plot(diffangCTRL30deg,'g')
            disp('Click around onset of FIRST ramp')
            pause
            r = ginput(2);
            tor(e,8) = max(CTRL30deg_data(r(1,1):r(2,1),1));
            vel(e,8) = mean(diffangCTRL30deg(r(1,1):r(2,1)));
            disp('Click around onset of SECOND ramp')
            pause
            q = ginput(2);
            tor(e+1,8) = max(CTRL30deg_data(q(1,1):q(2,1),1));
            vel(e+1,8) = mean(diffangCTRL30deg(q(1,1):q(2,1)));
            clear r q
            
            disp('Completed CTRL30')
            
        end
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 0 deg/s
        deg0_file = [filedir deg0 '.csv'];
        deg0_data4 = csvread([deg0_file],2,1);
        
        cf = 20; % center frequency
        le=size(deg0_data4,1);
        dt=1/sampling_rate;
        T=le*dt;
        df=1/T;
        frequency=(0:le/2)*df;
        [deg0_data2,bb]=WaveletFilter(deg0_data4',dt,cf,mode,type); % low pass
        deg0_data3 = deg0_data2';
        
        if z == 1
            deg0_data(:,1) = ((deg0_data3(:,1)*144.72)+2.3368);
            deg0_data(:,2) = ((deg0_data3(:,2)*34.634)+90);
            diffangdeg0 = diff(deg0_data(:,2))/0.001; % Superfluous - do not need
        else
            deg0_data(:,1) = ((deg0_data3(:,1)*144.72)+2.3368)*-1;
            deg0_data(:,2) = ((deg0_data3(:,2)*34.634)+90)*-1;
            diffangdeg0 = diff(deg0_data(:,2))/0.001; % Superfluous - do not need
        end
        
        % Pick out the ramp and determine peak torque for each trial
        figure(1)
        clf
        hold on
        plot(deg0_data)
        disp('Click during of FIRST MVC')
        pause
        r = ginput(2);
        tor(e,1) = max(deg0_data(r(1,1):r(2,1),1));
        vel(e,1) = 0;
        disp('Click during of SECOND MVC')
        pause
        q = ginput(2);
        tor(e+1,1) = max(deg0_data(q(1,1):q(2,1),1));
        vel(e+1,1) = 0;
        clear r q
        
        disp('Completed 0deg/s')

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 30 deg/s
        deg30_file = [filedir deg30 '.csv'];        
        deg30_data4 = csvread([deg30_file],2,1);
        
        [deg30_data2,bb]=WaveletFilter(deg30_data4',dt,cf,mode,type);
        deg30_data3 = deg30_data2';
        
        if z == 1
            deg30_data(:,1) = ((deg30_data3(:,1)*144.72)+2.3368);
            deg30_data(:,2) = ((deg30_data3(:,2)*34.634)+90);
            diffangdeg30 = diff(deg30_data(:,2))./.001;
        else
            deg30_data(:,1) = ((deg30_data3(:,1)*144.72)+2.3368)*-1;
            deg30_data(:,2) = ((deg30_data3(:,2)*34.634)+90)*-1;
            diffangdeg30 = diff(deg30_data(:,2))./.001;
        end
        
        % Pick out the ramp and determine peak torque for each trial
        figure(1)
        clf
        hold on
        plot(deg30_data)
        plot(diffangdeg30,'g')
        disp('Click around onset of FIRST ramp')
        pause
        r = ginput(2);
        tor(e,2) = max(deg30_data(r(1,1):r(2,1),1));
        vel(e,2) = mean(diffangdeg30(r(1,1):r(2,1)));
        disp('Click around onset of SECOND ramp')
        pause
        q = ginput(2);
        tor(e+1,2) = max(deg30_data(q(1,1):q(2,1),1));
        vel(e+1,2) = mean(diffangdeg30(q(1,1):q(2,1)));
        clear r q
        
        disp('Completed 30deg/s')
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 60 deg/s
        deg60_file = [filedir deg60 '.csv'];
        deg60_data4 = csvread([deg60_file],2,1);
        
        [deg60_data2,bb]=WaveletFilter(deg60_data4',dt,cf,mode,type);
        deg60_data3 = deg60_data2';
        
        if z == 1
        deg60_data(:,1) = ((deg60_data3(:,1)*144.72)+2.3368);
        deg60_data(:,2) = ((deg60_data3(:,2)*34.634)+90);
        diffangdeg60 = diff(deg60_data(:,2))/.001;
        else
        deg60_data(:,1) = ((deg60_data3(:,1)*144.72)+2.3368)*-1;
        deg60_data(:,2) = ((deg60_data3(:,2)*34.634)+90)*-1;
        diffangdeg60 = diff(deg60_data(:,2))/.001;
        end
        
        % Pick out the ramp and determine peak torque for each trial
        figure(1)
        clf
        hold on
        plot(deg60_data)
        plot(diffangdeg60,'g')
        disp('Click around onset of FIRST ramp')
        pause
        r = ginput(2);
        tor(e,3) = max(deg60_data(r(1,1):r(2,1),1));
        vel(e,3) = mean(diffangdeg60(r(1,1):r(2,1)));
        disp('Click around onset of SECOND ramp')
        pause
        q = ginput(2);
        tor(e+1,3) = max(deg60_data(q(1,1):q(2,1),1));
        vel(e+1,3) = mean(diffangdeg60(q(1,1):q(2,1)));
        clear r q
        
        disp('Completed 60deg/s')
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 90 deg/s
        deg90_file = [filedir deg90 '.csv'];
        deg90_data4 = csvread([deg90_file],2,1);
        
        [deg90_data2,bb]=WaveletFilter(deg90_data4',dt,cf,mode,type);
        deg90_data3 = deg90_data2';
        
        if z == 1
            deg90_data(:,1) = ((deg90_data3(:,1)*144.72)+2.3368);
            deg90_data(:,2) = ((deg90_data3(:,2)*34.634)+90);
            diffangdeg90 = diff(deg90_data(:,2))/.001;
        else
            deg90_data(:,1) = ((deg90_data3(:,1)*144.72)+2.3368)*-1;
            deg90_data(:,2) = ((deg90_data3(:,2)*34.634)+90)*-1;
            diffangdeg90 = diff(deg90_data(:,2))/.001;
        end
        
        % Pick out the ramp and determine peak torque for each trial
        figure(1)
        clf
        hold on
        plot(deg90_data)
        plot(diffangdeg90,'g')
        disp('Click around onset of FIRST ramp')
        pause
        r = ginput(2);
        tor(e,4) = max(deg90_data(r(1,1):r(2,1),1));
        vel(e,4) = mean(diffangdeg90(r(1,1):r(2,1)));
        disp('Click around onset of SECOND ramp')
        pause
        q = ginput(2);
        tor(e+1,4) = max(deg90_data(q(1,1):q(2,1),1));
        vel(e+1,4) = mean(diffangdeg90(q(1,1):q(2,1)));
        clear r q
        
        disp('Completed 90deg/s')
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 120 deg/s  JUST LONG ENOUGH TO DETERMINE PEAK TORQUE OVER PERIOD
        % OF VELOCITY
        deg120_file = [filedir deg120 '.csv'];
        deg120_data4 = csvread([deg120_file],2,1);
        
        [deg120_data2,bb]=WaveletFilter(deg120_data4',dt,cf,mode,type);
        deg120_data3 = deg120_data2';
        
        if z == 1
            deg120_data(:,1) = ((deg120_data3(:,1)*144.72)+2.3368);
            deg120_data(:,2) = ((deg120_data3(:,2)*34.634)+90);
            diffangdeg120 = diff(deg120_data(:,2))/.001;
        else
            deg120_data(:,1) = ((deg120_data3(:,1)*144.72)+2.3368)*-1;
            deg120_data(:,2) = ((deg120_data3(:,2)*34.634)+90)*-1;
            diffangdeg120 = diff(deg120_data(:,2))/.001;
        end
        
        % Pick out the ramp and determine peak torque for each trial
        figure(1)
        clf
        hold on
        plot(deg120_data)
        plot(diffangdeg120,'g')
        disp('Click around onset of FIRST ramp')
        pause
        r = ginput(2);
        tor(e,5) = max(deg120_data(r(1,1):r(2,1),1));
        vel(e,5) = mean(diffangdeg120(r(1,1):r(2,1)));
        disp('Click around onset of SECOND ramp')
        pause
        q = ginput(2);
        tor(e+1,5) = max(deg120_data(q(1,1):q(2,1),1));
        vel(e+1,5) = mean(diffangdeg120(q(1,1):q(2,1)));
        clear r q
        
        disp('Completed 120deg/s')
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %300 deg/s
        % NEED TO DETERMINE PEAK VEL, NOTE VEL AND INDEX, THEN FIND TORQUE
        % ALSO NEED DIFFERENT CUTOFF FREQUENCY TO FIND ACTUAL VEL
        deg300_file = [filedir deg300 '.csv'];
        deg300_data4 = csvread([deg300_file],2,1);
        
        cf = 45; %center frequency
        [deg300_data2,bb]=WaveletFilter(deg300_data4',dt,cf,mode,type);
        deg300_data3 = deg300_data2';
        
        if z == 1
            deg300_data(:,1) = ((deg300_data3(:,1)*144.72)+2.3368);
            deg300_data(:,2) = ((deg300_data3(:,2)*34.634)+90);
            diffangdeg300 = diff(deg300_data(:,2))/.001;
        else
            deg300_data(:,1) = ((deg300_data3(:,1)*144.72)+2.3368)*-1;
            deg300_data(:,2) = ((deg300_data3(:,2)*34.634)+90)*-1;
            diffangdeg300 = diff(deg300_data(:,2))/.001;
        end
        
        % Pick out the ramp and determine peak torque for each trial
        figure(1)
        clf
        hold on
        plot(deg300_data)
        plot(diffangdeg300,'g')
        disp('Click around onset of FIRST PEAK vel')
        pause
        r = ginput(2);
        [pks,locs] = findpeaks(deg300_data(r(1,1):r(2,1),1));
        temp300 = diffangdeg300(r(1,1):r(2,1),1);
        tor(e,6) = pks;
        vel(e,6) = temp300(locs);
        [pks,locs] = findpeaks(deg300_data(r(1,1):r(2,1),1),'sortstr','descend');
        temp300 = deg300_data(r(1,1):r(2,1),1);
        temp300vel=diffangdeg300(r(1,1):r(2,1),1);
        [Maxval300,indx300]=max(temp300);
        tor(e,6) = Maxval300;
        vel(e,6) = temp300vel(indx300);
        [pks,locs] = findpeaks(diffangdeg300(r(1,1):r(2,1),1),'sortstr','descend');
        temp300 = deg300_data(r(1,1):r(2,1),1);
        tor(e,6) = temp300(locs(1,1),1);
        tor(e,6) = max(deg300_data(r(1,1):r(2,1),1));
        vel(e,6) = mean(diffangdeg300(r(1,1):r(2,1)));
        disp('Click around onset of SECOND PEAK vel')
        pause
        q = ginput(2);
        temp300b = deg300_data(q(1,1):q(2,1),1);
        temp300velb=diffangdeg300(q(1,1):q(2,1),1);
        [Maxval300b,indx300b]=max(temp300b);
        tor(e+1,6) = Maxval300b;
        vel(e+1,6) = temp300velb(indx300b);
        [pks2,locs2] = findpeaks(deg300_data(q(1,1):q(2,1),1));
        temp300 = diffangdeg300(q(1,1):q(2,1),1);
        tor(e+1,6) = pks2;
        vel(e+1,6) = temp300(locs2);
        clear r q
        
        disp('Completed 300deg/s')
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 500 deg/s
        % NEED TO DETERMINE PEAK VEL, NOTE VEL AND INDEX, THEN FIND TORQUE
        % ALSO NEED DIFFERENT CUTOFF FREQUENCY TO FIND ACTUAL VEL
        deg500_file = [filedir deg500 '.csv'];
        deg500_data4 = csvread([deg500_file],2,1);
        
        cf = 70; %center frequency
        [deg500_data2,bb]=WaveletFilter(deg500_data4',dt,cf,mode,type);
        deg500_data3 = deg500_data2';
        
        if z == 1
            deg500_data(:,1) = ((deg500_data3(:,1)*144.72)+2.3368);
            deg500_data(:,2) = ((deg500_data3(:,2)*34.634)+90);
            diffangdeg500 = diff(deg500_data(:,2))/.001;
        else
            deg500_data(:,1) = ((deg500_data3(:,1)*144.72)+2.3368)*-1;
            deg500_data(:,2) = ((deg500_data3(:,2)*34.634)+90)*-1;
            diffangdeg500 = diff(deg500_data(:,2))/.001;
        end
        
        % Pick out the ramp and determine peak torque for each trial
        figure(1)
        clf
        hold on
        plot(deg500_data)
        plot(diffangdeg500,'g')
        disp('Click around onset of FIRST PEAK vel')
        pause
        r = ginput(2);
        
        temp500 = deg500_data(r(1,1):r(2,1),1);
        temp500vel=diffangdeg500(r(1,1):r(2,1),1);
        [Maxval500,indx500]=max(temp500);
        tor(e,7) = Maxval500;
        vel(e,7) = temp500vel(indx500);
        disp('Click around onset of SECOND PEAK vel')
        pause
        q = ginput(2);
        temp500b = deg500_data(q(1,1):q(2,1),1);
        temp500velb=diffangdeg500(q(1,1):q(2,1),1);
        [Maxval500b,indx500b]=max(temp500b);
        tor(e+1,7) = Maxval500b;
        vel(e+1,7) = temp500velb(indx500b);        
        
        [pks,locs] = findpeaks(deg500_data(r(1,1):r(2,1),1));
        temp500 = diffangdeg500(r(1,1):r(2,1),1);
        tor(e,7) = pks;
        vel(e,7) = temp500(locs);
        disp('Click around onset of SECOND PEAK vel')
        pause
        q = ginput(2);
        [pks2,locs2] = findpeaks(deg500_data(q(1,1):q(2,1),1));
        temp500b = diffangdeg500(q(1,1):q(2,1),1);
        tor(e+1,7) = pks2;
        vel(e+1,7) = temp500b(locs2);
        clear r q
        
        disp('Completed 500deg/s')
        
        if d == 1
%             eval(['Per4madata.active.CTRL.PX' num2str(sss) '.torque = tor;']);
%             eval(['Per4madata.active.CTRL.PX' num2str(sss) '.velocity = vel;']);
%             save('Per4madata.mat','Per4madata')
            disp('Completed CTRL')
            clear tor vel
        elseif d == 2
%             eval(['Per4madata.active.UA.PX' num2str(sss) '.torque = tor;']);
%             eval(['Per4madata.active.UA.PX' num2str(sss) '.velocity = vel;']);
%             save('Per4madata.mat','Per4madata')
            disp('Completed UA')
            clear tor vel
        else
%             eval(['Per4madata.active.XA.PX' num2str(sss) '.torque = tor;']);
%             eval(['Per4madata.active.XA.PX' num2str(sss) '.velocity = vel;']);
%             save('Per4madata.mat','Per4madata')
            disp('Completed XA')
            clear tor vel
        end
        
        
    clear d
        
    end

end
        




% %%   ARM WINGATE   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% load('Per4madata.mat')
% 
% % Active labels
% WG = '_WG';
% ctrl = 'CTRL';
% ua = 'UA';
% xa = 'XA';
% ch=2; % #channels: 1=armerg, 2=accel
% thresh = 2;
% 
% 
% for ssss = 6
%     
%     subj = sprintf('PX0%d',(ssss));
%     
%     for h = 1:3
%         if h == 1
%             file = [directory datafolder subj '\' ctrl '\' subj '_' ctrl WG '.emg'];
%         elseif h == 2
%             file = [directory datafolder subj '\' ua '\' subj '_' ua WG '.emg'];
%         elseif h ==3
%             file = [directory datafolder subj '\XACT\' subj '_' xa WG '.emg'];
%         end
%         
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % Read in .emg file
%         A = fopen(file, 'r'); %r stands for open file to read
%         B = fread(A,'float'); %float is float point numbers 32 bits
%         
%         [totalcal, to] = size(B);
%         m = totalcal/(ch+1);
%         data = reshape(B, ch+1, m); % all one file, so must split 
%         data = data';
%         ST = fclose(A);
%         armergtest(:,1) = data(:,2); %CHANGE
%         armergtest(:,2) = data(:,1); %CHANGE
%         
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % find onset trigger
%         figure(1)
%         clf
%         hold on
%         plot(armergtest)
%         prompt = 'Does it have a START trigger?\n Press 0 for a start trigger, 1 for only a stop trigger, and 2 for no trigger at all';
%         c = input(prompt) %Type in how many steps to process
%         
%         if c == 0 % Meaning it does have a start trigger
%             
%             onset = find(armergtest>5,1);
%             fullwingate = armergtest(onset:(onset+71999),2);
%             
%             % Determine rate of crank
%             % First find when voltage drops, only onset, record index
%             k = 1;
%             for i = 1:71999
%                 difference = fullwingate(i,1)-fullwingate(i+1,1);
%                 if difference > thresh
%                     index(k,1) = i;
%                     k = k+1;
%                 end
%             end
%             
%             % Find frame# diff between indices
%             for l = 1:length(index)-1;
%                 timebtwn(l,1) = index(l+1)-index(l);
%             end
%             
%             % Convert frame# to time in min (2400 f/s * 60 s/min)
%             wgrate = 144000./timebtwn;
%             
% %             % Then do this for the first and last 5 sec of wingate
% %             first5 = fullwingate(1:12000);
% %             last5 = fullwingate(60001:72000);
%             
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         elseif c == 1 %if it has no start trigger but has an end trigger
%             offset = find(armergtest>5,1);
%             fullwingate = armergtest((offset-72000):offset,2);
%             
%             % Determine rate of crank
%             % First find when voltage drops, only onset, record index
%             k = 1;
%             for i = 1:71999
%                 difference = fullwingate(i,1)-fullwingate(i+1,1);
%                 if difference > thresh
%                     index(k,1) = i;
%                     k = k+1;
%                 end
%             end
%             
%             % Find frame# diff between indices
%             for l = 1:length(index)-1;
%                 timebtwn(l,1) = index(l+1)-index(l);
%             end
%             
%             % Convert frame# to time in min (2400 f/s * 60 s/min)
%             wgrate = 144000./timebtwn;
%             
% %             % Then do this for the first and last 5 sec of wingate
% %             first5 = fullwingate(1:12000);
% %             last5 = fullwingate(60001:72000);
%             
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         else %you're in fuck city, bail yourself out
%             disp('Click where you think the wingate begins, slightly before they ramp it up')
%             k = ginput;
%             fullwingate = armergtest(k:(k+71999),2);
%             
%             % Determine rate of crank
%             % First find when voltage drops, only onset, record index
%             k = 1;
%             for i = 1:71999
%                 difference = fullwingate(i,1)-fullwingate(i+1,1);
%                 if difference > thresh
%                     index(k,1) = i;
%                     k = k+1;
%                 end
%             end
%             
%             % Find frame# diff between indices
%             for l = 1:length(index)-1;
%                 timebtwn(l,1) = index(l+1)-index(l);
%             end
%             
%             % Convert frame# to time in min (2400 f/s * 60 s/min)
%             wgrate = 144000./timebtwn;
%             
% %             % Then do this for the first and last 5 sec of wingate
% %             first5 = fullwingate(1:12000);
% %             last5 = fullwingate(60001:72000);
%             
%         end
%         
%         if h == 1
%             eval(['Per4madata.wingate.CTRL.PX0' num2str(ssss) '.wingate = fullwingate;']);
%             eval(['Per4madata.wingate.CTRL.PX0' num2str(ssss) '.rpm = wgrate;']);
%             save('Per4madata.mat','Per4madata')
%             disp('Completed CTRL Wingate')
%             clear tor vel
%         elseif h == 2
%             eval(['Per4madata.wingate.UA.PX0' num2str(ssss) '.wingate = fullwingate;']);
%             eval(['Per4madata.wingate.UA.PX0' num2str(ssss) '.rpm = wgrate;']);
%             save('Per4madata.mat','Per4madata')
%             disp('Completed UA Wingate')
%             clear tor vel
%         else
%             eval(['Per4madata.wingate.XA.PX0' num2str(ssss) '.wingate = fullwingate;']);
%             eval(['Per4madata.wingate.XA.PX0' num2str(ssss) '.rpm = wgrate;']);
%             save('Per4madata.mat','Per4madata')
%             disp('Completed XA Wingate')
%             clear tor vel
%         end
%         
%         clear armergtest data fullwingate wgrate l k i c timebtwn
%     end
%     
% end

% save('Per4madata.mat','Per4madata')








