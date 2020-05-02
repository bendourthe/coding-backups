clear all; close all;
%% DESCRIPTION:
%
% This code was designed to enable the processing of the ARM WINGATE data
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
%   - Raw arm ergometer data (1 CTRL, 1 UA, 1 XA)
%
%% Output:
%   - Excel file with all the Wingate parameters 
%
%% Dependencies:
%   Functions:
%       - None
  

%% ID and DIRECTORY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Subject ID
    ID = 17;     % Type subject ID number (2 digits after PX)
    if ID < 10
        subjectID = sprintf('PX0%d',(ID));
    else
        subjectID = sprintf('PX%d',(ID));
    end

    % Current directory (where the data are stored)
    cd 'C:\Users\bdour\OneDrive\Work\Calgary\Projects\Per4ma - X-ACT\Data Collection\Data\SHIRT\Valid\'
    directory = 'C:\Users\bdour\OneDrive\Work\Calgary\Projects\Per4ma - X-ACT\Data Collection\Data\SHIRT\Valid\';

    % Data folder
    datafolder = 'Wingate';

    % Labels
    WG = '_WG';
    ctrl = '_CTRL';
    ua = '_UA';
    xa = '_XA';    


%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Manual or Auto mode
    mode = 2;   % Type 1 for AUTO mode, Type 2 for MANUAL mode

    % Sampling rate (3 values in case the sampling rate was different
    % between trials)
    rate1 = 2400;           % in Hz
    rate2 = 2400;           % in Hz
    rate3 = 2400;           % in Hz
    
    % Wingate time (3 values in case the recording time was not enough to
    % include the full Wingate)
    wg_time1 = 30;          % in sec
    wg_time2 = 30;          % in sec
    wg_time3 = 30;          % in sec

    % Number of channels
    ch1 = 2;
    ch2 = 2;
    ch3 = 2;
    
    % Threshold used to determine location of each rotation
    thresh = 2;
    
    % Ergometer
    res = 4;            % Resistance of the ergometer
    d = 2*pi*0.14;      % Distance per rotation (= perimeter of the wheel)
    

%% CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 6
    
    path = [directory subjectID '\' datafolder];
    
    for h = 1:3
        if h == 1
            file = [subjectID WG ctrl '.emg'];
            rate = rate1;
            ch = ch1;
            wg_time = wg_time1;
            data_time = wg_time*rate-1; % How many data points for 30secs
        elseif h == 2
            file = [subjectID WG ua '.emg'];
            rate = rate2;
            ch = ch2;
            wg_time = wg_time2;
            data_time = wg_time*rate-1; % How many data points for 30secs
        elseif h ==3
            file = [subjectID WG xa '.emg'];
            rate = rate3;
            ch = ch3;
            wg_time = wg_time3;
            data_time = wg_time*rate-1; % How many data points for 30secs
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Read in .emg file
        data = read_data(path,file,ch);
        
        pulse = data(:,2);
        ergo = data(:,1);
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Defines the time scale according to the frequency rate
        time = 0+1/rate:1/rate:length(data)/rate;

        if mode == 1    % AUTO mode
            
            % Find onset trigger
            figure
            clf
            hold on
            plot(time,pulse,'r')
            plot(time,ergo,'b')
            disp('Press 0 if you see NO trigger')
            disp('Press 1 if you only see a START trigger (or both START and STOP triggers)')
            disp('Press 2 if you only see a STOP trigger')
            prompt = 'Type answer and PRESS ENTER: ';
            c = input(prompt)

            if c == 0       % NO trigger
                disp('Click where you think the wingate begins')
                disp('i.e. slightly before the rate seem to increase')
                onset = ginput(1)*rate;
                timePoints = [onset(1)/rate;(onset(1)+data_time)/rate];
                fullwingate = ergo(onset(1):onset(1)+data_time);            

            elseif c == 1   % START trigger only
                onset = find(pulse>2);
                timePoints = [onset(1)/rate;(onset(1)+data_time)/rate];
                fullwingate = ergo(onset(1):onset(1)+data_time);

            elseif c == 2    % END trigger only
                [M,onset] = max(pulse);
                timePoints = [(onset-data_time)/rate;onset/rate];
                fullwingate = ergo(onset-data_time:onset);
            end
            
        elseif mode == 2   % MANUAL mode
            
            % Find onset trigger
            figure
            clf
            hold on
            plot(time,pulse,'r')
            plot(time,ergo,'b')
            disp('Click where you think the wingate begins, then when you think it ends')
            onset = ginput(2)*rate;
            timePoints = [onset(1,1)/rate;onset(2,1)/rate];
            fullwingate = ergo(onset(1,1):onset(2,1));
            data_time = onset(2,1)-onset(1,1);
            wg_time = data_time/rate;
            
        end
        
        % Plots the corresponding Wingate for validation
        figure
        clf
        hold on
        plot(time,pulse,'r')
        plot(time,ergo,'b')
        plot(timePoints(1),max(ergo),'g*','MarkerSize',10)
        plot(timePoints(2),max(ergo),'g*','MarkerSize',10)
        disp('PRESS ENTER if ramps are correctly located - PRESS Ctr+C in the command window if not')
        pause
            
            % Determine location of each rotation during Wingate
            k = 1;
            for i = 1:data_time
                difference = fullwingate(i)-fullwingate(i+1);
                if difference > thresh
                    locR(k) = i;
                    k = k+1;
                end
            end
            
            % Find how long each rotation is (in secs)
            for j = 1:length(locR)-1
                rps(j,1) = (locR(j+1)-locR(j))/rate;
            end
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % RPM

                % Create the accumulative vector the see when each rotation
                % happens (first would be at 0, last at 30 for a 30 secs
                % Wingate)
                for l = 1:length(rps)
                    rps_sum(l)=sum(rps(1:l));
                end
                
                % RPM for every 5secs block
                for m = 1:wg_time/5
                    
                    if m == 1
                        idx = find(rps_sum < 5);                        % Find the location of the first 5 secs landmark
                        rpm(m) = 60/mean(rps(1:idx(end)));              % Mean RPM for the first 5secs of the Wingate
                    else
                        idx = find(rps_sum < 5*m & rps_sum > 5*(m-1));  % Find the location of the each 5 secs block
                        rpm(m) = 60/mean(rps(idx(1):idx(end)));         % Mean RPM for every 5secs of the Wingate
                    end
                end
                
                % Mean RPM for the full Wingate (30secs)
                mrpm(h) = 60/mean(rps); 
                
                % Peak RPM (first 5secs of the Wingate)
                prpm(h) = rpm(1);
                
                % Min RPM (last 5secs of the Wingate)
                lrpm(h) = rpm(end); 
                                
            % Power Output
            mpo(h) = (res * mrpm(h) * d)/5;   % Mean Power Output (over the full Wingate)
            ppo(h) = (res * rpm(1) * d)/5;          % Peak Power Output (first 5secs)
            lpo(h) = (res * rpm(end) * d)/5;        % Lowest Power Output (last 5secs)
            
            % Total work during Wingate
            tw(h) = sum((res * rpm * d)/5);
            
            % Fatigue Index
            fi(h) = (ppo(h)-lpo(h))/ppo(h)*100;
    
            clear data pulse ergo time fullwingate locR rps rps_sum
    end
    
end
            
          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = [mpo' ppo' lpo' mrpm' prpm' lrpm' fi' tw'];
Results = [M(1,:) M(2,:) M(3,:)];
csvwrite(['Results_' subjectID WG '.csv'], Results);
    
close all           
        