function [emg_1, emg_2, emg_3, emg_4, emg_5] = import_emg(pathName,subjectID,...
    numSessions,sessions,recTime,chanNum)
%% DESCRIPTION:
%
% [emg_1, emg_2, emg_3, emg_4, emg_5] = import_emg(pathName,subjectID,
%                                       numSessions,sessions,recTime,chanNum)
%
% Imports all EMG signals recorded during different walking trials for with
% different conditions.
%
%      Note: with current code, only works with recordings of 5 mins each,
%      up to 60 mins recorded per trial, up to 5 different conditions
%      studied (originally written for the Dr. Scholl's project)
%
%   Author: Benjamin Dourthe
%
%   Last update: June 7th, 2018
%
%% Input: 
%   - pathName: path towards the folder where the data is stored
%   - subjectID: ID of the corresponding subject (e.g. '10')
%   - numSessions: number of session/condition studied
%   - sessions: character (char) array corresponding to the number of
%   session/condition studied (e.g. for 3 sessions/conditions, sessions =
%   3x1 char array with '1', '2', and '3' in each row)
%   - recTime: total recording time for each trial (in mins)
%   - chanNum: number of channels used during each recording
%
%% Output:
%   - emg_n: reconstructed EMG data for the nth session/condition
%            size = n x number of channels + 1 (for time)
%            with n = recording time (mins) x 60 (to secs) x 2400 (sampling
%            frequency)
%
%% Dependencies:
%   Functions:
%       - None
%   Files:
%       - None


%% FUNCTION  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filenames

    emg1 = [strcat('\',sessions,'_',subjectID,'_emg_1.emg')];
    emg2 = [strcat('\',sessions,'_',subjectID,'_emg_2.emg')];
    emg3 = [strcat('\',sessions,'_',subjectID,'_emg_3.emg')];
    emg4 = [strcat('\',sessions,'_',subjectID,'_emg_4.emg')];
    emg5 = [strcat('\',sessions,'_',subjectID,'_emg_5.emg')];
    emg6 = [strcat('\',sessions,'_',subjectID,'_emg_6.emg')];
    emg7 = [strcat('\',sessions,'_',subjectID,'_emg_7.emg')];
    emg8 = [strcat('\',sessions,'_',subjectID,'_emg_8.emg')];
    emg9 = [strcat('\',sessions,'_',subjectID,'_emg_9.emg')];
    emg10 = [strcat('\',sessions,'_',subjectID,'_emg_10.emg')];
    emg11 = [strcat('\',sessions,'_',subjectID,'_emg_11.emg')];
    emg12 = [strcat('\',sessions,'_',subjectID,'_emg_12.emg')];

% Data reconstruction (assembles each recording to form one large signal
% for each session/condition)

    for i=1:numSessions
        if recTime >= 5
            data1 = read_data(pathName,emg1(i,:),chanNum);
        end
        if recTime >= 10
            data2 = read_data(pathName,emg2(i,:),chanNum);
        end
        if recTime >= 15
            data3 = read_data(pathName,emg3(i,:),chanNum);
        end
        if recTime >= 20
            data4 = read_data(pathName,emg4(i,:),chanNum);
        end
        if recTime >= 25
            data5 = read_data(pathName,emg5(i,:),chanNum);
        end
        if recTime >= 30
            data6 = read_data(pathName,emg6(i,:),chanNum);
        end
        if recTime >= 35
            data7 = read_data(pathName,emg7(i,:),chanNum);
        end
        if recTime >= 40
            data8 = read_data(pathName,emg8(i,:),chanNum);
        end
        if recTime >= 45
            data9 = read_data(pathName,emg9(i,:),chanNum);
        end
        if recTime >= 50
            data10 = read_data(pathName,emg10(i,:),chanNum);
        end
        if recTime >= 55
            data11 = read_data(pathName,emg11(i,:),chanNum);
        end
        if recTime >= 60
            data12 = read_data(pathName,emg12(i,:),chanNum);
        end
        
        % Condition 1
        if i==1
            if recTime == 5
                emg_1 = [data1(2:end,:)];
            elseif recTime == 10
                emg_1 = [data1(2:end,:);data2(2:end,:)];
            elseif recTime == 15
                emg_1 = [data1(2:end,:);data2(2:end,:);data3(2:end,:)];
            elseif recTime == 20
                emg_1 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:)];
            elseif recTime == 25
                emg_1 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:);data5(2:end,:)];
            elseif recTime == 30
                emg_1 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:);data5(2:end,:);...
                data6(2:end,:)];
            elseif recTime == 35
                emg_1 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:);data5(2:end,:);...
                data6(2:end,:);data7(2:end,:)];
            elseif recTime == 40
                emg_1 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:);data5(2:end,:);...
                data6(2:end,:);data7(2:end,:);data8(2:end,:)];
            elseif recTime == 45
                emg_1 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:);data5(2:end,:);...
                data6(2:end,:);data7(2:end,:);data8(2:end,:);data9(2:end,:)];
            elseif recTime == 50
                emg_1 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:);data5(2:end,:);...
                data6(2:end,:);data7(2:end,:);data8(2:end,:);data9(2:end,:);data10(2:end,:)];
            elseif recTime == 55
                emg_1 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:);data5(2:end,:);...
                data6(2:end,:);data7(2:end,:);data8(2:end,:);data9(2:end,:);data10(2:end,:);...
                data11(2:end,:)];
            elseif recTime == 60
                emg_1 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:);data5(2:end,:);...
                data6(2:end,:);data7(2:end,:);data8(2:end,:);data9(2:end,:);data10(2:end,:);...
                data11(2:end,:);data12(2:end,:)];
            end
        
        % Condition 2
        elseif i==2
            if recTime == 5
                emg_2 = [data1(2:end,:)];
            elseif recTime == 10
                emg_2 = [data1(2:end,:);data2(2:end,:)];
            elseif recTime == 15
                emg_2 = [data1(2:end,:);data2(2:end,:);data3(2:end,:)];
            elseif recTime == 20
                emg_2 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:)];
            elseif recTime == 25
                emg_2 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:);data5(2:end,:)];
            elseif recTime == 30
                emg_2 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:);data5(2:end,:);...
                data6(2:end,:)];
            elseif recTime == 35
                emg_2 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:);data5(2:end,:);...
                data6(2:end,:);data7(2:end,:)];
            elseif recTime == 40
                emg_2 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:);data5(2:end,:);...
                data6(2:end,:);data7(2:end,:);data8(2:end,:)];
            elseif recTime == 45
                emg_2 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:);data5(2:end,:);...
                data6(2:end,:);data7(2:end,:);data8(2:end,:);data9(2:end,:)];
            elseif recTime == 50
                emg_2 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:);data5(2:end,:);...
                data6(2:end,:);data7(2:end,:);data8(2:end,:);data9(2:end,:);data10(2:end,:)];
            elseif recTime == 55
                emg_2 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:);data5(2:end,:);...
                data6(2:end,:);data7(2:end,:);data8(2:end,:);data9(2:end,:);data10(2:end,:);...
                data11(2:end,:)];
            elseif recTime == 60
                emg_2 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:);data5(2:end,:);...
                data6(2:end,:);data7(2:end,:);data8(2:end,:);data9(2:end,:);data10(2:end,:);...
                data11(2:end,:),data12(2:end,:)];
            end            
        
        % Condition 3
        elseif i==3
            if recTime == 5
                emg_3 = [data1(2:end,:)];
            elseif recTime == 10
                emg_3 = [data1(2:end,:);data2(2:end,:)];
            elseif recTime == 15
                emg_3 = [data1(2:end,:);data2(2:end,:);data3(2:end,:)];
            elseif recTime == 20
                emg_3 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:)];
            elseif recTime == 25
                emg_3 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:);data5(2:end,:)];
            elseif recTime == 30
                emg_3 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:);data5(2:end,:);...
                data6(2:end,:)];
            elseif recTime == 35
                emg_3 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:);data5(2:end,:);...
                data6(2:end,:);data7(2:end,:)];
            elseif recTime == 40
                emg_3 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:);data5(2:end,:);...
                data6(2:end,:);data7(2:end,:);data8(2:end,:)];
            elseif recTime == 45
                emg_3 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:);data5(2:end,:);...
                data6(2:end,:);data7(2:end,:);data8(2:end,:);data9(2:end,:)];
            elseif recTime == 50
                emg_3 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:);data5(2:end,:);...
                data6(2:end,:);data7(2:end,:);data8(2:end,:);data9(2:end,:);data10(2:end,:)];
            elseif recTime == 55
                emg_3 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:);data5(2:end,:);...
                data6(2:end,:);data7(2:end,:);data8(2:end,:);data9(2:end,:);data10(2:end,:);...
                data11(2:end,:)];
            elseif recTime == 60
                emg_3 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:);data5(2:end,:);...
                data6(2:end,:);data7(2:end,:);data8(2:end,:);data9(2:end,:);data10(2:end,:);...
                data11(2:end,:),data12(2:end,:)];
            end
        
        % Condition 4
        elseif i==4
            if recTime == 5
                emg_4 = [data1(2:end,:)];
            elseif recTime == 10
                emg_4 = [data1(2:end,:);data2(2:end,:)];
            elseif recTime == 15
                emg_4 = [data1(2:end,:);data2(2:end,:);data3(2:end,:)];
            elseif recTime == 20
                emg_4 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:)];
            elseif recTime == 25
                emg_4 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:);data5(2:end,:)];
            elseif recTime == 30
                emg_4 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:);data5(2:end,:);...
                data6(2:end,:)];
            elseif recTime == 35
                emg_4 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:);data5(2:end,:);...
                data6(2:end,:);data7(2:end,:)];
            elseif recTime == 40
                emg_4 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:);data5(2:end,:);...
                data6(2:end,:);data7(2:end,:);data8(2:end,:)];
            elseif recTime == 45
                emg_4 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:);data5(2:end,:);...
                data6(2:end,:);data7(2:end,:);data8(2:end,:);data9(2:end,:)];
            elseif recTime == 50
                emg_4 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:);data5(2:end,:);...
                data6(2:end,:);data7(2:end,:);data8(2:end,:);data9(2:end,:);data10(2:end,:)];
            elseif recTime == 55
                emg_4 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:);data5(2:end,:);...
                data6(2:end,:);data7(2:end,:);data8(2:end,:);data9(2:end,:);data10(2:end,:);...
                data11(2:end,:)];
            elseif recTime == 60
                emg_4 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:);data5(2:end,:);...
                data6(2:end,:);data7(2:end,:);data8(2:end,:);data9(2:end,:);data10(2:end,:);...
                data11(2:end,:),data12(2:end,:)];
            end
            
        % Condition 5
        elseif i==5
            if recTime == 5
                emg_5 = [data1(2:end,:)];
            elseif recTime == 10
                emg_5 = [data1(2:end,:);data2(2:end,:)];
            elseif recTime == 15
                emg_5 = [data1(2:end,:);data2(2:end,:);data3(2:end,:)];
            elseif recTime == 20
                emg_5 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:)];
            elseif recTime == 25
                emg_5 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:);data5(2:end,:)];
            elseif recTime == 30
                emg_5 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:);data5(2:end,:);...
                data6(2:end,:)];
            elseif recTime == 35
                emg_5 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:);data5(2:end,:);...
                data6(2:end,:);data7(2:end,:)];
            elseif recTime == 40
                emg_5 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:);data5(2:end,:);...
                data6(2:end,:);data7(2:end,:);data8(2:end,:)];
            elseif recTime == 45
                emg_5 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:);data5(2:end,:);...
                data6(2:end,:);data7(2:end,:);data8(2:end,:);data9(2:end,:)];
            elseif recTime == 50
                emg_5 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:);data5(2:end,:);...
                data6(2:end,:);data7(2:end,:);data8(2:end,:);data9(2:end,:);data10(2:end,:)];
            elseif recTime == 55
                emg_5 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:);data5(2:end,:);...
                data6(2:end,:);data7(2:end,:);data8(2:end,:);data9(2:end,:);data10(2:end,:);...
                data11(2:end,:)];
            elseif recTime == 60
                emg_5 = [data1(2:end,:);data2(2:end,:);data3(2:end,:);data4(2:end,:);data5(2:end,:);...
                data6(2:end,:);data7(2:end,:);data8(2:end,:);data9(2:end,:);data10(2:end,:);...
                data11(2:end,:),data12(2:end,:)];
            end
        end   
        
        clear data1 data2 data3 data4 data5 data6 data7 data8 data9 data10 data11 data12
    
    end