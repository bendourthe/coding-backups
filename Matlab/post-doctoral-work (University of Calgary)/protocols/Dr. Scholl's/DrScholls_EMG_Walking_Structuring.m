%% DESCRIPTION:
%
% This code was designed to enable the formatting of the walking EMG data 
% obtained via the Dr. Scholl's Research Project:
%
%   Project Title: Quantifying the effects of different insole
%   configurations on fatigue
%
%   Author: Benjamin Dourthe
%
%   Last update: October 14th, 2018
%
%% Input: 
%   - Processed walking EMG data (Walk1.csv, Walk2.csv and Walk3.csv)
%
%% Output:
%   - Walk.csv: Structured EMG data from the intervention
%       The rows of each file are as followed:
%           1. Absolute change from Baseline to Post-intervention
%           2. Relative change from Baseline to Post-intervention
%       The columns of each file are as followed:
%           Each series of five columns represent the results for the
%           different conditions in this order: CTRL, CF, BPR, MG, REP
%           (zeros are added when the corresponding condition wasn't
%           tested)
%
%% Dependencies:
%
%   Files:
%       - Processed walking EMG data (Walk1.csv, Walk2.csv and Walk3.csv)
%
%   Functions:
%       - None

%% Initialization

    clear ; close all; clc
    
%% Subject ID   

    subjectID = '028';
    % Define third session: type 3 for bpr, 4 for mg, 5 for repeated
    condition3 = 4;
    
%% Directory

    pathName = ['F:\Dr. Scholls\Phase 3\' subjectID '\Results'];
    cd(pathName);
    
%% Execution  

    % Sessions
    sessions = num2str([1:3].','%01d');
    
    % Filenames
    walk = [strcat('\',subjectID,'_Walk',sessions,'.csv')];
    
    % Load data
    data1 = csvread([pathName walk(1,:)]);
    data2 = csvread([pathName walk(2,:)]);
    data3 = csvread([pathName walk(3,:)]);
    
    % Calculate relative changes
    data1r = (data1(end,:)-data1(3,:))*100./data1(3,:);
    data2r = (data2(end,:)-data2(3,:))*100./data2(3,:);
    data3r = (data3(end,:)-data3(3,:))*100./data3(3,:);
    
    % Created output matrix with desired structure
    if condition3 == 3
        dataEMGstruc = ...
            [data1(1,5) data2(1,5) data3(1,5) 0 0 data1(1,6) data2(1,6) data3(1,6) 0 0 ...
             data1(1,7) data2(1,7) data3(1,7) 0 0 data1(1,8) data2(1,8) data3(1,8) 0 0 ...
             data1(1,1) data2(1,1) data3(1,1) 0 0 data1(1,2) data2(1,2) data3(1,2) 0 0 ...
             data1(1,3) data2(1,3) data3(1,3) 0 0 data1(1,4) data2(1,4) data3(1,4) 0 0; ...
             data1r(1,5) data2r(1,5) data3r(1,5) 0 0 data1r(1,6) data2r(1,6) data3r(1,6) 0 0 ...
             data1r(1,7) data2r(1,7) data3r(1,7) 0 0 data1r(1,8) data2r(1,8) data3r(1,8) 0 0 ...
             data1r(1,1) data2r(1,1) data3r(1,1) 0 0 data1r(1,2) data2r(1,2) data3r(1,2) 0 0 ...
             data1r(1,3) data2r(1,3) data3r(1,3) 0 0 data1r(1,4) data2r(1,4) data3r(1,4) 0 0];
    elseif condition3 == 4
        dataEMGstruc = ...
            [data1(1,5) data2(1,5) 0 data3(1,5) 0 data1(1,6) data2(1,6) 0 data3(1,6) 0 ...
             data1(1,7) data2(1,7) 0 data3(1,7) 0 data1(1,8) data2(1,8) 0 data3(1,8) 0 ...
             data1(1,1) data2(1,1) 0 data3(1,1) 0 data1(1,2) data2(1,2) 0 data3(1,2) 0 ...
             data1(1,3) data2(1,3) 0 data3(1,3) 0 data1(1,4) data2(1,4) 0 data3(1,4) 0; ...
             data1r(1,5) data2r(1,5) 0 data3r(1,5) 0 data1r(1,6) data2r(1,6) 0 data3r(1,6) 0 ...
             data1r(1,7) data2r(1,7) 0 data3r(1,7) 0 data1r(1,8) data2r(1,8) 0 data3r(1,8) 0 ...
             data1r(1,1) data2r(1,1) 0 data3r(1,1) 0 data1r(1,2) data2r(1,2) 0 data3r(1,2) 0 ...
             data1r(1,3) data2r(1,3) 0 data3r(1,3) 0 data1r(1,4) data2r(1,4) 0 data3r(1,4) 0];
    elseif condition3 == 5
        dataEMGstruc = ...
            [data1(1,5) data2(1,5) 0 0 data3(1,5) data1(1,6) data2(1,6) 0 0 data3(1,6) ...
             data1(1,7) data2(1,7) 0 0 data3(1,7) data1(1,8) data2(1,8) 0 0 data3(1,8) ...
             data1(1,1) data2(1,1) 0 0 data3(1,1) data1(1,2) data2(1,2) 0 0 data3(1,2) ...
             data1(1,3) data2(1,3) 0 0 data3(1,3) data1(1,4) data2(1,4) 0 0 data3(1,4); ...
             data1r(1,5) data2r(1,5) 0 0 data3r(1,5) data1r(1,6) data2r(1,6) 0 0 data3r(1,6) ...
             data1r(1,7) data2r(1,7) 0 0 data3r(1,7) data1r(1,8) data2r(1,8) 0 0 data3r(1,8) ...
             data1r(1,1) data2r(1,1) 0 0 data3r(1,1) data1r(1,2) data2r(1,2) 0 0 data3r(1,2) ...
             data1r(1,3) data2r(1,3) 0 0 data3r(1,3) data1r(1,4) data2r(1,4) 0 0 data3r(1,4)];
    end
    
    % Export results
    xlswrite([subjectID '_Walk.xlsx'],dataEMGstruc);
    
    disp('     *****     DONE WITH STRUCTURING WALKING EMG DATA     *****     ')