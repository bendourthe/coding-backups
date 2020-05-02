clear all; close all;
%% DESCRIPTION:
%
% Imports the CSV document including the raw data from all the
% questionnaires of the Dr. Scholl's study and isolate the data from a
% specified subject to fit the corresponding Excel template
%
%   Author: Benjamin Dourthe
%
%   Last update: Dec. 11th, 2017
%
%% Input: 
%   - csvData: csv document including the raw data of every questionnaire
%   of the Dr. Scholl's study
%
%% Output:
%   - CSV document with the data of a single subject

%% SCRIPT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Subject ID
    subjectID = 101;
    
% Current directory (where the data are stored for one participant)
cd 'C:\Users\bdour\OneDrive\Work\Calgary\Projects\Dr. Scholls\Data\';
pathName = 'C:\Users\bdour\OneDrive\Work\Calgary\Projects\Dr. Scholls\Data\';

% Load data
    data = xlsread(strcat('Questionnaires (Raw Data).xlsx'));
    