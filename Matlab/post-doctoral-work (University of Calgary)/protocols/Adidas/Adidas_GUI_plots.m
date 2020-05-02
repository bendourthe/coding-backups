%% DESCRIPTION:
%
% This code was designed to enable the visualization and analysis of the
% data provided via the Athlete Assessment project using a GUI
%
%   Project Title: 
%
%   Author: Benjamin Dourthe
%
%   Last update: May 29th, 2018
%
%% Input: 
%   - Excel Master Sheet (.xlsx)
%
%% Output:
%   - 
%
%% Dependencies:
%
%   Files:
%       - 2017 - 12 Athlete Assessment - all data - Ben
%
%   Functions:
%       - uigetvariables.m


%% Initialization

    clear ; close all; clc
    
%% Directory

    cd 'C:\Users\bdour\OneDrive\Work\Calgary\Projects\Adidas\';
    
%% Filename

    fileName = '2017 - 12 Athlete Assessment - all data - Ben.xlsx';
    
%% Number of subjects

    nSub = 72;

%% Execution

    %% Import data
    data.HP = xlsread(fileName,1);
    data.Rise = xlsread(fileName,2);
    data.Ant = xlsread(fileName,3);
    data.FSens = xlsread(fileName,4);
    data.YBal = xlsread(fileName,5);
    data.GStr = xlsread(fileName,6);
    data.Cur = xlsread(fileName,7);
    data.FScan = xlsread(fileName,8);
    data.Kin = xlsread(fileName,9);
    
    %% Grouping
        % Identify what subject belongs to what group
        
    
    %% Generate relevant variable
        % Line 1: Generate variable
        % Line 2: Find missing and/or NaN data and save the subect ID in
        % NaN structure
        % Line 3: Remove NaN values from variable
        
        % Note: to facilitate the manual selection of variables to plot,
        % give variable an intuitive name
        
        % History & Preference
    
            % Hours on feet per day
                % Group 1
                HoursOnFeetPerDay_group1 = data.HP(group1,3);
                NaN.HoursOnFeetPerDay_group1 = find(isnan(HoursOnFeetPerDay_group1)==1);
                HoursOnFeetPerDay_group1(NaN.HoursOnFeetPerDay_group1) = [];
                % Group 2
                HoursOnFeetPerDay_group2 = data.HP(group2,3);
                NaN.HoursOnFeetPerDay_group2 = find(isnan(HoursOnFeetPerDay_group2)==1);
                HoursOnFeetPerDay_group2(NaN.HoursOnFeetPerDay_group2) = [];
                % Group 3
                HoursOnFeetPerDay_group3 = data.HP(group3,3);
                NaN.HoursOnFeetPerDay_group3 = find(isnan(HoursOnFeetPerDay_group3)==1);
                HoursOnFeetPerDay_group3(NaN.HoursOnFeetPerDay_group3) = [];
                % Group 4
                HoursOnFeetPerDay_group4 = data.HP(group4,3);
                NaN.HoursOnFeetPerDay_group4 = find(isnan(HoursOnFeetPerDay_group4)==1);
                HoursOnFeetPerDay_group4(NaN.HoursOnFeetPerDay_group4) = [];
            
        % Rise
        
        % Anthropometrics
        
        % Foot Sensitivity
        
        % Y-Balance
        
        % Grip Strength
        
        % Currex
        
        % Foot Scan
        
        % Kinect
            
    %% Select variables to plot
    plotVar = uigetvariables({'Variable 1','Variable 2'}, ...
        'Introduction',['Select variables you wish to plot']);
    
    %% Plot selected data
    
        % Generate a matrix with the mean and SD of each selected variable
        A = [mean(plotVar{1,1}) mean(plotVar{1,2})];
        B = [std(plotVar{1,1}) std(plotVar{1,2})];
    
        % Generate barplots using selected variables
        fig = figure;
        hold on
        b = bar(1:2,A);
        e = errorbar(1:2,A,B,'.');
        
        % Figure properties
            % Figure background color
            fig.Color = 'w';
            % Axes properties
            ax = gca;
            ax.FontSize = 20;
            ax.FontWeight = 'bold';
            ax.LineWidth = 2;
            ax.XTick = [];
            ax.XTickLabel = [];
            ax.TickDir = 'out'
            % Barplot properties
            b.FaceColor = 'w';
            b.LineWidth = 2;
            b.BarWidth = 0.5;
            % Errorbar properties
            e.Color = 'b';
            e.LineWidth = 2;


    
    
    
