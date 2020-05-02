%% DESCRIPTION:
%
% This code was designed to enable the automated classification of subjects
% based on the subjective ranking of 9 different shoe conditions
%
%   Project Title: Adidas Athlete Assessment Running
%
%   Author: Benjamin Dourthe
%
%   Last update: May 3rd, 2018
%
%% Input: 
%   - Ranking data (.xlsx)
%
%% Output:
%   - Ranking classification (.csv)
%
%% Dependencies:
%
%   Files:
%       - Ranking data 9.xlsx)
%
%   Functions:
%       - None


%% Initialization

cd = 'D:\Work';
pathName = 'D:\Work';

%% Directory

data = xlsread('Ranking raw.xlsx',1);

%% Execution  

    % Remove the empty column between assessments
    data(:,4) = [];

    % Remove the empty rows between subjects
    for i=1:size(data,1)
        test = isnan(data(i,:));
        if sum(test) == 6
            idx(i) = i;
        else
            idx(i) = 0;
        end
    end
    idx(idx==0)=[];
    data(idx,:) = [];

    % Seperate the data from both assessments
    data1 = data(:,1:3);
    data2 = data(:,4:6);

    % Reshape data into a single row
    rowData1 = reshape(data1',1,size(data1,1)*size(data1,2));
    rowData2 = reshape(data2',1,size(data2,1)*size(data2,2));

    % Reshape data (ROW = subject; COLUMN = condition)
    subData1 = reshape(rowData1,9,size(data,1)/3)';
    subData2 = reshape(rowData2,9,size(data,1)/3)';

    % Attribute zeros for each condition that was not selected
    idx1 = find(isnan(subData1)==1);
    idx2 = find(isnan(subData2)==1);

    % Attribute the value 5 for each non-selected condition (blank on Excel)
    subData1(idx1) = 5;
    subData2(idx2) = 5;

    % Sum the ranking of both assessments
    sumData = subData1 + subData2;

%% insert functions here

% Option 1 Average Ranking
classData = AverageRanking(subData1, subData2);
% 
% % Option 2 Buy Yes or No
% classData = BuyYesNo(sumData);

% % Option 3 Equation of Ranking and Rating
% classData = RankingLikertEquation(cd,subData1,subData2,sumData);

%% Export

    xlswrite('Classification_N=81.xlsx',classData);
    
disp('     *****     DONE WITH CLASSIFICATION AND EXPORT     *****     ')