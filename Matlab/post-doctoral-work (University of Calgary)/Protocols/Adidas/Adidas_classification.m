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

cd 'C:\Users\bdour\OneDrive\Work\Calgary\Projects\Adidas\';
pathName = 'C:\Users\bdour\OneDrive\Work\Calgary\Projects\Adidas';

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

    % Identify cases where one assessment ranked a shoe as most favourite 
    % and the other assessment ranked the same shoe as 3rd least favourite 
    % (any combination of 1 and 7)
    idx17 = find(subData1==1 & subData2==7 | subData1==7 & subData2==1);
    
    % Attribute the value 9 for each 1&7 case (note: change that value to 8 
    % if you wish to include these cases in the grouping)
    sumData(idx17) = 9;
    
    % Identify cases when one assessment ranked a shoe as most favourite
    % and the other assessment ranked the same shoe as least favourite
    % (any combination of 1 and 9 - shows the subject cannot rank)
    [row19,col19] = find(subData1==1 & subData2==9 | subData1==9 & subData2==1);
    
    % Attribute the value 9 for the whole row of each 1&9 case
    % (i.e. unclassifiable)
    sumData(row19,:) = 9;

%% Outcome

    % Classify the data and attribute the value 1 for each group that the
    % corresponding subject belong to
    classData = sumData;
    classData(classData<9) = 1;
    classData(classData>=9) = 0;
    
%% Export

    xlswrite('Classification.xlsx',classData);
    
disp('     *****     DONE WITH CLASSIFICATION AND EXPORT     *****     ')