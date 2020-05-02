function classData = AverageRanking(subData1, subData2)    

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