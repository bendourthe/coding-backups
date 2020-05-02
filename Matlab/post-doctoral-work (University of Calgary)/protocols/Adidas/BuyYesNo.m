function classData = BuyYesNo(sumData)    

    %% Outcome

    % Classify the data and attribute the value 1 where both answers were
    % 'Yes' and 2 where at least one answer was 'No'
    classData = sumData;
    classData(classData<3) = 1;
    classData(classData>=3) = 0;