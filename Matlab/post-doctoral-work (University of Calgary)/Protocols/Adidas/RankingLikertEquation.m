function classData = RankingLikertEquation(cd,subData1,subData2,sumData)  

    rating = xlsread('Rating raw.xlsx',1);
    
%% Execution  

    % Remove the empty column between assessments
    rating(:,4) = [];

    % Remove the empty rows between subjects
    for i=1:size(rating,1)
        test = isnan(rating(i,:));
        if sum(test) == 6
            idx(i) = i;
        else
            idx(i) = 0;
        end
    end
    idx(idx==0)=[];
    rating(idx,:) = [];

    % Seperate the data from both assessments
    rating1 = rating(:,1:3);
    rating2 = rating(:,4:6);

    % Reshape data into a single row
    rowRating1 = reshape(rating1',1,size(rating1,1)*size(rating1,2));
    rowRating2 = reshape(rating2',1,size(rating2,1)*size(rating2,2));

    % Reshape data (ROW = subject; COLUMN = condition)
    subRating1 = reshape(rowRating1,9,size(rating,1)/3)';
    subRating2 = reshape(rowRating2,9,size(rating,1)/3)';
    
%     % Invert Likert Scale
%     invRating1 = (subRating1-11)*-1;
%     invRating2 = (subRating2-11)*-1;
    
    % Normalize Rating to used scale ranke
    rangeRating1 = max(subRating1') - min(subRating1');
    rangeRating2 = max(subRating2') - min(subRating2');
    
    rangeRating1 = rangeRating1';
    rangeRating2 = rangeRating2';
    
    normRating1 = abs((subRating1-(max(subRating1'))')./rangeRating1);
    normRating2 = abs((subRating2-(max(subRating2'))')./rangeRating2);

    % Sum the ranking of both assessments
    sumRating = (normRating1 + normRating2)/2;
    
    % Build Equation
    equation = (sumRating.*sumData);
    
    % Exclude cases where range is too small
    % range should be greater/equal the average median range
    range1 = median(rangeRating1,2);
    range2 = median(rangeRating2,2);
    
    mean_range1 = mean(range1);
    mean_range2 = mean(range2);

    for le = 1:size(equation,1)
        if rangeRating1(le) <= mean_range1
            equation(le,:) = nan;
        end
        if rangeRating2(le) <= mean_range2
            equation(le,:) = nan;
        end
    end
    
    %% Outcome

    % Set threshold as mean of median result for each of the 9 shoes.
    thresholds = nanmedian(equation,1);
    mean_threshold = mean(thresholds);
    
    % Classify the data and attribute the value 1 where the result is below
    % the threshold and 0 where the result is above the threshold
    classData = equation;
    classData(classData<mean_threshold) = 1;
    classData(classData>=mean_threshold) = 0;