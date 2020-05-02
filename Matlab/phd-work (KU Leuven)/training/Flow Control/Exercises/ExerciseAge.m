%% ExerciseAge
clear all, 
close all, 
clc
%% Step 1: Determine the age of a person on a reference date
% The birthdate of a person is given by three integers bd, bm, by
% representing the day, month and year of birth. Make a script using the 
% if-statement which determines the age in years of a person at a certain 
% reference date. The reference date is represented by three integers rd, rm,
% ry (day, month, year). Remember to take into account that someone may have 
% had a birthday before the reference date. Test the script for multiple dates.

% Create three variables to define the birthdate (eg. bd, bm, by)
bd = 7;     % birth day
bm = 2;     % birth month
by = 1995;  % birth year
% Create three variables to define the referencedate (eg. rd, rm, by)
rd = 23;    % reference day
rm = 9;     % reference month
ry = 2013;  % reference year

% Computation of the age. Remember to take into account that someone may
% have had a birthday before the reference date. To solve this problem use 
% an IF statement. If the reference date has not been passed during the reference 
% year decrease the age with 1. Make sure no answer or zero is created if the
% reference day is prior to the birthday.
age = ry - by;
if age < 0
    disp('Not possible reference day prior to birthday')
    age = [];
elseif( bm > rm || ( bm == rm && bd > rd ) )
    % If the reference date has not been passed this year decrease the age with 1
    age = age - 1;
end;

%% Step 2: Birthday database 1000 subjects
% Import the data from the Excel file named 'Birthdays.xls'. This file
% contains birthday data of 1000 subjects. Extract the day, month and year of
% birth using string compare command ('strcmp').
[data,colheaders] = xlsread('Data\Birthdays.xls');
day = data(:,strcmp('day',colheaders));
month = data(:,strcmp('month',colheaders));
year = data(:,strcmp('year',colheaders));

% Copy and adapt the previous code to calculate the age of 1000 subjects. 
% To do so add a FOR statement to execute the previous code for all the 
% subjects. 
ages = zeros(size(year));
for i = 1:length(year)
    ages(i) = ry - year(i);
    if ages(i) < 0
        disp(['Not possible reference day prior to birthday for subject ', num2str(i)])
        ages(i) = 0;
    elseif( month(i) > rm || ( month(i) == rm && day(i) > rd ) )
        % If the reference date has not been passed this year decrease the age with 1
        ages(i) = ages(i) - 1;
    end;
end

%% Step 3: Population groups
% To visualize the age distribution for our birthday data we have to create
% age classes. For this exercise we want to define four age groups 
% (0-25, 26-50, 51-75, 76+). Visualize the distribution over the groups using 
% a bar plot. Save the group number and the age in an excel file 
% as a separated column. For all the age groups we want to know the number 
% of subjects included in the age group. Use a FOR loop to execute the 
% statements for all the subjects. To create the groups use the 
% SWITCH statement. For the definition of the case expressions use a cell 
% array. 
% (It is also possible to create the age groups using an IF statement. As
% and exercise you could program this as well)

% Construct a zeros-vector to store the group number for every subject.
population_group = zeros(size(year));
% Construct a zeros-vector to include the number of subjects included in
% the group.
group = zeros(1,4);
% Use a FOR loop over all subjects. Define the groups using SWITCH statement. 
% There are four groups in total numbered 1, 2, 3, and 4 for 
% 0-25, 26-50, 51-75, 76+ respectively.
for i = 1:length(year)
    switch ages(i)
        case num2cell(0:25)
            group(1) = group(1)+1;
            population_group(i) = 1;
        case num2cell(26:50)
            group(2) = group(2)+1;
            population_group(i) = 2;
        case num2cell(51:75)
            group(3) = group(3)+1;
            population_group(i) = 3;
        otherwise
            group(4) = group(4)+1;
            population_group(i) = 4;
    end
end

% Optional exercise: Use a FOR loop over all subjects. Define the age groups 
% using IF-ELSEIF-ELSE statement. 
% There are four groups in total numbered 1, 2, 3, and 4 for 
% 0-25,26-50,51-75,76+ respectively.
population_group1 = zeros(size(year));
group1 = zeros(1,4);
for i = 1:length(year)
    if ages(i) <= 25
        group1(1) = group1(1)+1;
        population_group1(i) = 1;
    elseif ages(i) <= 50
        group1(2) = group1(2)+1;
        population_group1(i) = 2;
    elseif ages(i) <= 75
        group1(3) = group1(3)+1;
        population_group1(i) = 3;
    else
        group1(4) = group1(4)+1;
        population_group1(i) = 4;
    end
end

% Construct a bar plot to visualize the distribution of the ages over the 
% four groups. Add a title and labels to the axes of the bar plot. 
figure, 
bar(group)
set(gca, 'XTickLabel', {'0-25'; '26-50';'51-75';'76+'}) % Command to change the labels of the ticks on the x-axis
xlabel('Age groups')
ylabel('Number of subjects')
title('Age distribution of 1000 subjects')

% Add the constructed vector containing the age group and the computed ages 
% to the rest of the imported data, and add the column headers 'Age' and 
% 'AgeGroup'
colheaders{end+1} = 'Age';
colheaders{end+1} = 'AgeGroup';
data = [data, ages, population_group];

% Write everything in a new excel file, called 'Birthdays_group.xls'
xlswrite('Data\Birthdays_group.xls', colheaders, 'A1:E1')
xlswrite('Data\Birthdays_group.xls', data, 'A2:E1001')
