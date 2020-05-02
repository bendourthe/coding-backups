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
...
% Create three variables to define the referencedate (eg. rd, rm, by)
...

% Computation of the age. Remember to take into account that someone may
% have had a birthday before the reference date. To solve this problem use 
% an IF statement. If the reference date has not been passed during the reference 
% year decrease the age with 1. Make sure no answer or zero is created if the
% reference day is prior to the birthday.
...

%% Step 2: Birthday database 1000 subjects
% Import the data from the Excel file named 'Birthdays.xls'. This file
% contains birthday data of 1000 subjects. Extract the day, month and year of
% birth using string compare command ('strcmp').
...

% Copy and adapt the previous code to calculate the age of 1000 subjects. 
% To do so add a FOR statement to execute the previous code for all the 
% subjects. 
...

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
...
% Construct a zeros-vector to include the number of subjects included in
% the group.
...
% Use a FOR loop over all subjects. Define the groups using SWITCH statement. 
% There are four groups in total numbered 1, 2, 3, and 4 for 
% 0-25, 26-50, 51-75, 76+ respectively.
...

% Optional exercise: Use a FOR loop over all subjects. Define the age groups 
% using IF-ELSEIF-ELSE statement. 
% There are four groups in total numbered 1, 2, 3, and 4 for 
% 0-25,26-50,51-75,76+ respectively.
...

% Construct a bar plot to visualize the distribution of the ages over the 
% four groups. Add a title and labels to the axes of the bar plot. 
figure, 
...
set(gca, 'XTickLabel', {'0-25'; '26-50';'51-75';'76+'}) % Command to change the labels of the ticks on the x-axis
...

% Add the constructed vector containing the age group and the computed ages 
% to the rest of the imported data, and add the column headers 'Age' and 
% 'AgeGroup'
...

% Write everything in a new excel file, called 'Birthdays_group.xls'
...
