%% Exercise on the classification of patients based on the BMI
clc
close all
clear all

%% Step 1: Open the excel file which contains data on the height, weight and
% BMI of 100 subjects. Make a scatter plot of the weight vs. the height of
% the person. Make a distinction between male and female, which is defined
% by the gender column (1 for male, 2 for female). How many males and
% females where included in the study?

% Import the data from the excel file named Height_weight_BMI_stats.xlsx,
% which has been made in exercise session on data import and export
...

% Extract the data from the Gender, Height (cm), Weight (kg) and BMI column
% using logical indexing and string compare command ('strcmp')
...

% find the indexes corresponding to males and those corresponding to the
% females, using the FIND command (1 for male, 2 for female)
...

% How many males and females were there included in the study?
...

% Make a scatter plot of the weight against the height. Use a different
% marker symbol and marker color (filled) for the males and females. Add
% also a meaningful title, labels and a legend to the plot
...

% Save the plot under the name weight_vs_height.jpg
...

%% Step 2: Classify the subjects according the BMI. The classes are:
% underweight (<18.5), Normal range (18.5-25), overweight (25-30) and obese
% (>30). Save the number of classes in the excel file as a separated
% column. 

% Construction of a zero-vector classes
...

% FOR loop over all subjects
% Classification using IF. Use classes numbered 1, 2, 3 and 4
... 

% Add the constructed vector to the rest of the imported data, and add a
% column header 'Classes'
...

% Write everything in a new excel file, called Height_weight_BMI_classes.xlsx
...

%% Step 3: Construct a stacked bar plot, two separated one for the females 
% and males.
% The height of the bar plot indicates the percentage of male and females
% in the different BMI classes. Use therefor the classes as determined in
% previous step.

% Construct a zeros-matrix for the bar plot with two rows (two groups: male 
% and female) and 4 columns (4 BMI classes - underweight, normal range, 
% overweight and obese)
... 

% Nested FOR loop to run over all gender and all classes. For each
% combination (gender and classes), determine the number of corresponding
% subjects (using the sum and relational operators). Fill the bar plot
% matrix with this number of subjects. 
...

% Normalize the barplot matrix using the number of females and males, such
% that the total of all classes in 100%.
...

% Plot the matrix using the barplot function, using a stacked layout. If
% everything is alright, the largest Y-value will be 100%. Add a legend to
% the plot, which states the different classes. Look to the given example of
% this barplot.
figure,
...
set(gca, 'XTickLabel', {'Male'; 'Female'}) % Command to change the labels of the ticks on the x-axis, to 'Male' and 'Female'
...

% Save the figure under the name BMI_classification.png
...