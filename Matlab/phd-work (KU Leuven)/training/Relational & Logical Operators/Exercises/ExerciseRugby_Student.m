%% Physiological characteristics of specific individual positions in junior 
% rugby league players
% Name:
% Date:

clear all
close all
clc

%% Step 1: Import a Excel file called 'RugbyPlayers.xlsx'
% Load the data form Excel file with name 'RugbyPlayers.xlsx'
...

% Compare if the number of columns and rows are the same for the data and 
% the headers. Use relational operators and the 'isequal' function. What is 
% the difference?
...

% Extract the data about the 'Position' of all the  players. To select the 
% correct column find 'Position' in the first row of the header file using 
% string compare and logical indexing.
...

%% Questions
% To solve the questions use relational and logical operators, string 
% compare and the find functions.
% 1) Calculate the average Height, Weight and VO2max of the players at 
% different positions. Select the position of players whose height, weight 
% as well as VO2max are above average. 
% The names of the columns are 'height (cm)', 'Body mass (kg)' and 
% 'VO2 max', respectively. Use 'strcmp' and logical indexing to select the
% data. Make sure you select the correct column!

% Height + average Height 
...

% Weight + average Weight
...

% VO2max + average VO2max
...

% Player positions (relational and logical operators + find function)
...

% 2) Use the weight en the height to calculate the BMI, with formula 
% weight/height², with weight in kg and height in meters. Is there a 
% player position in which the BMI index is above 27?

% BMI
...

% 3) Calculate the average speed in km/h of individual players over the 10m, 
% 20m and 40m sprint and calculate the average speed of all players over the 
% different distances. The data gives the time in seconds over a certain 
% distance. To convert the speeds from m/s to km/h use the following 
% relationship: km/h = 3.6*m/s. 
% Select the position of the players whose speed is lower than the average 
% speed of all the players over the different distances. 

% Average speed of individual players over 10, 20 and 40m sprint 
% + average speed of all players over the different distances.
...

% Player positions (relational and logical operators + find function)
...

% 4) Select the position of the players who jump-up at least 25% of their 
% own height in a vertical jump exercise.

% Jump + percentage jump with respect to height players
...

% Player positions (relational and logical operators + find function)
...