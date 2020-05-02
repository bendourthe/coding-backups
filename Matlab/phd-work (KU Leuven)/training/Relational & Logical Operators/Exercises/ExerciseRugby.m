%% Physiological characteristics of specific individual positions in junior 
% rugby league players
clear all,close all, clc

%% Step 1: Import a Excel file called 'RugbyPlayers.xlsx'
% Load the data form Excel file with name 'RugbyPlayers.xlsx'
[data, headers] = xlsread('Data/RugbyPlayers.xlsx');

% Compare if the number of columns and rows are the same for the data and 
% the headers. Use relational operators and the 'isequal' function. What is 
% the difference?
sizeEqual = (size(data) == size(headers));
sizeEqual2 = isequal(size(data),size(headers));
[rowData, colData] = size(data);
[rowHeader, colHeader] = size(headers);

% Extract the data about the 'Position' of all the  players. To select the 
% correct column find 'Position' in the first row of the header file using 
% string compare and logical indexing.  
Position = headers(2:end,strcmp('Position',headers(1,:)));

%% Questions
% To solve the questions use relational and logical operators, string 
% compare and the find functions.
% 1) Calculate the average Height, Weight and VO2max of the players at 
% different positions. Select the position of players whose height, weight 
% as well as VO2max are above average. 
% The names of the columns are 'height (cm)', 'Body mass (kg)' and 
% 'VO2 max', respectively. Use 'strcmp' and logical indexing to select the
% data. Make sure you select the correct column! 
Height = data(:,strcmp('Height (cm)',headers(1,:))-1);
aHeight = mean(Height);

Weight = data(:,strcmp('Body mass (kg)',headers(1,:))-1);
aWeight = mean(Weight);

VO2 = data(:,strcmp('VO2 max',headers(1,:))-1);
aVO2 = mean(VO2);

indQ1 = find((Height > aHeight) & (Weight > aWeight) & (VO2 > aVO2));
positionQ1 = Position(indQ1);

% 2) Use the weight en the height to calculate the BMI, with formula 
% weight/height², with weight in kg and height in meters. Is there a 
% player position in which the BMI index is above 27?
BMI = Weight./((Height/100).^2);
indQ2 = find(BMI > 27);
positionQ2 = Position(indQ2);

% 3) Calculate the average speed in km/h of individual players over the 10m, 
% 20m and 40m sprint and calculate the average speed of all players over the 
% different distances. The data gives the time in seconds over a certain 
% distance. To convert the speeds from m/s to km/h use the following 
% relationship: km/h = 3.6*m/s. 
% Select the position of the players whose speed is lower than the average 
% speed of all the players over the different distances.  
v_10 = 10./data(:,strcmp('10m (s)',headers(1,:))-1).*3.6;
v_20 = 20./data(:,strcmp('20m (s)',headers(1,:))-1).*3.6;
v_40 = 40./data(:,strcmp('40m (s)',headers(1,:))-1).*3.6;
aSpeed = [v_10,v_20,v_40];
aSpeedPlayers = mean(aSpeed);
indQ3 = find(aSpeed(:,1) < aSpeedPlayers(1) & ...
    aSpeed(:,2) < aSpeedPlayers(2) & ...
    aSpeed(:,3) < aSpeedPlayers(3));
positionQ3 = Position(indQ3);

% 4) Select the position of the players who jump-up at least 25% of their 
% own height in a vertical jump exercise.
Jump = data(:,strcmp('Vertical jump (cm)',headers(1,:)));
jump_perc = Jump.*100./Height;
indQ4 = find(jump_perc >= 25);
positionQ4 = Position(indQ4);