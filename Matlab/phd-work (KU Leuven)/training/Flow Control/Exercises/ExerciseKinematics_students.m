%% Exercise Kinematics
% Visualizing kinematic data using for and switch statements
clc
close all
clear all

%% Step 1:
% Open the file 'KS_motion_trial62.mot' with the kinematics values for the 
% simulated motion and plot the ankle, hip flexion, hip adduction and the 
% knee angles in four different subplots on 1 figure. Instead of doing this 
% as performed in exercise session II on plotting, try to do this as 
% efficient as possible, by using the FOR in combination with the SWITCH 
% command. Only one time the subplot command may be written.

% Import data from 'KS_motion_trial62.mot'
...
% Use the FOR statement to create the subplots. The SWITCH command is used 
% to determine the name of the column header to plot. For this dataset we 
% are interested in the following column names: 
% 'hip_flexion_r', 'hip_adduction_r', 'ankle_angle_r', 'knee_angle_r'.
%   o Find the column number by searching the matching column header in the 
%     imported column header data
%   o Create the subplot and plot the data

...

%% Step 2:
% Adapt the previous code, such that a title is added to each subplot, 
% which state the name of the variable (Hip flexion, Hip adduction, Knee 
% flexion and ankle flexion). Use an IF statement to add the xlabel 
% ('Time (s)') to the lower subplots and an ylabel ('Angle (°)') to the most 
% outer left subplots. the output is the same figure as in the Plotting 
% exercise session, but now implemented more efficient.

...% copy previous solution + add titles and axes labels %

%% Step 3:
% Repeat this workflow to open a second data file and plot the same 
% kinematic values on the corresponding subplot. Use a different color for 
% the different files. Implement this also as efficient as possible, using
% a nested for loop. Only the number in the filename is changed for the
% second file (KS_motion_trial69.mot, 69 instead of 62). Use this
% information to efficiently implement the FOR loop. Use the SWITCH command 
% to change the filename and the color of the line plot. 

...% copy previous solution + add new data file%
