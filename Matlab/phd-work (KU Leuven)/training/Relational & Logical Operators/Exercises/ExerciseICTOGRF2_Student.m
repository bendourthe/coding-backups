%% GRF detect IC/TO using relational operators and find command. 
% Compare the IC/TO with the moments of IC and TO found in 
% 'ExerciseICTOGRF'.
clc
close all
clear all

%% Step 1: Import a txt-file called TreadmillGRF.mot.
% Import the txt-file called TreadmillGRF.mot. This file contains measured 
% GRF during a measurement of a walk on a treadmill. Plot the vertical force 
% of the left and right leg. Plot the time on the x-axis. Add meaningful 
% labels to the axis, legend and the title of the plot. Use a line-width of 
% 1.5 for the curves and a front-size of 14 for the title and labels. Plot 
% the forces (left - black, right – magenta), and change the line-style.
% The vertical forces for the left and right leg are indicated with a 
% colheader '1_ground_force_vy' and 'ground_force_vy', respectively. The 
% colheader indicating the time is 'time'.  
  
% Import data
...

% Extract the time, vertical force of the right and left leg. Find the 
% column numbers in the data struct using logical indexing and string 
% compare 'strcmp'.
...

% Plot the forces (left - black, right - magenta. Change the line-style and 
% add label, title and a legend. Adapt the xlimit if necessary.
...

%% Step 2: Automatic initial contact and toe-off detection

% Detect the Initial Contact (IC) and Toe Off (TO) of the gait cycles 
% automatically. 
% Create a logical array to detect the moments in which there is contact with 
% the force plate is true, no contact is false. Use a threshold of 0.5 to 
% avoid selecting the noise in the signals.
...

% Determine the indices (points in the matrix) of initial contact (IC). 
% These are the indices for which the difference with the next logical 
% number is +1. Use the 'diff' and 'find' commands.
...

% Determine the indices (points in the matrix) toe off (TO). These are the 
% indices for which the difference with the next logical number is -1. 
% Check the definition of 'diff' function.
...
    
% Plot the IC and TO on the previous graph using marker symbols. Use
% different symbols to separate IC and TO. Use a different color for the 
% left and right foot. Adapt the legend. 
...

% Inspect the figure carefully. Are there any errors in detecting the IC
% and TO times? Use the interactive handles of the plot. 

% Put the IC and TO times in a structured array ICTO2
...

% Save the structured array in a mat-file.
...

%% Step 3: Compare the IC and TO times with manual selected timings

% Load the saved IC and TO times of the plotting exercise in practical 
% session II, in which you defined the IC and TO manually. This information 
% is saved in a structured array named ICTO
...

% Draw vertical lines on the instances of IC and TO. Use a different
% line style for IC and TO separately. Use the same color as the markers 
% of the left and right foot. 
ylimit = get(gca,'Ylim');   % Get the y-limit of the plot
ylim(ylimit);
...
    
% Adapt the xlimit of the figure to display the previous selected gait
% cycles using the xlim function.
...

% Save the created figure in a png-format and in a fig-format. 
...

% Inspect the figure carefully, are the manual selected IC and TO moments
% corresponding with the automatically determined times? Use the find
% command and relational and logical operators to check automatically.
...