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
data = importdata('Data/TreadmillGRF2.mot');

% Extract the time, vertical force of the right and left leg. Find the 
% column numbers in the data struct using logical indexing and string 
% compare 'strcmp'.
time = data.data(:,strcmp('time',data.colheaders));
GRF_left = data.data(:,strcmp('1_ground_force_vy',data.colheaders));
GRF_right = data.data(:,strcmp('ground_force_vy',data.colheaders));

% Plot the forces (left - black, right - magenta. Change the line-style and 
% add label, title and a legend. Adapt the xlimit if necessary.
figure, hold on
plot(time, GRF_right, 'm', 'LineWidth', 1.5);
plot(time, GRF_left, 'k', 'LineWidth', 1.5);
ylabel('Force (N)','FontSize',14); 
xlabel('Time (s)','FontSize',14);
title('Vertical force during treadmill walking', 'FontSize', 14)
legend({'Right', 'Left'})
range = [64 67]; xlim(range)


%% Step 2: Automatic initial contact and toe-off detection 

% Detect the Initial Contact (IC) and Toe Off (TO) of the gait cycles 
% automatically. 
% Create a logical array to detect the moments in which there is contact with 
% the force plate is true, no contact is false. Use a threshold of 0.5 to 
% avoid selecting the noise in the signals.
threshold = 0.5;
logicGRF_left = abs(GRF_left) > threshold;
logicGRF_right = abs(GRF_right) > threshold;

% Determine the indices (points in the matrix) of initial contact (IC). 
% These are the indices for which the difference with the next logical 
% number is +1. Use the 'diff' and 'find' commands.
indIC_left = find(diff(logicGRF_left)==1);
%indIC_left2 = find((logicGRF_left(2:end)-logicGRF_left(1:end-1))==1);
indIC_right = find(diff(logicGRF_right)==1);
%indIC_right2 = find((logicGRF_right(2:end)-logicGRF_right(1:end-1))==1);

% Determine the indices (points in the matrix) toe off (TO). These are the 
% indices for which the difference with the next logical number is -1. 
% Check the definition of 'diff' function.
indTO_left = find(diff(logicGRF_left)==-1)+1;
%indTO_left2 = find((logicGRF_left(2:end)-logicGRF_left(1:end-1))==-1)+1;
indTO_right = find(diff(logicGRF_right)==-1)+1;
%indTO_right2 = find((logicGRF_right(2:end)-logicGRF_right(1:end-1))==-1)+1;

% Plot the IC and TO on the previous graph using marker symbols. Use
% different symbols to separate IC and TO. Use a different color for the 
% left and right foot. Adapt the legend. 

plot(time(indIC_left),GRF_left(indIC_left),'or','MarkerSize',10)
plot(time(indIC_right),GRF_right(indIC_right),'ob','MarkerSize',10)
plot(time(indTO_left),GRF_left(indTO_left),'*r','MarkerSize',10)
plot(time(indTO_right),GRF_right(indTO_right),'*b','MarkerSize',10)
legend({'Right', 'Left','IC_L','IC_R','TO_L','TO_R'})

% Inspect the figure carefully. Are there any errors in detecting the IC
% and TO times? 

% Put the IC and TO times in a structured array ICTO2
IC.Left = time(indIC_left);
IC.Right = time(indIC_right);
TO.Left = time(indTO_left);
TO.Right = time(indTO_right);

ICTO2.IC = IC;
ICTO2.TO = TO;

% Save the structured array in a mat-file.
save('Data/ICTO2', 'ICTO2')

%% Step 3: Compare the IC and TO times with manual selected timings

% Load the saved IC and TO times of the plotting exercise in practical 
% session II, in which you defined the IC and TO manually. This information 
% is saved in a structured array named ICTO
load('Data/ICTO');

% Draw vertical lines on the instances of IC and TO. Use a different
% line style for IC and TO separately. Use the same color as the markers 
% of the left and right foot. 
ylimit = get(gca,'Ylim');   % Get the y-limit of the plot
%ylimit = [0 800];
plot([ICTO.IC.Left; ICTO.IC.Left], [ylimit; ylimit]', '--r');
%line([ICTO.IC.Left; ICTO.IC.Left],[ylimit; ylimit]','Color','r','LineStyle','--');
plot([ICTO.IC.Right; ICTO.IC.Right], [ylimit; ylimit]', '--b');
%line([ICTO.IC.Right; ICTO.IC.Right],[ylimit; ylimit]','Color','b','LineStyle','--');

plot([ICTO.TO.Left; ICTO.TO.Left], [ylimit; ylimit]', ':r');
%line([ICTO.TO.Left; ICTO.TO.Left],[ylimit; ylimit]','Color','r','LineStyle',':');
plot([ICTO.TO.Right; ICTO.TO.Right], [ylimit; ylimit]', ':b');
%line([ICTO.TO.Right; ICTO.TO.Right],[ylimit; ylimit]','Color','b','LineStyle',':');

% Adapt the xlimit of the figure to display the previous selected gait
% cycles using the xlim function
xlimit = [64 67];
%xlim(xlimit)

% Save the created figure in a png-format and in a fig-format. 
saveas(gcf, 'Data/GRF_ICTO_detection2.png');
saveas(gcf, 'Data/GRF_ICTO_detection2','fig');

% Inspect the figure carefully, are the manual selected IC and TO moments
% corresponding with the automatically determined times? Use the find
% command and relational and logical operators to check automatically.
ind_checkIC_Left = find((ICTO2.IC.Left == ICTO.IC.Left(1))|(ICTO2.IC.Left == ICTO.IC.Left(2)));
ind_checkIC_Right = find((ICTO2.IC.Right == ICTO.IC.Right(1))|(ICTO2.IC.Right == ICTO.IC.Right(2)));
ind_checkTO_Left = find((ICTO2.TO.Left == ICTO.TO.Left(1))|(ICTO2.TO.Left == ICTO.TO.Left(2)));
ind_checkTO_Right = find((ICTO2.TO.Right == ICTO.TO.Right(1))|(ICTO2.TO.Right == ICTO.TO.Right(2)));