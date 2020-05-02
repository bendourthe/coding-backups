%% Plotting plantar pressure from a dynamical plantar pressure measurement

clear all
clc
close all

%% Sensor Information
Total_width = 64;
Total_height = 64;
length_active_area = 0.488;
width_active_area = 0.325;
delta_y = length_active_area/Total_height; % 64 scanning lines
delta_x = width_active_area/Total_width; % 64 scanning lines

%% Step 1: Load the data from a dynamical pressure measurement, stored in 
% 'Jos_4 - RS_L - Dynamic Roll off.xls'. In this file a matrix is given for
% each time frame of the measurements. This matrix contains the pressure
% values for each sensor in the matrix. You can compare it with the data
% from an exercise from Practical Session 2. Read each frame and save this
% information in a large 3D-matrix: the 2 first dimensions for the the
% dimensions of the matrix, and the third one is the frame dimension. Do
% this as efficient as possible using a while loop...

file_name= 'Jos_4 - RS_L - Dynamic Roll off';

% ...


%% Step 2: Use the 3D-matrix from the previous step to plot each frame. 
% Use the knowledge from Practical session 2. Plot the frame on the same
% figure (without hold on), such that feels like a movie, displaying the
% plantar pressure.

% ...

%% Step 3: Copy the code from step 1 and step 2 and let it run for 2
% different files ('Jos_4 - RS_L - Dynamic Roll off.xls', 
% 'Jos_4 - BMe_R - Dynamic Roll off.xls'). Do this also as efficient as 
% possible.

file_names= {'Jos_4 - RS_L - Dynamic Roll off', 'Jos_4 - BMe_R - Dynamic Roll off'};
    
% ....

%% Step 4: EXTRA: add some code to capture a movie from the plots made. 
% Therefore you need to you the VideoWriter object. Try to find more
% information in the help file.

% ....
% 
% WriterObj = VideoWriter([data_folder, file_name, '.avi']);
% WriterObj.FrameRate = fps;
% open(WriterObj)
%
% figure
% frame = getframe(gca, rectangle);
% writeVideo(WriterObj, frame);
% close(WriterObj);
       