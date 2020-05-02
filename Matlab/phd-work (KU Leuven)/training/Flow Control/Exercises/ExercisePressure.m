%% Plotting plantar pressure from a dynamical plantar pressure measurement

clear all
clc
close all

% file_name = 'ExercisePressure.m';
% function_dir = which(file_name);
% function_dir = function_dir(1:end-length(file_name));

data_folder = 'Data/'; 

file_name_pressure= 'Jos_4 - RS_L - Dynamic Roll off';

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

stop = 0;
frame_nr = 0;
begin_row = 10;

[DATA,TXT,RAW] = xlsread([data_folder, file_name_pressure]);

frame_height = DATA(7,1);
frame_width = DATA(5,1);

while stop == 0
    frame_nr = frame_nr + 1;
    
    begin_row = begin_row + frame_height + 2;
    
    if begin_row > length(DATA(:,1))
        stop = 1;
    else
        NUM = DATA(begin_row:begin_row+frame_height-1,:);
        frames(:,:, frame_nr) = NUM;
    end
end

%% Step 2: Use the 3D-matrix from the previous step to plot each frame. 
% Use the knowledge from Practical session 2. Plot the frame on the same
% figure (without hold on), such that feels like a movie, displaying the
% plantar pressure.

figure;
for frame_nr = 1:1:length(frames(1,1,:))
    x_plot = [0:1:size(frames(:,:,frame_nr), 2)-1]*delta_x;
    y_plot = [0:1:size(frames(:,:,frame_nr), 1)-1]*delta_y;
    
    surf(x_plot, y_plot, frames(:,:,frame_nr))
    view(0, 90)
    grid off
    axis off
    axis image
    
    drawnow
end

%% Step 3: Copy the code from step 1 and step 2 and let it run for 2
% different files. Do this also as efficient as possible.

file_names= {'Jos_4 - RS_L - Dynamic Roll off', 'Jos_4 - BMe_R - Dynamic Roll off'};
    
for file_nr = 1:1:length(file_names)
    if file_nr >1;
        clear frames;
    end
    figure;
    stop = 0;
    frame_nr = 0;
    begin_row = 10;
    
    file_name = file_names{file_nr};
    
    [DATA,TXT,RAW] = xlsread([data_folder, file_name]);
    
    frame_height = DATA(7,1);
    frame_width = DATA(5,1);
    
    while stop == 0
        frame_nr = frame_nr + 1;
        
        begin_row = begin_row + frame_height + 2;
               
        if begin_row > length(DATA(:,1))
            stop = 1;
        else
            NUM = DATA(begin_row:begin_row+frame_height-1,:);
            frames(:,:, frame_nr) = NUM;
        end
    end
    for frame_nr = 1:1:length(frames(1,1,:))
        x_plot = [0:1:size(frames(:,:,frame_nr), 2)-1]*delta_x;
        y_plot = [0:1:size(frames(:,:,frame_nr), 1)-1]*delta_y;
        
        surf(x_plot, y_plot, frames(:,:,frame_nr))
        view(0, 90)
        grid off
        axis off
        axis image
        
        drawnow
    end
end

%% Step 4: EXTRA: add some code to capture a movie from the plots made. 
% Therefore you need to you the VideoWriter object. Try to find more
% information in the help file.

file_names= {'Jos_4 - RS_L - Dynamic Roll off', 'Jos_4 - BMe_R - Dynamic Roll off'};
    
for file_nr = 1:1:length(file_names)
    figure;
    clear frames
    stop = 0;
    frame_nr = 0;
    begin_row = 10;
    
    file_name = file_names{file_nr};
    
    [DATA,TXT,RAW] = xlsread([data_folder, file_name]);
    
    frame_height = DATA(7,1);
    frame_width = DATA(5,1);
    
    while stop == 0
        frame_nr = frame_nr + 1;
        
        begin_row = begin_row + frame_height + 2;
               
        if begin_row > length(DATA(:,1))
            stop = 1;
        else
            NUM = DATA(begin_row:begin_row+frame_height-1,:);
            frames(:,:, frame_nr) = NUM;
        end
    end
    fps = 1/0.01; % frames per second
    WriterObj = VideoWriter([data_folder, file_name, '.avi']);
    WriterObj.FrameRate = fps;
    open(WriterObj)
    for frame_nr = 1:1:length(frames(1,1,:))
        x_plot = [0:1:size(frames(:,:,frame_nr), 2)-1]*delta_x;
        y_plot = [0:1:size(frames(:,:,frame_nr), 1)-1]*delta_y;
        
        surf(x_plot, y_plot, frames(:,:,frame_nr))
        set(gcf,'Renderer','zbuffer');
        view(0, 90)
        grid off
        axis off
        axis image
        
        drawnow
        
        if frame_nr == 1
            frame = getframe;
            rectangle = size(frame.cdata);
            rectangle = [0 0 rectangle(2)-1 rectangle(1)-1];
        end
        
        frame = getframe(gca, rectangle);
        writeVideo(WriterObj, frame);
    end
    close(WriterObj);
end