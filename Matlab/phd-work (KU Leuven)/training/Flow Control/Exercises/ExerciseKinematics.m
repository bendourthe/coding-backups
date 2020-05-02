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

% Import data from KS_motion_trial62.mot
[data] = importdata('Data/KS_motion_trial62.mot');
% FOR loop over de 4 different subplots
for subplot_nr=1:4
    % SWITCH to determine the the name of the columnheader
    % ('hip_flexion_r', 'hip_adduction_r', 'ankle_angle_r', 'knee_angle_r')
    switch subplot_nr
        case 1
            string = 'hip_flexion_r';
        case 2
            string = 'hip_adduction_r';
        case 3
            string = 'ankle_angle_r';
        case 4
            string = 'knee_angle_r';
    end
    
    % Find the column number by searching the matching column header in the 
    % imported column header data
    index = find(strcmp(data.colheaders, string));
    % Plot the signal on the subplot
    subplot(2, 2, subplot_nr)
    hold on; plot(data.data(:,1), data.data(:, index), 'r');
end


%% Step 2:
% Adapt the previous code, such that a title is added to each subplot, 
% which state the name of the variable (Hip flexion, Hip adduction, Knee 
% flexion and ankle flexion). Use an IF statement to add the xlabel 
% ('Time (s)') to the lower subplots and an ylabel ('Angle (°)') to the most 
% outer left subplots. the output is the same figure as in the Plotting 
% exercise session, but now implemented more efficient.

[data] = importdata('Data/KS_motion_trial62.mot');
for subplot_nr=1:1:4
    switch subplot_nr
        case 1
            string = 'hip_flexion_r';
            title_string = 'Hip flexion';
        case 2
            string = 'hip_adduction_r';
            title_string = 'Hip adduction';
        case 3
            string = 'ankle_angle_r';
            title_string = 'Ankle flexion';
        case 4
            string = 'knee_angle_r';
            title_string = 'Knee flexion';
    end
    
    index = find(strcmp(data.colheaders, string));
    subplot(2, 2, subplot_nr)
    hold on; plot(data.data(:,1), data.data(:, index), 'r');
    title(title_string)
    if subplot_nr == 1
        ylabel('Angle (°)');
    elseif subplot_nr == 3
        ylabel('Angle (°)');
        xlabel('Time (s)');
    elseif subplot_nr == 4
        xlabel('Time (s)');
    end
end

%% Step 3:
% Repeat this workflow to open a second data file and plot the same 
% kinematic values on the corresponding subplot. Use a different color for 
% the different files. Implement this also as efficient as possible, using
% a nested for loop. Only the number in the filename is changed for the
% second file (KS_motion_trial69.mot, 69 instead of 62). Use this
% information to efficiently implement the FOR loop. Use the SWITCH command 
% to change the filename and the color of the line plot. 

% FOR loop over the different files
for situation = 1:2
    % SWITCH to determine the filename
    switch situation
        case 1
            file_nr = 62;
            color = 'r';
        case 2
            file_nr = 69;
            color = 'b';
    end
    
    % Import the data from txt-file with the name filename as specified
    % above
    [data] = importdata(['Data/KS_motion_trial', num2str(file_nr), '.mot']);
    % FOR loop over the different subplots (can be copied from Step 2)
    for subplot_nr=1:1:4
        switch subplot_nr
            case 1
                string = 'hip_flexion_r';
                title_string = 'Hip flexion';
            case 2
                string = 'hip_adduction_r';
                title_string = 'Hip adduction';
            case 3
                string = 'ankle_angle_r';
                title_string = 'Ankle flexion';
            case 4
                string = 'knee_angle_r';
                title_string = 'Knee flexion';
        end
        
        index = find(strcmp(data.colheaders, string));
        subplot(2, 2, subplot_nr)
        hold on; plot(data.data(:,1), data.data(:, index), color);
        title(title_string)
        if subplot_nr == 1
            ylabel('Angle (°)');
        elseif subplot_nr == 3
            ylabel('Angle (°)');
            xlabel('Time (s)');
        elseif subplot_nr == 4
            xlabel('Time (s)');
        end
    end
end

