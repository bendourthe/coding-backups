% 0. clear workspace and figures (not required)
clear all
close all

%1. load file
load('knee_moment_sts.mat')

%2. visualise
figure()                % open an new figure
plot(one_trial_moment); % you will learn more about this function in the next practical session

%3. peak flexion and extension
[peak_extension, peak_extension_percentage] = max(one_trial_moment);
[peak_flexion, peak_flexion_percentage] = min(one_trial_moment);
peak_flexion = peak_flexion.*-1;
% multiply by -1 because flexion is in the negative direction

%4. normalise to body mass
one_trial_moment = one_trial_moment./60;

%5. plot multiple trials
figure()    % open an new figure
plot(multiple_trials_moments); % you will learn more about this function in the next practical session


%6. max values + location
max_values_trials = max(multiple_trials_moments);
[peak_value, trial_peak_value] = max(max_values_trials); % second output argument max function 

%7. standard deviation max values
std_max_values = std(max_values_trials);

% 8 mean and std moments time
mean_moment = mean(multiple_trials_moments,2);% second dimension to calculate the mean on each percentage of gait cycle (mean along row)
std_moment = std(multiple_trials_moments,0,2);
