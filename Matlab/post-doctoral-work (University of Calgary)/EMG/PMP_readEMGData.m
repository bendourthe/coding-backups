clc
clearvars
close all

%% load data from excel sheet
[~, ~, rawData] = xlsread('C:\Users\Fabian\Desktop\Master\CALGARY\Intern\PreferedMovementPath\Mizuno EMG good trials for Fabian.xls', 1);
subjects2use = xlsread('C:\Users\Fabian\Desktop\Master\CALGARY\Intern\PreferedMovementPath\Mizuno EMG good trials for Fabian.xls', 2);

%% find the idices of subjects to use in rawData
% convert subject numbers to numeric array
raw_idx = [NaN; cell2mat(rawData(2:end,2))];
% match subjects to use to index in rawData
[~, subject_idx] = ismember(subjects2use, raw_idx, 'R2012a');

%% clear data of unsed subjects
header = rawData(1,:);
info = rawData(subject_idx,:);

clearvars -except header info

%% create struct with subject info and emg data
data = struct;
for i = 1 : length(info)
    
    this_subject = info{i,1};
    % header of subject info
    data.(this_subject).header = header;
    % subject info
    data.(this_subject).info = info(i,:);
    
    % iterate over conditions (eg. shoe + running/ walking)
    % and load emg data of good trials
    for ii = 9 : 16
        condition = data.(this_subject).header{ii};
        % remove whitespaces
        condition = regexprep(condition,'\W','');
        
        % subject number
        sub = data.(this_subject).info{1};
        % subject gender and age
        genderAge = data.(this_subject).info{3};
        % current condition 
        cond = data.(this_subject).header{ii};
        whitespace_idx = find(cond == ' ');
        if 2 == length(whitespace_idx)
            pace = lower(cond(whitespace_idx(2) + 1 : end));
            shoe = lower(cond(whitespace_idx(1) + 1 : whitespace_idx(2) - 1));
        else
            pace = lower(cond(whitespace_idx(1) + 1 : end));
            shoe = lower(cond(1 : whitespace_idx(1) - 1));
        end
        % good trials
        trials = data.(this_subject).info{ii};
        trials = str2num(trials);
        
        % iterate over trials in this condition and load emg data
        for iii = 1 : length(trials)
            if iii < 10
                this_trial = ['t0' num2str(iii)];
            else
                this_trial = ['t' num2str(iii)];
            end
            % create fileName
            fileName = [sub '_' genderAge '_' pace '_' shoe '_' num2str(trials(iii)) '.log'];
            % create pathName
            pathName = ['C:\Users\Fabian\Desktop\Master\CALGARY\Intern\MIZUNO\EMG\' sub '_' genderAge];
            % load number of channels
            numChannels = read_channels(pathName, fileName);
            % change ending of fileName
            fileName(end-3 : end) = '.emg';
            % load emg data
            emgData = read_data(pathName, fileName, numChannels);
            % save emg data in struct
            data.(this_subject).(condition).(this_trial) = emgData;
            
            clearvars this_trial fileName pathName numChannels emgData
        end
        clearvars iii trials shoe pace whitespace_idx cond condition genderAge sub 
    end
    clearvars ii this_subject
end
clearvars i info header