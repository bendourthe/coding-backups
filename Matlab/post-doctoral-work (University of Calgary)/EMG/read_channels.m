function out = read_channels(pathname,filename);

logname = [filename(1:end-3) 'log'];             % check if a logfile exists for desired .emg file
logname = fullfile(pathname,logname);
if(~exist(logname, 'file'))
    errordlg(['File: '' ' logname, ''' does not exist'])
    data_matrix = 0;
    parameters = 0;
    logfile_cells = 0;
    return
end
log_file_id = fopen(logname, 'r');                  % open the logfile

i = 1;
while ~feof(log_file_id)                            % read the whole logfile.
    logfile_cells{i} = {fgets(log_file_id)};        % save to logfile_cells variable
    i = i + 1;
end
fclose(log_file_id);
clear log_file_id;

for i = 1:length(logfile_cells)                     % convert the logfile to characters
    logfile_cells_disp = char(logfile_cells{i});
    %     disp(logfile_cells_disp(1:end-1));
end

logfile_cells_disp = char(logfile_cells{2});
columns = str2num(logfile_cells_disp(1:end-1));
channels = columns-1;
out = channels;
end