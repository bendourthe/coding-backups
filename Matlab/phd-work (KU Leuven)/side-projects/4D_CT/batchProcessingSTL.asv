 %% Batch processing STL files

%% C&C
close all;
clear all;

%% load data

% loading STL's from folder into struct "dataFolder"
folderFiles = 'C:\Users\u0078973\Desktop\4D 24-7-2013\STL Callibration Beads\001';
dataFolder = dir(folderFiles);

% Selecting correct files
startFile = 3; % place in dataFolder where bead A is located
endFile = 7; % place in dataFolder where bead E is located
counterFiles = 0;

for iFileNames = startFile : endFile;
    counterFiles = counterFiles + 1;
    correctFiles(counterFiles) = dataFolder(iFileNames);
end

for iFilename = 1 : 5
fpath = fullfile(pathname, iFilename);
end

% counter = 10;
% for iCorrectFiles = 1 : 5;
    
    fid = fopen('C:\Users\u0078973\Desktop\4D 24-7-2013\STL Callibration Beads\001', 'r');
    header = fread(fid, 84, '*char').';
    loc = regexpi(header, 'solid.*(facet)\s*normal', 'tokenExtents');
    counter = counter + 1;
% end

 
