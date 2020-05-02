%
% %%%%%%%%%%%%%  PEDAR ANALYSIS  %%%%%%%%%%%%%%%%%%%%%%%%% June 21, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This program will load raw(.asc) and force(.fgt) PEDAR files to anaylize
% COP (ML & AP), peak force, and partitians the insole to different regions

% Created by Chris Lam, using code modified by Colin Firminger
% 'HPL_data_analysis_v12'

% Start Fresh
close all;clc;clear all;


%% Allocate Folder Directories


% Allocate directory
directory = '\\CINCO\students\BennoGroup\MMohr\Running Style\';

% Pedar file path
pedarfolder = 'Pedar\';

% Raw Data file path
rawdatafolder = 'RawPedarforProcessing\';

% Save Folder
savefolder  = '\\CINCO\students\BennoGroup\MMohr\Running Style\Pedar\ChrisLamProcessing\';

% Raw Pressure .mat folder
rawpressurefolder = 'C:\Users\chris.lam1\Documents\MATLAB\Running Style Raw Pressure\';

% Processed Data
processedfolder = '\\CINCO\students\BennoGroup\MMohr\Running Style\Pedar\ChrisLamProcessing\COMPLETE\';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Processing Raw .fgt and .asc pedar files

%Read in insole area files
insoleareas = 'Pedar insole cell areas.xlsx';
insolearea_VScells = xlsread([directory pedarfolder insoleareas],'Sheet1','B2:B100');
insolearea_WScells = xlsread([directory pedarfolder insoleareas],'Sheet1','C2:C100');
insolearea_XScells = xlsread([directory pedarfolder insoleareas],'Sheet1','D2:D100');
insolearea_YScells = xlsread([directory pedarfolder insoleareas],'Sheet1','E2:E100');

%Read in centroid data for each Pedar insole (VS, WS, XS, YS) in units of mm
centroidfile = 'Centroid Calculation_All Insoles.xlsx';
V_centroids = xlsread([directory pedarfolder centroidfile],'VS Insole','O2:P100');
W_centroids = xlsread([directory pedarfolder centroidfile],'WS Insole','O2:P100');
X_centroids = xlsread([directory pedarfolder centroidfile],'XS Insole','O2:P100');
Y_centroids = xlsread([directory pedarfolder centroidfile],'YS Insole','O2:P100');
V_cent_max = [0.0701 0.2307]; %m
W_cent_max = [0.0748 0.2422]; %m
X_cent_max = [0.0782 0.257]; %m
Y_cent_max = [0.0828 0.2687]; %m


% Force Step Threshold
thresh = 100; %newtons

% Force Step Detection Threshold
threshold = 1000; %newtons

%%% Subject loop %%%

for s = 1:40 %40 subjects
    
    % For naming convention
    if s < 10
        subjnum = sprintf('Runner00%d',(s));
        subjnum1 = sprintf('r00%d',(s));
    else
        subjnum = sprintf('Runner0%d',(s));
        subjnum1 = sprintf('r0%d',(s));
    end
    
    
    % Load testing sheet to determine correct insole size from data collection sheet
    subjquestionnaire = [directory '\Data Collection\20160815_Running Style Questionnaire Results.xlsx'];
    
    % Find exact cell to determine the insole
    aaa = (2*s)+1;
    shoecell = sprintf('BE%d',(aaa));
    sheet = 'Sujective Questionaire Results';
    [var,insole] = xlsread(subjquestionnaire,sheet,shoecell);
    insole = cell2mat(insole);
    
    %Setting up correct insole size
    if strcmpi(insole,'beige')
        insolearea_total = 0.013931; %m^2
        insolearea_cells = insolearea_VScells;
        subj_cent = V_centroids;
        max_cent = V_cent_max;
        insole_length = 239.1; %mm
        insole_width = 75.6; %mm
    elseif strcmpi(insole,'blue')
        insolearea_total = 0.015542; %m^2
        insolearea_cells = insolearea_WScells;
        subj_cent = W_centroids;
        max_cent = W_cent_max;
        insole_length = 251.6; %mm
        insole_width = 81.2; %mm
    elseif strcmpi(insole,'black')
        insolearea_total = 0.017208; %m^2
        insolearea_cells = insolearea_XScells;
        subj_cent = X_centroids;
        max_cent = X_cent_max;
        insole_length = 267.2; %mm
        insole_width = 85; %mm
    elseif strcmpi(insole,'red')
        insolearea_total = 0.019113; %m^2
        insolearea_cells = insolearea_YScells;
        subj_cent = Y_centroids;
        max_cent = Y_cent_max;
        insole_length = 279.1; %mm
        insole_width = 89.9; %mm
    end
    
    %Creating v, w, x and y centroid matrices
    subj_centx = subj_cent(:,1)';
    subj_centy = subj_cent(:,2)';
    cent_matx = repmat(subj_centx,100,1); %
    cent_maty = repmat(subj_centy,100,1);
    
    
    %%%% Condition Loop %%%%
    
    % 5 shoe conditions
    for u = 1:5
        condnum = sprintf('_c%d',(u));
        condfolder = [subjnum1 condnum '_SD_'];
        errorfoldername = ['r001' condnum '_SD_'];
        
        % Generate array to hold 15 individual steps for only 2 trials
        steparrayforce = zeros(100,1); % 100 rows(norm data points), 31 col(steps)
        steparraycopxcalc = zeros(100,1); % 100 rows(norm data points), 31 col(steps) calculated through asc
        steparraycopycalc = zeros(100,1); % 100 rows(norm data points), 31 col(steps) calculated through asc
        ascholdingarray = zeros(100,99); %create array to hold individual steps to multiply by centroid
        ascholdingarray2 = zeros(100,99); %create array to hold individual steps to multiply by centroid
        pressurematrix = zeros(100,99,3); %matrix to hold raw pressure data for all 99 sensors
        indipressure = zeros(100,99); %array for avg pressure
        
        y = 0;
        
        % This section is only for the 4 or 5 conditions that have 2 trials
        if strcmpi(condfolder,'r007_c3_SD_')|| strcmpi(condfolder,'r035_c5_SD_')|| strcmpi(condfolder,'r030_c4_SD_')|| strcmpi(condfolder,'r028_c2_SD_')||strcmpi(condfolder,'r021_c1_SD_')|| strcmpi(condfolder,'r016_c2_SD_')|| strcmpi(condfolder,'r013_c1_SD_')|| strcmpi(condfolder,'r012_c1_SD_')|| strcmpi(condfolder,'r003_c3_SD_')|| strcmpi(condfolder,'r010_c1_SD_')|| strcmpi(condfolder,'r009_c2_SD_')|| strcmpi(condfolder,'r009_c3_SD_')|| strcmpi(condfolder,'r004_c5_SD_')|| strcmpi(condfolder,'r024_c4_SD_')||strcmpi(condfolder,'r002_c4_SD_')||strcmpi(condfolder,'r003_c1_SD_')|| strcmpi(condfolder,'r003_c2_SD_')|| strcmpi(condfolder,'r004_c4_SD_')|| strcmpi(condfolder,'r006_c3_SD_')
            
            %%%% Trial loop %%%%
            
            % 2 running trials per condition
            for t = 1:2
                
                % files are labelled 2 and 3
                tt = sprintf('_%d',(t+1));
                
                % Load Files
                filedirectory = [directory rawdatafolder subjnum '\' errorfoldername '\' condfolder tt];
                
                tempcondfold = [condfolder tt];
                
                % Reading in each .fgt pedar file (force data)
                pedar_fgt_filename = [filedirectory '.fgt'];
                % Reading pedar force data for the whole insole at each timestep from fgt file
                pedar_fgt_dataraw = dlmread([pedar_fgt_filename],'\t',9,5);
                pedar_fgt_data = pedar_fgt_dataraw(:,1);
                
                for i = 1:size(pedar_fgt_data)
                    if pedar_fgt_data(i) < 25
                        pedar_fgt_data(i) = 0;
                    end
                end
                
                
                plot(pedar_fgt_data(:,1))
                disp('Click onset and offset (x2)')
                f = ginput(2);
                newfgtdata = pedar_fgt_data(f(1,1):f(2,1),:); %clip the fgt data to boundaries set by cursor
                newfgtdataoriginal = newfgtdata;
                
                clear pedar_fgt_dataf pedar_fgt_datax pedar_fgt_datay;
                
                
                % Reading in each .asc pedar file (force data)
                pedar_asc_filename = [filedirectory '.asc'];
                % Reading pedar force data for the whole insole at each timestep from fgt file
                pedar_asc_data = dlmread([pedar_asc_filename],'\t',10,100); %'/t' tab delimited, starting at row 10, column 100
                if size(pedar_asc_data,2) == 100
                    pedar_asc_data(:,100) = [];
                end
                
                % Clip .asc data to window determined by force data
                newascdata = pedar_asc_data(f(1,1):f(2,1),:);
                newascdataoriginal = newascdata;
                
                % Issue with Pedar equipment, activates to 20-30 Hz
                % Resetting all pressure cells measuring less than 30 kPa to 0
                for i = 1:size(newascdata,1)
                    for j = 1:size(newascdata,2)
                        if newascdata(i,j) < 30
                            newascdata(i,j) = 0;
                        end
                    end
                end
                
                
                %Finding the peak of each step
                [peaks,locs] = findpeaks(newfgtdata,'MinPeakDistance',100,'MinPeakHeight',threshold);
                
                x = size(peaks,1);
                
                
                for h = 1:x
                    w = y+h;
                    
                    stepon = find(newfgtdata(:,1)>thresh,1); %Find first non-zero instance
                    newfgtdata2 = newfgtdata(stepon:end,:); %Clip data to start at first non-zero instance
                    stepoff2 = find(newfgtdata2(:,1)<thresh,1); %Find first zero at end of step
                    stepoff = stepoff2+1;
                    stepoff3 = stepoff2-1;
                    stepoff4 = stepoff2-2; %These 2 lines are technical to get right dimensions
                    stepoffNEXT = stepoff2+30;
                    step = newfgtdata2(1:stepoff,:); %creates array containing ONLY the step
                    
                    % Now clip the ASC data in teh same way
                    newascdata2 = newascdata(stepon:end,:);
                    stepasc = newascdata2(1:stepoff,:);
                    
                    % Make sure the step is more than 20 frames otherwise it is an artifact
                    r = size(step,1);
                    
                    % Interpolate Force data
                    step(r,:) = 0;
                    %%%%%
                    x0 = 0:(stepoff2); %Determine size of step array
                    x1 = 0:r/99:r; %100 data points
                    steparrayforce(:,w) = interp1(x0,step(:,1),x1,'cubic'); %interpolate force then place in array
                    
                    newfgtdata = newfgtdata2(stepoffNEXT:end,:); %Clips dataset to start AFTER previous step
                    
                    
                    % Interpolate and convert ASC data to Force
                    for l = 1:99
                        ascholdingarray(:,l) = interp1(x0,stepasc(:,l),x1, 'cubic');
                        ascholdingarray2(:,l) = ascholdingarray(:,l).*insolearea_cells(l).*1000;
                    end
                    
                    pressurematrix(:,:,w) = ascholdingarray2;
                    
                    % ASC processing to generate COP
                    pedar_fsum = squeeze(sum(ascholdingarray2,2));
                    force_matrix = ascholdingarray2; %Force matrix for whole insole
                    cop_x = sum((cent_matx.*force_matrix),2)./pedar_fsum./1000; %x-coord COP for whole insole (m)
                    cop_y = sum((cent_maty.*force_matrix),2)./pedar_fsum./1000; %y-coord COP for whole insole (m)
                    steparraycopxcalc(:,w) = cop_x(:); % place COP trace in array
                    steparraycopycalc(:,w) = cop_y(:); % place COP trace in array
                    
                    newascdata = newascdata2(stepoffNEXT:end,:); %Clips dataset to start AFTER previous step
                    
                end
                
                y = x;
                
            end %End of 2 Running trials Condition Loop
            
            % make all NaN values = 0
            steparraycopxcalc(isnan(steparraycopxcalc))=0;
            steparraycopycalc(isnan(steparraycopycalc))=0;
            % multiply by 1000 to get mm
            steparraycopxcalc = steparraycopxcalc*1000;
            steparraycopycalc = steparraycopycalc*1000;
            
            % Determine size of dataset and create empty arrays
            [aa,bb] = size(steparraycopxcalc);
            tempnormX = zeros(aa,bb+1);
            tempnormY = zeros(aa,bb+1);
            
            % Go through every cell to divide by length and width of insole
            for cc = 1:aa
                for dd = 1:bb
                    tempnormX(cc,dd) = steparraycopxcalc(cc,dd)./insole_width;
                    tempnormY(cc,dd) = steparraycopycalc(cc,dd)./insole_length;
                end
            end
            
            % Multiply all values by 100 to get to percentage
            tempnormX = tempnormX.*100;
            tempnormY = tempnormY.*100;
            
            
            clear steparrayforcetemp steparraycopxtemp steparraycopytemp steparraycopxcalctemp steparraycopycalctemp;
            
            % Now average all traces to get avg trace per condition
            steparrayforce(:,y+1) = mean(steparrayforce,2);
            steparraycopxcalc(:,y+1) = mean(steparraycopxcalc,2);
            steparraycopycalc(:,y+1) = mean(steparraycopycalc,2);
            tempnormX(:,bb+1) = mean(tempnormX,2);
            tempnormY(:,bb+1) = mean(tempnormY,2);
            indipressure(:,:) = mean(pressurematrix,3);
            
            xlswrite([savefolder subjnum condnum '.xlsx'],steparrayforce,'Force')
            xlswrite([savefolder subjnum condnum '.xlsx'],steparraycopxcalc(5:85,:),'COPx')
            xlswrite([savefolder subjnum condnum '.xlsx'],steparraycopycalc(5:85,:),'COPy')
            xlswrite([savefolder subjnum condnum '.xlsx'],tempnormX(5:85,:),'COPx_norm')
            xlswrite([savefolder subjnum condnum '.xlsx'],tempnormY(5:85,:),'COPy_norm')
            xlswrite([savefolder subjnum condnum '.xlsx'],indipressure,'Pressure')
            
            save([subjnum condnum], 'pressurematrix')
            
            disp([subjnum condnum ' is complete'])
            
            clear steparrayforce steparraycopx steparraycopy steparraycopxcalc steparraycopycalc
            
            
        else % this is for the rest of the trials which contain 3 trials per condition
            
            %%%% Trial loop %%%%
            
            y = 0;
            
            % 2 running trials per condition
            for t = 1:3
                
                % files are labelled 2 and 3
                tt = sprintf('_%d',(t+1));
                
                % Load Files
                filedirectory = [directory rawdatafolder subjnum '\' errorfoldername '\' condfolder tt];
                
                tempcondfold = [condfolder tt];
                
                % Reading in each .fgt pedar file (force data)
                pedar_fgt_filename = [filedirectory '.fgt'];
                % Reading pedar force data for the whole insole at each timestep from fgt file
                pedar_fgt_dataraw = dlmread([pedar_fgt_filename],'\t',9,5);
                pedar_fgt_data = pedar_fgt_dataraw(:,1);
                
                for i = 1:size(pedar_fgt_data)
                    if pedar_fgt_data(i) < 25
                        pedar_fgt_data(i) = 0;
                    end
                end
                
                plot(pedar_fgt_data(:,1))
                disp('Click onset and offset (x2)')
                f = ginput(2);
                newfgtdata = pedar_fgt_data(f(1,1):f(2,1),:); %clip the fgt data to boundaries set by cursor
                newfgtdataoriginal = newfgtdata;
                
                clear pedar_fgt_dataf pedar_fgt_datax pedar_fgt_datay;
                
                
                % Reading in each .asc pedar file (force data)
                pedar_asc_filename = [filedirectory '.asc'];
                % Reading pedar force data for the whole insole at each timestep from fgt file
                pedar_asc_data = dlmread([pedar_asc_filename],'\t',10,100); %'/t' tab delimited, starting at row 10, column 100
                if size(pedar_asc_data,2) == 100
                    pedar_asc_data(:,100) = [];
                end
                
                % Clip .asc data to window determined by force data
                newascdata = pedar_asc_data(f(1,1):f(2,1),:);
                newascdataoriginal = newascdata;
                
                % Issue with Pedar equipment, activates to 20-30 Hz
                % Resetting all pressure cells measuring less than 30 kPa to 0
                for i = 1:size(newascdata,1)
                    for j = 1:size(newascdata,2)
                        if newascdata(i,j) < 30
                            newascdata(i,j) = 0;
                        end
                    end
                end
                
                [peaks,locs] = findpeaks(newfgtdata,'MinPeakDistance',100,'MinPeakHeight',threshold);
                
                x = size(peaks,1);
                
                
                for h = 1:x
                    w = y+h;
                    
                    stepon = find(newfgtdata(:,1)>thresh,1); %Find first non-zero instance
                    newfgtdata2 = newfgtdata(stepon:end,:); %Clip data to start at first non-zero instance
                    stepoff2 = find(newfgtdata2(:,1)<thresh,1); %Find first zero at end of step
                    stepoff = stepoff2+1;
                    stepoff3 = stepoff2-1;
                    stepoff4 = stepoff2-2; %These 2 lines are technical to get right dimensions
                    stepoffNEXT = stepoff2+30;
                    step = newfgtdata2(1:stepoff,:); %creates array containing ONLY the step
                    
                    % Now clip the ASC data in teh same way
                    newascdata2 = newascdata(stepon:end,:);
                    stepasc = newascdata2(1:stepoff,:);
                    
                    % Make sure the step is more than 20 frames otherwise it is an artifact
                    r = size(step,1);
                    
                    % Interpolate Force data
                    step(r,:) = 0;
                    %%%%%
                    x0 = 0:(stepoff2); %Determine size of step array
                    x1 = 0:r/99:r; %100 data points
                    % for l = 1:99
                    steparrayforce(:,w) = interp1(x0,step(:,1),x1,'cubic'); %interpolate force then place in array
                    %                     end
                    
                    newfgtdata = newfgtdata2(stepoffNEXT:end,:); %Clips dataset to start AFTER previous step
                    
                    
                    % Interpolate and convert ASC data to Force
                    for l = 1:99
                        ascholdingarray(:,l) = interp1(x0,stepasc(:,l),x1, 'cubic');
                        ascholdingarray2(:,l) = ascholdingarray(:,l).*insolearea_cells(l).*1000;
                    end
                    
                    pressurematrix(:,:,w) = ascholdingarray2;
                    
                    % ASC processing to generate COP
                    pedar_fsum = squeeze(sum(ascholdingarray2,2));
                    force_matrix = ascholdingarray2; %Force matrix for whole insole
                    cop_x = sum((cent_matx.*force_matrix),2)./pedar_fsum./1000; %x-coord COP for whole insole (m)
                    cop_y = sum((cent_maty.*force_matrix),2)./pedar_fsum./1000; %y-coord COP for whole insole (m)
                    steparraycopxcalc(:,w) = cop_x(:); % place COP trace in array
                    steparraycopycalc(:,w) = cop_y(:); % palce COP trace in array
                    
                    newascdata = newascdata2(stepoffNEXT:end,:); %Clips dataset to start AFTER previous step
                    
                end
                
                y = size(steparrayforce,2);
                
            end %End of 2 Running trials Condition Loop
            
            % make all NaN values = 0
            steparraycopxcalc(isnan(steparraycopxcalc))=0;
            steparraycopycalc(isnan(steparraycopycalc))=0;
            % multiply by 1000 to get mm
            steparraycopxcalc = steparraycopxcalc*1000;
            steparraycopycalc = steparraycopycalc*1000;
            
            % Determine size of dataset and create empty arrays
            [aa,bb] = size(steparraycopxcalc);
            tempnormX = zeros(aa,bb+1);
            tempnormY = zeros(aa,bb+1);
            
            % Go through every cell to divide by length and width of insole
            for cc = 1:aa
                for dd = 1:bb
                    tempnormX(cc,dd) = steparraycopxcalc(cc,dd)./insole_width;
                    tempnormY(cc,dd) = steparraycopycalc(cc,dd)./insole_length;
                end
            end
            
            % Multiply all values by 100 to get to percentage
            tempnormX = tempnormX.*100;
            tempnormY = tempnormY.*100;
            
            
            clear steparrayforcetemp steparraycopxtemp steparraycopytemp steparraycopxcalctemp steparraycopycalctemp;
            
            % Now average all traces to get avg trace per condition
            steparrayforce(:,y+1) = mean(steparrayforce,2);
            steparraycopxcalc(:,y+1) = mean(steparraycopxcalc,2);
            steparraycopycalc(:,y+1) = mean(steparraycopycalc,2);
            tempnormX(:,bb+1) = mean(tempnormX,2);
            tempnormY(:,bb+1) = mean(tempnormY,2);
            indipressure(:,:) = mean(pressurematrix,3);
            
            xlswrite([savefolder subjnum condnum '.xlsx'],steparrayforce,'Force')
            xlswrite([savefolder subjnum condnum '.xlsx'],steparraycopxcalc(5:85,:),'COPx')
            xlswrite([savefolder subjnum condnum '.xlsx'],steparraycopycalc(5:85,:),'COPy')
            xlswrite([savefolder subjnum condnum '.xlsx'],tempnormX(5:85,:),'COPx_norm')
            xlswrite([savefolder subjnum condnum '.xlsx'],tempnormY(5:85,:),'COPy_norm')
            xlswrite([savefolder subjnum condnum '.xlsx'],indipressure,'Pressure')
            
            save([subjnum condnum], 'pressurematrix')
            
            disp([subjnum condnum ' is complete'])
            
            clear steparrayforce steparraycopx steparraycopy steparraycopxcalc steparraycopycalc
            
        end %End of Subject Loop
        
    end
end



%% Process .mat pressure files for SVM (INCOMPLETE)

svmfilepref = zeros(1,160); %create an array to hold data for SVM anaylsis (pref)
svmfileleast = zeros(1,160); %create an array to hold data for SVM anaylsis (least)
y = 0; %start ticker to keep track of array location (pref)
b = 0; %start ticker to keep track of array location (least)


%%% Subject Loop %%%

for s = 1:40 %40 subjects
    
    % For naming convention
    if s < 10
        subjnum = sprintf('Runner00%d',(s));
    else
        subjnum = sprintf('Runner0%d',(s));
    end
    
    %     % Load testing sheet to determine correct insole size from data collection sheet
    %     CKL_Datasheet = ['\\CINCO\students\BennoGroup\CLam\Running Style_CKL\Lam_ComfortRankings.xlsx'];
    %
    %     order_rank_VAS_master = xlsread(CKL_Datasheet,'comfort and fit rankings','C3:G122');
    %
    %     % Find exact cell to determine shoe order
    %     bbb = (3*(s-1))+1;
    %     ccc = (3*(s-1))+2;
    %     ddd = (3*(s-1))+3;
    
    
    %%% Condition Loop %%%
    
    for u = 1:5 %5 conditions
        condnum = sprintf('_c%d',(u));
        filename = [subjnum condnum];
        
        % PREFERRED
        if strcmpi(filename,'Runner004_c5')|| strcmpi(filename,'Runner011_c5')|| strcmpi(filename,'Runner012_c5')||...
                strcmpi(filename,'Runner013_c5')|| strcmpi(filename,'Runner014_c1')|| strcmpi(filename,'Runner022_c5')||...
                strcmpi(filename,'Runner024_c4')|| strcmpi(filename,'Runner027_c4')|| strcmpi(filename,'Runner028_c2')||...
                strcmpi(filename,'Runner034_c4')|| strcmpi(filename,'Runner035_c4')|| strcmpi(filename,'Runner037_c5')||...
                strcmpi(filename,'Runner009_c5')|| strcmpi(filename,'Runner016_c5')
            
            
            
            
            % Load raw pressure matrix
            file = [rawpressurefolder '\' filename '.mat'];
            load(file); %NOTE: File can be accessed through matrix name: 'pressurematrix'
            [p,q,r] = size(pressurematrix);
            
            if r<30
                display('Error, less than 30 steps')
                
            else
                for h = 1:30
                    z = y+h;
                    svmfilepref(z,:) = pressurematrix(5:85,:).'; %??? INCOMPLETE
                    labelpref(z,1) = prefcount;
                    labelpref(z,2) = 1;
                    
                    
                end
                y = size(svmfilepref,1);
                
            end
            
            
            
            
            
            
            % LEAST PREFERRED
        elseif strcmpi(filename,'Runner003_c3')|| strcmpi(filename,'Runner005_c4')|| strcmpi(filename,'Runner017_c1')||...
                strcmpi(filename,'Runner023_c1')|| strcmpi(filename,'Runner025_c5')|| strcmpi(filename,'Runner008_c3')||...
                strcmpi(filename,'Runner015_c4')|| strcmpi(filename,'Runner018_c5')|| strcmpi(filename,'Runner020_c4')||...
                strcmpi(filename,'Runner021_c2')|| strcmpi(filename,'Runner032_c4')|| strcmpi(filename,'Runner038_c3')
            
            % Load raw pressure matrix
            file = [rawpressurefolder '\' filename '.mat'];
            
            
            
        end
        
        
        
    end
    
end


%% Normalize Force to Bodyweight - extract BW from quiet stance trial


% Start Fresh
close all;clc;clear all;

% Allocate directory
directory  = '\\CINCO\students\BennoGroup\MMohr\Running Style';
quietstancedirectory = '\\CINCO\students\BennoGroup\MMohr\Running Style\RawPedarforProcessing\';
datafiles  = '\\CINCO\students\BennoGroup\MMohr\Running Style\Pedar\ChrisLamProcessing\Complete\';


for s = 1:40 %40 subjects
    
    % For naming convention
    if s < 10
        subjnum = sprintf('Runner00%d',(s));
        subjnum1 = sprintf('r00%d',(s));
    else
        subjnum = sprintf('Runner0%d',(s));
        subjnum1 = sprintf('r0%d',(s));
    end
    
%     % Load testing sheet to determine correct insole size from data collection sheet
%     subjquestionnaire = [directory '\Data Collection\20160815_Running Style Questionnaire Results.xlsx'];
%     
%     % Find exact cell to determine the weight
%     aaa = (2*s)+1;
%     subj_weight = sprintf('AZ%d',(aaa));
%     sheet = 'Sujective Questionaire Results';
%     [var,weight] = xlsread(subjquestionnaire,sheet,subj_weight);
%     BW = var*9.81;
        
    QS_filename = [quietstancedirectory subjnum '\' subjnum1 '_qs.fgt'];
    QS_file = dlmread(QS_filename,'\t',9,5);
    QS_data = QS_file(:,1);
    
    plot(QS_data) %plot to pick out onsets and offsets of quiet stance
    disp('Click onset and offset')
    ff = ginput(2);
    weightdata = QS_data(ff(1,1):ff(2,1)); %clip the fgt data to boundaries set by cursor
    
    BW = mean(weightdata); %take an average, do not divide by accel as we want it in N
    
    
    % Normalize Force to bodyweight
    for u = 1:5
        condnum = sprintf('_c%d',(u));
        filename = [subjnum condnum];
        file = [datafiles subjnum '\' filename '.xlsx'];
        Force_sheet = 'Force';
        
        % Load file to extract COPx and COPy
        temparrayforce = xlsread(file, Force_sheet);
        
        % Determine size of dataset and create empty arrays
        [a,b] = size(temparrayforce);
        tempnormBW = zeros(a,b);
        
        % Go through every cell to divide by bodyweight
        for c = 1:a
            for d = 1:b
                tempnormBW(c,d) = temparrayforce(c,d)./BW;
            end
        end
        
        % Write page to excel file
        xlswrite(file,tempnormBW,'Force_Norm')
        
    end
    disp([subjnum ' is complete'])
    
end

close


%% Process Excel data to analyze COP velocity, and displacement and peak force


% Array Allocation

% Arrays for saving Excel Summary files
Force_results_array = zeros(1,25);
Force_trace_array = zeros(100,1);
COPx_results_array = zeros(1,25);
X_Slope_array = zeros(1,25);
COPx_trace_array = zeros(81,1);
COPy_results_array = zeros(1,25);
Y_Slope_array = zeros(1,25);
COPy_trace_array = zeros(81,1);


%%% Subject Loop %%%

for s = 1:40 %40 subjects
    
    % For naming convention
    if s < 10
        subjnum = sprintf('Runner00%d',(s));
    else
        subjnum = sprintf('Runner0%d',(s));
    end
    
    
    %%% Condition Loop %%%
    
    for u = 1:5 %5 conditions
        condnum = sprintf('_c%d',(u));
        filename = [subjnum condnum];
        file = [processedfolder subjnum '\' filename '.xlsx'];
        % Pages within each file
        Force_sheet = 'Force_Norm';
        COPx_sheet = 'COPx_norm';
        COPy_sheet = 'COPy_norm';
        
        % Load Files
        forcearray = xlsread(file,Force_sheet);
        copxarray = xlsread(file,COPx_sheet);
        copyarray = xlsread(file,COPy_sheet);
        
        numsteps = size(forcearray,2);
        
        % Create matrix to hold all processed data per subject/condition
        raw_results_matrix = zeros(30,5,5);
        
        
        %%% Force Processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [force_max,force_max_index] = max(forcearray,[],1);
        forceslope_12 = forcearray(12,:);
        forceslope_45 = forcearray(45,:);
        forceslope_75 = forcearray(75,:);
        Force_results = [forceslope_12' forceslope_45' forceslope_75' force_max' force_max_index'];
        Force_results = Force_results(1:30,:);
        forcemean = mean(Force_results,1);
        Force_results = cat(1,Force_results,forcemean);
        
        %Want to get all 30 traces, but only use 30, take an avg and stdev
        Force_trace = forcearray(:,1:30);
        forcetracemean = mean(Force_trace,2);
        forcetracestdev = std(Force_trace');
        Force_trace = cat(2,Force_trace,forcetracemean,forcetracestdev');
        
        
        %%% COPx Processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        COPx_max = max(copxarray,[],1);
        COPx_min = min(copxarray,[],1);
        COPx_ranges = COPx_max-COPx_min;
        COPx_results = [COPx_max' COPx_min' COPx_ranges'];
        COPx_results = COPx_results(1:30,:);
        xmean2 = mean(COPx_results,1);
        COPx_results = cat(1,COPx_results,xmean2);

        gg = diff(copxarray);
        [xslope_max,x_max_index] = max(gg,[],1);
        x_max_index = x_max_index+4;
        xslope_12 = gg(8,:);
        xslope_45 = gg(41,:);
        xslope_75 = gg(71,:);
        X_Slope = [xslope_12' xslope_45' xslope_75' xslope_max' x_max_index'];
        X_Slope = X_Slope(1:30,:);
        xmean = mean(X_Slope,1);
        X_Slope = cat(1,X_Slope,xmean);
        
        %Want to get all 30 traces, but only use 30, take an avg and stdev
        COPx_trace = copxarray(:,1:30);
        copxtracemean = mean(COPx_trace,2);
        copxtracestdev = std(COPx_trace');
        COPx_trace = cat(2,COPx_trace,copxtracemean,copxtracestdev');
        
        
        %%% COPy Processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        COPy_max = max(copyarray,[],1);
        COPy_min = min(copyarray,[],1);
        COPy_ranges = COPy_max-COPy_min;
        COPy_results = [COPy_max' COPy_min' COPy_ranges'];
        COPy_results = COPy_results(1:30,:);
        ymean2 = mean(COPy_results,1);
        COPy_results = cat(1,COPy_results,ymean2);
        
        hh = diff(copyarray); %differentiate the normalized displacement to get slope
        [yslope_max,y_max_index] = max(hh,[],1);
        y_max_index = y_max_index+4; %since the window of COP is 5:85
        yslope_12 = gg(8,:); %corrected for window of 5:85
        yslope_45 = gg(41,:);
        yslope_75 = gg(71,:);
        Y_Slope = [yslope_12' yslope_45' yslope_75' yslope_max' y_max_index'];
        Y_Slope = Y_Slope(1:30,:);
        ymean = mean(Y_Slope,1);
        Y_Slope = cat(1,Y_Slope,ymean);
        
        %Want to get all 30 traces, but only use 30, take an avg and stdev
        COPy_trace = copyarray(:,1:30);
        copytracemean = mean(COPy_trace,2);
        copytracestdev = std(COPy_trace');
        COPy_trace = cat(2,COPy_trace,copytracemean,copytracestdev');
        
        
        %%% FILE SAVING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Setup onsets and offsets for array
        on5 = ((u-1)*5)+1;
        off5 = on5+4;
        on3 = ((u-1)*3)+1;
        off3 = on3+2;
        on2 = ((u-1)*2)+1;
        off2 = on2+1;
        
        
        % Array compiling data for all subjects and shoe conditions
        % Individual values per shoe condition per subject
        Force_results_array(s,on5:off5) = Force_results(31,:);
        COPx_results_array(s,on3:off3) = COPx_results(31,:);
        X_Slope_array(s,on5:off5) = X_Slope(31,:);
        COPy_results_array(s,on3:off3) = COPy_results(31,:);
        Y_Slope_array(s,on5:off5) = Y_Slope(31,:);
        
        % Whole traces per shoe condition per subject
        Force_trace_array(:,on2:off2) = Force_trace(:,31:32);
        COPx_trace_array(:,on2:off2) = COPx_trace(:,31:32);
        COPy_trace_array(:,on2:off2) = COPy_trace(:,31:32);
        
        
        % Matrices compiling data fro all raw results to be saved as .mat
        % Compile then save all raw calculated results
        raw_results_matrix(:,:,1) = Force_results(1:30,:);
        raw_results_matrix(:,1:3,2) = COPx_results(1:30,:);
        raw_results_matrix(:,:,3) = X_Slope(1:30,:);
        raw_results_matrix(:,1:3,4) = COPy_results(1:30,:);
        raw_results_matrix(:,:,5) = Y_Slope(1:30,:);
        
        rawname = ' Raw_Results';
        save([subjnum condnum rawname], 'raw_results_matrix')
        
    end
    
    % Save Trace files for each subject
    xlswrite([processedfolder subjnum '\' subjnum 'processed_traces.xlsx'],Force_trace_array,'Force_traces')
    xlswrite([processedfolder subjnum '\' subjnum 'processed_traces.xlsx'],COPx_trace_array,'COPx_traces')
    xlswrite([processedfolder subjnum '\' subjnum 'processed_traces.xlsx'],COPy_trace_array,'COPy_traces')
    
end

% Save Results Arrays with all subjects
xlswrite([savefolder 'Processed Results1.xlsx'],Force_results_array,'Force')
xlswrite([savefolder 'Processed Results1.xlsx'],COPx_results_array,'COPx_disp')
xlswrite([savefolder 'Processed Results1.xlsx'],X_Slope_array,'COPx_slope')
xlswrite([savefolder 'Processed Results1.xlsx'],COPy_results_array,'COPy_disp')
xlswrite([savefolder 'Processed Results1.xlsx'],Y_Slope_array,'COPy_slope')




%% Reorder results to match shoe rather than condition


% Start Fresh
% close all;clc;clear all;

% Allocate directory
groupingfile = '\\CINCO\students\BennoGroup\MMohr\Running Style\Data Collection\Lam_Comfort_Rankings_Grouping.xlsx';
processeddatafile  = '\\CINCO\students\BennoGroup\MMohr\Running Style\Pedar\ChrisLamProcessing\cop_disp_slope.xlsx';

% Pages within each file
CondOrder_sheet = 'Subject Trial Order';
Processed_data_sheet = 'Sheet2OUTOFORDER';

% Array Allocation
% Shoe condition order
condorder = xlsread(groupingfile,CondOrder_sheet,'B2:F41');
% Generate arrays for processed data
processed_force = xlsread(processeddatafile,Processed_data_sheet,'B3:Z42');
processed_X_disp = xlsread(processeddatafile,Processed_data_sheet,'AC3:AQ42');
processed_X_slope = xlsread(processeddatafile,Processed_data_sheet,'AT3:BR42');
processed_Y_disp = xlsread(processeddatafile,Processed_data_sheet,'BU3:CI42');
processed_Y_slope = xlsread(processeddatafile,Processed_data_sheet,'CL3:DJ42');
processed_force2 = zeros(40,25);
processed_X_disp2 = zeros(40,15);
processed_X_slope2 = zeros(40,25);
processed_Y_disp2 = zeros(40,15);
processed_Y_slope2 = zeros(40,25);


for s = 1:40 %40 subjects    
    for c = 1:5
        
        a = condorder(s,c);
        a5 = ((a-1)*5)+1;
        a3 = ((a-1)*3)+1;
        c5 = ((c-1)*5)+1;
        c3 = ((c-1)*3)+1;
        
        processed_force2(s,c5:c5+4) = processed_force(s,a5:a5+4);
        processed_X_disp2(s,c3:c3+2) = processed_X_disp(s,a3:a3+2);
        processed_X_slope2(s,c5:c5+4) = processed_X_slope(s,a5:a5+4);
        processed_Y_disp2(s,c3:c3+2) = processed_Y_disp(s,a3:a3+2);
        processed_Y_slope2(s,c5:c5+4) = processed_Y_slope(s,a5:a5+4);
        
    end
end
    
xlswrite([savefolder 'Processed Results_ord.xlsx'],processed_force2,'Force')
xlswrite([savefolder 'Processed Results_ord.xlsx'],processed_X_disp2,'COPx_disp')
xlswrite([savefolder 'Processed Results_ord.xlsx'],processed_X_slope2,'COPx_slope')
xlswrite([savefolder 'Processed Results_ord.xlsx'],processed_Y_disp2,'COPy_disp')
xlswrite([savefolder 'Processed Results_ord.xlsx'],processed_Y_slope2,'COPy_slope')



%% Analyzing based on Grouping
    
    
% Start Fresh
close all;clc;clear all;

% Allocate directory
directory  = '\\CINCO\students\BennoGroup\MMohr\Running Style\Pedar\ChrisLamProcessing\';
groupingfile = '\\CINCO\students\BennoGroup\MMohr\Running Style\Data Collection\Lam_Comfort_Rankings_Grouping.xlsx';
processeddatafile  = '\\CINCO\students\BennoGroup\MMohr\Running Style\Pedar\ChrisLamProcessing\cop_disp_slope.xlsx';

% Pages within each file
% CondOrder_sheet = 'Subject Trial Order';
Grouping_sheet = 'Grouping2';
Processed_data_sheet = 'ResultsREDUCED';


% Load grouping lists

% Grouping within subject
group_subject = xlsread(groupingfile,Grouping_sheet,'AN6:AO45');
% Grouping based on VAS
group_VAS_ALL = xlsread(groupingfile,Grouping_sheet,'A7:D56');
group_VAS_G  = xlsread(groupingfile,Grouping_sheet,'F8:G17');
group_VAS_Z  = xlsread(groupingfile,Grouping_sheet,'I8:J17');
group_VAS_B  = xlsread(groupingfile,Grouping_sheet,'F25:G34');
group_VAS_S  = xlsread(groupingfile,Grouping_sheet,'I25:J34');
group_VAS_P  = xlsread(groupingfile,Grouping_sheet,'F42:G51');
% Grouping based on Comfort Ranking
group_rank = xlsread(groupingfile,Grouping_sheet,'N7:W28');
% Grouping based on 'yes/no'
group_YN_avg = xlsread(groupingfile,Grouping_sheet,'AA7:AD23');
group_YN_testday = xlsread(groupingfile,Grouping_sheet,'AF7:AI39');

% Generate arrays for processed data
processed_X = xlsread(processeddatafile,Processed_data_sheet,'B3:P42');
processed_Y = xlsread(processeddatafile,Processed_data_sheet,'S3:AG42');

%%%%%%%% Gender

% Load Subject list of gender
gender_sheet = 'Groupingshoealphabet';
genderlist = xlsread(groupingfile,gender_sheet,'AA85:AB124');
malelist = xlsread(groupingfile,gender_sheet,'AD85:AE104');
femalelist = xlsread(groupingfile,gender_sheet,'AF85:AG104');

% Within subject grouping
withinsubjf = xlsread(groupingfile,Grouping_sheet,'AL50:AN69');
withinsubjm = xlsread(groupingfile,Grouping_sheet,'AP50:AR69');
% VAS grouping
% ALL
VAS_All_f = xlsread(groupingfile,Grouping_sheet,'A68:D92');
VAS_All_m = xlsread(groupingfile,Grouping_sheet,'A98:D122');
% By Shoe
VAS_G_f = xlsread(groupingfile,Grouping_sheet,'A129:B138');
VAS_Z_f = xlsread(groupingfile,Grouping_sheet,'D129:E138');
VAS_B_f = xlsread(groupingfile,Grouping_sheet,'G129:H138');
VAS_S_f = xlsread(groupingfile,Grouping_sheet,'J129:K138');
VAS_P_f = xlsread(groupingfile,Grouping_sheet,'M129:N138');
VAS_G_m = xlsread(groupingfile,Grouping_sheet,'A146:B155');
VAS_Z_m = xlsread(groupingfile,Grouping_sheet,'D146:E155');
VAS_B_m = xlsread(groupingfile,Grouping_sheet,'G146:H155');
VAS_S_m = xlsread(groupingfile,Grouping_sheet,'J146:K155');
VAS_P_m = xlsread(groupingfile,Grouping_sheet,'M146:N155');
% Comfort Rank Grouping
Top2_rank_f = zeros(15,10);
Top2_rank_f_temp = xlsread(groupingfile,Grouping_sheet,'N41:W51');
Top2_rank_m = xlsread(groupingfile,Grouping_sheet,'N57:W71');
Top2_rank_f(1:11,:) = Top2_rank_f_temp;
%%%%%%%%

% Matrix for holding all COP data
COPmatrix = zeros(6,5,40);


% Create matrix of data
for s = 1:40 %40 subjects
    for c = 1:5
        cc = ((c-1)*3)+1;
        COPmatrix(1:3,c,s) = processed_X(s,cc:cc+2);
        COPmatrix(4:6,c,s) = processed_Y(s,cc:cc+2);
    end
end
    
save('COPmatrix')

clear s c cc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Within subject grouping

group_subject_array = zeros(40,12);

for a = 1:40
    
    b = group_subject(a,1);
    bb = group_subject(a,2);
    group_subject_array(a,1:6) = COPmatrix(:,b,a);
    group_subject_array(a,7:12) = COPmatrix(:,bb,a);
    
end

clear a b bb
% group_subject_array is final array

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% VAS grouping

% ALL

group_VAS_ALL_array = zeros(50,12);

for d = 1:50
    
    e = group_VAS_ALL(d,1); %subject (pref)
    f = group_VAS_ALL(d,2); %shoe (pref)
    ee = group_VAS_ALL(d,3); %subject (nonpref)
    ff = group_VAS_ALL(d,4); %shoe (nonpref)
    group_VAS_ALL_array(d,1:6) = COPmatrix(:,f,e);
    group_VAS_ALL_array(d,7:12) = COPmatrix(:,ff,ee);
    
end

% By Shoe

group_VAS_G_array = zeros(10,12);
group_VAS_Z_array = zeros(10,12);
group_VAS_B_array = zeros(10,12);
group_VAS_S_array = zeros(10,12);
group_VAS_P_array = zeros(10,12);

for g = 1:10
    
   h = group_VAS_G(g,1);
   hh = group_VAS_G(g,2);
   i = group_VAS_Z(g,1);
   ii = group_VAS_Z(g,2);
   j = group_VAS_B(g,1);
   jj = group_VAS_B(g,2);
   k = group_VAS_S(g,1);
   kk = group_VAS_S(g,2);
   l = group_VAS_P(g,1);
   ll = group_VAS_P(g,2);
   group_VAS_G_array(g,1:6) = COPmatrix(:,1,h);
   group_VAS_G_array(g,7:12) = COPmatrix(:,1,hh);
   group_VAS_Z_array(g,1:6) = COPmatrix(:,2,i);
   group_VAS_Z_array(g,7:12) = COPmatrix(:,2,ii);
   group_VAS_B_array(g,1:6) = COPmatrix(:,3,j);
   group_VAS_B_array(g,7:12) = COPmatrix(:,3,jj);
   group_VAS_S_array(g,1:6) = COPmatrix(:,4,k);
   group_VAS_S_array(g,7:12) = COPmatrix(:,4,kk);
   group_VAS_P_array(g,1:6) = COPmatrix(:,5,l);
   group_VAS_P_array(g,7:12) = COPmatrix(:,5,ll);
    
end

group_VAS_indiv_array = cat(1,group_VAS_G_array,group_VAS_Z_array,group_VAS_B_array,group_VAS_S_array,group_VAS_P_array);

clear d e ee f ff g h hh i ii j jj k kk l ll group_VAS_G_array group_VAS_Z_array group_VAS_B_array group_VAS_S_array group_VAS_P_array
% group_VAS_ALL_array and group_VAS_indiv_array are the final arrays

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Comfort Rank Grouping

group_rank_G_array = zeros(22,12);
group_rank_Z_array = zeros(17,12);
group_rank_B_array = zeros(21,12);
group_rank_S_array = zeros(16,12);
group_rank_P_array = zeros(22,12);

group_rank(isnan(group_rank)) = 0 ;

for t = 1:22
    
    v = group_rank(t,1);
    vv = group_rank(t,6);
    w = group_rank(t,2);
    ww = group_rank(t,7);
    x = group_rank(t,3);
    xx = group_rank(t,8);
    y = group_rank(t,4);
    yy = group_rank(t,9);
    z = group_rank(t,5);
    zz = group_rank(t,10);
    
    if v == 0
    else
        group_rank_G_array(t,1:6) = COPmatrix(:,1,v);
    end
    if vv == 0
    else
        group_rank_G_array(t,7:12) = COPmatrix(:,1,vv);
    end

        if w == 0
    else
        group_rank_Z_array(t,1:6) = COPmatrix(:,2,w);
    end
    if ww == 0
    else
        group_rank_Z_array(t,7:12) = COPmatrix(:,2,ww);
    end

        if x == 0
    else
        group_rank_B_array(t,1:6) = COPmatrix(:,3,x);
    end
    if xx == 0
    else
        group_rank_B_array(t,7:12) = COPmatrix(:,3,xx);
    end

        if y == 0
    else
        group_rank_S_array(t,1:6) = COPmatrix(:,4,y);
    end
    if yy == 0
    else
        group_rank_S_array(t,7:12) = COPmatrix(:,4,yy);
    end

    if z == 0
    else
        group_rank_P_array(t,1:6) = COPmatrix(:,5,z);
    end
    if zz == 0
    else
        group_rank_P_array(t,7:12) = COPmatrix(:,5,zz);
    end

end

group_rank_array = cat(1,group_rank_G_array,group_rank_Z_array,group_rank_B_array,group_rank_S_array,group_rank_P_array);

clear t v vv w ww x xx y yy z zz group_rank_G_array group_rank_Z_array group_rank_B_array group_rank_S_array group_rank_P_array
% group_rank_array is final array

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Yes/No Grouping

group_YN_avg_array = zeros(17,12);
group_YN_testday_array = zeros(33,12);

for m = 1:17
    
    n = group_YN_avg(m,1); %subject pref
    nn = group_YN_avg(m,2); %shoe pref
    o = group_YN_avg(m,3); %subject no
    oo = group_YN_avg(m,4); %shoe no
    group_YN_avg_array(m,1:6) = COPmatrix(:,nn,n);
    group_YN_avg_array(m,7:12) = COPmatrix(:,oo,o);
    
end

for p = 1:33
    
    q = group_YN_testday(p,1); %subject pref
    qq = group_YN_testday(p,2); %shoe pref
    r = group_YN_testday(p,3); %subject no
    rr = group_YN_testday(p,4); %shoe no
    group_YN_testday_array(p,1:6) = COPmatrix(:,qq,q);
    group_YN_testday_array(p,7:12) = COPmatrix(:,rr,r);
    
end

clear m n nn o oo p q qq r rr
% group_YN_avg_array and group_YN_testday_array are the final arrays

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Between Shoes 

var_Xdisp_array = squeeze(COPmatrix(1,:,:));
var_Xslo_array = squeeze(COPmatrix(2,:,:));
var_Xind_array = squeeze(COPmatrix(3,:,:));
var_Ydisp_array = squeeze(COPmatrix(4,:,:));
var_Yslo_array = squeeze(COPmatrix(5,:,:));
var_Yind_array = squeeze(COPmatrix(6,:,:));

var_Xdisp_array = var_Xdisp_array';
var_Xslo_array = var_Xslo_array';
var_Xind_array = var_Xind_array';
var_Ydisp_array = var_Ydisp_array';
var_Yslo_array = var_Yslo_array';
var_Yind_array = var_Yind_array';

betweenshoes_array = cat(2,var_Xdisp_array,var_Xslo_array,var_Xind_array,var_Ydisp_array,var_Yslo_array,var_Yind_array);
clear var_Xdisp_array var_Xslo_array var_Xind_array var_Ydisp_array var_Yslo_array var_Yind_array

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sex Differences

% Load Subject list of gender
% gender_sheet = 'Groupingshoealphabet';
% genderlist = xlsread(groupingfile,gender_sheet,'AA85:AB124');
% malelist = xlsread(groupingfile,gender_sheet,'AD85:AE104');
% femalelist = xlsread(groupingfile,gender_sheet,'AF85:AG104');

% Within subject grouping
withinsubjf_array = zeros(20,12);
withinsubjm_array = zeros(20,12);

for a = 1:20
    
    b = withinsubjf(a,1);
    bb = withinsubjf(a,2);
    bbb = withinsubjf(a,3);
    c = withinsubjm(a,1);
    cc = withinsubjm(a,2);
    ccc = withinsubjm(a,3);
    withinsubjf_array(a,1:6) = COPmatrix(:,bb,b);
    withinsubjf_array(a,7:12) = COPmatrix(:,bbb,b);
    withinsubjm_array(a,1:6) = COPmatrix(:,cc,c);
    withinsubjm_array(a,7:12) = COPmatrix(:,ccc,c);
    
end

withinsubj_gender = cat(2,withinsubjf_array,withinsubjm_array);
clear a b bb bbb c cc ccc withinshoef_array withinshoem_array
% withinsubj_gender is the final array

%%%%%

% VAS grouping
% ALL
VAS_All_f_array = zeros(25,12);
VAS_All_m_array = zeros(25,12);

for d = 1:25
    
    e = VAS_All_f(d,1); %subject (pref)
    f = VAS_All_f(d,2); %shoe (pref)
    ee = VAS_All_f(d,3); %subject (nonpref)
    ff = VAS_All_f(d,4); %shoe (nonpref)
    eee = VAS_All_m(d,1); %subject (pref)
    fff = VAS_All_m(d,2); %shoe (pref)
    eeee = VAS_All_m(d,3); %subject (nonpref)
    ffff = VAS_All_m(d,4); %shoe (nonpref)
    VAS_All_f_array(d,1:6) = COPmatrix(:,f,e);
    VAS_All_f_array(d,7:12) = COPmatrix(:,ff,ee);
    VAS_All_m_array(d,1:6) = COPmatrix(:,fff,eee);
    VAS_All_m_array(d,7:12) = COPmatrix(:,ffff,eeee);
    
end

VAS_All_gender = cat(2,VAS_All_f_array,VAS_All_m_array);
clear d e ee eee eeee f ff fff ffff VAS_All_f_array VAS_All_m_array
% VAS_All_gender is the final array


% By Shoe
VAS_G_f_array = zeros(10,12);
VAS_Z_f_array = zeros(10,12);
VAS_B_f_array = zeros(10,12);
VAS_S_f_array = zeros(10,12);
VAS_P_f_array = zeros(10,12);
VAS_G_m_array = zeros(10,12);
VAS_Z_m_array = zeros(10,12);
VAS_B_m_array = zeros(10,12);
VAS_S_m_array = zeros(10,12);
VAS_P_m_array = zeros(10,12);

for g = 1:10
    
   h = VAS_G_f(g,1);
   hh = VAS_G_f(g,2);
   hhh = VAS_G_m(g,1);
   hhhh = VAS_G_m(g,2);
   i = VAS_Z_f(g,1);
   ii = VAS_Z_f(g,2);
   iii = VAS_Z_m(g,1);
   iiii = VAS_Z_m(g,2);
   j = VAS_B_f(g,1);
   jj = VAS_B_f(g,2);
   jjj = VAS_B_m(g,1);
   jjjj = VAS_B_m(g,2);
   k = VAS_S_f(g,1);
   kk = VAS_S_f(g,2);
   kkk = VAS_S_m(g,1);
   kkkk = VAS_S_m(g,2);
   l = VAS_P_f(g,1);
   ll = VAS_P_f(g,2);
   lll = VAS_P_m(g,1);
   llll = VAS_P_m(g,2);
   VAS_G_f_array(g,1:6) = COPmatrix(:,1,h);
   VAS_G_f_array(g,7:12) = COPmatrix(:,1,hh);
   VAS_G_m_array(g,1:6) = COPmatrix(:,1,hhh);
   VAS_G_m_array(g,7:12) = COPmatrix(:,1,hhhh);
   VAS_Z_f_array(g,1:6) = COPmatrix(:,2,i);
   VAS_Z_f_array(g,7:12) = COPmatrix(:,2,ii);
   VAS_Z_m_array(g,1:6) = COPmatrix(:,2,iii);
   VAS_Z_m_array(g,7:12) = COPmatrix(:,2,iiii);
   VAS_B_f_array(g,1:6) = COPmatrix(:,3,j);
   VAS_B_f_array(g,7:12) = COPmatrix(:,3,jj);
   VAS_B_m_array(g,1:6) = COPmatrix(:,3,jjj);
   VAS_B_m_array(g,7:12) = COPmatrix(:,3,jjjj);
   VAS_S_f_array(g,1:6) = COPmatrix(:,4,k);
   VAS_S_f_array(g,7:12) = COPmatrix(:,4,kk);
   VAS_S_m_array(g,1:6) = COPmatrix(:,4,kkk);
   VAS_S_m_array(g,7:12) = COPmatrix(:,4,kkkk);
   VAS_P_f_array(g,1:6) = COPmatrix(:,5,l);
   VAS_P_f_array(g,7:12) = COPmatrix(:,5,ll);
   VAS_P_m_array(g,1:6) = COPmatrix(:,5,lll);
   VAS_P_m_array(g,7:12) = COPmatrix(:,5,llll);
    
end

VAS_indiv_female = cat(2,VAS_G_f_array,VAS_Z_f_array,VAS_B_f_array,VAS_S_f_array,VAS_P_f_array);
VAS_indiv_male = cat(2,VAS_G_m_array,VAS_Z_m_array,VAS_B_m_array,VAS_S_m_array,VAS_P_m_array);
VAS_indiv_gender = cat(1,VAS_indiv_female,VAS_indiv_male);
clear g h hh hhh hhhh i ii iii iiii j jj jjj jjjj k kk kkk kkkk l ll lll llll VAS_G_f_array VAS_Z_f_array VAS_B_f_array VAS_S_f_array VAS_P_f_array VAS_G_m_array VAS_Z_m_array VAS_B_m_array VAS_S_m_array VAS_P_m_array
% VAS_indiv_gender is the final array

%%%%%

% Comfort Rank Grouping
group_rank_G_f_array = zeros(15,12);
group_rank_Z_f_array = zeros(15,12);
group_rank_B_f_array = zeros(15,12);
group_rank_S_f_array = zeros(15,12);
group_rank_P_f_array = zeros(15,12);
group_rank_G_m_array = zeros(15,12);
group_rank_Z_m_array = zeros(15,12);
group_rank_B_m_array = zeros(15,12);
group_rank_S_m_array = zeros(15,12);
group_rank_P_m_array = zeros(15,12);

Top2_rank_f(isnan(Top2_rank_f)) = 0 ;
Top2_rank_m(isnan(Top2_rank_m)) = 0 ;

for t = 1:15
    
    v = Top2_rank_f(t,1);
    vv = Top2_rank_f(t,6);
    vvv = Top2_rank_m(t,1);
    vvvv = Top2_rank_m(t,6);
    w = Top2_rank_f(t,2);
    ww = Top2_rank_f(t,7);
    www = Top2_rank_m(t,2);
    wwww = Top2_rank_m(t,7);
    x = Top2_rank_f(t,3);
    xx = Top2_rank_f(t,8);
    xxx = Top2_rank_m(t,3);
    xxxx = Top2_rank_m(t,8);
    y = Top2_rank_f(t,4);
    yy = Top2_rank_f(t,9);
    yyy = Top2_rank_m(t,4);
    yyyy = Top2_rank_m(t,9);
    z = Top2_rank_f(t,5);
    zz = Top2_rank_f(t,10);
    zzz = Top2_rank_m(t,5);
    zzzz = Top2_rank_m(t,10);
    
    if v == 0
    else
        group_rank_G_f_array(t,1:6) = COPmatrix(:,1,v);
    end
    if vv == 0
    else
        group_rank_G_f_array(t,7:12) = COPmatrix(:,1,vv);
    end
    if vvv == 0
    else
        group_rank_G_m_array(t,1:6) = COPmatrix(:,1,vvv);
    end
    if vvvv == 0
    else
        group_rank_G_m_array(t,7:12) = COPmatrix(:,1,vvvv);
    end
%%%
    if w == 0
    else
        group_rank_Z_f_array(t,1:6) = COPmatrix(:,2,w);
    end
    if ww == 0
    else
        group_rank_Z_f_array(t,7:12) = COPmatrix(:,2,ww);
    end
    if www == 0
    else
        group_rank_Z_m_array(t,1:6) = COPmatrix(:,2,www);
    end
    if wwww == 0
    else
        group_rank_Z_m_array(t,7:12) = COPmatrix(:,2,wwww);
    end
%%%
    if x == 0
    else
        group_rank_B_f_array(t,1:6) = COPmatrix(:,3,x);
    end
    if xx == 0
    else
        group_rank_B_f_array(t,7:12) = COPmatrix(:,3,xx);
    end
    if xxx == 0
    else
        group_rank_B_m_array(t,1:6) = COPmatrix(:,3,xxx);
    end
    if xxxx == 0
    else
        group_rank_B_m_array(t,7:12) = COPmatrix(:,3,xxxx);
    end
%%%
    if y == 0
    else
        group_rank_S_f_array(t,1:6) = COPmatrix(:,4,y);
    end
    if yy == 0
    else
        group_rank_S_f_array(t,7:12) = COPmatrix(:,4,yy);
    end
    if yyy == 0
    else
        group_rank_S_m_array(t,1:6) = COPmatrix(:,4,yyy);
    end
    if yyyy == 0
    else
        group_rank_S_m_array(t,7:12) = COPmatrix(:,4,yyyy);
    end
%%%
    if z == 0
    else
        group_rank_P_f_array(t,1:6) = COPmatrix(:,5,z);
    end
    if zz == 0
    else
        group_rank_P_f_array(t,7:12) = COPmatrix(:,5,zz);
    end
    if zzz == 0
    else
        group_rank_P_m_array(t,1:6) = COPmatrix(:,5,zzz);
    end
    if zzzz == 0
    else
        group_rank_P_m_array(t,7:12) = COPmatrix(:,5,zzzz);
    end
    
end

Top2_rank_f_array = cat(2,group_rank_G_f_array,group_rank_Z_f_array,group_rank_B_f_array,group_rank_S_f_array,group_rank_P_f_array);
Top2_rank_m_array = cat(2,group_rank_G_m_array,group_rank_Z_m_array,group_rank_B_m_array,group_rank_S_m_array,group_rank_P_m_array);
Top2_rank_gender = cat(1,Top2_rank_f_array,Top2_rank_m_array);
clear t v vv vvv vvvv w ww www wwww x xx xxx xxxx y yy yyy yyyy z zz zzz zzzz group_rank_G_f_array group_rank_Z_f_array group_rank_B_f_array group_rank_S_f_array group_rank_P_f_array group_rank_G_m_array group_rank_Z_m_array group_rank_B_m_array group_rank_S_m_array group_rank_P_m_array Top2_rank_f_array Top2_rank_m_array
% Top2_rank_gender is the final array


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xlswrite([directory 'Comfort_Comparisons.xlsx'],group_VAS_ALL_array,'VAS_All')
xlswrite([directory 'Comfort_Comparisons.xlsx'],group_VAS_indiv_array,'VAS_Indiv')
xlswrite([directory 'Comfort_Comparisons.xlsx'],group_rank_array,'Comfort_Ranking')
xlswrite([directory 'Comfort_Comparisons.xlsx'],group_YN_avg_array,'YesNo_Avg')
xlswrite([directory 'Comfort_Comparisons.xlsx'],group_YN_testday_array,'YesNo_TestDay')
xlswrite([directory 'Comfort_Comparisons.xlsx'],betweenshoes_array,'Between_Shoes')
xlswrite([directory 'Comfort_Comparisons.xlsx'],withinsubj_gender,'withinsubj_gender')
xlswrite([directory 'Comfort_Comparisons.xlsx'],VAS_All_gender,'VAS_all_gender')
xlswrite([directory 'Comfort_Comparisons.xlsx'],VAS_indiv_gender,'VAS_indiv_gender')
xlswrite([directory 'Comfort_Comparisons.xlsx'],Top2_rank_gender,'rank_gender')


% all variables to put into excel: withinshoe_gender VAS_indiv_gender VAS_All_gender Top2_rank_gender
