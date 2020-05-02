function [] = autoAlignAllMimicsOutputs(CTVolume, bone_files, Side, flippedScan)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% autoAlignAllMimicsOutput take the raw WRL files output from Mimics and 
% aligns them with theoriginal image files. In order to do so, this code 
% attempts to determine where mimics has placed the wrists
%
%    autoAlignMimicsOutput(subjectDir, Side) - will run against the subject
%      stored in 'subjectDir', and the left or right wrist depending upon
%      the side specified. In this mode, the program will ask the user if
%      the original scans had to be flipped. This is in relation to how
%      |dicom_to_mri| originally setup the MRI image. The radius (forearm),
%      must be in the first slice, but if it is not, then this must be
%      corrected using rotateMRI(scan,'frontback'), to make it so. This
%      changes the orientation of the scan, so must be accounted for.
%
%    autoAlignMimicsOutput(subjectDir, Side, flippedScan) - Same as above,
%      however, this time you specify if the scan was flipped on the
%      command line. 1 - says it was flipped (like answering 'N'), 0 - says
%      it was not flipped (like answering 'Y')
%
%      Side - 'L' or 'R' to signify left or right. To create left wrists,
%      mimics negates the Y values, but CS expects that the X values will
%      be negated.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% % Side = upper(Side);
% % rightWrist = (upper(Side)=='R');

if (exist('flippedScan','var')~=1)
    reply = input('Did the neutral scan originally come with the forearm in the first slice? Y/N [Y]: ', 's');
    if (isempty(reply)), reply = 'Y';
   
    end;
    
    if (upper(reply(1))=='N')
        flippedScan = 1;
    else
        flippedScan = 0;
    end;
    
else 
    flippedScan = 0;
end;

% Get Subject Information and Correct position info
try
%     warning off Wrist:NoStack;
%     warning off Wrist:NoIV;
%     warning off Wrist:NoInertia;
%     info = getSubjectFilePathE(subjectDir, Side);
    
    scaninfo = getSubjectScanInformationF(CTVolume);
catch
    warning('Wrist:DicomNotFound','Unable to locate dicom files. Unknown subject orientation. Continuing in default CT space.');
end
wrml_files=fullfile(bone_files, 'WRL.files');
 % find the neutral radius file
radfiles = dir(fullfile(wrml_files,'*rad*.wrl'));
if (length(radfiles)~=1)
    error('AutoAlignMimics:UnknownRadius','Unable to locate radius wrl file: needed for aligning output');
end;
radWRLfname = fullfile(wrml_files,radfiles(1).name);
% radWRLfname = fullfile(subjectDir,neutralSeries,'WRL.files',sprintf('rad.wrl',neutralSeries(2:end)));
[pts conn] = read_vrml_fastF(radWRLfname);

% determine if radius was export from Mimics as a left or right wrist.
% Mimics defaults to negating y-values if it was setup as a left wrist.
if (min(pts(:,2)) > 0 && Side=='L')
    error('Wrist:WrongSide', 'Wrist specified as a Left wrist, but mimics export looks like a right wrist (y values > 0)');
elseif (max(pts(:,2)) < 0 && Side=='R')
    error('Wrist:WrongSide', 'Wrist specified as a Right wrist, but mimics export looks like a left wrist (y values < 0)');
end;
%okay, so we look like we are on the right side.

% try and locate the bottom of the radius stack file, one end of it should
% be very close to either the first or last slice..... hopefully the first
% slice, or there are problems....
minRad = min(pts(:,3));
maxRad = max(pts(:,3));
minFirst = abs(minRad - scaninfo.firstSliceLocation);
minLast = abs(minRad - scaninfo.lastSliceLocation);
maxFirst = abs(maxRad - scaninfo.firstSliceLocation);
maxLast = abs(maxRad - scaninfo.lastSliceLocation);

if (minFirst < minLast && minFirst < maxFirst && minFirst < maxLast)
    % then the min of the radius is closest to the first slice
    radLocation = scaninfo.firstSliceLocation;
    mcLocation = scaninfo.lastSliceLocation;
    
elseif (minLast < maxFirst && minLast < maxLast)
    % then the min of the radius is closest to the last slice
    radLocation = scaninfo.lastSliceLocation;
    mcLocation = scaninfo.firstSliceLocation;
    
elseif (maxFirst < maxLast)
    % then the max of the radius is closest to the first slice
    radLocation = scaninfo.firstSliceLocation;
    mcLocation = scaninfo.lastSliceLocation;
    
else
    % then the max of the radius is closest to the last slice...
    radLocation = scaninfo.lastSliceLocation;
    mcLocation = scaninfo.firstSliceLocation;
end;


% Need to identify orientation of dicom relative to mri scan.... this is
% harder then it looks; so many options....

% loop through each of the old WRL.files
%moved_wrml_files = fullfile(info.neutralSeriesPath, 'WRL.files');
% ivFiles = fullfile(info.neutralSeriesPath, 'IV.files');
% stackFiles = fullfile(info.neutralSeriesPath, 'Stack.files');
% warning off MATLAB:MKDIR:DirectoryExists;

moved_wrml_files=fullfile(bone_files, 'movedWRL.files');
iv_files=fullfile(bone_files, 'IV.files');
stack_files=fullfile(bone_files, 'Stack.files');
mkdir(fullfile(bone_files, 'movedWRL.files'));
copyfile('meshconv.exe',moved_wrml_files)
mkdir(fullfile(bone_files, 'Stack.files'))
mkdir(fullfile(bone_files, 'IV.files'))

files = dir(fullfile(wrml_files,'*.wrl'));
for i=1:length(files),
    name = files(i).name;

     parameter = regexpi(name, '^.+_([a-zA-Z0-9]+) \d+_\d{3}\.wrl$','tokens');
     shortName = lower(parameter{1}{1});

         IVname=[shortName '15' Side '.iv'];
         WRLname=[shortName '15' Side '.wrl'];
    sname=[shortName '15' Side '.stack'];
    startFname = fullfile(wrml_files, name);
    endFname_iv = fullfile(iv_files, IVname);
    endFname_wrl = fullfile(moved_wrml_files, WRLname);
    stackFname=fullfile(stack_files, sname);
%     stackFname = fullfile(stackFiles, sprintf('%s%s.stack',shortName,info.neutralSeries(2:end)));
%     fprintf(1,'%s\t=> %s\n',startFname,endFname);
    
    % build the required RT
    RT = [eye(3); 0 0 0];
    RTx = [eye(3); 0 0 0];
    RTy = [eye(3); 0 0 0];
    RTz = [eye(3); 0 0 0];
    if (flippedScan)
        if (mcLocation < radLocation)
            if (rightWrist)
                fprintf('Case 1');
               % error('Wrist:IncompleteMode','This combination has not yet been programmed'); 
                 RT=RT;

            else
                fprintf('Case 2');
              RT(4,1:3) = [-scaninfo.x_RealSize scaninfo.y_RealSize -radLocation];
               RTz(1:3,1:3)=[cos(pi) -sin(pi) 0; sin(pi) cos(pi) 0; 0 0 1];
              %  RT=RT;
            end;
        else
              if (rightWrist)  
                  fprintf('Case 3');
                  RTy(1:3,1:3)=[cos(pi) 0 sin(pi);0 1 0; -sin(pi) 0 cos(pi)];
                 RTz(1:3,1:3)=[cos(pi) -sin(pi) 0; sin(pi) cos(pi) 0; 0 0 1];
                %  RTx(1:3,1:3)=[1 0 0; 0 cos(pi) -sin(pi); 0 sin(pi) cos(pi)];
                   RT(4,1:3) = [0 scaninfo.y_RealSize maxRad];
         %regular
       %RT(4,1:3) = [scaninfo.x_RealSize scaninfo.y_RealSize -radLocation];
       
       %RTz(1:3,1:3)=[cos(pi) -sin(pi) 0; sin(pi) cos(pi) 0; 0 0 1];
           %%x-rot
          % RTx(1:3,1:3)=[1 0 0; 0 cos(pi) -sin(pi); 0 sin(pi) cos(pi)];
           %%y-rot
         % RTy(1:3,1:3)=[cos(pi) 0 sin(pi);0 1 0; -sin(pi) 0 cos(pi)];
           
              else
                  fprintf('Case 4');
                 RT(4,1:3) = [-scaninfo.x_RealSize scaninfo.y_RealSize -radLocation];
                RT(4,1:3)=[0 0 radLocation];
              end
        end;
    else
        if (mcLocation < radLocation)
            if (rightWrist) 
                fprintf('Case 5');
                %error('Wrist:IncompleteMode','This combination has not yet been programmed');
                %%z-rot
%                  RTz(1:3,1:3)=[cos(pi) -sin(pi) 0; sin(pi) cos(pi) 0; 0 0 1];
                %%x-rot
             RTx(1:3,1:3)=[1 0 0; 0 cos(pi) -sin(pi); 0 sin(pi) cos(pi)];
                %%y-rot
%                  RTy(1:3,1:3)=[cos(pi) 0 sin(pi);0 1 0; -sin(pi) 0 cos(pi)];
            RT(4,1:3) = [0 scaninfo.y_RealSize radLocation-maxFirst]; %P:\Data\CMCOA\NORMALS\BN00126\E15294
          %RT(4,1:3) = [scaninfo.x_RealSize/2 scaninfo.y_RealSize/2 radLocation];
%             RT=RT;
            
            else
                fprintf('Case 6');
               %   RTy(1:3,1:3)=[cos(pi) 0 sin(pi);0 1 0; -sin(pi) 0 cos(pi)];
              %   RTz(1:3,1:3)=[cos(pi) -sin(pi) 0; sin(pi) cos(pi) 0; 0 0 1];
                 RTx(1:3,1:3)=[1 0 0; 0 cos(pi) -sin(pi); 0 sin(pi) cos(pi)];
                 RT(4,1:3) = [-scaninfo.x_RealSize/2 scaninfo.y_RealSize/2 radLocation];
            end
        else
            fprintf('Case 7');
           %error('Wrist:IncompleteMode','This combination has not yet been programmed');
                           %%z-rot
%                  RTz(1:3,1:3)=[cos(pi) -sin(pi) 0; sin(pi) cos(pi) 0; 0 0 1];
                %%x-rot
             RTx(1:3,1:3)=[1 0 0; 0 cos(pi) -sin(pi); 0 sin(pi) cos(pi)];
                %%y-rot
%                 RTy(1:3,1:3)=[cos(pi) 0 sin(pi);0 1 0; -sin(pi) 0 cos(pi)];
% RT=RT;
%   RT(4,1:3) = [-127	94.5 -minRad];
        end;        
    end;
    
    % convert from wrl => iv
   RT=RT;
   RTx=RTx;
   RTy=RTy;
   RTz=RTz;
   
   transforms = [RTx;RTy;RTz;RT];
   align_points_to_IV = fullfile(moved_wrml_files, 'align_points_to_IV.dat');
   dlmwrite(align_points_to_IV, transforms, '\t');
   
   transformIVF(startFname, RTx, endFname_iv);
    transformIVF(endFname_iv, RTy, endFname_iv);
     transformIVF(endFname_iv, RTz, endFname_iv); 
    transformIVF(endFname_iv, RT, endFname_iv);
    
   transformIVF(startFname, RTx, endFname_wrl);
    transformIVF(endFname_wrl, RTy, endFname_wrl);
     transformIVF(endFname_wrl, RTz, endFname_wrl); 
    transformIVF(endFname_wrl, RT, endFname_wrl);
%      % convert from iv => stack
    iv2stack_fastF(endFname_wrl, stackFname);
end

cd(moved_wrml_files)
files = dir('*.wrl');
for i=1:length(files),
    name = files(i).name;
     system(['meshconv.exe -c stl -tri ' name]);   
end
end
