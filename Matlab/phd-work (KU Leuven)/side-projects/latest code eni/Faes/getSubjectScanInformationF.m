function [subjectInfo] = getSubjectScanInformationF(CTVolume)
%%
% dicomFname = 'C:\Functional\E01424\DICOMS\E01424_15\E01424_01_001.dicom';

% subjectDir = 'P:\Data\CMCOA\Pilot\E14168';
% series = '15'
%folders = dir(fullfile(subjectDir,sprintf('S*%s',side)));
% parameter = regexp(subjectDir,'^(.+)\\(E\d{5})[\s\\]*$','tokens');
% baseDir = parameter{1}{1};
% subject = parameter{1}{2};
% 
% dicomFname = fullfile(subjectDir,'DICOMS',sprintf('%s_%s',subject,series),sprintf('%s_%s_002.dicom',subject, series));
% scanNumber = series;
% if (exist(dicomFname,'file')~=2),
 %s   dicomFname = fullfile(subjectDir,'DICOMS',sprintf('%s_%s',subject, series),sprintf('%s_02_002.dicom',subject));
%     scanNumber = 1;
% end;
% 
% if (exist(dicomFname,'file')~=2),
%     error('Invalid dicom file specified: %s\n',dicomFname);
% end;
d=dir(CTVolume);
info = dicominfo(fullfile(CTVolume,d(3,1).name));

subjectInfo.radiusLocation = info.SliceLocation;
subjectInfo.firstSliceLocation = info.SliceLocation;
%subjectInfo.firstSliceFilename = dicomFname;
subjectInfo.x_voxel = info.PixelSpacing(1);
subjectInfo.y_voxel = info.PixelSpacing(2);
subjectInfo.z_voxel = info.SliceThickness;

subjectInfo.x_size = double(info.Width);
subjectInfo.y_size = double(info.Height);

subjectInfo.x_RealSize = double(info.Width) * info.PixelSpacing(1);
subjectInfo.y_RealSize = double(info.Height) * info.PixelSpacing(2);

subjectInfo.PatientPosition = info.PatientPosition;

% try and find the last slice in the image
% pattern = fullfile(subjectDir,'DICOMS',sprintf('%s_%s',subject, series),'*.dicom');
% files = dir(pattern);
% 
% regpattern = sprintf('^%s_%02d_(\\d{3})\\.dicom$',subject,scanNumber);
% lastIndex = 1;
% 
% for i=1:length(files),
%     tokens = regexpi(files(i).name,regpattern,'tokens');
%     if (size(tokens,1)==1)
%         currentIndex = str2double(tokens{1}{1});
%         if (currentIndex > lastIndex)
%             lastFilename =files(i).name;
%             lastIndex = currentIndex;
%         end;
%     end;    
% end;
% lastFilename =files(i).name;
% dicomFname = fullfile(subjectDir,'DICOMS',sprintf('%s_%s',subject, series),lastFilename);
info = dicominfo(fullfile(CTVolume,d(size(d,1),1).name));
subjectInfo.lastSliceLocation = info.SliceLocation;
%subjectInfo.lastSliceFilename = dicomFname;
