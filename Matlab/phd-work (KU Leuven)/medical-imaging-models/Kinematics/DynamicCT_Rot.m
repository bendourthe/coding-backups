clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DynamicCT_Rot: allows the use of the kinematics algorithm on multiple
% frames by using the function KinFunctionRad
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input: 
%           - stl files of each bone (static, use for a better definition
%               of the coordinate system + dynamic)
%           - number of frames
% 
% Output:
%           - creates 3 text files which gather the rotation angles
%           calculated for each bone (each line corresponds to one set of
%           rotation angle (e.g. between frame 1 and 2) and each column
%           corresponds to a specific direction (column 1: X, column 2: Y,
%           column 3: Z)
%
% Dependencies:        
%               num2str
%               STL_ReadFile.m 
%                   TRI_RemoveInvalidTriangles.m
%                   TRI_RemoveBadlyConnectedTriangles.m
%               txt2mat
%               PlacePoints3
%               GUI_PlotShells
%               GUI_VisualisationUI
%               arrow
%               KinFunctionRad
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Selection of the current directory (where the STL files are)
cd 'H:\Data Faes\Patient_07_Post\STLs\4D_AutoSeg_Output_Flex\stls\'
dir = 'H:\Data Faes\Patient_07_Post\STLs\4D_AutoSeg_Output_Flex\stls\';

% Name of the text files where the numerical data will be saved
radius_axis_available = 0; % loads data if available (=1), else, generate data
    filenameAX = 'RadAxes_Ext_Flex.txt';
    fileAX = [dir filenameAX];
    
% Number of frames
num_frames = 3;

% FILES

% Static
    V_MC1_stat = 'Patient_07_Neutral_MC1_reg.stl';
    V_Trap_stat = 'Patient_07_Neutral_Trap_reg.stl';
    V_Scaph_stat = 'Patient_07_Neutral_Scap_reg.stl';
    Rad_stat = 'Patient_07_Neutral_Rad_reg.stl';
% Dynamic
    MC1 = num2str([1:num_frames].','Patient_07_4D_Flex_Frame25_MC1_SE000000_Frame%04d.stl');
    Trap = num2str([1:num_frames].','Patient_07_4D_Flex_Frame25_Trap_SE000000_Frame%04d.stl');
    Scap = num2str([1:num_frames].','Patient_07_4D_Flex_Frame25_Scap_SE000000_Frame%04d.stl');
    Rad = num2str([1:num_frames].','Patient_07_4D_Flex_Frame25_Rad_SE000000_Frame%04d.stl');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Import STL files (static scan and frame 1)
    [F_MC1_stat, V_MC1_stat] =  STL_ReadFile([dir V_MC1_stat],true);
    [F_Trap_stat, V_Trap_stat] =  STL_ReadFile([dir V_Trap_stat],true);
    [F_Scaph_stat, V_Scaph_stat] =  STL_ReadFile([dir V_Scaph_stat],true);
    [F_Rad_stat, V_Rad_stat] =  STL_ReadFile([dir Rad_stat],true);
    
    [F_MC1_f1, V_MC1_f1] =  STL_ReadFile([dir MC1(1,:)],true);
    [F_Trap_f1, V_Trap_f1] =  STL_ReadFile([dir Trap(1,:)],true);
    [F_Scaph_f1, V_Scaph_f1] =  STL_ReadFile([dir Scap(1,:)],true);
    [F_Rad_f1, V_Rad_f1] =  STL_ReadFile([dir Rad(1,:)],true);
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CS on radius:
%   origin: lowest point on the distal border of the ulnar notch
%   y: straight down (few mm) on the proximal border of the ulnar notch
%   z: on the tip of the radial styloid
%   x: Perpendicual to the z-axis to form a right-handed CS (calculated)

if(radius_axis_available)
    CSRad = txt2mat(fileAX);
    
    O = CSRad(1,:);
    Y = CSRad(2,:);
    Z = CSRad(3,:);
else
    [ad_curve] = PlacePoints3({F_Rad_stat}, {V_Rad_stat}, {F_Rad_stat}, {V_Rad_stat}, 'Select origin, point on Y axis and point on Z axis (in that order!). End by pressing enter and exiting the figure');
    close all;
    % Retrieves coordinates of the points:
    O = ad_curve.points{end,1}(1,:);
    Y = ad_curve.points{end,1}(2,:);
    Z = ad_curve.points{end,1}(3,:);
    CSRad = [O;Y;Z];
    
    % Creates save file
    fid = fopen(fileAX,'w+');
    fprintf(fid,'%5.8f\t%5.8f\t%5.8f\r\n',CSRad);
    fclose(fid);
end

Origin = O;
Yaxis = (Y-O)/norm(Y-O);
Zaxis = (Z-O)/norm(Z-O);
Xaxis = cross(Yaxis,Zaxis)/norm(cross(Yaxis,Zaxis));
Zaxis = cross(Xaxis,Yaxis)/norm(cross(Xaxis,Yaxis));

% Converges to radius CS
RotMat = [Xaxis', Yaxis', Zaxis'];
loc2glob = [RotMat, Origin'; 0 0 0 1]; 
   % These tranformation matrices convert from local to global CS.
   % Invert them to go from global to local.
glob2loc = [RotMat', -(RotMat')*Origin'; 0 0 0 1];
R_glob2loc = glob2loc(1:3,1:3);
T_glob2loc = glob2loc(1:3,4);

% Transforms radius data to radius CS
V_Rad_stat_cs = R_glob2loc*V_Rad_stat.' + repmat(T_glob2loc,1,size(V_Rad_stat,1));
V_Rad_stat_cs = V_Rad_stat_cs.';
V_Rad_f1_cs = R_glob2loc*V_Rad_f1.' + repmat(T_glob2loc,1,size(V_Rad_f1,1));
V_Rad_f1_cs = V_Rad_f1_cs.';

if(~radius_axis_available)
% Plots the radius with the new coordinate system for validation
    h0 = figure;
    [obj, li, ax] = GUI_PlotShells(h0, {F_Rad_stat; F_Rad_f1},...
        {V_Rad_stat_cs; V_Rad_f1_cs},{ones(size(V_Rad_stat_cs,1),1),ones(size(V_Rad_f1_cs,1),1)});
    GUI_VisualisationUI(ancestor(ax, 'figure'), true, ax, true, true);
        
        % Plots the coordinate system and its origin
    hold on
    arrow([0;0;0],[30;0;0],10,70,30,5,'EdgeColor','r','FaceColor','r');
    arrow([0;0;0],[0;30;0],10,70,30,5,'EdgeColor','g','FaceColor','g');
    arrow([0;0;0],[0;0;30],10,70,30,5,'EdgeColor','b','FaceColor','b');
    set(h0,'numbertitle','off','name','Is the coordinate system correct? PRESS Continue if yes, CLOSE the figure otherwise.');
    uicontrol('Position',[20 20 200 40],'String','Continue',...
              'Callback','uiresume(gcbf)');
    uiwait(gcf);
    close(h0);
end  
        % Re-plot the radius coordinate system without buttons
        h0 = figure;
        [obj, li, ax] = GUI_PlotShells(h0, {F_Rad_stat; F_Rad_f1},...
            {V_Rad_stat_cs; V_Rad_f1_cs},{ones(size(V_Rad_stat_cs,1),1),ones(size(V_Rad_f1_cs,1),1)});
        GUI_VisualisationUI(ancestor(ax, 'figure'), true, ax, true, true);
            % Plots the coordinate system and its origin
        hold on
        arrow([0;0;0],[30;0;0],10,70,30,5,'EdgeColor','r','FaceColor','r');
        arrow([0;0;0],[0;30;0],10,70,30,5,'EdgeColor','g','FaceColor','g');
        arrow([0;0;0],[0;0;30],10,70,30,5,'EdgeColor','b','FaceColor','b');
        set(h0,'numbertitle','off','name','Radius coordinate system');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Automatic calculation of rotation angles for each frames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:num_frames-1;
    [Hom_Rad,Hom_Scaph,Hom_Trap,Hom_MC1,Rot_Scaph(i,:),Rot_Trap(i,:),Rot_MC1(i,:)] = ...
    KinFunctionRad(dir,CSRad,Rad(i,:),Rad(i+1,:),Scap(i,:),Scap(i+1,:), ...
                    Trap(i,:),Trap(i+1,:),MC1(i,:),MC1(i+1,:),i);
end


    % Save rotation angles calculated for the MC1
    save('Rotation angle MC1.txt','Rot_MC1','-ascii')
    % Save rotation angles calculated for the Trapezium
    save('Rotation angle Trap.txt','Rot_Trap','-ascii')
    % Save rotation angles calculated for the Scaphoid
    save('Rotation angle Scaph.txt','Rot_Trap','-ascii')