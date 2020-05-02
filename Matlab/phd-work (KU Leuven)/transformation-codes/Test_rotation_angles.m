clear all; close all;
% Function: Define the radius coordinate system (CS), perform a 4-points
% registration based on an ICP algorithm to calculate the transformation
% matrices from a position 1 to a position 2, calculates the helical axis
% of each joint and plots them, then finally calculates the angle between
% each helical axis and the (Y-Z) plane of the radius CS as well as the
% coordinates of the intersection between those two components.
%
% Dependencies:        
%               txt2mat.m
%               STL_ReadFile.m
%               PlacePoints3.m
%               RegisterBones.m
%                   TRI_RemoveInvalidTriangles.m
%                       DeleteUnreferencedElements.m% 
%                   TRI_RemoveBadlyConnectedTriangles.m
%                       TRI_Edges.m
%                       StackEqualElementIndices.m
%                           RunLengthDecode.m
%                           IncrementalRuns.m
%                       TRI_Normals.m
%                           VectorNorms.m
%                       DeleteUnreferencedElements.m
%                   TRI_Areas.m
%                   KIN_DirectedSpringRotation.m
%               GUI_PlotShells.m
%               IntersectLineAndPlane.m 
%                   DistanceFromVertexToPlane.m 
%                       TRI_Normals.m	
%               VectorNorms.m
%               screw.m
%               FindPivotPoint
%               DistanceFromVertexToLine
%               arrow.m
% 
% Input: STl files of the radius, scaphoid, trapezium and first metacarpal
%       (MC1).
%
% Output: txt-files with transformations, interia axes and coordinate
% systems. Figures showing if the efficiency of the registration, the
% helical axis for each joint, the different inclination angle and the
% intersection of each helical axis with the (X-Y) plane of the radius CS
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Controle panel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Folder with files
dir = 'C:\Users\u0096636\Documents\Professional\PhD KU Leuven\Data\CAD Tests\Rot_Cube\'; 
data_available = 0; % loads data if available (=1), else, generate data
    filename = 'RotTransData.txt';
    file = [dir filename];
radius_axis_available = 0; % loads data if available (=1), else, generate data
    filenameAX = 'RadAxes.txt';
    fileAX = [dir filenameAX];
rot_angle_available = 0; % loads data if available (=1), else, generate data
    filenameRotAngle = 'RotAngle_pzyx.txt';
    fileRotAngle = [dir filenameRotAngle];

% FILES

% First position (referenced as P1)
Rad_P1 = 'cube.stl';

% Second position (referenced as P2)
Rad_P2 = 'Rot_m45zyx.stl';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CS on radius:
%   origin: lowest point on the distal border of the ulnar notch
%   y: straight down (few mm) on the proximal border of the ulnar notch
%   z: on the tip of the radial styloid
%   x: Perpendicual to the z-axis to form a right-handed CS (calculated)

if(radius_axis_available)
    A = txt2mat(fileAX);
    
    O = A(1,:);
    Y = A(2,:);
    Z = A(3,:);
else
    [F_Rad_P1, V_Rad_P1] =  STL_ReadFile([dir Rad_P1],true);
    [ad_curve] = PlacePoints3({F_Rad_P1}, {V_Rad_P1}, {F_Rad_P1}, {V_Rad_P1}, 'Select origin, point on Y axis and point on Z axis (in that order!). End by pressing enter and exiting the figure');
    close all;
    % Retrieves coordinates of the points:
    O = ad_curve.points{end,1}(1,:);
    Y = ad_curve.points{end,1}(2,:);
    Z = ad_curve.points{end,1}(3,:);
    
    % Creates save file
    fid = fopen(fileAX,'w+');
    fprintf(fid,'%5.8f\t%5.8f\t%5.8f\r\n',[O,Y,Z]);
    fclose(fid);
end

Origin = O;
Yaxis = (Y-O)/norm(Y-O);
Zaxis = (Z-O)/norm(Z-O);
Xaxis = cross(Yaxis,Zaxis)/norm(cross(Yaxis,Zaxis));
Zaxis = cross(Xaxis,Yaxis)/norm(cross(Xaxis,Yaxis));

% Converges to radius CS
RotMat = [Xaxis', Yaxis', Zaxis'];
loc2glob = [RotMat, Origin';
     0 0 0 1]; 
   % These tranformation matrices convert from local to global CS.
   % Invert them to go from global to local.
glob2loc = [RotMat', -(RotMat')*Origin'; 0 0 0 1];
R_glob2loc = glob2loc(1:3,1:3);
T_glob2loc = glob2loc(1:3,4);

% Transformes data to radius CS
V_Rad_P1 = R_glob2loc*V_Rad_P1.' + repmat(T_glob2loc,1,size(V_Rad_P1,1));
V_Rad_P1 = V_Rad_P1.';

% Plots the radius with the new coordinate system for validation
    h0 = figure;
    [obj, li, ax] = GUI_PlotShells(h0, {F_Rad_P1},...
        {V_Rad_P1},...
        {ones(size(V_Rad_P1,1),1)});
    GUI_VisualisationUI(ancestor(ax, 'figure'), true, ax, true, true);
        
        % Plots the coordinate system and its origin
    hold on
    arrow([0;0;0],[30;0;0],10,70,30,5,'EdgeColor','r','FaceColor','r');
    arrow([0;0;0],[0;30;0],10,70,30,5,'EdgeColor','g','FaceColor','g');
    arrow([0;0;0],[0;0;30],10,70,30,5,'EdgeColor','b','FaceColor','b');
    
        % Plots 2 buttons: 'Bad CS' if the coordinate system looks wrong,
        % this will stop the code, and 'Continue' if the coordinate system
        % looks right, will continue to run the rest of the code
    b1 = uicontrol('Position',[20 20 200 40],'String','Continue',...
              'Callback','uiresume(gcbf)');
    b2 = uicontrol('Position',[20 60 200 40],'String','Bad CS',...
              'Callback','close(gcbf)');
    uiwait(gcf);
    close(h0);

% Searches for radius fragment
if(data_available)
    % Loads data
    A = txt2mat(file);
        % Loads all the STL files
    [F_Rad_P1, V_Rad_P1] =  STL_ReadFile([dir Rad_P1],true);
    [F_Rad_P2, V_Rad_P2] =  STL_ReadFile([dir Rad_P2],true);

        % Associates the different parts of the RotTransData file to the
        % corresponding rotation and translation matrices for each bone
    R_Rad = A(1:3,1:3);
    T_Rad = A(1:3,4);
    H_Rad = [R_Rad, T_Rad; 0 0 0 1];
else
    
    % Loads data
       fid = fopen(file,'w+');
        % Loads all the STL files
    [F_Rad_P1, V_Rad_P1] =  STL_ReadFile([dir Rad_P1],true);
    [F_Rad_P2, V_Rad_P2] =  STL_ReadFile([dir Rad_P2],true);
    
    % Registers radius in position 2 with the radius in position 1
    [R_Rad, T_Rad, ERROR_Rad] = RegisterBones(F_Rad_P2, V_Rad_P2,F_Rad_P1, V_Rad_P1,1);
    H_Rad = [R_Rad, T_Rad; 0 0 0 1];
    
    % Saves results to file
    W = reshape([R_Rad, T_Rad].',1,[]);
    fprintf(fid,'%5.8f\t%5.8f\t%5.8f\t%5.8f\r\n',W); 
end

% Calculation of the rotation angles along each axis of the radius
% coordinate system for each bone. Each rotation angle is expressed in the
% radius coordinate system relatively to the more proximal bone.

    % Radius
    
    Rot_rad_x = asin(R_Rad(3,2)/cos(asin(-R_Rad(3,1))))*180/pi;
    Rot_rad_y = asin(-R_Rad(3,1))*180/pi;
    Rot_rad_z = asin(R_Rad(2,1)/cos(asin(-R_Rad(3,1))))*180/pi;
    Rot_rad = [Rot_rad_x,Rot_rad_y,Rot_rad_z];
    
    % Saves to file
        fid2 = fopen(fileRotAngle,'w+');
        fprintf(fid2,'%5.8f\t%5.8f\t%5.8f\r\n',Rot_rad);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%