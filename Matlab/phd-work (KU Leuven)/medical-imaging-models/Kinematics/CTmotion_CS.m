clear all; close all;
% Function: Select 3 landmarks on the radius to define a new coordinate
% system and plots the result on a 3D model.
%
% Dependencies:        
%               txt2mat.m
%               STL_ReadFile.m
%               PlacePoints3.m
%               arrow.m
%               GUI_PlotShells.m
%               arrow.m
%
% Input: STL file of the radius.
%
% Output: Figure of the radius with the corresponding coordinate system
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Controle panel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Folder with files
dir = 'D:\PhD KU Leuven\Data\Study_Brown_Rig\New Segmentation\Scan1\STL\Priscilla\'; 

% FILES

% First position (referenced as P1)
Rad_P1 = 'SCAN1_ADD_rad.stl';
% Second position (referenced as P2)
Rad_P2 = 'SCAN1_ABD_rad.stl';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CS on radius:
%   origin: lowest point on the distal border of the ulnar notch
%   y: straight down (few mm) on the proximal border of the ulnar notch
%   z: on the tip of the radial styloid
%   x: Perpendicual to the z-axis to form a right-handed CS (calculated)

    [F_Rad_P1, V_Rad_P1] =  STL_ReadFile([dir Rad_P1],true);
    [ad_curve] = PlacePoints3({F_Rad_P1}, {V_Rad_P1}, {F_Rad_P1}, {V_Rad_P1}, 'Select origin, point on Y axis and point on Z axis (in that order!). End by pressing enter and exiting the figure');
    close all;
    % Retrieves coordinates of the points:
    O = ad_curve.points{end,1}(1,:);
    Y = ad_curve.points{end,1}(2,:);
    Z = ad_curve.points{end,1}(3,:);
    
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

% Transformes all data to radius CS (incl. ransforming dynamic scan to static CS)
V_Rad_P1 = R_glob2loc*V_Rad_P1.' + repmat(T_glob2loc,1,size(V_Rad_P1,1));
V_Rad_P1 = V_Rad_P1.';

% Plots STl files and the new coordinate system
    h5 = figure;
    [obj, li, ax] = GUI_PlotShells(h5, {F_Rad_P1},...
        {V_Rad_P1},...
        {ones(size(V_Rad_P1,1),1)});
        
        % Plots the coordinate system and its origin
    hold on
    arrow([0;0;0],[30;0;0],10,70,30,5,'EdgeColor','r','FaceColor','r');
    arrow([0;0;0],[0;30;0],10,70,30,5,'EdgeColor','g','FaceColor','g');
    arrow([0;0;0],[0;0;30],10,70,30,5,'EdgeColor','b','FaceColor','b');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%