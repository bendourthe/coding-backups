clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Proximity_Dyn:
% - Plots a color map (previously calculated with the Proximity code) on a
% single STL. The corresponding image can be saved with a single
% orientation and used to make an animation showing the evolution of the
% proximity patterns during dynamic scanning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dependencies:        
%               STL_ReadFile.m 
%                   TRI_RemoveInvalidTriangles.m
%                   TRI_RemoveBadlyConnectedTriangles.m
%               GUI_PlotShells
% Input: 
%           MC1: STL file corresponding to first metacarpal
%           Trap: STL file corresponding to trapezium
% 
% Output:
%           Figure 1: Shows the proximity pattern (full range) of the MC1
%           Figure 2: Shows the proximity pattern (full range) of the Trap   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Selection of the current directory (where the STL files are)
dir = 'C:\Users\u0096636\Documents\Professional\PhD KU Leuven\Data\Forward Model\STL\';

% FILES
Bone1 = 'SCAN3_N_mc1.stl';
Bone2 = 'SCAN3_N_trap.stl';
filenameColorMC1 = 'Color_Trap.txt'; % Imports the corresponding color map
    fileColorMC1 = [dir filenameColorMC1];
filenameColorTrap = 'Color_Trap.txt'; % Imports the corresponding color map
    fileColorTrap = [dir filenameColorTrap];
    
% Read STL files
[F_Bone1, V_Bone1] =  STL_ReadFile([dir Bone1],true);
[F_Bone2, V_Bone2] =  STL_ReadFile([dir Bone2],true);

% Load the color matrices
colrstl1 = txt2mat(fileColorMC1);    
colrstl2 = txt2mat(fileColorTrap);

% Sets the range of the color map
dist_max = 3;       % in mm

% Plots the MC1 with the color map representing the proximity pattern
proximity1 = figure; 
    [obj, li, ax] = GUI_PlotShells(proximity1, {F_Bone1}, {V_Bone1},...
            {ones(size(V_Bone1,1),1)},[0,0,1]);
box off
view(-160,-45);
set(0,'defaultfigurecolor',[1 1 1])
gca.Color = [1 1 1];
hold on
patch('Faces',F_Bone1,'Vertices',V_Bone1, ...
    'FaceColor','interp', ...
    'FaceVertexCData',colrstl1, ...
    'EdgeColor', 'interp', ...
    'EdgeAlpha', 0, ...
    'CDataMapping', 'scaled',...
    'AmbientStrength', 0.4, ...
    'DiffuseStrength', 0.8, ...
    'SpecularStrength', 0.2, ...
    'SpecularColorReflectance', 0.5, ...
    'FaceLighting', 'gouraud')
hold on
step_h = (dist_max-1)/8;
scale_h = [dist_max,dist_max-step_h,dist_max-2*step_h,dist_max-3*step_h,...
    dist_max-4*step_h,dist_max-5*step_h,dist_max-6*step_h,dist_max-7*step_h,dist_max-8*step_h];
h=colorbar('YTickLabel',{scale_h});
set(h,'ylim',[0.1 0.9]);
title = get(h,'Title');
titleString = 'Distance (mm)';
set(title ,'String',titleString,'FontWeight','bold');
set(gcf,'numbertitle','off','name','MC1 - Proximity pattern (Full range)');

% Plots the trapezium with the color map representing the proximity
% pattern
proximity2 = figure; 
    [obj, li, ax] = GUI_PlotShells(proximity2, {F_Bone2}, {V_Bone2},...
            {ones(size(V_Bone2,1),1)},[0,0,1]);
box off
view(20,45);
set(0,'defaultfigurecolor',[1 1 1])
hold on
patch('Faces',F_Bone2,'Vertices',V_Bone2, ...
    'FaceColor','interp', ...
    'FaceVertexCData',colrstl2, ...
    'EdgeColor', 'interp', ...
    'EdgeAlpha', 0, ...
    'CDataMapping', 'scaled',...
    'AmbientStrength', 0.4, ...
    'DiffuseStrength', 0.8, ...
    'SpecularStrength', 0.2, ...
    'SpecularColorReflectance', 0.5, ...
    'FaceLighting', 'gouraud')
hold on
step_h = (dist_max-1)/8;
scale_h = [dist_max,dist_max-step_h,dist_max-2*step_h,dist_max-3*step_h,...
    dist_max-4*step_h,dist_max-5*step_h,dist_max-6*step_h,dist_max-7*step_h,dist_max-8*step_h];
h=colorbar('YTickLabel',{scale_h});
set(h,'ylim',[0.1 0.9]);
title = get(h,'Title');
titleString = 'Distance (mm)';
set(title ,'String',titleString,'FontWeight','bold');
set(gcf,'numbertitle','off','name','Trapezium - Proximity pattern (Full range)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%