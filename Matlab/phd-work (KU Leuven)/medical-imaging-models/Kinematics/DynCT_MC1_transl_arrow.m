clear all; close all;
% Function: Plots the trapezium and shows how the MC1 is translating
% between frames by drawing arrows (each representing a translation
% vector between 2 frames) which originates centroid of the MC1
%
% Dependencies:        
%               txt2mat.m
%               STL_ReadFile.m
%               GUI_PlotShells.m
%               arrow.m
% 
% Input: STl files of the trapezium and first metacarpal(MC1) (static scan)
%        Text files with the transformation matrices calculated between
%        each frame
%
% Output: Trapezium with the MC1 path represented by arrows.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Control panel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Folder with files
dir = 'C:\Users\u0096636\Documents\Professional\PhD KU Leuven\Data\Study_Faes_4D_CT\STL\registration\'; 

% Text files containing the transformation matrices (create as many copies
% as there are text files)
trans_f12 = 'RotTransData_f1f9.txt';
f12 = [dir trans_f12];
trans_f23 = 'RotTransData_f9f16.txt';
f23 = [dir trans_f23];
trans_f34 = 'RotTransData_f16f25.txt';
f34 = [dir trans_f34];

% Import the STL of the trapezium and MC1 from the static scan
Trap = '611L_reg1_trp.stl';
MC1 = '611L_reg1_mc1.stl';

% Read STL files
[F_Trap, V_Trap] =  STL_ReadFile([dir Trap],true);
[F_MC1, V_MC1] =  STL_ReadFile([dir MC1],true);

% Proximity distance
d_prox = 1.5; % in mm

% Import the transformation matrices text file
T_f12 = txt2mat(f12);
T_f23 = txt2mat(f23);
T_f34 = txt2mat(f34);

% Calculation of the matrix Di which contains the Euclidean distances between
% each vertices of each bone
Di = pdist2(V_MC1,V_Trap);
        
% Save the indexes of the vertices of each bone which are located
% in the proximity zone
[Di_MC1,Di_Trap] = find(Di<d_prox);
idxDef = [Di_MC1,Di_Trap];
idxDef(:,any(idxDef==0,1)) = []; % Removes all the zero in the idxDef Matrix
idx_MC1 = unique(idxDef(:,1));
idx_Trap = unique(idxDef(:,2));
        
% Calculates the coordinates of the centroid of the contact area
O = mean(V_MC1(idx_MC1,:));

% Plots the trapezium from the static scan with arrows representing each
% translation vector calculated between each frame (origin: centroid of the
% MC1)

h = figure;
    [obj, li, ax] = GUI_PlotShells(h, {F_Trap; F_MC1}, {V_Trap; V_MC1},{ones(size(V_Trap,1),1), ones(size(V_MC1,1),1)});
hold on
    plot3(O(1),O(2),O(3),'Marker','o','MarkerSize',10,'Color','k','MarkerFaceColor','k');
    arrow([O(1,1);O(1,2);O(1,3)],[T_f12(4,4);T_f12(5,4);T_f12(6,4)],10,70,30,5);
    arrow([T_f12(4,4);T_f12(5,4);T_f12(6,4)],[T_f23(4,4);T_f23(5,4);T_f23(6,4)],10,70,30,5);
    arrow([T_f23(4,4);T_f23(5,4);T_f23(6,4)],[T_f34(4,4);T_f34(5,4);T_f34(6,4)],10,70,30,5);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%