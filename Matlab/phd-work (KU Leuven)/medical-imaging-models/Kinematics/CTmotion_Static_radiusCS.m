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
% Output: txt-files with transformations, intertia axes and coordinate
% systems. Figures showing if the efficiency of the registration, the
% helical axis for each joint, the different inclination angle and the
% intersection of each helical axis with the (X-Y) plane of the radius CS
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Controle panel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Folder with files
dir = 'H:\Data\Study_Brown_Rig\Segmentation\Scan3\STL\'; 
data_available = 0; % loads data if available (=1), else, generate data
    filename = 'RotTransData_Ext_Flex.txt';
    file = [dir filename];
radius_axis_available = 0; % loads data if available (=1), else, generate data
    filenameAX = 'RadAxes_Ext_Flex.txt';
    fileAX = [dir filenameAX];
inertial_axis_available = 0; % loads data if available (=1), else, generate data
    filenameInertialAX = 'InertialAx_Ext_Flex.txt';
    fileInertialAX = [dir filenameInertialAX];
helical_axis_available = 0; % loads data if available (=1), else, generate data
    filenameHelicalAX = 'HelicalAx_Ext_Flex.txt';
    fileHelicalAX = [dir filenameHelicalAX];
rot_angle_available = 0; % loads data if available (=1), else, generate data
    filenameRotAngle = 'RotAngle_Ext_Flex.txt';
    fileRotAngle = [dir filenameRotAngle];

% FILES
    % Note: the code works so that P2 is registered on P1, so if you want
    % to study the motion Ext-Flex (for instance), you should put Ext as P2
    % and Flex as P1 (and vice-versa if you want to study Flex-Ext motion)
% First position (referenced as P1)
MC1_P1 = 'SCAN3_FLEX_mc1.stl';
Rad_P1 = 'SCAN3_FLEX_rad.stl';
Scaph_P1 = 'SCAN3_FLEX_scaph.stl';
Trap_P1 = 'SCAN3_FLEX_trap.stl';

% Second position (referenced as P2)
MC1_P2 = 'SCAN3_EXT_mc1_reg.stl';
Rad_P2 = 'SCAN3_EXT_rad.stl';
Scaph_P2 = 'SCAN3_EXT_scaph.stl';
Trap_P2 = 'SCAN3_EXT_trap.stl';

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

if(radius_axis_available)
else
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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Searches for radius fragment
if(data_available)
    % Loads data
    A = txt2mat(file);
        % Loads all the STL files
    [F_MC1_P1, V_MC1_P1] =  STL_ReadFile([dir MC1_P1],true);
    [F_Rad_P1, V_Rad_P1] =  STL_ReadFile([dir Rad_P1],true);
    [F_Scaph_P1, V_Scaph_P1] =  STL_ReadFile([dir Scaph_P1],true);
    [F_Trap_P1, V_Trap_P1] =  STL_ReadFile([dir Trap_P1],true);
    [F_MC1_P2, V_MC1_P2] =  STL_ReadFile([dir MC1_P2]',true);
    [F_Rad_P2, V_Rad_P2] =  STL_ReadFile([dir Rad_P2],true);
    [F_Scaph_P2, V_Scaph_P2] =  STL_ReadFile([dir Scaph_P2],true);    
    [F_Trap_P2, V_Trap_P2] =  STL_ReadFile([dir Trap_P2],true);
        % Associates the different parts of the RotTransData file to the
        % corresponding rotation and translation matrices for each bone
    R_Rad = A(1:3,1:3);
    T_Rad = A(1:3,4);
    H_Rad = [R_Rad, T_Rad; 0 0 0 1];
    R_MC1 = A(4:6,1:3);
    T_MC1 = A(4:6,4);
    H_MC1 = [R_MC1, T_MC1; 0 0 0 1];
    R_Scaph = A(10:12,1:3);
    T_Scaph = A(10:12,4);
    H_Scaph = [R_Scaph, T_Scaph; 0 0 0 1];
    R_Trap = A(7:9,1:3);
    T_Trap = A(7:9,4);
    H_Trap = [R_Trap, T_Trap; 0 0 0 1];
    
    ERROR_Rad = A(13,1);
    ERROR_MC1 = A(14,1);
    ERROR_Trap = A(15,1);
    ERROR_Scaph = A(16,1);    
else
    % Loads data
       fid = fopen(file,'w+');
        % Loads all the STL files
    [F_MC1_P1, V_MC1_P1] =  STL_ReadFile([dir MC1_P1],true);
    [F_Rad_P1, V_Rad_P1] =  STL_ReadFile([dir Rad_P1],true);
    [F_Scaph_P1, V_Scaph_P1] =  STL_ReadFile([dir Scaph_P1],true);
    [F_Trap_P1, V_Trap_P1] =  STL_ReadFile([dir Trap_P1],true);
    [F_MC1_P2, V_MC1_P2] =  STL_ReadFile([dir MC1_P2]',true);
    [F_Rad_P2, V_Rad_P2] =  STL_ReadFile([dir Rad_P2],true);
    [F_Scaph_P2, V_Scaph_P2] =  STL_ReadFile([dir Scaph_P2],true);    
    [F_Trap_P2, V_Trap_P2] =  STL_ReadFile([dir Trap_P2],true);
    
    % Registers radius in position 2 with the radius in position 1
    [R_Rad, T_Rad, ERROR_Rad] = RegisterBones(F_Rad_P2, V_Rad_P2,F_Rad_P1, V_Rad_P1,1);
    H_Rad = [R_Rad, T_Rad; 0 0 0 1];
    
    % Saves results to file
    W = reshape([R_Rad, T_Rad].',1,[]);
    fprintf(fid,'%5.8f\t%5.8f\t%5.8f\t%5.8f\r\n',W); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Transformes all data to radius CS
V_Rad_P1 = R_glob2loc*V_Rad_P1.' + repmat(T_glob2loc,1,size(V_Rad_P1,1));
V_Rad_P1 = V_Rad_P1.';
MC1_P1 = (glob2loc * [V_MC1_P1.' ; ones(1,size(V_MC1_P1,1))]).';
MC1_P1 = MC1_P1(:,1:3);
Scaph_P1 = (glob2loc * [V_Scaph_P1.' ; ones(1,size(V_Scaph_P1,1))]).';
Scaph_P1 = Scaph_P1(:,1:3);
Trap_P1 = (glob2loc * [V_Trap_P1.' ; ones(1,size(V_Trap_P1,1))]).';
Trap_P1 = Trap_P1(:,1:3);

V_Rad_P2 = (glob2loc * H_Rad * [V_Rad_P2.' ; ones(1,size(V_Rad_P2,1))]).';
V_Rad_P2 = V_Rad_P2(:,1:3);
MC1_P2 = (glob2loc * H_Rad * [V_MC1_P2.' ; ones(1,size(V_MC1_P2,1))]).';
MC1_P2 = MC1_P2(:,1:3);
Scaph_P2 = (glob2loc * H_Rad * [V_Scaph_P2.' ; ones(1,size(V_Scaph_P2,1))]).';
Scaph_P2 = Scaph_P2(:,1:3);
Trap_P2 = (glob2loc * H_Rad * [V_Trap_P2.' ; ones(1,size(V_Trap_P2,1))]).';
Trap_P2 = Trap_P2(:,1:3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculates centroids motion
% Scaphoid
cent_Scaph_P1 = mean(Scaph_P1);
cent_Scaph_P2 = mean(Scaph_P2);
% Trapezium
cent_Trap_P1 = mean(Trap_P1);
cent_Trap_P2 = mean(Trap_P2);

% Calculates centroids displacement (not realy necessary)
disp_Trap_P2P1 = cent_Trap_P1 - cent_Trap_P2;
disp_Scaph_P2P1 = cent_Scaph_P1 - cent_Scaph_P2;
dist_TrapScaph_P1 = cent_Scaph_P1 - cent_Trap_P1;
dist_TrapScaph_P2 = cent_Scaph_P2 - cent_Trap_P2;

% Calculates principle axes of inertia

if(inertial_axis_available)
    % Loads data from file
    I = txt2mat(fileInertialAX);
    R_b_Trap = I(1:3,1:3).';
    cent_Trap_P1 = I(1:3,4).';
    R_b_Sca = I(4:6,1:3).';
    cent_Scaph_P1 = I(4:6,4).';
    Trap2Rad = [R_b_Trap.' cent_Trap_P1'; 0 0 0 1];
    Rad2Trap = [R_b_Trap -R_b_Trap*cent_Trap_P1'; 0 0 0 1];
    Sca2Rad = [R_b_Sca.' cent_Scaph_P1'; 0 0 0 1];
    Rad2Sca = [R_b_Sca -R_b_Sca*cent_Scaph_P1'; 0 0 0 1];
else
    % Calculates inertia axes
    [weights_Trap, normals_Trap] = TRI_Areas(F_Trap_P1, Trap_P1, true);
    [~, R_b_Trap, loc_Trap] = KIN_DirectedSpringRotation(TRI_Centroids(F_Trap_P1, Trap_P1), normals_Trap, weights_Trap);
    Trap2Rad = [R_b_Trap.' cent_Trap_P1'; 0 0 0 1];
    Rad2Trap = [R_b_Trap -R_b_Trap*cent_Trap_P1'; 0 0 0 1];

    [weights_Sca, normals_Sca] = TRI_Areas(F_Scaph_P1, Scaph_P1, true);
    [~, R_b_Sca, loc_Sca] = KIN_DirectedSpringRotation(TRI_Centroids(F_Scaph_P1, Scaph_P1), normals_Sca, weights_Sca);
    Sca2Rad = [R_b_Sca.' cent_Scaph_P1'; 0 0 0 1];
    Rad2Sca = [R_b_Sca -R_b_Sca*cent_Scaph_P1'; 0 0 0 1];
    
        % Saves to file
        fid2 = fopen(fileInertialAX,'w+');
        W = reshape([R_b_Trap.' cent_Trap_P1'].',1,[]);
        fprintf(fid2,'%5.8f\t%5.8f\t%5.8f\t%5.8f\r\n',W); 
        W = reshape([R_b_Sca.' cent_Scaph_P1'].',1,[]);
        fprintf(fid2,'%5.8f\t%5.8f\t%5.8f\t%5.8f\r\n',W);
        fclose(fid2);
end
    % Plots bone fragments and inertia axes
    h123 = figure;
    [obj, li, ax] = GUI_PlotShells(h123, {F_Trap_P1,F_Scaph_P1},...
        {Trap_P1,Scaph_P1},{ones(size(Trap_P1,1),1),ones(size(Scaph_P1,1),1)});
    box off
    hold on
    arrow(cent_Trap_P1,cent_Trap_P1+20*R_b_Trap(1,:),5,70,30,'EdgeColor','r','FaceColor','r');
    arrow(cent_Trap_P1,cent_Trap_P1+20*R_b_Trap(2,:),5,70,30,'EdgeColor','g','FaceColor','g');
    arrow(cent_Trap_P1,cent_Trap_P1+20*R_b_Trap(3,:),5,70,30,'EdgeColor','b','FaceColor','b');
    arrow(cent_Scaph_P1,cent_Scaph_P1+20*R_b_Sca(1,:),5,70,30,'EdgeColor','r','FaceColor','r');
    arrow(cent_Scaph_P1,cent_Scaph_P1+20*R_b_Sca(2,:),5,70,30,'EdgeColor','g','FaceColor','g');
    arrow(cent_Scaph_P1,cent_Scaph_P1+20*R_b_Sca(3,:),5,70,30,'EdgeColor','b','FaceColor','b');
    hold off

% Transformes data to local CS. Once transformed the ICP starts.

% Scaphoid expressed in the radius CS
if(~data_available)
    [R_Scaph, T_Scaph, ERROR_Scaph] = RegisterBones(F_Scaph_P2, Scaph_P2, F_Scaph_P1, Scaph_P1,0);
end
Hom_Scaph = [R_Scaph, T_Scaph; 0 0 0 1];
inv_Hom_Sca = [R_Scaph.' -R_Scaph.'*T_Scaph; 0 0 0 1];

% Trapezium expressed in the scaphoid CS
Trap_P1_loc = (Rad2Sca * [Trap_P1.' ; ones(1,size(Trap_P1,1))]).';
Trap_P1_loc = -Trap_P1_loc(:,1:3);
Trap_P2_loc = (Rad2Sca * Hom_Scaph * [Trap_P2.' ; ones(1,size(Trap_P2,1))]).';
Trap_P2_loc = -Trap_P2_loc(:,1:3);
if(~data_available)
    [R_Trap, T_Trap, ERROR_Trap] = RegisterBones(F_Trap_P2, Trap_P2_loc, F_Trap_P1, Trap_P1_loc,0);
end
Hom_Trap = [R_Trap, T_Trap; 0 0 0 1];
inv_Hom_Trap = [R_Trap.' -R_Trap.'*T_Trap; 0 0 0 1];

% MC1 expressed in the trapezium CS.
% (the transformations are necessary to calculate the motion of MC1
% relative to the trapezium)
MC1_P1_loc = (Rad2Trap * [MC1_P1.' ; ones(1,size(MC1_P1,1))]).';
MC1_P1_loc = -MC1_P1_loc(:,1:3);
MC1_P2_loc = (Rad2Trap * Sca2Rad * Hom_Trap * Rad2Sca * Hom_Scaph *[MC1_P2.' ; ones(1,size(MC1_P2,1))]).';
MC1_P2_loc = -MC1_P2_loc(:,1:3);
if(~data_available)
    [R_MC1, T_MC1, ERROR_MC1] = RegisterBones(F_MC1_P2, MC1_P2_loc, F_MC1_P1, MC1_P1_loc,1);
end
Hom_MC1 = [R_MC1, T_MC1; 0 0 0 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Saving the gathered info (so you don't have to re-run the ICP)
if(~data_available)
    % Transformation of MC1
    W = reshape([R_MC1, T_MC1].',1,[]);
    fprintf(fid,'%5.8f\t%5.8f\t%5.8f\t%5.8f\r\n',W);
    % Transformation of Trapezium
    W = reshape([R_Trap, T_Trap].',1,[]);
    fprintf(fid,'%5.8f\t%5.8f\t%5.8f\t%5.8f\r\n',W);
    % Transformation of Scaphoid
    W = reshape([R_Scaph, T_Scaph].',1,[]);
    fprintf(fid,'%5.8f\t%5.8f\t%5.8f\t%5.8f\r\n',W);
    fprintf(fid,'%5.8f\r\n%5.8f\r\n%5.8f\r\n%5.8f\r\n',ERROR_Rad,ERROR_MC1,ERROR_Trap,ERROR_Scaph);
    fclose(fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Checks and tests

% Plots each motion relative to the radius CS.

    h8 = figure;
    [obj, li, ax] = GUI_PlotShells(h8, {F_Rad_P1;F_Rad_P2;F_Trap_P2;F_Trap_P1;F_MC1_P2;F_MC1_P1;F_Scaph_P2;F_Scaph_P1},...
        {V_Rad_P1;V_Rad_P2;Trap_P2;Trap_P1;MC1_P2;MC1_P1;Scaph_P2;Scaph_P1},...
        {ones(size(V_Rad_P1,1),1),ones(size(V_Rad_P2,1),1),ones(size(Trap_P2,1),1),ones(size(Trap_P1,1),1),ones(size(MC1_P2,1),1),ones(size(MC1_P1,1),1),ones(size(Scaph_P2,1),1),ones(size(Scaph_P1,1),1)});
    hold on
    arrow([0;0;0],[10;0;0],5,70,30);
    arrow([0;0;0],[0;10;0],5,70,30);
    arrow([0;0;0],[0;0;10],5,70,30);
    hold off;

% test: transform the first segments so it should match-up with the
% second one.
    h1 = figure;
    [obj, li, ax] = GUI_PlotShells(h1, {F_Rad_P1;F_Rad_P2}, {V_Rad_P1;V_Rad_P2},...
            {ones(size(V_Rad_P1,1),1),ones(size(V_Rad_P2,1),1)});
    box off
    
    Scaph_P2 = (Hom_Scaph * [Scaph_P2.' ; ones(1,size(Scaph_P2,1))]).';
    Scaph_P2 = Scaph_P2(:,1:3);
    h2 = figure;
    [obj, li, ax] = GUI_PlotShells(h2, {F_Scaph_P2;F_Scaph_P1}, {Scaph_P2;Scaph_P1},...
            {ones(size(Scaph_P2,1),1),ones(size(Scaph_P1,1),1)});
    box off
        
    Trap1_loc2 = (Hom_Trap * [Trap_P2_loc.' ; ones(1,size(Trap_P2_loc,1))]).';
    Trap1_loc2 = Trap1_loc2(:,1:3);
    [h3] = figure;
    [obj, li, ax] = GUI_PlotShells(h3, {F_Trap_P2;F_Trap_P1}, {Trap1_loc2;Trap_P1_loc},...
            {ones(size(Trap_P2,1),1),ones(size(Trap_P1,1),1)});
    box off
        
    MC1_loc2 = (Hom_MC1 * [MC1_P2_loc.' ; ones(1,size(MC1_P2_loc,1))]).';
    MC1_loc2 = MC1_loc2(:,1:3);
    h4 = figure;
    [obj, li, ax] = GUI_PlotShells(h4, {F_MC1_P2;F_MC1_P1}, {MC1_loc2;MC1_P1_loc},...
            {ones(size(MC1_P2,1),1),ones(size(MC1_P1,1),1)});
    box off  
    
%   Calculates the helical axes from position 1 to position 2
    % MC1 (plus relative to the radius)
[n_MC1,point_MC1,phi_MC1,t_MC1] = screw(Hom_MC1);
n_MC1_tmp = (Trap2Rad * [n_MC1 ; 1]); n_MC1 = n_MC1_tmp(1:3)/norm(n_MC1_tmp(1:3));
point_MC1_tmp = (Trap2Rad * [point_MC1 ; 1]); point_MC1 = point_MC1_tmp(1:3);
Hel_MC1_p = point_MC1.';
Hel_MC1_dir = n_MC1.';

    % Trapezium (plus relative to the radius)
[n_Trap,point_Trap,phi_Trap,t_Trap] = screw(Hom_Trap);
n_Trap_tmp = (Sca2Rad * [n_Trap ; 1]); n_Trap = n_Trap_tmp(1:3)/norm(n_Trap_tmp(1:3));
point_Trap_tmp = (Sca2Rad * [point_Trap ; 1]); point_Trap = point_Trap_tmp(1:3);
Hel_Trap_p = point_Trap.';
Hel_Trap_dir = n_Trap.';

    % Scaphoid
[n_Scaph,point_Scaph,phi_Scaph,t_Scaph] = screw(Hom_Scaph);
Hel_Scaph_p = point_Scaph.';
Hel_Scaph_dir = n_Scaph.';

% Plots STl files and the helical axes
    h5 = figure;
    [obj, li, ax] = GUI_PlotShells(h5, {F_Rad_P1;F_Trap_P1;F_MC1_P1;F_Scaph_P1},...
        {V_Rad_P1;Trap_P1;MC1_P1;Scaph_P1},...
        {ones(size(V_Rad_P1,1),1),ones(size(Trap_P1,1),1),ones(size(MC1_P1,1),1),ones(size(Scaph_P1,1),1)});
        
        % Plots the coordinate system and its origin
    hold on
    arrow([0;0;0],[10;0;0],5,70,30);
    arrow([0;0;0],[0;10;0],5,70,30);
    arrow([0;0;0],[0;0;10],5,70,30);
        
        % Plots the helical axis for each joint
    arrow(point_MC1-200*n_MC1,point_MC1+200*n_MC1,5,70,30,'EdgeColor','b','FaceColor','b','LineWidth',4);
    arrow(point_Trap-200*n_Trap,point_Trap+200*n_Trap,5,70,30,'EdgeColor','g','FaceColor','g','LineWidth',4);
    arrow(point_Scaph-200*n_Scaph,point_Scaph+200*n_Scaph,5,70,30,'EdgeColor','m','FaceColor','r','LineWidth',4);
        
        % Saves to file
    fid4 = fopen(fileHelicalAX,'w+');
    fprintf(fid4,'%5.8f\t%5.8f\t%5.8f\r\n',Hel_MC1_p,Hel_MC1_dir,Hel_Trap_p,Hel_Trap_dir,...
        Hel_Scaph_p,Hel_Scaph_dir);
    
% Definition of the horizontal plane defined on the radius
    % Equation of a plane: ax + by + cz + d = 0 (P)
        % To find the coefficient a, b and c we need 3 points
        % Here, we can use the points selected to define the radius
        % coordinate system, respectively called O, Y and z (see above)
        % We can find the equation of the plane by solving this equation:
        % OP.(OY x YZ) = 0     (with P(x,y,z))
        
        O = [0 0 0];
        Yaxis = [0 1 0];
        Zaxis = [0 0 1];
        OY = Yaxis - O;
        YZ = Zaxis - Yaxis;
        OYxYZ = [OY(1,2)*YZ(1,3)-YZ(1,2)*OY(1,3);...
            -(OY(1,1)*YZ(1,3)-YZ(1,1)*OY(1,3));...
            OY(1,1)*YZ(1,2)-YZ(1,1)*OY(1,2)];
        P = O;
        a = OYxYZ(1,1);
        b = OYxYZ(2,1);
        c = OYxYZ(3,1);
        d = a*P(1,1) + b*P(1,2) + c*P(1,3);

% Intersection between each helical axis and the plane OYZ (on the radius)
    % First, we define the parametric equation of each helical axis
        % xi = ai.t + pxi
        % yi = bi.t + pyi
        % zi = ci.t + pzi
        
        % With:
        % ai: n_MC(1,i) | bi: n_MC(2,i) | ci: n_MC(3,i)
        % pxi: point_MC(1,i) | pyi: point_MC(2,i) | pzi: point_MC(3,i)
        
        % To find the intersection between the plane (P) and the helical
        % axis, we need to solve: P(x(t),y(t),z(t)) = 0, then replace t in
        % the parametric equation to have the (x,y,z) coordinates of the
        % intersection point

% Scaphoid
    ai_Scaph = n_Scaph(1);
    bi_Scaph = n_Scaph(2);
    ci_Scaph = n_Scaph(3);
    px_Scaph = point_Scaph(1);
    py_Scaph = point_Scaph(2);
    pz_Scaph = point_Scaph(3);
    syms t
    solt_Scaph = double(solve(a*((ai_Scaph)*t+px_Scaph)+...
        b*((bi_Scaph)*t+py_Scaph)+c*((ci_Scaph)*t+pz_Scaph)+d==0,t));
    Px_Scaph = ai_Scaph*solt_Scaph+px_Scaph;
    Py_Scaph = bi_Scaph*solt_Scaph+py_Scaph;
    Pz_Scaph = ci_Scaph*solt_Scaph+pz_Scaph;

fig_Scaph = figure;
Inter_Scaph = scatter(Py_Scaph,Pz_Scaph,'fill','m');
grid on
yL_Scaph = xlim;
zL_Scaph = ylim;
hold on
plot(yL_Scaph, [0 0], 'k-') % Show Y axis = 0 line
plot([0 0], zL_Scaph, 'k-') % Show Z axis = 0 line
xlabel('Y axis');
ylabel('Z axis')
title('Scaphoid: Intersection of the helical axis with the YZ plane of the radius CS');

% Trapezium
    ai_Trap = n_Trap(1);
    bi_Trap = n_Trap(2);
    ci_Trap = n_Trap(3);
    px_Trap = point_Trap(1);
    py_Trap = point_Trap(2);
    pz_Trap = point_Trap(3);
    syms t
    solt_Trap = double(solve(a*((ai_Trap)*t+px_Trap)+...
        b*((bi_Trap)*t+py_Trap)+c*((ci_Trap)*t+pz_Trap)+d==0,t));
    Px_Trap = ai_Trap*solt_Trap+px_Trap;
    Py_Trap = bi_Trap*solt_Trap+py_Trap;
    Pz_Trap = ci_Trap*solt_Trap+pz_Trap;

firg_Trap = figure;
Inter_Trap = scatter(Py_Trap,Pz_Trap,'fill','g');
grid on
yL_Trap = xlim;
zL_Trap = ylim;
hold on
plot(yL_Trap, [0 0], 'k-') % Show Y axis = 0 line
plot([0 0], zL_Trap, 'k-') % Show Z axis = 0 line
xlabel('Y axis');
ylabel('Z axis')
title('Trapezium: Intersection of the helical axis with the YZ plane of the radius CS');

% MC1
    ai_MC1 = n_MC1(1);
    bi_MC1 = n_MC1(2);
    ci_MC1 = n_MC1(3);
    px_MC1 = point_MC1(1);
    py_MC1 = point_MC1(2);
    pz_MC1 = point_MC1(3);
    syms t
    solt_MC1 = double(solve(a*((ai_MC1)*t+px_MC1)+...
        b*((bi_MC1)*t+py_MC1)+c*((ci_MC1)*t+pz_MC1)+d==0,t));
    Px_MC1 = ai_MC1*solt_MC1+px_MC1;
    Py_MC1 = bi_MC1*solt_MC1+py_MC1;
    Pz_MC1 = ci_MC1*solt_MC1+pz_MC1;

fig_MC1 = figure;
Inter_MC1 = scatter(Py_MC1,Pz_MC1,'fill','b');
grid on
yL_MC1 = xlim;
zL_MC1 = ylim;
hold on
plot(yL_MC1, [0 0], 'k-') % Show Y axis = 0 line
plot([0 0], zL_MC1, 'k-') % Show Z axis = 0 line
xlabel('Y axis');
ylabel('Z axis')
title('MC1: Intersection of the helical axis with the YZ plane of the radius CS');

% Calculation of the angle between the helical axis and the YZ plane of the
% radius CS:
    % We can calculate this angle from the coordinates of the vector normal
    % to the plane (here: the x-axis = (1,0,0)) and the vector defining the 
    % direction of the helical axis (n_bone)

% Scaphoid
alpha_Scaph = asin(abs(n_Scaph(1))/(sqrt(n_Scaph(1)^2+n_Scaph(2)^2+n_Scaph(3)^2)))*180/pi

% Trapezium
alpha_Trap = asin(abs(n_Trap(1))/(sqrt(n_Trap(1)^2+n_Trap(2)^2+n_Trap(3)^2)))*180/pi

% MC1
alpha_MC1 = asin(abs(n_MC1(1))/(sqrt(n_MC1(1)^2+n_MC1(2)^2+n_MC1(3)^2)))*180/pi

% Gather all graphs together

Final_Graph = figure;
subplot(1,3,1);
Inter_Scaph = scatter(Py_Scaph,Pz_Scaph,'fill','r');
grid on
hold on
plot(yL_Scaph, [0 0], 'k-') % Show Y axis = 0 line
plot([0 0], zL_Scaph, 'k-') % Show Z axis = 0 line
xlabel('Y axis');
ylabel('Z axis')
title('Scaphoid: Intersection of the helical axis with the YZ plane of the radius CS');

subplot(1,3,2);
Inter_Trap = scatter(Py_Trap,Pz_Trap,'fill','g');
grid on
hold on
plot(yL_Trap, [0 0], 'k-') % Show Y axis = 0 line
plot([0 0], zL_Trap, 'k-') % Show Z axis = 0 line
xlabel('Y axis');
ylabel('Z axis')
title('Trapezium: Intersection of the helical axis with the YZ plane of the radius CS');

subplot(1,3,3);
Inter_MC1 = scatter(Py_MC1,Pz_MC1,'fill','b');
grid on
hold on
plot(yL_MC1, [0 0], 'k-') % Show Y axis = 0 line
plot([0 0], zL_MC1, 'k-') % Show Z axis = 0 line
xlabel('Y axis');
ylabel('Z axis')
title('MC1: Intersection of the helical axis with the YZ plane of the radius CS');

% Calculation of the rotation angles along each axis of the radius
% coordinate system for each bone. Each rotation angle is expressed in the
% radius coordinate system relatively to the more proximal bone.

    % MC1
    
    Rot_MC1_x = asin(R_MC1(3,2)/cos(asin(-R_MC1(3,1))))*180/pi;
    Rot_MC1_y = asin(-R_MC1(3,1))*180/pi;
    Rot_MC1_z = asin(R_MC1(2,1)/cos(asin(-R_MC1(3,1))))*180/pi;
    Rot_MC1 = [Rot_MC1_x,Rot_MC1_y,Rot_MC1_z];
    
    % Trapezium
    
    Rot_Trap_x = asin(R_Trap(3,2)/cos(asin(-R_Trap(3,1))))*180/pi;
    Rot_Trap_y = asin(-R_Trap(3,1))*180/pi;
    Rot_Trap_z = asin(R_Trap(2,1)/cos(asin(-R_Trap(3,1))))*180/pi;
    Rot_Trap = [Rot_Trap_x,Rot_Trap_y,Rot_Trap_z];
    
    % Scaph
    
    Rot_Scaph_x = asin(R_Scaph(3,2)/cos(asin(-R_Scaph(3,1))))*180/pi;
    Rot_Scaph_y = asin(-R_Scaph(3,1))*180/pi;
    Rot_Scaph_z = asin(R_Scaph(2,1)/cos(asin(-R_Scaph(3,1))))*180/pi;
    Rot_Scaph = [Rot_Scaph_x,Rot_Scaph_y,Rot_Scaph_z];
    
    % Saves to file
        fid5 = fopen(fileRotAngle,'w+');
        fprintf(fid5,'%5.8f\t%5.8f\t%5.8f\r\n',Rot_MC1,Rot_Trap,Rot_Scaph);