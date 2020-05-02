clear all; close all;
% Function: Define the trapezium coordinate system (CS), perform a 4-points
% registration based on an ICP algorithm to calculate the transformation
% matrices from a position 1 to a position 2, calculates the helical axis
% of the joint and plot it.
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
% Input: STl files of the trapezium and first metacarpal(MC1).
%
% Output: txt-files with transformations and coordinate systems. Figures 
% showing the efficiency of the registration and the helical axis.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Controle panel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Folder with files
dir = 'H:\Data\Study_Brown_Rig\Segmentation\Scan20\STL\'; 
data_available = 0; % loads data if available (=1), else, generate data
    filename = 'RotTransData_N_P80.txt';
    file = [dir filename];
trap_cs_available = 0;          % loads data if available (=1), else, generate data
    filenameCS_trap = 'TrapAxes_N_P80.txt';
    fileCS_trap = [dir filenameCS_trap];
helical_axis_available = 0; % loads data if available (=1), else, generate data
    filenameHelicalAX = 'HelicalAx_N_P80.txt';
    fileHelicalAX = [dir filenameHelicalAX];
rot_angle_available = 0; % loads data if available (=1), else, generate data
    filenameRotAngle = 'RotAngle_N_P80.txt';
    fileRotAngle = [dir filenameRotAngle];

% FILES
    % Note: the code works so that P2 is registered on P1, so if you want
    % to study the motion Ext-Flex (for instance), you should put Ext as P2
    % and Flex as P1 (and vice-versa if you want to study Flex-Ext motion)
% First position (referenced as P1)
MC1_P1 = 'SCAN20_P80_mc1.stl';
Trap_P1 = 'SCAN20_P80_trap.stl';

% Second position (referenced as P2)
MC1_P2 = 'SCAN20_N_mc1.stl';
Trap_P2 = 'SCAN20_N_trap.stl';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CS on trapezium:
%   origin: midway between the two landmarks defining the y-axis.
%   y: extends from the exact mid-point of the central ridge of the 
%       trapezial saddle to the center of the junction of the trapezium,
%       scaphoid and trapezoid.
%   x: runs in a dorsal-to-volar direction along a line perpendicular to
%       the central ridge of the trapezium and passes through the mid-point
%       of the dorsal surface to the proximal volar pole of the tubercle of
%       the trapezium.
%   z: is perpendicular to the X and Y axes and nearly parallel to the
%       central ridge of the trapezial metacarpal surface.

if(trap_cs_available)
    A = txt2mat(fileCS_trap);
    
    X1 = A(1,:);    % mid-point of the central ridge of the trapezial saddle
    X2 = A(2,:);    % center of the junction of the trapezium, scaphoid and trapezoid
    Y = A(3,:);     % center of the distal articular surface, facing the center of the base of MC1
else
    [F_Trap_P1, V_Trap_P1] =  STL_ReadFile([dir Trap_P1],true);
    [ad_curve] = PlacePoints3({F_Trap_P1}, {V_Trap_P1}, {F_Trap_P1}, {V_Trap_P1}, 'Select 3 landmarks (X1, X2 and Y, in that order!). End by pressing enter and exiting the figure');
    close all;
    % Retrieves coordinates of the points:
    X1 = ad_curve.points{end,1}(1,:);
    X2 = ad_curve.points{end,1}(2,:);
    Y = ad_curve.points{end,1}(3,:);
    
    % Save in a txt file
    fid = fopen(fileCS_trap,'w+');
    fprintf(fid,'%5.8f\t%5.8f\t%5.8f\r\n',[X1,X2,Y]);
    fclose(fid);
end

O_trap = (X1+X2)/2;
Xaxis_trap = (X2-X1)/norm(X2-X1);
Yaxis_trap = (O_trap- Y)/norm(O_trap-Y);
Zaxis_trap = cross(Xaxis_trap,Yaxis_trap)/norm(cross(Xaxis_trap,Yaxis_trap));
Yaxis_trap = cross(Zaxis_trap,Xaxis_trap)/norm(cross(Zaxis_trap,Xaxis_trap));
cs_trap = [Xaxis_trap; Yaxis_trap; Zaxis_trap];
% Converges to trapezium CS
RotMat_trap = [Xaxis_trap', Yaxis_trap', Zaxis_trap'];
loc2glob_trap = [RotMat_trap, O_trap';
     0 0 0 1]; 
   % These tranformation matrices convert from local to global CS.
   % Invert them to go from global to local.
glob2loc = [RotMat_trap', -(RotMat_trap')*O_trap'; 0 0 0 1];
R_glob2loc = glob2loc(1:3,1:3);
T_glob2loc = glob2loc(1:3,4);

if(trap_cs_available)
else
% Transformes data to trapezium CS
V_Trap_P1 = R_glob2loc*V_Trap_P1.' + repmat(T_glob2loc,1,size(V_Trap_P1,1));
V_Trap_P1 = V_Trap_P1.';

% Plots the trapezium with the new coordinate system for validation
    h0 = figure;
    [obj, li, ax] = GUI_PlotShells(h0, {F_Trap_P1},...
        {V_Trap_P1},...
        {ones(size(V_Trap_P1,1),1)});
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

% Plots the trapezium with the new coordinate system
    h0 = figure;
    [obj, li, ax] = GUI_PlotShells(h0, {F_Trap_P1},...
        {V_Trap_P1},...
        {ones(size(V_Trap_P1,1),1)});
    GUI_VisualisationUI(ancestor(ax, 'figure'), true, ax, true, true);
        
        % Plots the coordinate system and its origin
    hold on
    arrow([0;0;0],[30;0;0],10,70,30,5,'EdgeColor','r','FaceColor','r');
    arrow([0;0;0],[0;30;0],10,70,30,5,'EdgeColor','g','FaceColor','g');
    arrow([0;0;0],[0;0;30],10,70,30,5,'EdgeColor','b','FaceColor','b');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Searches for trapezium fragment
if(data_available)
    % Loads data
    A = txt2mat(file);
        % Loads all the STL files
    [F_MC1_P1, V_MC1_P1] =  STL_ReadFile([dir MC1_P1],true);
    [F_Trap_P1, V_Trap_P1] =  STL_ReadFile([dir Trap_P1],true);
    [F_MC1_P2, V_MC1_P2] =  STL_ReadFile([dir MC1_P2]',true);
    [F_Trap_P2, V_Trap_P2] =  STL_ReadFile([dir Trap_P2],true);

        % Associates the different parts of the RotTransData file to the
        % corresponding rotation and translation matrices for each bone
    R_Trap = A(1:3,1:3);
    T_Trap = A(1:3,4);
    H_Trap = [R_Trap, T_Trap; 0 0 0 1];
    R_MC1 = A(4:6,1:3);
    T_MC1 = A(4:6,4);
    H_MC1 = [R_MC1, T_MC1; 0 0 0 1];
    
else
    % Loads data
       fid = fopen(file,'w+');
        % Loads all the STL files
    [F_MC1_P1, V_MC1_P1] =  STL_ReadFile([dir MC1_P1],true);
    [F_Trap_P1, V_Trap_P1] =  STL_ReadFile([dir Trap_P1],true);
    [F_MC1_P2, V_MC1_P2] =  STL_ReadFile([dir MC1_P2]',true);
    [F_Trap_P2, V_Trap_P2] =  STL_ReadFile([dir Trap_P2],true);
    
    % Registers trapezium in position 2 with the trapezium in position 1
    [R_Trap, T_Trap, ERROR_Trap] = RegisterBones(F_Trap_P2, V_Trap_P2, F_Trap_P1, V_Trap_P1,1);
    H_Trap = [R_Trap, T_Trap; 0 0 0 1];
    
    % Saves results to file
    W = reshape([R_Trap, T_Trap].',1,[]);
    fprintf(fid,'%5.8f\t%5.8f\t%5.8f\t%5.8f\r\n',W); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Transformes all data to trapezium CS
V_Trap_P1 = R_glob2loc*V_Trap_P1.' + repmat(T_glob2loc,1,size(V_Trap_P1,1));
V_Trap_P1 = V_Trap_P1.';
MC1_P1 = (glob2loc * [V_MC1_P1.' ; ones(1,size(V_MC1_P1,1))]).';
MC1_P1 = MC1_P1(:,1:3);

V_Trap_P2 = (glob2loc * H_Trap * [V_Trap_P2.' ; ones(1,size(V_Trap_P2,1))]).';
V_Trap_P2 = V_Trap_P2(:,1:3);
MC1_P2 = (glob2loc * H_Trap * [V_MC1_P2.' ; ones(1,size(V_MC1_P2,1))]).';
MC1_P2 = MC1_P2(:,1:3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Transformes data to local CS. Once transformed the ICP starts.

% MC1 expressed in the trapezium CS.
if(~data_available)
    [R_MC1, T_MC1, ERROR_MC1] = RegisterBones(F_MC1_P2, MC1_P2, F_MC1_P1, MC1_P1,1);
end
Hom_MC1 = [R_MC1, T_MC1; 0 0 0 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Saving the gathered info (so you don't have to re-run the ICP)
if(~data_available)
    % Transformation of MC1
    W = reshape([R_MC1, T_MC1].',1,[]);
    fprintf(fid,'%5.8f\t%5.8f\t%5.8f\t%5.8f\r\n',W);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Checks and tests

% Plots each motion relative to the trapezium CS.

    h8 = figure;
    [obj, li, ax] = GUI_PlotShells(h8, {F_Trap_P1;F_Trap_P2;F_MC1_P2;F_MC1_P1},...
        {V_Trap_P1;V_Trap_P2;MC1_P2;MC1_P1},...
        {ones(size(V_Trap_P1,1),1),ones(size(V_Trap_P2,1),1),ones(size(MC1_P2,1),1),ones(size(MC1_P1,1),1)});
    hold on
    arrow([0;0;0],[10;0;0],5,70,30);
    arrow([0;0;0],[0;10;0],5,70,30);
    arrow([0;0;0],[0;0;10],5,70,30);
    hold off;

% test: transform the first segments so it should match-up with the
% second one.
    h1 = figure;
    [obj, li, ax] = GUI_PlotShells(h1, {F_Trap_P1;F_Trap_P2}, {V_Trap_P1;V_Trap_P2},...
            {ones(size(V_Trap_P1,1),1),ones(size(V_Trap_P2,1),1)});
    box off
        
    MC1_loc2 = (Hom_MC1 * [MC1_P2.' ; ones(1,size(MC1_P2,1))]).';
    MC1_loc2 = MC1_loc2(:,1:3);
    h4 = figure;
    [obj, li, ax] = GUI_PlotShells(h4, {F_MC1_P2;F_MC1_P1}, {MC1_loc2;MC1_P1},...
            {ones(size(MC1_P2,1),1),ones(size(MC1_P1,1),1)});
    box off  
    
%   Calculates the helical axes from position 1 to position 2
    % MC1
[n_MC1,point_MC1,phi_MC1,t_MC1] = screw(Hom_MC1);
Hel_MC1_p = point_MC1.';
Hel_MC1_dir = n_MC1.';

% Plots STl files and the helical axes
    h5 = figure;
    [obj, li, ax] = GUI_PlotShells(h5, {F_Trap_P1;F_MC1_P1},...
        {V_Trap_P1;MC1_P1},...
        {ones(size(V_Trap_P1,1),1),ones(size(MC1_P1,1),1)});
        
        % Plots the coordinate system and its origin
    hold on
    arrow([0;0;0],[10;0;0],5,70,30);
    arrow([0;0;0],[0;10;0],5,70,30);
    arrow([0;0;0],[0;0;10],5,70,30);
        
        % Plots the helical axis for each joint
    arrow(point_MC1-200*n_MC1,point_MC1+200*n_MC1,5,70,30,'EdgeColor','b','FaceColor','b','LineWidth',4);
        
        % Saves to file
    fid4 = fopen(fileHelicalAX,'w+');
    fprintf(fid4,'%5.8f\t%5.8f\t%5.8f\r\n',Hel_MC1_p,Hel_MC1_dir);
    

% Calculation of the rotation angles along each axis of the radius
% coordinate system for each bone. Each rotation angle is expressed in the
% radius coordinate system relatively to the more proximal bone.

    % MC1
    
    Rot_MC1_x = asin(R_MC1(3,2)/cos(asin(-R_MC1(3,1))))*180/pi;
    Rot_MC1_y = asin(-R_MC1(3,1))*180/pi;
    Rot_MC1_z = asin(R_MC1(2,1)/cos(asin(-R_MC1(3,1))))*180/pi;
    Rot_MC1 = [Rot_MC1_x,Rot_MC1_y,Rot_MC1_z];
    
    % Saves to file
        fid5 = fopen(fileRotAngle,'w+');
        fprintf(fid5,'%5.8f\t%5.8f\t%5.8f\r\n',Rot_MC1);