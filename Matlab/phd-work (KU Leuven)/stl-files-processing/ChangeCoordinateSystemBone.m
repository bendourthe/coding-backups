clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ChangeCoordinateSystemBone: User can manually select a new coordinate
%                               system for a selected STL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dependencies:        
%               STL_ReadFile.m 
%                   TRI_RemoveInvalidTriangles.m
%                   TRI_RemoveBadlyConnectedTriangles.m
%               vectarrow.m
% Input: 
%           Bone: STL file corresponding to the selected bone
% Output:
%           Plotting of the selecting bone in the new coordinate system       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Selection of the current directory (where the STL files are)
dir = 'D:\PhD KU Leuven\Data\Study_Microarchitecture_Trapezium\raw STLs\';

% FILE
Bone = '4.1_microCT_Declercq Viviane.stl';

% Read STL file
[F_Bone, V_Bone] =  STL_ReadFile([dir Bone],true);

% Plot the STL
stl = figure;
    [obj, li, ax] = GUI_PlotShells(stl, {F_Bone}, {V_Bone},...
            {ones(size(V_Bone,1),1)});
        
% Rotate the STL in order to have the best view for the points selection

rotate3d on
disp('Rotate the figure (5 seconds!)')    % 
pause(5) % Gives 5 seconds to rotate the figure in the desired orientation
rotate3d off

% Select the anatomical landmarks that will be used for the definition of
% the new coordinate system

    % Origin
disp('Select the first point - origin (dorsal beak for MC1)')
k1 = waitforbuttonpress; % k=0 if mouse pressed | k=1 if keyboard
if k1==0
    P = get(gca,'CurrentPoint');
    O = [P(1,1);P(1,2);P(1,3)];
end

    % Y direction
disp('Select the second point - y direction (volar beak for MC1)')
k2 = waitforbuttonpress;
if k2==0
    P = get(gca,'CurrentPoint');
    Y = [P(1,1);P(1,2);P(1,3)];
end

    % X direction
disp('Select the third point - x direction (dorsal beak for MC1) - Then PRESS ENTER')
k3 = waitforbuttonpress;
if k3==0
    P = get(gca,'CurrentPoint');
    X = [P(1,1);P(1,2);P(1,3)];
else
end

% Registration of the new coordinate system

Origin = O;
Yaxis = (Y-O)/norm(Y-O);
Xaxis = (X-O)/norm(X-O);
Zaxis = cross(Xaxis,Yaxis)/norm(cross(Xaxis,Yaxis));
Xaxis = cross(Yaxis,Zaxis)/norm(cross(Yaxis,Zaxis));

X=[Xaxis(1)+O(1);Xaxis(2)+O(2);Xaxis(3)+O(3)];
Y=[Yaxis(1)+O(1);Yaxis(2)+O(2);Yaxis(3)+O(3)];
Z=[Zaxis(1)+O(1);Zaxis(2)+O(2);Zaxis(3)+O(3)];

% Plot the STL with the new coordinate system
k4 = waitforbuttonpress;
if k4==1
    close all
    stl = figure;
    [obj, li, ax] = GUI_PlotShells(stl, {F_Bone}, {V_Bone},...
        {ones(size(V_Bone,1),1)});
    rotate3d on
    hold on
    plot3(O(1),O(2),O(3),'Marker','o','MarkerSize',5,'Color','b','MarkerFaceColor','b');
    hold on
    vectarrow(O,X)
    hold on
    vectarrow(O,Y)
    hold on
    vectarrow(O,Z)
end

% Modification of the V matrix in the new coordinate system

for i=1:1:length(V_Bone)
    V_Bone_new(i,1)=V_Bone(i,1)-O(1);
    V_Bone_new(i,2)=V_Bone(i,2)-O(2);
    V_Bone_new(i,3)=V_Bone(i,3)-O(3);
end