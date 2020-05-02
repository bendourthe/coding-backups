% Selection of the current directory (where the STL files are)
dir = 'D:\PhD KU Leuven\Data\Study_Brown_Rig\New Segmentation\Scan20\STL\';

% FILES
Bone1 = 'N_mc1_s_0.4_1.stl';
Bone2 = 'N_trap_s_0.4_1.stl';

% Read STL files
[F_Bone1, V_Bone1] =  STL_ReadFile([dir Bone1],true);
[F_Bone2, V_Bone2] =  STL_ReadFile([dir Bone2],true);

% Calculation of the shortest distance between the two STL selected
[idx1,D1]=knnsearch(V_Bone1,V_Bone2);
[C1,idxD1]=min(D1);
Shortest_Distance=C1

% Plot the two STL selected with an arrow joining the two points
%       calculated as the closest neighbors

[idx2,D2]=knnsearch(V_Bone2,V_Bone1);
[C2,idxD2]=min(D2);

P1=V_Bone1(idx1(idxD1),:);
P2=V_Bone2(idx2(idxD2),:);
SD=((P2(1)-P1(1))^2+(P2(2)-P1(2))^2+(P2(3)-P1(3))^2)^(1/2);

stl = figure;
    [obj, li, ax] = GUI_PlotShells(stl, {F_Bone2}, {V_Bone2},...
            {ones(size(V_Bone2,1),1)},[0.4375 0.75 0.9375]);
hold on