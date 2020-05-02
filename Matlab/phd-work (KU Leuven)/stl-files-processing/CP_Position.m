% Selection of the current directory (where the STL files are)
dir = 'C:\Users\u0096636\Documents\Professional\PhD KU Leuven\Data\STL dynamic Faes\WS\';

% FILES
Bone1 = '611L_phase2_mc1_WS.stl';
Bone2 = '611L_phase2_trp_WS.stl';

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

stl = figure;
    [obj, li, ax] = GUI_PlotShells(stl, {F_Bone1}, {V_Bone1},...
            {ones(size(V_Bone1,1),1)});
hold on
CSmax1=plot3(P1(1),P1(2),P1(3),'Marker','o','MarkerSize',10,'Color','k','MarkerFaceColor','k');