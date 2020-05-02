%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shortest_Distance_STL: Calculates the shortest distance between two STL                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dependencies:        
%               STL_ReadFile.m 
%                   TRI_RemoveInvalidTriangles.m
%                   TRI_RemoveBadlyConnectedTriangles.m
% Input: 
%           Bone1: STL file corresponding to the first bone
%           Bone2: STL file corresponding to second bone
% 
% Output:
%           Shortest_Distance: Shortest distance calculated between the two
%                               STL files selected (mm)
%           Plotting of the two STL selected with a vector to represent
%                               the shortest distance between those two
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Selection of the current directory (where the STL files are)
dir = 'C:\Users\u0096636\Documents\Professional\PhD KU Leuven\Data\CAD Tests\';

% FILES
Bone1 = 'Cube1.stl';
Bone2 = 'Cube2.stl';

% Read STL files
[F_Bone1, V_Bone1] =  STL_ReadFile([dir Bone1],true);
[F_Bone2, V_Bone2] =  STL_ReadFile([dir Bone2],true);

% Calculation of the shortest distance between the two STL selected
[idx1,D1]=knnsearch(V_Bone1,V_Bone2);
[C1,idxD1]=min(D1);
Shortest_Distance=C1

% Plotting of the two STL selected and the arrow joining the two points
%       calculated as the closest neighbors

[idx2,D2]=knnsearch(V_Bone2,V_Bone1);
[C2,idxD2]=min(D2);

P1=V_Bone1(idx1(idxD1),:);
P2=V_Bone2(idx2(idxD2),:);
distance=((P2(1)-P1(1))^2+(P2(2)-P1(2))^2+(P2(3)-P1(3))^2)^(1/2);        
        
h8 = figure;
    [obj, li, ax] = GUI_PlotShells(h8, {F_Bone1;F_Bone2}, {V_Bone1;V_Bone2},...
            {ones(size(V_Bone1,1),1),ones(size(V_Bone2,1),1)});
hold on
vectarrow(P1,P2);

D=pdist2(V_Bone1,V_Bone2)