%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shortest_Distance_STL_MC1_Trap: Calculates the shortest distance between
%                   the first metacarpal and the trapezium for each frame
%                   of a 4D CT dataset and plot those results as a frame
%                   function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dependencies:        
%               STL_ReadFile.m 
%                   TRI_RemoveInvalidTriangles.m
%                   TRI_RemoveBadlyConnectedTriangles.m
% Input: 
%           MC1_X: STL file corresponding to the Xth frame of the MC1
%           Trap_X: STL file corresponding to the Xth frame of the Trap
% Output:
%           Evolution of the Shortest Distance between the first metacarpal
%           and the trapezium for each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Selection of the current directory (where the STL files are)
dir = 'C:\Users\u0096636\Documents\Professional\PhD KU Leuven\Data\Static scans\STL\';

% FILES

% MC1 files:
MC1_1 = 'MC1_Extension.stl';
MC1_2 = 'MC1_Neutral.stl';
MC1_3 = 'Mc1_Flexion.stl';

% Trap files:
Trap_1 = 'Trap_Extension.stl';
Trap_2 = 'Trap_Neutral.stl';
Trap_3 = 'Trap_Flexion.stl';

% Read STL files

% MC1 files:
[F_MC1_1, V_MC1_1] =  STL_ReadFile([dir MC1_1],true);
[F_MC1_2, V_MC1_2] =  STL_ReadFile([dir MC1_2],true);
[F_MC1_3, V_MC1_3] =  STL_ReadFile([dir MC1_3],true);

% Trap files:
[F_Trap_1, V_Trap_1] =  STL_ReadFile([dir Trap_1],true);
[F_Trap_2, V_Trap_2] =  STL_ReadFile([dir Trap_2],true);
[F_Trap_3, V_Trap_3] =  STL_ReadFile([dir Trap_3],true);

% Calculation of the shortest distance between MC1 and Trap for each frame

[idx1,D1]=knnsearch(V_MC1_1,V_Trap_1);
SD1=min(D1);
[idx2,D2]=knnsearch(V_MC1_2,V_Trap_2);
SD2=min(D2);
[idx3,D3]=knnsearch(V_MC1_3,V_Trap_3);
SD3=min(D3);


SD=[SD1;SD2;SD3];
plot(SD);
title('Evolution of the shortest distance between MC1 and Trap for each frame');
xlabel('frames');
ylabel('Shortest Distance (mm)');


