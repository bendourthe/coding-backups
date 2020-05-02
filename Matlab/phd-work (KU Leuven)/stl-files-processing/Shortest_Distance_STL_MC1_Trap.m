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
dir = 'C:\Users\u0096636\Documents\Professional\PhD KU Leuven\Data\STL dynamic Faes\';

% FILES

% MC1 files:
MC1_1 = '611L_phase1_mc1.stl';
MC1_2 = '611L_phase2_mc1.stl';
MC1_3 = '611L_phase3_mc1.stl';
MC1_4 = '611L_phase4_mc1.stl';
MC1_5 = '611L_phase6_mc1.stl';
MC1_6 = '611L_phase9_mc1.stl';
MC1_7 = '611L_phase10_mc1.stl';
MC1_8 = '611L_phase13_mc1.stl';
MC1_9 = '611L_phase16_mc1.stl';
MC1_10 = '611L_phase17_mc1.stl';
MC1_11 = '611L_phase19_mc1.stl';
MC1_12 = '611L_phase22_mc1.stl';
MC1_13 = '611L_phase23_mc1.stl';
MC1_14 = '611L_phase24_mc1.stl';
MC1_15 = '611L_phase26_mc1.stl';
MC1_16 = '611L_phase27_mc1.stl';
MC1_17 = '611L_phase29_mc1.stl';
MC1_18 = '611L_phase33_mc1.stl';
MC1_19 = '611L_phase36_mc1.stl';

% Trap files:
Trap_1 = '611L_phase1_trp.stl';
Trap_2 = '611L_phase2_trp.stl';
Trap_3 = '611L_phase3_trp.stl';
Trap_4 = '611L_phase4_trp.stl';
Trap_5 = '611L_phase6_trp.stl';
Trap_6 = '611L_phase9_trp.stl';
Trap_7 = '611L_phase10_trp.stl';
Trap_8 = '611L_phase13_trp.stl';
Trap_9 = '611L_phase16_trp.stl';
Trap_10 = '611L_phase17_trp.stl';
Trap_11 = '611L_phase19_trp.stl';
Trap_12 = '611L_phase22_trp.stl';
Trap_13 = '611L_phase23_trp.stl';
Trap_14 = '611L_phase24_trp.stl';
Trap_15 = '611L_phase26_trp.stl';
Trap_16 = '611L_phase27_trp.stl';
Trap_17 = '611L_phase29_trp.stl';
Trap_18 = '611L_phase33_trp.stl';
Trap_19 = '611L_phase36_trp.stl';

% Read STL files

% MC1 files:
[F_MC1_1, V_MC1_1] =  STL_ReadFile([dir MC1_1],true);
[F_MC1_2, V_MC1_2] =  STL_ReadFile([dir MC1_2],true);
[F_MC1_3, V_MC1_3] =  STL_ReadFile([dir MC1_3],true);
[F_MC1_4, V_MC1_4] =  STL_ReadFile([dir MC1_4],true);
[F_MC1_5, V_MC1_5] =  STL_ReadFile([dir MC1_5],true);
[F_MC1_6, V_MC1_6] =  STL_ReadFile([dir MC1_6],true);
[F_MC1_7, V_MC1_7] =  STL_ReadFile([dir MC1_7],true);
[F_MC1_8, V_MC1_8] =  STL_ReadFile([dir MC1_8],true);
[F_MC1_9, V_MC1_9] =  STL_ReadFile([dir MC1_9],true);
[F_MC1_10, V_MC1_10] =  STL_ReadFile([dir MC1_10],true);
[F_MC1_11, V_MC1_11] =  STL_ReadFile([dir MC1_11],true);
[F_MC1_12, V_MC1_12] =  STL_ReadFile([dir MC1_12],true);
[F_MC1_13, V_MC1_13] =  STL_ReadFile([dir MC1_13],true);
[F_MC1_14, V_MC1_14] =  STL_ReadFile([dir MC1_14],true);
[F_MC1_15, V_MC1_15] =  STL_ReadFile([dir MC1_15],true);
[F_MC1_16, V_MC1_16] =  STL_ReadFile([dir MC1_16],true);
[F_MC1_17, V_MC1_17] =  STL_ReadFile([dir MC1_17],true);
[F_MC1_18, V_MC1_18] =  STL_ReadFile([dir MC1_18],true);
[F_MC1_19, V_MC1_19] =  STL_ReadFile([dir MC1_19],true);

% Trap files:
[F_Trap_1, V_Trap_1] =  STL_ReadFile([dir Trap_1],true);
[F_Trap_2, V_Trap_2] =  STL_ReadFile([dir Trap_2],true);
[F_Trap_3, V_Trap_3] =  STL_ReadFile([dir Trap_3],true);
[F_Trap_4, V_Trap_4] =  STL_ReadFile([dir Trap_4],true);
[F_Trap_5, V_Trap_5] =  STL_ReadFile([dir Trap_5],true);
[F_Trap_6, V_Trap_6] =  STL_ReadFile([dir Trap_6],true);
[F_Trap_7, V_Trap_7] =  STL_ReadFile([dir Trap_7],true);
[F_Trap_8, V_Trap_8] =  STL_ReadFile([dir Trap_8],true);
[F_Trap_9, V_Trap_9] =  STL_ReadFile([dir Trap_9],true);
[F_Trap_10, V_Trap_10] =  STL_ReadFile([dir Trap_10],true);
[F_Trap_11, V_Trap_11] =  STL_ReadFile([dir Trap_11],true);
[F_Trap_12, V_Trap_12] =  STL_ReadFile([dir Trap_12],true);
[F_Trap_13, V_Trap_13] =  STL_ReadFile([dir Trap_13],true);
[F_Trap_14, V_Trap_14] =  STL_ReadFile([dir Trap_14],true);
[F_Trap_15, V_Trap_15] =  STL_ReadFile([dir Trap_15],true);
[F_Trap_16, V_Trap_16] =  STL_ReadFile([dir Trap_16],true);
[F_Trap_17, V_Trap_17] =  STL_ReadFile([dir Trap_17],true);
[F_Trap_18, V_Trap_18] =  STL_ReadFile([dir Trap_18],true);
[F_Trap_19, V_Trap_19] =  STL_ReadFile([dir Trap_19],true);

% Calculation of the shortest distance between MC1 and Trap for each frame

[idx1,D1]=knnsearch(V_MC1_1,V_Trap_1);
SD1=min(D1);
[idx2,D2]=knnsearch(V_MC1_2,V_Trap_2);
SD2=min(D2);
[idx3,D3]=knnsearch(V_MC1_3,V_Trap_3);
SD3=min(D3);
[idx4,D4]=knnsearch(V_MC1_4,V_Trap_4);
SD4=min(D4);
[idx5,D5]=knnsearch(V_MC1_5,V_Trap_5);
SD5=min(D5);
[idx6,D6]=knnsearch(V_MC1_6,V_Trap_6);
SD6=min(D6);
[idx7,D7]=knnsearch(V_MC1_7,V_Trap_7);
SD7=min(D7);
[idx8,D8]=knnsearch(V_MC1_8,V_Trap_8);
SD8=min(D8);
[idx9,D9]=knnsearch(V_MC1_9,V_Trap_9);
SD9=min(D9);
[idx10,D10]=knnsearch(V_MC1_10,V_Trap_10);
SD10=min(D10);
[idx11,D11]=knnsearch(V_MC1_11,V_Trap_11);
SD11=min(D11);
[idx12,D12]=knnsearch(V_MC1_12,V_Trap_12);
SD12=min(D12);
[idx13,D13]=knnsearch(V_MC1_13,V_Trap_13);
SD13=min(D13);
[idx14,D14]=knnsearch(V_MC1_14,V_Trap_14);
SD14=min(D14);
[idx15,D15]=knnsearch(V_MC1_15,V_Trap_15);
SD15=min(D15);
[idx16,D16]=knnsearch(V_MC1_16,V_Trap_16);
SD16=min(D16);
[idx17,D17]=knnsearch(V_MC1_17,V_Trap_17);
SD17=min(D17);
[idx18,D18]=knnsearch(V_MC1_18,V_Trap_18);
SD18=min(D18);
[idx19,D19]=knnsearch(V_MC1_19,V_Trap_19);
SD19=min(D19);

SD=[SD1;SD2;SD3;SD4;SD5;SD6;SD7;SD8;SD9;SD10;SD11;SD12;SD13;SD14;SD15;SD16;SD17;SD18;SD19];
plot(SD);
title('Evolution of the shortest distance between MC1 and Trap for each frame');
xlabel('frames');
ylabel('Shortest Distance (mm)');


