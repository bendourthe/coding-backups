%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vol_Error_STL: Calculates the absolute volume difference between two
%                segmented bones and the volume error which can be linked
%                to the segmentation error                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dependencies:        
%               STL_ReadFile.m 
%                   TRI_RemoveInvalidTriangles.m
%                   TRI_RemoveBadlyConnectedTriangles.m
%               stlVolume.m
% 
% Input: 
%           Bone_Ref: STL file corresponding to the bone selected as the
%                       reference for error calculation (mostly: static case)
%           Bone_Segm: STL file corresponding to the segmented bone that we
%                       want to compare to the reference bone
% 
% Output:
%           VolDiff: Absolute volume difference between the two selected
%                       bones (in mm3)
%           VolError: Volume error calculated over the two selected bones
%                       (in %)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Selection of the current directory (where the STL files are)
dir = 'C:\Users\u0096636\Documents\Professional\PhD KU Leuven\Data\STL dynamic Faes\WS\';

% FILES
Bone_Ref = '611L_static_trp_WS.stl';
Bone_Segm1 = '611L_phase1_trp_WS.stl';
Bone_Segm2 = '611L_phase2_trp_WS.stl';
Bone_Segm3 = '611L_phase3_trp_WS.stl';
Bone_Segm4 = '611L_phase4_trp_WS.stl';
Bone_Segm5 = '611L_phase6_trp_WS.stl';
Bone_Segm6 = '611L_phase9_trp_WS.stl';
Bone_Segm7 = '611L_phase10_trp_WS.stl';
Bone_Segm8 = '611L_phase13_trp_WS.stl';
Bone_Segm9 = '611L_phase16_trp_WS.stl';
Bone_Segm10 = '611L_phase17_trp_WS.stl';
Bone_Segm11 = '611L_phase19_trp_WS.stl';
Bone_Segm12 = '611L_phase22_trp_WS.stl';
Bone_Segm13 = '611L_phase23_trp_WS.stl';
Bone_Segm14 = '611L_phase24_trp_WS.stl';
Bone_Segm15 = '611L_phase26_trp_WS.stl';
Bone_Segm16 = '611L_phase27_trp_WS.stl';
Bone_Segm17 = '611L_phase29_trp_WS.stl';
Bone_Segm18 = '611L_phase33_trp_WS.stl';
Bone_Segm19 = '611L_phase36_trp_WS.stl';

% Read STL files
[F_Bone_Ref, V_Bone_Ref] =  STL_ReadFile([dir Bone_Ref],true);
[F_Bone_Segm1, V_Bone_Segm1] =  STL_ReadFile([dir Bone_Segm1],true);
[F_Bone_Segm2, V_Bone_Segm2] =  STL_ReadFile([dir Bone_Segm2],true);
[F_Bone_Segm3, V_Bone_Segm3] =  STL_ReadFile([dir Bone_Segm3],true);
[F_Bone_Segm4, V_Bone_Segm4] =  STL_ReadFile([dir Bone_Segm4],true);
[F_Bone_Segm5, V_Bone_Segm5] =  STL_ReadFile([dir Bone_Segm5],true);
[F_Bone_Segm6, V_Bone_Segm6] =  STL_ReadFile([dir Bone_Segm6],true);
[F_Bone_Segm7, V_Bone_Segm7] =  STL_ReadFile([dir Bone_Segm7],true);
[F_Bone_Segm8, V_Bone_Segm8] =  STL_ReadFile([dir Bone_Segm8],true);
[F_Bone_Segm9, V_Bone_Segm9] =  STL_ReadFile([dir Bone_Segm9],true);
[F_Bone_Segm10, V_Bone_Segm10] =  STL_ReadFile([dir Bone_Segm10],true);
[F_Bone_Segm11, V_Bone_Segm11] =  STL_ReadFile([dir Bone_Segm11],true);
[F_Bone_Segm12, V_Bone_Segm12] =  STL_ReadFile([dir Bone_Segm12],true);
[F_Bone_Segm13, V_Bone_Segm13] =  STL_ReadFile([dir Bone_Segm13],true);
[F_Bone_Segm14, V_Bone_Segm14] =  STL_ReadFile([dir Bone_Segm14],true);
[F_Bone_Segm15, V_Bone_Segm15] =  STL_ReadFile([dir Bone_Segm15],true);
[F_Bone_Segm16, V_Bone_Segm16] =  STL_ReadFile([dir Bone_Segm16],true);
[F_Bone_Segm17, V_Bone_Segm17] =  STL_ReadFile([dir Bone_Segm17],true);
[F_Bone_Segm18, V_Bone_Segm18] =  STL_ReadFile([dir Bone_Segm18],true);
[F_Bone_Segm19, V_Bone_Segm19] =  STL_ReadFile([dir Bone_Segm19],true);

% Calculation of Total Volume and Total Area of each STL

[totalVolume_Bone_Ref,totalArea_Bone_Ref] = stlVolume(V_Bone_Ref,F_Bone_Ref);
[totalVolume_Bone_Segm1,totalArea_Bone_Segm1] = stlVolume(V_Bone_Segm1,F_Bone_Segm1);
[totalVolume_Bone_Segm2,totalArea_Bone_Segm2] = stlVolume(V_Bone_Segm2,F_Bone_Segm2);
[totalVolume_Bone_Segm3,totalArea_Bone_Segm3] = stlVolume(V_Bone_Segm3,F_Bone_Segm3);
[totalVolume_Bone_Segm4,totalArea_Bone_Segm4] = stlVolume(V_Bone_Segm4,F_Bone_Segm4);
[totalVolume_Bone_Segm5,totalArea_Bone_Segm5] = stlVolume(V_Bone_Segm5,F_Bone_Segm5);
[totalVolume_Bone_Segm6,totalArea_Bone_Segm6] = stlVolume(V_Bone_Segm6,F_Bone_Segm6);
[totalVolume_Bone_Segm7,totalArea_Bone_Segm7] = stlVolume(V_Bone_Segm7,F_Bone_Segm7);
[totalVolume_Bone_Segm8,totalArea_Bone_Segm8] = stlVolume(V_Bone_Segm8,F_Bone_Segm8);
[totalVolume_Bone_Segm9,totalArea_Bone_Segm9] = stlVolume(V_Bone_Segm9,F_Bone_Segm9);
[totalVolume_Bone_Segm10,totalArea_Bone_Segm10] = stlVolume(V_Bone_Segm10,F_Bone_Segm10);
[totalVolume_Bone_Segm11,totalArea_Bone_Segm11] = stlVolume(V_Bone_Segm11,F_Bone_Segm11);
[totalVolume_Bone_Segm12,totalArea_Bone_Segm12] = stlVolume(V_Bone_Segm12,F_Bone_Segm12);
[totalVolume_Bone_Segm13,totalArea_Bone_Segm13] = stlVolume(V_Bone_Segm13,F_Bone_Segm13);
[totalVolume_Bone_Segm14,totalArea_Bone_Segm14] = stlVolume(V_Bone_Segm14,F_Bone_Segm14);
[totalVolume_Bone_Segm15,totalArea_Bone_Segm15] = stlVolume(V_Bone_Segm15,F_Bone_Segm15);
[totalVolume_Bone_Segm16,totalArea_Bone_Segm16] = stlVolume(V_Bone_Segm16,F_Bone_Segm16);
[totalVolume_Bone_Segm17,totalArea_Bone_Segm17] = stlVolume(V_Bone_Segm17,F_Bone_Segm17);
[totalVolume_Bone_Segm18,totalArea_Bone_Segm18] = stlVolume(V_Bone_Segm18,F_Bone_Segm18);
[totalVolume_Bone_Segm19,totalArea_Bone_Segm19] = stlVolume(V_Bone_Segm19,F_Bone_Segm19);

% Calculation of the absolute Volume Difference between each bone (mm3)
VolDiff1=abs(totalVolume_Bone_Ref-totalVolume_Bone_Segm1);
VolDiff2=abs(totalVolume_Bone_Ref-totalVolume_Bone_Segm2);
VolDiff3=abs(totalVolume_Bone_Ref-totalVolume_Bone_Segm3);
VolDiff4=abs(totalVolume_Bone_Ref-totalVolume_Bone_Segm4);
VolDiff5=abs(totalVolume_Bone_Ref-totalVolume_Bone_Segm5);
VolDiff6=abs(totalVolume_Bone_Ref-totalVolume_Bone_Segm6);
VolDiff7=abs(totalVolume_Bone_Ref-totalVolume_Bone_Segm7);
VolDiff8=abs(totalVolume_Bone_Ref-totalVolume_Bone_Segm8);
VolDiff9=abs(totalVolume_Bone_Ref-totalVolume_Bone_Segm9);
VolDiff10=abs(totalVolume_Bone_Ref-totalVolume_Bone_Segm10);
VolDiff11=abs(totalVolume_Bone_Ref-totalVolume_Bone_Segm11);
VolDiff12=abs(totalVolume_Bone_Ref-totalVolume_Bone_Segm12);
VolDiff13=abs(totalVolume_Bone_Ref-totalVolume_Bone_Segm13);
VolDiff14=abs(totalVolume_Bone_Ref-totalVolume_Bone_Segm14);
VolDiff15=abs(totalVolume_Bone_Ref-totalVolume_Bone_Segm15);
VolDiff16=abs(totalVolume_Bone_Ref-totalVolume_Bone_Segm16);
VolDiff17=abs(totalVolume_Bone_Ref-totalVolume_Bone_Segm17);
VolDiff18=abs(totalVolume_Bone_Ref-totalVolume_Bone_Segm18);
VolDiff19=abs(totalVolume_Bone_Ref-totalVolume_Bone_Segm19);
VolDiff=[VolDiff1,VolDiff2,VolDiff3,VolDiff4,VolDiff5,VolDiff6,VolDiff7,VolDiff8,VolDiff9,VolDiff10,VolDiff11,VolDiff12,VolDiff13,VolDiff14,VolDiff15,VolDiff16,VolDiff17,VolDiff18,VolDiff19];

% Calculation of the Volume Error (%)
VolError1=abs(((totalVolume_Bone_Ref-totalVolume_Bone_Segm1)/totalVolume_Bone_Ref)*100);
VolError2=abs(((totalVolume_Bone_Ref-totalVolume_Bone_Segm2)/totalVolume_Bone_Ref)*100);
VolError3=abs(((totalVolume_Bone_Ref-totalVolume_Bone_Segm3)/totalVolume_Bone_Ref)*100);
VolError4=abs(((totalVolume_Bone_Ref-totalVolume_Bone_Segm4)/totalVolume_Bone_Ref)*100);
VolError5=abs(((totalVolume_Bone_Ref-totalVolume_Bone_Segm5)/totalVolume_Bone_Ref)*100);
VolError6=abs(((totalVolume_Bone_Ref-totalVolume_Bone_Segm6)/totalVolume_Bone_Ref)*100);
VolError7=abs(((totalVolume_Bone_Ref-totalVolume_Bone_Segm7)/totalVolume_Bone_Ref)*100);
VolError8=abs(((totalVolume_Bone_Ref-totalVolume_Bone_Segm8)/totalVolume_Bone_Ref)*100);
VolError9=abs(((totalVolume_Bone_Ref-totalVolume_Bone_Segm9)/totalVolume_Bone_Ref)*100);
VolError10=abs(((totalVolume_Bone_Ref-totalVolume_Bone_Segm10)/totalVolume_Bone_Ref)*100);
VolError11=abs(((totalVolume_Bone_Ref-totalVolume_Bone_Segm11)/totalVolume_Bone_Ref)*100);
VolError12=abs(((totalVolume_Bone_Ref-totalVolume_Bone_Segm12)/totalVolume_Bone_Ref)*100);
VolError13=abs(((totalVolume_Bone_Ref-totalVolume_Bone_Segm13)/totalVolume_Bone_Ref)*100);
VolError14=abs(((totalVolume_Bone_Ref-totalVolume_Bone_Segm14)/totalVolume_Bone_Ref)*100);
VolError15=abs(((totalVolume_Bone_Ref-totalVolume_Bone_Segm15)/totalVolume_Bone_Ref)*100);
VolError16=abs(((totalVolume_Bone_Ref-totalVolume_Bone_Segm16)/totalVolume_Bone_Ref)*100);
VolError17=abs(((totalVolume_Bone_Ref-totalVolume_Bone_Segm17)/totalVolume_Bone_Ref)*100);
VolError18=abs(((totalVolume_Bone_Ref-totalVolume_Bone_Segm18)/totalVolume_Bone_Ref)*100);
VolError19=abs(((totalVolume_Bone_Ref-totalVolume_Bone_Segm19)/totalVolume_Bone_Ref)*100);
VolError=[VolError1,VolError2,VolError3,VolError4,VolError5,VolError6,VolError7,VolError8,VolError9,VolError10,VolError11,VolError12,VolError13,VolError14,VolError15,VolError16,VolError17,VolError18,VolError19];

% Plot of each Volume Difference and Volume Error value in a two bar graphs
subplot(1,2,1);
bar(VolDiff,'gr');
title('Volume Difference');
xlabel('frames');
ylabel('Volume (mm3)');

subplot(1,2,2);
bar(VolError,'r');
title('Volume Error');
xlabel('frames');
ylabel('Error (%)');