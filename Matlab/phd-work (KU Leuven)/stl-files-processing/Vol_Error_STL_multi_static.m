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
dir = 'C:\Users\u0096636\Documents\Professional\PhD KU Leuven\Data\Study_Brown_Rig\Scan1\STL\';

% FILES
Bone_Ref = 'N_MC1_S-1-0.4.stl';
Bone_Segm1 = 'P80_MC1_S-1-0.4.stl';
Bone_Segm2 = 'J80_MC1_S-1-0.4.stl';

% Read STL files
[F_Bone_Ref, V_Bone_Ref] =  STL_ReadFile([dir Bone_Ref],true);
[F_Bone_Segm1, V_Bone_Segm1] =  STL_ReadFile([dir Bone_Segm1],true);
[F_Bone_Segm2, V_Bone_Segm2] =  STL_ReadFile([dir Bone_Segm2],true);

% Calculation of Total Volume and Total Area of each STL

[totalVolume_Bone_Ref,totalArea_Bone_Ref] = stlVolume(V_Bone_Ref,F_Bone_Ref);
[totalVolume_Bone_Segm1,totalArea_Bone_Segm1] = stlVolume(V_Bone_Segm1,F_Bone_Segm1);
[totalVolume_Bone_Segm2,totalArea_Bone_Segm2] = stlVolume(V_Bone_Segm2,F_Bone_Segm2);

% Calculation of the absolute Volume Difference between each bone (mm3)
VolDiff1=abs(totalVolume_Bone_Ref-totalVolume_Bone_Segm1);
VolDiff2=abs(totalVolume_Bone_Ref-totalVolume_Bone_Segm2);
VolDiff=[VolDiff1,VolDiff2];

% Calculation of the Volume Error (%)
VolError1=abs(((totalVolume_Bone_Ref-totalVolume_Bone_Segm1)/totalVolume_Bone_Ref)*100);
VolError2=abs(((totalVolume_Bone_Ref-totalVolume_Bone_Segm2)/totalVolume_Bone_Ref)*100);
VolError=[VolError1,VolError2];

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