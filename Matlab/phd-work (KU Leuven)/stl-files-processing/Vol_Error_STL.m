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
dir = 'C:\Users\u0096636\Documents\Professional\PhD KU Leuven\Data\CAD tests\';

% FILES
Bone_Ref = 'Sphere1.stl';
Bone_Segm = 'Sphere2.stl';

% Read STL files
[F_Bone_Ref, V_Bone_Ref] =  STL_ReadFile([dir Bone_Ref],true);
[F_Bone_Segm, V_Bone_Segm] =  STL_ReadFile([dir Bone_Segm],true);

% Calculation of Total Volume and Total Area of each STL

[totalVolume_Bone_Ref,totalArea_Bone_Ref] = stlVolume(V_Bone_Ref,F_Bone_Ref);
[totalVolume_Bone_Segm,totalArea_Bone_Segm] = stlVolume(V_Bone_Segm,F_Bone_Segm);

% Calculation of the absolute Volume Difference between each bone (mm3)
VolDiff=abs(totalVolume_Bone_Ref-totalVolume_Bone_Segm)

% Calculation of the Volume Error (%)
VolError=abs(((totalVolume_Bone_Ref-totalVolume_Bone_Segm)/totalVolume_Bone_Ref)*100)