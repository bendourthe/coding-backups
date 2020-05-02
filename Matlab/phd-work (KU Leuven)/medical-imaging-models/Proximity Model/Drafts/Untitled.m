% Selection of the current directory (where the STL files are)
dir = 'D:\PhD KU Leuven\Data\Study_Microarchitecture_Trapezium\raw STLs\';

% FILE
Bone = '4.1_microCT_Declercq Viviane.stl';

% Read STL file
[F_Bone, V_Bone] =  STL_ReadFile([dir Bone],true);