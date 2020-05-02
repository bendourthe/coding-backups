%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                               MatchBones                                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Functtion: Matches two pieces of bone by from manually selected STL
% files
%
% Dependencies:
%               Icp.m
%               STL_ReadFile.m 
%                   TRI_RemoveInvalidTriangles.m
%                       DeleteUnreferencedElements.m
% 
%                   TRI_RemoveBadlyConnectedTriangles.m
%                       TRI_Edges.m
%                       StackEqualElementIndices.m
%                           RunLengthDecode.m
%                           IncrementalRuns.m
%                       TRI_Normals.m
%                           VectorNorms.m
%                       DeleteUnreferencedElements.m
%                 GUI_Plotshells
%                     TRI_Merge.m
%                     GUI_DistributeColors.m
%                 
% 
%
% Input: Manual STL's
%
% Output: Plot of two overlaying bones
% 
% Created by: Frederik Van Eeghem, 2012
% FK, 17/04/2013: Added dependencies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%--------------------------------------------------------------------------
% Load data
%--------------------------------------------------------------------------
[F1, V1, fpath1, fpos1] =  STL_ReadFile;
[F2, V2, fpath2, fpos2] =  STL_ReadFile;


%--------------------------------------------------------------------------
% Plot bones
%--------------------------------------------------------------------------
h1 = figure;
[obj, li, ax] = GUI_PlotShells(h1,F1, V1, ones(size(V1,1),1));
h2 = figure;
[obj, li, ax] = GUI_PlotShells(h2,F2, V2, ones(size(V2,1),1));

%--------------------------------------------------------------------------
% Perform iterative closest point
%--------------------------------------------------------------------------
% Resample
V1_temp = V1; % (1:10:end,:);
V2_temp = V2; % (1:10:end,:);

% ICP
[R, T] = icp(V1_temp.',V2_temp.');


%--------------------------------------------------------------------------
% Plot both
%--------------------------------------------------------------------------
V2 = R*V2.'+repmat(T,1,size(V2,1)); V2 = V2.';
[obj, li, ax] = GUI_PlotShells({F1;F2}, {V1;V2},...
    {ones(size(V1,1),1),ones(size(V2,1),1)});