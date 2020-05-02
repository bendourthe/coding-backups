function [R, T, ERROR_final] = RegisterBones(Pinit,Pend,F_end,V_end,F1,V1,landmarks)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: uses ICP to match a fragment of bone onto a larger piece
%           of bone (automatically, knowing landmarks)
% 
% Input:    Pinit, Pend: data from initial landmark selection (frame 1)
%           F_end, V_end: stl data (from the smalles bone piece)
%
% Output:   R: rotationmatrix
%           T: translattionvector
% 
% Dependencies:
%   ICP_RegisterDatasets 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 5
    landmarks = 1;
end

% if 1 STL is only a bone fragment of another SLT that containt the whole
% bone:
if(landmarks) 
        
    % Marked landmarks (for initial positioning)
    points_init_V1 = Pinit;
    points_init_V_end = Pend;
    
    % ICP on initial landmark positioning 
    [R_init, T_init, ERROR_init] = icp(points_init_V1.',points_init_V_end.',20,4,2,1e-5);
    V_end = R_init*V_end.'+repmat(T_init,1,size(V_end,1)); V_end = V_end.';
    HomTot = [R_init, T_init; 0 0 0 1];


    % Final registration using ICP:
    % During each itteration pieces of the mesh from the large bone that ar not
    % found on the small bone fragment are discarded. This untill there is 
    % a best fit. 
    
    [R_icp,T_icp,~,ERROR_final] = ICP_RegisterDatasets(V_end,V1);
    R_icp = R_icp'; % inverteren (R is orthogonaal)
    Hom_icp = [R_icp, T_icp'; 0 0 0 1];
    HomTot = Hom_icp * HomTot;
    
else
    % Calculating ICP 
    V1_temp = V1(1:10:end,:);
    V_end_temp = V_end(1:10:end,:);
    [R, T, ERROR_final] = icp(V1_temp.',V_end_temp.',100,4,2,1e-4);
    ERROR_final = ERROR_final/size(V_end_temp,1);
    HomTot = [R, T; 0 0 0 1];
end

% Determining R and T:
R = HomTot(1:3,1:3);
T = HomTot(1:3,4);