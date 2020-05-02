function [Hom_Rad_2to1,Hom_Scaph_2to1,Hom_Trap_2to1,Hom_MC1_2to1, ...
                Rot_Scaph,Rot_Trap,Rot_MC1] = KinFunctionRad(dir,CSRad,...
                    Rad1,Rad2,Scaph1,Scaph2,Trap1,Trap2,M1,M2,frame)
% Function: Automatic kinematics algorithm which:
%               - imports the radius coordinate system (CS)
%               - calculates the transformations of each bone from position
%               1 to position 2 (always relative to the previous bone)
%               - calculates the rotation angles using that Tait-Brian
%               convention
%
% Dependencies:        
%               STL_ReadFile.m
%               RegisterBones.m
%                   PlacePoints3.m
%                   TRI_RemoveInvalidTriangles.m
%                       DeleteUnreferencedElements.m% 
%                   TRI_RemoveBadlyConnectedTriangles.m
%                       TRI_Edges.m
%                       StackEqualElementIndices.m
%                           RunLengthDecode.m
%                           IncrementalRuns.m
%                       TRI_Normals.m
%                           VectorNorms.m
%                       DeleteUnreferencedElements.m
%               TRI_Areas.m
%               TRI_Centroids
%               KIN_DirectedSpringRotation.m
%               GUI_PlotShells.m
%               arrow.m
% 
% Inputs:
%       - directory
%       - coordinate system (radius)
%       - stl files of each bones (frame 1 and frame 2, respect order!)
%       - number of frames
%
% Output: 
%       - transformation matrices (frame 1 to frame 2)
%       - rotation angles (expressed between frames, relatively to the
%       previous bone)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import CS on radius:
%   origin: lowest point on the distal border of the ulnar notch
%   y: straight down (few mm) on the proximal border of the ulnar notch
%   z: on the tip of the radial styloid
%   x: perpendicual to the z-axis to form a right-handed CS (calculated)
    
O = CSRad(1,:);
Y = CSRad(2,:);
Z = CSRad(3,:);

Origin = O;
Yaxis = (Y-O)/norm(Y-O);
Zaxis = (Z-O)/norm(Z-O);
Xaxis = cross(Yaxis,Zaxis)/norm(cross(Yaxis,Zaxis));
Zaxis = cross(Xaxis,Yaxis)/norm(cross(Xaxis,Yaxis));

% Converges to radius CS
RotMat = [Xaxis', Yaxis', Zaxis'];
loc2glob = [RotMat, Origin'; 0 0 0 1]; 
   % These tranformation matrices convert from local to global CS.
   % Invert them to go from global to local.
glob2loc = [RotMat', -(RotMat')*Origin'; 0 0 0 1];
R_glob2loc = glob2loc(1:3,1:3);
T_glob2loc = glob2loc(1:3,4);

    % Loads all the STL files
    [F_MC1_P1, V_MC1_P1] =  STL_ReadFile([dir M1],true);
    [F_Rad_P1, V_Rad_P1] =  STL_ReadFile([dir Rad1],true);
    [F_Scaph_P1, V_Scaph_P1] =  STL_ReadFile([dir Scaph1],true);
    [F_Trap_P1, V_Trap_P1] =  STL_ReadFile([dir Trap1],true);
    
    [F_MC1_P2, V_MC1_P2] =  STL_ReadFile([dir M2]',true);
    [F_Rad_P2, V_Rad_P2] =  STL_ReadFile([dir Rad2],true);
    [F_Scaph_P2, V_Scaph_P2] =  STL_ReadFile([dir Scaph2],true);    
    [F_Trap_P2, V_Trap_P2] =  STL_ReadFile([dir Trap2],true);
    
    % Transforms radius data to radius CS
    V_Rad_P1_cs = R_glob2loc*V_Rad_P1.' + repmat(T_glob2loc,1,size(V_Rad_P1,1));
    V_Rad_P1_cs = V_Rad_P1_cs.';
    V_Rad_P2_cs = R_glob2loc*V_Rad_P2.' + repmat(T_glob2loc,1,size(V_Rad_P2,1));
    V_Rad_P2_cs = V_Rad_P2_cs.';
    

    % Registers radius (position 2 to position 1)
    [R_Rad_2to1, T_Rad_2to1, ERROR_Rad_2to1] = ...
                RegisterBones(F_Rad_P2, V_Rad_P2,F_Rad_P1, V_Rad_P1,1);
    Hom_Rad_2to1 = [R_Rad_2to1, T_Rad_2to1; 0 0 0 1];
    
    % Transforms all data to radius CS
    V_MC1_P1_cs = (glob2loc * [V_MC1_P1.' ; ones(1,size(V_MC1_P1,1))]).';
    V_MC1_P1_cs = V_MC1_P1_cs(:,1:3);
    V_Scaph_P1_cs = (glob2loc * [V_Scaph_P1.' ; ones(1,size(V_Scaph_P1,1))]).';
    V_Scaph_P1_cs = V_Scaph_P1_cs(:,1:3);
    V_Trap_P1_cs = (glob2loc * [V_Trap_P1.' ; ones(1,size(V_Trap_P1,1))]).';
    V_Trap_P1_cs = V_Trap_P1_cs(:,1:3);

    V_MC1_P2_cs = (glob2loc * [V_MC1_P2.' ; ones(1,size(V_MC1_P2,1))]).';
    V_MC1_P2_cs = V_MC1_P2_cs(:,1:3);
    V_Scaph_P2_cs = (glob2loc * [V_Scaph_P2.' ; ones(1,size(V_Scaph_P2,1))]).';
    V_Scaph_P2_cs = V_Scaph_P2_cs(:,1:3);
    V_Trap_P2_cs = (glob2loc * [V_Trap_P2.' ; ones(1,size(V_Trap_P2,1))]).';
    V_Trap_P2_cs = V_Trap_P2_cs(:,1:3);
    
    % Calculates centroids motion
        % Scaphoid
        cent_Scaph_P1 = mean(V_Scaph_P1_cs);
        % Trapezium
        cent_Trap_P1 = mean(V_Trap_P1_cs);

        % Calculates inertia axes
        [weights_Trap, normals_Trap] = TRI_Areas(F_Trap_P1, V_Trap_P1_cs, true);
        [~, R_b_Trap, loc_Trap] = KIN_DirectedSpringRotation(TRI_Centroids(F_Trap_P1, V_Trap_P1_cs), normals_Trap, weights_Trap);
        Trap2Rad = [R_b_Trap.' cent_Trap_P1'; 0 0 0 1];
        Rad2Trap = [R_b_Trap -R_b_Trap*cent_Trap_P1'; 0 0 0 1];

        [weights_Sca, normals_Sca] = TRI_Areas(F_Scaph_P1, V_Scaph_P1_cs, true);
        [~, R_b_Sca, loc_Sca] = KIN_DirectedSpringRotation(TRI_Centroids(F_Scaph_P1, V_Scaph_P1_cs), normals_Sca, weights_Sca);
        Sca2Rad = [R_b_Sca.' cent_Scaph_P1'; 0 0 0 1];
        Rad2Sca = [R_b_Sca -R_b_Sca*cent_Scaph_P1'; 0 0 0 1];
    
    % Scaphoid expressed in the radius CS
    [R_Scaph_2to1, T_Scaph_2to1, ERROR_Scaph_2to1] = ...
                RegisterBones(F_Scaph_P2, V_Scaph_P2_cs, F_Scaph_P1, V_Scaph_P1_cs,0);
    Hom_Scaph_2to1 = [R_Scaph_2to1, T_Scaph_2to1; 0 0 0 1];

    % Trapezium expressed in the scaphoid CS
    V_Trap_P1_R2S = (Rad2Sca * [V_Trap_P1_cs.' ; ones(1,size(V_Trap_P1_cs,1))]).';
    V_Trap_P1_R2S = -V_Trap_P1_R2S(:,1:3);
    V_Trap_P2_R2S = (Rad2Sca * Hom_Scaph_2to1 * [V_Trap_P2_cs.' ; ones(1,size(V_Trap_P2_cs,1))]).';
    V_Trap_P2_R2S = -V_Trap_P2_R2S(:,1:3);

    [R_Trap_2to1, T_Trap_2to1, ERROR_Trap_2to1] = ...
                RegisterBones(F_Trap_P2, V_Trap_P2_R2S, F_Trap_P1, V_Trap_P1_R2S,0);
    Hom_Trap_2to1 = [R_Trap_2to1, T_Trap_2to1; 0 0 0 1];

    % MC1 expressed in the trapezium CS.
    % (the transformations are necessary to calculate the motion of MC1
    % relative to the trapezium)
    V_MC1_P1_R2T = (Rad2Trap * [V_MC1_P1_cs.' ; ones(1,size(V_MC1_P1_cs,1))]).';
    V_MC1_P1_R2T = -V_MC1_P1_R2T(:,1:3);
    V_MC1_P2_R2T = (Rad2Trap * Sca2Rad * Hom_Trap_2to1 * Rad2Sca * Hom_Scaph_2to1 *[V_MC1_P2_cs.' ; ones(1,size(V_MC1_P2_cs,1))]).';
    V_MC1_P2_R2T = -V_MC1_P2_R2T(:,1:3);

    [R_MC1_2to1, T_MC1_2to1, ERROR_MC1_2to1] = ...
                RegisterBones(F_MC1_P2, V_MC1_P2_R2T, F_MC1_P1, V_MC1_P1_R2T,1);
    Hom_MC1_2to1 = [R_MC1_2to1, T_MC1_2to1; 0 0 0 1];
    
    
    % Check the accuracy of the registration by plotting each bone on its
    % registered counterpart
        % Radius
    hrad = figure;
    [obj, li, ax] = GUI_PlotShells(hrad, {F_Rad_P1;F_Rad_P2}, {V_Rad_P1_cs;V_Rad_P2_cs},...
            {ones(size(V_Rad_P1_cs,1),1),ones(size(V_Rad_P2_cs,1),1)});
    box off
    set(hrad,'numbertitle','off','name',sprintf('Radius registration #%d: PRESS Continue if OK, CLOSE figure if not',frame));
    set(hrad,'menubar','figure');
    uicontrol('Position',[20 20 200 40],'String','Continue',...
              'Callback','uiresume(gcbf)');
    uiwait(gcf);
    close(hrad);
    
        % Scaphoid
    V_Scaph_P2_r = (Hom_Scaph_2to1 * [V_Scaph_P2_cs.' ; ones(1,size(V_Scaph_P2_cs,1))]).';
    V_Scaph_P2_r = V_Scaph_P2_r(:,1:3);
    hscaph = figure;
    [obj, li, ax] = GUI_PlotShells(hscaph, {F_Scaph_P2;F_Scaph_P1}, {V_Scaph_P2_r;V_Scaph_P1_cs},...
            {ones(size(V_Scaph_P2_r,1),1),ones(size(V_Scaph_P1_cs,1),1)});
    box off
    set(hscaph,'numbertitle','off','name',sprintf('Scaphoid registration #%d: PRESS Continue if OK, CLOSE figure if not',frame));
    set(hscaph,'menubar','figure');
    % Plots 2 buttons: 'Bad registration' if the registration looks wrong,
    % this will stop the code, and 'Continue' if the registration
    % looks right, will continue to run the rest of the code
    uicontrol('Position',[20 20 200 40],'String','Continue',...
              'Callback','uiresume(gcbf)');
    uiwait(gcf);
    close(hscaph);
    
        % Trapezium
    V_Trap_P2_r = (Hom_Trap_2to1 * [V_Trap_P2_R2S.' ; ones(1,size(V_Trap_P2_cs,1))]).';
    V_Trap_P2_r = V_Trap_P2_r(:,1:3);
    htrap = figure;
    [obj, li, ax] = GUI_PlotShells(htrap, {F_Trap_P2;F_Trap_P1}, {V_Trap_P2_r;V_Trap_P1_R2S},...
            {ones(size(V_Trap_P2_r,1),1),ones(size(V_Trap_P1_cs,1),1)});
    box off
    set(htrap,'numbertitle','off','name',sprintf('Trapezium registration #%d: PRESS Continue if OK, CLOSE figure if not',frame));
    set(htrap,'menubar','figure');
    % Plots 2 buttons: 'Bad registration' if the registration looks wrong,
    % this will stop the code, and 'Continue' if the registration
    % looks right, will continue to run the rest of the code
    uicontrol('Position',[20 20 200 40],'String','Continue',...
              'Callback','uiresume(gcbf)');
    uiwait(gcf);
    close(htrap);
        
    V_MC1_P2_r = (Hom_MC1_2to1 * [V_MC1_P2_R2T.' ; ones(1,size(V_MC1_P2_R2T,1))]).';
    V_MC1_P2_r = V_MC1_P2_r(:,1:3);
    hmc1 = figure;
    [obj, li, ax] = GUI_PlotShells(hmc1, {F_MC1_P2;F_MC1_P1}, {V_MC1_P2_r;V_MC1_P1_R2T},...
            {ones(size(V_MC1_P2_r,1),1),ones(size(V_MC1_P1_R2T,1),1)});
    box off
    set(hmc1,'numbertitle','off','name',sprintf('MC1 registration #%d: PRESS Continue if OK, CLOSE figure if not',frame));
    set(hmc1,'menubar','figure');
    % Plots 2 buttons: 'Bad registration' if the registration looks wrong,
    % this will stop the code, and 'Continue' if the registration
    % looks right, will continue to run the rest of the code
    uicontrol('Position',[20 20 200 40],'String','Continue',...
              'Callback','uiresume(gcbf)');
    uiwait(gcf);
    close(hmc1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of the rotation angles along each axis of the radius
% coordinate system for each bone. Each rotation angle is expressed in the
% radius coordinate system relatively to the more proximal bone.

    % MC1
    
    Rot_MC1_x = asin(R_MC1_2to1(3,2)/cos(asin(-R_MC1_2to1(3,1))))*180/pi;
    Rot_MC1_y = asin(-R_MC1_2to1(3,1))*180/pi;
    Rot_MC1_z = asin(R_MC1_2to1(2,1)/cos(asin(-R_MC1_2to1(3,1))))*180/pi;
    Rot_MC1 = [Rot_MC1_x,Rot_MC1_y,Rot_MC1_z];
    
    % Trapezium
    
    Rot_Trap_x = asin(R_Trap_2to1(3,2)/cos(asin(-R_Trap_2to1(3,1))))*180/pi;
    Rot_Trap_y = asin(-R_Trap_2to1(3,1))*180/pi;
    Rot_Trap_z = asin(R_Trap_2to1(2,1)/cos(asin(-R_Trap_2to1(3,1))))*180/pi;
    Rot_Trap = [Rot_Trap_x,Rot_Trap_y,Rot_Trap_z];
    
    % Scaph
    
    Rot_Scaph_x = asin(R_Scaph_2to1(3,2)/cos(asin(-R_Scaph_2to1(3,1))))*180/pi;
    Rot_Scaph_y = asin(-R_Scaph_2to1(3,1))*180/pi;
    Rot_Scaph_z = asin(R_Scaph_2to1(2,1)/cos(asin(-R_Scaph_2to1(3,1))))*180/pi;
    Rot_Scaph = [Rot_Scaph_x,Rot_Scaph_y,Rot_Scaph_z];
    