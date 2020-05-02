clear all; close all;
% Function: Matching bone(fragements) onto bone from static scn using ICP
% and manual landmark registration.
%
% Dependencies:        
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
% 
%
% Input: txt-files with rotations and translation data from each scan/frame
% relative tot the static scan. txt-files with inertia axes. OR stl files.
%
% Output:txt-files with transformations, interia axes and coordinate
% systems.
% 
% FK, 17/4/2013, Added description and dependencies.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Controle panel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir = 'C:\mimics\Medical Imaging\4D CT opposition Faes\4D CT 24-6-2013\STL\'; % Folder with files
data_available = 0; % loading data if available (=1), else, generate data
    filename = 'RotTransData_tovStat20.txt';
    file = [dir filename];
radius_axis_available = 0; % loading data if available (=1), else, generate data
    filenameAX = 'RadAxes_tovStat20.txt';
    fileAX = [dir filenameAX];
inertial_axis_available = 0; % loading data if available (=1), else, generate data
filenameInertialAX = 'InertialAx20.txt';
fileInertialAX = [dir filenameInertialAX];


% FILES

% *_static is the static scan data
Rad_frag = '611L_phase20_rad.stl';
Rad_static = '611L_static_rad.stl';
MC_frag = '611L_phase20_mc1.stl';
MC_static = '611L_static_mc1.stl';
Sca1 = '611L_phase20_sca.stl';
Sca_static = '611L_static_sca.stl';
Trap1 = '611L_phase20_trp.stl';
Trap_static = '611L_static_trp.stl';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CS on radius:
%   origin: medial tuberculum at the level of the styloids.
%   y: towards the tuberositas radii
%   z: towards radial styloid (if corrected for orthoganal CS)
%   x: Perpendicual to the z-axis to form a right-handed CS.

if(radius_axis_available)
    A = txt2mat(fileAX);
    
    O = A(1,:);
    Y = A(2,:);
    Z = A(3,:);
else
    [F_Rad, V_Rad] =  STL_ReadFile([dir Rad_static],true);
    [ad_curve] = PlacePoints3({F_Rad}, {V_Rad}, {F_Rad}, {V_Rad}, 'Select origin, point on Y axis and point on Z axis (in that order!). End by pressing enter and exiting the figure');
    close all;
    % Retrieving coördinates of the points:
    O = ad_curve.points{end,1}(1,:);
    Y = ad_curve.points{end,1}(2,:);
    Z = ad_curve.points{end,1}(3,:);
    
    % Creating save file
    fid = fopen(fileAX,'w+');
    fprintf(fid,'%5.8f\t%5.8f\t%5.8f\r\n',[O,Y,Z]);
    fclose(fid);
end

Origin = O;
Yaxis = (Y-O)/norm(Y-O);
Zaxis = (Z-O)/norm(Z-O);
Xaxis = cross(Yaxis,Zaxis)/norm(cross(Yaxis,Zaxis));
Zaxis = cross(Xaxis,Yaxis)/norm(cross(Xaxis,Yaxis));

% Converging to radius CS
RotMat = [Xaxis', Yaxis', Zaxis'];
loc2glob = [RotMat, Origin';
     0 0 0 1]; 
   % These tranformation matrices convert from local to global CS.
   % Invert them to go from global to local.
glob2loc = [RotMat', -(RotMat')*Origin'; 0 0 0 1];
R_glob2loc = glob2loc(1:3,1:3);
T_glob2loc = glob2loc(1:3,4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Searching for radius fragment
if(data_available)
    % Load data
    A = txt2mat(file);
    
    [F_Rad_end_init, V_Rad_end_init] =  STL_ReadFile([dir Rad_frag],true);
    [F_Rad, V_Rad] =  STL_ReadFile([dir Rad_static],true);
    [F_MC_end_init, V_MC_end_init] =  STL_ReadFile([dir MC_frag]',true);
    [F_MC, V_MC] =  STL_ReadFile([dir MC_static],true);
    [F_Trap_init, V_Trap_init] =  STL_ReadFile([dir Trap1],true);
    [F_Trap, V_Trap] =  STL_ReadFile([dir Trap_static],true);
    [F_Sca_init, V_Sca_init] =  STL_ReadFile([dir Sca1],true);
    [F_Sca, V_Sca] =  STL_ReadFile([dir Sca_static],true);
    
    R_Rad = A(1:3,1:3);
    T_Rad = A(1:3,4);
    H_Rad = [R_Rad, T_Rad; 0 0 0 1];
    R_MC = A(4:6,1:3);
    T_MC = A(4:6,4);
    R_Trap = A(7:9,1:3);
    T_Trap = A(7:9,4);
    R_Sca = A(10:12,1:3);
    T_Sca = A(10:12,4);
    
    ERROR_Rad = A(13,1);
    ERROR_MC = A(14,1);
    ERROR_Trap = A(15,1);
    ERROR_Sca = A(16,1);

    
    
else
    %   loading data
       fid = fopen(file,'w+');
        
    [F_Rad_end_init, V_Rad_end_init] =  STL_ReadFile([dir Rad_frag],true);
    [F_Rad, V_Rad] =  STL_ReadFile([dir Rad_static],true);
    [F_MC_end_init, V_MC_end_init] =  STL_ReadFile([dir MC_frag]',true);
    [F_MC, V_MC] =  STL_ReadFile([dir MC_static],true);
    [F_Trap_init, V_Trap_init] =  STL_ReadFile([dir Trap1],true);
    [F_Trap, V_Trap] =  STL_ReadFile([dir Trap_static],true);
    [F_Sca_init, V_Sca_init] =  STL_ReadFile([dir Sca1],true);
    [F_Sca, V_Sca] =  STL_ReadFile([dir Sca_static],true);
    
    %   Registering radius fragment in static CT scan
    [R_Rad, T_Rad, ERROR_Rad] = RegisterBones(F_Rad_end_init, V_Rad_end_init,F_Rad, V_Rad,1);
    H_Rad = [R_Rad, T_Rad; 0 0 0 1];
    
    
    % saving results to file
    W = reshape([R_Rad, T_Rad].',1,[]);
    fprintf(fid,'%5.8f\t%5.8f\t%5.8f\t%5.8f\r\n',W); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Transforming all data to radius CS (incl. ransforming dynamic scan to static CS)
Sca1 = (glob2loc * H_Rad * [V_Sca_init.' ; ones(1,size(V_Sca_init,1))]).';
Sca1 = Sca1(:,1:3);
Sca2 = (glob2loc * [V_Sca.' ; ones(1,size(V_Sca,1))]).';
Sca2 = Sca2(:,1:3);
V_Rad = R_glob2loc*V_Rad.' + repmat(T_glob2loc,1,size(V_Rad,1)); V_Rad = V_Rad.';
V_Rad_end = (glob2loc * H_Rad * [V_Rad_end_init.' ; ones(1,size(V_Rad_end_init,1))]).';
V_Rad_end = V_Rad_end(:,1:3);
Trap1 = (glob2loc * H_Rad * [V_Trap_init.' ; ones(1,size(V_Trap_init,1))]).';
Trap1 = Trap1(:,1:3);
Trap2 = (glob2loc * [V_Trap.' ; ones(1,size(V_Trap,1))]).';
Trap2 = Trap2(:,1:3);
MC1 = (glob2loc * H_Rad * [V_MC_end_init.' ; ones(1,size(V_MC_end_init,1))]).';
MC1 = MC1(:,1:3);
MC2 = (glob2loc * [V_MC.' ; ones(1,size(V_MC,1))]).';
MC2 = MC2(:,1:3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculating centroid motion
% Trapezium
cent_Trap1 = mean(Trap1);
cent_Trap2 = mean(Trap2);
% Scaphoïd
cent_Sca1 = mean(Sca1);
cent_Sca2 = mean(Sca2);
% Calculating centroid displacement (not realy necessary)
disp_Trap12 = cent_Trap2 - cent_Trap1;
disp_Sca12 = cent_Sca2 - cent_Sca1;
dist_TrapSca1 = cent_Sca1 - cent_Trap1;
dist_TrapSca2 = cent_Sca2 - cent_Trap2;

% Calculating principle axes of inertia. 

if(inertial_axis_available)
    % Loading data from file
    I = txt2mat(fileInertialAX);
    R_b_Trap = I(1:3,1:3).';
    cent_Trap2 = I(1:3,4).';
    R_b_Sca = I(4:6,1:3).';
    cent_Sca2 = I(4:6,4).';
    Trap2Rad = [R_b_Trap.' cent_Trap2'; 0 0 0 1];
    Rad2Trap = [R_b_Trap -R_b_Trap*cent_Trap2'; 0 0 0 1];
    Sca2Rad = [R_b_Sca.' cent_Sca2'; 0 0 0 1];
    Rad2Sca = [R_b_Sca -R_b_Sca*cent_Sca2'; 0 0 0 1];
else
    % Calculating interia axes
    [weights_Trap, normals_Trap] = TRI_Areas(F_Trap, Trap2, true);
    [~, R_b_Trap, loc_Trap] = KIN_DirectedSpringRotation(TRI_Centroids(F_Trap, Trap2), normals_Trap, weights_Trap);
    Trap2Rad = [R_b_Trap.' cent_Trap2'; 0 0 0 1];
    Rad2Trap = [R_b_Trap -R_b_Trap*cent_Trap2'; 0 0 0 1];

    [weights_Sca, normals_Sca] = TRI_Areas(F_Sca, Sca2, true);
    [~, R_b_Sca, loc_Sca] = KIN_DirectedSpringRotation(TRI_Centroids(F_Sca, Sca2), normals_Sca, weights_Sca);
    Sca2Rad = [R_b_Sca.' cent_Sca2'; 0 0 0 1];
    Rad2Sca = [R_b_Sca -R_b_Sca*cent_Sca2'; 0 0 0 1];
    
        % save to file
        fid2 = fopen(fileInertialAX,'w+');
        W = reshape([R_b_Trap.' cent_Trap2'].',1,[]);
        fprintf(fid2,'%5.8f\t%5.8f\t%5.8f\t%5.8f\r\n',W); 
        W = reshape([R_b_Sca.' cent_Sca2'].',1,[]);
        fprintf(fid2,'%5.8f\t%5.8f\t%5.8f\t%5.8f\r\n',W);
        fclose(fid2);
end
    % Plotting bone fragments and inertia axes
    h123 = figure;
    [obj, li, ax] = GUI_PlotShells(h123, {F_Trap,F_Sca}, {Trap2,Sca2},...
                {ones(size(Trap2,1),1),ones(size(Sca2,1),1)});
    hold on
    arrow(cent_Trap2,cent_Trap2+20*R_b_Trap(1,:),5,70,30,'EdgeColor','r','FaceColor','r');
    arrow(cent_Trap2,cent_Trap2+20*R_b_Trap(2,:),5,70,30,'EdgeColor','g','FaceColor','g');
    arrow(cent_Trap2,cent_Trap2+20*R_b_Trap(3,:),5,70,30,'EdgeColor','b','FaceColor','b');
    arrow(cent_Sca2,cent_Sca2+20*R_b_Sca(1,:),5,70,30,'EdgeColor','r','FaceColor','r');
    arrow(cent_Sca2,cent_Sca2+20*R_b_Sca(2,:),5,70,30,'EdgeColor','g','FaceColor','g');
    arrow(cent_Sca2,cent_Sca2+20*R_b_Sca(3,:),5,70,30,'EdgeColor','b','FaceColor','b');
    hold off

% Transforming data to lokal (radial) CS. Once transformed the ICP starts.
% The CS are defined in the static scan (coded as: Trap2,Sca2, etc.)

% Sca expressed in relation tot the radius CS
% Sca ifv lokaal RadiusAssenstelsel is al ok
if(~data_available)
    [R_Sca, T_Sca, ERROR_Sca] = RegisterBones(F_Sca_init, Sca1, F_Sca, Sca2,0);
end
Hom_Sca = [R_Sca, T_Sca; 0 0 0 1];
inv_Hom_Sca = [R_Sca.' -R_Sca.'*T_Sca; 0 0 0 1];

% Trap expressed in relation tot the Scaphoid CS
Trap1_loc = (Rad2Sca * Hom_Sca * [Trap1.' ; ones(1,size(Trap1,1))]).';
Trap1_loc = Trap1_loc(:,1:3);
Trap2_loc = (Rad2Sca * [Trap2.' ; ones(1,size(Trap2,1))]).';
Trap2_loc = Trap2_loc(:,1:3);
if(~data_available)
    [R_Trap, T_Trap, ERROR_Trap] = RegisterBones(F_Trap_init, Trap1_loc, F_Trap, Trap2_loc,0);
end
Hom_Trap = [R_Trap, T_Trap; 0 0 0 1];
inv_Hom_Trap = [R_Trap.' -R_Trap.'*T_Trap; 0 0 0 1];

% MC1 expressed in relation tot the trapezium CS.
% (as before, the transformations are necessary to calculare the motions of
% MC1 in relation to the static trapezium)

MC1_loc = (Rad2Trap * Sca2Rad * Hom_Trap * Rad2Sca * Hom_Sca *[MC1.' ; ones(1,size(MC1,1))]).';
MC1_loc = MC1_loc(:,1:3);
MC2_loc = (Rad2Trap * [MC2.' ; ones(1,size(MC2,1))]).';
MC2_loc = MC2_loc(:,1:3);
if(~data_available)
    [R_MC, T_MC, ERROR_MC] = RegisterBones(F_MC_end_init, MC1_loc, F_MC, MC2_loc,1);
end
Hom_MC = [R_MC, T_MC; 0 0 0 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Saving the gathered info (so you don't have to re-run the ICP)
if(~data_available)
    % Transformation of MC
    W = reshape([R_MC, T_MC].',1,[]);
    fprintf(fid,'%5.8f\t%5.8f\t%5.8f\t%5.8f\r\n',W);
    % Transformation of Trapezium
    W = reshape([R_Trap, T_Trap].',1,[]);
    fprintf(fid,'%5.8f\t%5.8f\t%5.8f\t%5.8f\r\n',W);
    % Transformation of Scaphoïd
    W = reshape([R_Sca, T_Sca].',1,[]);
    fprintf(fid,'%5.8f\t%5.8f\t%5.8f\t%5.8f\r\n',W);
    fprintf(fid,'%5.8f\r\n%5.8f\r\n%5.8f\r\n%5.8f\r\n',ERROR_Rad,ERROR_MC,ERROR_Trap,ERROR_Sca);
    fclose(fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Checks and testing

% Plotting all movement in relation to the radius CS.

    h8 = figure;
    [obj, li, ax] = GUI_PlotShells(h8, {F_Rad;F_Rad_end_init;F_Trap_init;F_Trap;F_MC_end_init;F_MC;F_Sca_init;F_Sca}, {V_Rad;V_Rad_end;Trap1;Trap2;MC1;MC2;Sca1;Sca2},...
            {ones(size(V_Rad,1),1),ones(size(V_Rad_end,1),1),ones(size(Trap1,1),1),ones(size(Trap2,1),1),ones(size(MC1,1),1),ones(size(MC2,1),1),ones(size(Sca1,1),1),ones(size(Sca2,1),1)});
    hold on
    arrow([0;0;0],[10;0;0],5,70,30);
    arrow([0;0;0],[0;10;0],5,70,30);
    arrow([0;0;0],[0;0;10],5,70,30);
    hold off;

% test: transforming the first segments so it should match-up with the
% second one.


    Trap1_loc2 = (Hom_Trap * [Trap1_loc.' ; ones(1,size(Trap1_loc,1))]).';
    Trap1_loc2 = Trap1_loc2(:,1:3);
    [h2] = figure
    h2 = figure;
    [obj, li, ax] = GUI_PlotShells(h2, {F_Trap_init;F_Trap}, {Trap1_loc2;Trap2_loc},...
            {ones(size(Trap1,1),1),ones(size(Trap2,1),1)});
        
       
    MC1_loc2 = (Hom_MC * [MC1_loc.' ; ones(1,size(MC1_loc,1))]).';
    MC1_loc2 = MC1_loc2(:,1:3);
    h4 = figure;
    [obj, li, ax] = GUI_PlotShells(h4, {F_MC_end_init;F_MC}, {MC1_loc2;MC2_loc},...
            {ones(size(MC1,1),1),ones(size(MC2,1),1)});
        
    Sca1 = (Hom_Sca * [Sca1.' ; ones(1,size(Sca1,1))]).';
    Sca1 = Sca1(:,1:3);
    h2 = figure;
    [obj, li, ax] = GUI_PlotShells(h2, {F_Sca_init;F_Sca}, {Sca1;Sca2},...
            {ones(size(Sca1,1),1),ones(size(Sca2,1),1)});
        
    h3 = figure;
    [obj, li, ax] = GUI_PlotShells(h3, {F_Rad;F_Rad_end_init}, {V_Rad;V_Rad_end},...
            {ones(size(V_Rad,1),1),ones(size(V_Rad_end,1),1)});