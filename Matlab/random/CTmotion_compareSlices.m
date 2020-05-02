clear all; close all;
% Function: Combining the data aqcuired in
% "CTmotion_plusplus_relToStatic.m" and calculating helical axes.
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
% relative tot the static scan. txt-files with inertia axes.
%
% Output: plot with helical axces of the dynamic scan (within the static
% radias CS)
% 
% FK, 17/4/2013, Added description and dependencies.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Controle panel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir = 'C:\Users\u0078973\Documents\KULAK Lokaal\PhD\Matlab\CTmotion\2012\STLs\STL\STL hand\';
    filename1 = 'RotTransData_tovStat4.txt';
    file1 = [dir filename1];
    filename2 = 'RotTransData_tovStat8.txt';
    file2 = [dir filename2];
    filename3 = 'RotTransData_tovStat12.txt';
    file3 = [dir filename3];
    filename4 = 'RotTransData_tovStat16.txt';
    file4 = [dir filename4];
    filename5 = 'RotTransData_tovStat19.txt';
    file5 = [dir filename5];
filenameInertialAX = 'InertialAx.txt';
fileInertialAX = [dir filenameInertialAX];
filenameAX = 'RadAxes_tovStat.txt';
fileAX = [dir filenameAX];
fig = 1; %1 to plot figure, else 0

Rad_static = 'Radius_statisch.stl';
MC_static = 'SE00_Static_MC1.stl';
Sca_static = 'SE00_Static_Scaph.stl';
Trap_static = 'SE00_Static_Trap.stl';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loading data and allocate variables
A1 = txt2mat(file1);
   
    R_Rad1 = A1(1:3,1:3);
    T_Rad1 = A1(1:3,4);
    H_Rad1 = [R_Rad1, T_Rad1; 0 0 0 1];
    R_MC1 = A1(4:6,1:3);
    T_MC1 = A1(4:6,4);
    R_Trap1 = A1(7:9,1:3);
    T_Trap1 = A1(7:9,4);
    R_Sca1 = A1(10:12,1:3);
    T_Sca1 = A1(10:12,4);
    
    ERROR_Rad1 = A1(13,1);
    ERROR_MC1 = A1(14,1);
    ERROR_Trap1 = A1(15,1);
    ERROR_Sca1 = A1(16,1);
    
A2 = txt2mat(file2);
   
    R_Rad2 = A2(1:3,1:3);
    T_Rad2 = A2(1:3,4);
    H_Rad2 = [R_Rad2, T_Rad2; 0 0 0 1];
    R_MC2 = A2(4:6,1:3);
    T_MC2 = A2(4:6,4);
    R_Trap2 = A2(7:9,1:3);
    T_Trap2 = A2(7:9,4);
    R_Sca2 = A2(10:12,1:3);
    T_Sca2 = A2(10:12,4);
    
    ERROR_Rad2 = A2(13,1);
    ERROR_MC2 = A2(14,1);
    ERROR_Trap2 = A2(15,1);
    ERROR_Sca2 = A2(16,1);
    
A3 = txt2mat(file3);
   
    R_Rad3 = A3(1:3,1:3);
    T_Rad3 = A3(1:3,4);
    H_Rad3 = [R_Rad3, T_Rad3; 0 0 0 1];
    R_MC3 = A3(4:6,1:3);
    T_MC3 = A3(4:6,4);
    R_Trap3 = A3(7:9,1:3);
    T_Trap3 = A3(7:9,4);
    R_Sca3 = A3(10:12,1:3);
    T_Sca3 = A3(10:12,4);
    
    ERROR_Rad3 = A3(13,1);
    ERROR_MC3 = A3(14,1);
    ERROR_Trap3 = A3(15,1);
    ERROR_Sca3 = A3(16,1);

A4 = txt2mat(file4);
   
    R_Rad4 = A4(1:3,1:3);
    T_Rad4 = A4(1:3,4);
    H_Rad4 = [R_Rad4, T_Rad4; 0 0 0 1];
    R_MC4 = A4(4:6,1:3);
    T_MC4 = A4(4:6,4);
    R_Trap4 = A4(7:9,1:3);
    T_Trap4 = A4(7:9,4);
    R_Sca4 = A4(10:12,1:3);
    T_Sca4 = A4(10:12,4);
    
    ERROR_Rad4 = A4(13,1);
    ERROR_MC4 = A4(14,1);
    ERROR_Trap4 = A4(15,1);
    ERROR_Sca4 = A4(16,1);
    
    
A5 = txt2mat(file5);
   
    R_Rad5 = A5(1:3,1:3);
    T_Rad5 = A5(1:3,4);
    H_Rad5 = [R_Rad5, T_Rad5; 0 0 0 1];
    R_MC5 = A5(4:6,1:3);
    T_MC5 = A5(4:6,4);
    R_Trap5 = A5(7:9,1:3);
    T_Trap5 = A5(7:9,4);
    R_Sca5 = A5(10:12,1:3);
    T_Sca5 = A5(10:12,4);
    
    ERROR_Rad5 = A5(13,1);
    ERROR_MC5 = A5(14,1);
    ERROR_Trap5 = A5(15,1);
    ERROR_Sca5 = A5(16,1);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% All transformations are in the coordinate system of the underlying 
% bone fragment of the static scan. Grouping transformations eg; frame 1 ->
% static, static -> frame 2 to get the transition from frame 1 to frame 2.

% MC
    H_MC_1to2 = [R_MC2.', -R_MC2.'*T_MC2; 0 0 0 1] * [R_MC1, T_MC1; 0 0 0 1];
    H_MC_2to3 = [R_MC3.', -R_MC3.'*T_MC3; 0 0 0 1] * [R_MC2, T_MC2; 0 0 0 1];
    H_MC_3to4 = [R_MC4.', -R_MC4.'*T_MC4; 0 0 0 1] * [R_MC3, T_MC3; 0 0 0 1];
    H_MC_4to5 = [R_MC5.', -R_MC5.'*T_MC5; 0 0 0 1] * [R_MC4, T_MC4; 0 0 0 1];
    
% Trap
    H_Trap_1to2 = [R_Trap2.', -R_Trap2.'*T_Trap2; 0 0 0 1] * [R_Trap1, T_Trap1; 0 0 0 1];
    H_Trap_2to3 = [R_Trap3.', -R_Trap3.'*T_Trap3; 0 0 0 1] * [R_Trap2, T_Trap2; 0 0 0 1];
    H_Trap_3to4 = [R_Trap4.', -R_Trap4.'*T_Trap4; 0 0 0 1] * [R_Trap3, T_Trap3; 0 0 0 1];
    H_Trap_4to5 = [R_Trap5.', -R_Trap5.'*T_Trap5; 0 0 0 1] * [R_Trap4, T_Trap4; 0 0 0 1];
    
% Sca
    H_Sca_1to2 = [R_Sca2.', -R_Sca2.'*T_Sca2; 0 0 0 1] * [R_Sca1, T_Sca1; 0 0 0 1];
    H_Sca_2to3 = [R_Sca3.', -R_Sca3.'*T_Sca3; 0 0 0 1] * [R_Sca2, T_Sca2; 0 0 0 1];
    H_Sca_3to4 = [R_Sca4.', -R_Sca4.'*T_Sca4; 0 0 0 1] * [R_Sca3, T_Sca3; 0 0 0 1];
    H_Sca_4to5 = [R_Sca5.', -R_Sca5.'*T_Sca5; 0 0 0 1] * [R_Sca4, T_Sca4; 0 0 0 1];
    
% Linking to radius

Ax = txt2mat(fileInertialAX);

Trap2Rad = [Ax(1:3,1:3) Ax(1:3,4); 0 0 0 1];
Sca2Rad = [Ax(4:6,1:3) Ax(4:6,4); 0 0 0 1];
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating helical axes


[n_MC1,point_MC1,phi_MC1,t_MC1] = screw(H_MC_1to2);
[n_MC2,point_MC2,phi_MC2,t_MC2] = screw(H_MC_2to3);
[n_MC3,point_MC3,phi_MC3,t_MC3] = screw(H_MC_3to4);
[n_MC4,point_MC4,phi_MC4,t_MC4] = screw(H_MC_4to5);
n_MC1 = (Trap2Rad * [n_MC1 ; 1]); n_MC1 = n_MC1(1:3)/norm(n_MC1(1:3));
point_MC1 = (Trap2Rad * [point_MC1 ; 1]); point_MC1 = point_MC1(1:3);
n_MC2 = (Trap2Rad * [n_MC2 ; 1]); n_MC2 = n_MC2(1:3)/norm(n_MC2(1:3));
point_MC2 = (Trap2Rad * [point_MC2 ; 1]); point_MC2 = point_MC2(1:3);
n_MC3 = (Trap2Rad * [n_MC3 ; 1]); n_MC3 = n_MC3(1:3)/norm(n_MC3(1:3));
point_MC3 = (Trap2Rad * [point_MC3 ; 1]); point_MC3 = point_MC3(1:3);
n_MC4 = (Trap2Rad * [n_MC4 ; 1]); n_MC4 = n_MC4(1:3)/norm(n_MC4(1:3));
point_MC4 = (Trap2Rad * [point_MC4 ; 1]); point_MC4 = point_MC4(1:3);

[n_Trap1,point_Trap1,phi_Trap1,t_Trap1] = screw(H_Trap_1to2);
[n_Trap2,point_Trap2,phi_Trap2,t_Trap2] = screw(H_Trap_2to3);
[n_Trap3,point_Trap3,phi_Trap3,t_Trap3] = screw(H_Trap_3to4);
[n_Trap4,point_Trap4,phi_Trap4,t_Trap4] = screw(H_Trap_4to5);
n_Trap1 = (Sca2Rad * [n_Trap1 ; 1]); n_Trap1 = n_Trap1(1:3)/norm(n_Trap1(1:3));
point_Trap1 = (Sca2Rad * [point_Trap1 ; 1]); point_Trap1 = point_Trap1(1:3);
n_Trap2 = (Sca2Rad * [n_Trap2 ; 1]); n_Trap2 = n_Trap2(1:3)/norm(n_Trap2(1:3));
point_Trap2 = (Sca2Rad * [point_Trap2 ; 1]); point_Trap2 = point_Trap2(1:3);
n_Trap3 = (Sca2Rad * [n_Trap3 ; 1]); n_Trap3 = n_Trap3(1:3)/norm(n_Trap3(1:3));
point_Trap3 = (Sca2Rad * [point_Trap3 ; 1]); point_Trap3 = point_Trap3(1:3);
n_Trap4 = (Sca2Rad * [n_Trap4 ; 1]); n_Trap4 = n_Trap4(1:3)/norm(n_Trap4(1:3));
point_Trap4 = (Sca2Rad * [point_Trap4 ; 1]); point_Trap4 = point_Trap4(1:3);

[n_Sca1,point_Sca1,phi_Sca1,t_Sca1] = screw(H_Sca_1to2);
[n_Sca2,point_Sca2,phi_Sca2,t_Sca2] = screw(H_Sca_2to3);
[n_Sca3,point_Sca3,phi_Sca3,t_Sca3] = screw(H_Sca_3to4);
[n_Sca4,point_Sca4,phi_Sca4,t_Sca4] = screw(H_Sca_4to5);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot to visualy check things
if(fig)
    % if fig = 1
    % loading data(stl's and radius CS)
    [F_Rad, V_Rad] =  STL_ReadFile([dir Rad_static],true);
    [F_MC, V_MC] =  STL_ReadFile([dir MC_static],true);
    [F_Trap, V_Trap] =  STL_ReadFile([dir Trap_static],true);
    [F_Sca, V_Sca] =  STL_ReadFile([dir Sca_static],true);
    A = txt2mat(fileAX);
    O = A(1,:);
    Y = A(2,:);
    Z = A(3,:);
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
    
    % Transforming data to radius CS (to plot one large figure)
    Sca2 = (glob2loc * [V_Sca.' ; ones(1,size(V_Sca,1))]).';
    Sca2 = Sca2(:,1:3);
    V_Rad = R_glob2loc*V_Rad.' + repmat(T_glob2loc,1,size(V_Rad,1)); V_Rad = V_Rad.';
    Trap2 = (glob2loc * [V_Trap.' ; ones(1,size(V_Trap,1))]).';
    Trap2 = Trap2(:,1:3);
    MC2 = (glob2loc * [V_MC.' ; ones(1,size(V_MC,1))]).';
    MC2 = MC2(:,1:3);
    
    % Plot stl's and helical axes
    h8 = figure;
    [obj, li, ax] = GUI_PlotShells(h8, {F_Rad;F_Trap;F_MC;F_Sca}, {V_Rad;Trap2;MC2;Sca2},...
            {ones(size(V_Rad,1),1),ones(size(Trap2,1),1),ones(size(MC2,1),1),ones(size(Sca2,1),1)});
    hold on
    arrow([0;0;0],[10;0;0],5,70,30);
    arrow([0;0;0],[0;10;0],5,70,30);
    arrow([0;0;0],[0;0;10],5,70,30);
    arrow(point_Sca1-200*n_Sca1,point_Sca1+200*n_Sca1,5,70,30,'EdgeColor','r','FaceColor','r');
    arrow(point_Trap1-200*n_Trap1,point_Trap1+200*n_Trap1,5,70,30,'EdgeColor','g','FaceColor','g');
    arrow(point_MC1-200*n_MC1,point_MC1+200*n_MC1,5,70,30,'EdgeColor','b','FaceColor','b');
    arrow(point_Sca2-200*n_Sca2,point_Sca2+200*n_Sca2,5,70,30,'EdgeColor','r','FaceColor','r');
    arrow(point_Trap2-200*n_Trap2,point_Trap2+200*n_Trap2,5,70,30,'EdgeColor','g','FaceColor','g');
    arrow(point_MC2-200*n_MC2,point_MC2+200*n_MC2,5,70,30,'EdgeColor','b','FaceColor','b');
    arrow(point_Sca3-200*n_Sca3,point_Sca3+200*n_Sca3,5,70,30,'EdgeColor','r','FaceColor','r');
    arrow(point_Trap3-200*n_Trap3,point_Trap3+200*n_Trap3,5,70,30,'EdgeColor','g','FaceColor','g');
    arrow(point_MC3-200*n_MC3,point_MC3+200*n_MC3,5,70,30,'EdgeColor','b','FaceColor','b');
    arrow(point_Sca4-200*n_Sca4,point_Sca4+200*n_Sca4,5,70,30,'EdgeColor','r','FaceColor','r');
    arrow(point_Trap4-200*n_Trap4,point_Trap4+200*n_Trap4,5,70,30,'EdgeColor','g','FaceColor','g');
    arrow(point_MC4-200*n_MC4,point_MC4+200*n_MC4,5,70,30,'EdgeColor','b','FaceColor','b');
    hold off;
end