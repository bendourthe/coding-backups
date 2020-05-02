clear all; close all;
% Function: Plots the helical axis calculated with CTmotion_Static during
% Ab-Ad and Ex-Fl motion on the same figure for comparison and calculates
% the angles between the helical axis of each position.
%
% Input: STl files of the radius, scaphoid, trapezium and first metacarpal
%       (MC1) from the ABD scan.
%
% Outputs: 
%       -figure showing the 4 bones and the helical axis of the MC1,
%       trapezium and scaphoid
%       -figure showing the helical axis projection in the XY-plane
%       -figure showing the helical axis projection in the YZ-plane
%       -figure showing the helical axis projection in the XZ-plane
%       -displays the values of the angles between each helical axis (3D
%       and 2D projections)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Controle panel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Folder with files
dir = 'F:\Profesional\PhD KU Leuven\Data\Study_Brown_Rig\New Segmentation\Scan30\STL\'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILES
% ABD scan
MC1 = 'SCAN30_ABD_mc1.stl';
Rad = 'SCAN30_ABD_rad_reg.stl';
Scaph = 'SCAN30_ABD_scaph.stl';
Trap = 'SCAN30_ABD_trap.stl';
% EXT scan
Rad2 = 'SCAN30_EXT_rad.stl';

% Helical axis shift coefficient
s = 0.5;

        % Loads all the STL files
    [F_Rad, V_Rad] =  STL_ReadFile([dir Rad],true);
    [F_Scaph, V_Scaph] =  STL_ReadFile([dir Scaph],true);
    [F_Trap, V_Trap] =  STL_ReadFile([dir Trap],true);
    [F_MC1, V_MC1] =  STL_ReadFile([dir MC1],true);
    
    [F_Rad2, V_Rad2] =  STL_ReadFile([dir Rad2],true);
    
% Data
%   Coordinate system (Abd-Add)
radius_axis_available = 1; % loads data if available (=1), else, generate data
    filenameAX = 'RadAxes_Abd_Add.txt';
    fileAX = [dir filenameAX];
%   Helical axis    
helical_axis_Ab_Ad_available = 1; % loads data if available (=1), else, error
    filenameHelicalAX_Ab_Ad = 'HelicalAx_Abd_Add.txt';
    fileHelicalAX_Ab_Ad = [dir filenameHelicalAX_Ab_Ad];
helical_axis_Ex_Fl_available = 1; % loads data if available (=1), else, error
    filenameHelicalAX_Ex_Fl = 'HelicalAx_Ext_Flex.txt';
    fileHelicalAX_Ex_Fl = [dir filenameHelicalAX_Ex_Fl];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% Load CS data to make sure everything is expressed within the same CS
if(radius_axis_available)
    A = txt2mat(fileAX);
    
    O = A(1,:);
    Y = A(2,:);
    Z = A(3,:);
end
Origin = O;
Yaxis = (Y-O)/norm(Y-O);
Zaxis = (Z-O)/norm(Z-O);
Xaxis = cross(Yaxis,Zaxis)/norm(cross(Yaxis,Zaxis));
Zaxis = cross(Xaxis,Yaxis)/norm(cross(Xaxis,Yaxis));

% Converges to radius CS
RotMat = [Xaxis', Yaxis', Zaxis'];
loc2glob = [RotMat, Origin';
     0 0 0 1]; 
   % These tranformation matrices convert from local to global CS.
   % Invert them to go from global to local.
glob2loc = [RotMat', -(RotMat')*Origin'; 0 0 0 1];
R_glob2loc = glob2loc(1:3,1:3);
T_glob2loc = glob2loc(1:3,4);

    % Registers radius in position 2 with the radius in position 1
    [R_Rad, T_Rad, ERROR_Rad] = RegisterBones(F_Rad2, V_Rad2,F_Rad, V_Rad,1);
    H_Rad = [R_Rad, T_Rad; 0 0 0 1];
    
    % Transformes all data to radius CS
V_Rad = R_glob2loc*V_Rad.' + repmat(T_glob2loc,1,size(V_Rad,1));
V_Rad = V_Rad.';
V_MC1 = (glob2loc * [V_MC1.' ; ones(1,size(V_MC1,1))]).';
V_MC1 = V_MC1(:,1:3);
V_Scaph = (glob2loc * [V_Scaph.' ; ones(1,size(V_Scaph,1))]).';
V_Scaph = V_Scaph(:,1:3);
V_Trap = (glob2loc * [V_Trap.' ; ones(1,size(V_Trap,1))]).';
V_Trap = V_Trap(:,1:3);

V_Rad2 = (glob2loc * H_Rad * [V_Rad2.' ; ones(1,size(V_Rad2,1))]).';
V_Rad2 = V_Rad2(:,1:3);

% test: transform the first segments so it should match-up with the
% second one.
    h1 = figure;
    [obj, li, ax] = GUI_PlotShells(h1, {F_Rad;F_Rad2}, {V_Rad;V_Rad2},...
            {ones(size(V_Rad,1),1),ones(size(V_Rad2,1),1)});
    box off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Helical Axis
if(helical_axis_Ab_Ad_available)
    A_Ab_Ad = txt2mat(fileHelicalAX_Ab_Ad);
    
    p_Ab_Ad_mc1 = A_Ab_Ad(1,:);
    n_Ab_Ad_mc1 = A_Ab_Ad(2,:);
    p_Ab_Ad_trap = A_Ab_Ad(3,:);
    n_Ab_Ad_trap = A_Ab_Ad(4,:);
    p_Ab_Ad_scaph = A_Ab_Ad(5,:);
    n_Ab_Ad_scaph = A_Ab_Ad(6,:);
end
if(helical_axis_Ex_Fl_available)
    A_Ex_Fl = txt2mat(fileHelicalAX_Ex_Fl);
    
    p_Ex_Fl_mc1 = A_Ex_Fl(1,:);
    n_Ex_Fl_mc1 = A_Ex_Fl(2,:);
    p_Ex_Fl_trap = A_Ex_Fl(3,:);
    n_Ex_Fl_trap = A_Ex_Fl(4,:);
    p_Ex_Fl_scaph = A_Ex_Fl(5,:);
    n_Ex_Fl_scaph = A_Ex_Fl(6,:);
end

p_Ex_Fl_mc1 = R_glob2loc*p_Ex_Fl_mc1.' + repmat(T_glob2loc,1,size(p_Ex_Fl_mc1,1));
p_Ex_Fl_mc1 = p_Ex_Fl_mc1.';
n_Ex_Fl_mc1 = R_glob2loc*n_Ex_Fl_mc1.' + repmat(T_glob2loc,1,size(n_Ex_Fl_mc1,1));
n_Ex_Fl_mc1 = n_Ex_Fl_mc1.';
p_Ex_Fl_trap = R_glob2loc*p_Ex_Fl_trap.' + repmat(T_glob2loc,1,size(p_Ex_Fl_trap,1));
p_Ex_Fl_trap = p_Ex_Fl_trap.';
n_Ex_Fl_trap = R_glob2loc*n_Ex_Fl_trap.' + repmat(T_glob2loc,1,size(n_Ex_Fl_trap,1));
n_Ex_Fl_trap = n_Ex_Fl_trap.';
p_Ex_Fl_scaph = R_glob2loc*p_Ex_Fl_scaph.' + repmat(T_glob2loc,1,size(p_Ex_Fl_scaph,1));
p_Ex_Fl_scaph = p_Ex_Fl_scaph.';
n_Ex_Fl_scaph = R_glob2loc*n_Ex_Fl_scaph.' + repmat(T_glob2loc,1,size(n_Ex_Fl_scaph,1));
n_Ex_Fl_scaph = n_Ex_Fl_scaph.';

% Plots STl files and the helical axes
    h = figure;
    [obj, li, ax] = GUI_PlotShells(h, {F_Rad;F_Scaph;F_Trap;F_MC1},...
        {V_Rad;V_Scaph;V_Trap;V_MC1},...
        {ones(size(V_Rad,1),1),ones(size(V_Scaph,1),1),ones(size(V_Trap,1),1),ones(size(V_MC1,1),1)});
    hold on
    arrow([0;0;0],[30;0;0],10,70,30,5,'EdgeColor','r','FaceColor','r');
    arrow([0;0;0],[0;30;0],10,70,30,5,'EdgeColor','g','FaceColor','g');
    arrow([0;0;0],[0;0;30],10,70,30,5,'EdgeColor','b','FaceColor','b');
    
        % Plots the helical axis for each joint
    arrow(p_Ab_Ad_mc1-200*n_Ab_Ad_mc1,p_Ab_Ad_mc1+200*n_Ab_Ad_mc1,5,70,30,'EdgeColor','r','FaceColor','r','LineWidth',4);
    arrow(p_Ex_Fl_mc1-200*n_Ex_Fl_mc1,p_Ex_Fl_mc1+200*n_Ex_Fl_mc1,5,70,30,'EdgeColor','m','FaceColor','m','LineWidth',4);
    arrow(p_Ab_Ad_trap-200*n_Ab_Ad_trap,p_Ab_Ad_trap+200*n_Ab_Ad_trap,5,70,30,'EdgeColor','b','FaceColor','b','LineWidth',4);
    arrow(p_Ex_Fl_trap-200*n_Ex_Fl_trap,p_Ex_Fl_trap+200*n_Ex_Fl_trap,5,70,30,'EdgeColor','c','FaceColor','c','LineWidth',4);
    arrow(p_Ab_Ad_scaph-200*n_Ab_Ad_scaph,p_Ab_Ad_scaph+200*n_Ab_Ad_scaph,5,70,30,'EdgeColor','y','FaceColor','y','LineWidth',4);
    arrow(p_Ex_Fl_scaph-200*n_Ex_Fl_scaph,p_Ex_Fl_scaph+200*n_Ex_Fl_scaph,5,70,30,'EdgeColor','g','FaceColor','g','LineWidth',4);

    title('Coordinate System (CS) and Helical Axis (HA) of bone chain (n=4)');
    legend('Radius','Scaphoid','Trapezium','MC1','CS: X-axis','CS: Y-axis','CS: Z-axis','HA: MC1 - Ab-Ad',...
        'HA: MC1 - Ex-Fl','HA: Trap - Ab-Ad','HA: Trap - Ex-Fl','HA: Scaph - Ab-Ad','HA: Scaph - Ex-Fl','Location','SouthEastOutside')

    % Plots the helical axis in 2D for angle observation
        % XY-Plane
        n_aa_xy_mc1 = [n_Ab_Ad_mc1(2)/sqrt(n_Ab_Ad_mc1(1)^2+n_Ab_Ad_mc1(2)^2) n_Ab_Ad_mc1(1)/sqrt(n_Ab_Ad_mc1(1)^2+n_Ab_Ad_mc1(2)^2)];
        n_ef_xy_mc1 = [n_Ex_Fl_mc1(2)/sqrt(n_Ex_Fl_mc1(1)^2+n_Ex_Fl_mc1(2)^2) n_Ex_Fl_mc1(1)/sqrt(n_Ex_Fl_mc1(1)^2+n_Ex_Fl_mc1(2)^2)];
        n_aa_xy_trap = [n_Ab_Ad_trap(2)/sqrt(n_Ab_Ad_trap(1)^2+n_Ab_Ad_trap(2)^2) n_Ab_Ad_trap(1)/sqrt(n_Ab_Ad_trap(1)^2+n_Ab_Ad_trap(2)^2)];
        n_ef_xy_trap = [n_Ex_Fl_trap(2)/sqrt(n_Ex_Fl_trap(1)^2+n_Ex_Fl_trap(2)^2) n_Ex_Fl_trap(1)/sqrt(n_Ex_Fl_trap(1)^2+n_Ex_Fl_trap(2)^2)];
        n_aa_xy_scaph = [n_Ab_Ad_scaph(2)/sqrt(n_Ab_Ad_scaph(1)^2+n_Ab_Ad_scaph(2)^2) n_Ab_Ad_scaph(1)/sqrt(n_Ab_Ad_scaph(1)^2+n_Ab_Ad_scaph(2)^2)];
        n_ef_xy_scaph = [n_Ex_Fl_scaph(2)/sqrt(n_Ex_Fl_scaph(1)^2+n_Ex_Fl_scaph(2)^2) n_Ex_Fl_scaph(1)/sqrt(n_Ex_Fl_scaph(1)^2+n_Ex_Fl_scaph(2)^2)];
        
        h_xy_mc1 = figure;
        arrow(-n_aa_xy_mc1,n_aa_xy_mc1,5,70,30,'EdgeColor','r','FaceColor','r','LineWidth',4);
        arrow(-n_ef_xy_mc1,n_ef_xy_mc1,5,70,30,'EdgeColor','m','FaceColor','m','LineWidth',4);
        arrow(s-n_aa_xy_trap,s+n_aa_xy_trap,5,70,30,'EdgeColor','b','FaceColor','b','LineWidth',4);
        arrow(s-n_ef_xy_trap,s+n_ef_xy_trap,5,70,30,'EdgeColor','c','FaceColor','c','LineWidth',4);
        arrow(2*s-n_aa_xy_scaph,2*s+n_aa_xy_scaph,5,70,30,'EdgeColor','y','FaceColor','y','LineWidth',4);
        arrow(2*s-n_ef_xy_scaph,2*s+n_ef_xy_scaph,5,70,30,'EdgeColor','g','FaceColor','g','LineWidth',4);
            xlabel('X axis');
            ylabel('Y axis');
            title('Helical Axis (HA) projection (XY-Plane)');
            legend('HA: MC1 - Ab-Ad','HA: MC1 - Ex-Fl','HA: Trap - Ab-Ad','HA: Trap - Ex-Fl','HA: Scaph - Ab-Ad','HA: Scaph - Ex-Fl','Location','NorthEastOutside')     
        
        % YZ-Plane
        n_aa_yz_mc1 = [n_Ab_Ad_mc1(3)/sqrt(n_Ab_Ad_mc1(2)^2+n_Ab_Ad_mc1(3)^2) n_Ab_Ad_mc1(2)/sqrt(n_Ab_Ad_mc1(2)^2+n_Ab_Ad_mc1(3)^2)];
        n_ef_yz_mc1 = [n_Ex_Fl_mc1(3)/sqrt(n_Ex_Fl_mc1(2)^2+n_Ex_Fl_mc1(3)^2) n_Ex_Fl_mc1(2)/sqrt(n_Ex_Fl_mc1(2)^2+n_Ex_Fl_mc1(3)^2)];
        n_aa_yz_trap = [n_Ab_Ad_trap(3)/sqrt(n_Ab_Ad_trap(2)^2+n_Ab_Ad_trap(3)^2) n_Ab_Ad_trap(2)/sqrt(n_Ab_Ad_trap(2)^2+n_Ab_Ad_trap(3)^2)];
        n_ef_yz_trap = [n_Ex_Fl_trap(3)/sqrt(n_Ex_Fl_trap(2)^2+n_Ex_Fl_trap(3)^2) n_Ex_Fl_trap(2)/sqrt(n_Ex_Fl_trap(2)^2+n_Ex_Fl_trap(3)^2)];
        n_aa_yz_scaph = [n_Ab_Ad_scaph(3)/sqrt(n_Ab_Ad_scaph(2)^2+n_Ab_Ad_scaph(3)^2) n_Ab_Ad_scaph(2)/sqrt(n_Ab_Ad_scaph(2)^2+n_Ab_Ad_scaph(3)^2)];
        n_ef_yz_scaph = [n_Ex_Fl_scaph(3)/sqrt(n_Ex_Fl_scaph(2)^2+n_Ex_Fl_scaph(3)^2) n_Ex_Fl_scaph(2)/sqrt(n_Ex_Fl_scaph(2)^2+n_Ex_Fl_scaph(3)^2)];
        
        h_yz_mc1 = figure;
        arrow(-n_aa_yz_mc1,n_aa_yz_mc1,5,70,30,'EdgeColor','r','FaceColor','r','LineWidth',4);
        arrow(-n_ef_yz_mc1,n_ef_yz_mc1,5,70,30,'EdgeColor','m','FaceColor','m','LineWidth',4);
        arrow(s-n_aa_yz_trap,s+n_aa_yz_trap,5,70,30,'EdgeColor','b','FaceColor','b','LineWidth',4);
        arrow(s-n_ef_yz_trap,s+n_ef_yz_trap,5,70,30,'EdgeColor','c','FaceColor','c','LineWidth',4);
        arrow(2*s-n_aa_yz_scaph,2*s+n_aa_yz_scaph,5,70,30,'EdgeColor','y','FaceColor','y','LineWidth',4);
        arrow(2*s-n_ef_yz_scaph,2*s+n_ef_yz_scaph,5,70,30,'EdgeColor','g','FaceColor','g','LineWidth',4);
            xlabel('Y axis');
            ylabel('Z axis');
            title('Helical Axis (HA) projection (YZ-Plane)');
            legend('HA: MC1 - Ab-Ad','HA: MC1 - Ex-Fl','HA: Trap - Ab-Ad','HA: Trap - Ex-Fl','HA: Scaph - Ab-Ad','HA: Scaph - Ex-Fl','Location','NorthEastOutside')
        
        % XZ-Plane
        n_aa_xz_mc1 = [n_Ab_Ad_mc1(3)/sqrt(n_Ab_Ad_mc1(1)^2+n_Ab_Ad_mc1(3)^2) n_Ab_Ad_mc1(1)/sqrt(n_Ab_Ad_mc1(1)^2+n_Ab_Ad_mc1(3)^2)];
        n_ef_xz_mc1 = [n_Ex_Fl_mc1(3)/sqrt(n_Ex_Fl_mc1(1)^2+n_Ex_Fl_mc1(3)^2) n_Ex_Fl_mc1(1)/sqrt(n_Ex_Fl_mc1(1)^2+n_Ex_Fl_mc1(3)^2)];
        n_aa_xz_trap = [n_Ab_Ad_trap(3)/sqrt(n_Ab_Ad_trap(1)^2+n_Ab_Ad_trap(3)^2) n_Ab_Ad_trap(1)/sqrt(n_Ab_Ad_trap(1)^2+n_Ab_Ad_trap(3)^2)];
        n_ef_xz_trap = [n_Ex_Fl_trap(3)/sqrt(n_Ex_Fl_trap(1)^2+n_Ex_Fl_trap(3)^2) n_Ex_Fl_trap(1)/sqrt(n_Ex_Fl_trap(1)^2+n_Ex_Fl_trap(3)^2)];
        n_aa_xz_scaph = [n_Ab_Ad_scaph(3)/sqrt(n_Ab_Ad_scaph(1)^2+n_Ab_Ad_scaph(3)^2) n_Ab_Ad_scaph(1)/sqrt(n_Ab_Ad_scaph(1)^2+n_Ab_Ad_scaph(3)^2)];
        n_ef_xz_scaph = [n_Ex_Fl_scaph(3)/sqrt(n_Ex_Fl_scaph(1)^2+n_Ex_Fl_scaph(3)^2) n_Ex_Fl_scaph(1)/sqrt(n_Ex_Fl_scaph(1)^2+n_Ex_Fl_scaph(3)^2)];
        
        h_xz_mc1 = figure;
        arrow(-n_aa_xz_mc1,n_aa_xz_mc1,5,70,30,'EdgeColor','r','FaceColor','r','LineWidth',4);
        arrow(-n_ef_xz_mc1,n_ef_xz_mc1,5,70,30,'EdgeColor','m','FaceColor','m','LineWidth',4);
        arrow(s-n_aa_xz_trap,s+n_aa_xz_trap,5,70,30,'EdgeColor','b','FaceColor','b','LineWidth',4);
        arrow(s-n_ef_xz_trap,s+n_ef_xz_trap,5,70,30,'EdgeColor','c','FaceColor','c','LineWidth',4);
        arrow(2*s-n_aa_xz_scaph,2*s+n_aa_xz_scaph,5,70,30,'EdgeColor','y','FaceColor','y','LineWidth',4);
        arrow(2*s-n_ef_xz_scaph,2*s+n_ef_xz_scaph,5,70,30,'EdgeColor','g','FaceColor','g','LineWidth',4);
            xlabel('X axis');
            ylabel('Z axis');
            title('Helical Axis (HA) projection (XZ-Plane)');
            legend('HA: MC1 - Ab-Ad','HA: MC1 - Ex-Fl','HA: Trap - Ab-Ad','HA: Trap - Ex-Fl','HA: Scaph - Ab-Ad','HA: Scaph - Ex-Fl','Location','NorthEastOutside')
        
    % Calculation of the angle between each helical axis

    % MC1
    angle_MC1 = acos((n_Ab_Ad_mc1(1)*n_Ex_Fl_mc1(1)+n_Ab_Ad_mc1(2)*n_Ex_Fl_mc1(2)+n_Ab_Ad_mc1(3)*n_Ex_Fl_mc1(3))/(sqrt(n_Ab_Ad_mc1(1)^2+n_Ab_Ad_mc1(2)^2+n_Ab_Ad_mc1(3)^2)*sqrt(n_Ex_Fl_mc1(1)^2+n_Ex_Fl_mc1(2)^2+n_Ex_Fl_mc1(3)^2)))*180/pi
    mc1_XY = abs(atan(n_Ab_Ad_mc1(2)/n_Ab_Ad_mc1(1))-atan(n_Ex_Fl_mc1(2)/n_Ex_Fl_mc1(1)))*180/pi
    mc1_YZ = abs(atan(n_Ab_Ad_mc1(3)/n_Ab_Ad_mc1(2))-atan(n_Ex_Fl_mc1(3)/n_Ex_Fl_mc1(2)))*180/pi
    mc1_XZ = abs(atan(n_Ab_Ad_mc1(3)/n_Ab_Ad_mc1(1))-atan(n_Ex_Fl_mc1(3)/n_Ex_Fl_mc1(1)))*180/pi
    
    % Trapezium
    angle_Trap = acos((n_Ab_Ad_trap(1)*n_Ex_Fl_trap(1)+n_Ab_Ad_trap(2)*n_Ex_Fl_trap(2)+n_Ab_Ad_trap(3)*n_Ex_Fl_trap(3))/(sqrt(n_Ab_Ad_trap(1)^2+n_Ab_Ad_trap(2)^2+n_Ab_Ad_trap(3)^2)*sqrt(n_Ex_Fl_trap(1)^2+n_Ex_Fl_trap(2)^2+n_Ex_Fl_trap(3)^2)))*180/pi
    trap_XY = abs(atan(n_Ab_Ad_trap(2)/n_Ab_Ad_trap(1))-atan(n_Ex_Fl_trap(2)/n_Ex_Fl_trap(1)))*180/pi
    trap_YZ = abs(atan(n_Ab_Ad_trap(3)/n_Ab_Ad_trap(2))-atan(n_Ex_Fl_trap(3)/n_Ex_Fl_trap(2)))*180/pi
    trap_XZ = abs(atan(n_Ab_Ad_trap(3)/n_Ab_Ad_trap(1))-atan(n_Ex_Fl_trap(3)/n_Ex_Fl_trap(1)))*180/pi
    
    % Scaphoid
    angle_Scaph = acos((n_Ab_Ad_scaph(1)*n_Ex_Fl_scaph(1)+n_Ab_Ad_scaph(2)*n_Ex_Fl_scaph(2)+n_Ab_Ad_scaph(3)*n_Ex_Fl_scaph(3))/(sqrt(n_Ab_Ad_scaph(1)^2+n_Ab_Ad_scaph(2)^2+n_Ab_Ad_scaph(3)^2)*sqrt(n_Ex_Fl_scaph(1)^2+n_Ex_Fl_scaph(2)^2+n_Ex_Fl_scaph(3)^2)))*180/pi
    scaph_XY = abs(atan(n_Ab_Ad_scaph(2)/n_Ab_Ad_scaph(1))-atan(n_Ex_Fl_scaph(2)/n_Ex_Fl_scaph(1)))*180/pi
    scaph_YZ = abs(atan(n_Ab_Ad_scaph(3)/n_Ab_Ad_scaph(2))-atan(n_Ex_Fl_scaph(3)/n_Ex_Fl_scaph(2)))*180/pi
    scaph_XZ = abs(atan(n_Ab_Ad_scaph(3)/n_Ab_Ad_scaph(1))-atan(n_Ex_Fl_scaph(3)/n_Ex_Fl_scaph(1)))*180/pi
    