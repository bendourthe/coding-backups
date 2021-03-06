clear all; close all;

% Function: Uses the txt-files with transformations and coordinate systems 
% generated by "CTmotion_plusplus_relToStatic", a plot with helical axes
% and average helical axes.
% 
% Dependencies: 
%   txt2mat.m
%   IntersectLineAndPlane.m 
%       DistanceFromVertexToPlane.m 
%           TRI_Normals.m	
%   VectorNorms.m
%   screw.m
%   FindPivotPoint
%   DistanceFromVertexToLine
%       VectorNorms.m
% 
% Created by: Frederik Van Eeghem, 2012
% FK, 17/4/2013: added dependencies

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Controle panel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine which files (phases/frames) are available
available = [1,2,3,4,5,6,8,9,10,12,13,15];

% Directery of txt-files with transformations and interial axes 
dir = 'C:\Users\u0096636\Documents\Professional\PhD KU Leuven\Data\Study_Faes_4D_CT\STL\';
    filename{1} = 'RotTransData_tovStat1.txt';
    filename{2} = 'RotTransData_tovStat2.txt';
    filename{3} = 'RotTransData_tovStat3.txt';
    filename{4} = 'RotTransData_tovStat4.txt';
    filename{5} = 'RotTransData_tovStat5.txt';
    filename{6} = 'RotTransData_tovStat6.txt';
    filename{7} = 'RotTransData_tovStat8.txt';
    filename{8} = 'RotTransData_tovStat9.txt';
    filename{9} = 'RotTransData_tovStat10.txt';
    filename{10} = 'RotTransData_tovStat12.txt';
    filename{11} = 'RotTransData_tovStat13.txt';
    filename{12} = 'RotTransData_tovStat15.txt';
    filename{13} = 'RotTransData_tovStat16.txt';
    filename{14} = 'RotTransData_tovStat17.txt';
    filename{15} = 'RotTransData_tovStat19.txt';

filenameInertialAX = 'InertialAx1.txt';
fileInertialAX = [dir filenameInertialAX];
filenameAX = 'RadAxes_tovStat20.txt';
fileAX = [dir filenameAX];
fig = 1; %1 to plot figure, else 0

% STL's filenames from each bone (static scan) 
Rad_static = '611L_static_rad.stl';
MC_static = '611L_static_mc1.stl';
Sca_static = '611L_static_sca.stl';
Trap_static = '611L_static_trp.stl';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------%
% Loading data %
%--------------%

% preallocating variables 

R_Rad = cell(size(available));
T_Rad = cell(size(available));
H_Rad = cell(size(available));
R_MC = cell(size(available));
T_MC = cell(size(available));
R_Trap = cell(size(available));
T_Trap = cell(size(available));
R_Sca = cell(size(available));
T_Sca = cell(size(available));

ERROR_Rad = cell(size(available));
ERROR_MC = cell(size(available));
ERROR_Trap = cell(size(available));
ERROR_Sca = cell(size(available));

% Loading data
for nr = available
    A = txt2mat([dir filename{nr}]);
    
    R_Rad{nr} = A(1:3,1:3);
    T_Rad{nr} = A(1:3,4);
    H_Rad{nr} = [R_Rad{nr}, T_Rad{nr}; 0 0 0 1];
    R_MC{nr} = A(4:6,1:3);
    T_MC{nr} = A(4:6,4);
    R_Trap{nr} = A(7:9,1:3);
    T_Trap{nr} = A(7:9,4);
    R_Sca{nr} = A(10:12,1:3);
    T_Sca{nr} = A(10:12,4);
    
    ERROR_Rad{nr} = A(13,1);
    ERROR_MC{nr} = A(14,1);
    ERROR_Trap{nr} = A(15,1);
    ERROR_Sca{nr} = A(16,1);
    
end

%---------------%
% Transforming %
%---------------%
% All transformations are in the coordinate system of the underlying 
% bone fragment of the static scan. Grouping transformations eg; frame 1 ->
% static, static -> frame 2 to get the transition from frame 1 to frame 2.

% Alle transformaties in assenstelsel van onderliggend botfragement van
% statische scan. Steeds transformaties samenstellen: 
% vb frame 1 -> statisch en dan statisch -> frame2 om zo de overgang 
% frame 1 -> frame2 te bekomen

% preallocating for speed
H_MC = cell(1,length(available)-1);
H_Trap = cell(1,length(available)-1);
H_Sca = cell(1,length(available)-1);

for nr = 1:length(available)-1
    H_MC{nr} = [R_MC{available(nr+1)}.', -R_MC{available(nr+1)}.'*T_MC{available(nr+1)}; 0 0 0 1] * [R_MC{available(nr)}, T_MC{available(nr)}; 0 0 0 1];
    H_Trap{nr} = [R_Trap{available(nr+1)}.', -R_Trap{available(nr+1)}.'*T_Trap{available(nr+1)}; 0 0 0 1] * [R_Trap{available(nr)}, T_Trap{available(nr)}; 0 0 0 1];
    H_Sca{nr} = [R_Sca{available(nr+1)}.', -R_Sca{available(nr+1)}.'*T_Sca{available(nr+1)}; 0 0 0 1] * [R_Sca{available(nr)}, T_Sca{available(nr)}; 0 0 0 1];
end

% Linking to radius
Ax = txt2mat(fileInertialAX);
Trap2Rad = [Ax(1:3,1:3) Ax(1:3,4); 0 0 0 1];
Sca2Rad = [Ax(4:6,1:3) Ax(4:6,4); 0 0 0 1];
    
%------------------------%
% Calculating helical axes %
%------------------------%

% preallocating for speed
n_MC = zeros(3,length(available)-1);
point_MC = zeros(3,length(available)-1);
phi_MC = zeros(1,length(available)-1);
t_MC = zeros(1,length(available)-1);
n_Trap = zeros(3,length(available)-1);
point_Trap = zeros(3,length(available)-1);
phi_Trap = zeros(1,length(available)-1);
t_Trap = zeros(1,length(available)-1);
n_Sca = zeros(3,length(available)-1);
point_Sca = zeros(3,length(available)-1);
phi_Sca = zeros(1,length(available)-1);
t_Sca = zeros(1,length(available)-1);

% Helical axis between frames
for nr = 1:length(available)-1
    [n_MC(:,nr),point_MC(:,nr),phi_MC(:,nr),t_MC(:,nr)] = screw(H_MC{nr},1);
    n_MC_tmp = (Trap2Rad * [n_MC(:,nr) ; 1]); n_MC(:,nr) = n_MC_tmp(1:3)/norm(n_MC_tmp(1:3));
    point_MC_tmp = (Trap2Rad * [point_MC(:,nr) ; 1]); point_MC(:,nr) = point_MC_tmp(1:3);
    [n_Trap(:,nr),point_Trap(:,nr),phi_Trap(:,nr),t_Trap(:,nr)] = screw(H_Trap{nr});
    n_Trap_tmp = (Sca2Rad * [n_Trap(:,nr) ; 1]); n_Trap(:,nr) = n_Trap_tmp(1:3)/norm(n_Trap_tmp(1:3));
    point_Trap_tmp = (Sca2Rad * [point_Trap(:,nr) ; 1]); point_Trap(:,nr) = point_Trap_tmp(1:3);
    [n_Sca(:,nr),point_Sca(:,nr),phi_Sca(:,nr),t_Sca(:,nr)] = screw(H_Sca{nr});
end

% Average helixal axis based on 2 points

lines = zeros((2*(length(available)-1))*3,3);
lines(1:2:end-1,:) = [point_Sca'; point_Trap'; point_MC'];
lines(2:2:end,:) = [(point_Sca+n_Sca)'; (point_Trap+n_Trap)'; (point_MC+n_MC)'];
planes = [-20 0 0; -20 1 0; -20 0 1; 30 0 0; 30 1 0; 30 0 1];
vertices = IntersectLineAndPlane(lines, planes);
average = [mean(vertices(1:(length(available)-1),:,1),1);
            mean(vertices(1:(length(available)-1),:,2),1);
            mean(vertices(length(available):2*(length(available)-1),:,1),1);
            mean(vertices(length(available):2*(length(available)-1),:,2),1);
            mean(vertices((2*(length(available)-1)+1):3*(length(available)-1),:,1),1);
            mean(vertices((2*(length(available)-1)+1):3*(length(available)-1),:,2),1)];
        
% Calculationg PivotPoint (could take a while)
[PivotPoint_Sca, ssd_Sca] = FindPivotPoint(point_Sca',n_Sca',[0 0 0]);
[PivotPoint_Trap, ssd_Trap] = FindPivotPoint(point_Trap',n_Trap',[0 0 0]);
[PivotPoint_MC, ssd_MC] = FindPivotPoint(point_MC',n_MC',[0 0 0]);
PivotPoints = [PivotPoint_Sca; PivotPoint_Trap; PivotPoint_MC];


%---------%
% Plotting %
%---------%

if(fig)
    % if fig = 1
    % loading data (stl's and radius CS)
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
    
    % Transforming data to radius CS (to make one large plot)
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
        
        % Plots the coordinate system and its origin
    hold on
    arrow([0;0;0],[10;0;0],5,70,30);
    arrow([0;0;0],[0;10;0],5,70,30);
    arrow([0;0;0],[0;0;10],5,70,30);
        
        % Plots the helical axis of each frame
    for nr = 1:length(available)-1
    arrow(point_Sca(:,nr)-200*n_Sca(:,nr),point_Sca(:,nr)+200*n_Sca(:,nr),5,70,30,'EdgeColor','r','FaceColor','r');
    arrow(point_Trap(:,nr)-200*n_Trap(:,nr),point_Trap(:,nr)+200*n_Trap(:,nr),5,70,30,'EdgeColor','g','FaceColor','g');
    arrow(point_MC(:,nr)-200*n_MC(:,nr),point_MC(:,nr)+200*n_MC(:,nr),5,70,30,'EdgeColor','b','FaceColor','b');
    end
        % Plots the average helical axis
    line(average(1:2,1),average(1:2,2),average(1:2,3),'Color','r','LineWidth',4)
    line(average(3:4,1),average(3:4,2),average(3:4,3),'Color','g','LineWidth',4)
    line(average(5:6,1),average(5:6,2),average(5:6,3),'Color','b','LineWidth',4)
    scatter3(PivotPoints(:,1),PivotPoints(:,2),PivotPoints(:,3),50,'black','filled');
    hold off;
end

% Definition of the horizontal plane defined on the radius
    % Equation of a plane: ax + by + cz + d = 0 (P)
        % To find the coefficient a, b and c we need 3 points
        % Here, we can use the points selected to define the radius
        % coordinate system, respectively called O, Y and z (see above)
        % We can find the equation of the plane by solving this equation:
        % OP.(OY x YZ) = 0     (with P(x,y,z))
        
        O = [0 0 0];
        Yaxis = [0 1 0];
        Zaxis = [0 0 1];
        OY = Yaxis - O;
        YZ = Zaxis - Yaxis;
        OYxYZ = [OY(1,2)*YZ(1,3)-YZ(1,2)*OY(1,3);...
            -(OY(1,1)*YZ(1,3)-YZ(1,1)*OY(1,3));...
            OY(1,1)*YZ(1,2)-YZ(1,1)*OY(1,2)];
        P = O;
        a = OYxYZ(1,1);
        b = OYxYZ(2,1);
        c = OYxYZ(3,1);
        d = a*P(1,1) + b*P(1,2) + c*P(1,3);

% Intersection between each helical axis and the plane OYZ (on the radius)
    % First, we define the parametric equation of each helical axis
        % xi = ai.t + pxi
        % yi = bi.t + pyi
        % zi = ci.t + pzi
        
        % With:
        % ai: n_MC(1,i) | bi: n_MC(2,i) | ci: n_MC(3,i)
        % pxi: point_MC(1,i) | pyi: point_MC(2,i) | pzi: point_MC(3,i)
        
        % To find the intersection between the plane (P) and the helical
        % axis, we need to solve: P(x(t),y(t),z(t)) = 0, then replace t in
        % the parametric equation to have the (x,y,z) coordinates of the
        % intersection point

% Scaphoid
for nr = 1:length(available)-1
    ai_Sca(nr) = n_Sca(1,nr);
    bi_Sca(nr) = n_Sca(2,nr);
    ci_Sca(nr) = n_Sca(3,nr);
    px_Sca(nr) = point_Sca(1,nr);
    py_Sca(nr) = point_Sca(2,nr);
    pz_Sca(nr) = point_Sca(3,nr);
    syms t
    solt_Sca(nr) = double(solve(a*((ai_Sca(nr))*t+px_Sca(nr))+...
        b*((bi_Sca(nr))*t+py_Sca(nr))+c*((ci_Sca(nr))*t+pz_Sca(nr))+d==0,t));
    Px_Sca(nr) = ai_Sca(nr)*solt_Sca(nr)+px_Sca(nr);
    Py_Sca(nr) = bi_Sca(nr)*solt_Sca(nr)+py_Sca(nr);
    Pz_Sca(nr) = ci_Sca(nr)*solt_Sca(nr)+pz_Sca(nr);
end
fig_Sca = figure;
Inter_Sca = scatter(Py_Sca,Pz_Sca,'fill','r');
grid on
yL_Sca = xlim;
zL_Sca = ylim;
hold on
plot(yL_Sca, [0 0], 'k-') % Show Y axis = 0 line
plot([0 0], zL_Sca, 'k-') % Show Z axis = 0 line
xlabel('Y axis');
ylabel('Z axis')
title('Scaphoid: Intersection of the helical axis with the YZ plane of the radius CS');

% Trapezium
for nr = 1:length(available)-1
    ai_Trap(nr) = n_Trap(1,nr);
    bi_Trap(nr) = n_Trap(2,nr);
    ci_Trap(nr) = n_Trap(3,nr);
    px_Trap(nr) = point_Trap(1,nr);
    py_Trap(nr) = point_Trap(2,nr);
    pz_Trap(nr) = point_Trap(3,nr);
    syms t
    solt_Trap(nr) = double(solve(a*((ai_Trap(nr))*t+px_Trap(nr))+...
        b*((bi_Trap(nr))*t+py_Trap(nr))+c*((ci_Trap(nr))*t+pz_Trap(nr))+d==0,t));
    Px_Trap(nr) = ai_Trap(nr)*solt_Trap(nr)+px_Trap(nr);
    Py_Trap(nr) = bi_Trap(nr)*solt_Trap(nr)+py_Trap(nr);
    Pz_Trap(nr) = ci_Trap(nr)*solt_Trap(nr)+pz_Trap(nr);
end
firg_Trap = figure;
Inter_Trap = scatter(Py_Trap,Pz_Trap,'fill','g');
grid on
yL_Trap = xlim;
zL_Trap = ylim;
hold on
plot(yL_Trap, [0 0], 'k-') % Show Y axis = 0 line
plot([0 0], zL_Trap, 'k-') % Show Z axis = 0 line
xlabel('Y axis');
ylabel('Z axis')
title('Trapezium: Intersection of the helical axis with the YZ plane of the radius CS');

% MC1
for nr = 1:length(available)-1
    ai_MC(nr) = n_MC(1,nr);
    bi_MC(nr) = n_MC(2,nr);
    ci_MC(nr) = n_MC(3,nr);
    px_MC(nr) = point_MC(1,nr);
    py_MC(nr) = point_MC(2,nr);
    pz_MC(nr) = point_MC(3,nr);
    syms t
    solt_MC(nr) = double(solve(a*((ai_MC(nr))*t+px_MC(nr))+...
        b*((bi_MC(nr))*t+py_MC(nr))+c*((ci_MC(nr))*t+pz_MC(nr))+d==0,t));
    Px_MC(nr) = ai_MC(nr)*solt_MC(nr)+px_MC(nr);
    Py_MC(nr) = bi_MC(nr)*solt_MC(nr)+py_MC(nr);
    Pz_MC(nr) = ci_MC(nr)*solt_MC(nr)+pz_MC(nr);
end
fig_MC = figure;
Inter_MC = scatter(Py_MC,Pz_MC,'fill','b');
grid on
yL_MC = xlim;
zL_MC = ylim;
hold on
plot(yL_MC, [0 0], 'k-') % Show Y axis = 0 line
plot([0 0], zL_MC, 'k-') % Show Z axis = 0 line
xlabel('Y axis');
ylabel('Z axis')
title('MC1: Intersection of the helical axis with the YZ plane of the radius CS');

% Calculation of the angle between th helical axis and the YZ plane of the
% radius CS:
    % We can calculate this angle from the coordinates of the vector normal
    % to the plane (here: the x-axis = (1,0,0)) and the vector defining the 
    % direction of the helical axis (n_bone)

for nr = 1:length(available)-1
    frame_nr(nr) = nr;
end

% Scaphoid
for nr = 1:length(available)-1
    alpha_Sca(nr) = asin(abs(n_Sca(1,nr))/(sqrt(n_Sca(1,nr)^2+n_Sca(2,nr)^2+n_Sca(3,nr)^2)))*180/pi;
end
figalpha_Sca = figure;
angle_Sca = scatter(frame_nr,alpha_Sca,'fill','r');
grid on
axL_Sca = xlim;
fyL_Sca = ylim;
xlim([0,axL_Sca(2)]);
ylim([0,round(fyL_Sca(2))]);
xlabel('frames');
ylabel('angle (deg)')
title('Scaphoid: Angle between helical axis and the YZ plane of the radius CS');

% Trapezium
for nr = 1:length(available)-1
    alpha_Trap(nr) = asin(abs(n_Trap(1,nr))/(sqrt(n_Trap(1,nr)^2+n_Trap(2,nr)^2+n_Trap(3,nr)^2)))*180/pi;
end
figalpha_Trap = figure;
angle_Trap = scatter(frame_nr,alpha_Trap,'fill','g');
grid on
axL_Trap = xlim;
fyL_Trap = ylim;
xlim([0,axL_Trap(2)]);
ylim([0,round(fyL_Trap(2))]);
xlabel('frames');
ylabel('angle (deg)')
title('Trapexium: Angle between helical axis and the YZ plane of the radius CS');

% MC1
for nr = 1:length(available)-1
    alpha_MC(nr) = asin(abs(n_MC(1,nr))/(sqrt(n_MC(1,nr)^2+n_MC(2,nr)^2+n_MC(3,nr)^2)))*180/pi;
end
figalpha_MC = figure;
angle_MC = scatter(frame_nr,alpha_MC,'fill','b');
grid on
axL_MC = xlim;
fyL_MC = ylim;
xlim([0,axL_MC(2)]);
ylim([0,round(fyL_MC(2))]);
xlabel('frames');
ylabel('angle (deg)')
title('MC1: Angle between helical axis and the YZ plane of the radius CS');

% Gather all graphs together

Final_Graph = figure;
subplot(2,3,1);
Inter_Sca = scatter(Py_Sca,Pz_Sca,'fill','r');
grid on
hold on
plot(yL_Sca, [0 0], 'k-') % Show Y axis = 0 line
plot([0 0], zL_Sca, 'k-') % Show Z axis = 0 line
xlabel('Y axis');
ylabel('Z axis')
title('Scaphoid: Intersection of the helical axis with the YZ plane of the radius CS');

subplot(2,3,2);
Inter_Trap = scatter(Py_Trap,Pz_Trap,'fill','g');
grid on
hold on
plot(yL_Trap, [0 0], 'k-') % Show Y axis = 0 line
plot([0 0], zL_Trap, 'k-') % Show Z axis = 0 line
xlabel('Y axis');
ylabel('Z axis')
title('Trapezium: Intersection of the helical axis with the YZ plane of the radius CS');

subplot(2,3,3);
Inter_MC = scatter(Py_MC,Pz_MC,'fill','b');
grid on
hold on
plot(yL_MC, [0 0], 'k-') % Show Y axis = 0 line
plot([0 0], zL_MC, 'k-') % Show Z axis = 0 line
xlabel('Y axis');
ylabel('Z axis')
title('MC1: Intersection of the helical axis with the YZ plane of the radius CS');

subplot(2,3,4);
angle_Sca = scatter(frame_nr,alpha_Sca,'fill','r');
grid on
xlim([0,axL_Sca(2)]);
ylim([0,round(fyL_Sca(2))]);
xlabel('frames');
ylabel('angle (deg)')
title('Scaphoid: Angle between helical axis and the YZ plane of the radius CS');

subplot(2,3,5);
angle_Trap = scatter(frame_nr,alpha_Trap,'fill','g');
grid on
xlim([0,axL_Trap(2)]);
ylim([0,round(fyL_Trap(2))]);
xlabel('frames');
ylabel('angle (deg)')
title('Trapexium: Angle between helical axis and the YZ plane of the radius CS');

subplot(2,3,6);
angle_MC = scatter(frame_nr,alpha_MC,'fill','b');
grid on
xlim([0,axL_MC(2)]);
ylim([0,round(fyL_MC(2))]);
xlabel('frames');
ylabel('angle (deg)')
title('MC1: Angle between helical axis and the YZ plane of the radius CS');