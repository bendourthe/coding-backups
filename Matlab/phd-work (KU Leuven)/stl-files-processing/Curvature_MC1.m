clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Curvature_MC1: Estimates the curvature of the articular surface of the
%                   base of the first metacarpal (MC1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dependencies:        
%               STL_ReadFile.m 
%                   TRI_RemoveInvalidTriangles.m
%                   TRI_RemoveBadlyConnectedTriangles.m
%               vectarrow.m
% Input: 
%           MC1: STL file corresponding to the first metacarpal
% Output:
%           Estimation od the curvature of the articular surface of the
%               base of MC1 according to three critical points manually
%               selected
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Selection of the current directory (where the STL files are)
dir = 'C:\Users\u0096636\Documents\Professional\PhD KU Leuven\Data\Study_Brown_Rig\Scan18\STL\';

% FILE
MC1 = 'N_Trap_S-1-0.4.stl';

% Read STL file
[F_MC1, V_MC1] =  STL_ReadFile([dir MC1],true);

% Plot the STL
stl = figure;
    [obj, li, ax] = GUI_PlotShells(stl, {F_MC1}, {V_MC1},...
            {ones(size(V_MC1,1),1)});
        
% Rotate the STL in order to have the best view for the points selection

rotate3d on
cprintf('Keywords','Rotate the figure')
cprintf('Errors',' WARNING: you have 5 seconds!')    % 
pause(5) % Gives 5 seconds to rotate the figure in the desired orientation
rotate3d off

% Select 3 relevant points that will estimate the curvature along the
% articular surface of MC1

    % Point 1
disp(' ')
cprintf('Keywords','Select the first point (P1 - dorsal beak)')
k1 = waitforbuttonpress; % k=0 if mouse pressed | k=1 if keyboard
if k1==0
    P = get(gca,'CurrentPoint');
    P1 = [P(1,1);P(1,2);P(1,3)];
    hold on
    plot3(P1(1),P1(2),P1(3),'Marker','o','MarkerSize',5,'Color','r','MarkerFaceColor','r')
end

    % Point 2
disp(' ')
cprintf('Keywords','Select the second point (P2 - volar beak)')
k2 = waitforbuttonpress;
if k2==0
    P = get(gca,'CurrentPoint');
    P2 = [P(1,1);P(1,2);P(1,3)];
    hold on
    plot3(P2(1),P2(2),P2(3),'Marker','o','MarkerSize',5,'Color','r','MarkerFaceColor','r')
end

    % Point 3
disp(' ')
cprintf('Keywords','Select the third point (between P1 and P2, where it seems the deepest)')
k3 = waitforbuttonpress;
if k3==0
    P = get(gca,'CurrentPoint');
    P3 = [P(1,1);P(1,2);P(1,3)];
    hold on
    plot3(P3(1),P3(2),P3(3),'Marker','o','MarkerSize',5,'Color','r','MarkerFaceColor','r')
    pause(1);
else
end
close all

% Calculation of the center of the circle passing by our 3 selected points

    % To find the equation of a circle passing by 3 points, we need 4
    % equations:
       % (1) (x1-x0)^2+(y1-y0)^2+(z1-z0)^2 = r0^2
       % (2) (x2-x0)^2+(y2-y0)^2+(z2-z0)^2 = r0^2
       % (3) (x3-x0)^2+(y3-y0)^2+(z3-z0)^2 = r0^2
    % With:
       % P1=(x1,y1,z1) | P2=(x2,y2,z2) | P3=(x3,y3,z3)
       % C=(x0,y0,z0) -> center of the circle and r0 its radius
    % For the 4th equation, we can say that the point C belongs to the
    % plane defined by our 3 previous points, hence:
       % Equation of the plane:
          % P1P2=[x2-x1;y2-y1;z2-z1]
          % P1P3=[x3-x1;y3-y1;z3-z1]
          % System:
             % x = x1+(x2-x1)*l+(x3-x1)u
             % y = y1+(y2-y1)*l+(y3-y1)u
             % z = z1+(z2-z1)*l+(z3-z1)u
          % Matrix:
             % P = [(x-x1) (x2-x1) (x3-x1);(y-y1) (y2-y1) (y3-y1);(z-z1)
             %     (z2-z1) (z3-z1)]
       % (4) det(P) = 0
    % To obtain the equation of our circle and hence, its radius (used to
    % estimate our curvature), we need to solve this system:
    
x1=P1(1);
y1=P1(2);
z1=P1(3);

x2=P2(1);
y2=P2(2);
z2=P2(3);

x3=P3(1);
y3=P3(2);
z3=P3(3);

syms x0 y0 z0 r0
S = solve((x1-x0)^2+(y1-y0)^2+(z1-z0)^2 == r0^2,(x2-x0)^2+(y2-y0)^2+(z2-z0)^2 == r0^2,(x3-x0)^2+(y3-y0)^2+(z3-z0)^2 == r0^2,(x0-x1)*((y2-y1)*(z3-z1)-(z2-z1)*(y3-y1))-(y0-y1)*((x2-x1)*(z3-z1)-(z2-z1)*(x3-x1))+(z0-z1)*((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1))==0);
S = [S.x0 S.y0 S.z0 S.r0];

x0=S(1);
y0=S(2);
z0=S(3);
xo=vpa(x0,4);
yo=vpa(y0,4);
zo=vpa(z0,4);
C=[xo yo zo];
disp(' ')
radius=abs(vpa(S(4),4))

% Estimation of the curvature of the base of MC1 according to the selection
% of 3 critical points

Curvature=vpa(1/radius,4)