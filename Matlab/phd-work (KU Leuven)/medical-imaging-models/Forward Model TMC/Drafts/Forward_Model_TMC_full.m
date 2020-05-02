clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward_Model_TMC:
%       This model is designed to estimate the mechanical response of the
% articular cartilage and the main ligaments stabilizing the TMC joint when
% a specific joint loading is applied to the MC1. The aims of this model 
% are (1) to understand how muscles control the kinematics of the thumb and
% how joint loading can affect its motion, (2) to provide a new simulation 
% tool to understand the occurrence mechanism of osteoarthritis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dependencies:        
%               STL_ReadFile.m
%                   TRI_RemoveInvalidTriangles.m
%                   TRI_RemoveBadlyConnectedTriangles.m
%               PlacePoints3
%               GUI_PlotShells
%               GUI_VisualisationUI
%               polyhedronCentroid
% Input: 
%           Trap: STL file corresponding to the trapezium
%           MC1: STL file corresponding to the first metacarpal
%           Cartilage properties
%           Ligement properties
%           Orientation and position of MC1 (Trap attached to CS origin)
%           Origin and Insertion points of ligaments
%           Point of application and direction of the resultant force
% 
% Output:
%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Selection of the current directory (where the STL files are)
dir = 'C:\Users\u0096636\Documents\Professional\PhD KU Leuven\International Mobility\USC\Forward Model\STL\';
trap_cs_available = 1;          % loads data if available (=1), else, generate data
    filenameCS_trap = 'trap_cs.txt';
    fileCS_trap = [dir filenameCS_trap];
mc1_cs_available = 1;           % loads data if available (=1), else, generate data
    filenameCS_mc1 = 'mc1_cs.txt';
    fileCS_mc1 = [dir filenameCS_mc1];
origin_lig_available = 1;       % loads data if available (=1), else, generate data
    filenameO_lig = 'origin_lig.txt';
    fileO_lig = [dir filenameO_lig];
insertion_lig_available = 1;    % loads data if available (=1), else, generate data
    filenameI_lig = 'insertion_lig.txt';
    fileI_lig = [dir filenameI_lig];
loading_available = 1;    % loads data if available (=1), else, generate data
    filename_loading = 'loading.txt';
    fileloading = [dir filename_loading];
    
% FILES
Trap = 'SCAN3_N_trap.stl';
MC1 = 'SCAN3_N_mc1.stl';
[F_trap, V_trap] =  STL_ReadFile([dir Trap],true);
[F_mc1, V_mc1] =  STL_ReadFile([dir MC1],true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL PARAMETERS

    % Time frame
tmax = 10; % in ms (equivalent to number of iteration)

    % Bone properties
        %MC1
rho = 500; % Density (in kg/m3)
M_mc1 = 0.002; % in kg
% Principal moments of inertia (1: x; 2: y; 3: z), in kg.mm2
I1 = 0.04241;
I2 = 0.413;
I3 = 0.4198;
   
    % Position & Orientation of MC1
x_translation = 0;
y_translation = 0;
z_translation = 0;
x_angle = 0;
y_angle = 0;
z_angle = 0;

    % Ligament properties
CSA_SAOL = 27.1; % average cross-sectional area of SAOL (in mm2)
CSA_dAOL = 10.9; % average cross-sectional area of dAOL (in mm2)
CSA_DRL = 14.9; % average cross-sectional area of DRL (in mm2)
CSA_POL = 29.6; % average cross-sectional area of POL (in mm2)
CSA_UCL = 8.9; % average cross-sectional area of UCL (in mm2)
C1 = 1; % in MPa
C2 = 1; % in MPa
C3 = 33; % in MPa
nu1 = 1; % in MPa.s
nu2 = 22; % in MPa.s

    % Cartilage properties
HA = 0.80;            % Average Aggregate Modulus for a human cartilage (MPa)
d0 = 0.90;            % Initial solid content (material constant)
a0 = 0.20;            % Average fluid-to-solid true density ratio
d1 = d0*(a0 + 1)/(a0 + d0); % Material constant
Ttrap = 0.8;          % Average trapezium cartilage thickness (mm)
Tmc1 = 0.7;           % Average first metacarpal cartilage thickness (mm)
Ttot = Ttrap + Tmc1;  % Total joint cartilage thickness (mm)

    % Loading conditions    
        % Define the magnitude of the force vector:
magnitude_x = 0; % in N
magnitude_y = 10;
magnitude_z = 0;
        % Define the orientation of the force
resultant_x_angle = 0; % in degrees
resultant_y_angle = 0; % in degrees
resultant_z_angle = 0; % in degrees

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COORDINATE SYSTEMS AND ALIGNMENT

% CS on trapezium:
%   origin: midway between the two landmarks defining the y-axis.
%   y: extends from the exact mid-point of the central ridge of the 
%       trapezial saddle to the center of the junction of the trapezium,
%       scaphoid and trapezoid.
%   x: runs in a dorsal-to-volar direction along a line perpendicular to
%       the central ridge of the trapezium and passes through the mid-point
%       of the dorsal surface to the proximal volar pole of the tubercle of
%       the trapezium.
%   z: is perpendicular to the X and Y axes and nearly parallel to the
%       central ridge of the trapezial metacarpal surface.

if(trap_cs_available)
    A = txt2mat(fileCS_trap);
    
    Y1 = A(1,:);    % mid-point of the central ridge of the trapezial saddle
    Y2 = A(2,:);    % center of the junction of the trapezium, scaphoid and trapezoid
    X = A(3,:);     % center of the distal articular surface, facing the center of the base of MC1
else
    [ad_curve] = PlacePoints3({F_trap}, {V_trap}, {F_trap}, {V_trap}, 'Select 3 landmarks (Y1, Y2 and X, in that order!). End by pressing enter and exiting the figure');
    close all;
    % Retrieves coordinates of the points:
    Y1 = ad_curve.points{end,1}(1,:);
    Y2 = ad_curve.points{end,1}(2,:);
    X = ad_curve.points{end,1}(3,:);
    
    % Save in a txt file
    fid = fopen(fileCS_trap,'w+');
    fprintf(fid,'%5.8f\t%5.8f\t%5.8f\r\n',[Y1,Y2,X]);
    fclose(fid);
end

O_trap = (Y1+Y2)/2;
Yaxis_trap = (Y2-Y1)/norm(Y2-Y1);
Xaxis_trap = (X-O_trap)/norm(X-O_trap);
Zaxis_trap = cross(Xaxis_trap,Yaxis_trap)/norm(cross(Xaxis_trap,Yaxis_trap));
Xaxis_trap = cross(Yaxis_trap,Zaxis_trap)/norm(cross(Yaxis_trap,Zaxis_trap));
cs_trap = [Xaxis_trap; Yaxis_trap; Zaxis_trap];
% Converges to trapezium CS
RotMat_trap = [Xaxis_trap', Yaxis_trap', Zaxis_trap'];
loc2glob_trap = [RotMat_trap, O_trap';
     0 0 0 1]; 
   % These tranformation matrices convert from local to global CS.
   % Invert them to go from global to local.
glob2loc_trap = [RotMat_trap', -(RotMat_trap')*O_trap'; 0 0 0 1];
R_glob2loc_trap = glob2loc_trap(1:3,1:3);
T_glob2loc_trap = glob2loc_trap(1:3,4);

% Transformes data to trapezium CS
V_trap = R_glob2loc_trap*V_trap.' + repmat(T_glob2loc_trap,1,size(V_trap,1));
V_trap = V_trap.';

% Plots the trapezium with the new coordinate system for validation
    f1 = figure;
    [obj, li, ax] = GUI_PlotShells(f1, {F_trap},{V_trap},...
        {ones(size(V_trap,1),1)});
    GUI_VisualisationUI(ancestor(ax, 'figure'), true, ax, true, true);
    box off
        
        % Plots the coordinate system and its origin
    hold on
    arrow([0;0;0],[15;0;0],10,70,30,5,'EdgeColor','r','FaceColor','r');
    arrow([0;0;0],[0;15;0],10,70,30,5,'EdgeColor','g','FaceColor','g');
    arrow([0;0;0],[0;0;15],10,70,30,5,'EdgeColor','b','FaceColor','b');
    view([180,90,45]);
    
        % Plots 2 buttons: 'Bad CS' if the coordinate system looks wrong,
        % this will stop the code, and 'Continue' if the coordinate system
        % looks right, will continue to run the rest of the code
    b1 = uicontrol('Position',[20 20 200 40],'String','Continue',...
              'Callback','uiresume(gcbf)');
    b2 = uicontrol('Position',[20 60 200 40],'String','Bad CS',...
              'Callback','close(gcbf)');
    uiwait(gcf);
    close(f1);

% CS on MC1:
%   origin: midway between the base and head of each metacarpal.
%   x: line parallel to a line from the center of the distal head of the 
%       metacarpal to the midpoint of the base of the metacarpal.
%   y: the x- and y-axis will form a sagittal plane that splits the
%       metacarpal into mirror images.
%   z: the common line perpendicular to the x- and y-axis.

if(mc1_cs_available)
    B = txt2mat(fileCS_mc1);
    
    X1 = B(1,:);    % center of the distal head of MC1
    X2 = B(2,:);    % center of the base of MC1
    Y = B(3,:);     % tip of the palmar beak
    % Convert mc1 into trapezium CS
    V_mc1 = R_glob2loc_trap*V_mc1.' + repmat(T_glob2loc_trap,1,size(V_mc1,1));
    V_mc1 = V_mc1.';
else
    % Convert mc1 into trapezium CS
    V_mc1 = R_glob2loc_trap*V_mc1.' + repmat(T_glob2loc_trap,1,size(V_mc1,1));
    V_mc1 = V_mc1.';
    [ad_curve] = PlacePoints3({F_mc1}, {V_mc1}, {F_mc1}, {V_mc1}, 'Select 3 landmarks (X1, X2 and Y, in that order!). End by pressing enter and exiting the figure');
    close all;
    % Retrieves coordinates of the points:
    X1 = ad_curve.points{end,1}(1,:);
    X2 = ad_curve.points{end,1}(2,:);
    Y = ad_curve.points{end,1}(3,:);
    
    % Save in a txt file
    fid = fopen(fileCS_mc1,'w+');
    fprintf(fid,'%5.8f\t%5.8f\t%5.8f\r\n',[X1,X2,Y]);
    fclose(fid);
end

O_mc1 = (X1+X2)/2;
Xaxis_mc1 = (X2-X1)/norm(X2-X1);
Yaxis_mc1 = (Y-O_mc1)/norm(Y-O_mc1);
Zaxis_mc1 = cross(Xaxis_mc1,Yaxis_mc1)/norm(cross(Xaxis_mc1,Yaxis_mc1));
Yaxis_mc1 = cross(Zaxis_mc1,Xaxis_mc1)/norm(cross(Zaxis_mc1,Xaxis_mc1));
cs_mc1 = [Xaxis_mc1; Yaxis_mc1; Zaxis_mc1];
% Converges to mc1 CS
RotMat_mc1 = [Xaxis_mc1', Yaxis_mc1', Zaxis_mc1'];
loc2glob_mc1 = [RotMat_mc1, O_mc1';
     0 0 0 1]; 
   % These tranformation matrices convert from local to global CS.
   % Invert them to go from global to local.
glob2loc_mc1 = [RotMat_mc1', -(RotMat_mc1')*O_mc1'; 0 0 0 1];
R_glob2loc_mc1 = glob2loc_mc1(1:3,1:3);
T_glob2loc_mc1 = glob2loc_mc1(1:3,4);

% Transformes data to mc1 CS for CS visualization
V0_mc1 = R_glob2loc_mc1*V_mc1.' + repmat(T_glob2loc_mc1,1,size(V_mc1,1));
V0_mc1 = V0_mc1.';

% Plots the mc1 with the new coordinate system for validation
    f2 = figure;
    [obj, li, ax] = GUI_PlotShells(f2, {F_mc1},{V0_mc1},...
        {ones(size(V0_mc1,1),1)});
    box off
    GUI_VisualisationUI(ancestor(ax, 'figure'), true, ax, true, true);
        
        % Plots the coordinate system and its origin
    hold on
    arrow([0;0;0],[40;0;0],10,70,30,5,'EdgeColor','r','FaceColor','r');
    arrow([0;0;0],[0;15;0],10,70,30,5,'EdgeColor','g','FaceColor','g');
    arrow([0;0;0],[0;0;15],10,70,30,5,'EdgeColor','b','FaceColor','b');
    view([180,90,45]);
        % Plots 2 buttons: 'Bad CS' if the coordinate system looks wrong,
        % this will stop the code, and 'Continue' if the coordinate system
        % looks right, will continue to run the rest of the code
    b1 = uicontrol('Position',[20 20 200 40],'String','Continue',...
              'Callback','uiresume(gcbf)');
    b2 = uicontrol('Position',[20 60 200 40],'String','Bad CS',...
              'Callback','close(gcbf)');
    uiwait(gcf);
    close(f2);

% Plots the trapezium with the new coordinate system
    f1 = figure;
    [obj, li, ax] = GUI_PlotShells(f1, {F_trap},{V_trap},...
        {ones(size(V_trap,1),1)});
    box off
    view([180,90,45]);
        % Plots the coordinate system and its origin
    hold on
    arrow([0;0;0],[15;0;0],10,70,30,5,'EdgeColor','r','FaceColor','r');
    arrow([0;0;0],[0;15;0],10,70,30,5,'EdgeColor','g','FaceColor','g');
    arrow([0;0;0],[0;0;15],10,70,30,5,'EdgeColor','b','FaceColor','b');
    set(gcf,'numbertitle','off','name','Coordinate system of the trapezium');

% Plots the mc1 with the new coordinate system for validation
    f2 = figure;
    [obj, li, ax] = GUI_PlotShells(f2, {F_mc1},{V0_mc1},...
        {ones(size(V0_mc1,1),1)});
    box off
    view([180,90,45]);
        % Plots the coordinate system and its origin
    hold on
    arrow([0;0;0],[40;0;0],10,70,30,5,'EdgeColor','r','FaceColor','r');
    arrow([0;0;0],[0;15;0],10,70,30,5,'EdgeColor','g','FaceColor','g');
    arrow([0;0;0],[0;0;15],10,70,30,5,'EdgeColor','b','FaceColor','b');
    set(gcf,'numbertitle','off','name','Coordinate system of the MC1');
    
% Plots the two bones in the trapezium CS
f3 = figure;
    [obj, li, ax] = GUI_PlotShells(f3, {F_trap;F_mc1},{V_trap;V_mc1},...
        {ones(size(V_trap,1),1);ones(size(V_mc1,1),1)});
    box off
    view([45,45,45]);
    % Plots the trapezium CS
    hold on
    arrow([0;0;0],[15;0;0],10,70,30,5,'EdgeColor','r','FaceColor','r');
    arrow([0;0;0],[0;15;0],10,70,30,5,'EdgeColor','g','FaceColor','g');
    arrow([0;0;0],[0;0;15],10,70,30,5,'EdgeColor','b','FaceColor','b');
    set(gcf,'numbertitle','off','name','Original system');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ORIENTATION OF MC1

    % Define the translation matrix
    Tr_mc1 = [x_translation; y_translation; z_translation];
    % Define rotation matrix
    a = x_angle*pi/180;
    b = y_angle*pi/180;
    c = z_angle*pi/180;

    Rot_mc1 = [cos(c)*cos(b), cos(c)*sin(b)*sin(a)-sin(c)*cos(a), cos(c)*sin(b)*cos(a)+sin(c)*sin(a);...
        sin(c)*cos(b), sin(c)*sin(b)*sin(a)+cos(c)*cos(a), sin(c)*sin(b)*cos(a)-cos(c)*sin(a);
        -sin(b), cos(b)*sin(a), cos(b)*cos(a)];

    % Change the orientation of MC1 according to the rotations defined
    Vor_mc1 = Rot_mc1*V_mc1.' + repmat(Tr_mc1,1,size(V_mc1,1));
    Vor_mc1 = Vor_mc1.';

% Plots the two bones with the new position and orientation of mc1
f4 = figure;
    [obj, li, ax] = GUI_PlotShells(f4, {F_trap;F_mc1},{V_trap;Vor_mc1},...
        {ones(size(V_trap,1),1);ones(size(Vor_mc1,1),1)});
    box off
    view([45,45,45]);
    % Plots the trapezium CS
    hold on
    arrow([0;0;0],[15;0;0],10,70,30,5,'EdgeColor','r','FaceColor','r');
    arrow([0;0;0],[0;15;0],10,70,30,5,'EdgeColor','g','FaceColor','g');
    arrow([0;0;0],[0;0;15],10,70,30,5,'EdgeColor','b','FaceColor','b');
    set(gcf,'numbertitle','off','name','System after modifying MC1 orientation');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LIGAMENTS ORIGIN

if(origin_lig_available)
    C = txt2mat(fileO_lig);
    
    O_SAOL = C(1,:);
    O_dAOL = C(2,:);
    O_DRL = C(3,:);
    O_POL = C(4,:);
    O_UCL = C(5,:);
    idxDO1 = C(6,1);
    idxDO2 = C(7,1);
    idxDO3 = C(8,1);
    idxDO4 = C(9,1);
    idxDO5 = C(10,1);
    O_SAOL_i = V_trap(C(6,1),:);
    O_dAOL_i = V_trap(C(7,1),:);
    O_DRL_i = V_trap(C(8,1),:);
    O_POL_i = V_trap(C(9,1),:);
    O_UCL_i = V_trap(C(10,1),:);
else
    [ad_curve] = PlacePoints3({F_trap}, {V_trap}, {F_trap}, {V_trap}, 'Select the origin of the SAOL, dAOL, DRL, POL, UCL, IN THIS ORDER. End by pressing enter and exiting the figure');
    close;
    % Retrieves coordinates of the points:
    O_SAOL = ad_curve.points{end,1}(1,:);
    O_dAOL = ad_curve.points{end,1}(2,:);
    O_DRL = ad_curve.points{end,1}(3,:);
    O_POL = ad_curve.points{end,1}(4,:);
    O_UCL = ad_curve.points{end,1}(5,:);
    
    % Find the closest vertice
    [idxO1,DO1] = knnsearch(O_SAOL,V_trap);
    [CO1,idxDO1] = min(DO1);
    [idxO2,DO2] = knnsearch(O_dAOL,V_trap);
    [CO2,idxDO2] = min(DO2);
    [idxO3,DO3] = knnsearch(O_DRL,V_trap);
    [CO3,idxDO3] = min(DO3);
    [idxO4,DO4] = knnsearch(O_POL,V_trap);
    [CO4,idxDO4] = min(DO4);
    [idxO5,DO5] = knnsearch(O_UCL,V_trap);
    [CO5,idxDO5] = min(DO5);
    
    O_SAOL_i = V_trap(idxDO1,:);
    O_dAOL_i = V_trap(idxDO2,:);
    O_DRL_i = V_trap(idxDO3,:);
    O_POL_i = V_trap(idxDO4,:);
    O_UCL_i = V_trap(idxDO5,:);
    
    idxO = [idxDO1 0 0; idxDO2 0 0; idxDO3 0 0; idxDO4 0 0; idxDO5 0 0];
    
    % Save in a txt file
    fid = fopen(fileO_lig,'w+');
    fprintf(fid,'%5.8f\t%5.8f\t%5.8f\r\n',[O_SAOL,O_dAOL,O_DRL,O_POL,O_UCL,...
        idxO(1,:),idxO(2,:),idxO(3,:),idxO(4,:),idxO(5,:)]);
    fclose(fid);
end
OPs = [O_SAOL_i;O_dAOL_i;O_DRL_i;O_POL_i;O_UCL_i];

% LIGAMENTS INSERTION

if(insertion_lig_available)
    D = txt2mat(fileI_lig);
    
    I_SAOL = D(1,:);
    I_dAOL = D(2,:);
    I_DRL = D(3,:);
    I_POL = D(4,:);
    I_UCL = D(5,:);
    idxDI1 = D(6,1);
    idxDI2 = D(7,1);
    idxDI3 = D(8,1);
    idxDI4 = D(9,1);
    idxDI5 = D(10,1);
    I_SAOL_i = Vor_mc1(D(6,1),:);
    I_dAOL_i = Vor_mc1(D(7,1),:);
    I_DRL_i = Vor_mc1(D(8,1),:);
    I_POL_i = Vor_mc1(D(9,1),:);
    I_UCL_i = Vor_mc1(D(10,1),:);

else
    [ad_curve] = PlacePoints3({F_mc1}, {Vor_mc1}, {F_mc1}, {Vor_mc1}, 'Select the insertion of the SAOL, dAOL, DRL, POL, UCL, IN THIS ORDER. End by pressing enter and exiting the figure');
    close;
    % Retrieves coordinates of the points:
    I_SAOL = ad_curve.points{end,1}(1,:);
    I_dAOL = ad_curve.points{end,1}(2,:);
    I_DRL = ad_curve.points{end,1}(3,:);
    I_POL = ad_curve.points{end,1}(4,:);
    I_UCL = ad_curve.points{end,1}(5,:);
    
    % Find the closest vertice
    [idxI1,DI1] = knnsearch(I_SAOL,Vor_mc1);
    [CI1,idxDI1] = min(DI1);
    [idxI2,DI2] = knnsearch(I_dAOL,Vor_mc1);
    [CI2,idxDI2] = min(DI2);
    [idxI3,DI3] = knnsearch(I_DRL,Vor_mc1);
    [CI3,idxDI3] = min(DI3);
    [idxI4,DI4] = knnsearch(I_POL,Vor_mc1);
    [CI4,idxDI4] = min(DI4);
    [idxI5,DI5] = knnsearch(I_UCL,Vor_mc1);
    [CI5,idxDI5] = min(DI5);
    
    I_SAOL_i = Vor_mc1(idxDI1,:);
    I_dAOL_i = Vor_mc1(idxDI2,:);
    I_DRL_i = Vor_mc1(idxDI3,:);
    I_POL_i = Vor_mc1(idxDI4,:);
    I_UCL_i = Vor_mc1(idxDI5,:);
    
    idxI = [idxDI1 0 0; idxDI2 0 0; idxDI3 0 0; idxDI4 0 0; idxDI5 0 0];
    
    % Save in a txt file
    fid = fopen(fileI_lig,'w+');
    fprintf(fid,'%5.8f\t%5.8f\t%5.8f\r\n',[I_SAOL,I_dAOL,I_DRL,I_POL,I_UCL,...
        idxI(1,:),idxI(2,:),idxI(3,:),idxI(4,:),idxI(5,:)]);
    fclose(fid);
end
IPs = [I_SAOL_i;I_dAOL_i;I_DRL_i;I_POL_i;I_UCL_i];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% JOINT LOADING

% Choose the point of application of the Resultant force
if(loading_available)
    E = txt2mat(fileloading);
    
    PointA = E(1,:);
    idxDPA = E(2,1);
    PointA_idx = Vor_mc1(idxDPA,:);
else
    [ad_curve] = PlacePoints3({F_mc1}, {Vor_mc1}, {F_mc1}, {Vor_mc1}, 'Select the point of application of the Resultant force. End by pressing enter and exiting the figure');
    close;

    % Retrieves coordinates of the points:
    PointA = ad_curve.points{end,1}(1,:);
    
    % Find the closest vertice
    [idxPA,DPA] = knnsearch(PointA,Vor_mc1);
    [CPA,idxDPA] = min(DPA);
    
    PointA_idx = Vor_mc1(idxDPA,:);
    idxPA = [idxDPA 0 0];
    
    % Save in a txt file
    fid = fopen(fileloading,'w+');
    fprintf(fid,'%5.8f\t%5.8f\t%5.8f\r\n',[PointA,idxPA]);
    fclose(fid);
end
    
    % Default position of the second point:
    Vorigin = [PointA_idx(1)-magnitude_x;PointA_idx(2)-magnitude_y;PointA_idx(3)-magnitude_z];
    
    % Define the orientation of the force
    a1 = resultant_x_angle*pi/180;
    b1 = resultant_y_angle*pi/180;
    c1 = resultant_z_angle*pi/180;

Rot_resultant = [cos(c1)*cos(b1), cos(c1)*sin(b1)*sin(a1)-sin(c1)*cos(a1), cos(c1)*sin(b1)*cos(a1)+sin(c1)*sin(a1);...
    sin(c1)*cos(b1), sin(c1)*sin(b1)*sin(a1)+cos(c1)*cos(a1), sin(c1)*sin(b1)*cos(a1)-cos(c1)*sin(a1);
    -sin(b1), cos(b1)*sin(a1), cos(b1)*cos(a1)];

    % Change the default orientation:    
    Vorigin1 = (Rot_resultant*Vorigin).';
    
    % Definition of the resultant force as a vector    
    FR = [PointA_idx(1)-Vorigin1(1);PointA_idx(2)-Vorigin1(2);PointA_idx(3)-Vorigin1(3)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FINAL SYSTEM

    % Plots the final system with the joint CS, the different ligaments and
    % the resultant force vector in the chosen orientation:
    
    f6 = figure;
    [obj, li, ax] = GUI_PlotShells(f6, {F_trap;F_mc1},{V_trap;Vor_mc1},...
        {ones(size(V_trap,1),1);ones(size(Vor_mc1,1),1)});
    box off
    view([45,45,45]);
    hold on
    arrow([Vorigin1(1);Vorigin1(2);Vorigin1(3)],[PointA_idx(1);PointA_idx(2);PointA_idx(3)]...
        ,1,70,30,5,'EdgeColor','r','FaceColor','r');
    plot3(PointA_idx(1),PointA_idx(2),PointA_idx(3),'Marker','o','MarkerSize',10,'Color',...
        'k','MarkerFaceColor','k');
    hold on
    for i=1:5;
        arrow([OPs(i,1);OPs(i,2);OPs(i,3)],[IPs(i,1);IPs(i,2);IPs(i,3)],...
        10,70,30,5,'EdgeColor',[1 0.5 0.2],'FaceColor',[1 0.5 0.2])
        plot3(OPs(i,1),OPs(i,2),OPs(i,3),'Marker','o','MarkerSize',10,'Color',...
        'r','MarkerFaceColor','r');
        plot3(IPs(i,1),IPs(i,2),IPs(i,3),'Marker','o','MarkerSize',10,'Color',...
        'b','MarkerFaceColor','b');
    end
    set(gcf,'numbertitle','off','name','System in its INITIAL configuration');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATION

% INITIAL CONDITIONS

    % Calculation of the shortest distance between the two bones selected
    [idx1,D1] = knnsearch(V_trap,Vor_mc1);
    [minD1,idxD1] = min(D1);
    Shortest_Distance = minD1;
    
    % Definition of the point of contact (CP) which will be our center of
    % rotation at each time step
    CP = V_trap(idxD1,:);

    % Transposition of inertia from center of gravity to point of contact
    CG = polyhedronCentroid(Vor_mc1,F_mc1); % center of gravity of MC1
    dcg_cr = sqrt((CP(1)-CG(1))^2+(CP(2)-CG(2))^2+(CP(3)-CG(3))^2); % distance from CG to CR
    I1cr = I1 - M_mc1*dcg_cr;
    I2cr = I2 - M_mc1*dcg_cr;
    I3cr = I3 - M_mc1*dcg_cr;
    
    % Calculation of the resultant torque (N.mm)
    MA = [PointA_idx(1) - CP(1);PointA_idx(2) - CP(2); PointA_idx(3) - CP(3)]; % Moment arm
    TR = [MA(2)*FR(3) - MA(3)*FR(2);MA(3)*FR(1) - MA(1)*FR(3);MA(1)*FR(2) - MA(2)*FR(1)];
    
    % Calculation of the initial length of the ligaments
    linit_SAOL = norm(I_SAOL_i-O_SAOL_i);
    linit_dAOL = norm(I_dAOL_i-O_dAOL_i);
    linit_DRL = norm(I_DRL_i-O_DRL_i);
    linit_POL = norm(I_POL_i-O_POL_i);
    linit_UCL = norm(I_UCL_i-O_UCL_i);
    
    % Calculation of the initial contact stress tensor (cartilage)
        % Calculation of the matrix Di which contains the Euclidean distances between
        % each vertices of each bone
        Di = pdist2(V_trap,Vor_mc1);
        
        % Save the indexes of the vertices of each bone which are located
        % in the deformation zone
        [Di_trap,Di_mc1] = find(Di<Ttot);
        idxDef = [Di_trap,Di_mc1];
        idxDef(:,any(idxDef==0,1)) = []; % Removes all the zero in the idxD Matrix
        idx_trap = unique(idxDef(:,1));
        idx_mc1 = unique(idxDef(:,2));

        % We can now test to which triangle these vertices correspond in
        % the F matrix and save their indexes
        test_F = ismember(F_trap,idx_trap);
        idx_F_trap = find(test_F(:,1)==1 & test_F(:,2)==1 & test_F(:,3)==1);
            % This idx_F_trap matrix is a Nx1 matrix, where N is the amount
            % of triangles located in the deformation zone, and each value
            % in this matrix corresponds to the index of a line in the F
            % matrix that defines one of these triangles
        for i=1:length(idx_F_trap)
            CA_F_trap(i,:) = F_trap(idx_F_trap(i),:);
        end
    
                % Each triangle is composed of 3 points called a,b and c
Xa1 = V_trap(CA_F_trap(:,1),1);
Ya1 = V_trap(CA_F_trap(:,1),2);
Za1 = V_trap(CA_F_trap(:,1),3);

Xb1 = V_trap(CA_F_trap(:,2),1);
Yb1 = V_trap(CA_F_trap(:,2),2);
Zb1 = V_trap(CA_F_trap(:,2),3);

Xc1 = V_trap(CA_F_trap(:,3),1);
Yc1 = V_trap(CA_F_trap(:,3),2);
Zc1 = V_trap(CA_F_trap(:,3),3);

        % Now we need to calculate the coordinates of the 2 vectors ab and ac
for i=1:length(CA_F_trap)
    Xab1(i) = Xb1(i)-Xa1(i);
    Yab1(i) = Yb1(i)-Ya1(i);
    Zab1(i) = Zb1(i)-Za1(i);
    Xac1(i) = Xc1(i)-Xa1(i);
    Yac1(i) = Yc1(i)-Ya1(i);
    Zac1(i) = Zc1(i)-Za1(i);
end

        % Now we need to calculate the surface area of each triangle and sum
        % all of them to obtain to total contact area
for i = 1:length(CA_F_trap)
    S_F_trap(i) = (1/2)*((Yab1(i)*Zac1(i)-Zab1(i)*Yac1(i))^2+(Zab1(i)*Xac1(i)-Xab1(i)*Zac1(i))^2+(Xab1(i)*Yac1(i)-Yab1(i)*Xac1(i))^2)^(1/2);
end
S_F_trap = S_F_trap.';
PCA_trap = sum(S_F_trap(:));

        % We need now to calculate the corresponding stress value for each
        % vertex which appears in the F matrix, so that we can calculate
        % the average stress value for the whole triangle and convert it
        % into force

for i=1:length(CA_F_trap)
    % for each point P1, P2 and P3 of the F matrix of the trapezium, we
    % need to know the index of the closest point in the V matrix of the
    % MC1
    Di_trap_mc1_P1 = pdist2(V_trap(CA_F_trap(i,1),:),Vor_mc1);
    [minDiP1(i),idxP1(i)] = min(Di_trap_mc1_P1);
    Di_trap_mc1_P2 = pdist2(V_trap(CA_F_trap(i,2),:),Vor_mc1);
    [minDiP2(i),idxP2(i)] = min(Di_trap_mc1_P2);
    Di_trap_mc1_P3 = pdist2(V_trap(CA_F_trap(i,3),:),Vor_mc1);
    [minDiP3(i),idxP3(i)] = min(Di_trap_mc1_P3);
    % distance between each vertices for P1, P2 and P3 of the F matrix
    DiP1(i) = ((V_trap(CA_F_trap(i,1),1)-Vor_mc1(idxP1(i),1)).^2+...
        (V_trap(CA_F_trap(i,1),2)-Vor_mc1(idxP1(i),2)).^2+...
        (V_trap(CA_F_trap(i,1),3)-Vor_mc1(idxP1(i),3)).^2).^(1/2);
    DiP2(i) = ((V_trap(CA_F_trap(i,2),1)-Vor_mc1(idxP2(i),1)).^2+...
        (V_trap(CA_F_trap(i,2),2)-Vor_mc1(idxP2(i),2)).^2+...
        (V_trap(CA_F_trap(i,2),3)-Vor_mc1(idxP2(i),3)).^2).^(1/2);
    DiP3(i) = ((V_trap(CA_F_trap(i,3),1)-Vor_mc1(idxP3(i),1)).^2+...
        (V_trap(CA_F_trap(i,3),2)-Vor_mc1(idxP3(i),2)).^2+...
        (V_trap(CA_F_trap(i,3),3)-Vor_mc1(idxP3(i),3)).^2).^(1/2);
    % deformation in MC1 cartilage layer for P1, P2 and P3 of the F matrix
    Defmc1_P1(i) = (Ttot-DiP1(i))/2/Tmc1+1;
    Defmc1_P2(i) = (Ttot-DiP2(i))/2/Tmc1+1;
    Defmc1_P3(i) = (Ttot-DiP3(i))/2/Tmc1+1;
    % deformation in MC1 cartilage layer for P1, P2 and P3 of the F matrix
    Deftrap_P1(i) = (Ttot-DiP1(i))/2/Ttrap+1;
    Deftrap_P2(i) = (Ttot-DiP2(i))/2/Ttrap+1;
    Deftrap_P3(i) = (Ttot-DiP3(i))/2/Ttrap+1;
    % stress at P1
    stressP1_mc1(i) = 1/4*HA*(1 + d1*(Defmc1_P1(i) - 1))*1/Defmc1_P1(i)*(Defmc1_P1(i).^2 - 1/(Defmc1_P1(i).^2));
    stressP1_trap(i) = 1/4*HA*(1 + d1*(Deftrap_P1(i) - 1))*1/Deftrap_P1(i)*(Deftrap_P1(i).^2 - 1/(Deftrap_P1(i).^2));
    stressP1(i) = stressP1_mc1(i) + stressP1_trap(i);
    % stress at P2
    stressP2_mc1(i) = 1/4*HA*(1 + d1*(Defmc1_P2(i) - 1))*1/Defmc1_P2(i)*(Defmc1_P2(i).^2 - 1/(Defmc1_P2(i).^2));
    stressP2_trap(i) = 1/4*HA*(1 + d1*(Deftrap_P2(i) - 1))*1/Deftrap_P2(i)*(Deftrap_P2(i).^2 - 1/(Deftrap_P2(i).^2));
    stressP2(i) = stressP2_mc1(i) + stressP2_trap(i);
    % stress at P3
    stressP3_mc1(i) = 1/4*HA*(1 + d1*(Defmc1_P3(i) - 1))*1/Defmc1_P3(i)*(Defmc1_P3(i).^2 - 1/(Defmc1_P3(i).^2));
    stressP3_trap(i) = 1/4*HA*(1 + d1*(Deftrap_P3(i) - 1))*1/Deftrap_P3(i)*(Deftrap_P3(i).^2 - 1/(Deftrap_P3(i).^2));
    stressP3(i) = stressP3_mc1(i) + stressP3_trap(i);
end
    % create a matrix (same size than F) which contains the stress values
    % instead of the indexes of points P1, P2 and P3
    CS_F_trap = [stressP1 stressP2 stressP3];
    
    % calculate the average stress for each set of 3 points (each row of
    % the F matrix)
    for i=1:length(CA_F_trap)
        AV_CS_trap(i) = (stressP1(i)+stressP2(i)+stressP3(i))/3;
    end
    AV_CS_trap = AV_CS_trap.';
    CSmax = max(AV_CS_trap);
    CSav = mean(AV_CS_trap);
    
    % calculate the reaction force matrix for each triangle of the F matrix
    for i=1:length(CA_F_trap)
        % origin
        P1(i,:) = V_trap(CA_F_trap(i,1),:);
        P2(i,:) = V_trap(CA_F_trap(i,2),:);
        P3(i,:) = V_trap(CA_F_trap(i,3),:);
        triangle = [P1(i,:);P2(i,:);P3(i,:)];
        OF(i,:) = [(P1(i,1)+P2(i,1)+P3(i,1))/3 (P1(i,2)+P2(i,2)+P3(i,2))/3 (P1(i,3)+P2(i,3)+P3(i,3))/3];
        % direction (unit vector)
        P1P2(i,:) = (V_trap(CA_F_trap(i,2),:) - V_trap(CA_F_trap(i,1),:))/norm(V_trap(CA_F_trap(i,2),:) - V_trap(CA_F_trap(i,1),:));
        P1P3(i,:) = (V_trap(CA_F_trap(i,3),:) - V_trap(CA_F_trap(i,1),:))/norm(V_trap(CA_F_trap(i,3),:) - V_trap(CA_F_trap(i,1),:));
        dir_F(i,:) = cross(P1P2(i,:),P1P3(i,:))/norm(cross(P1P2(i,:),P1P3(i,:)));
        % magnitude
        magn_cart(i) = AV_CS_trap(i)*S_F_trap(i);
        % matrix containing the corresponding vector forces
        F_cart(i,:) = dir_F(i,:)*magn_cart(i);
    end
    magn_cart = magn_cart.';

    % calculate the corresponding reaction torque matrix (N.mm)
    for i=1:length(CA_F_trap)
        MA_cart(i,:) = [CP(1)-OF(i,1) CP(2)-OF(i,2) CP(3)-OF(i,3)];
        T_cart(i,:) = [MA_cart(i,2)*F_cart(i,3)-MA_cart(i,3)*F_cart(i,2)...
            MA_cart(i,3)*F_cart(i,1)-MA_cart(i,1)*F_cart(i,3)...
            MA_cart(i,1)*F_cart(i,2)-MA_cart(i,2)*F_cart(i,1)];
    end
    
    % total torque due to cartilage deformation (N.mm)
    Ttot_cart = sum(T_cart).';
    
    % total initial torque at t=0
    Tinit = TR + Ttot_cart;
    
% TIME ANALYSIS

    % The initial torque calculated before generates a specific angular
    % velocity around the center of rotation. This angular velocity is a
    % function of time. We need to constantly recalculate the total torque
    % at each time step to correct the angular velocity for each
    % corresponding position of the MC1, until equilibrium is reached

    Vi_mc1 = Vor_mc1;    
for i=1:tmax*10
        t(i)=i/10000;    % time (s)
        dt = 1/10000;    % time step delta t (s)
        
        % LOADING
            % New coordinates for point of application
            PointA_idxi = Vi_mc1(idxDPA,:);            
            % Default position of the second point:
            Vorigini = [PointA_idxi(1)-magnitude_x;PointA_idxi(2)-magnitude_y;PointA_idxi(3)-magnitude_z];
            
            % Definition of the resultant force as a vector    
            FR = [PointA_idxi(1)-Vorigini(1) PointA_idxi(2)-Vorigini(2) PointA_idxi(3)-Vorigini(3)];
        
        % Calculation of the shortest distance between the two bones
        % selected in new position
        [idx_cart,D_cart] = knnsearch(V_trap,Vi_mc1);
        [minD_cart,idxD_cart] = min(D_cart);
        SD(i) = minD_cart;
        
        % Definition of the point of contact (CP) which will be our center of
        % rotation at each time step
        CP(i,:) = V_trap(idxD_cart,:);
        
        % Transposition of inertia from center of gravity to point of contact
        CG = polyhedronCentroid(Vi_mc1,F_mc1); % center of gravity of MC1
        dcg_cr = sqrt((CP(i,1)-CG(1))^2+(CP(i,2)-CG(2))^2+(CP(i,3)-CG(3))^2); % distance from CG to CR
        I1cr = I1 - M_mc1*dcg_cr;
        I2cr = I2 - M_mc1*dcg_cr;
        I3cr = I3 - M_mc1*dcg_cr;
    
        % Calculation of the resultant torque (N.mm)
        MA = [PointA_idxi(1) - CP(i,1); PointA_idxi(2) - CP(i,2); PointA_idxi(3) - CP(i,3)]; % Moment arm
        TRi = [MA(2)*FR(3) - MA(3)*FR(2); MA(3)*FR(1) - MA(1)*FR(3); MA(1)*FR(2) - MA(2)*FR(1)];        
        TRi(i,1) = TRi(1);
        TRi(i,2) = TRi(2);
        TRi(i,3) = TRi(3);
                
        % LIGAMENTS    
            % new insertion coordinates
            I_SAOL_i = Vi_mc1(idxDI1,:);
            I_dAOL_i = Vi_mc1(idxDI2,:);
            I_DRL_i = Vi_mc1(idxDI3,:);
            I_POL_i = Vi_mc1(idxDI4,:);
            I_UCL_i = Vi_mc1(idxDI5,:);
            IPsi = [I_SAOL_i;I_dAOL_i;I_DRL_i;I_POL_i;I_UCL_i];                
            % new length of the ligaments
            l_SAOL(i) = norm(I_SAOL_i-O_SAOL_i);
            l_dAOL(i) = norm(I_dAOL_i-O_dAOL_i);
            l_DRL(i) = norm(I_DRL_i-O_DRL_i);
            l_POL(i) = norm(I_POL_i-O_POL_i);
            l_UCL(i) = norm(I_UCL_i-O_UCL_i);
            % stretch of the ligaments
            stretch_SAOL(i)=l_SAOL(i)/linit_SAOL;
            stretch_dAOL(i)=l_dAOL(i)/linit_dAOL;
            stretch_DRL(i)=l_DRL(i)/linit_DRL;
            stretch_POL(i)=l_POL(i)/linit_POL;
            stretch_UCL(i)=l_UCL(i)/linit_UCL;
            % elastic stress in the stretched ligaments
                % SAOL
            if stretch_SAOL(i) > 1
                se33_SAOL(i) = 2*(C1*(1-1/stretch_SAOL(i).^3)+C2*exp(C3*...
                    (stretch_SAOL(i).^2-1).^2)*(stretch_SAOL(i).^2-1));
            else
                se33_SAOL(i) = 0;
            end
                % dAOL
            if stretch_dAOL(i) > 1
                se33_dAOL(i) = 2*(C1*(1-1/stretch_dAOL(i).^3)+C2*exp(C3*...
                    (stretch_dAOL(i).^2-1).^2)*(stretch_dAOL(i).^2-1));
            else
                se33_dAOL(i) = 0;
            end
                % DRL
            if stretch_DRL(i) > 1
                se33_DRL(i) = 2*(C1*(1-1/stretch_DRL(i).^3)+C2*exp(C3*...
                    (stretch_DRL(i).^2-1).^2)*(stretch_DRL(i).^2-1));
            else
                se33_DRL(i) = 0;
            end
                % POL
            if stretch_POL(i) > 1
                se33_POL(i) = 2*(C1*(1-1/stretch_POL(i).^3)+C2*exp(C3*...
                    (stretch_POL(i).^2-1).^2)*(stretch_POL(i).^2-1));
            else
                se33_POL(i) = 0;
            end
                % UCL
            if stretch_UCL(i) > 1
                se33_UCL(i) = 2*(C1*(1-1/stretch_UCL(i).^3)+C2*exp(C3*...
                    (stretch_UCL(i).^2-1).^2)*(stretch_UCL(i).^2-1));
            else
                se33_UCL(i) = 0;
            end
            % viscious stress in the stretched ligaments
                % SAOL
            if stretch_SAOL(i) > 1
                sv33_SAOL(i) = 2*(nu1*(stretch_SAOL(i).^2+2/stretch_SAOL(i)-3)*...
                    (1+1/(2*stretch_SAOL(i).^6))*stretch_SAOL(i)/dt+nu2*...
                    (stretch_SAOL(i).^2-1)^2*stretch_SAOL(i)/dt);
            else
                sv33_SAOL(i) = 0;
            end
                % dAOL
            if stretch_dAOL(i) > 1
                sv33_dAOL(i) = 2*(nu1*(stretch_dAOL(i).^2+2/stretch_dAOL(i)-3)*...
                    (1+1/(2*stretch_dAOL(i).^6))*stretch_dAOL(i)/dt+nu2*...
                    (stretch_dAOL(i).^2-1)^2*stretch_dAOL(i)/dt);
            else
                sv33_dAOL(i) = 0;
            end
                % DRL
            if stretch_DRL(i) > 1
                sv33_DRL(i) = 2*(nu1*(stretch_DRL(i).^2+2/stretch_DRL(i)-3)*...
                    (1+1/(2*stretch_DRL(i).^6))*stretch_DRL(i)/dt+nu2*...
                    (stretch_DRL(i).^2-1)^2*stretch_DRL(i)/dt);
            else
                sv33_DRL(i) = 0;
            end
                % POL
            if stretch_POL(i) > 1
                sv33_POL(i) = 2*(nu1*(stretch_POL(i).^2+2/stretch_POL(i)-3)*...
                    (1+1/(2*stretch_POL(i).^6))*stretch_POL(i)/dt+nu2*...
                    (stretch_POL(i).^2-1)^2*stretch_POL(i)/dt);
            else
                sv33_POL(i) = 0;
            end
                % UCL
            if stretch_UCL(i) > 1
                sv33_UCL(i) = 2*(nu1*(stretch_UCL(i).^2+2/stretch_UCL(i)-3)*...
                    (1+1/(2*stretch_UCL(i).^6))*stretch_UCL(i)/dt+nu2*...
                    (stretch_UCL(i).^2-1)^2*stretch_UCL(i)/dt);
            else
                sv33_UCL(i) = 0;
            end
            % total stress in the ligaments
            stot_SAOL(i) = se33_SAOL(i) + sv33_SAOL(i);
            stot_dAOL(i) = se33_dAOL(i) + sv33_dAOL(i);
            stot_DRL(i) = se33_DRL(i) + sv33_DRL(i);
            stot_POL(i) = se33_POL(i) + sv33_POL(i);
            stot_UCL(i) = se33_UCL(i) + sv33_UCL(i);
            % forces in the ligaments
                % magnitudes
            magn_SAOL(i) = stot_SAOL(i)*CSA_SAOL;
            magn_dAOL(i) = stot_dAOL(i)*CSA_dAOL;
            magn_DRL(i) = stot_DRL(i)*CSA_DRL;
            magn_POL(i) = stot_POL(i)*CSA_POL;
            magn_UCL(i) = stot_UCL(i)*CSA_UCL;
                % directions
            dir_SAOL(i,:) = (O_SAOL_i - I_SAOL_i)/norm(O_SAOL_i - I_SAOL_i);
            dir_dAOL(i,:) = (O_dAOL_i - I_dAOL_i)/norm(O_dAOL_i - I_dAOL_i);
            dir_DRL(i,:) = (O_DRL_i - I_DRL_i)/norm(O_DRL_i - I_DRL_i);
            dir_POL(i,:) = (O_POL_i - I_POL_i)/norm(O_POL_i - I_POL_i);
            dir_UCL(i,:) = (O_UCL_i - I_UCL_i)/norm(O_UCL_i - I_UCL_i);
                % force vectors
            F_SAOL(i,:) = dir_SAOL(i,:)*magn_SAOL(i);
            F_dAOL(i,:) = dir_dAOL(i,:)*magn_dAOL(i);
            F_DRL(i,:) = dir_DRL(i,:)*magn_DRL(i);
            F_POL(i,:) = dir_POL(i,:)*magn_POL(i);
            F_UCL(i,:) = dir_UCL(i,:)*magn_UCL(i);
            % torques created by the ligaments
                % SAOL
            MA_SAOL(i,:) = CP(i,:)-I_SAOL_i;
            T_SAOL(i,:) = [MA_SAOL(i,2)*F_SAOL(i,3)-MA_SAOL(i,3)*F_SAOL(i,2)...
            MA_SAOL(i,3)*F_SAOL(i,1)-MA_SAOL(i,1)*F_SAOL(i,3)...
            MA_SAOL(i,1)*F_SAOL(i,2)-MA_SAOL(i,2)*F_SAOL(i,1)];
                % dAOL
            MA_dAOL(i,:) = CP(i,:)-I_dAOL_i;
            T_dAOL(i,:) = [MA_dAOL(i,2)*F_dAOL(i,3)-MA_dAOL(i,3)*F_dAOL(i,2)...
            MA_dAOL(i,3)*F_dAOL(i,1)-MA_dAOL(i,1)*F_dAOL(i,3)...
            MA_dAOL(i,1)*F_dAOL(i,2)-MA_dAOL(i,2)*F_dAOL(i,1)];
                % DRL
            MA_DRL(i,:) = CP(i,:)-I_DRL_i;
            T_DRL(i,:) = [MA_DRL(i,2)*F_DRL(i,3)-MA_DRL(i,3)*F_DRL(i,2)...
            MA_DRL(i,3)*F_DRL(i,1)-MA_DRL(i,1)*F_DRL(i,3)...
            MA_DRL(i,1)*F_DRL(i,2)-MA_DRL(i,2)*F_DRL(i,1)];
                % POL
            MA_POL(i,:) = CP(i,:)-I_POL_i;
            T_POL(i,:) = [MA_POL(i,2)*F_POL(i,3)-MA_POL(i,3)*F_POL(i,2)...
            MA_POL(i,3)*F_POL(i,1)-MA_POL(i,1)*F_POL(i,3)...
            MA_POL(i,1)*F_POL(i,2)-MA_POL(i,2)*F_POL(i,1)];
                % UCL
            MA_UCL(i,:) = CP(i,:)-I_UCL_i;
            T_UCL(i,:) = [MA_UCL(i,2)*F_UCL(i,3)-MA_UCL(i,3)*F_UCL(i,2)...
            MA_UCL(i,3)*F_UCL(i,1)-MA_UCL(i,1)*F_UCL(i,3)...
            MA_UCL(i,1)*F_UCL(i,2)-MA_UCL(i,2)*F_UCL(i,1)];
    % Total torque from ligaments stretch (N.mm)
    Ttot_lig(i,:) = T_SAOL(i,:) + T_dAOL(i,:) + T_DRL(i,:) + T_POL(i,:) + T_UCL(i,:);
        
        % CARTILAGE
            % Calculation of the matrix Di which contains the Euclidean distances between
            % each vertices of each bone
            Di = pdist2(V_trap,Vi_mc1);
        
            % Save the indexes of the vertices of each bone which are located
            % in the deformation zone
            [Di_trap,Di_mc1] = find(Di<Ttot);
            idxDef = [Di_trap,Di_mc1];
            idxDef(:,any(idxDef==0,1)) = []; % Removes all the zero in the idxD Matrix
            idx_trap = unique(idxDef(:,1));
            idx_mc1 = unique(idxDef(:,2));

            % We can now test to which triangle these vertices correspond in
            % the F matrix and save their indexes
            test_F = ismember(F_trap,idx_trap);
            idx_F_trap = find(test_F(:,1)==1 & test_F(:,2)==1 & test_F(:,3)==1);
                % This idx_F_trap matrix is a Nx1 matrix, where N is the amount
                % of triangles located in the deformation zone, and each value
                % in this matrix corresponds to the index of a line in the F
                % matrix that defines one of these triangles
            for j=1:length(idx_F_trap)
                CA_F_trap(i,:) = F_trap(idx_F_trap(j),:);
            end
    
                % Each triangle is composed of 3 points called a,b and c
                Xa1 = V_trap(CA_F_trap(:,1),1);
                Ya1 = V_trap(CA_F_trap(:,1),2);
                Za1 = V_trap(CA_F_trap(:,1),3);

                Xb1 = V_trap(CA_F_trap(:,2),1);
                Yb1 = V_trap(CA_F_trap(:,2),2);
                Zb1 = V_trap(CA_F_trap(:,2),3);

                Xc1 = V_trap(CA_F_trap(:,3),1);
                Yc1 = V_trap(CA_F_trap(:,3),2);
                Zc1 = V_trap(CA_F_trap(:,3),3);

            % Now we need to calculate the coordinates of the 2 vectors ab and ac
                for j=1:length(CA_F_trap)
                    Xab1(j) = Xb1(j)-Xa1(j);
                    Yab1(j) = Yb1(j)-Ya1(j);
                    Zab1(j) = Zb1(j)-Za1(j);
                    Xac1(j) = Xc1(j)-Xa1(j);
                    Yac1(j) = Yc1(j)-Ya1(j);
                    Zac1(j) = Zc1(j)-Za1(j);
                end

            % Now we need to calculate the surface area of each triangle and sum
            % all of them to obtain to total contact area
                for j = 1:length(CA_F_trap)
                    S_F_trap(j) = (1/2)*((Yab1(j)*Zac1(j)-Zab1(j)*Yac1(j))^2+...
                        (Zab1(j)*Xac1(j)-Xab1(j)*Zac1(j))^2+...
                        (Xab1(j)*Yac1(j)-Yab1(j)*Xac1(j))^2)^(1/2);
                end
                S_F_trap = S_F_trap.';
                PCA_trap(i) = sum(S_F_trap(:));

            % We need now to calculate the corresponding stress value for each
            % vertex which appears in the F matrix, so that we can calculate
            % the average stress value for the whole triangle and convert it
            % into force

            for j=1:length(CA_F_trap)
                % for each point P1, P2 and P3 of the F matrix of the trapezium, we
                % need to know the index of the closest point in the V matrix of the
                % MC1
                Di_trap_mc1_P1 = pdist2(V_trap(CA_F_trap(j,1),:),Vi_mc1);
                [minDiP1(j),idxP1(j)] = min(Di_trap_mc1_P1);
                Di_trap_mc1_P2 = pdist2(V_trap(CA_F_trap(j,2),:),Vi_mc1);
                [minDiP2(j),idxP2(j)] = min(Di_trap_mc1_P2);
                Di_trap_mc1_P3 = pdist2(V_trap(CA_F_trap(j,3),:),Vi_mc1);
                [minDiP3(j),idxP3(j)] = min(Di_trap_mc1_P3);
                % distance between each vertices for P1, P2 and P3 of the F matrix
                DiP1(j) = ((V_trap(CA_F_trap(j,1),1)-Vi_mc1(idxP1(j),1)).^2+...
                    (V_trap(CA_F_trap(j,1),2)-Vi_mc1(idxP1(j),2)).^2+...
                    (V_trap(CA_F_trap(j,1),3)-Vi_mc1(idxP1(j),3)).^2).^(1/2);
                DiP2(j) = ((V_trap(CA_F_trap(j,2),1)-Vi_mc1(idxP2(j),1)).^2+...
                    (V_trap(CA_F_trap(j,2),2)-Vi_mc1(idxP2(j),2)).^2+...
                    (V_trap(CA_F_trap(j,2),3)-Vi_mc1(idxP2(j),3)).^2).^(1/2);
                DiP3(j) = ((V_trap(CA_F_trap(j,3),1)-Vi_mc1(idxP3(j),1)).^2+...
                    (V_trap(CA_F_trap(j,3),2)-Vi_mc1(idxP3(j),2)).^2+...
                    (V_trap(CA_F_trap(j,3),3)-Vi_mc1(idxP3(j),3)).^2).^(1/2);
                % deformation in MC1 cartilage layer for P1, P2 and P3 of the F matrix
                Defmc1_P1(j) = (Ttot-DiP1(j))/2/Tmc1+1;
                Defmc1_P2(j) = (Ttot-DiP2(j))/2/Tmc1+1;
                Defmc1_P3(j) = (Ttot-DiP3(j))/2/Tmc1+1;
                % deformation in MC1 cartilage layer for P1, P2 and P3 of the F matrix
                Deftrap_P1(j) = (Ttot-DiP1(j))/2/Ttrap+1;
                Deftrap_P2(j) = (Ttot-DiP2(j))/2/Ttrap+1;
                Deftrap_P3(j) = (Ttot-DiP3(j))/2/Ttrap+1;
                % stress at P1
                stressP1_mc1(j) = 1/4*HA*(1 + d1*(Defmc1_P1(j) - 1))*1/Defmc1_P1(j)*(Defmc1_P1(j).^2 - 1/(Defmc1_P1(j).^2));
                stressP1_trap(j) = 1/4*HA*(1 + d1*(Deftrap_P1(j) - 1))*1/Deftrap_P1(j)*(Deftrap_P1(j).^2 - 1/(Deftrap_P1(j).^2));
                stressP1(j) = stressP1_mc1(j) + stressP1_trap(j);
                % stress at P2
                stressP2_mc1(j) = 1/4*HA*(1 + d1*(Defmc1_P2(j) - 1))*1/Defmc1_P2(j)*(Defmc1_P2(j).^2 - 1/(Defmc1_P2(j).^2));
                stressP2_trap(j) = 1/4*HA*(1 + d1*(Deftrap_P2(j) - 1))*1/Deftrap_P2(j)*(Deftrap_P2(j).^2 - 1/(Deftrap_P2(j).^2));
                stressP2(j) = stressP2_mc1(j) + stressP2_trap(j);
                % stress at P3
                stressP3_mc1(j) = 1/4*HA*(1 + d1*(Defmc1_P3(j) - 1))*1/Defmc1_P3(j)*(Defmc1_P3(j).^2 - 1/(Defmc1_P3(j).^2));
                stressP3_trap(j) = 1/4*HA*(1 + d1*(Deftrap_P3(j) - 1))*1/Deftrap_P3(j)*(Deftrap_P3(j).^2 - 1/(Deftrap_P3(j).^2));
                stressP3(j) = stressP3_mc1(j) + stressP3_trap(j);
            end
            % create a matrix (same size than F) which contains the stress values
            % instead of the indexes of points P1, P2 and P3
            CS_F_trap = [stressP1 stressP2 stressP3];
    
            % calculate the average stress for each set of 3 points (each row of
            % the F matrix)
            for j=1:length(CA_F_trap)
                AV_CS_trap(j) = (stressP1(j)+stressP2(j)+stressP3(j))/3;
            end
            AV_CS_trap = AV_CS_trap.';
            CSmax(i) = max(AV_CS_trap);
            CSav(i) = mean(AV_CS_trap);    
        
            % calculate the reaction force matrix for each triangle of the new F matrix
                for j=1:length(CA_F_trap)
                    % origin
                    Pi1(j,:) = V_trap(CA_F_trap(j,1),:);
                    Pi2(j,:) = V_trap(CA_F_trap(j,2),:);
                    Pi3(j,:) = V_trap(CA_F_trap(j,3),:);
                    triangle = [Pi1(j,:);Pi2(j,:);Pi3(j,:)];
                    OF(j,:) = [(Pi1(j,1)+Pi2(j,1)+Pi3(j,1))/3 (Pi1(j,2)+Pi2(j,2)+Pi3(j,2))/3 (Pi1(j,3)+Pi2(j,3)+Pi3(j,3))/3];
                    % direction (unit vector)
                    Pi1Pi2(j,:) = (V_trap(CA_F_trap(j,2),:) - V_trap(CA_F_trap(j,1),:))/norm(V_trap(CA_F_trap(j,2),:) - V_trap(CA_F_trap(j,1),:));
                    Pi1Pi3(j,:) = (V_trap(CA_F_trap(j,3),:) - V_trap(CA_F_trap(j,1),:))/norm(V_trap(CA_F_trap(j,3),:) - V_trap(CA_F_trap(j,1),:));
                    dir_Fi(j,:) = cross(Pi1Pi2(j,:),Pi1Pi3(j,:))/norm(cross(Pi1Pi2(j,:),Pi1Pi3(j,:)));
                    % magnitude
                    magn_carti(j) = AV_CS_trap(j)*S_F_trap(j);
                    % matrix containing the corresponding vector forces
                    F_cart(j,:) = dir_F(j,:)*magn_cart(j);

                % calculate the corresponding reaction torque matrix (N.mm)
                    MA_cart(j,:) = [CP(i,1)-OF(j,1) CP(i,2)-OF(j,2) CP(i,3)-OF(j,3)];
                    T_cart(j,:) = [MA_cart(j,2)*F_cart(j,3)-MA_cart(j,3)*F_cart(j,2)...
                    MA_cart(j,3)*F_cart(j,1)-MA_cart(j,1)*F_cart(j,3)...
                    MA_cart(j,1)*F_cart(j,2)-MA_cart(j,2)*F_cart(j,1)];
                end
    % Total torque due to cartilage deformation (N.mm)
        Ttot_carti(i,:) = [sum(T_cart(:,1)) sum(T_cart(:,2)) sum(T_cart(:,3))];

    % TOTAL TORQUE AT TIME T (N.mm)
    Tt(i,1) = TRi(i,1) + Ttot_lig(i,1) + Ttot_carti(i,1);
    Tt(i,2) = TRi(i,2) + Ttot_lig(i,2) + Ttot_carti(i,2);
    Tt(i,3) = TRi(i,3) + Ttot_lig(i,3) + Ttot_carti(i,3);
    
    % ANGULAR VELOCITY (rad)
    Rot_x(i) = Tt(i,1)*10^-3/I1cr*t(i);
    Rot_y(i) = Tt(i,2)*10^-3/I2cr*t(i);
    Rot_z(i) = Tt(i,3)*10^-3/I3cr*t(i);
    % (degrees)
    Rot_deg_x(i) = Rot_x(i)*180/pi;
    Rot_deg_y(i) = Rot_y(i)*180/pi;
    Rot_deg_z(i) = Rot_z(i)*180/pi;
    
    % ROTATION MATRIX FOR NEW POSITION OF MC1 DUE TO TOTAL TORQUE
    Roti_mc1 = [cos(Rot_z(i))*cos(Rot_y(i)),...
        cos(Rot_z(i))*sin(Rot_y(i))*sin(Rot_x(i))-sin(Rot_z(i))*cos(Rot_x(i)),...
        cos(Rot_z(i))*sin(Rot_y(i))*cos(Rot_x(i))+sin(Rot_z(i))*sin(Rot_x(i));...
        sin(Rot_z(i))*cos(Rot_y(i)),...
        sin(Rot_z(i))*sin(Rot_y(i))*sin(Rot_x(i))+cos(Rot_z(i))*cos(Rot_x(i)),...
        sin(Rot_z(i))*sin(Rot_y(i))*cos(Rot_x(i))-cos(Rot_z(i))*sin(Rot_x(i));...
        -sin(Rot_y(i)), cos(Rot_y(i))*sin(Rot_x(i)), cos(Rot_y(i))*cos(Rot_x(i))];
    
    % ALIGNMENT OF THE CENTER OF ROTATION AND THE ORIGIN
    Trans_CP = -CP(i,:).';
    Vi_mc1 = Vi_mc1.' + repmat(Trans_CP,1,size(Vi_mc1,1));
    Vi_mc1 = Vi_mc1.';    
    
    % ROTATION OF MC1 DUE TO TOTAL TORQUE
    Vi_mc1 = Roti_mc1*Vi_mc1.';
    Vi_mc1 = Vi_mc1.';
    
    % NEW POSITION OF MC1 DUE TO TOTAL TORQUE
    Vi_mc1 = Vi_mc1.' - repmat(Trans_CP,1,size(Vi_mc1,1));
    Vi_mc1 = Vi_mc1.';
   
    % CONDITION TO STOP
        % if the total torque gets close to 0 (equilibrium reached)
    if abs(Tt(i,1)) < 0.1 && abs(Tt(i,2)) < 0.1 && abs(Tt(i,3)) < 0.1
        break
    end
        % if the total torque gets too high (system getting unstable)
    if abs(Tt(i,1)) > 1000 || abs(Tt(i,2)) > 1000 || abs(Tt(i,3)) > 1000
        break
    end
end

% FINAL COORDINATES FOR LIGAMENTS INSERTION AND CENTER OF ROTATION       
    % FINAL INSERTIONS COORDINATES FOR LIGAMENTS
    I_SAOL_i = Vi_mc1(idxDI1,:);
    I_dAOL_i = Vi_mc1(idxDI2,:);
    I_DRL_i = Vi_mc1(idxDI3,:);
    I_POL_i = Vi_mc1(idxDI4,:);
    I_UCL_i = Vi_mc1(idxDI5,:);
    IPsi = [I_SAOL_i;I_dAOL_i;I_DRL_i;I_POL_i;I_UCL_i];
    
% FINAL COORDINATES FOR FORCE APPLICATION AND DIRECTION
    PointA_idxi = Vi_mc1(idxDPA,:);            
    Vorigini = [PointA_idxi(1)-magnitude_x;PointA_idxi(2)-magnitude_y;PointA_idxi(3)-magnitude_z];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FINAL RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plots the system in its final position    
    f10 = figure;
    [obj, li, ax] = GUI_PlotShells(f10, {F_trap;F_mc1;F_mc1},{V_trap;Vi_mc1;Vor_mc1},...
        {ones(size(V_trap,1),1);ones(size(Vi_mc1,1),1);ones(size(Vor_mc1,1),1)});
    box off
    view([45,45,45]);
    hold on
    arrow([Vorigini(1);Vorigini(2);Vorigini(3)],[PointA_idxi(1);PointA_idxi(2);PointA_idxi(3)]...
        ,1,70,30,5,'EdgeColor','r','FaceColor','r');
    plot3(PointA_idxi(1),PointA_idxi(2),PointA_idxi(3),'Marker','o','MarkerSize',10,'Color',...
        'k','MarkerFaceColor','k');
    hold on
    for i=1:5;
        arrow([OPs(i,1);OPs(i,2);OPs(i,3)],[IPsi(i,1);IPsi(i,2);IPsi(i,3)],...
        10,70,30,5,'EdgeColor',[1 0.5 0.2],'FaceColor',[1 0.5 0.2])
        plot3(OPs(i,1),OPs(i,2),OPs(i,3),'Marker','o','MarkerSize',10,'Color',...
        'r','MarkerFaceColor','r');
        plot3(IPsi(i,1),IPsi(i,2),IPsi(i,3),'Marker','o','MarkerSize',10,'Color',...
        'b','MarkerFaceColor','b');
    end
    % Plots the trapezium CS
    hold on
    arrow([0;0;0],[15;0;0],10,70,30,5,'EdgeColor','r','FaceColor','r');
    arrow([0;0;0],[0;15;0],10,70,30,5,'EdgeColor','g','FaceColor','g');
    arrow([0;0;0],[0;0;15],10,70,30,5,'EdgeColor','b','FaceColor','b');
    set(gcf,'numbertitle','off','name','System in its FINAL configuration (purple: MC1 initial position; green: MC1 final position)');

% Plots the trapezium with a color map representing to stress distribution
% in the initial stage of the system
colormap_initial = figure; 
    [obj, li, ax] = GUI_PlotShells(colormap_initial, {F_trap}, {V_trap},...
            {ones(size(V_trap,1),1)},[0,0,1]);
box off
view(90,0);
for i=1:length(V_trap)
    colrstl_init(i,:)=[0,0,1];
end

hold on
for j=1:length(CA_F_trap)
    % for each point P1, P2 and P3 of the F matrix of the trapezium, we
    % need to know the index of the closest point in the V matrix of the
    % MC1
    Di_trap_mc1_P1_init = pdist2(V_trap(CA_F_trap(j,1),:),Vor_mc1);
    [minDiP1(j),idxP1_init(j)] = min(Di_trap_mc1_P1_init);
    Di_trap_mc1_P2_init = pdist2(V_trap(CA_F_trap(j,2),:),Vor_mc1);
    [minDiP2(j),idxP2_init(j)] = min(Di_trap_mc1_P2_init);
    Di_trap_mc1_P3_init = pdist2(V_trap(CA_F_trap(j,3),:),Vor_mc1);
    [minDiP3(j),idxP3_init(j)] = min(Di_trap_mc1_P3_init);
    % distance between each vertices for P1, P2 and P3 of the F matrix
    DiP1_init(j) = ((V_trap(CA_F_trap(j,1),1)-Vor_mc1(idxP1_init(j),1)).^2+...
        (V_trap(CA_F_trap(j,1),2)-Vor_mc1(idxP1_init(j),2)).^2+...
        (V_trap(CA_F_trap(j,1),3)-Vor_mc1(idxP1_init(j),3)).^2).^(1/2);
    DiP2_init(j) = ((V_trap(CA_F_trap(j,2),1)-Vor_mc1(idxP2_init(j),1)).^2+...
        (V_trap(CA_F_trap(j,2),2)-Vor_mc1(idxP2_init(j),2)).^2+...
        (V_trap(CA_F_trap(j,2),3)-Vor_mc1(idxP2_init(j),3)).^2).^(1/2);
    DiP3_init(j) = ((V_trap(CA_F_trap(j,3),1)-Vor_mc1(idxP3_init(j),1)).^2+...
        (V_trap(CA_F_trap(j,3),2)-Vor_mc1(idxP3_init(j),2)).^2+...
        (V_trap(CA_F_trap(j,3),3)-Vor_mc1(idxP3_init(j),3)).^2).^(1/2);
    % deformation in MC1 cartilage layer for P1, P2 and P3 of the F matrix
    Defmc1_P1_init(j) = (Ttot-DiP1_init(j))/2/Tmc1+1;
    Defmc1_P2_init(j) = (Ttot-DiP2_init(j))/2/Tmc1+1;
    Defmc1_P3_init(j) = (Ttot-DiP3_init(j))/2/Tmc1+1;
    % deformation in MC1 cartilage layer for P1, P2 and P3 of the F matrix
    Deftrap_P1_init(j) = (Ttot-DiP1_init(j))/2/Ttrap+1;
    Deftrap_P2_init(j) = (Ttot-DiP2_init(j))/2/Ttrap+1;
    Deftrap_P3_init(j) = (Ttot-DiP3_init(j))/2/Ttrap+1;
    % stress at P1
    stressP1_mc1_init(j) = 1/4*HA*(1 + d1*(Defmc1_P1_init(j) - 1))*1/Defmc1_P1_init(j)*(Defmc1_P1_init(j).^2 - 1/(Defmc1_P1_init(j).^2));
    stressP1_trap_init(j) = 1/4*HA*(1 + d1*(Deftrap_P1_init(j) - 1))*1/Deftrap_P1_init(j)*(Deftrap_P1_init(j).^2 - 1/(Deftrap_P1_init(j).^2));
    stressP1_init(j) = stressP1_mc1_init(j) + stressP1_trap_init(j);
    % stress at P2
    stressP2_mc1_init(j) = 1/4*HA*(1 + d1*(Defmc1_P2_init(j) - 1))*1/Defmc1_P2_init(j)*(Defmc1_P2_init(j).^2 - 1/(Defmc1_P2_init(j).^2));
    stressP2_trap_init(j) = 1/4*HA*(1 + d1*(Deftrap_P2_init(j) - 1))*1/Deftrap_P2_init(j)*(Deftrap_P2_init(j).^2 - 1/(Deftrap_P2_init(j).^2));
    stressP2_init(j) = stressP2_mc1_init(j) + stressP2_trap_init(j);
    % stress at P3
    stressP3_mc1_init(j) = 1/4*HA*(1 + d1*(Defmc1_P3_init(j) - 1))*1/Defmc1_P3_init(j)*(Defmc1_P3_init(j).^2 - 1/(Defmc1_P3_init(j).^2));
    stressP3_trap_init(j) = 1/4*HA*(1 + d1*(Deftrap_P3_init(j) - 1))*1/Deftrap_P3_init(j)*(Deftrap_P3_init(j).^2 - 1/(Deftrap_P3_init(j).^2));
    stressP3_init(j) = stressP3_mc1_init(j) + stressP3_trap_init(j);
       
        colrstl_init(CA_F_trap(j,1),:)= ...
        [-((stressP1_init(j)-HA)^2)/((HA)^2)+1,...
        -((stressP1_init(j)-HA/2)^2)/((-HA/2)^2)+1,...
        -(stressP1_init(j)^2)/((HA)^2)+1];
    if stressP1_init(j) >= HA
       colrstl_init(CA_F_trap(j,1),:) = [1,0,0];
    end
        colrstl_init(CA_F_trap(j,2),:)= ...
        [-((stressP2_init(j)-HA)^2)/((HA)^2)+1,...
        -((stressP2_init(j)-HA/2)^2)/((-HA/2)^2)+1,...
        -(stressP2_init(j)^2)/((HA)^2)+1];
    if stressP2_init(j) >= HA
       colrstl_init(CA_F_trap(j,2),:) = [1,0,0];
    end
            colrstl_init(CA_F_trap(j,3),:)= ...
        [-((stressP3_init(j)-HA)^2)/((HA)^2)+1,...
        -((stressP3_init(j)-HA/2)^2)/((-HA/2)^2)+1,...
        -(stressP3_init(j)^2)/((HA)^2)+1];
    if stressP3_init(j) >= HA
       colrstl_init(CA_F_trap(j,3),:) = [1,0,0];
    end
end
hold on
patch('Faces',F_trap,'Vertices',V_trap, ...
    'FaceColor','interp', ...
    'FaceVertexCData',colrstl_init, ...
    'EdgeColor', 'interp', ...
    'EdgeAlpha', 0, ...
    'CDataMapping', 'scaled',...
    'AmbientStrength', 0.4, ...
    'DiffuseStrength', 0.8, ...
    'SpecularStrength', 0.2, ...
    'SpecularColorReflectance', 0.5, ...
    'FaceLighting', 'gouraud')
hold on
h=colorbar('YTickLabel',{'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8 = HA or more'});
set(h,'ylim',[0.1 0.9]);
title = get(h,'Title');
titleString = 'MPa';
set(title ,'String',titleString,'FontWeight','bold');
set(gcf,'numbertitle','off','name','INITIAL stress distribution on the trapezium');
    
% Plots the trapezium with a color map representing to stress distribution
% in the final stage of the system
colormap_final = figure; 
    [obj, li, ax] = GUI_PlotShells(colormap_final, {F_trap}, {V_trap},...
            {ones(size(V_trap,1),1)},[0,0,1]);
box off
view(90,0);
for i=1:length(V_trap)
    colrstl(i,:)=[0,0,1];
end

hold on
for j=1:length(CA_F_trap)
    % for each point P1, P2 and P3 of the F matrix of the trapezium, we
    % need to know the index of the closest point in the V matrix of the
    % MC1
    Di_trap_mc1_P1 = pdist2(V_trap(CA_F_trap(j,1),:),Vi_mc1);
    [minDiP1(j),idxP1(j)] = min(Di_trap_mc1_P1);
    Di_trap_mc1_P2 = pdist2(V_trap(CA_F_trap(j,2),:),Vi_mc1);
    [minDiP2(j),idxP2(j)] = min(Di_trap_mc1_P2);
    Di_trap_mc1_P3 = pdist2(V_trap(CA_F_trap(j,3),:),Vi_mc1);
    [minDiP3(j),idxP3(j)] = min(Di_trap_mc1_P3);
    % distance between each vertices for P1, P2 and P3 of the F matrix
    DiP1(j) = ((V_trap(CA_F_trap(j,1),1)-Vi_mc1(idxP1(j),1)).^2+...
        (V_trap(CA_F_trap(j,1),2)-Vi_mc1(idxP1(j),2)).^2+...
        (V_trap(CA_F_trap(j,1),3)-Vi_mc1(idxP1(j),3)).^2).^(1/2);
    DiP2(j) = ((V_trap(CA_F_trap(j,2),1)-Vi_mc1(idxP2(j),1)).^2+...
        (V_trap(CA_F_trap(j,2),2)-Vi_mc1(idxP2(j),2)).^2+...
        (V_trap(CA_F_trap(j,2),3)-Vi_mc1(idxP2(j),3)).^2).^(1/2);
    DiP3(j) = ((V_trap(CA_F_trap(j,3),1)-Vi_mc1(idxP3(j),1)).^2+...
        (V_trap(CA_F_trap(j,3),2)-Vi_mc1(idxP3(j),2)).^2+...
        (V_trap(CA_F_trap(j,3),3)-Vi_mc1(idxP3(j),3)).^2).^(1/2);
    % deformation in MC1 cartilage layer for P1, P2 and P3 of the F matrix
    Defmc1_P1(j) = (Ttot-DiP1(j))/2/Tmc1+1;
    Defmc1_P2(j) = (Ttot-DiP2(j))/2/Tmc1+1;
    Defmc1_P3(j) = (Ttot-DiP3(j))/2/Tmc1+1;
    % deformation in MC1 cartilage layer for P1, P2 and P3 of the F matrix
    Deftrap_P1(j) = (Ttot-DiP1(j))/2/Ttrap+1;
    Deftrap_P2(j) = (Ttot-DiP2(j))/2/Ttrap+1;
    Deftrap_P3(j) = (Ttot-DiP3(j))/2/Ttrap+1;
    % stress at P1
    stressP1_mc1(j) = 1/4*HA*(1 + d1*(Defmc1_P1(j) - 1))*1/Defmc1_P1(j)*(Defmc1_P1(j).^2 - 1/(Defmc1_P1(j).^2));
    stressP1_trap(j) = 1/4*HA*(1 + d1*(Deftrap_P1(j) - 1))*1/Deftrap_P1(j)*(Deftrap_P1(j).^2 - 1/(Deftrap_P1(j).^2));
    stressP1(j) = stressP1_mc1(j) + stressP1_trap(j);
    % stress at P2
    stressP2_mc1(j) = 1/4*HA*(1 + d1*(Defmc1_P2(j) - 1))*1/Defmc1_P2(j)*(Defmc1_P2(j).^2 - 1/(Defmc1_P2(j).^2));
    stressP2_trap(j) = 1/4*HA*(1 + d1*(Deftrap_P2(j) - 1))*1/Deftrap_P2(j)*(Deftrap_P2(j).^2 - 1/(Deftrap_P2(j).^2));
    stressP2(j) = stressP2_mc1(j) + stressP2_trap(j);
    % stress at P3
    stressP3_mc1(j) = 1/4*HA*(1 + d1*(Defmc1_P3(j) - 1))*1/Defmc1_P3(j)*(Defmc1_P3(j).^2 - 1/(Defmc1_P3(j).^2));
    stressP3_trap(j) = 1/4*HA*(1 + d1*(Deftrap_P3(j) - 1))*1/Deftrap_P3(j)*(Deftrap_P3(j).^2 - 1/(Deftrap_P3(j).^2));
    stressP3(j) = stressP3_mc1(j) + stressP3_trap(j);
       
        colrstl(CA_F_trap(j,1),:)= ...
        [-((stressP1(j)-HA)^2)/((HA)^2)+1,...
        -((stressP1(j)-HA/2)^2)/((-HA/2)^2)+1,...
        -(stressP1(j)^2)/((HA)^2)+1];
    if stressP1(j) >= HA
       colrstl(CA_F_trap(j,1),:) = [1,0,0];
    end
        colrstl(CA_F_trap(j,2),:)= ...
        [-((stressP2(j)-HA)^2)/((HA)^2)+1,...
        -((stressP2(j)-HA/2)^2)/((-HA/2)^2)+1,...
        -(stressP2(j)^2)/((HA)^2)+1];
    if stressP2(j) >= HA
       colrstl(CA_F_trap(j,2),:) = [1,0,0];
    end
            colrstl(CA_F_trap(j,3),:)= ...
        [-((stressP3(j)-HA)^2)/((HA)^2)+1,...
        -((stressP3(j)-HA/2)^2)/((-HA/2)^2)+1,...
        -(stressP3(j)^2)/((HA)^2)+1];
    if stressP3(j) >= HA
       colrstl(CA_F_trap(j,3),:) = [1,0,0];
    end
end
hold on
patch('Faces',F_trap,'Vertices',V_trap, ...
    'FaceColor','interp', ...
    'FaceVertexCData',colrstl, ...
    'EdgeColor', 'interp', ...
    'EdgeAlpha', 0, ...
    'CDataMapping', 'scaled',...
    'AmbientStrength', 0.4, ...
    'DiffuseStrength', 0.8, ...
    'SpecularStrength', 0.2, ...
    'SpecularColorReflectance', 0.5, ...
    'FaceLighting', 'gouraud')
hold on
h=colorbar('YTickLabel',{'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8 = HA or more'});
set(h,'ylim',[0.1 0.9]);
title = get(h,'Title');
titleString = 'MPa';
set(title ,'String',titleString,'FontWeight','bold');
set(gcf,'numbertitle','off','name','FINAL stress distribution on the trapezium');

% Plots the evolution of the tension in the ligaments through time
f7=figure;
plot(t,magn_SAOL,'r',t,magn_dAOL,'g',t,magn_DRL,'b',t,magn_POL,'m',t,magn_UCL,'c');
xlabel('time (ms)');
ylabel('Tension in the ligaments (N)')
legend('SAOL','dAOL','DRL','POL','UCL','Location','NorthEastOutside')
set(gcf,'numbertitle','off','name','Evolution of the tension in the ligaments');

% Plots the evolution of the maximal and average contact stress through time
f8=figure;
plot(t,CSmax,'r',t,CSav,'g');
xlabel('time (ms)');
ylabel('Contact stress (MPa)')
legend('CSmax','CSav','Location','NorthEastOutside')
set(gcf,'numbertitle','off','name','Evolution of the maximal and average contact stress');

% Plots the evolution of the total torque through time
f9=figure;
plot(t,Tt(:,1),'r',t,Tt(:,2),'g',t,Tt(:,3),'b',t,t*0,'k');
xlabel('time (ms)');
ylabel('Torque (N.mm)')
legend('Ttot-x','Ttot-y','Ttot-z','Location','NorthEastOutside')
set(gcf,'numbertitle','off','name','Evolution of the total torque');

for i=1:length(Tt)
    Tin(i,1)=TR(1);
    Tin(i,2)=TR(2);
    Tin(i,3)=TR(3);
end
% Plots the evolution of the torque on the x-axis through time
f10=figure;
plot(t,Tin(:,1),'r',t,Ttot_lig(:,1),'g',t,Ttot_carti(:,1),'b',t,t*0,'k');
xlabel('time (ms)');
ylabel('Tx (N.mm)')
legend('Tin-x','Tlig-x','Tcart-x','Location','NorthEastOutside')
set(gcf,'numbertitle','off','name','Evolution of the torque on X');

% Plots the evolution of the torque on the y-axis through time
f11=figure;
plot(t,Tin(:,2),'r',t,Ttot_lig(:,2),'g',t,Ttot_carti(:,2),'b',t,t*0,'k');
xlabel('time (ms)');
ylabel('Ty (N.mm)')
legend('Tin-y','Tlig-y','Tcart-y','Location','NorthEastOutside')
set(gcf,'numbertitle','off','name','Evolution of the torque on Y');

% Plots the evolution of the torque on the y-axis through time
f12=figure;
plot(t,Tin(:,3),'r',t,Ttot_lig(:,3),'g',t,Ttot_carti(:,3),'b',t,t*0,'k');
xlabel('time (ms)');
ylabel('Tz (N.mm)')
legend('Tin-z','Tlig-z','Tcart-z','Location','NorthEastOutside')
set(gcf,'numbertitle','off','name','Evolution of the torque on Z');

% Plots the evolution of the projected contact area through time
f13=figure;
plot(t,PCA_trap,'b');
xlabel('time (ms)');
ylabel('Projected contact area (mm2)')
set(gcf,'numbertitle','off','name','Evolution of the projected contact area');

% Plots the evolution of the rotation angles
f14=figure;
plot(t,Rot_deg_x,'r',t,Rot_deg_y,'g',t,Rot_deg_z,'b',t,t*0,'k');
xlabel('time (ms)');
ylabel('Rotation angle (degrees)')
legend('Rot_x','Rot_y','Rot_z','Location','NorthEastOutside')
set(gcf,'numbertitle','off','name','Evolution of the rotation angles');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%