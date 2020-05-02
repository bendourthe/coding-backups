clear all; close all;
% Function: Define the trapezium coordinate system (CS), perform a 4-points
% registration based on an ICP algorithm to calculate the transformation
% matrices from a position 2 to a position 1, calculates the helical axis
% of the joint and show it on a figure.
% Also calculates the shortest distance between the two bones of the joint
% for both frames and show the proximity patterns on the articular surfaces
% of each bone in both frames.
%
% Dependencies:        
%               txt2mat.m
%               STL_ReadFile.m
%               PlacePoints3.m
%               RegisterBones.m
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
%                   TRI_Areas.m
%                   KIN_DirectedSpringRotation.m
%               GUI_PlotShells.m
%               IntersectLineAndPlane.m 
%                   DistanceFromVertexToPlane.m 
%                       TRI_Normals.m	
%               VectorNorms.m
%               screw.m
%               FindPivotPoint
%               DistanceFromVertexToLine
%               arrow.m
% 
% Input: STl files of the trapezium and first metacarpal(MC1).
%
% Output: txt-files with transformations and coordinate systems. Figures 
% showing the efficiency of the registration and the helical axis. Figures
% showing the proximity patterns.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Control panel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Folder with files
dir = 'H:\Data\Study_Brown_Rig\Segmentation\Scan27\STL\'; 
data_available = 0; % loads data if available (=1), else, generate data
    filename = 'RotTransData_f1f2.txt';
    file = [dir filename];
trap_cs_available = 0;          % loads data if available (=1), else, generate data
    filenameCS_trap = 'TrapAxes_f1f2.txt';
    fileCS_trap = [dir filenameCS_trap];
helical_axis_available = 0; % loads data if available (=1), else, generate data
    filenameHelicalAX = 'HelicalAx_f1f2.txt';
    fileHelicalAX = [dir filenameHelicalAX];
rot_angle_available = 0; % loads data if available (=1), else, generate data
    filenameRotAngle = 'RotAngle_f1f2.txt';
    fileRotAngle = [dir filenameRotAngle];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Note: the code works so that P2 is registered on P1, so if you want
    % to study the motion from Frame 1 to Frame 2, you should put Frame 1
    % as P2 and Frame 2 as P1 (and vice-versa)
% First position (referenced as P1)
MC1_P1 = 'SCAN27_ADD_mc1.stl';
Trap_P1 = 'SCAN27_ADD_trap.stl';

% Second position (referenced as P2)
MC1_P2 = 'SCAN27_ABD_mc1.stl';
Trap_P2 = 'SCAN27_ABD_trap.stl';

% Setting the maximal distance that the code will consider between the 2
% bones and the distance that we consider critical (underneath which we can
% expect high amount of stress due to high compression)
dist_max = 3;       % in mm
dist_prox =   1.5;    % in mm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COORDINATE SYSTEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    
    X1 = A(1,:);    % mid-point of the central ridge of the trapezial saddle
    X2 = A(2,:);    % center of the junction of the trapezium, scaphoid and trapezoid
    Y = A(3,:);     % center of the distal articular surface, facing the center of the base of MC1
else
    [F_Trap_P1, V_Trap_P1] =  STL_ReadFile([dir Trap_P1],true);
    [ad_curve] = PlacePoints3({F_Trap_P1}, {V_Trap_P1}, {F_Trap_P1}, {V_Trap_P1}, 'Select 3 landmarks (X1, X2 and Y, in that order!). End by pressing enter and exiting the figure');
    close all;
    % Retrieves coordinates of the points:
    X1 = ad_curve.points{end,1}(1,:);
    X2 = ad_curve.points{end,1}(2,:);
    Y = ad_curve.points{end,1}(3,:);
    
    % Save in a txt file
    fid = fopen(fileCS_trap,'w+');
    fprintf(fid,'%5.8f\t%5.8f\t%5.8f\r\n',[X1,X2,Y]);
    fclose(fid);
end

O_trap = (X1+X2)/2;
Xaxis_trap = (X2-X1)/norm(X2-X1);
Yaxis_trap = (O_trap- Y)/norm(O_trap-Y);
Zaxis_trap = cross(Xaxis_trap,Yaxis_trap)/norm(cross(Xaxis_trap,Yaxis_trap));
Yaxis_trap = cross(Zaxis_trap,Xaxis_trap)/norm(cross(Zaxis_trap,Xaxis_trap));
cs_trap = [Xaxis_trap; Yaxis_trap; Zaxis_trap];
% Converges to trapezium CS
RotMat_trap = [Xaxis_trap', Yaxis_trap', Zaxis_trap'];
loc2glob_trap = [RotMat_trap, O_trap';
     0 0 0 1]; 
   % These tranformation matrices convert from local to global CS.
   % Invert them to go from global to local.
glob2loc = [RotMat_trap', -(RotMat_trap')*O_trap'; 0 0 0 1];
R_glob2loc = glob2loc(1:3,1:3);
T_glob2loc = glob2loc(1:3,4);

if(trap_cs_available)
else
% Transformes data to trapezium CS
V_Trap_P1 = R_glob2loc*V_Trap_P1.' + repmat(T_glob2loc,1,size(V_Trap_P1,1));
V_Trap_P1 = V_Trap_P1.';

% Plots the trapezium with the new coordinate system for validation
    h0 = figure;
    [obj, li, ax] = GUI_PlotShells(h0, {F_Trap_P1},...
        {V_Trap_P1},...
        {ones(size(V_Trap_P1,1),1)});
    GUI_VisualisationUI(ancestor(ax, 'figure'), true, ax, true, true);
        
        % Plots the coordinate system and its origin
    hold on
    arrow([0;0;0],[30;0;0],10,70,30,5,'EdgeColor','r','FaceColor','r');
    arrow([0;0;0],[0;30;0],10,70,30,5,'EdgeColor','g','FaceColor','g');
    arrow([0;0;0],[0;0;30],10,70,30,5,'EdgeColor','b','FaceColor','b');
    
        % Plots 2 buttons: 'Bad CS' if the coordinate system looks wrong,
        % this will stop the code, and 'Continue' if the coordinate system
        % looks right, will continue to run the rest of the code
    b1 = uicontrol('Position',[20 20 200 40],'String','Continue',...
              'Callback','uiresume(gcbf)');
    b2 = uicontrol('Position',[20 60 200 40],'String','Bad CS',...
              'Callback','close(gcbf)');
    uiwait(gcf);
    close(h0);
end

% Plots the trapezium with the new coordinate system
[F_Trap_P1, V_Trap_P1] =  STL_ReadFile([dir Trap_P1],true);
% Transformes data to trapezium CS
V_Trap_P1 = R_glob2loc*V_Trap_P1.' + repmat(T_glob2loc,1,size(V_Trap_P1,1));
V_Trap_P1 = V_Trap_P1.';
    h0 = figure;
    [obj, li, ax] = GUI_PlotShells(h0, {F_Trap_P1},...
        {V_Trap_P1},...
        {ones(size(V_Trap_P1,1),1)});
    GUI_VisualisationUI(ancestor(ax, 'figure'), true, ax, true, true);
        
        % Plots the coordinate system and its origin
    hold on
    arrow([0;0;0],[30;0;0],10,70,30,5,'EdgeColor','r','FaceColor','r');
    arrow([0;0;0],[0;30;0],10,70,30,5,'EdgeColor','g','FaceColor','g');
    arrow([0;0;0],[0;0;30],10,70,30,5,'EdgeColor','b','FaceColor','b');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KINEMATICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Searches for trapezium fragment
if(data_available)
    % Loads data
    A = txt2mat(file);
        % Loads all the STL files
    [F_MC1_P1, V_MC1_P1] =  STL_ReadFile([dir MC1_P1],true);
    [F_Trap_P1, V_Trap_P1] =  STL_ReadFile([dir Trap_P1],true);
    [F_MC1_P2, V_MC1_P2] =  STL_ReadFile([dir MC1_P2]',true);
    [F_Trap_P2, V_Trap_P2] =  STL_ReadFile([dir Trap_P2],true);

        % Associates the different parts of the RotTransData file to the
        % corresponding rotation and translation matrices for each bone
    R_Trap = A(1:3,1:3);
    T_Trap = A(1:3,4);
    H_Trap = [R_Trap, T_Trap; 0 0 0 1];
    R_MC1 = A(4:6,1:3);
    T_MC1 = A(4:6,4);
    H_MC1 = [R_MC1, T_MC1; 0 0 0 1];
    
else
    % Loads data
       fid = fopen(file,'w+');
        % Loads all the STL files
    [F_MC1_P1, V_MC1_P1] =  STL_ReadFile([dir MC1_P1],true);
    [F_Trap_P1, V_Trap_P1] =  STL_ReadFile([dir Trap_P1],true);
    [F_MC1_P2, V_MC1_P2] =  STL_ReadFile([dir MC1_P2]',true);
    [F_Trap_P2, V_Trap_P2] =  STL_ReadFile([dir Trap_P2],true);
    
    % Registers trapezium in position 2 with the trapezium in position 1
    [R_Trap, T_Trap, ERROR_Trap] = RegisterBones(F_Trap_P2, V_Trap_P2, F_Trap_P1, V_Trap_P1,1);
    H_Trap = [R_Trap, T_Trap; 0 0 0 1];
    
    % Saves results to file
    W = reshape([R_Trap, T_Trap].',1,[]);
    fprintf(fid,'%5.8f\t%5.8f\t%5.8f\t%5.8f\r\n',W); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Transformes all data to trapezium CS
V_Trap_P1 = R_glob2loc*V_Trap_P1.' + repmat(T_glob2loc,1,size(V_Trap_P1,1));
V_Trap_P1 = V_Trap_P1.';
V_MC1_P1 = (glob2loc * [V_MC1_P1.' ; ones(1,size(V_MC1_P1,1))]).';
V_MC1_P1 = V_MC1_P1(:,1:3);

V_Trap_P2 = (glob2loc * H_Trap * [V_Trap_P2.' ; ones(1,size(V_Trap_P2,1))]).';
V_Trap_P2 = V_Trap_P2(:,1:3);
V_MC1_P2 = (glob2loc * H_Trap * [V_MC1_P2.' ; ones(1,size(V_MC1_P2,1))]).';
V_MC1_P2 = V_MC1_P2(:,1:3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Transformes data to local CS. Once transformed the ICP starts.

% MC1 expressed in the trapezium CS.
if(~data_available)
    [R_MC1, T_MC1, ERROR_MC1] = RegisterBones(F_MC1_P2, V_MC1_P2, F_MC1_P1, V_MC1_P1,1);
end
Hom_MC1 = [R_MC1, T_MC1; 0 0 0 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Saving the gathered info (so you don't have to re-run the ICP)
if(~data_available)
    % Transformation of MC1
    W = reshape([R_MC1, T_MC1].',1,[]);
    fprintf(fid,'%5.8f\t%5.8f\t%5.8f\t%5.8f\r\n',W);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Checks and tests

% Plots each motion relative to the trapezium CS.

    h8 = figure;
    [obj, li, ax] = GUI_PlotShells(h8, {F_Trap_P1;F_Trap_P2;F_MC1_P2;F_MC1_P1},...
        {V_Trap_P1;V_Trap_P2;V_MC1_P2;V_MC1_P1},...
        {ones(size(V_Trap_P1,1),1),ones(size(V_Trap_P2,1),1),ones(size(V_MC1_P2,1),1),ones(size(V_MC1_P1,1),1)});
    hold on
    arrow([0;0;0],[10;0;0],5,70,30);
    arrow([0;0;0],[0;10;0],5,70,30);
    arrow([0;0;0],[0;0;10],5,70,30);
    hold off;

% test: transform the first segments so it should match-up with the
% second one.
    h1 = figure;
    [obj, li, ax] = GUI_PlotShells(h1, {F_Trap_P1;F_Trap_P2}, {V_Trap_P1;V_Trap_P2},...
            {ones(size(V_Trap_P1,1),1),ones(size(V_Trap_P2,1),1)});
    box off
        
    MC1_loc2 = (Hom_MC1 * [V_MC1_P2.' ; ones(1,size(V_MC1_P2,1))]).';
    MC1_loc2 = MC1_loc2(:,1:3);
    h4 = figure;
    [obj, li, ax] = GUI_PlotShells(h4, {F_MC1_P2;F_MC1_P1}, {MC1_loc2;V_MC1_P1},...
            {ones(size(V_MC1_P2,1),1),ones(size(V_MC1_P1,1),1)});
    box off  
    
%   Calculates the helical axes from position 1 to position 2
    % MC1
[n_MC1,point_MC1,phi_MC1,t_MC1] = screw(Hom_MC1);
Hel_MC1_p = point_MC1.';
Hel_MC1_dir = n_MC1.';

% Plots STl files and the helical axes
    h5 = figure;
    [obj, li, ax] = GUI_PlotShells(h5, {F_Trap_P1;F_MC1_P1},...
        {V_Trap_P1;V_MC1_P1},...
        {ones(size(V_Trap_P1,1),1),ones(size(V_MC1_P1,1),1)});
        
        % Plots the coordinate system and its origin
    hold on
    arrow([0;0;0],[10;0;0],5,70,30);
    arrow([0;0;0],[0;10;0],5,70,30);
    arrow([0;0;0],[0;0;10],5,70,30);
        
        % Plots the helical axis for each joint
    arrow(point_MC1-200*n_MC1,point_MC1+200*n_MC1,5,70,30,'EdgeColor','b','FaceColor','b','LineWidth',4);
        
        % Saves to file
    fid4 = fopen(fileHelicalAX,'w+');
    fprintf(fid4,'%5.8f\t%5.8f\t%5.8f\r\n',Hel_MC1_p,Hel_MC1_dir);
    

% Calculation of the rotation angles along each axis of the radius
% coordinate system for each bone. Each rotation angle is expressed in the
% radius coordinate system relatively to the more proximal bone.

    % MC1
    
    Rot_MC1_x = asin(R_MC1(3,2)/cos(asin(-R_MC1(3,1))))*180/pi;
    Rot_MC1_y = asin(-R_MC1(3,1))*180/pi;
    Rot_MC1_z = asin(R_MC1(2,1)/cos(asin(-R_MC1(3,1))))*180/pi;
    Rot_MC1 = [Rot_MC1_x,Rot_MC1_y,Rot_MC1_z];
    
    % Saves to file
        fid5 = fopen(fileRotAngle,'w+');
        fprintf(fid5,'%5.8f\t%5.8f\t%5.8f\r\n',Rot_MC1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROXIMITY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FRAME 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MC1:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculation of the shortest distance between the two STL selected
[idx1,D1] = knnsearch(V_MC1_P1,V_Trap_P1);
[C1,idxD1] = min(D1);
dist_min_P1 = C1

% Calculation of the matrix D which contains the Euclidean distances between
%       each points of each STL file
D = pdist2(V_MC1_P1,V_Trap_P1);

% Define the condition to take only one part of the stl surface into
% consideration for the following steps (zone of interest)
Def1 = D < dist_max.';

%     Now we calculate the matrix idxD which corresponds to a 2xn matrix.
%     The first row gives all the indexes of the points located in the
%     zone of interest of the V_Bone1 matrix. The second row gives all the
%     indexes of the points located in the zone of interest of the V_Bone2
%     matrix.
[idxrD,idxcD] = find(D<dist_max);
idxD = [idxrD,idxcD];
idxD(:,any(idxD==0,1)) = []; % Removes all the zero in the idxD Matrix
idxD1 = unique(idxD(:,1));
idxD2 = unique(idxD(:,2));

%     Then we isolate the coordinates of each point of each STL file
%     located in the zone of interest in two different matrices.
for i = 1:length(idxD1)
    CA_V_Bone1(i,:) = V_MC1_P1(idxD1(i),:);
end
for i = 1:length(idxD2)
    CA_V_Bone2(i,:) = V_Trap_P1(idxD2(i),:);
end

%     Then we assign each point from the first STL with the corresponding
%     one from the second STL
for i = 1:length(CA_V_Bone1)
    DCA = pdist2(CA_V_Bone1(i,:),CA_V_Bone2);
    [minDCA(i),idxCA(i)] = min(DCA);
end

% Define the condition to consider the zone of proximity (red zone)
Defprox1 = D < dist_prox.';

%     Now we calculate the matrix idxD which corresponds to a 2xn matrix.
%     The first row gives all the indexes of the points located in the
%     zone of interest of the V_Bone1 matrix. The second row gives all the
%     indexes of the points located in the zone of interest of the V_Bone2
%     matrix.
[idxrDprox1,idxcDprox1] = find(D<dist_prox);
idxDprox1 = [idxrDprox1,idxcDprox1];
idxDprox1(:,any(idxDprox1==0,1)) = []; % Removes all the zero in the idxD Matrix
idxD1prox = unique(idxDprox1(:,1));
idxD2prox = unique(idxDprox1(:,2));

%     Then we isolate the coordinates of each point of each STL file
%     located in the zone of interest in two different matrices.
for i = 1:length(idxD1prox)
    CA_V_Bone1prox(i,:) = V_MC1_P1(idxD1prox(i),:);
end
for i = 1:length(idxD2prox)
    CA_V_Bone2prox(i,:) = V_Trap_P1(idxD2prox(i),:);
end

%     Then we assign each point from the first STL with the corresponding
%     one from the second STL
for i = 1:length(CA_V_Bone1prox)
    DCAprox1 = pdist2(CA_V_Bone1prox(i,:),CA_V_Bone2prox);
    [minDCAprox1(i),idxCAprox1(i)] = min(DCAprox1);
end

        % We can now test to which triangle these vertices correspond in
        % the F matrix and save their indexes
        test_F1 = ismember(F_Bone1,idxD1prox);
        idx_F1 = find(test_F1(:,1)==1 & test_F1(:,2)==1 & test_F1(:,3)==1);
            % This idx_F_trap matrix is a Nx1 matrix, where N is the amount
            % of triangles located in the deformation zone, and each value
            % in this matrix corresponds to the index of a line in the F
            % matrix that defines one of these triangles
        for i=1:length(idx_F1)
            CA_F_Bone1prox(i,:) = F_Bone1(idx_F1(i),:);
        end
        % We can now test to which triangle these vertices correspond in
        % the F matrix and save their indexes
        test_F2 = ismember(F_Bone2,idxD2prox);
        idx_F2 = find(test_F2(:,1)==1 & test_F2(:,2)==1 & test_F2(:,3)==1);
            % This idx_F_trap matrix is a Nx1 matrix, where N is the amount
            % of triangles located in the deformation zone, and each value
            % in this matrix corresponds to the index of a line in the F
            % matrix that defines one of these triangles
        for i=1:length(idx_F2)
            CA_F_Bone2prox(i,:) = F_Bone2(idx_F2(i),:);
        end

%       Each triangle is composed of 3 points called a,b and c
Xa1 = V_MC1_P1(CA_F_Bone1prox(:,1),1);
Ya1 = V_MC1_P1(CA_F_Bone1prox(:,1),2);
Za1 = V_MC1_P1(CA_F_Bone1prox(:,1),3);

Xb1 = V_MC1_P1(CA_F_Bone1prox(:,2),1);
Yb1 = V_MC1_P1(CA_F_Bone1prox(:,2),2);
Zb1 = V_MC1_P1(CA_F_Bone1prox(:,2),3);

Xc1 = V_MC1_P1(CA_F_Bone1prox(:,3),1);
Yc1 = V_MC1_P1(CA_F_Bone1prox(:,3),2);
Zc1 = V_MC1_P1(CA_F_Bone1prox(:,3),3);

%       Now we need to calculate the coordinates of the 2 vectors ab and ac
for i=1:length(Xa1)
    Xab1(i) = Xb1(i)-Xa1(i);
    Yab1(i) = Yb1(i)-Ya1(i);
    Zab1(i) = Zb1(i)-Za1(i);
    Xac1(i) = Xc1(i)-Xa1(i);
    Yac1(i) = Yc1(i)-Ya1(i);
    Zac1(i) = Zc1(i)-Za1(i);
end

%       Now we need to calculate the surface area of each triangle and sum
%       all of them to obtain to total contact area
for i = 1:length(Xab1)
    S1(i) = (1/2)*((Yab1(i)*Zac1(i)-Zab1(i)*Yac1(i))^2+(Zab1(i)*Xac1(i)-Xab1(i)*Zac1(i))^2+(Xab1(i)*Yac1(i)-Yab1(i)*Xac1(i))^2)^(1/2);
end
PA_MC1_P1 = sum(S1(:))

%     Now we can plot those points located in the compression area directly
%     on the STL files in order to visualize the contact area
XB1 = CA_V_Bone1(:,1);
YB1 = CA_V_Bone1(:,2);
ZB1 = CA_V_Bone1(:,3);
XB2 = CA_V_Bone2(:,1);
YB2 = CA_V_Bone2(:,2);
ZB2 = CA_V_Bone2(:,3);

% Plots the MC1 with the color map representing the proximity
% pattern
area1 = figure; 
    [obj, li, ax] = GUI_PlotShells(area1, {F_MC1_P1}, {V_MC1_P1},...
            {ones(size(V_MC1_P1,1),1)},[0,0,1]);
box off
view(-135,60);
for i=1:length(V_MC1_P1)
    colrstl1(i,:)=[0,0,1];
end

hold on
for i=1:length(idxD1)
    Dist1(i)=((XB1(i)-XB2(idxCA(i)))^2+(YB1(i)-YB2(idxCA(i)))^2+(ZB1(i)-ZB2(idxCA(i)))^2)^(1/2);
    if Dist1(i) <= 1
       colrstl1(idxD1(i),:) = [1,0,0];
    elseif Dist1(i) > 1   
    colrstl1(idxD1(i),:)= ...
        [-((0.5*(Dist1(i)+1)-1)^2)+1,...
        -((1.5*(Dist1(i)-1)-dist_max/2)^2)/((dist_max/2)^2)+1,...
        -((0.5*(Dist1(i)-1)-1)^2)+1];
    end
end
hold on
patch('Faces',F_MC1_P1,'Vertices',V_MC1_P1, ...
    'FaceColor','interp', ...
    'FaceVertexCData',colrstl1, ...
    'EdgeColor', 'interp', ...
    'EdgeAlpha', 0, ...
    'CDataMapping', 'scaled',...
    'AmbientStrength', 0.4, ...
    'DiffuseStrength', 0.8, ...
    'SpecularStrength', 0.2, ...
    'SpecularColorReflectance', 0.5, ...
    'FaceLighting', 'gouraud')
hold on
step_h = (dist_max-1)/8;
scale_h = [dist_max,dist_max-step_h,dist_max-2*step_h,dist_max-3*step_h,...
    dist_max-4*step_h,dist_max-5*step_h,dist_max-6*step_h,dist_max-7*step_h,dist_max-8*step_h];
h=colorbar('YTickLabel',{scale_h});
set(h,'ylim',[0.1 0.9]);
title = get(h,'Title');
titleString = 'Distance (mm)';
set(title ,'String',titleString,'FontWeight','bold');
set(gcf,'numbertitle','off','name','MC1 - Proximity pattern (Position 1)');

% Trapezium:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read STL files
[F_Trap_P12, V_Trap_P12] =  STL_ReadFile([dir Trap_P1],true);
[F_MC1_P12, V_MC1_P12] =  STL_ReadFile([dir MC1_P1],true);

% Calculation of the matrix D which contains the Euclidean distances between
%       each points of each STL file
D3 = pdist2(V_Trap_P12,V_MC1_P12);

% Define the condition to take only one part of the stl surface into
% consideration for the following steps (zone of interest)
Def3 = D3 < dist_max.';

%     Now we calculate the matrix idxD which corresponds to a 2xn matrix.
%     The first row gives all the indexes of the points located in the
%     zone of interest of the V_Bone1 matrix. The second row gives all the
%     indexes of the points located in the zone of interest of the V_Bone2
%     matrix.
[idxrDt,idxcDt] = find(D3<dist_max);
idxDt = [idxrDt,idxcDt];
idxDt(:,any(idxDt==0,1)) = []; % Removes all the zero in the idxD Matrix
idxD3 = unique(idxDt(:,1));
idxD4 = unique(idxDt(:,2));

%     Then we isolate the coordinates of each point of each STL file
%     located in the zone of interest in two different matrices.
for i = 1:length(idxD3)
    CA_V_Bone3(i,:) = V_Trap_P12(idxD3(i),:);
end
for i = 1:length(idxD4)
    CA_V_Bone4(i,:) = V_MC1_P12(idxD4(i),:);
end

%     Then we assign each point from the first STL with the corresponding
%     one from the second STL
for i = 1:length(CA_V_Bone3)
    DCAt = pdist2(CA_V_Bone3(i,:),CA_V_Bone4);
    [minDCAt(i),idxCAt(i)] = min(DCAt);
end

% Define the condition to consider the zone of proximity (red zone)
Def3prox = D3 < dist_prox.';

%     Now we calculate the matrix idxD which corresponds to a 2xn matrix.
%     The first row gives all the indexes of the points located in the
%     zone of interest of the V_Bone1 matrix. The second row gives all the
%     indexes of the points located in the zone of interest of the V_Bone2
%     matrix.
[idxrDtprox,idxcDtprox] = find(D3<dist_prox);
idxDtprox = [idxrDtprox,idxcDtprox];
idxDtprox(:,any(idxDtprox==0,1)) = []; % Removes all the zero in the idxD Matrix
idxD3prox = unique(idxDtprox(:,1));
idxD4prox = unique(idxDtprox(:,2));

%     Then we isolate the coordinates of each point of each STL file
%     located in the zone of interest in two different matrices.
for i = 1:length(idxD3prox)
    CA_V_Bone3prox(i,:) = V_Trap_P12(idxD3prox(i),:);
end
for i = 1:length(idxD4prox)
    CA_V_Bone4prox(i,:) = V_MC1_P12(idxD4prox(i),:);
end

%     Then we assign each point from the first STL with the corresponding
%     one from the second STL
for i = 1:length(CA_V_Bone3prox)
    DCAtprox = pdist2(CA_V_Bone3prox(i,:),CA_V_Bone4prox);
    [minDCAtprox(i),idxCAprox2(i)] = min(DCAtprox);
end

        % We can now test to which triangle these vertices correspond in
        % the F matrix and save their indexes
        test_F3 = ismember(F_Bone3,idxD3prox);
        idx_F3 = find(test_F3(:,1)==1 & test_F3(:,2)==1 & test_F3(:,3)==1);
            % This idx_F_trap matrix is a Nx1 matrix, where N is the amount
            % of triangles located in the deformation zone, and each value
            % in this matrix corresponds to the index of a line in the F
            % matrix that defines one of these triangles
        for i=1:length(idx_F3)
            CA_F_Bone3prox(i,:) = F_Bone3(idx_F3(i),:);
        end
        % We can now test to which triangle these vertices correspond in
        % the F matrix and save their indexes
        test_F4 = ismember(F_Bone4,idxD4prox);
        idx_F4 = find(test_F4(:,1)==1 & test_F4(:,2)==1 & test_F4(:,3)==1);
            % This idx_F_trap matrix is a Nx1 matrix, where N is the amount
            % of triangles located in the deformation zone, and each value
            % in this matrix corresponds to the index of a line in the F
            % matrix that defines one of these triangles
        for i=1:length(idx_F4)
            CA_F_Bone4prox(i,:) = F_Bone4(idx_F4(i),:);
        end

%       Each triangle is composed of 3 points called a,b and c
Xa3 = V_Trap_P12(CA_F_Bone3prox(:,1),1);
Ya3 = V_Trap_P12(CA_F_Bone3prox(:,1),2);
Za3 = V_Trap_P12(CA_F_Bone3prox(:,1),3);

Xb3 = V_Trap_P12(CA_F_Bone3prox(:,2),1);
Yb3 = V_Trap_P12(CA_F_Bone3prox(:,2),2);
Zb3 = V_Trap_P12(CA_F_Bone3prox(:,2),3);

Xc3 = V_Trap_P12(CA_F_Bone3prox(:,3),1);
Yc3 = V_Trap_P12(CA_F_Bone3prox(:,3),2);
Zc3 = V_Trap_P12(CA_F_Bone3prox(:,3),3);

%       Now we need to calculate the coordinates of the 2 vectors ab and ac
for i=1:length(Xa3)
    Xab3(i) = Xb3(i)-Xa3(i);
    Yab3(i) = Yb3(i)-Ya3(i);
    Zab3(i) = Zb3(i)-Za3(i);
    Xac3(i) = Xc3(i)-Xa3(i);
    Yac3(i) = Yc3(i)-Ya3(i);
    Zac3(i) = Zc3(i)-Za3(i);
end

%       Now we need to calculate the surface area of each triangle and sum
%       all of them to obtain to total contact area
for i = 1:length(Xab3)
    S3(i) = (1/2)*((Yab3(i)*Zac3(i)-Zab3(i)*Yac3(i))^2+(Zab3(i)*Xac3(i)-Xab3(i)*Zac3(i))^2+(Xab3(i)*Yac3(i)-Yab3(i)*Xac3(i))^2)^(1/2);
end
PA_Trap_P1 = sum(S3(:))

%     Now we can plot those points located in the compression area directly
%     on the STL files in order to visualize the contact area
XB3 = CA_V_Bone3(:,1);
YB3 = CA_V_Bone3(:,2);
ZB3 = CA_V_Bone3(:,3);
XB4 = CA_V_Bone4(:,1);
YB4 = CA_V_Bone4(:,2);
ZB4 = CA_V_Bone4(:,3);

% Plots the trapezium with the color map representing the proximity
% pattern
area2 = figure; 
    [obj, li, ax] = GUI_PlotShells(area2, {F_Trap_P12}, {V_Trap_P12},...
            {ones(size(V_Trap_P12,1),1)},[0,0,1]);
box off
view(0,-45);
for i=1:length(V_Trap_P12)
    colrstl2(i,:)=[0,0,1];
end

hold on
for i=1:length(idxD3)
    Dist2(i)=((XB3(i)-XB4(idxCAt(i)))^2+(YB3(i)-YB4(idxCAt(i)))^2+(ZB3(i)-ZB4(idxCAt(i)))^2)^(1/2);
    if Dist2(i) <= 1
       colrstl2(idxD3(i),:) = [1,0,0];
    elseif Dist2(i) > 1   
    colrstl2(idxD3(i),:)= ...
        [-((0.5*(Dist2(i)+1)-1)^2)+1,...
        -((1.5*(Dist2(i)-1)-dist_max/2)^2)/((dist_max/2)^2)+1,...
        -((0.5*(Dist2(i)-1)-1)^2)+1];
    end
end
hold on
patch('Faces',F_Trap_P12,'Vertices',V_Trap_P12, ...
    'FaceColor','interp', ...
    'FaceVertexCData',colrstl2, ...
    'EdgeColor', 'interp', ...
    'EdgeAlpha', 0, ...
    'CDataMapping', 'scaled',...
    'AmbientStrength', 0.4, ...
    'DiffuseStrength', 0.8, ...
    'SpecularStrength', 0.2, ...
    'SpecularColorReflectance', 0.5, ...
    'FaceLighting', 'gouraud')
hold on
step_h = (dist_max-1)/8;
scale_h = [dist_max,dist_max-step_h,dist_max-2*step_h,dist_max-3*step_h,...
    dist_max-4*step_h,dist_max-5*step_h,dist_max-6*step_h,dist_max-7*step_h,dist_max-8*step_h];
h=colorbar('YTickLabel',{scale_h});
set(h,'ylim',[0.1 0.9]);
title = get(h,'Title');
titleString = 'Distance (mm)';
set(title ,'String',titleString,'FontWeight','bold');
set(gcf,'numbertitle','off','name','Trapezium - Proximity pattern (Position 1)');


% FRAME 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MC1:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculation of the shortest distance between the two STL selected
[idx1b,D1b] = knnsearch(V_MC1_P2,V_Trap_P2);
[C1b,idxD1b] = min(D1b);
dist_min_P2 = C1b

% Calculation of the matrix D which contains the Euclidean distances between
%       each points of each STL file
Db = pdist2(V_MC1_P2,V_Trap_P2);

% Define the condition to take only one part of the stl surface into
% consideration for the following steps (zone of interest)
Def1b = Db < dist_max.';

%     Now we calculate the matrix idxD which corresponds to a 2xn matrix.
%     The first row gives all the indexes of the points located in the
%     zone of interest of the V_Bone1 matrix. The second row gives all the
%     indexes of the points located in the zone of interest of the V_Bone2
%     matrix.
[idxrDb,idxcDb] = find(Db<dist_max);
idxDb = [idxrDb,idxcDb];
idxDb(:,any(idxDb==0,1)) = []; % Removes all the zero in the idxD Matrix
idxD1b = unique(idxDb(:,1));
idxD2b = unique(idxDb(:,2));

%     Then we isolate the coordinates of each point of each STL file
%     located in the zone of interest in two different matrices.
for i = 1:length(idxD1b)
    CA_V_Bone1b(i,:) = V_MC1_P2(idxD1b(i),:);
end
for i = 1:length(idxD2b)
    CA_V_Bone2b(i,:) = V_Trap_P2(idxD2b(i),:);
end

%     Then we assign each point from the first STL with the corresponding
%     one from the second STL
for i = 1:length(CA_V_Bone1b)
    DCAb = pdist2(CA_V_Bone1b(i,:),CA_V_Bone2b);
    [minDCAb(i),idxCAb(i)] = min(DCAb);
end

% Define the condition to consider the zone of proximity (red zone)
Defprox1b = Db < dist_prox.';

%     Now we calculate the matrix idxD which corresponds to a 2xn matrix.
%     The first row gives all the indexes of the points located in the
%     zone of interest of the V_Bone1 matrix. The second row gives all the
%     indexes of the points located in the zone of interest of the V_Bone2
%     matrix.
[idxrDprox1b,idxcDprox1b] = find(Db<dist_prox);
idxDprox1b = [idxrDprox1b,idxcDprox1b];
idxDprox1b(:,any(idxDprox1b==0,1)) = []; % Removes all the zero in the idxD Matrix
idxD1proxb = unique(idxDprox1b(:,1));
idxD2proxb = unique(idxDprox1b(:,2));

%     Then we isolate the coordinates of each point of each STL file
%     located in the zone of interest in two different matrices.
for i = 1:length(idxD1proxb)
    CA_V_Bone1proxb(i,:) = V_MC1_P2(idxD1proxb(i),:);
end
for i = 1:length(idxD2proxb)
    CA_V_Bone2proxb(i,:) = V_Trap_P2(idxD2proxb(i),:);
end

%     Then we assign each point from the first STL with the corresponding
%     one from the second STL
for i = 1:length(CA_V_Bone1proxb)
    DCAprox1b = pdist2(CA_V_Bone1proxb(i,:),CA_V_Bone2proxb);
    [minDCAprox1b(i),idxCAprox1b(i)] = min(DCAprox1b);
end

% Now we can calculate the proximity area by calculating each
% triangle's area formed by the points located in the proximity zone
for i = 1:length(idxD1proxb)
    CA_F_Bone1proxb(i,:) = F_MC1_P2(idxD1proxb(i),:);
end
for i = 1:length(idxD2proxb)
    CA_F_Bone2proxb(i,:) = F_Trap_P2(idxD2proxb(i),:);
end

%       Each triangle is composed of 3 points called a,b and c
Xa1b = V_MC1_P2(CA_F_Bone1proxb(:,1),1);
Ya1b = V_MC1_P2(CA_F_Bone1proxb(:,1),2);
Za1b = V_MC1_P2(CA_F_Bone1proxb(:,1),3);

Xb1b = V_MC1_P2(CA_F_Bone1proxb(:,2),1);
Yb1b = V_MC1_P2(CA_F_Bone1proxb(:,2),2);
Zb1b = V_MC1_P2(CA_F_Bone1proxb(:,2),3);

Xc1b = V_MC1_P2(CA_F_Bone1proxb(:,3),1);
Yc1b = V_MC1_P2(CA_F_Bone1proxb(:,3),2);
Zc1b = V_MC1_P2(CA_F_Bone1proxb(:,3),3);

%       Now we need to calculate the coordinates of the 2 vectors ab and ac
for i=1:length(Xa1b)
    Xab1b(i) = Xb1b(i)-Xa1b(i);
    Yab1b(i) = Yb1b(i)-Ya1b(i);
    Zab1b(i) = Zb1b(i)-Za1b(i);
    Xac1b(i) = Xc1b(i)-Xa1b(i);
    Yac1b(i) = Yc1b(i)-Ya1b(i);
    Zac1b(i) = Zc1b(i)-Za1b(i);
end

%       Now we need to calculate the surface area of each triangle and sum
%       all of them to obtain to total contact area
for i = 1:length(Xab1b)
    S1b(i) = (1/2)*((Yab1b(i)*Zac1b(i)-Zab1b(i)*Yac1b(i))^2+(Zab1b(i)*Xac1b(i)-Xab1b(i)*Zac1b(i))^2+(Xab1b(i)*Yac1b(i)-Yab1b(i)*Xac1b(i))^2)^(1/2);
end
PA_MC1_P2 = sum(S1b(:))

%     Now we can plot those points located in the compression area directly
%     on the STL files in order to visualize the contact area
XB1b = CA_V_Bone1b(:,1);
YB1b = CA_V_Bone1b(:,2);
ZB1b = CA_V_Bone1b(:,3);
XB2b = CA_V_Bone2b(:,1);
YB2b = CA_V_Bone2b(:,2);
ZB2b = CA_V_Bone2b(:,3);

% Plots the MC1 with the color map representing the proximity
% pattern
area1b = figure; 
    [obj, li, ax] = GUI_PlotShells(area1b, {F_MC1_P2}, {V_MC1_P2},...
            {ones(size(V_MC1_P2,1),1)},[0,0,1]);
box off
view(-135,60);
for i=1:length(V_MC1_P2)
    colrstl1b(i,:)=[0,0,1];
end

hold on
for i=1:length(idxD1b)
    Dist1b(i)=((XB1b(i)-XB2b(idxCAb(i)))^2+(YB1b(i)-YB2b(idxCAb(i)))^2+(ZB1b(i)-ZB2b(idxCAb(i)))^2)^(1/2);
    if Dist1b(i) <= 1
       colrstl1b(idxD1b(i),:) = [1,0,0];
    elseif Dist1b(i) > 1   
    colrstl1b(idxD1b(i),:)= ...
        [-((0.5*(Dist1b(i)+1)-1)^2)+1,...
        -((1.5*(Dist1b(i)-1)-dist_max/2)^2)/((dist_max/2)^2)+1,...
        -((0.5*(Dist1b(i)-1)-1)^2)+1];
    end
end
hold on
patch('Faces',F_MC1_P2,'Vertices',V_MC1_P2, ...
    'FaceColor','interp', ...
    'FaceVertexCData',colrstl1b, ...
    'EdgeColor', 'interp', ...
    'EdgeAlpha', 0, ...
    'CDataMapping', 'scaled',...
    'AmbientStrength', 0.4, ...
    'DiffuseStrength', 0.8, ...
    'SpecularStrength', 0.2, ...
    'SpecularColorReflectance', 0.5, ...
    'FaceLighting', 'gouraud')
hold on
step_h = (dist_max-1)/8;
scale_h = [dist_max,dist_max-step_h,dist_max-2*step_h,dist_max-3*step_h,...
    dist_max-4*step_h,dist_max-5*step_h,dist_max-6*step_h,dist_max-7*step_h,dist_max-8*step_h];
h=colorbar('YTickLabel',{scale_h});
set(h,'ylim',[0.1 0.9]);
title = get(h,'Title');
titleString = 'Distance (mm)';
set(title ,'String',titleString,'FontWeight','bold');
set(gcf,'numbertitle','off','name','MC1 - Proximity pattern (Frame 2)');

% Trapezium:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read STL files
[F_Trap_P22, V_Trap_P22] =  STL_ReadFile([dir Trap_P1],true);
[F_MC1_P22, V_MC1_P22] =  STL_ReadFile([dir MC1_P1],true);

% Calculation of the matrix D which contains the Euclidean distances between
%       each points of each STL file
D3b = pdist2(V_Trap_P22,V_MC1_P22);

% Define the condition to take only one part of the stl surface into
% consideration for the following steps (zone of interest)
Def3b = D3b < dist_max.';

%     Now we calculate the matrix idxD which corresponds to a 2xn matrix.
%     The first row gives all the indexes of the points located in the
%     zone of interest of the V_Bone1 matrix. The second row gives all the
%     indexes of the points located in the zone of interest of the V_Bone2
%     matrix.
[idxrDtb,idxcDtb] = find(D3b<dist_max);
idxDtb = [idxrDtb,idxcDtb];
idxDtb(:,any(idxDtb==0,1)) = []; % Removes all the zero in the idxD Matrix
idxD3b = unique(idxDtb(:,1));
idxD4b = unique(idxDtb(:,2));

%     Then we isolate the coordinates of each point of each STL file
%     located in the zone of interest in two different matrices.
for i = 1:length(idxD3b)
    CA_V_Bone3b(i,:) = V_Trap_P22(idxD3b(i),:);
end
for i = 1:length(idxD4b)
    CA_V_Bone4b(i,:) = V_MC1_P22(idxD4b(i),:);
end

%     Then we assign each point from the first STL with the corresponding
%     one from the second STL
for i = 1:length(CA_V_Bone3b)
    DCAtb = pdist2(CA_V_Bone3b(i,:),CA_V_Bone4b);
    [minDCAtb(i),idxCAtb(i)] = min(DCAtb);
end

% Define the condition to consider the zone of proximity (red zone)
Def3proxb = D3b < dist_prox.';

%     Now we calculate the matrix idxD which corresponds to a 2xn matrix.
%     The first row gives all the indexes of the points located in the
%     zone of interest of the V_Bone1 matrix. The second row gives all the
%     indexes of the points located in the zone of interest of the V_Bone2
%     matrix.
[idxrDtproxb,idxcDtproxb] = find(D3b<dist_prox);
idxDtproxb = [idxrDtproxb,idxcDtproxb];
idxDtproxb(:,any(idxDtproxb==0,1)) = []; % Removes all the zero in the idxD Matrix
idxD3proxb = unique(idxDtproxb(:,1));
idxD4proxb = unique(idxDtproxb(:,2));

%     Then we isolate the coordinates of each point of each STL file
%     located in the zone of interest in two different matrices.
for i = 1:length(idxD3proxb)
    CA_V_Bone3proxb(i,:) = V_Trap_P22(idxD3proxb(i),:);
end
for i = 1:length(idxD4proxb)
    CA_V_Bone4proxb(i,:) = V_MC1_P22(idxD4proxb(i),:);
end

%     Then we assign each point from the first STL with the corresponding
%     one from the second STL
for i = 1:length(CA_V_Bone3proxb)
    DCAtproxb = pdist2(CA_V_Bone3proxb(i,:),CA_V_Bone4proxb);
    [minDCAtproxb(i),idxCAprox2b(i)] = min(DCAtproxb);
end

% Now we can calculate the contact surface area by calculating each
% triangle's area formed by the points located in the contact zone
for i = 1:length(idxD3proxb)
    CA_F_Bone3proxb(i,:) = F_Trap_P22(idxD3proxb(i),:);
end
for i = 1:length(idxD4proxb)
    CA_F_Bone4proxb(i,:) = F_MC1_P22(idxD4proxb(i),:);
end

%       Each triangle is composed of 3 points called a,b and c
Xa3b = V_Trap_P22(CA_F_Bone3proxb(:,1),1);
Ya3b = V_Trap_P22(CA_F_Bone3proxb(:,1),2);
Za3b = V_Trap_P22(CA_F_Bone3proxb(:,1),3);

Xb3b = V_Trap_P22(CA_F_Bone3proxb(:,2),1);
Yb3b = V_Trap_P22(CA_F_Bone3proxb(:,2),2);
Zb3b = V_Trap_P22(CA_F_Bone3proxb(:,2),3);

Xc3b = V_Trap_P22(CA_F_Bone3proxb(:,3),1);
Yc3b = V_Trap_P22(CA_F_Bone3proxb(:,3),2);
Zc3b = V_Trap_P22(CA_F_Bone3proxb(:,3),3);

%       Now we need to calculate the coordinates of the 2 vectors ab and ac
for i=1:length(Xa3b)
    Xab3b(i) = Xb3b(i)-Xa3b(i);
    Yab3b(i) = Yb3b(i)-Ya3b(i);
    Zab3b(i) = Zb3b(i)-Za3b(i);
    Xac3b(i) = Xc3b(i)-Xa3b(i);
    Yac3b(i) = Yc3b(i)-Ya3b(i);
    Zac3b(i) = Zc3b(i)-Za3b(i);
end

%       Now we need to calculate the surface area of each triangle and sum
%       all of them to obtain to total contact area
for i = 1:length(Xab3b)
    S3b(i) = (1/2)*((Yab3b(i)*Zac3b(i)-Zab3b(i)*Yac3b(i))^2+(Zab3b(i)*Xac3b(i)-Xab3b(i)*Zac3b(i))^2+(Xab3b(i)*Yac3b(i)-Yab3b(i)*Xac3b(i))^2)^(1/2);
end
PA_Trap_P2 = sum(S3b(:))

%     Now we can plot those points located in the compression area directly
%     on the STL files in order to visualize the contact area
XB3b = CA_V_Bone3b(:,1);
YB3b = CA_V_Bone3b(:,2);
ZB3b = CA_V_Bone3b(:,3);
XB4b = CA_V_Bone4b(:,1);
YB4b = CA_V_Bone4b(:,2);
ZB4b = CA_V_Bone4b(:,3);

% Plots the trapezium with the color map representing the proximity
% pattern
area2b = figure; 
    [obj, li, ax] = GUI_PlotShells(area2b, {F_Trap_P22}, {V_Trap_P22},...
            {ones(size(V_Trap_P22,1),1)},[0,0,1]);
box off
view(0,-45);
for i=1:length(V_Trap_P22)
    colrstl2b(i,:)=[0,0,1];
end

hold on
for i=1:length(idxD3b)
    Dist2b(i)=((XB3b(i)-XB4b(idxCAtb(i)))^2+(YB3b(i)-YB4b(idxCAtb(i)))^2+(ZB3b(i)-ZB4b(idxCAtb(i)))^2)^(1/2);
    if Dist2b(i) <= 1
       colrstl2b(idxD3b(i),:) = [1,0,0];
    elseif Dist2b(i) > 1   
    colrstl2b(idxD3b(i),:)= ...
        [-((0.5*(Dist2b(i)+1)-1)^2)+1,...
        -((1.5*(Dist2b(i)-1)-dist_max/2)^2)/((dist_max/2)^2)+1,...
        -((0.5*(Dist2b(i)-1)-1)^2)+1];
    end
end
hold on
patch('Faces',F_Trap_P22,'Vertices',V_Trap_P22, ...
    'FaceColor','interp', ...
    'FaceVertexCData',colrstl2b, ...
    'EdgeColor', 'interp', ...
    'EdgeAlpha', 0, ...
    'CDataMapping', 'scaled',...
    'AmbientStrength', 0.4, ...
    'DiffuseStrength', 0.8, ...
    'SpecularStrength', 0.2, ...
    'SpecularColorReflectance', 0.5, ...
    'FaceLighting', 'gouraud')
hold on
step_h = (dist_max-1)/8;
scale_h = [dist_max,dist_max-step_h,dist_max-2*step_h,dist_max-3*step_h,...
    dist_max-4*step_h,dist_max-5*step_h,dist_max-6*step_h,dist_max-7*step_h,dist_max-8*step_h];
h=colorbar('YTickLabel',{scale_h});
set(h,'ylim',[0.1 0.9]);
title = get(h,'Title');
titleString = 'Distance (mm)';
set(title ,'String',titleString,'FontWeight','bold');
set(gcf,'numbertitle','off','name','Trapezium - Proximity pattern (Position 2)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%