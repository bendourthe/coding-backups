clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ColorMap_Dist:
% - Calculates the shortest distance between two STL
% - Plots each STL with an arrow joining the two closest neighbours
% - Plots each STL with a colormap indicating the relative distance between
% each articular surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dependencies:        
%               STL_ReadFile.m 
%                   TRI_RemoveInvalidTriangles.m
%                   TRI_RemoveBadlyConnectedTriangles.m
%               vectarrow.m
%               GUI_PlotShells
% Input: 
%           Bone1: STL file corresponding to the first bone
%           Bone2: STL file corresponding to second bone
% 
% Output:
%           Shortest_Distance: Shortest distance calculated between the two
%                               STL files selected (mm)
%           CSmax: Maximal contact stress value (MPa)
%           Articular_AreaMC1: Contact area of the first metacarpal (mm2)
%           Articular_AreaTrap: Contact area of the trapezium (mm2)
%           Figure 1: Shows the 2 bones with an arrow linking the 2 closest
%                     neighbours
%           Figure 2: Shows the colormap of the 1st bone
%           Figure 3: Shows the colormap of the 2nd bone     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Selection of the current directory (where the STL files are)
dir = 'F:\PhD KU Leuven\Data\Study_Brown_Rig\New Segmentation\Scan18\STL\';

% FILES
Bone1 = 'SCAN18_J80_mc1.stl';
Bone2 = 'SCAN18_J80_trap.stl';

% Setting the maximal distance that the code will consider between the 2
% bones and the distance that we consider critical (underneath which we can
% expect high amount of stress due to high compression)
dist_max = 3;       % in mm
dist_crit =   0.044938;    % in mm

% MC1:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read STL files
[F_Bone1, V_Bone1] =  STL_ReadFile([dir Bone1],true);
[F_Bone2, V_Bone2] =  STL_ReadFile([dir Bone2],true);

% Calculation of the shortest distance between the two STL selected
[idx1,D1] = knnsearch(V_Bone1,V_Bone2);
[C1,idxD1] = min(D1);
Shortest_Distance = C1

% Plot the two STL selected with an arrow joining the two points
%       calculated as the closest neighbors
[idx2,D2] = knnsearch(V_Bone2,V_Bone1);
[C2,idxD2] = min(D2);

P1 = V_Bone1(idx1(idxD1),:);
P2 = V_Bone2(idx2(idxD2),:);
SD = ((P2(1)-P1(1))^2+(P2(2)-P1(2))^2+(P2(3)-P1(3))^2)^(1/2);

stl = figure;
    [obj, li, ax] = GUI_PlotShells(stl, {F_Bone1;F_Bone2}, {V_Bone1;V_Bone2},...
            {ones(size(V_Bone1,1),1),ones(size(V_Bone2,1),1)});
box off

% Calculation of the matrix D which contains the Euclidean distances between
%       each points of each STL file
D = pdist2(V_Bone1,V_Bone2);

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
    CA_V_Bone1(i,:) = V_Bone1(idxD1(i),:);
end
for i = 1:length(idxD2)
    CA_V_Bone2(i,:) = V_Bone2(idxD2(i),:);
end

%     Then we assign each point from the first STL with the corresponding
%     one from the second STL
for i = 1:length(CA_V_Bone1)
    DCA = pdist2(CA_V_Bone1(i,:),CA_V_Bone2);
    [minDCA(i),idxCA(i)] = min(DCA);
end

% Now we can calculate the contact surface area by calculating each
% triangle's area formed by the points located in the contact zone
for i = 1:length(idxD1)
    CA_F_Bone1(i,:) = F_Bone1(idxD1(i),:);
end
for i = 1:length(idxD2)
    CA_F_Bone2(i,:) = F_Bone2(idxD2(i),:);
end

% Define the condition to consider the zone of proximity (red zone)
Defprox1 = D < dist_crit.';

%     Now we calculate the matrix idxD which corresponds to a 2xn matrix.
%     The first row gives all the indexes of the points located in the
%     zone of interest of the V_Bone1 matrix. The second row gives all the
%     indexes of the points located in the zone of interest of the V_Bone2
%     matrix.
[idxrDprox1,idxcDprox1] = find(D<dist_max);
idxDprox1 = [idxrDprox1,idxcDprox1];
idxDprox1(:,any(idxDprox1==0,1)) = []; % Removes all the zero in the idxD Matrix
idxD1prox = unique(idxDprox1(:,1));
idxD2prox = unique(idxDprox1(:,2));

%     Then we isolate the coordinates of each point of each STL file
%     located in the zone of interest in two different matrices.
for i = 1:length(idxD1prox)
    CA_V_Bone1prox(i,:) = V_Bone1(idxD1prox(i),:);
end
for i = 1:length(idxD2prox)
    CA_V_Bone2prox(i,:) = V_Bone2(idxD2prox(i),:);
end

%     Then we assign each point from the first STL with the corresponding
%     one from the second STL
for i = 1:length(CA_V_Bone1prox)
    DCAprox1 = pdist2(CA_V_Bone1prox(i,:),CA_V_Bone2prox);
    [minDCAprox1(i),idxCAprox1(i)] = min(DCAprox1);
end

% Now we can calculate the proximity area by calculating each
% triangle's area formed by the points located in the proximity zone
for i = 1:length(idxD1prox)
    CA_F_Bone1prox(i,:) = F_Bone1(idxD1prox(i),:);
end
for i = 1:length(idxD2prox)
    CA_F_Bone2prox(i,:) = F_Bone2(idxD2prox(i),:);
end

%       Each triangle is composed of 3 points called a,b and c
Xa1 = V_Bone1(CA_F_Bone1prox(:,1),1);
Ya1 = V_Bone1(CA_F_Bone1prox(:,1),2);
Za1 = V_Bone1(CA_F_Bone1prox(:,1),3);

Xb1 = V_Bone1(CA_F_Bone1prox(:,2),1);
Yb1 = V_Bone1(CA_F_Bone1prox(:,2),2);
Zb1 = V_Bone1(CA_F_Bone1prox(:,2),3);

Xc1 = V_Bone1(CA_F_Bone1prox(:,3),1);
Yc1 = V_Bone1(CA_F_Bone1prox(:,3),2);
Zc1 = V_Bone1(CA_F_Bone1prox(:,3),3);

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

Articular_AreaMC1 = sum(S1(:))

%     Now we can plot those points located in the compression area directly
%     on the STL files in order to visualize the contact area
XB1 = CA_V_Bone1(:,1);
YB1 = CA_V_Bone1(:,2);
ZB1 = CA_V_Bone1(:,3);
XB2 = CA_V_Bone2(:,1);
YB2 = CA_V_Bone2(:,2);
ZB2 = CA_V_Bone2(:,3);

% Plots the STL Bone1 with the color map representing the stress
% distribution in the compression zone
area1 = figure; 
    [obj, li, ax] = GUI_PlotShells(area1, {F_Bone1}, {V_Bone1},...
            {ones(size(V_Bone1,1),1)},[0,0,1]);
box off

for i=1:length(V_Bone1)
    colrstl1(i,:)=[0,0,1];
end

hold on
for i=1:length(idxD1)
    Dist1(i)=((XB1(i)-XB2(idxCA(i)))^2+(YB1(i)-YB2(idxCA(i)))^2+(ZB1(i)-ZB2(idxCA(i)))^2)^(1/2);
       colrstl1(idxD1(i),:)= ...
        [-(Dist1(i)^2)/(dist_max^2)+1,...
        -((Dist1(i)-(dist_max+dist_crit)/2)^2)/(((dist_max-dist_crit)/2)^2)+1,...
        -((Dist1(i)-(dist_max+dist_crit))^2)/((dist_max)^2)+1];
end
hold on
patch('Faces',F_Bone1,'Vertices',V_Bone1, ...
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
step_h = (dist_max-dist_crit)/8;
scale_h = [dist_max,dist_max-step_h,dist_max-2*step_h,dist_max-3*step_h,...
    dist_max-4*step_h,dist_max-5*step_h,dist_max-6*step_h,dist_max-7*step_h,dist_max-8*step_h];
h=colorbar('YTickLabel',{scale_h});
set(h,'ylim',[0.1 0.9]);
title = get(h,'Title');
titleString = 'Distance (mm)';
set(title ,'String',titleString,'FontWeight','bold');

areaprox = figure; 
    [obj, li, ax] = GUI_PlotShells(areaprox, {F_Bone1}, {V_Bone1},...
            {ones(size(V_Bone1,1),1)},[0,0,1]);
box off

%   Plots the MC1 with only the critical proximity area

XB1prox = CA_V_Bone1prox(:,1);
YB1prox = CA_V_Bone1prox(:,2);
ZB1prox = CA_V_Bone1prox(:,3);
XB2prox = CA_V_Bone2prox(:,1);
YB2prox = CA_V_Bone2prox(:,2);
ZB2prox = CA_V_Bone2prox(:,3);

for i=1:length(V_Bone1)
    colrstlprox1(i,:)=[0,0,1];
end

hold on
for i=1:length(idxD1prox)
    Distprox1(i)=((XB1prox(i)-XB2prox(idxCAprox1(i)))^2+(YB1prox(i)-YB2prox(idxCAprox1(i)))^2 ...
    +(ZB1prox(i)-ZB2prox(idxCAprox1(i)))^2)^(1/2);
       colrstlprox1(idxD1prox(i),:)= ...
        [-(Distprox1(i)^2)/(dist_max^2)+1,...
        -((Distprox1(i)-(dist_max+dist_crit)/2)^2)/(((dist_max-dist_crit)/2)^2)+1,...
        -((Distprox1(i)-(dist_max+dist_crit))^2)/((dist_max)^2)+1];
end
for i=1:length(idxD1prox)
    colrstlprox1(idxD1prox(i),:)=[1,0,0];
end
hold on
patch('Faces',F_Bone1,'Vertices',V_Bone1, ...
    'FaceColor','interp', ...
    'FaceVertexCData',colrstlprox1, ...
    'EdgeColor', 'interp', ...
    'EdgeAlpha', 0, ...
    'CDataMapping', 'scaled',...
    'AmbientStrength', 0.4, ...
    'DiffuseStrength', 0.8, ...
    'SpecularStrength', 0.2, ...
    'SpecularColorReflectance', 0.5, ...
    'FaceLighting', 'gouraud')
hold on
step_h = (dist_max-dist_crit)/8;
scale_h = [dist_max,dist_max-step_h,dist_max-2*step_h,dist_max-3*step_h,...
    dist_max-4*step_h,dist_max-5*step_h,dist_max-6*step_h,dist_max-7*step_h,dist_max-8*step_h];
h=colorbar('YTickLabel',{scale_h});
set(h,'ylim',[0.1 0.9]);
title = get(h,'Title');
titleString = 'Distance (mm)';
set(title ,'String',titleString,'FontWeight','bold');

% Trapezium:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read STL files
[F_Bone3, V_Bone3] =  STL_ReadFile([dir Bone2],true);
[F_Bone4, V_Bone4] =  STL_ReadFile([dir Bone1],true);

% Calculation of the shortest distance between the two STL selected
[idx3,D3] = knnsearch(V_Bone3,V_Bone4);
[C3,idxD3] = min(D3);

% Plot the two STL selected with an arrow joining the two points
%       calculated as the closest neighbors
[idx4,D4] = knnsearch(V_Bone4,V_Bone3);
[C4,idxD4] = min(D4);

% Calculation of the matrix D which contains the Euclidean distances between
%       each points of each STL file
D3 = pdist2(V_Bone3,V_Bone4);

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
    CA_V_Bone3(i,:) = V_Bone3(idxD3(i),:);
end
for i = 1:length(idxD4)
    CA_V_Bone4(i,:) = V_Bone4(idxD4(i),:);
end

%     Then we assign each point from the first STL with the corresponding
%     one from the second STL
for i = 1:length(CA_V_Bone3)
    DCAt = pdist2(CA_V_Bone3(i,:),CA_V_Bone4);
    [minDCAt(i),idxCAt(i)] = min(DCAt);
end

% Now we can calculate the contact surface area by calculating each
% triangle's area formed by the points located in the contact zone
for i = 1:length(idxD3)
    CA_F_Bone3(i,:) = F_Bone3(idxD3(i),:);
end
for i = 1:length(idxD4)
    CA_F_Bone4(i,:) = F_Bone4(idxD4(i),:);
end

% Define the condition to consider the zone of proximity (red zone)
Def3prox = D3 < 1.5*dist_crit.';

%     Now we calculate the matrix idxD which corresponds to a 2xn matrix.
%     The first row gives all the indexes of the points located in the
%     zone of interest of the V_Bone1 matrix. The second row gives all the
%     indexes of the points located in the zone of interest of the V_Bone2
%     matrix.
[idxrDtprox,idxcDtprox] = find(D3<dist_max);
idxDtprox = [idxrDtprox,idxcDtprox];
idxDtprox(:,any(idxDtprox==0,1)) = []; % Removes all the zero in the idxD Matrix
idxD3prox = unique(idxDtprox(:,1));
idxD4prox = unique(idxDtprox(:,2));

%     Then we isolate the coordinates of each point of each STL file
%     located in the zone of interest in two different matrices.
for i = 1:length(idxD3prox)
    CA_V_Bone3prox(i,:) = V_Bone3(idxD3prox(i),:);
end
for i = 1:length(idxD4prox)
    CA_V_Bone4prox(i,:) = V_Bone4(idxD4prox(i),:);
end

%     Then we assign each point from the first STL with the corresponding
%     one from the second STL
for i = 1:length(CA_V_Bone3prox)
    DCAtprox = pdist2(CA_V_Bone3prox(i,:),CA_V_Bone4prox);
    [minDCAtprox(i),idxCAtprox(i)] = min(DCAtprox);
end

% Now we can calculate the contact surface area by calculating each
% triangle's area formed by the points located in the contact zone
for i = 1:length(idxD3prox)
    CA_F_Bone3prox(i,:) = F_Bone3(idxD3prox(i),:);
end
for i = 1:length(idxD4prox)
    CA_F_Bone4prox(i,:) = F_Bone4(idxD4prox(i),:);
end

%       Each triangle is composed of 3 points called a,b and c
Xa3 = V_Bone3(CA_F_Bone3prox(:,1),1);
Ya3 = V_Bone3(CA_F_Bone3prox(:,1),2);
Za3 = V_Bone3(CA_F_Bone3prox(:,1),3);

Xb3 = V_Bone3(CA_F_Bone3prox(:,2),1);
Yb3 = V_Bone3(CA_F_Bone3prox(:,2),2);
Zb3 = V_Bone3(CA_F_Bone3prox(:,2),3);

Xc3 = V_Bone3(CA_F_Bone3prox(:,3),1);
Yc3 = V_Bone3(CA_F_Bone3prox(:,3),2);
Zc3 = V_Bone3(CA_F_Bone3prox(:,3),3);

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

Articular_AreaTrap = sum(S3(:))

%     Now we can plot those points located in the compression area directly
%     on the STL files in order to visualize the contact area
XB3 = CA_V_Bone3(:,1);
YB3 = CA_V_Bone3(:,2);
ZB3 = CA_V_Bone3(:,3);
XB4 = CA_V_Bone4(:,1);
YB4 = CA_V_Bone4(:,2);
ZB4 = CA_V_Bone4(:,3);

% Plots the STL Bone2 with the color map representing the stress
% distribution in the compression zone
area2 = figure; 
    [obj, li, ax] = GUI_PlotShells(area2, {F_Bone3}, {V_Bone3},...
            {ones(size(V_Bone3,1),1)},[0,0,1]);
box off

for i=1:length(V_Bone3)
    colrstl2(i,:)=[0,0,1];
end

hold on
for i=1:length(idxD3)
    Dist2(i)=((XB3(i)-XB4(idxCAt(i)))^2+(YB3(i)-YB4(idxCAt(i)))^2+(ZB3(i)-ZB4(idxCAt(i)))^2)^(1/2);
       colrstl2(idxD3(i),:)= ...
        [-(Dist2(i)^2)/(dist_max^2)+1,...
        -((Dist2(i)-(dist_max+dist_crit)/2)^2)/(((dist_max-dist_crit)/2)^2)+1,...
        -((Dist2(i)-(dist_max+dist_crit))^2)/((dist_max)^2)+1];
end
hold on
patch('Faces',F_Bone3,'Vertices',V_Bone3, ...
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
step_h = (dist_max-dist_crit)/8;
scale_h = [dist_max,dist_max-step_h,dist_max-2*step_h,dist_max-3*step_h,...
    dist_max-4*step_h,dist_max-5*step_h,dist_max-6*step_h,dist_max-7*step_h,dist_max-8*step_h];
h=colorbar('YTickLabel',{scale_h});
set(h,'ylim',[0.1 0.9]);
title = get(h,'Title');
titleString = 'Distance (mm)';
set(title ,'String',titleString,'FontWeight','bold');

%   Plots the trapezium with only the critical proximity area
areaprox2 = figure; 
    [obj, li, ax] = GUI_PlotShells(areaprox2, {F_Bone3}, {V_Bone3},...
            {ones(size(V_Bone3,1),1)},[0,0,1]);
box off

for i=1:length(V_Bone3)
    colrstlprox2(i,:)=[0,0,1];
end

hold on
for i=1:length(idxD3prox)
       colrstlprox2(idxD3prox(i),:)= [1,0,0];
end
hold on
patch('Faces',F_Bone3,'Vertices',V_Bone3, ...
    'FaceColor','interp', ...
    'FaceVertexCData',colrstlprox2, ...
    'EdgeColor', 'interp', ...
    'EdgeAlpha', 0, ...
    'CDataMapping', 'scaled',...
    'AmbientStrength', 0.4, ...
    'DiffuseStrength', 0.8, ...
    'SpecularStrength', 0.2, ...
    'SpecularColorReflectance', 0.5, ...
    'FaceLighting', 'gouraud')
hold on
step_h = (dist_max-dist_crit)/8;
scale_h = [dist_max,dist_max-step_h,dist_max-2*step_h,dist_max-3*step_h,...
    dist_max-4*step_h,dist_max-5*step_h,dist_max-6*step_h,dist_max-7*step_h,dist_max-8*step_h];
h=colorbar('YTickLabel',{scale_h});
set(h,'ylim',[0.1 0.9]);
title = get(h,'Title');
titleString = 'Distance (mm)';
set(title ,'String',titleString,'FontWeight','bold');