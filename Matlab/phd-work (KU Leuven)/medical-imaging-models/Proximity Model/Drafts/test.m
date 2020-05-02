% Selection of the current directory (where the STL files are)
dir = 'D:\PhD KU Leuven\Data\Study_Brown_Rig\New Segmentation\Scan4\STL\Priscilla\';

% FILES
Bone1 = 'SCAN4_N_mc1.stl';
Bone2 = 'SCAN4_N_trap.stl';

% Setting the maximal distance that the code will consider between the 2
% bones and the distance that we consider critical (underneath which we can
% expect high amount of stress due to high compression)
dist_max = 3;       % in mm
dist_crit =    1.0063;    % in mm

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

%       Each triangle is composed of 3 points called a,b and c
Xa1 = V_Bone1(CA_F_Bone1(:,1),1);
Ya1 = V_Bone1(CA_F_Bone1(:,1),2);
Za1 = V_Bone1(CA_F_Bone1(:,1),3);

Xb1 = V_Bone1(CA_F_Bone1(:,2),1);
Yb1 = V_Bone1(CA_F_Bone1(:,2),2);
Zb1 = V_Bone1(CA_F_Bone1(:,2),3);

Xc1 = V_Bone1(CA_F_Bone1(:,3),1);
Yc1 = V_Bone1(CA_F_Bone1(:,3),2);
Zc1 = V_Bone1(CA_F_Bone1(:,3),3);

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
for i=1:min(length(idxD1),length(idxD2))
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
step_h = (3-dist_crit)/8;
scale_h = [3,3-step_h,3-2*step_h,3-3*step_h,3-4*step_h,3-5*step_h,3-6*step_h,3-7*step_h,3-8*step_h];
h=colorbar('YTickLabel',{scale_h});
set(h,'ylim',[0.1 0.9]);
title = get(h,'Title');
titleString = 'Distance (mm)';
set(title ,'String',titleString,'FontWeight','bold');