clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Proximity_Model:
% - Calculates the shortest distance between two STL
% - Plots each STL with an arrow joining the two closest neighbours
% - Plots each STL with a colormap indicating the relative distance between
% each articular surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dependencies:        
%               STL_ReadFile.m 
%                   TRI_RemoveInvalidTriangles.m
%                   TRI_RemoveBadlyConnectedTriangles.m
%               GUI_PlotShells
% Input: 
%           Bone1: STL file corresponding to the first bone
%           Bone2: STL file corresponding to second bone
% 
% Output:
%           dist_min: Shortest distance calculated between the two STL
%                   files selected (mm)
%           PA_MC1: Proximity area of the first metacarpal (mm2)
%           PA_Trap: Proximity area of the trapezium (mm2)
%           dist_av: Average distance calculated between the two STL
%                   files selected (mm)
%           Figure 1: Shows the two bones
%           Figure 2: Shows the proximity pattern (full range) of the MC1
%           Figure 3: Shows the proximity area of the MC1 (d < 1.5)
%           Figure 4: Shows the proximity pattern (full range) of the Trap 
%           Figure 5: Shows the proximity area of the Trap (d < 1.5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Selection of the current directory (where the STL files are)
dir = 'C:\Users\u0096636\Documents\Professional\PhD KU Leuven\Data\Study_Microarchitecture_Trapezium\Mimics Projects\CBCT\';

% FILES
Bone1 = '7.7_mc1_red.stl';
Bone2 = '7.7_trap_red.stl';
    
% Setting the maximal distance that the code will consider between the 2
% bones and the distance that we consider critical (underneath which we can
% expect high amount of stress due to high compression)
dist_max = 3;       % in mm
dist_prox =   1.5;    % in mm

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
    CA_V_Bone3prox(i,:) = V_Bone3(idxD3prox(i),:);
end
for i = 1:length(idxD4prox)
    CA_V_Bone4prox(i,:) = V_Bone4(idxD4prox(i),:);
end

%     Then we assign each point from the first STL with the corresponding
%     one from the second STL
for i = 1:length(CA_V_Bone3prox)
    DCAtprox = pdist2(CA_V_Bone3prox(i,:),CA_V_Bone4prox);
    [minDCAtprox(i),idxCAprox2(i)] = min(DCAtprox);
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

PA_Trap = sum(S3(:))

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
    [obj, li, ax] = GUI_PlotShells(area2, {F_Bone3}, {V_Bone3},...
            {ones(size(V_Bone3,1),1)},[0,0,1]);
box off
view(0,-45);
for i=1:length(V_Bone3)
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
step_h = (dist_max-1)/8;
scale_h = [dist_max,dist_max-step_h,dist_max-2*step_h,dist_max-3*step_h,...
    dist_max-4*step_h,dist_max-5*step_h,dist_max-6*step_h,dist_max-7*step_h,dist_max-8*step_h];
h=colorbar('YTickLabel',{scale_h});
set(h,'ylim',[0.1 0.9]);
title = get(h,'Title');
titleString = 'Distance (mm)';
set(title ,'String',titleString,'FontWeight','bold');
set(gcf,'numbertitle','off','name','Trapezium - Proximity pattern (Full range)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%