clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contact_Stress_TMC:
% - Calculates the shortest distance between two STL
% - Plots the two STL together
% - Calculates the contact stress values for the cartilage area under
%   compression according to a finite deformation biphasic theory
% - Calculates the maximal and average contact stress for the whole contact
%   area
% - Calculates the projected contact areas (PCA)
% - Plots each STL and shows the stress distribution on the contact areas
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
%           CSav: Average contact stress value (MPa)
%           PCA_MC1: Contact area of the first metacarpal (mm2)
%           PCA_Trap: Contact area of the trapezium (mm2)
%           Figure 1: Shows the 2 bones
%           Figure 2: Shows the stress distribution of the 1st bone
%           Figure 3: Shows the stress distribution of the 2nd bone        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Selection of the current directory (where the STL files are)
dir = 'C:\Users\u0096636\Documents\Professional\PhD KU Leuven\Data\Forward Model\STL\';

% FILES
Bone1 = 'SCAN3_N_mc1.stl';
Bone2 = 'SCAN3_N_trap.stl';

% Cartilage properties
HA = 0.80;            % Average Aggregate Modulus for a human cartilage (MPa)
d0 = 0.90;            % Initial solid content (material constant)
a0 = 0.20;            % Average fluid-to-solid true density ratio
d1 = d0*(a0 + 1)/(a0 + d0); % Material constant
Ttrap = 0.8;          % Average trapezium cartilage thickness (mm)
Tmc1 = 0.7;           % Average first metacarpal cartilage thickness (mm)
Ttot = Ttrap + Tmc1;  % Total joint cartilage thickness (mm)

% MC1:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read STL files
[F_Bone1, V_Bone1] =  STL_ReadFile([dir Bone1],true);
[F_Bone2, V_Bone2] =  STL_ReadFile([dir Bone2],true);

% Calculation of the shortest distance between the two STL selected
[idx1,D1] = knnsearch(V_Bone1,V_Bone2);
[C1,idxD1] = min(D1);
dmin = C1

% Calculation of the matrix D which contains the Euclidean distances between
%       each points of each STL file
D = pdist2(V_Bone1,V_Bone2);

% Calculation of the contact stresses values corresponding to the area where
% the cartillage is under compression (where the deformation is positive)
Def_mc1 = (D(D < Ttot).' - Ttot)/2/Tmc1 + 1;
for i=1:length(Def_mc1)
    if Def_mc1(i) >1
        Def_mc1(i) = 1;
    end
    CS_mc1(i) = 1/4*HA*(1 + d1*(Def_mc1(i) - 1))*1/Def_mc1(i)*(Def_mc1(i).^2 - 1/(Def_mc1(i).^2));
end
Def_trap = (D(D < Ttot).' - Ttot)/2/Ttrap + 1;
for i=1:length(Def_trap)
    if Def_trap(i) >1
        Def_trap(i) = 1;
    end
    CS_trap(i) = 1/4*HA*(1 + d1*(Def_trap(i) - 1))*1/Def_trap(i)*(Def_trap(i).^2 - 1/(Def_trap(i).^2));
end
for i=1:length(Def_mc1)
    CS(i) = CS_mc1(i) + CS_trap(i);
end
CS(CS==0) = []; % Removes all the zero in the CS Matrix
CSmax = max(abs(CS))
CSav = mean(abs(CS))

% Now that we have all the contact stresses values in an 1xn matrix, we
% need to link those values to the corresponding points of the two STL

%     We first calculate the matrix idxD which corresponds to a 2xn matrix.
%     The first row gives all the indexes of the points located in the
%     compression area in the V_Bone1 matrix. The second row gives all the
%     indexes of the points located in the compression area in the V_Bone2
%     matrix.
[idxrD,idxcD] = find(D<Ttot);
idxD = [idxrD,idxcD];
idxD(:,any(idxD==0,1)) = []; % Removes all the zero in the idxD Matrix
idxD1 = unique(idxD(:,1));
idxD2 = unique(idxD(:,2));

%     Then we isolate the coordinates of each point of each STL file
%     located in the compression area in two different matrices.
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

        % We can now test to which triangle these vertices correspond in
        % the F matrix and save their indexes
        test_F1 = ismember(F_Bone1,idxD1);
        idx_F1 = find(test_F1(:,1)==1 & test_F1(:,2)==1 & test_F1(:,3)==1);
            % This idx_F_trap matrix is a Nx1 matrix, where N is the amount
            % of triangles located in the deformation zone, and each value
            % in this matrix corresponds to the index of a line in the F
            % matrix that defines one of these triangles
        for i=1:length(idx_F1)
            CA_F_Bone1(i,:) = F_Bone1(idx_F1(i),:);
        end
        % We can now test to which triangle these vertices correspond in
        % the F matrix and save their indexes
        test_F2 = ismember(F_Bone2,idxD2);
        idx_F2 = find(test_F2(:,1)==1 & test_F2(:,2)==1 & test_F2(:,3)==1);
            % This idx_F_trap matrix is a Nx1 matrix, where N is the amount
            % of triangles located in the deformation zone, and each value
            % in this matrix corresponds to the index of a line in the F
            % matrix that defines one of these triangles
        for i=1:length(idx_F2)
            CA_F_Bone2(i,:) = F_Bone2(idx_F2(i),:);
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

PCA_MC1 = sum(S1(:))

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
view(200,45)

for i=1:length(V_Bone1)
    colrstl1(i,:)=[0,0,1];
end

hold on
for i=1:min(length(idxD1),length(idxD2))
    dist_m1(i) = ((XB1(i)-XB2(idxCA(i)))^2+(YB1(i)-YB2(idxCA(i)))^2+(ZB1(i)-ZB2(idxCA(i)))^2)^(1/2);
    Def_mc1(i) = (dist_m1(i) - Ttot)/2/Tmc1+1;
    if Def_mc1(i) >1
        Def_mc1(i) = 1;
    end
    Def_trap(i) = (dist_m1(i) - Ttot)/2/Ttrap+1;
    if Def_trap(i) >1
        Def_trap(i) = 1;
    end
    stress1_mc1(i) = 1/4*HA*(1 + d1*(Def_mc1(i) - 1))*1/Def_mc1(i)*(Def_mc1(i).^2 - 1/(Def_mc1(i).^2));
    stress1_trap(i) = 1/4*HA*(1 + d1*(Def_trap(i) - 1))*1/Def_trap(i)*(Def_trap(i).^2 - 1/(Def_trap(i).^2));
    stress1(i) = stress1_mc1(i) + stress1_trap(i);
    
    if abs(stress1(i)) >= HA
       colrstl1(idxD1(i),:) = [1,0,0];
    elseif abs(stress1(i)) < HA
        colrstl1(idxD1(i),:)= ...
        [-((abs(stress1(i))-HA)^2)/((HA)^2)+1,...
        -((abs(stress1(i))-HA/2)^2)/((-HA/2)^2)+1,...
        -(abs(stress1(i))^2)/((HA)^2)+1];    
    end
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
h=colorbar('YTickLabel',{'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8 = HA or more'});
set(h,'ylim',[0.1 0.9]);
title = get(h,'Title');
titleString = 'MPa';
set(title ,'String',titleString,'FontWeight','bold');
set(gcf,'numbertitle','off','name','Stress distribution in the MC1');

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

% Calculation of the contact stresses values corresponding to the area where
% the cartillage is under compression (where the deformation is positive)
Def3_mc1 = (D3(D3 < Ttot).' - Ttot)/2/Tmc1 + 1;
for i=1:length(Def3_mc1)
    if Def3_mc1(i) >1
        Def3_mc1(i) = 1;
    end
    CS3_mc1(i) = 1/4*HA*(1 + d1*(Def3_mc1(i) - 1))*1/Def3_mc1(i)*(Def3_mc1(i).^2 - 1/(Def3_mc1(i).^2));
end
Def3_trap = (D3(D3 < Ttot).' - Ttot)/2/Ttrap + 1;
for i=1:length(Def3_trap)
    if Def3_trap(i) >1
        Def3_trap(i) = 1;
    end
    CS3_trap(i) = 1/4*HA*(1 + d1*(Def3_trap(i) - 1))*1/Def3_trap(i)*(Def3_trap(i).^2 - 1/(Def3_trap(i).^2));
end
for i=1:length(Def3_mc1)
    CS3(i) = CS3_mc1(i) + CS3_trap(i);
end
CS3(CS3==0) = []; % Removes all the zero in the CS Matrix

% Now that we have all the contact stresses values in an 1xn matrix, we
% need to link those values to the corresponding points of the two STL

%     We first calculate the matrix idxD which corresponds to a 2xn matrix.
%     The first row gives all the indexes of the points located in the
%     compression area in the V_Bone1 matrix. The second row gives all the
%     indexes of the points located in the compression area in the V_Bone2
%     matrix.
[idxrDt,idxcDt] = find(D3<Ttot);
idxDt = [idxrDt,idxcDt];
idxDt(:,any(idxDt==0,1)) = []; % Removes all the zero in the idxD Matrix
idxD3 = unique(idxDt(:,1));
idxD4 = unique(idxDt(:,2));

%     Then we isolate the coordinates of each point of each STL file
%     located in the compression area in two different matrices.
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

% We can now test to which triangle these vertices correspond in
        % the F matrix and save their indexes
        test_F3 = ismember(F_Bone3,idxD3);
        idx_F3 = find(test_F3(:,1)==1 & test_F3(:,2)==1 & test_F3(:,3)==1);
            % This idx_F_trap matrix is a Nx1 matrix, where N is the amount
            % of triangles located in the deformation zone, and each value
            % in this matrix corresponds to the index of a line in the F
            % matrix that defines one of these triangles
        for i=1:length(idx_F3)
            CA_F_Bone3(i,:) = F_Bone3(idx_F3(i),:);
        end
        % We can now test to which triangle these vertices correspond in
        % the F matrix and save their indexes
        test_F4 = ismember(F_Bone4,idxD4);
        idx_F4 = find(test_F4(:,1)==1 & test_F4(:,2)==1 & test_F4(:,3)==1);
            % This idx_F_trap matrix is a Nx1 matrix, where N is the amount
            % of triangles located in the deformation zone, and each value
            % in this matrix corresponds to the index of a line in the F
            % matrix that defines one of these triangles
        for i=1:length(idx_F4)
            CA_F_Bone4(i,:) = F_Bone4(idx_F4(i),:);
        end

%       Each triangle is composed of 3 points called a,b and c
Xa3 = V_Bone3(CA_F_Bone3(:,1),1);
Ya3 = V_Bone3(CA_F_Bone3(:,1),2);
Za3 = V_Bone3(CA_F_Bone3(:,1),3);

Xb3 = V_Bone3(CA_F_Bone3(:,2),1);
Yb3 = V_Bone3(CA_F_Bone3(:,2),2);
Zb3 = V_Bone3(CA_F_Bone3(:,2),3);

Xc3 = V_Bone3(CA_F_Bone3(:,3),1);
Yc3 = V_Bone3(CA_F_Bone3(:,3),2);
Zc3 = V_Bone3(CA_F_Bone3(:,3),3);

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

PCA_Trap = sum(S3(:))

%     Now we can plot those points located in the compression area directly
%     on the STL files in order to visualize the contact area
XB3 = CA_V_Bone3(:,1);
YB3 = CA_V_Bone3(:,2);
ZB3 = CA_V_Bone3(:,3);
XB4 = CA_V_Bone4(:,1);
YB4 = CA_V_Bone4(:,2);
ZB4 = CA_V_Bone4(:,3);

% Plots the STL Bone1 with the color map representing the stress
% distribution in the compression zone
area2 = figure; 
    [obj, li, ax] = GUI_PlotShells(area2, {F_Bone3}, {V_Bone3},...
            {ones(size(V_Bone3,1),1)},[0,0,1]);
box off
view(0,-60)

for i=1:length(V_Bone3)
    colrstl2(i,:)=[0,0,1];
end

hold on
for i=1:min(length(idxD3),length(idxD4))
    dist_m2(i) = ((XB3(i)-XB4(idxCAt(i)))^2+(YB3(i)-YB4(idxCAt(i)))^2+(ZB3(i)-ZB4(idxCAt(i)))^2)^(1/2);
    Def3_mc1(i) = (dist_m2(i) - Ttot)/2/Tmc1+1;
    Def3_trap(i) = (dist_m2(i) - Ttot)/2/Ttrap+1;
    stress2_mc1(i) = 1/4*HA*(1 + d1*(Def3_mc1(i) - 1))*1/Def3_mc1(i)*(Def3_mc1(i).^2 - 1/(Def3_mc1(i).^2));
    stress2_trap(i) = 1/4*HA*(1 + d1*(Def3_trap(i) - 1))*1/Def3_trap(i)*(Def3_trap(i).^2 - 1/(Def3_trap(i).^2));
    stress2(i) = stress2_mc1(i) + stress2_trap(i);
    
    if abs(stress2(i)) >= HA
       colrstl2(idxD3(i),:) = [1,0,0];
    elseif abs(stress2(i)) < HA
        colrstl2(idxD3(i),:)= ...
        [-((abs(stress2(i))-HA)^2)/((HA)^2)+1,...
        -((abs(stress2(i))-HA/2)^2)/((-HA/2)^2)+1,...
        -(abs(stress2(i))^2)/((HA)^2)+1];    
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
h=colorbar('YTickLabel',{'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8 = HA or more'});
set(h,'ylim',[0.1 0.9]);
title = get(h,'Title');
titleString = 'MPa';
set(title ,'String',titleString,'FontWeight','bold');
set(gcf,'numbertitle','off','name','Stress distribution in the trapezium');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%