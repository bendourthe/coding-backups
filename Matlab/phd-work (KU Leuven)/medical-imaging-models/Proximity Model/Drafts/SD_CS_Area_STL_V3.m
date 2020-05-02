clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SD_CS_Area_STL_V2:
% - Calculates the shortest distance between two STL
% - Calculates the macimal contact stress value
% - Plots each STL with an arrow joining the two closest neighbours
% - Calculates the contact stress values for the cartilage area under
%   compression according to a finite deformation biphasic theory
% - Calculates the average contact stress for the whole contact area
% - Calculates the contact areas
% - Plots each STL and shows the stress distribution on the contact areas
% - Quantify the stress distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dependencies:        
%               STL_ReadFile.m 
%                   TRI_RemoveInvalidTriangles.m
%                   TRI_RemoveBadlyConnectedTriangles.m
%               vectarrow.m
% Input: 
%           Bone1: STL file corresponding to the first bone
%           Bone2: STL file corresponding to second bone
% 
% Output:
%           Shortest_Distance: Shortest distance calculated between the two
%                               STL files selected (mm)
%           CSmax: Maximal contact stress value (MPa)
%           Contact_AreaMC1: Contact area of the first metacarpal (mm2)
%           Contact_AreaTrap: Contact area of the trapezium (mm2)
%           Figure 1: Shows the 2 bones with an arrow linking the 2 closest
%                     neighbours
%           Figure 2: Shows the stress distribution of the 1st bone
%           Figure 3: Shows the stress distribution of the 2nd bone
%           Figure 4: Quantify the stress distribution         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Selection of the current directory (where the STL files are)
dir = 'D:\PhD KU Leuven\Data\Study_Brown_Rig\New Segmentation\Scan17\STL\';

% FILES
Bone1 = 'N_mc1_s_0.4_1.stl';
Bone2 = 'N_trap_s_0.4_1.stl';

% Cartilage properties
HA = 0.80;            % Average Aggregate Modulus for a human cartilage (MPa)
d0 = 0.90;            % Initial solid content (material constant)
a0 = 0.20;            % Average fluid-to-solid true density ratio
d1 = d0*(a0 + 1)/(a0 + d0); % Material constant
Ttrap = 0.8;          % Average trapezium cartilage thickness (mm)
Tmc1 = 0.7;           % Average first metacarpal cartilage thickness (mm)
Ttot = Ttrap + Tmc1;  % Total joint cartilage thickness (mm)

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
hold on
vectarrow(P1,P2);

% Calculation of the matrix D which contains the Euclidean distances between
%       each points of each STL file
D = pdist2(V_Bone1,V_Bone2);

% Calculation of the contact stresses values corresponding to the area where
% the cartillage is under compression (where the deformation is positive)
Def = (Ttot - D(D < Ttot).') + 1;
for i=1:length(Def)
    CS(i) = 1/4*HA*(1 + d1*(Def(i) - 1))*1/Def(i)*(Def(i).^2 - 1/(Def(i).^2));
end
CS(CS==0) = []; % Removes all the zero in the CS Matrix
CSmax = max(CS)
CSav = mean(CS)

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

Xa2 = V_Bone2(CA_F_Bone2(:,1),1);
Ya2 = V_Bone2(CA_F_Bone2(:,1),2);
Za2 = V_Bone2(CA_F_Bone2(:,1),3);

Xb2 = V_Bone2(CA_F_Bone2(:,2),1);
Yb2 = V_Bone2(CA_F_Bone2(:,2),2);
Zb2 = V_Bone2(CA_F_Bone2(:,2),3);

Xc2 = V_Bone2(CA_F_Bone2(:,3),1);
Yc2 = V_Bone2(CA_F_Bone2(:,3),2);
Zc2 = V_Bone2(CA_F_Bone2(:,3),3);

%       Now we need to calculate the coordinates of the 2 vectors ab and ac
for i=1:length(Xa1)
    Xab1(i) = Xb1(i)-Xa1(i);
    Yab1(i) = Yb1(i)-Ya1(i);
    Zab1(i) = Zb1(i)-Za1(i);
    Xac1(i) = Xc1(i)-Xa1(i);
    Yac1(i) = Yc1(i)-Ya1(i);
    Zac1(i) = Zc1(i)-Za1(i);
end

for i=1:length(Xa2)
    Xab2(i) = Xb2(i)-Xa2(i);
    Yab2(i) = Yb2(i)-Ya2(i);
    Zab2(i) = Zb2(i)-Za2(i);
    Xac2(i) = Xc2(i)-Xa2(i);
    Yac2(i) = Yc2(i)-Ya2(i);
    Zac2(i) = Zc2(i)-Za2(i);
end

%       Now we need to calculate the surface area of each triangle and sum
%       all of them to obtain to total contact area
for i = 1:length(Xab1)
    S1(i) = (1/2)*((Yab1(i)*Zac1(i)-Zab1(i)*Yac1(i))^2+(Zab1(i)*Xac1(i)-Xab1(i)*Zac1(i))^2+(Xab1(i)*Yac1(i)-Yab1(i)*Xac1(i))^2)^(1/2);
end

for i = 1:length(Xab2)
    S2(i) = (1/2)*((Yab2(i)*Zac2(i)-Zab2(i)*Yac2(i))^2+(Zab2(i)*Xac2(i)-Xab2(i)*Zac2(i))^2+(Xab2(i)*Yac2(i)-Yab2(i)*Xac2(i))^2)^(1/2);
end

Contact_AreaMC1 = sum(S1(:))
Contact_AreaTrap = sum(S2(:))

%     Now we can plot those points located in the compression area directly
%     on the STL files in order to visualize the contact area
XB1 = CA_V_Bone1(:,1);
YB1 = CA_V_Bone1(:,2);
ZB1 = CA_V_Bone1(:,3);
XB2 = CA_V_Bone2(:,1);
YB2 = CA_V_Bone2(:,2);
ZB2 = CA_V_Bone2(:,3);

mblue1=[0 0 0];
mcyan1=[0 0 0];
mgreen1=[0 0 0];
myellow1=[0 0 0];
mred1=[0 0 0];

for j=1:length(F_Bone1)
   colr(j,:)=[j/length(F_Bone1),0,0]; 
end

area1 = figure; % Plots the STL Bone1 and the corresponding points
    [obj, li, ax] = GUI_PlotShells(area1, {F_Bone1}, {V_Bone1},...
            {ones(size(V_Bone1,1),1)});

hold on
CSmax1 = plot3(P1(1),P1(2),P1(3),'Marker','o','MarkerSize',10,'Color','k','MarkerFaceColor','k');

hold on
for i=1:length(XB1)
    Def(i)=(Ttot-((XB1(i)-XB2(idxCA(i)))^2+(YB1(i)-YB2(idxCA(i)))^2+(ZB1(i)-ZB2(idxCA(i)))^2)^(1/2))+1;
    if (1/4*HA*(1 + d1*(Def(i) - 1))*1/Def(i)*(Def(i).^2 - 1/Def(i).^2)) < 0.25*HA
        mblue1(i)=XB1(i);
        hold on
        blue1 = scatter3(XB1(i),YB1(i),ZB1(i),'fill','b');
    end
    if (1/4*HA*(1 + d1*(Def(i) - 1))*1/Def(i)*(Def(i).^2 - 1/Def(i).^2)) >= 0.25*HA && (1/4*HA*(1 + d1*(Def(i) - 1))*1/Def(i)*(Def(i).^2 - 1/Def(i).^2)) < 0.50*HA
        mcyan1(i)=XB1(i);
        hold on
        cyan1 = scatter3(XB1(i),YB1(i),ZB1(i),'fill','c');
    end
    if (1/4*HA*(1 + d1*(Def(i) - 1))*1/Def(i)*(Def(i).^2 - 1/Def(i).^2)) >= 0.50*HA && (1/4*HA*(1 + d1*(Def(i) - 1))*1/Def(i)*(Def(i).^2 - 1/Def(i).^2)) < 0.75*HA
        mgreen1(i)=XB1(i);
        hold on
        green1 = scatter3(XB1(i),YB1(i),ZB1(i),'fill','g');
    end
    if (1/4*HA*(1 + d1*(Def(i) - 1))*1/Def(i)*(Def(i).^2 - 1/Def(i).^2)) >= 0.75*HA && (1/4*HA*(1 + d1*(Def(i) - 1))*1/Def(i)*(Def(i).^2 - 1/Def(i).^2)) < HA
        myellow1(i)=XB1(i);
        hold on
        yellow1 = scatter3(XB1(i),YB1(i),ZB1(i),'fill','y');
    end
    if (1/4*HA*(1 + d1*(Def(i) - 1))*1/Def(i)*(Def(i).^2 - 1/Def(i).^2)) >= HA
        mred1(i)=XB1(i);
        hold on
        red1 = scatter3(XB1(i),YB1(i),ZB1(i),'fill','r');
    end
end
hold on
blue1=plot3(0,0,0,'Marker','o','Color','b','MarkerFaceColor','b');
cyan1=plot3(0,0,0,'Marker','o','Color','c','MarkerFaceColor','c');
green1=plot3(0,0,0,'Marker','o','Color','g','MarkerFaceColor','g');
yellow1=plot3(0,0,0,'Marker','o','Color','y','MarkerFaceColor','y');
red1=plot3(0,0,0,'Marker','o','Color','r','MarkerFaceColor','r');
legend([CSmax1,red1,yellow1,green1,cyan1,blue1],['CSmax = ',num2str(CSmax),' MPa'],['> HA = ',num2str(HA),' MPa'],['',num2str(HA*0.75),' < CS < ',num2str(HA),' MPa'],['',num2str(HA*0.5),' < CS < ',num2str(HA*0.75),' MPa'],['',num2str(HA*0.25),' < CS < ',num2str(HA*0.5),' MPa'],['0 < CS < ',num2str(HA*0.25),' MPa'],'Location','NorthEastOutside')

area2 = figure; % Plots the STL Bone2 and the corresponding points
    [obj, li, ax] = GUI_PlotShells(area2, {F_Bone2}, {V_Bone2},...
            {ones(size(V_Bone2,1),1)});

hold on
CSmax2 = plot3(P2(1),P2(2),P2(3),'Marker','o','MarkerSize',10,'Color','k','MarkerFaceColor','k');

hold on
for i=1:length(XB1)
    Def(i)=(Ttot-((XB1(i)-XB2(idxCA(i)))^2+(YB1(i)-YB2(idxCA(i)))^2+(ZB1(i)-ZB2(idxCA(i)))^2)^(1/2))+1;
    if (1/4*HA*(1 + d1*(Def(i) - 1))*1/Def(i)*(Def(i).^2 - 1/Def(i).^2)) < 0.20*max(CS)
        mblue2(i)=XB2(idxCA(i));
        hold on
        blue2 = scatter3(XB2(idxCA(i)),YB2(idxCA(i)),ZB2(idxCA(i)),'fill','b');
    end
    if (1/4*HA*(1 + d1*(Def(i) - 1))*1/Def(i)*(Def(i).^2 - 1/Def(i).^2)) >= 0.25*HA && (1/4*HA*(1 + d1*(Def(i) - 1))*1/Def(i)*(Def(i).^2 - 1/Def(i).^2)) < 0.50*HA
        mcyan2(i)=XB2(idxCA(i));
        hold on
        cyan2 = scatter3(XB2(idxCA(i)),YB2(idxCA(i)),ZB2(idxCA(i)),'fill','c');
    end
    if (1/4*HA*(1 + d1*(Def(i) - 1))*1/Def(i)*(Def(i).^2 - 1/Def(i).^2)) >= 0.50*HA && (1/4*HA*(1 + d1*(Def(i) - 1))*1/Def(i)*(Def(i).^2 - 1/Def(i).^2)) < 0.75*HA
        mgreen2(i)=XB2(idxCA(i));
        hold on
        green2 = scatter3(XB2(idxCA(i)),YB2(idxCA(i)),ZB2(idxCA(i)),'fill','g');
    end
    if (1/4*HA*(1 + d1*(Def(i) - 1))*1/Def(i)*(Def(i).^2 - 1/Def(i).^2)) >= 0.75*HA && (1/4*HA*(1 + d1*(Def(i) - 1))*1/Def(i)*(Def(i).^2 - 1/Def(i).^2)) < HA
        myellow2(i)=XB2(idxCA(i));
        hold on
        yellow2 = scatter3(XB2(idxCA(i)),YB2(idxCA(i)),ZB2(idxCA(i)),'fill','y');
    end
    if (1/4*HA*(1 + d1*(Def(i) - 1))*1/Def(i)*(Def(i).^2 - 1/Def(i).^2)) >= HA
        mred2(i)=XB2(idxCA(i));
        hold on
        red2 = scatter3(XB2(idxCA(i)),YB2(idxCA(i)),ZB2(idxCA(i)),'fill','r');
    end
end
blue2=plot3(0,0,0,'Marker','o','Color','b','MarkerFaceColor','b');
cyan2=plot3(0,0,0,'Marker','o','Color','c','MarkerFaceColor','c');
green2=plot3(0,0,0,'Marker','o','Color','g','MarkerFaceColor','g');
yellow2=plot3(0,0,0,'Marker','o','Color','y','MarkerFaceColor','y');
red2=plot3(0,0,0,'Marker','o','Color','r','MarkerFaceColor','r');
legend([CSmax2,red2,yellow2,green2,cyan2,blue2],['CSmax = ',num2str(CSmax),' MPa'],['> HA = ',num2str(HA),' MPa'],['',num2str(HA*0.75),' < CS < ',num2str(HA),' MPa'],['',num2str(HA*0.5),' < CS < ',num2str(HA*0.75),' MPa'],['',num2str(HA*0.25),' < CS < ',num2str(HA*0.5),' MPa'],['0 < CS < ',num2str(HA*0.25),' MPa'],'Location','NorthEastOutside')

% Calculation and representation of the percentage of the contact area that
% represents each color

%   MC1
mblue1(mblue1==0) = [];
mcyan1(mcyan1==0) = [];
mgreen1(mgreen1==0) = [];
myellow1(myellow1==0) = [];
mred1(mred1==0) = [];

pblue1 = length(mblue1)*100/length(XB1);
pcyan1 = length(mcyan1)*100/length(XB1);
pgreen1 = length(mgreen1)*100/length(XB1);
pyellow1 = length(myellow1)*100/length(XB1);
pred1 = length(mred1)*100/length(XB1);

barplot1 = figure;

Pmc = [pblue1,pcyan1,pgreen1,pyellow1,pred1];
Pmc(2,:) = 0; % add a row of zero to allow the stacked bar plot
sb1 = bar(Pmc,'stacked');
ax = axis;ax(2) = 1.6; % change the x axis limit
axis(ax)
Col=findobj(gca,'type','patch');
myCol = [1 0 0
    1 1 0
    0 1 0
    0 1 1
    0 0 1];
for n = 1 : length(Col) 
    set(Col(n),'facecolor',myCol(n,:));
end
colormap(myCol)
ylabel('% of the contact area')
title('Repartition of the contact stress over the contact area')
legend(sb1,['',num2str(pred1,2),' %'],['',num2str(pyellow1,2),' %'],['',num2str(pgreen1,2),' %'],['',num2str(pcyan1,2),' %'],['',num2str(pblue1,2),' %'],'Location','NorthEastOutside')

% There is no need to run the same calculation for the trapezium, since
% each point in the contact area of MC1 is associated to one of the
% trapezium, hence, the total number of points is similar for each
% bone and the percentage obtained for each should be the same (which is
% different from the calculation of the contact area, which doesn't only
% depends on the number of points, but also on the shape of the surface)
