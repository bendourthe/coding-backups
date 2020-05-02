%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                               CalcCentroid                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function: Calculated the distance between centroids of a 3D mesh
%
% Dependencies:      STL_ReadFile.m 
%                    TRI_RemoveInvalidTriangles.m
%                       DeleteUnreferencedElements.m
%                    TRI_RemoveBadlyConnectedTriangles.m
%                       TRI_Edges.m
%                       StackEqualElementIndices.m
%                           RunLengthDecode.m
%                           IncrementalRuns.m
%                       TRI_Normals.m
%                           VectorNorms.m
%                       DeleteUnreferencedElements.m% 
%
% Input: STL's
%
% Output: distance between two centroids
% 
% Created by: Faes Kerkhof
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear & Close
close all;
clear all;

%% Load data

% nrBeads = 8; % set nr of bead you want to analyse
% names = dir('C:\Users\u0078973\Desktop\Static CT beads STD\*.stl');   

FrameNr = 25; % set patient ID nr
PathStr = 'C:\Users\u0078973\Desktop\Static CT beads STD\'; % set patient nr and STL folder
Files = dir(fullfile(PathStr,'*.stl')); 

for k=1:length(Files)
    [F{k}, V{k}, fpath{k}, fpos{k}]=STL_ReadFile(fullfile(PathStr,Files(k).name));
end

% for z=1:length(Files)
%     Vmat(z) = cell2mat(V(:,:));
% end

V1 = cell2mat(V(1,1));
V2 = cell2mat(V(1,2));
V3 = cell2mat(V(1,3));
V4 = cell2mat(V(1,4));
V5 = cell2mat(V(1,5));
V6 = cell2mat(V(1,6));
V7 = cell2mat(V(1,7));
V8 = cell2mat(V(1,8));
% [F1, V1, fpath1, fpos1] =  STL_ReadFile;
% [F2, V2, fpath2, fpos2] =  STL_ReadFile;
% [F3, V3, fpath3, fpos3] =  STL_ReadFile;
% [F4, V4, fpath4, fpos4] =  STL_ReadFile;
% [F5, V5, fpath5, fpos5] =  STL_ReadFile;

% Calculating centroid
for i = 1:length(Files)
    CenV{i} = mean(V{1,i});
end

CenV1 = cell2mat(CenV(1,1));
CenV2 = cell2mat(CenV(1,2));
CenV3 = cell2mat(CenV(1,3));
CenV4 = cell2mat(CenV(1,4));
CenV5 = cell2mat(CenV(1,5));
CenV6 = cell2mat(CenV(1,6));
CenV7 = cell2mat(CenV(1,7));
CenV8 = cell2mat(CenV(1,8));
% CenV1 = mean(V1); % Define the 2 centroids. Centroid 1 at (x1, y1, z1) and centroid 2 at (x2, y2, z2).
% CenV2 = mean(V2); % Make sure they are in the same coordinate system.
% CenV3 = mean(V3);
% CenV4 = mean(V4);
% CenV5 = mean(V5);

%% Plot centroid
hold on 
scatter3(V1(:,1), V1(:,2), V1(:,3),'m')
scatter3(CenV1(1,1), CenV1(1,2), CenV1(1,3),'filled') % visualy check if centroid is correct

scatter3(V2(:,1), V2(:,2), V2(:,3),'c')
scatter3(CenV2(1,1), CenV2(1,2), CenV2(1,3),'filled') % visualy check if centroid is correct

scatter3(V3(:,1), V3(:,2), V3(:,3),'y')
scatter3(CenV3(1,1), CenV3(1,2), CenV3(1,3),'filled')

scatter3(V4(:,1), V4(:,2), V4(:,3),'+b')
scatter3(CenV4(1,1), CenV4(1,2), CenV4(1,3),'filled')

scatter3(V5(:,1), V5(:,2), V5(:,3),'dr')
scatter3(CenV5(1,1), CenV5(1,2), CenV5(1,3),'filled')

scatter3(V6(:,1), V6(:,2), V6(:,3),'y')
scatter3(CenV6(1,1), CenV6(1,2), CenV6(1,3),'filled')

scatter3(V7(:,1), V7(:,2), V7(:,3),'+b')
scatter3(CenV7(1,1), CenV7(1,2), CenV7(1,3),'filled')

scatter3(V8(:,1), V8(:,2), V8(:,3),'dr')
scatter3(CenV8(1,1), CenV8(1,2), CenV8(1,3),'filled')

hold off

%% Calculating distance between centroids  
% n = 5; % number of beads
% for i = 1:n
    
% C(i) = [CenV(i)(1,1), CenV(i)(1,2), CenV(i)(1,3)]; % coordinate of
% centroid 1

% end


C1 = [CenV1(1,1), CenV1(1,2), CenV1(1,3)]; % coordinate of centroid 1
C2 = [CenV2(1,1), CenV2(1,2), CenV2(1,3)]; % coordinate of centroid 2
C3 = [CenV3(1,1), CenV3(1,2), CenV3(1,3)]; % etc.
C4 = [CenV4(1,1), CenV4(1,2), CenV4(1,3)];
C5 = [CenV5(1,1), CenV5(1,2), CenV5(1,3)];
C6 = [CenV6(1,1), CenV6(1,2), CenV6(1,3)]; % etc.
C7 = [CenV7(1,1), CenV7(1,2), CenV7(1,3)];
C8 = [CenV8(1,1), CenV8(1,2), CenV8(1,3)];

xd1 = C2(1,1)- C1(1,1); % distances in x,y and z direction
yd1 = C2(1,2)- C1(1,2);
zd1 = C2(1,3)- C1(1,3);

Distance1 = sqrt(xd1^2 + yd1^2 + zd1^2); % sqrt is not beneficial for fast code

xd2 = C3(1,1)- C1(1,1); 
yd2 = C3(1,2)- C1(1,2);
zd2 = C3(1,3)- C1(1,3);
 
Distance2 = sqrt(xd2^2 + yd2^2 + zd2^2);


xd3 = C4(1,1)- C1(1,1); 
yd3 = C4(1,2)- C1(1,2);
zd3 = C4(1,3)- C1(1,3);
 
Distance3 = sqrt(xd3^2 + yd3^2 + zd3^2);


xd4 = C5(1,1)- C1(1,1); 
yd4 = C5(1,2)- C1(1,2);
zd4 = C5(1,3)- C1(1,3);
 
Distance4 = sqrt(xd4^2 + yd4^2 + zd4^2);

xd5 = C6(1,1)- C1(1,1); 
yd5 = C6(1,2)- C1(1,2);
zd5 = C6(1,3)- C1(1,3);
 
Distance5 = sqrt(xd5^2 + yd5^2 + zd5^2);

xd6 = C7(1,1)- C1(1,1); 
yd6 = C7(1,2)- C1(1,2);
zd6 = C7(1,3)- C1(1,3);
 
Distance6 = sqrt(xd6^2 + yd6^2 + zd6^2);

xd7 = C8(1,1)- C1(1,1); 
yd7 = C8(1,2)- C1(1,2);
zd7 = C8(1,3)- C1(1,3);
 
Distance7 = sqrt(xd7^2 + yd7^2 + zd7^2);


xd8 = C2(1,1)- C3(1,1); 
yd8 = C2(1,2)- C3(1,2);
zd8 = C2(1,3)- C3(1,3);
 
Distance8 = sqrt(xd8^2 + yd8^2 + zd8^2);


xd9 = C2(1,1)- C4(1,1); 
yd9 = C2(1,2)- C4(1,2);
zd9 = C2(1,3)- C4(1,3);
 
Distance9 = sqrt(xd9^2 + yd9^2 + zd9^2);

xd10 = C2(1,1)- C5(1,1); 
yd10 = C2(1,2)- C5(1,2);
zd10 = C2(1,3)- C5(1,3);
 
Distance10 = sqrt(xd10^2 + yd10^2 + zd10^2);

xd11 = C2(1,1)- C6(1,1); 
yd11 = C2(1,2)- C6(1,2);
zd11 = C2(1,3)- C6(1,3);
 
Distance11 = sqrt(xd11^2 + yd11^2 + zd11^2);

xd12 = C2(1,1)- C7(1,1); 
yd12 = C2(1,2)- C7(1,2);
zd12 = C2(1,3)- C7(1,3);
 
Distance12 = sqrt(xd12^2 + yd12^2 + zd12^2);

    
xd13 = C2(1,1)- C8(1,1); 
yd13 = C2(1,2)- C8(1,2);
zd13 = C2(1,3)- C8(1,3);
 
Distance13 = sqrt(xd13^2 + yd13^2 + zd13^2);   
    
xd14 = C3(1,1)- C4(1,1); 
yd14 = C3(1,2)- C4(1,2);
zd14 = C3(1,3)- C4(1,3);
 
Distance14 = sqrt(xd14^2 + yd14^2 + zd14^2);
 
xd15 = C3(1,1)- C5(1,1); 
yd15 = C3(1,2)- C5(1,2);
zd15 = C3(1,3)- C5(1,3);
 
Distance15 = sqrt(xd15^2 + yd15^2 + zd15^2);
 
xd16 = C3(1,1)- C6(1,1); 
yd16 = C3(1,2)- C6(1,2);
zd16 = C3(1,3)- C6(1,3);
 
Distance16 = sqrt(xd16^2 + yd16^2 + zd16^2);
 
xd17 = C3(1,1)- C7(1,1); 
yd17 = C3(1,2)- C7(1,2);
zd17 = C3(1,3)- C7(1,3);
 
Distance17 = sqrt(xd17^2 + yd17^2 + zd17^2);
 
    
xd18 = C3(1,1)- C8(1,1); 
yd18 = C3(1,2)- C8(1,2);
zd18 = C3(1,3)- C8(1,3);
 
Distance18 = sqrt(xd18^2 + yd18^2 + zd18^2);   
    
xd20 = C4(1,1)- C5(1,1); 
yd20 = C4(1,2)- C5(1,2);
zd20 = C4(1,3)- C5(1,3);
 
Distance20 = sqrt(xd20^2 + yd20^2 + zd20^2);
 
xd21 = C4(1,1)- C6(1,1); 
yd21 = C4(1,2)- C6(1,2);
zd21 = C4(1,3)- C6(1,3);
 
Distance21 = sqrt(xd21^2 + yd21^2 + zd21^2);
 
xd22 = C4(1,1)- C7(1,1); 
yd22 = C4(1,2)- C7(1,2);
zd22 = C4(1,3)- C7(1,3);
 
Distance22 = sqrt(xd22^2 + yd22^2 + zd22^2);
 
    
xd23 = C4(1,1)- C8(1,1); 
yd23 = C4(1,2)- C8(1,2);
zd23 = C4(1,3)- C8(1,3);
 
Distance23 = sqrt(xd23^2 + yd23^2 + zd23^2);   
    

xd24 = C5(1,1)- C6(1,1); 
yd24 = C5(1,2)- C6(1,2);
zd24 = C5(1,3)- C6(1,3);
 
Distance24 = sqrt(xd24^2 + yd24^2 + zd24^2);
 
xd25 = C5(1,1)- C7(1,1); 
yd25 = C5(1,2)- C7(1,2);
zd25 = C5(1,3)- C7(1,3);
 
Distance25 = sqrt(xd25^2 + yd25^2 + zd25^2);
 
    
xd26 = C5(1,1)- C8(1,1); 
yd26 = C5(1,2)- C8(1,2);
zd26 = C5(1,3)- C8(1,3);
 
Distance26 = sqrt(xd26^2 + yd26^2 + zd26^2);   
    


xd27 = C6(1,1)- C7(1,1); 
yd27 = C6(1,2)- C7(1,2);
zd27 = C6(1,3)- C7(1,3);
 
Distance27 = sqrt(xd27^2 + yd27^2 + zd27^2);
 
    
xd28 = C6(1,1)- C8(1,1); 
yd28 = C6(1,2)- C8(1,2);
zd28 = C6(1,3)- C8(1,3);
 
Distance28 = sqrt(xd28^2 + yd28^2 + zd28^2);   
    
    
    
AllDistance =   [Distance1...
Distance2...
Distance3...
Distance4...
Distance5...
Distance6...
Distance7...
Distance8...
Distance9...
Distance10...
Distance11...
Distance12...
Distance13...
Distance14...
Distance15...
Distance16...
Distance17...
Distance18...
Distance20...
Distance21...
Distance22...
Distance23...
Distance24...
Distance25...
Distance26...
Distance27...
Distance28];

%% Write all calculated distances to excel file
filename = [PathStr 'FrameNr_' num2str(FrameNr) '.xlsx'];
xlswrite(filename,AllDistance);


