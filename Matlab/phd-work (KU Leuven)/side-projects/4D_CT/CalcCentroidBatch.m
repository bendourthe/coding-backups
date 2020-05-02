clc%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                               CalcCentroidBatch                              %
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
% Input: folder with 5 STL's
%
% Output: distance between all centroids
%
% Created by:   Maarten Afschrift & Faes Kerkhof
%              
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear & Close
close all;
clear all;
clc

% add the path
startpath='C:\Users\u0078973\Desktop\4D 24-7-2013\';
% addpath(fullfile(startpath,'Iterative closest point')); % Why this line?

% locate the different folders
dir_all=dir(fullfile(startpath,'STL Callibration Beads'));
color_sel={'m','c','y','+b','dr','k','g'}; % setting the collors for the beads. 

% counter for output
count=1;
% loop over folders
for i=1:length(dir_all)
    
    if dir_all(i).isdir  && ~strcmp(dir_all(i).name,'.') &&    ~strcmp(dir_all(i).name,'..') % excluding the cells from the array that contain '.' or '...'
        % open the files
        files_dir=dir(fullfile(startpath,'STL Callibration Beads',num2str(dir_all(i).name),'*.stl')); % Uses the info from the cell 'name' in the struct 'dir_all' to name 'file_dir"
        % loop over files
        n_files=length(files_dir); % nr of bead in the dir files_dir
        figure()
        hold on
        
        % pre allocate centroid_pos
        centroid_pos=zeros(n_files,3);
        
        for j=1:n_files
            % load the file
            file_path=fullfile(startpath,'STL Callibration Beads',num2str(dir_all(i).name),files_dir(j).name);
            [data(j).F1, data(j).V1, data(j).fpath1, data(j).fpos1] =  STL_ReadFile(file_path); % places the output from STL_ReadFile from each bead ('j').
            
            % scatter plot
            scatter3(data(j).V1(:,1), data(j).V1(:,2), data(j).V1(:,3),color_sel{j})
            scatter3(mean(data(j).V1(:,1)),mean(data(j).V1(:,2)), mean(data(j).V1(:,3)),'filled') % visualy check if centroid is correct
            
            %
            centroid_pos(j,:)=mean(data(j).V1);
            
        end
        hold off
        key=input('do you want to analyse these files: N or Y','s'); % check visualy if the centroids and spheres make sence.
        
        if strcmp(key,'Y') || strcmp(key,'y')
            % calculate distance
            %            output(count).distance=calculate_distance_centroids(data);
            % pre allocate data_3d
            data_3d=zeros(n_files,3,n_files);
            
            for j=1:n_files
                data_3d(:,:,j)=[centroid_pos(j,1)-centroid_pos(:,1) centroid_pos(j,2)-centroid_pos(:,2) centroid_pos(j,3)-centroid_pos(:,3)];
            end
            data_3d_dist=sqrt(sum(data_3d.^2,2));            
            for j=1:n_files
                data(j).distance=data_3d_dist(:,j);
            end
            
            
            
            % counter
            count=count+1;
            
        else
            
        end
        
        clear data n_files
        
    end
end




% 
% 
% 
% % Load data
% [F1, V1, fpath1, fpos1] =  STL_ReadFile;
% [F2, V2, fpath2, fpos2] =  STL_ReadFile;
% [F3, V3, fpath3, fpos3] =  STL_ReadFile;
% [F4, V4, fpath4, fpos4] =  STL_ReadFile;
% [F5, V5, fpath5, fpos5] =  STL_ReadFile;
% 
% % Calculating centroid
% CenV1 = mean(V1); % Define the 2 centroids. Centroid 1 at (x1, y1, z1) and centroid 2 at (x2, y2, z2).
% CenV2 = mean(V2); % Make sure they are in the same coordinate system.
% CenV3 = mean(V3);
% CenV4 = mean(V4);
% CenV5 = mean(V5);
% 
% % Plot centroid
% hold on
% scatter3(V1(:,1), V1(:,2), V1(:,3),'m')
% scatter3(CenV1(1,1), CenV1(1,2), CenV1(1,3),'filled') % visualy check if centroid is correct
% 
% scatter3(V2(:,1), V2(:,2), V2(:,3),'c')
% scatter3(CenV2(1,1), CenV2(1,2), CenV2(1,3),'filled') % visualy check if centroid is correct
% 
% scatter3(V3(:,1), V3(:,2), V3(:,3),'y')
% scatter3(CenV3(1,1), CenV3(1,2), CenV3(1,3),'filled')
% 
% scatter3(V4(:,1), V4(:,2), V4(:,3),'+b')
% scatter3(CenV4(1,1), CenV4(1,2), CenV4(1,3),'filled')
% 
% scatter3(V5(:,1), V5(:,2), V5(:,3),'dr')
% scatter3(CenV5(1,1), CenV5(1,2), CenV5(1,3),'filled')
% 
% hold off
% 
% % Calculating distance between centroids
% % n = 5; % number of beads
% % for i = 1:n
% 
% % C(i) = [CenV(i)(1,1), CenV(i)(1,2), CenV(i)(1,3)]; % coordinate of
% % centroid 1
% 
% % end
% 
% C1 = [CenV1(1,1), CenV1(1,2), CenV1(1,3)]; % coordinate of centroid 1
% C2 = [CenV2(1,1), CenV2(1,2), CenV2(1,3)]; % coordinate of centroid 2
% C3 = [CenV3(1,1), CenV3(1,2), CenV3(1,3)]; % etc.
% C4 = [CenV4(1,1), CenV4(1,2), CenV4(1,3)];
% C5 = [CenV5(1,1), CenV5(1,2), CenV5(1,3)];
% 
% xd1 = C2(1,1)- C1(1,1); % distances in x,y and z direction
% yd1 = C2(1,2)- C1(1,2);
% zd1 = C2(1,3)- C1(1,3);
% 
% Distance1 = sqrt(xd1^2 + yd1^2 + zd1^2); % sqrt is not beneficial for fast code
% 
% xd2 = C3(1,1)- C1(1,1);
% yd2 = C3(1,2)- C1(1,2);
% zd2 = C3(1,3)- C1(1,3);
% 
% Distance2 = sqrt(xd2^2 + yd2^2 + zd2^2);
% 
% 
% xd3 = C4(1,1)- C1(1,1);
% yd3 = C4(1,2)- C1(1,2);
% zd3 = C4(1,3)- C1(1,3);
% 
% Distance3 = sqrt(xd3^2 + yd3^2 + zd3^2);
% 
% 
% xd4 = C5(1,1)- C1(1,1);
% yd4 = C5(1,2)- C1(1,2);
% zd4 = C5(1,3)- C1(1,3);
% 
% Distance4 = sqrt(xd4^2 + yd4^2 + zd4^2);
% 
% xd5 = C2(1,1)- C3(1,1);
% yd5 = C2(1,2)- C3(1,2);
% zd5 = C2(1,3)- C3(1,3);
% 
% Distance5 = sqrt(xd5^2 + yd5^2 + zd5^2);
% 
% 
% xd6 = C2(1,1)- C4(1,1);
% yd6 = C2(1,2)- C4(1,2);
% zd6 = C2(1,3)- C4(1,3);
% 
% Distance6 = sqrt(xd6^2 + yd6^2 + zd6^2);
% 
% 
% xd7 = C2(1,1)- C5(1,1);
% yd7 = C2(1,2)- C5(1,2);
% zd7 = C2(1,3)- C5(1,3);
% 
% Distance7 = sqrt(xd7^2 + yd7^2 + zd7^2);
% 
% 
% xd8 = C3(1,1)- C4(1,1);
% yd8 = C3(1,2)- C4(1,2);
% zd8 = C3(1,3)- C4(1,3);
% 
% Distance8 = sqrt(xd8^2 + yd8^2 + zd8^2);
% 
% 
% xd9 = C3(1,1)- C5(1,1);
% yd9 = C3(1,2)- C5(1,2);
% zd9 = C3(1,3)- C5(1,3);
% 
% Distance9 = sqrt(xd9^2 + yd9^2 + zd9^2);
% 
% xd10 = C5(1,1)- C4(1,1);
% yd10 = C5(1,2)- C4(1,2);
% zd10 = C5(1,3)- C4(1,3);
% 
% Distance10 = sqrt(xd10^2 + yd10^2 + zd10^2);
% 
% 
% 
% 
% 
% 
% 









