clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Curvature_Bones:
%   Calculates the mean, max, and min curvatures in the two principal
%   directions of a surface (e.g. articular surface). Uses a 3D-screen
%   which:
%       - scans the surface along both directions
%       - calculates the best circle fit in each frame
%       - calculates the curvature for each fram as the inverse of the
%       radius of each circle
%       - calculates the mean, max and min curvature along both directions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dependencies:        
%               STL_ReadFile.m 
%               GUI_PlotShells
%               PlacePoints3
%               CircleFitByPratt
%               circle_3D
% Input: 
%           Surf: STL file corresponding to the surface to study (e.g.
%           articular surface of a bone)
%           n_frame: number of frames the screen will scan
% 
% Output:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Selection of the current directory (where the STL files are)
dir = 'H:\Data\Study_Brown_Rig\Segmentation\Articular_Surfaces\';
% Available data
cs_available = 0;          % loads data if available (=1), else, generate data
    filenamecs = 'cs_scan26_trap.txt';
    filecs = [dir filenamecs];
    
% Number of screening for curvature calculation (i.e. how many frames will
% be created to calculate the curvature -> the higher this value is, the
% more the mean, max and min curvatures will be accurate)
n_frame = 50;

% Bone:
mc1 = 0;    % = 1 if studying the MC1, = 0 otherwise
trap = 1;   % = 1 if studying the trapezium, = 0 otherwise

% FILES
Bone = 'SCAN3_N_mc1.stl';
Surf = 'SCAN3_N_mc1_AS.stl';
    
% Read STL
[F, V] =  STL_ReadFile([dir Bone],true);
[F1, V1] =  STL_ReadFile([dir Surf],true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coordinate System (cs) on surface:
%   origin: centroid of the surface (calculated automatically).
%   x: in the dorsal-volar direction
%   y: in the distal-proximal direction%   
%   z: in the ulnar-radial direction

if(cs_available)
    A = txt2mat(filecs);
    
    X1 = A(1,:);    % dorsal tip (MC1) or notch (trapezium)
    X2 = A(2,:);    % volar tip (or palmar beak, MC1) or volar notch (trapezium)
    Y = A(3,:);     % deepest (MC1) or highest point located along the dorsal-volar axis
else
    [ad_curve] = PlacePoints3({F}, {V}, {F}, {V}, 'Select 3 landmarks (dorsal tip/notch, volar tip/notch and center, in that order!). End by pressing enter and exiting the figure');
    close all;
    % Retrieves coordinates of the points:
    X1 = ad_curve.points{end,1}(1,:);
    X2 = ad_curve.points{end,1}(2,:);
    Y = ad_curve.points{end,1}(3,:);
    
    % Save in a txt file
    fid = fopen(filecs,'w+');
    fprintf(fid,'%5.8f\t%5.8f\t%5.8f\r\n',[X1,X2,Y]);
    fclose(fid);
end

centr = mean(V1);
for i=1:length(V1)
    dif = V1(i,:) - centr;
end
[D, idx_o] = min(dif);
O = (X1+X2)/2;
Xax = (X2-X1)/norm(X2-X1);
Yax = (O-Y)/norm(O-Y);
Zax = cross(Xax,Yax)/norm(cross(Xax,Yax));
Yax = cross(Zax,Xax)/norm(cross(Zax,Xax));
cs_trap = [Xax; Yax; Zax];
O = mean(V1);

% Converges to new CS
RotMat = [Xax', Yax', Zax'];

   % These tranformation matrices convert from local to global CS.
   % Invert them to go from global to local.
glob2loc = [RotMat', -(RotMat')*O'; 0 0 0 1];
R_glob2loc = glob2loc(1:3,1:3);
T_glob2loc = glob2loc(1:3,4);

% Transformes data to new CS
    % Whole bone
    V = R_glob2loc*V.' + repmat(T_glob2loc,1,size(V,1));
    V = V.';
    % Surface only
    V1 = R_glob2loc*V1.' + repmat(T_glob2loc,1,size(V1,1));
    V1 = V1.';

if(cs_available)
else
% Plots the surface with the new coordinate system for validation
    h0 = figure;
    [obj, li, ax] = GUI_PlotShells(h0, {F;F1},{V;V1},{ones(size(V,1),1);ones(size(V1,1),1)});
    GUI_VisualisationUI(ancestor(ax, 'figure'), true, ax, true, true);
        
        % Plots the coordinate system and its origin
    hold on
    box off
    arrow([0;0;0],[30;0;0],10,70,30,5,'EdgeColor','r','FaceColor','r');
    arrow([0;0;0],[0;30;0],10,70,30,5,'EdgeColor','g','FaceColor','g');
    arrow([0;0;0],[0;0;30],10,70,30,5,'EdgeColor','b','FaceColor','b');
    set(gcf,'numbertitle','off','name','New Coordinate System on Surface');
    
        % Plots 2 buttons: 'Bad CS' if the coordinate system looks wrong,
        % this will stop the code, and 'Continue' if the coordinate system
        % looks right, will continue to run the rest of the code
    ans=1;
    b1 = uicontrol('Position',[20 20 200 40],'String','Continue',...
              'Callback','uiresume(gcbf)');
    b2 = uicontrol('Position',[20 60 200 40],'String','Bad CS',...
              'Callback','BreakOnClick(ans)');
    uiwait(gcf);
    if ans == 0
        close(h0)
        error('Incorrect coordinate system. Please try again!');
        break
    end
end

% Plots the mc1/trapezium with the new coordinate system
    h0 = figure;
    [obj, li, ax] = GUI_PlotShells(h0, {F;F1},{V;V1},{ones(size(V,1),1);ones(size(V1,1),1)});
    box off    
        % Plots the coordinate system and its origin
    hold on
    arrow([0;0;0],[30;0;0],10,70,30,5,'EdgeColor','r','FaceColor','r');
    arrow([0;0;0],[0;30;0],10,70,30,5,'EdgeColor','g','FaceColor','g');
    arrow([0;0;0],[0;0;30],10,70,30,5,'EdgeColor','b','FaceColor','b');
    set(gcf,'numbertitle','off','name','New Coordinate System on Surface');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the ranges along which the screen will scan the surface
M1 = max(V1(:,1));
m1 = min(V1(:,1));
M3 = max(V1(:,3));
m3 = min(V1(:,3));
r_dv = M1 - m1; % from the most dorsal (u) to the most volar (r) point
r_ur = M3 - m3; % from the most ulnar (u) to the most radial (r) point
s_ur = r_dv/n_frame; % define the step size in the dorsal/volar direction
s_dv = r_ur/n_frame; % define the step size in the ulnar/radial direction

% Ulnar-Radial Curvature
ur_curv = figure;
    [obj, li, ax] = GUI_PlotShells(ur_curv, {F1},{V1},{ones(size(V1,1),1)});
    hold on
    box off
for i=1:n_frame
    idx_ur = find(V1(:,1) < M1 - (i-1)*s_ur & V1(:,1) > M1 - i*s_ur);
    temp_ur = V1(idx_ur, 2:3);
    if numel(temp_ur) < 3; % cannot calculate a circle fit with less than 3 points
        curv_ur(i) = NaN;
    else
        fit_ur = CircleFitByPratt(temp_ur);
        Y_ur(i) = fit_ur(1);
        Z_ur(i) = fit_ur(2);
        X_ur(i) = (M1 - (i-1)*s_ur);
        R_ur(i) = fit_ur(3);
        if Y_ur(i) > 0
            curv_ur(i) = -1/R_ur(i);
        else
            curv_ur(i) = 1/R_ur(i);
        end
    end
    if mc1 == 1 && trap == 0;
    elseif mc1 == 0 && trap == 1;
            curv_ur(i) = - curv_ur(i);
    end
end
curv_ur = curv_ur(~isnan(curv_ur)); % Remove NaN values (usually around the
                                    % extremities of the surface where there
                                    % is only one point)
hold on
for i=2:numel(curv_ur)-1        
    h_circle = circle_3D(R_ur(i), [X_ur(i) Y_ur(i) Z_ur(i)], [1 0 0]);
end
set(gcf,'numbertitle','off','name','Set of fitting circles for ulnar-radial curvature');


% Dorsal-Volar Curvature
dv_curv = figure;
    [obj, li, ax] = GUI_PlotShells(dv_curv, {F1},{V1},{ones(size(V1,1),1)});
    hold on
    box off
for i=1:n_frame
    idx_dv = find(V1(:,3) < M3 - (i-1)*s_dv & V1(:,3) > M3 - i*s_dv);
    temp_dv = V1(idx_dv, 1:2);
    if numel(temp_dv) < 3; % cannot calculate a circle fit with less than 3 points
        curv_dv(i) = NaN;
    else
        fit_dv = CircleFitByPratt(temp_dv);
        X_dv(i) = fit_dv(1);
        Y_dv(i) = fit_dv(2);
        Z_dv(i) = (M3 - (i-1)*s_dv);
        R_dv(i) = fit_dv(3);
        if Y_dv(i) > 0
            curv_dv(i) = -1/R_dv(i);
        else
            curv_dv(i) = 1/R_dv(i);
        end
    end
    if mc1 == 1 && trap == 0;
    elseif mc1 == 0 && trap == 1;
            curv_dv(i) = - curv_dv(i);
    end
end
curv_dv = curv_dv(~isnan(curv_dv)); % Remove NaN values (usually around the
                                    % extremities of the surface where there
                                    % is only one point)
hold on
for i=2:numel(curv_dv)-1        
    h_circle = circle_3D(R_dv(i), [X_dv(i) Y_dv(i) Z_dv(i)], [0 0 1]);
end
set(gcf,'numbertitle','off','name','Set of fitting circles for dorsal-volar curvature');

% Final Results
ulnar_radial_mean_max_min = [mean(curv_ur) max(curv_ur) min(curv_ur)]
dorsal_volar_mean_max_min = [mean(curv_dv) max(curv_dv) min(curv_dv)]
        


