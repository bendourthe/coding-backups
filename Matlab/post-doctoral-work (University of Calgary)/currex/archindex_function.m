function [AILeft,AIRight] = archindex_function(data)
%% DESCRIPTION:
%
% archIdx = archindex_function(data)
%
% Calculates the average foot print from static Currex data (video files),
% then divided the foot print into 3 areas (forefoot, midfoot and heel) and
% calculates the Arch Index (= number of active cells in the arch divided
% by the overall number of active cells in the foot
%
%   Author: Sasa Cigoja
%
%   Last update: Dec. 7th, 2017
%
%% Input: 
%   - data: static video file recorded by Currex
%
%% Output:
%   - AI: Arch index for the corresponding static trial (left and right)

%% Function

data=data(:,2:3); % gets rid of the frames and pressure columns
data=unique(data,'rows','sorted'); % finds all unique rows taking x- and y-coordinates of the active cells in account
fig = figure;
maxfig(fig,1);
plot(data(:,1),data(:,2),'*') % plots all active cells
[x y]=getpts; % select 4 points from the pressure plot in the order described below
% 1) left heel
% 2) left MTP
% 3) right heel
% 4) right MTP
close all

ind=find(data(:,1)<21); % finds all cells with an x-coordinate lower than 21 -> indicating the left foot
ind=(ind(end));
left=data(1:ind,:); % determines all active sensors of left foot
right=data(ind+1:end,:); % determines all active sensors of right foot

%% sort coordinates
heel.left(1,1)=round(x(1,1)); % x-coordinate of left heel
heel.left(1,2)=round(y(1,1)); % y-coordinate of left heel
MTP.left(1,1)=round(x(2,1)); % x-coordinate of left MTP
MTP.left(1,2)=round(y(2,1)); % y-coordinate of left MTP

heel.right(1,1)=round(x(3,1)); % x-coordinate of right heel
heel.right(1,2)=round(y(3,1)); % y-coordinate of right heel
MTP.right(1,1)=round(x(4,1)); % x-coordinate of right MTP
MTP.right(1,2)=round(y(4,1)); % y-coordinate of right MTP

footlength.left=heel.left(1,2)-MTP.left(1,2); % left footlength in currex coordinates
footlength.right=heel.left(1,2)-MTP.right(1,2); % right footlength in currex coordinates

archstart.left=heel.left(1,2)-round(footlength.left/3); % the beginning of the arch is defined as subtracting a third of the footlength starting from the heel
archstart.right=heel.right(1,2)-round(footlength.right/3); % same for the right foot

archend.left=archstart.left-round(footlength.left/3); % the end of the arch is defined as subtracting a third of the footlength starting from the arch
archend.right=archstart.right-round(footlength.right/3); % same for the right foot

%%
left=sort(left,'descend'); % sorts the active cells of a foot in descending order; focus is on the second column -> y-coordinates
right=sort(right,'descend');

ind=find(left(:,2)<MTP.left(1,2)); % finds all active cells whose y-coordinate is smaller than the selected MTP -> those cells represent the toes
left(ind,:)=[]; % deletes the cell toes

ind=find(right(:,2)<MTP.right(1,2)); % same for the right foot
right(ind,:)=[];

cells.act.left=length(left(:,2)); % determines all active cells of the respective foot
cells.act.right=length(right(:,2));

heel.act.left=find(left(:,2)>=archstart.left); % finds all active cells in the heel whose y-coordinates are greater than the beginning of the arch
heel.act.left=length(heel.act.left);

heel.act.right=find(right(:,2)>=archstart.right); 
heel.act.right=length(heel.act.right);

arch.act.left=find(left(:,2)>=archend.left); % finds all active cells in the arch whose y-coordinates are greater than the end of the arch
arch.act.left=length(arch.act.left)-heel.act.left;

arch.act.right=find(right(:,2)>=archend.right);
arch.act.right=length(arch.act.right)-heel.act.right;

forefoot.act.left=find(left(:,2)>=MTP.left(1,2)); % finds all active cells in the forefoot whose y-coordinates are greater than the selected end of the MTP
forefoot.act.left=length(forefoot.act.left)-heel.act.left-arch.act.left;

forefoot.act.right=find(right(:,2)>=MTP.right(1,2));
forefoot.act.right=length(forefoot.act.right)-heel.act.right-arch.act.right;

AILeft=arch.act.left/cells.act.left; % calculates Arch Index based on the number of active cells in the arch divided by the overall number of active cells in the foot
AIRight=arch.act.right/cells.act.right;