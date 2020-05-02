function AI = archindexDyn_function(data)
%% DESCRIPTION:
%
% AI = archindexDyn_function(data)
%
% Calculates the average foot print from dynamic Currex data (mean files),
% then divided the foot print into 3 areas (forefoot, midfoot and heel) and
% calculates the Arch Index (= number of active cells in the arch divided
% by the overall number of active cells in the foot)
%
%   Author: Sasa Cigoja
%
%   Last update: May 1st, 2017
%
%% Input: 
%   - data: dynamic mean file recorded by Currex
%
%% Output:
%   - AI: Arch index for the corresponding static trial (left or right)

%% Function

data=data(:,2:3); % gets rid of the frames and pressure columns
data=unique(data,'rows','sorted'); % finds all unique rows taking x- and y-coordinates of the active cells in account
fig = figure;
plot(data(:,1),data(:,2),'*') % plots all active cells
set(gcf,'position',[500,70,300,550])
[x y]=getpts; % select 2 points from the pressure plot in the order described below
% 1) heel
% 2) MTP
close all

foot=data;

%% sort coordinates
heel.foot(1,1)=round(x(1,1)); % x-coordinate of heel
heel.foot(1,2)=round(y(1,1)); % y-coordinate of heel
MTP.foot(1,1)=round(x(2,1)); % x-coordinate of MTP
MTP.foot(1,2)=round(y(2,1)); % y-coordinate of MTP

footlength.foot=heel.foot(1,2)-MTP.foot(1,2); % footlength in currex coordinates

archstart.foot=heel.foot(1,2)-round(footlength.foot/3); % the beginning of the arch is defined as subtracting a third of the footlength starting from the heel

archend.foot=archstart.foot-round(footlength.foot/3); % the end of the arch is defined as subtracting a third of the footlength starting from the arch

%%
foot=sort(foot,'descend'); % sorts the active cells of a foot in descending order; focus is on the second column -> y-coordinates

ind=find(foot(:,2)<MTP.foot(1,2)); % finds all active cells whose y-coordinate is smaller than the selected MTP -> those cells represent the toes
foot(ind,:)=[]; % deletes the cell toes

cells.act.foot=length(foot(:,2)); % determines all active cells of the respective foot

heel.act.foot=find(foot(:,2)>=archstart.foot); % finds all active cells in the heel whose y-coordinates are greater than the beginning of the arch
heel.act.foot=length(heel.act.foot);

arch.act.foot=find(foot(:,2)>=archend.foot); % finds all active cells in the arch whose y-coordinates are greater than the end of the arch
arch.act.foot=length(arch.act.foot)-heel.act.foot;

forefoot.act.foot=find(foot(:,2)>=MTP.foot(1,2)); % finds all active cells in the forefoot whose y-coordinates are greater than the selected end of the MTP
forefoot.act.foot=length(forefoot.act.foot)-heel.act.foot-arch.act.foot;

AI=arch.act.foot/cells.act.foot; % calculates Arch Index based on the number of active cells in the arch divided by the overall number of active cells in the foot