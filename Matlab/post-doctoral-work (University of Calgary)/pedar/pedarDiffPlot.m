function [] = pedarDiffPlot(cellSheet, areaSheet, insole, cond1, cond2)
%PEDARDIFFPLOT creates a difference plot for for the pedar insole system
%
%   written by Fabian Hoitz, October 16th 2017
%
%   This function creates a difference plot for data created by the Pedar
%   insole system. Each subplot displayes the right insole
%   pressure cells with their load in kPa.
%   Necessary additional files:
%       colormapDifference.txt
%
%   INPUT:
%       Path to the excel file containing the cell coordinates
%           (ask Jordyn if you don't have it)
%       Path to the excel file containing the area of each cell
%           (ask Jordyn if you don't have it)
%       Insole condition
%       Path to data for first trial
%           (.asc file created by pedar)
%       Path to data for second trial
%           (.asc file created by pedar)
%
%   OUTPUT:
%       None
%
%   EXAMPLE:
%       cellsheet = '...\Centroid Calculation_All Insoles_B.xlsx';
%       areasheet = '...\Pedar_Insole_area_B.xlsx';
%       insole = 'YS Insole';
%       cond1 = '...\DrScholls_SD_1.asc';
%       cond2 = '...\DrScholls_SD_2.asc';
%       pedarDiffPlot(cellsheet, areasheet, insole, cond1, cond2)

%% coordinates for all 99 cells
num = xlsread(cellSheet,insole);

xy_coordinates = num(:,15:16);
xy_left = [xy_coordinates(:,1) * -1, xy_coordinates(:,2)];
xy_right = [xy_coordinates(:,1) + 200, xy_coordinates(:,2)];
clearvars num xy_coordinates

%% area of each cell
num = xlsread(areaSheet, insole);

area = num(3:end, 2) * 250000;
clearvars num

%% load first condition
insole_data = importdata(cond1);
data_kpa_left1 = insole_data.data(2:end,1 : 99);
data_kpa_right1 = insole_data.data(2:end, 100 : 198);
kpa_left1 = mean(insole_data.data(2:end,1 : 99));
kpa_right1 = mean(insole_data.data(2:end, 100 : 198));
clearvars insole_data

%% load second condition
insole_data = importdata(cond2);
data_kpa_left2 = insole_data.data(2:end,1 : 99);
data_kpa_right2 = insole_data.data(2:end, 100 : 198);
kpa_left2 = mean(insole_data.data(2:end,1 : 99));
kpa_right2 = mean(insole_data.data(2:end, 100 : 198)) * 1.26;
clearvars insole_data

%% calculate difference
kpa_left3 = kpa_left1 - kpa_left2;
kpa_right3 = kpa_right1 - kpa_right2;

%% define what side was measured (or the dominant side if both feet were measured)
if mean(mean([kpa_left1;kpa_left2])) > mean(mean([kpa_right1;kpa_right2]))
    kpa_dom = 'kpa_left';
else
    kpa_dom = 'kpa_right';
end

%% plotting section
max_diff = max([max(abs(kpa_left3)) max(abs(kpa_right3))]);
max_cond = max([max(kpa_left1) max(kpa_left2) max(kpa_right1) max(kpa_right2)]);

figure
% set(gcf,'color', [19/255 55/255 160/255]);
load('colormapDifference.txt');
colormap(colormapDifference);

for i = 1 : 3
    
    subplot(1,3,i)
    
    [x, y] = deal(xy_right(:,1), xy_right(:,2));
    z = eval([kpa_dom num2str(i)]);
    
    [xq, yq] = meshgrid(linspace(0, 300, 1000), linspace(0, 300, 1000));
    zq = griddata(x,y,z,xq,yq);
    
    contourf(xq,yq,zq,100,'Linestyle','None');
    colorbar('Ticks', [])
    
    xlim([110 200])
    
    if i == 3
        caxis([max_diff * -1 max_diff]);
    else
        caxis([0 max_cond]);
        colormap(gca, 'jet')
    end
    
    axis off
end
end

