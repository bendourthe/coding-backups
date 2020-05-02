clearvars
clc

%% load insole metrics
load(['\\CINCO\students\BennoGroup\Fabian\2018\' ...
    'Insoles - A visual demonstration of insole effects\'...
    'pedar_insole_metrics.mat']);

%% load peda data for insole 1 (CTRL)
data_c1a = importdata(['\\CINCO\students\BennoGroup\Benjamin_Dourthe\' ...
    'Dr. Scholls\Data Collection\Protocol V1\101\'...
    '1_101_ped_1.asc']);
data_c1a = data_c1a.data(2:end,100:end);

data_c1b = importdata(['\\CINCO\students\BennoGroup\Benjamin_Dourthe\' ...
    'Dr. Scholls\Data Collection\Protocol V1\101\'...
    '1_101_ped_2.asc']);
data_c1b = data_c1b.data(2:end,100:end);

%% load peda data for insole 2 (BPR)
data_c2a = importdata(['\\CINCO\students\BennoGroup\Benjamin_Dourthe\' ...
    'Dr. Scholls\Data Collection\Protocol V1\101\'...
    '2_101_ped_1.asc']);
data_c2a = data_c2a.data(2:end,100:end);

data_c2b = importdata(['\\CINCO\students\BennoGroup\Benjamin_Dourthe\' ...
    'Dr. Scholls\Data Collection\Protocol V1\101\'...
    '2_101_ped_2.asc']);
data_c2b = data_c2b.data(2:end,100:end);

%% load peda data for insole 3 (MG)
data_c3a = importdata(['\\CINCO\students\BennoGroup\Benjamin_Dourthe\' ...
    'Dr. Scholls\Data Collection\Protocol V1\101\'...
    '3_101_ped_1.asc']);
data_c3a = data_c3a.data(2:end,100:end);

data_c3b = importdata(['\\CINCO\students\BennoGroup\Benjamin_Dourthe\' ...
    'Dr. Scholls\Data Collection\Protocol V1\101\'...
    '3_101_ped_2.asc']);
data_c3b = data_c3b.data(2:end,100:end);

%% load peda data for insole 4 (CF)
data_c4a = importdata(['\\CINCO\students\BennoGroup\Benjamin_Dourthe\' ...
    'Dr. Scholls\Data Collection\Protocol V1\101\'...
    '4_101_ped_1.asc']);
data_c4a = data_c4a.data(2:end,100:end);

data_c4b = importdata(['\\CINCO\students\BennoGroup\Benjamin_Dourthe\' ...
    'Dr. Scholls\Data Collection\Protocol V1\101\'...
    '4_101_ped_2.asc']);
data_c4b = data_c4b.data(2:end,100:end);

%% prepare figure for plot
[x, y] = deal(insoles.YS(:,1) + 100, insoles.YS(:,2));
[xq, yq] = meshgrid(linspace(15, 100, 1000), linspace(0, 280, 1000));
z = zeros(1,99);
zq = griddata(x,y,z,xq,yq);

zq(zq == 0) = 1;

figure
colormap bone
subplot(1,3,2:3)
hold on
contourf(xq,yq,zq,1,'Linestyle','None');

clearvars x y xq yq z zq

%% processing
for cond = 1 : 8
    
    %% switch insole conditions
    switch cond
        case 1
            data = data_c1a;
        case 2
            data = data_c1b;
        case 3
            data = data_c2a;
        case 4
            data = data_c2b;
        case 5
            data = data_c3a;
        case 6
            data = data_c3b;
        case 7
            data = data_c4a;
        case 8
            data = data_c4b;
    end
    
    %% CoP calculation
    [cop_x, cop_y] = COP_trajectory_Pedar(data, insoles.YS(:, 3), insoles.YS(:, 1:2));
    CoP = [cop_x, cop_y];
    clearvars cop_x cop_y
    
    %% step detection
    % # of active cells per frame
    active_cells = WaveletFilter(smooth(sum(sign(data),2), 10)', 1/200, 0.5, 3, 1);
    
    % create binary signal to help identify steps
    bin_sig = zeros(length(active_cells), 1);
    
    % set frames with more than the mean amount of active cells to value 1
    bin_sig(active_cells > mean(active_cells)) = 1;
    
    % count sequences of ones (i.e. individual steps)
    [steps, total_steps] = bwlabel(bin_sig);
    
    % interpolate data for each step to 100 samples (ignore first and last step)
    for i = 2 : total_steps - 1
        
        curr_step = CoP(steps == i, :);
        
        CoP_stepwise(:,:,i-1) = interp1(curr_step, linspace(1, length(curr_step), 100));
        
    end
    
    clearvars curr_step i bin_sig active_cells steps total_steps
    
    %% calculate mean and std
    mean = mean(CoP_stepwise,3);
    std = std(CoP_stepwise,0, 3);
    
    %% add CoP trajectory to plot
    if cond <= 2
        h.(['p' num2str(cond, '%02d')]) = shadedErrorBar(mean(:,1)+100,...
            mean(:,2), std, {'-r', 'Linewidth', 1},1);
        %         plot(mean(:, 1) + 100, mean(:, 2), 'r')
    elseif cond > 2 && cond < 5
        h.(['p' num2str(cond, '%02d')]) = shadedErrorBar(mean(:,1)+100,...
            mean(:,2), std, {'-y', 'Linewidth', 1},1);
        %         plot(mean(:, 1) + 100, mean(:, 2), 'y')
    elseif cond > 4 && cond < 7
        h.(['p' num2str(cond, '%02d')]) = shadedErrorBar(mean(:,1)+100,...
            mean(:,2), std, {'-c', 'Linewidth', 1},1);
        %         plot(mean(:, 1) + 100, mean(:, 2), 'c')
    elseif cond >= 7
        h.(['p' num2str(cond, '%02d')]) = shadedErrorBar(mean(:,1)+100,...
            mean(:,2), std, {'-g', 'Linewidth', 1},1);
        %         plot(mean(:, 1) + 100, mean(:, 2), 'g')
    end
    
    clearvars mean std
end

%% beauty patches for figure
% box off
% xticks([25 35 45 55 65 75 85 95])
% xticklabels({'0', '10', '20', '30', '40', '50', '60', '70'})
% yticks([60 110 160 210 260])
% yticklabels({'50', '100', '150', '200', '250'})
% ax = gca;
% ax.FontSize = 15;
% ax.FontWeight = 'bold';
% ax.LineWidth = 3;
% ax.TickDir = 'out';
% ylabel({'Anterior/ Posterior';' displacement [mm]'}, 'FontWeight', 'bold', 'FontSize', 15)
% xlabel('Medial/ Lateral displacement [mm]', 'FontWeight', 'bold', 'FontSize', 15)
% rotateYLabel
% legend([h1.mainLine, h2.mainLine], 'Custom Fit', 'Back Pain Relief', 'Location', 'southwest')
