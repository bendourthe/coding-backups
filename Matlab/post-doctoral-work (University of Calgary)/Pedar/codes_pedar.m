%% setup
load(['\\CINCO\students\BennoGroup\Fabian\2018\' ...
    'Insoles - A visual demonstration of insole effects\'...
    'pedar_insole_metrics.mat']);

data_c1 = importdata(['\\CINCO\students\BennoGroup\Fabian\2018\' ...
    'Insoles - A visual demonstration of insole effects\'...
    '3_1_01_ped.asc']);
data_c1 = data_c1.data(2:end,1:99);

data_c2 = importdata(['\\CINCO\students\BennoGroup\Fabian\2018\' ...
    'Insoles - A visual demonstration of insole effects\'...
    '4_1_01_ped.asc']);
data_c2 = data_c2.data(2:end,1:99);

%% CoP calculation
for cond = 1 : 2
    
    switch cond
        case 1
            [cop_x, cop_y] = COP_trajectory_Pedar(data_c1, insoles.WS(:, 3), insoles.WS(:, 1:2));
            CoP_c1 = [cop_x, cop_y];
        case 2
            [cop_x, cop_y] = COP_trajectory_Pedar(data_c2, insoles.WS(:, 3), insoles.WS(:, 1:2));
            CoP_c2 = [cop_x, cop_y];
    end
end
clearvars cop_x cop_y cond

%% step detection
for cond = 1 : 2
    % count # of active cells per frame
    switch cond
        case 1
            active_cells = WaveletFilter(smooth(sum(sign(data_c1),2), 10)', 1/200, 0.5, 3, 1);
        case 2
            active_cells = WaveletFilter(smooth(sum(sign(data_c2),2), 10)', 1/200, 0.5, 3, 1);
    end
    
    % create binary signal to help identify steps
    bin_sig = zeros(length(active_cells), 1);
    
    % set frames with more than 10 active cells to value 1
    bin_sig(active_cells > mean(active_cells)) = 1;
    
    switch cond
        case 1
            % count sequences of ones (i.e. individual steps)
            [c1_steps, c1_total_steps] = bwlabel(bin_sig);
        case 2
            % count sequences of ones (i.e. individual steps)
            [c2_steps, c2_total_steps] = bwlabel(bin_sig);
    end
end
clearvars curr_data i bin_sig active_cells cond

%% interpolate CoP trajectory to 100 samples
for i = 1 : 2
    switch i
        case 1
            for ii = 2 : c1_total_steps - 1
                
                curr_step_data = CoP_c1(c1_steps == ii, :);
                               
                CoP_c1_all(:,:,ii-1) = interp1(curr_step_data, linspace(1, length(curr_step_data), 100));
                
            end
        case 2
            for ii = 2 : c2_total_steps - 1
                
                curr_step_data = CoP_c2(c2_steps == ii, :);
                
                CoP_c2_all(:,:,ii-1) = interp1(curr_step_data, linspace(1, length(curr_step_data), 100));
                
            end
    end
end
clearvars ii curr_step_data i

%% mean plot
c1_mean = mean(CoP_c1_all,3);
c1_std = std(CoP_c1_all,0, 3);

c2_mean = mean(CoP_c2_all,3);
c2_std = std(CoP_c2_all,0, 3);

[x, y] = deal(insoles.WS(:,1) + 100, insoles.WS(:,2));
[xq, yq] = meshgrid(linspace(10, 110, 1000), linspace(0, 250, 1000));
z = zeros(1,99);
zq = griddata(x,y,z,xq,yq);

zq(zq == 0) = 1;

colormap bone
subplot(1,3,2:3)
contourf(xq,yq,zq,1,'Linestyle','None');
hold on
h1 = shadedErrorBar(c1_mean(:,1)+100, c1_mean(:,2), c1_std, {'-r', 'Linewidth', 2},1);
h2 = shadedErrorBar(c2_mean(:,1)+100, c2_mean(:,2), c2_std, {'-b', 'Linewidth', 2},1);

box off
xlim([25 95])
xticks([25 35 45 55 65 75 85 95])
xticklabels({'0', '10', '20', '30', '40', '50', '60', '70'})
ylim([10 245])
yticks([60 110 160 210 260])
yticklabels({'50', '100', '150', '200', '250'})
ax = gca;
ax.FontSize = 15;
ax.FontWeight = 'bold';
ax.LineWidth = 3;
ax.TickDir = 'out';
ylabel({'Anterior/ Posterior';' displacement [mm]'}, 'FontWeight', 'bold', 'FontSize', 15)
xlabel('Medial/ Lateral displacement [mm]', 'FontWeight', 'bold', 'FontSize', 15)
rotateYLabel
legend([h1.mainLine, h2.mainLine], 'Custom Fit', 'Back Pain Relief', 'Location', 'southwest')

clearvars c1_mean c2_mean c1_std c2_std x y xq yq z zq
%% Pedar diff plot (percentage wise, i.e. @5%, 10%, 15%, etc.)
% close all
% interpolate pedar data of each step to 100 frames
for i = 1 : 2
    switch i
        case 1
            for ii = 2 : c1_total_steps - 1
                
                curr_step_data = data_c1(c1_steps == ii, :);
                               
                data_c1_all(:,:,ii-1) = interp1(curr_step_data, linspace(1, length(curr_step_data), 100));
                
            end
        case 2
            for ii = 2 : c2_total_steps - 1
                
                curr_step_data = data_c2(c2_steps == ii, :);
                
                data_c2_all(:,:,ii-1) = interp1(curr_step_data, linspace(1, length(curr_step_data), 100));
                
            end
    end
end
clearvars ii curr_step_data i

% build mean for each frame
data_c1_mean = nanmean(data_c1_all, 3);
data_c2_mean = nanmean(data_c2_all, 3);

writerObj = VideoWriter('insoleDifference1.avi');
writerObj.FrameRate = 10;
open(writerObj)

for i = 1 : size(data_c1_mean, 1)
    
    pedarDiffPlot_v2(insoles, data_c1_mean(i,:), data_c2_mean(i,:), 'insole', 'WS')
    subplot(1,3,2)
    title(['@ ' num2str(i) '% of stance phase'], 'FontSize', 20)
    
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    set(gcf, 'color', [.6 .6 .6])
    
    frame = getframe(gcf);
    writeVideo(writerObj, frame);
    close(gcf)
end

close(writerObj)
%% Anitmation
% [x, y] = deal(insoles.WS(:,1) + 100, insoles.WS(:,2));
% [xq, yq] = meshgrid(linspace(0, 250, 50), linspace(0, 250, 50));
% 
% for i = 1 : length(data_c1)
%     
%     disp(i)
%     z = data_c1(i,:);
%     zq = griddata(x,y,z,xq,yq);
%     
%     contourf(xq,yq,zq,50,'Linestyle','None');
%     drawnow
%     
% end