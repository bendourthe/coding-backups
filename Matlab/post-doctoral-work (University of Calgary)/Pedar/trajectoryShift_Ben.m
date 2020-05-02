%% setup
clearvars
close all
clc

subject = 2;
insole = 'YS';

%% load relevant data
% insole metrics
load(['\\CINCO\students\BennoGroup\Fabian\2018\' ...
    'Insoles - A visual demonstration of insole effects\'...
    'pedar_insole_metrics.mat']);

% insole data
path = fullfile('\\CINCO\students\BennoGroup\Benjamin_Dourthe\Pedar', ...
    num2str(subject, '%03d'));

ctrl_fatigued = importdata(fullfile(path, ['1_' num2str(subject, '%03d') ...
    '_pedSD_2.asc']));
ctrl_unfatigued = importdata(fullfile(path, ['1_' num2str(subject, '%03d') ...
    '_pedSD_1.asc']));
cf_fatigued = importdata(fullfile(path, ['2_' num2str(subject, '%03d') ...
    '_pedSD_2.asc']));
cf_unfatigued = importdata(fullfile(path, ['2_' num2str(subject, '%03d') ...
    '_pedSD_1.asc']));

figureName = fullfile(path, ['plot' num2str(subject, '%03d') ...
    '.fig']);

%% determine active 99 cells per condition and calculate COP trajectory
for i = 1 : 4
    % switch conditions
    switch i
        case 1
            c = ctrl_fatigued;
        case 2
            c = ctrl_unfatigued;
        case 3
            c = cf_fatigued;
        case 4
            c = cf_unfatigued;
    end
    
    % determine which side has higher pressure values (i.e. which insole was used)
    sumLeft = sum(sum(c.data(2:end, 1:99)));
    sumRight = sum(sum(c.data(2:end, 100:end)));
    
    if sumLeft > sumRight
        c = c.data(2:end, 1:99);
        side = 'left';
    elseif sumLeft < sumRight
        c = c.data(2:end, 100:end);
        side = 'right';
    end
    clearvars sumLeft sumRight
    
    % calculate COP trajectory
    switch insole 
        case 'XS'
            in = insoles.XS;
        case 'YS'
            in = insoles.YS;
    end
    
    [cop_x, cop_y] = COP_trajectory_Pedar(c, in(:, 3), in(:, 1:2));
    CoP = [cop_x, cop_y];
    clearvars cop_x cop_y
    
    % step detection
    active_cells = c;
    active_cells(active_cells ~= 0) = 1;
    active_cells = sum(active_cells, 2);
    active_cells(active_cells <= 15) = 0;
    
    bin_sig = zeros(length(active_cells), 1);
    bin_sig(active_cells == 0) = 1;
    
    [L, NUM] = bwlabel(bin_sig);
    for ii = 1 : NUM
        idx = find(L == ii);
        if length(idx) <= 15
            bin_sig(idx) = 0;
        end
    end
    
    [steps, stepN] = bwlabel(~bin_sig);
    
    clearvars ii idx active_cells bin_sig L NUM
    
    % interpolate data for random 50 steps to 100 samples
    for ii = i : stepN
        
        curr_step = CoP(steps == ii, :);
        curr_step(~any(curr_step, 2), :) = [];
        
        CoP100(:,:,ii) = interp1(curr_step, linspace(1, length(curr_step), 100));
        
    end
    clearvars curr_step ii
    
    % delete 2d parts that are all 0
    idx0 = logical(squeeze(sum(~any(CoP100))));
    CoP100(:,:, idx0) = [];
    N = size(CoP100, 3);
    
    switch i
        case 1
            ctrl_fatigued = struct('data', c, 'side', side, ...
                'CoP', CoP100, 'N', N);
        case 2
            ctrl_unfatigued = struct('data', c, 'side', side, ...
                'CoP', CoP100, 'N', N);
        case 3
            cf_fatigued = struct('data', c, 'side', side, ...
                'CoP', CoP100, 'N', N);
        case 4
            cf_unfatigued = struct('data', c, 'side', side, ...
                'CoP', CoP100, 'N', N);
    end
    
    clearvars c CoP100 step50 steps side CoP stepN idx0 N
end

%% plotting
% insole outline
[x, y] = deal(in(:,1) + 100, in(:,2));
[xq, yq] = meshgrid(linspace(0, 100, 1000), linspace(0, 270, 1000));
z = zeros(1,99);
zq = griddata(x,y,z,xq,yq);

zq(zq == 0) = 1;

figure
colormap bone
subplot(1,3,2:3)
hold on
contourf(xq,yq,zq,1,'Linestyle','None');

clearvars x y xq yq z zq

% cop trajectories
mean = mean(ctrl_unfatigued.CoP, 3);
std = std(ctrl_unfatigued.CoP, 0, 3);
h1 = shadedErrorBar(mean(:,1)+100, mean(:,2), std, {'-', 'Linewidth', 2, 'Color', [.7 .7 .7]},1);
clearvars mean std

mean = mean(cf_unfatigued.CoP, 3);
std = std(cf_unfatigued.CoP, 0, 3);
h2 = shadedErrorBar(mean(:,1)+100, mean(:,2), std, {'-g', 'Linewidth', 2},1);
clearvars mean std

mean = mean(ctrl_fatigued.CoP, 3);
std = std(ctrl_fatigued.CoP, 0, 3);
h3 = shadedErrorBar(mean(:,1)+100, mean(:,2), std, {'--', 'Linewidth', 2, 'Color', [.7 .7 .7]},1);
clearvars mean std

mean = mean(cf_fatigued.CoP, 3);
std = std(cf_fatigued.CoP, 0, 3);
h4 = shadedErrorBar(mean(:,1)+100, mean(:,2), std, {'--g', 'Linewidth', 2},1);
clearvars mean std

title(['Subject ' num2str(subject)])
box off
xticks([20 30 40 50 60 70 80 90])
xticklabels({'0', '10', '20', '30', '40', '50', '60', '70'})
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
legend([h1.mainLine, h2.mainLine, h3.mainLine, h4.mainLine], ...
    ['CTRL - No Fatigue (N = ' num2str(ctrl_unfatigued.N) ')'], ...
    ['CF - No Fatigue (N = ' num2str(cf_unfatigued.N) ')'], ...
    ['CTRL - Fatigue (N = ' num2str(ctrl_fatigued.N) ')'], ...
    ['CF - Fatigue (N = ' num2str(cf_unfatigued.N) ')'], ...
    'Location', 'southwest')
savefig(gcf, figureName)