function [] = MMPplot( signals_c1, signals_c2)
%MMPplot_int creates a difference plot for emg signals of two conditions
%
%   written by Fabian Hoitz september 9th 2017
%
%   This function creates a multi muscle patter (MMP) plot where the first
%   coloumn shows the emg signals of all muscles in the first conditon,
%   the second coloumn the emg singnals of all muscles in the second
%   condition, and the third colum the difference between those two for
%   each muscle.
%   It is not recommended to have more than 6 muscles displayed at the same
%   time.
%   Necessary additonal files:
%       subplot1.m
%       colormapDifference.txt
%
%   INPUT:
%       emg signals of condition 1
%           (cell array where each cell contains a matrix that holds x
%           wavelets and their intensities per time point)
%       emg signals of condition 2 - same length as first one
%           (cell array where each cell contains a matrix that holds x
%           wavelets and their intensities per time point)
%   OUTPUT:
%       None
%
%   EXAMPLE:
%       MMPplot_int(signals_con1, signals_con2);

num_muscles = length(signals_c1);

for i = 1 : num_muscles
    max_c1(i) = max(max(signals_c1{i}, [], 2)');
    max_c2(i) = max(max(signals_c2{i}, [], 2)');
    max_max(i) = max([max_c1(i), max_c2(i)]);
    difference{i} = signals_c1{i} - signals_c2{i};
end

% create figure, define subplots, set background color and load/set colormap
figure
subplot1(num_muscles, 3, 'Gap', [0.04 0.0], 'YTickL', 'None', 'XTickL', 'None')
set(gcf,'color', [19/255 55/255 160/255]);
load('colormapDifference.txt');
colormap(colormapDifference);

% indices of subplots for each condition
c1_idx = 1:3:num_muscles*3;
c2_idx = 2:3:num_muscles*3;
c3_idx = 3:3:num_muscles*3;

% iterate over muscles
for i = 1 : num_muscles
    % iterate over conditions
    for ii = 1 : 3
        
        switch ii % switch input depending on condition
            case 1
                con_idx = c1_idx(i);
                signal = signals_c1{i};
            case 2
                con_idx = c2_idx(i);
                signal = signals_c2{i};
            case 3
                con_idx = c3_idx(i);
                signal = difference{i};
        end
        
        subplot1(con_idx); % put focus to current subplot
        contourf(signal, 20,'LineStyle', 'none'); % plot
        set(gca,'XColor',[1 1 1],'YColor',[1 1 1]);
        if ii == 3
            caxis([max(max(difference{i}))*-1 max(max(difference{i}))]);
        else
            caxis([0 max_max(i)]);
            colormap(gca, 'jet')
        end
        hold on;
        line([length(signal)/2 length(signal)/2],get(gca, 'YLim'),'color','r','linewidth',2);
        hold off;
    end
end

% annotation('textbox',[0.01 0.9 0 0],'String','Flexor','FontSize',30, 'Color', 'w');
% annotation('textbox',[0.01 0.725 0 0],'String','Extensor','FontSize',30, 'Color', 'w');
% annotation('textbox',[0.01 0.550 0 0],'String','Triceps','FontSize',30, 'Color', 'w');
% annotation('textbox',[0.01 0.375 0 0],'String','Deltoid','FontSize',30, 'Color', 'w');
% annotation('textbox',[0.01 0.200 0 0],'String','Gastroc','FontSize',30, 'Color', 'w');
% annotation('textbox',[0.2 1 0 0],'String','Hit','FontSize',30, 'Color', 'w');
% annotation('textbox',[0.475 1 0 0],'String','Miss','FontSize',30, 'Color', 'w');
% annotation('textbox',[0.74 1 0 0],'String','Difference','FontSize',30, 'Color', 'w');
end

