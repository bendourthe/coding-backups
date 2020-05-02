function [] = pedarDiffPlot_v2(insoleData, cond1, cond2, varargin)
%% DESCRIPTION:
%
% [outputVariables] = FunctionName(inputVariables, varargin)
%
% Quick description of what the function does
%
%   Author:
%
%   Last update: Feb. 2nd, 2018
%
%% Input:
%   - List of input variables (if possible add their size into brackets;
%   e.g. matrix size m x n)
%
%% Output:
%   - List of output variables (if possible add their size into brackets;
%   e.g. matrix size m x n)
%
%% Dependencies:
%   Functions:
%       - List of .m functions that are necessary to run the code (e.g.
%       WaveletFilter.m)
%   Files:
%       - List of files that are necessary to run the code (e.g. Centroid
%       Calculation_All Insoles_B.xlsx)
%
%% Set default inputs parameters for name-value pairs
% adjust according to your code/ needs
defaults = struct(...
    'insole', 'YS' ...  % insole condition (i.e. WS, YS, etc.)
    );

% match passed name-value pairs to default pairs
params = inputHandler(defaults, varargin, 'pedarDiffPlot_v2');
clearvars defaults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% coordinates and area of cells
if ~any([strcmp('YS', params.insole), strcmp('XS', params.insole), ...
        strcmp('VS', params.insole), strcmp('WS', params.insole)])
    
    error('specified insole could not be found')
    
end

switch params.insole
    case 'YS'
        xy_coordinates = insoleData.YS(:, 1:2);
        area = insoleData.YS(:, 3) * 250000;
    case 'WS'
        xy_coordinates = insoleData.WS(:, 1:2);
        area = insoleData.WS(:, 3) * 250000;
    case 'XS'
        xy_coordinates = insoleData.XS(:, 1:2);
        area = insoleData.XS(:, 3) * 250000;
    case 'VS'
        xy_coordinates = insoleData.VS(:, 1:2);
        area = insoleData.VS(:, 3) * 250000;
end

xy_left = [xy_coordinates(:,1) * -1, xy_coordinates(:,2)];
xy_right = [xy_coordinates(:,1) + 200, xy_coordinates(:,2)];
clearvars xy_coordinates

%% calculate difference
cond3 = cond1 - cond2;

%% plotting
max_diff = max(abs(cond3));
max_cond = max([cond1, cond2]);

figure

if exist('colormapDifference.txt', 'file') ~= 2
    error('colormapdifference.txt could not be found')
else
    load('colormapdifference.txt')
    colormap(colormapdifference)
end

for i = 1 : 3
    
    subplot(1,3,i)
    
    [x, y] = deal(xy_right(:,1), xy_right(:,2));
    z = eval(['cond' num2str(i)]);
    
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
    
    switch i 
        case 1
            annotation('textbox', [0.19 0.51 0.3 0.3], 'String', 'Insole A', ...
                'FitBoxToText', 'on', 'LineStyle', 'none', 'FontSize', 18, 'FontWeight', 'bold')
        case 2
            annotation('textbox', [0.47 0.51 0.3 0.3], 'String', 'Insole B', ...
                'FitBoxToText', 'on', 'LineStyle', 'none', 'FontSize', 18, 'FontWeight', 'bold')
        case 3
            annotation('textbox', [0.74 0.51 0.3 0.3], 'String', 'Difference', ...
                'FitBoxToText', 'on', 'LineStyle', 'none', 'FontSize', 18, 'FontWeight', 'bold')
    end
    
    axis off
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


%% InputHandler Function - DO NOT CHANGE
% function to simplify the way name/ value pairs get handled when passed
% as optional inputs
function [defaults] = inputHandler(defaults,inputs,functionName)
%INPUTHANDLER can be used to handle inputs for functions that use multiple
%   name/value pairs as input.
%
%   Programmed by Fabian Hoitz, April 2018
%
%   Input:
%       defaults (structur holding default values of higher-order function)
%       inputs (varargin variable of higher-order function)
%       functionName (string containing name of higher-order function)
%
%   Output:
%       Structure containing the input parameters of higher-order function
%
%   Example:
%       parameters = inputHandler(defaults, varargin, 'functionName');

% Names of parameters
paramNames = fieldnames(defaults);

% Number of arguments passed 
nArgs = length(inputs);

% warning if name-value pairs are not within default structure
if round(nArgs/2) ~= nArgs/2
    error([functionName ' needs propertyName/propertyValue pairs'])
end

for pair = reshape(inputs,2,[]) % pair is {propName;propValue}
    
    inpName = pair{1};
    %     inpName = lower(pair{1}); % make case insensitive
    
    if any(strcmp(inpName,paramNames))
        % overwrite options. If you want you can test for the right class here
        % Also, if you find out that there is an option you keep getting wrong,
        % you can use "if strcmp(inpName,'problemOption'),testMore,end"-statements
        defaults.(inpName) = pair{2};
    else
        error('%s is not a recognized parameter name',inpName)
    end
end
end