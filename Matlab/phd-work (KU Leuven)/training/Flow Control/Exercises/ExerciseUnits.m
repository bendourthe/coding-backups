%% Metric unit converter
clear all
close all
clc
%% Part 1: contverter
% Write a program to convert metric units to centimeters. Use the SWITCH 
% command to define different cases e.g. mm, cm, dm, m, km. Test your code
% for different units and values.
x = 3.0;        % numeric variable
units = 'mm';   % string variable
switch units
    case{'millimeter','mm'}% multiple matches
        y = x/10;          % converts to centimeters
        disp([num2str(x),' ',units,' converted to cm is: ',num2str(y)])
    case{'cm','centimeter'}
        y = x;
        disp([num2str(x),' ',units,' converted to cm is: ',num2str(y)])
    case{'dm','decimeter'} 
        y = x*10;           
        disp([num2str(x),' ',units,' converted to cm is: ',num2str(y)])
    case{'m','meter'}
        y = x*100;         
        disp([num2str(x),' ',units,' converted to cm is: ',num2str(y)])
    case{'km','kilometer'}
        y = x*100000;
        disp([num2str(x),' ',units,' converted to cm is: ',num2str(y)])
    otherwise
        disp(['unknown units: ',units])
        y = NaN;  % not a number
end

%% Part 2: Extension with more data
% Copy the previous code and adapt it to use this program for converting
% arrays of numbers and units. Use a FOR loop.
X = [0.0200, 3.45, 447, 8.80, 32.7, 6.5, 1.6747, 36, 3245, 163];
Units = {'m','dm','mm','cm','centimeter','meter','decimeter','millimeter','mm','cm'};

% Create a zero-matrix with the same size as X.
Y = zeros(size(X));
for i = 1:length(X);
    switch Units{i}
        case{'millimeter','mm'}	% multiple matches
            Y(i) = X(i)/10;     % converts to centimeters
            disp([num2str(X(i)),' ',Units{i},' converted to cm is: ',num2str(Y(i))])
        case{'cm','centimeter'}
            Y(i) = X(i);
            disp([num2str(X(i)),' ',Units{i},' converted to cm is: ',num2str(Y(i))])
        case{'dm','decimeter'} 
            Y(i) = X(i)*10;           
            disp([num2str(X(i)),' ',Units{i},' converted to cm is: ',num2str(Y(i))])
        case{'m','meter'}
            Y(i) = X(i)*100;        
            disp([num2str(X(i)),' ',Units{i},' converted to cm is: ',num2str(Y(i))])
        case{'km','kilometer'}
            Y(i) = X(i)*100000;
            disp([num2str(X(i)),' ',Units{i},' converted to cm is: ',num2str(Y(i))])
        otherwise
            disp(['unknown units: ',Units{i}])
            Y(i) = NaN;  % not a number
    end
end