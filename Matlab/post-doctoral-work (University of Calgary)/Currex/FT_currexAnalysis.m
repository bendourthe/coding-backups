clc
clearvars

subject = 'Fabiand';
filePath = 'C:\Users\Fabian\Desktop\Master\CALGARY\Intern\Basketball\PilotData\Currex\';

data = csvread([filePath subject '\dVideoRight4.csv']);

%% Footline determination methods
%   1 - mean of 1/3 lateral/ medial points 
%   2 - mean of most anterior/ posterior points
%   3 - median of 1/3 anterior/ posterior points
%   4 - most intense of 1/3 anterior/ posterior points
pointEstimation = 4;

%% Animation
colormap jet
[~, startInd] = ismember(1:data(1,1), data(2:end,1), 'R2012a');
for i = 1 : data(1,1)
    if i == data(1,1)
        thisData = data(startInd(i)+1:end,2:4);
    else
        thisData = data(startInd(i)+1:startInd(i+1),2:4);
    end
    scatter(thisData(:,1),thisData(:,2), 25, thisData(:,3), 'filled');
    axis([0 64 0 64])
    line([32 32],[0 64])    %tLine
    drawnow
end

%% Foot wise angle assessment
% Get unique data points (e.g. Footprint)
uniqueData = zeros(size(unique(data(2:end,2:3),'rows'),1),3);
[uniqueData, ia, ic] = unique(data(2:end,2:3),'rows');
for i = 1 : size(ia,1)
    thisDataSet = find(ic == i) + 1;
    uniqueData(i,3) = mean(data(thisDataSet,4));
end
% uniqueData = uniqueData(uniqueData(:,3) > 100,:);

% cluster points to aquire left and right foot
index = kmeans(uniqueData(:,1:2), 2);
Foot1 = uniqueData(index == 1,:);
Foot2 = uniqueData(index == 2,:);

%% Different ways to find foot line
switch pointEstimation
    case 1
        
        F1totdist = max(Foot1(:,2)) - min(Foot1(:,2));
        F1top = Foot1(Foot1(:,2) >= (max(Foot1(:,2)) - F1totdist/3),:);
        F1bot = Foot1(Foot1(:,2) < (min(Foot1(:,2)) + F1totdist/3),:);
        
        F2totdist = max(Foot2(:,2)) - min(Foot2(:,2));
        F2top = Foot2(Foot2(:,2) >= (max(Foot2(:,2)) - F2totdist/3),:);
        F2bot = Foot2(Foot2(:,2) < (min(Foot2(:,2)) + F2totdist/3),:);
        
        f1p1 = [mean(F1top(:,1)), ...
            mean(F1top(:,2))];
        f1p2 = [mean(F1bot(:,1)), ...
            mean(F1bot(:,2))];
        
        f2p1 = [mean(F2top(:,1)), ...
            mean(F2top(:,2))];
        f2p2 = [mean(F2bot(:,1)), ...
            mean(F2bot(:,2))];

    case 2
        
        f1p1 = [mean(Foot1(find(Foot1(:,2) == max(Foot1(:,2))),1)), ...
            max(Foot1(:,2))];
        f1p2 = [mean(Foot1(find(Foot1(:,2) == min(Foot1(:,2))),1)), ...
            min(Foot1(:,2))];
        
        f2p1 = [mean(Foot2(find(Foot2(:,2) == max(Foot2(:,2))),1)), ...
            max(Foot2(:,2))];
        f2p2 = [mean(Foot2(find(Foot2(:,2) == min(Foot2(:,2))),1)), ...
            min(Foot2(:,2))]; 
        
    case 3
        
        F1totdist = max(Foot1(:,2)) - min(Foot1(:,2));
        F1top = Foot1(Foot1(:,2) >= (max(Foot1(:,2)) - F1totdist/3),:);
        F1bot = Foot1(Foot1(:,2) < (min(Foot1(:,2)) + F1totdist/3),:);
        
        F2totdist = max(Foot2(:,2)) - min(Foot2(:,2));
        F2top = Foot2(Foot2(:,2) >= (max(Foot2(:,2)) - F2totdist/3),:);
        F2bot = Foot2(Foot2(:,2) < (min(Foot2(:,2)) + F2totdist/3),:);
        
        f1p1 = [median(F1top(:,1)), ...
            median(F1top(:,2))];
        f1p2 = [median(F1bot(:,1)), ...
            median(F1bot(:,2))];
        
        f2p1 = [median(F2top(:,1)), ...
            median(F2top(:,2))];
        f2p2 = [median(F2bot(:,1)), ...
            median(F2bot(:,2))];
        
    case 4
        
        F1totdist = max(Foot1(:,2)) - min(Foot1(:,2));
        F1top = Foot1(Foot1(:,2) >= (max(Foot1(:,2)) - F1totdist/3),:);
        F1bot = Foot1(Foot1(:,2) < (min(Foot1(:,2)) + F1totdist/3),:);
        
        F2totdist = max(Foot2(:,2)) - min(Foot2(:,2));
        F2top = Foot2(Foot2(:,2) >= (max(Foot2(:,2)) - F2totdist/3),:);
        F2bot = Foot2(Foot2(:,2) < (min(Foot2(:,2)) + F2totdist/3),:);
        
        [~, maxF1top] = max(F1top(:,3));
        [~, maxF1bot] = max(F1bot(:,3));
        [~, maxF2top] = max(F2top(:,3));
        [~, maxF2bot] = max(F2bot(:,3));
        
        f1p1 = [F1top(maxF1top,1), ...
            F1top(maxF1top,2)];
        f1p2 = [F1bot(maxF1bot,1), ...
            F1bot(maxF1bot,2)];
        
        f2p1 = [F2top(maxF2top,1), ...
            F2top(maxF2top,2)];
        f2p2 = [F2bot(maxF2bot,1), ...
            F2bot(maxF2bot,2)];
end
%% Determine foot line
tl1 = [32, 0];
tl2 = [32, 64];

f1 = f1p2 - f1p1;
f2 = f2p2 - f2p1;
tl = tl2 - tl1;

angleF1 = 180 - rad2deg(acos(dot(f1, tl)/sqrt(sum(f1.^2) * sum(tl.^2))))
angleF2 = 180 - rad2deg(acos(dot(f2, tl)/sqrt(sum(f2.^2) * sum(tl.^2))))

figure
colormap jet
% scatter(Foot1(:,1),Foot1(:,2), 25, [1 0 0])
scatter(Foot1(:,1),Foot1(:,2), 25, Foot1(:,3), '*');
hold on
% scatter(Foot2(:,1),Foot2(:,2), 25, [0 0.5 0])
scatter(Foot2(:,1),Foot2(:,2), 25, Foot2(:,3), 'filled');
line([tl1(1) tl2(1)],[tl1(2) tl2(2)], 'linestyle', '--', 'LineWidth', 1, 'color', 'k')    %tLine
line([f1p1(1) f1p2(1)], [f1p1(2) f1p2(2)], 'LineWidth', 4, 'color', 'r')
text(f1p1(1) - 5, f1p1(2) + 10, num2str(angleF1))
line([f2p1(1) f2p2(1)], [f2p1(2) f2p2(2)], 'LineWidth', 4, 'color', 'r')
text(f2p1(1) - 5, f2p1(2) + 10, num2str(angleF2))
axis([0 64 0 64])


