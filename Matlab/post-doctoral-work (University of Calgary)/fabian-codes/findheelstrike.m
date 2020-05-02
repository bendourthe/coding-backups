function [ heelstrikes ] = findheelstrike( signal )
%findheelstrike detect all heelstrikes in acceleration data
%   detects the index of heelstirkes (rise in acceleration due to impact)
%   in data collected from accelerometer mounted at the heel of foot.
% INPUT:
%   signal (1d vector of acceration data)
% OUTPUT:
%   heelstrikes (indices of heelstrikes detected in signal)
signal = smooth(signal,7);
baseline = mean(signal(1 : 2000));
zeroMeanSig = signal - baseline;
abs0MS = abs(zeroMeanSig);

[~, peaks] = findpeaks(abs0MS, ...
    'MinPeakHeight',0.02, 'MinPeakDistance', 1300,'MinPeakProminence', 0.03);

heelstrikes = peaks;

% % check if it is a double peak
% for i = 1 : length(peaks)
%     
%     if peaks(i) - 50 <= 0
%         continue
%     end
%     
%     [~, secondPeak] = findpeaks(flip(abs0MS(peaks(i) - 50 : peaks(i))), ...
%         'MinPeakProminence', 0.01);
%     
%     if length(secondPeak) > 1
%         secondPeak = secondPeak(end);
%         peaks(i) = peaks(i) - secondPeak;
%     elseif length(secondPeak) == 1
%         peaks(i) = peaks(i) - secondPeak;
%     end
% end
% 
% heelstrikes = arrayfun(@(x)...
%     x - find(flip(diff(abs0MS(1 : x-1))) <= 0, 1, 'first'), ...
%     peaks);
% 
% hold on
% figure
% findpeaks(abs0MS, ...
%     'MinPeakHeight',0.02, 'MinPeakDistance', 1300, 'MinPeakProminence', 0.03);

% plot(zeroMeanSig)
% plot(abs(diff(abs0MS)))
% for i = 1 : length(heelstrikes)
%     line([heelstrikes(i) heelstrikes(i)], get(gca, 'YLim'),'color','g');
% end

end

