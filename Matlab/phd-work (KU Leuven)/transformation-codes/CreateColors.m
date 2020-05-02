function [colors] = CreateColors(N,S,V)
%Creates N colors equally spaced in hue
if nargin < 3
    V = 0.75;
    if nargin < 2
        S = 0.75;
    end
end
HSV = zeros(N,3);
HSV(:,1) = (1:N)'/N;
HSV(:,2) = S*ones(N,1);
HSV(:,3) = V*ones(N,1);
colors = hsv2rgb(HSV);
end

