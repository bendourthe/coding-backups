
point1 = [3,5];
point2 = [6,8];
testpoint = [4,10];

% A line running through two points can be represented by y = ax + b
% The slope of the line is represented by a = (y2-y1)/(x2-x1)
a = (point2(2)-point1(2))/(point2(1)-point1(1));
% The y-cordinate where he line intersects the y-axis is represented by b 
b = point1(2) - ( a * point1(1) );

% Enter the x-coordinate of the testpoint in the parameterization of the
% line
yTestpointOnLine = a * testpoint(1) + b;
% Test if testpoint is above, on, or below the line
if( yTestpointOnLine < testpoint(2) )
    disp('ABOVE');
elseif( yTestpointOnLine == testpoint(2) )
    disp('ON');
else
    disp('BELOW');
end;