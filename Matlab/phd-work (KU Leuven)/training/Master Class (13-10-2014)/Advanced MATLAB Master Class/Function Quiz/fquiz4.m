function y = fquiz4(x)

% Copyright 2007 The MathWorks, Inc.

%% What kind of function is MAKELINE?
a = 3;
b = 5;
y = makeline(x,a,b);


function y = makeline(x,a,b)
% MAKELINE  Make a line
y = a*x + b;
