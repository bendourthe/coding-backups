function testInPlace

% Copyright 2007 The MathWorks, Inc.

%% Show in-place function calling

%% Make some data
% n=50*2^20;
% x=randn(n,1); % 400MB variable (3G switch)
% n = 35*2^20;
% x = randn(n,1); % 280MB variable (1G RAM and 1.5G swap + fineprint driver)
% n = 48*2^20;  % 3G switch, 1G RAM, fineprint
n = 38*2^20;
x = randn(n,1);
disp 'create initial var'
pause(5)

%% Call regular function with same LHS
x = myfunc(x); 
disp 'not in place but no new var'
pause

%% Call inplace function with same LHS
x = myfuncInPlace(x);
disp 'in place'
pause

% %% Call regular function with different LHS
% y = myfunc(x); 
% disp 'not in place and new var'
% pause(4)
% 
% %% Call inplace function with same LHS
% % change to new LHS and get an error
% x = myfuncInPlace(x);
% disp 'in place'
% pause(4)
% disp 'finished'
% 
% 
