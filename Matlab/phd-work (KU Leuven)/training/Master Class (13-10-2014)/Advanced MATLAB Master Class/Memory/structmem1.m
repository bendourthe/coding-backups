%% Example: Managing Structure Memory 

% Copyright 2007 The MathWorks, Inc.

clear all, close all

%% Create a structure
% (3000*3000 elts) * (8 bytes/elt) = 72,000,000 B ~= 70MB
% Could be interesting to bring up task manager here
s.A = rand(3000,3000);

%%
s.B = rand(3000,3000);

%% Copy the structure
% What will happen to memory?
% Watch task manager here
sNew = s;       


%% Modify the structure
% Can you explain what happened to memory?
s.A(1,1) = 17;   
