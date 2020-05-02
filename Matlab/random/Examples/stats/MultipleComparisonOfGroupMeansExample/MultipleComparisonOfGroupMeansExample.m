%% Multiple Comparison of Group Means
%%
% Load the sample data.

% Copyright 2015 The MathWorks, Inc.

load carsmall
%%
% Perform a one-way analysis of variance (ANOVA) to see if there is any
% difference between the mileage of the cars by origin.
[p,t,stats] = anova1(MPG,Origin,'off');
%%
% Perform a multiple comparison of the group means.
[c,m,h,nms] = multcompare(stats);
%%
% |multcompare| displays the estimates with comparison intervals around
% them. You can click the graphs of each country to compare its mean to
% those of other countries.
%%
% Now display the mean estimates and the standard errors with the
% corresponding group names.
[nms num2cell(m)]


