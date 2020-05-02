%% Show New categorical Datatype
clear variables
load patients SelfAssessedHealthStatus
%%
open SelfAssessedHealthStatus
%%
% cell array of strings, but only a few values
whos
%%
HealthStatus = categorical(SelfAssessedHealthStatus);
categories(HealthStatus)
%%
whos
%%
% can still read the values even though storage is compact
open HealthStatus
