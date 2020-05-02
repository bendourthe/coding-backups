clear all; close all; clc

cd 'C:\Users\bdour\OneDrive\Work\Calgary\Projects\Dr. Scholls\Results\';
pathName = 'C:\Users\bdour\OneDrive\Work\Calgary\Projects\Dr. Scholls\Results';

data = xlsread('Dr.Scholls_Raw_Data.xlsx',1);

Condition = data(1,2:6);
Subject = data(2:end,1);
LeftAI = data(2:end,2:6);
RightAI = data(2:end,7:11);
MaxT = data(2:end,12:16);
PercVC = data(2:end,17:21);
P2P = data(2:end,22:26);
MeanT = data(2:end,27:31);
MeanA = data(2:end,32:36);
MeanF = data(2:end,37:41);

[pLAI,~,statsLAI] = anova1(LeftAI,Condition,'off');
[resultsLAI,meansLAI] = multcompare(statsLAI,'CType','bonferroni')

[pRAI,~,statsRAI] = anova1(RightAI,Condition,'off');
[resultsRAI,meansRAI] = multcompare(statsRAI,'CType','bonferroni')

[pMT,~,statsMT] = anova1(MaxT,Condition,'off');
[resultsMT,meansMT] = multcompare(statsMT,'CType','bonferroni')

[pPVC,~,statsPVC] = anova1(PercVC,Condition,'off');
[resultsPVC,meansPVC] = multcompare(statsPVC,'CType','bonferroni')

[pP2P,~,statsP2P] = anova1(P2P,Condition,'off');
[resultsP2P,meansP2P] = multcompare(statsP2P,'CType','bonferroni')

[pMeT,~,statsMeT] = anova1(MeanT,Condition,'off');
[resultsMeT,meansMeT] = multcompare(statsMeT,'CType','bonferroni')

[pMA,~,statsMA] = anova1(MeanA,Condition,'off');
[resultsMA,meansMA] = multcompare(statsMA,'CType','bonferroni')

[pMF,~,statsMF] = anova1(MeanF,Condition,'off');
[resultsMF,meansMF] = multcompare(statsMF,'CType','bonferroni')