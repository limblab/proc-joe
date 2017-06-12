%% generate fake data from two poisson distributions 
clear;

T = 1;
numSamples = 1000;
binSize = T/numSamples;
changePoint = 500;
lambda1 = 20;
lambda2 = 80;
lambda3 = 20;
t = (binSize/2:binSize:T-binSize/2)';
data = zeros(numSamples,1);
data(1:changePoint,1) = poissrnd(lambda1,size(data(1:changePoint,1),1),1);
data(changePoint+1:changePoint+20,1) = poissrnd(lambda2,size(data(changePoint+1:changePoint+20,1),1),1);
data(changePoint+21:end,1) = poissrnd(lambda3,size(data(changePoint+21:end,1),1),1);

% do binary splitting to get change points
alpha = 0.05;
changePointIdx = findChangePointsMLE(t,data,alpha);

%%
findchangepts(data,'Statistic','mean','MaxNumChanges',2)

%%
load('SpindleStimData_5ms.mat')

t = bE(1:end-1);
data = floor(bC*20);
% findchangepts(data,'Statistic','mean','MaxNumChanges',2)
changePoints = findChangePointsMLE(t',data');
[~,cIdx] = max(L);
% cIdx