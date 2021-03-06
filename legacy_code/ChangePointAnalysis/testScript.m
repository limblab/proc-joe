%% generate fake data from two poisson distributions 
clear;

T = 100;
numSamples = 110;
binSize = T/numSamples;
changePoint = 50;
lambda1 = 3;
lambda2 = 1;
lambda3 = 1;
t = (0:binSize:T-binSize)';
data = zeros(numSamples,1);
data(1:changePoint,1) = poissrnd(lambda1,size(data(1:changePoint,1),1),1);
data(changePoint+1:end,1) = poissrnd(lambda2,size(data(changePoint+1:end,1),1),1);
% data(changePoint+101:end,1) = poissrnd(lambda3,size(data(changePoint+101:end,1),1),1);

plot(t,data)
% do binary splitting to get change points
alpha = 0.05;
changePointIdx = findChangePointsMLE(t,data,alpha);
% figure;
% findchangepts(data,'MaxNumChanges',2)

%%
load('SpindleStimData_10ms.mat')

t = bE(1:floor(numel(bE)/2));
data = floor(bC(1:floor(numel(bE)/2)));
% findchangepts(data,'Statistic','mean','MaxNumChanges',2)
[changePoints,confInter] = findChangePointsMLE(t',data',0.05);

%% 






