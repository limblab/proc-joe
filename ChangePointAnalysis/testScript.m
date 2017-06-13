%% generate fake data from two poisson distributions 
clear;

T = 100;
numSamples = 1000;
binSize = T/numSamples;
changePoint = 500;
lambda1 = 2;
lambda2 = 10;
lambda3 = 2;
t = (0:binSize:T-binSize)';
data = zeros(numSamples,1);
data(1:changePoint,1) = poissrnd(lambda1,size(data(1:changePoint,1),1),1);
data(changePoint+1:changePoint+100,1) = poissrnd(lambda2,size(data(changePoint+1:changePoint+100,1),1),1);
data(changePoint+101:end,1) = poissrnd(lambda3,size(data(changePoint+101:end,1),1),1);

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






