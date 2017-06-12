function [ theta0,theta1,tau,L ] = computePoissonEstimatorsAndLikelihoods(t,data)
%COMPUTEPOISSONESTIMATORS Summary of this function goes here
%   Detailed explanation goes here
t = t-t(1); % shift t so that t(1) = 0
binSize = t(2)-t(1);
numSamples = size(data,1);

theta0 = zeros(numSamples,1);
theta1 = zeros(numSamples,1);
tau = zeros(numSamples,1);
L = zeros(numSamples,1);

for i = 3:numel(t)-2
    % estimators
    theta0(i,1) = sum(data(1:i-1,1))/(t(i));
    theta1(i,1) = sum(data(i+1:end,1))/(t(end)-t(i+1)+binSize);
    tau(i,1) = t(i) + (data(i,1) - theta1(i,1))/(theta0(i,1)-theta1(i,1));
    if(tau(i,1) < t(i))
        tau(i,1) = t(i);
    elseif(tau(i,1) > t(i)+binSize)
        tau(i,1) = t(i)+binSize;
    else
        disp('hi')
    end
    % likelihood
    L(i,1) = computePoissonLikelihood(t,data,theta0(i,1),theta1(i,1),tau(i,1),i);

end
L(isnan(L) | L==0) = mean(L(~(isnan(L)|L==0)));



end

