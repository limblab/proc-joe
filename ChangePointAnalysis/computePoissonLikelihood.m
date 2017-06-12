function [ L ] = computePoissonLikelihood( t,data,theta0,theta1,tau,tauIdx)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

pTau = t(tauIdx+1)-tau;
L = -tau*theta0 + sum(data(1:tauIdx-1,1))*log(theta0) ...
        - (t(end)-(tau+1))*theta1 + sum(data(tauIdx+1:end,1))*log(theta1) ...
        - pTau*theta0-(1-pTau)*theta1+data(tauIdx,1)*log(pTau*theta0+(1-pTau)*theta1) ...
        - 0;%log(sum(factorial(data)));

end

