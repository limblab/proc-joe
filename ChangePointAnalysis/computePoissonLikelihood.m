function [ L ] = computePoissonLikelihood( t,data,theta0,theta1,tau,tauIdx)

binSize = t(2)-t(1);
pTau = (tau-t(tauIdx))/binSize;
theta0 = theta0*binSize; % has to do with the factorial bugging out later
theta1 = theta1*binSize;

% L = -tau*theta0 + sum(data(1:tauIdx-1,1))*log(theta0) ...
%         - (t(end)-(tau+1))*theta1 + sum(data(tauIdx+1:end,1))*log(theta1) ...
%         - pTau*theta0-(1-pTau)*theta1+data(tauIdx,1)*log(pTau*theta0+(1-pTau)*theta1) ...
%         - 0;%sum(log(factorial(floor(data))));

L = -(t(tauIdx))*theta0 + sum(data(1:tauIdx-1,1))*log(theta0) ...
        - (t(end)-t(tauIdx+1))*theta1 + sum(data(tauIdx+1:end,1))*log(theta1) ...
        - pTau*theta0-(1-pTau)*theta1+data(tauIdx,1)*log(pTau*theta0+(1-pTau)*theta1) ...
        - sumLogFactorial(data);

end

function [out] = sumLogFactorial(data)

out = 0;
for i = 1:numel(data)
    out = out + sum(log(1:floor(data(i))));
end

end