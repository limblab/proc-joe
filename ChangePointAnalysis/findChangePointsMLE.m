function [ changePoint, confInter] = findChangePointsMLE(t,data,alpha)
% based on West and Ogden 1997

changePoint = [];
confInter = [];

numSamples = numel(data);
binSize = t(2)-t(1);
thetaNoChange = mean(data)/binSize;

% compute L0 -- model with no changepoints
L0 = computePoissonLikelihood(t,data,thetaNoChange,thetaNoChange,t(floor(numSamples/2)),floor(numSamples/2));

% compute L -- likelihoods 
[theta0,theta1,tau,L] = computePoissonEstimatorsAndLikelihoods(t,data);

% find best L
[Lmax,maxIdx] = max(L);
% determine test statistic
D = 2*(Lmax - L0);

% degrees of freedom = dalt - dnull = number of change points in each
dof = 3-1; % for this one, its easy?

% chi-squared test for significance
pVal = 1-chi2cdf(D,dof);

if(pVal <= alpha/numSamples) % if we pass the test
    % update changePoints
    changePoint = maxIdx;
    % do conf interval for this by finding highest density region
    Lcut = max(L); % start from highes point
    stepSize = Lcut*0.0001;
    while(computeDensity(L,Lcut) < alpha)
        Lcut = Lcut - stepSize;  
    end
    Lcut = Lcut + stepSize;
    % get confidence interval from Lcut 
    confInter = getConfidenceIntervalCPA(L,Lcut);
end



end

