function [ changePoints ] = findChangePointsMLE(t,data,alpha)
numSamples = numel(data);
thetaNoChange = mean(data)*numSamples;

% compute L0 -- model with no changepoints
L0 = computePoissonLikelihood(t,data,thetaNoChange,thetaNoChange,t(floor(numSamples/2)),floor(numSamples/2));

% compute L -- likelihoods 
[~,~,~,L] = computePoissonEstimatorsAndLikelihoods(t,data);

% find best L
[Lmax,maxIdx] = max(L);
% determine test statistic
D = 2*(log(Lmax) - log(L0));

% degrees of freedom = dalt - dnull = number of change points in each
dof = 1-0; % for this one, its easy

% chi-squared test for significance
pVal = 1-chi2cdf(D,dof);

if(pVal < alpha) % if we pass the test
    % update changePoints
    changePoints = maxIdx;

    % check for more splits by doing binary splitting?
    
end
% 
% [LRatioMax,maxIdx] = max(LRatio);
% LRatio(isnan(LRatio) | LRatio == 0) = 0;
% LRatio(LRatio < 0) = 0;
% LRatio = LRatio/sum(LRatio);
% LRatio = LRatio - mean(LRatio);
% changePoints = maxIdx;
% if(LRatioMax > 0.005)
%     % set changePointIdx
%     changePointIdx = maxIdx;
%     % check for split
%     tLow = t(1:changePointIdx-1,1);
%     tHigh = t(changePointIdx+1:end,1);
%     dataLow = data(1:changePointIdx-1,1);
%     dataHigh = data(changePointIdx+1:end,1);
%     
%     outHigh = findChangePointsMLE(tHigh,dataHigh);
%     outLow = findChangePointsMLE(tLow,dataLow);
%     if(~isempty(outHigh))
%         changePointIdx(end+1:end+numel(outHigh))=outHigh;
%     end
%     if(~isempty(outLow))
%         changePointIdx(end+1:end+numel(outLow))=outLow;
%     end
% end



end

