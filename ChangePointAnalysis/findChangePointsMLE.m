function [ changePointIdx,L,L0 ] = findChangePointsMLE(t,data)
changePointIdx = [];
numSamples = numel(data);
thetaNoChange = mean(data)*numSamples;
L0 = computePoissonLikelihood(t,data,thetaNoChange,thetaNoChange,t(floor(numSamples/2)),floor(numSamples/2));
[~,~,~,L] = computePoissonEstimatorsAndLikelihoods(t,data);
LRatio = -2*(log(L0)-log(L));
[LRatioMax,maxIdx] = max(LRatio);
LRatio(isnan(LRatio) | LRatio == 0) = 0;
LRatio(LRatio < 0) = 0;
LRatio = LRatio/sum(LRatio);
LRatio = LRatio - mean(LRatio);
changePointIdx = maxIdx;
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

