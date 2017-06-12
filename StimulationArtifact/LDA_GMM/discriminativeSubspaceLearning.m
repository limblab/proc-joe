%% Noise-robust unsupervised spike sorting based on discriminative subspace learning with
% outlier handling, Mohammad Reza Keshtkaran and Zhi Yang

% build X -- data of spike waveforms
chan =1;
wf=nevData.waveforms;
wfMask = nevData.elec == chan;
wf = wf(wfMask,:);

X = wf';

% init K, K', d, counter;
K = 1; Kprime = K;
d = 20;
counter = 0;
while Kprime >= K && counter < 6
    K = K+1;
    Kprime = K;
    % do LDA-GMM on X to remove outliers and get Y
    Y = LDA_GMM(X,K,d);
    % do k-means on Y to obtain cluster labels L
    [L] = kmeans(Y',K);
    % build M -- matrix of means
    M = [];
    for k = 1:K
        clusterMask = L==k;
        clusterData = Y(:,L);
        M(:,k) = mean(clusterData,2);
    end
    
    % check if stopping
    for iter1 = 1:K
        for iter2 = 1:iter1-1
            X1 = Y(:,L==iter1);
            X2 = Y(:,L==iter2);
            mean1 = M(:,iter1);
            mean2 = M(:,iter2);
            X1proj = (mean1-mean2)'*X1;
            X2proj = (mean1-mean2)'*X2;
            X12proj = [X1proj,X2proj];
            [~,~,adstat1,~] = adtest(X1proj);
            [~,~,adstat2,~] = adtest(X2proj);
            [~,~,adstat12,~] = adtest(X12proj);
            if adstat1 > adstat12 || adstat2 > adstat12 % stop
                Kprime = K-1;
            end
        end
    end
    counter = counter+1;
end
