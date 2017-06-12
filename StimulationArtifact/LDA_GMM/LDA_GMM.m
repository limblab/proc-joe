function [ Y ] = LDA_GMM(X, K, d)
% X -- m x n matrix or spikes (m samples per spike)
epsilonWhitening = 0.0001;
threshold = 0.001;
counter = 1;
vold = 0;
vdiff = threshold + 1;
%% initialize W -- m x d matrix 
if d >= K
    d = K;
end

m = size(X,1);
numIters = 10;
R = floor(size(X,2)/2); % use half of the spikes
idx = 1:1:size(X,2);
bestW = [];
vBestW = -10000;
for i = 1:numIters
    Xmask=datasample(idx,R,'Replace',false);
    S = X(:,Xmask);
    if mod(i,2)==0 % even
        W0 = rand(m,d);
    else % odd
        W0 = pca(S');
        W0 = W0(:,1:d);
    end
    [W,L] = LDA_KM(S,W0,d,epsilonWhitening);
    % build M
    M = [];
    for k = 1:K
        clusterMask = L==k;
        clusterData = S(:,L);
        M(:,k) = mean(clusterData,2);
    end
    % reformat L
    newL = zeros(numel(L),K);
    for k = 1:K
        newLmask = L==k;
        newL(:,k) = newLmask;
    end
    L=newL;
    Sb = M*L'*L*M';
    Sw = (S-M*L')*(S-M*L')';
    v = trace(W'*Sb*W)/trace(W'*Sw*W);
    % pick best W based on v -- trace ratio
    if i == 1 || v > vBestW
        bestW = W;
        vBestW = v;
    end
end

while vdiff > threshold && counter < 20
    %% Y = W'X
    Y = W'*X;

    %% whiten Y
    Y = whiten(Y,epsilonWhitening);

    %% do GMM with K+1 components to obtain L
    GMModel = fitgmdist(Y',K+1,'RegularizationValue',0.1);
    L = cluster(GMModel,Y');
    %% remove outliers from L
    for k = 1:K
        if sum(L==k) < 10
            Lmask = L==k;
            L(Lmask) = -1;
        end
    end
    Xmask = L~=-1;
    Xuse = X(:,Xmask);
    L = L(Xmask);
    %% do LDA to update W -- needs to be an m x d
    W = LDA(Xuse',L)';
    % remove constants and add necessary columns
    W = W(2:end,:);
    if size(W,2) ~= d
        numAdd = d-size(W,2);
        W(:,end+1:end+numAdd) = 0;
    end
    %% calculate new v
    % build M
    M = [];
    for k = 1:K
        clusterMask = L==k;
        clusterData = Xuse(:,L);
        M(:,k) = mean(clusterData,2);
    end
    % reformat L
    newL = zeros(numel(L),K);
    for k = 1:K
        newLmask = L==k;
        newL(:,k) = newLmask;
    end
    L=newL;
    Sb = M*L'*L*M';
    Sw = (Xuse-M*L')*(Xuse-M*L')';
    vnew = trace(W'*Sb*W)/trace(W'*Sw*W);
    %% if change in v is small, end
    if counter == 1 
        vold = vnew;
        vdiff = threshold+1;
    else
        vdiff = abs(vnew-vold);
        vold = vnew;
    end

    counter = counter+1;
end

end

function [W,L] = LDA_KM(X,W0,K,epsilon)
W = W0; L = []; 
i=0;
imax = 5;
while i<imax
    % project
    Y=W'*X;
    Lnew = L;
    % whiten
    Y = whiten(Y,epsilon);
    % K-means
    [L] = kmeans(Y',K);
    % LDA using L to get W
    W = LDA(X',L)';
    % remove constants
    W = W(2:end,:);
    % repeat
    i=i+1;
end

end