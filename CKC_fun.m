
%% step 0 -- take two pulse trains and mix them
N = 2; % num sources
M = 5; % num electrodes
L = 10; % length of a source pulses
K = 50; % large
numSamples = 2000;
source1Wave = 10*ones(L,1);
source2Wave = 20*ones(L,1);
noise = 4;
    
signal = noise*rand(numSamples,M) - noise/2;
mixMatrix = rand(M,N)*10;
source1 = 0;
source2 = 0;
source1Idx = 1;
source2Idx = 1;
for s = 1:numSamples
    if(source1)
        if(source1Idx == L)
            source1 = 0;
        end
        for i = 1:M
            signal(s,i) = signal(s,i) + mixMatrix(i,1)*source1Wave(source1Idx);
        end
        source1Idx = source1Idx + 1;
    elseif(~source1 && rand(1) > 0.999)
        source1 = 1;
        source1Idx = 1;
    end
    if(source2)
        if(source2Idx == L)
            source2 = 0;
        end
        for i = 1:M
            signal(s,i) = signal(s,i) + mixMatrix(i,2)*source2Wave(source2Idx);
        end
        source2Idx = source2Idx + 1;
    elseif(~source2 && rand(1) > 0.999)
        source2 = 1;
        source2Idx = 1;
    end
end


x = signal';

%% step 1
gamma = zeros(size(x,2),1);
n=K+1;
% for n = 1:size(x,2)
    % same observation repeated
%     xbar = x(:,n);
%     xbar = repmat(xbar,1,K);
%     xbar = reshape(xbar',1,numel(xbar))';
    % observations back in time
    if(n >= K)
        xbar = x(:,n-K+1:n);
        xbar = reshape(xbar',1,numel(xbar))';
    end
    
%     Rx = pinv(mahal(xbar,xbar)*mahal(xbar,xbar)');
%     covXbar = xbar*xbar'; % this line is wrong. should be a numSamp x numSamp matrix
%     gamma(n,1) = (xbar)'*Rx*(xbar);
% end


% noiseVar = % smallest eigenvalues
% noiseThreshold = noiseVar*norm(inv(cov(gamma)),1);