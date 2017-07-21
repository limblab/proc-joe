
%% step 0 -- take two pulse trains and add them together?
N = 2; % num sources
M = 5; % num electrodes
L = 10; % length of a source pulses
K = 100; % large
numSamples = 20000;
source1Wave = 10*ones(L,1);
source2Wave = 20*ones(L,1);
noise = 4;
    
signal = noise*rand(numSamples,M) - noise/2;
for i = 1:M
    mix = rand(2,1)*10;
    source1 = 0;
    source2 = 0;
    source1Idx = 1;
    source2Idx = 1;
    for s = 1:numSamples
        if(source1)
            if(source1Idx == L)
                source1 = 0;
            end
            signal(s,i) = signal(s,i) + mix(1)*source1Wave(source1Idx);
            source1Idx = source1Idx + 1;
        elseif(~source1 && rand(1) > 0.999)
            source1 = 1;
            source1Idx = 1;
        end
        if(source2)
            if(source2Idx == L)
                source2 = 0;
            end
            signal(s,i) = signal(s,i) + mix(2)*source2Wave(source2Idx);
            source2Idx = source2Idx + 1;
        elseif(~source2 && rand(1) > 0.999)
            source2 = 1;
            source2Idx = 1;
        end
    end
    
end

x = signal';

%% step 1
gamma = zeros(size(x,2),1);
for n = 1:size(x,2)
    xbar = x(:,n);
    xbar = repmat(xbar,1,K);
    xbar = reshape(xbar',1,numel(xbar))';
    Rx = sum(xbar*xbar');
%     covXbar = xbar*xbar'; % this line is wrong. should be a numSamp x numSamp matrix
    gamma(n,1) = (xbar)'*inv(covXbar)*(xbar);
end


% noiseVar = % smallest eigenvalues
% noiseThreshold = noiseVar*norm(inv(cov(gamma)),1);