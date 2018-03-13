%% messing around with modelling the artifact
load('artifact_lowGain')

figure; plot(squeeze(artifact(18,1,1:200))')
a=squeeze(artifact(18,1,1:200))';

% optimize parameters to fit all artifacts -- use sum of max xcorr as
% metric. we want maximum
numTrials = 1000;
scores = zeros(numTrials,1);
downsampleFactor = 333;

R = zeros(numTrials,3);
C = zeros(numTrials,2);

R(1,:) = [10,10,10]*1000;
C(1,:) = [0.000001,0.001];  

stimParams = [0.2,0.053,1]; % in ms
dt = 0.0001; % in ms
T = 2/3*10; % in ms
for i = 1:numTrials % 100 steps
    [~,uc,~] = runArtifactModel(R(i,:), C(i,:), dt, T, stimParams);
    uc = downsample(uc,downsampleFactor);
    for j = 1:size(a,1)  
        % match in time and amplitude
        [r,lag] = xcorr(uc(1:200),a(j,:));
        [~,lagIdx] = max(r);
        uc = circshift(uc(1:200),-1*lag(lagIdx));
        uc = uc.*(uc(:)>0)./max(uc)*abs(max(a(j,:))) + uc.*(uc(:)<0)./abs(min(uc))*abs(min(a(j,:)));
        % compute SSE
        scores(i,1) = scores(i,1) + sum((uc(104:118,1)-a(j,104:118)').^2);
    end
    
    % update parameters R(i+1) = ... C(i+1) = ....
    % randomly select a parameter, change that based on gradient (run the
    % simulation again)
    R(i+1,:) = R(1,:) + rand(1,3)*100 - 50;
    C(i+1,:) = C(1,:) + rand(1,2).*[0.00000001 0.000001] - [0.00000001 0.000001]/2;
    paramRC = ceil(rand()*2);
    
    
end

[~,idx] = min(scores);
[t,uc,ui] = runArtifactModel(R(idx,:), C(idx,:), dt, T, stimParams);
t = downsample(t,downsampleFactor);
uc = downsample(uc,downsampleFactor);
[r,lag] = xcorr(uc(1:200),a(j,:));
[~,lagIdx] = max(r);
uc = circshift(uc(1:200),-1*lag(lagIdx));
uc = uc.*(uc(:)>0)./max(uc)*abs(max(a(j,:))) + uc.*(uc(:)<0)./abs(min(uc))*abs(min(a(j,:)));

idx=1;
hold on
plot(uc);

figure;
plot(uc-a')
