threshold = 30;
sensitivity = 10; % more = less sensitive
x = 1:1:100;
y = 1./(1+exp(-(x-threshold)/sensitivity));


numSteps = 300;
stepReset = 25;
amp = zeros(numSteps,1);
detect = zeros(numSteps,1);
stepSize = [5,5]; % if detect, if not detect
n = 2;

step = 1;

numRight = 0;
amp(1) = 100; % max amp

for s = 1:numSteps-1
    if(mod(step,stepReset) == 0)
        amp(s+1) = amp(1); % reset amp to max amp
        numRight = 0;
        detect(s) = rand() < 1/(1+exp(-(amp(s)-threshold)/sensitivity));
        step = 1;
    elseif(rand() < 1/(1+exp(-(amp(s)-threshold)/sensitivity)))
        numRight = numRight + 1;
        if(numRight == n)
           amp(s+1) = max(0,amp(s) - stepSize(1));
           numRight = 0;
           step = step + 1;
        else
            amp(s+1) = amp(s);
        end
        detect(s) = 1;
    else
        amp(s+1) = min(100,amp(s) + stepSize(2));
        detect(s) = 0;
        step = step + 1;
    end
end

mean_detect = zeros(100,1);
num_points = zeros(100,1);
for a = 1:1:100
    if(~isnan(mean(detect(amp==a))))
        mean_detect(a) = mean(detect(amp==a));
        num_points(a) = sum(amp==a);
    end
end
% figure()
% subplot(3,1,1)
% plot(1:1:100,mean_detect)
% hold on
% plot(x,y)
% subplot(3,1,2)
% plot(1:1:100,num_points)
% subplot(3,1,3)
% plot(1:numel(amp),amp)

figure()
plot(y(num_points~=0),num_points(num_points~=0),'.','markersize',12)