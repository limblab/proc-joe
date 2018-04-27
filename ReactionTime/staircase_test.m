threshold = 30;
sensitivity = 10; % more = less sensitive
x = 1:1:100;
y = 1./(1+exp(-(x-threshold)/sensitivity));


numSteps = 280;
amp = zeros(numSteps,1);
detect = zeros(numSteps,1);
stepSize = [10,20]; % if detect, if not detect
n = 1;

numRight = 0;
amp(1) = 100; % max amp

for s = 1:numSteps-1
    if(rand() < 1/(1+exp(-(amp(s)-threshold)/sensitivity)))
        numRight = numRight + 1;
        if(numRight == n)
           amp(s+1) = max(1,amp(s) - stepSize(1));
           numRight = 0;
        else
            amp(s+1) = amp(s);
        end
        detect(s) = 1;
    else
        amp(s+1) = min(100,amp(s) + stepSize(2));
        detect(s) = 0;
    end
end
figure()
plot(1:numSteps,amp)
mean_detect = zeros(100,1);
num_points = zeros(100,1);
for a = 1:1:100
    if(~isnan(mean(detect(amp==a))))
        mean_detect(a) = mean(detect(amp==a));
        num_points(a) = sum(amp==a);
    end
end
figure()
subplot(2,1,1)
plot(1:1:100,mean_detect)
hold on
plot(x,y)
subplot(2,1,2)
plot(1:1:100,num_points)