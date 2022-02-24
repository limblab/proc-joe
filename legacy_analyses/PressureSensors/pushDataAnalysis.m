%% load in all data
fileList = dir('*incremental*');
data = zeros(1000,numel(fileList));

for f = 1:numel(fileList)
    dataTemp = readtable(fileList(f).name);
    data(1:numel(dataTemp),f) = dataTemp{:,1};
end
%% plot data
plot(data(:,2))

%% we want the peak data - there are 5 presses at 3 magnitues, find those peaks
peakData = zeros(15,numel(fileList));

for f = 1:numel(fileList)
    baseline = 1035;%mean(data(1:10,f))+10;
    baselineCrosses = find(data(:,f) > baseline);
    pushGroups = find(diff(baselineCrosses) > 3);
    for pd = 1:15
        if(pd == 1)
            peakData(pd,f) = max(data(baselineCrosses(1):baselineCrosses(pushGroups(pd)),f));
        elseif(pd == 15)
            peakData(pd,f) = max(data(baselineCrosses(pushGroups(pd-1)):baselineCrosses(end),f));
        else
            peakData(pd,f) = max(data(baselineCrosses(pushGroups(pd-1)):baselineCrosses(pushGroups(pd)),f));
        end
    end
end


%% mean and variance for each mineral oil amount for each push condition

means = zeros(3,numel(fileList));
vars = zeros(3,numel(fileList));
for f = 1:numel(fileList)
    for push = 1:3
        means(push,f) = mean(peakData(1+(push-1)*5:5+(push-1)*5,f));
        vars(push,f) = var(peakData(1+(push-1)*5:5+(push-1)*5,f));
    end
end


colors = {'r','b',[0,0.5,0],'k'};
b = bar(means - 990); % subtract away baseline
for f = 1:numel(fileList)
    b(f).FaceColor = colors{f};
end
f = gcf;
f.Children.XTickLabel = {'light','medium','heavy'};
ylabel('Sensor output from baseline')
formatForLee(f)
f.Children.FontSize = 16;

%% make an errorbar plot
figure();
xOff = [-0.15,-0.075,0.075,0.15];
for f = 1:numel(fileList)
    for p = 1:3
        errorbar(xOff(f)+p,means(p,f),sqrt(vars(p,f)),'color',colors{f});
        hold on
    end
end
