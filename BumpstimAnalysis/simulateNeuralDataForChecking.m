fr = 20; % Hz
dt = 1/10000; % s
nBins = 1000; % 100 ms spike train
numTrials = 2300;
train = rand(numTrials, nBins) < fr*dt;

% raster sorted by first stim time
x = [];
y = [];
for i = 1:numTrials
    if(isempty(x))
        x = find(train(i,:)==1);
        y = i*ones(1,numel(find(train(i,:)==1)));
    else
        x = [x, find(train(i,:)==1)];
        y = [y, i*ones(1,numel(find(train(i,:)==1)))];
    end
end
%
optsPlot = [];
optsPlot.SORT_DATA = 'postStimuliTime';
optsSave = [];
plotRasterLIB(x,y,optsPlot,optsSave);
