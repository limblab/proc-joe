%% this assumes a cds has been loaded by analyze3000PulseData

%% get bin counts and bin edges across array
% make opts struct for future use
opts.MAP_FILENAME = mapFileName;
opts.PRE_TIME = 10/1000;
opts.POST_TIME = 90/1000;
opts.BIN_SIZE = 0.0002;
opts.WAVEFORM_IDX = [1:numel(unique(cds.waveforms.waveSent))];
opts.CHANNEL_IDX = [1:numel(unique(cds.waveforms.chanSent))];
opts.BIN_DATA = 1;
opts.SPIKE_DATA = 1;

[arrayData] = getCountsAcrossArray(cds,opts);

%% get probability of response for all units

opts.WINDOW = [1,5]/1000;
opts.AUTOMATIC_WINDOW = 0; % not implemented currently
opts.SUBTRACT_BASELINE = 0;

[arrayData] = getProbabilityOfResponse(arrayData,opts);


%% all of this needs to be rewritten with what changed above

%% look at the distributions (poisson, normal, etc) of the stim response

opts.WINDOW = [60,64];
opts.BASELINE_WINDOW_IDX = [];
opts.PROJECT_TO_SLOPE_1 = 0;
opts.COLOR_MARKERS_RESPONSE = 1;
opts.PLOT_FIT_LINE = 1;
opts.PLOT_LOG = 0;
opts.MIN_RATIO = 1;
opts.MAX_RATIO = 4;

opts.FIGURE_SAVE = 1;
opts.FIGURE_DIR = folderpath;
opts.FIGURE_NAME = 'Chips_20171026_distributionAnalysis';

arrayDataFits = getDistributionsOfStimResponse(arrayData,opts);


%% plot data with best fit line going through (0,0)
idx = 1;
figure
plot(arrayData.means{idx},arrayData.vars{idx},'.',...
                'markersize',markerSize)
% find best fit line going through (0,0)
[F,gof] = fit(reshape(arrayData.means{idx},1,[])',reshape(arrayData.vars{idx},1,[])','a*x');

xFit = [0,max(max(arrayData.means{idx}))];
yFit = F.a*xFit;
hold on
plot(xFit,yFit,'k')
%% plots stim and non-stim region with colors based on ratio
markerSize = 10;
colormap = jet;
colors = colormap;
maxRatio = 4;
minRatio = 0;
figure();
hold on
for i = 1:size(arrayData.means{idx},1)
    for j = 1:size(arrayData.means{idx},2)
        if(arrayData.means{1}(i,j) > 0 && arrayData.means{2}(i,j) > 0)
            % get color based on change in firing rate between idx 1 and idx 2
            r = arrayData.means{1}(i,j)/arrayData.means{2}(i,j);
            r = min(maxRatio,max(minRatio,r));
            colorIdx = floor(size(colors,1)*(r-minRatio)/(maxRatio-minRatio));
            % plot baseline point
            plot(arrayData.means{2}(i,j),arrayData.vars{2}(i,j),'.',...
                'markersize',markerSize,'color',colors(colorIdx,:))

            % plot stim point
            plot(arrayData.means{1}(i,j),arrayData.vars{1}(i,j),'.',...
                'markersize',markerSize,'color',colors(colorIdx,:))

            % plot arrow between 2
            p1=[arrayData.means{2}(i,j),arrayData.vars{2}(i,j)];
            p2=[arrayData.means{1}(i,j),arrayData.vars{1}(i,j)];
            dp = p2-p1;
            quiver(p1(1),p1(2),dp(1),dp(2),'MaxHeadSize',0.01,'color',colors(colorIdx,:))

        end
    end
end

% plot best fit line from before
xFit = [0,max(max(arrayData.means{idx}))];
yFit = F.a*xFit;
hold on
plot(xFit,yFit,'k')
% ylim([0,5.5E-6])
% xlim([0,8E-3])
%% generate spike sequences with random mean firing rates, bin and compute variance, then plot
numFr = 100;
numRuns = 3000;
dt = 2/30000;
maxTime = 60/1000;
fr = 70*rand(numFr,1)+5;
spikes = rand(numFr,numRuns,floor(maxTime/dt)) < fr*dt;

binnedSpikes = sum(spikes,2)/numRuns;
meanSpikes = mean(binnedSpikes,3);
varSpikes = var(binnedSpikes,0,3);
figure();
plot(meanSpikes,varSpikes,'.','markersize',10)