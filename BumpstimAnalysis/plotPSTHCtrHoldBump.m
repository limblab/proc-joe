function [] = plotPSTHCtrHoldBump(cds,neuronNumber,ctrHoldBumpNoStimMask,ctrHoldBumpStimMask,varargin)
% plot rasters for ctrHold bumps for the provided neuronNumber nn
bumpDirectionColors = {'m','r','b',[0,0.5,0]};
rpLineLength = 0.33;
verticalLine = 0;
markerSize = 6;
stimsPerBump = 10;
trialLength = 0.3; % 200 ms -- 100ms of stim (10 stims * 10ms per stim)
preTime = 0.15;
postTime = 0.25;
plotTitle = 1;
binSize = 0.005; % 5ms bins
gaussianSmooth = 0;
gaussianSmooth_std = binSize;
makeFigure = 1;
makeLegend = 1;
makeAxisLabels = 1;
for i = 1:2:size(varargin,2)
    switch varargin{i}
        case 'stimsPerBump'
            stimsPerBump = varargin{i+1};
        case 'trialLength'
            times = varargin{i+1};
            preTime = times(1);
            postTime = times(2);
        case 'plotTitle'
            plotTitle = varargin{i+1};
        case 'gaussianSmooth'
            gaussianSmooth = varargin{i+1};
        case 'gaussianStd'
            gaussianSmooth_std = varargin{i+1};
        case 'binSize'
            binSize = varargin{i+1};
            gaussianSmooth_std = binSize;
        case 'makeFigure'
            makeFigure = varargin{i+1};
        case 'makeLegend'
            makeLegend = varargin{i+1};
        case 'makeAxisLabels'
            makeAxisLabels = varargin{i+1};
    end
end

if(makeFigure)
    close all
    figure(1);
    figure(2);
    figure(3);
    figure(4);
end
hold on
if(plotTitle)
    suptitle(num2str(neuronNumber));
end


for cond = 1:2 % 1 = no stim, bump, 2 = stim bump
    if(cond==1)
        trialTable = cds.trials(ctrHoldBumpNoStimMask,:);
    else
        trialTable = cds.trials(ctrHoldBumpStimMask,:);
    end
    yLblStr = 'Firing rate (Hz)';
    xLblStr = 'Time after bump start (ms)';
    
    trialOffset = 0;
    bumpDirections = unique(trialTable(:,:).bumpDir);
    for bd = 1:numel(bumpDirections)
        spikesPlot = [];
        figure(bd)
        hold on
        bumpDir = bumpDirections(bd);
        trialTableBumpDirMask = trialTable(:,:).bumpDir == bumpDir;
        trialTableBumpDir = trialTable(trialTableBumpDirMask,:);
%         trialTableBumpDir = trialTable;
        for trial = 1:size(trialTableBumpDir,1) % for each trial
            % plot trial related data
            if(cond == 1)
                trialStart = trialTableBumpDir(trial,:).bumpTime;
            elseif(cond == 2) % stim data
                % find set of 10 times in cds.stimOn
                stimStart = cds.stimOn(1:stimsPerBump:end);
                stimOffset = stimStart - trialTableBumpDir(trial,:).startTime;
                stimOffset(stimOffset < 0) = 10000;
                [~,stimStartIdx] = min(stimOffset);
                stimOnIdx = (stimStartIdx-1)*stimsPerBump + 1;
                stimOnTimes = cds.stimOn(stimOnIdx:stimOnIdx+stimsPerBump-1) - cds.stimOn(stimOnIdx);
                trialStart = cds.stimOn(stimOnIdx);
                stimOffTimes = stimOnTimes + 1/1000; % blank for 1ms
            end

            % get spike data
            spikeMask = cds.units(neuronNumber).spikes.ts  > trialStart - preTime & ...
                    cds.units(neuronNumber).spikes.ts < trialStart + postTime;
            spikesPlot = [spikesPlot;cds.units(neuronNumber).spikes.ts(spikeMask)-trialStart];
                      
        end
        
        
        % bin spike data
        bE = -preTime:binSize:postTime;
        [bC,bE] = histcounts(spikesPlot,bE);
        bE = (bE(1:end-1) + (bE(2)-bE(1))/2)*1000;
        bC = bC/binSize/size(trialTableBumpDir,1);
        
        if(gaussianSmooth)
            kernel_width = ceil(3*gaussianSmooth_std/binSize);
            kernel = normpdf(-kernel_width*binSize: ...
                binSize: ...
                kernel_width*binSize,...
                0, gaussianSmooth_std); 
            normalizer = conv(kernel,ones(1,length(bE)));
            smoothed_fr_inter = conv(kernel,bC)./normalizer;
            smoothed_fr = smoothed_fr_inter(kernel_width+1:end-kernel_width);
            bC = smoothed_fr;
        end
        
        if(cond==1) % no stim
            linestyle = '-';
            c='k';
        else % stim
            linestyle = '-';
            c='r';
        end
        
        if(cond==1)
            plot(bE,bC,'linestyle',linestyle,'color','k','linewidth',2)
        else
            plot(bE,bC,'linestyle',linestyle,'color',bumpDirectionColors{bd},'linewidth',2)
        end
        xlim([-preTime*1000,postTime*1000])   
        
        formatForLee(gcf);
        if(makeAxisLabels)
            if(bumpDir == 180)
                ylabel(yLblStr,'fontsize',22)
            elseif(bumpDir == 90)
                xlabel(xLblStr,'fontsize',22)
            end
        end
    end

end

% make ylims the same for all 4 plots
yLimMax = 0;
for fig = 1:4
    f = figure(fig);
    if(f.Children(1).YLim(2) > yLimMax)
        yLimMax = f.Children(1).YLim(2);
    end
end

for fig = 1:4
    f = figure(fig);
    f.Children(1).YLim = [0,yLimMax];
end

if(makeLegend)
    l=legend('Bump, No Stim','Bump and Stim');
    set(l,'box','off','location','northwest');
end


end

