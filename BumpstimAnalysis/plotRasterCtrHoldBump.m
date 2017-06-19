function [] = plotRasterCtrHoldBump(cds,neuronNumber,ctrHoldBumpNoStimMask,ctrHoldBumpStimMask,varargin)
% plot rasters for ctrHold bumps for the provided neuronNumber nn
colorRect = 'r';
rpLineLength = 0.33;
verticalLine = 0;
markerSize = 6;
stimsPerBump = 10;
trialLength = 0.3; % 200 ms -- 100ms of stim (10 stims * 10ms per stim)
preTime = 0.15;
postTime = 0.25;
plotArtifactData = 0;
spikesArtifact = [];
plotTitle = 1;

for i = 1:2:size(varargin,2)
    switch varargin{i}
        case 'stimsPerBump'
            stimsPerBump = varargin{i+1};
        case 'trialLength'
            times = varargin{i+1};
            preTime = times(1);
            postTime = times(2);
        case 'plotArtifactData'
            plotArtifactData = varargin{i+1};
        case 'plotTitle'
            plotTitle = varargin{i+1};
    end
end

figure('Position',[1123 324 639 608]);
if(plotTitle)
    suptitle(num2str(neuronNumber));
end
% plot a set of zoomed in rasters
for cond = 1:2 % 1 = no stim, bump, 2 = stim bump
    subplot(2,1,cond)
    hold on

    if(cond==1)
        trialTable = cds.trials(ctrHoldBumpNoStimMask,:);
        yLblStr = 'Bump, No Stim Trials';
    else
        trialTable = cds.trials(ctrHoldBumpStimMask,:);
        yLblStr = 'Bump, Stim Trials';
        xLblStr = 'Time after bump start (ms)';
    end

    trialOffset = 0;
    bumpDirections = unique(trialTable(:,:).bumpDir);
    for bd = 1:numel(bumpDirections)
        bumpDir = bumpDirections(bd);
        trialTableBumpDirMask = trialTable(:,:).bumpDir == bumpDir;
        trialTableBumpDir = trialTable(trialTableBumpDirMask,:);
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
                for st = 1:numel(stimOnTimes)
                    dtStim = stimOffTimes(st) - stimOnTimes(st);
                    hold on
                    pRect = [stimOnTimes(st)*1000, trial+trialOffset-0.5, dtStim, 1]; % [x y dx dy]
                    r=rectangle('Position',pRect,'FaceColor',colorRect,'EdgeColor',colorRect); 
                end
            end

            % get spike data
            spikeMask = cds.units(neuronNumber).spikes.ts  > trialStart - preTime & ...
                    cds.units(neuronNumber).spikes.ts < trialStart + postTime;
            spikesPlot = cds.units(neuronNumber).spikes.ts(spikeMask) - trialStart;
            
            if(cond == 2 && plotArtifactData == 1) % save spike data for plotting
                spikesArtifact = [spikesArtifact;spikesPlot+trialStart];
            end
            % plot spike data
            if(verticalLine)
                plot([spikesPlot*1000,spikesPlot*1000]', ...
                    [(trial+trialOffset)*ones(length(spikesPlot),1)-rpLineLength/2,...
                    (trial+trialOffset)*ones(length(spikesPlot),1)+rpLineLength/2]','k','linewidth',1.5)      
            else
                plot(spikesPlot*1000,(trial+trialOffset)*ones(length(spikesPlot),1),'k.','markerSize',markerSize)
            end            
        end
        % plot horizontal line to denote bump directions
        trialOffset = trialOffset + size(trialTableBumpDir,1) + 1;
        if(bd ~= numel(bumpDirections))
            plot([-preTime*1500,postTime*1500],[trialOffset, trialOffset],'b','linewidth',1.5)
            trialOffset = trialOffset + 1;
        end
    end
    
    ylim([-3,trialOffset]+rpLineLength/2)
    xlim([-preTime*1000,postTime*1000])
    ylabel(yLblStr)
    if(cond == 2)
        xlabel(xLblStr);
    end
    formatForLee(gcf);
    ax=gca;
    ax.YTick = [0;trialOffset];
    ax.YMinorTick = 'off';
    ax.YTickLabel = {num2str(0),num2str(size(trialTable,1))};
end

% plot raw waveforms around the artifact data
if(plotArtifactData == 1)
    figure();
    rowSubplot = 1;
    colSubplot = 1;
    numSubplots = rowSubplot*colSubplot;
    wavesPerSubplot = 50;%ceil(size(spikesArtifact,1)/numSubplots);
    xData = (0:1:numel(cds.rawData.waveforms(1,:))-1)/30000*1000;
    subplotNum = 1;
    for waveTs = 1:wavesPerSubplot
        subplot(rowSubplot,colSubplot,subplotNum);
        hold on
        waveIdx = find(cds.units(neuronNumber).spikes.ts == spikesArtifact(waveTs));
        rawIdx = getRawDataIdx(cds.units(neuronNumber).spikes.ts(waveIdx),cds.units(neuronNumber).chan,...
                cds.rawData.ts,cds.rawData.elec);
        if(rawIdx ~= -1)
            plot(xData,cds.rawData.waveforms(rawIdx,:))
            ylim([-500,500])
        end    
        if(mod(waveTs,wavesPerSubplot) == 0)
            subplotNum = subplotNum+1;
        end
        
    end
    
    if(rowSubplot*colSubplot == 1)
        xlabel('Time (ms)')
        ylabel('Voltage (\muV)')
        formatForLee(gcf);
    end
end

end

