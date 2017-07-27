function [] = plotPSTHStim(cds,neuronNumber,varargin)
preTime = 10/1000; % in seconds
postTime = 20/1000; % in seconds
plotAllArtifacts = 0;
spikesArtifact = [];
plotTitle = 1;
titleToPlot = '';
stimsPerBump = -1;
bumpTask = 0;
makeFigure = 1;
numWaveformTypes = 1;
waveformsSentExist = any(isfield(cds,'waveforms'));
waveformsMakeSubplots = 0;
waveformTypesPlot = 1;
if(waveformsSentExist)
    waveformTypesPlot = 1:1:numel(unique(cds.waveforms.waveSent));
end
plotLine = 0;
lineColor = 'k';
lineWidth = 2;
binSize = 0.0002;

% plot waves near artifact stuff
plotSpikeWaveforms = 0;
timeAfterStimRawNoStim = 20/1000;
timeAfterStimRawArtifact = 5/1000;

% plot artifact sample
plotArtifacts = 0;
maxArtifactsPerPlot = 5;
rowSubplotArtifact = 4;
colSubplotArtifact = 5;

stimElectrode = -1;
chansPlot = 1;
alignWaves = 1;

% save stuff
saveFigures = 0;
figDir = '';
figPrefix = '';

plotAllOnOneFigure = 0;
makeLegend = 0;
legStr = '';
plotStimOn = 0;

plotOnlyStimuliWithResponse = 0;
plotAllStimuli = 1;
plotOnlyStimuliWithoutResponse = 0;
%% deal with varagin
for i = 1:2:size(varargin,2)
    switch varargin{i}
        case 'bumpTask'
            bumpTask = varargin{i+1};
        case 'stimsPerBump'
            stimsPerBump = varargin{i+1};
        case 'preTime'
            preTime = varargin{i+1};
        case 'postTime'
            postTime = varargin{i+1};
        case 'plotSpikeWaveforms'
            plotSpikeWaveforms = varargin{i+1};
        case 'timeAfterStimRawNoStim'
            timeAfterStimRawNoStim = varargin{i+1};
        case 'timeAfterStimRawArtifact'
            timeAfterStimRawArtifact = varargin{i+1};
        case 'plotTitle'
            plotTitle = varargin{i+1};
        case 'title'
            titleToPlot = varargin{i+1};
        case 'plotArtifacts'
            plotArtifacts = varargin{i+1};
        case 'makeFigure'
            makeFigure = varargin{i+1};
        case 'makeSubplots'
            waveformsMakeSubplots = varargin{i+1};
        case 'waveformTypes'
            waveformTypesPlot = varargin{i+1};
        case 'maxArtifactsPerPlot'
            maxArtifactsPerPlot = varargin{i+1};
        case 'colSubplotArtifact'
            colSubplotArtifact = varargin{i+1};
        case 'rowSubplotArtifact'
            rowSubplotArtifact = varargin{i+1};
        case 'binSize'
            binSize = varargin{i+1};
        case 'chans'
            chansPlot = varargin{i+1};
        case 'stimElectrode'
            stimElectrode = varargin{i+1};
        case 'saveFigures'
            saveFigures = varargin{i+1};
        case 'figDir'
            figDir = varargin{i+1};
        case 'figPrefix'
            figPrefix = varargin{i+1};
        case 'plotLine'
            plotLine = varargin{i+1};
        case 'lineColor'
            lineColor = varargin{i+1};
        case 'lineWidth'
            lineWidth = varargin{i+1};
        case 'plotAllOnOneFigure'
            plotAllOnOneFigure = varargin{i+1};
        case 'makeLegend'
            makeLegend = varargin{i+1};
        case 'legendString'
            legStr = varargin{i+1};
        case 'plotStimOn'
            plotStimOn = varargin{i+1};
        case 'plotStimuliGroup'
            if(strcmpi(varargin{i+1},'responsive'))
                plotOnlyStimuliWithResponse = 1;
                plotAllStimuli = 0;
                plotOnlyStimuliWithoutResponse = 0;
            elseif(strcmpi(varargin{i+1},'nonResponsive'))
                plotOnlyStimuliWithResponse = 0;
                plotAllStimuli = 0;
                plotOnlyStimuliWithoutResponse = 1;
            else % default case
                plotOnlyStimuliWithResponse = 0;
                plotAllStimuli = 1;
                plotOnlyStimuliWithoutResponse = 0;
            end
        
    end
end

%% extract number of waveform types if applicable
if(waveformsSentExist)
    numWaveformTypes = numel(unique(cds.waveforms.waveSent));
end

if(any(isfield(cds.waveforms,'chanSent')))
    numChans = numel(unique(cds.waveforms.chanSent));
    chanList = unique(cds.waveforms.chanSent);
else
    chanList = stimElectrode;
    numChans = 1;
end

%% extract data
for c = 1:numChans
    for i = 1:numWaveformTypes
        spikeTimeData{c,i} = [];
    end
end
stimNum = zeros(numChans,numWaveformTypes);
if(bumpTask) % plot things aligned to bump times and whatnot
    % write this later :D
else % get data after stimulations
    for st = 1:numel(cds.stimOn)
        spikeMask = cds.units(neuronNumber).spikes.ts > cds.stimOn(st)-preTime & cds.units(neuronNumber).spikes.ts < cds.stimOn(st)+postTime;
        spikesPlot = (cds.units(neuronNumber).spikes.ts(spikeMask) - cds.stimOn(st));
       
        if(plotOnlyStimuliWithResponse) % check if spikesPlot has a spike within a specified amount of time -- if not, remove
            mask = spikesPlot > 0 & spikesPlot < min(5/1000,timeAfterStimRawArtifact);
            if(sum(mask) == 0)
                spikesPlot = [];
            end
        elseif(plotOnlyStimuliWithoutResponse)
            mask = spikesPlot > 0 & spikesPlot < min(5/1000,timeAfterStimRawArtifact);
            if(sum(mask) ~= 0)
                spikesPlot = [];
            end
        end
        
        if(alignWaves && ~isempty(spikesPlot))
            spikesPlot = spikesPlot + 4/30000;
        end
        numWaves = sum(spikeMask==1);
        if(waveformsSentExist)
            if(any(isfield(cds.waveforms,'chanSent')))
                chanNumber = find(unique(cds.waveforms.chanSent)==cds.waveforms.chanSent(st));
                stimNum(chanNumber,cds.waveforms.waveSent(st)) = stimNum(chanNumber,cds.waveforms.waveSent(st)) + 1;
                if(~isempty(spikesPlot))
                    spikeTimeData{chanNumber,cds.waveforms.waveSent(st)}(end+1:end+numWaves) = spikesPlot';
                end
            else
                stimNum(1,cds.waveforms.waveSent(st)) = stimNum(1,cds.waveforms.waveSent(st)) + 1;
                if(~isempty(spikesPlot))
                    spikeTimeData{1,cds.waveforms.waveSent(st)}(end+1:end+numWaves) = spikesPlot';
                end
            end
        else
            if(any(isfield(cds.waveforms,'chanSent')))
                chanNumber = find(unique(cds.waveforms.chanSent)==cds.waveforms.chanSent(st));
                stimNum(chanNumber,1) = stimNum(chanNumber,1) + 1;
                spikeTimeData{chanNumber,1}(end+1:end+numWaves) = spikesPlot';
            else
                stimNum(1,1) = stimNum(1,1) + 1;
                spikeTimeData{1,1}(end+1:end+numWaves) = spikesPlot';
            end
        end
    end
end
%% bin and plot data - stimuli data is y-axis. spike data is x-axis
maxYLim = 0;
% find y maximum y limit
for chan = chansPlot
    for fig = waveformTypesPlot
        bE = -preTime:binSize:postTime;
        [bC,~] = histcounts(spikeTimeData{chan,fig},bE);
        bC = bC/binSize/size(spikeTimeData{chan,fig},2);
        maxYLim = max(max(bC)*1.1,maxYLim);
    end
end
for chan = chansPlot
    for fig = waveformTypesPlot
        % deals with making figure
        if(makeFigure)
            if(plotAllOnOneFigure && ((numel(waveformTypesPlot)>1 && fig == 1) || (numel(waveformTypesPlot) == 1 && chan == 1)))
                figure();
            elseif(plotAllOnOneFigure && ~((numel(waveformTypesPlot)>1 && fig == 1) || (numel(waveformTypesPlot) == 1 && chan == 1)))
                % do not make figure
            elseif(waveformsSentExist && waveformsMakeSubplots && fig == 1)
                figure(); % make figure for the subplots
                subplot(numWaveformTypes,1,fig);
            elseif(waveformsSentExist && waveformsMakeSubplots)
                subplot(numWaveformTypes,1,fig); % subplots
            elseif(waveformsSentExist && ~waveformsMakeSubplots)
                figure(); % figure for each waveform type
            elseif(fig==1)
                figure(); % figure for a single waveform type
            end
        end

        % bin data
        bE = -preTime:binSize:postTime;
        [bC,bE] = histcounts(spikeTimeData{chan,fig},bE);
        bE = bE*1000;
        bC = bC/binSize/size(spikeTimeData{chan,fig},2);
        % plot actual data now
        if(plotLine)
            if(iscell(lineColor))
                lc = lineColor{chan*fig};
            else
                lc = lineColor;
            end
            plot(bE(1:end-1)+(bE(2)-bE(1))/2,bC,'color',lc,'linewidth',lineWidth);
        else
            bar(bE(1:end-1)+(bE(2)-bE(1))/2,bC)
        end
        if(plotAllOnOneFigure)
            hold on
        end
        
        if(makeLegend && ~isempty(legStr))
            l=legend(legStr);
            set(l,'box','off');
        end
        if(plotStimOn)
            hold on
            plot([0,0],[-1000,maxYLim*10],'r','linewidth',1.5)
        end
        % clean up graph
        ylim([0,maxYLim])
        xlim([-preTime*1000,postTime*1000])
        ylabel('Average Firing Rate (spikes/s)')
        if(~waveformsSentExist || ~waveformsMakeSubplots || fig == numWaveformTypes)
            xlabel('Time after stimulation onset (ms)')  
        end
        % deals with title requests
        if(plotTitle)
            if(strcmp(titleToPlot,'') == 0)
                title(titleToPlot);
            elseif(plotAllOnOneFigure)
                % no title
            elseif(numChans > 1 && waveformsSentExist)
                title(strcat('Stim Chan: ',num2str(chanList(chan)),' Wave: ',num2str(fig)));
            elseif(numChans == 1 && waveformsSentExist)
                title(strcat('Stim Chan: ',num2str(stimElectrode),' Wave: ',num2str(fig)));
            else
                
            end
        end
        
        formatForLee(gcf);
        
        if(~waveformsMakeSubplots && saveFigures && ~plotAllOnOneFigure)
            fname = strcat(figPrefix,'nn',num2str(neuronNumber),'_chan',num2str(cds.units(neuronNumber).chan),'_stimChan',num2str(chanList(chan)),'_waveNum',num2str(fig),'_PSTH');
            saveFiguresLab(gcf,figDir,fname);
        end
    end
    
    if(waveformsMakeSubplots && saveFigures && ~plotAllOnOneFigure)
        fname = strcat(figPrefix,'nn',num2str(neuronNumber),'_chan',num2str(cds.units(neuronNumber).chan),'_stimChan',num2str(chanList(chan)),'_waveNum',num2str(fig),'_PSTH');
        saveFiguresLab(gcf,figDir,fname);
    end
end
if(saveFigures && plotAllOnOneFigure)
    fname = strcat(figPrefix,'nn',num2str(neuronNumber),'_chan',num2str(cds.units(neuronNumber).chan),'_PSTHall');
    saveFiguresLab(gcf,figDir,fname);
end



end

%% old code below
% % plot a set of zoomed in rasters
% for cond = 1:2 % 1 = no stim, bump, 2 = stim bump
%     subplot(2,1,cond)
%     hold on
% 
%     if(cond==1)
%         trialTable = cds.trials(ctrHoldBumpNoStimMask,:);
%         yLblStr = 'Bump, No Stim Trials';
%     else
%         trialTable = cds.trials(ctrHoldBumpStimMask,:);
%         yLblStr = 'Bump, Stim Trials';
%         xLblStr = 'Time after bump start (ms)';
%     end
% 
%     trialOffset = 0;
%     bumpDirections = unique(trialTable(:,:).bumpDir);
%     for bd = 1:numel(bumpDirections)
%         bumpDir = bumpDirections(bd);
%         trialTableBumpDirMask = trialTable(:,:).bumpDir == bumpDir;
%         trialTableBumpDir = trialTable(trialTableBumpDirMask,:);
%         for trial = 1:size(trialTableBumpDir,1) % for each trial
%             % plot trial related data
%             if(cond == 1)
%                 trialStart = trialTableBumpDir(trial,:).bumpTime;
%             elseif(cond == 2) % stim data
%                 % find set of 10 times in cds.stimOn
%                 stimStart = cds.stimOn(1:stimsPerBump:end);
%                 stimOffset = stimStart - trialTableBumpDir(trial,:).startTime;
%                 stimOffset(stimOffset < 0) = 10000;
%                 [~,stimStartIdx] = min(stimOffset);
%                 stimOnIdx = (stimStartIdx-1)*stimsPerBump + 1;
%                 stimOnTimes = cds.stimOn(stimOnIdx:stimOnIdx+stimsPerBump-1) - cds.stimOn(stimOnIdx);
%                 trialStart = cds.stimOn(stimOnIdx);
%                 stimOffTimes = stimOnTimes + 1/1000; % blank for 1ms
%                 for st = 1:numel(stimOnTimes)
%                     dtStim = stimOffTimes(st) - stimOnTimes(st);
%                     hold on
%                     pRect = [stimOnTimes(st)*1000, trial+trialOffset-0.5, dtStim, 1]; % [x y dx dy]
%                     r=rectangle('Position',pRect,'FaceColor',colorRect,'EdgeColor',colorRect); 
%                 end
%             end
% 
%             % get spike data
%             spikeMask = cds.units(neuronNumber).spikes.ts  > trialStart - preTime & ...
%                     cds.units(neuronNumber).spikes.ts < trialStart + postTime;
%             spikesPlot = cds.units(neuronNumber).spikes.ts(spikeMask) - trialStart;
%             
%             if(cond == 2 && plotSpikesNearArtifacts == 1) % save spike data for plotting
%                 spikesArtifact = [spikesArtifact;spikesPlot+trialStart];
%             end
%             % plot spike data
%             if(verticalLine)
%                 plot([spikesPlot*1000,spikesPlot*1000]', ...
%                     [(trial+trialOffset)*ones(length(spikesPlot),1)-rpLineLength/2,...
%                     (trial+trialOffset)*ones(length(spikesPlot),1)+rpLineLength/2]','k','linewidth',1.5)      
%             else
%                 plot(spikesPlot*1000,(trial+trialOffset)*ones(length(spikesPlot),1),'k.','markerSize',markerSize)
%             end            
%         end
%         % plot horizontal line to denote bump directions
%         trialOffset = trialOffset + size(trialTableBumpDir,1) + 1;
%         if(bd ~= numel(bumpDirections))
%             plot([-preTime*1500,postTime*1500],[trialOffset, trialOffset],'b','linewidth',1.5)
%             trialOffset = trialOffset + 1;
%         end
%     end
%     
%     ylim([-3,trialOffset]+rpLineLength/2)
%     xlim([-preTime*1000,postTime*1000])
%     ylabel(yLblStr)
%     if(cond == 2)
%         xlabel(xLblStr);
%     end
%     formatForLee(gcf);
%     ax=gca;
%     ax.YTick = [0;trialOffset];
%     ax.YMinorTick = 'off';
%     ax.YTickLabel = {num2str(0),num2str(size(trialTable,1))};
% end