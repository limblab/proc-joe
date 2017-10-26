function [  ] = plotArtifactsStim( cds, neuronNumber, chanNum, figNum, varargin )
makeFigure = 1;
plotTitle = 1;
titleToPlot = '';
maxArtifactsPerPlot = 20;
rowSubplot = 4; % must be an even integer
colSubplot = 5;
waveformsSentExist = any(isfield(cds,'waveforms'));
plotFiltered = 0;
[bFilter,aFilter] = butter(6,[500]/(30000/2),'high');
timeAfterStim = 5/1000;
plotArtifactsSeparated = 1;
randomSample = 1;
stimElectrode = -1;
templateSubtract = 0;

% save stuff
saveFigures = 0;
figDir = '';
figPrefix = '';

rangeProvided = 0;
plotXRange = [];
for i = 1:2:size(varargin,2)
    switch varargin{i}
        case 'plotFiltered'
            plotFiltered = varargin{i+1};
        case 'makeFigure'
            makeFigure = varargin{i+1};
        case 'plotTitle'
            plotTitle = varargin{i+1};
        case 'title'
            titleToPlot = varargin{i+1};
        case 'maxArtifactsPerPlot'
            maxArtifactsPerPlot = varargin{i+1};
        case 'rowSubplot'
            rowSubplotTemp = varargin{i+1};
            if(mod(rowSubplotTemp,2)==0)
                rowSubplot = rowSubplotTemp;
            else
                rowSubplot = max(2,rowSubplotTemp-1);
            end
        case 'colSubplot'
            colSubplot = varargin{i+1};
        case 'timeAfterStim'
            timeAfterStim = varargin{i+1};
        case 'plotArtifactsSeparated'
            plotArtifactsSeparated = varargin{i+1};
        case 'randomSample'
            randomSample = varargin{i+1};
        case 'stimElectrode'
            stimElectrode = varargin{i+1};
        case 'saveFigures'
            saveFigures = varargin{i+1};
        case 'figDir'
            figDir = varargin{i+1};
        case 'figPrefix'
            figPrefix = varargin{i+1};
        case 'plotXRange'
            rangeProvided = 1;
            plotXRange = varargin{i+1};
            if(plotXRange(2) < 21)
                plotXRange(2) = 21;
            end
        case 'templateSubtract'
            templateSubtract = varargin{i+1};
    end
end

if(saveFigures && strcmp(figDir,''))
    saveFigures = 0;
end

if(waveformsSentExist && any(isfield(cds.waveforms,'chanSent')))
    chanList = unique(cds.waveforms.chanSent);
    numChans = numel(chanList);
else
    chanList = stimElectrode;
    numChans = 1;
end

neuronChan = cds.units(neuronNumber).chan;
spikes = cds.units(neuronNumber).spikes.ts;
numArtifacts = size(cds.artifactData.artifact,1);
xData = ((0:1:size(cds.artifactData.artifact,3)-1)-30)/30000*1000;

if(plotArtifactsSeparated)
    numArtifactsPerCond = maxArtifactsPerPlot*rowSubplot/2*colSubplot;
    maxArtCond = 2;
else
    numArtifactsPerCond = maxArtifactsPerPlot*rowSubplot*colSubplot;
    maxArtCond = 1;
end

if(makeFigure)
    figure();   
end


% [coeff,score] = pca(squeeze(cds.artifactData.artifact(:,chan,:)));
% coeff(:,1:3) = 0;
% artifact=coeff*score';
for artCond = 1:maxArtCond
    artifactsPlot = [];
    
    if(plotArtifactsSeparated)
        if(artCond == 1) % plot sample of artifacts with spike afterwards
            for art = 1:numArtifacts
                artifactsMask = spikes >= cds.artifactData.t(art) & spikes <= cds.artifactData.t(art) + timeAfterStim;
                if(sum(artifactsMask)>0 && ((waveformsSentExist && cds.waveforms.waveSent(art) == figNum) || ~waveformsSentExist) && (~any(isfield(cds.waveforms,'chanSent')) || cds.waveforms.chanSent(art)==chanList(chanNum)))
                    artifactsPlot(end+1,1) = art;
                end
            end
        else % plot sample of artifacts without spike afterwards 
            for art = 1:numArtifacts
                artifactsMask = spikes >= cds.artifactData.t(art) & spikes <= cds.artifactData.t(art) + timeAfterStim;

                if(sum(artifactsMask)==0 && cds.artifactData.t(art) ~= 0 && ((waveformsSentExist && cds.waveforms.waveSent(art) == figNum) || ~waveformsSentExist) && (~any(isfield(cds.waveforms,'chanSent')) || cds.waveforms.chanSent(art)==chanList(chanNum)))
                    artifactsPlot(end+1,1) = art;
                end
            end
        end
    else
        for art = 1:numArtifacts
            if(((waveformsSentExist && cds.waveforms.waveSent(art) == figNum) || ~waveformsSentExist) && (~any(isfield(cds.waveforms,'chanSent')) || cds.waveforms.chanSent(art)==chanList(chanNum)))
                artifactsPlot(end+1,1) = art;
            end
        end
    end

    % grab random sample from artifactsPlot
    if(numel(artifactsPlot) > numArtifactsPerCond)
        if(randomSample)
           artifactsPlot = datasample(artifactsPlot,numArtifactsPerCond,'Replace',false);
        else
           artifactsPlot = artifactsPlot(end-numArtifactsPerCond:end,1);
        end
    end

    artCount = 1;
    if(plotArtifactsSeparated)
        for sub = 1:rowSubplot/2*colSubplot
            subplot(rowSubplot,colSubplot, sub+(rowSubplot/2*colSubplot)*(artCond-1))

            for p = 1:maxArtifactsPerPlot
                if(artCount <= numel(artifactsPlot))
                    plotStimDataArtifactStim(cds,artifactsPlot,artCount,neuronChan,rangeProvided,plotXRange,plotFiltered,...
                        bFilter,aFilter,templateSubtract,waveformsSentExist,figNum,chanList,chanNum)

                end
                hold on
                artCount = artCount+1;
            end
%             ylim([-400 400])
%             xlabel('Time after stimulation onset (ms)')
%             ylabel('Voltage (\muV)')
            formatForLee(gcf)
        end
    else
        for sub = 1:rowSubplot*colSubplot
            subplot(rowSubplot,colSubplot, sub)
            for p = 1:maxArtifactsPerPlot
                if(artCount <= numel(artifactsPlot))
                    plotStimDataArtifactStim(cds,artifactsPlot,artCount,neuronChan,rangeProvided,plotXRange,plotFiltered,...
                        bFilter,aFilter,templateSubtract,waveformsSentExist,figNum,chanList,chanNum)
                end
                hold on
                artCount = artCount+1;
%                 ylim([-400,400])
            end
%             plot((stimDataPlot(1:end-200,:) - mean(stimDataPlot(1:end-200,:),2)))

%             ylim([-400 400])
%             xlabel('Time after stimulation onset (ms)')
%             ylabel('Voltage (\muV)')
            formatForLee(gcf)
        end
    end

end

if(plotTitle)
    if(strcmp(titleToPlot,'') == 0)
        suptitle(titleToPlot);
    elseif(numChans > 1 && waveformsSentExist)
        suptitle(strcat('Stim Chan: ',num2str(chanList(chanNum)),' Wave: ',num2str(figNum)));
    elseif(numChans == 1 && waveformsSentExist)
        suptitle(strcat('Stim Chan: ',num2str(stimElectrode),' Wave: ',num2str(figNum)));
    else

    end
end
   
if(saveFigures)
    if(plotFiltered)
        fname = strcat(figPrefix,'nn',num2str(neuronNumber),'_chan',num2str(cds.units(neuronNumber).chan),'_stimChan',num2str(chanList(chanNum)),'_waveNum',num2str(figNum),'_artifactDataFiltered');
    else
        fname = strcat(figPrefix,'nn',num2str(neuronNumber),'_chan',num2str(cds.units(neuronNumber).chan),'_stimChan',num2str(chanList(chanNum)),'_waveNum',num2str(figNum),'_artifactDataRaw');
    end
    saveFiguresLIB(gcf,figDir,fname);
end


end

function [] = plotStimDataArtifactStim(cds,artifactsPlot,artCount,neuronChan,rangeProvided,plotXRange,plotFiltered,bFilter,aFilter,templateSubtract,waveformsSentExist,figNum,chanList,chanNum)

if(plotFiltered && ~any(isfield(cds.artifactData,'artifactProcessed')))
    if(rangeProvided)
        stimData = squeeze(cds.artifactData.artifact(artifactsPlot(artCount),neuronChan,plotXRange(1):plotXRange(2)));
    else
        stimData = squeeze(cds.artifactData.artifact(artifactsPlot(artCount),neuronChan,:));
    end
    stimData = [stimData;mean(stimData(end-20:end))*ones(200,1)];
    stimData = fliplr(filter(bFilter,aFilter,fliplr(stimData')))';
%     if(templateSubtract)
%         templateIdx = [];
%         for a = 1:size(cds.stimOn,1)
%             if(((waveformsSentExist && cds.waveforms.waveSent(a) == figNum) || ~waveformsSentExist) && (~any(isfield(cds.waveforms,'chanSent')) || cds.waveforms.chanSent(a)==chanList(chanNum)))
%                 templateIdx(end+1,1) = a;
%             end
%         end
%         templateData = squeeze(cds.artifactData.artifact(templateIdx,neuronChan,plotXRange(1):plotXRange(2)));
%         templateData = [templateData,mean(templateData(:,end-20:end),2).*ones(size(templateData,1),200)];
%         templateData = fliplr(filter(bFilter,aFilter,fliplr(templateData)')');
%         template = mean(templateData);
%         stimData(1:numel(template),1) = stimData(1:numel(template),1) - template';
%     end
    plot((0:1:(numel(stimData)-201))/30,stimData(1:end-200)/0.254)
%     [coeff,score] = pca(squeeze(cds.artifactData.artifact(artifactsPlot(artCount),:,:))');
%     coeff(:,1:4) = 0;
%     artifact=coeff*score';
%     stimData = artifact(chan,:)';
% 
%     stimData = artifact(:,artifactsPlot(artCount));
% 
%     stimData = [stimData;zeros(200,1)];
%     stimDataPlot = fliplr(filter(bFilter,aFilter,fliplr(stimData')))';
%     plot((0:1:(numel(stimDataPlot)-201))/30,stimDataPlot(1:end-200))
elseif(plotFiltered && any(isfield(cds.artifactData,'artifactProcessed')))
    if(rangeProvided)
        plot((0:1:(numel(squeeze(cds.artifactData.artifactProcessed(artifactsPlot(artCount),neuronChan,plotXRange(1):plotXRange(2))))-1))/30,...
            squeeze(cds.artifactData.artifactProcessed(artifactsPlot(artCount),neuronChan,plotXRange(1):plotXRange(2)))/0.254)
    else
        plot((0:1:(numel(squeeze(cds.artifactData.artifactProcessed(artifactsPlot(artCount),neuronChan,:)))-1))/30,...
            squeeze(cds.artifactData.artifactProcessed(artifactsPlot(artCount),neuronChan,:))/0.254)
    end
else    
    if(rangeProvided)
        plot((0:1:(numel(squeeze(cds.artifactData.artifact(artifactsPlot(artCount),neuronChan,plotXRange(1):plotXRange(2))))-1))/30,...
            squeeze(cds.artifactData.artifact(artifactsPlot(artCount),neuronChan,plotXRange(1):plotXRange(2)))/0.254)
    else
        plot((0:1:(numel(squeeze(cds.artifactData.artifact(artifactsPlot(artCount),neuronChan,:)))-1))/30,...
            squeeze(cds.artifactData.artifact(artifactsPlot(artCount),neuronChan,:))/0.254)
    end
end



end




