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

% save stuff
saveFigures = 0;
figDir = '';
figPrefix = '';

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
numArtifacts = numel(cds.stimOn);
xData = ((0:1:size(cds.artifactData.artifact,3)-1)-30)/30000*1000;

if(plotArtifactsSeparated)
    numArtifactsPerCond = maxArtifactsPerPlot*rowSubplot/2*colSubplot;
    maxArtCond = 2;
else
    numArtifactsPerCond = maxArtifactsPerPlot*rowSubplot*colSubplot;
    maxArtCond = 1;
end

if(makeFigure && ~plotArtifactsSeparated)
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
                artifactsMask = spikes >= cds.artifactData.t(art) & spikes <= cds.artifactData.t(art) + (size(cds.artifactData.artifact,3)-30)/30000;

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
           artifactsPlot = artifactsPlot(1:numArtifactsPerCond,1);
        end
    end

    artCount = 1;
    if(plotArtifactsSeparated)
        for sub = 1:rowSubplot/2*colSubplot
            subplot(rowSubplot,colSubplot, sub+(rowSubplot/2*colSubplot)*(artCond-1))

            for p = 1:maxArtifactsPerPlot
                if(artCount <= numel(artifactsPlot))
                    if(plotFiltered)
                        stimData = squeeze(cds.artifactData.artifact(artifactsPlot(artCount),neuronChan,:));
                        stimData = [stimData;zeros(200,1)];
                        stimDataPlot = fliplr(filter(bFilter,aFilter,fliplr(stimData')))';
                        plot((0:1:(numel(stimDataPlot)-201))/30,stimDataPlot(1:end-200)/0.254)
%                         [coeff,score] = pca(squeeze(cds.artifactData.artifact(artifactsPlot(artCount),:,:))');
%                         coeff(:,1:4) = 0;
%                         artifact=coeff*score';
%                         stimData = artifact(chan,:)';
                                             
%                         stimData = artifact(:,artifactsPlot(artCount));
%                         
%                         stimData = [stimData;zeros(200,1)];
%                         stimDataPlot = fliplr(filter(bFilter,aFilter,fliplr(stimData')))';
%                         plot((0:1:(numel(stimDataPlot)-201))/30,stimDataPlot(1:end-200))   
                    else
                        plot((0:1:(numel(squeeze(cds.artifactData.artifact(artifactsPlot(artCount),neuronChan,:)))-1))/30,...
                            squeeze(cds.artifactData.artifact(artifactsPlot(artCount),neuronChan,:))/0.254)
                    end
                end
                hold on
                artCount = artCount+1;
            end
            ylim([-1200 1200])
            xlabel('Time after stimulation onset (ms)')
            ylabel('Voltage (\muV)')
            formatForLee(gcf)
        end
    else
        for sub = 1:rowSubplot*colSubplot
%             subplot(rowSubplot,colSubplot, sub)
            for p = 1:maxArtifactsPerPlot
                if(artCount <= numel(artifactsPlot))
                    if(plotFiltered)
                        stimData = squeeze(cds.artifactData.artifact(artifactsPlot(artCount),neuronChan,:));
                        stimData = [stimData;zeros(200,1)];
                        stimDataPlot = fliplr(filter(bFilter,aFilter,fliplr(stimData')))';
                        plot((0:1:(numel(stimDataPlot)-201))/30,stimDataPlot(1:end-200)/0.254)
%                         [coeff,score] = pca(squeeze(cds.artifactData.artifact(artifactsPlot(artCount),:,:))');
%                         coeff(:,1:4) = 0;
%                         artifact=coeff*score';
%                         stimData = artifact(chan,:)';
% %                         
%                         stimData = artifact(:,artifactsPlot(artCount));
%                         
%                         stimData = [stimData;zeros(200,1)];
%                         stimDataPlot = fliplr(filter(bFilter,aFilter,fliplr(stimData')))';
%                         plot((0:1:(numel(stimDataPlot)-201))/30,stimDataPlot(1:end-200))
                    else
                        plot((0:1:(numel(squeeze(cds.artifactData.artifact(artifactsPlot(artCount),neuronChan,:)))-1))/30,...
                            squeeze(cds.artifactData.artifact(artifactsPlot(artCount),neuronChan,:))/0.254)
                    end
                end
                hold on
                artCount = artCount+1;
            end
            ylim([-1200 1200])
            xlabel('Time after stimulation onset (ms)')
            ylabel('Voltage (\muV)')
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
    saveFiguresLab(gcf,figDir,fname);
end


end

