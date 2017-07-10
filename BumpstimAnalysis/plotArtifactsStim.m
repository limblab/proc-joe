function [  ] = plotArtifactsStim( cds, neuronNumber, figNum, varargin )
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
    end
end

chan = cds.units(neuronNumber).chan;
spikes = cds.units(neuronNumber).spikes.ts;
numArtifacts = numel(cds.stimOn);
xData = ((0:1:size(cds.artifactData.artifact,3)-1)-30)/30000*1000;

if(plotArtifactsSeparated)
    numArtifactsPerCond = maxArtifactsPerPlot*rowSubplot/2*colSubplot;
else
    numArtifactsPerCond = maxArtifactsPerPlot*rowSubplot*colSubplot;
end

if(makeFigure)
    figure();   
end


if(plotArtifactsSeparated)
    maxArtCond = 2;
else
    maxArtCond = 1;
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
                if(sum(artifactsMask)>0 && ((waveformsSentExist && cds.waveforms.waveSent(art) == figNum) || ~waveformsSentExist))
                    artifactsPlot(end+1,1) = art;
                end
            end
        else % plot sample of artifacts without spike afterwards 
            for art = 1:numArtifacts
                artifactsMask = spikes >= cds.artifactData.t(art) & spikes <= cds.artifactData.t(art) + (size(cds.artifactData.artifact,3)-30)/30000;

                if(sum(artifactsMask)==0 && cds.artifactData.t(art) ~= 0)
                    artifactsPlot(end+1,1) = art;
                end
            end
        end
    else
        for art = 1:numArtifacts
            if(((waveformsSentExist && cds.waveforms.waveSent(art) == figNum) || ~waveformsSentExist))
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
                        stimData = squeeze(cds.artifactData.artifact(artifactsPlot(artCount),chan,:));
                        stimData = [stimData;zeros(200,1)];
                        stimDataPlot = fliplr(filter(bFilter,aFilter,fliplr(stimData')))';
                        plot((0:1:(numel(stimDataPlot)-201))/30,stimDataPlot(1:end-200))
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
                        plot((0:1:(numel(squeeze(cds.artifactData.artifact(artifactsPlot(artCount),chan,:)))-1))/30,...
                            squeeze(cds.artifactData.artifact(artifactsPlot(artCount),chan,:)))
                    end
                end
                hold on
                artCount = artCount+1;
            end
            ylim([-500 500])
            formatForLee(gcf)
        end
    else
        for sub = 1:rowSubplot*colSubplot
%             subplot(rowSubplot,colSubplot, sub)
            figure
            for p = 1:maxArtifactsPerPlot
                if(artCount <= numel(artifactsPlot))
                    if(plotFiltered)
                        stimData = squeeze(cds.artifactData.artifact(artifactsPlot(artCount),chan,:));
                        stimData = [stimData;zeros(200,1)];
                        stimDataPlot = fliplr(filter(bFilter,aFilter,fliplr(stimData')))';
                        plot((0:1:(numel(stimDataPlot)-201))/30,stimDataPlot(1:end-200))
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
                        plot((0:1:(numel(squeeze(cds.artifactData.artifact(artifactsPlot(artCount),chan,:)))-1))/30,...
                            squeeze(cds.artifactData.artifact(artifactsPlot(artCount),chan,:)))
                    end
                end
                hold on
                artCount = artCount+1;
            end
            ylim([-500 500])
            formatForLee(gcf)
        end
    end

end

if(plotTitle)
    if(strcmp(titleToPlot,'') == 0)
        suptitle(titleToPlot);
    elseif(waveformsSentExist)
        suptitle(num2str(figNum));
    else
        suptitle(num2str(neuronNumber));
    end
end
    

end

