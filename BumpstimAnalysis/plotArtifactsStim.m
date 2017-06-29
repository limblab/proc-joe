function [  ] = plotArtifactsStim( cds, neuronNumber, figNum, varargin )
makeFigure = 1;
plotTitle = 1;
titleToPlot = '';
maxArtifactsPerPlot = 5;
rowSubplot = 4; % must be an even integer
colSubplot = 5;
waveformsSentExist = any(isfield(cds,'waveforms'));
for i = 1:2:size(varargin,2)
    switch varargin{i}
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
    end
end

chan = cds.units(neuronNumber).chan;
spikes = cds.units(neuronNumber).spikes.ts;
numArtifacts = size(cds.artifactData.artifact,1);
xData = ((0:1:size(cds.artifactData.artifact,3)-1)-30)/30000*1000;

numArtifactsPerCond = maxArtifactsPerPlot*rowSubplot/2*colSubplot;
if(makeFigure)
    figure();   
end

for artCond = 1:2
    artifactsPlot = [];
    if(artCond == 1) % plot sample of artifacts with spike afterwards
        for art = 1:numArtifacts
            artifactsMask = spikes >= cds.artifactData.t(art) & spikes <= cds.artifactData.t(art) + (size(cds.artifactData.artifact,3)-30)/30000;
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

    % grab random sample from artifactsPlot
    if(numel(artifactsPlot) > numArtifactsPerCond)
        artifactsPlot = datasample(artifactsPlot,numArtifactsPerCond,'Replace',false);
    end

    artCount = 1;
    for sub = 1:rowSubplot/2*colSubplot
        subplot(rowSubplot,colSubplot, sub+(rowSubplot/2*colSubplot)*(artCond-1))

        for p = 1:maxArtifactsPerPlot
            if(artCount <= numel(artifactsPlot))
                plot(squeeze(squeeze(cds.artifactData.artifact(artifactsPlot(artCount),chan,:))))
            end
            hold on
            artCount = artCount+1;
        end
        ylim([-500 500])
        formatForLee(gcf)
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

