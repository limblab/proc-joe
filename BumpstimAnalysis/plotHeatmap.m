function [ output_args ] = plotHeatmap( cds,mapFilename, varargin )
% plot PSTH for the whole array

plotLine = 1;
binSize = 0.0002;
waveformTypes = 1:numel(cds.waveforms.parameters);
chans = 1:numel(unique(cds.waveforms.chanSent));
lineColor = 'k';
lineWidth = 1;
plotAllOnOneFigure = 0;
plotStimOn = 1;
plotStimChan = 1;
plotProbabilityText = 0;
stimPreTime = 0;
stimPostTime = 5/1000;
baselinePreTime = -10/1000;
baselinePostTime = -2/1000;
saveFigures = 0;
for i = 1:2:numel(varargin)
    switch varargin{i}
        case 'plotLine'
            plotLine = varargin{i+1};
        case 'lineColor'
            lineColor = varargin{i+1};
        case 'plotAllOnOneFigure'
            plotAllOnOneFigure = varargin{i+1};
        case 'stimPreTime'
            stimPreTime = varargin{i+1};
        case 'stimPostTime'
            stimPostTime = varargin{i+1};
        case 'baselinePreTime'
            baselinePreTime = varargin{i+1};
        case 'baselinePostTime'
            baselinePostTime = varargin{i+1};
        case 'binSize'
            binSize = varargin{i+1};
        case 'waveformTypes'
            waveformTypes = varargin{i+1};
        case 'chans'
            chans = varargin{i+1};
        case 'plotStimOn'
            plotStimOn = varargin{i+1};
        case 'plotStimChan'
            plotStimChan = varargin{i+1};
        case 'plotProbabilityText'
            plotProbabilityText = varargin{i+1};
        case 'saveFigures'
            saveFigures = varargin{i+1};
        case 'figDir'
            figDir = varargin{i+1};
        case 'figPrefix'
            figPrefix = varargin{i+1};
    end
end

mapData=loadMapFile(mapFilename);

%establish trackers for information that is common to all files:
posList=[];
posList = [mapData.row,mapData.col];
plottedHere = zeros(10,10);

% extract data and put in matrix 
dataPost = nan(numel(chans),numel(waveformTypes),10,10);
dataPre = nan(numel(chans),numel(waveformTypes),10,10);
        
for nn = 1:size(cds.units,2)
    if(cds.units(nn).ID~=0 && cds.units(nn).ID~=255)

        posIdx=find(mapData.chan==cds.units(nn).chan);
        eRow=posList(posIdx,1);
        eRow = 11 - eRow;
        eCol=posList(posIdx,2);
%                 h=subplot(10,10,10*(eRow-1)+eCol);
        plottedHere(eRow,eCol) = 1;

        [bE,bC,bCVar] = plotPSTHStim(cds,nn,'binSize',binSize,'makeFigure',0,'makeSubplots',0,'plotTitle',0,'waveformTypes',waveformTypes,...
                'chans',chans,'preTime',-1*(baselinePreTime-5/1000),'postTime',stimPostTime+5/1000,'lineColor',lineColor,'plotLine',plotLine,'plotStimOn',plotStimOn,...
                'lineWidth',lineWidth,'noPlot',1);
        
        for c = 1:numel(chans)
            for wave = 1:numel(waveformTypes)
                binIdx = [max(find(bE{chans(c),waveformTypes(wave)}(:) <= stimPreTime*1000)), min(find(bE{chans(c),waveformTypes(wave)}(:) >= stimPostTime*1000));...
                    max(find(bE{chans(c),waveformTypes(wave)}(:) <= baselinePreTime*1000)), min(find(bE{chans(c),waveformTypes(wave)}(:) >= baselinePostTime*1000))];

                posIdx=find(mapData.chan==cds.units(nn).chan);
                eRow=posList(posIdx,1);
%                 eRow = 11 - eRow;
                eCol=posList(posIdx,2);
                dataPost(c,wave,eRow,eCol) = mean(bC{chans(c),waveformTypes(wave)}(binIdx(1,1):binIdx(1,2)));
                dataPre(c,wave,eRow,eCol) = mean(bC{chans(c),waveformTypes(wave)}(binIdx(2,1):binIdx(2,2)));
                dataPostVar(c,wave,eRow,eCol) = bCVar{chans(c),waveformTypes(wave)}(2,2);
                dataPreVar(c,wave,eRow,eCol) = bCVar{chans(c),waveformTypes(wave)}(2,1);
            end
        end
    end
end

% get colors
colormap jet;
colors = colormap;

%% plot data for each channel
% ratios of dataPost and dataPre
dataAll = dataPost./dataPre;
dataAll(dataAll > 8) = 8;
% t-test for significance, plot significance (remove non-significant data)
% pVal = 0.05/96/numel(chans)/numel(waveformTypes);
% 
% numStims = zeros(numel(chans),numel(waveformTypes));
% chanList = unique(cds.waveforms.chanSent);
% for c  = 1:numel(chans)
%     for w = 1:numel(waveformTypes)
%         numStims(c,w) = sum(cds.waveforms.waveSent == w & cds.waveforms.chanSent == chanList(c));
%     end
% end
% sAdjusted = sqrt(dataPostVar./numStims + dataPreVar./numStims);
% tAll = (dataPost - dataPre)./sAdjusted;
% dfAll = (sAdjusted.^2)./(((dataPostVar./numStims).^2)./(numStims-1) + ((dataPreVar./numStims).^2)./(numStims-1));
% tStatAll = tcdf(tAll,dfAll,'upper');
% dataAll = NaN(size(tStatAll));
% dataAll(tStatAll(:,:,:,:) < pVal) = tStatAll(tStatAll(:,:,:,:) < pVal);
% dataAll(dataAll(:,:,:,:) < 1) = nan;
% dataAll(4,1,3,5) = 8;
maxData = max(max(max(max(dataAll))));
minData = min(min(min(min(dataAll))));
minData = 0;
for c = 1:numel(chans)
    for wave = 1:numel(waveformTypes)
        data = squeeze(dataAll(c,wave,:,:));       
        % make heatmap
        figure;        
        % create ROWSxCOLS grid
        surface(zeros(11,11))
        ax = gca;
        ax.Children(1).LineWidth = 1.5;
        colormap([1,1,1]);
        f=gcf;
        ax = gca;
        ax.YTickLabel = {};
        ax.XTickLabel = {};
        hold on
        for rowElec = 1:size(data,1)
            for colElec = 1:size(data,2)
                if(isnan(data(rowElec,colElec)))
%                     plot([colElec,colElec+1],[rowElec,rowElec+1],'k','linewidth',2);
%                     plot([colElec+1,colElec],[rowElec,rowElec+1],'k','linewidth',2);
                else
                    colorIdx = floor(size(colors,1)*((data(rowElec,colElec)-minData)/(maxData-minData)));
                    if(colorIdx < 1)
                        colorIdx = 1;
                    elseif(colorIdx > size(colors,1))
                        colorIdx = size(colors,1);
                    end
                    rectangle('Position',[colElec,rowElec,1,1],'EdgeColor','k',...
                            'FaceColor',colors(colorIdx,:),'linewidth',0.1);
                end
            end
        end
        
        % plot stim chan
        stimChan = unique(cds.waveforms.chanSent);
        stimChan = stimChan(chans(c));
        posIdx=find(mapData.chan==stimChan);
        eRow=posList(posIdx,1);
%         eRow = 11 - eRow;
        eCol=posList(posIdx,2);
        
        % plot box
        plot([eCol,eCol+1,eCol+1,eCol,eCol],[eRow,eRow,eRow+1,eRow+1,eRow],'m','linewidth',3)
        % plot magenta x
        if(isnan(data(eRow,eCol)))
            plot([eCol,eCol+1],[eRow,eRow+1],'m','linewidth',3);
            plot([eCol+1,eCol],[eRow,eRow+1],'m','linewidth',3);
        end
        set(gca,'Visible','off')
        
        if(saveFigures)
            fname = strcat(figPrefix,'_stimChan',num2str(stimChan),'_waveNum',num2str(wave),'_heatmap');
            saveFiguresLIB(gcf,figDir,fname);
        end
    end
end


%% make colorbar
figure
b=colorbar;
colormap jet;
set(gca,'Visible','off');
b.FontSize = 14;
b.Ticks = [0,0.25,0.5,0.75,1.0];
maxDataRound = round(maxData,1);
minDataRound = round(minData,1);
b.TickLabels = {};
for i = b.Ticks
    
    if(i==b.Ticks(end))
        b.TickLabels{end+1,1} = strcat('>',num2str(i*(maxDataRound-minDataRound) + minDataRound));
    else
        b.TickLabels{end+1,1} = num2str(i*(maxDataRound-minDataRound) + minDataRound);
    end
end
if(saveFigures)
    fname = strcat(figPrefix,'_heatmapScale');
    saveFiguresLIB(gcf,figDir,fname);
end


end

