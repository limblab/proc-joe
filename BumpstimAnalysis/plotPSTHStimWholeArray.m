function [ output_args ] = plotPSTHStimWholeArray( cds,mapFilename, varargin )
% plot PSTH for the whole array

plotLine = 1;
preTime = 10/1000;
postTime = 30/1000;
binSize = 0.0002;
waveformTypes = 1:numel(cds.waveforms.parameters);
chans = 1:numel(unique(cds.waveforms.chanSent));
for i = 1:2:numel(varargin)
    switch varargin{i}
        case 'plotLine'
            plotLine = varargin{i+1};
        case 'preTime'
            preTime = varargin{i+1};
        case 'postTime'
            postTime = varargin{i+1};
        case 'binSize'
            binSize = varargin{i+1};
        case 'waveformTypes'
            waveformTypes = varargin{i+1};
        case 'chans'
            chans = varargin{i+1};
    end
end

mapData=loadMapFile(mapFilename);

%establish trackers for information that is common to all files:
posList=[];

figure();

for nn = 1:size(cds.units,2)
    if(cds.units(nn).ID~=0 && cds.units(nn).ID~=255)
        posIdx=find(strcmp(mapData.chan,cds.units(nn).chan));
        eRow=posList(posIdx,1);
        eCol=posList(posIdx,2);
        h=subplot(10,10,10*(eRow-1)+eCol);
        
        plotPSTHStim(cds,nn,'binSize',binSize,'makeFigure',0,'makeSubplots',0,'plotTitle',0,'waveformTypes',waveformTypes,...
                'chans',[1:1:numel(unique(cds.waveforms.chanSent))],'preTime',preTime,'postTime',postTime)
    end
end

end

