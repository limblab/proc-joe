function [  ] = plotArrayMap( cds,neuronNumber,mapFileName, varargin )
% this function plots things on a 10x10 grid of the array. Namely the
% stimulated and recording channel, though other functionality may be added
% in the future

% electrodea are passed as a channel -- not 'elec##'

numRows = 10+1;
numCols = 10+1;
makeFigure = 1;
userColorMap = [1,1,1]; % white is default
stimElectrode = [];
stimElectrodeColor = 'k';
stimElectrodeLabel = 'string';
recordingElectrode = [];
recordingElectrodeColor = 'k';
recordingElectrodeLabel = 'string';

% save stuff
saveFigures = 0;
figDir = '';
figPrefix = '';

%% deal with varagin
for i = 1:2:size(varargin,2)
    switch varargin{i}
        case 'numRows'
            numRows = varargin{i+1}+1;
        case 'numCols'
            numCols = varargin{i+1}+1;
        case 'makeFigure'
            makeFigure = varargin{i+1};
        case 'colormap'
            userColorMap = varargin{i+1};
        case 'stimElectrode'
            stimElectrode = varargin{i+1};
        case 'stimElectrodeColor'
            stimElectrodeColor = varargin{i+1};
        case 'stimElectrodeLabel'
            stimElectrodeLabel = varargin{i+1};
        case 'recordingElectrode'
            recordingElectrode = varargin{i+1};
        case 'recordingElectrodeColor'
            recordingElectrodeColor = varargin{i+1};
        case 'recordingElectrodeLabel'
            recordingElectrodeLabel = varargin{i+1};
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

%% load map file
arrayMap=loadMapFile(mapFileName);

%% make figure if necessary
if(makeFigure)
    figure();
end

%% remove axis labels
f=gcf;
ax = gca;
ax.YTickLabel = {};
ax.XTickLabel = {};

%% create ROWSxCOLS grid
surface(zeros(numRows,numCols))
colormap(userColorMap);

%% label stim electrode -- currently using an 'S' though this can be changed maybe later
if(~isempty(stimElectrode))
    for stElec = 1:numel(stimElectrode)
        arrayMapIdx = find(arrayMap.chan == stimElectrode(stElec));
        if(~isempty(arrayMapIdx))
            rowElec = arrayMap.row(arrayMapIdx);
            colElec = arrayMap.col(arrayMapIdx);
            if(strcmpi(stimElectrodeLabel,'string'))
                if(numel(stimElectrodeColor) > 1)
                    text(colElec+0.25,rowElec+0.5,'S','FontWeight','bold','fontsize',20,'color',stimElectrodeColor{stElec});
                else
                    text(colElec+0.25,rowElec+0.5,'S','FontWeight','bold','fontsize',20,'color',stimElectrodeColor);
                end
            elseif(strcmpi(stimElectrodeLabel,'box'))
                if(numel(stimElectrodeColor) > 1)
                    rectangle('Position',[colElec,rowElec,1,1],'EdgeColor',stimElectrodeColor{stElec},...
                        'FaceColor',stimElectrodeColor,'linewidth',0.1);
                else
                    rectangle('Position',[colElec,rowElec,1,1],'EdgeColor',stimElectrodeColor,...
                        'FaceColor',stimElectrodeColor,'linewidth',0.1);
                end
            end
        end
    end
end
%% label recording electrode
if(~isempty(recordingElectrode))
    for recElec = 1:numel(recordingElectrode)
        arrayMapIdx = find(arrayMap.chan == recordingElectrode(recElec));
        if(~isempty(arrayMapIdx))
            rowElec = arrayMap.row(arrayMapIdx);
            colElec = arrayMap.col(arrayMapIdx);
            if(strcmpi(recordingElectrodeLabel,'string'))
                if(numel(recordingElectrodeColor)>1)
                    text(colElec+0.25,rowElec+0.5,'R','FontWeight','bold','fontsize',20,'color',recordingElectrodeColor{recElec});
                else
                    text(colElec+0.25,rowElec+0.5,'R','FontWeight','bold','fontsize',20,'color',recordingElectrodeColor);
                end
            elseif(strcmpi(recordingElectrodeLabel,'box'))
                rectangle('Position',[colElec,rowElec,1,1],'EdgeColor',recordingElectrodeColor,...
                    'FaceColor',recordingElectrodeColor,'linewidth',0.1);
                if(numel(recordingElectrodeColor)>1)
                    rectangle('Position',[colElec,rowElec,1,1],'EdgeColor',recordingElectrodeColor,...
                        'FaceColor',recordingElectrodeColor{recElec},'linewidth',0.1);
                else
                    rectangle('Position',[colElec,rowElec,1,1],'EdgeColor',recordingElectrodeColor,...
                        'FaceColor',recordingElectrodeColor,'linewidth',0.1);
                end
            end
        end
    end
end

ax = gca;
set(ax,'Visible','off')
if(saveFigures)
    fname = strcat(figPrefix,'nn',num2str(neuronNumber),'_chan',num2str(cds.units(neuronNumber).chan),'_arrayMap');
    saveFiguresLIB(gcf,figDir,fname);
end

end

