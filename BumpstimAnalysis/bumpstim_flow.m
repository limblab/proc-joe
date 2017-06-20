%% load in processed data
pwd = cd;
folderpath= 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\Han_processed\20170614\';
mapFile = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
cd(folderpath)
filelist = dir('*processed*');
if(numel(filelist) > 0)
    load(filelist(1).name);
else
    disp('No processed data file. Process data before trying to analyze');
end
cd(pwd);
if(cds.processed == 0)
    disp('This needs to be processed offline more')
end
unitsRemove = [];
%% make sure to include the stimulation electrode and mapfile data
cds.stimElectrode = 'elec64';
mapData=loadMapFile(mapFile);
cds.posList = [mapData.row(:) mapData.col(:)];
cds.elecList = mapData.label;
cds.chanList = mapData.chan;


%% for every unit -- plot 500 raw waveforms for verification. remove if too not
% neural like
unitsRemove = [];
unitsResort = [];
numWaves = 500;
for nn = 1:numel(cds.units)
    if(cds.units(nn).ID ~= 0 && cds.units(nn).ID ~= 255)
        figure
        suptitle(num2str(nn))
        subplotNum = 1;
        subplot(10,10,subplotNum)
        hold on
        for wave = 1:numWaves
            waveIdx = floor(rand()*size(cds.units(nn).spikes.ts,1)) + 1;
            rawIdx = getRawDataIdx(cds.units(nn).spikes.ts(waveIdx),cds.units(nn).chan,...
                cds.rawData.ts,cds.rawData.elec);
            if(rawIdx ~= -1)
                plot(cds.rawData.waveforms(rawIdx,:))
                ylim([-400,400])
            end
            if(mod(wave,5)==0 && wave ~= 500) % edge case
                subplotNum = subplotNum+1;
                subplot(10,10,subplotNum)
                hold on
            end
        end
        inp = -1;
        while(inp~=0 && inp~=1 && inp~=2)
            inp = input('Enter a number: 0 to discard, 1 to keep, 2 to re-sort: ');
            if(~isnumeric(inp))
                inp = -1;
            end
        end
        
        if(inp==0)
            unitsRemove(end+1,1) = nn;
        elseif(inp==2)
            unitsResort(end+1,1) = nn;
        end
    end
end

%% for every unit -- use inverseWaveform to get neural waveform. save all
% of those. If does not look like neuron, remove from units
cds.predictedNeuralWaves = [];
cds.predictedNeuralWaves.neuronNumber = [];
cds.predictedNeuralWaves.waves = [];
figure();
plotCount = 1;
maxPlots = 25;
for nn = 1:numel(cds.units)
    if(cds.units(nn).ID ~= 0 && cds.units(nn).ID ~= 255)
        subplot(5,5,plotCount)
        wave = mean(cds.units(nn).spikes{:,2:end});
        [neuralWave,numGuesses] = runGradDescentCOBumpstim(wave,cds.filter.b,cds.filter.a);
        cds.predictedNeuralWaves.neuronNumber(end+1,1) = nn;
        cds.predictedNeuralWaves.waves(end+1,:) = neuralWave(numGuesses,:);
    
        plot(neuralWave(numGuesses,:))
        title(num2str(nn));
        plotCount = plotCount + 1;
        if(plotCount > maxPlots)
            figure();
            plotCount = 1;
        end
    end
end

%% remove units
if(cds.processed == 0)
    cds.unitsRaw = cds.units;
    cds.unitsRemove = unitsRemove;
    unitsRemove = unique(unitsRemove);
    for unit = numel(unitsRemove):-1:1
        cds.units(unitsRemove(unit)) = [];
    end
end
cds.processed = 1;

%% get row and col in array, then save
cds.unitPos = [];
for nn = 1:numel(cds.units)
    if(cds.units(nn).ID ~= 0 && cds.units(nn).ID ~= 255)
        posIdx = find(strcmpi(cds.elecList,cds.units(nn).label));
        cds.unitPos(end+1,:) = [nn,cds.posList(posIdx,:)];
    end
end
cds.stimElecPos = [];
posIdx = find(strcmpi(cds.elecList,cds.stimElectrode));
cds.stimElecPos(1,:) = cds.posList(posIdx,:);
save([folderpath filelist(1).name],'cds','-v7.3')


%% for every unit -- plot rasters. 
figDir = folderpath;
filename = strcat(cds.meta.monkey,'_20170614_COBumpstim_Raster_neuronNumber');

for nn = 1:numel(cds.units)
    if(cds.units(nn).ID ~= 0 && cds.units(nn).ID ~= 255)
        generateCOBumpstimRaster(cds,nn,'stimsPerBump',cds.stimsPerBump,'trialLength',[0.1,0.2],'plotTitle',0);
%         saveFigure(gcf,figDir,[filename,num2str(nn)]);
%         close all
    end
end

%% look at single units more closely -- plot data around artifact for confirmation
nn = 120;
if(cds.units(nn).ID ~= 0 && cds.units(nn).ID ~= 255)
    % plot all artifacts to see if "spike" is an artifact
    generateCOBumpstimRaster(cds,nn,'stimsPerBump',cds.stimsPerBump,'trialLength',[0.1,0.2],...
        'plotSpikesNearArtifacts',0,'plotTitle',0,'plotAllArtifacts',1);
end


%% for every unit, do a PSTH 
figDir = folderpath;
filename = strcat(cds.meta.monkey,'_20170614_COBumpstim_PSTH_neuronNumber');
for nn = 1:numel(cds.units)
    if(cds.units(nn).ID ~= 0 && cds.units(nn).ID ~= 255)
        generateCOBumpstimPSTH(cds,nn,'stimsPerBump',cds.stimsPerBump,'trialLength',[0.1,0.2],'plotTitle',0,'gaussianSmooth',0,'makeFigure',1);
        saveFigure(gcf,figDir,[filename,num2str(nn)]);
        close all
    end
end

%% look at single PSTH
nn = 120;
if(cds.units(nn).ID ~= 0 && cds.units(nn).ID ~= 255)
    generateCOBumpstimPSTH(cds,nn,'stimsPerBump',cds.stimsPerBump,'trialLength',[0.1,0.2],'plotTitle',0,'gaussianSmooth',0,'makeFigure',1,...
        'makeLegend',0,'makeAxisLabels',1);
end

%% plot PSTH on a 10x10 array grid
figure();
subplotsUsed = [];
for nn = 1:numel(cds.units)
    if(cds.units(nn).ID ~= 0 && cds.units(nn).ID ~= 255) % plot the PSTH
        posIdx = find(cds.unitPos(:,1)==nn);
        eRow=cds.unitPos(posIdx,2);
        eCol=cds.unitPos(posIdx,3);
        h=subplot(10,10,10*(eRow-1)+eCol);
        subplotsUsed(end+1,1) = 10*(eRow-1)+eCol;
        generateCOBumpstimPSTH(cds,nn,'stimsPerBump',cds.stimsPerBump,'trialLength',[0.1,0.2],...
            'plotTitle',0,'gaussianSmooth',0,'makeFigure',0,'makeLegend',0,'makeAxisLabels',0);  
        axis tight%keeps from padding time, we will set y axis below:
        set(gca,'XTickLabel',[])
        set(gca,'YTickLabel',[])
        set(gca,'YMinorTick','off')
        set(gca,'YTick',[])
        set(gca,'XTick',[])
    end
end

% for the stim channel -- if used, plot purple outline, otherwise purple X
subplotStimChannel = 10*(cds.stimElecPos(1)-1)+cds.stimElecPos(2);
h=subplot(10,10,subplotStimChannel);
hold on
if(numel(find(subplotsUsed == sp))==0) % not used
    badMarkX=[0,1];
    badMarkY=[0,1];
    plot(badMarkX,badMarkY,'mp-','lineWidth',3)
    hold on
    badMarkX=[0,1];
    badMarkY=[1,0];
    plot(badMarkX,badMarkY,'mp-','lineWidth',3)
    subplotsUsed(end+1,1) = subplotStimChannel;
else
    ax=gca;
    stimBoxX=[ax.XLim(1),ax.XLim(2),ax.XLim(2),ax.XLim(1),ax.XLim(1)];
    stimBoxY=[ax.YLim(1),ax.YLim(1),ax.YLim(2),ax.YLim(2),ax.YLim(1)];
    plot(stimBoxX,stimBoxY,'mp-','linewidth',3)
end

axis tight%keeps from padding time, we will set y axis below:
set(gca,'XTickLabel',[])
set(gca,'YTickLabel',[])
set(gca,'YMinorTick','off')
set(gca,'YTick',[])
set(gca,'XTick',[])
set(gca,'color','none')
axis off

% for all channels not used, plot a black X
for sp = 1:100
    if(numel(find(subplotsUsed == sp))==0)
        h=subplot(10,10,sp);
        badMarkX=[0,1];
        badMarkY=[0,1];
        plot(badMarkX,badMarkY,'kp-','lineWidth',3)
        hold on
        badMarkX=[0,1];
        badMarkY=[1,0];
        plot(badMarkX,badMarkY,'kp-','lineWidth',3)
        axis tight%keeps from padding time, we will set y axis below:
        set(gca,'XTickLabel',[])
        set(gca,'YTickLabel',[])
        set(gca,'YMinorTick','off')
        set(gca,'YTick',[])
        set(gca,'XTick',[])
        set(gca,'color','none')
        axis off
    end
end




%%% things to do
% 2. PSTH? -- smoothing almost thing
% 3. Spatial (on the array)
