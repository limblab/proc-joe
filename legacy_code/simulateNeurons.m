%% this function takes a cds with artifactData, adds artificial neurons to
% channels without units (determined by cds or by user settings), filters
% and thresholds those artifact snippets and generates a NEV file. Then,
% the simulated neurons can be sorted, cds generated, and a reliability
% metric can be determined.

load('C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\Chips_20171026\Chips_20171026_CObump_chan_all_processed');
load('C:\Users\Joseph\Desktop\Lab\GIT\proc-joe\ProcessStimArtifactData\neuronData_testing.mat');
opts.MAP_FILE_NAME = 'R:\limblab\lab_folder\Animal-Miscellany\Chips_12H1\map_files\left S1\SN 6251-001455.cmp';

%% settings and run all code to generate a NEV file
opts.CHANS_USE = [];
opts.PACKET_WIDTH = 104;
opts.WAVE_AMPLITUDE = 100;
opts.MAX_CHANS_USE = 10;
opts.THRESHOLD_MULT = -4;
opts.IDX_ADD_NEURONS = [25:60];
opts.NEURONS_PER_IDX = 300;
[opts.B, opts.A] =butter(2,250/(30000/2),'high');

% get a random set of channels to use

if(isempty(opts.CHANS_USE)) % find channels based on unit data
    opts.CHANS_USE = 1:96;
    chansWithUnits = [];
    for i = 1:size(cds.units,2)
        if(cds.units(i).ID == 1)
            chansWithUnits(end+1) = cds.units(i).chan;
        end
    end
    chansWithUnits = unique(chansWithUnits);
    opts.CHANS_USE = setdiff(opts.CHANS_USE,chansWithUnits);
    opts.CHANS_USE = sort(datasample(opts.CHANS_USE,opts.MAX_CHANS_USE,'replace',false));
end

% add artificial neurons to artifact data on specified channels
artifact = cds.artifactData.artifact(cds.waveforms.waveSent((0:size(cds.artifactData.artifact,1)-1)*5 + 1)==5,opts.CHANS_USE,:);
artifactData_neurons = addArtificialNeurons(artifact, neuronData, opts);

% get threshold data
% neuralData = thresholdData(artifactData_neurons.artifact,opts);

% % write nev file
% filename = 'Chips_20171026_neuronSimulation';
% writeNEV(neuralData, opts.PACKET_WIDTH, filename, opts.MAP_FILE_NAME, '')
% save(strcat(filename,'_inputData'),'opts');

%% generate get spike times form nev file
NEV = openNEV(strcat(filename,'-s.nev'),'nosave');
spikeData = NEV.Data.Spikes;
spikeData.TimeStamp = spikeData.TimeStamp/30000;
%% process sorted spikes to get a percent found at each time step
% whole number is the artifact number in artifactData_neurons
% neuron idx stores the idx for that artifact number

foundStruct = getPercentFound(spikeData,artifactData_neurons,opts);

%% plot example artifacts with spikes at specified idxs
saveFigures = 0;
folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\Chips_20171026\SimulationStudy\';

plotFiltered = 1;
idxsPlot = 30:3:opts.IDX_ADD_NEURONS(end);
chan = 5;
figure();
% artifactIdxAll = [];
for idx = idxsPlot
    artifactIdx = find(artifactData_neurons.neuronIdx == idx);
    artifactIdx = datasample(artifactIdx,1,'replace',false);
%     artifactIdxAll = [artifactIdxAll;artifactIdx];
    xData = ((1:1:size(artifactData_neurons.artifact,3))-1)/30;
    if(plotFiltered)
        plot(xData,fliplr(filter(opts.B,opts.A,fliplr(squeeze(artifactData_neurons.artifact(artifactIdx,chan,:))')))','linewidth',1.5)
    else
        plot(xData,squeeze(artifactData_neurons.artifact(artifactIdx,chan,:))','linewidth',2)
    end
    
    hold on
end

ylim([-1000,1000])
xlim([0,3])
xlabel('Time after stimulation onset (ms)')
ylabel('Voltage (\muV)')
set(gca,'fontsize',16)
formatForLee(gcf)
if(saveFigures)
    if(~plotFiltered)
        saveFiguresLIB(gcf,folderpath,'Chips_20171026_exampleArtifactsWithSimulatedSpikes_chan50_stimChan50_unfiltered');
    else
        saveFiguresLIB(gcf,folderpath,'Chips_20171026_exampleArtifactsWithSimulatedSpikes_filtered');
    end
end
%% plot example artifacts and percent found for a single channel
colors = parula(size(foundStruct.percentCorrect,1));
for ch = 1:size(foundStruct.percentCorrect,1)
    figure
    ax(1)=subplot(2,1,1); % artifact plot
    plot(xData(1:60),squeeze(artifactData_neurons.artifact(1:10:100,ch,1:60))')
    ylabel('Voltage (\muV)')
    formatForLee(gcf)
    ax(2) = subplot(2,1,2); % percent found
    plot(opts.IDX_ADD_NEURONS/30,foundStruct.percentCorrect(ch,:),'linewidth',2,'color',colors(ch,:))
    linkaxes([ax(1),ax(2)],'x');
    ylabel('Percent recovered')
    xlabel('Time after stimulation onset (ms)')
    formatForLee(gcf)
    if(saveFigures)
        saveFiguresLIB(gcf,folderpath,strcat('Chips_20171026_artifactPercentRecovered_channelIdx',num2str(ch)));
    end
end



%% plot all percent correct 
colors = parula(size(foundStruct.percentCorrect,1));
figure
for c = 1:size(foundStruct.percentCorrect,1)
    plot(opts.IDX_ADD_NEURONS/30,foundStruct.percentCorrect(c,:),'linewidth',2,'color',colors(c,:))
    hold on
end
xlabel('Time after stimulation onset (ms)')
ylabel('Percent Correct')
set(gca,'fontsize',16)
formatForLee(gcf)

if(saveFigures)
    saveFiguresLIB(gcf,folderpath,'Chips_20171026_percentCorrect');
end
