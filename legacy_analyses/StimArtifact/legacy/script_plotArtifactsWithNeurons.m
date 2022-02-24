%% load in artifact data and cds to get spike times
    inputData.folderpath = 'C:\Users\jts3256\Desktop\Han_stim_data\Han_20190821_stimrec\chan60\';
    inputData.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
%     inputData.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\left S1 20190205\SN 6251-002087.cmp';


    folderpath = inputData.folderpath; % rest of code uses folderpath currently...may have switched this, not 100% certain

    inputData.task='taskCObump';
    inputData.ranBy='ranByJoseph'; 
    inputData.array1='arrayLeftS1'; 
    inputData.monkey='monkeyHan';
    inputData.labnum = 6;

    pwd=cd;
    cd(inputData.folderpath)
    fileList = dirSorted('*spikesExtracted.nev*');
    
    fileNumber = 1;
    
    cds = commonDataStructure();
    cd(inputData.folderpath)
    cds.file2cds([inputData.folderpath fileList(fileNumber).name],inputData.task,inputData.ranBy,...
        inputData.monkey,inputData.labnum,inputData.array1,inputData.mapFileName); % DO NOT USE RECOVER PRE SYNC, currently this shifts the units and the analog signal differently

    inputDataFileList = dir('*inputData*');
    outputDataFileList = dir('*outputData*');
    
    load(inputDataFileList(1).name);
    load(outputDataFileList(1).name);

%% make plot with artifacts with neurons and artifacts without neurons

    chan_rec = 60;
    unit_idx = find([cds.units.chan] == chan_rec & [cds.units.ID] == 1);
    wave_idx = 8;
    
    window = [0,5]/1000; % s
    
    
    artifactData = outputData.artifactData;
    waveforms = outputData.waveforms;

    % get artifacts with spikes after
    artifact_spike_idx = [];
    spike_idx = [];
    for a = 1:numel(artifactData.t)
        spike_mask = cds.units(unit_idx).spikes.ts > artifactData.t(a) + window(1) & cds.units(unit_idx).spikes.ts < artifactData.t(a) + window(2);
        if(sum(spike_mask) > 0 && ~isempty(wave_idx) && outputData.waveforms.waveSent(a) == wave_idx)
            artifact_spike_idx = [artifact_spike_idx; a];
            spike_idx = [spike_idx;find(spike_mask)];
        end
    end
    
    artifact_no_spike_idx = setdiff(1:numel(artifactData.t),artifact_spike_idx);
    
    %
    num_plot = 10;
    plot_filtered = 1;
    y_limits = 250*[-1,1];
    x_data = ((1:size(artifactData.artifact,3))-1-inputData.presample)/30;% - ...

    % plot artifact with spikes afterwards
    figure()
    ax1=subplot(2,1,1);
    [~,idx_use] = datasample(artifact_spike_idx,min(num_plot,numel(artifact_spike_idx)),'Replace',false);
    artifact_chan = chan_rec;
    if(artifact_chan > size(artifactData.artifact,2)) % in case I only dealt with the stimulated channel
        artifact_chan = 1;
    end
    if(plot_filtered)
        plot(x_data,acausalFilter(squeeze(artifactData.artifact(idx_use,artifact_chan,:))'))
    else
        plot(x_data,squeeze(artifactData.artifact(idx_use,artifact_chan,:))');
    end
    xlim([-0.3,5]);
    formatForLee(gcf)
    ylabel('Voltage (\muV)');
    set(gca,'fontsize',14)
    ylim(y_limits);
    
    % plot artifact withOUT spikes afterwards
    ax1=subplot(2,1,2);
    [~,idx_use] = datasample(artifact_no_spike_idx,min(num_plot,numel(artifact_no_spike_idx)),'Replace',false);
    if(plot_filtered)
        plot(x_data,acausalFilter(squeeze(artifactData.artifact(idx_use,artifact_chan,:))'))
    else
        plot(x_data,squeeze(artifactData.artifact(idx_use,artifact_chan,:))');
    end
    xlim([-0.3,5]);
    formatForLee(gcf)
    ylabel('Voltage (\muV)');
    set(gca,'fontsize',14)
    ylim(y_limits);

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    