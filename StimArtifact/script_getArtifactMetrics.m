%% get artifact data in a nice format and get metrics for each stimulation chan and recording chan
% from outputData files through StimRecProcessing scripts
% also need inputData
    artifactData = outputData.artifactData;
    waveforms = outputData.waveforms;
    
    unique_chan_list = uniquecell(waveforms.chanSent);
    single_chan_mask = cellfun(@numel,unique_chan_list) == 1;
    mean_time_off_rails = zeros(numel(unique_chan_list),numel(waveforms.parameters),96);
    for chan_sent_idx = 1:numel(unique_chan_list)
        for wave_sent = 1:numel(waveforms.parameters)
            artifact_idx = find(waveforms.waveSent == wave_sent & ...
                checkChanListEquality(waveforms.chanSent, unique_chan_list{chan_sent_idx}));

            artifact = artifactData.artifact(artifact_idx,:,:);
            
            [metrics_out] = getArtifactMetrics(artifact);
            mean_time_off_rails(chan_sent_idx,wave_sent,:)  = mean(metrics_out.idx_off_rails - inputData.presample,'omitnan')/30; % in ms
        end
    end

    mean_time_off_rails = mean_time_off_rails(single_chan_mask==1,:,:);
%% for each recording channel, plot time off rails
    f=figure();
    f.Name = '20190218_timeOffRails';
    map_file_name = 'R:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\left S1 20190205\SN 6251-002087.cmp';      
%     map_file_name = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
    map_data = loadMapFile(map_file_name);
    
    % plot time off rails for each stim electrode against distance from stim
    % chan
    ax=gca;
    for chan_stim_idx = 1:size(mean_time_off_rails,1)
        for chan_rec_idx = 1:size(mean_time_off_rails,3)
            for wave = 1%:size(mean_time_off_rails,2)
                map_data_chan_stim = find(map_data.chan == unique_chan_list{chan_stim_idx});
                map_data_chan_rec = find(map_data.chan == chan_rec_idx);
                dist = 400*sqrt((map_data.row(map_data_chan_stim)-map_data.row(map_data_chan_rec))^2 + (map_data.col(map_data_chan_stim)-map_data.col(map_data_chan_rec))^2);
                plot(dist,squeeze(mean_time_off_rails(chan_stim_idx,wave,chan_rec_idx)),'.','markersize',16,'color',getColorFromList(1,chan_stim_idx))
                hold on
            end
        end
    end
    ax.YLim(1) = 0;
    formatForLee(gcf);
    xlabel('Distance (\mum)');
    ylabel('Time off rails (ms)');
    set(gca,'fontsize',14);
%% make artifact_keep_mask
    artifact_keep_mask = ones(size(mean_time_off_rails,3),1);
    for chan = 1:size(mean_time_off_rails,3)
        if(any(squeeze(mean_time_off_rails(:,end,chan)) > 1.75))
            artifact_keep_mask(chan) = 0;
        end
    end

%% put mean_time_off_rails into arrayData

    for arr_idx = 1:numel(arrayData)
        arrayData{arr_idx}.mean_time_off_rail = squeeze(mean_time_off_rails(:,:,arrayData{arr_idx}.CHAN_REC));
    end