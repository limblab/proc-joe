%% get artifact data in a nice format and get metrics for each stimulation chan and recording chan
% from outputData files through StimRecProcessing scripts
% also need inputData
    artifactData = outputData.artifactData;
    waveforms = outputData.waveforms;
    
    unique_chan_list = uniquecell(waveforms.chanSent);
    mean_time_off_rails = zeros(numel(unique_chan_list),numel(waveforms.parameters),96);
    for chan_sent_idx = 1:numel(unique_chan_list)
        for wave_sent = 1:numel(waveforms.parameters)
            artifact_idx = find(waveforms.waveSent == wave_sent & ...
                checkChanListEquality(waveforms.chanSent, unique_chan_list{chan_sent_idx}));

            artifact = artifactData.artifact(artifact_idx,:,:);
            
            [metrics_out] = getArtifactMetrics(artifact);
            mean_time_off_rails(chan_sent_idx,wave_sent,:)  = mean(metrics_out.idx_off_rails - inputData.presample)/30; % in ms
        end
    end

%% for each recording channel, plot time off rails
    f=figure();
    f.Name = '20190218_timeOffRails';

    ax = gca;
    for chan = 1:size(mean_time_off_rails,3)
        for wave = 1%:size(mean_time_off_rails,2)
            plot(chan,squeeze(mean_time_off_rails(:,wave,chan)),'.','markersize',16,'color',getColorFromList(1,wave))
            hold on
        end
    end
    ax.YLim(1) = 0;
    formatForLee(gcf);
    xlabel('Elec');
    ylabel('Time off rails (ms)');
    set(gca,'fontsize',14);
%% make artifact_keep_mask
    artifact_keep_mask = ones(size(mean_time_off_rails,3),1);
    for chan = 1:size(mean_time_off_rails,3)
        if(any(squeeze(mean_time_off_rails(:,end,chan)) > 1.75))
            artifact_keep_mask(chan) = 0;
        end
    end

