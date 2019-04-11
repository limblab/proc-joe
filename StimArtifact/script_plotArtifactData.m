%% get artifact data in a nice format
% from outputData files through StimRecProcessing scripts
% also need inputData

    waveSent = 1;
    chanSent = [6];
    chanRec = 1;
    
    artifactData = outputData.artifactData;
    waveforms = outputData.waveforms;

%     artifact_idx = find(waveforms.waveSent == waveSent & checkChanListEquality(waveforms.chanSent, chanSent));
    artifact_idx = find(ones(size(outputData.waveforms.chanSent))==1);
    % x_data is in ms and from stimulation offset
    x_data = ((1:size(artifactData.artifact,3))-1-inputData.presample)/30;% - ...
%         waveforms.parameters(waveSent).pWidth1/1000 - waveforms.parameters(waveSent).pWidth2/1000 - 53/1000;

    % plot artifact and filtered artifact
    num_plot = 10;
    template = median(squeeze(artifactData.artifact(artifact_idx,chanRec,:)));
    % plot raw artifact
    figure()
    ax1=subplot(2,1,1);
    [~,idx_use] = datasample(artifact_idx,min(num_plot,numel(artifact_idx)),'Replace',false);
    plot(x_data,squeeze(artifactData.artifact(idx_use,chanRec,:))')
    xlim([-0.3,5]);
    formatForLee(gcf)
    ylabel('Voltage (\muV)');
    set(gca,'fontsize',14)

    % plot artifact after acausal filtering
    ax2=subplot(2,1,2);
    plot(x_data,acausalFilter(squeeze(artifactData.artifact(idx_use,chanRec,:))'-template'))
    xlim([-0.3,5]);
    ylim([-250,250])
    formatForLee(gcf)
    xlabel('Time after stimulation offset (ms)');
    ylabel('Voltage (\muV)');
    set(gca,'fontsize',14)
    linkaxes([ax1,ax2],'x');
    
    
%% simulate neurons on top of the artifact (demonstration)
% load in desired neuron type (neuronMeanWave is used)
    mean_latency = 0.5; % ms
    std_latency = 0.2; % ms
    percent_with_spikes = 0.5; % percentage of artifacts with spikes
    base_artifact_idx = 50;
    
    % get pretty single artifact example
    base_artifact = squeeze(artifactData.artifact(artifact_idx(base_artifact_idx),chanRec,:));
    time_after_add = [0.3:0.2:1.1,2]; % in ms
    idx_add = floor(time_after_add*30) + inputData.presample + ...
        floor((waveforms.parameters(waveSent).pWidth1/1000 + waveforms.parameters(waveSent).pWidth2/1000 + 53/1000)*30);
    
    simulated_artifact_single = zeros(size(base_artifact,1),numel(idx_add));
    for i = 1:numel(idx_add)
        simulated_artifact_single(:,i) = base_artifact;
        simulated_artifact_single(idx_add(i):idx_add(i)+numel(neuronMeanWave)-1,i) = ...
            base_artifact(idx_add(i):idx_add(i)+numel(neuronMeanWave)-1) + neuronMeanWave'*0.5;
        
    end
    
    % add fake neurons to artifactData.artifact before building template to
    % simulate an actual process (assume normal distribution)
    simulated_artifact_all = squeeze(artifactData.artifact(artifact_idx,chanRec,:))';
    artifacts_with_spikes = randperm(size(simulated_artifact_all,2),floor(size(simulated_artifact_all,2)*percent_with_spikes));
    
    spike_latency = normrnd(mean_latency,std_latency,size(artifacts_with_spikes));
    % resample spike latencies below 0
    counter = 0;
    while(sum(spike_latency < 0) > 0 && counter < 1000)
        if(mod(counter,10) == 0)
            warning('spike latencies below 0 detected, resampling');
        end
        spike_latency(spike_latency < 0) = normrnd(mean_latency,std_latency,sum(spike_latency < 0),1);
    end
    
    % place spikes into simulated_artifact_all
    for i = 1:numel(artifacts_with_spikes)
        idx_latency = floor(spike_latency(i)*30 + inputData.presample + (waveforms.parameters(waveSent).pWidth1/1000 ...
            + waveforms.parameters(waveSent).pWidth2/1000 + 53/1000)*30);
        simulated_artifact_all(idx_latency:idx_latency+numel(neuronMeanWave)-1,artifacts_with_spikes(i)) = ...
            simulated_artifact_all(idx_latency:idx_latency+numel(neuronMeanWave)-1,artifacts_with_spikes(i)) + neuronMeanWave'*0.5;
        
    end
    
    % build median template for each artifact?
%     for i = 1:size(simulated_artifact_all,2)
%         
%         
%     end
    template = median(simulated_artifact_all,2);

    
    % plot nice example without template subtraction to show problem. Then
    % plot real examples after template subtraction with template
    figure();
    ax1=subplot(2,1,1)
    plot(x_data,acausalFilter(simulated_artifact_single))
    hold on
    plot(x_data,acausalFilter(base_artifact),'k','linewidth',2)
    xlim([-0.3,3]);
    ylim([-1000,1000])
    formatForLee(gcf)
    ylabel('Voltage (\muV)');
    set(gca,'fontsize',14)
    
%     subplot(2,1,2)
%     plot(x_data,acausalFilter(simulated_artifact_single - template));
%     xlim([-0.3,3]);
%     ylim([-1000,1000])
%     formatForLee(gcf)
%     ylabel('Voltage (\muV)');
%     set(gca,'fontsize',14)
    
    ax2=subplot(2,1,2)
    plot(x_data,acausalFilter(simulated_artifact_all(:,artifacts_with_spikes(1:10:end))-template))
    hold on
    plot(x_data,acausalFilter(template),'k','linewidth',2)
    xlim([-0.3,3]);
    ylim([-1000,1000])
    formatForLee(gcf)
    ylabel('Voltage (\muV)');
    set(gca,'fontsize',14)
    
    linkaxes([ax1,ax2],'xy');
    
%% full simulation given a method (lol)
% this will randomly give neurons to artifacts at random times, then 






































