%% plot artifact data
% from outputData files through StimRecProcessing scripts
% also need inputData

    waveSent = 1;
    chanSent = 36;

    artifactData = outputData.artifactData;
    waveforms = outputData.waveforms;

    artifact_idx = find(waveforms.waveSent == waveSent & checkChanListEquality(waveforms.chanSent, chanSent));

    % x_data is in ms and from stimulation offset
    x_data = ((1:size(artifactData.artifact,3))-1-inputData.presample)/30 - ...
        waveforms.parameters(waveSent).pWidth1/1000 - waveforms.parameters(waveSent).pWidth2/1000 - 53/1000;

    
    % plot raw artifact
    figure()
    subplot(2,1,1)
    plot(x_data,squeeze(artifactData.artifact(artifact_idx(1:10:end),chanSent,:))')
    xlim([-0.3,3]);
    formatForLee(gcf)
    ylabel('Voltage (\muV)');
    set(gca,'fontsize',14)
    
    % plot artifact after acausal filtering
    subplot(2,1,2)
    plot(x_data,acausalFilter(squeeze(artifactData.artifact(artifact_idx(1:10:end),chanSent,:))'))
    xlim([-0.3,3]);
    formatForLee(gcf)
    xlabel('Time after stimulation offset (ms)');
    ylabel('Voltage (\muV)');
    set(gca,'fontsize',14)
    
    





































