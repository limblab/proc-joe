figure;

stimChan = 85;
art =1;
stimstep = 20;
% x = 1:1:numel(outputData.artifactData(art).artifact(stimChan,1,90:260));
% plot(x/30000*1000,1/1000*squeeze(outputData.artifactData(art).artifact(stimChan,1:stimstep:end,90:260)),'linewidth',1)
x = 1:1:numel(outputData.artifactData(art).artifact(stimChan,1,1:end));
plot(x/30000*1000,1/1000*squeeze(outputData.artifactData(art).artifact(stimChan,1:1:5,1:end)),'linewidth',1)

xlabel('Time (ms)')
ylabel('Voltage (mV)')

title('Chan 85 - high gain')
for i = 1:numel(squeeze(outputData.artifactData(art).artifact(stimChan,1:1:5,:)))
    legendStr{i} = strcat('Stimulation ', num2str((1) + stimStep*(i-1)));
end
l=legend(legendStr);
set(l,'box','off');
formatForLee(gcf);