function [  ] = plotAverageWaveform(cds,neuronNumber)
%


spikesTable = cds.units(neuronNumber).spikes;
wavesWave = spikesTable{:,:};
if(size(wavesWave,1)>1)
    aveWavesWave = mean(wavesWave);
else
    aveWavesWave = wavesWave;
end
   
plot(aveWavesWave,'k','linewidth',2);
hold on
plot(aveWavesWave + std(wavesWave),'--k','linewidth',2);
plot(aveWavesWave - std(wavesWave),'--k','linewidth',2);
ylabel('Voltage (\muV)');
xlabel('Wave Point');
formatForLee(gcf);


end

