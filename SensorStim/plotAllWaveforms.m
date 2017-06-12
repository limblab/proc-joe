function [  ] = plotAllWaveforms( cds,neuronNumber )
%
spikesTable = cds.units(neuronNumber).spikes;
wavesWave = spikesTable{:,:};
  
plot(wavesWave','k','linewidth',1);
ylabel('Voltage (\muV)');
xlabel('Wave Point');
formatForLee(gcf);

end

