figure()
spike_times_mean = spike_times_use;
plot((bin_edges(1:end-1) + mode(diff(bin_edges)))*1000,spike_times_mean,'k','linewidth',2)
xlabel('Time after stimulation offset (ms)');
ylabel('Total spike count');
formatForLee(gcf);
set(gca,'fontsize',14)
xlim([-7.5,7.5])