load('D:\Lab\Data\stim_ephys_paper\artifact_analysis\sim_study\Han_Duncan_20210120_duke_sim_plot');
duke_bin_edges = bin_edges;
duke_prop_rec = prop_recovered;
load('D:\Lab\Data\stim_ephys_paper\artifact_analysis\sim_study\Han_Duncan_20210120_blackrock_sim_plot');
black_bin_edges = bin_edges;
black_prop_rec = prop_recovered;

%%

f=figure(); hold on;

% amps_plot = [2,5,8,9];
amps_plot = [1:1:9];
colors = inferno(numel(amps_plot)+1);

% plot duke
col_cnt = 1;
for i_amp = amps_plot
%     plot(duke_bin_edges(1:end-1)+mode(diff(duke_bin_edges))/2,duke_prop_rec(i_amp,:),...
%         '-','color',colors(col_cnt,:),'linewidth',2,'markersize',16,'marker','.');

    % plot blackrock
    plot(black_bin_edges(1:end-1)+mode(diff(black_bin_edges))/2,black_prop_rec(i_amp,:),...
        '--','color',colors(col_cnt,:),'linewidth',2,'markersize',6,'marker','s');
    col_cnt = col_cnt + 1;
end

formatForLee(gcf);
set(gca,'fontsize',14);
xlabel('Time after stim offset (ms)');
ylabel('Proportion spikes recovered');

ylim([0,1]);
xlim([0,5]);



