    load('D:\Lab\Data\stim_ephys_paper\artifact_analysis\sim_study\Han_Duncan_20210120_duke_sim_plot');
    duke_bin_edges = bin_edges;
    duke_prop_rec = prop_recovered;
    load('D:\Lab\Data\stim_ephys_paper\artifact_analysis\sim_study\Han_Duncan_20210120_blackrock_sim_plot');
    black_bin_edges = bin_edges;
    black_prop_rec = prop_recovered;

%%

    f=figure(); hold on;

    amps_plot = [2,5,8,9];
%     amps_plot = [1:1:9];
    colors = inferno(numel(amps_plot)+1);

    % plot duke
    col_cnt = 1;
    for i_amp = amps_plot
        plot(duke_bin_edges(1:end-1)+mode(diff(duke_bin_edges))/2,duke_prop_rec(i_amp,:),...
            '-','color',colors(col_cnt,:),'linewidth',2,'markersize',16,'marker','.');

        % plot blackrock
%         plot(black_bin_edges(1:end-1)+mode(diff(black_bin_edges))/2,black_prop_rec(i_amp,:),...
%             '--','color',colors(col_cnt,:),'linewidth',2,'markersize',6,'marker','s');
        col_cnt = col_cnt + 1;
    end

    formatForLee(gcf);
    set(gca,'fontsize',14);
    xlabel('Time after stim offset (ms)');
    ylabel('Proportion spikes recovered');

    ylim([0,1]);
    xlim([0,5]);

%% statistics on this data -- load in some new ones
    load('D:\Lab\Data\stim_ephys_paper\artifact_analysis\sim_study\Han_Duncan_20210120_duke_sim_study_small_save');
    duke_spikes_found_ts = spikes_found_ts_list;
    duke_spikes_found_chan = spikes_found_chan_list;
    duke_spikes_total_ts = spikes_total_ts_list;
    duke_spikes_total_chan = spikes_total_chan_list;

    load('D:\Lab\Data\stim_ephys_paper\artifact_analysis\sim_study\Han_Duncan_20210120_blackrock_sim_study_small_save');
    black_spikes_found_ts = spikes_found_ts_list;
    black_spikes_found_chan = spikes_found_chan_list;
    black_spikes_total_ts = spikes_total_ts_list;
    black_spikes_total_chan = spikes_total_chan_list;

    amps = amp_list;
%% treat each channel as a repitition, and do stats comparing amplitude and duke or blackrock
    % get percent recovered in each bin for each channel and duke/blackrock
    bin_size = 0.1; % ms
    
    wave_length = 0.453;
    bin_edges = (0.25:bin_size:7);
    prop_recovered = nan(numel(duke_spikes_found_ts),numel(bin_edges)-1);
    
    
    bin_list = []; amp_list = []; is_duke = []; recov_list = []; chan_list = [];
    for i_amp = 1:numel(duke_spikes_found_ts)
        for i = 1:2
            switch i
                case 1
                    spikes_found_ts = duke_spikes_found_ts{i_amp}/30 - wave_length;
                    spikes_total_ts = duke_spikes_total_ts{i_amp}/30 - wave_length;
                    chan_found = duke_spikes_found_chan{i_amp};
                    chan_total = duke_spikes_total_chan{i_amp};
                case 2
                    spikes_found_ts = black_spikes_found_ts{i_amp}/30 - wave_length;
                    spikes_total_ts = black_spikes_total_ts{i_amp}/30 - wave_length;
                    chan_found = black_spikes_found_chan{i_amp};
                    chan_total = black_spikes_total_chan{i_amp};
            end
            
        % convert to time post stim offset (or onset), then bin and get
        % counts
            unique_chans = unique(chan_total);
            for i_chan = 1:numel(unique_chans)
                temp_found_counts = histcounts(spikes_found_ts(chan_found==i_chan),bin_edges);
                temp_total_counts = histcounts(spikes_total_ts(chan_total==i_chan),bin_edges);

                temp_recov = temp_found_counts./temp_total_counts;
                % store data in lists
                recov_list = [recov_list;temp_recov'];
                is_duke = [is_duke;(i==1)*ones(numel(temp_recov),1)];
                amp_list = [amp_list; amps(i_amp)*ones(numel(temp_recov),1)];
                bin_list = [bin_list; bin_edges(1:end-1)' + mode(diff(bin_edges))/2];
                chan_list = [chan_list; i_chan*ones(numel(temp_recov),1)];
            end
        
        end

    end
   
recov_tbl = table(bin_list,amp_list,categorical(is_duke),recov_list,categorical(chan_list),'VariableNames',{'t','amp','is_duke','recov','chan'});

recov_tbl(recov_tbl.recov <= 0,:) = [];

mdl_spec = 'recov~ amp*is_duke + t*is_duke + chan';
% mdl_spec = 'recov~t + amp + is_duke + chan + t:amp + t:is_duke + t:is_duke:amp';
recov_mdl = fitglm(recov_tbl,mdl_spec,'distribution','binomial','link','logit')


%% get time "recovered" to shade raster plots -- run previous two sections first
    unique_bins = unique(bin_list);

    mean_recov_duke = zeros(numel(amps),numel(unique_bins));
    mean_recov_black = zeros(size(mean_recov_duke));
    % smooth with bins before and after
    for i_amp = 1:numel(amps)
        for i_bin = 1:numel(unique_bins)
            bins_keep = [i_bin-1,i_bin,i_bin+1];
            bins_keep(bins_keep < 1 | bins_keep > numel(unique_bins)) = [];
            is_bin_mask = any(bin_list == unique_bins(bins_keep)',2);
            
            is_amp_mask = amp_list == amps(i_amp);
            is_duke_mask = is_duke;
            mean_recov_duke(i_amp,i_bin) = mean(recov_list(is_bin_mask & is_amp_mask & is_duke_mask));
            mean_recov_black(i_amp,i_bin) = mean(recov_list(is_bin_mask & is_amp_mask & ~is_duke_mask));

        end
    end
    
    prop_thresh = 0.5;
    duke_blank_times = zeros(size(black_prop_rec,1),2);
    black_blank_times = zeros(size(duke_blank_times));
    
    for i_amp = 1:size(black_prop_rec,1)
        duke_blank_times(i_amp,:) = [-0.453, unique_bins(find(mean_recov_duke(i_amp,:) > prop_thresh,1,'first'))];
        
        black_blank_times(i_amp,:) = [-0.453, unique_bins(find(mean_recov_black(i_amp,:) > prop_thresh,1,'first'))];
    end
    
    plot(duke_blank_times(:,2))















    