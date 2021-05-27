%% get experimental data
    exp_input_data.home_computer = 1;
    [exp_data] = getExperimentDoublePulseData(exp_input_data);
    exp_array_data = exp_data.array_data;
    
    % sanitize experimental array data because field names and shape are
    % not consistent
    for i_unit = 1:numel(exp_array_data)
        exp_array_data{i_unit}.spikeTrialTimes = exp_array_data{i_unit}.spikeTrialTimes';
        exp_array_data{i_unit}.trial_num = exp_array_data{i_unit}.trial_num';
        exp_array_data{i_unit}.numStims = exp_array_data{i_unit}.num_stims';
        exp_array_data{i_unit}.binCounts = exp_array_data{i_unit}.binCounts';
        exp_array_data{i_unit}.binEdges = exp_array_data{i_unit}.binEdges';
        exp_array_data{i_unit}.stimData = exp_array_data{i_unit}.trial_num;
    end
    
    % get baseline FR and adjust spike timing
    exp_array_data = adjustArrayDataSpikeTimes(exp_array_data, 0.453/1000); % stim pulse length
    exp_array_data = getBaselineFiringRate(exp_array_data,[-80,-5]/1000); % window relative to stim onset


% santize array data
for i_unit = 1:numel(exp_array_data)
    exp_array_data{i_unit}.spikeTrialTimes = exp_array_data{i_unit}.spikeTrialTimes';
    exp_array_data{i_unit}.trial_num = exp_array_data{i_unit}.trial_num';
    exp_array_data{i_unit}.num_stims = exp_array_data{i_unit}.num_stims';
    exp_array_data{i_unit}.binCounts = exp_array_data{i_unit}.binCounts';
    exp_array_data{i_unit}.binEdges = exp_array_data{i_unit}.binEdges';
    exp_array_data{i_unit}.stimData = exp_array_data{i_unit}.stimData';
end

%%
num_stims = [];
idx_keep = [3:10,28:31];
for i=idx_keep
    num_stims = [num_stims,exp_array_data{i}.numStims([2,3,4,5])];
end


%% plot raster of example neuron(s)

    raster_input_data = [];
    raster_input_data.x_lim = [-100,480]; % ms
    raster_input_data.cond_list = [2,3,4,5]; % 
    raster_input_data.marker_style = '.'; % line is the correct way, but much slower
    
    for exp_idx = 1:numel(exp_array_data) % 6 was used for paper figures
        raster_input_data.is_model = 0;
        if(numel(exp_array_data{exp_idx}.num_stims) == 8)
            f=plotModelExpDoublePulseRaster(exp_array_data{exp_idx},raster_input_data);
            f.Name = num2str(exp_idx);
        end
    end

    

    
%% rebound excitation stats 
% duration, percent of cells
    rebound_input_data.cond_list = [2,3,4,5];
    rebound_input_data.bin_window = [-200,500]./1000;
    rebound_input_data.pre_window = [-80,-10]./1000;
    rebound_input_data.post_window = [10,250]./1000;
    rebound_input_data.blank_time = [-10,10]/1000; %s
    
    rebound_input_data.bin_size = 5/1000; % s
    rebound_input_data.kernel_length = 1;
    rebound_input_data.num_consec_bins = 6;
    [rebound_data] = getReboundExcitationWrapper(exp_array_data,rebound_input_data);
    
    
%% plot rebound excitation data -- 

    % plot duration across frequencies (with lines between same cells)
    f=figure('Position',[1403 522 395 420]);
    f.Name = 'Han_duncan_short_trains_rebound_excitation';
    suptitle('Short Trains')
    subplot(2,1,2); hold on;
    freq_list=[179,94,49,20];
    for i_cell = 1:size(rebound_data.is_rebound,1)
        plot(freq_list+7*rand()-3.5,rebound_data.rebound_dur(i_cell,:)*1000,'-k','marker','o','markersize',6);
    end
    xlim([0,200])
    ylim([0,300]);
    formatForLee(gcf);
    xlabel('Frequency (Hz)');
    ylabel('Duration (ms)');
    set(gca,'fontsize',14);
    ax=gca;
    ax.XTick = sort(freq_list);
    ax.XMinorTick = 'off';


    % also plot % of cells with rebound across frequencies
    subplot(2,1,1); hold on;
    frac_data = nan(size(freq_list));
    for i_freq = 1:size(rebound_data.is_rebound,2)
        frac_data(i_freq) = sum(rebound_data.is_rebound(:,i_freq),'omitnan')/sum(~isnan(rebound_data.is_rebound(:,i_freq)));
    end
    
    bar(freq_list,frac_data,'EdgeColor','k','FaceColor','k');
    xlim([0,200])
    ylim([0,1]);
    formatForLee(gcf);
    xlabel('Frequency (Hz)');
    ylabel('Fraction cells');
    set(gca,'fontsize',14);
    ax=gca;
    ax.XTick = sort(freq_list);
    ax.XMinorTick = 'off';
    
    
    
%% compute inhibition duration and plot across train frequencies

    inhib_input_data = [];
    inhib_input_data.pre_window = [-80,-10]/1000; % s 
    inhib_input_data.post_window = [0,260]/1000; % s
    inhib_input_data.bin_window = [inhib_input_data.pre_window(1),490/1000];
    inhib_input_data.max_time_start = 150/1000; % s
    inhib_input_data.bin_size = 5/1000; % s
    inhib_input_data.kernel_length = 1;
    inhib_input_data.blank_time = [-5,10]/1000; % s
    
    inhib_input_data.num_consec_bins = 2;
    inhib_input_data.cond_list = [1,2,3,4,5]; % single, 180Hz, 90 Hz, 50Hz, 20 Hz
    
    exp_inhib_data = getInhibitionDurationDoublePulseWrapper(exp_array_data,inhib_input_data);

%% plot inhibition data for short trains

    % plot inhibition duration across frequencies (with lines between same cells)
    f=figure('Position',[1403 522 395 420]);
    f.Name = 'Han_duncan_short_trains_inhibition';
    suptitle('Short Trains')
    subplot(2,1,2); hold on;
    freq_list=[179,94,49,20];
    for i_cell = 1:size(exp_inhib_data.is_inhib,1)
        plot(freq_list+7*rand()-3.5,exp_inhib_data.inhib_dur(i_cell,2:end)*1000,'-k','marker','o','markersize',6);
    end
    xlim([0,200])
    ylim([0,300]);
    formatForLee(gcf);
    xlabel('Frequency (Hz)');
    ylabel('Duration (ms)');
    set(gca,'fontsize',14);
    ax=gca;
    ax.XTick = sort(freq_list);
    ax.XMinorTick = 'off';

    % also plot % of cells with inhibition response
    subplot(2,1,1); hold on;
    frac_data = nan(size(freq_list)-1);
    for i_freq = 2:size(exp_inhib_data.is_inhib,2)
        frac_data(i_freq-1) = sum(exp_inhib_data.is_inhib(:,i_freq),'omitnan')/sum(~isnan(exp_inhib_data.is_inhib(:,i_freq)));
    end
    
    bar(freq_list,frac_data,'EdgeColor','k','FaceColor','k');
    xlim([0,200])
    ylim([0,1]);
    formatForLee(gcf);
    xlabel('Frequency (Hz)');
    ylabel('Fraction cells');
    set(gca,'fontsize',14);
    ax=gca;
    ax.XTick = sort(freq_list);
    ax.XMinorTick = 'off';
    
%% inhib dur stats
    freq_list=[179,94,49,20];
    
    freqs = []; inhib_durs = []; cell_id = [];
    for i_cell = 1:size(exp_inhib_data.is_inhib,1)
        num_add = sum(exp_inhib_data.is_inhib(i_cell,2:end));
        cell_id = [cell_id; i_cell+zeros(num_add,1)];
        freqs = [freqs; freq_list(exp_inhib_data.is_inhib(i_cell,2:end)==1)'];
        inhib_durs = [inhib_durs; exp_inhib_data.inhib_dur(i_cell,find(exp_inhib_data.is_inhib(i_cell,2:end))+1)'];
    end
    inhib_durs = inhib_durs*1000;
    
    inhib_tbl = table(freqs, inhib_durs, categorical(cell_id), 'VariableNames',{'freq','dur','cell'});
    mdl_spec = 'dur~freq + cell';
    
    inhib_mdl = fitlm(inhib_tbl,mdl_spec)
    
%% percent of baseline blocked during stimulus
    blocking_input_data.cond_list = [2,3,4,5];
    blocking_input_data.blank_time = [5]/1000; %s
    blocking_input_data.max_train_time = 0.2; % s
    
    [rebound_data] = getBaselineActivityBlocked(exp_array_data,blocking_input_data);


    
%% get decay rate for each neuron across stim frequencies
    decay_input_data.cond_list = [2,3,4,5];
    decay_input_data.spike_window = [1,7]; % s,
    [decay_data] = getDecayShortTrains(exp_array_data,decay_input_data);
    
    
%% plot decay rate and rebound data together to see if there is a relationship

    figure();
    boxplot_params = [];
    boxplot_params.use_same_color_for_all=1;
    for i_freq = 1%size(rebound_data.is_rebound,2)
        boxplot_params.master_color = getColorFromList(1,0);
        boxplot_wrapper(0,decay_data.slope(rebound_data.is_rebound(:,i_freq)==1,i_freq),boxplot_params)
        boxplot_params.master_color = getColorFromList(1,1);
        boxplot_wrapper(1,decay_data.slope(rebound_data.is_rebound(:,i_freq)==0,i_freq),boxplot_params)
    end


%% sanity check code
figure(); hold on
idx = 27;

for i=1:size(exp_inhib_data.inhib_dur,2)
    x_data = inhib_input_data.bin_window(1):inhib_input_data.bin_size:inhib_input_data.bin_window(2);
    x_data = x_data(1:end-1) + mode(diff(x_data))/2;
    plot(x_data,squeeze(exp_inhib_data.filtered_PSTH(idx,i,:))','color',getColorFromList(1,i-1),'linewidth',1.5)
end
plot([x_data(1),x_data(end)],exp_inhib_data.threshold(idx,1)+[0,0],'k--','linewidth',1)

exp_inhib_data.inhib_dur(idx,:)