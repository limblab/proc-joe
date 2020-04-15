%% model response at different baseline firing rates

    window = [0,2]; % ms

    bin_idx = [find(ICMS_SinglePulse_Response{2,2} > window(1),1,'first'),...
        find(ICMS_SinglePulse_Response{2,2} > window(2),1,'first')];

    amp_list = [];
    response_amp = [];

    for a = 2:size(ICMS_SinglePulse_Response,1)
        amp_list(a-1,1) = ICMS_SinglePulse_Response{a,1};
        for baseline_idx = 3:5
            response_amp(a-1,baseline_idx-2) = sum(ICMS_SinglePulse_Response{a,baseline_idx}(bin_idx(1):bin_idx(2)));
        end
    end

    f=figure(); hold on
    f.Name = 'Model_singlePulse_baselineFreq';
    linestyles = {'--','-.',':'};
    for line_idx = 3:-1:1
        plot(amp_list,response_amp(:,line_idx),linestyles{line_idx},'linewidth',2,'color',getColorFromList(1,1))
    end
    formatForLee(gcf);
    set(gca,'fontsize',14)
    xlabel('Amplitude (\muA)');
    ylabel('Num Spikes/stims/neurons');
    l=legend('0Hz','10Hz','20Hz');
    set(l,'box','off','location','best')

%% experimental and model data plotted together (low amp version and a high amp version)

    baseline_idx = 3; % corresponds to 20 Hz
    window = [0,7];
    model_bin_idx = [find(ICMS_SinglePulse_Response{2,2} > window(1),1,'first'),...
        find(ICMS_SinglePulse_Response{2,2} > window(2),1,'first')];
    bin_centers = bin_edges(1:end-1) + mode(diff(bin_edges))/2;
    
    exp_bin_idx = [find(bin_centers > window(1),1,'first'),...
        find(bin_centers > window(2),1,'first')];
    
%     amps_plot = [5,10];
%     amps_plot = [20,30];
    amps_plot = [40,50];
    
%     amps_plot = [5,10,20];
%     amps_plot = [30,40,50];

    colors = inferno(numel(amps_plot)+1);
    
    model_line = '--'; model_marker = 's';
    model_offset = [-0.1,-0.01,0.08];
    exp_line = '-'; exp_marker = '.';
    exp_offset = [-0.08,0.01,0.1]; 
    
    % get baseline FR
    baseline_count = [];
    for unit = 1:numel(spikesStruct)
        baseline_count(unit) = spikesStruct{unit}.mean_baseline_fr*mode(diff(bin_edges))/1000;
    end
    mean_baseline_count = mean(baseline_count);
    
    % get std of bin counts
    std_bin_count = zeros(numel(bin_edges)-1,numel(amps_plot));
    
    for a = 1:numel(amps_plot)
        bin_count_amp = [];
        for unit = 1:numel(spikesStruct)
            amp_idx = find(amps_plot(a) == spikesStruct{unit}.amp,1,'first');
            if(~isempty(amp_idx))
                bin_count_amp(:,end+1) = histcounts(spikesStruct{unit}.spike_times_post_stim{amp_idx},bin_edges/1000)/arrayData{unit}.numStims(amp_idx);
            end
        end
        std_bin_count(:,a) = std(bin_count_amp,0,2);%/sqrt(num_responsive(amp_idx));
    end
    
    %%
    
    figure(); hold on
    color_counter = 1;
    for model_amp_idx = 2:size(ICMS_SinglePulse_Response,1)
        if(any(ICMS_SinglePulse_Response{model_amp_idx,1} == amps_plot))
            exp_amp_idx = find(amp_list == ICMS_SinglePulse_Response{model_amp_idx,1});
            
            model_x = ICMS_SinglePulse_Response{model_amp_idx,2}(model_bin_idx(1):model_bin_idx(2)) + model_offset(color_counter);
            model_y = ICMS_SinglePulse_Response{model_amp_idx,baseline_idx}(model_bin_idx(1):model_bin_idx(2)) - 0.02;
            exp_x = bin_centers(exp_bin_idx(1):exp_bin_idx(2)) + exp_offset(color_counter);
            exp_y = bin_counts{exp_amp_idx}(exp_bin_idx(1):exp_bin_idx(2))/num_responsive(exp_amp_idx) - mean_baseline_count;
            exp_bars = std_bin_count(exp_bin_idx(1):exp_bin_idx(2),color_counter);
            
            errorbar(exp_x,exp_y,exp_bars,'linestyle',exp_line,'marker',exp_marker,'markersize',26,...
                'linewidth',1,'color',colors(color_counter,:));
            plot(model_x,model_y,'linestyle',model_line,'marker',model_marker,'markersize',12,...
                'linewidth',1,'color',colors(color_counter,:));
       
            color_counter = color_counter + 1;
        end
    end

    formatForLee(gcf)
    set(gca,'fontsize',14)
    xlabel('Time post stim (ms)');
    ylabel('Num spikes/stims/neurons');
        
%% plot model and experiment counts in certain bins
% plot experiment in [1,3], model in [0,1],[1,3],[0,4]

    baseline_idx = 3; % corresponds to 20 Hz
%     exp_window = [1,3];
%     model_window = [1,3; 0,1; 0,3];
    exp_window = [3,7];
    model_window = [3,7];
    
    exp_color = getColorFromList(1,1);
    model_color = [getColorFromList(1,1); getColorFromList(1,0); getColorFromList(1,2)];
    exp_markers = {'.'}; exp_marker_size = 26;
    model_markers = {'s','s','s'}; model_marker_size = [10,6,10];
    
    amps_plot = [5,10,15,20,25,30,40,50,100];
    
    % get baseline FR
    baseline_count = [];
    for unit = 1:numel(spikesStruct)
        baseline_count(unit) = spikesStruct{unit}.mean_baseline_fr*mode(diff(bin_edges))/1000;
    end
    mean_baseline_count = mean(baseline_count);
    
    

    
    bin_centers = bin_edges(1:end-1) + mode(diff(bin_edges))/2;
   
    
    figure(); hold on
    % experiment window first
    for i_window = 1:size(exp_window,1)
        exp_bin_idx = [find(bin_centers > exp_window(i_window,1),1,'first'),...
            find(bin_centers > exp_window(i_window,2),1,'first')-1];

        % get std of bin counts
        std_bin_count = zeros(1,numel(amps_plot));

        for a = 1:numel(amps_plot)
            bin_count_amp = [];
            for unit = 1:numel(spikesStruct)
                amp_idx = find(amps_plot(a) == spikesStruct{unit}.amp,1,'first');
                if(~isempty(amp_idx))
                    bin_count_amp(:,end+1) = sum(histcounts(spikesStruct{unit}.spike_times_post_stim{amp_idx},bin_edges(exp_bin_idx(1):exp_bin_idx(2)+1)/1000)/arrayData{unit}.numStims(amp_idx));
                end
            end
            std_bin_count(a) = std(bin_count_amp);%/sqrt(num_responsive(amp_idx));
        end
        
        for i_amp = 1:numel(amps_plot)
            exp_x = amp_list(i_amp)-0.5;
            exp_y = sum(bin_counts{i_amp}(exp_bin_idx(1):exp_bin_idx(2))/num_responsive(i_amp) - mean_baseline_count);
            exp_bars = std_bin_count(i_amp);
            
            errorbar(exp_x,exp_y,exp_bars,'linestyle','-','marker',exp_markers{i_window},'markersize',exp_marker_size(i_window),...
                'linewidth',1,'color',exp_color(i_window,:)); hold on
        end
    end
    
    % model windows next
    for i_window = 1:size(model_window,1)
        model_bin_idx = [find(ICMS_SinglePulse_Response{2,2} > model_window(i_window,1),1,'first'),...
        find(ICMS_SinglePulse_Response{2,2} > model_window(i_window,2),1,'first')-1];
    
        for i_amp = 1:numel(amps_plot)
            model_amp_idx = find([ICMS_SinglePulse_Response{2:end,1}] == amps_plot(i_amp)) + 1;
            
            model_x = amp_list(i_amp)+0.5;
            model_y = sum(ICMS_SinglePulse_Response{model_amp_idx,baseline_idx}(model_bin_idx(1):model_bin_idx(2)) - 0.01);
            exp_bars = std_bin_count(i_amp);
            
            plot(model_x,model_y,'linestyle','-','marker',model_markers{i_window},'markersize',model_marker_size(i_window),...
                'linewidth',1,'color',model_color(i_window,:)); hold on
        end
    end

    formatForLee(gcf); 
    xlabel('Amplitude (\muA)');
    ylabel('Num spikes/stim/neuron');
    set(gca,'fontsize',14)
