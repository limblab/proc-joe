%% make array data using analyzeStimData or load one in here

%% for each unit, get number of spikes in a window and then match pairs of electrodes
    map_file = 'R:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\right S1 20180919\SN 6251-001804.cmp';
    map_data = loadMapFile(map_file);
    
    spike_window_time = [1,5]; % in ms
    baseline_window_time = [-10,-1]; % in ms
    
    num_spikes_above_baseline = {};
    for arr_idx = 1:numel(arrayData)
        % collapse spiking data into 1 number (number of spikes per stimulation
        % in a provided window (spike_window_time)
        spikes_window_idx = [find(spike_window_time(1) > arrayData{1}.bE{1,1},1,'last') ...
            find(spike_window_time(2) > arrayData{1}.bE{1,1},1,'last')];
        baseline_window_idx = [find(baseline_window_time(1) > arrayData{1}.bE{1,1},1,'last') ...
            find(baseline_window_time(2) > arrayData{1}.bE{1,1},1,'last')];
        num_spikes_above_baseline{arr_idx} = zeros(size(arrayData{arr_idx}.stimData));
        for i = 1:size(num_spikes_above_baseline{arr_idx},1)
            for j = 1:size(num_spikes_above_baseline{arr_idx},2)
                num_spikes_above_baseline{arr_idx}(i,j) = sum(arrayData{arr_idx}.bC{i,j}(spikes_window_idx(1):spikes_window_idx(2))) - ...
                    mean(arrayData{arr_idx}.bC{i,j}(baseline_window_idx(1):baseline_window_idx(2)))*(spikes_window_idx(2)-spikes_window_idx(1)+1);
            end
        end

    end

% for each unit, match electrode pairs and individual electrode responses
    match_data = {};
    
    for arr_idx = 1:numel(num_spikes_above_baseline)
        match_data{arr_idx}.together = [];
        match_data{arr_idx}.individual = [];
        match_data{arr_idx}.independence = [];
        match_data{arr_idx}.wave = [];
        match_data{arr_idx}.pos = [];
        for i = 1:size(num_spikes_above_baseline{arr_idx},1)
            for j = 1:size(num_spikes_above_baseline{arr_idx},2)
                % if not a single electrode, make independence prediction
                % and store plot data
                if(numel(arrayData{arr_idx}.CHAN_LIST{i}) > 1)
                    elec_idx = [];
                    for chan_idx = 1:numel(arrayData{arr_idx}.CHAN_LIST{i})
                        for list_idx = 1:numel(arrayData{arr_idx}.CHAN_LIST)
                            if(numel(arrayData{arr_idx}.CHAN_LIST{list_idx}) == 1 && arrayData{arr_idx}.CHAN_LIST{i}(chan_idx) == arrayData{arr_idx}.CHAN_LIST{list_idx})
                                elec_idx(end+1) = list_idx;
                            end
                        end                        
                    end
                    
                    match_data{arr_idx}.wave(end+1,1) = j;
                    match_data{arr_idx}.together(end+1,1) = num_spikes_above_baseline{arr_idx}(i,j);
                    match_data{arr_idx}.individual(end+1,:) = num_spikes_above_baseline{arr_idx}(elec_idx,j)';
                    match_data{arr_idx}.independence(end+1,1) = sum(num_spikes_above_baseline{arr_idx}(elec_idx,j)) - prod(num_spikes_above_baseline{arr_idx}(elec_idx,j));
                    idx_elec1 = find(map_data.chan == elec_idx(1));
                    idx_elec2 = find(map_data.chan == elec_idx(2));
                    match_data{arr_idx}.pos(end+1,:) = [map_data.row(idx_elec1),map_data.col(idx_elec1),...
                        map_data.row(idx_elec2),map_data.col(idx_elec2)]; % row1, col1, row2, col2
                    
                end % end if
            end
        end
        
    end

%% plot observed vs independence
    f=figure(); % for all on one plot

    for arr_idx = 1:numel(match_data)
%         f = figure(); % for individual plots
        for j = 1:size(num_spikes_above_baseline{arr_idx},2)
            plot(match_data{arr_idx}.independence(match_data{arr_idx}.wave == j),match_data{arr_idx}.together(match_data{arr_idx}.wave == j),'.','markersize',16,'color',getColorFromList(1,j))
            ax = gca;
            hold on
        end
        plot([0,1],[0,1],'r--','linewidth',1.5)
        xlim([0,1]);
        ylim([0,1]);
        xlabel('Independence Prediction');
        ylabel('Observed');
        formatForLee(gcf);
        set(gca,'fontsize',14);
    end
    
%% plot gain (together vs best individual)
    f=figure(); % for all on one plot

    for arr_idx = 1:numel(match_data)
%         f = figure(); % for individual plots
        for j = 1:size(num_spikes_above_baseline{arr_idx},2)
            plot(max(match_data{arr_idx}.individual(match_data{arr_idx}.wave == j,:),[],2),match_data{arr_idx}.together(match_data{arr_idx}.wave == j),'.','markersize',16,'color',getColorFromList(1,j))
            hold on
        end
        ax = gca;

        plot([0,1],[0,1],'r--','linewidth',1.5)
        xlim([0,1]);
        ylim([0,1]);
        xlabel('Best individual');
        ylabel('Two channels');
        formatForLee(gcf);
        set(gca,'fontsize',14);
    end
    
%% plot gain with distance from neuron each electrode on each axis (gain = color)
    f=figure(); % for all on one plot
    colors = inferno;
    max_gain = 0.1;
    min_gain = 0;
    for arr_idx = 1:numel(match_data)
%         f = figure(); % for individual plots
        for idx = 1:numel(match_data{arr_idx}.together)
            gain = match_data{arr_idx}.together(idx) - max(match_data{arr_idx}.individual(idx,:),[],2);
            gain = min(max(gain,min_gain+eps),max_gain); % force between bounds
            
            color_to_plot = colors(ceil((gain-min_gain)/(max_gain-min_gain)*size(colors,1)),:);
%             color_to_plot = colors(ceil(match_data{arr_idx}.together(idx)*size(colors,1)),:);

            dist_1 = 0.4*sqrt((match_data{arr_idx}.pos(idx,1)-arrayData{arr_idx}.ROW).^2 + ...
                (match_data{arr_idx}.pos(idx,2)-arrayData{arr_idx}.COL).^2);
            
            dist_2 = 0.4*sqrt((match_data{arr_idx}.pos(idx,2)-arrayData{arr_idx}.ROW).^2 + ...
                (match_data{arr_idx}.pos(idx,4)-arrayData{arr_idx}.COL).^2);
            
            plot(dist_1,dist_2,'.','markersize',16,'color',color_to_plot)
            hold on
        end
        ax = gca;

        xlabel('Distance to 1st electrode (mm)');
        ylabel('Distance to 2nd electrode (mm)');
        formatForLee(gcf);
        set(gca,'fontsize',14);
        
        b = colorbar;
        colormap(inferno);
    end

%% plot gain vs distance between two electrodes
    f=figure(); % for all on one plot

    for arr_idx = 1:numel(match_data)
%         f = figure(); % for individual plots
        for j = 1:size(num_spikes_above_baseline{arr_idx},2)
            gains = match_data{arr_idx}.together(match_data{arr_idx}.wave==j) - max(match_data{arr_idx}.individual(match_data{arr_idx}.wave==j,:),[],2);
            distances = sqrt((match_data{arr_idx}.pos(match_data{arr_idx}.wave==j,1)-match_data{arr_idx}.pos(match_data{arr_idx}.wave==j,3)).^2 + ...
                (match_data{arr_idx}.pos(match_data{arr_idx}.wave==j,2)-match_data{arr_idx}.pos(match_data{arr_idx}.wave==j,4)).^2);
            
            plot(distances,gains,'.','markersize',16,'color',getColorFromList(1,j))
            hold on
        end
        ax = gca;
        xlabel('Distance between electrodes');
        ylabel('Gain');
        formatForLee(gcf);
        set(gca,'fontsize',14);
    end


    
%% plot gain as a color with P(elec1) and P(elec2) on the axes
    f=figure(); % for all on one plot
    colors = inferno;
    max_gain = 0.1;
    min_gain = 0;
    gains = [];
    for arr_idx = 1:numel(match_data)
%         f = figure(); % for individual plots
        for idx = 1:numel(match_data{arr_idx}.together)
            gain = match_data{arr_idx}.together(idx) - max(match_data{arr_idx}.individual(idx,:),[],2);
            gains(end+1) = gain;
            gain = min(max(gain,min_gain+eps),max_gain);
            color_to_plot = colors(ceil((gain-min_gain)/(max_gain-min_gain)*size(colors,1)),:);
            plot(max(match_data{arr_idx}.individual(idx,:),[],2),min(match_data{arr_idx}.individual(idx,:),[],2),'.','markersize',16,'color',color_to_plot)
            hold on
        end
        ax = gca;

        plot([0,1],[0,1],'r--','linewidth',1.5)
        xlim([0,1]);
        ylim([0,1]);
        xlabel('Best individual');
        ylabel('Second best individual');
        formatForLee(gcf);
        set(gca,'fontsize',14);
        
        b = colorbar;
        colormap(inferno);
    end