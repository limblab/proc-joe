%% build match data with bootstrapping for all neurons in array data
% and repeat
    params.monkey_name = 'Duncan';
    params.bootstrap = 0;
    params.num_boot = 1;
    params.array_name = 'LeftS1';
    params.post_stim_window_size = 3;
    
% remove indices from arrayData that have bad artifacts
    arrayDataOld = arrayData;
    if(exist('artifact_keep_mask'))
        keep_idx = ones(size(arrayData));
        for arr_idx = 1:numel(arrayData)
            keep_idx(arr_idx) = artifact_keep_mask(arrayData{arr_idx}.CHAN_REC);
        end
        arrayData = arrayData(keep_idx==1);
    else
        warning('using all channels without looking at artifact settling');
    end


%% get match data and bootstrap (if params.bootstrap)
    match_data = {};
    for arr_idx = 1:numel(arrayData)
        % compute metrics (prob spike, independent prob)
        match_data{arr_idx} = buildMatchData(arrayData{arr_idx},params);
    end
    disp('done');
    
    

%% plot observed vs independence, with confidence bounds (95% based on bootstrap)
    min_prob = 0;
    f=figure(); % for all on one plot
    f.Name = [params.monkey_name,'_',params.array_name,'_ObservedVsIndependence'];
    for arr_idx = 1:numel(match_data) % for each neuron
%         disp(match_data{arr_idx}.chan_rec);
%         f = figure(); % for individual plots
        for cond = 1:size(match_data{arr_idx}.independence,1)
            % get conf bounds (95% bootstrap)
            if(sum(match_data{arr_idx}.individual(cond,:) > min_prob) >= size(match_data{arr_idx}.individual,2))
                independence_bounds = [prctile(squeeze(match_data{arr_idx}.independence(cond,:,:)),0.025),prctile(squeeze(match_data{arr_idx}.independence(cond,:,:)),0.975)];
                together_bounds = [prctile(squeeze(match_data{arr_idx}.together(cond,:,:)),0.025),prctile(squeeze(match_data{arr_idx}.together(cond,:,:)),0.975)];

                % plot with errorbars, color each amplitude differently
                errorbar(squeeze(match_data{arr_idx}.independence(cond,:,1)),squeeze(match_data{arr_idx}.together(cond,:,1)),...
                    abs(squeeze(match_data{arr_idx}.together(cond,:,1) - together_bounds(1))),abs(squeeze(match_data{arr_idx}.together(cond,:,1) - together_bounds(2))),...
                    abs(squeeze(match_data{arr_idx}.independence(cond,:,1) - independence_bounds(1))),abs(squeeze(match_data{arr_idx}.independence(cond,:,1) - independence_bounds(2))),...
                    '.','markersize',16,'color',getColorFromList(1,match_data{arr_idx}.wave(cond)-1))

                hold on
%             if(match_data{arr_idx}.independence(cond) > 0 && match_data{arr_idx}.together(cond) < 0 && match_data{arr_idx}.wave(cond) == 1)
%                 disp(arr_idx)
%             end
            end
        end

    end
    
    plot([0,1],[0,1],'k--','linewidth',1.5)
    xlim([0,1]);
    ylim([0,1]);
    xlabel('Independence Prediction');
    ylabel('Observed');
    formatForLee(gcf);
    set(gca,'fontsize',14);
    
    
    
%% histogram of differences between observed and independence

    f=figure(); % for all on one plot
    f.Name = [params.monkey_name,'_',params.array_name,'_ObservedVsIndependence'];
    difference_data = [];
    for arr_idx = 1:numel(match_data) % for each neuron
        % get difference between independence and together, also store wave
        independence_mask = squeeze(match_data{arr_idx}.independence(:,:,1)) > 0.; % remove times where neurons don't respond 
        difference_data = [difference_data; squeeze(match_data{arr_idx}.together(independence_mask,:,1)) - squeeze(match_data{arr_idx}.independence(independence_mask,:,1)), squeeze(match_data{arr_idx}.wave(independence_mask))];

        
    end
    
    histogram(difference_data(difference_data(:,2) == 3),[-1:0.025:1])
    formatForLee(gcf);
    set(gca,'fontsize',14);
 
    
    
%% plot observed - independence as a function of the worst elec prob of response
    f=figure(); % for all on one plot
    f.Name = [params.monkey_name,'_',params.array_name,'_dependenceVsWorstResponse'];

    normalize_y_data = 0;
    
    x_data_all = cell(numel(unique(match_data{arr_idx}.wave)),1);
    y_data_all = cell(numel(unique(match_data{arr_idx}.wave)),1);
    for arr_idx = 1:numel(match_data)
%         f = figure(); % for individual plots
        for j = 1:numel(unique(match_data{arr_idx}.wave))
            plot_mask = match_data{arr_idx}.together(:,:,1) > 0.025;
            x_data = (min(match_data{arr_idx}.individual(match_data{arr_idx}.wave == j & plot_mask,:,1),[],2));
            y_data = match_data{arr_idx}.together(match_data{arr_idx}.wave == j & plot_mask,:,1) - ...
                match_data{arr_idx}.independence(match_data{arr_idx}.wave == j & plot_mask,:,1);
            
            
            plot(x_data,y_data,...
                '.','markersize',16,'color',getColorFromList(1,j-1))
            
            x_data_all{j} = [x_data_all{j}; x_data];
            y_data_all{j} = [y_data_all{j}; y_data];
            
            ax = gca;
            hold on
        end
        xlabel('Worst electrode strength');
        ylabel('Observed - Independence');
        formatForLee(gcf);
        set(gca,'fontsize',14);
    end

%     FITS = {};
%     gof = {};
%     for j = 1:numel(x_data_all)
%         [FITS{j},gof{j}] = fit(x_data_all{j}(x_data_all{j} < 0.8),y_data_all{j}(x_data_all{j} < 0.8),'a*x+b');
%         plot([0,0.6],feval(FITS{j},[0,0.6]),'--','color',getColorFromList(1,j-1),'linewidth',2)
%     end
    
%% plot observed-independence as a function of the distance between electrodes and the neuron
    f=figure(); % for all on one plot
    f.Name = [params.monkey_name,'_',params.array_name,'_dependenceVsDistanceElectrodes'];
        
    distance_data_all = cell(numel(unique(match_data{arr_idx}.wave)),1);
    difference_data_all = cell(numel(unique(match_data{arr_idx}.wave)),1);
    for arr_idx = 1:numel(match_data)
        neuron_pos = match_data{arr_idx}.chan_rec_pos;
        
        for j = 1:size(match_data{arr_idx}.wave,1)
            stim_chan_pos_list = match_data{arr_idx}.pos{j};
            
            
            distance_data = 400*sqrt((stim_chan_pos_list(:,1)-neuron_pos(:,1)).^2 + (stim_chan_pos_list(:,2)-neuron_pos(:,2)).^2); % distances
            distance_data = mean(distance_data);
            
            difference_data = squeeze(match_data{arr_idx}.together(j,:,1)) - squeeze(match_data{arr_idx}.independence(j,:,1));
            
            distance_data_all{match_data{arr_idx}.wave(j)} = [distance_data_all{match_data{arr_idx}.wave(j)}; distance_data];
            difference_data_all{match_data{arr_idx}.wave(j)} = [difference_data_all{match_data{arr_idx}.wave(j)}; difference_data];
            
            plot(distance_data,difference_data,'.','markersize',16,'color',getColorFromList(1,match_data{arr_idx}.wave(j)-1))
            
            ax = gca;
            hold on
        end
    end
    xlabel('Mean distance between neuron and electrodes (\mum)');
    ylabel('Observed - Independence');
    formatForLee(gcf);
    set(gca,'fontsize',14);

    
%% FOR TWO ELECTRODES : plot observed-independence as a function of the distance between each electrode and the unit
    wave_to_plot = 1;

    f=figure(); % for all on one plot
    f.Name = [params.monkey_name,'_',params.array_name,'_dependence_distanceNeuronAndElectrode_100uA'];
    
    colors = inferno;
    max_gain = 0.3;
    min_gain = 0;
    gain_all = [];
    dist_1_all = [];
    dist_2_all = [];
    for arr_idx = 1:numel(match_data)
%         f = figure(); % for individual plots
        neuron_pos = match_data{arr_idx}.chan_rec_pos;
        
        for idx = 1:size(match_data{arr_idx}.together,1)
            if(min(match_data{arr_idx}.individual(idx,:,1)) > -10.00 && match_data{arr_idx}.wave(idx) == wave_to_plot)
                gain = match_data{arr_idx}.independence(idx,:,1) - match_data{arr_idx}.together(idx,:,1);
                gain_all(end+1) = gain;
                gain = min(max(gain,min_gain+eps),max_gain); % force between bounds

                color_to_plot = colors(ceil((gain-min_gain)/(max_gain-min_gain)*size(colors,1)),:);
                
                stim_chan_pos_list = match_data{arr_idx}.pos{idx};

                dist_1 = 400*sqrt((stim_chan_pos_list(1,1)-neuron_pos(1,1)).^2 + ...
                        (stim_chan_pos_list(1,2)-neuron_pos(1,2)).^2);
                dist_2 = 400*sqrt((stim_chan_pos_list(2,1)-neuron_pos(1,1)).^2 + ...
                        (stim_chan_pos_list(2,2)-neuron_pos(1,2)).^2);

    
                dist_1_all(end+1) = dist_1;
                dist_2_all(end+1) = dist_2;
                
                dist_1 = dist_1 + rand(size(dist_1))*80 - 40;
                dist_2 = dist_2 + rand(size(dist_2))*80 - 40;
                
                plot(dist_1,dist_2,'.','markersize',16,'color',color_to_plot)
                hold on
            end
        end
        
    end
    
    ax = gca;
    xlabel('Distance to elec 1 (\mum)');
    ylabel('Distance to elec 2 (\mum)');
    formatForLee(gcf);
    set(gca,'fontsize',14);
    xlim([0,5500])
    ylim([0,5500])
    b = colorbar;
    b.Label.String = 'Observed-Independence';
    b.TickDirection = 'out';
    colormap(flip(inferno));
    b.Label.FontSize = 14;
    tick_labels = flip(linspace(min_gain,max_gain,numel(b.TickLabels)));
    for t = 1:numel(b.TickLabels)
        if(t == 1)
            b.TickLabels{t} = ['< -',num2str(max_gain)];
        elseif(t == numel(b.TickLabels))
            b.TickLabels{t} = ['> ',num2str(min_gain)];                
        else
            b.TickLabels{t} = -tick_labels(t);
        end
    end
    ax.Position(3:4) = ax.Position(3:4)*0.975;

    
    tbl = table(dist_1_all',dist_2_all',gain_all','VariableNames',{'dist_1','dist_2','gain'});
    lm = fitlm(tbl,'gain~dist_1*dist_2')
    

%% FOR TWO ELECTRODES : plot amplitude of one electrode vs amplitude of second electrode
    prob_list = [];
    wave_list = [];
    x_edges = [-0.1:0.05:1];
    y_edges = x_edges;
    for i_unit = 1:numel(match_data)
        prob_list = [prob_list; match_data{i_unit}.individual];
        wave_list = [wave_list; match_data{i_unit}.wave];
    end
    
    [n_his, x_edges, y_edges] = histcounts2(prob_list(:,1),prob_list(:,2),x_edges,y_edges);
    n_his = log10(n_his);
    x_center = x_edges(1:end-1) + mode(diff(x_edges))/2;
    y_center = y_edges(1:end-1) + mode(diff(y_edges))/2;
    
    imagesc(x_center, y_center, n_his)
    colormap(inferno);
    b=colorbar;
    set(gca,'YDir','normal');
    formatForLee(gcf);
    xlabel('Response to elec 1');
    ylabel('Response to elec 2');
    set(gca,'fontsize',14);
    
    b.Label.String = 'log_1_0(Count)';
    
%% plot observed vs. independence as a sum of total neural activity
    f=figure(); % for all on one plot
    f.Name = [params.monkey_name,'_',params.array_name,'_popResponse'];

    chan_list = unique(match_data{1}.chans,'rows');
    num_pairs = size(chan_list,1);
    mean_total_spikes = zeros(num_pairs,numel(unique(match_data{1}.wave)));
    mean_independence_spikes = zeros(num_pairs,numel(unique(match_data{1}.wave)));
    match_data_idx = 1;
    
    for pair_idx = 1:num_pairs
        for j = 1:numel(unique(match_data{1}.wave))
            for arr_idx = 1:numel(match_data)
                if(match_data{arr_idx}.independence(match_data_idx) > 0)
                    mean_total_spikes(pair_idx,j) = mean_total_spikes(pair_idx,j) + match_data{arr_idx}.together(match_data_idx);
                    mean_independence_spikes(pair_idx,j) = mean_independence_spikes(pair_idx,j) + match_data{arr_idx}.independence(match_data_idx);
                end
            end
            plot(mean_independence_spikes(pair_idx,j),mean_total_spikes(pair_idx,j),'.','markersize',16,'color',getColorFromList(1,j-1))
%             plot(j,mean_total_spikes(pair_idx,j)./mean_independence_spikes(pair_idx,j),'.','markersize',16,'color',getColorFromList(1,j-1))
            hold on
            match_data_idx = match_data_idx + 1;
        end
    end
    
    plot([0,15],[0,15],'r--','linewidth',1.5)
    
    xlabel('Pop predicted (spikes/stim)');
    ylabel('Pop observed (spikes/stim)');
    formatForLee(gcf);
    set(gca,'fontsize',14);
    
%% plot observed vs. independence as a ratio of total neural activity 
    f=figure(); % for all on one plot
    f.Name = [params.monkey_name,'_',params.array_name,'_popResponse'];

    chan_list = unique(match_data{1}.chans,'rows');
    num_pairs = size(chan_list,1);
    mean_total_spikes = zeros(num_pairs,numel(unique(match_data{1}.wave)));
    mean_independence_spikes = zeros(num_pairs,numel(unique(match_data{1}.wave)));
    match_data_idx = 1;
    
    for pair_idx = 1:num_pairs
        for j = 1:numel(unique(match_data{1}.wave))
            for arr_idx = 1:numel(match_data)
                if(match_data{arr_idx}.independence(match_data_idx) > 0)
                    mean_total_spikes(pair_idx,j) = mean_total_spikes(pair_idx,j) + match_data{arr_idx}.together(match_data_idx);
                    mean_independence_spikes(pair_idx,j) = mean_independence_spikes(pair_idx,j) + match_data{arr_idx}.independence(match_data_idx);
                end
            end
            plot(j,mean_total_spikes(pair_idx,j)./mean_independence_spikes(pair_idx,j),'.','markersize',16,'color',getColorFromList(1,j-1))
            hold on
            match_data_idx = match_data_idx + 1;
        end
    end
    
    plot([0,15],[0,15],'r--','linewidth',1.5)
    
    xlabel('Pop predicted (spikes/stim)');
    ylabel('Pop observed (spikes/stim)');
    formatForLee(gcf);
    set(gca,'fontsize',14);

%% separate arrayData into stimulations on single electrodes and stimulation on many electrodes
arrayData_single = {};
arrayData_multi = {};
for arr_idx = 1:numel(arrayData)
    % get idx of single channel
    single_chan_mask = (cellfun(@numel,arrayData{1}.CHAN_LIST)==1);
    arrayData_single{arr_idx} = arrayData{arr_idx}; 
    arrayData_multi{arr_idx} = arrayData{arr_idx};
    
    % 
    
    % go through each field and separate if needed
    for f = fieldnames(arrayData{arr_idx})'
       if(size(arrayData{arr_idx}.(f{1}),1) == numel(single_chan_mask)) 
            arrayData_single{arr_idx}.(f{1}) = arrayData{arr_idx}.(f{1})(single_chan_mask,:);
            arrayData_multi{arr_idx}.(f{1}) = arrayData{arr_idx}.(f{1})(~single_chan_mask,:);
       elseif(strcmpi(f{1},'CHAN_SENT') || strcmpi(f{1},'WAVEFORM_SENT')) % deal with CHAN_SENT and WAVEFORM_SENT
           chan_sent_mask = cellfun(@numel,arrayData{arr_idx}.CHAN_SENT)==1;
           arrayData_single{arr_idx}.(f{1}) = arrayData{arr_idx}.(f{1})(chan_sent_mask);
           arrayData_multi{arr_idx}.(f{1}) = arrayData{arr_idx}.(f{1})(~chan_sent_mask);
           
       end
    
        
    end
    
    
    arrayData_single{arr_idx}.CHAN_SENT = [arrayData_single{arr_idx}.CHAN_SENT{:}]';
    arrayData_single{arr_idx}.CHAN_LIST = [arrayData_single{arr_idx}.CHAN_LIST{:}]';
end
    
    



%% plot PSTH's with observed and predicted
    arr_idx = 1; % 4, 26
    wave = 4;
    chans_to_plot = [21,42];
    window = [-5,10]; % ms
    
    optsPlot.NUM_PLOTS = [];
    optsPlot.BAR_STYLE = 'line';
    optsPlot.LINE_WIDTH = 2;
    optsSave = [];
    
    x_data = [];
    y_data = [];
    % get individual responses
    for chan = chans_to_plot
        chan_list_idx = find(checkChanListEquality(arrayData{arr_idx}.CHAN_LIST,chan));
        window_idx = [find(arrayData{arr_idx}.bE{1,1} > window(1),1,'first'), find(arrayData{arr_idx}.bE{1,1} > window(2),1,'first')];
        x_data(end+1,:) = arrayData{arr_idx}.bE{chan_list_idx,wave}(window_idx(1):window_idx(2)) + mode(diff(arrayData{arr_idx}.bE{chan_list_idx,wave}))/2;
        y_data(end+1,:) = arrayData{arr_idx}.bC{chan_list_idx,wave}(window_idx(1):window_idx(2));
    end

    % get multichannel stim response
    chan_list_idx = find(checkChanListEquality(arrayData{arr_idx}.CHAN_LIST,chans_to_plot));
    window_idx = [find(arrayData{arr_idx}.bE{1,1} > window(1),1,'first'), find(arrayData{arr_idx}.bE{1,1} > window(2),1,'first')];
    x_data(end+1,:) = arrayData{arr_idx}.bE{chan_list_idx,wave}(window_idx(1):window_idx(2)) + mode(diff(arrayData{arr_idx}.bE{chan_list_idx,wave}))/2;
    y_data(end+1,:) = arrayData{arr_idx}.bC{chan_list_idx,wave}(window_idx(1):window_idx(2));
    
    
    % get summation of individual respones (independence)
    x_data(end+1,:) = x_data(end,:);
    y_data(end+1,:) = y_data(1,:) + y_data(2,:) - y_data(1,:).*y_data(2,:);
    
    % plot
    plotPSTHLIB(x_data',y_data',optsPlot,optsSave);
    l=legend('elec1','elec2','elec1 and elec2','independence');
    set(l,'box','off','location','best');
    ylabel('Num spikes per stimulation');
    xlabel('Time after stimulation onset (ms)');
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    