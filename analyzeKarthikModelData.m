%% cocatenate data across stim iterations and make a struct for each neuron type
% for each neuron, we want spike times relative to stim onset within a
% window around each stimulation
    amps = [5,10,15,20,25,30,40,50,100];
    
    wave_length = 0.453;
    ICMS_window = [-250,250];
    ICMS_freq = 2;
    ICMS_times = 2000:1000/ICMS_freq:(10000-1); % should be 16
    
    neuron_struct_all = {};
    iteration_list = []; neuron_num_list = [];
    dist_list = [];
    data_folder = 'D:\Lab\Data\StimArtifact\singlePulseModelExperiment\data\'; % must have file_sep after
    neuron_types = {'L23','L4','L5','L6'};
    neuron_type_list = [];
    % load soma coordinate file
    for i_type = 1:numel(neuron_types)
        soma_file = dir([data_folder,neuron_types{i_type},'_soma*']);
        load([data_folder,soma_file(1).name]);

        % go through each iteration and amplitude and collect data
        for i_amp = 1:numel(amps)
            iteration_file = dir([data_folder,neuron_types{i_type},'_amp_',num2str(amps(i_amp)),'uA*']);
            for i_iteration = 1:numel(iteration_file)
                load([data_folder,iteration_file(i_iteration).name]);

                % for each data entry (neuron)
                for i_neuron = 1:numel(data)
                    % see if neuron exists already
                    neuron_idx = find(neuron_type_list == i_type & iteration_list == i_iteration & neuron_num_list == i_neuron);
                    if(isempty(neuron_idx))
                        % add entry
                        neuron_idx = numel(neuron_struct_all) + 1;
                        iteration_list(end+1) = i_iteration; neuron_num_list(end+1) = i_neuron; neuron_type_list(end+1) = i_type;

                        neuron_struct_all{neuron_idx}.coordinates = iteration(i_iteration).soma_coord(i_neuron,:);
                        neuron_struct_all{neuron_idx}.dist = norm(neuron_struct_all{neuron_idx}.coordinates,2); 
                        dist_list(end+1) = neuron_struct_all{neuron_idx}.dist; 
                        neuron_struct_all{neuron_idx}.type = neuron_types{i_type};
                        
                        
                        neuron_struct_all{neuron_idx}.spike_times = {}; % spike times relative to stimulation
                        neuron_struct_all{neuron_idx}.stim_idx = {}; % which stimulation each spike time happened for
                        neuron_struct_all{neuron_idx}.amps = []; % amplitude list
                        neuron_struct_all{neuron_idx}.baseline_fr = []; % baseline fr per amp
                    end

                    spike_times = []; stim_idx = [];
                    % get spike times relative to stim onset
                    for i_stim = 1:numel(ICMS_times)
                        spike_mask = data(i_neuron).times > ICMS_times(i_stim) + ICMS_window(1) & ...
                            data(i_neuron).times < ICMS_times(i_stim) + ICMS_window(2);

                        spike_times = [spike_times, data(i_neuron).times(spike_mask) - ICMS_times(i_stim)];
                        stim_idx = [stim_idx, i_stim*ones(1,sum(spike_mask))];
                    end

                    % update entry with spike times for new amplitude
                    neuron_struct_all{neuron_idx}.amps(end+1,1) = amps(i_amp);
                    neuron_struct_all{neuron_idx}.spike_times{end+1,1} = spike_times;
                    neuron_struct_all{neuron_idx}.stim_idx{end+1,1} = stim_idx;
                    neuron_struct_all{neuron_idx}.baseline_fr(end+1,1) = sum(spike_times < 0)/(abs(ICMS_window(1)/1000)*numel(ICMS_times));
                end

            end
        end
    end

    firing_rate_list = [];
    for i_neuron = 1:numel(neuron_struct_all)
        firing_rate_list(i_neuron) = mean(neuron_struct_all{i_neuron}.baseline_fr);
    end
    

%% plot activation vs distance for the model

    
    bin_edges_model = [0,3]; 
    dist_range = [0,200];
    wave_length = 0;
    bin_counts_model = getProbSpikeKarthikModel(neuron_struct_all,bin_edges_model,dist_range,numel(ICMS_times),amps,wave_length,0); % don't do mean

    neuron_mask = squeeze(bin_counts_model(9,1,:)) > -1000;
    neuron_struct = neuron_struct_all(neuron_mask == 1);
    
    bin_counts_model = getProbSpikeKarthikModel(neuron_struct,bin_edges_model,dist_range,numel(ICMS_times),amps,wave_length,0); % don't do mean

    figure()
    for amp_idx = 1:9
        subplot(3,3,amp_idx)
%         densityplot(dist_list(neuron_mask==1),squeeze(bin_counts_model(amp_idx,1,:)),'edges',{[dist_range(1):20:dist_range(end)],[0:0.05:1]});
        plot(dist_list(neuron_mask==1),squeeze(bin_counts_model(amp_idx,1,:)),'.')
        ylim([-0.2,1.2]);
        formatForLee(gcf)
        if(amp_idx == 1)
            xlabel('Distance (\mum)')
            ylabel('Response amp')
        end
    end
%% plot mean FR at different bin points for model neurons at specific distances
% compare to recorded neurons (arrayData)
    bin_size = 1; min_bin = 0; max_bin = 7;
    dist_range = [0,200];
    adj_wave_time = 0;
    
    bin_edges_model = [min_bin:bin_size:max_bin];
    
    bin_counts_model = zeros(numel(amps),numel(bin_edges_model)-1);
    num_neurons_model = zeros(numel(amps),1);
    
    bin_counts_exp = zeros(numel(amps),numel(bin_edges_model)-1);
    num_neurons_exp = zeros(numel(amps),1);
    
    % bin model data
    if(adj_wave_time)
        wave_length = 0.453;
    else
        wave_length = 0;
    end
    do_mean = 1;
    bin_counts_model = getProbSpikeKarthikModel(neuron_struct,bin_edges_model,dist_range,numel(ICMS_times),amps,wave_length,do_mean);
    
    for i_amp = 1:numel(amps)
        % bin experimental data
        for i_neuron = 1:numel(arrayData)
            amp_idx = find([arrayData{i_neuron}.STIM_PARAMETERS.amp1] == amps(i_amp),1,'first');
            amp_100uA_idx = find(spikesStruct{i_neuron}.amp == 100,1,'first');
            if(~isempty(amp_idx))
                baseline_count = spikesStruct{i_neuron}.mean_baseline_fr*bin_size/1000;
                bin_counts_exp(i_amp,:) = bin_counts_exp(i_amp,:) +  ... 
                    histcounts(arrayData{i_neuron}.spikeTrialTimes{amp_idx} - 0.453/1000,bin_edges_model/1000)/arrayData{i_neuron}.numStims(amp_idx) - baseline_count;
                num_neurons_exp(i_amp) = num_neurons_exp(i_amp) + 1;
            end
        end
        bin_counts_exp(i_amp,:) = bin_counts_exp(i_amp,:)./num_neurons_exp(i_amp);
    end
    
    
    % plot bin counts for stim and model on same figures
    y_limits = [-0.1,0.4]; x_limits = [0,7];
    amp_idx = [1,2,4,6,7,8];
    subplot_idx = [1,1,2,2,3,3,];
    color_idx = [0,1,2,3,4,5,6,7,8];
    x_data = bin_edges_model(1:end-1) + mode(diff(bin_edges_model))/2;
    figure('Position',[2014 462 1713 420]);
    for i = 1:numel(amp_idx)
        subplot(1,max(subplot_idx),subplot_idx(i)); hold on;
        % plot model (dashed line, squares)
        plot(x_data,bin_counts_model(amp_idx(i),:),'s--','color',getColorFromList(1,color_idx(i)),'markersize',6)
        % plot experimental (solid line, dot)
        plot(x_data,bin_counts_exp(amp_idx(i),:),'.-','color',getColorFromList(1,color_idx(i)),'markersize',18)
        
        % plot no recording zone
        
        ylim(y_limits);
        xlim(x_limits);
        
        formatForLee(gca)
        set(gca,'fontsize',14)
        xlabel('Time after stim offset (ms)')
        ylabel('Spikes per stim per neuron')
    end
    
%% plot activation across neurons for model and experiment

    bin_edges_model = [0,1.5;1.5,3]; 
    dist_range = [0,200];
    wave_length = 0;
    model_offset = [-1.5,-0.75,0];
    model_colors = [getColorFromList(1,1); getColorFromList(1,2); getColorFromList(1,3)];
    bin_edges_exp = [1.5,3];
    exp_offset = [0.75,1.5]; exp_colors = [getColorFromList(1,0); getColorFromList(1,4)];
        
    figure(); subplot(1,2,1); hold on;
    % make experiment plots
    for i_bin = 1:size(bin_edges_exp,1)
        bin_counts_exp = nan(numel(amps),numel(arrayData));
        num_neurons_exp = zeros(numel(amps),1);
        for i_amp = 1:numel(amps)
            % bin experimental data
            for i_neuron = 1:16%numel(arrayData)
                amp_idx = find([arrayData{i_neuron}.STIM_PARAMETERS.amp1] == amps(i_amp),1,'first');
                amp_100uA_idx = find(spikesStruct{i_neuron}.amp == 100,1,'first');
                if(~isempty(amp_idx))
                    baseline_count = spikesStruct{i_neuron}.mean_baseline_fr*bin_size/1000;
                    bin_counts_exp(i_amp,i_neuron) = histcounts(arrayData{i_neuron}.spikeTrialTimes{amp_idx} - 0.453/1000,bin_edges_exp(i_bin,:)/1000)/arrayData{i_neuron}.numStims(amp_idx) - baseline_count;
                    num_neurons_exp(i_amp) = num_neurons_exp(i_amp) + 1;
                end
            end
            bin_counts_exp(i_amp,:) = bin_counts_exp(i_amp,:);
        end
        
        errorbar(amps+exp_offset(i_bin),mean(bin_counts_exp,2,'omitnan'),std(bin_counts_exp,0,2,'omitnan'),'linestyle','none',...
            'marker','.','markersize',20,'linewidth',1.5,'color',exp_colors(i_bin,:)); 
        
    end

    
    % make model plots

    subplot(1,2,2); hold on
    for i_bin = 1:size(bin_edges_model,1)
        bin_counts_model = squeeze(getProbSpikeKarthikModel(neuron_struct,bin_edges_model(i_bin,:),dist_range,numel(ICMS_times),amps,wave_length,0)); % don't do mean

        errorbar(amps+model_offset(i_bin),mean(bin_counts_model,2,'omitnan'),std(bin_counts_model,0,2,'omitnan'),'linestyle','none',...
            'marker','s','markersize',8,'linewidth',1.5,'color',model_colors(i_bin,:)); 
    end
    
% %% plot raster for a given neuron
% 
%     neuron_idx = 12;
%     num_stims = numel(ICMS_times);
%     
%     figure(); hold on
%     stim_idx = 1;
%     for i_amp = 1:numel(neuron_struct{neuron_idx}.amps)
%         for i_stim = 1:num_stims
%             spike_mask = neuron_struct{neuron_idx}.stim_idx{i_amp} == i_stim;
%             
%             plot(neuron_struct{neuron_idx}.spike_times{i_amp}(spike_mask),stim_idx*ones(1,sum(spike_mask)),'k.')
%             
%             stim_idx = stim_idx + 1;
%         end
%         
%         if(i_amp < numel(neuron_struct{neuron_idx}.amps))
%             plot(ICMS_window,stim_idx-0.5 + [0,0],'r--')
%         end
%         
%     end
%     
%     xlim(ICMS_window)
%     ylim([0.5,stim_idx-0.5])
%     



