% get model data -- space constant DONT RUN, LOAD .MAT FILE INSTAED
    mdl_input_data = [];
    mdl_input_data.folderpath = 'C:\Users\Joseph Sombeck\Box\Miller-Grill_S1-stim\ModelData\SpaceConstant\';
%     mdl_input_data.folderpath = 'C:\Users\Joseph Sombeck\Box\Miller-Grill_S1-stim\ModelData\StimChannelResponse\no_Intrinsic_Activity\';
    mdl_input_data.diam_list = [1,2,3];
    mdl_input_data.amp_list = [15,30,50,100];
    mdl_input_data.get_axon_dendrite_locs = 0;
    mdl_input_data.cell_id_list = [1,6,11,16,21];
    mdl_input_data.num_clones = 5;
    mdl_input_data.stim_times = (201); % ms
%     mdl_input_data.stim_times = (201)+(0:1:9)*500; % ms
    mdl_input_data.wave_length = 0.453; % ms
    mdl_input_data.stim_window = [-100,400]; % ms around stim
    mdl_input_data.get_IPIs = 0;
    mdl_input_data.gaba_ratio_list = [1];
    % cell_id=6:10 %L23 PC, clones 1-5
    % cell_id=11:15 %L4 LBC, clones 1-5
    % cell_id=16:20 %L5 PC, clones 1-5
    % cell_id=21:25 %L6 PC, clones 1-5
    
%     mdl_input_data.get_synapses = 1;
%     [mdl_data_all,mdl_syn_array_data_all,mdl_syn_mask_data_all] = getModelStimChannelData(mdl_input_data);    
%     mdl_syn_array_data = mdl_syn_array_data_all; mdl_syn_mask_data = mdl_syn_mask_data_all;
    mdl_input_data.temporal_data = 0;
    
    mdl_input_data.get_synapses = 0;
    [mdl_data_all,mdl_array_data_all,mdl_mask_data_all] = getModelStimChannelData(mdl_input_data);  
    mdl_array_data = mdl_array_data_all; mdl_mask_data = mdl_mask_data_all;
        
%% get experiment data -- space constant
    input_data.home_computer = 1;
    exp_array_data = getExperimentSpaceConstantData(input_data);
    exp_array_data = adjustArrayDataSpikeTimes(exp_array_data, 0.453/1000); % stim pulse length
    exp_array_data = getBaselineFiringRate(exp_array_data,[-25,-5]/1000); % window relative to stim onset
    

%% resample neurons based on proportions observed in rat cortex
%     num_sample = 200;
    [mdl_array_data,mdl_mask_data] = resampleModelData(mdl_array_data,mdl_mask_data);
%     [mdl_syn_array_data,mdl_syn_mask_data] = resampleModelData(mdl_syn_array_data_all,mdl_syn_mask_data_all);
    
%% get response data for both model and experiment
    resp_input_data = []; 
    resp_input_data.home_computer = 1;
    resp_input_data.monkey_list = {'Han','Duncan'};
    
    resp_input_data.sub_baseline = 1; % doesn't apply to model
    resp_input_data.max_lat_distance = 200; % max distance to count latency of a spike
    
%     resp_input_data.is_model = 0;
%     resp_input_data.spike_window = [1,5]/1000; % s
%     [exp_resp_data] = getResponseAndDistanceData(exp_array_data, resp_input_data);

    resp_input_data.is_model = 1;
    resp_input_data.spike_window = [-5,5]/1000; % s
    [mdl_resp_data] = getResponseAndDistanceData(mdl_array_data, resp_input_data);

%     resp_input_data.is_model = 1;
%     resp_input_data.spike_window = [0,5]/1000; % s
%     [mdl_syn_resp_data] = getResponseAndDistanceData(mdl_syn_array_data, resp_input_data);
    
    % adjust response amplitude because there's no intrinsic activity
%     mdl_syn_resp_data.response_amp = mdl_syn_resp_data.response_amp*0.5;
    
    % get experiment mean baseline FR
%     exp_baseline_fr = 0;
%     for i_unit = 1:numel(exp_array_data)
%         exp_baseline_fr = exp_baseline_fr + exp_array_data{i_unit}.baseline_fr;
%     end
%     exp_baseline_fr = exp_baseline_fr/numel(exp_array_data);
    
%% 
    diam = 2;
    amp_list = unique(mdl_resp_data.amp);
    color_list = inferno(numel(amp_list)+2);
    color_list(2,:) = [];
    
    bin_edges = [0:100:2000];
%     bin_edges = [0:10:150];
    [~,~,bin_idx] = histcounts(mdl_resp_data.dist_from_stim,bin_edges);
    resp_amp_bin = zeros(numel(amp_list),numel(bin_edges)-1);
    num_neurons_bin = zeros(numel(amp_list),numel(bin_edges)-1);
    
    f=figure();
    fits = {}; gof = {}; 
    bin_vol = 4/3*pi*(bin_edges(2:end).^3 - bin_edges(1:end-1).^3);
    for i_amp = 1:numel(amp_list)
        mask = mdl_resp_data.diam == diam & mdl_resp_data.amp == amp_list(i_amp);

        for i_bin = 1:size(resp_amp_bin,2)
            resp_amp_bin(i_amp,i_bin) = sum(mdl_resp_data.response_amp(mask & bin_idx == i_bin));
            num_neurons_bin(i_amp,i_bin) = sum(mask & bin_idx == i_bin);
        end
        
        density_bin = num_neurons_bin(i_amp,:)./bin_vol;
        ideal_density = mean(density_bin(bin_edges(1:end-1) > 200));
        correction_factor = density_bin./ideal_density;
        
        
        x=bin_edges(1:end-1)+mode(diff(bin_edges))/2;
        y=resp_amp_bin(i_amp,:);
%         fit_mask = ~isnan(y);
%         
% %         fits{i_amp} = fit(x(fit_mask)', y(fit_mask)','a*exp(-1*x/b)','lower',[0,0],'upper',[2,2000],'StartPoint',[1,500]);
%         [fits{i_amp},gof{i_amp}] = fit(x(fit_mask)', y(fit_mask)','1-1/(1+exp(-a*(x-b)))','lower',[-1,0],'upper',[1,2000],'StartPoint',[0.1,100]);
        [fits{i_amp},gof{i_amp}] = fit(mdl_resp_data.dist_from_stim(mask),mdl_resp_data.response_amp(mask),'1-1/(1+exp(-a*(x-b)))','lower',[-1,0],'upper',[1,2000],'StartPoint',[0.1,100]);
%        
        plot(x,resp_amp_bin(i_amp,:)./num_neurons_bin(i_amp,:),'color',color_list(i_amp,:),'marker','.','linestyle','none','markersize',30); hold on
        plot([0,x],feval(fits{i_amp},[0,x]),'color',color_list(i_amp,:),'linestyle','--','linewidth',2);
    
    end
    
    formatForLee(gcf);
    xlabel('Distance from stim (\mum)')
    ylabel('Proportion of neurons activated')
    set(gca,'fontsize',14);

%% using fits of distance and probability activated, estimate number of spikes evoked at different distances
bin_edges = 0:1:2000;
dBin = mode(diff(bin_edges));
bin_centers = bin_edges(1:end-1)+dBin/2;
density_factor = 1; % just scales result

bin_vol = (bin_edges(2:end) - bin_edges(1:end-1));
figure(); hold on;
for i_amp = 1:numel(fits)
    prob_act = feval(fits{i_amp},bin_centers);
    expected_num_spikes = prob_act.*bin_vol';
    plot(bin_centers,expected_num_spikes);
end

    
    
%% plot response vs distance -- use bins for experiment and model, model + synapses
% across all neurons for experiment
% across all clones for model
    max_dist = 1000*ceil(max([mdl_resp_data.dist_from_stim])/1000);
    bin_size = 100; % um
    bin_edges = 200:bin_size:max_dist;
    bin_centers = bin_edges(1:end-1) + mode(diff(bin_edges))/2;
%     [~,~,exp_bin_idx] = histcounts(exp_resp_data.dist_from_stim,bin_edges);
    [~,~,mdl_bin_idx] = histcounts(mdl_resp_data.dist_from_stim,bin_edges);
%     [~,~,mdl_syn_bin_idx] = histcounts(mdl_syn_resp_data.dist_from_stim,bin_edges);
    fit_bin_min = 0;
    
    mdl_diam = 2; % 2 is monkey diameter
    
    figure(); hold on;
    % experiment and model activations across all neurons/cell types/clones
    % (one diameter for model)
    fits = {}; gofs = {};
    ax_list = [];
    for i_type = 2:2
%         ax_list(end+1,1) = subplot(1,2,i_type); hold on;
        if(i_type == 1)
            resp_data = exp_resp_data;
            bin_idx = exp_bin_idx;
            type_mask = ones(size(exp_resp_data.response_amp));
            marker = '.';
            markersize = 20;
            linestyle = '-';
        elseif(i_type == 2) % model, no synapses
            resp_data = mdl_resp_data;
            bin_idx = mdl_bin_idx;
            type_mask = mdl_resp_data.diam == mdl_diam;
            marker = 's';
            markersize = 9;
            linestyle = '-';
        else % model synapses
            resp_data = mdl_syn_resp_data;
            bin_idx = mdl_syn_bin_idx;
            type_mask = mdl_syn_resp_data.diam == mdl_diam & mdl_syn_resp_data.cell_id == 21;
            marker = 'o';
            markersize = 9;
            linestyle = '-';
        end
        
        unique_amps = unique(resp_data.amp);

        for i_amp = 1:numel(unique_amps)
            ax_list(i_amp) = subplot(1,4,i_amp); hold on;
            mask = resp_data.amp == unique_amps(i_amp) & resp_data.dist_from_stim > 0 & type_mask;
            resp_amp = []; bin_val = [];
            for i_bin = 1:numel(bin_centers)
                if(sum(mask & bin_idx == i_bin) > 0)
                    resp_amp(end+1,1) = mean(resp_data.response_amp(mask & bin_idx == i_bin));
                    bin_val(end+1,1) = bin_centers(i_bin);
                end
            end

            plot(bin_val,resp_amp,'linestyle','none','marker',marker,'color',getColorFromList(1,i_type-1),'markersize',markersize);

            % get exponential fit
            [fits{i_type,i_amp},gofs{i_type,i_amp}] = fit(bin_val(bin_val >= fit_bin_min),resp_amp(bin_val >= fit_bin_min),'a*exp(-x/b)+c','Lower',[0,0,-0.4],'Upper',[5,10000,0.4],'StartPoint',[0.03,-0.005,0]);
            plot(bin_centers(1):50:bin_centers(end),feval(fits{i_type,i_amp},bin_centers(1):50:bin_centers(end)),linestyle,'color',getColorFromList(1,i_type-1),'linewidth',1.25)
        
            formatForLee(gcf);
            set(gca,'fontsize',14);
            xlabel('Dist from stim elec (\mum)');
            ylabel('Num spikes per stim per neuron');
        end
    end    
    
    linkaxes(ax_list,'xy');
    ylim([-0.05,0.5])
    
    xlim([0,4500])
    
    
    
%% plot latency of evoked spikes in model with and without synapses
    bin_edges = 0:0.25:5; % in ms, convert latency data to ms in for loop
    bin_centers = bin_edges(1:end-1) + mode(diff(bin_edges))/2;
    
    figure();
    diam = 2;
    ax_list = [];
    for i_amp = 1:numel(mdl_input_data.amp_list)
        ax_list(end+1) = subplot(1,4,i_amp); hold on;
        
        for i_type = 1:3
            if(i_type == 1)
                data = exp_resp_data;
                amp_list = [15,30,60,100];
                lat_mask = ones(size(exp_resp_data.latency));
                type_mask = ones(size(exp_resp_data.response_amp));
                temp_baseline_fr = exp_baseline_fr;
                color_use = 'k';
            elseif(i_type == 2)
                data = mdl_resp_data;
                amp_list = mdl_input_data.amp_list;
                lat_mask = mdl_resp_data.latency_diam == diam;
                type_mask = mdl_resp_data.diam == diam;
                temp_baseline_fr = 0;
                color_use = getColorFromList(1,1);
            else
                data = mdl_syn_resp_data;
                amp_list = mdl_input_data.amp_list;
                lat_mask = mdl_syn_resp_data.latency_diam == diam;
                type_mask = mdl_syn_resp_data.diam == diam;
                temp_baseline_fr = 0;
                color_use = getColorFromList(1,0);
            end
            
            lat_data = 1000*data.latency(data.latency_amp == i_amp & lat_mask);
            num_cells = sum(data.amp == amp_list(i_amp) & data.dist_from_stim < resp_input_data.max_lat_distance & type_mask);
            bin_data = histcounts(lat_data,bin_edges)/num_cells/mean(data.latency_num_stims(data.latency_amp == i_amp & lat_mask));
            bin_data = bin_data - temp_baseline_fr*mode(diff(bin_edges))/1000;
            plot(bin_centers,bin_data,'color',color_use,'linewidth',2);
        end
        
        formatForLee(gcf);
        set(gca,'fontsize',14);
        xlabel('Time after stim offset (ms)');
        ylabel('Num spikes per cell per stim');
    end
    
    linkaxes(ax_list,'xy');
%% compute space constants for each condition in model and experiment
%     max_dist = 1000*ceil(max([exp_resp_data.dist_from_stim; mdl_resp_data.dist_from_stim])/1000);
    max_dist = 2000;
    bin_size = 200; % um
    bin_edges = 200:bin_size:max_dist;

    mdl_space_constant_data = zeros(numel(mdl_input_data.diam_list),numel(mdl_input_data.cell_id_list),...
        mdl_input_data.num_clones,numel(mdl_input_data.amp_list),2); % with and without synapses
    mdl_space_constant_rsquare = size(mdl_space_constant_data);
    
    for i_diam = 1:numel(mdl_input_data.diam_list)
        for i_cell = 1:numel(mdl_input_data.cell_id_list)
            for i_clone = 1:mdl_input_data.num_clones
                for i_amp = 1:numel(mdl_input_data.amp_list)
                    mask = mdl_resp_data.diam == mdl_input_data.diam_list(i_diam) & ...
                        mdl_resp_data.cell_id == mdl_input_data.cell_id_list(i_cell) & ...
                        mdl_resp_data.clone_num == i_clone & ...
                        mdl_resp_data.amp == mdl_input_data.amp_list(i_amp);
                    
                    [mdl_space_constant_data(i_diam,i_cell,i_clone,i_amp,1),mdl_space_constant_rsquare(i_diam,i_cell,i_clone,i_amp,1)] = ...
                        getSpaceConstant(mdl_resp_data.response_amp(mask),mdl_resp_data.dist_from_stim(mask),bin_edges);
                    
                    
                    mask = mdl_syn_resp_data.diam == mdl_input_data.diam_list(i_diam) & ...
                        mdl_syn_resp_data.cell_id == mdl_input_data.cell_id_list(i_cell) & ...
                        mdl_syn_resp_data.clone_num == i_clone & ...
                        mdl_syn_resp_data.amp == mdl_input_data.amp_list(i_amp);
                    
                    [mdl_space_constant_data(i_diam,i_cell,i_clone,i_amp,2),mdl_space_constant_rsquare(i_diam,i_cell,i_clone,i_amp,2)] = ...
                        getSpaceConstant(mdl_syn_resp_data.response_amp(mask),mdl_syn_resp_data.dist_from_stim(mask),bin_edges);
                end
            end
        end
    end    
% get experiment space constant data
    % go through each unique channel and monkey combo and get a space constant
    unique_chan_monk = unique([exp_resp_data.monkey, exp_resp_data.chan_stim],'rows');
    unique_amps = unique(exp_resp_data.amp);
    exp_space_constant_data = nan(size(unique_chan_monk,1),numel(unique_amps));
    exp_space_constant_rsquare = nan(size(exp_space_constant_data));
    
    for i_cond = 1:size(unique_chan_monk,1)
        for i_amp = 1:numel(unique_amps)
            mask = exp_resp_data.monkey == unique_chan_monk(i_cond,1) & exp_resp_data.chan_stim == unique_chan_monk(i_cond,2) & ...
                exp_resp_data.amp == unique_amps(i_amp);
            
            [exp_space_constant_data(i_cond,i_amp),exp_space_constant_rsquare(i_cond,i_amp)] = getSpaceConstant(exp_resp_data.response_amp(mask),exp_resp_data.dist_from_stim(mask),bin_edges);
        end
    end
    
    
%% plot space constants -- model at each diameter (with and without synapses) and experiment on one plot
    plot_violin = 0; % else boxplot
    min_rsquare = 0.1;
    
    figure();
    ax_list = [];
    for i_amp = 1:4
        ax_list(i_amp) = subplot(1,4,i_amp); hold on;
        data_all = [];
        rsquare_all = [];
        group_all = [];
        % diam 1, 2 ,3 with and without synapses
        
        for i_type = 1:3
            num_data = size(mdl_space_constant_data,2)*size(mdl_space_constant_data,3);
            if(plot_violin)
                data_all(end+1:end+num_data,1) = reshape(mdl_space_constant_data(i_type,:,:,i_amp,1),num_data,1); % no synapses
                rsquare_all(end+1:end+num_data,1) = reshape(mdl_space_constant_rsquare(i_type,:,:,i_amp,1),num_data,1);
                group_all(end+1:end+num_data,1) = (2*i_type - 1)*ones(num_data,1);

                data_all(end+1:end+num_data,1) = reshape(mdl_space_constant_data(i_type,:,:,i_amp,2),num_data,1); % synapses
                group_all(end+1:end+num_data,1) = (2*i_type)*ones(num_data,1);
                rsquare_all(end+1:end+num_data,1) = reshape(mdl_space_constant_rsquare(i_type,:,:,i_amp,2),num_data,1);
            else % boxplot
                boxplot_params = [];
                boxplot_params.use_same_color_for_all = 1; % omits median color
                boxplot_params.master_color = getColorFromList(1,1);
                boxplot_params.median_color = 'k';
                boxplot_params.box_width = 0.25;
        
                % no synapse
                data = reshape(mdl_space_constant_data(i_type,:,:,i_amp,1),num_data,1);
                rsquare = reshape(mdl_space_constant_rsquare(i_type,:,:,i_amp,1),num_data,1);
                
                boxplot_wrapper(i_type-0.15,data(rsquare > min_rsquare),boxplot_params);
                
                % synapse
                boxplot_params.master_color = getColorFromList(1,0);
                data = reshape(mdl_space_constant_data(i_type,:,:,i_amp,2),num_data,1);
                rsquare = reshape(mdl_space_constant_rsquare(i_type,:,:,i_amp,2),num_data,1);
                
                boxplot_wrapper(i_type+0.15,data(rsquare>min_rsquare),boxplot_params);
            end
        end
        
        % experiment
        if(plot_violin)
            num_data = size(exp_space_constant_data,1);
            data_all(end+1:end+num_data,1) = exp_space_constant_data(:,i_amp);
            rsquare_all(end+1:end+num_data,1) = exp_space_constant_rsquare(:,i_amp);
            group_all(end+1:end+num_data,1) = 7;
            
            violin_plot(data_all,group_all);
        else
            boxplot_params = [];
            boxplot_params.use_same_color_for_all = 1; % omits median color
            boxplot_params.master_color = 'k';
            boxplot_params.median_color = 'k';
            boxplot_params.box_width = 0.25;
            
            boxplot_wrapper(4, exp_space_constant_data(exp_space_constant_rsquare(:,i_amp) > min_rsquare,i_amp),boxplot_params);
        end
        
        % format plot
        formatForLee(gcf);
        set(gca,'fontsize',14);
        xlim([0.5,4.5])
        if(i_amp == 1)
            ylabel('Space Constant (\mum)');
        end
    end

    linkaxes(ax_list,'xy');


    








