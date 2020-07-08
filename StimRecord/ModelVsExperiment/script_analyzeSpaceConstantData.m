% get model data -- space constant DONT RUN, LOAD .MAT FILE INSTAED
    mdl_input_data = [];
    mdl_input_data.folderpath = 'C:\Users\Joseph Sombeck\Box\Miller-Grill_S1-stim\ModelData\SpaceConstant\';
    mdl_input_data.diam_list = [1,2,3];
    mdl_input_data.amp_list = [15,30,50,100];
    mdl_input_data.get_axon_dendrite_locs = 0;
    mdl_input_data.cell_id_list = [1,6,11,16,21];
    mdl_input_data.num_clones = 5;
    mdl_input_data.stim_times = (201) - 0.453; % ms
    mdl_input_data.wave_length = 0.453; % ms
    mdl_input_data.stim_window = [-100,400]; % ms around stim
    mdl_input_data.get_IPIs = 0;
    
    % cell_id=6:10 %L23 PC, clones 1-5
    % cell_id=11:15 %L4 LBC, clones 1-5
    % cell_id=16:20 %L5 PC, clones 1-5
    % cell_id=21:25 %L6 PC, clones 1-5
    
    mdl_input_data.get_synapses = 1;
    [mdl_data_all,mdl_syn_array_data,mdl_syn_mask_data] = getModelStimChannelData(mdl_input_data);    
    mdl_input_data.get_synapses = 0;
    [mdl_data_all,mdl_array_data,mdl_mask_data] = getModelStimChannelData(mdl_input_data);    
    
%% get experiment data -- space constant
    input_data.home_computer = 1;
    exp_array_data = getExperimentSpaceConstantData(input_data);
    exp_array_data = adjustArrayDataSpikeTimes(exp_array_data, 0.453/1000); % stim pulse length
    exp_array_data = getBaselineFiringRate(exp_array_data,[-25,-5]/1000); % window relative to stim onset
    

%% get response data for both model and experiment
    resp_input_data = []; 
    resp_input_data.home_computer = input_data.home_computer;
    resp_input_data.monkey_list = {'Han','Duncan'};
    
    
    resp_input_data.sub_baseline = 1; % doesn't apply to model
    
    resp_input_data.is_model = 0;
    resp_input_data.spike_window = [1,5]/1000; % s
    [exp_resp_data] = getResponseAndDistanceData(exp_array_data, resp_input_data);

    resp_input_data.is_model = 1;
    resp_input_data.spike_window = [0,5]/1000; % s
    [mdl_resp_data] = getResponseAndDistanceData(mdl_array_data, resp_input_data);

    resp_input_data.is_model = 1;
    resp_input_data.spike_window = [0,5]/1000; % s
    [mdl_syn_resp_data] = getResponseAndDistanceData(mdl_syn_array_data, resp_input_data);
    
    % adjust response amplitude because there's no intrinsic activity
    mdl_syn_resp_data.response_amp = mdl_syn_resp_data.response_amp*0.5;
    mdl_resp_data.response_amp = mdl_resp_data.response_amp*0.5;
    
%% plot response vs distance -- use bins for experiment and model, model + synapses
% across all neurons for experiment
% across all clones for model
    max_dist = 1000*ceil(max([exp_resp_data.dist_from_stim; mdl_resp_data.dist_from_stim])/1000);
    bin_size = 400; % um
    bin_edges = 200:bin_size:max_dist;
    bin_centers = bin_edges(1:end-1) + mode(diff(bin_edges))/2;
    [~,~,exp_bin_idx] = histcounts(exp_resp_data.dist_from_stim,bin_edges);
    [~,~,mdl_bin_idx] = histcounts(mdl_resp_data.dist_from_stim,bin_edges);
    [~,~,mdl_syn_bin_idx] = histcounts(mdl_syn_resp_data.dist_from_stim,bin_edges);
    fit_bin_min = 375;
    
    mdl_diam = 2; % 2 is monkey diameter
    
    
    figure('Position',[2097 497 1640 420]); hold on;
    % experiment and model activations across all neurons/cell types/clones
    % (one diameter for model)
    fits = {}; gofs = {};
    ax_list = [];
    for i_type = 1:3
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
            type_mask = mdl_syn_resp_data.diam == mdl_diam;
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
    
    xlim([0,3000])
    
    
%% compare number of evoked spikes in model with and without synapses

    % 1 plot per diameter, 4 amplitudes and each cell type per plot
    figure();
    ax_list = [];
    for i_diam = 1:numel(mdl_input_data.diam_list)
        ax_list(i_diam) = subplot(1,3,i_diam); hold on;
        num_spikes = nan(numel(mdl_input_data.cell_id_list), numel(mdl_input_data.amp_list),2); % no synapses, synapses
        for i_cell_id = 1:numel(mdl_input_data.cell_id_list)
            for i_amp = 1:numel(mdl_input_data.amp_list)
                mask = mdl_resp_data.diam == mdl_input_data.diam_list(i_diam) & mdl_resp_data.amp == mdl_input_data.amp_list(i_amp) & ...
                    mdl_resp_data.cell_id == mdl_input_data.cell_id_list(i_cell_id);
                num_spikes(i_cell_id,i_amp,1) = sum(mdl_resp_data.response_amp(mask));
                
                mask = mdl_syn_resp_data.diam == mdl_input_data.diam_list(i_diam) & mdl_syn_resp_data.amp == mdl_input_data.amp_list(i_amp) & ...
                    mdl_syn_resp_data.cell_id == mdl_input_data.cell_id_list(i_cell_id);
                num_spikes(i_cell_id,i_amp,2) = sum(mdl_syn_resp_data.response_amp(mask));
            end
            
            % plot data for cell id
            plot(mdl_input_data.amp_list,squeeze(num_spikes(i_cell_id,:,1)),'-','color',getColorFromList(1,i_cell_id-1),'linewidth',2)
            plot(mdl_input_data.amp_list,squeeze(num_spikes(i_cell_id,:,2)),'--','color',getColorFromList(1,i_cell_id-1),'linewidth',2)
        end
        
    end

    linkaxes(ax_list,'xy');
    
    

    
%% compute space constants for each condition in model and experiment
    max_dist = 1000*ceil(max([exp_resp_data.dist_from_stim; mdl_resp_data.dist_from_stim])/1000);
    bin_size = 100; % um
    bin_edges = 0:bin_size:max_dist;

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
            data_all(end+1:end+num_data,1) = reshape(mdl_space_constant_data(i_type,:,:,i_amp,1),num_data,1); % no synapses
            rsquare_all(end+1:end+num_data,1) = reshape(mdl_space_constant_rsquare(i_type,:,:,i_amp,1),num_data,1);
            group_all(end+1:end+num_data,1) = (2*i_type - 1)*ones(num_data,1);
            
            data_all(end+1:end+num_data,1) = reshape(mdl_space_constant_data(i_type,:,:,i_amp,2),num_data,1); % synapses
            group_all(end+1:end+num_data,1) = (2*i_type)*ones(num_data,1);
            rsquare_all(end+1:end+num_data,1) = reshape(mdl_space_constant_rsquare(i_type,:,:,i_amp,2),num_data,1);
        end
        
        % experiment
        num_data = size(exp_space_constant_data,1);
        data_all(end+1:end+num_data,1) = exp_space_constant_data(:,i_amp);
        rsquare_all(end+1:end+num_data,1) = exp_space_constant_rsquare(:,i_amp);
        group_all(end+1:end+num_data,1) = 7;
        
        
        % rsquare_mask
        rsquare_mask = rsquare_all > 0.25 & data_all < 10000;
        boxplot(data_all(rsquare_mask),group_all(rsquare_mask));
    end

    linkaxes(ax_list,'xy');


    








