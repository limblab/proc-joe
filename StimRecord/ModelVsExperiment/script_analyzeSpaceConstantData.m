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
    
    resp_input_data.spike_window = [0,5]/1000; % s
    resp_input_data.sub_baseline = 1; % doesn't apply to model
    
    resp_input_data.is_model = 0;
    [exp_resp_data] = getResponseAndDistanceData(exp_array_data, resp_input_data);

    resp_input_data.is_model = 1;
    [mdl_resp_data] = getResponseAndDistanceData(mdl_array_data, resp_input_data);

%% plot response vs distance -- use small bins for experiment and model
% across all neurons for experiment
% across all clones for model
    max_dist = 1000*ceil(max([exp_resp_data.dist_from_stim; mdl_resp_data.dist_from_stim])/1000);
    bin_size = 100; % um
    bin_edges = 0:bin_size:max_dist;
    bin_centers = bin_edges(1:end-1) + mode(diff(bin_edges))/2;
    [~,~,exp_bin_idx] = histcounts(exp_resp_data.dist_from_stim,bin_edges);
    [~,~,mdl_bin_idx] = histcounts(mdl_resp_data.dist_from_stim,bin_edges);
    mdl_diam = 2; % 2 is monkey
    
    figure(); hold on;
    % experiment and model activations across all neurons/cell types/clones
    % (one diameter for model)
    fits = {}; gofs = {};
    ax_list = [];
    for i_type = 1:2
%         ax_list(end+1,1) = subplot(1,2,i_type); hold on;
        if(i_type == 1)
            resp_data = exp_resp_data;
            bin_idx = exp_bin_idx;
            type_mask = ones(size(exp_resp_data.response_amp));
            marker = '.';
            markersize = 20;
            linestyle = '--';
        else
            resp_data = mdl_resp_data;
            bin_idx = mdl_bin_idx;
            type_mask = mdl_resp_data.diam == mdl_diam;
            marker = 's';
            markersize = 7;
            linestyle = ':';
        end
        unique_amps = unique(resp_data.amp);
        for i_amp = 1:2:numel(unique_amps)
            mask = resp_data.amp == unique_amps(i_amp) & resp_data.dist_from_stim > 0 & type_mask;
            resp_amp = []; bin_val = [];
            for i_bin = 1:numel(bin_centers)
                if(sum(mask & bin_idx == i_bin) > 0)
                    resp_amp(end+1,1) = mean(resp_data.response_amp(mask & bin_idx == i_bin));
                    bin_val(end+1,1) = bin_centers(i_bin);
                end
            end

            plot(bin_val,resp_amp,'linestyle','none','marker',marker,'color',getColorFromList(1,i_amp-1),'markersize',markersize);

            % get exponential fit
            [fits{i_type,i_amp},gofs{i_type,i_amp}] = fit(bin_val,resp_amp,'a*exp(-x/b)+c','Lower',[0,0,-0.2],'Upper',[1,8000,0.2],'StartPoint',[0.03,-0.005,0]);
            plot(bin_centers,feval(fits{i_type,i_amp},bin_centers'),linestyle,'color',getColorFromList(1,i_amp-1),'linewidth',1.25)
        end
    end    
    
%     linkaxes(ax_list,'xy');
    
%% compute space constants for each condition in experiment and model

%% plot space constants

