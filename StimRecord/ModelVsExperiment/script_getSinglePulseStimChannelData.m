% get model data -- stim channel response
    mdl_input_data = [];
    mdl_input_data.folderpath = 'C:\Users\Joseph Sombeck\Box\Miller-Grill_S1-stim\ModelData\StimChannelResponse\with_Intrinsic_Activity\';
    mdl_input_data.diam_list = [1,2,3];
    mdl_input_data.amp_list = [5,10,15,20,25,30,40,50,100];
    mdl_input_data.get_axon_dendrite_locs = 0;
    mdl_input_data.cell_id_list = [6,11,16,21];
    mdl_input_data.stim_times = 200:500:(200+500*9); % ms
    mdl_input_data.wave_length = 0.453; % ms
    mdl_input_data.stim_window = [-100,400]; % ms around stim
    
    % cell_id=6:10 %L23 PC, clones 1-5
    % cell_id=11:15 %L4 LBC, clones 1-5
    % cell_id=16:20 %L5 PC, clones 1-5
    % cell_id=21:25 %L6 PC, clones 1-5

    [mdl_data_all,mdl_array_data,mdl_mask_data] = getModelStimChannelData(mdl_input_data);
    mdl_dists = getModelDistances(mdl_array_data);
    
    
%% get experimental data -- stim channel response
    exp_input_data.home_computer = 1;
    [exp_data] = getExperimentStimChannelData(exp_input_data);
    exp_array_data = exp_data.array_data;
    
%% Spike times after stimulation, activation threshold
    activation_input_data.spike_window = [0,4]/1000;
    activation_input_data.remove_intrinsic = 1;
    
    activation_input_data.is_model = 0;
    exp_threshold_data = getActivationThreshold(exp_array_data,activation_input_data);
    
    activation_input_data.is_model = 1;
    mdl_threshold_data = getActivationThreshold(mdl_array_data,activation_input_data);
    
% plot activation threshold data -- NEED TO ORDER BY DIAMETER
% boxplots across all cell types for each diameter and experiment
% boxplots for each cell type and diameter
    figure();
    subplot(2,2,1) % across cell types, for each diameter and experiment
    for i_diam = 1:numel(mdl_input_data.diam_list) % model data
        boxplot_params = [];
        boxplot_params.use_same_color_for_all = 1;
        boxplot_params.master_color = getColorFromList(1,i_diam-1);
        
        data = mdl_threshold_data.thresholds(mdl_threshold_data.is_responsive & mask_data.diam == i_diam);
        boxplot_wrapper(i_diam, data, boxplot_params);
    end
    % experimental data
    boxplot_params = [];
    boxplot_params.use_same_color_for_all = 1;
    boxplot_params.master_color = 'k';

    data = exp_threshold_data.thresholds(exp_threshold_data.is_responsive==1);
    boxplot_wrapper(4, data, boxplot_params);
    
    % plot percent responsive all cell types for each diameter
    subplot(2,2,2); hold on
    for i_diam = 1:numel(mdl_input_data.diam_list) % model data
        keep_mask = mask_data.diam == mdl_input_data.diam_list(i_diam);
        perc_resp = sum(mdl_threshold_data.is_responsive(keep_mask))/sum(keep_mask);
            
        b=bar(i_diam, perc_resp);
        b.FaceColor = getColorFromList(1,i_diam-1);
    end
    % experimental data
    perc_resp = sum(exp_threshold_data.is_responsive)/numel(exp_threshold_data.is_responsive);
    b=bar(4, perc_resp);
    b.FaceColor = 'k';
    
    % boxplot for each cell type grouped by diameter
    subplot(2,2,3) 
    x_pos = 1;
    for i_diam = 1:numel(mdl_input_data.diam_list) % model data
        for i_cell = 1:numel(mdl_input_data.cell_id_list)
            boxplot_params = [];
            boxplot_params.use_same_color_for_all = 1;
            boxplot_params.master_color = getColorFromList(1,i_cell+2);

            data = mdl_threshold_data.thresholds(mdl_threshold_data.is_responsive & ...
                mask_data.diam == mdl_input_data.diam_list(i_diam) & mask_data.cell_id == mdl_input_data.cell_id_list(i_cell));
            
            boxplot_wrapper(x_pos, data, boxplot_params);
            x_pos = x_pos + 1;
        end
        if(i_diam < numel(mdl_input_data.diam_list))
            x_pos = x_pos + 2;
            plot([x_pos-1.5,x_pos-1.5],[0,100],'--','color',[0.5,0.5,0.5])
        end
    end
    xlim([0,x_pos])
    
    % plot percent responsive for each cell type grouped by diameter
    subplot(2,2,4); hold on
    x_pos = 1;
    for i_diam = 1:numel(mdl_input_data.diam_list) % model data
        for i_cell = 1:numel(mdl_input_data.cell_id_list)
            keep_mask = mask_data.diam == mdl_input_data.diam_list(i_diam) & mask_data.cell_id == mdl_input_data.cell_id_list(i_cell);
            perc_resp = sum(mdl_threshold_data.is_responsive(keep_mask))/sum(keep_mask);
            
            b=bar(x_pos, perc_resp);
            b.FaceColor = getColorFromList(1,i_cell+2);
            x_pos = x_pos + 1;
        end
        if(i_diam < numel(mdl_input_data.diam_list))
            x_pos = x_pos + 2;
            plot([x_pos-1.5,x_pos-1.5],[0,1],'--','color',[0.5,0.5,0.5])
        end
    end
    xlim([0,x_pos])
    
    
%% get binned spike times after stim
    