% get model data -- stim channel responses for paired pulse data
    mdl_input_data = [];
    mdl_input_data.folderpath = 'C:\Users\Joseph Sombeck\Box\Miller-Grill_S1-stim\ModelData\TemporalResponse\SinglePulse\';
    mdl_input_data.amp_list = [15,30,50,100];
    mdl_input_data.IPI_list = [-1];
    mdl_input_data.gaba_ratio_list = [0.2,0.4,0.6,0.8,1.0];
    mdl_input_data.get_axon_dendrite_locs = 0;
    mdl_input_data.cell_id_list = [16];
    mdl_input_data.num_clones = 1;
    mdl_input_data.stim_times = (1000:500:(1000+500*9)); % ms, just the first one
    mdl_input_data.wave_length = 0.; % ms
    mdl_input_data.stim_window = [-200,600]; % ms around stim
    mdl_input_data.get_IPIs = 1;
    mdl_input_data.get_synapses = 0;
    
    mdl_input_data.is_train = 0;
    mdl_input_data.is_single_pulse = 1;
    % cell_id=6:10 %L23 PC, clones 1-5
    % cell_id=11:15 %L4 LBC, clones 1-5
    % cell_id=16:20 %L5 PC, clones 1-5
    % cell_id=21:25 %L6 PC, clones 1-5

    [mdl_data_all,mdl_array_data_all,mdl_mask_data_all] = getModelDoublePulseData(mdl_input_data);
    mdl_array_data = mdl_array_data_all; mdl_mask_data = mdl_mask_data_all; % use same proportion of cells as default, can resample below
    
    mdl_dists = getModelDistances(mdl_array_data);


%% get experimental data -- stim channel response
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
    exp_array_data = rebinArrayData(exp_array_data,1);

%% plot raster of example neuron(s)

    raster_input_data = [];
    raster_input_data.x_lim = [-100,500]; % ms
    raster_input_data.cond_list = [1,7,8,6];
    raster_input_data.marker_style = 'line'; % line is the correct way, but much slower
    
%     for exp_idx = 6 % 6
%         raster_input_data.is_model = 0;
%         plotModelExpDoublePulseRaster(exp_array_data{exp_idx},raster_input_data);
%     end
    for mdl_idx = [1:10]
        raster_input_data.is_model = 1;
        raster_input_data.cond_list = [1,2,3];
        plotModelExpDoublePulseRaster(mdl_array_data{mdl_idx},raster_input_data);
    end
    
    
    
%% compare response to first pulse to response to second pulse
    resp_input_data.spike_window = [-0,5]; % ms post stim to count spikes
    
    [mdl_array_data] = getResponseToEachPulse(mdl_array_data,resp_input_data);
    exp_array_data = getResponseToEachPulse(exp_array_data,resp_input_data);
    
    % 6,7,8 = 200, 10 20 ms for exp
    % 1,2,3 = 10,20,300 ms for model
%% plot response to first pulse against response to second pulse for each IPI

    figure(); hold on;
    for i_unit = 1:50%numel(mdl_array_data)
        plot(mdl_array_data{i_unit}.mean_resp_each_pulse{1}(1),mdl_array_data{i_unit}.mean_resp_each_pulse{1}(2),'r.')
%         plot(mdl_array_data{i_unit}.mean_resp_each_pulse{2}(1),mdl_array_data{i_unit}.mean_resp_each_pulse{2}(2),'b.')
%         plot(mdl_array_data{i_unit}.mean_resp_each_pulse{3}(1),mdl_array_data{i_unit}.mean_resp_each_pulse{3}(2),'k.')
    end
    
    for i_unit = 1:numel(exp_array_data)
        if(numel(exp_array_data{i_unit}.mean_resp_each_pulse) > 6)
            plot(exp_array_data{i_unit}.mean_resp_each_pulse{7}(1),exp_array_data{i_unit}.mean_resp_each_pulse{7}(2),'r.')
%             plot(exp_array_data{i_unit}.mean_resp_each_pulse{8}(1),exp_array_data{i_unit}.mean_resp_each_pulse{8}(2),'b.')
%             plot(exp_array_data{i_unit}.mean_resp_each_pulse{6}(1),exp_array_data{i_unit}.mean_resp_each_pulse{6}(2),'k.')
        end
    end
    
    plot([0,300],[0,300],'k--')
%% compute inhibition duration and plot across amplitudes

    inhib_input_data = [];
    inhib_input_data.pre_window = [-80,-15]/1000; % s 
    inhib_input_data.post_window = [0,220]/1000; % s
    inhib_input_data.bin_window = [inhib_input_data.pre_window(1),490/1000];
    inhib_input_data.max_time_start = 60/1000; % s
    inhib_input_data.bin_size = 5/1000; % s
    inhib_input_data.kernel_length = 2;
    inhib_input_data.blank_time = [-5,10]/1000; % s
    
    inhib_input_data.num_consec_bins = 2;
    inhib_input_data.IPI_list = [-1,200,10,20];
    inhib_input_data.cond_list = [1,6,6,7,8]; % single, 200 ms (first pulse), 200 ms (second pulse), 10 ms, 20 ms
    
    
    exp_inhib_data = getInhibitionDurationDoublePulseWrapper(exp_array_data,inhib_input_data);

  
    f=figure(); hold on
    subplot(1,2,1); hold on
    plot(exp_inhib_data.inhib_dur(:,1),exp_inhib_data.inhib_dur(:,3),'.','markersize',20)
    plot(exp_inhib_data.inhib_dur(:,1),exp_inhib_data.inhib_dur(:,4),'.','markersize',20)
    plot(exp_inhib_data.inhib_dur(:,1),exp_inhib_data.inhib_dur(:,5),'.','markersize',20)
    plot([0,0.25],[0,0.25],'k--')
    
    subplot(1,2,2); hold on;
    plot(exp_inhib_data.inhib_dur(:,1),exp_inhib_data.inhib_dur(:,2),'.','markersize',20)
    plot(exp_inhib_data.inhib_dur(:,2),exp_inhib_data.inhib_dur(:,3),'.','markersize',20)
    plot([0,0.25],[0,0.25],'k--')


    
%% sanity check code
    figure(); hold on
    idx = 20;

    for i=1:size(exp_inhib_data.inhib_dur,2)
        x_data = inhib_input_data.bin_window(1):inhib_input_data.bin_size:inhib_input_data.bin_window(2);
        x_data = x_data(1:end-1) + mode(diff(x_data))/2;
        plot(x_data,squeeze(exp_inhib_data.filtered_PSTH(idx,i,:))','color',getColorFromList(1,i-1),'linewidth',1.5)
    end
    plot([x_data(1),x_data(end)],exp_inhib_data.threshold(idx,1)+[0,0],'k--','linewidth',1)

    exp_inhib_data.inhib_dur(idx,:)