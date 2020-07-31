%% santize array data
for i_unit = 1:numel(exp_array_data)
    exp_array_data{i_unit}.spikeTrialTimes = exp_array_data{i_unit}.spikeTrialTimes';
    exp_array_data{i_unit}.trial_num = exp_array_data{i_unit}.trial_num';
    exp_array_data{i_unit}.num_stims = exp_array_data{i_unit}.num_stims';
    exp_array_data{i_unit}.binCounts = exp_array_data{i_unit}.binCounts';
    exp_array_data{i_unit}.binEdges = exp_array_data{i_unit}.binEdges';
    exp_array_data{i_unit}.stimData = exp_array_data{i_unit}.stimData';
end


%% compute inhibition duration and plot across amplitudes

    inhib_input_data = [];
    inhib_input_data.pre_window = [-80,-5]/1000; % s 
    inhib_input_data.post_window = [0,220]/1000; % s
    inhib_input_data.bin_window = [inhib_input_data.pre_window(1),inhib_input_data.post_window(2)+10/1000];
    inhib_input_data.max_time_start = 40/1000; % s
    inhib_input_data.bin_size = 5/1000; % s
    inhib_input_data.kernel_length = 2;
    inhib_input_data.blank_time = [0,10]/1000; % s
    
    inhib_input_data.num_consec_bins = 2;
    inhib_input_data.IPI_list = [-1,200,10,20];
    inhib_input_data.cond_list = [1,6,7,8];
    
    
    exp_inhib_data = getInhibitionDurationDoublePulseWrapper(exp_array_data,inhib_input_data);

  