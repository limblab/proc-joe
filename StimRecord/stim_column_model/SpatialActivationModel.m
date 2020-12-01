%% makes spatial activation prediction for a Histed brain and a Stoney brain
% so that analysis from Karthik's biophysical model can be put in
% perspective

% all we do is generate a set of neurons in a cortical column, then choose
% which neurons are activated based on rules (either Histed or Stoney
% rule)

    column_size = [2000,400,400]/1000; % mm, height, width, depth
    num_neurons = 1000;

    neuron_locs = generateCorticalColumn(column_size, num_neurons);

% sample activated population at different amplitudes using either a Stoney
% rule (square root activation based on distance) or a Histed rule (fill in
% volume, but volume size does not change)
    % stoney rule
    stoney_input_data = [];
    stoney_input_data.rule = 'Stoney';
    stoney_input_data.amp = 100;
    stoney_input_data.k =  1292; % uA/(mm^2)
    stoney_input_data.locs = neuron_locs;
    
    stoney_act_mask = simulateActivatedPopulation(stoney_input_data);

    % histed rule
    histed_input_data = [];
    histed_input_data.rule = 'Histed';
    histed_input_data.amp = 100;
    histed_input_data.max_amp = 100;
    histed_input_data.min_amp = 0;
    histed_input_data.min_p = 0;
    histed_input_data.max_p = 0.4;
    histed_input_data.locs = neuron_locs;
    
    histed_act_mask = simulateActivatedPopulation(histed_input_data);
    
% look at metrics for both data sets



