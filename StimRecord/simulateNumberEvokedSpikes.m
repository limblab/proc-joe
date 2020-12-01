function [output_data] = simulateNumberEvokedSpikes(input_data)
% this function simulates the number of evoked spikes in a volume based on
% an activation function and cell density. 

% input_data contains the following:
%     activation_func -- probability of activation for differnt distances (mm)
%     min_r; % mm
%     max_r % mm
%     dr % mm
%     density % cells/0.001 mm^3...


    r_data = input_data.min_r:input_data.dr:input_data.max_r;


    activation_data = max(input_data.activation_func(r_data),0); % no activation below 0
    
    vol_data = ((r_data+input_data.dr).^3) - (r_data.^3);

    slice_data = 4.3*pi.*vol_data.*activation_data.*(input_data.density/0.001);

    num_evoked_spikes = sum(slice_data);
    
    % package outputs
    output_data.num_evoked_spikes = num_evoked_spikes;
    output_data.slice_data = slice_data;
    output_data.r_data = r_data;
    
end