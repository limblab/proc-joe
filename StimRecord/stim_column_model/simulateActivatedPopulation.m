function [act_mask] = simulateActivatedPopulation(input_data)

% this simulates activated population (input_data.neuron_locs) based on a
% specified set of rules

    if(strcmpi(input_data.rule,'Stoney')==1)
        % probability of activation depends on distance from stim elec
        % (0,0,0). If within max radius: activated, else: not activated.
        
        max_r = sqrt(input_data.amp/input_data.k);
        r_data = sqrt(sum(input_data.locs.^2,2));
        
        act_mask = r_data < max_r;
        
    elseif(strcmpi(input_data.rule,'Histed')==1)
        % probability of activation depends on amplitude -- higher amp
        % means all neurons have higher probability of activation
       
        
        p = (input_data.amp-input_data.min_amp)/(input_data.max_amp-input_data.min_amp)*(input_data.max_p-input_data.min_p) + input_data.min_p;
        p = min(1,max(p,0)); % restric p to [0,1]
        
        act_mask = rand(size(input_data.locs,1),1) < p;
        
    else
        error('rule not implemented (yet)');
        
    end


end