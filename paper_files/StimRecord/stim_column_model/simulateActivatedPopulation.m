function [act_mask] = simulateActivatedPopulation(input_data)

% this simulates activated population (input_data.neuron_locs) based on a
% specified set of rules

for i_amp = 1:numel(input_data.amp)
    amp = input_data.amp(i_amp);

    if(strcmpi(input_data.rule,'Stoney')==1)
        % probability of activation depends on distance from stim elec
        % (0,0,0). If within max radius: activated, else: not activated.
        
        max_r = sqrt(amp/input_data.k);
        r_data = sqrt(sum(input_data.locs.^2,2));
        
        act_mask(:,i_amp) = r_data < max_r;
        
    elseif(strcmpi(input_data.rule,'Histed')==1)
        % probability of activation depends on amplitude -- higher amp
        % means all neurons have higher probability of activation
       
        
        p = (amp-input_data.min_amp)/(input_data.max_amp-input_data.min_amp)*(input_data.max_p-input_data.min_p) + input_data.min_p;
        p = min(1,max(p,0)); % restric p to [0,1]
        
        act_mask(:,i_amp) = rand(size(input_data.locs,1),1) < p;
        
    else
        error('rule not implemented (yet)');
        
    end


end