function [ output_data ] = getBiomimeticPattern( td, input_data )
% gets a biomimetic pattern that stimulates on the number of electrodes
% specified in input_data (divided by 2) based on the  PD of those
% electrodes and the dir specified in input data. The cos(PD-dir) determines
% the relative frequency or amplitude for each electrode. Nonbiomimetic patterns 
% are developed by shuffling the relative magnitudes across the number of
% electrodes

 % randomly choose well tuned channels in half of the circe
    well_tuned_chan_idx = find(input_data.well_tuned & cos(angleDiff(input_data.pd_elec,input_data.dir_model)) > 0);
    elec_idx = well_tuned_chan_idx(randperm(numel(well_tuned_chan_idx),min(numel(well_tuned_chan_idx),input_data.num_elec/2)));
    
    % randomly choose well tuned channels in the other half of the circle
    well_tuned_chan_idx = find(input_data.well_tuned & cos(angleDiff(input_data.pd_elec,input_data.dir_model)) <= 0);
    elec_idx = [elec_idx; well_tuned_chan_idx(randperm(numel(well_tuned_chan_idx),min(numel(well_tuned_chan_idx),input_data.num_elec/2)))];
    
    
    output_data.chan_num = td(1).LeftS1_unit_guide(elec_idx,1);
     
     
    output_data.biomimetic_freq_norm = cos(angleDiff(input_data.pd_elec(elec_idx),input_data.dir_model,1));
    output_data.biomimetic_freq_norm(output_data.biomimetic_freq_norm < 0) = 0;
     
    
    
    % adjust frequencies to match potential frequencies, if provided
    if(isfield(input_data,'freq_all_norm'))
        for i = 1:numel(output_data.biomimetic_freq_norm)
            if(output_data.biomimetic_freq_norm(i) == 1)
                output_data.biomimetic_freq_norm(i) = input_data.freq_all_norm(end);
            elseif(output_data.biomimetic_freq_norm(i) ~= 0)
                idx_small = find(output_data.biomimetic_freq_norm(i) > input_data.freq_all_norm,1,'last');
                dist = [output_data.biomimetic_freq_norm(i) - input_data.freq_all_norm(idx_small:idx_small+1)];
                if(isempty(idx_small))
                    output_data.biomimetic_freq_norm(i) = input_data.freq_all_norm(1);
                elseif(dist(1) < dist(2))
                    output_data.biomimetic_freq_norm(i) = input_data.freq_all_norm(idx_small);
                else
                    output_data.biomimetic_freq_norm(i) = input_data.freq_all_norm(idx_small+1);
                end
            end
        end
    end
    
    output_data.nonbiomimetic_freq_norm = output_data.biomimetic_freq_norm(randperm(numel(output_data.biomimetic_freq_norm)));
    output_data.dir = input_data.dir_model;
    
end

