function [output_data] = makeNeurometricCurves(td,psych_data)

    % function gets data and plots firing rate per trial as a function of
    % bump direction relative to target axis
    % psych data is assumed to be of a single axis.
    
    % trial by num neurons matrix
    spike_data = zeros(numel(psych_data.trial_ids),size(td(1).LeftS1_unit_guide,1));
   
    % trial by 1 matrix
    bump_dirs = zeros(numel(psych_data.trial_ids),1);
    reached_0_deg = zeros(numel(psych_data.trial_ids),1);
    % build spike_data matrix and get list of bump dirs for each trial
    for trial_counter = 1:numel(psych_data.trial_ids) % for each trial
        trial_id = psych_data.trial_ids(trial_counter);
        
        % find entry in td
        td_idx = find(trial_id == [td.trial_id]);
        
        % build window based on bump time
        window = td(td_idx).idx_bumpTime + [0,400]; 
        
        % get firing rate for each neuron for that trial
        if(~isempty(td_idx))
            spike_data(trial_counter,:) = sum(td(td_idx).LeftS1_spikes(window(1):window(2),:),1)/(diff(window)*td(td_idx).bin_size);
        end
        
        % get bump dir for this trial (relative to target axis, 0-180
        % degrees)
%         if(td(td_idx).bumpDir > 180)
%             bump_dirs(trial_counter,1) = -1*(td(td_idx).bumpDir-360);
%         else
            bump_dirs(trial_counter,1) = td(td_idx).bumpDir;
%         end
        
        % determine if animal reached to the "0-deg" target
        
        % get if 0 deg target was answer or not
        if(bump_dirs(trial_counter,1) > 90)
            reached_0_deg(trial_counter,1) = 0;
        else
            reached_0_deg(trial_counter,1) = 1;
        end
        % flip if monkey failed
        if(td(td_idx).result == 'F')
            reached_0_deg(trial_counter,1) = ~reached_0_deg(trial_counter,1);
        end
        
    end

    % plot firing rate vs bump dir
    for neuron_idx = 1:size(spike_data,2)
        figure();
        plot(bump_dirs(reached_0_deg == 1),spike_data(reached_0_deg == 1,neuron_idx),'.','markersize',20);
        hold on
        plot(bump_dirs(reached_0_deg == 0),spike_data(reached_0_deg == 0,neuron_idx),'x','markersize',20);
    end
    
end