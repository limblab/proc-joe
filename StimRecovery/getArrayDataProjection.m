function [proj_data,arrayData] = getArrayDataProjection(arrayData, dpca_info, bin_size, sqrt_transform,do_pot_null)

    % outputs projections onto dPCs during baseline. 
    
    % rebin array data
    arrayData = rebinArrayData(arrayData,bin_size);
    
    % get baseline spike matrix
    num_bins = numel(arrayData{1}.binEdges{1}) - 1;
    num_stims = arrayData{1}.numStims(1);
    spikes = zeros(numel(arrayData),num_stims, num_bins);
    
    for i_arr = 1:numel(arrayData)
        for i_stim = 1:num_stims
            spikes(i_arr,i_stim,:) = histcounts(arrayData{i_arr}.spikeTrialTimes{1}(arrayData{i_arr}.stimData{1} == i_stim),...
                arrayData{i_arr}.binEdges{1}/1000);
        end
    end
    
    if(sqrt_transform)
        spikes = sqrt(spikes);
    end    
    
    
    % project baseline data onto dPCs
    if(do_pot_null == 0)
        proj_data = zeros(dpca_info.num_comps,num_stims,num_bins);
        for i_stim = 1:num_stims
            spikes_2D = squeeze(spikes(:,i_stim,:));
            spikes_2D = spikes_2D - dpca_info.mu';
            proj_data_temp_2D = dpca_info.W'*spikes_2D;
            proj_data(:,i_stim,:) = reshape(proj_data_temp_2D,size(proj_data_temp_2D,1),1,size(proj_data_temp_2D,2));
        end
    else
        % project spikes down, then project spikes onto potent and null
        % spaces
        pot_proj = zeros(size(dpca_info.V_potent,2),num_stims,num_bins);
        null_proj = zeros(size(dpca_info.V_null,2),num_stims,num_bins);
        num_dims = size(dpca_info.V_potent,1);
        for i_stim = 1:num_stims
            spikes_2D = squeeze(spikes(:,i_stim,:));
            spikes_2D = spikes_2D - dpca_info.mu_in';
            proj_data_temp_2D = dpca_info.w_in(:,1:num_dims)'*spikes_2D;
            pot_data_temp = dpca_info.V_potent'*proj_data_temp_2D;
            null_data_temp = dpca_info.V_null'*proj_data_temp_2D;
            
            pot_proj(:,i_stim,:) = reshape(pot_data_temp,size(pot_data_temp,1),1,size(pot_data_temp,2));
            null_proj(:,i_stim,:) = reshape(null_data_temp,size(null_data_temp,1),1,size(null_data_temp,2));
        end
        
        proj_data.pot_proj = pot_proj;
        proj_data.null_proj = null_proj;
    end
       
    
    
end