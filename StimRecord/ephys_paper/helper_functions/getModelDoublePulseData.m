function [mdl_data_all, array_data, mask_data] = getModelDoublePulseData(input_data)
% input data contains f
    %folderpath: data directory
    %IPI_list : inter-pulse interval in ms
    %amp_list : amplitude of stimulation
    %get_axon_dendrite_locs : bool, do you want locations of all axon or
            %dendrites?
    %cell_id_list : 6,11,16,21 corresponds to type of cell
    %num_clones : number of clones per cell type
    %gaba_ratio_list : list of gaba_B/gaba_A param
    
    curr_fpath = cd;
    had_noise_file = 0;
    
    gaba_list_all = []; amp_all = []; cell_id_all = []; clone_all = []; IPI_all = [];
    mdl_data_all = [];
    mdl_data_all.IPIs = [];
    cond_counter = 1;
    for i_gaba = 1:numel(input_data.gaba_ratio_list)
        for i_cell_id = 1:numel(input_data.cell_id_list)
            for i_clone = 1:input_data.num_clones
                for i_amp = 1:numel(input_data.amp_list)
                    for i_IPI = 1:numel(input_data.IPI_list)
                        gaba_list_all(end+1,1) = input_data.gaba_ratio_list(i_gaba);
                        amp_all(end+1,1) = input_data.amp_list(i_amp);
                        cell_id_all(end+1,1) = input_data.cell_id_list(i_cell_id);
                        clone_all(end+1,1) = i_clone;

                        % make coord data
                        coord_data = [];
                        coord_data.gaba_ratio = input_data.gaba_ratio_list(i_gaba);
                        coord_data.amp = input_data.amp_list(i_amp);
                        coord_data.IPI = input_data.IPI_list(i_IPI);
                        coord_data.cell_id = input_data.cell_id_list(i_cell_id);
                        coord_data.clone_id = i_clone;
                        coord_data.get_axon_dendrite_locs = input_data.get_axon_dendrite_locs;
                        coord_data.folderpath = input_data.folderpath;

                        % get location and spike data
                        cd(input_data.folderpath)
                        location_data = getTrueCoordinates(coord_data);

                        data_soma = [];
                        if(input_data.is_train == 1)
                            load([input_data.folderpath,'data_soma',num2str(coord_data.cell_id + coord_data.clone_id - 1),'_amp_',num2str(input_data.amp_list(i_amp)),'uA_freq_',...
                                num2str(1000/input_data.IPI_list(i_IPI)),'Hz_gabab_',num2str(input_data.gaba_ratio_list(i_gaba)),'_syn_act_100%.mat']);
                        else
                            load([input_data.folderpath,'data_soma',num2str(coord_data.cell_id + coord_data.clone_id - 1),'_amp_',num2str(input_data.amp_list(i_amp)),'uA_pp_',...
                                num2str(input_data.IPI_list(i_IPI)),'ms_gabab_',num2str(input_data.gaba_ratio_list(i_gaba)),'_syn_act_100%.mat']);
                        end
                        
                        if(exist('PoissonNoise.mat') > 0)
                            load('PoissonNoise.mat');
                            had_noise_file = 1;
                        else
                            noise = [];
                        end
                        mdl_data_all(cond_counter).soma = location_data.soma;
                        mdl_data_all(cond_counter).spike_data = data_soma;
                        mdl_data_all(cond_counter).gaba_ratio = coord_data.gaba_ratio;
                        mdl_data_all(cond_counter).amp = coord_data.amp;
                        mdl_data_all(cond_counter).IPI = coord_data.IPI;
                        mdl_data_all(cond_counter).cell_id = coord_data.cell_id;
                        mdl_data_all(cond_counter).folderpath = input_data.folderpath;
                        mdl_data_all(cond_counter).get_axon_dendrite_locs = coord_data.get_axon_dendrite_locs;
                        mdl_data_all(cond_counter).noise = noise;

                        % get IPI distribution for each condition combined....
                        if(input_data.get_IPIs)
                            for i_soma = 1:numel(data_soma)
                                if(~isempty(data_soma(i_soma).times))
                                    mdl_data_all(cond_counter).IPIs(end+1:end+numel(data_soma(i_soma).times)-1) = diff(sort(data_soma(i_soma).times/1000));
                                end
                            end
                        end
                        % update counter
                        cond_counter = cond_counter + 1;
                    end
                end
            end
        end
    end    
    
    
    
    
    % make an array_data struct with each cell -- cell_id and soma idx
    % defines a cell
    array_data = {};
    mask_data.gaba_ratio = [];
    mask_data.IPI = [];
    mask_data.cell_id = [];
    mask_data.clone = [];
    
    cell_counter = 1;
    % for each cell type and diameter
    for i_gaba = 1:numel(input_data.gaba_ratio_list)
        for i_amp = 1:numel(input_data.amp_list)
            for i_cell_id = 1:numel(input_data.cell_id_list)
                for i_clone = 1:input_data.num_clones
                    IPI_idx = find(gaba_list_all == input_data.gaba_ratio_list(i_gaba) & ...
                        cell_id_all == input_data.cell_id_list(i_cell_id) & clone_all == i_clone & ...
                        amp_all == input_data.amp_list(i_amp));
                    
                    % for each cell, then go through IPIs
                    for i_soma = 1:numel(mdl_data_all(IPI_idx(1)).spike_data)
                        % setup masks which may be useful later
                        mask_data.gaba_ratio(cell_counter,1) = input_data.gaba_ratio_list(i_gaba);
                        mask_data.cell_id(cell_counter,1) = input_data.cell_id_list(i_cell_id);
                        mask_data.clone(cell_counter,1) = i_clone;
                        mask_data.amp(cell_counter,1) = input_data.amp_list(i_amp);

                        % populate array_data
                        array_data{cell_counter}.stimData = cell(1,numel(IPI_idx));
                        array_data{cell_counter}.spikeTrialTimes = cell(1,numel(IPI_idx));
                        array_data{cell_counter}.intrinsic_near_stim = cell(1,numel(IPI_idx));
                        array_data{cell_counter}.baseline_fr = 0;
                        for i_IPI = 1:numel(IPI_idx)
                            % bin spikes around each stimulation based on a
                            % provided window. Store spikes in spikeTrialTimes and
                            % stim num in stimData
                            % spike time is relative to stim offset
                            % also need to determine if intrinsically generated
                            % spike was near (some window) stimulation. Label as
                            % such
                            array_data{cell_counter}.intrinsic_near_stim{i_IPI} = zeros(numel(input_data.stim_times),1);
                            for i_stim = 1:numel(input_data.stim_times)
                                % determine if intrinsic noise occurs 0-5 ms around
                                % stim
                                if(had_noise_file)
                                    noise_to_stim = noise(i_soma).times-input_data.stim_times(i_stim);
                                    noise_to_stim(noise_to_stim < 0) = 10000; % set negative entries as large
                                    array_data{cell_counter}.intrinsic_near_stim{1,i_IPI}(i_stim) = min(noise_to_stim) < 5;
                                    array_data{cell_counter}.baseline_fr = 1000/(mean(diff(noise(i_soma).times))); % Hz
                                end
                                % get spike times around stim time
                                spike_mask = mdl_data_all(IPI_idx(i_IPI)).spike_data(i_soma).times > input_data.stim_times(i_stim)+input_data.stim_window(1) & ...
                                    mdl_data_all(IPI_idx(i_IPI)).spike_data(i_soma).times <= input_data.stim_times(i_stim)+input_data.stim_window(2);

                                array_data{cell_counter}.spikeTrialTimes{i_IPI}(end+1:end+sum(spike_mask),1) = mdl_data_all(IPI_idx(i_IPI)).spike_data(i_soma).times(spike_mask==1) - (input_data.stim_times(i_stim)+input_data.wave_length);

                                array_data{cell_counter}.stimData{i_IPI}(end+1:end+sum(spike_mask),1) = i_stim+zeros(sum(spike_mask),1);
                                
                                array_data{cell_counter}.PULSE_TIMES{i_IPI,1}{i_stim} = [0,input_data.IPI_list(i_IPI)/1000]; % done here to match experimental array data
                            end
                            % convert times to s
                            array_data{cell_counter}.spikeTrialTimes{i_IPI} = array_data{cell_counter}.spikeTrialTimes{i_IPI}/1000;

                            array_data{cell_counter}.STIM_PARAMETERS(i_IPI).amp1 = [mdl_data_all(IPI_idx(i_IPI)).amp];
                            
                            
                        end % end amp

                        % meta data : STIM_PARAMETERS, WAVEFORM_LIST, ID, location,
                        % monkey, numStims, diameter, cell_id
                        array_data{cell_counter}.monkey = 'model';
                        array_data{cell_counter}.loc = mdl_data_all(IPI_idx(1)).soma(i_soma).coord;
                        array_data{cell_counter}.waveform_list = 1:1:numel(IPI_idx);

                        array_data{cell_counter}.numStims = numel(input_data.stim_times)+zeros(1,numel(IPI_idx));
                        array_data{cell_counter}.diam = input_data.gaba_ratio_list(i_gaba);
                        array_data{cell_counter}.cell_id = input_data.cell_id_list(i_cell_id);
                        array_data{cell_counter}.clone_num = i_clone;

                        
                        % update counter after each soma
                        cell_counter = cell_counter + 1;
                    end % end soma
                end
            end
        end
    end
    
       
    % go back to original folder
    cd(curr_fpath);

end




function [ output_data ] = getTrueCoordinates( input_data )
% input data contains f
    %diam : 1, 2, 3 for diameter
    %amp : amplitude of stimulation
    %get_axon_dendrite_locs : bool, do you want locations of all axon or
            %dendrites?
    %cell_id : 6,11,16,21 corresponds to type of cell
    
% outputs location data for each cell
    cell_id = input_data.cell_id + input_data.clone_id - 1;
    
    load('realx.dat') %x-coordinate (in um) of randomly seeded soma within shpherical volume
    load('realy.dat') %y-coordinate (in um) of randomly seeded soma within shpherical volume
    load('realz.dat') %z-coordinate (in um) of randomly seeded soma within shpherical volume

    load('realang.dat') %Angle (in radian) of random rotation of neuron in azimuthal direction 

    intx=load(['intx_' num2str(cell_id) '.dat']); %x-coordinates (in um) of all neural section (dendrite, soma and axon) of a single neuron 
    inty=load(['inty_' num2str(cell_id) '.dat']); %y-coordinates (in um) of all neural section (dendrite, soma and axon) of a single neuron
    intz=load(['intz_' num2str(cell_id) '.dat']); %z-coordinates (in um) of all neural section (dendrite, soma and axon) of a single neuron

    soma_coord=load(['soma_coord_' num2str(cell_id) '.dat']); %x,y,z coordinates (in um) of somatic section (of a single neuron) 

    %To get true location of somatic section of each neuron within the spherical
    %volume
    soma = [];
    for k=1:length(realx)

        net_ad_x=((realx(k)+soma_coord(1))*cos(realang(k)))-((realz(k)+soma_coord(3))*sin(realang(k)))-realx(k);
        net_ad_y=realy(k)+soma_coord(2)-realy(k);
        net_ad_z=((realx(k)+soma_coord(1))*sin(realang(k)))+((realz(k)+soma_coord(3))*cos(realang(k)))-realz(k);

        x1=((realx(k)+soma_coord(1))*cos(realang(k)))-((realz(k)+soma_coord(3))*sin(realang(k)));
        y1=realy(k)+soma_coord(2);
        z1=((realx(k)+soma_coord(1))*sin(realang(k)))+((realz(k)+soma_coord(3))*cos(realang(k)));

        soma(k).coord(1)=x1-net_ad_x;
        soma(k).coord(2)=y1-net_ad_y;
        soma(k).coord(3)=z1-net_ad_z;
        soma(k).meta = 'xyz';

    end
    
    % package outputs
    output_data = [];
    
    output_data.soma = soma;
    output_data.folderpath = input_data.folderpath;
    output_data.amp = input_data.amp;
    output_data.cell_id = input_data.cell_id;
    output_data.get_axon_dendrite_locs = input_data.get_axon_dendrite_locs;
    
    
end

