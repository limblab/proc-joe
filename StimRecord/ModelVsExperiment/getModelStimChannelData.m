function [mdl_data_all, array_data, mask_data] = getModelStimChannelData(input_data)
% input data contains f
    %folderpath: data directory
    %diam_list : 1, 2, 3 for diameter
    %amp_list : amplitude of stimulation
    %get_axon_dendrite_locs : bool, do you want locations of all axon or
            %dendrites?
    %cell_id_list : 6,11,16,21 corresponds to type of cell
    curr_fpath = cd;
    had_noise_file = 0;
    
    diam_all = []; amp_all = []; cell_id_all = []; clone_all = [];
    mdl_data_all = [];
    mdl_data_all.IPIs = [];
    cond_counter = 1;
    for i_diam = 1:numel(input_data.diam_list)
        for i_cell_id = 1:numel(input_data.cell_id_list)
            for i_clone = 1:input_data.num_clones
                for i_amp = 1:numel(input_data.amp_list)
                    diam_all(end+1,1) = input_data.diam_list(i_diam);
                    amp_all(end+1,1) = input_data.amp_list(i_amp);
                    cell_id_all(end+1,1) = input_data.cell_id_list(i_cell_id);
                    clone_all(end+1,1) = i_clone;
                    
                    % make coord data
                    coord_data = [];
                    coord_data.diam = input_data.diam_list(i_diam);
                    coord_data.amp = input_data.amp_list(i_amp);
                    coord_data.cell_id = input_data.cell_id_list(i_cell_id);
                    coord_data.clone_id = i_clone;
                    coord_data.get_axon_dendrite_locs = input_data.get_axon_dendrite_locs;
                    coord_data.folderpath = input_data.folderpath;

                    % get location and spike data
                    data_path = [coord_data.folderpath, 'Diam',num2str(coord_data.diam),'_',num2str(coord_data.amp),'uA\'];
                    cd(data_path);
                    location_data = getTrueCoordinates(coord_data);
                    
                    data_soma = [];
                    try
                        if(input_data.get_synapses)
                            load(['data_soma',num2str(coord_data.cell_id + coord_data.clone_id - 1),'_syn']);
                        else
                            load(['data_soma',num2str(coord_data.cell_id + coord_data.clone_id - 1)]);
                        end
                    end
                    if(exist('PoissonNoise.mat') > 0)
                        load('PoissonNoise.mat');
                        had_noise_file = 1;
                    else
                        noise = [];
                    end
                    mdl_data_all(cond_counter).soma = location_data.soma;
                    mdl_data_all(cond_counter).axon = location_data.axon;
                    mdl_data_all(cond_counter).loc_all = location_data.loc_all;
                    mdl_data_all(cond_counter).spike_data = data_soma;
                    mdl_data_all(cond_counter).diam = coord_data.diam;
                    mdl_data_all(cond_counter).amp = coord_data.amp;
                    mdl_data_all(cond_counter).cell_id = coord_data.cell_id;
                    mdl_data_all(cond_counter).folderpath = data_path;
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
    
    
    
    
    % make an array_data struct with each cell -- cell_id and soma idx
    % defines a cell
    array_data = {};
    mask_data.diam = [];
    mask_data.cell_id = [];
    mask_data.clone = [];
    
    cell_counter = 1;
    % for each cell type and diameter
    for i_diam = 1:numel(input_data.diam_list)
        for i_cell_id = 1:numel(input_data.cell_id_list)
            for i_clone = 1:input_data.num_clones
                amp_idx = find(diam_all == input_data.diam_list(i_diam) & cell_id_all == input_data.cell_id_list(i_cell_id) & clone_all == i_clone);
                % for each cell, then go through amplitudes
                for i_soma = 1:numel(mdl_data_all(amp_idx(1)).spike_data)
                    % setup masks which may be useful later
                    mask_data.diam(cell_counter,1) = input_data.diam_list(i_diam);
                    mask_data.cell_id(cell_counter,1) = input_data.cell_id_list(i_cell_id);
                    mask_data.clone(cell_counter,1) = i_clone;
                    
                    % populate array_data
                    array_data{cell_counter}.stimData = cell(1,numel(amp_idx));
                    array_data{cell_counter}.spikeTrialTimes = cell(1,numel(amp_idx));
                    array_data{cell_counter}.intrinsic_near_stim = cell(1,numel(amp_idx));
                    array_data{cell_counter}.baseline_fr = 0;
                    for i_amp = 1:numel(amp_idx)
                        % bin spikes around each stimulation based on a
                        % provided window. Store spikes in spikeTrialTimes and
                        % stim num in stimData
                        % spike time is relative to stim offset
                        % also need to determine if intrinsically generated
                        % spike was near (some window) stimulation. Label as
                        % such
                        array_data{cell_counter}.intrinsic_near_stim{i_amp} = zeros(numel(input_data.stim_times),1);
                        for i_stim = 1:numel(input_data.stim_times)
                            % determine if intrinsic noise occurs 0-5 ms around
                            % stim
                            if(had_noise_file)
                                noise_to_stim = noise(i_soma).times-input_data.stim_times(i_stim);
                                noise_to_stim(noise_to_stim < 0) = 10000; % set negative entries as large
                                array_data{cell_counter}.intrinsic_near_stim{1,i_amp}(i_stim) = min(noise_to_stim) < 5;
                                array_data{cell_counter}.baseline_fr = 1000/(mean(diff(noise(i_soma).times))); % Hz
                            end
                            % get spike times around stim time
                            spike_mask = mdl_data_all(amp_idx(i_amp)).spike_data(i_soma).times > input_data.stim_times(i_stim)+input_data.stim_window(1) & ...
                                mdl_data_all(amp_idx(i_amp)).spike_data(i_soma).times <= input_data.stim_times(i_stim)+input_data.stim_window(2);

                            array_data{cell_counter}.spikeTrialTimes{i_amp}(end+1:end+sum(spike_mask),1) = mdl_data_all(amp_idx(i_amp)).spike_data(i_soma).times(spike_mask==1) - (input_data.stim_times(i_stim)+input_data.wave_length);

                            array_data{cell_counter}.stimData{i_amp}(end+1:end+sum(spike_mask),1) = i_stim+zeros(sum(spike_mask),1);
                        end
                        % convert times to s
                        array_data{cell_counter}.spikeTrialTimes{i_amp} = array_data{cell_counter}.spikeTrialTimes{i_amp}/1000;

                        array_data{cell_counter}.STIM_PARAMETERS(i_amp).amp1 = [mdl_data_all(amp_idx(i_amp)).amp];
                    end % end amp

                    % meta data : STIM_PARAMETERS, WAVEFORM_LIST, ID, location,
                    % monkey, numStims, diameter, cell_id
                    array_data{cell_counter}.monkey = 'model';
                    array_data{cell_counter}.loc = mdl_data_all(amp_idx(1)).soma(i_soma).coord;
                    array_data{cell_counter}.waveform_list = 1:1:numel(amp_idx);

                    array_data{cell_counter}.numStims = numel(input_data.stim_times)+zeros(1,numel(amp_idx));
                    array_data{cell_counter}.diam = input_data.diam_list(i_diam);
                    array_data{cell_counter}.cell_id = input_data.cell_id_list(i_cell_id);
                    array_data{cell_counter}.clone_num = i_clone;
                    
                    % update counter after each soma
                    cell_counter = cell_counter + 1;
                end % end soma
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
    
    load('realx.dat') %x-coordinate (in um) of randomly seeded soma within spherical volume
    load('realy.dat') %y-coordinate (in um) of randomly seeded soma within spherical volume
    load('realz.dat') %z-coordinate (in um) of randomly seeded soma within spherical volume

    load('realang.dat') %Angle (in radian) of random rotation of neuron in azimuthal direction 

    soma_coord=load(['soma_coord_' num2str(cell_id) '.dat']); %x,y,z coordinates (in um) of somatic section (of a single neuron) 

    %To get true location of somatic section of each neuron within the spherical
    %volume
    soma = []; axon = []; loc_all = [];
    
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

    if(input_data.get_axon_dendrite_locs==1)
        intx=load(['intx_' num2str(cell_id) '.dat']); %x-coordinates (in um) of all neural section (dendrite, soma and axon) of a single neuron 
        inty=load(['inty_' num2str(cell_id) '.dat']); %y-coordinates (in um) of all neural section (dendrite, soma and axon) of a single neuron
        intz=load(['intz_' num2str(cell_id) '.dat']); %z-coordinates (in um) of all neural section (dendrite, soma and axon) of a single neuron


        x_axon=load(['x_axon_' num2str(cell_id) '.dat']); %x-coordinates (in um) of all axonal sections (of a single neuron)
        y_axon=load(['y_axon_' num2str(cell_id) '.dat']); %y-coordinates (in um) of all axonal sections (of a single neuron)
        z_axon=load(['z_axon_' num2str(cell_id) '.dat']); %z-coordinates (in um) of all axonal sections (of a single neuron)


    %     To get true location of axon section of each neuron within the spherical
    %     volume
        for k=1:length(realx)

            net_ad_x=((realx(k)+soma_coord(1))*cos(realang(k)))-((realz(k)+soma_coord(3))*sin(realang(k)))-realx(k);
            net_ad_y=realy(k)+soma_coord(2)-realy(k);
            net_ad_z=((realx(k)+soma_coord(1))*sin(realang(k)))+((realz(k)+soma_coord(3))*cos(realang(k)))-realz(k);

            for i=1:length(x_axon)

                x1=((realx(k)+x_axon(i))*cos(realang(k)))-((realz(k)+z_axon(i))*sin(realang(k)));
                y1=realy(k)+y_axon(i);
                z1=((realx(k)+x_axon(i))*sin(realang(k)))+((realz(k)+z_axon(i))*cos(realang(k)));

                axon(k).coord(i,1)=x1-net_ad_x;
                axon(k).coord(i,2)=y1-net_ad_y;
                axon(k).coord(i,3)=z1-net_ad_z;
            end
        end



    %     To get true location of all neural sections (i.e. dendrite, axon and soma) of each neuron within the spherical
    %     volume
        for k=1:length(realx)

            net_ad_x=((realx(k)+soma_coord(1))*cos(realang(k)))-((realz(k)+soma_coord(3))*sin(realang(k)))-realx(k);
            net_ad_y=realy(k)+soma_coord(2)-realy(k);
            net_ad_z=((realx(k)+soma_coord(1))*sin(realang(k)))+((realz(k)+soma_coord(3))*cos(realang(k)))-realz(k);

            for i=1:length(intx)

                x1=((realx(k)+intx(i))*cos(realang(k)))-((realz(k)+intz(i))*sin(realang(k)));
                y1=realy(k)+inty(i);
                z1=((realx(k)+intx(i))*sin(realang(k)))+((realz(k)+intz(i))*cos(realang(k)));

                loc_all(k).coord(i,1)=x1-net_ad_x;
                loc_all(k).coord(i,2)=y1-net_ad_y;
                loc_all(k).coord(i,3)=z1-net_ad_z;
            end
        end
    end
    
    % package outputs
    output_data = [];
    
    output_data.soma = soma;
    output_data.axon = axon;
    output_data.loc_all = loc_all;
    output_data.diam = input_data.diam;
    output_data.folderpath = input_data.folderpath;
    output_data.amp = input_data.amp;
    output_data.cell_id = input_data.cell_id;
    output_data.get_axon_dendrite_locs = input_data.get_axon_dendrite_locs;
    
    
end

