%% makes spatial activation prediction for a Histed brain and a Stoney brain
% so that analysis from Karthik's biophysical model can be put in
% perspective

% all we do is generate a set of neurons in a cortical column, then choose
% which neurons are activated based on rules (either Histed or Stoney
% rule)

    column_size = [2000,2000,2000]/1000; % mm, height, width, depth
    num_neurons = 1000000;

    neuron_locs = generateCorticalColumn(column_size, num_neurons);
    neuron_dists = sqrt(sum(neuron_locs.^2,2));
    neuron_angles = atan2(neuron_locs(:,3),neuron_locs(:,2));
% sample activated population at different amplitudes using either a Stoney
% rule (square root activation based on distance) or a Histed rule (fill in
% volume, but volume size does not change)
    % stoney rule
    stoney_input_data = [];
    stoney_input_data.rule = 'Stoney';
    stoney_input_data.amp = [10:10:100];
    stoney_input_data.k =  1292; % uA/(mm^2)
    stoney_input_data.locs = neuron_locs;
    
    stoney_act_mask = simulateActivatedPopulation(stoney_input_data);

    % histed rule
    histed_input_data = [];
    histed_input_data.rule = 'Histed';
    histed_input_data.amp = [10:10:100];
    histed_input_data.max_amp = 100;
    histed_input_data.min_amp = 0;
    histed_input_data.min_p = 0;
    histed_input_data.max_p = 0.4;
    histed_input_data.locs = neuron_locs;
    
    histed_act_mask = simulateActivatedPopulation(histed_input_data);
    
%% look at metrics/plots for both data sets

%% plot distance vs activation for each case
    figure()
    for i = 1:numel(histed_input_data.amp)
        ax=subplot(ceil(numel(histed_input_data.amp)/2),2,i,polaraxes);
        polarplot(ax,neuron_angles(histed_act_mask(:,i)),neuron_dists(histed_act_mask(:,i)),'.')
        rlim([0,max(neuron_dists)])
    end
    
    figure()
    for i = 1:numel(stoney_input_data.amp)
        ax=subplot(ceil(numel(stoney_input_data.amp)/2),2,i,polaraxes);
        polarplot(ax,neuron_angles(stoney_act_mask(:,i)),neuron_dists(stoney_act_mask(:,i)),'.')
        rlim([0,max(neuron_dists)])
    end
    
%% look at density of activated neurons at different radiuses across
% amplitudes
    
    
    r_list = linspace(0.1,1,6);
    color_list = inferno(numel(r_list)+1);
    
    stoney_dens = zeros(numel(stoney_input_data.amp),numel(r_list));
    histed_dens = zeros(numel(histed_input_data.amp),numel(r_list));
    
    for i_r = 1:numel(r_list)
        % computer density of activated cells within each radius       
        in_radius_mask = neuron_dists < r_list(i_r);
        
        for type = 1:2
            if(type==1)
                cell_act_mask = in_radius_mask & stoney_act_mask;
            else
                cell_act_mask = in_radius_mask & histed_act_mask;
            end
            
            num_cells_act = sum(cell_act_mask,1);
            
            if(type==1)
                stoney_dens(:,i_r) = num_cells_act./(4/3*pi*(r_list(i_r).^3));
%                 stoney_dens(:,i_r) = num_cells_act./sum(in_radius_mask);
            else
                histed_dens(:,i_r) = num_cells_act./(4/3*pi*(r_list(i_r).^3));
%                 histed_dens(:,i_r) = num_cells_act./sum(in_radius_mask);
            end
        end
    end

    
    figure(); hold on;
    set(gca,'ColorOrder',color_list,'ColorOrderIndex',1)
    plot(stoney_input_data.amp,stoney_dens,'--')
    set(gca,'ColorOrderIndex',1);
    plot(histed_input_data.amp,histed_dens,'-')
    
    