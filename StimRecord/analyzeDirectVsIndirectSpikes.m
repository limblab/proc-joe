%% estimates number of direct and indirectly evoked spikes at different amps. 
% This is based on Stoney's results, my non-stimulated channel data, and
% some reasonable assumptions

% load in exponential fit from my non-stimulated channel data and also set
% Stoney's parameters (k), as well as cell density parameter

    load('D:\Lab\Data\StimArtifact\Han_Duncan_spread_fits.mat');
    k = 1292; % uA/mm^2, 1292 is from Stoney 1968, could use other params
    cell_density = 75; % cells/0.001 mm^3, from Tehovnik 1995 
    max_r = 5; % mm, computational load depends on this value
    dr = 0.0001; % mm, computational load depends on this value

% estimate number of direct and indirect spikes in a volume
% estimating spikes in a volume requires a probability of activation
% function. For direct spikes, this will be sigmoidal. For indirect spikes,
% this will be exponential (provided by data)


% to estimate direct spikes, use square root relationship to get 50% point,
% and other sigmoidal parameters for rest of distances....Need to assume
% sigmoidal parameters
    amps = [15,30,60,100];
    sim_input_data = []; direct_data = [];  
    sim_input_data.min_r = 0; % for direct spikes
    sim_input_data.max_r = max_r; % mm
    sim_input_data.dr = dr; % mm
    sim_input_data.density = cell_density; % cells/0.001 mm^3...
    
    a = 20; % a controls how steep the function is at the 50% point (bigger = steeper)

    for i_amp = 1:numel(amps)
        % get activation function based on Stoney's results and sigmoid
        % params....
        
        b = (amps(i_amp)/k)^(1/2); % b is the 50% point from Stoney's results, in mm
        
        a_func = @(x) 1-(1./(1+exp(-a*(x-b))));
        
        sim_input_data.activation_func = a_func;
       
        temp_data = simulateNumberEvokedSpikes(sim_input_data);
        if(isempty(direct_data))
            direct_data = temp_data;
        else
            direct_data(i_amp) = temp_data;
        end
    end


    
% for indirect spikes, use exponential fit and a larger min_radius
    indirect_data = [];
    sim_input_data = [];    
    sim_input_data.min_r = 0.2; % mm
    sim_input_data.max_r = 5; % mm
    sim_input_data.dr = dr; % mm
    sim_input_data.density = cell_density; % cells/0.001 mm^3...
    
    for i_amp = 1:numel(exp_fits)
        % get activation function based on Stoney's results and sigmoid
        % params....
                
        a_func = @(x) exp_fits{i_amp}.a*exp(-x./(exp_fits{i_amp}.b/1000)) + exp_fits{i_amp}.c;
        sim_input_data.activation_func = a_func;
       
        temp_data = simulateNumberEvokedSpikes(sim_input_data);
        if(isempty(indirect_data))
            indirect_data = temp_data;
        else
            indirect_data(i_amp) = temp_data;
        end
    end

    

%%



















